#!/usr/bin/env python3
"""Shared parsers and writers for orthology-program group-table adapters."""
from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path

from msap_io import atomic_output, check_output_path

MISSING_GROUP_VALUES = {"", "*", "-", "NA", "N/A"}


def taxon_from_fasta_label(label: str) -> str:
    """Match the conventional taxon label used by FASTA-per-species tools."""
    name = Path(label.strip()).name
    return name.rsplit(".", 1)[0] if "." in name else name


def split_genes(value: str, *, source: str, group_id: str, taxon: str) -> list[str]:
    value = value.strip()
    if value in MISSING_GROUP_VALUES:
        return []
    genes = [gene.strip() for gene in value.split(",")]
    if any(not gene for gene in genes):
        raise ValueError(f"{source}: {group_id}, {taxon}: malformed comma-separated gene list")
    if len(set(genes)) != len(genes):
        raise ValueError(f"{source}: {group_id}, {taxon}: duplicated gene ID")
    return genes


def write_standard_table(output: Path, taxa: list[str], groups: list[tuple[str, dict[str, list[str]]]], inputs=()) -> None:
    if not taxa or len(set(taxa)) != len(taxa):
        raise ValueError("Taxon names are empty or duplicated after parsing input headers")
    seen = set()
    for group_id, members in groups:
        if not group_id or group_id in seen:
            raise ValueError(f"Duplicate or empty orthogroup ID: {group_id!r}")
        seen.add(group_id)
        if set(members).difference(taxa):
            raise ValueError(f"{group_id}: member table contains an unknown taxon")
    check_output_path(output, inputs=inputs)
    with atomic_output(output, inputs=inputs) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Orthogroup", *taxa])
        for group_id, members in groups:
            writer.writerow([group_id, *[", ".join(members.get(taxon, [])) for taxon in taxa]])


def parse_matrix_table(path: Path, *, metadata_columns: int, header_comment: bool, group_prefix: str) -> tuple[list[str], list[tuple[str, dict[str, list[str]]]]]:
    """Parse a table whose final columns are per-species comma-separated IDs."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Empty group table: {path}")
        if header_comment:
            header[0] = header[0].lstrip("#").strip()
        if len(header) <= metadata_columns:
            raise ValueError(f"{path}: expected {metadata_columns} metadata columns followed by taxon columns")
        taxa = [taxon_from_fasta_label(value) for value in header[metadata_columns:]]
        groups = []
        for line_number, row in enumerate(reader, 2):
            if not row or not any(value.strip() for value in row) or row[0].lstrip().startswith("#"):
                continue
            if len(row) != len(header):
                raise ValueError(f"{path}:{line_number}: expected {len(header)} tab-separated columns, found {len(row)}")
            group_id = row[0].strip() if metadata_columns == 1 else f"{group_prefix}{len(groups) + 1:07d}"
            members = {
                taxon: split_genes(value, source=str(path), group_id=group_id, taxon=taxon)
                for taxon, value in zip(taxa, row[metadata_columns:])
            }
            groups.append((group_id, members))
    if not groups:
        raise ValueError(f"No orthogroups found in {path}")
    return taxa, groups


def parse_prefixed_groups(path: Path, *, separator: str, strip_prefix: bool) -> tuple[list[str], list[tuple[str, dict[str, list[str]]]]]:
    """Parse ``GROUP: taxon|gene ...`` lines used by OrthoMCL and OMA."""
    taxa = OrderedDict()
    groups = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                raise ValueError(f"{path}:{line_number}: expected 'GROUP: member member ...'")
            group_id, members_text = line.split(":", 1)
            group_id = group_id.strip()
            members = members_text.split()
            if not group_id or not members:
                raise ValueError(f"{path}:{line_number}: group ID or members are empty")
            table = {}
            for member in members:
                if separator not in member:
                    raise ValueError(f"{path}:{line_number}: member lacks separator {separator!r}: {member}")
                taxon, gene = member.split(separator, 1)
                if not taxon or not gene:
                    raise ValueError(f"{path}:{line_number}: invalid taxon/gene member: {member}")
                taxa.setdefault(taxon, None)
                table.setdefault(taxon, []).append(gene if strip_prefix else member)
            for taxon, genes in table.items():
                if len(set(genes)) != len(genes):
                    raise ValueError(f"{path}:{line_number}: duplicate gene in {group_id}, {taxon}")
            groups.append((group_id, table))
    if not groups:
        raise ValueError(f"No orthogroups found in {path}")
    return list(taxa), groups


def parse_prefixed_tab_groups(path: Path, *, separator: str, strip_prefix: bool) -> tuple[list[str], list[tuple[str, dict[str, list[str]]]]]:
    """Parse ``GROUP<TAB>taxon|gene<TAB>...`` rows used by OMA groups."""
    taxa = OrderedDict()
    groups = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, 1):
            if not row or not any(value.strip() for value in row) or row[0].lstrip().startswith("#"):
                continue
            group_id, members = row[0].strip(), [member.strip() for member in row[1:] if member.strip()]
            if not group_id or not members:
                raise ValueError(f"{path}:{line_number}: group ID or members are empty")
            table = {}
            for member in members:
                if separator not in member:
                    raise ValueError(f"{path}:{line_number}: member lacks separator {separator!r}: {member}")
                taxon, gene = member.split(separator, 1)
                if not taxon or not gene:
                    raise ValueError(f"{path}:{line_number}: invalid taxon/gene member: {member}")
                taxa.setdefault(taxon, None)
                table.setdefault(taxon, []).append(gene if strip_prefix else member)
            for taxon, genes in table.items():
                if len(set(genes)) != len(genes):
                    raise ValueError(f"{path}:{line_number}: duplicate gene in {group_id}, {taxon}")
            groups.append((group_id, table))
    if not groups:
        raise ValueError(f"No orthogroups found in {path}")
    return list(taxa), groups
