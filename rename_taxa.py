#!/usr/bin/env python3
"""Rename FASTA sequence IDs using a mapping table and safe-name rules."""
from __future__ import annotations
import argparse, csv, re, sys
from pathlib import Path
from Bio import SeqIO
from msap_io import atomic_output, check_output_path

def normalized_taxon_id(identifier):
    return identifier.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")


def read_mapping(mapping_file):
    mapping = {}
    with open(mapping_file, encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError("Mapping table requires two tab-separated columns: old_id and new_id")
            old_id, new_id = normalized_taxon_id(row[0]), normalized_taxon_id(row[1])
            if old_id in mapping and mapping[old_id] != new_id:
                raise ValueError(f"Conflicting mapping after OrthoFinder normalization: {row[0]}")
            mapping[old_id] = new_id
    return mapping


def rename_fasta(infile, outfile, mapping_file=None, map_output=None, allow_duplicate_ids=False,
                 mapping=None):
    if mapping_file is not None and mapping is not None:
        raise ValueError("Provide either mapping_file or mapping, not both")
    inputs = [infile] + ([mapping_file] if mapping_file else [])
    check_output_path(outfile, inputs)
    check_output_path(map_output, inputs + [outfile])
    records = list(SeqIO.parse(infile, "fasta"))
    if not records:
        raise ValueError("No sequences found in input FASTA.")
    mapping = read_mapping(mapping_file) if mapping_file else (mapping or {})
    result, seen, map_rows = [], set(), []
    for record in records:
        old = normalized_taxon_id(record.id)
        new = mapping.get(old, old)

        if not new: raise ValueError(f"{old}: renamed ID is empty")
        if any(char.isspace() for char in new):
            raise ValueError(f'{old}: FASTA IDs must not contain whitespace, even with --keep-special: {new!r}')
        if new in seen and not allow_duplicate_ids:
            raise ValueError(f"Duplicate ID after renaming: {new}")
        seen.add(new); record.id = new; record.name = new; record.description = ""
        result.append(record); map_rows.append((old, new))
    with atomic_output(outfile, inputs=[infile]) as handle: SeqIO.write(result, handle, "fasta")
    if map_output:
        with atomic_output(map_output) as handle:
            handle.write("original_id\trenamed_id\n")
            handle.writelines(f"{old}\t{new}\n" for old, new in map_rows)

def main(argv=None):
    parser = argparse.ArgumentParser(description="Rename FASTA sequence IDs for phylogenetic software.", add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="""Examples:
  rename_taxa.py -i input.fasta -o renamed.fasta -m taxa.tsv
  rename_taxa.py -i input.fasta -o renamed.fasta --map-output renamed.tsv

Mapping format: tab-separated original_id and renamed_id, one pair per line.
Special characters are replaced with underscores by default.""")
    required=parser.add_argument_group("required arguments"); optional=parser.add_argument_group("renaming options")
    required.add_argument("-i","--input",required=True,metavar="FASTA",help="Input FASTA file.")
    required.add_argument("-o","--output",required=True,metavar="FASTA",help="Output FASTA file.")
    optional.add_argument("-m","--mapping","--map",metavar="TSV",help="Two-column ID mapping table.")
    optional.add_argument("--map-output",metavar="TSV",help="Write original-to-renamed IDs.")
    optional.add_argument("--allow-duplicate-ids", action="store_true",
                        help="Write duplicate renamed IDs for a downstream validator; normally duplicates are an error.")
    optional.add_argument("--keep-special", action="store_true", help="Deprecated compatibility option; OrthoFinder substitutions still apply.")
    optional.add_argument("-h","--help",action="help",help="Show this help message and exit.")
    optional.add_argument("-v","--version",action="version",version="v1.0.0",help="Show program's version number and exit.")
    args=parser.parse_args(argv)
    try: rename_fasta(args.input,args.output,args.mapping,args.map_output,args.allow_duplicate_ids)
    except (OSError,ValueError) as error: print(f"rename_taxa: {error}",file=sys.stderr); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
