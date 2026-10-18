#!/usr/bin/env python3
"""Build an offline English HTML report from one PhyloPrep output directory.

This file is independently runnable with Python's standard library. No source
FASTA, mapping file, plotting package, network access or analysis rerun is needed.
"""

import argparse
from array import array
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile

__version__ = "1.0.2"
UNAVAILABLE = "Unavailable"
BLOCK = 8190
ALIGNMENT_SUFFIX = re.compile(
    r"\.(?:mafft|muscle|prank|clustalw2|macse)\.(?:codon|prot|nucl)"
    r"(?:\.trimmed)?\.aln(?:\.fasta)?$")
MATRIX_NAMES = {"supermatrix.fasta", "fourfold.fasta"} | {
    f"codon{i}{suffix}.supermatrix.fasta" for i, suffix in ((1, "st"), (2, "nd"), (3, "rd"))}


def basename(value):
    return str(value).replace("\\", "/").rsplit("/", 1)[-1]


def locus(value):
    return ALIGNMENT_SUFFIX.sub("", basename(value))


def partition_id(value):
    name = re.sub(r"[^A-Za-z0-9_]", "_", locus(value))
    return name if name and name[0].isalpha() else "gene_" + name


def rows(path, required=()):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"Missing report columns in {path.name}: {', '.join(required)}")
        yield from reader


def path_list(path):
    if not path.is_file():
        return None
    with path.open(encoding="utf-8-sig") as handle:
        return {line.strip() for line in handle if line.strip() and not line.startswith("#")}


def read_json(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def visible_files(directory):
    """Do not descend into publication backups, staging trees or private caches."""
    for current, dirs, files in os.walk(directory):
        dirs[:] = sorted(d for d in dirs if not d.startswith(".") and not Path(current, d).is_symlink())
        for name in sorted(files):
            path = Path(current, name)
            if not name.startswith(".") and not path.is_symlink():
                yield path


def managed(path):
    manifest = path.parent / ".phyloprep-products.json"
    if manifest.is_file():
        names = read_json(manifest)
        if not isinstance(names, list) or any(not isinstance(n, str) for n in names):
            raise ValueError(f"Invalid product manifest: {manifest}")
        return path.name in names
    return True


class FastaSpool:
    """Bounded parsing even for unwrapped matrices; only IDs/offsets stay in RAM."""

    def __init__(self, source, directory):
        self.file = tempfile.TemporaryFile(dir=directory)
        self.records = []
        self.length = None
        identifier, start, length = None, 0, 0
        seen = set()

        def finish():
            if identifier is None:
                return
            if not length or (self.length is not None and self.length != length):
                raise ValueError(f"Empty or unequal-length sequence: {identifier}")
            self.length = length
            self.records.append((identifier, start))

        try:
            with source.open("rb") as handle:
                line_start = True
                while True:
                    part = handle.readline(65536)
                    if not part:
                        break
                    if line_start and part.startswith(b">"):
                        finish()
                        if not part.endswith(b"\n"):
                            raise ValueError("FASTA header exceeds 64 KiB")
                        fields = part[1:].decode("utf-8").split()
                        if not fields or fields[0] in seen:
                            raise ValueError("Empty or duplicate FASTA identifier")
                        identifier = fields[0]
                        seen.add(identifier)
                        start, length = self.file.tell(), 0
                    else:
                        chunk = b"".join(part.split()).upper()
                        if chunk:
                            if identifier is None:
                                raise ValueError("Sequence before FASTA header")
                            if not chunk.isascii():
                                raise ValueError("Non-ASCII sequence symbols")
                            self.file.write(chunk)
                            length += len(chunk)
                    line_start = part.endswith(b"\n")
            finish()
            if not self.records:
                raise ValueError("No sequences in matrix")
            self.file.flush()
        except BaseException:
            self.file.close()
            raise

    def close(self):
        self.file.close()


def matrix_stats(path, seqtype, temporary):
    spool = FastaSpool(path, temporary)
    try:
        if seqtype == "codon" and spool.length % 3:
            raise ValueError("Full codon matrix length is not a multiple of three")
        alphabet = b"ACDEFGHIKLMNPQRSTVWY" if seqtype == "prot" else b"ACGT"
        indices = {symbol: i for i, symbol in enumerate(alphabet)}
        characters = Counter()
        sites = Counter()
        for start in range(0, spool.length, BLOCK):
            width = min(BLOCK, spool.length - start)
            counts = [array("I", [0]) * width for _ in alphabet]
            for _, offset in spool.records:
                spool.file.seek(offset + start)
                chunk = spool.file.read(width)
                characters.update(chunk)
                for column, symbol in enumerate(chunk):
                    index = indices.get(symbol)
                    if index is not None:
                        counts[index][column] += 1
            for column in zip(*counts):
                states = sum(n > 0 for n in column)
                sites["unresolved_only" if states == 0 else "constant" if states == 1 else "variable"] += 1
                sites["pis"] += sum(n >= 2 for n in column) >= 2
                sites["singleton"] += states == 2 and min(n for n in column if n) == 1
        total = spool.length * len(spool.records)
        resolved = sum(characters[symbol] for symbol in alphabet)
        gaps, missing = characters[ord("-")], characters[ord("?")]
        return dict(taxa=[name for name, _ in spool.records], length=spool.length,
                    unit="aa" if seqtype == "prot" else "nt", cells=total,
                    complete_codons=spool.length // 3 if seqtype == "codon" and spool.length % 3 == 0 else None,
                    gap_characters=gaps, missing_characters=missing,
                    ambiguous_characters=total - resolved - gaps - missing,
                    gap_ratio=gaps / total, missing_ratio=missing / total,
                    ambiguous_ratio=(total - resolved - gaps - missing) / total,
                    gc_ratio=(characters[ord("G")] + characters[ord("C")]) / resolved
                    if seqtype != "prot" and resolved else None,
                    **{key: sites[key] for key in ("constant", "variable", "pis", "singleton", "unresolved_only")})
    finally:
        spool.close()


def partitions(path, length):
    source = path.with_name("supermatrix_part_iqtree.txt")
    if not source.is_file() or not managed(source):
        return None
    result, previous = {}, 0
    with source.open(encoding="utf-8") as handle:
        for line in handle:
            match = re.match(r"\s*charset\s+(\S+)\s*=\s*(\d+)\s*-\s*(\d+)\s*;", line, re.I)
            if match:
                name, start, end = match.group(1), int(match.group(2)), int(match.group(3))
                if name in result or start != previous + 1 or end < start:
                    raise ValueError("Non-contiguous or duplicate gene partitions")
                result[name], previous = (start, end), end
    if not result or previous != length:
        raise ValueError("Gene partition coordinates do not span the full matrix")
    return result


def padding_stats(path, parts, taxa, db):
    source = path.with_suffix(".missing.tsv")
    if not source.is_file() or not managed(source) or parts is None:
        return None
    ids = set(taxa)
    db.execute("DELETE FROM padding")
    for row in rows(source, ("locus", "taxon", "start", "end", "length", "fill_character")):
        name, taxon = row["locus"], row["taxon"]
        start, end, length = (int(row[col]) for col in ("start", "end", "length"))
        if (name not in parts or taxon not in ids or (start, end) != parts[name]
                or length != end - start + 1 or row["fill_character"] not in {"?", "-"}):
            raise ValueError("Missing-block metadata does not match matrix partitions/taxa")
        try:
            db.execute("INSERT INTO padding VALUES (?, ?, ?)", (name, taxon, length))
        except sqlite3.IntegrityError as error:
            raise ValueError("Duplicate missing taxon/locus block") from error
    counts = {taxon: (n, length) for taxon, n, length in db.execute(
        "SELECT taxon, count(*), sum(length) FROM padding GROUP BY taxon")}
    return dict(blocks=sum(n for n, _ in counts.values()),
                characters=sum(length for _, length in counts.values()),
                coverage=[dict(taxon=taxon, observed_loci=len(parts) - counts.get(taxon, (0, 0))[0],
                               total_loci=len(parts), artificial_characters=counts.get(taxon, (0, 0))[1])
                          for taxon in taxa], source=source.name)


def active_branch(group, parameters):
    """Current-run metadata excludes saved branches from previous workflow modes."""
    workflow = parameters.get('seqtype')
    if workflow not in {'nucl', 'prot', 'codon', 'pseudogene'}:
        return True
    parts = Path(group).parts
    if (workflow == 'pseudogene') != ('pseudogene' in parts):
        return False
    if parameters.get('notrim') and parts and parts[-1] == 'trimmed':
        return False
    representations = set(parts) & {'nucl', 'prot', 'codon'}
    if workflow in {'codon', 'pseudogene'}:
        return not representations - {'codon', 'prot'}
    return 'paired' not in parts and (not representations or workflow in representations)


def qc_summary(root, warnings, parameters=None):
    result = []
    for path in visible_files(root / "qc"):
        paired = path.name == "loci.details.tsv"
        if path.name != "alignment.details.tsv" and not paired:
            continue
        if not active_branch(path.parent.relative_to(root / 'qc'), parameters or {}):
            continue
        passed_list = path_list(path.with_name("alignment.pass.tsv")) if not paired else None
        failed_list = path_list(path.with_name("alignment.fail.tsv")) if not paired else None
        passed = {basename(p) for p in passed_list} if passed_list is not None else None
        failed = {basename(p) for p in failed_list} if failed_list is not None else None
        retained = path_list(path.with_name("independent.retained.tsv")) if not paired else None
        kept = {basename(p) for p in retained} if retained is not None else None
        hits, combinations, multiplicity = Counter(), Counter(), Counter()
        counts = Counter()
        seen = set()
        # Paired details carry the otherwise unavailable independent PIS thresholds.
        thresholds = {}
        if not paired:
            relative = path.parent.relative_to(root / "qc")
            pair_path = root / "qc"
            if "pseudogene" in relative.parts:
                pair_path /= "pseudogene"
            pair_path = pair_path / "paired" / path.parent.name / "loci.details.tsv"
            if pair_path.is_file():
                first = next(rows(pair_path), None)
                if first:
                    thresholds = {"codon": int(first["Minimum_codon_PIS"]),
                                  "prot": int(first["Minimum_protein_PIS"])}
        for row in rows(path):
            key = row["Locus"] if paired else basename(row["File"])
            if key in seen:
                warnings.append(f"Duplicate QC row ignored: {path.relative_to(root)} / {key}")
                continue
            seen.add(key)
            reasons = {r for r in row.get("Reasons" if paired else "Issues", "").split(";") if r and r != "NA"}
            if paired:
                status = row.get("Paired_status")
            elif kept is not None:
                status = "pass" if key in kept else "fail"
                threshold = thresholds.get(row.get("Seqtype"))
                if threshold is not None and int(row["Parsimony_informative_sites"]) < threshold:
                    reasons.add("below_independent_PIS_threshold")
            elif passed is not None and key in passed:
                status = "pass"
            elif failed is not None and key in failed:
                status = "fail"
            else:
                status = "not_evaluated"
            if status not in {"pass", "fail"}:
                status = "not_evaluated"
            counts[status] += 1
            if status == "fail":
                reasons = reasons or {"failure_reason_not_recorded"}
                hits.update(reasons)
                combinations[" + ".join(sorted(reasons))] += 1
                multiplicity[str(len(reasons))] += 1
        result.append(dict(branch=path.parent.relative_to(root / "qc").as_posix(),
                           evaluated=counts["pass"] + counts["fail"], passed=counts["pass"],
                           failed=counts["fail"], not_evaluated=counts["not_evaluated"],
                           reason_hits=dict(hits), failure_combinations=dict(combinations),
                           reasons_per_locus=dict(multiplicity),
                           selection="Paired" if paired else "Independent" if kept is not None else "Recorded QC",
                           source=path.relative_to(root).as_posix()))
    return result


def audit_summary(root, db):
    path = root / "reports" / "sequence_changes.tsv"
    if not path.is_file():
        return None
    db.execute("CREATE TABLE events (digest TEXT PRIMARY KEY, taxon TEXT, sequence_key TEXT, event TEXT, action TEXT, nt INTEGER, aa INTEGER, unknown INTEGER)")
    duplicates = 0
    for row in rows(path, ("input_file", "taxa_id", "event", "action")):
        digest = hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()
        taxon = row.get("mapped_taxa_id")
        if not taxon or taxon == "NA":
            taxon = row.get("taxa_id") or "Unassigned"
        if taxon == "NA":
            taxon = "Unassigned"
        sequence = json.dumps([row.get("input_file"), row.get("taxa_id")])
        # Source coordinate spans measure affected characters, not output length.
        lengths = []
        for field in ("source_nt_positions", "source_aa_positions"):
            value = row.get(field, "")
            if re.fullmatch(r"\d+(?:-\d+)?", value):
                ends = list(map(int, value.split("-")))
                lengths.append(ends[-1] - ends[0] + 1 if ends[0] > 0 and ends[-1] >= ends[0] else None)
            else:
                lengths.append(None)
        cursor = db.execute("INSERT OR IGNORE INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                            (digest, taxon, sequence if taxon != "Unassigned" else None,
                             row.get("event", "Unrecorded"), row.get("action", "Unrecorded"),
                             lengths[0], lengths[1], int(all(n is None for n in lengths))))
        duplicates += cursor.rowcount == 0
    taxa = [dict(taxon=t, events=n, affected_sequences=s, nt_characters=nt,
                 aa_characters=aa, events_without_coordinates=u)
            for t, n, s, nt, aa, u in db.execute(
                "SELECT taxon, count(*), count(DISTINCT sequence_key), sum(nt), sum(aa), sum(unknown) "
                "FROM events GROUP BY taxon ORDER BY count(*) DESC, taxon")]
    kinds = [dict(event=e, action=a, count=n) for e, a, n in db.execute(
        "SELECT event, action, count(*) FROM events GROUP BY event, action ORDER BY count(*) DESC, event, action")]
    by_taxon_event = [dict(taxon=t, event=e, count=n) for t, e, n in db.execute(
        "SELECT taxon, event, count(*) FROM events GROUP BY taxon, event ORDER BY taxon, event")]
    return dict(source="reports/sequence_changes.tsv", duplicates_ignored=duplicates,
                events=sum(t["events"] for t in taxa), taxa=taxa, types=kinds,
                by_taxon_event=by_taxon_event)


def processing_overview(root, manifest):
    listed = path_list(root / "inputs.pathlist")
    total = len(listed) if listed is not None else None
    inputs = manifest.get("inputs")
    if total is None and isinstance(inputs, list):
        if any("input_key" in entry for entry in inputs):
            total = sum(entry.get("input_key") == "inputs" and not entry.get("name", "").endswith(".pathlist")
                        for entry in inputs)
        else:
            # Older manifests did not identify mapping/list auxiliaries.
            total = None
    lists = list((root / "alignments" / "path-lists").glob("all.*.pathlist"))
    successful = set()
    for path in lists:
        successful.update(locus(p) for p in path_list(path) or ())
    # Lists are authoritative when present; scanning files could include stale loci.
    if not lists:
        for path in visible_files(root / "alignments"):
            if ALIGNMENT_SUFFIX.search(path.name) and ".trimmed.aln" not in path.name:
                successful.add(locus(path.name))
    return dict(input_loci=total, aligned_loci=len(successful) if lists or successful else None,
                failed_alignments=None, reused_loci=manifest.get("reused_loci"),
                alignment_count_source="Published alignment path lists" if lists else "Present alignment files (inferred)" if successful else UNAVAILABLE)


def inclusion(path):
    source = path.with_suffix(".report.tsv")
    if not source.is_file() or not managed(source):
        return None
    kept, skipped = {}, 0
    for row in rows(source):
        if row.get("status") in {"included", "padded"}:
            key = partition_id(row["gene"])
            if key in kept:
                raise ValueError("Duplicate included locus in concatenation report")
            kept[key] = row["gene"]
        elif row.get("status") == "skipped":
            skipped += 1
    return dict(kept=kept, skipped=skipped)


def build_report(root, run_metadata=None):
    root = Path(root).resolve()
    if not root.is_dir() or not any((root / name).is_dir() for name in ("alignments", "qc", "matrices", "reports")):
        raise ValueError("Input must be a PhyloPrep final output directory")
    warnings = []
    manifest_path = root / "phyloprep-run.json"
    manifest = run_metadata if run_metadata is not None else read_json(manifest_path) if manifest_path.is_file() else {}
    if not isinstance(manifest, dict):
        raise ValueError("The run manifest must contain a JSON object")
    parameters = manifest.get('parameters') or {}
    report = dict(schema="phyloprep-report-v1", report_version=__version__,
                  generated_at=datetime.now(timezone.utc).isoformat(), run_name=root.name,
                  overview=processing_overview(root, manifest), qc=qc_summary(root, warnings, parameters),
                  matrices=[], trimming_comparisons=[], warnings=warnings,
                  provenance={key: manifest[key] for key in (
                      "created_at", "module", "script_version", "source_sha256", "parameters", "equivalent_cli",
                      "runtime", "reused_loci", "matrix_padding", "engine_metadata") if key in manifest},
                  provenance_source=("Current PhyloPrep workflow" if run_metadata is not None
                                     else "phyloprep-run.json") if manifest else None)
    if manifest:
        report["provenance"]["checksum_inventory"] = {
            key: dict(files=len(manifest[key]), sha256_recorded=sum(bool(item.get("sha256")) for item in manifest[key]))
            for key in ("inputs", "outputs") if isinstance(manifest.get(key), list)}
    # All temporary sequence/event data stay on the filesystem beside the output.
    with tempfile.TemporaryDirectory(prefix=".phyloprep-report-", dir=root) as temporary:
        db = sqlite3.connect(str(Path(temporary) / "events.sqlite"))
        try:
            db.execute("PRAGMA cache_size=-2048")
            db.execute("PRAGMA temp_store=FILE")
            db.execute("CREATE TABLE padding (locus TEXT, taxon TEXT, length INTEGER, PRIMARY KEY(locus,taxon))")
            report["audit"] = audit_summary(root, db)
            groups = {}
            for path in visible_files(root / "matrices"):
                if path.name not in MATRIX_NAMES or not managed(path):
                    continue
                relative = path.relative_to(root / "matrices")
                if not active_branch(relative.parent, parameters):
                    continue
                seqtype = next((p for p in relative.parts[:-1] if p in {"nucl", "prot", "codon"}), None)
                if not seqtype:
                    warnings.append(f"Unknown sequence type; matrix not evaluated: {relative}")
                    continue
                kind = "full" if path.name == "supermatrix.fasta" else path.name.split(".")[0]
                actual_type = seqtype if kind == "full" else "nucl"
                key = relative.as_posix()
                print(f"[PhyloPrep report] Reading {key}", flush=True)
                entry = dict(matrix=key, group=relative.parent.as_posix(), kind=kind, seqtype=actual_type)
                try:
                    entry.update(matrix_stats(path, actual_type, temporary))
                    entry["status"] = "evaluated"
                except (ValueError, OSError) as error:
                    entry.update(status="not_evaluated", error=str(error))
                    warnings.append(f"{key}: {error}")
                    report["matrices"].append(entry)
                    continue
                parts = included = padding = None
                if kind == "full":
                    try:
                        parts = partitions(path, entry["length"])
                        included = inclusion(path)
                        if parts is not None and included is not None and set(parts) != set(included["kept"]):
                            raise ValueError("Concatenation report and gene partitions disagree")
                        padding = padding_stats(path, parts, entry["taxa"], db)
                    except (ValueError, KeyError, sqlite3.Error) as error:
                        warnings.append(f"{key}: {error}; locus/coverage metadata not evaluated")
                        parts = included = padding = None
                    entry.update(loci=len(parts) if parts is not None else len(included["kept"]) if included else None,
                                 skipped_after_qc=included["skipped"] if included is not None else None,
                                 padding=padding)
                    groups[entry["group"]] = (entry, parts, included)
                else:
                    entry.update(loci=None, padding=None)
                report["matrices"].append(entry)
            # Derived positions share full-matrix loci; fourfold contributors are not recorded.
            for entry in report["matrices"]:
                if entry["kind"].startswith("codon") and entry["group"] in groups:
                    entry["loci"] = groups[entry["group"]][0]["loci"]
            for group, (raw, raw_parts, _) in groups.items():
                if group.split("/")[-1] != "raw":
                    continue
                trimmed_group = group.rsplit("/", 1)[0] + "/trimmed"
                if trimmed_group not in groups:
                    continue
                trimmed, trimmed_parts, _ = groups[trimmed_group]
                comparison = dict(raw=group, trimmed=trimmed_group, common_loci=None,
                                  raw_common_columns=None, trimmed_common_columns=None)
                if raw_parts is not None and trimmed_parts is not None:
                    common = set(raw_parts) & set(trimmed_parts)
                    comparison.update(common_loci=len(common), raw_only=len(set(raw_parts) - common),
                                      trimmed_only=len(set(trimmed_parts) - common),
                                      raw_common_columns=sum(raw_parts[k][1] - raw_parts[k][0] + 1 for k in common),
                                      trimmed_common_columns=sum(trimmed_parts[k][1] - trimmed_parts[k][0] + 1 for k in common),
                                      unit=raw["unit"])
                report["trimming_comparisons"].append(comparison)
        finally:
            db.close()
    return report


def esc(value):
    return html.escape(str(value), quote=True)


def display(value):
    if value is None:
        return UNAVAILABLE
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def pct(value):
    return UNAVAILABLE if value is None else f"{100 * value:.2f}%"


def table(headers, data, title="Table", row_keys=None):
    body = []
    for index, row in enumerate(data):
        key = f' data-matrix-key="{esc(row_keys[index])}"' if row_keys is not None else ''
        body.append(f'<tr{key}>' + ''.join(
            f'<td data-label="{esc(header)}">{esc(display(value))}</td>'
            for header, value in zip(headers, row)) + '</tr>')
    if not body:
        return '<p class="empty">No records available for this table.</p>'
    return f'<div class="table-scroll data-table" data-title="{esc(title)}"><table><thead><tr>' + ''.join(
        f'<th scope="col">{esc(h)}</th>' for h in headers) + '</tr></thead><tbody>' + ''.join(body) + '</tbody></table></div>'


REASON_LABELS = {
    'duplicate_ids': 'Duplicate sequence identifiers', 'length_mismatch': 'Unequal sequence lengths',
    'invalid_characters': 'Invalid sequence symbols', 'empty_sequence': 'Empty sequence',
    'all_gap_columns': 'Columns containing only gaps', 'frame_error': 'Invalid reading frame or codon gaps',
    'terminal_stop': 'Unexpected terminal stop', 'internal_stop': 'Internal stop detected',
    'short_alignment': 'Alignment below minimum length', 'gap_ratio': 'Excessive alignment gaps',
    'N_ratio': 'Excessive nucleotide ambiguity', 'X_ratio': 'Excessive protein ambiguity',
    'parsimony_informative_sites': 'Insufficient parsimony-informative sites',
    'below_independent_PIS_threshold': 'Below independent PIS threshold',
    'codon_qc_failed': 'Codon QC failed', 'protein_qc_failed': 'Protein QC failed',
    'pis_policy_failed': 'Paired PIS rule not met', 'failure_reason_not_recorded': 'Failure reason not recorded',
}


def reason_label(code):
    return REASON_LABELS.get(code, code.replace('_', ' ').capitalize())


def combination_label(value):
    return ' + '.join(reason_label(code) for code in value.split(' + '))


def group_label(group):
    parts = group.split('/')
    labels = ['Pseudogene'] if 'pseudogene' in parts else []
    labels.append('Paired QC' if 'paired' in parts else 'Independent QC')
    labels.extend({'nucl': 'Nucleotide', 'prot': 'Protein', 'codon': 'Codon'}[part]
                  for part in parts if part in {'nucl', 'prot', 'codon'})
    labels.extend(part.title() for part in parts if part in {'raw', 'trimmed'})
    return ' · '.join(labels)


def matrix_label(matrix):
    return group_label(matrix['group']) + (' · ' + derived_name(matrix['kind']) if matrix['kind'] != 'full' else '')


def event_label(code):
    names = {'terminal_stop': 'Terminal stop processing', 'protein_terminal_stop': 'Protein terminal stop processing',
             'terminal_partial_codon': 'Incomplete terminal codon removal', 'internal_stop': 'Internal stop masking',
             'protein_internal_stop': 'Protein internal stop masking', 'nucleotide_ambiguity': 'Nucleotide ambiguity normalization',
             'protein_unsupported_residue': 'Non-standard amino-acid normalization', 'uracil_to_thymine': 'Uracil-to-thymine conversion',
             'validation_error': 'Sequence validation event', 'invalid_frame': 'Invalid reading frame'}
    return names.get(code, code.replace('_', ' ').capitalize())


def event_category(code):
    if code in {'terminal_stop', 'protein_terminal_stop', 'terminal_partial_codon'}:
        return 'Terminal cleanup'
    if code in {'internal_stop', 'protein_internal_stop'}:
        return 'Internal stops'
    if code in {'nucleotide_ambiguity', 'protein_unsupported_residue', 'uracil_to_thymine'}:
        return 'Symbol normalization'
    if 'frame' in code:
        return 'Frame-related events'
    if 'renam' in code or 'identifier' in code:
        return 'Identifier changes'
    return 'Other recorded events'


def readable_time(value):
    try:
        parsed = datetime.fromisoformat(value.replace('Z', '+00:00'))
        return parsed.strftime('%d %b %Y, %H:%M') + (' UTC' if parsed.utcoffset() is not None and parsed.utcoffset().total_seconds() == 0 else parsed.strftime(' %z'))
    except (ValueError, AttributeError):
        return value


def setting_value(key, value, parameters):
    mode = parameters.get('seqtype')
    if key in {'min_nucl_pis', 'min_protein_pis'} and value is None:
        if (key == 'min_nucl_pis' and mode == 'prot') or (key == 'min_protein_pis' and mode == 'nucl'):
            return 'Not applicable'
        return '0 (default)'
    if key == 'min_length' and value is not None:
        return f'{value} nt' if mode == 'nucl' else f'{value} aa' if mode == 'prot' else f'{value} codons / aa'
    if isinstance(value, (list, dict)):
        return json.dumps(value)
    return 'Not explicitly set' if value is None else value


def bar_chart(title, values, denominator=None, limit=12, labeler=None):
    ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
    if not ordered:
        return '<p class="empty">No recorded observations.</p>'
    omitted = max(0, len(ordered) - limit)
    ordered = ordered[:limit]
    width, height = 900, 48 + 48 * len(ordered)
    maximum = max(1, denominator if denominator is not None else max(v for _, v in ordered))
    parts = [f'<svg class="desktop-chart" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" font-family="sans-serif" role="img" aria-label="{esc(title)}">',
             f'<title>{esc(title)}</title>']
    for i, (label, value) in enumerate(ordered):
        y = 24 + i * 48
        readable = labeler(label) if labeler else label
        short = readable if len(readable) <= 48 else readable[:45] + "…"
        suffix = f" ({value / denominator:.1%})" if denominator else ""
        parts.extend([f'<text x="0" y="{y + 14}" font-size="12" fill="#334155">{esc(short)}</text>',
                      f'<rect x="365" y="{y}" width="{value / maximum * 400:.1f}" height="24" rx="4" fill="#197c80"><title>{esc(readable)} ({esc(label)}): {value:,}{esc(suffix)}</title></rect>',
                      f'<text x="775" y="{y + 17}" font-size="12" fill="#0f172a">{value:,}{esc(suffix)}</text>'])
    parts.append('</svg>')
    parts.append('<div class="mobile-bars">')
    for label, value in ordered:
        readable = labeler(label) if labeler else label
        suffix = f' ({value / denominator:.1%})' if denominator else ''
        parts.append(f'<div class="bar-item"><div class="bar-caption"><span>{esc(readable)}</span><strong>{value:,}{esc(suffix)}</strong></div>'
                     f'<div class="bar-track"><span style="width:{value / maximum * 100:.3f}%"></span></div></div>')
    parts.append('</div>')
    if omitted:
        parts.append(f'<p class="muted">Top {limit} categories shown; {omitted} additional categories are included in the tables or summary JSON.</p>')
    return ''.join(parts)


def stacked_chart(title, series, legend, limit=12):
    shown = series[:limit]
    if not shown:
        return '<p class="empty">No recorded observations.</p>'
    fragments = ['<div class="legend">' + ''.join(f'<span><i style="background:{color}"></i>{esc(label)}</span>' for label, color in legend) + '</div>']
    svg = [f'<svg class="export-chart" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 {80 + len(shown) * 62}" font-family="sans-serif" role="img"><title>{esc(title)}</title>']
    legend_x = 0
    for label, color in legend:
        svg.append(f'<rect x="{legend_x}" y="4" width="12" height="12" fill="{color}"/><text x="{legend_x + 18}" y="15" font-size="12">{esc(label)}</text>')
        legend_x += 40 + 7 * len(label)
    for index, (label, values) in enumerate(shown):
        total = sum(values)
        labels = ', '.join(f'{legend[i][0]}: {count:,}' for i, count in enumerate(values))
        fragments.append(f'<div class="stack-row"><div class="bar-caption"><strong>{esc(label)}</strong><span>{esc(labels)}</span></div><div class="stack-track" role="img" aria-label="{esc(label + "; " + labels)}">')
        x, y = 0, 60 + index * 62
        svg.append(f'<text x="0" y="{y - 10}" font-size="13">{esc(label)} — {esc(labels)}</text>')
        for i, count in enumerate(values):
            width = count / total * 100 if total else 0
            if count:
                fragments.append(f'<span style="width:{width:.3f}%;background:{legend[i][1]}"><span class="sr-only">{esc(legend[i][0])}: {count:,}</span></span>')
                svg.append(f'<rect x="{x:.2f}" y="{y}" width="{width * 9:.2f}" height="20" fill="{legend[i][1]}"><title>{esc(legend[i][0])}: {count:,}</title></rect>')
            x += width * 9
        fragments.append('</div></div>')
    svg.append('</svg>')
    fragments.extend(svg)
    if len(series) > len(shown):
        fragments.append(f'<p class="muted">Showing {len(shown)} of {len(series):,} rows. Complete values remain in the table and summary JSON.</p>')
    return ''.join(fragments)


def coverage_chart(matrices):
    full = [m for m in matrices if m['kind'] == 'full' and m.get('status') == 'evaluated']
    if not any(m.get('padding') is not None for m in full):
        return '<p class="empty">Unavailable: matching gene partitions and artificial missing-block tables are required.</p>'
    lookup = {m['group']: {r['taxon']: r for r in m['padding']['coverage']} if m.get('padding') is not None else {} for m in full}
    taxa = sorted({t for m in full for t in m['taxa']})
    shown = taxa[:50]
    width, height = 230 + len(full) * 150, 110 + 30 * len(shown)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" style="width:{width}px;max-width:none" viewBox="0 0 {width} {height}" font-family="sans-serif" role="img" aria-label="Taxon locus occupancy"><title>Taxon locus occupancy</title>']
    for i, matrix in enumerate(full):
        lines = group_label(matrix['group']).split(' · ')
        for line, label in enumerate(lines):
            svg.append(f'<text x="{240 + i * 150}" y="{16 + 15 * line}" font-size="11" fill="#334155">{esc(label)}</text>')
    for j, taxon in enumerate(shown):
        y = 82 + 30 * j
        svg.append(f'<text x="0" y="{y + 16}" font-size="12" fill="#334155">{esc(taxon[:30])}<title>{esc(taxon)}</title></text>')
        for i, matrix in enumerate(full):
            row = lookup[matrix["group"]].get(taxon)
            ratio = row["observed_loci"] / row["total_loci"] if row and row["total_loci"] else None
            color = f'hsl(182 60% {94 - 60 * ratio:.1f}%)' if ratio is not None else '#e2e8f0'
            text = f'{ratio:.0%}' if ratio is not None else 'N/R' if taxon in matrix['taxa'] else 'Absent'
            description = f'{taxon}; {matrix["group"]}; {row["observed_loci"]}/{row["total_loci"]} loci' if row else f'{taxon}: occupancy not recorded' if taxon in matrix['taxa'] else f'{taxon}: not present in this matrix'
            svg.append(f'<rect x="{230 + i * 150}" y="{y}" width="140" height="24" rx="3" fill="{color}"><title>{esc(description)}</title></rect>')
            svg.append(f'<text x="{280 + i * 150}" y="{y + 16}" font-size="11" fill="{"white" if ratio is not None and ratio > .7 else "#0f172a"}">{text}</text>')
    svg.append('</svg>')
    if len(taxa) > len(shown):
        svg.append(f'<p class="muted">First {len(shown)} of {len(taxa):,} taxa shown. All values are available in the summary JSON.</p>')
    return ''.join(svg)


def matrix_identity(matrix):
    """Keep raw/trimmed, sequence representations and selection routes explicit."""
    parts = matrix["group"].split("/")
    sequence = next((part for part in parts if part in {"codon", "prot", "nucl"}), matrix["seqtype"])
    mode = "Pseudogene" if "pseudogene" in parts else "Standard"
    selection = "Paired QC" if "paired" in parts else "Independent QC"
    return mode, selection, sequence, parts[-1]


def matrix_anchor(matrix):
    return "matrix-" + hashlib.sha256(matrix["matrix"].encode()).hexdigest()[:16]


def representation_name(sequence):
    return {"codon": "Codon (all positions)", "prot": "Protein", "nucl": "Nucleotide"}.get(sequence, sequence)


def derived_name(kind):
    return {"codon1st": "Codon position 1", "codon2nd": "Codon position 2",
            "codon3rd": "Codon position 3", "fourfold": "Fourfold-degenerate sites"}.get(kind, kind)


def matrix_catalog(matrices):
    primary = [matrix for matrix in matrices if matrix["kind"] == "full"]
    if not primary:
        return '<p class="empty">No primary concatenated matrices were found in the saved outputs.</p>'
    entries = {matrix_identity(matrix): matrix for matrix in primary}
    categories = {identity[:3] for identity in entries}
    # A codon output has a protein counterpart, even if that counterpart was not saved.
    categories |= {(mode, selection, "prot") for mode, selection, sequence in categories if sequence == "codon"}
    sequence_order = {"codon": 0, "prot": 1, "nucl": 2}
    ordered = sorted(categories, key=lambda item: (item[0], item[1], sequence_order.get(item[2], 3)))
    show_mode = len({category[0] for category in categories}) > 1 or any(category[0] == "Pseudogene" for category in categories)
    fragments = ['<div class="table-scroll matrix-catalog"><table><thead><tr>']
    for heading in (["Workflow"] if show_mode else []) + ["Selection", "Sequence type", "Raw", "Trimmed"]:
        fragments.append(f'<th scope="col">{heading}</th>')
    fragments.append('</tr></thead><tbody>')
    for mode, selection, sequence in ordered:
        fragments.append('<tr>')
        if show_mode:
            fragments.append(f'<td class="catalog-mode">{mode}</td>')
        fragments.append(f'<td class="catalog-selection">{selection}</td><th scope="row">{esc(representation_name(sequence))}</th>')
        for variant in ("raw", "trimmed"):
            entry = entries.get((mode, selection, sequence, variant))
            if entry is None:
                fragments.append(f'<td data-label="{variant.title()}"><span class="muted">Not found in saved outputs</span></td>')
                continue
            evaluated = entry.get("status") == "evaluated"
            label = "View matrix" if evaluated else "Not evaluated — view details"
            accessible = f'{mode} / {selection} / {representation_name(sequence)} / {variant.title()}: {label}'
            fragments.append(f'<td data-label="{variant.title()}"><a class="matrix-link" href="#{matrix_anchor(entry)}" aria-label="{esc(accessible)}">'
                             + ('<span aria-hidden="true">✓ </span>' if evaluated else '') + esc(label) + '</a>')
            if evaluated:
                fragments.append(f'<span class="matrix-meta">{esc(display(entry.get("loci")))} loci · '
                                 f'{len(entry["taxa"]):,} taxa · {entry["length"]:,} {esc(entry["unit"])}</span>')
            fragments.append('</td>')
        fragments.append('</tr>')
    fragments.append('</tbody></table></div>')
    return ''.join(fragments)


def matrix_details(matrix, matrices):
    mode, selection, sequence, variant = matrix_identity(matrix)
    title = f'{selection} · {representation_name(sequence)} · {variant.title()}'
    if mode == "Pseudogene":
        title = mode + " · " + title
    fragments = [f'<details class="matrix-detail" id="{matrix_anchor(matrix)}"><summary>{esc(title)}</summary>']
    fragments.append(f'<p class="muted">FASTA path relative to the output directory:</p><p><code>{esc("matrices/" + matrix["matrix"])}</code></p>')
    if matrix.get("status") != "evaluated":
        fragments.append(f'<p class="note warning">Not evaluated: {esc(matrix.get("error", "Statistics are unavailable."))}</p>')
    else:
        fragments.append(table(["Measure", "Value"], [
            ("Actually included loci", matrix.get("loci")), ("Taxa", len(matrix["taxa"])),
            ("Alignment length", f'{matrix["length"]:,} {matrix["unit"]}'),
            ("Complete codons", matrix.get("complete_codons") if sequence == "codon" else "Not applicable"),
            ("Skipped after QC", matrix.get("skipped_after_qc")),
            ("Artificial missing blocks", matrix["padding"]["blocks"] if matrix.get("padding") else None),
            ("Artificially filled characters", matrix["padding"]["characters"] if matrix.get("padding") else None)]))
    derived = [entry for entry in matrices if entry["group"] == matrix["group"] and entry["kind"] != "full"]
    if derived:
        fragments.append('<h3>Datasets derived from this codon matrix</h3>')
        fragments.append(table(["Dataset", "Source loci", "Length", "FASTA path relative to output directory"], [
            (derived_name(entry["kind"]), entry.get("loci"),
             f'{entry["length"]:,} {entry["unit"]}' if entry.get("status") == "evaluated" else "Not evaluated",
             "matrices/" + entry["matrix"]) for entry in derived]))
        fragments.append('<p class="muted">These are subsets of the parent codon matrix. Fourfold source-locus counts remain unavailable when no contribution map was recorded.</p>')
    fragments.append('</details>')
    return ''.join(fragments)


STYLE = """
:root{color-scheme:light;font-family:system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#172b3a;background:#f3f6f8;font-size:15px}
*{box-sizing:border-box}body{margin:0}main{max-width:1340px;margin:auto;padding:36px 30px}header{border-top:5px solid #197c80;padding:25px 0 18px}h1{font-size:36px;letter-spacing:-.8px;margin:8px 0 10px}h2{font-size:24px;margin:0 0 14px}h3{font-size:18px;margin:24px 0 10px}p{line-height:1.7;max-width:100ch}.lede{font-size:18px;color:#445c69;margin:8px 0 14px}a{color:#12676b}nav{display:flex;flex-wrap:wrap;gap:10px 20px;margin:23px 0 16px;font-size:14px}.eyebrow{color:#197c80;font-weight:700;letter-spacing:1.8px;font-size:12px}.muted,.empty{color:#526875;font-size:13px}.empty{padding:10px 0}.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px;margin:20px 0}.card,section{background:white;border:1px solid #dce5e9;border-radius:11px;padding:26px;margin-bottom:23px}.cards .card{margin:0;padding:19px 23px}.metric{font-size:31px;font-weight:700;margin-top:5px}.label{font-size:13px;color:#526875}.note{border-left:4px solid #197c80;background:#eef8f7;padding:13px 18px;font-size:14px;line-height:1.7;margin:18px 0}.warning{border-left-color:#b76620;background:#fff7ed}.note ul{margin:5px 0;padding-left:20px}.section-number{color:#197c80;font-size:11px;letter-spacing:1.4px;margin-bottom:7px}details{margin:13px 0}summary{cursor:pointer;font-weight:600;line-height:1.6}details[open]>summary{margin-bottom:14px}details.explanation{font-size:14px}.qc-panel{padding:18px;background:#fafcfc;border:1px solid #dce5e9;border-radius:8px}.chart-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:28px}.chart-grid h3{margin-top:6px}.chart-grid .desktop-chart{display:none}.chart-grid .mobile-bars{display:block}
.table-scroll{overflow:auto;max-height:550px;margin:12px 0 20px}table{border-collapse:collapse;width:100%;font-size:13px;line-height:1.5;font-variant-numeric:tabular-nums}th,td{text-align:left;border-bottom:1px solid #e5ecef;padding:11px 12px;vertical-align:top}th{position:sticky;top:0;background:#f1f6f7;white-space:nowrap}tbody tr:nth-child(even){background:#fbfcfd}td:first-child{max-width:330px;min-width:160px;overflow-wrap:anywhere}tr[hidden]{display:none!important}.table-tools,.filter-bar{display:flex;flex-wrap:wrap;gap:12px;align-items:end;margin:15px 0 8px}.control{display:flex;flex-direction:column;gap:6px;font-size:12px;font-weight:600;color:#445c69}.control input,.control select{font:inherit;font-size:13px;font-weight:400;color:#172b3a;background:white;border:1px solid #b7cbd0;border-radius:6px;padding:8px 10px;min-height:36px;max-width:100%}.table-tools .control:first-child{flex:1;min-width:150px}.table-tools .row-count{font-size:12px;color:#526875;padding-bottom:10px}.filter-bar{background:#f4f8f9;padding:15px;border-radius:8px}.filter-bar .control{flex:1;min-width:130px}.table-tools button{margin:0 0 1px}.header-sort{font:inherit;font-weight:600;text-align:left;border:0;background:transparent;padding:0;margin:0;color:#172b3a;cursor:pointer}.header-sort:hover{background:transparent;text-decoration:underline}th[aria-sort=ascending] .header-sort::after{content:' ↑'}th[aria-sort=descending] .header-sort::after{content:' ↓'}
.chart{overflow:auto;margin:16px 0}.chart svg{display:block;width:100%;max-width:100%;height:auto}.chart .export-chart{display:none}.mobile-bars{display:none}.bar-item{margin:15px 0 20px}.bar-caption{display:flex;justify-content:space-between;gap:14px;align-items:baseline;font-size:13px;line-height:1.6;margin-bottom:7px}.bar-caption strong:last-child{white-space:nowrap}.bar-caption>span:first-child{overflow-wrap:anywhere}.bar-track{height:13px;background:#e9f0f2;border-radius:4px;overflow:hidden}.bar-track>span{display:block;height:100%;background:#197c80;border-radius:4px}.stack-row{margin:18px 0 22px}.stack-row .bar-caption{flex-wrap:wrap;gap:3px 15px}.stack-row .bar-caption>span{color:#526875;font-size:12px}.stack-track{height:22px;background:#e9f0f2;border-radius:5px;display:flex;overflow:hidden}.stack-track>span{display:block;height:100%}.legend{display:flex;flex-wrap:wrap;gap:10px 20px;font-size:12px;color:#445c69;line-height:1.6;margin:12px 0}.legend>span{display:flex;align-items:center;gap:7px}.legend i{display:inline-block;width:13px;height:13px;border:1px solid #c6d3d8;border-radius:3px}.heatmap{max-height:620px}.chart.heatmap svg{max-width:none}.trim-card{border:1px solid #dce5e9;border-radius:8px;padding:17px 20px;margin:15px 0}.trim-card p{margin:7px 0 12px;font-size:13px}.trim-bar{display:grid;grid-template-columns:58px minmax(0,1fr) 120px;align-items:center;gap:10px;font-size:12px;margin:10px 0}.trim-bar strong{text-align:right}.trim-bar .bar-track{height:16px}
button{border:1px solid #b7cbd0;background:white;color:#12676b;border-radius:6px;padding:8px 12px;cursor:pointer;font:inherit;font-size:12px;margin:8px 8px 8px 0}button:hover{background:#eef8f7}button:focus-visible,a:focus-visible,summary:focus-visible,input:focus-visible,select:focus-visible{outline:2px solid #197c80;outline-offset:3px}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f6f8;padding:16px;border-radius:6px;font-size:12px}.matrix-meta{display:block;color:#526875;font-size:12px;margin-top:7px;line-height:1.6}.matrix-catalog{max-height:none;margin:18px 0}.matrix-catalog th[scope=row]{position:static;white-space:normal;background:white}.matrix-detail{scroll-margin-top:18px;padding:14px;border:1px solid #dce5e9;border-radius:7px}.matrix-detail:target{border-color:#197c80;background:#f5fbfa}.matrix-detail code{overflow-wrap:anywhere;font-size:12px}.matrix-link{font-weight:600}footer{font-size:12px;color:#526875;padding:10px 0 25px;line-height:1.7}.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
@media(max-width:850px){.chart-grid{grid-template-columns:1fr;gap:15px}}
@media(max-width:700px){main{padding:14px}h1{font-size:29px}h2{font-size:22px}section{padding:19px}.cards .card{padding:16px}.metric{font-size:25px;overflow-wrap:anywhere}.chart .desktop-chart{display:none}.mobile-bars{display:block}.control{width:100%}.table-tools .control:first-child{width:100%}.table-tools .row-count{padding-bottom:0}.table-scroll.data-table{overflow:visible;max-height:none}.data-table table,.data-table tbody{display:block}.data-table thead{display:none}.data-table tbody tr{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border:1px solid #dce5e9;border-radius:7px;margin:12px 0;padding:7px;background:white}.data-table td{border:0;min-width:0;padding:9px;font-size:13px;overflow-wrap:anywhere}.data-table td:first-child{grid-column:1/-1;max-width:none;font-weight:600;background:#f4f8f9;border-radius:4px}.data-table td::before{content:attr(data-label);display:block;font-size:11px;color:#526875;font-weight:500;margin-bottom:5px}.matrix-catalog{overflow:visible}.matrix-catalog table,.matrix-catalog tbody{display:block}.matrix-catalog thead{display:none}.matrix-catalog tbody tr{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border:1px solid #dce5e9;border-radius:7px;margin-bottom:14px;padding:4px;background:white}.matrix-catalog .catalog-mode,.matrix-catalog .catalog-selection,.matrix-catalog th[scope=row]{grid-column:1/-1;border:0;padding-bottom:4px}.matrix-catalog .catalog-selection{font-weight:600;color:#526875}.matrix-catalog th[scope=row]{padding-top:2px;padding-bottom:8px}.matrix-catalog td[data-label]{border:0;padding:10px;min-width:0}.matrix-catalog td[data-label]::before{content:attr(data-label);display:block;font-size:12px;color:#526875;font-weight:600;margin-bottom:8px}.matrix-catalog .matrix-meta{overflow-wrap:anywhere}.trim-card{padding:14px}.trim-bar{grid-template-columns:50px minmax(0,1fr) 95px;font-size:11px;gap:7px}}
@media print{body{background:white;print-color-adjust:exact;-webkit-print-color-adjust:exact}main{padding:0;max-width:none}nav,button:not(.header-sort),.table-tools,.filter-bar,#qc-branch{display:none!important}section{border:0;padding:16px 0;margin-bottom:18px}h2,h3,summary{break-after:avoid}.table-scroll{max-height:none;overflow:visible}th{position:static;white-space:normal}table{font-size:11px}th,td{padding:6px}td:first-child{min-width:0;max-width:none}.qc-panel[hidden]{display:block!important}.heatmap{max-height:none;overflow:visible}.heatmap svg{max-width:100%!important;width:100%!important}.trim-card,.stack-row{break-inside:avoid}.chart .export-chart{display:none}}
"""

REPORT_SCRIPT = r"""
function saveBlob(blob,name){const a=document.createElement('a');const u=URL.createObjectURL(blob);a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),60000);}
const report=JSON.parse(document.getElementById('report-data').textContent);
function csvCell(value){let text=String(value??'');if(/^[=+@]/.test(text)||(/^\-/.test(text)&&!/^\-\d[\d,.]*(%| nt| aa)?$/.test(text)))text="'"+text;return '"'+text.replaceAll('"','""')+'"';}
function downloadCsv(rows,name){saveBlob(new Blob(['\uFEFF'+rows.map(row=>row.map(csvCell).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'}),name);}
document.getElementById('download-summary').addEventListener('click',()=>saveBlob(new Blob([document.getElementById('report-data').textContent],{type:'application/json'}),'phyloprep-report-summary.json'));
document.querySelectorAll('time[datetime]').forEach(el=>{const date=new Date(el.dateTime);if(!Number.isNaN(date.getTime()))el.textContent=new Intl.DateTimeFormat('en-GB',{year:'numeric',month:'short',day:'2-digit',hour:'2-digit',minute:'2-digit',timeZoneName:'short'}).format(date);});
function openDetails(id){const target=document.getElementById(id);if(!target||target.tagName!=='DETAILS')return;let parent=target;while(parent){parent.open=true;parent=parent.parentElement?.closest('details');}target.hidden=false;if(target.classList.contains('qc-panel')){document.getElementById('qc-branch').value=id;showQc();}target.querySelector('summary').focus({preventScroll:true});target.scrollIntoView({block:'start'});}
document.querySelectorAll('a.matrix-link').forEach(a=>a.addEventListener('click',()=>openDetails(a.hash.slice(1))));window.addEventListener('hashchange',()=>openDetails(location.hash.slice(1)));
const qcSelect=document.getElementById('qc-branch');function showQc(){if(!qcSelect)return;document.querySelectorAll('.qc-panel').forEach(panel=>{panel.hidden=panel.id!==qcSelect.value;if(!panel.hidden)panel.open=true;});}qcSelect?.addEventListener('change',showQc);showQc();if(location.hash)openDetails(location.hash.slice(1));
const matrices=new Map(report.matrices.map(m=>[m.matrix,m]));const filters=document.getElementById('matrix-filters');const tableStates=[];
function matrixMatches(key){if(!key)return true;const m=matrices.get(key);if(!m)return false;const values=Object.fromEntries(Array.from(filters.querySelectorAll('select'),el=>[el.dataset.field,el.value]));const parts=m.group.split('/');return(values.kind==='all'||(values.kind==='full'?m.kind==='full':m.kind!=='full'))&&(values.selection==='all'||values.selection===(parts.includes('paired')?'paired':'independent'))&&(values.representation==='all'||parts.includes(values.representation))&&(values.variant==='all'||parts.at(-1)===values.variant);}
function applyRows(state){const query=state.search.value.trim().toLowerCase();let visible=0;for(const row of state.rows){row.hidden=!matrixMatches(row.dataset.matrixKey)||!row.textContent.toLowerCase().includes(query);if(!row.hidden)visible++;}state.count.textContent=visible+' of '+state.rows.length+' rows';}
function sortValue(value){const text=value.trim();if(/^-?\d[\d,]*(\.\d+)?(%| nt| aa)?$/.test(text))return Number.parseFloat(text.replaceAll(',',''));return text.toLowerCase();}
function sortRows(state,index,direction){const ascending=direction==='asc';const positions=new Map(state.rows.map((row,i)=>[row,i]));const ordered=[...state.rows].sort((a,b)=>{if(index<0)return positions.get(a)-positions.get(b);const av=sortValue(a.cells[index].textContent),bv=sortValue(b.cells[index].textContent);const unknown=v=>typeof v==='string'&&/^(unavailable|not applicable|not recorded)/.test(v);if(unknown(av)!==unknown(bv))return unknown(av)?1:-1;const comparison=typeof av==='number'&&typeof bv==='number'?av-bv:String(av).localeCompare(String(bv),undefined,{numeric:true});return comparison*(ascending?1:-1);});ordered.forEach(row=>state.body.appendChild(row));state.headers.forEach((h,i)=>h.setAttribute('aria-sort',i===index?(ascending?'ascending':'descending'):'none'));state.sort.value=index<0?'original':index+':'+direction;}
document.querySelectorAll('.data-table').forEach((container,index)=>{const table=container.querySelector('table'),body=table.tBodies[0],rows=Array.from(body.rows),headers=Array.from(table.tHead.rows[0].cells),title=container.dataset.title;const tools=document.createElement('div');tools.className='table-tools';const searchLabel=document.createElement('label');searchLabel.className='control';searchLabel.textContent='Search table';const search=document.createElement('input');search.type='search';search.placeholder='Search rows…';search.setAttribute('aria-label','Search '+title);searchLabel.appendChild(search);const sortLabel=document.createElement('label');sortLabel.className='control';sortLabel.textContent='Sort rows';const sort=document.createElement('select');sort.setAttribute('aria-label','Sort '+title);sort.add(new Option('Original order','original'));headers.forEach((header,i)=>{for(const direction of ['asc','desc'])sort.add(new Option(header.textContent+' '+(direction==='asc'?'↑':'↓'),i+':'+direction));});sortLabel.appendChild(sort);const button=document.createElement('button');button.textContent='Download visible CSV';const count=document.createElement('span');count.className='row-count';count.setAttribute('role','status');tools.append(searchLabel,sortLabel,button,count);container.before(tools);const state={container,body,rows,headers,search,sort,count};tableStates.push(state);search.addEventListener('input',()=>applyRows(state));sort.addEventListener('change',()=>{if(sort.value==='original')sortRows(state,-1,'asc');else{const [column,direction]=sort.value.split(':');sortRows(state,Number(column),direction);}});headers.forEach((header,i)=>{const label=header.textContent;const control=document.createElement('button');control.className='header-sort';control.textContent=label;control.addEventListener('click',()=>sortRows(state,i,header.getAttribute('aria-sort')==='ascending'?'desc':'asc'));header.textContent='';header.appendChild(control);});button.addEventListener('click',()=>{const visible=Array.from(body.rows).filter(row=>!row.hidden);const withPaths=visible.some(row=>row.dataset.matrixKey);const exportRows=[headers.map(h=>h.textContent).concat(withPaths?['FASTA path']:[]),...visible.map(row=>Array.from(row.cells,c=>c.textContent).concat(withPaths?['matrices/'+(row.dataset.matrixKey||'')]:[]))];downloadCsv(exportRows,'phyloprep-'+(title||String(index)).toLowerCase().replace(/[^a-z0-9]+/g,'-').slice(0,90)+'.csv');});applyRows(state);});
const scope=document.createElement('p');scope.id='matrix-scope';scope.className='muted';scope.setAttribute('role','status');filters.after(scope);function refreshMatrixFilters(){tableStates.forEach(applyRows);scope.textContent='Showing: '+Array.from(filters.querySelectorAll('select'),el=>el.selectedOptions[0].textContent).join(' · ')+'.';}filters.addEventListener('change',refreshMatrixFilters);refreshMatrixFilters();
document.querySelectorAll('.chart svg').forEach((svg,i)=>{const b=document.createElement('button');b.textContent='Save chart as SVG';b.addEventListener('click',()=>{const copy=svg.cloneNode(true);copy.removeAttribute('class');copy.removeAttribute('style');saveBlob(new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'}),'phyloprep-chart-'+(i+1)+'.svg');});svg.parentElement.appendChild(b);});
document.getElementById('download-coverage').addEventListener('click',()=>{const rows=[['Matrix','Taxon','Observed loci','Included loci','Occupancy (%)','Status']];for(const m of report.matrices.filter(m=>m.kind==='full'&&m.status==='evaluated')){const coverage=new Map((m.padding?.coverage||[]).map(row=>[row.taxon,row]));for(const taxon of m.taxa){const row=coverage.get(taxon);rows.push([m.group,taxon,row?.observed_loci??'',row?.total_loci??'',row&&row.total_loci?(100*row.observed_loci/row.total_loci).toFixed(2):'',row?'Recorded':'Not recorded']);}}downloadCsv(rows,'phyloprep-taxon-occupancy.csv');});
let printState=[];window.addEventListener('beforeprint',()=>{printState=Array.from(document.querySelectorAll('details'),element=>({element,open:element.open,hidden:element.hidden}));printState.forEach(({element})=>{element.open=true;element.hidden=false;});});window.addEventListener('afterprint',()=>{printState.forEach(({element,open,hidden})=>{element.open=open;element.hidden=hidden;});printState=[];});
"""


def render_html(report):
    overview = report['overview']
    matrices = report['matrices']
    full = [m for m in matrices if m['kind'] == 'full']
    evaluated = [m for m in matrices if m.get('status') == 'evaluated']
    generated = report['generated_at']
    fragments = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
                 '<meta name="viewport" content="width=device-width,initial-scale=1">',
                 '<title>PhyloPrep Analysis Report</title><style>' + STYLE + '</style></head><body><main>',
                 '<header><div class="eyebrow">PHYLOPREP · RESULTS &amp; PROVENANCE</div><h1>PhyloPrep Analysis Report</h1>',
                 '<p class="lede">Review locus retention, matrix composition and sequence changes.</p>',
                 f'<p class="muted">Run: <strong>{esc(report["run_name"])}</strong><br>Generated <time datetime="{esc(generated)}">{esc(readable_time(generated))}</time> · Report v{__version__}</p>',
                 '<nav aria-label="Report sections"><a href="#matrix-catalog">Matrix catalog</a><a href="#overview">Processing</a><a href="#qc">QC</a><a href="#matrices">Statistics &amp; trimming</a><a href="#events">Sequence changes</a><a href="#coverage">Coverage</a><a href="#methods">Definitions &amp; provenance</a></nav>',
                 '<button id="download-summary">Download summary JSON</button><button onclick="window.print()">Print / save PDF</button></header>']
    if report['warnings']:
        fragments.append('<aside class="note warning"><strong>Data notes</strong><ul>' + ''.join(f'<li>{esc(w)}</li>' for w in report['warnings'][:3]) + '</ul>'
                         + ('<a href="#availability">View all data notes</a>' if len(report['warnings']) > 3 else '') + '</aside>')
    fragments.append('<div class="cards">')
    for label, value in (('Input loci', overview['input_loci']), ('Aligned loci', overview['aligned_loci'])):
        fragments.append(f'<div class="card"><div class="label">{label}</div><div class="metric">{esc(display(value))}</div></div>')
    fragments.append('</div><section id="matrix-catalog"><h2>Which matrices were produced?</h2><p>Select a matrix to open its locus count, taxon count, length and file path. Codon-position and fourfold datasets are listed with their parent codon matrix.</p>')
    fragments.append(matrix_catalog(matrices))
    fragments.append('<details class="explanation"><summary>How to read this catalog</summary><p><strong>Independent QC</strong> selects codon and protein loci separately. <strong>Paired QC</strong> retains the same locus set for both. The resulting matrices can be identical when their retained sets match.</p><p><strong>Codon (all positions)</strong> includes all three positions of the retained codons; it does not mean all input loci passed QC. FASTA, PHYLIP and NEXUS are formats of the same dataset.</p></details></section>')
    fragments.append('<section id="overview"><div class="section-number">01 / PROCESSING</div><h2>From input loci to retained datasets</h2><p>QC and concatenation make separate decisions for each branch. Counts across representations or branches must not be added to estimate unique input loci.</p>')
    overview_rows = [('Input loci', overview['input_loci']), ('Successfully aligned loci', overview['aligned_loci'])]
    unavailable = []
    for key, label in (('failed_alignments', 'Recorded alignment failures'), ('reused_loci', 'Reused aligned loci')):
        if overview[key] is None:
            unavailable.append((label, 'Not recorded in the available run records'))
        else:
            overview_rows.append((label, overview[key]))
    fragments.append(table(['Processing measure', 'Loci'], overview_rows, 'Locus processing'))
    fragments.append('<p class="muted">QC pass counts and actual concatenation inclusion are shown separately below. Input minus aligned is not assumed to equal failed alignments.</p></section>')
    fragments.append('<section id="qc"><div class="section-number">02 / QUALITY CONTROL</div><h2>Which loci passed QC?</h2>')
    legend = [('Passed', '#197c80'), ('Failed', '#b76620'), ('Not evaluated', '#b5c2cc')]
    qcs = report['qc']
    def qc_label(q):
        label = group_label(q['branch'])
        return label.replace('Independent QC', 'Recorded QC') if q['selection'] == 'Recorded QC' else label
    fragments.append('<div class="chart">' + stacked_chart('QC outcomes by branch', [(qc_label(q), [q['passed'], q['failed'], q['not_evaluated']]) for q in qcs], legend, limit=30) + '</div>')
    fragments.append('<details><summary>Exact QC counts</summary>' + table(['QC branch', 'Passed', 'Unique failed loci', 'Not evaluated'], [(qc_label(q), q['passed'], q['failed'], q['not_evaluated']) for q in qcs], 'QC outcomes') + '</details>')
    fragments.append('<p class="note"><strong>One failed locus can have several reasons.</strong> Reason hits overlap, so their sum can exceed the number of failed loci and their percentages can exceed 100%. Complete failure combinations count each failed locus once. No primary reason is assigned.</p>')
    if qcs:
        fragments.append('<label class="control">Inspect a QC branch<select id="qc-branch">' + ''.join(f'<option value="qc-panel-{i}">{esc(qc_label(q))}</option>' for i, q in enumerate(qcs)) + '</select></label>')
    else:
        fragments.append('<p class="empty">QC records are unavailable.</p>')
    for index, q in enumerate(qcs):
        fragments.append(f'<details class="qc-panel" id="qc-panel-{index}"' + (' open' if index == 0 else '') + f'><summary>{esc(qc_label(q))} — {q["failed"]:,} unique failed loci</summary>')
        if not q['failed']:
            fragments.append('<p class="note">No failures were recorded in this branch.</p>')
        else:
            combinations = Counter(q['failure_combinations'])
            top = dict(combinations.most_common(10))
            rest = sum(combinations.values()) - sum(top.values())
            if rest:
                top['Other complete reason combinations'] = rest
            fragments.append('<div class="chart-grid"><div><h3>Overlapping reason hits</h3><p class="muted">A locus may appear in several bars. Percentages use unique failed loci as the denominator.</p><div class="chart">' + bar_chart('Overlapping QC reason hits', q['reason_hits'], q['failed'], labeler=reason_label) + '</div></div>')
            fragments.append('<div><h3>Complete failure combinations</h3><p class="muted">Each locus appears once. Counts sum to the unique failed-locus total.</p><div class="chart">' + bar_chart('Complete QC failure combinations', top, q['failed'], labeler=combination_label) + '</div></div></div>')
            fragments.append('<details><summary>Exact reasons, combinations and original codes</summary>')
            fragments.append(table(['Reason', 'Original code', 'Loci with reason', '% of failed loci'], [(reason_label(k), k, v, pct(v / q['failed'])) for k, v in sorted(q['reason_hits'].items())], qc_label(q) + ' reason hits'))
            fragments.append(table(['Complete combination', 'Original codes', 'Unique loci'], [(combination_label(k), k, v) for k, v in sorted(combinations.items())], qc_label(q) + ' failure combinations'))
            fragments.append('</details>')
        fragments.append('</details>')
    fragments.append('<p class="muted">Paired QC reasons summarize codon/protein QC outcomes and paired PIS decisions. Inspect the representation branches for individual conditions. Recorded QC tables may precede independent PIS filtering in older outputs.</p></section>')
    fragments.append('<section id="matrices"><div class="section-number">03 / MATRIX STATISTICS</div><h2>What is in each final matrix?</h2><p>Values are recomputed from the final FASTA files. The filters below apply to all three tables. Start with primary matrices, then select derived datasets when needed.</p>')
    fragments.append('<div class="filter-bar" id="matrix-filters"><label class="control">Datasets<select data-field="kind"><option value="full">Primary matrices</option><option value="derived">Derived datasets</option><option value="all">All datasets</option></select></label><label class="control">Selection<select data-field="selection"><option value="all">All selections</option><option value="independent">Independent QC</option><option value="paired">Paired QC</option></select></label><label class="control">Representation<select data-field="representation"><option value="all">All representations</option><option value="codon">Codon</option><option value="prot">Protein</option><option value="nucl">Nucleotide</option></select></label><label class="control">Variant<select data-field="variant"><option value="all">Raw and trimmed</option><option value="raw">Raw</option><option value="trimmed">Trimmed</option></select></label></div>')
    keys = [m['matrix'] for m in evaluated]
    fragments.append('<h3>Matrix composition</h3>' + table(['Matrix', 'Source loci', 'Taxa', 'Length', 'Complete codons'], [(matrix_label(m), m.get('loci'), len(m['taxa']), f'{m["length"]:,} {m["unit"]}', m['complete_codons'] if m['seqtype'] == 'codon' else 'Not applicable') for m in evaluated], 'Matrix composition', keys))
    fragments.append('<h3>Variation and phylogenetic information</h3><p class="muted">PIS = parsimony-informative sites. Only resolved character states contribute to these site counts.</p>' + table(['Matrix', 'Variable sites', 'PIS', 'PIS (%)', 'Singleton sites', 'Constant sites', 'No resolved states'], [(matrix_label(m), m['variable'], m['pis'], pct(m['pis'] / m['length']), m['singleton'], m['constant'], m['unresolved_only']) for m in evaluated], 'Matrix variation', keys))
    fragments.append('<h3>Gaps, missing data and ambiguity</h3><p class="muted"><strong>−</strong> = alignment gap; <strong>?</strong> = missing character; <strong>N/X</strong> = unresolved nucleotide/protein character. Only missing-block records identify artificial padding. GC uses resolved A/C/G/T; it does not apply to proteins. These are final character-cell percentages; codon QC ambiguity checks use codon-level units.</p>' + table(['Matrix', 'Gap (%)', 'Missing ? (%)', 'Ambiguous (%)', 'GC (%)'], [(matrix_label(m), pct(m['gap_ratio']), pct(m['missing_ratio']), pct(m['ambiguous_ratio']), 'Not applicable' if m['gc_ratio'] is None else pct(m['gc_ratio'])) for m in evaluated], 'Matrix missingness and composition', keys))
    fragments.append('<details><summary>Matrix files, inclusion details and derived datasets</summary>')
    for m in full:
        fragments.append(matrix_details(m, matrices))
    fragments.append('</details><h3>What changed after trimming?</h3><p>Length comparisons use only loci included in both raw and trimmed matrices of the same representation and selection. The retained locus sets of the final matrices can still differ.</p>')
    comparisons = report['trimming_comparisons']
    for c in comparisons:
        raw, trimmed = c['raw_common_columns'], c['trimmed_common_columns']
        fragments.append(f'<div class="trim-card"><strong>{esc(group_label(c["raw"].rsplit("/", 1)[0]))}</strong>')
        if raw is None or trimmed is None:
            fragments.append('<p class="empty">Unavailable: matching gene partition records are needed.</p>')
        elif not c['common_loci']:
            fragments.append('<p class="empty">No shared retained loci; a matched length comparison is not applicable.</p>')
        else:
            fragments.append(f'<p>{c["common_loci"]:,} shared loci · {c.get("raw_only", 0):,} raw-only · {c.get("trimmed_only", 0):,} trimmed-only</p>')
            maximum = max(raw, trimmed, 1)
            for label, length, color in [('Raw', raw, '#91b9bb'), ('Trimmed', trimmed, '#197c80')]:
                fragments.append(f'<div class="trim-bar"><span>{label}</span><div class="bar-track"><span style="width:{100 * length / maximum:.3f}%;background:{color}"></span></div><strong>{length:,} {esc(c.get("unit", "columns"))}</strong></div>')
            text = 'No columns removed from the shared locus set.' if raw == trimmed else f'{pct((raw - trimmed) / raw)} fewer columns in the shared locus set.' if raw else 'No raw columns were available for a percentage comparison.'
            fragments.append(f'<p class="muted">{esc(text)}</p>')
        fragments.append('</div>')
    if not comparisons:
        fragments.append('<p class="empty">No matching raw/trimmed matrix pairs were available.</p>')
    fragments.append('<details><summary>Exact matched-locus comparison values</summary>' + table(['Representation / selection', 'Shared loci', 'Raw-only loci', 'Trimmed-only loci', 'Raw shared-set columns', 'Trimmed shared-set columns', 'Unit'], [(group_label(c['raw'].rsplit('/', 1)[0]), c['common_loci'], c.get('raw_only'), c.get('trimmed_only'), c['raw_common_columns'], c['trimmed_common_columns'], c.get('unit')) for c in comparisons], 'Matched trimming comparison') + '</details><p class="muted">Bars are scaled within each comparison; nucleotide and amino-acid lengths have different units. Matched loci do not necessarily imply identical taxon sets. Paired refers to codon/protein pairing within a variant, not raw/trimmed pairing.</p></section>')
    fragments.append('<section id="events"><div class="section-number">04 / SEQUENCE CHANGES</div><h2>Which sequence changes were recorded?</h2>')
    audit = report.get('audit')
    if audit is None:
        fragments.append('<p class="empty">Unavailable: reports/sequence_changes.tsv was not found. Per-locus copies are not added together because that could double-count events.</p>')
    else:
        fragments.append(f'<p><strong>{audit["events"]:,} distinct recorded events</strong>, counted once from the canonical audit report. Routine normalization and cleanup are included; an event is not itself a QC failure. Exact duplicate rows ignored: {audit["duplicates_ignored"]:,}.</p>')
        if 'by_taxon_event' in audit:
            categories = ['Terminal cleanup', 'Symbol normalization', 'Internal stops', 'Frame-related events', 'Identifier changes', 'Other recorded events']
            palette = ['#91b9bb', '#197c80', '#b76620', '#725899', '#4b77a1', '#9cabb3']
            counts = {}
            for entry in audit['by_taxon_event']:
                values = counts.setdefault(entry['taxon'], [0] * len(categories))
                values[categories.index(event_category(entry['event']))] += entry['count']
            active = [i for i in range(len(categories)) if any(values[i] for values in counts.values())]
            series = [(taxon['taxon'], [counts.get(taxon['taxon'], [0] * len(categories))[i] for i in active]) for taxon in audit['taxa']]
            fragments.append('<div class="chart">' + stacked_chart('Sequence-change events by taxon and category', series, [(categories[i], palette[i]) for i in active]) + '</div><p class="muted">Bar segments show each taxon\'s event composition; exact counts are printed above the bars.</p>')
        else:
            fragments.append('<div class="chart">' + bar_chart('Sequence-change events by taxon', {t['taxon']: t['events'] for t in audit['taxa']}) + '</div>')
            fragments.append('<p class="muted">Taxon-by-event categories were not recorded in this older summary. The recorded event-type totals are shown separately.</p><details><summary>Recorded event-type totals</summary><div class="chart">' + bar_chart('Recorded event types', {entry['event'] + ' / ' + entry['action']: entry['count'] for entry in audit['types']}, labeler=lambda key: event_label(key.split(' / ')[0]) + ' / ' + key.split(' / ')[1]) + '</div></details>')
        fragments.append(table(['Taxon / identifier', 'Events', 'Affected sequences', 'NT character-event spans', 'AA character-event spans', 'Events without coordinates'], [(t['taxon'], t['events'], t['affected_sequences'] if t['taxon'] != 'Unassigned' else None, t['nt_characters'], t['aa_characters'], t['events_without_coordinates']) for t in audit['taxa']], 'Sequence changes by taxon'))
        fragments.append('<details><summary>Event types, actions and original codes</summary>' + table(['Event description', 'Original code', 'Action', 'Events'], [(event_label(t['event']), t['event'], t['action'], t['count']) for t in audit['types']], 'Sequence event types') + '</details><p class="muted">Taxa use mapped IDs when recorded; otherwise original sequence IDs are shown. Character-event spans sum source-coordinate lengths and may overlap. They are not unique edited positions. NT and AA spans are separate units. Event rates relative to original sequence length are not inferred.</p>')
    fragments.append('</section><section id="coverage"><div class="section-number">05 / COVERAGE</div><h2>How complete is taxon coverage?</h2><p>Occupancy is the fraction of included loci for which a taxon was not artificially padded. An observed sequence containing gaps or ambiguity still counts as an observed locus.</p><div class="legend"><span><i style="background:hsl(182 60% 94%)"></i>0%</span><span><i style="background:hsl(182 60% 64%)"></i>50%</span><span><i style="background:hsl(182 60% 34%)"></i>100%</span><span><i style="background:#e2e8f0"></i>N/R = not recorded; Absent = taxon not in matrix</span></div><div class="chart heatmap">' + coverage_chart(matrices) + '</div><button id="download-coverage">Download occupancy CSV (all taxa)</button></section>')
    fragments.append('<section id="methods"><div class="section-number">06 / DEFINITIONS &amp; PROVENANCE</div><h2>How to interpret these values</h2>')
    definitions = [('Locus', 'An input alignment unit, often a gene or gene fragment.'), ('Resolved states', 'A/C/G/T for nucleotides; the 20 standard amino acids for proteins.'), ('Variable site', 'A column containing at least two resolved states.'), ('Parsimony-informative site (PIS)', 'A column with at least two resolved states, each observed at least twice.'), ('Constant / no resolved states', 'Constant means exactly one resolved state. Columns with no resolved states are counted separately.'), ('Singleton site', 'Exactly two resolved states, with at least one observed once.'), ('Gap / missing / ambiguity percentages', 'Fractions of all taxon-by-column cells. Gaps and unresolved symbols do not contribute resolved states.'), ('GC (%)', 'G+C divided by resolved A+C+G+T. Not applicable to proteins or matrices with no resolved bases.'), ('Unavailable', 'A required record is absent or could not be validated. This never means zero.'), ('Fourfold source-locus count', 'Unavailable unless a contributing-locus map was recorded; the parent matrix count is not substituted.')]
    fragments.append('<details><summary>Metric definitions</summary>' + table(['Metric', 'Definition'], definitions, 'Metric definitions') + '</details>')
    fragments.append('<details id="availability"><summary>Data availability and interpretation limits</summary><p>Header-only missing-block tables establish zero artificial blocks. Missing tables do not. Original absolute input paths are not followed. Recorded checksums are not reverified by this report.</p>' + table(['Optional record', 'Availability'], unavailable, 'Optional records') + ''.join(f'<p class="note warning">{esc(w)}</p>' for w in report['warnings']) + '<p>These metrics describe processing and filtering decisions. They do not establish orthology, biological alignment accuracy or suitability for every downstream method.</p></details>')
    if report['provenance_source']:
        prov = report['provenance']
        params = prov.get('parameters') or {}
        friendly = {'genetic_code': 'NCBI genetic-code table', 'seqtype': 'Sequence mode', 'align_software': 'Aligner', 'min_length': 'Minimum alignment length', 'min_nucl_pis': 'Minimum nucleotide PIS', 'min_protein_pis': 'Minimum protein PIS', 'paired_pis_policy': 'Paired PIS rule', 'missing_taxa': 'Missing-taxon policy', 'fourfold_missing_policy': 'Fourfold missing-codon policy', 'G': 'Trimming: maximum gap proportion', 'N': 'Trimming: maximum N proportion', 'X': 'Trimming: maximum X proportion', 'notrim': 'Skip trimming', 'allow_internal_stop': 'Allow internal stops'}
        fragments.append('<details><summary>Recorded run settings</summary>' + table(['Setting', 'Recorded value / effective default'], [(friendly.get(k, k.replace('_', ' ').capitalize()), setting_value(k, v, params)) for k, v in params.items()], 'Run settings') + '<details><summary>Full recorded provenance (JSON)</summary><pre>' + esc(json.dumps(prov, indent=2, ensure_ascii=False)) + '</pre></details></details>')
    else:
        fragments.append('<p class="empty">A portable run manifest was not available. Run parameters, versions, checksums and execution time are not inferred. Report generation time is not analysis runtime.</p>')
    payload = json.dumps(report, ensure_ascii=True, separators=(',', ':')).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    fragments.append('</section><footer>Generated locally with Python standard-library tools. This self-contained report makes no network requests. SVG charts and visible table rows can be downloaded.</footer>')
    fragments.append(f'<script type="application/json" id="report-data">{payload}</script><script>' + REPORT_SCRIPT + '</script></main></body></html>')
    return ''.join(fragments)


def write_report(root, run_metadata=None):
    """Atomically publish a report; raise errors for the caller to handle."""
    root = Path(root).resolve()
    report = build_report(root, run_metadata=run_metadata)
    output = root / "phyloprep-report.html"
    content = render_html(report)
    fd, temporary = tempfile.mkstemp(prefix=".phyloprep-report-", suffix=".html", dir=root)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
        os.replace(temporary, output)
    finally:
        Path(temporary).unlink(missing_ok=True)
    print(f"[PhyloPrep report] Saved {output}", flush=True)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path, help="Final PhyloPrep output directory")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    try:
        write_report(args.input)
        return 0
    except (ValueError, OSError, KeyError, csv.Error, sqlite3.Error) as error:
        print(f"PhyloPrep report: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
