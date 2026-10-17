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

__version__ = "1.0.0"
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


def qc_summary(root, warnings):
    result = []
    for path in visible_files(root / "qc"):
        paired = path.name == "loci.details.tsv"
        if path.name != "alignment.details.tsv" and not paired:
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
    return dict(source="reports/sequence_changes.tsv", duplicates_ignored=duplicates,
                events=sum(t["events"] for t in taxa), taxa=taxa, types=kinds)


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


def build_report(root):
    root = Path(root).resolve()
    if not root.is_dir() or not any((root / name).is_dir() for name in ("alignments", "qc", "matrices", "reports")):
        raise ValueError("Input must be a PhyloPrep final output directory")
    warnings = []
    manifest_path = root / "phyloprep-run.json"
    manifest = read_json(manifest_path) if manifest_path.is_file() else {}
    if not isinstance(manifest, dict):
        raise ValueError("The run manifest must contain a JSON object")
    report = dict(schema="phyloprep-report-v1", report_version=__version__,
                  generated_at=datetime.now(timezone.utc).isoformat(), run_name=root.name,
                  overview=processing_overview(root, manifest), qc=qc_summary(root, warnings),
                  matrices=[], trimming_comparisons=[], warnings=warnings,
                  provenance={key: manifest[key] for key in (
                      "created_at", "module", "script_version", "source_sha256", "parameters", "equivalent_cli",
                      "runtime", "reused_loci", "matrix_padding", "engine_metadata") if key in manifest},
                  provenance_source="phyloprep-run.json" if manifest else None)
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
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def pct(value):
    return UNAVAILABLE if value is None else f"{100 * value:.2f}%"


def table(headers, data):
    return '<div class="table-scroll"><table><thead><tr>' + ''.join(
        f'<th scope="col">{esc(h)}</th>' for h in headers) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join(f'<td>{esc(display(v))}</td>' for v in row) + '</tr>' for row in data
    ) + '</tbody></table></div>'


def bar_chart(title, values, denominator=None, limit=12):
    ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
    if not ordered:
        return '<p class="empty">No recorded observations.</p>'
    omitted = max(0, len(ordered) - limit)
    ordered = ordered[:limit]
    width, height = 900, 48 + 48 * len(ordered)
    maximum = max(1, max(v for _, v in ordered))
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" font-family="sans-serif" role="img" aria-label="{esc(title)}">',
             f'<title>{esc(title)}</title>']
    for i, (label, value) in enumerate(ordered):
        y = 24 + i * 48
        short = label if len(label) <= 48 else label[:45] + "…"
        suffix = f" ({value / denominator:.1%})" if denominator else ""
        parts.extend([f'<text x="0" y="{y + 14}" font-size="12" fill="#334155">{esc(short)}</text>',
                      f'<rect x="365" y="{y}" width="{max(2, value / maximum * 400):.1f}" height="24" rx="4" fill="#197c80"><title>{esc(label)}: {value:,}{esc(suffix)}</title></rect>',
                      f'<text x="775" y="{y + 17}" font-size="12" fill="#0f172a">{value:,}{esc(suffix)}</text>'])
    parts.append('</svg>')
    if omitted:
        parts.append(f'<p class="muted">Top {limit} categories shown; {omitted} additional categories are included in the tables or summary JSON.</p>')
    return ''.join(parts)


def coverage_chart(matrices):
    full = [m for m in matrices if m.get("padding") is not None]
    if not full:
        return '<p class="empty">Unavailable: matching gene partitions and artificial missing-block tables are required.</p>'
    lookup = {m["group"]: {r["taxon"]: r for r in m["padding"]["coverage"]} for m in full}
    taxa = sorted({t for group in lookup.values() for t in group})
    shown = taxa[:50]
    width, height = 300 + len(full) * 140, 100 + 28 * len(shown)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" font-family="sans-serif" role="img" aria-label="Taxon locus occupancy"><title>Taxon locus occupancy</title>']
    for i, matrix in enumerate(full):
        svg.append(f'<text x="{310 + i * 140}" y="18" font-size="11" fill="#334155">{esc(matrix["group"][:23])}<title>{esc(matrix["group"])}</title></text>')
    for j, taxon in enumerate(shown):
        y = 45 + 28 * j
        svg.append(f'<text x="0" y="{y + 16}" font-size="12" fill="#334155">{esc(taxon[:43])}</text>')
        for i, matrix in enumerate(full):
            row = lookup[matrix["group"]].get(taxon)
            ratio = row["observed_loci"] / row["total_loci"] if row and row["total_loci"] else None
            color = f'hsl(182 60% {94 - 60 * ratio:.1f}%)' if ratio is not None else '#e2e8f0'
            text = f'{ratio:.0%}' if ratio is not None else 'N/A'
            description = f'{taxon}; {matrix["group"]}; {row["observed_loci"]}/{row["total_loci"]} loci' if row else f'{taxon}: not present in this matrix'
            svg.append(f'<rect x="{300 + i * 140}" y="{y}" width="130" height="23" rx="3" fill="{color}"><title>{esc(description)}</title></rect>')
            svg.append(f'<text x="{348 + i * 140}" y="{y + 16}" font-size="11" fill="{"white" if ratio is not None and ratio > .7 else "#0f172a"}">{text}</text>')
    svg.append('</svg>')
    if len(taxa) > len(shown):
        svg.append(f'<p class="muted">First {len(shown)} of {len(taxa):,} taxa shown. All values are available in the summary JSON.</p>')
    return ''.join(svg)


STYLE = """
:root{color-scheme:light;font-family:system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#172b3a;background:#f3f6f8}
*{box-sizing:border-box}body{margin:0}main{max-width:1420px;margin:auto;padding:40px 36px}header{border-top:5px solid #197c80;padding:26px 0 18px}
h1{font-size:38px;letter-spacing:-1px;margin:5px 0}h2{font-size:24px;margin:0 0 15px}h3{font-size:17px;margin:22px 0 8px}
p{line-height:1.65}a{color:#12676b}nav{display:flex;flex-wrap:wrap;gap:18px;margin:18px 0}.eyebrow{color:#197c80;font-weight:700;letter-spacing:2px;font-size:12px}
.muted,.empty{color:#566a79;font-size:14px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin:24px 0}.card,section{background:white;border:1px solid #dde5ea;border-radius:10px;padding:24px;margin-bottom:22px}
.cards .card{margin:0;padding:20px}.metric{font-size:29px;font-weight:700;margin-top:8px}.label{font-size:13px;color:#566a79}.note{border-left:4px solid #197c80;background:#eef8f7;padding:12px 18px;font-size:14px}.warning{border-left-color:#bd7125;background:#fff7ed}
.table-scroll{overflow:auto;max-height:530px}table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}th,td{text-align:left;border-bottom:1px solid #e6edf0;padding:11px 13px;vertical-align:top}th{position:sticky;top:0;background:#f1f6f7;white-space:nowrap}td:first-child{max-width:350px;overflow-wrap:anywhere}
.chart{overflow:auto;margin:16px 0}.chart svg{display:block;min-width:720px;max-width:100%;height:auto}button{border:1px solid #bbced0;background:#fff;color:#12676b;border-radius:6px;padding:8px 12px;cursor:pointer;font:inherit;font-size:12px;margin:8px 8px 8px 0}button:hover{background:#eef8f7}
details{margin:12px 0}summary{cursor:pointer;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f6f8;padding:16px;border-radius:6px;font-size:12px}.section-number{color:#197c80;font-size:12px;letter-spacing:1px}footer{font-size:12px;color:#566a79;padding:10px 0 25px}
@media(max-width:700px){main{padding:15px}h1{font-size:29px}section{padding:18px}.metric{font-size:24px}}
@media print{body{background:white}main{padding:0;max-width:none}nav,button{display:none}section{break-inside:avoid;border:0;padding:12px 0}.table-scroll{max-height:none;overflow:visible}th{position:static}.chart svg{min-width:0}}
"""


def render_html(report):
    overview = report["overview"]
    matrices = report["matrices"]
    full = [m for m in matrices if m["kind"] == "full"]
    evaluated = [m for m in matrices if m.get("status") == "evaluated"]
    fragments = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
                 '<meta name="viewport" content="width=device-width,initial-scale=1">',
                 '<title>PhyloPrep analysis report</title><style>' + STYLE + '</style></head><body><main>',
                 '<header><div class="eyebrow">PHYLOPREP · ANALYSIS REPORT</div><h1>From loci to analysis-ready matrices</h1>',
                 f'<p class="muted">Run: {esc(report["run_name"])} · Generated {esc(report["generated_at"])} · Report v{__version__}</p>',
                 '<nav><a href="#overview">Overview</a><a href="#qc">Quality control</a><a href="#matrices">Matrices</a><a href="#events">Sequence changes</a><a href="#coverage">Taxon coverage</a><a href="#methods">Methods &amp; provenance</a></nav>',
                 '<button id="download-summary">Download summary JSON</button><button onclick="window.print()">Print / save PDF</button></header>',
                 '<div class="cards">']
    for label, value in (("Input loci", overview["input_loci"]), ("Aligned loci", overview["aligned_loci"]),
                         ("Full concatenated matrices", len(full)), ("Evaluated matrix datasets", len(evaluated))):
        fragments.append(f'<div class="card"><div class="label">{label}</div><div class="metric">{esc(display(value))}</div></div>')
    fragments.append('</div><section id="overview"><div class="section-number">01 / PROCESSING</div><h2>Locus processing and retention</h2>')
    fragments.append('<p>Each QC branch and each matrix has its own locus set. Counts across raw/trimmed, codon/protein, and independent/paired outputs must not be added to estimate unique input loci.</p>')
    fragments.append(table(["Measure", "Value"], [("Input loci", overview["input_loci"]), ("Successfully aligned loci", overview["aligned_loci"]),
                                                  ("Recorded alignment failures", overview["failed_alignments"]), ("Reused aligned loci", overview["reused_loci"]),
                                                  ("Alignment count source", overview["alignment_count_source"])]))
    fragments.append('<p class="muted">Input minus aligned is not interpreted as failures: skipped, incomplete, unrecorded, and failed tasks cannot always be distinguished.</p>')
    fragments.append('<div class="chart">' + bar_chart("Actual concatenation inclusion by matrix", {m["group"]: m["loci"] for m in full if m.get("loci") is not None}) + '</div>')
    fragments.append('</section><section id="qc"><div class="section-number">02 / QUALITY CONTROL</div><h2>QC outcomes and overlapping exclusion reasons</h2>')
    fragments.append('<p class="note">A locus can violate several QC conditions. Unique failed loci are counted once per branch, while reason hits overlap. Their sum may exceed the failed-locus total, and reason percentages may sum to more than 100%. No arbitrary primary reason is assigned.</p><p class="muted">Paired reason labels summarize codon/protein QC outcomes and paired PIS decisions. Individual QC conditions are shown in the corresponding representation branches.</p>')
    fragments.append(table(["Branch", "Selection", "Passed", "Unique failed", "Not evaluated"],
                           [(q["branch"], q["selection"], q["passed"], q["failed"], q["not_evaluated"]) for q in report["qc"]]))
    if not report["qc"]:
        fragments.append('<p class="empty">QC records are unavailable.</p>')
    for q in report["qc"]:
        fragments.append(f'<details><summary>{esc(q["branch"])} — {q["failed"]:,} unique failed loci</summary>')
        fragments.append('<h3>Overlapping reason hits</h3><p class="muted">Denominator: unique failed loci in this branch.</p><div class="chart">' + bar_chart("Overlapping QC reason hits", q["reason_hits"], q["failed"]) + '</div>')
        fragments.append(table(["Reason", "Loci with reason", "% of unique failed loci"], [(k, v, pct(v / q["failed"]) if q["failed"] else "Not applicable") for k, v in sorted(q["reason_hits"].items())]))
        combinations = Counter(q["failure_combinations"])
        top = dict(combinations.most_common(10))
        rest = sum(combinations.values()) - sum(top.values())
        if rest:
            top["Other complete reason combinations"] = rest
        fragments.append('<h3>Mutually exclusive failure combinations</h3><p class="muted">Each failed locus belongs to exactly one complete combination. Counts sum to the unique failed-locus total.</p><div class="chart">' + bar_chart("Complete QC failure combinations", top, q["failed"]) + '</div>')
        fragments.append(table(["Complete combination", "Unique loci"], sorted(combinations.items())))
        fragments.append('</details>')
    fragments.append('</section><section id="matrices"><div class="section-number">03 / FINAL DATASETS</div><h2>Statistics of the actual final matrices</h2>')
    fragments.append('<p>Statistics are recomputed from each final FASTA matrix. QC pass counts do not imply inclusion: concatenation can subsequently skip loci with missing taxa. Codon-position lengths are nucleotide columns; full codon matrices also report complete codons.</p>')
    fragments.append(table(["Matrix", "Taxa", "Source loci", "Length", "Unit", "Complete codons", "Variable", "PIS", "Constant", "Unresolved-only", "Gap −", "Missing ?", "Ambiguous", "GC"],
                           [(m["matrix"], len(m["taxa"]), m.get("loci"), m["length"], m["unit"], m["complete_codons"], m["variable"], m["pis"], m["constant"], m["unresolved_only"], pct(m["gap_ratio"]), pct(m["missing_ratio"]), pct(m["ambiguous_ratio"]), "Not applicable" if m["seqtype"] == "prot" else pct(m["gc_ratio"])) for m in evaluated]))
    for m in full:
        if m.get("status") == "evaluated":
            fragments.append(f'<details><summary>{esc(m["group"])} — inclusion and missing blocks</summary>' + table(["Measure", "Value"], [
                ("Actually included loci", m.get("loci")), ("Skipped after QC", m.get("skipped_after_qc")),
                ("Artificial missing blocks", m["padding"]["blocks"] if m.get("padding") else None),
                ("Artificially filled characters", m["padding"]["characters"] if m.get("padding") else None)]) + '</details>')
    fragments.append('<p class="muted">Fourfold-degenerate source-locus counts are unavailable unless a contribution map is recorded; the parent matrix locus count is not substituted. Literal ? counts include any existing missing markers; only .missing.tsv identifies artificial padding.</p><h3>Raw versus trimmed: matched locus sets</h3>')
    fragments.append(table(["Representation / selection", "Common included loci", "Raw-only loci", "Trimmed-only loci", "Raw common-set columns", "Trimmed common-set columns", "Unit"],
                           [(c["raw"].rsplit("/", 1)[0], c["common_loci"], c.get("raw_only"), c.get("trimmed_only"), c["raw_common_columns"], c["trimmed_common_columns"], c.get("unit")) for c in report["trimming_comparisons"]]))
    fragments.append('<p class="muted">Trimming comparisons use the intersection of actually included gene partitions within each representation and selection. This controls locus composition for length comparisons; it does not enforce identical taxon sets. Paired refers to codon/protein pairing within a variant, not raw/trimmed pairing.</p></section>')
    fragments.append('<section id="events"><div class="section-number">04 / SEQUENCE PROVENANCE</div><h2>Sequence-change events by taxon</h2>')
    audit = report.get("audit")
    if audit is None:
        fragments.append('<p class="empty">Unavailable: the canonical reports/sequence_changes.tsv was not found. Per-locus copies are not combined because they could double-count events.</p>')
    else:
        fragments.append(f'<p>{audit["events"]:,} distinct recorded events. Exact duplicate audit rows ignored: {audit["duplicates_ignored"]:,}. Taxa use mapped identifiers when recorded, otherwise original sequence identifiers.</p>')
        fragments.append('<div class="chart">' + bar_chart("Sequence-change events by taxon", {t["taxon"]: t["events"] for t in audit["taxa"]}) + '</div>')
        fragments.append(table(["Taxon / identifier", "Events", "Affected sequences", "NT character-event spans", "AA character-event spans", "Events without coordinates"],
                               [(t["taxon"], t["events"], t["affected_sequences"] if t["taxon"] != "Unassigned" else None, t["nt_characters"], t["aa_characters"], t["events_without_coordinates"]) for t in audit["taxa"]]))
        fragments.append('<h3>Event types and actions</h3>' + table(["Event", "Action", "Count"], [(t["event"], t["action"], t["count"]) for t in audit["types"]]))
        fragments.append('<p class="muted">Source audit rows are counted once, not once per raw/trimmed or codon/protein output. Terminal cleanup and internal-stop events retain separate event types. Character-event spans sum coordinate lengths and can overlap; they are not unique edited positions. NT and AA spans are never added. Identifier changes remain separate event types. Original-sequence denominators and retained/excluded event attribution are not inferred from aligned sequences.</p>')
    fragments.append('</section><section id="coverage"><div class="section-number">05 / TAXON COVERAGE</div><h2>Observed loci versus artificial missing blocks</h2>')
    fragments.append('<p>Occupancy is the fraction of included loci for which a taxon was not artificially padded. An observed sequence containing gaps or ambiguous characters remains an observed locus. Gray cells indicate taxa absent from that matrix.</p><div class="chart">' + coverage_chart(matrices) + '</div></section>')
    fragments.append('<section id="methods"><div class="section-number">06 / INTERPRETATION</div><h2>Definitions, availability and provenance</h2>')
    fragments.append('<p>Resolved states are A/C/G/T for nucleotide columns and the 20 standard amino acids for protein columns. A variable column has at least two resolved states; a parsimony-informative site (PIS) has at least two states each observed at least twice. Constant columns have exactly one resolved state. Columns with no resolved states are reported separately. Singleton columns have exactly two resolved states with at least one occurring once. Gaps, ? and unresolved symbols do not contribute states.</p>')
    fragments.append('<p>Gap, missing and ambiguity percentages use all taxon-by-column cells as denominator. Ambiguity includes N for nucleotides, X for proteins, and any other unresolved symbols. GC uses resolved A/C/G/T only. Counts are descriptive and do not establish orthology or biological alignment accuracy.</p>')
    fragments.append('<p>Unavailable means the needed records were not present or could not be validated; it never means zero. A header-only missing-block table establishes zero artificial blocks. QC tables marked Recorded QC can precede independent PIS selection in older outputs. The report reads published product manifests when present and does not follow original absolute input paths. Provenance checksums are recorded values, not newly verified by this report; per-file input/output hashes remain in the run manifest.</p>')
    if report["provenance_source"]:
        fragments.append(f'<details><summary>Recorded run parameters and environment ({esc(report["provenance_source"])})</summary><pre>{esc(json.dumps(report["provenance"], indent=2, ensure_ascii=False))}</pre></details>')
    else:
        fragments.append('<p class="empty">Run parameters, versions, checksums and execution timings: not recorded in a portable run manifest. Report generation time is not analysis runtime.</p>')
    if report["warnings"]:
        fragments.append('<h3>Data availability warnings</h3><ul class="note warning">' + ''.join(f'<li>{esc(w)}</li>' for w in report["warnings"]) + '</ul>')
    payload = json.dumps(report, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    fragments.extend(['</section><footer>Generated locally with Python standard-library tools. This self-contained report makes no network requests. Embedded SVG charts can be saved individually.</footer>',
                      f'<script type="application/json" id="report-data">{payload}</script>',
                      """<script>
function saveBlob(blob,name){const a=document.createElement('a');const u=URL.createObjectURL(blob);a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),60000);}
document.getElementById('download-summary').addEventListener('click',()=>saveBlob(new Blob([document.getElementById('report-data').textContent],{type:'application/json'}),'phyloprep-report-summary.json'));
document.querySelectorAll('.chart svg').forEach((svg,i)=>{const b=document.createElement('button');b.textContent='Save chart as SVG';b.addEventListener('click',()=>saveBlob(new Blob([new XMLSerializer().serializeToString(svg)],{type:'image/svg+xml'}),'phyloprep-chart-'+(i+1)+'.svg'));svg.parentElement.appendChild(b);});
</script></main></body></html>"""])
    return ''.join(fragments)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path, help="Final PhyloPrep output directory")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    try:
        report = build_report(args.input)
        output = args.input.resolve() / "phyloprep-report.html"
        content = render_html(report)
        fd, temporary = tempfile.mkstemp(prefix=".phyloprep-report-", suffix=".html", dir=output.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(content)
            os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
        print(f"[PhyloPrep report] Saved {output}", flush=True)
        return 0
    except (ValueError, OSError, KeyError, csv.Error, sqlite3.Error) as error:
        print(f"PhyloPrep report: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
