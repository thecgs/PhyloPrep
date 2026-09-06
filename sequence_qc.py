#!/usr/bin/env python3
"""Check individual FASTA sequences and write a QC report."""

from __future__ import annotations

import argparse
import html
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from msap_io import atomic_output


DNA_CHARS = set("ACGTURYSWKMBDHVN-?")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*?-" )


def _records(path: str):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError("No sequences found in input FASTA.")
    return records


def sequence_qc(path: str, seqtype: str, genetic_code: int = 1,
                    min_length: int = 1) -> dict:
    """Return summary, sequence and column findings for one FASTA file."""
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"Input FASTA file does not exist: {path}")
    records = _records(path)
    actual_type = seqtype
    allowed = DNA_CHARS if actual_type in {"nucl", "codon"} else PROTEIN_CHARS
    ids = [record.id for record in records]
    duplicate_ids = sorted(identifier for identifier, count in Counter(ids).items() if count > 1)
    lengths = [len(record.seq) for record in records]
    aligned = len(set(lengths)) == 1
    alignment_length = lengths[0] if aligned else 0
    sequence_rows = []
    all_invalid = []
    for record in records:
        sequence = str(record.seq).upper()
        invalid = sorted(set(sequence) - allowed)
        all_invalid.extend(invalid)
        row = {
            "id": record.id,
            "length": len(sequence),
            "gap_count": sequence.count("-"),
            "gap_ratio": sequence.count("-") / len(sequence) if sequence else 0.0,
            "n_count": sequence.count("N"),
            "n_ratio": sequence.count("N") / len(sequence) if sequence else 0.0,
            "x_count": sequence.count("X"),
            "x_ratio": sequence.count("X") / len(sequence) if sequence else 0.0,
            "invalid": "".join(invalid),
            "short": len(sequence) < min_length,
        }
        if actual_type == "codon":
            row["frame_ok"] = len(sequence) % 3 == 0
            row["terminal_stop"] = False
            row["internal_stop_count"] = 0
            if row["frame_ok"] and sequence and "-" not in sequence:
                table = CodonTable.unambiguous_dna_by_id[genetic_code]
                codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
                row["terminal_stop"] = codons[-1] if codons[-1] in table.stop_codons else "NA"
                hits = [f"{i + 1}-{codon}" for i, codon in enumerate(codons[:-1]) if codon in table.stop_codons]
                row["internal_stop_count"] = ";".join(hits) if hits else "NA"
        sequence_rows.append(row)

    full_gap_columns = []
    column_gap_ratios = []
    if aligned and alignment_length:
        columns = [str(record.seq).upper() for record in records]
        for index in range(alignment_length):
            column = [sequence[index] for sequence in columns]
            gap_ratio = column.count("-") / len(column)
            column_gap_ratios.append(gap_ratio)
            if gap_ratio == 1.0:
                full_gap_columns.append(index + 1)

    valid_sites = alignment_length - len(full_gap_columns) if aligned else 0
    if actual_type == "codon" and aligned and alignment_length % 3 == 0:
        valid_sites = sum(
            1 for start in range(0, alignment_length - 2, 3)
            if "-" not in "".join(str(record.seq).upper()[start:start + 3] for record in records)
        )
    issues = []
    if duplicate_ids:
        issues.append("duplicate_ids")
    if not aligned:
        issues.append("length_mismatch")
    if all_invalid:
        issues.append("invalid_characters")
    if any(not row["length"] for row in sequence_rows):
        issues.append("empty_sequence")
    if any(row["short"] for row in sequence_rows):
        issues.append("short_sequence")
    if full_gap_columns:
        issues.append("all_gap_columns")
    if actual_type == "codon":
        if any(not row["frame_ok"] for row in sequence_rows):
            issues.append("frame_error")
        if any(row["terminal_stop"] != "NA" or row["internal_stop_count"] != "NA" for row in sequence_rows):
            issues.append("stop_codon")
    return {
        "file": str(source), "seqtype": actual_type, "sequence_count": len(records),
        "aligned": aligned, "alignment_length": alignment_length,
        "valid_sites": valid_sites, "full_gap_columns": full_gap_columns,
        "max_gap_ratio": max(column_gap_ratios, default=0.0),
        "issues": issues, "sequences": sequence_rows,
    }


SUMMARY_COLUMNS = ["file", "seqtype", "sequence_count", "aligned", "alignment_length",
                   "valid_sites", "full_gap_columns", "max_gap_ratio", "issues"]
SEQUENCE_COLUMNS = ["id", "length", "gap_count", "gap_ratio", "n_count", "n_ratio",
                    "x_count", "x_ratio", "invalid", "short", "frame_ok",
                    "terminal_stop", "internal_stop_count"]


def write_tsv(results, output):
    with atomic_output(output) as handle:
        protein = results[0]["seqtype"] == "prot"
        codon = results[0]["seqtype"] == "codon"
        columns = ["File", "Taxa_ID", "Sequence_length", "X_count" if protein else "N_count", "Gap_count"]
        if codon:
            columns += ["Terminal_stop", "Internal_stop_count"]
        handle.write("\t".join(columns) + "\n")
        for result in results:
            for row in result["sequences"]:
                terminal = "NA"
                internal = "NA"
                if result["seqtype"] == "codon":
                    terminal = row["terminal_stop"]
                    internal = row["internal_stop_count"]
                    if len(str(next(r.seq for r in _records(result["file"]) if r.id == row["id"]))) % 3:
                        terminal = "not_multiple_of_3"
                values = [result["file"], row["id"], row["length"],
                          row["x_count"] if protein else row["n_count"], row["gap_count"]]
                if codon:
                    values += [terminal, internal]
                handle.write("\t".join(map(str, values)) + "\n")


def write_html(results, output):
    protein = results[0]["seqtype"] == "prot"
    codon = results[0]["seqtype"] == "codon"
    columns = ["File", "Taxa_ID", "Sequence_length", "X_count" if protein else "N_count", "Gap_count"]
    if codon:
        columns += ["Terminal_stop", "Internal_stop_count"]
    rows = []
    for result in results:
        for row in result["sequences"]:
            values = [result["file"], row["id"], row["length"], row["x_count"] if protein else row["n_count"], row["gap_count"]]
            if codon:
                values += [row.get("terminal_stop", "NA"), row.get("internal_stop_count", "NA")]
            if codon and len(str(next(r.seq for r in _records(result["file"]) if r.id == row["id"]))) % 3:
                values[-2] = "not_multiple_of_3"
                values[-1] = "NA"
            rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in values) + "</tr>")
    body = "<table><thead><tr>" + "".join(f"<th>{html.escape(key)}</th>" for key in columns) + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    with atomic_output(output) as handle:
        handle.write("<!doctype html><meta charset='utf-8'><title>Alignment QC</title>" + body)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="""Check individual FASTA sequences and write a QC report.

The sequence type must be specified explicitly. Codon mode additionally checks
reading frame and terminal/internal stop codons.""",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  sequence_qc.py -i alignment.fasta -st codon -o qc/alignment.tsv
  sequence_qc.py -i alignments/*.fasta -st nucl -o qc/alignment.html --format html
  sequence_qc.py -i proteins.fasta -st prot -o qc/proteins.tsv --min-length 50

Notes:
  * -st/--seqtype is required and must be codon, prot, or nucl.
  * TSV is the default report format; .html and .htm outputs select HTML.
  * -g/--genetic-code uses NCBI genetic-code table IDs (default: 1).""",
    )
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("quality-control options")
    required.add_argument("-i", "--input", nargs="+", required=True, metavar="FASTA",
                          help="Input FASTA files.")
    optional.add_argument("-o", "--output", default="-", metavar="REPORT",
                          help="Output TSV or HTML report (default: stdout).")
    required.add_argument("-st", "--seqtype", required=True,
                          choices=["nucl", "prot", "codon"],
                          help="Sequence type: nucl, prot, or codon.")
    optional.add_argument("-g", "--genetic-code", type=int, default=1, metavar="TABLE",
                          help="NCBI genetic-code table for codon mode (default: 1).")
    optional.add_argument("--min-length", type=int, default=1, metavar="N",
                          help="Flag sequences shorter than N (default: 1).")
    optional.add_argument("--format", choices=["tsv", "html"], default=None,
                          help="Report format; inferred from output extension when omitted.")
    optional.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    optional.add_argument("-v", "--version", action="version", version="sequence_qc v1.00",
                          help="Show program's version number and exit.")
    args = parser.parse_args(argv)
    if args.min_length < 0:
        parser.error("--min-length must be non-negative")
    try:
        results = [sequence_qc(path, args.seqtype, args.genetic_code, args.min_length) for path in args.input]
        for result in results:
            if result["issues"]:
                print(f"Warning: {result['file']}: {', '.join(result['issues'])}", file=sys.stderr)
        output = None if args.output == "-" else Path(args.output)
        report_format = args.format or ("html" if output and output.suffix.lower() in {".html", ".htm"} else "tsv")
        if report_format == "html":
            write_html(results, output)
        else:
            write_tsv(results, output)
    except (OSError, ValueError, KeyError) as error:
        print(f"sequence_qc: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
