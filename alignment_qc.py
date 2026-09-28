#!/usr/bin/env python3
"""Check individual FASTA sequences and write a QC report."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from msap_io import atomic_output, check_output_path, is_stop_codon
from msap_io import expand_input_paths


def length_units(length, seqtype):
    """Return the length in the units used by the --min-length threshold."""
    return length // 3 if seqtype == "codon" else length

def alignment_stats(path, seqtype="nucl"):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"{path}: no sequences found")
    seqs = [str(r.seq).upper() if seqtype == 'prot' else str(r.seq).upper().replace('U', 'T') for r in records]
    lengths = [len(s) for s in seqs]
    aligned = len(set(lengths)) == 1
    nchar = lengths[0] if aligned else 0
    ambiguity = "X" if seqtype == "prot" else "N"
    if seqtype == "codon":
        # A codon containing N translates to one ambiguous protein residue.
        # Count it once, rather than counting its individual N bases, so the
        # codon --max-n-ratio threshold corresponds to the protein X ratio.
        ambiguity_count = sum(
            "N" in sequence[start:start + 3]
            for sequence in seqs
            for start in range(0, len(sequence) - 2, 3)
        )
        ambiguity_total = sum(len(sequence) // 3 for sequence in seqs)
    else:
        ambiguity_count = sum(s.count(ambiguity) for s in seqs)
        ambiguity_total = sum(lengths)
    result = {"file": str(path), "seqtype": seqtype, "sequence_count": len(seqs),
              "alignment_length": nchar, "gap_ratio": sum(s.count("-") for s in seqs) / (sum(lengths) or 1),
              "variable_sites": 0, "parsimony_informative_sites": 0,
              "distinct_patterns": 0, "singleton_sites": 0, "constant_sites": 0,
              ("X_ratio" if seqtype == "prot" else "N_ratio"): ambiguity_count / (ambiguity_total or 1)}
    if aligned:
        patterns = [tuple(s[i] for s in seqs) for i in range(nchar)]
        # Count only resolved states; ambiguity is not an additional allele.
        alphabet = set("ACDEFGHIKLMNPQRSTVWYOU") if seqtype == "prot" else set("ACGT")
        for col in patterns:
            valid = [x if seqtype == "prot" else x.replace("U", "T") for x in col]
            valid = [x for x in valid if x in alphabet]
            states = set(valid)
            if len(states) <= 1: result["constant_sites"] += 1
            else:
                result["variable_sites"] += 1
                if sum(valid.count(x) >= 2 for x in states) >= 2: result["parsimony_informative_sites"] += 1
                elif any(valid.count(x) == 1 for x in states): result["singleton_sites"] += 1
        result["distinct_patterns"] = len(set(patterns))
    return result, aligned


DNA_CHARS = set("ACGTURYSWKMBDHVN-?")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*?-" )


def _records(path: str):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError("No sequences found in input FASTA.")
    return records


def alignment_qc(path: str, seqtype: str, genetic_code: int = 1,
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
            # In codon mode a threshold of N means N complete codons, so it
            # has the same biological unit as the translated protein matrix.
            "short": length_units(len(sequence), actual_type) < min_length,
        }
        if actual_type == "codon":
            codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
            row["frame_ok"] = len(sequence) % 3 == 0 and all(
                "-" not in codon or codon == "---" for codon in codons)
            row["terminal_stop"] = "NA"
            row["internal_stop_count"] = "NA"
            if len(sequence) % 3 == 0 and sequence:
                occupied = [i for i, codon in enumerate(codons) if codon != '---']
                last = occupied[-1] if occupied else -1
                row["terminal_stop"] = codons[last] if last >= 0 and is_stop_codon(codons[last], genetic_code) else "NA"
                hits = [f"{i + 1}-{codon}" for i, codon in enumerate(codons) if i != last and is_stop_codon(codon, genetic_code)]
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

def validate_report_paths(prefix, inputs):
    outputs = [str(prefix) + suffix for suffix in (".pass.tsv", ".fail.tsv", ".details.tsv")]
    for i, output in enumerate(outputs):
        check_output_path(output, list(inputs) + outputs[:i])


def write_qc_reports(results, passed, failed, prefix, seqtype, allow_terminal, allow_internal):
    validate_report_paths(prefix, [result["file"] for result in results])
    summary = ["File", "Seqtype", "Sequence_count", "Alignment_length",
               "N_ratio" if seqtype != "prot" else "X_ratio", "Gap_ratio",
               "Variable_sites", "Parsimony_informative_sites", "Distinct_patterns",
               "Singleton_sites", "Constant_sites"]
    def rows(items):
        for r in items:
            s = r["alignment_stats"]
            yield [str(Path(r["file"]).resolve()), s["seqtype"], s["sequence_count"],
                   s["alignment_length"], s["N_ratio" if seqtype != "prot" else "X_ratio"],
                   s["gap_ratio"], s["variable_sites"], s["parsimony_informative_sites"],
                   s["distinct_patterns"], s["singleton_sites"], s["constant_sites"]]
    for suffix, items in ((".pass.tsv", passed), (".fail.tsv", failed)):
        with atomic_output(str(prefix) + suffix) as handle:
            w = csv.writer(handle, delimiter="\t"); w.writerows(([row[0]] for row in rows(items)))
    detail = summary + (["Terminal_stop", "Internal_stop_count"] if seqtype == "codon" else []) + ["Issues"]
    with atomic_output(str(prefix) + ".details.tsv") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(detail)
        for r in results:
            for row in rows([r]):
                extra = []
                if seqtype == "codon":
                    stops = [(x["id"], x["terminal_stop"]) for x in r["sequences"] if x["terminal_stop"] != "NA"]
                    internal = [(x["id"], x["internal_stop_count"]) for x in r["sequences"] if x["internal_stop_count"] != "NA"]
                    extra = [";".join(f"{a}-{b}" for a,b in stops) or "NA", ";".join(f"{a}-inter-{b}" for a,b in internal) or "NA"]
                w.writerow(row + extra + [";".join(r["reasons"])])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="""Check FASTA alignment matrices and write pass/fail QC reports.

The sequence type must be specified explicitly. Codon mode additionally checks
reading frame and terminal/internal stop codons.""",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  alignment_qc.py -i alignment.fasta -st codon -prefix qc/alignment
  alignment_qc.py -i alignments/*.fasta -st nucl -prefix qc/alignment
  alignment_qc.py -i proteins.fasta -st prot -prefix qc/proteins --min-length 50

Notes:
  * -st/--seqtype is required and must be codon, prot, or nucl.
  * --min-length is applied to every sequence and to the aligned matrix length.
    Its unit is nucleotides for nucl, amino acids for prot, and complete codons
    for codon. Thus --min-length 100 applies the same 100-residue threshold to
    paired codon and protein matrices in phyloprep.
  * --max-n-ratio applies to nucl/codon matrices; --max-x-ratio applies to prot.
    In nucl mode it is the fraction of N bases. In codon mode it is the
    fraction of complete codons containing at least one N, matching the
    fraction of X residues in the translated protein matrix.
    Ratios are calculated over the complete matrix: characters for nucl/prot,
    and complete codons for codon ambiguity; --max-gap-ratio is always the
    corresponding character-level gap ratio.
  * --min-parsimony-informative-sites requires the whole matrix to contain at
    least N parsimony-informative sites; use 0 to disable this requirement.
    Site statistics count resolved states only (DNA A/C/G/T, U treated as T;
    protein residues excluding B/J/X/Z, gaps, unknowns and stop symbols).
  * Codon gaps must occupy complete triplets (---); partial triplet gaps fail
    with frame_error. This check cannot be disabled by stop-codon options.
  * Stop checks normalize U to T and recognize IUPAC codons only when every
    possible expansion is a stop in the selected genetic code (e.g. TAR in 1).
    Dual-coding codons also encoding amino acids are not definite stops in an
    alignment, even at its end; complete-CDS preprocessing handles terminal stops.
  * In codon mode, terminal and internal stop codons fail QC by default.
    --allow-terminal-stop and --allow-internal-stop independently permit them.
  * Output files are PREFIX.pass.tsv, PREFIX.fail.tsv and PREFIX.details.tsv.
    The first two contain absolute input paths without headers; details.tsv
    contains the complete statistics table.
  * -g/--genetic-code uses NCBI genetic-code table IDs (default: 1).""",
    )
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("quality-control options")
    required.add_argument("-i", "--input", nargs="+", required=True, metavar="FASTA|LIST",
                          help="FASTA files and/or path-list files (one FASTA path per line).")
    required.add_argument("-p", "-prefix", "--prefix", required=True, metavar="PREFIX",
                          help="Output prefix; writes .pass.tsv, .fail.tsv and .details.tsv.")
    required.add_argument("-st", "--seqtype", required=True,
                          choices=["nucl", "prot", "codon"],
                          help="Sequence type: nucl, prot, or codon.")
    optional.add_argument("-g", "--genetic-code", type=int, default=1, metavar="TABLE",
                          help="NCBI genetic-code table for codon mode (default: 1).")
    optional.add_argument("--min-length", type=int, default=1, metavar="N",
                          help="Minimum nt (nucl), aa (prot), or complete codons (codon) for every sequence and matrix (default: 1).")
    optional.add_argument("--max-n-ratio", type=float, default=1.0, metavar="RATIO",
                          help="Maximum N-base ratio (nucl) or N-containing codon ratio (codon) (default: 1).")
    optional.add_argument("--max-x-ratio", type=float, default=1.0, metavar="RATIO",
                          help="Maximum X ratio for protein matrices (default: 1).")
    optional.add_argument("--max-gap-ratio", type=float, default=1.0, metavar="RATIO",
                          help="Maximum gap ratio for the complete matrix (default: 1).")
    optional.add_argument("--min-parsimony-informative-sites", type=int, default=0, metavar="N",
                          help="Minimum parsimony-informative sites in the matrix (default: 0).")
    optional.add_argument("--allow-terminal-stop", action="store_true",
                          help="Allow terminal stop codons in codon mode.")
    optional.add_argument("--allow-internal-stop", action="store_true",
                          help="Allow internal stop codons in codon mode.")
    optional.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    optional.add_argument("-v", "--version", action="version", version="v1.0.0",
                          help="Show program's version number and exit.")
    args = parser.parse_args(argv)
    if args.min_length < 0 or args.min_parsimony_informative_sites < 0 or any(
            value < 0 or value > 1 for value in (args.max_n_ratio, args.max_x_ratio, args.max_gap_ratio)):
        parser.error("length thresholds must be non-negative and ratios must be between 0 and 1")
    try:
        paths = expand_input_paths(args.input)
        validate_report_paths(args.prefix, [*args.input, *paths])
        results, passed, failed = [], [], []
        for path in paths:
            result = alignment_qc(path, args.seqtype, args.genetic_code, args.min_length)
            stat, aligned = alignment_stats(path, args.seqtype)
            result["alignment_stats"] = stat
            reasons = [x for x in result["issues"] if x != "stop_codon"]
            if args.seqtype == "codon":
                has_terminal = any(x["terminal_stop"] != "NA" for x in result["sequences"])
                has_internal = any(x["internal_stop_count"] != "NA" for x in result["sequences"])
                if (has_terminal and not args.allow_terminal_stop) or (has_internal and not args.allow_internal_stop):
                    reasons.append("stop_codon")
            if not aligned: reasons.append("length_mismatch")
            if length_units(stat["alignment_length"], args.seqtype) < args.min_length:
                reasons.append("short_alignment")
            if stat["gap_ratio"] > args.max_gap_ratio: reasons.append("gap_ratio")
            amb_key = "X_ratio" if args.seqtype == "prot" else "N_ratio"
            if stat[amb_key] > (args.max_x_ratio if args.seqtype == "prot" else args.max_n_ratio): reasons.append(amb_key)
            if stat["parsimony_informative_sites"] < args.min_parsimony_informative_sites: reasons.append("parsimony_informative_sites")
            result["pass"] = not reasons; result["reasons"] = sorted(set(reasons))
            (passed if result["pass"] else failed).append(result); results.append(result)
        write_qc_reports(results, passed, failed, args.prefix, args.seqtype, args.allow_terminal_stop, args.allow_internal_stop)
    except (OSError, ValueError, KeyError) as error:
        print(f"alignment_qc: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
