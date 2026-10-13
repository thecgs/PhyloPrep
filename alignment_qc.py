#!/usr/bin/env python3
"""Check individual FASTA sequences and write a QC report."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from msap_io import atomic_output, check_output_path, is_stop_codon
from msap_io import expand_input_paths


def length_units(length, seqtype):
    """Return the length in the units used by the --min-length threshold."""
    return length // 3 if seqtype == "codon" else length

def alignment_stats(path, seqtype="nucl", records=None):
    records = list(records) if records is not None else list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"{path}: no sequences found")
    seqs = [str(r.seq).upper() if seqtype == 'prot' else str(r.seq).upper().replace('U', 'T') for r in records]
    lengths = [len(s) for s in seqs]
    aligned = len(set(lengths)) == 1
    nchar = lengths[0] if aligned else 0
    ambiguity = "X" if seqtype == "prot" else "N"
    if seqtype == "codon":
        # Count each N-containing codon once rather than counting individual N
        # bases. This keeps the threshold in codon units. It is intentionally
        # independent of protein X: ambiguous codons such as GCN can translate
        # deterministically and therefore need not produce X.
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


def _audit_locus_key(path: str, source_input=False) -> str:
    """Map source/audit/alignment filenames from one locus to a stable key."""
    name = Path(path).name
    for suffix in (".macse_original.NT.fasta", ".macse_original.AA.fasta",
                   ".macse_original.NT.fa", ".macse_original.AA.fa"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    if source_input:
        return Path(name).stem
    match = re.fullmatch(
        r"(.+)\.(?:mafft|muscle|prank|clustalw2|macse)\.(?:nucl|codon|prot)(?:\.trimmed)?\.aln",
        name,
    )
    return match.group(1) if match else Path(name).stem


def load_sequence_changes(paths):
    """Index sequence-change provenance by locus and final taxon ID."""
    indexed = defaultdict(lambda: defaultdict(lambda: {"terminal": [], "internal": []}))
    sources_by_locus = defaultdict(set)
    for path in paths or []:
        source = Path(path)
        if not source.is_file():
            raise ValueError(f"Sequence-change report does not exist: {path}")
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"input_file", "taxa_id", "event", "source_nt_positions", "original"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f"{path}: not a PhyloPrep sequence_changes.tsv report")
            for row in reader:
                event_name = row.get("event", "")
                if event_name not in {"terminal_stop", "internal_stop",
                                     "protein_terminal_stop", "protein_internal_stop"}:
                    continue
                input_file = row.get("input_file", "")
                locus = _audit_locus_key(input_file, source_input=True)
                if input_file:
                    sources_by_locus[locus].add(str(Path(input_file).resolve()))
                taxon = row.get("mapped_taxa_id")
                if not taxon or taxon == "NA":
                    taxon = row.get("taxa_id", "NA")
                category = "internal" if "internal_stop" in event_name else "terminal"
                position = row.get("source_nt_positions", "NA")
                if position in {"", "NA"}:
                    position = row.get("source_aa_positions", "NA")
                original = row.get("original", "NA")
                indexed[locus][taxon][category].append(f"{position}-{original}")
    collisions = {key: values for key, values in sources_by_locus.items() if len(values) > 1}
    if collisions:
        details = "; ".join(
            f"{key}: {', '.join(sorted(values))}" for key, values in sorted(collisions.items())
        )
        raise ValueError(
            "Cannot match --sequence-changes unambiguously because multiple source loci "
            f"share the same basename-derived locus key ({details}). Use unique locus basenames."
        )
    return indexed


def _audit_for_alignment(path, audit_index):
    return audit_index.get(_audit_locus_key(path), {}) if audit_index else {}


def alignment_qc(path: str, seqtype: str, genetic_code: int = 1,
                    audit_index=None) -> dict:
    """Return summary, sequence and column findings for one FASTA file."""
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"Input FASTA file does not exist: {path}")
    records = _records(path)
    actual_type = seqtype
    stat, aligned = alignment_stats(path, actual_type, records=records)
    audit_rows = _audit_for_alignment(path, audit_index) if actual_type in {"codon", "prot"} else {}
    allowed = DNA_CHARS if actual_type in {"nucl", "codon"} else PROTEIN_CHARS
    ids = [record.id for record in records]
    duplicate_ids = sorted(identifier for identifier, count in Counter(ids).items() if count > 1)
    lengths = [len(record.seq) for record in records]
    alignment_length = stat["alignment_length"]
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
            "audited_terminal_stop": ";".join(audit_rows.get(record.id, {}).get("terminal", [])) or "NA",
            "audited_internal_stop": ";".join(audit_rows.get(record.id, {}).get("internal", [])) or "NA",
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
        elif actual_type == "prot":
            occupied = [i for i, residue in enumerate(sequence) if residue != "-"]
            last = occupied[-1] if occupied else -1
            row["terminal_stop"] = "*" if last >= 0 and sequence[last] == "*" else "NA"
            hits = [f"{i + 1}-*" for i, residue in enumerate(sequence)
                    if i != last and residue == "*"]
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
            if all("-" not in sequence[start:start + 3] for sequence in columns)
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
    if full_gap_columns:
        issues.append("all_gap_columns")
    if actual_type == "codon":
        if any(not row["frame_ok"] for row in sequence_rows):
            issues.append("frame_error")
    if actual_type in {"codon", "prot"} and any(
            row.get("terminal_stop", "NA") != "NA" or
            row.get("internal_stop_count", "NA") != "NA"
            for row in sequence_rows):
        issues.append("stop_codon")
    return {
        "file": str(source), "seqtype": actual_type, "sequence_count": len(records),
        "aligned": aligned, "alignment_length": alignment_length,
        "valid_sites": valid_sites, "full_gap_columns": full_gap_columns,
        "max_gap_ratio": max(column_gap_ratios, default=0.0),
        "issues": issues, "sequences": sequence_rows,
        "alignment_stats": stat,
        "audited_terminal_stop": any(row["audited_terminal_stop"] != "NA" for row in sequence_rows),
        "audited_internal_stop": any(row["audited_internal_stop"] != "NA" for row in sequence_rows),
    }


def validate_report_paths(prefix, inputs):
    outputs = [str(prefix) + suffix for suffix in (".pass.tsv", ".fail.tsv", ".details.tsv")]
    for i, output in enumerate(outputs):
        check_output_path(output, list(inputs) + outputs[:i])


def write_qc_reports(results, passed, failed, prefix, seqtype):
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
    detail = summary + (["Terminal_stop", "Internal_stop_count"] if seqtype in {"codon", "prot"} else []) + ["Audited_terminal_stop", "Audited_internal_stop", "Issues"]
    with atomic_output(str(prefix) + ".details.tsv") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(detail)
        for r in results:
            for row in rows([r]):
                extra = []
                if seqtype in {"codon", "prot"}:
                    stops = [(x["id"], x.get("terminal_stop", "NA")) for x in r["sequences"] if x.get("terminal_stop", "NA") != "NA"]
                    internal = [(x["id"], x.get("internal_stop_count", "NA")) for x in r["sequences"] if x.get("internal_stop_count", "NA") != "NA"]
                    extra = [";".join(f"{a}-{b}" for a,b in stops) or "NA", ";".join(f"{a}-inter-{b}" for a,b in internal) or "NA"]
                audited_terminal = [(x["id"], x["audited_terminal_stop"]) for x in r["sequences"] if x["audited_terminal_stop"] != "NA"]
                audited_internal = [(x["id"], x["audited_internal_stop"]) for x in r["sequences"] if x["audited_internal_stop"] != "NA"]
                extra += [
                    ";".join(f"{a}-{b}" for a, b in audited_terminal) or "NA",
                    ";".join(f"{a}-{b}" for a, b in audited_internal) or "NA",
                ]
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
  * --min-length is applied to the alignment length. Its unit is nucleotides
    for nucl, amino acids for prot, and complete codons for codon. Unequal
    sequence lengths already fail as length_mismatch, so no second per-sequence
    length gate is applied.
  * --max-n-ratio applies to nucl/codon matrices; --max-x-ratio applies to prot.
    In nucl mode it is the fraction of N bases. In codon mode it is the
    fraction of complete codons containing at least one N. This need not equal
    the protein X ratio because some ambiguous codons translate deterministically.
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
    Dual-coding codons also encoding amino acids are not definite stops in a
    final alignment, even when they occupy its last non-gap position. Source-CDS
    preprocessing uses terminal context and records removed terminal stops.
  * In codon and protein modes, stop symbols observed directly in the supplied
    alignment fail QC by default. --allow-terminal-stop permits a terminal stop;
    --allow-internal-stop permits internal stops. Protein direct-stop detection
    uses '*'; custom stop symbols should be normalized before QC.
  * --sequence-changes accepts one or more PhyloPrep sequence_changes.tsv reports.
    Masked internal stops recorded by preprocessing are restored as QC provenance
    for both codon and paired protein alignments. They fail QC by default and are
    permitted only by --allow-internal-stop. Terminal stops that MSAP removed from
    complete CDS inputs are reported as audited provenance but do not fail QC.
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
    optional.add_argument("--sequence-changes", nargs="+", metavar="TSV",
                          help="Optional PhyloPrep sequence_changes.tsv report(s). Masked internal stops are restored as QC provenance for codon/protein alignments.")
    optional.add_argument("--allow-terminal-stop", action="store_true",
                          help="Allow terminal stops observed directly in codon or protein alignments. MSAP-removed terminal stops are provenance only.")
    optional.add_argument("--allow-internal-stop", action="store_true",
                          help="Allow internal stops observed directly in codon/protein alignments or recorded as masked in --sequence-changes.")
    optional.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    optional.add_argument("-v", "--version", action="version", version="v1.0.0",
                          help="Show program's version number and exit.")
    args = parser.parse_args(argv)
    if args.min_length < 0 or args.min_parsimony_informative_sites < 0 or any(
            value < 0 or value > 1 for value in (args.max_n_ratio, args.max_x_ratio, args.max_gap_ratio)):
        parser.error("length thresholds must be non-negative and ratios must be between 0 and 1")
    try:
        paths = expand_input_paths(args.input)
        audit_index = load_sequence_changes(args.sequence_changes)
        validate_report_paths(args.prefix, [*args.input, *paths, *(args.sequence_changes or [])])
        results, passed, failed = [], [], []
        for path in paths:
            result = alignment_qc(path, args.seqtype, args.genetic_code, audit_index)
            stat = result["alignment_stats"]
            reasons = [x for x in result["issues"] if x != "stop_codon"]
            observed_terminal = args.seqtype in {"codon", "prot"} and any(
                x.get("terminal_stop", "NA") != "NA" for x in result["sequences"])
            observed_internal = args.seqtype in {"codon", "prot"} and any(
                x.get("internal_stop_count", "NA") != "NA" for x in result["sequences"])
            audited_internal = args.seqtype in {"codon", "prot"} and any(
                x["audited_internal_stop"] != "NA" for x in result["sequences"])
            if observed_terminal and not args.allow_terminal_stop:
                reasons.append("terminal_stop")
            if (observed_internal or audited_internal) and not args.allow_internal_stop:
                reasons.append("internal_stop")
            if length_units(stat["alignment_length"], args.seqtype) < args.min_length:
                reasons.append("short_alignment")
            if stat["gap_ratio"] > args.max_gap_ratio: reasons.append("gap_ratio")
            amb_key = "X_ratio" if args.seqtype == "prot" else "N_ratio"
            if stat[amb_key] > (args.max_x_ratio if args.seqtype == "prot" else args.max_n_ratio): reasons.append(amb_key)
            if stat["parsimony_informative_sites"] < args.min_parsimony_informative_sites: reasons.append("parsimony_informative_sites")
            result["pass"] = not reasons; result["reasons"] = sorted(set(reasons))
            (passed if result["pass"] else failed).append(result); results.append(result)
        write_qc_reports(results, passed, failed, args.prefix, args.seqtype)
    except (OSError, ValueError, KeyError) as error:
        print(f"alignment_qc: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
