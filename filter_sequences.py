#!/usr/bin/env python3
"""Filter FASTA sequences by length, ambiguity, and taxon lists."""
from __future__ import annotations
import argparse, sys
from Bio import SeqIO
from msap_io import atomic_output

def _read_ids(path):
    with open(path, encoding="utf-8") as handle:
        return {line.strip().split()[0] for line in handle if line.strip() and not line.startswith("#")}

def _remove_empty_sites(records, seqtype):
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        return 0
    length = lengths.pop()
    if seqtype == "codon" and length % 3:
        raise ValueError("codon alignment length is not divisible by 3")
    if seqtype == "codon":
        kept = []
        for start in range(0, length, 3):
            codons = [str(r.seq)[start:start + 3].upper() for r in records]
            if not all(c in {"NNN", "---"} for c in codons):
                kept.extend(range(start, start + 3))
    else:
        missing = {"N", "-"} if seqtype == "nucl" else {"X", "-"}
        kept = [i for i in range(length) if not all(str(r.seq)[i].upper() in missing for r in records)]
    removed = length - len(kept)
    if removed:
        for r in records:
            r.seq = r.seq.__class__("".join(str(r.seq)[i] for i in kept))
    return removed

def filter_fasta(infile, outfile, seqtype, min_length=0, max_gap_ratio=1.0,
                 max_n_ratio=1.0, max_x_ratio=1.0, keep_ids=None,
                 remove_ids=None):
    records=list(SeqIO.parse(infile,"fasta"));
    if not records: raise ValueError("No sequences found in input FASTA.")
    if keep_ids is not None and not keep_ids:
        raise ValueError("--keep-taxa list is empty; refusing to retain all taxa")
    remove_ids=remove_ids or set(); selected=[]
    for record in records:
        sequence=str(record.seq).upper(); length=len(sequence)
        if (length < min_length or (sequence.count("-")/length if length else 1)>max_gap_ratio or
                (seqtype in {"nucl", "codon"} and (sequence.count("N")/length if length else 1)>max_n_ratio) or
                (seqtype == "prot" and (sequence.count("X")/length if length else 1)>max_x_ratio)): continue
        if keep_ids is not None and record.id not in keep_ids: continue
        if record.id in remove_ids: continue
        selected.append(record)
    if not selected: raise ValueError("No sequences remain after filtering.")
    if (keep_ids or remove_ids) and len(selected) < 2:
        raise ValueError("taxon filtering retained fewer than two sequences; refusing to write a matrix")
    removed_sites = _remove_empty_sites(selected, seqtype)
    if any(not record.seq for record in selected):
        raise ValueError("No alignment sites remain after filtering; no output was written")
    with atomic_output(outfile,inputs=[infile]) as handle: SeqIO.write(selected,handle,"fasta")
    return {"input":len(records),"output":len(selected), "removed":len(records)-len(selected), "sites_removed":removed_sites}

def main(argv=None):
    parser=argparse.ArgumentParser(description="Filter complete FASTA sequences (taxa) by quality and taxon lists.",add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter,epilog="""Examples:
  filter_sequences.py -i alignment.fasta -o filtered.fasta -st nucl --min-length 100
  filter_sequences.py -i alignment.fasta -o filtered.fasta -st nucl --max-gap-ratio 0.5 --max-n-ratio 0.1
  filter_sequences.py -i proteins.fasta -o filtered.fasta -st prot --max-x-ratio 0.2
  filter_sequences.py -i alignment.fasta -o filtered.fasta -st codon --keep-taxa taxa.txt
  filter_sequences.py -i alignment.fasta -o filtered.fasta -st nucl --remove-taxa taxa.txt

Important:
  This script removes whole sequences/taxa. N, gap, and X ratios are calculated
  within each complete sequence. After taxa filtering, columns containing only
  missing states are removed: N/- for nucl, X/- for prot, and NNN/--- codons
  for codon. Zero retained sites is an error and leaves existing output intact.
  Other site-level trimming remains the job of trimAlnSeq.py.""")
    required=parser.add_argument_group("required arguments"); optional=parser.add_argument_group("filtering options")
    required.add_argument("-i","--input",required=True,metavar="FASTA",help="Input FASTA file."); required.add_argument("-o","--output",required=True,metavar="FASTA",help="Filtered FASTA file."); required.add_argument("-st","--seqtype","--seq-type",required=True,choices=["nucl","prot","codon"],help="Sequence type: nucl, prot, or codon.")
    optional.add_argument("--min-length",type=int,default=0,metavar="N",help="Remove sequences shorter than N."); optional.add_argument("--max-gap-ratio",type=float,default=1.0,metavar="RATIO",help="Maximum per-sequence gap ratio (default: 1)."); optional.add_argument("--max-n-ratio",type=float,default=1.0,metavar="RATIO",help="Maximum per-sequence N ratio (default: 1)."); optional.add_argument("--max-x-ratio",type=float,default=1.0,metavar="RATIO",help="Maximum per-sequence X ratio in protein mode (default: 1).")
    optional.add_argument("--keep-taxa", metavar="TXT", help="Taxon ID list file; retain only listed IDs.")
    optional.add_argument("--remove-taxa", metavar="TXT", help="Taxon ID list file; remove listed IDs.")
    optional.add_argument("-h","--help",action="help",help="Show this help message and exit."); optional.add_argument("-v","--version",action="version",version="v1.0.0",help="Show program's version number and exit.")
    args=parser.parse_args(argv)
    if args.min_length<0 or not 0<=args.max_gap_ratio<=1 or not 0<=args.max_n_ratio<=1 or not 0<=args.max_x_ratio<=1: parser.error("length must be non-negative and ratios must be between 0 and 1")
    try:
        keep_ids = _read_ids(args.keep_taxa) if args.keep_taxa else None
        remove_ids = _read_ids(args.remove_taxa) if args.remove_taxa else None
        input_ids = {record.id for record in SeqIO.parse(args.input, "fasta")}
        for label, ids in (("--keep-taxa", keep_ids), ("--remove-taxa", remove_ids)):
            if ids:
                absent = sorted(ids - input_ids)
                if absent:
                    print(f"filter_sequences: warning: {label} contains IDs absent from input: "
                          + ", ".join(absent), file=sys.stderr)
        if keep_ids is not None and remove_ids is not None:
            overlap = sorted(keep_ids & remove_ids)
            if overlap:
                raise ValueError("--keep-taxa and --remove-taxa overlap: " + ", ".join(overlap))
        filter_fasta(args.input, args.output, args.seqtype, args.min_length,
                     args.max_gap_ratio, args.max_n_ratio, args.max_x_ratio,
                     keep_ids, remove_ids)
    except (OSError,ValueError) as error: print(f"filter_sequences: {error}",file=sys.stderr); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
