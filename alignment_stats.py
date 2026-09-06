#!/usr/bin/env python3
"""Summarize FASTA alignment dimensions, missing data, and variable sites."""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
from Bio import SeqIO
from msap_io import atomic_output

def stats(path, seqtype="nucl"):
    records=list(SeqIO.parse(path,"fasta"));
    if not records: raise ValueError(f"{path}: no sequences found")
    seqs=[str(r.seq).upper() for r in records]; lengths=[len(s) for s in seqs]; aligned=len(set(lengths))==1; nchar=lengths[0] if aligned else 0
    ambiguity_char = "X" if seqtype == "prot" else "N"
    ambiguity_ratio = sum(s.count(ambiguity_char) for s in seqs)/(sum(lengths) or 1)
    variable=parsimony=singleton=constant=distinct_patterns=0
    if aligned:
        patterns = []
        for i in range(nchar):
            raw_col = tuple(s[i] for s in seqs)
            patterns.append(raw_col)
            missing = "-?N" if seqtype in {"nucl", "codon"} else "-?"
            col=[base for base in raw_col if base not in missing]
            states=set(col)
            if len(states) <= 1:
                constant += 1
            else:
                variable+=1
                if len(col)>=2 and sum(col.count(x)>=2 for x in states)>=2: parsimony+=1
                elif len(col) and any(col.count(x) == 1 for x in states): singleton += 1
        distinct_patterns = len(set(patterns))
    result = {"file":str(path),"seqtype":seqtype,"sequence_count":len(seqs),
              "alignment_length":nchar,"gap_ratio":sum(s.count("-") for s in seqs)/(sum(lengths) or 1),
              "variable_sites":variable,"parsimony_informative_sites":parsimony,
              "distinct_patterns":distinct_patterns,"singleton_sites":singleton,
              "constant_sites":constant}
    result["N_ratio" if seqtype in {"nucl", "codon"} else "X_ratio"] = ambiguity_ratio
    return result, aligned

def main(argv=None):
    parser=argparse.ArgumentParser(description="Summarize FASTA alignments for phylogenetic QC.",add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter,epilog="""Examples:
  alignment_stats.py -i alignments/*.fasta -o qc/alignment_stats.tsv -st codon
  alignment_stats.py -i trimmed/*.fasta -o qc/stats.tsv -st nucl""")
    required=parser.add_argument_group("required arguments"); optional=parser.add_argument_group("statistics options")
    required.add_argument("-i","--input",nargs="+",required=True,metavar="FASTA",help="Input FASTA alignment files."); required.add_argument("-st","--seqtype",required=True,choices=["nucl","prot","codon"],help="Sequence type: nucl, prot, or codon.")
    optional.add_argument("-o","--output",default="-",metavar="TSV",help="Output TSV report (default: stdout).")
    optional.add_argument("-h","--help",action="help",help="Show this help message and exit."); optional.add_argument("-v","--version",action="version",version="alignment_stats v1.00",help="Show program's version number and exit.")
    args=parser.parse_args(argv)
    try:
        checked=[stats(path,args.seqtype) for path in args.input]
        rows=[row for row, _aligned in checked]
        for path, (_row, aligned) in zip(args.input, checked):
            if not aligned:
                print(f"Warning: {path} is not an alignment matrix; sequence lengths differ and data are not aligned.", file=sys.stderr)
        keys=["File", "Seqtype", "Sequence_count", "Alignment_length",
              "N_ratio" if args.seqtype in {"nucl", "codon"} else "X_ratio",
              "Gap_ratio", "Variable_sites", "Parsimony_informative_sites",
              "Distinct_patterns", "Singleton_sites", "Constant_sites"]
        aliases = {"File": "file", "Seqtype": "seqtype", "Sequence_count": "sequence_count",
                   "Alignment_length": "alignment_length", "N_ratio": "N_ratio", "X_ratio": "X_ratio",
                   "Gap_ratio": "gap_ratio", "Variable_sites": "variable_sites",
                   "Parsimony_informative_sites": "parsimony_informative_sites",
                   "Distinct_patterns": "distinct_patterns", "Singleton_sites": "singleton_sites",
                   "Constant_sites": "constant_sites"}
        rows = [{key: row[aliases[key]] for key in keys} for row in rows]
        with atomic_output(None if args.output == "-" else args.output) as handle:
            writer=csv.DictWriter(handle,fieldnames=keys,delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    except (OSError,ValueError) as error: print(f"alignment_stats: {error}",file=sys.stderr); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
