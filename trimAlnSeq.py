#!/usr/bin/env python
# coding: utf-8

import sys
import gzip
import argparse
from Bio import SeqIO
from msap_io import atomic_output, check_output_path


def open_input(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def open_output(filename):
    if filename == "-":
        return sys.stdout
    elif filename.endswith(".gz"):
        return gzip.open(filename, "wt")
    else:
        return open(filename, "w")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trim fasta format alignment sequence. For codon mode, filtering is performed per codon triplet.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Date:2026/06/21 Author:Guisen Chen Email:thecgs001@foxmail.com"
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-i", "--input",
        metavar="str",
        required=True,
        help="Input FASTA alignment, optionally gzip-compressed; use - for stdin."
    )

    optional.add_argument(
        "-o", "--output",
        metavar="str",
        default="-",
        help="Output FASTA file, optionally .gz; - or omitted means stdout. Zero retained sites is an error."
    )

    optional.add_argument(
        "-st", "--seqtype",
        metavar="str",
        default="codon",
        choices=["codon", "prot", "nucl"],
        help="Sequence type: codon, prot or nucl. default=codon"
    )

    optional.add_argument(
        "-G", "--G",
        default=0.2,
        type=float,
        metavar="float",
        help="Maximum gap ratio per site. range 0-1. default=0.2"
    )

    optional.add_argument(
        "-N", "--N",
        default=0.2,
        type=float,
        metavar="float",
        help="Maximum N ratio per site for nucleotide/codon. default=0.2"
    )

    optional.add_argument(
        "-X", "--X",
        default=0.2,
        type=float,
        metavar="float",
        help="Maximum X ratio per site for protein. default=0.2"
    )

    optional.add_argument(
        "-ra", "--removeall",
        action="store_true",
        help="Remove invariant columns."
    )

    optional.add_argument(
        "-h", "--help",
        action="help",
        help="Show help message and exit."
    )

    optional.add_argument(
        "-v", "--version",
        action="version",
        version="v2.02"
    )

    args = parser.parse_args()

    infile = args.input
    outfile = args.output
    st = args.seqtype
    G = args.G
    N = args.N
    X = args.X
    ra = args.removeall

    if not (0 <= G <= 1):
        raise ValueError("-G must be between 0 and 1")

    if not (0 <= N <= 1):
        raise ValueError("-N must be between 0 and 1")

    if not (0 <= X <= 1):
        raise ValueError("-X must be between 0 and 1")

    # --------------------------------------------------
    # Read alignment
    # --------------------------------------------------

    handle = open_input(infile)

    ids = []
    seqs = []

    expected_len = None

    for record in SeqIO.parse(handle, "fasta"):
        seq = str(record.seq).upper()

        if expected_len is None:
            expected_len = len(seq)
        elif len(seq) != expected_len:
            raise ValueError(
                f"Alignment length mismatch: {record.id} has length "
                f"{len(seq)}, expected {expected_len}"
            )

        ids.append(record.id)
        seqs.append(seq)

    if handle is not sys.stdin:
        handle.close()

    seqnum = len(seqs)

    if seqnum == 0:
        sys.exit("No sequences found.")

    if st == "codon" and expected_len % 3 != 0:
        raise ValueError(
            f"Codon alignment length ({expected_len}) "
            "is not divisible by 3."
        )

    # --------------------------------------------------
    # Determine columns to retain
    # --------------------------------------------------

    keep = []

    if st == "codon":

        nsites = expected_len // 3

        for site in range(nsites):
            start = site * 3
            end = start + 3

            gap_count = 0
            n_count = 0

            first = None
            invariant = True

            for seq in seqs:
                codon = seq[start:end]

                if first is None:
                    first = codon
                elif codon != first:
                    invariant = False

                # Gap and N thresholds are independent, including partial gaps.
                if "-" in codon:
                    gap_count += 1
                if "N" in codon:
                    n_count += 1

            if gap_count / seqnum > G:
                continue

            if n_count / seqnum > N:
                continue

            if ra and invariant:
                continue

            keep.append((start, end))

    else:

        nsites = expected_len

        for site in range(nsites):

            gap_count = 0
            bad_count = 0

            first = seqs[0][site]
            invariant = True

            for seq in seqs:
                base = seq[site]

                if base != first:
                    invariant = False

                if base == "-":
                    gap_count += 1

                if st == "prot":
                    if base == "X":
                        bad_count += 1
                else:
                    if base == "N":
                        bad_count += 1

            if gap_count / seqnum > G:
                continue

            if st == "prot":
                if bad_count / seqnum > X:
                    continue
            else:
                if bad_count / seqnum > N:
                    continue

            if ra and invariant:
                continue

            keep.append(site)

    # --------------------------------------------------
    # Output
    # --------------------------------------------------

    if not keep:
        sys.exit("No alignment sites remain after trimming; no output was written.")

    destination = None if outfile == "-" else outfile
    sources = [] if infile == "-" else [infile]
    with atomic_output(destination, inputs=sources) as out:
        if st == "codon":
            for name, seq in zip(ids, seqs):
                trimmed = "".join(seq[start:end] for start, end in keep)
                out.write(f">{name}\n{trimmed}\n")
        else:
            for name, seq in zip(ids, seqs):
                trimmed = "".join(seq[i] for i in keep)
                out.write(f">{name}\n{trimmed}\n")
