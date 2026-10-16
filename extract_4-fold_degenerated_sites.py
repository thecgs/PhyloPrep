#!/usr/bin/env python
# coding: utf-8

import argparse
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict
from msap_io import atomic_output

def get_four_fold_codons(genetic_code=1):
    """
    Return a set containing all four-fold degenerate codons.
    """
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]

    # amino acid -> {first_two_bases: count}
    aa_prefix_count = defaultdict(lambda: defaultdict(int))

    for codon, aa in codon_table.forward_table.items():
        aa_prefix_count[aa][codon[:2]] += 1

    four_fold_codons = set()

    for codon, aa in codon_table.forward_table.items():
        if aa_prefix_count[aa][codon[:2]] == 4:
            four_fold_codons.add(codon)

    return four_fold_codons
    
def extract_four_fold_site(infile, outfile, genetic_code=1):
    """
    Extract four-fold degenerate codons from a codon-aligned FASTA.

    This version:
      - avoids dict/Counter/list creation inside the codon loop
      - uses list.append() instead of repeated string concatenation
      - stops checking a site immediately when one sequence fails
      - converts Bio.Seq objects to str only once
    """

    ids = []
    seqs = []

    for record in SeqIO.parse(infile, "fasta"):
        ids.append(record.id)
        seqs.append(str(record.seq).upper().replace('U', 'T'))

    if not seqs:
        raise ValueError("No sequences found in input FASTA.")

    nseq = len(seqs)
    matrix_length = len(seqs[0])

    # Check alignment lengths once
    for seq in seqs[1:]:
        if len(seq) != matrix_length:
            raise ValueError("Sequences have different lengths.")

    if matrix_length % 3 != 0:
        raise ValueError(
            f"Alignment length ({matrix_length}) is not divisible by 3."
        )

    four_fold_codons = get_four_fold_codons(genetic_code)

    # One list for each sequence
    output = [[] for _ in range(nseq)]

    for pos in range(0, matrix_length, 3):
        #print(str(round((pos+3)/matrix_length*100, 2)) + " %", end='\r')
        
        # Use the first sequence as reference
        codon0 = seqs[0][pos:pos + 3]

        # If reference itself is not four-fold, reject immediately
        if codon0 not in four_fold_codons:
            continue

        prefix = codon0[:2]

        valid = True

        # Start from sequence 2
        for j in range(1, nseq):
            codon = seqs[j][pos:pos + 3]

            # This checks both:
            # 1. same first two bases
            # 2. valid four-fold codon
            if codon[:2] != prefix or codon not in four_fold_codons:
                valid = False
                break

        if not valid:
            continue

        # Site passed all checks.
        # Append complete codon.
        for j in range(nseq):
            #output[j].append(seqs[j][pos:pos + 3])
            output[j].append(seqs[j][pos + 2])

    with atomic_output(outfile, inputs=[infile]) as out:
        for seq_id, fragments in zip(ids, output):
            out.write(f">{seq_id}\n")
            out.write("".join(fragments))
            out.write("\n")
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract third positions from four-fold degenerate codon sites in a codon FASTA alignment.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  extract_4-fold_degenerated_sites.py -i codon.aln.fasta -o fourfold.fasta -g 1

The genetic-code table uses NCBI IDs (default: 1).""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='CODON_FASTA', help='Input codon alignment in FASTA format.', required=True)
    required.add_argument('-o', '--output', metavar='FASTA', help='Output four-fold-site alignment in FASTA format.', required=True)
    optional.add_argument('-g', '--genetic_code', '--genetic-code', metavar='TABLE', default=1, type=int,
                          help='NCBI genetic-code table (default: 1).')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0', help="Show program's version number and exit.")
    args = parser.parse_args()
    extract_four_fold_site(infile=args.input, outfile=args.output, genetic_code=args.genetic_code)
