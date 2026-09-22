#!/usr/bin/env python
# coding: utf-8

import argparse
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import Counter, defaultdict
from msap_io import atomic_output

#def get_four_fold_codons(genetic_code=1):
#    """
#    Return a set containing all four-fold degenerate codons.
#    """
#    
#    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]
#    
#    #print('The condon table used:', codon_table, file=sys.stderr)
#    
#    
#    
#    four_fold_codons = {}
#    acids2codons = defaultdict(list)
#    
#    for codon in codon_table.forward_table:
#        acids2codons[codon_table.forward_table[codon]].append(codon)
#        
#    for acid in acids2codons:
#        tmp = defaultdict(list)
#        for codon in acids2codons[acid]:
#            tmp[codon[0:2]].append(codon)
#        for i in tmp:
#            if len(tmp[i]) == 4:
#                for c in tmp[i]:
#                    four_fold_codons.setdefault(c, acid)
#    
#    return four_fold_codons
 
#def extract_four_fold_site(infile, outfile, genetic_code):
#    """
#    Extract four-fold degenerate codons from a codon-aligned FASTA.
#
#    This version:
#      - avoids dict/Counter/list creation inside the codon loop
#      - uses list.append() instead of repeated string concatenation
#      - stops checking a site immediately when one sequence fails
#      - converts Bio.Seq objects to str only once
#    """
#    
#    alignment_seqences = {}
#    for record in SeqIO.parse(infile, 'fasta'):
#        alignment_seqences.setdefault(record.id, record.seq.upper())
#    
#    four_fold_codons = get_four_fold_codons(genetic_code)
#    alignment_seqences_four_fold_site = defaultdict(list)
#    matrix_length = len(record.seq)
#    for i in range(0, matrix_length, 3):
#        print(str(round((i+3)/matrix_length*100, 2)) + " %", end='\r')
#        codons = {}
#        for geneID in alignment_seqences:
#            codons[geneID] = alignment_seqences[geneID][i: i+3]
#
#        if len(Counter([c[0:2] for c in codons.values()])) == 1:
#            codons_bool = [True if c in four_fold_codons else False for c in codons.values()]
#            if all(codons_bool):
#                for geneID in alignment_seqences:
#                    alignment_seqences_four_fold_site[geneID].append(codons[geneID][2])
#                    
#    out = open(outfile, 'w')
#    for k in alignment_seqences_four_fold_site:
#        seq = ''.join(alignment_seqences_four_fold_site[k])
#        print(f'>{k}\n{seq}', file=out)
#    out.close()
#    return None 

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
    optional.add_argument('-v', '--version', action='version', version='extract_4-fold_degenerated_sites v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    extract_four_fold_site(infile=args.input, outfile=args.output, genetic_code=args.genetic_code)
