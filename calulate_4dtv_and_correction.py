#!/usr/bin/env python
# coding: utf-8

import os
import sys
import math
import argparse
import tempfile
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict

def get_four_fold_codons(genetic_code=1):
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]
    print('The condon table used:', codon_table, file=sys.stderr)
    four_fold_codons = {}
    acids2codons = defaultdict(list)
    for codon in codon_table.forward_table:
        acids2codons[codon_table.forward_table[codon]].append(codon)
        
    for acid in acids2codons:
        tmp = defaultdict(list)
        for codon in acids2codons[acid]:
            tmp[codon[0:2]].append(codon)
        for i in tmp:
            if len(tmp[i]) == 4:
                for c in tmp[i]:
                    four_fold_codons.setdefault(c, acid)
    
    return four_fold_codons

def get_len(seq):
    if len(seq)%3 == 0:
        l = len(seq)
    elif len(seq)%3 == 1:
        l = len(seq) - 1
    else:
        l = len(seq) - 2
    return l

def is_transversion(base1, base2):
    transversion = {'A': 'TC', 'G': 'TC', 'C': 'AG', 'T':'AG'}
    if base2 in transversion[base1]:
        res = True
    else:
        res = False
    return res

def calculate_4DTV_correction(infile, four_fold_codons, error=sys.stderr):
    if not os.path.isfile(infile):
        raise ValueError(f"Input file does not exist: {infile}")
    file_prefix = os.path.splitext(os.path.basename(infile))[0]
    records = list(SeqIO.parse(infile, 'fasta'))
    if len(records) != 2:
        raise ValueError(
            f"{infile}: expected exactly two sequences, found {len(records)}"
        )
    if records[0].id == records[1].id:
        raise ValueError(f"{infile}: sequence IDs must be different")

    seq1_name = records[0].id
    seq2_name = records[1].id
    seq1 = str(records[0].seq).upper()
    seq2 = str(records[1].seq).upper()
    if not seq1 or not seq2:
        raise ValueError(f"{infile}: sequences must not be empty")
    if len(seq1) != len(seq2):
        raise ValueError(
            f"{infile}: sequence lengths differ ({len(seq1)} and {len(seq2)})"
        )
    if len(seq1) % 3 != 0:
        raise ValueError(
            f"{infile}: alignment length {len(seq1)} is not divisible by 3"
        )
    
    seq = ''
    fourfold_sites_total_number = 0
    fourfold_sites_transversion_number = 0
    for i in range(0, get_len(seq1), 3):
        codon1 = seq1[i: i+3]
        codon2 = seq2[i: i+3]
        if (codon1 in four_fold_codons) and (codon2 in four_fold_codons) and (codon1[0:2] == codon2[0:2]):
            fourfold_sites_total_number += 1
            seq += codon1[2]
            seq += codon2[2]
            if is_transversion(codon1[2], codon2[2]):
                fourfold_sites_transversion_number += 1

    if fourfold_sites_total_number == 0:
        return (
            file_prefix, seq1_name, seq2_name, "NA", "NA",
            fourfold_sites_total_number, fourfold_sites_transversion_number,
        )

    raw_4dtv = fourfold_sites_transversion_number / fourfold_sites_total_number

    A = 0.5*seq.count('A')/fourfold_sites_total_number
    C = 0.5*seq.count('C')/fourfold_sites_total_number
    G = 0.5*seq.count('G')/fourfold_sites_total_number
    T = 0.5*seq.count('T')/fourfold_sites_total_number
    Y = 0.5*(seq.count('T') + seq.count('C'))/fourfold_sites_total_number
    R = 0.5*(seq.count('A') + seq.count('G'))/fourfold_sites_total_number

    if (A != 0) and (C != 0) and (G != 0) and (T != 0) and (Y != 0) and (R != 0):
        if (1-raw_4dtv*(T*C*R/Y+A*G*Y/R)/(2*(T*C*R+A*G*Y))) == 0:
            corrected_4dtv = "NA"
        else:
            a= -1*math.log(1-raw_4dtv*(T*C*R/Y+A*G*Y/R)/(2*(T*C*R+A*G*Y)))
            if (1-raw_4dtv/(2*Y*R) > 0):
                b=-1*math.log(1-raw_4dtv/(2*Y*R))
                corrected_4dtv=2*a*(T*C/Y+A*G/R)-2*b*(T*C*R/Y+A*G*Y/R-Y*R)
            else:
                corrected_4dtv = "NA"
    else:
        corrected_4dtv = 'NA'
    return file_prefix, seq1_name, seq2_name, corrected_4dtv, raw_4dtv, fourfold_sites_total_number, fourfold_sites_transversion_number

def main(infiles, outfile, genetic_code=1):
    four_fold_codons = get_four_fold_codons(genetic_code=genetic_code)
    # Calculate every input before touching the destination file. This keeps
    # an existing result intact when one alignment is invalid.
    results = []
    for infile in infiles:
        results.append(calculate_4DTV_correction(infile, four_fold_codons))

    output_dir = os.path.dirname(os.path.abspath(outfile)) or "."
    os.makedirs(output_dir, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=output_dir,
            prefix=".4dtv-", suffix=".tmp", delete=False
        ) as out:
            temporary_name = out.name
            print('Input file prefix\tSequence1 name\tSequence2 name\tcorrected_4dtv\traw_4dtv\tfourfold_sites_total_number\tfourfold_sites_transversion_number', file=out)
            for result in results:
                print(*result, sep='\t', file=out)
        os.replace(temporary_name, outfile)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)
    return None
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate raw and HKY-corrected 4DTV from codon FASTA alignments.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  calulate_4dtv_and_correction.py -i gene1.fasta gene2.fasta -o 4dtv.tsv -g 1

The genetic-code table uses NCBI IDs (default: 1).""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='CODON_FASTA', help='Input codon alignments in FASTA format.', required=True, nargs='+')
    required.add_argument('-o', '--output', metavar='TSV', help='Output TSV file.', required=True)
    optional.add_argument('-g', '--genetic_code', '--genetic-code', metavar='TABLE', default=1, type=int,
                          help='NCBI genetic-code table (default: 1).')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='calulate_4dtv_and_correction v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    main(infiles=args.input, outfile=args.output, genetic_code=args.genetic_code)
