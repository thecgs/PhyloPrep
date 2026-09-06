#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from msap_io import atomic_output, check_output_path, read_fasta_alignment

def get_prefix(infile):
    return os.path.splitext(os.path.basename(infile))[0]

def split(infile):

    records = read_fasta_alignment(infile, codon=True)
    outputs = [f"codon{position}.{get_prefix(infile)}.fasta" for position in ("1st", "2nd", "3rd")]
    for outfile in outputs:
        check_output_path(outfile, [infile])
    for offset, outfile in enumerate(outputs):
        with atomic_output(outfile, inputs=[infile]) as out:
            for record in records:
                print(f">{record.id}\n{record.seq[offset::3]}", file=out)
    
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Split a codon FASTA alignment into first, second, and third positions.',
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  split_codon_seqence_alignment.py -i codon.alignment.fasta')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='CODON_FASTA',
                          help='Input codon alignment in FASTA format.', required=True)
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='split_codon_seqence_alignment v1.00',
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    split(infile=args.input)
