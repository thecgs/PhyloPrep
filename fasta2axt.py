#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from msap_io import atomic_output, read_fasta_alignment

def fasta2axt(inputfile, output=None):
    records = read_fasta_alignment(inputfile, count=2)
    with atomic_output(output, inputs=[inputfile]) as axt_out:
        print(">" + "|".join(record.id for record in records), file=axt_out)
        for record in records:
            print(record.seq, file=axt_out)
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert exactly two aligned FASTA sequences to simplified Ka/Ks AXT format.', add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  fasta2axt.py -i alignment.fasta -o alignment.axt')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA', required=True,
                          help='Input FASTA alignment.')
    optional.add_argument('-o', '--output', metavar='AXT',
                          help='Output AXT file (default: stdout).', default=None)
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='fasta2axt v1.00',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    fasta2axt(inputfile=args.input, output=args.output)
