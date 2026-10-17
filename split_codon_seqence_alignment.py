#!/usr/bin/env python
# coding: utf-8

import os
import argparse
from contextlib import ExitStack
from msap_matrix import AlignmentSpool
from msap_io import atomic_output, check_output_path

def get_prefix(infile):
    return os.path.splitext(os.path.basename(infile))[0]

def split(infile):

    outputs = [f"codon{position}.{get_prefix(infile)}.fasta" for position in ("1st", "2nd", "3rd")]
    for outfile in outputs:
        check_output_path(outfile, [infile])
    with AlignmentSpool(infile, codon=True, directory=os.getcwd()) as matrix, ExitStack() as stack:
        handles = [stack.enter_context(atomic_output(path, inputs=[infile])) for path in outputs]
        for index, (identifier, _, _) in enumerate(matrix.records):
            for handle in handles:
                handle.write(f'>{identifier}\n')
            for chunk in matrix.chunks(index):
                for offset, handle in enumerate(handles):
                    handle.write(chunk[offset::3])
            for handle in handles:
                handle.write('\n')

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
    optional.add_argument('-v', '--version', action='version', version='v1.0.0',
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    split(infile=args.input)
