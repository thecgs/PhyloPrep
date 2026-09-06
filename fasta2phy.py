#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from msap_io import atomic_output, read_fasta_alignment

def fasta2phy(infile, outfile):

    records = read_fasta_alignment(infile)
    alignment = MultipleSeqAlignment([record.upper() for record in records])
    with atomic_output(outfile, inputs=[infile]) as out:
        AlignIO.write(alignment, out, "phylip-relaxed")
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a FASTA alignment to relaxed PHYLIP format.', add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  fasta2phy.py -i alignment.fasta -o alignment.phy')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA', help='Input FASTA alignment.', required=True)
    optional.add_argument('-o', '--output', metavar='PHY', help='Output PHYLIP file (default: stdout).', default=None)
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='fasta2phy v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    fasta2phy(infile=args.input, outfile=args.output)
