#!/usr/bin/env python
# coding: utf-8

import argparse, sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from msap_io import atomic_output, read_fasta_alignment


def fasta2nex(infile, outfile=None, seqtype="DNA"):
    records = read_fasta_alignment(infile)
    for record in records:
        record.annotations["molecule_type"] = seqtype
    with atomic_output(outfile, inputs=[infile]) as handle:
        AlignIO.write(MultipleSeqAlignment(records), handle, "nexus")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a FASTA alignment to NEXUS format.", add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  fasta2nex.py -i alignment.fasta -o alignment.nex -st DNA")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA', required=True, help='Input FASTA alignment.')
    optional.add_argument('-o', '--output', metavar='NEXUS', help='Output NEXUS file (default: stdout).', default=None)
    optional.add_argument('-st', '--seqtype', '--seq-type', default='DNA', choices=["DNA", "RNA", "protein"], help='Alignment sequence type (default: DNA).')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='fasta2nex v2.00', help="Show program's version number and exit.")
    args = parser.parse_args()

    fasta2nex(args.input, args.output, args.seqtype)
