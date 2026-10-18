#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path
from msap_io import atomic_output, check_output_path
from fasta2nex import prepare_records

def write_phylip(matrix, out):
    """Write an already prepared alignment without parsing or normalizing it again."""
    out.write(f'{len(matrix.records)} {matrix.length}\n')
    for index, (identifier, _, _) in enumerate(matrix.records):
        out.write(f'{identifier}  ')
        for chunk in matrix.chunks(index):
            out.write(chunk)
        out.write('\n')

def fasta2phy(infile, outfile, seqtype='DNA', changes_report=None):

    check_output_path(outfile, inputs=[infile])
    matrix, changes = prepare_records(infile, seqtype, Path(outfile).absolute().parent if outfile else None)
    try:
        if changes:
            if changes_report is None:
                if outfile is None:
                    raise ValueError('Use --changes-report when normalizing input while writing PHYLIP to stdout')
                changes_report = str(outfile) + '.sequence_changes.tsv'
            check_output_path(changes_report, inputs=[infile, outfile])
        with atomic_output(outfile, inputs=[infile]) as out:
            write_phylip(matrix, out)
        if changes:
            changes.publish(changes_report)
    finally:
        matrix.close()
        changes.close()
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
    optional.add_argument('-st', '--seqtype', default='DNA', choices=['DNA', 'protein'], help='Preserves missing-data ? and gaps; DNA masks ambiguity as N, protein as X.')
    optional.add_argument('--changes-report', metavar='TSV', help='Audit TSV; defaults to OUTPUT.sequence_changes.tsv when edits occur.')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0', help="Show program's version number and exit.")
    args = parser.parse_args()
    fasta2phy(args.input, args.output, args.seqtype, args.changes_report)
