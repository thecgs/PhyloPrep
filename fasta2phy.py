#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from msap_io import atomic_output, read_fasta_alignment, check_output_path
from sequence_audit import normalize_export_records, write_report

def fasta2phy(infile, outfile, seqtype='DNA', changes_report=None):

    records = read_fasta_alignment(infile)
    changes, protein = normalize_export_records(records, infile, seqtype)
    if changes:
        if changes_report is None:
            if outfile is None:
                raise ValueError('Use --changes-report when normalizing input while writing PHYLIP to stdout')
            changes_report = str(outfile) + '.sequence_changes.tsv'
        check_output_path(changes_report, inputs=[infile, outfile])
    with atomic_output(outfile, inputs=[infile]) as out:
        out.write(f'{len(records)} {len(records[0].seq)}\n')
        for record in records:
            out.write(f'{record.id}  {str(record.seq).upper()}\n')
    if changes:
        write_report(changes_report, changes, protein=protein)
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
    optional.add_argument('-st', '--seqtype', default='DNA', choices=['DNA', 'protein'], help='DNA masks ambiguity as N; protein masks ambiguity as X.')
    optional.add_argument('--changes-report', metavar='TSV', help='Audit TSV; defaults to OUTPUT.sequence_changes.tsv when edits occur.')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='fasta2phy v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    fasta2phy(args.input, args.output, args.seqtype, args.changes_report)
