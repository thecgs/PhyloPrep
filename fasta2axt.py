#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from msap_io import atomic_output, read_fasta_alignment, check_output_path
from sequence_audit import normalize_export_records, write_report

def fasta2axt(inputfile, output=None, seqtype='DNA', changes_report=None):
    records = read_fasta_alignment(inputfile, count=2)
    changes, protein = normalize_export_records(records, inputfile, seqtype)
    if changes:
        if changes_report is None:
            if output is None:
                raise ValueError('Use --changes-report when normalizing input while writing AXT to stdout')
            changes_report = str(output) + '.sequence_changes.tsv'
        check_output_path(changes_report, inputs=[inputfile, output])
    with atomic_output(output, inputs=[inputfile]) as axt_out:
        print(">" + "|".join(record.id for record in records), file=axt_out)
        for record in records:
            print(record.seq, file=axt_out)
    if changes:
        write_report(changes_report, changes, protein=protein)
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
    optional.add_argument('-st', '--seqtype', default='DNA', choices=['DNA', 'protein'],
                          help='DNA masks ambiguity as N; protein masks ambiguity as X.')
    optional.add_argument('--changes-report', metavar='TSV',
                          help='Audit TSV; defaults to OUTPUT.sequence_changes.tsv when edits occur.')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    fasta2axt(args.input, args.output, args.seqtype, args.changes_report)
