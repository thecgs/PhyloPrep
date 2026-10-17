#!/usr/bin/env python
# coding: utf-8

"""Convert aligned FASTA to interleaved NEXUS without using Biopython.

The Biopython NEXUS reader builds interleaved alignments by repeatedly
concatenating short blocks. That makes its own round-trip check impractical
for multi-million-column supermatrices. This module writes a conventional
interleaved matrix from a temporary disk spool in bounded pages instead.
"""

import argparse
import csv
import tempfile
from pathlib import Path
from msap_matrix import AlignmentSpool

from msap_io import atomic_output, check_output_path
from sequence_audit import FIELDS, event, normalize_nucleotides


_NEXUS_PUNCTUATION = set("()[]{}\\,;:=*'\"`+-<>")
_NEXUS_WHITESPACE = set(" \t\n")
_LINE_WIDTH = 70


def nexus_name(identifier):
    """Quote an identifier when required by the NEXUS token grammar."""
    escaped = identifier.replace("'", "''")
    if set(escaped).intersection(_NEXUS_PUNCTUATION | _NEXUS_WHITESPACE):
        return f"'{escaped}'"
    return escaped


def normalize_sequence(source, seqtype, infile, identifier, changes, offset=0):
    if seqtype == 'DNA':
        sequence = normalize_nucleotides(source)
        for index, (original, replacement) in enumerate(zip(source, sequence), offset + 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier,
                                     'uracil_to_thymine' if original.upper() == 'U' else 'nucleotide_ambiguity',
                                     'NA', source_nt_positions=f'{index}-{index}',
                                     original=original, replacement=replacement, action='replace'))
    elif seqtype == 'RNA':
        sequence = normalize_nucleotides(source).replace('T', 'U')
        for index, (original, replacement) in enumerate(zip(source, sequence), offset + 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier,
                                     'thymine_to_uracil' if original.upper() == 'T' else 'nucleotide_ambiguity',
                                     'NA', source_nt_positions=f'{index}-{index}',
                                     original=original, replacement=replacement, action='replace'))
    else:
        sequence = ''.join(char.upper() if char.upper() in 'ACDEFGHIKLMNPQRSTVWYX-' else 'X'
                           for char in source)
        for index, (original, replacement) in enumerate(zip(source, sequence), offset + 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier, 'protein_unsupported_residue', 'NA',
                                     source_aa_positions=f'{index}-{index}', original=original,
                                     replacement='X', action='replace'))
    return sequence


class ChangeLog:
    """Spool export edits too; ambiguity-heavy matrices can have many edits."""
    def __init__(self, directory, protein):
        self.file = tempfile.TemporaryFile(mode='w+t', encoding='utf-8', newline='',
                                          dir=directory)
        self.count = 0
        fields = FIELDS + (['source_aa_positions'] if protein else [])
        self.writer = csv.DictWriter(self.file, fieldnames=fields, delimiter='\t', restval='NA')
        self.writer.writeheader()

    def append(self, row):
        self.writer.writerow(row)
        self.count += 1

    def __bool__(self):
        return bool(self.count)

    def publish(self, path):
        self.file.seek(0)
        with atomic_output(path) as handle:
            while True:
                chunk = self.file.read(65536)
                if not chunk:
                    break
                handle.write(chunk)

    def close(self):
        self.file.close()


def prepare_records(infile, seqtype, directory=None):
    changes = ChangeLog(directory, seqtype == 'protein')
    try:
        matrix = AlignmentSpool(infile, directory=directory, transform=lambda identifier, chunk, offset:
                                normalize_sequence(chunk, seqtype, infile, identifier, changes, offset))
        return matrix, changes
    except BaseException:
        changes.close()
        raise


def write_nexus(matrix, handle, seqtype):
    datatype = {'DNA': 'dna', 'RNA': 'rna', 'protein': 'protein'}[seqtype]
    handle.write('#NEXUS\n\nbegin data;\n')
    handle.write(f'  dimensions ntax={len(matrix.records)} nchar={matrix.length};\n')
    handle.write(f'  format datatype={datatype} missing=? gap=- interleave;\n  matrix\n')
    labels = [nexus_name(identifier) for identifier, _, _ in matrix.records]
    label_width = max(map(len, labels))
    positions = [record[1] for record in matrix.records]
    # Reuse a small page per taxon instead of seeking disk for every 70 bases.
    page_width = _LINE_WIDTH * 64
    for start in range(0, matrix.length, page_width):
        size = min(page_width, matrix.length - start)
        pages = []
        for index in range(len(labels)):
            matrix.file.seek(positions[index])
            pages.append(matrix.file.read(size))
            positions[index] = matrix.file.tell()
        for offset in range(0, size, _LINE_WIDTH):
            for label, page in zip(labels, pages):
                handle.write(f'{label:<{label_width}} {page[offset:offset + _LINE_WIDTH]}\n')
            handle.write('\n')
    handle.write('  ;\nend;\n')


def fasta2nex(infile, outfile=None, seqtype='DNA', changes_report=None):
    """Convert an aligned FASTA file to NEXUS using bounded memory."""
    check_output_path(outfile, inputs=[infile])
    matrix, changes = prepare_records(infile, seqtype, Path(outfile).absolute().parent if outfile else None)
    try:
        if changes:
            if changes_report is None:
                if outfile is None:
                    raise ValueError('Use --changes-report when normalizing input while writing NEXUS to stdout')
                changes_report = str(outfile) + '.sequence_changes.tsv'
            check_output_path(changes_report, inputs=[infile, outfile])
        with atomic_output(outfile, inputs=[infile]) as handle:
            write_nexus(matrix, handle, seqtype)
        if changes:
            changes.publish(changes_report)
    finally:
        matrix.close()
        changes.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a FASTA alignment to NEXUS using a streaming writer.', add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  fasta2nex.py -i alignment.fasta -o alignment.nex -st DNA')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA', required=True, help='Input FASTA alignment.')
    optional.add_argument('-o', '--output', metavar='NEXUS', help='Output NEXUS file (default: stdout).', default=None)
    optional.add_argument('-st', '--seqtype', '--seq-type', default='DNA', choices=['DNA', 'RNA', 'protein'], help='Alignment type: DNA normalizes U/IUPAC ambiguity to T/N; protein masks ambiguity as X (default: DNA).')
    optional.add_argument('--changes-report', metavar='TSV', help='Audit TSV for normalization edits; defaults to OUTPUT.sequence_changes.tsv when edits occur.')
    optional.add_argument('-h', '--help', action='help', help='Show program help and exit.')
    optional.add_argument('-v', '--version', action='version', version='v1.0.0', help='Show version and exit.')
    args = parser.parse_args()
    fasta2nex(args.input, args.output, args.seqtype, args.changes_report)
