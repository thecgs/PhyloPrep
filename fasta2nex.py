#!/usr/bin/env python
# coding: utf-8

"""Convert aligned FASTA to interleaved NEXUS without using Biopython.

The Biopython NEXUS reader builds interleaved alignments by repeatedly
concatenating short blocks. That makes its own round-trip check impractical
for multi-million-column supermatrices. This module writes a conventional
interleaved matrix from temporary per-taxon sequence spools instead.
"""

import argparse
import tempfile

from msap_io import atomic_output, check_output_path
from sequence_audit import event, normalize_nucleotides, write_report


_NEXUS_PUNCTUATION = set("()[]{}\\,;:=*'\"`+-<>")
_NEXUS_WHITESPACE = set(" \t\n")
_LINE_WIDTH = 70


def fasta_records(infile):
    """Yield ``(id, sequence)`` records with SeqIO-compatible FASTA IDs."""
    identifier = None
    sequence = []
    with open(infile, encoding='utf-8') as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if identifier is not None:
                    yield identifier, ''.join(sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError('FASTA record has an empty identifier')
                identifier = header.split()[0]
                sequence = []
            elif identifier is None:
                raise ValueError('FASTA sequence data appears before its first header')
            else:
                sequence.append(''.join(line.split()))
    if identifier is not None:
        yield identifier, ''.join(sequence)


def nexus_name(identifier):
    """Quote an identifier when required by the NEXUS token grammar."""
    escaped = identifier.replace("'", "''")
    if set(escaped).intersection(_NEXUS_PUNCTUATION | _NEXUS_WHITESPACE):
        return f"'{escaped}'"
    return escaped


def normalize_sequence(source, seqtype, infile, identifier, changes):
    if seqtype == 'DNA':
        sequence = normalize_nucleotides(source)
        for index, (original, replacement) in enumerate(zip(source, sequence), 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier,
                                     'uracil_to_thymine' if original.upper() == 'U' else 'nucleotide_ambiguity',
                                     'NA', source_nt_positions=f'{index}-{index}',
                                     original=original, replacement=replacement, action='replace'))
    elif seqtype == 'RNA':
        sequence = normalize_nucleotides(source).replace('T', 'U')
        for index, (original, replacement) in enumerate(zip(source, sequence), 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier,
                                     'thymine_to_uracil' if original.upper() == 'T' else 'nucleotide_ambiguity',
                                     'NA', source_nt_positions=f'{index}-{index}',
                                     original=original, replacement=replacement, action='replace'))
    else:
        sequence = ''.join(char.upper() if char.upper() in 'ACDEFGHIKLMNPQRSTVWYX-' else 'X'
                           for char in source)
        for index, (original, replacement) in enumerate(zip(source, sequence), 1):
            if original.upper() != replacement:
                changes.append(event(infile, identifier, 'protein_unsupported_residue', 'NA',
                                     source_aa_positions=f'{index}-{index}', original=original,
                                     replacement='X', action='replace'))
    return sequence


def prepare_records(infile, seqtype):
    """Validate FASTA and spool normalized sequences for interleaved output."""
    ids = set()
    changes = []
    records = []
    count = 0
    length = None
    for identifier, source in fasta_records(infile):
        if identifier in ids:
            raise ValueError(f'Duplicate sequence ID: {identifier}')
        ids.add(identifier)
        sequence = normalize_sequence(source, seqtype, infile, identifier, changes)
        if not sequence:
            raise ValueError(f'{identifier}: sequence is empty')
        if length is None:
            length = len(sequence)
        elif len(sequence) != length:
                raise ValueError(f'Alignment length mismatch: {identifier} has length {len(sequence)}, expected {length}')
        spool = tempfile.SpooledTemporaryFile(max_size=8 * 1024 * 1024, mode='w+t', encoding='utf-8')
        spool.write(sequence)
        spool.seek(0)
        records.append((identifier, spool))
        count += 1
    if not count:
        raise ValueError('No sequences found in input FASTA.')
    return records, count, length, changes


def write_nexus(records, handle, seqtype, count, length):
    """Write a standards-compliant, wrapped interleaved NEXUS matrix."""
    datatype = {'DNA': 'dna', 'RNA': 'rna', 'protein': 'protein'}[seqtype]
    handle.write('#NEXUS\n\nbegin data;\n')
    handle.write(f'  dimensions ntax={count} nchar={length};\n')
    handle.write(f'  format datatype={datatype} missing=? gap=- interleave;\n  matrix\n')
    labels = [(nexus_name(identifier), spool) for identifier, spool in records]
    label_width = max(len(label) for label, _ in labels)
    for start in range(0, length, _LINE_WIDTH):
        for label, spool in labels:
            handle.write(f'{label:<{label_width}} {spool.read(_LINE_WIDTH)}\n')
        handle.write('\n')
    handle.write('  ;\nend;\n')


def fasta2nex(infile, outfile=None, seqtype='DNA', changes_report=None):
    """Convert an aligned FASTA file to NEXUS using bounded memory."""
    check_output_path(outfile, inputs=[infile])
    records, count, length, changes = prepare_records(infile, seqtype)
    if changes:
        if changes_report is None:
            if outfile is None:
                raise ValueError('Use --changes-report when normalizing input while writing NEXUS to stdout')
            changes_report = str(outfile) + '.sequence_changes.tsv'
        check_output_path(changes_report, inputs=[infile, outfile])
    try:
        with atomic_output(outfile, inputs=[infile]) as handle:
            write_nexus(records, handle, seqtype, count, length)
        if changes:
            write_report(changes_report, changes, protein=seqtype == 'protein')
    finally:
        for _, spool in records:
            spool.close()


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
