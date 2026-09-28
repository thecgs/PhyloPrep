#!/usr/bin/env python
# coding: utf-8

import argparse, sys
import io
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from msap_io import atomic_output, read_fasta_alignment, check_output_path
from sequence_audit import event, normalize_nucleotides, write_report


def fasta2nex(infile, outfile=None, seqtype="DNA", changes_report=None):
    records = read_fasta_alignment(infile)
    changes = []
    for record in records:
        if seqtype == 'DNA':
            source = str(record.seq)
            sequence = normalize_nucleotides(source)
            for index, (original, replacement) in enumerate(zip(source, sequence), 1):
                if original.upper() != replacement:
                    changes.append(event(infile, record.id,
                                         'uracil_to_thymine' if original.upper() == 'U' else 'nucleotide_ambiguity',
                                         'NA', source_nt_positions=f'{index}-{index}',
                                         original=original, replacement=replacement, action='replace'))
            record.seq = Seq(sequence)
        elif seqtype == 'RNA':
            source = str(record.seq)
            sequence = normalize_nucleotides(source).replace('T', 'U')
            for index, (original, replacement) in enumerate(zip(source, sequence), 1):
                if original.upper() != replacement:
                    changes.append(event(infile, record.id,
                                         'thymine_to_uracil' if original.upper() == 'T' else 'nucleotide_ambiguity',
                                         'NA', source_nt_positions=f'{index}-{index}',
                                         original=original, replacement=replacement, action='replace'))
            record.seq = Seq(sequence)
        else:
            source = str(record.seq)
            sequence = ''.join(char.upper() if char.upper() in 'ACDEFGHIKLMNPQRSTVWYX-' else 'X'
                               for char in source)
            for index, (original, replacement) in enumerate(zip(source, sequence), 1):
                if original.upper() != replacement:
                    changes.append(event(infile, record.id, 'protein_unsupported_residue', 'NA',
                                         source_aa_positions=f'{index}-{index}', original=original,
                                         replacement='X', action='replace'))
            record.seq = Seq(sequence)
        record.annotations["molecule_type"] = seqtype
    buffer = io.StringIO()
    AlignIO.write(MultipleSeqAlignment(records), buffer, "nexus")
    text = buffer.getvalue()
    try:
        restored = list(AlignIO.read(io.StringIO(text), 'nexus'))
    except Exception as error:
        raise ValueError('NEXUS round-trip validation failed; unsupported characters must be normalized before export') from error
    if [(r.id, str(r.seq).upper()) for r in restored] != [(r.id, str(r.seq).upper()) for r in records]:
        raise ValueError('NEXUS round-trip changed taxon IDs or sequences')
    if changes:
        if changes_report is None:
            if outfile is None:
                raise ValueError('Use --changes-report when normalizing input while writing NEXUS to stdout')
            changes_report = str(outfile) + '.sequence_changes.tsv'
        check_output_path(changes_report, inputs=[infile, outfile])
    with atomic_output(outfile, inputs=[infile]) as handle:
        handle.write(text)
    if changes:
        write_report(changes_report, changes, protein=seqtype == 'protein')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a FASTA alignment to NEXUS format.", add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  fasta2nex.py -i alignment.fasta -o alignment.nex -st DNA")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA', required=True, help='Input FASTA alignment.')
    optional.add_argument('-o', '--output', metavar='NEXUS', help='Output NEXUS file (default: stdout).', default=None)
    optional.add_argument('-st', '--seqtype', '--seq-type', default='DNA', choices=["DNA", "RNA", "protein"], help='Alignment type: DNA normalizes U/IUPAC ambiguity to T/N; protein masks ambiguity as X (default: DNA).')
    optional.add_argument('--changes-report', metavar='TSV', help='Audit TSV for normalization edits; defaults to OUTPUT.sequence_changes.tsv when edits occur.')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0', help="Show program's version number and exit.")
    args = parser.parse_args()

    fasta2nex(args.input, args.output, args.seqtype, args.changes_report)
