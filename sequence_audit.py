"""Shared sequence-change provenance and MACSE workflow support.

Coordinates are 1-based and inclusive within input_file. Codon preprocessing
uses input sequence positions; MACSE export uses native NT alignment columns.
Protein edits use source_aa_positions; source_nt_positions is NA for these rows.
"""
from pathlib import Path
import csv
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.Data import CodonTable
from msap_io import atomic_output, read_fasta_alignment, is_stop_codon, normalize_dna

FIELDS = ['input_file', 'taxa_id', 'event', 'source_nt_positions',
          'original', 'replacement', 'action', 'genetic_code']


def write_report(path, rows, protein=False):
    rows = list(rows)
    fields = FIELDS + (['source_aa_positions'] if protein or any('source_aa_positions' in r for r in rows) else [])
    with atomic_output(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', restval='NA')
        writer.writeheader()
        writer.writerows(rows)


def event(infile, taxon, kind, code, **values):
    row = dict.fromkeys(FIELDS, 'NA')
    row.update(input_file=str(Path(infile).resolve()), taxa_id=taxon,
               event=kind, genetic_code=code)
    row.update(values)
    return row


def rejected_input(infile, prefix, mode, code, reason):
    rows = []
    if Path(infile).is_file():
        for record in SeqIO.parse(infile, 'fasta'):
            if mode == 'codon' and len(record.seq) % 3:
                rows.append(event(infile, record.id, 'invalid_frame', code, action='reject'))
    rows.append(event(infile, 'NA', 'validation_error', code, action='reject'))
    write_report(prefix + '.sequence_changes.tsv', rows)


def preprocess_codon(infile, outfile, code, report):
    records = list(SeqIO.parse(infile, 'fasta'))
    rows, outputs = [], []
    for record in records:
        original_sequence = str(record.seq)
        remainder = len(original_sequence) % 3
        sequence = original_sequence
        if remainder:
            start = len(original_sequence) - remainder
            removed = original_sequence[start:]
            rows.append(event(infile, record.id, 'terminal_partial_codon', code,
                              source_nt_positions=f'{start + 1}-{len(original_sequence)}',
                              original=removed, replacement='-', action='delete'))
            print(
                f'WARNING: {record.id}: CDS length {len(original_sequence)} is not divisible by 3; '
                f'trimmed {remainder} trailing nucleotide(s) at positions {start + 1}-{len(original_sequence)}. '
                f'See {report}.', file=sys.stderr,
            )
            sequence = original_sequence[:start]
        if not sequence:
            rows.append(event(infile, record.id, 'validation_error', code, action='reject'))
            write_report(report, rows)
            raise ValueError(f'{record.id}: CDS has no complete codon after trailing-base trimming')
        fragments = []
        for i in range(0, len(sequence), 3):
            # Detect a definite stop before collapsing IUPAC ambiguity.  For
            # example TAR is certainly a stop under the standard code, while
            # TAN is not; this preserves the established stop policy.
            source_codon = normalize_dna(sequence[i:i+3])
            codon = normalize_nucleotides(source_codon)
            replacement = codon
            terminal = i + 3 == len(sequence)
            remove_terminal = terminal and is_stop_codon(source_codon, code, terminal=True)
            if is_stop_codon(source_codon, code, terminal=terminal):
                replacement = '' if remove_terminal else 'NNN'
                rows.append(event(infile, record.id,
                                  'terminal_stop' if terminal else 'internal_stop', code,
                                  source_nt_positions=f'{i+1}-{i+3}', original=original_sequence[i:i+3],
                                  replacement=replacement or '-', action='delete' if remove_terminal else 'replace'))
            else:
                rows.extend(nucleotide_changes(infile, record.id, original_sequence[i:i+3], code, offset=i))
            fragments.append(replacement)
        outputs.append((record.id, ''.join(fragments)))
    with atomic_output(outfile, inputs=[infile]) as handle:
        for identifier, seq in outputs:
            handle.write(f'>{identifier}\n{seq}\n')
    write_report(report, rows)


def rna_changes(infile, taxon, sequence, code='NA', offset=0):
    return [event(infile, taxon, 'uracil_to_thymine', code,
                  source_nt_positions=f'{offset+i+1}-{offset+i+1}',
                  original=base, replacement='T', action='replace')
            for i, base in enumerate(sequence) if base.upper() == 'U']


def nucleotide_changes(infile, taxon, sequence, code='NA', offset=0):
    """Audit conversion to the canonical DNA alphabet A/C/G/T/N."""
    changes = rna_changes(infile, taxon, sequence, code, offset)
    for i, base in enumerate(sequence):
        if base.upper() in 'RYSWKMBDHV?':
            changes.append(event(infile, taxon, 'nucleotide_ambiguity', code,
                                 source_nt_positions=f'{offset+i+1}-{offset+i+1}',
                                 original=base, replacement='N', action='replace'))
    return changes


def normalize_nucleotides(sequence):
    """Return DNA using only A/C/G/T/N and gaps; U is normalized to T."""
    return normalize_dna(sequence).translate(str.maketrans({base: 'N' for base in 'RYSWKMBDHV?'}))


def normalize_export_records(records, infile, seqtype):
    """Canonicalize aligned records for an interchange format and audit edits."""
    rows = []
    protein = seqtype.lower() in {'protein', 'prot', 'aa'}
    for record in records:
        source = str(record.seq)
        if protein:
            normalized = ''.join(c.upper() if c.upper() in 'ACDEFGHIKLMNPQRSTVWYX-' else 'X'
                                 for c in source)
        else:
            normalized = normalize_nucleotides(source)
        for index, (original, replacement) in enumerate(zip(source, normalized), 1):
            if original.upper() != replacement:
                kind = ('protein_unsupported_residue' if protein else
                        ('uracil_to_thymine' if original.upper() == 'U' else 'nucleotide_ambiguity'))
                rows.append(event(infile, record.id, kind, 'NA',
                                  source_aa_positions=f'{index}-{index}' if protein else 'NA',
                                  source_nt_positions='NA' if protein else f'{index}-{index}',
                                  original=original, replacement=replacement, action='replace'))
        record.seq = record.seq.__class__(normalized)
    return rows, protein


def preprocess_nucl(infile, outfile, report):
    rows = []
    with atomic_output(outfile, inputs=[infile]) as handle:
        for record in SeqIO.parse(infile, 'fasta'):
            rows.extend(nucleotide_changes(infile, record.id, str(record.seq)))
            handle.write(f'>{record.id}\n{normalize_nucleotides(record.seq)}\n')
    write_report(report, rows)


def preprocess_protein(infile, outfile, report, stop_symbol='*'):
    """Mask internal stops and remove trailing stops, preserving source coordinates."""
    rows, outputs = [], []
    for record in SeqIO.parse(infile, 'fasta'):
        source = str(record.seq)
        end = len(source)
        while end and source[end-1] in (stop_symbol, '-'):
            end -= 1
        if not end:
            write_report(report, [event(infile, record.id, 'validation_error', 'NA',
                                       action='reject', source_aa_positions='NA')], protein=True)
            raise ValueError(f'{record.id}: protein has no residues after terminal-stop removal')
        cleaned = []
        for i, char in enumerate(source):
            if char == stop_symbol:
                terminal = i >= end
                replacement = '' if terminal else 'X'
                rows.append(event(infile, record.id,
                                  'protein_terminal_stop' if terminal else 'protein_internal_stop', 'NA',
                                  source_aa_positions=f'{i+1}-{i+1}', original=char,
                                  replacement=replacement or '-', action='delete' if terminal else 'replace'))
                cleaned.append(replacement)
            elif char.upper() not in 'ACDEFGHIKLMNPQRSTVWYX-':
                rows.append(event(infile, record.id, 'protein_unsupported_residue', 'NA',
                                  source_aa_positions=f'{i+1}-{i+1}', original=char,
                                  replacement='X', action='replace'))
                cleaned.append('X')
            else:
                cleaned.append(char.upper())
        outputs.append((record.id, ''.join(cleaned)))
    with atomic_output(outfile, inputs=[infile]) as handle:
        for identifier, sequence in outputs:
            handle.write(f'>{identifier}\n{sequence}\n')
    write_report(report, rows, protein=True)
    unsupported = sum(row['event'] == 'protein_unsupported_residue' for row in rows)
    if unsupported:
        print(f'prot: masked {unsupported} ambiguous/non-standard residue(s) as X for alignment and format compatibility; see {report}', file=sys.stderr)


def add_protein_options(group):
    group.add_argument('--protein-stop-symbol', '--stop-symbol', default='*', metavar='CHAR',
                       help="prot only: stop character (default: '*'); trailing stops are removed, internal stops become X. Use '.' for dot-encoded stops.")


def add_macse_options(group):
    group.add_argument('--macse-jar', default=str(Path(__file__).with_name('macse_v2.07.jar')),
                       help='MACSE jar used by pseudogene mode (default: bundled macse_v2.07.jar).')
    group.add_argument('--macse-memory', default='2g', help='Java maximum heap per batch job (default: 2g).')


def normalize_workflow(args):
    symbol = getattr(args, 'protein_stop_symbol', '*')
    if len(symbol) != 1 or not symbol.isprintable() or symbol.isspace() or symbol.isalnum() or symbol in '- >':
        raise ValueError('--protein-stop-symbol must be one non-letter, non-digit symbol other than a gap or FASTA delimiter')
    if args.seqtype != 'prot' and symbol != '*':
        raise ValueError('--protein-stop-symbol is only supported with -st prot')
    if args.seqtype == 'pseudogene':
        if args.align_software not in (None, 'macse'):
            raise ValueError('pseudogene mode requires MACSE; remove -s or use -s macse')
        args.align_software = 'macse'
        args.macse_jar = str(Path(args.macse_jar).resolve())
        if not Path(args.macse_jar).is_file():
            raise ValueError(f'MACSE jar does not exist: {args.macse_jar}')
        if not re.fullmatch(r'[1-9][0-9]*[mMgG]', args.macse_memory):
            raise ValueError('--macse-memory must be a positive heap size such as 2g or 512m')
    else:
        if args.align_software == 'macse':
            raise ValueError('-s macse requires -st pseudogene')
        args.align_software = args.align_software or 'mafft'


def macse_align(args, infile, prefix, java):
    original_nt = prefix + '.macse_original.NT.fasta'
    original_aa = prefix + '.macse_original.AA.fasta'
    nt, aa = prefix + '.macse.codon.aln', prefix + '.macse.prot.aln'
    base = [java, '-Xmx' + args.macse_memory, '-jar', args.macse_jar]
    subprocess.run([*base, '-prog', 'alignSequences', '-seq', infile,
                    '-gc_def', str(args.genetic_code), '-out_NT', original_nt, '-out_AA', original_aa], check=True)
    originals = read_fasta_alignment(original_nt, codon=True)
    sources = {r.id: str(r.seq).upper() for r in SeqIO.parse(infile, 'fasta')}
    if set(sources) != {r.id for r in originals}:
        raise ValueError('MACSE changed the input taxon set')
    kinds = {}
    for record in originals:
        seq = str(record.seq).upper()
        if ''.join(c for c in seq if c not in '-!') != sources[record.id]:
            raise ValueError(f'{record.id}: cannot map MACSE alignment to input sequence')
        occupied = [i for i in range(0, len(seq), 3) if seq[i:i+3] != '---']
        for i in occupied:
            codon = seq[i:i+3]
            if '!' in codon:
                kind = 'external_frameshift' if i in (occupied[0], occupied[-1]) else 'internal_frameshift'
            elif is_stop_codon(codon, args.genetic_code):
                kind = 'terminal_stop' if i == occupied[-1] else 'internal_stop'
            else:
                continue
            kinds[(record.id, i)] = kind
    # Keep columns fixed so native, standardized NT and AA coordinates correspond exactly.
    subprocess.run([*base, '-prog', 'exportAlignment', '-align', original_nt,
                    '-gc_def', str(args.genetic_code), '-codonForInternalStop', 'NNN',
                    '-codonForFinalStop', 'NNN', '-codonForInternalFS', 'NNN',
                    '-codonForExternalFS', 'NNN', '-charForRemainingFS', 'N',
                    '-keep_gap_only_sites_ON', '-out_NT', nt, '-out_AA', aa], check=True)
    exported = {r.id: normalize_nucleotides(str(r.seq)) for r in read_fasta_alignment(nt, codon=True)}
    proteins = {r.id: str(r.seq).upper() for r in read_fasta_alignment(aa)}
    changes = []
    canonical_proteins = {}
    for record in originals:
        before, after = str(record.seq).upper(), exported[record.id]
        protein = list(proteins[record.id])
        # A non-standard MACSE amino-acid symbol cannot safely be represented
        # in all downstream formats.  Mask its codon as unknown as well, so
        # the paired NT and AA matrices carry the same uncertainty.
        for aa_index, residue in enumerate(protein):
            if residue not in 'ACDEFGHIKLMNPQRSTVWYX-':
                start = 3 * aa_index
                after = after[:start] + 'NNN' + after[start + 3:]
                protein[aa_index] = 'X'
                changes.append(event(original_nt, record.id, 'protein_unsupported_residue', args.genetic_code,
                                     source_nt_positions=f'{start+1}-{start+3}', original=residue,
                                     replacement='X', action='replace',
                                     source_aa_positions=f'{aa_index+1}-{aa_index+1}'))
        exported[record.id] = after
        canonical_proteins[record.id] = ''.join(protein)
        if len(before) != len(after) or len(after) != 3 * len(proteins[record.id]):
            raise ValueError('MACSE export changed alignment dimensions')
        if set(after) - set('ACGTN-') or set(canonical_proteins[record.id]) - set('ACDEFGHIKLMNPQRSTVWYX-'):
            raise ValueError('MACSE export retains unsupported characters')
        for i in range(0, len(before), 3):
            if before[i:i+3] != after[i:i+3]:
                changes.append(event(original_nt, record.id, kinds.get((record.id, i), 'export_change'), args.genetic_code,
                                     source_nt_positions=f'{i+1}-{i+3}',
                                     original=before[i:i+3],
                                     replacement=after[i:i+3], action='replace'))
    # Replace MACSE exports only after their dimensions and symbols have been
    # checked, keeping the public matrices canonical and paired.
    with atomic_output(nt, inputs=[original_nt]) as handle:
        for record in originals:
            handle.write(f'>{record.id}\n{exported[record.id]}\n')
    with atomic_output(aa, inputs=[original_aa]) as handle:
        for record in originals:
            handle.write(f'>{record.id}\n{canonical_proteins[record.id]}\n')
    write_report(prefix + '.sequence_changes.tsv', changes)
    return nt, aa


def merge_reports(paths, output):
    rows = []
    protein = False
    for path in paths:
        if not Path(path).is_file():
            continue
        with open(path, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            has_aa = 'source_aa_positions' in (reader.fieldnames or [])
            protein |= has_aa
            for row in reader:
                rows.append({field: row.get(field, 'NA') for field in FIELDS + (['source_aa_positions'] if has_aa else [])})
    write_report(output, rows, protein=protein)


def relocate_change_report(report, source_dir, destination_dir):
    """Update references to native matrices before batch staging is promoted."""
    if not Path(report).is_file():
        return
    with open(report, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        protein = 'source_aa_positions' in (reader.fieldnames or [])
    for row in rows:
        path = Path(row['input_file'])
        if path.parent.resolve() == Path(source_dir).resolve():
            row['input_file'] = str((Path(destination_dir) / path.name).resolve())
    write_report(report, rows, protein=protein)
