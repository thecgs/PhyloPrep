"""Change provenance, true MACSE export semantics, and workflow boundaries."""
import csv
import shutil
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import pytest
from Bio import AlignIO, SeqIO

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from sequence_audit import FIELDS, preprocess_codon, macse_align
from sequence_audit import event, write_report
from alignment_qc import alignment_qc
from phyloprep import relabel_change_report


def rows(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def cli(script, *args, cwd, env=None):
    return subprocess.run([sys.executable, str(REPO / script), *map(str, args)],
                          cwd=cwd, env=env, capture_output=True, text=True)


@pytest.mark.parametrize('code,expected,events', [
    (1, 'ATGNNNGCT', [('internal_stop', '4-6', 'tga', 'NNN'), ('terminal_stop', '10-12', 'taa', '-')]),
    (2, 'ATGTGAGCT', [('terminal_stop', '10-12', 'taa', '-')]),
])
def test_codon_exact_events_and_genetic_code(tmp_path, code, expected, events):
    source, output, report = [tmp_path / name for name in ('gene.fa', 'clean.fa', 'gene.sequence_changes.tsv')]
    source.write_text('>a\natgtgagcttaa\n>b\nATGGCTGCTGCT\n')
    preprocess_codon(source, output, code, report)
    assert str(next(SeqIO.parse(output, 'fasta')).seq) == expected
    assert [(r['event'], r['source_nt_positions'], r['original'], r['replacement']) for r in rows(report)] == events
    assert all(r['taxa_id'] == 'a' and r['input_file'] == str(source) for r in rows(report))
    assert source.read_text().startswith('>a\natgtga')


def test_partial_terminal_codon_is_trimmed_and_audited(tmp_path, capsys):
    source, output, report = [tmp_path / name for name in ('gene.fa', 'clean.fa', 'gene.sequence_changes.tsv')]
    source.write_text('>a\nATGG\n>b\nATGAA\n')
    output.write_text('previous output')
    preprocess_codon(source, output, 1, report)
    assert {record.id: str(record.seq) for record in SeqIO.parse(output, 'fasta')} == {'a': 'ATG', 'b': 'ATG'}
    findings = rows(report)
    assert [(row['taxa_id'], row['event'], row['source_nt_positions'], row['original'], row['replacement'], row['action'])
            for row in findings] == [
        ('a', 'terminal_partial_codon', '4-4', 'G', '-', 'delete'),
        ('b', 'terminal_partial_codon', '4-5', 'AA', '-', 'delete'),
    ]
    assert 'WARNING: a: CDS length 4 is not divisible by 3; trimmed 1 trailing nucleotide(s)' in capsys.readouterr().err


def test_stops_detected_despite_gaps(tmp_path):
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nATG---TGAGCTTAA---\n>b\nATGGCTGCTGCTGCT---\n')
    row = alignment_qc(str(source), 'codon')['sequences'][0]
    assert row['terminal_stop'] == 'TAA'
    assert row['internal_stop_count'] == '3-TGA'


def test_batch_trims_terminal_partial_codons_and_keeps_report(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'broken.fa'
    source.write_text('>a\nATGG\n>b\nATGAA\n')
    result = cli('MSAP_batch.py', '-i', source, '-o', tmp_path/'out', '-t', '1', cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr
    report = rows(tmp_path/'out/reports/sequence_changes.tsv')
    assert sum(r['event'] == 'terminal_partial_codon' for r in report) == 2
    assert (tmp_path/'out/broken.mafft.codon.aln').exists()


def test_pipeline_records_only_prepared_file_and_taxa_in_failure_report(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'broken.fa'
    source.write_text('>original_a\nATGG\n>original_b\nATGAA\n')
    mapping = tmp_path / 'taxa.tsv'
    mapping.write_text('original_a\tSpecies_A\noriginal_b\tSpecies_B\n')
    result = cli('phyloprep.py', '-i', source, '-m', mapping, '-o', tmp_path/'out', '-t', '1', cwd=tmp_path, env=env)
    assert result.returncode != 0
    report = [r for r in rows(tmp_path/'out/reports/sequence_changes.tsv') if r['event'] == 'terminal_partial_codon']
    assert {r['taxa_id'] for r in report} == {'Species_A', 'Species_B'}
    assert all(r['input_file'] == str(tmp_path/'out/renamed/broken.fa') for r in report)
    assert all(list(r) == FIELDS for r in report)


def test_pipeline_maps_special_ids_and_skips_terminal_stop_only_locus(tmp_path, fake_aligner):
    _, env = fake_aligner
    skipped, retained = tmp_path / 'skipped.fa', tmp_path / 'retained.fa'
    skipped.write_text('>original:a\nTAA\n>original,b\nTAA\n')
    retained.write_text('>original:a\nGCT\n>original,b\nGCC\n')
    mapping = tmp_path / 'taxa.tsv'
    mapping.write_text('original:a\tSpecies_A\noriginal,b\tSpecies_B\n')
    output = tmp_path / 'out'
    result = cli('phyloprep.py', '-i', skipped, retained, '-m', mapping, '-o', output,
                 '-t', '1', '--notrim', cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'WARNING: original:a, original,b: no coding codons remain' in result.stderr
    report = [r for r in rows(output/'reports/sequence_changes.tsv') if r['event'] == 'terminal_stop']
    assert {r['taxa_id'] for r in report} == {'Species_A', 'Species_B'}
    assert all(r['input_file'] == str(output/'renamed/skipped.fa') for r in report)
    assert not (output/'alignments/skipped.mafft.codon.aln').exists()
    assert {r.id for r in SeqIO.parse(output/'matrices/codon/raw/supermatrix.fasta', 'fasta')} == {'Species_A', 'Species_B'}


def test_relabel_change_report_maps_macse_native_audit_rows(tmp_path):
    source, target = tmp_path/'source.fa', tmp_path/'renamed.fa'
    native, report = tmp_path/'gene.macse_original.NT.fasta', tmp_path/'changes.tsv'
    source.write_text('>original:a\nATG\n')
    target.write_text('>Species_A\nATG\n')
    native.write_text('>original:a\nATG\n')
    write_report(report, [event(native, 'original:a', 'internal_frameshift', 1, action='replace')])
    relabel_change_report(report, {
        str(source.resolve()): (target.resolve(), {'original_a': 'Species_A'}),
    })
    finding = rows(report)[0]
    assert finding['taxa_id'] == 'Species_A'
    # Native MACSE alignment columns remain the provenance coordinate system.
    assert finding['input_file'] == str(native.resolve())


@pytest.mark.skipif(not shutil.which('java') or not (REPO/'macse_v2.07.jar').exists(), reason='local MACSE/Java unavailable')
def test_macse_export_reports_exact_changes_and_placeholder_coordinates(tmp_path, monkeypatch):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nATGTGAGCTAA\n>b\nATGGCTGCTGCTTAA\n')
    native_nt = '>a\nATGTGA---GC!TAA\n>b\nATGGCTGCTGCTTAA\n'
    real_run = subprocess.run
    def controlled_alignment(command, **kwargs):
        if 'alignSequences' in command:
            Path(command[command.index('-out_NT')+1]).write_text(native_nt)
            Path(command[command.index('-out_AA')+1]).write_text('>a\nM*-!*\n>b\nMAAA*\n')
            return subprocess.CompletedProcess(command, 0)
        return real_run(command, **kwargs)  # Actual MACSE 2.07 export, not a simulated replacement.
    monkeypatch.setattr('sequence_audit.subprocess.run', controlled_alignment)
    args = Namespace(macse_memory='512m', macse_jar=str(REPO/'macse_v2.07.jar'), genetic_code=1)
    prefix = str(tmp_path/'gene')
    nt, aa = macse_align(args, str(source), prefix, shutil.which('java'))
    assert {r.id: str(r.seq) for r in SeqIO.parse(nt, 'fasta')} == {'a':'ATGNNN---NNNNNN', 'b':'ATGGCTGCTGCTNNN'}
    assert '!' not in Path(nt).read_text()
    changes = rows(prefix+'.sequence_changes.tsv')
    frame = next(r for r in changes if r['event'] == 'internal_frameshift')
    assert frame['source_nt_positions'] == '10-12'
    assert frame['input_file'] == prefix + '.macse_original.NT.fasta'
    assert list(frame) == ['input_file', 'taxa_id', 'event', 'source_nt_positions', 'original', 'replacement', 'action', 'genetic_code']
    assert (frame['original'], frame['replacement']) == ('GC!', 'NNN')
    assert {r['event'] for r in changes} == {'internal_frameshift', 'terminal_stop', 'internal_stop'}
    assert not Path(prefix+'.sequence_events.tsv').exists()
    assert not Path(prefix+'.macse.source_positions.tsv').exists()
    assert Path(prefix+'.macse_original.NT.fasta').read_text() == native_nt
    assert len(AlignIO.read(nt, 'fasta')[0].seq) == 3*len(AlignIO.read(aa, 'fasta')[0].seq)


@pytest.mark.skipif(not shutil.which('java') or not (REPO/'macse_v2.07.jar').exists(), reason='local MACSE/Java unavailable')
def test_real_pseudogene_pipeline_resume_and_default_products(tmp_path):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nATGGCTGGTGCTGGTTAA\n>b\nATGGCCGGAGCCGGATAA\n>c\nATGGCGGGCAGCGGGCTAA\n')
    out = tmp_path/'out'
    arguments = ['-i', source, '-st', 'pseudogene', '-o', out, '-t', '1', '--macse-memory', '512m']
    result = cli('phyloprep.py', *arguments, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    for variant in ('raw', 'trimmed'):
        for kind, datatype in [('codon', 'dna'), ('prot', 'protein')]:
            folder = out/'matrices/pseudogene'/kind/variant
            for extension, fmt in [('fasta','fasta'),('phy','phylip-relaxed'),('nex','nexus')]:
                alignment = AlignIO.read(folder/f'supermatrix.{extension}', fmt)
                assert len(alignment) == 3
            assert f'datatype={datatype}' in (folder/'supermatrix.nex').read_text().lower()
        nt = out/'matrices/pseudogene/codon'/variant
        for stem in ('codon1st.supermatrix', 'codon2nd.supermatrix', 'codon3rd.supermatrix', 'fourfold'):
            for extension in ('fasta', 'phy', 'nex'):
                assert (nt/f'{stem}.{extension}').exists()
    assert 'do not establish' in result.stdout
    for report in (out/'reports/sequence_changes.tsv', out/'alignments/reports/sequence_changes.tsv',
                   out/'alignments/gene.sequence_changes.tsv'):
        changes = rows(report)
        assert changes
        for row in changes:
            assert row['input_file'] == str(out/'alignments/gene.macse_original.NT.fasta')
            native = {r.id: str(r.seq) for r in SeqIO.parse(row['input_file'], 'fasta')}
            start, end = map(int, row['source_nt_positions'].split('-'))
            assert native[row['taxa_id']][start-1:end] == row['original']
    for pattern in ('*.sequence_events.tsv', 'sequence_events.tsv', '*.source_positions.tsv', '*.columns.tsv', 'taxa_provenance.tsv'):
        assert not list(out.rglob(pattern))
    # The default rerun retains derived matrices and resumes the expensive alignment.
    result = cli('phyloprep.py', *arguments, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'skipped (already complete)' in result.stdout
    assert list((out/'matrices').rglob('codon1st*'))
    assert list((out/'matrices').rglob('fourfold*'))
    assert not list((out/'matrices/pseudogene/prot').rglob('codon*'))


def test_explicit_incompatible_aligner_is_rejected(tmp_path):
    result = cli('MSAP.py', '-i', 'unused.fa', '-st', 'pseudogene', '-s', 'mafft', cwd=tmp_path)
    assert result.returncode != 0
    assert 'pseudogene mode requires MACSE' in result.stderr


def test_ambiguous_terminal_stop_retains_previous_codon_behavior(tmp_path):
    source, output, report = [tmp_path/name for name in ('gene.fa', 'clean.fa', 'gene.sequence_changes.tsv')]
    source.write_text('>a\nATGGCTTAR\n>b\nATGGCTGCT\n')
    preprocess_codon(source, output, 1, report)
    assert str(next(SeqIO.parse(output, 'fasta')).seq) == 'ATGGCT'
    assert [(r['event'], r['source_nt_positions'], r['original']) for r in rows(report)] == [('terminal_stop','7-9','TAR')]


def test_iupac_nucleotides_are_masked_as_n_with_audit_rows(tmp_path):
    source, output, report = [tmp_path/name for name in ('gene.fa', 'clean.fa', 'gene.sequence_changes.tsv')]
    source.write_text('>a\nATGRYS\n>b\nATGGCT\n')
    preprocess_codon(source, output, 1, report)
    assert '>a\nATGNNN\n' in output.read_text()
    changes = rows(report)
    assert [(row['source_nt_positions'], row['original'], row['replacement']) for row in changes] == [
        ('4-4', 'R', 'N'), ('5-5', 'Y', 'N'), ('6-6', 'S', 'N')]


def test_trimming_preserves_result_without_column_report(tmp_path):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nATGNNNGCT\n>b\nATGNNNGCC\n')
    output = tmp_path/'trim.fa'
    result = cli('trimAlnSeq.py', '-i', source, '-o', output, '-st', 'codon', cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert {r.id: str(r.seq) for r in SeqIO.parse(output, 'fasta')} == {'a': 'ATGGCT', 'b': 'ATGGCC'}
    assert not list(tmp_path.glob('*.columns.tsv'))


@pytest.mark.skipif(not shutil.which('java') or not (REPO/'macse_v2.07.jar').exists(), reason='local MACSE/Java unavailable')
def test_pseudogene_empty_fourfold_is_explained_and_skipped(tmp_path):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nATGATG\n>b\nATGATG\n')
    out = tmp_path/'out'
    result = cli('phyloprep.py', '-i', source, '-st', 'pseudogene', '--notrim', '-o', out, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    folder = out/'matrices/pseudogene/codon/raw'
    assert 'no usable four-fold sites; skipping FASTA/PHY/NEX' in result.stdout
    assert (folder/'fourfold.skipped.txt').exists()
    assert not (folder/'fourfold.fasta').exists()
    assert (folder/'codon3rd.supermatrix.nex').exists()
    assert (out/'matrices/pseudogene/prot/raw/supermatrix.phy').exists()


@pytest.mark.parametrize('script', ['MSAP.py', 'MSAP_batch.py', 'phyloprep.py'])
def test_help_includes_pseudogene_example_without_optional_switches(tmp_path, script):
    result = cli(script, '--help', cwd=tmp_path)
    assert result.returncode == 0
    assert '-st pseudogene --macse-jar macse_v2.07.jar' in result.stdout
    assert '--split-codons' not in result.stdout
    assert '--extract-fourfold' not in result.stdout
