import csv
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sequence_audit import preprocess_protein, merge_reports


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


@pytest.mark.parametrize('symbol', ['*', '.'])
def test_protein_stop_edits_and_exact_coordinates(tmp_path, symbol):
    source, output, report = [tmp_path/n for n in ('input.fa', 'clean.fa', 'changes.tsv')]
    original = f'>a\nmA{symbol}G{symbol}{symbol}--\n>b\nM-AAG\n'
    source.write_text(original)
    preprocess_protein(source, output, report, symbol)
    assert {r.id: str(r.seq) for r in SeqIO.parse(output, 'fasta')} == {'a': 'MAXG--', 'b': 'M-AAG'}
    edits = rows(report)
    assert [(r['source_aa_positions'], r['action'], r['replacement']) for r in edits] == [
        ('3-3', 'replace', 'X'), ('5-5', 'delete', '-'), ('6-6', 'delete', '-')]
    sources = {r.id: str(r.seq) for r in SeqIO.parse(source, 'fasta')}
    for row in edits:
        start, end = map(int, row['source_aa_positions'].split('-'))
        assert sources[row['taxa_id']][start-1:end] == row['original'] == symbol
        assert row['source_nt_positions'] == row['genetic_code'] == 'NA'
        assert row['input_file'] == str(source)
    assert source.read_text() == original


@pytest.mark.parametrize('sequence', ['***', '*--*', '---'])
def test_empty_protein_after_cleaning_rejected(tmp_path, sequence):
    source, output, report = [tmp_path/n for n in ('input.fa', 'clean.fa', 'changes.tsv')]
    source.write_text(f'>a\n{sequence}\n>b\nMAAG\n')
    output.write_text('previous result')
    with pytest.raises(ValueError, match='no residues'):
        preprocess_protein(source, output, report)
    assert output.read_text() == 'previous result'
    assert rows(report)[0]['action'] == 'reject'


@pytest.mark.parametrize('script', ['MSAP.py', 'MSAP_batch.py', 'phyloprep.py'])
def test_custom_stop_through_all_entrypoints(tmp_path, fake_aligner, script):
    _, env = fake_aligner
    source = tmp_path/'protein.fa'
    source.write_text('>a\nMA.G..\n>b\nMAAG.\n')
    output = tmp_path/'out'
    args = [sys.executable, str(ROOT/script), '-i', str(source), '-st', 'prot',
            '--stop-symbol', '.', '--notrim', '-t', '1']
    if script != 'MSAP.py':
        args += ['-o', str(output)]
    result = subprocess.run(args, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    report = tmp_path/'protein.sequence_changes.tsv' if script == 'MSAP.py' else output/'reports/sequence_changes.tsv'
    edits = rows(report)
    assert len(edits) == 4
    assert all(r['input_file'] == str(source) and r['original'] == '.' for r in edits)
    assert all(r['source_aa_positions'] != 'NA' for r in edits)
    alignment_dir = tmp_path if script == 'MSAP.py' else (output/'alignments' if script == 'phyloprep.py' else output)
    assert {r.id: str(r.seq) for r in SeqIO.parse(alignment_dir/'protein.mafft.prot.aln', 'fasta')} == {'a': 'MAXG', 'b': 'MAAG'}
    assert not (alignment_dir/'protein.protein.cleaned.fasta').exists()
    if script != 'MSAP.py':
        repeated = subprocess.run(args, cwd=tmp_path, env=env, capture_output=True, text=True)
        assert repeated.returncode == 0 and 'already complete' in repeated.stdout


def test_header_only_protein_report_survives_merge(tmp_path):
    source, output, report = [tmp_path/n for n in ('input.fa', 'clean.fa', 'changes.tsv')]
    source.write_text('>a\nMAGG\n>b\nMAAG\n')
    preprocess_protein(source, output, report)
    merged = tmp_path/'merged.tsv'
    merge_reports([report], merged)
    assert merged.read_text() == report.read_text()
    assert 'source_aa_positions' in merged.read_text()


def test_all_ambiguous_protein_symbols_are_masked_as_x(tmp_path):
    source, output, report = [tmp_path/n for n in ('input.fa', 'clean.fa', 'changes.tsv')]
    source.write_text('>a\nMBZJOU?X\n>b\nMAAAAAAX\n')
    preprocess_protein(source, output, report)
    assert {r.id: str(r.seq) for r in SeqIO.parse(output, 'fasta')}['a'] == 'MXXXXXXX'
    edits = rows(report)
    assert [r['original'] for r in edits] == list('BZJOU?')
    assert all(r['replacement'] == 'X' and r['action'] == 'replace' for r in edits)


def test_stop_option_invalidates_resume_and_trimming_uses_cleaned_protein(tmp_path, fake_aligner):
    _, env = fake_aligner
    source, output = tmp_path/'protein.fa', tmp_path/'out'
    source.write_text('>a\nMAAG\n>b\nMAAG\n')
    args = [sys.executable, str(ROOT/'MSAP_batch.py'), '-i', str(source),
            '-o', str(output), '-st', 'prot', '-t', '1']
    for extra in ([], ['--protein-stop-symbol', '.']):
        result = subprocess.run(args + extra, cwd=tmp_path, env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert 'already complete' not in result.stdout
    assert Path(env['MSAP_TEST_LOG']).read_text().count('run') == 2
    source.write_text('>a\nMA.G..\n>b\nMAAG.\n')
    result = subprocess.run(args + ['--protein-stop-symbol', '.'], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert {str(r.seq) for r in SeqIO.parse(output/'protein.mafft.prot.trimal.aln', 'fasta')} == {'MAG'}
    assert len(rows(output/'reports/sequence_changes.tsv')) == 4


@pytest.mark.parametrize('symbol', ['X', '-', 'AA', '1', '>'])
def test_invalid_stop_symbol_rejected(tmp_path, symbol):
    result = subprocess.run([sys.executable, str(ROOT/'MSAP.py'), '-i', 'unused.fa',
                             '-st', 'prot', '--protein-stop-symbol', symbol],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0 and '--protein-stop-symbol' in result.stderr
