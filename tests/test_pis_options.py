"""Verify representation-specific PIS filtering and reject ineffective options."""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('mode,flag,sequence', [
    ('nucl', '--min-nucl-pis', 'AC'),
    ('codon', '--min-nucl-pis', 'GCT'),
    ('prot', '--min-protein-pis', 'A'),
])
def test_standalone_pis_threshold(tmp_path, mode, flag, sequence):
    source = tmp_path / 'gene.fa'
    source.write_text(''.join(f'>{name}\n{sequence}\n' for name in 'abcd'))
    for threshold, passed in [(0, True), (1, False)]:
        prefix = tmp_path / f'qc-{threshold}'
        result = subprocess.run([sys.executable, str(REPO / 'alignment_qc.py'),
                                 '-i', str(source), '-st', mode, '-p', str(prefix),
                                 flag, str(threshold)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert (str(source) in prefix.with_suffix('.pass.tsv').read_text()) == passed


@pytest.mark.parametrize('script', ['phyloprep.py', 'alignment_qc.py'])
@pytest.mark.parametrize('mode,flag', [
    ('nucl', '--min-protein-pis'), ('prot', '--min-nucl-pis'),
    ('codon', '--min-parsimony-informative-sites'), ('codon', '--min-codon-pis'),
])
def test_invalid_pis_options_are_rejected(tmp_path, script, mode, flag):
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCT\n>b\nGCT\n')
    result = subprocess.run([sys.executable, str(REPO / script), '-i', str(source),
                             '-st', mode, flag, '0', *(['-p', str(tmp_path / 'qc')] if script == 'alignment_qc.py' else [])], capture_output=True, text=True)
    assert result.returncode != 0
    assert flag in result.stderr


@pytest.mark.parametrize('script', ['phyloprep.py', 'alignment_qc.py'])
@pytest.mark.parametrize('mode,flag', [('nucl', '--min-nucl-pis'), ('prot', '--min-protein-pis')])
def test_negative_pis_threshold_is_rejected(tmp_path, script, mode, flag):
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCT\n>b\nGCT\n')
    result = subprocess.run([sys.executable, str(REPO / script), '-i', str(source),
                             '-st', mode, flag, '-1', *(['-p', str(tmp_path / 'qc')] if script == 'alignment_qc.py' else [])], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'non-negative' in result.stderr


@pytest.mark.parametrize('mode,flag,sequence', [('nucl', '--min-nucl-pis', 'AC'), ('prot', '--min-protein-pis', 'A')])
def test_single_representation_pipeline_pis(tmp_path, fake_aligner, mode, flag, sequence):
    _, env = fake_aligner
    source = tmp_path / 'gene.fa'
    source.write_text(''.join(f'>{name}\n{sequence}\n' for name in 'abcd'))
    for threshold, passes in [(0, True), (1, False)]:
        output = tmp_path / f'results-{threshold}'
        result = subprocess.run([sys.executable, str(REPO / 'phyloprep.py'),
                                 '-i', str(source), '-st', mode, '--notrim', '-t', '1',
                                 '-o', str(output), flag, str(threshold)],
                                env=env, capture_output=True, text=True)
        assert (result.returncode == 0) == passes, result.stdout + result.stderr
        assert (output / f'matrices/{mode}/raw/supermatrix.fasta').exists() == passes
