"""Terminal-stop permission belongs to standalone QC, not preprocessing workflows."""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('mode,sequence', [('codon', 'GCTTAA'), ('prot', 'A*')])
def test_standalone_terminal_stop_permission_preserves_input(tmp_path, mode, sequence):
    source = tmp_path / 'gene.fa'
    original = f'>a\n{sequence}\n>b\n{sequence}\n'
    source.write_text(original)
    for allow in [False, True]:
        prefix = tmp_path / f'qc-{allow}'
        result = subprocess.run([sys.executable, str(REPO / 'alignment_qc.py'),
                                 '-i', str(source), '-p', str(prefix), '-st', mode,
                                 *(['--allow-terminal-stop'] if allow else [])],
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert (str(source) in prefix.with_suffix('.pass.tsv').read_text()) == allow
        assert source.read_text() == original


@pytest.mark.parametrize('mode', ['codon', 'prot', 'pseudogene'])
def test_complete_workflow_rejects_terminal_stop_switch(tmp_path, mode):
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCTTAA\n>b\nGCTTAA\n')
    result = subprocess.run([sys.executable, str(REPO / 'phyloprep.py'),
                             '-i', str(source), '-st', mode, '--allow-terminal-stop'],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert 'unrecognized arguments: --allow-terminal-stop' in result.stderr
