import importlib.util
import json
from pathlib import Path

import phyloprep

ROOT = Path(__file__).resolve().parents[1]


def browser_module():
    spec = importlib.util.spec_from_file_location('output_manifest_browser', ROOT / 'site/browser_runner.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_manifest_does_not_read_output_contents(tmp_path, monkeypatch):
    browser = browser_module()
    (tmp_path / 'matrix.fasta').write_text('>a\nAC\n')
    hidden = tmp_path / '.staging'
    hidden.mkdir()
    (hidden / 'private.fasta').write_text('private')
    monkeypatch.setattr(Path, 'read_text', lambda *args, **kwargs: (_ for _ in ()).throw(MemoryError('contents must not be read')))
    files = browser._visible_output_files(tmp_path, manifest_only=True)
    assert files == [dict(name='matrix.fasta', path=str(tmp_path / 'matrix.fasta'), size=6)]


def test_finalize_serializes_metadata_not_large_matrix(tmp_path, monkeypatch):
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    def variant(args, results, stage, variant):
        target = stage / 'prot' / variant
        target.mkdir(parents=True)
        (target / 'supermatrix.fasta').write_bytes(b'>a\n' + b'A' * (8 * 1024 * 1024) + b'\n')
    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    raw = browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=True), str(tmp_path / 'out'), str(ROOT), manifest_only=True)
    assert len(raw) < 10000
    result = json.loads(raw)
    assert result['ok']
    matrix = next(file for file in result['files'] if file['name'].endswith('supermatrix.fasta'))
    assert matrix['size'] > 8 * 1024 * 1024
    assert 'content' not in matrix
    assert Path(matrix['path']).exists()


def test_bounded_logs_retain_diagnostics_without_full_history():
    browser = browser_module()
    log = browser.BoundedLog(limit=100)
    for _ in range(10000):
        log.write('old message\n')
    log.write('last diagnostic')
    assert log.length <= 100
    assert log.getvalue().startswith('[Earlier log output truncated]\n')
    assert log.getvalue().endswith('last diagnostic')
    log = browser.BoundedLog()
    for _ in range(10000):
        log.write('x')
    assert len(log.parts) <= 1024
    log.write('A' * 1000000)
    assert log.length == 65536


def test_finalize_memory_failure_uses_prebuilt_response(tmp_path, monkeypatch):
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    def fail(*args, **kwargs):
        raise MemoryError('simulated memory exhaustion')
    monkeypatch.setattr(phyloprep, 'prepare_variant', fail)
    monkeypatch.setattr(browser.BoundedLog, 'getvalue', fail)
    monkeypatch.setattr(browser.traceback, 'print_exc', fail)
    monkeypatch.setattr(browser.json, 'dumps', fail)
    raw = browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=True), str(tmp_path / 'out'), str(ROOT), manifest_only=True)
    assert raw is browser._MEMORY_ERROR_RESULT
    result = json.loads(raw)
    assert not result['ok']
    assert 'Python/Pyodide' in result['error']


def test_incremental_finalize_exports_qc_without_publishing_staged_matrices(tmp_path, monkeypatch):
    browser = browser_module()
    import phyloprep
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    seen = []

    def variant(args, root, stage, name):
        assert source.exists()
        if name == 'trimmed':
            assert (root / 'qc/raw/alignment.details.tsv').exists()
        report = root / 'qc' / name / 'alignment.details.tsv'
        report.parent.mkdir(parents=True)
        report.write_text('locus\tpass\ngene\tyes\n')
        matrix = stage / 'prot' / name / 'supermatrix.fasta'
        matrix.parent.mkdir(parents=True)
        matrix.write_text('>a\nMA\n>b\nMA\n')
        seen.append(name)

    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT),
        manifest_only=True, incremental=True)
    try:
        for expected in ('raw', 'trimmed'):
            checkpoint = json.loads(browser.advance_phyloprep_finalize(steps))
            assert checkpoint['pending']
            assert seen[-1] == expected
            assert all(file['name'].startswith(('qc/', 'reports/')) for file in checkpoint['files'])
            assert all(Path(file['path']).exists() for file in checkpoint['files'])
            assert not (tmp_path / 'out/matrices').exists()
        completed = json.loads(browser.advance_phyloprep_finalize(steps))
        assert completed['ok']
        assert any(file['name'].startswith('matrices/') for file in completed['files'])
        assert source.exists()
    finally:
        steps.close()


def test_closing_incremental_finalize_restores_engine_and_streams(tmp_path, monkeypatch):
    import sys
    import phyloprep
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    original_run = phyloprep.run
    original_stdout, original_stderr = sys.stdout, sys.stderr
    monkeypatch.setattr(phyloprep, 'prepare_variant', lambda *args: None)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT),
        manifest_only=True, incremental=True)
    try:
        next(steps)
    finally:
        steps.close()
    assert phyloprep.run is original_run
    assert sys.stdout is original_stdout
    assert sys.stderr is original_stderr
    assert not list((tmp_path / 'out').glob('.matrices-*'))
    assert not (tmp_path / 'out/matrices').exists()
