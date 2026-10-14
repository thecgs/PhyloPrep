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
