"""Fingerprint work follows necessary I/O and stays local to the current run."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

import phyloprep
import rename_taxa
from msap_io import file_digest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('newline', ['\n', '\r\n'])
def test_rename_hashes_exact_bytes_during_existing_io(tmp_path, monkeypatch, newline):
    source, target = tmp_path / 'input.fa', tmp_path / 'renamed.fa'
    raw = newline.join(['>gène description', 'ATGGCT', '>b', 'ATGGCC', '']).encode()
    source.write_bytes(raw)
    def forbidden(*args, **kwargs):
        raise AssertionError('Renaming must not reread source/output to hash them')
    monkeypatch.setattr(rename_taxa, 'file_signature', forbidden)
    fingerprints = {}
    rename_taxa.rename_fasta(source, target, mapping={'gène': 'Espèce'}, fingerprints=fingerprints)
    assert fingerprints == {'source': hashlib.sha256(raw).hexdigest(),
                            'target': hashlib.sha256(target.read_bytes()).hexdigest()}
    assert target.read_text().startswith('>Espèce\nATGGCT\n')


def test_mapping_snapshot_hashes_bytes_while_parsing(tmp_path, monkeypatch):
    source = tmp_path / 'map.tsv'
    raw = b'# comment\r\na\tSpecies_A\r\nb\tSpecies_B\r\n'
    source.write_bytes(raw)
    monkeypatch.setattr(rename_taxa, 'file_signature', lambda *_: pytest.fail('Extra mapping read'))
    mapping, signature = rename_taxa.read_mapping_snapshot(source)
    assert mapping == {'a': 'Species_A', 'b': 'Species_B'}
    assert signature['sha256'] == hashlib.sha256(raw).hexdigest()
    assert signature['size'] == len(raw)


def test_fingerprint_cache_checks_metadata_and_is_not_persistent(tmp_path, monkeypatch):
    path = tmp_path / 'file.fa'
    path.write_bytes(b'AAA')
    original = rename_taxa.file_signature
    calls = []
    def signature(path):
        calls.append(str(path))
        return original(path)
    monkeypatch.setattr(rename_taxa, 'file_signature', signature)
    cache = rename_taxa.FileFingerprintCache()
    first = cache.signature(path)
    assert cache.signature(path) == first
    assert len(calls) == 1
    # Restore mtime to ensure that ctime also participates in cache invalidation.
    before = path.stat()
    path.write_bytes(b'CCC')
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert cache.signature(path)['sha256'] == hashlib.sha256(b'CCC').hexdigest()
    assert len(calls) == 2
    rename_taxa.FileFingerprintCache().signature(path)
    assert len(calls) == 3, 'A fresh run must verify file contents again'


def test_cli_rename_checkpoints_keep_real_source_and_output_hashes(tmp_path):
    source, mapping = tmp_path / 'gene.fa', tmp_path / 'mapping.tsv'
    source.write_text('>a\nATGGCT\n>b\nATGGCC\n')
    mapping.write_text('a\tSpecies_A\nb\tSpecies_B\n')
    parsed, signature = rename_taxa.read_mapping_snapshot(mapping)
    args = SimpleNamespace(mapping=mapping, mapping_signature=signature,
                           taxon_mapping=parsed, no_resume=False)
    root = tmp_path / 'out'
    target, = phyloprep.rename_alignment_matrices(args, [source], root, 'codon', 'raw')
    state = json.loads(phyloprep.rename_state_path(target.parent, source).read_text())
    assert state['source']['sha256'] == file_digest(source)
    assert state['outputs']['alignment']['sha256'] == file_digest(target)
    assert not args._fingerprints.entries.get(str(source.resolve()))
    assert not args._fingerprints.entries.get(str(target.resolve()))
    # Across-run resume must reject a modified target, even at unchanged size.
    target.write_text(target.read_text().replace('ATGGCT', 'ATGCCC'))
    again = SimpleNamespace(mapping=mapping, mapping_signature=signature,
                            taxon_mapping=parsed, no_resume=False)
    phyloprep.rename_alignment_matrices(again, [source], root, 'codon', 'raw')
    assert 'ATGCCC' not in target.read_text()


def test_browser_rename_reuses_verified_source_and_export_digest(tmp_path, monkeypatch):
    source, mapping = tmp_path / 'gene.fa', tmp_path / 'map.tsv'
    source.write_text('>a\nATGGCT\n>b\nATGGCC\n')
    mapping.write_text('a\tA\nb\tB\n')
    parsed, signature = rename_taxa.read_mapping_snapshot(mapping)
    cache = rename_taxa.FileFingerprintCache()
    cache.remember(source, file_digest(source))
    args = SimpleNamespace(mapping=mapping, mapping_signature=signature,
                           taxon_mapping=parsed, no_resume=True,
                           _browser_runtime=True, _fingerprints=cache)
    root = tmp_path / 'out'
    target, = phyloprep.rename_alignment_matrices(args, [source], root, 'codon', 'raw')
    assert not list(root.rglob('.phyloprep-rename-*.json'))
    spec = importlib.util.spec_from_file_location('fingerprint_browser', ROOT / 'site/browser_runner.py')
    browser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(browser)
    monkeypatch.setattr(rename_taxa, 'file_signature', lambda *_: pytest.fail('Export rehashed renamed data'))
    for _ in range(2):
        files = browser._export_output_files(root, [dict(name='renamed.fa', path=str(target), size=target.stat().st_size)], fingerprints=cache)
        assert files[0]['sha256'] == hashlib.sha256(target.read_bytes()).hexdigest()
