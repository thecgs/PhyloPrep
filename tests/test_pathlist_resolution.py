from pathlib import Path

import pytest

from msap_io import expand_input_paths


def test_relocated_list_uses_same_directory_basename(tmp_path):
    fasta = tmp_path / 'gene.fa'
    fasta.write_text('>a\nAC\n')
    paths = tmp_path / 'orthogroups.pathlist'
    paths.write_text('# comment\n/old/machine/data/gene.fa\n; comment\n')
    assert expand_input_paths([paths]) == [str(fasta.resolve())]


def test_relative_path_is_resolved_from_list_directory(tmp_path, monkeypatch):
    folder = tmp_path / 'data'
    folder.mkdir()
    fasta = folder / 'gene.fa'
    fasta.write_text('>a\nAC\n')
    paths = tmp_path / 'orthogroups.pathlist'
    paths.write_text('data/gene.fa\n')
    monkeypatch.chdir(folder)
    assert expand_input_paths([paths]) == [str(fasta.resolve())]


def test_existing_absolute_reference_wins_over_sibling(tmp_path):
    original = tmp_path / 'original'
    original.mkdir()
    fasta = original / 'gene.fa'
    fasta.write_text('>a\nAC\n')
    (tmp_path / 'gene.fa').write_text('>a\nTT\n')
    paths = tmp_path / 'orthogroups.pathlist'
    paths.write_text(str(fasta) + '\n')
    assert expand_input_paths([paths]) == [str(fasta.resolve())]


def test_missing_reference_names_list_and_entry(tmp_path):
    paths = tmp_path / 'orthogroups.pathlist'
    paths.write_text('/missing/gene.fa\n')
    with pytest.raises(ValueError, match='Missing FASTA referenced by.*orthogroups.pathlist.*gene.fa'):
        expand_input_paths([paths])


def test_windows_reference_can_fall_back_to_sibling(tmp_path):
    fasta = tmp_path / 'gene.fa'
    fasta.write_text('>a\nAC\n')
    paths = tmp_path / 'orthogroups.pathlist'
    paths.write_text('C:\\old\\gene.fa\n')
    assert expand_input_paths([paths]) == [str(fasta.resolve())]
