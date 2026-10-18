"""Large unwrapped supermatrices must not be retained as Python record lists."""
import csv
import importlib.util
import random
from pathlib import Path
import tracemalloc

import pytest

from fasta2nex import fasta2nex
from fasta2phy import fasta2phy
from msap_matrix import AlignmentSpool
from split_codon_seqence_alignment import split

spec = importlib.util.spec_from_file_location('fourfold', Path(__file__).parents[1] / 'extract_4-fold_degenerated_sites.py')
fourfold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fourfold)


@pytest.mark.parametrize('operation', ['split', 'fourfold', 'phy', 'nex'])
def test_large_unwrapped_matrix_is_bounded(tmp_path, monkeypatch, operation):
    source = tmp_path / 'matrix.fa'
    with source.open('w') as handle:
        for taxon in range(10):
            handle.write(f'>taxon{taxon}\n')
            handle.write('GCT' * 100001 + '\n')
    monkeypatch.chdir(tmp_path)
    output = tmp_path / ('out.' + operation)
    tracemalloc.start()
    try:
        if operation == 'split':
            split(source)
        elif operation == 'fourfold':
            fourfold.extract_four_fold_site(source, output)
        else:
            {'phy': fasta2phy, 'nex': fasta2nex}[operation](source, output)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 2_500_000, (operation, peak)
    if operation == 'split':
        for position, base in zip(('1st', '2nd', '3rd'), 'GCT'):
            with AlignmentSpool(tmp_path / f'codon{position}.matrix.fasta') as matrix:
                assert matrix.length == 100001
                assert ''.join(matrix.chunks(0)) == base * 100001
    elif operation == 'fourfold':
        with AlignmentSpool(output) as matrix:
            assert matrix.length == 100001
            assert ''.join(matrix.chunks(9)) == 'T' * 100001
    else:
        assert output.stat().st_size > 3_000_000
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(
        ['matrix.fa'] + ([f'codon{pos}.matrix.fasta' for pos in ('1st', '2nd', '3rd')]
                         if operation == 'split' else [output.name]))


@pytest.mark.parametrize('converter', [fasta2phy, fasta2nex])
def test_export_audit_coordinates_cross_chunk_boundary(tmp_path, converter):
    source, output = tmp_path / 'in.fa', tmp_path / 'out'
    source.write_text('>a\n' + 'A' * 65534 + 'U?RAC\n')
    converter(source, output)
    with Path(str(output) + '.sequence_changes.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [row['source_nt_positions'] for row in rows] == ['65535-65535', '65537-65537']


@pytest.mark.parametrize('content', ['>a\nGCT\n>a\nGCT\n', '>a\nGCT\n>b\nGC\n', '>a\n\n'])
def test_invalid_matrix_keeps_existing_output(tmp_path, content):
    source, output = tmp_path / 'in.fa', tmp_path / 'out'
    source.write_text(content)
    output.write_text('previous')
    for operation in (fasta2phy, fasta2nex, fourfold.extract_four_fold_site):
        with pytest.raises(ValueError):
            operation(source, output)
        assert output.read_text() == 'previous'


@pytest.mark.parametrize('code', [1, 2, 5])
def test_fourfold_matches_column_definition_across_blocks(tmp_path, code):
    rng = random.Random(17)
    fragments = [[], [], []]
    expected = [[], [], []]
    eligible = fourfold.get_four_fold_codons(code)
    for _ in range(22001):
        prefix = rng.choice(['GC', 'CG', 'AG', 'TC', 'AT'])
        codons = [prefix + rng.choice('ACGT') for _ in fragments]
        if rng.random() < 0.2:
            codons[1] = rng.choice(['NNN', '---', 'GCT', 'ATG'])
        for row, codon in zip(fragments, codons):
            row.append(codon)
        if all(codon in eligible and codon[:2] == codons[0][:2] for codon in codons):
            for row, codon in zip(expected, codons):
                row.append(codon[2])
    sequences = [''.join(row) for row in fragments]
    source, output = tmp_path / 'in.fa', tmp_path / 'fourfold.fa'
    source.write_text(''.join(f'>t{i}\n{seq}\n' for i, seq in enumerate(sequences)))
    fourfold.extract_four_fold_site(source, output, code)
    with AlignmentSpool(output) as matrix:
        assert [''.join(matrix.chunks(i)) for i in range(3)] == [''.join(row) for row in expected]
