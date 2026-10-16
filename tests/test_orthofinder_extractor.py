import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'extract_orthofinder_orthogroups.py'


def run(tmp_path, *args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], cwd=tmp_path,
                          text=True, capture_output=True)


def write_inputs(tmp_path):
    table = tmp_path / 'Orthogroups.tsv'
    table.write_text(
        'Orthogroup\tHuman\tMouse\tFish\n'
        'OG0001\tHuman_gene:1\tMouse_gene:2\tFish_gene:3\n'
        'OG0002\tHuman_gene:4, Human_gene:5\tMouse_gene:4\tFish_gene:8\n'
        'OG0003\tHuman_gene:6\t\tFish_gene:7\n'
    )
    cds = tmp_path / 'cds'
    cds.mkdir()
    (cds / 'Human.fa').write_text('>gene:1 original description\nATG\n>gene:4\nATG\n>gene:5\nATG\n>gene:6\nATG\n')
    (cds / 'Mouse.fasta').write_text('>gene:2\nATC\n>gene:4\nATC\n')
    (cds / 'Fish.pep').write_text('>gene:3\nATT\n>gene:8\nATT\n>gene:7\nATT\n')
    return table, cds


def test_default_ids_match_orthofinder_without_x_and_map_is_written(tmp_path):
    table, cds = write_inputs(tmp_path)
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    records = list(SeqIO.parse(tmp_path / 'out/OG0001.fasta', 'fasta'))
    assert [record.id for record in records] == ['Human_gene_1', 'Mouse_gene_2', 'Fish_gene_3']
    assert (tmp_path / 'out/map.tsv').read_text().splitlines() == [
        'Human_gene_1\tHuman', 'Mouse_gene_2\tMouse', 'Fish_gene_3\tFish',
    ]
    assert (tmp_path / 'out/orthogroups.pathlist').read_text().splitlines() == [
        str((tmp_path / 'out/OG0001.fasta').resolve()),
    ]
    assert not (tmp_path / 'out/OG0002.fasta').exists()


def test_x_writes_unprefixed_ids_and_map_to_custom_path(tmp_path):
    table, cds = write_inputs(tmp_path)
    mapping = tmp_path / 'gene_to_taxon.tsv'
    pathlist = tmp_path / 'phyloprep_inputs.txt'
    result = run(tmp_path, '-og', table, '-i', cds, '-X', '--map-output', mapping,
                 '--pathlist-output', pathlist, '-o', tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    assert [record.id for record in SeqIO.parse(tmp_path / 'out/OG0001.fasta', 'fasta')] == ['gene_1', 'gene_2', 'gene_3']
    assert mapping.read_text().splitlines() == ['gene_1\tHuman', 'gene_2\tMouse', 'gene_3\tFish']
    assert pathlist.read_text().splitlines() == [str((tmp_path / 'out/OG0001.fasta').resolve())]
    assert not (tmp_path / 'out/map.tsv').exists()
    assert not (tmp_path / 'out/orthogroups.pathlist').exists()


def test_unprefixed_table_gets_taxon_prefix_without_extractor_x(tmp_path):
    table, cds = write_inputs(tmp_path)
    table.write_text(
        'Orthogroup\tHuman\tMouse\tFish\n'
        'OG0001\tgene:1\tgene:2\tgene:3\n'
    )
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    assert [record.id for record in SeqIO.parse(tmp_path / 'out/OG0001.fasta', 'fasta')] == [
        'Human_gene_1', 'Mouse_gene_2', 'Fish_gene_3',
    ]
    assert (tmp_path / 'out/map.tsv').read_text().splitlines() == [
        'Human_gene_1\tHuman', 'Mouse_gene_2\tMouse', 'Fish_gene_3\tFish',
    ]


def test_mixed_orthofinder_id_styles_require_an_explicit_nonstandard_table_fix(tmp_path):
    table, cds = write_inputs(tmp_path)
    table.write_text(
        'Orthogroup\tHuman\tMouse\tFish\n'
        'OG0001\tHuman_gene:1\tgene:2\tFish_gene:3\n'
    )
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out')
    assert result.returncode == 2
    assert 'mixes taxon-prefixed and unprefixed gene IDs' in result.stderr


def test_keep_drop_and_min_taxa_define_occupancy_after_taxon_filtering(tmp_path):
    table, cds = write_inputs(tmp_path)
    keep = run(tmp_path, '-og', table, '-i', cds, '--keep-taxa', 'Human', 'Fish', '-o', tmp_path / 'keep')
    assert keep.returncode == 0, keep.stderr
    assert (tmp_path / 'keep/OG0001.fasta').exists()
    assert (tmp_path / 'keep/OG0003.fasta').exists()
    drop = run(tmp_path, '-og', table, '-i', cds, '--drop-taxa', 'Human', '-o', tmp_path / 'drop')
    assert drop.returncode == 0, drop.stderr
    assert (tmp_path / 'drop/OG0002.fasta').exists()
    partial = run(tmp_path, '-og', table, '-i', cds, '--min-taxa', '2', '-o', tmp_path / 'partial')
    assert partial.returncode == 0, partial.stderr
    assert {path.stem for path in (tmp_path / 'partial').glob('OG*.fasta')} == {'OG0001', 'OG0002', 'OG0003'}
    assert [record.id for record in SeqIO.parse(tmp_path / 'partial/OG0002.fasta', 'fasta')] == ['Mouse_gene_4', 'Fish_gene_8']
    invalid = run(tmp_path, '-og', table, '-i', cds, '--min-taxa', '4')
    assert invalid.returncode == 2 and '--min-taxa must be between 2 and 3' in invalid.stderr


def test_sequence_directory_uses_orthofinder_filename_stem_rules(tmp_path):
    table, cds = write_inputs(tmp_path)
    (cds / 'Mouse.fasta').unlink()
    (cds / 'Mouse.fna').write_text('>gene:2\nATC\n')
    result = run(tmp_path, '-og', table, '-i', cds)
    assert result.returncode == 2
    assert 'Sequence FASTA is missing for selected taxa: Mouse' in result.stderr


def test_f_alias_accepts_protein_fasta_directory(tmp_path):
    table, _ = write_inputs(tmp_path)
    proteins = tmp_path / 'proteins'
    proteins.mkdir()
    (proteins / 'Human.faa').write_text('>gene:1\nMK\n')
    (proteins / 'Mouse.faa').write_text('>gene:2\nMA\n')
    (proteins / 'Fish.faa').write_text('>gene:3\nMV\n')
    result = run(tmp_path, '-og', table, '-f', proteins, '-o', tmp_path / 'protein_ogs')
    assert result.returncode == 0, result.stderr
    assert [str(record.seq) for record in SeqIO.parse(tmp_path / 'protein_ogs/OG0001.fasta', 'fasta')] == ['MK', 'MA', 'MV']


def test_keep_drop_overlap_errors_and_all_copies_preserves_paralogs(tmp_path):
    table, cds = write_inputs(tmp_path)
    overlap = run(tmp_path, '-og', table, '-i', cds,
                  '--keep-taxa', 'Human', 'Mouse', '--drop-taxa', 'Mouse')
    assert overlap.returncode == 2
    assert '--keep-taxa and --drop-taxa select the same taxa: Mouse' in overlap.stderr
    result = run(tmp_path, '-og', table, '-i', cds,
                 '--copy-mode', 'all-copies', '--min-taxa', '3', '-o', tmp_path / 'all_copies')
    assert result.returncode == 0, result.stderr
    records = list(SeqIO.parse(tmp_path / 'all_copies/OG0002.fasta', 'fasta'))
    assert [record.id for record in records] == [
        'Human_gene_4', 'Human_gene_5', 'Mouse_gene_4', 'Fish_gene_8',
    ]
    assert (tmp_path / 'all_copies/OG0003.fasta').exists() is False


@pytest.mark.parametrize('prefixed', [False, True])
@pytest.mark.parametrize('x', [False, True])
def test_export_prefix_policy_is_independent_of_table_style(tmp_path, prefixed, x):
    table, cds = write_inputs(tmp_path)
    values = ['Human_gene:1', 'Mouse_gene:2', 'Fish_gene:3'] if prefixed else ['gene:1', 'gene:2', 'gene:3']
    table.write_text('Orthogroup\tHuman\tMouse\tFish\nOG0001\t' + '\t'.join(values) + '\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output, *(['-X'] if x else []))
    assert result.returncode == 0, result.stderr
    records = list(SeqIO.parse(output / 'OG0001.fasta', 'fasta'))
    expected = ['gene_1', 'gene_2', 'gene_3'] if x else ['Human_gene_1', 'Mouse_gene_2', 'Fish_gene_3']
    assert [r.id for r in records] == expected
    assert [str(r.seq) for r in records] == ['ATG', 'ATC', 'ATT']
    assert (output / 'map.tsv').read_text().splitlines() == [f'{identifier}\t{taxon}' for identifier, taxon in zip(expected, ['Human', 'Mouse', 'Fish'])]


def test_default_prefixes_disambiguate_shared_ids_but_x_rejects(tmp_path):
    cds = tmp_path / 'cds'
    cds.mkdir()
    for taxon in ['Human', 'Mouse']:
        (cds / f'{taxon}.fa').write_text('>gene1\nATG\n')
    table = tmp_path / 'Orthogroups.tsv'
    table.write_text('Orthogroup\tHuman\tMouse\nOG1\tgene1\tgene1\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output)
    assert result.returncode == 0, result.stderr
    assert [r.id for r in SeqIO.parse(output / 'OG1.fasta', 'fasta')] == ['Human_gene1', 'Mouse_gene1']
    rejected = tmp_path / 'rejected'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', rejected, '-X')
    assert result.returncode != 0
    assert 'output IDs are not unique' in result.stderr
    assert not (rejected / 'OG1.fasta').exists()


def test_x_rejects_ambiguous_map_ids_across_different_orthogroups(tmp_path):
    cds = tmp_path / 'cds'
    cds.mkdir()
    (cds / 'Human.fa').write_text('>gene1\nATG\n')
    (cds / 'Mouse.fa').write_text('>mouse1\nATG\n>mouse2\nATG\n')
    (cds / 'Fish.fa').write_text('>gene1\nATG\n')
    table = tmp_path / 'Orthogroups.tsv'
    table.write_text('Orthogroup\tHuman\tMouse\tFish\nOG1\tgene1\tmouse1\t\nOG2\t\tmouse2\tgene1\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output, '--min-taxa', '2', '-X')
    assert result.returncode != 0
    assert 'maps to multiple taxa across orthogroups' in result.stderr
    assert not (output / 'map.tsv').exists()
    assert not list(output.glob('*.fasta'))


def test_longest_copy_retains_one_per_taxon_and_records_candidates(tmp_path):
    table, cds = write_inputs(tmp_path)
    (cds / 'Human.fa').write_text('>gene:1\nATG\n>gene:4\nATG\n>gene:5\nATGAAA\n>gene:6\nATG\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output,
                 '--copy-mode', 'all-copies', '--longest-copy')
    assert result.returncode == 0, result.stderr
    records = list(SeqIO.parse(output / 'OG0002.fasta', 'fasta'))
    assert [r.id for r in records] == ['Human_gene_5', 'Mouse_gene_4', 'Fish_gene_8']
    assert str(records[0].seq) == 'ATGAAA'
    mapping = (output / 'map.tsv').read_text()
    assert 'Human_gene_5\tHuman\n' in mapping
    assert 'Human_gene_4\tHuman\n' not in mapping
    assert len((output / 'orthogroups.pathlist').read_text().splitlines()) == 2
    selection = (output / 'copy_selection.tsv').read_text()
    assert selection.startswith('orthogroup\ttaxon\tgene_id\tlength\teffective_length\tsequence_type\tselected\n')
    assert 'OG0002\tHuman\tgene_4\t3\t3\tcodon\tno\n' in selection
    assert 'OG0002\tHuman\tgene_5\t6\t6\tcodon\tyes\n' in selection


def test_longest_copy_ties_are_independent_of_cell_order(tmp_path):
    table, cds = write_inputs(tmp_path)
    for index, genes in enumerate(['Human_gene:4, Human_gene:5', 'Human_gene:5, Human_gene:4']):
        contents = table.read_text()
        contents = contents.replace('Human_gene:4, Human_gene:5', genes).replace('Human_gene:5, Human_gene:4', genes)
        table.write_text(contents)
        out = tmp_path / f'out{index}'
        result = run(tmp_path, '-og', table, '-i', cds, '-o', out,
                     '--copy-mode', 'all-copies', '--longest-copy')
        assert result.returncode == 0, result.stderr
        assert next(SeqIO.parse(out / 'OG0002.fasta', 'fasta')).id == 'Human_gene_4'


def test_longest_copy_requires_all_copies(tmp_path):
    table, cds = write_inputs(tmp_path)
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out', '--longest-copy')
    assert result.returncode == 2
    assert '--longest-copy requires --copy-mode all-copies' in result.stderr
    assert not (tmp_path / 'out').exists()


def test_longest_copy_does_not_hide_missing_candidate_sequence(tmp_path):
    table, cds = write_inputs(tmp_path)
    (cds / 'Human.fa').write_text('>gene:1\nATG\n>gene:4\nATGAAAAAA\n>gene:6\nATG\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output,
                 '--copy-mode', 'all-copies', '--longest-copy')
    assert result.returncode == 2
    assert "sequence ID 'gene_5'" in result.stderr
    assert not list(output.glob('*.fasta'))
    assert not (output / 'copy_selection.tsv').exists()


def test_relaxed_single_copy_warns_about_downstream_missing_taxa(tmp_path):
    table, cds = write_inputs(tmp_path)
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out',
                 '--copy-mode', 'single-copy', '--min-taxa', 2)
    assert result.returncode == 0, result.stderr
    assert '--missing-taxa pad-gaps' in result.stderr
    assert '3 selected taxa' in result.stderr
    assert 'Multi-copy taxa are omitted' in result.stderr
    assert 'skip-gene' in result.stderr


def test_full_occupancy_does_not_warn_about_padding(tmp_path):
    table, cds = write_inputs(tmp_path)
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    assert '--missing-taxa pad-gaps' not in result.stderr


def test_longest_copy_warns_about_paralog_selection(tmp_path):
    table, cds = write_inputs(tmp_path)
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out',
                 '--copy-mode', 'all-copies', '--longest-copy')
    assert result.returncode == 0, result.stderr
    assert 'Different species may retain different paralogs' in result.stderr
    assert 'does not establish single-copy orthologs' in result.stderr
    assert 'copy_selection.tsv' in result.stderr


@pytest.mark.parametrize('sequence_type,first,second,first_score,second_score', [
    ('nucl', 'aCgTuNRYWSKMBDHVX?*-.', 'ACGTAC', 5, 6),
    ('prot', 'aCDEFGHIKLMNPQRSTVWYnbzjxuo?*-.', 'ACDEFGHIKLMNPQRSTVWYNN', 21, 22),
    ('nucl', 'NNNNNNNNNN', 'ACG', 0, 3),
    ('codon', 'NNNNNNNNNN', 'ACG', 0, 3),
    ('pseudogene', 'NNNNNNNNNN', 'ACG', 0, 3),
    ('prot', 'XXXXXXXXXX', 'NNN', 0, 3),
])
def test_longest_copy_counts_only_non_ambiguous_characters(tmp_path, sequence_type, first, second, first_score, second_score):
    import csv
    table, cds = write_inputs(tmp_path)
    (cds / 'Human.fa').write_text(f'>gene:1\nATG\n>gene:4\n{first}\n>gene:5\n{second}\n>gene:6\nATG\n')
    output = tmp_path / 'out'
    result = run(tmp_path, '-og', table, '-i', cds, '-o', output,
                 '--copy-mode', 'all-copies', '--longest-copy', '-st', sequence_type)
    assert result.returncode == 0, result.stderr
    records = list(SeqIO.parse(output / 'OG0002.fasta', 'fasta'))
    assert records[0].id == 'Human_gene_5'
    assert str(records[0].seq) == second
    with (output / 'copy_selection.tsv').open() as handle:
        rows = {r['gene_id']: r for r in csv.DictReader(handle, delimiter='\t') if r['orthogroup'] == 'OG0002' and r['taxon'] == 'Human'}
    assert rows['gene_4']['length'] == str(len(first))
    assert rows['gene_4']['effective_length'] == str(first_score)
    assert rows['gene_5']['effective_length'] == str(second_score)
    assert rows['gene_5']['sequence_type'] == sequence_type


@pytest.mark.parametrize('first,second', [('ACG', 'ACGNNNN'), ('NN', 'NNNNNN')])
def test_longest_copy_effective_length_ties_ignore_raw_length(tmp_path, first, second):
    table, cds = write_inputs(tmp_path)
    (cds / 'Human.fa').write_text(f'>gene:1\nATG\n>gene:4\n{first}\n>gene:5\n{second}\n>gene:6\nATG\n')
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out',
                 '--copy-mode', 'all-copies', '--longest-copy')
    assert result.returncode == 0, result.stderr
    record = next(SeqIO.parse(tmp_path / 'out' / 'OG0002.fasta', 'fasta'))
    assert record.id == 'Human_gene_4'
    assert str(record.seq) == first


@pytest.mark.parametrize('flag', ['-st', '--seqtype', '--seq-type'])
def test_longest_copy_sequence_type_aliases(tmp_path, flag):
    table, cds = write_inputs(tmp_path)
    (cds / 'Human.fa').write_text('>gene:1\nATG\n>gene:4\nXXXXX\n>gene:5\nNNN\n>gene:6\nATG\n')
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out',
                 '--copy-mode', 'all-copies', '--longest-copy', flag, 'prot')
    assert result.returncode == 0, result.stderr
    assert next(SeqIO.parse(tmp_path / 'out' / 'OG0002.fasta', 'fasta')).id == 'Human_gene_5'


@pytest.mark.parametrize('copy_mode,longest_copy,expected_human', [
    ('single-copy', False, {'gene_1'}),
    ('all-copies', False, {'gene_1', 'gene_4', 'gene_5'}),
    ('all-copies', True, {'gene_1', 'gene_4', 'gene_5'}),
])
def test_only_retained_og_candidates_are_loaded(tmp_path, monkeypatch, copy_mode, longest_copy, expected_human):
    import extract_orthofinder_orthogroups as extractor
    table, cds = write_inputs(tmp_path)
    with (cds / 'Human.fa').open('a') as handle:
        handle.write('>unused\n' + 'A' * 100000 + '\n')
    calls = {}
    original = extractor.load_sequences

    def load(path, wanted_ids=None):
        result = original(path, wanted_ids)
        calls[path.stem] = (set(wanted_ids), set(result))
        return result

    monkeypatch.setattr(extractor, 'load_sequences', load)
    extractor.extract_orthologs(table, cds, tmp_path / 'out', copy_mode=copy_mode, longest_copy=longest_copy)
    assert calls['Human'] == (expected_human, expected_human)
    assert len(calls) == 3
    assert (tmp_path / 'out/map.tsv').is_file()


def test_filtered_sequence_loading_keeps_existing_duplicate_validation(tmp_path):
    from extract_orthofinder_orthogroups import load_sequences
    source = tmp_path / 'source.fa'
    source.write_text('>wanted\nATG\n>unused:1\nAAA\n>unused_1\nCCC\n')
    with pytest.raises(ValueError, match='duplicate OrthoFinder-normalized sequence ID'):
        load_sequences(source, {'wanted'})


def test_filtered_sequence_loading_distinguishes_no_matches_from_empty_file(tmp_path):
    from extract_orthofinder_orthogroups import load_sequences
    source = tmp_path / 'source.fa'
    source.write_text('>other\nATG\n')
    assert load_sequences(source, {'missing'}) == {}
    source.write_text('')
    with pytest.raises(ValueError, match='No sequences found'):
        load_sequences(source, {'missing'})
