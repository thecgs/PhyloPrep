import subprocess
import sys
from pathlib import Path

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


def test_orthofinder_x_table_is_detected_and_preserved_without_extractor_x(tmp_path):
    table, cds = write_inputs(tmp_path)
    table.write_text(
        'Orthogroup\tHuman\tMouse\tFish\n'
        'OG0001\tgene:1\tgene:2\tgene:3\n'
    )
    result = run(tmp_path, '-og', table, '-i', cds, '-o', tmp_path / 'out')
    assert result.returncode == 0, result.stderr
    assert [record.id for record in SeqIO.parse(tmp_path / 'out/OG0001.fasta', 'fasta')] == [
        'gene_1', 'gene_2', 'gene_3',
    ]
    assert (tmp_path / 'out/map.tsv').read_text().splitlines() == [
        'gene_1\tHuman', 'gene_2\tMouse', 'gene_3\tFish',
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
