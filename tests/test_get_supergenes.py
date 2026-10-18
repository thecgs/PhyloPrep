import subprocess
import sys
import tracemalloc
from pathlib import Path
import pytest


SCRIPT = Path(__file__).parents[1] / "get_supergenes.py"


def write_alignment(path, records):
    path.write_text("".join(f">{taxon}\n{sequence}\n" for taxon, sequence in records))


def run_supergenes(tmp_path, missing_taxa=None):
    complete_gene = tmp_path / "complete.fasta"
    incomplete_gene = tmp_path / "incomplete.fasta"
    prefix = tmp_path / (missing_taxa or 'default')
    write_alignment(
        complete_gene,
        [("sp1", "AAA"), ("sp2", "AAA"), ("sp3", "AAA")],
    )
    write_alignment(incomplete_gene, [("sp1", "CCC"), ("sp2", "CCC")])

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(complete_gene),
            str(incomplete_gene),
            "-p",
            str(prefix),
            *(["--missing-taxa", missing_taxa] if missing_taxa else []),
        ],
        check=True,
    )
    return prefix


def test_skip_gene_omits_a_gene_missing_any_taxon(tmp_path):
    prefix = run_supergenes(tmp_path, "skip-gene")

    assert (prefix.parent / f"{prefix.name}.fasta").read_text() == (
        ">sp1\nAAA\n>sp2\nAAA\n>sp3\nAAA\n"
    )
    report = (prefix.parent / f"{prefix.name}.report.tsv").read_text()
    assert "incomplete.fasta\tskipped\tmissing taxa: sp3" in report
    assert (prefix.parent / f"{prefix.name}.phy").read_text() == (
        "3 3\nsp1 AAA\nsp2 AAA\nsp3 AAA\n"
    )
    partition_config = (prefix.parent / f"{prefix.name}_partition_finder.cfg.txt").read_text()
    assert f"alignment = {prefix.name}.phy;" in partition_config


@pytest.mark.parametrize('policy', [None, 'pad-gaps'])
def test_pad_gaps_preserves_a_gene_missing_a_taxon(tmp_path, policy):
    prefix = run_supergenes(tmp_path, policy)

    assert (prefix.parent / f"{prefix.name}.fasta").read_text() == (
        ">sp1\nAAACCC\n>sp2\nAAACCC\n>sp3\nAAA???\n"
    )
    report = (prefix.parent / f"{prefix.name}.report.tsv").read_text()
    assert "incomplete.fasta\tpadded\tmissing taxa: sp3" in report
    assert (prefix.parent / f"{prefix.name}.phy").read_text() == (
        "3 6\nsp1 AAACCC\nsp2 AAACCC\nsp3 AAA???\n"
    )


def test_skip_gene_writes_a_report_when_no_gene_can_be_included(tmp_path):
    first_gene = tmp_path / "first.fasta"
    second_gene = tmp_path / "second.fasta"
    prefix = tmp_path / "all_skipped"
    write_alignment(first_gene, [("sp1", "AAA")])
    write_alignment(second_gene, [("sp2", "CCC")])

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(first_gene),
            str(second_gene),
            "-p",
            str(prefix),
            "--missing-taxa", "skip-gene",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    report = (prefix.parent / f"{prefix.name}.report.tsv").read_text()
    assert "first.fasta\tskipped\tmissing taxa: sp2" in report
    assert "second.fasta\tskipped\tmissing taxa: sp1" in report


def test_rejects_colliding_partition_names(tmp_path):
    first_gene = tmp_path / "gene-1.fasta"
    second_gene = tmp_path / "gene.1.fasta"
    prefix = tmp_path / "collision"
    write_alignment(first_gene, [("sp1", "AAA"), ("sp2", "AAA")])
    write_alignment(second_gene, [("sp1", "CCC"), ("sp2", "CCC")])

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(first_gene),
            str(second_gene),
            "-p",
            str(prefix),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Duplicate partition names after filename normalization" in result.stderr
    assert "gene_1_fasta" in result.stderr


def test_trimmed_alignment_suffix_preserves_locus_names():
    from get_supergenes import tidy_name
    for aligner in ('mafft', 'muscle', 'prank', 'clustalw2', 'macse'):
        for seqtype in ('nucl', 'prot', 'codon'):
            for variant in ('aln.fasta', 'trimmed.aln.fasta'):
                assert tidy_name(f'/work/ND5.{aligner}.{seqtype}.{variant}') == 'ND5'


def test_pipeline_can_defer_phylip_output(tmp_path):
    source = tmp_path / 'gene.mafft.nucl.aln.fasta'
    write_alignment(source, [('a', 'ACGT'), ('b', 'ACGT')])
    prefix = tmp_path / 'supermatrix'
    result = subprocess.run([sys.executable, str(SCRIPT), '-i', str(source),
                             '-p', str(prefix), '--no-phylip'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path/'supermatrix.fasta').exists()
    assert not (tmp_path/'supermatrix.phy').exists()
    assert 'alignment = supermatrix.phy;' in (tmp_path/'supermatrix_partition_finder.cfg.txt').read_text()


def test_concatenation_peak_does_not_scale_with_total_matrix(tmp_path):
    from get_supergenes import get_supergenes, read_alignment
    inputs = []
    for index in range(30):
        source = tmp_path / f'gene{index}.fa'
        write_alignment(source, [('a', 'A' * 100000), ('b', 'C' * 100000)])
        inputs.append(str(source))
    read_alignment(inputs[0])  # Warm the parser before measuring this operation.
    tracemalloc.start()
    try:
        get_supergenes(inputs, str(tmp_path / 'matrix'))
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    # The matrix contains six million characters; the old implementation kept
    # all locus sequences and joined matrices simultaneously (>12 MB).
    assert peak < 3 * 1024 * 1024, f'Unexpected Python allocation peak: {peak}'
    from Bio import SeqIO
    records = list(SeqIO.parse(tmp_path / 'matrix.fasta', 'fasta'))
    assert [(r.id, len(r.seq)) for r in records] == [('a', 3000000), ('b', 3000000)]
    assert not list(tmp_path.glob('.supergenes-*'))


def test_spool_handles_many_taxa_and_padding(tmp_path):
    from get_supergenes import get_supergenes
    records = [(f'sp{i}', 'AAA') for i in range(40)]
    a, b = tmp_path / 'a.fa', tmp_path / 'b.fa'
    write_alignment(a, records)
    write_alignment(b, [('sp0', 'CCC')])
    prefix = str(tmp_path / 'matrix')
    get_supergenes([str(a), str(b)], prefix, missing_taxa='pad-gaps')
    fasta = Path(prefix + '.fasta').read_text()
    assert '>sp0\nAAACCC\n' in fasta
    assert '>sp39\nAAA???\n' in fasta
    assert Path(prefix + '.phy').read_text().startswith('40 6\nsp0 AAACCC\n')
    assert 'b_fa = 4-6;' in Path(prefix + '_part_iqtree.txt').read_text()


@pytest.mark.parametrize('content,reason', [
    ('>a\nAAA\n>b\nCC\n', 'has length 2, expected 3'),
    ('>a\nAAA\n>a\nCCC\n', 'duplicate sequence ID'),
    ('', 'empty alignment'),
    ('>a\nAA\n>b\nCC\n', 'alignment length differs from the first pass'),
])
def test_invalid_second_pass_preserves_old_matrix_and_cleans_spools(tmp_path, monkeypatch, content, reason):
    import get_supergenes as module
    source = tmp_path / 'a.fa'
    write_alignment(source, [('a', 'AAA'), ('b', 'CCC')])
    prefix = str(tmp_path / 'matrix')
    Path(prefix + '.fasta').write_text('old matrix')
    original = module.read_alignment
    def changing(path):
        source.write_text(content)
        return original(path)
    monkeypatch.setattr(module, 'read_alignment', changing)
    with pytest.raises(ValueError, match=reason):
        module.get_supergenes([str(source)], prefix)
    assert Path(prefix + '.fasta').read_text() == 'old matrix'
    assert not list(tmp_path.glob('.supergenes-*'))


@pytest.mark.parametrize('policy', ['pad-gaps', 'skip-gene'])
def test_two_sequence_passes_without_content_hashing(tmp_path, monkeypatch, policy):
    import get_supergenes as module
    import msap_io
    from collections import Counter
    sources = [tmp_path / 'first.fa', tmp_path / 'second.fa']
    for source, sequence in zip(sources, ['AAA', 'CCC']):
        write_alignment(source, [('a', sequence), ('b', sequence)])
    reads = Counter()
    original = module.SeqIO.parse
    def parse(path, *args, **kwargs):
        reads[Path(path)] += 1
        return original(path, *args, **kwargs)
    def unexpected_digest(*args, **kwargs):
        raise AssertionError('Concatenation must not hash alignment contents')
    monkeypatch.setattr(module.SeqIO, 'parse', parse)
    monkeypatch.setattr(msap_io, 'file_digest', unexpected_digest)
    prefix = str(tmp_path / 'matrix')
    module.get_supergenes([str(source) for source in sources], prefix, missing_taxa=policy)
    assert reads == Counter({source: 2 for source in sources})
    assert Path(prefix + '.fasta').read_text() == '>a\nAAACCC\n>b\nAAACCC\n'


@pytest.mark.parametrize('content,reason', [
    ('>a\nAAA\n>a\nCCC\n', 'duplicate sequence ID'),
    ('>a\nAAA\n>b\nCC\n', 'has length 2, expected 3'),
])
def test_scan_preserves_structural_checks(tmp_path, content, reason):
    from get_supergenes import get_supergenes
    source = tmp_path / 'invalid.fa'
    source.write_text(content)
    with pytest.raises(ValueError, match=reason):
        get_supergenes([str(source)], str(tmp_path / 'matrix'))
    assert not (tmp_path / 'matrix.fasta').exists()
