"""Exercise the full CLI chain with a deterministic external aligner."""
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import AlignIO

REPO = Path(__file__).resolve().parents[1]


def invoke(tmp_path, env, *args):
    return subprocess.run(
        [sys.executable, str(REPO / "phyloprep.py"), *map(str, args)],
        cwd=tmp_path, env=env, capture_output=True, text=True,
    )


def check_formats(directory, stem, expected):
    for extension, fmt in (("fasta", "fasta"), ("phy", "phylip-relaxed"), ("nex", "nexus")):
        records = AlignIO.read(directory / f"{stem}.{extension}", fmt)
        assert {r.id: str(r.seq) for r in records} == expected


@pytest.mark.parametrize("mapping_alias", [None, "-m", "--mapping", "--map"])
def test_codon_pipeline_and_resume(tmp_path, fake_aligner, mapping_alias):
    _, env = fake_aligner
    source = tmp_path / "gene.fa"
    source.write_text(">a\nGCTGGT\n>b\nGCCGGA\n")
    mapping = tmp_path / "taxa.tsv"
    mapping.write_text("a\tSpecies_A\nb\tSpecies_B\n")
    inputs = tmp_path / "genes.list"
    inputs.write_text(str(source) + "\n")
    output = tmp_path / "results"
    args = ["-i", inputs, "-o", output, "-t", "1", "-g", "1"]
    if mapping_alias:
        args.extend([mapping_alias, mapping])
    result = invoke(tmp_path, env, *args)
    assert result.returncode == 0, result.stdout + result.stderr
    a, b = ("Species_A", "Species_B") if mapping_alias else ("a", "b")
    if mapping_alias:
        raw_alignment = output / "alignments/gene.mafft.codon.aln"
        assert [record.id for record in AlignIO.read(raw_alignment, "fasta")] == ["a", "b"]
        renamed_alignment = output / "renamed/codon/raw/gene.mafft.codon.aln"
        assert [record.id for record in AlignIO.read(renamed_alignment, "fasta")] == ["Species_A", "Species_B"]
    for variant in ("raw", "trimmed"):
        directory = output / "matrices" / "codon" / variant
        check_formats(directory, "supermatrix", {a: "GCTGGT", b: "GCCGGA"})
        check_formats(directory, "codon1st.supermatrix", {a: "GG", b: "GG"})
        check_formats(directory, "codon2nd.supermatrix", {a: "CG", b: "CG"})
        check_formats(directory, "codon3rd.supermatrix", {a: "TT", b: "CA"})
        check_formats(directory, "fourfold", {a: "TT", b: "CA"})
        protein_dir = output / "matrices" / "prot" / variant
        check_formats(protein_dir, "supermatrix", {a: "AG", b: "AG"})
        assert "datatype=protein" in (protein_dir / "supermatrix.nex").read_text().lower()
        assert not list(protein_dir.glob("codon*"))
        suffix = ".aln" if variant == "raw" else ".trimal.aln"
        for seqtype in ("codon", "prot"):
            paths = (output / "qc" / seqtype / variant / "alignment.pass.tsv").read_text()
            assert f".{seqtype}{suffix}" in paths
    assert source.read_text() == ">a\nGCTGGT\n>b\nGCCGGA\n"
    result = invoke(tmp_path, env, *args)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "skipped (already complete)" in result.stdout
    if mapping_alias:
        assert "rename gene.mafft.codon.aln: skipped (already complete)" in result.stdout
    assert Path(env["MSAP_TEST_LOG"]).read_text().splitlines() == ["run"]


def test_unmapped_gene_ids_that_cannot_concatenate_keep_a_report_and_stop(tmp_path, fake_aligner):
    _, env = fake_aligner
    first, second = tmp_path / "gene1.fa", tmp_path / "gene2.fa"
    first.write_text(">gene_a\nACGT\n>gene_b\nACGA\n")
    second.write_text(">gene_c\nACGT\n>gene_d\nACGA\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", first, second, "-st", "nucl", "--notrim", "-o", output, "-t", "1")
    assert result.returncode == 2
    assert "provide phyloprep.py -m TAXA.tsv" in result.stderr
    report = output / "reports/concatenation/nucl/raw/original_ids.report.tsv"
    assert report.exists()
    assert report.read_text().count("\tskipped\tmissing taxa:") == 2
    assert not (output / "matrices").exists()


def test_mapping_that_collapses_paralogs_is_reported_by_qc(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / "all_copies.fa"
    source.write_text(">gene_1\nACGT\n>gene_2\nACGA\n")
    mapping = tmp_path / "taxa.tsv"
    mapping.write_text("gene_1\tSpecies_A\ngene_2\tSpecies_A\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", source, "-st", "nucl", "--notrim", "-m", mapping,
                    "-o", output, "-t", "1")
    assert result.returncode == 2
    assert "No alignments passed QC" in result.stderr
    assert "ASTRAL-Pro3" in result.stderr
    assert (output / "alignments/all_copies.mafft.nucl.aln").exists()
    details = (output / "qc/nucl/raw/alignment.details.tsv").read_text()
    assert "duplicate_ids" in details


def test_rename_checkpoint_is_invalidated_when_mapping_changes(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / "gene.fa"
    source.write_text(">a\nGCTGGT\n>b\nGCCGGA\n")
    mapping = tmp_path / "taxa.tsv"
    mapping.write_text("a\tSpecies_A\nb\tSpecies_B\n")
    output = tmp_path / "results"
    args = ["-i", source, "-m", mapping, "-o", output, "--notrim", "-t", "1"]
    first = invoke(tmp_path, env, *args)
    assert first.returncode == 0, first.stdout + first.stderr
    mapping.write_text("a\tSpecies_C\nb\tSpecies_D\n")
    second = invoke(tmp_path, env, *args)
    assert second.returncode == 0, second.stdout + second.stderr
    assert "rename gene.mafft.codon.aln: skipped" not in second.stdout
    check_formats(output / "matrices/codon/raw", "supermatrix",
                  {"Species_C": "GCTGGT", "Species_D": "GCCGGA"})
    assert Path(env["MSAP_TEST_LOG"]).read_text().splitlines() == ["run"]


def test_codon_iupac_is_masked_before_translation_and_export(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'ambiguous.fa'
    source.write_text('>a\nATGWTAGCTGGT\n>b\nATGATAGCCGGA\n')
    output = tmp_path / 'results'
    result = invoke(tmp_path, env, '-i', source, '-o', output, '--notrim', '-t', '1')
    assert result.returncode == 0, result.stdout + result.stderr
    codon = output / 'matrices/codon/raw/supermatrix.fasta'
    protein = output / 'matrices/prot/raw/supermatrix.fasta'
    assert 'ATGNTAGCTGGT' in codon.read_text()
    assert 'MXAG' in protein.read_text()
    report = (output / 'reports/sequence_changes.tsv').read_text()
    assert '\tnucleotide_ambiguity\t4-4\tW\tN\treplace\t1\n' in report


@pytest.mark.parametrize("seqtype,good,bad", [("nucl", "ACGT", "NNNN"), ("prot", "ACDE", "XXXX")])
def test_qc_filters_genes_and_non_codon_modes(tmp_path, fake_aligner, seqtype, good, bad):
    _, env = fake_aligner
    sources = [tmp_path / "good.fa", tmp_path / "bad.fa"]
    for source, sequence in zip(sources, (good, bad)):
        source.write_text(f">a\n{sequence}\n>b\n{sequence}\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", *sources, "-o", output, "-st", seqtype,
                    "--notrim", "-t", "1", "--max-n-ratio", "0.1", "--max-x-ratio", "0.1")
    assert result.returncode == 0, result.stdout + result.stderr
    check_formats(output / "matrices" / seqtype / "raw", "supermatrix", {"a": good, "b": good})
    assert "bad.mafft" in (output / "qc" / seqtype / "raw/alignment.fail.tsv").read_text()
    assert "good.mafft" in (output / "qc" / seqtype / "raw/alignment.pass.tsv").read_text()
    assert not list((output / "matrices" / seqtype / "raw").glob("codon*"))
    assert not (output / "matrices" / seqtype / "trimmed").exists()
    if seqtype == "prot":
        assert "datatype=protein" in (output / "matrices/prot/raw/supermatrix.nex").read_text().lower()


@pytest.mark.parametrize("failure", ["alignment", "qc", "fourfold"])
def test_failure_does_not_publish_matrices(tmp_path, fake_aligner, failure):
    _, env = fake_aligner
    source = tmp_path / "gene.fa"
    source.write_text(">a\nATGATG\n>b\nATGATG\n")
    output = tmp_path / "results"
    extra = []
    if failure == "alignment":
        env["MSAP_TEST_FAIL"] = "1"
    elif failure == "qc":
        extra = ["--min-length", "100"]
    result = invoke(tmp_path, env, "-i", source, "-o", output, "-t", "1", "--notrim", *extra)
    assert result.returncode != 0
    assert not (output / "matrices").exists()
    if failure == "alignment":
        assert not (output / "qc").exists()
    elif failure == "qc":
        assert "No alignments passed QC" in result.stderr
    else:
        assert "four-fold" in result.stderr


def test_raw_and_trimmed_are_distinct_and_concatenate_multiple_genes(tmp_path, fake_aligner):
    _, env = fake_aligner
    sources = [tmp_path / "gene1.fa", tmp_path / "gene2.fa"]
    sources[0].write_text(">a\nGCTNNNGGT\n>b\nGCCGCAGGA\n")
    sources[1].write_text(">a\nGGT\n>b\nGGA\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", *sources, "-o", output, "-t", "1")
    assert result.returncode == 0, result.stdout + result.stderr
    for seqtype, raw, trimmed in (
        ("codon", {"a": "GCTNNNGGTGGT", "b": "GCCGCAGGAGGA"},
         {"a": "GCTGGTGGT", "b": "GCCGGAGGA"}),
        ("prot", {"a": "AXGG", "b": "AAGG"}, {"a": "AGG", "b": "AGG"}),
    ):
        check_formats(output / "matrices" / seqtype / "raw", "supermatrix", raw)
        check_formats(output / "matrices" / seqtype / "trimmed", "supermatrix", trimmed)
    for variant in ("raw", "trimmed"):
        check_formats(output / "matrices/codon" / variant, "fourfold", {"a": "TTT", "b": "CAA"})


def test_codon_notrim_exports_protein_without_trimmed_groups(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / "gene.fa"
    source.write_text(">a\nGCTGGT\n>b\nGCCGGA\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", source, "-o", output, "-t", "1", "--notrim")
    assert result.returncode == 0, result.stdout + result.stderr
    check_formats(output / "matrices/prot/raw", "supermatrix", {"a": "AG", "b": "AG"})
    check_formats(output / "matrices/codon/raw", "supermatrix", {"a": "GCTGGT", "b": "GCCGGA"})
    assert not list(output.glob("matrices/*/trimmed"))
    assert not list(output.glob("qc/*/trimmed"))


@pytest.mark.parametrize("seqtype,ambiguous", [("nucl", "N"), ("prot", "X")])
def test_qc_is_independent_for_raw_and_trimmed(tmp_path, fake_aligner, seqtype, ambiguous):
    _, env = fake_aligner
    good, repaired = tmp_path / "good.fa", tmp_path / "repaired.fa"
    good.write_text(">a\nACGT\n>b\nACGT\n")
    repaired.write_text(f">a\n{ambiguous}AAA\n>b\n{ambiguous}AAA\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", good, repaired, "-o", output, "-st", seqtype,
                    "-t", "1", "--max-n-ratio", "0.1", "--max-x-ratio", "0.1")
    assert result.returncode == 0, result.stdout + result.stderr
    check_formats(output / "matrices" / seqtype / "raw", "supermatrix", {"a": "ACGT", "b": "ACGT"})
    check_formats(output / "matrices" / seqtype / "trimmed", "supermatrix", {"a": "ACGTAAA", "b": "ACGTAAA"})
    assert "repaired.mafft" in (output / "qc" / seqtype / "raw/alignment.fail.tsv").read_text()
    assert "repaired.mafft" in (output / "qc" / seqtype / "trimmed/alignment.pass.tsv").read_text()


def test_late_protein_qc_failure_does_not_publish_codon_matrices(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / "gene.fa"
    source.write_text(">a\nGCTGGT\n>b\nGCCGGA\n")
    output = tmp_path / "results"
    result = invoke(tmp_path, env, "-i", source, "-o", output, "-t", "1",
                    "--notrim", "--min-length", "6")
    assert result.returncode != 0
    assert "No alignments passed QC for prot/raw" in result.stderr
    assert (output / "qc/codon/raw/alignment.pass.tsv").read_text().strip()
    assert not (output / "matrices").exists()
