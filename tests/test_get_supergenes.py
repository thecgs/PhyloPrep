import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "get_supergenes.py"


def write_alignment(path, records):
    path.write_text("".join(f">{taxon}\n{sequence}\n" for taxon, sequence in records))


def run_supergenes(tmp_path, missing_taxa):
    complete_gene = tmp_path / "complete.fasta"
    incomplete_gene = tmp_path / "incomplete.fasta"
    prefix = tmp_path / missing_taxa
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
            "--missing-taxa",
            missing_taxa,
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
    partition_config = (prefix.parent / f"{prefix.name}_partition_finder.cfg").read_text()
    assert f"alignment = {prefix.name}.phy;" in partition_config


def test_pad_gaps_preserves_a_gene_missing_a_taxon(tmp_path):
    prefix = run_supergenes(tmp_path, "pad-gaps")

    assert (prefix.parent / f"{prefix.name}.fasta").read_text() == (
        ">sp1\nAAACCC\n>sp2\nAAACCC\n>sp3\nAAA---\n"
    )
    report = (prefix.parent / f"{prefix.name}.report.tsv").read_text()
    assert "incomplete.fasta\tpadded\tmissing taxa: sp3" in report
    assert (prefix.parent / f"{prefix.name}.phy").read_text() == (
        "3 6\nsp1 AAACCC\nsp2 AAACCC\nsp3 AAA---\n"
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
