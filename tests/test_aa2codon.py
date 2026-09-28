import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "AA2Codon.py"


def run_aa2codon(tmp_path, cds, protein, output_name="output.fasta", extra_args=()):
    cds_file = tmp_path / "cds.fasta"
    protein_file = tmp_path / "protein.fasta"
    output_file = tmp_path / output_name
    cds_file.write_text(cds)
    protein_file.write_text(protein)
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-c",
            str(cds_file),
            "-p",
            str(protein_file),
            "-o",
            str(output_file),
            *extra_args,
        ],
        capture_output=True,
        text=True,
    )
    return result, output_file


def test_converts_a_gapped_protein_alignment_after_validation(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nATGAAA\n",
        ">sp1\nM-K\n",
        extra_args=("-g", "1"),
    )

    assert result.returncode == 0, result.stderr
    assert output_file.read_text() == ">sp1\nATG---AAA\n"


def test_standard_code_does_not_treat_gtg_as_an_initiator(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nGTGAAA\n",
        ">sp1\nV-K\n",
        extra_args=("-g", "1"),
    )

    assert result.returncode == 0, result.stderr
    assert output_file.read_text() == ">sp1\nGTG---AAA\n"


def test_alternative_start_codon_matches_msap_initiator_methionine(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nGTGAAA\n",
        ">sp1\nM-K\n",
        extra_args=("-g", "11"),
    )

    assert result.returncode == 0, result.stderr
    assert output_file.read_text() == ">sp1\nGTG---AAA\n"


def test_rejects_mismatched_id_sets_without_creating_output(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nATG\n",
        ">sp2\nM\n",
    )

    assert result.returncode != 0
    assert "CDS/protein alignment ID mismatch" in result.stderr
    assert not output_file.exists()


def test_rejects_residue_and_codon_count_mismatch_without_creating_output(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nATGAAA\n",
        ">sp1\nM-\n",
    )

    assert result.returncode != 0
    assert "protein alignment has 1 non-gap residues, but CDS has 2 codons" in result.stderr
    assert not output_file.exists()


def test_rejects_translation_mismatch_without_creating_output(tmp_path):
    result, output_file = run_aa2codon(
        tmp_path,
        ">sp1\nATG\n",
        ">sp1\nK\n",
        extra_args=("-g", "1"),
    )

    assert result.returncode != 0
    assert "translates to M, but the protein alignment contains K" in result.stderr
    assert not output_file.exists()
