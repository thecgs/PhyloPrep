import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "calulate_4dtv_and_correction.py"


def run_4dtv(tmp_path, fasta, output_name="result.tsv"):
    input_file = tmp_path / "alignment.fasta"
    output_file = tmp_path / output_name
    input_file.write_text(fasta)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(input_file), "-o", str(output_file)],
        capture_output=True,
        text=True,
    )
    return result, output_file


def test_no_fourfold_sites_returns_na(tmp_path):
    result, output_file = run_4dtv(tmp_path, ">seq1\nATG\n>seq2\nATG\n")

    assert result.returncode == 0, result.stderr
    row = output_file.read_text().splitlines()[1].split("\t")
    assert row[3:5] == ["NA", "NA"]
    assert row[5:7] == ["0", "0"]


def test_rejects_fewer_than_two_sequences_without_overwriting_output(tmp_path):
    output = tmp_path / "result.tsv"
    output.write_text("old result\n")
    result, output_file = run_4dtv(tmp_path, ">seq1\nGCT\n")

    assert result.returncode != 0
    assert "expected exactly two sequences" in result.stderr
    assert output_file.read_text() == "old result\n"


def test_rejects_alignment_length_not_divisible_by_three(tmp_path):
    result, output_file = run_4dtv(tmp_path, ">seq1\nGCTA\n>seq2\nGCTA\n")

    assert result.returncode != 0
    assert "not divisible by 3" in result.stderr
    assert not output_file.exists()
