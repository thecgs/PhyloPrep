import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "sequence_qc.py"


def run(tmp_path, fasta, *extra):
    source = tmp_path / "alignment.fasta"
    output = tmp_path / "report.tsv"
    source.write_text(fasta)
    result = subprocess.run([sys.executable, str(SCRIPT), "-i", str(source), "-o", str(output), *extra], capture_output=True, text=True)
    return result, output


def test_reports_alignment_and_full_gap_columns(tmp_path):
    result, output = run(tmp_path, ">a\nACT-\n>b\nAC--\n", "-st", "nucl")
    assert result.returncode == 0, result.stderr
    text = output.read_text()
    assert "File\tTaxa_ID\tSequence_length\tN_count\tGap_count" in text
    assert "Terminal_stop" not in text.splitlines()[0]
    assert "\ta\t4\t0\t1" in text


def test_reports_codon_frame_and_stop(tmp_path):
    result, output = run(tmp_path, ">a\nATGTAA\n>b\nATGTAA\n", "-st", "codon")
    assert result.returncode == 0, result.stderr
    assert "\tTAA\tNA" in output.read_text()


def test_protein_report_has_x_without_n(tmp_path):
    result, output = run(tmp_path, ">a\nACX\n>b\nACD\n", "-st", "prot")
    assert result.returncode == 0, result.stderr
    header = output.read_text().splitlines()[0]
    assert header == "File\tTaxa_ID\tSequence_length\tX_count\tGap_count"
    assert "N_count" not in header


def test_html_output_and_invalid_input(tmp_path):
    result, output = run(tmp_path, ">a\nAC!\n>b\nACG\n", "-st", "nucl", "--format", "html")
    assert result.returncode == 0
    assert "<table>" in output.read_text()
