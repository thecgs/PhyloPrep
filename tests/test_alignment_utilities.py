import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import AlignIO


REPO = Path(__file__).parents[1]
CONVERTERS = ["fasta2phy.py", "fasta2nex.py", "fasta2axt.py"]


def run(script, *args, cwd):
    return subprocess.run([sys.executable, str(REPO / script), *map(str, args)],
                          cwd=cwd, capture_output=True, text=True)


@pytest.mark.parametrize("script", CONVERTERS)
@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_conversion_cannot_destroy_input(tmp_path, script, alias):
    source = tmp_path / "in.fasta"
    original = ">a\nACT\n>b\nACG\n"
    source.write_text(original)
    output = source if alias == "same" else tmp_path / "alias"
    if alias == "symlink":
        output.symlink_to(source)
    elif alias == "hardlink":
        os.link(source, output)
    result = run(script, "-i", source, "-o", output, cwd=tmp_path)
    assert result.returncode != 0
    assert "different files" in result.stderr
    assert source.read_text() == original


@pytest.mark.parametrize("script", CONVERTERS)
@pytest.mark.parametrize("fasta", ["", ">a\n\n>b\n\n", ">a\nACT\n>b\nAC\n", ">a\nACT\n>a\nACG\n"])
def test_invalid_conversion_preserves_existing_result(tmp_path, script, fasta):
    source, output = tmp_path / "in.fasta", tmp_path / "result"
    source.write_text(fasta)
    output.write_text("previous result\n")
    result = run(script, "-i", source, "-o", output, cwd=tmp_path)
    assert result.returncode != 0
    assert output.read_text() == "previous result\n"


@pytest.mark.parametrize("script,fmt", [("fasta2phy.py", "phylip-relaxed"), ("fasta2nex.py", "nexus")])
def test_conversion_roundtrip(tmp_path, script, fmt):
    source, output = tmp_path / "in.fasta", tmp_path / "result"
    source.write_text(">long_taxon_name_1\nACGT-\n>long_taxon_name_2\nAC-TN\n")
    result = run(script, "-i", source, "-o", output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    alignment = AlignIO.read(output, fmt)
    assert [(r.id, str(r.seq)) for r in alignment] == [
        ("long_taxon_name_1", "ACGT-"), ("long_taxon_name_2", "AC-TN")]


def test_axt_retains_kaks_dialect_and_requires_pair(tmp_path):
    source = tmp_path / "in.fasta"
    source.write_text(">a\nATG---\n>b\nATGAAA\n")
    result = run("fasta2axt.py", "-i", source, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ">a|b\nATG---\nATGAAA\n"
    source.write_text(source.read_text() + ">c\nATGCCC\n")
    result = run("fasta2axt.py", "-i", source, cwd=tmp_path)
    assert result.returncode != 0
    assert "exactly 2" in result.stderr


@pytest.mark.parametrize("fasta", [">a\nATGA\n>b\nATGA\n", ">a\nATG\n>b\nATGAAA\n", ""])
def test_invalid_codon_split_preserves_all_outputs(tmp_path, fasta):
    source = tmp_path / "in.fasta"
    source.write_text(fasta)
    outputs = [tmp_path / f"codon{pos}.in.fasta" for pos in ("1st", "2nd", "3rd")]
    for output in outputs:
        output.write_text("previous result\n")
    result = run("split_codon_seqence_alignment.py", "-i", source, cwd=tmp_path)
    assert result.returncode != 0
    assert all(output.read_text() == "previous result\n" for output in outputs)


def test_codon_split_retains_positions(tmp_path):
    source = tmp_path / "in.fasta"
    source.write_text(">a\nATGGCT\n>b\nATG---\n")
    result = run("split_codon_seqence_alignment.py", "-i", source, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "codon1st.in.fasta").read_text() == ">a\nAG\n>b\nA-\n"
    assert (tmp_path / "codon2nd.in.fasta").read_text() == ">a\nTC\n>b\nT-\n"
    assert (tmp_path / "codon3rd.in.fasta").read_text() == ">a\nGT\n>b\nG-\n"
