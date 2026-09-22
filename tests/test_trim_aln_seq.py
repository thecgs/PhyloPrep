import subprocess
import sys
import os
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "trimAlnSeq.py"
MSAP_SCRIPT = Path(__file__).parents[1] / "MSAP.py"
sys.path.insert(0, str(MSAP_SCRIPT.parent))
import MSAP


def run_trim(tmp_path, seqtype, sequences):
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(
        "".join(f">seq{i}\n{sequence}\n" for i, sequence in enumerate(sequences))
    )
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "-st",
            seqtype,
        ],
        check=True,
    )
    return output_file.read_text()


@pytest.mark.parametrize(
    ("seqtype", "ambiguous_base"),
    [
        ("nucl", "-"),
        ("nucl", "N"),
        ("prot", "X"),
    ],
)
def test_default_threshold_retains_a_site_at_twenty_percent(
    tmp_path, seqtype, ambiguous_base
):
    """The default 0.2 threshold retains a site observed in one of five sequences."""
    output = run_trim(
        tmp_path,
        seqtype,
        ["ACT", f"{ambiguous_base}CT", "ACT", "ACT", "ACT"],
    )

    assert output == "".join(
        f">seq{i}\n{sequence}\n"
        for i, sequence in enumerate(
            ["ACT", f"{ambiguous_base}CT", "ACT", "ACT", "ACT"]
        )
    )


def test_codon_site_with_gap_and_n_is_filtered_by_n_threshold(tmp_path):
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"
    input_file.write_text(">seq0\nATG\n>seq1\nN--\n")

    output_file.write_text("previous matrix\n")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "-st",
            "codon",
            "-G",
            "1",
            "-N",
            "0",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "No alignment sites remain" in result.stderr
    assert output_file.read_text() == "previous matrix\n"


def test_msap_can_use_trimal_with_a_path_containing_spaces(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    mafft = bin_dir / "mafft"
    mafft.write_text(
        "#!/bin/sh\n"
        "for arg do input=\"$arg\"; done\n"
        "cat \"$input\"\n"
    )
    trimal = bin_dir / "trimal"
    trimal.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' \"$@\" > trimal.args\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  case \"$1\" in\n"
        "    -in) infile=\"$2\"; shift 2 ;;\n"
        "    -out) outfile=\"$2\"; shift 2 ;;\n"
        "    *) shift ;;\n"
        "  esac\n"
        "done\n"
        "cp \"$infile\" \"$outfile\"\n"
        "printf '#ColumnsMap\\t0,1,2\\n'\n"
    )
    mafft.chmod(0o755)
    trimal.chmod(0o755)

    input_file = tmp_path / "input with spaces.fasta"
    input_file.write_text(">seq1\nACT\n>seq2\nACT\n")
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}

    subprocess.run(
        [
            sys.executable,
            str(MSAP_SCRIPT),
            "-i",
            str(input_file),
            "-st",
            "nucl",
            "--trim-software",
            "trimal",
        ],
        check=True,
        cwd=tmp_path,
        env=env,
    )

    output = tmp_path / "input with spaces.mafft.nucl.trimal.aln"
    assert output.read_text() == input_file.read_text()
    assert "-automated1" in (tmp_path / "trimal.args").read_text().splitlines()

    subprocess.run(
        [
            sys.executable,
            str(MSAP_SCRIPT),
            "-i",
            str(input_file),
            "-st",
            "nucl",
            "--trim-software",
            "trimal",
            "--trimal-args",
            "-gt",
            "0.9",
        ],
        check=True,
        cwd=tmp_path,
        env=env,
    )

    trimal_args = (tmp_path / "trimal.args").read_text().splitlines()
    assert trimal_args[-2:] == ["-gt", "0.9"]
    assert "-automated1" not in trimal_args

    trimal.unlink()
    untrimmed_input = tmp_path / "untrimmed input.fasta"
    untrimmed_input.write_text(">seq1\nACT\n>seq2\nACT\n")
    subprocess.run(
        [
            sys.executable,
            str(MSAP_SCRIPT),
            "-i",
            str(untrimmed_input),
            "-st",
            "nucl",
            "--trim-software",
            "trimal",
            "--notrim",
        ],
        check=True,
        cwd=tmp_path,
        env=env,
    )

    assert (tmp_path / "untrimmed input.mafft.nucl.aln").read_text() == untrimmed_input.read_text()
    assert not (tmp_path / "untrimmed input.mafft.nucl.trimal.aln").exists()


@pytest.mark.parametrize(
    ("seqtype", "fasta", "error"),
    [
        ("nucl", ">sp1\nACT\n>sp1\nACT\n", "Duplicate sequence ID: sp1"),
    ],
)
def test_msap_validates_input_before_checking_external_dependencies(
    tmp_path, seqtype, fasta, error
):
    input_file = tmp_path / "invalid.fasta"
    input_file.write_text(fasta)

    result = subprocess.run(
        [sys.executable, str(MSAP_SCRIPT), "-i", str(input_file), "-st", seqtype],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert error in result.stderr


def _fake_trimal(path, output_fasta, columns):
    """Write a minimal trimAl replacement for codon-transfer tests."""
    path.write_text(
        "#!/bin/sh\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  case \"$1\" in\n"
        "    -out) outfile=\"$2\"; shift 2 ;;\n"
        "    *) shift ;;\n"
        "  esac\n"
        "done\n"
        f"cat > \"$outfile\" <<'EOF'\n{output_fasta}EOF\n"
        f"printf '#ColumnsMap\\t{columns}\\n'\n"
    )
    path.chmod(0o755)


def test_trimal_codon_columns_are_transferred_as_whole_triplets(tmp_path):
    protein = tmp_path / "protein.aln"
    codon = tmp_path / "codon.aln"
    protein_out = tmp_path / "protein.trim.aln"
    codon_out = tmp_path / "codon.trim.aln"
    trimal = tmp_path / "trimal"
    protein.write_text(">a\nMKA\n>b\nMKA\n")
    codon.write_text(">a\nATGAAAGCT\n>b\nATGAAGGCC\n")
    _fake_trimal(trimal, ">a\nMA\n>b\nMA\n", "0,2")

    MSAP.trim_codon_with_trimal(
        str(trimal), [], str(protein), str(codon), str(protein_out), str(codon_out)
    )

    assert protein_out.read_text() == ">a\nMA\n>b\nMA\n"
    assert codon_out.read_text() == ">a\nATGGCT\n>b\nATGGCC\n"
    assert all(len(line) % 3 == 0 for line in codon_out.read_text().splitlines()
               if line and not line.startswith(">"))


def test_trimal_codon_rejects_changed_sequence_set_without_writing_outputs(tmp_path):
    protein = tmp_path / "protein.aln"
    codon = tmp_path / "codon.aln"
    protein_out = tmp_path / "protein.trim.aln"
    codon_out = tmp_path / "codon.trim.aln"
    trimal = tmp_path / "trimal"
    protein.write_text(">a\nMKA\n>b\nMKA\n")
    codon.write_text(">a\nATGAAAGCT\n>b\nATGAAGGCC\n")
    protein_out.write_text("previous protein\n")
    codon_out.write_text("previous codon\n")
    _fake_trimal(trimal, ">a\nMA\n", "0,2")

    with pytest.raises(ValueError, match="changed the sequence set"):
        MSAP.trim_codon_with_trimal(
            str(trimal), [], str(protein), str(codon), str(protein_out), str(codon_out)
        )

    assert protein_out.read_text() == "previous protein\n"
    assert codon_out.read_text() == "previous codon\n"


@pytest.mark.parametrize("option", ["-in", "-out", "-fasta", "-backtrans", "-colnumbering"])
def test_trimal_codon_rejects_arguments_that_override_its_mapping_contract(tmp_path, option):
    protein = tmp_path / "protein.aln"
    codon = tmp_path / "codon.aln"
    protein.write_text(">a\nMKA\n>b\nMKA\n")
    codon.write_text(">a\nATGAAAGCT\n>b\nATGAAGGCC\n")

    with pytest.raises(ValueError, match="managed automatically"):
        MSAP.trim_codon_with_trimal(
            "unused", [option], str(protein), str(codon),
            str(tmp_path / "protein.out"), str(tmp_path / "codon.out")
        )
