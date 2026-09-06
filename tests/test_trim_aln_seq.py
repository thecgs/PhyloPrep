import subprocess
import sys
import os
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "trimAlnSeq.py"
MSAP_SCRIPT = Path(__file__).parents[1] / "MSAP.py"


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

    subprocess.run(
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
        check=True,
    )

    assert output_file.read_text() == ">seq0\n\n>seq1\n\n"


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
        ("codon", ">sp1\nATGA\n>sp2\nATG\n", "CDS length 4 is not divisible by 3"),
    ],
)
def test_msap_validates_input_before_checking_external_dependencies(
    tmp_path, seqtype, fasta, error
):
    input_file = tmp_path / "invalid.fasta"
    input_file.write_text(fasta)

    result = subprocess.run(
        [sys.executable, str(MSAP_SCRIPT), "-i", str(input_file), "-st", seqtype],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert error in result.stderr
