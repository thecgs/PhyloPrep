"""Shared validation and safe output for alignment utilities."""

import os
import sys
import tempfile
import gzip
from contextlib import contextmanager
from pathlib import Path

from Bio import SeqIO

def expand_input_paths(items):
    """Expand FASTA paths and arbitrary-extension text path lists."""
    paths = []
    for item in items:
        path = Path(item)
        if path.is_file():
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            first = next((line.strip() for line in lines if line.strip() and not line.lstrip().startswith(("#", ";"))), "")
            if not first.startswith(">"):
                paths.extend(line.strip() for line in lines
                             if line.strip() and not line.lstrip().startswith(("#", ";")))
                continue
        paths.append(str(path))
    return paths


def read_fasta_alignment(infile, codon=False, count=None):
    records = list(SeqIO.parse(infile, "fasta"))
    if not records:
        raise ValueError("No sequences found in input FASTA.")
    if count is not None and len(records) != count:
        raise ValueError(f"Expected exactly {count} sequences, found {len(records)}.")
    length = len(records[0].seq)
    ids = set()
    for record in records:
        if record.id in ids:
            raise ValueError(f"Duplicate sequence ID: {record.id}")
        ids.add(record.id)
        if not record.seq:
            raise ValueError(f"{record.id}: sequence is empty")
        if len(record.seq) != length:
            raise ValueError(
                f"Alignment length mismatch: {record.id} has length "
                f"{len(record.seq)}, expected {length}"
            )
    if codon and length % 3:
        raise ValueError(f"Codon alignment length ({length}) is not divisible by 3.")
    return records


def check_output_path(outfile, inputs=()):
    if outfile is None:
        return
    output = Path(outfile)
    output.parent.mkdir(parents=True, exist_ok=True)
    for infile in inputs:
        source = Path(infile)
        if output.resolve() == source.resolve() or (
            output.exists() and source.exists() and os.path.samefile(output, source)
        ):
            raise ValueError(f"Input and output must be different files: {outfile}")


@contextmanager
def atomic_output(outfile, inputs=()):
    """Keep existing files intact on validation/writer errors; never close stdout."""
    check_output_path(outfile, inputs)
    if outfile is None:
        yield sys.stdout
        return
    output = Path(outfile)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=output.absolute().parent,
                                         prefix=".msap-", suffix=".tmp", delete=False) as raw:
            temporary = raw.name
        opener = gzip.open if output.name.endswith(".gz") else open
        with opener(temporary, "wt", encoding="utf-8") as handle:
            yield handle
        os.replace(temporary, output)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
