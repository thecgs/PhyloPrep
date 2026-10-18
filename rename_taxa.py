#!/usr/bin/env python3
"""Rename FASTA sequence IDs using a mapping table and safe-name rules."""
from __future__ import annotations
import argparse, csv, sys
import contextlib, hashlib, io, os
from pathlib import Path
from Bio import SeqIO
from msap_io import atomic_output, check_output_path, file_signature

def normalized_taxon_id(identifier):
    return identifier.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")


def file_stamp(path):
    path = Path(path)
    stat = path.stat()
    return (str(path.resolve()), stat.st_size, stat.st_mtime_ns,
            stat.st_ctime_ns, stat.st_dev, stat.st_ino)


class FileFingerprintCache:
    """Reuse verified digests only within one run while file metadata matches."""
    def __init__(self):
        self.entries = {}

    def peek(self, path):
        stamp = file_stamp(path)
        entry = self.entries.get(stamp[0])
        return entry[1] if entry and entry[0] == stamp else None

    def remember(self, path, sha256):
        stamp = file_stamp(path)
        signature = dict(path=stamp[0], size=stamp[1], mtime_ns=stamp[2], sha256=sha256)
        self.entries[stamp[0]] = (stamp, signature)
        return signature

    def signature(self, path):
        cached = self.peek(path)
        if cached is not None:
            return cached
        before = file_stamp(path)
        signature = file_signature(path)
        if file_stamp(path) != before:
            raise ValueError(f"File changed while fingerprinting: {path}")
        self.entries[before[0]] = (before, signature)
        return signature

    def forget(self, path):
        self.entries.pop(str(Path(path).resolve()), None)

    def clear(self):
        self.entries.clear()


class _DigestReader(io.RawIOBase):
    def __init__(self, source, digest):
        self.source, self.digest = source, digest

    def readable(self):
        return True

    def readinto(self, buffer):
        count = self.source.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count


@contextlib.contextmanager
def _hashed_text_input(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as raw, io.TextIOWrapper(
            io.BufferedReader(_DigestReader(raw, digest)), encoding='utf-8') as source:
        yield source, digest


class _DigestWriter:
    def __init__(self, target):
        self.target, self.digest = target, hashlib.sha256()

    def write(self, text):
        count = self.target.write(text)
        # atomic_output uses UTF-8 with native newline translation.
        self.digest.update(text.replace('\n', os.linesep).encode('utf-8'))
        return count


def _read_mapping_rows(handle):
    mapping = {}
    for row in csv.reader(handle, delimiter="\t"):
        if not row or row[0].startswith("#"):
            continue
        if len(row) < 2:
            raise ValueError("Mapping table requires two tab-separated columns: old_id and new_id")
        old_id, new_id = normalized_taxon_id(row[0]), normalized_taxon_id(row[1])
        if old_id in mapping and mapping[old_id] != new_id:
            raise ValueError(f"Conflicting mapping after OrthoFinder normalization: {row[0]}")
        mapping[old_id] = new_id
    return mapping


def read_mapping(mapping_file):
    with open(mapping_file, encoding="utf-8", newline="") as handle:
        return _read_mapping_rows(handle)


def read_mapping_snapshot(mapping_file):
    """Hash the exact mapping bytes while parsing, without extra file reads."""
    before = file_stamp(mapping_file)
    with _hashed_text_input(mapping_file) as (source, digest):
        mapping = _read_mapping_rows(source)
    if file_stamp(mapping_file) != before:
        raise ValueError("Mapping changed while being read; rerun with a stable mapping file")
    signature = dict(path=before[0], size=before[1], mtime_ns=before[2], sha256=digest.hexdigest())
    return mapping, signature


def rename_fasta(infile, outfile, mapping_file=None, map_output=None, allow_duplicate_ids=False,
                 mapping=None, fingerprints=None, input_sha256=None):
    if mapping_file is not None and mapping is not None:
        raise ValueError("Provide either mapping_file or mapping, not both")
    inputs = [infile] + ([mapping_file] if mapping_file else [])
    check_output_path(outfile, inputs)
    check_output_path(map_output, inputs + [outfile])
    mapping = read_mapping(mapping_file) if mapping_file else (mapping or {})
    seen = set()
    map_rows = [] if map_output else None

    def renamed_records(source):
        for record in SeqIO.parse(source, "fasta"):
            old = normalized_taxon_id(record.id)
            new = mapping.get(old, old)
            if not new:
                raise ValueError(f"{old}: renamed ID is empty")
            if any(char.isspace() for char in new):
                raise ValueError(f'{old}: FASTA IDs must not contain whitespace, even with --keep-special: {new!r}')
            if new in seen and not allow_duplicate_ids:
                raise ValueError(f"Duplicate ID after renaming: {new}")
            seen.add(new)
            record.id = record.name = new
            record.description = ""
            if map_rows is not None:
                map_rows.append((old, new))
            yield record

    with contextlib.ExitStack() as stack:
        if fingerprints is not None and input_sha256 is None:
            source, input_digest = stack.enter_context(_hashed_text_input(infile))
        else:
            source = stack.enter_context(open(infile))
            input_digest = None
        handle = stack.enter_context(atomic_output(outfile, inputs=[infile]))
        writer = _DigestWriter(handle) if fingerprints is not None and not str(outfile).endswith('.gz') else handle
        if not SeqIO.write(renamed_records(source), writer, "fasta"):
            raise ValueError("No sequences found in input FASTA.")
    if fingerprints is not None:
        fingerprints['source'] = input_sha256 or input_digest.hexdigest()
        fingerprints['target'] = writer.digest.hexdigest() if isinstance(writer, _DigestWriter) else file_signature(outfile)['sha256']
    if map_output:
        with atomic_output(map_output) as handle:
            handle.write("original_id\trenamed_id\n")
            handle.writelines(f"{old}\t{new}\n" for old, new in map_rows)

def main(argv=None):
    parser = argparse.ArgumentParser(description="Rename FASTA sequence IDs for phylogenetic software.", add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="""Examples:
  rename_taxa.py -i input.fasta -o renamed.fasta -m taxa.tsv
  rename_taxa.py -i input.fasta -o renamed.fasta --map-output renamed.tsv

Mapping format: tab-separated original_id and renamed_id, one pair per line.
Special characters are replaced with underscores by default.""")
    required=parser.add_argument_group("required arguments"); optional=parser.add_argument_group("renaming options")
    required.add_argument("-i","--input",required=True,metavar="FASTA",help="Input FASTA file.")
    required.add_argument("-o","--output",required=True,metavar="FASTA",help="Output FASTA file.")
    optional.add_argument("-m","--mapping","--map",metavar="TSV",help="Two-column ID mapping table.")
    optional.add_argument("--map-output",metavar="TSV",help="Write original-to-renamed IDs.")
    optional.add_argument("--allow-duplicate-ids", action="store_true",
                        help="Write duplicate renamed IDs for a downstream validator; normally duplicates are an error.")
    optional.add_argument("--keep-special", action="store_true", help="Deprecated compatibility option; OrthoFinder substitutions still apply.")
    optional.add_argument("-h","--help",action="help",help="Show this help message and exit.")
    optional.add_argument("-v","--version",action="version",version="v1.0.0",help="Show program's version number and exit.")
    args=parser.parse_args(argv)
    try: rename_fasta(args.input,args.output,args.mapping,args.map_output,args.allow_duplicate_ids)
    except (OSError,ValueError) as error: print(f"rename_taxa: {error}",file=sys.stderr); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
