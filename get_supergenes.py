#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import re
import tempfile
from pathlib import Path
from collections import OrderedDict
from contextlib import ExitStack
from Bio import SeqIO
from msap_io import atomic_output, expand_input_paths, check_output_path, NoGenesIncludedError, file_signature

def tidy_name(file):
    # Strip only the current PhyloPrep alignment suffix.
    file = re.sub(
        r'\.(?:mafft|muscle|clustalw2|prank|macse)\.(?:nucl|codon|prot)(?:\.trimmed)?\.aln$',
        '', os.path.basename(file),
    ).replace('.', '_').replace('-', '_')

    file = re.sub(r'[^A-Za-z0-9_]', '_', file)
    if not file or not file[0].isalpha():
        file = 'gene_' + file
    return file


def read_alignment(file):
    """Read one FASTA alignment and validate its IDs and sequence lengths."""
    sequences = {}
    expected_length = None
    for record in SeqIO.parse(file, "fasta"):
        if record.id in sequences:
            raise ValueError(f"{file}: duplicate sequence ID: {record.id}")
        if expected_length is None:
            expected_length = len(record.seq)
        if len(record.seq) != expected_length:
            raise ValueError(f"{file}: {record.id} has length {len(record.seq)}, expected {expected_length}")
        sequences[record.id] = str(record.seq)
    if expected_length is None:
        return None, "empty alignment"
    if expected_length == 0:
        return None, "zero-column alignment"
    return sequences, None


def scan_alignment(file):
    """Validate one record at a time, retaining identifiers but no sequences."""
    taxa = []
    seen = set()
    length = None
    for record in SeqIO.parse(file, "fasta"):
        if record.id in seen:
            raise ValueError(f"{file}: duplicate sequence ID: {record.id}")
        seen.add(record.id)
        taxa.append(record.id)
        if length is None:
            length = len(record.seq)
        if len(record.seq) != length:
            raise ValueError(f"{file}: {record.id} has length {len(record.seq)}, expected {length}")
    return taxa, length


class TaxonSpool:
    """Append taxon fragments to temporary files with a bounded handle cache."""
    chunk_size = 65536

    def __init__(self, directory, taxa):
        self.paths = {taxon: Path(directory) / str(i) for i, taxon in enumerate(taxa)}
        self.handles = OrderedDict()

    def append(self, taxon, sequence, length):
        handle = self.handles.pop(taxon, None)
        if handle is None:
            if len(self.handles) >= 32:
                _, oldest = self.handles.popitem(last=False)
                oldest.close()
            handle = self.paths[taxon].open('a', encoding='utf-8')
        self.handles[taxon] = handle
        if sequence is not None:
            handle.write(sequence)
        else:
            # Do not allocate an entire missing-gene sequence.
            while length:
                count = min(length, self.chunk_size)
                handle.write('-' * count)
                length -= count

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()

    def write_outputs(self, prefix, taxa, length, write_phylip):
        self.close()
        with ExitStack() as stack:
            fasta = stack.enter_context(atomic_output(prefix + '.fasta'))
            phylip = stack.enter_context(atomic_output(prefix + '.phy')) if write_phylip else None
            if phylip:
                phylip.write(f'{len(taxa)} {length}\n')
            for taxon in taxa:
                fasta.write(f'>{taxon}\n')
                if phylip:
                    phylip.write(taxon + ' ')
                with self.paths[taxon].open(encoding='utf-8') as source:
                    while True:
                        chunk = source.read(self.chunk_size)
                        if not chunk:
                            break
                        fasta.write(chunk)
                        if phylip:
                            phylip.write(chunk)
                fasta.write('\n')
                if phylip:
                    phylip.write('\n')
                # In Pyodide these files also occupy memory; release each one
                # as its final outputs are written instead of keeping all spools.
                self.paths[taxon].unlink()


def get_supergenes(infiles, prefix, missing_taxa="skip-gene", input_lists=(), write_phylip=True):
    """Create a supermatrix while handling taxa missing from individual genes."""
    infiles = list(infiles)
    destinations = [prefix + suffix for suffix in
                    (".report.tsv", "_partition_finder.cfg", "_part_iqtree.txt", ".fasta")]
    if write_phylip:
        destinations.append(prefix + ".phy")
    for i, destination in enumerate(destinations):
        check_output_path(destination, infiles + list(input_lists) + destinations[:i])
    alignments = []
    taxa_order = []
    taxa_seen = set()
    report_rows = []

    # First pass retains only file metadata and the global taxon order.
    for file in infiles:
        if os.path.getsize(file) == 0:
            report_rows.append((file, "skipped", "empty file"))
            continue
        signature = file_signature(file)
        taxa, length = scan_alignment(file)
        if signature != file_signature(file):
            raise ValueError(f"Input changed while scanning: {file}")
        if not taxa or not length:
            report_rows.append((file, "skipped", "empty alignment" if not taxa else "zero-column alignment"))
            continue
        for taxon in taxa:
            if taxon not in taxa_seen:
                taxa_order.append(taxon)
                taxa_seen.add(taxon)
        alignments.append((file, length, signature))

    # A single locus is loaded at a time in the second pass. Spools use actual
    # disk on the CLI and reduce duplicated Python strings in browser MEMFS.
    with tempfile.TemporaryDirectory(prefix='.supergenes-', dir=Path(prefix).absolute().parent) as directory:
        spool = TaxonSpool(directory, taxa_order)
        try:
            selected_alignments = []
            for file, alignment_length, signature in alignments:
                if signature != file_signature(file):
                    raise ValueError(f"Input changed before concatenation: {file}")
                sequences, error = read_alignment(file)
                if error or signature != file_signature(file):
                    raise ValueError(f"Input changed during concatenation: {file}")
                missing = [taxon for taxon in taxa_order if taxon not in sequences]
                if missing and missing_taxa == "skip-gene":
                    reason = "missing taxa: " + ", ".join(missing)
                    report_rows.append((file, "skipped", reason))
                    print(f"Skipping {file}: {reason}", file=sys.stderr)
                    del sequences
                    continue
                report_rows.append((file, "padded" if missing else "included",
                                    "missing taxa: " + ", ".join(missing) if missing else ""))
                selected_alignments.append((file, alignment_length))
                for taxon in taxa_order:
                    spool.append(taxon, sequences.get(taxon), alignment_length)
                del sequences
            spool.close()
            _write_supermatrix_products(prefix, report_rows, selected_alignments,
                                        taxa_order, spool, write_phylip)
        finally:
            spool.close()
    return None


def _write_supermatrix_products(prefix, report_rows, selected_alignments, taxa_order, spool, write_phylip):
    report_file = prefix + ".report.tsv"
    with atomic_output(report_file) as report:
        print("gene\tstatus\tdetails", file=report)
        for file, status, details in report_rows:
            print(f"{file}\t{status}\t{details}", file=report)

    if not selected_alignments:
        raise NoGenesIncludedError(
            "No genes were included in the supermatrix. "
            f"See {report_file}; use --missing-taxa pad-gaps to retain genes "
            "with missing taxa."
        )

    partition_names = {}
    for file, _ in selected_alignments:
        partition_names.setdefault(tidy_name(file), []).append(file)
    duplicate_partitions = {
        name: files for name, files in partition_names.items() if len(files) > 1
    }
    if duplicate_partitions:
        details = "; ".join(
            f"{name}: {', '.join(files)}"
            for name, files in duplicate_partitions.items()
        )
        raise ValueError(
            "Duplicate partition names after filename normalization: " + details
        )

    partition_path = prefix + "_partition_finder.cfg"
    iqtree_path = prefix + '_part_iqtree.txt'
    #part_raxml = open('part_raxml.txt', 'w')

    with atomic_output(partition_path) as partition, atomic_output(iqtree_path) as part_iqtree:
      print(f"""## ALIGNMENT FILE ##
alignment = {os.path.basename(prefix)+".phy"};

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = mrbayes;

## MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]    
""", file=partition)
    
      print("#nexus\nbegin sets;", file=part_iqtree)
    
      for n, (file, alignment_length) in enumerate(selected_alignments):
        if n==0:
            start_pos = 1
            end_pos = alignment_length
        else:
            start_pos = end_pos + 1
            end_pos = end_pos + alignment_length
        
        print(f"{tidy_name(file)} = {start_pos}-{end_pos};", file=partition)
        #print(f"DNA, {tidy_name(file)} = {start_pos}-{end_pos}", file=part_raxml)
        print(f"    charset {tidy_name(file)} = {start_pos}-{end_pos};", file=part_iqtree)
      included_files = [file for file, _ in selected_alignments]
    # Keep the empty model prefix for the user's IQ-TREE -m MFP workflow.
      print(f"    charpartition my_genes = :{', :'.join([tidy_name(file) for file in included_files])};\nend;", file=part_iqtree)
      print(f"#partition my_genes = {len(included_files)} : {', '.join([tidy_name(file) for file in included_files])};\n#end;", file=partition)
      print("""
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;    
""", file=partition)

    supermatrix_length = sum(length for _, length in selected_alignments)
    spool.write_outputs(prefix, taxa_order, supermatrix_length, write_phylip)
    return None
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build a supermatrix and partition files from FASTA alignments.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  get_supergenes.py -i gene1.fasta gene2.fasta -p supermatrix
  get_supergenes.py -i gene1.fasta gene2.fasta -p supermatrix \\
      --missing-taxa pad-gaps""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='FASTA|LIST',
                          help='FASTA files and/or path-list files (one FASTA path per line).', nargs='+', required=True)
    required.add_argument('-p', '--prefix', metavar='PREFIX', required=True,
                          help='Prefix for supermatrix, partition, and report files.')
    optional.add_argument('--missing-taxa', choices=['skip-gene', 'pad-gaps'],
                          default='skip-gene',
                          help='How to handle a taxon missing from a gene (default: skip-gene).')
    optional.add_argument('--no-phylip', action='store_true',
                          help='Skip PHYLIP output when a downstream converter will generate it.')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0',
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    try:
        get_supergenes(infiles=expand_input_paths(args.input), prefix=args.prefix, missing_taxa=args.missing_taxa, input_lists=args.input, write_phylip=not args.no_phylip)
    except (OSError, ValueError) as error:
        print(f"get_supergenes: {error}", file=sys.stderr)
        sys.exit(3 if isinstance(error, NoGenesIncludedError) else 2)
    
