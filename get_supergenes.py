#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
from Bio import SeqIO
from msap_io import atomic_output

def tidy_name(file):
    file = os.path.basename(file).replace('.', '_').replace('-', '_')
    file = file.replace("_mafft_nucl_trimal_aln", "")
    file = file.replace("_mafft_codon_trimal_aln", "")
    file = file.replace("_mafft_prot_trimal_aln", "")
    file = file.replace("_mafft_nucl_aln", "")
    file = file.replace("_mafft_codon_aln", "")
    file = file.replace("_mafft_prot_aln", "")

    file = file.replace("_muscle_nucl_trimal_aln", "")
    file = file.replace("_muscle_codon_trimal_aln", "")
    file = file.replace("_muscle_prot_trimal_aln", "")
    file = file.replace("_muscle_nucl_aln", "")
    file = file.replace("_muscle_codon_aln", "")
    file = file.replace("_muscle_prot_aln", "")

    file = file.replace("_clustalw2_nucl_trimal_aln", "")
    file = file.replace("_clustalw2_codon_trimal_aln", "")
    file = file.replace("_clustalw2_prot_trimal_aln", "")
    file = file.replace("_clustalw2_nucl_aln", "")
    file = file.replace("_clustalw2_codon_aln", "")
    file = file.replace("_clustalw2_prot_aln", "")

    file = file.replace("_prank_nucl_trimal_aln", "")
    file = file.replace("_prank_codon_trimal_aln", "")
    file = file.replace("_prank_prot_trimal_aln", "")
    file = file.replace("_prank_nucl_aln", "")
    file = file.replace("_prank_codon_aln", "")
    file = file.replace("_prank_prot_aln", "")

    file = file.replace("_mafft_codon_trimal_fasta", "")
    file = file.replace("_muscle_codon_trimal_fasta", "")
    file = file.replace("_clustalw2_codon_trimal_fasta", "")
    file = file.replace("_prank_codon_trimal_fasta", "")

    return file


def read_alignment(file):
    """Read one FASTA alignment and validate its IDs and sequence lengths."""
    records = list(SeqIO.parse(file, "fasta"))
    if not records:
        return None, "empty alignment"

    sequences = {}
    expected_length = len(records[0].seq)
    for record in records:
        if record.id in sequences:
            raise ValueError(f"{file}: duplicate sequence ID: {record.id}")
        if len(record.seq) != expected_length:
            raise ValueError(
                f"{file}: {record.id} has length {len(record.seq)}, "
                f"expected {expected_length}"
            )
        sequences[record.id] = str(record.seq)
    if expected_length == 0:
        return None, "zero-column alignment"
    return sequences, None


def get_supergenes(infiles, prefix, missing_taxa="skip-gene"):
    """Create a supermatrix while handling taxa missing from individual genes."""
    alignments = []
    taxa_order = []
    taxa_seen = set()
    report_rows = []

    for file in infiles:
        if os.path.getsize(file) == 0:
            report_rows.append((file, "skipped", "empty file"))
            continue

        sequences, error = read_alignment(file)
        if error:
            report_rows.append((file, "skipped", error))
            continue

        for taxon in sequences:
            if taxon not in taxa_seen:
                taxa_order.append(taxon)
                taxa_seen.add(taxon)
        alignments.append((file, sequences))

    selected_alignments = []
    supergenes = {taxon: [] for taxon in taxa_order}

    for file, sequences in alignments:
        missing = [taxon for taxon in taxa_order if taxon not in sequences]
        if missing and missing_taxa == "skip-gene":
            reason = "missing taxa: " + ", ".join(missing)
            report_rows.append((file, "skipped", reason))
            print(f"Skipping {file}: {reason}", file=sys.stderr)
            continue

        alignment_length = len(next(iter(sequences.values())))
        if missing:
            report_rows.append((file, "padded", "missing taxa: " + ", ".join(missing)))
        else:
            report_rows.append((file, "included", ""))

        selected_alignments.append((file, alignment_length))
        for taxon in taxa_order:
            supergenes[taxon].append(sequences.get(taxon, "-" * alignment_length))

    report_file = prefix + ".report.tsv"
    with atomic_output(report_file) as report:
        print("gene\tstatus\tdetails", file=report)
        for file, status, details in report_rows:
            print(f"{file}\t{status}\t{details}", file=report)

    if not selected_alignments:
        raise ValueError(
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

    #part_raxml.close()
    supermatrix = {
        species: "".join(supergenes[species]) for species in taxa_order
    }
    supermatrix_length = len(supermatrix[taxa_order[0]])

    with atomic_output(prefix+".fasta") as out:
        for species in taxa_order:
            print(f'>{species}\n{supermatrix[species]}', file=out)

    with atomic_output(prefix+".phy") as out:
        print(f"{len(taxa_order)} {supermatrix_length}", file=out)
        for species in taxa_order:
            print(f"{species} {supermatrix[species]}", file=out)
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
    required.add_argument('-i', '--input', metavar='FASTA',
                          help='Input FASTA alignments.', nargs='+', required=True)
    required.add_argument('-p', '--prefix', metavar='PREFIX', required=True,
                          help='Prefix for supermatrix, partition, and report files.')
    optional.add_argument('--missing-taxa', choices=['skip-gene', 'pad-gaps'],
                          default='skip-gene',
                          help='How to handle a taxon missing from a gene (default: skip-gene).')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='get_supergenes v1.00',
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    get_supergenes(infiles=args.input, prefix=args.prefix, missing_taxa=args.missing_taxa)
    
