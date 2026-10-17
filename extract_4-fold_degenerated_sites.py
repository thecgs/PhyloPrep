#!/usr/bin/env python
# coding: utf-8

import argparse
import tempfile
from pathlib import Path
from msap_matrix import AlignmentSpool, CHUNK_SIZE
from Bio.Data import CodonTable
from collections import defaultdict
from msap_io import atomic_output, check_output_path

def get_four_fold_codons(genetic_code=1):
    """
    Return a set containing all four-fold degenerate codons.
    """
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]

    # amino acid -> {first_two_bases: count}
    aa_prefix_count = defaultdict(lambda: defaultdict(int))

    for codon, aa in codon_table.forward_table.items():
        aa_prefix_count[aa][codon[:2]] += 1

    four_fold_codons = set()

    for codon, aa in codon_table.forward_table.items():
        if aa_prefix_count[aa][codon[:2]] == 4:
            four_fold_codons.add(codon)

    return four_fold_codons
    
def extract_four_fold_site(infile, outfile, genetic_code=1):
    """Intersect qualifying codons in bounded blocks, then stream third bases."""
    check_output_path(outfile, [infile])
    directory = Path(outfile).absolute().parent
    four_fold_codons = get_four_fold_codons(genetic_code)
    with AlignmentSpool(infile, codon=True, directory=directory,
                        transform=lambda _, chunk, offset: chunk.upper().replace('U', 'T')) as matrix, \
            tempfile.TemporaryFile(mode='w+b', dir=directory) as mask:
        positions = [record[1] for record in matrix.records]
        for start in range(0, matrix.length, CHUNK_SIZE):
            size = min(CHUNK_SIZE, matrix.length - start)
            reference = None
            keep = bytearray(size // 3)
            for index in range(len(matrix.records)):
                matrix.file.seek(positions[index])
                chunk = matrix.file.read(size)
                positions[index] = matrix.file.tell()
                if reference is None:
                    reference = chunk
                    for site, pos in enumerate(range(0, size, 3)):
                        keep[site] = chunk[pos:pos + 3] in four_fold_codons
                else:
                    for site, pos in enumerate(range(0, size, 3)):
                        if keep[site]:
                            codon = chunk[pos:pos + 3]
                            if codon not in four_fold_codons or codon[:2] != reference[pos:pos + 2]:
                                keep[site] = 0
            mask.write(keep)
        with atomic_output(outfile, inputs=[infile]) as out:
            for index, (identifier, _, _) in enumerate(matrix.records):
                out.write(f'>{identifier}\n')
                mask.seek(0)
                for chunk in matrix.chunks(index):
                    keep = mask.read(len(chunk) // 3)
                    out.write(''.join(chunk[site * 3 + 2] for site, retain in enumerate(keep) if retain))
                out.write('\n')
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract third positions from four-fold degenerate codon sites in a codon FASTA alignment.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  extract_4-fold_degenerated_sites.py -i codon.aln.fasta -o fourfold.fasta -g 1

The genetic-code table uses NCBI IDs (default: 1).""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='CODON_FASTA', help='Input codon alignment in FASTA format.', required=True)
    required.add_argument('-o', '--output', metavar='FASTA', help='Output four-fold-site alignment in FASTA format.', required=True)
    optional.add_argument('-g', '--genetic_code', '--genetic-code', metavar='TABLE', default=1, type=int,
                          help='NCBI genetic-code table (default: 1).')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0', help="Show program's version number and exit.")
    args = parser.parse_args()
    extract_four_fold_site(infile=args.input, outfile=args.output, genetic_code=args.genetic_code)
