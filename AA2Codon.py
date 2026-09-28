#!/usr/bin/env python
# coding: utf-8

import os, sys
import textwrap
import argparse
from Bio import SeqIO, Seq
from Bio.Data import CodonTable
from msap_io import atomic_output, normalize_dna

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Back-translate a protein FASTA alignment to a codon alignment.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  AA2Codon.py -c CDS.fasta -p protein.aln.fasta -o codon.aln.fasta
  AA2Codon.py -c CDS.fasta -p protein.aln.fasta -g 1""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-c', '--cds', metavar='CDS_FASTA', required=True,
                          help='Input CDS FASTA file; RNA U is normalized to DNA T.')
    required.add_argument('-p', '--protaln', '--protein-alignment', metavar='PROTEIN_FASTA', required=True,
                          help='Input protein alignment in FASTA format.')
    optional.add_argument('-o', '--out', '--output', metavar='CODON_FASTA', default=None,
                          help='Output codon alignment in FASTA format (default: stdout).')
    optional.add_argument('-g', '--genetic_code', '--genetic-code', metavar='TABLE', default=None,
                          type=int, help='NCBI genetic-code table used for validation; recognized initial start codons validate as M (default: disabled).')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0',
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    cdsfile=args.cds
    protalnfile = args.protaln
    outfile = args.out
    genetic_code = args.genetic_code
    
    mapping = {}
    for record in SeqIO.parse(cdsfile, 'fasta'):
        if record.id in mapping:
            raise ValueError(f"Duplicate CDS ID: {record.id}")
        if len(record.seq) % 3 != 0:
            raise ValueError(
                f"{record.id}: CDS length {len(record.seq)} is not divisible by 3"
            )
        mapping[record.id] = textwrap.wrap(normalize_dna(record.seq), 3)

    if not mapping:
        raise ValueError("No CDS sequences found.")

    protein_records = list(SeqIO.parse(protalnfile, 'fasta'))
    if not protein_records:
        raise ValueError("No protein-alignment sequences found.")
    protein_ids = set()
    table = CodonTable.unambiguous_dna_by_id[genetic_code] if genetic_code is not None else None
    for record in protein_records:
        if not record.seq:
            raise ValueError(f"{record.id}: protein alignment sequence is empty")
        if len(record.seq) != len(protein_records[0].seq):
            raise ValueError(f"{record.id}: protein alignment length mismatch")
        if record.id in protein_ids:
            raise ValueError(f"Duplicate protein-alignment ID: {record.id}")
        protein_ids.add(record.id)

    cds_ids = set(mapping)
    if cds_ids != protein_ids:
        missing_from_protein = sorted(cds_ids - protein_ids)
        missing_from_cds = sorted(protein_ids - cds_ids)
        details = []
        if missing_from_protein:
            details.append("missing from protein alignment: " + ", ".join(missing_from_protein))
        if missing_from_cds:
            details.append("missing from CDS: " + ", ".join(missing_from_cds))
        raise ValueError("CDS/protein alignment ID mismatch: " + "; ".join(details))

    for record in protein_records:
        protein_sequence = str(record.seq).upper()
        residue_count = sum(amino_acid != "-" for amino_acid in protein_sequence)
        codons = mapping[record.id]
        codon_count = len(codons)
        if residue_count != codon_count:
            raise ValueError(
                f"{record.id}: protein alignment has {residue_count} non-gap residues, "
                f"but CDS has {codon_count} codons"
            )
        if genetic_code is not None:
            codon_index = 0
            for amino_acid in protein_sequence:
                if amino_acid == "-":
                    continue
                codon = codons[codon_index]
                translated = str(Seq.Seq(codon).translate(table=genetic_code))
                # Match MSAP.translate_seq(): a codon recognized as an
                # initiator by the selected table is represented as M only
                # at the first translated position.
                expected = "M" if codon_index == 0 and codon in table.start_codons else translated
                if expected != amino_acid:
                    raise ValueError(
                        f"{record.id}: codon index {codon_index} ({codon}) translates to "
                        f"{expected}, but the protein alignment contains {amino_acid}"
                    )
                codon_index += 1

    with atomic_output(outfile, inputs=[cdsfile, protalnfile]) as out:
        for record in protein_records:
            print(">"+record.id, file=out)
            s = []
            idx = 0
            for i in record.seq.upper():
                if i != '-':
                    codon = mapping[record.id][idx]
                    s.append(codon)
                    idx += 1
                else:
                    s.append('---')
            print("".join(s), file=out)
