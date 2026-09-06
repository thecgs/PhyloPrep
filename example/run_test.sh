#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# Only the original 13 genes belong in the input list, including on reruns.
genes=(ATPase6 ATPase8 COX1 COX2 COX3 Cytb ND1 ND2 ND3 ND4 ND4L ND5 ND6)
inputs=()
for gene in "${genes[@]}"; do
    inputs+=("${PWD}/${gene}.fasta")
done

mkdir -p run-results
cd run-results
python ../../MSAP_batch.py -i "${inputs[@]}" -g 2 -G 0.2 -N 0.2 -X 0.2 -st codon -j 13 -t 1 -o alignments
python ../../get_supergenes.py -i alignments/*.mafft.codon.aln -p supermatrix_mafft_aln
python ../../extract_4-fold_degenerated_sites.py -i supermatrix_mafft_aln.fasta -o four_fold_site.fasta -g 2
python ../../split_codon_seqence_alignment.py -i supermatrix_mafft_aln.fasta
