# Reproducibility tutorial

This directory contains two reproducible phylogenomics examples, both run with the `phyloprep.sif` container image. 

| Example | Input | Purpose |
| --- | --- | --- |
| Mitochondrial CDS | `00.data/mitochondrial_data/` | Run the complete codon workflow from per-gene FASTA files. |
| Nuclear single-copy orthologs | `00.data/nuclear_data/` | Infer orthogroups, extract one-to-one CDS orthologs, and construct species matrices. |

Precomputed output directories are included as reference results. To avoid
overwriting them, use the review-specific output names shown below.

## Requirements

Run all commands from this `example/` directory in a POSIX shell. Either of
the following environments is sufficient:

- A SingularityCE or Apptainer image named `phyloprep-v1.0.0.sif` in this
  directory. Replace `singularity` with `apptainer` if needed.
- A local installation with Python, Biopython, MAFFT, trimAl, and (for the
  nuclear example) OrthoFinder available on `PATH`.

Confirm the command-line interface before running an analysis:

```bash
singularity exec ./phyloprep-v1.0.0.sif phyloprep.py --help
```

For a local installation, replace the prefix above with `python phyloprep.py`.

## 1. Mitochondrial CDS workflow

The directory `00.data/mitochondrial_data/` contains one coding-sequence
FASTA file per locus and `map.tsv`. The mapping table has two tab-separated
columns:

```text
input_sequence_id    species_id
```

The first column must exactly match FASTA IDs. Mapping occurs after alignment,
so raw inputs and the audit trail retain the original IDs; final matrices use
the mapped species IDs.

Run the example into a new directory:

```bash
singularity exec ./phyloprep-v1.0.0.sif phyloprep.py \
  -i ./00.data/mitochondrial_data/*.fasta \
  -m ./00.data/mitochondrial_data/map.tsv \
  -st codon -g 2 -s mafft -ts trimAlnSeq \
  -j 4 -t 1 \
  -o review_mitochondrial_results
```

`-g 2` selects the vertebrate mitochondrial genetic code. `-j` is the number
of loci processed concurrently, and `-t` is the number of alignment threads
per locus; choose values appropriate for the review machine. The supplied
`run_example_mitochondrial.sh` contains an equivalent command.

Expected primary outputs include:

```text
review_mitochondrial_results/
├── matrices/codon/raw/supermatrix.{fasta,phy,nex}
├── matrices/codon/trimmed/supermatrix.{fasta,phy,nex}
├── matrices/codon/trimmed/codon{1st,2nd,3rd}.supermatrix.fasta
├── matrices/codon/trimmed/fourfold.fasta
├── matrices/prot/{raw,trimmed}/supermatrix.fasta
├── qc/codon/{raw,trimmed}/alignment.{pass,fail,details}.tsv
└── reports/sequence_changes.tsv
```

Raw and trimmed matrices are evaluated independently. `sequence_changes.tsv`
records sequence normalization and validation events with the original input
file and original FASTA ID; when `-m` is supplied, it also includes
`mapped_taxa_id`.

## 2. Nuclear single-copy ortholog workflow

This example starts from per-taxon CDS and protein FASTA files under
`00.data/nuclear_data/`. It has three stages. Each stage may be inspected or
run independently.

### 2.1 Infer orthogroups from proteins

```bash
singularity exec ./phyloprep-v1.0.0.sif orthofinder \
  -f 00.data/nuclear_data/pep/ -M msa -S diamond -T fasttree \
  -o review_nuclear_orthofinder
```

Use the `Orthogroups.tsv` produced within the resulting `Results_*/`
directory in the next step. The wildcard avoids assuming a tool-generated
timestamped directory name.

### 2.2 Extract one-to-one CDS orthologs

```bash
singularity exec ./phyloprep-v1.0.0.sif extract_orthofinder_orthogroups.py \
  -i 00.data/nuclear_data/cds/ \
  -og review_nuclear_orthofinder/Results_*/Orthogroups/Orthogroups.tsv \
  -o review_nuclear_orthologs
```

This writes one FASTA file per retained orthogroup, plus:

- `review_nuclear_orthologs/orthogroups.pathlist`: absolute paths to the
  exported FASTA files;
- `review_nuclear_orthologs/map.tsv`: FASTA ID to species ID mapping.

The paths in the generated list are valid while the directory remains in its
current location. If it is moved, rerun this extraction step.

### 2.3 Align, QC, concatenate, and export matrices

```bash
singularity exec ./phyloprep-v1.0.0.sif phyloprep.py \
  -i review_nuclear_orthologs/orthogroups.pathlist \
  -m review_nuclear_orthologs/map.tsv \
  -st codon -g 1 -s mafft -ts trimal \
  -j 4 -t 1 \
  -o review_nuclear_results
```

`-g 1` selects the standard genetic code. In codon mode with `-ts trimal`,
PhyloPrep trims the paired protein alignment and transfers retained columns as
complete codon triplets. It does not trim nucleotide positions independently.

## 3. Review checks

After either workflow, the following checks give a compact assessment of the
result and its audit trail:

```bash
find review_mitochondrial_results/matrices -name '*.nex' -type f | sort
head review_mitochondrial_results/qc/codon/trimmed/alignment.pass.tsv
head review_mitochondrial_results/reports/sequence_changes.tsv
```

For codon matrices, sequence lengths in `supermatrix.fasta` are divisible by
three. The position-specific matrices and four-fold matrix are derived only
after a supermatrix is successfully produced. A normal codon workflow stops if
no usable four-fold sites remain; pseudogene mode instead writes a
`fourfold.skipped.txt` marker and retains its other matrices.

## 4. Optional modes

Protein inputs can be processed without trimming:

```bash
singularity exec ./phyloprep-v1.0.0.sif phyloprep.py \
  -i proteins.pathlist -st prot --notrim -o review_protein_results
```

For frameshift- or stop-tolerant coding sequences, use the pseudogene mode
with MACSE available in the image:

```bash
singularity exec ./phyloprep-v1.0.0.sif phyloprep.py \
  -i pseudogene.pathlist -st pseudogene \
  --macse-jar macse_v2.07.jar -j 2 -o review_pseudogene_results
```

## 5. Using your own data

For mitochondrial or other preassembled CDS datasets, follow Example 1: use one FASTA file per gene and provide a two-column `map.tsv` (`sequence_ID<TAB>species_ID`). For nuclear data, follow Example 2: extract CDS and proteins from each species' genome and annotation, use OrthoFinder to select single-copy orthologs, and supply the resulting `orthogroups.pathlist` and `map.tsv` to PhyloPrep. See the repository-level `README.md` for complete details about input requirements, QC thresholds, and output files.

## 6. Reproducibility notes

- Inputs are never modified in place.
- Output publication is atomic: incomplete conversions do not replace an
  existing final matrix.
- Completed per-locus work is reused when inputs and relevant settings are
  unchanged. Add `--no-resume` to recompute all loci.
