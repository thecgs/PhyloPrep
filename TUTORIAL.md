# PhyloPrep Containerized Example Tutorial

This directory contains two reproducible phylogenomics examples, both run with the [`phyloprep.sif`](https://zenodo.org/api/records/22949441/draft/files/phyloprep.sif/content) container image. After downloading and unpacking the example archive, keep the following directory layout:

```text
my_path/
├── phyloprep.sif
├── benchmark_mitochondrial_genes/
│   ├── data/
│   └── run.sh
└── benchmark_nuclear_genes/
    ├── extract_sequence_from_gff3.py
    ├── get_longest_transcript_gff3.py
    └── run.sh
```

> **Zenodo:** The container image and both complete examples will be deposited on Zenodo. `[Zenodo DOI: https://doi.org/10.5281/zenodo.22949440]`

## Prerequisites

Use Linux or another POSIX environment that supports Singularity/Apptainer.
You will need:

- SingularityCE or Apptainer;
- `bash`;
- for the nuclear-gene example, host installations of `wget`, `bgzip`
  (htslib), and `grep`;
- access to the EBI Ensembl FTP server for the first nuclear-genome download.

If using Apptainer, replace `singularity` with `apptainer` throughout this tutorial.

```bash
singularity exec phyloprep.sif phyloprep.py --help
```

If this command prints the help message, the image is ready to use. It includes
PhyloPrep, MAFFT, trimAl, OrthoFinder 2.5.5, and their required dependencies.

## Example 1: mitochondrial protein-coding genes

`benchmark_mitochondrial_genes/data/` contains CDS FASTA files for 13 vertebrate mitochondrial protein-coding genes, together with `map.tsv`. The first column of this mapping file is the original sequence ID and the second is the species name. PhyloPrep uses it to rename sequences by species after alignment and then construct species-level supermatrices.

Run the complete workflow:

```bash
cd benchmark_mitochondrial_genes
bash run.sh
```

The core command in `run.sh` is:

```bash
singularity exec ../phyloprep.sif phyloprep.py \
  -i data/*.fasta -m data/map.tsv \
  -j 13 -t 1 -ts trimal -st codon -s mafft -g 2
```

Here, `-st codon` selects codon-aware CDS processing; `-g 2` selects the vertebrate mitochondrial genetic code; `-s mafft` selects MAFFT; and `-ts trimal` enables trimAl. The workflow processes the 13 genes concurrently with `-j 13`, while allocating one thread per alignment with `-t 1`. Adjust these two settings for the available computational resources.

Key output files are written below `phyloprep-results/`:

- `matrices/codon/{raw,trimmed}/supermatrix.{fasta,phy,nex}`: concatenated CDS
  supermatrices;
- `matrices/codon/{raw,trimmed}/codon{1st,2nd,3rd}.supermatrix.*`: matrices
  for the three codon positions;
- `matrices/codon/{raw,trimmed}/fourfold.*`: four-fold degenerate sites;
- `matrices/prot/{raw,trimmed}/supermatrix.*`: translated protein
  supermatrices;
- `qc/` and `reports/`: quality-control and sequence-standardization records.

## Example 2: nuclear single-copy orthologs

This example uses zebrafish, western clawed frog, mouse, chicken, and human. Its workflow is: download reference genomes and GFF3 annotations → retain the longest transcript per gene → extract CDS and protein sequences → infer orthogroups with OrthoFinder → extract one-to-one ortholog CDS sequences → build alignments and supermatrices with PhyloPrep.

```bash
cd benchmark_nuclear_genes
bash run.sh
```

The script checks for each downloaded genome and GFF3 file by name and skips the download when it is already present. The first run downloads large genome files, and OrthoFinder plus the subsequent alignments may require substantial time. The script uses `mkdir -p cds pep`, so existing directories do not stop the workflow; rerun behavior for later stages depends on their outputs and the corresponding tool's resume behavior.

### Workflow stages and intermediate results

1. **Longest transcripts and sequence extraction:**
   `get_longest_transcript_gff3.py` selects the longest transcript from the non-mitochondrial annotations of each GFF3 file (records beginning with `MT` are excluded). `extract_sequence_from_gff3.py` uses the genome FASTA to extract CDS (`-t CDS`) and protein sequences (`-t prot`). Both scripts are provided with this example and originate from [QuickProt](https://github.com/thecgs/quickprot).

2. **Ortholog inference:** 

   OrthoFinder in the container analyses the five proteomes in `pep/` and writes its results to `OrthoFinder/`.

3. **One-to-one orthologs:** 

   `extract_orthofinder_orthogroups.py` combines the OrthoFinder `Orthogroups.tsv` file with `cds/` and writes `one-to-one_orthologs/`. Its `orthogroups.pathlist` is the PhyloPrep input list, while `map.tsv` maps transcript IDs to species names.

4. **Codon alignment and concatenation:** 

   The final `phyloprep.py` command uses the standard genetic code (`-g 1`) and writes raw and trimmed CDS, protein, codon-position, and four-fold-degenerate matrices, together with QC reports, under `phyloprep-results/`.

The final command in this example is:

```bash
singularity exec ../phyloprep.sif phyloprep.py \
  -i one-to-one_orthologs/orthogroups.pathlist \
  -m one-to-one_orthologs/map.tsv \
  -j 20 -t 1 -ts trimal -st codon -s mafft -g 1
```

To rerun one stage from scratch, use a separate working copy and remove only that stage's output directory before rerunning its command. Do not casually remove input data from the original example archive. PhyloPrep normally reuses completed results whose inputs have not changed.

## Using your own data

For mitochondrial or other preassembled CDS datasets, follow Example 1: use one FASTA file per gene and provide a two-column `map.tsv` (`sequence_ID<TAB>species_ID`). For nuclear data, follow Example 2: extract CDS and proteins from each species' genome and annotation, use OrthoFinder to select single-copy orthologs, and supply the resulting `orthogroups.pathlist` and `map.tsv` to PhyloPrep. See the repository-level `README.md` for complete details about input requirements, QC thresholds, and output files.
