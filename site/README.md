# PhyloPrep Browser Edition

This directory is a buildless static site designed for GitHub Pages. It uses
Pyodide to run the repository's original selected Python scripts locally in a
visitor's browser; no FASTA content is posted to a server by this page.

## Deploy with GitHub Pages

1. Push this repository to a GitHub repository you control.
2. Open **Settings → Pages** and select **GitHub Actions** as the source.
3. Push to `main`, or use **Actions → Deploy PhyloPrep Browser Edition → Run workflow**.

The included workflow generates `site/engine/` from the original source files,
then uploads `site/`. For a local preview, run `node site/build-engine.mjs`
from the repository root before starting an HTTP server for `site/`.

The browser edition runs the listed standalone preparation modules through the
same original scripts used by the CLI, including quality control,
back-translation, four-fold-site tools, 4DTV, supermatrix construction,
filtering, renaming, trimming, and conversion.

The browser also exposes the original single-locus `MSAP.py` workflow for
`codon`, `prot`, and `nucl` input. `MSAP.py` itself is not modified: the browser
runs it in Pyodide until it requests MAFFT or MUSCLE, delegates that command to
Biowasm MUSCLE 5.1.0 or optional MAFFT E-INS-i 7.520 (`tbfast` + `dvtditr`)/WebAssembly, and then resumes the same
original workflow with the result. Back-translation and the default
`trimAlnSeq.py` trimming remain the original bundled Python scripts. The
browser implementation uses one alignment worker. The Biowasm MAFFT package
does not provide an upstream `mafft.js` wrapper; the browser runs the compiled
`tbfast` + `dvtditr` programs directly instead. PRANK, ClustalW2, trimAl, MACSE, and pseudogene mode
remain CLI-only.

FASTA content remains local to the browser. Pyodide and Biowasm assets are
downloaded from their public CDNs on first use, but sequence data are not sent
to PhyloPrep or a PhyloPrep server.


## Upstream browser workflows

The browser catalog also exposes `MSAP_batch.py` and `phyloprep.py`.
Because browsers do not provide the Linux process-management primitives used by
the CLI orchestrators (process groups, `fcntl` locks, and normal child
processes), the web edition uses a browser-safe compatibility layer rather than
modifying those core scripts.

`MSAP_batch.py` browser mode processes loci sequentially with the existing
MSAP/Pyodide + Biowasm MAFFT/MUSCLE workflow, then creates batch-style path lists and
a merged sequence-change report. Those validated path lists are also the browser finalizer's only source of downstream alignment membership, matching the CLI contract. `phyloprep.py` browser mode reuses those
alignments and invokes the repository's provenance-aware paired codon/protein QC,
taxon-renaming, supermatrix, codon-position, four-fold and FASTA/PHYLIP/NEXUS
modules inside Pyodide.

Current browser constraints are explicit in the parameter UI: MAFFT E-INS-i
7.520 is the default; MUSCLE 5.1.0 remains available for explicit selection. The browser uses one
alignment worker, one locus at a time, and
`trimAlnSeq.py` trimming. Pseudogene/MACSE, PRANK, ClustalW2 and trimAl remain
CLI-only.

### Stop provenance and paired codon/protein QC

In browser `phyloprep.py`, codon and protein alignments from the same CDS locus are treated as a pair. The browser uses the merged `reports/sequence_changes.tsv` automatically when running QC. Internal stops that MSAP already masked are therefore still visible to QC as provenance events in both the codon and protein representations. They fail the pair by default and are allowed only when `--allow-internal-stop` is selected; the masked sequence state is never restored. Terminal stops normally removed from complete CDS inputs are reported but do not fail QC. The complete workflow has no `--allow-terminal-stop` switch: unexpected directly observed stops in normal codon/protein outputs fail QC. Standalone `alignment_qc.py` retains this option for externally supplied alignments, permitting observed terminal stops without modifying sequences.

The paired workflow also keeps codon and protein site structure synchronized. Browser codon trimming uses the same paired site mask as the CLI: gap, codon-N, and protein-X thresholds are evaluated for the same amino-acid/codon site, and rejected sites are removed from both representations together. For the additional paired datasets, both sides must pass non-PIS representation-specific QC. PIS is evaluated separately for codon and protein, and the browser exposes a **Paired PIS policy** selector with `both`, `either`, `codon`, and `protein`. Optional **Minimum codon PIS** and **Minimum protein PIS** fields let the two representations use different thresholds; blank fields fall back to the general minimum-PIS value. The additional matrices in `matrices/paired/<codon|prot>/<raw|trimmed>/` use the same retained locus set. Independent matrices in `matrices/<codon|prot>/<raw|trimmed>/` use each side's own QC and PIS threshold; they include all independently passing loci. Representation-specific pass/fail reports are preserved, while final retained manifests are written as `qc/paired/<raw|trimmed>/codon.retained.tsv` and `prot.retained.tsv`. Pair decisions, both PIS values, both thresholds, selected policy, and final status are written under `qc/paired/<raw|trimmed>/loci.details.tsv`.

The standalone browser `alignment_qc.py` tool accepts an optional `sequence_changes.tsv` input. Without it, directly visible stops are detected from codons or protein `*`; with it, masked internal-stop provenance is restored for QC. Terminal and internal observed-stop controls are therefore available for codon and standalone protein QC.


## Orthology tools in the browser

The web catalog includes the OrthoFinder orthogroup extractor plus adapters for
OMA, OrthoMCL, Proteinortho and SonicParanoid. The OrthoFinder extractor accepts
multiple per-taxon FASTA files in the browser and reconstructs the CLI's
sequence-directory layout internally before running the original
`extract_orthofinder_orthogroups.py`.

`orthogroup_table.py` remains bundled only as an internal shared library for
the orthology adapters and is not shown as a standalone web tool.

The browser also exposes `split_codon_seqence_alignment.py` to write first-,
second- and third-codon-position FASTA matrices. The single-locus `MSAP.py`
entry is intentionally hidden from the web catalog because `MSAP_batch.py`
covers the same browser-facing workflow, while `MSAP.py` remains bundled
internally as the implementation used by the batch adapter.


## Workspace safeguards and reproducibility

The browser workspace now performs lightweight input preflight checks before a
run, groups and searches the tool catalog, hides parameters that do not apply to
the selected sequence type, and supports direct handoff between common
orthology-preparation steps.

Each successful run adds a `phyloprep-run.json` manifest containing the
selected module, parameters, equivalent CLI command, runtime provenance, and
SHA256 hashes for input and output content. Pyodide run directories are removed
after each worker request so repeated batch runs do not accumulate virtual
filesystem workspaces.

MAFFT E-INS-i and MUSCLE are executed in a dedicated alignment Web Worker. Cancelling a
run terminates that worker and therefore interrupts the active WebAssembly
alignment rather than only cancelling the surrounding Pyodide step.

The Pages workflow runs syntax and browser-bundle smoke checks before building
and deploying the static site.


## Multi-file browser input

For `phyloprep.py`, `MSAP_batch.py`, and the existing multi-alignment tools,
selecting more than one FASTA file automatically switches the browser to a
path-list input model. The equivalent CLI command therefore uses one generated
`.pathlist` instead of expanding every FASTA path on the command line.
`MSAP_batch.py` retains the generated list as `path-lists/inputs.pathlist`;
the browser-adapted `phyloprep.py` result retains `inputs.pathlist`, matching
the top-level pipeline's CLI organization.


## Repeated browser alignments

The Biowasm MUSCLE package is configured with Aioli `reinit: true`. Batch workflows reuse
the alignment worker across loci, but the MUSCLE WebAssembly module itself is
reinitialized after each invocation so argument-parser and other module-global
state cannot leak from one FASTA alignment into the next.


MAFFT E-INS-i uses a stricter isolation model than MUSCLE. Each locus creates a
fresh Aioli instance containing `tbfast` and `dvtditr`; that instance is
closed after the two-stage alignment completes. This preserves the shared
`/shared/data/pre` file within a locus while preventing MAFFT module-global
state from leaking into the next FASTA in a batch.

Independent retained lists are written as `qc/<codon|prot>/<raw|trimmed>/independent.retained.tsv`. QC runs once for each representation/variant. Empty datasets produce `skipped.txt`, while other passing datasets are retained. Identical independent/paired selections reuse matrix products.

### Local pathlists and large selections

Select a `.pathlist` and its input folder, or select the referenced FASTAs alongside it. Browsers require explicit directory authorization; a pathlist cannot grant filesystem access. The resolver checks referenced paths within authorized inputs, then same-directory basenames. Missing and ambiguous entries report the list and filename. Only referenced files enter a list-driven run. RUN reads at most four files concurrently; preflight samples at most 100 FASTAs and displays the resolved count. CLI utilities use the shared resolver with the same basename fallback for moved folders.

Browser batch progress shows completed/total loci with elapsed time and estimated remaining time. Python QC and matrix construction stream stage updates through the existing worker bridge. The interface refreshes elapsed time while a stage is active and stops timers on result, error or cancellation. Stage counts are not compute-time percentages; progress requires no additional libraries.

The web catalog and selected-tool details load each script version from `engine/manifest.json`. The build reads the scripts’ own `--version` declarations and records their source SHA256 hashes; version strings are not duplicated in the interface. `phyloprep-run.json` records the version and hashes of the engine bundle actually used. Biowasm package versions come from the alignment-worker configuration and are labeled as package versions; a program’s own banner can report a different version. A missing declaration or failed manifest request is shown explicitly.

### Self-hosted 4 GiB aligners

Pages builds now compile the three alignment programs with a 4 GiB memory-growth ceiling and load them from `aligners/` on the same site. Build them before running `build-engine.mjs` for local previews. See [build instructions and limits](../docs/browser-aligner-builds.md). The ceiling is not a memory reservation or a guarantee of successful alignment; the browser may provide less. Aligner failures stop immediately, and diagnostics are specific to MAFFT or MUSCLE.

The OrthoFinder extractor offers an optional **Keep copy with greatest non-ambiguous length** checkbox in `all-copies` mode (`--longest-copy`). It selects separately within each OG and taxon, using non-ambiguous character counts and gene-ID order to break effective-length ties. Set **Sequence type for copy selection** using `-st/--seqtype/--seq-type`: `codon` (default), `nucl` and `pseudogene` all count A/C/G/T/U; `prot` counts the 20 standard amino acids, including N. These match the CLI choices and default. Copy selection does not translate, group codons or invoke MACSE. Ambiguity symbols, gaps and stops do not count; raw length does not break ties, and exported sequences remain unchanged. The sequence type is not inferred automatically. Switching to single-copy omits the option from both execution and the CLI preview. `copy_selection.tsv` audits candidates and choices with original `length`, `effective_length` and `sequence_type`. Longest-copy selection does not infer orthology or identify paralogs.

Extractor policy notices are shown before execution and in successful results: a relaxed single-copy taxon threshold prompts `--missing-taxa pad-gaps` for downstream retention of incomplete loci, and longest-copy selection warns about different paralogs across species. The selected-taxon count reflects keep/drop filters. Warnings are also recorded in `phyloprep-run.json`; no downstream option is changed automatically.

The OG extractor scans each selected taxon FASTA once and retains Python sequence contents only for genes in retained OGs (all eligible candidates when longest-copy selection is enabled). Source files still occupy browser staging storage; this optimization does not remove browser memory limits.

**Minimum occupied taxa** (`--min-taxa`) is an integer taxon count, not a percentage or copy count: `single-copy` counts taxa with exactly one listed gene; `all-copies` counts taxa with at least one listed gene, once per taxon, including with `--longest-copy`. It measures OG membership, not sequence quality. The default is all taxa remaining after keep/drop filtering.

In `all-copies` mode, **Maximum copies per taxon** (`--max-copies-per-taxon`) defaults visibly to **inf (unlimited)**. The field accepts `inf` or a positive integer. A positive integer excludes an entire OG if any selected taxon exceeds the inclusive limit, after keep/drop filtering and before longest-copy selection. It neither truncates copies nor omits taxa. Switching to single-copy hides the field and omits it from execution and the command preview. This limit does not implement KinFin's target single-copy fraction.

### Bounded concatenation and diagnostics

Concatenation now uses two passes, temporary taxon spools and chunked output rather than retaining every locus plus complete joined matrices in Python. Spools are released taxon by taxon during final output. Browser MEMFS still consumes memory for temporary and final files. Captured stdout/stderr each retain only a tail of up to 65,536 characters (and at most 1,024 fragments), with an explicit truncation marker. Browser entry points return a pre-serialized short diagnostic on MemoryError, avoiding traceback/log collection and another JSON allocation while memory is exhausted. These changes reduce duplication; they do not guarantee that an arbitrarily large dataset fits in Pyodide.

Run progress is mirrored in the main workspace log as well as the stage label. Input preparation reports resolved path counts, read file counts and elapsed time before worker execution; output saving reports completed file counts and the current file. Pathlist folder roots are indexed once instead of rescanned for each entry.
