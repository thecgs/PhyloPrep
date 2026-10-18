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
a merged sequence-change report. In the full `phyloprep.py` browser workflow, the
batch adapter validates the input set once and forwards alignments plus individual
audit reports; it skips standalone summaries and unused preprocessing intermediates.
The Python finalizer generates the path lists and merged report once. Its path
lists remain the only source of downstream alignment membership, matching the CLI
contract. It invokes the repository's provenance-aware paired codon/protein QC,
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

### Fingerprint reuse during processing

Saved-input SHA-256 verification runs while copying chunks into the analysis workspace, rather than reading the same Blob once to hash and again to stage it. Checksum mismatches reject and remove the incomplete staged input. The verified source digest is reused during renaming. Mapping fingerprints are calculated while parsing the mapping; renamed-output fingerprints are calculated while writing FASTA and reused for checkpoint/export metadata. The worker still checks output bytes during chunk transfer, which adds no separate disk-reading pass and detects incomplete temporary-disk reads. Pipeline alignments already confirmed saved and verified during staging are not retransferred after QC; temporary copies are released only after the main thread confirms the saved output is still registered.

Unchanged output fingerprints are reused within the current branch using file size, modification/change times and file identity. Branch caches are cleared after downstream export, and separate runs validate disk content again. The browser no longer writes hidden CLI rename-state files that browser recovery never consumes; its saved stage catalogs and run-manifest recovery remain the browser checkpoint mechanisms. CLI rename-state files and source/output provenance digests are retained.

### Reused matrix preparation and rename pathlists

Each final FASTA is prepared and normalized once for both PHYLIP and NEXUS export. The writers share the bounded disk spool, IDs, alignment length and normalization-event log; fourfold output is validated in that preparation instead of being parsed separately first. Paired-QC locus counts and equality flags are retained for matrix selection rather than rereading retained lists.

Each rename directory also exports one current-run pathlist (`all.TYPE.aln.fasta.pathlist` for raw, `all.TYPE.trimmed.aln.fasta.pathlist` for trimmed), with filenames relative to its own folder. The lists include all renamed inputs to QC, not just QC-passing loci, and work with the CLI after the complete folder is saved. Complete rename-checkpoint hits regenerate the same lists from the current locus collection.

### Bounded concatenation and diagnostics

Concatenation uses two sequence-reading passes, temporary taxon spools and chunked output rather than retaining every locus plus complete joined matrices in Python. The concatenator retains duplicate-ID and equal-length checks but does not repeatedly hash inputs; it assumes inputs remain unchanged during processing. Checkpoint, saved-input read-back and output-provenance checksums are unchanged. Spools are released taxon by taxon during final output. The temporary disk workspace described below avoids retaining working-file contents in MEMFS when supported; the fallback still uses MEMFS. Captured stdout/stderr each retain only a tail of up to 65,536 characters (and at most 1,024 fragments), with an explicit truncation marker. Browser entry points return a pre-serialized short diagnostic on MemoryError, avoiding traceback/log collection and another JSON allocation while memory is exhausted. These changes reduce duplication; they do not guarantee that an arbitrarily large dataset fits in Pyodide.

Detailed run progress is shown in the main workspace log. Batch alignment and concatenation share the main workflow’s structured progress format. Concatenation streams scanning, concatenating and writing phases with locus/taxon counts and elapsed time; failed concatenation is never labelled complete. The status beside Repository and the Workflow status label show only the broad phase (Checking input, Processing locally or Saving output), plus completion, cancellation or errors; locus names, aligner steps and elapsed-time ticks do not change these labels. Input preparation reports resolved path counts, read file counts and elapsed time before worker execution; output saving reports completed file counts and the current file. Pathlist folder roots are indexed once instead of rescanned for each entry.


### Incremental saving to an authorized folder

When an output folder is selected, completed locus alignments are queued for local saving during the alignment loop. In the full pipeline, their exported paths are `alignments/<filename>`; internal batch paths and MEMFS files remain unchanged. Completed QC/report branches are exported between raw/trimmed processing steps. Final matrices are exported only after normal staging, validation and publication. No internal files needed by later steps are deleted.

The per-run save queue caches shared directory promises and permits at most four concurrent file writes. A file counts as saved only after its stream closes successfully. Repeated identical outputs are skipped; changed content at the same path is written after any preceding write finishes. A save failure pauses new writes while active writes settle. Analysis results remain in the workspace; **Save remaining to folder** retries failed/pending files without rewriting confirmed successes. Changing the destination starts a separate save session for that folder. Cancelled/failed analyses retain already produced files for inspection/retry; these are intermediate outputs, not a completed run.

Local-save counts are displayed separately from analysis progress. Completion waits for all pending saves; incomplete saving produces a warning rather than a success claim. Without an authorized folder, results remain in the browser workspace for preview and manual download; no automatic downloads occur. This overlaps analysis and saving but retains data needed for downstream processing/preview, and does not by itself remove the browser's memory limit.

### PartitionFinder configuration filenames

CLI and browser outputs both use `*_partition_finder.cfg.txt`. The final `.txt` extension avoids browser restrictions on creating `.cfg` files in an authorized output folder. Before using the configuration with PartitionFinder, rename it to `partition_finder.cfg` in the analysis directory. Its contents and matrix references are unchanged.

### Disk-backed saved outputs

With an authorized output folder, outputs use shared storage references. Only after a writable stream closes successfully and a disk snapshot is obtained is the corresponding text removed from the result list and save queue. Preview reads a bounded slice; downloads, ZIP export, tool handoff and copying to a different folder read the saved file on demand without permanently restoring its text. Files must remain accessible and unchanged; size/timestamp changes are reported, and alignment read-back into the analysis also verifies its SHA256.

The full pipeline reloads saved alignments for the current raw/trimmed branch only. Each renamed alignment is exported as soon as that locus is renamed, before QC/concatenation; its OPFS working copy remains available for downstream reads. Branch completion exports QC/reports and releases eligible working alignments after confirmed local saves. The worker removes each eligible MEMFS file only after receiving confirmation of a successful local save. Renamed outputs remain in the exported file manifest and can be previewed/downloaded from disk. Failed saves retain their text and analysis files for retry. Without an output folder, the existing memory-backed behavior is retained. Final matrices stay staged until normal publication; this optimization does not remove Pyodide's memory limit or all input/matrix memory use.

### Per-locus temporary disk workspace

With an authorized output folder, the full pipeline uses a synchronous OPFS-backed working filesystem when supported. Python keeps the same paths and ordinary file interfaces, but original/renamed alignments, QC working files, concatenation spools and staged matrices are read and written in 4 KiB disk blocks instead of accumulating file contents in MEMFS. Each written file is flushed when closed; later QC and concatenation read the current locus as needed. `rename_taxa.py` also streams records rather than retaining the entire locus in a list; ID sets and optional ID-report rows still occupy memory.

This temporary working cache is browser-managed OPFS storage, separate from the user-selected output folder. Final outputs continue to use the existing validated save/publication flow and the same names. The cache is removed after success/failure; a later run reclaims abandoned caches after cancellation, while exclusive file locks protect active workers/tabs. Disk quota/I/O failures stop analysis with a temporary-disk diagnostic. If the synchronous API is unavailable, the existing branch-level memory mode is retained and the log reports the fallback. The run manifest records `workspace_storage`; `cached_file_bytes` describes working-file buffers only, not total Python/browser RAM.

This trades additional local disk I/O and quota usage for lower retained memory. Active-locus sequence strings, filesystem metadata and Python objects still consume memory; neither OPFS nor streaming renaming raises the WebAssembly heap limit. Renamed files reach the authorized output folder one locus at a time, while their working copies are disk-backed from the first write.

### Bounded matrix derivatives and format conversion

Codon-position splitting, fourfold-site extraction and PHYLIP/NEXUS conversion use a shared disk spool and bounded sequence blocks, including FASTA records written on one very long line. They do not materialize the complete supermatrix as Python sequence lists. Fourfold extraction intersects eligible columns in blocks before emitting third positions; export normalization edits are spooled rather than retained as a list. NEXUS retains a small page per taxon for interleaving. Temporary files are placed beside outputs so they use the OPFS workspace during disk-backed finalization; CLI stdout conversion uses the system temporary directory. Fourfold output validation also uses the bounded reader. These changes reduce matrix-sized allocations without increasing the WebAssembly heap limit.

### Alignment filenames

CLI and browser alignment outputs use the explicit FASTA suffix: `GENE.ALIGNER.TYPE.aln.fasta` for untrimmed alignments and `GENE.ALIGNER.TYPE.trimmed.aln.fasta` for trimmed alignments. Batch pathlists are named `all.TYPE.aln.fasta.pathlist` and `all.TYPE.trimmed.aln.fasta.pathlist`. The file contents remain FASTA; QC, provenance matching, paired trimming and concatenation use the same naming convention.

### Portable exports and bounded tool results

All Python tool results, including OG extraction and per-locus MSAP completion, return file metadata and use acknowledged 1 MiB output chunks instead of serializing all contents into one JSON response. Full-pipeline checkpoints fingerprint files with SHA256; unchanged files are reused in the final results instead of being transmitted and decoded again. Changed reports are retransmitted.

Exported output pathlists and QC pass/fail/retained lists use paths relative to their own directory. Export copies do not alter engine paths used for QC or concatenation. Lists referring to original input FASTAs use basenames because the originals are not bundled; place those original files alongside the list before CLI reuse. Batch alignment lists use `../<alignment filename>`. Download the complete output folder to preserve these relative paths.

Starting another run clears the previous result list and preview immediately. Duplicate input-stem checking uses a set, and the CLI indexes existing output filenames once before checking prefix ownership.

Manual bulk saving holds a fixed writer and file snapshot until completion. It disables Run, folder selection and other export controls, restoring them on success or failure. Worker crashes, undecodable messages and incomplete results discard the affected worker; retries initialize a fresh instance, and late errors from a discarded instance cannot interrupt a new run. Large browser TSV merges append rows without argument expansion. Python sequence-change report merging reads headers first and streams records into an atomic output, preserving mixed nucleotide/protein columns and the previous completed report if reading fails.

### Verified sources and fixed run provenance

The build emits immutable `engine/by-sha256/<hash>/<name>` source assets, including `browser_runner.py`. Before installing or executing Python code, the worker verifies the manifest bundle digest and every downloaded source SHA256. A mismatched cached source is revalidated once against the network; a persistent mismatch stops execution with a reload/retry diagnostic. The service worker revalidates the engine manifest against the network, with a cached fallback when offline; cached source contents must still pass verification. `no-store` requests bypass the cache.

Each run captures parameters, output naming and its equivalent CLI command before reading inputs. Input, parameter and tool-selection controls remain locked until analysis/saving ends; previous disabled states are retained. Output sizes and SHA256 hashes computed by Python survive worker transfer, batching and final result assembly and are reused in `phyloprep-run.json`. Files without fingerprints are encoded and hashed once, with the computed metadata reused afterward.

Browser incremental finalization leaves progress contexts, redirected output streams and engine overrides before yielding files to JavaScript. Each resumed Python call opens and closes its own analysis scope, so a different Pyodide task context cannot invalidate progress tokens. Early generator closure also restores temporary engine state without publishing staged matrices. Checkbox controls display their CLI option alongside the label, including `--longest-copy`, `-X`, `--allow-internal-stop`, `-n, --notrim` and optional translation validation via `-g`.

### Lazy inputs, binary exports and alignment reuse

Preflight samples at most 100 files, at most 1 MiB per file and 8 MiB total; larger inputs are explicitly deferred to full engine validation. Inputs remain File/Blob references after a small header check; they are staged in 1 MiB chunks only when needed. Finalization uses source names without reloading raw FASTA contents. Saved-file checksums and input fingerprints use incremental SHA-256. Outputs larger than 1 MiB remain binary Blobs during transfer instead of being decoded and joined into large strings; successful saves release their Blob references. Batch report merging uses Blob slices, and mapped audit reports are rewritten row by row. Empty workers are terminated after a workflow; the pipeline also releases alignment and Python workers between alignment and finalization because WebAssembly heaps cannot shrink. Recreating a worker reloads its runtime, using available browser caches.

With an authorized output folder, browser workflows first look for complete stage checkpoints in `reports/stage-checkpoints/`. A matching `alignment.json` skips the alignment-stage loop, including per-locus checkpoint and output prechecks. The fast match compares input names, sizes and modification times, alignment settings, relevant Python source checksums and aligner build metadata. It restores saved input fingerprints and creates lazy references to output files; it does not open each saved alignment. Inputs without reliable File modification times use the per-locus fallback. Changing a file while preserving both its size and modification time cannot be detected by this metadata fast match; use `--no-resume` when input metadata is not trustworthy.

Without a matching complete alignment checkpoint, loci are checked individually using input content and saved-output SHA-256 hashes. Per-locus records under `reports/alignment-cache/` support interrupted runs; a compatible prior `phyloprep-run.json` can also supply fingerprints. A complete stage checkpoint is written only after every required alignment and audit report has been saved. Before recomputing a stage, its old completion marker is removed. Plain alignment files without matching provenance are not accepted by filename alone.

Renaming uses the same stage-first approach, with separate `rename-raw.json` and `rename-trimmed.json` checkpoints. These bind the source alignment fingerprints, mapping content and renaming implementation. A complete match bypasses per-file rename prechecks and reuses the saved targets. If the locus set differs, matching entries in the recorded catalog are checked individually. Without a rename checkpoint, a compatible prior `phyloprep-run.json` can supply source, mapping and renamed-output fingerprints for per-file verification; entries without trustworthy metadata are renamed again. Each rename checkpoint is saved after all targets in that branch have been saved, before QC and matrix construction, so a subsequent QC failure does not discard completed renaming. QC thresholds and missing-taxon policies do not invalidate alignment or rename checkpoints; mapping changes invalidate renaming only. QC and matrices always run again.

Reused alignments and renamed targets are read from the chosen folder when downstream processing needs them, with size and SHA-256 validation at that point. This required reading is distinct from the skipped precheck loop. Missing or changed saved data on a complete-checkpoint path stops the run with a diagnostic; select **Recompute existing alignments** (`--no-resume`) to rebuild alignment, trimming and renaming. Output preview/download also reads saved data on demand. Keep the folder accessible and unchanged during analysis. These checkpoints resume completed tasks, not an interrupted aligner operation, and are separate from CLI state files. Temporary storage and Blob processing reduce duplicate data but do not remove WebAssembly memory limits or browser storage quotas.

### Missing taxa and padding provenance

The default `--missing-taxa pad-gaps` policy fills absent taxon sequences with `?`, not `-`; the dropdown explains the fill character. Original alignment gaps remain `-`. The shared Python exporters preserve `?` in FASTA/PHYLIP/NEXUS and codon-position derivatives, so browser and CLI outputs agree. Each full supermatrix includes `supermatrix.missing.tsv`, a sparse table of artificially filled taxon/locus intervals (1-based inclusive full-matrix coordinates, nt for nucleotide/codon data and aa for protein data). It is published and saved with the matrices, including independent/paired and raw/trimmed branches. Empty padding sets produce a header-only table. QC details show a separate `Missing_ratio`; provenance requires the missing-block table rather than counting all `?` as newly generated. Source-input normalization and fourfold rules for ambiguous and gapped observed codons remain unchanged. `phyloprep-run.json` records the padding symbols and metadata suffix. Alignment reuse still verifies relevant source checksums; downstream matrices are rebuilt.

`pad-gaps` is the default in both the full pipeline and standalone concatenation tool; select `skip-gene` explicitly for complete taxon coverage. When a padded codon matrix has no sites satisfying the selected fourfold rule, only that derived dataset is skipped, with `fourfold.skipped.txt` explaining why; other matrices remain available. Taxon mapping is still required when FASTA headers are gene IDs.

The pipeline and standalone fourfold extractor expose `--fourfold-missing-policy`: **Ignore complete missing codons (???)** (`ignore-missing`, default) or **Require all taxa to be evaluable** (`require-all`). The default requires at least two observed taxa with qualifying codons and matching first two bases, preserves missing output bases as `?`, and never ignores NNN, gaps, partial missing codons or conflicting observed codons. The setting is included in CLI command previews and run/report provenance. The pipeline hides this control in nucleotide/protein modes.
