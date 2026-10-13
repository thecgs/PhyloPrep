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
a merged sequence-change report. `phyloprep.py` browser mode reuses those
alignments and invokes the repository's original QC, taxon-renaming,
supermatrix, codon-position, four-fold and FASTA/PHYLIP/NEXUS modules inside
Pyodide.

Current browser constraints are explicit in the parameter UI: MUSCLE 5.1.0 is
the recommended default; MAFFT E-INS-i 7.520 is optional. The browser uses one
alignment worker, one locus at a time, and
`trimAlnSeq.py` trimming. Pseudogene/MACSE, PRANK, ClustalW2 and trimAl remain
CLI-only.


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
