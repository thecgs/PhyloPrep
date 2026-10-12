import { readFile } from "node:fs/promises";

const requiredFiles = [
  "site/index.html",
  "site/app.js",
  "site/pyodide-worker.js",
  "site/muscle-worker.js",
  "site/browser_runner.py",
  "site/build-engine.mjs",
];

const contents = Object.fromEntries(await Promise.all(requiredFiles.map(async (path) => [path, await readFile(path, "utf8")])));

const app = contents["site/app.js"];
const worker = contents["site/pyodide-worker.js"];
const build = contents["site/build-engine.mjs"];
const html = contents["site/index.html"];

const visibleModules = [
  "phyloprep.py",
  "MSAP_batch.py",
  "extract_orthofinder_orthogroups.py",
  "oma_to_orthogroups.py",
  "orthomcl_to_orthogroups.py",
  "proteinortho_to_orthogroups.py",
  "sonicparanoid_to_orthogroups.py",
  "split_codon_seqence_alignment.py",
  "AA2Codon.py",
  "alignment_qc.py",
  "calulate_4dtv_and_correction.py",
  "extract_4-fold_degenerated_sites.py",
  "get_supergenes.py",
  "filter_sequences.py",
  "rename_taxa.py",
  "trimAlnSeq.py",
  "fasta2axt.py",
  "fasta2nex.py",
  "fasta2phy.py",
];

for (const script of visibleModules) {
  if (!app.includes(`script: "${script}"`)) throw new Error(`Missing browser module: ${script}`);
}
if (app.includes('script: "MSAP.py", title:')) throw new Error("MSAP.py must remain hidden from the browser catalog.");
if (app.includes('script: "orthogroup_table.py", title:')) throw new Error("orthogroup_table.py must remain an internal dependency.");
for (const internal of ["MSAP.py", "orthogroup_table.py"]) {
  if (!build.includes(`"${internal}"`)) throw new Error(`Missing required internal engine source: ${internal}`);
}
for (const selector of ["script-search", "preflight-summary", "use-next-step"]) {
  if (!html.includes(`id="${selector}"`)) throw new Error(`Missing workspace UI element: ${selector}`);
}
if (!worker.includes("removeTree(runtime, stage.runRoot)")) throw new Error("Pyodide run workspace cleanup is missing.");
if (!app.includes("phyloprep-run.json")) throw new Error("Run provenance manifest is missing.");
if (!app.includes('new Worker("muscle-worker.js")')) throw new Error("Dedicated cancellable MUSCLE worker is missing.");

console.log(`PhyloPrep Web smoke checks passed for ${visibleModules.length} visible modules.`);
