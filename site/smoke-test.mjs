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


if (!app.includes("function rawParameterValues()")) throw new Error("Raw runtime parameter collection is missing.");
if (!app.includes("const optionValues = rawParameterValues();")) throw new Error("Browser runtime must use raw UI parameter values.");
if (app.includes("protein_stop_symbol: cliQuote(")) throw new Error("Protein stop symbol must not be CLI-quoted in runtime options.");
if (!app.includes("protein_stop_symbol_cli: cliQuote(")) throw new Error("Protein stop symbol CLI quoting is missing.");
if (!app.includes("taxon_separator_cli: cliQuote(")) throw new Error("Taxon separator CLI quoting is missing.");


const muscleWorker = contents["site/muscle-worker.js"];
if (!muscleWorker.includes("await cli.cat(outputPath)")) throw new Error("MUSCLE output must be read directly from the Aioli filesystem.");
if (muscleWorker.includes("cli.download(outputPath)")) throw new Error("Worker-scoped Aioli blob URLs must not be used for MUSCLE output.");


if (!app.includes('script: "phyloprep.py"') || !app.includes('script: "MSAP_batch.py"')) throw new Error("Batch workflow modules are missing.");
const phyloprepBlock = app.slice(app.indexOf('script: "phyloprep.py"'), app.indexOf('script: "MSAP_batch.py"'));
const batchBlock = app.slice(app.indexOf('script: "MSAP_batch.py"'), app.indexOf('script: "extract_orthofinder_orthogroups.py"'));
if (!phyloprepBlock.includes('multiple: true, pathList: true')) throw new Error("phyloprep.py multi-file input must use browser pathlist mode.");
if (!batchBlock.includes('multiple: true, pathList: true')) throw new Error("MSAP_batch.py multi-file input must use browser pathlist mode.");
if (!app.includes('name: "path-lists/inputs.pathlist"')) throw new Error("MSAP_batch browser input pathlist output is missing.");
if (!app.includes('name: "inputs.pathlist"')) throw new Error("PhyloPrep browser input pathlist output is missing.");

if (!muscleWorker.includes('cli.exec("muscle -align "')) throw new Error("MUSCLE must use Aioli single-string exec form.");
if (muscleWorker.includes('cli.exec("muscle", [')) throw new Error("Two-argument Aioli exec form must not be used for MUSCLE.");
