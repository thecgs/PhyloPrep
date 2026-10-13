import { readFile } from "node:fs/promises";

const requiredFiles = [
  "site/index.html",
  "site/app.js",
  "site/pyodide-worker.js",
  "site/alignment-worker.js",
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
if (!app.includes('new Worker("alignment-worker.js")')) throw new Error("Dedicated cancellable alignment worker is missing.");

console.log(`PhyloPrep Web smoke checks passed for ${visibleModules.length} visible modules.`);


if (!app.includes("function rawParameterValues()")) throw new Error("Raw runtime parameter collection is missing.");
if (!app.includes("const optionValues = rawParameterValues();")) throw new Error("Browser runtime must use raw UI parameter values.");
if (app.includes("protein_stop_symbol: cliQuote(")) throw new Error("Protein stop symbol must not be CLI-quoted in runtime options.");
if (!app.includes("protein_stop_symbol_cli: cliQuote(")) throw new Error("Protein stop symbol CLI quoting is missing.");
if (!app.includes("taxon_separator_cli: cliQuote(")) throw new Error("Taxon separator CLI quoting is missing.");


const alignmentWorker = contents["site/alignment-worker.js"];
if (!alignmentWorker.includes("await cli.cat(outputPath)")) throw new Error("MUSCLE output must be read directly from the Aioli filesystem.");
if (alignmentWorker.includes("cli.download(outputPath)")) throw new Error("Worker-scoped Aioli blob URLs must not be used for MUSCLE output.");


if (!app.includes('script: "phyloprep.py"') || !app.includes('script: "MSAP_batch.py"')) throw new Error("Batch workflow modules are missing.");
const phyloprepBlock = app.slice(app.indexOf('script: "phyloprep.py"'), app.indexOf('script: "MSAP_batch.py"'));
const batchBlock = app.slice(app.indexOf('script: "MSAP_batch.py"'), app.indexOf('script: "extract_orthofinder_orthogroups.py"'));
if (!phyloprepBlock.includes('multiple: true, pathList: true')) throw new Error("phyloprep.py multi-file input must use browser pathlist mode.");
if (!batchBlock.includes('multiple: true, pathList: true')) throw new Error("MSAP_batch.py multi-file input must use browser pathlist mode.");
if (!app.includes('name: "path-lists/inputs.pathlist"')) throw new Error("MSAP_batch browser input pathlist output is missing.");
if (!app.includes('name: "inputs.pathlist"')) throw new Error("PhyloPrep browser input pathlist output is missing.");

if (!alignmentWorker.includes('cli.exec("muscle -align "')) throw new Error("MUSCLE must use Aioli single-string exec form.");
if (alignmentWorker.includes('cli.exec("muscle", [')) throw new Error("Two-argument Aioli exec form must not be used for MUSCLE.");

if (!app.includes('["muscle", "MUSCLE 5.1.0 (recommended)"]')) throw new Error("MUSCLE must be the recommended browser default.");
if (!app.includes('["mafft", "MAFFT E-INS-i 7.520 (Biowasm)"]')) throw new Error("MAFFT E-INS-i 7.520 must be exposed as an optional browser aligner.");
if (!worker.includes('option(values.align_software, "muscle")')) throw new Error("MUSCLE must be the browser fallback aligner.");
if (!alignmentWorker.includes('program: "tbfast"') || !alignmentWorker.includes('program: "dvtditr"')) throw new Error("Biowasm MAFFT core programs are missing from the alignment worker.");
if (alignmentWorker.includes('"mafft/7.520"')) throw new Error("Biowasm MAFFT has no mafft.js wrapper; load tbfast/dvtditr explicitly.");

if (!app.includes("CLI -s mafft runs upstream mafft --auto")) throw new Error("The Web/CLI MAFFT difference must be disclosed in the command preview.");

const missingTaxaOptions = 'options: [["skip-gene", "skip-gene"], ["pad-gaps", "pad-gaps"]]';
if (app.split(missingTaxaOptions).length - 1 < 2) {
  throw new Error("phyloprep.py and get_supergenes.py must use identical skip-gene/pad-gaps Missing taxa policy labels.");
}

if (!alignmentWorker.includes('{ tool: "muscle", version: "5.1.0", program: "muscle", reinit: true }')) {
  throw new Error("MUSCLE must reinitialize after each Biowasm invocation so multi-locus runs do not reuse parser/global state.");
}

if (!alignmentWorker.includes("async function createMafftCli()")) {
  throw new Error("MAFFT must create a fresh Aioli instance for each locus.");
}
if (!alignmentWorker.includes('if (tool === "mafft") {\n        try { await cli.close(); } catch {}')) {
  throw new Error("MAFFT Aioli instances must be closed after each locus.");
}
if (alignmentWorker.includes('let cliPromise;')) {
  throw new Error("MAFFT and MUSCLE must not share one persistent Aioli instance.");
}

if (!alignmentWorker.includes('"base/1.0.0"')) {
  throw new Error("Aioli MUSCLE reinit requires base/1.0.0 as the persistent filesystem module.");
}
if (!worker.includes("function pythonJson(value)")) {
  throw new Error("Pyodide worker must convert JavaScript objects to Python through json.loads().");
}
if (!worker.includes("run_phyloprep_finalize(${pythonJson(manifest)}, ${pythonJson(request.options)}")) {
  throw new Error("PhyloPrep finalizer must not inject JSON booleans directly as Python literals.");
}

if (!app.includes("No map.tsv selected and no FASTA IDs are shared across all loci")) {
  throw new Error("PhyloPrep preflight must warn when mapping is absent and locus IDs do not overlap.");
}
if (!app.includes("if (result.partial)")) {
  throw new Error("Recoverable PhyloPrep finalizer failures must render diagnostic outputs instead of a raw traceback.");
}

if (!html.includes('id="clear-input-files"')) {
  throw new Error("Input selection needs a one-click Clear all control.");
}
if (!app.includes("function clearAllInputFiles()") || !app.includes('elements.clearInputFiles.addEventListener("click", clearAllInputFiles)')) {
  throw new Error("Clear all must remove every selected input file and refresh the workspace.");
}

if (!app.includes("function normalizeCatalogText(value)") || !app.includes("function catalogSearchText(module)")) {
  throw new Error("Tool search must normalize punctuation and search across module metadata.");
}
if (!app.includes('elements.scriptSearch.addEventListener("search", refreshCatalogSearch)')) {
  throw new Error("Tool search must respond to the native search-field clear action.");
}

const styles = await readFile("site/styles.css", "utf8");
if (!/\.output-files-list[\s\S]*?resize:\s*vertical/.test(styles)) {
  throw new Error("Output files panel must be vertically resizable.");
}
if (!/\.output-preview\s*\{[\s\S]*?resize:\s*vertical/.test(styles)) {
  throw new Error("Output preview panel must be vertically resizable.");
}
