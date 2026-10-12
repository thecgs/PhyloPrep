const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.2/full/";
const ENGINE_ROOT = "/engine";
let runtimePromise;
let engineManifest;

function safeRelativePath(value, fallback) {
  const pieces = String(value || fallback).replaceAll("\\", "/").split("/")
    .filter((piece) => piece && piece !== "." && piece !== "..");
  return pieces.join("/") || fallback;
}

function uniqueFileName(name, used) {
  const safeName = safeRelativePath(name, "input.fasta").split("/").at(-1);
  const dot = safeName.lastIndexOf(".");
  const stem = dot > 0 ? safeName.slice(0, dot) : safeName;
  const extension = dot > 0 ? safeName.slice(dot) : "";
  let candidate = safeName;
  let number = 2;
  while (used.has(candidate)) candidate = `${stem}-${number++}${extension}`;
  used.add(candidate);
  return candidate;
}

function isFasta(file) {
  return file.content.trimStart().startsWith(">");
}

async function fetchText(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path}.`);
  return response.text();
}

async function installOriginalEngine(runtime) {
  const manifestResponse = await fetch("engine/manifest.json", { cache: "no-cache" });
  if (!manifestResponse.ok) {
    throw new Error("The original-script bundle is missing. Run node site/build-engine.mjs before local preview.");
  }
  engineManifest = await manifestResponse.json();
  if (!Array.isArray(engineManifest.files) || !engineManifest.files.length) {
    throw new Error("The original-script bundle manifest is invalid.");
  }
  runtime.FS.mkdirTree(ENGINE_ROOT);
  const files = await Promise.all(engineManifest.files.map(async ({ name }) => [name, await fetchText(`engine/${name}`)]));
  files.forEach(([name, source]) => runtime.FS.writeFile(`${ENGINE_ROOT}/${name}`, source));
  runtime.FS.writeFile("/browser_runner.py", await fetchText("browser_runner.py"));
  await runtime.runPythonAsync("exec(compile(open('/browser_runner.py', encoding='utf-8').read(), '/browser_runner.py', 'exec'), globals())");
}

async function ensureRuntime(runId) {
  if (!runtimePromise) {
    runtimePromise = (async () => {
      self.postMessage({ type: "stage", runId, stage: "runtime", label: "Loading Pyodide and original scripts" });
      importScripts(`${PYODIDE_INDEX_URL}pyodide.js`);
      const runtime = await loadPyodide({ indexURL: PYODIDE_INDEX_URL });
      await runtime.loadPackage("biopython");
      await installOriginalEngine(runtime);
      return runtime;
    })();
  }
  return runtimePromise;
}

function stageBrowserInputs(runtime, request) {
  const runRoot = `/work/runs/run-${request.runId}`;
  const inputRoot = `${runRoot}/input`;
  const optionsRoot = `${runRoot}/options`;
  const outputRoot = `${runRoot}/output`;
  runtime.FS.mkdirTree(inputRoot);
  runtime.FS.mkdirTree(optionsRoot);
  runtime.FS.mkdirTree(outputRoot);
  const names = new Set();
  const files = request.files.map((file) => {
    const name = uniqueFileName(file.name, names);
    const path = `${inputRoot}/${name}`;
    runtime.FS.writeFile(path, file.content);
    return { ...file, name, path, fasta: isFasta(file) };
  });
  const fastaFiles = files.filter((file) => file.fasta);
  let pathList;
  if (request.usesPathList && fastaFiles.length > 1) {
    pathList = `${optionsRoot}/${request.script.replace(/\.py$/i, "")}.pathlist`;
    runtime.FS.writeFile(pathList, `${fastaFiles.map((file) => file.path).join("\n")}\n`);
  }
  const writeTaxaList = (key) => {
    const values = request.options[key];
    if (!Array.isArray(values) || !values.length) return undefined;
    const path = `${optionsRoot}/${key.replaceAll("_", "-")}.txt`;
    runtime.FS.writeFile(path, `${values.join("\n")}\n`);
    return path;
  };
  return {
    files,
    fastaFiles,
    plainFiles: files.filter((file) => !file.fasta),
    pathList,
    keepTaxa: writeTaxaList("keep_taxa"),
    removeTaxa: writeTaxaList("remove_taxa"),
    optionsRoot,
    outputRoot,
  };
}

function outputPath(stage, outputName, fallback) {
  return `${stage.outputRoot}/${safeRelativePath(outputName, fallback)}`;
}

function inputFile(stage, index = 0) {
  const file = stage.fastaFiles[index] || stage.files[index];
  if (!file) throw new Error("No FASTA input file was staged.");
  return file.path;
}

function listInputs(stage) {
  if (!stage.fastaFiles.length) throw new Error("No FASTA input file was staged.");
  return stage.pathList ? [stage.pathList] : stage.fastaFiles.map((file) => file.path);
}

function option(value, fallback = "") {
  return value === undefined || value === null || value === "" ? fallback : String(value);
}

function originalArguments(request, stage) {
  const values = request.options;
  const output = outputPath(stage, request.outputName, "output.txt");
  switch (request.script) {
    case "AA2Codon.py":
      return ["-c", inputFile(stage, 0), "-p", inputFile(stage, 1), "-o", output, "-g", option(values.genetic_code, "1")];
    case "alignment_qc.py":
      return ["-i", ...listInputs(stage), "-p", outputPath(stage, request.outputName, "phyloprep-qc"), "-st", option(values.seqtype, "nucl"), "-g", option(values.genetic_code, "1"), "--min-length", option(values.min_length, "1"), "--max-n-ratio", option(values.max_n_ratio, "1"), "--max-x-ratio", option(values.max_x_ratio, "1"), "--max-gap-ratio", option(values.max_gap_ratio, "1"), "--min-parsimony-informative-sites", option(values.min_parsimony_informative_sites, "0"), ...(values.allow_terminal_stop ? ["--allow-terminal-stop"] : []), ...(values.allow_internal_stop ? ["--allow-internal-stop"] : [])];
    case "calulate_4dtv_and_correction.py":
      return ["-i", ...listInputs(stage), "-o", output, "-g", option(values.genetic_code, "1")];
    case "extract_4-fold_degenerated_sites.py":
      return ["-i", inputFile(stage), "-o", output, "-g", option(values.genetic_code, "1")];
    case "get_supergenes.py":
      return ["-i", ...listInputs(stage), "-p", outputPath(stage, request.outputName, "supergene"), "--missing-taxa", option(values.missing_taxa, "skip-gene")];
    case "filter_sequences.py":
      return ["-i", inputFile(stage), "-o", output, "-st", option(values.seqtype, "nucl"), "--min-length", option(values.min_length, "0"), "--max-gap-ratio", option(values.max_gap_ratio, "1"), "--max-n-ratio", option(values.max_n_ratio, "1"), "--max-x-ratio", option(values.max_x_ratio, "1"), ...(stage.keepTaxa ? ["--keep-taxa", stage.keepTaxa] : []), ...(stage.removeTaxa ? ["--remove-taxa", stage.removeTaxa] : [])];
    case "rename_taxa.py": {
      const mapping = stage.plainFiles[0]?.path;
      const mapOutput = values.map_output ? outputPath(stage, `${safeRelativePath(request.outputName, "renamed.fasta").replace(/\.[^.]+$/, "")}.map.tsv`, "renamed.map.tsv") : undefined;
      return ["-i", inputFile(stage), "-o", output, ...(mapping ? ["-m", mapping] : []), ...(mapOutput ? ["--map-output", mapOutput] : []), ...(values.allow_duplicate_ids ? ["--allow-duplicate-ids"] : []), ...(values.keep_special ? ["--keep-special"] : [])];
    }
    case "trimAlnSeq.py":
      return ["-i", inputFile(stage), "-o", output, "-st", option(values.seqtype, "codon"), "-G", option(values.gap_ratio, "0.2"), "-N", option(values.n_ratio, "0.2"), "-X", option(values.x_ratio, "0.2"), ...(values.removeall ? ["-ra"] : [])];
    case "fasta2axt.py":
    case "fasta2nex.py":
    case "fasta2phy.py":
      return ["-i", inputFile(stage), "-o", output, "-st", option(values.seqtype, "DNA")];
    default:
      throw new Error(`Unsupported original module: ${request.script}`);
  }
}


function msapArguments(request, stage) {
  const values = request.options;
  return [
    "-i", inputFile(stage),
    "-st", option(values.seqtype, "codon"),
    "-s", "muscle",
    "-t", "1",
    "-g", option(values.genetic_code, "1"),
    "-ts", "trimAlnSeq",
    "-G", option(values.gap_ratio, "0.2"),
    "-N", option(values.n_ratio, "0.2"),
    "-X", option(values.x_ratio, "0.2"),
    ...(values.seqtype === "prot" ? ["--protein-stop-symbol", option(values.protein_stop_symbol, "*")] : []),
    ...(values.notrim ? ["-n"] : []),
  ];
}

self.addEventListener("message", async (event) => {
  const request = event.data;
  if (!["run", "msap-probe", "msap-complete", "phyloprep-finalize"].includes(request.type)) return;
  const { runId } = request;
  try {
    const runtime = await ensureRuntime(runId);
    const stage = stageBrowserInputs(runtime, request);

    if (request.type === "msap-probe") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Preparing MSAP.py MUSCLE input" });
      const args = msapArguments(request, stage);
      const raw = await runtime.runPythonAsync(`run_msap_probe(${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)})`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "MSAP.py could not prepare the MUSCLE input.");
      self.postMessage({ type: "result", runId, result });
      return;
    }

    if (request.type === "msap-complete") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Finishing original MSAP.py workflow" });
      const alignmentPath = `${stage.optionsRoot}/biowasm-muscle.aln.fasta`;
      runtime.FS.writeFile(alignmentPath, request.alignmentContent || "");
      const args = msapArguments(request, stage);
      const raw = await runtime.runPythonAsync(`run_msap_complete(${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)}, ${JSON.stringify(alignmentPath)})`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "The original MSAP.py workflow could not complete.");
      self.postMessage({
        type: "result",
        runId,
        result: {
          ...result,
          provenance: `Untouched MSAP.py + Biowasm MUSCLE 5.1.0 · source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
        },
      });
      return;
    }

    if (request.type === "phyloprep-finalize") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Running PhyloPrep QC and matrix construction" });
      const manifest = stage.files.map(({ name, path, role, inputKey, fasta }) => ({ name, path, role, inputKey, fasta }));
      const raw = await runtime.runPythonAsync(`run_phyloprep_finalize(${JSON.stringify(manifest)}, ${JSON.stringify(request.options)}, ${JSON.stringify(stage.outputRoot)})`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "The browser PhyloPrep pipeline could not complete.");
      self.postMessage({
        type: "result",
        runId,
        result: {
          ...result,
          provenance: `Browser phyloprep.py adapter · original downstream modules · source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
        },
      });
      return;
    }

    self.postMessage({ type: "stage", runId, stage: "process", label: "Running original script locally" });
    const args = originalArguments(request, stage);
    const raw = await runtime.runPythonAsync(`run_original_module(${JSON.stringify(request.script)}, ${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)})`);
    const result = JSON.parse(raw);
    if (!result.ok) throw new Error(result.error || "The original script could not complete.");
    self.postMessage({
      type: "result",
      runId,
      result: {
        ...result,
        provenance: `Original source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
      },
    });
  } catch (error) {
    self.postMessage({ type: "error", runId, message: error.message || String(error) });
  }
});
