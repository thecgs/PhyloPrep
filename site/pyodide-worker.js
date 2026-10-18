const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.2/full/";
const ENGINE_ROOT = "/engine";
let runtimePromise;
let engineManifest;

const outputAcknowledgements = new Map();
let outputSequence = 0;

function acknowledgeOutput(runId, message, timeout, label) {
  const sequence = ++outputSequence;
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      outputAcknowledgements.delete(sequence);
      reject(new Error(`${label}: ${message.name}`));
    }, timeout);
    outputAcknowledgements.set(sequence, {runId, resolve: ack => {clearTimeout(timer); resolve(ack);}});
    self.postMessage({...message, runId, sequence});
  });
}

async function streamOutputFiles(runtime, runId, files, transferId, waitForSave = false, transferred = new Map()) {
  if (transferId && files.some(file => file.release)) {
    self.postMessage({type: 'stage', runId, stage: 'process',
      label: 'QC and matrix construction for this branch finished; checking saved outputs and clearing temporary copies'});
  }
  for (let index = 0; index < files.length; index++) {
    const file = files[index];
    if (!file.path) continue;
    if (file.sha256 && transferred.get(file.name) === file.sha256) {
      const message = { type: "output-file-reuse", index, name: file.name, transferId,
        release: Boolean(file.release), fileCount: files.length };
      if (file.release && waitForSave) {
        const ack = await acknowledgeOutput(runId, message, 300000, 'Local output reuse acknowledgement timed out');
        if (ack?.saved && runtime.FS.analyzePath(file.path).exists) runtime.FS.unlink(file.path);
      } else self.postMessage({...message, runId});
      continue;
    }
    if (file.sha256 && !globalThis.PhyloPrepHash) importScripts('browser-hash.js');
    const hash = file.sha256 ? new PhyloPrepHash.SHA256() : null;
    const stream = runtime.FS.open(file.path, "r");
    try {
      let position = 0;
      while (position < file.size) {
        const bytes = new Uint8Array(Math.min(1024 * 1024, file.size - position));
        const read = runtime.FS.read(stream, bytes, 0, bytes.length, position);
        if (!read) throw new Error(`Incomplete output transfer: ${file.name}`);
        const chunk = bytes.subarray(0, read);
        hash?.update(chunk);
        const sequence = ++outputSequence;
        await new Promise((resolve, reject) => {
          const timer = setTimeout(() => {
            outputAcknowledgements.delete(sequence);
            reject(new Error(`Output transfer acknowledgement timed out: ${file.name}`));
          }, 30000);
          outputAcknowledgements.set(sequence, { runId, resolve: () => { clearTimeout(timer); resolve(); } });
          self.postMessage({ type: "output-chunk", runId, sequence, index, name: file.name,
            fileCount: files.length, transferId, bytes: chunk }, [bytes.buffer]);
        });
        position += read;
      }
      if (hash && hash.hex() !== file.sha256) throw new Error(`Output checksum mismatch while reading browser temporary disk: ${file.name}. The incomplete output was not saved.`);
      const end = { type: "output-file-end", runId, index, name: file.name, transferId, size: file.size, sha256: file.sha256 };
      let acknowledgement;
      if (waitForSave) {
        acknowledgement = await acknowledgeOutput(runId, end, 300000, 'Local output queue acknowledgement timed out');
      } else self.postMessage(end);
      if (file.sha256) transferred.set(file.name, file.sha256);
      if (file.release && acknowledgement?.saved) runtime.FS.unlink(file.path);
    } finally {
      runtime.FS.close(stream);
    }
  }
}

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
  return file.fasta ?? String(file.content || "").trimStart().startsWith(">");
}

async function hashSourceBytes(bytes) {
  const digest = await crypto.subtle.digest("SHA-256", bytes);
  return [...new Uint8Array(digest)].map(byte => byte.toString(16).padStart(2, "0")).join("");
}

async function verifiedEngineSource(file) {
  if (!file || !/^[0-9a-f]{64}$/.test(file.sha256) || !file.url) {
    throw new Error("The engine source metadata is incomplete. Reload the page to load a matching build.");
  }
  for (let attempt = 0; attempt < 2; attempt++) {
    const response = await fetch(`engine/${file.url}`, attempt ? { cache: "no-cache" } : undefined);
    if (!response.ok) throw new Error(`Could not load verified source: ${file.name}. Reload the page and retry.`);
    const bytes = await response.arrayBuffer();
    if (await hashSourceBytes(bytes) === file.sha256) return new TextDecoder().decode(bytes);
  }
  throw new Error(`Source checksum mismatch: ${file.name}. Reload the page to refresh the engine; analysis was not started.`);
}

async function installOriginalEngine(runtime) {
  const manifestResponse = await fetch("engine/manifest.json", { cache: "no-cache" });
  if (!manifestResponse.ok) {
    throw new Error("The original-script bundle is missing. Run node site/build-engine.mjs before local preview.");
  }
  const manifest = await manifestResponse.json();
  if (!Array.isArray(manifest.files) || !manifest.files.length || !manifest.runner) {
    throw new Error("The original-script bundle manifest is incomplete. Reload the page to load a matching build.");
  }
  const records = [...manifest.files, manifest.runner];
  const bundleBytes = new TextEncoder().encode(records.map(({name, sha256}) => `${name}:${sha256}`).join("\n"));
  if (await hashSourceBytes(bundleBytes) !== manifest.bundleSha256) {
    throw new Error("The engine bundle checksum is invalid. Reload the page and retry.");
  }
  // Verify every source, including the browser adapter, before installing/executing any.
  const sources = await Promise.all(records.map(verifiedEngineSource));
  runtime.FS.mkdirTree(ENGINE_ROOT);
  manifest.files.forEach((file, index) => runtime.FS.writeFile(`${ENGINE_ROOT}/${file.name}`, sources[index]));
  runtime.FS.writeFile("/browser_runner.py", sources.at(-1));
  await runtime.runPythonAsync("exec(compile(open('/browser_runner.py', encoding='utf-8').read(), '/browser_runner.py', 'exec'), globals())");
  engineManifest = manifest;
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

async function readDiskInput(file) {
  const saved = file.disk;
  const source = await saved.handle.getFile();
  if (source.size !== saved.size || source.lastModified !== saved.lastModified) {
    throw new Error(`Saved input changed on disk: ${file.name}. Restore the original file before retrying.`);
  }
  if (!globalThis.PhyloPrepHash) importScripts('browser-hash.js');
  if (file.sha256 && await PhyloPrepHash.blob(source) !== file.sha256) throw new Error(`Saved input checksum mismatch: ${file.name}`);
  return source;
}

async function writeBlobInput(runtime, path, source) {
  const stream=runtime.FS.open(path,'w');
  try {
    for(let position=0;position<source.size;position+=1024*1024) {
      const bytes=new Uint8Array(await source.slice(position,position+1024*1024).arrayBuffer());
      let done=0;
      while(done<bytes.length) {
        const count=runtime.FS.write(stream,bytes,done,bytes.length-done,position+done);
        if(!count)throw Error(`Incomplete input write: ${path}`);
        done+=count;
      }
    }
  } finally {runtime.FS.close(stream);}
}

async function stageBrowserInputs(runtime, request) {
  const runRoot = `/work/runs/run-${request.runId}`;
  const inputRoot = `${runRoot}/input`;
  const optionsRoot = `${runRoot}/options`;
  const outputRoot = `${runRoot}/output`;
  runtime.FS.mkdirTree(inputRoot);
  runtime.FS.mkdirTree(optionsRoot);
  runtime.FS.mkdirTree(outputRoot);
  const names = new Set();
  const files = [];
  for (const file of request.files) {
    const name = uniqueFileName(file.name, names);
    const path = `${inputRoot}/${name}`;
    const deferred = request.incrementalSave && request.type === "phyloprep-finalize" && file.role === "batch-output" && file.disk && name.endsWith('.aln.fasta');
    if (request.type === 'phyloprep-finalize' && file.role === 'source') {
      // Finalization only needs source names; alignments and audit reports are supplied separately.
      runtime.FS.writeFile(path,'');
    } else if (!deferred) {
      if(file.disk || file.blob) await writeBlobInput(runtime,path,file.disk ? await readDiskInput(file) : file.blob);
      else runtime.FS.writeFile(path,file.content);
    }
    const staged = { ...file, name, path, on_disk: Boolean(deferred), fasta: file.fasta ?? (file.disk ? /\.(fa|fas|fasta)$/i.test(name) : isFasta(file)) };
    delete staged.content;
    delete staged.blob;
    delete file.content;
    delete file.blob;
    files.push(staged);
  }
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
    runRoot,
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

function removeTree(runtime, path) {
  try {
    const entries = runtime.FS.readdir(path).filter((name) => name !== "." && name !== "..");
    entries.forEach((name) => {
      const child = `${path}/${name}`;
      const stat = runtime.FS.stat(child);
      if (runtime.FS.isDir(stat.mode)) removeTree(runtime, child);
      else runtime.FS.unlink(child);
    });
    runtime.FS.rmdir(path);
  } catch {
    // A cancelled/failed run may already have a partially removed workspace.
  }
}

function pythonJson(value) {
  const json = JSON.stringify(value);
  return `__import__("json").loads(${JSON.stringify(json)})`;
}

function outputPath(stage, outputName, fallback) {
  return `${stage.outputRoot}/${safeRelativePath(outputName, fallback)}`;
}

function inputFile(stage, index = 0) {
  const file = stage.fastaFiles[index] || stage.files[index];
  if (!file) throw new Error("No FASTA input file was staged.");
  return file.path;
}

function fileByKey(stage, key, index = 0) {
  const matches = stage.files.filter((file) => file.inputKey === key);
  const file = matches[index];
  if (!file) throw new Error(`Missing staged input: ${key}`);
  return file;
}

function orthofinderArguments(runtime, request, stage) {
  const values = request.options;
  const orthogroups = fileByKey(stage, "orthogroups").path;
  const sequences = stage.files.filter((file) => file.inputKey === "sequences" && file.fasta);
  if (!sequences.length) throw new Error("No per-taxon FASTA files were staged.");
  const sequenceDir = `${stage.optionsRoot}/orthofinder-sequences`;
  runtime.FS.mkdirTree(sequenceDir);
  sequences.forEach((file) => {
    runtime.FS.rename(file.path, `${sequenceDir}/${file.name}`);
    file.path = `${sequenceDir}/${file.name}`;
  });
  const outdir = outputPath(stage, request.outputName, "one-to-one_orthologs");
  const args = [
    "-og", orthogroups,
    "-i", sequenceDir,
    "-o", outdir,
    "-m", option(values.copy_mode, "single-copy"),
    "--orthofinder-id-style", option(values.orthofinder_id_style, "auto"),
  ];
  const keep = String(values.keep_taxa_text || "").split(/[\s,]+/).filter(Boolean);
  const drop = String(values.drop_taxa_text || "").split(/[\s,]+/).filter(Boolean);
  if (keep.length) args.push("-k", ...keep);
  if (drop.length) args.push("-d", ...drop);
  if (values.min_taxa) args.push("--min-taxa", String(values.min_taxa));
  if (values.copy_mode === "all-copies" && values.max_copies_per_taxon !== undefined && values.max_copies_per_taxon !== "") args.push("--max-copies-per-taxon", String(values.max_copies_per_taxon));
  if (values.copy_mode === "all-copies" && values.longest_copy) args.push("--longest-copy", "-st", option(values.seqtype, "codon"));
  if (values.no_species_prefix) args.push("-X");
  const mapName = safeRelativePath(values.map_output_name || "map.tsv", "map.tsv").split("/").at(-1);
  const pathlistName = safeRelativePath(values.pathlist_output_name || "orthogroups.pathlist", "orthogroups.pathlist").split("/").at(-1);
  args.push("-mo", `${outdir}/${mapName}`);
  args.push("--pathlist-output", `${outdir}/${pathlistName}`);
  return args;
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
      return ["-c", inputFile(stage, 0), "-p", inputFile(stage, 1), "-o", output, ...(values.validate_translation ? ["-g", option(values.genetic_code, "1")] : [])];
    case "alignment_qc.py": {
      const changes = stage.files.find((file) => file.inputKey === "sequence_changes");
      return ["-i", ...listInputs(stage), "-p", outputPath(stage, request.outputName, "phyloprep-qc"), "-st", option(values.seqtype, "nucl"), "-g", option(values.genetic_code, "1"), ...(changes ? ["--sequence-changes", changes.path] : []), "--min-length", option(values.min_length, "1"), "--max-n-ratio", option(values.max_n_ratio, "1"), "--max-x-ratio", option(values.max_x_ratio, "1"), "--max-gap-ratio", option(values.max_gap_ratio, "1"), ...(values.seqtype === "prot" ? ["--min-protein-pis", option(values.min_protein_pis, "0")] : ["--min-nucl-pis", option(values.min_nucl_pis, "0")]), ...(values.allow_terminal_stop ? ["--allow-terminal-stop"] : []), ...(values.allow_internal_stop ? ["--allow-internal-stop"] : [])];
    }
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
      return ["-i", inputFile(stage), "-o", output, ...(mapping ? ["-m", mapping] : []), ...(mapOutput ? ["--map-output", mapOutput] : []), ...(values.allow_duplicate_ids ? ["--allow-duplicate-ids"] : [])];
    }
    case "trimAlnSeq.py":
      return ["-i", inputFile(stage), "-o", output, "-st", option(values.seqtype, "codon"), "-G", option(values.gap_ratio, "0.2"), "-N", option(values.n_ratio, "0.2"), "-X", option(values.x_ratio, "0.2"), ...(values.removeall ? ["-ra"] : [])];
    case "oma_to_orthogroups.py":
      return ["-i", fileByKey(stage, "input").path, "-o", output, "--taxon-separator", option(values.taxon_separator, "|"), ...(values.strip_taxon_prefix ? ["--strip-taxon-prefix"] : [])];
    case "orthomcl_to_orthogroups.py":
      return ["-i", fileByKey(stage, "input").path, "-o", output, "--taxon-separator", option(values.taxon_separator, "|"), ...(values.strip_taxon_prefix ? ["--strip-taxon-prefix"] : [])];
    case "proteinortho_to_orthogroups.py":
    case "sonicparanoid_to_orthogroups.py":
      return ["-i", fileByKey(stage, "input").path, "-o", output];
    case "split_codon_seqence_alignment.py":
      return ["-i", fileByKey(stage, "input").path];
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
    "-s", option(values.align_software, "mafft"),
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
  if (request.type === "output-chunk-ack") {
    const pending = outputAcknowledgements.get(request.sequence);
    if (pending?.runId === request.runId) {
      outputAcknowledgements.delete(request.sequence);
      pending.resolve(request);
    }
    return;
  }
  if (!["run", "msap-probe", "msap-complete", "phyloprep-finalize", "orthofinder-extract"].includes(request.type)) return;
  const { runId } = request;
  let runtime;
  let stage;
  let diskWorkspace;
  const transferred = new Map();
  try {
    runtime = await ensureRuntime(runId);
    if(request.expectedBundle && request.expectedBundle!==engineManifest.bundleSha256)throw Error("The engine build changed. Reload the page before reusing saved alignments.");
    if (request.incrementalSave) {
      if (!globalThis.PhyloPrepDiskFS) importScripts('browser-disk-fs.js');
      diskWorkspace = await PhyloPrepDiskFS.open(runtime, `/work/runs/run-${runId}`);
      self.postMessage({type:'stage',runId,stage:'process',label:diskWorkspace
        ? 'Using temporary disk for per-locus alignment processing'
        : 'Temporary disk is unavailable; using branch-level memory storage'});
    }
    stage = await stageBrowserInputs(runtime, request);

    if (request.type === "msap-probe") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Preparing MSAP.py alignment input" });
      const args = msapArguments(request, stage);
      const raw = await runtime.runPythonAsync(`run_msap_probe(${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)})`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "MSAP.py could not prepare the alignment input.");
      self.postMessage({ type: "result", runId, result });
      return;
    }

    if (request.type === "msap-complete") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Finishing original MSAP.py workflow" });
      const alignmentPath = `${stage.optionsRoot}/biowasm-alignment.aln.fasta`;
      runtime.FS.writeFile(alignmentPath, request.alignmentContent || "");
      delete request.alignmentContent;
      const args = msapArguments(request, stage);
      const raw = await runtime.runPythonAsync(`run_msap_complete(${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)}, ${JSON.stringify(alignmentPath)}, manifest_only=True)`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "The original MSAP.py workflow could not complete.");
      await streamOutputFiles(runtime, runId, result.files || [], undefined, request.incrementalSave, transferred);
      self.postMessage({
        type: "result",
        runId,
        result: {
          ...result,
          engineMetadata: engineManifest,
          provenance: `Untouched MSAP.py + ${request.options.align_software === "mafft" ? "MAFFT E-INS-i" : "MUSCLE"} (Biowasm package ${engineManifest.aligners[request.options.align_software || "mafft"].packageVersion}) · source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
        },
      });
      return;
    }

    if (request.type === "phyloprep-finalize") {
      // SQLite is optional report infrastructure. A package/network failure
      // must never stop alignment QC, matrix publication or output download.
      if (!request.options.no_report) {
        try { await runtime.loadPackage("sqlite3"); }
        catch (error) { request.options = { ...request.options, report_dependency_error:
          `Could not load the browser report dependency sqlite3: ${error.message || error}` }; }
      }
      self.postMessage({ type: "stage", runId, stage: "process", label: "Preparing downstream alignment processing" });
      const manifest = stage.files.map(({ name, path, role, inputKey, fasta, on_disk }) => ({ name, path, role, inputKey, fasta, on_disk }));
      runtime.globals.set("_phyloprep_progress_callback", (eventJson) => {
        self.postMessage({ type: "progress", runId, progress: JSON.parse(eventJson) });
      });
      let raw;
      let stepsStarted = false;
      try {
        if (request.incrementalSave) {
          await runtime.runPythonAsync(`_phyloprep_steps = _phyloprep_finalize_steps(${pythonJson(manifest)}, ${pythonJson(request.options)}, ${JSON.stringify(stage.outputRoot)}, progress_callback=_phyloprep_progress_callback, manifest_only=True, incremental=True, disk_backed=True)`);
          stepsStarted = true;
          let checkpoint = 0;
          while (true) {
            raw = await runtime.runPythonAsync("advance_phyloprep_finalize(_phyloprep_steps)");
            const step = JSON.parse(raw);
            if (!step.pending) break;
            if (step.load_variant) {
              for (const file of stage.files) {
                if (!file.on_disk || file.name.endsWith('.trimmed.aln.fasta') !== (step.load_variant === 'trimmed')) continue;
                const target = `${stage.outputRoot}/alignments/${file.name}`;
                await writeBlobInput(runtime,target,await readDiskInput(file));
              }
            }
            await streamOutputFiles(runtime, runId, step.files || [], `checkpoint-${checkpoint++}`, true, transferred);
          }
        } else {
          raw = await runtime.runPythonAsync(`run_phyloprep_finalize(${pythonJson(manifest)}, ${pythonJson(request.options)}, ${JSON.stringify(stage.outputRoot)}, progress_callback=_phyloprep_progress_callback, manifest_only=True)`);
        }
      } finally {
        if (stepsStarted) {
          await runtime.runPythonAsync("_phyloprep_steps.close(); del _phyloprep_steps");
        }
        runtime.globals.delete("_phyloprep_progress_callback");
      }
      const result = JSON.parse(raw);
      if (!result.ok && !result.recoverable) {
        throw new Error(result.user_error || result.error || "The browser PhyloPrep pipeline could not complete.");
      }
      await streamOutputFiles(runtime, runId, result.files || [], undefined, request.incrementalSave, transferred);
      self.postMessage({
        type: "result",
        runId,
        result: {
          ...result,
          engineMetadata: engineManifest,
          partial: Boolean(result.recoverable),
          workspace_storage: diskWorkspace ? diskWorkspace.stats() : {kind:'memfs'},
          provenance: `Browser phyloprep.py adapter · original downstream modules · source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
        },
      });
      return;
    }

    if (request.type === "orthofinder-extract") {
      self.postMessage({ type: "stage", runId, stage: "process", label: "Extracting OrthoFinder orthogroups" });
      const args = orthofinderArguments(runtime, request, stage);
      const raw = await runtime.runPythonAsync(`run_original_module("extract_orthofinder_orthogroups.py", ${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)}, manifest_only=True)`);
      const result = JSON.parse(raw);
      if (!result.ok) throw new Error(result.error || "OrthoFinder orthogroup extraction failed.");
      await streamOutputFiles(runtime, runId, result.files || [], undefined, request.incrementalSave, transferred);
      self.postMessage({
        type: "result",
        runId,
        result: {
          ...result,
          engineMetadata: engineManifest,
          warnings: String(result.stderr || "").split(/\r?\n/).filter(line => line.startsWith("WARNING:")),
          provenance: `Original extract_orthofinder_orthogroups.py · source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
        },
      });
      return;
    }

    self.postMessage({ type: "stage", runId, stage: "process", label: "Running original script locally" });
    const args = originalArguments(request, stage);
    const raw = await runtime.runPythonAsync(`run_original_module(${JSON.stringify(request.script)}, ${JSON.stringify(args)}, ${JSON.stringify(stage.outputRoot)}, manifest_only=True)`);
    const result = JSON.parse(raw);
    if (!result.ok) throw new Error(result.error || "The original script could not complete.");
    await streamOutputFiles(runtime, runId, result.files || [], undefined, request.incrementalSave, transferred);
    self.postMessage({
      type: "result",
      runId,
      result: {
        ...result,
        engineMetadata: engineManifest,
        provenance: `Original source bundle ${engineManifest.bundleSha256.slice(0, 12)}`,
      },
    });
  } catch (error) {
    self.postMessage({ type: "error", runId, message: diskWorkspace?.failure || error.message || String(error) });
  } finally {
    if (runtime && stage?.runRoot) removeTree(runtime, stage.runRoot);
    else if (runtime && diskWorkspace) removeTree(runtime, diskWorkspace.path);
    if (diskWorkspace) {
      try {await diskWorkspace.close();}
      catch(error) {console.warn(`Could not remove temporary disk cache: ${error.message}`);}
    }
  }
});
