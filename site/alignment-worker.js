const AIOLI_URL = "https://biowasm.com/cdn/v3/aioli.js";
let muscleCliPromise;
let aioliLoaded = false;

function localAligner(config) {
  return { ...config, urlPrefix: new URL(`aligners/${config.tool}/${config.version}`, self.location.href).href };
}

function ensureAioliLoaded() {
  if (!aioliLoaded) {
    importScripts(AIOLI_URL);
    if (typeof globalThis.Aioli !== "function") throw new Error("Biowasm Aioli did not initialize.");
    aioliLoaded = true;
  }
}

async function ensureMuscleCli() {
  if (!muscleCliPromise) {
    muscleCliPromise = (async () => {
      ensureAioliLoaded();
      return await new globalThis.Aioli([
        "base/1.0.0",
        localAligner({ tool: "muscle", version: "5.1.0", program: "muscle", reinit: true }),
      ], { printInterleaved: false });
    })();
  }
  return muscleCliPromise;
}

async function createMafftCli() {
  ensureAioliLoaded();
  return await new globalThis.Aioli([
    localAligner({ tool: "mafft", version: "7.520", program: "tbfast", reinit: false }),
    localAligner({ tool: "mafft", version: "7.520", program: "dvtditr", reinit: false }),
  ], { printInterleaved: false });
}

function sequenceStats(text) {
  if(text && typeof text==='object')return text;
  let count=0,longest=0,length=0,lineStart=true,header=false;
  for(const char of String(text || '')){
    if(lineStart && char==='>'){longest=Math.max(longest,length);length=0;count++;header=true;}
    else if(!header && char!=='\n' && char!=='\r' && char!==' ' && char!=='\t')length++;
    if(char==='\n' || char==='\r'){lineStart=true;header=false;}else lineStart=false;
  }
  return {count,longest:Math.max(longest,length)};
}

function memoryFailureMessage(tool, diagnostic, sequenceText) {
  if (!/myalloc(?:64)?\([^\n]*failed|out of memory|cannot enlarge memory|allocation (?:failed|error)|memory access out of bounds|Array buffer allocation failed/i.test(diagnostic)) return null;
  const {count,longest}=sequenceStats(sequenceText);
  const limit = /(?:the )?limit is (\d+) bytes/i.exec(diagnostic);
  const limitNote = limit ? ` The runtime reported a limit of ${(Number(limit[1]) / 1073741824).toFixed(2)} GiB.` : "";
  const alternative = tool === "MUSCLE"
    ? " Alternatively, explicitly select MAFFT to try its browser implementation. The MUSCLE RAM banner is a build-reported value, not a measurement of your computer's RAM."
    : " MAFFT has already failed; selecting MAFFT again will not increase the available memory.";
  return `${tool} could not allocate memory in the browser/WebAssembly runtime (${count} sequences; longest ${longest} residues/bases).` +
    limitNote + " Long sequences can require large alignment matrices even with few taxa. The self-hosted aligner builds permit growth up to 4 GiB per module, but the browser may provide less and the task may require more. " +
    "Run this dataset with the native PhyloPrep CLI." + alternative +
    " The aligner was not switched automatically.\n\nOriginal diagnostic:\n" + diagnostic;
}

function checkExecution(execution, tool, sequenceText) {
  const diagnostic = [execution?.stderr, execution?.stdout].filter(Boolean).join("\n");
  const memoryError = memoryFailureMessage(tool, diagnostic, sequenceText);
  if (memoryError) throw new Error(memoryError);
  // Aioli 3 returns captured streams rather than a process exit-code field.
  if (/---Fatal error---|Allocation error|program exited \(with status: [1-9]\d*\)|Aborted\(|RuntimeError:/i.test(diagnostic)) {
    throw new Error(`${tool} failed. Original diagnostic:\n${diagnostic}`);
  }
}

function validateAlignment(alignment, tool, execution) {
  if (typeof alignment === "string" && alignment.trimStart().startsWith(">")) return alignment;
  const stderr = String(execution?.stderr || "").trim();
  const stdout = String(execution?.stdout || "").trim();
  throw new Error(stderr || stdout || `${tool} did not produce a readable FASTA alignment.`);
}

self.addEventListener("message", async (event) => {
  const request = event.data;
  request.sequenceStats=sequenceStats(request.sequenceText);
  if (request.type !== "align") return;
  const { runId, tool, seqtype } = request;
  let cli;
  let inputPath;
  let outputPath;
  try {
    cli = tool === "mafft" ? await createMafftCli() : await ensureMuscleCli();
    const token = `phyloprep-${Date.now()}-${runId}`;
    const mounted = await cli.mount([{ name: `${token}.fasta`, data: request.sequenceText }]);
    delete request.sequenceText;
    inputPath = mounted?.[0];
    if (!inputPath) throw new Error("Biowasm could not mount the alignment input.");

    if (tool === "mafft") {
      // Biowasm exposes MAFFT's compiled core programs (tbfast + dvtditr),
      // not the upstream shell wrapper named "mafft". Follow the supported
      // two-stage Biowasm recipe and force the same sequence type selected in
      // PhyloPrep (-P protein/codon, -D nucleotide).
      const seqFlag = seqtype === "nucl" ? "-D" : "-P";
      const prePath = "/shared/data/pre";
      const tbfast = await cli.exec(
        "tbfast _ -u 0.0 -l 2.7 -C 0 -b 62 -g 0.0 -f -2.00 -Q 100.0 -h 0.0 " +
        "-O -6.00 -E -0.000 -N -Z _ -+ 16 -W 0.00001 -V -1.53 -s 0.0 -O -C 0 " +
        "-b 62 -f -1.53 -Q 100.0 -h 0.000 -l 2.7 -X 0.1 " + seqFlag + " -i " + inputPath
      );
      checkExecution(tbfast, "MAFFT", request.sequenceStats);
      // Do not refine a partial or missing pre file after a failed first stage.
      validateAlignment(await cli.cat(prePath), "MAFFT", tbfast);
      const refine = await cli.exec(
        "dvtditr -W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F -l 2.7 -z 50 -b 62 " +
        "-f -1.53 -Q 100.0 -h 0.000 -I 16 -X 0.1 " + seqFlag +
        " -p BAATARI2 -K 0 -i " + prePath
      );
      checkExecution(refine, "MAFFT", request.sequenceStats);
      let alignment = await cli.cat(prePath);
      alignment = validateAlignment(alignment, "MAFFT", {
        stdout: refine?.stdout || tbfast?.stdout,
        stderr: [tbfast?.stderr, refine?.stderr].filter(Boolean).join("\n"),
      });
      self.postMessage({ type: "result", runId, alignment });
      return;
    }

    if (tool === "muscle") {
      outputPath = `/shared/data/${token}.aligned.fasta`;
      const execution = await cli.exec("muscle -align " + inputPath + " -output " + outputPath + " -threads 1");
      checkExecution(execution, "MUSCLE", request.sequenceStats);
      let alignment = await cli.cat(outputPath);
      if (typeof alignment !== "string") {
        const blob = await cli.downloadBlob(outputPath);
        if (blob instanceof Blob) alignment = await blob.text();
      }
      alignment = validateAlignment(alignment, "MUSCLE", execution);
      self.postMessage({ type: "result", runId, alignment });
      return;
    }

    throw new Error(`Unsupported browser aligner: ${tool}`);
  } catch (error) {
    const diagnostic = error.message || String(error);
    const message = diagnostic.includes("Original diagnostic:") ? diagnostic
      : memoryFailureMessage(tool === "muscle" ? "MUSCLE" : "MAFFT", diagnostic, request.sequenceStats) || diagnostic;
    self.postMessage({ type: "error", runId, message });
    if (tool === "muscle") {
      try { await cli?.close(); } catch {}
      muscleCliPromise = undefined;
    }
  } finally {
    if (cli) {
      try { if (outputPath) await cli.fs.unlink(outputPath); } catch {}
      try { if (tool === "mafft") await cli.fs.unlink("/shared/data/pre"); } catch {}
      try { if (inputPath) await cli.fs.unlink(inputPath); } catch {}
      if (tool === "mafft") {
        try { await cli.close(); } catch {}
      }
    }
  }
});
