const AIOLI_URL = "https://biowasm.com/cdn/v3/aioli.js";
let muscleCliPromise;
let aioliLoaded = false;

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
        { tool: "muscle", version: "5.1.0", program: "muscle", reinit: true },
      ], { printInterleaved: false });
    })();
  }
  return muscleCliPromise;
}

async function createMafftCli() {
  ensureAioliLoaded();
  return await new globalThis.Aioli([
    { tool: "mafft", version: "7.520", program: "tbfast", reinit: false },
    { tool: "mafft", version: "7.520", program: "dvtditr", reinit: false },
  ], { printInterleaved: false });
}

function validateAlignment(alignment, tool, execution) {
  if (typeof alignment === "string" && alignment.trimStart().startsWith(">")) return alignment;
  const stderr = String(execution?.stderr || "").trim();
  const stdout = String(execution?.stdout || "").trim();
  throw new Error(stderr || stdout || `${tool} did not produce a readable FASTA alignment.`);
}

self.addEventListener("message", async (event) => {
  const request = event.data;
  if (request.type !== "align") return;
  const { runId, tool, seqtype } = request;
  let cli;
  let inputPath;
  let outputPath;
  try {
    cli = tool === "mafft" ? await createMafftCli() : await ensureMuscleCli();
    const token = `phyloprep-${Date.now()}-${runId}`;
    const mounted = await cli.mount([{ name: `${token}.fasta`, data: request.sequenceText }]);
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
      const refine = await cli.exec(
        "dvtditr -W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F -l 2.7 -z 50 -b 62 " +
        "-f -1.53 -Q 100.0 -h 0.000 -I 16 -X 0.1 " + seqFlag +
        " -p BAATARI2 -K 0 -i " + prePath
      );
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
    self.postMessage({ type: "error", runId, message: error.message || String(error) });
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
