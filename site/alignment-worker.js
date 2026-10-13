const AIOLI_URL = "https://biowasm.com/cdn/v3/aioli.js";
let cliPromise;

async function ensureCli() {
  if (!cliPromise) {
    cliPromise = (async () => {
      importScripts(AIOLI_URL);
      if (typeof globalThis.Aioli !== "function") throw new Error("Biowasm Aioli did not initialize.");
      return await new globalThis.Aioli(["mafft/7.520", "muscle/5.1.0"], { printInterleaved: false });
    })();
  }
  return cliPromise;
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
    cli = await ensureCli();
    const token = `phyloprep-${Date.now()}-${runId}`;
    const mounted = await cli.mount([{ name: `${token}.fasta`, data: request.sequenceText }]);
    inputPath = mounted?.[0];
    if (!inputPath) throw new Error("Biowasm could not mount the alignment input.");

    if (tool === "mafft") {
      const mode = seqtype === "nucl" ? "--nuc" : "--amino";
      const execution = await cli.exec(`mafft ${mode} --thread 1 --quiet --auto ${inputPath}`);
      const alignment = validateAlignment(execution?.stdout, "MAFFT", execution);
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
      try { if (inputPath) await cli.fs.unlink(inputPath); } catch {}
    }
  }
});
