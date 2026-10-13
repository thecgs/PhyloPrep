const AIOLI_URL = "https://biowasm.com/cdn/v3/aioli.js";
let cliPromise;

async function ensureCli() {
  if (!cliPromise) {
    cliPromise = (async () => {
      importScripts(AIOLI_URL);
      if (typeof globalThis.Aioli !== "function") throw new Error("Biowasm Aioli did not initialize.");
      return await new globalThis.Aioli(["muscle/5.1.0"], { printInterleaved: false });
    })();
  }
  return cliPromise;
}

self.addEventListener("message", async (event) => {
  const request = event.data;
  if (request.type !== "align") return;
  const { runId } = request;
  let cli;
  let inputPath;
  let outputPath;
  try {
    cli = await ensureCli();
    const token = `phyloprep-${Date.now()}-${runId}`;
    const mounted = await cli.mount([{ name: `${token}.fasta`, data: request.sequenceText }]);
    inputPath = mounted?.[0];
    if (!inputPath) throw new Error("Biowasm could not mount the MUSCLE input.");
    outputPath = `/shared/data/${token}.aligned.fasta`;
    const execution = await cli.exec("muscle", ["-align", inputPath, "-output", outputPath, "-threads", "1"]);

    // Aioli.download() creates a blob URL inside Aioli's own internal worker.
    // This file is already running in a second Web Worker, so reading that
    // worker-scoped URL with fetch() is not reliable. MUSCLE output is text;
    // read it directly from Aioli's shared Emscripten filesystem instead.
    let alignment = await cli.cat(outputPath);
    if (typeof alignment !== "string") {
      const blob = await cli.downloadBlob(outputPath);
      if (blob instanceof Blob) alignment = await blob.text();
    }
    if (typeof alignment !== "string" || !alignment.trimStart().startsWith(">")) {
      const stderr = String(execution?.stderr || "").trim();
      const stdout = String(execution?.stdout || "").trim();
      const detail = stderr || stdout || `MUSCLE did not create a readable FASTA file at ${outputPath}.`;
      throw new Error(detail);
    }
    self.postMessage({ type: "result", runId, alignment });
  } catch (error) {
    self.postMessage({ type: "error", runId, message: error.message || String(error) });
  } finally {
    if (cli) {
      try { if (outputPath) await cli.fs.unlink(outputPath); } catch {}
      try { if (inputPath) await cli.fs.unlink(inputPath); } catch {}
    }
  }
});
