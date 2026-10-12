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
  let url;
  try {
    cli = await ensureCli();
    const token = `phyloprep-${Date.now()}-${runId}`;
    const mounted = await cli.mount([{ name: `${token}.fasta`, data: request.sequenceText }]);
    inputPath = mounted?.[0];
    if (!inputPath) throw new Error("Biowasm could not mount the MUSCLE input.");
    outputPath = `/shared/data/${token}.aligned.fasta`;
    const execution = await cli.exec("muscle", ["-align", inputPath, "-output", outputPath, "-threads", "1"]);
    url = await cli.download(outputPath);
    const response = await fetch(url);
    if (!response.ok) throw new Error("Biowasm could not read the MUSCLE output.");
    const alignment = await response.text();
    if (!alignment.trimStart().startsWith(">")) {
      const detail = execution?.stderr || execution?.stdout || "MUSCLE did not produce a FASTA alignment.";
      throw new Error(String(detail).trim());
    }
    self.postMessage({ type: "result", runId, alignment });
  } catch (error) {
    self.postMessage({ type: "error", runId, message: error.message || String(error) });
  } finally {
    if (url) URL.revokeObjectURL(url);
    if (cli) {
      try { if (outputPath) await cli.fs.unlink(outputPath); } catch {}
      try { if (inputPath) await cli.fs.unlink(inputPath); } catch {}
    }
  }
});
