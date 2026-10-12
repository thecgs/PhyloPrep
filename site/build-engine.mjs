import { copyFile, mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const siteDirectory = dirname(fileURLToPath(import.meta.url));
const sourceDirectory = dirname(siteDirectory);
const engineDirectory = join(siteDirectory, "engine");
const sourceFiles = [
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
  "msap_io.py",
  "sequence_audit.py",
];

await mkdir(engineDirectory, { recursive: true });
const files = [];
for (const name of sourceFiles) {
  const source = join(sourceDirectory, name);
  const contents = await readFile(source);
  await copyFile(source, join(engineDirectory, name));
  files.push({ name, sha256: createHash("sha256").update(contents).digest("hex") });
}

const bundleSha256 = createHash("sha256")
  .update(files.map(({ name, sha256 }) => `${name}:${sha256}`).join("\n"))
  .digest("hex");
await writeFile(
  join(engineDirectory, "manifest.json"),
  `${JSON.stringify({ files, bundleSha256 }, null, 2)}\n`,
);

console.log(`Bundled ${files.length} original Python sources (${bundleSha256.slice(0, 12)}).`);
