import { copyFile, mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { spawnSync } from "node:child_process";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const siteDirectory = dirname(fileURLToPath(import.meta.url));
const sourceDirectory = dirname(siteDirectory);
const engineDirectory = join(siteDirectory, "engine");
const sourceFiles = [
  "MSAP.py",
  "MSAP_batch.py",
  "phyloprep.py",
  "extract_orthofinder_orthogroups.py",
  "oma_to_orthogroups.py",
  "orthogroup_table.py",
  "orthomcl_to_orthogroups.py",
  "proteinortho_to_orthogroups.py",
  "sonicparanoid_to_orthogroups.py",
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
  "msap_progress.py",
  "sequence_audit.py",
  "split_codon_seqence_alignment.py",
];

await mkdir(engineDirectory, { recursive: true });
const files = [];
for (const name of sourceFiles) {
  const source = join(sourceDirectory, name);
  const contents = await readFile(source);
  await copyFile(source, join(engineDirectory, name));
  files.push({ name, sha256: createHash("sha256").update(contents).digest("hex") });
}

const metadataProcess = spawnSync(process.env.PYTHON || "python3", [join(siteDirectory, "script_metadata.py"), ...sourceFiles.map(name => join(engineDirectory, name))], {
  encoding: "utf8", cwd: sourceDirectory,
});
if (metadataProcess.status !== 0) throw new Error(`Could not read script versions: ${metadataProcess.stderr}`);
const metadata = JSON.parse(metadataProcess.stdout);
for (const file of files) file.version = metadata[file.name]?.version || null;
const alignmentWorkerSource = await readFile(join(siteDirectory, "alignment-worker.js"), "utf8");
const aligners = {};
for (const match of alignmentWorkerSource.matchAll(/tool:\s*"(mafft|muscle)",\s*version:\s*"([^"]+)"/g)) {
  aligners[match[1]] = { packageVersion: match[2], provider: "Biowasm" };
}

const bundleSha256 = createHash("sha256")
  .update(files.map(({ name, sha256 }) => `${name}:${sha256}`).join("\n"))
  .digest("hex");
await writeFile(
  join(engineDirectory, "manifest.json"),
  `${JSON.stringify({ files, bundleSha256, aligners }, null, 2)}\n`,
);

console.log(`Bundled ${files.length} original Python sources (${bundleSha256.slice(0, 12)}).`);
