import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import vm from "node:vm";

const source = await readFile("site/app.js", "utf8");
const helpers = source.slice(source.indexOf("function scriptMetadata("), source.indexOf("function refreshVersionDisplay("));
const context = vm.createContext({ versionManifest: undefined, versionLoadFailed: false });
vm.runInContext(helpers, context);
assert.equal(vm.runInContext('scriptVersionLabel("tool.py")', context), "Loading version…");
context.versionManifest = { files: [{ name: "tool.py", version: "v9.8.7" }], aligners: { muscle: { packageVersion: "42.0" } } };
assert.equal(vm.runInContext('scriptVersionLabel("tool.py")', context), "v9.8.7");
assert.equal(vm.runInContext('scriptVersionLabel("unknown.py")', context), "Version not declared");
assert.equal(vm.runInContext('externalAlignerLabel("muscle")', context), "MUSCLE (Biowasm package 42.0)");
context.versionManifest = undefined;
context.versionLoadFailed = true;
assert.equal(vm.runInContext('scriptVersionLabel("tool.py")', context), "Version unavailable");
const manifest = JSON.parse(await readFile("site/engine/manifest.json", "utf8"));
for (const match of source.matchAll(/script: "([^"]+\.py)", title:/g)) {
  assert.ok(manifest.files.find(file => file.name === match[1])?.version, `Missing version: ${match[1]}`);
}
const worker = await readFile("site/alignment-worker.js", "utf8");
for (const match of worker.matchAll(/tool:\s*"(mafft|muscle)",\s*version:\s*"([^"]+)"/g)) {
  assert.equal(manifest.aligners[match[1]].packageVersion, match[2]);
}
console.log("Dynamic script and Biowasm package version checks passed.");
