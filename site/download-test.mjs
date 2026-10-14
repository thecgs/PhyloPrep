import { readFile } from 'node:fs/promises';
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const clicked = [], created = [];
let attached = false, builds = 0;
const context = vm.createContext({ TextEncoder, URL: {
  createObjectURL(blob) { const url = URL.createObjectURL(blob); created.push(url); return url; },
}, document: {
  body: { append() { attached = true; } },
  createElement() { return { click() { assert.equal(attached, true); clicked.push(this.href); }, remove() { attached = false; } }; },
} });
vm.runInContext(await readFile('site/browser-downloads.js', 'utf8'), context);
const { requestDownload, zipParts, writeZip } = context.PhyloPrepDownloads;
const owner = {};
const makeBlob = () => { builds++; return new Blob(['complete output']); };
requestDownload(makeBlob, 'result.txt', owner);
requestDownload(makeBlob, 'result.txt', owner);
assert.equal(builds, 1); assert.equal(clicked[0], clicked[1]);
await new Promise(resolve => setTimeout(resolve, 1200));
assert.equal(await (await fetch(clicked[0])).text(), 'complete output');
for (const url of created) URL.revokeObjectURL(url);
const files = [{ name: 'matrices/测试.fasta', content: '>a\nACGT\n' }, { name: 'empty.txt', content: '' }];
const parts = [];
let closed = false;
await writeZip(files, { async write(bytes) { parts.push(bytes); }, async close() { closed = true; }, async abort() { assert.fail('Unexpected archive abort'); } });
assert.equal(closed, true);
const zip = Buffer.concat(parts.map(part => Buffer.from(part)));
const result = spawnSync('python3', ['-c', 'import io,sys,zipfile; z=zipfile.ZipFile(io.BytesIO(sys.stdin.buffer.read())); assert z.testzip() is None; assert z.read("matrices/测试.fasta")==b">a\\nACGT\\n"; assert z.read("empty.txt")==b""'], { input: zip });
assert.equal(result.status, 0, result.stderr?.toString());
assert.throws(() => Array.from(zipParts(new Array(0xffff))), /file-count limit/);
assert.throws(() => Array.from(zipParts([{ name: 'a'.repeat(65536), content: '' }])), /ZIP32 limits/);
let aborted = false;
await assert.rejects(writeZip(files, { async write() { throw new Error('disk full'); }, async close() { assert.fail('Should not close'); }, async abort() { aborted = true; } }), /disk full/);
assert.equal(aborted, true);
console.log('Download lifetime/retry, ZIP round-trip, and failed-write cleanup checks passed.');
