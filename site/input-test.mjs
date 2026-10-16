import { readFile } from "node:fs/promises";
import assert from "node:assert/strict";
import vm from "node:vm";
const context = vm.createContext({});
vm.runInContext(await readFile('site/browser-inputs.js', 'utf8'), context);
const { readBounded, resolveFiles, readText, pathlistFolders } = context.PhyloPrepInputs;
const file = (name, content, path = '') => ({ name, webkitRelativePath: path, text: async () => content });
const list = file('orthogroups.pathlist', '# ignored\n/old/a.fa\nb.fa\n; ignored\n');
const a = file('a.fa', '>a\nAC\n'), b = file('b.fa', '>a\nTT\n'), unrelated = file('extra.fa', '>a\nGG\n');
assert.deepEqual(Array.from(await resolveFiles([list, a, b, unrelated])).map(f => f.name), ['a.fa', 'b.fa']);
assert.deepEqual(Array.from(await resolveFiles([list], { sibling: async entry => entry.endsWith('a.fa') ? a : b })).map(f => f.name), ['a.fa', 'b.fa']);
await assert.rejects(resolveFiles([list]), /Missing FASTA.*a.fa.*input folder/);
await assert.rejects(resolveFiles([list, a, file('a.fa', '>b\nAC\n')]), /Ambiguous FASTA/);
await assert.rejects(readText({ name: 'unreadable.fa', text: async () => { throw new Error('permission revoked'); } }), /unreadable.fa.*Reselect/);
const rootFile = file('a.fa', '>a\nAC\n', 'folder/a.fa');
const nestedFile = file('a.fa', '>a\nTT\n', 'folder/sub/a.fa');
assert.equal((await resolveFiles([file('genes.pathlist', '/old/a.fa'), rootFile, nestedFile]))[0], rootFile);
assert.equal((await resolveFiles([file('genes.pathlist', 'sub/a.fa'), rootFile, nestedFile]))[0], nestedFile);
let active = 0, maxActive = 0;
const many = Array.from({ length: 7644 }, (_, i) => i);
const result = await readBounded(many, async i => {
  active++; maxActive = Math.max(maxActive, active);
  await Promise.resolve();
  active--; return i;
});
assert.equal(maxActive, 4);
assert.equal(result.length, 7644);
assert.equal(result[7643], 7643);
console.log('Browser pathlist and 7,644-file bounded-read checks passed.');

const parsed = pathlistFolders('# comment\n/home/data/a.fa\n/home/data/b.fa\nC:\\input\\c.fa\nd.fa\n; comment\n');
assert.equal(parsed.count, 4);
assert.deepEqual(Array.from(parsed.folders), ['/home/data', 'C:/input', 'beside the pathlist']);


const readEvents = [];
const reported = await readBounded([a, b], async f => f.text(), 2, event => readEvents.push(event));
assert.equal(reported.length, 2);
assert.deepEqual(readEvents.map(e => e.completed), [0, 1, 2]);
assert.ok(readEvents.every(e => e.total === 2));
const resolveEvents = [];
await resolveFiles([list, a, b], { onProgress: event => resolveEvents.push(event) });
assert.deepEqual(resolveEvents.map(e => e.completed), [0, 1, 2]);
assert.equal(resolveEvents.at(-1).name, 'orthogroups.pathlist');
// Do not keep reporting reads that complete after another read failed.
const failureEvents = [];
let finishLate;
const late = new Promise(resolve => { finishLate = resolve; });
await assert.rejects(readBounded([a, b], async f => {
  if (f === a) throw new Error('read failed');
  await late;
  return 'late';
}, 2, event => failureEvents.push(event)), /read failed/);
finishLate();
await Promise.resolve(); await Promise.resolve();
assert.equal(failureEvents.length, 1);
