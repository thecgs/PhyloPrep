import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const source = await readFile('site/alignment-worker.js', 'utf8');
async function scenario(tool, responses, alignment = '>a\nAAAA\n>b\nAAAA\n') {
  let handler, cats = 0, closes = 0;
  const posted = [], calls = [], configs = [];
  class Aioli {
    constructor(config) { configs.push(config); }
    async mount() { return ['/shared/input.fa']; }
    async exec(command) { calls.push(command); return responses.shift() || {}; }
    async cat() { cats++; return alignment; }
    async close() { closes++; }
    fs = { unlink: async () => {} };
  }
  const context = vm.createContext({ Aioli, URL, Blob, importScripts() {}, self: {
    location: { href: 'https://example.org/PhyloPrep/alignment-worker.js' },
    addEventListener(type, callback) { handler = callback; }, postMessage(message) { posted.push(message); },
  } });
  vm.runInContext(source, context);
  await handler({ data: { type: 'align', runId: 1, tool, seqtype: 'prot', sequenceText: '>a\nAAAA\n>b\nAAAAA\n' } });
  return { posted, calls, configs, cats, closes, context };
}
const muscle = await scenario('muscle', [{stderr:'---Fatal error---\nmyalloc64(4, 3.25e+08) failed'}]);
assert.equal(muscle.cats, 0);
assert.equal(muscle.closes, 1);
assert.equal(muscle.posted[0].type, 'error');
assert.match(muscle.posted[0].message, /2 sequences; longest 5/);
assert.match(muscle.posted[0].message, /explicitly select MAFFT/);
assert.match(muscle.posted[0].message, /MUSCLE RAM banner/);
assert.match(muscle.posted[0].message, /4 GiB/);
assert.equal(vm.runInContext('muscleCliPromise', muscle.context), undefined);
assert.equal(vm.runInContext('memoryFailureMessage("MUSCLE", "unrelated input error", "")', muscle.context), null);
const oom = 'Cannot enlarge memory, asked to go up to 2147495936 bytes, but the limit is 2147483648 bytes!\nAllocation error(10148 x 10148 doublemtx)';
const mafft = await scenario('mafft', [{stderr:oom}]);
assert.equal(mafft.calls.length, 1, 'No refinement after tbfast failure');
assert.equal(mafft.cats, 0);
assert.equal(mafft.closes, 1);
assert.equal(mafft.posted[0].type, 'error');
assert.match(mafft.posted[0].message, /2.00 GiB/);
assert.doesNotMatch(mafft.posted[0].message, /explicitly select MAFFT|MUSCLE|100Mb/);
assert.match(mafft.posted[0].message, /native PhyloPrep CLI/);
const refine = await scenario('mafft', [{}, {stderr:'Allocation error(100 x 100 doublemtx)'}]);
assert.equal(refine.calls.length, 2);
assert.equal(refine.posted[0].type, 'error', 'Stale pre must not hide failed refinement');
const nonmemory = await scenario('mafft', [{stderr:'program exited (with status: 1)'}]);
assert.equal(nonmemory.calls.length, 1);
assert.equal(nonmemory.posted[0].type, 'error');
const good = await scenario('mafft', [{}, {}]);
assert.equal(good.posted[0].type, 'result');
for (const result of [muscle, mafft]) {
  for (const config of result.configs.flat().filter(x => typeof x === 'object')) {
    assert.match(config.urlPrefix, /^https:\/\/example.org\/PhyloPrep\/aligners\//);
  }
}
console.log('Aligner-specific memory messages, failure handling, cleanup and self-hosted URLs passed.');
