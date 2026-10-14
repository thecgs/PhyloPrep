import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const source = await readFile('site/alignment-worker.js', 'utf8');
let handler;
const posted = [];
let catCalled = false;
class Aioli {
  async mount() { return ['/shared/input.fa']; }
  async exec() { return { stderr: '---Fatal error---\nmyalloc64(4, 5.27e+07) failed' }; }
  async cat() { catCalled = true; return ''; }
  fs = { unlink: async () => {} };
}
const context = vm.createContext({ Aioli, importScripts() {}, self: {
  addEventListener(type, callback) { handler = callback; },
  postMessage(message) { posted.push(message); },
}, Blob });
vm.runInContext(source, context);
await handler({ data: { type: 'align', runId: 1, tool: 'muscle', seqtype: 'prot', sequenceText: '>a\nAAAA\n>b\nAAAAA\n' } });
assert.equal(catCalled, false);
assert.equal(posted[0].type, 'error');
assert.match(posted[0].message, /browser\/WebAssembly runtime/);
assert.match(posted[0].message, /2 sequences; longest 5/);
assert.match(posted[0].message, /native PhyloPrep CLI/);
assert.match(posted[0].message, /Original diagnostic/);
assert.equal(vm.runInContext('memoryFailureMessage("MUSCLE", "unrelated input error", ">a\\nAA")', context), null);
console.log('Browser alignment memory-error checks passed.');
