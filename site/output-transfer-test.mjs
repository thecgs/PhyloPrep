import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const text = 'A'.repeat(1024 * 1024 - 1) + 'é' + 'B'.repeat(1024 * 1024) + '\n';
const data = new TextEncoder().encode(text);
let handler, open = 0, closed = 0, outstanding = 0, maxOutstanding = 0;
const decoder = new TextDecoder(), parts = [], chunks = [];
const context = vm.createContext({ setTimeout, clearTimeout, self: {
  addEventListener(type, callback) { handler = callback; },
  postMessage(message, transfer = []) {
    const received = structuredClone(message, { transfer });
    if (received.type === 'output-chunk') {
      outstanding++; maxOutstanding = Math.max(maxOutstanding, outstanding);
      chunks.push(received.bytes.length);
      parts.push(decoder.decode(received.bytes, { stream: true }));
      queueMicrotask(() => { outstanding--; void handler({ data: { type: 'output-chunk-ack', runId: received.runId, sequence: received.sequence } }); });
    } else if (received.type === 'output-file-end') {
      parts.push(decoder.decode());
      if (received.sequence !== undefined) queueMicrotask(() => void handler({data:{type:'output-chunk-ack',runId:received.runId,sequence:received.sequence}}));
    }
  },
} });
vm.runInContext(await readFile('site/pyodide-worker.js', 'utf8'), context);
context.runtime = { FS: {
  open() { open++; return {}; },
  read(stream, bytes, offset, length, position) { bytes.set(data.subarray(position, position + length)); return length; },
  close() { closed++; },
} };
context.files = [{ name: 'matrix.fasta', path: '/matrix.fasta', size: data.length }, { name: 'empty.txt', path: '/empty.txt', size: 0 }];
await vm.runInContext('streamOutputFiles(runtime, 9, files)', context);
assert.equal(parts.join(''), text);
assert.equal(maxOutstanding, 1);
assert.ok(chunks.every(size => size <= 1024 * 1024));
assert.equal(open, 2); assert.equal(closed, 2);
assert.equal(vm.runInContext('outputAcknowledgements.size', context), 0);
await vm.runInContext('streamOutputFiles(runtime, 10, files, "checkpoint-0", true)', context);
assert.equal(vm.runInContext('outputAcknowledgements.size', context), 0);
assert.equal(open, 4); assert.equal(closed, 4);
console.log('Acknowledged output-chunk transfer and UTF-8 boundary checks passed.');
