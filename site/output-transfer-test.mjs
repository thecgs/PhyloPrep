import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const text = 'A'.repeat(1024 * 1024 - 1) + 'é' + 'B'.repeat(1024 * 1024) + '\n';
const data = new TextEncoder().encode(text);
let handler, open = 0, closed = 0, outstanding = 0, maxOutstanding = 0;
const decoder = new TextDecoder(), parts = [], chunks = [], reused = [];
const context = vm.createContext({ setTimeout, clearTimeout, self: {
  addEventListener(type, callback) { handler = callback; },
  postMessage(message, transfer = []) {
    const received = structuredClone(message, { transfer });
    if (received.type === 'output-chunk') {
      outstanding++; maxOutstanding = Math.max(maxOutstanding, outstanding);
      chunks.push(received.bytes.length);
      parts.push(decoder.decode(received.bytes, { stream: true }));
      queueMicrotask(() => { outstanding--; void handler({ data: { type: 'output-chunk-ack', runId: received.runId, sequence: received.sequence } }); });
    } else if (received.type === 'output-file-reuse') {
      reused.push(received);
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

context.transferred = new Map();
context.files.forEach((file, i) => {file.sha256 = `hash-${i}`;});
await vm.runInContext('streamOutputFiles(runtime, 11, files, "checkpoint", true, transferred)', context);
const readsBeforeReuse = open;
await vm.runInContext('streamOutputFiles(runtime, 11, files, undefined, true, transferred)', context);
assert.equal(open, readsBeforeReuse, 'Unchanged checkpoint outputs must not be read/transferred again');
assert.equal(reused.length, 2);
assert.deepEqual(reused.map(file => file.index), [0, 1]);
context.files[0].sha256 = 'changed';
await vm.runInContext('streamOutputFiles(runtime, 11, files, undefined, true, transferred)', context);
assert.equal(open, readsBeforeReuse + 1, 'Changed output must be transferred again');
console.log('Checkpoint/final output reuse and changed-content retransmission passed.');

// Actual main-thread receiver reuses checkpoint contents and defers MSAP exports.
const app = await readFile('site/app.js', 'utf8');
let receiver;
const published = [];
const receiveContext = vm.createContext({
  TextDecoder, clearInterval, Promise,
  activeWorkerRun: {id:42, streamedOutputs:new Map()}, pyodideWorker:undefined,
  producedFiles:new Map([['qc/pass.tsv', {name:'qc/pass.tsv',content:'gene.fasta\n'}]]),
  elements:{cancelRun:{},runtimeStatus:{}}, setRunStage:()=>{},
  publishProducedFiles: async files => {published.push(...files);},
  Worker: class {
    constructor() {this.handlers={};receiver=this;}
    addEventListener(type, callback) {this.handlers[type]=callback;}
    postMessage() {}
  },
});
vm.runInContext(app.slice(app.indexOf('function createPyodideWorker('),app.indexOf('function runInWorker(')),receiveContext);
vm.runInContext('createPyodideWorker()',receiveContext);
receiver.handlers.message({data:{type:'output-file-reuse',runId:42,index:0,name:'qc/pass.tsv'}});
assert.equal(receiveContext.activeWorkerRun.streamedOutputs.get('final:0').content,'gene.fasta\n');
receiveContext.activeWorkerRun.deferPublish=true;
receiver.handlers.message({data:{type:'output-chunk',runId:42,index:1,name:'gene.prot.aln.fasta',bytes:new TextEncoder().encode('>a\nAA\n')}});
receiver.handlers.message({data:{type:'output-file-end',runId:42,index:1,name:'gene.prot.aln.fasta',size:6,sha256:'b'.repeat(64)}});
await Promise.resolve();
assert.equal(published.length,0,'MSAP outputs must wait for the batch adapter to assign alignments/ paths');
assert.equal(receiveContext.activeWorkerRun.streamedOutputs.get('final:1').content,'>a\nAA\n');
assert.equal(receiveContext.activeWorkerRun.streamedOutputs.get('final:1').size,6);
assert.equal(receiveContext.activeWorkerRun.streamedOutputs.get('final:1').sha256,'b'.repeat(64));
let finished;
receiveContext.activeWorkerRun.resolve=result=>{finished=result;};
receiver.handlers.message({data:{type:'result',runId:42,result:{files:[
  {name:'qc/pass.tsv',path:'/qc/pass.tsv',size:11,sha256:'a'.repeat(64)},
  {name:'gene.prot.aln.fasta',path:'/gene.prot.aln.fasta',size:6,sha256:'b'.repeat(64)},
]}}});
assert.equal(finished.files[1].sha256,'b'.repeat(64));
assert.equal(finished.files[1].size,6);
console.log('Main-thread checkpoint reuse and deferred MSAP export paths passed.');

// A failed worker is terminated; retry gets a new instance, and stale events are ignored.
let terminated=0, rejected=[], created=0;
receiveContext.Worker=class {
  constructor(){this.handlers={};this.sent=[];receiver=this;created++;}
  addEventListener(type,callback){this.handlers[type]=callback;}
  postMessage(message){this.sent.push(message);}
  terminate(){terminated++;}
};
receiveContext.pyodideWorker=undefined;
vm.runInContext(app.slice(app.indexOf('function runInWorker('),app.indexOf('function cancelActiveRun(')),receiveContext);
receiveContext.nextRunId=100;receiveContext.outputSaveSession=undefined;
let pending=vm.runInContext('runInWorker({script:"test.py"})',receiveContext);
pending.catch(error=>rejected.push(error.message));
const dead=receiver;
dead.handlers.error({message:'worker crashed'});
await pending.catch(()=>{});
assert.equal(receiveContext.pyodideWorker,undefined);
assert.equal(receiveContext.activeWorkerRun,undefined);
assert.equal(terminated,1);
pending=vm.runInContext('runInWorker({script:"retry.py"})',receiveContext);
pending.catch(error=>rejected.push(error.message));
assert.equal(created,2);
dead.handlers.error({message:'late old error'});
assert.ok(receiveContext.activeWorkerRun,'An old worker must not reject a new run');
receiver.handlers.messageerror({});await pending.catch(()=>{});
assert.equal(terminated,2);
assert.equal(rejected.length,2);

pending=vm.runInContext('runInWorker({script:"incomplete.py"})',receiveContext);
pending.catch(()=>{});
receiver.handlers.message({data:{type:'result',runId:103,result:{files:[{name:'missing.fa',path:'/missing.fa'}]}}});
await pending.catch(()=>{});
assert.equal(terminated,3,'An incomplete result must also discard its worker');
assert.equal(receiveContext.activeWorkerRun,undefined);
console.log('Worker crash/decode/incomplete-output cleanup, retry recreation and stale-event isolation passed.');

receiveContext.Worker.prototype.postMessage=()=>{throw new Error('cannot clone request');};
pending=vm.runInContext('runInWorker({script:"clone-failure.py"})',receiveContext);
await assert.rejects(pending,/cannot clone request/);
assert.equal(terminated,4);
assert.equal(receiveContext.pyodideWorker,undefined);
assert.equal(receiveContext.activeWorkerRun,undefined);

assert.equal(receiveContext.activeWorkerRun, undefined);
