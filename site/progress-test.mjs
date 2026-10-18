import { PhyloPrepSave } from './test-output-store.mjs';
import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const app = await readFile('site/app.js', 'utf8');
const source = (start, end) => app.slice(app.indexOf(start), app.indexOf(end, app.indexOf(start)));
const messages = [];
const timers = new Set();
let clock = 0;
const stageClasses = [1,2,3].map(stage => ({dataset:{stage:String(stage)}, classes:new Set(),
  get classList(){return {add:(...names)=>names.forEach(name=>this.classes.add(name)),
    remove:(...names)=>names.forEach(name=>this.classes.delete(name))};}}));
const elements = {
  command:{textContent:"test command"},run: { disabled: true }, cancelRun: {}, runtimeStatus: {}, runStageLabel: {},
  runStages: { querySelectorAll: () => stageClasses }, outputName: { value: '' }, useNextStep: {},
  log: { set textContent(value) { messages.push(value); }, get textContent() { return messages.at(-1); } },
};
const context = vm.createContext({releaseIdleWorkers:()=>{}, PhyloPrepSave, Blob, TextDecoder,
  elements, performance: { now: () => clock += 101 }, Map,
  setInterval: callback => { timers.add(callback); return callback; },
  clearInterval: callback => timers.delete(callback),
});
vm.runInContext(source('function setLog(', 'function selectModule('), context);
vm.runInContext(source('function setRunStage(', 'function renderEmptyOutputFiles('), context);
vm.runInContext(source('function showLocusProgress(', 'async function runMsapBatchBrowserWorkflow('), context);
vm.runInContext('showLocusProgress(2, 10, 0, "Aligning gene.fa")', context);
assert.match(elements.log.textContent, /2\/10 loci complete.*Aligning gene.fa/);
assert.equal(elements.runStageLabel.textContent, 'Processing locally');
assert.equal(elements.runtimeStatus.textContent, 'Processing locally');
vm.runInContext('showLocusProgress(3, 10, 0, "Aligning next.fa")', context);
assert.match(elements.log.textContent, /3\/10 loci complete.*Aligning next.fa/);
assert.equal(elements.runStageLabel.textContent, 'Processing locally');
assert.equal(elements.runtimeStatus.textContent, 'Processing locally');
vm.runInContext('setLog("specific failure", "error"); setRunStage(2, "Failed", true)', context);
assert.equal(elements.log.textContent, 'specific failure');
assert.equal(elements.runStageLabel.textContent, 'Failed');
vm.runInContext('setRunStage(2, "Run cancelled")', context);
assert.equal(elements.runStageLabel.textContent, 'Run cancelled');

let saved = 0;
Object.assign(context, {
  activeModule: () => ({script:'tool.py', inputFields:[], defaultOutput:'results'}),
  browserInputFiles: () => [], inputControls: () => [], rawParameterValues: () => ({}),
  setRunControls: running => { elements.run.disabled = running; },
  preflightGeneration: 0, resolvedFolderFiles: new Map(), generatedFiles: [],
  runInWorker: async () => ({files:[{name:'a.fa',content:'A'},{name:'b.fa',content:'B'}], summary:'Completed'}),
  buildRunManifest: async () => '{}',
  setConfigurationLocked:()=>{}, localSaveBusy:false, outputDirectoryHandle: undefined, outputSaveSession: undefined, producedFiles: new Map(),
  publishProducedFiles: async files => { saved += files.length; },
  savedDestination: () => 'browser workspace',
  renderGeneratedFiles: () => {}, transferTarget: () => undefined,
});
vm.runInContext(source('async function runBrowserModule(', 'function refreshCatalogSearch('), context);
await vm.runInContext('runBrowserModule()', context);
assert.equal(saved, 3);
assert.ok(messages.includes('Waiting for local file saves to finish'));
assert.equal(elements.log.textContent, 'Completed');
assert.equal(elements.runStageLabel.textContent, 'Completed');
assert.equal(timers.size, 0, 'Input heartbeat must stop before execution/completion');
context.runInWorker = async () => { throw new Error('worker failure'); };
await vm.runInContext('runBrowserModule()', context);
assert.equal(elements.log.textContent, 'worker failure');
assert.equal(timers.size, 0, 'Input heartbeat must stop on errors');

let worker;
context.Worker = class {
  constructor() { this.handlers = {}; worker = this; }
  addEventListener(type, handler) { this.handlers[type] = handler; }
  postMessage() {}
};
context.pyodideWorker = undefined;
context.activeWorkerRun = {id:1,streamedOutputs:new Map()};
elements.run.disabled = true;
vm.runInContext(source('function createPyodideWorker(', 'function runInWorker('), context);
vm.runInContext('createPyodideWorker()', context);
worker.handlers.message({data:{type:'stage',runId:1,stage:'runtime',label:'Loading scripts'}});
assert.equal(elements.log.textContent, 'Loading scripts');
worker.handlers.message({data:{type:'progress',runId:1,progress:{completed:1,total:3,unit:'stages',label:'raw: QC',elapsed_seconds:5}}});
assert.match(elements.log.textContent, /1\/3 stages.*raw: QC/);
assert.equal(elements.runStageLabel.textContent, 'Processing locally');
assert.equal(elements.runtimeStatus.textContent, 'Processing locally');
assert.equal(timers.size,1);
worker.handlers.message({data:{type:'stage',runId:1,stage:'process',label:'Checking saved branch outputs'}});
assert.equal(timers.size,0,'A new stage must stop the previous QC heartbeat');
worker.handlers.message({data:{type:'output-chunk',runId:1,transferId:'rename-1',index:0,name:'renamed/codon/raw/gene.aln.fasta',bytes:new Uint8Array([65])}});
assert.equal(elements.runStageLabel.textContent,'Processing locally');
assert.ok(stageClasses[1].classes.has('active'));
assert.ok(!stageClasses[1].classes.has('complete'),'Intermediate saves must not mark analysis complete');
assert.match(elements.log.textContent,/intermediate output \(analysis continues\)/);
context.producedFiles.set('renamed/codon/raw/gene.aln.fasta',{name:'renamed/codon/raw/gene.aln.fasta'});
worker.handlers.message({data:{type:'output-file-reuse',runId:1,transferId:'branch-end',release:true,index:0,fileCount:1,name:'renamed/codon/raw/gene.aln.fasta'}});
assert.match(elements.log.textContent,/Reusing saved output.*checking temporary-file release/);
assert.equal(elements.runStageLabel.textContent,'Processing locally');
context.activeWorkerRun.deferPublish=true;
worker.handlers.message({data:{type:'output-chunk',runId:1,index:1,name:'gene.prot.aln.fasta',bytes:new Uint8Array([65])}});
assert.equal(elements.runStageLabel.textContent,'Processing locally','Per-locus MSAP output remains part of analysis');
context.activeWorkerRun.deferPublish=false;
worker.handlers.message({data:{type:'output-chunk',runId:1,index:2,name:'matrices/codon/raw/supermatrix.fasta',bytes:new Uint8Array([65])}});
assert.equal(elements.runStageLabel.textContent,'Saving output');
assert.ok(stageClasses[1].classes.has('complete'));
assert.match(elements.log.textContent,/Receiving final output/);
console.log('Workspace progress, output-save counts, worker messages and heartbeat cleanup checks passed.');
