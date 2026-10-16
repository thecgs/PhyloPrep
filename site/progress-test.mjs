import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const app = await readFile('site/app.js', 'utf8');
const source = (start, end) => app.slice(app.indexOf(start), app.indexOf(end, app.indexOf(start)));
const messages = [];
const timers = new Set();
let clock = 0;
const elements = {
  run: { disabled: true }, cancelRun: {}, runtimeStatus: {}, runStageLabel: {},
  runStages: { querySelectorAll: () => [] }, outputName: { value: '' }, useNextStep: {},
  log: { set textContent(value) { messages.push(value); }, get textContent() { return messages.at(-1); } },
};
const context = vm.createContext({
  elements, performance: { now: () => clock += 101 }, Map,
  setInterval: callback => { timers.add(callback); return callback; },
  clearInterval: callback => timers.delete(callback),
});
vm.runInContext(source('function setLog(', 'function selectModule('), context);
vm.runInContext(source('function setRunStage(', 'function renderEmptyOutputFiles('), context);
vm.runInContext(source('function showLocusProgress(', 'async function runMsapBatchBrowserWorkflow('), context);
vm.runInContext('showLocusProgress(2, 10, 0, "Aligning gene.fa")', context);
assert.match(elements.log.textContent, /2\/10 loci complete.*Aligning gene.fa/);
vm.runInContext('setLog("specific failure", "error"); setRunStage(2, "Failed", true)', context);
assert.equal(elements.log.textContent, 'specific failure');

let saved = 0;
Object.assign(context, {
  activeModule: () => ({script:'tool.py', inputFields:[], defaultOutput:'results'}),
  browserInputFiles: () => [], inputControls: () => [], rawParameterValues: () => ({}),
  setRunControls: running => { elements.run.disabled = running; },
  preflightGeneration: 0, resolvedFolderFiles: new Map(), generatedFiles: [],
  runInWorker: async () => ({files:[{name:'a.fa',content:'A'},{name:'b.fa',content:'B'}], summary:'Completed'}),
  buildRunManifest: async () => '{}',
  saveGeneratedFile: async () => {
    assert.match(elements.log.textContent, new RegExp(`Saving output: ${saved}/3 files`));
    saved++;
    return 'saved';
  },
  renderGeneratedFiles: () => {}, transferTarget: () => undefined,
});
vm.runInContext(source('async function runBrowserModule(', 'function refreshCatalogSearch('), context);
await vm.runInContext('runBrowserModule()', context);
assert.equal(saved, 3);
assert.ok(messages.some(m => /Saving output: 3\/3 files/.test(m)));
assert.equal(elements.log.textContent, 'Completed');
assert.equal(timers.size, 0, 'Input heartbeat must stop before execution/completion');
context.runInWorker = async () => { throw new Error('worker failure'); };
await vm.runInContext('runBrowserModule()', context);
assert.equal(elements.log.textContent, 'worker failure');
assert.equal(timers.size, 0, 'Input heartbeat must stop on errors');

let worker;
context.Worker = class {
  constructor() { this.handlers = {}; worker = this; }
  addEventListener(type, handler) { this.handlers[type] = handler; }
};
context.pyodideWorker = undefined;
context.activeWorkerRun = {id:1};
elements.run.disabled = true;
vm.runInContext(source('function createPyodideWorker(', 'function runInWorker('), context);
vm.runInContext('createPyodideWorker()', context);
worker.handlers.message({data:{type:'stage',runId:1,stage:'runtime',label:'Loading scripts'}});
assert.equal(elements.log.textContent, 'Loading scripts');
worker.handlers.message({data:{type:'progress',runId:1,progress:{completed:1,total:3,unit:'stages',label:'raw: QC',elapsed_seconds:5}}});
assert.match(elements.log.textContent, /1\/3 stages.*raw: QC/);
console.log('Workspace progress, output-save counts, worker messages and heartbeat cleanup checks passed.');
