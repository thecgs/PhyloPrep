import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const scope = vm.createContext({});
vm.runInContext(await readFile('site/browser-save.js', 'utf8'), scope);
const { createWriter } = scope.PhyloPrepSave;
function fakeFolder() {
  const disk = new Map(), calls = [], busy = new Set();
  let active = 0, peak = 0, failName;
  const directory = path => ({
    async getDirectoryHandle(name) { calls.push(`dir:${path}/${name}`); await Promise.resolve(); return directory(`${path}/${name}`); },
    async getFileHandle(name) {
      const key = `${path}/${name}`; calls.push(`file:${key}`);
      return { async createWritable() {
        assert.ok(!busy.has(key), `Concurrent overwrite of ${key}`); busy.add(key);
        active++; peak = Math.max(peak, active); let text;
        return {
          async write(value) { text = value; await Promise.resolve(); },
          async close() {
            if (failName === name) throw new Error('disk full');
            disk.set(key, text); calls.push(`closed:${key}`); busy.delete(key); active--;
          },
          async abort() { calls.push(`abort:${key}`); busy.delete(key); active--; },
        };
      } };
    },
  });
  return { root: directory(''), disk, calls, get peak() { return peak; }, fail(name) { failName = name; } };
}
const folder = fakeFolder(), writer = createWriter(folder.root);
const files = Array.from({length:100}, (_,i) => ({name:`alignments/codon/gene${i}.aln.fasta`,content:`seq${i}`}));
for (const file of files) writer.save(file);
assert.notEqual(writer.stats().saved, 100, 'Must wait for close before reporting completion');
await writer.settle();
assert.equal(writer.stats().saved, 100);
assert.equal(folder.peak, 4);
assert.equal(folder.calls.filter(c => c.startsWith('dir:')).length, 2, 'Directory creation must be shared across concurrent saves');
for (const file of files) writer.save({...file});
await writer.settle();
assert.equal(folder.calls.filter(c => c.startsWith('closed:')).length, 100, 'Identical results must not be saved twice');
writer.save({name:'matrices/result.fa',content:'first'});
writer.save({name:'matrices/result.fa',content:'second'});
await writer.settle();
assert.equal(folder.disk.get('/matrices/result.fa'), 'second');
assert.ok(writer.saved({name:'matrices/result.fa',content:'second'}));
assert.ok(!writer.saved({name:'matrices/result.fa',content:'first'}));

const broken = fakeFolder(); broken.fail('gene0.aln.fasta');
const retry = createWriter(broken.root);
files.forEach(file => retry.save(file));
await retry.settle();
assert.ok(retry.error);
assert.ok(retry.stats().failed > 0);
assert.ok(retry.stats().pending > 0);
assert.ok(broken.calls.some(c => c.startsWith('abort:')));
const successes = broken.calls.filter(c => c.startsWith('closed:'));
broken.fail(undefined);
await retry.retry();
files.forEach(file => retry.save(file));
await retry.settle();
assert.equal(retry.stats().saved, 100);
for (const call of successes) assert.equal(broken.calls.filter(c => c === call).length, 1);

// Integration: canonical pipeline paths are available before the next locus starts.
const app = await readFile('site/app.js','utf8');
const slice = (start,end) => app.slice(app.indexOf(start),app.indexOf(end,app.indexOf(start)));
const published = new Map(); let processed = 0;
const context = vm.createContext({
  performance:{now:()=>0}, setInterval:()=>1, clearInterval:()=>{},
  showLocusProgress:()=>{}, requireUniqueStems:()=>{}, selectedInputPathList:()=>'',
  externalAlignerLabel:()=> 'MAFFT', batchPathLists:()=>[], mergeTsvContents:()=>'',
  publishProducedFiles:async files => {for (const file of files) published.set(file.name,file.content);},
  runMsapBrowserWorkflow:async request => {
    if (processed) assert.equal(published.get('alignments/first.codon.aln.fasta'),'aligned');
    processed++;
    return {files:[{name:`${request.files[0].name.split('.')[0]}.codon.aln.fasta`,content:'aligned'}]};
  },
});
vm.runInContext(slice('async function runMsapBatchBrowserWorkflow(', 'async function runPhyloPrepBrowserWorkflow('),context);
context.request={pipelineOutput:true,options:{},files:[{inputKey:'inputs',name:'first.fa',content:'>a\nATG'},{inputKey:'inputs',name:'second.fa',content:'>a\nATG'}]};
const batch=await vm.runInContext('runMsapBatchBrowserWorkflow(request)',context);
assert.ok(batch.files.some(file => file.name === 'first.codon.aln.fasta'), 'Internal batch names must stay unchanged');
assert.ok(published.has('alignments/second.codon.aln.fasta'));
console.log('Four-writer bound, directory caching, deduplication, overwrite serialization, failure retry and per-locus export checks passed.');

// Exercise the actual app saving helpers, completion warnings and manual retry.
let logs = [];
const ui = {
  run:{disabled:false},cancelRun:{},runtimeStatus:{},saveProgress:{},outputName:{value:''},useNextStep:{},downloadAll:{},
};
const exported = [{name:'a.fa',content:'A'},{name:'b.fa',content:'B'}];
const integration = vm.createContext({
  elements:ui,performance:{now:()=>Date.now()},PhyloPrepSave:scope.PhyloPrepSave,
  outputDirectoryHandle:undefined,outputSaveSession:undefined,producedFiles:new Map(),generatedFiles:[],
  resolvedFolderFiles:new Map(),preflightGeneration:0,
  setInterval:()=>1,clearInterval:()=>{},
  activeModule:()=>({script:'tool.py',inputFields:[],defaultOutput:'results'}),
  browserInputFiles:()=>[],inputControls:()=>[],rawParameterValues:()=>({}),
  setRunControls:running=>{ui.run.disabled=running;},
  setRunStage:()=>{},setLog:(text,status)=>logs.push({text,status}),
  runInWorker:async()=>({files:exported,summary:'Completed'}),buildRunManifest:async()=>'{}',
  renderGeneratedFiles:()=>{},transferTarget:()=>undefined,
});
vm.runInContext(slice('function createOutputSaveSession(', 'function renderGeneratedFiles('),integration);
vm.runInContext(slice('async function runBrowserModule(', 'function refreshCatalogSearch('),integration);
vm.runInContext(slice('async function downloadAllOutputs(', 'function createAlignmentWorker('),integration);
await vm.runInContext('runBrowserModule()',integration);
assert.equal(integration.generatedFiles.length,3);
assert.ok(integration.generatedFiles.every(file=>file.destination==='browser workspace'));
assert.equal(integration.outputSaveSession,undefined, 'No folder must mean no automatic writes/downloads');

const local = fakeFolder(); integration.outputDirectoryHandle=local.root;
await vm.runInContext('runBrowserModule()',integration);
assert.equal(local.disk.size,3);
assert.equal(logs.at(-1).status,'success');
const writeCount=local.calls.filter(c=>c.startsWith('closed:')).length;
await vm.runInContext('downloadAllOutputs()',integration);
assert.equal(local.calls.filter(c=>c.startsWith('closed:')).length,writeCount);

const full = fakeFolder();full.fail('a.fa');integration.outputDirectoryHandle=full.root;
await vm.runInContext('runBrowserModule()',integration);
assert.equal(logs.at(-1).status,'warning');
assert.match(logs.at(-1).text,/local saving stopped/);
assert.equal(integration.generatedFiles.length,3);
full.fail(undefined);
await vm.runInContext('downloadAllOutputs()',integration);
assert.equal(full.disk.size,3);
assert.ok(integration.generatedFiles.every(file=>file.destination==='selected local folder'));
const other=fakeFolder();integration.outputDirectoryHandle=other.root;
await vm.runInContext('downloadAllOutputs()',integration);
assert.equal(other.disk.size,3,'A new destination needs its own save session');

integration.runInWorker=async()=>{
  await integration.publishProducedFiles([{name:'alignments/partial.aln.fasta',content:'ATG'}]);
  throw new Error('analysis interrupted');
};
await vm.runInContext('runBrowserModule()',integration);
assert.equal(logs.at(-1).text,'analysis interrupted');
assert.equal(integration.generatedFiles.length,1,'Produced files must survive analysis errors');
assert.equal(other.disk.get('/alignments/partial.aln.fasta'),'ATG');
console.log('No-folder behavior, full workflow save completion, paused-save warning, manual retry, destination changes and partial-result retention passed.');
