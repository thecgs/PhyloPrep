import { PhyloPrepSave } from './test-output-store.mjs';
import { readFile } from 'node:fs/promises';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const scope = vm.createContext({ PhyloPrepSave, Blob, });
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
          async write(value) { text = typeof value === "string" ? value : await value.text(); await Promise.resolve(); },
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
const context = vm.createContext({ PhyloPrepSave, Blob,
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
  command:{textContent:"test command"},run:{disabled:false},cancelRun:{},runtimeStatus:{},saveProgress:{},outputName:{value:''},useNextStep:{},downloadAll:{},chooseFolder:{},downloadSelected:{},
};
const exported = [{name:'a.fa',content:'A'},{name:'b.fa',content:'B'}];
const integration = vm.createContext({ PhyloPrepSave, Blob,
  elements:ui,performance:{now:()=>Date.now()},PhyloPrepSave:scope.PhyloPrepSave,
  setConfigurationLocked:()=>{}, localSaveBusy:false, outputDirectoryHandle:undefined,outputSaveSession:undefined,producedFiles:new Map(),generatedFiles:[],
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

// Portable batch lists and a linear duplicate-stem check use actual app helpers.
const names = vm.createContext({ PhyloPrepSave, Blob, });
vm.runInContext(slice('function fileStem(', 'function mergeTsvContents('), names);
vm.runInContext(slice('function batchPathLists(', 'function showLocusProgress('), names);
names.inputs = Array.from({length:60000}, (_,i) => ({name:`gene${i}.fasta`}));
vm.runInContext('requireUniqueStems(inputs)', names);
names.inputs.push({name:'gene0.fa'});
assert.throws(() => vm.runInContext('requireUniqueStems(inputs)', names), /duplicate stem: gene0/);
const portable = vm.runInContext('batchPathLists([{name:"gene.mafft.codon.aln.fasta"}])', names);
assert.equal(portable[0].content, '../gene.mafft.codon.aln.fasta\n');
let cleared = false;
integration.renderGeneratedFiles = () => { if (!integration.generatedFiles.length) cleared = true; };
integration.runInWorker = async () => {throw new Error('early failure');};
await vm.runInContext('runBrowserModule()', integration);
assert.ok(cleared, 'Starting a failing run must clear the previous result list');
console.log('Portable list naming, 60,000-input stem checking and stale-result clearing passed.');

// A suspended save must keep its original session even if globals change.
let releaseSave;
const saveGate = new Promise(resolve => {releaseSave=resolve;});
const sessionWrites=[];
const sessionA={root:'A',save:file=>sessionWrites.push(`A:${file.name}`),ready:()=>saveGate,settle:async()=>{},saved:()=>true};
const sessionB={root:'B',save:file=>sessionWrites.push(`B:${file.name}`),ready:async()=>{},settle:async()=>{},saved:()=>true};
integration.outputDirectoryHandle='A';integration.outputSaveSession=sessionA;
integration.snapshotFiles=[{name:'old1.fa',content:'old'},{name:'old2.fa',content:'old'}];
const isolated=vm.runInContext('saveOutputCollection(snapshotFiles)',integration);
integration.outputSaveSession=sessionB;integration.outputDirectoryHandle='B';releaseSave();
assert.equal(await isolated,sessionA);
assert.deepEqual(sessionWrites,['A:old1.fa','A:old2.fa']);

// Manual retry locks all conflicting controls until writing settles.
let releaseRetry;
const retryGate=new Promise(resolve=>{releaseRetry=resolve;});
integration.outputDirectoryHandle='locked';
integration.outputSaveSession={root:'locked',retry:()=>retryGate,save:()=>{},ready:async()=>{},settle:async()=>{},saved:()=>true};
integration.generatedFiles=exported;
const retrying=vm.runInContext('downloadAllOutputs()',integration);
assert.ok(integration.localSaveBusy);
for(const control of [ui.run,ui.chooseFolder,ui.downloadAll,ui.downloadSelected]) assert.ok(control.disabled);
let analyses=0;
integration.activeModule=()=>{analyses++;return {script:'tool.py',inputFields:[]};};
await vm.runInContext('runBrowserModule()',integration);
assert.equal(analyses,0,'No analysis may start during manual saving');
vm.runInContext(slice('async function chooseOutputFolder(', 'function cancellationError('),integration);
await vm.runInContext('chooseOutputFolder()',integration); // Guard must return before opening any picker.
releaseRetry();await retrying;
assert.equal(integration.localSaveBusy,false);
for(const control of [ui.run,ui.chooseFolder,ui.downloadAll,ui.downloadSelected]) assert.equal(control.disabled,false);

vm.runInContext(slice('function mergeTsvContents(', 'function selectedInputPathList('),names);
names.reports=[{content:'header\n'+'first\n'.repeat(150000)},{content:'header\nlast\n'}];
const largeReport=vm.runInContext('mergeTsvContents(reports)',names);
assert.equal(largeReport.split('\n').length,150003);
assert.ok(largeReport.startsWith('header\nfirst\n'));
assert.ok(largeReport.endsWith('last\n'));
console.log('Save-session isolation, retry control locking and 150,001-row report merging passed.');

// Actual workflow must capture its options/output/CLI before asynchronous processing.
const rawOptions={genetic_code:1,taxa:['original']};let submitted;
ui.command.textContent='original command -g 1';ui.outputName.value='original-output';
integration.outputDirectoryHandle=undefined;integration.outputSaveSession=undefined;
integration.activeModule=()=>({script:'tool.py',inputFields:[],defaultOutput:'results'});
integration.rawParameterValues=()=>rawOptions;
integration.runInWorker=async request=>{
  submitted=request;ui.command.textContent='changed command -g 2';ui.outputName.value='changed-output';
  rawOptions.genetic_code=2;rawOptions.taxa.push('changed');
  return {files:exported};
};
integration.buildRunManifest=async(module,request)=>{
  assert.equal(request.equivalentCli,'original command -g 1');
  assert.equal(request.options.genetic_code,1);assert.equal(request.options.taxa.length,1);
  assert.equal(request.outputName,'original-output');return '{}';
};
await vm.runInContext('runBrowserModule()',integration);
assert.equal(submitted.options.genetic_code,1);
assert.equal(logs.at(-1).status,'success');
console.log('Workflow start snapshots survive subsequent UI/option changes.');

// Full pipeline validates once and lets Python build its own lists/merged report.
let validations=0, listBuilds=0, reportMerges=0, finalized=0, aligned=0;
let finalRequest;
const adapter=vm.createContext({ PhyloPrepSave, Blob,
  performance:{now:()=>0},setInterval:()=>1,clearInterval:()=>{},
  showLocusProgress:()=>{},setRunStage:()=>{},elements:{runtimeStatus:{}},
  externalAlignerLabel:()=> 'MAFFT',publishProducedFiles:async()=>{},
  runMsapBrowserWorkflow:async request=>{
    aligned++;const stem=request.files[0].name.split('.')[0];
    return {files:[
      {name:`${stem}.mafft.codon.aln.fasta`,content:'>a\nATG\n',size:7,sha256:'a'.repeat(64)},
      {name:`${stem}.sequence_changes.tsv`,content:`header\n${stem}\n`},
      {name:`${stem}.cds.fasta`,content:'unused intermediate'},
    ]};
  },
  runInWorker:async request=>{finalized++;finalRequest=request;return {files:[{name:'matrices/result.fasta',content:'matrix'}]};},
});
vm.runInContext(slice('function fileStem(', 'function showLocusProgress('),adapter);
const validateInputs=adapter.requireUniqueStems,buildLists=adapter.batchPathLists,mergeReports=adapter.mergeTsvContents;
adapter.requireUniqueStems=files=>{validations++;return validateInputs(files);};
adapter.batchPathLists=files=>{listBuilds++;return buildLists(files);};
adapter.mergeTsvContents=files=>{reportMerges++;return mergeReports(files);};
vm.runInContext(slice('async function runMsapBatchBrowserWorkflow(', 'async function sha256Text('),adapter);
adapter.request={options:{},files:[
  {inputKey:'inputs',name:'first.fa',content:'>a\nATG'},
  {inputKey:'inputs',name:'second.fa',content:'>a\nATG'},
  {inputKey:'mapping',name:'map.tsv',content:'a\tTaxon\n'},
]};
const pipeline=await vm.runInContext('runPhyloPrepBrowserWorkflow(request)',adapter);
assert.equal(validations,1);assert.equal(listBuilds,0);assert.equal(reportMerges,0);assert.equal(finalized,1);
assert.equal(finalRequest.files.filter(file=>file.role==='source').length,2);
assert.equal(finalRequest.files.filter(file=>file.role==='mapping').length,1);
assert.ok(finalRequest.files.filter(file=>file.role==='batch-output').every(file=>!file.name.includes('/') && !file.name.endsWith('.cds.fasta')));
assert.equal(finalRequest.files.find(file=>file.name==='first.mafft.codon.aln.fasta').sha256,'a'.repeat(64));
assert.ok(pipeline.files.some(file=>file.name==='inputs.pathlist'));
validations=0;
const standalone=await vm.runInContext('runMsapBatchBrowserWorkflow(request)',adapter);
assert.equal(validations,1);assert.equal(listBuilds,1);assert.equal(reportMerges,1);
assert.ok(standalone.files.some(file=>file.name==='path-lists/all.codon.aln.fasta.pathlist'));
assert.ok(standalone.files.some(file=>file.name==='first.cds.fasta'));
assert.equal(standalone.files.find(file=>file.name==='reports/sequence_changes.tsv').content,'header\nfirst\nsecond\n');
adapter.request.files=[{inputKey:'inputs',name:'first.fa',content:'>a\nATG'},{inputKey:'inputs',name:'first.fasta',content:'>a\nATG'}];
const alignedBeforeFailure=aligned;
await assert.rejects(vm.runInContext('runPhyloPrepBrowserWorkflow(request)',adapter),/duplicate stem/);
assert.equal(aligned,alignedBeforeFailure);
adapter.request.files=[];
await assert.rejects(vm.runInContext('runPhyloPrepBrowserWorkflow(request)',adapter),/Select at least one FASTA/);
assert.equal(aligned,alignedBeforeFailure);
console.log('Single input validation, pipeline-only product retention and standalone list/report outputs passed.');
