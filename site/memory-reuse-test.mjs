import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash,randomBytes} from 'node:crypto';
import vm from 'node:vm';

const scope=vm.createContext({Blob,File,TextEncoder,TextDecoder,Uint8Array});
for(const name of ['browser-hash.js','browser-save.js','browser-reuse.js'])vm.runInContext(await readFile('site/'+name,'utf8'),scope);
const {PhyloPrepHash:hash,PhyloPrepSave:save,PhyloPrepReuse:reuse}=scope;
const digest=bytes=>createHash('sha256').update(bytes).digest('hex');
for(const size of [0,1,55,56,63,64,65,1024,1048579]){
  const bytes=randomBytes(size), incremental=new hash.SHA256();
  for(let i=0;i<size;i+=17)incremental.update(bytes.subarray(i,i+17));
  assert.equal(incremental.hex(),digest(bytes));
}
const large=new Blob(['ATG'.repeat(800000)]);let maxSlice=0;
assert.equal(await hash.blob({size:large.size,slice(start,end){maxSlice=Math.max(maxSlice,end-start);return large.slice(start,end);}}),digest(await large.bytes()));
assert.ok(maxSlice<=1048576);

function folder(){
  const files=new Map(),directories=new Set(['']);let time=0,writes=0,reads=0;
  function directory(prefix){return {
    async removeEntry(name){if(!files.delete(prefix+name))throw Object.assign(Error('not found'),{name:'NotFoundError'});},
    async getDirectoryHandle(name,{create=false}={}){
      const key=prefix+name+'/';if(!directories.has(key)&&!create)throw Object.assign(Error('not found'),{name:'NotFoundError'});directories.add(key);return directory(key);
    },
    async getFileHandle(name,{create=false}={}){
      const key=prefix+name;if(!files.has(key)&&!create)throw Error('not found');
      return {async getFile(){reads++;const item=files.get(key);if(!item)throw Error('not found');return new File([item.blob],name,{lastModified:item.time});},
        async createWritable(){let blob;return {async write(value){blob=new Blob([value]);},async close(){files.set(key,{blob,time:++time});writes++;},async abort(){}};}};
    },
  };}
  return {root:directory(''),files,get writes(){return writes;},get reads(){return reads;}};
}
const disk=folder(),writer=save.createWriter(disk.root);
const engine={bundleSha256:'bundle-a',files:['MSAP.py','AA2Codon.py','trimAlnSeq.py','sequence_audit.py','msap_io.py','rename_taxa.py'].map(name=>({name,sha256:digest(name)})),aligners:{mafft:{packageVersion:'test',maximumMemoryBytes:4294967296}}};
const input={name:'gene.fa',blob:new File(['>a\nATGGCT\n'],'gene.fa')};
async function fingerprint(file){const source=await save.readBlob(file);return {name:file.name,bytes:source.size,sha256:await hash.blob(source)};}
const options={seqtype:'codon',notrim:true};
const store=reuse.create(disk.root,engine,options,{fingerprint});
assert.equal(await store.restore(input),null);
assert.equal(await store.restoreBatch([input]),null);
const products=[];
for(const name of ['gene.mafft.codon.aln.fasta','gene.mafft.prot.aln.fasta','gene.sequence_changes.tsv']){
  const file=save.record({name,content:'>a\nATGGCT\n'});
  Object.assign(file,{size:10,sha256:digest(file.content)});
  const output=save.record({...file,name:'alignments/'+name});writer.save(output);products.push(file);
}
await writer.settle();
const checkpoint=await store.checkpoint(input,{engineMetadata:engine},products);
assert.ok(checkpoint);writer.save(checkpoint);await writer.settle();
const restored=await store.restore(input);assert.equal(restored.files.length,3);
assert.ok(restored.files.every(file=>file.content===undefined&&file.storage.disk));
const nextWriter=save.createWriter(disk.root),before=disk.writes;
for(const file of restored.files){file.name='alignments/'+file.name;nextWriter.adopt(file);nextWriter.save(file);}
await nextWriter.settle();assert.equal(disk.writes,before,'Verified reuse must not rewrite alignment files');
assert.equal(nextWriter.stats().saved,3);
assert.ok(await reuse.create(disk.root,engine,{...options,min_nucl_pis:50,missing_taxa:'pad-gaps',fourfold_missing_policy:'require-all'},{fingerprint}).restore(input),'QC changes do not invalidate alignment');
assert.ok(await reuse.create(disk.root,{...engine,bundleSha256:'new-qc-only-build'},options,{fingerprint}).restore(input),'Unchanged alignment sources survive downstream-only build changes');
for(const changed of [{genetic_code:2},{gap_ratio:0.5},{notrim:false},{seqtype:'prot'},{align_software:'muscle'}]){
  assert.equal(await reuse.create(disk.root,engine,{...options,...changed},{fingerprint}).restore(input),null);
}
assert.equal(await store.restore({...input,blob:new Blob(['>a\nATGGCC\n'])}),null);
const altered=structuredClone(engine);altered.files[0].sha256='changed';
assert.equal(await reuse.create(disk.root,altered,options,{fingerprint}).restore(input),null);
const path='alignments/gene.mafft.codon.aln.fasta',original=disk.files.get(path);
disk.files.set(path,{...original,blob:new Blob(['>a\nATGGCC\n'])});
assert.equal(await store.restore(input),null,'Same-size corrupt alignments must fail SHA validation');
disk.files.set(path,original);
// A complete stage hit reads only its manifest, not locus caches or saved products.
const batchCheckpoint=await store.checkpointBatch([input],products);
assert.ok(batchCheckpoint);writer.save(batchCheckpoint);await writer.settle();
const reads=disk.reads;
const fast=await store.restoreBatch([input]);
assert.equal(disk.reads-reads,1);assert.equal(fast.files.length,3);
assert.ok(fast.files.every(file=>file.storage.disk.root));
assert.equal(await save.readText(fast.files[0]),'>a\nATGGCT\n');
assert.equal(await store.restoreBatch([{...input,blob:new File(['>a\nATGGCT\n'],'gene.fa',{lastModified:input.blob.lastModified+1})}]),null);
assert.equal(await store.restoreBatch([{...input,blob:new Blob(['>a\nATGGCT\n'])}]),null);
assert.equal(await reuse.create(disk.root,engine,{...options,genetic_code:2},{fingerprint}).restoreBatch([input]),null);
disk.files.delete(path);
assert.ok(await store.restoreBatch([input]),'Stage hits defer output reads until actual use');
await assert.rejects(()=>save.readText(fast.files.find(file=>file.name.includes('.codon.'))),/not found/);
disk.files.set(path,original);
await store.clearBatch();assert.equal(await store.restoreBatch([input]),null);

// Renaming binds mapping content, source alignment hashes and renamer sources.
const aligned=products.filter(file=>file.name.endsWith('.aln.fasta'));
const mapping={name:'mapping.tsv',bytes:10,sha256:digest('mapping')};
const rename=reuse.createRename(disk.root,engine,mapping)('raw');
assert.deepEqual(Array.from((await rename.restore(aligned)).files),[]);
const renamed=aligned.map(file=>save.record({name:`renamed/${file.name.includes('.codon.')?'codon':'prot'}/raw/${file.name}`,content:'>taxon\nAA\n',size:10,sha256:digest('>taxon\nAA\n')}));
for(const file of renamed)writer.save(file);await writer.settle();
const renameCheckpoint=rename.checkpoint(aligned,renamed);assert.ok(renameCheckpoint);writer.save(renameCheckpoint);await writer.settle();
const renameReads=disk.reads, cachedRename=await rename.restore(aligned);
assert.equal(cachedRename.complete,true);assert.equal(disk.reads-renameReads,1);
assert.equal(cachedRename.files.length,2);
const changedAlignment=[{...aligned[0],sha256:digest('new alignment')},aligned[1]];
const partial=await rename.restore(changedAlignment);
assert.equal(partial.complete,false);assert.equal(partial.files.length,1);assert.ok(partial.files[0].name.includes('/prot/'));
assert.equal((await reuse.createRename(disk.root,engine,{...mapping,sha256:digest('new map')})('raw').restore(aligned)).files.length,0);
const newRenamer=structuredClone(engine);newRenamer.files.find(file=>file.name==='rename_taxa.py').sha256=digest('new implementation');
assert.equal((await reuse.createRename(disk.root,newRenamer,mapping)('raw').restore(aligned)).files.length,0);
await rename.clear();assert.equal((await rename.restore(aligned)).complete,false);
writer.save({name:'phyloprep-run.json',content:JSON.stringify({engine_metadata:engine,inputs:[{...mapping,input_key:'mapping'}],outputs:[
  ...aligned.map(file=>({name:'alignments/'+file.name,bytes:file.size,sha256:file.sha256})),
  ...renamed.map(file=>({name:file.name,bytes:file.size,sha256:file.sha256}))]})});await writer.settle();
const oldRunRename=await reuse.createRename(disk.root,engine,mapping)('raw').restore(aligned);
assert.equal(oldRunRename.complete,false,'Legacy manifests must use per-file checks, not the complete-stage fast path');
assert.equal(oldRunRename.files.length,2);


const cachePath=checkpoint.name;disk.files.delete(cachePath);
const previous={engine_metadata:engine,parameters:options,inputs:[await fingerprint(input)],outputs:products.map(file=>({name:'alignments/'+file.name,bytes:file.size,sha256:file.sha256}))};
writer.save({name:'phyloprep-run.json',content:JSON.stringify(previous)});await writer.settle();
assert.ok(await reuse.create(disk.root,engine,options,{fingerprint}).restore(input),'Prior run manifests can bootstrap reuse');
disk.files.delete(path);assert.equal(await reuse.create(disk.root,engine,options,{fingerprint}).restore(input),null);

const blobOutput=save.record({name:'large.fasta',blob:large});
assert.equal((await save.readBlob(blobOutput)).size,large.size);assert.equal(save.wireFile(blobOutput).blob.size,large.size);
writer.save(blobOutput);await writer.settle();assert.equal(blobOutput.storage.blob,undefined);assert.ok(blobOutput.storage.disk);
console.log('Bounded SHA-256, Blob release, verified reuse, no-rewrite adoption, cache invalidation and prior-manifest reuse passed.');

const app=await readFile('site/app.js','utf8');
let released=0;
Object.assign(scope,{activeAlignmentRun:undefined,activeWorkerRun:undefined,
  alignmentWorker:{terminate(){released++;}},pyodideWorker:{terminate(){released++;}}});
vm.runInContext(app.slice(app.indexOf('function releaseIdleWorkers('),app.indexOf('async function runMsapBrowserWorkflow(')),scope);
vm.runInContext('releaseIdleWorkers({pyodide:true})',scope);assert.equal(released,2);
scope.activeWorkerRun={};scope.pyodideWorker={terminate(){released++;}};
vm.runInContext('releaseIdleWorkers({pyodide:true})',scope);assert.equal(released,2,'Active workers must survive idle cleanup');
