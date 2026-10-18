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
  const files=new Map(),directories=new Set(['']);let time=0,writes=0;
  function directory(prefix){return {
    async getDirectoryHandle(name,{create=false}={}){
      const key=prefix+name+'/';if(!directories.has(key)&&!create)throw Error('not found');directories.add(key);return directory(key);
    },
    async getFileHandle(name,{create=false}={}){
      const key=prefix+name;if(!files.has(key)&&!create)throw Error('not found');
      return {async getFile(){const item=files.get(key);if(!item)throw Error('not found');return new File([item.blob],name,{lastModified:item.time});},
        async createWritable(){let blob;return {async write(value){blob=new Blob([value]);},async close(){files.set(key,{blob,time:++time});writes++;},async abort(){}};}};
    },
  };}
  return {root:directory(''),files,get writes(){return writes;}};
}
const disk=folder(),writer=save.createWriter(disk.root);
const engine={bundleSha256:'bundle-a',files:['MSAP.py','AA2Codon.py','trimAlnSeq.py','sequence_audit.py','msap_io.py'].map(name=>({name,sha256:digest(name)})),aligners:{mafft:{packageVersion:'test',maximumMemoryBytes:4294967296}}};
const input={name:'gene.fa',blob:new File(['>a\nATGGCT\n'],'gene.fa')};
async function fingerprint(file){const source=await save.readBlob(file);return {name:file.name,bytes:source.size,sha256:await hash.blob(source)};}
const options={seqtype:'codon',notrim:true};
const store=reuse.create(disk.root,engine,options,{fingerprint});
assert.equal(await store.restore(input),null);
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
