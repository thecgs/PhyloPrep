import { PhyloPrepSave } from './test-output-store.mjs';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash, webcrypto} from 'node:crypto';
import vm from 'node:vm';
const app=await readFile('site/app.js','utf8');
const worker=await readFile('site/pyodide-worker.js','utf8');
const sha=value=>createHash('sha256').update(value).digest('hex');
const manifest=JSON.parse(await readFile('site/engine/manifest.json','utf8'));
const records=[...manifest.files,manifest.runner];
assert.equal(manifest.bundleSha256,sha(records.map(({name,sha256})=>`${name}:${sha256}`).join('\n')));
for(const file of records){
  assert.equal(sha(await readFile(`site/engine/${file.url}`)),file.sha256);
  assert.ok(file.url.includes(`/`+file.sha256+`/`));
}
const installed=new Map();let executed=0, corrupt, repairOnRefresh=false;
const loader=vm.createContext({ PhyloPrepSave, Blob, crypto:webcrypto,TextEncoder,TextDecoder,Uint8Array,
  setTimeout,clearTimeout,self:{addEventListener:()=>{}},
  fetch:async (url,options)=>url==='engine/manifest.json'?new Response(JSON.stringify(manifest)):
    new Response(corrupt && url.endsWith(`/${corrupt}`) && (!repairOnRefresh || options?.cache!=='no-cache')?'mismatched cached source':await readFile(`site/${url}`)),
  runtime:{FS:{mkdirTree:()=>{},writeFile:(path,value)=>installed.set(path,value)},runPythonAsync:async()=>{executed++;}},
});
vm.runInContext(worker,loader);
await vm.runInContext('installOriginalEngine(runtime)',loader);
assert.equal(installed.size,records.length);assert.equal(executed,1);
for(const name of [manifest.files[0].name,'browser_runner.py']){
  installed.clear();executed=0;corrupt=name;
  await assert.rejects(vm.runInContext('installOriginalEngine(runtime)',loader),/Source checksum mismatch/);
  assert.equal(installed.size,0,'No source may be installed before all checksums pass');
  assert.equal(executed,0);
}
installed.clear();executed=0;repairOnRefresh=true;
await vm.runInContext('installOriginalEngine(runtime)',loader);
assert.equal(executed,1,'A corrupt cached source must recover after a successful network revalidation');
assert.equal(installed.size,records.length);
repairOnRefresh=false;corrupt=undefined;
const validBundle=manifest.bundleSha256;manifest.bundleSha256='0'.repeat(64);
await assert.rejects(vm.runInContext('installOriginalEngine(runtime)',loader),/bundle checksum is invalid/);
manifest.bundleSha256=validBundle;

// Manifest refresh must not return a stale cache entry when the network has a newer build.
let fetchHandler, responsePromise, fetches=0, offline=false;
const old=new Response('old manifest');
const sw=vm.createContext({ PhyloPrepSave, Blob, URL,self:{location:{origin:'https://test.example'},addEventListener:(type,callback)=>{if(type==='fetch')fetchHandler=callback;}},
  caches:{match:async()=>old.clone(),open:async()=>({put:async()=>{}})},
  fetch:async()=>{fetches++;if(offline)throw new Error('offline');return new Response('new manifest');},
});
vm.runInContext(await readFile('site/sw.js','utf8'),sw);
fetchHandler({request:new Request('https://test.example/engine/manifest.json',{cache:'no-cache'}),respondWith:p=>{responsePromise=p;}});
assert.equal(await (await responsePromise).text(),'new manifest');assert.equal(fetches,1);
offline=true;
fetchHandler({request:new Request('https://test.example/engine/manifest.json',{cache:'no-cache'}),respondWith:p=>{responsePromise=p;}});
assert.equal(await (await responsePromise).text(),'old manifest','Offline fallback must remain available for verified cached sources');

let encodes=0,digests=0;
const fingerprints=vm.createContext({ PhyloPrepSave, Blob, TextEncoder:class extends TextEncoder{encode(value){encodes++;return super.encode(value);}},
  crypto:{subtle:{digest:async(...args)=>{digests++;return webcrypto.subtle.digest(...args);}}},Uint8Array,
  elements:{command:{textContent:'changed command -g 2'}},versionManifest:undefined,
  navigator:{userAgent:'test'},implementationProfile:()=>({}),externalAlignerLabel:()=>'',
});
vm.runInContext(await readFile('site/browser-hash.js','utf8'),fingerprints);
let blobHashes=0;const hashBlob=fingerprints.PhyloPrepHash.blob;fingerprints.PhyloPrepHash.blob=async file=>{blobHashes++;return hashBlob(file);};
vm.runInContext(app.slice(app.indexOf('async function sha256Text('),app.indexOf('function transferTarget(')),fingerprints);
fingerprints.module={script:'phyloprep.py'};
fingerprints.request={equivalentCli:'original command -g 1',options:{genetic_code:1},files:[]};
fingerprints.result={};fingerprints.outputs=[{name:'matrix.fasta',content:'A'.repeat(8*1024*1024),size:8*1024*1024,sha256:'a'.repeat(64)}];
let record=JSON.parse(await vm.runInContext('buildRunManifest(module,request,result,outputs)',fingerprints));
assert.equal(record.equivalent_cli,'original command -g 1');assert.equal(record.parameters.genetic_code,1);
assert.equal(encodes,0,'A previously fingerprinted matrix must not be encoded again');assert.equal(digests,0);
fingerprints.outputs=[{name:'unicode.fasta',content:'é😀\n'}];
record=JSON.parse(await vm.runInContext('buildRunManifest(module,request,result,outputs)',fingerprints));
assert.equal(encodes,0);assert.equal(digests,0);assert.equal(blobHashes,1);
assert.equal(record.outputs[0].bytes,Buffer.byteLength('é😀\n'));
assert.equal(record.outputs[0].sha256,sha('é😀\n'));
await vm.runInContext('buildRunManifest(module,request,result,outputs)',fingerprints);
assert.equal(blobHashes,1,'Fallback fingerprints must be reused too');

const roots=[{},{},{}];
const config=vm.createContext({ PhyloPrepSave, Blob, configurationDisabled:new Map(),elements:{parameterList:roots[0],inputFileFields:roots[1],scriptList:roots[2],outputName:{disabled:false},clearInputFiles:{disabled:true}}});
vm.runInContext(app.slice(app.indexOf('function setConfigurationLocked('),app.indexOf('function setRunControls(')),config);
vm.runInContext('setConfigurationLocked(true)',config);
assert.ok(roots.every(root=>root.inert));assert.equal(config.elements.outputName.disabled,true);
vm.runInContext('setConfigurationLocked(false)',config);
assert.ok(roots.every(root=>!root.inert));assert.equal(config.elements.outputName.disabled,false);
assert.equal(config.elements.clearInputFiles.disabled,true,'Previously disabled controls must stay disabled');
console.log('Source/adapter/bundle verification, manifest refresh/offline cache, frozen provenance and single-pass fingerprints passed.');
