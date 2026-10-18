import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import vm from 'node:vm';
import { PhyloPrepSave } from './test-output-store.mjs';

function diskFolder() {
  const files = new Map(); let clock = 0, failClose = false, reads = 0, truncateWrites = false;
  const folder = prefix => ({
    async getDirectoryHandle(name) { return folder(`${prefix}${name}/`); },
    async getFileHandle(name) {
      const key = prefix + name;
      return {
        async getFile() {
          reads++; const item = files.get(key);
          if (!item) throw new Error('File removed');
          return new File([item.text], name, { lastModified: item.time });
        },
        async createWritable() {
          let text;
          return {
            async write(value) { text = typeof value === 'string' ? value : await value.text(); },
            async close() { if (failClose) throw new Error('disk full'); files.set(key, {text:truncateWrites ? text.slice(0,-1) : text, time:++clock}); },
            async abort() {},
          };
        },
      };
    },
  });
  return {root:folder(''),files,get reads(){return reads;},fail(value){failClose=value;},truncate(value){truncateWrites=value;}};
}
const disk=diskFolder(), writer=PhyloPrepSave.createWriter(disk.root);
const file=PhyloPrepSave.record({name:'renamed/codon/raw/gene.aln.fasta',content:'>taxon\nATG\n'});
const alias=PhyloPrepSave.record({...file});
await writer.save(file);
assert.equal(file.content,undefined); assert.equal(alias.content,undefined);
assert.ok(writer.saved(alias));
assert.equal(await PhyloPrepSave.readText(alias),'>taxon\nATG\n');
assert.equal(file.content,undefined,'Read-back must not repopulate a permanent text cache');
assert.ok(!Object.hasOwn(PhyloPrepSave.wireFile(alias),'storage'));
assert.ok(PhyloPrepSave.wireFile(alias).disk.handle);
const reads=disk.reads;
await writer.save(alias);assert.equal(disk.reads,reads,'Confirmed saves must be deduplicated');
const other=diskFolder(), otherWriter=PhyloPrepSave.createWriter(other.root);
await otherWriter.save(alias);
assert.equal(other.files.get(file.name).text,'>taxon\nATG\n','Changing destination must read the saved source');
assert.equal(await PhyloPrepSave.readText(file),'>taxon\nATG\n');
const changed=PhyloPrepSave.record({...file,content:'new content'});
assert.notEqual(changed.storage,file.storage,'New content must detach from an old shared disk reference');
assert.equal(await PhyloPrepSave.readText(file),'>taxon\nATG\n');

const bad=diskFolder(), retry=PhyloPrepSave.createWriter(bad.root);
const pending=PhyloPrepSave.record({name:'matrix.fasta',content:'original'});
bad.fail(true);await retry.save(pending);
assert.equal(pending.content,'original');assert.equal(retry.saved(pending),false);
bad.fail(false);await retry.retry();await retry.save(pending);
assert.equal(pending.content,undefined);assert.equal(await PhyloPrepSave.readText(pending),'original');
bad.files.set('matrix.fasta',{text:'modified',time:100});
await assert.rejects(PhyloPrepSave.readText(pending),/changed on disk/);
const shortDisk=diskFolder(), shortWriter=PhyloPrepSave.createWriter(shortDisk.root);
const shortOutput=PhyloPrepSave.record({name:'renamed.fasta',content:'>a\nATG\n'});
shortDisk.truncate(true);await shortWriter.save(shortOutput);
assert.match(shortWriter.error.message,/Saved file size mismatch/);
assert.equal(shortWriter.saved(shortOutput),false);
assert.equal(shortOutput.content,'>a\nATG\n','A truncated save must not release the source');
shortDisk.truncate(false);await shortWriter.retry();await shortWriter.save(shortOutput);
assert.ok(shortWriter.saved(shortOutput));

const zipContext=vm.createContext({TextEncoder});
vm.runInContext(await readFile('site/browser-downloads.js','utf8'),zipContext);
let archive=[];
await zipContext.PhyloPrepDownloads.writeZip([file],{write:async bytes=>archive.push(bytes),close:async()=>{}},PhyloPrepSave.readBlob);
assert.ok(archive.length>0);assert.equal(file.content,undefined);
const whole=new Uint8Array(await new Blob(archive).arrayBuffer());
assert.ok(new TextDecoder().decode(whole).includes('>taxon\nATG\n'));

// Read-back staging verifies a saved alignment before it enters the analysis FS.
const worker=vm.createContext({TextEncoder,TextDecoder,Uint8Array,self:{addEventListener(){}},crypto:(await import('node:crypto')).webcrypto});
vm.runInContext(await readFile('site/browser-hash.js','utf8'),worker);
vm.runInContext(await readFile('site/pyodide-worker.js','utf8'),worker);
worker.file=PhyloPrepSave.wireFile(file);
assert.equal(await (await vm.runInContext('readDiskInput(file)',worker)).text(),'>taxon\nATG\n');
worker.file.sha256='0'.repeat(64);
await assert.rejects(vm.runInContext('readDiskInput(file)',worker),/checksum mismatch/);

// A reused renamed alignment is released only after confirmed disk saving.
let messages=[], saved=true, readsFromFS=0;
const backing=new Map(), released=[];
worker.runtime={FS:{
 open(path){if(!backing.has(path))throw Error('missing file');return {path};},
 close(){},
 read(stream,bytes,offset,length,position){readsFromFS++;const data=backing.get(stream.path);const count=Math.min(length,data.length-position);bytes.set(data.subarray(position,position+count),offset);return count;},
 analyzePath:path=>({exists:backing.has(path)}),
 unlink(path){assert.ok(backing.has(path));backing.delete(path);released.push(path);},
}};
Object.assign(worker,{setTimeout,clearTimeout});
worker.self.postMessage=(message,transfer=[])=>{
 messages.push(structuredClone(message,{transfer}));
 if(message.sequence!==undefined){
  worker.ack={sequence:message.sequence,saved};
  vm.runInContext('outputAcknowledgements.get(ack.sequence).resolve(ack); outputAcknowledgements.delete(ack.sequence)',worker);
 }
};
worker.files=[{name:'renamed/gene.aln.fasta',path:'/renamed/gene',size:10,sha256:(await import('node:crypto')).createHash('sha256').update('>a\nATGGCT\n').digest('hex')}];
backing.set('/renamed/gene',new TextEncoder().encode('>a\nATGGCT\n'));
await vm.runInContext('globalThis.transferred = new Map(); streamOutputFiles(runtime,1,files,"rename",true,transferred)',worker);
assert.equal(released.length,0,'The first transfer must keep the working copy for QC');
worker.files[0].release=true;const beforeReads=readsFromFS;
saved=false;
await vm.runInContext('streamOutputFiles(runtime,1,files,"after-QC",true,transferred)',worker);
assert.equal(released.length,0,'A failed save must preserve the working file');
saved=true;
await vm.runInContext('streamOutputFiles(runtime,1,files,"after-QC",true,transferred)',worker);
assert.deepEqual(released,['/renamed/gene']);
assert.equal(readsFromFS,beforeReads,'Reuse release must not retransmit the alignment');
await vm.runInContext('streamOutputFiles(runtime,1,files,undefined,true,transferred)',worker);
assert.equal(released.length,1,'Final archived outputs must not be unlinked twice');

// Same-size corruption must fail before a file-end/save acknowledgement.
messages=[];backing.set('/renamed/gene',new TextEncoder().encode('>a\nATGGCC\n'));
await assert.rejects(vm.runInContext('streamOutputFiles(runtime,2,files,"corrupt",true,new Map())',worker),/Output checksum mismatch.*renamed\/gene.aln.fasta/);
assert.ok(!messages.some(message=>message.type==='output-file-end'));
assert.ok(backing.has('/renamed/gene'),'A corrupt output must not be released');
console.log('Disk-backed release, shared aliases, lazy read-back, destination changes, failed-save retry, ZIP and changed-file checks passed.');
