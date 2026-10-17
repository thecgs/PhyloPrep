import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import vm from 'node:vm';
import { PhyloPrepSave } from './test-output-store.mjs';

function diskFolder() {
  const files = new Map(); let clock = 0, failClose = false, reads = 0;
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
            async close() { if (failClose) throw new Error('disk full'); files.set(key, {text, time:++clock}); },
            async abort() {},
          };
        },
      };
    },
  });
  return {root:folder(''),files,get reads(){return reads;},fail(value){failClose=value;}};
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
console.log('Disk-backed release, shared aliases, lazy read-back, destination changes, failed-save retry, ZIP and changed-file checks passed.');
