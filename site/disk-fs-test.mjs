import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import vm from 'node:vm';
const scope=vm.createContext({Uint8Array,Map});
vm.runInContext(await readFile('site/browser-disk-fs.js','utf8'),scope);
const {filesystem,open}=scope.PhyloPrepDiskFS;
assert.equal(await open({},'/work'),null,'Unavailable OPFS must retain the existing fallback');
const codes={ENOSPC:51,EIO:29,EINVAL:28};
let data=new Uint8Array(0),fail=false,flushes=0,partial=Infinity;
const access={
 write(view,{at}) {
  if(fail) {const e=new Error('quota exceeded');e.name='QuotaExceededError';throw e;}
  if(data.length<at+view.length){const next=new Uint8Array(Math.max(at+view.length,data.length*2,4096));next.set(data);data=next;}
  const length=Math.min(view.length,partial);data.set(view.subarray(0,length),at);return length;
 },
 read(view,{at}) {const length=Math.min(view.length,partial);view.set(data.subarray(at,at+length));return length;},
 flush(){flushes++;},
};
const FS={
 isDir:mode=>(mode&0xf000)===0x4000,isFile:mode=>(mode&0xf000)===0x8000,
 ErrnoError:class extends Error{constructor(errno){super(String(errno));this.errno=errno;}},
 filesystems:{},
};
const memory={
 mount(){return this.createNode(null,'/',0x4000);},
 createNode(parent,name,mode){
  const node={name,parent,mode,usedBytes:0,contents:FS.isDir(mode)?{}:null};
  if(parent)parent.contents[name]=node;
  node.node_ops={
   setattr(node,attr){if(attr.size!==undefined)throw Error('MEMFS data allocation was reached');Object.assign(node,attr);},
   unlink(parent,name){delete parent.contents[name];},
   rename(source,parent,name){delete source.parent.contents[source.name];parent.contents[name]=source;source.parent=parent;source.name=name;},
  };
  node.stream_ops={llseek(stream,offset,whence){return offset+(whence===2?stream.node.usedBytes:whence===1?stream.position:0);}};
  return node;
 },
};
FS.filesystems.MEMFS=memory;
const disk=filesystem(FS,access,codes),root=disk.mount({});
const create=name=>root.node_ops.mknod(root,name,0x8000,0);
const write=(node,bytes,position=0)=>node.stream_ops.write({node},bytes,0,bytes.length,position);
const read=(node,size=node.usedBytes)=>{const bytes=new Uint8Array(size);node.stream_ops.read({node},bytes,0,size,0);return bytes;};
const loci=[];
for(let i=0;i<30;i++){
 const node=create(`gene${i}.fasta`),bytes=new Uint8Array(100000).fill(i+1),stream={node};
 node.stream_ops.open(stream);write(node,bytes);node.stream_ops.close(stream);loci.push(node);
 assert.equal(node.contents,null,'Completed loci must not retain file bytes in MEMFS');
 assert.equal(disk.stats().cached_file_bytes,0);
}
assert.ok(flushes>=30,'Each closed locus must be flushed');
for(let i=0;i<loci.length;i++)assert.deepEqual(read(loci[i]),new Uint8Array(100000).fill(i+1));
const hole=create('hole');write(hole,new Uint8Array([7,8]),8197);
assert.deepEqual(read(hole).slice(0,8197),new Uint8Array(8197));
hole.node_ops.setattr(hole,{size:8198});hole.node_ops.setattr(hole,{size:8200});
assert.deepEqual(read(hole).slice(8197),new Uint8Array([7,0,0]));
const before=disk.stats().files;
root.node_ops.rename(loci[0],root,'gene1.fasta');
assert.equal(disk.stats().files,before-1,'Atomic replacement must reclaim overwritten disk blocks');
assert.deepEqual(read(root.contents['gene1.fasta']),new Uint8Array(100000).fill(1));
const stream={node:loci[2]};loci[2].stream_ops.open(stream);
root.node_ops.unlink(root,'gene2.fasta');
assert.deepEqual(read(loci[2]),new Uint8Array(100000).fill(3),'Unlink must preserve an open reader');
loci[2].stream_ops.close(stream);
const short=create('reused');write(short,new Uint8Array([1]));short.node_ops.setattr(short,{size:4096});
assert.deepEqual(read(short).slice(1),new Uint8Array(4095),'Reused blocks must not expose old locus bytes');
partial=7;
const partialNode=create('partial');write(partialNode,new Uint8Array(31).fill(12));
assert.deepEqual(read(partialNode),new Uint8Array(31).fill(12),'Short disk reads/writes must complete without losing bytes');
partial=Infinity;
fail=true;
assert.throws(()=>write(create('quota'),new Uint8Array([1])),error=>error.errno===codes.ENOSPC);
assert.match(disk.failure,/temporary disk is full/);
console.log('Per-locus disk flush/read-back, no MEMFS file buffers, sparse/truncated data, replacement, open-unlink and quota checks passed.');
