// Synchronous OPFS backing for the worker's ordinary Python file operations.
// Only directory/block metadata stays in JS; file contents never accumulate in MEMFS.
(() => {
  const BLOCK = 4096;
  function filesystem(FS, access, codes) {
    const memory = FS.filesystems.MEMFS, free = [], zero = new Uint8Array(BLOCK);
    let next = 0, files = 0, bytes = 0, peak = 0, reads = 0, writes = 0, failure;
    function io(callback) {
      try { return callback(); }
      catch (error) {
        failure = error.name === 'QuotaExceededError'
          ? 'Browser temporary disk is full. Free browser/site storage or use the native CLI.'
          : `Browser temporary disk I/O failed: ${error.message || error}`;
        throw new FS.ErrnoError(error.name === 'QuotaExceededError' ? codes.ENOSPC : codes.EIO);
      }
    }
    function writeAll(view, at) {
      let done = 0;
      while (done < view.byteLength) {
        const count = io(() => access.write(view.subarray(done), {at:at + done}));
        if (!count) { failure = 'Browser temporary disk could not complete a write.'; throw new FS.ErrnoError(codes.ENOSPC); }
        done += count;
      }
    }
    function resize(node, size) {
      if (size < 0) throw new FS.ErrnoError(codes.EINVAL);
      if (size < node.usedBytes) {
        for (const [index, block] of node.diskBlocks) {
          if (index * BLOCK >= size) { free.push(block); node.diskBlocks.delete(index); }
        }
        const remainder = size % BLOCK, block = node.diskBlocks.get(Math.floor(size / BLOCK));
        if (remainder && block !== undefined) writeAll(zero.subarray(remainder), block * BLOCK + remainder);
      }
      bytes += size - node.usedBytes; peak = Math.max(peak, bytes);
      node.usedBytes = size; node.timestamp = Date.now(); node.diskDirty = true;
    }
    function discard(node) {
      if (!node.diskBlocks || node.diskGone) return;
      node.diskUnlinked = true;
      if (node.diskOpen) return;
      for (const block of node.diskBlocks.values()) free.push(block);
      node.diskBlocks.clear(); bytes -= node.usedBytes; files--; node.diskGone = true;
    }
    function patch(node) {
      const original = node.node_ops;
      if (FS.isDir(node.mode)) {
        node.node_ops = {...original,
          mknod(parent, name, mode, dev) { return patch(memory.createNode(parent, name, mode, dev)); },
          unlink(parent, name) { const child=parent.contents[name]; original.unlink(parent,name); discard(child); },
          rename(source, parent, name) {
            const replaced = parent.contents[name];
            original.rename(source,parent,name);
            if (replaced && replaced !== source) discard(replaced);
          },
        };
      } else if (FS.isFile(node.mode)) {
        files++; node.diskBlocks = new Map(); node.diskOpen = 0;
        node.node_ops = {...original, setattr(target, attributes) {
          const {size, ...rest} = attributes;
          original.setattr(target,rest);
          if (size !== undefined) resize(target,size);
        }};
        const previous = node.stream_ops;
        node.stream_ops = {
          open(stream) { stream.node.diskOpen++; },
          close(stream) {
            const target=stream.node;
            try { if (target.diskDirty) { io(()=>access.flush()); target.diskDirty=false; } }
            finally { target.diskOpen--; if (target.diskUnlinked) discard(target); }
          },
          llseek: previous.llseek,
          allocate(stream, offset, length) { if (offset+length > stream.node.usedBytes) resize(stream.node,offset+length); },
          read(stream, buffer, offset, length, position) {
            const target=stream.node, count=Math.max(0,Math.min(length,target.usedBytes-position));
            let done=0;
            while (done<count) {
              const index=Math.floor((position+done)/BLOCK), within=(position+done)%BLOCK;
              const chunk=Math.min(count-done,BLOCK-within), block=target.diskBlocks.get(index);
              const view=buffer.subarray(offset+done,offset+done+chunk);
              if (block === undefined) view.fill(0);
              else {
                let got=0;
                while(got<chunk) {
                  const n=io(()=>access.read(view.subarray(got),{at:block*BLOCK+within+got}));
                  if (!n) { failure='Browser temporary disk returned incomplete data.'; throw new FS.ErrnoError(codes.EIO); }
                  got+=n;
                }
              }
              done+=chunk;
            }
            reads+=count; return count;
          },
          write(stream, buffer, offset, length, position) {
            const target=stream.node; let done=0;
            while(done<length) {
              const index=Math.floor((position+done)/BLOCK), within=(position+done)%BLOCK;
              const chunk=Math.min(length-done,BLOCK-within);
              let block=target.diskBlocks.get(index);
              if (block === undefined) {
                block=free.length ? free.pop() : next++;
                try { writeAll(zero,block*BLOCK); }
                catch(error) { free.push(block); throw error; }
                target.diskBlocks.set(index,block);
              }
              writeAll(buffer.subarray(offset+done,offset+done+chunk),block*BLOCK+within);
              done+=chunk;
            }
            if (position+length > target.usedBytes) resize(target,position+length);
            target.timestamp=Date.now(); target.diskDirty=true;
            writes+=length; return length;
          },
        };
      }
      return node;
    }
    return {
      mount(mount) { return patch(memory.mount(mount)); },
      stats() {return {kind:'opfs',cached_file_bytes:0,live_file_bytes:bytes,peak_file_bytes:peak,files,bytes_read:reads,bytes_written:writes};},
      get failure() {return failure;},
    };
  }
  async function open(runtime, path) {
    if (!globalThis.navigator?.storage?.getDirectory) return null;
    let directory, name, handle, file;
    const create = async () => {
      directory = await (await navigator.storage.getDirectory()).getDirectoryHandle('phyloprep-working-cache',{create:true});
      // Active workers hold an exclusive sync-access lock. Reclaim only unlocked
      // orphan files, e.g. after cancellation. Serialize creation against cleanup.
      if (navigator.locks) for await (const [oldName, old] of directory.entries()) {
        if (old.kind!=='file' || !oldName.startsWith('workspace-')) continue;
        try { const unused=await old.createSyncAccessHandle(); unused.close(); await directory.removeEntry(oldName); }
        catch {} // Another tab/worker can still own this cache.
      }
      name=`workspace-${crypto.randomUUID()}`;
      file=await directory.getFileHandle(name,{create:true});
      if (!file.createSyncAccessHandle) {await directory.removeEntry(name);return null;}
      handle=await file.createSyncAccessHandle();
      return handle;
    };
    try {
      const result = navigator.locks ? await navigator.locks.request('phyloprep-cache-maintenance',create) : await create();
      if (!result) return null;
    } catch(error) {
      if (file && !handle) {try {await directory.removeEntry(name);} catch {}}
      if (['NotAllowedError','SecurityError','NotSupportedError'].includes(error.name)) return null;
      throw new Error(`Could not create browser temporary disk: ${error.message}`);
    }
    const FS=runtime.FS, backend=filesystem(FS,handle,runtime.ERRNO_CODES);
    try {FS.mkdirTree(path);FS.mount(backend,{},path);}
    catch(error) {handle.close();await directory.removeEntry(name);throw error;}
    return {
      path, stats:backend.stats, get failure(){return backend.failure;},
      async close() {
        try {FS.unmount(path);FS.rmdir(path);}
        finally {try {handle.close();} finally {await directory.removeEntry(name);}}
      },
    };
  }
  globalThis.PhyloPrepDiskFS={open,filesystem};
})();
