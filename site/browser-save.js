// Shared output records: confirmed disk writes can release their in-memory text.
(() => {
  function record(file) {
    const descriptor = Object.getOwnPropertyDescriptor(file, 'content');
    if (!file.storage || (descriptor && 'value' in descriptor && descriptor.value !== file.storage.content)) {
      file.storage = { content: file.content, ...(file.blob ? {blob:file.blob} : {}) };
    }
    if (file.storage.blob) delete file.blob;
    if (Object.prototype.hasOwnProperty.call(file, 'content')) delete file.content;
    Object.defineProperty(file, 'content', {
      configurable: true, enumerable: false,
      get() { return this.storage.content; },
      set(value) { this.storage = { content: value }; },
    });
    return file;
  }
  async function readBlob(file) {
    const storage = file.storage;
    if (typeof (storage ? storage.content : file.content) === 'string') {
      return new Blob([storage ? storage.content : file.content]);
    }
    if (storage?.blob || file.blob) return storage?.blob || file.blob;
    const disk = storage?.disk || file.disk;
    if (!disk) throw new Error(`Output content is unavailable: ${file.name}`);
    const source = await disk.handle.getFile();
    if (source.size !== disk.size || source.lastModified !== disk.lastModified) {
      throw new Error(`Saved file changed on disk: ${file.name}. Restore the original file before retrying.`);
    }
    return source;
  }
  async function readText(file) { return (await readBlob(file)).text(); }
  function wireFile(file) {
    const { storage, ...metadata } = file;
    const content = storage ? storage.content : file.content;
    return typeof content === 'string' ? { ...metadata, content } : storage?.blob || file.blob
      ? {...metadata, blob:storage?.blob || file.blob} : { ...metadata, disk: storage?.disk || file.disk };
  }
  function createWriter(root, { concurrency = 4, onProgress = () => {} } = {}) {
    const directories = new Map([['', Promise.resolve(root)]]);
    const latest = new Map(), pending = new Set(), queue = [];
    let active = 0, failure = null;
    const counts = { saved: 0, failed: 0, pending: 0, writing: 0 };
    const bucket = status => ['saved', 'failed', 'writing'].includes(status) ? status : 'pending';
    function status(item, value) {
      if (latest.get(item.name) === item) {
        counts[bucket(item.status)]--; counts[bucket(value)]++;
      }
      item.status = value;
    }
    const stats = () => ({ total: latest.size, ...counts });
    const notify = (settled = false) => onProgress({ ...stats(), settled });
    function directory(path) {
      if (directories.has(path)) return directories.get(path);
      const split = path.lastIndexOf('/');
      const parent = split < 0 ? '' : path.slice(0, split);
      const name = split < 0 ? path : path.slice(split + 1);
      const promise = directory(parent).then(handle => handle.getDirectoryHandle(name, { create: true }));
      directories.set(path, promise);
      promise.catch(() => { if (directories.get(path) === promise) directories.delete(path); });
      return promise;
    }
    async function write(item) {
      const split = item.name.lastIndexOf('/');
      const handle = await directory(split < 0 ? '' : item.name.slice(0, split));
      const file = await handle.getFileHandle(item.name.slice(split + 1), { create: true });
      let stream;
      try {
        stream = await file.createWritable();
        const source = await readBlob(item.file);
        await stream.write(source);
        await stream.close();
        // Do not discard text unless close AND a readable disk snapshot succeed.
        if (file.getFile) {
          const snapshot = await file.getFile();
          if (snapshot.size !== source.size) throw new Error(`Saved file size mismatch: expected ${source.size} bytes, found ${snapshot.size}`);
          item.storage.disk = { handle: file, size: snapshot.size, lastModified: snapshot.lastModified };
          delete item.storage.content;
          delete item.storage.blob;
        }
      } catch (error) {
        if (stream) { try { await stream.abort(); } catch {} }
        throw error;
      }
    }
    function pump() {
      if (failure) {
        for (const item of queue.splice(0)) { status(item, 'pending'); item.finish(item); }
        notify();
        return;
      }
      while (active < concurrency && queue.length) {
        const item = queue.shift();
        active++;
        status(item, 'writing');
        void write(item).then(() => { status(item, 'saved'); }, error => {
          status(item, 'failed'); item.error = error;
          failure ||= new Error(`Could not save ${item.name}: ${error.message || error}`, { cause: error });
        }).finally(() => {
          active--; item.finish(item); notify(); pump();
        });
      }
      notify();
    }
    function matches(item, file) {
      return Boolean(item && (item.storage === file.storage ||
        (typeof file.content === 'string' && item.storage.content === file.content) ||
        (file.sha256 && file.sha256 === item.file.sha256 && file.size === item.file.size)));
    }
    function save(file) {
      record(file);
      const name = file.name.replaceAll('\\', '/');
      const previous = latest.get(name);
      if (matches(previous, file) && ['saved', 'queued', 'writing', 'waiting'].includes(previous.status)) {
        file.storage = previous.storage;
        return previous.promise;
      }
      const item = { name, file, storage: file.storage, status: 'waiting' };
      item.promise = new Promise(resolve => { item.finish = resolve; });
      if (previous) counts[bucket(previous.status)]--;
      latest.set(name, item); counts.pending++;
      pending.add(item.promise);
      void item.promise.then(() => pending.delete(item.promise));
      // Serialize replacements of the same path, even across a retry/version change.
      const enqueue = () => {
        if (failure) { status(item, 'pending'); item.finish(item); notify(); return; }
        status(item, 'queued'); queue.push(item); pump();
      };
      if (previous && ['waiting', 'queued', 'writing'].includes(previous.status)) void previous.promise.then(enqueue);
      else enqueue();
      return item.promise;
    }
    function adopt(file) {
      record(file);
      const name=file.name.replaceAll('\\','/');
      if(latest.has(name) || !file.storage.disk)return;
      const item={name,file,storage:file.storage,status:'saved'};
      item.promise=Promise.resolve(item);latest.set(name,item);counts.saved++;notify();
    }
    async function settle() { while (pending.size) await Promise.all([...pending]); notify(true); return stats(); }
    async function ready(limit = 8) { while (!failure && pending.size > limit) await Promise.race([...pending]); }
    function state(file) {
      const item = latest.get(file.name.replaceAll('\\', '/'));
      return matches(item, file) ? item.status : 'pending';
    }
    function saved(file) { return state(file) === 'saved'; }
    async function retry() { await settle(); failure = null; }
    return { root, save, adopt, settle, ready, saved, state, retry, stats, get error() { return failure; } };
  }
  globalThis.PhyloPrepSave = { createWriter, record, readBlob, readText, wireFile };
})();
