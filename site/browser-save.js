// A per-run local-directory writer. Files remain available to the analysis engine.
(() => {
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
        await stream.write(item.content);
        await stream.close();
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
    function save(file) {
      const name = file.name.replaceAll('\\', '/');
      const previous = latest.get(name);
      if (previous?.content === file.content && ['saved', 'queued', 'writing', 'waiting'].includes(previous.status)) {
        previous.content = file.content;
        return previous.promise;
      }
      const item = { name, content: file.content, status: 'waiting' };
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
    async function settle() { while (pending.size) await Promise.all([...pending]); notify(true); return stats(); }
    async function ready(limit = 8) { while (!failure && pending.size > limit) await Promise.race([...pending]); }
    function state(file) { const item = latest.get(file.name.replaceAll('\\', '/')); return item?.content === file.content ? item.status : 'pending'; }
    function saved(file) { return state(file) === 'saved'; }
    async function retry() { await settle(); failure = null; }
    return { root, save, settle, ready, saved, state, retry, stats, get error() { return failure; } };
  }
  globalThis.PhyloPrepSave = { createWriter };
})();
