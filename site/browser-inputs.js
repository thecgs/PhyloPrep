// Shared, testable local-file loading. Browser directory access is user-authorized.
(() => {
  const normalize = path => String(path).replaceAll("\\", "/").replace(/^\.\//, "");
  const basename = path => normalize(path).split("/").at(-1);
  const isList = file => /\.(pathlist|list)$/i.test(file.name);
  async function readText(file) {
    try { return await file.text(); }
    catch (error) {
      throw new Error(`Cannot read ${file.name}: ${error.message}. Reselect this file or authorize its input folder and retry.`, { cause: error });
    }
  }
  async function readBounded(files, read = readText, concurrency = 4) {
    const results = new Array(files.length);
    let next = 0, stopped = false;
    await Promise.all(Array.from({ length: Math.min(concurrency, files.length) }, async () => {
      while (!stopped && next < files.length) {
        const index = next++;
        try { results[index] = await read(files[index]); }
        catch (error) { stopped = true; throw error; }
      }
    }));
    return results;
  }
  async function resolveFiles(files, { sibling, read = readText } = {}) {
    const lists = files.filter(isList);
    if (!lists.length) return files;
    const available = files.filter(file => !isList(file));
    const exact = new Map(), byBase = new Map();
    for (const file of available) {
      const keyPath = normalize(file.webkitRelativePath || file.name);
      if (exact.has(keyPath) && exact.get(keyPath) !== file) throw new Error(`Ambiguous FASTA path: ${keyPath}. Select the matching folder.`);
      exact.set(keyPath, file);
      const key = basename(file.name);
      if (!byBase.has(key)) byBase.set(key, []);
      byBase.get(key).push(file);
    }
    const result = [], seen = new Set();
    for (const list of lists) {
      const content = await read(list);
      for (const line of content.split(/\r?\n/)) {
        const entry = line.trim();
        if (!entry || /^[#;]/.test(entry)) continue;
        const listPath = normalize(list.webkitRelativePath || list.name);
        const roots = new Set(available.filter(file => file.webkitRelativePath).map(file => normalize(file.webkitRelativePath).split("/")[0]));
        const parent = listPath.includes("/") ? listPath.slice(0, listPath.lastIndexOf("/") + 1)
          : roots.size === 1 ? `${[...roots][0]}/` : "";
        let file = exact.get(parent + normalize(entry)) || exact.get(normalize(entry));
        // A folder authorization can resolve subpaths and the same-directory fallback.
        if (!file && sibling) file = await sibling(entry, list);
        if (!file) file = exact.get(parent + basename(entry));
        if (!file) {
          const matches = byBase.get(basename(entry)) || [];
          if (matches.length > 1) throw new Error(`Ambiguous FASTA basename in ${list.name}: ${entry}. Select the matching folder.`);
          file = matches[0];
        }
        if (!file) throw new Error(`Missing FASTA in ${list.name}: ${entry}. Select its input folder or the referenced FASTA files; the same-folder basename ${basename(entry)} was also checked.`);
        if (!seen.has(file)) { result.push(file); seen.add(file); }
      }
    }
    if (!result.length) throw new Error("The selected pathlist contains no FASTA paths.");
    return result;
  }
  globalThis.PhyloPrepInputs = { readText, readBounded, resolveFiles, normalize, basename };
})();
