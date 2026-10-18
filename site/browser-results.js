(() => {
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
  const category = name => name === 'phyloprep-run.json' ? 'run' : name === 'phyloprep-report.html' || name.startsWith('reports/') ? 'reports' : name.startsWith('matrices/') ? 'matrices' : name.startsWith('qc/') ? 'qc' : 'other';
  const rank = { run: 0, reports: 1, matrices: 2, qc: 3, other: 4 };
  function createEntries(files) {
    return files.map((file, index) => ({ file, index, category: category(file.name) }))
      .sort((a, b) => rank[a.category] - rank[b.category] || collator.compare(a.file.name, b.file.name) || a.index - b.index);
  }
  function pageEntries(entries, { folder = 'all', query = '', page = 0, pageSize = 200 } = {}) {
    const text = query.trim().toLowerCase();
    const filtered = entries.filter(entry => (folder === 'all' || entry.category === folder) && (!text || entry.file.name.toLowerCase().includes(text)));
    const pages = Math.max(1, Math.ceil(filtered.length / pageSize));
    page = Math.max(0, Math.min(page, pages - 1));
    const start = page * pageSize;
    return { entries: filtered.slice(start, start + pageSize), total: filtered.length, pages, page, start };
  }
  globalThis.PhyloPrepResults = { createEntries, pageEntries };
})();
