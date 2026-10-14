(() => {
  function taxaFromHeader(text) {
    const line = text.replace(/^\uFEFF/, '').split(/\r?\n/).find(line => line.trim());
    if (!line) return [];
    const values = [], parts = [];
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const char = line[i];
      if (char === '"') {
        if (quoted && line[i + 1] === '"') { parts.push('"'); i++; }
        else quoted = !quoted;
      } else if (char === '\t' && !quoted) { values.push(parts.join('').trim()); parts.length = 0; }
      else parts.push(char);
    }
    values.push(parts.join('').trim());
    return values[0] === 'Orthogroup' ? values.slice(1) : [];
  }
  function selectedTaxonCount(taxa, keepText, dropText) {
    const names = text => new Set(String(text || '').split(/[\s,]+/).filter(Boolean));
    const keep = names(keepText), drop = names(dropText);
    return taxa.filter(taxon => (!keep.size || keep.has(taxon)) && !drop.has(taxon)).length;
  }
  globalThis.PhyloPrepOrthogroupDefaults = { taxaFromHeader, selectedTaxonCount };
})();
