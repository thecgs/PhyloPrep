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
  function policyWarnings(copyMode, minimum, total, longestCopy) {
    const messages = [];
    if (copyMode === 'single-copy' && Number.isFinite(total) && total > 0 && minimum >= 2 && minimum < total) {
      messages.push(`This single-copy threshold permits missing taxa (${minimum} of ${total}). Taxa with multiple copies are omitted. To retain incomplete loci during downstream concatenation, use --missing-taxa pad-gaps; pad-gaps is the default. Explicit --missing-taxa skip-gene excludes loci missing a target taxon.`);
    }
    if (copyMode === 'all-copies' && longestCopy) {
      messages.push('Longest-copy selection chooses representatives by non-ambiguous length, not by inferred orthology. Different species may retain different paralogs, and long chimeric or fusion sequences may be selected. These are not established single-copy orthologs. Inspect copy_selection.tsv and independent orthology evidence before concatenation.');
    }
    return messages;
  }
  globalThis.PhyloPrepOrthogroupDefaults = { taxaFromHeader, selectedTaxonCount, policyWarnings };
})();
