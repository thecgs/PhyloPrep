import { readFile } from 'node:fs/promises';
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import vm from 'node:vm';
const app = await readFile('site/app.js', 'utf8');
const context = vm.createContext({ document: { querySelector: () => null } });
vm.runInContext(app.slice(0, app.indexOf('let activeModuleIndex')), context);
const modules = vm.runInContext('modules', context);
const captured = spawnSync('python3', ['tests/inspect_cli_options.py', ...modules.map(module => module.script)], { encoding: 'utf8' });
assert.equal(captured.status, 0, captured.stderr);
const cli = JSON.parse(captured.stdout);
for (const module of modules) {
  const actions = cli[module.script];
  assert.ok(Array.isArray(actions), `${module.script}: CLI parser must be inspectable`);
  for (const field of module.inputFields) {
    for (const flag of field.flag.match(/--?[A-Za-z][A-Za-z0-9_-]*/g) || []) {
      assert.ok(actions.some(action => action.flags.includes(flag)), `${module.script}: unknown input flag ${flag}`);
    }
  }
  for (const parameter of module.parameters) {
    const flag = parameter.flag.split(/[ ,]/)[0];
    const action = actions.find(action => action.flags.includes(flag));
    assert.ok(action, `${module.script}: unknown browser flag ${flag}`);
    if (parameter.token === 'thread' && module.browserWorkflow) {
      assert.equal(parameter.value, '1');
      assert.match(parameter.description, /fixed|compatibility/i);
    } else if (parameter.type === 'select' && action.default != null) {
      assert.equal(parameter.value, String(action.default), `${module.script} ${flag}: mismatched dropdown default`);
    }
    if (parameter.type === 'select' && action.choices) {
      for (const [value] of parameter.options) assert.ok(action.choices.includes(value), `${module.script} ${flag}: unsupported ${value}`);
    }
    if (action.default != null && ['number', 'text'].includes(parameter.type)) {
      const expected = String(action.default);
      assert.equal(Number.isFinite(Number(expected)) ? Number(parameter.value) : parameter.value,
                   Number.isFinite(Number(expected)) ? Number(expected) : expected,
                   `${module.script} ${flag}: mismatched fixed default`);
      assert.notEqual(parameter.value, '', `${module.script} ${flag}: blank fixed default`);
    }
  }
}
const extractor = modules.find(module => module.script === 'extract_orthofinder_orthogroups.py');
assert.deepEqual(Array.from(extractor.parameters.find(p => p.token === 'copy_mode').options, pair => Array.from(pair)), [['single-copy', 'single-copy'], ['all-copies', 'all-copies']]);
assert.ok(!modules.find(m => m.script === 'rename_taxa.py').parameters.some(p => p.flag === '--keep-special'));
const workerContext = vm.createContext({ self: { addEventListener() {} } });
vm.runInContext(await readFile('site/pyodide-worker.js', 'utf8'), workerContext);
workerContext.stage = { fastaFiles: [{path:'/cds.fa'}, {path:'/prot.fa'}], outputRoot:'/out' };
workerContext.request = { script:'AA2Codon.py', options:{ genetic_code:'2' }, outputName:'codon.fa' };
assert.ok(!vm.runInContext('originalArguments(request, stage)', workerContext).includes('-g'));
workerContext.request.options.validate_translation = true;
const args = vm.runInContext('originalArguments(request, stage)', workerContext);
assert.equal(args[args.indexOf('-g') + 1], '2');
const defaultsContext = vm.createContext({});
vm.runInContext(await readFile('site/orthogroup-defaults.js', 'utf8'), defaultsContext);
const { taxaFromHeader, selectedTaxonCount } = defaultsContext.PhyloPrepOrthogroupDefaults;
const taxa = taxaFromHeader('Orthogroup\tHuman\tMouse\tFish\nOG1\tg1\tg2\tg3');
assert.equal(selectedTaxonCount(taxa, '', ''), 3);
assert.equal(selectedTaxonCount(taxa, 'Human Mouse', ''), 2);
assert.equal(selectedTaxonCount(taxa, '', 'Fish'), 2);
assert.equal(taxaFromHeader('"Orthogroup"\t"Human"\t"Mouse"\n').length, 2);
console.log('19-tool CLI option/default audit and dynamic-default checks passed.');
const rawSource = app.slice(app.indexOf('function rawParameterValues()'), app.indexOf('function valuesForCommand()'));
vm.runInContext(rawSource, context);
const controls = new Map(extractor.parameters.map(parameter => [parameter.token || parameter.flag.replace(/^--?/, '').replaceAll('-', '_'), {
  type: parameter.type, value: parameter.value, checked: parameter.value === true,
}]));
controls.get('min_taxa').value = '3'; controls.get('min_taxa').dataset = { autoDefault:'true' };
context.parameterToken = parameter => parameter.token || parameter.flag.replace(/^--?/, '').replaceAll('-', '_');
context.activeModule = () => extractor;
context.testPanel = { querySelector(selector) { return controls.get(selector.match(/data-parameter="([^"]+)"/)?.[1]); } };
vm.runInContext('elements.parameterList = testPanel', context);
assert.equal(vm.runInContext('rawParameterValues().min_taxa', context), undefined);
controls.get('min_taxa').dataset.autoDefault = 'false'; controls.get('min_taxa').value = '2';
assert.equal(vm.runInContext('rawParameterValues().min_taxa', context), '2');
controls.get('max_copies_per_taxon').value = '2';
controls.get('longest_copy').checked = true;
controls.get('copy_mode').value = 'single-copy';
assert.equal(vm.runInContext('rawParameterValues().longest_copy', context), undefined);
assert.equal(vm.runInContext('rawParameterValues().max_copies_per_taxon', context), undefined);
controls.get('copy_mode').value = 'all-copies';
assert.equal(vm.runInContext('rawParameterValues().longest_copy', context), true);
assert.equal(vm.runInContext('rawParameterValues().max_copies_per_taxon', context), '2');
workerContext.runtime = { FS: { mkdirTree() {}, writeFile() {} } };
workerContext.stage = { optionsRoot:'/opts', outputRoot:'/out', files:[
  {inputKey:'orthogroups',path:'/og.tsv'},
  {inputKey:'sequences',name:'Human.fa',fasta:true,content:'>g\nAAA'},
] };
workerContext.request = { options:{copy_mode:'all-copies',longest_copy:true}, outputName:'results' };
assert.ok(vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('--longest-copy'));
assert.ok(vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('codon'));
workerContext.request.options.max_copies_per_taxon = '2';
const copyArgs = vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext);
assert.equal(copyArgs[copyArgs.indexOf('--max-copies-per-taxon') + 1], '2');
workerContext.request.options.seqtype = 'prot';
assert.ok(vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('prot'));
workerContext.request.options.longest_copy = false;
assert.ok(vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('--max-copies-per-taxon'));
workerContext.request.options.copy_mode = 'single-copy';
assert.ok(!vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('--longest-copy'));
assert.ok(!vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('-st'));
assert.ok(!vm.runInContext('orthofinderArguments(runtime, request, stage)', workerContext).includes('--max-copies-per-taxon'));
const policyWarnings = defaultsContext.PhyloPrepOrthogroupDefaults.policyWarnings;
assert.match(policyWarnings('single-copy', 80, 100, false).join(' '), /--missing-taxa pad-gaps/);
assert.equal(policyWarnings('single-copy', 80, 80, false).length, 0);
assert.equal(policyWarnings('single-copy', 80, null, false).length, 0);
assert.match(policyWarnings('all-copies', 100, 100, true).join(' '), /different paralogs/);
assert.equal(policyWarnings('all-copies', 100, 100, false).length, 0);
assert.ok(app.includes('...(result.warnings || [])'), 'Successful extraction must display policy warnings');
