// Real Pyodide + FileSystemFileHandle integration; no aligner or dialogs required.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { chromium } = require('playwright');
const root = path.resolve('site');
const server = http.createServer((req,res) => {
  const file = path.resolve(root, '.' + new URL(req.url,'http://localhost').pathname);
  if (!file.startsWith(root + path.sep) || !fs.existsSync(file) || !fs.statSync(file).isFile()) {res.writeHead(404).end();return;}
  res.setHeader('Content-Type',file.endsWith('.js') ? 'application/javascript' : file.endsWith('.json') ? 'application/json' : 'text/html');
  fs.createReadStream(file).pipe(res);
});
(async()=>{
 await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
 const browser=await chromium.launch({args:['--no-sandbox']});
 const timeout=setTimeout(()=>void browser.close(),180000);
 try{
  const context=await browser.newContext({serviceWorkers:"block"});
  if(process.env.PHYLOPREP_TEST_PYODIDE_DIR) await context.route('https://cdn.jsdelivr.net/pyodide/**',async route=>{
    const file=path.join(process.env.PHYLOPREP_TEST_PYODIDE_DIR,path.basename(new URL(route.request().url()).pathname));
    if(fs.existsSync(file))await route.fulfill({path:file,headers:{'access-control-allow-origin':'*'},contentType:file.endsWith('.js')?'application/javascript':file.endsWith('.wasm')?'application/wasm':undefined});
    else await route.continue();
  });
  const page=await context.newPage();
  page.on('console',message=>{if(message.type()==='error')console.error(message.text());});
  await page.goto(`http://127.0.0.1:${server.address().port}/index.html`);
  const result=await page.evaluate(async()=>{
    const folder=await (await navigator.storage.getDirectory()).getDirectoryHandle('outputs',{create:true});
    outputDirectoryHandle=folder; outputSaveSession=createOutputSaveSession(folder);producedFiles=new Map();
    const dna='>gene_a\nATGGCT\n>gene_b\nATGGCT\n', protein='>gene_a\nMA\n>gene_b\nMA\n';
    const batch=[];
    for(const variant of ['', '.trimmed']) for(const type of ['codon','prot']){
      const name=`gene.mafft.${type}${variant}.aln.fasta`;
      const file=PhyloPrepSave.record({name:`alignments/${name}`,content:type==='codon'?dna:protein});
      await fileFingerprint(file);await publishProducedFiles([file]);await outputSaveSession.settle();
      if(file.content!==undefined||!file.storage.disk)throw Error('Saved alignment text was retained');
      batch.push(PhyloPrepSave.record({...file,name,role:'batch-output'}));
    }
    const result=await runInWorker({type:'phyloprep-finalize',script:'phyloprep.py',files:[
      {name:'gene.fa',role:'source',content:dna},
      {name:'map.tsv',role:'mapping',content:'gene_a\tsp_a\ngene_b\tsp_b\n'},...batch],
      options:{seqtype:'codon',genetic_code:1},outputName:'',usesPathList:false});
    await outputSaveSession.settle();
    if(outputSaveSession.error)throw outputSaveSession.error;
    if(result.files.some(file=>file.content!==undefined||!file.storage.disk))throw Error('Final result text was retained');
    const renamed=result.files.filter(file=>file.name.startsWith('renamed/'));
    const matrix=result.files.find(file=>file.name==='matrices/codon/raw/supermatrix.fasta');
    const paired=result.files.find(file=>file.name==='matrices/paired/codon/raw/supermatrix.fasta');
    const expected='>sp_a\nATGGCT\n>sp_b\nATGGCT\n';
    if(await PhyloPrepSave.readText(matrix)!==expected||await PhyloPrepSave.readText(paired)!==expected)throw Error('Matrix changed');
    generatedFiles=result.files;renderGeneratedFiles();selectedGeneratedFileIndex=result.files.indexOf(matrix);
    await previewGeneratedFile(matrix);
    if(elements.preview.textContent!==expected||matrix.content!==undefined)throw Error('Lazy preview failed');
    const destination=await folder.getFileHandle('download-copy.fasta',{create:true});
    window.showSaveFilePicker=async()=>destination;
    await downloadGeneratedFile(matrix);
    if(await (await destination.getFile()).text()!==expected)throw Error('Lazy download failed');
    const second=await (await navigator.storage.getDirectory()).getDirectoryHandle('second',{create:true});
    outputDirectoryHandle=second;
    await saveOutputCollection(result.files,true);
    if(outputSaveSession.error)throw outputSaveSession.error;
    const reread=await PhyloPrepSave.readText(matrix);
    return {files:result.files.length,renamed:renamed.length,matrix:reread,allReleased:result.files.every(file=>file.content===undefined)};
  });
  assert.equal(result.renamed,4);assert.ok(result.files>10);assert.ok(result.allReleased);
  console.log(`Real Pyodide disk-backed raw/trimmed QC and paired/independent matrices passed (${result.files} outputs, ${result.renamed} renamed alignments); preview, download and destination-copy read-back passed.`);
 }finally{clearTimeout(timeout);await browser.close();server.close();}
})().catch(error=>{console.error(error);server.close();process.exitCode=1;});
