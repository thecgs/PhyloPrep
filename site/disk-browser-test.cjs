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
  await page.evaluate(async()=>{
    const helper=new URL('browser-disk-fs.js',location.href).href;
    const script=`importScripts(${JSON.stringify(helper)});
      let cache;
      const runtime={FS:{filesystems:{MEMFS:{}},mkdirTree(){},mount(){},unmount(){},rmdir(){}},ERRNO_CODES:{}};
      onmessage=async event=>{try {
        if(event.data==='open') cache=await PhyloPrepDiskFS.open(runtime,'/work/test');
        else await cache.close();
        postMessage({ok:true});
      }catch(error){postMessage({error:String(error)});}};`;
    const url=URL.createObjectURL(new Blob([script],{type:'application/javascript'}));
    const workers=[new Worker(url),new Worker(url),new Worker(url)];
    const ask=(worker,command)=>new Promise((resolve,reject)=>{
      const timer=setTimeout(()=>reject(Error('Cache lifecycle test timed out')),10000);
      worker.onmessage=event=>{clearTimeout(timer);event.data.error?reject(Error(event.data.error)):resolve();};
      worker.onerror=event=>{clearTimeout(timer);reject(Error(event.message));};worker.postMessage(command);
    });
    const directory=await (await navigator.storage.getDirectory()).getDirectoryHandle('phyloprep-working-cache',{create:true});
    const count=async()=>{let total=0;for await(const [name] of directory.entries())if(name.startsWith('workspace-'))total++;return total;};
    try{
      await ask(workers[0],'open');await ask(workers[1],'open');
      if(await count()!==2)throw Error('An active worker cache was deleted');
      workers[0].terminate();await new Promise(resolve=>setTimeout(resolve,100));
      await ask(workers[2],'open');
      if(await count()!==2)throw Error('An abandoned cache was not reclaimed');
      await ask(workers[1],'close');await ask(workers[2],'close');
      if(await count()!==0)throw Error('Finished caches were not removed');
    }finally{workers.forEach(worker=>worker.terminate());URL.revokeObjectURL(url);}
  });
  console.log('Real OPFS active-worker isolation, cancelled-cache recovery and completed-cache cleanup passed.');
  const quotaProbe = await page.evaluate(async()=>{
    const helper=new URL('browser-disk-fs.js',location.href).href;
    const script=`importScripts(${JSON.stringify(helper)});
    onmessage=async()=>{
      const root=await navigator.storage.getDirectory();
      const name='quota-probe-'+crypto.randomUUID();
      const handle=await(await root.getFileHandle(name,{create:true})).createSyncAccessHandle();
      let result;
      try {
        const codes={ENOSPC:51,EIO:29,EINVAL:28};
        const FS={isDir:mode=>mode===0x4000,isFile:mode=>mode===0x8000,
          ErrnoError:class extends Error{constructor(errno){super(String(errno));this.errno=errno;}},filesystems:{}};
        function createNode(parent,name,mode){
          const node={parent,name,mode,usedBytes:0,contents:mode===0x4000?{}:null,
            node_ops:{setattr:Object.assign,unlink(parent,name){delete parent.contents[name];}},stream_ops:{}};
          if(parent)parent.contents[name]=node;return node;
        }
        FS.filesystems.MEMFS={createNode,mount:()=>createNode(null,'/',0x4000)};
        let forceQuota=true,returned;
        const quota=(await navigator.storage.estimate()).quota;
        const access={getSize:()=>handle.getSize(),truncate:size=>handle.truncate(size),flush:()=>handle.flush(),
          read:(view,options)=>handle.read(view,options),
          write(view,{at}){returned=handle.write(view,{at:forceQuota?quota+65536+at:at});return returned;}};
        const backend=PhyloPrepDiskFS.filesystem(FS,access,codes),directory=backend.mount({});
        const node=directory.node_ops.mknod(directory,'fourfold-spool',0x8000);
        let rejected=false;
        try{node.stream_ops.write({node},new Uint8Array([7]),0,1,0);}
        catch(error){rejected=[codes.ENOSPC,codes.EIO].includes(error.errno);}
        if(!rejected||node.usedBytes||node.diskBlocks.size)throw Error('Quota failure created a phantom data block');
        const failure=backend.failure;
        forceQuota=false;
        const bytes=new Uint8Array(100000).fill(42);
        node.stream_ops.open({node});node.stream_ops.dup({node});
        node.stream_ops.write({node},bytes,0,bytes.length,0);
        directory.node_ops.unlink(directory,'fourfold-spool');node.stream_ops.close({node});
        const copy=new Uint8Array(bytes.length);
        node.stream_ops.read({node},copy,0,copy.length,0);
        if(!copy.every(value=>value===42))throw Error('Open duplicate lost spool data');
        node.stream_ops.close({node});
        if(handle.getSize()!==0)throw Error('Closed spool did not return backing-file quota');
        result={ok:true,failure,quota};
      }catch(error){result={error:String(error)};}
      finally{handle.close();await root.removeEntry(name);}
      postMessage(result);
    };`;
    const url=URL.createObjectURL(new Blob([script],{type:'application/javascript'}));
    const worker=new Worker(url);
    try{return await new Promise((resolve,reject)=>{
      worker.onmessage=event=>resolve(event.data);worker.onerror=event=>reject(Error(event.message));worker.postMessage('start');
    });}finally{worker.terminate();URL.revokeObjectURL(url);}
  });
  assert.equal(quotaProbe.ok,true,JSON.stringify(quotaProbe));
  console.log('Real OPFS quota failure rejected before read-back; duplicate spool lifetime and physical tail-quota reclamation passed.');
  const result=await page.evaluate(async()=>{
    const folder=await (await navigator.storage.getDirectory()).getDirectoryHandle('outputs',{create:true});
    outputDirectoryHandle=folder; outputSaveSession=createOutputSaveSession(folder);producedFiles=new Map();
    // Concatenation crosses the bounded-reader block boundary (>65,535 bases).
    const unit='ATGGCT'.repeat(2000), amino='MA'.repeat(2000);
    const dna=`>gene_a\n${unit}\n>gene_b\n${unit}\n`, protein=`>gene_a\n${amino}\n>gene_b\n${amino}\n`;
    const batch=[], sources=[];
    const loci=6;
    for(let locus=0;locus<loci;locus++) {
    sources.push({name:`gene${locus}.fa`,role:'source',content:dna});
    for(const variant of ['', '.trimmed']) for(const type of ['codon','prot']){
      const name=`gene${locus}.mafft.${type}${variant}.aln.fasta`;
      const file=PhyloPrepSave.record({name:`alignments/${name}`,content:type==='codon'?dna:protein});
      await fileFingerprint(file);await publishProducedFiles([file]);await outputSaveSession.settle();
      if(file.content!==undefined||!file.storage.disk)throw Error('Saved alignment text was retained');
      batch.push(PhyloPrepSave.record({...file,name,role:'batch-output'}));
    }
    }
    const result=await runInWorker({type:'phyloprep-finalize',script:'phyloprep.py',files:[
      ...sources,
      {name:'map.tsv',role:'mapping',content:'gene_a\tsp_a\ngene_b\tsp_b\n'},...batch],
      options:{seqtype:'codon',genetic_code:1},outputName:'',usesPathList:false});
    await outputSaveSession.settle();
    if(outputSaveSession.error)throw outputSaveSession.error;
    const html=result.files.find(file=>file.name==='phyloprep-report.html');
    if(result.report?.status!=='generated'||!html)throw Error('Automatic HTML report failed: '+JSON.stringify(result.warnings));
    const reportText=await PhyloPrepSave.readText(html);
    const reportData=JSON.parse(reportText.match(/id="report-data">(.*?)<\/script>/s)[1]);
    if(reportData.overview.input_loci!==loci||reportData.matrices.length!==24)throw Error('Automatic report metadata incomplete');
    if(result.files.some(file=>file.content!==undefined||!file.storage.disk))throw Error('Final result text was retained');
    if(result.workspace_storage?.kind!=='opfs'||result.workspace_storage.cached_file_bytes!==0)throw Error('Locus files were not disk-backed');
    const renamed=result.files.filter(file=>file.name.startsWith('renamed/'));
    const matrix=result.files.find(file=>file.name==='matrices/codon/raw/supermatrix.fasta');
    const paired=result.files.find(file=>file.name==='matrices/paired/codon/raw/supermatrix.fasta');
    const expected=`>sp_a\n${unit.repeat(loci)}\n>sp_b\n${unit.repeat(loci)}\n`;
    if(await PhyloPrepSave.readText(matrix)!==expected||await PhyloPrepSave.readText(paired)!==expected)throw Error('Matrix changed');
    const fourfold=result.files.find(file=>file.name==='matrices/codon/raw/fourfold.fasta');
    const third=result.files.find(file=>file.name==='matrices/codon/raw/codon3rd.supermatrix.fasta');
    if(await PhyloPrepSave.readText(fourfold)!==`>sp_a\n${'T'.repeat(2000*loci)}\n>sp_b\n${'T'.repeat(2000*loci)}\n`)throw Error('Fourfold matrix changed across blocks');
    if(await PhyloPrepSave.readText(third)!==`>sp_a\n${'GT'.repeat(2000*loci)}\n>sp_b\n${'GT'.repeat(2000*loci)}\n`)throw Error('Codon-position matrix changed across blocks');
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
    return {files:result.files.length,renamed:renamed.length,matrix:reread,allReleased:result.files.every(file=>file.content===undefined),storage:result.workspace_storage,loci};
  });
  assert.equal(result.renamed,4*result.loci);assert.ok(result.files>10);assert.ok(result.allReleased);
  console.log(`Real Pyodide disk-backed raw/trimmed QC and paired/independent matrices passed (${result.files} outputs, ${result.renamed} renamed alignments); preview, download and destination-copy read-back passed.`);
 const extra=await page.evaluate(async()=>{
   const root=await (await navigator.storage.getDirectory()).getDirectoryHandle('reuse-tests',{create:true});
   outputDirectoryHandle=root;outputSaveSession=createOutputSaveSession(root);producedFiles=new Map();
   const engine=await (await fetch('engine/manifest.json')).json();
   const original=runMsapBrowserWorkflow;let calls=0;
   runMsapBrowserWorkflow=async()=>{
     calls++;
     const files=[{name:'cached.mafft.codon.aln.fasta',content:'>g1\nATGGCT\n>g2\nATGGCT\n'},
       {name:'cached.mafft.prot.aln.fasta',content:'>g1\nMA\n>g2\nMA\n'},
       {name:'cached.sequence_changes.tsv',content:'header\n'}];
     for(const file of files)await fileFingerprint(file);
     return {files,engineMetadata:engine};
   };
   try {
     const request={pipelineOutput:true,options:{seqtype:'codon',notrim:true},files:[
       {name:'cached.fa',inputKey:'inputs',fasta:true,blob:new File(['>g1\nATGGCT\n>g2\nATGGCT\n'],'cached.fa')}]};
     await runMsapBatchBrowserWorkflow(request);await outputSaveSession.settle();
     const alignments=await root.getDirectoryHandle('alignments');
     const saved=await alignments.getFileHandle('cached.mafft.codon.aln.fasta');
     const before=(await saved.getFile()).lastModified;
     outputSaveSession=createOutputSaveSession(root);producedFiles=new Map();
     const reused=await runMsapBatchBrowserWorkflow({...request,options:{...request.options,min_nucl_pis:5}});
     await outputSaveSession.settle();
     if(calls!==1 || reused.reusedLoci!==1)throw Error('Matching locus was recomputed');
     if((await saved.getFile()).lastModified!==before)throw Error('Reused alignment was unnecessarily rewritten');
     outputSaveSession=createOutputSaveSession(root);
     await runMsapBatchBrowserWorkflow({...request,options:{...request.options,no_resume:true}});await outputSaveSession.settle();
     if(calls!==2)throw Error('--no-resume did not recompute');
     const stream=await saved.createWritable();await stream.write('>g1\nATGGCC\n>g2\nATGGCC\n');await stream.close();
     outputSaveSession=createOutputSaveSession(root);
     await runMsapBatchBrowserWorkflow(request);await outputSaveSession.settle();
     if(calls!==3)throw Error('Changed saved alignment was incorrectly reused');
   }finally{runMsapBrowserWorkflow=original;}
   outputSaveSession=createOutputSaveSession(root);producedFiles=new Map();
   const sequence='GCT'.repeat(220000), input=new File([`>a\n${sequence}\n>b\n${sequence}\n`],'large.fa');
   const converted=await runInWorker({script:'fasta2phy.py',files:[{name:input.name,blob:input,fasta:true}],
     options:{seqtype:'DNA'},outputName:'large.phy',usesPathList:false});
   await outputSaveSession.settle();
   const matrix=converted.files.find(file=>file.name==='large.phy');
   if(matrix.content!==undefined || matrix.storage.blob || !matrix.storage.disk)throw Error('Large binary output was retained');
   if((await PhyloPrepSave.readBlob(matrix)).size<1300000)throw Error('Large output truncated');
   const extracted=await runInWorker({type:'orthofinder-extract',script:'extract_orthofinder_orthogroups.py',
     files:[{name:'Orthogroups.tsv',inputKey:'orthogroups',blob:new File(['Orthogroup\tsp_a\tsp_b\nOG1\tg1\tg2\n'],'Orthogroups.tsv'),fasta:false},
       ...['sp_a','sp_b'].map((name,i)=>({name:name+'.fa',inputKey:'sequences',blob:new File([`>g${i+1}\nATGGCT\n`],name+'.fa'),fasta:true}))],
     options:{copy_mode:'single-copy',orthofinder_id_style:'unprefixed'},outputName:'ogs',usesPathList:false});
   await outputSaveSession.settle();
   const og=extracted.files.find(file=>/OG1\.(fa|fasta)$/.test(file.name));
   if(!og || !(await PhyloPrepSave.readText(og)).includes('>sp_a_g1'))throw Error('Blob-backed OG input staging failed');
   const e2e=await (await navigator.storage.getDirectory()).getDirectoryHandle('pipeline-reuse',{create:true});
   outputDirectoryHandle=e2e;outputSaveSession=createOutputSaveSession(e2e);producedFiles=new Map();
   const align=runBrowserAlignment;let aligned=0;
   // Identity alignment for identical sequences; original MSAP/Python stages remain real.
   runBrowserAlignment=async sequence=>{aligned++;return sequence;};
   try{
     const request={options:{seqtype:'codon',notrim:true},files:[
       ...[0,1].map(i=>({name:`e2e${i}.fa`,inputKey:'inputs',fasta:true,blob:new File([i===0?'>g1\nATGGCT\n>g2\nATGGCT\n':'>g1\nATGGCT\n>g3\nATGGCT\n'],`e2e${i}.fa`)})),
       {name:'mapping.tsv',inputKey:'mapping',fasta:false,blob:new File(['g1\tsp_a\ng2\tsp_b\ng3\tsp_c\n'],'mapping.tsv')}]};
     await runPhyloPrepBrowserWorkflow(request);await outputSaveSession.settle();
     if(aligned!==2)throw Error('End-to-end alignment did not run');
     outputSaveSession=createOutputSaveSession(e2e);producedFiles=new Map();
     const rerun=await runPhyloPrepBrowserWorkflow(request);await outputSaveSession.settle();
     if(rerun.report?.status!=='generated'||!rerun.files.some(file=>file.name==='phyloprep-report.html'))throw Error('Cached pipeline did not regenerate report');
     if(aligned!==2 || rerun.reused_loci!==2)throw Error('End-to-end cached pipeline recomputed alignment');
     const matrix=rerun.files.find(file=>file.name==='matrices/codon/raw/supermatrix.fasta');
     if(await PhyloPrepSave.readText(matrix)!=='>sp_a\nATGGCTATGGCT\n>sp_b\nATGGCT??????\n>sp_c\n??????ATGGCT\n')throw Error('Cached pipeline matrix changed');
     for(const extension of ['phy','nex']){
       const file=rerun.files.find(file=>file.name===`matrices/codon/raw/supermatrix.${extension}`);
       if(!(await PhyloPrepSave.readText(file)).includes('??????'))throw Error('Padding marker lost during export');
     }
     const missing=rerun.files.find(file=>file.name==='matrices/codon/raw/supermatrix.missing.tsv');
     if(!(await PhyloPrepSave.readText(missing)).includes('e2e1\tsp_b\t7\t12\t6\t?'))throw Error('Missing block provenance lost');
     const fourfold=rerun.files.find(file=>file.name==='matrices/codon/raw/fourfold.fasta');
     if(await PhyloPrepSave.readText(fourfold)!=='>sp_a\nTT\n>sp_b\nT?\n>sp_c\n?T\n')throw Error('Default observed-taxon fourfold selection failed');
     const html=rerun.files.find(file=>file.name==='phyloprep-report.html');
     const report=JSON.parse((await PhyloPrepSave.readText(html)).match(/id="report-data">(.*?)<\/script>/s)[1]);
     if(report.provenance.parameters.fourfold_missing_policy!=='ignore-missing')throw Error('Fourfold rule missing from report provenance');
     const strictFolder=await e2e.getDirectoryHandle('strict-fourfold',{create:true});
     outputDirectoryHandle=strictFolder;outputSaveSession=createOutputSaveSession(strictFolder);producedFiles=new Map();
     const strict=await runPhyloPrepBrowserWorkflow({...request,options:{...request.options,fourfold_missing_policy:'require-all'}});
     await outputSaveSession.settle();
     const skipped=strict.files.find(file=>file.name==='matrices/codon/raw/fourfold.skipped.txt');
     if(!skipped || !(await PhyloPrepSave.readText(skipped)).includes('all-taxon rule'))throw Error('Strict fourfold strategy did not skip an empty dataset');
     if(strict.files.some(file=>file.name==='matrices/codon/raw/fourfold.fasta'))throw Error('Strict selection published empty fourfold FASTA');
   }finally{runBrowserAlignment=align;releaseIdleWorkers({pyodide:true});}
   return {calls,size:matrix.size};
 });
 console.log(`Real browser verified reuse, --no-resume, changed-file recomputation, no-rewrite adoption, Blob input/output (${extra.size} bytes) and OG extraction passed.`);
 await context.route('**/sqlite3-*.whl',route=>route.abort('failed'));
 const reportFailure=await page.evaluate(async()=>{
   const run=async no_report=>{
     releaseIdleWorkers({pyodide:true});
     outputDirectoryHandle=undefined;outputSaveSession=undefined;producedFiles=new Map();
     return runInWorker({type:'phyloprep-finalize',script:'phyloprep.py',files:[
       {name:'report-failure.fa',role:'source',content:'>a\nACGT\n>b\nACGA\n'},
       {name:'report-failure.mafft.nucl.aln.fasta',role:'batch-output',content:'>a\nACGT\n>b\nACGA\n'}],
       options:{seqtype:'nucl',notrim:true,no_report},outputName:'',usesPathList:false});
   };
   try{
     const failed=await run(false);
     if(!failed.ok||failed.report?.status!=='failed'||!failed.warnings?.some(w=>w.includes('sqlite3')))throw Error('Report package failure affected analysis');
     if(!failed.files.some(file=>file.name==='matrices/nucl/raw/supermatrix.fasta')||failed.files.some(file=>file.name==='phyloprep-report.html'))throw Error('Report failure lost matrices or returned stale HTML');
     const skipped=await run(true);
     if(!skipped.ok||skipped.report?.status!=='skipped'||skipped.warnings.length)throw Error('--no-report attempted report package loading');
     return {failed:failed.report.status,skipped:skipped.report.status};
   }finally{releaseIdleWorkers({pyodide:true});}
 });
 console.log('Real browser report package failure retained matrices; --no-report skipped package loading: '+JSON.stringify(reportFailure));
 }finally{clearTimeout(timeout);await browser.close();server.close();}
})().catch(error=>{console.error(error);server.close();process.exitCode=1;});
