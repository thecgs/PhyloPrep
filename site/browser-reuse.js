// Reuse completed loci only when inputs, alignment settings, engine and outputs match.
(() => {
  function settings(options) {
    const number=(key,fallback)=>Number(options[key] ?? fallback);
    return {seqtype:options.seqtype || 'codon',align_software:options.align_software || 'mafft',
      genetic_code:number('genetic_code',1),gap_ratio:number('gap_ratio',0.2),
      n_ratio:number('n_ratio',0.2),x_ratio:number('x_ratio',0.2),notrim:Boolean(options.notrim),
      protein_stop_symbol:options.protein_stop_symbol || '*'};
  }
  function engineKey(engine,tool) {
    const names=['MSAP.py','AA2Codon.py','trimAlnSeq.py','sequence_audit.py','msap_io.py'];
    const sources=names.map(name=>[name,engine?.files?.find(file=>file.name===name)?.sha256]);
    if(sources.some(([,sha])=>!sha))return null;
    return JSON.stringify({sources,aligner:engine.aligners?.[tool] || null});
  }
  function create(root,engine,options,{pipeline=true,fingerprint}={}) {
    if(!root || !engine?.bundleSha256 || !fingerprint)return null;
    const configuration=settings(options);
    const alignmentEngine=engineKey(engine,configuration.align_software);
    if(!alignmentEngine)return null;
    const signature=JSON.stringify({configuration,alignmentEngine});
    let legacy;
    const key=input=>'reports/alignment-cache/'+new PhyloPrepHash.SHA256().update(new TextEncoder().encode(input.name)).hex()+'.json';
    async function handle(path,create=false) {
      const parts=path.split('/');
      if(parts.some(part=>!part || part==='.' || part==='..' || part.includes('\\')))throw Error('Invalid cache path');
      let directory=root;
      for(const part of parts.slice(0,-1))directory=await directory.getDirectoryHandle(part,{create});
      return directory.getFileHandle(parts.at(-1),{create});
    }
    function names(input) {
      const leaf=input.name.split('/').at(-1), dot=leaf.lastIndexOf('.');
      const stem=dot>0?leaf.slice(0,dot):leaf, tool=configuration.align_software;
      const types=configuration.seqtype==='codon'?['codon','prot']:[configuration.seqtype];
      const files=[`${stem}.sequence_changes.tsv`];
      for(const type of types)for(const variant of configuration.notrim?['']:['','.trimmed'])files.push(`${stem}.${tool}.${type}${variant}.aln.fasta`);
      return files;
    }
    async function previousRun(input) {
      if(legacy===undefined) {
        try {legacy=JSON.parse(await (await (await handle('phyloprep-run.json')).getFile()).text());}
        catch {legacy=null;}
        if(legacy && (engineKey(legacy.engine_metadata,configuration.align_software)!==alignmentEngine ||
          JSON.stringify(settings(legacy.parameters || {}))!==JSON.stringify(configuration) ||
          JSON.stringify(legacy.engine_metadata?.aligners?.[configuration.align_software] || null)!==JSON.stringify(engine.aligners?.[configuration.align_software] || null)))legacy=null;
        if(legacy){legacy.inputIndex=new Map((legacy.inputs||[]).map(item=>[item.name,item]));legacy.outputIndex=new Map((legacy.outputs||[]).map(item=>[item.name,item]));}
      }
      if(!legacy)return null;
      const source=legacy.inputIndex.get(input.name);
      if(!source)return null;
      const files=names(input).map(name=>{
        const item=legacy.outputIndex.get(pipeline?'alignments/'+name:name);
        return item && {name,path:item.name,size:item.bytes,sha256:item.sha256};
      });
      return files.every(Boolean)?{schema:1,signature,input:{name:source.name,bytes:source.bytes,sha256:source.sha256},files}:null;
    }
    async function restore(input) {
      try {
        let state;
        try {state=JSON.parse(await (await (await handle(key(input))).getFile()).text());}catch {state=await previousRun(input);}
        const source=await fingerprint(input);
        if(state?.schema!==1 || state.signature!==signature || JSON.stringify(state.input)!==JSON.stringify(source))return null;
        const expected=names(input);
        if(state.files?.length!==expected.length || new Set(state.files.map(item=>item.name)).size!==expected.length || state.files.some(item=>!expected.includes(item.name)))return null;
        const files=[];
        for(const item of state.files) {
          const expectedPath=pipeline?'alignments/'+item.name:item.name;
          if(item.path!==expectedPath || !/^[0-9a-f]{64}$/.test(item.sha256))return null;
          const saved=await handle(item.path), snapshot=await saved.getFile();
          if(snapshot.size!==item.size || await PhyloPrepHash.blob(snapshot)!==item.sha256)return null;
          files.push(PhyloPrepSave.record({name:item.name,size:item.size,sha256:item.sha256,
            storage:{disk:{handle:saved,size:snapshot.size,lastModified:snapshot.lastModified}}}));
        }
        return {files,reused:true,engineMetadata:engine,originBundle:state.origin_bundle || legacy?.engine_metadata?.bundleSha256,summary:`Reused verified alignment: ${input.name}`};
      }catch{return null;}
    }
    async function checkpoint(input,result,products) {
      if(result.engineMetadata?.bundleSha256!==engine.bundleSha256)return null;
      const expected=names(input), files=[];
      for(const name of expected) {
        const product=products.find(file=>file.name===name);
        if(!product?.storage?.disk || !product.sha256)return null;
        files.push({name,path:pipeline?'alignments/'+name:name,size:product.size,sha256:product.sha256});
      }
      return {name:key(input),content:JSON.stringify({schema:1,signature,origin_bundle:result.originBundle || result.engineMetadata.bundleSha256,input:await fingerprint(input),files})+'\n'};
    }
    const stage = catalog(root, 'alignment', signature);
    function identity(input) {
      const file=input.storage?.blob || input.blob;
      return Number.isSafeInteger(file?.lastModified) ? {name:input.name,size:file.size,lastModified:file.lastModified} : null;
    }
    async function restoreBatch(inputs) {
      const identities=inputs.map(identity);
      if(identities.some(item=>!item))return null;
      const state=await stage.read();
      if(!stage.matches(state,identities))return null;
      const expected=inputs.flatMap(input=>names(input));
      if(state.files.length!==expected.length || state.files.some((file,index)=>file.name!==(pipeline?'alignments/':'')+expected[index]))return null;
      if(state.inputs?.length!==inputs.length)return null;
      for(let i=0;i<inputs.length;i++){
        if(!state.inputs[i] || state.inputs[i].name!==inputs[i].name || state.inputs[i].bytes!==identities[i].size || !validHash(state.inputs[i].sha256))return null;
      }
      inputs.forEach((input,index)=>Object.assign(input,{size:state.inputs[index].bytes,sha256:state.inputs[index].sha256}));
      return {files:stage.records(state).map(file=>{file.name=file.name.split('/').at(-1);return file;}),engineMetadata:engine,reused:true};
    }
    async function checkpointBatch(inputs,products) {
      const identities=inputs.map(identity);
      if(identities.some(item=>!item))return null;
      const index=new Map(products.map(file=>[file.name,file]));
      const ordered=inputs.flatMap(input=>names(input)).map(name=>index.get(name));
      if(ordered.some(file=>!file))return null;
      const files=ordered.map(file=>({...file,name:(pipeline?'alignments/':'')+file.name}));
      const result=stage.checkpoint(identities,files);
      if(!result)return null;
      const state=JSON.parse(result.content);state.inputs=await Promise.all(inputs.map(fingerprint));
      result.content=JSON.stringify(state)+'\n';return result;
    }
    return {restore,checkpoint,restoreBatch,checkpointBatch,clearBatch:stage.clear};
  }
  const validHash=value=>/^[0-9a-f]{64}$/.test(value || '');
  // Stage manifests contain metadata only. Open saved data lazily when it is used.
  function catalog(root, name, signature) {
    const path='reports/stage-checkpoints/'+name+'.json';
    async function read() {
      try {return JSON.parse(await (await (await PhyloPrepSave.diskHandle({root,path})).getFile()).text());}
      catch {return null;}
    }
    function matches(state, sources) {
      return state?.schema===1 && state.complete===true && state.signature===signature &&
        JSON.stringify(state.sources)===JSON.stringify(sources) && Array.isArray(state.files) &&
        state.files.every(file=>typeof file?.name==='string' && validHash(file.sha256) && Number.isSafeInteger(file.size) && file.size>=0 &&
          !file.name.split('/').some(part=>!part || part==='.' || part==='..' || part.includes('\\'))) &&
        new Set(state.files.map(file=>file.name)).size===state.files.length;
    }
    function records(state) {
      return state.files.map(file=>PhyloPrepSave.record({...file,storage:{disk:{root,path:file.name,size:file.size,lastModified:file.lastModified}}}));
    }
    function checkpoint(sources, files) {
      if(files.some(file=>!file.storage?.disk || !validHash(file.sha256) || !Number.isSafeInteger(file.size) || file.size<0))return null;
      return {name:path,content:JSON.stringify({schema:1,complete:true,signature,sources,
        files:files.map(file=>({name:file.name,size:file.size,sha256:file.sha256,lastModified:file.storage.disk.lastModified}))})+'\n'};
    }
    async function clear() {
      try {const reports=await root.getDirectoryHandle('reports');const directory=await reports.getDirectoryHandle('stage-checkpoints');await directory.removeEntry(name+'.json');}
      catch(error){if(error.name!=='NotFoundError')throw error;}
    }
    return {read,matches,records,checkpoint,clear};
  }
  function createRename(root,engine,mapping) {
    const sources=['rename_taxa.py','msap_io.py'].map(name=>[name,engine?.files?.find(file=>file.name===name)?.sha256]);
    if(!root || !mapping || sources.some(([,sha])=>!sha))return null;
    const signature=JSON.stringify({sources,mapping});
    let previous;
    async function previousRun(){
      if(previous===undefined)previous=(async()=>{
        try{
          const run=JSON.parse(await (await (await PhyloPrepSave.diskHandle({root,path:'phyloprep-run.json'})).getFile()).text());
          const input=run.inputs?.find(item=>item.input_key==='mapping');
          if(!input || input.name!==mapping.name || input.bytes!==mapping.bytes || input.sha256!==mapping.sha256 ||
            sources.some(([name,sha])=>run.engine_metadata?.files?.find(file=>file.name===name)?.sha256!==sha))return null;
          return new Map((run.outputs || []).map(file=>[file.name,file]));
        }catch{return null;}
      })();
      return previous;
    }
    return variant=>{
      const stage=catalog(root,'rename-'+variant,signature);
      const expected=files=>files.map(file=>({name:file.name,size:file.size,sha256:file.sha256}));
      // Derive the alphabet from the stable .TYPE[.trimmed].aln.fasta suffix.
      const targetName=file=>{const match=file.name.match(/\.(codon|prot|nucl)(?:\.trimmed)?\.aln\.fasta$/);return match && `renamed/${match[1]}/${variant}/${file.name}`;};
      async function restore(files) {
        let state=await stage.read();
        const wanted=expected(files);
        let legacy=false;
        if(!state){
          const run=await previousRun();
          if(run){
            legacy=true;
            const available=files.filter(file=>run.has('alignments/'+file.name) && run.has(targetName(file)));
            state={schema:1,complete:true,signature,
              sources:available.map(file=>{const item=run.get('alignments/'+file.name);return {name:file.name,size:item.bytes,sha256:item.sha256};}),
              files:available.map(file=>{const item=run.get(targetName(file));return {name:item.name,size:item.bytes,sha256:item.sha256};})};
          }
        }
        if(!legacy && stage.matches(state,wanted) && state.files.length===files.length && state.files.every((file,i)=>file.name===targetName(files[i])))return {files:stage.records(state),complete:true};
        // A different locus set can still reuse matching individual entries.
        if(!state || !Array.isArray(state.sources) || !stage.matches(state,state.sources))return {files:[],complete:false};
        const saved=new Map(state.files.map(file=>[file.name,file])), previous=new Map(state.sources.map(file=>[file.name,file]));
        const restored=[];
        for(const file of files){
          const old=previous.get(file.name), item=saved.get(targetName(file));
          if(JSON.stringify(old)!==JSON.stringify(expected([file])[0]) || !item || !validHash(item.sha256))continue;
          try{
            const record=stage.records({files:[item]})[0], blob=await PhyloPrepSave.readBlob(record);
            if(await PhyloPrepHash.blob(blob)===item.sha256)restored.push(record);
          }catch {} // Missing or changed files are renamed again.
        }
        return {files:restored,complete:false};
      }
      return {restore,clear:stage.clear,checkpoint:(files,outputs)=>{
        const index=typeof outputs.get==='function'?outputs:new Map(outputs.map(file=>[file.name,file]));
        const selected=files.map(file=>index.get(targetName(file)));
        return selected.every(Boolean)?stage.checkpoint(expected(files),selected):null;
      }};
    };
  }
  globalThis.PhyloPrepReuse={create,settings,createRename};
})();
