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
    return {restore,checkpoint};
  }
  globalThis.PhyloPrepReuse={create,settings};
})();
