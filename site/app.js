const inputFolderHandles = new Map();
const inputFolderFiles = new Map();
const resolvedFolderFiles = new Map();
let preflightGeneration = 0;
let pathlistNoticeGeneration = 0;
const elements = {
  runtimeStatus: document.querySelector("#runtime-status"),
  scriptList: document.querySelector("#script-list"),
  selectedName: document.querySelector("#selected-name"),
  selectedSummary: document.querySelector("#selected-summary"),
  parameterList: document.querySelector("#parameter-list"),
  workspaceScript: document.querySelector("#workspace-script"),
  workspaceMode: document.querySelector("#workspace-mode"),
  chooseFolder: document.querySelector("#choose-output-folder"),
  outputFolder: document.querySelector("#output-folder-display"),
  inputFileFields: document.querySelector("#input-file-fields"),
  inputFilesList: document.querySelector("#input-files-list"),
  clearInputFiles: document.querySelector("#clear-input-files"),
  inputFilesHelp: document.querySelector("#input-files-help"),
  outputNameLabel: document.querySelector("#output-name-label"),
  outputName: document.querySelector("#output-name"),
  outputNameHelp: document.querySelector("#output-name-help"),
  outputNameSection: document.querySelector("#output-name-section"),
  run: document.querySelector("#run-browser-module"),
  cancelRun: document.querySelector("#cancel-run"),
  log: document.querySelector("#workspace-log"),
  implementationNote: document.querySelector("#implementation-note"),
  outputFilesList: document.querySelector("#output-files-list"),
  resultCount: document.querySelector("#result-count"),
  resultControls: document.querySelector("#result-controls"),
  resultFolder: document.querySelector("#result-folder"),
  resultSearch: document.querySelector("#result-search"),
  resultPrev: document.querySelector("#result-prev"),
  resultNext: document.querySelector("#result-next"),
  resultPage: document.querySelector("#result-page"),
  downloadSelected: document.querySelector("#download-selected"),
  downloadAll: document.querySelector("#download-all"),
  preview: document.querySelector("#output-preview"),
  command: document.querySelector("#command-example"),
  copyCommand: document.querySelector("#copy-command"),
  runStages: document.querySelector("#run-stages"),
  runStageLabel: document.querySelector("#run-stage-label"),
  scriptSearch: document.querySelector("#script-search"),
  preflightStatus: document.querySelector("#preflight-status"),
  preflightSummary: document.querySelector("#preflight-summary"),
  useNextStep: document.querySelector("#use-next-step"),
};

const geneticCodeOptions = [
  ["1", "1. The Standard Code"],
  ["2", "2. The Vertebrate Mitochondrial Code"],
  ["3", "3. The Yeast Mitochondrial Code"],
  ["4", "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code"],
  ["5", "5. The Invertebrate Mitochondrial Code"],
  ["6", "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code"],
  ["9", "9. The Echinoderm and Flatworm Mitochondrial Code"],
  ["10", "10. The Euplotid Nuclear Code"],
  ["11", "11. The Bacterial, Archaeal and Plant Plastid Code"],
  ["12", "12. The Alternative Yeast Nuclear Code"],
  ["13", "13. The Ascidian Mitochondrial Code"],
  ["14", "14. The Alternative Flatworm Mitochondrial Code"],
  ["15", "15. Blepharisma Nuclear Code"],
  ["16", "16. Chlorophycean Mitochondrial Code"],
  ["21", "21. Trematode Mitochondrial Code"],
  ["22", "22. Scenedesmus obliquus Mitochondrial Code"],
  ["23", "23. Thraustochytrium Mitochondrial Code"],
  ["24", "24. Rhabdopleuridae Mitochondrial Code"],
  ["25", "25. Candidate Division SR1 and Gracilibacteria Code"],
  ["26", "26. Pachysolen tannophilus Nuclear Code"],
  ["27", "27. Karyorelict Nuclear Code"],
  ["28", "28. Condylostoma Nuclear Code"],
  ["29", "29. Mesodinium Nuclear Code"],
  ["30", "30. Peritrich Nuclear Code"],
  ["31", "31. Blastocrithidia Nuclear Code"],
  ["32", "32. Balanophoraceae Plastid Code"],
  ["33", "33. Cephalodiscidae Mitochondrial UAA-Tyr Code"],
];

function geneticCodeParameter() {
  return {
    flag: "-g TABLE",
    label: "Genetic-code table",
    description: "Valid NCBI genetic-code table ID. The browser accepts only codes in the DDBJ/NCBI list.",
    type: "select",
    value: "1",
    options: geneticCodeOptions,
    token: "genetic_code",
    reference: "https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes",
  };
}

function proteinStopParameter() {
  return {
    flag: "--protein-stop-symbol",
    label: "Protein stop symbol",
    description: "Protein mode only. Trailing occurrences are removed and internal occurrences become X. Default: *.",
    type: "text",
    value: "*",
    token: "protein_stop_symbol",
  };
}

const modules = [
  {
    script: "phyloprep.py", title: "Complete phylogenomic preparation pipeline", category: "Pipeline", browser: true, browserWorkflow: "phyloprep", inputListKey: "inputs", cliOutputDir: "phyloprep-results",
    summary: "Aligns multiple loci with browser MAFFT or MUSCLE, applies PhyloPrep QC and optional taxon renaming, builds supermatrices, and exports downstream matrices.",
    inputs: ["One or more input FASTA files", "Optional two-column taxon mapping table"],
    inputFields: [
      { key: "inputs", label: "Input FASTA files", flag: "-i, --input", multiple: true, pathList: true },
      { key: "mapping", label: "Taxon mapping table", flag: "-m, --mapping", accept: ".tsv,.txt,.csv", required: false },
    ],
    outputs: ["Per-locus alignments and audit reports", "QC reports", "Optional renamed alignments", "Independent and paired raw/trimmed supermatrices", "Codon-position and four-fold matrices", "FASTA/PHYLIP/NEXUS exports"],
    outputLabel: "Pipeline output", outputHelp: "Browser output is collected as downloadable files or written into the selected local folder.", defaultOutput: "", hideOutputName: true,
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Browser pipeline supports codon, protein, and nucleotide modes. Pseudogene/MACSE remains CLI-only.", type: "select", value: "codon", options: [["codon", "Codon / CDS"], ["prot", "Protein"], ["nucl", "Nucleotide"]], token: "seqtype" },
      geneticCodeParameter(),
      { flag: "--missing-taxa", label: "Missing taxa policy", description: "Skip genes missing taxa, or pad missing taxa with gaps during concatenation.", type: "select", value: "skip-gene", options: [["skip-gene", "skip-gene"], ["pad-gaps", "pad-gaps"]], token: "missing_taxa" },
      { flag: "-s", label: "Alignment software", description: "Runs locally through Biowasm/WebAssembly. MAFFT E-INS-i is the browser default through Biowasm's compiled tbfast + dvtditr pipeline; CLI -s mafft continues to use upstream mafft --auto and can differ.", type: "select", value: "mafft", options: [["mafft", "MAFFT E-INS-i 7.520 (Biowasm, default)"], ["muscle", "MUSCLE 5.1.0"]], token: "align_software" },
      { flag: "-t", label: "Threads per locus", description: "Browser fixed to one thread for GitHub Pages/WebAssembly compatibility.", type: "select", value: "1", options: [["1", "1 (browser fixed)"]], token: "thread" },
      { flag: "-j", label: "Concurrent loci", description: "Browser edition processes loci sequentially to keep memory predictable.", type: "select", value: "1", options: [["1", "1 (browser sequential)"]], token: "jobs" },
      { flag: "-ts", label: "Trimming software", description: "Browser fixed to PhyloPrep's bundled trimAlnSeq.py.", type: "select", value: "trimAlnSeq", options: [["trimAlnSeq", "trimAlnSeq.py (browser fixed)"]], token: "trim_software" },
      { flag: "-n", label: "Skip trimming", description: "Process raw alignments only.", type: "checkbox", value: false, token: "notrim" },
      proteinStopParameter(),
      { flag: "-G", label: "Maximum gap proportion", description: "Maximum gap ratio per site during trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "gap_ratio" },
      { flag: "-N", label: "Maximum N proportion", description: "Maximum N ratio per site during nucleotide/codon trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "n_ratio" },
      { flag: "-X", label: "Maximum X proportion", description: "Maximum X ratio per site during protein trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "x_ratio" },
      { flag: "--min-length", label: "Minimum QC length", description: "Minimum sequence/alignment length; codon mode uses complete codons.", type: "number", value: "1", token: "min_length" },
      { flag: "--max-n-ratio", label: "Maximum QC N ratio", description: "Maximum N-containing-codon ratio for codon QC (or N-base ratio for nucleotide QC). Codon/protein loci are retained as pairs.", type: "number", value: "1", step: "0.01", token: "max_n_ratio" },
      { flag: "--max-x-ratio", label: "Maximum QC X ratio", description: "Maximum X ratio for protein QC. In codon workflow the final retained locus set is synchronized with codon QC.", type: "number", value: "1", step: "0.01", token: "max_x_ratio" },
      { flag: "--max-gap-ratio", label: "Maximum QC gap ratio", description: "Maximum gap ratio allowed by alignment QC. Range: 0–1.", type: "number", value: "1", step: "0.01", token: "max_gap_ratio" },
      { flag: "--min-nucl-pis", label: "Minimum nucleotide PIS", description: "Minimum nucleotide-column PIS per locus in nucleotide/codon modes. Default: 0.", type: "number", value: "0", token: "min_nucl_pis" },
      { flag: "--min-protein-pis", label: "Minimum protein PIS", description: "Minimum amino-acid PIS per locus in protein/codon modes. Default: 0.", type: "number", value: "0", token: "min_protein_pis" },
      { flag: "--paired-pis-policy", label: "Paired PIS policy", description: "Codon workflow only. Decide how the minimum PIS threshold is applied to the paired codon/protein locus: both must pass, either may pass, or only codon/protein PIS controls retention. Final codon and protein matrices still use the same locus set.", type: "select", value: "both", options: [["both", "both — codon AND protein must pass"], ["either", "either — codon OR protein may pass"], ["codon", "codon — use codon PIS only"], ["protein", "protein — use protein PIS only"]], token: "paired_pis_policy" },
      { flag: "--allow-internal-stop", label: "Allow internal stop", description: "Retain a protein locus or codon/protein locus pair when an internal stop is observed directly or recorded by MSAP as masked to NNN/X. The masked NNN/X is never restored.", type: "checkbox", value: false, token: "allow_internal_stop" },
    ],
    command: "python phyloprep.py -i {inputs} -o phyloprep-results{mapping_option} -st {seqtype} -s {align_software} -t 1 -j 1 -g {genetic_code} -ts trimAlnSeq -G {gap_ratio} -N {n_ratio} -X {x_ratio} --missing-taxa {missing_taxa} --min-length {min_length} --max-n-ratio {max_n_ratio} --max-x-ratio {max_x_ratio} --max-gap-ratio {max_gap_ratio}{min_nucl_pis_option}{min_protein_pis_option} --paired-pis-policy {paired_pis_policy} --protein-stop-symbol {protein_stop_symbol_cli}{allow_internal_stop}{notrim_option}",
  },
  {
    script: "MSAP_batch.py", title: "Batch alignment workflow", category: "Batch", browser: true, browserWorkflow: "batch", inputListKey: "inputs", cliOutputDir: "msap-results",
    summary: "Runs the MSAP workflow over multiple FASTA loci in the browser using Biowasm MAFFT or MUSCLE and aggregates path lists and sequence-change reports.",
    inputs: ["One or more input FASTA files"],
    inputFields: [{ key: "inputs", label: "Input FASTA files", flag: "-i, --input", multiple: true, pathList: true }],
    outputs: ["Per-locus preprocessing and audit files", "MAFFT or MUSCLE alignments", "Optional trimmed alignments", "Batch path lists", "Merged sequence-change report"],
    outputLabel: "Batch output", outputHelp: "Browser output is collected as downloadable files or written into the selected local folder.", defaultOutput: "", hideOutputName: true,
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Browser batch supports codon, protein, and nucleotide modes. Pseudogene/MACSE remains CLI-only.", type: "select", value: "codon", options: [["codon", "Codon / CDS"], ["prot", "Protein"], ["nucl", "Nucleotide"]], token: "seqtype" },
      geneticCodeParameter(),
      { flag: "-s", label: "Alignment software", description: "Runs locally through Biowasm/WebAssembly. MAFFT E-INS-i is the browser default through Biowasm's compiled tbfast + dvtditr pipeline; CLI -s mafft continues to use upstream mafft --auto and can differ.", type: "select", value: "mafft", options: [["mafft", "MAFFT E-INS-i 7.520 (Biowasm, default)"], ["muscle", "MUSCLE 5.1.0"]], token: "align_software" },
      { flag: "-t", label: "Threads per locus", description: "Browser fixed to one thread.", type: "select", value: "1", options: [["1", "1 (browser fixed)"]], token: "thread" },
      { flag: "-j", label: "Concurrent loci", description: "Browser edition processes loci sequentially.", type: "select", value: "1", options: [["1", "1 (browser sequential)"]], token: "jobs" },
      { flag: "-ts", label: "Trimming software", description: "Browser fixed to PhyloPrep's bundled trimAlnSeq.py.", type: "select", value: "trimAlnSeq", options: [["trimAlnSeq", "trimAlnSeq.py (browser fixed)"]], token: "trim_software" },
      { flag: "-n", label: "Skip trimming", description: "Keep raw alignments only.", type: "checkbox", value: false, token: "notrim" },
      proteinStopParameter(),
      { flag: "-G", label: "Maximum gap proportion", description: "Maximum gap ratio per site during trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "gap_ratio" },
      { flag: "-N", label: "Maximum N proportion", description: "Maximum N ratio per site during nucleotide/codon trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "n_ratio" },
      { flag: "-X", label: "Maximum X proportion", description: "Maximum X ratio per site during protein trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "x_ratio" },
    ],
    command: "python MSAP_batch.py -i {inputs} -o msap-results -st {seqtype} -s {align_software} -t 1 -j 1 -g {genetic_code} -ts trimAlnSeq -G {gap_ratio} -N {n_ratio} -X {x_ratio} --protein-stop-symbol {protein_stop_symbol_cli}{notrim_option}",
  },
  {
    script: "extract_orthofinder_orthogroups.py", title: "Extract OrthoFinder orthogroups", category: "Orthology", browser: true, browserWorkflow: "orthofinder-extract",
    summary: "Extracts selected single-copy or multi-copy orthogroups from Orthogroups.tsv using one FASTA file per taxon.",
    inputs: ["OrthoFinder Orthogroups.tsv", "Sequence FASTA files, one per taxon"],
    inputFields: [
      { key: "orthogroups", label: "Orthogroups.tsv", flag: "-og, --orthogroups", accept: ".tsv,.txt" },
      { key: "sequences", label: "Per-taxon FASTA files", flag: "-i, -f, --sequence-dir", multiple: true, accept: ".fa,.faa,.fasta,.fas,.pep" },
    ],
    outputs: ["One FASTA per retained orthogroup", "map.tsv", "orthogroups.pathlist"],
    outputLabel: "Output directory name", outputHelp: "Equivalent to -o/--outdir. Browser-selected FASTA files are staged as a temporary sequence directory.", defaultOutput: "one-to-one_orthologs",
    parameters: [
      { flag: "-m, --copy-mode", label: "Copy mode", description: "single-copy requires exactly one gene in each counted taxon; all-copies retains every listed copy.", type: "select", value: "single-copy", options: [["single-copy", "Single-copy"], ["all-copies", "All copies"]], token: "copy_mode" },
      { flag: "--min-taxa", label: "Minimum occupied taxa", description: "Leave blank to require all selected taxa.", type: "number", value: "", token: "min_taxa" },
      { flag: "-k, --keep-taxa", label: "Keep taxa", description: "Optional taxon names separated by commas or whitespace.", type: "text", value: "", token: "keep_taxa_text" },
      { flag: "-d, --drop-taxa", label: "Drop taxa", description: "Optional taxon names separated by commas or whitespace.", type: "text", value: "", token: "drop_taxa_text" },
      { flag: "-X", label: "Write unprefixed gene IDs", description: "Match OrthoFinder -X-style exported IDs.", type: "checkbox", value: false, token: "no_species_prefix" },
      { flag: "--orthofinder-id-style", label: "Orthogroup ID style", description: "Auto-detect normal taxon-prefixed or OrthoFinder -X unprefixed table IDs.", type: "select", value: "auto", options: [["auto", "Auto-detect"], ["prefixed", "Prefixed"], ["unprefixed", "Unprefixed"]], token: "orthofinder_id_style" },
      { flag: "-mo, --map-output", label: "Map file name", description: "Written inside the selected output directory.", type: "text", value: "map.tsv", token: "map_output_name" },
      { flag: "--pathlist-output", label: "Pathlist file name", description: "Written inside the selected output directory.", type: "text", value: "orthogroups.pathlist", token: "pathlist_output_name" },
    ],
    command: "python extract_orthofinder_orthogroups.py -og {orthogroups} -i SEQUENCE_DIR -o {output}{orthofinder_keep_taxa_option}{orthofinder_drop_taxa_option} -m {copy_mode}{min_taxa_option}{no_species_prefix_option} --orthofinder-id-style {orthofinder_id_style} -mo {map_output_name} --pathlist-output {pathlist_output_name}",
  },
  {
    script: "oma_to_orthogroups.py", title: "Convert OMA orthogroups", category: "Orthology", browser: true,
    summary: "Converts OMA Standalone OrthologousGroups.txt into the PhyloPrep orthogroup-by-taxon TSV standard.",
    inputs: ["OMA OrthologousGroups.txt"],
    inputFields: [{ key: "input", label: "OMA OrthologousGroups.txt", flag: "-i, --input", accept: ".txt,.tsv" }],
    outputs: ["Standard Orthogroup TSV"],
    outputLabel: "Output TSV", outputHelp: "Standard orthogroup-by-taxon table.", defaultOutput: "oma.orthogroups.tsv",
    parameters: [
      { flag: "--taxon-separator", label: "Taxon separator", description: "Genome/gene separator used by OMA members.", type: "text", value: "|", token: "taxon_separator" },
      { flag: "--strip-taxon-prefix", label: "Strip taxon prefix", description: "Write gene IDs without the leading genome and separator.", type: "checkbox", value: false, token: "strip_taxon_prefix" },
    ],
    command: "python oma_to_orthogroups.py -i {input} -o {output} --taxon-separator {taxon_separator_cli}{strip_taxon_prefix_option}",
  },
  {
    script: "orthomcl_to_orthogroups.py", title: "Convert OrthoMCL orthogroups", category: "Orthology", browser: true,
    summary: "Converts OrthoMCL groups.txt into the PhyloPrep orthogroup-by-taxon TSV standard.",
    inputs: ["OrthoMCL groups.txt"],
    inputFields: [{ key: "input", label: "OrthoMCL groups.txt", flag: "-i, --input", accept: ".txt,.tsv" }],
    outputs: ["Standard Orthogroup TSV"],
    outputLabel: "Output TSV", outputHelp: "Standard orthogroup-by-taxon table.", defaultOutput: "orthomcl.orthogroups.tsv",
    parameters: [
      { flag: "--taxon-separator", label: "Taxon separator", description: "Taxon/gene separator used by OrthoMCL members.", type: "text", value: "|", token: "taxon_separator" },
      { flag: "--strip-taxon-prefix", label: "Strip taxon prefix", description: "Write gene IDs without the leading taxon and separator.", type: "checkbox", value: false, token: "strip_taxon_prefix" },
    ],
    command: "python orthomcl_to_orthogroups.py -i {input} -o {output} --taxon-separator {taxon_separator_cli}{strip_taxon_prefix_option}",
  },
  {
    script: "proteinortho_to_orthogroups.py", title: "Convert Proteinortho groups", category: "Orthology", browser: true,
    summary: "Converts Proteinortho .proteinortho/.poff output into the PhyloPrep orthogroup-by-taxon TSV standard.",
    inputs: ["Proteinortho result table"],
    inputFields: [{ key: "input", label: "Proteinortho table", flag: "-i, --input", accept: ".proteinortho,.poff,.tsv,.txt" }],
    outputs: ["Standard Orthogroup TSV"],
    outputLabel: "Output TSV", outputHelp: "Standard orthogroup-by-taxon table.", defaultOutput: "proteinortho.orthogroups.tsv",
    parameters: [],
    command: "python proteinortho_to_orthogroups.py -i {input} -o {output}",
  },
  {
    script: "sonicparanoid_to_orthogroups.py", title: "Convert SonicParanoid groups", category: "Orthology", browser: true,
    summary: "Converts SonicParanoid flat.ortholog_groups.tsv into the PhyloPrep orthogroup-by-taxon TSV standard.",
    inputs: ["SonicParanoid flat.ortholog_groups.tsv"],
    inputFields: [{ key: "input", label: "SonicParanoid group table", flag: "-i, --input", accept: ".tsv,.txt" }],
    outputs: ["Standard Orthogroup TSV"],
    outputLabel: "Output TSV", outputHelp: "Standard orthogroup-by-taxon table.", defaultOutput: "sonicparanoid.orthogroups.tsv",
    parameters: [],
    command: "python sonicparanoid_to_orthogroups.py -i {input} -o {output}",
  },
  {
    script: "split_codon_seqence_alignment.py", title: "Split codon positions", category: "Alignment", browser: true,
    summary: "Splits a codon alignment into first-, second-, and third-position FASTA matrices.",
    inputs: ["Codon alignment FASTA"],
    inputFields: [{ key: "input", label: "Codon alignment FASTA", flag: "-i, --input", accept: ".fa,.fas,.fasta,.aln" }],
    outputs: ["codon1st.*.fasta", "codon2nd.*.fasta", "codon3rd.*.fasta"],
    outputLabel: "Output naming", outputHelp: "The original script derives all three output names from the input filename.", defaultOutput: "", hideOutputName: true,
    parameters: [],
    command: "python split_codon_seqence_alignment.py -i {input}",
  },
  {
    script: "AA2Codon.py", title: "Protein-alignment back-translation", category: "Back-translation", browser: true,
    summary: "Maps protein-alignment positions back to matching CDS codons. Sequence IDs must match exactly.",
    inputs: ["CDS FASTA (-c, --cds)", "Protein alignment FASTA (-p, --protaln)"],
    inputFields: [{ key: "cds", label: "CDS FASTA", flag: "-c, --cds" }, { key: "protaln", label: "Protein alignment FASTA", flag: "-p, --protaln" }],
    outputs: ["Codon alignment FASTA (-o, --out; stdout when omitted)"],
    outputLabel: "Output file name", outputHelp: "Codon alignment FASTA written by -o.", defaultOutput: "codon.aln.fasta",
    parameters: [geneticCodeParameter()],
    command: "python AA2Codon.py -c {cds} -p {protaln} -o {output} -g {genetic_code}",
  },
  {
    script: "alignment_qc.py", title: "Alignment quality control", category: "QC", browser: true, inputListKey: "inputs",
    summary: "Checks FASTA alignments locally and can combine observed alignment features with PhyloPrep sequence-change provenance, including masked internal stops.",
    inputs: ["FASTA files and/or pathlist files (one FASTA path per line) (-i, --input)", "Optional PhyloPrep sequence_changes.tsv provenance", "Required report prefix (-p, --prefix)", "Sequence type (-st, --seqtype)"],
    inputFields: [
      { key: "inputs", label: "Alignment input files", flag: "-i, --input", multiple: true, pathList: true },
      { key: "sequence_changes", label: "Sequence changes report", flag: "--sequence-changes", accept: ".tsv,.txt", required: false },
    ],
    outputs: ["PREFIX.pass.tsv", "PREFIX.fail.tsv", "PREFIX.details.tsv"],
    outputLabel: "Report prefix", outputHelp: "The command-line script creates its QC reports from this prefix.", defaultOutput: "phyloprep-qc",
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Expected sequence type for the command-line quality-control workflow.", type: "select", value: "nucl", options: [["nucl", "Nucleotide"], ["prot", "Protein"], ["codon", "Codon"]], token: "seqtype" },
      geneticCodeParameter(),
      { flag: "--min-length", label: "Minimum length", description: "Minimum alignment or sequence length. Default: 1.", type: "number", value: "1" },
      { flag: "--max-n-ratio", label: "Maximum N ratio", description: "Maximum N-base ratio for nucleotide matrices or N-containing codon ratio in codon matrices. Range: 0–1.", type: "number", value: "1", step: "0.01", token: "max_n_ratio" },
      { flag: "--max-x-ratio", label: "Maximum X ratio", description: "Maximum X ratio for protein matrices. Range: 0–1.", type: "number", value: "1", step: "0.01", token: "max_x_ratio" },
      { flag: "--max-gap-ratio", label: "Maximum gap ratio", description: "Proportion allowed before the alignment is rejected. Default: 1.", type: "number", value: "1", step: "0.01" },
      { flag: "--min-nucl-pis", label: "Minimum nucleotide PIS", description: "Minimum nucleotide-column PIS per locus in nucleotide/codon modes. Default: 0.", type: "number", value: "0", token: "min_nucl_pis" },
      { flag: "--min-protein-pis", label: "Minimum protein PIS", description: "Minimum amino-acid PIS per locus in protein/codon modes. Default: 0.", type: "number", value: "0", token: "min_protein_pis" },
      { flag: "--allow-terminal-stop", label: "Allow observed terminal stop", description: "Allow a terminal stop directly present in a codon alignment or a terminal * in a protein alignment. Stops recorded as removed in sequence_changes.tsv remain provenance only.", type: "checkbox", value: false, token: "allow_terminal_stop" },
      { flag: "--allow-internal-stop", label: "Allow internal stop", description: "Allow an internal stop observed in a codon/protein alignment or recorded as masked in sequence_changes.tsv. Protein direct-stop detection uses *; masking is preserved.", type: "checkbox", value: false, token: "allow_internal_stop" },
    ],
    command: "python alignment_qc.py -i {inputs} -p {prefix} -st {seqtype} -g {genetic_code}{sequence_changes_option} --min-length {min_length} --max-n-ratio {max_n_ratio} --max-x-ratio {max_x_ratio} --max-gap-ratio {max_gap_ratio}{min_nucl_pis_option}{min_protein_pis_option}{allow_terminal_stop}{allow_internal_stop}",
  },
  {
    script: "calulate_4dtv_and_correction.py", title: "4DTV calculation and correction", category: "Divergence", browser: true,
    summary: "Calculates fourfold-degenerate transversion values and applies the script's correction method.",
    inputs: ["One or more codon-alignment FASTA files (-i, --input)"], outputs: ["4DTV result table (-o, --out)"],
    inputFields: [{ key: "inputs", label: "Alignment input files", flag: "-i, --input", multiple: true, pathList: true }],
    outputLabel: "Output file name", outputHelp: "Tabular 4DTV output written by -o.", defaultOutput: "4dtv.tsv",
    parameters: [geneticCodeParameter()],
    command: "python calulate_4dtv_and_correction.py -i {inputs} -o {output} -g {genetic_code}",
  },
  {
    script: "extract_4-fold_degenerated_sites.py", title: "Extract four-fold degenerate sites", category: "Divergence", browser: true,
    summary: "Extracts four-fold-degenerate positions from codon-aligned FASTA records.",
    inputs: ["Codon-alignment FASTA (-i, --input)"], outputs: ["Four-fold-site FASTA (-o, --out)"],
    inputFields: [{ key: "input", label: "Alignment input files", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "FASTA file containing extracted sites.", defaultOutput: "fourfold-sites.fasta",
    parameters: [geneticCodeParameter()],
    command: "python extract_4-fold_degenerated_sites.py -i {input} -o {output} -g {genetic_code}",
  },
  {
    script: "get_supergenes.py", title: "Build supergenes", category: "Matrix", browser: true,
    summary: "Combines compatible single-copy alignments into concatenated supermatrices.",
    inputs: ["One or more alignment FASTA files (-i, --input)"], outputs: ["PREFIX.supergene.fasta", "Partition information and occupancy reports"],
    inputFields: [{ key: "inputs", label: "Alignment input files", flag: "-i, --input", multiple: true, pathList: true }],
    outputLabel: "Matrix prefix", outputHelp: "The command-line script derives supermatrix files from this prefix.", defaultOutput: "supergene",
    parameters: [{ flag: "--missing-taxa", label: "Missing taxa policy", description: "How records absent from a gene alignment are handled.", type: "select", value: "skip-gene", options: [["skip-gene", "skip-gene"], ["pad-gaps", "pad-gaps"]], token: "missing_taxa" }],
    command: "python get_supergenes.py -i {inputs} -p {prefix} --missing-taxa {missing_taxa}",
  },
  {
    script: "filter_sequences.py", title: "Filter sequences", category: "Filtering", browser: true,
    summary: "Filters aligned sequence data using length, ambiguity, gap, and taxa rules.",
    inputs: ["Alignment FASTA (-i, --input)", "Required sequence type (-st, --seqtype)"], outputs: ["Filtered alignment FASTA (-o, --out)"],
    inputFields: [{ key: "input", label: "Alignment FASTA", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "Filtered FASTA written by -o.", defaultOutput: "filtered.fasta",
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Required sequence type for filtering.", type: "select", value: "nucl", options: [["nucl", "Nucleotide"], ["prot", "Protein"], ["codon", "Codon"]], token: "seqtype" },
      { flag: "--min-length", label: "Minimum length", description: "Minimum retained sequence length.", type: "number", value: "0", token: "min_length" },
      { flag: "--max-gap-ratio", label: "Maximum gap ratio", description: "Maximum allowed gap fraction.", type: "number", value: "1", step: "0.01" },
      { flag: "--max-n-ratio", label: "Maximum N ratio", description: "Maximum allowed N fraction in nucleotide or codon sequences.", type: "number", value: "1", step: "0.01", token: "max_n_ratio" },
      { flag: "--max-x-ratio", label: "Maximum X ratio", description: "Maximum allowed X fraction in protein sequences.", type: "number", value: "1", step: "0.01", token: "max_x_ratio" },
      { flag: "--keep-taxa", label: "Keep taxa", description: "Taxon IDs are parsed automatically from the selected input FASTA.", type: "taxon-select", value: [], token: "keep_taxa" },
      { flag: "--remove-taxa", label: "Remove taxa", description: "Taxon IDs are parsed automatically from the selected input FASTA.", type: "taxon-select", value: [], token: "remove_taxa" },
    ],
    command: "python filter_sequences.py -i {input} -o {output} -st {seqtype} --min-length {min_length} --max-gap-ratio {max_gap_ratio} --max-n-ratio {max_n_ratio} --max-x-ratio {max_x_ratio}{keep_taxa_option}{remove_taxa_option}",
  },
  {
    script: "rename_taxa.py", title: "Rename taxa", category: "Identifiers", browser: true,
    summary: "Renames FASTA record IDs using a mapping file, with safeguards for duplicate identifiers.",
    inputs: ["FASTA input (-i, --input)", "Taxon mapping file (-m, --map)"], outputs: ["Renamed FASTA (-o, --out)", "Optional map audit (--map-output)"],
    inputFields: [{ key: "input", label: "Input FASTA", flag: "-i, --input" }, { key: "mapping", label: "Taxon mapping table", flag: "-m, --map", accept: ".tsv,.txt,.csv" }],
    outputLabel: "Output file name", outputHelp: "Renamed FASTA written by -o.", defaultOutput: "renamed.fasta",
    parameters: [
      { flag: "--map-output", label: "Write map audit", description: "Save the resolved mapping table beside the output.", type: "checkbox", value: false, token: "map_output" },
      { flag: "--allow-duplicate-ids", label: "Allow duplicate renamed IDs", description: "Write duplicated IDs for downstream validation; normally they are rejected.", type: "checkbox", value: false, token: "allow_duplicate_ids" },
      { flag: "--keep-special", label: "Keep special IDs", description: "Preserve special record identifiers when supported by the script.", type: "checkbox", value: false, token: "keep_special" },
    ],
    command: "python rename_taxa.py -i {input} -o {output}{mapping_option}{map_output_option}{allow_duplicate_ids_option}{keep_special_option}",
  },
  {
    script: "trimAlnSeq.py", title: "Trim alignment sequences", category: "Trimming", browser: true,
    summary: "Removes low-quality columns or sequences from an alignment under explicit thresholds.",
    inputs: ["Alignment FASTA (-i, --input)", "Required sequence type (-st, --seqtype)"], outputs: ["Trimmed alignment FASTA (-o, --out; optional)"],
    inputFields: [{ key: "input", label: "Alignment FASTA", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "Trimmed FASTA written by -o when provided.", defaultOutput: "trimmed.aln",
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Alignment sequence type.", type: "select", value: "codon", options: [["codon", "Codon"], ["prot", "Protein"], ["nucl", "Nucleotide"]], token: "seqtype" },
      { flag: "-G", label: "Maximum gap proportion", description: "Column gap threshold for trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "gap_ratio" },
      { flag: "-N", label: "Maximum N proportion", description: "Column N threshold for nucleotide or codon alignments. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "n_ratio" },
      { flag: "-X", label: "Maximum X proportion", description: "Column X threshold for protein alignments. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "x_ratio" },
      { flag: "-ra", label: "Remove invariant columns", description: "Remove alignment columns with one resolved state.", type: "checkbox", value: false, token: "removeall" },
    ],
    command: "python trimAlnSeq.py -i {input} -o {output} -st {seqtype} -G {gap_ratio} -N {n_ratio} -X {x_ratio}{removeall}",
  },
  {
    script: "fasta2axt.py", title: "Convert FASTA to AXT", category: "Convert", browser: true,
    summary: "Converts exactly two equal-length aligned FASTA records into AXT format locally.",
    inputs: ["Two-sequence aligned FASTA (-i, --input)"], outputs: ["AXT file (-o, --out)"],
    inputFields: [{ key: "input", label: "Two-sequence aligned FASTA", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "AXT output generated locally in the browser.", defaultOutput: "alignment.axt",
    parameters: [{ flag: "-st", label: "Sequence type", description: "Sequence type used by the converter.", type: "select", value: "DNA", options: [["DNA", "DNA"], ["protein", "Protein"]], token: "seqtype" }], command: "python fasta2axt.py -i {input} -o {output} -st {seqtype}",
  },
  {
    script: "fasta2nex.py", title: "Convert FASTA to NEXUS", category: "Convert", browser: true,
    summary: "Converts an aligned FASTA matrix to NEXUS format locally.",
    inputs: ["Aligned FASTA (-i, --input)"], outputs: ["NEXUS file (-o, --out)"],
    inputFields: [{ key: "input", label: "Aligned FASTA", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "NEXUS output generated locally in the browser.", defaultOutput: "alignment.nex",
    parameters: [{ flag: "-st", label: "Sequence type", description: "Alignment type used by the converter.", type: "select", value: "DNA", options: [["DNA", "DNA"], ["RNA", "RNA"], ["protein", "Protein"]], token: "seqtype" }], command: "python fasta2nex.py -i {input} -o {output} -st {seqtype}",
  },
  {
    script: "fasta2phy.py", title: "Convert FASTA to PHYLIP", category: "Convert", browser: true,
    summary: "Converts an aligned FASTA matrix to PHYLIP format locally.",
    inputs: ["Aligned FASTA (-i, --input)"], outputs: ["PHYLIP file (-o, --out)"],
    inputFields: [{ key: "input", label: "Aligned FASTA", flag: "-i, --input" }],
    outputLabel: "Output file name", outputHelp: "PHYLIP output generated locally in the browser.", defaultOutput: "alignment.phy",
    parameters: [{ flag: "-st", label: "Sequence type", description: "Sequence type used by the converter.", type: "select", value: "DNA", options: [["DNA", "DNA"], ["protein", "Protein"]], token: "seqtype" }], command: "python fasta2phy.py -i {input} -o {output} -st {seqtype}",
  },
];

let activeModuleIndex = 0;
let outputDirectoryHandle;
let generatedFiles = [];
let selectedGeneratedFileIndex = 0;
let resultEntries = [];
let resultPageNumber = 0;
let pyodideWorker;
let activeWorkerRun;
let nextRunId = 0;
let selectedInputFiles = new Map();
let alignmentWorker;
let activeAlignmentRun;
let alignmentRunCounter = 0;
let catalogQuery = "";
let pendingTransfer;

function activeModule() { return modules[activeModuleIndex]; }

function implementationProfile(module) {
  if (module.browserWorkflow === "phyloprep") {
    return {
      key: "compatible",
      label: "Browser PhyloPrep pipeline",
      note: "The browser orchestrates the repository's original preparation modules locally: MSAP preprocessing/back-translation/trimming in Pyodide, MUSCLE 5.1.0 or optional MAFFT E-INS-i 7.520 in Biowasm, then provenance-aware paired codon/protein QC, renaming, concatenation, codon splitting, four-fold extraction and format-conversion code. Linux process-management features are replaced by a browser-safe sequential adapter.",
    };
  }
  if (module.browserWorkflow === "batch") {
    return {
      key: "compatible",
      label: "Browser batch adapter",
      note: "Processes loci sequentially in the browser with the same MSAP sequence workflow and Biowasm MUSCLE 5.1.0 or MAFFT E-INS-i 7.520, then generates batch-style path lists and a merged sequence-change report. Linux process groups, file locks and concurrent jobs are intentionally not used in the browser.",
    };
  }
  if (module.fullWorkflow) {
    return {
      key: "compatible",
      label: "Original MSAP.py + Biowasm aligner",
      note: "Runs the repository's untouched MSAP.py in Pyodide. Its MAFFT or MUSCLE subprocess call is fulfilled locally by Biowasm/WebAssembly; bundled AA2Codon.py and trimAlnSeq.py continue inside the same browser Python runtime. No sequence file is uploaded to PhyloPrep.",
    };
  }
  return {
    key: "original",
    label: "Original script in Pyodide",
    note: `Runs the repository's original ${module.script} source locally in the browser, with the same CLI implementation and bundled dependencies.`,
  };
}

function normalizeCatalogText(value) {
  return String(value || "")
    .toLowerCase()
    .replace(/\.py\b/g, "")
    .replace(/[^a-z0-9]+/g, " ")
    .trim();
}

function catalogSearchText(module) {
  return normalizeCatalogText([
    module.script,
    module.title,
    module.category,
    module.summary,
    ...(module.inputs || []),
    ...(module.outputs || []),
    ...(module.parameters || []).flatMap((parameter) => [parameter.label, parameter.flag, parameter.description]),
  ].join(" "));
}

function renderCatalog() {
  elements.scriptList.replaceChildren();
  const query = normalizeCatalogText(catalogQuery);
  const tokens = query.split(/\s+/).filter(Boolean);
  const matched = modules
    .map((module, index) => ({ module, index }))
    .filter(({ module }) => {
      if (!tokens.length) return true;
      const haystack = catalogSearchText(module);
      const compact = haystack.replace(/\s+/g, "");
      return tokens.every((token) => haystack.includes(token) || compact.includes(token));
    });
  if (!matched.length) {
    const empty = document.createElement("p");
    empty.className = "catalog-empty";
    empty.textContent = "No matching tools.";
    elements.scriptList.append(empty);
    return;
  }
  const groups = new Map();
  matched.forEach((entry) => {
    const category = entry.module.category || "Other";
    if (!groups.has(category)) groups.set(category, []);
    groups.get(category).push(entry);
  });
  groups.forEach((entries, category) => {
    const group = document.createElement("details");
    group.className = "catalog-group";
    group.open = true;
    const summary = document.createElement("summary");
    summary.textContent = `${category} · ${entries.length}`;
    group.append(summary);
    entries.forEach(({ module, index }) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = `script-button${index === activeModuleIndex ? " active" : ""}`;
      button.setAttribute("role", "tab");
      button.setAttribute("aria-selected", String(index === activeModuleIndex));
      const name = document.createElement("span");
      name.className = "script-name";
      name.textContent = module.script;
      const title = document.createElement("span");
      title.className = "script-category";
      title.textContent = module.title;
      button.append(name, title);
      button.addEventListener("click", () => selectModule(index));
      group.append(button);
    });
    elements.scriptList.append(group);
  });
}

function parameterToken(parameter) {
  return parameter.token || parameter.flag.replace(/[^a-z0-9]+/gi, "_").replace(/^_|_$/g, "").toLowerCase();
}

function renderParameterControl(parameter) {
  const card = document.createElement("section");
  card.className = "parameter-card";
  card.dataset.parameterToken = parameterToken(parameter);
  const id = `parameter-${parameterToken(parameter)}`;
  if (parameter.type === "checkbox") {
    const label = document.createElement("label");
    label.className = "checkbox-row";
    const control = document.createElement("input");
    control.id = id;
    control.type = "checkbox";
    control.checked = Boolean(parameter.value);
    control.dataset.parameter = parameterToken(parameter);
    const text = document.createElement("span");
    text.textContent = parameter.label;
    label.append(control, text);
    card.append(label);
  } else {
    const label = document.createElement("label");
    label.htmlFor = id;
    const labelText = document.createElement("span");
    labelText.textContent = parameter.label;
    const flag = document.createElement("code");
    flag.textContent = parameter.flag;
    label.append(labelText, flag);
    let control;
    if (parameter.type === "taxon-select") {
      control = document.createElement("select");
      control.multiple = true;
      control.size = 6;
      control.disabled = true;
      const placeholder = document.createElement("option");
      placeholder.textContent = "Select an input FASTA file first";
      placeholder.disabled = true;
      control.append(placeholder);
    } else if (parameter.type === "select") {
      control = document.createElement("select");
      parameter.options.forEach(([value, text]) => {
        const option = document.createElement("option");
        option.value = value;
        option.textContent = text;
        option.selected = value === parameter.value;
        control.append(option);
      });
    } else {
      control = document.createElement("input");
      control.type = parameter.type || "text";
      control.value = parameter.value ?? "";
      if (parameter.step) control.step = parameter.step;
    }
    control.id = id;
    control.dataset.parameter = parameterToken(parameter);
    const description = document.createElement("p");
    description.textContent = parameter.description;
    if (parameter.reference) {
      const reference = document.createElement("a");
      reference.href = parameter.reference;
      reference.target = "_blank";
      reference.rel = "noreferrer";
      reference.textContent = " Open NCBI reference ↗";
      description.append(reference);
    }
    card.append(label, description, control);
  }
  return card;
}

function renderParameters(module) {
  elements.parameterList.replaceChildren();
  if (!module.parameters.length) {
    const message = document.createElement("p");
    message.className = "no-parameters";
    message.textContent = "This converter has no additional script parameters.";
    elements.parameterList.append(message);
  } else {
    const priority = (parameter) => {
      const token = parameterToken(parameter);
      if (token === "seqtype") return 0;
      if (token === "genetic_code") return 1;
      return 2;
    };
    const sortedParameters = module.parameters
      .map((parameter, index) => ({ parameter, index }))
      .sort((left, right) => priority(left.parameter) - priority(right.parameter) || left.index - right.index)
      .map(({ parameter }) => parameter);
    const coreTokens = new Set(["seqtype", "genetic_code", "missing_taxa", "keep_taxa", "remove_taxa"]);
    const coreParameters = sortedParameters.filter((parameter) => coreTokens.has(parameterToken(parameter)));
    const advancedParameters = sortedParameters.filter((parameter) => !coreTokens.has(parameterToken(parameter)));

    if (coreParameters.length) {
      const group = document.createElement("section");
      group.className = "parameter-group";
      const title = document.createElement("p");
      title.className = "parameter-group-title";
      title.textContent = "Core options";
      group.append(title, ...coreParameters.map(renderParameterControl));
      elements.parameterList.append(group);
    }
    if (advancedParameters.length) {
      const group = document.createElement("details");
      group.className = "advanced-parameters";
      const summary = document.createElement("summary");
      summary.textContent = `Advanced options (${advancedParameters.length})`;
      const body = document.createElement("div");
      body.className = "advanced-parameters-body";
      body.append(...advancedParameters.map(renderParameterControl));
      group.append(summary, body);
      elements.parameterList.append(group);
    }
    elements.parameterList.querySelectorAll("input, select").forEach((control) => {
      const refresh = () => { updateContextualParameters(); updateCommand(); void updatePreflight(); };
      control.addEventListener("input", refresh);
      control.addEventListener("change", refresh);
    });
    updateContextualParameters();
  }
}

function updateContextualParameters() {
  const seqtype = elements.parameterList.querySelector('[data-parameter="seqtype"]')?.value;
  const notrim = elements.parameterList.querySelector('[data-parameter="notrim"]')?.checked;
  const script = activeModule().script;
  const hide = (token, hidden) => {
    const card = elements.parameterList.querySelector(`[data-parameter-token="${token}"]`);
    if (card) card.dataset.contextHidden = String(Boolean(hidden));
  };
  hide("genetic_code", seqtype === "prot" || seqtype === "nucl");
  hide("protein_stop_symbol", seqtype !== "prot");
  hide("n_ratio", seqtype === "prot" || Boolean(notrim));
  hide("x_ratio", seqtype === "nucl" || Boolean(notrim));
  hide("gap_ratio", Boolean(notrim));
  hide("max_n_ratio", seqtype === "prot");
  hide("max_x_ratio", seqtype === "nucl" || (script === "alignment_qc.py" && seqtype === "codon"));
  hide("allow_terminal_stop", script !== "alignment_qc.py" || !["codon", "prot"].includes(seqtype));
  hide("allow_internal_stop", seqtype === "nucl");
  const pairedPis = script === "phyloprep.py" && seqtype === "codon";
  hide("paired_pis_policy", !pairedPis);
  hide("min_nucl_pis", seqtype === "prot");
  hide("min_protein_pis", !(seqtype === "prot" || pairedPis));
}

function inputControls() {
  return [...elements.inputFileFields.querySelectorAll(".module-input")];
}

function filesForInput(control) {
  return selectedInputFiles.get(control.dataset.inputKey) || [];
}

function selectedInputNames() {
  return inputControls().flatMap((control) => filesForInput(control).map((file) => file.name));
}

function browserInputFiles() {
  return inputControls().flatMap((control) => filesForInput(control));
}

function fileIdentity(file) {
  return `${file.webkitRelativePath || file.name}\u0000${file.size}\u0000${file.lastModified}`;
}

function rememberInputFiles(control) {
  const selected = [...control.files];
  const existing = filesForInput(control);
  const files = control.multiple
    ? [...new Map([...existing, ...selected].map(file => [fileIdentity(file), file])).values()]
    : selected.slice(0, 1);
  selectedInputFiles.set(control.dataset.inputKey, files);
  control.value = "";
  void renderInputFiles();
}

function removeInputFile(inputKey, index) {
  const remaining = [...(selectedInputFiles.get(inputKey) || [])];
  remaining.splice(index, 1);
  selectedInputFiles.set(inputKey, remaining);
  void renderInputFiles();
}

function clearAllInputFiles() {
  inputControls().forEach((control) => {
    selectedInputFiles.set(control.dataset.inputKey, []);
    control.value = "";
  });
  void renderInputFiles();
}

function virtualInputPath(name) {
  return `/work/input/${name}`;
}

function pathListPath(module = activeModule()) {
  return `/work/pathlists/${module.script.replace(/\.py$/i, "")}.pathlist`;
}

function renderInputFileFields(module) {
  elements.inputFileFields.replaceChildren();
  selectedInputFiles = new Map();
  inputFolderHandles.clear();
  inputFolderFiles.clear();
  module.inputFields.forEach((field) => {
    const section = document.createElement("section");
    section.className = "workspace-section input-file-field";
    const heading = document.createElement("div");
    heading.className = "field-heading";
    const label = document.createElement("label");
    const id = `module-input-${field.key}`;
    label.htmlFor = id;
    label.textContent = field.label;
    const flag = document.createElement("span");
    flag.className = "required-label";
    flag.textContent = field.flag;
    if (field.required === false) flag.textContent += " · optional";
    heading.append(label, flag);
    const picker = document.createElement("label");
    picker.className = "folder-picker input-file-picker";
    picker.htmlFor = id;
    const control = document.createElement("input");
    control.className = "module-input input-file-control";
    control.id = id;
    control.dataset.inputKey = field.key;
    control.dataset.inputLabel = field.label;
    control.dataset.inputRequired = String(field.required !== false);
    control.type = "file";
    control.multiple = Boolean(field.multiple);
    const acceptedFiles = field.accept || ".fa,.fas,.fasta,.aln,.phy,.phylip,.nex,.nexus,.txt";
    control.accept = field.pathList ? `${acceptedFiles},.pathlist,.list` : acceptedFiles;
    control.addEventListener("change", () => rememberInputFiles(control));
    const display = document.createElement("span");
    display.className = "input-file-display";
    display.textContent = field.multiple ? "Select local files" : "Select a local file";
    picker.append(control, display);
    section.append(heading, picker);
    if (field.pathList) {
      const folderButton = document.createElement("button");
      folderButton.type = "button";
      folderButton.className = "folder-picker";
      folderButton.dataset.pathlistFolderKey = field.key;
      folderButton.textContent = "Select pathlist input folder";
      folderButton.addEventListener("click", async () => {
        try {
          if ("showDirectoryPicker" in window) {
            const handle = await window.showDirectoryPicker({ mode: "read" });
            inputFolderHandles.set(field.key, handle);
            folderButton.textContent = `Input folder: ${handle.name}`;
            if (!filesForInput(control).length) {
              const found = [];
              for await (const [name, entry] of handle.entries()) {
                if (entry.kind === "file" && /\.(pathlist|list|fa|faa|fas|fasta|aln|pep)$/i.test(name)) found.push(await entry.getFile());
              }
              const lists = found.filter(file => /\.(pathlist|list)$/i.test(file.name));
              selectedInputFiles.set(field.key, lists.length ? lists : found);
            }
            await renderInputFiles();
          } else {
            directoryControl.click();
          }
        } catch (error) {
          if (error.name !== "AbortError") setLog(`Cannot access input folder: ${error.message}`, "error");
        }
      });
      const directoryControl = document.createElement("input");
      directoryControl.type = "file";
      directoryControl.multiple = true;
      directoryControl.setAttribute("webkitdirectory", "");
      directoryControl.hidden = true;
      directoryControl.addEventListener("change", () => {
        const files = [...directoryControl.files];
        inputFolderFiles.set(field.key, files);
        if (!filesForInput(control).length) {
          const lists = files.filter(file => /\.(pathlist|list)$/i.test(file.name));
          selectedInputFiles.set(field.key, lists.length ? lists : files.filter(file => /\.(fa|faa|fas|fasta|aln|pep)$/i.test(file.name)));
        }
        folderButton.textContent = "Input folder selected";
        void renderInputFiles();
      });
      const notice = document.createElement("p");
      notice.dataset.pathlistNoticeKey = field.key;
      notice.setAttribute("role", "status");
      notice.hidden = true;
      section.append(folderButton, directoryControl, notice);
    }
    elements.inputFileFields.append(section);
  });
}

function cliQuote(value) {
  const text = String(value);
  return /^[A-Za-z0-9_./:=+,-]+$/u.test(text) ? text : JSON.stringify(text);
}

function rawParameterValues() {
  const module = activeModule();
  const parameterValues = {};
  module.parameters.forEach((parameter) => {
    const control = elements.parameterList.querySelector(`[data-parameter="${parameterToken(parameter)}"]`);
    parameterValues[parameterToken(parameter)] = control?.type === "checkbox"
      ? control.checked
      : control?.multiple
        ? [...control.selectedOptions].map((option) => option.value)
        : (control?.value || parameter.value || "");
  });
  if (parameterValues.seqtype === "prot") delete parameterValues.min_nucl_pis;
  if (parameterValues.seqtype === "nucl" || (module.script === "alignment_qc.py" && parameterValues.seqtype === "codon")) {
    delete parameterValues.min_protein_pis;
  }
  return parameterValues;
}

function valuesForCommand() {
  const module = activeModule();
  const allInputNames = selectedInputNames();
  const listControl = module.inputListKey
    ? inputControls().find((control) => control.dataset.inputKey === module.inputListKey)
    : undefined;
  const commandInputNames = listControl
    ? filesForInput(listControl).map((file) => file.name)
    : allInputNames;
  const parameterValues = rawParameterValues();
  const rawOutput = elements.outputName.value.trim() || module.defaultOutput;
  const outputPath = outputDirectoryHandle ? `${outputDirectoryHandle.name}/${rawOutput}` : rawOutput;
  const usePathList = module.inputFields.some((field) => field.pathList) && commandInputNames.length > 1;
  const inputPaths = commandInputNames.map((name) => cliQuote(virtualInputPath(name)));
  const mappingFile = selectedInputFiles.get("mapping")?.[0];
  const mappingName = mappingFile?.name || "taxa-map.tsv";
  const mapStem = rawOutput.replace(/\.[^.]+$/, "") || "renamed";
  const mapOutputPath = outputDirectoryHandle ? `${outputDirectoryHandle.name}/${mapStem}.map.tsv` : `${mapStem}.map.tsv`;
  return {
    ...parameterValues,
    input: inputPaths[0] || "/work/input/input.fasta",
    inputs: usePathList ? pathListPath(module) : (inputPaths.length ? inputPaths.join(" ") : "/work/input/input.fasta"),
    orthogroups: selectedInputFiles.get("orthogroups")?.[0] ? cliQuote(virtualInputPath(selectedInputFiles.get("orthogroups")[0].name)) : "/work/input/Orthogroups.tsv",
    cds: inputPaths[0] || "/work/input/cds.fasta",
    protaln: inputPaths[1] || "/work/input/protein.aln.fasta",
    output: cliQuote(outputPath),
    output_dir: module.cliOutputDir || "results",
    prefix: cliQuote(outputPath),
    seqtype: parameterValues.seqtype || "nucl",
    genetic_code: parameterValues.genetic_code || "1",
    protein_stop_symbol: parameterValues.protein_stop_symbol || "*",
    protein_stop_symbol_cli: cliQuote(parameterValues.protein_stop_symbol || "*"),
    taxon_separator_cli: cliQuote(parameterValues.taxon_separator || "|"),
    allow_terminal_stop: parameterValues.allow_terminal_stop ? " --allow-terminal-stop" : "",
    allow_internal_stop: parameterValues.allow_internal_stop ? " --allow-internal-stop" : "",
    removeall: parameterValues.removeall ? " -ra" : "",
    keep_taxa_option: parameterValues.keep_taxa?.length ? " --keep-taxa /work/path-lists/filter_sequences.keep-taxa.txt" : "",
    remove_taxa_option: parameterValues.remove_taxa?.length ? " --remove-taxa /work/path-lists/filter_sequences.remove-taxa.txt" : "",
    mapping_option: mappingFile ? ` -m ${cliQuote(virtualInputPath(mappingName))}` : "",
    sequence_changes_option: selectedInputFiles.get("sequence_changes")?.[0]
      ? ` --sequence-changes ${cliQuote(virtualInputPath(selectedInputFiles.get("sequence_changes")[0].name))}`
      : "",
    min_nucl_pis_option: parameterValues.seqtype !== "prot" && parameterValues.min_nucl_pis !== undefined && parameterValues.min_nucl_pis !== ""
      ? ` --min-nucl-pis ${parameterValues.min_nucl_pis}`
      : "",
    min_protein_pis_option: ["prot", "codon", "pseudogene"].includes(parameterValues.seqtype) && parameterValues.min_protein_pis !== undefined && parameterValues.min_protein_pis !== ""
      ? ` --min-protein-pis ${parameterValues.min_protein_pis}`
      : "",
    map_output_option: parameterValues.map_output ? ` --map-output ${cliQuote(mapOutputPath)}` : "",
    allow_duplicate_ids_option: parameterValues.allow_duplicate_ids ? " --allow-duplicate-ids" : "",
    keep_special_option: parameterValues.keep_special ? " --keep-special" : "",
    notrim_option: parameterValues.notrim ? " -n" : "",
    strip_taxon_prefix_option: parameterValues.strip_taxon_prefix ? " --strip-taxon-prefix" : "",
    no_species_prefix_option: parameterValues.no_species_prefix ? " -X" : "",
    min_taxa_option: parameterValues.min_taxa ? ` --min-taxa ${parameterValues.min_taxa}` : "",
    orthofinder_keep_taxa_option: parameterValues.keep_taxa_text
      ? ` -k ${parameterValues.keep_taxa_text.split(/[\s,]+/).filter(Boolean).map(cliQuote).join(" ")}`
      : "",
    orthofinder_drop_taxa_option: parameterValues.drop_taxa_text
      ? ` -d ${parameterValues.drop_taxa_text.split(/[\s,]+/).filter(Boolean).map(cliQuote).join(" ")}`
      : "",
  };
}

function updateCommand() {
  const values = valuesForCommand();
  let command = activeModule().command.replace(/\{([a-z_]+)\}/g, (_, key) => values[key] ?? `{${key}}`);
  if (["phyloprep.py", "MSAP_batch.py"].includes(activeModule().script) && values.align_software === "mafft") {
    command += "\n# Web note: MAFFT E-INS-i 7.520 (Biowasm) is selected here; CLI -s mafft runs upstream mafft --auto and may produce different alignments.";
  }
  elements.command.textContent = command;
}

function setRunStage(stage = 0, label = "Ready to run", failed = false) {
  elements.runStageLabel.textContent = label;
  [...elements.runStages.querySelectorAll(".run-stage")].forEach((item) => {
    const number = Number(item.dataset.stage);
    item.classList.remove("active", "complete", "error");
    if (stage === 4 || (stage > 0 && number < stage)) item.classList.add("complete");
    if (stage === number) item.classList.add(failed ? "error" : "active");
    if (stage === 4) item.classList.add("complete");
  });
}

function renderEmptyOutputFiles() {
  elements.resultControls.hidden = true;
  resultEntries = [];
  resultPageNumber = 0;
  elements.outputFilesList.replaceChildren();
  elements.downloadSelected.hidden = true;
  elements.downloadAll.hidden = true;
  const message = document.createElement("p");
  message.className = "empty-state";
  message.textContent = "Generated files will appear here after the workflow runs.";
  elements.outputFilesList.append(message);
}

function resetGeneratedOutput() {
  generatedFiles = [];
  renderEmptyOutputFiles();
  elements.resultCount.textContent = "No output yet";
  elements.preview.textContent = "Run the selected module to inspect generated output here.";
  elements.useNextStep.hidden = true;
  setRunStage();
}

function setLog(message, status = "") {
  elements.log.textContent = message;
  elements.log.className = `workspace-log${status ? ` ${status}` : ""}`;
}

function selectModule(index) {
  activeModuleIndex = index;
  const module = activeModule();
  elements.selectedName.textContent = module.title;
  elements.selectedSummary.textContent = module.summary;
  elements.workspaceScript.textContent = module.script;
  const profile = implementationProfile(module);
  elements.workspaceMode.textContent = profile.label;
  elements.workspaceMode.className = `mode-badge ${profile.key}`;
  elements.implementationNote.textContent = profile.note;
  elements.outputNameSection.hidden = Boolean(module.hideOutputName);
  elements.outputNameLabel.textContent = module.outputLabel;
  elements.outputName.value = module.defaultOutput;
  elements.outputNameHelp.textContent = module.outputHelp;
  elements.inputFilesHelp.textContent = module.inputFields.some((field) => field.pathList)
    ? "Select FASTA files or a pathlist file. Multiple selected files are merged into a browser pathlist."
    : "Select the files required by this module.";
  renderInputFileFields(module);
  elements.inputFilesList.textContent = "No files selected.";
  elements.clearInputFiles.hidden = true;
  if (pendingTransfer?.script === module.script) {
    selectedInputFiles = pendingTransfer.inputs;
    pendingTransfer = undefined;
    void renderInputFiles();
  }
  renderCatalog();
  renderParameters(module);
  resetGeneratedOutput();
  setLog("Select the required input file or files, then run the local workflow.");
  updateCommand();
}

async function refreshTaxonSelectors() {
  const selectors = [...elements.parameterList.querySelectorAll('[data-parameter="keep_taxa"], [data-parameter="remove_taxa"]')];
  if (!selectors.length) return;
  const file = browserInputFiles()[0];
  let taxa = [];
  if (file) {
    try {
      taxa = (await file.text()).split(/\r?\n/)
        .filter((line) => line.startsWith(">"))
        .map((line) => line.slice(1).trim().split(/\s+/, 1)[0])
        .filter(Boolean);
      taxa = [...new Set(taxa)];
    } catch {
      taxa = [];
    }
  }
  selectors.forEach((selector) => {
    const previouslySelected = new Set([...selector.selectedOptions].map((option) => option.value));
    selector.replaceChildren();
    if (!taxa.length) {
      const placeholder = document.createElement("option");
      placeholder.textContent = "Select an input FASTA file first";
      placeholder.disabled = true;
      selector.append(placeholder);
      selector.disabled = true;
      return;
    }
    taxa.forEach((taxon) => {
      const option = document.createElement("option");
      option.value = taxon;
      option.textContent = taxon;
      option.selected = previouslySelected.has(taxon);
      selector.append(option);
    });
    selector.disabled = false;
  });
}

async function updatePathlistFolderNotices() {
  const generation = ++pathlistNoticeGeneration;
  for (const control of inputControls()) {
    const key = control.dataset.inputKey;
    const section = control.closest(".input-file-field");
    const notice = section.querySelector("[data-pathlist-notice-key]");
    const button = section.querySelector("[data-pathlist-folder-key]");
    if (!notice || !button) continue;
    const selected = filesForInput(control);
    const lists = selected.filter(file => /\.(pathlist|list)$/i.test(file.name));
    notice.hidden = !lists.length;
    if (!lists.length) {
      if (!inputFolderHandles.has(key) && !inputFolderFiles.has(key)) button.textContent = "Select pathlist input folder";
      continue;
    }
    notice.textContent = "Reading pathlist references…";
    try {
      const texts = await PhyloPrepInputs.readBounded(lists);
      if (generation !== pathlistNoticeGeneration) return;
      const parsed = texts.map(PhyloPrepInputs.pathlistFolders);
      const folders = [...new Set(parsed.flatMap(item => item.folders))];
      const count = parsed.reduce((total, item) => total + item.count, 0);
      const folderSummary = folders.slice(0, 3).join("; ") + (folders.length > 3 ? `; and ${folders.length - 3} more` : "");
      const references = `Pathlist contains ${count} reference(s). Recorded folders: ${folderSummary || "none"}. `;
      if (inputFolderHandles.has(key) || inputFolderFiles.has(key)) {
        notice.textContent = references + "Input folder selected. Preflight will check referenced files and same-folder basenames.";
      } else {
        let needsFolder = false;
        try { await PhyloPrepInputs.resolveFiles(selected); }
        catch { needsFolder = true; }
        if (generation !== pathlistNoticeGeneration) return;
        if (needsFolder) {
          button.textContent = "Authorize pathlist input folder";
          notice.textContent = references + "Click Authorize pathlist input folder and select the folder containing these FASTAs. If the data were moved, select their current folder (usually beside the pathlist). Browser permission requires your confirmation; the listed paths do not grant access.";
        } else {
          button.textContent = "Select pathlist input folder";
          notice.textContent = references + "Referenced FASTAs are already selected; no additional folder access is needed.";
        }
      }
    } catch (error) {
      if (generation !== pathlistNoticeGeneration) return;
      notice.textContent = error.message;
    }
  }
}

async function renderInputFiles() {
  const controls = inputControls();
  const names = selectedInputNames();
  elements.inputFilesList.replaceChildren();
  controls.forEach((control) => {
    const selected = filesForInput(control);
    const display = control.closest(".input-file-picker")?.querySelector(".input-file-display");
    if (!display) return;
    display.textContent = selected.length
      ? selected.length > 20 ? `${selected.length} files selected` : selected.map((file) => file.name).join(", ")
      : (control.multiple ? "Select local files" : "Select a local file");
    display.title = selected.length > 20 ? `${selected.length} files selected` : selected.map((file) => file.name).join(", ");
  });
  elements.clearInputFiles.hidden = names.length === 0;
  if (!names.length) {
    elements.inputFilesList.textContent = "No files selected.";
  } else {
    controls.forEach((control) => {
      filesForInput(control).slice(0, 100).forEach((file, index) => {
        const chip = document.createElement("span");
        chip.className = "file-chip";
        chip.title = `${control.dataset.inputLabel}: ${file.name}`;
        const fileName = document.createElement("span");
        fileName.className = "file-chip-name";
        fileName.textContent = `${control.dataset.inputLabel}: ${file.name}`;
        const remove = document.createElement("button");
        remove.className = "file-chip-remove";
        remove.type = "button";
        remove.setAttribute("aria-label", `Remove ${file.name}`);
        remove.title = `Remove ${file.name}`;
        remove.textContent = "×";
        remove.addEventListener("click", () => removeInputFile(control.dataset.inputKey, index));
        chip.append(fileName, remove);
        elements.inputFilesList.append(chip);
      });
    });
  }
  if (names.length > 100) {
    const note = document.createElement("span");
    note.textContent = "Showing the first 100 files per input; all selected files will be processed.";
    elements.inputFilesList.append(note);
  }
  const module = activeModule();
  const pathLists = browserInputFiles().filter((file) => /\.(pathlist|list)$/i.test(file.name));
  const fastaFiles = browserInputFiles().filter((file) => !/\.(pathlist|list)$/i.test(file.name));
  if (module.inputFields.some((field) => field.pathList) && names.length > 1) {
    const browserNote = pathLists.length && !fastaFiles.length
      ? " Select the pathlist input folder or its referenced FASTA files."
      : "";
    elements.inputFilesHelp.textContent = `Multiple selections are merged into a browser pathlist for the CLI preview: ${pathListPath(module)}.${browserNote}`;
  } else if (module.inputFields.some((field) => field.pathList) && pathLists.length) {
    elements.inputFilesHelp.textContent = "A .pathlist file is selected. Select its input folder or referenced FASTA files. Missing paths fall back to basenames beside the list.";
  } else {
    elements.inputFilesHelp.textContent = "Select the files required by this module."
  }
  await updatePathlistFolderNotices();
  await refreshTaxonSelectors();
  await updatePreflight();
  resetGeneratedOutput();
  updateCommand();
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / 1024 ** 2).toFixed(1)} MB`;
}

async function updatePreflight() {
  const generation = ++preflightGeneration;
  let files = browserInputFiles();
  if (!files.length) {
    elements.preflightStatus.textContent = "Waiting for input";
    elements.preflightSummary.textContent = "Select input files to estimate browser load and validate basic FASTA structure.";
    return;
  }
  try {
    const expanded = [];
    for (const control of inputControls()) {
      if (generation !== preflightGeneration) return;
      const field = activeModule().inputFields.find(field => field.key === control.dataset.inputKey);
      const selected = filesForInput(control);
      if (field?.pathList && selected.some(file => /\.(pathlist|list)$/i.test(file.name))) {
        expanded.push(...await PhyloPrepInputs.resolveFiles([...selected, ...(inputFolderFiles.get(field.key) || [])], {
          sibling: async (entry, list) => {
            if (generation !== preflightGeneration) throw new Error("Preflight cancelled");
            return readPathlistSibling(field.key, entry, list);
          },
        }));
      } else expanded.push(...selected);
    }
    if (generation !== preflightGeneration) return;
    files = expanded;
  } catch (error) {
    if (generation !== preflightGeneration) return;
    elements.preflightStatus.textContent = "Input resolution needed";
    elements.preflightSummary.textContent = error.message;
    return;
  }
  const totalBytes = files.reduce((sum, file) => sum + (file.size || 0), 0);
  let fastaFiles = 0;
  let sequences = 0;
  let duplicateIds = 0;
  let malformedCodon = 0;
  const fastaIdSets = [];
  const sampled = files.filter(file => /\.(fa|faa|fas|fasta|aln|pep)$/i.test(file.name)).slice(0, 100);
  for (const file of sampled) {
    if (generation !== preflightGeneration) return;
    let text;
    try { text = await PhyloPrepInputs.readText(file); }
    catch (error) {
      if (generation !== preflightGeneration) return;
      elements.preflightStatus.textContent = "File unreadable";
      elements.preflightSummary.textContent = error.message;
      return;
    }
    if (generation !== preflightGeneration) return;
    if (!text.trimStart().startsWith(">")) continue;
    fastaFiles += 1;
    const ids = text.split(/\r?\n/).filter((line) => line.startsWith(">")).map((line) => line.slice(1).trim().split(/\s+/, 1)[0]).filter(Boolean);
    sequences += ids.length;
    duplicateIds += ids.length - new Set(ids).size;
    fastaIdSets.push(new Set(ids));
    if ((rawParameterValues().seqtype || "nucl") === "codon") {
      const chunks = text.split(/^>/m).slice(1);
      malformedCodon += chunks.filter((chunk) => {
        const lines = chunk.split(/\r?\n/);
        const seq = lines.slice(1).join("").replace(/\s+/g, "").replace(/-/g, "");
        return seq && seq.length % 3 !== 0;
      }).length;
    }
  }
  const problems = [];
  if (duplicateIds) problems.push(`<span class="preflight-error">${duplicateIds} duplicate FASTA ID(s)</span>`);
  if (malformedCodon) problems.push(`<span class="preflight-warning">${malformedCodon} sequence(s) are not complete codon multiples</span>`);
  const mappingProvided = Boolean(selectedInputFiles.get("mapping")?.length);
  if (activeModule().script === "phyloprep.py" && !mappingProvided && fastaIdSets.length > 1) {
    const sharedIds = [...fastaIdSets[0]].filter((id) => fastaIdSets.slice(1).every((set) => set.has(id)));
    if (!sharedIds.length) {
      problems.push('<span class="preflight-warning">No map.tsv selected and no FASTA IDs are shared across all loci. If IDs are locus-specific (for example GenBank/CDS IDs), provide FASTA_ID→Taxon mapping before building a supermatrix.</span>');
    }
  }
  const level = totalBytes > 100 * 1024 ** 2 ? "Large — CLI recommended" : totalBytes > 25 * 1024 ** 2 ? "Moderate browser workload" : "Browser-friendly";
  elements.preflightStatus.textContent = level;
  elements.preflightSummary.innerHTML = [
    `<strong>${files.length}</strong> file(s)`,
    files.length > 100 ? "Preflight samples up to 100 FASTA files; RUN processes all inputs" : "",
    `<strong>${formatBytes(totalBytes)}</strong>`,
    fastaFiles ? `<strong>${fastaFiles}</strong> FASTA file(s) · <strong>${sequences}</strong> sequence(s)` : "",
    problems.length ? problems.join(" · ") : '<span class="preflight-ok">Basic structure checks passed</span>',
  ].filter(Boolean).join(" · ");
}

async function chooseOutputFolder() {
  if (!("showDirectoryPicker" in window)) {
    setLog("This browser does not support local folder selection. Output will use its configured download folder.", "error");
    return;
  }
  try {
    outputDirectoryHandle = await window.showDirectoryPicker({ mode: "readwrite" });
    elements.outputFolder.textContent = `Selected local folder: ${outputDirectoryHandle.name}`;
    elements.downloadAll.textContent = "Save all to folder";
    updateCommand();
    setLog("Output folder selected. Generated files will be written there.", "success");
  } catch (error) {
    if (error.name !== "AbortError") setLog("Unable to select the output folder. Browser download folder remains in use.", "error");
  }
}

function cancellationError() {
  const error = new Error("Run cancelled by user.");
  error.code = "CANCELLED";
  return error;
}

function setRunControls(running) {
  elements.run.disabled = running;
  elements.run.textContent = running ? "Running…" : "Run";
  elements.cancelRun.hidden = !running;
}

function createPyodideWorker() {
  const worker = new Worker("pyodide-worker.js");
  worker.addEventListener("message", (event) => {
    const message = event.data;
    if (!activeWorkerRun || message.runId !== activeWorkerRun.id) return;
    if (message.type === "output-chunk" || message.type === "output-file-end") {
      clearInterval(activeWorkerRun.progressTimer);
      const outputs = activeWorkerRun.streamedOutputs;
      let file = outputs.get(message.index);
      if (!file) {
        file = { decoder: new TextDecoder(), parts: [], name: message.name };
        outputs.set(message.index, file);
      }
      if (message.type === "output-chunk") {
        file.parts.push(file.decoder.decode(message.bytes, { stream: true }));
        setRunStage(3, `Receiving output ${message.index + 1}/${message.fileCount}: ${message.name}`);
        worker.postMessage({ type: "output-chunk-ack", runId: message.runId, sequence: message.sequence });
      } else {
        file.parts.push(file.decoder.decode());
        file.content = file.parts.join("");
        file.parts = [];
        file.decoder = undefined;
      }
      return;
    }
    if (message.type === "stage") {
      if (message.stage === "runtime") {
        elements.runtimeStatus.textContent = `${message.label}…`;
        setRunStage(2, message.label);
      } else if (message.stage === "process") {
        elements.runtimeStatus.textContent = "Local worker running";
        setRunStage(2, message.label);
      }
      return;
    }
    if (message.type === "progress") {
      const progress = message.progress;
      const received = performance.now();
      const render = () => {
        const elapsed = progress.elapsed_seconds + (performance.now() - received) / 1000;
        const label = `${progress.completed}/${progress.total} ${progress.unit} · ${progress.label} · ${Math.round(elapsed)}s elapsed`;
        elements.runtimeStatus.textContent = label;
        setRunStage(2, label);
      };
      clearInterval(activeWorkerRun.progressTimer);
      render();
      activeWorkerRun.progressTimer = setInterval(render, 1000);
      return;
    }
    const run = activeWorkerRun;
    clearInterval(run.progressTimer);
    activeWorkerRun = undefined;
    elements.cancelRun.hidden = true;
    if (message.type === "result") {
      elements.runtimeStatus.textContent = "Local worker ready";
      try {
        if (message.result.files?.some(file => file.path)) {
          message.result.files = message.result.files.map((file, index) => {
            const streamed = run.streamedOutputs.get(index);
            if (!streamed || streamed.content === undefined) throw new Error(`Incomplete output received: ${file.name}`);
            return { name: file.name, content: streamed.content };
          });
        }
        run.resolve(message.result);
      } catch (error) {
        run.reject(error);
      }
    } else if (message.type === "error") {
      run.reject(new Error(message.message || "The local workflow could not be completed."));
    }
  });
  worker.addEventListener("error", (event) => {
    if (!activeWorkerRun) return;
    const run = activeWorkerRun;
    clearInterval(run.progressTimer);
    activeWorkerRun = undefined;
    elements.cancelRun.hidden = true;
    run.reject(new Error(event.message || "The local worker stopped unexpectedly."));
  });
  pyodideWorker = worker;
  return worker;
}

function runInWorker(payload) {
  const worker = pyodideWorker || createPyodideWorker();
  const runId = ++nextRunId;
  return new Promise((resolve, reject) => {
    activeWorkerRun = { id: runId, resolve, reject, streamedOutputs: new Map() };
    worker.postMessage({ type: payload.type || "run", runId, ...payload });
  });
}

function cancelActiveRun() {
  if (cancelAlignmentRun()) {
    elements.cancelRun.hidden = true;
    elements.runtimeStatus.textContent = "Alignment stopped";
    return;
  }
  if (!activeWorkerRun) return;
  const run = activeWorkerRun;
  clearInterval(run.progressTimer);
  activeWorkerRun = undefined;
  pyodideWorker?.terminate();
  pyodideWorker = undefined;
  elements.cancelRun.hidden = true;
  elements.runtimeStatus.textContent = "Local worker stopped";
  run.reject(cancellationError());
}

async function saveGeneratedFile(file) {
  if (outputDirectoryHandle) {
    const parts = file.name.replaceAll("\\", "/").split("/").filter(Boolean);
    let directory = outputDirectoryHandle;
    for (const part of parts.slice(0, -1)) {
      directory = await directory.getDirectoryHandle(part, { create: true });
    }
    const fileHandle = await directory.getFileHandle(parts.at(-1) || "output.txt", { create: true });
    const writable = await fileHandle.createWritable();
    await writable.write(file.content);
    await writable.close();
    return "selected local folder";
  }
  return "browser workspace";
}

function renderGeneratedFiles() {
  resultEntries = PhyloPrepResults.createEntries(generatedFiles);
  resultPageNumber = 0;
  elements.resultFolder.value = "all";
  elements.resultSearch.value = "";
  elements.resultControls.hidden = !generatedFiles.length;
  elements.resultCount.textContent = `${generatedFiles.length.toLocaleString()} output files`;
  elements.downloadSelected.hidden = !generatedFiles.length;
  elements.downloadAll.hidden = generatedFiles.length < 2;
  selectedGeneratedFileIndex = resultEntries[0]?.index ?? 0;
  renderGeneratedPage();
}

function renderGeneratedPage() {
  elements.outputFilesList.replaceChildren();
  const result = PhyloPrepResults.pageEntries(resultEntries, {
    folder: elements.resultFolder.value, query: elements.resultSearch.value, page: resultPageNumber,
  });
  resultPageNumber = result.page;
  elements.resultPrev.disabled = result.page === 0;
  elements.resultNext.disabled = result.page + 1 >= result.pages;
  elements.resultPage.textContent = result.total ? `${result.start + 1}–${result.start + result.entries.length} of ${result.total.toLocaleString()} matching files · page ${result.page + 1}/${result.pages}` : "No matching files";
  if (!result.entries.some(entry => entry.index === selectedGeneratedFileIndex) && result.entries.length) {
    selectedGeneratedFileIndex = result.entries[0].index;
  }
  result.entries.forEach(({ file, index }) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `output-file-button${index === selectedGeneratedFileIndex ? " active" : ""}`;
    const fileName = document.createElement("span");
    fileName.className = "output-file-name";
    fileName.textContent = file.name;
    const type = document.createElement("span");
    type.className = "output-file-type";
    const extension = file.name.split(".").pop()?.toLowerCase() || "file";
    type.textContent = ({ fa: "FASTA", fas: "FASTA", fasta: "FASTA", tsv: "TSV", nex: "NEXUS", nexus: "NEXUS", phy: "PHYLIP", phylip: "PHYLIP", axt: "AXT", txt: "TEXT" })[extension] || extension.toUpperCase();
    const destination = document.createElement("span");
    destination.className = "output-file-destination";
    destination.textContent = file.destination;
    button.append(fileName, type, destination);
    button.setAttribute("aria-label", `${file.name}, ${type.textContent}, saved to ${file.destination}`);
    button.title = file.name;
    button.addEventListener("click", () => {
      selectedGeneratedFileIndex = index;
      elements.preview.textContent = previewOutput(file.content);
      elements.outputFilesList.querySelectorAll("button").forEach((item) => item.classList.remove("active"));
      button.classList.add("active");
    });
    elements.outputFilesList.append(button);
  });
  elements.preview.textContent = result.entries.length ? previewOutput(generatedFiles[selectedGeneratedFileIndex]?.content || "") : "No files match the current filter.";
  elements.downloadSelected.hidden = !result.entries.length;
}

function previewOutput(content) {
  const limit = 200000;
  return content.length > limit ? content.slice(0, limit) + "\n[Preview truncated; download contains the complete file.]" : content;
}

async function downloadGeneratedFile(file) {
  const name = file.name.replaceAll("\\", "/").split("/").at(-1);
  if ("showSaveFilePicker" in window) {
    const handle = await window.showSaveFilePicker({ suggestedName: name });
    const writable = await handle.createWritable();
    try {
      await writable.write(file.content);
      await writable.close();
    } catch (error) {
      try { await writable.abort(); } catch {}
      throw error;
    }
    return "saved";
  }
  PhyloPrepDownloads.requestDownload(() => new Blob([file.content], { type: "text/plain;charset=utf-8" }), name, file);
  return "requested";
}

async function downloadSelectedOutput() {
  const file = generatedFiles[selectedGeneratedFileIndex];
  if (!file) return;
  try {
    const status = await downloadGeneratedFile(file);
    setLog(status === "saved" ? `Saved ${file.name}.` : `Download requested for ${file.name}; check your browser downloads for completion.`, "success");
  } catch (error) {
    setLog(error.name === "AbortError" ? "Save cancelled; results remain available." : `Unable to save ${file.name}: ${error.message}. Results remain available for retry.`, error.name === "AbortError" ? "warning" : "error");
  }
}

async function copyEquivalentCommand() {
  const command = elements.command.textContent;
  try {
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(command);
    } else {
      const temporary = document.createElement("textarea");
      temporary.value = command;
      document.body.append(temporary);
      temporary.select();
      document.execCommand("copy");
      temporary.remove();
    }
    elements.copyCommand.textContent = "Copied";
    window.setTimeout(() => { elements.copyCommand.textContent = "Copy"; }, 1400);
  } catch {
    elements.copyCommand.textContent = "Copy failed";
    window.setTimeout(() => { elements.copyCommand.textContent = "Copy"; }, 1800);
  }
}

async function downloadAllOutputs() {
  if (generatedFiles.length < 2) return;
  const stem = (elements.outputName.value.trim() || activeModule().script.replace(/\.py$/i, "outputs")).replace(/\.[^.]+$/, "").replace(/[^a-z0-9._-]+/gi, "-");
  try {
    if (outputDirectoryHandle) {
      for (const [index, file] of generatedFiles.entries()) {
        setLog(`Saving file ${index + 1}/${generatedFiles.length}: ${file.name}`, "info");
        await saveGeneratedFile(file);
      }
      setLog(`Saved ${generatedFiles.length} files to the selected output folder.`, "success");
      return;
    }
    PhyloPrepDownloads.checkZipCount(generatedFiles);
    const name = `${stem || "phyloprep"}-outputs.zip`;
    if ("showSaveFilePicker" in window) {
      const handle = await window.showSaveFilePicker({ suggestedName: name });
      const writable = await handle.createWritable();
      setLog("Writing ZIP archive to the selected file…", "info");
      await PhyloPrepDownloads.writeZip(generatedFiles, writable);
      setLog(`Saved ${generatedFiles.length} files as a ZIP archive.`, "success");
    } else {
      PhyloPrepDownloads.requestDownload(() => new Blob([...PhyloPrepDownloads.zipParts(generatedFiles)], { type: "application/zip" }), name, generatedFiles);
      setLog(`ZIP download requested for ${generatedFiles.length} files; check your browser downloads for completion.`, "success");
    }
  } catch (error) {
    setLog(error.name === "AbortError" ? "Save cancelled; results remain available." : `Unable to save outputs: ${error.message}. Results remain available for retry.`, error.name === "AbortError" ? "warning" : "error");
  }
}


function createAlignmentWorker() {
  const worker = new Worker("alignment-worker.js");
  worker.addEventListener("message", (event) => {
    if (!activeAlignmentRun || event.data.runId !== activeAlignmentRun.id) return;
    const run = activeAlignmentRun;
    activeAlignmentRun = undefined;
    if (event.data.type === "result") run.resolve(event.data.alignment);
    else run.reject(new Error(event.data.message || "Alignment failed."));
  });
  worker.addEventListener("error", (event) => {
    if (!activeAlignmentRun) return;
    const run = activeAlignmentRun;
    activeAlignmentRun = undefined;
    run.reject(new Error(event.message || "Alignment worker stopped unexpectedly."));
  });
  alignmentWorker = worker;
  return worker;
}

function runBrowserAlignment(sequenceText, tool, seqtype) {
  const worker = alignmentWorker || createAlignmentWorker();
  const runId = ++alignmentRunCounter;
  return new Promise((resolve, reject) => {
    activeAlignmentRun = { id: runId, resolve, reject };
    worker.postMessage({ type: "align", runId, sequenceText, tool, seqtype });
  });
}

function cancelAlignmentRun() {
  if (!activeAlignmentRun) return false;
  const run = activeAlignmentRun;
  activeAlignmentRun = undefined;
  alignmentWorker?.terminate();
  alignmentWorker = undefined;
  run.reject(cancellationError());
  return true;
}


async function runMsapBrowserWorkflow(request) {
  setRunStage(2, "Preparing input with original MSAP.py");
  elements.runtimeStatus.textContent = "Preparing MSAP.py workflow";
  const probe = await runInWorker({ type: "msap-probe", ...request });
  if (!probe.alignmentInput) throw new Error("MSAP.py did not provide an input for the selected browser aligner.");

  const alignerLabel = probe.tool === "mafft" ? "MAFFT E-INS-i 7.520 (Biowasm)" : "MUSCLE 5.1.0";
  elements.cancelRun.hidden = false;
  setRunStage(2, `Running Biowasm ${alignerLabel}`);
  elements.runtimeStatus.textContent = `Biowasm ${alignerLabel} running`;
  const alignmentContent = await runBrowserAlignment(probe.alignmentInput, probe.tool, request.options.seqtype || "codon");

  setRunStage(2, "Finishing original MSAP.py workflow");
  elements.runtimeStatus.textContent = "Finishing MSAP.py workflow";
  elements.cancelRun.hidden = false;
  return runInWorker({ type: "msap-complete", ...request, alignmentContent });
}


function fileStem(name) {
  const leaf = String(name).replaceAll("\\", "/").split("/").at(-1) || "";
  const dot = leaf.lastIndexOf(".");
  return dot > 0 ? leaf.slice(0, dot) : leaf;
}

function requireUniqueStems(files) {
  const stems = files.map((file) => fileStem(file.name));
  const duplicate = stems.find((stem, index) => stems.indexOf(stem) !== index);
  if (duplicate) throw new Error(`Input FASTA basenames must be unique; duplicate stem: ${duplicate}`);
}

function mergeTsvContents(files) {
  let header;
  const rows = [];
  files.forEach((file) => {
    const lines = String(file.content || "").trim().split(/\r?\n/).filter(Boolean);
    if (!lines.length) return;
    if (!header) header = lines[0];
    rows.push(...lines.slice(1));
  });
  return header ? `${header}\n${rows.join("\n")}${rows.length ? "\n" : ""}` : "";
}

function selectedInputPathList(inputs) {
  return `${inputs.map((file) => virtualInputPath(file.name)).join("\n")}\n`;
}

function batchPathLists(files) {
  const suffixes = ["nucl.aln", "prot.aln", "codon.aln", "nucl.trimmed.aln", "prot.trimmed.aln", "codon.trimmed.aln"];
  const outputs = [];
  suffixes.forEach((suffix) => {
    const names = files.filter((file) => file.name.endsWith(`.${suffix}`)).map((file) => file.name).sort();
    if (names.length) outputs.push({ name: `path-lists/all.${suffix}.pathlist`, content: `${names.join("\n")}\n` });
  });
  return outputs;
}

function showLocusProgress(completed, total, started, label) {
  const elapsed = (performance.now() - started) / 1000;
  const eta = completed ? ` · estimated remaining ${Math.round(elapsed * (total - completed) / completed)}s` : "";
  const text = `${completed}/${total} loci complete · ${label} · ${Math.round(elapsed)}s elapsed${eta}`;
  elements.runtimeStatus.textContent = text;
  setRunStage(2, text);
}

async function runMsapBatchBrowserWorkflow(request) {
  const inputs = request.files.filter((file) => file.inputKey === "inputs" && String(file.content || "").trimStart().startsWith(">"));
  if (!inputs.length) throw new Error("Select at least one FASTA input.");
  requireUniqueStems(inputs);
  const outputs = [];
  if (inputs.length > 1) {
    outputs.push({
      name: "path-lists/inputs.pathlist",
      content: selectedInputPathList(inputs),
    });
  }
  const progressStarted = performance.now();
  for (let index = 0; index < inputs.length; index += 1) {
    const input = inputs[index];
    showLocusProgress(index, inputs.length, progressStarted, `Aligning ${input.name}`);
    const heartbeat = setInterval(() => showLocusProgress(index, inputs.length, progressStarted, `Aligning ${input.name}`), 1000);
    let result;
    try {
      result = await runMsapBrowserWorkflow({
        ...request, script: "MSAP.py", files: [input], outputName: "", usesPathList: false,
      });
    } finally {
      clearInterval(heartbeat);
    }
    showLocusProgress(index + 1, inputs.length, progressStarted, input.name);
    outputs.push(...result.files.map((file) => ({ name: file.name, content: file.content })));
  }
  outputs.push(...batchPathLists(outputs));
  const changeReports = outputs.filter((file) => file.name.endsWith(".sequence_changes.tsv") && !file.name.includes("/"));
  const mergedReport = mergeTsvContents(changeReports);
  if (mergedReport) outputs.push({ name: "reports/sequence_changes.tsv", content: mergedReport });
  return {
    files: outputs,
    summary: `Processed ${inputs.length} loci with the browser MSAP batch workflow.`,
    provenance: `Browser MSAP_batch adapter · untouched MSAP.py + Biowasm ${request.options.align_software === "mafft" ? "MAFFT E-INS-i 7.520 (Biowasm)" : "MUSCLE 5.1.0"}`,
  };
}

async function runPhyloPrepBrowserWorkflow(request) {
  const sources = request.files.filter((file) => file.inputKey === "inputs" && String(file.content || "").trimStart().startsWith(">"));
  if (!sources.length) throw new Error("Select at least one FASTA input.");
  requireUniqueStems(sources);
  const mapping = request.files.find((file) => file.inputKey === "mapping");
  const batch = await runMsapBatchBrowserWorkflow({ ...request, files: sources });
  const batchProducts = batch.files
    .filter((file) => !file.name.includes("/") && (file.name.endsWith(".aln") || file.name.endsWith(".sequence_changes.tsv")))
    .map((file) => ({ ...file, role: "batch-output" }));
  const finalFiles = [
    ...sources.map((file) => ({ ...file, role: "source" })),
    ...(mapping ? [{ ...mapping, role: "mapping" }] : []),
    ...batchProducts,
  ];
  setRunStage(2, "Running PhyloPrep QC and matrix construction");
  elements.runtimeStatus.textContent = "PhyloPrep matrix construction";
  const result = await runInWorker({
    type: "phyloprep-finalize",
    script: "phyloprep.py",
    files: finalFiles,
    options: request.options,
    outputName: "",
    usesPathList: false,
  });
  if (sources.length > 1) {
    result.files.push({
      name: "inputs.pathlist",
      content: selectedInputPathList(sources),
    });
  }
  return result;
}

async function sha256Text(text) {
  const digest = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
  return [...new Uint8Array(digest)].map((byte) => byte.toString(16).padStart(2, "0")).join("");
}

async function buildRunManifest(module, request, result, outputs) {
  const inputManifest = [];
  for (const file of request.files) {
    inputManifest.push({
      name: file.name,
      input_key: file.inputKey,
      bytes: new TextEncoder().encode(file.content).length,
      sha256: await sha256Text(file.content),
    });
  }
  const outputManifest = [];
  for (const file of outputs) {
    outputManifest.push({
      name: file.name,
      bytes: new TextEncoder().encode(file.content).length,
      sha256: await sha256Text(file.content),
    });
  }
  return JSON.stringify({
    schema: "phyloprep-web-run-v1",
    created_at: new Date().toISOString(),
    module: module.script,
    implementation: implementationProfile(module),
    equivalent_cli: elements.command.textContent,
    parameters: request.options,
    inputs: inputManifest,
    outputs: outputManifest,
    provenance: result.provenance || null,
    runtime: {
      pyodide: "0.27.2",
      aligner: ["phyloprep.py", "MSAP_batch.py"].includes(module.script)
        ? (request.options.align_software === "mafft" ? "Biowasm MAFFT E-INS-i 7.520" : "Biowasm MUSCLE 5.1.0")
        : null,
      user_agent: navigator.userAgent,
    },
  }, null, 2) + "\n";
}

function transferTarget(module) {
  if (["oma_to_orthogroups.py", "orthomcl_to_orthogroups.py", "proteinortho_to_orthogroups.py", "sonicparanoid_to_orthogroups.py"].includes(module.script)) {
    return "extract_orthofinder_orthogroups.py";
  }
  if (module.script === "extract_orthofinder_orthogroups.py") return "phyloprep.py";
  return null;
}

function prepareNextStep() {
  const target = transferTarget(activeModule());
  if (!target || !generatedFiles.length) return;
  const inputs = new Map();
  if (target === "extract_orthofinder_orthogroups.py") {
    const table = generatedFiles.find((file) => /\.tsv$/i.test(file.name) && file.name !== "phyloprep-run.json");
    if (!table) return;
    inputs.set("orthogroups", [new File([table.content], table.name, { type: "text/tab-separated-values" })]);
  } else if (target === "phyloprep.py") {
    const fastas = generatedFiles.filter((file) => /\.(fa|fas|fasta)$/i.test(file.name));
    const mapping = generatedFiles.find((file) => /(^|\/)map\.tsv$/i.test(file.name));
    if (!fastas.length) return;
    inputs.set("inputs", fastas.map((file) => new File([file.content], file.name.split("/").at(-1), { type: "text/plain" })));
    if (mapping) inputs.set("mapping", [new File([mapping.content], "map.tsv", { type: "text/tab-separated-values" })]);
  }
  const index = modules.findIndex((module) => module.script === target);
  if (index < 0) return;
  pendingTransfer = { script: target, inputs };
  selectModule(index);
}

async function readPathlistSibling(inputKey, entry, list) {
  const root = inputFolderHandles.get(inputKey);
  if (!root) return undefined;
  const normalized = PhyloPrepInputs.normalize(entry);
  const paths = [];
  if (!normalized.startsWith("/") && !/^[A-Za-z]:/.test(normalized)) paths.push(normalized);
  paths.push(PhyloPrepInputs.basename(entry));
  for (const path of [...new Set(paths)]) {
    const parts = path.split("/").filter(part => part && part !== ".");
    if (parts.includes("..")) continue;
    try {
      let directory = root;
      for (const part of parts.slice(0, -1)) directory = await directory.getDirectoryHandle(part);
      const cacheKey = `${inputKey}\u0000${parts.join("/")}`;
      if (resolvedFolderFiles.has(cacheKey)) return resolvedFolderFiles.get(cacheKey);
      const file = await (await directory.getFileHandle(parts.at(-1))).getFile();
      resolvedFolderFiles.set(cacheKey, file);
      return file;
    } catch (error) {
      if (error.name !== "NotFoundError" && error.name !== "TypeMismatchError") {
        throw new Error(`Cannot read ${entry} from input folder: ${error.message}. Reselect the input folder.`);
      }
    }
  }
  return undefined;
}

async function runBrowserModule() {
  const module = activeModule();
  const files = browserInputFiles();
  const missingInput = inputControls().find((control) => control.dataset.inputRequired !== "false" && !filesForInput(control).length);
  if (missingInput) {
    setLog(`Select ${missingInput.dataset.inputLabel} before running.`, "error");
    setRunStage(1, "Input required", true);
    missingInput.focus();
    return;
  }
  setRunControls(true);
  elements.cancelRun.hidden = true;
  let stage = 2;
  setRunStage(1, "Reading input files");
  setLog("Preparing input for the local worker…");
  try {
    ++preflightGeneration;
    resolvedFolderFiles.clear();
    const filePayload = [];
    for (const control of inputControls()) {
      const field = module.inputFields.find(field => field.key === control.dataset.inputKey);
      let inputs = filesForInput(control);
      if (field?.pathList) {
        const available = [...inputs, ...(inputFolderFiles.get(field.key) || [])];
        inputs = await PhyloPrepInputs.resolveFiles(available, {
          sibling: (entry, list) => readPathlistSibling(field.key, entry, list),
        });
      }
      const contents = await PhyloPrepInputs.readBounded(inputs, async file => {
        try { return await PhyloPrepInputs.readText(file); }
        catch (error) {
          if (!field?.pathList || !inputFolderHandles.has(field.key)) throw error;
          const fresh = await readPathlistSibling(field.key, file.name, file);
          if (!fresh) throw error;
          return PhyloPrepInputs.readText(fresh);
        }
      });
      for (let index = 0; index < inputs.length; index++) {
        filePayload.push({ name: inputs[index].name, content: contents[index], inputKey: control.dataset.inputKey });
      }
    }
    const optionValues = rawParameterValues();
    setRunStage(stage, "Loading local workflow");
    elements.cancelRun.hidden = false;
    const request = {
      script: module.script,
      files: filePayload,
      options: optionValues,
      outputName: elements.outputName.value.trim() || module.defaultOutput,
      usesPathList: module.inputFields.some((field) => field.pathList),
    };
    const result = module.browserWorkflow === "phyloprep"
      ? await runPhyloPrepBrowserWorkflow(request)
      : module.browserWorkflow === "batch"
        ? await runMsapBatchBrowserWorkflow(request)
        : module.browserWorkflow === "orthofinder-extract"
          ? await runInWorker({ type: "orthofinder-extract", ...request })
          : await runInWorker(request);
    stage = 3;
    setRunStage(stage, "Saving output");
    generatedFiles = [];
    const resultOutputs = result.files.map((output) => ({ name: output.name, content: output.content }));
    const manifest = await buildRunManifest(module, request, result, resultOutputs);
    resultOutputs.push({ name: "phyloprep-run.json", content: manifest });
    for (const output of resultOutputs) {
      const generated = { name: output.name, content: output.content };
      generated.destination = await saveGeneratedFile(generated);
      generatedFiles.push(generated);
    }
    renderGeneratedFiles();
    elements.useNextStep.hidden = !transferTarget(module);
    if (result.partial) {
      setRunStage(3, "Needs map.tsv or missing-taxa adjustment", true);
      setLog(result.user_error || result.summary || "PhyloPrep needs input adjustment before a supermatrix can be completed.", "warning");
    } else {
      setRunStage(4, "Completed");
      setLog([result.summary || `Created ${generatedFiles.length} output file(s).`, result.provenance].filter(Boolean).join(" "), "success");
    }
  } catch (error) {
    const cancelled = error.code === "CANCELLED";
    setRunStage(stage, cancelled ? "Run cancelled" : "Workflow needs attention", !cancelled);
    setLog(cancelled ? "Run cancelled. You can change inputs or parameters and try again." : (error.message || String(error)), cancelled ? "" : "error");
  } finally {
    setRunControls(false);
  }
}

function refreshCatalogSearch() {
  catalogQuery = elements.scriptSearch.value;
  renderCatalog();
}

elements.scriptSearch.addEventListener("input", refreshCatalogSearch);
elements.scriptSearch.addEventListener("search", refreshCatalogSearch);
elements.useNextStep.addEventListener("click", prepareNextStep);
elements.clearInputFiles.addEventListener("click", clearAllInputFiles);
elements.chooseFolder.addEventListener("click", chooseOutputFolder);
elements.outputName.addEventListener("input", updateCommand);
elements.run.addEventListener("click", runBrowserModule);
elements.cancelRun.addEventListener("click", cancelActiveRun);
elements.copyCommand.addEventListener("click", copyEquivalentCommand);
elements.downloadSelected.addEventListener("click", downloadSelectedOutput);
elements.downloadAll.addEventListener("click", downloadAllOutputs);
for (const control of [elements.resultFolder, elements.resultSearch]) {
  control.addEventListener("input", () => { resultPageNumber = 0; renderGeneratedPage(); });
}
elements.resultPrev.addEventListener("click", () => { resultPageNumber--; renderGeneratedPage(); });
elements.resultNext.addEventListener("click", () => { resultPageNumber++; renderGeneratedPage(); });

if ("serviceWorker" in navigator) {
  window.addEventListener("load", () => {
    navigator.serviceWorker.register("./sw.js").catch(() => {
      // The site remains fully usable when a browser or local preview blocks service workers.
    });
  });
}

selectModule(activeModuleIndex);
