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
  downloadSelected: document.querySelector("#download-selected"),
  downloadAll: document.querySelector("#download-all"),
  preview: document.querySelector("#output-preview"),
  command: document.querySelector("#command-example"),
  copyCommand: document.querySelector("#copy-command"),
  runStages: document.querySelector("#run-stages"),
  runStageLabel: document.querySelector("#run-stage-label"),
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

const modules = [
  {
    script: "MSAP.py", title: "Complete alignment workflow", category: "Workflow", browser: true, fullWorkflow: true,
    summary: "Runs the untouched single-locus MSAP.py workflow locally with MUSCLE 5.1.0 supplied by Biowasm/WebAssembly.",
    inputs: ["Input FASTA (-i, --input)"],
    inputFields: [{ key: "input", label: "Input FASTA", flag: "-i, --input" }],
    outputs: ["Preprocessed FASTA and audit report", "MUSCLE alignment", "Codon back-translation for CDS input", "Optional PhyloPrep trimming outputs"],
    outputLabel: "Output naming", outputHelp: "MSAP.py derives output names from the input FASTA stem.", defaultOutput: "", hideOutputName: true,
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Browser workflow supports codon, protein, and nucleotide modes. Pseudogene/MACSE remains CLI-only.", type: "select", value: "codon", options: [["codon", "Codon / CDS"], ["prot", "Protein"], ["nucl", "Nucleotide"]], token: "seqtype" },
      geneticCodeParameter(),
      { flag: "-n", label: "Skip trimming", description: "Keep the raw MUSCLE alignment(s) without running trimAlnSeq.py.", type: "checkbox", value: false, token: "notrim" },
      { flag: "-G", label: "Maximum gap proportion", description: "Maximum gap ratio per site for PhyloPrep's bundled trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "gap_ratio" },
      { flag: "-N", label: "Maximum N proportion", description: "Maximum N ratio per site for nucleotide/codon trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "n_ratio" },
      { flag: "-X", label: "Maximum X proportion", description: "Maximum X ratio per site for protein trimming. Range: 0–1.", type: "number", value: "0.2", step: "0.01", token: "x_ratio" },
    ],
    command: "python MSAP.py -i {input} -st {seqtype} -s muscle -t 1 -g {genetic_code} -ts trimAlnSeq -G {gap_ratio} -N {n_ratio} -X {x_ratio}{notrim_option}",
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
    script: "alignment_qc.py", title: "Alignment quality control", category: "QC", browser: true,
    summary: "Checks FASTA alignments locally and creates pass, fail, and detailed QC reports.",
    inputs: ["FASTA files and/or pathlist files (one FASTA path per line) (-i, --input)", "Required report prefix (-p, --prefix)", "Sequence type (-st, --seqtype)"],
    inputFields: [{ key: "inputs", label: "Alignment input files", flag: "-i, --input", multiple: true, pathList: true }],
    outputs: ["PREFIX.pass.tsv", "PREFIX.fail.tsv", "PREFIX.details.tsv"],
    outputLabel: "Report prefix", outputHelp: "The command-line script creates its QC reports from this prefix.", defaultOutput: "phyloprep-qc",
    parameters: [
      { flag: "-st", label: "Sequence type", description: "Expected sequence type for the command-line quality-control workflow.", type: "select", value: "nucl", options: [["nucl", "Nucleotide"], ["prot", "Protein"], ["codon", "Codon"]], token: "seqtype" },
      geneticCodeParameter(),
      { flag: "--min-length", label: "Minimum length", description: "Minimum alignment or sequence length. Default: 1.", type: "number", value: "1" },
      { flag: "--max-n-ratio", label: "Maximum N ratio", description: "Maximum N-base ratio for nucleotide matrices or N-containing codon ratio in codon matrices. Range: 0–1.", type: "number", value: "1", step: "0.01", token: "max_n_ratio" },
      { flag: "--max-x-ratio", label: "Maximum X ratio", description: "Maximum X ratio for protein matrices. Range: 0–1.", type: "number", value: "1", step: "0.01", token: "max_x_ratio" },
      { flag: "--max-gap-ratio", label: "Maximum gap ratio", description: "Proportion allowed before the alignment is rejected. Default: 1.", type: "number", value: "1", step: "0.01" },
      { flag: "--min-parsimony-informative-sites", label: "Minimum parsimony-informative sites", description: "Minimum count required in the complete alignment matrix. Default: 0.", type: "number", value: "0", token: "min_parsimony_informative_sites" },
      { flag: "--allow-terminal-stop", label: "Allow terminal stop codons", description: "Permit terminal stop codons in codon-mode QC.", type: "checkbox", value: false, token: "allow_terminal_stop" },
      { flag: "--allow-internal-stop", label: "Allow internal stop codons", description: "Permit internal stop codons in codon-mode QC.", type: "checkbox", value: false, token: "allow_internal_stop" },
    ],
    command: "python alignment_qc.py -i {inputs} -p {prefix} -st {seqtype} -g {genetic_code} --min-length {min_length} --max-n-ratio {max_n_ratio} --max-x-ratio {max_x_ratio} --max-gap-ratio {max_gap_ratio} --min-parsimony-informative-sites {min_parsimony_informative_sites}{allow_terminal_stop}{allow_internal_stop}",
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
    parameters: [{ flag: "--missing-taxa", label: "Missing taxa policy", description: "How records absent from a gene alignment are handled.", type: "select", value: "skip-gene", options: [["skip-gene", "Skip gene"], ["pad-gaps", "Pad with gaps"]], token: "missing_taxa" }],
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
    outputLabel: "Output file name", outputHelp: "Trimmed FASTA written by -o when provided.", defaultOutput: "trimmed.fasta",
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
let pyodideWorker;
let activeWorkerRun;
let nextRunId = 0;
let selectedInputFiles = new Map();
let muscleCliPromise;
let muscleRunCounter = 0;

function activeModule() { return modules[activeModuleIndex]; }

function implementationProfile(module) {
  if (module.fullWorkflow) {
    return {
      key: "compatible",
      label: "Original MSAP.py + Biowasm MUSCLE",
      note: "Runs the repository's untouched MSAP.py in Pyodide. Its MUSCLE subprocess call is fulfilled locally by Biowasm MUSCLE 5.1.0; bundled AA2Codon.py and trimAlnSeq.py continue inside the same browser Python runtime. No sequence file is uploaded to PhyloPrep.",
    };
  }
  return {
    key: "original",
    label: "Original script in Pyodide",
    note: `Runs the repository's original ${module.script} source locally in the browser, with the same CLI implementation and bundled dependencies.`,
  };
}

function renderCatalog() {
  elements.scriptList.replaceChildren();
  modules.forEach((module, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `script-button${index === activeModuleIndex ? " active" : ""}`;
    button.setAttribute("role", "tab");
    button.setAttribute("aria-selected", String(index === activeModuleIndex));
    const name = document.createElement("span");
    name.className = "script-name";
    name.textContent = module.script;
    const category = document.createElement("span");
    category.className = "script-category";
    category.textContent = module.category;
    button.append(name, category);
    button.addEventListener("click", () => selectModule(index));
    elements.scriptList.append(button);
  });
}

function parameterToken(parameter) {
  return parameter.token || parameter.flag.replace(/[^a-z0-9]+/gi, "_").replace(/^_|_$/g, "").toLowerCase();
}

function renderParameterControl(parameter) {
  const card = document.createElement("section");
  card.className = "parameter-card";
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
      control.addEventListener("input", updateCommand);
      control.addEventListener("change", updateCommand);
    });
  }
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
  return `${file.name}\u0000${file.size}\u0000${file.lastModified}`;
}

function rememberInputFiles(control) {
  const selected = [...control.files];
  const existing = filesForInput(control);
  const files = control.multiple
    ? [...existing, ...selected].filter((file, index, all) => all.findIndex((candidate) => fileIdentity(candidate) === fileIdentity(file)) === index)
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

function virtualInputPath(name) {
  return `/work/input/${name}`;
}

function pathListPath(module = activeModule()) {
  return `/work/pathlists/${module.script.replace(/\.py$/i, "")}.pathlist`;
}

function renderInputFileFields(module) {
  elements.inputFileFields.replaceChildren();
  selectedInputFiles = new Map();
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
    heading.append(label, flag);
    const picker = document.createElement("label");
    picker.className = "folder-picker input-file-picker";
    picker.htmlFor = id;
    const control = document.createElement("input");
    control.className = "module-input input-file-control";
    control.id = id;
    control.dataset.inputKey = field.key;
    control.dataset.inputLabel = field.label;
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
    elements.inputFileFields.append(section);
  });
}

function cliQuote(value) {
  const text = String(value);
  return /[\s"']/u.test(text) ? `"${text.replaceAll("\\", "\\\\").replaceAll("\"", "\\\"")}"` : text;
}

function valuesForCommand() {
  const module = activeModule();
  const inputNames = selectedInputNames();
  const parameterValues = {};
  module.parameters.forEach((parameter) => {
    const control = elements.parameterList.querySelector(`[data-parameter="${parameterToken(parameter)}"]`);
    parameterValues[parameterToken(parameter)] = control?.type === "checkbox"
      ? control.checked
      : control?.multiple
        ? [...control.selectedOptions].map((option) => option.value)
        : (control?.value || parameter.value || "");
  });
  const rawOutput = elements.outputName.value.trim() || module.defaultOutput;
  const outputPath = outputDirectoryHandle ? `${outputDirectoryHandle.name}/${rawOutput}` : rawOutput;
  const usePathList = module.inputFields.some((field) => field.pathList) && inputNames.length > 1;
  const inputPaths = inputNames.map((name) => cliQuote(virtualInputPath(name)));
  const mappingName = inputNames.find((name, index) => index > 0 && /\.(tsv|txt|csv)$/i.test(name)) || inputNames[1] || "taxa-map.tsv";
  const mapStem = rawOutput.replace(/\.[^.]+$/, "") || "renamed";
  const mapOutputPath = outputDirectoryHandle ? `${outputDirectoryHandle.name}/${mapStem}.map.tsv` : `${mapStem}.map.tsv`;
  return {
    ...parameterValues,
    input: inputPaths[0] || "/work/input/input.fasta",
    inputs: usePathList ? pathListPath(module) : (inputPaths.length ? inputPaths.join(" ") : "/work/input/input.fasta"),
    cds: inputPaths[0] || "/work/input/cds.fasta",
    protaln: inputPaths[1] || "/work/input/protein.aln.fasta",
    output: cliQuote(outputPath),
    prefix: cliQuote(outputPath),
    seqtype: parameterValues.seqtype || "nucl",
    genetic_code: parameterValues.genetic_code || "1",
    allow_terminal_stop: parameterValues.allow_terminal_stop ? " --allow-terminal-stop" : "",
    allow_internal_stop: parameterValues.allow_internal_stop ? " --allow-internal-stop" : "",
    removeall: parameterValues.removeall ? " -ra" : "",
    keep_taxa_option: parameterValues.keep_taxa?.length ? " --keep-taxa /work/path-lists/filter_sequences.keep-taxa.txt" : "",
    remove_taxa_option: parameterValues.remove_taxa?.length ? " --remove-taxa /work/path-lists/filter_sequences.remove-taxa.txt" : "",
    mapping_option: inputNames.length > 1 ? ` -m ${cliQuote(virtualInputPath(mappingName))}` : "",
    map_output_option: parameterValues.map_output ? ` --map-output ${cliQuote(mapOutputPath)}` : "",
    allow_duplicate_ids_option: parameterValues.allow_duplicate_ids ? " --allow-duplicate-ids" : "",
    keep_special_option: parameterValues.keep_special ? " --keep-special" : "",
    notrim_option: parameterValues.notrim ? " -n" : "",
  };
}

function updateCommand() {
  const values = valuesForCommand();
  elements.command.textContent = activeModule().command.replace(/\{([a-z_]+)\}/g, (_, key) => values[key] ?? `{${key}}`);
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

async function renderInputFiles() {
  const controls = inputControls();
  const names = selectedInputNames();
  elements.inputFilesList.replaceChildren();
  controls.forEach((control) => {
    const selected = filesForInput(control);
    const display = control.closest(".input-file-picker")?.querySelector(".input-file-display");
    if (!display) return;
    display.textContent = selected.length
      ? selected.map((file) => file.name).join(", ")
      : (control.multiple ? "Select local files" : "Select a local file");
    display.title = selected.map((file) => file.name).join(", ");
  });
  if (!names.length) {
    elements.inputFilesList.textContent = "No files selected.";
  } else {
    controls.forEach((control) => {
      filesForInput(control).forEach((file, index) => {
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
  const module = activeModule();
  const pathLists = browserInputFiles().filter((file) => /\.(pathlist|list)$/i.test(file.name));
  const fastaFiles = browserInputFiles().filter((file) => !/\.(pathlist|list)$/i.test(file.name));
  if (module.inputFields.some((field) => field.pathList) && names.length > 1) {
    const browserNote = pathLists.length && !fastaFiles.length
      ? " Upload the referenced FASTA files as well to run them in the browser."
      : "";
    elements.inputFilesHelp.textContent = `Multiple selections are merged into a browser pathlist for the CLI preview: ${pathListPath(module)}.${browserNote}`;
  } else if (module.inputFields.some((field) => field.pathList) && pathLists.length) {
    elements.inputFilesHelp.textContent = "A .pathlist file is selected. Upload its referenced FASTA files as well for browser execution.";
  } else {
    elements.inputFilesHelp.textContent = "Select the files required by this module."
  }
  await refreshTaxonSelectors();
  resetGeneratedOutput();
  updateCommand();
}

async function chooseOutputFolder() {
  if (!("showDirectoryPicker" in window)) {
    setLog("This browser does not support local folder selection. Output will use its configured download folder.", "error");
    return;
  }
  try {
    outputDirectoryHandle = await window.showDirectoryPicker({ mode: "readwrite" });
    elements.outputFolder.textContent = `Selected local folder: ${outputDirectoryHandle.name}`;
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
    const run = activeWorkerRun;
    activeWorkerRun = undefined;
    elements.cancelRun.hidden = true;
    if (message.type === "result") {
      elements.runtimeStatus.textContent = "Local worker ready";
      run.resolve(message.result);
    } else if (message.type === "error") {
      run.reject(new Error(message.message || "The local workflow could not be completed."));
    }
  });
  worker.addEventListener("error", (event) => {
    if (!activeWorkerRun) return;
    const run = activeWorkerRun;
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
    activeWorkerRun = { id: runId, resolve, reject };
    worker.postMessage({ type: payload.type || "run", runId, ...payload });
  });
}

function cancelActiveRun() {
  if (!activeWorkerRun) return;
  const run = activeWorkerRun;
  activeWorkerRun = undefined;
  pyodideWorker?.terminate();
  pyodideWorker = undefined;
  elements.cancelRun.hidden = true;
  elements.runtimeStatus.textContent = "Local worker stopped";
  run.reject(cancellationError());
}

async function saveGeneratedFile(file) {
  if (outputDirectoryHandle) {
    const fileHandle = await outputDirectoryHandle.getFileHandle(file.name, { create: true });
    const writable = await fileHandle.createWritable();
    await writable.write(file.content);
    await writable.close();
    return "selected local folder";
  }
  return "browser workspace";
}

function renderGeneratedFiles() {
  elements.outputFilesList.replaceChildren();
  selectedGeneratedFileIndex = 0;
  elements.resultCount.textContent = generatedFiles.length === 1 ? "1 output file" : `${generatedFiles.length} output files`;
  elements.downloadSelected.hidden = !generatedFiles.length;
  elements.downloadAll.hidden = generatedFiles.length < 2;
  generatedFiles.forEach((file, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `output-file-button${index === 0 ? " active" : ""}`;
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
      elements.preview.textContent = file.content;
      elements.outputFilesList.querySelectorAll("button").forEach((item) => item.classList.remove("active"));
      button.classList.add("active");
    });
    elements.outputFilesList.append(button);
  });
  elements.preview.textContent = generatedFiles[0]?.content || "";
}

function downloadGeneratedFile(file) {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(new Blob([file.content], { type: "text/plain;charset=utf-8" }));
  link.download = file.name;
  link.click();
  window.setTimeout(() => URL.revokeObjectURL(link.href), 1000);
}

function downloadSelectedOutput() {
  const file = generatedFiles[selectedGeneratedFileIndex];
  if (!file) return;
  downloadGeneratedFile(file);
  setLog(`Downloaded ${file.name}.`, "success");
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

let crcTable;

function crc32(bytes) {
  if (!crcTable) {
    crcTable = Array.from({ length: 256 }, (_, index) => {
      let value = index;
      for (let bit = 0; bit < 8; bit += 1) value = value & 1 ? (value >>> 1) ^ 0xedb88320 : value >>> 1;
      return value >>> 0;
    });
  }
  let value = 0xffffffff;
  bytes.forEach((byte) => { value = (value >>> 8) ^ crcTable[(value ^ byte) & 0xff]; });
  return (value ^ 0xffffffff) >>> 0;
}

function makeZip(files) {
  const encoder = new TextEncoder();
  const chunks = [];
  const centralDirectory = [];
  let offset = 0;
  files.forEach((file) => {
    const name = encoder.encode(file.name);
    const content = encoder.encode(file.content);
    const checksum = crc32(content);
    const local = new Uint8Array(30 + name.length);
    const localView = new DataView(local.buffer);
    localView.setUint32(0, 0x04034b50, true);
    localView.setUint16(4, 20, true);
    localView.setUint32(14, checksum, true);
    localView.setUint32(18, content.length, true);
    localView.setUint32(22, content.length, true);
    localView.setUint16(26, name.length, true);
    local.set(name, 30);
    chunks.push(local, content);

    const central = new Uint8Array(46 + name.length);
    const centralView = new DataView(central.buffer);
    centralView.setUint32(0, 0x02014b50, true);
    centralView.setUint16(4, 20, true);
    centralView.setUint16(6, 20, true);
    centralView.setUint32(16, checksum, true);
    centralView.setUint32(20, content.length, true);
    centralView.setUint32(24, content.length, true);
    centralView.setUint16(28, name.length, true);
    centralView.setUint32(42, offset, true);
    central.set(name, 46);
    centralDirectory.push(central);
    offset += local.length + content.length;
  });
  const centralOffset = offset;
  const centralSize = centralDirectory.reduce((size, item) => size + item.length, 0);
  const end = new Uint8Array(22);
  const endView = new DataView(end.buffer);
  endView.setUint32(0, 0x06054b50, true);
  endView.setUint16(8, files.length, true);
  endView.setUint16(10, files.length, true);
  endView.setUint32(12, centralSize, true);
  endView.setUint32(16, centralOffset, true);
  return new Blob([...chunks, ...centralDirectory, end], { type: "application/zip" });
}

function downloadAllOutputs() {
  if (generatedFiles.length < 2) return;
  const stem = (elements.outputName.value.trim() || activeModule().script.replace(/\.py$/i, "outputs")).replace(/\.[^.]+$/, "").replace(/[^a-z0-9._-]+/gi, "-");
  const link = document.createElement("a");
  link.href = URL.createObjectURL(makeZip(generatedFiles));
  link.download = `${stem || "phyloprep"}-outputs.zip`;
  link.click();
  window.setTimeout(() => URL.revokeObjectURL(link.href), 1000);
  setLog(`Downloaded ${generatedFiles.length} output files as a ZIP archive.`, "success");
}


function loadExternalScript(src) {
  return new Promise((resolve, reject) => {
    const existing = [...document.scripts].find((script) => script.src === src);
    if (existing?.dataset.loaded === "true") {
      resolve();
      return;
    }
    const script = existing || document.createElement("script");
    const onLoad = () => {
      script.dataset.loaded = "true";
      resolve();
    };
    const onError = () => reject(new Error(`Could not load ${src}.`));
    script.addEventListener("load", onLoad, { once: true });
    script.addEventListener("error", onError, { once: true });
    if (!existing) {
      script.src = src;
      script.async = true;
      document.head.append(script);
    }
  });
}

async function ensureMuscleCli() {
  if (!muscleCliPromise) {
    muscleCliPromise = (async () => {
      await loadExternalScript("https://biowasm.com/cdn/v3/aioli.js");
      if (typeof globalThis.Aioli !== "function") throw new Error("Biowasm Aioli did not initialize.");
      return await new globalThis.Aioli(["muscle/5.1.0"], { printInterleaved: false });
    })();
  }
  return muscleCliPromise;
}

async function runMuscleAlignment(sequenceText) {
  const cli = await ensureMuscleCli();
  const token = `phyloprep-${Date.now()}-${++muscleRunCounter}`;
  const mounted = await cli.mount([{ name: `${token}.fasta`, data: sequenceText }]);
  const inputPath = mounted?.[0];
  if (!inputPath) throw new Error("Biowasm could not mount the MUSCLE input.");
  const outputPath = `/shared/data/${token}.aligned.fasta`;
  let execution;
  let url;
  try {
    execution = await cli.exec("muscle", ["-align", inputPath, "-output", outputPath, "-threads", "1"]);
    url = await cli.download(outputPath);
    const alignment = await fetch(url).then((response) => response.text());
    if (!alignment.trimStart().startsWith(">")) {
      const detail = execution?.stderr || execution?.stdout || "MUSCLE did not produce a FASTA alignment.";
      throw new Error(String(detail).trim());
    }
    return alignment;
  } finally {
    if (url) URL.revokeObjectURL(url);
    try { await cli.fs.unlink(outputPath); } catch {}
    try { await cli.fs.unlink(inputPath); } catch {}
  }
}

async function runMsapBrowserWorkflow(request) {
  setRunStage(2, "Preparing input with original MSAP.py");
  elements.runtimeStatus.textContent = "Preparing MSAP.py workflow";
  const probe = await runInWorker({ type: "msap-probe", ...request });
  if (!probe.alignmentInput) throw new Error("MSAP.py did not provide an alignment input for MUSCLE.");

  elements.cancelRun.hidden = true;
  setRunStage(2, "Running Biowasm MUSCLE 5.1.0");
  elements.runtimeStatus.textContent = "Biowasm MUSCLE running";
  const alignmentContent = await runMuscleAlignment(probe.alignmentInput);

  setRunStage(2, "Finishing original MSAP.py workflow");
  elements.runtimeStatus.textContent = "Finishing MSAP.py workflow";
  elements.cancelRun.hidden = false;
  return runInWorker({ type: "msap-complete", ...request, alignmentContent });
}

async function runBrowserModule() {
  const module = activeModule();
  const files = browserInputFiles();
  const missingInput = inputControls().find((control) => !filesForInput(control).length);
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
    const filePayload = await Promise.all(files.map(async (file) => ({ name: file.name, content: await file.text() })));
    const optionValues = valuesForCommand();
    setRunStage(stage, "Loading local workflow");
    elements.cancelRun.hidden = false;
    const request = {
      script: module.script,
      files: filePayload,
      options: optionValues,
      outputName: elements.outputName.value.trim() || module.defaultOutput,
      usesPathList: module.inputFields.some((field) => field.pathList),
    };
    const result = module.fullWorkflow
      ? await runMsapBrowserWorkflow(request)
      : await runInWorker(request);
    stage = 3;
    setRunStage(stage, "Saving output");
    generatedFiles = [];
    for (const output of result.files) {
      const generated = { name: output.name, content: output.content };
      generated.destination = await saveGeneratedFile(generated);
      generatedFiles.push(generated);
    }
    renderGeneratedFiles();
    setRunStage(4, "Completed");
    setLog([result.summary || `Created ${generatedFiles.length} output file(s).`, result.provenance].filter(Boolean).join(" "), "success");
  } catch (error) {
    const cancelled = error.code === "CANCELLED";
    setRunStage(stage, cancelled ? "Run cancelled" : "Workflow needs attention", !cancelled);
    setLog(cancelled ? "Run cancelled. You can change inputs or parameters and try again." : (error.message || String(error)), cancelled ? "" : "error");
  } finally {
    setRunControls(false);
  }
}

elements.chooseFolder.addEventListener("click", chooseOutputFolder);
elements.outputName.addEventListener("input", updateCommand);
elements.run.addEventListener("click", runBrowserModule);
elements.cancelRun.addEventListener("click", cancelActiveRun);
elements.copyCommand.addEventListener("click", copyEquivalentCommand);
elements.downloadSelected.addEventListener("click", downloadSelectedOutput);
elements.downloadAll.addEventListener("click", downloadAllOutputs);

if ("serviceWorker" in navigator) {
  window.addEventListener("load", () => {
    navigator.serviceWorker.register("./sw.js").catch(() => {
      // The site remains fully usable when a browser or local preview blocks service workers.
    });
  });
}

selectModule(activeModuleIndex);
