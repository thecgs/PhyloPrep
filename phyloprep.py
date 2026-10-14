#!/usr/bin/env python3
"""Prepare phylogenetic matrices through the existing MSAP command-line tools."""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile

from sequence_audit import normalize_workflow, merge_reports
from MSAP_batch import build_parser as batch_parser
from msap_io import (atomic_output, expand_input_paths, read_fasta_alignment,
                     file_digest, file_signature, read_json, NoGenesIncludedError)
from rename_taxa import rename_fasta, read_mapping_snapshot, normalized_taxon_id

SCRIPTS = Path(__file__).resolve().parent
QC_VALUES = {
    "min_length": (int, 1),
    "max_n_ratio": (float, 1.0),
    "max_x_ratio": (float, 1.0),
    "max_gap_ratio": (float, 1.0),
}
RENAME_STATE_VERSION = 2


def build_parser():
    parser = batch_parser()
    parser.description = "Align, quality-check, optionally rename, concatenate and export phylogenetic matrices."
    parser.set_defaults(output_dir="phyloprep-results")
    parser.epilog = """Examples:
  phyloprep.py -i genes.pathlist -m taxa.tsv -st codon -t 4 -j 2 -o results
  phyloprep.py -i gene1.fa gene2.fa --min-length 90 --max-gap-ratio 0.2
  phyloprep.py -i proteins.fa -st prot --notrim
  phyloprep.py -i genes.pathlist -st pseudogene --macse-jar macse_v2.07.jar -j 2 -o results

All MSAP_batch.py options are supported. --trimal-args must be last.
nucl/codon accept RNA: U is normalized to T and IUPAC ambiguity codes to N,
with audit rows; outputs use the canonical DNA alphabet. prot removes trailing
stops and masks internal stops, ambiguous and non-standard residues as X.
Use --protein-stop-symbol '.' for dot-encoded stops (default: '*'); edits are audited.
Raw and trimmed alignments are QC-checked once, then independently and jointly filtered;
--notrim processes raw alignments only. Codon mode also exports protein
supermatrices. Independent matrices use each representation's own QC/PIS passing loci;
additional paired matrices use the selected paired PIS policy and shared loci.
Paired outputs are under matrices/paired/ (matrices/pseudogene/paired/ for pseudogenes).
Codon matrices additionally undergo position splitting and
four-fold extraction. FASTA, PHYLIP and NEXUS outputs are written under
matrices/<seqtype>/<raw|trimmed>/. Non-PIS QC is evaluated separately for codon
and protein. Both sides must pass non-PIS QC only for the paired matrices;
paired PIS retention is then controlled by
--paired-pis-policy with independent codon/protein thresholds when requested.
Pseudogene mode selects MACSE and writes matrices/pseudogene/<codon|prot>/<raw|trimmed>/.
It preserves native alignments and change reports, and automatically generates
codon-position and four-fold matrices from inferred coding positions. These do
not establish functional coding or synonymous-site neutrality. Empty four-fold
results are skipped with a message. See reports/ for audit TSVs.
No loci in any independent/paired dataset, or no four-fold sites in a nonempty
normal codon dataset, is an error; final
matrices are published only after all requested conversions succeed.
"""
    for action in parser._actions:
        if isinstance(action, argparse._VersionAction):
            action.version = "v1.0.0"
        if action.dest == "output_dir":
            action.help = "Pipeline result directory (default: %(default)s)."
    pipeline = parser.add_argument_group("pipeline options")
    pipeline.add_argument("-m", "--mapping", "--map", metavar="TSV",
                          help="Rename all alignment IDs using this two-column mapping table before QC and concatenation.")
    pipeline.add_argument("--missing-taxa", choices=["skip-gene", "pad-gaps"],
                          default="skip-gene", help="Missing-taxon policy (default: skip-gene).")
    qc = parser.add_argument_group("alignment QC options (same defaults as alignment_qc.py)")
    for name, (kind, default) in QC_VALUES.items():
        unit = " (nt for nucl, aa for prot, complete codons for codon)" if name == "min_length" else ""
        qc.add_argument("--" + name.replace("_", "-"), type=kind, default=default,
                        help=f"{name.replace('_', ' ')}{unit} (default: %(default)s).")
    qc.add_argument("--allow-internal-stop", action="store_true",
                    help="Retain codon/protein locus pairs with internal stops detected directly or recorded as masked by MSAP preprocessing. The stop remains masked as NNN/X.")
    qc.add_argument("--paired-pis-policy", choices=["both", "either", "codon", "protein"],
                    default="both",
                    help="For additional paired matrices, combine the effective codon/protein PIS thresholds: both representations must pass (both; default), either may pass, or only codon/protein PIS controls pair retention.")
    qc.add_argument("--min-nucl-pis", type=int, default=None, metavar="N",
                    help="Minimum per-locus nucleotide-column PIS in nucl/codon/pseudogene modes, for independent and paired matrices (default: 0).")
    qc.add_argument("--min-protein-pis", type=int, default=None, metavar="N",
                    help="Minimum per-locus amino-acid PIS in prot/codon/pseudogene modes, for independent and paired matrices (default: 0).")
    return parser


def validate_pipeline_options(args):
    """Validate PhyloPrep-owned QC and paired-retention options."""
    qc_ratios = [args.max_n_ratio, args.max_x_ratio, args.max_gap_ratio]
    if any(not 0 <= value <= 1 for value in qc_ratios):
        raise ValueError("QC ratios must be between 0 and 1")
    if (args.min_length < 0 or
            (args.min_nucl_pis is not None and args.min_nucl_pis < 0) or
            (args.min_protein_pis is not None and args.min_protein_pis < 0)):
        raise ValueError("QC length/site thresholds must be non-negative")
    if args.paired_pis_policy not in {"both", "either", "codon", "protein"}:
        raise ValueError("Invalid paired PIS policy")
    if args.seqtype == "prot" and args.min_nucl_pis is not None:
        raise ValueError("--min-nucl-pis is not applicable to prot mode")
    if args.seqtype == "nucl" and args.min_protein_pis is not None:
        raise ValueError("--min-protein-pis is not applicable to nucl mode")
    if args.seqtype not in {"codon", "pseudogene"} and args.paired_pis_policy != "both":
        raise ValueError("--paired-pis-policy requires a paired codon/protein workflow")


def run(script, *arguments, cwd=None):
    command = [sys.executable, str(SCRIPTS / script), *map(str, arguments)]
    print("[phyloprep] " + shlex.join(command), flush=True)
    # Keep the batch alive long enough to clean up its independent job groups.
    process = None
    interrupted = False

    def cancel(_signum, _frame):
        nonlocal interrupted
        interrupted = True
        if process is not None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    previous = {sig: signal.signal(sig, cancel) for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        process = subprocess.Popen(command, cwd=cwd, start_new_session=True)
        if interrupted:
            cancel(None, None)
        while True:
            try:
                returncode = process.wait(timeout=0.2)
                break
            except subprocess.TimeoutExpired:
                if interrupted and script != "MSAP_batch.py":
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    returncode = process.wait()
                    break
        if interrupted:
            raise KeyboardInterrupt
        if returncode:
            raise subprocess.CalledProcessError(returncode, command)
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def write_paths(path, paths):
    with atomic_output(path) as handle:
        handle.writelines(f"{p}\n" for p in paths)


def rename_state_path(directory, source):
    key = hashlib.sha256(str(Path(source).resolve()).encode()).hexdigest()[:16]
    return Path(directory) / f".phyloprep-rename-{key}.json"


def valid_rename_state(state, source, mapping, target, mapping_signature=None, script_digest=None):
    try:
        if Path(target).is_symlink():
            return False
        expected = {
            "version": RENAME_STATE_VERSION,
            "status": "complete",
            "source": file_signature(source),
            "mapping": mapping_signature if mapping_signature is not None else file_signature(mapping),
            "rename_taxa_sha256": script_digest if script_digest is not None else file_digest(SCRIPTS / "rename_taxa.py"),
            "outputs": {"alignment": file_signature(target)},
        }
    except OSError:
        return False
    return state == expected


def write_rename_state(path, source_signature, mapping_signature, script_digest, target):
    state = {
        "version": RENAME_STATE_VERSION,
        "status": "complete",
        "source": source_signature,
        "mapping": mapping_signature,
        "rename_taxa_sha256": script_digest,
        "outputs": {"alignment": file_signature(target)},
    }
    with atomic_output(path) as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")


def batch_arguments(args, input_list, output):
    result = ["-i", str(input_list), "-o", str(output)]
    if args.seqtype == 'prot':
        result.extend(['--protein-stop-symbol', args.protein_stop_symbol])
    for flag, name in (("-t", "thread"), ("-j", "jobs"), ("-s", "align_software"),
                       ("-st", "seqtype"), ("-g", "genetic_code"), ("-ts", "trim_software"),
                       ("-G", "G"), ("-N", "N"), ("-X", "X")):
        result.extend([flag, str(getattr(args, name))])
    for name in ("notrim", "no_resume"):
        if getattr(args, name):
            result.append("--" + name.replace("_", "-"))
    if args.seqtype == "pseudogene":
        result.extend(["--macse-jar", args.macse_jar, "--macse-memory", args.macse_memory])
    if args.trimal_args:
        result.extend(["--trimal-args", *args.trimal_args])
    return result


def recover_matrix_publication(root):
    """Recover a directory switch interrupted before the new tree was installed."""
    destination = root / "matrices"
    backup = root / ".phyloprep-matrices-backup"
    if backup.exists():
        if destination.exists():
            shutil.rmtree(backup)
        else:
            os.replace(backup, destination)


def publish_matrices(stage_root, root, variants, block_signals=True):
    """Build a complete replacement tree, then switch it with rollback.

    Old groups and unmanaged files are preserved. The backup lives outside the
    temporary directory so a failed rollback or forced process exit cannot
    cause automatic temporary-directory cleanup to delete the old results.
    The caller holds the pipeline lock throughout publication and recovery.
    """
    recover_matrix_publication(root)
    destination = root / "matrices"
    backup = root / ".phyloprep-matrices-backup"
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=root) as temporary:
        candidate = Path(temporary) / "matrices"
        if destination.exists():
            shutil.copytree(destination, candidate)
        else:
            candidate.mkdir()
        # All validation and file writes happen privately, before changing any result.
        groups = sorted(path for path in stage_root.rglob('*')
                        if path.is_dir() and path.name in variants)
        for stage in groups:
            target = candidate / stage.relative_to(stage_root)
            target.mkdir(parents=True, exist_ok=True)
            manifest = target / '.phyloprep-products.json'
            previous = json.loads(manifest.read_text()) if manifest.exists() else []
            if not isinstance(previous, list) or any(
                    not isinstance(name, str) or not name or name in {'.', '..'} or
                    Path(name).name != name for name in previous):
                raise ValueError(f"Invalid matrix product manifest: {manifest}")
            products = list(stage.iterdir())
            names = [product.name for product in products]
            for product in products:
                shutil.copy2(product, target / product.name)
            for name in previous:
                if name not in names:
                    (target / name).unlink(missing_ok=True)
            with atomic_output(manifest) as handle:
                json.dump(names, handle)

        # Defer Ctrl-C/SIGTERM until the tree is either committed or restored.
        old_mask = (signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
                    if block_signals else None)
        moved_old = False
        try:
            if destination.exists():
                os.replace(destination, backup)
                moved_old = True
            try:
                os.replace(candidate, destination)
            except BaseException:
                if moved_old:
                    os.replace(backup, destination)
                raise
        finally:
            if block_signals:
                signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        if moved_old:
            # Cleanup failure must not turn a successful commit into a failed run.
            try:
                shutil.rmtree(backup)
            except OSError as error:
                print(f"[phyloprep] Matrices committed; retained backup {backup}: {error}",
                      file=sys.stderr)


def rename_alignment_matrices(args, paths, root, seqtype, variant):
    """Rename all aligned matrices after MSA has used unique gene IDs."""
    renamed_dir = root / "renamed" / seqtype / variant
    renamed_dir.mkdir(parents=True, exist_ok=True)
    renamed = []
    mapping_signature = getattr(args, "mapping_signature", None)
    if mapping_signature is None:
        args.taxon_mapping, mapping_signature = read_mapping_snapshot(args.mapping)
        args.mapping_signature = mapping_signature
    if file_signature(args.mapping) != mapping_signature:
        raise ValueError("Mapping changed after it was read; rerun before renaming alignments")
    script_digest = file_digest(SCRIPTS / "rename_taxa.py")
    for source in paths:
        target = renamed_dir / source.name
        state_path = rename_state_path(renamed_dir, source)
        if not args.no_resume and valid_rename_state(
                read_json(state_path), source, args.mapping, target, mapping_signature, script_digest):
            print(f"[phyloprep] rename {source.name}: skipped (already complete)", flush=True)
            renamed.append(target)
            continue
        source_signature = file_signature(source)
        try:
            rename_fasta(source, target, allow_duplicate_ids=True, mapping=args.taxon_mapping)
        except ValueError as error:
            raise ValueError(
                f"ID mapping cannot rename {source} to unique taxon IDs. "
                "Each QC-passing alignment must contain at most one sequence for each mapped taxon. "
                "Multi-copy orthogroups should retain gene IDs for gene-tree methods such as ASTRAL-Pro3."
            ) from error
        if file_signature(source) != source_signature:
            raise ValueError(
                f"Rename inputs changed while processing {source}; outputs were not checkpointed. Rerun the pipeline."
            )
        write_rename_state(state_path, source_signature, mapping_signature, script_digest, target)
        renamed.append(target)
    if (file_signature(args.mapping) != mapping_signature or
            file_digest(SCRIPTS / "rename_taxa.py") != script_digest):
        raise ValueError("Rename mapping or implementation changed during processing; rerun the pipeline")
    return renamed


def preserve_concat_report(stage, root, seqtype, variant, label):
    """Keep the diagnostic report when private matrix construction fails."""
    report = stage / "supermatrix.report.tsv"
    if not report.exists():
        return None
    destination = root / "reports" / "concatenation" / seqtype / variant / f"{label}.report.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(report, destination)
    return destination


def annotate_change_report(path, mapping):
    """Append mapped taxon IDs without changing primary audit provenance."""
    if not path.is_file():
        return
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    fields = list(reader.fieldnames or [])
    if "mapped_taxa_id" not in fields:
        fields.append("mapped_taxa_id")
    for row in rows:
        normalized = normalized_taxon_id(row["taxa_id"])
        row["mapped_taxa_id"] = mapping.get(normalized, normalized)
    with atomic_output(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows(rows)


def prepare(args, inputs, root):
    recover_matrix_publication(root)
    if args.seqtype == 'pseudogene':
        print('[phyloprep] Pseudogene mode: frameshift/stop codons are masked for analysis. '
              'Codon-position and four-fold matrices are generated automatically from inferred coding positions; '
              'they do not establish that sequences are functional proteins or that sites are neutral synonymous sites. '
              'Positions containing masked or ambiguous codons are excluded by the four-fold extractor. '
              'Audited stop/frameshift events are reported but do not reject pseudogene locus pairs. '
              'If no usable four-fold sites remain, that matrix is skipped while other outputs are retained.', flush=True)
    taxon_mapping = args.taxon_mapping if args.mapping else None
    input_list = root / "inputs.pathlist"
    write_paths(input_list, inputs)
    alignments = root / "alignments"
    report_files = {name: alignments / 'reports' / f'{name}.tsv'
                    for name in ('sequence_changes',)}
    previous_stats = {name: path.stat().st_mtime_ns if path.exists() else None
                      for name, path in report_files.items()}
    try:
        run("MSAP_batch.py", *batch_arguments(args, input_list, alignments))
    finally:
        for name, path in report_files.items():
            updated = path.exists() and path.stat().st_mtime_ns != previous_stats[name]
            merge_reports([path] if updated else [], root / 'reports' / f'{name}.tsv')
            if taxon_mapping is not None:
                annotate_change_report(root / 'reports' / f'{name}.tsv', taxon_mapping)

    variants = ("raw",) if args.notrim else ("raw", "trimmed")
    with tempfile.TemporaryDirectory(prefix=".matrices-", dir=root) as temporary:
        stage_root = Path(temporary)
        for variant in variants:
            prepare_variant(args, root, stage_root, variant)
        publish_matrices(stage_root, root, variants)

    print(f"[phyloprep] Complete: {root / 'matrices'}", flush=True)


def _group_stage(stage_root, workflow, seqtype, variant, paired=False):
    base = Path(stage_root)
    if workflow == "pseudogene":
        base = base / "pseudogene"
    if paired:
        base = base / "paired"
    return base / seqtype / variant


def _alignment_locus_key(path, args, seqtype, variant):
    name = Path(path).name
    suffix = f".{args.align_software}.{seqtype}"
    suffix += ".aln" if variant == "raw" else ".trimmed.aln"
    return name[:-len(suffix)] if name.endswith(suffix) else Path(name).stem


def run_group_qc(args, root, stage, seqtype, variant):
    """Run one representation-specific QC pass and return its pairing context."""
    print(f"[phyloprep] QC {seqtype}/{variant}", flush=True)
    list_suffix = f"{seqtype}.aln" if variant == "raw" else f"{seqtype}.trimmed.aln"
    batch_list = root / "alignments" / "path-lists" / f"all.{list_suffix}.pathlist"
    if not batch_list.is_file():
        raise ValueError(
            f"Validated MSAP batch path list is missing for {seqtype}/{variant}: {batch_list}"
        )
    selected = [
        Path(line.strip()).resolve()
        for line in batch_list.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith(("#", ";"))
    ]
    if not selected:
        raise ValueError(f"No validated completed alignments are available for {seqtype}/{variant}")

    qc_dir = root / "qc" / ("pseudogene" if args.seqtype == "pseudogene" else "") / seqtype / variant
    qc_paths = selected
    if args.mapping:
        qc_paths = rename_alignment_matrices(args, selected, root, seqtype, variant)
    qc_dir.mkdir(parents=True, exist_ok=True)
    qc_inputs = qc_dir / "inputs.pathlist"
    write_paths(qc_inputs, qc_paths)
    prefix = qc_dir / "alignment"
    qc_args = ["-i", qc_inputs, "-st", seqtype, "-g", args.genetic_code, "-p", prefix]
    for name in QC_VALUES:
        value = getattr(args, name)
        qc_args.extend(["--" + name.replace("_", "-"), value])
    pis_flag = "--min-protein-pis" if seqtype == "prot" else "--min-nucl-pis"
    threshold = args.min_protein_pis if seqtype == "prot" else args.min_nucl_pis
    # Paired workflows combine PIS decisions after computing non-PIS QC once.
    qc_args.extend([pis_flag, 0 if args.seqtype in {"codon", "pseudogene"} else (threshold or 0)])

    changes = root / "reports" / "sequence_changes.tsv"
    if seqtype in {"codon", "prot"} and changes.is_file():
        qc_args.extend(["--sequence-changes", changes])

    allow_internal = args.allow_internal_stop or args.seqtype == "pseudogene"
    if args.seqtype == "pseudogene":
        qc_args.append("--allow-terminal-stop")
    if allow_internal:
        qc_args.append("--allow-internal-stop")

    run("alignment_qc.py", *qc_args)
    return {
        "stage": Path(stage),
        "seqtype": seqtype,
        "variant": variant,
        "selected": selected,
        "qc_paths": qc_paths,
        "prefix": prefix,
        "passed": Path(str(prefix) + ".pass.tsv"),
        "failed": Path(str(prefix) + ".fail.tsv"),
    }


def _pass_map(context, args):
    result = {}
    text = context["passed"].read_text().strip()
    if not text:
        return result
    for line in text.splitlines():
        path = Path(line.strip())
        result[_alignment_locus_key(path, args, context["seqtype"], context["variant"])] = path
    return result


def _context_map(context, args):
    return {
        _alignment_locus_key(path, args, context["seqtype"], context["variant"]): Path(path)
        for path in context["qc_paths"]
    }


def _pis_map(context, args):
    """Reuse PIS already computed by alignment_qc.py instead of rereading alignments."""
    details = Path(str(context["prefix"]) + ".details.tsv")
    values = {}
    with details.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"File", "Parsimony_informative_sites"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"{details}: missing PIS columns from alignment_qc output")
        for row in reader:
            path = Path(row["File"])
            locus = _alignment_locus_key(
                path, args, context["seqtype"], context["variant"]
            )
            values[locus] = int(row["Parsimony_informative_sites"])
    expected = set(_context_map(context, args))
    if set(values) != expected:
        missing = sorted(expected - set(values))
        unexpected = sorted(set(values) - expected)
        raise ValueError(
            f"{details}: QC/PIS locus mismatch; missing={missing or 'none'}, "
            f"unexpected={unexpected or 'none'}"
        )
    return values

def _paired_pis_thresholds(args):
    return (args.min_nucl_pis or 0, args.min_protein_pis or 0)


def _paired_pis_decision(policy, codon_pis, protein_pis, codon_threshold, protein_threshold):
    codon_ok = codon_pis >= codon_threshold
    protein_ok = protein_pis >= protein_threshold
    if policy == "both":
        return codon_ok and protein_ok
    if policy == "either":
        return codon_ok or protein_ok
    if policy == "codon":
        return codon_ok
    if policy == "protein":
        return protein_ok
    raise ValueError(f"Unknown paired PIS policy: {policy}")


def synchronize_paired_qc(args, root, variant, codon_context, protein_context):
    """Keep one locus set while applying non-PIS QC and a user-selected paired PIS policy."""
    codon_all = _context_map(codon_context, args)
    protein_all = _context_map(protein_context, args)
    codon_pass = _pass_map(codon_context, args)
    protein_pass = _pass_map(protein_context, args)
    codon_pis = _pis_map(codon_context, args)
    protein_pis = _pis_map(protein_context, args)
    codon_threshold, protein_threshold = _paired_pis_thresholds(args)
    policy = args.paired_pis_policy

    non_pis_common = set(codon_pass) & set(protein_pass)
    common = sorted(
        locus for locus in non_pis_common
        if _paired_pis_decision(
            policy, codon_pis[locus], protein_pis[locus],
            codon_threshold, protein_threshold
        )
    )

    paired_dir = root / "qc" / ("pseudogene" if args.seqtype == "pseudogene" else "") / "paired" / variant
    paired_dir.mkdir(parents=True, exist_ok=True)
    report = paired_dir / "loci.details.tsv"
    with atomic_output(report) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "Locus", "Codon_file", "Protein_file",
            "Codon_nonPIS_QC", "Protein_nonPIS_QC",
            "Codon_PIS", "Protein_PIS", "Minimum_codon_PIS", "Minimum_protein_PIS",
            "PIS_policy", "PIS_status",
            "Paired_status", "Reasons",
        ])
        for locus in sorted(set(codon_all) | set(protein_all)):
            cp = locus in codon_pass
            pp = locus in protein_pass
            cp_pis = codon_pis.get(locus, "NA")
            pp_pis = protein_pis.get(locus, "NA")
            pis_ok = (
                cp_pis != "NA" and pp_pis != "NA" and
                _paired_pis_decision(
                    policy, cp_pis, pp_pis, codon_threshold, protein_threshold
                )
            )
            reasons = []
            if not cp:
                reasons.append("codon_qc_failed")
            if not pp:
                reasons.append("protein_qc_failed")
            if not pis_ok:
                reasons.append("pis_policy_failed")
            paired_ok = cp and pp and pis_ok
            writer.writerow([
                locus,
                str(codon_all.get(locus, "NA")),
                str(protein_all.get(locus, "NA")),
                "pass" if cp else "fail",
                "pass" if pp else "fail",
                cp_pis,
                pp_pis,
                codon_threshold,
                protein_threshold,
                policy,
                "pass" if pis_ok else "fail",
                "pass" if paired_ok else "fail",
                ";".join(reasons) or "NA",
            ])

    for context, passing, pis, threshold in (
            (codon_context, codon_pass, codon_pis, codon_threshold),
            (protein_context, protein_pass, protein_pis, protein_threshold)):
        independent_path = context["prefix"].parent / "independent.retained.tsv"
        independent = sorted(locus for locus in passing if pis[locus] >= threshold)
        write_paths(independent_path, [passing[locus].resolve() for locus in independent])
        context["independent_retained"] = independent_path

    for context, passing in ((codon_context, codon_pass), (protein_context, protein_pass)):
        retained_path = paired_dir / f"{context['seqtype']}.retained.tsv"
        with atomic_output(retained_path) as handle:
            for locus in common:
                handle.write(str(passing[locus].resolve()) + "\n")
        context["paired_retained"] = retained_path

    print(
        f"[phyloprep] Paired QC {variant}: retained {len(common)} locus/loci "
        f"(PIS policy={policy}, codon minimum={codon_threshold}, "
        f"protein minimum={protein_threshold})",
        flush=True,
    )
    return common


def _require_passing_group(args, context):
    if context["passed"].read_text().strip():
        return
    seqtype, variant = context["seqtype"], context["variant"]
    prefix = context["prefix"]
    if args.mapping:
        raise ValueError(
            f"No alignments passed QC for {seqtype}/{variant}; see {prefix}.details.tsv. "
            "Mapped IDs must be unique within every alignment. Multi-copy orthogroups with several genes "
            "mapped to one taxon require gene-tree methods such as ASTRAL-Pro3, not a species supermatrix."
        )
    raise ValueError(f"No alignments passed QC for {seqtype}/{variant}; see {prefix}.details.tsv")


def build_group_matrices(args, root, context):
    """Concatenate one already-QC-filtered representation and export derived matrices."""
    stage = context["stage"]
    seqtype = context["seqtype"]
    variant = context["variant"]
    passed = context.get("retained", context["passed"])
    stage.mkdir(parents=True, exist_ok=True)

    supermatrix = stage / "supermatrix.fasta"
    try:
        run("get_supergenes.py", "-i", passed, "-p", stage / "supermatrix",
            "--missing-taxa", args.missing_taxa, "--no-phylip")
    except subprocess.CalledProcessError as error:
        if error.returncode != 3:
            raise ValueError(
                f"Concatenation failed for {seqtype}/{variant} (exit {error.returncode}); "
                "see the original get_supergenes error above."
            ) from error
        label = "renamed" if args.mapping else "original_ids"
        if context.get("paired"):
            label = "paired_" + label
        report = preserve_concat_report(stage, root, seqtype, variant, label)
        if not args.mapping:
            location = f" See {report}." if report else ""
            raise NoGenesIncludedError(
                f"No {seqtype}/{variant} genes can be concatenated with the original gene IDs.{location} "
                "If these IDs represent the same taxa across loci, provide phyloprep.py -m TAXA.tsv "
                "to rename all alignment matrices before QC and concatenation."
            ) from error
        location = f" See {report}." if report else ""
        raise NoGenesIncludedError(
            f"No {seqtype}/{variant} genes can be concatenated after ID mapping.{location} "
            "Check that every intended gene ID is mapped and that each alignment has one unique ID per taxon. "
            "Multi-copy orthogroups with several genes mapped to one taxon require gene-tree methods such as ASTRAL-Pro3, not a species supermatrix."
        ) from error

    matrices = [supermatrix]
    if seqtype == 'codon':
        run("split_codon_seqence_alignment.py", "-i", supermatrix, cwd=stage)
        matrices.extend(stage / f"codon{position}.supermatrix.fasta"
                        for position in ("1st", "2nd", "3rd"))
        fourfold = stage / "fourfold.fasta"
        run("extract_4-fold_degenerated_sites.py", "-i", supermatrix,
            "-o", fourfold, "-g", args.genetic_code)
        try:
            read_fasta_alignment(fourfold)
        except ValueError as error:
            if args.seqtype != 'pseudogene':
                raise ValueError(f"No usable four-fold matrix for {seqtype}/{variant} "
                                 "(possibly zero qualifying sites); final matrices were not published") from error
            print(f'[phyloprep] {seqtype}/{variant}: no usable four-fold sites; skipping FASTA/PHY/NEX for this matrix.', flush=True)
            fourfold.unlink()
            (stage / 'fourfold.skipped.txt').write_text('No usable four-fold sites; no matrix exported.\n')
        else:
            matrices.append(fourfold)

    for matrix in matrices:
        run("fasta2phy.py", "-i", matrix, "-o", matrix.with_suffix(".phy"),
            "-st", "protein" if seqtype == "prot" else "DNA")
        run("fasta2nex.py", "-i", matrix, "-o", matrix.with_suffix(".nex"),
            "-st", "protein" if seqtype == "prot" else "DNA")


class NoAlignmentsPassedQCError(ValueError):
    """Expected QC rejection with diagnostic files available for inspection."""


def no_passing_paired_qc_error(args, variant, contexts):
    message = f"No codon/protein alignments passed independent or paired QC for {variant}."
    details = Path(str(contexts["codon"]["prefix"]) + ".details.tsv")
    with details.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    stopped = sum("internal_stop" in row.get("Issues", "").split(";") for row in rows)
    if stopped:
        message += (
            f" {stopped} of {len(rows)} codon loci were rejected for internal stops "
            f"(including stops masked during preprocessing). Current genetic-code table: -g {args.genetic_code}. "
            "Check that -g matches the source CDS; vertebrate mitochondrial CDS typically use table 2, "
            "rather than the standard table 1. Internal stops can also reflect sequence or annotation errors. "
            "Correct the genetic-code setting and rerun from the original CDS if it is wrong; "
            "do not use --allow-internal-stop as a substitute for correcting it."
        )
    message += " Inspect alignment.details.tsv and reports/sequence_changes.tsv for rejection reasons."
    return NoAlignmentsPassedQCError(message)


def prepare_variant(args, root, stage_root, variant):
    """Prepare one raw/trimmed variant, pairing codon and protein when applicable."""
    if args.seqtype in {"codon", "pseudogene"}:
        contexts = {}
        for seqtype in ("codon", "prot"):
            stage = _group_stage(stage_root, args.seqtype, seqtype, variant)
            contexts[seqtype] = run_group_qc(args, root, stage, seqtype, variant)
        paired_loci = synchronize_paired_qc(
            args, root, variant, contexts["codon"], contexts["prot"])
        if not paired_loci and not any(
                context["independent_retained"].read_text().strip()
                for context in contexts.values()):
            raise no_passing_paired_qc_error(args, variant, contexts)
        built = False
        for seqtype in ("codon", "prot"):
            context = contexts[seqtype]
            for paired in (False, True):
                selected = context["paired_retained" if paired else "independent_retained"]
                stage = _group_stage(stage_root, args.seqtype, seqtype, variant, paired=paired)
                if not selected.read_text().strip():
                    stage.mkdir(parents=True, exist_ok=True)
                    (stage / "skipped.txt").write_text(
                        "No loci passed " + ("paired" if paired else "independent") + " QC.\n")
                    continue
                independent_stage = context["stage"]
                if (paired and selected.read_text() == context["independent_retained"].read_text()
                        and (independent_stage / "supermatrix.fasta").is_file()):
                    shutil.copytree(independent_stage, stage)
                    built = True
                    continue
                try:
                    build_group_matrices(args, root, {**context, "stage": stage,
                                                     "retained": selected, "paired": paired})
                except NoGenesIncludedError as error:
                    # An empty missing-taxon selection does not invalidate other datasets.
                    (stage / "skipped.txt").write_text(str(error) + "\n")
                else:
                    built = True
        if not built:
            raise NoGenesIncludedError(f"No independent or paired matrices could be built for {variant}")
        return

    stage = _group_stage(stage_root, args.seqtype, args.seqtype, variant)
    context = run_group_qc(args, root, stage, args.seqtype, variant)
    _require_passing_group(args, context)
    build_group_matrices(args, root, context)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        normalize_workflow(args)
        validate_pipeline_options(args)
        inputs = [Path(p).resolve() for p in expand_input_paths(args.input)]
        if not inputs or any(not p.is_file() for p in inputs):
            raise ValueError("Provide at least one existing input FASTA")
        if len({p.stem for p in inputs}) != len(inputs):
            raise ValueError("Input files must be unique and have unique basenames")
        if args.mapping:
            args.mapping = Path(args.mapping).resolve()
            if not args.mapping.is_file():
                raise ValueError(f"Mapping table does not exist: {args.mapping}")
            args.taxon_mapping, args.mapping_signature = read_mapping_snapshot(args.mapping)
        else:
            args.taxon_mapping = None
        root = Path(args.output_dir).resolve()
        # Pipeline-owned subdirectories must not contain source inputs.
        for source in [*inputs, *([args.mapping] if args.mapping else [])]:
            if source == root or root in source.parents:
                raise ValueError("Input FASTAs and mapping must be outside the output directory")
        root.mkdir(parents=True, exist_ok=True)
        with (root / ".phyloprep.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError("Another phyloprep process is using this output directory") from None
            prepare(args, inputs, root)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"phyloprep: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("phyloprep: interrupted", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
