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

from Bio.Data import CodonTable
from sequence_audit import normalize_workflow, merge_reports, write_report
from MSAP_batch import build_parser as batch_parser
from msap_io import atomic_output, expand_input_paths, read_fasta_alignment

SCRIPTS = Path(__file__).resolve().parent
QC_VALUES = {
    "min_length": (int, 1),
    "max_n_ratio": (float, 1.0),
    "max_x_ratio": (float, 1.0),
    "max_gap_ratio": (float, 1.0),
    "min_parsimony_informative_sites": (int, 0),
}
RENAME_STATE_VERSION = 1


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
Raw and trimmed alignments are QC-checked and concatenated independently;
--notrim processes raw alignments only. Codon mode also exports protein
supermatrices. Codon matrices additionally undergo position splitting and
four-fold extraction. FASTA, PHYLIP and NEXUS outputs are written under
matrices/<seqtype>/<raw|trimmed>/. QC thresholds apply independently to each
group; length and site counts use nucleotides for codon and amino acids for prot.
Pseudogene mode selects MACSE and writes matrices/pseudogene/<codon|prot>/<raw|trimmed>/.
It preserves native alignments and change reports, and automatically generates
codon-position and four-fold matrices from inferred coding positions. These do
not establish functional coding or synonymous-site neutrality. Empty four-fold
results are skipped with a message. See reports/ for audit TSVs.
No passing genes or no four-fold sites in normal codon mode is an error; final
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
    for name in ("allow-terminal-stop", "allow-internal-stop"):
        qc.add_argument("--" + name, action="store_true", help=name.replace("-", " ") + ".")
    return parser


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


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_signature(path):
    path = Path(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": file_digest(path)}


def read_json(path):
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def rename_state_path(directory, source):
    key = hashlib.sha256(str(Path(source).resolve()).encode()).hexdigest()[:16]
    return Path(directory) / f".phyloprep-rename-{key}.json"


def valid_rename_state(state, source, mapping, target, map_report):
    try:
        if Path(target).is_symlink() or Path(map_report).is_symlink():
            return False
        expected = {
            "version": RENAME_STATE_VERSION,
            "status": "complete",
            "source": file_signature(source),
            "mapping": file_signature(mapping),
            "rename_taxa_sha256": file_digest(SCRIPTS / "rename_taxa.py"),
            "outputs": {"alignment": file_signature(target), "map_report": file_signature(map_report)},
        }
    except OSError:
        return False
    return state == expected


def write_rename_state(path, source_signature, mapping_signature, script_digest, target, map_report):
    state = {
        "version": RENAME_STATE_VERSION,
        "status": "complete",
        "source": source_signature,
        "mapping": mapping_signature,
        "rename_taxa_sha256": script_digest,
        "outputs": {"alignment": file_signature(target), "map_report": file_signature(map_report)},
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


def publish_matrices(stage_root, root, variants):
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
        old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
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
    report_dir = root / "reports" / "renamed_ids" / seqtype / variant
    renamed_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    renamed = []
    for source in paths:
        target = renamed_dir / source.name
        map_report = report_dir / f"{source.name}.tsv"
        state_path = rename_state_path(renamed_dir, source)
        if not args.no_resume and valid_rename_state(
                read_json(state_path), source, args.mapping, target, map_report):
            print(f"[phyloprep] rename {source.name}: skipped (already complete)", flush=True)
            renamed.append(target)
            continue
        source_signature = file_signature(source)
        mapping_signature = file_signature(args.mapping)
        script_digest = file_digest(SCRIPTS / "rename_taxa.py")
        try:
            run("rename_taxa.py", "-i", source, "-o", target, "-m", args.mapping,
                "--map-output", map_report, "--allow-duplicate-ids")
        except subprocess.CalledProcessError as error:
            raise ValueError(
                f"ID mapping cannot rename {source} to unique taxon IDs. "
                "Each QC-passing alignment must contain at most one sequence for each mapped taxon. "
                "Multi-copy orthogroups should retain gene IDs for gene-tree methods such as ASTRAL-Pro3."
            ) from error
        if (file_signature(source) != source_signature or
                file_signature(args.mapping) != mapping_signature or
                file_digest(SCRIPTS / "rename_taxa.py") != script_digest):
            raise ValueError(
                f"Rename inputs changed while processing {source}; outputs were not checkpointed. Rerun the pipeline."
            )
        write_rename_state(state_path, source_signature, mapping_signature, script_digest, target, map_report)
        renamed.append(target)
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


def prepare_report_inputs(args, inputs, root):
    """Create mapped source copies used solely as audit-report provenance.

    Alignment still receives original IDs, which must remain unique even when
    several genes map to one taxon.  Sequence-change events, however, are
    meaningful to users under the supplied taxon names, including when a
    malformed input prevents an alignment from being produced.
    """
    if not args.mapping:
        return {}
    renamed_dir = root / "renamed"
    report_dir = root / "reports" / "renamed_ids" / "inputs"
    renamed_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    prepared = {}
    for source in inputs:
        target = renamed_dir / source.name
        map_report = report_dir / f"{source.name}.tsv"
        run("rename_taxa.py", "-i", source, "-o", target, "-m", args.mapping,
            "--map-output", map_report, "--allow-duplicate-ids")
        with map_report.open(newline="") as handle:
            mapping = {row["original_id"]: row["renamed_id"]
                       for row in csv.DictReader(handle, delimiter="\t")}
        prepared[str(source.resolve())] = (target.resolve(), mapping)
    return prepared


def relabel_change_report(path, prepared_inputs):
    """Apply taxon mapping to batch audit rows and point at mapped source copies."""
    if not prepared_inputs or not path.is_file():
        return
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        protein = "source_aa_positions" in (reader.fieldnames or [])
    changed = False
    global_mapping = {}
    for _target, mapping in prepared_inputs.values():
        global_mapping.update(mapping)
    for row in rows:
        prepared = prepared_inputs.get(str(Path(row["input_file"]).resolve()))
        normalized_id = row["taxa_id"].replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")
        renamed_id = global_mapping.get(normalized_id)
        if renamed_id is not None:
            row["taxa_id"] = renamed_id
            changed = True
        if prepared is not None:
            target, _mapping = prepared
            row["input_file"] = str(target)
            changed = True
    if changed:
        write_report(path, rows, protein=protein)


def prepare(args, inputs, root):
    recover_matrix_publication(root)
    if args.seqtype == 'pseudogene':
        print('[phyloprep] Pseudogene mode: frameshift/stop codons are masked for analysis. '
              'Codon-position and four-fold matrices are generated automatically from inferred coding positions; '
              'they do not establish that sequences are functional proteins or that sites are neutral synonymous sites. '
              'Positions containing masked or ambiguous codons are excluded by the four-fold extractor. '
              'If no usable four-fold sites remain, that matrix is skipped while other outputs are retained.', flush=True)
    prepared_inputs = prepare_report_inputs(args, inputs, root)
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
            # An early dependency/argument failure may leave an older batch report.
            updated = path.exists() and path.stat().st_mtime_ns != previous_stats[name]
            merge_reports([path] if updated else [], root / 'reports' / f'{name}.tsv')
            relabel_change_report(root / 'reports' / f'{name}.tsv', prepared_inputs)

    seqtypes = ("codon", "prot") if args.seqtype in {"codon", "pseudogene"} else (args.seqtype,)
    variants = ("raw",) if args.notrim else ("raw", "trimmed")
    # Finish every group before publishing any final matrices.
    with tempfile.TemporaryDirectory(prefix=".matrices-", dir=root) as temporary:
        stage_root = Path(temporary)
        for seqtype in seqtypes:
            for variant in variants:
                stage = stage_root / ("pseudogene" if args.seqtype == "pseudogene" else "") / seqtype / variant
                stage.mkdir(parents=True)
                prepare_group(args, inputs, root, stage, seqtype, variant)
        publish_matrices(stage_root, root, variants)

    print(f"[phyloprep] Complete: {root / 'matrices'}", flush=True)


def prepare_group(args, inputs, root, stage, seqtype, variant):
    """QC and concatenate one sequence type and trimming state independently."""
    print(f"[phyloprep] Processing {seqtype}/{variant}", flush=True)
    # Select only this invocation's inputs, never historical output globs.
    suffix = f".{args.align_software}.{seqtype}"
    suffix += ".aln" if variant == "raw" else ".trimal.aln"
    selected = [root / "alignments" / (source.stem + suffix) for source in inputs]
    selected = [path for path in selected if path.is_file()]
    if not selected:
        raise ValueError(f"No completed alignments are available for {seqtype}/{variant}")
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
        qc_args.extend(["--" + name.replace("_", "-"), getattr(args, name)])
    for name in ("allow_terminal_stop", "allow_internal_stop"):
        if getattr(args, name):
            qc_args.append("--" + name.replace("_", "-"))
    run("alignment_qc.py", *qc_args)
    passed = Path(str(prefix) + ".pass.tsv")
    if not passed.read_text().strip():
        if args.mapping:
            raise ValueError(
                f"No alignments passed QC for {seqtype}/{variant}; see {prefix}.details.tsv. "
                "Mapped IDs must be unique within every alignment. Multi-copy orthogroups with several genes "
                "mapped to one taxon require gene-tree methods such as ASTRAL-Pro3, not a species supermatrix."
            )
        raise ValueError(f"No alignments passed QC for {seqtype}/{variant}; see {prefix}.details.tsv")

    supermatrix = stage / "supermatrix.fasta"
    try:
        run("get_supergenes.py", "-i", passed, "-p", stage / "supermatrix",
            "--missing-taxa", args.missing_taxa)
    except subprocess.CalledProcessError as error:
        label = "renamed" if args.mapping else "original_ids"
        report = preserve_concat_report(stage, root, seqtype, variant, label)
        if not args.mapping:
            location = f" See {report}." if report else ""
            raise ValueError(
                f"No {seqtype}/{variant} genes can be concatenated with the original gene IDs.{location} "
                "If these IDs represent the same taxa across loci, provide phyloprep.py -m TAXA.tsv "
                "to rename all alignment matrices before QC and concatenation."
            ) from error
        location = f" See {report}." if report else ""
        raise ValueError(
            f"No {seqtype}/{variant} genes can be concatenated after ID mapping.{location} "
            "Check that every intended gene ID is mapped and that each alignment has one unique ID per taxon. "
            "Multi-copy orthogroups with several genes mapped to one taxon require gene-tree methods such as ASTRAL-Pro3, not a species supermatrix."
        ) from error
    matrices = [supermatrix]
    if seqtype == 'codon':
        run("split_codon_seqence_alignment.py", "-i", supermatrix, cwd=stage)
        matrices.extend(stage / f"codon{position}.supermatrix.fasta"
                        for position in ("1st", "2nd", "3rd"))
    if seqtype == 'codon':
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


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        normalize_workflow(args)
        if args.jobs < 1 or args.thread < 1:
            raise ValueError("--jobs and --thread must be positive")
        ratios = [args.G, args.N, args.X, args.max_n_ratio, args.max_x_ratio, args.max_gap_ratio]
        if any(not 0 <= value <= 1 for value in ratios):
            raise ValueError("Trimming and QC ratios must be between 0 and 1")
        if args.min_length < 0 or args.min_parsimony_informative_sites < 0:
            raise ValueError("QC length/site thresholds must be non-negative")
        if not args.notrim and args.trimal_args and args.trim_software != "trimal":
            raise ValueError("--trimal-args requires --trim-software trimal")
        if args.genetic_code not in CodonTable.unambiguous_dna_by_id:
            raise ValueError("Unknown genetic-code table")
        inputs = [Path(p).resolve() for p in expand_input_paths(args.input)]
        if not inputs or any(not p.is_file() for p in inputs):
            raise ValueError("Provide at least one existing input FASTA")
        if len({p.stem for p in inputs}) != len(inputs):
            raise ValueError("Input files must be unique and have unique basenames")
        if args.mapping:
            args.mapping = Path(args.mapping).resolve()
            if not args.mapping.is_file():
                raise ValueError(f"Mapping table does not exist: {args.mapping}")
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
