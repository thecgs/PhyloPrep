#!/usr/bin/env python3
"""Run MSAP.py on many FASTA files with resumable, interrupt-safe jobs.

Each input is executed in a private staging directory.  Its files become
visible in the result directory only after MSAP.py exits successfully.
"""

from __future__ import annotations

import argparse
import fcntl
import concurrent.futures
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import tempfile
from pathlib import Path
from typing import Any

from msap_io import atomic_output, read_fasta_alignment, expand_input_paths
from sequence_audit import add_macse_options, add_protein_options, normalize_workflow, merge_reports, relocate_change_report


STOP = False
STATE_VERSION = 2


def _stop_handler(signum: int, _frame: Any) -> None:
    global STOP
    STOP = True
    print(f"Received {signal.Signals(signum).name}; terminating active task process groups...", file=sys.stderr)


def _input_key(path: Path) -> str:
    digest = hashlib.sha256(str(path.resolve()).encode()).hexdigest()[:10]
    return f"{path.stem}_{digest}"


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": _digest(path)}


def _read_json(path: Path) -> dict[str, Any]:
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
        return state if isinstance(state, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    with atomic_output(path) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _configuration(args: argparse.Namespace) -> dict[str, Any]:
    from MSAP import check_dependencies
    import Bio
    softwares = ["java" if args.seqtype == "pseudogene" else args.align_software]
    if not args.notrim and args.trim_software == "trimal":
        softwares.append("trimal")
    binaries = check_dependencies(softwares)
    if args.seqtype == "pseudogene":
        binaries["macse_jar"] = args.macse_jar
    scripts = ["MSAP.py", "MSAP_batch.py", "AA2Codon.py", "trimAlnSeq.py", "msap_io.py", "sequence_audit.py"]
    options = ("thread", "align_software", "seqtype", "genetic_code", "notrim",
               "trim_software", "G", "N", "X", "trimal_args", "macse_jar", "macse_memory", "protein_stop_symbol")
    return {"options": {key: getattr(args, key) for key in options},
            "scripts": {name: _digest(Path(__file__).with_name(name)) for name in scripts},
            "tools": {name: _signature(Path(path)) for name, path in binaries.items()},
            "python": sys.version, "biopython": Bio.__version__}


def _claim_prefixes(inputs: list[Path], root: Path) -> None:
    """Called under the output-directory lock, before any jobs mutate files."""
    owner_file = root / ".msap-batch-owners.json"
    owners = _read_json(owner_file)
    if owner_file.exists():
        try:
            owners = json.loads(owner_file.read_text(encoding='utf-8'))
            if not isinstance(owners, dict) or any(not isinstance(v, str) for v in owners.values()):
                raise ValueError('Invalid ownership data')
        except (OSError, ValueError) as error:
            raise ValueError("Cannot read output-prefix ownership; use a new output directory") from error
    # Recover ownership from older checkpoints, including unfinished jobs.
    for state_file in root.glob(".msap-batch-state-*.json"):
        state = _read_json(state_file)
        source = state.get("input", {}).get("path")
        if not source:
            raise ValueError(f"Cannot determine output ownership from {state_file.name}")
        key = state_file.name[len(".msap-batch-state-"):-len(".json")]
        stem = key.rsplit("_", 1)[0]
        if stem in owners and owners[stem] != source:
            raise ValueError(f"Conflicting historical owners for prefix {stem}; use a new output directory")
        owners[stem] = source
    for path in inputs:
        source = str(path.resolve())
        owner = owners.get(path.stem)
        if owner is not None and owner != source:
            raise ValueError(f"Output prefix {path.stem} already belongs to {owner}; use another output directory")
        if owner is None and any(p.name.startswith(path.stem + ".") for p in root.iterdir()):
            raise ValueError(f"Unmanaged outputs already use prefix {path.stem}; use another output directory")
        owners[path.stem] = source
    _atomic_json(owner_file, owners)


def _build_command(args: argparse.Namespace, infile: Path) -> list[str]:
    command = [sys.executable, str(Path(__file__).resolve().with_name("MSAP.py")), "-i", str(infile.resolve()),
               "-t", str(args.thread), "-s", args.align_software, "-st", args.seqtype,
               "-g", str(args.genetic_code), "-ts", args.trim_software,
               "-G", str(args.G), "-N", str(args.N), "-X", str(args.X)]
    if args.seqtype == 'prot':
        command.extend(['--protein-stop-symbol', args.protein_stop_symbol])
    if args.seqtype == "pseudogene":
        command.extend(["--macse-jar", args.macse_jar, "--macse-memory", args.macse_memory])
    if args.notrim:
        command.append("--notrim")
    if args.trimal_args:
        command.extend(["--trimal-args", *args.trimal_args])
    return command


def _expected_alignment_outputs(args: argparse.Namespace, infile: Path, stage: Path) -> list[Path]:
    suffixes = ["prot", "codon"] if args.seqtype in {"codon", "pseudogene"} else [args.seqtype]
    outputs = [stage / f"{infile.stem}.{args.align_software}.{suffix}.aln" for suffix in suffixes]
    if not args.notrim:
        outputs.extend(stage / f"{infile.stem}.{args.align_software}.{suffix}.trimal.aln" for suffix in suffixes)
    return outputs


def _manifest_valid(manifest: dict[str, Any], root: Path, expected: list[Path]) -> bool:
    if not isinstance(manifest, dict) or not manifest:
        return False
    if not all(path.name in manifest for path in expected):
        return False
    for name, metadata in manifest.items():
        if Path(name).name != name or not isinstance(metadata, dict):
            return False
        path = root / name
        try:
            if (not path.is_file() or path.is_symlink() or
                    path.stat().st_size != metadata.get("size") or
                    _digest(path) != metadata.get("sha256")):
                return False
        except OSError:
            return False
    return True


def _write_path_lists(root: Path, args: argparse.Namespace, inputs: list[Path]) -> None:
    """Write absolute-path lists for each alignment type produced by the batch."""
    list_dir = root / "path-lists"
    list_dir.mkdir(parents=True, exist_ok=True)
    suffixes = ["prot.aln", "codon.aln"] if args.seqtype in {"codon", "pseudogene"} else [f"{args.seqtype}.aln"]
    if not args.notrim:
        suffixes += [suffix.replace(".aln", ".trimal.aln") for suffix in list(suffixes)]
    for kind in ('nucl', 'prot', 'codon'):
        for variant in ('aln', 'trimal.aln'):
            suffix = f'{kind}.{variant}'
            if suffix not in suffixes:
                (list_dir / f'all.{suffix}.pathlist').unlink(missing_ok=True)
    valid_outputs = set()
    for source in inputs:
        state = _read_json(root / f".msap-batch-state-{_input_key(source)}.json")
        expected = _expected_alignment_outputs(args, source, root)
        try:
            valid = (state.get("version") == STATE_VERSION and
                     state.get("status") == "complete" and
                     state.get("input") == _signature(source) and
                     state.get("configuration") == args.run_configuration and
                     _manifest_valid(state.get("outputs"), root, expected))
        except OSError:
            valid = False
        if valid:
            valid_outputs.update(expected)
    for suffix in suffixes:
        files = sorted(path for path in valid_outputs if path.name.endswith("." + suffix))
        list_path = list_dir / f"all.{suffix}.pathlist"
        with atomic_output(list_path) as handle:
            handle.write("".join(f"{path.resolve()}\n" for path in files))


def _terminate_group(process: subprocess.Popen) -> None:
    # The child has its own session; descendants share its process group.
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        process.poll()  # Reap the direct child while waiting for descendants.
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.05)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def _copy_replace(source: Path, destination: Path) -> None:
    """Keep source intact, and expose only complete destination files."""
    with tempfile.TemporaryDirectory(prefix='.msap-publish-', dir=destination.parent) as temporary:
        candidate = Path(temporary) / destination.name
        shutil.copy2(source, candidate)
        os.replace(candidate, destination)


def _recover_publication(backup: Path, root: Path) -> None:
    journal_path = backup / 'journal.json'
    if not journal_path.exists():
        # Preparation did not reach the first modification of public files.
        shutil.rmtree(backup)
        return
    journal = json.loads(journal_path.read_text(encoding='utf-8'))
    if journal.get('status') in {'rolled_back', 'committed'}:
        if journal['status'] == 'committed':
            stage = root / '.msap-batch-staging' / backup.name
            if stage.exists():
                shutil.rmtree(stage)
        shutil.rmtree(backup)
        return
    names = journal['names']
    if any(Path(name).name != name for name in names):
        raise ValueError('Invalid batch publication journal')
    state_file = root / f'.msap-batch-state-{backup.name}.json'
    state = _read_json(state_file)
    committed = (state.get('status') == 'complete' and
                 state.get('outputs') == journal['outputs'] and
                 _manifest_valid(journal['outputs'], root, []))
    if not committed:
        original_files = journal.get('original_files')
        if original_files is None:
            # Older journals cannot distinguish deleted backups from absent originals.
            raise ValueError(f'Legacy publication journal requires manual recovery: {backup}')
        if not isinstance(original_files, list) or not set(original_files) <= set(names):
            raise ValueError('Invalid original-file inventory in publication journal')
        for name in original_files:
            saved = backup / 'files' / name
            if not saved.is_file() or saved.is_symlink():
                raise ValueError(f'Missing rollback copy; public outputs left intact: {saved}')
        for name in names:
            saved = backup / 'files' / name
            if name in original_files:
                _copy_replace(saved, root / name)
            else:
                (root / name).unlink(missing_ok=True)
        _atomic_json(state_file, journal['running_state'])
        _atomic_json(journal_path, {**journal, 'status': 'rolled_back'})
    else:
        _atomic_json(journal_path, {**journal, 'status': 'committed'})
        stage = root / '.msap-batch-staging' / backup.name
        if stage.exists():
            shutil.rmtree(stage)
    shutil.rmtree(backup)


def _publish_task(stage, root, state_file, running, report_name):
    """Publish with durable rollback copies; original staging reports stay valid."""
    produced = list(stage.iterdir())
    if any(not path.is_file() or path.is_symlink() for path in produced):
        raise RuntimeError('Unexpected non-file output in task staging directory')
    backup = root / '.msap-batch-backups' / stage.name
    backup.mkdir(parents=True)
    try:
        with tempfile.TemporaryDirectory(prefix='.msap-publish-', dir=root) as temporary:
            candidate = Path(temporary)
            for path in produced:
                shutil.copy2(path, candidate / path.name)
            relocate_change_report(candidate / report_name, stage, root)
            manifest = {p.name: {'size': p.stat().st_size, 'sha256': _digest(p)} for p in candidate.iterdir()}
            names = sorted(set(manifest) | set(running.get('outputs', {})))
            saved = backup / 'files'
            saved.mkdir()
            original_files = []
            for name in names:
                if Path(name).name != name:
                    raise ValueError('Invalid previous output name')
                destination = root / name
                if destination.is_symlink() or (destination.exists() and not destination.is_file()):
                    raise ValueError(f'Cannot replace non-regular output: {destination}')
                if destination.exists():
                    shutil.copy2(destination, saved / name)
                    original_files.append(name)
            _atomic_json(backup / 'journal.json', {'names': names, 'outputs': manifest,
                         'running_state': running, 'original_files': original_files, 'status': 'pending'})
            for name in manifest:
                if STOP:
                    raise KeyboardInterrupt
                os.replace(candidate / name, root / name)
            for name in names:
                if name not in manifest:
                    (root / name).unlink(missing_ok=True)
            if STOP:
                raise KeyboardInterrupt
            _atomic_json(state_file, {**running, 'status': 'complete', 'outputs': manifest, 'finished': time.time()})
    except BaseException:
        _recover_publication(backup, root)
        raise
    # Recovery recognizes the committed checkpoint and only removes backups/staging.
    _recover_publication(backup, root)


def _run_one(args: argparse.Namespace, infile: Path, root: Path) -> tuple[str, str]:
    key = _input_key(infile)
    state_file = root / f".msap-batch-state-{key}.json"
    signature = _signature(infile)
    state = _read_json(state_file)
    expected = _expected_alignment_outputs(args, infile, root)
    if (not args.no_resume and state.get("version") == STATE_VERSION and
            state.get("status") == "complete" and state.get("input") == signature and
            state.get("configuration") == args.run_configuration and
            _manifest_valid(state.get("outputs"), root, expected)):
        return key, "skipped (already complete)"
    if STOP:
        return key, "cancelled before start"

    staging_root = root / ".msap-batch-staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    stage = staging_root / key
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    previous_outputs = state.get("outputs", {})
    running = {"version": STATE_VERSION, "status": "running", "input": signature,
               "configuration": args.run_configuration, "outputs": previous_outputs,
               "started": time.time()}
    # Invalidate the old completion marker before rerunning or promoting files.
    _atomic_json(state_file, running)
    command = _build_command(args, infile)
    process = subprocess.Popen(command, cwd=stage, start_new_session=True)
    try:
        while process.poll() is None:
            if STOP:
                raise KeyboardInterrupt
            time.sleep(0.1)
        if STOP:
            raise KeyboardInterrupt
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)
        for path in _expected_alignment_outputs(args, infile, stage):
            read_fasta_alignment(path, codon=path.name.endswith((".codon.aln", ".codon.trimal.aln")))
        if _signature(infile) != signature:
            raise RuntimeError("Input changed during alignment; results were not promoted")
    except BaseException:
        _terminate_group(process)
        raise

    _publish_task(stage, root, state_file, running, f'{infile.stem}.sequence_changes.tsv')
    return key, "completed"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run MSAP.py on multiple FASTA files concurrently with validated results and resume support.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""MSAP.py-compatible examples:
  MSAP_batch.py -i CDS.fasta -st codon -g 1
  MSAP_batch.py -i genes.pathlist -st pseudogene --macse-jar macse_v2.07.jar -j 2 -o results
  MSAP_batch.py -i protein1.fasta protein2.fasta -st prot -s muscle -t 4 -j 2
  MSAP_batch.py -i 16S.fasta rbcL.fasta -st nucl --notrim
  MSAP_batch.py -i gene1.fasta gene2.fasta -st codon --trim-software trimal \\
      --trimal-args -automated1

Batch-only options:
  -j/--jobs controls how many input files run concurrently.
  -t/--thread controls MAFFT/MUSCLE threads; MACSE, PRANK and ClustalW2 do not use it.
  nucl/codon accept RNA, normalize U to T and IUPAC ambiguity codes to N, and audit the changes.
  codon trims a trailing 1-2-base excess with a warning and terminal_partial_codon audit event.
  prot removes trailing stops and masks internal stops, ambiguous and non-standard residues as X,
  with recorded amino-acid positions.
  --protein-stop-symbol '.' selects dot-encoded stops (default: '*').
  Input accepts individual FASTA paths, multiple FASTA paths, a path-list file
  (one FASTA path per line, any extension), or any mixture of these forms.
  Path lists include only validated, completed inputs from the current run;
  failed or interrupted replacements are excluded even if old files remain.
  After completion, all path lists are written under the output directory's
  path-lists/ subdirectory (for example, path-lists/all.codon.aln.pathlist).
  Results are written directly into <output-dir> only after success; incomplete work remains in
  .msap-batch-staging and is restarted on the next run.""",
    )
    required = parser.add_argument_group("MSAP.py required arguments")
    workflow = parser.add_argument_group("MSAP.py workflow options")
    trim_seq = parser.add_argument_group("MSAP.py trimAlnSeq.py options")
    trim_al = parser.add_argument_group("MSAP.py trimAl options")
    batch = parser.add_argument_group("batch-only options")
    required.add_argument("-i", "--input", nargs="+", required=True, metavar="FASTA|LIST", help="FASTA files and/or path-list files (one FASTA path per line).")
    workflow.add_argument("-t", "--thread", type=int, default=os.cpu_count(), help="MAFFT/MUSCLE threads per input; unused by MACSE/PRANK/ClustalW2 (default: %(default)s).")
    workflow.add_argument("-s", "--align_software", "--align-software", default=None, choices=["mafft", "muscle", "prank", "clustalw2", "macse"], help="Alignment program (default: mafft; pseudogene requires macse).")
    workflow.add_argument("-n", "--notrim", action="store_true", help="Skip alignment trimming.")
    workflow.add_argument("-ts", "--trim-software", default="trimAlnSeq", choices=["trimAlnSeq", "trimal"], help="Trimming program (default: trimAlnSeq).")
    trim_seq.add_argument("-G", "--G", "--trimAlnSeq-G", dest="G", type=float, default=0.2, help="Maximum gap ratio (default: 0.2).")
    trim_seq.add_argument("-N", "--N", "--trimAlnSeq-N", dest="N", type=float, default=0.2, help="Maximum N ratio (default: 0.2).")
    trim_seq.add_argument("-X", "--X", "--trimAlnSeq-X", dest="X", type=float, default=0.2, help="Maximum X ratio (default: 0.2).")
    trim_al.add_argument("--trimal-args", nargs=argparse.REMAINDER, default=[], help="Arguments passed to trimAl; must be final. In codon/pseudogene mode, do not pass -in, -out, -fasta, -backtrans or -colnumbering.")
    workflow.add_argument("-st", "--seqtype", "--seq-type", default="codon", choices=["codon", "prot", "nucl", "pseudogene"], help="Input sequence type (default: codon).")
    workflow.add_argument("-g", "--genetic_code", "--genetic-code", type=int, default=1, help="NCBI genetic-code table (default: 1).")
    batch.add_argument("-o", "--output-dir", default="msap-results", metavar="DIR", help="Result directory (default: msap-results).")
    batch.add_argument("-j", "--jobs", type=int, default=1, help="Input files processed simultaneously (default: 1).")
    batch.add_argument("--no-resume", action="store_true", help="Ignore completed task records and rerun all inputs.")
    parser.add_argument("-v", "--version", action="version", version="MSAP_batch 1.0")
    add_macse_options(workflow)
    add_protein_options(workflow)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        normalize_workflow(args)
    except ValueError as error:
        parser.error(str(error))
    if args.jobs < 1 or args.thread < 1:
        parser.error("--jobs and --thread must be positive integers")
    if not args.notrim and args.trimal_args and args.trim_software != "trimal":
        parser.error("--trimal-args requires --trim-software trimal")
    if any(not 0 <= value <= 1 for value in (args.G, args.N, args.X)):
        parser.error("-G, -N, and -X must be between 0 and 1")
    try:
        inputs = [Path(item).resolve() for item in expand_input_paths(args.input)]
    except ValueError as error:
        parser.error(str(error))
    if not inputs:
        parser.error('No input FASTA files were supplied')
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        parser.error("input file(s) do not exist: " + ", ".join(missing))
    if len({str(path.resolve()) for path in inputs}) != len(inputs):
        parser.error("input files must be unique")
    stems = [path.stem for path in inputs]
    if len(set(stems)) != len(stems):
        parser.error("input files must have unique basenames when using a flat output directory")

    global STOP
    STOP = False
    signal.signal(signal.SIGINT, _stop_handler)
    signal.signal(signal.SIGTERM, _stop_handler)
    root = Path(args.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    args.run_configuration = _configuration(args)
    with (root / ".msap-batch.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("another MSAP batch is already using this output directory")
        try:
            for backup in sorted((root / '.msap-batch-backups').glob('*')):
                _recover_publication(backup, root)
            _claim_prefixes(inputs, root)
        except ValueError as error:
            parser.error(str(error))
        failures = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(_run_one, args, infile, root): infile for infile in inputs}
            try:
                for future in concurrent.futures.as_completed(futures):
                    infile = futures[future]
                    try:
                        key, status = future.result()
                        print(f"[{status}] {infile} -> {root}")
                    except KeyboardInterrupt:
                        STOP = True
                        failures += 1
                    except Exception as error:  # report all tasks, then return non-zero
                        failures += 1
                        print(f"[failed] {infile}: {error}", file=sys.stderr)
            except KeyboardInterrupt:
                STOP = True
                for future in futures:
                    future.cancel()
                print("Cancellation requested; active jobs will not be promoted to results.", file=sys.stderr)
                failures += 1
        _write_path_lists(root, args, inputs)
        # Failed jobs retain diagnostics in staging; never summarize stale successful reports.
        for report_name in ('sequence_changes',):
            report_paths = []
            for source in inputs:
                stage = root / '.msap-batch-staging' / _input_key(source)
                directory = stage if stage.exists() else root
                report_paths.append(directory / f'{source.stem}.{report_name}.tsv')
            merge_reports(report_paths, root / 'reports' / f'{report_name}.tsv')
        return 130 if STOP else (1 if failures else 0)


if __name__ == "__main__":
    raise SystemExit(main())
