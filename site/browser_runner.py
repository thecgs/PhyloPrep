"""Execute a bundled PhyloPrep CLI script inside Pyodide without reimplementing it."""

from __future__ import annotations

import contextlib
import io
import importlib
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import traceback
import types
from collections import deque
from functools import wraps


_MEMORY_ERROR_RESULT = ('{"ok":false,"error":"PhyloPrep exhausted memory in the Python/Pyodide runtime. '
                        'Try fewer loci or the native CLI. Raising the MAFFT/MUSCLE heap limit does not '
                        'increase Pyodide memory. The original diagnostic may be unavailable.",'
                        '"stdout":"","stderr":""}')


class BoundedLog(io.TextIOBase):
    """Keep only a bounded diagnostic tail, including for many skipped loci."""
    def __init__(self, limit=65536):
        super().__init__()
        self.limit = limit
        self.parts = deque()
        self.length = 0
        self.truncated = False

    def write(self, text):
        count = len(text)
        if not count:
            return 0
        if count >= self.limit:
            self.parts.clear()
            self.parts.append(text[-self.limit:])
            self.truncated = self.truncated or count > self.limit or self.length > 0
            self.length = self.limit
        else:
            self.parts.append(text)
            self.length += count
            while self.length > self.limit:
                first = self.parts.popleft()
                excess = self.length - self.limit
                if len(first) > excess:
                    self.parts.appendleft(first[excess:])
                    self.length -= excess
                else:
                    self.length -= len(first)
                self.truncated = True
            while len(self.parts) > 1024:
                self.length -= len(self.parts.popleft())
                self.truncated = True
        return count

    def getvalue(self):
        return ('[Earlier log output truncated]\n' if self.truncated else '') + ''.join(self.parts)


def memory_safe_result(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except MemoryError:
            # Already serialized at import time: do not format a traceback,
            # inspect files or allocate another JSON payload under pressure.
            return _MEMORY_ERROR_RESULT
    return wrapped


def _exit_code(value):
    if value is None:
        return 0
    if isinstance(value, int):
        return value
    return 1


def _output_files(output_root: Path, manifest_only=False, path_lists=()):
    if manifest_only:
        return _export_output_files(output_root, _visible_output_files(output_root, manifest_only=True), path_lists=path_lists)
    files = []
    for path in sorted(output_root.rglob("*")):
        if path.is_file():
            files.append({
                "name": path.relative_to(output_root).as_posix(),
                "content": path.read_text(encoding="utf-8", errors="replace"),
            })
    return files


@memory_safe_result
def run_original_module(module, argv, output_root, engine_root="/engine", manifest_only=False):
    """Run an untouched source script as its CLI entry point and return its files."""
    engine = Path(engine_root).resolve()
    script = engine / module
    output = Path(output_root).resolve()
    if not script.is_file():
        return json.dumps({"ok": False, "error": f"Original script not found: {module}"})

    output.mkdir(parents=True, exist_ok=True)
    old_argv, old_cwd = sys.argv, Path.cwd()
    stdout, stderr = BoundedLog(), BoundedLog()
    code = 0
    try:
        if str(engine) not in sys.path:
            sys.path.insert(0, str(engine))
        sys.argv = [str(script), *argv]
        os.chdir(output)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                runpy.run_path(str(script), run_name="__main__")
            except SystemExit as error:
                code = _exit_code(error.code)
    except Exception:
        code = 1
        traceback.print_exc(file=stderr)
    finally:
        sys.argv = old_argv
        os.chdir(old_cwd)

    message = stderr.getvalue().strip() or stdout.getvalue().strip()
    if code:
        return json.dumps({
            "ok": False,
            "error": message or f"{module} exited with status {code}.",
            "stdout": stdout.getvalue(),
            "stderr": stderr.getvalue(),
        })
    list_paths = []
    for index, argument in enumerate(argv):
        if str(argument) == '--pathlist-output' and index + 1 < len(argv):
            list_paths.append(argv[index + 1])
        elif str(argument).startswith('--pathlist-output='):
            list_paths.append(str(argument).split('=', 1)[1])
    list_paths = [Path(value) if Path(value).is_absolute() else output / value for value in list_paths]
    files = _output_files(output, manifest_only=manifest_only, path_lists=list_paths)
    if not files:
        return json.dumps({"ok": False, "error": f"{module} did not create an output file."})
    return json.dumps({
        "ok": True,
        "files": files,
        "summary": f"Created {len(files)} output file(s) with the original {module} script.",
        "stdout": stdout.getvalue(),
        "stderr": stderr.getvalue(),
    })


class _BrowserAlignmentRequest(Exception):
    """Signal that the browser must satisfy one external alignment command."""

    def __init__(self, command):
        super().__init__("Alignment execution delegated to the browser")
        self.command = [str(item) for item in command]


def _browser_program_name(command):
    if not isinstance(command, (list, tuple)) or not command:
        return ""
    return os.path.basename(str(command[0]))


def _browser_prepare_executable(name):
    bindir = Path("/web-bin")
    bindir.mkdir(parents=True, exist_ok=True)
    executable = bindir / name
    if not executable.exists():
        executable.write_text("# browser placeholder\n", encoding="utf-8")
    os.chmod(executable, 0o755)
    paths = os.environ.get("PATH", os.defpath).split(os.pathsep)
    if str(bindir) not in paths:
        os.environ["PATH"] = os.pathsep.join([str(bindir), *paths])


def _browser_python_child(command, kwargs):
    script = Path(str(command[1]))
    old_argv = sys.argv
    code = 0
    try:
        sys.argv = [str(script), *[str(item) for item in command[2:]]]
        try:
            runpy.run_path(str(script), run_name="__main__")
        except SystemExit as error:
            code = _exit_code(error.code)
    finally:
        sys.argv = old_argv
    if code and kwargs.get("check"):
        raise subprocess.CalledProcessError(code, command)
    return subprocess.CompletedProcess(
        command,
        code,
        stdout="" if kwargs.get("stdout") == subprocess.PIPE else None,
        stderr="" if kwargs.get("stderr") == subprocess.PIPE else None,
    )


def _browser_execute_msap(argv, output_root, handler, engine_root="/engine"):
    engine = Path(engine_root).resolve()
    script = engine / "MSAP.py"
    output = Path(output_root).resolve()
    if not script.is_file():
        raise FileNotFoundError("Original MSAP.py is missing from the browser engine.")
    output.mkdir(parents=True, exist_ok=True)
    _browser_prepare_executable("muscle")
    _browser_prepare_executable("mafft")

    old_argv, old_cwd = sys.argv, Path.cwd()
    original_run = subprocess.run
    stdout, stderr = BoundedLog(), BoundedLog()
    code = 0
    try:
        if str(engine) not in sys.path:
            sys.path.insert(0, str(engine))
        sys.argv = [str(script), *argv]
        os.chdir(output)
        subprocess.run = handler
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                runpy.run_path(str(script), run_name="__main__")
            except SystemExit as error:
                code = _exit_code(error.code)
    finally:
        subprocess.run = original_run
        sys.argv = old_argv
        os.chdir(old_cwd)
    return code, stdout.getvalue(), stderr.getvalue()


@memory_safe_result
def run_msap_probe(argv, output_root, engine_root="/engine"):
    """Run untouched MSAP.py until it requests a browser-supported aligner."""
    def handler(command, *args, **kwargs):
        program = _browser_program_name(command)
        if program in {"muscle", "mafft"}:
            raise _BrowserAlignmentRequest(command)
        raise RuntimeError(
            "The browser MSAP workflow encountered an unsupported external command: "
            + " ".join(str(item) for item in command)
        )

    try:
        code, stdout, stderr = _browser_execute_msap(argv, output_root, handler, engine_root)
    except _BrowserAlignmentRequest as request:
        command = request.command
        program = _browser_program_name(command)
        try:
            if program == "muscle":
                input_path = command[command.index("-align") + 1]
                expected_output = Path(command[command.index("-output") + 1]).name
                version = "5.1.0"
            elif program == "mafft":
                input_path = command[-1]
                expected_output = None
                version = "7.520"
            else:
                raise ValueError(program)
        except (ValueError, IndexError):
            return json.dumps({"ok": False, "error": f"Could not interpret the {program} command generated by MSAP.py."})
        path = Path(input_path)
        if not path.is_file():
            return json.dumps({"ok": False, "error": f"MSAP.py prepared no {program} input at {input_path}."})
        return json.dumps({
            "ok": True,
            "tool": program,
            "version": version,
            "alignmentInput": path.read_text(encoding="utf-8", errors="replace"),
            "inputName": path.name,
            "expectedOutput": expected_output,
            "command": command,
        })
    except Exception:
        captured = BoundedLog()
        traceback.print_exc(file=captured)
        return json.dumps({"ok": False, "error": captured.getvalue().strip()})

    message = stderr.strip() or stdout.strip()
    return json.dumps({
        "ok": False,
        "error": message or f"MSAP.py exited with status {code} before requesting a supported browser aligner.",
    })


@memory_safe_result
def run_msap_complete(argv, output_root, alignment_file, engine_root="/engine", manifest_only=False):
    """Rerun untouched MSAP.py and satisfy its browser-supported alignment call."""
    alignment = Path(alignment_file)
    if not alignment.is_file():
        return json.dumps({"ok": False, "error": "The browser alignment is missing."})

    def handler(command, *args, **kwargs):
        program = _browser_program_name(command)
        if program == "muscle":
            command = [str(item) for item in command]
            try:
                output_path = command[command.index("-output") + 1]
            except (ValueError, IndexError):
                raise RuntimeError("Could not interpret the MUSCLE output path generated by MSAP.py.")
            shutil.copyfile(alignment, output_path)
            return subprocess.CompletedProcess(command, 0)
        if program == "mafft":
            destination = kwargs.get("stdout")
            if destination is None or not hasattr(destination, "write"):
                raise RuntimeError("MSAP.py did not provide a writable stdout destination for MAFFT.")
            with alignment.open(encoding="utf-8", errors="replace") as source:
                shutil.copyfileobj(source, destination, length=65536)
            destination.flush()
            return subprocess.CompletedProcess(command, 0)
        if isinstance(command, (list, tuple)) and len(command) > 1 and str(command[1]).endswith(".py"):
            return _browser_python_child(command, kwargs)
        raise RuntimeError(
            "This browser workflow supports MAFFT/MUSCLE plus PhyloPrep's bundled Python helpers only; "
            f"unsupported external command: {' '.join(str(item) for item in command)}"
        )

    try:
        code, stdout, stderr = _browser_execute_msap(argv, output_root, handler, engine_root)
    except Exception:
        captured = BoundedLog()
        traceback.print_exc(file=captured)
        return json.dumps({"ok": False, "error": captured.getvalue().strip()})
    if code:
        return json.dumps({
            "ok": False,
            "error": stderr.strip() or stdout.strip() or f"MSAP.py exited with status {code}.",
            "stdout": stdout,
            "stderr": stderr,
        })
    files = _output_files(Path(output_root).resolve(), manifest_only=manifest_only)
    if not files:
        return json.dumps({"ok": False, "error": "MSAP.py did not create output files."})
    return json.dumps({
        "ok": True,
        "files": files,
        "summary": f"Created {len(files)} output file(s) with the untouched MSAP.py workflow and browser alignment.",
        "stdout": stdout,
        "stderr": stderr,
    })



def _browser_install_fcntl_stub():
    """Allow importing CLI orchestrators when Pyodide lacks POSIX fcntl."""
    try:
        import fcntl  # noqa: F401
        return
    except ImportError:
        module = types.ModuleType("fcntl")
        module.LOCK_EX = 2
        module.LOCK_NB = 4
        module.LOCK_UN = 8
        module.flock = lambda *_args, **_kwargs: None
        sys.modules["fcntl"] = module


def _browser_run_engine_script(engine, script, arguments, cwd=None):
    """Run one bundled Python CLI module in the current Pyodide interpreter."""
    path = Path(engine) / script
    if not path.is_file():
        raise FileNotFoundError(f"Bundled script not found: {script}")
    old_argv, old_cwd = sys.argv, Path.cwd()
    code = 0
    try:
        sys.argv = [str(path), *[str(value) for value in arguments]]
        if cwd is not None:
            os.chdir(cwd)
        try:
            runpy.run_path(str(path), run_name="__main__")
        except SystemExit as error:
            code = _exit_code(error.code)
    finally:
        sys.argv = old_argv
        os.chdir(old_cwd)
    if code:
        raise subprocess.CalledProcessError(code, [sys.executable, str(path), *map(str, arguments)])


def _visible_output_files(root, manifest_only=False):
    root = Path(root)
    files = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if not path.is_file() or any(part.startswith(".") for part in relative.parts):
            continue
        file = {"name": relative.as_posix()}
        if manifest_only:
            file.update(path=str(path), size=path.stat().st_size)
        else:
            file["content"] = path.read_text(encoding="utf-8", errors="replace")
        files.append(file)
    return files


def _export_output_files(root, files, path_lists=()):
    """Prepare portable list copies; never change paths used by the engine."""
    root = Path(root).resolve()
    list_paths = {Path(path).resolve() for path in path_lists}
    for item in files:
        path = Path(item['path'])
        if path.resolve() in list_paths or path.suffix in {'.pathlist', '.list'} or path.name.endswith(('.pass.tsv', '.fail.tsv', '.retained.tsv')):
            exported = root / '.browser-export' / item['name']
            exported.parent.mkdir(parents=True, exist_ok=True)
            with path.open(encoding='utf-8') as source, exported.open('w', encoding='utf-8') as target:
                for line in source:
                    entry = line.strip()
                    candidate = Path(entry)
                    if entry and not entry.startswith(('#', ';')):
                        if not candidate.is_absolute():
                            candidate = path.parent / candidate
                        candidate = candidate.resolve()
                        if candidate.is_relative_to(root):
                            # Completed branches may already be safely stored on disk.
                            entry = Path(os.path.relpath(candidate, path.parent)).as_posix()
                        elif candidate.is_file():
                            # Inputs are not bundled: place the original FASTA
                            # alongside this exported list before CLI reuse.
                            entry = candidate.name
                    target.write(entry + '\n')
            path = exported
            item.update(path=str(path), size=path.stat().st_size)
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        item['sha256'] = digest.hexdigest()
    return files


def _phyloprep_finalize_steps(file_manifest, options, output_root, engine_root="/engine", progress_callback=None, manifest_only=False, incremental=False, disk_backed=False):
    """Finish the browser pipeline with the repository's original downstream modules."""
    engine = Path(engine_root).resolve()
    root = Path(output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    stdout, stderr = BoundedLog(), BoundedLog()
    mapping_supplied = False
    phyloprep = None
    archived = {}

    def current_outputs():
        files = _visible_output_files(root, manifest_only=manifest_only)
        if manifest_only:
            files = _export_output_files(root, files)
            return list({**archived, **{item['name']: item for item in files}}.values())
        return files

    try:
        if str(engine) not in sys.path:
            sys.path.insert(0, str(engine))
        _browser_install_fcntl_stub()
        phyloprep = importlib.import_module("phyloprep")
        sequence_audit = importlib.import_module("sequence_audit")
        rename_taxa = importlib.import_module("rename_taxa")

        sources = [Path(item["path"]) for item in file_manifest if item.get("role") == "source"]
        batch_outputs = [item for item in file_manifest if item.get("role") == "batch-output"]
        mapping_items = [item for item in file_manifest if item.get("role") == "mapping"]
        if not sources:
            raise ValueError("No source FASTA files were supplied to the PhyloPrep browser finalizer.")
        phyloprep.write_paths(root / 'inputs.pathlist', [source.name for source in sources])
        stems = [path.stem for path in sources]
        if len(set(stems)) != len(stems):
            raise ValueError("Input FASTA files must have unique basenames.")

        alignments = root / "alignments"
        alignments.mkdir(parents=True, exist_ok=True)
        change_reports = []
        alignment_targets = []
        for item in batch_outputs:
            source = Path(item["path"])
            target = alignments / Path(item["name"]).name
            if source.exists():
                if disk_backed:
                    source.replace(target)
                else:
                    shutil.copyfile(source, target)
            elif not (disk_backed and item.get("on_disk")):
                raise FileNotFoundError(source)
            alignment_targets.append(target)
            if target.name.endswith(".sequence_changes.tsv"):
                change_reports.append(target)

        path_lists = alignments / "path-lists"
        path_lists.mkdir(parents=True, exist_ok=True)
        for kind in ("nucl", "prot", "codon"):
            for variant, suffix in (("raw", f".{kind}.aln.fasta"), ("trimmed", f".{kind}.trimmed.aln.fasta")):
                files = sorted(
                    path.resolve() for path in alignment_targets
                    if path.name.endswith(suffix)
                )
                if files:
                    list_name = f"all.{kind}.aln.fasta.pathlist" if variant == "raw" else f"all.{kind}.trimmed.aln.fasta.pathlist"
                    phyloprep.write_paths(path_lists / list_name, files)

        reports = root / "reports"
        reports.mkdir(parents=True, exist_ok=True)
        merged_changes = reports / "sequence_changes.tsv"
        sequence_audit.merge_reports(change_reports, merged_changes)

        mapping = Path(mapping_items[0]["path"]) if mapping_items else None
        mapping_supplied = mapping is not None
        taxon_mapping, mapping_signature = (rename_taxa.read_mapping_snapshot(mapping)
                                            if mapping else (None, None))
        if mapping:
            phyloprep.annotate_change_report(merged_changes, taxon_mapping)

        def number(name, default, kind=float):
            value = options.get(name, default)
            return kind(value)

        def optional_number(name, kind=int):
            value = options.get(name)
            if value in (None, ""):
                return None
            return kind(value)

        args = types.SimpleNamespace(
            seqtype=str(options.get("seqtype", "codon")),
            align_software=str(options.get("align_software", "mafft")),
            notrim=bool(options.get("notrim", False)),
            genetic_code=int(options.get("genetic_code", 1)),
            G=number("gap_ratio", 0.2),
            N=number("n_ratio", 0.2),
            X=number("x_ratio", 0.2),
            mapping=mapping,
            taxon_mapping=taxon_mapping,
            mapping_signature=mapping_signature,
            no_resume=True,
            no_report=bool(options.get('no_report', False)),
            _report_dependency_error=options.get('report_dependency_error'),
            _report_context=options.get('report_context'),
            _browser_runtime=True,
            missing_taxa=str(options.get("missing_taxa", "skip-gene")),
            min_length=number("min_length", 1, int),
            max_n_ratio=number("max_n_ratio", 1.0),
            max_x_ratio=number("max_x_ratio", 1.0),
            max_gap_ratio=number("max_gap_ratio", 1.0),
            min_nucl_pis=optional_number("min_nucl_pis", int),
            min_protein_pis=optional_number("min_protein_pis", int),
            paired_pis_policy=str(options.get("paired_pis_policy", "both")),
            allow_internal_stop=bool(options.get("allow_internal_stop", False)),
        )
        phyloprep.validate_pipeline_options(args)
        from msap_progress import ProgressReporter, progress_context, report_phase
        variants = ("raw",) if args.notrim else ("raw", "trimmed")
        args._progress = ProgressReporter('PhyloPrep', len(variants) + 1 + (not args.no_report), unit='stages', callback=progress_callback)

        def browser_run(script, *arguments, cwd=None):
            print("[phyloprep-web] " + " ".join([script, *map(str, arguments)]), flush=True)
            return _browser_run_engine_script(engine, script, arguments, cwd=cwd)

        @contextlib.contextmanager
        def analysis_step():
            original_run = phyloprep.run
            phyloprep.run = browser_run
            try:
                with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr), progress_context(args._progress):
                    yield
            finally:
                phyloprep.run = original_run

        with tempfile.TemporaryDirectory(prefix=".matrices-", dir=root) as temporary:
            matrix_root = Path(temporary)
            for index, variant in enumerate(variants):
                if disk_backed:
                    yield {"load_variant": variant, "files": []}
                if incremental and mapping:
                    args._prepared_renamed = {}
                    kinds = ('codon', 'prot') if args.seqtype == 'codon' else (args.seqtype,)
                    for kind in kinds:
                        paths = phyloprep.alignment_paths(args, root, kind, variant)
                        renamed = []
                        iterator = phyloprep.iter_renamed_alignments(args, paths, root, kind, variant)
                        while True:
                            with analysis_step():
                                try:
                                    target = next(iterator)
                                except StopIteration:
                                    break
                                report_phase(args, f'{variant}/{kind}: renamed {target.name}')
                            renamed.append(target)
                            item = dict(name=target.relative_to(root).as_posix(), path=str(target), size=target.stat().st_size)
                            checkpoint = _export_output_files(root, [item])
                            archived[item['name']] = dict(item)
                            # Keep the OPFS working copy for QC; saving is now per locus.
                            yield checkpoint
                        args._prepared_renamed[(kind, variant)] = renamed
                # Each runPythonAsync call can have a different Context. Leave
                # progress/log/engine scopes before yielding files to JavaScript.
                with analysis_step():
                    report_phase(args, f'{variant}: QC and matrix preparation')
                    phyloprep.prepare_variant(args, root, matrix_root, variant)
                    args._progress.update(f'{variant}: complete', completed=index + 1)
                if incremental:
                    # Staged matrices stay private until publication succeeds.
                    checkpoint = []
                    for folder in ('qc', 'reports'):
                        for item in _visible_output_files(root / folder, manifest_only=True):
                            item['name'] = folder + '/' + item['name']
                            checkpoint.append(item)
                    if disk_backed:
                        for folder in ('alignments', 'renamed'):
                            for item in _visible_output_files(root / folder, manifest_only=True):
                                if not item['name'].endswith('.aln.fasta'):
                                    continue
                                trimmed = item['name'].endswith('.trimmed.aln.fasta')
                                if trimmed != (variant == 'trimmed'):
                                    continue
                                item.update(name=folder + '/' + item['name'], release=True)
                                checkpoint.append(item)
                    checkpoint = _export_output_files(root, checkpoint)
                    archived.update((item['name'], dict(item)) for item in checkpoint)
                    yield checkpoint
            with analysis_step():
                report_phase(args, 'Publishing matrices')
                phyloprep.publish_matrices(matrix_root, root, variants, block_signals=False)
                args._progress.update('Matrices published', completed=len(variants) + 1)

        with analysis_step():
            if not args.no_report:
                report_phase(args, 'Generating HTML report')
            report_result = phyloprep.generate_pipeline_report(args, root)
            args._progress.update('Complete', completed=args._progress.total)

        files = current_outputs()
        if not files:
            raise ValueError("PhyloPrep produced no downloadable output files.")
        return json.dumps({
            "ok": True,
            "files": files,
            "summary": f"Created {len(files)} files with the browser PhyloPrep pipeline.",
            "report": report_result,
            "warnings": [report_result['warning']] if report_result['status'] == 'failed' else [],
            "stdout": stdout.getvalue(),
            "stderr": stderr.getvalue(),
        })
    except MemoryError:
        return _MEMORY_ERROR_RESULT
    except Exception as error:
        if phyloprep is not None and isinstance(error, phyloprep.NoAlignmentsPassedQCError):
            return json.dumps({
                "ok": False,
                "recoverable": True,
                "files": current_outputs(),
                "user_error": str(error),
                "summary": "No loci passed QC; diagnostic reports were preserved.",
                "stdout": stdout.getvalue(),
                "stderr": stderr.getvalue(),
            })
        traceback.print_exc(file=stderr)
        technical_error = stderr.getvalue().strip() or "The browser PhyloPrep pipeline failed."
        if phyloprep is not None and isinstance(error, phyloprep.NoGenesIncludedError):
            files = current_outputs()
            if mapping_supplied:
                user_error = (
                    "No supermatrix could be built with missing-taxa=skip-gene because no locus retained the full "
                    "taxon set. Inspect supermatrix.report.tsv. If taxa are genuinely absent, use pad-gaps; if the "
                    "mapping is incorrect, correct map.tsv and rerun."
                )
            else:
                user_error = (
                    "No supermatrix could be built because no map.tsv was provided and the FASTA IDs do not resolve "
                    "to a shared taxon set across loci. This commonly happens with locus-specific GenBank/CDS IDs. "
                    "Provide a two-column map.tsv (FASTA_ID<TAB>Taxon) and rerun. Use pad-gaps only for genuinely "
                    "missing taxa, not to replace taxon mapping. Diagnostic reports and intermediate outputs are "
                    "available below."
                )
            return json.dumps({
                "ok": False,
                "recoverable": True,
                "files": files,
                "user_error": user_error,
                "summary": "PhyloPrep stopped before supermatrix completion; diagnostic files were preserved.",
                "technical_error": technical_error,
                "stdout": stdout.getvalue(),
                "stderr": stderr.getvalue(),
            })
        return json.dumps({
            "ok": False,
            "error": technical_error,
            "stdout": stdout.getvalue(),
            "stderr": stderr.getvalue(),
        })


@memory_safe_result
def advance_phyloprep_finalize(steps):
    """One synchronous analysis step, followed by bounded asynchronous JS export."""
    try:
        files = next(steps)
    except StopIteration as completed:
        return completed.value
    return json.dumps({"pending": True, **files} if isinstance(files, dict)
                      else {"pending": True, "files": files})


@memory_safe_result
def run_phyloprep_finalize(file_manifest, options, output_root, engine_root="/engine", progress_callback=None, manifest_only=False):
    """Existing synchronous interface; incremental exports are browser opt-in."""
    steps = _phyloprep_finalize_steps(file_manifest, options, output_root, engine_root,
                                     progress_callback, manifest_only)
    try:
        while True:
            next(steps)
    except StopIteration as completed:
        return completed.value
    finally:
        steps.close()
