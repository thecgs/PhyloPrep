"""Execute a bundled PhyloPrep CLI script inside Pyodide without reimplementing it."""

from __future__ import annotations

import contextlib
import io
import importlib
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import traceback
import types


def _exit_code(value):
    if value is None:
        return 0
    if isinstance(value, int):
        return value
    return 1


def _output_files(output_root: Path):
    files = []
    for path in sorted(output_root.rglob("*")):
        if path.is_file():
            files.append({
                "name": path.relative_to(output_root).as_posix(),
                "content": path.read_text(encoding="utf-8", errors="replace"),
            })
    return files


def run_original_module(module, argv, output_root, engine_root="/engine"):
    """Run an untouched source script as its CLI entry point and return its files."""
    engine = Path(engine_root).resolve()
    script = engine / module
    output = Path(output_root).resolve()
    if not script.is_file():
        return json.dumps({"ok": False, "error": f"Original script not found: {module}"})

    output.mkdir(parents=True, exist_ok=True)
    old_argv, old_cwd = sys.argv, Path.cwd()
    stdout, stderr = io.StringIO(), io.StringIO()
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
    files = _output_files(output)
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
    stdout, stderr = io.StringIO(), io.StringIO()
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
        captured = io.StringIO()
        traceback.print_exc(file=captured)
        return json.dumps({"ok": False, "error": captured.getvalue().strip()})

    message = stderr.strip() or stdout.strip()
    return json.dumps({
        "ok": False,
        "error": message or f"MSAP.py exited with status {code} before requesting a supported browser aligner.",
    })


def run_msap_complete(argv, output_root, alignment_file, engine_root="/engine"):
    """Rerun untouched MSAP.py and satisfy its browser-supported alignment call."""
    alignment = Path(alignment_file)
    if not alignment.is_file():
        return json.dumps({"ok": False, "error": "The browser alignment is missing."})
    alignment_text = alignment.read_text(encoding="utf-8", errors="replace")

    def handler(command, *args, **kwargs):
        program = _browser_program_name(command)
        if program == "muscle":
            command = [str(item) for item in command]
            try:
                output_path = command[command.index("-output") + 1]
            except (ValueError, IndexError):
                raise RuntimeError("Could not interpret the MUSCLE output path generated by MSAP.py.")
            Path(output_path).write_text(alignment_text, encoding="utf-8")
            return subprocess.CompletedProcess(command, 0)
        if program == "mafft":
            destination = kwargs.get("stdout")
            if destination is None or not hasattr(destination, "write"):
                raise RuntimeError("MSAP.py did not provide a writable stdout destination for MAFFT.")
            destination.write(alignment_text)
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
        captured = io.StringIO()
        traceback.print_exc(file=captured)
        return json.dumps({"ok": False, "error": captured.getvalue().strip()})
    if code:
        return json.dumps({
            "ok": False,
            "error": stderr.strip() or stdout.strip() or f"MSAP.py exited with status {code}.",
            "stdout": stdout,
            "stderr": stderr,
        })
    files = _output_files(Path(output_root).resolve())
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


def _visible_output_files(root):
    return [
        file for file in _output_files(Path(root))
        if not any(part.startswith(".") for part in Path(file["name"]).parts)
    ]


def run_phyloprep_finalize(file_manifest, options, output_root, engine_root="/engine"):
    """Finish the browser pipeline with the repository's original downstream modules."""
    engine = Path(engine_root).resolve()
    root = Path(output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    stdout, stderr = io.StringIO(), io.StringIO()
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
        stems = [path.stem for path in sources]
        if len(set(stems)) != len(stems):
            raise ValueError("Input FASTA files must have unique basenames.")

        alignments = root / "alignments"
        alignments.mkdir(parents=True, exist_ok=True)
        change_reports = []
        for item in batch_outputs:
            source = Path(item["path"])
            target = alignments / Path(item["name"]).name
            shutil.copyfile(source, target)
            if target.name.endswith(".sequence_changes.tsv"):
                change_reports.append(target)

        reports = root / "reports"
        reports.mkdir(parents=True, exist_ok=True)
        merged_changes = reports / "sequence_changes.tsv"
        sequence_audit.merge_reports(change_reports, merged_changes)

        mapping = Path(mapping_items[0]["path"]) if mapping_items else None
        taxon_mapping = rename_taxa.read_mapping(mapping) if mapping else None
        if mapping:
            phyloprep.annotate_change_report(merged_changes, taxon_mapping)

        def number(name, default, kind=float):
            value = options.get(name, default)
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
            no_resume=True,
            missing_taxa=str(options.get("missing_taxa", "skip-gene")),
            min_length=number("min_length", 1, int),
            max_n_ratio=number("max_n_ratio", 1.0),
            max_x_ratio=number("max_x_ratio", 1.0),
            max_gap_ratio=number("max_gap_ratio", 1.0),
            min_parsimony_informative_sites=number("min_parsimony_informative_sites", 0, int),
            allow_terminal_stop=bool(options.get("allow_terminal_stop", False)),
            allow_internal_stop=bool(options.get("allow_internal_stop", False)),
        )
        ratios = [args.G, args.N, args.X, args.max_n_ratio, args.max_x_ratio, args.max_gap_ratio]
        if any(value < 0 or value > 1 for value in ratios):
            raise ValueError("Trimming and QC ratios must be between 0 and 1.")
        if args.min_length < 0 or args.min_parsimony_informative_sites < 0:
            raise ValueError("QC length/site thresholds must be non-negative.")

        original_run = phyloprep.run

        def browser_run(script, *arguments, cwd=None):
            print("[phyloprep-web] " + " ".join([script, *map(str, arguments)]), flush=True)
            return _browser_run_engine_script(engine, script, arguments, cwd=cwd)

        phyloprep.run = browser_run
        try:
            seqtypes = ("codon", "prot") if args.seqtype == "codon" else (args.seqtype,)
            variants = ("raw",) if args.notrim else ("raw", "trimmed")
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                for seqtype in seqtypes:
                    for variant in variants:
                        stage = root / "matrices" / seqtype / variant
                        stage.mkdir(parents=True, exist_ok=True)
                        phyloprep.prepare_group(args, sources, root, stage, seqtype, variant)
        finally:
            phyloprep.run = original_run

        files = _visible_output_files(root)
        if not files:
            raise ValueError("PhyloPrep produced no downloadable output files.")
        return json.dumps({
            "ok": True,
            "files": files,
            "summary": f"Created {len(files)} files with the browser PhyloPrep pipeline.",
            "stdout": stdout.getvalue(),
            "stderr": stderr.getvalue(),
        })
    except Exception:
        traceback.print_exc(file=stderr)
        return json.dumps({
            "ok": False,
            "error": stderr.getvalue().strip() or "The browser PhyloPrep pipeline failed.",
            "stdout": stdout.getvalue(),
            "stderr": stderr.getvalue(),
        })
