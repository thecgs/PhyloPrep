"""Execute a bundled PhyloPrep CLI script inside Pyodide without reimplementing it."""

from __future__ import annotations

import contextlib
import io
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import traceback


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


class _BrowserMuscleRequest(Exception):
    """Signal that the browser must satisfy one MUSCLE command."""

    def __init__(self, command):
        super().__init__("MUSCLE execution delegated to the browser")
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
    """Run untouched MSAP.py until it requests MUSCLE and return that input."""
    def handler(command, *args, **kwargs):
        if _browser_program_name(command) == "muscle":
            raise _BrowserMuscleRequest(command)
        raise RuntimeError(
            "The browser MSAP workflow encountered an external command before MUSCLE: "
            + " ".join(str(item) for item in command)
        )

    try:
        code, stdout, stderr = _browser_execute_msap(argv, output_root, handler, engine_root)
    except _BrowserMuscleRequest as request:
        command = request.command
        try:
            input_path = command[command.index("-align") + 1]
            output_path = command[command.index("-output") + 1]
        except (ValueError, IndexError):
            return json.dumps({"ok": False, "error": "Could not interpret the MUSCLE command generated by MSAP.py."})
        path = Path(input_path)
        if not path.is_file():
            return json.dumps({"ok": False, "error": f"MSAP.py prepared no MUSCLE input at {input_path}."})
        return json.dumps({
            "ok": True,
            "tool": "muscle",
            "version": "5.1.0",
            "alignmentInput": path.read_text(encoding="utf-8", errors="replace"),
            "inputName": path.name,
            "expectedOutput": Path(output_path).name,
            "command": command,
        })
    except Exception:
        captured = io.StringIO()
        traceback.print_exc(file=captured)
        return json.dumps({"ok": False, "error": captured.getvalue().strip()})

    message = stderr.strip() or stdout.strip()
    return json.dumps({
        "ok": False,
        "error": message or f"MSAP.py exited with status {code} before requesting MUSCLE.",
    })


def run_msap_complete(argv, output_root, alignment_file, engine_root="/engine"):
    """Rerun untouched MSAP.py and satisfy its MUSCLE call with a browser result."""
    alignment = Path(alignment_file)
    if not alignment.is_file():
        return json.dumps({"ok": False, "error": "The browser MUSCLE alignment is missing."})

    def handler(command, *args, **kwargs):
        program = _browser_program_name(command)
        if program == "muscle":
            command = [str(item) for item in command]
            try:
                output_path = command[command.index("-output") + 1]
            except (ValueError, IndexError):
                raise RuntimeError("Could not interpret the MUSCLE output path generated by MSAP.py.")
            Path(output_path).write_text(
                alignment.read_text(encoding="utf-8", errors="replace"),
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(command, 0)
        if isinstance(command, (list, tuple)) and len(command) > 1 and str(command[1]).endswith(".py"):
            return _browser_python_child(command, kwargs)
        raise RuntimeError(
            "This browser workflow supports MUSCLE plus PhyloPrep's bundled Python helpers only; "
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
        "summary": f"Created {len(files)} output file(s) with the untouched MSAP.py workflow and browser MUSCLE.",
        "stdout": stdout,
        "stderr": stderr,
    })
