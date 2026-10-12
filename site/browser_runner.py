"""Execute a bundled PhyloPrep CLI script inside Pyodide without reimplementing it."""

from __future__ import annotations

import contextlib
import io
import json
import os
from pathlib import Path
import runpy
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
