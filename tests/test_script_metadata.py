import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("script_metadata", ROOT / "site/script_metadata.py")
metadata = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metadata)


def test_capture_version_without_running_work(tmp_path):
    script = tmp_path / "tool.py"
    script.write_text("import argparse\np = argparse.ArgumentParser()\np.add_argument('--version', action='version', version='v9.8.7')\np.parse_args()\nraise RuntimeError('work must not execute')\n")
    original = argparse.ArgumentParser.parse_args
    assert metadata.collect_script_metadata([script])[script.name]["version"] == "v9.8.7"
    assert argparse.ArgumentParser.parse_args is original


def test_restore_parser_after_import_failure(tmp_path):
    script = tmp_path / "broken.py"
    script.write_text("raise RuntimeError('broken import')\n")
    original = argparse.ArgumentParser.parse_args
    with pytest.raises(RuntimeError, match="broken import"):
        metadata.collect_script_metadata([script])
    assert argparse.ArgumentParser.parse_args is original


@pytest.mark.parametrize("script", ["oma_to_orthogroups.py", "orthomcl_to_orthogroups.py", "proteinortho_to_orthogroups.py", "sonicparanoid_to_orthogroups.py"])
def test_adapter_version_matches_cli(script):
    declared = metadata.collect_script_metadata([ROOT / script])[script]["version"]
    actual = subprocess.run([sys.executable, str(ROOT / script), "--version"], check=True, capture_output=True, text=True)
    assert actual.stdout.strip() == declared
