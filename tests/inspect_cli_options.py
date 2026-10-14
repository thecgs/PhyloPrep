"""Expose production CLI metadata to the browser parity tests."""
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('script_metadata', ROOT / 'site/script_metadata.py')
metadata = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metadata)
result = metadata.collect_script_metadata([ROOT / name for name in sys.argv[1:]])
print(json.dumps({name: entry['actions'] for name, entry in result.items()}, default=str))
