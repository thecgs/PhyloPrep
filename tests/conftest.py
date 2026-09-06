import os
import sys
from pathlib import Path

import pytest


@pytest.fixture
def fake_aligner(tmp_path):
    """A controllable executable; exercises real CLI/process/file boundaries."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    executable = bin_dir / "mafft"
    executable.write_text(
        f"#!{sys.executable}\n"
        "import os, sys, time, signal\n"
        "from pathlib import Path\n"
        "with open(os.environ['MSAP_TEST_LOG'], 'a') as log: log.write('run\\n')\n"
        "if os.environ.get('MSAP_TEST_SLEEP'):\n"
        "    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "    Path(os.environ['MSAP_TEST_PID']).write_text(str(os.getpid()))\n"
        "    time.sleep(60)\n"
        "if os.environ.get('MSAP_TEST_FAIL'): sys.exit(7)\n"
        "if not os.environ.get('MSAP_TEST_EMPTY'):\n"
        "    print(Path(sys.argv[-1]).read_text(), end='')\n"
    )
    executable.chmod(0o755)
    env = {**os.environ, "PATH": str(bin_dir) + os.pathsep + os.environ.get("PATH", ""),
           "MSAP_TEST_LOG": str(tmp_path / "aligner.log"),
           "MSAP_TEST_PID": str(tmp_path / "aligner.pid")}
    return executable, env
