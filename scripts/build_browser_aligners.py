"""Rebuild pinned Biowasm aligners with a 4 GiB wasm32 heap ceiling."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BIOWASM_REV = '4afa546aa207172355fdf2f639167ea78ace539a'
EMSDK_VERSION = '2.0.25'
MAXIMUM_MEMORY = 4294967296
TOOLS = {
    'mafft': ('7.520', 'v7.520', '52b59f064c600da59bca8233736418fb8bb35d5e', ['tbfast', 'dvtditr']),
    'muscle': ('5.1.0', '5.1.0', '6c601163998616bb88991931e443c645858e162c', ['muscle']),
}


def run(*args, cwd=None, **kwargs):
    return subprocess.run(args, cwd=cwd, check=True, **kwargs)


def build(work_dir, output):
    if not shutil.which('emcc'):
        raise SystemExit('emcc is required; activate Emscripten 2.0.25 or use the Pages CI build.')
    compiler = subprocess.check_output(['emcc', '--version'], text=True)
    if EMSDK_VERSION not in compiler.splitlines()[0]:
        raise SystemExit(f'Expected Emscripten {EMSDK_VERSION}, got {compiler.splitlines()[0]}')
    work_dir.mkdir(parents=True, exist_ok=True)
    source = work_dir / 'biowasm'
    if source.exists():
        raise SystemExit(f'Use an empty build directory; refusing to overwrite {source}')
    run('git', 'clone', '--no-checkout', 'https://github.com/biowasm/biowasm.git', str(source))
    run('git', 'checkout', '--detach', BIOWASM_REV, cwd=source)
    run('git', 'submodule', 'update', '--init', 'tools/mafft/src', 'tools/muscle/src', cwd=source)
    shared = source / 'bin/shared.sh'
    contents = shared.read_text()
    needle = '-s ALLOW_MEMORY_GROWTH=1'
    if contents.count(needle) != 1:
        raise RuntimeError('Unexpected upstream memory flags')
    shared.write_text(contents.replace(needle, needle + '\n\t-s MAXIMUM_MEMORY=4294967296'))
    output.mkdir(parents=True, exist_ok=True)
    manifest = {'biowasmRevision': BIOWASM_REV, 'emscripten': EMSDK_VERSION,
                'maximumMemoryBytes': MAXIMUM_MEMORY, 'allowMemoryGrowth': True, 'tools': {}}
    env = dict(os.environ, MAKEFLAGS='-j2')
    for tool, (version, tag, revision, programs) in TOOLS.items():
        tool_dir = source / 'tools' / tool
        src = tool_dir / 'src'
        run('git', 'checkout', '--detach', revision, cwd=src)
        run('git', 'apply', str(tool_dir / 'patches' / f'{tag}.patch'), cwd=src)
        (tool_dir / 'build').mkdir(exist_ok=True)
        run('bash', '-e', '-c', f'. bin/shared.sh; cd tools/{tool}/src; bash -e ../compile.sh', cwd=source, env=env)
        dest = output / tool / version
        dest.mkdir(parents=True, exist_ok=True)
        files = []
        for program in programs:
            for ext in ['js', 'wasm']:
                original = tool_dir / 'build' / f'{program}.{ext}'
                if not original.is_file() or not original.stat().st_size:
                    raise RuntimeError(f'Missing build output: {original}')
                shutil.copy2(original, dest / original.name)
                files.append({'name': original.name, 'sha256': hashlib.sha256(original.read_bytes()).hexdigest()})
        for license_name in ['license', 'license.extensions', 'LICENSE']:
            if (src / license_name).exists():
                shutil.copy2(src / license_name, dest / license_name)
        # Ship the exact patches and source revisions alongside the binaries.
        shutil.copy2(tool_dir / 'patches' / f'{tag}.patch', dest / 'biowasm.patch')
        manifest['tools'][tool] = {'packageVersion': version, 'sourceRevision': revision, 'files': files}
    shutil.copy2(source / 'LICENSE', output / 'BIOWASM-LICENSE')
    # A matching small program proves that the compiler/runtime really crosses 2 GiB.
    probe = work_dir / 'memory_probe.c'
    probe.write_text('''#include <stdio.h>
#include <stdlib.h>
#include <emscripten/heap.h>
int main(void) {
    size_t n = (size_t)1100 * 1024 * 1024;
    volatile unsigned char *p = malloc(n);
    volatile unsigned char *q = malloc(n);
    if (!p || !q) { puts("PROBE_FAILED"); return 1; }
    p[0] = 17; p[n-1] = 19; q[0] = 23; q[n-1] = 29;
    printf("HEAP_BYTES=%zu END=%u\\n", emscripten_get_heap_size(), (unsigned)q[n-1]);
    free((void *)q); free((void *)p); return 0;
}
''')
    dest = output / 'probe' / '1'
    dest.mkdir(parents=True, exist_ok=True)
    run('emcc', str(probe), '-O2', '-sALLOW_MEMORY_GROWTH=1', '-sMAXIMUM_MEMORY=4294967296',
        '-sMODULARIZE=1', '-sINVOKE_RUN=0', '-sFORCE_FILESYSTEM=1',
        '-sEXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS"]',
        '-sENVIRONMENT=web,worker', '-lworkerfs.js', '-lproxyfs.js', '-o', str(dest / 'probe.js'))
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Built 3 aligner programs with maximum memory {MAXIMUM_MEMORY} bytes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, default=ROOT / '.cache/browser-aligners')
    parser.add_argument('--output', type=Path, default=ROOT / 'site/aligners')
    args = parser.parse_args()
    build(args.work_dir.resolve(), args.output.resolve())
