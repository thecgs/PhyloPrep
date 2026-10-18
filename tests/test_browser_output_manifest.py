import importlib.util
import json
import pytest
from pathlib import Path

import phyloprep

ROOT = Path(__file__).resolve().parents[1]


def browser_module():
    spec = importlib.util.spec_from_file_location('output_manifest_browser', ROOT / 'site/browser_runner.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_manifest_does_not_read_output_contents(tmp_path, monkeypatch):
    browser = browser_module()
    (tmp_path / 'matrix.fasta').write_text('>a\nAC\n')
    hidden = tmp_path / '.staging'
    hidden.mkdir()
    (hidden / 'private.fasta').write_text('private')
    monkeypatch.setattr(Path, 'read_text', lambda *args, **kwargs: (_ for _ in ()).throw(MemoryError('contents must not be read')))
    files = browser._visible_output_files(tmp_path, manifest_only=True)
    assert files == [dict(name='matrix.fasta', path=str(tmp_path / 'matrix.fasta'), size=6)]


def test_finalize_serializes_metadata_not_large_matrix(tmp_path, monkeypatch):
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    def variant(args, results, stage, variant):
        target = stage / 'prot' / variant
        target.mkdir(parents=True)
        (target / 'supermatrix.fasta').write_bytes(b'>a\n' + b'A' * (8 * 1024 * 1024) + b'\n')
    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    raw = browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=True), str(tmp_path / 'out'), str(ROOT), manifest_only=True)
    assert len(raw) < 10000
    result = json.loads(raw)
    assert result['ok']
    matrix = next(file for file in result['files'] if file['name'].endswith('supermatrix.fasta'))
    assert matrix['size'] > 8 * 1024 * 1024
    assert 'content' not in matrix
    assert Path(matrix['path']).exists()


def test_bounded_logs_retain_diagnostics_without_full_history():
    browser = browser_module()
    log = browser.BoundedLog(limit=100)
    for _ in range(10000):
        log.write('old message\n')
    log.write('last diagnostic')
    assert log.length <= 100
    assert log.getvalue().startswith('[Earlier log output truncated]\n')
    assert log.getvalue().endswith('last diagnostic')
    log = browser.BoundedLog()
    for _ in range(10000):
        log.write('x')
    assert len(log.parts) <= 1024
    log.write('A' * 1000000)
    assert log.length == 65536


def test_finalize_memory_failure_uses_prebuilt_response(tmp_path, monkeypatch):
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    def fail(*args, **kwargs):
        raise MemoryError('simulated memory exhaustion')
    monkeypatch.setattr(phyloprep, 'prepare_variant', fail)
    monkeypatch.setattr(browser.BoundedLog, 'getvalue', fail)
    monkeypatch.setattr(browser.traceback, 'print_exc', fail)
    monkeypatch.setattr(browser.json, 'dumps', fail)
    raw = browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=True), str(tmp_path / 'out'), str(ROOT), manifest_only=True)
    assert raw is browser._MEMORY_ERROR_RESULT
    result = json.loads(raw)
    assert not result['ok']
    assert 'Python/Pyodide' in result['error']


def test_incremental_finalize_exports_qc_without_publishing_staged_matrices(tmp_path, monkeypatch):
    browser = browser_module()
    import phyloprep
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    seen = []

    def variant(args, root, stage, name):
        assert source.exists()
        if name == 'trimmed':
            assert (root / 'qc/raw/alignment.details.tsv').exists()
        report = root / 'qc' / name / 'alignment.details.tsv'
        report.parent.mkdir(parents=True)
        report.write_text('locus\tpass\ngene\tyes\n')
        matrix = stage / 'prot' / name / 'supermatrix.fasta'
        matrix.parent.mkdir(parents=True)
        matrix.write_text('>a\nMA\n>b\nMA\n')
        seen.append(name)

    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT),
        manifest_only=True, incremental=True)
    try:
        for expected in ('raw', 'trimmed'):
            checkpoint = json.loads(browser.advance_phyloprep_finalize(steps))
            assert checkpoint['pending']
            assert seen[-1] == expected
            assert all(file['name'].startswith(('qc/', 'reports/')) for file in checkpoint['files'])
            assert all(Path(file['path']).exists() for file in checkpoint['files'])
            assert not (tmp_path / 'out/matrices').exists()
        completed = json.loads(browser.advance_phyloprep_finalize(steps))
        assert completed['ok']
        assert any(file['name'].startswith('matrices/') for file in completed['files'])
        assert source.exists()
    finally:
        steps.close()


def test_closing_incremental_finalize_restores_engine_and_streams(tmp_path, monkeypatch):
    import sys
    import phyloprep
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    original_run = phyloprep.run
    original_stdout, original_stderr = sys.stdout, sys.stderr
    monkeypatch.setattr(phyloprep, 'prepare_variant', lambda *args: None)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT),
        manifest_only=True, incremental=True)
    try:
        next(steps)
    finally:
        steps.close()
    assert phyloprep.run is original_run
    assert sys.stdout is original_stdout
    assert sys.stderr is original_stderr
    assert not list((tmp_path / 'out').glob('.matrices-*'))
    assert not (tmp_path / 'out/matrices').exists()


def test_portable_export_lists_resolve_after_moving_without_changing_engine_files(tmp_path):
    import shutil
    from msap_io import expand_input_paths
    browser = browser_module()
    root = tmp_path / 'run'
    alignment = root / 'alignments/gene.mafft.codon.aln.fasta'
    alignment.parent.mkdir(parents=True)
    alignment.write_text('>a\nATG\n>b\nATG\n')
    for name in ('alignments/path-lists/all.codon.aln.fasta.pathlist',
                 'qc/codon/raw/alignment.pass.tsv', 'qc/paired/raw/codon.retained.tsv'):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(str(alignment) + '\n')
    files = browser._export_output_files(root, browser._visible_output_files(root, manifest_only=True))
    saved = tmp_path / 'saved'
    for item in files:
        target = saved / item['name']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(item['path'], target)
    assert (root / 'qc/codon/raw/alignment.pass.tsv').read_text() == str(alignment) + '\n'
    shutil.rmtree(root)
    for item in files:
        if item['name'].endswith(('.pathlist', '.pass.tsv', '.retained.tsv')):
            assert expand_input_paths([saved / item['name']]) == [str(saved / 'alignments/gene.mafft.codon.aln.fasta')]


def test_independent_tool_manifest_avoids_serializing_large_output(tmp_path):
    browser = browser_module()
    engine = tmp_path / 'engine'
    engine.mkdir()
    (engine / 'large.py').write_text("from pathlib import Path\nPath('matrix.fasta').write_bytes(b'A' * (8 * 1024 * 1024))\n")
    result = browser.run_original_module('large.py', [], str(tmp_path / 'out'), str(engine), manifest_only=True)
    assert len(result) < 2000
    item = json.loads(result)['files'][0]
    assert item['size'] == 8 * 1024 * 1024
    assert 'content' not in item
    assert len(item['sha256']) == 64


def test_export_fingerprints_change_when_report_content_changes(tmp_path):
    browser = browser_module()
    path = tmp_path / 'report.tsv'
    path.write_text('first')
    first = browser._export_output_files(tmp_path, browser._visible_output_files(tmp_path, manifest_only=True))[0]['sha256']
    path.write_text('other')  # Same size must still be detected.
    second = browser._export_output_files(tmp_path, browser._visible_output_files(tmp_path, manifest_only=True))[0]['sha256']
    assert first != second


def test_input_only_export_list_uses_basename_without_bundling_source(tmp_path):
    browser = browser_module()
    source = tmp_path / 'original.fa'
    source.write_text('>a\nATG\n')
    root = tmp_path / 'out'
    root.mkdir()
    path = root / 'inputs.pathlist'
    path.write_text(str(source) + '\n')
    files = browser._export_output_files(root, browser._visible_output_files(root, manifest_only=True))
    assert Path(files[0]['path']).read_text() == 'original.fa\n'
    assert path.read_text() == str(source) + '\n'
    assert len(files) == 1


def test_explicit_custom_pathlist_extension_is_portable(tmp_path):
    browser = browser_module()
    engine = tmp_path / 'engine'
    engine.mkdir()
    (engine / 'custom.py').write_text("from pathlib import Path\nimport sys\np=Path('gene.fasta')\np.write_text('>a\\nATG\\n')\nPath(sys.argv[2]).write_text(str(p.resolve()) + '\\n')\n")
    root = tmp_path / 'out'
    result = json.loads(browser.run_original_module('custom.py', ['--pathlist-output', 'files.txt'], str(root), str(engine), manifest_only=True))
    assert result['ok']
    item = next(file for file in result['files'] if file['name'] == 'files.txt')
    assert Path(item['path']).read_text() == 'gene.fasta\n'
    assert (root / 'files.txt').read_text() == str(root / 'gene.fasta') + '\n'


@pytest.mark.parametrize("notrim", [False, True])
def test_incremental_finalize_resumes_in_different_contexts(tmp_path, monkeypatch, notrim):
    from contextvars import Context
    import sys
    import msap_progress
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    original_run = phyloprep.run
    original_streams = (sys.stdout, sys.stderr)
    events = []

    def variant(args, root, stage, name):
        assert msap_progress._CURRENT.get() is args._progress
        assert phyloprep.run is not original_run
        report = root / 'qc' / name / 'alignment.details.tsv'
        report.parent.mkdir(parents=True)
        report.write_text('locus\tpass\ngene\tyes\n')
        matrix = stage / 'prot' / name / 'supermatrix.fasta'
        matrix.parent.mkdir(parents=True)
        matrix.write_text('>a\nMA\n>b\nMA\n')

    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=notrim), str(tmp_path / 'out'), str(ROOT),
        progress_callback=lambda text: events.append(json.loads(text)), manifest_only=True, incremental=True)
    try:
        for pending in ((True, False) if notrim else (True, True, False)):
            context = Context()
            marker = object()
            context.run(msap_progress._CURRENT.set, marker)
            result = json.loads(context.run(browser.advance_phyloprep_finalize, steps))
            assert context.run(msap_progress._CURRENT.get) is marker
            assert (sys.stdout, sys.stderr) == original_streams
            assert phyloprep.run is original_run
            assert bool(result.get('pending')) is pending
        assert result['ok'], result
        matrices = [file for file in result['files'] if file['name'].startswith('matrices/')]
        assert matrices and all(Path(file['path']).is_file() for file in matrices)
        assert events[-1]['completed'] == events[-1]['total'] == (3 if notrim else 4)
    finally:
        Context().run(steps.close)


def test_incremental_finalize_closes_in_another_context(tmp_path, monkeypatch):
    from contextvars import Context
    import sys
    import msap_progress
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    original_run = phyloprep.run
    original_streams = (sys.stdout, sys.stderr)
    monkeypatch.setattr(phyloprep, 'prepare_variant', lambda *args: None)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT), manifest_only=True, incremental=True)
    first, closing = Context(), Context()
    marker = object()
    first.run(msap_progress._CURRENT.set, marker)
    first.run(next, steps)
    assert first.run(msap_progress._CURRENT.get) is marker
    assert phyloprep.run is original_run
    assert (sys.stdout, sys.stderr) == original_streams
    closing.run(steps.close)
    assert not list((tmp_path / 'out').glob('.matrices-*'))
    assert not (tmp_path / 'out/matrices').exists()
    assert phyloprep.run is original_run
    assert (sys.stdout, sys.stderr) == original_streams


def test_incremental_finalize_failure_restores_calling_context(tmp_path, monkeypatch):
    from contextvars import Context
    import sys
    import msap_progress
    browser = browser_module()
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    original_run = phyloprep.run
    original_streams = (sys.stdout, sys.stderr)

    def fail(args, *others):
        assert msap_progress._CURRENT.get() is args._progress
        raise ValueError('simulated QC failure')

    monkeypatch.setattr(phyloprep, 'prepare_variant', fail)
    steps = browser._phyloprep_finalize_steps(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot'), str(tmp_path / 'out'), str(ROOT), manifest_only=True, incremental=True)
    context = Context()
    marker = object()
    context.run(msap_progress._CURRENT.set, marker)
    result = json.loads(context.run(browser.advance_phyloprep_finalize, steps))
    assert not result['ok']
    assert 'simulated QC failure' in result['error']
    assert context.run(msap_progress._CURRENT.get) is marker
    assert phyloprep.run is original_run
    assert (sys.stdout, sys.stderr) == original_streams
    assert not list((tmp_path / 'out').glob('.matrices-*'))


@pytest.mark.parametrize('fail_trim', [False, True])
def test_disk_backed_branches_release_original_and_renamed_after_use(tmp_path, monkeypatch, fail_trim):
    """Run real QC/concatenation, simulating successful checkpoint disk writes."""
    browser = browser_module()
    source = tmp_path / 'gene.fa'
    source.write_text('>gene_a\nATGGCT\n>gene_b\nATGGCT\n')
    mapping = tmp_path / 'map.tsv'
    mapping.write_text('gene_a\tsp_a\ngene_b\tsp_b\n')
    root = tmp_path / 'out'
    products = {
        'gene.mafft.codon.aln.fasta': '>gene_a\nATGGCT\n>gene_b\nATGGCT\n',
        'gene.mafft.prot.aln.fasta': '>gene_a\nMA\n>gene_b\nMA\n',
        'gene.mafft.codon.trimmed.aln.fasta': '>gene_a\nATGGCT\n>gene_b\nATGGCT\n',
        'gene.mafft.prot.trimmed.aln.fasta': '>gene_a\nMA\n>gene_b\nMA\n',
    }
    manifest = [dict(role='source', path=str(source), name=source.name),
                dict(role='mapping', path=str(mapping), name=mapping.name)]
    manifest += [dict(role='batch-output', path=str(tmp_path / name), name=name, on_disk=True)
                 for name in products]
    steps = browser._phyloprep_finalize_steps(manifest, dict(seqtype='codon'), str(root),
        str(ROOT), manifest_only=True, incremental=True, disk_backed=True)
    saved = {}
    if fail_trim:
        original_prepare = phyloprep.prepare_variant
        def prepare(args, root, stage, variant):
            if variant == 'trimmed':
                raise phyloprep.NoAlignmentsPassedQCError('trimmed QC rejected all loci')
            return original_prepare(args, root, stage, variant)
        monkeypatch.setattr(phyloprep, 'prepare_variant', prepare)
    try:
        for variant in ('raw', 'trimmed'):
            loading = json.loads(browser.advance_phyloprep_finalize(steps))
            assert loading['pending'] and loading['load_variant'] == variant
            for name, content in products.items():
                if ('.trimmed.aln.fasta' in name) == (variant == 'trimmed'):
                    (root / 'alignments' / name).write_text(content)
            checkpoint = json.loads(browser.advance_phyloprep_finalize(steps))
            early = []
            while checkpoint.get('pending') and checkpoint.get('files') and all(
                    file['name'].startswith('renamed/') and not file.get('release')
                    for file in checkpoint['files']):
                assert len(checkpoint['files']) == 1
                item = checkpoint['files'][0]
                assert Path(item['path']).exists()
                assert not (root / 'matrices').exists()
                saved[item['name']] = Path(item['path']).read_bytes()
                early.append(item['name'])
                checkpoint = json.loads(browser.advance_phyloprep_finalize(steps))
            assert len(early) == 2  # Each renamed locus is exported before QC/concatenation.
            if fail_trim and variant == 'trimmed':
                assert checkpoint['recoverable']
                names = {file['name'] for file in checkpoint['files']}
                assert 'renamed/codon/raw/gene.mafft.codon.aln.fasta' in names
                assert not (root / 'matrices').exists()
                return
            assert checkpoint.get('pending'), checkpoint
            releasable = [file for file in checkpoint['files'] if file.get('release')]
            assert len(releasable) == 4  # original and renamed codon/protein
            assert all(('.trimmed.aln.fasta' in file['name']) == (variant == 'trimmed')
                       for file in releasable)
            for file in checkpoint['files']:
                saved[file['name']] = Path(file['path']).read_bytes()
            for file in releasable:
                Path(file['path']).unlink()
            assert not list((root / 'renamed').rglob('*.aln.fasta'))
        result = json.loads(browser.advance_phyloprep_finalize(steps))
        assert result['ok'], result
        files = {file['name']: file for file in result['files']}
        assert all('alignments/' + name in files for name in products)
        assert len([name for name in files if name.startswith('renamed/')]) == 4
        assert (root / 'matrices/codon/raw/supermatrix.fasta').read_text() == '>sp_a\nATGGCT\n>sp_b\nATGGCT\n'
        assert (root / 'matrices/paired/codon/raw/supermatrix.fasta').read_text() == '>sp_a\nATGGCT\n>sp_b\nATGGCT\n'
        for name in ('qc/codon/raw/alignment.pass.tsv', 'qc/codon/trimmed/alignment.pass.tsv'):
            assert str(root).encode() not in Path(files[name]['path']).read_bytes()
        assert source.read_text().startswith('>gene_a')
    finally:
        steps.close()


@pytest.mark.parametrize('reuse_protein', [True, False])
def test_rename_checkpoint_reuses_saved_targets_before_qc(tmp_path, monkeypatch, reuse_protein):
    """Cached targets skip rename_fasta; missing entries alone are recomputed."""
    import hashlib
    browser = browser_module()
    source = tmp_path / 'gene.fa'
    source.write_text('>g1\nATGGCT\n>g2\nATGGCT\n')
    mapping = tmp_path / 'mapping.tsv'
    mapping.write_text('g1\tsp_a\ng2\tsp_b\n')
    root = tmp_path / 'out'
    manifest = [dict(role='source', name=source.name, path=str(source)),
                dict(role='mapping', name=mapping.name, path=str(mapping))]
    products = {'gene.mafft.codon.aln.fasta': '>g1\nATGGCT\n>g2\nATGGCT\n',
                'gene.mafft.prot.aln.fasta': '>g1\nMA\n>g2\nMA\n'}
    for name, content in products.items():
        manifest.append(dict(role='batch-output', name=name, path=str(tmp_path / name), on_disk=True))
        if '.prot.' in name and not reuse_protein:
            continue
        kind = 'codon' if '.codon.' in name else 'prot'
        relative = f'renamed/{kind}/raw/{name}'
        text = content.replace('g1', 'sp_a').replace('g2', 'sp_b')
        manifest.append(dict(role='renamed-output', name='cached-'+name, rename_path=relative,
                             path=str(tmp_path / ('cached-'+name)), on_disk=True,
                             sha256=hashlib.sha256(text.encode()).hexdigest()))
    calls = []
    original = phyloprep.rename_fasta
    def rename(source, *args, **kwargs):
        calls.append(source.name)
        return original(source, *args, **kwargs)
    monkeypatch.setattr(phyloprep, 'rename_fasta', rename)
    steps = browser._phyloprep_finalize_steps(manifest,
        dict(seqtype='codon', notrim=True, no_report=True, rename_checkpoints=True),
        str(root), str(ROOT), manifest_only=True, incremental=True, disk_backed=True)
    first = json.loads(browser.advance_phyloprep_finalize(steps))
    assert first['load_variant'] == 'raw'
    for name, content in products.items():
        (root / 'alignments' / name).write_text(content)
    for item in manifest:
        if item['role'] == 'renamed-output':
            target = root / item['rename_path']
            target.parent.mkdir(parents=True, exist_ok=True)
            name = target.name
            target.write_text(products[name].replace('g1', 'sp_a').replace('g2', 'sp_b'))
    checkpoints = []
    while True:
        step = json.loads(browser.advance_phyloprep_finalize(steps))
        if not step.get('pending'):
            result = step
            break
        checkpoints.append(step)
        if step.get('rename_checkpoint'):
            assert not (root / 'matrices').exists(), 'Rename checkpoint must precede QC/matrices'
    assert result['ok'], result
    assert calls == ([] if reuse_protein else ['gene.mafft.prot.aln.fasta'])
    assert sum(step.get('rename_checkpoint') == 'raw' for step in checkpoints) == 1
    assert (root / 'matrices/codon/raw/supermatrix.fasta').read_text() == '>sp_a\nATGGCT\n>sp_b\nATGGCT\n'
    assert len([item for item in result['files'] if item['name'].startswith('renamed/')]) == 2
