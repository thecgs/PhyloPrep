"""Batch pathlists reuse confirmed task outcomes without another validation pass."""
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import MSAP_batch as batch


def test_pathlist_writer_does_not_read_or_validate_files(tmp_path, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError('pathlist generation must not revalidate task files')

    for name in ('_read_json', '_signature', '_digest', '_manifest_valid'):
        monkeypatch.setattr(batch, name, unexpected)
    # No file access is needed: task results are the authority for this run.
    outputs = [tmp_path / 'b.mafft.nucl.aln.fasta', tmp_path / 'a.mafft.nucl.aln.fasta']
    batch._write_path_lists(tmp_path, SimpleNamespace(seqtype='nucl', notrim=True), outputs * 2)
    pathlist = tmp_path / 'path-lists/all.nucl.aln.fasta.pathlist'
    assert pathlist.read_text().splitlines() == [str(p.resolve()) for p in sorted(outputs)]


def test_batch_reuses_fresh_and_resumed_results(tmp_path, fake_aligner, monkeypatch):
    _, env = fake_aligner
    for name, value in env.items():
        monkeypatch.setenv(name, value)
    monkeypatch.setattr(batch.signal, 'signal', lambda *args: None)
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nACGT\n>b\nACGA\n')
    root = tmp_path / 'results'
    argv = ['-i', str(source), '-o', str(root), '-st', 'nucl', '--notrim', '-t', '1']
    signatures = Counter()
    original = batch._signature

    def signature(path):
        signatures[Path(path)] += 1
        return original(path)

    monkeypatch.setattr(batch, '_signature', signature)
    assert batch.main(argv) == 0
    # Protect the live alignment with before/after checks, but no third read for lists.
    assert signatures[source] == 2
    pathlist = root / 'path-lists/all.nucl.aln.fasta.pathlist'
    expected = str(root / 'gene.mafft.nucl.aln.fasta') + '\n'
    assert pathlist.read_text() == expected
    signatures.clear()
    assert batch.main(argv) == 0
    assert signatures[source] == 1
    assert pathlist.read_text() == expected
    assert Path(env['MSAP_TEST_LOG']).read_text().splitlines() == ['run']

    # A corrupt resumed output must still be detected and rebuilt.
    (root / 'gene.mafft.nucl.aln.fasta').write_text('>a\nTTTT\n>b\nTTTA\n')
    assert batch.main(argv) == 0
    assert Path(env['MSAP_TEST_LOG']).read_text().splitlines() == ['run', 'run']
    monkeypatch.setenv('MSAP_TEST_FAIL', '1')
    assert batch.main(argv + ['--no-resume']) == 1
    assert pathlist.read_text() == ''
    assert (root / 'gene.mafft.nucl.aln.fasta').exists()


def test_cancelled_batch_lists_only_confirmed_success(tmp_path, monkeypatch):
    monkeypatch.setattr(batch, 'STOP', False)
    monkeypatch.setattr(batch.signal, 'signal', lambda *args: None)
    monkeypatch.setattr(batch, '_configuration', lambda args: {})
    sources = [tmp_path / f'{name}.fa' for name in ('complete', 'cancelled', 'skipped', 'failed')]
    for source in sources:
        source.write_text('>a\nACGT\n')
    root = tmp_path / 'results'
    output = root / 'complete.mafft.nucl.aln.fasta'

    def run_one(args, source, root):
        if source.stem == 'complete':
            output.write_text('>a\nACGT\n')
            return batch.TaskResult('complete', 'completed', (output,))
        # Stale public files for other tasks must not make it into the list.
        (root / f'{source.stem}.mafft.nucl.aln.fasta').write_text('>old\nACGT\n')
        if source.stem == 'failed':
            raise RuntimeError('failed publication')
        return batch.TaskResult(source.stem, 'cancelled before start' if source.stem == 'cancelled'
                                else 'skipped (no coding residues)')

    def interrupt_wait(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(batch, '_run_one', run_one)
    # Interrupt before the main loop observes results; executor shutdown joins
    # active jobs, and the final collection still retains confirmed successes.
    monkeypatch.setattr(batch.concurrent.futures, 'wait', interrupt_wait)
    assert batch.main(['-i', *map(str, sources), '-o', str(root), '-st', 'nucl',
                       '--notrim', '-j', '4', '-t', '1']) == 130
    assert (root / 'path-lists/all.nucl.aln.fasta.pathlist').read_text() == str(output) + '\n'
