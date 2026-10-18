import io
import json
import pytest

from msap_progress import ProgressReporter, heartbeat, progress_context


def test_progress_callback_and_throttled_heartbeat():
    now = [10.0]
    events = []
    stream = io.StringIO()
    reporter = ProgressReporter('MSAP_batch', 4, 'loci', callback=events.append,
                                clock=lambda: now[0], stream=stream)
    reporter.update('Aligning')
    now[0] += 2
    reporter.tick()
    assert len(events) == 1
    now[0] += 4
    reporter.tick()
    assert len(events) == 2
    reporter.update('gene1.fa', completed=1, status='failed')
    event = json.loads(events[-1])
    assert event['completed'] == 1 and event['total'] == 4
    assert event['counts'] == {'failed': 1}
    assert event['eta_seconds'] == 18
    assert '1/4 loci' in stream.getvalue()
    assert 'failed=1' in stream.getvalue()


def test_stage_progress_has_no_time_percentage_or_eta():
    events = []
    reporter = ProgressReporter('PhyloPrep', 3, 'stages', callback=events.append, stream=io.StringIO())
    with progress_context(reporter):
        reporter.update('raw QC', completed=1)
        heartbeat()
    event = json.loads(events[-1])
    assert event['unit'] == 'stages'
    assert event['eta_seconds'] is None
    count = len(events)
    heartbeat()
    assert len(events) == count


@pytest.mark.parametrize('no_report', [False, True])
def test_browser_python_progress_streams_before_result(tmp_path, monkeypatch, no_report):
    import importlib.util
    from pathlib import Path
    import phyloprep
    from msap_progress import report_phase
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('progress_browser_runner', root / 'site/browser_runner.py')
    browser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(browser)
    source = tmp_path / 'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    events = []
    def variant(args, results, stage, variant):
        report_phase(args, 'Building test matrix')
        assert events and json.loads(events[-1])['label'] == 'Building test matrix'
        target = stage / 'prot' / variant
        target.mkdir(parents=True)
        (target / 'supermatrix.fasta').write_text('>a\nMA\n>b\nMA\n')
    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    result = json.loads(browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=True, no_report=no_report), str(tmp_path / 'out'), str(root), progress_callback=events.append))
    assert result['ok']
    final = json.loads(events[-1])
    assert final['completed'] == final['total'] == (2 if no_report else 3)
    assert json.loads(events[-1])['label'] == 'Complete'


def test_batch_progress_is_visible_before_alignment_finishes(tmp_path, fake_aligner):
    import subprocess
    import sys
    from pathlib import Path
    executable, env = fake_aligner
    executable.write_text(executable.read_text().replace('import os, sys, time, signal', 'import os, sys, time, signal\ntime.sleep(1.5)'))
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    root = Path(__file__).resolve().parents[1]
    process = subprocess.Popen([sys.executable, str(root / 'MSAP_batch.py'), '-i', str(source),
                                '-o', str(tmp_path / 'out'), '--notrim', '-t', '1'],
                               env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        first = process.stderr.readline()
        assert '0/1 loci' in first
        assert process.poll() is None
        stdout, stderr = process.communicate(timeout=20)
        assert process.returncode == 0, stdout + stderr
        assert '1/1 loci' in stderr and 'completed=1' in stderr
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate()


def test_batch_failed_tasks_are_counted_without_success_label(tmp_path, fake_aligner):
    import subprocess
    import sys
    from pathlib import Path
    _, env = fake_aligner
    env['MSAP_TEST_FAIL'] = '1'
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([sys.executable, str(root / 'MSAP_batch.py'), '-i', str(source),
                             '-o', str(tmp_path / 'out'), '--notrim', '-t', '1'],
                            env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'failed=1' in result.stderr
    assert 'Finished with failures' in result.stderr
    assert '| Complete |' not in result.stderr


def test_throttling_retains_all_task_counts():
    now = [0.0]
    events = []
    reporter = ProgressReporter('test', 3, 'loci', clock=lambda: now[0], callback=events.append, stream=io.StringIO())
    reporter.update('start')
    reporter.update('first', completed=1, status='included', force=False)
    reporter.update('second', completed=2, status='included', force=False)
    assert len(events) == 1
    reporter.update('third', completed=3, status='skipped')
    assert json.loads(events[-1])['counts'] == {'included': 2, 'skipped': 1}


def test_concatenation_streams_nested_phases_and_failure(tmp_path):
    from get_supergenes import get_supergenes
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nAAA\n>b\nCCC\n')
    events = []
    parent = ProgressReporter('PhyloPrep', 4, 'stages', callback=events.append, stream=io.StringIO())
    parent.update('raw/nucl: concatenating loci')
    with progress_context(parent):
        get_supergenes([str(source)], str(tmp_path / 'matrix'))
    child = [json.loads(event) for event in events if json.loads(event)['scope'] == 'get_supergenes']
    assert child and all('raw/nucl: concatenating loci' in event['label'] for event in child)
    for label in ['Scanning loci', 'Concatenating loci', 'Writing supermatrix', 'Complete']:
        assert any(label in event['label'] for event in child)
    assert parent.completed == 0
    source.write_text('>a\nAAA\n>a\nCCC\n')
    with progress_context(parent), pytest.raises(ValueError, match='duplicate'):
        get_supergenes([str(source)], str(tmp_path / 'bad'))
    assert 'Failed' in json.loads(events[-1])['label']
