import csv
from types import SimpleNamespace

import phyloprep


def test_internal_stop_rejections_suggest_verifying_genetic_code(tmp_path):
    prefix = tmp_path / 'alignment'
    with prefix.with_suffix('.details.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['Issues'])
        writer.writerows([['internal_stop;N_ratio'], ['internal_stop'], ['short_alignment']])
    error = phyloprep.no_passing_paired_qc_error(SimpleNamespace(genetic_code=1), 'raw', {'codon': {'prefix': prefix}})
    assert isinstance(error, phyloprep.NoAlignmentsPassedQCError)
    assert '2 of 3 codon loci' in str(error)
    assert '-g 1' in str(error)
    assert 'table 2' in str(error)
    assert 'original CDS' in str(error)


def test_pis_only_rejection_does_not_suggest_genetic_code(tmp_path):
    prefix = tmp_path / 'alignment'
    prefix.with_suffix('.details.tsv').write_text('Issues\n""\n')
    error = phyloprep.no_passing_paired_qc_error(SimpleNamespace(genetic_code=1), 'trimmed', {'codon': {'prefix': prefix}})
    assert 'trimmed' in str(error)
    assert 'Inspect alignment.details.tsv' in str(error)
    assert 'Check that -g' not in str(error)


def test_wrong_table_returns_actionable_qc_error(tmp_path, fake_aligner):
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    _, env = fake_aligner
    source = tmp_path / 'mito.fa'
    source.write_text('>a\nGCTTGAGCT\n>b\nGCCTGAGCC\n')
    for table in (1, 2):
        output = tmp_path / f'code-{table}'
        result = subprocess.run([sys.executable, str(root / 'phyloprep.py'), '-i', str(source),
                                 '-o', str(output), '-g', str(table), '--notrim', '-t', '1'],
                                env=env, capture_output=True, text=True)
        if table == 1:
            assert result.returncode != 0
            assert '1 of 1 codon loci' in result.stderr
            assert 'Current genetic-code table: -g 1' in result.stderr
            assert 'Traceback' not in result.stderr
            assert (output / 'qc/codon/raw/alignment.details.tsv').exists()
        else:
            assert result.returncode == 0, result.stdout + result.stderr
            assert (output / 'matrices/codon/raw/supermatrix.fasta').exists()


def test_browser_qc_rejection_is_readable_and_keeps_reports(tmp_path, monkeypatch):
    import importlib.util
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('hint_browser_runner', root / 'site/browser_runner.py')
    browser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(browser)
    source = tmp_path / 'source.fa'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    output = tmp_path / 'out'
    def reject(args, results, stage, variant):
        report = results / 'qc/codon/raw/alignment.details.tsv'
        report.parent.mkdir(parents=True)
        report.write_text('Issues\ninternal_stop\n')
        raise phyloprep.NoAlignmentsPassedQCError('Check that -g matches the source CDS.')
    monkeypatch.setattr(phyloprep, 'prepare_variant', reject)
    result = json.loads(browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='codon', notrim=True), str(output), str(root)))
    assert not result['ok'] and result['recoverable']
    assert 'Check that -g' in result['user_error']
    assert 'Traceback' not in result['stderr']
    assert result['files']
    assert not (output / 'matrices').exists()
