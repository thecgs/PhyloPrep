"""Automatic reporting runs after publication and cannot invalidate analysis."""
import importlib.util
import json
from pathlib import Path
import re

import pytest

import phyloprep
import phyloprep_report

ROOT = Path(__file__).resolve().parents[1]


def payload(path):
    return json.loads(re.search(r'id="report-data">(.*?)</script>', path.read_text(), re.S).group(1))


def cli_run(tmp_path, monkeypatch, fake_aligner, *options):
    _, env = fake_aligner
    for name in ('PATH', 'MSAP_TEST_LOG', 'MSAP_TEST_PID'):
        monkeypatch.setenv(name, env[name])
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nACGT\n>b\nACGA\n')
    out = tmp_path / 'results'
    status = phyloprep.main(['-i', str(source), '-o', str(out), '-st', 'nucl',
                             '--notrim', '-t', '1', *options])
    return status, out


def test_cli_default_report_follows_published_matrices(tmp_path, monkeypatch, fake_aligner):
    original = phyloprep_report.write_report
    seen = []

    def inspect(root, **kwargs):
        assert (root / 'matrices/nucl/raw/supermatrix.fasta').is_file()
        assert not list(root.glob('.matrices-*'))
        seen.append(root)
        return original(root, **kwargs)

    monkeypatch.setattr(phyloprep_report, 'write_report', inspect)
    status, out = cli_run(tmp_path, monkeypatch, fake_aligner)
    assert status == 0 and seen == [out]
    data = payload(out / 'phyloprep-report.html')
    assert data['overview']['input_loci'] == 1
    assert data['provenance']['parameters']['seqtype'] == 'nucl'
    assert data['matrices'][0]['length'] == 4


@pytest.mark.parametrize('error', [ValueError('broken report table'), OSError('report disk full'),
                                  MemoryError('report heap exhausted'), ImportError('report package unavailable')])
def test_cli_report_failure_is_nonfatal_and_removes_old_html(tmp_path, monkeypatch, fake_aligner, capsys, error):
    out = tmp_path / 'results'
    out.mkdir()
    (out / 'phyloprep-report.html').write_text('OLD RUN')

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(phyloprep_report, 'write_report', fail)
    status, out = cli_run(tmp_path, monkeypatch, fake_aligner)
    assert status == 0
    assert (out / 'matrices/nucl/raw/supermatrix.fasta').read_text() == '>a\nACGT\n>b\nACGA\n'
    assert (out / 'qc/nucl/raw/alignment.details.tsv').exists()
    assert not (out / 'phyloprep-report.html').exists()
    message = capsys.readouterr().err
    assert type(error).__name__ in message and str(error) in message
    assert 'Completed analysis outputs have been retained' in message
    assert 'Retry: python phyloprep_report.py -i' in message


def test_cli_opt_out_does_not_call_report_and_removes_previous_html(tmp_path, monkeypatch, fake_aligner):
    out = tmp_path / 'results'
    out.mkdir()
    (out / 'phyloprep-report.html').write_text('OLD RUN')
    monkeypatch.setattr(phyloprep_report, 'write_report', lambda *a, **kw: pytest.fail('report should be skipped'))
    status, out = cli_run(tmp_path, monkeypatch, fake_aligner, '--no-report')
    assert status == 0
    assert (out / 'matrices/nucl/raw/supermatrix.fasta').exists()
    assert not (out / 'phyloprep-report.html').exists()


def test_failed_analysis_never_generates_or_invalidates_report(tmp_path, monkeypatch, fake_aligner):
    out = tmp_path / 'results'
    out.mkdir()
    (out / 'phyloprep-report.html').write_text('OLD PUBLISHED RUN')
    def fail(*args):
        raise ValueError('QC failed before publication')
    monkeypatch.setattr(phyloprep, 'prepare_variant', fail)
    monkeypatch.setattr(phyloprep_report, 'write_report', lambda *a, **kw: pytest.fail('analysis did not complete'))
    status, out = cli_run(tmp_path, monkeypatch, fake_aligner)
    assert status == 2
    assert (out / 'phyloprep-report.html').read_text() == 'OLD PUBLISHED RUN'


def browser_run(tmp_path, options):
    spec = importlib.util.spec_from_file_location('report_browser', ROOT / 'site/browser_runner.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nACGT\n>b\nACGA\n')
    alignment = tmp_path / 'gene.mafft.nucl.aln.fasta'
    alignment.write_text(source.read_text())
    out = tmp_path / 'results'
    result = json.loads(module.run_phyloprep_finalize([
        dict(role='source', name=source.name, path=str(source)),
        dict(role='batch-output', name=alignment.name, path=str(alignment))],
        dict(seqtype='nucl', notrim=True, **options), str(out), str(ROOT), manifest_only=True))
    return result, out


def test_browser_default_report_is_downloadable_and_contains_current_metadata(tmp_path):
    result, out = browser_run(tmp_path, dict(report_context={'reused_loci': 1}))
    assert result['ok'] and result['report']['status'] == 'generated'
    assert not result['warnings']
    assert any(file['name'] == 'phyloprep-report.html' and 'content' not in file for file in result['files'])
    data = payload(out / 'phyloprep-report.html')
    assert data['overview']['input_loci'] == 1 and data['overview']['reused_loci'] == 1
    assert data['provenance']['parameters']['seqtype'] == 'nucl'
    assert 'report_context' not in data['provenance']['parameters']


@pytest.mark.parametrize('error', [OSError('report storage full'), MemoryError('report memory full')])
def test_browser_report_failure_preserves_success_and_downloads(tmp_path, monkeypatch, error):
    def fail(*args, **kwargs):
        raise error
    monkeypatch.setattr(phyloprep_report, 'write_report', fail)
    result, out = browser_run(tmp_path, {})
    assert result['ok'] and result['report']['status'] == 'failed'
    assert str(error) in result['warnings'][0] and type(error).__name__ in result['stderr']
    assert any(file['name'] == 'matrices/nucl/raw/supermatrix.fasta' for file in result['files'])
    assert not any(file['name'] == 'phyloprep-report.html' for file in result['files'])
    assert (out / 'matrices/nucl/raw/supermatrix.nex').exists()


def test_browser_report_package_failure_and_opt_out(tmp_path):
    result, out = browser_run(tmp_path, {'report_dependency_error': 'SQLite download failed: network offline'})
    assert result['ok'] and result['report']['status'] == 'failed'
    assert 'SQLite download failed: network offline' in result['warnings'][0]
    other = tmp_path / 'skipped'
    other.mkdir()
    result, out = browser_run(other, {'no_report': True, 'report_dependency_error': 'must be ignored'})
    assert result['ok'] and result['report']['status'] == 'skipped'
    assert not result['warnings']
    assert not (out / 'phyloprep-report.html').exists()


def test_auto_metadata_does_not_reuse_previous_run_manifest(tmp_path):
    root = tmp_path / 'output'
    root.mkdir()
    (root / 'reports').mkdir()
    (root / 'phyloprep-run.json').write_text('{"parameters":{"genetic_code":99}}')
    phyloprep_report.write_report(root, run_metadata={'parameters': {'genetic_code': 2}})
    assert payload(root / 'phyloprep-report.html')['provenance']['parameters']['genetic_code'] == 2


def test_report_current_scope_excludes_saved_other_modes_and_trimmed_branches(tmp_path):
    for group, sequence in [('nucl/raw', 'ACGT'), ('nucl/trimmed', 'ACGT'), ('prot/raw', 'MACD'),
                            ('paired/prot/raw', 'MACD'), ('pseudogene/codon/raw', 'GCTGCT')]:
        folder = tmp_path / 'matrices' / group
        folder.mkdir(parents=True)
        (folder / 'supermatrix.fasta').write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    data = phyloprep_report.build_report(tmp_path, run_metadata={'parameters': {'seqtype': 'nucl', 'notrim': True}})
    assert [matrix['matrix'] for matrix in data['matrices']] == ['nucl/raw/supermatrix.fasta']
    assert (tmp_path / 'matrices/prot/raw/supermatrix.fasta').exists()
