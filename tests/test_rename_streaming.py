import tracemalloc

import pytest
from rename_taxa import rename_fasta


def test_rename_does_not_retain_all_sequence_records(tmp_path):
    source, target = tmp_path / 'locus.fasta', tmp_path / 'renamed.fasta'
    sequence = 'A' * 250000
    with source.open('w') as handle:
        for index in range(50):
            handle.write(f'>gene_{index}\n{sequence}\n')
    # Warm the parser/writer before measuring just this operation.
    warm = tmp_path / 'warm.fasta'
    warm.write_text('>gene\nAAA\n')
    rename_fasta(warm, tmp_path / 'warm-renamed.fasta')
    tracemalloc.start()
    try:
        rename_fasta(source, target)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 4 * 1024 * 1024, peak
    assert target.stat().st_size > 12_500_000


@pytest.mark.parametrize('content', ['', '>gene\nAA\n>gene\nCC\n'])
def test_streaming_failure_keeps_existing_output(tmp_path, content):
    source, target = tmp_path / 'locus.fasta', tmp_path / 'renamed.fasta'
    source.write_text(content)
    target.write_text('previous result\n')
    with pytest.raises(ValueError):
        rename_fasta(source, target)
    assert target.read_text() == 'previous result\n'
    assert source.read_text() == content


def test_mapping_audit_does_not_retain_all_report_rows(tmp_path):
    from phyloprep import annotate_change_report
    report = tmp_path / 'sequence_changes.tsv'
    with report.open('w') as handle:
        handle.write('taxa_id\tevent\toriginal\n')
        for _ in range(30000):
            handle.write('gene\tnucleotide_ambiguity\tR\n')
    tracemalloc.start()
    try:
        annotate_change_report(report, {'gene': 'species'})
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 2 * 1024 * 1024, peak
    with report.open() as handle:
        assert handle.readline().strip().endswith('mapped_taxa_id')
        assert all(line.strip().endswith('\tspecies') for line in handle)
