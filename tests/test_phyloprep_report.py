"""Report semantics, portability, bounded reads and actual workflow integration."""
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tracemalloc

import pytest

import phyloprep_report as report

REPO = Path(__file__).resolve().parents[1]


def tsv(path, fields, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def matrix(root, group='codon/raw', sequences=None, parts=None, missing=(), included=None):
    directory = root / 'matrices' / group
    directory.mkdir(parents=True, exist_ok=True)
    sequences = sequences or {'a': 'ACG???', 'b': 'ACG---', 'c': 'ATGNNN', 'd': 'ATG???'}
    (directory / 'supermatrix.fasta').write_text(''.join(f'>{t}\n{s}\n' for t, s in sequences.items()))
    parts = parts or [('gene1', 1, 3), ('gene2', 4, 6)]
    (directory / 'supermatrix_part_iqtree.txt').write_text('#nexus\nbegin sets;\n' + ''.join(
        f' charset {name} = {start}-{end};\n' for name, start, end in parts) + 'end;\n')
    tsv(directory / 'supermatrix.missing.tsv', ['locus', 'taxon', 'start', 'end', 'length', 'fill_character'],
        [dict(locus=name, taxon=taxon, start=start, end=end, length=end - start + 1, fill_character='?')
         for name, taxon, start, end in missing])
    tsv(directory / 'supermatrix.report.tsv', ['gene', 'status', 'details'], included if included is not None else
        [dict(gene=f'/unavailable/{name}.mafft.codon.aln.fasta', status='included', details='') for name, _, _ in parts])
    return directory


def qc(root, branch='nucl/raw', issues=('gap_ratio;short_alignment', 'gap_ratio;short_alignment', 'N_ratio', '')):
    directory = root / 'qc' / branch
    tsv(directory / 'alignment.details.tsv', ['File', 'Seqtype', 'Parsimony_informative_sites', 'Issues'],
        [dict(File=f'/old/{i}.mafft.nucl.aln.fasta', Seqtype='nucl', Parsimony_informative_sites=0, Issues=reason)
         for i, reason in enumerate(issues)])
    for status, selected in [('pass', [i for i, r in enumerate(issues) if not r]),
                             ('fail', [i for i, r in enumerate(issues) if r])]:
        (directory / f'alignment.{status}.tsv').write_text(''.join(f'/old/{i}.mafft.nucl.aln.fasta\n' for i in selected))
    return directory


def test_overlapping_reasons_and_exclusive_combinations(tmp_path):
    qc(tmp_path)
    q = report.build_report(tmp_path)['qc'][0]
    assert q['failed'] == 3
    assert q['passed'] == 1
    assert sum(q['reason_hits'].values()) == 5
    assert sum(q['failure_combinations'].values()) == 3
    assert q['failure_combinations'] == {'gap_ratio + short_alignment': 2, 'N_ratio': 1}
    assert q['reasons_per_locus'] == {'2': 2, '1': 1}


def test_missing_qc_status_is_not_invented_pass(tmp_path):
    directory = qc(tmp_path)
    (directory / 'alignment.pass.tsv').unlink()
    (directory / 'alignment.fail.tsv').unlink()
    q = report.build_report(tmp_path)['qc'][0]
    assert q['passed'] == q['failed'] == 0
    assert q['not_evaluated'] == 4


def test_mapping_and_pathlist_are_not_counted_as_input_loci(tmp_path):
    manifest = dict(inputs=[{'name': 'gene.fa', 'input_key': 'inputs'},
                            {'name': 'taxa.tsv', 'input_key': 'mapping'},
                            {'name': 'genes.pathlist', 'input_key': 'inputs'}])
    assert report.processing_overview(tmp_path, manifest)['input_loci'] == 1
    (tmp_path / 'inputs.pathlist').write_text('/source/gene1.fa\n/source/gene2.fa\n')
    assert report.processing_overview(tmp_path, manifest)['input_loci'] == 2


def test_final_matrix_resolved_statistics_and_artificial_padding(tmp_path):
    matrix(tmp_path, missing=[('gene2', 'a', 4, 6)])
    result = report.build_report(tmp_path)['matrices'][0]
    assert result['length'] == 6 and result['complete_codons'] == 2
    assert result['pis'] == result['variable'] == 1
    assert result['constant'] == 2 and result['unresolved_only'] == 3
    assert result['gap_characters'] == 3
    assert result['missing_characters'] == 6
    assert result['ambiguous_characters'] == 3
    assert result['gc_ratio'] == pytest.approx(6 / 12)
    assert result['padding']['characters'] == 3  # Not all six literal '?' characters.
    assert result['padding']['blocks'] == 1
    coverage = {r['taxon']: r['observed_loci'] for r in result['padding']['coverage']}
    assert coverage == {'a': 1, 'b': 2, 'c': 2, 'd': 2}  # Existing gaps/N are observed.


def test_missing_table_absent_vs_known_zero(tmp_path):
    directory = matrix(tmp_path)
    assert report.build_report(tmp_path)['matrices'][0]['padding']['blocks'] == 0
    (directory / 'supermatrix.missing.tsv').unlink()
    assert report.build_report(tmp_path)['matrices'][0]['padding'] is None


def test_empty_file_is_not_a_known_zero_missing_table(tmp_path):
    directory = matrix(tmp_path)
    (directory / 'supermatrix.missing.tsv').write_text('')
    result = report.build_report(tmp_path)
    assert result['matrices'][0]['padding'] is None
    assert result['warnings']


def test_protein_states_have_no_gc_and_n_is_resolved(tmp_path):
    directory = matrix(tmp_path, 'prot/raw', sequences={'a': 'ANNX-?', 'b': 'ANNX-?'},
                       parts=[('gene1', 1, 3), ('gene2', 4, 6)])
    result = report.build_report(tmp_path)['matrices'][0]
    assert result['unit'] == 'aa'
    assert result['gc_ratio'] is result['complete_codons'] is None
    assert result['ambiguous_characters'] == 2
    assert result['constant'] == 3
    assert result['unresolved_only'] == 3


def test_actual_skipped_locus_report(tmp_path):
    included = [dict(gene=f'/old/{n}.mafft.codon.aln.fasta', status=s, details='')
                for n, s in [('gene1', 'included'), ('gene2', 'padded'), ('gene3', 'skipped')]]
    matrix(tmp_path, included=included)
    result = report.build_report(tmp_path)['matrices'][0]
    assert result['loci'] == 2
    assert result['skipped_after_qc'] == 1


def test_invalid_padding_is_unavailable(tmp_path):
    matrix(tmp_path, missing=[('gene2', 'not_a_taxon', 4, 6)])
    result = report.build_report(tmp_path)
    assert result['matrices'][0]['padding'] is None
    assert result['warnings']


def test_actual_inclusion_not_qc_pass_and_matched_trimming(tmp_path):
    qc(tmp_path, issues=('', '', '', ''))
    matrix(tmp_path, sequences={'a': 'ACGACGACG', 'b': 'ATGATGATG'}, parts=[('g1', 1, 3), ('g2', 4, 9)])
    matrix(tmp_path, 'codon/trimmed', sequences={'a': 'ACGACG', 'b': 'ATGATG'}, parts=[('g2', 1, 3), ('g3', 4, 6)])
    result = report.build_report(tmp_path)
    assert result['qc'][0]['passed'] == 4
    assert result['matrices'][0]['loci'] == 2
    comparison = result['trimming_comparisons'][0]
    assert comparison['common_loci'] == comparison['raw_only'] == comparison['trimmed_only'] == 1
    assert comparison['raw_common_columns'] == 6
    assert comparison['trimmed_common_columns'] == 3


def test_events_canonical_mapped_and_deduplicated(tmp_path):
    fields = ['input_file', 'taxa_id', 'mapped_taxa_id', 'event', 'action', 'source_nt_positions', 'source_aa_positions']
    event = dict(input_file='gene.fa', taxa_id='id1', mapped_taxa_id='Taxon1', event='internal_stop',
                 action='replace', source_nt_positions='4-6', source_aa_positions='NA')
    protein = dict(input_file='protein.fa', taxa_id='id1', mapped_taxa_id='Taxon1', event='protein_internal_stop',
                   action='replace', source_nt_positions='NA', source_aa_positions='4-4')
    tsv(tmp_path / 'reports/sequence_changes.tsv', fields, [event, event, protein])
    tsv(tmp_path / 'alignments/extra.sequence_changes.tsv', fields, [event] * 5)
    audit = report.build_report(tmp_path)['audit']
    assert audit['events'] == 2 and audit['duplicates_ignored'] == 1
    assert audit['taxa'][0] == dict(taxon='Taxon1', events=2, affected_sequences=2,
                                   nt_characters=3, aa_characters=1, events_without_coordinates=0)


def test_published_manifest_excludes_unmanaged_matrix(tmp_path):
    directory = matrix(tmp_path)
    (directory / 'fourfold.fasta').write_text('>a\nA\n')
    (directory / '.phyloprep-products.json').write_text(json.dumps([
        'supermatrix.fasta', 'supermatrix_part_iqtree.txt', 'supermatrix.missing.tsv', 'supermatrix.report.tsv']))
    hidden = tmp_path / 'matrices/.backup/codon/raw'
    hidden.mkdir(parents=True)
    (hidden / 'supermatrix.fasta').write_text('not a matrix')
    result = report.build_report(tmp_path)
    assert len(result['matrices']) == 1


def test_derived_loci_fourfold_not_assumed(tmp_path):
    directory = matrix(tmp_path)
    (directory / 'codon1st.supermatrix.fasta').write_text('>a\nA?\n>b\nA-\n')
    (directory / 'fourfold.fasta').write_text('>a\nC\n>b\nT\n')
    result = {m['kind']: m for m in report.build_report(tmp_path)['matrices']}
    assert result['codon1st']['loci'] == 2
    assert result['fourfold']['loci'] is None
    assert result['codon1st']['complete_codons'] is None


def test_portable_standalone_cli_and_safe_english_html(tmp_path):
    # Filesystem cannot contain slash in a name; the hostile taxon covers that case.
    root = tmp_path / 'run-<unsafe>'
    matrix(root, sequences={'</script><script>alert(1)</script>': 'ACG???', 'b': 'ATGNNN'})
    standalone = tmp_path / 'standalone'
    standalone.mkdir()
    shutil.copy(REPO / 'phyloprep_report.py', standalone)
    completed = subprocess.run([sys.executable, '-I', str(standalone / 'phyloprep_report.py'), '-i', str(root)],
                               text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    output = (root / 'phyloprep-report.html').read_text()
    assert '<html lang="en">' in output
    assert '<script>alert(1)</script>' not in output
    assert 'network requests' in output
    assert 'src="http' not in output and 'href="http' not in output
    assert 'Unique failed' in output and 'Save chart as SVG' in output
    payload = re.search(r'id="report-data">(.*?)</script>', output).group(1)
    assert json.loads(payload)['matrices'][0]['taxa'][0] == '</script><script>alert(1)</script>'
    assert not list(root.glob('.phyloprep-report-*'))


def test_long_unwrapped_matrix_reads_bounded(tmp_path, monkeypatch):
    source = tmp_path / 'long.fasta'
    sequence = 'ACGT' * 300_000
    source.write_text('>a\n' + sequence + '\n>b\n' + sequence + '\n')
    del sequence
    original = Path.open

    class Guard:
        def __init__(self, handle):
            self.handle = handle

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.handle.close()

        def readline(self, size=-1):
            assert 0 < size <= 65536
            return self.handle.readline(size)

    def guarded(path, *args, **kwargs):
        handle = original(path, *args, **kwargs)
        return Guard(handle) if path == source else handle

    monkeypatch.setattr(Path, 'open', guarded)
    tracemalloc.start()
    try:
        # Parsing stays small even when individual unwrapped lines exceed a MiB.
        spool = report.FastaSpool(source, tmp_path)
        assert spool.length == 1_200_000
        spool.close()
        assert tracemalloc.get_traced_memory()[1] < 1_000_000
    finally:
        tracemalloc.stop()


def test_actual_codon_pipeline_and_independent_pis(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCTGGT\n>b\nGCTGGT\n>c\nGCCGGA\n>d\nGCCGGA\n')
    root = tmp_path / 'results'
    completed = subprocess.run([sys.executable, str(REPO / 'phyloprep.py'), '-i', str(source),
                                '-o', str(root), '-t', '1', '--min-nucl-pis', '1'],
                               cwd=tmp_path, env=env, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    result = report.build_report(root)
    assert not result['warnings']
    assert result['overview']['input_loci'] == result['overview']['aligned_loci'] == 1
    assert len(result['matrices']) == 24
    assert all(m['loci'] == 1 for m in result['matrices'] if m['kind'] != 'fourfold')
    assert all(m['loci'] is None for m in result['matrices'] if m['kind'] == 'fourfold')
    assert all(q['passed'] == 1 for q in result['qc'])
    assert result['trimming_comparisons'][0]['common_loci'] == 1


def test_independent_threshold_failure_not_lost(tmp_path):
    directory = qc(tmp_path, branch='codon/raw', issues=('',))
    tsv(directory / 'alignment.details.tsv', ['File', 'Seqtype', 'Parsimony_informative_sites', 'Issues'],
        [dict(File='/old/0.mafft.codon.aln.fasta', Seqtype='codon', Parsimony_informative_sites=1, Issues='')])
    (directory / 'independent.retained.tsv').write_text('')
    tsv(tmp_path / 'qc/paired/raw/loci.details.tsv',
        ['Locus', 'Minimum_codon_PIS', 'Minimum_protein_PIS', 'Paired_status', 'Reasons'],
        [dict(Locus='0', Minimum_codon_PIS=2, Minimum_protein_PIS=0, Paired_status='fail', Reasons='pis_policy_failed')])
    result = report.build_report(tmp_path)
    independent = next(q for q in result['qc'] if q['branch'] == 'codon/raw')
    assert independent['failed'] == 1
    assert independent['reason_hits'] == {'below_independent_PIS_threshold': 1}


def test_invalid_folder_does_not_create_output(tmp_path):
    assert report.main(['-i', str(tmp_path)]) == 1
    assert not (tmp_path / 'phyloprep-report.html').exists()
