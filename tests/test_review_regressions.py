"""Regression cases from the repository-wide review."""
import csv
import importlib.util
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from alignment_qc import alignment_stats
from alignment_qc import alignment_qc
from filter_sequences import filter_fasta
from get_supergenes import get_supergenes
from rename_taxa import rename_fasta


def cli(script, tmp_path, *args, env=None):
    return subprocess.run([sys.executable, str(REPO / script), *map(str, args)],
                          cwd=tmp_path, env=env, capture_output=True, text=True, timeout=45)


@pytest.mark.parametrize('seqtype,states,expected', [
    ('prot', 'AAXX', 0), ('prot', 'AABB', 0), ('prot', 'AALL', 1),
    ('nucl', 'AARR', 0), ('nucl', 'AANN', 0), ('nucl', 'AACC', 1),
    ('nucl', 'TTUU', 0),
])
def test_pis_counts_resolved_states(tmp_path, seqtype, states, expected):
    source = tmp_path / 'input.fa'
    source.write_text(''.join(f'>t{i}\n{s}\n' for i, s in enumerate(states)))
    assert alignment_stats(source, seqtype)[0]['parsimony_informative_sites'] == expected


@pytest.mark.parametrize('seqtype', ['nucl', 'prot', 'codon'])
def test_qc_detail_columns_and_issue(tmp_path, seqtype):
    source = tmp_path / 'input.fa'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    prefix = tmp_path / 'qc'
    result = cli('alignment_qc.py', tmp_path, '-i', source, '-st', seqtype,
                 '-prefix', prefix, '--min-length', '100')
    with prefix.with_suffix('.details.tsv').open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert len(rows) == 2, result.stderr
    assert len(rows[0]) == len(rows[1])
    assert rows[0][-1] == 'Issues'
    assert 'short_alignment' in rows[1][-1].split(';')


def test_codon_min_length_counts_complete_triplets(tmp_path):
    source = tmp_path / 'input.fa'
    source.write_text('>a\nGCTGGT\n>b\nGCCGGA\n')
    assert 'short_sequence' not in alignment_qc(source, 'codon', min_length=2)['issues']
    assert 'short_sequence' in alignment_qc(source, 'codon', min_length=3)['issues']


def test_codon_n_ratio_counts_ambiguous_triplets(tmp_path):
    source = tmp_path / 'input.fa'
    source.write_text('>a\nNAAGCT\n>b\nAAAGCC\n')
    # One of four codons is ambiguous. Counting bases would incorrectly
    # report 1/12 rather than the protein-comparable 1/4.
    assert alignment_stats(source, 'codon')[0]['N_ratio'] == 0.25


@pytest.mark.parametrize('alias', ['same', 'symlink', 'hardlink'])
@pytest.mark.parametrize('utility', ['supergenes', 'rename', '4dtv'])
def test_output_collision_preserves_inputs(tmp_path, alias, utility):
    source = tmp_path / 'input.fasta'
    original = '>a\nGCT\n>b\nGCC\n'
    source.write_text(original)
    destination = source if alias == 'same' else tmp_path / 'alias.fasta'
    if alias == 'symlink': destination.symlink_to(source)
    if alias == 'hardlink': os.link(source, destination)
    with pytest.raises(ValueError, match='different files'):
        if utility == 'supergenes':
            prefix = str(destination)[:-len('.fasta')]
            get_supergenes([str(source)], prefix)
        elif utility == 'rename':
            rename_fasta(source, tmp_path / 'renamed.fa', map_output=destination)
        else:
            spec = importlib.util.spec_from_file_location('fourdtv', REPO / 'calulate_4dtv_and_correction.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.main([str(source)], str(destination))
    assert source.read_text() == original
    assert not (tmp_path / 'renamed.fa').exists()


def test_rename_checks_all_destinations_before_writing(tmp_path):
    source, mapping, output = (tmp_path / name for name in ('in.fa', 'map.tsv', 'out.fa'))
    source.write_text('>a\nACG\n')
    mapping.write_text('a\tb\n')
    output.write_text('previous result')
    for destination in (mapping, output):
        with pytest.raises(ValueError, match='different files'):
            rename_fasta(source, output, mapping, map_output=destination)
        assert output.read_text() == 'previous result'
        assert mapping.read_text() == 'a\tb\n'


def test_empty_keep_list_preserves_output(tmp_path):
    source, output = tmp_path / 'in.fa', tmp_path / 'out.fa'
    source.write_text('>a\nACG\n>b\nACC\n')
    output.write_text('previous result')
    with pytest.raises(ValueError, match='list is empty'):
        filter_fasta(source, output, 'nucl', keep_ids=set())
    assert output.read_text() == 'previous result'


@pytest.mark.parametrize('seqtype,sequence,flag', [
    ('prot', 'MAGAGG', '--amino'), ('codon', 'ATGGCTGGT', '--amino'), ('nucl', 'ACGTAC', '--nuc')])
def test_mafft_alphabet_is_explicit(tmp_path, fake_aligner, seqtype, sequence, flag):
    executable, env = fake_aligner
    executable.write_text(executable.read_text().replace(
        'import os, sys, time, signal',
        f'import os, sys, time, signal\nassert {flag!r} in sys.argv'))
    source = tmp_path / 'input.fa'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    result = cli('MSAP.py', tmp_path, '-i', source, '-st', seqtype, '--notrim', '-t', 1, env=env)
    assert result.returncode == 0, result.stderr


def test_failed_rerun_excludes_stale_matrix(tmp_path, fake_aligner):
    _, env = fake_aligner
    source, output = tmp_path / 'gene.fa', tmp_path / 'results'
    source.write_text('>a\nACGT\n>b\nACGA\n')
    args = ['-i', source, '-o', output, '-st', 'nucl', '--notrim', '-t', 1]
    first = cli('MSAP_batch.py', tmp_path, *args, env=env)
    assert first.returncode == 0, first.stderr
    paths = output / 'path-lists/all.nucl.aln.pathlist'
    assert paths.read_text().strip()
    source.write_text('>a\nTTTT\n>b\nTTTA\n')
    result = cli('MSAP_batch.py', tmp_path, *args, env={**env, 'MSAP_TEST_FAIL': '1'})
    assert result.returncode != 0
    assert paths.read_text() == ''
    # Failed replacement keeps the old matrix for recovery, but never lists it.
    assert 'ACGT' in (output / 'gene.mafft.nucl.aln').read_text()


@pytest.mark.parametrize('interrupt', [signal.SIGINT, signal.SIGTERM])
def test_pipeline_interrupt_reaps_aligner(tmp_path, fake_aligner, interrupt):
    _, env = fake_aligner
    env['MSAP_TEST_SLEEP'] = '1'
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nACGT\n>b\nACGA\n')
    child_pid = None
    process = subprocess.Popen([sys.executable, str(REPO / 'phyloprep.py'), '-i', str(source),
                                '-o', str(tmp_path / 'out'), '-st', 'nucl', '--notrim', '-t', '1'],
                               cwd=tmp_path, env=env, start_new_session=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        deadline = time.monotonic() + 20
        pidfile = Path(env['MSAP_TEST_PID'])
        while not pidfile.exists() and time.monotonic() < deadline:
            assert process.poll() is None
            time.sleep(0.05)
        child_pid = int(pidfile.read_text())
        process.send_signal(interrupt)
        stdout, stderr = process.communicate(timeout=20)
        assert process.returncode == 130, stdout + stderr
        # An orphan zombie may await init reaping; it must not still be running.
        stat = Path(f'/proc/{child_pid}/stat')
        assert not stat.exists() or stat.read_text().split()[2] == 'Z'
    finally:
        if child_pid:
            try: os.killpg(os.getpgid(child_pid), signal.SIGKILL)
            except ProcessLookupError: pass
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
        process.communicate(timeout=10)
