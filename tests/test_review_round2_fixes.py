"""Stop interpretation, output preservation, and matrix publication boundaries."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from alignment_qc import alignment_qc
from msap_io import is_stop_codon
from sequence_audit import preprocess_codon
import phyloprep


def cli(tmp_path, script, *args, data=None, env=None):
    return subprocess.run([sys.executable, str(REPO / script), *map(str, args)],
                          cwd=tmp_path, input=data, capture_output=True, text=True,
                          env=env, timeout=30)


@pytest.mark.parametrize('codon,code,expected', [
    ('TAR', 1, True), ('TRA', 1, True), ('UGA', 1, True), ('uar', 1, True),
    ('TAN', 1, False), ('TRR', 1, False), ('NNN', 1, False), ('---', 1, False),
    ('TGA', 2, False), ('AGR', 2, True), ('AGR', 1, False), ('TAR', 6, False)])
def test_definite_stops_only(codon, code, expected):
    assert is_stop_codon(codon, code) is expected


def test_ambiguous_stops_preprocessing_and_report(tmp_path):
    source, output, report = (tmp_path / n for n in ('in.fa', 'out.fa', 'changes.tsv'))
    source.write_text('>a\nATGTARGCTTAR\n>b\nATGTANGCTGCT\n')
    preprocess_codon(source, output, 1, report)
    assert output.read_text() == '>a\nATGNNNGCT\n>b\nATGTANGCTGCT\n'
    with report.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [(r['event'], r['source_nt_positions'], r['original'], r['replacement']) for r in rows] == [
        ('internal_stop', '4-6', 'TAR', 'NNN'), ('terminal_stop', '10-12', 'TAR', '-')]


@pytest.mark.parametrize('sequence,code,stops', [
    ('ATGTARGCT', 1, '2-TAR'), ('AUGUGAGCU', 1, '2-UGA'),
    ('ATGTANGCT', 1, 'NA'), ('ATGTGAGCT', 2, 'NA')])
def test_qc_uses_same_stop_rules(tmp_path, sequence, code, stops):
    source = tmp_path / 'in.fa'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    result = alignment_qc(source, 'codon', code)
    assert result['sequences'][0]['internal_stop_count'] == stops
    assert ('stop_codon' in result['issues']) == (stops != 'NA')


def test_msap_internal_ambiguous_stop_backtranslates(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nATGTARGCT\n>b\nATGGCTGCT\n')
    result = cli(tmp_path, 'MSAP.py', '-i', source, '-st', 'codon', '--notrim', '-t', 1, env=env)
    assert result.returncode == 0, result.stderr
    assert 'ATGNNNGCT' in (tmp_path / 'gene.mafft.codon.aln').read_text()


@pytest.mark.parametrize('suffix', ['pass', 'fail', 'details'])
@pytest.mark.parametrize('alias', ['same', 'symlink', 'hardlink'])
def test_qc_protects_inputs_before_any_report_write(tmp_path, suffix, alias):
    source = tmp_path / ('qc.' + suffix + '.tsv' if alias == 'same' else 'in.fa')
    content = '>a\nGCT\n>b\nGCC\n'
    source.write_text(content)
    for name in ('pass', 'fail', 'details'):
        path = tmp_path / f'qc.{name}.tsv'
        if name == suffix:
            if alias == 'symlink': path.symlink_to(source)
            elif alias == 'hardlink': os.link(source, path)
        else:
            path.write_text('old ' + name)
    result = cli(tmp_path, 'alignment_qc.py', '-i', source, '-st', 'codon', '-p', tmp_path / 'qc')
    assert result.returncode != 0
    assert 'different files' in result.stderr
    assert source.read_text() == content
    for name in ('pass', 'fail', 'details'):
        if name != suffix:
            assert (tmp_path / f'qc.{name}.tsv').read_text() == 'old ' + name


def test_qc_protects_input_pathlist(tmp_path):
    source, listing = tmp_path / 'in.fa', tmp_path / 'qc.pass.tsv'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    listing.write_text(str(source) + '\n')
    result = cli(tmp_path, 'alignment_qc.py', '-i', listing, '-st', 'codon', '-p', tmp_path / 'qc')
    assert result.returncode != 0
    assert listing.read_text() == str(source) + '\n'
    assert not (tmp_path / 'qc.details.tsv').exists()


@pytest.mark.parametrize('stdin', [False, True])
@pytest.mark.parametrize('explicit_stdout', [False, True])
def test_trim_standard_streams(tmp_path, stdin, explicit_stdout):
    content = '>a\nGCT\n>b\nGCC\n'
    source = tmp_path / 'in.fa'; source.write_text(content)
    args = ['-i', '-' if stdin else source]
    if explicit_stdout: args += ['-o', '-']
    result = cli(tmp_path, 'trimAlnSeq.py', *args, data=content if stdin else None)
    assert result.returncode == 0, result.stderr
    assert result.stdout == content
    assert not (tmp_path / '-').exists()


@pytest.mark.parametrize('script', ['filter_sequences.py', 'trimAlnSeq.py'])
@pytest.mark.parametrize('seqtype,sequence', [('nucl', 'NNN'), ('prot', 'XXX'), ('codon', 'NNN')])
def test_zero_columns_preserve_existing_output(tmp_path, script, seqtype, sequence):
    source, output = tmp_path / 'in.fa', tmp_path / 'out.fa'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    output.write_text('previous matrix')
    result = cli(tmp_path, script, '-i', source, '-o', output, '-st', seqtype)
    assert result.returncode != 0
    assert 'No alignment sites remain' in result.stderr
    assert output.read_text() == 'previous matrix'


@pytest.mark.parametrize('sequence,valid', [('AT-GCT', False), ('A--GCT', False),
    ('AT-G-TACT', False), ('---GCT', True), ('ATGGCT', True), ('NNNGCT', True)])
def test_codon_qc_triplet_gaps(tmp_path, sequence, valid):
    source = tmp_path / 'in.fa'
    source.write_text(f'>a\n{sequence}\n>b\n' + 'GCT' * (len(sequence) // 3) + '\n')
    result = alignment_qc(source, 'codon')
    assert result['sequences'][0]['frame_ok'] is valid
    assert ('frame_error' in result['issues']) is not valid


def trees(tmp_path):
    root, stage = tmp_path / 'out', tmp_path / 'stage'
    root.mkdir(); stage.mkdir()
    for kind in ('codon', 'prot'):
        old = root / 'matrices' / kind / 'raw'; old.mkdir(parents=True)
        new = stage / kind / 'raw'; new.mkdir(parents=True)
        (old / 'supermatrix.fasta').write_text('old ' + kind)
        (old / 'obsolete.phy').write_text('obsolete')
        (old / 'user.txt').write_text('user data')
        (old / '.phyloprep-products.json').write_text(json.dumps(['supermatrix.fasta', 'obsolete.phy']))
        (new / 'supermatrix.fasta').write_text('new ' + kind)
    return root, stage


def snapshot(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob('*') if p.is_file()}


def test_invalid_later_manifest_preserves_every_group(tmp_path):
    root, stage = trees(tmp_path)
    (root / 'matrices/prot/raw/.phyloprep-products.json').write_text('{')
    before = snapshot(root / 'matrices')
    with pytest.raises(ValueError):
        phyloprep.publish_matrices(stage, root, ('raw',))
    assert snapshot(root / 'matrices') == before


@pytest.mark.parametrize('failure', [OSError, KeyboardInterrupt])
def test_failed_directory_switch_rolls_back(tmp_path, monkeypatch, failure):
    root, stage = trees(tmp_path)
    before = snapshot(root / 'matrices')
    replace = os.replace
    def fail_install(source, target):
        if Path(source).name == 'matrices' and Path(source).parent.name.startswith('.publish-'):
            raise failure('injected install failure')
        return replace(source, target)
    monkeypatch.setattr(phyloprep.os, 'replace', fail_install)
    with pytest.raises(failure):
        phyloprep.publish_matrices(stage, root, ('raw',))
    assert snapshot(root / 'matrices') == before
    assert not (root / '.phyloprep-matrices-backup').exists()


def test_success_preserves_unmanaged_files_and_other_groups(tmp_path):
    root, stage = trees(tmp_path)
    other = root / 'matrices/nucl/raw'; other.mkdir(parents=True)
    (other / 'supermatrix.fasta').write_text('other run')
    phyloprep.publish_matrices(stage, root, ('raw',))
    for kind in ('codon', 'prot'):
        group = root / 'matrices' / kind / 'raw'
        assert (group / 'supermatrix.fasta').read_text() == 'new ' + kind
        assert (group / 'user.txt').read_text() == 'user data'
        assert not (group / 'obsolete.phy').exists()
    assert (other / 'supermatrix.fasta').read_text() == 'other run'
    assert not (root / '.phyloprep-matrices-backup').exists()


def test_recovers_interrupted_switch(tmp_path):
    root, _ = trees(tmp_path)
    before = snapshot(root / 'matrices')
    os.replace(root / 'matrices', root / '.phyloprep-matrices-backup')
    phyloprep.recover_matrix_publication(root)
    assert snapshot(root / 'matrices') == before


def test_failed_rollback_retains_recoverable_backup(tmp_path, monkeypatch):
    root, stage = trees(tmp_path)
    before = snapshot(root / 'matrices')
    replace = os.replace
    def fail_install_and_restore(source, target):
        if Path(target) == root / 'matrices':
            raise OSError('injected filesystem failure')
        return replace(source, target)
    with monkeypatch.context() as patch:
        patch.setattr(phyloprep.os, 'replace', fail_install_and_restore)
        with pytest.raises(OSError):
            phyloprep.publish_matrices(stage, root, ('raw',))
    assert snapshot(root / '.phyloprep-matrices-backup') == before
    phyloprep.recover_matrix_publication(root)
    assert snapshot(root / 'matrices') == before


def test_first_publication_failure_leaves_no_partial_tree(tmp_path, monkeypatch):
    root, stage = tmp_path / 'out', tmp_path / 'stage'
    root.mkdir()
    group = stage / 'nucl/raw'; group.mkdir(parents=True)
    (group / 'supermatrix.fasta').write_text('new matrix')
    replace = os.replace
    def fail_install(source, target):
        if Path(target) == root / 'matrices':
            raise OSError('injected install failure')
        return replace(source, target)
    monkeypatch.setattr(phyloprep.os, 'replace', fail_install)
    with pytest.raises(OSError):
        phyloprep.publish_matrices(stage, root, ('raw',))
    assert not (root / 'matrices').exists()


def test_trim_zero_columns_writes_nothing_to_stdout(tmp_path):
    result = cli(tmp_path, 'trimAlnSeq.py', '-i', '-', data='>a\nNNN\n>b\nNNN\n')
    assert result.returncode != 0
    assert result.stdout == ''
    assert not (tmp_path / '-').exists()
