import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import MSAP_batch as batch
from alignment_qc import alignment_qc
from sequence_audit import preprocess_codon
from msap_io import is_stop_codon


def cli(script, *args, cwd, env=None):
    return subprocess.run([sys.executable, str(ROOT/script), *map(str, args)],
                          cwd=cwd, env=env, capture_output=True, text=True)


@pytest.mark.parametrize('code,codon', [(27, 'TGA'), (28, 'TAA'), (28, 'TGA'), (31, 'TAA')])
def test_dual_coding_retained_internally_and_not_rejected_after_cleaning(tmp_path, code, codon):
    source, output, report = [tmp_path/n for n in ('input.fa', 'clean.fa', 'changes.tsv')]
    source.write_text(f'>a\nATG{codon}{codon}\n>b\nATG{codon}{codon}\n')
    preprocess_codon(source, output, code, report)
    assert {str(r.seq) for r in SeqIO.parse(output, 'fasta')} == {'ATG'+codon}
    with report.open() as handle:
        changes = list(csv.DictReader(handle, delimiter='\t'))
    assert len(changes) == 2 and all(r['event'] == 'terminal_stop' for r in changes)
    assert not is_stop_codon(codon, code)
    assert is_stop_codon(codon, code, terminal=True)
    assert 'stop_codon' not in alignment_qc(output, 'codon', code)['issues']


def test_code27_real_pipeline_preserves_internal_tryptophan(tmp_path):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nATGTGAGCTTGA\n>b\nATGTGAGCCTGA\n')
    out = tmp_path/'out'
    result = cli('phyloprep.py', '-i', source, '-o', out, '-g', 27, '-t', 1, '--notrim', cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    matrix = AlignIO.read(out/'matrices/prot/raw/supermatrix.fasta', 'fasta')
    assert {str(r.seq) for r in matrix} == {'MWA'}


def test_real_protein_pipeline_masks_and_audits_special_residues(tmp_path):
    source = tmp_path/'protein.fa'
    source.write_text('>a\nMAjouG\n>b\nMAAALG\n')
    out = tmp_path/'out'
    result = cli('phyloprep.py', '-i', source, '-o', out, '-st', 'prot', '-t', 1, '--notrim', cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    matrices = []
    for extension, fmt in [('fasta', 'fasta'), ('phy', 'phylip-relaxed'), ('nex', 'nexus')]:
        matrices.append({r.id: str(r.seq) for r in AlignIO.read(out/f'matrices/prot/raw/supermatrix.{extension}', fmt)})
    assert matrices[0] == matrices[1] == matrices[2]
    assert matrices[0]['a'].replace('-', '') == 'MAXXXG'
    with (out/'reports/sequence_changes.tsv').open() as handle:
        changes = list(csv.DictReader(handle, delimiter='\t'))
    assert [(r['original'], r['source_aa_positions']) for r in changes] == [('j', '3-3'), ('o', '4-4'), ('u', '5-5')]
    assert all(r['event'] == 'protein_unsupported_residue' and r['replacement'] == 'X' for r in changes)


@pytest.mark.parametrize('symbol', ['J', 'O', 'U'])
def test_nexus_masks_unhandled_protein_symbols_with_audit(tmp_path, symbol):
    source, output = tmp_path/'protein.fa', tmp_path/'result.nex'
    source.write_text(f'>a\nMA{symbol}G\n>b\nMAAG\n')
    result = cli('fasta2nex.py', '-i', source, '-o', output, '-st', 'protein', cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert {r.id: str(r.seq) for r in AlignIO.read(output, 'nexus')}['a'] == 'MAXG'
    report = output.with_name(output.name + '.sequence_changes.tsv')
    rows = list(csv.DictReader(report.open(), delimiter='\t'))
    assert [(row['original'], row['replacement'], row['source_aa_positions']) for row in rows] == [(symbol, 'X', '3-3')]


def test_rename_whitespace_rejected_without_overwriting(tmp_path):
    source, mapping, output = [tmp_path/n for n in ('input.fa', 'map.tsv', 'out.fa')]
    source.write_text('>a\nACG\n>b\nACG\n')
    mapping.write_text('a\tSpecies one\nb\tSpecies two\n')
    output.write_text('previous result')
    result = cli('rename_taxa.py', '-i', source, '-m', mapping, '-o', output, '--keep-special', cwd=tmp_path)
    assert result.returncode != 0 and 'whitespace' in result.stderr
    assert output.read_text() == 'previous result'


@pytest.mark.parametrize('script', ['get_supergenes.py', 'calulate_4dtv_and_correction.py'])
@pytest.mark.parametrize('link', ['same', 'symlink', 'hardlink'])
def test_path_lists_are_protected(tmp_path, script, link):
    source = tmp_path/'gene.fa'
    source.write_text('>a\nGCT\n>b\nGCC\n')
    output = tmp_path/('out.fasta' if script == 'get_supergenes.py' else 'out.tsv')
    manifest = output if link == 'same' else tmp_path/'inputs.list'
    original = str(source)+'\n'
    manifest.write_text(original)
    if link == 'symlink':
        output.symlink_to(manifest)
    elif link == 'hardlink':
        os.link(manifest, output)
    options = ['-p', tmp_path/'out'] if script == 'get_supergenes.py' else ['-o', output]
    result = cli(script, '-i', manifest, *options, cwd=tmp_path)
    assert result.returncode != 0 and 'different files' in result.stderr
    assert manifest.read_text() == output.read_text() == original
    assert not (tmp_path/'out.report.tsv').exists()


@pytest.mark.parametrize('committed', [False, True])
def test_recovery_after_partial_backup_cleanup(tmp_path, monkeypatch, committed):
    key = 'gene_test'
    stage = tmp_path/'.msap-batch-staging'/key
    stage.mkdir(parents=True)
    name = 'gene.mafft.prot.aln'
    (tmp_path/name).write_text('OLD')
    (stage/name).write_text('NEW')
    state = tmp_path/f'.msap-batch-state-{key}.json'
    running = {'status': 'running', 'outputs': {name: {}}}
    batch._atomic_json(state, running)
    backup = tmp_path/'.msap-batch-backups'/key
    real_replace, real_rmtree = batch.os.replace, batch.shutil.rmtree
    failed = False

    def replace(source, dest):
        nonlocal failed
        if not committed and Path(dest) == state and not failed:
            failed = True
            raise OSError('checkpoint failure')
        return real_replace(source, dest)

    def cleanup(path, *args, **kwargs):
        if Path(path) == backup:
            real_rmtree(backup/'files')
            raise OSError('cleanup interrupted')
        return real_rmtree(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(batch.os, 'replace', replace)
        patch.setattr(batch.shutil, 'rmtree', cleanup)
        with pytest.raises(OSError, match='cleanup'):
            batch._publish_task(stage, tmp_path, state, running, 'gene.sequence_changes.tsv')
    expected = 'NEW' if committed else 'OLD'
    assert (tmp_path/name).read_text() == expected
    batch._recover_publication(backup, tmp_path)
    assert (tmp_path/name).read_text() == expected
    assert not backup.exists()


@pytest.mark.parametrize('legacy', [False, True])
def test_missing_backup_is_never_interpreted_as_permission_to_delete(tmp_path, legacy):
    backup = tmp_path/'.msap-batch-backups/gene_test'
    backup.mkdir(parents=True)
    output = tmp_path/'old.aln'
    output.write_text('KEEP')
    journal = {'names': ['old.aln'], 'outputs': {}, 'running_state': {}}
    if not legacy:
        journal['original_files'] = ['old.aln']
    (backup/'journal.json').write_text(json.dumps(journal))
    with pytest.raises(ValueError):
        batch._recover_publication(backup, tmp_path)
    assert output.read_text() == 'KEEP'
