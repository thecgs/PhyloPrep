"""Release regressions: RNA, publication recovery, and lossless identifiers."""
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import AlignIO, SeqIO
from Bio.Nexus import Nexus

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import MSAP_batch as batch
from alignment_qc import alignment_stats
from get_supergenes import tidy_name
from sequence_audit import preprocess_codon, write_report, event


def cli(script, *args, cwd, env=None):
    return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)],
                          cwd=cwd, env=env, capture_output=True, text=True)


def records(path, fmt='fasta'):
    return {r.id: str(r.seq) for r in AlignIO.read(path, fmt)}


def check_report(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        source = {r.id: str(r.seq) for r in SeqIO.parse(row['input_file'], 'fasta')}
        start, end = map(int, row['source_nt_positions'].split('-'))
        assert source[row['taxa_id']][start-1:end] == row['original']
    return rows


@pytest.mark.parametrize('mode', ['nucl', 'codon'])
def test_rna_pipeline_raw_trimmed_and_conversions(tmp_path, fake_aligner, mode):
    _, env = fake_aligner
    source = tmp_path / 'rna.fa'
    text = '>a\nguggcuugagcuuaa\n>b\nGUGGCCUGAGCCUAA\n'
    source.write_text(text)
    out = tmp_path / 'out'
    result = cli('phyloprep.py', '-i', source, '-st', mode, '-t', 1, '-o', out, cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    products = list((out / 'matrices').rglob('*.fasta'))
    assert products
    for path in products:
        matrix = records(path)
        assert matrix == records(path.with_suffix('.phy'), 'phylip-relaxed')
        assert matrix == records(path.with_suffix('.nex'), 'nexus')
        if '/prot/' not in str(path):
            assert all('U' not in seq for seq in matrix.values())
    assert any('/raw/' in str(p) for p in products)
    assert any('/trimmed/' in str(p) for p in products)
    if mode == 'codon':
        assert any(p.name == 'fourfold.fasta' for p in products)
        assert any(p.name == 'codon3rd.supermatrix.fasta' for p in products)
    rows = check_report(out / 'reports/sequence_changes.tsv')
    assert any(r['event'] == 'uracil_to_thymine' for r in rows)
    assert source.read_text() == text


def test_rna_stop_events_do_not_overlap_normalization(tmp_path):
    source, output, report = [tmp_path / n for n in ('rna.fa', 'dna.fa', 'changes.tsv')]
    source.write_text('>a\ngugugagcuuaa\n>b\nGUGUGAGCUUAA\n')
    preprocess_codon(source, output, 1, report)
    assert records(output) == {'a': 'GTGNNNGCT', 'b': 'GTGNNNGCT'}
    rows = check_report(report)
    for taxon in ('a', 'b'):
        positions = []
        for row in rows:
            if row['taxa_id'] == taxon:
                start, end = map(int, row['source_nt_positions'].split('-'))
                positions.extend(range(start, end+1))
        assert len(positions) == len(set(positions))


def test_rna_backtranslation_initiator_and_fourfold(tmp_path):
    cds, protein, output = [tmp_path / n for n in ('rna.fa', 'aa.fa', 'codon.fa')]
    cds.write_text('>a\nguggcu\n>b\nGUGGCC\n')
    protein.write_text('>a\nVA\n>b\nVA\n')
    result = cli('AA2Codon.py', '-c', cds, '-p', protein, '-g', 1, '-o', output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert records(output) == {'a': 'GTGGCT', 'b': 'GTGGCC'}
    fourfold = tmp_path / 'fourfold.fa'
    result = cli('extract_4-fold_degenerated_sites.py', '-i', cds, '-o', fourfold, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert records(fourfold) == {'a': 'GT', 'b': 'GC'}
    statistics = tmp_path/'4dtv.tsv'
    result = cli('calulate_4dtv_and_correction.py', '-i', cds, '-o', statistics, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert statistics.read_text().splitlines()[-1].split('\t')[5:7] == ['2', '0']


def test_qc_t_and_u_have_identical_patterns(tmp_path):
    source = tmp_path / 'rna.fa'
    source.write_text('>a\nUT\n>b\nTU\n')
    stats, aligned = alignment_stats(source, 'nucl')
    assert aligned and stats['distinct_patterns'] == 1 and stats['variable_sites'] == 0
    stats, _ = alignment_stats(source, 'prot')
    assert stats['distinct_patterns'] == 2


@pytest.mark.parametrize('mode,flag,sequence', [('prot', '-protein', 'AGAG'), ('codon', '-protein', 'GCTGGT'), ('nucl', '-DNA', 'ACGU')])
def test_prank_explicit_alphabet(tmp_path, fake_aligner, mode, flag, sequence):
    executable, env = fake_aligner
    prank = executable.with_name('prank')
    prank.write_text(f'''#!{sys.executable}
import sys
from pathlib import Path
assert {flag!r} in sys.argv
source = next(a[3:] for a in sys.argv if a.startswith('-d='))
output = next(a[3:] for a in sys.argv if a.startswith('-o='))
Path(output + '.best.fas').write_text(Path(source).read_text())
''')
    prank.chmod(0o755)
    source = tmp_path / 'gene.fa'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    result = cli('MSAP.py', '-i', source, '-st', mode, '-s', 'prank', '--notrim', cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr


def test_protein_stops_and_u_preprocessed(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'protein.fa'
    source.write_text('>a\nMA*G\n>b\nMAAG\n')
    result = cli('MSAP.py', '-i', source, '-st', 'prot', '--notrim', cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr
    assert records(tmp_path / 'protein.mafft.prot.aln')['a'] == 'MAXG'
    assert Path(env['MSAP_TEST_LOG']).exists()
    source.write_text('>a\nMAUG\n>b\nMAAG\n')
    result = cli('MSAP.py', '-i', source, '-st', 'prot', '--notrim', cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr
    assert records(tmp_path / 'protein.mafft.prot.aln')['a'] == 'MAXG'


def test_aligner_cannot_silently_change_residues(tmp_path, fake_aligner):
    executable, env = fake_aligner
    executable.write_text(f'#!{sys.executable}\nprint(">a\\nMA-G\\n>b\\nMAAG")\n')
    source = tmp_path / 'protein.fa'
    source.write_text('>a\nMAUG\n>b\nMAAG\n')
    result = cli('MSAP_batch.py', '-i', source, '-st', 'prot', '--notrim', '-o', tmp_path/'out', cwd=tmp_path, env=env)
    assert result.returncode != 0 and 'non-gap sequence content' in result.stderr
    assert not (tmp_path/'out/protein.mafft.prot.aln').exists()


@pytest.mark.parametrize('text', ['', '# empty list\n; comment\n'])
def test_empty_input_cannot_create_batch_state(tmp_path, text):
    source = tmp_path / 'empty.list'
    source.write_text(text)
    out = tmp_path/'out'
    result = cli('MSAP_batch.py', '-i', source, '-o', out, cwd=tmp_path)
    assert result.returncode != 0 and 'Empty input' in result.stderr
    assert not out.exists()


def test_valid_empty_ownership_can_be_reused(tmp_path):
    (tmp_path / '.msap-batch-owners.json').write_text('{}')
    source = tmp_path/'gene.fa'
    batch._claim_prefixes([source], tmp_path)
    assert json.loads((tmp_path / '.msap-batch-owners.json').read_text()) == {'gene': str(source)}


def test_mode_switch_retires_managed_pathlists(tmp_path, fake_aligner):
    _, env = fake_aligner
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nGCTGGT\n>b\nGCCGGA\n')
    out = tmp_path/'out'
    result = cli('MSAP_batch.py', '-i', source, '-o', out, '-t', 1, cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr
    custom = out/'path-lists/custom.pathlist'
    custom.write_text('user data')
    result = cli('MSAP_batch.py', '-i', source, '-o', out, '-st', 'nucl', '--notrim', '-t', 1, cwd=tmp_path, env=env)
    assert result.returncode == 0, result.stderr
    assert {p.name for p in custom.parent.iterdir()} == {'custom.pathlist', 'all.nucl.aln.pathlist'}
    assert custom.read_text() == 'user data'


def test_identifiers_preserved_across_formats(tmp_path):
    ids = ['a[1]', 'a1', 'a:b;c,(d)', 'a|b']
    source = tmp_path/'gene.fa'
    source.write_text(''.join(f'>{name}\nACGT\n' for name in ids))
    for script, extension, fmt in [('fasta2phy.py', 'phy', 'phylip-relaxed'), ('fasta2nex.py', 'nex', 'nexus')]:
        output = source.with_suffix('.'+extension)
        result = cli(script, '-i', source, '-o', output, cwd=tmp_path)
        assert result.returncode == 0, result.stderr
        assert records(output, fmt) == records(source)


def test_partition_identifiers_are_safe_and_parseable(tmp_path):
    names = ['gene 1.fa', '9;[gene].fa']
    for name in names:
        assert re.fullmatch('[A-Za-z][A-Za-z0-9_]*', tidy_name(name))
    source = tmp_path/names[0]
    source.write_text('>a\nACG\n>b\nACG\n')
    result = cli('get_supergenes.py', '-i', source, '-p', tmp_path/'out', cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    partition = next(tmp_path.glob('out*iqtree*'))
    # IQ-TREE's empty model selector is an extension to the NEXUS grammar.
    text = re.sub(r'charpartition[^;]+;', '', partition.read_text())
    text = text.replace('#nexus', '#nexus\nbegin data; dimensions ntax=2 nchar=3; '
                        'format datatype=dna; matrix\na ACG\nb ACG\n; end;')
    assert Nexus.Nexus(text).charsets == {'gene_1_fa': [0, 1, 2]}


@pytest.mark.parametrize('seqtype,source_seq,expected', [('DNA', 'aUg', 'ATG'), ('RNA', 'aTg', 'AUG'), ('protein', 'MAXG', 'MAXG')])
def test_nexus_uses_selected_alphabet(tmp_path, seqtype, source_seq, expected):
    source, output = tmp_path/'in.fa', tmp_path/'out.nex'
    source.write_text(f'>a\n{source_seq}\n>b\n{source_seq}\n')
    result = cli('fasta2nex.py', '-i', source, '-o', output, '-st', seqtype, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert records(output, 'nexus') == {'a': expected, 'b': expected}


@pytest.mark.parametrize('failure', ['native', 'report', 'checkpoint', 'crash'])
def test_publication_failure_preserves_reports_and_recovers(tmp_path, monkeypatch, failure):
    root = tmp_path/'out'
    root.mkdir()
    key = 'gene_123'
    stage = root/'.msap-batch-staging'/key
    stage.mkdir(parents=True)
    native = 'gene.macse_original.NT.fasta'
    report = 'gene.sequence_changes.tsv'
    (root/native).write_text('>a\nGCT\n')
    (stage/native).write_text('>a\nTGA\n')
    write_report(root/report, [event(root/native, 'a', 'export_change', 1,
                 source_nt_positions='1-3', original='GCT', replacement='NNN', action='replace')])
    write_report(stage/report, [event(stage/native, 'a', 'internal_stop', 1,
                 source_nt_positions='1-3', original='TGA', replacement='NNN', action='replace')])
    obsolete = root/'old.aln'
    obsolete.write_text('old alignment')
    state = root/f'.msap-batch-state-{key}.json'
    running = {'status': 'running', 'outputs': {native: {}, report: {}, 'old.aln': {}}}
    batch._atomic_json(state, running)
    real_replace = batch.os.replace
    real_recover = batch._recover_publication
    target = state if failure in {'checkpoint', 'crash'} else root/(report if failure == 'report' else native)
    failed = False

    def fail_once(source, destination):
        nonlocal failed
        if Path(destination) == target and not failed:
            failed = True
            raise OSError('simulated publication failure')
        return real_replace(source, destination)

    monkeypatch.setattr(batch.os, 'replace', fail_once)
    if failure == 'crash':
        monkeypatch.setattr(batch, '_recover_publication', lambda *args: None)
    with pytest.raises(OSError, match='simulated'):
        batch._publish_task(stage, root, state, running, report)
    if failure == 'crash':
        assert (root/native).read_text() == '>a\nTGA\n'
        assert not obsolete.exists()
        real_recover(root/'.msap-batch-backups'/key, root)
    assert check_report(stage/report)[0]['original'] == 'TGA'
    assert check_report(root/report)[0]['original'] == 'GCT'
    assert obsolete.read_text() == 'old alignment'
    assert json.loads(state.read_text())['status'] == 'running'
    monkeypatch.setattr(batch, '_recover_publication', real_recover)
    batch._publish_task(stage, root, state, running, report)
    assert check_report(root/report)[0]['original'] == 'TGA'
    assert not stage.exists() and not obsolete.exists()
    assert json.loads(state.read_text())['status'] == 'complete'
