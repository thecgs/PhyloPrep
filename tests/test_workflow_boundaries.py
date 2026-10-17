"""Regression coverage for mapping/provenance and matrix publication boundaries."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import alignment_qc
import phyloprep
from rename_taxa import read_mapping_snapshot
from sequence_audit import event, write_report

ROOT = Path(__file__).resolve().parents[1]


def test_mapping_changed_after_read_cannot_be_checkpointed(tmp_path):
    mapping = tmp_path/'map.tsv'
    mapping.write_text('a\tOLD_A\nb\tOLD_B\n')
    parsed, signature = read_mapping_snapshot(mapping)
    args = SimpleNamespace(mapping=mapping, taxon_mapping=parsed,
                           mapping_signature=signature, no_resume=False)
    mapping.write_text('a\tNEW_A\nb\tNEW_B\n')
    source = tmp_path/'gene.mafft.prot.aln.fasta'
    source.write_text('>a\nMA\n>b\nMA\n')
    with pytest.raises(ValueError, match='Mapping changed after'):
        phyloprep.rename_alignment_matrices(args, [source], tmp_path/'out', 'prot', 'raw')
    assert not list((tmp_path/'out').rglob('*.aln.fasta'))
    assert not list((tmp_path/'out').rglob('*.json'))


@pytest.mark.parametrize('suffix', ['NT.fasta', 'AA.fasta', 'NT.fa', 'AA.fa'])
def test_macse_like_source_filename_preserves_stop_provenance(tmp_path, suffix):
    source = tmp_path/f'gene.macse_original.{suffix}'
    alignment = tmp_path/f'{source.stem}.mafft.codon.aln.fasta'
    alignment.write_text('>a\nATGNNNATG\n>b\nATGAAAATG\n')
    report = tmp_path/'changes.tsv'
    write_report(report, [event(source, 'a', 'internal_stop', 1,
                               source_nt_positions='4-6', original='TAA')])
    audit = alignment_qc.load_sequence_changes([report])
    assert alignment_qc.alignment_qc(str(alignment), 'codon', audit_index=audit)['audited_internal_stop']


def test_explicit_macse_intermediate_provenance_still_matches(tmp_path):
    source = tmp_path/'gene.macse_original.NT.fasta'
    alignment = tmp_path/'gene.macse.codon.aln.fasta'
    alignment.write_text('>a\nATGNNNATG\n>b\nATGAAAATG\n')
    report = tmp_path/'changes.tsv'
    write_report(report, [event(source, 'a', 'internal_stop', 1,
                               source_nt_positions='4-6', original='TAA', input_kind='macse_alignment')])
    audit = alignment_qc.load_sequence_changes([report])
    assert alignment_qc.alignment_qc(str(alignment), 'codon', audit_index=audit)['audited_internal_stop']


@pytest.mark.parametrize('flag', ['--max-x-ratio', '--max-n-ratio', '--max-gap-ratio'])
def test_standalone_qc_rejects_nan_before_output(tmp_path, flag):
    source = tmp_path/'input.fa'
    source.write_text('>a\nXXXX\n>b\nXXXX\n')
    with pytest.raises(SystemExit) as error:
        alignment_qc.main(['-i', str(source), '-st', 'prot', '-p', str(tmp_path/'qc'), flag, 'nan'])
    assert error.value.code == 2
    assert not list(tmp_path.glob('qc.*'))


def test_partition_collision_not_misdiagnosed_as_mapping_error(tmp_path, capfd):
    files = []
    for name in ('gene-a', 'gene_a'):
        path = tmp_path/f'{name}.mafft.prot.aln.fasta'
        path.write_text('>a\nMA\n>b\nMA\n')
        files.append(path)
    listing = tmp_path/'inputs.pathlist'
    phyloprep.write_paths(listing, files)
    context = dict(stage=tmp_path/'stage', seqtype='prot', variant='raw', passed=listing)
    with pytest.raises(ValueError, match='Concatenation failed') as error:
        phyloprep.build_group_matrices(SimpleNamespace(mapping=None, missing_taxa='skip-gene'), tmp_path, context)
    assert 'provide phyloprep.py -m' not in str(error.value)
    assert 'Duplicate partition names' in capfd.readouterr().err


@pytest.mark.parametrize('fail_trimmed', [False, True])
def test_browser_publishes_requested_variants_together(tmp_path, monkeypatch, fail_trimmed):
    spec = importlib.util.spec_from_file_location('boundary_browser_runner', ROOT/'site/browser_runner.py')
    browser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(browser)
    source = tmp_path/'source.fa'
    source.write_text('>a\nMA\n>b\nMA\n')
    root = tmp_path/'out'
    old = root/'matrices/prot/raw/supermatrix.fasta'
    old.parent.mkdir(parents=True)
    old.write_text('previous successful matrix')

    def variant(args, root, stage, variant):
        if variant == 'trimmed' and fail_trimmed:
            raise phyloprep.NoGenesIncludedError('No loci retained')
        destination = stage/'prot'/variant
        destination.mkdir(parents=True)
        (destination/'supermatrix.fasta').write_text('>a\nMA\n>b\nMA\n')

    monkeypatch.setattr(phyloprep, 'prepare_variant', variant)
    result = json.loads(browser.run_phyloprep_finalize(
        [dict(role='source', path=str(source), name=source.name)],
        dict(seqtype='prot', notrim=False), str(root), str(ROOT)))
    if fail_trimmed:
        assert not result['ok'] and result['recoverable']
        assert old.read_text() == 'previous successful matrix'
        assert not (root/'matrices/prot/trimmed').exists()
    else:
        assert result['ok']
        assert old.read_text().startswith('>a')
        assert (root/'matrices/prot/trimmed/supermatrix.fasta').is_file()
    assert not list(root.glob('.matrices-*'))


def test_identical_independent_and_paired_sets_reuse_matrix_products(tmp_path, monkeypatch):
    from sequence_audit import normalize_workflow
    root = tmp_path/'results'
    alignments = root/'alignments'
    lists = alignments/'path-lists'
    lists.mkdir(parents=True)
    for kind, sequence in [('codon', 'GCT'), ('prot', 'A')]:
        path = alignments/f'gene.mafft.{kind}.aln.fasta'
        path.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
        phyloprep.write_paths(lists/f'all.{kind}.aln.fasta.pathlist', [path])
    args = phyloprep.build_parser().parse_args(['-i', 'unused.fa', '--notrim', '-t', '1'])
    normalize_workflow(args)
    args.taxon_mapping = None
    calls = []
    build = phyloprep.build_group_matrices
    def record(*arguments):
        calls.append(arguments[2]['stage'])
        return build(*arguments)
    monkeypatch.setattr(phyloprep, 'build_group_matrices', record)
    stage = tmp_path/'stage'
    phyloprep.prepare_variant(args, root, stage, 'raw')
    assert len(calls) == 2
    for kind in ('codon', 'prot'):
        assert (stage/kind/'raw/supermatrix.fasta').read_bytes() == (stage/'paired'/kind/'raw/supermatrix.fasta').read_bytes()
