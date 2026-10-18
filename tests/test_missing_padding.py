"""Missing loci stay distinct from gaps throughout concatenation/export."""
import csv
import importlib.util
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest
from Bio import SeqIO

from alignment_qc import alignment_qc, main as qc_main
from fasta2axt import fasta2axt
from fasta2nex import fasta2nex
from fasta2phy import fasta2phy
from get_supergenes import get_supergenes
from sequence_audit import normalize_nucleotides
from split_codon_seqence_alignment import split

spec = importlib.util.spec_from_file_location('fourfold', Path(__file__).parents[1] / 'extract_4-fold_degenerated_sites.py')
fourfold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fourfold)


def sequences(path, fmt='fasta'):
    return {record.id: str(record.seq) for record in SeqIO.parse(path, fmt)}


@pytest.mark.parametrize('seqtype', ['DNA', 'protein'])
def test_export_preserves_missing_without_masking_or_auditing(tmp_path, seqtype):
    source = tmp_path / 'in.fa'
    sequence = 'AXN-???' if seqtype == 'protein' else 'ACN-???'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    expected = {'a': sequence, 'b': sequence}
    for converter, suffix, fmt in [(fasta2phy, 'phy', 'phylip-relaxed'), (fasta2nex, 'nex', 'nexus')]:
        output = tmp_path / ('matrix.' + suffix)
        converter(source, output, seqtype)
        assert sequences(output, fmt) == expected
        assert not Path(str(output) + '.sequence_changes.tsv').exists()
    axt = tmp_path / 'matrix.axt'
    fasta2axt(source, axt, seqtype)
    assert axt.read_text().splitlines()[1:] == [sequence, sequence]
    assert not Path(str(axt) + '.sequence_changes.tsv').exists()


def test_padding_coordinates_and_codon_derivatives(tmp_path, monkeypatch):
    complete, partial = tmp_path / 'complete.fa', tmp_path / 'partial.fa'
    complete.write_text('>a\nGCT---GCC\n>b\nGCTGCTGCC\n')
    partial.write_text('>a\nGCTGCC\n')
    prefix = str(tmp_path / 'matrix')
    get_supergenes([str(complete), str(partial)], prefix, missing_taxa='pad-gaps')
    source = Path(prefix + '.fasta')
    assert sequences(source) == {'a': 'GCT---GCCGCTGCC', 'b': 'GCTGCTGCC??????'}
    assert sequences(Path(prefix + '.phy'), 'phylip-relaxed') == sequences(source)
    with Path(prefix + '.missing.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert rows == [dict(locus='partial_fa', taxon='b', start='10', end='15', length='6', fill_character='?')]
    monkeypatch.chdir(tmp_path)
    split(source)
    for position, expected in [('1st', 'G-GGG'), ('2nd', 'C-CCC'), ('3rd', 'T-CTC')]:
        output = tmp_path / f'codon{position}.matrix.fasta'
        assert sequences(output)['a'] == expected
        assert sequences(output)['b'].endswith('??')
        for converter, suffix, fmt in [(fasta2phy, 'phy', 'phylip-relaxed'), (fasta2nex, 'nex', 'nexus')]:
            converted = output.with_suffix('.' + suffix)
            converter(output, converted)
            assert sequences(converted, fmt) == sequences(output)
    extracted = tmp_path / 'fourfold.fa'
    fourfold.extract_four_fold_site(source, extracted)
    # The complete, ungapped GC-family columns pass; absent taxa still exclude a site.
    assert sequences(extracted) == {'a': 'TC', 'b': 'TC'}


def test_empty_missing_table_and_existing_unknowns_are_not_padding(tmp_path):
    source = tmp_path / 'gene.fa'
    source.write_text('>a\nAC???-\n>b\nAC???-\n')
    prefix = str(tmp_path / 'matrix')
    get_supergenes([str(source)], prefix, missing_taxa='pad-gaps')
    with Path(prefix + '.missing.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == []
    assert sequences(Path(prefix + '.fasta')) == sequences(source)
    assert normalize_nucleotides('U?R') == 'TNN'
    assert normalize_nucleotides('U?R', preserve_missing=True) == 'T?N'


@pytest.mark.parametrize('mode,sequence', [('nucl', 'ACN-??'), ('prot', 'ANX-??'), ('codon', 'GCT---NNN???')])
def test_qc_reports_missing_separately_from_gap_and_ambiguity(tmp_path, mode, sequence):
    source = tmp_path / 'matrix.fa'
    source.write_text(f'>a\n{sequence}\n>b\n{sequence}\n')
    result = alignment_qc(str(source), mode)
    stats = result['alignment_stats']
    assert stats['missing_ratio'] == sequence.count('?') / len(sequence)
    assert stats['gap_ratio'] == sequence.count('-') / len(sequence)
    assert result['sequences'][0]['missing_count'] == sequence.count('?')
    assert 'invalid_characters' not in result['issues']


def test_qc_details_exposes_separate_missing_ratio(tmp_path):
    source, prefix = tmp_path / 'matrix.fa', tmp_path / 'qc'
    source.write_text('>a\nGCT---NNN???\n>b\nGCTGCCGCTGCT\n')
    qc_main(['-i', str(source), '-st', 'codon', '-p', str(prefix)])
    with Path(str(prefix) + '.details.tsv').open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert float(row['Missing_ratio']) == 0.125
    assert float(row['Gap_ratio']) == 0.125
    assert float(row['N_ratio']) == 0.125


CODEML = os.environ.get('PHYLOPREP_TEST_CODEML') or shutil.which('codeml')


@pytest.mark.skipif(not CODEML, reason='Optional PAML executable not installed')
@pytest.mark.parametrize('cleandata', [0, 1])
def test_codeml_missing_and_gap_padding_have_equal_likelihood(tmp_path, cleandata):
    data = {
        'a': 'ATGGCTGCCGCAAAGAAAGGTGGCTTTTTCCTGCTACCTCCGTCAGACTGACCTGGT',
        'b': 'ATGGCCGCCGCAAAGAAAGGTGGTTTTTTCCTGCTACCTCCGTCAGACTGACCTGGT',
        'c': 'ATGGCTGCCGCAAAAAAAGGTGGCTTTTTCCTGCTACCTCCGTCAGACTGACCTGGC',
        'd': 'ATGGCTGCCGCAAAGAAAGGTGGCTTTTTCCTGCTACCTCCGTCAGACTGACCTGGA',
    }
    likelihoods = []
    for name, marker in [('gap', '-'), ('missing', '?')]:
        directory = tmp_path / name
        directory.mkdir()
        sequences = {**data, 'd': data['d'][:-9] + marker * 9}
        source = directory / 'input.fa'
        source.write_text(''.join(f'>{taxon}\n{sequence}\n' for taxon, sequence in sequences.items()))
        fasta2phy(source, directory / 'input.phy')
        (directory / 'tree.nwk').write_text('4 1\n((a:0.1,b:0.1):0.1,c:0.1,d:0.1);\n')
        (directory / 'codeml.ctl').write_text(
            'seqfile = input.phy\ntreefile = tree.nwk\noutfile = result.txt\n'
            'noisy = 0\nverbose = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\n'
            'clock = 0\nmodel = 0\nNSsites = 0\nfix_kappa = 1\nkappa = 2\n'
            f'fix_omega = 1\nomega = 0.5\nfix_blength = 2\ncleandata = {cleandata}\n')
        result = subprocess.run([str(CODEML)], cwd=directory, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        match = re.search(r'lnL\([^\n]*?:\s*(-?\d+\.\d+)', (directory / 'result.txt').read_text())
        assert match, result.stdout
        likelihoods.append(float(match.group(1)))
    assert likelihoods[0] == likelihoods[1]
