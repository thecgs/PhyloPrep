import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(tmp_path, script, *args):
    return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)], cwd=tmp_path,
                          text=True, capture_output=True)


def table(path):
    with path.open() as handle:
        return list(csv.reader(handle, delimiter='\t'))


def test_proteinortho_converter_writes_standard_table(tmp_path):
    source, output = tmp_path / 'groups.proteinortho', tmp_path / 'groups.tsv'
    source.write_text('# Species\tGenes\tAlg.-Conn.\tHuman.faa\tMouse.fa\n2\t3\t0.5\tH1,H2\tM1\n')
    result = run(tmp_path, 'proteinortho_to_orthogroups.py', '-i', source, '-o', output)
    assert result.returncode == 0, result.stderr
    assert table(output) == [['Orthogroup', 'Human', 'Mouse'], ['POG0000001', 'H1, H2', 'M1']]


def test_sonicparanoid_converter_writes_standard_table(tmp_path):
    source, output = tmp_path / 'flat.ortholog_groups.tsv', tmp_path / 'groups.tsv'
    source.write_text('group_id\tHuman.faa\tMouse.pep\nSP0001\tH1\tM1,M2\n')
    result = run(tmp_path, 'sonicparanoid_to_orthogroups.py', '-i', source, '-o', output)
    assert result.returncode == 0, result.stderr
    assert table(output) == [['Orthogroup', 'Human', 'Mouse'], ['SP0001', 'H1', 'M1, M2']]


def test_orthomcl_and_oma_converters_accept_their_group_layouts(tmp_path):
    orthomcl, oma = tmp_path / 'groups.txt', tmp_path / 'OrthologousGroups.txt'
    orthomcl_out, oma_out = tmp_path / 'orthomcl.tsv', tmp_path / 'oma.tsv'
    orthomcl.write_text('ORTHOMCL1(3 genes,2 taxa): Human|H1 Mouse|M1 Mouse|M2\n')
    oma.write_text('OMA0001\tHuman|H1\tMouse|M1\n')
    result = run(tmp_path, 'orthomcl_to_orthogroups.py', '-i', orthomcl, '-o', orthomcl_out,
                 '--strip-taxon-prefix')
    assert result.returncode == 0, result.stderr
    assert table(orthomcl_out) == [['Orthogroup', 'Human', 'Mouse'], ['ORTHOMCL1(3 genes,2 taxa)', 'H1', 'M1, M2']]
    result = run(tmp_path, 'oma_to_orthogroups.py', '-i', oma, '-o', oma_out, '--strip-taxon-prefix')
    assert result.returncode == 0, result.stderr
    assert table(oma_out) == [['Orthogroup', 'Human', 'Mouse'], ['OMA0001', 'H1', 'M1']]
