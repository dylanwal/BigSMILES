import pytest

import bigsmiles


test_molecules = [
"CC1=C(C(=NC=C1C=O)C)O"
]


@pytest.mark.parametrize("molecule", test_molecules)
def test_whole_system_molecule(molecule: str):
    result = bigsmiles.BigSMILES(molecule)
    assert str(result) == molecule


test_polymers = [
    "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
]


@pytest.mark.parametrize("polymer", test_polymers)
def test_whole_system(polymer: str):
    result = bigsmiles.BigSMILES(polymer)
    assert str(result) == polymer
