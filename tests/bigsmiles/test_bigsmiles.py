import pytest

import bigsmiles


test_molecules = [
    "CC1=C(C(=NC=C1C=O)C)O",
    "C1=CC(C(C(=C1)C(=O)O)O)O",
    "C(C(=O)COP(=O)(O)O)N",
    "C1(C(C(C(C(C1O)O)OP(=O)(O)O)O)O)O",
    "C1=CC(=C(C=C1Cl)Cl)Cl",
    "CCCCCC(=O)C=CC1C(CC(=O)C1CCCCCCC(=O)O)O",
    "C1=CC(=C(C(=C1)O)O)C(=O)O",
    "CSCCC(=O)C(=COP(=O)(O)O)O",
    "C(C(C(COP(=O)(O)O)O)O)C(=O)C(=O)O",
    "CC(CC1=CC(=C(C=C1)O)O)(C(=O)OC)N",
    "C1=CC=C(C=C1)S(=O)(=O)NNC2=NC(=NC(=N2)Cl)Cl"
]


@pytest.mark.parametrize("molecule", test_molecules)
def test_whole_system_molecule(molecule: str):
    result = bigsmiles.BigSMILES(molecule)
    assert str(result) == molecule


test_polymers = [
    "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]",
    "{[>][$]CC[$],[$]CC(CC)[$][<]}",
    "{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}",
    "{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}",
    "C{[$][$]CC[$],[$]CC(CC)[$][$]}",
    "CC{[>][<]CC(C)[>][<]}CC(C)=C",  # explicit end groups
    "{[][$]CC[$],[$]CC(CC)[$][]}",  # implicit end groups
    "{[]C([$])C([$])CC[]}",  # test end groups in middle

    # From BCPD
    "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][$]}[H]",
    "{[][$]CC(c1cc(C(=O)Oc2ccc(OCCCC)cc2)ccc1(C(=O)Oc3ccc(OCCCC)cc3))[$][$]}{[>][<]Si(C)(C)O[>][]}",
    # "{[][<]CCO[>][<]}{[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][]}",
    # "{[][$]C\C=C/C[$],[$]C\C=C\C[$],[$]CC(C=C)[$][$]}{[>][<][Si](C)(C)O[>][]}",
    # "{[][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}{[$][$]CC(c1ccccc1)[$][]}",
    "{[][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][]}",
    "COCCOCCO{[>][<]CCO[>][<]}{[>][<]C(CC)CO[>][<]}[H]",
    "{[][<]CCO[>][<]}{[>][<]C(CC)CO[>][]}",
    "{[][$]CCCC[$],[$]CC(CC)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][]}",
    # "CCC(C){[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}{[>][<]CCO[>][<]}[H]",
    "{[][<]C(C)C(=O)O[>][<]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][]}",
    # "{[][$]CC(c1ccccc1)[$][$]}{[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][]}",
    # "{[][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}{[$][$]CC(c1ccccc1)[$][]}",
    # "[H]{[>][>]CCO[<][<]}{[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}C(C)CC",
    "{[][$]CC(c1ccccc1)[$],[$]CC(c1ccc(S(=O)(=O)O)cc1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][]}",
    "{[][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][]}",
    "{[][$]CC(C1CCCCC1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][]}",

]


@pytest.mark.parametrize("polymer", test_polymers)
def test_whole_system(polymer: str):
    result = bigsmiles.BigSMILES(polymer)
    assert str(result) == polymer


validation_cases = [
    ["CCCCC1"],  # Ring not closed
    ["CCCC("],  # branch not closed
    ["CCCC)"],  # branch not started
    ["((CC))"],   # no double branch
    ["C((C)C)"],   # no branch right away
    ["C(C)(C)(C)(C)C"],  # break bond limit of carbon

    # bigsmiles
    ["CCC,C"],  # stochastic seperator
    ["CC}CC"],  # stochastic object no start
    ["CC{CC"],   # stochastic object no end
    ["{CC}"],
    ["{[]CC[]}"],
    ["{[][$]CC[]}"],
    ["{[][>]CC[$][]}"],
    ["{[][>]CC[>][]}"],
    ["{[][>]CC[>][]}CC"],  # implicit and explict end groups same time
    ["CC{[<][>]CC[>][]}CC"],  # implicit and explict end groups same time
    ["{[][>]CC[>];[<]C[]}"],  # only one end group provided
    ["{[][>]CC[>];[$]C[]}"],
    ["{[$1][$]CC[$][$1]}"],  # index don't match
    ["{[$][<]CC[>][$]"],  # end group bonding descriptor don't match stochastic fragment
]


@pytest.mark.parametrize("polymer", validation_cases)
def test_validation(polymer: list):
    with pytest.raises(bigsmiles.BigSMILESError) as e:
        bigsmiles.BigSMILES(polymer[0])


"""

N1=NN=C[N]1 
[CH2]C=C 
[CH]C=C 
[CH2]CCCC

test for radicals 

"""