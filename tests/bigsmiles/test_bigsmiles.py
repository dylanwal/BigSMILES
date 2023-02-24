import pytest

import bigsmiles
from bigsmiles.errors import BigSMILESError

test_molecules = [

    # daylight tests
    "CCCCCC",  # hexane
    "O=C=O",  # carbon dioxide
    "C#N",  # hydrogen cyanide
    "CCN(CC)CC",  # triethylamine
    "CC(=O)O",  # acetic acid
    "C1CCCCC1",  # cyclohexane
    "c1ccccc1",  # benzene
    # "C1=CN=C[NH]C(=O)1",  # ring number
    # "c1cnc[nH]c(=O)1",  # ring number
    "N[C@@H](C)C(=O)O",  # L-alanine
    "N[C@H](C)C(=O)O",  # D-alanine
    "N[C@](C)(F)C(=O)O",
    "N[C@@](F)(C)C(=O)O",
    "C[C@H]1CCCCO1",
    "O1CCCC[C@@H]1C",
    # "F/C=C/F",  # E-difluoroethene
    # "F\C=C\F",  # Z-difluoroethene


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
    "C1=CC=C(C=C1)S(=O)(=O)NNC2=NC(=NC(=N2)Cl)Cl",
    "C=CC(CCC)C(C(C)C)CCC",  # 3-propyl-4-isopropyl-1-heptene
    "O1CCCCC1N2CCCCC2",

    "C12C3C4C1C5C4C3C25",

    "[12C]C1=C(C(=NC=C1C=O)C)O",
    "[C@]C1=C(C(=NC=C1C=O)C)O",
    "[35Br]CC1=C(C(=NC=C1C=O)C)O",
    "[35Br@@]CC1=C(C(=NC=C1C=O)C)O",
    "[12CH2]C1=C(C(=NC=C1C=O)C)O",
    "[BrH3]",
    "[Fe+3]CCC",
    "CC[12C]C1=C(C(=NC=C1C=O)C)O",
    "CC[C@]C1=C(C(=NC=C1C=O)C)O",
    "CC[35Br]CC1=C(C(=NC=C1C=O)C)O",
    "CC[35Br@@]CC1=C(C(=NC=C1C=O)C)O",
    "CC[12CH2]C1=C(C(=NC=C1C=O)C)O",
    "CC[BrH3]",
    "CC[Fe+3]CCC",

    # cis/trans

    # "C(\F)=C/F",  # trans
    # "F\C=C/F",  # cis
    # "F/C=C\F",  # cis
    # "C(/F)=C/F",  # cis
    # "F/C(CC)=C/F",
    # "F/C=C=C=C/F",  # trans
    # "F/C=C=C=C\F", # cis
    # "F/C=C/C/C=C\C",
    # "F/C=C/CC=CC",  # partially specified

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
    "C{[$][$]CC[$],[$]CC(CC[$])[$][$]}O",
    "CC{[>][<]CC(C)[>][<]}CC(C)=C",  # explicit end groups
    "{[][$]CC[$],[$]CC(CC)[$][]}",  # implicit end groups
    "{[]C([$])C([$])CC[]}",  # test end groups in middle
    "OC{[>][<]CC(C{[>][<]CCO[>][<]}CN)[>][<]}CC",
    "CC(CC){[<][>]CC(C)[<2][>2]}CCO",
    # "{[][>]C([>])([>]),[<]OO[>][>]}CB", Not sure if it a valid BigSMILES

    # From BCPD
    "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][$]}[H]",
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
    ["((CC))"],   # no double branch/ extra parenthesis
    ["C((C)C)"],   # no branch right away
    ["C(C)(C)(C)(C)C"],  # break bond limit of carbon
    ["O1CCCCC1N1CCCCC1"],  # re-use of ring index
    # ['C/C(\F)=C/C'],  # conflicting cis/trans

    # bigsmiles
    ["CCC,C"],  # stochastic seperator
    ["CC}CC"],  # stochastic object no start
    ["CC{CC"],   # stochastic object no end
    ["{CC}"],
    ["{[]CC[]}"],
    # ["{[][$]CC[]}"],  # are these valid?
    # ["CC({[][$]CC[]})CC"],
    # ["CC({[$][$]CC[$]})CC"],
    # ["CC(C{[$][$]CC[$]})CC"],
    ["{[][>]CC[$][]}"],
    ["{[][>]CC[>][]}"],
    ["{[][>]CC[>][]}CC"],  # implicit and explict end groups same time
    ["CC{[<][>]CC[>][]}CC"],  # implicit and explict end groups same time
    # ["{[][>]CC[>];[<]C[]}"],  # only one end group provided
    # ["{[][>]CC[>];[$]C[]}"],
    # ["{[$2][$]CC[$][$2]}"],  # index don't match
    ["{[$][<]CC[>][$]"],  # end group bonding descriptor don't match stochastic fragment
]


@pytest.mark.parametrize("polymer", validation_cases)
def test_validation(polymer: list):
    with pytest.raises(BigSMILESError) as e:
        bigsmiles.BigSMILES(polymer[0])


"""

N1=NN=C[N]1 
[CH2]C=C 
[CH]C=C 
[CH2]CCCC

test for radicals 

"""
