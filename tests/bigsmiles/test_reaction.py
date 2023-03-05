import logging

import pytest

import bigsmiles
from bigsmiles.errors import BigSMILESError

bigsmiles.Config.show_aromatic_bond = False

test_rxn = [
    # no agent
    ["C=CCBr >> C=CCI", [[bigsmiles.BigSMILES("C=CCBr")], [], [bigsmiles.BigSMILES("C=CCI")]]],
    ["CC(=O)[OH]>>CC(=O)OCC", [[bigsmiles.BigSMILES("CC(=O)[OH]")], [], [bigsmiles.BigSMILES("CC(=O)OCC")]]],

    # agent
    ["CC(=O)[OH]>>CC(=O)OCC", [[bigsmiles.BigSMILES("CC(=O)[OH]")], [], [bigsmiles.BigSMILES("CC(=O)OCC")]]],

    # polymer
    ["C=C(C(=O)OC)>CCCCSC(=S)SC(C)C(=O)O,CS(=O)C>{[>][<]CC(C(=O)OC)[>][<]}",  # RAFT of MA in DMSO
     [[bigsmiles.BigSMILES("C=C(C(=O)OC)")],
      [bigsmiles.BigSMILES("CCCCSC(=S)SC(C)C(=O)O"), bigsmiles.BigSMILES("CS(=O)C")],
      [bigsmiles.BigSMILES("{[>][<]CC(C(=O)OC)[>][<]}")]]
     ],
]


@pytest.mark.parametrize("case", test_rxn)
def test_whole_system_molecule(case: str):
    input_, answer = case
    result = bigsmiles.Reaction(input_)
    # reactants
    for chem, ans in zip(result.reactants, answer[0]):
        assert chem == ans
    # agents
    for chem, ans in zip(result.agents, answer[1]):
        assert chem == ans
    # products
    for chem, ans in zip(result.products, answer[2]):
        assert chem == ans
