
import pytest

import bigsmiles


bigsmiles.Config.show_aromatic_bond = False

test_rxn = [
    # no agent
    ["C=CCBr >> C=CCI", [[bigsmiles.BigSMILES("C=CCBr")], [], [bigsmiles.BigSMILES("C=CCI")]]],
    ["CC(=O)[OH]>>CC(=O)OCC", [[bigsmiles.BigSMILES("CC(=O)[OH]")], [], [bigsmiles.BigSMILES("CC(=O)OCC")]]],

    # agent
    ["C1=CC=CC=C1C(O)C(C)NC>SO(Cl)Cl>C1=CC=CC=C1C(Cl)C(C)NC",
     [[bigsmiles.BigSMILES("C1=CC=CC=C1C(O)C(C)NC")],
      [bigsmiles.BigSMILES("SO(Cl)Cl")],
      [bigsmiles.BigSMILES("C1=CC=CC=C1C(Cl)C(C)NC")]
      ]
     ],
    ["C1=CC=CC=C1C(Cl)C(C)NC>[Pd],[H][H]>C1=CC=CC=C1CC(C)NC",
     [[bigsmiles.BigSMILES("C1=CC=CC=C1C(Cl)C(C)NC")],
      [bigsmiles.BigSMILES("[Pd]"), bigsmiles.BigSMILES("[H][H]")],
      [bigsmiles.BigSMILES("C1=CC=CC=C1CC(C)NC")]
      ]
     ],


    # polymer
    ["C=C(C(=O)OC)>CCCCSC(=S)SC(C)C(=O)O,CS(=O)C>{[>][<]CC(C(=O)OC)[>][<]}",  # RAFT of MA in DMSO
     [[bigsmiles.BigSMILES("C=C(C(=O)OC)")],
      [bigsmiles.BigSMILES("CCCCSC(=S)SC(C)C(=O)O"), bigsmiles.BigSMILES("CS(=O)C")],
      [bigsmiles.BigSMILES("{[>][<]CC(C(=O)OC)[>][<]}")]]
     ],
    ["C=Cc1ccccc1.C[CH-](.[Li+])CC>Cc1ccccc1>CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]",
     [
         [bigsmiles.BigSMILES("C=Cc1ccccc1"), bigsmiles.BigSMILES("C[CH-](.[Li+])CC")],
         [bigsmiles.BigSMILES("Cc1ccccc1")],
         [bigsmiles.BigSMILES("CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]")]
     ]
     ],  # anionic of styrene
    ["OCCO,OC(=O)c1ccc(cc1)C(=O)O>>{[][<]OCCO[<],[>]C(=O)c1ccc(cc1)C(=O)[>],[>][H],[<]O[]},O",
     [
         [bigsmiles.BigSMILES("OCCO"), bigsmiles.BigSMILES("OC(=O)c1ccc(cc1)C(=O)O")],
         [],
         [bigsmiles.BigSMILES("{[][<]OCCO[<],[>]C(=O)c1ccc(cc1)C(=O)[>],[>][H],[<]O[]}"), bigsmiles.BigSMILES("O")]
     ]
     ],  # PET synthesis

    # atom index
    ["[OH:1][CH2:2][CH2:3][OH:4],OC(=O)c1ccc(cc1)C(=O)O>>"
     "{[][<][O:1][CH2:2][CH2:3][O:4][<],[>]C(=O)c1ccc(cc1)C(=O)[>],[>][H],[<]O[]},O",
     [
         [bigsmiles.BigSMILES("[OH:1][CH2:2][CH2:3][OH:4]"), bigsmiles.BigSMILES("OC(=O)c1ccc(cc1)C(=O)O")],
         [],
         [bigsmiles.BigSMILES("{[][<][O:1][CH2:2][CH2:3][O:4][<],[>]C(=O)c1ccc(cc1)C(=O)[>],[>][H],[<]O[]}"),
          bigsmiles.BigSMILES("O")]
     ]
     ],  # PET synthesis
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
