import logging

import pytest

import bigsmiles
from bigsmiles.errors import BigSMILESError

bigsmiles.Config.show_aromatic_bond = False


test_molecules = [
    # simple
    "C",
    "U",
    "Zn",
    "CC",
    "CCCCCC",  # hexane
    "O=C=O",  # carbon dioxide
    "CC#CC",  # 2-butyne
    "C#N",  # hydrogen cyanide
    "CCN(CC)CC",  # triethylamine
    "CC(=O)O",  # acetic acid
    "C1CCCCC1",  # cyclohexane

    # rings
    "C=1CCCCC1",
    "C12C3C4C1C5C4C3C25",
    "O1CCCCC1N2CCCCC2",
    "C12(CCCCC1)CCCCC2",  # spiro[5.5]undecane

    # aromatic
    "c1ccccc1",  # benzene
    "n1ccccc1",  # pyridine
    "Cn1cccc1",  # Methyl-pyrrole
    "[nH]1cccc1",  # 1H-pyrrole
    "O=c1[nH]cccc1",  # 2-pyridone
    "Oc1ncccc1",  # 2-pyridinol
    "C1=CC=CC(CCC2)=C12",  # indane
    "c1ccc2CCCc2c1",  # indane
    "C1=CC=PC=C1",  # Phosphorine
    "c1ccp(=O)cc1",  # Phosphorine oxide
    "O=P1=CC=CC=C1",  # Phosphorine oxide
    "O=p1ccccc1",  # Phosphorine oxide
    "c1ccccc1-c2ccccc2",
    "C1(C=C2)=C3C2=CC=C3C=C1",  # Acepentalene


    # cis/trans
    # "F/C=C/F",  # E-difluoroethene
    # "F\C=C\F",  # Z-difluoroethene
    # "C(\F)=C/F",  # trans
    # "F\C=C/F",  # cis
    # "F/C=C\F",  # cis
    # "C(/F)=C/F",  # cis
    # "F/C(CC)=C/F",
    # "F/C=C=C=C/F",  # trans
    # "F/C=C=C=C\F", # cis
    # "F/C=C/C/C=C\C",
    # "F/C=C/CC=CC",  # partially specified

    # charges
    "[Fe+]",
    "[Fe+2]",
    "[OH-]",
    "[OH3+]",
    "[O-][n+]1ccccc1",  # Pyridine-N-oxide
    "[O-][N+]1=CC=CC=C1",

    # chiral
    "N[C@@H](C)C(=O)O",  # L-alanine
    "N[C@H](C)C(=O)O",  # D-alanine
    "N[C@](C)(F)C(=O)O",
    "N[C@@](F)(C)C(=O)O",
    "C[C@H]1CCCCO1",
    "O1CCCC[C@@H]1C",
    "N[P@](C(O)=O)C",

    # disconnected structures
    "C1.CCCCC1",  # hexane
    "[Na+].[Cl-]",  # sodium chloride
    "[Na+].[O-]c1ccccc1",  # sodium phenoxide
    "c1cc([Na+].[O-])ccc1",  # sodium phenoxide
    "[NH4+].[NH4+].[O-]S(=O)(=O)[S-]",  # diammonium thiosulfate
    "c1c2c3c4cc1.Br2.Cl3.Cl4",  # 1-bromo-2,3-dichlorobenzene

    # atom indexing
    "[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH3:6]",  # hexane
    "[CH2:1]1[CH2:2][CH:3]=[CH:4][CH2:5][CH2:6]1",  # cylcohexene
    "[CH3:1]CC[CH2:2][C:3]",  # hexane
    "[CH2:1]=[CH:2][CH2:1][CH2:3][CH:4](C)[CH3:3]",
    "[OH:1][C:2](=[O:3])[CH:4]([CH3:5])[OH:6]",  # lactic acid

    # metal complex
    "[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl",

    # extended atom
    "[C]",
    "[Fe]",
    "[12C]",
    "[238U]",
    "[13CH4]",
    "[12C@@H2+]",
    "[38Br-]",
    "[34SH2+2]",
    "[CH4]",
    "[H][CH2][H]",
    "[H][CH]([H])[H]",
    "[ClH]",
    "[2H][CH2]C",

    # complex / random
    "OS(=O)(=S)O",
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
    "[12C]C1=C(C(=NC=C1C=O)C)O",
    "[C@]C1=C(C(=NC=C1C=O)C)O",
    "[35Br]CC1=C(C(=NC=C1C=O)C)O",
    "[35Br@@]CC1=C(C(=NC=C1C=O)C)O",
    "[12CH2]C1=C(C(=NC=C1C=O)C)O",
    "[BrH3]",
    "[Fe+3]CCC",
    "CC[12C]C1=C(C(=NC=C1C=O)C)O",
    "CC[C@]C1=C(C(=NC=C1C=O)C)O",
    "CC[12CH2]C1=C(C(=NC=C1C=O)C)O",
    "CC[BrH3]",
    "CC[Fe+3]CCC",

]


@pytest.mark.parametrize("molecule", test_molecules)
def test_whole_system_molecule(molecule: str):
    result = bigsmiles.BigSMILES(molecule)
    assert str(result) == molecule


test_polymers = [
    # homopolymer
    "[H]{[$][$]CC[$][$]}[H]",  # poly(ethylene)
    "[H]{[>][<]CC[>][<]}[H]",
    "[H]{[$][$]CC(C)[$][$]}[H]",  # poly(propylene) random orientation
    "[H]{[>][<]CC(C)[>][<]}[H]",
    "[H]{[>][<]CC(C1=CC=CC=C1)[>][<]}[H]",  # PS head-to-tail
    "[H]{[>][<]CC(C1)=CC=CC=C1[>][<]}[H]",
    "[H]{[>][<]NCCCCCC(=O)[>][<]}O",  # nylon 6
    "[H]{[>][<]OC(C)C(=O)[>][<]}O",  # PLA
    "O{[>][<]C(=O)C(C)O[>][<]}[H]",
    "[H]{[>][<]CCO[>][<]}[H]",  # polyethylene oxide
    "[H]{[>][<]CC(C)O[>][<]}[H]",  # polypropylene oxide
    "[H]{[>][<]CC(C#N)[>][<]}[H]",  # polyacrylonitrile
    "[H]{[>][<]CC=CC[>][<]}[H]",  # polybutadiene
    # r"[H]{[<][>]CC\=\CC[<][>]}[H]",  # polybutadiene (trans)
    # r"[H]{[<][>]CC\=/CC[<][>]}[H]",  # polybutadiene (cis)
    # with branches
    "[H]{[<][$]CC(C)(CC)[$][>]}[H]",
    "[H]{[<][$]CC(CCC)[$][>]}[H]",
    "[H]{[<][$]CC(C(CC)C)[$][>]}[H]",
    "CC{[>][<]CC(C)[>][<]}CC(C)=C",
    # with rings
    "[H]{[<][$]C1CC(C1)[$][>]}[H]",  # 4 membered rings
    "[H]{[<][$]C1CCC(C1)[$][>]}[H]",  # 5 membered rings (3, 2)
    "[H]{[<][$]C1CC(CC1)[$][>]}[H]",  # 5 membered rings
    "[H]{[<][$]C1(CCCCC1)[$][>]}[H]",  # 6 membered rings
    # multiple bonds
    "C={[$][$]=CC=[$][$]}=C",  # polyacetylene (traditional)
    "C={[$][$]=C[$2][$2]}[H]",  # polyacetylene (minimal)
    "C#{[$][$]#CC1=CC=C(C=C1)C#[$][$]}#C",

    # copolymers
    "C{[$][$]CC[$],[$]CC(CC[$])[$][$]}O",

    # implicit end groups
    "{[][$]CC[$],[$]CC(CC)[$][]}",
    "{[][$]=CC=[$][]}",

    # mixed end groups
    "{[][>]C([>])[>],[<]C[>][>]}C",
    "[H]{[<][>]OCCO[>],[<]C(=O)C1=CC=C(C=C1)C(=O)[<],[>]O[]}",  # PET
    "{[][>]NCCCCCCN[>],[<]C(=O)CCCCC(=O)[<],[>]Cl,[<][H][]}",  # nylon 6,6

    # test end groups in middle
    "{[]C([$])C([$])CC[]}",
    "[H]{[$]C([$])C([$])CC[$]}[H]",

    # nested linear
    "OC{[>][<]CC(C{[>][<]CCO[>][<]}CN)[>][<]}CC",

    # bonding descriptor index
    "CC(CC){[<][>]CC(C)[<2][>2]}CCO",

    # alternating
    "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]",
    "[H]{[>][<]CC(C1=CC=CC=C1)[<],[>]C2C(C(=O)OC2=O)[>][<]}[H]",  # PS-alt-poly(maleic anhydride)

    # block copolymer
    "[H]{[$][$]CC[$][$]}{[$][$]CC[$][$]}[H]",
    '[H]{[$][$]CC[$][$]}{[$][$]CC(C)[$][$]}[H]',  # poly(ethylene)-block-poly(propylene)
    "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][$]}[H]",
    '[H]{[>][<]OCC[>][<]}O{[>][<]C(=O)C(C)O[>][<]}[H]',  # PEG-Oxygen-PLA
    'O{[>][<]CCO[>][<]}{[>][<]C(=O)C(C)O[>][<]}[H]',  # Oxygen - PEG-block-PLA
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
    "{[][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][]}",
    "{[][$]CC(C1CCCCC1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C2CCCCC2)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C3CCCCC3)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C4CCCCC4)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C5CCCCC5)[$][$]}"
    "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C6CCCCC6)[$][]}",
    # multi bonds
    "C={[$][$]=CC=[$][$]}={[$][$]=CC=[$][$]}=C",
    "C={[$][$]=CCCCCCCC=[$],[$]=CC1CCC(C1)C=[$][$]}=C",  # poly(cyclooctene)-rand-poly(noroborene)

    # star
    'CC({[$][$]CC[$][$]}[H])CCC',
    'CC({[$][$]CC[$][$]}[H])({[$][$]CC[$][$]}[H])CCC',
    'CC({[$][$]CC[$][$]}[H])({[$][$]CC[$][$]}O)CCC',
    'CC({[$][$]CC[$][$]}N)({[$][$]CC[$][$]}O)CCC',

    # # graft/bottlebrush
    "[H]{[>][<]CC([>2])[>],[<2]CC[>2],[<2][H][<]}[H]",  # bottlebrush
    "[H]{[>][<]CC([>2])[>],[<]CC(C)[>],[<2]CC[>2],[<2][H][<]}[H]",  # graft copolymerization
    "[H]{[>][<]CC([>2])[>],[<2]CC[>2],[<2][>3],[<3]CC(C)[>3],[<3][H][<]}[H]",  # graft diblock brushes
    "C={[>][<]=CC([>2])=[>],[<2]CC[>2],[<2][>3],[<3]CC(C)[>3],[<3][H][<]}=C",  # graft diblock brushes (double bond)
    "C={[>][<]=CC(=[>2])=[>],[<2]=CC=[>2],[<2]=[>3],[<3]=CC(C)=[>3],[<3]=C[<]}=C",  # graft diblock brushes(double bond)
    "[H]{[>][<]CC([>2])[>],[<2]CC(C)[>2],[<2]CC[>2],[<2][H][<]}[H]",  # graft random copolymer brushes
    # nesting
    "[H]{[>][<]CC({[>][<]CC[>][<]}[H])[>][<]}[H]",  # bottlebrush

    # hyperbranched
    "[H]{[$][$]CC([$])[$][$]}[H]",  # hyperbranch Poly(ethylene)
    "[H]{[$][$]CC([$])[$],[$]CC[$][$]}[H]",  # hyperbranch Poly(ethylene)

    # rings/polycyclic (stochastic)
    "C1{[$][$]CC[$][$]}1",
    "C({[>][<]CCO[>][<]}1)CCCCCC1",
    "C1{[>][<]CCO[>][<]}CO1",  # PEG rings
    "C({[>][<]CCO[>][<]}1)CCO1",  # PEG rings wierd written
    "C({[>][<]CCO[>][<]}1)({[>][<]CCO[>][<]}2)OCC12",  # dual PEG rings
    "C({[>][<]CCO[>][<]}1)({[>][<]CCO[>][<]}2)({[>][<]CCO[>][<]}3)OCC123",  # triple PEG rings
    # multiple bonds
    "C(={[>][<]=CC=[>][<]}=1)CCCCCC=1",
    "C=1CCC(={[>][<]=CC=[>][<]}=1)CCC",
    # aromatic notation
    "c1c{[$][$]cc[$][$]}1",  # cyclic polyacetylene

    # networks
    "{[][$]CC=CC[$],[$]CC([<])C([>])C[$],[>]S[<],[$]C[]}",  # poly(1,4-butadiene) rubber vulcanized
    "[H]{[>][<]CCN([>])[<][<]}[H]",  # branched polyethylenimine

    # ladder


]


@pytest.mark.parametrize("polymer", test_polymers)
def test_whole_system(polymer: str):
    result = bigsmiles.BigSMILES(polymer)
    assert str(result) == polymer


cases_with_changes = [
    # hydrogens
    ["[ClH1]", "[ClH]"],

    # charges
    ["[OH-1]", "[OH-]"],
    ["[Fe++]", "[Fe+2]"],
    ["[12C@@H2+1]", "[12C@@H2+]"],
    ["[O-1][n+1]1ccccc1", "[O-][n+]1ccccc1"],

    # branching
    ["C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C",
     "C(CCCCCCCCCCCCCCCCCCCC)C"],

    # rings
    ["C=1CCCCC=1", "C=1CCCCC1"],
    ["C1CCCCC=1", "C=1CCCCC1"],
    ["C1=CN=C[NH]C(=O)1", "C1=CN=C[NH]C1=O"],  # ring number position change
    ["c1cnc[nH]c(=O)1", "c1cnc[nH]c1=O"],  # ring number position change
    ["C(={[>][<]=CC=[>][<]}12)CCCCCC12", "C(={[>][<]=CC=[>][<]}=1)CCCCCC=1"],
    ["C2(C=C3)=C1C3=CC=C1C=C2",  'C1(C=C2)=C3C2=CC=C3C=C1'],  # Acepentalene

    # polymers
    ["{[][>]C([>])([>]),[<]C[>][>]}C", "{[][>]C([>])[>],[<]C[>][>]}C"],  # drop () around bond descriptor
    ["[H]{[>][<]CC(C1=CC=CC=C1)[<],[>]C2C(C(=O)OC2(=O))[>][<]}[H]",
     "[H]{[>][<]CC(C1=CC=CC=C1)[<],[>]C2C(C(=O)OC2=O)[>][<]}[H]"],  # drop () around =O
    ["[H]{[>1][<]CC([>2])[>],[<2]CC[>2],[<2][H][<1]}[H]", "[H]{[>][<]CC([>2])[>],[<2]CC[>2],[<2][H][<]}[H]"],  # drop 1
    ["C(={[>][<]=CC=[>][<]}=1)CCCCCC1", "C(={[>][<]=CC=[>][<]}=1)CCCCCC=1"],  # add double bond to second ring
    ["C(={[>][<]=CC=[>][<]}1)CCCCCC=1", "C(={[>][<]=CC=[>][<]}=1)CCCCCC=1"],  # add double bond to second ring
    ["C1CCC(={[>][<]=CC=[>][<]}=1)CCC", "C=1CCC(={[>][<]=CC=[>][<]}=1)CCC"],  # add double bond to second ring
    ["C=1CCC(={[>][<]=CC=[>][<]}1)CCC", "C=1CCC(={[>][<]=CC=[>][<]}=1)CCC"],  # add double bond to second ring

    # random
    ["C1.C1", "CC"],  # disconnect and ring notation cancel out
]


@pytest.mark.parametrize("case", cases_with_changes)
def test_whole_system_with_changes(case: list):
    input_, answer = case
    result = bigsmiles.BigSMILES(input_)
    assert str(result) == answer


ring_fix_cases = [
    ["O1CCCCC1N1CCCCC1", "O1CCCCC1N2CCCCC2"],  # re-use of ring index],  # re-use of ring index
    ["[H]{[>][<]CC(C1=CC=CC=C1)[<],[>]C1C(C(=O)OC1=O)[>][<]}[H]",
     "[H]{[>][<]CC(C1=CC=CC=C1)[<],[>]C2C(C(=O)OC2=O)[>][<]}[H]"],
    ["{[][$]CC(c1ccccc1)[$],[$]CC(c1ccc(S(=O)(=O)O)cc1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][]}",
     "{[][$]CC(c1ccccc1)[$],[$]CC(c2ccc(S(=O)(=O)O)cc2)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][]}"],
    ["{[][$]CC(C1CCCCC1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C2CCCCC2)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][]}",
     "{[][$]CC(C1CCCCC1)[$][$]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C2CCCCC2)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C3CCCCC3)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C4CCCCC4)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C5CCCCC5)[$][$]}"
     "{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[$][$]CC(C6CCCCC6)[$][]}"
     ]

]


@pytest.mark.parametrize("case", ring_fix_cases)
def test_ring_fix(caplog, case: list):
    input_, answer = case
    with caplog.at_level(logging.WARNING):
        result = bigsmiles.BigSMILES(input_)
    assert str(result) == answer
    assert "Duplicate ring index detected" in caplog.text


cases_add_explicit_hydrogens = [
    ["{[$][$]CC[$][$]}", "[H]{[$][$]CC[$][$]}[H]"],
    ["O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}", "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"],
    ["{[$][$]CC[$],[$]CC(CC)[$][$]}", "[H]{[$][$]CC[$],[$]CC(CC)[$][$]}[H]"],
    ["{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}", "[H]{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"],
    ["{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}", "[H]{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}[H]"],
    ["C{[$][$]CC[$],[$]CC(CC)[$][$]}", "C{[$][$]CC[$],[$]CC(CC)[$][$]}[H]"],
]


@pytest.mark.parametrize("case", cases_add_explicit_hydrogens)
def test_add_explicit_hydrogens(case: list):
    input_, answer = case
    result = bigsmiles.BigSMILES(input_)
    assert str(result) == answer


cases_incomplete_valance = [
    ["N1=NN=C[N]1", "N1=NN=C[N]1"],  # under-saturated nitrogen
    ["CCC[12C]CCC", "CCC[12C]CCC"],  # under-saturated carbon
    ["[CH2]CCCC", "[CH2]CCCC"],  # under-saturated carbon
    ["[CH]C=C", "[CH]C=C"],  # under-saturated carbon
    ["[CH2]C=C", "[CH2]C=C"],  # under-saturated carbon
    ["CC[35Br]CC1=C(C(=NC=C1C=O)C)O", "CC[35Br]CC1=C(C(=NC=C1C=O)C)O"],  # over-saturated Br
    ["CC[35Br@@]CC1=C(C(=NC=C1C=O)C)O", "CC[35Br@@]CC1=C(C(=NC=C1C=O)C)O"],  # over-saturated Br
    ["C(C)(C)(C)(C)C", "C(C)(C)(C)(C)C"],  # over-saturated C
    ["O=n1ccccc1", "O=n1ccccc1"],  # Pyridine-N-oxide
]


@pytest.mark.parametrize("case", cases_incomplete_valance)
def test_incomplete_valance(caplog, case: list):
    input_, answer = case
    with caplog.at_level(logging.WARNING):
        result = bigsmiles.BigSMILES(input_)
    assert str(result) == answer
    assert "Incomplete valence detected" in caplog.text


validation_cases = [
    ["DJW"],  # not an element
    ["[C"],
    ["[C]]"],
    ["[]"],
    ["[.]"],
    ["[rings=]"],
    ["CCCCC1"],  # Ring not closed
    ["CCCC("],  # branch not closed
    ["CCCC)"],  # branch not started
    ["((CC))"],  # no double branch/ extra parenthesis
    ["C((C)C)"],  # no branch right away

    # ['C/C(\F)=C/C'],  # conflicting cis/trans

    # bigsmiles
    ["CCC,C"],  # stochastic seperator
    ["CC}CC"],  # stochastic object no start
    ["CC{CC"],  # stochastic object no end
    ["{CC}"],
    ["{[]CC[]}"],
    # ["{[][$]CC[]}"],  # are these valid?
    ["CC({[][$]CC[]})CC"],
    # ["CC({[$][$]CC[$]})CC"],
    # ["CC(C{[$][$]CC[$]})CC"],
    # ["{[>][$]CC[$],[$]CC(CC)[$][<]}"],
    ["{[][>]CC[$][]}"],
    ["{[][>]CC[>][]}"],
    ["{[][>]CC[>][]}CC"],  # implicit and explict end groups same time
    ["CC{[<][>]CC[>][]}CC"],  # implicit and explict end groups same time
    # ["{[][>]CC[>];[<]C[]}"],  # only one end group provided
    # ["{[][>]CC[>];[$]C[]}"],
    # ["{[$2][$]CC[$][$2]}"],  # index don't match
    ["{[$][<]CC[>][$]"],  # end group bonding descriptor don't match stochastic fragment
    ["C={[$][$]=C[$][$]}[H]"],  # multiple bond types to [$]
    ["C{[>][<]=CC([>2])=[>],[<2]CC[>2],[<2][>3],[<3]CC(C)[>3],[<3][H][<]}C"],  # multiple bond types to [>]
    ["{[$][$]=CC=[$][$]}"],  # end group must have double bond
    ["C={[$][$]=CC=[$][$]}"],  # end group must have double bond
    ["{[$][$]=CC=[$][$]}=C"],  # end group must have double bond
    ["1{[$][$]CC[$][$]}1"],  # ring id can't start
    ["={[$][$]CC[$][$]}1"],  # ring id can't start
    ["[$]{[$][$]CC[$][$]}1"],  # ring id can't start
    ["C1{[$][$]CC[$][$]}1C"],  # last carbon has nothing to bond to
    ["C{[>][<]=CC=([>2])=[>],[<2]CC[>2],[<2][>3],[<3]CC(C)[>3],[<3][H][<]}C"],  # first repeat extra =
    ["C{[>][<]=CC=([>2])[>],[<2]CC[>2],[<2][>3],[<3]CC(C)[>3],[<3][H][<]}C"],  # first repeat wrong = location

    # ["C{[$][$]C.C[$][$]}C"]

]


@pytest.mark.parametrize("polymer", validation_cases)
def test_validation(polymer: list):
    with pytest.raises(BigSMILESError) as _:
        bigsmiles.BigSMILES(polymer[0])


# case_for_conjugated = [
#     ["c1ccccc1", {0: True, 1: True, 2: True, 3: True, 4: True, 5: True}],  # benzene
#     ["c1cnccc1", {0: True, 1: True, 2: True, 3: True, 4: True, 5: True}],  # pyridine
#     ["c1cocc1", {0: True, 1: True, 2: True, 3: True, 4: True, 5: True}],  # furan
#     ["CCc1ccccc1", {0: False, 1: False, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True}],  # ethyl benzene
#     ["C=Cc1ccccc1", {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True}],  # styrene
#     ["C=CC=C", {0: True, 1: True, 2: True}],  # butadiene
#     ["C=CC(C)=C", {0: True, 1: True, 2: False, 3: True}],  # isoprene
#     ["C=COC=C", {0: True, 1: True, 2: True, 3: True}],
#     ["CC=COC=C", {0: False, 1: True, 2: True, 3: True, 4: True}],
#     ["C=C[C+]C=C", {0: True, 1: True, 2: True, 3: True}],
#     ["C=C[C-]C=C", {0: True, 1: True, 2: True, 3: True}],
#     ["c1ccc(cc1)/C=C/C=O",  {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True}]
# c1ccc2cc3ccccc3cc2c1  # anthacene
#
#     # negative
#     ["C=CCC=C", {0: False, 1: False, 2: False, 3: False}],
# ]
#
#
# @pytest.mark.parametrize("case", case_for_conjugated)
# def test_incomplete_conjugated(case: list):
#     input_, answer = case
#     result = bigsmiles.BigSMILES(input_)
#     for k, v in answer.items():
#         assert result.bonds[k].conjugated == v


#    "OCCn2c(=N)n(CCOc1ccc(Cl)cc1Cl)c3ccccc23", Too many bonds to n, c, n
    # "OCCN2C(=N)N(CCOC1=CC=C(Cl)C=C1Cl)C3=CC=CC=C23",