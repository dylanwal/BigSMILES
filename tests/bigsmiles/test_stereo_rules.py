
import pytest

import bigsmiles


cases = [
    ["F/C=C/F", {1: "E"}],  # E-difluoroethene
    [r"F\C=C\F", {1: "E"}],  # E-difluoroethene
    [r"F\C=C/F", {1: "Z"}],  # Z-difluoroethene
    [r"F/C=C\F", {1: "Z"}],  # Z-difluoroethene
    [r"C(\F)=C/F", {1: "E"}],  #
    ["F/C(CC)=C/F", {3: "E"}],
    ["CC(F)/C=C/F", {3: "E"}],
    [r"CCC(\F)=C/F", {3: "E"}],
    [r"CCC(\F)=C(/F)CC", {3: "E"}],
    [r"CCC(\F)=C(\F)CC", {3: "Z"}],
    [r"CCCCC(\CCCCF)=C(\CCCCF)CCCC", {9: "Z"}],
    [r"CCCC(C)C(\C(C)CCCF)=C(\CCCCF)CCCC", {11: "Z"}],
    [r"C(CCCC)C(\C(CCCF)C)=C(\CCCCF)CCCC", {11: "Z"}],
    ["C/C(CC)=C(C(OCC)C(CF)C)/C(C(CC)C)OCC", {3: "E"}],
    ["C/C(CC)=C(C(CCC)C(CF)C)/C(C(CC)C)=C/CC", {3: "Z", 17: "E"}],
    ["F/C(C)=C(C(C(CC)F)CC(CCl)CBr)/C(C(Cl)CC)CC(CF)CBr", {2: "E"}],
    ["F/C(C)=C(C(C(CC)F)OC(CCl)NBr)/C(C(CC)Cl)OC(CF)NBr", {2: "E"}],
    ["F/C(C)=C(C(C(F)=C)OC(CCl)NBr)/C(CCl)OC(CF)NBr", {2: "E"}],
    ["F/C(C)=C(C(C(F)=C)OC(CCl)NBr)/C(C(C)(C)F)OC(CF)NBr", {2: "Z"}],
    [r"F/C=C(CC(O)(O)OCC)\CC(OCC)=O", {2: "Z"}],
    [r"F/C=C(CC(O)(OC)OCC)\CC(OCC)=O", {2: "E"}],  # chemdraw E, rdkit Z
    ["F/C=C(C/C(OCC)=N/C)/CC(NC)(OCC)NCC", {2, "Z"}],
    ["F/C=C(C/C(OC)=N/CC)/CC(NC)(OCC)N(C)C", {2: "Z"}],  # chemdraw E, rdkit Z
    [r"F/C(C)=C(C(OC(NBr)CF)C(C)(F)C)\C(OC(C)C)C(F)=C", {2: "E"}],
    [r"F/C=C(CC(OC)(N(C)CC)N)\C/C(OC)=N/CC", {2: "E", 22: "Z"}],  # chemdraw E, rdkit Z (rdkit wrong)

    # rings and bridges
    ["FC1C2C3C/C(CCC3C2C1)=C4CC5CCC(F)C5C/4", {9: "Z"}],

    # conjugated
    [r"F/C=C/C/C=C\C", {1: "E", 4: "Z"}],
    ["F/C=C/CC=CC", {1: "E", 4: None}],  # partially specified

    # Cumulene
    ["F/C=C=C=C/F", {1: None, 2: "E", 3: None}],
    [r"F/C=C=C=C\F", {1: None, 2: "Z", 3: None}],
    [r"F/C(CC)=C=C=C\F", {3: None, 4: "Z", 5: None}],
    [r"F/C(CC)=C=C=C(CC)\F", {3: None, 4: "Z", 5: None}],
    [r"F/C(CC)=C=C=C(\CC)F", {3: None, 4: "E", 5: None}],

    # double bonds
    [r"C(=C\F)(/C=C\C)\C=C\C", {1: "E", 5: "Z", 9: "E"}],
    [r"C/1\2=C\F.C\1=C\C.C/2=C\C", {1: "E", 5: "Z", 9: "E"}],
    [r"C/C=C/C(/C=C\C)=C(/C=C/C)\C=C/C", {1: "E", 4: "Z", 6: "E", 8: "E", 11: "Z"}],

]


# @pytest.mark.parametrize("case", cases)
# def test_ez_cases(case: str):
#     input_, answer = case
#     result = bigsmiles.BigSMILES(input_)
#     for k in answer:
#         assert result.bonds[k].double_bond_ez == answer[k]


cases_errors = [
    # missing '/'
    "CC/C=CCC",
    "CC/C=C",
    r"CC/C(\F)=C/F",
    r"F/C=C(/C=C/C)/C=C\C",  # two stereo symbols for first double bond right side
    # terminal H
    "CC/C=C/",
]

