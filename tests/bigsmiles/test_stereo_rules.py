
import pytest

import bigsmiles


cases = [
    ["F/C=C/F", {1: "E"}],  # E-difluoroethene
    ["F\C=C\F", {1: "E"}],  # E-difluoroethene
    ["F\C=C/F", {1: "Z"}],  # Z-difluoroethene
    ["F/C=C\F", {1: "Z"}],  # Z-difluoroethene
    ["C(\F)=C/F", {1: "E"}],  #
    ["F/C(CC)=C/F", {3: "E"}],
    ["CC(F)/C=C/F", {3: "E"}],
    ["CCC(\F)=C/F", {3: "E"}],
    ["CCC(\F)=C(/F)CC", {3: "E"}],
    ["CCC(\F)=C(\F)CC", {3: "Z"}],
    ["CCCCC(\CCCCF)=C(\CCCCF)CCCC", {9: "Z"}],
    ["CCCC(C)C(\C(C)CCCF)=C(\CCCCF)CCCC", {11: "Z"}],
    ["C(CCCC)C(\C(CCCF)C)=C(\CCCCF)CCCC", {11: "Z"}],
    ["FC1C2C3C/C(CCC3C2C1)=C4CC5CCC(F)C5C/4", {9: "Z"}],  # rings and bridges

    # conjugated
    ["F/C=C/C/C=C\C", {1: "E", 4: "Z"}],
    ["F/C=C/CC=CC", {1: "E", 4: None}],  # partially specified

    # Cumulene
    ["F/C=C=C=C/F", {1: None, 2: "E", 3: None}],
    ["F/C=C=C=C\F", {1: None, 2: "Z", 3: None}],
    ["F/C(CC)=C=C=C\F", {3: None, 4: "Z", 5: None}],
    ["F/C(CC)=C=C=C(CC)\F", {3: None, 4: "Z", 5: None}],
    ["F/C(CC)=C=C=C(\CC)F", {3: None, 4: "E", 5: None}],

]


@pytest.mark.parametrize("case", cases)
def test_ez_cases(case: str):
    input_, answer = case
    result = bigsmiles.BigSMILES(input_)
    for k in answer:
        assert result.bonds[k].double_bond_ez == answer[k]


cases_errors = [
    # missing '/'
    "CC/C=CCC",
    "CC/C=C",
    "CC/C(\F)=C/F",

    # terminal H
    "CC/C=C/",
]

