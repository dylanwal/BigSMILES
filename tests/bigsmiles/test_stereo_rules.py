
import pytest

import bigsmiles


cases = [
    ["F/C=C/F", {2: "E"}],  # E-difluoroethene
    ["F\C=C\F", {2: "E"}],  # E-difluoroethene
    ["F\C=C/F", {2: "Z"}],  # Z-difluoroethene
    ["F/C=C\F", {2: "Z"}],  # Z-difluoroethene
    ["C(\F)=C/F", {2: "E"}],  #
    ["F/C(CC)=C/F", {4: "E"}],

    ["F/C=C=C=C/F", {2: "E"}],
    ["F/C=C=C=C\F", {2: "Z"}],
    ["F/C=C/C/C=C\C", {2: "Z"}],
    ["F/C=C/CC=CC", {2: "Z"}],  # partially specified
    ["CC(F)/C=C/F", {2: "Z"}],
    ["CCC(\F)=C/F", {2: "Z"}],  #
    ["CCC(\F)=C(/F)CC", {2: "Z"}],
]


@pytest.mark.parametrize("case", cases)
def test_ez_cases(case: str):
    input_, answer = case
    result = bigsmiles.BigSMILES(input_)
    for k in answer:
        assert result.bonds[k].double_bond_ez == answer[k]


cases_errors = [
"CC/C(\F)=C/F",
]