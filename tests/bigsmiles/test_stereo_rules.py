
import pytest

import bigsmiles


cases = [
    ["F/C=C/F", {5: "E"}],  # E-difluoroethene
    ["F\C=C\F", {5: "Z"}],  # Z-difluoroethene
    ["C(\F)=C/F", {5: "Z"}],  # trans
    ["F\C=C/F", {5: "Z"}],  # cis
    ["F/C=C\F", {5: "Z"}],  # cis
    ["F/C(CC)=C/F", {5: "Z"}],
    ["F/C=C=C=C/F", {5: "Z"}],  # trans
    ["F/C=C=C=C\F", {5: "Z"}], # cis
    ["F/C=C/C/C=C\C", {5: "Z"}],
    ["F/C=C/CC=CC", {5: "Z"}],  # partially specified
    ["CC(F)/C=C/F", {5: "Z"}],
    ["CCC(\F)=C/F", {5: "Z"}],  #
    ["CCC(\F)=C(/F)CC", {5: "Z"}],
]


@pytest.mark.parametrize("case", cases)
def test_whole_system_molecule(case: str):
    input_, answer = case
    result = bigsmiles.Reaction(input_)
    for chem, ans in zip(result.reactants, answer[0]):
        assert chem == ans


cases_errors = [
"CC/C(\F)=C/F",
]