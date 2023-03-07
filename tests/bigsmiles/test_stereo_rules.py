
import pytest

import bigsmiles


cases = [
    ["F/C=C/F", {5: "E"}],  # E-difluoroethene
    ["F\C=C\F", {5: "Z"}],  # Z-difluoroethene
    "C(\F)=C/F",  # trans
    "F\C=C/F",  # cis
    "F/C=C\F",  # cis
    "F/C(CC)=C/F",
    "F/C=C=C=C/F",  # trans
    "F/C=C=C=C\F", # cis
    "F/C=C/C/C=C\C",
    "F/C=C/CC=CC",  # partially specified
    "CC(F)/C=C/F",
    "CCC(\F)=C/F",  #
    "CCC(\F)=C(/F)CC",
    "CC/C(\F)=C/F",  # error
]


@pytest.mark.parametrize("case", cases)
def test_whole_system_molecule(case: str):
    input_, answer = case
    result = bigsmiles.Reaction(input_)
    for chem, ans in zip(result.reactants, answer[0]):
        assert chem == ans