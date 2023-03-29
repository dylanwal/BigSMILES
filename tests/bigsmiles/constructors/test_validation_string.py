import pytest

import bigsmiles.errors as errors
import bigsmiles.constructors.validation.validation_string as validation

ring_renumbering_cases = [
    ["O1CCCCC1N1CCCCC1", "O1CCCCC1N2CCCCC2"],
    ["O1CC1CC1C1N1CCCCC1", "O1CC1CC2C2N3CCCCC3"],
    ["[12C]C1=C(C(=NC=C1C=O)C)O", "[12C]C1=C(C(=NC=C1C=O)C)O"],  # no change
    ["C1C(C)([CH7])[N+7]([7C-4])C1(=O)C14C3C2C1C5C2C3C45C(=O)C69C8[$1]C7C6C%10C7C8C9%10[>7]",
     "C1C(C)([CH7])[N+7]([7C-4])C1(=O)C%114C3C2C%11C5C2C3C45C(=O)C69C8[$1]C7C6C%10C7C8C9%10[>7]"]
]


@pytest.mark.parametrize("case", ring_renumbering_cases)
def test_ring_renumbering(case):
    input_, answer = case
    result = validation.validate_ring_numbering(input_, True)
    assert result == answer


ring_renumbering_cases_error = [
    ["O1CCCCC1N1CCCCC1", "O1CCCCC1N2CCCCC2"],
    ["O1CC1CC1C1N1CCCCC1", "O1CC1CC2C2N3CCCCC3"],
    ["C1C(C)([CH7])[N+7]([7C-4])C1(=O)C14C3C2C1C5C2C3C45C(=O)C69C8[$1]C7C6C%10C7C8C9%10[>7]",
     "C1C(C)([CH7])[N+7]([7C-4])C1(=O)C%114C3C2C%11C5C2C3C45C(=O)C69C8[$1]C7C6C%10C7C8C9%10[>7]"]
]


@pytest.mark.parametrize("case", ring_renumbering_cases_error)
def test_ring_renumbering_error(case):
    input_, answer = case
    with pytest.raises(errors.BigSMILESError):
        validation.validate_ring_numbering(input_, False)
