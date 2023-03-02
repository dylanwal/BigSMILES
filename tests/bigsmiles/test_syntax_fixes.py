import pytest

import bigsmiles

bigsmiles.Config.show_aromatic_bond = False

branch_cases = [
    # Branch syntax fixes
    ['CC(CC)(CC)', 'CC(CC)CC'],
    ['CC(CC(CC))', 'CCCCCC'],
    ['CCCCCC()', 'CCCCCC'],
    ["C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C", "C(CCCCCCCCCCCCCCCCCCCC)C"],
    ["{[][$]CC(c1cc(C(=O)Oc2ccc(OCCCC)cc2)ccc1(C(=O)Oc3ccc(OCCCC)cc3))[$][$]}{[>][<]Si(C)(C)O[>][]}",
     "{[][$]CC(c1cc(C(=O)Oc2ccc(OCCCC)cc2)ccc1C(=O)Oc3ccc(OCCCC)cc3)[$][$]}{[>][<]Si(C)(C)O[>][]}"],
]


@pytest.mark.parametrize("case", branch_cases)
def test_branch_syntax_fixes(case):
    """
    mainly testing the following functions:
    * bigsmiles.syntax_fixes.remove_unnecessary_branch_symbols() --> remove unnecessary branches
    * bigsmiles.constructor.BigSMILESConstructor.close_branch() --> remove empty branches
    """
    result = bigsmiles.BigSMILES(case[0])
    assert str(result) == case[1]


ring_cases = [
    # Branch syntax fixes
    ['C12CCCC12C', 'C=1CCCC1C'],
    ['C2CCCC2C', 'C1CCCC1C'],
    ['CC(CC2CCC2)CC', 'CC(CC1CCC1)CC'],
    ['C%10CCCC%10C', 'C1CCCC1C'],
    ['C1C2C3C4C5C6C7C8C9C%10C%11OC1C2C3C4C5C6C7C8C9C%10C%11', 'C1C2C3C4C5C6C7C8C9C%10C%11OC1C2C3C4C5C6C7C8C9C%10C%11'],
]


@pytest.mark.parametrize("case", ring_cases)
def test_ring_syntax_fixes(case):
    """
    mainly testing the following functions:
    * bigsmiles.constructor.BigSMILESConstructor.add_ring() --> two rings goes to a double bond
    * bigsmiles.constructor.BigSMILESConstructor.add_ring_from_atoms() --> two rings goes to a double bond
    """
    result = bigsmiles.BigSMILES(case[0])
    assert str(result) == case[1]


