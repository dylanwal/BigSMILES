import pytest

from bigsmiles.errors import TokenizeError
from bigsmiles.tokenizer import TokenKind, tokenize

token_tests = [
    # Atoms
    ["Zn", [TokenKind.Atom]],
    ["C", [TokenKind.Atom]],
    ["CCC", [TokenKind.Atom, TokenKind.Atom, TokenKind.Atom]],
    ["c", [TokenKind.Aromatic]],
    ["[C@]", [TokenKind.AtomExtend]],
    ["[C@@]", [TokenKind.AtomExtend]],
    ["[Fe+]", [TokenKind.AtomExtend]],
    ["[Fe++]", [TokenKind.AtomExtend]],
    ["[Fe+2]", [TokenKind.AtomExtend]],
    ["[OH3+]", [TokenKind.AtomExtend]],

    # Rings
    ["C1C", [TokenKind.Atom, TokenKind.Ring, TokenKind.Atom]],
    ["C7C", [TokenKind.Atom, TokenKind.Ring, TokenKind.Atom]],
    ["C77C", [TokenKind.Atom, TokenKind.Ring, TokenKind.Ring, TokenKind.Atom]],
    ["C%77C", [TokenKind.Atom, TokenKind.Ring2, TokenKind.Atom]],

    # Bonds
    ["C=C", [TokenKind.Atom, TokenKind.Bond, TokenKind.Atom]],
    ["C#N", [TokenKind.Atom, TokenKind.Bond, TokenKind.Atom]],
    ["C/C=C/C", [TokenKind.Atom, TokenKind.BondEZ, TokenKind.Atom, TokenKind.Bond, TokenKind.Atom,
                 TokenKind.BondEZ, TokenKind.Atom]],

    # Branch
    ["(", [TokenKind.BranchStart]],
    [")", [TokenKind.BranchEnd]],
    ["C(C)C", [TokenKind.Atom, TokenKind.BranchStart, TokenKind.Atom, TokenKind.BranchEnd, TokenKind.Atom]],

    # Mixture
    [".", [TokenKind.Mix]],
    ["C.C", [TokenKind.Atom, TokenKind.Mix, TokenKind.Atom]],

    # Reaction
    ["C>CC>>CCC", [TokenKind.Atom, TokenKind.Rxn, TokenKind.Atom, TokenKind.Atom, TokenKind.Rxn, TokenKind.Atom,
                   TokenKind.Atom, TokenKind.Atom]],

    # Bonding Descriptor
    ["[>]", [TokenKind.BondDescriptor]],
    ["[<]", [TokenKind.BondDescriptor]],
    ["[$]", [TokenKind.BondDescriptor]],
    ["[>1]", [TokenKind.BondDescriptor]],
    ["[>12]", [TokenKind.BondDescriptor]],
    ["[$]CCC(=[$])C[$]", [TokenKind.BondDescriptor, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom,
                      TokenKind.BranchStart, TokenKind.Bond, TokenKind.BondDescriptor, TokenKind.BranchEnd,
                      TokenKind.Atom, TokenKind.BondDescriptor]],
    ["[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]",
     [TokenKind.AtomExtend, TokenKind.Atom, TokenKind.StochasticStart, TokenKind.BondDescriptor,
      TokenKind.BondDescriptor, TokenKind.Atom, TokenKind.BranchStart, TokenKind.Bond, TokenKind.Atom,
      TokenKind.BranchEnd, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom,
      TokenKind.BranchStart, TokenKind.Bond, TokenKind.Atom, TokenKind.BranchEnd, TokenKind.BondDescriptor,
      TokenKind.StochasticSeperator, TokenKind.BondDescriptor, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom,
      TokenKind.Atom, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom, TokenKind.BondDescriptor,
      TokenKind.BondDescriptor, TokenKind.StochasticEnd, TokenKind.AtomExtend]],

    # ladder
    ["[$1[$2]1]", [TokenKind.BondDescriptorLadder]],
    ["{[]CCC([$1[$1]1])CC(OC[$1[$2]1])(CC[$1[$1]2])OC[$1[$2]2][]}",
     [TokenKind.StochasticStart, TokenKind.ImplictEndGroup, TokenKind.Atom, TokenKind.Atom, TokenKind.Atom,
      TokenKind.BranchStart, TokenKind.BondDescriptorLadder, TokenKind.BranchEnd, TokenKind.Atom, TokenKind.Atom,
      TokenKind.BranchStart, TokenKind.Atom, TokenKind.Atom, TokenKind.BondDescriptorLadder, TokenKind.BranchEnd,
      TokenKind.BranchStart, TokenKind.Atom, TokenKind.Atom, TokenKind.BondDescriptorLadder, TokenKind.BranchEnd,
      TokenKind.Atom, TokenKind.Atom, TokenKind.BondDescriptorLadder, TokenKind.ImplictEndGroup,
      TokenKind.StochasticEnd]],

    # stochastic symbols
    ["C,C", [TokenKind.Atom, TokenKind.StochasticSeperator, TokenKind.Atom]],
    ["C;C", [TokenKind.Atom, TokenKind.StochasticSeperator, TokenKind.Atom]],
    ["C{C", [TokenKind.Atom, TokenKind.StochasticStart, TokenKind.Atom]],
    ["C}C", [TokenKind.Atom, TokenKind.StochasticEnd, TokenKind.Atom]],

]


@pytest.mark.parametrize("test_case", token_tests)
def test_tokenizer(test_case: list):
    results = tokenize(test_case[0])
    assert all([result.kind == answer for result, answer in zip(results, test_case[1])])


negative_token_tests = [
    "E",  # not an element
    "%",  # not valid on its own
    "[C<]",  # not valid
]


@pytest.mark.parametrize("test_case", negative_token_tests)
def test_tokenizer_error(test_case: str):
    with pytest.raises(TokenizeError) as p:
        tokenize(test_case)
