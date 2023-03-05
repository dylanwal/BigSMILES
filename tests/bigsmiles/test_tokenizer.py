import pytest

from bigsmiles.errors import TokenizeError
from bigsmiles.constructors.tokenizer import TokenKind, tokenize, Token, tokenize_bonding_descriptor, \
    tokenize_atom_symbol

cases_for_atom_pattern = [
    ["C", {"element": "C", "isotope": None, "stereo": None, "hydrogens": None, "charge": 0, "class_": None}],
    ["Zn", {"element": "Zn", "isotope": None, "stereo": None, "hydrogens": None, "charge": 0, "class_": None}],
    ["c", {"element": "c", "isotope": None, "stereo": None, "hydrogens": None, "charge": 0, "class_": None}],

    # isotope
    ["[13C]", {"element": "C", "isotope": 13, "stereo": None, "hydrogens": 0, "charge": 0, "class_": None}],
    ["[238U]", {"element": "U", "isotope": 238, "stereo": None, "hydrogens": 0, "charge": 0, "class_": None}],

    # stereo
    ["[C@]", {"element": "C", "isotope": None, "stereo": "@", "hydrogens": 0, "charge": 0, "class_": None}],
    ["[C@@]", {"element": "C", "isotope": None, "stereo": "@@", "hydrogens": 0, "charge": 0, "class_": None}],

    # hydrogens
    ["[H]", {"element": "H", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 0, "class_": None}],
    ["[C]", {"element": "C", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 0, "class_": None}],
    ["[CH]", {"element": "C", "isotope": None, "stereo": None, "hydrogens": 1, "charge": 0, "class_": None}],
    ["[CH2]", {"element": "C", "isotope": None, "stereo": None, "hydrogens": 2, "charge": 0, "class_": None}],
    ["[CH3]", {"element": "C", "isotope": None, "stereo": None, "hydrogens": 3, "charge": 0, "class_": None}],
    ["[ZnH3]", {"element": "Zn", "isotope": None, "stereo": None, "hydrogens": 3, "charge": 0, "class_": None}],

    # charge
    ["[Fe+]", {"element": "Fe", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 1, "class_": None}],
    ["[Fe-]", {"element": "Fe", "isotope": None, "stereo": None, "hydrogens": 0, "charge": -1, "class_": None}],
    ["[Fe+2]", {"element": "Fe", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 2, "class_": None}],
    ["[Fe+++]", {"element": "Fe", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 3, "class_": None}],

    # atom class_
    ["[CH3:1]", {"element": "C", "isotope": None, "stereo": None, "hydrogens": 3, "charge": 0, "class_": 1}],
    ["[Zn:53]", {"element": "Zn", "isotope": None, "stereo": None, "hydrogens": 0, "charge": 0, "class_": 53}],

    # aromatic_elements
    ["[13cH:1]", {"element": "c", "isotope": 13, "stereo": None, "hydrogens": 1, "charge": 0, "class_": 1}],

]


@pytest.mark.parametrize("case", cases_for_atom_pattern)
def test_atom_tokenizer(case: list):
    input_, answer = case
    results = tokenize_atom_symbol(input_)
    for key in answer:
        assert results[key] == answer[key]


cases_for_bonding_descriptor = [
    [">", [">", 1]],
    ["<", ["<", 1]],
    ["$", ["$", 1]],
    [">2", [">", 2]],
    [">22", [">", 22]],

    ["[>]", [">", 1]],
    ["[<]", ["<", 1]],
    ["[$]", ["$", 1]],
    ["[>2]", [">", 2]],
    ["[>22]", [">", 22]],
]


@pytest.mark.parametrize("case", cases_for_bonding_descriptor)
def test_bonding_descriptor_tokenizer(case: list[str, list]):
    input_, answer = case
    results = tokenize_bonding_descriptor(input_)
    for result, ans in zip(results, answer):
        assert result == ans


cases_for_full_tokenizer = [
    # Atoms
    ["Zn", [Token(TokenKind.Atom, "Zn")]],
    ["C", [Token(TokenKind.Atom, "C")]],
    ["CCC", [Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C")]],
    ["c", [Token(TokenKind.Aromatic, "c")]],
    ["[C@]", [Token(TokenKind.AtomExtend, "[C@]")]],
    ["[C@@]", [Token(TokenKind.AtomExtend, "[C@@]")]],
    ["[Fe+]", [Token(TokenKind.AtomExtend, "[Fe+]")]],
    ["[Fe++]", [Token(TokenKind.AtomExtend, "[Fe++]")]],
    ["[Fe+2]", [Token(TokenKind.AtomExtend, "[Fe+2]")]],
    ["[OH3+]", [Token(TokenKind.AtomExtend, "[OH3+]")]],

    # Rings
    ["C1C", [Token(TokenKind.Atom, "C"), Token(TokenKind.Ring, "1"), Token(TokenKind.Atom, "C")]],
    ["C7C", [Token(TokenKind.Atom, "C"), Token(TokenKind.Ring, "7"), Token(TokenKind.Atom, "C")]],
    ["C77C",
     [Token(TokenKind.Atom, "C"), Token(TokenKind.Ring, "7"), Token(TokenKind.Ring, "7"), Token(TokenKind.Atom, "C")]],
    ["C%77C", [Token(TokenKind.Atom, "C"), Token(TokenKind.Ring2, "%77"), Token(TokenKind.Atom, "C")]],

    # Bonds
    ["C=C", [Token(TokenKind.Atom, "C"), Token(TokenKind.Bond, "="), Token(TokenKind.Atom, "C")]],
    ["C#N", [Token(TokenKind.Atom, "C"), Token(TokenKind.Bond, "#"), Token(TokenKind.Atom, "N")]],
    ["C/C=C/C", [Token(TokenKind.Atom, "C"), Token(TokenKind.BondEZ, "/"), Token(TokenKind.Atom, "C"),
                 Token(TokenKind.Bond, "="), Token(TokenKind.Atom, "C"), Token(TokenKind.BondEZ, "/"),
                 Token(TokenKind.Atom, "C")]],

    # Branch
    ["(", [Token(TokenKind.BranchStart, "(")]],
    [")", [Token(TokenKind.BranchEnd, ")")]],
    ["C(C)C", [Token(TokenKind.Atom, "C"), Token(TokenKind.BranchStart, "("), Token(TokenKind.Atom, "C"),
               Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Atom, "C")]],

    # Mixture
    [".", [Token(TokenKind.Disconnected, ".")]],
    ["C.C", [Token(TokenKind.Atom, "C"), Token(TokenKind.Disconnected, "."), Token(TokenKind.Atom, "C")]],

    # Reaction
    ["C>CC>CCC", [Token(TokenKind.Atom, "C"), Token(TokenKind.Rxn, ">"), Token(TokenKind.Atom, "C"),
                  Token(TokenKind.Atom, "C"), Token(TokenKind.Rxn, ">"), Token(TokenKind.Atom, "C"),
                  Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C")]],

    # Bonding Descriptor
    ["[>]", [Token(TokenKind.BondDescriptor, "[>]")]],
    ["[<]", [Token(TokenKind.BondDescriptor, "[<]")]],
    ["[$]", [Token(TokenKind.BondDescriptor, "[$]")]],
    ["[>1]", [Token(TokenKind.BondDescriptor, "[>1]")]],
    ["[>12]", [Token(TokenKind.BondDescriptor, "[>12]")]],
    ["[$]CCC(=[$])C[$]",
     [Token(TokenKind.BondDescriptor, "[$]"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.BranchStart, "("), Token(TokenKind.Bond, "="),
      Token(TokenKind.BondDescriptor, "[$]"), Token(TokenKind.BranchEnd, ")"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.BondDescriptor, "[$]")]],
    ["[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]",
     [Token(TokenKind.AtomExtend, "[H]"), Token(TokenKind.Atom, "O"), Token(TokenKind.StochasticStart, "{"),
      Token(TokenKind.BondDescriptor, "[>]"), Token(TokenKind.BondDescriptor, "[<]"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.BranchStart, "("), Token(TokenKind.Bond, "="), Token(TokenKind.Atom, "O"),
      Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.BranchStart, "("), Token(TokenKind.Bond, "="), Token(TokenKind.Atom, "O"),
      Token(TokenKind.BranchEnd, ")"), Token(TokenKind.BondDescriptor, "[<]"),
      Token(TokenKind.StochasticSeperator, ","),
      Token(TokenKind.BondDescriptor, "[>]"), Token(TokenKind.Atom, "N"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "N"), Token(TokenKind.BondDescriptor, "[>]"),
      Token(TokenKind.BondDescriptor, "[<]"), Token(TokenKind.StochasticEnd, "}"), Token(TokenKind.AtomExtend, "[H]")]],

    # ladder
    ["[$1[$2]1]", [Token(TokenKind.BondDescriptorLadder, "[$1[$2]1]")]],
    ["{[]CCC([$1[$1]1])CC(OC[$1[$2]1])(CC[$1[$1]2])OC[$1[$2]2][]}",
     [Token(TokenKind.StochasticStart, "{"), Token(TokenKind.ImplictEndGroup, "[]"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"), Token(TokenKind.BranchStart, "("),
      Token(TokenKind.BondDescriptorLadder, "[$1[$1]1]"), Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.BranchStart, "("), Token(TokenKind.Atom, "O"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.BondDescriptorLadder, "[$1[$2]1]"), Token(TokenKind.BranchEnd, ")"),
      Token(TokenKind.BranchStart, "("), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.BondDescriptorLadder, "[$1[$1]2]"), Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Atom, "O"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.BondDescriptorLadder, "[$1[$2]2]"),
      Token(TokenKind.ImplictEndGroup, "[]"), Token(TokenKind.StochasticEnd, "}")]],

    # stochastic symbols
    ["C,C", [Token(TokenKind.Atom, "C"), Token(TokenKind.StochasticSeperator, ","), Token(TokenKind.Atom, "C")]],
    ["C;C", [Token(TokenKind.Atom, "C"), Token(TokenKind.StochasticSeperator, ";"), Token(TokenKind.Atom, "C")]],
    ["C{C", [Token(TokenKind.Atom, "C"), Token(TokenKind.StochasticStart, "{"), Token(TokenKind.Atom, "C")]],
    ["C}C", [Token(TokenKind.Atom, "C"), Token(TokenKind.StochasticEnd, "}"), Token(TokenKind.Atom, "C")]],

    # Atom class
    ["[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH3:6]",
     [Token(TokenKind.AtomExtend, "[CH3:1]"), Token(TokenKind.AtomExtend, "[CH2:2]"),
      Token(TokenKind.AtomExtend, "[CH2:3]"), Token(TokenKind.AtomExtend, "[CH2:4]"),
      Token(TokenKind.AtomExtend, "[CH2:5]"), Token(TokenKind.AtomExtend, "[CH3:6]")]],
    ["[CH3:1]CC[CH2:2][CH2:3][CH3:4]",
     [Token(TokenKind.AtomExtend, "[CH3:1]"), Token(TokenKind.Atom, "C"),
      Token(TokenKind.Atom, "C"), Token(TokenKind.AtomExtend, "[CH2:2]"),
      Token(TokenKind.AtomExtend, "[CH2:3]"), Token(TokenKind.AtomExtend, "[CH3:4]")]],

    # aromatic_elements
    ["c", [Token(TokenKind.Aromatic, "c")]],
    ["c1ccccc1", [Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1"), Token(TokenKind.Aromatic, "c"),
                  Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"),
                  Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1")]],
    ["c1[13cH]cc[c+]c1", [Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1"),
                          Token(TokenKind.AtomExtend, "[13cH]"), Token(TokenKind.Aromatic, "c"),
                          Token(TokenKind.Aromatic, "c"), Token(TokenKind.AtomExtend, "[c+]"),
                          Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1")]],
    ["c1ccccc1-c2ccccc2", [
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "1"), Token(TokenKind.Bond, "-"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "2"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Ring, "2")
    ]],
    ["OCCn2c(=N)nc3ccccc23", [
        Token(TokenKind.Atom, "O"), Token(TokenKind.Atom, "C"), Token(TokenKind.Atom, "C"),
        Token(TokenKind.Aromatic, "n"), Token(TokenKind.Ring, "2"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.BranchStart, "("), Token(TokenKind.Bond, "="), Token(TokenKind.Atom, "N"),
        Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Aromatic, "n"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Ring, "3"), Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"), Token(TokenKind.Aromatic, "c"),
        Token(TokenKind.Ring, "2"), Token(TokenKind.Ring, "3")

    ]],

]


@pytest.mark.parametrize("case", cases_for_full_tokenizer)
def test_tokenizer(case: list):
    input_, answer = case
    results = tokenize(input_)
    for result, ans in zip(results, answer):
        assert result == ans


negative_token_tests = [
    "E",  # not an element
    "%",  # not valid on its own
    "[C<]",  # not valid
    "[CH2:]",
    "[:3]",
    "[12]"
]


@pytest.mark.parametrize("test_case", negative_token_tests)
def test_tokenizer_error(test_case: str):
    with pytest.raises(TokenizeError) as p:
        tokenize(test_case)
