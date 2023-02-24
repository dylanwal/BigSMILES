"""

Code for tokenizing a BigSMILES string.

"""

import enum
import re

from bigsmiles.errors import TokenizeError
from bigsmiles.config import Config


class TokenKind(enum.Enum):
    Bond = 0
    Atom = 1
    Aromatic = 2
    AtomExtend = 3
    BranchStart = 4
    BranchEnd = 5
    Ring = 6
    Ring2 = 7
    BondEZ = 8
    Mix = 9
    Rxn = 10
    BondDescriptor = 11
    StochasticSeperator = 12
    StochasticStart = 13
    StochasticEnd = 14
    ImplictEndGroup = 15
    BondDescriptorLadder = 16


_isotope_pattern = r'(?P<isotope>[\d]{1,3})?'
_element_pattern = r'(?P<element>' + "|".join(Config.elements_ordered) + "|".join(Config.aromatic) + '{1})'
_stereo_pattern = r'(?P<stereo>@{1,2})?'
_hydrogen_pattern = r'(?P<hcount>H[\d]?)?'
_charge_pattern = r'(?P<charge>[-|\+]{1,3}[\d]?)?'
atom_pattern = r"(?:\[)" + _isotope_pattern + _element_pattern + _stereo_pattern + \
               _hydrogen_pattern + _charge_pattern + r"(?:\])"

token_specification = [
    # order in the list is important; regex stops at first match
    (TokenKind.Bond.name, r'[=|#]'),
    (TokenKind.Atom.name, "|".join(Config.elements_ordered)),
    (TokenKind.Aromatic.name, "|".join(Config.aromatic)),
    (TokenKind.AtomExtend.name, atom_pattern),  # Atom in brackets
    (TokenKind.BranchStart.name, r'\('),
    (TokenKind.BranchEnd.name, r'\)'),
    (TokenKind.Ring.name, r'[\d]{1}'),
    (TokenKind.Ring2.name, r'%[\d]{2}'),  # ring_id with two-digit numbers
    (TokenKind.BondEZ.name, r'/|\\'),  # cis trans
    (TokenKind.Mix.name, r"\."),  # mixture
    (TokenKind.Rxn.name, r">>|>"),  # reaction -->

    (TokenKind.BondDescriptorLadder.name, r"\[[$<>][\d]\[[$<>][\d]?\][\d]?\]"),  # Ladder
    (TokenKind.BondDescriptor.name, r"\[[$<>][\d]?[\d]?\]"),
    (TokenKind.StochasticSeperator.name, r",|;"),
    (TokenKind.StochasticStart.name, r'\{'),
    (TokenKind.StochasticEnd.name, r'\}'),
    (TokenKind.ImplictEndGroup.name, r'\[\]'),

    ('SKIP', r'[ \t]+'),  # Skip over spaces and tabs
    ('MISMATCH', r'.'),  # Any other character
]

tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)

tok_regex = re.compile(tok_regex)


class Token:
    __slots__ = ("kind", "value")

    def __init__(self, kind: TokenKind, value: str):
        self.kind = kind
        self.value = value

    def __str__(self):
        return f"{self.kind}: {self.value}"


def tokenize(text: str) -> list[Token]:
    """

    Tokenizes a bigSMILES string.

    Parameters
    ----------
    text: str
        BigSMILES string

    Returns
    -------
    result: list[Token]
        BigSMILES as a token list

    """
    result = []
    for match in re.finditer(tok_regex, text.strip()):
        kind = match.lastgroup

        value = match.group()
        if kind == 'SKIP':
            continue
        elif kind == 'MISMATCH':
            raise TokenizeError(f'Invalid symbol (or group of symbols). (starting with {value!r}; '
                                f'index: {match.span()[0]})'
                                f'\n{text}' + "\n" + " " * match.span()[0] + "^(and forward)")

        result.append(
            Token(TokenKind[kind], value)
        )

    return result
