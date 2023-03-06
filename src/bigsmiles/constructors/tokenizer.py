"""

Code for tokenizing a BigSMILES string.

"""

import enum
import re

from bigsmiles.errors import TokenizeError
import bigsmiles.reference_data.chemical_data as chemical_data


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
    Disconnected = 9
    Rxn = 10
    BondDescriptor = 11
    StochasticSeperator = 12
    StochasticStart = 13
    StochasticEnd = 14
    ImplictEndGroup = 15
    BondDescriptorLadder = 16


_isotope_pattern = r'(?P<isotope>[\d]{1,3})?'
_element_pattern = r'(?P<element>' + "|".join(chemical_data.aromatic_elements) + "|" \
                   + "|".join(chemical_data.elements_ordered) + ')'
_stereo_pattern = r'(?P<stereo>@{1,2})?'
_hydrogen_pattern = r'(?P<hydrogens>H[\d]?)?'
_charge_pattern = r'(?P<charge>[-|\+]{1,3}[\d]?)?'
_class_pattern = r'(?P<class_>:\d{1,3})?'
atom_pattern = r"(?:\[)" + _isotope_pattern + _element_pattern + _stereo_pattern + \
               _hydrogen_pattern + _charge_pattern + _class_pattern + r"(?:\])"

token_specification = [
    # order in the list is important; regex stops at first match
    (TokenKind.Bond.name, r'[-|=|#|$]'),
    (TokenKind.Atom.name, "|".join(chemical_data.organic_ordered)),
    (TokenKind.Aromatic.name, "|".join(chemical_data.aromatic_elements)),
    (TokenKind.AtomExtend.name, atom_pattern),  # Atom in brackets
    (TokenKind.BranchStart.name, r'\('),
    (TokenKind.BranchEnd.name, r'\)'),
    (TokenKind.Ring.name, r'[\d]{1}'),
    (TokenKind.Ring2.name, r'%[\d]{2}'),  # ring_id with two-digit numbers
    (TokenKind.BondEZ.name, r'/|\\'),  # cis trans
    (TokenKind.Disconnected.name, r"\."),  # mixture

    (TokenKind.BondDescriptorLadder.name, r"\[[$<>][\d]\[[$<>][\d]?\][\d]?\]"),  # Ladder
    (TokenKind.BondDescriptor.name, r"\[[$<>][\d]?[\d]?\]"),
    (TokenKind.StochasticSeperator.name, r",|;"),
    (TokenKind.StochasticStart.name, r'\{'),
    (TokenKind.StochasticEnd.name, r'\}'),
    (TokenKind.ImplictEndGroup.name, r'\[\]'),

    (TokenKind.Rxn.name, r">>|>"),  # reaction -->

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
        return f"{self.kind.name}: {self.value}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self.kind != other.kind:
            return False
        if self.value != other.value:
            return False
        return True


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
    for match in re.finditer(tok_regex, text.replace(" ", "")):
        kind = match.lastgroup

        value = match.group()
        if kind == 'SKIP':
            continue
        elif kind == 'MISMATCH':
            if text[:2] in chemical_data.element_symbols:
                raise TokenizeError(f"Invalid symbol. If the element is not in the list below it must be in []."
                                    f"\nNon-bracket elements: {chemical_data.organic_elements}"
                                    f"\n(starting with {value!r}; index: {match.span()[0]})"
                                f'\n{text}' + "\n" + " " * match.span()[0] + "^(and forward)")

            raise TokenizeError(f'Invalid symbol (or group of symbols). (starting with {value!r}; '
                                f'index: {match.span()[0]})'
                                f'\n{text}' + "\n" + " " * match.span()[0] + "^(and forward)")

        result.append(
            Token(TokenKind[kind], value)
        )

    return result


ATOM_PATTERN = re.compile(atom_pattern)


def tokenize_atom_symbol(symbol: str) -> dict:
    if symbol in chemical_data.elements_aromatic:
        return {"element": symbol, "isotope": None, "stereo": None, "hydrogens": None, "charge": 0, "class_": None}

    try:
        results = ATOM_PATTERN.match(symbol).groupdict()

        if results["isotope"] is not None:
            results["isotope"] = int(results["isotope"])

        if results["hydrogens"] is None:
            results["hydrogens"] = 0
        else:
            if results["hydrogens"][-1].isdigit():
                results["hydrogens"] = int(results["hydrogens"][-1])
            else:
                results["hydrogens"] = 1

        if results["charge"] is None:
            results["charge"] = 0
        elif results["charge"] == "+":
            results["charge"] = 1
        elif results["charge"] == "-":
            results["charge"] = -1
        elif results["charge"].count("+") > 1:
            results["charge"] = results["charge"].count("+")
        elif results["charge"].count("-") > 1:
            results["charge"] = -1 * results["charge"].count("-")
        else:
            results["charge"] = int(results["charge"].replace("+", ""))

        if results["class_"] is not None:
            results["class_"] = int(results["class_"][1:])  # index 0 is ":", so skip it

        return results
    except AttributeError:
        raise TokenizeError(f"Issue tokenizing atom: {symbol}")


DEFAULT_BONDING_DESCRIPTOR_INDEX = 1


def tokenize_bonding_descriptor(symbol: str) -> tuple[str, int]:
    symbol = symbol.replace("[", "").replace("]", "")
    if not symbol:
        return symbol, DEFAULT_BONDING_DESCRIPTOR_INDEX

    if symbol[-1].isdigit():
        return symbol[0], int(symbol[1:])

    return symbol, DEFAULT_BONDING_DESCRIPTOR_INDEX
