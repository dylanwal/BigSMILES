import re
from bigsmiles.config import Config


ATOM_PATTERN = re.compile(
    r"(?P<isotope>\d{1,3})?" +
    r'(?P<element>' + "|".join(Config.elements_ordered) + "|".join(Config.aromatic) + ')' +
    r"(?P<stereo>@{0,2})?" +
    r'(?P<hydrogens>H\d?)?' +
    r'(?P<charge>[-|+]{1,3}\d?)?'
)


def atom_symbol_to_attributes(symbol: str) -> dict:
    return re.match(ATOM_PATTERN, symbol).groupdict()


def process_bonding_descriptor_symbol(symbol: str) -> tuple[str, int]:
    """
    Process bonding descriptor into symbol and index

    Parameters
    ----------
    symbol: str
        bonding descriptor symbol
        '[>1]'

    Returns
    -------
    symbol: str
        bonding descriptor symbol
        '<', '>', '$', ''
    index: int
        bonding descriptor index

    """
    symbol = symbol.replace("[", "").replace("]", "")
    if not symbol:
        return symbol, 0

    if symbol[-1].isdigit():
        return symbol[0], int(symbol[-1])

    return symbol, 0

