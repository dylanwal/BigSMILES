from __future__ import annotations
import enum
import re

from bigsmiles.config import Config


class AtomChirality(enum.Enum):
    none = ""
    R = "@"
    S = "@@"


def get_charge_number(symbol_text: str, charge_symbol: str) -> tuple[str, int]:
    count = symbol_text.count(charge_symbol)
    if count > 1:
        # format like Fe++
        return symbol_text[:-count], count

    index_positive_symbol = symbol_text.find(charge_symbol)
    if index_positive_symbol == len(symbol_text) - 1:
        # format like Fe+
        return symbol_text[:-1], count

    # format like Fe+2
    return symbol_text[:-2], int(symbol_text[-1])


def get_charge(symbol_text: str) -> tuple[str, int]:
    if "+" in symbol_text:
        return get_charge_number(symbol_text, "+")
    elif "-" in symbol_text:
        symbol_text, count = get_charge_number(symbol_text, "-")
        return symbol_text, -1 * count

    return symbol_text, 0


def get_chirality(symbol_text: str) -> tuple[str, AtomChirality]:
    result = re.search(r'@{1,2}', symbol_text)
    if result is None:
        return symbol_text, AtomChirality.none

    return symbol_text.replace("@", ""), AtomChirality(result.group())


def get_isotope(symbol_text: str) -> tuple[str, int | None]:
    if symbol_text[0].isdigit():
        isotope = re.search(r'\d+', symbol_text).group()
        return symbol_text.lstrip(isotope), int(isotope)

    return symbol_text, None


def get_hydrogens(symbol_text: str) -> tuple[str, int | None]:
    if symbol_text[0] == "H":
        return "H", None

    result = re.search(r'H\d?', symbol_text)
    if result is None:
        return symbol_text, None

    if result.group()[-1].isdigit():
        count = int(result.group()[-1])
    else:
        count = 1

    return symbol_text[:result.span()[0]] + symbol_text[result.span()[1]:], count


def process_atom_symbol(symbol_text: str) -> tuple[str, tuple[int] | None, AtomChirality, int, int | None, int | None]:
    """
    Process atom symbol into symbol, valence, chiral, charge, isotope, hydrogens

    Parameters
    ----------
    symbol_text: str
        atom symbol
        'C', '[235U]', '	[OH3+]'

    Returns
    -------
    symbol_text: str
        Atom symbol
    valance: tuple[int]
        Accepted valance electrons - From 'Config'
    chiral: AtomChirality
        chirality
    charge: int
        electron charge
    isotope: int
        isotope number
    hydrogens: int
        explict hydrogens

    """
    if "[" in symbol_text:
        symbol_text = symbol_text.replace('[', "").replace(']', '')
        symbol_text, isotope = get_isotope(symbol_text)
        symbol_text, chiral = get_chirality(symbol_text)
        symbol_text, hydrogens = get_hydrogens(symbol_text)
        symbol_text, charge = get_charge(symbol_text)
    else:
        charge = 0
        chiral = AtomChirality.none
        isotope = None
        hydrogens = None

    if symbol_text in Config.atoms_with_valence:
        valence = Config.atoms_with_valence[symbol_text]["valence"]
    else:
        valence = None

    return symbol_text, valence, chiral, charge, isotope, hydrogens


class Atom:
    __slots__ = ["id_", "valance_possible", "chiral", "charge", "isotope", "_hydrogens", "symbol", "valance",
                 "organic", "bonds"]
    _tree_print_repr = True

    def __init__(self, symbol: str, id_: int = None):
        self.id_ = id_
        self.symbol, self.valance_possible, self.chiral, self.charge, self.isotope, self._hydrogens = \
            process_atom_symbol(symbol)
        self.valance = self.valance_possible[0] if self.valance_possible is not None else None
        self.organic = True if self.symbol in Config.organics else False
        self.bonds = []

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return str(self) + "{" + str(self.id_) + "}"

    @property
    def bonds_available(self) -> int | None:
        if self.valance is None:
            return None
        return self.valance - self.number_of_bonds

    @property
    def number_of_bonds(self) -> int:
        hydrogen_count = 0
        if self._hydrogens is not None:
            hydrogen_count = self._hydrogens
        return sum((bond.type_.value for bond in self.bonds)) + hydrogen_count

    @property
    def ring_indexes(self) -> list[int]:
        ring_index = []
        for bond in self.bonds:
            if bond.ring_id is not None:
                ring_index.append(bond.ring_id)

        return ring_index

    def to_string(self, show_hydrogens: bool = False) -> str:
        """
        Construct Atom symbol

        Parameters
        ----------
        show_hydrogens: bool
            add implicit hydrogens to string

        Returns
        -------
        text: str
            atom symbol

        """
        text = self.symbol
        bracket_flag = False

        if self.symbol == "H":
            bracket_flag = True

        if self.chiral is not AtomChirality.none:
            text += self.chiral.value
            if self.bonds_available == 1:
                show_hydrogens = True
            bracket_flag = True

        if (show_hydrogens and self.bonds_available is not None and self.bonds_available > 0) or \
                self._hydrogens is not None:
            if self._hydrogens is not None:
                text += f"H{self._hydrogens if self._hydrogens > 1 else ''}"
            else:
                text += f"H{self.bonds_available if self.bonds_available > 1 else ''}"

            bracket_flag = True

        if self.charge > 0:
            text += f"+{self.charge if self.charge > 1 else ''}"
            bracket_flag = True
        elif self.charge < 0:
            text += f"{self.charge}"
            bracket_flag = True

        if self.isotope is not None:
            text = str(self.isotope) + text
            bracket_flag = True

        if bracket_flag:
            text = "[" + text + "]"

        return text + "".join((str(id_) for id_ in self.ring_indexes))


class BondType(enum.Enum):
    single = 1
    double = 2
    triple = 3


bond_mapping = {
    "": BondType.single,
    "=": BondType.double,
    "#": BondType.triple
}

#  Not implemented yet
#
# class BondConfigTypes(enum.Enum):
#     none = ""
#     E = "E"
#     Z = "Z"
#
#
# class BondConfig:
#     __slots__ = ["type_", "left_up", "left_down", "right_up", "right_down"]
#
#     def __init__(self,
#                  left_up: Atom | BondDescriptorAtom | Branch | StochasticObject,
#                  left_down: Atom | BondDescriptorAtom | Branch | StochasticObject,
#                  right_up: Atom | BondDescriptorAtom | Branch | StochasticObject,
#                  right_down: Atom | BondDescriptorAtom | Branch | StochasticObject,
#                  ):
#         self.type_ = type_


class Bond:
    __slots__ = ["id_", "type_", "symbol", "atom1", "atom2", "ring_id"]
    _tree_print_repr = True

    def __init__(self,
                 symbol: str,
                 atom1: Atom | BondDescriptorAtom | StochasticObject,
                 atom2: Atom | BondDescriptorAtom | StochasticObject | None = None,
                 id_: int = None,
                 ring_id: int = None
                 ):
        self.id_ = id_
        self.type_ = bond_mapping[symbol]
        self.symbol = symbol
        self.atom1 = atom1
        self.atom2 = atom2
        self.ring_id = ring_id

    def __str__(self):
        text = self.symbol
        if Config.color_output:
            return Config.colors(text, "Blue")

        return text

    def __repr__(self):
        bond_repr_symbols = {
            Atom: "A",
            BondDescriptorAtom: "BD",
            StochasticObject: "SO",
        }
        text = self.symbol + "{" + str(self.id_) + f"|{bond_repr_symbols[type(self.atom1)]}{self.atom1.id_}-"
        if self.atom2 is not None:
            text += f">{bond_repr_symbols[type(self.atom2)]}{self.atom2.id_}"
        text += "}"
        if Config.color_output:
            return Config.colors(text, "Blue")

        return text

    def __iter__(self):
        return iter((self.atom1, self.atom2))

    @property
    def bond_order(self) -> int:
        return bond_mapping[self.symbol].value


class BondDescriptorTypes(enum.Enum):
    Left = "<"
    Right = ">"
    Dollar = "$"
    Implicit = ""


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


class BondDescriptor:
    __slots__ = ["symbol", "type_", "index_", "instances", "stochastic_object", "_text"]

    def __init__(self, stochastic_object: StochasticObject, symbol: str):
        self._text = symbol  # here for fast string checks
        self.symbol, self.index_ = process_bonding_descriptor_symbol(symbol)
        self.type_ = BondDescriptorTypes(self.symbol)
        self.instances = []
        self.stochastic_object = stochastic_object

    def __str__(self):
        if self.index_ == 0 and not Config.show_bond_descriptor_zero_index:
            index = ""
        else:
            index = self.index_

        return "[" + self.symbol + str(index) + "]"

    def __repr__(self):
        return str(self)

    def is_pair(self, bd: BondDescriptor) -> bool:
        """ Returns true if symbols are <> and index match. """
        if self.index_ == bd.index_:
            if (self.type_ is BondDescriptorTypes.Left and bd.type_ is BondDescriptorTypes.Right) or \
                    (self.type_ is BondDescriptorTypes.Right and bd.type_ is BondDescriptorTypes.Left):
                return True

        return False


class BondDescriptorAtom:
    _tree_print_repr = True
    __slots__ = ["descriptor", "id_", "bond"]

    def __init__(self, bond_descriptor: BondDescriptor, id_: int = None):
        self.descriptor = bond_descriptor
        bond_descriptor.instances.append(self)
        self.id_ = id_
        self.bond = None

    def __str__(self):
        text = str(self.descriptor)
        if Config.color_output:
            return Config.colors(text, "Green")

        return text

    def __repr__(self):
        text = str(self.descriptor) + "{" + str(self.id_) + "}"

        if Config.color_output:
            return Config.colors(text, "Green")

        return text


class Branch:
    _tree_print_repr = False
    __slots__ = ["nodes", "id_", "parent"]

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None):
        self.nodes: list[Atom | Bond | Branch | StochasticObject | BondDescriptorAtom] = []
        self.id_ = id_
        self.parent = parent

    def __str__(self):
        return "(" + "".join((str(node) for node in self.nodes)) + ")"

    def __repr__(self):
        return "({" + str(self.id_) + "} " + " ".join((repr(node) for node in self.nodes)) + " {" + str(self.id_) + \
               "})"

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class StochasticFragment:
    _tree_print_repr = False
    __slots__ = ["nodes", "id_", "parent", "bonding_descriptors", "rings"]

    def __init__(self, parent: StochasticObject, id_: int = None):
        self.nodes: list[Atom | Bond | BondDescriptorAtom | Branch | StochasticObject] = []
        self.rings: list[Bond] = []
        self.bonding_descriptors: list[BondDescriptor] = []
        self.id_ = id_
        self.parent = parent

    def __str__(self):
        return "".join((str(node) for node in self.nodes))

    def __repr__(self):
        return " ".join((repr(node) for node in self.nodes))

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class StochasticObject:
    _tree_print_repr = False
    __slots__ = ["nodes", "bonding_descriptors", "id_", "parent", "end_group_left", "end_group_right",
                 "bond_left", "bond_right"]

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None):
        self.nodes: list[StochasticFragment] = []
        self.bonding_descriptors: list[BondDescriptor] = []
        self.end_group_left: BondDescriptorAtom | None = None
        self.end_group_right: BondDescriptorAtom | None = None
        self.id_ = id_
        self.bond_left: Bond | None = None
        self.bond_right: Bond | None = None
        self.parent = parent

    def __str__(self):
        if Config.color_output:
            fragments = Config.colors(",", "Red").join((str(node) for node in self.nodes))
            return Config.colors("{", "Red") + str(self.end_group_left) + fragments + str(self.end_group_right) + \
                   Config.colors("}", "Red")

        fragments = ",".join((str(node) for node in self.nodes))
        return "{" + str(self.end_group_left) + fragments + str(self.end_group_right) + "}"

    def __repr__(self):
        if Config.color_output:
            fragments = Config.colors(",", "Red").join((repr(node) for node in self.nodes))
            return Config.colors("{", "Red") + repr(self.end_group_left) + fragments + repr(self.end_group_right) + \
                   Config.colors("}", "Red")

        fragments = ",".join((repr(node) for node in self.nodes))
        return "{" + repr(self.end_group_left) + fragments + repr(self.end_group_right) + "}"

    @property
    def implicit_endgroups(self) -> bool:
        """ Returns true if one or more are implicit """
        return True if self.end_group_left.descriptor.type_ is BondDescriptorTypes.Implicit or \
                       self.end_group_right.descriptor.type_ is BondDescriptorTypes.Implicit else False

    @property
    def in_stochastic_object(self) -> bool:
        return True

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class BigSMILES:
    _tree_print_repr = False
    __slots__ = ["nodes", "atoms", "bonds", "rings", "input_text", "_tokens"]

    def __init__(self, input_text: str = None):
        self.nodes: list[Atom | Bond | StochasticObject | Branch] = []
        self.atoms: list[Atom] = []  # does count bonds in sub-objects
        self.bonds: list[Bond] = []  # does count bonds in sub-objects
        self.rings: list[Bond] = []  # does not count rings in sub-objects

        # process input string
        self.input_text = input_text
        self._tokens = []

        if input_text:
            from bigsmiles.create_parse_tree import create_parse_tree
            create_parse_tree(self)

    def __str__(self):
        return "".join((str(node) for node in self.nodes))

    def __repr__(self):
        return " ".join((repr(node) for node in self.nodes))

    def __getitem__(self, obj: int | slice) -> Atom | list[Atom]:
        return self.atoms[obj]

    def __iter__(self) -> list[Atom]:
        return self.atoms

    @property
    def in_stochastic_object(self) -> bool:
        return False

    @property
    def root(self) -> BigSMILES:
        return self

    def print_tree(self, show_object_label: bool = True, print_repr: bool = False):
        """
        prints a tree representation of the parsed bigsmiles

        Parameters
        ----------
        show_object_label: bool
            show object labels
        print_repr: bool
            use repr() instead of str()

        """
        from bigsmiles.tree_to_string import tree_to_string  # here to avoid circular imports
        print(tree_to_string(self, show_object_label, print_repr))

    def graph(self):
        from bigsmiles.nx_graph.create_nx_graph import create_nx_graph  # here to avoid circular imports and optional package requirements
        return create_nx_graph(self)
