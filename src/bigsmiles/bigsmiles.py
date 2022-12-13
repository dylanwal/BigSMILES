from __future__ import annotations
import enum
import re

from bigsmiles.config import Config
from bigsmiles.tokenizer import Token, tokenize, TokenKind


class BigSMILESError(Exception):
    def __init__(self, text: str, optional: str = None):
        if optional is not None:
            text += "(" + optional + ")"
        super().__init__(text)


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
    result = re.search(r'H[\d]?', symbol_text)
    if result is None:
        return symbol_text, None

    if result.group()[-1].isdigit():
        count = int(result.group()[-1])
    else:
        count = 1

    return symbol_text[:result.span()[0]] + symbol_text[result.span()[1]:], count


def check_atom_symbol(symbol_text: str) -> tuple[str, tuple[int] | None, AtomChirality, int, int | None, int | None]:
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
    _tree_print_label = True

    def __init__(self, symbol: str, id_: int = None):
        self.id_ = id_
        self.symbol, self.valance_possible, self.chiral, self.charge, self.isotope, self._hydrogens = \
            check_atom_symbol(symbol)
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
            if bond.ring is not None:
                ring_index.append(bond.ring)

        return ring_index

    def _add_bond(self, bond: Bond):
        if self.bonds_available is not None and self.bonds_available == 0:
            if self.valance == self.valance_possible[-1] or self.valance_possible[-1] + bond.value:
                raise BigSMILESError("Too many bonds trying to be made.", repr(self))

            for val in self.valance_possible:
                if val > self.valance:
                    self.valance = val
                    break

        self.bonds.append(bond)

    def to_string(self, show_hydrogens: bool = False) -> str:
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


class Bond:
    __slots__ = ["id_", "type_", "symbol", "_atom1", "_atom2", "ring"]
    _tree_print_label = True

    def __init__(self,
                 symbol: str,
                 atom1: Atom | BondDescriptor | StochasticObject,
                 atom2: Atom | BondDescriptor | StochasticObject | None = None,
                 ring: int = None,
                 id_: int = None
                 ):
        self.id_ = id_
        self.type_ = bond_mapping[symbol]
        self.symbol = symbol
        self._atom1 = None
        self._atom2 = None
        self.atom1 = atom1
        self.atom2 = atom2
        self.ring = ring

    def __str__(self):
        text = self.symbol
        if Config.color_output:
            return Config.colors(text, "Blue")

        return text

    def __repr__(self):
        bond_repr_symbols = {
            Atom: "A",
            BondDescriptor: "BD",
            StochasticObject: "SO",
        }
        text = self.symbol + "{" + str(self.id_) + f"|{bond_repr_symbols[type(self.atom1)]}{self.atom1.id_}-" \
                                                   f">{bond_repr_symbols[type(self.atom2)]}{self.atom2.id_}" + "}"
        if Config.color_output:
            return Config.colors(text, "Blue")

        return text

    @property
    def value(self) -> int:
        return self.type_.value

    @property
    def atom1(self) -> Atom:
        return self._atom1

    @atom1.setter
    def atom1(self, atom: Atom | BondDescriptor):
        atom._add_bond(self)
        self._atom1 = atom

    @property
    def atom2(self) -> Atom:
        return self._atom2

    @atom2.setter
    def atom2(self, atom: Atom | BondDescriptor | None):
        if atom is not None:
            atom._add_bond(self)
            self._atom2 = atom


class BondDescriptorTypes(enum.Enum):
    Left = "<"
    Right = ">"
    Dollar = "$"


def process_bonding_descriptor_symbol(symbol: str) -> tuple[str, int]:
    symbol = symbol.replace("[", "").replace("]", "")
    if symbol[-1].isdigit():
        return symbol[0], int(symbol[-1])

    return symbol, 0


class BondDescriptor:
    __slots__ = ["id_", "symbol", "type_", "index_", "bond"]
    _tree_print_label = True

    def __init__(self, symbol: str, id_: int = None):
        self.id_ = id_
        self.symbol, self.index_ = process_bonding_descriptor_symbol(symbol)
        self.type_ = BondDescriptorTypes(self.symbol)
        self.bond = None

    def __str__(self):
        if self.index_ == 0 and not Config.show_bond_descriptor_zero_index:
            index = ""
        else:
            index = self.index_
        text = "[" + self.symbol + str(index) + "]"
        if Config.color_output:
            return Config.colors(text, "Green")

        return text

    def __repr__(self):
        text = "[" + str(self.symbol) + "]" + "{" + str(self.id_) + "}"

        if Config.color_output:
            return Config.colors(text, "Green")

        return text

    def _add_bond(self, bond: Bond):
        self.bond = bond


class GetParentID:
    _tree_print_label = False

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch | StochasticObject):
        self.parent = parent
        if isinstance(parent, BigSMILES):
            self.root = parent
        else:
            self.root = parent.root

    def _get_atom_id(self) -> int:
        return self.root._get_atom_id()

    def _get_bond_id(self) -> int:
        return self.root._get_bond_id()

    def _get_bond_descriptor_id(self) -> int:
        return self.root._get_bond_descriptor_id()

    def _get_branch_id(self) -> int:
        return self.root._get_branch_id()

    def _get_stochastic_fragment_id(self) -> int:
        return self.root._get_stochastic_fragment_id()

    def _get_stochastic_object_id(self) -> int:
        return self.root._get_stochastic_object_id()


class Branch(GetParentID):
    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None):
        self.nodes: list[Atom | Bond | Branch | StochasticObject | BondDescriptor] = []
        self.id_ = id_
        super().__init__(parent)

    def __str__(self):
        return "(" + "".join((str(node) for node in self.nodes)) + ")"

    def __repr__(self):
        return "({" + str(self.id_) + "} " + " ".join((repr(node) for node in self.nodes)) + " {" + str(self.id_) + \
               "})"

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object

    def _get_ring(self, id_: int) -> Bond | None:
        return self.root._get_ring(id_)

    def _add_ring(self, ring: Bond):
        self.root.rings.append(ring)


class StochasticFragment(GetParentID):

    def __init__(self, parent: StochasticObject, id_: int = None):
        self.nodes: list[Atom | Bond | BondDescriptor | Branch | StochasticObject] = []
        self.id_ = id_
        super().__init__(parent)

    def __str__(self):
        return "".join((str(node) for node in self.nodes))

    def __repr__(self):
        return " ".join((repr(node) for node in self.nodes))

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object


class StochasticObject(GetParentID):

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None):
        self.nodes: list[StochasticFragment] = []
        self.end_group_left = None
        self.end_group_right = None
        self.id_ = id_
        self.bond_left = None
        self.bond_right = None
        super().__init__(parent)

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
        return False  # TODO:

    @property
    def in_stochastic_object(self) -> bool:
        return True

    def _add_bond(self, bond: Bond):
        if self.bond_left is None:
            self.bond_left = bond
        else:
            self.bond_right = bond


class BigSMILES:
    _tree_print_label = False

    def __init__(self, input_text: str):
        self.nodes: list[Atom | Bond | StochasticObject | Branch] = []
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []
        self.rings: list[Bond] = []

        # setup id counters
        self._atom_counter = 0
        self._bond_counter = 0
        self._bond_descriptor = 0
        self._branch_counter = 0
        self._stochastic_fragment = 0
        self._stochastic_object = 0

        # process input string
        self.input_text = input_text
        self._tokens = tokenize(input_text)

        import bigsmiles.tokens_to_objects
        bigsmiles.tokens_to_objects.tokens_to_objects(self, self._tokens)

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

    def _get_atom_id(self) -> int:
        self._atom_counter += 1
        return self._atom_counter - 1

    def _get_bond_id(self) -> int:
        self._bond_counter += 1
        return self._bond_counter - 1

    def _get_bond_descriptor_id(self) -> int:
        self._bond_descriptor += 1
        return self._bond_descriptor - 1

    def _get_branch_id(self) -> int:
        self._branch_counter += 1
        return self._branch_counter - 1

    def _get_stochastic_fragment_id(self) -> int:
        self._stochastic_object += 1
        return self._stochastic_object - 1

    def _get_stochastic_object_id(self) -> int:
        self._stochastic_object += 1
        return self._stochastic_object - 1

    def _get_ring(self, id_: int) -> Bond | None:
        for ring in self.rings:
            if ring.ring == id_:
                return ring

        return None

    def _add_ring(self, ring: Bond):
        self.rings.append(ring)

    def add_atom(self):
        pass

    def add_bond(self):
        pass

    def add_bond_atom_pair(self):
        pass

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
        from bigsmiles.tree_to_string import tree_to_string
        print(tree_to_string(self, show_object_label, print_repr))
