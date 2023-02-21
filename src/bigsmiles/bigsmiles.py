from __future__ import annotations

from bigsmiles.config import Config


class Atom:
    __slots__ = ["id_", "element", "isotope", "stereo", "hydrogens", "charge", "valence",
                 "organic", "bonds", "possible_valence", "_default_valance", "__dict__"]
    _tree_print_repr = True

    def __init__(self,
                 id_: int,
                 element: str,
                 isotope: int | None = None,
                 stereo: str = '',
                 hydrogens: int = 0,
                 charge: int = 0,
                 valence: int = None,
                 **kwargs
                 ):
        self.id_ = id_
        self.element = element
        self.isotope = isotope
        self.stereo = stereo
        self.hydrogens = hydrogens
        self.charge = charge
        self.possible_valence: tuple[int] = Config.get_atom_possible_valence(element)
        if valence is None:
            self._default_valance: bool = True
            self.valence = self.possible_valence[0]
        else:
            self._default_valance: bool = False
            self.valence = valence
        self.organic = True if element in Config.organics else False

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

        self.bonds = []

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return str(self) + "{" + str(self.id_) + "}"

    @property
    def bond_capacity(self) -> int:
        return self.valence - self.hydrogens + self.charge

    @property
    def number_of_bonds(self) -> int:
        return sum((bond.bond_order for bond in self.bonds)) + self.hydrogens

    @property
    def bonds_available(self) -> int:
        return self.bond_capacity - self.number_of_bonds

    @property
    def ring_indexes(self) -> list[int]:
        ring_index = []
        for bond in self.bonds:
            if bond.ring_id is not None:
                ring_index.append(bond.ring_id)

        return ring_index

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False) -> str:
        text = self.element
        bracket_flag = False

        if self.element == "H":
            bracket_flag = True

        if self.stereo:
            text += self.stereo
            bracket_flag = True

        if show_hydrogens or self.hydrogens > 0:
            if self.hydrogens > 0 or self.bonds_available != 0:
                text += f"H{self.bonds_available + self.hydrogens}"
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

        if self.ring_indexes:
            for bond in self.bonds:
                if bond.ring_id is not None:
                    if bond.ring_id > 9:
                        ring_id = "%" + str(bond.ring_id)
                    else:
                        ring_id = str(bond.ring_id)
                    text += bond.symbol + ring_id

        return text

    def _increase_valance(self, requested_valance_increase: int = 0) -> bool:
        " If default valance, try to increase when an additional bond added that needs it. "
        if not self._default_valance:
            return False

        for valence in self.possible_valence:
            if valence > self.valence:
                old_valance = self.valence
                self.valence = valence
                if self.bonds_available < requested_valance_increase:
                    # increase was not enough for requested valance increase; set valence back to original value
                    self.valence = old_valance
                    continue

                return True

        return False


bond_mapping = {
    "": 1,
    "=": 2,
    "#": 3
}


class Bond:
    __slots__ = ["id_", "symbol", "atom1", "atom2", "ring_id", "__dict__"]
    _tree_print_repr = True

    def __init__(self,
                 symbol: str,
                 atom1: Atom | BondDescriptorAtom | StochasticObject,
                 atom2: Atom | BondDescriptorAtom | StochasticObject | None = None,
                 id_: int = None,
                 ring_id: int = None
                 ):
        self.id_ = id_
        self.symbol = symbol
        self.atom1 = atom1
        self.atom2 = atom2
        self.ring_id = ring_id

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        return Config.add_color(self.symbol, 'Blue', skip_color)

    def __iter__(self):
        return iter((self.atom1, self.atom2))

    @property
    def bond_order(self) -> int:
        return bond_mapping[self.symbol]

    @bond_order.setter
    def bond_order(self, count: int):
        for k, v in bond_mapping.items():
            if v == count:
                self.symbol = k
                return

        raise ValueError(f"Invalid bond_order. \nGiven: {count} \n Acceptable: {bond_mapping}")


class BondDescriptor:
    __slots__ = ["descriptor", "index_", "instances", "stochastic_object", "__dict__"]

    def __init__(self, stochastic_object: StochasticObject, descriptor: str, index_: int):
        self.descriptor = descriptor
        self.index_ = index_
        self.instances = []
        self.stochastic_object = stochastic_object

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        if self.index_ == 1 and not Config.show_bond_descriptor_zero_index:
            index = ""
        else:
            index = self.index_
        return "[" + self.descriptor + str(index) + "]"

    @property
    def symbol(self) -> str:
        return self.descriptor + str(self.index_)

    def is_pair(self, bd: BondDescriptor) -> bool:
        """ Returns true if symbols are <> and index match. """
        if self.index_ == bd.index_:
            if (self.descriptor == "<" and bd.descriptor == ">") or \
                    (self.descriptor == ">" and bd.descriptor == "<"):
                return True

        return False


class BondDescriptorAtom:
    _tree_print_repr = True
    __slots__ = ["descriptor", "id_", "bond", "__dict__"]

    def __init__(self, bond_descriptor: BondDescriptor, id_: int = None):
        self.descriptor = bond_descriptor
        bond_descriptor.instances.append(self)
        self.id_ = id_
        self.bond = None

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        return Config.add_color(self.descriptor.to_string(show_hydrogens, print_repr), 'Green', skip_color)


class Branch:
    _tree_print_repr = False
    __slots__ = ["nodes", "id_", "parent", "__dict__"]

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None):
        self.nodes: list[Atom | Bond | Branch | StochasticObject | BondDescriptorAtom] = []
        self.id_ = id_
        self.parent = parent

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        text = "".join(
            node.to_string(show_hydrogens, print_repr) for node in self.nodes
        )

        return "(" + text + ")"

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class StochasticFragment:
    _tree_print_repr = False
    __slots__ = ["nodes", "id_", "parent", "bonding_descriptors", "rings", "__dict__"]

    def __init__(self, parent: StochasticObject, id_: int = None):
        self.nodes: list[Atom | Bond | BondDescriptorAtom | Branch | StochasticObject] = []
        self.rings: list[Bond] = []
        self.bonding_descriptors: list[BondDescriptor] = []
        self.id_ = id_
        self.parent = parent

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        return "".join(node.to_string(show_hydrogens, print_repr, skip_color) for node in self.nodes)

    @property
    def in_stochastic_object(self) -> bool:
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class StochasticObject:
    _tree_print_repr = False
    __slots__ = ["nodes", "bonding_descriptors", "id_", "parent", "end_group_left", "end_group_right",
                 "bond_left", "bond_right", "__dict__"]

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
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        text = Config.add_color("{", "Red", skip_color)
        text += self.end_group_left.to_string(show_hydrogens, print_repr, skip_color)
        text += ",".join(node.to_string(show_hydrogens, print_repr, skip_color) for node in self.nodes)
        text += self.end_group_right.to_string(show_hydrogens, print_repr, skip_color)
        return text + Config.add_color("}", "Red", skip_color)

    @property
    def implicit_endgroups(self) -> bool:
        """ Returns true if one or more are implicit """
        return True if self.end_group_left.descriptor.descriptor == "" or \
                       self.end_group_right.descriptor.descriptor == "" else False

    @property
    def in_stochastic_object(self) -> bool:
        return True

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


def contains_stochastic_object(nodes: Atom | Bond | StochasticObject | Branch):
    for node in nodes:
        if isinstance(node, StochasticObject):
            return True
        if isinstance(node, Branch):
            result = contains_stochastic_object(node.nodes)
            if result:
                return True

    return False


class BigSMILES:
    _tree_print_repr = False
    __slots__ = ["nodes", "atoms", "bonds", "rings", "__dict__"]

    def __init__(self, input_text: str = None):
        self.nodes: list[Atom | Bond | StochasticObject | Branch] = []
        self.atoms: list[Atom] = []  # includes atoms in sub-objects
        self.bonds: list[Bond] = []  # includes bonds in sub-objects
        self.rings: list[Bond] = []  # does not include rings in sub-objects

        # parse input string
        if input_text:
            from bigsmiles.parse_bigsmiles_str import parse_bigsmiles_str
            parse_bigsmiles_str(input_text, self)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __getitem__(self, obj: int | slice) -> Atom | list[Atom]:
        return self.atoms[obj]

    def __iter__(self) -> list[Atom]:
        return self.atoms

    def __bool__(self):
        if self.nodes:
            return True
        else:
            return False

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False) -> str:
        return "".join(node.to_string(show_hydrogens, print_repr, skip_color) for node in self.nodes)

    @property
    def in_stochastic_object(self) -> bool:
        return False

    @property
    def contains_stochastic_object(self) -> bool:
        return contains_stochastic_object(self.nodes)

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
        from bigsmiles.graph.create_graph import create_graph  # here to avoid circular imports and optional package requirements
        return create_graph(self)