from __future__ import annotations
import enum

from bigsmiles.config import Config


def to_repr(obj) -> str:
    return repr(obj)


def to_str(obj) -> str:
    return str(obj)


class Atom:
    __slots__ = ["id_", "element", "bond_symbol", "isotope", "stereo", "hydrogens", "charge", "valence",
                 "organic", "bonds"]
    _tree_print_repr = True

    def __init__(self,
                 id_: int,
                 element: str,
                 isotope: int | None = None,
                 stereo: str = '',
                 hydrogens: int = 0,
                 charge: int = 0,
                 valance: int = None,
                 **kwargs
                 ):
        self.id_ = id_
        self.element = element
        self.isotope = isotope
        self.stereo = stereo
        self.hydrogens = hydrogens
        self.charge = charge
        self.valence = valance
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
    def bond_capacity(self) -> int | None:
        if self.valence is None:
            return None
        return self.valence - self.hydrogens + self.charge

    @property
    def number_of_bonds(self) -> int:
        return sum((bond.type_.value for bond in self.bonds)) + self.hydrogens

    @property
    def bonds_available(self) -> int | None:
        if self.valence is None:
            return None
        return self.bond_capacity - self.number_of_bonds

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
        text = self.element
        bracket_flag = False

        if self.element == "H":
            bracket_flag = True

        if self.stereo:
            text += self.stereo
            bracket_flag = True

        if (show_hydrogens and self.bonds_available is not None and self.bonds_available > 0) or \
                self.hydrogens is not None:
            if self.hydrogens is not None:
                text += f"H{self.hydrogens if self.hydrogens > 1 else ''}"
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
            return Config.add_color(text, "Blue")

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
            return Config.add_color(text, "Blue")

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


class BondDescriptor:
    __slots__ = ["descriptor", "type_", "index_", "instances", "stochastic_object"]

    def __init__(self, stochastic_object: StochasticObject, descriptor: str, index_: int):
        self.descriptor = descriptor
        self.index_ = index_
        self.type_ = BondDescriptorTypes(self.descriptor)
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

    @property
    def symbol(self) -> str:
        return self.descriptor + str(self.index_)

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
            return Config.add_color(text, "Green")

        return text

    def __repr__(self):
        text = str(self.descriptor) + "{" + str(self.id_) + "}"

        if Config.color_output:
            return Config.add_color(text, "Green")

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

    # def __str__(self):
    #     return self.to_string()
    #
    # def __repr__(self):
    #     return self.to_string(print_repr=True, skip_color=True)
    #
    # def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
    #     text = "".join(
    #         Config.add_color(node.to_string(show_hydrogens, print_repr), 'Red', skip_color) for node in self.nodes
    #     )
    #
    #     return text

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
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        text = Config.add_color("{", "Red", skip_color)
        text += self.end_group_left.to_string(show_hydrogens, print_repr, skip_color)
        text += "".join(
            Config.add_color(node.to_string(show_hydrogens, print_repr, skip_color), 'Red', skip_color) for node in self.nodes
        )
        text += self.end_group_right.to_string(show_hydrogens, print_repr, skip_color)

        return text + Config.add_color("}", "Red", skip_color)

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
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __getitem__(self, obj: int | slice) -> Atom | list[Atom]:
        return self.atoms[obj]

    def __iter__(self) -> list[Atom]:
        return self.atoms

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False) -> str:
        return "".join(node.to_string(show_hydrogens, print_repr, skip_color) for node in self.nodes)

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
