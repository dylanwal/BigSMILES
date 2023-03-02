from __future__ import annotations
import logging

from bigsmiles.errors import ConstructorError
import bigsmiles.chemical_data as chemical_data
from bigsmiles.config import Config


class Atom:
    __slots__ = ["id_", "element", "isotope", "stereo", "hydrogens", "charge", "valence", "_valene_warning_raised"
                 "organic", "_bonds", "possible_valence", "_default_valence", "__dict__", "parent"]
    _tree_print_repr = True

    def __init__(self,
                 id_: int,
                 element: str,
                 isotope: int | None = None,
                 stereo: str = '',
                 hydrogens: int | None = None,
                 charge: int = 0,
                 valence: int = None,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs
                 ):
        self.id_ = id_
        self.element = element
        self.isotope = isotope
        self.stereo = stereo
        self.hydrogens = hydrogens
        self.charge = charge
        self.possible_valence: tuple[int] = chemical_data.get_atom_possible_valence(element)
        if valence is None:
            self._default_valence: bool = True
            self.valence = self.possible_valence[0]
        else:
            self._default_valence: bool = False
            self.valence = valence
        self.organic = True if element in chemical_data.organics else False

        self.parent = parent

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

        self._bonds = []
        self._valene_warning_raised = False

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return str(self) + "{" + str(self.id_) + "}"

    @property
    def bonds(self) -> list[Bond]:
        return self._bonds

    @bonds.setter
    def bonds(self, bonds: list[Bond]):
        self._bonds = bonds
        # Validation for available valence
        if self.bonds_available < 0:
            if not self._increase_valence(self.bonds_available):
                logging.error(f"Too many bonds trying to be made. {str(self)}")

    @property
    def implicit_hydrogens(self) -> int:
        """ Number of implicit hydrogens. zero if explict hydrogens specified. """
        return self.bonds_available if self.hydrogens is None else 0

    @property
    def bond_capacity(self) -> int:
        """ Total capacity of the atom. """
        return self.valence + self.charge

    @property
    def number_of_bonds(self) -> int:
        """ Number of bonds already formed. Not including implicit hydrogens. """
        num_bonds = sum((bond.bond_order for bond in self.bonds))
        if self.hydrogens is not None:
            num_bonds += self.hydrogens
        return num_bonds

    @property
    def bonds_available(self) -> int:
        """ Number bonds that remain open for bonding. Will reduce implicit hydrogen count. """
        bonds_available = self.bond_capacity - self.number_of_bonds
        if bonds_available < 0:
            return 0
        return bonds_available

    @property
    def full_valence(self) -> bool:
        """ Returns true if the atoms valence is full"""
        if self.valence - self.number_of_bonds - self.implicit_hydrogens == 0:
            return True

        return False

    @property
    def ring_indexes(self) -> list[int]:
        ring_index = []
        for bond in self.bonds:
            if bond.ring_id is not None:
                ring_index.append(bond.ring_id)

        return ring_index

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False) -> str:
        text = self._to_string(show_hydrogens)
        if not self.full_valence and not self._valene_warning_raised:
            self._valene_warning_raised = True
            logging.warning(f"Incomplete valence detected on atom: {text}")

        return text

    def _to_string(self, show_hydrogens: bool = False) -> str:
        text = self.element
        bracket_flag = False

        if self.element == "H":
            bracket_flag = True

        if self.stereo:
            text += self.stereo
            bracket_flag = True

        if show_hydrogens or self.hydrogens is not None:
            if self.hydrogens is not None:
                if self.hydrogens == 1:
                    text += f"H"
                elif self.hydrogens > 1:
                    text += f"H{self.hydrogens}"

                bracket_flag = True
            else:
                if self.implicit_hydrogens == 1:
                    text += f"H"
                    bracket_flag = True
                if self.implicit_hydrogens > 1:
                    text += f"H{self.implicit_hydrogens}"
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

    def _increase_valence(self, requested_valence_increase: int = 0) -> bool:
        """ If default valence, try to increase when an additional bond added that needs it. And does the increase!! """
        if not self._default_valence:
            return False

        for valence in self.possible_valence:
            if valence > self.valence:
                old_valence = self.valence
                self.valence = valence
                if self.bonds_available < requested_valence_increase:
                    # increase was not enough for requested valence increase; set valence back to original value
                    self.valence = old_valence
                    continue

                return True

        return False

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


bond_mapping = {
    None: 0,
    "": 1,
    "=": 2,
    "#": 3
}


class Bond:
    __slots__ = ["id_", "symbol", "atom1", "atom2", "ring_id", "__dict__", "parent"]
    _tree_print_repr = True

    def __init__(self,
                 symbol: str,
                 atom1: Atom | BondDescriptorAtom | StochasticObject,
                 atom2: Atom | BondDescriptorAtom | StochasticObject | None = None,
                 id_: int = None,
                 ring_id: int = None,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs
                 ):
        self.id_ = id_
        self.symbol = symbol
        self.atom1 = atom1
        self.atom2 = atom2
        self.ring_id = ring_id
        self.parent = parent

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        return Config.add_color(self.symbol, 'Blue', skip_color)

    def __iter__(self):
        return iter((self.atom1, self.atom2))

    def __reversed__(self):
        return iter((self.atom2, self.atom1))

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

    @property
    def root(self) -> BigSMILES:
        return self.parent.root


class BondDescriptor:
    __slots__ = ["parent", "descriptor", "index_", "_instances", "_instances_up_to_date", "_bond_symbol", "__dict__"]

    def __init__(self, parent: StochasticObject, descriptor: str, index_: int, **kwargs):
        self.parent = parent
        self.descriptor = descriptor
        self.index_ = index_

        self._instances_up_to_date: bool = True
        self._instances: list[BondDescriptorAtom] = []
        self._bond_symbol: str | None = None

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

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

    @property
    def root(self) -> BigSMILES:
        return self.parent.root

    @property
    def instances(self) -> list[BondDescriptorAtom]:
        return self._instances

    @instances.setter
    def instances(self, instances: list[BondDescriptor]):
        self._instances = instances
        self._instances_up_to_date = False

    @property
    def bond_symbol(self) -> str:
        if not self._instances_up_to_date:
            self._get_bond_symbol()
            self._instances_up_to_date = True

        return self._bond_symbol

    def _get_bond_symbol(self):
        for instance_ in self.instances:
            if instance_.bond_symbol is not None:
                if self._bond_symbol is not None and self._bond_symbol != instance_.bond_symbol:
                    logging.warning("Multiple bond orders to same bonding descriptor.")
                else:
                    self._bond_symbol = instance_.bond_symbol

    @property
    def bond_order(self) -> int:
        return bond_mapping[self.bond_symbol]

    @property
    def implicit(self) -> bool:
        if self.descriptor == "":
            return True
        return False


class BondDescriptorAtom:
    _tree_print_repr = True
    __slots__ = ["descriptor", "id_", "_bond", "__dict__", "parent"]

    def __init__(self,
                 bond_descriptor: BondDescriptor,
                 id_: int = None,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs):
        self.descriptor = bond_descriptor
        bond_descriptor.instances += [self]
        self.id_ = id_
        self._bond = None  # Added after construction
        self.parent = parent

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        return Config.add_color(self.descriptor.to_string(show_hydrogens, print_repr), 'Green', skip_color)

    @property
    def bond(self) -> Bond:
        return self._bond

    @bond.setter
    def bond(self, bond: Bond):
        if self.bond is not None:
            raise ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond = bond

    @property
    def root(self) -> BigSMILES:
        return self.parent.root

    @property
    def bond_symbol(self) -> str | None:
        if self.bond is None:
            return None
        return self.bond.symbol

    @property
    def bond_order(self) -> 0:
        if self.bond is None:
            return 0
        return self.bond.bond_order


class Branch:
    _tree_print_repr = False
    __slots__ = ["nodes", "id_", "parent", "__dict__"]

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None, **kwargs):
        self.nodes: list[Atom | Bond | Branch | StochasticObject | BondDescriptorAtom] = []
        self.id_ = id_
        self.parent = parent

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

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

    def _get_id(self) -> int:
        return self.root._get_id()


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

    def _get_id(self) -> int:
        return self.root._get_id()


class StochasticObject:
    _tree_print_repr = False
    __slots__ = ["nodes", "bonding_descriptors", "id_", "parent", "bd_left", "bd_right",
                 "_bond_left", "_bond_right", "__dict__"]

    def __init__(self, parent: BigSMILES | StochasticFragment | Branch, id_: int = None, **kwargs):
        self.nodes: list[StochasticFragment] = []
        self.bonding_descriptors: list[BondDescriptor] = []
        self.bd_left: BondDescriptor | None = None
        self.bd_right: BondDescriptor | None = None
        self.id_ = id_
        self._bond_left: Bond | None = None
        self._bond_right: Bond | None = None
        self.parent = parent

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    @property
    def bond_left(self) -> Bond | None:
        return self._bond_left

    @bond_left.setter
    def bond_left(self, bond: Bond):
        if self._bond_left is not None:
            raise ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond_left = bond

    @property
    def bond_right(self) -> Bond | None:
        return self._bond_right

    @bond_right.setter
    def bond_right(self, bond: Bond):
        if self._bond_right is not None:
            raise ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond_right = bond

    def to_string(self, show_hydrogens: bool = False, print_repr: bool = False, skip_color: bool = False):
        text = ""
        text += Config.add_color("{", "Red", skip_color)
        text += self.bd_left.to_string(show_hydrogens, print_repr, skip_color) if self.bd_left is not None else ""
        text += ",".join(node.to_string(show_hydrogens, print_repr, skip_color) for node in self.nodes)
        text += self.bd_right.to_string(show_hydrogens, print_repr, skip_color) if self.bd_right is not None else ""
        text += Config.add_color("}", "Red", skip_color)
        if self.bond_right is not None and self.bond_right.ring_id is not None:
            text += self.bond_right.symbol + str(self.bond_right.ring_id)
        return text

    @property
    def implicit_endgroups(self) -> bool:
        """ Returns true if one or more are implicit """
        if (self.bd_left is not None and self.bd_left.implicit) or \
                (self.bd_right is not None and self.bd_right.implicit):
            return True

        return False

    @property
    def in_stochastic_object(self) -> bool:
        return True

    @property
    def root(self) -> BigSMILES:
        return self.parent.root

    @property
    def bonds(self) -> list[Bond]:
        bonds = []
        if self.bond_left is not None:
            bonds.append(self.bond_left)
        if self.bond_right is not None:
            bonds.append(self.bond_right)
        return bonds

    def _get_id(self) -> int:
        return self.root._get_id()


def contains_stochastic_object(nodes: list[Atom, Bond, Branch, StochasticObject]):
    """ recursive search for a stochastic_object. """
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
        self.atoms: list[Atom | BondDescriptorAtom] = []  # includes atoms in sub-objects
        self.bonds: list[Bond] = []  # includes bonds in sub-objects
        self.rings: list[Bond] = []  # does not include rings in sub-objects

        self.__id = 0

        # parse input string
        if input_text:
            # import here to avoid circular imports
            from bigsmiles.constructors.parse_bigsmiles_str import parse_bigsmiles_str
            parse_bigsmiles_str(input_text, self)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __getitem__(self, obj: int | slice) -> Atom | list[Atom]:
        return self.atoms[obj]

    def __iter__(self):
        return iter(self.atoms)

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

    def _get_id(self) -> int:
        id_ = self.__id
        self.__id += 1
        return id_

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


# types
has_node_attr = BigSMILES | Branch | StochasticObject | StochasticFragment
has_ring_attr = BigSMILES | StochasticFragment
has_parent_attr = Branch | StochasticObject | StochasticFragment | Bond | BondDescriptorAtom | Atom
has_root_attr = BigSMILES | has_parent_attr
