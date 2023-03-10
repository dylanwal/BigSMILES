"""

The following contains the following data structures:
* Atom
* Bond
* Branch
* Bond descriptor
* Bond descriptor atom
* Stochastic fragment
* Stochastic object
* BigSMILES

"""

from __future__ import annotations
import logging

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.config import Config

_conjugated_warning = True


class Atom:
    """
    this class represents an atom
    """
    __slots__ = ["id_", "symbol", "isotope", "stereo", "hydrogens", "charge", "class_", "organic",
                 "aromatic", "valence", "possible_valence", "_default_valence", "_valene_warning_raised", "_bonds",
                 "parent", "__dict__"]
    _eq_attr = ("id_", "symbol", "isotope", "stereo", "hydrogens", "charge", "valence", "aromatic")

    def __init__(self,
                 id_: int,
                 symbol: str,
                 isotope: int | None = None,
                 stereo: str | None = None,
                 hydrogens: int | None = None,
                 charge: int = 0,
                 valence: int | float | None = None,
                 class_: int | None = None,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs
                 ):
        """
        Attributes
        ----------
        id_: int
            id of atom (id is limited to atoms). Range: [1, inf]
        symbol: str
            symbol symbol (e.g., H, C, O, Zn)
        isotope: int | None
            isotope (e.g., [13C])
        stereo: str | None
            stereochemistry [None, "@", "@@"] (None means not defined)
        hydrogens: int | None
            number of explict hydrogens  (e.g., [CH2])
        charge: int
            symbol charge  (e.g., [Fe+3])
        valence: int
            The capacity to form bonds with other atoms
        class_: int | None
            index of class (e.g., [C:1] class_ = 1)
        parent: BigSMILES | Branch | StochasticFragment | None
            the owner of the atom
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
        self.id_ = id_
        self.symbol = symbol.capitalize()
        self.isotope = isotope
        self.stereo = stereo
        self.hydrogens = hydrogens
        self.charge = charge
        self.class_ = class_
        self.organic = True if self.symbol in chemical_data.organic_elements else False
        # TODO: calculate aromatic
        self.aromatic = True if symbol in chemical_data.aromatic_elements else False
        self.possible_valence: tuple[int] = chemical_data.atom_valences[self.symbol]
        if valence is None:
            self._default_valence: bool = True
            self.valence = self.possible_valence[0]
        else:
            self._default_valence: bool = False
            self.valence = valence

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

    def __eq__(self, other: Atom):
        """
        equality is based on the following parameters:

        * id_
        * symbol
        * isotope
        * stereo
        * hydrogens
        * charge
        * valence
        * aromatic
        * bonds ~ bond.id_
        """
        if not isinstance(other, Atom):
            return False

        # check attributes
        for attr in self._eq_attr:
            if getattr(self, attr) != getattr(other, attr):
                return False

        # check atoms
        for bond, other_bond in zip(self.bonds, other.bonds):
            if bond.id_ != other_bond.id_:
                return False

        return True

    @property
    def bonds(self) -> list[Bond]:
        """ list of bonds (excludes implicit hydrogens) """
        return self._bonds

    @bonds.setter
    def bonds(self, bonds: list[Bond]):
        self._bonds = bonds
        # Validation for available valence
        if self.number_of_bonds > self.bond_capacity:
            if not self._increase_valence(self.number_of_bonds - self.bond_capacity):
                logging.error(f"Too many bonds trying to be made. {str(self)}")

    def delete_bond(self, delete_bond: Bond | list[Bond]):
        """
        delete bond or bonds

        Parameters
        ----------
        delete_bond: Bond | list[Bond]

        """
        if not isinstance(delete_bond, list):
            delete_bond = list(delete_bond)

        for bond in delete_bond:
            try:
                self.bonds.remove(bond)
            except ValueError:
                raise ValueError("Bond not connect to atom so it can't be deleted."
                                 f"\n Atom:{self.details}\nBond attempting to remove: {bond.details}")

    @property
    def implicit_hydrogens(self) -> float:
        """ number of implicit hydrogens. (zero if explict hydrogens have been specified) """
        return self.bonds_available if self.hydrogens is None else 0

    @property
    def bond_capacity(self) -> float:
        """ total capacity of the atom account for charge """
        return self.valence + self.charge

    @property
    def number_of_bonds(self) -> float:
        """ number of bonds already formed (not including implicit hydrogens; but including explicit hydrogens) """
        num_bonds = sum((bond.bond_order for bond in self.bonds))
        if self.hydrogens is not None:
            num_bonds += self.hydrogens
        return num_bonds

    @property
    def bonds_available(self) -> float:
        """ number bonds that remain open for bonding """
        bonds_available = self.bond_capacity - self.number_of_bonds
        if bonds_available < 0:
            return 0
        return bonds_available

    @property
    def full_valence(self) -> bool:
        """ returns true if the atom valence is full """
        if self.valence - self.number_of_bonds - self.implicit_hydrogens - abs(self.charge) == 0:
            return True

        return False

    @property
    def ring_indexes(self) -> list[int]:
        """ list of ring index """
        ring_index = []
        for bond in self.bonds:
            if bond.ring_id is not None:
                ring_index.append(bond.ring_id)

        return ring_index

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        for attr in self._eq_attr:
            text += f"{attr}: {getattr(self, attr)}, "
        return text[:-2] + "}"

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        text = self._to_string(show_hydrogens, show_atom_index)
        if not self.full_valence and not self._valene_warning_raised:
            self._valene_warning_raised = True
            logging.warning(f"Incomplete valence detected on atom: {text}")

        return text

    def _to_string(self, show_hydrogens: bool = False, show_atom_index: bool = True) -> str:
        text = ""
        bracket_flag = False

        if self.symbol == "H":
            bracket_flag = True

        if self.isotope is not None:
            text += str(self.isotope)
            bracket_flag = True

        if self.aromatic:
            text += self.symbol.lower()
        else:
            text += self.symbol

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
            text += f"-{self.charge if self.charge > 1 else ''}"
            bracket_flag = True

        if self.class_ is not None and show_atom_index:
            text += ":" + str(self.class_)

        if bracket_flag:
            text = "[" + text + "]"

        if self.ring_indexes:
            # get ring ids
            rings = []
            for bond in self.bonds:
                if bond.ring_id is not None:
                    rings.append([bond.ring_id, bond])

            # sort to numerical order
            rings.sort(key=lambda x: x[0])
            for ring_id, bond in rings:
                ring_text = ""
                if bond.symbol != ":":  # don't show aromatic bonds in ring index
                    if not self.root.contains_stochastic_object and not Config.show_multi_bonds_on_both_ring_index \
                            and bond.atom2 == self:
                        pass
                    else:
                        ring_text += bond.symbol
                if ring_id > 9:
                    ring_text += "%"

                text += ring_text + str(ring_id)

        return text

    def _increase_valence(self, requested_valence_increase: float = 0) -> bool:
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
        """ the owner at the top of the parent tree """
        return self.parent.root


class Bond:
    """
    this class represents a bond
    """
    __slots__ = ["id_", "symbol", "atom1", "atom2", "ring_id", "parent", "__dict__", "double_bond_stereo"]
    _eq_attr = ("id_", "symbol", "ring_id")

    def __init__(self,
                 id_: int,
                 symbol: str,
                 atom1: Atom | BondDescriptorAtom | StochasticObject,
                 atom2: Atom | BondDescriptorAtom | StochasticObject | None = None,
                 ring_id: int | None = None,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs
                 ):
        """

        Parameters
        ----------
        id_: int
            id of bond (id is limited to bonds). Range: [1, inf]
        symbol: str
            bond symbol (e.g. '', '=', '#')
        atom1: Atom | BondDescriptorAtom | StochasticObject
            atom that the bond starts at
        atom2: Atom | BondDescriptorAtom | StochasticObject | None
            atom that the bond ends at
        ring_id: int | None
            if the bond is closing a ring, this is the index for that ring
        parent: BigSMILES | Branch | StochasticFragment | None
            the owner of the bond
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
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

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        if self.symbol == ":" and not Config.show_aromatic_bond:
            return ""
        return Config.add_color(self.symbol, 'Blue', skip_color)

    def __iter__(self):
        """ iterates over Atoms """
        return iter((self.atom1, self.atom2))

    def __reversed__(self):
        return iter((self.atom2, self.atom1))

    def __eq__(self, other: Bond):
        """
        equality is based on the following parameters:

        * id_
        * symbol
        * ring_id
        * atom1.id_
        * atom2.id_
         """
        if not isinstance(other, Bond):
            return False

        # check attributes
        for attr in self._eq_attr:
            if getattr(self, attr) != getattr(other, attr):
                return False

        # check atoms
        if self.atom1.id_ != other.atom1.id_ or self.atom2.id_ != other.atom2.id_:
            return False

        return True

    def delete(self):
        """ deletes bond and removes it from atoms and root """
        self.atom1.delete_bond(self)
        self.atom2.delete_bond(self)
        self.parent.nodes.remove(self)
        if self.parent is not self.root:
            self.root.bonds.remove(self)
        del self

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        for attr in self._eq_attr:
            text += f"{attr}: {getattr(self, attr)}, "
        text += f"bond_order: {self.bond_order}, "
        text += "atoms: " + repr(self.atom1) + " <--> " + repr(self.atom2)
        return text + "}"

    @property
    def bond_order(self) -> int:
        """  number of covalent bonds represented """
        return chemical_data.bond_mapping[self.symbol]

    @bond_order.setter
    def bond_order(self, count: int):
        for k, v in chemical_data.bond_mapping.items():
            if v == count:
                self.symbol = k
                return

        raise ValueError(f"Invalid bond_order. \nGiven: {count} \n Acceptable: {chemical_data.bond_mapping}")

    @property
    def root(self) -> BigSMILES:
        """ the owner at the top of the parent tree """
        return self.parent.root

    @property
    def double_bond_ez(self) -> str | None:
        """ returns labels for geometric isomers of alkenes ["E", "Z", None] """
        return None
        if self.symbol != "=":
            return None

        # check if atoms have bonds with "/" or "\"
        left_stereo_bond = len([True for bond in self.atom1.bonds if bond.symbol in chemical_data.stereo_bonds])
        right_stereo_bond = len([True for bond in self.atom2.bonds if bond.symbol in chemical_data.stereo_bonds])

        if left_stereo_bond == 0 and right_stereo_bond == 0:
            return None
        if left_stereo_bond == 1 and right_stereo_bond == 1:
            from bigsmiles.data_structures.stereo_rules import get_double_bond_ez
            return get_double_bond_ez(self)

        raise errors.BigSMILESError("Only one double bond stereochemistry detected.")

    @property
    def aromatic(self) -> bool:
        """ Limited accuracy """
        if self.symbol == ":":
            return True
        return False

    # @property
    # def conjugated(self) -> bool:
    #     # TODO: small molecules and bridge across stochastic object and stochastic fragments
    #     # only for
    #
    #     global _conjugated_warning  # so it is only raised once
    #     if _conjugated_warning and self.root.contains_stochastic_object:
    #
    #         logging.warning("May be wrong for bonds bridging stochastic objects and stochastic fragments.")
    #         _conjugated_warning = False
    #
    #     if self.bond_order <= 1:
    #         return False
    #
    #     if self.aromatic:
    #         return True
    #
    #     # a depth of two graph traversal looking for a double or triple bond
    #     for atom in self:
    #         for bond in atom.bonds:
    #             if bond.id_ == self.id_:
    #                 continue  # skip self
    #             for atom_ in bond:
    #                 if atom_.id_ == atom.id_:
    #                     continue  # skip if past atom
    #                 for bond_ in atom_.bonds:
    #                     if bond_.bond_order > 1:
    #                         return True
    #
    #     return False


class BondDescriptor:
    """
    this class represents a bonding descriptor
    """
    __slots__ = ["parent", "descriptor", "index_", "_instances", "_instances_up_to_date", "bond_symbol", "__dict__"]
    _eq_attr = ("descriptor", "index_", "bond_symbol")

    def __init__(self,
                 parent: StochasticObject,
                 descriptor: str,
                 index_: int,
                 bond_symbol: str | None,
                 **kwargs
                 ):
        """

        Parameters
        ----------
        parent: StochasticObject
            the owner of the bonding descriptor
        descriptor: str
            descriptor symbol [$, >, <]
        index_: int
            descriptor index [1, inf]
        bond_symbol: str | None
            bond symbol
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
        self.parent = parent
        self.descriptor = descriptor
        self.index_ = index_

        self._instances_up_to_date: bool = True
        self._instances: list[BondDescriptorAtom] = []
        self.bond_symbol: str | None = bond_symbol

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __eq__(self, other: BondDescriptor):
        """
        equality is based on the following parameters:

        * descriptor
        * index_
        * bond_symbol
        """
        if not isinstance(other, BondDescriptor):
            return False

        for attr in self._eq_attr:
            if getattr(self, attr) != getattr(other, attr):
                return False

        return True

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        if self.index_ == 1 and not Config.show_bond_descriptor_one_index:
            index = ""
        else:
            index = self.index_
        return "[" + self.descriptor + str(index) + "]"

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        for attr in self._eq_attr:
            text += f"{attr}: {getattr(self, attr)}, "
        return text[:-1] + "}"

    @property
    def symbol(self) -> str:
        """ descriptor + index """
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
        """ the owner at the top of the parent tree """
        return self.parent.root

    @property
    def instances(self) -> list[BondDescriptorAtom]:
        """ list of BondDescriptorAtom that are based on this bonding descriptor """
        return self._instances

    @instances.setter
    def instances(self, instances: list[BondDescriptor]):
        self._instances = instances
        self._instances_up_to_date = False

    @property
    def bond_order(self) -> int | None:
        """ number of covalent bonds represented """
        return chemical_data.bond_mapping[self.bond_symbol]

    @property
    def implicit(self) -> bool:
        """ returns true is it is an implicit bonding descriptor """
        if self.descriptor == "":
            return True
        return False

    @property
    def aromatic(self) -> bool:
        """ limited accuracy """
        return self.bond_symbol == ":"


class BondDescriptorAtom:
    """
    this class represents a bonding descriptor and has the ability to bond
    it should be thought of as an atom with only on bonding site that can be any bond (single, double triple)
    """
    __slots__ = ["descriptor", "id_", "_bond", "__dict__", "parent"]

    def __init__(self,
                 id_: int,
                 bond_descriptor: BondDescriptor,
                 parent: BigSMILES | Branch | StochasticFragment | None = None,
                 **kwargs
                 ):
        """
        Parameters
        ----------
        id_: int
            id of bonding descriptor atom (id is limited to bonding descriptor atoms). Range: [1, inf]
        bond_descriptor: BondDescriptor
            bonding descriptor that this atom represents
        parent: BigSMILES | Branch | StochasticFragment | None
            the owner of the bonding descriptor atom
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
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

    def __eq__(self, other: BondDescriptorAtom):
        """
        equality is based on the following parameters:

        * id_
        * descriptor
        * bond.id_
        """
        if not isinstance(other, BondDescriptorAtom):
            return False

        if self.id_ != other.id_:
            return False

        if self.descriptor != other.descriptor:
            return False

        if self.bond is None:
            if other.bond is not None:
                return False
        else:
            if other.bond is None or self.id_ != other.id_:
                return False

        return True

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        return Config.add_color(self.descriptor.to_string(show_hydrogens, print_repr), 'Green', skip_color)

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        text += f"id_: {self.id_}"
        return text + "}"

    @property
    def bond(self) -> Bond:
        """ bond to bonding descriptor atom """
        return self._bond

    @bond.setter
    def bond(self, bond: Bond):
        if self.bond is not None:
            raise errors.ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond = bond

    @property
    def root(self) -> BigSMILES:
        """ the owner at the top of the parent tree """
        return self.parent.root

    @property
    def aromatic(self) -> bool:
        """ limited accuracy """
        return self.descriptor.aromatic


class Branch:
    """
    this class represents a branch
    """
    __slots__ = ["nodes", "id_", "parent", "__dict__"]

    def __init__(self, id_: int, parent: BigSMILES | StochasticFragment | Branch, **kwargs):
        """
        Parameters
        ----------
        id_: int
            id of bond (id is limited to bonds). Range: [1, inf]
        parent: BigSMILES | StochasticFragment | Branch
            owner of the branch
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
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

    def __eq__(self, other: Branch):
        """ equality is determining all nodes and rings are equal """
        if not isinstance(other, Branch):
            return False

        for node, other_node in zip(self.nodes, other.nodes):
            if node != other_node:
                return False

        return True

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        text = "".join(node.to_string(show_hydrogens, show_atom_index, print_repr, skip_color) for node in self.nodes)
        return "(" + text + ")"

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        text += f"id_: {self.id_}, "
        text += f"num_nodes: {len(self.nodes)}"
        return text + "}"

    @property
    def in_stochastic_object(self) -> bool:
        """ returns that if on the tree up to the root there is a stochastic object parent """
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        """ the owner at the top of the parent tree """
        return self.parent.root

    def _get_id(self, node_type) -> int:
        return self.root._get_id(node_type)  # noqa


class StochasticFragment:
    """
    this class represents a stochastic fragment.

    stochastic fragments are only found within stochastic objects
    """
    __slots__ = ["nodes", "id_", "parent", "bonding_descriptors", "rings", "__dict__"]

    def __init__(self, id_: int, parent: StochasticObject, **kwargs):
        """
        Parameters
        ----------
        id_: int
            id of bond (id is limited to bonds). Range: [1, inf]
        parent
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
        self.id_ = id_
        self.parent = parent
        self.nodes: list[Atom | Bond | BondDescriptorAtom | Branch | StochasticObject] = []
        self.rings: list[Bond] = []
        self.bonding_descriptors: list[BondDescriptor] = []

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __eq__(self, other: StochasticFragment):
        """ equality is determining all nodes and rings are equal """
        if not isinstance(other, StochasticFragment):
            return False

        for node, other_node in zip(self.nodes, other.nodes):
            if node != other_node:
                return False

        for ring, other_ring in zip(self.rings, other.rings):
            if ring != other_ring:
                return False

        return True

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        return "".join(node.to_string(show_hydrogens, show_atom_index, print_repr, skip_color) for node in self.nodes)

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        text += f"id_: {self.id_}, "
        text += f"num_nodes: {len(self.nodes)}"
        return text + "}"

    @property
    def in_stochastic_object(self) -> bool:
        """ returns that if on the tree up to the root there is a stochastic object parent """
        return self.parent.in_stochastic_object

    @property
    def root(self) -> BigSMILES:
        """ the owner at the top of the parent tree """
        return self.parent.root

    def _get_id(self, node_type) -> int:
        return self.root._get_id(node_type)  # noqa


class StochasticObject:
    """
    this class is represents a stochastic object
    """
    __slots__ = ["nodes", "bonding_descriptors", "id_", "parent", "bd_left", "bd_right",
                 "_bond_left", "_bond_right", "__dict__"]

    def __init__(self, id_: int, parent: BigSMILES | StochasticFragment | Branch, **kwargs):
        """
        Parameters
        ----------
        id_: int
            id of bond (id is limited to bonds). Range: [1, inf]
        parent: BigSMILES | StochasticFragment | Branch
            the owner of the atom
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
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

    def __eq__(self, other: StochasticObject):
        """ equality is determining all nodes and bonding descriptors are equal """
        if not isinstance(other, StochasticObject):
            return False

        for node, other_node in zip(self.nodes, other.nodes):
            if node != other_node:
                return False

        if self.bd_left != other.bd_left:
            return False

        if self.bd_right != other.bd_right:
            return False

        return True

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        text += f"id_: {self.id_}, "
        text += f"num_nodes: {len(self.nodes)}"
        return text + "}"

    @property
    def bond_left(self) -> Bond | None:
        """ bond leaving the left of the stochastic object """
        return self._bond_left

    @bond_left.setter
    def bond_left(self, bond: Bond):
        if self._bond_left is not None:
            raise errors.ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond_left = bond

    @property
    def bond_right(self) -> Bond | None:
        """ bond leaving the right of the stochastic object """
        return self._bond_right

    @bond_right.setter
    def bond_right(self, bond: Bond):
        if self._bond_right is not None:
            raise errors.ConstructorError(f"Trying to make a bond to {self} when it already has a bond.")
        self._bond_right = bond

    @property
    def aromatic(self) -> bool:
        """ limited accuracy """
        if self.bd_left.bond_order == 1.5:
            # not checking right bond descriptor, maybe we should? but during construction it is not defined
            return True
        return False

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        text = ""
        text += Config.add_color("{", "Red", skip_color)
        if self.bd_left is not None:
            text += self.bd_left.to_string(show_hydrogens, show_atom_index, print_repr, skip_color)
        text += ",".join(node.to_string(show_hydrogens, show_atom_index, print_repr, skip_color) for node in self.nodes)
        if self.bd_right is not None:
            text += self.bd_right.to_string(show_hydrogens, show_atom_index, print_repr, skip_color)
        text += Config.add_color("}", "Red", skip_color)
        if self.bond_right is not None and self.bond_right.ring_id is not None:
            if self.bond_right.symbol != ":":
                text += self.bond_right.symbol
            text += str(self.bond_right.ring_id)
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
        """ returns True as it is the stochastic object """
        return True

    @property
    def root(self) -> BigSMILES:
        """ the owner at the top of the parent tree """
        return self.parent.root

    @property
    def bonds(self) -> list[Bond]:
        """ bonds entering the left and right """
        bonds = []
        if self.bond_left is not None:
            bonds.append(self.bond_left)
        if self.bond_right is not None:
            bonds.append(self.bond_right)
        return bonds

    def _get_id(self, node_type) -> int:
        return self.root._get_id(node_type)  # noqa


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
    """
    this class is used to represent a BigSMILES.

    BigSMILES is a superset of SMILES; so any this class handles SMILES as well

    Examples
    --------
    bigsmiles("C1CCCCC1")

    bigsmiles("CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]")

    """
    __slots__ = ["nodes", "atoms", "bonds", "rings", "__dict__"]

    def __init__(self, text: str | None = None, **kwargs):
        """
        Parameters
        ----------
        text: str | None
            BigSMILES or SMILES string. String will be processed automatically. (pass 'None' if using 'constructor')
        kwargs:
            any additional keyword arguments are accepted and set as additional attributes
        """
        self.nodes: list[Atom | Bond | StochasticObject | Branch] = []
        self.atoms: list[Atom | BondDescriptorAtom] = []  # includes atoms in sub-objects
        self.bonds: list[Bond] = []  # includes bonds in sub-objects
        self.rings: list[Bond] = []  # does not include rings in sub-objects

        self._ids_ = {Atom: 0, Bond: 0, BondDescriptorAtom: 0, StochasticFragment: 0, StochasticObject: 0, Branch: 0}

        # parse input string
        if text:
            # import here to avoid circular imports
            from bigsmiles.constructors.constructor_tokens import parse_bigsmiles_str
            parse_bigsmiles_str(text, self)

        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string(print_repr=True, skip_color=True)

    def __getitem__(self, obj: int | slice) -> Atom | list[Atom]:
        """ get items by node index or slice """
        return self.nodes[obj]

    def __iter__(self):
        """ iterator over nodes """
        return iter(self.nodes)

    def __bool__(self):
        """ returns True if BigSMILES has nodes """
        if self.nodes:
            return True
        else:
            return False

    def __eq__(self, other: BigSMILES):
        """ equality is determining all nodes and rings are equal """
        if not isinstance(other, BigSMILES):
            return False

        for node, other_node in zip(self.nodes, other.nodes):
            if node != other_node:
                return False

        for ring, other_ring in zip(self.rings, other.rings):
            if ring != other_ring:
                return False

        return True

    def to_string(self,
                  show_hydrogens: bool = False,
                  show_atom_index: bool = True,
                  print_repr: bool = False,
                  skip_color: bool = False
                  ) -> str:
        return "".join(node.to_string(show_hydrogens, show_atom_index, print_repr, skip_color) for node in self.nodes)

    @property
    def details(self) -> str:
        """ long string representation """
        text = str(self) + "  {"
        text += f"num_nodes: {len(self.nodes)}, "
        text += f"num_atoms: {len(self.atoms)}, "
        text += f"num_bonds: {len(self.bonds)}"
        return text + "}"

    @property
    def in_stochastic_object(self) -> bool:
        """ returns False as BigSMILES is always a root node and can't be contained in anything """
        return False

    @property
    def contains_stochastic_object(self) -> bool:
        """ returns True if a stochastic object is present in the BigSMILES (anywhere) """
        return contains_stochastic_object(self.nodes)

    @property
    def has_disconnect(self) -> bool:
        """ returns True if '.' (disconnect symbol) is present in the BigSMILES (anywhere) """
        for bond in self.bonds:
            if bond.symbol == ".":
                return True
        return False

    @property
    def root(self) -> BigSMILES:
        """ root will return itself as BigSMILES object is always the root node """
        return self

    def _get_id(self, node_type) -> int:
        self._ids_[node_type] += 1
        return self._ids_[node_type]


# types  (for python >3.10)
# has_node_attr = BigSMILES | Branch | StochasticObject | StochasticFragment
# has_ring_attr = BigSMILES | StochasticFragment
# has_parent_attr = Branch | StochasticObject | StochasticFragment | Bond | BondDescriptorAtom | Atom
# has_root_attr = BigSMILES | has_parent_attr

# types (for python <=3.9)
import typing

bigsmiles_obj = [Atom, Bond, Branch, BondDescriptorAtom, StochasticObject, StochasticFragment, BigSMILES]
has_node_attr = typing.Union[BigSMILES, Branch, StochasticObject, StochasticFragment]
has_ring_attr = typing.Union[BigSMILES, StochasticObject]
has_parent_attr = typing.Union[Branch, StochasticObject, StochasticFragment, Bond, BondDescriptorAtom, Atom]
has_root_attr = typing.Union[BigSMILES, has_parent_attr]
