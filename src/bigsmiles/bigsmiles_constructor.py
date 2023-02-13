import enum
from functools import wraps

from bigsmiles.errors import BigSMILESError
from bigsmiles.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, StochasticObject, \
    BigSMILES, BondDescriptorAtom
from bigsmiles.validation import run_validation


class ConstructorStates(enum.Enum):
    start = 0
    atom = 1
    bond_descriptor = 2
    branch_start = 3
    branch = 4
    ring = 5
    stochastic_fragment = 6
    stochastic_object_start = 7
    stochastic_object_end = 8
    stochastic_object = 9


def check_atom_for_making_bond(bond: Bond, atom: Atom):
    """
    Check to see if atom can accept bond.

    Parameters
    ----------
    bond: Bond
        Bond to add to Atom
    atom: Atom
        Atom to receive Bond

    """
    if atom.bonds_available is None:
        pass
    elif atom.bonds_available == 0:
        if atom.valance == atom.valance_possible[-1] or atom.valance_possible[-1] < atom.valance + bond.type_.value:
            raise BigSMILESError("Too many bonds trying to be made.", repr(atom))

        for val in bond.atom1.valance_possible:
            if val > atom.valance:
                atom.valance = val
                break

    atom.bonds.append(bond)


def add_bond_to_connected_objects(bond: Bond):
    """ Adds bonds to Atoms, and BondDescriptorAtom. """
    for obj in bond:
        if isinstance(obj, Atom):
            check_atom_for_making_bond(bond, obj)
        elif isinstance(obj, BondDescriptorAtom):
            obj.bond = bond
        elif isinstance(obj, StochasticObject):  # left bond add with
            obj.bond_right = bond


def add_bond_descriptor_to_stochastic_fragment(stoch_frag: StochasticFragment, loop_obj: Branch = None):
    if loop_obj is None:
        loop_obj = stoch_frag

    for obj in loop_obj.nodes:
        if isinstance(obj, BondDescriptorAtom) and obj.descriptor not in stoch_frag.bonding_descriptors:
            stoch_frag.bonding_descriptors.append(obj.descriptor)
        if isinstance(obj, Branch):
            add_bond_descriptor_to_stochastic_fragment(stoch_frag, obj)  # recursive


def in_stochastic_object(func):
    """ Decorator to ensure function call only occurs within stochastic object."""
    @wraps(func)
    def _in_stochastic_object(*args, **kwargs):
        if not args[0].stack[-1].in_stochastic_object:
            raise BigSMILESError("Must in stochastic object.")
        return func(*args, **kwargs)

    return _in_stochastic_object


class BigSMILESConstructor:

    def __init__(self, obj: BigSMILES):
        self.bigsmiles = obj

        # setup id counters
        self._atom_counter = 0
        self._bond_counter = 0
        self._bond_descriptor_atom = 0
        self._branch_counter = 0
        self._stochastic_fragment = 0
        self._stochastic_object = 0
        self._ring_counter = 1  # start ring counter at 1

        self.state = ConstructorStates.start
        self.stack: list[BigSMILES | StochasticObject | StochasticFragment | Branch] = [self.bigsmiles]

    def _get_atom_id(self) -> int:
        self._atom_counter += 1
        return self._atom_counter - 1

    def _get_bond_id(self) -> int:
        self._bond_counter += 1
        return self._bond_counter - 1

    def _get_ring_id(self) -> int:
        self._ring_counter += 1
        return self._ring_counter - 1

    def _get_bond_descriptor_atom_id(self) -> int:
        self._bond_descriptor_atom += 1
        return self._bond_descriptor_atom - 1

    def _get_branch_id(self) -> int:
        self._branch_counter += 1
        return self._branch_counter - 1

    def _get_stochastic_fragment_id(self) -> int:
        self._stochastic_object += 1
        return self._stochastic_object - 1

    def _get_stochastic_object_id(self) -> int:
        self._stochastic_object += 1
        return self._stochastic_object - 1

    def add_ring(self, ring_id: int) -> Bond:
        # check if ring_id already exists
        ring_parent = self._get_ring_parent()
        for ring in ring_parent.rings:
            if ring.ring_id == ring_id:
                if ring.atom2 is not None:
                    raise BigSMILESError(f"Ring already formed for ring id {ring_id}.")
                ring.atom2 = self._get_prior(self.stack[-1], (Atom,))
                add_bond_to_connected_objects(ring)
                return ring

        # Make new ring_id
        bond = Bond("", self._get_prior(self.stack[-1], (Atom,)), None, self._get_bond_id(), ring_id)
        self.bigsmiles.bonds.append(bond)
        ring_parent.rings.append(bond)
        # add_bond_to_connected_objects(bond)

        return bond

    def add_ring_from_atoms(self, atom1: Atom, atom2: Atom):
        ring_parent = self._get_ring_parent()

        # Make new ring_id
        bond = Bond("", atom1, atom2, self._get_bond_id(), self._get_ring_id())
        self.bigsmiles.bonds.append(bond)
        ring_parent.rings.append(bond)
        add_bond_to_connected_objects(bond)

    def _get_ring_parent(self) -> BigSMILES | StochasticFragment:
        for node in reversed(self.stack):
            if hasattr(node, "rings"):
                return node

    def add_atom(self, symbol: str) -> Atom:
        atom = Atom(symbol, self._get_atom_id())
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = ConstructorStates.atom
        return atom

    def add_bond(self, bond_symbol: str, atom1: Atom | BondDescriptorAtom, atom2: Atom | BondDescriptorAtom | None) -> Bond:
        bond = Bond(bond_symbol, atom1, atom2, self._get_bond_id())
        self.stack[-1].nodes.append(bond)
        self.bigsmiles.bonds.append(bond)
        add_bond_to_connected_objects(bond)
        return bond

    def _get_bonding_descriptor_atom(self, bd_symbol: str) -> BondDescriptorAtom:
        bd = self._get_bonding_descriptor(bd_symbol)
        return BondDescriptorAtom(bd, self._get_bond_descriptor_atom_id())

    def _get_bonding_descriptor(self, bd_symbol: str) -> BondDescriptor:
        stoch_obj = self._get_current_stochastic_object()

        # check if bd in there already
        for bd in stoch_obj.bonding_descriptors:
            if bd_symbol == bd._text:
                return bd

        # create new bonding descriptor
        new_bd = BondDescriptor(stoch_obj, bd_symbol)
        stoch_obj.bonding_descriptors.append(new_bd)
        return new_bd

    def _get_current_stochastic_object(self) -> StochasticObject:
        for obj in self.stack:
            if isinstance(obj, StochasticObject):
                return obj

        raise BigSMILESError("Coding error.")

    @in_stochastic_object
    def add_bonding_descriptor(self, bd_symbol: str) -> BondDescriptorAtom:
        """ [<], [>], [$], [$1], [>2], ... """
        bd_atom = self._get_bonding_descriptor_atom(bd_symbol)
        self.stack[-1].nodes.append(bd_atom)

        self.state = ConstructorStates.bond_descriptor
        return bd_atom

    def add_bond_atom_pair(self, bond_symbol: str, atom_symbol: str) -> tuple[Bond, Atom]:
        atom = Atom(atom_symbol, self._get_atom_id())
        prior_atom = self._get_prior(self.stack[-1], (Atom, BondDescriptorAtom, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, atom)
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = ConstructorStates.atom
        return bond, atom

    @in_stochastic_object
    def add_bond_bonding_descriptor_pair(self, bond_symbol: str, bd_symbol: str) -> tuple[Bond, BondDescriptorAtom]:
        bd_atom = self._get_bonding_descriptor_atom(bd_symbol)
        prior_atom = self._get_prior(self.stack[-1], (Atom, BondDescriptorAtom, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, bd_atom)
        self.stack[-1].nodes.append(bd_atom)

        self.state = ConstructorStates.bond_descriptor
        return bond, bd_atom

    def open_branch(self):
        if isinstance(self.stack[-1], Branch) and len(self.stack[-1].nodes) == 0:
            raise BigSMILESError("BigSMILES string or branch can't start with branch")

        branch = Branch(self.stack[-1], self._get_branch_id())
        self.stack[-1].nodes.append(branch)
        self.stack.append(branch)
        self.state = ConstructorStates.branch_start

    def close_branch(self):
        if not isinstance(self.stack[-1], Branch):
            raise BigSMILESError("Error closing branch. Possible issues: "
                                 "\n\tClosing a branch with another intermediate node started."
                                 "\n\tNo starting branch symbol.")

        self.stack.pop()

    def open_stochastic_object(self, bd_symbol: str) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        new_bd = BondDescriptor(stoch_obj, bd_symbol)
        stoch_obj.bonding_descriptors.append(new_bd)
        stoch_obj.end_group_left = BondDescriptorAtom(new_bd, self._get_bond_descriptor_atom_id())
        self.stack[-1].nodes.append(stoch_obj)
        self.stack.append(stoch_obj)

        # open stochastic_fragment
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = ConstructorStates.stochastic_fragment
        return stoch_frag

    def open_stochastic_object_with_bond(self, bond_symbol: str, bd_symbol: str) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        new_bd = BondDescriptor(stoch_obj, bd_symbol)
        stoch_obj.bonding_descriptors.append(new_bd)
        stoch_obj.end_group_left = BondDescriptorAtom(new_bd, self._get_bond_descriptor_atom_id())

        prior_atom = self._get_prior(self.stack[-1], (Atom, BondDescriptor, StochasticObject))
        bond = Bond(bond_symbol, prior_atom, stoch_obj, self._get_bond_id())
        stoch_obj.bond_left = bond
        self.stack[-1].nodes.append(bond)

        self.stack[-1].nodes.append(stoch_obj)
        self.stack.append(stoch_obj)

        return self.open_stochastic_fragment()

    def close_stochastic_object(self, bd_symbol: str):
        if not isinstance(self.stack[-1], StochasticObject):
            raise BigSMILESError("Error closing StochasticObject. Possible issues: "
                                 "\n\tClosing a StochasticObject with another intermediate node started."
                                 "\n\tNo starting StochasticObject symbol.")

        self.stack[-1].end_group_right = self._get_bonding_descriptor_atom(bd_symbol)
        self.stack.pop()

    @in_stochastic_object
    def open_stochastic_fragment(self) -> StochasticFragment:
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = ConstructorStates.stochastic_fragment
        return stoch_frag

    @in_stochastic_object
    def close_open_stochastic_fragment(self) -> StochasticFragment:
        if not isinstance(self.stack[-1], StochasticFragment):
            raise BigSMILESError("Stochastic seperator can only follow a stochastic fragments.")

        self.close_stochastic_fragment()
        return self.open_stochastic_fragment()

    @in_stochastic_object
    def close_stochastic_fragment(self):
        if not isinstance(self.stack[-1], StochasticFragment):
            raise BigSMILESError("Error StochasticFragment branch. Possible issues: "
                                 "\n\tClosing a branch with another intermediate node started.")

        add_bond_descriptor_to_stochastic_fragment(self.stack[-1])

        # check for at-least one bonding descriptor
        if not self.stack[-1].bonding_descriptors:
            raise BigSMILESError(f"No bonding descriptor in the stochastic fragment. ({repr(self.stack[-1])})")

        self.stack.pop()

    def _get_prior(self, obj: BigSMILES | StochasticObject | StochasticFragment | Branch, types_: tuple,
                   flag: bool = False) -> Atom:
        if obj.nodes:
            for node in reversed(obj.nodes):
                if isinstance(node, types_):
                    return node

            raise BigSMILESError("Bug in the code")

        elif not flag:
            return self._get_prior(self.stack[-2], types_, flag=True)

        raise BigSMILESError("Bond attempted to be made to that has nothing to bond back to.")

    def final_validation(self):
        """ Run various validations. """
        if len(self.stack) != 1:
            raise BigSMILESError(f"{type(self.stack[-1])} is missing closing symbol.")

        run_validation(self.bigsmiles)
