import enum

from bigsmiles.errors import BigSMILESError
from bigsmiles.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, StochasticObject, \
    BigSMILES, BondDescriptorAtom
from bigsmiles.validation import post_validation
from bigsmiles.syntax_fixes import run_syntax_fixes


class States(enum.Enum):
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


def add_bond_to_connected_objects(bond: Bond):
    """ Adds bonds to Atoms, and BondDescriptorAtom. """
    for obj in bond:
        if isinstance(obj, Atom):
            if obj.bonds_available < bond.bond_order:
                if obj._increase_valance(bond.bond_order):
                    continue
                raise BigSMILESError("Too many bonds trying to be made.", str(obj))
            obj.bonds.append(bond)
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


def get_common_bond(atom1: Atom, atom2: Atom) -> Bond | None:
    # Check if bond already between two atoms
    for bond in atom1.bonds:
        if bond in atom2.bonds:
            return bond

    return None  # No common bond found


def remove_partial_ring(parent: BigSMILES | Branch | StochasticFragment, ring_id: int):
    for i, ring in enumerate(parent.rings):
        if ring.ring_id == ring_id:
            break
    else:
        return

    bond = parent.rings.pop(i)
    parent.bonds.remove(bond)
    if hasattr(parent, "parent"):
        parent.parent.bonds.remove(bond)


class BigSMILESConstructor:

    def __init__(self, bigsmiles: BigSMILES = None, syntax_fix: bool = True):
        self.bigsmiles = bigsmiles if bigsmiles is not None else BigSMILES()

        self.syntax_fix = syntax_fix  # only with context manager

        # setup id counters
        self._atom_counter = 0
        self._bond_counter = 0
        self._ring_counter = 1  # start ring counter at 1
        self._bond_descriptor_atom = 0
        self._branch_counter = 0
        self._stochastic_fragment_counter = 0
        self._stochastic_object_counter = 0

        self.state = States.start
        self.stack: list[BigSMILES | StochasticObject | StochasticFragment | Branch] = [self.bigsmiles]

    def exit(self):
        if len(self.stack) != 1:
            raise BigSMILESError(f"{type(self.stack[-1])} is missing closing symbol.")
        if self.syntax_fix:
            run_syntax_fixes(self.bigsmiles)

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
        self._stochastic_object_counter += 1
        return self._stochastic_object_counter - 1

    def _get_stochastic_object_id(self) -> int:
        self._stochastic_object_counter += 1
        return self._stochastic_object_counter - 1

    def add_ring(self, ring_id: int, bond_symbol: str = "") -> Bond:
        # check if ring_id already exists
        ring_parent = self._get_ring_parent()
        for ring in ring_parent.rings:
            if ring.ring_id == ring_id:
                if ring.atom2 is not None:
                    raise BigSMILESError(f"Ring already formed for ring id {ring_id}.")
                atom2 = self.get_prior(self.stack[-1], (Atom,))

                bond = get_common_bond(ring.atom1, atom2)
                if bond:  # Check if ring already between two atoms and just increase bond order
                    bond.bond_order += 1
                    remove_partial_ring(ring_parent, ring_id)
                else:
                    ring.atom2 = atom2
                    add_bond_to_connected_objects(ring)
                return ring

        # Make new ring_id
        bond = Bond(bond_symbol, self.get_prior(self.stack[-1], (Atom,)), None, self._get_bond_id(), ring_id)
        self.bigsmiles.bonds.append(bond)
        ring_parent.rings.append(bond)

        return bond

    def add_ring_from_atoms(self, atom1: Atom, atom2: Atom, bond_symbol: str = ""):
        ring_parent = self._get_ring_parent()

        bond = get_common_bond(atom1, atom2)
        if bond:  # Check if ring already between two atoms and just increase bond order
            bond.bond_order += 1
        else:
            # Make new ring_id
            bond = Bond(bond_symbol, atom1, atom2, self._get_bond_id(), self._get_ring_id())
            self.bigsmiles.bonds.append(bond)
            ring_parent.rings.append(bond)
            add_bond_to_connected_objects(bond)

    def _get_ring_parent(self) -> BigSMILES | StochasticFragment:
        for node in reversed(self.stack):
            if hasattr(node, "rings"):
                return node

    def add_atom(self,
                 element: str,
                 isotope: int | None = None,
                 stereo: str = '',
                 hcount: int = 0,
                 charge: int = 0,
                 valance: int = None,
                 **kwargs
                 ) -> Atom:
        atom = Atom(self._get_atom_id(), element, isotope, stereo, hcount, charge, valance, **kwargs)
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = States.atom
        return atom

    def add_bond(self,
                 bond_symbol: str,
                 atom1: Atom | BondDescriptorAtom,
                 atom2: Atom | BondDescriptorAtom | StochasticObject | None
                 ) -> Bond:
        bond = Bond(bond_symbol, atom1, atom2, self._get_bond_id())
        self.stack[-1].nodes.append(bond)
        self.bigsmiles.bonds.append(bond)
        add_bond_to_connected_objects(bond)
        return bond

    def _get_bonding_descriptor_atom(self, descriptor: str, index_: int) -> BondDescriptorAtom:
        bd = self._get_bonding_descriptor(descriptor, index_)
        return BondDescriptorAtom(bd, self._get_bond_descriptor_atom_id())

    def _get_bonding_descriptor(self, descriptor: str, index_: int) -> BondDescriptor:
        stoch_obj = self._get_current_stochastic_object()

        # check if bd in there already
        for bd in stoch_obj.bonding_descriptors:
            if descriptor == bd.descriptor and index_ == bd.index_:
                return bd

        # create new bonding descriptor
        new_bd = BondDescriptor(stoch_obj, descriptor, index_)
        stoch_obj.bonding_descriptors.append(new_bd)
        return new_bd

    def _get_current_stochastic_object(self) -> StochasticObject:
        for obj in self.stack:
            if isinstance(obj, StochasticObject):
                return obj

        raise BigSMILESError("Coding error.")

    def add_bonding_descriptor(self, descriptor: str, index_: int) -> BondDescriptorAtom:
        """ [<], [>], [$], [$1], [>2], ... """
        bd_atom = self._get_bonding_descriptor_atom(descriptor, index_)
        self.stack[-1].nodes.append(bd_atom)

        self.state = States.bond_descriptor
        return bd_atom

    def add_bond_atom_pair(self,
                           bond_symbol: str,
                           element: str,
                           isotope: int | None = None,
                           stereo: str = '',
                           hcount: int = 0,
                           charge: int = 0,
                           valance: int = None,
                           **kwargs
                           ) -> tuple[Bond, Atom]:
        atom = Atom(self._get_atom_id(), element, isotope, stereo, hcount, charge, valance, **kwargs)
        prior_atom = self.get_prior(self.stack[-1], (Atom, BondDescriptorAtom, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, atom)
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = States.atom
        return bond, atom

    def add_bond_bonding_descriptor_pair(self, bond_symbol: str, descriptor: str, index_: int) \
            -> tuple[Bond, BondDescriptorAtom]:
        bd_atom = self._get_bonding_descriptor_atom(descriptor, index_)
        prior_atom = self.get_prior(self.stack[-1], (Atom, BondDescriptorAtom, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, bd_atom)
        self.stack[-1].nodes.append(bd_atom)

        self.state = States.bond_descriptor
        return bond, bd_atom

    def open_branch(self):
        if isinstance(self.stack[-1], Branch) and len(self.stack[-1].nodes) == 0:
            raise BigSMILESError("BigSMILES string or branch can't start with branch")

        branch = Branch(self.stack[-1], self._get_branch_id())
        self.stack[-1].nodes.append(branch)
        self.stack.append(branch)
        self.state = States.branch_start

    def close_branch(self):
        if not isinstance(self.stack[-1], Branch):
            raise BigSMILESError("Error closing branch. Possible issues: "
                                 "\n\tClosing a branch with another intermediate node started."
                                 "\n\tNo starting branch symbol.")

        branch = self.stack.pop()

        # check if branch is empty; if so remove it
        if len(branch.nodes) == 0:
            self.stack[-1].nodes.pop()
            self._branch_counter -= 1

    def open_stochastic_object(self, descriptor: str, index_: int) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        new_bd = BondDescriptor(stoch_obj, descriptor, index_)
        stoch_obj.bonding_descriptors.append(new_bd)
        stoch_obj.end_group_left = BondDescriptorAtom(new_bd, self._get_bond_descriptor_atom_id())
        self.stack[-1].nodes.append(stoch_obj)
        self.stack.append(stoch_obj)

        # open stochastic_fragment
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = States.stochastic_fragment
        return stoch_frag

    def open_stochastic_object_with_bond(self, bond_symbol: str, descriptor: str, index_: int) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        new_bd = BondDescriptor(stoch_obj, descriptor, index_)
        stoch_obj.bonding_descriptors.append(new_bd)
        stoch_obj.end_group_left = BondDescriptorAtom(new_bd, self._get_bond_descriptor_atom_id())

        prior_atom = self.get_prior(self.stack[-1], (Atom, BondDescriptor, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, stoch_obj)
        stoch_obj.bond_left = bond

        self.stack[-1].nodes.append(stoch_obj)
        self.stack.append(stoch_obj)

        return self.open_stochastic_fragment()

    def close_stochastic_object(self, descriptor: str, index_: int):
        if not isinstance(self.stack[-1], StochasticObject):
            raise BigSMILESError("Error closing StochasticObject. Possible issues: "
                                 "\n\tClosing a StochasticObject with another intermediate node started."
                                 "\n\tNo starting StochasticObject symbol.")

        self.stack[-1].end_group_right = self._get_bonding_descriptor_atom(descriptor, index_)
        self.stack.pop()

    def open_stochastic_fragment(self) -> StochasticFragment:
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = States.stochastic_fragment
        return stoch_frag

    def close_open_stochastic_fragment(self) -> StochasticFragment:
        if not isinstance(self.stack[-1], StochasticFragment):
            raise BigSMILESError("Stochastic seperator can only follow a stochastic fragments.")

        self.close_stochastic_fragment()
        return self.open_stochastic_fragment()

    def close_stochastic_fragment(self):
        if not isinstance(self.stack[-1], StochasticFragment):
            raise BigSMILESError("Error StochasticFragment branch. Possible issues: "
                                 "\n\tClosing a branch with another intermediate node started.")

        add_bond_descriptor_to_stochastic_fragment(self.stack[-1])

        # check for at-least one bonding descriptor
        if not self.stack[-1].bonding_descriptors:
            raise BigSMILESError(f"No bonding descriptor in the stochastic fragment. ({repr(self.stack[-1])})")

        self.stack.pop()

    def add_fragment(self, bond_symbol: str | None, bigsmiles: BigSMILES):
        if bond_symbol is not None:
            if not isinstance(bigsmiles.nodes[0], Atom | StochasticObject):
                raise BigSMILESError("First node must be an 'Atom' to added fragment.")

            atom = bigsmiles.nodes[0]
            prior_atom = self.get_prior(self.stack[-1], (Atom, BondDescriptorAtom, StochasticObject))
            self.add_bond(bond_symbol, prior_atom, atom)
            # self.stack[-1].nodes.append(atom)
            # self.bigsmiles.atoms.append(atom)

        # append bigsmiles
        # re-direct 'parent' to new bigsmiles object
        for node in bigsmiles.nodes:
            if hasattr(node, 'parent'):
                node.parent = self.bigsmiles
        # re-numbering
        self.re_number_rings(bigsmiles, set())
        self.re_number_nodes(bigsmiles)

        # append fragment
        self.bigsmiles.nodes += bigsmiles.nodes
        self.bigsmiles.atoms += bigsmiles.atoms
        self.bigsmiles.bonds += bigsmiles.bonds
        self.bigsmiles.rings += bigsmiles.rings

        # self.stack[-1].nodes.append(atom)
        self.state = States.atom

    def re_number_rings(self, obj, seen: set[Bond]):
        """ Recursive renumbering of rings. """
        for node in obj.nodes:
            if isinstance(node, Bond) and node.ring_id is not None and node not in seen:
                node.ring_id = self._get_ring_id()
                seen.add(node)
            if hasattr(node, 'nodes'):
                self.re_number_rings(node, seen)

    def re_number_nodes(self, obj):
        """ Recursive renumbering 'id_'. """
        for node in obj.nodes:
            if isinstance(node, Atom):
                node.id_ = self._get_atom_id()
            elif isinstance(node, Bond):
                node.id_ = self._get_bond_id()
            elif isinstance(node, BondDescriptorAtom):
                node.id_ = self._get_bond_descriptor_atom_id()
            elif isinstance(node, Branch):
                node.id_ = self._get_branch_id()
                self.re_number_nodes(node)
            elif isinstance(node, StochasticFragment):
                node.id_ = self._get_stochastic_fragment_id()
                self.re_number_nodes(node)
            elif isinstance(node, StochasticObject):
                node.id_ = self._get_stochastic_object_id()
                self.re_number_nodes(node)

    def get_prior(self, obj: BigSMILES | StochasticObject | StochasticFragment | Branch, types_: tuple,
                   flag: bool = False) -> Atom:
        if obj.nodes:
            for node in reversed(obj.nodes):  # no atom can ever have more than 10 bonds
                if isinstance(node, types_):
                    return node

            raise BigSMILESError("Bug in the code")

        elif not flag:
            return self.get_prior(self.stack[-2], types_, flag=True)

        raise BigSMILESError("Bond attempted to be made to that has nothing to bond back to.")

    def run_validation(self):
        """ Run various validations. """
        post_validation(self.bigsmiles)
