import enum
from functools import wraps

from bigsmiles.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, StochasticObject, \
    BigSMILES


class BigSMILESError(Exception):
    def __init__(self, text: str, optional: str = None):
        if optional is not None:
            text += "(" + optional + ")"
        super().__init__(text)


class States(enum.Enum):
    start = enum.auto()
    atom = enum.auto()
    bond_descriptor = enum.auto()
    branch_start = enum.auto()
    branch = enum.auto()
    stochastic_fragment = enum.auto()
    stochastic_object_start = enum.auto()
    stochastic_object_end = enum.auto()
    stochastic_object = enum.auto()


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
    """ Adds bonds to Atoms, and BondDescriptor. """
    for obj in bond:
        if isinstance(obj, Atom):
            check_atom_for_making_bond(bond, obj)
        elif isinstance(obj, BondDescriptor):
            obj.bond = bond


def in_stochastic_object(func):
    """ Decorator to ensure function call only occurs within stochastic object."""
    @wraps(func)
    def _in_stochastic_object(*args, **kwargs):
        if not args[0].stack[-1].in_stochastic_object:
            raise BigSMILESError("Must be in stochastic object.")
        return func(*args, **kwargs)

    return _in_stochastic_object


class BigSMILESConstructor:

    def __init__(self, obj: BigSMILES):
        self.bigsmiles = obj

        # setup id counters
        self._atom_counter = 0
        self._bond_counter = 0
        self._bond_descriptor = 0
        self._branch_counter = 0
        self._stochastic_fragment = 0
        self._stochastic_object = 0

        self.state = States.start
        self.stack: list[BigSMILES | StochasticObject | StochasticFragment | Branch] = [self.bigsmiles]

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

    def add_ring(self, ring_id: int) -> Bond:
        # check if ring already exists
        for ring in self.bigsmiles.rings:
            if ring.ring == ring_id:
                if ring.atom2 is not None:
                    raise BigSMILESError("Ring already formed")
                ring.atom2 = self._get_prior(self.stack[-1], (Atom,))
                return ring

        # Make new ring
        bond = Bond("", self._get_prior(self.stack[-1], (Atom,)), None, self._get_bond_id(), ring_id)
        self.bigsmiles.nodes.append(bond)
        self.bigsmiles.bonds.append(bond)
        add_bond_to_connected_objects(bond)

        return bond

    def add_atom(self, symbol: str) -> Atom:
        atom = Atom(symbol, self._get_atom_id())
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = States.atom
        return atom

    def add_bond(self, bond_symbol: str, atom1: Atom | BondDescriptor, atom2: Atom | BondDescriptor | None) -> Bond:
        bond = Bond(bond_symbol, atom1, atom2, self._get_bond_id())
        self.stack[-1].nodes.append(bond)
        self.bigsmiles.bonds.append(bond)
        add_bond_to_connected_objects(bond)
        return bond

    @in_stochastic_object
    def add_bonding_descriptor(self, bd_symbol: str) -> BondDescriptor:
        """ [<], [>], [$], [$1], [>2], ... """
        bd = BondDescriptor(bd_symbol, self._get_bond_descriptor_id())
        self.stack[-1].nodes.append(bd)

        self.state = States.bond_descriptor
        return bd

    def add_bond_atom_pair(self, bond_symbol: str, atom_symbol: str) -> tuple[Bond, Atom]:
        atom = Atom(atom_symbol, self._get_atom_id())
        prior_atom = self._get_prior(self.stack[-1], (Atom, BondDescriptor, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, atom)
        self.stack[-1].nodes.append(atom)
        self.bigsmiles.atoms.append(atom)

        self.state = States.atom
        return bond, atom

    @in_stochastic_object
    def add_bond_bonding_descriptor_pair(self, bond_symbol: str, bd_symbol: str) -> tuple[Bond, BondDescriptor]:
        bd = BondDescriptor(bd_symbol, self._get_bond_descriptor_id())
        prior_atom = self._get_prior(self.stack[-1], (Atom, BondDescriptor, StochasticObject))
        bond = self.add_bond(bond_symbol, prior_atom, bd)
        self.stack[-1].nodes.append(bd)

        self.state = States.bond_descriptor
        return bond, bd

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
        self.stack.pop()

    def open_stochastic_object(self, bd_symbol: str) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        bd = BondDescriptor(bd_symbol, self._get_bond_descriptor_id())
        stoch_obj.end_group_left = bd
        self.stack[-1].nodes.append(stoch_obj)
        self.stack.append(stoch_obj)

        # open stochastic_fragment
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = States.stochastic_fragment
        return stoch_frag

    def open_stochastic_object_with_bond(self, bond_symbol: str, bd_symbol: str) -> StochasticFragment:
        stoch_obj = StochasticObject(self.stack[-1], self._get_stochastic_object_id())
        bd = BondDescriptor(bd_symbol, self._get_bond_descriptor_id())
        stoch_obj.end_group_left = bd

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

        self.stack[-1].end_group_right = BondDescriptor(bd_symbol, self._get_bond_descriptor_id())
        self.stack.pop()

    @in_stochastic_object
    def open_stochastic_fragment(self) -> StochasticFragment:
        stoch_frag = StochasticFragment(self.stack[-1], self._get_stochastic_fragment_id())
        self.stack[-1].nodes.append(stoch_frag)
        self.stack.append(stoch_frag)
        self.state = States.stochastic_fragment
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

        self.stack.pop()

    def _get_prior(self, obj: BigSMILES | StochasticObject | StochasticFragment | Branch, types_: tuple,
                   flag: bool = False) -> Atom | None:
        if obj.nodes:
            if isinstance(obj.nodes[-1], types_):
                return obj.nodes[-1]
            if isinstance(obj.nodes[-1], Branch):
                if isinstance(obj.nodes[-2], types_):
                    return obj.nodes[-2]
                if isinstance(obj.nodes[-2], Branch):
                    if isinstance(obj.nodes[-3], types_):
                        return obj.nodes[-3]
        elif hasattr(obj, "parent") and not flag:
            return self._get_prior(self.stack[-2], types_, flag=True)

        return None
