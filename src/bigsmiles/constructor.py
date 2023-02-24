from bigsmiles.errors import BigSMILESError
from bigsmiles.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, StochasticObject, \
    BigSMILES, BondDescriptorAtom
from bigsmiles.validation import post_validation
from bigsmiles.syntax_fixes import run_syntax_fixes

ParentType = BigSMILES | Branch | StochasticFragment | StochasticObject


## Utility functions for construction ##
######################################################################################################################
def add_bond_to_connected_objects(bond: Bond):
    """ Adds bonds to Atoms, and BondDescriptorAtom. """
    for obj in bond:
        if isinstance(obj, Atom):

            # Validation for available valence
            if obj.bonds_available < bond.bond_order:
                if obj._increase_valence(bond.bond_order):
                    continue
                raise BigSMILESError("Too many bonds trying to be made.", str(obj))

            obj.bonds.append(bond)
        elif isinstance(obj, BondDescriptorAtom):
            obj.bond = bond
        elif isinstance(obj, StochasticObject):  # left bond add with
            obj.bond_right = bond


def add_bond_descriptor_to_stochastic_fragment(stoch_frag: StochasticFragment, loop_obj: Branch = None):
    """ Add bonding descriptor to stoch_frag.bonding_descriptors. Recursive """
    if loop_obj is None:
        loop_obj = stoch_frag

    for obj in loop_obj.nodes:
        if isinstance(obj, BondDescriptorAtom) and obj.descriptor not in stoch_frag.bonding_descriptors:
            stoch_frag.bonding_descriptors.append(obj.descriptor)
        if isinstance(obj, Branch):
            add_bond_descriptor_to_stochastic_fragment(stoch_frag, obj)  # recursive


def get_common_bond(atom1: Atom, atom2: Atom) -> Bond | None:
    """ Checks if the two atoms have a bond in common. """
    # Check if bond already between two atoms
    for bond in atom1.bonds:
        if bond in atom2.bonds:
            return bond

    return None  # No common bond found


def remove_partial_ring(parent: BigSMILES | Branch | StochasticFragment, ring_id: int):
    """ Given Ring id, remove ring. """
    for i, ring in enumerate(parent.rings):
        if ring.ring_id == ring_id:
            bond = parent.rings.pop(i)
            parent.bonds.remove(bond)
            if hasattr(parent, "parent"):
                parent.parent.bonds.remove(bond)
            return


def insert_obj_into_list(list_: list, obj, object_to_insert):
    """ Insert 'object_to_insert' in a list right behind. """
    index_ = list_.index(obj)
    for i, item in enumerate(list_[index_:]):
        if isinstance(item, Branch):  # put it at the end of branches
            continue

        list_.insert(index_ + 1 + i, object_to_insert)
        break


def get_prior(obj: BigSMILES | StochasticObject | StochasticFragment | Branch, types_: tuple,
              flag: bool = False) -> Atom:
    if obj.nodes:
        for node in reversed(obj.nodes):  # no atom can ever have more than 10 bonds
            if type(node) in types_:
                return node

        raise BigSMILESError("Trying to make a bond to a prior atom that isn't there.??")

    elif not flag:
        return get_prior(obj.parent, types_, flag=True)

    raise BigSMILESError("Bond attempted to be made to that has nothing to bond back to.")


def get_ring_parent(parent) -> BigSMILES | StochasticFragment:
    """ Recursively find a parent that has the 'ring' attribute. only BigSMILES and Stochastic Fragments do. """
    if hasattr(parent, "rings"):
        return parent

    return get_ring_parent(parent.parent)


def run_validation(parent):
    """ Run various validations. """
    post_validation(parent)


def exit_(parent):
    if not isinstance(parent, BigSMILES):
        raise BigSMILESError(f"{type(parent)} is missing closing symbol.")
    run_syntax_fixes(parent)


## Functions for constructing a bigsmiles sequentially or step by step ##
#######################################################################################################################
def add_atom(
        parent: ParentType,
        element: str,
        isotope: int | None = None,
        stereo: str = '',
        hcount: int = 0,
        charge: int = 0,
        valance: int = None,
        **kwargs
) -> ParentType:
    atom = Atom(parent._get_id(), element, isotope, stereo, hcount, charge, valance,
                parent=parent, **kwargs)
    parent.nodes.append(atom)
    parent.root.atoms.append(atom)

    return parent


def add_bond(
        parent: ParentType,
        bond_symbol: str,
        atom1: Atom | BondDescriptorAtom,
        atom2: Atom | BondDescriptorAtom | StochasticObject | None,
        **kwargs
) -> ParentType:
    bond = Bond(bond_symbol, atom1, atom2, parent._get_id(), parent=parent, **kwargs)
    parent.nodes.append(bond)
    parent.root.bonds.append(bond)
    add_bond_to_connected_objects(bond)

    return parent


def add_ring(parent: ParentType, ring_id: int, bond_symbol: str = "", **kwargs) -> ParentType:
    # check if ring_id already exists
    ring_parent = get_ring_parent(parent)
    for ring in ring_parent.rings:
        if ring.ring_id == ring_id:
            if ring.atom2 is not None:
                raise BigSMILESError(f"Ring already formed for ring id {ring_id}.")
            atom2 = get_prior(parent, (Atom,))

            bond = get_common_bond(ring.atom1, atom2)
            if bond is None:
                # Add atom2
                ring.atom2 = atom2
                add_bond_to_connected_objects(ring)
            else:
                # Check if ring already between two atoms and just increase bond order
                bond.bond_order += 1
                remove_partial_ring(ring_parent, ring_id)

            return parent

    # Make new ring_id
    bond = Bond(bond_symbol, get_prior(parent, (Atom,)), None, parent._get_id(), ring_id,
                parent=parent, **kwargs)
    parent.root.bonds.append(bond)
    ring_parent.rings.append(bond)

    return parent


def add_ring_from_atoms(parent: ParentType, atom1: Atom, atom2: Atom, bond_symbol: str = "", **kwargs) \
        -> ParentType:
    ring_parent = get_ring_parent(parent)
    bond = get_common_bond(atom1, atom2)  # Check if ring already between two atoms
    if bond is None:
        # Make new ring_id
        bond = Bond(bond_symbol, atom1, atom2, parent._get_id(), len(parent.root.rings) + 1,
                    parent=parent, **kwargs)
        parent.root.bonds.append(bond)
        ring_parent.rings.append(bond)
        add_bond_to_connected_objects(bond)
    else:
        # increase bond order of existing bond
        bond.bond_order += 1

    return parent


def _get_bonding_descriptor_atom(parent, descriptor: str, index_: int, **kwargs) -> BondDescriptorAtom:
    bd = _get_bonding_descriptor(parent, descriptor, index_)
    return BondDescriptorAtom(bd, parent._get_id(), parent=parent, **kwargs)


def _get_bonding_descriptor(parent, descriptor: str, index_: int) -> BondDescriptor:
    stoch_obj = _get_current_stochastic_object(parent)

    # check if bd in there already
    for bd in stoch_obj.bonding_descriptors:
        if descriptor == bd.descriptor and index_ == bd.index_:
            return bd

    # create new bonding descriptor
    new_bd = BondDescriptor(stoch_obj, descriptor, index_)
    stoch_obj.bonding_descriptors.append(new_bd)
    return new_bd


def _get_current_stochastic_object(parent) -> StochasticObject:
    """ Recursive"""
    if isinstance(parent, StochasticObject):
        return parent

    if hasattr(parent, "parent"):
        return _get_current_stochastic_object(parent.parent)
    else:
        raise BigSMILESError("No stochastic object found.")


def add_bonding_descriptor(parent, descriptor: str, index_: int, **kwargs) -> BondDescriptorAtom:
    """ [<], [>], [$], [$1], [>2], ... """
    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_, **kwargs)
    parent.nodes.append(bd_atom)
    parent.root.atoms.append(bd_atom)
    return parent


def add_bond_atom_pair(
        parent,
        bond_symbol: str,
        element: str,
        isotope: int | None = None,
        stereo: str = '',
        hcount: int = 0,
        charge: int = 0,
        valance: int = None,
        atom_kwargs: dict = None,
        bond_kwargs: dict = None
) -> tuple[Bond, Atom]:
    atom_kwargs = atom_kwargs if atom_kwargs else {}
    bond_kwargs = bond_kwargs if bond_kwargs else {}

    atom = Atom(parent._get_id(), element, isotope, stereo, hcount, charge, valance,
                parent=parent, **atom_kwargs)
    prior_atom = get_prior(parent, (Atom, BondDescriptorAtom, StochasticObject))
    add_bond(parent, bond_symbol, prior_atom, atom, **bond_kwargs)
    parent.nodes.append(atom)
    parent.root.atoms.append(atom)

    return parent


def add_bond_bonding_descriptor_pair(
        parent,
        bond_symbol: str,
        descriptor: str,
        index_: int,
        bond_descriptor_kwargs: dict = None,
        bond_kwargs: dict = None
):
    bond_descriptor_kwargs = bond_descriptor_kwargs if bond_descriptor_kwargs else {}
    bond_kwargs = bond_kwargs if bond_kwargs else {}

    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_, **bond_descriptor_kwargs)
    prior_atom = get_prior(parent, (Atom, BondDescriptorAtom, StochasticObject))
    add_bond(parent, bond_symbol, prior_atom, bd_atom, **bond_kwargs)
    parent.nodes.append(bd_atom)
    parent.root.atoms.append(bd_atom)
    return parent


def open_branch(parent, **kwargs):
    if len(parent.nodes) == 0 and hasattr(parent, 'parent') and isinstance(parent.parent, Branch):
        raise BigSMILESError("BigSMILES string or branch can't start with branch")

    branch = Branch(parent, parent._get_id(), **kwargs)
    parent.nodes.append(branch)

    return branch


def close_branch(parent):
    if not isinstance(parent, Branch):
        raise BigSMILESError("Error closing branch. Possible issues: "
                             "\n\tClosing a branch with another intermediate node started."
                             "\n\tNo starting branch symbol.")

    # check if branch is empty; if so remove it
    if len(parent.nodes) == 0:
        parent.parent.nodes.pop()

    return parent.parent


def open_stochastic_object(parent, descriptor: str, index_: int, **kwargs) -> StochasticFragment:
    stoch_obj = StochasticObject(parent, parent._get_id())
    new_bd = BondDescriptor(stoch_obj, descriptor, index_)
    stoch_obj.bonding_descriptors.append(new_bd)
    stoch_obj.end_group_left = BondDescriptorAtom(new_bd, parent._get_id(),
                                                  parent=parent, **kwargs)
    parent.nodes.append(stoch_obj)

    fragment = StochasticFragment(stoch_obj, stoch_obj._get_id())
    stoch_obj.nodes.append(fragment)

    return fragment


def open_stochastic_object_with_bond(parent, bond_symbol: str, descriptor: str, index_: int, **kwargs) \
        -> StochasticFragment:
    stoch_obj = StochasticObject(parent, parent._get_id())
    new_bd = BondDescriptor(stoch_obj, descriptor, index_)
    stoch_obj.bonding_descriptors.append(new_bd)
    stoch_obj.end_group_left = BondDescriptorAtom(new_bd, parent._get_id(),
                                                  parent=parent, **kwargs)

    prior_atom = get_prior(parent, (Atom, BondDescriptor, StochasticObject))
    bond = add_bond(parent, bond_symbol, prior_atom, stoch_obj)
    stoch_obj.bond_left = bond

    parent.nodes.append(stoch_obj)

    return open_stochastic_fragment(stoch_obj)


def close_stochastic_object(parent, descriptor: str, index_: int):
    if not isinstance(parent, StochasticObject):
        raise BigSMILESError("Error closing StochasticObject. Possible issues: "
                             "\n\tClosing a StochasticObject with another intermediate node started."
                             "\n\tNo starting StochasticObject symbol.")

    parent.end_group_right = _get_bonding_descriptor_atom(parent, descriptor, index_)

    # TODO: add [H] and bond if no atom before or after
    return parent.parent


def open_stochastic_fragment(parent) -> StochasticFragment:
    stoch_frag = StochasticFragment(parent, parent._get_id())
    parent.nodes.append(stoch_frag)
    return stoch_frag


def close_open_stochastic_fragment(parent) -> StochasticFragment:
    if not isinstance(parent, StochasticFragment):
        raise BigSMILESError("Stochastic seperator can only follow a stochastic fragments.")

    parent = close_stochastic_fragment(parent)
    return open_stochastic_fragment(parent)


def close_stochastic_fragment(parent):
    if not isinstance(parent, StochasticFragment):
        raise BigSMILESError("Error StochasticFragment branch. Possible issues: "
                             "\n\tClosing a branch with another intermediate node started.")

    add_bond_descriptor_to_stochastic_fragment(parent)

    # check for at-least one bonding descriptor
    if not parent.bonding_descriptors:
        raise BigSMILESError(f"No bonding descriptor in the stochastic fragment. ({repr(parent)})")

    return parent.parent


## functions for building BigSMILES in chunks ##
###################################################################################################################
def append_bigsmiles_fragment(parent, bond_symbol: str | None, bigsmiles_: BigSMILES):
    if bond_symbol is not None:
        if not isinstance(bigsmiles_.nodes[0], Atom | StochasticObject):
            raise BigSMILESError("First node must be an 'Atom' to added fragment.")

        atom = bigsmiles_.nodes[0]
        prior_atom = get_prior(parent, (Atom, BondDescriptorAtom, StochasticObject))
        add_bond(parent, bond_symbol, prior_atom, atom)
        # stack[-1].nodes.append(atom)
        # bigsmiles.atoms.append(atom)

    # append bigsmiles
    # re-direct 'parent' to new bigsmiles object
    for node in bigsmiles_.nodes:
        if hasattr(node, 'parent'):
            node.parent = parent
    # re-numbering
    # re_number_rings(parent, set())
    # re_number_nodes(bigsmiles_)

    # append fragment
    parent.nodes += bigsmiles_.nodes
    parent.root.atoms += bigsmiles_.atoms
    parent.root.bonds += bigsmiles_.bonds
    parent.root.rings += bigsmiles_.rings


def attach_bigsmiles_branch(parent, bond_symbol: str | None, bigsmiles_: BigSMILES, index_: int):
    if index_ > len(parent.root.atoms) - 1:
        raise BigSMILESError(f"Branch can't be added. 'index' outside atom list."
                             f"\n current smiles: {parent.root} \n desired atom index: {index_}")

    # re-direct 'parent' to new bigsmiles object
    for node in bigsmiles_.nodes:
        if hasattr(node, 'parent'):
            node.parent = parent
    # # re-numbering
    # re_number_rings(parent, set())
    # re_number_nodes(bigsmiles_)

    atom = parent.root.atoms[index_]

    # make branch
    branch = Branch(parent, parent._get_id())
    branch.nodes = bigsmiles_.nodes
    parent.root.atoms += bigsmiles_.atoms
    parent.root.bonds += bigsmiles_.bonds

    # make bond
    bond = Bond(bond_symbol, atom, branch.nodes[0], parent._get_id(), parent=parent)
    branch.nodes.insert(0, bond)
    parent.root.bonds.append(bond)
    add_bond_to_connected_objects(bond)

    # insert branch by put it after all the other branches
    insert_obj_into_list(parent.root.nodes, atom, branch)


# def re_number_rings(obj, seen: set[Bond]):
#     """ Recursive renumbering of rings. """
#     for node in obj.nodes:
#         if isinstance(node, Bond) and node.ring_id is not None and node not in seen:
#             obj.root.rings.append(node)
#             node.ring_id = len(obj.root.rings)
#             seen.add(node)
#         if hasattr(node, 'nodes'):
#             re_number_rings(node, seen)
#
#
# def re_number_nodes(obj):
#     """ Recursive renumbering 'id_'. """
#     for node in obj.nodes:
#         if isinstance(node, Atom):
#             node.id_ = parent._get_id()
#         elif isinstance(node, Bond):
#             node.id_ = parent._get_id()
#         elif isinstance(node, BondDescriptorAtom):
#             node.id_ = parent._get_id()
#         elif isinstance(node, Branch):
#             node.id_ = parent._get_id()
#             re_number_nodes(node)
#         elif isinstance(node, StochasticFragment):
#             node.id_ = parent._get_id()
#             re_number_nodes(node)
#         elif isinstance(node, StochasticObject):
#             node.id_ = parent._get_id()
#             re_number_nodes(node)


def add_bond_bonding_descriptor_via_index(
        parent,
        bond_symbol: str,
        descriptor: str,
        index_: int,
        prior_atom: Atom
):
    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_)
    bond = add_bond(parent, bond_symbol, prior_atom, bd_atom)

    if prior_atom == prior_atom.parent.nodes[0]:
        parent.nodes.insert(0, bond)
        parent.bonds.insert(0, bond)
        parent.nodes.insert(0, bd_atom)
        parent.atoms.insert(0, bd_atom)
    elif prior_atom == prior_atom.parent.nodes[-1]:
        parent.nodes.append(bond)
        parent.bonds.append(bond)
        parent.nodes.append(bd_atom)
        parent.atoms.append(bd_atom)
    else:
        branch = Branch(parent, parent._get_id())
        branch.nodes = [bond, bd_atom]
        branch.atoms = [bd_atom]
        branch.bond = [bond]
        insert_obj_into_list(parent.nodes, prior_atom, branch)

    insert_obj_into_list(parent.root.bonds, prior_atom.bonds[0], bond)
    insert_obj_into_list(parent.root.atoms, prior_atom, bd_atom)


# def add_stochastic_object(
#         parent: BigSMILES | Branch,
#         left_descriptor: str,
#         left_index: int,
#         right_descriptor: str,
#         right_index: int
# ):
#     stoch_obj = StochasticObject(parent, parent._get_id())
#
#     # create bonding descriptors
#     left_bd = BondDescriptor(stoch_obj, left_descriptor, left_index)
#     right_bd = BondDescriptor(stoch_obj, right_descriptor, right_index)
#     stoch_obj.bonding_descriptors.append(left_bd)
#     stoch_obj.bonding_descriptors.append(right_bd)
#
#     if not stoch_obj.implicit_endgroups:
#         try:
#             prior_atom = get_prior(parent, (Atom, BondDescriptor, StochasticObject))
#         except BigSMILESError:
#             # no prior atoms, it must be the first atom
#             if isinstance(parent, BigSMILES) and:
#                 prior_atom = Atom(parent._get_id(), "H", None, "", 0, 0, 1, parent)
#                 parent.atoms.append(prior_atom)
#
#         bond = add_bond(parent, "", prior_atom, stoch_obj)
#         stoch_obj.bond_left = bond
#
#     parent.nodes.append(stoch_obj)
#
#     return stoch_obj


__all__ = ["ParentType", "add_atom", "add_bond", "add_ring", "add_ring_from_atoms", "add_bonding_descriptor",
           "add_bond_atom_pair", "add_bond_bonding_descriptor_pair", "open_branch", "close_branch",
           "open_stochastic_object", "open_stochastic_object_with_bond", "close_stochastic_object",
           "open_stochastic_fragment", "close_open_stochastic_fragment", "close_stochastic_fragment",
           "exit_", "run_validation"
           ]
