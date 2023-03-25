"""
Constructor:

This file contains functions that can be called to construct a BigSMILES object.

There are two approaches:
* build it step by step;  where you open branch add atoms then close branch
* build it in chunks;  where you build the branch in a separate BigSMILES object and attach it.

"""
from __future__ import annotations
import logging

import bigsmiles.errors as errors
import bigsmiles.reference_data.chemical_data as chemical_data
from bigsmiles.data_structures.bigsmiles import Atom, Bond, BondDescriptor, Branch, StochasticFragment, \
    StochasticObject, BigSMILES, BondDescriptorAtom, has_node_attr, has_ring_attr, has_parent_attr
import bigsmiles.constructors.validation.validation_bigsmiles_obj as validation_bigsmiles_obj


## Utility functions for construction ## noqa
######################################################################################################################
def add_bond_to_connected_objects(bond: Bond):
    """
    Adds Bond to Atoms.bonds attribute
    also works on BondDescriptorAtom and stochastic objects
    """
    if bond.ring_id is not None and isinstance(bond.atom2, StochasticObject):
        atoms = reversed(bond)
    else:
        atoms = iter(bond)

    for i, obj in enumerate(atoms):
        if isinstance(obj, Atom):
            obj.bonds += [bond]
        elif isinstance(obj, BondDescriptorAtom):
            obj.bond = bond
        elif isinstance(obj, StochasticObject):
            if i == 0:
                obj.bond_right = bond
            else:
                obj.bond_left = bond


def add_bond_descriptor_to_stochastic_fragment(stoch_frag: StochasticFragment, loop_obj: Branch = None):
    """ Add bonding descriptor to stoch_frag.bonding_descriptors. Recursive. """
    if loop_obj is None:
        loop_obj = stoch_frag

    for obj in loop_obj.nodes:
        if isinstance(obj, BondDescriptorAtom) and obj.descriptor not in stoch_frag.bonding_descriptors:
            stoch_frag.bonding_descriptors.append(obj.descriptor)
        if isinstance(obj, Branch):
            add_bond_descriptor_to_stochastic_fragment(stoch_frag, obj)  # recursive


def get_common_bond(atom1: Atom, atom2: Atom | StochasticObject) -> Bond | None:
    """ Checks if the two atoms have a bond in common. """
    if isinstance(atom2, StochasticObject):
        if atom2.bond_right is not None and atom2.bond_right in atom1.bonds:
            raise errors.ConstructorError("An atom can't have two bonds to the same stochastic object.")
    else:
        for bond in atom1.bonds:
            if bond in atom2.bonds:
                return bond

    return None  # No common bond found


def remove_partial_ring(parent: has_ring_attr, ring_id: int):
    """ Given a ring_id, remove ring from parent. """
    for i, ring in enumerate(parent.rings):
        if ring.ring_id == ring_id:
            bond = parent.rings.pop(i)
            parent.bonds.remove(bond)
            if hasattr(parent, "parent"):
                parent.parent.bonds.remove(bond)
            break


def insert_obj_into_list(list_: list, obj, object_to_insert):
    """ Insert 'object_to_insert' in a list right behind 'obj'. """
    index_ = list_.index(obj) + 1
    for i, item in enumerate(list_[index_:]):
        if isinstance(item, Branch):  # put it at the end of existing branches
            continue

        list_.insert(index_ + i, object_to_insert)
        break


def get_prior(obj: has_node_attr, types_: tuple, flag: bool = False) -> Atom:
    """ Get prior """
    if obj.nodes:
        for node in reversed(obj.nodes):  # no atom can ever have more than 10 bonds
            if type(node) in types_:
                return node

        raise errors.ConstructorError("Trying to make a bond to a prior atom that isn't there.??")

    elif not flag:
        return get_prior(obj.parent, types_, flag=True)

    raise errors.ConstructorError("Bond attempted to be made to that has nothing to bond back to.")


def get_ring_parent(parent: has_parent_attr) -> BigSMILES | StochasticFragment:
    """ Recursively find a parent that has the 'ring' attribute. only BigSMILES and Stochastic Fragments do. """
    if hasattr(parent, "rings"):
        return parent

    return get_ring_parent(parent.parent)


def add_explict_hydrogen_to_stochastic_objects(bigsmiles_: BigSMILES):
    """ {[$][$]CC[$][$]} --> [H]{[$][$]CC[$][$]}[H] """
    if isinstance(bigsmiles_.nodes[0], StochasticObject) and not bigsmiles_.nodes[0].implicit_endgroups:
        if bigsmiles_.nodes[0].bd_left.bond_order > 1:
            raise errors.ConstructorError(
                f"Invalid BigSMILES. A double bond coming out of a stochastic object requires an "
                f"explict endgroup defined. \n Current Bigsmiles: {bigsmiles_}")

        insert_atom_and_bond(bigsmiles_, None, "", "H")

    if isinstance(bigsmiles_.nodes[-1], StochasticObject) and not bigsmiles_.nodes[-1].implicit_endgroups and \
            bigsmiles_.nodes[-1].bond_right is None:
        if bigsmiles_.nodes[-1].bd_right.bond_order > 1:
            raise errors.ConstructorError(
                f"Invalid BigSMILES. A double bond coming out of a stochastic object requires an "
                f"explict endgroup defined. \n Current Bigsmiles: {bigsmiles_}")

        add_bond_atom_pair(bigsmiles_, "", "H")


def exit_construction(parent: BigSMILES, syntax_fix: bool = True, validation: bool = True):
    """
    runs several methods to fix any syntax issues and runs more extensive validation

    Parameters
    ----------
    parent:
        BigSMILES to be checked
    syntax_fix:
        True to perform syntax fixes
    validation:
        True to run validation
    """
    if not isinstance(parent, BigSMILES):
        raise errors.ConstructorError(f"{type(parent)} is missing closing symbol.")

    if not parent:
        return parent
    if syntax_fix:
        add_explict_hydrogen_to_stochastic_objects(parent)
        validation_bigsmiles_obj.run_syntax_fixes(parent)
    if validation:
        validation_bigsmiles_obj.run_validation(parent)


## Functions for constructing a bigsmiles sequentially or step by step ## noqa
#######################################################################################################################
def add_atom(
        parent: has_node_attr,
        symbol: str,
        isotope: int | None = None,
        stereo: str | None = None,
        hydrogens: int | None = None,
        charge: int = 0,
        valence: int | None = None,
        class_: int | None = None,
        **kwargs
) -> has_node_attr:
    """
    appends an `Atom` to the end of the `parent`

    Parameters
    ----------
    parent:
    symbol:
        symbol symbol (e.g., H, C, O, Zn)
    isotope:
        isotope (e.g., [13C])
    stereo:
        stereochemistry [None, "@", "@@"] (None means not defined)
    hydrogens:
        number of explict hydrogens  (e.g., [CH2])
    charge:
        symbol charge  (e.g., [Fe+3])
    valence:
        The capacity to form bonds with other atoms
    class_:
        index of class (e.g., [C:1] class_ = 1)
    kwargs:
        any additional keyword arguments are accepted and set as additional attributes to the Atom

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    atom = Atom(parent._get_id(Atom), symbol, isotope, stereo, hydrogens, charge, valence, class_,
                parent=parent, **kwargs)
    parent.nodes.append(atom)
    parent.root.atoms.append(atom)

    return parent


def add_bond(
        parent: has_node_attr,
        bond_symbol: str | None,
        atom1: Atom | BondDescriptorAtom,
        atom2: Atom | BondDescriptorAtom | StochasticObject | None,
        **kwargs
) -> has_node_attr:
    """
    add 'Bond' between two 'Atoms' (typically not be used directly)

    Parameters
    ----------
    parent:
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    atom1:
        atom the bond starts at
    atom2:
        atom the bond ends at
    kwargs:
        any additional keyword arguments are accepted and set as additional attributes to the Bond

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    bond = Bond(parent._get_id(Bond), bond_symbol, atom1, atom2, parent=parent, **kwargs)
    parent.nodes.append(bond)
    parent.root.bonds.append(bond)
    add_bond_to_connected_objects(bond)

    return parent


def add_ring_by_index(parent: has_node_attr, ring_id: int, bond_symbol: str | None, **kwargs) -> has_node_attr:
    """
    add ring by index; rings are just 'Bond'

    !!! warning

        This function needs to be called twice to work; once at each time the parsing hits a ring index.

        'exit_construction()' will check to confirm that there is no uncompleted rings

    Parameters
    ----------
    parent:
    ring_id:
        index of ring
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    kwargs:
        any additional keyword arguments are accepted and set as additional attributes to the Bond

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    # check if ring_id already exists
    ring_parent = get_ring_parent(parent)
    for ring in ring_parent.rings:
        if ring.ring_id == ring_id:
            if ring.atom2 is not None:
                raise errors.ConstructorError(f"Ring already formed for ring id {ring_id}.")
            atom2 = get_prior(parent, (Atom, StochasticObject))

            # validation
            bond = get_common_bond(ring.atom1, atom2)
            if bond:
                logging.warning("Duplicate ring detected and merged into one with a higher bond order. ")
                bond.bond_order += chemical_data.bond_mapping[bond_symbol]
                remove_partial_ring(ring_parent, ring_id)
                return parent

            # Add atom2
            # ring.atom1, ring.atom2 = atom2, ring.atom1
            ring.atom2 = atom2

            # this is to handle that the multi-bond symbol (-, #) could be on either or both ring index
            if chemical_data.bond_mapping[bond_symbol] > ring.bond_order:
                ring.symbol = bond_symbol

            # aromatic_elements ring closure
            if ring.atom1.aromatic and ring.atom2.aromatic:
                ring.symbol = ":"

            add_bond_to_connected_objects(ring)
            return parent

    # Make new ring_id
    bond = Bond(parent._get_id(Bond), bond_symbol, get_prior(parent, (Atom, StochasticObject)), None, ring_id,
                parent=parent, **kwargs)
    parent.root.bonds.append(bond)
    ring_parent.rings.append(bond)

    return parent


def add_ring_from_atoms(parent: has_node_attr, atom1: Atom, atom2: Atom, bond_symbol: str | None = "", **kwargs) \
        -> has_node_attr:
    """
    add ring from with both 'Atoms'; rings are just 'Bond'

    Parameters
    ----------
    parent:
    atom1:
        atom the ring starts at
    atom2:
        atom the ring ends at
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    kwargs
        any additional keyword arguments are accepted and set as additional attributes to the Bond

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    ring_parent = get_ring_parent(parent)
    bond = get_common_bond(atom1, atom2)  # Check if ring already between two atoms
    if bond is None:
        # Make new ring_id
        bond = Bond(parent._get_id(Bond), bond_symbol, atom1, atom2, len(parent.root.rings) + 1,
                    parent=parent, **kwargs)
        parent.root.bonds.append(bond)
        ring_parent.rings.append(bond)
        add_bond_to_connected_objects(bond)
    else:
        # increase bond order of existing bond
        bond.bond_order += 1

    return parent


def _get_bonding_descriptor_atom(parent, descriptor: str, index_: int, bond_symbol: str | None, **kwargs) \
        -> BondDescriptorAtom:
    stoch_obj = _get_current_stochastic_object(parent)
    bd = _get_bonding_descriptor(stoch_obj, descriptor, index_, bond_symbol)
    return BondDescriptorAtom(parent._get_id(BondDescriptorAtom), bd, parent=parent, **kwargs)


def _get_bonding_descriptor(stoch_obj: StochasticObject, descriptor: str, index_: int, bond_symbol: str | None) \
        -> BondDescriptor:
    # check if bd in there already
    for bd in stoch_obj.bonding_descriptors:
        if descriptor == bd.descriptor and index_ == bd.index_:
            if bond_symbol is not None and bd.bond_symbol != bond_symbol:
                raise errors.ConstructorError("Multiple bond orders to same bonding descriptor.")
                # logging.warning("Multiple bond orders to same bonding descriptor.")
            return bd

    # create new bonding descriptor
    new_bd = BondDescriptor(stoch_obj, descriptor, index_, bond_symbol)
    stoch_obj.bonding_descriptors.append(new_bd)
    return new_bd


def _get_current_stochastic_object(parent) -> StochasticObject:
    """ Recursive"""
    if isinstance(parent, StochasticObject):
        return parent

    if hasattr(parent, "parent"):
        return _get_current_stochastic_object(parent.parent)
    else:
        raise errors.ConstructorError("No stochastic object found.")


def add_bonding_descriptor_atom(parent, descriptor: str, index_: int, bond_symbol: str | None = None, **kwargs) \
        -> has_node_attr:
    """
    appends 'BondDescriptorAtom' to the end of the `parent`

    Parameters
    ----------
    parent:
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    kwargs
        any additional keyword arguments are accepted and set as additional attributes to the BondDescriptorAtom

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_, bond_symbol, **kwargs)
    parent.nodes.append(bd_atom)
    parent.root.atoms.append(bd_atom)
    return parent


def add_bond_atom_pair(
        parent: has_node_attr,
        bond_symbol: str | None,
        symbol: str,
        isotope: int | None = None,
        stereo: str | None = None,
        hydrogens: int | None = None,
        charge: int = 0,
        valence: int = None,
        class_: int | None = None,
        kwargs_atom: dict = None,
        kwargs_bond: dict = None
) -> has_node_attr:
    """
    appends a `Bond` followed by an `Atom`

    Parameters
    ----------
    parent:
        parent
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    symbol:
        symbol symbol (e.g., H, C, O, Zn)
    isotope:
        isotope (e.g., [13C])
    stereo:
        stereochemistry [None, "@", "@@"] (None means not defined)
    hydrogens:
        number of explict hydrogens  (e.g., [CH2])
    charge:
        symbol charge  (e.g., [Fe+3])
    valence:
        The capacity to form bonds with other atoms
    class_:
        index of class (e.g., [C:1] class_ = 1)
    kwargs_atom:
        any additional keyword arguments are accepted and set as additional attributes to the Atom
    kwargs_bond
        any additional keyword arguments are accepted and set as additional attributes to the Bond

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    kwargs_atom = kwargs_atom if kwargs_atom is not None else {}
    kwargs_bond = kwargs_bond if kwargs_bond is not None else {}

    atom = Atom(parent._get_id(Atom), symbol, isotope, stereo, hydrogens, charge, valence, class_,
                parent=parent, **kwargs_atom)
    prior_atom = get_prior(parent, (Atom, BondDescriptorAtom, StochasticObject))
    add_bond(parent, bond_symbol, prior_atom, atom, **kwargs_bond)
    parent.nodes.append(atom)
    parent.root.atoms.append(atom)

    return parent


def add_bond_bonding_descriptor_pair(
        parent: StochasticFragment | Branch,
        bond_symbol: str | None,
        descriptor: str,
        index_: int,
        kwargs_bond_descriptor: dict = None,
        kwargs_bond: dict = None
) -> StochasticFragment | Branch:
    """
    appends a `Bond` followed by an `BondingDescriptorAtom`

    Parameters
    ----------
    parent:
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    kwargs_bond_descriptor:
        any additional keyword arguments are accepted and set as additional attributes to the BondDescriptorAtom
    kwargs_bond:
        any additional keyword arguments are accepted and set as additional attributes to the Bond

    Returns
    -------
    parent:
        returns the parent passed to it
    """
    kwargs_bond_descriptor = kwargs_bond_descriptor if kwargs_bond_descriptor is not None else {}
    kwargs_bond = kwargs_bond if kwargs_bond is not None else {}

    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_, bond_symbol, **kwargs_bond_descriptor)
    prior_atom = get_prior(parent, (Atom, BondDescriptorAtom, StochasticObject))
    add_bond(parent, bond_symbol, prior_atom, bd_atom, **kwargs_bond)
    parent.nodes.append(bd_atom)
    parent.root.atoms.append(bd_atom)
    return parent


def open_branch(parent: has_node_attr, **kwargs) -> Branch:
    """
    open branch

    !!! warning

        close_branch() must be called at some later point

    Parameters
    ----------
    parent:
    kwargs:
        any additional keyword arguments are accepted and set as additional attributes to the Branch

    Returns
    -------
    parent:
        returns Branch
    """
    if len(parent.nodes) == 0 and hasattr(parent, 'parent') and isinstance(parent.parent, Branch):
        raise errors.ConstructorError("BigSMILES string or branch can't start with branch")

    branch = Branch(parent._get_id(Branch), parent, **kwargs)
    parent.nodes.append(branch)

    return branch


def close_branch(parent: Branch) -> has_node_attr:
    """
    closes branch and appends it to `parent.parent`

    Parameters
    ----------
    parent:

    Returns
    -------
    parent:
        returns parent.parent

    """
    if not isinstance(parent, Branch):
        raise errors.ConstructorError("Error closing branch. Possible issues: "
                                      "\n\tClosing a branch with another intermediate node started."
                                      "\n\tNo starting branch symbol.")

    # check if branch is empty; if so remove it
    if len(parent.nodes) == 0:
        parent.parent.nodes.pop()

    return parent.parent


def open_stochastic_object(parent: has_node_attr, descriptor: str, index_: int, **kwargs) -> StochasticObject:
    """
    open `StochasticObject`

    !!! warning

        close_stochastic_object() must be called at a later point

    Parameters
    ----------
    parent:
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    kwargs
        any additional keyword arguments are accepted and set as additional attributes to the BondDescriptorAtom

    Returns
    -------
    parent:
        returns StochasticObject
    """
    stoch_obj = StochasticObject(parent._get_id(StochasticObject), parent, **kwargs)
    stoch_obj.bd_left = _get_bonding_descriptor(stoch_obj, descriptor, index_, None)
    parent.nodes.append(stoch_obj)

    return stoch_obj


def open_stochastic_object_fragment(parent: has_node_attr, descriptor: str, index_: int, **kwargs) \
        -> StochasticFragment:
    """
    opens `StochasticObject` and `StochasticFragment`

    !!! warning

        close_stochastic_fragment() must be called at a later point

        close_stochastic_object() must be called at a later point

    Parameters
    ----------
    parent:
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    kwargs
        any additional keyword arguments are accepted and set as additional attributes to the BondDescriptorAtom

    Returns
    -------
    parent:
        returns StochasticFragment
    """
    stoch_obj = open_stochastic_object(parent, descriptor, index_, **kwargs)

    return open_stochastic_fragment(stoch_obj)


def open_stochastic_object_fragment_with_bond(parent: has_node_attr, bond_symbol: str | None, descriptor: str,
                                              index_: int, **kwargs) -> StochasticFragment:
    """
    opens `StochasticObject` and `StochasticFragment`

    !!! warning

        close_stochastic_fragment() must be called at a later point

        close_stochastic_object() must be called at a later point

    Parameters
    ----------
    parent:
    bond_symbol:
        bond symbol ("", "=", "#", etc.)
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    kwargs
        any additional keyword arguments are accepted and set as additional attributes to the BondDescriptorAtom

    Returns
    -------
    parent:
        returns StochasticFragment
    """
    stoch_obj = StochasticObject(parent._get_id(StochasticObject), parent, **kwargs)
    stoch_obj.bd_left = _get_bonding_descriptor(stoch_obj, descriptor, index_, bond_symbol)

    # bond made without 'add_bond' function to ensure it's added to 'bond_left'
    prior_atom = get_prior(parent, (Atom, StochasticObject))
    bond = Bond(parent._get_id(Bond), bond_symbol, prior_atom, stoch_obj, parent=parent)
    parent.nodes.append(bond)
    parent.root.bonds.append(bond)
    add_bond_to_connected_objects(bond)
    parent.nodes.append(stoch_obj)

    return open_stochastic_fragment(stoch_obj)


def close_stochastic_object(parent: StochasticObject, descriptor: str, index_: int, bond_symbol: str | None) \
        -> has_node_attr:
    """
    closes `StochasticObject`

    Parameters
    ----------
    parent:
    descriptor:
        bonding descriptor symbol [<, >, $]
    index_:
        bonding descriptor index [1, inf]
    bond_symbol:
        bond symbol ("", "=", "#", etc.)

    Returns
    -------
    parent:
        returns parent.parent

    Raises
    ------
    ConstructorError
        Error closing StochasticObject. Possible issues:

        * Closing a StochasticObject with another intermediate node started.
        * No starting StochasticObject symbol.
    """
    if not isinstance(parent, StochasticObject):
        raise errors.ConstructorError("Error closing StochasticObject. Possible issues: "
                                      "\n\tClosing a StochasticObject with another intermediate node started."
                                      "\n\tNo starting StochasticObject symbol.")

    parent.bd_right = _get_bonding_descriptor(parent, descriptor, index_, bond_symbol)
    return parent.parent


def open_stochastic_fragment(parent: StochasticObject) -> StochasticFragment:
    """
    open `StochasticFragment`

    !!! warning

        close_stochastic_fragment() must be called at a later point


    Parameters
    ----------
    parent:

    Returns
    -------
    parent:
        returns StochasticFragment
    """
    stoch_frag = StochasticFragment(parent._get_id(StochasticFragment), parent)
    parent.nodes.append(stoch_frag)
    return stoch_frag


def close_open_stochastic_fragment(parent: StochasticFragment) -> StochasticFragment:
    """

    closes `StochasticFragment` and opens a new `StochasticFragment`

    Parameters
    ----------
    parent:

    Returns
    -------
    parent:
        returns a new StochasticFragment
    """
    if not isinstance(parent, StochasticFragment):
        raise errors.ConstructorError("Stochastic seperator can only follow a stochastic fragments.")

    parent = close_stochastic_fragment(parent)
    return open_stochastic_fragment(parent)


def close_stochastic_fragment(parent: StochasticFragment) -> StochasticObject:
    """

    Parameters
    ----------
    parent:

    Returns
    -------
    parent:
        returns parent.parent
    """
    if not isinstance(parent, StochasticFragment):
        raise errors.ConstructorError("Error StochasticFragment branch. Possible issues: "
                                      "\n\tClosing a branch with another intermediate node started.")

    add_bond_descriptor_to_stochastic_fragment(parent)

    # check for at-least one bonding descriptor
    if not parent.bonding_descriptors:
        raise errors.ConstructorError(f"No bonding descriptor in the stochastic fragment. ({repr(parent)})")

    return parent.parent


## functions for building BigSMILES in chunks ## noqa
###################################################################################################################
def append_bigsmiles_fragment(parent, bigsmiles_: BigSMILES, bond_symbol: str | None, **kwargs) -> BigSMILES:
    """

    Parameters
    ----------
    parent
    bigsmiles_
    bond_symbol
    kwargs

    Returns
    -------

    """
    if not isinstance(bigsmiles_.nodes[0], (Atom, StochasticObject)):
        raise errors.ConstructorError("First node must be an 'Atom' or 'StochasticObject' to added fragment.")
    if not parent:
        return bigsmiles_

    atom = bigsmiles_.nodes[0]
    prior_atom = get_prior(parent, (Atom, StochasticObject))
    add_bond(parent, bond_symbol, prior_atom, atom, **kwargs)

    set_new_parent(parent, bigsmiles_)

    # append fragment
    parent.nodes += bigsmiles_.nodes
    parent.root.atoms += bigsmiles_.atoms
    parent.root.bonds += bigsmiles_.bonds
    parent.root.rings += bigsmiles_.rings

    return parent


def attach_bigsmiles_branch(parent, bond_symbol: str | None | None, bigsmiles_: BigSMILES, index_: int,
                            kwargs_bond: dict = None, kwargs_branch: dict = None):
    """

    Parameters
    ----------
    parent
    bond_symbol
    bigsmiles_
    index_
    kwargs_bond
    kwargs_branch

    Returns
    -------

    """
    kwargs_bond = kwargs_bond if kwargs_bond is not None else {}
    kwargs_branch = kwargs_branch if kwargs_branch is not None else {}

    if index_ > len(parent.root.atoms) - 1:
        raise errors.ConstructorError(f"Branch can't be added. 'index' outside atom list."
                                      f"\n current smiles: {parent.root} \n desired atom index: {index_}")

    set_new_parent(parent, bigsmiles_)

    # make branch
    branch = Branch(parent._get_id(Branch), parent, **kwargs_branch)
    branch.nodes = bigsmiles_.nodes
    parent.root.atoms += bigsmiles_.atoms
    parent.root.bonds += bigsmiles_.bonds

    # make bond
    atom = parent.root.atoms[index_]
    bond = Bond(parent._get_id(Bond), bond_symbol, atom, branch.nodes[0], parent=parent, **kwargs_bond)
    branch.nodes.insert(0, bond)
    parent.root.bonds.append(bond)
    add_bond_to_connected_objects(bond)

    # insert branch (put it after all the other branches on the atom)
    insert_obj_into_list(parent.root.nodes, atom, branch)

    return parent


def insert_atom_and_bond(
        parent: has_node_attr,
        prior_atom: Atom | BondDescriptorAtom | Branch | StochasticObject | None,
        bond_symbol: str | None,
        element: str,
        isotope: int | None = None,
        stereo: str = '',
        hydrogens: int = 0,
        charge: int = 0,
        valance: int = None,
        class_: int | None = None,
        kwargs_atom: dict = None,
        kwargs_bond: dict = None
) -> has_node_attr:
    """

    Parameters
    ----------
    parent
    prior_atom
    bond_symbol
    element
    isotope
    stereo
    hydrogens
    charge
    valance
    class_
    kwargs_atom
    kwargs_bond

    Returns
    -------

    """
    kwargs_atom = kwargs_atom if kwargs_atom is not None else {}
    kwargs_bond = kwargs_bond if kwargs_bond is not None else {}

    atom = Atom(parent._get_id(Atom), element, isotope, stereo, hydrogens, charge, valance, class_,
                parent=parent, **kwargs_atom)

    if prior_atom is None and isinstance(parent, BigSMILES):  # add to beginning
        prior_atom = parent.nodes[0]
        bond = Bond(parent._get_id(Bond), bond_symbol, atom, prior_atom, parent=parent, **kwargs_bond)
        add_bond_to_connected_objects(bond)
        parent.nodes.insert(0, bond)
        parent.root.bonds.insert(0, bond)
        parent.nodes.insert(0, atom)
        parent.root.atoms.insert(0, atom)

    elif prior_atom == prior_atom.parent.nodes[-1]:  # add to end
        bond = Bond(parent._get_id(Bond), bond_symbol, parent.atoms[-1], atom, parent=parent, **kwargs_bond)
        add_bond_to_connected_objects(bond)
        parent.nodes.append(bond)
        parent.root.bonds.append(bond)
        parent.nodes.append(atom)
        parent.root.atoms.append(atom)

    validation_bigsmiles_obj.re_number_node_ids(parent.root)
    return parent


def insert_atom_into_bond(
        parent: has_node_attr,
        bond_to_insert_into: Bond,
        element: str,
        isotope: int | None = None,
        stereo: str = '',
        hydrogens: int = 0,
        charge: int = 0,
        valance: int = None,
        class_: int | None = None,
        kwargs_atom: dict = None,
) -> has_node_attr:
    """

    Parameters
    ----------
    parent
    bond_to_insert_into
    element
    isotope
    stereo
    hydrogens
    charge
    valance
    class_
    kwargs_atom

    Returns
    -------

    """
    kwargs_atom = kwargs_atom if kwargs_atom is not None else {}

    # create atom
    atom = Atom(parent._get_id(Atom), element, isotope, stereo, hydrogens, charge, valance, class_,
                parent=parent, **kwargs_atom)

    # create bonds
    new_bond1 = Bond(parent._get_id(Bond), bond_to_insert_into.symbol, bond_to_insert_into.atom1, atom, parent=parent)
    new_bond2 = Bond(parent._get_id(Bond), bond_to_insert_into.symbol, atom, bond_to_insert_into.atom2, parent=parent)

    # remove old bond and add new bonds and atoms to parent
    bond_index = parent.nodes.index(bond_to_insert_into)
    parent.nodes.insert(bond_index + 1, new_bond1)
    parent.nodes.insert(bond_index + 2, atom)
    parent.nodes.insert(bond_index + 3, new_bond2)
    parent.nodes.pop(bond_index)

    # remove old bond and add new bonds and atoms to root
    bond_index = parent.root.bonds.index(bond_to_insert_into)
    parent.root.bonds.insert(bond_index + 1, new_bond1)
    parent.root.bonds.insert(bond_index + 2, new_bond2)
    parent.root.bonds.pop(bond_index)
    atom_index = parent.root.atoms.index(bond_to_insert_into.atom1)
    parent.root.atoms.insert(atom_index + 1, atom)

    # remove old bond from atoms and add new bonds to atoms
    bond_to_insert_into.atom1.bonds.remove(bond_to_insert_into)
    bond_to_insert_into.atom2.bonds.remove(bond_to_insert_into)
    add_bond_to_connected_objects(new_bond1)
    add_bond_to_connected_objects(new_bond2)

    return parent


def add_bonding_descriptor_bond_via_index(
        parent: StochasticFragment,
        bond_symbol: str | None,
        descriptor: str,
        index_: int,
        prior_atom: Atom,
        kwargs_bond: dict = None,
):
    """

    Parameters
    ----------
    parent
    bond_symbol
    descriptor
    index_
    prior_atom
    kwargs_bond

    Returns
    -------

    """
    kwargs_bond = kwargs_bond if kwargs_bond is not None else {}

    bd_atom = _get_bonding_descriptor_atom(parent, descriptor, index_, bond_symbol)
    bond = Bond(parent._get_id(Bond), bond_symbol, prior_atom, bd_atom, parent=parent, **kwargs_bond)
    add_bond_to_connected_objects(bond)

    if prior_atom == prior_atom.parent.nodes[0]:  # if first atom
        parent.nodes.insert(0, bond)
        parent.root.bonds.insert(0, bond)
        parent.nodes.insert(0, bd_atom)
        parent.root.atoms.insert(0, bd_atom)
    elif prior_atom == prior_atom.parent.nodes[-1] or \
            all([False for node in parent.nodes[parent.nodes.index(prior_atom) + 1:] if not isinstance(node, Branch)]):
        # if last atom or last atom followed by only branches
        parent.nodes.append(bond)
        parent.root.bonds.append(bond)
        parent.nodes.append(bd_atom)
        parent.root.atoms.append(bd_atom)
    else:
        branch = Branch(parent._get_id(Branch), parent)
        branch.nodes = [bond, bd_atom]
        insert_obj_into_list(parent.nodes, prior_atom, branch)
        insert_obj_into_list(parent.root.bonds, prior_atom.bonds[0], bond)
        insert_obj_into_list(parent.root.atoms, prior_atom, bd_atom)


def add_stochastic_fragment(stoch_object: StochasticObject, stoch_fragment: StochasticFragment) -> StochasticObject:
    """

    Parameters
    ----------
    stoch_object
    stoch_fragment

    Returns
    -------

    """
    stoch_object.nodes.append(stoch_fragment)
    stoch_fragment.parent = stoch_object

    for node in stoch_fragment.nodes:
        if isinstance(node, Atom):
            stoch_object.root.atoms.append(node)
        elif isinstance(node, Bond):
            stoch_object.root.bonds.append(node)
        elif isinstance(node, BondDescriptorAtom):
            # re-direct descriptors to the one in the stoch_object
            for descriptor in stoch_object.bonding_descriptors:
                if descriptor == node.descriptor:
                    node.descriptor = descriptor
                    continue

            # if descriptor not found already in stochastic object, add it
            stoch_object.bonding_descriptors.append(node.descriptor)

    return stoch_object


def add_bigsmiles_as_stochastic_fragment(stoch_obj: StochasticObject, bigsmiles_: BigSMILES):
    """

    Parameters
    ----------
    stoch_obj
    bigsmiles_

    Returns
    -------

    """
    stoch_fragment = StochasticFragment(stoch_obj.root._get_id(StochasticFragment), stoch_obj)
    stoch_fragment.nodes = bigsmiles_.nodes
    stoch_fragment.rings = bigsmiles_.rings
    stoch_fragment.bonding_descriptors = \
        list({node.descriptor for node in bigsmiles_.nodes if isinstance(node, BondDescriptorAtom)})

    return add_stochastic_fragment(stoch_obj, stoch_fragment)


def set_new_parent(new_parent, obj):
    """ re-direct 'parent' to new bigsmiles object """
    for node in obj.nodes:
        if hasattr(node, 'parent'):
            node.parent = new_parent


__all__ = ["has_node_attr", "add_atom", "add_bond", "add_ring_by_index", "add_ring_from_atoms", "add_bonding_descriptor_atom",
           "add_bond_atom_pair", "add_bond_bonding_descriptor_pair", "open_branch", "close_branch",
           "open_stochastic_object", "open_stochastic_object_fragment", "open_stochastic_object_fragment_with_bond",
           "close_stochastic_object", "open_stochastic_fragment", "close_open_stochastic_fragment",
           "close_stochastic_fragment", "exit_construction", "get_prior"
           ]
