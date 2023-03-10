"""

Extend existing construction functions for string input

"""
from __future__ import annotations
from functools import wraps

from bigsmiles.errors import BigSMILESError
from bigsmiles.constructors.tokenizer import tokenize_atom_symbol, tokenize_bonding_descriptor
from bigsmiles.data_structures.bigsmiles import StochasticObject, StochasticFragment, BigSMILES, Branch
from bigsmiles.constructors.constructor import *


def in_stochastic_object(func):
    """ Decorator to ensure function call only occurs within stochastic object."""

    @wraps(func)
    def _in_stochastic_object(*args, **kwargs):
        if not args[0].in_stochastic_object:
            raise BigSMILESError("Must in stochastic object.")
        return func(*args, **kwargs)

    return _in_stochastic_object


## Extend existing construction functions for string input ##
#######################################################################################################################
def add_atom_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment, symbol: str) \
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
    appends an `Atom` to the end of the `parent`

    Parameters
    ----------
    parent:
        parent Atom will be added to
    symbol:
        atom symbol - will be parsed into symbol, isotope, charge, etc.

    Returns
    -------
    parent:
        returns same parent that was provided

    """
    return add_atom(parent, **tokenize_atom_symbol(symbol))


@in_stochastic_object
def add_bonding_descriptor_atom_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                                    bd_symbol: str, bond_symbol: str | None = None)\
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
    appends an `BondDescriptorAtom` to the end of the `parent`

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
        parent that the BondDescriptorAtom will be added to
    bd_symbol: str
        bonding descriptor symbol - will be parsed into symbol, index
    bond_symbol: str
        bond symbol ("", "=", "#", etc.)

    Returns
    -------
    parent:
        returns same parent that was provided

    """
    return add_bonding_descriptor_atom(parent, *tokenize_bonding_descriptor(bd_symbol), bond_symbol=bond_symbol)


def add_bond_atom_pair_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                           bond_symbol: str | None, atom_symbol: str) \
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
     appends an `Bond` followed by an 'Atom' to the end of the `parent`

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
        parent that the Bond and Atom will be added to
    bond_symbol: str
        bond symbol ("", "=", "#", etc.)
    atom_symbol: str
        atom symbol - will be parsed into symbol, isotope, charge, etc.

    Returns
    -------
    parent:
        returns same parent that was provided

    """
    return add_bond_atom_pair(parent, bond_symbol, **tokenize_atom_symbol(atom_symbol))


@in_stochastic_object
def add_bond_bonding_descriptor_pair_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                                         bond_symbol: str | None, bd_symbol: str)\
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
    appends an `Bond` followed by an 'BondDescriptorAtom' to the end of the `parent`

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
        parent that the Bond and BondDescriptorAtom will be added to
    bond_symbol: str
        bond symbol ("", "=", "#", etc.)
    bd_symbol: str
        bonding descriptor symbol - will be parsed into symbol, index

    Returns
    -------
    parent:
        returns same parent that was provided

    """
    return add_bond_bonding_descriptor_pair(parent, bond_symbol, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_fragment_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                                        bd_symbol: str) -> StochasticFragment:
    """
    opens 'StochasticObject' followed by opening 'StochasticFragment'

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
    bd_symbol: str
        bonding descriptor symbol - will be parsed into symbol, index

    Returns
    -------
    parent: StochasticFragment
        returns parent is stochastic fragment

    """
    return open_stochastic_object_fragment(parent, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                               bd_symbol: str) -> StochasticObject:
    """
    opens 'StochasticObject'

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
    bd_symbol: str
        bonding descriptor symbol - will be parsed into symbol, index

    Returns
    -------
    parent: StochasticFragment
        returns parent is stochastic object

    """
    return open_stochastic_object(parent, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_fragment_with_bond_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment,
                                                  bond_symbol: str | None, bd_symbol: str) \
        -> StochasticFragment:
    """
    opens 'StochasticObject' followed by opening 'StochasticFragment'

    Parameters
    ----------
    parent: BigSMILES, Branch, StochasticFragment
        parent of StochasticObject
    bd_symbol: str
        bonding descriptor symbol - will be parsed into symbol, index
    bond_symbol: str
        bond symbol ("", "=", "#", etc.)

    Returns
    -------
    parent: StochasticFragment
        returns parent is StochasticFragment

    """
    return open_stochastic_object_fragment_with_bond(parent, bond_symbol, *tokenize_bonding_descriptor(bd_symbol))


def close_stochastic_object_str(parent: BigSMILES | Branch | StochasticObject | StochasticFragment, bd_symbol: str,
                                bond_symbol: str | None) \
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
    closes 'StochasticObject' and append it to `parent`

    Parameters
    ----------
    parent: StochasticObject
        parent of StochasticObject
    bd_symbol: str
        bonding descriptor symbol (right side) - will be parsed into symbol, index
    bond_symbol: str
        bond symbol ("", "=", "#", etc.)

    Returns
    -------
    parent.parent: BigSMILES, Branch, StochasticFragment
        returns parent of StochasticObject
    """
    return close_stochastic_object(parent, *tokenize_bonding_descriptor(bd_symbol), bond_symbol=bond_symbol)


@in_stochastic_object
def open_stochastic_fragment(parent) -> StochasticFragment:
    """
    opens 'StochasticFragment'

    Parameters
    ----------
    parent: StochasticObject
        parent of StochasticFragment

    Returns
    -------
    parent: StochasticFragment
        returns the freshly opened StochasticFragment

    """
    return open_stochastic_fragment(parent)


@in_stochastic_object
def close_open_stochastic_fragment_str(parent) -> StochasticFragment:
    """
    closes `StochasticFragment` and appends it to 'parent.parent' and opens a new `StochasticFragment`

    Parameters
    ----------
    parent: StochasticFragment
        StochasticFragment to be closed

    Returns
    -------
    parent: StochasticFragment
        returns the freshly opened StochasticFragment

    """
    return close_open_stochastic_fragment(parent)


@in_stochastic_object
def close_stochastic_fragment_str(parent: StochasticFragment) \
        -> BigSMILES | Branch | StochasticObject | StochasticFragment:
    """
    closes `StochasticFragment` and appends it to 'parent.parent'

    Parameters
    ----------
    parent: StochasticFragment
        StochasticFragment to be closed

    Returns
    -------
    parent.parent: StochasticObject
        returns the parent.parent

    """
    return close_stochastic_fragment(parent)
