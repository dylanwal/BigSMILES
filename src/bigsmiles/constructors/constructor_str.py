"""

Extend existing construction functions for string input

"""

from functools import wraps

from bigsmiles.errors import BigSMILESError
from bigsmiles.constructors.tokenizer import tokenize_atom_symbol, tokenize_bonding_descriptor
from bigsmiles.data_structures.bigsmiles import StochasticObject, StochasticFragment
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
def add_atom_str(parent: has_node_attr, symbol: str) -> has_node_attr:
    return add_atom(parent, **tokenize_atom_symbol(symbol))


@in_stochastic_object
def add_bonding_descriptor_str(parent: has_node_attr, bd_symbol: str) -> has_node_attr:
    return add_bonding_descriptor(parent, *tokenize_bonding_descriptor(bd_symbol))


def add_bond_atom_pair_str(parent: has_node_attr, bond_symbol: str | None, atom_symbol: str) -> has_node_attr:
    return add_bond_atom_pair(parent, bond_symbol, **tokenize_atom_symbol(atom_symbol))


@in_stochastic_object
def add_bond_bonding_descriptor_pair_str(parent: has_node_attr, bond_symbol: str | None, bd_symbol: str)\
        -> has_node_attr:
    return add_bond_bonding_descriptor_pair(parent, bond_symbol, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_fragment_str(parent: has_node_attr, bd_symbol: str) -> StochasticFragment:
    return open_stochastic_object_fragment(parent, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_str(parent: has_node_attr, bd_symbol: str) -> StochasticObject:
    return open_stochastic_object(parent, *tokenize_bonding_descriptor(bd_symbol))


def open_stochastic_object_with_bond_str(parent: has_node_attr, bond_symbol: str | None,
                                         bd_symbol: str) -> has_node_attr:
    return open_stochastic_object_with_bond(parent, bond_symbol, *tokenize_bonding_descriptor(bd_symbol))


def close_stochastic_object_str(parent: has_node_attr, bd_symbol: str) -> has_node_attr:
    return close_stochastic_object(parent, *tokenize_bonding_descriptor(bd_symbol))


@in_stochastic_object
def open_stochastic_fragment(parent) -> has_node_attr:
    return open_stochastic_fragment(parent)


@in_stochastic_object
def close_open_stochastic_fragment_str(parent) -> has_node_attr:
    return close_open_stochastic_fragment(parent)


@in_stochastic_object
def close_stochastic_fragment_str(parent) -> has_node_attr:
    return close_stochastic_fragment(parent)
