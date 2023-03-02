import re
from functools import wraps

import bigsmiles.chemical_data as chemical_data
from bigsmiles.errors import BigSMILESError
from bigsmiles.constructors.tokenizer import atom_pattern
from bigsmiles.data_structures.bigsmiles import StochasticObject, StochasticFragment
from bigsmiles.constructors.constructor import *


ATOM_PATTERN = re.compile(atom_pattern)


def atom_symbol_to_attributes(symbol: str) -> dict:
    if symbol in chemical_data.elements_aromatic:
        return {"element": symbol}

    try:
        results = ATOM_PATTERN.match(symbol).groupdict()

        if results["hydrogens"] is None:
            results["hydrogens"] = 0
        else:
            if results["hydrogens"][-1].isdigit():
                results["hydrogens"] = int(results["hydrogens"][-1])
            else:
                results["hydrogens"] = 1

        if results["charge"] is None:
            results["charge"] = 0
        else:
            results["charge"] = int(results["charge"])

        return results
    except AttributeError:
        raise BigSMILESError(f"Issue parsing atom: {symbol}")


DEFAULT_BONDING_DESCRIPTOR_INDEX = 1


def process_bonding_descriptor_symbol(symbol: str) -> tuple[str, int]:
    symbol = symbol.replace("[", "").replace("]", "")
    if not symbol:
        return symbol, DEFAULT_BONDING_DESCRIPTOR_INDEX

    if symbol[-1].isdigit():
        return symbol[0], int(symbol[-1])

    return symbol, DEFAULT_BONDING_DESCRIPTOR_INDEX


def in_stochastic_object(func):
    """ Decorator to ensure function call only occurs within stochastic object."""

    @wraps(func)
    def _in_stochastic_object(*args, **kwargs):
        if not args[0].in_stochastic_object:
            raise BigSMILESError("Must in stochastic object.")
        return func(*args, **kwargs)

    return _in_stochastic_object


## Extend existing construction functions for string inputs ##
#######################################################################################################################
def add_atom_str(parent: has_node_attr, symbol: str) -> has_node_attr:
    return add_atom(parent, **atom_symbol_to_attributes(symbol))


@in_stochastic_object
def add_bonding_descriptor_str(parent: has_node_attr, bd_symbol: str) -> has_node_attr:
    return add_bonding_descriptor(parent, *process_bonding_descriptor_symbol(bd_symbol))


def add_bond_atom_pair_str(parent: has_node_attr, bond_symbol: str,
                           atom_symbol: str) -> has_node_attr:
    return add_bond_atom_pair(parent, bond_symbol, **atom_symbol_to_attributes(atom_symbol))


@in_stochastic_object
def add_bond_bonding_descriptor_pair_str(parent: has_node_attr, bond_symbol: str,
                                         bd_symbol: str) -> has_node_attr:
    return add_bond_bonding_descriptor_pair(parent, bond_symbol, *process_bonding_descriptor_symbol(bd_symbol))


def open_stochastic_object_fragment_str(parent: has_node_attr, bd_symbol: str) -> StochasticFragment:
    return open_stochastic_object_fragment(parent, *process_bonding_descriptor_symbol(bd_symbol))


def open_stochastic_object_str(parent: has_node_attr, bd_symbol: str) -> StochasticObject:
    return open_stochastic_object(parent, *process_bonding_descriptor_symbol(bd_symbol))


def open_stochastic_object_with_bond_str(parent: has_node_attr, bond_symbol: str,
                                         bd_symbol: str) -> has_node_attr:
    return open_stochastic_object_with_bond(parent, bond_symbol, *process_bonding_descriptor_symbol(bd_symbol))


def close_stochastic_object_str(parent: has_node_attr, bd_symbol: str) -> has_node_attr:
    return close_stochastic_object(parent, *process_bonding_descriptor_symbol(bd_symbol))


@in_stochastic_object
def open_stochastic_fragment(parent) -> has_node_attr:
    return open_stochastic_fragment(parent)


@in_stochastic_object
def close_open_stochastic_fragment_str(parent) -> has_node_attr:
    return close_open_stochastic_fragment(parent)


@in_stochastic_object
def close_stochastic_fragment_str(parent) -> has_node_attr:
    return close_stochastic_fragment(parent)
