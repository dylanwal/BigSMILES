from __future__ import annotations
import re

import bigsmiles.errors as errors
from bigsmiles.data_structures.bigsmiles import BigSMILES, Bond, Atom


rxn_split_pattern = re.compile(r"(?<!\[)(>>|>)")


def parse_reaction(text: str) -> tuple[list[BigSMILES], list[BigSMILES], list[BigSMILES]]:
    """
    Main entry point to parsing reaction SMILES.

    !!! info "Regular Expression Pattern"

        pattern: r"(?<!\[)(>>|>)"

    Parameters
    ----------
    text: str
        reaction BigSMILES string

    Returns
    -------
    reactants: list[BigSMILES]
        reactants
    agents: list[BigSMILES]
        agents
    products: list[BigSMILES]
        products
    """
    result = re.split(rxn_split_pattern, text.replace(" ", ""))

    try:
        if ">>" in result:
            return parse_no_agent_reaction(result)
        return parse_agent_reaction(result)
    except errors.BigSMILESError as e:
        raise errors.BigSMILESError(f"Parsing failed on '{text}'.", e) from e


def parse_no_agent_reaction(text_list: list[str]) -> tuple[list[BigSMILES], list[BigSMILES], list[BigSMILES]]:
    """ the pattern should be [reactants, >>, products] """
    if len(text_list) != 3:
        raise errors.TokenizeError("Invalid BigSMILES reaction. Too many or few '>' or '>>' detected. "
                                   f"\nDetected chunks: {len(text_list)} (expected chunks: 3)")
    if ">" in text_list:
        raise errors.TokenizeError("Invalid BigSMILES reaction. Reactions must follow one of these patterns: "
                                   "\n\t 'reactants >> products'\n\t 'reactants > agents > products' ")
    if text_list[1] != ">>":
        raise errors.TokenizeError("Invalid BigSMILES reaction. More than one '>>' detected.")

    reactants = process_chemical_block(text_list[0])
    products = process_chemical_block(text_list[2])

    return reactants, [], products


def parse_agent_reaction(text_list: list[str]) -> tuple[list[BigSMILES], list[BigSMILES], list[BigSMILES]]:
    """ the pattern should be [reactants, >, agents, >, products] """
    if len(text_list) != 5:
        raise errors.TokenizeError("Invalid BigSMILES reaction. Too many or few '>' or '>>' detected. "
                                   f"\nDetected chunks: {len(text_list)} (expected chunks: 5)")
    if text_list[1] != ">" and text_list[3] != ">":
        raise errors.TokenizeError("Invalid BigSMILES reaction. "
                                   "The pattern should be: 'reactants > agents > products'.")

    reactants = process_chemical_block(text_list[0])
    agents = process_chemical_block(text_list[2])
    products = process_chemical_block(text_list[4])

    return reactants, agents, products


def process_chemical_block(text: str) -> list[BigSMILES]:
    text_list = comma_split(text)

    chemicals = []
    for chunk in text_list:
        chemical = BigSMILES(chunk)
        chemical = split_chemical(chemical)
        chemicals += chemical

    return chemicals


def comma_split(text: str) -> list[str]:
    """ splits on commas only at the top level (don't hit on commas in stochastic objects) """
    flag = 0
    buffer = ""
    result = []
    for char_ in text:
        if char_ == "{":
            flag += 1
        elif char_ == "}":
            flag -= 1
        elif char_ == "," and flag == 0:
            result.append(buffer)
            buffer = ""
            continue

        buffer += char_

    if buffer:
        result.append(buffer)

    return result


## Split chemicals with '.' disconnect notation ## noqa
#######################################################################################################################
def split_chemical(bigsmiles_: BigSMILES) -> list[BigSMILES]:
    """
    Tries to split a BigSMILES at disconnects '.' if it thinks they are different molecules/polymers.

    If it can't determine, it won't do anything.

    Parameters
    ----------
    bigsmiles_: BigSMILES
        bigsmiles which may or may not be made of multiple molecules/polymers seperated by disconnect notation '.'

    Returns
    -------
    results: list[BigSMILES]
        list of BigSMILES of individual molecules/polymers

    """
    if not bigsmiles_.has_disconnect:
        return [bigsmiles_]

    split_index = determine_split_index(bigsmiles_)
    if not split_index:
        return [bigsmiles_]

    return split_bigsmiles(bigsmiles_, split_index, True)


def determine_split_index(bigsmiles_: BigSMILES) -> list[int]:
    disconnect_index = []
    for i, node in enumerate(bigsmiles_.nodes):
        if isinstance(node, Bond) and node.bond_order == 0:
            disconnect_index.append(i)

    split_index = []  # list[int]
    buffer = 0
    can_split = True
    for index_ in disconnect_index:
        buffer_fragment = 0
        # check fragment
        for node in bigsmiles_.nodes[:index_]:
            if find_bridging_bond(node, bigsmiles_.nodes[index_].id_):
                can_split = False
                break

            if isinstance(node, Atom) and node.charge != 0:
                buffer_fragment += node.charge

        buffer += buffer_fragment
        if buffer == 0 and can_split:
            split_index.append(index_)

    return split_index


def find_bridging_bond(node, cut_off_index: int) -> bool:
    """
    check bonds for an atom that has an id_ greater than cut off index.
    cut off index is the disconnect symbol index so atom index is greater than the bond bridges across a disconnect
    """
    if not hasattr(node, "bonds"):
        return False

    for bond in node.bonds:
        if bond.bond_order == 0:
            continue
        for atom in bond:
            if atom.id_ > cut_off_index:
                return True  # found a bond that bridge across disconnect

    return False


def split_bigsmiles(bigsmiles_: BigSMILES, split_index: list[int], delete_index: bool = False) -> list[BigSMILES]:
    """ Brake bigsmiles up into multiple BigSMILES. """
    if not split_index:
        return [bigsmiles_]
    split_index.sort()

    result: list[BigSMILES] = []
    for index_ in reversed(split_index):
        new_bigsmiles = BigSMILES()

        # delete bond
        if delete_index:
            bigsmiles_.nodes[index_].delete()

        # move nodes
        new_bigsmiles.nodes += bigsmiles_.nodes[index_:]
        del bigsmiles_.nodes[index_:]

        # move atoms
        atom_index = get_index(new_bigsmiles.nodes, bigsmiles_.atoms, Atom)
        if atom_index is not None:
            new_bigsmiles.root.atoms += bigsmiles_.atoms[atom_index:]
            del bigsmiles_.atoms[atom_index:]

        # move bonds
        bond_index = get_index(new_bigsmiles.nodes, bigsmiles_.bonds, Bond)
        if bond_index is not None:
            new_bigsmiles.root.bonds += bigsmiles_.bonds[bond_index:]
            del bigsmiles_.bonds[bond_index:]

        # move rings
        ring_index = get_ring_index(new_bigsmiles, bigsmiles_)
        if ring_index is not None:
            new_bigsmiles.root.rings += bigsmiles_.rings[index_:]
            del bigsmiles_.rings[index_:]

        set_new_parent(new_bigsmiles, new_bigsmiles)
        re_number_node_ids(new_bigsmiles)
        result.append(new_bigsmiles)

    result.append(bigsmiles_)
    return list(reversed(result))


def set_new_parent(new_parent, obj):
    """ re-direct 'parent' to new bigsmiles object; only need to do the first layer """
    for node in obj.nodes:
        if hasattr(node, 'parent'):
            node.parent = new_parent


def get_index(node_list, other_list, type_) -> int | None:
    for node in node_list:
        if isinstance(node, type_):
            break
    else:
        return None

    for i, obj in enumerate(other_list):
        if obj == node:
            return i

    return None


def get_ring_index(new_bigsmiles: BigSMILES, old_bigsmiles: BigSMILES) -> int | None:
    for i, ring in enumerate(old_bigsmiles.rings):
        if ring.atom1 in new_bigsmiles.atoms:
            return i

    return None


def re_number_node_ids(obj, id_: int = 0):
    """ Recursive renumbering 'id_'. """
    root = obj.root
    for node in obj.nodes:
        node.id_ = root._get_id(type(node))
        if hasattr(node, 'nodes'):
            re_number_node_ids(node, id_)
