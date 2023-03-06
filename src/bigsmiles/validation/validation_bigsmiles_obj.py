"""

Remaining validation done on the bigsmile object

Note:
* tokenization should take care of symbol validation
* Bigsmiles_constructor should take care of most bond and valance issues

This is mainly to catch additional leftover validation.

"""
from __future__ import annotations

from bigsmiles.data_structures.bigsmiles import BigSMILES, StochasticObject, Branch, BondDescriptor, Bond, has_node_attr
from bigsmiles.errors import ValidationError


def run_syntax_fixes(bigsmiles: BigSMILES):
    remove_unnecessary_branch_symbols(bigsmiles)
    check_ring_numbers(bigsmiles)
    re_number_node_ids(bigsmiles)


def remove_unnecessary_branch_symbols(obj: has_node_attr):
    """
    Checks if branch is notation is need; if not needed remove it
    Example: CCC(C(CC)) --> CCCCCC
    recursive function
    """
    for node in obj.nodes:
        # recursive call
        if hasattr(node, "nodes"):
            remove_unnecessary_branch_symbols(node)

    # do check
    if isinstance(obj.nodes[-1], Branch):
        # checks if branch is end of nodes; if it is then it's not needed
        branch = obj.nodes.pop()
        obj.nodes += branch.nodes


def check_ring_numbers(obj: has_node_attr):
    """ Renumber rings. """
    rings = _get_rings(obj, [])

    if not rings:
        return

    for i, ring in enumerate(rings):
        ring.ring_id = i + 1


def _get_rings(obj: has_node_attr, ring_list: list[Bond]) -> list[Bond]:
    """ Recursive get rings. """
    for node in obj.nodes:
        if hasattr(node, 'nodes'):
            _get_rings(node, ring_list)

    if hasattr(obj, 'rings'):
        ring_list += obj.rings
    return ring_list


def re_number_node_ids(obj: has_node_attr, id_: int = 0):
    """ Recursive renumbering 'id_'. """
    for node in obj.nodes:
        node.id_ = id_
        id_ += 1
        if hasattr(node, 'nodes'):
            re_number_node_ids(node, id_)


## Error checks ##
#######################################################################################################################
def run_validation(bigsmiles: BigSMILES):
    """ Main entry point for post-construction validation. """
    check_ring_closure(bigsmiles)
    check_bonding_descriptors(bigsmiles)
    check_implicit_endgroups_ends(bigsmiles)


def check_ring_closure(bigsmiles: BigSMILES):
    """ Checks to see that all rings have been closed that were started. """
    for ring in bigsmiles.rings:
        if ring.atom2 is None:
            raise ValidationError(f"Ring opened, but not closed. (ring id: {ring.ring_id})")


def check_bonding_descriptors(bigsmiles: BigSMILES | StochasticObject | Branch):
    """

    rules:
    if [$] must be 2 or more
    if [>] there must be [<]
    if [<] there must be [>]
    if end is implicit; there must be single bonding units
    if ends are explict they must both have a matching pair in the stochastic fragments



    Parameters
    ----------
    bigsmiles

    Returns
    -------

    """
    for node in bigsmiles.nodes:
        if isinstance(node, StochasticObject):
            _do_check_bonding_descriptors(node)

        # check nested stochastic objects
        if hasattr(node, "nodes"):
            check_bonding_descriptors(node)


def _do_check_bonding_descriptors(stoch_obj: StochasticObject):
    for bd in stoch_obj.bonding_descriptors:
        a = bd.bond_symbol  # runs bond order check

        # if [$], [$0], [$1], ... must be 2 or more
        if bd.descriptor == "$":
            if len(bd.instances) < 1:
                raise ValidationError("[$] type bonding descriptors require more than one instances in a string.")

        # if [>] there must be [<] and if [<] there must be [>]
        if bd.descriptor in ("<", ">"):
            bd_pair = _find_complementing_bonding_descriptor(stoch_obj, bd)
            if bd_pair is None:
                raise ValidationError(f'{bd} complementary partner not found.')


def _find_complementing_bonding_descriptor(stoch_obj: StochasticObject, bond_descr: BondDescriptor) \
        -> BondDescriptor | None:
    for bd in stoch_obj.bonding_descriptors:
        if bond_descr.is_pair(bd):
            return bd

    return None


def check_implicit_endgroups_ends(obj: has_node_attr, parent_obj: has_node_attr = None):
    # if end is implicit; there must be single bonding units
    for i, node in enumerate(obj.nodes):
        if isinstance(node, StochasticObject):
            if node.bd_left.implicit:
                if i != 0:
                    # nothing allowed to the left
                    raise ValidationError("With a the left end group implicit, "
                                          "there should be nothing to the left of the stochastic object.")
                if parent_obj is not None:
                    # if it isn't BigSMILES it can't have a left implicit end group
                    raise ValidationError("Implicit left end group not allowed within interior.")

            if node.bd_right.implicit:

                if i != len(obj.nodes) - 1:
                    # nothing allowed to the right
                    raise ValidationError("With a implicit right end group, there should not be anything to the "
                                          "right of it.")

        if hasattr(node, "nodes"):
            check_implicit_endgroups_ends(node, obj)
