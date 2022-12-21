
"""

Remaining validation done on the bigsmile object

Note:
* tokenization should take care of symbol validation
* Bigsmiles_constructor should take care of most bond and valance issues

This is mainly to catch additional leftover validation.

"""

from bigsmiles.errors import BigSMILESError
from bigsmiles.bigsmiles import BigSMILES, BondDescriptorTypes, StochasticObject, StochasticFragment, \
    BondDescriptor, Branch


def run_validation(bigsmiles: BigSMILES):
    """ Main entry point for validation. """
    check_ring_closure(bigsmiles)
    check_bonding_descriptors(bigsmiles)
    check_implicit_endgroups_ends(bigsmiles)


def check_ring_closure(bigsmiles: BigSMILES):
    """ Checks to see that all rings have been closed that were started. """
    for ring in bigsmiles.rings:
        if ring.atom2 is None:
            raise BigSMILESError(f"Ring opened, but not closed. (ring id: {ring.ring_id})")


def check_bonding_descriptors(bigsmiles: BigSMILES | StochasticObject):
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
            do_check_bonding_descriptors(node)

            # check nested stochastic objects
            check_bonding_descriptors(node)


def do_check_bonding_descriptors(stoch_obj: StochasticObject):

    for bd in stoch_obj.bonding_descriptors:
        # if [$], [$0], [$1], ... must be 2 or more
        if bd.type_ is BondDescriptorTypes.Dollar:
            if len(bd.instances) <= 1:
                raise BigSMILESError("[$] type bonding descriptors require more than one instances in a string.")

        # if [>] there must be [<] and if [<] there must be [>]
        if bd.type_ in (BondDescriptorTypes.Left, BondDescriptorTypes.Right):
            bd_pair = find_complementing_bonding_descriptor(stoch_obj, bd)
            if bd_pair is None:
                raise BigSMILESError(f'{bd} complementary partner not found.')


def find_complementing_bonding_descriptor(stoch_obj: StochasticObject, bond_descr: BondDescriptor) \
        -> BondDescriptor | None:
    for bd in stoch_obj.bonding_descriptors:
        if bond_descr.is_pair(bd):
            return bd

    return None


def check_implicit_endgroups_ends(obj: BigSMILES | StochasticObject | StochasticFragment | Branch,
                             parent_obj: BigSMILES | StochasticObject | StochasticFragment | Branch = None):
    # if end is implicit; there must be single bonding units
    for i, node in enumerate(obj.nodes):
        if isinstance(node, StochasticObject):
            if node.end_group_left.descriptor.type_ is BondDescriptorTypes.Implicit:
                if i != 0:
                    # nothing allowed to the left
                    raise BigSMILESError("With a the left end group implicit, "
                                         "there should be nothing to the left of the stochastic object.")
                if parent_obj is not None:
                    # if it isn't BigSMILES it can't have a left implicit end group
                    raise BigSMILESError("Implicit left end group not allowed within interior.")

            if node.end_group_right.descriptor.type_ is BondDescriptorTypes.Implicit:

                if i != len(obj.nodes)-1:
                    # nothing allowed to the right
                    raise BigSMILESError("With a implicit right end group, there should not be anything to the "
                                         "right of it.")

        if hasattr(node, "nodes"):
            check_implicit_endgroups_ends(node, obj)

