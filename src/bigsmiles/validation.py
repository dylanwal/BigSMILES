
"""

Remaining validation done on the bigsmile object

Note:
* tokenization should take care of symbol validation
* Bigsmiles_constructor should take care of most bond and valance issues

This is mainly to catch additional leftover validation.

"""

from bigsmiles.errors import BigSMILESError
from bigsmiles.bigsmiles import BigSMILES, BondDescriptorTypes, StochasticObject, StochasticFragment


def run_validation(bigsmiles: BigSMILES):
    """ Main entry point for validation. """
    check_ring_closure(bigsmiles)
    check_bonding_descriptors(bigsmiles)


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

        # if [>] there must be [<]
        if bd.type_ is BondDescriptorTypes.Left:
            pass

        # if [<] there must be [>]