import pkg_resources
__version__ = pkg_resources.require("bigsmiles")[0].version

from bigsmiles.config import Config
from bigsmiles.errors import ERRORS
from bigsmiles.data_structures.bigsmiles import Atom, Bond, Branch, BondDescriptor, BondDescriptorAtom, \
    StochasticObject, StochasticFragment, BigSMILES
from bigsmiles.data_structures.reaction import Reaction
