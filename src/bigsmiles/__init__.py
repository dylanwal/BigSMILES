import pkg_resources  # part of setuptools
__version__ = pkg_resources.require("bigsmiles")[0].version

from bigsmiles.config import Config
from bigsmiles.tokenizer import BigSMILESTokenizeError
from bigsmiles.bigsmiles_constructor import BigSMILESError
from bigsmiles.bigsmiles import Atom, Bond, Branch, BondDescriptor, StochasticObject, StochasticFragment, BigSMILES, \
    AtomChirality, BondDescriptorTypes, BondType

errors = [
    BigSMILESError,
    BigSMILESTokenizeError
]

enums = [
    BondType,
    BondDescriptorTypes,
    AtomChirality
]
