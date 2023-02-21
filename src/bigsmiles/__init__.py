import pkg_resources
__version__ = pkg_resources.require("bigsmiles")[0].version


from bigsmiles.config import Config
from bigsmiles.errors import BigSMILESError
from bigsmiles.tokenizer import BigSMILESTokenizeError
from bigsmiles.bigsmiles import Atom, Bond, Branch, BondDescriptor, StochasticObject, StochasticFragment, BigSMILES

errors = [
    BigSMILESError,
    BigSMILESTokenizeError
]

