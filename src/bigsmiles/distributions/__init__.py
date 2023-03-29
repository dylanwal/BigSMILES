"""
Optional package

The structure of distributions class (inheritance) are:

Distribution
├── ContinuousDistribution
│    ├── LogNormal
│    ├── SchulzZimm
│    ├── Gaussian
│    ├── Uniform
│    └── Custom
├── DiscreteDistribution
     ├── Poisson
     └── FlorySchulz

"""

try:
    # optional packages need for this feature
    import scipy
    import numpy
except ImportError:
    raise ImportError("'Scipy' is an optional package and needs to be installed. "
                      "\n'pip install scipy' or 'conda install scipy' ")  # scipy will install numpy too

from bigsmiles.distributions.continuous import LogNormal, SchulzZimm, Gaussian, Uniform, CustomDistribution
from bigsmiles.distributions.discrete import FlorySchulz, Poisson
from bigsmiles.distributions.plot import plot_x_i, plot_x_i_pmd, plot_w_i, plot_x_i_cdf, plot_w_i_cdf

distributions = {
    "flory_schulz": FlorySchulz,
    "schulz-zimm": SchulzZimm,
    "log_normal": LogNormal,
    "poisson": Poisson,
    "gauss": Gaussian,
    "uniform": Uniform,
    "custom": CustomDistribution
}
