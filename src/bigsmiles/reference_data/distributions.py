"""
Optional

"""

import abc
import math
from ast import literal_eval as make_tuple


try:
    # optional packages need for this feature
    from scipy import stats  # noqa
    import numpy as np  # noqa
except ImportError:
    raise ImportError("'Scipy' is an optional package and needs to be installed. "
                      "\n'pip install scipy' or 'conda install scipy' ")

import bigsmiles.errors as errors


_GLOBAL_RNG = np.random.default_rng()


class Distribution(abc.ABC):
    """
    Generic class to generate molecular weight distribution.
    """
    def __init__(self, Mn: int | float, D: int | float = 1.2):
        """
        Initialize the generic distribution.
        """
        self._Mn = Mn
        self._D = D
        self._distribution = None
        self._monomer_MW = None
        self._create_distribution()

    @property
    def Mn(self) -> float | int:
        """ number average MW """
        return self._Mn

    @property
    def Mw(self) -> float | int:
        """ weight average MW """
        return self._D * self._Mn

    @property
    def D(self) -> float | int:
        """ molecular weight dispersity """
        return self._D

    @property
    def N(self) -> float | int | None:
        """ Chain length """
        if self._monomer_MW is None:
            return None

        return self._Mn / self._monomer_MW

    @property
    def monomer_MW(self) -> float | int:
        # TODO: have it figure out monomer MW from bigSMILES
        return self.monomer_MW

    @property
    def skew(self) -> int | float:
        return self._distribution.stats('s')

    @property
    def kurtosis(self) -> int | float:
        return self._distribution.stats('K')

    @property
    def pdf(self, x) -> np.ndarray:
        return self._distribution.pdf()

    @property
    def cdf(self, x) -> np.ndarray:
        return self._distribution.cdf()

    @abc.abstractmethod
    def _create_distribution(self):
        """ create distribution"""

    def draw_mw(self, rng: np.random.Generator = None):
        """
        Draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        rng:
             Numpy random number generator for the generation of numbers.

        """
        if rng is None:
            rng = _GLOBAL_RNG
        return self._distribution.rvs(random_state=rng)

    def prob_mw(self, mw: int | float):
        """
        Calculate the probability that this mw was from this distribution.

        Parameters:
        ----------
        mw:
             Integer heavy atom molecular weight.
        """
        if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
            return self._distribution.cdf(mw.value) - self._distribution.cdf(mw.previous)

        return self._distribution.pdf(mw)

    def plot(self):
        try:
            import plotly.graph_objs as go
        except ImportError:
            raise ImportError("'plotly' is an optional package and needs to be installed. "
                              "\n'pip install plotly'")

        return go.Scatter(x=self.pdf[:, 0], y=self.pdf[:, 1])


class FlorySchulz(Distribution):
    """
    Flory-Schulz distribution of molecular weights for geometrically distributed chain lengths.
    :math:`W_a(N) = a^2 N (1-a)^M`
    where :math:`0<a<1` is the experimentally determined constant of remaining monomers and
    :math:`k` is the chain length.
    The textual representation of this distribution is: `flory_schulz(a)`

    ideal step-growth polymerization process, One-Parameter Functions

    [wiki](https://en.wikipedia.org/wiki/Flory%E2%80%93Schulz_distribution)
    """

    class flory_schulz_gen(stats.rv_discrete):
        """ Flory Schulz distribution """

        def _pmf(self, k, a):
            return a ** 2 * k * (1 - a) ** (k - 1)

    def __init__(self, conversion: float | int = None, N: int = None):
        """
        Initialization of Flory-Schulz distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `flory_schulz`.
        """
        if N is None:
            N = 1/(1-conversion)
        else:
            conversion = 1-1/N

        D = 1 + conversion

        super().__init__()

    def _create_distribution(self):
        self._distribution = self.flory_schulz_gen(name="Flory-Schulz")

    def generate_string(self, extension):
        if extension:
            return f"|flory_schulz({self._a})|"
        return ""

    @property
    def generable(self):
        return True

    def draw_mw(self, rng=None):
        if rng is None:
            rng = _GLOBAL_RNG
        return self._distribution.rvs(a=self._a, random_state=rng)

    def prob_mw(self, mw):
        if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
            return self._distribution.cdf(mw.value, a=self._a) - self._distribution.cdf(
                mw.previous, a=self._a
            )
        return self._distribution.pmf(int(mw), a=self._a)


class SchulzZimm(Distribution):
    """
    Schulz-Zimm distributions is effective at modeling SEC traces for polymers with Ɖ ranging from 1 to 2

    """

    def __init__(self, Mn: int | float, D: float):
        """

        """
        super().__init__(Mn, D)

    def _create_distribution(self):
        self.std = math.sqrt(math.log(self.D))
        self.mean = math.log(self.Mn) - self.std**2/2
        self._distribution = stats.lognorm(s=self.std, loc=self.mean)

    def generate_string(self, extension):
        if extension:
            return f"|schulz_zimm({self.Mn}, {self.D})|"
        return ""


class LogNormal(Distribution):
    """
    Similar to Schulz-Zimm and commonly used to fit experimental molecular weight distributions

    """

    def __init__(self, Mn: int | float, D: float):
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        Mn: M
        """
        super().__init__(Mn, D)

    def _create_distribution(self):
        self.std = math.sqrt(math.log(self.D))
        self.mean = math.log(self.Mn) - self.std**2/2
        self._distribution = stats.lognorm(s=self.std, loc=self.mean)

    def generate_string(self, extension):
        if extension:
            return f"|lognorm({self.Mn}, {self.D})|"
        return ""


class Poisson(Distribution):
    """
   Poisson's distribution are only effective at describing the distributions of anionic
polymerizations at low MW, as there are few non-idealities. Poisson distributions tend to be too
narrow for most living polymerizations, which has motivated the implantation of Schulz-Zimm
distributions and log-normal distributions.

    """

    def __init__(self, Mn: int | float, monomer_MW: int | float):
        """

        Parameters:
        ----------
        Mn:

        monomer_MW:


        """
        D = 1 + monomer_MW/Mn
        super().__init__(Mn, D)

    def _create_distribution(self):
        self._distribution = stats.poisson(mu=self.N)

    def generate_string(self, extension):
        if extension:
            return f"|poisson({self.Mn}, {self.D})|"
        return ""


class Gauss(Distribution):
    """
    Gauss distribution of molecular weights for geometrically distributed chain lengths.
    :math:`G_{\\sigma,\\mu}(N) = 1/\\sqrt{\\sigma^2 2\\pi} \\exp(-1/2 (x-\\mu^2)/\\sigma`
    where :math:`\\mu` is the mean and :math:`\\sigma^2` the variance.
    The textual representation of this distribution is: `gauss(\\mu, \\sigma)`
    """

    def __init__(self, raw_text):
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `gauss`.
        """
        super().__init__(raw_text)

        if not self._raw_text.startswith("gauss"):
            raise RuntimeError(
                f"Attempt to initialize Gaussian distribution from text '{raw_text}' that does not start with 'gauss'"
            )

        self._mu, self._sigma = make_tuple(self._raw_text[len("gauss"):])
        self._mu = float(self._mu)
        self._sigma = float(self._sigma)
        self._distribution = stats.norm(loc=self._mu, scale=self._sigma)

    def generate_string(self, extension):
        if extension:
            return f"|gauss({self._mu}, {self._sigma})|"
        return ""

    @property
    def generable(self):
        return True

    def prob_mw(self, mw):
        if self._sigma < 1e-6 and abs(self._mu - mw) < 1e-6:
            return 1.0
        return super().prob_mw(mw)


class Uniform(Distribution):
    """
    Uniform distribution of different lengths, usually useful for short chains.
    The textual representation of this distribution is: `uniform(low, high)`
    """

    def __init__(self, raw_text):
        """
        Initialization of Uniform distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `gauss`.
        """
        super().__init__(raw_text)

        if not self._raw_text.startswith("uniform"):
            raise RuntimeError(
                f"Attempt to initialize Uniform distribution from text '{raw_text}' that does not start with 'uniform'"
            )

        self._low, self._high = make_tuple(self._raw_text[len("uniform"):])
        self._low = int(self._low)
        self._high = int(self._high)
        self._distribution = stats.uniform(loc=self._low, scale=(self._high - self._low))

    def generate_string(self, extension):
        if extension:
            return f"|uniform({self._low}, {self._high})|"
        return ""

    @property
    def generable(self):
        return True


distributions = {
    "flory_schulz": FlorySchulz,
    "schulz-zimm": SchulzZimm,
    "log_normal": LogNormal,
    "poisson": Poisson,
    "gauss": Gauss,
    "uniform": Uniform,
}


def get_distribution(distribution_text):
    try:
        return distributions[distribution_text]
    except KeyError:
        raise errors.DistributionError(f"Unknown distribution type {distribution_text}.")
