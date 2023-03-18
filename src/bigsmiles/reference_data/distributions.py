"""
Optional

"""
from __future__ import annotations

import abc
import math
from ast import literal_eval as make_tuple


try:
    # optional packages need for this feature
    from scipy import stats  # noqa
    import numpy as np  # noqa
except ImportError:
    raise ImportError("'Scipy' is an optional package and needs to be installed. "
                      "\n'pip install scipy' or 'conda install scipy' ")  # scipy will install numpy too

import bigsmiles.errors as errors


_GLOBAL_RNG = np.random.default_rng()


def cal_Mn_D_from_xi(mw_i: np.ndarray, x_i: np.ndarray) -> tuple[float, float]:
    """ calculate Mn and D from xi vs MW data (MW goes low to high) """
    mw_n = np.sum(x_i * mw_i)
    mw_w = np.sum(x_i * mw_i**2) / mw_n
    mw_d = mw_w / mw_n
    return mw_n, mw_d


class Distribution(abc.ABC):
    """
    Generic class to generate molecular weight distribution.
    """
    def __init__(self, distribution):
        """
        Initialize the generic distribution.
        """

        self._Mn = None
        self._D = None
        self._monomer_MW = None
        self._distribution = distribution

    def __str__(self):
        return self._create_label()

    def __bool__(self):
        return self.__str__()

    @abc.abstractmethod
    def _create_label(self):
        """ Create string for bigSMILES. """

    @property
    def Mn(self) -> float | int:  # noqa
        """ number average MW """
        return self._Mn

    @property
    def Mw(self) -> float | int:  # noqa
        """ weight average MW """
        return self._D * self._Mn

    @property
    def D(self) -> float | int:  # noqa
        """ molecular weight dispersity """
        return self._D

    @property
    def N(self) -> float | int | None:  # noqa
        """ Chain length """
        if self._monomer_MW is None:
            return None

        return self._Mn / self._monomer_MW

    @property
    def monomer_MW(self) -> float | int | None:  # noqa
        # TODO: have it figure out monomer MW from bigSMILES
        return self._monomer_MW

    @property
    def mean(self) -> int | float:
        return np.trapz(x=self.mw_i, y=self.mw_i*self.xi)

    @property
    def std(self) -> int | float:
        return np.sqrt(np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 2))

    @property
    def skew(self) -> int | float:
        """
        skew (3 moment)

        Returns
        -------

        Notes
        -----

        * Symmetric: Values between -0.5 to 0.5
        * Moderated Skewed: Values between -1 and -0.5 or between 0.5 and 1
        * Highly Skewed: Values less than -1 or greater than 1

        """
        return np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 3) / self.mw_std ** 3

    @property
    def kurtosis(self) -> int | float:
        """
        kurtosis (4 moment)

        Returns
        -------

        Notes
        -----

        * Mesokurtic: normal distribution
        * Leptokurtic: kurtosis > 3; the distribution has fatter tails and a sharper peak
        * Platykurtic: kurtosis < -3; the distribution has a lower and wider peak and thinner tails

        """
        return (np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 4) / self.mw_std ** 4)- 3

    def pdf(self, x: np.ndarray | None = None, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """ probability density function """
        if x is None:
            x = np.logspace(2, 6, 1000)  # 100 to 1_000_000

        return x, self._distribution.pdf(x, **kwargs)

    def cdf(self, x: np.ndarray | None = None, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """ cumulative distribution function """
        if x is None:
            x = np.linspace(0, 1, 1000)

        return x, self._distribution.cdf(x, **kwargs)

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
        # if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
        #     return self._distribution.cdf(mw.value) - self._distribution.cdf(mw.previous)

        return self._distribution.pdf(mw)

    def plot_pdf(self, log_scale: bool = True):
        try:
            import plotly.graph_objs as go
        except ImportError:
            raise ImportError("'plotly' is an optional package and needs to be installed. "
                              "\n'pip install plotly'")

        x, pdf = self.pdf()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=pdf, mode='lines'))

        if log_scale:
            fig.update_xaxes(type="log")

        return fig


class LogNormal(Distribution):
    """
    Similar to Schulz-Zimm and commonly used to fit experimental molecular weight distributions


    Notes
    -----
    The probability density function for `lognorm` is:
    .. math::
        f(x, s) = \frac{1}{s x \sqrt{2\pi}}
                  \exp\left(-\frac{\log^2(x)}{2s^2}\right)
    for :math:`x > 0`, :math:`s > 0`.

    """
    label = "lognorm"

    class log_normal_dist(stats.rv_continuous):
        """ Flory Schulz distribution """
        def __init__(self, Mn, D):
            super().__init__(momtype=0, a=0, b=1_000_000)
            self.D = D
            self.Mn = Mn

        def _pdf(self, x_):
            pre_factor = 1 / (x_ * math.sqrt(2*math.pi*math.log(self.D)))
            exp_ = np.exp(-1 * (np.log(x_/self.Mn) + math.log(self.D)/2)**2 / (2 * math.log(self.D)))
            return pre_factor * exp_

    def __init__(self, Mn: int | float, D: int | float):  # noqa
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        Mn: M
        """
        distribution = self.log_normal_dist(Mn, D)
        super().__init__(distribution)
        self._Mn = Mn
        self._D = D

    def _create_label(self):
        return f"|{self.label}({self.Mn}, {self.D})|"


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


class CustomDistribution(Distribution):
    """
    Distribution from data
    """

    def __init__(self, pdf: np.ndarray):
        Mn = 100
        D = 100
        super().__init__(Mn, D)

    def _get_dis(self):
        kernel = stats.gaussian_kde(data)
        class rv(stats.rv_continuous):
            def _cdf(self, x):
                return kernel.integrate_box_1d(-numpy.Inf, x)
        return rv()


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
