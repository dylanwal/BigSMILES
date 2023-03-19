"""
Optional package


"""
from __future__ import annotations
import abc
import logging
import math

try:
    # optional packages need for this feature
    from scipy.integrate import cumulative_trapezoid
    from scipy.special import gamma
    from scipy.optimize import minimize_scalar, root_scalar
    import numpy as np
except ImportError:
    raise ImportError("'Scipy' is an optional package and needs to be installed. "
                      "\n'pip install scipy' or 'conda install scipy' ")  # scipy will install numpy too

import bigsmiles.errors as errors


_GLOBAL_RNG = np.random.default_rng()


def calculate_Mn_D_from_xi(mw_i: np.ndarray, x_i: np.ndarray) -> tuple[float, float]:  # noqa
    """
    calculate Mn and D from xi (molar probability distribution function) vs MW_i

    Parameters
    ----------
    mw_i:
        array containing molecular weight of polymer 'i' ( i is degree of polymerization, chain length)
    x_i:
        array containing mole fraction of polymer 'i' ( i is degree of polymerization, chain length)

    Returns
    -------
    mw_n: float
        number average molecular weight
    mw_d: float
        molecular weight dispersity

    """
    mw_n = np.trapz(x=mw_i, y=mw_i * x_i)
    mw_w = np.trapz(x=mw_i, y=x_i * mw_i**2) / mw_n
    mw_d = mw_w / mw_n
    return mw_n, mw_d


def generate_cdf(pdf, npts: int = 1000, x_range: tuple[float, float] = (0, 1)) -> tuple[np.ndarray, np.ndarray]:
    """
    calculate cumulative distribution function (cdf)

    Parameters
    ----------
    pdf:
        probability distribution function
    npts:
        number of points for cdf
    x_range:
        range for which the func, will be evaluated over

    Returns
    -------
    x:
        x position of cdf
    cdf:
        cdf

    """
    x = np.linspace(x_range[0], x_range[1], npts)
    y = pdf(x)

    y_norm = y / np.trapz(y, x)
    cdf = cumulative_trapezoid(y_norm, x)
    cdf = np.insert(cdf, 0, 0)

    # deal with numerical round off errors
    cdf, index_ = np.unique(cdf, return_index=True)
    x = x[index_]
    x[-1] = x_range[1]

    return x, cdf


class Distribution(abc.ABC):
    """
    Generic class for molecular weight distributions
    """

    def __init__(self, repeat_MW: int | float | None = None): # noqa
        """
        Parameters
        ----------
        repeat_MW:
            repeat unit MW
        """
        self._monomer_MW = repeat_MW
        self._Mn: int | float | None = None
        self._D: int | float | None = None
        self._N: int | float | None = None
        self._mean: int | float | None = None
        self._std: int | float | None = None
        self._max: int | float | None = None

        self._pdf: np.ndarray | None = None
        self._x: np.ndarray | None = None

    def __str__(self):
        return self._create_label()

    def __bool__(self):
        return self.__str__()

    @abc.abstractmethod
    def _create_label(self):
        """ create string for bigSMILES """

    @abc.abstractmethod
    def _compute_pdf(self, x: int | float | np.ndarray) -> np.ndarray:
        """ computes probability distribution function """

    @property
    def x(self) -> np.ndarray:
        """ either molecular weight or chain length (x-axis for pdf) """
        if self._x is None:
            self._x = self._get_x()

        return self._x

    @property
    def Mn(self) -> float | int:  # noqa
        """ number average MW """
        if self._Mn is None:
            self._Mn, self._D = calculate_Mn_D_from_xi(self._pdf(self.x), self.x)

        return self._Mn

    @property
    def Mw(self) -> float | int:  # noqa
        """ weight average MW """
        return self._D * self._Mn

    @property
    def D(self) -> float | int:  # noqa
        """ molecular weight dispersity """
        if self._D is None:
            self._Mn, self._D = calculate_Mn_D_from_xi(self._pdf(self.x), self.x)

        return self._D

    @property
    def N(self) -> float | int | None:  # noqa
        """ Chain length """
        if self.N is None:
            if self.Mn and self._monomer_MW:
                self._N = self.Mn / self._monomer_MW

        return self._N

    @property
    def monomer_MW(self) -> float | int | None:  # noqa
        # TODO: have it figure out monomer MW from bigSMILES
        return self._monomer_MW

    @property
    def max(self):
        """ max molecular weight """
        if self._max is None:
            max_result = minimize_scalar(lambda x_: -1 * self.pdf(x_), bounds=[-10, 10_000_000])
            self._max = max_result.x

        return self._max

    @property
    def mean(self) -> int | float:
        """ mean molecular weight """
        if self._mean is None:
            self._mean = np.trapz(x=self.x, y=self.x*self.pdf(self.x))

        return self._mean

    @property
    def std(self) -> int | float:
        """ standard deviation molecular weight """
        if self.std is None:
            self._std = np.sqrt(np.trapz(x=self.x, y=self.pdf(self.x)*(self.x-self.mean)**2))

        return self._std

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
        return np.trapz(x=self.x, y=self.pdf(self.x) * (self.x - self.mean)**3) / self.std**3

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
        return (np.trapz(x=self.x, y=self.pdf(self.x)*(self.x-self.mean)**4)/self.std**4)-3

    def pdf(self, x: int | float | np.ndarray = None) -> np.ndarray:
        """ probability distribution function (number/molar) """
        if x is None:
            x = self.x
        if self._pdf is not None and x == self._x:
            return self._pdf

        return self._compute_pdf(x)

    def cdf(self) -> tuple[np.ndarray, np.ndarray]:
        """ cumulative distribution function """
        return generate_cdf(self.pdf)

    def draw_mw(self, n: int = 1):
        """
        draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        rng:
             Numpy random number generator for the generation of numbers.

        """
        rnd = np.random.random((n,))
        x, cdf = self.cdf()

        return np.interp(rnd, cdf, x)

    def prob_mw(self, mw: int | float):
        """
        calculate the probability that this mw was from this distribution.

        Parameters:
        ----------
        mw:
             Integer heavy atom molecular weight.
        """
        # if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
        #     return self._distribution.cdf(mw.value) - self._distribution.cdf(mw.previous)

        return self.pdf(mw)

    def _get_x(self, mw_or_n: bool = True) -> np.ndarray:
        """
        attempt to determine a good x

        Parameters
        ----------
        mw_or_n:
            True: MW; False N

        Returns
        -------
        x:
            array of Mn or N

        Notes
        -----

        * attempts to find a good x by first finding the peak max, then finding the lower bound.
        It will then create a high density of points around the peak and low density everywhere else

        """
        if mw_or_n:
            # mw
            peak_max_range = (-100, 10_000_000)
            low_cutoff = 1000
            low_cutoff_max_mw = math.log10(30_000)
            x_min = math.log10(300)
            x_max = math.log10(10_000_000)
        else:
            # chain length
            peak_max_range = (-2, 10_000)
            low_cutoff = 50
            low_cutoff_max_mw = math.log10(400)
            x_min = 1
            x_max = math.log10(10_000)

        max_result = minimize_scalar(lambda x_: -1 * self.pdf(x_), bounds=peak_max_range)
        max_MW, max_y = max_result.x, -1 * max_result.fun  # -1 is to undo the lambda function one line above
        self._max = max_MW

        # low MW
        if max_MW < low_cutoff:
            return np.logspace(1, low_cutoff_max_mw, 10_000)

        # high MW
        min_ = root_scalar(lambda x_: self.pdf(x_) - (max_y * 0.03), bracket=[50, max_MW - 10])
        min_MW = min_.root

        low_connection = math.log10(min_MW - 10)
        lower_connection = math.log10(min_MW)
        up_connection = math.log10(max_MW + 2 * (max_MW - min_MW))
        upper_connection = math.log10(max_MW + 2 * (max_MW - min_MW) + 100)

        x = np.concatenate(
            (
                np.logspace(x_min, low_connection, 500),
                np.logspace(lower_connection, up_connection, 9000),
                np.logspace(upper_connection, x_max, 500)
                            )
            )

        return x

    def details(self) -> str:
        """ long string representation """
        text = self.__str__()
        text += f"(Mn: {self.Mn}, D: {self.D}, N: {self.N}, Monomer_MW: {self.monomer_MW})"
        text += f"\n(mean: {self.mean}, std: {self.std}, skew: {self.skew}, kurtosis: {self.kurtosis})"

        return text

    def plot_pdf(self, log_scale: bool = True, normalize: bool = True, fig=None):
        """
        for quick visualization of distribution

        Parameters
        ----------
        log_scale:
            True: x-axis log-scale; False: x-axis is linear scale
        normalize:
            True: max set 1
        fig:
            plotly figure

        Returns
        -------
        fig:
            plotly figure

        """
        try:
            import plotly.graph_objs as go
        except ImportError:
            raise ImportError("'plotly' is an optional package and needs to be installed. "
                              "\n'pip install plotly'")

        if fig is None:  # create figure is one is not provided
            fig = go.Figure()

        if normalize:
            y = self.pdf(self.x)/self.pdf(self.max)
        else:
            y = self.pdf(self.x)

        fig.add_trace(go.Scatter(x=self.x, y=y, mode='lines', name=str(self)))

        if log_scale:
            fig.update_xaxes(type="log")

        return fig


class LogNormal(Distribution):
    r"""
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

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        Mn: M
        """
        super().__init__(repeat_MW)
        self._Mn = Mn
        self._D = D

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            if x == 0:
                return 0
            pre_factor = 1 / (x * math.sqrt(2*math.pi*math.log(self.D)))
        else:
            mask = x > 0  # mask is to avoid ZeroDivisionError
            pre_factor = np.zeros_like(x)
            pre_factor[mask] = 1 / (x[mask] * math.sqrt(2*math.pi*math.log(self.D)))

        exp_ = np.exp(-1 * (np.log(x/self.Mn) + math.log(self.D)/2)**2 / (2 * math.log(self.D)))
        return pre_factor * exp_


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
    label = "flory_schulz"

    def __init__(self, conversion: float | int, repeat_MW: int | float | None = None):  # noqa
        """
        Initialization of Flory-Schulz distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `flory_schulz`.
        """
        super().__init__(repeat_MW)
        self.conversion = conversion
        self.a = 1 - self.conversion
        self._D = 1 + conversion
        self._N = 1/(1-conversion)

    def _compute_pdf(self, x: int | float | np.ndarray) -> np.ndarray:
        return self.a ** 2 * x * (1 - self.a) ** (x - 1)

    @property
    def mean(self) -> int | float:
        return 2/self.a - 1

    @property
    def std(self) -> int | float:
        return math.sqrt((2-2*self.a)/self.a**2)

    @property
    def skew(self) -> int | float:
        return (2-self.a)/math.sqrt(2-2*self.a)

    @property
    def kurtosis(self) -> int | float:
        return (self.a * (self.a - 6) + 6) / (2 - 2 * self.a)

    def _create_label(self):
        return f"|{self.label}({self.conversion:.2f})|"


class SchulzZimm(Distribution):
    """
    Schulz-Zimm distributions is effective at modeling SEC traces for polymers with Ɖ ranging from 1 to 2

    """
    label = "schulz_zimm"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """

        """
        super().__init__(repeat_MW)
        self._Mn = Mn
        self._D = D
        self.z = 1 / (D - 1)

    def _compute_pdf(self, x: int | float | np.ndarray) -> np.ndarray:
        return self.z**(self.z+1)/gamma(self.z) * x**(self.z-1)/self.Mn**self.z * np.exp(-self.z * x / self.Mn)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"


class Poisson(Distribution):
    """
   Poisson's distribution are only effective at describing the distributions of anionic
polymerizations at low MW, as there are few non-idealities. Poisson distributions tend to be too
narrow for most living polymerizations, which has motivated the implantation of Schulz-Zimm
distributions and log-normal distributions.

    """
    label = "poisson"
    default_x = np.linspace(1, 10_000, 10_000)

    def __init__(self, N: int | float, repeat_MW: int | float | None = None):  # noqa
        """

        Parameters:
        ----------
        Mn:

        repeat_MW:


        """
        super().__init__(repeat_MW)
        self._N = N
        self._D = 1 + 1/N

    def _compute_pdf(self, x: int | float | np.ndarray) -> np.ndarray:
        return (self.N**x * math.exp(-self.N)) / math.factorial(x)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f})|"


class Gauss(Distribution):
    """
    Gauss distribution of molecular weights for geometrically distributed chain lengths.
    :math:`G_{\\sigma,\\mu}(N) = 1/\\sqrt{\\sigma^2 2\\pi} \\exp(-1/2 (x-\\mu^2)/\\sigma`
    where :math:`\\mu` is the mean and :math:`\\sigma^2` the variance.
    The textual representation of this distribution is: `gauss(\\mu, \\sigma)`
    """
    label = "gauss"

    def __init__(self, mean: int | float, std: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `gauss`.
        """
        super().__init__(repeat_MW)
        self._mean = mean
        self._std = std

    def _compute_pdf(self, x: int | float | np.ndarray) -> np.ndarray:
        return 1 / (self.std * math.sqrt(2*math.pi)) * np.exp(-0.5 * ((x - self.mean)/self.std)**2)

    def _create_label(self):
        return f"|{self.label}({self.mean:.0f},{self.std:.2f})|"

    def prob_mw(self, mw):
        # if self.std < 1e-6 and abs(self.mean - mw) < 1e-6:
        #     return 1.0
        return super().prob_mw(mw)


class Uniform(Distribution):
    """
    Uniform distribution of different lengths, usually useful for short chains.
    The textual representation of this distribution is: `uniform(low, high)`
    """
    label = "uniform"

    def __init__(self, low: int | float, high: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Initialization of Uniform distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `gauss`.
        """
        super().__init__(repeat_MW)
        self.low = low
        self.high = high

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            return 1 if self.low < x < self.high else 0

        low_index = np.argmin(np.abs(x-self.low))
        high_index = np.argmin(np.abs(x-self.high))
        output = np.zeros_like(x)
        output[low_index: high_index] = 1
        return output

    def _create_label(self):
        return f"|{self.label}({self.low:.0f},{self.high:.0f})|"


class CustomDistribution(Distribution):
    """
    Distribution from data
    """
    label = "custom"

    def __init__(self, x: np.ndarray, pdf: np.ndarray, repeat_MW: int | float | None = None):  # noqa
        super().__init__(repeat_MW)

        if np.any(x < 0):
            logging.warning("Negative values detected in 'CustomDistribution`. Results may be inaccurate.")

        self._x = x
        self.pdf = pdf

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        return np.interp(x, self.x, self.pdf)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"


distributions = {
    "flory_schulz": FlorySchulz,
    "schulz-zimm": SchulzZimm,
    "log_normal": LogNormal,
    "poisson": Poisson,
    "gauss": Gauss,
    "uniform": Uniform,
    "custom": CustomDistribution
}


def get_distribution(distribution_text):
    try:
        return distributions[distribution_text]
    except KeyError:
        raise errors.DistributionError(f"Unknown distribution type {distribution_text}.")
