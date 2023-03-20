"""
Optional package


"""
from __future__ import annotations
import abc
import logging
import math
from functools import wraps

try:
    # optional packages need for this feature
    from scipy.integrate import cumulative_trapezoid
    from scipy.special import gamma, factorial, gammaln, xlogy
    from scipy.optimize import minimize, root_scalar
    import numpy as np
except ImportError:
    raise ImportError("'Scipy' is an optional package and needs to be installed. "
                      "\n'pip install scipy' or 'conda install scipy' ")  # scipy will install numpy too

import bigsmiles.errors as errors

_GLOBAL_RNG = np.random.default_rng(0)


def check_for_repeat_MW(func):
    @wraps(func)
    def _check_for_repeat_MW(*args, **kwargs):
        self = args[0]
        if hasattr(self, "repeat_MW"):
            if self.repeat_MW is None:
                # logging.warning("'repeat_MW' needs to be defined.")
                return None

        return func(*args, **kwargs)

    return _check_for_repeat_MW


class Distribution(abc.ABC):
    """
    Generic class for molecular weight distributions
    """

    def __init__(self, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        repeat_MW:
            repeat unit MW
        """
        self._repeat_MW = repeat_MW
        self._Mn: int | float | None = None
        self._D: int | float | None = None
        self._N: int | float | None = None
        self._mean: int | float | None = None
        self._std: int | float | None = None
        self._peak_mw: int | float | None = None
        self._peak_N: int | float | None = None

        self._pdf: np.ndarray | None = None
        self._mw_i: np.ndarray | None = None
        self._N_i: np.ndarray | None = None

    def __str__(self):
        return self._create_label()

    def __bool__(self):
        return self.__str__()

    @abc.abstractmethod
    def _create_label(self):
        """ create string for bigSMILES """

    @property
    def repeat_MW(self) -> float | int | None:  # noqa
        # TODO: have it figure out monomer MW from bigSMILES
        return self._repeat_MW


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
    mw_w = np.trapz(x=mw_i, y=x_i * mw_i ** 2) / mw_n
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


class DistributionContinuous(Distribution, abc.ABC):
    """
    Generic class for molecular weight distributions
    """

    def __init__(self, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        repeat_MW:
            repeat unit MW
        """
        super().__init__(repeat_MW)
        self._pdf: np.ndarray | None = None

    @abc.abstractmethod
    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        """ computes probability distribution function """

    @property
    def mw_i(self) -> np.ndarray:
        """ molecular weight array (x-axis for pdf) """
        if self._mw_i is None:
            self._mw_i = self._get_mw_i()

        return self._mw_i

    @property
    @check_for_repeat_MW
    def N_i(self) -> np.ndarray | None:  # noqa
        """ chain length array (x-axis for pdf) """
        if self._N_i is None:
            self._N_i = self.mw_i / self.repeat_MW

        return self._N_i

    @property
    def Mn(self) -> float | int:  # noqa
        """ number average MW """
        if self._Mn is None:
            self._Mn, self._D = calculate_Mn_D_from_xi(self.mw_i, self.pdf(self.mw_i))

        return self._Mn

    @property
    def Mw(self) -> float | int:  # noqa
        """ weight average MW """
        return self._D * self._Mn

    @property
    def D(self) -> float | int:  # noqa
        """ molecular weight dispersity """
        if self._D is None:
            self._Mn, self._D = calculate_Mn_D_from_xi(self.mw_i, self.pdf(self.mw_i))

        return self._D

    @property
    @check_for_repeat_MW
    def N(self) -> float | int | None:  # noqa
        """ average chain length """
        if self._N is None:
            self._N = self.Mn / self._repeat_MW

        return self._N

    @property
    def peak_mw(self) -> float | int | None:
        """ peak molecular weight """
        if self._peak_mw is None:
            self._peak_mw = self._get_peak_mw()

        return self._peak_mw

    @property
    @check_for_repeat_MW
    def peak_N(self) -> float | int | None:  # noqa
        """ peak chain length """
        return self.peak_mw / self.repeat_MW

    @property
    def mean(self) -> int | float:
        """ mean molecular weight """
        return self._Mn

    @property
    def std(self) -> int | float:
        """ standard deviation molecular weight """
        if self._std is None:
            self._std = np.sqrt(np.trapz(x=self.mw_i, y=self.pdf(self.mw_i) * (self.mw_i - self.mean) ** 2))

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
        return np.trapz(x=self.mw_i, y=self.pdf(self.mw_i) * (self.mw_i - self.mean) ** 3) / self.std ** 3

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
        return (np.trapz(x=self.mw_i, y=self.pdf(self.mw_i) * (self.mw_i - self.mean) ** 4) / self.std ** 4) - 3

    def pdf(self, mw_i: int | float | np.ndarray = None) -> np.ndarray | None:
        """ probability distribution function (number/molar) """
        if mw_i is None:
            if self._pdf is None:
                self._pdf = self._compute_pdf(self.mw_i)
            return self._pdf

        return self._compute_pdf(mw_i)

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

    def _get_peak_mw(self) -> int | float:
        # get distribution max
        if self._Mn is None:
            x = np.logspace(2, 6, 100)
            index_ = np.argmax(self.pdf(x))
            x0 = np.max(x[index_])
        else:
            x0 = np.array((self.Mn,))
        scaler = 10_000/self.pdf(x0)  # minimizer works bes with large values
        max_result = minimize(lambda x_: -scaler * self.pdf(x_), x0=x0, bounds=((10, 10_000_000),))

        return max_result.x[0]

    def _get_mw_i(self) -> np.ndarray:
        """
        attempt to determine a good MW_i

        Returns
        -------
        x:
            array of MW_i

        Notes
        -----

        * attempts to find a good x by first finding the peak max, then finding the lower bound.
        It will then create a high density of points around the peak and low density everywhere else

        """
        max_MW = self.peak_mw
        max_y = self.pdf(self.peak_mw)

        # for low MW
        if max_MW < 1000:
            return np.logspace(1, math.log10(30_000), 10_000)

        # for med MW
        if self.pdf(0) > max_y * 0.03:
            x = np.concatenate(
                (
                    np.logspace(1, math.log10(max_MW), 5000),
                    np.logspace(math.log10(max_MW + 10), math.log10(max_MW * 5), 4500),
                    np.logspace(math.log10(max_MW * 5 + 100), math.log10(max_MW * 20), 500)
                )
            )
            return x

        # for high MW
        min_ = root_scalar(lambda x_: self.pdf(x_) - (max_y * 0.03), bracket=[0, max_MW - 10])
        min_MW = max(min_.root, 500)

        low_connection = math.log10(min_MW)
        lower_connection = math.log10(min_MW + 10)
        up_connection = math.log10(max_MW + 2 * (max_MW - min_MW))
        upper_connection = math.log10(max_MW + 2 * (max_MW - min_MW) + 100)

        x = np.concatenate(
            (
                np.logspace(1, low_connection, 500),
                np.logspace(lower_connection, up_connection, 9000),
                np.logspace(upper_connection, 7, 500)
            )
        )

        return x

    def details(self) -> str:
        """ long string representation """
        text = self.__str__()
        text += f"(Mn: {self.Mn}, D: {self.D}, N: {self.N}, repeat_MW: {self.repeat_MW})"
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
            y = self.pdf(self.mw_i) / self.pdf(self.peak_mw)
        else:
            y = self.pdf(self.mw_i)

        fig.add_trace(go.Scatter(x=self.mw_i, y=y, mode='lines', name=str(self)))

        if log_scale:
            fig.update_xaxes(type="log")

        if hasattr(self, "_plot_range"):
            fig.update_xaxes(range=[math.log10(i) for i in self._plot_range])

        return fig


class LogNormal(DistributionContinuous):
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
        Mn:
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
            pre_factor = 1 / (x * math.sqrt(2 * math.pi * math.log(self.D)))
        else:
            mask = x > 0  # mask is to avoid ZeroDivisionError
            pre_factor = np.zeros_like(x)
            pre_factor[mask] = 1 / (x[mask] * math.sqrt(2 * math.pi * math.log(self.D)))

        exp_ = np.exp(-1 * (np.log(x / self.Mn) + math.log(self.D) / 2) ** 2 / (2 * math.log(self.D)))
        return pre_factor * exp_


class SchulzZimm(DistributionContinuous):
    """
    Schulz-Zimm distributions is effective at modeling SEC traces for polymers with Ɖ ranging from 1 to 1.8

    * prone to numerical round-off errors; log normal provides similar behavior without numerical issues.

    """
    label = "schulz_zimm"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """

        """
        super().__init__(repeat_MW)
        self._Mn = Mn

        if D > 1.99:
            raise ValueError("'D' must be between (1, 1.99) to avoid numerical issues")

        self.z = 1 / (D - 1)

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            if x == 0:
                return 0
            return self.z ** (self.z + 1) / gamma(self.z + 1) * x ** (self.z - 1) / self.Mn ** self.z * np.exp(
                -self.z * x / self.Mn)

        mask = x > 0  # mask is to avoid Zero gamma
        pdf = np.zeros_like(x)
        pdf[mask] = (self.z ** (self.z + 1) / gamma(self.z + 1)) * (
                x[mask] ** (self.z - 1) / self.Mn ** self.z) * np.exp(-self.z * x[mask] / self.Mn)

        return pdf

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"


class Gauss(DistributionContinuous):
    """
    Gauss distribution of molecular weights for geometrically distributed chain lengths.
    :math:`G_{\\sigma,\\mu}(N) = 1/\\sqrt{\\sigma^2 2\\pi} \\exp(-1/2 (x-\\mu^2)/\\sigma`
    where :math:`\\mu` is the mean and :math:`\\sigma^2` the variance.
    The textual representation of this distribution is: `gauss(\\mu, \\sigma)`

    Only useful for low dispersity samples (D <1.1); use log-normal instead.
    """
    label = "gauss"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Initialization of Gaussian distribution object.
        Arguments:
        ----------
        raw_text: str
             Text representation of the distribution.
             Has to start with `gauss`.
        """
        super().__init__(repeat_MW)
        self._Mn = Mn
        self._D = D
        self._std = Mn * math.sqrt(D - 1)

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        return 1 / (self.std * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((x - self.mean) / self.std) ** 2)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"

    def prob_mw(self, mw):
        # if self.std < 1e-6 and abs(self.mean - mw) < 1e-6:
        #     return 1.0
        return super().prob_mw(mw)


class Uniform(DistributionContinuous):
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
        self._Mn = (high + low) / 2
        self._peak_mw = self._Mn

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            return 1 if self.low < x < self.high else 0

        low_index = np.argmin(np.abs(x - self.low))
        high_index = np.argmin(np.abs(x - self.high))
        output = np.zeros_like(x)
        output[low_index: high_index] = 1
        return output

    def _create_label(self):
        return f"|{self.label}({self.low:.0f},{self.high:.0f})|"

    def _get_mw_i(self) -> np.ndarray:
        return np.array([10, self.low - 1, self.low, self.high - 1, self.high, 10_000_000])


def normalize_pdf(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """ make area = 1"""
    area = np.trapz(x=x, y=y)
    return y / area


class CustomDistribution(DistributionContinuous):
    """
    Distribution from data
    """
    label = "custom"

    def __init__(self, mw_i: np.ndarray, pdf: np.ndarray, repeat_MW: int | float | None = None):  # noqa
        super().__init__(repeat_MW)

        if np.any(mw_i < 0):
            logging.warning("Negative values detected in 'CustomDistribution`. Results may be inaccurate.")

        if mw_i[0] > mw_i[-1]:
            mw_i = np.flip(mw_i)
            pdf = np.flip(pdf)

        self._raw_mw_i = mw_i
        self._raw_pdf = normalize_pdf(mw_i, pdf)
        self._plot_range = (min(mw_i), max(mw_i))

    def _compute_pdf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        return np.interp(x, self._raw_mw_i, self._raw_pdf, left=0, right=0)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"


## Discrete Distribution ##
#######################################################################################################################
def calculate_Mn_D_from_ni(N_i: np.ndarray, x_i: np.ndarray, repeat_mw: int | float) -> tuple[float, float]:  # noqa
    """
    calculate Mn and D from xi (molar probability distribution function) vs MW_i

    Parameters
    ----------
    N_i:
        array containing polymer chain length 'i'
    x_i:
        array containing count fraction of polymer 'i'

    Returns
    -------
    mw_n: float
        number average molecular weight
    mw_d: float
        molecular weight dispersity

    """
    mw_n = np.sum(N_i * x_i)
    mw_w = np.sum(x_i * N_i ** 2) / mw_n
    mw_d = mw_w / mw_n
    return repeat_mw * mw_n, mw_d


class DistributionDiscrete(Distribution, abc.ABC):
    """
    Generic class for molecular weight distributions
    """

    def __init__(self, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        repeat_MW:
            repeat unit MW
        """
        super().__init__(repeat_MW)
        self._pmf: np.ndarray | None = None
        self._peak_N: int | None = None

    @abc.abstractmethod
    def _compute_pmf(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        """ computes probability mass function """

    @property
    @check_for_repeat_MW
    def mw_i(self) -> np.ndarray | None:
        """ molecular weight array (x-axis for pdf) """
        if self._mw_i is None:
            self._mw_i = self.N_i * self.repeat_MW

        return self._mw_i

    @property
    def N_i(self) -> np.ndarray:  # noqa
        """ chain length array (x-axis for pdf) """
        if self._N_i is None:
            self._N_i = np.linspace(0, 9999, 10_000, dtype="uint")  # self._get_N_i()

        return self._N_i

    @property
    @check_for_repeat_MW
    def Mn(self) -> float | int | None:  # noqa
        """ number average MW """
        if self._Mn is None:
            self._Mn, self._D = calculate_Mn_D_from_ni(self.N_i, self.pmf(self.N_i), self.repeat_MW)

        return self._Mn

    @property
    @check_for_repeat_MW
    def Mw(self) -> float | int | None:  # noqa
        """ weight average MW """
        return self._D * self._Mn

    @property
    @check_for_repeat_MW
    def D(self) -> float | int | None:  # noqa
        """ molecular weight dispersity """
        if self._D is None:
            self._Mn, self._D = calculate_Mn_D_from_ni(self.N_i, self.pmf(self.N_i), self.repeat_MW)

        return self._D

    @property
    def N(self) -> float | int | None:  # noqa
        """ average chain length """
        if self._N is None:
            self._N = self.Mn / self._repeat_MW

        return self._N

    @property
    def peak_N(self) -> float | int | None:  # noqa
        """ peak chain length """
        if self._peak_N is None:
            self._peak_N = self.N_i[np.argmax(self.pmf())]

        return self._peak_N

    @property
    @check_for_repeat_MW
    def peak_mw(self) -> float | int | None:
        """ peak molecular weight """
        return self.peak_N * self.repeat_MW

    @check_for_repeat_MW
    def pdf(self, mw_i: int | float | np.ndarray = None) -> np.ndarray | None:
        """ probability distribution function (number/molar) """
        if mw_i is None:
            return self.pmf() * self.repeat_MW

        mw_i = mw_i // self.repeat_MW
        return self._compute_pmf(mw_i) * self.repeat_MW

    @check_for_repeat_MW
    def cdf(self) -> tuple[np.ndarray, np.ndarray]:
        """ cumulative distribution function """
        return generate_cdf(self.pdf)

    def pmf(self, N_i: int | np.ndarray = None) -> np.ndarray | None:  # noqa
        """ probability distribution function (number/molar) """
        if N_i is None:
            if self._pmf is None:
                self._pmf = self._compute_pmf(self.N_i)
            return self._pmf

        return self._compute_pmf(N_i)

    @check_for_repeat_MW
    def draw_mw(self, n: int = 1):
        """
        draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        rng:
             Numpy random number generator for the generation of numbers.

        """
        return self.draw_N(n) * self.repeat_MW

    def draw_N(self, n: int = 1) -> int | np.ndarray:  # noqa
        """
        draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        rng:
             Numpy random number generator for the generation of numbers.

        """
        return np.random.choice(self.N_i, n, p=self.pmf())

    @check_for_repeat_MW
    def prob_mw(self, mw: int | float) -> float:
        """
        calculate the probability that this mw was from this distribution.

        Parameters:
        ----------
        mw:
             Integer heavy atom molecular weight.
        """
        # if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
        #     return self._distribution.cdf(mw.value) - self._distribution.cdf(mw.previous)

        return self.prob_N(mw // self.repeat_MW)

    def prob_N(self, N: int) -> float:  # noqa
        """
        calculate the probability that this mw was from this distribution.

        Parameters:
        ----------
        mw:
             Integer heavy atom molecular weight.
        """
        # if isinstance(mw, bigsmiles_gen.mol_prob.RememberAdd):
        #     return self._distribution.cdf(mw.value) - self._distribution.cdf(mw.previous)

        return self.pmf(N)

    def details(self) -> str:
        """ long string representation """
        text = self.__str__()
        text += f"(Mn: {self.Mn}, D: {self.D}, N: {self.N}, repeat_MW: {self.repeat_MW})"

        return text

    def plot_pmf(self, log_scale: bool = True, normalize: bool = True, fig=None):
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
            y = self.pmf(self.N_i) / self.pmf(self.peak_N)
        else:
            y = self.pmf(self.N_i)

        fig.add_trace(go.Scatter(x=self.N_i, y=y, mode='lines', name=str(self)))

        if log_scale:
            fig.update_xaxes(type="log")

        if hasattr(self, "_plot_range"):
            fig.update_xaxes(range=[math.log10(i) for i in self._plot_range])

        return fig

    @check_for_repeat_MW
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
            y = self.pdf(self.mw_i) / self.pdf(self.peak_mw)
        else:
            y = self.pdf(self.mw_i)

        fig.add_trace(go.Scatter(x=self.mw_i, y=y, mode='lines', name=str(self)))

        if log_scale:
            fig.update_xaxes(type="log")

        if hasattr(self, "_plot_range"):
            fig.update_xaxes(range=[math.log10(i) for i in self._plot_range])

        return fig


class FlorySchulz(DistributionDiscrete):
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
        self._N = 2 / self.a - 1

    def _compute_pmf(self, N_i: int | np.ndarray) -> int | float | np.ndarray:  # noqa
        return self.a ** 2 * N_i * (1 - self.a) ** (N_i - 1)

    @property
    def N_std(self) -> int | float:  # noqa
        return math.sqrt((2 - 2 * self.a) / self.a ** 2)

    @property
    def N_skew(self) -> int | float:  # noqa
        return (2 - self.a) / math.sqrt(2 - 2 * self.a)

    @property
    def N_kurtosis(self) -> int | float:  # noqa
        return (self.a * (self.a - 6) + 6) / (2 - 2 * self.a)

    def _create_label(self):
        return f"|{self.label}({self.conversion:.4f})|"


class Poisson(DistributionDiscrete):
    """
    Poisson's distribution are only effective at describing the distributions of anionic
    polymerizations at low MW, as there are few non-idealities. Poisson distributions tend to be too
    narrow for most living polymerizations, which has motivated the implantation of Schulz-Zimm
    distributions and log-normal distributions.

    """
    label = "poisson"

    def __init__(self, N: int | float, repeat_MW: int | float | None = None):  # noqa
        """

        Parameters:
        ----------
        Mn:

        repeat_MW:


        """
        super().__init__(repeat_MW)
        self._N = N
        # self._D = 1 + 1/N

    def _compute_pmf(self, N_i: int | np.ndarray) -> int | float | np.ndarray:
        """

        Parameters
        ----------
        N_i

        Returns
        -------

        Notes
        -----

        * computing directly causes numerical roundoff/overflow issues with factorial (e.g.;
        (self.N ** N_i * math.exp(-self.N)) / factorial(N_i, exact=True))

        """
        return np.exp(xlogy(N_i, self.N) - gammaln(N_i + 1) - self.N)

    def _create_label(self):
        return f"|{self.label}({self.N:.0f})|"


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
