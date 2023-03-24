from __future__ import annotations
import abc
import logging
import math
from functools import wraps

from scipy.special import gamma
from scipy.optimize import minimize, root_scalar
import numpy as np

from bigsmiles.distributions.base import Distribution
import bigsmiles.distributions.utils as utils


class DistributionContinuous(Distribution, abc.ABC):
    """
    Generic class for continuous molecular weight distributions

    !!! info

        See [Distribution]() for attributes.

    """

    @abc.abstractmethod
    def _compute_x_i(self, mw_i: int | float | np.ndarray) -> int | float | np.ndarray:
        """ computes mole fraction of molecular weight as a function of $mw_i$ """

    def _compute_mw_i(self):
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
        max_y = self._compute_x_i(self.peak_mw)

        # for low MW
        if max_MW < 1000:
            self._mw_i = np.logspace(1, math.log10(30_000), 10_000)
            return

        # for med MW
        if self._compute_x_i(0) > max_y * 0.03:
            x = np.concatenate(
                (
                    np.logspace(1, math.log10(max_MW), 5000),
                    np.logspace(math.log10(max_MW + 10), math.log10(max_MW * 5), 4500),
                    np.logspace(math.log10(max_MW * 5 + 100), math.log10(max_MW * 20), 500)
                )
            )
            self._mw_i = x
            return

        # for high MW
        min_ = root_scalar(lambda x_: self._compute_x_i(x_) - (max_y * 0.03), bracket=[0, max_MW - 10])
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

        self._mw_i = x

    @utils.check_for_repeat_MW
    def _compute_N_i(self):
        n_i = self.mw_i / self.repeat_MW
        self._N_i = n_i.astype(int)

    def _compute_Mn(self):
        self._Mn, self._D = utils.compute_Mn_D_from_xi_mw_i(self.mw_i, self.x_i(self.mw_i))

    def _compute_D(self):
        self._Mn, self._D = utils.compute_Mn_D_from_xi_mw_i(self.mw_i, self.x_i(self.mw_i))

    @utils.check_for_repeat_MW
    def _compute_N(self):
        self._N = self.Mn / self._repeat_MW

    @utils.check_for_repeat_MW
    def _compute_peak_N(self):
        self._peak_N = self.peak_mw / self.repeat_MW

    def _compute_peak_mw(self):
        if self._Mn is None:
            x = np.logspace(2, 6, 100)
            index_ = np.argmax(self.x_i(x))
            x0 = np.max(x[index_])
        else:
            x0 = np.array((self.Mn,))
        scaler = 1_000 / self.x_i(x0)  # minimizer works best with large values
        max_result = minimize(lambda x_: -scaler * self.x_i(x_), x0=x0, bounds=((10, self._mw_i_max),))

        self._peak_mw = max_result.x[0]

    @utils.check_for_repeat_MW
    def _compute_x_i_pmd(self, N_i: int | float | np.ndarray) -> int | float | np.ndarray:  # noqa
        if self._x_i_pmd is None:
            x_i = self.x_i(np.linspace(0, np.max(self.N_i), np.max(self.N_i)+1) * self.repeat_MW)
            self._x_i_pmd = x_i / np.sum(x_i)

        return self._x_i_pmd[N_i]

    @wraps(Distribution.draw_mw)
    def draw_mw(self, n: int = 1, rng=None) -> int | float | np.ndarray:
        if rng is None:
            rng = self._DEFAULT_RNG

        rnd = rng.random((n,))
        x, cdf = self.x_i_cdf()

        return np.interp(rnd, cdf, x)

    @wraps(Distribution.draw_N)
    @utils.check_for_repeat_MW
    def draw_N(self, n: int = 1) -> int | np.ndarray:  # noqa
        return np.random.choice(self.N_i, n, p=self.x_i_pmd())


class LogNormal(DistributionContinuous):
    r"""
    > Log normal distributions model narrow molecular weight dispersities (D = 1-2) well. They are among the
    > most popular distribution to model living polymerization.
    >> See [Wikipida](https://en.wikipedia.org/wiki/Log-normal_distribution) for more information.

    The probability density function is:

    $$
        P(x) = \frac{1}{x\sigma_{LN}\sqrt{2\pi}}exp\Bigg(-\frac{(ln(x)-\mu_{LN})^2}{2\sigma^2}\Bigg)
    $$

    * $P(x)$: probability density function
    * $x$: random variable
    * $\mu_{LN}$: log normal mean
    * $\sigma_{LN}$: log normal standard deviation


    To put it into polymer specific terminology:

    $$
        x_i(mw_i) = \frac{1}{mw_i\sqrt{2\pi ln(Đ))}} exp\Bigg(-\frac{\Big(ln\big(\frac{mw_i}{M_n}\Big)+
        \frac{Đ}{2}\Big)^2}{2\sigma^2}\Bigg)
         ~~~~~~~~~~~
        w_i(mw_i) = \frac{1}{M_n\sqrt{2\pi ln(Đ))}} exp\Bigg(-\frac{\Big(ln\big(\frac{mw_i}{M_n}\Big)+
        \frac{Đ}{2}\Big)^2}{2\sigma^2}\Bigg)
    $$

    * $x_i(mw_i)$: mole fraction of molecular weight $mw_i$
    * $w_i(mw_i)$: weight fraction of molecular weight $mw_i$
    * $mw_i$: molecular weight of polymer 'i' units long
    * $M_n$: number average molecular weight
    * $Đ$: molecular weight dispersity

    $$\mu_{LN}=ln(M_{n})-\frac{\sigma_{LN}^2}{2} ~~~~~~~~~~~ \sigma_{LN}=\sqrt{ln(D)} ~~~~~~~~~~~ x=mw_i$$

    ---

    """
    label = "lognorm"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        Mn:
            number average molecular weight
        D:
            molecular weight dispersity
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self._Mn = Mn
        self._D = D

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"

    def _compute_x_i(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            if x == 0:
                return 0
            return 1 / (x * math.sqrt(2 * math.pi * math.log(self.D))) * \
                   np.exp(-1 * (np.log(x / self.Mn) + math.log(self.D) / 2) ** 2 / (2 * math.log(self.D)))

        mask = x > 0  # mask is to avoid ZeroDivisionError
        result = np.zeros_like(x, dtype='f')
        result[mask] = 1 / (x[mask] * math.sqrt(2 * math.pi * math.log(self.D))) * \
                       np.exp(-1 * (np.log(x[mask] / self.Mn) + math.log(self.D) / 2) ** 2 / (2 * math.log(self.D)))

        return result


class SchulzZimm(DistributionContinuous):
    r"""
    > Schulz-Zimm distributions models narrow molecular weight dispersities (D = 1-1.8) well.
    > Schulz-Zimm distributions are prone to numerical round-off errors; log normal provides similar behavior
    > without numerical issues.
    >> See [Wikipida](https://en.wikipedia.org/wiki/Schulz-Zimm_distribution) for more information.

    The probability density function is:

    $$
        x_i(mw_i) = \frac{z^{z+1}}{\Gamma(z+1)}\frac{mw_i^{z-1}}{M_n^z}exp\Big(-\frac{zmw_i}{M_n}\Big)
         ~~~~~~~~~~~
        w_i(mw_i) = \frac{z^{z+1}}{\Gamma(z+1)}\frac{mw_i^{z}}{M_n^{z+1}}exp\Big(-\frac{zmw_i}{M_n}\Big)
    $$

    * $x_i(mw_i)$: mole fraction of molecular weight $mw_i$
    * $w_i(mw_i)$: mole fraction of molecular weight $mw_i$
    * $mw_i$: molecular weight of polymer 'i' units long
    * $M_n$: number average molecular weight
    * $z$: parameter related to molecular weight dispersity

    $$ Đ = \frac{z+1}{z} $$

    """
    label = "schulz_zimm"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        Mn:
            number average molecular weight
        D:
            molecular weight dispersity
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self._Mn = Mn

        if D > 1.99:
            raise ValueError("'D' must be between (1, 1.99) to avoid numerical issues")

        self.z = 1 / (D - 1)

    def _compute_x_i(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            if x == 0:
                return 0
            return self.z ** (self.z + 1) / gamma(self.z + 1) * x ** (self.z - 1) / self.Mn ** self.z * np.exp(
                -self.z * x / self.Mn)

        mask = x > 0  # mask is to avoid Zero gamma
        pdf = np.zeros_like(x, dtype='f')
        pdf[mask] = (self.z ** (self.z + 1) / gamma(self.z + 1)) * (
                x[mask] ** (self.z - 1) / self.Mn ** self.z) * np.exp(-self.z * x[mask] / self.Mn)

        return pdf

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"


class Gaussian(DistributionContinuous):
    r"""
    > Gaussian distribution of molecular weights for geometrically distributed chain lengths.
    > Only useful for low dispersity samples (D <1.1); use log-normal instead.
    >> See [Wikipida](https://en.wikipedia.org/wiki/Normal_distribution) for more information.

    The probability density function is:

    $$
        P(x) = \frac{1}{\sigma\sqrt{2\pi}} exp\Big(-\Big(\frac{x-\mu}{4\sigma}\Big)^2\Big)
    $$

    * $P(x)$: probability density function
    * $x$: random variable
    * $\mu$: mean
    * $\sigma$: standard deviation

    $$  \mu = M_n ~~~~~~~~~~~ \sigma =M_n\sqrt{Đ-1} $$

    """
    label = "gauss"

    def __init__(self, Mn: int | float, D: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        Mn:
            number average molecular weight
        D:
            molecular weight dispersity
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self._Mn = Mn
        self._D = D
        self._std_mw = Mn * math.sqrt(D - 1)

    def _compute_x_i(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        return 1 / (self.std_mw * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((x - self.Mn) / self.std_mw) ** 2)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"

    # def prob_mw(self, mw):
    #     # if self.std < 1e-6 and abs(self.mean - mw) < 1e-6:
    #     #     return 1.0
    #     return super().prob_mw(mw)


class Uniform(DistributionContinuous):
    r"""
    > Uniform distribution of different lengths, usually useful for short chains.

    The probability density function is:

    $$
        x_i(mw_i) = \begin{cases} 0 &mw_i < low \\ 1 & low \le mw_i \le high \\ 0 &mw_i > high \end{cases}
    $$

    * $x_i(mw_i)$: mole fraction of molecular weight $mw_i$
    * $mw_i$: molecular weight of polymer 'i' units long
    * $low$: low molecular weight bound
    * $high$: high molecular weight bound

    """
    label = "uniform"

    def __init__(self, low: int | float, high: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        low:
            low molecular weight bound
        high:
            high molecular weight bound
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self.low = low
        self.high = high
        self._Mn = (high + low) / 2
        self._peak_mw = self._Mn

    def _compute_x_i(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        if isinstance(x, int) or isinstance(x, float):
            return 1 if self.low < x < self.high else 0

        low_index = np.argmin(np.abs(x - self.low))
        high_index = np.argmin(np.abs(x - self.high))
        output = np.zeros_like(x, 'f')
        output[low_index: high_index] = 1
        return output

    def _create_label(self):
        return f"|{self.label}({self.low:.0f},{self.high:.0f})|"

    def _get_mw_i(self) -> np.ndarray:
        return np.array([10, self.low - 1, self.low, self.high - 1, self.high, 10_000_000])


class CustomDistribution(DistributionContinuous):
    r"""
    Provide your own molecule weight distribution.

    """
    label = "custom"

    def __init__(self,
                 mw_i: np.ndarray,
                 x_i: np.ndarray = None,
                 w_i: np.ndarray = None,
                 repeat_MW: int | float | None = None  # noqa
                 ):
        """
        Parameters
        ----------
        mw_i:
            molecular weight of polymer 'i' units long
        x_i:
            mole fraction of molecular weight $mw_i$ (!! provide one x_i or w_i !!)
        w_i:
            weight fraction of molecular weight $mw_i$ (!! provide one x_i or w_i !!)
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)

        self._raw_mw_i = None
        self._raw_x_i = None

        self._process_inputs(mw_i, x_i, w_i)
        self._plot_range = (min(mw_i), max(mw_i))

    def _compute_x_i(self, x: int | float | np.ndarray) -> int | float | np.ndarray:
        return np.interp(x, self._raw_mw_i, self._raw_x_i, left=0, right=0)

    def _create_label(self):
        return f"|{self.label}({self.Mn:.0f},{self.D:.2f})|"

    def _process_inputs(self, mw_i: np.ndarray | None, x_i: np.ndarray | None, w_i: np.ndarray | None):
        if x_i is None and w_i is None or x_i is not None and w_i is not None:
            raise ValueError("Only provide one 'x_i' or 'w_i'. ")

        flip_flag = False
        if mw_i[0] > mw_i[-1]:
            mw_i = np.flip(mw_i)

        if x_i is not None:
            if flip_flag:
                x_i = np.flip(x_i)
            if np.any(x_i < 0):
                logging.warning("Negative values detected in 'CustomDistribution`. Results may be inaccurate.")
        elif w_i is not None:
            if flip_flag:
                w_i = np.flip(w_i)
            if np.any(w_i < 0):
                logging.warning("Negative values detected in 'CustomDistribution`. Results may be inaccurate.")

            self._Mn, self._D = utils.compute_Mn_D_from_wi_mw_i(mw_i, w_i)
            x_i = utils.wi_to_xi(mw_i, w_i, self._Mn)

        self._raw_x_i = x_i
        self._raw_mw_i = mw_i
