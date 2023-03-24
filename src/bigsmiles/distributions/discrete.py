"""
Discrete Distribution

"""
import abc
from functools import wraps

import numpy as np
from scipy.special import gammaln, xlogy

from bigsmiles.distributions.base import Distribution
import bigsmiles.distributions.utils as utils


class DistributionDiscrete(Distribution, abc.ABC):
    """
    Generic class for molecular weight distributions

    !!! info

        See [Distribution]() for attributes.

    """

    @abc.abstractmethod
    def _compute_x_i_pmd(self, N_i: int | np.ndarray) -> int | np.ndarray:  # noqa
        """ computes probability mass function """

    @utils.check_for_repeat_MW
    def _compute_mw_i(self) -> np.ndarray:
        return self.N_i * self.repeat_MW

    def _compute_N_i(self) -> np.ndarray:
        return np.linspace(0, self._N_i_max - 1, self._N_i_max, dtype="uint")

    @utils.check_for_repeat_MW
    def _compute_Mn(self):
        self._Mn = self.repeat_MW * self.N

    def _compute_D(self):
        self._D = np.sum(self.x_i(self.N_i) * self.N_i ** 2) / self.N

    def _compute_N(self):
        self._N = np.sum(self.N_i * self.x_i(self.N_i))

    def _compute_peak_N(self):  # noqa
        self._peak_N = self.N_i[np.argmax(self.x_i_pmd())]

    @utils.check_for_repeat_MW
    def _compute_peak_mw(self):
        """ peak molecular weight """
        self._peak_mw = self.peak_N * self.repeat_MW

    @utils.check_for_repeat_MW
    def _compute_std_mw(self):
        super()._compute_std_mw()

    @utils.check_for_repeat_MW
    def _compute_skew_mw(self):
        super()._compute_skew_mw()

    @utils.check_for_repeat_MW
    def _compute_kurtosis_mw(self):
        super()._compute_kurtosis_mw()

    @utils.check_for_repeat_MW
    def _compute_x_i(self, mw_i: int | float | np.ndarray) -> int | float | np.ndarray:
        if self._x_i is None:
            self._x_i_pmd = utils.normalize_distribution(self.N_i, self.x_i_pmd(self.N_i))

        return np.interp(mw_i, self._x_i_pmd, self.mw_i, 0, 0)

    @wraps(Distribution.draw_mw)
    @utils.check_for_repeat_MW
    def draw_mw(self, n: int = 1):
        return self.draw_N(n) * self.repeat_MW

    @wraps(Distribution.draw_N)
    def draw_N(self, n: int = 1) -> int | np.ndarray:  # noqa
        return np.random.choice(self.N_i, n, p=self.x_i_pmd())


class FlorySchulz(DistributionDiscrete):
    r"""
    > Flory-Schulz distribution model ideal step-growth polymerization.
    >> See [Wikipida](https://en.wikipedia.org/wiki/Flory%E2%80%93Schulz_distribution) for more information.

    The probability mass function is:

    $$
        x_i(N_i) = a^2 N_i (1-a)^(N_i -1)
    $$

    * $x_i(N_i)$: mole fraction of chain length $N_i$
    * $N_i$: chain length of polymer 'i' units long
    * $a$: fraction of remaining monomer ($a \epsilon [0,1]$) $a=1-conversion$

    """
    label = "flory_schulz"

    def __init__(self, conversion: float | int, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        conversion:
            monomer conversion
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self.conversion = conversion
        self.a = 1 - self.conversion
        self._D = 1 + conversion
        self._N = 2 / self.a - 1

    def _compute_x_i_pmd(self, N_i: int | np.ndarray) -> int | float | np.ndarray:  # noqa
        return self.a ** 2 * N_i * (1 - self.a) ** (N_i - 1)

    # @property
    # def N_std(self) -> int | float:  # noqa
    #     return math.sqrt((2 - 2 * self.a) / self.a ** 2)
    #
    # @property
    # def N_skew(self) -> int | float:  # noqa
    #     return (2 - self.a) / math.sqrt(2 - 2 * self.a)
    #
    # @property
    # def N_kurtosis(self) -> int | float:  # noqa
    #     return (self.a * (self.a - 6) + 6) / (2 - 2 * self.a)

    def _create_label(self):
        return f"|{self.label}({self.conversion:.4f})|"


class Poisson(DistributionDiscrete):
    r"""
    > Poisson's distribution model ideal living polymerization.
    > Poisson distributions tend to be to narrow for most living polymerization's (except for anionic polymerization)
    > and log-normal distributions provide better modeling (for ATRP, RAFT, NMP, etc.).
    >> See [Wikipida](https://en.wikipedia.org/wiki/Poisson_distribution) for more information.

    The probability mass function is:

    $$
        x_i(N_i) = \frac{N^N_i exp(-N)}{N_i!} \approx exp(N_i ln(N) - ln(\Gamma + 1)-N)
    $$

    * $x_i(N_i)$: mole fraction of chain length $N_i$
    * $N_i$: chain length of polymer 'i' units long
    * $N$: average chain length

    $$ Đ = 1 +\frac{1}{N} $$

    Note: approximation used to avoid roundoff/overflow issues with factorial.
    See for more information on the approximation:

    * DOI: https://doi.org/10.1090/S0025-5718-63-99188-9
    * William H Press, Saul A Teukolsky, William T Vetterling, and Brian P Flannery.
    Numerical recipes: The art of scientific computing. Cambridge University Press, third edition, 2007

    """
    label = "poisson"

    def __init__(self, N: int | float, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        N:
            chain length
        repeat_MW:
            repeat unit molecular weight
        """
        super().__init__(repeat_MW)
        self._N = N
        self._D = 1 + 1 / N

    def _compute_x_i_pmd(self, N_i: int | np.ndarray) -> int | float | np.ndarray:  # noqa
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
