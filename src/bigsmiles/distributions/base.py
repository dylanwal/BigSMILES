from __future__ import annotations
import abc

import numpy as np

import bigsmiles.distributions.utils as utils

# TODO: add numerical bound of validity for all distirbutions; when distributions cut off at Zero cause issues


class Distribution(abc.ABC):
    """
    Generic class for molecular weight distributions
    """
    _DEFAULT_RNG = np.random.default_rng(0)
    _N_i_max = 100_000
    _mw_i_max = 10_000_000

    def __init__(self, repeat_MW: int | float | None = None):  # noqa
        """
        Parameters
        ----------
        repeat_MW:
            repeat unit MW
        """
        self._repeat_MW = repeat_MW

        # the goal is not compute unless needed;
        # and once computed store value so doesn't need to be computed again
        # molecular weight
        self._Mn: int | float | None = None
        self._D: int | float | None = None
        self._mean_mw: int | float | None = None
        self._std_mw: int | float | None = None
        self._skew_mw: int | float | None = None
        self._kurtosis_mw: int | float | None = None
        self._peak_mw: int | float | None = None
        # chain length
        self._N: int | float | None = None
        self._peak_N: int | float | None = None

        # numpy arrays
        self._mw_i: np.ndarray | None = None
        self._N_i: np.ndarray | None = None
        self._x_i: np.ndarray | None = None
        self._x_i_pmd: np.ndarray | None = None
        self._w_i: np.ndarray | None = None
        self._w_i_pmd: np.ndarray | None = None

    def __str__(self):
        return self._create_label()

    def __bool__(self):
        return self.__str__()

    @abc.abstractmethod
    def _create_label(self):
        """ create string for bigSMILES """

    def details(self) -> str:
        """ long string representation """
        text = self.__str__()
        text += f"(Mn: {self.Mn:.0f}, D: {self.D:.3f}, N: {self.N:.0f}, repeat_MW: {self.repeat_MW:.2f})"

        return text

    @property
    def repeat_MW(self) -> float | int | None:  # noqa
        """ repeat unit MW """
        return self._repeat_MW

    @repeat_MW.setter
    def repeat_MW(self, repeat_MW: float | int):  # noqa
        if (isinstance(repeat_MW, int) or isinstance(repeat_MW, float)) and repeat_MW > 0:
            self._repeat_MW = repeat_MW
        else:
            raise ValueError(f"'repeat_MW' must be a positive integer or float. Given: {repeat_MW}")

    @property
    def mw_i(self) -> np.ndarray:
        """ molecular weight array (x-axis) """
        if self._mw_i is None:
            self._compute_mw_i()

        return self._mw_i

    @property
    def N_i(self) -> np.ndarray | None:  # noqa
        """ chain length array (x-axis)"""
        if self._N_i is None:
            self._compute_N_i()

        return self._N_i

    @property
    def Mn(self) -> float | int:  # noqa
        """ number average MW """
        if self._Mn is None:
            self._compute_Mn()

        return self._Mn

    @property
    def Mw(self) -> float | int:  # noqa
        """ weight average MW """
        return self._D * self._Mn

    @property
    def D(self) -> float | int:  # noqa
        """ molecular weight dispersity """
        if self._D is None:
            self._compute_D()

        return self._D

    @property
    def N(self) -> float | int | None:  # noqa
        """ average chain length """
        if self._N is None:
            self._compute_N()

        return self._N

    @property
    def peak_mw(self) -> float | int | None:
        """ molecular weight at distribution max value """
        if self._peak_mw is None:
            self._compute_peak_mw()

        return self._peak_mw

    @property
    def peak_N(self) -> float | int | None:  # noqa
        """ chain length at distribution max value """
        if self._N is None:
            self._compute_peak_N()

        return self._N

    @property
    def mean_mw(self) -> int | float:
        """ mean molecular weight """
        return self.Mn

    @property
    def std_mw(self) -> int | float:
        """ standard deviation molecular weight """
        if self._std_mw is None:
            self._compute_std_mw()

        return self._std_mw

    @property
    def skew_mw(self) -> int | float:
        """
        skew (3 moment of a distribution)

        Notes
        -----

        * Symmetric: Values between -0.5 to 0.5
        * Moderated Skewed: Values between -1 and -0.5 or between 0.5 and 1
        * Highly Skewed: Values less than -1 or greater than 1

        """
        if self._skew_mw is None:
            self._compute_skew_mw()

        return self._skew_mw

    @property
    def kurtosis_mw(self) -> int | float:
        """
        kurtosis (4 moment of a distribution)

        Notes
        -----

        * Mesokurtic: normal distribution
        * Leptokurtic: kurtosis > 3; the distribution has fatter tails and a sharper peak
        * Platykurtic: kurtosis < -3; the distribution has a lower and wider peak and thinner tails

        """
        if self._kurtosis_mw is None:
            self._compute_kurtosis_mw()

        return self._kurtosis_mw

    def x_i(self, mw_i: int | float | np.ndarray = None) -> int | float | np.ndarray:
        """ mole fraction of molecular weight as a function of $mw_i$ """
        if mw_i is None:
            if self._x_i is None:
                self._x_i = self._compute_x_i(self.mw_i)
            return self._x_i

        return self._compute_x_i(mw_i)

    def x_i_cdf(self) -> tuple[np.ndarray, np.ndarray]:
        """
        cumulative mole fraction of molecular weight as a function of $mw_i$

        Returns
        -------
        x:
            x position of cdf
        cdf:
            cdf
        """
        return utils.compute_cdf_x(self.x_i(self.mw_i), self.mw_i)

    def w_i(self, mw_i: int | float | np.ndarray = None) -> np.ndarray | None:
        """
        weight fraction of molecular weight as a function of $mw_i$

        """
        if self._w_i is None:
            self._w_i = utils.xi_to_wi(self.mw_i, self.x_i(self.mw_i), self.Mn)
        if mw_i is None:
            return self._w_i

        return np.interp(mw_i, self.mw_i, self._w_i)

    def w_i_cdf(self, n: int = None) -> np.ndarray | None:
        """
        cumulative weight fraction of molecular weight as a function of $mw_i$

        Parameters
        ----------
        n:
            number of points

        Returns
        -------
        x:
            x position of cdf
        cdf:
            cdf
        """
        return utils.compute_cdf_x(self.w_i(self.mw_i), self.mw_i)

    def x_i_pmd(self, N_i: int | np.ndarray = None) -> int | np.ndarray:  # noqa
        """ mole fraction of molecular weight as a function of $N_i$  (Probability mass distribution) """
        if N_i is None:
            N_i = self.N_i
        elif isinstance(N_i, np.ndarray) and (N_i.dtype.kind == 'i' or N_i.dtype.kind == 'u'):
            if np.max(N_i) > np.max(self.N_i):
                raise ValueError(f'Probability mass distribution limit exceeded. Max value: {self._N_i_max}')
        elif isinstance(N_i, int):
            pass
        else:
            raise TypeError("'x_i_pmd' accepts only integer or numpy array dtype 'i'.")

        return self._compute_x_i_pmd(N_i)

    def draw_mw(self, n: int = 1, rng=None) -> int | float | np.ndarray:
        """
        draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        n:
             Numpy random number generator for the generation of numbers.
        rng:
            Provide your own random number generator;
            [np.random.default_rng()](https://numpy.org/doc/stable/reference/random/generator.html)

        """

    def draw_N(self, n: int = 1, rng: np.random.Generator = None) -> int | np.ndarray:  # noqa
        """
        draw a sample from the molecular weight distribution.

        Parameters:
        ----------
        n:
            number of points
        rng:
             Numpy random number generator for the generation of numbers.

        """

    @abc.abstractmethod
    def _compute_mw_i(self):
        ...

    @abc.abstractmethod
    def _compute_N_i(self):  # noqa
        ...

    @abc.abstractmethod
    def _compute_Mn(self):  # noqa
        ...

    @abc.abstractmethod
    def _compute_D(self):  # noqa
        ...

    @abc.abstractmethod
    def _compute_N(self):  # noqa
        ...

    @abc.abstractmethod
    def _compute_peak_N(self): # noqa
        ...

    @abc.abstractmethod
    def _compute_peak_mw(self):
        ...

    def _compute_std_mw(self):
        self._std_mw = np.sqrt(np.trapz(x=self.mw_i, y=self.x_i(self.mw_i) * (self.mw_i - self.Mn) ** 2))

    def _compute_skew_mw(self):
        self._skew_mw = np.trapz(x=self.mw_i, y=self.x_i(self.mw_i) * (self.mw_i - self.Mn) ** 3) / self.std_mw ** 3

    def _compute_kurtosis_mw(self):
        self._kurtosis_mw = (np.trapz(
            x=self.mw_i,
            y=self.x_i(self.mw_i) * (self.mw_i - self.Mn) ** 4
        ) / self.std_mw ** 4) - 3

    @abc.abstractmethod
    def _compute_x_i(self, mw_i: int | float | np.ndarray) -> int | float | np.ndarray:
        ...

    @abc.abstractmethod
    def _compute_x_i_pmd(self, N_i: int | float | np.ndarray) -> int | float | np.ndarray: # noqa
        ...
