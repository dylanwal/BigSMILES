from __future__ import annotations
from functools import wraps

import numpy as np
from scipy.integrate import cumulative_trapezoid


def check_for_repeat_MW(func):  # noqa
    @wraps(func)
    def _check_for_repeat_MW(*args, **kwargs):  # noqa
        self = args[0]
        if hasattr(self, "repeat_MW"):
            if self.repeat_MW is None:
                # logging.warning("'repeat_MW' needs to be defined.")
                return None

        return func(*args, **kwargs)

    return _check_for_repeat_MW


def compute_Mn_D_from_xi_mw_i(mw_i: np.ndarray, x_i: np.ndarray) -> tuple[float, float]:  # noqa
    """
    calculate Mn and D from x_i vs MW_i

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


def compute_Mn_D_from_wi_mw_i(mw_i: np.ndarray, w_i: np.ndarray) -> tuple[float, float]:  # noqa
    """
    calculate Mn and D from w_i vs MW_i

    Parameters
    ----------
    mw_i:
        array containing molecular weight of polymer 'i' (i is degree of polymerization, chain length)
    w_i:
        array containing weight fraction of polymer 'i' (i is degree of polymerization, chain length)

    Returns
    -------
    mw_n: float
        number average molecular weight
    mw_d: float
        molecular weight dispersity
    """
    # avoid dividing by zero
    mask = mw_i > 1
    mw_i = mw_i[mask]
    w_i = w_i[mask]

    mw_n = np.sum(w_i) / np.sum(w_i/mw_i)
    mw_w = np.sum(w_i * mw_i) / np.sum(w_i)
    mw_d = mw_w / mw_n
    return mw_n, mw_d


def compute_Mn_D_from_xi_ni(N_i: np.ndarray, x_i: np.ndarray, repeat_mw: int | float) -> tuple[float, float]:  # noqa
    """
    calculate Mn and D from xi  vs N_i

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


def wi_to_xi(mw_i: np.ndarray, w_i: np.ndarray, Mn: int | float) -> np.ndarray:  # noqa
    return w_i*Mn/mw_i


def xi_to_wi(mw_i: np.ndarray, x_i: np.ndarray, Mn: int | float) -> np.ndarray:  # noqa
    return x_i*mw_i/Mn


def compute_cdf(pdf, npts: int = 10_000, x_range: tuple[float, float] = (0, 1)) -> tuple[np.ndarray, np.ndarray]:
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


def compute_cdf_x(pdf: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    calculate cumulative distribution function (cdf)

    Parameters
    ----------
    pdf:
        probability distribution function
    x:
        points of evaluation !!! will be reduced !!!

    Returns
    -------
    x:
        x position of cdf
    cdf:
        cdf

    """
    y_norm = pdf / np.trapz(pdf, x)
    cdf = cumulative_trapezoid(y_norm, x)
    cdf = np.insert(cdf, 0, 0)

    # deal with numerical round off errors
    cdf, index_ = np.unique(cdf, return_index=True)
    x = x[index_]

    return x, cdf


def normalize_distribution(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """ normalize a distribution to area = 1"""
    area = np.trapz(x=x, y=y)
    return y / area


def cutoff(y: np.ndarray, cutoff_: int | float) -> np.ndarray:
    mask = y < (cutoff_ * np.max(y))
    y[mask] = 0
    return y
