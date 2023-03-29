import pytest

import numpy as np

import bigsmiles.distributions as distributions
import bigsmiles.distributions.utils as utils

repeat_unit = 104.14

cases_flory_schulz = [0.1, 0.8, 0.99, 0.999]


@pytest.mark.parametrize("conv", cases_flory_schulz)
def test_flory_schulz(conv: float):
    dis = distributions.FlorySchulz(conversion=conv, repeat_MW=repeat_unit)

    a = 1 - conv
    assert np.isclose(dis.N, 2 / a - 1)
    assert np.isclose(dis.std_mw/dis.repeat_MW, np.sqrt((2 - 2 * a) / a**2))
    assert np.isclose(dis.skew_mw, (2 - a) / np.sqrt(2 - 2 * a))


cases_poisson = [10, 100, 1000, 8000]


@pytest.mark.parametrize("N", cases_poisson)
def test_cases_poisson(N: float):
    dis = distributions.Poisson(N=N, repeat_MW=repeat_unit)
    cal_Mn, cal_D = utils.compute_Mn_D_from_xi_mw_i(dis.mw_i, dis.x_i())
    assert np.isclose(cal_Mn/dis.repeat_MW, N, rtol=1.e-3)
    assert np.isclose(cal_D, 1 + 1 / N, rtol=1.e-3)

