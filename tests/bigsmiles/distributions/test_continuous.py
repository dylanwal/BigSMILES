import pytest

import numpy as np

import bigsmiles.distributions as distributions
import bigsmiles.distributions.utils as utils

repeat_unit = 104.14

cases_log_normal = \
    [
        # Mn, D
        [1432.324, 1.011],
        [1432.324, 1.12],
        [1432.324, 1.773],
        [12345.2, 1.01],
        [12345.2, 1.13],
        [12345.2, 1.764],
        [259_093.23, 1.015],
        [259_093.23, 1.178],
        [259_093.23, 1.81]
    ]


@pytest.mark.parametrize("case", cases_log_normal)
def test_log_normal(case):
    Mn, D = case
    dis = distributions.LogNormal(Mn, D, repeat_MW=repeat_unit)

    cal_Mn, cal_D = utils.compute_Mn_D_from_xi_mw_i(dis.mw_i, dis.x_i())
    assert np.isclose(cal_Mn, Mn, rtol=1.e-3)
    assert np.isclose(cal_D, D, rtol=1.e-3)


cases_schulz_zimm = \
    [
        # Mn, D
        [1432.324, 1.0267],
        [1432.324, 1.12],
        [1432.324, 1.773],
        [12345.2, 1.0334],
        [12345.2, 1.13],
        [12345.2, 1.764],
        [259_093.23, 1.0255],
        [259_093.23, 1.178],
        [259_093.23, 1.81]
    ]


@pytest.mark.parametrize("case", cases_schulz_zimm)
def test_schulz_zimm(case):
    Mn, D = case
    dis = distributions.SchulzZimm(Mn, D, repeat_MW=repeat_unit)

    cal_Mn, cal_D = utils.compute_Mn_D_from_xi_mw_i(dis.mw_i, dis.x_i())
    assert np.isclose(cal_Mn, Mn, rtol=1.e-3)
    assert np.isclose(cal_D, D, rtol=1.e-3)


cases_gaussian = \
    [
        # Mn, D
        [1432.324, 1.011],
        [1432.324, 1.12],
        [1432.324, 1.243],
        [12345.2, 1.01],
        [12345.2, 1.13],
        [12345.2, 1.24],
        [259_093.23, 1.015],
        [259_093.23, 1.178],
        [259_093.23, 1.21]
    ]


@pytest.mark.parametrize("case", cases_gaussian)
def test_gaussian(case):
    Mn, D = case
    dis = distributions.Gaussian(Mn, D, repeat_MW=repeat_unit)

    cal_Mn, cal_D = utils.compute_Mn_D_from_xi_mw_i(dis.mw_i, dis.x_i())
    assert np.isclose(cal_Mn, Mn, rtol=1.e-2)
    assert np.isclose(cal_D, D, rtol=1.e-2)


cases_uniform = \
    [
        # [Mn low, Mn_high], [Mn, D]answer
        [1012.1, 5032.3, 3021, 1.147],
        [1012.1, 10032.3, 5522, 1.22],
        [1012.1, 45032.3, 23000, 1.305]
    ]


@pytest.mark.parametrize("case", cases_uniform)
def test_uniform(case):
    low_MW, high_MW, Mn, D = case
    dis = distributions.Uniform(low_MW, high_MW, repeat_MW=repeat_unit)

    cal_Mn, cal_D = utils.compute_Mn_D_from_xi_mw_i(dis.mw_i, dis.x_i())
    assert np.isclose(cal_Mn, Mn, rtol=1.e-2)
    assert np.isclose(cal_D, D, rtol=1.e-2)

