import pytest

import numpy as np

import bigsmiles.distributions.utils as utils


def test_cutoff():
    cutoff_ = 0.5
    output = utils.cutoff(np.array([1, 2, 3, 4, 5]), cutoff_)
    np.testing.assert_array_equal(output, np.array([0, 0, 3, 4, 5]))

