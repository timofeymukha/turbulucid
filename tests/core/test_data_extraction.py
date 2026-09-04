# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from os import path

import numpy as np
import pytest
from numpy.testing import assert_allclose

import turbulucid
from turbulucid.core.data_extraction import (
    dist,
    isoline,
    profile_along_line,
    sample_by_plane,
)


def test_dist_orthogonal():
    casePath = path.join(turbulucid.__path__[0], "datasets",
                         "test_case_block", "averaged.vtm")
    case = turbulucid.Case(casePath)

    distance = dist(case, "bottomWall")
    assert_allclose(distance, 0.25)

    distance = dist(case, "topWall")
    assert_allclose(distance, 0.25)

    distance = dist(case, "inlet")
    assert_allclose(distance, 1/6, rtol=1e-5)

    distance = dist(case, "outlet")
    assert_allclose(distance, 1/6, rtol=1e-5)


@pytest.mark.parametrize(
    ("p1", "p2", "error"),
    [
        ((0, 0), (0, 0), ValueError),
        ((0,), (1, 1), ValueError),
        ((0, np.nan), (1, 1), ValueError),
        (("bad", 0), (1, 1), TypeError),
    ],
)
def test_profile_along_line_validates_points(block_case, p1, p2, error):
    with pytest.raises(error):
        profile_along_line(block_case, p1, p2)


def test_profile_along_line_validates_flags(block_case):
    with pytest.raises(TypeError, match="correctDistance"):
        profile_along_line(block_case, (0, 0), (1, 1), correctDistance=1)
    with pytest.raises(TypeError, match="excludeBoundaries"):
        profile_along_line(block_case, (0, 0), (1, 1), excludeBoundaries=1)


@pytest.mark.parametrize(
    ("resolution", "error"),
    [
        ((10,), ValueError),
        ((10, 10, 10), ValueError),
        ((1, 10), ValueError),
        ((10.0, 10), TypeError),
        (10, TypeError),
    ],
)
def test_sample_by_plane_validates_resolution(block_case, resolution, error):
    with pytest.raises(error):
        sample_by_plane(block_case, resolution)


def test_dist_validates_corrected_flag(block_case):
    with pytest.raises(TypeError, match="corrected"):
        dist(block_case, "inlet", corrected=1)


def test_isoline_validation_and_empty_result(block_case):
    with pytest.raises(ValueError, match="not present"):
        isoline(block_case, "missing", 1)
    with pytest.raises(ValueError, match="scalar"):
        isoline(block_case, "vectorField", 1)
    with pytest.raises(ValueError, match="finite"):
        isoline(block_case, "scalarField", np.inf)

    result = isoline(block_case, "scalarField", 1e30)
    assert result.shape == (0, 2)
