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
    edge_lengths,
    isoline,
    normals,
    profile_along_line,
    sample_by_plane,
    sort_indices,
    tangents,
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


def test_profile_along_line_samples_cells_and_boundaries(block_case):
    """A horizontal cut through the lower cell row of the 3x2 block case."""
    block_case["ramp"] = block_case.cellCentres[:, 0].copy()

    distance, data = profile_along_line(block_case, (0.0, 0.25), (1.0, 0.25))

    # Three cell centres at x = 1/6, 1/2, 5/6, plus the inlet and outlet
    # faces at x = 0 and x = 1.
    assert_allclose(distance, [0.0, 1/6, 0.5, 5/6, 1.0], rtol=1e-5)
    assert_allclose(data["ramp"], [1/6, 1/6, 0.5, 5/6, 5/6], rtol=1e-5)


def test_profile_along_line_can_exclude_boundaries(block_case):
    block_case["ramp"] = block_case.cellCentres[:, 0].copy()

    distance, data = profile_along_line(
        block_case, (0.0, 0.25), (1.0, 0.25), excludeBoundaries=True)

    assert_allclose(distance, [1/6, 0.5, 5/6], rtol=1e-5)
    assert_allclose(data["ramp"], [1/6, 0.5, 5/6], rtol=1e-5)


def test_profile_along_line_is_sorted_regardless_of_direction(block_case):
    forward = profile_along_line(block_case, (0.0, 0.25), (1.0, 0.25))[0]
    backward = profile_along_line(block_case, (1.0, 0.25), (0.0, 0.25))[0]

    assert np.all(np.diff(forward) >= 0)
    assert np.all(np.diff(backward) >= 0)
    assert_allclose(forward, backward, rtol=1e-5)


def test_tangents_and_normals_point_the_right_way(block_case):
    """Outward normals on an axis-aligned box are the unit axis vectors."""
    expected = {
        "inlet": [-1.0, 0.0],
        "outlet": [1.0, 0.0],
        "bottomWall": [0.0, -1.0],
        "topWall": [0.0, 1.0],
    }
    for boundary, outward in expected.items():
        computed = normals(block_case, boundary)
        assert_allclose(computed, np.tile(outward, (computed.shape[0], 1)),
                        atol=1e-12)

        # A tangent is perpendicular to its normal and of unit length.
        tangent = tangents(block_case, boundary)
        assert_allclose(np.einsum("ij,ij->i", tangent, computed), 0.0,
                        atol=1e-12)
        assert_allclose(np.linalg.norm(tangent, axis=1), 1.0)


def test_edge_lengths_match_the_mesh_spacing(block_case):
    assert_allclose(edge_lengths(block_case, "inlet"), 0.5)
    assert_allclose(edge_lengths(block_case, "outlet"), 0.5)
    assert_allclose(edge_lengths(block_case, "bottomWall"), 1/3, rtol=1e-5)
    assert_allclose(edge_lengths(block_case, "topWall"), 1/3, rtol=1e-5)


def test_sort_indices_orders_along_the_requested_axis(block_case):
    ordering = sort_indices(block_case, "bottomWall", "x")
    points = block_case.boundary_data("bottomWall")[0][ordering]

    assert np.all(np.diff(points[:, 0]) > 0)

    with pytest.raises(ValueError, match="axis should be x or y"):
        sort_indices(block_case, "bottomWall", "z")


def test_sample_by_plane_covers_the_geometry(block_case):
    points, data = sample_by_plane(block_case, (5, 7))

    assert points.shape == (35, 2)
    assert "vtkValidPointMask" in data
    assert data["scalarField"].shape == (35,)
    # The plane is laid over the bounding box, so every point is inside.
    assert np.all(data["vtkValidPointMask"] > 0)


def test_isoline_follows_a_known_contour(block_case):
    block_case["ramp"] = block_case.cellCentres[:, 0].copy()

    line = isoline(block_case, "ramp", 0.5)

    assert line.shape[1] == 2
    assert line.shape[0] > 0
    assert_allclose(line[:, 0], 0.5, atol=1e-9)
