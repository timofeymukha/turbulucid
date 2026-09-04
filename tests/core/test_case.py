import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from turbulucid import Case


def test_native_case_structure(block_case):
    assert block_case.bounds == (0.0, 1.0, 0.0, 1.0)
    assert block_case.cellCentres.shape == (6, 2)
    assert block_case.fields == [
        "scalarField",
        "symtensorField",
        "tensorField",
        "vectorField",
    ]
    assert block_case.boundaries == [
        "inlet",
        "outlet",
        "bottomWall",
        "topWall",
    ]


def test_field_add_overwrite_and_delete(block_case):
    first = np.arange(block_case.cellCentres.shape[0], dtype=float)
    second = first + 10

    block_case["derived"] = first
    assert_array_equal(block_case["derived"], first)
    assert block_case.fields.count("derived") == 1

    block_case["derived"] = second
    assert_array_equal(block_case["derived"], second)
    assert block_case.fields.count("derived") == 1

    for boundary in block_case.boundaries:
        boundary_values = block_case.boundary_cell_data(boundary)[1]["derived"]
        assert np.all(np.isin(boundary_values, second))

    del block_case["derived"]
    assert "derived" not in block_case.fields


def test_boundary_cell_data_returns_independent_arrays(block_case):
    _, first = block_case.boundary_cell_data("inlet")
    expected = first["scalarField"].copy()
    first["scalarField"][:] = -1

    _, second = block_case.boundary_cell_data("inlet")
    assert_array_equal(second["scalarField"], expected)


@pytest.mark.parametrize("method", ["boundary_data", "boundary_cell_data"])
def test_boundary_sort_validation(block_case, method):
    with pytest.raises(ValueError, match="sort should"):
        getattr(block_case, method)("inlet", sort="z")


def test_native_round_trip(block_case, tmp_path):
    output = tmp_path / "round_trip.vtm"
    block_case.write(str(output))

    restored = Case(str(output))
    assert restored.boundaries == block_case.boundaries
    assert restored.fields == block_case.fields
    assert_allclose(restored.cellCentres, block_case.cellCentres)
    for field in block_case.fields:
        assert_allclose(restored[field], block_case[field])
