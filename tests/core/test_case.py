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


def test_case_properties_follow_direct_vtk_changes(block_case):
    values = np.arange(block_case.vtkData.GetNumberOfCells(), dtype=float)
    block_case.vtkData.CellData.append(values, "directField")

    assert "directField" in block_case.fields
    assert_array_equal(block_case["directField"], values)

    block_case.vtkData.CellData["scalarField"][:] = 42
    for boundary in block_case.boundaries:
        boundary_values = block_case.boundary_cell_data(boundary)[1]
        assert_array_equal(
            boundary_values["scalarField"],
            np.full(boundary_values["scalarField"].shape, 42),
        )

    original_centres = block_case.cellCentres
    block_case.vtkData.Points[:, 0] += 2
    assert_allclose(block_case.cellCentres[:, 0], original_centres[:, 0] + 2)
    assert block_case.bounds[:2] == pytest.approx((2.0, 3.0))


def test_case_round_trips_tensor_fields(block_case):
    n_cells = block_case.vtkData.GetNumberOfCells()
    values = np.arange(n_cells * 9, dtype=float).reshape(n_cells, 3, 3)

    block_case["derivedTensor"] = values

    assert_array_equal(block_case["derivedTensor"], values)


@pytest.mark.parametrize(
    ("name", "values", "error"),
    [
        (1, np.ones(6), TypeError),
        ("", np.ones(6), ValueError),
        ("badLength", np.ones(2), ValueError),
        ("notNumeric", np.array(["x"] * 6), TypeError),
        ("complex", np.full(6, 1 + 2j), TypeError),
        ("notFinite", np.full(6, np.nan), ValueError),
        ("badTensor", np.ones((6, 2, 2)), ValueError),
    ],
)
def test_field_assignment_validation(block_case, name, values, error):
    with pytest.raises(error):
        block_case[name] = values


def test_transform_validation(block_case):
    with pytest.raises(ValueError, match="nonzero"):
        block_case.scale(0, 1)
    with pytest.raises(ValueError, match="finite"):
        block_case.translate(np.inf, 0)
    with pytest.raises(TypeError, match="scalar"):
        block_case.rotate([90])


def test_boundary_name_validation(block_case):
    with pytest.raises(ValueError, match="not present"):
        block_case.boundary_data("missing")
    with pytest.raises(TypeError, match="string"):
        block_case.boundary_cell_data(1)


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
