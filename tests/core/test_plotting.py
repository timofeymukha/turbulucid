import numpy as np
import pytest

from turbulucid.core.plotting import (
    plot_field,
    plot_streamlines,
    plot_vectors,
)


def test_plot_field_accepts_field_name_without_mutating_case(block_case):
    block_case["temp"] = np.arange(block_case.cellCentres.shape[0], dtype=float)
    expected = block_case["temp"]
    fields = block_case.fields.copy()

    result = plot_field(block_case, "scalarField", colorbar=False)

    assert result is not None
    assert block_case.fields == fields
    np.testing.assert_array_equal(block_case["temp"], expected)


def test_plot_field_accepts_array_without_mutating_case(block_case):
    fields = block_case.fields.copy()

    result = plot_field(block_case, block_case["scalarField"], colorbar=False)

    assert result is not None
    assert block_case.fields == fields
    assert "__turbulucid_plot_data__" not in block_case.vtkData.CellData.keys()


def test_plot_field_rejects_wrong_array_length(block_case):
    with pytest.raises(ValueError, match="dimensionality"):
        plot_field(block_case, np.ones(2), colorbar=False)


def test_sampled_vectors_support_color_field(block_case):
    result = plot_vectors(
        block_case,
        "vectorField",
        colorField="scalarField",
        sampleByPlane=True,
        planeResolution=(8, 9),
        plotBoundaries=False,
    )

    assert result is not None


def test_streamlines_use_requested_resolution(block_case):
    result = plot_streamlines(
        block_case,
        "vectorField",
        planeResolution=(8, 9),
        plotBoundaries=False,
    )

    assert result is not None


def test_colored_streamlines_return_plot_object(block_case):
    result = plot_streamlines(
        block_case,
        "vectorField",
        colorField="scalarField",
        planeResolution=(8, 9),
        plotBoundaries=False,
    )

    assert result is not None
