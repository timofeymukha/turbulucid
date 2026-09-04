import matplotlib.pyplot as plt
import numpy as np
import pytest

from turbulucid.core.data_extraction import sample_by_plane
from turbulucid.core.plotting import (
    _line_segments,
    _polygon_vertices,
    plot_boundaries,
    plot_contour,
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


def test_boundaries_default_to_black(block_case):
    collection = plot_boundaries(block_case)

    np.testing.assert_array_equal(collection.get_color(), [[0.0, 0.0, 0.0, 1.0]])


@pytest.mark.parametrize("keyword", ["color", "colors"])
@pytest.mark.parametrize("function", [plot_boundaries, plot_contour])
def test_line_plots_honour_either_colour_keyword(block_case, function, keyword):
    """LineCollection accepts both spellings; neither may be overridden."""
    arguments = {keyword: "red"}
    if function is plot_contour:
        collection = function(block_case, "scalarField", 2.5, **arguments)
    else:
        collection = function(block_case, **arguments)

    np.testing.assert_array_equal(collection.get_color(), [[1.0, 0.0, 0.0, 1.0]])


def test_plot_vectors_normalizes_without_mutating_the_input(block_case):
    field = block_case["vectorField"].astype(float)
    original = field.copy()

    result = plot_vectors(block_case, field, normalize=True, plotBoundaries=False)

    np.testing.assert_array_equal(field, original)
    lengths = np.hypot(result.U, result.V)
    np.testing.assert_allclose(lengths, 1.0)


def test_plot_vectors_normalize_leaves_zero_vectors_alone(block_case):
    field = np.zeros((block_case.vtkData.GetNumberOfCells(), 3))
    field[0] = [3.0, 4.0, 0.0]

    result = plot_vectors(block_case, field, normalize=True, plotBoundaries=False)

    lengths = np.hypot(result.U, result.V)
    assert lengths[0] == pytest.approx(1.0)
    np.testing.assert_allclose(lengths[1:], 0.0)


@pytest.mark.parametrize(
    ("field", "value", "error", "message"),
    [
        ("missing", 1.0, ValueError, "not present"),
        ("vectorField", 1.0, ValueError, "scalar"),
        ("scalarField", np.inf, ValueError, "finite"),
        (1, 1.0, TypeError, "string"),
        ("scalarField", "high", TypeError, "real scalar"),
    ],
)
def test_plot_contour_validates_its_request(block_case, field, value, error, message):
    with pytest.raises(error, match=message):
        plot_contour(block_case, field, value)


def test_streamline_colour_array_matches_the_named_field(block_case):
    """An array colorField must be interpreted in sample_by_plane's order."""
    block_case["ramp"] = block_case.cellCentres[:, 0].copy()
    resolution = (8, 9)
    _, sampled = sample_by_plane(block_case, resolution)

    from_name = plot_streamlines(
        block_case,
        "vectorField",
        colorField="ramp",
        planeResolution=resolution,
        plotBoundaries=False,
    )
    expected = np.asarray(from_name.lines.get_array())
    plt.close("all")

    from_array = plot_streamlines(
        block_case,
        "vectorField",
        colorField=sampled["ramp"],
        planeResolution=resolution,
        plotBoundaries=False,
    )

    np.testing.assert_allclose(np.asarray(from_array.lines.get_array()), expected)
    # A ramp across a non-square grid: the two reshape orders really differ.
    assert np.nanmin(expected) < np.nanmax(expected)


def test_streamline_colour_array_length_is_checked(block_case):
    with pytest.raises(ValueError, match="one value per sampling point"):
        plot_streamlines(
            block_case,
            "vectorField",
            colorField=np.ones(5),
            planeResolution=(8, 9),
            plotBoundaries=False,
        )


def test_plot_field_clips_to_the_requested_limits(block_case):
    """Narrowing xlim must drop the cells that fall outside it."""
    full = plot_field(block_case, "scalarField", colorbar=False)
    fullCount = len(full.get_paths())
    plt.close("all")

    clipped = plot_field(
        block_case, "scalarField", xlim=[0.0, 0.4], colorbar=False)

    assert len(clipped.get_paths()) < fullCount
    np.testing.assert_allclose(clipped.axes.get_xlim(), [0.0, 0.4])


def test_plot_field_without_limits_covers_the_whole_case(block_case):
    collection = plot_field(block_case, "scalarField", colorbar=False)

    assert len(collection.get_paths()) == block_case.vtkData.GetNumberOfCells()
    np.testing.assert_allclose(collection.axes.get_xlim(), block_case.xlim)
    np.testing.assert_allclose(collection.axes.get_ylim(), block_case.ylim)


@pytest.mark.parametrize("limits", [[0.0], [0.0, 1.0, 2.0], [1.0, 0.0],
                                    [0.0, np.nan]])
@pytest.mark.parametrize("axis", ["xlim", "ylim"])
def test_plot_field_validates_limits(block_case, axis, limits):
    with pytest.raises(ValueError, match=axis):
        plot_field(block_case, "scalarField", colorbar=False, **{axis: limits})


@pytest.mark.parametrize(
    "function",
    [plot_boundaries, plot_field, plot_vectors, plot_streamlines, plot_contour],
)
def test_scaling_factors_must_be_positive(block_case, function):
    arguments = {"plot_field": ("scalarField",), "plot_vectors": ("vectorField",),
                 "plot_streamlines": ("vectorField",),
                 "plot_contour": ("scalarField", 2.5)}
    positional = arguments.get(function.__name__, ())

    with pytest.raises(ValueError, match="Scaling factors"):
        function(block_case, *positional, scaleX=0)


def test_plot_field_scales_the_geometry(block_case):
    collection = plot_field(block_case, "scalarField", scaleX=2, colorbar=False)

    np.testing.assert_allclose(collection.axes.get_xlim(), block_case.xlim/2)


# ---------------------------------------------------------------------------
# Geometry extraction. plot_field reads the cell array in one vectorised go
# rather than walking cells; these pin it to the per-cell reference.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(("scaleX", "scaleY"), [(1, 1), (2, 0.5)])
def test_polygon_vertices_match_the_per_cell_walk(
        block_case, polygon_reference, scaleX, scaleY):
    expected = polygon_reference(block_case.vtkData, scaleX, scaleY)

    vertices = _polygon_vertices(block_case.vtkData, scaleX, scaleY)

    assert len(vertices) == len(expected)
    for computed, reference in zip(vertices, expected, strict=True):
        np.testing.assert_allclose(computed, reference)


def test_uniform_meshes_use_the_array_fast_path(block_case):
    """An all-quad mesh must come back as one (nCells, 4, 2) array.

    PolyCollection is markedly faster for that than for a list, so this
    is a performance contract, not just a detail.

    """
    vertices = _polygon_vertices(block_case.vtkData, 1, 1)

    assert isinstance(vertices, np.ndarray)
    assert vertices.shape == (block_case.vtkData.GetNumberOfCells(), 4, 2)


def test_mixed_meshes_fall_back_to_a_list(mixed_cell_case, polygon_reference):
    expected = polygon_reference(mixed_cell_case.vtkData)

    vertices = _polygon_vertices(mixed_cell_case.vtkData, 1, 1)

    assert not isinstance(vertices, np.ndarray)
    assert [len(polygon) for polygon in vertices] == [4, 3, 5]
    for computed, reference in zip(vertices, expected, strict=True):
        np.testing.assert_allclose(computed, reference)


def test_plot_field_handles_mixed_cell_sizes(mixed_cell_case):
    collection = plot_field(mixed_cell_case, "p", colorbar=False,
                            plotBoundaries=False)

    # Matplotlib closes each polygon, so every path gains one vertex.
    counts = sorted(len(path.vertices) for path in collection.get_paths())
    assert counts == [4, 5, 6]
    np.testing.assert_allclose(collection.get_array(), [0.0, 1.0, 2.0])


def test_line_segments_match_the_per_cell_walk(block_case):
    block = block_case.extract_block_by_name("bottomWall")

    segments = _line_segments(block, 2, 0.5)

    assert segments.shape == (block.GetNumberOfCells(), 2, 2)
    for cellId in range(block.GetNumberOfCells()):
        points = block_case.extract_block_by_name(
            "bottomWall").GetCell(cellId).GetPoints()
        expected = np.array([points.GetPoint(0)[:2],
                             points.GetPoint(1)[:2]])/[2, 0.5]
        np.testing.assert_allclose(segments[cellId], expected)


# ---------------------------------------------------------------------------
# The collections are added with autolim=False and the data limits supplied
# directly, so the limits have to be checked explicitly.
# ---------------------------------------------------------------------------

def test_plot_field_sets_the_data_limits(block_case):
    collection = plot_field(block_case, "scalarField", colorbar=False,
                            plotBoundaries=False)

    np.testing.assert_allclose(
        collection.axes.dataLim.get_points(), [[0.0, 0.0], [1.0, 1.0]])


def test_plot_field_data_limits_follow_clipping_and_scaling(block_case):
    collection = plot_field(block_case, "scalarField", xlim=[0.2, 0.7],
                            scaleX=2, colorbar=False, plotBoundaries=False)

    limits = collection.axes.dataLim.get_points()
    np.testing.assert_allclose(limits[:, 0], [0.1, 0.35])
    np.testing.assert_allclose(limits[:, 1], [0.0, 1.0])


def test_plot_boundaries_sets_the_data_limits(block_case):
    collection = plot_boundaries(block_case)

    np.testing.assert_allclose(
        collection.axes.dataLim.get_points(), [[0.0, 0.0], [1.0, 1.0]])


def test_plot_contour_sets_the_data_limits(block_case):
    # scalarField is a uniform 2.5, which has no interior contour, so
    # contour a field that actually varies across the geometry.
    block_case["ramp"] = block_case.cellCentres[:, 0].copy()

    collection = plot_contour(block_case, "ramp", 0.5)

    limits = collection.axes.dataLim.get_points()
    np.testing.assert_allclose(limits[:, 0], [0.5, 0.5])
    np.testing.assert_allclose(limits[:, 1], [0.0, 1.0])


def test_an_empty_contour_leaves_the_data_limits_alone(block_case):
    """A contour that matches nothing must not poison the limits."""
    collection = plot_contour(block_case, "scalarField", 1e30)

    assert collection.get_segments() == []
    # Matplotlib's "no data yet" bbox is [[inf, inf], [-inf, -inf]].
    assert np.all(np.isinf(collection.axes.dataLim.get_points()))


def test_autoscale_after_plotting_still_frames_the_geometry(block_case):
    """A user re-enabling autoscale must not get an empty view."""
    collection = plot_field(block_case, "scalarField", colorbar=False,
                            plotBoundaries=False)
    ax = collection.axes

    ax.autoscale()

    assert ax.get_xlim()[0] <= 0.0 and ax.get_xlim()[1] >= 1.0
    assert ax.get_ylim()[0] <= 0.0 and ax.get_ylim()[1] >= 1.0
