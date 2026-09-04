# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import matplotlib.pyplot as plt
import numpy as np
import vtkmodules
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits import axes_grid1
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData

from .data_extraction import _validate_contour_request, sample_by_plane

__all__ = ["plot_boundaries", "plot_vectors", "plot_streamlines", "plot_field",
           "add_colorbar", "plot_contour"]


def _colour_requested(kwargs):
    """Whether the user already asked for a specific line colour.

    LineCollection accepts both spellings, so checking only one of them
    silently overrides the other.

    """
    return any(key in kwargs for key in ("color", "colors"))


def _axis_limits(limits, name):
    """Validate a user-supplied pair of axis limits."""
    limits = np.asarray(limits, dtype=float)
    if limits.shape != (2,):
        raise ValueError(f"{name} must contain exactly two values.")
    if not np.all(np.isfinite(limits)):
        raise ValueError(f"{name} values must be finite.")
    if limits[0] >= limits[1]:
        raise ValueError(f"{name} must be increasing.")
    return limits


def _update_datalim(ax, points):
    """Grow the axes data limits to cover a set of xy points."""
    if len(points) == 0:
        return
    ax.update_datalim([points.min(axis=0), points.max(axis=0)])


def _as_polydata(data):
    """Unwrap a dataset_adapter wrapper, if that is what we were given."""
    return getattr(data, "VTKObject", data)


def _scaled_points(polyData, scaleX, scaleY):
    """The xy point coordinates of a polydata, scaled, as one array."""
    vtkPoints = polyData.GetPoints()
    if vtkPoints is None or vtkPoints.GetNumberOfPoints() == 0:
        return np.empty((0, 2))
    return vtk_to_numpy(vtkPoints.GetData())[:, :2]/[scaleX, scaleY]


def _cell_array_topology(cellArray):
    """Return (offsets, connectivity) of a vtkCellArray as NumPy arrays.

    Returns None if this build of VTK predates the offset/connectivity
    representation, so that callers can fall back to walking the cells.

    """
    try:
        offsets = cellArray.GetOffsetsArray()
        connectivity = cellArray.GetConnectivityArray()
    except AttributeError:
        return None
    if offsets is None or connectivity is None:
        return None
    return vtk_to_numpy(offsets), vtk_to_numpy(connectivity)


def _polygon_vertices(data, scaleX, scaleY):
    """Return the polygons of a polydata as scaled vertex arrays.

    When every polygon has the same number of vertices -- the usual
    all-quad or all-triangle cut plane -- a single (nCells, nVertices, 2)
    array is returned, because PolyCollection has a much faster path for
    that than for a list of separate arrays. Otherwise a list of
    (nVertices, 2) arrays is returned.

    """
    polyData = _as_polydata(data)
    nCells = polyData.GetNumberOfCells()
    if nCells == 0:
        return []

    topology = None
    if polyData.GetNumberOfPolys() == nCells:
        # Only take the vectorised route when the polys are the whole
        # story: GetCell() walks verts, lines, polys and strips in turn,
        # and the cell data is ordered to match.
        topology = _cell_array_topology(polyData.GetPolys())

    if topology is None:
        return _polygon_vertices_by_cell(polyData, scaleX, scaleY)

    offsets, connectivity = topology
    points = _scaled_points(polyData, scaleX, scaleY)
    vertices = points[connectivity]

    sizes = np.diff(offsets)
    if np.all(sizes == sizes[0]):
        return vertices.reshape(nCells, int(sizes[0]), 2)
    return np.split(vertices, offsets[1:-1])


def _polygon_vertices_by_cell(polyData, scaleX, scaleY):
    """Per-cell fallback for datasets of mixed cell types."""
    polygons = []
    for c in range(polyData.GetNumberOfCells()):
        cellPoints = polyData.GetCell(c).GetPoints()
        vertices = np.array(
            [cellPoints.GetPoint(i)[:2]
             for i in range(cellPoints.GetNumberOfPoints())])
        polygons.append(vertices/[scaleX, scaleY])
    return polygons


def _line_segments(data, scaleX, scaleY):
    """Collect the two-point cells of a polydata as scaled xy segments.

    Returns an (nSegments, 2, 2) array.

    """
    polyData = _as_polydata(data)
    nCells = polyData.GetNumberOfCells()
    if nCells == 0:
        return np.empty((0, 2, 2))

    topology = None
    if polyData.GetNumberOfLines() == nCells:
        topology = _cell_array_topology(polyData.GetLines())

    if topology is None:
        return _line_segments_by_cell(polyData, scaleX, scaleY)

    offsets, connectivity = topology
    points = _scaled_points(polyData, scaleX, scaleY)

    # Anything that is not a plain two-point segment is skipped.
    starts = offsets[:-1][np.diff(offsets) == 2]
    if starts.size == 0:
        return np.empty((0, 2, 2))
    return np.stack(
        (points[connectivity[starts]], points[connectivity[starts + 1]]),
        axis=1,
    )


def _line_segments_by_cell(polyData, scaleX, scaleY):
    """Per-cell fallback for datasets of mixed cell types."""
    segments = []
    for c in range(polyData.GetNumberOfCells()):
        cell = polyData.GetCell(c)
        if cell.GetNumberOfPoints() != 2:
            continue
        cellPoints = cell.GetPoints()
        segments.append((cellPoints.GetPoint(0)[:2], cellPoints.GetPoint(1)[:2]))
    if not segments:
        return np.empty((0, 2, 2))
    return np.array(segments)/[scaleX, scaleY]


def _add_line_collection(segments, case, scaleX, scaleY, kwargs):
    """Add line segments to the current axes, framed on the geometry."""
    collection = LineCollection(segments, **kwargs)
    if not _colour_requested(kwargs):
        collection.set_color("Black")

    ax = plt.gca()
    # autolim would walk every path to find the extent, which we already
    # know. The limits are set explicitly below, so feed the data limits
    # in directly instead.
    ax.add_collection(collection, autolim=False)
    _update_datalim(ax, np.asarray(segments).reshape(-1, 2))
    ax.set_xlim(case.xlim/scaleX)
    ax.set_ylim(case.ylim/scaleY)
    ax.set_aspect('equal')

    return collection


def _temporary_field_name(case):
    """Return a field name that cannot collide with user data."""
    name = "__turbulucid_plot_data__"
    while name in case.fields:
        name += "_"
    return name


def add_colorbar(data, aspect=20, padFraction=0.5, **kwargs):
    """Add a vertical colorbar to an image plot.

    Parameters
    ----------
    data
        The data with a .axes attribute.
    aspect : float, optional
        The ratio between the height and the width of the colorbar.
    padFraction : float, optional
        The horizontal distance between the figure and the colorbar
        as a fraction of the width of the colorbar.


    Returns
    -------
    colorbar
        The colorbar object

    """
    divider = axes_grid1.make_axes_locatable(data.axes)
    width = axes_grid1.axes_size.AxesY(data.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(padFraction, width)
    currentAx = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(currentAx)

    return data.axes.figure.colorbar(data, cax=cax, **kwargs)


def plot_boundaries(case, scaleX=1, scaleY=1, **kwargs):
    """Plot the boundaries the domain.

    Parameters
    ----------
    case : Case
        The case to draw the boundaries for.
    scaleX : float, optional
        A scaling factor for the abscissa.
    scaleY : float, optional
        A scaling factor for the ordinate.
    **kwargs
        Additional options to pass to LineCollection constructor.

    Raises
    ------
    ValueError
        If one or both scaling factors are non-positive.

    Returns
    -------
    LineCollection
        Collection of line segments defining the boundary.

    """
    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    segments = [_line_segments(case.extract_block_by_name(boundary),
                               scaleX, scaleY)
                for boundary in case.boundaries]
    segments = (np.concatenate(segments) if segments
                else np.empty((0, 2, 2)))

    return _add_line_collection(segments, case, scaleX, scaleY, kwargs)


def plot_vectors(case, field, colorField=None,
                 normalize=False, scaleX=1, scaleY=1,
                 sampleByPlane=False, planeResolution=None,
                 plotBoundaries=True,
                 **kwargs):
    """Plot a vector field.

    This function wraps pyplot.quiver. See that function's documentation
    for additional keyword argumnets.

    Parameters
    ----------
    case : Case
        The case that the vector field belongs to.
    field : str or ndarray
        Either a string with the name of the field as found in the case
        or an ndarray with the data.
    colorField : string or ndarray
        Data used to colour the vectors, either name of the field or
        an array. Does not work with sampleByPlane unless you carefully
        craft the array to fit the resampled data size.
    normalize : bool, optional
        Whether to normalize the the length of the vectors.
        Default is False.
    scaleX : float, optional
        A scaling factor for the abscissa.
    scaleY : float, optional
        A scaling factor for the ordinate.
    sampleByPlane : bool, optional
        Instead of using the cell-centre coordinates use points equally
        distributed over a plane overlayed on the
        geometry.
    planeResolution : 2-tuple, optional
        Only needed in case sampleByPlane is True. Sets the amount of
        sampling points in the x and y directions.
    plotBoundaries : bool, optional
        Whether to plot the boundary of the geometry as a black line.
    **kwargs
        Additional arguments to be passed to pyplot.quiver.

    Raises
    ------
    TypeError
        If field is neither a string or and ndarray.
    ValueError
        If the data to be plotted has less dimensions than two.
        If one or both scaling factors are non-positive.

    Returns
    -------
    Quiver
        As returned by pyplot.quiver.

    """
    pointsX = np.copy(case.cellCentres[:, 0])
    pointsY = np.copy(case.cellCentres[:, 1])

    temporaryField = None
    if isinstance(field, str):
        data = case[field]
    elif isinstance(field, (vtkmodules.numpy_interface.dataset_adapter.VTKArray,
                            np.ndarray)):
        data = np.asarray(field)
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got " + str(type(field)))

    if np.ndim(data) < 2:
        raise ValueError("The selected field appears to be a scalar!")

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")

    if sampleByPlane:
        if planeResolution is None:
            planeResolution = [50, 50]

        if not isinstance(field, str):
            temporaryField = _temporary_field_name(case)
            case[temporaryField] = field

        try:
            points, sampledData = sample_by_plane(case, planeResolution)
        finally:
            if temporaryField is not None:
                del case[temporaryField]

        pointsX = points[:, 0]
        pointsY = points[:, 1]

        if isinstance(field, str):
            data = sampledData[field]
        else:
            data = np.copy(sampledData[temporaryField])

        validPointsIdx = sampledData['vtkValidPointMask']
        data = np.ma.array(data)
        data[np.nonzero(1 - validPointsIdx), 0] = np.ma.masked
        data[np.nonzero(1 - validPointsIdx), 1] = np.ma.masked

    if normalize:
        norms = np.linalg.norm(data[:, [0, 1]], axis=1)
        # Out of place: data may alias an array owned by the caller.
        data = data/np.where(norms == 0, 1, norms)[:, np.newaxis]

    if colorField is None:
        return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                          data[:, 1], **kwargs)

    if sampleByPlane:
        if isinstance(colorField, str):
            colorData = sampledData[colorField]
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], colorData, **kwargs)
        else:
            colorData = np.ma.masked_array(colorField, mask=1 - validPointsIdx)
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], colorData, **kwargs)
    else:
        if isinstance(colorField, str):
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], case[colorField], **kwargs)
        else:
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], colorField, **kwargs)


def plot_streamlines(case, field, colorField=None,
                     scaleX=1, scaleY=1,
                     planeResolution=None,
                     plotBoundaries=True,
                     **kwargs):
    """Produce a streamline plot.

    This function wraps pyplot.streamplot. See that functions
    documentation for additional customization parameters.

    Parameters
    ----------
    case : Case
        The case that the vector field used for the streamlines belongs
        to.
    field : str or ndarray
        The vector field used for computing the streamlines. Either a
        string with the name of the field as found in the case or an
        ndarray with the data.
    colorField : string or ndarray
        Data used to colour the streamlines, either the name of a field
        in the case or an array. An array must hold one value per
        sampling point, i.e. planeResolution[0]*planeResolution[1]
        values, ordered as returned by sample_by_plane.
    scaleX : float, optional
        A scaling factor for the abscissa.
    scaleY : float, optional
        A scaling factor for the ordinate.
    planeResolution : 2-tuple, optional
        Sets the amount of sampling points in the x and y directions.
        Note that this does not control the density of the streamlines.
    plotBoundaries : bool, optional
        Whether to plot the boundary of the geometry as a black line.
    **kwargs
        Additional arguments to be passed to pyplot.streamplot.

    Raises
    ------
    TypeError
        If field is neither a string or and ndarray.
    ValueError
        If the data to be plotted has less dimensions than two.
        If one or both scaling factors are non-positive.
        If colorField is an array of the wrong length.

    Returns
    -------
    StreamPlotSet
        As returned by pyplot.streamplot.

    """
    temporaryField = None
    if isinstance(field, str):
        data = case[field]
    elif isinstance(field, (vtkmodules.numpy_interface.dataset_adapter.VTKArray,
                            np.ndarray)):
        data = np.asarray(field)
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got " + str(type(field)))

    if np.ndim(data) < 2:
        raise ValueError("The selected field appears to be a scalar!")

    if planeResolution is None:
        planeResolution = (50, 50)

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    if not isinstance(field, str):
        temporaryField = _temporary_field_name(case)
        case[temporaryField] = field

    try:
        points, sampledData = sample_by_plane(case, planeResolution)
    finally:
        if temporaryField is not None:
            del case[temporaryField]

    pointsX = points[:, 0]
    pointsY = points[:, 1]

    if isinstance(field, str):
        data = sampledData[field]
    else:
        data = np.copy(sampledData[temporaryField])

    nRows, nCols = planeResolution
    pointsX = pointsX.reshape(nCols, nRows)[:, 0]
    pointsY = pointsY.reshape(nCols, nRows)[0, :]
    dataX = data[:, 0].reshape((nRows, nCols), order='F')
    dataY = data[:, 1].reshape((nRows, nCols), order='F')

    validPoints = sampledData['vtkValidPointMask'].reshape(
        (nRows, nCols), order='F').astype(bool)
    dataX = np.ma.masked_where(~validPoints, dataX)
    dataY = np.ma.masked_where(~validPoints, dataY)

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")

    if colorField is None:
        return plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY,
                              **kwargs)
    else:
        if isinstance(colorField, str):
            colorData = sampledData[colorField]
        else:
            colorData = np.asarray(colorField)
            if colorData.shape != (nRows*nCols,):
                raise ValueError(
                    "colorField must contain one value per sampling point, "
                    f"i.e. {nRows*nCols} values, got {colorData.size}.")
        # Same ordering as the sampled vector components above.
        colorData = colorData.reshape((nRows, nCols), order='F')
        colorData = np.ma.masked_where(~validPoints, colorData)
        return plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY,
                              color=colorData, **kwargs)


def plot_field(case, field, scaleX=1, scaleY=1, xlim=None, ylim=None, plotBoundaries=True,
               colorbar=True, **kwargs):
    """Plot a field.

    This function uses a matplotlib PolyCollection to compose the
    plot. Additional customization parameters can be passed to the
    constructor of the PolyCollection via kwargs. In particular,
    cmap can be used to set the colormap and edgecolor to color the
    edges of the cells.

    Parameters
    ----------
    case : Case
        The case that the vector field used for the streamlines belongs
        to.
    field : str or ndarray
        The scalar field that will be plotted. Either a  string with
        the name of the field as found in the case or an ndarray with
        the data.
    scaleX : float, optional
        A scaling factor for the abscissa.
    scaleY : float, optional
        A scaling factor for the ordinate.
    xlim : list with two elements
        Limits for the plotted data in the x direction. Defaults to no limits.
    ylim : list with two elements
        Limits for the plotted data in the y direction. Defaults to no limits.
    plotBoundaries : bool, optional
        Whether to plot the boundary of the geometry as a black line.
    colorbar : bool, optional
        Whether to add a vertical colorbar to the right of the plot.
    **kwargs
        Additional arguments to be passed to PolyCollection constructor.

    Raises
    ------
    TypeError
        If field is neither a string or and ndarray.
    ValueError
        If the field to be plotted has more dimensions than one.
        If one or both scaling factors are non-positive.
        If xlim or ylim does not contain exactly two values.

    Returns
    -------
    PolyCollection
        The collection of polygons defining the cells.

    """
    from vtkmodules.vtkCommonDataModel import vtkBox
    from vtkmodules.vtkFiltersCore import vtkClipPolyData

    xlim = case.xlim if xlim is None else _axis_limits(xlim, "xlim")
    ylim = case.ylim if ylim is None else _axis_limits(ylim, "ylim")

    if isinstance(field, str):
        fieldName = field
        data = case[field]
        plotData = case.vtkData.VTKObject
    elif isinstance(field, (vtkmodules.numpy_interface.dataset_adapter.VTKArray,
                            np.ndarray)):
        fieldName = "__turbulucid_plot_data__"
        data = np.asarray(field)
        if data.ndim == 0 or data.shape[0] != case.vtkData.GetNumberOfCells():
            raise ValueError("The dimensionality of the provided field "
                             "does not match that of the case.")
        plotData = case.vtkData.VTKObject.NewInstance()
        plotData.ShallowCopy(case.vtkData.VTKObject)
        dsa.WrapDataObject(plotData).CellData.append(data, fieldName)
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got " + str(type(field)))

    if np.ndim(data) > 1:
        raise ValueError("The selected field appears to not be a scalar!")

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    # init to case.vtkData in case we do not need clipping
    clippedData = dsa.WrapDataObject(plotData)

    if np.any(xlim - case.xlim) or np.any(ylim - case.ylim):
        clipper = vtkClipPolyData()
        box = vtkBox()
        box.SetBounds(xlim[0], xlim[1], ylim[0], ylim[1], -1, 1)
        clipper.SetClipFunction(box)
        clipper.SetInputData(plotData)
        clipper.SetInsideOut(1)
        clipper.Update()
        clippedData = dsa.WrapDataObject(clipper.GetOutput())

    polys = _polygon_vertices(clippedData, scaleX, scaleY)

    polyCollection = PolyCollection(polys, **kwargs)

    if "edgecolor" not in kwargs:
        polyCollection.set_edgecolor("face")
    data = np.copy(vtk_to_numpy(clippedData.GetCellData()[fieldName]))
    polyCollection.set_array(data)

    ax = plt.gca()
    ax.add_collection(polyCollection, autolim=False)
    dataBounds = _as_polydata(clippedData).GetBounds()
    _update_datalim(ax, np.array([
        [dataBounds[0]/scaleX, dataBounds[2]/scaleY],
        [dataBounds[1]/scaleX, dataBounds[3]/scaleY],
    ]))

    if colorbar:
        add_colorbar(polyCollection)

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY)

    ax.set_xlim(xlim/scaleX)
    ax.set_ylim(ylim/scaleY)
    ax.set_aspect('equal')

    return polyCollection


def plot_contour(case, field, value, scaleX=1, scaleY=1, **kwargs):
    """Plot a contour plot of a scalar field.

    The cell data is first interpolated to points in order to
    use vtkContourFilter to extract the contour line.
    The kwargs are passed to the constructor of a LineCollection
    and can be used to customize the plotted line, e.g. its colour.

    Parameters
    ----------
    case : Case
        The case to draw the boundaries for.
    field : string
        The field to extract the contour from.
    value : float
        The value associated with the contour.
    scaleX : float, optional
        A scaling factor for the abscissa.
    scaleY : float, optional
        A scaling factor for the ordinate.
    **kwargs
        Additional options to pass to pyplot.tricontour.

    Raises
    ------
    TypeError
        If ``field`` is not a string or ``value`` is not a real scalar.
    ValueError
        If one or both scaling factors are non-positive.
        If the field is missing or non-scalar, or the value is non-finite.

    Returns
    -------
    LineCollection
        Collection of line segments defining the contour line.

    """
    from vtkmodules.vtkFiltersCore import vtkContourFilter

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    value = _validate_contour_request(case, field, value, "plot_contour")

    toPoint = vtkCellDataToPointData()
    toPoint.SetInputData(case.vtkData.VTKObject)
    toPoint.Update()
    pointData = toPoint.GetOutput()
    pointData.GetPointData().SetActiveScalars(field)

    contour = vtkContourFilter()
    contour.SetInputData(pointData)
    contour.SetValue(0, value)
    contour.Update()
    contour = dsa.WrapDataObject(contour.GetOutput())

    segments = _line_segments(contour, scaleX, scaleY)
    return _add_line_collection(segments, case, scaleX, scaleY, kwargs)
