# This file is part of turbulucid
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from mpl_toolkits import axes_grid1
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from .data_extraction import sample_by_plane


__all__ = ["plot_boundaries", "plot_vectors", "plot_streamlines", "plot_field",
           "add_colorbar"]


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
        Additional options to pass to pyplot.tricontour.

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

    ax = plt.gca()
    segments = []
    for boundary in case.boundaries:
        block = case.extract_block_by_name(boundary)
        for c in range(block.GetNumberOfCells()):
            point0 = block.GetCell(c).GetPoints().GetPoint(0)[:2]
            point1 = block.GetCell(c).GetPoints().GetPoint(1)[:2]

            point0 = np.array(point0)/[scaleX, scaleY]
            point1 = np.array(point1)/[scaleX, scaleY]
            segments.append((point0, point1))
    collection = LineCollection(segments, **kwargs)
    if "color" not in kwargs:
        collection.set_color("Black")

    ax.add_collection(collection)
    ax.set_xlim(case.xlim/scaleX)
    ax.set_ylim(case.ylim/scaleY)
    ax.set_aspect('equal')

    return collection


def plot_vectors(case, field, color=None,
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
    color : ndarray
        Data used to colour the vectors.
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

    if type(field) == str:
        data = case[field]
    elif ((type(field) == vtk.numpy_interface.dataset_adapter.VTKArray) or
          (type(field) == np.ndarray)):
        case['temp'] = field
        data = case['temp']
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got "+str(type(field)))

    if np.ndim(data) < 2:
        raise ValueError("The selected field appears to be a scalar!")

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")

    if sampleByPlane:
        plane = vtk.vtkPlaneSource()
        if planeResolution is not None:
            plane.SetResolution(planeResolution[0], planeResolution[1])
        else:
            plane.SetResolution(50, 50)

        points, sampledData = sample_by_plane(case, planeResolution)

        pointsX = points[:, 0]
        pointsY = points[:, 1]

        if type(field) == str:
            data = sampledData[field]
        else:
            data = np.copy(sampledData['temp'])
            case.__delitem__('temp')

        validPointsIdx = sampledData['vtkValidPointMask']
        data = np.ma.array(data)
        data[np.nonzero(1 - validPointsIdx), 0] = np.ma.masked
        data[np.nonzero(1 - validPointsIdx), 1] = np.ma.masked

    if normalize:
        norms = np.linalg.norm(data, axis=1)
        for j in range(data.shape[1]):
            for i in range(data.shape[0]):
                if norms[i] != 0:
                    data[i, j] /= norms[i]

    if color is None:
        return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                          data[:, 1], **kwargs)

    if sampleByPlane:
        if type(color) == str:
            colorData = sampledData[color].reshape(planeResolution[1]+1,
                                                   planeResolution[0]+1)
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], colorData, **kwargs)
        else:
            colorData = np.ma.masked_array(color, mask=1-validPointsIdx)
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], colorData, **kwargs)
    else:
        if type(color) == str:
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], case[color], **kwargs)
        else:
            return plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0],
                              data[:, 1], color, **kwargs)


def plot_streamlines(case, field, color=None,
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
    color : ndarray
        Data used to colour the vectors.
    normalize : bool, optional
        Whether to normalize the the length of the vectors.
        Default is False.
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
    StreamPlotSet
        As returned by pyplot.streamplot.

    """
    if type(field) == str:
        data = case[field]
    elif ((type(field) == vtk.numpy_interface.dataset_adapter.VTKArray) or
          (type(field) == np.ndarray)):
        case['temp'] = field
        data = case['temp']
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got "+str(type(field)))

    if np.ndim(data) < 2:
        raise ValueError("The selected field appears to be a scalar!")

    plane = vtk.vtkPlaneSource()
    if planeResolution is None:
        planeResolution = (50, 50)

    points, sampledData = sample_by_plane(case, planeResolution)
    pointsX = points[:, 0]
    pointsY = points[:, 1]

    if type(field) == str:
        data = sampledData[field]
    else:
        data = np.copy(sampledData['temp'])
        case.__delitem__('temp')

    pointsX = pointsX.reshape(planeResolution[1]+1, planeResolution[0]+1)[:, 0]
    pointsY = pointsY.reshape(planeResolution[1]+1, planeResolution[0]+1)[0, :]
    dataX = data[:, 0].reshape((planeResolution[0]+1, planeResolution[1]+1),
                               order='F')
    dataY = data[:, 1].reshape((planeResolution[0]+1, planeResolution[1]+1),
                               order='F')

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")

    if color is None:
        return plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY,
                              **kwargs)
    else:
        if type(color) == str:
            colorData = sampledData[color].reshape(planeResolution[1]+1,
                                                   planeResolution[0]+1)
        else:
            colorData = color
        plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY,
                       color=colorData, **kwargs)


def plot_field(case, field, scaleX=1, scaleY=1, plotBoundaries=True,
               colorbar=True, **kwargs):

    """Plot a field.

    This function uses a matplotlib PatchCollection to compose the
    plot. Additional customization parameters can be passed to the
    constructor of the PatchCollection via kwargs. In particular,
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
    plotBoundaries : bool, optional
        Whether to plot the boundary of the geometry as a black line.
    colorbar : bool, optional
        Whether to add a vertical colorbar to the right of the plot.
    **kwargs
        Additional arguments to be passed to PatchCollection constructor.

    Raises
    ------
    TypeError
        If field is neither a string or and ndarray.
    ValueError
        If the field to be plotted has more dimensions than one.
        If one or both scaling factors are non-positive.

    Returns
    -------
    PatchCollection
        The collection of polygons defining the cells.

    """
    if type(field) == str:
        data = case[field]
    elif ((type(field) == vtk.numpy_interface.dataset_adapter.VTKArray) or
          (type(field) == np.ndarray)):
        case['temp'] = field
        data = case['temp']
    else:
        raise TypeError("field should be a name of an existing field or an"
                        " array of values. Got "+str(type(field)))

    if np.ndim(data) > 1:
        raise ValueError("The selected field appears to not be a scalar!")

    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    polys = []
    for i in range(case.vtkData.GetNumberOfCells()):
        cell = case.vtkData.GetCell(i)
        nPoints = cell.GetNumberOfPoints()
        points = np.zeros((nPoints, 2))
        for pointI in range(nPoints):
            points[pointI, :] = cell.GetPoints().GetPoint(pointI)[:2]
            points[pointI, :] /= [scaleX, scaleY]

        polys.append(Polygon(points))

    patchCollection = PatchCollection(polys, **kwargs)

    if "edgecolor" not in kwargs:
        patchCollection.set_edgecolor("face")
    patchCollection.set_array(data)

    ax = plt.gca()
    ax.add_collection(patchCollection)
    ax.set_xlim(case.xlim/scaleX)
    ax.set_ylim(case.ylim/scaleY)
    ax.set_aspect('equal')

    if colorbar:
        add_colorbar(patchCollection)

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY)
    return patchCollection
