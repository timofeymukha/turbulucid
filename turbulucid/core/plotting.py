from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from matplotlib import tri
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

__all__ = ["plot_boundaries", "plot_vectors", "plot_streamlines"]


def plot_boundaries(case, scaleX=1, scaleY=1, **kwargs):
    """Plot the boundaries the domain.

    The function defines a field that is 1 at the boundary and zero
    elsewhere. A contour plot of this field is the produced using
    pyplot.tricontour. See the documentation of this function for
    additional customization parameters.

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
            If on or both scaling factors are non-positive.

    """
    if (scaleX <= 0) or (scaleY <= 0):
        raise ValueError("Scaling factors must be positive.")

    pointsX = np.copy(case.cellCentres[:, 0])
    pointsY = np.copy(case.cellCentres[:, 1])

    nEdgePoints = 0
    for i in case.boundaries:
        points = case.boundary_data(i)[0]
        nEdgePoints += points.shape[0]
        pointsX = np.append(pointsX, points[:, 0])
        pointsY = np.append(pointsY, points[:, 1])

    triang = tri.Triangulation(pointsX/scaleX, pointsY/scaleY)
    z = np.zeros(pointsX.shape[0])
    z[-nEdgePoints:] = 1

    plt.tricontour(triang, z, levels=[1], **kwargs)


def plot_vectors(case, field, color=None,
                 normalize=False, scaleX=1, scaleY=1,
                 sampleByPlane=False, planeResolution=None,
                 plotBoundaries=True,
                 **kwargs):
    """Plot a vector field.

    This function wraps pyplot.quiver. See that function's
    documentation for additional customisation parameters.

    Parameters
    ----------
        case : Case
            The case that the vector field belongs to.
        field : str or ndarray
            Either a string with the name of the field as found in the
            case or an ndarray with the data.
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
            Instead of using the cell-centre coordinates use points
            equally distributed over a plane overlayed on the
            geometry.
        planeResolution : 2-tuple, optional
            Only needed in case sampleByPlane is True. Sets the
            amount of sampling points in the x and y directions.
        **kwargs
            Additional arguments to be passed to pyplot.quiver.

    Raises
    ------
        TypeError
            If field is neither a string or and ndarray.
        ValueError
            If the data to be plotted has less dimensions than two.
            If one or both scaling factors are non-positive.

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

    if sampleByPlane:
        plane = vtk.vtkPlaneSource()
        if planeResolution is not None:
            plane.SetResolution(planeResolution[0], planeResolution[1])
        else:
            plane.SetResolution(50, 50)

        minX = np.min(case.cellCentres[:, 0])
        maxX = np.max(case.cellCentres[:, 0])
        minY = np.min(case.cellCentres[:, 1])
        maxY = np.max(case.cellCentres[:, 1])

        plane.SetOrigin(minX, minY, case.zValue)
        plane.SetPoint1(minX, maxY, case.zValue)
        plane.SetPoint2(maxX, minY, case.zValue)
        plane.Update()

        probeFilter = vtk.vtkProbeFilter()
        probeFilter.SetSourceData(case.vtkData.VTKObject)
        probeFilter.SetInputConnection(plane.GetOutputPort())
        probeFilter.Update()

        probeData = dsa.WrapDataObject(probeFilter.GetOutput())
        pointsX = probeData.Points[:, 0]
        pointsY = probeData.Points[:, 1]
        if type(field) == str:
            data = probeData.PointData[field]
        else:
            data = np.copy(probeData.PointData['temp'])
            case.__delitem__('temp')

        validPointsIdx = probeData.PointData['vtkValidPointMask']
        validPointsIdx = np.nonzero(validPointsIdx)
        pointsX = pointsX[validPointsIdx]
        pointsY = pointsY[validPointsIdx]
        data = data[validPointsIdx]

    if normalize:
        norms = np.linalg.norm(data, axis=1)
        idx = np.nonzero(norms)[0]
        for i in range(data.shape[1]):
            data[idx, i] /= norms

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")


    if (color == None) or sampleByPlane:
        plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0], data[:, 1],
                   **kwargs)
    else:
        plt.quiver(pointsX/scaleX, pointsY/scaleY, data[:, 0], data[:, 1],
                   color, **kwargs)


def plot_streamlines(case, field, color=None,
                     scaleX=1, scaleY=1,
                     planeResolution=None,
                     plotBoundaries=True,
                     **kwargs):

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
    plane.SetResolution(planeResolution[0], planeResolution[1])

    plane.SetOrigin(case.bounds[0], case.bounds[2], case.zValue)
    plane.SetPoint1(case.bounds[0], case.bounds[3], case.zValue)
    plane.SetPoint2(case.bounds[1], case.bounds[2], case.zValue)
    plane.Update()

    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(case.vtkData.VTKObject)
    probeFilter.SetInputConnection(plane.GetOutputPort())
    probeFilter.Update()

    probeData = dsa.WrapDataObject(probeFilter.GetOutput())
    pointsX = probeData.Points[:, 0]
    pointsY = probeData.Points[:, 1]

    if type(field) == str:
        data = probeData.PointData[field]
    else:
        data = np.copy(probeData.PointData['temp'])
        case.__delitem__('temp')

    pointsX = pointsX.reshape(planeResolution[1]+1, planeResolution[0]+1)[:, 0]
    pointsY = pointsY.reshape(planeResolution[1]+1, planeResolution[0]+1)[0, :]
    dataX = data[:, 0].reshape((planeResolution[0]+1, planeResolution[1]+1),
                               order='F')
    dataY = data[:, 1].reshape((planeResolution[0]+1, planeResolution[1]+1),
                               order='F')

    if plotBoundaries:
        plot_boundaries(case, scaleX=scaleX, scaleY=scaleY, colors="Black")

    if color == None:
        plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY, **kwargs)
    else:
        plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY, **kwargs)
        #color = color.reshape((planeResolution[0]+1, planeResolution[1]+1),
        #                       order='F')
        #plt.streamplot(pointsX/scaleX, pointsY/scaleY, dataX, dataY,
        #               color=color, **kwargs)

