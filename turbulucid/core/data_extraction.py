from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import interp1d
from collections import OrderedDict

__all__ = ["interpolate_dataset", "profile_along_line", "tangents", "normals",
           "dist"]


def profile_along_line(case, p1, p2, correctDistance=False,
                       excludeBoundaries=False, debug=False):
    """Extract a linear profile from the data.

    Returns the data itself and the its location in a coordinate
    system tangential to the line.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    p1 : 2-tuple
        Start point defining the line.
    p2 : 2-tuple
        End point defining the line.
    correctDistance : bool
        Whether to place the origo of the coordinate system at p1
        or at the point where the the line first intersects the
        geometry.
    excludeBoundaries : bool
        Exclude boundary data when extracting the profile, default
        is False.
    debug : bool
        Debug switch, if True, additional output is given.

    """
    # Convert point to ndarrays
    zValue = case.zValue

    # Add the z-value to the points.
    p1 = np.append(np.array(p1), zValue)
    p2 = np.append(np.array(p2), zValue)

    # Compute the plane-normal as a cross-product
    unit = (p2 - p1)/np.linalg.norm(p2 - p1)
    geomNormal = np.array([0, 0, 1])
    planeNormal = np.cross(unit, geomNormal)

    # Define the cutting plane
    plane = vtk.vtkPlane()
    plane.SetNormal(planeNormal[0], planeNormal[1], planeNormal[2])
    plane.SetOrigin(p1[0], p1[1], p1[2])

    # Create the cutter and extract the sampled data
    planeCut = vtk.vtkCutter()
    planeCut.SetInputData(case.vtkData.VTKObject)
    planeCut.GenerateTrianglesOff()
    planeCut.SetCutFunction(plane)
    planeCut.Update()
    cutData = dsa.WrapDataObject(planeCut.GetOutput())

    # Create the cell-centres on the cut -- that is where the data is
    cCenters = vtk.vtkCellCenters()
    cCenters.SetInputData(planeCut.GetOutput())
    cCenters.Update()
    cCentersData = dsa.WrapDataObject(cCenters.GetOutput())

    # Grab the data and its coordinates
    coords = np.copy(dsa.WrapDataObject(cCenters.GetOutput()).Points)
    data = OrderedDict()

    for field in cCentersData.PointData.keys():
        data[field] = cCentersData.PointData[field]

    # Extract data from the boundaries
    if not excludeBoundaries:
        boundaryCC = vtk.vtkCellCenters()
        for boundary in case.boundaries:

            block = case.extract_block_by_name(boundary)
            planeCut.SetInputData(block)
            planeCut.Update()

            boundaryCC.SetInputData(planeCut.GetOutput())
            boundaryCC.Update()
            boundaryCCData = dsa.WrapDataObject(boundaryCC.GetOutput())

            if boundaryCCData.Points is not None:
                coords = np.append(coords, boundaryCCData.Points, axis=0)

                for field in data.keys():
                    data[field] = np.append(data[field],
                                            boundaryCCData.PointData[field],
                                            axis=0)

    validIds = []

    l = np.linalg.norm(p2 - p1)

    for i in range(coords.shape[0]):
        distP1 = np.linalg.norm(coords[i, :] - p1)
        distP2 = np.linalg.norm(coords[i, :] - p2)
        if (distP1 <= l) and (distP2 <= l):
            validIds.append(i)

    coords = coords[validIds, :]
    for field in data.keys():
        if np.ndim(data[field]) == 1:
            data[field] = data[field][validIds]
        else:
            data[field] = data[field][validIds, :]

    nCells = coords.shape[0]
    distance = np.zeros(nCells)

    # Compute the distance from p1 to each coordinate
    for i in range(nCells):
        distance[i] = np.linalg.norm(coords[i, :] - p1)

    # Sort
    order = np.argsort(distance)
    distance = distance[order]

    # Sort the data and convert to numpy arrays
    dataNumpy = {}
    for key in data.keys():
        if np.ndim(data[key]) > 1:
            dataNumpy[key] = np.array(data[key])[order, :]
        else:
            dataNumpy[key] = np.array(data[key])[order]

    # Find the point (not cell-center!) closest to p1, get correction
    planeCut.SetInputData(case.vtkData.VTKObject)
    planeCut.Update()

    shiftPointId = cutData.VTKObject.FindPoint(p1)
    shiftPoint = cutData.Points[shiftPointId, :]
    correction = np.linalg.norm(shiftPoint - p1)

    # Correct the distance values
    if correctDistance:
        distance -= correction

    if debug:
        print("p1 ", p1)
        print("p2 ", p2)
        print("Plane normal", planeNormal)
        print("Point id closest to p1 ", shiftPointId)
        print("The point closest to p1", shiftPoint)
        print("The distance correction is", correction)

        return distance, dataNumpy, coords

    return distance, dataNumpy


def tangents(case, name):
    """ Compute unit tangent vectors for a given boundary.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    name : str
        The name of the boundary.

    Returns
    -------
    ndarray
        The tangent vectors

    """
    block = case.extract_block_by_name(name)
    nCells = block.GetNumberOfCells()

    tangents = np.zeros([nCells, 3])

    for cellI in range(nCells):
        cell = block.GetCell(cellI)
        point0 = np.array(cell.GetPoints().GetPoint(0))
        point1 = np.array(cell.GetPoints().GetPoint(1))
        tangents[cellI, :] = (point1 - point0)/np.linalg.norm(point1 - point0)

    return tangents


def normals(case, name):
    """Compute outward unit normal vectors for a given boundary.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    name : str
        The name of the boundary.

    Returns
    -------
    ndarray
        The outward normal vectors

    """
    t = tangents(case, name)

    n = np.zeros(t.shape)
    z = np.array([0, 0, 1])

    for i in range(n.shape[0]):
        n[i] = np.cross(t[i, :], z)

    boundaryCoords = case.boundary_data(name)[0]
    cellCoords = case.boundary_cell_data(name)[0]
    d = cellCoords - boundaryCoords

    for i in range(n.shape[0]):
        n[i] = np.cross(t[i, :], z)
        # Ensure it is outward
        if np.dot(n[i], d[i]) > 0:
            n[i] = -n[i]

    return n


def dist(case, name, corrected=True):
    """Compute the distances between the boundary and the adjacent
    cell-centre.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    name : str
        The name of the boundary.
    corrected : bool
        If true, projects the distance between face center and cell
        center onto the wall-normal direction.

    Returns
    -------
    ndarray
        The value of the distance for each adjacent cell.

    """
    boundaryCoords = case.boundary_data(name)[0]
    cellCoords = case.boundary_cell_data(name)[0]

    d = cellCoords - boundaryCoords

    if not corrected:
        return np.linalg.norm(d, axis=1)
    else:
        n = normals(case, name)
        dNormal = np.zeros(d.shape[0])

        for i in range(dNormal.size):
            dNormal[i] = np.linalg.norm(np.dot(d[i, :], n[i, :])*n[i, :])

        return dNormal


def interpolate_dataset(dataset, value, xAxis, yAxis):
    dataFile = h5py.File(dataset)

    availVals = np.sort(np.array(map(float, dataFile.keys()), dtype=np.int64))
    if value < availVals.min() or value > availVals.max():
        raise ValueError("Error: desired value outside available interpolation"
                         " range")

    argMin = np.argmin(np.abs(availVals - value))

    if availVals[argMin] < value:
        valLeft = str(availVals[argMin])
        valRight = str(availVals[argMin+1])
    else:
        valLeft = str(availVals[argMin-1])
        valRight = str(availVals[argMin])

    interpLeft = interp1d(dataFile[valLeft][xAxis],
                          dataFile[valLeft][yAxis])

    interpRight = interp1d(dataFile[valRight][xAxis],
                           dataFile[valRight][yAxis])

    interpX = dataFile[valRight][xAxis][:]

    interpY = np.zeros(interpX.shape)

    for i in range(interpX.size):
        interpY[i] = interp1d(np.array([float(valLeft), float(valRight)]),
                              np.array([interpLeft(interpX[i]),
                                        interpRight(interpX[i])]))(value)

    dataFile.close()
    return interpX, interpY
