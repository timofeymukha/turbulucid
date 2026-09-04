# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from collections import OrderedDict
from numbers import Integral, Real

import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa

try:
    from vtkmodules.vtkFiltersGeneral import vtkCellCenters
except ImportError:
    from vtkmodules.vtkFiltersCore import vtkCellCenters

__all__ = ["profile_along_line", "tangents", "normals",
           "dist", "sort_indices", "sample_by_plane", "edge_lengths",
           "isoline"]


def _point_2d(point, name):
    try:
        point = np.asarray(point, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must contain two real coordinates.") from error
    if point.shape != (2,):
        raise ValueError(f"{name} must contain exactly two coordinates.")
    if not np.all(np.isfinite(point)):
        raise ValueError(f"{name} coordinates must be finite.")
    return point


def _boolean(value, name):
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a boolean.")
    return bool(value)


def _plane_resolution(resolution):
    try:
        values = tuple(resolution)
    except TypeError as error:
        raise TypeError("resolution must contain two integer values.") from error
    if len(values) != 2:
        raise ValueError("resolution must contain exactly two values.")
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
           for value in values):
        raise TypeError("resolution values must be integers.")
    if any(value < 2 for value in values):
        raise ValueError("resolution values must be at least 2.")
    return tuple(int(value) for value in values)


def _validate_contour_request(case, field, value, caller):
    """Validate a request for a contour of a scalar field.

    Shared by :func:`isoline` and :func:`turbulucid.plot_contour` so that
    both reject a missing or non-scalar field instead of quietly producing
    an empty contour.

    Parameters
    ----------
    case : Case
        The case the field belongs to.
    field : str
        The name of the field to contour.
    value : float
        The value associated with the contour.
    caller : str
        The name of the calling function, used in the error message.

    Returns
    -------
    float
        The validated contour value.

    Raises
    ------
    TypeError
        If ``field`` is not a string or ``value`` is not a real scalar.
    ValueError
        If the field is missing or non-scalar, or the value is non-finite.

    """
    if not isinstance(field, str):
        raise TypeError("field must be a string.")
    if field not in case.fields:
        raise ValueError(f"Field {field} not present in the case.")
    if case[field].ndim != 1:
        raise ValueError(f"{caller} requires a scalar field.")
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError("value must be a real scalar.")
    value = float(value)
    if not np.isfinite(value):
        raise ValueError("value must be finite.")
    return value


def profile_along_line(case, p1, p2, correctDistance=False,
                       excludeBoundaries=False):
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

    Returns
    -------
    (ndarray, dictionary)
        The first element of the tuple is an array of coordinates along
        the line where the data is located. The second element is a
        dictionary with the case's fields as keys and arrays of values
        of these fields as values.

    Raises
    ------
    TypeError
        If coordinates or boolean options have incompatible types.
    ValueError
        If either point is invalid or the line has zero length.

    """
    from vtkmodules.vtkCommonDataModel import vtkPlane
    from vtkmodules.vtkFiltersCore import vtkCutter

    p1 = _point_2d(p1, "p1")
    p2 = _point_2d(p2, "p2")
    correctDistance = _boolean(correctDistance, "correctDistance")
    excludeBoundaries = _boolean(excludeBoundaries, "excludeBoundaries")
    if np.array_equal(p1, p2):
        raise ValueError("p1 and p2 must define a nonzero-length line.")
    if case.vtkData.GetNumberOfPoints() == 0:
        raise ValueError("The case contains no geometry to sample.")

    p1 = np.append(p1, case.vtkData.Points[0, 2])
    p2 = np.append(p2, case.vtkData.Points[0, 2])

    # Compute the plane-normal as a cross-product
    unit = (p2 - p1)/np.linalg.norm(p2 - p1)
    geomNormal = np.array([0, 0, 1])
    planeNormal = np.cross(unit, geomNormal)

    # Define the cutting plane
    plane = vtkPlane()
    plane.SetNormal(planeNormal[0], planeNormal[1], planeNormal[2])
    plane.SetOrigin(p1[0], p1[1], p1[2])

    # Create the cutter and extract the sampled data
    planeCut = vtkCutter()
    planeCut.SetInputData(case.vtkData.VTKObject)
    planeCut.GenerateTrianglesOff()
    planeCut.SetCutFunction(plane)
    planeCut.Update()
    cutData = dsa.WrapDataObject(planeCut.GetOutput())

    # Create the cell-centres on the cut -- that is where the data is
    cCenters = vtkCellCenters()
    cCenters.SetInputData(planeCut.GetOutput())
    cCenters.Update()
    cCentersData = dsa.WrapDataObject(cCenters.GetOutput())

    # Grab the data and its coordinates
    centrePoints = dsa.WrapDataObject(cCenters.GetOutput()).Points
    if centrePoints is None:
        coords = np.empty((0, 3))
    else:
        coords = np.array(centrePoints)
    data = OrderedDict()

    for field in cCentersData.PointData.keys():
        data[field] = cCentersData.PointData[field]

    # Extract data from the boundaries
    if not excludeBoundaries:
        boundaryCC = vtkCellCenters()
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

    lineLength = np.linalg.norm(p2 - p1)

    for i in range(coords.shape[0]):
        distP1 = np.linalg.norm(coords[i, :] - p1)
        distP2 = np.linalg.norm(coords[i, :] - p2)
        if (distP1 <= lineLength) and (distP2 <= lineLength):
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

    # Correct the distance values
    if correctDistance:
        # Find the point (not cell-center!) closest to p1, get correction
        planeCut.SetInputData(case.vtkData.VTKObject)
        planeCut.Update()

        if cutData.Points is None or cutData.GetNumberOfPoints() == 0:
            raise ValueError("The line does not intersect the case geometry.")
        shiftPointId = cutData.VTKObject.FindPoint(p1)
        if shiftPointId < 0:
            raise ValueError("Could not locate the line intersection point.")
        shiftPoint = cutData.Points[shiftPointId, :]
        correction = np.linalg.norm(shiftPoint - p1)
        distance -= correction

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

    Raises
    ------
    ValueError
        If the named boundary contains a zero-length edge.

    """
    block = case.extract_block_by_name(name)
    nCells = block.GetNumberOfCells()

    tangents = np.zeros([nCells, 2])

    for cellI in range(nCells):
        cell = block.GetCell(cellI)
        point0 = np.array(cell.GetPoints().GetPoint(0))[:2]
        point1 = np.array(cell.GetPoints().GetPoint(1))[:2]
        edge = point1 - point0
        length = np.linalg.norm(edge)
        if length == 0:
            raise ValueError(f"Boundary {name} contains a zero-length edge.")
        tangents[cellI, :] = edge/length

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

    boundaryCoords = case.boundary_data(name)[0]
    cellCoords = case.boundary_cell_data(name)[0]
    d = cellCoords - boundaryCoords

    for i in range(n.shape[0]):

        n[i] = np.cross(np.append(t[i, :], 0), z)[:2]

        # Ensure normal is outward
        if np.dot(n[i], d[i]) > 0:
            n[i] = -n[i]

    return n


def dist(case, name, corrected=True, sort=None):
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
    sort : {None, 'x', 'y'}, optional
        Whether to sort the data along a coordinate. Use 'x' and
        'y' to sort along x and y, respectively. Default is no
        sorting.

    Returns
    -------
    ndarray
        The value of the distance for each adjacent cell.

    Raises
    ------
    TypeError
        If ``corrected`` is not a boolean.

    """
    corrected = _boolean(corrected, "corrected")
    boundaryCoords = case.boundary_data(name)[0]
    cellCoords = case.boundary_cell_data(name)[0]

    if sort is not None:
        idx = sort_indices(case, name, sort)
    else:
        idx = np.arange(0, boundaryCoords.shape[0], 1, dtype=np.int32)

    d = cellCoords - boundaryCoords
    if not corrected:
        return np.linalg.norm(d, axis=1)[idx]
    else:
        n = normals(case, name)
        dNormal = np.zeros(d.shape[0])

        for i in range(dNormal.size):
            dNormal[i] = np.linalg.norm(np.dot(d[i, :], n[i, :])*n[i, :])

        return dNormal[idx]


def edge_lengths(case, name):
    """Compute the lengths of boundary edges.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    name : str
        The name of the boundary.

    Returns
    -------
    ndarray
        The lengths of each edge of the boundary.

    """
    sizes = np.zeros(case.boundary_data(name)[0].shape[0])

    block = case.extract_block_by_name(name)
    for c in range(block.GetNumberOfCells()):
        point0 = block.GetCell(c).GetPoints().GetPoint(0)[:2]
        point1 = block.GetCell(c).GetPoints().GetPoint(1)[:2]

        sizes[c] = np.sqrt((point1[0] - point0[0])**2 +
                           (point1[1] - point0[1])**2)

    return sizes


def sort_indices(case, name, axis):
    """Compute indices sorting the values on a boundary along an axis.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    name : str
        The name of the boundary.
    axis : {'x', 'y'}
        Axis to sort along.

    Returns
    -------
    ndarray
        The  array of indices.

    """

    if axis not in {"x", "y"}:
        raise ValueError("axis should be x or y.")

    blockData = case.extract_block_by_name(name)

    cCenters = vtkCellCenters()
    cCenters.SetInputData(blockData)
    cCenters.Update()

    points = np.array(dsa.WrapDataObject(cCenters.GetOutput()).Points)

    if axis == "x":
        return np.argsort(points[:, 0])
    return np.argsort(points[:, 1])


def sample_by_plane(case, resolution):
    """Sample the field values by a Cartesian grid.

    Parameters
    ----------
    case : Case
        The case to extract data from.
    resolution : pair
        The resolution of the plane, (number of rows, number of cols).

    Returns
    -------
    (ndarray, dictionary)
        The first element of the tuple is an array of coordinates of
        the points on the plane. The second is a dictionary with the
        case's fields as keys and arrays of values of these fields as
        values.

    Raises
    ------
    TypeError
        If the resolution values are not integers.
    ValueError
        If resolution does not contain two values of at least two.

    """
    from vtkmodules.vtkFiltersSources import vtkPlaneSource
    from vtkmodules.vtkFiltersCore import vtkProbeFilter

    resolution = _plane_resolution(resolution)
    plane = vtkPlaneSource()
    plane.SetResolution(resolution[0] - 1, resolution[1] - 1)

    smallDy = (case.bounds[3] - case.bounds[2])/10000
    smallDx = (case.bounds[1] - case.bounds[0])/10000

    plane.SetOrigin(case.bounds[0] - smallDx, case.bounds[2] - smallDy, 0)
    plane.SetPoint1(case.bounds[0] - smallDx, case.bounds[3] + smallDy, 0)
    plane.SetPoint2(case.bounds[1] + smallDx, case.bounds[2] - smallDy, 0)

    plane.Update()

    probeFilter = vtkProbeFilter()
    probeFilter.SetSourceData(case.vtkData.VTKObject)
    probeFilter.SetInputConnection(plane.GetOutputPort())
    probeFilter.Update()

    probeData = dsa.WrapDataObject(probeFilter.GetOutput())
    points = probeData.Points[:, [0, 1]]

    # validPointsIdx = probeData.PointData['vtkValidPointMask']
    # validPointsIdx = np.nonzero(validPointsIdx)
    # points = points[validPointsIdx, :]

    data = {}
    for key in probeData.PointData.keys():
        if probeData.PointData[key].ndim == 1:
            # data[key] = np.array(probeData.PointData[key][validPointsIdx])
            data[key] = np.array(probeData.PointData[key])
        else:
            # data[key] = np.array(probeData.PointData[key][validPointsIdx, :])
            data[key] = np.array(probeData.PointData[key])

    return points, data


def isoline(case, field, value):
    """Extract an isoline of a scalar field.

    The cell data is first interpolated to points in order to
    use vtkContourFilter to extract the isoline.

    Parameters
    ----------
    case : Case
        The case to draw the boundaries for.
    field : string
        The field to extract the contour from.
    value : float
        The value associated with the contour.

    Returns
    -------
    ndarray
        Points defining the isoline.

    Raises
    ------
    TypeError
        If ``field`` is not a string or ``value`` is not a real scalar.
    ValueError
        If the field is missing or non-scalar, or the value is non-finite.

    """
    from vtkmodules.vtkFiltersCore import vtkCellDataToPointData, vtkContourFilter

    value = _validate_contour_request(case, field, value, "isoline")

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

    if contour.Points is None:
        return np.empty((0, 2))
    return np.array(contour.Points[:, :2])
