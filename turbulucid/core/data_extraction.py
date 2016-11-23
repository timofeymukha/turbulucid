from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
from scipy.interpolate import interp1d

__all__ = ["interpolate_dataset", "profile_along_line"]


def profile_along_line(case, p1, p2):
    # Convert point to ndarrays
    p1 = np.array(p1)
    p2 = np.array(p2)

    unit = (p2 - p1)/np.linalg.norm(p2 - p1)

    geomNormal = np.array([0, 0, 1])
    planeNormal = np.cross(unit, geomNormal)

    plane = vtk.vtkPlane()
    plane.SetNormal(planeNormal[0], planeNormal[1], planeNormal[2])
    plane.SetOrigin(p1[0], p1[1], p1[2])

    planeCut = vtk.vtkCutter()
    planeCut.SetInputData(case.reader.GetOutput())
    planeCut.GenerateTrianglesOff()
    planeCut.SetCutFunction(plane)
    planeCut.Update()
    cutData = dsa.WrapDataObject(planeCut.GetOutput())

    cCenters = vtk.vtkCellCenters()
    cCenters.SetInputData(planeCut.GetOutput())
    cCenters.Update()

    coords = dsa.WrapDataObject(cCenters.GetOutput()).Points
    data = cutData.CellData

    nCells = cutData.GetNumberOfCells()

    distance = np.zeros(nCells)

    for i in range(nCells):
        distance[i] = np.linalg.norm(coords[i, :] - p1)
    return [distance, data, coords]


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
