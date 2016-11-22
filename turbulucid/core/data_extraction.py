from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from scipy.interpolate import interp1d

__all__ = ["profile_along_gridline", "interpolate_dataset", "profile_along_line"]


def profile_along_gridline(field, point, direction="y", index=-1):
    """Return the values along a path following a grid-line 
    Equivalent to a line-profile on a rectangular mesh.
    
    Parameters
    ----------
        field -- the field as a dictionary with X, Y and Z
        point -- the point from where to start the path 
        direction -- which direction to move in (x or y), default is "y"
        index -- instead of point, choose the array index directly
    """

    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    if (direction == "y"):
        if( index == -1):
            idx = np.argmin(abs(X[0,:]-point))
            actual_x = X[0,idx]

            if (actual_x != point):
                print("Note: using point "+str(actual_x)+" instead of "+str(point))
        else:
            idx = index
        values = V[:,idx]
        coords = Y[:,idx]
    elif (direction == "x"):
        if( index == -1):
            idx = np.argmin(abs(Y[:,0]-point))
            actual_y = Y[idx,0]

            if (actual_y != point):
                print("Note: using point = "+str(actual_y)+" instead of "+str(point))
        else:
            idx = index
        values = V[idx,:]
        coords = X[idx,:]

    return [coords, values]


def profile_along_line(case, field, startPoint, endPoint):
    # Convert point to ndarrays
    startPoint = np.array(startPoint)
    endPoint = np.array(endPoint)

    # Create a locator to check which cells are intersected by the line
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(case.reader.GetOutput())
    locator.BuildLocator()
    idList = vtk.vtkIdList()
    locator.FindCellsAlongLine(startPoint, endPoint, 0.0, idList)

    nCells = idList.GetNumberOfIds()
    nDims = case.vtkData.CellData[field][0].size

    # unit vector point along the line
    direction = (endPoint - startPoint)/np.linalg.norm(endPoint-startPoint)

    coords = np.zeros((nCells, 3))
    distance = np.zeros((nCells, 1))
    data = np.zeros((nCells, nDims))

    for i in range(nCells):
        coords[i, :] = case.cellCentres[idList.GetId(i)]
        data[i, :] = case.vtkData.CellData[field][idList.GetId(i)]

        # Project the vector connecting the starting point and the cell-center
        # to get the distance
        distance[i] = np.sum((coords[i,:] - startPoint)*direction)
    return [coords, distance, data]


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
