from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
from scipy.interpolate import interp1d

__all__ = ["profile_along_gridline", "interpolate_dataset"]


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
