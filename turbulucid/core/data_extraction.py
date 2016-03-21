import numpy as np

__all__ = ["profile_along_gridline"]


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
                print "Note: using point "+str(actual_x)+" instead of "+str(point)
        else:
            idx = index
        values = V[:,idx]
        coords = Y[:,idx]
    elif (direction == "x"):
        if( index == -1):
            idx = np.argmin(abs(Y[:,0]-point))
            actual_y = Y[idx,0]

            if (actual_y != point):
                print "Note: using point = "+str(actual_y)+" instead of "+str(point)
        else:
            idx = index
        values = V[idx,:]
        coords = X[idx,:]

    return [coords, values]

