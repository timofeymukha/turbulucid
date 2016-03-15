import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from os import path

def chauhan_U_inner(yPlus, kappa=0.384, a=-10.361):
    alpha = (-1/kappa - a)/2
    beta = np.sqrt(-2*a*alpha-alpha**2)
    R = np.sqrt(alpha**2 + beta**2)
    return  1/kappa*np.log(-(yPlus-a)/a) + R**2/(a*(4*alpha-a))* \
           ((4*alpha+a)*np.log(-a/R*np.sqrt((yPlus-alpha)**2 + beta**2)/ \
           (yPlus - a)) + alpha/beta*(4*alpha + 5*a)*(np.arctan((yPlus-alpha)/beta) \
           + np.arctan(alpha/beta)))    

def chauhan_U_inner_mod(yPlus, kappa=0.384, a=-10.361):
    return chauhan_U_inner(yPlus, kappa, a) + 1/2.85*np.exp(-np.log(yPlus/30)**2)

def chauhan_wake(eta, Pi, a2=132.8410, a3=-166.2041, a4=71.9114):
    nom1 = 1 - np.exp(-0.25*(5*a2 + 6*a3 + 7*a4)*eta**4 + a2*eta**5 + a3*eta**6 + a4*eta**7)
    nom2 = 1 - 0.5/Pi*np.log(eta)
    denom = 1 - np.exp(-0.25*(a2+2*a3+3*a4))
    return nom1*nom2/denom

def chauhan_U_composite(yPlus, eta, Pi, kappa=0.384, a=-10.361):
    return chauhan_U_inner_mod(yPlus, kappa, a) + 2*Pi/kappa*chauhan_wake(eta, Pi)

def epsilon_ReT(y, ReTau, kappa, APlus):
    return (0.5*np.sqrt(1 + kappa**2*ReTau**2/9*(1 - (y - 1)**2)**2*(1 + 2*(y-1)**2)**2*(1 - np.exp(-y*ReTau/APlus))**2) - 0.5)

def zpgtbl_momentum_thickness(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    x = X[0,:]
    theta = np.zeros(x.size) 
    for i in range(x.size):
        yIdxTop = np.argmin(np.abs(Y[:,i]-yLoc[i]))
        y = Y[0:yIdxTop,i]
        theta[i] = np.trapz(V[0:yIdxTop,i]/u0[i]*(1-V[0:yIdxTop,i]/u0[i]), y)     

    return [x, theta]    

def zpgtbl_delta_star(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    x = X[0,:]

    deltaStar = np.zeros(x.size) 

    for i in xrange(x.size):
        yIdxTop = np.argmin(np.abs(Y[:,i]-yLoc[i]))
        y =Y[0:yIdxTop,i]
        deltaStar[i] = np.trapz((1-V[0:yIdxTop,i]/u0[i]), y)


    return [x, deltaStar]

def zpgtbl_delta_99(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    x = X[0,:]

    delta99 = np.zeros(x.shape)

    for i in xrange(x.size):
        #interp = interp1d(Y[:,i], V[:,i])
        #newY = np.linspace(Y[0,i], yLoc[i], 10000)
        #newV = interp(newY)
        for j in range(V[:,i].size):
            if V[j, i] >= 0.99*u0[i]:
                delta99[i] = Y[j, i]
                break
    return delta99


def plot_Cp(field, pRefPoint, wall, uRef, **kwargs):

    if (wall == "top"):
        [x,v] = profile_along_gridline(field, 10, direction='x')
    else: 
        [x,v] = profile_along_gridline(field, 0, direction='x')
    
    v = v-v[np.argmin(abs(x[:]-pRefPoint))]
    uRef = uRef[np.argmin(abs(x[:]-pRefPoint))]
    print uRef

    plt.plot(x/l_ramp,v/(0.5*uRef**2), **kwargs)
    plt.xlim([-2,8])
    plt.xlabel('x/l')
    plt.ylabel('Cp')
    plt.legend(loc=4)
    plt.grid()

def delta_99(field):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    x = X[0,:]
    delta99 = np.zeros(x.size) 
    u0 = np.zeros(x.size) 
    for i in range(x.size):
        # Index of the point in the middle
        yMidIdx = np.argmin(abs(Y[:,i]-0.5*(Y[0,i]+Y[-1,i]))) 
        # Freestream velocity
        u0 = V[yMidIdx,i]
        delta99[i] = Y[np.argmin(abs((u0-V[0:yMidIdx,i])/u0-0.01)),i]     

    return [x, delta99]    

def delta_star_top(field, relDiff=0.001):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    x = X[0,:]
    deltaStar = np.zeros(x.size) 
    deltaStarOld = np.copy(deltaStar)

    u0 = np.zeros(x.size) 
    convY = np.zeros(x.size) 
    yLoc = np.zeros(x.size) 

    for i in xrange(x.size):
        maxYidx = np.min(np.argwhere(V[:,i] == np.max(V[:,i])))

        for yIdx in xrange(Y.shape[0]-10,maxYidx,-1):
            y = 1000*Y[yIdx:Y.shape[0],i]
            u0[i] = V[yIdx,i]
            deltaStar[i] = np.trapz((1-V[yIdx:Y.shape[0],i]/u0[i]), y)
            yLoc[i] = Y[yIdx,i]

            refDiffComp = np.abs(deltaStar[i]-deltaStarOld[i])/deltaStarOld[i]
            if ( refDiffComp < relDiff):
                break

            deltaStarOld[i] = deltaStar[i]

    return [x, deltaStar, u0, yLoc]    

def delta_star_bottom(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    x = X[0,:]

    deltaStar = np.zeros(x.size) 

    for i in xrange(x.size):
        yIdxTop = np.argmin(np.abs(Y[:,i]-yLoc[i]))
        y = 1000*Y[0:yIdxTop,i]
        deltaStar[i] = np.trapz((1-V[0:yIdxTop,i]/u0[i]), y)


    return [x, deltaStar]


def theta(field):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    x = X[0,:]
    theta = np.zeros(x.size) 
    u0 = np.zeros(x.size) 
    for i in range(x.size):
        # Index of the point in the middle
        yMidIdx = np.argmin(abs(Y[:,i]-0.5*(Y[0,i]+Y[-1,i]))) 
        # Freestream velocity
        u0[i] = V[yMidIdx,i]
        y=1000*(Y[0:yMidIdx,i]-Y[0,i])
        theta[i] = np.trapz(V[0:yMidIdx,i]/u0[i]*(1-V[0:yMidIdx,i]/u0[i]), y)     

    return [x, theta]    

def Re_theta(field, l, nu):
    
    [x, fieldTheta] = theta(field)

    return [x, fieldTheta*l/nu]    

def u_0(field, wall="bottom", relDiff=0.01):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    u0Top = delta_star_top(field)[2]
    yLocTop = delta_star_top(field)[3]
    x = X[0,:]

    if (wall == "top"):
        u0 = u0Top
        yLoc = yLocTop
    
    else:
        u0 = np.zeros(x.size)
        yLoc = np.zeros(x.size)

        for i in xrange(x.size):
            yIdxTop = np.argmin(np.abs(Y[:,i]-yLocTop[i]))

            for j in xrange(yIdxTop, 0, -1):
                if (np.abs(V[j,i]-u0Top[i])/u0Top[i] > relDiff):
                    u0[i] = V[j-1,i]
                    yLoc[i] = Y[j-1,i]
                    break
    

    return [u0, yLoc]

def u0_from_grad(field, gradField, wall="bottom", relDiff=0.05, sanity=2000, start=0.02):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    grad = np.abs(gradField["V"])

    x = X[0,:]

    u0 = np.zeros(x.size) 
    yLoc = np.zeros(x.size) 

    if (wall == "top"):
        for i in xrange(x.size):
            startIdx = np.argmin(np.abs(Y[:,i]-(Y[-1,i]-start)))

            for yIdx in xrange(startIdx,0,-1):
                if (grad[yIdx,i] > sanity):
                    continue

                refDiffComp = np.abs(grad[yIdx,i]-grad[yIdx+1,i])/grad[yIdx+1,i]

                if (grad[yIdx,i] < sanity):
                    u0[i] = V[yIdx,i]
                    yLoc[i] = Y[yIdx,i]
                    break
    else:
        for i in xrange(x.size):
            startIdx = np.argmin(np.abs(Y[:,i]-(Y[0,i]+start)))

            for yIdx in xrange(startIdx,Y.shape[0]):
                if (grad[yIdx,i] > sanity):
                    continue

                refDiffComp = np.abs(grad[yIdx,i]-grad[yIdx-1,i])/grad[yIdx-1,i]

                if (grad[yIdx,i] < sanity):
                    u0[i] = V[yIdx,i]
                    yLoc[i] = Y[yIdx,i]
                    break

    return [u0, yLoc]

#Flat plate TBL estimates
def estimate_cf_from_ReX(ReX):
    return 0.0282*pow(ReX, -2.0/13)

def estimate_cf_from_ReDelta99(ReDelta99):
    return 0.0203*pow(ReDelta99, -2.0/11)

def estimate_ReDelta99_from_ReX(ReX):
    return 0.1635*pow(ReX, 11.0/13)

def estimate_ReTau_from_ReX(ReX):
    return 0.0194*pow(ReX, 10.0/13)

def estimate_ReTheta_from_ReX(ReX):
    return 0.0167*pow(ReX, 11.0/13)+378.43

def estimate_cf_from_ReTheta(ReTheta):
    return 0.0134*pow(ReTheta - 378.43, -2./11)



def field_as_meshgrid(field, nrows):
    """ LEGACY """
    X = field.x.as_matrix().reshape((nrows,-1),order='F')
    Y = field.y.as_matrix().reshape((nrows,-1),order='F')
    V = field.iloc[:,-1].as_matrix().reshape((nrows,-1),order='F')

    return [X, Y, V]

def profile_along_gridline(field, point, direction="y", index=-1):
    """Return the values along a path following a grid-line 
    Equivalent to a line-profile on a rectangular mesh
    
    Arguments:
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



def profile_values(x, field, nrows):
    """LEGACY, use profile_along_gridline()
    Return the values of a field along y for a given x
    
    Arguments:
        x -- the value of x
        field -- the field as a pandas DataFrame
        nrows -- the number of cells along y
    """

    X = field.x.as_matrix().reshape((nrows,-1),order='F')
    Y = field.y.as_matrix().reshape((nrows,-1),order='F')
    V = field.iloc[:,-1].as_matrix().reshape((nrows,-1),order='F')

    return profile_along_gridline({'X':X, 'Y':Y, 'V':V}, x)


def load_field(path, nrows):
    """LEGACY
       Load the field from file and do proper reordering
    
    Arguments:
        path -- string with path to the file
        nrows -- number of cells along y
    """

    field = pd.DataFrame.from_csv(path, sep=' ').sort(['x'])
    X = field.x.as_matrix().reshape((nrows,-1),order='F')
    Y = field.y.as_matrix().reshape((nrows,-1),order='F')
    V = field.iloc[:,-1].as_matrix().reshape((nrows,-1),order='F')
    
    idx = np.argsort(Y, axis=0)

    for i in range(Y.shape[1]):
        X[:,i] = X[:,i][idx[:,i]]
        Y[:,i] = Y[:,i][idx[:,i]]
        V[:,i] = V[:,i][idx[:,i]]
    
    field['x'] = X.flatten(order='F')
    field['y'] = Y.flatten(order='F')
    field[field.columns[-1]] = V.flatten(order='F')

    return field

def plot_field(field, cmap="RdBu_r", aspect="equal"):
    """ Plot the 2d field as a contour plot."""

    pcm = plt.pcolormesh(field["X"], field["Y"], field["V"], cmap='RdBu_r')
    pcm.axes.set_aspect(aspect)
    plt.xlim([np.min(field["X"][:,0]), np.max(field["X"][:,-1])])
    plt.ylim([np.min(field["Y"][0,:]), np.max(field["Y"][-1,:])])
    plt.colorbar(orientation="horizontal")

def open_case(path):
    case = SolutionDirectory(path, parallel=False, archive=None)    
    blockData = block_data(case)
    return {'case':case, 'blockData':blockData}

def block_data(case):

    # Open the blockMeshDict
    blockMeshDict = ParsedParameterFile(case.polyMeshDir()+"/blockMeshDict")

    # Read in the blocks data
    blockData = blockMeshDict["blocks"]

    # Rearrange the data into a 2d structure, 1 tuple per block
    blockData = zip(*[iter(blockData)]*5)

    # Number of blocks
    nBlocks = len(blockData)

    # Block connectivity matrix 
    connectivityMatrix = -1*np.ones((nBlocks, 4), dtype=np.int)


    # Matrix with the faces of the blocks
    faces = []
    for i in range(nBlocks):
        vertices = blockData[i][1]
        facesI = []
        facesI.append([vertices[0], vertices[3], vertices[7], vertices[4]])
        facesI.append([vertices[1], vertices[2], vertices[6], vertices[5]])
        facesI.append([vertices[7], vertices[6], vertices[2], vertices[3]])
        facesI.append([vertices[4], vertices[5], vertices[1], vertices[0]])
        faces.append(facesI)

    # Search through blocks to find shared faces and
    # establish connectivity
    for blockI in range(nBlocks):
        for blockJ in range(nBlocks):
            if (blockI != blockJ):
                for faceI in range(4):
                    if faces[blockI][faceI] in faces[blockJ]:
                        connectivityMatrix[blockI][faceI] = blockJ


    # Fill in the sizes of the blocks
    sizesMatrix = np.zeros((nBlocks, 3), dtype=np.int)

    for blockI in range(nBlocks):
        sizesMatrix[blockI][0] = blockData[blockI][2][0]
        sizesMatrix[blockI][1] = blockData[blockI][2][1]
        sizesMatrix[blockI][2] = blockData[blockI][2][0]* \
                                    blockData[blockI][2][1]
        
    return {'connectivity':connectivityMatrix, 'sizes':sizesMatrix}


def open_field(case, fieldName, time, utilName = "averageAlongAxis"):

    blockData = case["blockData"]
    case = case["case"]
    # Form the path and load in the field as a database
    fieldPath = path.join(case.systemDir(),"..","postProcessing", utilName, \
            str(time),fieldName)
    field = pd.DataFrame.from_csv(fieldPath, sep=' ')
    
    blockSizes = blockData["sizes"]
    blockConnectivity = blockData["connectivity"]
    nBlocks = blockSizes.shape[0]

    # Find the block with no neighbour on the left and bottom 
    startBlock = -1
    for blockI in range(nBlocks):
        if ((blockConnectivity[blockI][0] == -1) and \
            (blockConnectivity[blockI][3] == -1)):
            startBlock = blockI
    
    if (startBlock == -1):
        print "Error in open_field(): found no starting block"

    # Start index for each block
    startIndices = [0]
    for i in range(nBlocks-1):
        startIndices.append(startIndices[i]+blockSizes[i][2])

    # Load the data from the blocks into separate matrices
    blockFields = []
    for blockI in range(nBlocks):
        startIdx = startIndices[blockI]
        X = field[startIdx:startIdx+blockSizes[blockI][2]].as_matrix()[:,0] \
                .reshape(blockSizes[blockI][1],-1)
        Y = field[startIdx:startIdx+blockSizes[blockI][2]].as_matrix()[:,1] \
                .reshape(blockSizes[blockI][1],-1)
        V = field[startIdx:startIdx+blockSizes[blockI][2]].as_matrix()[:,-1] \
                .reshape(blockSizes[blockI][1],-1)
        blockFields.append((X,Y,V))

    # Block assembly process. This assumes that we can map the data
    # to a rectangular grid

    # Start with the block in the lower-left
    currentBlock = startBlock
    rowStartBlock = startBlock

    while (rowStartBlock != -1):
        currentBlock = rowStartBlock
        nextBlock = blockConnectivity[currentBlock,1]
        
        # Variables containing strips along "x"
        rowX = blockFields[currentBlock][0]
        rowY = blockFields[currentBlock][1]
        rowV = blockFields[currentBlock][2]
        while (nextBlock != -1):
            rowX = np.concatenate((rowX, blockFields[nextBlock][0]),axis=1)
            rowY = np.concatenate((rowY, blockFields[nextBlock][1]),axis=1)
            rowV = np.concatenate((rowV, blockFields[nextBlock][2]),axis=1)
            currentBlock = nextBlock
            nextBlock = blockConnectivity[nextBlock,1]

        # Concatenate in the y direction
        if( rowStartBlock == startBlock):
            X = rowX
            Y = rowY
            V = rowV
        else:
            X = np.concatenate((X, rowX), axis=0)
            Y = np.concatenate((Y, rowY), axis=0)
            V = np.concatenate((V, rowV), axis=0)

        rowStartBlock = blockConnectivity[rowStartBlock,2]

    startIdx = startIdx + blockSizes[-1][2]

    # Add boundary values at the top and bottom
    # Assumes there are only two patches
    for i in range(2):
        patchX = field[startIdx:startIdx+X[0,:].size].as_matrix()[:,0] \
                .reshape(1,-1)
        patchY = field[startIdx:startIdx+X[0,:].size].as_matrix()[:,1] \
                .reshape(1,-1)
        patchV = field[startIdx:startIdx+X[0,:].size].as_matrix()[:,-1] \
                .reshape(1,-1)
        if (patchY[0,0] < Y[0,0]):
            X = np.concatenate((patchX, X), axis = 0)
            Y = np.concatenate((patchY, Y), axis = 0)
            V = np.concatenate((patchV, V), axis = 0)
        else:
            X = np.concatenate((X, patchX), axis = 0)
            Y = np.concatenate((Y, patchY), axis = 0)
            V = np.concatenate((V, patchV), axis = 0)
        
        startIdx = startIdx + X[0,:].size

    return {'X':X, 'Y':Y, 'V':V}
