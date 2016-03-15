from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from os import path
import numpy as np
import pandas as pd

__all__ = ["open_case", "open_field"]

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

