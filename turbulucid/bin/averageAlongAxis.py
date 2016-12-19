""" 
ASSMPTIONS
    * The averaging direction is the z axis
    * The geometry is "fit in" between to x-y planes.
"""

from __future__ import print_function
from __future__ import division
import os
import argparse
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import *
from collections import OrderedDict
from mpi4py import MPI
from sys import exit


def get_block_names(blocks):
    """ Return the names of the blocks in a given multiblock dataset

    Parameters
    ----------
        blocks : vtkMultiBlockDataSet
            The dataset with the blocks.

    Returns
    -------
        names : list of str
            A list with the names of the blocks.

    """
    names = []
    for i in range(blocks.GetNumberOfBlocks()):
        blockName = blocks.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
        names.append(blockName)

    return names


def get_block_index(blocks, name):
    """Get the index of the block by name.
    
    Parameters
    ----------
        blocks : vtkMultiBlockDataSet
            The dataset with the blocks.
        name : str
            The name of the block that is sought.
    
    Returns
    -------
        int
            The index of the sought block.

    """
    
    number = -1
    for i in range(blocks.GetNumberOfBlocks()):
        if (blocks.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME()) ==
            name):
            number = i
            break

    if number == -1:
        raise NameError("No block named "+name+" found")
    
    return number


def read(casePath):

    # Check that paths are valid
    if not os.path.exists(casePath):
        raise ValueError("Provided path to .foam file invalid!")

    # Case reader
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(casePath)
    reader.Update()

    reader.CreateCellToPointOff()
    reader.DisableAllPointArrays()
    reader.EnableAllPatchArrays()
    reader.Update()

    reader.SetTimeValue(vtk_to_numpy(reader.GetTimeValues())[-1])
    reader.Update()
    reader.UpdateInformation()

    return reader


def print_progress(current, total, freq=10., tabLevel=0):
    """Prints the number of percents complete with given frequency."""
    tabs = ""

    for i in range(tabLevel):
        tabs += "    "

    if np.mod(current, total/freq) == 0:
        print(tabs+"Done about "+str(int(current/total*100.))+"%")
              

def config_to_dict(configPath):
    """Parse a config file to a dictionary."""
    configFile = open(configPath, mode='r')

    configDict = {}

    for line in configFile:
        if (line[0] == '#') or (line == '\n'):
            continue
        configDict[line.split()[0]] = line.split()[1]

    configFile.close()
    return configDict


def minimal_length(patchData):
    """ Compute minimal length as sqrt of smallest area on a patch."""
    
    areaFilter = vtk.vtkMeshQuality()
    areaFilter.SetInputData(patchData) 
    areaFilter.SetTriangleQualityMeasureToArea()
    areaFilter.SetQuadQualityMeasureToArea()
    areaFilter.Update()
    area = dsa.WrapDataObject(areaFilter.GetOutput())
    area = area.CellData["Quality"]

    return np.sqrt(np.min(area))


def average_internal_field_data(block, patchData, nSamples):

    blockData = block.GetCellData()
    bounds = block.GetBounds()
    smallDz = (bounds[5] - bounds[2])/10000

    patchCellCenters = vtk.vtkCellCenters()
    patchCellCenters.SetInputData(patchData)
    patchCellCenters.Update()

    patchCellCenters = patchCellCenters.GetOutput()
    nSeedPoints = patchCellCenters.GetNumberOfPoints()

    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(block)
    patchCellData = patchData.GetCellData()

    averageFields = OrderedDict()

    nFields = blockData.GetNumberOfArrays()

    for field in range(nFields):
        name = blockData.GetArrayName(field)
        nCols = blockData.GetArray(field).GetNumberOfComponents()
        averageFields[name] = np.zeros( (nSeedPoints, nCols))

    print("Sampling and averaging internal field")
    for seed in range(int(nSeedPoints)):
        print_progress(seed, nSeedPoints)

        seedPoint = patchCellCenters.GetPoint(seed)
        line.SetResolution(nSamples)
        line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
        line.SetPoint2(seedPoint[0], seedPoint[1], bounds[5]-smallDz)
        line.Update()

        probeFilter.SetInputConnection(line.GetOutputPort())
        probeFilter.Update()

        probeData = dsa.WrapDataObject(probeFilter.GetOutput())

        for field in averageFields:
            averageFields[field][seed] = np.mean(probeData.PointData[field],
                                                 axis=0)

    wrappedPatchData = dsa.WrapDataObject(patchData)
    for field in averageFields:
        fieldI = averageFields[field]
        nComp = fieldI.shape[1]

        if nComp == 1:  # scalar
            wrappedPatchData.CellData[field][:] = fieldI[:, 0]
        else:
            wrappedPatchData.CellData[field][:, :] = fieldI[:, :]


def average_patch_data(data, patchPolys, nSamples, bounds):
    print("Averaging data for patches")

    patchAveragedFields = OrderedDict()

    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()

    smallDz = (bounds[5] - bounds[4])/1000.

    for boundary in patchPolys:
        print("Patch "+boundary)

        polyI = patchPolys[boundary]
        zero_out_arrays(polyI)

        blockNumber = get_block_index(data.GetBlock(1), boundary)
        patchBlock = data.GetBlock(1).GetBlock(blockNumber)
        patchBlockData = patchBlock.GetCellData()
        nSeedPoints = polyI.GetNumberOfPoints()

        nFields = patchBlockData.GetNumberOfArrays()

        for field in range(nFields):
            name = patchBlockData.GetArrayName(field)
            nCols = patchBlockData.GetArray(field).GetNumberOfComponents()
            patchAveragedFields[name] = np.zeros((nSeedPoints, nCols))

        probeFilter.SetSourceData(patchBlock)

        for seed in range(nSeedPoints):

            seedPoint = polyI.GetPoint(seed)
            line.SetResolution(nSamples)
            line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
            line.SetPoint2((seedPoint[0], seedPoint[1], bounds[5]-smallDz))
            line.Update()

            probeFilter.SetInputConnection(line.GetOutputPort())
            probeFilter.Update()

            probeData = dsa.WrapDataObject(probeFilter.GetOutput())

            for field in patchAveragedFields:
                patchAveragedFields[field][seed] = \
                    np.mean(probeData.PointData[field], axis=0)

        wrappedPoly = dsa.WrapDataObject(polyI)

        for field in patchAveragedFields:
            fieldI = patchAveragedFields[field]
            nComp = fieldI.shape[1]
            if nComp == 1:  # scalar
                wrappedPoly.PointData[field][:] = fieldI[:, 0]
            elif nComp == 3:  # vector
                wrappedPoly.PointData[field][:, :] = fieldI[:, :]


def colour_boundary_points(patchData, boundaryPoints):
    pointColouring = np.zeros((patchData.GetNumberOfPoints(),
                               len(boundaryPoints)),
                              dtype=np.int64)
    pointColouring[:] = 0
    for boundary in boundaryPoints:
        for pointI in range(boundaryPoints[boundary].shape[0]):
            pointId = patchData.FindPoint(
                boundaryPoints[boundary][pointI, :])
            pointColouring[pointId, boundaryPoints.keys().index(boundary)] = 1

    return pointColouring


def mark_boundary_cells(patchData, patchPolys, boundaryPoints):
    print("Marking boundary cells.")

    # For each point on a given boundary find the id of this point in the
    # patch data. Then set the colouring to 1 in at row number with found
    # id and column number according to boundary ordering
    pointColouring = colour_boundary_points(patchData, boundaryPoints)

    boundaryCellsConn = {}

    for field in range(len(patchPolys)):
        key = list(patchPolys.keys())[field]
        boundaryCellsConn[key] = -1*np.ones(patchPolys[key].GetNumberOfPoints())

    for cellI in range(patchData.GetNumberOfCells()):
        print_progress(cellI, patchData.GetNumberOfCells(), 5, 1)
        cellPointsIdsVtk = vtk.vtkIdList()
        patchData.GetCellPoints(cellI, cellPointsIdsVtk)
        cellPointsIds = np.zeros(cellPointsIdsVtk.GetNumberOfIds(),
                                 dtype=np.int64)
        for i in range(cellPointsIds.size):
            cellPointsIds[i] = cellPointsIdsVtk.GetId(i)

        for boundary in boundaryPoints:
            colIdx = boundaryPoints.keys().index(boundary)
            boundaryPointIds = np.where(pointColouring[:,colIdx] == 1)[0]
            intersection = np.intersect1d(cellPointsIds, boundaryPointIds)
            if intersection.size >= 2:
                point1 = np.array(patchData.GetPoint(intersection[0]))
                point2 = np.array(patchData.GetPoint(intersection[1]))
                midPoint = 0.5*(point1 + point2)
                idOnPatchPoly = patchPolys[boundary].FindPoint(midPoint)
                boundaryCellsConn[boundary][idOnPatchPoly] = cellI

    # Check if some connectivity is not establishes
    for key in boundaryCellsConn:
        if np.any(boundaryCellsConn[key] == -1):
            print("ERROR: some connectivity not established for boundary "+key)

    for key in boundaryCellsConn:
        toVtk = numpy_to_vtk(boundaryCellsConn[key])
        toVtk.SetName(key)
        patchData.GetFieldData().AddArray(toVtk)


def assemble_multiblock(patchData, patchPolys):
    print("Assembling multi-block structure")
    multiBlock = vtk.vtkMultiBlockDataSet()
    multiBlock.SetNumberOfBlocks(len(patchPolys) + 1)
    multiBlock.SetBlock(0, patchData)
    multiBlock.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(),
                                  "internalField")

    i = 1
    for boundary in patchPolys:
        multiBlock.SetBlock(i, patchPolys[boundary])
        multiBlock.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), boundary)
        i += 1

    return multiBlock


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script for averaging a 3D field along z.")

    parser.add_argument('--config',
                        type=str,
                        help='The config file.',
                        required=True)

    return parser


def zero_out_arrays(polyData):
    wrapped = dsa.WrapDataObject(polyData)

    for field in wrapped.PointData.keys():
        wrapped.PointData[field][:] = np.zeros(wrapped.PointData[field].shape)

    for field in wrapped.CellData.keys():
        wrapped.CellData[field][:] = np.zeros(wrapped.CellData[field].shape)





def main():

    parser = create_parser()
    args = parser.parse_args()

    config = config_to_dict(args.config)
    try:
        casePath = config["case"]
        seedPatchName = config["patch"]
        time = float(config["time"])
        nSamples = int(config["nSamples"])
        writePath = config["file"]
    except KeyError:
        print("ERROR: required parameter not specified in config file.")
        raise

    # Case reader
    print("Reading")
    reader = read(casePath)

    # Writer
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(writePath)

    caseData = reader.GetOutput()
    internalBlock = caseData.GetBlock(0)
    bounds = internalBlock.GetBounds()

    patchBlocks = caseData.GetBlock(1)
    patchBlockI = patchBlocks.GetBlock(get_block_index(patchBlocks,
                                                       seedPatchName))

    patchData = vtk.vtkPolyData()
    patchData.ShallowCopy(patchBlockI)

    zero_out_arrays(patchData)

    #average_internal_field_data(internalBlock, patchData, nSamples)

    boundaryNames = vtk.vtkStringArray()
    boundaryNames.SetName("boundaries")

    # Find seed points for the patches
    print("Finding seed points for the patches")

    patchFeatureEdgesFilter = vtk.vtkFeatureEdges()

    # Polydata with points containing boundary data
    # also used as seeds for sampling the patch data
    patchPolys = OrderedDict()

    # The boundary points extracted from the patches.
    # Needed for finding boundary cells and associating them with the boundary
    # names
    boundaryPoints = OrderedDict()

    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
    selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)

    #for field in range(nBoundaries-1):
    for patchName in get_block_names(patchBlocks):
        # Get the patch data
        patchBlockI = patchBlocks.GetBlock(get_block_index(patchBlocks,
                                                           patchName))

        # Extract feature edges
        patchFeatureEdgesFilter.SetInputData(patchBlockI)
        patchFeatureEdgesFilter.Update()
        patchFeatureEdgesData = patchFeatureEdgesFilter.GetOutput()

        boundaryIPoints = dsa.WrapDataObject(patchFeatureEdgesData).Points

        # Extract cell centers on the feature edges
        patchFeatureEdgesCCsFilter = vtk.vtkCellCenters()
        patchFeatureEdgesCCsFilter.SetInputData(patchFeatureEdgesData)
        patchFeatureEdgesCCsFilter.Update()

        patchFeatureEdgesCCs = dsa.WrapDataObject(
                patchFeatureEdgesCCsFilter.GetOutput())

        ccPoints = patchFeatureEdgesCCs.Points

        # The patch is an x-y plane
        if (np.all(boundaryIPoints[:, 2] == bounds[4]) or
            np.all(boundaryIPoints[:, 2] == bounds[5])):
                continue
        else:
            boundaryNames.InsertNextValue(patchName)

            # Select the points located at the boundary of the patch
            idx = ccPoints[:, 2] == patchData.GetPoint(0)[2]
            ids = np.where(idx == True)

            idx = boundaryIPoints[:, 2] == patchData.GetPoint(0)[2]
            ids = np.where(idx == True)

            cellIds = vtk.vtkIntArray()
            for i in range(patchFeatureEdgesData.GetNumberOfCells()):
                cellI = patchFeatureEdgesData.GetCell(i)
                pointIds = vtk.vtkIdList()
                patchFeatureEdgesData.GetCellPoints(i, pointIds)

                flag = 0
                for j in range(pointIds.GetNumberOfIds()):
                    if not idx[pointIds.GetId(j)]:
                        break
                    flag += 1

                if flag == 2:
                    cellIds.InsertNextValue(i)


            newPoly = vtk.vtkPolyData()
            newPoly.ShallowCopy(patchFeatureEdgesFilter.GetOutput())
            cellIds2 = vtk_to_numpy(cellIds)
            for i in range(newPoly.GetNumberOfCells()):
                if i not in cellIds2:
                    newPoly.DeleteCell(i)

            newPoly.RemoveDeletedCells()

            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(newPoly)
            cleaner.Update()

            newPoly.ShallowCopy(cleaner.GetOutput())



            selectionNode.SetSelectionList(numpy_to_vtk(ids[0]))

            #selectionNode.SetSelectionList(cellIds)

            selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
            #selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
            selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
            selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1)

            selection = vtk.vtkSelection()
            selection.AddNode(selectionNode)

            extractSelection = vtk.vtkExtractSelection()
            #extractSelection.SetInputConnection(0,
            #    patchFeatureEdgesCCsFilter.GetOutputPort())
            extractSelection.SetInputConnection(0,
                                                patchFeatureEdgesFilter.GetOutputPort())
            extractSelection.SetInputData(1, selection)
            extractSelection.Update()

            #newPoly = vtk.vtkPolyData()
            #newPoly.ShallowCopy(extractSelection.GetOutput())

            #newPoly.ShallowCopy(patchFeatureEdgesFilter.GetOutput())

            patchPolys[patchName] = newPoly

            #idx = boundaryIPoints[:, 2] == patchData.GetPoint(0)[2]
            boundaryPoints[patchName] = boundaryIPoints[idx, :]

    # Add the names of the boundaries as field data
    patchData.GetFieldData().AddArray(boundaryNames)

  #  average_patch_data(caseData, patchPolys, nSamples, bounds)

    minLength = minimal_length(patchData)
    tol = 0.001*minLength

   # mark_boundary_cells(patchData, patchPolys, boundaryPoints)

    multiBlock = assemble_multiblock(patchData, patchPolys)

    print("Writing")
    writer.SetInputData(multiBlock)
    writer.Write()


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    nProcs = comm.size

    print(nProcs)
    main()
