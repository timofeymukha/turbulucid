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
import matplotlib.pyplot as plt
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import *
from collections import OrderedDict
from sys import exit


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


def read_case(casePath):
    print(casePath)
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


def average_patch_data(data, patchPolys, nSamples, bounds):
    print("Averaging data for patches")

    patchAveragedFields = OrderedDict()

    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()

    smallDz = (bounds[5] - bounds[4])/1000.

    for boundary in patchPolys:
        print("Patch "+boundary)
        blockNumber = get_block_index(data.GetBlock(1), boundary)
        patchBlock = data.GetBlock(1).GetBlock(blockNumber)
        patchBlockData = patchBlock.GetCellData()
        nSeedPoints = patchPolys[boundary].GetNumberOfPoints()

        nFields = patchBlockData.GetNumberOfArrays()

        for field in range(nFields):
            name = patchBlockData.GetArrayName(field)
            nCols = patchBlockData.GetArray(field).GetNumberOfComponents()
            patchAveragedFields[name] = np.zeros((nSeedPoints, nCols))

        probeFilter.SetSourceData(patchBlock)

        for seed in range(nSeedPoints):
            #print_progress(seed, nSeedPoints, 5, 1)

            seedPoint = patchPolys[boundary].GetPoint(seed)
            line.SetResolution(nSamples)
            line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
            line.SetPoint2((seedPoint[0], seedPoint[1] ,bounds[5]-smallDz))
            line.Update()

            probeFilter.SetInputConnection(line.GetOutputPort())
            probeFilter.Update()

            probeData = dsa.WrapDataObject(probeFilter.GetOutput())

            for field in patchAveragedFields:
                patchAveragedFields[field][seed] = \
                    np.mean(probeData.PointData[field])

        for field in patchAveragedFields:
            nComp = patchAveragedFields[field].shape[1]
            if nComp == 1: #scalar
                patchBlockData.SetActiveScalars(field)
                patchBlockData.SetScalars(numpy_to_vtk(patchAveragedFields[field]))
            elif nComp == 3: #vector
                patchBlockData.SetActiveVectors(field)
                patchBlockData.SetVectors(numpy_to_vtk(patchAveragedFields[field]))
            elif nComp == 6: #symmetric tensor
                # Add three dummy components

                sixComp = patchAveragedFields[field]
                nineComp = np.column_stack((sixComp, np.zeros((sixComp.shape[0], 3))))
                #nineCompVTK = numpy_to_vtk(nineComp)
                #nineCompVTK.SetName(field)

                #patchBlockData.SetActiveTensors(field)
                #patchBlockData.SetTensors(nineCompVTK)
            elif nComp == 9: #tensor
                patchBlockData.SetActiveTensors(field)
                patchBlockData.SetTensors(numpy_to_vtk(patchAveragedFields[field]))


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


def main():
    parser = argparse.ArgumentParser(
            description="Script for averaging a 3D field along a chosen \
                            direction.")

    parser.add_argument('--config',
                        type=str,
                        help='The config file.',
                        required=False)

    parser.add_argument('-n',
                        type=int,
                        help='Number of samples.',
                        required=False)

    parser.add_argument('-c', '--case',
                        type=str,
                        help='Path to the OpenFOAM case. Default is cwd.',
                        default='.',
                        required=False)

    parser.add_argument('-p', '--patch',
                        type=str,
                        help='Path to the seed patch, saved as a vtk file.',
                        default=os.path.join(".", "postProcessing", "sufraces",
                                             "0", "seedPatch.vtk"),
                        required=False)

    parser.add_argument('-t', '--time',
                        type=float,
                        help='The time value for the fields to be averaged. \
                              Default is latest time.',
                        required=False)

    parser.add_argument('-f', '--file',
                        type=str,
                        help='The name of the output file.',
                        required=False)

    args = parser.parse_args()

    if args.config is not None:
        config = config_to_dict(args.config)
        try:
            path = config["case"]
            pathPatch = config["patch"]
            time = float(config["time"])
            nSamples = int(config["nSamples"])
            writePath = config["file"]
        except KeyError:
            print("ERROR: required parameter not specified in config file.")
            raise
    else:
        path = args.case
        pathPatch = args.patch
        nSamples = args.n
        time = args.time
        writePath = args.file

    print("Reading")
    # Case reader
    reader = read_case(path)

    # Seed patch reader
    patchReader = vtk.vtkPolyDataReader()
    patchReader.SetFileName(pathPatch)
    patchReader.Update()

    # Writer
    #writer = vtk.vtkPolyDataWriter()
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(writePath)

    data = reader.GetOutput()
    block = data.GetBlock(0)
    blockData = block.GetCellData()
    bounds = block.GetBounds()

    smallDz = (bounds[5] - bounds[2])/10000

    patchData = patchReader.GetOutput()
    nCells = patchData.GetNumberOfCells()

    patchCellCenters = vtk.vtkCellCenters()
    patchCellCenters.SetInputData(patchData)
    patchCellCenters.Update()

    patchCellCenters = patchCellCenters.GetOutput()
    nSeedPoints = patchCellCenters.GetNumberOfPoints()

    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(block)
    patchCellData = patchData.GetCellData()

    averageFields = []
    fieldNames = []
    nFields = blockData.GetNumberOfArrays()

    for field in range(nFields):
        name = blockData.GetArrayName(field)
        newAverageField = vtk.vtkFloatArray()
        newAverageField.SetName(name)
        newAverageField.SetNumberOfComponents(blockData.GetArray(field).GetNumberOfComponents())
        newAverageField.SetNumberOfTuples(nSeedPoints)
        averageFields.append(newAverageField)
        fieldNames.append(name)

    print("Sampling and averaging internal field")
    for seed in range(int(nSeedPoints)):
        print_progress(seed, nSeedPoints)

        seedPoint = patchCellCenters.GetPoint(seed)
        line.SetResolution(nSamples)
        line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
        line.SetPoint2((seedPoint[0], seedPoint[1] ,bounds[5]-smallDz))
        line.Update()

        probeFilter.SetInputConnection(line.GetOutputPort())
        probeFilter.Update()


        probeData = probeFilter.GetOutput().GetPointData()

        for fieldI in range(nFields):
            field = probeData.GetArray(fieldI)
            name = probeData.GetArray(fieldI).GetName()
            field = vtk_to_numpy(field)
            averageField = np.array(np.mean(field, axis=0), dtype=np.float32)
            idx = fieldNames.index(name)

            for comp in range(averageField.size):
                if averageField.size > 1:
                    averageFields[idx].SetComponent(seed, comp, averageField[comp])
                else:
                    averageFields[idx].SetValue(seed, averageField)

    # Add the averaged data to the patch
    for field in range(nFields):
        nComp = averageFields[field].GetNumberOfComponents()
        if nComp == 1:  # scalar
            patchCellData.SetActiveScalars(averageFields[field].GetName())
            patchCellData.SetScalars(averageFields[field])
        elif nComp == 3:  # vector
            patchCellData.SetActiveVectors(averageFields[field].GetName())
            patchCellData.SetVectors(averageFields[field])
        elif nComp == 6:  # symmetric tensor
            # Add three dummy components
            nineComp = vtk_to_numpy(averageFields[field])
            nineComp = np.column_stack((nineComp, np.zeros((nineComp.shape[0], 3))))
            nineCompVTK = numpy_to_vtk(nineComp)
            nineCompVTK.SetName(averageFields[field].GetName())

            patchCellData.SetActiveTensors(averageFields[field].GetName())
            patchCellData.SetTensors(nineCompVTK)
        elif nComp == 9:  # tensor
            patchCellData.SetActiveTensors(averageFields[field].GetName())
            patchCellData.SetTensors(averageFields[field])

    nBoundaries = reader.GetNumberOfPatchArrays()
    boundaryNames = vtk.vtkStringArray()
    # one is actually the internal field and 2 are x-y planes
    boundaryNames.SetNumberOfValues(nBoundaries - 3)
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

    for field in range(nBoundaries-1):
        # Get the patch data
        patchBlock = data.GetBlock(1).GetBlock(field)

        # Extract feature edges
        patchFeatureEdgesFilter.SetInputData(patchBlock)
        patchFeatureEdgesFilter.Update()
        patchFeatureEdgesData = patchFeatureEdgesFilter.GetOutput()

        boundaryIPoints = dsa.WrapDataObject(patchFeatureEdgesData).Points

        # Extract cell centers on the feature edges
        patchFeatureEdgesCCsFilter =  vtk.vtkCellCenters()
        patchFeatureEdgesCCsFilter.SetInputData(patchFeatureEdgesData)
        patchFeatureEdgesCCsFilter.Update()

        patchFeatureEdgesCCs = dsa.WrapDataObject(
                patchFeatureEdgesCCsFilter.GetOutput())
        ccPoints = patchFeatureEdgesCCs.Points

        # The patch is an x-y plane
        if (np.all(ccPoints[:,2] == bounds[4]) or
            np.all(ccPoints[:,2] == bounds[5])):
                continue
        else:
            nameNum = 0
            for j in range(boundaryNames.GetNumberOfValues()):
                if boundaryNames.GetValue(j) == '':
                    boundaryNames.SetValue(j, reader.GetPatchArrayName(field + 1))
                    nameNum = j
                    break

            idx = ccPoints[:, 2] == bounds[4]
            ids =  np.where(idx == True)
            selectionNode.SetSelectionList(numpy_to_vtk(ids[0]))

            selection = vtk.vtkSelection()
            selection.AddNode(selectionNode)

            extractSelection = vtk.vtkExtractSelection()
            extractSelection.SetInputConnection(0,
                patchFeatureEdgesCCsFilter.GetOutputPort())
            extractSelection.SetInputData(1, selection)
            extractSelection.Update()

            newPoly = vtk.vtkPolyData()
            newPoly.ShallowCopy(extractSelection.GetOutput())

            patchPolys[boundaryNames.GetValue(nameNum)] = newPoly

            idx = boundaryIPoints[:, 2] == patchData.GetPoint(0)[2]
            boundaryPoints[boundaryNames.GetValue(nameNum)] = boundaryIPoints[idx, :]

    # Add the names of the boundaries as field data
    patchData.GetFieldData().AddArray(boundaryNames)

    average_patch_data(data, patchPolys, nSamples, bounds)

    minLength = minimal_length(patchData)
    tol = 0.001*minLength

    # For each point on a given boundary find the id of this point in the
    # patch data. Then set the colouring to 1 in at row number with found
    # id and column number according to boundary ordering

    print("Marking boundary cells.")
    pointColouring =  colour_boundary_points(patchData, boundaryPoints)

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

    multiBlock = assemble_multiblock(patchData, patchPolys)

    print("Writing")
    writer.SetInputData(multiBlock)
    writer.Write()


if __name__ == '__main__':
    main()
