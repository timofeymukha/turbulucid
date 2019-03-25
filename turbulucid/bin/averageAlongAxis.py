# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

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
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import *
from collections import OrderedDict
from sys import exit


def get_block_names(blocks):
    """ Return the names of the blocks in a given multiblock dataset

    Parameters
    ----------
    blocks : vtkMultiBlockDataSet
        The dataset with the blocks.

    Returns
    -------
    list of str
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


def print_progress(current, total, freq=10., tabLevel=0):
    """Prints the number of percents complete with given
    frequency.

    """
    tabs = ""

    for i in range(tabLevel):
        tabs += "    "

    if np.mod(current, int(total/freq)) == 0:
        print(tabs+"Done about "+str(int(current/total*100.))+"%")


def read(casePath, time, debug=False):
    """Read the case from a given path to .foam file.

    Parameters
    ----------
    casePath : str
        The path to the .foam file.
    time : float
        The time step to load, default to latest time
    debug : bool
        Debug switch

    Returns
    -------
        The reader updated with the read case.

    Raises
    ------
    ValueError
        If the path is not valid.


    """
    # Check that paths are valid
    if not os.path.exists(casePath):
        raise ValueError("Provided path to .foam file invalid!")

    if debug:
        print("    Opening the case")
    # Case reader
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(casePath)
    reader.Update()

    if debug:
        print("    Changing reader parameters")
    reader.CreateCellToPointOff()
    reader.DisableAllPointArrays()
    reader.EnableAllPatchArrays()
    reader.DecomposePolyhedraOff()
    reader.Update()
    reader.UpdateInformation()

    info = reader.GetExecutive().GetOutputInformation(0)

    if debug:
        print("The available timesteps are", vtk_to_numpy(reader.GetTimeValues()))

    print(time)
    if time is None:
        print("Selecting the latest available time step")
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),
                vtk_to_numpy(reader.GetTimeValues())[-1])
    else:
        print("Selecting the time step", time)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), time)

    reader.Update()
    reader.UpdateInformation()

    return reader


def config_to_dict(configPath):
    """Parse a config file to a dictionary.

    """
    configFile = open(configPath, mode='r')

    configDict = {}

    for line in configFile:
        if (line[0] == '#') or (not line.strip()):
            continue
        configDict[line.split(None, 1)[0]] = line.split(None, 1)[1].rstrip()

    configFile.close()
    return configDict


def minimal_length(patchData):
    """ Compute minimal length as sqrt of smallest area on a patch.

    """
    areaFilter = vtk.vtkMeshQuality()
    areaFilter.SetInputData(patchData) 
    areaFilter.SetTriangleQualityMeasureToArea()
    areaFilter.SetQuadQualityMeasureToArea()
    areaFilter.Update()
    area = dsa.WrapDataObject(areaFilter.GetOutput())
    area = area.CellData["Quality"]

    return np.sqrt(np.min(area))


def get_cell_points(polyData, cellId):
    """Get the ids of the points of the cell with id cellId.

    Wraps the corresponding VTK function data.GetCellPoints().

    Parameters
    ----------
    polyData : vtkPolyData
        The data object with the cells.
    cellId : int
        The id of the cell.

    Returns
    -------
    list
        List of point ids.

    """
    cellPointsIdsVtk = vtk.vtkIdList()
    polyData.GetCellPoints(cellId, cellPointsIdsVtk)
    cellPointsIds = np.zeros(cellPointsIdsVtk.GetNumberOfIds(),
                             dtype=np.int64)
    for id in range(cellPointsIds.size):
        cellPointsIds[id] = cellPointsIdsVtk.GetId(id)

    return cellPointsIds


def average_internal_field_data(block, internalData, nSamples, debug, dry):

    blockCellData = block.GetCellData()
    bounds = block.GetBounds()
    smallDz = (bounds[5] - bounds[2])/10000

    if debug:
        print("    Computing cell centers of the seed patch")
    patchCellCenters = vtk.vtkCellCenters()
    patchCellCenters.SetInputData(internalData)
    patchCellCenters.Update()

    patchCellCenters = patchCellCenters.GetOutput()
    nSeedPoints = patchCellCenters.GetNumberOfPoints()

    if debug:
        print("    The number of seed points is", nSeedPoints)


    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(block)

    avrgFields = OrderedDict()

    nFields = blockCellData.GetNumberOfArrays()

    for field in range(nFields):
        name = blockCellData.GetArrayName(field)
        nCols = blockCellData.GetArray(field).GetNumberOfComponents()
        if debug:
            print("    Will average field", name, "with", nCols, "components")
        avrgFields[name] = np.zeros((nSeedPoints, nCols))

    if not dry:
        for seed in range(int(nSeedPoints)):
            print_progress(seed, nSeedPoints, tabLevel=1)

            seedPoint = patchCellCenters.GetPoint(seed)
            line.SetResolution(nSamples-1)
            line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
            line.SetPoint2(seedPoint[0], seedPoint[1], bounds[5]-smallDz)
            line.Update()

            probeFilter.SetInputConnection(line.GetOutputPort())
            probeFilter.Update()

            probeData = dsa.WrapDataObject(probeFilter.GetOutput()).PointData

            for field in avrgFields:
                if avrgFields[field].shape[1] == 9:  # a tensor
                    reshaped = probeData[field].reshape((nSamples, 9))
                    avrgFields[field][seed] = np.mean(reshaped, axis=0)
                else:
                    avrgFields[field][seed] = np.mean(probeData[field], axis=0)
    else:
        print("    This is a dry run, will not actually average")

    if debug:
        print("    Assigning sampled data to the internal field")
    wrappedPatchData = dsa.WrapDataObject(internalData)
    for field in avrgFields:
        print("        Assigning field", field)
        fieldI = avrgFields[field]
        nComp = fieldI.shape[1]

        if nComp == 1:  # scalar
            wrappedPatchData.CellData[field][:] = fieldI[:, 0]
        elif nComp == 9:  # tensor
            wrappedPatchData.CellData[field][:] = fieldI.reshape(nSeedPoints,
                                                                 3, 3)
        else:
            wrappedPatchData.CellData[field][:, :] = fieldI[:, :]


def average_patch_data(data, boundaryData, nSamples, bounds):

    avergFields = OrderedDict()

    line = vtk.vtkLineSource()
    probeFilter = vtk.vtkProbeFilter()

    cellCenters = vtk.vtkCellCenters()

    smallDz = (bounds[5] - bounds[4])/1000.

    for boundary in boundaryData:
        print("Patch "+boundary)

        polyI = boundaryData[boundary]
        zero_out_arrays(polyI)

        blockNumber = get_block_index(data.GetBlock(1), boundary)
        patchBlock = data.GetBlock(1).GetBlock(blockNumber)
        patchBlockData = patchBlock.GetCellData()
        nSeedPoints = polyI.GetNumberOfCells()

        cellCenters.SetInputData(polyI)
        cellCenters.Update()

        nFields = patchBlockData.GetNumberOfArrays()

        for field in range(nFields):
            name = patchBlockData.GetArrayName(field)
            nCols = patchBlockData.GetArray(field).GetNumberOfComponents()
            avergFields[name] = np.zeros((nSeedPoints, nCols))

        probeFilter.SetSourceData(patchBlock)

        for seed in range(nSeedPoints):

            seedPoint = cellCenters.GetOutput().GetPoint(seed)
            line.SetResolution(nSamples - 1)
            line.SetPoint1(seedPoint[0], seedPoint[1], bounds[4]+smallDz)
            line.SetPoint2((seedPoint[0], seedPoint[1], bounds[5]-smallDz))
            line.Update()

            probeFilter.SetInputConnection(line.GetOutputPort())
            probeFilter.Update()

            probeData = dsa.WrapDataObject(probeFilter.GetOutput()).PointData

            for field in avergFields:
                if avergFields[field].shape[1] == 9:  #tensor
                    reshaped = probeData[field].reshape((nSamples, 9))
                    avergFields[field][seed] = np.mean(reshaped, axis=0)
                else:
                    avergFields[field][seed] = np.mean(probeData[field],
                                                       axis=0)

        wrappedPoly = dsa.WrapDataObject(polyI)

        for field in avergFields:
            fieldI = avergFields[field]
            nComp = fieldI.shape[1]
            if nComp == 1:  # scalar
                wrappedPoly.CellData[field][:] = fieldI[:, 0]
            elif nComp == 9:  # tensor
                wrappedPoly.CellData[field][:, :] = fieldI.reshape(nSeedPoints,
                                                                   3, 3)
            else:
                wrappedPoly.CellData[field][:, :] = fieldI[:, :]


def get_point_ids(polyData, points):
    """Given a list of point coordinates return their ids in given data.

    Parameters
    ----------
    polyData : vtkPolyData
        The polydata for which the

    Returns
    -------
    ndarray
        Array with the ids of the points.

    """
    ids = np.zeros(points.shape, dtype=np.int64)

    for pointI in range(points.shape[0]):
        ids[pointI] = polyData.FindPoint(points[pointI, :])

    return ids


def get_closest_cell(point, internalData, debug):
    """For a given point, find the cell located closest to it.

    Based on vtkCell.EvaluatePosition.

    Parameters
    ---------
    point : triple
        The point for which to find the closest cell.
    internalData : polydata
        The polydata with the cells.
    debug : bool
        Debug switch

    Returns
    -------
        The id of the cell and the distance to it


    """
    distance = np.zeros(internalData.GetNumberOfCells())

    closestPoint = [0, 0, 0]
    subId = vtk.mutable(0)
    dist2 = vtk.mutable(0.0)
    pcoords = [0, 0, 0]
    weights = [0, 0, 0, 0]


    for i in range(distance.shape[0]):
        cellI = internalData.GetCell(i)
        found = cellI.EvaluatePosition(point, closestPoint, subId,
                                       pcoords, dist2, weights)
        distance[i] = dist2
        if found == -1:
            print("    WARNING: could not evaluate position for "
                  "cell", i)

    foundCellId = np.argmin(distance)

    return foundCellId, distance[foundCellId]


def mark_boundary_cells(patchData, patchPolys, debug):
    """Find the internal cell adjacent to each cell in the boundary
    data.

    Goes through each cell center for all boundary polydata and
    attempts to find the adjacent cell in the interalField.
    Creates a connectivity array and adds it as field data to the
    internal field polydata.

    """
    boundaryCellsConn = OrderedDict()
    cellCenters = vtk.vtkCellCenters()

    for boundary in patchPolys:
        boundaryCellsConn[boundary] = -1*np.ones(patchPolys[boundary].GetNumberOfCells(), dtype=np.int32)

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(patchData)
    locator.Update()

    for boundary in patchPolys:
        if debug:
            print("    Marking cells for patch", boundary)

        if debug:
            print("    Computing edge centers")
        polyI = patchPolys[boundary]
        cellCenters.SetInputData(polyI)
        cellCenters.Update()
        points = dsa.WrapDataObject(cellCenters.GetOutput()).Points
        if debug:
            print("    Found", points.shape[0], "centers")

        if debug:
            print("    Computing normals")

        normals = np.zeros(points.shape)
        for i in range(polyI.GetNumberOfCells()):
            cellI = polyI.GetCell(i)
            tan = np.array(cellI.GetPoints().GetPoint(1)) - np.array(cellI.GetPoints().GetPoint(0)[:])
            normals[i] = np.cross(tan, [0, 0, 1])
            normals[i] /= np.linalg.norm(normals[i])

        if debug:
            print("    The mean normal is ", np.mean(normals, axis=0))
        cell = vtk.vtkGenericCell()
        tol2 = 0.0
        pcoords = [0, 0, 0]
        weights = [0, 0, 0]

        for i in range(points.shape[0]):
            pointI = points[i, :]
            foundCellId = locator.FindCell(pointI, tol2, cell, pcoords, weights)

            # Attmept going along the normal
            if foundCellId == -1:
                print("    Failed to find adjacent cell for boundary point",
                      pointI, "on boundary", boundary)
                print("    Attempting to perturb location a long the normal")
                pointI += 1e-6*normals[i]
                foundCellId = locator.FindCell(pointI, tol2, cell, pcoords, weights)
            if foundCellId == -1:
                pointI -= 1e-6*normals[i]
                foundCellId = locator.FindCell(pointI, tol2, cell, pcoords, weights)


            boundaryCellsConn[boundary][i] = foundCellId

            if foundCellId == -1:
                print("    Attempting with slow algorithm based on minimum"
                      " distance")
                foundCellId, distance = get_closest_cell(pointI, patchData, debug)

                print("    Found cell with id", foundCellId, "located",
                      distance, "away.")
                boundaryCellsConn[boundary][i] = foundCellId

    if debug:
        print("    Assigning the connectivity lists as FieldData")
    for key in boundaryCellsConn:
        if np.any(boundaryCellsConn[key] == -1):
            print("ERROR: some connectivity not established for boundary "+key)

        wrappedData = dsa.WrapDataObject(patchData)
        wrappedData.FieldData.append(boundaryCellsConn[key], key)


def add_boundary_names_to_fielddata(polyData, boundaryData):
    """Add the names of the boundaries to the polydata as a field array.

    Parameters
    ----------
    polyData : vtkPolyData
        The polydata that the field data will be added to.
    boundaryData : dictionary
        Dictionary with the keys being the names of the boundaries.

    """
    boundaryNames = vtk.vtkStringArray()
    boundaryNames.SetName("boundaries")

    for boundary in boundaryData:
        boundaryNames.InsertNextValue(boundary)

    polyData.GetFieldData().AddArray(boundaryNames)


def create_boundary_polydata(patchBlocks, patchData, bounds, debug):

    patchFeatureEdgesFilter = vtk.vtkFeatureEdges()

    patchPolys = OrderedDict()

    for patchName in get_block_names(patchBlocks):
        if debug:
            print("    Extracting edges for patch", patchName)

        # Get the patch data
        patchBlockI = patchBlocks.GetBlock(get_block_index(patchBlocks,
                                                           patchName))

        # Extract feature edges
        patchFeatureEdgesFilter.SetInputData(patchBlockI)
        patchFeatureEdgesFilter.Update()
        patchFeatureEdgesData = patchFeatureEdgesFilter.GetOutput()

        boundaryIPoints = dsa.WrapDataObject(patchFeatureEdgesData).Points

        if debug:
            print("    Found", boundaryIPoints.shape[0], "edge points")

        b4Points = np.ones(boundaryIPoints.shape[0])*bounds[4]
        b5Points = np.ones(boundaryIPoints.shape[0])*bounds[5]
        # The patch is an x-y plane

        if (np.allclose(boundaryIPoints[:, 2], b4Points, atol=1e-6, rtol=1e-4) or
            np.allclose(boundaryIPoints[:, 2], b5Points, atol=1e-6, rtol=1e-4)):
            if debug:
                print("    The patch appears to lie in the x-y plane, ignoring it")
            continue
        else:

            # Select the points located at the boundary of the patch
            idx = np.abs(boundaryIPoints[:, 2] - patchData.GetPoint(0)[2]) < 3e-4

            if debug:
                print("    Found", np.sum(idx), "edge points with the same z-values as the seed patch")

            # Find ids of the cells in the x-y plane
            cellIds = vtk.vtkIntArray()
            for i in range(patchFeatureEdgesData.GetNumberOfCells()):

                # Get ids of the points of cell i
                pointIds = vtk.vtkIdList()
                patchFeatureEdgesData.GetCellPoints(i, pointIds)

                flag = 0
                for j in range(pointIds.GetNumberOfIds()):
                    if not idx[pointIds.GetId(j)]:
                        break
                    flag += 1

                if flag == 2:
                    cellIds.InsertNextValue(i)

            # Create the polydata and remove all cells but the ones in cellIds
            newPoly = vtk.vtkPolyData()
            newPoly.ShallowCopy(patchFeatureEdgesFilter.GetOutput())
            cellIds2 = vtk_to_numpy(cellIds)
            for i in range(newPoly.GetNumberOfCells()):
                if i not in cellIds2:
                    newPoly.DeleteCell(i)

            newPoly.RemoveDeletedCells()

            # Clean up redundant points
            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(newPoly)
            cleaner.Update()

            newPoly.ShallowCopy(cleaner.GetOutput())

            if debug:
                print("    Created polyData with", newPoly.GetNumberOfCells(), "cells")
            patchPolys[patchName] = newPoly
    return patchPolys


def assemble_multiblock(internalData, edgeDataDict):
    """Assemble the multiblock data set from the internal field and
    boundary data.

    Parameters
    ----------
    internalData : vtkPolyData
        The polydata with the internal field

    edgeDataDict : dictionary
        A dictionary with each entry being a vtkPolydata corresponding
        to one edge.

    Returns
    -------
    vtkMultiBlockDataSet
        The assembled dataset.

    """
    multiBlock = vtk.vtkMultiBlockDataSet()
    multiBlock.SetNumberOfBlocks(len(edgeDataDict) + 1)
    multiBlock.SetBlock(0, internalData)
    multiBlock.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(),
                                  "internalField")

    i = 1
    for boundary in edgeDataDict:
        multiBlock.SetBlock(i, edgeDataDict[boundary])
        multiBlock.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), boundary)
        i += 1

    return multiBlock


def create_parser():
    """Create a parse for the command line arguments.

    """
    parser = argparse.ArgumentParser(
        description="Script for averaging a 3D field along z.")

    parser.add_argument('--config',
                        type=str,
                        help='The config file.',
                        required=True)

    return parser


def zero_out_arrays(polyData):
    """Assign all point and cell arrays in given polyData to 0.

    Parameters
    ----------
    polyData : vtkPolyData
        The polyData.

    """
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
        nSamples = int(config["nSamples"])
        writePath = config["file"]
    except KeyError:
        print("ERROR: required parameter not specified in config file.")
        raise

    try:
        debug = bool(int(config["debug"]))
    except KeyError:
        debug = False
        pass

    try:
        time = float(config["time"])
    except KeyError:
        time = None
        pass

    try:
        dry = bool(int(config["dry"]))
    except KeyError:
        dry = False
        pass

    if debug:
        print("The debug switch is on")
        print("")
        print("The case path is "+casePath)
        print("The name of the seed patch is "+seedPatchName)
        print("The number of samples to be taken along z is "+str(nSamples))
        print("The produced filename will be "+writePath)

    # Case reader
    print("Reading")
    reader = read(casePath, time, debug)
    # Writer
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(writePath)

    caseData = reader.GetOutput()
    internalBlock = caseData.GetBlock(0)
    bounds = internalBlock.GetBounds()

    patchBlocks = caseData.GetBlock(1)
    seedPatchBlock = patchBlocks.GetBlock(get_block_index(patchBlocks,
                                                          seedPatchName))

    # The polyData for the 2d fields, copied from the seed patch
    internalData = vtk.vtkPolyData()
    internalData.ShallowCopy(seedPatchBlock)
    internalData.BuildLinks()

    print("Sampling and averaging internal field")
    zero_out_arrays(internalData)
    average_internal_field_data(internalBlock, internalData, nSamples, debug, dry)

    print("Creating boundary polyData")
    boundaryData = create_boundary_polydata(patchBlocks, internalData, bounds, debug)

    print("Averaging data for patches")
    average_patch_data(caseData, boundaryData, nSamples, bounds)
    add_boundary_names_to_fielddata(internalData, boundaryData)

    print("Marking boundary cells.")
    mark_boundary_cells(internalData, boundaryData, debug)

    print("Assembling multi-block structure")
    multiBlock = assemble_multiblock(internalData, boundaryData)

    print("Writing")
    writer.SetInputData(multiBlock)
    writer.Write()


if __name__ == '__main__':
    main()
