from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
from collections import OrderedDict
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
import abc

__all__ = ["Reader", "LegacyReader"]


class Reader():
    """Abstract base class for file readers."""
    __metaclass__ = abc.ABCMeta

    @property
    @abc.abstractmethod
    def vtkReader(self):
        """Return the associated VTK reader."""


class LegacyReader:
    """Reader for data in legacy VTK format, i.e. .vtk."""

    def __init__(self, filename):

        self._vtkReader = vtk.vtkPolyDataReader()
        self._filename = filename

        self._vtkReader.SetFileName(self._filename)
        self._vtkReader.Update()
        internalData = self.__transform()
        boundaryData = self.__extract_boundary_data(internalData)
        bDict = {'boundary': boundaryData}
        self.mark_boundary_cells(internalData, bDict)
        self._data = self.__assemble_multiblock_data(internalData, boundaryData)

    @property
    def vtkReader(self):
        return self._vtkReader

    @property
    def filename(self):
        return self._filename

    @property
    def data(self):
        return self._data

    def __compute_normal(self):
        vtkNormals = vtk.vtkPolyDataNormals()
        vtkNormals.ComputeCellNormalsOn()
        vtkNormals.SetInputConnection(self._vtkReader.GetOutputPort())
        vtkNormals.Update()
        normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
        meanNormal = np.mean(normals, axis=0)
        meanNormal /= np.linalg.norm(meanNormal)
        return meanNormal

    def __transform(self):
        transform = vtk.vtkTransform()

        meanNormal = self.__compute_normal()

        axis = np.cross(meanNormal, [0, 0, 1])
        angle = np.rad2deg(np.arccos(np.dot(meanNormal, [0, 0, 1])))

        transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
        transform.Update()

        filter = vtk.vtkTransformPolyDataFilter()
        filter.SetInputConnection(self.vtkReader.GetOutputPort())
        filter.SetTransform(transform)
        filter.Update()

        zValue = filter.GetOutput().GetPoint(0)[-1]

        transform = vtk.vtkTransform()
        transform.Translate(0, 0, -zValue)
        transform.Update()

        data = vtk.vtkTransformPolyDataFilter()
        data.SetInputConnection(filter.GetOutputPort())
        data.SetTransform(transform)
        data.Update()
        return data.GetOutput()

    def __extract_boundary_data(self, internalData):
        patchFeatureEdgesFilter = vtk.vtkFeatureEdges()
        patchFeatureEdgesFilter.FeatureEdgesOff()
        patchFeatureEdgesFilter.NonManifoldEdgesOff()
        patchFeatureEdgesFilter.ManifoldEdgesOff()

        patchFeatureEdgesFilter.SetInputData(internalData)
        patchFeatureEdgesFilter.Update()
        return patchFeatureEdgesFilter.GetOutput()

    def __assemble_multiblock_data(self, internalData, boundaryData):

        multiBlock = vtk.vtkMultiBlockDataSet()
        multiBlock.SetNumberOfBlocks(2)
        multiBlock.SetBlock(0, internalData)
        multiBlock.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(),
                                      "internalField")
        multiBlock.SetBlock(1, boundaryData)
        multiBlock.GetMetaData(1).Set(vtk.vtkCompositeDataSet.NAME(),
                                      "boundary")

        boundaryNames = vtk.vtkStringArray()
        boundaryNames.SetName("boundaries")
        boundaryNames.InsertNextValue("boundary")

        internalData.GetFieldData().AddArray(boundaryNames)
        return multiBlock

    def mark_boundary_cells(self, patchData, patchPolys):
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

            polyI = patchPolys[boundary]
            cellCenters.SetInputData(polyI)
            cellCenters.Update()

            points = dsa.WrapDataObject(cellCenters.GetOutput()).Points

            cell = vtk.vtkGenericCell()
            tol2 = 0.0
            pcoords = [0, 0, 0]
            weights = []

            for i in range(points.shape[0]):
                pointI = points[i, :]
                foundCellId = locator.FindCell(pointI, tol2, cell, pcoords, weights)
                boundaryCellsConn[boundary][i] = foundCellId

                if foundCellId == -1:
                    print("Failed to find adjacent cell for boundary point",
                          pointI, "on boundary", boundary)
                    print("    Attempting with slow algorithm based on minimum"
                          " distance")
                    foundCellId, distance = get_closest_cell(pointI, patchData)

                    print("    Found cell with id", foundCellId, "located",
                          distance, "away.")
                    boundaryCellsConn[boundary][i] = foundCellId

        for key in boundaryCellsConn:
            if np.any(boundaryCellsConn[key] == -1):
                print("ERROR: some connectivity not established for boundary "+key)









