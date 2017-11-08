from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from collections import OrderedDict
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
import abc
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy

__all__ = ["Reader", "LegacyReader", "mark_boundary_cells", "NativeReader"]


def mark_boundary_cells(internalData, boundaryDataDict):
    """Find the internal cell adjacent to each cell in the boundary
    data.

    Goes through each cell center for all boundary polydata and
    attempts to find the adjacent cell in the interalField.
    Creates a connectivity array and adds it as field data to the
    internal field polydata.

    """
    boundaryCellsConn = OrderedDict()
    cellCenters = vtk.vtkCellCenters()

    for boundary in boundaryDataDict:
        boundaryCellsConn[boundary] = -1*np.ones(boundaryDataDict[boundary].GetNumberOfCells(), dtype=np.int32)

    for boundary in boundaryDataDict:

        boundaryDataI = boundaryDataDict[boundary]

        boundaryCellsConn[boundary] = \
            vtk_to_numpy(boundaryDataI.GetAttributes(vtk.vtkDataObject.CELL).
                         GetPedigreeIds())

    for key in boundaryCellsConn:
        if np.any(boundaryCellsConn[key] == -1):
            print("ERROR: some connectivity not established for boundary "+key)
        wrappedData = dsa.WrapDataObject(internalData)
        wrappedData.FieldData.append(boundaryCellsConn[key], key)


class Reader():
    """Abstract base class for file readers."""
    __metaclass__ = abc.ABCMeta

    def __init__(self, fileName):
        if not os.path.exists(fileName):
            raise ValueError("ERROR: The file "+fileName+" does not exist")
        pass

    @abc.abstractproperty
    def vtkReader(self):
        pass

    @abc.abstractproperty
    def fileName(self):
        pass

    @abc.abstractproperty
    def data(self):
        pass

    def _transform(self):
        transform = vtk.vtkTransform()

        meanNormal = self._compute_normal()

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

    def _compute_normal(self):
        vtkNormals = vtk.vtkPolyDataNormals()
        vtkNormals.ComputeCellNormalsOn()
        vtkNormals.SetInputConnection(self._vtkReader.GetOutputPort())
        vtkNormals.Update()
        normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
        meanNormal = np.mean(normals, axis=0)
        meanNormal /= np.linalg.norm(meanNormal)
        return meanNormal

    def _extract_boundary_data(self, internalData):
        patchFeatureEdgesFilter = vtk.vtkFeatureEdges()
        patchFeatureEdgesFilter.FeatureEdgesOff()
        patchFeatureEdgesFilter.NonManifoldEdgesOff()
        patchFeatureEdgesFilter.ManifoldEdgesOff()

        patchFeatureEdgesFilter.SetInputData(internalData)
        patchFeatureEdgesFilter.Update()
        return patchFeatureEdgesFilter.GetOutput()

    def _assemble_multiblock_data(self, internalData, boundaryData):

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


class LegacyReader(Reader):
    """Reader for data in legacy VTK format, i.e. .vtk."""

    def __init__(self, filename):
        Reader.__init__(self, filename)

        self._vtkReader = vtk.vtkPolyDataReader()
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()
        internalData = self._transform()
        boundaryData = self._extract_boundary_data(internalData)
        bDict = {'boundary': boundaryData}
        mark_boundary_cells(internalData, bDict)
        self._data = self._assemble_multiblock_data(internalData, boundaryData)

    @property
    def vtkReader(self):
        return self._vtkReader

    @property
    def fileName(self):
        return self._fileName

    @property
    def data(self):
        return self._data


class NativeReader(Reader):
    """Reader for native turbulucid format."""

    def __init__(self, fileName):
        Reader.__init__(self, fileName)

        self._vtkReader = vtk.vtkXMLMultiBlockDataReader()
        self._fileName = fileName

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()
        self._data = self._vtkReader.GetOutput()

    @property
    def vtkReader(self):
        return self._vtkReader

    @property
    def fileName(self):
        return self._fileName

    @property
    def data(self):
        return self._data











