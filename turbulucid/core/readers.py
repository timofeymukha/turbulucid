# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import numpy as np
from collections import OrderedDict
import os
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersCore import vtkPointDataToCellData
import abc
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

__all__ = ["Reader", "LegacyReader", "mark_boundary_cells", "NativeReader",
           "XMLReader", "VTUReader"]


def mark_boundary_cells(internalData, boundaryDataDict):
    """Find the internal cell adjacent to each cell in the boundary
    data.

    Goes through each cell center for all boundary polydata and
    attempts to find the adjacent cell in the interalField.
    Creates a connectivity array and adds it as field data to the
    internal field polydata.

    """

    boundaryCellsConn = OrderedDict()

    for boundary in boundaryDataDict:
        boundaryDataI = boundaryDataDict[boundary]

        boundaryCellsConn[boundary] = \
            vtk_to_numpy(boundaryDataI.GetAttributes(vtkDataObject.CELL).
                         GetPedigreeIds())

    for key in boundaryCellsConn:
        wrappedData = dsa.WrapDataObject(internalData)
        wrappedData.FieldData.append(boundaryCellsConn[key], key)


class Reader():
    """Abstract base class for file readers."""
    __metaclass__ = abc.ABCMeta

    def __init__(self, fileName):
        if not os.path.exists(fileName):
            raise ValueError("ERROR: The file " + fileName + " does not exist")

    @abc.abstractproperty
    def vtkReader(self):
        pass

    @abc.abstractproperty
    def fileName(self):
        pass

    @abc.abstractproperty
    def data(self):
        pass

    def _clean(self, data):
        """Removes cells that have area less than 1e-10. Also runs the
        data through vtkCleanPolyData().

        """
        from vtkmodules.vtkFiltersVerdict import vtkMeshQuality
        from vtkmodules.vtkFiltersCore import vtkCleanPolyData

        data.BuildLinks()

        area = vtkMeshQuality()
        area.SetTriangleQualityMeasureToArea()
        area.SetInputData(data)
        area.Update()
        area = dsa.WrapDataObject(area.GetOutput()).CellData["Quality"]

        for i in range(data.GetNumberOfCells()):
            if area[i] < 1e-10:
                data.DeleteCell(i)

        data.RemoveDeletedCells()

        cleaner = vtkCleanPolyData()
        cleaner.SetInputData(data)
        cleaner.Update()
        return cleaner.GetOutput()

    def _transform(self, inputData):
        from vtkmodules.vtkCommonTransforms import vtkTransform
        from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

        transform = vtkTransform()
        meanNormal = self._compute_normal(inputData)

        axis = np.cross(meanNormal, [0, 0, 1])
        angle = np.rad2deg(np.arccos(np.dot(meanNormal, [0, 0, 1])))

        print(meanNormal, angle)

        # Do not rotate 180 degrees, no guarantee that it will be better
        # than no rotation at all
        if np.allclose([angle], [180]):
            angle = 0

        transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
        transform.Update()

        filter = vtkTransformPolyDataFilter()
        filter.SetInputData(inputData)
        filter.SetTransform(transform)
        filter.Update()

        data = filter.GetOutput()

        for i in range(data.GetNumberOfPoints()):
            p = data.GetPoints().GetPoint(i)
            x = p[0]
            y = p[1]
            data.GetPoints().SetPoint(i, [x, y, 0])

        #        zValue = filter.GetOutput().GetPoint(0)[-1]
        #        print(zValue)

        #        transform = vtkTransform()
        #        transform.Translate(0, 0, -zValue)
        #        transform.Update()

        #        data = vtkTransformPolyDataFilter()
        #        data.SetInputConnection(filter.GetOutputPort())
        #        data.SetTransform(transform)
        #        data.Update()
        return data

    def _compute_normal(self, inputData):
        from vtkmodules.vtkFiltersCore import vtkPolyDataNormals

        vtkNormals = vtkPolyDataNormals()
        vtkNormals.ComputeCellNormalsOn()
        vtkNormals.SetInputData(inputData)
        vtkNormals.Update()
        normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
        meanNormal = np.mean(normals, axis=0)
        meanNormal /= np.linalg.norm(meanNormal)
        return meanNormal

    def _extract_boundary_data(self, internalData):
        from vtkmodules.vtkFiltersCore import vtkFeatureEdges

        patchFeatureEdgesFilter = vtkFeatureEdges()
        patchFeatureEdgesFilter.FeatureEdgesOff()
        patchFeatureEdgesFilter.NonManifoldEdgesOff()
        patchFeatureEdgesFilter.ManifoldEdgesOff()

        patchFeatureEdgesFilter.SetInputData(internalData)
        patchFeatureEdgesFilter.Update()

        return patchFeatureEdgesFilter.GetOutput()

    def _assemble_multiblock_data(self, internalData, boundaryData):
        from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkCompositeDataSet
        from vtkmodules.vtkCommonCore import vtkStringArray

        multiBlock = vtkMultiBlockDataSet()
        multiBlock.SetNumberOfBlocks(2)
        multiBlock.SetBlock(0, internalData)
        multiBlock.GetMetaData(0).Set(vtkCompositeDataSet.NAME(),
                                      "internalField")
        multiBlock.SetBlock(1, boundaryData)
        multiBlock.GetMetaData(1).Set(vtkCompositeDataSet.NAME(),
                                      "boundary")

        boundaryNames = vtkStringArray()
        boundaryNames.SetName("boundaries")
        boundaryNames.InsertNextValue("boundary")

        internalData.GetFieldData().AddArray(boundaryNames)
        return multiBlock


class LegacyReader(Reader):
    """Reader for data in legacy VTK format, i.e. .vtk."""

    def __init__(self, filename, clean=False, pointData=False):
        from vtkmodules.vtkIOLegacy import vtkPolyDataReader

        Reader.__init__(self, filename)

        self._vtkReader = vtkPolyDataReader()
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()

        internalData = self._transform(self._vtkReader.GetOutput())
        if clean:
            internalData = self._clean(internalData)

        if pointData:
            interp = vtkPointDataToCellData()
            interp.SetInputData(internalData)
            interp.PassPointDataOff()
            interp.Update()
            internalData = interp.GetOutput()

        internalData.BuildLinks()

        n = internalData.GetNumberOfCells()
        pids = np.arange(n)

        internalData.GetAttributes(vtkDataObject.CELL).SetPedigreeIds(
            numpy_to_vtk(pids))

        boundaryData = self._extract_boundary_data(internalData)
        bDict = {'boundary': boundaryData}
        mark_boundary_cells(internalData, bDict)
        self._data = self._assemble_multiblock_data(internalData, boundaryData)

    @property
    def vtkReader(self):
        """The VTK reader for the data."""
        return self._vtkReader

    @property
    def fileName(self):
        """The path to the file with the data."""
        return self._fileName

    @property
    def data(self):
        """The read in data."""
        return self._data


class XMLReader(Reader):
    """Reader for data in XML polydata format, i.e. .vtp."""

    def __init__(self, filename, clean=False, pointData=False):
        Reader.__init__(self, filename)

        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader, vtkXMLStructuredGridReader
        from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter

        if ".vtu" in filename:
            self._vtkReader = vtkXMLUnstructuredGridReader()
        elif ".vtp" in filename:
            self._vtkReader = vtkXMLPolyDataReader()
        elif ".vts" in filename:
            self._vtkReader = vtkXMLStructuredGridReader()
        else:
            raise NotImplementedError
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()

        if (".vtu" in filename) or (".vts" in filename):
            polydata = vtkDataSetSurfaceFilter()
            polydata.SetInputData(self._vtkReader.GetOutput())
            polydata.Update()
            readData = polydata.GetOutput()
        else:
            readData = self._vtkReader.GetOutput()

        internalData = self._transform(readData)
        if clean:
            internalData = self._clean(internalData)

        if pointData:
            interp = vtkPointDataToCellData()
            interp.SetInputData(internalData)
            interp.PassPointDataOff()
            interp.Update()
            internalData = interp.GetOutput()

        internalData.BuildLinks()

        n = internalData.GetNumberOfCells()
        pids = np.arange(n)

        internalData.GetAttributes(vtkDataObject.CELL).SetPedigreeIds(
            numpy_to_vtk(pids))

        boundaryData = self._extract_boundary_data(internalData)
        bDict = {'boundary': boundaryData}
        mark_boundary_cells(internalData, bDict)
        self._data = self._assemble_multiblock_data(internalData, boundaryData)
    @property
    def vtkReader(self):
        """The VTK reader for the data."""
        return self._vtkReader

    @property
    def fileName(self):
        """The path to the file with the data."""
        return self._fileName

    @property
    def data(self):
        """The read in data."""
        return self._data

class VTUReader(Reader):
    """Reader for data in XML unstructured grid format, i.e. .vtu."""

    def __init__(self, filename, clean=False, pointData=False):
        Reader.__init__(self, filename)

        self._vtkReader = vtk.vtkXMLUnstructuredGridReader()
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()

        internalData = self._vtkReader.GetOutput()
        self._data = self._vtkReader.GetOutput()

    #        internalData = self._transform()
#        if clean:
#            internalData = self._clean(internalData)

        if pointData:
            interp = vtk.vtkPointDataToCellData()
            interp.SetInputData(internalData)
            interp.PassPointDataOff()
            interp.Update()
            internalData = interp.GetOutput()

        internalData.BuildLinks()

        n = internalData.GetNumberOfCells()
        pids = np.arange(n)

        internalData.GetAttributes(vtk.vtkDataObject.CELL).SetPedigreeIds(
            numpy_to_vtk(pids))

        boundaryData = self._extract_boundary_data(internalData)
        bDict = {'boundary': boundaryData}
        mark_boundary_cells(internalData, bDict)
        self._data = self._assemble_multiblock_data(internalData, boundaryData)
    @property
    def vtkReader(self):
        """The VTK reader for the data."""
        return self._vtkReader

    @property
    def fileName(self):
        """The path to the file with the data."""
        return self._fileName

    @property
    def data(self):
        """The read in data."""
        return self._data



class NativeReader(Reader):
    """Reader for native turbulucid format."""

    def __init__(self, fileName):
        from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataReader

        Reader.__init__(self, fileName)

        self._vtkReader = vtkXMLMultiBlockDataReader()
        self._fileName = fileName

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()
        self._data = self._vtkReader.GetOutput()

        for i in range(self._data.GetNumberOfBlocks()):
            points = dsa.WrapDataObject(self._data.GetBlock(i)).Points
            points[:, 2] = 0

    @property
    def vtkReader(self):
        """The VTK reader for the data."""
        return self._vtkReader

    @property
    def fileName(self):
        """The path to the file with the data."""
        return self._fileName

    @property
    def data(self):
        """The read in data."""
        return self._data

    def _compute_normal(self):
        from vtkmodules.vtkFiltersCore import vtkPolyDataNormals

        vtkNormals = vtkPolyDataNormals()
        vtkNormals.ComputeCellNormalsOn()
        vtkNormals.SetInputData(self._vtkReader.GetOutput().GetBlock(0))
        vtkNormals.Update()
        normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
        meanNormal = np.mean(normals, axis=0)
        meanNormal /= np.linalg.norm(meanNormal)
        return meanNormal
