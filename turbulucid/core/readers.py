# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import abc
import os
import warnings
from collections import OrderedDict

import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersCore import vtkPointDataToCellData

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

    wrappedData = dsa.WrapDataObject(internalData)
    for key in boundaryCellsConn:
        wrappedData.FieldData.append(boundaryCellsConn[key], key)


class Reader(abc.ABC):
    """Abstract base class for file readers."""

    def __init__(self, fileName):
        fileName = os.fspath(fileName)
        if not os.path.isfile(fileName):
            raise FileNotFoundError(f"No such file: {fileName}")

    @property
    @abc.abstractmethod
    def vtkReader(self):
        pass

    @property
    @abc.abstractmethod
    def fileName(self):
        pass

    @property
    @abc.abstractmethod
    def data(self):
        pass

    def _clean(self, data):
        """Removes cells that have area less than 1e-10. Also runs the
        data through vtkCleanPolyData().

        Only triangles and quadrilaterals are considered.  vtkMeshQuality
        reports NaN for general polygons, which never compares less than the
        threshold, so such cells are always kept.

        """
        from vtkmodules.vtkFiltersCore import vtkCleanPolyData
        from vtkmodules.vtkFiltersVerdict import vtkMeshQuality

        data.BuildLinks()

        area = vtkMeshQuality()
        area.SetTriangleQualityMeasureToArea()
        # Without this the quads are scored with the default measure, which
        # reports 1e30 for degenerate cells and so never triggers the check.
        area.SetQuadQualityMeasureToArea()
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

        targetNormal = np.array([0.0, 0.0, 1.0])
        axis = np.cross(meanNormal, targetNormal)
        sinAngle = np.linalg.norm(axis)
        cosAngle = np.clip(np.dot(meanNormal, targetNormal), -1.0, 1.0)

        if sinAngle > 1e-12:
            axis /= sinAngle
            angle = np.rad2deg(np.arctan2(sinAngle, cosAngle))
            transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
        transform.Update()

        transformFilter = vtkTransformPolyDataFilter()
        transformFilter.SetInputData(inputData)
        transformFilter.SetTransform(transform)
        transformFilter.Update()

        data = transformFilter.GetOutput()

        points = dsa.WrapDataObject(data).Points
        points[:, 2] = 0
        return data

    def _compute_normal(self, inputData):
        """Return the normal of a validated best-fit plane.

        Point-based plane fitting is independent of polygon winding, unlike
        averaging cell normals.  The validation prevents silently flattening
        genuinely three-dimensional or degenerate input.
        """
        vtkPoints = inputData.GetPoints()
        if vtkPoints is None or vtkPoints.GetNumberOfPoints() < 3:
            raise ValueError("The input must contain at least three points.")

        points = np.asarray(
            vtk_to_numpy(vtkPoints.GetData()), dtype=float
        )
        if not np.all(np.isfinite(points)):
            raise ValueError("The input geometry contains non-finite points.")

        centredPoints = points - np.mean(points, axis=0)
        _, singularValues, directions = np.linalg.svd(
            centredPoints, full_matrices=False
        )
        scale = singularValues[0]
        if scale == 0 or singularValues[1] <= scale * 1e-6:
            raise ValueError("The input geometry is degenerate or collinear.")

        normal = directions[-1]
        maxDistance = np.max(np.abs(centredPoints @ normal))
        if maxDistance > max(1e-10, scale * 1e-6):
            raise ValueError("The input geometry is not planar.")

        # SVD leaves the sign undetermined. Prefer the orientation requiring
        # at most a 90-degree rotation towards the positive z axis.
        if normal[2] < 0:
            normal = -normal
        elif np.isclose(normal[2], 0):
            dominant = np.argmax(np.abs(normal[:2]))
            if normal[dominant] < 0:
                normal = -normal

        return normal

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
        from vtkmodules.vtkCommonCore import vtkStringArray
        from vtkmodules.vtkCommonDataModel import (
            vtkCompositeDataSet,
            vtkMultiBlockDataSet,
        )

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

    def _build(self, readData, clean, pointData):
        """Turn raw polydata from a VTK reader into the multiblock layout.

        The steps are identical for every single-block input format, so
        the concrete readers only have to produce the polydata.

        Parameters
        ----------
        readData : vtkPolyData
            The geometry as it came out of the format-specific reader.
        clean : bool
            Whether to remove degenerate cells.
        pointData : bool
            Whether the data is point data that must be interpolated to
            the cells.

        Returns
        -------
        vtkMultiBlockDataSet
            The assembled dataset.

        """
        internalData = self._transform(readData)

        if clean:
            internalData = self._clean(internalData)

        if pointData:
            interpolator = vtkPointDataToCellData()
            interpolator.SetInputData(internalData)
            interpolator.PassPointDataOff()
            interpolator.Update()
            internalData = interpolator.GetOutput()

        internalData.BuildLinks()

        # The pedigree ids let the boundary blocks refer back to the
        # internal cell each boundary face belongs to.
        pedigreeIds = np.arange(internalData.GetNumberOfCells())
        internalData.GetAttributes(vtkDataObject.CELL).SetPedigreeIds(
            numpy_to_vtk(pedigreeIds, deep=True))

        boundaryData = self._extract_boundary_data(internalData)
        mark_boundary_cells(internalData, {"boundary": boundaryData})
        return self._assemble_multiblock_data(internalData, boundaryData)


class LegacyReader(Reader):
    """Reader for data in legacy VTK format, i.e. .vtk."""

    def __init__(self, filename, clean=False, pointData=False):
        from vtkmodules.vtkIOLegacy import vtkPolyDataReader

        super().__init__(filename)

        self._vtkReader = vtkPolyDataReader()
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()

        self._data = self._build(self._vtkReader.GetOutput(), clean, pointData)

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
    """Reader for XML VTK polydata, unstructured, and structured grids."""

    def __init__(self, filename, clean=False, pointData=False):
        super().__init__(filename)

        from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
        from vtkmodules.vtkIOXML import (
            vtkXMLPolyDataReader,
            vtkXMLStructuredGridReader,
            vtkXMLUnstructuredGridReader,
        )

        extension = os.path.splitext(filename)[1].lower()
        readerTypes = {
            ".vtp": vtkXMLPolyDataReader,
            ".vts": vtkXMLStructuredGridReader,
            ".vtu": vtkXMLUnstructuredGridReader,
        }
        try:
            self._vtkReader = readerTypes[extension]()
        except KeyError as error:
            raise ValueError(
                f"Unsupported XML VTK file extension: {extension or '<none>'}"
            ) from error
        self._fileName = filename

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()

        if extension in {".vtu", ".vts"}:
            polydata = vtkDataSetSurfaceFilter()
            polydata.SetInputData(self._vtkReader.GetOutput())
            polydata.Update()
            readData = polydata.GetOutput()
        else:
            readData = self._vtkReader.GetOutput()

        self._data = self._build(readData, clean, pointData)

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


class VTUReader(XMLReader):
    """Deprecated compatibility wrapper for :class:`XMLReader`."""

    def __init__(self, filename, clean=False, pointData=False):
        warnings.warn(
            "VTUReader is deprecated; use XMLReader instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(filename, clean=clean, pointData=pointData)


class NativeReader(Reader):
    """Reader for native turbulucid format."""

    def __init__(self, fileName):
        from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataReader

        super().__init__(fileName)

        self._vtkReader = vtkXMLMultiBlockDataReader()
        self._fileName = fileName

        self._vtkReader.SetFileName(self._fileName)
        self._vtkReader.Update()
        self._data = self._vtkReader.GetOutput()

        if self._data is None or self._data.GetNumberOfBlocks() == 0:
            raise ValueError(
                f"{self._fileName} does not contain any data blocks.")

        for i in range(self._data.GetNumberOfBlocks()):
            block = self._data.GetBlock(i)
            if block is None:
                raise ValueError(f"Block {i} of {self._fileName} is empty.")

            # The geometry is planar by construction, so flatten away any
            # round-off left in the stored z coordinates.
            points = dsa.WrapDataObject(block).Points
            if points is not None and len(points) > 0:
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
