# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import os

import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.numpy_support import numpy_to_vtk
from .readers import NativeReader, LegacyReader, XMLReader
import numpy as np

__all__ = ["Case"]


class Case:
    """A class representing a simulation case.

    """

    def __init__(self, fileName, clean=False, pointData=False):
        """
        Create Case from file.

        Parameters
        ----------
        fileName : str
            The file to be read in. Should be data in VTK format.

        clean : bool
            Whether to attempt to clean the data of redundant cells.

        """
        self.fileName = os.fspath(fileName)

        # Read in the data
        self._blockData = self.read(clean, pointData)
        self._validate_structure()

    @property
    def blockData(self):
        """vtkMultiBlockDataSet : the multiblock data assembled by the
        reader.

        """
        return self._blockData

    @property
    def vtkData(self):
        """wrapped PolyData : The actual data read by the reader."""

        return dsa.WrapDataObject(self._blockData.GetBlock(0))

    @property
    def cellCentres(self):
        """wrapped VTKArray : the cell centres of the read data """

        cellCentres = vtk.vtkCellCenters()
        cellCentres.SetInputData(self._blockData.GetBlock(0))
        cellCentres.Update()
        points = dsa.WrapDataObject(cellCentres.GetOutput()).GetPoints()
        return np.array(points[:, :2])

    @property
    def boundaries(self):
        """list : A list of names of the boundaries present the case."""

        return self._fill_boundary_list()

    @property
    def bounds(self):
        """tuple : (min(x), max(x), min(y), max(y))."""

        return self.vtkData.VTKObject.GetBounds()[:4]

    @property
    def fields(self):
        """list of str: The names of the fields present in the case."""

        return list(self.vtkData.CellData.keys())

    @property
    def xlim(self):
        """list of two floats: The x limits that cover the
        geometry of the case, plus small a margin.

        """
        return self._compute_plot_limits()[0]

    @property
    def ylim(self):
        """list of two floats: The y limits that cover the
        geometry of the case, plus a small margin.

        """
        return self._compute_plot_limits()[1]

    def _validate_structure(self):
        """Validate the multiblock structure required by the public API."""
        if self._blockData is None or self._blockData.GetNumberOfBlocks() == 0:
            raise ValueError("The case does not contain an internal data block.")
        if self._blockData.GetBlock(0) is None:
            raise ValueError("The case's internal data block is empty.")
        if "boundaries" not in self.vtkData.FieldData.keys():
            raise ValueError(
                "The internal data block has no 'boundaries' field metadata."
            )

        expectedBlocks = len(self._fill_boundary_list()) + 1
        if self._blockData.GetNumberOfBlocks() != expectedBlocks:
            raise ValueError(
                "The number of boundary blocks does not match the boundary metadata."
            )

    def _fill_boundary_list(self):
        fieldData = self.vtkData.FieldData['boundaries']
        boundaryList = []

        for i in range(fieldData.GetNumberOfValues()):
            boundaryList.append(fieldData.GetValue(i))

        return boundaryList

    def __getitem__(self, item):
        """Return a cell array by name.

        Parameters
        ----------
        item : string
            The name of the cell array.

        Returns
        -------
        ndarray
            Array of values of the requested field.

        """
        self._validate_field_name(item)
        if item not in self.fields:
            raise ValueError(f"Field {item} not present in the case.")

        return np.copy(np.array((self.vtkData.CellData[item])))

    def __setitem__(self, item, values):
        """Add another internal field to the case.

        Parameters
        ----------
        item : string
            The name of the cell array.
        values : ndarray
            The values of the field.

        """
        self._validate_field_name(item)
        values = np.asarray(values)
        if values.ndim == 0 or values.shape[0] != self.vtkData.GetNumberOfCells():
            raise ValueError("The dimensionality of the provided field "
                             "does not match that of the case.")
        if (not np.issubdtype(values.dtype, np.number) or
                np.issubdtype(values.dtype, np.complexfloating)):
            raise TypeError("Field values must be real numbers.")
        if not np.all(np.isfinite(values)):
            raise ValueError("Field values must be finite.")
        if values.ndim > 3 or (
                values.ndim == 3 and values.shape[1:] != (3, 3)):
            raise ValueError(
                "Fields must be scalar, component arrays, or 3-by-3 tensors."
            )
        if any(size == 0 for size in values.shape[1:]):
            raise ValueError("Field component dimensions must not be empty.")

        self._add_cell_array(self.vtkData.VTKObject, item, values)

        # Boundary values for derived fields are copied from adjacent cells.
        for boundary in self.boundaries:
            boundaryCellIds = np.asarray(
                self.vtkData.FieldData[boundary], dtype=np.intp
            )
            boundaryValues = values[boundaryCellIds, ...]
            block = self.extract_block_by_name(boundary)
            self._add_cell_array(block, item, boundaryValues)

    def __delitem__(self, item):
        """Delete an internal field form the case.

        Parameters
        ----------
        item : str
            Name of the field to delete.

        """
        self._validate_field_name(item)
        if item not in self.fields:
            raise ValueError(f"Field {item} not present in the case.")

        self.vtkData.VTKObject.GetCellData().RemoveArray(item)

        for boundary in self.boundaries:
            block = self.extract_block_by_name(boundary)
            block.GetCellData().RemoveArray(item)

    @staticmethod
    def _validate_field_name(name):
        if not isinstance(name, str):
            raise TypeError("Field name must be a string.")
        if not name:
            raise ValueError("Field name must not be empty.")

    @staticmethod
    def _add_cell_array(data, name, values):
        """Add a copied NumPy array to a VTK dataset's cell data."""
        if values.ndim == 3:
            # VTK stores tensor components in column-major matrix order.
            values = values.transpose(0, 2, 1).reshape(values.shape[0], 9)
        elif values.ndim == 2:
            values = values.reshape(values.shape[0], -1)
        values = np.ascontiguousarray(values)
        vtkValues = numpy_to_vtk(values, deep=True)
        vtkValues.SetName(name)
        data.GetCellData().AddArray(vtkValues)

    def _compute_plot_limits(self):
        """ Compute xlim and ylim."""

        minX = self.bounds[0]
        maxX = self.bounds[1]
        minY = self.bounds[2]
        maxY = self.bounds[3]

        marginX = (maxX - minX)/60
        marginY = (maxY - minY)/60

        return (np.array([minX - marginX, maxX + marginX]),
                np.array([minY - marginY, maxY + marginY]))

    def _transform(self, transform):
        """Transform the geometry according to a vtkTransform filter."""

        from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

        # Transform the internal field
        filter = vtkTransformPolyDataFilter()
        filter.SetInputData(self.blockData.GetBlock(0))
        filter.SetTransform(transform)
        filter.Update()

        self._blockData.SetBlock(0, filter.GetOutput())

        # Transform boundary data
        i = 1
        for boundary in self.boundaries:
            filter = vtkTransformPolyDataFilter()
            filter.SetTransform(transform)
            filter.SetInputData(self.blockData.GetBlock(i))
            filter.Update()
            self.blockData.SetBlock(i, filter.GetOutput())
            i += 1

    def translate(self, dx, dy):
        """Translate the geometry of the case.

        Parameters
        ----------
        dx : float
            The translation along the x axis.
        dy : float
            The translation along the y axis.

        """
        dx = self._finite_scalar(dx, "dx")
        dy = self._finite_scalar(dy, "dy")
        transform = vtk.vtkTransform()
        transform.Translate(dx, dy, 0)
        transform.Update()

        self._transform(transform)

    def scale(self, scaleX, scaleY):
        """Scale the geometry of the case.

        The coordinates get divided by the scaling factors.

        Parameters
        ----------
        scaleX : float
            The scaling factor along x.
        scaleY : float
            The scaling factor along y.

        """
        scaleX = self._finite_scalar(scaleX, "scaleX", nonzero=True)
        scaleY = self._finite_scalar(scaleY, "scaleY", nonzero=True)
        transform = vtk.vtkTransform()
        transform.Scale(1/scaleX, 1/scaleY, 1)
        transform.Update()
        self._transform(transform)

    def rotate(self, angle):
        """Rotate the geometry of the case around the z axis.

        Parameters
        ----------
        dx : angle
            Rotation angle in degrees.

        """
        angle = self._finite_scalar(angle, "angle")
        axis = [0, 0, 1]
        transform = vtk.vtkTransform()
        transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
        transform.Update()
        self._transform(transform)

    @staticmethod
    def _finite_scalar(value, name, nonzero=False):
        if not np.isscalar(value):
            raise TypeError(f"{name} must be a scalar.")
        try:
            value = float(value)
        except (TypeError, ValueError) as error:
            raise TypeError(f"{name} must be a real number.") from error
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite.")
        if nonzero and value == 0:
            raise ValueError(f"{name} must be nonzero.")
        return value

    def boundary_cell_data(self, boundary, sort=None):
        """Return cell-centre coordinates and data from cells adjacent
        to a specific boundary.

        Parameters
        ----------
        boundary : str
            The name of the boundary.
        sort : {None, 'x', 'y'}, optional
            Whether to sort the data along a coordinate. Use 'x' and
            'y' to sort along x and y, respectively. Default is no
            sorting.

        Returns
        -------
            Two ndarrays

        """
        self._validate_boundary(boundary)
        self._validate_sort(sort)
        cellIds = np.asarray(self.vtkData.FieldData[boundary], dtype=np.intp)
        points = self.cellCentres[cellIds, :]
        data = {
            field: self[field][cellIds, ...]
            for field in self.fields
        }

        if sort is None:
            return points, data
        if sort == "x":
            ind = np.argsort(points[:, 0])
        else:
            ind = np.argsort(points[:, 1])

        points = points[ind]

        for key in data:
            data[key] = data[key][ind, ...]

        return points, data

    def extract_block_by_name(self, name):
        """Extract a block from the case by a given name."""
        self._validate_boundary(name)
        return self._blockData.GetBlock(self.boundaries.index(name) + 1)

    def _validate_boundary(self, boundary):
        if not isinstance(boundary, str):
            raise TypeError("Boundary name must be a string.")
        if boundary not in self.boundaries:
            raise ValueError(f"Boundary {boundary} not present in the case.")

    @staticmethod
    def _validate_sort(sort):
        if sort not in {None, "x", "y"}:
            raise ValueError("sort should be 'x', 'y', or None.")

    def boundary_data(self, boundary, sort=None):
        """Return cell-center coordinates and data from a boundary.

        Parameters
        ----------
        boundary : str
            The name of the boundary.
        sort : str
            Whether to sort the data along a coordinate. Use "x" and
            "y" to sort along x and y, respectively. Default is no
            sorting.

        Returns
        -------
        Two ndarrays
            The coordinates of the boundary face centres.
            The corresponding data.
        """

        self._validate_sort(sort)
        blockData = self.extract_block_by_name(boundary)

        cCenters = vtk.vtkCellCenters()
        cCenters.SetInputData(blockData)
        cCenters.Update()

        points = np.array(dsa.WrapDataObject(cCenters.GetOutput()).Points)
        dataVTK = dsa.WrapDataObject(blockData).CellData

        data = {}
        for key in dataVTK.keys():
            data[key] = np.array(dataVTK[key])

        if sort is None:
            return points[:, [0, 1]], data
        if sort == "x":
            ind = np.argsort(points[:, 0])
        else:
            ind = np.argsort(points[:, 1])

        points = points[ind]

        for key in data:
            data[key] = data[key][ind]

        return points[:, [0, 1]], data

    def read(self, clean, pointData):
        """Read in the data from a file.

        Parameters
        ----------
        clean : bool
            Whether to attempt cleaning the case of degenerate cells upon
            read.
        pointData : bool
            Whether the file contains point data instead of cell data.
            Cell data will be computed by interpolation.

        Raises
        ------
        ValueError
            If the provided file does not exist.

        """
        fileName = self.fileName

        fileExt = os.path.splitext(fileName)[1].lower()

        if fileExt == ".vtm":
            reader = NativeReader(fileName)
            return reader.data
        elif fileExt == ".vtk":
            return LegacyReader(fileName, clean=clean,
                                pointData=pointData).data
        elif fileExt in [".vtu", ".vtp", ".vts"]:
            return XMLReader(fileName, clean=clean, pointData=pointData).data
        else:
            raise ValueError("Unsupported file format.", fileName, fileExt)

    def write(self, writePath):
        """Save the case to a .vtm format.

        Parameters
        ----------
        writePath : str
            The name of the file.

        Raises
        ------
        OSError
            If the data could not be written to the given path.

        """
        from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataWriter
        from vtkmodules.vtkCommonMisc import vtkErrorCode

        writePath = os.fspath(writePath)

        writer = vtkXMLMultiBlockDataWriter()
        writer.SetFileName(writePath)
        writer.SetInputData(self._blockData)
        writer.Write()

        # The writer does not raise, and its return value stays 1 even on
        # failure, so the error code is the only reliable status.
        errorCode = writer.GetErrorCode()
        if errorCode != vtkErrorCode.NoError:
            raise OSError(
                f"Could not write the case to {writePath}: "
                + vtkErrorCode.GetStringFromErrorCode(errorCode)
            )
