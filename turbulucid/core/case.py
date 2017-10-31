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

__all__ = ["Case"]


class Case:
    """A class representing a simulation case.

    """
    def __init__(self, fileName):
        """
        Create Case from file.

        Parameters
        ----------
        fileName : str
            The file to be read in. Should be data in VTK format.

        """

        self.fileName = fileName

        # Read in the data
        self._blockData = self.read(fileName)

        # Compute the cell-centres
        self._cellCentres = vtk.vtkCellCenters()
        self._cellCentres.SetInputData(self._blockData.GetBlock(0))
        self._cellCentres.Update()
        self._cellCentres =\
            dsa.WrapDataObject(self._cellCentres.GetOutput()).GetPoints()
        self._cellCentres = np.array(self._cellCentres[:, :2])

        self._vtkData = dsa.WrapDataObject(self._blockData.GetBlock(0))

        self._zValue = self._vtkData.Points[0, 2]

        self._boundaries = self._fill_boundary_list()

        self._bounds = self._vtkData.VTKObject.GetBounds()

        self._fields = self._vtkData.CellData.keys()

    def read(self, fileName):
        """Read in the data from a file.

        Parameters
        ----------
        fileName : str
            The path to the file with the data.

        Raises
        ------
        ValueError
            If the provided file does not exist.

        """
        if not os.path.exists(fileName):
            raise ValueError("ERROR: The file "+fileName+" does not exist")

        fileExt = os.path.splitext(fileName)[1]

        if fileExt == ".vtm":
            reader = vtk.vtkXMLMultiBlockDataReader()

        reader.SetFileName(fileName)
        reader.Update()

        return reader.GetOutput()

    @property
    def vtkData(self):
        """wrapped PolyData : The actual data read by the reader."""

        return self._vtkData

    @property
    def cellCentres(self):
        """wrapped VTKArray : the cell centres of the read data """

        return self._cellCentres

    @property
    def zValue(self):
        """float : the value of z for the considered geometry."""

        return self._zValue

    @property
    def boundaries(self):
        """list : A list of names of the boundaries present the case."""

        return self._boundaries

    @property
    def bounds(self):
        """tuple : (min(x), max(x), min(y), max(y), min(z), max(z))."""

        return self._bounds

    @property
    def fields(self):
        """list of str: The names of the fields present in the case."""

        return self._fields

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
        return np.copy(vtk_to_numpy(self.vtkData.CellData[item]))

    def __setitem__(self, item, values):
        """Add another internal field to the case.

        Parameters
        ----------
        item : string
            The name of the cell array.
        values : ndarray
            The values of the field.

        """
        if values.shape[0] != self[self.fields[0]].shape[0]:
            raise ValueError("The dimensionality of the provided field "
                             "does not match that of the case.")

        cellData = self._vtkData.VTKObject.GetCellData()
        valuesVtk = vtk.vtkDoubleArray()

        if np.ndim(values) > 1:
            valuesVtk.SetNumberOfComponents(values.shape[1])
            valuesVtk.SetNumberOfTuples(values.shape[0])
            for i in range(values.shape[0]):
                valuesVtk.SetTupleValue(i, values[i, :])
        else:
            valuesVtk.SetNumberOfComponents(1)
            valuesVtk.SetNumberOfValues(values.shape[0])
            for i in range(values.shape[0]):
                valuesVtk.SetValue(i, values[i])

        valuesVtk.SetName(item)

        cellData.AddArray(valuesVtk)
        self.fields.append(item)

    def __delitem__(self, item):
        """Delete an internal field form the case.

        Parameters
        ----------
        item : str
            Name of the field to delete.

        """
        self.vtkData.VTKObject.GetCellData().RemoveArray(item)
        self.fields.remove(item)

    def extract_boundary_cells(self, boundary):
        """Extract cells adjacent to a certain boundary.

        Parameters
        ----------
        boundary : str
            The name of the boundary.

        Returns
        -------
            vtkExtractSelection

        """
        if boundary not in self.vtkData.FieldData.keys():
            raise(NameError("No boundary named "+boundary))

        cellIds = self.vtkData.FieldData[boundary]
        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(numpy_to_vtk(cellIds))

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        extractSelection = vtk.vtkExtractSelection()

        extractSelection.SetInputData(0, self.vtkData.VTKObject)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        return extractSelection

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
        selection = self.extract_boundary_cells(boundary)
        cCenters = vtk.vtkCellCenters()
        cCenters.SetInputData(selection.GetOutput())
        cCenters.Update()

        points = np.array(dsa.WrapDataObject(cCenters.GetOutput()).Points)
        dataVTK = dsa.WrapDataObject(selection.GetOutput()).CellData

        data = {}
        for key in dataVTK.keys():
            data[key] = np.array(dataVTK[key])

        if sort is None:
            return points, data
        elif sort == "x":
            ind = np.argsort(points[:, 0])
        elif sort == "y":
            ind = np.argsort(points[:, 1])

        points = points[ind]

        for key in data:
            data[key] = data[key][ind]

        return points, data

    def extract_block_by_name(self, name):
        """Extract a block from the case by a given name."""

        return self._blockData.GetBlock(self.boundaries.index(name) + 1)

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
            return points, data
        elif sort == "x":
            ind = np.argsort(points[:, 0])
        elif sort == "y":
            ind = np.argsort(points[:, 1])

        points = points[ind]

        for key in data:
            data[key] = data[key][ind]

        return points, data
