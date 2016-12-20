from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
from collections import OrderedDict
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

__all__ = ["Case"]


class Case:
    """A class for representing the simulation case.

    Attributes
    ----------
    """

    def __init__(self, fileName):
        self.fileName = fileName

        if not os.path.exists(fileName):
            raise ValueError("ERROR: The file "+fileName+" does not exist")

        # Read the file
        self._reader = vtk.vtkXMLMultiBlockDataReader()

        self._reader.SetFileName(fileName)
        self._reader.Update()

        self._blockData = self._reader.GetOutput()

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


    @property
    def reader(self):
        """vtkPolyDataReader : A vtk reader for the stored data."""

        return self._reader

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
        """list : The fields present in the case."""

        return self._fields

    def _fill_boundary_list(self):
        fieldData = self.vtkData.FieldData['boundaries']
        boundaryList = []

        for i in range(fieldData.GetNumberOfValues()):
            boundaryList.append(fieldData.GetValue(i))

        return boundaryList

    def __getitem__(self, item):
        """Return a cell array of a given name.

        Parameters
        ----------
        item : string
            The name of the cell array.

        """
        # TODO Consider tensors
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

        # TODO Consider tensors
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

    def boundary_cell_data(self, boundary):
        """Return cell-centre coordinates and data from cells adjacent
        to a specific boundary.

        Parameters
        ----------
        boundary : str
            The name of the boundary.

        Returns
        -------
            Two ndarrays

        """
        selection = self.extract_boundary_cells(boundary)
        cCenters = vtk.vtkCellCenters()
        cCenters.SetInputData(selection.GetOutput())
        cCenters.Update()

        coords = dsa.WrapDataObject(cCenters.GetOutput()).Points
        data = dsa.WrapDataObject(selection.GetOutput()).CellData

        return coords, data

    def extract_block_by_name(self, name):
        """Extract a block from the case by a given name."""

        return self._blockData.GetBlock(self.boundaries.index(name) + 1)

    def boundary_data(self, boundary):
        """Return cell-center coordinates and data from a boundary.

        Parameters
        ----------
        boundary : str
            The name of the boundary.

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

        points = dsa.WrapDataObject(cCenters.GetOutput()).Points
        data = dsa.WrapDataObject(blockData).CellData

        return points, data
