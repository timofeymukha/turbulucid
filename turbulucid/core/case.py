from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

__all__ = ["Case"]


class Case:
    """A class for representing the simulation.

    Attributes
    ----------
    """

    def __init__(self, fileName):
        self.fileName = fileName

        if not os.path.exists(fileName):
            raise ValueError("ERROR: The file "+fileName+" does not exist")

        # Read the file
        self._reader = vtk.vtkPolyDataReader()
        self._reader.SetFileName(fileName)
        self._reader.Update()

        # Compute the cell-centres
        self._cellCentres = vtk.vtkCellCenters()
        self._cellCentres.SetInputConnection(self._reader.GetOutputPort())
        self._cellCentres.Update()
        self._cellCentres =\
            dsa.WrapDataObject(self._cellCentres.GetOutput()).GetPoints()

        self._vtkData = dsa.WrapDataObject(self._reader.GetOutput())

        self._zValue = self._vtkData.Points[0, 2]

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
        """vtkPolyDataReader : A vtk reader for the stored data."""

        return self._zValue

    def __getitem__(self, item):
        """Return a cell array of a given name.

        Parameters
        ----------
        item : string
            The name of the cell array.

        """
        return self.vtkData.CellData[item]

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
        extractSelection.SetInputConnection(0, self.reader.GetOutputPort())

        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        return extractSelection

    def boundary_cell_data(self, boundary):
        """Return cell-center coordinates and data from cells adjacent
        to a specific boundary.

        Parameters
        ----------
            boundary : str
                The name of the boundary.

        Returns
        -------


        """
        selection = self.extract_boundary_cells(boundary)
        cCenters = vtk.vtkCellCenters()
        cCenters.SetInputData(selection.GetOutput())
        cCenters.Update()

        coords = dsa.WrapDataObject(cCenters.GetOutput()).Points
        data = dsa.WrapDataObject(selection.GetOutput()).CellData

        return coords, data


