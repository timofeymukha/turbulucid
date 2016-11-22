from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
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
            print("ERROR: The file "+fileName+" does not exist")
            raise ValueError

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

    def __getitem__(self, item):
        """Return a cell array of a given name.

        Parameters
        ----------
        item : string
            The name of the cell array.

        """
        return self.vtkData.CellData[item]

