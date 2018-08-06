# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from .readers import NativeReader, LegacyReader, XMLReader

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
        self.fileName = fileName

        # Read in the data
        self._blockData = self.read(clean, pointData)

        # Compute the cell-centres
        self._cellCentres = vtk.vtkCellCenters()
        self._cellCentres.SetInputData(self._blockData.GetBlock(0))
        self._cellCentres.Update()
        self._cellCentres =\
            dsa.WrapDataObject(self._cellCentres.GetOutput()).GetPoints()
        self._cellCentres = np.array(self._cellCentres[:, :2])

        self._vtkData = dsa.WrapDataObject(self._blockData.GetBlock(0))

        self._boundaries = self._fill_boundary_list()

        self._bounds = self._vtkData.VTKObject.GetBounds()[:4]

        self._fields = self._vtkData.CellData.keys()

        plot_limits = self._compute_plot_limits()
        self._xlim = plot_limits[0]
        self._ylim = plot_limits[1]

    @property
    def blockData(self):
        """vtkMultiBlockDataSet : the multiblock data assembled by the
        reader.

        """
        return self._blockData

    @property
    def vtkData(self):
        """wrapped PolyData : The actual data read by the reader."""

        return self._vtkData

    @property
    def cellCentres(self):
        """wrapped VTKArray : the cell centres of the read data """

        return self._cellCentres

    @property
    def boundaries(self):
        """list : A list of names of the boundaries present the case."""

        return self._boundaries

    @property
    def bounds(self):
        """tuple : (min(x), max(x), min(y), max(y))."""

        return self._bounds

    @property
    def fields(self):
        """list of str: The names of the fields present in the case."""

        return self._fields

    @property
    def xlim(self):
        """list of two floats: The x limits that cover the
        geometry of the case, plus small a margin.

        """
        return self._xlim

    @property
    def ylim(self):
        """list of two floats: The y limits that cover the
        geometry of the case, plus a small margin.

        """
        return self._ylim

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
        if item not in self._fields:
            raise ValueError("Field " + item + " not present in the case.")

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

        self.fields.append(item)

        cellData = self._vtkData.VTKObject.GetCellData()
        valuesVtk = vtk.vtkDoubleArray()

        if np.ndim(values) > 1:
            valuesVtk.SetNumberOfComponents(values.shape[1])
            valuesVtk.SetNumberOfTuples(values.shape[0])
            for i in range(values.shape[0]):
                valuesVtk.SetTuple(i, values[i, :])
        else:
            valuesVtk.SetNumberOfComponents(1)
            valuesVtk.SetNumberOfValues(values.shape[0])
            for i in range(values.shape[0]):
                valuesVtk.SetValue(i, values[i])

        valuesVtk.SetName(item)

        cellData.AddArray(valuesVtk)

        # Add boundary data by copying from boundary cells
        for boundary in self.boundaries:
            block = self.extract_block_by_name(boundary)
            cellData = block.GetCellData()
            valuesVtk = vtk.vtkDoubleArray()

            nVals = self.boundary_cell_data(boundary)[0][:, 0].size
            bCellData = self.boundary_cell_data(boundary)[1][item]

            if np.ndim(values) > 1:
                valuesVtk.SetNumberOfComponents(values.shape[1])
                valuesVtk.SetNumberOfTuples(nVals)
                for i in range(nVals):
                    valuesVtk.SetTuple(i, bCellData[i, :])
            else:
                valuesVtk.SetNumberOfComponents(1)
                valuesVtk.SetNumberOfValues(nVals)
                for i in range(nVals):
                    valuesVtk.SetValue(i, bCellData[i])

            valuesVtk.SetName(item)
            cellData.AddArray(valuesVtk)

    def __delitem__(self, item):
        """Delete an internal field form the case.

        Parameters
        ----------
        item : str
            Name of the field to delete.

        """
        self.vtkData.VTKObject.GetCellData().RemoveArray(item)
        self.fields.remove(item)

        for boundary in self.boundaries:
            block = self.extract_block_by_name(boundary)
            block.GetCellData().RemoveArray(item)

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

        # Transform the internal field
        filter = vtk.vtkTransformPolyDataFilter()
        filter.SetInputData(self.blockData.GetBlock(0))
        filter.SetTransform(transform)
        filter.Update()

        self._blockData.SetBlock(0, filter.GetOutput())

        # Transform boundary data
        i = 1
        for boundary in self.boundaries:
            filter = vtk.vtkTransformPolyDataFilter()
            filter.SetTransform(transform)
            filter.SetInputData(self.blockData.GetBlock(i))
            filter.Update()
            self.blockData.SetBlock(i, filter.GetOutput())
            i += 1

        # Update attuributes
        self._cellCentres = vtk.vtkCellCenters()
        self._cellCentres.SetInputData(self.blockData.GetBlock(0))
        self._cellCentres.Update()
        self._cellCentres = \
            dsa.WrapDataObject(self._cellCentres.GetOutput()).GetPoints()
        self._cellCentres = np.array(self._cellCentres[:, :2])

        self._vtkData = dsa.WrapDataObject(self._blockData.GetBlock(0))

        self._bounds = self._vtkData.VTKObject.GetBounds()[:4]

        plot_limits = self._compute_plot_limits()
        self._xlim = plot_limits[0]
        self._ylim = plot_limits[1]

    def translate(self, dx, dy):
        """Translate the geometry of the case.

        Parameters
        ----------
        dx : float
            The translation along the x axis.
        dy : float
            The translation along the y axis.

        """
        transform = vtk.vtkTransform()
        transform.Translate(dx, dy, 0)
        transform.Update()

        self._transform(transform)

    def rotate(self, angle):
        """Rotate the geometry of the case around the z axis.

        Parameters
        ----------
        dx : angle
            Rotation angle in degrees.

        """
        axis = [0, 0, 1]
        transform = vtk.vtkTransform()
        transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
        transform.Update()
        self._transform(transform)

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
            return points[:, [0, 1]], data
        elif sort == "x":
            ind = np.argsort(points[:, 0])
        elif sort == "y":
            ind = np.argsort(points[:, 1])

        points = points[ind]

        for key in data:
            data[key] = data[key][ind]

        return points[:, [0, 1]], data

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
            return points[:, [0, 1]], data
        elif sort == "x":
            ind = np.argsort(points[:, 0])
        elif sort == "y":
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

        fileExt = os.path.splitext(fileName)[1]

        if fileExt == ".vtm":
            return NativeReader(fileName).data
        elif fileExt == ".vtk":
            return LegacyReader(fileName, clean=clean,
                                pointData=pointData).data
        elif (fileExt == ".vtu") or (fileExt == ".vtp"):
            return XMLReader(fileName, clean=clean, pointData=pointData).data
        else:
            raise ValueError("Unsupported file format.", fileName, fileExt)

    def write(self, writePath):
        """Save the case to a .vtm format.

        Parameters
        ----------
        writePath : str
            The name of the file.

        """
        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetFileName(writePath)
        writer.SetInputData(self._blockData)
        writer.Write()
