from __future__ import print_function
from __future__ import division
import pytest
import vtk
import turbulucid
from turbulucid.core.readers import *
from numpy.testing import assert_allclose


#def test_create_abstract_reader():
#    with pytest.raises(NotImplementedError):
#        Reader()

def create_single_cell():
    points = vtk.vtkPoints()
    points.InsertPoint(0, 0.0, 0.0, 0.0)
    points.InsertPoint(1, 1.0, 0.0, 0.0)
    points.InsertPoint(2, 1.0, 1.0, 0.0)
    points.InsertPoint(3, 0.0, 1.0, 0.0)

    strips = vtk.vtkCellArray()
    strips.InsertNextCell(4) # number of points
    strips.InsertCellPoint(0)
    strips.InsertCellPoint(1)
    strips.InsertCellPoint(2)
    strips.InsertCellPoint(3)

    data = vtk.vtkPolyData()
    data.SetPoints(points)
    data.SetPolys(strips)

    v = vtk.vtkDoubleArray()
    v.SetName("Pressure")
    v.InsertNextValue(2.7)

    data.GetCellData().SetScalars(v)
    return data


def test_legacy_reader(tmpdir):
    data = create_single_cell()

    writer = vtk.vtkPolyDataWriter()
    print(type(tmpdir))
    writer.SetFileName(tmpdir.join("test.vtk").strpath)
    writer.SetInputData(data)
    writer.Write()

    filename = writer.GetFileName()
    reader = LegacyReader(filename)
    readerData = reader.data

    assert(readerData.GetNumberOfBlocks() == 2)
    assert(readerData.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "internalField")
    assert(readerData.GetMetaData(1).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "boundary")

    vtkNormals = vtk.vtkPolyDataNormals()
    vtkNormals.ComputeCellNormalsOn()
    vtkNormals.SetInputData(readerData.GetBlock(0))
    vtkNormals.Update()
    normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
    meanNormal = np.mean(normals, axis=0)
    meanNormal /= np.linalg.norm(meanNormal)



