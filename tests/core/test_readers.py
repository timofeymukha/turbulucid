# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from __future__ import print_function
from __future__ import division
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
import turbulucid
from turbulucid.core.readers import *
import numpy as np
from numpy.testing import assert_allclose
import pytest


def create_single_cell(z, axis, angle):

    points = vtk.vtkPoints()
    points.InsertPoint(0, 0.0, 0.0, z)
    points.InsertPoint(1, 1.0, 0.0, z)
    points.InsertPoint(2, 1.0, 1.0, z)
    points.InsertPoint(3, 0.0, 1.0, z)

    strips = vtk.vtkCellArray()
    strips.InsertNextCell(4)
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

    transform = vtk.vtkTransform()

    transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
    transform.Update()

    filter = vtk.vtkTransformPolyDataFilter()
    filter.SetInputData(data)
    filter.SetTransform(transform)
    filter.Update()

    return filter.GetOutput()


def write_data(data, writerType, path):
    """Write data to a temporary directory and return the path to the file.

    """
    if writerType == "legacy":
        writer = vtk.vtkPolyDataWriter()
        format = "vtk"
    else:
        writer = vtk.vtkXMLPolyDataWriter()
        format = "vtu"

    filename = path.join("test." + format).strpath
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()

    return writer.GetFileName()


# Fixtures for testing different types of initial data
@pytest.fixture(params=[-1, 0, 1])
def different_z(request):
    return create_single_cell(request.param, [0, 1, 0], 0)


@pytest.fixture(params=[0, 45, 90, 180, 270, 360])
def different_angle(request):
    return create_single_cell(0, [0, 1, 0], request.param)


@pytest.fixture(params=[[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1],
                        [0, 1, 1], [0.5, 0.345, 0.12]])
def different_axis(request):
    return create_single_cell(0, request.param, 56.562)


# Test for the direction of the normal being [0, 0, 1]
def normal_direction(fixture, writer, tmpdir):
    data = fixture

    filename = write_data(data, writer, tmpdir)

    if writer == "legacy":
        reader = LegacyReader(filename)
    else:
        reader = XMLReader(filename)
    readerData = reader.data

    vtkNormals = vtk.vtkPolyDataNormals()
    vtkNormals.ComputeCellNormalsOn()
    vtkNormals.SetInputData(readerData.GetBlock(0))
    vtkNormals.Update()

    normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
    meanNormal = np.mean(normals, axis=0)
    meanNormal /= np.linalg.norm(meanNormal)
    assert_allclose(meanNormal, [0, 0, 1], rtol=1e-5, atol=1e-5)


def test_legacy_normal_direction_different_z(different_z, tmpdir):
    normal_direction(different_z, "legacy", tmpdir)


def test_legacy_normal_direction_different_angle(different_angle, tmpdir):
    normal_direction(different_angle, "legacy", tmpdir)


def test_legacy_normal_direction_different_axis(different_axis, tmpdir):
    normal_direction(different_axis, "legacy", tmpdir)


def test_xml_normal_direction_different_z(different_z, tmpdir):
    normal_direction(different_z, "xml", tmpdir)


def test_xml_normal_direction_different_angle(different_angle, tmpdir):
    normal_direction(different_angle, "xml", tmpdir)


def test_xml_normal_direction_different_axis(different_axis, tmpdir):
    normal_direction(different_axis, "xml", tmpdir)


# Test for z value of the resulting data being 0
def zvalue(fixture, writer, tmpdir):
    data = fixture

    filename = write_data(data, writer, tmpdir)

    if writer == "legacy":
        reader = LegacyReader(filename)
    else:
        reader = XMLReader(filename)

    readerData = reader.data

    for pointI in range(readerData.GetBlock(0).GetNumberOfPoints()):
        zValue = readerData.GetBlock(0).GetPoint(pointI)[-1]
        assert_allclose(zValue, [0], rtol=1e-5, atol=1e-5)


def test_legacy_zvalue_different_z(different_z, tmpdir):
    zvalue(different_z, "legacy", tmpdir)


def test_legacy_zvalue_different_angle(different_angle, tmpdir):
    zvalue(different_angle, "legacy", tmpdir)


def test_legacy_zvalue_different_axis(different_axis, tmpdir):
    zvalue(different_axis, "legacy", tmpdir)


def test_xml_zvalue_different_z(different_z, tmpdir):
    zvalue(different_z, "xml", tmpdir)


def test_xml_zvalue_different_angle(different_angle, tmpdir):
    zvalue(different_angle, "xml", tmpdir)


def test_xml_zvalue_different_axis(different_axis, tmpdir):
    zvalue(different_axis, "xml", tmpdir)


# Test multiblock structure
def test_legacy_block_structure(tmpdir):
    data = create_single_cell(0, [0, 1, 0], 0)
    filename = write_data(data, "legacy", tmpdir)
    reader = LegacyReader(filename)
    readerData = reader.data

    assert(readerData.GetNumberOfBlocks() == 2)
    assert(readerData.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "internalField")
    assert(readerData.GetMetaData(1).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "boundary")


def test_xml_block_structure(tmpdir):
    data = create_single_cell(0, [0, 1, 0], 0)
    filename = write_data(data, "xml", tmpdir)
    reader = XMLReader(filename)
    readerData = reader.data

    assert(readerData.GetNumberOfBlocks() == 2)
    assert(readerData.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "internalField")
    assert(readerData.GetMetaData(1).Get(vtk.vtkCompositeDataSet.NAME()) ==
           "boundary")


# Test boundary field data
def test_legacy_boundary_field_data(tmpdir):
    data = create_single_cell(0, [0, 1, 0], 0)
    filename = write_data(data, "legacy", tmpdir)
    reader = LegacyReader(filename)
    readerData = reader.data
    w = dsa.WrapDataObject(readerData.GetBlock(0))

    assert("boundaries" in w.FieldData.keys())
    assert(w.FieldData["boundaries"].GetNumberOfTuples() == 1)


def test_xml_boundary_field_data(tmpdir):
    data = create_single_cell(0, [0, 1, 0], 0)
    filename = write_data(data, "xml", tmpdir)
    reader = XMLReader(filename)
    readerData = reader.data
    w = dsa.WrapDataObject(readerData.GetBlock(0))

    assert("boundaries" in w.FieldData.keys())
    assert(w.FieldData["boundaries"].GetNumberOfTuples() == 1)
