# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import numpy as np
import pytest
from numpy.testing import assert_allclose
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray, vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkCompositeDataSet,
    vtkPolyData,
)
from vtkmodules.vtkCommonTransforms import vtkTransform

from turbulucid.core.readers import LegacyReader, Reader, VTUReader, XMLReader


def test_reader_base_class_is_abstract():
    with pytest.raises(TypeError, match="abstract"):
        Reader("unused")


def create_single_cell(z, axis, angle):
    from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

    points = vtkPoints()
    points.InsertPoint(0, 0.0, 0.0, z)
    points.InsertPoint(1, 1.0, 0.0, z)
    points.InsertPoint(2, 1.0, 1.0, z)
    points.InsertPoint(3, 0.0, 1.0, z)

    strips = vtkCellArray()
    strips.InsertNextCell(4)
    strips.InsertCellPoint(0)
    strips.InsertCellPoint(1)
    strips.InsertCellPoint(2)
    strips.InsertCellPoint(3)

    data = vtkPolyData()
    data.SetPoints(points)
    data.SetPolys(strips)

    v = vtkDoubleArray()
    v.SetName("Pressure")
    v.InsertNextValue(2.7)

    data.GetCellData().SetScalars(v)

    transform = vtkTransform()

    transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])
    transform.Update()

    filter = vtkTransformPolyDataFilter()
    filter.SetInputData(data)
    filter.SetTransform(transform)
    filter.Update()

    return filter.GetOutput()


def create_polydata(point_values, cells):
    points = vtkPoints()
    for index, point in enumerate(point_values):
        points.InsertPoint(index, *point)

    polygons = vtkCellArray()
    for cell in cells:
        polygons.InsertNextCell(len(cell))
        for point_id in cell:
            polygons.InsertCellPoint(point_id)

    data = vtkPolyData()
    data.SetPoints(points)
    data.SetPolys(polygons)

    values = vtkDoubleArray()
    values.SetName("Pressure")
    for index in range(len(cells)):
        values.InsertNextValue(float(index))
    data.GetCellData().SetScalars(values)
    return data


def write_data(data, writerType, path):
    """Write data to a temporary directory and return the path to the file.

    """
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    from vtkmodules.vtkIOLegacy import vtkPolyDataWriter

    if writerType == "legacy":
        writer = vtkPolyDataWriter()
        format = "vtk"
    else:
        writer = vtkXMLPolyDataWriter()
        format = "vtp"

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
    from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
    data = fixture

    filename = write_data(data, writer, tmpdir)

    if writer == "legacy":
        reader = LegacyReader(filename)
    else:
        reader = XMLReader(filename)
    readerData = reader.data

    vtkNormals = vtkPolyDataNormals()
    vtkNormals.ComputeCellNormalsOn()
    vtkNormals.SetInputData(readerData.GetBlock(0))
    vtkNormals.Update()

    normals = dsa.WrapDataObject(vtkNormals.GetOutput()).CellData["Normals"]
    meanNormal = np.mean(normals, axis=0)
    meanNormal /= np.linalg.norm(meanNormal)
    assert_allclose(meanNormal[:2], [0, 0], rtol=1e-5, atol=1e-5)

    # Can be both +1 and -1, but aligned with z
    assert_allclose([1], np.abs(meanNormal[-1]), rtol=1e-5, atol=1e-5)


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
    assert(readerData.GetMetaData(0).Get(vtkCompositeDataSet.NAME()) ==
           "internalField")
    assert(readerData.GetMetaData(1).Get(vtkCompositeDataSet.NAME()) ==
           "boundary")


def test_xml_block_structure(tmpdir):
    data = create_single_cell(0, [0, 1, 0], 0)
    filename = write_data(data, "xml", tmpdir)
    reader = XMLReader(filename)
    readerData = reader.data

    assert(readerData.GetNumberOfBlocks() == 2)
    assert(readerData.GetMetaData(0).Get(vtkCompositeDataSet.NAME()) ==
           "internalField")
    assert(readerData.GetMetaData(1).Get(vtkCompositeDataSet.NAME()) ==
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


def test_xml_unstructured_grid_reader(tmp_path):
    from vtkmodules.vtkFiltersCore import vtkAppendFilter
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter

    polydata = create_single_cell(0, [0, 1, 0], 0)
    converter = vtkAppendFilter()
    converter.AddInputData(polydata)
    converter.Update()

    filename = tmp_path / "test.vtu"
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(filename))
    writer.SetInputData(converter.GetOutput())
    assert writer.Write() == 1

    reader = XMLReader(str(filename))
    assert reader.data.GetBlock(0).GetNumberOfCells() == 1
    assert reader.data.GetBlock(0).GetNumberOfPoints() == 4

    with pytest.warns(DeprecationWarning, match="use XMLReader"):
        compatibility_reader = VTUReader(str(filename))
    assert compatibility_reader.data.GetBlock(0).GetNumberOfCells() == 1
    assert compatibility_reader.data.GetBlock(0).GetNumberOfPoints() == 4


def test_xml_reader_rejects_unsupported_extension(tmp_path):
    filename = tmp_path / "test.xml"
    filename.touch()

    with pytest.raises(ValueError, match="Unsupported XML VTK file extension"):
        XMLReader(str(filename))


def test_plane_fit_is_independent_of_cell_winding(tmpdir):
    from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

    data = create_polydata(
        [
            (0, 0, 0),
            (1, 0, 0),
            (2, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            (2, 1, 0),
        ],
        [
            (0, 1, 4, 3),
            (1, 4, 5, 2),  # Deliberately opposite winding.
        ],
    )
    transform = vtkTransform()
    transform.RotateWXYZ(57, 1, 0.3, 0.2)
    transformed = vtkTransformPolyDataFilter()
    transformed.SetInputData(data)
    transformed.SetTransform(transform)
    transformed.Update()

    filename = write_data(transformed.GetOutput(), "xml", tmpdir)
    reader = XMLReader(filename)
    points = dsa.WrapDataObject(reader.data.GetBlock(0)).Points

    assert_allclose(points[:, 2], 0, atol=1e-7)


def test_plane_fit_rejects_nonplanar_geometry(tmpdir):
    data = create_polydata(
        [(0, 0, 0), (1, 0, 0), (1, 1, 0.1), (0, 1, 0)],
        [(0, 1, 2, 3)],
    )
    filename = write_data(data, "xml", tmpdir)

    with pytest.raises(ValueError, match="not planar"):
        XMLReader(filename)


def test_plane_fit_rejects_collinear_geometry(tmpdir):
    data = create_polydata(
        [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)],
        [(0, 1, 2, 3)],
    )
    filename = write_data(data, "xml", tmpdir)

    with pytest.raises(ValueError, match="degenerate or collinear"):
        XMLReader(filename)


def cell_areas(data):
    """Return the area of every cell in a polydata."""
    from vtkmodules.vtkFiltersVerdict import vtkMeshQuality

    quality = vtkMeshQuality()
    quality.SetTriangleQualityMeasureToArea()
    quality.SetQuadQualityMeasureToArea()
    quality.SetInputData(data)
    quality.Update()
    return np.asarray(dsa.WrapDataObject(quality.GetOutput()).CellData["Quality"])


def degenerate_quad_mesh():
    """Two unit quads with a zero-area quad wedged between them."""
    return create_polydata(
        [
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            # Two coincident pairs, so this quad encloses no area at all.
            (2, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 0),
            (3, 0, 0), (4, 0, 0), (4, 1, 0), (3, 1, 0),
        ],
        [
            (0, 1, 2, 3),
            (4, 5, 6, 7),
            (8, 9, 10, 11),
        ],
    )


@pytest.mark.parametrize("writer", ["legacy", "xml"])
def test_clean_removes_degenerate_quads(writer, tmpdir):
    """vtkMeshQuality scores quads with its own measure, which must be area.

    With the default quad measure a degenerate cell scores 1e30 rather than
    0, so it never falls below the threshold and clean silently does
    nothing on the quad meshes that make up most cut planes.

    """
    filename = write_data(degenerate_quad_mesh(), writer, tmpdir)
    readerType = LegacyReader if writer == "legacy" else XMLReader

    uncleaned = readerType(filename, clean=False)
    assert uncleaned.data.GetBlock(0).GetNumberOfCells() == 3

    cleaned = readerType(filename, clean=True).data.GetBlock(0)
    assert cleaned.GetNumberOfCells() == 2

    assert_allclose(sorted(cell_areas(cleaned)), [1.0, 1.0])


@pytest.mark.parametrize("writer", ["legacy", "xml"])
def test_clean_keeps_a_healthy_quad_mesh_intact(writer, tmpdir):
    data = create_polydata(
        [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, 0), (1, 1, 0), (2, 1, 0)],
        [(0, 1, 4, 3), (1, 2, 5, 4)],
    )
    filename = write_data(data, writer, tmpdir)
    readerType = LegacyReader if writer == "legacy" else XMLReader

    assert readerType(filename, clean=True).data.GetBlock(0).GetNumberOfCells() == 2
