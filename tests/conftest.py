from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import turbulucid


def build_polydata(pointValues, cells, fieldName="Pressure", values=None):
    """A vtkPolyData with the given points, polygons and one cell array."""
    from vtkmodules.vtkCommonCore import vtkDoubleArray, vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData

    points = vtkPoints()
    for index, point in enumerate(pointValues):
        points.InsertPoint(index, *point)

    polygons = vtkCellArray()
    for cell in cells:
        polygons.InsertNextCell(len(cell))
        for pointId in cell:
            polygons.InsertCellPoint(pointId)

    data = vtkPolyData()
    data.SetPoints(points)
    data.SetPolys(polygons)

    if values is None:
        values = range(len(cells))
    array = vtkDoubleArray()
    array.SetName(fieldName)
    for value in values:
        array.InsertNextValue(float(value))
    data.GetCellData().SetScalars(array)
    return data


def write_polydata(data, path):
    """Write a polydata to a .vtp file and return the path as a string."""
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter

    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(data)
    assert writer.Write() == 1
    return str(path)


def polygon_vertices_by_cell(data, scaleX=1, scaleY=1):
    """Reference implementation: walk the cells one at a time.

    The plotting code takes a vectorised route through the cell array;
    this is what it has to agree with.

    """
    polygons = []
    for i in range(data.GetNumberOfCells()):
        cell = data.GetCell(i)
        nPoints = cell.GetNumberOfPoints()
        vertices = np.zeros((nPoints, 2))
        for pointI in range(nPoints):
            vertices[pointI, :] = cell.GetPoints().GetPoint(pointI)[:2]
            vertices[pointI, :] /= [scaleX, scaleY]
        polygons.append(vertices)
    return polygons


@pytest.fixture
def polygon_reference():
    """The per-cell reference walk, for comparing against the fast path."""
    return polygon_vertices_by_cell


@pytest.fixture
def block_case():
    """Return a fresh, small native case for core API tests."""
    case_path = (
        Path(turbulucid.__path__[0])
        / "datasets"
        / "test_case_block"
        / "averaged.vtm"
    )
    return turbulucid.Case(str(case_path))


@pytest.fixture
def mixed_cell_case(tmp_path):
    """A case whose cells are a quad, a triangle and a pentagon.

    The bundled cases are uniform quad meshes, so only this one exercises
    the mixed-vertex-count branch of the polygon extraction.

    """
    data = build_polydata(
        [
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (2, 0, 0), (1, 1, 0),
            (3, 0, 0), (3.2, 0.6, 0), (2.5, 1, 0), (1.8, 0.6, 0),
        ],
        [
            (0, 1, 2, 3),
            (1, 4, 5),
            (4, 6, 7, 8, 9),
        ],
        fieldName="p",
    )
    return turbulucid.Case(write_polydata(data, tmp_path / "mixed.vtp"))


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Keep plotting tests isolated from Matplotlib's global state."""
    yield
    plt.close("all")
