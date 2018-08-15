# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from __future__ import print_function
from __future__ import division
import pytest
from os import path
import vtk
import turbulucid
from turbulucid.bin.averageAlongAxis import *
from numpy.testing import assert_allclose


@pytest.fixture
def read_test_case_1():
    """Return the reader to the first test case."""

    casePath = path.join(turbulucid.__path__[0], "datasets",
                         "test_case_1", "test_case_1.foam")
    return read(casePath)


@pytest.fixture
def read_test_case_block():
    """Return the reader to the test case ."""

    casePath = path.join(turbulucid.__path__[0], "datasets",
                         "test_case_block", "test_case_block.foam")
    return read(casePath)


def test_get_block_names(read_test_case_1):
    blocks = read_test_case_1.GetOutput()

    names = get_block_names(blocks)
    patchNames = get_block_names(blocks.GetBlock(1))

    assert names == ['internalMesh', 'Patches']
    assert patchNames == ['left', 'right', 'top', 'inlet', 'botOrthoHex',
                          'botCurved', 'botPrism', 'outlet', 'botSkewedHex']


def test_get_block_index(read_test_case_1):
    blocks = read_test_case_1.GetOutput()
    patchBlocks = blocks.GetBlock(1)

    assert get_block_index(blocks, 'internalMesh') == 0
    assert get_block_index(blocks, 'Patches') == 1
    assert get_block_index(patchBlocks, 'left') == 0
    assert get_block_index(patchBlocks, 'outlet') == 7


def test_config_to_dict():
    configPath = path.join(turbulucid.__path__[0], "datasets",
                           "test_case_1", "averagingConfig")
    dict = config_to_dict(configPath)
    assert dict["case"] == "testpath"
    assert dict["patch"] == "left"
    assert dict["time"] == "0"
    assert dict["nSamples"] == "10"


# def test_average_internal_field_data(read_test_case_block):
#    blocks = read_test_case_block.GetOutput()
#    patchBlocks = blocks.GetBlock(1)

#    seedPatchName = "left"
#    seedPatchBlock = patchBlocks.GetBlock(get_block_index(patchBlocks,
#                                                          seedPatchName))

#    internalData = vtk.vtkPolyData()
#    internalData.ShallowCopy(seedPatchBlock)

#    average_internal_field_data(blocks.GetBlock(0), internalData, 10)

#    ccFilter = vtk.vtkCellCenters()
#    ccFilter.SetInputData(internalData)
#    ccFilter.Update()
#    cc = dsa.WrapDataObject(ccFilter.GetOutput()).Points

#    wrappedData = dsa.WrapDataObject(internalData).CellData
#    for i in range(cc.shape[0]):
#        assert_allclose(wrappedData["scalarField"][i], 2.5)
#        assert_allclose(wrappedData["vectorField"][i, 0], 2.5)
#        assert_allclose(wrappedData["vectorField"][i, 1], 2.5)
#        assert_allclose(wrappedData["vectorField"][i, 2], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 0], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 1], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 2], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 3], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 4], 2.5)
#        assert_allclose(wrappedData["symtensorField"][i, 5], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 0, 0], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 0, 1], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 0, 2], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 1, 0], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 1, 1], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 1, 2], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 2, 0], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 2, 1], 2.5)
#        assert_allclose(wrappedData["tensorField"][i, 2, 2], 2.5)


def test_create_boundary_polydata(read_test_case_block):
    blocks = read_test_case_block.GetOutput()
    patchBlocks = blocks.GetBlock(1)

    seedPatchName = "left"
    seedPatchBlock = patchBlocks.GetBlock(get_block_index(patchBlocks,
                                                          seedPatchName))

    internalData = vtk.vtkPolyData()
    internalData.ShallowCopy(seedPatchBlock)

    bounds = [0, 1, 0, 1, 0, 1]

    boundaryData = create_boundary_polydata(patchBlocks, internalData, bounds)

    correctKeys = ["inlet", "outlet", "bottomWall", "topWall"]
    correctPoints = {}
    correctPoints["inlet"] = [[0, 0, 0], [0, 0.5, 0], [0, 1, 0]]
    correctPoints["outlet"] = [[1, 0, 0], [1, 0.5, 0], [1, 1, 0]]
    correctPoints["bottomWall"] = [[0, 0, 0], [1/3, 0, 0], [2/3, 0., 0],
                                   [1, 0, 0]]

    correctPoints["topWall"] = [[0, 1, 0], [1/3, 1, 0], [2/3, 1., 0],
                                [1, 1, 0]]

    for boundary in boundaryData:
        assert boundary in correctKeys

        dataI = dsa.WrapDataObject(boundaryData[boundary])
        assert dataI.Points.shape[0] == len(correctPoints[boundary])

        flag = 0
#        for i in range(dataI.Points.shape[0]):
#            for j in range(dataI.Points.shape[0]):
#                print(np.any(dataI.Points[i, :] - correctPoints[boundary][j]))
#                if not np.any(dataI.Points[i, :] -
#                              correctPoints[boundary][j]):
#                    print(dataI.Points[i, :], j)
#                    flag = 1
#            assert flag == 1
#            flag = 0
