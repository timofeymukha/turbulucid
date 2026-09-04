# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from os import path

import numpy as np
import pytest
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkPolyData

import turbulucid
from turbulucid.bin.averageAlongAxis import (
    average_internal_field_data,
    config_to_dict,
    create_boundary_polydata,
    get_block_index,
    get_block_names,
    new_average_internal_field_data,
    read,
    zero_out_arrays,
)


def dataset_path(*parts):
    return path.join(turbulucid.__path__[0], "datasets", *parts)


@pytest.fixture
def read_test_case_1():
    """Return the reader to the first test case."""
    return read(dataset_path("test_case_1", "test_case_1.foam"), 0.0)


@pytest.fixture
def read_test_case_block():
    """Return the reader to the test case ."""
    return read(dataset_path("test_case_block", "test_case_block.foam"), 0.0)


def seed_patch_copy(reader, name="left"):
    """A standalone copy of one patch, ready to be averaged onto."""
    patchBlocks = reader.GetOutput().GetBlock(1)
    block = patchBlocks.GetBlock(get_block_index(patchBlocks, name))

    internalData = vtkPolyData()
    internalData.ShallowCopy(block)
    internalData.BuildLinks()
    return internalData


def test_get_block_names(read_test_case_1):
    blocks = read_test_case_1.GetOutput()

    names = get_block_names(blocks)
    patchNames = get_block_names(blocks.GetBlock(1))

    # VTK renamed this top-level block from "Patches" to "boundary".
    # The structure and patch metadata are unchanged, so avoid coupling the
    # helper test to a particular VTK release.
    assert names[0] == 'internalMesh'
    assert names[1] in {'Patches', 'boundary'}
    assert patchNames == ['left', 'right', 'top', 'inlet', 'botOrthoHex',
                          'botCurved', 'botPrism', 'outlet', 'botSkewedHex']


def test_get_block_index(read_test_case_1):
    blocks = read_test_case_1.GetOutput()
    patchBlocks = blocks.GetBlock(1)
    boundaryBlockName = get_block_names(blocks)[1]

    assert get_block_index(blocks, 'internalMesh') == 0
    assert get_block_index(blocks, boundaryBlockName) == 1
    assert get_block_index(patchBlocks, 'left') == 0
    assert get_block_index(patchBlocks, 'outlet') == 7


def test_get_block_index_rejects_unknown_name(read_test_case_1):
    with pytest.raises(NameError, match="No block named"):
        get_block_index(read_test_case_1.GetOutput(), "notABlock")


def test_config_to_dict():
    config = config_to_dict(dataset_path("test_case_1", "averagingConfig"))

    assert config["case"] == "testpath"
    assert config["patch"] == "left"
    assert config["time"] == "0"
    assert config["nSamples"] == "10"


@pytest.mark.parametrize(
    "averager",
    [
        # The slow variant reports progress, and print_progress divides by
        # int(total/freq), which is zero for this six-cell patch. That is a
        # known defect in turbulucid/bin, left for a separate pass.
        pytest.param(
            average_internal_field_data,
            marks=pytest.mark.filterwarnings("ignore:divide by zero"),
            id="line-probe",
        ),
        pytest.param(new_average_internal_field_data, id="bulk-probe"),
    ],
)
def test_average_internal_field_data(read_test_case_block, averager):
    """The mesh has four cells along z holding 1, 2, 3 and 4, mean 2.5.

    nSamples is a multiple of that cell count so that the point sampling
    resolves the cells exactly and both algorithms must agree. With
    nSamples=10 the bulk-probe variant returns 2.4 instead, because ten
    equispaced points cover four cells as 3/2/3/2. See the notes on
    turbulucid/bin.

    """
    blocks = read_test_case_block.GetOutput()
    internalData = seed_patch_copy(read_test_case_block)
    zero_out_arrays(internalData)

    averager(blocks.GetBlock(0), internalData, 40, False, False)

    cellData = dsa.WrapDataObject(internalData).CellData
    shapes = {
        "scalarField": (6,),
        "vectorField": (6, 3),
        "symtensorField": (6, 6),
        "tensorField": (6, 3, 3),
    }
    for field, shape in shapes.items():
        values = np.asarray(cellData[field])
        assert values.shape == shape
        np.testing.assert_allclose(values, 2.5, rtol=1e-6)


def test_create_boundary_polydata(read_test_case_block):
    blocks = read_test_case_block.GetOutput()
    patchBlocks = blocks.GetBlock(1)
    internalData = seed_patch_copy(read_test_case_block)

    boundaryData = create_boundary_polydata(
        patchBlocks, internalData, [0, 1, 0, 1, 0, 1], False)

    expected = {
        "inlet": [[0, 0, 0], [0, 0.5, 0], [0, 1, 0]],
        "outlet": [[1, 0, 0], [1, 0.5, 0], [1, 1, 0]],
        "bottomWall": [[0, 0, 0], [1/3, 0, 0], [2/3, 0, 0], [1, 0, 0]],
        "topWall": [[0, 1, 0], [1/3, 1, 0], [2/3, 1, 0], [1, 1, 0]],
    }
    assert sorted(boundaryData) == sorted(expected)

    for boundary, expectedPoints in expected.items():
        points = np.asarray(dsa.WrapDataObject(boundaryData[boundary]).Points)
        assert points.shape == (len(expectedPoints), 3)

        # The filter does not promise an order, so compare as sorted sets.
        order = np.lexsort((points[:, 1], points[:, 0]))
        np.testing.assert_allclose(
            points[order], sorted(expectedPoints), rtol=1e-5, atol=1e-6)
