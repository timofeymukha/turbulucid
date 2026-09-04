from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import turbulucid


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


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Keep plotting tests isolated from Matplotlib's global state."""
    yield
    plt.close("all")
