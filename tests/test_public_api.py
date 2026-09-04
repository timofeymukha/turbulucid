"""Guards on the shape of the installed package and its public names."""

import pytest

import turbulucid
from turbulucid import core

SUBMODULES = ["case", "plotting", "data_extraction", "quantities", "readers"]


def test_core_all_matches_the_submodules():
    """core.__all__ is a literal, so it can drift from the submodules."""
    expected = SUBMODULES + [
        name for module in SUBMODULES for name in getattr(core, module).__all__
    ]

    assert core.__all__ == expected


def test_top_level_all_re_exports_core():
    assert turbulucid.__all__ == ["core", *core.__all__]


@pytest.mark.parametrize("name", turbulucid.__all__)
def test_every_exported_name_resolves(name):
    assert hasattr(turbulucid, name)


def test_no_duplicate_exports():
    assert len(turbulucid.__all__) == len(set(turbulucid.__all__))


def test_private_helpers_are_not_exported():
    """Internal validators must not leak into the public namespace."""
    leaked = [name for name in turbulucid.__all__ if name.startswith("_")]

    assert leaked == []


def test_package_does_not_ship_docs_or_tests():
    """Namespace discovery used to install "docs" and "tests" top level."""
    from pathlib import Path

    # setuptools is not installed in a bare virtual environment.
    discovery = pytest.importorskip("setuptools.discovery")

    root = Path(turbulucid.__file__).parent.parent
    if not (root / "pyproject.toml").is_file():
        pytest.skip("not running from a source checkout")

    found = discovery.PEP420PackageFinder.find(
        where=str(root), include=["turbulucid*"])

    assert "docs" not in found
    assert "tests" not in found
    assert "turbulucid" in found
