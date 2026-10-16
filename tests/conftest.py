import importlib.resources
import json
from pathlib import Path

import pytest

from diffpy.structure import PDFFitStructure


@pytest.fixture
def user_filesystem(tmp_path):
    base_dir = Path(tmp_path)
    home_dir = base_dir / "home_dir"
    home_dir.mkdir(parents=True, exist_ok=True)
    cwd_dir = base_dir / "cwd_dir"
    cwd_dir.mkdir(parents=True, exist_ok=True)

    home_config_data = {"username": "home_username", "email": "home@email.com"}
    with open(home_dir / "diffpyconfig.json", "w") as f:
        json.dump(home_config_data, f)
    yield tmp_path


@pytest.fixture(scope="session")
def datafile():
    """Fixture to load a test data file from the testdata package
    directory."""

    def _datafile(filename):
        return importlib.resources.files("tests.testdata").joinpath(filename)

    return _datafile


@pytest.fixture(scope="session")
def as_list():
    def _as_list(values):
        """Unavailable uncertainties are None rather than an array."""
        return None if values is None else values.tolist()

    return _as_list


@pytest.fixture
def ni_pdffit_structure(datafile):
    """Fixture for the Ni structure from ni.cif with Biso = 1 on every
    atom."""
    structure = PDFFitStructure()
    structure.read(str(datafile("ni.cif")))
    for atom in structure:
        atom.Bisoequiv = 1
    return structure
