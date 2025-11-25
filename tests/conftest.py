import importlib.resources
import logging
from functools import lru_cache

import pytest

logger = logging.getLogger(__name__)


# diffpy.structure
@lru_cache()
def has_diffpy_structure():
    try:
        import diffpy.structure as m

        del m
        return True
    except ImportError:
        logger.warning(
            "Cannot import diffpy.structure, Structure tests skipped."
        )
        return False


# diffpy.srreal


@lru_cache()
def has_diffpy_srreal():
    try:
        import diffpy.srreal.pdfcalculator as m

        del m
        return True
    except ImportError:
        logger.warning("Cannot import diffpy.srreal, PDF tests skipped.")
        return False


@pytest.fixture(scope="session")
def diffpy_structure_available():
    return has_diffpy_structure()


@pytest.fixture(scope="session")
def diffpy_srreal_available():
    return has_diffpy_srreal()


@pytest.fixture(scope="session")
def datafile():
    """Fixture to load a test data file from the testdata package
    directory."""

    def _datafile(filename):
        return importlib.resources.files("tests.testdata").joinpath(filename)

    return _datafile
