"""Unit tests for __version__.py."""

import diffpy.cmipdf  # noqa


def test_package_version():
    """Ensure the package version is defined and not set to the initial
    placeholder."""
    assert hasattr(diffpy.cmipdf, "__version__")
    assert diffpy.cmipdf.__version__ != "0.0.0"
