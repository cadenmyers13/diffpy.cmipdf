"""Tests for the pdfcontribution module."""

import io
import pickle
import re

import numpy as np
import pytest

from diffpy.cmipdf import PDFContribution
from diffpy.srfit.exceptions import SrFitError
from diffpy.structure import Structure, load_structure


def test_set_qmax():
    contribution = PDFContribution("pdf")
    # C1: qmax is set to 21 before a structure is added.
    # Expected: The generator created for the structure gets qmax = 21.
    contribution.set_qmax(21)
    contribution.add_structure(Structure(), name="empty")
    actual_generator_qmax = contribution.empty.get_qmax()
    expected_generator_qmax = 21
    assert actual_generator_qmax == expected_generator_qmax
    # C2: qmax is set to 22 after the structure is added.
    # Expected: Both the contribution and its generator have qmax = 22.
    contribution.set_qmax(22)
    actual_qmaxes = (contribution.get_qmax(), contribution.empty.get_qmax())
    expected_qmaxes = (22, 22)
    assert actual_qmaxes == expected_qmaxes


# get_qmax looks for qmax in the contribution, then in its generators, then
# in its profile.
def test_get_qmax_not_set():
    # C1: qmax is not set anywhere.
    # Expected: Returns None.
    contribution = PDFContribution("pdf")
    actual_qmax = contribution.get_qmax()
    expected_qmax = None
    assert actual_qmax == expected_qmax


def test_get_qmax_from_contribution():
    # C1: qmax is set to 17 on the contribution.
    # Expected: Returns 17.
    contribution = PDFContribution("pdf")
    contribution.set_qmax(17)
    actual_qmax = contribution.get_qmax()
    expected_qmax = 17
    assert actual_qmax == expected_qmax


def test_get_qmax_from_generator():
    # C1: qmax is set to 18 on the contribution's generator only.
    # Expected: Returns 18.
    contribution = PDFContribution("pdf")
    contribution.add_structure(Structure(), name="empty")
    contribution.empty.set_qmax(18)
    actual_qmax = contribution.get_qmax()
    expected_qmax = 18
    assert actual_qmax == expected_qmax


def test_get_qmax_from_profile():
    # C1: qmax is set to 19 in the profile metadata only.
    # Expected: Returns 19.
    contribution = PDFContribution("pdf")
    contribution.profile.meta["qmax"] = 19
    actual_qmax = contribution.get_qmax()
    expected_qmax = 19
    assert actual_qmax == expected_qmax


def test_savetxt(datafile):
    # C1: A contribution with data on r = 0 to 10 is evaluated and saved.
    # Expected: The file has a header line and one line for each of the
    # 1000 points.
    contribution = PDFContribution("pdf")
    contribution.load_data(datafile("si-q27r60-xray.gr"))
    contribution.set_calculation_range(0, 10)
    contribution.add_structure(Structure(), name="empty")
    contribution.evaluate()
    file_object = io.BytesIO()
    contribution.savetxt(file_object)
    actual_text = file_object.getvalue().decode()
    actual_line_count = len(actual_text.strip().split("\n"))
    expected_line_count = 1001
    assert actual_line_count == expected_line_count


def test_savetxt_bad(datafile):
    # C1: A contribution is saved before it is evaluated.
    # Expected: Raises an SrFitError, since there is no calculated profile.
    contribution = PDFContribution("pdf")
    contribution.load_data(datafile("si-q27r60-xray.gr"))
    contribution.set_calculation_range(0, 10)
    contribution.add_structure(Structure(), name="empty")
    expected_msg = "ycalc is None"
    with pytest.raises(SrFitError, match=re.escape(expected_msg)):
        contribution.savetxt(io.BytesIO())


def _make_ni_contribution(datafile, uiso):
    contribution = PDFContribution("pdf")
    contribution.load_data(datafile("ni-q27r100-neutron.gr"))
    ni = load_structure(str(datafile("ni.cif")))
    ni.Uisoequiv = uiso
    contribution.add_structure(ni, name="ni")
    contribution.set_calculation_range(0, 10)
    return contribution


def test_pickling(datafile):
    contribution = _make_ni_contribution(datafile, uiso=0.003)
    unpickled_contribution = pickle.loads(pickle.dumps(contribution))
    # C1: A Ni contribution is pickled and unpickled.
    # Expected: The copy has the same residual as the original.
    actual_residual = unpickled_contribution.residual()
    expected_residual = contribution.residual()
    assert np.array_equal(actual_residual, expected_residual)
    # C2: Uiso of every atom in the copy is changed from 0.003 to 0.004.
    # Expected: The copy's residual equals that of a contribution built
    # with Uiso = 0.004.
    for parameter in unpickled_contribution.iterate_over_parameters("Uiso"):
        parameter.value = 0.004
    actual_residual = unpickled_contribution.residual()
    expected_residual = _make_ni_contribution(datafile, uiso=0.004).residual()
    assert np.allclose(actual_residual, expected_residual)
