import io
import pickle
from itertools import chain

import numpy as np
import pytest

from diffpy.cmipdf import PDFContribution
from diffpy.srfit.exceptions import SrFitError
from diffpy.structure import Structure, loadStructure


def test_set_qmax():
    """Check PDFContribution.setQmax()"""
    pc = PDFContribution("pdf")
    pc.set_qmax(21)
    pc.add_structure(Structure(), name="empty")
    assert 21 == pc.empty.get_qmax()
    pc.set_qmax(22)
    assert 22 == pc.get_qmax()
    assert 22 == pc.empty.get_qmax()
    return


def test_get_qmax():
    """Check PDFContribution.get_qmax()"""
    # cover all code branches in PDFContribution._get_meta_value
    # (1) contribution metadata
    pc1 = PDFContribution("pdf")
    assert pc1.get_qmax() is None
    pc1.set_qmax(17)
    assert 17 == pc1.get_qmax()
    # (2) contribution metadata
    pc2 = PDFContribution("pdf")
    pc2.add_structure(Structure(), name="empty")
    pc2.empty.set_qmax(18)
    assert 18 == pc2.get_qmax()
    # (3) profile metadata
    pc3 = PDFContribution("pdf")
    pc3.profile.meta["qmax"] = 19
    assert 19 == pc3.get_qmax()
    return


def test_savetxt(datafile):
    "check PDFContribution.savetxt()"

    pc = PDFContribution("pdf")
    pc.load_data(datafile("si-q27r60-xray.gr"))
    pc.set_calculation_range(0, 10)
    pc.add_structure(Structure(), name="empty")
    fp = io.BytesIO()
    with pytest.raises(SrFitError):
        pc.savetxt(fp)
    pc.evaluate()
    pc.savetxt(fp)
    txt = fp.getvalue().decode()
    nlines = len(txt.strip().split("\n"))
    assert 1001 == nlines
    return


def test_pickling(datafile):
    "validate PDFContribution.residual() after pickling."
    pc = PDFContribution("pdf")
    pc.load_data(datafile("ni-q27r100-neutron.gr"))
    ciffile = datafile("ni.cif")
    cif_path = str(ciffile)
    ni = loadStructure(cif_path)
    ni.Uisoequiv = 0.003
    pc.add_structure(ni, name="ni")
    pc.set_calculation_range(0, 10)
    pc2 = pickle.loads(pickle.dumps(pc))
    res0 = pc.residual()
    assert np.array_equal(res0, pc2.residual())
    for p in chain(
        pc.iterate_over_parameters("Uiso"), pc2.iterate_over_parameters("Uiso")
    ):
        p.value = 0.004
    res1 = pc.residual()
    assert not np.allclose(res0, res1)
    assert np.array_equal(res1, pc2.residual())
    return
