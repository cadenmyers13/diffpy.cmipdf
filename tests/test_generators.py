import numpy as np
import pytest

from diffpy.cmipdf import PDFGenerator
from diffpy.srreal.pdfcalculator import PDFCalculator
from diffpy.structure import PDFFitStructure


def testGenerator(datafile):
    qmax = 27.0
    gen = PDFGenerator()
    gen.set_scattering_type("N")
    assert "N" == gen.get_scattering_type()
    gen.set_qmax(qmax)
    assert qmax == pytest.approx(gen.get_qmax())

    structure = PDFFitStructure()
    ciffile = datafile("ni.cif")
    cif_path = str(ciffile)
    structure.read(cif_path)
    for i in range(4):
        structure[i].Bisoequiv = 1
    gen.set_structure(structure)

    calc = gen._calc
    # Test parameters
    for par in gen.iterate_over_parameters(recurse=False):
        pname = par.name
        defval = calc._getDoubleAttr(pname)
        assert defval == par.getValue()
        # Test setting values
        par.set_value(1.0)
        assert 1.0 == par.getValue()
        par.set_value(defval)
        assert defval == par.getValue()

    r = np.arange(0, 10, 0.1)
    y = gen(r)

    # Now create a reference PDF. Since the calculator is testing its
    # output, we just have to make sure we can calculate from the
    # PDFGenerator interface.

    calc = PDFCalculator()
    calc.rstep = r[1] - r[0]
    calc.rmin = r[0]
    calc.rmax = r[-1] + 0.5 * calc.rstep
    calc.qmax = qmax
    calc.setScatteringFactorTableByType("N")
    calc.eval(structure)
    yref = calc.pdf

    diff = y - yref
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)
    return


def test_set_qmin():
    """Verify qmin is propagated to the calculator object."""
    gen = PDFGenerator()
    assert 0 == gen.get_qmin()
    assert 0 == gen._calc.qmin
    gen.set_qmin(0.93)
    assert 0.93 == gen.get_qmin()
    assert 0.93 == gen._calc.qmin
    return
