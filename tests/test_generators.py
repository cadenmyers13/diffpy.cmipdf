"""Tests for the pdfgenerator module."""

import numpy as np
import pytest

from diffpy.cmipdf import DebyePDFGenerator, PDFGenerator
from diffpy.cmistructure import constrain_as_space_group
from diffpy.cmistructure.diffpyparset import DiffpyStructureParSet
from diffpy.srreal.pdfcalculator import PDFCalculator


# set_scattering_type stores the radiation type used for the calculation.
@pytest.mark.parametrize(
    "input_scattering_type, expected_scattering_type",
    [
        # C1: The scattering type is set to neutron.
        # Expected: get_scattering_type returns "N".
        ("N", "N"),
        # C2: The scattering type is set to x-ray.
        # Expected: get_scattering_type returns "X".
        ("X", "X"),
    ],
)
def test_set_scattering_type(input_scattering_type, expected_scattering_type):
    generator = PDFGenerator()
    generator.set_scattering_type(input_scattering_type)
    actual_scattering_type = generator.get_scattering_type()
    assert actual_scattering_type == expected_scattering_type


def test_set_qmax():
    # C1: qmax is set to 27 on a new generator.
    # Expected: get_qmax returns 27.
    generator = PDFGenerator()
    generator.set_qmax(27.0)
    actual_qmax = generator.get_qmax()
    expected_qmax = 27.0
    assert actual_qmax == pytest.approx(expected_qmax)


def test_set_qmin():
    generator = PDFGenerator()
    # C1: qmin is read from a new generator.
    # Expected: qmin is 0 on both the generator and its calculator.
    actual_qmins = (generator.get_qmin(), generator._calc.qmin)
    expected_qmins = (0, 0)
    assert actual_qmins == expected_qmins
    # C2: qmin is set to 0.93.
    # Expected: qmin is 0.93 on both the generator and its calculator.
    generator.set_qmin(0.93)
    actual_qmins = (generator.get_qmin(), generator._calc.qmin)
    expected_qmins = (0.93, 0.93)
    assert actual_qmins == expected_qmins


# set_structure adapts the structure with a diffpy.cmistructure parameter
# set, which the generator then manages as its phase.
@pytest.mark.parametrize(
    "input_generator_class",
    [
        # C1: A diffpy.structure Structure is set on a generator,
        # 1. PDFGenerator, 2. DebyePDFGenerator.
        # Expected: The phase is a diffpy.cmistructure DiffpyStructureParSet
        # that adapts the structure.
        PDFGenerator,
        DebyePDFGenerator,
    ],
)
def test_set_structure(ni_pdffit_structure, input_generator_class):
    generator = input_generator_class()
    generator.set_structure(ni_pdffit_structure)
    actual_phase = (type(generator.phase), generator.phase.structure)
    expected_phase = (DiffpyStructureParSet, ni_pdffit_structure)
    assert actual_phase == expected_phase


# set_structure_from_parset manages the given diffpy.cmistructure parameter
# set as the phase, and periodic sets whether it uses symmetry.
@pytest.mark.parametrize(
    "input_periodic, expected_using_symmetry",
    [
        # C1: The parameter set is set as periodic.
        # Expected: The generator manages that parameter set, its structure
        # is the adapted structure, and the parameter set uses symmetry.
        (True, True),
        # C2: The parameter set is set as not periodic.
        # Expected: The same, except that it does not use symmetry.
        (False, False),
    ],
)
def test_set_structure_from_parset(
    ni_pdffit_structure, input_periodic, expected_using_symmetry
):
    parameter_set = DiffpyStructureParSet("phase", ni_pdffit_structure)
    generator = PDFGenerator()
    generator.set_structure_from_parset(parameter_set, periodic=input_periodic)
    actual_state = (
        generator.phase is parameter_set,
        generator.structure is ni_pdffit_structure,
        parameter_set.using_symmetry(),
    )
    expected_state = (True, True, expected_using_symmetry)
    assert actual_state == expected_state


def test_constrain_as_space_group_accepts_generator(ni_pdffit_structure):
    # C1: A PDFGenerator holding face-centered cubic Ni is passed to
    # diffpy.cmistructure's constrain_as_space_group with Fm-3m.
    # Expected: Only the lattice length a of its phase is left free.
    generator = PDFGenerator()
    generator.set_structure(ni_pdffit_structure)
    space_group_parameters = constrain_as_space_group(generator, "Fm-3m")
    actual_free_lattice_names = [
        parameter.name for parameter in space_group_parameters.latpars
    ]
    expected_free_lattice_names = ["a"]
    assert actual_free_lattice_names == expected_free_lattice_names


def test_parameter_names(ni_pdffit_structure):
    # C1: A structure is set on the generator.
    # Expected: The generator holds one Parameter for each calculator
    # attribute it wraps.
    generator = PDFGenerator()
    generator.set_structure(ni_pdffit_structure)
    actual_parameter_names = [
        parameter.name
        for parameter in generator.iterate_over_parameters(recurse=False)
    ]
    expected_parameter_names = ["delta1", "delta2", "qbroad", "scale", "qdamp"]
    assert actual_parameter_names == expected_parameter_names


# Each generator Parameter wraps the calculator attribute of the same name.
@pytest.mark.parametrize(
    "input_parameter_name",
    [
        # C1: A wrapped calculator attribute is read and then set to 1,
        # 1. delta1, 2. delta2, 3. qbroad, 4. scale, 5. qdamp.
        # Expected: The Parameter starts at the calculator's value, and
        # setting it to 1 sets the calculator attribute to 1.
        "delta1",
        "delta2",
        "qbroad",
        "scale",
        "qdamp",
    ],
)
def test_parameter_wraps_calculator(ni_pdffit_structure, input_parameter_name):
    generator = PDFGenerator()
    generator.set_structure(ni_pdffit_structure)
    parameter = generator.get(input_parameter_name)
    actual_initial_value = parameter.get_value()
    expected_initial_value = generator._calc._getDoubleAttr(
        input_parameter_name
    )
    assert actual_initial_value == expected_initial_value
    parameter.set_value(1.0)
    actual_values = (
        parameter.get_value(),
        generator._calc._getDoubleAttr(input_parameter_name),
    )
    expected_values = (1.0, 1.0)
    assert actual_values == expected_values


def test_generator_matches_calculator(ni_pdffit_structure):
    # C1: The generator computes the neutron PDF of Ni with qmax = 27 on
    # r = 0, 0.1, ..., 9.9.
    # Expected: The PDF equals the one from a PDFCalculator with the same
    # settings.
    qmax = 27.0
    r = np.arange(0, 10, 0.1)
    generator = PDFGenerator()
    generator.set_scattering_type("N")
    generator.set_qmax(qmax)
    generator.set_structure(ni_pdffit_structure)
    actual_pdf = generator(r)
    calculator = PDFCalculator()
    calculator.rstep = r[1] - r[0]
    calculator.rmin = r[0]
    calculator.rmax = r[-1] + 0.5 * calculator.rstep
    calculator.qmax = qmax
    calculator.scatteringfactortable = "N"
    calculator.eval(ni_pdffit_structure)
    expected_pdf = calculator.pdf
    assert np.allclose(actual_pdf, expected_pdf)
