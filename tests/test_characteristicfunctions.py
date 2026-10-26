"""Tests for the nanoparticle characteristic functions."""

# FIXME: remove this line when `docformatter` fixes the blank line bug


import re
import unittest
import warnings

import numpy
import numpy.testing as npt
import pytest

import diffpy.cmipdf.characteristicfunctions as cf

# ----------------------------------------------------------------------------
# The characteristic functions below compute the attenuation of the PDF due
# to a finite nanoparticle shape. Each has a known value of 1 at r = 0 (the
# particle always overlaps itself at zero distance) and drops to 0 once r
# reaches the largest distance spanned by the shape.


@pytest.mark.parametrize(
    "input_r, input_particle_diameter, expected_characteristic_function",
    [
        # C1: r ranges over values inside and outside the particle.
        # Expected: the polynomial form 1 - 1.5x + 0.5x^3 (x = r /
        # particle_diameter) applies inside the particle, and the function
        # is zero at and beyond the particle diameter.
        (
            numpy.array([0.0, 5.0, 10.0, 15.0]),
            10.0,
            [1.0, 0.3125, 0.0, 0.0],
        ),
    ],
)
def test_spherical_particle(
    input_r, input_particle_diameter, expected_characteristic_function
):
    actual_characteristic_function = cf.spherical_particle(
        input_r, input_particle_diameter
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_equatorial_radius, input_polar_radius, "
    "expected_characteristic_function",
    [
        # C1: equatorial_radius equals polar_radius, i.e. a sphere.
        # Expected: matches spherical_particle for the equivalent diameter.
        (
            numpy.linspace(0, 25, 6),
            10.0,
            10.0,
            cf.spherical_particle(numpy.linspace(0, 25, 6), 20.0),
        ),
        # C2: polar_radius > equatorial_radius, i.e. a prolate spheroid.
        # Expected: values from the Lei et al. characteristic function.
        (
            numpy.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0]),
            10.0,
            20.0,
            [
                1.0,
                0.40106264900760513,
                0.054200238412168256,
                0.003650768214115807,
                0.0,
                0.0,
            ],
        ),
        # C3: polar_radius < equatorial_radius, i.e. an oblate spheroid.
        # Expected: values from the Lei et al. characteristic function.
        (
            numpy.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0]),
            20.0,
            10.0,
            [
                1.0,
                0.744181556741916,
                0.5061470768546105,
                0.30368052370886184,
                0.15456586067544859,
                0.0675583474738014,
            ],
        ),
    ],
)
def test_spheroidal_particle(
    input_r,
    input_equatorial_radius,
    input_polar_radius,
    expected_characteristic_function,
):
    actual_characteristic_function = cf.spheroidal_particle(
        input_r, input_equatorial_radius, input_polar_radius
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_particle_diameter, input_particle_diameter_sigma, "
    "expected_characteristic_function",
    [
        # C1: the lognormal distribution has no width.
        # Expected: matches spherical_particle for the same diameter.
        (
            numpy.array([0.0, 5.0, 10.0]),
            10.0,
            0.0,
            cf.spherical_particle(numpy.array([0.0, 5.0, 10.0]), 10.0),
        ),
        # C2: a genuine lognormal size distribution is applied.
        # Expected: values from the closed-form lognormal-averaged
        # spherical characteristic function.
        (
            numpy.array([1.0, 5.0, 10.0, 15.0]),
            10.0,
            2.0,
            [
                0.8617610662266728,
                0.362144891493816,
                0.039068792680198694,
                0.0009056141323687261,
            ],
        ),
    ],
)
def test_lognormal_spherical_particle(
    input_r,
    input_particle_diameter,
    input_particle_diameter_sigma,
    expected_characteristic_function,
):
    actual_characteristic_function = cf.lognormal_spherical_particle(
        input_r, input_particle_diameter, input_particle_diameter_sigma
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_sheet_thickness, expected_characteristic_function",
    [
        # C1: r is an array with values inside and outside the sheet
        # thickness.
        # Expected: 1 - 0.5*r/thickness inside the sheet, and
        # 0.5*thickness/r beyond it.
        (
            numpy.array([0.0, 2.0, 4.0, 8.0]),
            4.0,
            [1.0, 0.75, 0.5, 0.25],
        ),
        # C2: r is a scalar inside the sheet thickness.
        # Expected: 1 - 0.5*r/thickness.
        (2.0, 4.0, 0.75),
    ],
)
def test_sheet_particle(
    input_r, input_sheet_thickness, expected_characteristic_function
):
    actual_characteristic_function = cf.sheet_particle(
        input_r, input_sheet_thickness
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_radius, input_thickness, "
    "expected_characteristic_function",
    [
        # C1: r ranges from the shell center out past the outer radius.
        # Expected: values from the Lei et al. shell characteristic
        # function; the function is 1 at r = 0 and 0 once r exceeds the
        # outer radius (radius + thickness).
        (
            numpy.array([0.0, 5.0, 12.5, 15.0, 20.0, 30.0]),
            10.0,
            5.0,
            [
                1.0,
                0.4934210526315789,
                0.19736842105263158,
                0.16447368421052633,
                0.12335526315789473,
                0.0,
            ],
        ),
    ],
)
def test_shell_particle(
    input_r, input_radius, input_thickness, expected_characteristic_function
):
    actual_characteristic_function = cf.shell_particle(
        input_r, input_radius, input_thickness
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


# ----------------------------------------------------------------------------
# A refinement can step a shape parameter into a region that has no physical
# meaning, such as a negative diameter or thickness. The characteristic
# functions must not raise there, because that would abort a refinement that
# would otherwise recover on its own, so they return zero for every r instead.
# A zero return is flat in the shape parameter and so stalls the optimizer,
# which the user cannot see from the fit output alone; each function therefore
# also warns. The warning is issued only once per process for each distinct
# problem so that a refinement loop does not bury the user in repeated
# messages.
#
# Every such warning states which requirement was violated and then closes
# with the same explanation of the consequence and of how to avoid it. Only
# the requirement varies between cases.


@pytest.mark.parametrize(
    "function_name, input_args, expected_requirement",
    [
        # C1: the sphere diameter is zero.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "spherical_particle",
            (numpy.array([0.0, 5.0]), 0.0),
            "'particle_diameter' must be positive",
        ),
        # C2: the sphere diameter is negative.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), -10.0),
            "'particle_diameter' must be positive",
        ),
        # C3: the spheroid equatorial radius is zero, which the axis ratio
        # would otherwise be divided by.
        # Expected: zero everywhere, warning that the equatorial radius must
        # be positive, and no ZeroDivisionError.
        (
            "spheroidal_particle",
            (numpy.array([1.0, 5.0]), 0.0, 5.0),
            "'equatorial_radius' must be positive",
        ),
        # C4: the spheroid polar radius is negative.
        # Expected: zero everywhere, warning that the polar radius must be
        # positive.
        (
            "spheroidal_particle",
            (numpy.array([1.0, 5.0]), 10.0, -5.0),
            "'polar_radius' must be positive",
        ),
        # C5: the mean particle diameter of the lognormal distribution is
        # zero.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 2.0]), 0.0, 1.0),
            "'particle_diameter' must be positive",
        ),
        # C6: the width of the lognormal distribution is negative. A negative
        # width is meaningless rather than a limiting case, so it must not be
        # quietly treated as a distribution of zero width.
        # Expected: zero everywhere, warning that the width must not be
        # negative.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -2.0),
            "'particle_diameter_sigma' must not be negative",
        ),
        # C7: the mean particle diameter is too small for the requested
        # distribution width, which puts the underlying lognormal mean below
        # zero and invalidates the closed form.
        # Expected: zero everywhere, warning that the diameter is too small
        # for the width.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0]), 2.0, 4.0),
            "'particle_diameter' is too small for 'particle_diameter_sigma'",
        ),
        # C8: the sheet thickness is zero.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "sheet_particle",
            (numpy.array([1.0, 2.0]), 0.0),
            "'sheet_thickness' must be positive",
        ),
        # C9: the sheet thickness is negative.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "sheet_particle",
            (numpy.array([1.0, 2.0]), -1.0),
            "'sheet_thickness' must be positive",
        ),
        # C10: the shell thickness is zero, which would otherwise make the
        # normalizing denominator vanish and return an unattenuated one at
        # every r.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 0.0),
            "'thickness' must be positive",
        ),
        # C11: the shell thickness is negative, which would otherwise return
        # a smoothly varying branch with values below negative one that an
        # optimizer can mistake for a real solution.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -5.0),
            "'thickness' must be positive",
        ),
        # C12: the shell inner radius is negative.
        # Expected: zero everywhere, warning that the radius must not be
        # negative.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), -20.0, 5.0),
            "'radius' must not be negative",
        ),
    ],
)
def test_non_physical_input_returns_zero_and_warns(
    function_name,
    input_args,
    expected_requirement,
    reset_characteristic_function_warnings,
):
    function = getattr(cf, function_name)
    input_r = input_args[0]
    expected_remedy = (
        "The characteristic function was set to zero for every r, which "
        "flattens the fit residual so a refinement cannot recover on its own. "
        "Please use a physically meaningful starting value, or keep the "
        "parameter in range with FitRecipe.add_soft_bounds."
    )
    expected_msg = (
        f"In '{function_name}', {expected_requirement}. {expected_remedy}"
    )
    with pytest.warns(RuntimeWarning, match=re.escape(expected_msg)):
        actual_characteristic_function = function(*input_args)
    npt.assert_allclose(
        actual_characteristic_function, numpy.zeros_like(input_r)
    )


@pytest.mark.parametrize(
    "function_name, input_args, expected_characteristic_function",
    [
        # C1: a lognormal distribution of zero width is the sphere limit
        # rather than a mistake.
        # Expected: matches spherical_particle for the same diameter, with no
        # warning.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 0.0),
            cf.spherical_particle(numpy.array([1.0, 5.0, 10.0]), 10.0),
        ),
        # C2: a shell of zero inner radius is a solid sphere rather than a
        # mistake.
        # Expected: matches spherical_particle for twice the thickness, with
        # no warning.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 0.0, 5.0),
            cf.spherical_particle(numpy.array([1.0, 5.0, 10.0]), 10.0),
        ),
    ],
)
def test_limiting_input_does_not_warn(
    function_name,
    input_args,
    expected_characteristic_function,
    reset_characteristic_function_warnings,
):
    function = getattr(cf, function_name)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual_characteristic_function = function(*input_args)
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


def test_non_physical_input_warns_only_once(
    reset_characteristic_function_warnings,
):
    # A refinement evaluates the characteristic function on every iteration,
    # so a function that warns on each call would flood the user. The
    # suppression must not rely on the warnings registry, which callers can
    # disable, so this asserts under the "always" filter.
    with warnings.catch_warnings(record=True) as raised_warnings:
        warnings.simplefilter("always")
        for _ in range(5):
            cf.spherical_particle(numpy.array([1.0, 2.0]), -10.0)
    actual_warning_count = len(raised_warnings)
    expected_warning_count = 1
    assert actual_warning_count == expected_warning_count


if __name__ == "__main__":
    unittest.main()
