#!/usr/bin/env python
########################################################################
#
# (c) 2009 The Trustees of Columbia University in the City of New York.
# (c) 2026 Contributors to diffpy.cmipdf.
# All rights reserved.
#
# File coded by: Chris Farrow and members of the diffpy community.
#
# Originally developed in diffpy.srfit by the DANSE Diffraction group and
# Simon J. L. Billinge.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.cmipdf/graphs/contributors
#
# See LICENSE.rst and LICENSE_DANSE.rst for license information.
#
########################################################################
"""Example of fitting a crystal-like nanoparticle (nanocrystal) PDF.

This is an example of modeling the PDF from a nanocrystal as an
attenuated bulk PDF. This involves a crystal PDF calculation and a
spherical nanoparticle characteristic function. The equation we model is
Gnano(r) = f(r) * Gbulk(r), where f(r) is the nanoparticle
characteristic function for the nanoparticle shape. Functions for
calculating the characteristic function in the
diffpy.cmipdf.characteristicfunctions module.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from crystalpdf import scipy_optimize
from pyobjcryst.crystal import create_crystal_from_cif

from diffpy.cmipdf import PDFGenerator, PDFParser
from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
)


def make_recipe(ciffile, grdata):
    """Make a recipe to model a crystal-like nanoparticle PDF."""
    # Set up a PDF fit as has been done in other examples.
    pdfprofile = Profile()

    pdfparser = PDFParser()
    pdfparser.parse_file(grdata)
    pdfprofile.load_parsed_data(pdfparser)
    pdfprofile.set_calculation_range(xmin=0.1, xmax=20)

    pdfcontribution = FitContribution("pdf")
    pdfcontribution.set_profile(pdfprofile, xname="r")

    pdfgenerator = PDFGenerator("G")
    pdfgenerator.set_qmax(30.0)
    stru = create_crystal_from_cif(str(ciffile))
    pdfgenerator.set_structure(stru)
    pdfcontribution.add_profile_generator(pdfgenerator)

    # Register the nanoparticle shape factor.
    from diffpy.cmipdf.characteristicfunctions import spherical_particle

    pdfcontribution.register_function(spherical_particle, name="f")

    # Now we set up the fitting equation.
    pdfcontribution.set_equation("f * G")

    # Now make the recipe. Make sure we fit the characteristic function shape
    # parameters, in this case 'particle_diameter'.
    recipe = FitRecipe()
    recipe.add_contribution(pdfcontribution)

    phase = pdfgenerator.phase
    for par in phase.sgpars:
        recipe.add_variable(par)

    recipe.add_variable(pdfcontribution.particle_diameter, 20)
    recipe.add_variable(pdfgenerator.scale, 1)
    recipe.add_variable(pdfgenerator.delta2, 0)
    recipe.B11_0 = 0.1

    return recipe


def plot_results(recipe):
    """Plot the results contained within a refined FitRecipe."""
    r = recipe.pdf.profile.x
    g = recipe.pdf.profile.y

    # These two curves are not part of the standard observed/fit/diff plot
    # that plot_recipe produces, so we overlay them afterwards.
    gcryst = recipe.pdf.evaluate_equation("G")
    gcryst /= recipe.scale.value

    fr = recipe.pdf.evaluate_equation("f")
    fr *= max(g) / fr[0]

    fig, ax = recipe.plot_recipe(
        show=False,
        return_fig=True,
        xlabel=r"$r (\AA)$",
        ylabel=r"$G (\AA^{-2})$",
    )
    ax.plot(r, gcryst, "y--", label="G(r) Crystal")
    ax.plot(r, fr, "k--", label="f(r) calculated (scaled)")
    ax.legend(loc=1)

    plt.show()
    return


if __name__ == "__main__":

    ciffile = Path(__file__).parent / "data/pb.cif"
    grdata = Path(__file__).parent / "data/pb_100_qmin1.gr"

    recipe = make_recipe(ciffile, grdata)
    scipy_optimize(recipe)

    res = FitResults(recipe)
    res.print_results()

    plot_results(recipe)

# End of file
