#!/usr/bin/env python
##############################################################################
#
# (c) 2025 Simon Billinge.
# All rights reserved.
#
# File coded by: Caden Myers, Simon Billinge, and members of the Billinge
#                group.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.cmipdf/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""PDF profile generator.

The PDFGenerator class can take a diffpy.structure,
pyobjcryst.crystal.Crystal or pyobjcryst.molecule.Molecule object and
calculate the crystal PDF from it. The passed structure object is
wrapped in a StructureParameter set, which makes its attributes
refinable. See the class definition for more details and the examples
for its use.
"""

__all__ = ["PDFGenerator"]

from pathlib import Path

from diffpy.cmipdf.basepdfgenerator import BasePDFGenerator
from diffpy.srreal.pdfcalculator import PDFCalculator
from diffpy.structure import loadStructure


class PDFGenerator(BasePDFGenerator):
    """A class for calculating the PDF from a single crystal structure.

    This works with diffpy.structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes
    ----------
    _calc
        PDFCalculator instance for calculating the PDF
    _phase
        The structure ParameterSet used to calculate the profile.
    _lastr
        The last value of r over which the PDF was calculated. This is
        used to configure the calculator when r changes.

    Managed Parameters
    ------------------
    scale
        Scale factor
    delta1
        Linear peak broadening term
    delta2
        Quadratic peak broadening term
    qbroad
        Resolution peak broadening term
    qdamp
        Resolution peak dampening term
    Managed ParameterSets
        The structure ParameterSet (SrRealStructure instance) used to
        calculate the profile is named by the user.

    Usable Metadata
    ---------------
    stype
        The scattering type "X" for x-ray, "N" for neutron (see
        'set_scattering_type').
    qmax
        The maximum scattering vector used to generate the PDF (see
        set_qmax).
    qmin
        The minimum scattering vector used to generate the PDF (see
        set_qmin).
    scale
        See Managed Parameters.
    delta1
        See Managed Parameters.
    delta2
        See Managed Parameters.
    qbroad
        See Managed Parameters.
    qdamp
        See Managed Parameters.
    """

    def __init__(self, name="pdf"):
        """Initialize the generator."""
        BasePDFGenerator.__init__(self, name)
        self._set_calculator(PDFCalculator())
        return

    def generate_pdf_from_structure(
        self,
        structure,
        rmin=0,
        rmax=30,
        qmin=0.1,
        qmax=25.0,
        qdamp=0.03,
        qbroad=0.0,
        delta1=0.0,
        delta2=0.0,
        uiso=0.007,
    ):
        """Calculate the PDF from a structure and return G vs. r.

        This is a convenience method that allows the user to calculate
        the PDF from a
        structure without having to set up the calculator and
        structure ParameterSet
        manually. The structure can be passed as a path to a
        structure file or as a
        `diffpy.structure.Structure` object.

        Parameters
        ----------
        structure : Path, str, or diffpy.structure.Structure
            The structure to calculate the PDF from. Can be a path
            to a structure file or a diffpy.structure.Structure
            object.
        rmin : float, optional
            The minimum r value in Angstroms for the PDF (default 0).
        rmax : float, optional
            The maximum r value in Angstroms for the PDF (default 30).
        qmin : float, optional
            The minimum scattering vector used to generate the PDF
            (default 0.1).
        qmax : float, optional
            The maximum scattering vector used to generate the PDF
            (default 25.0).
        qdamp : float, optional
            The resolution dampening term to use in the PDF calculation
            (default 0.03).
        qbroad : float, optional
            The resolution broadening term to use in the PDF calculation
            (default 0.0).
        delta1 : float, optional
            The linear peak broadening term to use in the PDF calculation
            (default 0.0).
        delta2 : float, optional
            The quadratic peak broadening term to use in the PDF calculation
            (default 0.0).
        uiso : float, optional
            The isotropic atomic displacement parameter to use for all
            atoms in the structure (default 0.007).

        Returns
        -------
        r : numpy.ndarray
            The r values for the PDF in units of Angstroms.
        G : numpy.ndarray
            The G values for the PDF in units of 1/Angstrom^2.

        Example
        -------
        .. code-block:: python
            cif_path = "path/to/ni.cif"
            gen = PDFGenerator()
            r, g = gen.generate_pdf_from_structure(cifpath)
        """
        if isinstance(structure, Path):
            structure = loadStructure(str(structure))
        elif isinstance(structure, str):
            structure = loadStructure(structure)
        structure.Uisoequiv = uiso

        calc = PDFCalculator()
        calc.qmin = qmin
        calc.qmax = qmax
        calc.rstep = 0.01
        calc.rmin = rmin
        calc.rmax = rmax
        calc.qdamp = qdamp
        calc.qbroad = qbroad
        calc.delta1 = delta1
        calc.delta2 = delta2

        r, g = calc(structure)
        return r, g


# End class PDFGenerator
