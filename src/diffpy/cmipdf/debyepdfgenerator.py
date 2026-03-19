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
"""PDF profile generator using the Debye equation.

The DebyePDFGenerator class can take a diffpy.structure,
pyobjcryst.crystal.Crystal or pyobjcryst.molecule.Molecule object and
calculate the PDF from it. This generator is especially appropriate for
isolated scatterers, such as nanoparticles and molecules.
"""

__all__ = ["DebyePDFGenerator"]

from pathlib import Path

from diffpy.cmipdf.basepdfgenerator import BasePDFGenerator
from diffpy.srreal.pdfcalculator import DebyePDFCalculator
from diffpy.structure import loadStructure


class DebyePDFGenerator(BasePDFGenerator):
    """A class for calculating the PDF from an isolated scatterer.

    This works with diffpy.structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes
    ----------
    _calc
        DebyePDFCalculator instance for calculating the PDF
    _phase
        The structure ParameterSets used to calculate the profile.
    structure
        The structure objected adapted by _phase.
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
        The structure ParameterSet (SrRealStructure instance) used to calculate
        the profile is named by the user.

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

    def set_structure(self, structure, name="phase", periodic=False):
        """Set the structure that will be used to calculate the PDF.

        This creates a DiffpyStructureParSet, ObjCrystCrystalParSet or
        ObjCrystMoleculeParSet that adapts structure to a ParameterSet
        interface.
        See those classes (located in diffpy.srfit.structure) for how they are
        used. The resulting ParameterSet will be managed by this generator.

        Parameters
        ----------
        structure : Structure object
            The `diffpy.structure.Structure`, `pyobjcryst.crystal.Crystal` or
            `pyobjcryst.molecule.Molecule` instance.
        name : str, optional
            A name to give to the managed ParameterSet that adapts structure
            (default "phase").
        periodic : bool, optional
            The structure should be treated as periodic (default
            False). Note that some structures do not support
            periodicity, in which case this will have no effect on the
            PDF calculation.
        """
        return BasePDFGenerator.set_structure(self, structure, name, periodic)

    def set_structure_from_parset(self, parset, periodic=False):
        """Set the phase that will be used to calculate the PDF.

        Set the phase directly with a DiffpyStructureParSet,
        ObjCrystCrystalParSet or ObjCrystMoleculeParSet that adapts a structure
        object (from diffpy or pyobjcryst).  The passed ParameterSet will be
        managed by this generator.

        Parameters
        ----------
        parset : SrealParSet object
            The SrRealParSet that holds the structural information.
            This can be used to share the phase between multiple
            BasePDFGenerators, and have the changes in one reflect in
            another.
        periodic : bool, optional
            The structure should be treated as periodic (default True).
            Note that some structures do not support periodicity, in
            which case this will be ignored.
        """
        return BasePDFGenerator.set_structure_from_parset(
            self, parset, periodic
        )

    def __init__(self, name="pdf"):
        """Initialize the generator."""
        BasePDFGenerator.__init__(self, name)
        self._set_calculator(DebyePDFCalculator())
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
        uiso : float, optional
            The isotropic atomic displacement parameter to use for all
            atoms in the structure (default 0.007).
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
            xyz_path = "path/to/nanoparticle.xyz"
            gen = PDFGenerator()
            r, g = gen.generate_pdf_from_structure(xyz_path)
        """
        if isinstance(structure, Path):
            structure = loadStructure(str(structure))
        elif isinstance(structure, str):
            structure = loadStructure(structure)
        structure.Uisoequiv = uiso

        calc = DebyePDFCalculator()
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


# End class DebyePDFGenerator

# End of file
