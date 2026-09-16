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
"""Modules and classes that adapt structure representations to the
ParameterSet interface and automatic structure constraint generation
from space group information."""

from diffpy.cmipdf.structure.sgconstraints import constrain_as_space_group

__all__ = ["constrain_as_space_group", "struToParameterSet"]


def struToParameterSet(name, stru):
    """Creates a ParameterSet from an structure.

    This returns a ParameterSet adapted for the structure depending on its
    type.

    Parameters
    ----------
    stru
        a structure object known by this module
    name
        A name to give the structure.

    Raises TypeError if stru cannot be adapted
    """
    from diffpy.cmipdf.structure.diffpyparset import DiffpyStructureParSet

    if DiffpyStructureParSet.can_adapt(stru):
        return DiffpyStructureParSet(name, stru)

    from diffpy.cmipdf.structure.objcrystparset import ObjCrystCrystalParSet

    if ObjCrystCrystalParSet.can_adapt(stru):
        return ObjCrystCrystalParSet(name, stru)

    from diffpy.cmipdf.structure.objcrystparset import ObjCrystMoleculeParSet

    if ObjCrystMoleculeParSet.can_adapt(stru):
        return ObjCrystMoleculeParSet(name, stru)

    from diffpy.cmipdf.structure.cctbxparset import CCTBXCrystalParSet

    if CCTBXCrystalParSet.can_adapt(stru):
        return CCTBXCrystalParSet(name, stru)

    raise TypeError("Unadaptable structure format")


# silence pyflakes checker
assert constrain_as_space_group

# End of file
