=============
Release notes
=============

.. current developments

0.1.1
=====

**Added:**

* Examples in ``docs/examples``, ported from the ``diffpy.srfit`` PDF examples: ``crystalpdf``, ``simplepdf``, ``parallel_refinement``, ``crystalpdfobjcryst``, ``crystalpdftwophase``, ``simplepdftwophase``, ``crystalpdftwodata``, ``crystalpdfall``, ``nppdfcrystal``, ``nppdfobjcryst`` and ``coreshellnp``.
* ``diffpy.cmistructure`` as a dependency.

**Changed:**

* Replace the placeholder description in the README with a description of the package.
* Rename the parameters of the characteristic functions ``spherical_particle``, ``spheroidal_particle``, ``lognormal_spherical_particle``, ``sheet_particle``, and ``shell_particle`` to match ``diffpy.srfit.pdf.characteristicfunctions``.
* Issue a ``RuntimeWarning`` once per process when a characteristic function is given non-physical input and returns zero.
* Use the same license as ``diffpy.srfit``: ``LICENSE.rst``, the documentation License page and the ``License :: Free To Use But Restricted`` classifier now match it.

**Fixed:**

* Read the PDF metadata in the data file header in ``PDFContribution.load_data``, which used the generic ``ProfileParser`` and lost the scattering type and qmax from NOMAD files with a free-text instrument comment.
* The README now gives the command-line program as ``diffpy-cmipdf``, which is its installed name, instead of ``diffpy.cmipdf``.
* Require ``diffpy.srfit`` 3.4.0 or later. Older versions lack the methods diffpy.cmipdf calls, so installing it next to ``diffpy.srfit`` 3.3.1 failed at run time.
* Build ``diffpy.cmistructure`` parameter sets in ``set_structure`` and ``add_structure``, and accept them in ``set_structure_from_parset`` and ``add_structure_from_parset``. These used ``diffpy.srfit.structure`` and raised ``AttributeError: stru`` for a ``diffpy.cmistructure`` parameter set.

**Removed:**

* The scikit-package template pages (Getting started, the example table snippet and the logo) from the documentation.
* Support for Python 3.11. ``diffpy.srfit`` 3.4.0, which diffpy.cmipdf requires, needs Python 3.12 or later.
* Remove the non-functional ``SASCF`` calculator from ``characteristicfunctions``. SAS characteristic functions now live in ``diffpy.cmisas``.


0.1.0
=====

**Added:**

* Add source code from ``diffpy.srfit.pdf`` to ``diffpy.cmipdf``.
