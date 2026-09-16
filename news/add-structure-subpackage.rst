**Added:**

* ``diffpy.cmipdf.structure`` subpackage, migrated from ``diffpy.srfit.structure``.
  It adapts ``diffpy.structure``, ``pyobjcryst`` and ``cctbx`` structures to the
  ``ParameterSet`` interface and generates space group constraints. The migrated
  API is snake_case only; ``getValue`` and ``addParameter`` stay camelCase until
  core srfit renames them.

**Changed:**

* ``BasePDFGenerator`` now takes ``struToParameterSet`` from
  ``diffpy.cmipdf.structure`` instead of ``diffpy.srfit.structure``.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* ``src/diffpy/__init__.py`` now extends the ``diffpy`` namespace path, so
  ``diffpy.cmipdf`` no longer shadows the other ``diffpy`` packages.

**Security:**

* <news item>
