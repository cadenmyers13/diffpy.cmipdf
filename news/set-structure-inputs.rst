**Added:**

* Let ``PDFGenerator.set_structure`` and ``DebyePDFGenerator.set_structure`` take a diffpy.cmistructure ParameterSet, used under its own name, or the path to a structure file, as well as a structure object.

**Changed:**

* The ``name`` argument of ``set_structure`` defaults to ``None``, which means "phase" for a structure object or a path. Giving a ``name`` with a ParameterSet raises a ``ValueError``.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
