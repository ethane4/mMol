/*
 * main_metal.cpp — Stubs for the Metal-only build, where main.cpp and
 * main_appkit.mm are excluded. Provides:
 *   - _gScaleFactor global
 *   - MainAsPyList / MainFromPyList / MainSavingUnderWhileIdle stubs
 */

#include "pymol/platform/os_python.h"
#include "pymol/base/pymol_globals.h"

int _gScaleFactor = 1;

int MainSavingUnderWhileIdle(void) { return 0; }

PyObject* MainAsPyList(PyMOLGlobals* G) { Py_RETURN_NONE; }

int MainFromPyList(PyMOLGlobals* G, PyObject* list) {
  return 1; /* success (no-op) */
}
