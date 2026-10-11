/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

// Helpers shared by the cmd*.cc binding files; not part of any API.

#pragma once

#include "pymol/platform/os_python.h"
#include "pymol/ui/but_mode.h"

#define API_SETUP_PYMOL_GLOBALS G = _api_get_pymol_globals(self)
#define API_SETUP_ARGS(G, self, args, ...)                  \
  if (!PyArg_ParseTuple(args, __VA_ARGS__)) return nullptr; \
  G = _api_get_pymol_globals(self);                         \
  API_ASSERT(G);

/**
 * Reports an error but keeps going
 */
#define API_HANDLE_ERROR               \
  if (PyErr_Occurred()) PyErr_Print(); \
  fprintf(stderr, "API-Error: in %s line %d.\n", __FILE__, __LINE__);

/**
 * If `x` is false, raises CmdException("x")
 */
#define API_ASSERT(x)                                                         \
  if (!(x)) {                                                                 \
    if (!PyErr_Occurred())                                                    \
      PyErr_SetString(P_CmdException ? P_CmdException : PyExc_Exception, #x); \
    return nullptr;                                                           \
  }

// Defined in cmd.cc. A binding brackets its work in APIEnter* and APIExit*,
// holding the API lock: APIEnter releases the GIL and APIExit takes it back,
// the Blocked variants keep it, and the NotModal variants return false without
// entering while a modal draw is active. The rest build return values: None on
// success, CmdException on failure.
PyMOLGlobals* _api_get_pymol_globals(PyObject* self);
void APIEnter(PyMOLGlobals* G);
int APIEnterNotModal(PyMOLGlobals* G);
void APIEnterBlocked(PyMOLGlobals* G);
int APIEnterBlockedNotModal(PyMOLGlobals* G);
void APIExit(PyMOLGlobals* G);
void APIExitBlocked(PyMOLGlobals* G);
PyObject* APISuccess(void);
PyObject* APIFailure(void);
PyObject* APIFailure(PyMOLGlobals* G, const char* msg = nullptr);
PyObject* APIFailure(PyMOLGlobals* G, const pymol::Error& error);
PyObject* APIResultCode(int code);
PyObject* APIResultOk(int ok);
PyObject* APIResultOk(PyMOLGlobals* G, bool ok);
PyObject* APIAutoNone(PyObject* result);

/**
 * If `res` is true, return res.result().
 * Else raise CmdException(res.error()).
 */
template <typename T>
PyObject* APIResult(PyMOLGlobals* G, pymol::Result<T>& res) {
  if (res) return PConvToPyObject(res.result());
  return APIFailure(G, res.error());
}

/// The bindings of each cmd_*.cc file, ending in a sentinel; PyInit__cmd
/// joins them with those in cmd.cc.
extern const PyMethodDef kEditingMethods[];
extern const PyMethodDef kViewingMethods[];
extern const PyMethodDef kQueryingMethods[];
extern const PyMethodDef kIoMethods[];
extern const PyMethodDef kCreatingMethods[];
extern const PyMethodDef kSettingsMethods[];
