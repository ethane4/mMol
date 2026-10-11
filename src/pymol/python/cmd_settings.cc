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

// Python bindings for settings and feedback (pymol/setting.py, feedingback.py).

#include "pymol/base/err.h"
#include "pymol/executive/executive.h"
#include "pymol/executive/setting_effects.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/feedback.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdGetSettingUpdates(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int state;
  char* name;

  if (!PyArg_ParseTuple(args, "Osi", &self, &name, &state)) {
    API_HANDLE_ERROR;
    ok_raise(2);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  result = PConvToPyObject(SettingGetUpdateList(G, name, state));

  APIExitBlocked(G);
ok_except2:
  return (APIAutoNone(result));
}

static PyObject* CmdGetSettingIndices(PyObject* self, PyObject* args) {
  return SettingGetSettingIndices();
}

static PyObject* CmdFeedback(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1, i2, result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &i1, &i2);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    /* NO API Entry for performance,
     *feedback (MACRO) just accesses a safe global */
    result = Feedback(G, i1, i2);
  }
  return Py_BuildValue("i", result);
}

static PyObject* CmdSetFeedbackMask(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1, i2, i3;
  API_SETUP_ARGS(G, self, args, "Oiii", &self, &i1, &i2, &i3);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetFeedbackMask(G, i1, static_cast<unsigned char>(i2),
                           static_cast<unsigned char>(i3));
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdGetSettingType(PyObject*, PyObject* args) {
  int index;
  if (!PyArg_ParseTuple(args, "i", &index)) {
    return nullptr;
  }
  return PyLong_FromLong(SettingGetType(index));
}

static PyObject* CmdGetSettingOfType(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1, int2, int3;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Oisii", &self, &int1, &str1, &int2, &int3);
  APIEnterBlocked(G);
  // may set the Python error indicator and return nullptr
  auto result = ExecutiveGetSettingOfType(G, int1, str1, int2, int3);
  APIExitBlocked(G);
  return result;
}

static PyObject* CmdUnset(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int index;
  char* str;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oisiii", &self, &index, &str, &state, &quiet,
                 &updates);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveUnsetSetting(G, index, str, state, quiet, updates);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdUnsetBond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oissiii", &self, &index, &str3, &str4, &state,
                 &quiet, &updates);
  API_ASSERT(APIEnterNotModal(G));
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    ExecutiveUnsetBondSetting(G, index, tmpsele1->getName(),
                              tmpsele2->getName(), state, quiet, updates);
    return {};
  }();
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdSetBond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int index;
  PyObject* value;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "OiOssiii", &self, &index, &value, &str3, &str4,
                 &state, &quiet, &updates);
  API_ASSERT(APIEnterNotModal(G));
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    if (!ExecutiveSetBondSetting(G, index, value, tmpsele1->getName(),
                                 tmpsele2->getName(), state, quiet, updates)) {
      return pymol::Error();
    }
    return {};
  }();
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetBond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int index;
  char *str3, *str4;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "Oissiii", &self, &index, &str3, &str4, &state,
                 &quiet, &updates);
  APIEnterBlocked(G);
  // TODO move selection handling to Executive
  auto res = [&]() -> pymol::Result<PyObject*> {
    auto tmpsele1 = SelectorTmp::make(G, str3);
    p_return_if_error(tmpsele1);
    auto tmpsele2 = SelectorTmp::make(G, str4);
    p_return_if_error(tmpsele2);
    return ExecutiveGetBondSetting(G, index,
                                   (/* TODO */ char*)tmpsele1->getName(),
                                   tmpsele2->getName(), state, quiet, updates);
  }();
  APIExitBlocked(G);
  return APIResult(G, res);
}

const PyMethodDef kSettingsMethods[] = {
    {"feedback", CmdFeedback, METH_VARARGS},
    {"get_setting_of_type", CmdGetSettingOfType, METH_VARARGS},
    {"get_setting_type", CmdGetSettingType, METH_VARARGS},
    {"get_setting_updates", CmdGetSettingUpdates, METH_VARARGS},
    {"get_setting_indices", CmdGetSettingIndices, METH_VARARGS},
    {"set_bond", CmdSetBond, METH_VARARGS},
    {"get_bond", CmdGetBond, METH_VARARGS},
    {"set_feedback", CmdSetFeedbackMask, METH_VARARGS},
    {"unset", CmdUnset, METH_VARARGS},
    {"unset_bond", CmdUnsetBond, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
