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

// Python bindings for the editing commands (pymol/editing.py).

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule3.h"
#include "pymol/executive/executive.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/executive/editor.h"
#include "pymol/chem/sculpt_cache.h"
#include "pymol/base/feedback.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdFixChemistry(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str2, *str3;
  int quiet;
  int invalidate;
  API_SETUP_ARGS(G, self, args, "Ossii", &self, &str2, &str3, &invalidate,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFixChemistry(G, str2, str3, invalidate, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdGetOrigin(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  float origin[3];
  char* object;
  ok = PyArg_ParseTuple(args, "Os", &self, &object);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterBlockedNotModal(G))) {
    if ((!object) || (!object[0])) {
      SceneOriginGet(G, origin);
    } else {
      pymol::CObject* obj = ExecutiveFindObjectByName(G, object);
      if (!obj) {
        ok = false;
      } else {
        if (obj->TTTFlag) {
          origin[0] = -obj->TTT[12];
          origin[1] = -obj->TTT[13];
          origin[2] = -obj->TTT[14];
        } else {
          SceneOriginGet(G, origin); /* otherwise, return scene origin */
        }
      }
    }
    APIExitBlocked(G);
  }
  if (ok) {
    return (Py_BuildValue("fff", origin[0], origin[1], origin[2]));
  } else {
    return APIFailure();
  }
}

static PyObject* CmdSetSymmetry(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state;
  float a, b, c, alpha, beta, gamma;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Osiffffffsi", &self, &str1, &state, &a, &b, &c,
                 &alpha, &beta, &gamma, &str2, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSetSymmetry(G, str1, state, a, b, c, alpha, beta, gamma,
                                  str2, quiet);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdSmooth(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sele;
  int cycles, window, first, last, ends, quiet;
  float cutoff = -1;
  int pbc = true;
  API_SETUP_ARGS(G, self, args, "Osiiiiii|fi", &self, &sele, &cycles, &window,
                 &first, &last, &ends, &quiet, &cutoff, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSmooth(G, sele, cycles, window, first, last, ends,
                                quiet, cutoff, pbc);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetName(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  API_SETUP_ARGS(G, self, args, "Oss", &self, &str1, &str2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetName(G, str1, str2);
  APIExit(G);
  if (result) {
    return APISuccess();
  }
  return APIFailure(G, result.error());
}

static PyObject* CmdSculptPurge(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    SculptCachePurge(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdSculptDeactivate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  char* str1;
  ok = PyArg_ParseTuple(args, "Os", &self, &str1);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSculptDeactivate(G, str1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdSculptActivate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int int1, int2, int3;
  char* str1;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &int1, &int2, &int3);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveSculptActivate(G, str1, int1, int2, int3);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdSculptIterate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int int1, int2;
  char* str1;
  float total_strain = 0.0F;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    total_strain = ExecutiveSculptIterate(G, str1, int1, int2);
    APIExit(G);
  }
  return PyFloat_FromDouble((double)total_strain);
}

static PyObject* CmdSetObjectTTT(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float ttt[16];
  int quiet;
  char* name;
  int state;
  API_SETUP_ARGS(G, self, args, "Os(ffffffffffffffff)ii", &self, &name, &ttt[0],
                 &ttt[1], &ttt[2], &ttt[3], &ttt[4], &ttt[5], &ttt[6], &ttt[7],
                 &ttt[8], &ttt[9], &ttt[10], &ttt[11], &ttt[12], &ttt[13],
                 &ttt[14], &ttt[15], &state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveSetObjectTTT(G, name, ttt, state, quiet,
                            SettingGetGlobal_i(G, cSetting_movie_auto_store));
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTranslateObjectTTT(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float mov[3];
  char* name;
  API_SETUP_ARGS(G, self, args, "Os(fff)", &self, &name, &mov[0], &mov[1],
                 &mov[2]);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTranslateObjectTTT(
      G, name, mov, SettingGetGlobal_i(G, cSetting_movie_auto_store), true);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdCombineObjectTTT(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  PyObject* m;
  float ttt[16];
  API_SETUP_ARGS(G, self, args, "OsO", &self, &name, &m);
  if (PConvPyListToFloatArrayInPlace(m, ttt, 16) <= 0) {
    return APIFailure(G, "Bad Matrix");
  }
  API_ASSERT(APIEnterNotModal(G))
  auto result = ExecutiveCombineObjectTTT(G, name, ttt, false, -1);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMapSetBorder(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  float level;
  int state;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osfi", &self, &name, &level, &state);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMapSetBorder(G, name, level, state);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdMapSet(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *operands;
  int target_state, source_state, operator_;
  int zoom, quiet;

  API_SETUP_ARGS(G, self, args, "Osisiiii", &self, &name, &operator_, &operands,
                 &target_state, &source_state, &zoom, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapSet(G, name, operator_, operands, target_state,
                                source_state, zoom, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMapTrim(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *sele;
  int map_state, sele_state;
  float buffer;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossfiii", &self, &name, &sele, &buffer,
                 &map_state, &sele_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveMapTrim(G, name, sele, buffer, map_state, sele_state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMapDouble(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int state;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapDouble(G, name, state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMapHalve(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int state;
  int smooth;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &name, &state, &smooth);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapHalve(G, name, state, smooth);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTranslateAtom(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state, log, mode;
  float v[3];
  API_SETUP_ARGS(G, self, args, "Osfffiii", &self, &str1, v, v + 1, v + 2,
                 &state, &mode, &log);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTranslateAtom(G, str1, v, state, mode, log);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMatrixCopy(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *source_name, *target_name;
  int source_mode, target_mode;
  int source_state, target_state, target_undo;
  int log;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossiiiiiii", &self, &source_name, &target_name,
                 &source_mode, &target_mode, &source_state, &target_state,
                 &target_undo, &log, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveMatrixCopy(G, source_name, target_name, source_mode, target_mode,
                      source_state, target_state, target_undo, log, quiet);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdResetMatrix(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int mode;
  int state;
  int log;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osiiii", &self, &name, &mode, &state, &log,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveResetMatrix(G, name, mode, state, log, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTransformObject(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *sele;
  int state, log;
  PyObject* m;
  float matrix[16];
  int homo;
  API_SETUP_ARGS(G, self, args, "OsiOisi", &self, &name, &state, &m, &log,
                 &sele, &homo);
  if (PConvPyListToFloatArrayInPlace(m, matrix, 16) <= 0) {
    return APIFailure(G, "Bad Matrix");
  }
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTransformObjectSelection(G, name, state, sele, log,
                                                  matrix, homo, true);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTransformSelection(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sele;
  int state, log;
  int homo;
  PyObject* m;
  float ttt[16];
  API_SETUP_ARGS(G, self, args, "OsiOii", &self, &sele, &state, &m, &log,
                 &homo);
  if (PConvPyListToFloatArrayInPlace(m, ttt, 16) <= 0) {
    return APIFailure(G, "Bad Matrix");
  }
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveTransformSelection(G, state, sele, log, ttt, homo);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdGetEditorScheme(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    result = EditorGetScheme(G);
  }
  return (APIResultCode(result));
}

static PyObject* CmdSetTitle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osis", &self, &str1, &int1, &str2);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetTitle(G, str1, int1, str2);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdPushUndo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str0;
  int state;
  OrthoLineType s0 = "";
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str0, &state);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    if (str0[0]) ok = (SelectorGetTmp(G, str0, s0) >= 0);
    if (ok) ok = ExecutiveSaveUndo(G, s0, state);
    if (s0[0]) SelectorFreeTmp(G, s0);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdInvert(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &int1);
  API_ASSERT(APIEnterNotModal(G));
  auto res = EditorInvert(G, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdTorsion(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float float1;
  API_SETUP_ARGS(G, self, args, "Of", &self, &float1);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorTorsion(G, float1);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdUndo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ExecutiveUndo(G, int1); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdProtect(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveProtect(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSymmetryCopy(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *source_name, *target_name;
  int source_state, target_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossiii", &self, &source_name, &target_name,
                 &source_state, &target_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSymmetryCopy(G, source_name, target_name, source_state,
                                   target_state, quiet);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdBond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int order, mode;
  int quiet;
  const char* symop = "";
  API_SETUP_ARGS(G, self, args, "Ossiii|s", &self, &str1, &str2, &order, &mode,
                 &quiet, &symop);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveBond(G, str1, str2, order, mode, quiet, symop);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAddBond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int atm1, atm2;
  int order;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &oname, &atm1, &atm2, &order);
  APIEnterBlocked(G);

  auto result = ExecutiveAddBondByIndices(G, oname, atm1, atm2, order);

  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdRebond(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int state;
  int pbc = 0;
  API_SETUP_ARGS(G, self, args, "Osi|i", &self, &oname, &state, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveRebond(G, oname, state, pbc);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdRevalence(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *sele1, *sele2, *source;
  int source_state, target_state, reset;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssiiii", &self, &sele1, &sele2, &source,
                 &target_state, &source_state, &reset, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRevalence(G, sele1, sele2, source, target_state,
                                   source_state, reset, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdPBCUnwrap(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;
  int bymol = true;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &oname, &bymol);
  API_ASSERT(APIEnterNotModal(G));

  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    APIExit(G);
    return APIFailure(G, "cannot find object");
  }

  ObjectMoleculePBCUnwrap(*obj, bymol);

  APIExit(G);
  return APISuccess();
}

static PyObject* CmdPBCWrap(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;
  PyObject* pycenter = nullptr;
  API_SETUP_ARGS(G, self, args, "OsO", &self, &oname, &pycenter);

  std::vector<float> center;
  if (pycenter != Py_None) {
    API_ASSERT(PConvFromPyObject(G, pycenter, center) && center.size() == 3);
  }

  API_ASSERT(APIEnterNotModal(G));

  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    APIExit(G);
    return APIFailure(G, "cannot find object");
  }

  ObjectMoleculePBCWrap(*obj, center.empty() ? nullptr : center.data());

  APIExit(G);
  return APISuccess();
}

static PyObject* CmdVdwFit(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state1, state2, quiet;
  float buffer;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osisifi", &self, &str1, &state1, &str2, &state2,
                        &buffer, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveVdwFit(G, str1, state1, str2, state2, buffer, quiet);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdAlter(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int read_only, quiet;
  PyObject* space;
  API_SETUP_ARGS(G, self, args, "OssiiO", &self, &str1, &str2, &read_only,
                 &quiet, &space);
  API_ASSERT(APIEnterBlockedNotModal(G));
  pymol::Result<int> result{-1};
  if (read_only) {
    result = ExecutiveIterate(G, str1, str2, read_only, quiet, space);
  } else {
    result = ExecutiveIterate(G, str1, str2, read_only, quiet, space);
  }
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdAlterList(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int quiet;
  PyObject* space;
  PyObject* list;
  API_SETUP_ARGS(G, self, args, "OsOiO", &self, &str1, &list, &quiet, &space);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveIterateList(G, str1, list, false, quiet, space);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdAlterState(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state, read_only, quiet;
  PyObject* obj;
  API_SETUP_ARGS(G, self, args, "OissiiO", &self, &state, &str1, &str2,
                 &read_only, &quiet, &obj);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result =
      ExecutiveIterateState(G, state, str1, str2, read_only, quiet, obj);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdUpdate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int int1, int2;
  int matchmaker, quiet;
  API_SETUP_ARGS(G, self, args, "Ossiiii", &self, &str1, &str2, &int1, &int2,
                 &matchmaker, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveUpdateCmd(G, str1, str2, int1, int2, matchmaker, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetDihe(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2, *str3, *str4;
  float float1;
  int int1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossssfii", &self, &str1, &str2, &str3, &str4,
                 &float1, &int1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveSetDihe(G, str1, str2, str3, str4, float1, int1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdFlag(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int flag;
  int action;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oisii", &self, &flag, &str1, &action, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFlag(G, flag, str1, action, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetStateOrder(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* oname;
  PyObject* order;
  pymol::CObject* obj = nullptr;
  int* int_array = nullptr;
  int ok = false;

  if (!PyArg_ParseTuple(args, "OsO", &self, &oname, &order)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  ok_assert(1, PyList_Check(order));

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G && APIEnterNotModal(G));

  obj = ExecutiveFindObjectByName(G, oname);
  if (!obj || obj->type != cObjectMolecule) {
    ErrMessage(G, "SetStateOrder", "named object molecule not found.");
    ok_raise(2);
  }

  if (PConvPyListToIntArray(order, &int_array)) {
    int len = PyList_Size(order);

    PBlock(G);
    ok = ObjectMoleculeSetStateOrder((ObjectMolecule*)obj, int_array, len);
    PUnblock(G);

    FreeP(int_array);
  } else {
    ErrMessage(G, "SetStateOrder", "not an integer list.");
    ok_raise(2);
  }

  APIExit(G);
  return APIResultOk(ok);
ok_except2:
  APIExit(G);
ok_except1:
  return APIFailure();
}

static PyObject* CmdSort(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  API_SETUP_ARGS(G, self, args, "Os", &self, &name);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSort(G, name);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAssignSS(PyObject* self, PyObject* args)

/* EXPERIMENTAL */
{
  PyMOLGlobals* G = nullptr;
  int state, quiet;
  char *str1, *str2;
  int preserve;
  API_SETUP_ARGS(G, self, args, "Osisii", &self, &str1, &state, &str2,
                 &preserve, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveAssignSS(G, str1, state, str2, preserve, nullptr, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdReference(PyObject* self, PyObject* args) {
  int ok = false;
  PyMOLGlobals* G = nullptr;
  OrthoLineType s1;
  int action;
  char* sele1;
  int state;
  int quiet;
  ok = PyArg_ParseTuple(args, "Oisii", &self, &action, &sele1, &state, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp(G, sele1, s1) >= 0);
    if (ok) ok = ExecutiveReference(G, action, s1, state, quiet);
    SelectorFreeTmp(G, s1);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdHAdd(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int quiet;
  int state;
  int legacy;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &str1, &quiet, &state, &legacy);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveAddHydrogens(G, str1, quiet, state, legacy);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetObjectColor(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *color;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &name, &color, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetObjectColor(G, name, color, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdRemove(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRemoveAtoms(G, str1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdRemovePicked(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &i1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorRemove(G, i1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdHFill(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorHFill(G, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdHFix(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int quiet;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  {
    SelectorTmp2 s1(G, str1);
    result = EditorHFix(G, s1.getName(), quiet);
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdCycleValence(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorCycleValence(G, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdReplace(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1, i2;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osiisi", &self, &str1, &i1, &i2, &str2,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorReplace(G, str1, i1, i2, str2, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetGeometry(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int geom, valence;
  char* sele;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &sele, &geom, &valence);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSetGeometry(G, sele, geom, valence);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAttach(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1, i2;
  char* str1;
  int quiet;
  char* name;
  API_SETUP_ARGS(G, self, args, "Osiis", &self, &str1, &i1, &i2, &name, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = EditorAttach(G, str1, i1, i2, name, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdFuse(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int mode;
  int recolor;
  int move_flag;
  API_SETUP_ARGS(G, self, args, "Ossiii", &self, &str1, &str2, &mode, &recolor,
                 &move_flag);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFuse(G, str1, str2, mode, recolor, move_flag);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdUnpick(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  EditorInactivate(G);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdEdit(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str0, *str1, *str2, *str3;
  int pkresi, pkbond;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossssiii", &self, &str0, &str1, &str2, &str3,
                 &pkresi, &pkbond, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  if (!str0[0]) {
    EditorInactivate(G);
  } else {
    result = EditorSelect(G, str0, str1, str2, str3, pkresi, pkbond, quiet);
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdDrag(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str0;
  OrthoLineType s0 = "";
  int ok = false;
  int quiet;
  int mode;
  ok = PyArg_ParseTuple(args, "Osii", &self, &str0, &quiet, &mode);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = (SelectorGetTmp2(G, str0, s0) >= 0);
    if (ok) {
      ok = ExecutiveSetDrag(G, s0, quiet, mode);
      SelectorFreeTmp(G, s0);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdRename(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveRenameObjectAtoms(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetDiscrete(PyObject* self, PyObject* args) {
  const char* name;
  int discrete;
  bool status = false;

  if (!PyArg_ParseTuple(args, "Osi", &self, &name, &discrete)) {
    API_HANDLE_ERROR;
  } else {
    PyMOLGlobals* G = nullptr;
    API_SETUP_PYMOL_GLOBALS;

    if (G && APIEnterBlockedNotModal(G)) {
      ObjectMolecule* obj = ExecutiveFindObjectMoleculeByName(G, name);

      if (!obj) {
        PRINTFB(G, FB_Executive, FB_Errors)
        " Executive-Error: object '%s' not found.\n", name ENDFB(G);
      } else {
        status = ObjectMoleculeSetDiscrete(G, obj, discrete);
      }

      APIExitBlocked(G);
    }
  }

  return APIResultOk(status);
}

static PyObject* CmdLookAt(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* mobileObjName;
  const char* targetObjName;
  API_SETUP_ARGS(G, self, args, "Oss", &self, &targetObjName, &mobileObjName);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveLookAt(G, targetObjName, mobileObjName);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdMoveOnCurve(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* mobileObjName;
  const char* curveObjName;
  float t;
  API_SETUP_ARGS(G, self, args, "Ossf", &self, &mobileObjName, &curveObjName,
                 &t);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveMoveOnCurve(G, mobileObjName, curveObjName, t);
  APIExitBlocked(G);
  return APIResult(G, result);
}

const PyMethodDef kEditingMethods[] = {
    {"alter", CmdAlter, METH_VARARGS},
    {"alter_list", CmdAlterList, METH_VARARGS},
    {"alter_state", CmdAlterState, METH_VARARGS},
    {"attach", CmdAttach, METH_VARARGS},
    {"bond", CmdBond, METH_VARARGS},
    {"add_bond", CmdAddBond, METH_VARARGS},
    {"rebond", CmdRebond, METH_VARARGS},
    {"combine_object_ttt", CmdCombineObjectTTT, METH_VARARGS},
    {"cycle_valence", CmdCycleValence, METH_VARARGS},
    {"drag", CmdDrag, METH_VARARGS},
    {"edit", CmdEdit, METH_VARARGS},
    {"torsion", CmdTorsion, METH_VARARGS},
    {"fix_chemistry", CmdFixChemistry, METH_VARARGS},
    {"flag", CmdFlag, METH_VARARGS},
    {"dss", CmdAssignSS, METH_VARARGS},
    {"fuse", CmdFuse, METH_VARARGS},
    {"get_editor_scheme", CmdGetEditorScheme, METH_VARARGS},
    {"get_origin", CmdGetOrigin, METH_VARARGS},
    {"h_add", CmdHAdd, METH_VARARGS},
    {"h_fill", CmdHFill, METH_VARARGS},
    {"h_fix", CmdHFix, METH_VARARGS},
    {"invert", CmdInvert, METH_VARARGS},
    {"look_at", CmdLookAt, METH_VARARGS},
    {"map_double", CmdMapDouble, METH_VARARGS},
    {"map_halve", CmdMapHalve, METH_VARARGS},
    {"map_set", CmdMapSet, METH_VARARGS},
    {"map_set_border", CmdMapSetBorder, METH_VARARGS},
    {"map_trim", CmdMapTrim, METH_VARARGS},
    {"move_on_curve", CmdMoveOnCurve, METH_VARARGS},
    {"protect", CmdProtect, METH_VARARGS},
    {"push_undo", CmdPushUndo, METH_VARARGS},
    {"pbc_unwrap", CmdPBCUnwrap, METH_VARARGS},
    {"pbc_wrap", CmdPBCWrap, METH_VARARGS},
    {"reference", CmdReference, METH_VARARGS},
    {"remove", CmdRemove, METH_VARARGS},
    {"remove_picked", CmdRemovePicked, METH_VARARGS},
    {"rename", CmdRename, METH_VARARGS},
    {"replace", CmdReplace, METH_VARARGS},
    {"reset_matrix", CmdResetMatrix, METH_VARARGS},
    {"revalence", CmdRevalence, METH_VARARGS},
    {"sculpt_deactivate", CmdSculptDeactivate, METH_VARARGS},
    {"sculpt_activate", CmdSculptActivate, METH_VARARGS},
    {"sculpt_iterate", CmdSculptIterate, METH_VARARGS},
    {"sculpt_purge", CmdSculptPurge, METH_VARARGS},
    {"set_dihe", CmdSetDihe, METH_VARARGS},
    {"set_discrete", CmdSetDiscrete, METH_VARARGS},
    {"set_name", CmdSetName, METH_VARARGS},
    {"set_geometry", CmdSetGeometry, METH_VARARGS},
    {"set_object_ttt", CmdSetObjectTTT, METH_VARARGS},
    {"set_object_color", CmdSetObjectColor, METH_VARARGS},
    {"set_state_order", CmdSetStateOrder, METH_VARARGS},
    {"set_symmetry", CmdSetSymmetry, METH_VARARGS},
    {"set_title", CmdSetTitle, METH_VARARGS},
    {"smooth", CmdSmooth, METH_VARARGS},
    {"sort", CmdSort, METH_VARARGS},
    {"symmetry_copy", CmdSymmetryCopy, METH_VARARGS},
    {"matrix_copy", CmdMatrixCopy, METH_VARARGS},
    {"transform_object", CmdTransformObject, METH_VARARGS},
    {"transform_selection", CmdTransformSelection, METH_VARARGS},
    {"translate_atom", CmdTranslateAtom, METH_VARARGS},
    {"translate_object_ttt", CmdTranslateObjectTTT, METH_VARARGS},
    {"vdw_fit", CmdVdwFit, METH_VARARGS},
    {"undo", CmdUndo, METH_VARARGS},
    {"unpick", CmdUnpick, METH_VARARGS},
    {"update", CmdUpdate, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
