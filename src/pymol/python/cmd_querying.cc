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

// Python bindings for queries, properties and selections (pymol/querying.py,
// properties.py, selecting.py, monitoring.py).

#include "pymol/base/err.h"
#include "pymol/executive/executive.h"
#include "pymol/executive/executive_python.h"
#include "pymol/ray/scene_ray.h"
#include "pymol/scene/movie.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/executive/editor.h"
#include "pymol/ui/seq.h"
#include "pymol/base/feedback.h"
#include "pymol/io/cif_file.h"
#include "pymol/io/molecule_exporter.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdGetModalDraw(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int status = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    APIEnterBlocked(G);
    status = PyMOL_GetModalDraw(G->PyMOL);
    APIExitBlocked(G);
  }
  return APIResultCode(status);
}

static PyObject* CmdGetMovieLocked(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  return APIResultCode(MovieLocked(G));
}

static PyObject* CmdGetRawAlignment(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* name;
  int active_only;
  int state = 0;
  API_SETUP_ARGS(G, self, args, "Osi|i", &self, &name, &active_only, &state);
  pymol::Result<PyObject*> result;
  APIEnterBlocked(G);
  result = ExecutiveGetRawAlignment(G, name, active_only, state);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdGetVolumeField(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int state = 0;
  int ok = false;
  char* objName;
  short copy = 1;
  ok = PyArg_ParseTuple(args, "Os|ih", &self, &objName, &state, &copy);

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterBlockedNotModal(G))) {
    CField* field = ExecutiveGetVolumeField(G, objName, state);
    if (field) {
      result = FieldAsNumPyArray(field, copy);
    }
    APIExitBlocked(G);
  }

  if (!result) {
    return APIFailure();
  } else
    return result;
}

static PyObject* CmdGetVolumeHistogram(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* objName;
  float min_val = 0.f, max_val = 0.f;
  int n_points = 64;
  API_SETUP_ARGS(G, self, args, "Os|i(ff)", &self, &objName, &n_points,
                 &min_val, &max_val);
  API_ASSERT(APIEnterBlockedNotModal(G));

  auto res = ExecutiveGetHistogram(G, objName, n_points, min_val, max_val);

  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject* CmdGetSymmetry(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  float a, b, c, alpha, beta, gamma;
  int state;
  WordType sg;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);

  APIEnter(G);
  auto res = ExecutiveGetSymmetry(G, str1, state, &a, &b, &c, &alpha, &beta,
                                  &gamma, sg);
  APIExit(G);

  if (!res) {
    return APIFailure(G, res.error());
  }

  if (!res.result()) {
    return APIAutoNone(nullptr);
  }

  return Py_BuildValue("[fff fff s]", a, b, c, alpha, beta, gamma, sg);
}

static PyObject* CmdGetColor(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int mode;
  int a, nc, nvc;
  const float* rgb;
  int index;
  PyObject* result = nullptr;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &mode);
  APIEnterBlocked(G);
  {
    switch (mode) {
      case 0: /* by name or index, return floats */
        index = ColorGetIndex(G, name);
        if (index >= 0) {
          rgb = ColorGet(G, index);
          result = Py_BuildValue("fff", rgb[0], rgb[1], rgb[2]);
        }
        break;
      case 1: /* get color names with NO NUMBERS in their names */
        nc = ColorGetNColor(G);
        nvc = 0;
        for (a = 0; a < nc; a++) {
          if (ColorGetStatus(G, a) == 1) nvc++;
        }
        result = PyList_New(nvc);
        nvc = 0;
        for (a = 0; a < nc; a++) {
          if (ColorGetStatus(G, a) == 1) {
            const char* color_name = ColorGetName(G, a);
            if (!color_name) {
              color_name = "";
            }
            PyObject* tup = Py_BuildValue("si", color_name, a);
            PyList_SetItem(result, nvc++, tup);
          }
        }
        break;
      case 2: /* get all colors */
        nc = ColorGetNColor(G);
        nvc = 0;
        for (a = 0; a < nc; a++) {
          if (ColorGetStatus(G, a) != 0) nvc++;
        }
        result = PyList_New(nvc);
        nvc = 0;
        for (a = 0; a < nc; a++) {
          if (ColorGetStatus(G, a)) {
            const char* color_name = ColorGetName(G, a);
            if (!color_name) {
              color_name = "";
            }
            PyObject* tup = Py_BuildValue("si", color_name, a);
            PyList_SetItem(result, nvc++, tup);
          }
        }
        break;
      case 3: /* get a single color index */
        result = PyInt_FromLong(ColorGetIndex(G, name));
        break;
      case 4: /* by name or index, return floats including negative R for
                 special colors */
        index = ColorGetIndex(G, name);
        rgb = ColorGetSpecial(G, index);
        result = Py_BuildValue("fff", rgb[0], rgb[1], rgb[2]);
        break;
    }
    APIExitBlocked(G);
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetChains(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;

  char* str1;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &int1);
  APIEnter(G);
  auto res = ExecutiveGetChains(G, str1, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetRenderer(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *vendor = nullptr, *renderer = nullptr, *version = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneGetCardInfo(G, &vendor, &renderer, &version);
  APIExit(G);
  return Py_BuildValue("(sss)", vendor, renderer, version);
}

static PyObject* CmdGetPosition(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float v[3] = {0.0F, 0.0F, 0.0F};
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneGetCenter(G, v);
  APIExit(G);
  return PConvFloatArrayToPyList(v, 3);
}

static PyObject* CmdGetPhiPsi(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  PyObject* result = Py_None;
  pymol::vla<ObjectMolecule*> oVLA;
  pymol::vla<int> iVLA;
  pymol::vla<float> pVLA, sVLA;
  int l = 0;
  int a;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);
  APIEnter(G);
  {
    l = ExecutivePhiPsi(G, str1, &oVLA, &iVLA, &pVLA, &sVLA, state);
    APIExit(G);
    if (iVLA) {
      result = PyDict_New();
      for (a = 0; a < l; a++) {
        auto key =
            Py_BuildValue("si", oVLA[a]->Name, iVLA[a] + 1 /* 1-based */);
        auto value = Py_BuildValue("ff", pVLA[a], sVLA[a]);
        PyDict_SetItem(result, key, value);
        Py_DECREF(key);
        Py_DECREF(value);
      }
    } else {
      result = PyDict_New();
    }
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetCoordsAsNumPy(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state = 0;
  OrthoLineType s1;
  PyObject* result = nullptr;

  API_SETUP_ARGS(G, self, args, "Os|i", &self, &str1, &state);
  API_ASSERT(str1[0]);
  APIEnterBlocked(G);

  if (SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if (sele1 >= 0) {
      result = SelectorGetCoordsAsNumPy(G, sele1, state);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject* CmdGetCoordSetAsNumPy(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  CoordSet* cs;
  int state = 0;
  char* name;
  short copy = 1;

  API_SETUP_ARGS(G, self, args, "Os|ih", &self, &name, &state, &copy);

  ok_assert(2, name[0]);
  ok_assert(2, state >= 0);

  APIEnterBlocked(G);

  ok_assert(1, cs = ExecutiveGetCoordSet(G, name, state));
  result = CoordSetAsNumPyArray(cs, copy);

ok_except1:
  APIExitBlocked(G);
ok_except2:
  return (APIAutoNone(result));
}

static PyObject* CmdGetViewPort(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int width, height;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneGetWidthHeight(G, &width, &height);
  APIExit(G);
  return Py_BuildValue("ii", width, height);
}

static PyObject* CmdGetTitle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1;
  PyObject* result = Py_None;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &int1);
  APIEnter(G);
  const char* str2 = ExecutiveGetTitle(G, str1, int1);
  APIExit(G);
  if (str2) result = PyString_FromString(str2);
  return (APIAutoNone(result));
}

static PyObject* CmdGetArea(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, int2;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  APIEnter(G);
  auto res = ExecutiveGetArea(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetObjectSettings(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  pymol::CObject* obj = nullptr;
  PyObject* result = nullptr;
  const char* oname;
  int state = -1;

  if (!PyArg_ParseTuple(args, "Os|i", &self, &oname, &state)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnterBlocked(G);

  obj = ExecutiveFindObjectByName(G, oname);

  if (!obj) {
    ErrMessage(G, "GetObjectSettings", "named object not found.");
  } else {
    auto handle = obj->getSettingHandle(-1);

    if (state != -1) {
      auto handle_state = obj->getSettingHandle(state);

      // only accept handle if different from object-level settings
      handle = (handle_state == handle) ? nullptr : handle_state;
    }

    if (handle) {
      result = SettingAsPyList(handle->get(), true);
    }
  }

  APIExitBlocked(G);
ok_except1:
  return APIAutoNone(result);
}

static PyObject* CmdGetUnusedName(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* prefix = nullptr;
  int alwaysnumber = false;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osi", &self, &prefix, &alwaysnumber);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterBlockedNotModal(G))) {
    auto result =
        PConvToPyObject(ExecutiveGetUnusedName(G, prefix, alwaysnumber));
    APIExitBlocked(G);
    return result;
  } else {
    return APIResultOk(ok);
  }
}

static PyObject* CmdGetDragObjectName(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = Py_None;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    char* name = nullptr;
    APIEnter(G);
    {
      pymol::CObject* obj = EditorDragObject(G);
      if (obj) name = obj->Name;
    }
    APIExit(G);
    if (name)
      result = PyString_FromString(name);
    else
      result = PyString_FromString("");
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetLegalName(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  WordType name;
  char* str0;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str0);
  APIEnter(G);
  UtilNCopy(name, str0, sizeof(WordType));
  ObjectMakeValidName(G, name, true /* quiet */);
  APIExit(G);
  return PyString_FromString(name);
}

static PyObject* CmdPop(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int quiet;
  int result = 0;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossi", &self, &str1, &str2, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    result = ExecutivePop(G, str1, str2, quiet);
    APIExit(G);
  } else
    result = -1;
  return (APIResultCode(result));
}

static PyObject* CmdGetVRML(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ver;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &ver);
  {
    char* vla = nullptr;
    API_ASSERT(APIEnterNotModal(G));
    SceneRay(G, 0, 0, (ver == 1) ? 6 : 4, /* VRML1 or 2? */
             nullptr, &vla, 0.0F, 0.0F, false, nullptr, false, -1);
    APIExit(G);
    if (vla) {
      result = Py_BuildValue("s", vla);
    }
    VLAFreeP(vla);
  }
  return (APIAutoNone(result));
}

/**
 * Return a COLLADA string or None on failure
 */
static PyObject* CmdGetCOLLADA(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ver;
  char* vla = nullptr;

  API_SETUP_ARGS(G, self, args, "Oi", &self, &ver);
  API_ASSERT(APIEnterNotModal(G));

  SceneRay(G, 0, 0, 8, /* mode 8 = COLLADA */
           nullptr, &vla, 0.0F, 0.0F, false, nullptr, false, -1);
  APIExit(G);

  if (vla && vla[0]) {
    result = Py_BuildValue("s", vla);
  }

  VLAFreeP(vla);

  return APIAutoNone(result);
}

static PyObject* CmdGetIdtf(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *node = nullptr, *rsrc = nullptr;
    API_ASSERT(APIEnterNotModal(G));
    SceneRay(G, 0, 0, cSceneRay_MODE_IDTF, &node, &rsrc, 0.0F, 0.0F, false,
             nullptr, false, -1);
    APIExit(G);
    if (node && rsrc) {
      result = Py_BuildValue("(ss)", node, rsrc);
    }
    VLAFreeP(node);
    VLAFreeP(rsrc);
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetPovRay(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *header = nullptr, *geom = nullptr;
    API_ASSERT(APIEnterNotModal(G));
    SceneRay(G, 0, 0, 1, &header, &geom, 0.0F, 0.0F, false, nullptr, false, -1);
    APIExit(G);
    if (header && geom) {
      result = Py_BuildValue("(ss)", header, geom);
    }
    VLAFreeP(header);
    VLAFreeP(geom);
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetMtlObj(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    char *obj = nullptr, *mtl = nullptr;
    API_ASSERT(APIEnterNotModal(G));
    SceneRay(G, 0, 0, 5, &obj, &mtl, 0.0F, 0.0F, false, nullptr, false, -1);
    APIExit(G);
    if (obj && mtl) {
      result = Py_BuildValue("(ss)", mtl, obj);
    }
    VLAFreeP(obj);
    VLAFreeP(mtl);
  }
  return (APIAutoNone(result));
}

static PyObject* CmdOverlap(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state1, state2;
  float overlap = -1.0;
  float adjust;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossiif", &self, &str1, &str2, &state1, &state2,
                        &adjust);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    overlap = ExecutiveOverlap(G, str1, state1, str2, state2, adjust);
    APIExit(G);
  }
  if (overlap < 0.f) {
    return APIFailure(G);
  }
  return (Py_BuildValue("f", overlap));
}

static PyObject* CmdDist(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *str1, *str2;
  float cutoff;
  int labels, quiet;
  int mode, reset, state, zoom;
  int state1, state2;
  API_SETUP_ARGS(G, self, args, "Osssifiiiiiii", &self, &name, &str1, &str2,
                 &mode, &cutoff, &labels, &quiet, &reset, &state, &zoom,
                 &state1, &state2);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveDistance(G, name, str1, str2, mode, cutoff, labels, quiet,
                               reset, state, zoom, state1, state2);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdAngle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *str1, *str2, *str3;
  int labels, quiet;
  int mode;
  int reset, zoom;
  int state;
  int state1, state2, state3;
  API_SETUP_ARGS(G, self, args, "Ossssiiiiiiiii", &self, &name, &str1, &str2,
                 &str3, &mode, &labels, &reset, &zoom, &quiet, &state, &state1,
                 &state2, &state3);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveAngle(G, name, str1, str2, str3, mode, labels, reset,
                            zoom, quiet, state, state1, state2, state3);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdDihedral(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *str1, *str2, *str3, *str4;
  int labels, quiet;
  int mode;
  int reset, zoom;
  int state;
  API_SETUP_ARGS(G, self, args, "Osssssiiiiii", &self, &name, &str1, &str2,
                 &str3, &str4, &mode, &labels, &reset, &zoom, &quiet, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveDihedral(G, name, str1, str2, str3, str4, mode, labels,
                               reset, zoom, quiet, state);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdSelectList(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *sele_name;
  int quiet;
  int mode;
  int state;
  PyObject* list;
  API_SETUP_ARGS(G, self, args, "OssO!iii", &self, &sele_name, &str1,
                 &PyList_Type, &list, &state, &mode, &quiet);
  std::vector<int> int_array;
  API_ASSERT(PConvFromPyObject(G, list, int_array));
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSelectList(G, sele_name, str1, int_array.data(),
                                    int_array.size(), state, mode, quiet);
  SceneInvalidate(G);
  SeqDirty(G);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdReady(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    return (APIResultCode(G->Ready));
  } else {
    return (APIResultCode(0));
  }
}

static PyObject* CmdCountFrames(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  SceneCountFrames(G);
  int result = SceneGetNFrame(G, nullptr);
  APIExit(G);
  return (APIResultCode(result));
}

static PyObject* CmdGetMovieLength(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int result = 0;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  result = MovieGetLength(G);
  APIExit(G);
  return (APIResultCode(result));
}

static PyObject* CmdIndex(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int mode;
  pymol::vla<ObjectMolecule*> oVLA;
  pymol::vla<int> iVLA;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &mode);
  APIEnter(G);
  int l = ExecutiveIndex(G, str1, mode, &iVLA, &oVLA);
  APIExit(G);

  if (l == -1) {
    return APIFailure(G, "invalid selection");
  }

  PyObject* result = PyList_New(l);
  for (int a = 0; a < l; ++a) {
    auto tuple = Py_BuildValue("si", oVLA[a]->Name, iVLA[a] + 1 /* 1-based */);
    PyList_SetItem(result, a, tuple);
  }
  return result;
}

static PyObject* CmdFindPairs(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state1, state2;
  float cutoff;
  float angle;
  int mode;
  pymol::vla<int> iVLA;
  pymol::vla<ObjectMolecule*> oVLA;

  API_SETUP_ARGS(G, self, args, "Ossiiiff", &self, &str1, &str2, &state1,
                 &state2, &mode, &cutoff, &angle);
  APIEnter(G);
  auto res = ExecutivePairIndices(G, str1, str2, state1, state2, mode, cutoff,
                                  angle, &iVLA, &oVLA);
  APIExit(G);

  if (!res) {
    return APIFailure(G, res.error());
  }

  int l = res.result();
  PyObject* result = PyList_New(l);
  for (int a = 0; a < l; ++a) {
    auto tuple =
        Py_BuildValue("(si)(si)",  //
                      oVLA[a * 2]->Name, iVLA[a * 2] + 1 /* 1-based */,
                      oVLA[a * 2 + 1]->Name, iVLA[a * 2 + 1] + 1 /* 1-based */);
    PyList_SetItem(result, a, tuple);
  }
  return result;
}

static PyObject* CmdGetModel(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  char* ref_object;
  int ref_state;
  OrthoLineType s1;
  PyObject* result = nullptr;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Osisi", &self, &str1, &state, &ref_object,
                 &ref_state);
  if (!ref_object[0]) ref_object = nullptr;
  APIEnterBlocked(G);
  ok = (SelectorGetTmp(G, str1, s1) >= 0);
  if (ok)
    result = ExecutiveSeleToChemPyModel(G, s1, state, ref_object, ref_state);
  SelectorFreeTmp(G, s1);
  APIExitBlocked(G);
  if (!result) return APIFailure(G);
  return result;
}

static PyObject* CmdGetBonds(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  char* sele;
  int state;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &sele, &state);
  APIEnter(G);

  result = MoleculeExporterGetPyBonds(G, sele, state);
  APIExit(G);

  return APIAutoNone(result);
}

static PyObject* CmdGetAtomCoords(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &state, &quiet);
  APIEnter(G);
  auto result = ExecutiveGetAtomVertex(G, str1, state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdGetDistance(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int int1;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &int1);
  APIEnter(G);
  auto res = ExecutiveGetDistance(G, str1, str2, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetAngle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2, *str3;
  int int1;
  API_SETUP_ARGS(G, self, args, "Osssi", &self, &str1, &str2, &str3, &int1);
  APIEnter(G);
  auto res = ExecutiveGetAngle(G, str1, str2, str3, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetDihe(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2, *str3, *str4;
  int int1;
  API_SETUP_ARGS(G, self, args, "Ossssi", &self, &str1, &str2, &str3, &str4,
                 &int1);
  APIEnter(G);
  auto res = ExecutiveGetDihe(G, str1, str2, str3, str4, int1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetProgress(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int reset; /* TO DO */
  ok = PyArg_ParseTuple(args, "Oi", &self, &reset);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    if (G->Ready && (!SettingGetGlobal_b(G, cSetting_sculpting))) {
      /* assumes status is already locked */

      float result = -1.0F;
      float value = 0.0F, range = 1.0F;
      int ok = false;
      int int1;
      int offset;
      int progress[PYMOL_PROGRESS_SIZE];

      ok = PyArg_ParseTuple(args, "Oi", &self, &int1);
      if (ok) {
        if (PyMOL_GetBusy(G->PyMOL, false)) {
          PyMOL_GetProgress(G->PyMOL, progress, false);

          for (offset = PYMOL_PROGRESS_FAST; offset >= PYMOL_PROGRESS_SLOW;
               offset -= 2) {
            if (progress[offset + 1]) {
              float old_value = value;
              float old_range = range;

              range = (float)(progress[offset + 1]);
              value = (float)(progress[offset]);

              value += (1.0F / range) * (old_value / old_range);

              result = value / range;
            }
          }
        }
      }
      return (PyFloat_FromDouble((double)result));
    }
  }
  return (PyFloat_FromDouble(-1.0));
}

static PyObject* CmdGetMinMax(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float mn[3], mx[3];
  char* str1;
  int state;
  OrthoLineType s1;
  int flag;

  int ok = false;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &state);
  APIEnter(G);
  {
    ok = (SelectorGetTmp2(G, str1, s1) >= 0);
    flag = ExecutiveGetExtent(G, s1, mn, mx, true, state, false);
    SelectorFreeTmp(G, s1);
    APIExit(G);
    if (flag)
      return Py_BuildValue("[[fff],[fff]]", mn[0], mn[1], mn[2], mx[0], mx[1],
                           mx[2]);
  }
  return Py_BuildValue("[[fff],[fff]]", -0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
}

static PyObject* CmdGetObjectMatrix(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  char* name;
  double* history = nullptr;
  int ok = false;
  int found;
  int state;
  int incl_ttt = true;
  ok = PyArg_ParseTuple(args, "Osi|i", &self, &name, &state, &incl_ttt);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    APIEnter(G);
    found = ExecutiveGetObjectMatrix(G, name, state, &history, incl_ttt);
    APIExit(G);
    if (found) {
      if (history)
        result = Py_BuildValue("dddddddddddddddd", history[0], history[1],
                               history[2], history[3], history[4], history[5],
                               history[6], history[7], history[8], history[9],
                               history[10], history[11], history[12],
                               history[13], history[14], history[15]);
      else
        result =
            Py_BuildValue("dddddddddddddddd", 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    }
  }
  return APIAutoNone(result);
}

static PyObject* CmdGetObjectTTT(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  const char* name;
  int state, quiet;
  const float* ttt = nullptr;

  if (!PyArg_ParseTuple(args, "Osii", &self, &name, &state, &quiet)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnterBlocked(G);

  ExecutiveGetObjectTTT(G, name, &ttt, state, quiet);
  if (ttt) result = PConvFloatArrayToPyList(ttt, 16);

  APIExitBlocked(G);
ok_except1:
  return APIAutoNone(result);
}

static PyObject* CmdOnOff(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int state;
  int parents = 0;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &name, &state, &parents);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetObjVisib(G, name, state, parents);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdSelect(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *sname, *sele;
  int quiet;
  int state = 0;
  char* domain;
  int enable = -1;
  int merge = 0;
  API_SETUP_ARGS(G, self, args, "Ossiis|ii", &self, &sname, &sele, &quiet,
                 &state, &domain, &enable, &merge);
  API_ASSERT(APIEnterNotModal(G));

  auto seleargs = ExecutiveSelectPrepareArgs(G, sname, sele);
  auto res = ExecutiveSelect(G, seleargs, enable, quiet, merge, state, domain);

  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetObjectColorIndex(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
  auto result = ExecutiveGetObjectColorIndex(G, str1);
  APIExit(G);
  return (APIResultCode(result));
}

static PyObject* CmdSetProperty(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* value;
  char *sname, *propname;
  int proptype;
  int state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "OsOsiii", &self, &propname, &value, &sname,
                 &proptype, &state, &quiet);
  APIEnterBlocked(G);
  pymol::Result<bool> result;
  {
    auto r = ExecutiveSetPropertyForObject(G, propname, value, sname, state,
                                           proptype, quiet);
    if (r)
      result = true;
    else
      result = pymol::make_error(r.error().what());
  }
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdSetAtomProperty(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* value;
  char* sele;
  OrthoLineType s1;
  char* propname;
  int proptype;
  int state;
  int quiet;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "OsOsiii", &self, &propname, &value, &sele,
                 &proptype, &state, &quiet);

  pymol::Result<> result;
  {
    APIEnterBlocked(G);
    ok = (SelectorGetTmp(G, sele, s1) >= 0);
    if (ok) {
      ok = ExecutiveSetAtomPropertyForSelection(G, propname, value, s1, state,
                                                proptype, quiet);
      SelectorFreeTmp(G, s1);
    }
    if (!ok) {
      result = pymol::Error::QUIET;
    }
    APIExitBlocked(G);
  }
  return APIResult(G, result);
}

static PyObject* CmdGetProperty(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *sname, *propname;
  int state, quiet;
  PyObject* result = Py_None;
  int ok = false;
  ok =
      PyArg_ParseTuple(args, "Ozsii", &self, &propname, &sname, &state, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    APIEnterBlocked(G);
    result = ExecutiveGetPropertyForObject(G, propname, sname, state, quiet);
    APIExitBlocked(G);
  }
  return APIAutoNone(result);
}

static PyObject* CmdCountDiscrete(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  OrthoLineType s1;
  pymol::vla<ObjectMolecule*> list;
  int discrete = 0;

  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnterBlocked(G);
  ok_assert(2, SelectorGetTmp(G, str1, s1) >= 0);

  if ((list = ExecutiveGetObjectMoleculeVLA(G, s1))) {
    unsigned int i, size = VLAGetSize(list);
    for (i = 0; i < size; i++)
      if (list[i]->DiscreteFlag) discrete++;
  }

  SelectorFreeTmp(G, s1);
ok_except2:
  APIExitBlocked(G);
  return Py_BuildValue("i", discrete);
}

/**
 * Experimental - SUBJECT TO CHANGE
 */
static PyObject* CmdCifGetArray(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char *name, *key, *dtype = "";
  ObjectMolecule* obj;
  PyObject* ret = nullptr;

  ok_assert(1, PyArg_ParseTuple(args, "Oss|s", &self, &name, &key, &dtype));
  API_SETUP_PYMOL_GLOBALS;
  ok_assert(1, G);
  APIEnterBlocked(G);

  obj = ExecutiveFindObjectMoleculeByName(G, name);

  if (!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
    " Executive-Error: object '%s' not found.\n", name ENDFB(G);
  } else if (!obj->m_cifdata) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    " Executive-Warning: no cif data for object '%s'\n"
    " ! The 'cif_keepinmemory' setting needs to be set prior to loading a cif "
    "file.\n",
        name ENDFB(G);
  } else {
    auto* arr = obj->m_cifdata->get_arr(key);
    if (!arr) {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Details: key '%s' not in cif data for object '%s'.\n", key,
          name ENDFB(G);
    } else {
      switch (dtype[0]) {
        case 'i':
          ret = PConvToPyObject(arr->to_vector<int>());
          break;
        case 'f':
          ret = PConvToPyObject(arr->to_vector<double>());
          break;
        default:
          ret = PConvToPyObject(arr->to_vector<const char*>());
          break;
      }
    }
  }

  APIExitBlocked(G);
  return APIAutoNone(ret);
ok_except1:
  API_HANDLE_ERROR;
  return APIAutoNone(nullptr);
}

const PyMethodDef kQueryingMethods[] = {
    {"angle", CmdAngle, METH_VARARGS},
    {"cif_get_array", CmdCifGetArray, METH_VARARGS},
    {"count_frames", CmdCountFrames, METH_VARARGS},
    {"count_discrete", CmdCountDiscrete, METH_VARARGS},
    {"dihedral", CmdDihedral, METH_VARARGS},
    {"dist", CmdDist, METH_VARARGS},
    {"find_pairs", CmdFindPairs, METH_VARARGS},
    {"get_angle", CmdGetAngle, METH_VARARGS},
    {"get_area", CmdGetArea, METH_VARARGS},
    {"get_atom_coords", CmdGetAtomCoords, METH_VARARGS},
    {"get_chains", CmdGetChains, METH_VARARGS},
    {"get_collada", CmdGetCOLLADA, METH_VARARGS},
    {"get_color", CmdGetColor, METH_VARARGS},
    {"get_coords", CmdGetCoordsAsNumPy, METH_VARARGS},
    {"get_coordset", CmdGetCoordSetAsNumPy, METH_VARARGS},
    {"get_distance", CmdGetDistance, METH_VARARGS},
    {"get_dihe", CmdGetDihe, METH_VARARGS},
    {"get_drag_object_name", CmdGetDragObjectName, METH_VARARGS},
    {"get_idtf", CmdGetIdtf, METH_VARARGS},
    {"get_legal_name", CmdGetLegalName, METH_VARARGS},
    {"get_min_max", CmdGetMinMax, METH_VARARGS},
    {"get_mtl_obj", CmdGetMtlObj, METH_VARARGS},
    {"get_model", CmdGetModel, METH_VARARGS},
    {"get_property", CmdGetProperty, METH_VARARGS},
    {"get_bonds", CmdGetBonds, METH_VARARGS},
    {"get_modal_draw", CmdGetModalDraw, METH_VARARGS},
    {"get_movie_length", CmdGetMovieLength, METH_VARARGS},
    {"get_movie_locked", CmdGetMovieLocked, METH_VARARGS},
    {"get_object_color_index", CmdGetObjectColorIndex, METH_VARARGS},
    {"get_object_matrix", CmdGetObjectMatrix, METH_VARARGS},
    {"get_object_ttt", CmdGetObjectTTT, METH_VARARGS},
    {"get_object_settings", CmdGetObjectSettings, METH_VARARGS},
    {"get_position", CmdGetPosition, METH_VARARGS},
    {"get_povray", CmdGetPovRay, METH_VARARGS},
    {"get_progress", CmdGetProgress, METH_VARARGS},
    {"get_phipsi", CmdGetPhiPsi, METH_VARARGS},
    {"get_renderer", CmdGetRenderer, METH_VARARGS},
    {"get_raw_alignment", CmdGetRawAlignment, METH_VARARGS},
    {"get_symmetry", CmdGetSymmetry, METH_VARARGS},
    {"get_title", CmdGetTitle, METH_VARARGS},
    {"get_unused_name", CmdGetUnusedName, METH_VARARGS},
    {"get_viewport", CmdGetViewPort, METH_VARARGS},
    {"get_volume_field", CmdGetVolumeField, METH_VARARGS},
    {"get_volume_histogram", CmdGetVolumeHistogram, METH_VARARGS},
    {"get_vrml", CmdGetVRML, METH_VARARGS},
    {"index", CmdIndex, METH_VARARGS},
    {"onoff", CmdOnOff, METH_VARARGS},
    {"overlap", CmdOverlap, METH_VARARGS},
    {"pop", CmdPop, METH_VARARGS},
    {"ready", CmdReady, METH_VARARGS},
    {"select", CmdSelect, METH_VARARGS},
    {"select_list", CmdSelectList, METH_VARARGS},
    {"set_property", CmdSetProperty, METH_VARARGS},
    {"set_atom_property", CmdSetAtomProperty, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
