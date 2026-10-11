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

/*
   NOTICE:

   Important thread safety tip:

   PM operations which will ultimately call GLUT can only be called by
   the main GLUT thread (with some exceptions, such as the simple
   drawing operations which seem to be thread safe).

   Thus, pm.py needs to guard against direct invocation of certain _cmd
   (API) calls from outside threads [Examples: _cmd.png(..),
   _cmd.mpng(..) ].  Instead, it needs to hand them over to the main
   thread by way of a cmd.mdo(...)  statement.

   Note that current, most glut operations have been pushed into the
   main event and redraw loops to avoid these kinds of problems - so
   I'm not sure how important this really is anymore.

*/

/* TODO: Put in some exception handling and reporting for the
 * python calls, especially, PyArg_ParseTuple()
 */

#define PY_SSIZE_T_CLEAN

#include "pymol/base/pymol_options.h"
#include "pymol/executive/executive.h"
#include "pymol/executive/executive_python.h"
#include "pymol/scene/movie.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/control.h"
#include "pymol/ui/wizard.h"
#include "pymol/base/feedback.h"

#include "pymol/scene/metal_pick.h"
#include "pymol/python/cmd_internal.h"

#define tmpSele "_tmp"
#define tmpSele1 "_tmp1"
#define tmpSele2 "_tmp2"

static int flush_count = 0;

/*
 * C-level tests (not built)
 */
static PyObject* CmdTest2(PyObject*, PyObject*) {
  PyErr_SetString(PyExc_NotImplementedError, "C-level tests are not built");
  return nullptr;
}

/**
 * Start a headless singleton instance in the current thread.
 *
 * Unlike when calling `pymol.finish_launching()`, there is no event loop,
 * so animations, continuous sculpting and modal draw are not supported.
 *
 * After calling this, SingletonPyMOLGlobals will be available.
 */
static void launch_library_singleton() {
  PyRun_SimpleString(
      "import pymol.invocation, pymol2\n"
      "pymol.invocation.parse_args(['pymol', '-cqk'])\n"
      "pymol2.SingletonPyMOL().start()");
}

/// Creating non-singleton instances disables auto-library mode
static bool auto_library_mode_disabled = false;

/**
 * Get the PyMOLGlobals pointer from the `self` object (_self._COb in Python).
 *
 * If _COb is None, launch a headless singleton ("library mode").
 */
PyMOLGlobals* _api_get_pymol_globals(PyObject* self) {
  if (self == Py_None) {
    if (auto_library_mode_disabled) {
      PyErr_SetString(PyExc_RuntimeError, "Missing PyMOL instance");
      return nullptr;
    }

    launch_library_singleton();
    return SingletonPyMOLGlobals;
  }

  if (self && PyCapsule_CheckExact(self)) {
    auto G_handle =
        reinterpret_cast<PyMOLGlobals**>(PyCapsule_GetPointer(self, nullptr));
    if (G_handle) {
      return *G_handle;
    }
  }

  return nullptr;
}

/* NOTE: the glut_thread_keep_out variable can only be changed by the thread
   holding the API lock, therefore this is safe even through increment
   isn't (necessarily) atomic. */

void APIEnter(PyMOLGlobals* G) { /* assumes API is locked */
  PRINTFD(G, FB_API)
  " APIEnter-DEBUG: as thread %ld.\n", PyThread_get_thread_ident() ENDFD;

  if (G->Terminating) { /* try to bail */

    exit(0);
  }

  if (!PIsGlutThread()) G->P_inst->glut_thread_keep_out++;
  PUnblock(G);
}

int APIEnterNotModal(PyMOLGlobals* G) { /* assumes API is locked */
  if (PyMOL_GetModalDraw(G->PyMOL)) {
    return false;
  } else {
    APIEnter(G);
    return true;
  }
}

void APIEnterBlocked(PyMOLGlobals* G) { /* assumes API is locked */

  PRINTFD(G, FB_API)
  " APIEnterBlocked-DEBUG: as thread %ld.\n", PyThread_get_thread_ident() ENDFD;

  if (G->Terminating) { /* try to bail */

    exit(0);
  }

  if (!PIsGlutThread()) G->P_inst->glut_thread_keep_out++;
}

int APIEnterBlockedNotModal(PyMOLGlobals* G) { /* assumes API is locked */
  if (PyMOL_GetModalDraw(G->PyMOL)) {
    return false;
  } else {
    APIEnterBlocked(G);
    return true;
  }
}

void APIExit(PyMOLGlobals* G) { /* assumes API is locked */
  PBlock(G);
  if (!PIsGlutThread()) G->P_inst->glut_thread_keep_out--;
  PRINTFD(G, FB_API)
  " APIExit-DEBUG: as thread %ld.\n", PyThread_get_thread_ident() ENDFD;
}

void APIExitBlocked(PyMOLGlobals* G) { /* assumes API is locked */

  if (!PIsGlutThread()) G->P_inst->glut_thread_keep_out--;
  PRINTFD(G, FB_API)
  " APIExitBlocked-DEBUG: as thread %ld.\n", PyThread_get_thread_ident() ENDFD;
}

/**
 * Return None (pymol.constants.DEFAULT_SUCCESS)
 */
PyObject* APISuccess(void) { /* success returns None */
  return PConvAutoNone(Py_None);
}

/**
 * Return -1 (pymol.constants.DEFAULT_ERROR)
 */
PyObject* APIFailure(void) { /* returns -1: a general unspecified
                              * error */
  return (Py_BuildValue("i", -1));
}

/**
 * Raise CmdException(msg).
 */
PyObject* APIFailure(PyMOLGlobals* G, const char* msg) {
  if (msg) {
    PyErr_SetString(P_CmdException, msg);
  } else {
    PyErr_SetNone(P_CmdException);
  }
  return nullptr;
}

/**
 * Raise `error` as an exception.
 */
PyObject* APIFailure(PyMOLGlobals* G, const pymol::Error& error) {
  if (PyErr_Occurred()) {
    return nullptr;
  }

  PyObject* exc_type;
  switch (error.code()) {
    case pymol::Error::QUIET:
      exc_type = P_QuietException;
      break;
    case pymol::Error::MEMORY:
      exc_type = PyExc_MemoryError;
      break;
    case pymol::Error::INCENTIVE_ONLY:
      exc_type = P_IncentiveOnlyException;
      break;
    default:
      exc_type = P_CmdException;
  }

  PyErr_SetString(exc_type, error.what().c_str());
  return nullptr;
}

PyObject* APIResultCode(int code) { /* innteger result code
                                     * (could be a value, a
                                     * count, or a boolean) */
  return (Py_BuildValue("i", code));
}

/**
 * If `ok` is true, return None (DEFAULT_SUCCESS).
 * Else return -1 (DEFAULT_ERROR).
 */
PyObject* APIResultOk(int ok) {
  if (ok)
    return APISuccess();
  else
    return APIFailure();
}

/**
 * If `ok` is true, return None (DEFAULT_SUCCESS).
 * Else raise CmdException().
 */
PyObject* APIResultOk(PyMOLGlobals* G, bool ok) {
  if (ok) return APISuccess();
  return APIFailure(G);
}

static PyObject* APIIncRef(PyObject* result) {
  Py_INCREF(result);
  return (result);
}

PyObject* APIAutoNone(PyObject* result) { /* automatically owned Py_None */
  if (result == Py_None)
    Py_INCREF(result);
  else if (result == nullptr) {
    result = Py_None;
    Py_INCREF(result);
  }
  return (result);
}

static PyObject* CmdRayAntiThread(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayAntiThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayAntiThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject* CmdRayHashThread(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayHashThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayHashThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject* CmdRayTraceThread(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CRayThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  RayTraceThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject* CmdCoordSetUpdateThread(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CCoordSetUpdateThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  CoordSetUpdateThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject* CmdObjectUpdateThread(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* py_thread_info;
  API_SETUP_ARGS(G, self, args, "OO", &self, &py_thread_info);

  auto thread_info = reinterpret_cast<CObjectUpdateThreadInfo*>(
      PyCapsule_GetPointer(py_thread_info, nullptr));
  API_ASSERT(thread_info);

  PUnblock(G);
  SceneObjectUpdateThread(thread_info);
  PBlock(G);

  return APISuccess();
}

static PyObject* CmdFakeDrag(PyObject* self, PyObject* args) {
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
    PyMOL_NeedFakeDrag(G->PyMOL);
  }
  return APISuccess();
}

static PyObject* CmdReinitialize(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int what;
  char* object;
  API_SETUP_ARGS(G, self, args, "Ois", &self, &what, &object);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveReinitialize(G, what, object);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetVersion(PyObject* self, PyObject* args) {
  return Py_BuildValue("(sdiisi)", _PyMOL_VERSION, _PyMOL_VERSION_double,
                       _PyMOL_VERSION_int,
#ifdef _PyMOL_BUILD_DATE
                       _PyMOL_BUILD_DATE, _PYMOL_BUILD_GIT_SHA, 0
#else
                       0, "", 0
#endif
  );
}

static PyObject* CmdGetCapabilities(PyObject*, PyObject*) {
  static PyObject* caps = nullptr;

  if (!caps) {
    caps = PySet_New(nullptr);
    // VMD molfile plugins
    PySet_Add(caps, PConvToPyObject("vmdplugins"));
    // object and atom level properties
    PySet_Add(caps, PConvToPyObject("properties"));
  }

  Py_INCREF(caps);
  return caps;
}

/**
 * Screen-space atom pick for the Metal backend (see
 * src/pymol/scene/metal_pick.h).
 *
 * `cam` is a 19-float list: 3x3 row-major rotation, camera position, rotation
 * origin, tan_half, aspect, clip_front, clip_back. `objects` is a list naming
 * the objects to search, in tie-break order; names that aren't molecular
 * objects are skipped. Returns None for empty space, else the tuple
 * (d2, object, chain, resi, resn, segi, name, ndc_x, ndc_y).
 */
static PyObject* CmdMetalPick(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject *objects_py, *cam_py;
  int state, rep_mask, guide_rep_mask;
  float ndc_x, ndc_y, max_ndc2, cluster_ndc2;
  API_SETUP_ARGS(G, self, args, "OOOiffffii", &self, &objects_py, &cam_py,
                 &state, &ndc_x, &ndc_y, &max_ndc2, &cluster_ndc2, &rep_mask,
                 &guide_rep_mask);

  std::vector<std::string> names;
  std::vector<float> cam_v;
  API_ASSERT(PConvFromPyObject(G, objects_py, names));
  API_ASSERT(PConvFromPyObject(G, cam_py, cam_v));
  API_ASSERT(cam_v.size() == 19 || cam_v.size() == 20);

  MetalPickCamera cam;
  std::copy_n(cam_v.begin(), 9, cam.rot);
  std::copy_n(cam_v.begin() + 9, 3, cam.pos);
  std::copy_n(cam_v.begin() + 12, 3, cam.origin);
  cam.tan_half = cam_v[15];
  cam.aspect = cam_v[16];
  cam.clip_front = cam_v[17];
  cam.clip_back = cam_v[18];
  if (cam_v.size() == 20) cam.ortho_half = cam_v[19];

  API_ASSERT(APIEnterBlockedNotModal(G));

  std::vector<ObjectMolecule*> objects;
  objects.reserve(names.size());
  for (const auto& name : names) {
    if (auto* obj = ExecutiveFindObjectMoleculeByName(G, name.c_str()))
      objects.push_back(obj);
  }

  auto hit = MetalPickAtom(G, objects, state, cam, ndc_x, ndc_y, max_ndc2,
                           cluster_ndc2, rep_mask, guide_rep_mask);

  PyObject* result = nullptr;
  if (!hit.obj) {
    result = PConvAutoNone(nullptr);
  } else {
    const AtomInfoType* ai = hit.obj->AtomInfo + hit.atm;
    char resi[8];
    AtomResiFromResv(resi, sizeof(resi), ai);
    result =
        Py_BuildValue("(fssssssff)", hit.d2, hit.obj->Name,
                      LexStr(G, ai->chain), resi, LexStr(G, ai->resn),
                      LexStr(G, ai->segi), LexStr(G, ai->name), hit.sx, hit.sy);
  }

  APIExitBlocked(G);
  return result;
}

static PyObject* CmdGetType(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
  auto res = ExecutiveGetType(G, str1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdGetNames(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1, int2;
  char* str0;
  API_SETUP_ARGS(G, self, args, "Oiis", &self, &int1, &int2, &str0);
  APIEnter(G);
  auto res = ExecutiveGetNames(G, int1, int2, str0);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdInterrupt(PyObject* self, PyObject* args) {
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
  if (ok) {
    PyMOL_SetInterrupt(G->PyMOL, int1);
  }
  return APIResultOk(ok);
}

static PyObject* CmdFlushNow(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  ;
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && G->Ready) {
    /* only called by the GLUT thread with unlocked API, blocked interpreter */
    if (flush_count < 8) { /* prevent super-deep recursion */
      flush_count++;
      PFlushFast(G);
      flush_count--;
    } else {
      PRINTFB(G, FB_CCmd, FB_Warnings)
      " Cmd: PyMOL lagging behind API requests...\n" ENDFB(G);
    }
  }
  return APISuccess();
}

static PyObject* CmdWaitQueue(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  ;

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    int waiting = 0;

    /* called by non-GLUT thread with unlocked API, blocked interpreter */
    if (!G->Terminating) {
      if (APIEnterBlockedNotModal(G)) {
        if (OrthoCommandWaiting(G) || (flush_count > 1))
          waiting = 1; /* commands are waiting or we're in nested execution */
        APIExitBlocked(G);
      } else {
        waiting = 1; /* we're performing a "modal" task... */
      }
    } else {
      waiting = 1; /* we're shutting down... */
    }
    result = PyInt_FromLong(waiting);
  }
  return APIAutoNone(result);
}

static PyObject* CmdWaitDeferred(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  ;
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    if (!G->Terminating) {
      if (APIEnterBlockedNotModal(G)) {
        if (OrthoDeferredWaiting(G))
          result = PyInt_FromLong(1);
        else
          result = PyInt_FromLong(0);
        APIExitBlocked(G);
      }
    }
  }
  return APIAutoNone(result);
}

static PyObject* CmdPaste(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject *list, *str;
  const char* st;
  int l, a;
  int ok = false;
  ok = PyArg_ParseTuple(args, "OO", &self, &list);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    if (!list)
      ok = false;
    else if (!PyList_Check(list))
      ok = false;
    else {
      l = PyList_Size(list);
      for (a = 0; a < l; a++) {
        str = PyList_GetItem(list, a);
        if (str) {
          if (PyString_Check(str)) {
            st = PyString_AsString(str);
            if ((ok = APIEnterNotModal(G))) {
              OrthoPasteIn(G, st);
              if (a < (l - 1)) OrthoPasteIn(G, "\n");
              APIExit(G);
            }
          } else {
            ok = false;
          }
        }
      }
    }
  }
  return APIResultOk(ok);
}

static PyObject* CmdGetWizard(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  auto result = WizardGet(G);
  APIExit(G);
  PyObject* res = result ? result : Py_None;
  return APIIncRef(res);
}

static PyObject* CmdGetWizardStack(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = WizardGetStack(G);
  APIExitBlocked(G);
  return result;
}

static PyObject* CmdSetWizard(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* obj;
  int replace;
  API_SETUP_ARGS(G, self, args, "OOi", &self, &obj, &replace);
  API_ASSERT(APIEnterNotModal(G));
  if (!obj) {
    return APIFailure(G, "Invalid wizard.");
  }
  auto result = WizardSet(G, obj, replace);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetWizardStack(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* obj;
  API_SETUP_ARGS(G, self, args, "OO", &self, &obj);
  API_ASSERT(APIEnterNotModal(G));
  if (!obj) {
    return APIFailure(G, "Invalid wizard.");
  }
  auto result = WizardSetStack(G, obj);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdRefreshWizard(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  WizardRefresh(G);
  OrthoInvalidateDoDraw(G);
  OrthoDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdDirtyWizard(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  WizardDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdSplash(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int query;
  int result = 1;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &query);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (!query) {
    if (ok && (ok = APIEnterNotModal(G))) {
      OrthoSplash(G);
      APIExit(G);
    }
  } else {
  }
  return APIResultCode(result);
}

static PyObject* CmdCls(PyObject* self, PyObject* args) {
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
    OrthoClear(G);
    APIExit(G);
  }
  return APISuccess();
}

static void PyMOLGlobalsCapsuleDestructor(PyObject* self) {
  assert(self != Py_None);
  auto G = _api_get_pymol_globals(self);
  assert(G);
  PyMOL_Free(G->PyMOL);
}

/**
 * Create a `_COb` instance.
 *
 * @param pymol The `pymol` module or a `pymol2.PyMOL` instance
 * @param options Options like `pymol.invocation.options`
 * @param singleton (bool) Whether this instance should be registered as
 * the SingletonPyMOLGlobals (raises RuntimeError if it already exists)
 */
static PyObject* Cmd_New(PyObject* self, PyObject* args) {
  PyObject* pymol = nullptr; /* pymol object instance */
  PyObject* pyoptions = Py_None;
  int singleton = false;

  if (!PyArg_ParseTuple(args, "O|Op", &pymol, &pyoptions, &singleton)) {
    return nullptr;
  }

  if (singleton && SingletonPyMOLGlobals) {
    PyErr_SetString(PyExc_RuntimeError, "Singleton already exists");
    return nullptr;
  }

  CPyMOLOptions* options = PyMOLOptions_New();
  assert(options);

  if (pyoptions != Py_None) {
    PConvertOptions(options, pyoptions);
  }

  CPyMOL* I = PyMOL_NewWithOptions(options);
  PyMOLOptions_Free(options);

  if (!I) {
    PyErr_SetString(PyExc_Exception, "PyMOL_NewWithOptions failed");
    return nullptr;
  }

  PyMOLGlobals* G = PyMOL_GetGlobals(I);
  assert(G);

  if (singleton) {
    assert(!SingletonPyMOLGlobals);
    SingletonPyMOLGlobals = G;
  } else {
    // Creating a non-singleton instance disables auto-library mode
    auto_library_mode_disabled = true;
  }

  G->P_inst = pymol::calloc<CP_inst>(1);
  G->P_inst->obj = pymol;
  G->P_inst->dict = PyObject_GetAttrString(pymol, "__dict__");
  Py_DECREF(G->P_inst->dict);  // borrow reference

  PyObject* tmp = PyCapsule_New(I, nullptr, nullptr);
  PyObject_SetAttrString(pymol, "__pymol__", tmp);
  Py_DECREF(tmp);

  for (int a = 0; a < MAX_SAVED_THREAD; ++a) {
    G->P_inst->savedThread[a].id = -1;
  }

  return PyCapsule_New(PyMOL_GetGlobalsHandle(I), nullptr,
                       PyMOLGlobalsCapsuleDestructor);
}

static PyObject* Cmd_Start(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* cmd = nullptr;
  int ok = true;
  ok = PyArg_ParseTuple(args, "OO", &self, &cmd);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    G->P_inst->cmd = cmd;
    PyMOL_StartWithPython(G->PyMOL);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_Stop(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    PyMOL_Stop(G->PyMOL);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_Idle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  int result = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    result = PyMOL_Idle(G->PyMOL);
    PBlockAndUnlockAPI(G);
  }
  return APIResultCode(result);
}

static PyObject* Cmd_Reshape(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  int width, height, force;
  ok = PyArg_ParseTuple(args, "Oiii", &self, &width, &height, &force);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Reshape(G->PyMOL, width, height, force);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_GetRedisplay(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  int result = false;
  int reset;
  ok = PyArg_ParseTuple(args, "Oi", &self, &reset);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    result = PyMOL_GetRedisplay(G->PyMOL, reset);
    PBlockAndUnlockAPI(G);
  }
  return APIResultCode(result);
}

static PyObject* Cmd_Draw(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Draw(G->PyMOL);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_Button(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  int button, state;
  int x, y, modifiers;
  ok = PyArg_ParseTuple(args, "Oiiiii", &self, &button, &state, &x, &y,
                        &modifiers);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Button(G->PyMOL, button, state, x, y, modifiers);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_Drag(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = true;
  int x, y, modifiers;
  ok = PyArg_ParseTuple(args, "Oiii", &self, &x, &y, &modifiers);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G && G->PyMOL);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && PTryLockAPIAndUnblock(G)) {
    PyMOL_Drag(G->PyMOL, x, y, modifiers);
    PBlockAndUnlockAPI(G);
  }
  return APIResultOk(ok);
}

static PyObject* Cmd_Sdof(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float tx, ty, tz, rx, ry, rz;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Offffff", &self, &tx, &ty, &tz, &rx, &ry, &rz);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    PLockStatus(G);
    ControlSdofUpdate(G, tx, ty, tz, rx, ry, rz);
    PUnlockStatus(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdRunPyMOL(PyObject* self, PyObject* args) {
  PyErr_SetString(PyExc_NotImplementedError, "compile with --glut");
  return nullptr;

  return APISuccess();
}

static PyObject* CmdCountStates(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  OrthoLineType s1;
  int ok = false;
  int count = 0;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  APIEnter(G);
  {
    // A bare object name ("obj" / "?obj"): ask the object, without building a
    // temporary selection over every atom (~5 ms each with 300k atoms loaded,
    // paid per object by the panel poll).
    const char* name = str1[0] == '?' ? str1 + 1 : str1;
    if (auto* obj = ExecutiveFindObjectByName(G, name)) {
      count = obj->getNFrame();
      APIExit(G);
      return APIResultCode(count);
    }
  }
  ok = (SelectorGetTmp2(G, str1, s1) >= 0);
  count = ExecutiveCountStates(G, s1);
  if (count < 0) ok = false;
  SelectorFreeTmp(G, s1);
  APIExit(G);
  return ok ? APIResultCode(count) : APIFailure(G);
}

/* mMol: a hash of every enabled molecule's atom colors and cartoon_color. The
   sequence strip polls it to notice a recolor from anywhere; the same scan in
   Python cost ~30 ms per poll with 300k atoms loaded, this well under 1 ms. */
static PyObject* CmdColorFingerprint(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  unsigned long long h = 1469598103934665603ULL;
  auto mix = [&h](long long v) {
    h = (h ^ (unsigned long long)v) * 1099511628211ULL;
  };
  if (auto objs = ExecutiveFindObjectsByType(G, cObjectMolecule)) {
    for (size_t i = 0, n = VLAGetSize(objs); i < n; ++i) {
      auto obj = static_cast<ObjectMolecule*>(objs[i]);
      if (!obj->Enabled) continue;
      mix(obj->NAtom);
      for (int a = 0; a < obj->NAtom; ++a) mix(obj->AtomInfo[a].color);
      mix(SettingGet_i(G, obj->Setting.get(), nullptr, cSetting_cartoon_color));
    }
    VLAFreeP(objs);
  }
  APIExit(G);
  return PyLong_FromUnsignedLongLong(h);
}

/* mMol: whether any atom of object `name` has its own value for one of the
   settings `indices`. The panel's per-atom transparency badge asks every poll;
   answering "no" here spares it a Python callback per atom. */
static PyObject* CmdHasAtomSettings(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  PyObject* indices_py;
  API_SETUP_ARGS(G, self, args, "OsO", &self, &name, &indices_py);
  std::vector<int> indices;
  API_ASSERT(PConvFromPyObject(G, indices_py, indices));
  APIEnter(G);
  bool found = false;
  if (auto* obj = ExecutiveFindObjectMoleculeByName(G, name)) {
    for (int a = 0; a < obj->NAtom && !found; ++a) {
      const AtomInfoType* ai = obj->AtomInfo + a;
      if (!ai->has_setting) continue;
      for (int s : indices)
        if (SettingUniqueCheck(G, ai->unique_id, s)) {
          found = true;
          break;
        }
    }
  }
  APIExit(G);
  return PyBool_FromLong(found);
}

/**
 * Representations shown on a molecule's atoms, with their atom colors, in one
 * pass over the object (no selection): (mask, {rep: {color: count}}). `mask`
 * has bit 1 << rep for every rep any atom shows; the counts cover the atoms
 * showing that rep, guide atoms only for cartoon. The object card's poll used
 * a `count_atoms('obj & rep X')` per rep plus an `iterate` per shown rep, each
 * a scan of every atom in the session.
 */
static PyObject* CmdObjectRepColors(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  API_SETUP_ARGS(G, self, args, "Os", &self, &name);
  APIEnter(G);
  int mask = 0;
  // Colors in first-seen atom order, so the commonest color's ties resolve as
  // an iterate over the atoms would.
  std::map<int, std::vector<std::pair<int, int>>> counts;
  if (auto* obj = ExecutiveFindObjectMoleculeByName(G, name)) {
    for (int a = 0; a < obj->NAtom; ++a) {
      const AtomInfoType* ai = obj->AtomInfo + a;
      int const vis = ai->visRep & cRepBitmask;
      if (!vis) continue;
      mask |= vis;
      for (int r = 0; r < cRepCnt; ++r) {
        if (!(vis & (1 << r))) continue;
        if (r == cRepCartoon && !(ai->flags & cAtomFlag_guide)) continue;
        auto& cols = counts[r];
        auto it = std::find_if(cols.begin(), cols.end(),
                               [&](auto& c) { return c.first == ai->color; });
        if (it == cols.end())
          cols.emplace_back(ai->color, 1);
        else
          ++it->second;
      }
    }
  }
  APIExit(G);
  PyObject* dict = PyDict_New();
  for (auto& [rep, cols] : counts) {
    PyObject* sub = PyDict_New();
    for (auto& [color, n] : cols) {
      PyObject* k = PyLong_FromLong(color);
      PyObject* v = PyLong_FromLong(n);
      PyDict_SetItem(sub, k, v);
      Py_DECREF(k);
      Py_DECREF(v);
    }
    PyObject* k = PyLong_FromLong(rep);
    PyDict_SetItem(dict, k, sub);
    Py_DECREF(k);
    Py_DECREF(sub);
  }
  return Py_BuildValue("(iN)", mask, dict);
}

/**
 * 1-based indices of `name`'s atoms in residue `resi` of chain `chain` (an
 * empty chain matches any), for `select_list`. Hover picks use it to set the
 * preview without parsing a selection over every atom in the session.
 */
static PyObject* CmdResidueIndices(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *name, *chain, *resi;
  API_SETUP_ARGS(G, self, args, "Osss", &self, &name, &chain, &resi);
  APIEnter(G);
  std::vector<int> idx;
  if (auto* obj = ExecutiveFindObjectMoleculeByName(G, name)) {
    char buf[16];
    for (int a = 0; a < obj->NAtom; ++a) {
      const AtomInfoType* ai = obj->AtomInfo + a;
      if (chain[0] && strcmp(LexStr(G, ai->chain), chain) != 0) continue;
      AtomResiFromResv(buf, sizeof(buf), ai);
      if (strcmp(buf, resi) == 0) idx.push_back(a + 1);
    }
  }
  APIExit(G);
  return PConvToPyObject(idx);
}

/**
 * Sorted "object/chain/resi" keys of the residues named selection `sele`
 * touches, counting guide atoms and non-polymer heavy atoms (ligands, ions,
 * waters), as the sequence strip lists them. Membership is read per atom, with
 * no selection expression to evaluate over the whole session.
 */
static PyObject* CmdSelectionResidueKeys(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sele_name;
  API_SETUP_ARGS(G, self, args, "Os", &self, &sele_name);
  APIEnter(G);
  std::set<std::string> keys;
  SelectorID_t const sele = SelectorIndexByName(G, sele_name);
  if (sele >= 0) {
    ObjectMolecule* obj = nullptr;
    void* hidden = nullptr;
    char resi[16];
    while (ExecutiveIterateObjectMolecule(G, &obj, &hidden)) {
      for (int a = 0; a < obj->NAtom; ++a) {
        const AtomInfoType* ai = obj->AtomInfo + a;
        if (!(ai->flags & cAtomFlag_guide) &&
            ((ai->flags & cAtomFlag_polymer) || ai->isHydrogen()))
          continue;
        if (!SelectorIsMember(G, ai->selEntry, sele)) continue;
        AtomResiFromResv(resi, sizeof(resi), ai);
        keys.insert(std::string(obj->Name) + "/" + LexStr(G, ai->chain) + "/" +
                    resi);
      }
    }
  }
  APIExit(G);
  return PConvToPyObject(std::vector<std::string>(keys.begin(), keys.end()));
}

static PyObject* CmdIdentify(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int mode;
  PyObject* result = nullptr;
  pymol::vla<int> iVLA;
  pymol::vla<ObjectMolecule*> oVLA;

  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &mode);
  APIEnter(G);
  int l = ExecutiveIdentifyObjects(G, str1, mode, &iVLA, &oVLA);
  APIExit(G);

  if (l < 0) {
    return APIFailure(G, "invalid selection");
  }

  if (!iVLA) {
    result = PyList_New(0);
  } else if (!mode) {
    result = PConvIntVLAToPyList(iVLA);
  } else { /* object mode */
    result = PyList_New(l);
    for (int a = 0; a < l; ++a) {
      auto tuple = Py_BuildValue("si", oVLA[a]->Name, iVLA[a]);
      PyList_SetItem(result, a, tuple);
    }
  }

  return result;
}

static PyObject* CmdSystem(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int ok = false;
  int async;
  ok = PyArg_ParseTuple(args, "Osi", &self, &str1, &async);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    if (async) {
      PUnblock(G); /* free up PyMOL and the API */
      ok = system(str1);
      PBlock(G);
    } else if ((ok = APIEnterNotModal(G))) { /* keep PyMOL locked */
      ok = system(str1);
      APIExit(G);
    }
  }
  return APIResultOk(ok);
}

static PyObject* CmdGetFeedback(PyObject* self, PyObject* args) {
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
    if (G->Ready) {
      PyObject* result = nullptr;

      if (G->Terminating) { /* try to bail */
        /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h")
         */
        /* END PROPRIETARY CODE SEGMENT */
        exit(0);
      }

      /* ALLOWED DURING MODAL DRAWING */
      APIEnterBlocked(G);
      auto buffer = OrthoFeedbackOut(G, *G->Ortho);
      APIExitBlocked(G);
      if (!buffer.empty()) result = Py_BuildValue("s", buffer.c_str());
      return (APIAutoNone(result));
    }
  }
  return (APIAutoNone(nullptr));
}

static pymol::Result<PyObject*> _getObjectMoleculeNamePyList(PyMOLGlobals* G,
                                                             char const* sele) {
  if (!sele[0]) {
    return PyList_New(0);
  }
  auto tmpsele1 = SelectorTmp::make(G, sele);
  p_return_if_error(tmpsele1);
  auto const list = ExecutiveGetObjectMoleculeVLA(G, tmpsele1->getName());
  assert(list);
  auto const size = list.size();
  auto result = PyList_New(size);
  if (result) {
    for (unsigned a = 0; a < size; ++a) {
      PyList_SetItem(result, a, PyUnicode_FromString(list[a]->Name));
    }
  }
  return result;
}

static PyObject* CmdGetObjectList(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto res = _getObjectMoleculeNamePyList(G, str1);
  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject* CmdDo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int log;
  int echo;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &log, &echo);
  API_ASSERT(APIEnterNotModal(G));

  if (str1[0] != '_') { /* suppress internal call-backs */
    if (strncmp(str1, "cmd._", 5) && (strncmp(str1, "_cmd.", 5))) {
      if (echo) {
        OrthoAddOutput(G, "mMol>");
        OrthoAddOutput(G, str1);
        OrthoNewLine(G, nullptr, true);
      }
      if (!strncmp(str1, "PyMOL>", 6) || !strncmp(str1, "mMol>", 5)) {
        /* ignore pasted-in /PyMOL>\s?/ or /mMol>\s?/ for end-user convenience
         */
        str1 += (str1[0] == 'P') ? 6 : 5;
        if (str1[0] == ' ') str1++;
      }
      if (log)
        if (WordMatch(G, str1, "quit", true) == 0) /* don't log quit */
          PLog(G, str1, cPLog_pml);
    }
    PParse(G, str1);
  } else if (str1[1] == ' ') {
    /* "_ command" suppresses echoing of command, but it is still logged */
    if (log)
      if (WordMatch(G, str1 + 2, "quit", true) == 0) /* don't log quit */
        PLog(G, str1 + 2, cPLog_pml);
    PParse(G, str1 + 2);
  } else {
    PParse(G, str1);
  }
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdRefreshLater(PyObject* self, PyObject* args) {
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
    SceneInvalidate(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdRefresh(PyObject* self, PyObject* args) {
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
    SceneInvalidateCopy(G, false);
    ExecutiveDrawNow(G); /* TODO STATUS */
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdRefreshNow(PyObject* self, PyObject* args) {
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
    PyMOL_PushValidContext(G->PyMOL); /* we're trusting the caller on this... */
    SceneInvalidateCopy(G, false);
    ExecutiveDrawNow(G); /* TODO STATUS */
    PyMOL_PopValidContext(G->PyMOL);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdMPNG(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, int2, int3, int4, format, mode, quiet;
  int ok = false;
  int width = 0, height = 0;
  ok = PyArg_ParseTuple(args, "Osiiiiiiiii", &self, &str1, &int1, &int2, &int3,
                        &int4, &format, &mode, &quiet, &width, &height);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = MoviePNG(G, str1, SettingGetGlobal_b(G, cSetting_cache_frames), int1,
                  int2, int3, int4, format, mode, quiet, width, height);
    /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(G, ok);
}

static PyObject* CmdDelete(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  API_SETUP_ARGS(G, self, args, "Os", &self, &sname);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveDelete(G, sname);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdDeleteStates(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* objName;
  PyObject* statesPyList;
  API_SETUP_ARGS(G, self, args, "OsO", &self, &objName, &statesPyList);
  std::vector<int> states;
  PConvFromPyObject(G, statesPyList, states);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveDeleteStates(G, objName, states);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdQuit(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int code = EXIT_SUCCESS;
  ok = PyArg_ParseTuple(args, "O|i", &self, &code);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    if (!APIEnterNotModal(G)) { /* override modal state to enable quit */
      PyMOL_SetModalDraw(G->PyMOL, nullptr);
      APIEnter(G);
    }
    if (!G->Option->no_quit) {
      G->Terminating = true;
      PExit(G, code);
    } else {
      OrthoAddOutput(G, "Cmd-Error: cannot quit from within this context.\n");
    }
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdFinishObject(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;

  API_SETUP_ARGS(G, self, args, "Os", &self, &oname);
  API_ASSERT(APIEnterNotModal(G));

  auto origObj = ExecutiveFindObjectByName(G, oname);
  auto ok = bool(origObj);
  if (ok) {
    if (origObj->type == cObjectMolecule) {
      ObjectMoleculeUpdateIDNumbers((ObjectMolecule*)origObj);
      ObjectMoleculeUpdateNonbonded((ObjectMolecule*)origObj);
      origObj->invalidate(cRepAll, cRepInvAll, -1);
    }
    ExecutiveUpdateObjectSelection(G, origObj); /* TODO STATUS */
  }
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdLoadObject(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* oname;
  PyObject* model;
  int frame, type;
  int finish, discrete;
  int quiet;
  int zoom;
  API_SETUP_ARGS(G, self, args, "OsOiiiiii", &self, &oname, &model, &frame,
                 &type, &finish, &discrete, &quiet, &zoom);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveLoadObject(G, oname, model, frame, type, finish, discrete, quiet,
                      zoom);
  APIExit(G);
  return APISuccess();
}

static const PyMethodDef Cmd_methods[] = {
    {"_new", Cmd_New, METH_VARARGS},
    {"_start", Cmd_Start, METH_VARARGS},
    {"_stop", Cmd_Stop, METH_VARARGS},
    {"_idle", Cmd_Idle, METH_VARARGS},
    {"_reshape", Cmd_Reshape, METH_VARARGS},
    {"_getRedisplay", Cmd_GetRedisplay, METH_VARARGS},
    {"_draw", Cmd_Draw, METH_VARARGS},
    {"_button", Cmd_Button, METH_VARARGS},
    {"_drag", Cmd_Drag, METH_VARARGS},
    {"_sdof", Cmd_Sdof, METH_VARARGS},
    /*  {"cache",                 CmdCache,                METH_VARARGS }, */
    {"cls", CmdCls, METH_VARARGS},
    {"coordset_update_thread", CmdCoordSetUpdateThread, METH_VARARGS},
    {"count_states", CmdCountStates, METH_VARARGS},
    {"color_fingerprint", CmdColorFingerprint, METH_VARARGS},
    {"has_atom_settings", CmdHasAtomSettings, METH_VARARGS},
    {"object_rep_colors", CmdObjectRepColors, METH_VARARGS},
    {"residue_indices", CmdResidueIndices, METH_VARARGS},
    {"selection_residue_keys", CmdSelectionResidueKeys, METH_VARARGS},
    {"fake_drag", CmdFakeDrag, METH_VARARGS},
    {"delete", CmdDelete, METH_VARARGS},
    {"delete_states", CmdDeleteStates, METH_VARARGS},
    {"dirty_wizard", CmdDirtyWizard, METH_VARARGS},
    {"do", CmdDo, METH_VARARGS},
    {"finish_object", CmdFinishObject, METH_VARARGS},
    {"flush_now", CmdFlushNow, METH_VARARGS},
    {"get_feedback", CmdGetFeedback, METH_VARARGS},
    {"get_names", CmdGetNames, METH_VARARGS},
    {"get_object_list", CmdGetObjectList, METH_VARARGS},
    {"get_type", CmdGetType, METH_VARARGS},
    {"get_version", CmdGetVersion, METH_VARARGS},
    {"get_capabilities", CmdGetCapabilities, METH_NOARGS,
     "Get a set of compiled-in capabilities"},
    {"get_wizard", CmdGetWizard, METH_VARARGS},
    {"get_wizard_stack", CmdGetWizardStack, METH_VARARGS},
    {"identify", CmdIdentify, METH_VARARGS},
    {"interrupt", CmdInterrupt, METH_VARARGS},
    {"wait_deferred", CmdWaitDeferred, METH_VARARGS},
    {"wait_queue", CmdWaitQueue, METH_VARARGS},
    {"load_object", CmdLoadObject, METH_VARARGS},
    {"metal_pick", CmdMetalPick, METH_VARARGS},
    {"mpng_", CmdMPNG, METH_VARARGS},
    {"object_update_thread", CmdObjectUpdateThread, METH_VARARGS},
    {"paste", CmdPaste, METH_VARARGS},
    {"quit", CmdQuit, METH_VARARGS},
    {"ray_trace_thread", CmdRayTraceThread, METH_VARARGS},
    {"ray_hash_thread", CmdRayHashThread, METH_VARARGS},
    {"ray_anti_thread", CmdRayAntiThread, METH_VARARGS},
    {"refresh", CmdRefresh, METH_VARARGS},
    {"refresh_later", CmdRefreshLater, METH_VARARGS},
    {"refresh_now", CmdRefreshNow, METH_VARARGS},
    {"refresh_wizard", CmdRefreshWizard, METH_VARARGS},
    {"reinitialize", CmdReinitialize, METH_VARARGS},
    {"runpymol", CmdRunPyMOL, METH_VARARGS},
    {"set_wizard", CmdSetWizard, METH_VARARGS},
    {"set_wizard_stack", CmdSetWizardStack, METH_VARARGS},
    {"splash", CmdSplash, METH_VARARGS},
    {"system", CmdSystem, METH_VARARGS},
    {"test2", CmdTest2, METH_VARARGS},
    {NULL, nullptr} /* sentinel */
};

#ifdef __cplusplus
extern "C" {
#endif

PyMODINIT_FUNC PyInit__cmd(void) {
  // The bindings are split by topic into several tables; join them.
  static std::vector<PyMethodDef> methods;
  if (methods.empty()) {
    for (const PyMethodDef* table :
         {Cmd_methods, kEditingMethods, kViewingMethods, kQueryingMethods,
          kIoMethods, kCreatingMethods, kSettingsMethods}) {
      for (; table->ml_name; ++table) methods.push_back(*table);
    }
    methods.push_back({nullptr, nullptr});
  }
  static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT, "pymol._cmd",
                                         "DO NOT USE", -1, methods.data()};
  return PyModule_Create(&moduledef);
}

// Required for AxPyMOL
void init_cmd(void) {
  auto _cmd = PyInit__cmd();
  if (_cmd) {
    PyDict_SetItemString(PyImport_GetModuleDict(), "pymol._cmd", _cmd);
    Py_DECREF(_cmd);
  }
}

#ifdef __cplusplus
}
#endif
