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

// Python bindings for viewing, movies and controls (pymol/viewing.py,
// moving.py, controlling.py, keyboard.py).

#include "pymol/base/err.h"
#include "pymol/executive/executive.h"
#include "pymol/scene/movie.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/control.h"
#include "pymol/base/feedback.h"
#include "pymol/scene/movie_scene.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdDelColorection(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  PyObject* list;
  char* prefix;
  API_SETUP_ARGS(G, self, args, "OOs", &self, &list, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
  ok = SelectorColorectionFree(G, list, prefix);
  APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdSetColorection(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  char* prefix;
  PyObject* list;
  API_SETUP_ARGS(G, self, args, "OOs", &self, &list, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
  ok = SelectorColorectionApply(G, list, prefix);
  APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdGetColorection(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  const char* prefix;
  API_SETUP_ARGS(G, self, args, "Os", &self, &prefix);
  API_ASSERT(APIEnterBlockedNotModal(G));
  result = SelectorColorectionGet(G, prefix);
  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject* CmdGetVis(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterBlockedNotModal(G));
  result = ExecutiveGetVisAsPyDict(G);
  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject* CmdSetVis(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  PyObject* visDict;
  API_SETUP_ARGS(G, self, args, "OO", &self, &visDict);
  API_ASSERT(APIEnterBlockedNotModal(G));
  ok = ExecutiveSetVisFromPyDict(G, visDict);
  APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdSpectrum(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *expr, *prefix;
  float min, max;
  int digits, start, stop, byres;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossffiisiii", &self, &str1, &expr, &min, &max,
                 &start, &stop, &prefix, &digits, &byres, &quiet);
  API_ASSERT(APIEnterNotModal(G));

  auto res = ExecutiveSpectrum(G, str1, expr, min, max, start, stop, prefix,
                               digits, byres, quiet);

  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdMDump(PyObject* self, PyObject* args) {
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
    APIEnter(G);
    MovieDump(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdAccept(PyObject* self, PyObject* args) {
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
    MovieSetLock(G, false);
    PRINTFB(G, FB_Movie, FB_Actions)
    " Movie: Risk accepted by user.  Movie commands have been enabled.\n" ENDFB(
        G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdDecline(PyObject* self, PyObject* args) {
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
    MovieReset(G);
    PRINTFB(G, FB_Movie, FB_Actions)
    " Movie: Risk declined by user.  Movie commands have been deleted.\n" ENDFB(
        G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdScene(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;

  MovieSceneFuncArgs margs;
  const char *key, *action, *message = nullptr, *new_key = nullptr;
  const char* sele = "all";

  API_SETUP_ARGS(G, self, args, "Oss|zbbbbbfzbbs", &self, &key, &action,
                 &message, &margs.store_view, &margs.store_color,
                 &margs.store_active, &margs.store_rep, &margs.store_frame,
                 &margs.animate, &new_key, &margs.hand, &margs.quiet, &sele);
  API_ASSERT(APIEnterBlockedNotModal(G));

  margs.key = key;
  margs.action = action;
  margs.message = message ? message : "";
  margs.new_key = new_key ? new_key : "";
  margs.sele = sele;
  auto res = MovieSceneFunc(G, margs);
  APIExitBlocked(G);
  return APIResult(G, res);
}

static PyObject* CmdSceneOrder(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;

  const char* location;
  unsigned char sort;
  PyObject* pynames = nullptr;

  API_SETUP_ARGS(G, self, args, "OObs", &self, &pynames, &sort, &location);

  std::vector<std::string> names;
  API_ASSERT(PConvFromPyObject(G, pynames, names));

  API_ASSERT(APIEnterBlockedNotModal(G));

  auto result = MovieSceneOrder(G, std::move(names), sort, location);

  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdGetSceneOrder(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;

  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnterBlocked(G);

  result = PConvToPyObject(MovieSceneGetOrder(G));

  APIExitBlocked(G);
  return (APIAutoNone(result));
}

static PyObject* CmdGetSceneThumbnail(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* name;

  API_SETUP_ARGS(G, self, args, "Os", &self, &name);
  APIEnterBlocked(G);

  png_outbuf_t pngbuf = MovieSceneGetThumbnail(G, name);

  APIExitBlocked(G);

  return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(pngbuf.data()),
                                   pngbuf.size());
}

static PyObject* CmdGetSceneMessage(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* name;

  API_SETUP_ARGS(G, self, args, "Os", &self, &name);
  APIEnterBlocked(G);

  std::string message = MovieSceneGetMessage(G, name);

  APIExitBlocked(G);

  auto result = PyString_FromString(message.c_str());

  return APIAutoNone(result);
}

static PyObject* CmdSetSceneMessage(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* name;
  const char* message;

  API_SETUP_ARGS(G, self, args, "Oss", &self, &name, &message);
  APIEnterBlocked(G);

  // TODO : Needs to be processed by undo.
  auto result = MovieSceneSetMessage(G, name, message);

  APIExitBlocked(G);
  if (result) {
    return APISuccess();
  }
  return APIFailure(G, result.error());
}

static PyObject* CmdLoadPNG(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int ok = false;
  int quiet;
  int movie, stereo;
  API_SETUP_ARGS(G, self, args, "Osiii", &self, &str1, &movie, &stereo, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ok = SceneLoadPNG(G, str1, movie, stereo, quiet);
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdBackgroundColor(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  API_SETUP_ARGS(G, self, args, "Os", &self, &str1);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveBackgroundColor(G, str1);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdGetMoviePlaying(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  ;
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
    result = PyInt_FromLong(MoviePlaying(G));
  } else {
    API_HANDLE_ERROR;
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetView(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  SceneViewType view;
  API_SETUP_ARGS(G, self, args, "O", &self);
  {
    APIEnter(G);
    SceneGetView(G, view);
    APIExit(G);
    return (Py_BuildValue("(fffffffffffffffffffffffff)", view[0], view[1],
                          view[2], view[3],
                          /* 4x4 mat */
                          view[4], view[5], view[6], view[7], view[8], view[9],
                          view[10], view[11], view[12], view[13], view[14],
                          view[15], view[16], view[17], view[18],
                          /* pos */
                          view[19], view[20], view[21], /* origin */
                          view[22], view[23],           /* clip */
                          view[24]                      /* orthoscopic */
                          ));
  }
}

static PyObject* CmdSetView(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  SceneViewType view;
  int quiet;
  float animate;
  int hand;
  API_SETUP_ARGS(G, self, args, "O(fffffffffffffffffffffffff)ifi", &self,
                 &view[0], &view[1], &view[2], &view[3], /* 4x4 mat */
                 &view[4], &view[5], &view[6], &view[7], &view[8], &view[9],
                 &view[10], &view[11], &view[12], &view[13], &view[14],
                 &view[15], &view[16], &view[17], &view[18], /* pos */
                 &view[19], &view[20], &view[21],            /* origin */
                 &view[22], &view[23],                       /* clip */
                 &view[24],                                  /* orthoscopic */
                 &quiet, &animate, &hand);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetView(G, view, quiet, animate, hand); /* TODO STATUS */
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdGetState(PyObject* self, PyObject* args) {
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
    result = SceneGetState(G); /* shouldn't this be +1? */
  }
  return (APIResultCode(result));
}

static PyObject* CmdGetFrame(PyObject* self, PyObject* args) {
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
    result = SceneGetFrame(G) + 1;
  }
  return (APIResultCode(result));
}

static PyObject* CmdMask(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, quiet;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMask(G, str1, int1, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdButton(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1, i2;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oii", &self, &i1, &i2);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ButModeSet(G, i1, i2); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdLabel(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  ExecutiveLabel(G, str1, str2, quiet, cExecutiveLabelEvalOn);
  APIExitBlocked(G);
  if (PyErr_Occurred()) {
    return nullptr;
  }
  return APISuccess();
}

static PyObject* CmdLabel2(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveLabel(G, str1, str2, quiet, cExecutiveLabelEvalAlt);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdRecolor(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  cRep_t rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &rep);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveInvalidateRep(G, str1, rep, cRepInvColor);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdRebuild(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* str1;
  cRep_t rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &str1, &rep);
  API_ASSERT(APIEnterNotModal(G));

  pymol::Result<> res;
  if (WordMatchExact(G, str1, cKeywordAll, true)) {
    ExecutiveRebuildAll(G);
  } else {
    const cRepInv_t level = SettingGet<bool>(G, cSetting_defer_builds_mode)
                                ? cRepInvPurgeAll
                                : cRepInvAll;
    res = ExecutiveInvalidateRep(G, str1, rep, level);
  }

  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdResetRate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  ButModeResetRate(G);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdOrient(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  float animate;
  int quiet = false; /* TODO */
  API_SETUP_ARGS(G, self, args, "Osif", &self, &str1, &state, &animate);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveOrient(G, str1, state, animate, false, 0.0F, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdDirty(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  API_ASSERT(APIEnterNotModal(G));
  OrthoDirty(G);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdRock(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &int1);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ControlRock(G, int1);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdSetBusy(PyObject* self, PyObject* args) {
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
    PLockStatus(G);
    PyMOL_SetBusy(G->PyMOL, int1);
    PUnlockStatus(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdSetFrame(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int mode, frm;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &mode, &frm);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetFrame(G, mode, frm);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdFrame(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int frm, trigger;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &frm, &trigger);
  API_ASSERT(APIEnterNotModal(G));
  SceneSetFrame(G, trigger ? 4 : 0, frm);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdStereo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int i1;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &i1);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveStereo(G, i1);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdReset(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* obj;
  API_SETUP_ARGS(G, self, args, "Os", &self, &obj);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveReset(G, obj);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMDo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* cmd;
  int frame;
  int append;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oisi", &self, &frame, &cmd, &append);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    if (frame < 0) {
      if (frame == -1) {
        frame = SceneGetFrame(G);
      } else {
        frame = MovieGetLength(G) + 2 + frame;
        if (frame < 0) frame = 0;
      }
    }
    if (append) {
      MovieAppendCommand(G, frame, cmd);
    } else {
      MovieSetCommand(G, frame, cmd);
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdMPlay(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int cmd;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &cmd);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    MoviePlay(G, cmd);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdMMatrix(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int cmd;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Oi", &self, &cmd);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = MovieMatrix(G, cmd);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdMClear(PyObject* self, PyObject* args) {
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
    MovieClearImages(G);
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdMSet(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int start_from, freeze;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &start_from, &freeze);
  API_ASSERT(APIEnterNotModal(G));
  MovieSet(G, str1, start_from, freeze);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdMModify(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* object;
  int action, index, count, target, freeze, quiet;
  API_SETUP_ARGS(G, self, args, "Oiiiisii", &self, &action, &index, &count,
                 &target, &object, &freeze, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveMotionViewModify(G, static_cast<ViewElemAction>(action), index,
                                count, target, object, freeze, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMView(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int action, first, last, simple, wrap, window, cycles, quiet, state, autogen;
  float power, bias, linear, hand, scene_cut;
  char *object, *scene_name;
  ok = PyArg_ParseTuple(args, "Oiiiffifsiiiisfiii", &self, &action, &first,
                        &last, &power, &bias, &simple, &linear, &object, &wrap,
                        &hand, &window, &cycles, &scene_name, &scene_cut,
                        &quiet, &state, &autogen);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ok = ExecutiveMotionView(G, static_cast<MViewAction>(action), first, last,
                             power, bias, simple, linear, object, wrap, hand,
                             window, cycles, scene_name, scene_cut, state,
                             quiet, autogen);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdViewport(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int w, h;
  API_SETUP_ARGS(G, self, args, "Oii", &self, &w, &h);
  API_ASSERT(APIEnterNotModal(G));
  {
    {
      if (!(w < 1 && h < 1 && ExecutiveIsFullScreen(G))) {
        if (((w > 0) && (h <= 0)) || ((h > 0) && (w <= 0))) {
          int cw, ch;
          SceneGetWidthHeight(G, &cw, &ch);
          if (h <= 0) {
            h = (w * ch) / cw;
          }
          if (w <= 0) {
            w = (h * cw) / ch;
          }
        }

        if ((w > 0) && (h > 0)) {
          if (w < 10) w = 10;
          if (h < 10) h = 10;

          if (SettingGetGlobal_b(G, cSetting_internal_gui)) {
            w += DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_width));
          }

          if (SettingGetGlobal_i(G, cSetting_internal_feedback)) {
            h += (SettingGetGlobal_i(G, cSetting_internal_feedback) - 1) *
                     cOrthoLineHeight +
                 cOrthoBottomSceneMargin;
          }

          h += MovieGetPanelHeight(G);
        } else {
          w = -1;
          h = -1;
        }
      } else {
        w = 0;
        h = 0;
      }

      {
        PyMOL_NeedReshape(G->PyMOL, 2, 0, 0, w, h);
      }
      APIExit(G);
    }
  }
  return APISuccess();
}

static PyObject* CmdColor(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *color;
  int flags;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Ossii", &self, &color, &str1, &flags, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveColorFromSele(G, str1, color, flags, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdColorDef(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* color;
  float v[3];
  int mode;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osfffii", &self, &color, v, v + 1, v + 2,
                 &mode, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  ColorDef(G, color, v, mode, quiet);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdDraw(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1, int2;
  int quiet, antialias;

  API_SETUP_ARGS(G, self, args, "Oiiii", &self, &int1, &int2, &antialias,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));

  bool entire_window = false;

  if (antialias == -2) {
    // capture action
    entire_window = true;
    int1 = 0;
    int2 = 0;
    antialias = 0;
  }

  bool ok = ExecutiveDrawCmd(G, int1, int2, antialias, entire_window, quiet);
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdRay(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int w, h, mode;
  float angle, shift;
  int quiet;
  int antialias;
  API_SETUP_ARGS(G, self, args, "Oiiiffii", &self, &w, &h, &antialias, &angle,
                 &shift, &mode, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  {
    if (mode < 0) mode = SettingGetGlobal_i(G, cSetting_ray_default_renderer);
    ExecutiveRay(G, w, h, mode, angle, shift, quiet, false,
                 antialias); /* TODO STATUS */
    APIExit(G);
  }
  return APISuccess();
}

static PyObject* CmdClip(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  float dist;
  char* str1;
  int state;
  API_SETUP_ARGS(G, self, args, "Osfsi", &self, &sname, &dist, &str1, &state);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result;
  {
    SelectorTmp2 s1(G, str1);
    result = SceneClipFromMode(G, sname, dist, s1.getName(), state);
  }
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdGetClip(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  API_SETUP_ARGS(G, self, args, "O", &self);
  APIEnter(G);
  auto result = SceneGetClip(G);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMove(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  float dist;
  API_SETUP_ARGS(G, self, args, "Osf", &self, &sname, &dist);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMove(G, sname, dist);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTurn(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  float angle;
  API_SETUP_ARGS(G, self, args, "Osf", &self, &sname, &angle);
  API_ASSERT(APIEnterNotModal(G));
  SceneRotateAxis(G, angle, sname[0]);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdSet(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int index;
  PyObject* value;
  char* str;
  int state;
  int quiet;
  int updates;
  API_SETUP_ARGS(G, self, args, "OiOsiii", &self, &index, &value, &str, &state,
                 &quiet, &updates);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveSetSetting(G, index, value, str, state, quiet, updates);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdCartoon(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  int type;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &type);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveCartoon(G, type, sname);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdShowHide(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* sname;
  int rep;
  int state;
  API_SETUP_ARGS(G, self, args, "Osii", &self, &sname, &rep, &state);
  API_ASSERT(APIEnterNotModal(G));
  if (sname[0] == '@') {
    // DEPRECATED
    sname = cKeywordAll;
    rep = cRepBitmask;
  }
  pymol::Result<> res;
  {
    auto tmpsele1 = SelectorTmp2::make(G, sname);
    if (tmpsele1) {
      ExecutiveSetRepVisMask(G, tmpsele1->getName(), rep, state);
    } else {
      res = tmpsele1.error_move();
    }
  }
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdOnOffBySele(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  int onoff;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &onoff);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveSetOnOffBySele(G, sname, onoff);
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdToggle(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sname;
  int rep;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &sname, &rep);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveToggleRepVisib(G, sname, rep);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdFullScreen(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int flag = 0;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Oi", &self, &flag);
  API_ASSERT(APIEnterNotModal(G));
  ExecutiveFullScreen(G, flag); /* TODO STATUS */
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdOrigin(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *obj;
  float v[3];
  int state;
  API_SETUP_ARGS(G, self, args, "Oss(fff)i", &self, &str1, &obj, v, v + 1,
                 v + 2, &state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveOrigin(G, str1, 1, obj, v, state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdCenter(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  int origin;
  float animate;
  int quiet = false; /* TODO */
  API_SETUP_ARGS(G, self, args, "Osiif", &self, &str1, &state, &origin,
                 &animate);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> res;
  {
    auto tmpsele1 = SelectorTmp2::make(G, str1);
    if (tmpsele1) {
      ExecutiveCenter(G, tmpsele1->getName(), state, origin, animate, nullptr,
                      quiet);
    } else {
      res = tmpsele1.error_move();
    }
  }
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdZoom(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  float buffer;
  int state;
  int inclusive;
  float animate;
  int quiet = false; /* TODO */
  API_SETUP_ARGS(G, self, args, "Osfiif", &self, &str1, &buffer, &state,
                 &inclusive, &animate);
  API_ASSERT(APIEnterNotModal(G));
  {
    SelectorTmp2 s1(G, str1);
    ExecutiveWindowZoom(G, s1.getName(), buffer, state, inclusive, animate,
                        quiet);
  }
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdOrder(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int int1, int2;

  API_SETUP_ARGS(G, self, args, "Osii", &self, &str1, &int1, &int2);
  APIEnterNotModal(G);
  auto result = ExecutiveOrder(G, str1, int1, int2);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdScrollTo(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int i = 0, r = -1;

  if (!PyArg_ParseTuple(args, "Os|i", &self, &name, &i)) {
    API_HANDLE_ERROR;
    ok_raise(1);
  }

  ok_assert(1, name && name[0]);

  API_SETUP_PYMOL_GLOBALS;
  ok_assert(2, G && APIEnterBlockedNotModal(G));

  r = ExecutiveScrollTo(G, name, i);

ok_except2:
  APIExitBlocked(G);
ok_except1:
  return Py_BuildValue("i", r);
}

static PyObject* CmdWindow(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int int1, x, y, width, height;
  API_SETUP_ARGS(G, self, args, "Oiiiii", &self, &int1, &x, &y, &width,
                 &height);
  API_ASSERT(G->HaveGUI);
  API_ASSERT(APIEnterNotModal(G));
  {
    APIExit(G);
  }
  return APISuccess();
}

const PyMethodDef kViewingMethods[] = {
    {"accept", CmdAccept, METH_VARARGS},
    {"bg_color", CmdBackgroundColor, METH_VARARGS},
    {"button", CmdButton, METH_VARARGS},
    {"cartoon", CmdCartoon, METH_VARARGS},
    {"center", CmdCenter, METH_VARARGS},
    {"clip", CmdClip, METH_VARARGS},
    {"color", CmdColor, METH_VARARGS},
    {"colordef", CmdColorDef, METH_VARARGS},
    {"decline", CmdDecline, METH_VARARGS},
    {"del_colorection", CmdDelColorection, METH_VARARGS},
    {"dirty", CmdDirty, METH_VARARGS},
    {"draw", CmdDraw, METH_VARARGS},
    {"frame", CmdFrame, METH_VARARGS},
    {"full_screen", CmdFullScreen, METH_VARARGS},
    {"get_clip", CmdGetClip, METH_VARARGS},
    {"get_colorection", CmdGetColorection, METH_VARARGS},
    {"get_frame", CmdGetFrame, METH_VARARGS},
    {"get_movie_playing", CmdGetMoviePlaying, METH_VARARGS},
    {"get_state", CmdGetState, METH_VARARGS},
    {"get_view", CmdGetView, METH_VARARGS},
    {"get_vis", CmdGetVis, METH_VARARGS},
    {"label", CmdLabel, METH_VARARGS},
    {"label2", CmdLabel2, METH_VARARGS},
    {"load_png", CmdLoadPNG, METH_VARARGS},
    {"mask", CmdMask, METH_VARARGS},
    {"mclear", CmdMClear, METH_VARARGS},
    {"mdo", CmdMDo, METH_VARARGS},
    {"mdump", CmdMDump, METH_VARARGS},
    {"mmodify", CmdMModify, METH_VARARGS},
    {"move", CmdMove, METH_VARARGS},
    {"mset", CmdMSet, METH_VARARGS},
    {"mplay", CmdMPlay, METH_VARARGS},
    {"mmatrix", CmdMMatrix, METH_VARARGS},
    {"mview", CmdMView, METH_VARARGS},
    {"origin", CmdOrigin, METH_VARARGS},
    {"orient", CmdOrient, METH_VARARGS},
    {"onoff_by_sele", CmdOnOffBySele, METH_VARARGS},
    {"order", CmdOrder, METH_VARARGS},
    {"scrollto", CmdScrollTo, METH_VARARGS},
    {"rebuild", CmdRebuild, METH_VARARGS},
    {"recolor", CmdRecolor, METH_VARARGS},
    {"render", CmdRay, METH_VARARGS},
    {"reset", CmdReset, METH_VARARGS},
    {"reset_rate", CmdResetRate, METH_VARARGS},
    {"rock", CmdRock, METH_VARARGS},
    {"set", CmdSet, METH_VARARGS},
    {"scene", CmdScene, METH_VARARGS},
    {"scene_order", CmdSceneOrder, METH_VARARGS},
    {"get_scene_order", CmdGetSceneOrder, METH_VARARGS},
    {"get_scene_thumbnail", CmdGetSceneThumbnail, METH_VARARGS},
    {"get_scene_message", CmdGetSceneMessage, METH_VARARGS},
    {"set_scene_message", CmdSetSceneMessage, METH_VARARGS},
    {"set_busy", CmdSetBusy, METH_VARARGS},
    {"set_colorection", CmdSetColorection, METH_VARARGS},
    {"set_frame", CmdSetFrame, METH_VARARGS},
    {"set_view", CmdSetView, METH_VARARGS},
    {"set_vis", CmdSetVis, METH_VARARGS},
    {"showhide", CmdShowHide, METH_VARARGS},
    {"spectrum", CmdSpectrum, METH_VARARGS},
    {"stereo", CmdStereo, METH_VARARGS},
    {"toggle", CmdToggle, METH_VARARGS},
    {"turn", CmdTurn, METH_VARARGS},
    {"viewport", CmdViewport, METH_VARARGS},
    {"window", CmdWindow, METH_VARARGS},
    {"zoom", CmdZoom, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
