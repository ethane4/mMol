
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
#include "pymol/platform/os_python.h"

#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"
#include "pymol/platform/os_gl.h"
#include "pymol/base/err.h"
#include "pymol/app/main.h"

#include "pymol/platform/memory_debug.h"
#include "pymol/ui/ortho.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/pop_up.h"

#include "pymol/ui/wizard.h"
#include "pymol/scene/scene.h"

#include "pymol/executive/executive.h"
#include "pymol/ui/block.h"
#include "pymol/text/text.h"
#include "pymol/graphics/cgo.h"
#include "pymol/base/vla.h"
#include "pymol/type_traits.h"
#include "pymol/base/util2.h"
#include "pymol/graphics/immediate_helper.h"

#define cWizEventPick 1
#define cWizEventSelect 2
#define cWizEventKey 4
#define cWizEventSpecial 8
#define cWizEventScene 16
#define cWizEventState 32
#define cWizEventFrame 64
#define cWizEventDirty 128
#define cWizEventView 256
#define cWizEventPosition 512

using cWizEvent_t = int;

struct WizardLine {
  int type;
  WordType text;
  OrthoLineType code;
};

struct CWizard : public Block {
  std::vector<unique_PyObject_ptr_auto_gil> Wiz{};
  pymol::vla<WizardLine> Line{};
  ov_size NLine{0};
  int EventMask{0};
  int Dirty{};
  int LastUpdatedState{-1};
  int LastUpdatedFrame{-1};
  float LastUpdatedPosition[3]{};
  SceneViewType LastUpdatedView{};

  CWizard(PyMOLGlobals* G) : Block(G) {};

  bool isEventType(cWizEvent_t eventType) const noexcept;
};

void WizardDirty(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  I->Dirty = true;
  OrthoDirty(G);
}

int WizardUpdate(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  int result = false;

  if (OrthoGetDirty(G)) {
    WizardDoDirty(G);
  }

  {
    int frame = SettingGetGlobal_i(G, cSetting_frame);
    if (frame != I->LastUpdatedFrame) {
      I->LastUpdatedFrame = frame;
      WizardDoFrame(G);
    }
  }
  {
    int state = SettingGetGlobal_i(G, cSetting_state);
    if (state != I->LastUpdatedState) {
      I->LastUpdatedState = state;
      WizardDoState(G);
    }
  }
  WizardDoPosition(G, false);
  WizardDoView(G, false);
  if (I->Dirty) {
    WizardRefresh(G);
    I->Dirty = false;
    result = true;
  }
  return result;
}

void WizardPurgeStack(PyMOLGlobals* G) {
  pymol::pautoblock block(G);
  G->Wizard->Wiz.clear();
}

bool CWizard::isEventType(cWizEvent_t eventType) const noexcept {
  return EventMask & eventType;
}

/**
 * Calls wizard's method in python
 *
 * @param wiz wizard
 * @param funcName wizard method name
 * @param func C/C++ function pointer to call the method (see
 * src/pymol/python/p.cc)
 * @param fargs function arguments for func
 * @return result of func
 */

template <typename Func, typename... FuncArgs>
auto WizardCallPython(PyMOLGlobals* G, PyObject* wiz, const char* funcName,
                      Func&& func, FuncArgs&&... fargs)
    -> decltype(func(std::declval<PyObject*>(), std::declval<const char*>(),
                     fargs...)) {
  using result_t = decltype(func(std::declval<PyObject*>(),
                                 std::declval<const char*>(), fargs...));
  result_t result{};
  assert(wiz != nullptr);
  if (PyObject_HasAttrString(wiz, funcName)) {
    result = func(wiz, funcName, std::forward<FuncArgs>(fargs)...);
    PErrPrintIfOccurred(G);
  }
  return result;
}

/**
 * Run when user selects something with the mouse, in a wizard
 */
int WizardDoSelect(PyMOLGlobals* G, const char* name, int state) {
  CWizard* I = G->Wizard;

  /* if the event is a selection and we're listening for selections */
  if (!I->isEventType(cWizEventSelect)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  /* log if necessary */
  auto buf = pymol::string_format("cmd.get_wizard().do_select('''%s''')", name);
  PLog(G, buf, cPLog_pym);
  /* block and call (in Python) the wizard's do_select */
  pymol::pblock block(G);
  WizardCallPython(G, wiz, "do_pick_state", PTruthCallStr1i, state + 1);
  return WizardCallPython(G, wiz, "do_select", PTruthCallStr, name);
}
void WizardRefresh(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  char* vla = nullptr;
  pymol::pautoblock block(G);

  /* get the current prompt */
  auto wiz = WizardGet(G);
  if (wiz) {
    auto P_list =
        WizardCallPython(G, wiz, "get_prompt", PyObject_CallMethod, "");
    if (P_list) PConvPyListToStringVLA(P_list, &vla);
    Py_XDECREF(P_list);
  }

  OrthoSetWizardPrompt(G, vla);

  /* get the current panel list */

  I->NLine = 0;
  if (wiz) {
    I->EventMask = cWizEventPick + cWizEventSelect;

    auto i =
        WizardCallPython(G, wiz, "get_event_mask", PyObject_CallMethod, "");
    if (i) {
      if (!PConvPyIntToInt(i, &I->EventMask))
        I->EventMask = cWizEventPick + cWizEventSelect;
      Py_XDECREF(i);
    }

    auto P_list =
        WizardCallPython(G, wiz, "get_panel", PyObject_CallMethod, "");
    if (P_list) {
      if (PyList_Check(P_list)) {
        std::size_t ll = PyList_Size(P_list);
        I->Line.check(ll);
        for (std::size_t a = 0u; a < ll; a++) {
          /* fallback defaults */

          I->Line[a].text[0] = 0;
          I->Line[a].code[0] = 0;
          I->Line[a].type = 0;

          i = PyList_GetItem(P_list, a);
          if (PyList_Check(i))
            if (PyList_Size(i) > 2) {
              PConvPyObjectToInt(PyList_GetItem(i, 0), &I->Line[a].type);
              PConvPyObjectToStrMaxLen(PyList_GetItem(i, 1), I->Line[a].text,
                                       sizeof(WordType) - 1);
              PConvPyObjectToStrMaxLen(PyList_GetItem(i, 2), I->Line[a].code,
                                       sizeof(OrthoLineType) - 1);
            }
        }
        I->NLine = ll;
      }
      Py_XDECREF(P_list);
    }
  }
  if (I->NLine) {
    int LineHeight =
        DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
    OrthoReshapeWizard(G, LineHeight * I->NLine + 4);
  } else {
    OrthoReshapeWizard(G, 0);
  }
}
pymol::Result<> WizardSet(PyMOLGlobals* G, PyObject* wiz, bool replace) {
  CWizard* I = G->Wizard;
  pymol::pautoblock block(G);
  if ((!wiz) || (wiz == Py_None) || (!I->Wiz.empty() && replace)) {
    /* remove wizard from stack first */
    if (!I->Wiz.empty()) {
      auto old_wiz = std::move(I->Wiz.back());
      I->Wiz.pop_back();
      if (old_wiz) {
        /* then call cleanup, etc. */
        auto result = WizardCallPython(G, old_wiz.get(), "cleanup",
                                       PyObject_CallMethod, "");
        PXDecRef(result);
      }
    }
  }
  if (wiz && (wiz != Py_None)) { /* push */
    I->Wiz.emplace_back(PIncRef(wiz));
  }
  WizardRefresh(G);
  return {};
}
Block* WizardGetBlock(PyMOLGlobals* G) { return G->Wizard; }
int WizardDoPick(PyMOLGlobals* G, int bondFlag, int state) {
  /**
   * Run when user picks something
   */
  CWizard* I = G->Wizard;

  /* process the pick if it happened and we're listening for it */
  if (!I->isEventType(cWizEventPick)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  if (bondFlag)
    PLog(G, "cmd.get_wizard().do_pick(1)", cPLog_pym);
  else
    PLog(G, "cmd.get_wizard().do_pick(0)", cPLog_pym);

  pymol::pblock block(G);
  WizardCallPython(G, wiz, "do_pick_state", PTruthCallStr1i, state + 1);
  return WizardCallPython(G, wiz, "do_pick", PTruthCallStr1i, bondFlag);
}

int WizardDoKey(PyMOLGlobals* G, unsigned char k, int x, int y, int mod) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventKey)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  auto buffer = pymol::string_format("cmd.get_wizard().do_key(%d,%d,%d,%d)", k,
                                     x, y, mod);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_key", PTruthCallStr4i, k, x, y, mod);
}

int WizardDoPosition(PyMOLGlobals* G, int force) {
  CWizard* I = G->Wizard;
  int result = false;
  if (!I->isEventType(cWizEventPosition)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int changed = force;
  if (!changed) {
    float pos[3];
    SceneGetCenter(G, pos);
    changed = ((fabs(pos[0] - I->LastUpdatedPosition[0]) > R_SMALL4) ||
               (fabs(pos[1] - I->LastUpdatedPosition[1]) > R_SMALL4) ||
               (fabs(pos[2] - I->LastUpdatedPosition[2]) > R_SMALL4));
  }
  if (changed) {
    SceneGetCenter(G, I->LastUpdatedPosition);
    pymol::pblock block(G);
    result = WizardCallPython(G, wiz, "do_position", PTruthCallStr0);
  }
  return result;
}

int WizardDoView(PyMOLGlobals* G, int force) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventView)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int result = false;
  int changed = force;
  if (!changed) {
    SceneViewType view;
    SceneGetView(G, view);
    changed = !SceneViewEqual(view, I->LastUpdatedView);
  }
  if (changed) {
    SceneGetView(G, I->LastUpdatedView);
    pymol::pblock block(G);
    result = WizardCallPython(G, wiz, "do_view", PTruthCallStr0);
  }
  return result;
}

int WizardDoScene(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventScene)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  std::string buffer = "cmd.get_wizard().do_scene()";
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_scene", PTruthCallStr0);
}

int WizardDoDirty(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventDirty)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  std::string buffer = "cmd.get_wizard().do_dirty()";
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_dirty", PTruthCallStr0);
}

int WizardDoState(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventState)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int state = SettingGet<int>(G, cSetting_state);
  auto buffer = pymol::string_format("cmd.get_wizard().do_state(%d)", state);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_state", PTruthCallStr1i, state);
}

int WizardDoFrame(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  if (!I->isEventType(cWizEventFrame)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  int frame = SettingGet<int>(G, cSetting_frame) + 1;
  auto buffer = pymol::string_format("cmd.get_wizard().do_frame(%d)", frame);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_frame", PTruthCallStr1i, frame);
}

int WizardDoSpecial(PyMOLGlobals* G, int k, int x, int y, int mod) {
  CWizard* I = G->Wizard;

  if (!I->isEventType(cWizEventSpecial)) {
    return false;
  }
  auto wiz = WizardGet(G);
  if (!wiz) {
    return false;
  }
  auto buffer = pymol::string_format("cmd.get_wizard().do_special(%d,%d,%d,%d)",
                                     k, x, y, mod);
  PLog(G, buffer, cPLog_pym);
  pymol::pblock block(G);
  return WizardCallPython(G, wiz, "do_special", PTruthCallStr4i, k, x, y, mod);
}
PyObject* WizardGet(PyMOLGlobals* G) {
  auto I = G->Wizard;
  return I->Wiz.empty() ? nullptr : I->Wiz.back().get();
}
PyObject* WizardGetStack(PyMOLGlobals* G) {
  CWizard* I = G->Wizard;
  auto result = PyList_New(I->Wiz.size());
  for (std::size_t i = 0; i < I->Wiz.size(); i++) {
    auto wiz = I->Wiz[i].get();
    Py_INCREF(wiz);
    PyList_SetItem(result, i, wiz); /* steals ref */
  }
  return result;
}
pymol::Result<> WizardSetStack(PyMOLGlobals* G, PyObject* list) {
  CWizard* I = G->Wizard;

  if (list == nullptr || !PyList_Check(list)) {
    return pymol::make_error("Invalid list.");
  }
  WizardPurgeStack(G);
  auto size = PyList_Size(list);
  pymol::pautoblock block(G);
  for (int a = 0; a < size; a++) {
    auto wiz = PyList_GetItem(list, a);
    I->Wiz.emplace_back(PIncRef(wiz));
  }
  WizardRefresh(G);
  OrthoDirty(G);
  return {};
}
int WizardInit(PyMOLGlobals* G) {
  auto I = new CWizard(G);
  G->Wizard = I;
  I->active = true;
  I->TextColor[0] = 0.2F;
  I->TextColor[1] = 1.0F;
  I->TextColor[2] = 0.2F;

  OrthoAttach(G, I, cOrthoTool);

  I->Line.resize(1);
  return 1;
}
void WizardFree(PyMOLGlobals* G) {
  WizardPurgeStack(G);
  DeleteP(G->Wizard);
}
