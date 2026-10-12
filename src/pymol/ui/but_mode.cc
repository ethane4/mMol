
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
#include "pymol/platform/memory_debug.h"
#include "pymol/app/main.h"
#include "pymol/base/base.h"
#include "pymol/ui/but_mode.h"
#include "pymol/scene/scene.h"
#include "pymol/base/util.h"
#include "pymol/ui/ortho.h"
#include "pymol/settings/setting.h"
#include "pymol/python/p.h"
#include "pymol/text/text.h"
#include "pymol/ui/menu.h"
#include "pymol/graphics/cgo.h"
#include "pymol/scene/movie.h"

struct CButMode : public Block {
  int NCode{};
  int Mode[cButModeInputCount]{};
  int NBut{};
  float Rate{0.0f};
  float RateShown{0.0f};
  float Samples{0.0f}, Delay{0.0f};
  int DeferCnt{0};
  float DeferTime{0.0f};

  CButMode(PyMOLGlobals* G) : Block(G) {}
};
Block* ButModeGetBlock(PyMOLGlobals* G) {
  CButMode* I = G->ButMode;
  {
    return (I);
  }
}

int ButModeGet(PyMOLGlobals* G, int button) {
  CButMode* I = G->ButMode;
  if ((button >= 0) && (button < I->NBut)) {
    return I->Mode[button];
  }
  return 0;
}

void ButModeSet(PyMOLGlobals* G, int button, int action) {
  CButMode* I = G->ButMode;
  if ((button >= 0) && (button < I->NBut) && (action >= 0) &&
      (action < I->NCode)) {
    I->Mode[button] = action;
    OrthoDirty(G);
  }
}
void ButModeSetRate(PyMOLGlobals* G, float interval) {
  CButMode* I = G->ButMode;

  if (interval >= 0.001F) { /* sub-millisecond, defer... */
    if (I->DeferCnt) {
      interval = (interval + I->DeferTime) / (I->DeferCnt + 1);
      I->DeferCnt = 0;
      I->DeferTime = 0.0F;
    }
    I->Delay -= interval;
    if (interval < 1.0F) {
      I->Samples *= 0.95 * (1.0F - interval);
      I->Rate *= 0.95 * (1.0F - interval);
    } else {
      I->Samples = 0.0F;
      I->Rate = 0.0F;
    }

    I->Samples++;
    I->Rate += 1.0F / interval;
  } else {
    I->DeferCnt++;
    I->DeferTime += interval;
  }
}
void ButModeResetRate(PyMOLGlobals* G) {
  CButMode* I = G->ButMode;
  I->Samples = 0.0;
  I->Rate = 0.0;
  I->RateShown = 0;
  I->Delay = 0;
}
void ButModeFree(PyMOLGlobals* G) { DeleteP(G->ButMode); }
int ButModeInit(PyMOLGlobals* G) {
  CButMode* I = nullptr;
  if ((I = (G->ButMode = new CButMode(G)))) {
    int a;

    I->Rate = 0.0;
    I->Samples = 0.0;
    I->RateShown = 0.0;
    I->Delay = 0.0;
    I->DeferCnt = 0;
    I->DeferTime = 0.0F;
    I->NCode = cButModeCount;
    I->NBut = cButModeInputCount;

    for (a = 0; a < I->NBut; a++) {
      I->Mode[a] = -1;
    }

    I->active = true;

    I->TextColor[0] = 0.2F;
    I->TextColor[1] = 1.0F;
    I->TextColor[2] = 0.2F;

    OrthoAttach(G, I, cOrthoTool);
    return 1;
  } else
    return 0;
}
int ButModeCheckPossibleSingleClick(PyMOLGlobals* G, int button, int mod) {
  int click_button = -1;
  switch (button) {
    case P_GLUT_LEFT_BUTTON:
      click_button = P_GLUT_SINGLE_LEFT;
      break;
    case P_GLUT_MIDDLE_BUTTON:
      click_button = P_GLUT_SINGLE_MIDDLE;
      break;
    case P_GLUT_RIGHT_BUTTON:
      click_button = P_GLUT_SINGLE_RIGHT;
      break;
  }
  if (click_button < 0)
    return false;
  else
    return (ButModeTranslate(G, click_button, mod) >= 0);
}

int ButModeTranslate(PyMOLGlobals* G, int button, int mod) {
  int mode = cButModeNothing;
  CButMode* I = G->ButMode;
  switch (button) {
    case P_GLUT_LEFT_BUTTON:
      mode = 0;
      break;
    case P_GLUT_MIDDLE_BUTTON:
      mode = 1;
      break;
    case P_GLUT_RIGHT_BUTTON:
      mode = 2;
      break;
    case P_GLUT_BUTTON_SCROLL_FORWARD:
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      switch (mod) {
        case 0:
          mode = 12;
          break;
        case cOrthoSHIFT:
          mode = 13;
          break;
        case cOrthoCTRL:
          mode = 14;
          break;
        case (cOrthoCTRL + cOrthoSHIFT):
          mode = 15;
      }
      mod = 0;
      switch (I->Mode[mode]) {
        case cButModeScaleSlab:
          if (button == P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeScaleSlabExpand;
          } else {
            return cButModeScaleSlabShrink;
          }
          break;
        case cButModeMoveSlab:
          if (button == P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeMoveSlabForward;
          } else {
            return cButModeMoveSlabBackward;
          }
          break;
        case cButModeMoveSlabAndZoom:
          if (button == P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeMoveSlabAndZoomForward;
          } else {
            return cButModeMoveSlabAndZoomBackward;
          }
          break;
        case cButModeInvMoveSlabAndZoom:
          if (button != P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeMoveSlabAndZoomForward;
          } else {
            return cButModeMoveSlabAndZoomBackward;
          }
          break;
        case cButModeTransZ:
          if (button == P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeZoomForward;
          } else {
            return cButModeZoomBackward;
          }
          break;
        case cButModeInvTransZ:
          if (button != P_GLUT_BUTTON_SCROLL_FORWARD) {
            return cButModeZoomForward;
          } else {
            return cButModeZoomBackward;
          }
          break;
      }
      return -1;
      break;
    case P_GLUT_DOUBLE_LEFT:
    case P_GLUT_DOUBLE_MIDDLE:
    case P_GLUT_DOUBLE_RIGHT:
    case P_GLUT_SINGLE_LEFT:
    case P_GLUT_SINGLE_MIDDLE:
    case P_GLUT_SINGLE_RIGHT:
      switch (button) {
        case P_GLUT_DOUBLE_LEFT:
          mode = 16;
          break;
        case P_GLUT_DOUBLE_MIDDLE:
          mode = 17;
          break;
        case P_GLUT_DOUBLE_RIGHT:
          mode = 18;
          break;
        case P_GLUT_SINGLE_LEFT:
          mode = 19;
          break;
        case P_GLUT_SINGLE_MIDDLE:
          mode = 20;
          break;
        case P_GLUT_SINGLE_RIGHT:
          mode = 21;
          break;
      }
      switch (mod) {
        case cOrthoSHIFT:
          mode += 6;
          break;
        case cOrthoCTRL:
          mode += 12;
          break;
        case (cOrthoCTRL + cOrthoSHIFT):
          mode += 18;
          break;
        case cOrthoALT:
          mode += 24;
          break;
        case (cOrthoALT + cOrthoSHIFT):
          mode += 30;
          break;
        case (cOrthoALT + cOrthoCTRL):
          mode += 36;
          break;
        case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
          mode += 42;
          break;
      }
      mod = 0;
      break;
  }
  switch (mod) {
    case 0:
      break;
    case cOrthoSHIFT:
      mode += 3;
      break;
    case cOrthoCTRL:
      mode += 6;
      break;
    case (cOrthoCTRL + cOrthoSHIFT):
      mode += 9;
      break;
    case cOrthoALT:
      mode += 68;
      break;
    case (cOrthoALT + cOrthoSHIFT):
      mode += 71;
      break;
    case (cOrthoALT + cOrthoCTRL):
      mode += 74;
      break;
    case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
      mode += 77;
      break;
  }
  return (I->Mode[mode]);
}
