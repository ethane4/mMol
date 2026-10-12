
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

#include "pymol/app/main.h"
#include "pymol/base/base.h"
#include "pymol/ui/control.h"
#include "pymol/scene/scene.h"
#include "pymol/executive/executive.h"
#include "pymol/settings/setting.h"
#include "pymol/python/p.h"
#include "pymol/ui/seq.h"
#include "pymol/scene/movie.h"
#include "pymol/text/text.h"
#include "pymol/base/util.h"
#include "pymol/ui/ortho.h"
#include "pymol/graphics/cgo.h"
#include "pymol/graphics/immediate_helper.h"

#define SDOF_QUEUE_MASK 0x1F

struct CControl : public Block {
  /* not saved */

  int sdofActive{};
  double sdofLastIterTime{};
  int sdofMode{};
  float sdofTrans[3]{};
  float sdofRot[3]{};
  unsigned int sdofWroteTo{}, sdofReadFrom{}; /* queue synchronization fields */
  float sdofBuffer[(SDOF_QUEUE_MASK + 1) * 6]{};

  CControl(PyMOLGlobals* G) : Block(G) {};
};

int ControlSdofUpdate(PyMOLGlobals* G, float tx, float ty, float tz, float rx,
                      float ry, float rz) {
  /* may be called asynchronously anytime after CControl has been initialized */

  CControl* I = G->Control;
  if (I) {
    float tol = 0.0001;
    bool active = fabs(tx) >= tol || fabs(ty) >= tol || fabs(tz) >= tol ||
                  fabs(rx) >= tol || fabs(ry) >= tol || fabs(rz) >= tol;

    if (active) {
      int slot = (I->sdofWroteTo + 1) & SDOF_QUEUE_MASK;
      float* buffer = I->sdofBuffer + (6 * slot);

      buffer[0] = tx;
      buffer[1] = ty;
      buffer[2] = tz;
      buffer[3] = rx;
      buffer[4] = ry;
      buffer[5] = rz;

      I->sdofWroteTo = slot;

      {
        if (!I->sdofActive) {
          I->sdofLastIterTime = UtilGetSeconds(G);
        }
      }
      /*printf("wrote %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d %d
       * %d\n",tx,ty,tz,rx,ry,rz,I->sdofReadFrom,I->sdofWroteTo,slot); */
    }
    I->sdofActive = active;
  }
  return 1;
}

int ControlSdofIterate(PyMOLGlobals* G) {
  CControl* I = G->Control;
  if (I->sdofWroteTo != I->sdofReadFrom && I->sdofActive) {
    int slot = I->sdofWroteTo;

    /* new data available */

    float* buffer = I->sdofBuffer + (6 * slot);

    I->sdofTrans[0] = buffer[0];
    I->sdofTrans[1] = buffer[1];
    I->sdofTrans[2] = buffer[2];
    I->sdofRot[0] = buffer[3];
    I->sdofRot[1] = buffer[4];
    I->sdofRot[2] = buffer[5];

    I->sdofReadFrom = slot;

    double now = UtilGetSeconds(G);
    double delta = now - I->sdofLastIterTime;
    I->sdofLastIterTime = now;

    {
      /* suppress small amounts of combined motion using a truncated switching
       * function */

      float len_rot = length3f(I->sdofRot);
      float len_trans = length3f(I->sdofTrans);
      float *dom, *sub;
      if (len_rot > len_trans) {
        dom = &len_rot;
        sub = &len_trans;
      } else {
        dom = &len_trans;
        sub = &len_rot;
      }

      {
        float expo = 2.0F;
        float trnc = 0.05F;
        float sub_fxn = (*sub) / (*dom);
        if (sub_fxn < trnc) {
          sub_fxn = 0.0F;
        } else if (sub_fxn < 0.5F) {
          sub_fxn = (sub_fxn - trnc) / (0.5F - trnc);
          sub_fxn = (float)pow(sub_fxn, expo);
        } else {
          sub_fxn = 1.0F - (float)pow(1.0F - sub_fxn, expo);
        }
        *dom = 1.0F;
        *sub = sub_fxn;
        scale3f(I->sdofTrans, len_trans, I->sdofTrans);
        scale3f(I->sdofRot, len_rot, I->sdofRot);
      }
    }
    /* translate */
    SceneTranslateScaled(G, delta * I->sdofTrans[0], -delta * I->sdofTrans[1],
                         -delta * I->sdofTrans[2], I->sdofMode);

    /* rotate */
    SceneRotateScaled(G, 2.0F * delta * I->sdofRot[0],
                      -2.0F * delta * I->sdofRot[1],
                      -2.0F * delta * I->sdofRot[2], I->sdofMode);

    SceneDirty(G);
  }
  return 1;
}

int ControlRocking(PyMOLGlobals* G) {
  if (G->Interrupt) {
    SettingSetGlobal_b(G, cSetting_rock, false);
  }
  return SettingGetGlobal_b(G, cSetting_rock);
}

Block* ControlGetBlock(PyMOLGlobals* G) {
  CControl* I = G->Control;
  {
    return (I);
  }
}
int ControlIdling(PyMOLGlobals* G) {
  CControl* I = G->Control;
  return (I->sdofActive || MoviePlaying(G) ||
          SettingGetGlobal_b(G, cSetting_rock) ||
          SettingGetGlobal_b(G, cSetting_sculpting));
}
void ControlFree(PyMOLGlobals* G) { DeleteP(G->Control); }
pymol::Result<bool> ControlRock(PyMOLGlobals* G, int mode) {
  switch (mode) {
    case -2:
      break;
    case -1:
      SettingSetGlobal_b(G, cSetting_rock,
                         !SettingGetGlobal_b(G, cSetting_rock));
      if (SettingGetGlobal_b(G, cSetting_rock)) {
        SceneRestartSweepTimer(G);
      }
      break;
    case 0:
      SettingSetGlobal_b(G, cSetting_rock, false);
      break;
    case 1:
      SettingSetGlobal_b(G, cSetting_rock, true);
      SceneRestartSweepTimer(G);
      break;
  }
  if (mode != -2) {
    SceneRestartFrameTimer(G);
    OrthoDirty(G);
  }
  return SettingGetGlobal_b(G, cSetting_rock);
}
int ControlInit(PyMOLGlobals* G) {
  CControl* I = nullptr;
  if ((I = (G->Control = new CControl(G)))) {
    I->active = true;
    I->TextColor[0] = 1.0;
    I->TextColor[1] = 0.75;
    I->TextColor[2] = 0.75;
    OrthoAttach(G, I, cOrthoTool);
    return 1;
  } else
    return 0;
}
