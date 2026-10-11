
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
#include "pymol/graphics/immediate_helper.h"

#include "pymol/app/main.h"
#include "pymol/base/feedback.h"
#include "pymol/reps/rep.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/objects/coord_set.h"
#include "pymol/python/p.h"
#include "pymol/base/util.h"
#include "pymol/scene/scene.h"
/**
 * Delete this instance and return a newly created instance.
 *
 * If creating the new instance fails (e.g. empty rep):
 * - Don't delete this instance
 * - Set cs->Active[rep] to false
 *
 * @todo Use move semantics instead? Return false on failure?
 */
Rep* Rep::rebuild() {
  assert(cs);
  assert(fNew);
  Rep* tmp = fNew(cs, getState());
  if (tmp) {
    tmp->fNew = fNew;
    delete this;
    return tmp;
  }

  // nothing returned -- visibility is zero...
  // keep the old object around, but inactive (TODO why?)
  cs->Active[type()] = false;
  return this;
}
/**
 * Rebuild if necessary (according to invalidation status). Returns either this
 * instance, or deletes this instance and returns a new one, or nullptr if the
 * rep became empty/inactive.
 */
Rep* Rep::update() {
  assert(cs);

  if (MaxInvalid == cRepInvNone) {
    return this;
  }

  auto I = this;
  auto const rep = type();
  auto const* cs_ = cs;

  assert(cs_->Active[rep]);

  if (MaxInvalid == cRepInvPick) {
    switch (rep) {
      case cRepLine:
      case cRepCyl:
      case cRepRibbon:
      case cRepNonbonded:
        // TODO Is this needed? How about other pickable reps, like:
        // - cartoon
        // - spheres
        // - surface with pick_surface=on
        MaxInvalid = cRepInvRep;
    }
  }

  if (MaxInvalid < cRepInvColor) {
    // nothing to do
  } else if (MaxInvalid == cRepInvColor) {
    I = recolor();
  } else if (MaxInvalid > cRepInvVisib || !sameVis()) {
    I = rebuild();
  } else if (!sameColor()) {
    I = recolor();
  }

  if (!cs_->Active[rep]) {
    delete I;
    return nullptr;
  }

  if (I) {
    I->MaxInvalid = cRepInvNone;
  }

  return I;
}
/**
 * Request that the rep gets updated. Update happens on next scene redraw.
 *
 * @param level Invalidate to at least this level
 */
void Rep::invalidate(cRepInv_t level) {
  auto I = this;
  SceneInvalidatePicking(
      I->G);  // for now, if anything invalidated, then invalidate picking
  if (level > I->MaxInvalid) I->MaxInvalid = level;
}

/**
 * Get the visRep mask according to auto_show_* settings
 */
cRepBitmask_t RepGetAutoShowMask(PyMOLGlobals* G) {
  cRepBitmask_t mask = 0;
  if (SettingGetGlobal_b(G, cSetting_auto_show_lines)) mask |= cRepLineBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_spheres)) mask |= cRepSphereBit;
  if (SettingGetGlobal_b(G, cSetting_auto_show_nonbonded))
    mask |= cRepNonbondedBit;
  return mask;
}
/**
 * Derived classes should overrride this method.
 *
 * TODO: What is this default implementation useful for? Debugging?
 */
void Rep::render(RenderInfo* info) {
  if (G->HaveGUI && G->ValidContext) {
    /* TODO */
  }
}
Rep::Rep(pymol::CObject* obj_, int state) : G(obj_->G), obj(obj_) {
  context.object = obj_;
  context.state = state;
}

Rep::Rep(CoordSet* cs_, int state) : Rep(cs_->Obj, state) { cs = cs_; }

Rep::~Rep() { FreeP(P); }

RepIterator::RepIterator(PyMOLGlobals* G, int rep_) {
  if (rep_ < 0) {
    end = cRepCnt;
    rep = -1;
  } else {
    end = rep_ + 1;
    rep = rep_ - 1;
  }
}
