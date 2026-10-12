/**
 * Dashed-line reps of a DistSet shared by angles and dihedrals: V holds N
 * dash endpoints (in pairs), drawn as lines or cylinders and ray traced as
 * sausages or flat-capped cylinders.
 */

#pragma once

#include "pymol/base/vla.h"
#include "pymol/reps/rep.h"

struct CGO;
struct DistSet;

struct RepDashLines : Rep {
  RepDashLines(pymol::CObject* obj, int state, int color_setting)
      : Rep(obj, state), color_setting(color_setting) {}
  ~RepDashLines() override;

  void render(RenderInfo* info) override;

  pymol::vla<float> V;
  int N = 0;
  DistSet* ds = nullptr;
  float linewidth, radius;
  CGO* shaderCGO = nullptr;
  int color_setting;  // cSetting_angle_color or cSetting_dihedral_color
};
