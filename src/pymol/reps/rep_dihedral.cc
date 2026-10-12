
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
#include "pymol/reps/rep_dihedral.h"
#include "pymol/reps/rep_dash_internal.h"
#include "pymol/objects/object_dist.h"
#include "pymol/color/color.h"
#include "pymol/scene/scene.h"
#include "pymol/app/main.h"
#include "pymol/geometry/vector.h"
#include "pymol/settings/setting.h"
#include "pymol/objects/pymol_object.h"
#include "pymol/graphics/shader_mgr.h"
#include "pymol/graphics/cgo.h"
#include "pymol/objects/coord_set.h"

struct RepDihedral : RepDashLines {
  RepDihedral(pymol::CObject* obj, int state)
      : RepDashLines(obj, state, cSetting_dihedral_color) {}
  cRep_t type() const override { return cRepDihedral; }
};

Rep* RepDihedralNew(DistSet* ds, int state) {
  PyMOLGlobals* G = ds->G;
  int a;
  int n;
  float* v;
  float dash_len, dash_gap, dash_sum;
  int ok = true;

  if (!ok || !ds->NDihedralIndex) {
    return (nullptr);
  }

  auto I = new RepDihedral(ds->Obj, state);

  if (ds && ds->getNRep() > cRepDihedral && ds->Rep[cRepDihedral])
    I->cs = ds->Rep[cRepDihedral]->cs;

  dash_len =
      SettingGet_f(G, nullptr, ds->Obj->Setting.get(), cSetting_dash_length);
  dash_gap =
      SettingGet_f(G, nullptr, ds->Obj->Setting.get(), cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if (dash_sum < R_SMALL4) dash_sum = 0.5;

  I->ds = ds;

  n = 0;
  if (ds->NDihedralIndex) {
    float *v1, *v2, *v3, *v4, *v5;

    float d12[3], d32[3], d43[3], n12[3], n32[3], n43[3];
    float p12[3], p43[3], np12[3], np43[3], v12[3], v43[3];
    float s12[3], s43[3];

    float a32[3];
    float l1, l2;
    float d3[3], n1[3], n3[3], x[3], y[3];
    float radius, length, angle, phase, pos;
    float dihedral_size = SettingGet_f(G, nullptr, ds->Obj->Setting.get(),
                                       cSetting_dihedral_size);

    I->V.resize(ds->NDihedralIndex * 10);
    CHECKOK(ok, I->V);
    for (a = 0; ok && a < ds->NDihedralIndex; a = a + 6) {
      v1 = ds->DihedralCoord + 3 * a;
      v2 = v1 + 3;
      v3 = v1 + 6;
      v4 = v1 + 9;
      v5 = v1 + 12;

      subtract3f(v1, v2, d12);
      subtract3f(v3, v2, d32);
      subtract3f(v4, v3, d43);

      normalize23f(d12, n12);
      normalize23f(d32, n32);
      normalize23f(d43, n43);

      remove_component3f(d12, n32, p12);
      remove_component3f(d43, n32, p43);

      average3f(v2, v3, a32);

      l1 = (float)length3f(p12);
      l2 = (float)length3f(p43);

      if (l1 > l2)
        radius = l2;
      else
        radius = l1;
      radius *= dihedral_size;

      normalize23f(p12, np12);
      normalize23f(p43, np43);

      scale3f(np12, radius, v12);
      scale3f(np43, radius, v43);

      extrapolate3f(v12, n12, s12);
      add3f(s12, v2, s12);
      extrapolate3f(v43, n43, s43);
      add3f(s43, v3, s43);

      add3f(a32, v12, v12);
      add3f(a32, v43, v43);

      angle = get_angle3f(p12, p43);

      normalize23f(p12, n1);

      remove_component3f(p43, n1, d3);

      if (length3f(d3) < R_SMALL8) {
        d3[0] = 1.0F;
        d3[1] = 0.0F;
        d3[2] = 0.0F;
      } else {
        normalize23f(d3, n3);
      }

      scale3f(n1, radius, x);
      scale3f(n3, radius, y);

      I->V.check((n * 3) + 5);
      CHECKOK(ok, I->V);
      if (ok) {
        v = I->V + n * 3;
        copy3f(v12, v);
        v += 3;
        copy3f(a32, v);
        n += 2;
      }
      if (ok) I->V.check((n * 3) + 5);
      CHECKOK(ok, I->V);
      if (ok) {
        v = I->V + n * 3;
        copy3f(v43, v);
        v += 3;
        copy3f(a32, v);
        n += 2;
      }

      if (ok && v5[0] != 0.0F) { /* line 1 flag */

        I->V.check((n * 3) + 5);
        CHECKOK(ok, I->V);
        if (ok) {
          v = I->V + n * 3;
          copy3f(v1, v);
          v += 3;
          copy3f(v2, v);
          n += 2;
        }
      }

      if (ok && v5[1] != 0.0F) { /* line 2 flag */
        I->V.check((n * 3) + 5);
        CHECKOK(ok, I->V);
        if (ok) {
          v = I->V + n * 3;
          copy3f(v3, v);
          v += 3;
          copy3f(v2, v);
          n += 2;
        }
      }

      if (ok && v5[2] != 0.0F) { /* line 3 flag */
        I->V.check((n * 3) + 5);
        CHECKOK(ok, I->V);
        if (ok) {
          v = I->V + n * 3;
          copy3f(v3, v);
          v += 3;
          copy3f(v4, v);
          n += 2;
        }
      }

      /* now we have a relevant orthogonal axes */

      length = (float)(angle * radius * 2);

      /* figure out dash/gap phasing that will lead to nicely space dashes and
       * gaps */

      phase = dash_sum - (float)fmod(length / 2 + (dash_gap / 2), dash_sum);
      pos = -phase;

      if (length > R_SMALL4) {
        float vx[3], vy[3];
        float cur_angle;
        float cons_pos1, cons_pos2;

        while (ok && pos < length) {
          I->V.check((n * 3) + 5);
          CHECKOK(ok, I->V);
          if (ok) {
            cons_pos1 = pos;
            if (cons_pos1 < 0.0F) cons_pos1 = 0.0F;
            cons_pos2 = pos + dash_len;
            if (cons_pos2 > length) cons_pos2 = length;
          }
          if (ok && cons_pos1 < cons_pos2) {
            cur_angle = angle * cons_pos1 / length;

            v = I->V + n * 3;
            scale3f(x, (float)cos(cur_angle), vx);
            scale3f(y, (float)sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(a32, v, v);

            cur_angle = angle * cons_pos2 / length;

            v += 3;
            scale3f(x, (float)cos(cur_angle), vx);
            scale3f(y, (float)sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(a32, v, v);

            n += 2;
          }
          pos += dash_sum;
        }
      }
    }
    if (ok) I->V.resize(n * 3);
    CHECKOK(ok, I->V);
    I->N = n;
  }
  if (!ok) {
    delete I;
    I = nullptr;
  }
  return (Rep*)I;
}
