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
-* Cameron Mura
-*
-*
Z* -------------------------------------------------------------------
*/

#include "pymol/base/err.h"
#include "pymol/reps/rep_cartoon.h"
#include "pymol/base/feedback.h"

#include "pymol/select/atom_iterators.h"
#include "pymol/reps/rep_cartoon_internal.h"

enum {
  CARTOON_CYLINDRICAL_HELICES_CURVED = 1,
  CARTOON_CYLINDRICAL_HELICES_STRAIGHT = 2
};

#include "pymol/executive/executive.h"

RepCartoon::~RepCartoon() {
  auto I = this;
  assert(I->ray != I->preshader);
  CGOFree(I->preshader);
  CGOFree(I->ray);
  CGOFree(I->std);
  FreeP(I->LastVisib);
}

/**
 * CGOAddTwoSidedBackfaceSpecialOps: this function takes in a CGO,
 * and outputs a CGO with that CGO wrapped with the two operations:
 * enabling and disabling back faces if two_sided_lighting is not set
 * (i.e., the CGOSpecial operations, see below)
 *
 */
static CGO* CGOAddTwoSidedBackfaceSpecialOps(CGO* cgo) {
  if (!CGOHasOperations(cgo)) {
    return cgo;
  }

  auto G = cgo->G;
  CGO* tmpCGO = CGONew(G);
  CGOSpecial(tmpCGO, ENABLE_BACK_FACES_IF_NOT_TWO_SIDED);
  tmpCGO->free_append(cgo);
  CGOSpecial(tmpCGO, DISABLE_BACK_FACES_IF_NOT_TWO_SIDED);
  CGOStop(tmpCGO);
  return tmpCGO;
}

static int RepCartoonCGOGenerate(RepCartoon* I, RenderInfo* info) {
  PyMOLGlobals* G = I->G;
  int ok = true;

  int use_shaders;
  float alpha =
      1.0F - SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(),
                          cSetting_cartoon_transparency);

  bool const hasAlpha = alpha < 1 || [](RepCartoon* I) {
    for (CoordSetAtomIterator iter(I->cs); iter.next();) {
      auto ai = iter.getAtomInfo();
      if ((ai->visRep & cRepCartoonBit) &&
          AtomSettingGetWD(I->G, ai, cSetting_cartoon_transparency, 0.0f) > 0) {
        return true;
      }
    }
    return false;
  }(I);

  I->setHasTransparency(hasAlpha);

  use_shaders = SettingGetGlobal_b(G, cSetting_use_shaders) &&
                SettingGetGlobal_b(G, cSetting_cartoon_use_shader);

  assert(!I->std);

  if (use_shaders) {
    if (hasAlpha && (SettingGetGlobal_i(G, cSetting_transparency_mode) != 3)) {
      // some transparency
      std::unique_ptr<CGO> simplified(CGOSimplify(I->preshader));
      std::unique_ptr<CGO> optimized(
          CGOOptimizeToVBOIndexedWithColorEmbedTransparentInfo(
              simplified.get(), 0, nullptr, true));

      CGO* tmp2CGO = CGONew(G);
      CGOEnable(tmp2CGO, GL_BACK_FACE_CULLING);
      tmp2CGO->move_append(std::move(*optimized));
      CGODisable(tmp2CGO, GL_BACK_FACE_CULLING);
      CGOStop(tmp2CGO);

      I->std = tmp2CGO;
    } else {
      CGO* leftOverCGOBorrowed = nullptr;
      std::unique_ptr<CGO> leftOverCGOManaged;
      std::unique_ptr<CGO> convertcgo(CGONew(G));

      leftOverCGOBorrowed = I->preshader;

      bool has_spheres_to_optimize =
          CGOHasSphereOperations(leftOverCGOBorrowed);
      if (has_spheres_to_optimize) {
        /* Optimize spheres and putting them into convertcgo as a shader CGO
         * operation */
        std::unique_ptr<CGO> leftOverAfterSpheresCGO(CGONew(G));
        std::unique_ptr<CGO> sphereVBOs(CGOOptimizeSpheresToVBONonIndexed(
            leftOverCGOBorrowed, 0, true, leftOverAfterSpheresCGO.get()));
        if (sphereVBOs) {
          convertcgo->move_append(std::move(*sphereVBOs));
          leftOverCGOManaged = std::move(leftOverAfterSpheresCGO);
          leftOverCGOBorrowed = leftOverCGOManaged.get();
        }
      }

      /* For the rest of the primitives that exist, simplify them into Geometry
       * (should probably be no more, but do this anyway) */
      std::unique_ptr<CGO> leftOverCGOSimplified(
          CGOSimplify(leftOverCGOBorrowed));
      if (leftOverCGOSimplified) {
        std::unique_ptr<CGO> optimized(
            CGOOptimizeToVBONotIndexed(leftOverCGOSimplified.get()));
        if (optimized) {
          convertcgo->move_append(std::move(*optimized));
        }
      }
      I->std = CGOAddTwoSidedBackfaceSpecialOps(convertcgo.release());
    }
    I->std->use_shader = true;
  } else {
    if (ok) {
      auto simplifiedCGO = CGOSimplify(I->preshader, 0);
      if (hasAlpha) {
        auto convertedCGO = CGOConvertTrianglesToAlpha(simplifiedCGO);
        CGOFree(simplifiedCGO);
        I->std = convertedCGO;
        if (I->std) {
          I->std->render_alpha =
              1;  // CGORenderGL should call CGOSetZVector/CGORenderGLAlpha only
        }
      } else {
        I->std = simplifiedCGO;
        CHECKOK(ok, I->std);
      }
      if (I->std) {
        I->std = CGOAddTwoSidedBackfaceSpecialOps(I->std);
      }
    }
  }
  I->disposePreshaderCGO();

  return ok;
}

void RepCartoon::render(RenderInfo* info) {
  auto I = this;

  if (info->ray) {
    CGO* raycgo = I->ray ? I->ray : I->preshader;

    if (raycgo && !CGORenderRay(raycgo, info->ray, info, nullptr, nullptr,
                                I->cs->Setting.get(), I->obj->Setting.get())) {
      PRINTFB(G, FB_RepCartoon, FB_Warnings)
      " %s-Warning: ray rendering failed\n", __func__ ENDFB(G);
      CGOFree(I->ray);
    }
    return;
  }

  if (G->HaveGUI && G->ValidContext) {
    if (I->preshader) {
      assert(!I->std);

      bool ok = RepCartoonCGOGenerate(I, info);

      if (!ok) {
        I->disposePreshaderCGO();
        I->invalidate(cRepInvPurge);
        I->cs->Active[cRepCartoon] = false;
      }
    }

    // Generation is what discovers transparency, and it runs in whichever pass
    // reaches the rep first (opaque). Drawing transparent geometry there loses
    // the OIT look for that one frame, which strobes while a transparency
    // slider is dragged; CoordSet draws it in the transparent pass instead.
    if (info->pass == RenderPass::Opaque && I->hasTransparency()) return;

    if (I->std && CGOHasOperations(I->std)) {
      assert(!I->preshader);

      CGORender(I->std, nullptr, I->cs->Setting.get(), I->obj->Setting.get(),
                info, I);
    }
  }
}

static void GenerateRepCartoonDrawDebugLineAlongPath(CGO* cgo, int nAt,
                                                     float* pv) {
  float *v, *v1, *v2, *v3, *v4;
  float t0[3], t1[3];
  int a;
  CGOColor(cgo, 1.0, 1.0, 1.0);
  v1 = nullptr;
  v2 = nullptr;
  v3 = nullptr;
  v4 = nullptr;
  v = pv;
  if (nAt > 1) {
    CGOBegin(cgo, GL_LINE_STRIP);
    for (a = 0; a < nAt; a++) {
      v4 = v3;
      v3 = v2;
      v2 = v1;
      v1 = v;
      if (v1 && v2 && v3 && v4) {
        add3f(v1, v4, t0);
        add3f(v2, v3, t1);
        /*            scale3f(t0,0.2024,t0);
                      scale3f(t1,0.2727,t1); */
        scale3f(t0, 0.2130F, t0);
        scale3f(t1, 0.2870F, t1);

        add3f(t0, t1, t0);
        CGOVertexv(cgo, t0);
      }
      v += 3;
    }
    CGOEnd(cgo);
  }
}

static int GenerateRepCartoonDrawDebugNormals(CGO* cgo, CExtrude* ex, int n_p) {
  int b, ok;
  float *v, *vn;
  float t0[3];
  ok = CGOColor(cgo, 0.0, 1.0, 0.0);
  v = ex->p;
  vn = ex->n + 3;
  if (ok) ok &= CGOBegin(cgo, GL_LINES);
  for (b = 0; ok && b < n_p; b++) {
    ok &= CGOVertexv(cgo, v);
    add3f(v, vn, t0);
    if (ok) ok &= CGOVertexv(cgo, t0);
    v += 3;
    vn += 9;
  }
  if (ok) ok &= CGOEnd(cgo);
  return ok;
}

static int GenerateRepCartoonDrawDebugOrient(CGO* cgo, int nAt, float* pv,
                                             float* pvo, float* tv) {
  int ok, a;
  float *v1, *v2, *v3;
  float t0[3];
  ok = CGOColor(cgo, 1.0, 1.0, 1.0);
  if (ok) ok &= CGOBegin(cgo, GL_LINES);
  v1 = pv;
  v2 = pvo;
  v3 = tv;
  for (a = 0; ok && a < nAt; a++) {
    ok &= CGOVertexv(cgo, v1);
    if (ok) {
      add3f(v1, v2, t0);
      add3f(v2, t0, t0);
      ok &= CGOVertexv(cgo, t0);
    }
    if (ok) {
      subtract3f(v1, v3, t0);
      ok &= CGOVertexv(cgo, t0);
    }
    if (ok) {
      add3f(v1, v3, t0);
      ok &= CGOVertexv(cgo, t0);
    }
    v1 += 3;
    v2 += 3;
    v3 += 3;
  }
  if (ok) ok &= CGOEnd(cgo);
  return ok;
}

static int GenerateRepCartoonProcessCylindricalHelices(
    PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs, CGO* cgo, CExtrude* ex,
    int nAt, int* seg, float* pv, float* tv, float* pvo, const CCInOut* cc,
    int* at, float* dl, int cartoon_color, int discrete_colors,
    float loop_radius, const float objAlpha) {
  int ok = true;
  int n_p, n_pm1, n_pm2;
  const float* v0;
  float *v, *v1, *v2, *vo, *d;
  float* valpha;
  float* vc = nullptr;
  int atom_index1, atom_index2, *s, *atp, a, cur_car;
  unsigned* vi;
  int last_color, uniform_color;
  bool hasAtomLevelTrans = false;
  int contFlag, extrudeFlag;
  int b, c1, c2;
  float *h_start = nullptr, *h_end = nullptr;
  float t0[3], t1[3], t2[3], t3[3];
  float t4[3];
  float helix_radius;
  CGOPickColor(cgo, 0, cPickableNoPick);
  helix_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_cartoon_helix_radius);

  /* this is confusing because we're borrowing Extrude's arrays
   * for convenient storage, but not actually calling Extrude */
  n_p = 0;
  v = ex->p;
  vc = ex->c;
  valpha = ex->alpha;
  vi = ex->i;

  last_color = -1;
  uniform_color = true;

  v1 = pv; /* points */
  v2 = tv; /* tangents */
  vo = pvo;
  d = dl;
  s = seg;
  atp = at; /* cs index pointer */
  a = 0;
  contFlag = true;
  cur_car = cCartoon_skip;
  extrudeFlag = false;

  while (contFlag) {
    if ((*cc) != cur_car) { /* new cartoon type */
      if (n_p) {            /* any cartoon points? */
        extrudeFlag = true;
      } else {
        cur_car = *(cc); /* now: go ahead and switch cartoons */
        n_p = 0;
        v = ex->p;
        vc = ex->c;
        valpha = ex->alpha;
        vi = ex->i;
        last_color = -1;
        uniform_color = true;
      }
    }
    if (a && !extrudeFlag) {
      if ((*s) != *(s - 1)) { /* new segment */
        if (n_p) {            /* any cartoon points? */
          extrudeFlag = true;
        } else {
          n_p = 0;
          v = ex->p;
          vc = ex->c;
          valpha = ex->alpha;
          vi = ex->i;
          last_color = -1;
          uniform_color = true;
        }
      }
    }
    if (!extrudeFlag) {
      if ((a < (nAt - 1)) &&
          (*s == *(s + 1))) { /* working in the same segment... */
        AtomInfoType *ai1, *ai2;
        atom_index1 = cs->IdxToAtm[*atp];
        atom_index2 = cs->IdxToAtm[*(atp + 1)];
        ai1 = obj->AtomInfo + atom_index1;
        ai2 = obj->AtomInfo + atom_index2;

        c1 = AtomSettingGetWD(G, ai1, cSetting_cartoon_color, cartoon_color);
        c2 = AtomSettingGetWD(G, ai2, cSetting_cartoon_color, cartoon_color);

        float alpha1 =
            1.0f - AtomSettingGetWD(G, ai1, cSetting_cartoon_transparency,
                                    1.0f - objAlpha);
        float alpha2 =
            1.0f - AtomSettingGetWD(G, ai2, cSetting_cartoon_transparency,
                                    1.0f - objAlpha);
        if (!hasAtomLevelTrans && (alpha1 != objAlpha || alpha2 != objAlpha)) {
          hasAtomLevelTrans = true;
        }

        if (c1 < 0) c1 = ai1->color;
        if (c2 < 0) c2 = ai2->color;

        if ((*(cc) == *(cc + 1)) && (c1 != c2)) uniform_color = false;
        if (last_color >= 0) {
          if (c1 != last_color) uniform_color = false;
        }
        last_color = c1;

        v0 = ColorGet(G, c1);
        *(vc++) = *(v0++);
        *(vc++) = *(v0++);
        *(vc++) = *(v0++);
        *(valpha++) = alpha1;
        *(vi++) = ai1->masked ? -1 : atom_index1;

        v0 = ColorGet(G, c2); /* kludge */
        *(vc) = *(v0++);
        *(vc + 1) = *(v0++);
        *(vc + 2) = *(v0++);
        *(valpha) = alpha2;
        *(vi) = ai2->masked ? -1 : atom_index2;
      } else {
        vc += 3; /* part of kludge */
        valpha++;
        vi++;
      }
      if (cur_car == cCartoon_skip_helix) {
        if (!n_p) {
          h_start = v1;
          h_end = v1;
        } else {
          h_end = v1;
        }
        copy3f(
            v1,
            v); /* just store coordinates until we have a complete cylinder */
        v += 3;
        n_p++;
      }
      v1 += 3;
      v2 += 3;
      vo += 3;
      d++;
      atp++;
      s++;
      cc++;
    }
    a++;
    if (a == nAt) {
      contFlag = false;
      if (n_p) extrudeFlag = true;
    }
    if (extrudeFlag) { /* generate cylinder */
      if (n_p > 1) {
        atom_index1 = cs->IdxToAtm[*(atp - 1)];

        c1 = AtomSettingGetWD(G, obj->AtomInfo + atom_index1,
                              cSetting_cartoon_color, cartoon_color);

        auto ai1 = obj->AtomInfo + atom_index1;
        auto ai2 = obj->AtomInfo + atom_index2;

        float alpha1 =
            1.0f - AtomSettingGetWD(G, ai1, cSetting_cartoon_transparency,
                                    1.0f - objAlpha);
        float alpha2 =
            1.0f - AtomSettingGetWD(G, ai2, cSetting_cartoon_transparency,
                                    1.0f - objAlpha);
        if (!hasAtomLevelTrans && (alpha1 != objAlpha || alpha2 != objAlpha)) {
          hasAtomLevelTrans = true;
        }
        if (c1 < 0) c1 = (obj->AtomInfo + atom_index1)->color;

        if (n_p < 5) {
          copy3f(ex->p, t3);
          copy3f(v - 3, t4);
        } else {
          add3f(ex->p, ex->p + 9, t0);
          add3f(ex->p + 3, ex->p + 6, t1);
          scale3f(t0, 0.2130F, t0);
          scale3f(t1, 0.2870F, t1);
          add3f(t0, t1, t3);
          add3f(v - 3, v - 12, t0);
          add3f(v - 6, v - 9, t1);
          scale3f(t0, 0.2130F, t0);
          scale3f(t1, 0.2870F, t1);
          add3f(t0, t1, t4);

          /* extend helix to line up with CA */
          subtract3f(t4, t3, t0);
          normalize3f(t0);
          subtract3f(v - 3, t3, t1);
          project3f(t1, t0, t4);
          add3f(t3, t4, t4);
          invert3f(t0);
          subtract3f(ex->p, t4, t1);
          project3f(t1, t0, t3);
          add3f(t3, t4, t3);

          /* relocate terminal CA to touch helix, if necessary */
          if (h_start && h_end) {
            float f0;
            subtract3f(h_start, t3, t0);
            f0 = helix_radius - loop_radius * 2;
            if (length3f(t0) > f0) {
              normalize3f(t0);
              scale3f(t0, f0, t1);
              add3f(t1, t3, h_start);
            }
            subtract3f(h_end, t4, t0);
            if (length3f(t0) > f0) {
              normalize3f(t0);
              scale3f(t0, f0, t1);
              add3f(t1, t4, h_end);
            }
          }
        }
        /* push helix out a tad to consume loop */
        subtract3f(t4, t3, t0);
        normalize3f(t0);
        scale3f(t0, loop_radius * 2, t0);
        add3f(t0, t4, t4);
        invert3f(t0);
        add3f(t0, t3, t3);
        if (uniform_color && !hasAtomLevelTrans) {
          cgo->add<cgo::draw::cylinder>(t3, t4, helix_radius, ex->c, ex->c);
        } else {
          subtract3f(t4, t3, t0);
          n_pm1 = n_p - 1;
          n_pm2 = n_p - 2;
          for (b = 0; ok && b < n_pm1; b++) {
            scale3f(t0, ((float)b) / n_pm1, t1);
            scale3f(t0, ((float)b + 1) / n_pm1, t2);

            add3f(t3, t1, t1);
            add3f(t3, t2, t2);
            if (hasAtomLevelTrans) {
              cgo->add<cgo::draw::custom_cylinder_alpha>(
                  t1, t2, helix_radius, ex->c + (b * 3), ex->c + (b + 1) * 3,
                  ex->alpha[b], ex->alpha[b + 1],
                  (float)(b ? cCylCap::None : cCylCap::Flat),
                  (float)(b == n_pm2 ? cCylCap::Flat : cCylCap::None));
            } else {
              cgo->add<cgo::draw::custom_cylinder>(
                  t1, t2, helix_radius, ex->c + (b * 3), ex->c + (b + 1) * 3,
                  (float)(b ? cCylCap::None : cCylCap::Flat),
                  (float)(b == n_pm2 ? cCylCap::Flat : cCylCap::None));
            }
          }
        }
      }
      a--; /* undo above... */
      extrudeFlag = false;
      n_p = 0;
      v = ex->p;
      vc = ex->c;
      valpha = ex->alpha;
      vi = ex->i;
      uniform_color = true;
      last_color = -1;
    }
  }
  return ok;
}

/**
 * skip > dash > loop
 */
inline int prioritize(int a, int b) {
  if (a == b) return a;
  if (a == cCartoon_skip || b == cCartoon_skip) return cCartoon_skip;
  if (a == cCartoon_dash || b == cCartoon_dash) return cCartoon_dash;
  return cCartoon_loop;
}

/* CheckExtrudeContigFlags: compares current cartoon type and next to figure out
   whether to extrude, truncate, and/or whether its contiguous */
static int CheckExtrudeContigFlags(int nAt, int n_p, int a,
                                   int* cur_car,  // in-out
                                   const CCInOut* cc, const int* segptr,
                                   int* /* contigFlag */, int* extrudeFlag) {
  int truncate = false;
  int next_car = cCartoon_skip;

  if (a < (nAt - 1) && (*segptr) == *(segptr + 1)) {
    next_car = prioritize(cc[0].getCCOut(), cc[1].getCCIn());
  }

  if (*cur_car != next_car) { /* new cartoon type */
    if (n_p) { /* any cartoon points? then extrude, otherwise truncate */
      *extrudeFlag = true;
    } else {
      *cur_car = next_car; /* no: go ahead and switch cartoons */
      truncate = true;
    }
  }

  return truncate;
}

static void ComputeCartoonAtomColors(PyMOLGlobals* G, ObjectMolecule* obj,
                                     CoordSet* cs, int* nuc_flag,
                                     int atom_index1, int atom_index2, int* c1a,
                                     int* c2a, int* atp, const CCInOut* cc,
                                     int cur_car, int cartoon_color,
                                     float& alpha1, float& alpha2,
                                     int nucleic_color, int discrete_colors,
                                     int n_p, int contigFlag) {
  int c1, c2;

  if (nucleic_color >= 0 && (nuc_flag[*atp] || nuc_flag[*(atp + 1)])) {
    c1 = c2 = nucleic_color;
  } else {
    c1 = c2 = cartoon_color;
  }

  auto ai1 = obj->AtomInfo + atom_index1;
  auto ai2 = obj->AtomInfo + atom_index2;

  AtomSettingGetIfDefined(G, ai1, cSetting_cartoon_color, &c1);
  AtomSettingGetIfDefined(G, ai2, cSetting_cartoon_color, &c2);
  alpha1 = 1.0f - AtomSettingGetWD(G, ai1, cSetting_cartoon_transparency,
                                   1.0f - alpha1);
  alpha2 = 1.0f - AtomSettingGetWD(G, ai2, cSetting_cartoon_transparency,
                                   1.0f - alpha2);

  if (c1 < 0) c1 = ai1->color;
  if (c2 < 0) c2 = ai2->color;

  if (discrete_colors) {
    int next_car = *(cc + 1);
    // end of loop or dash segment
    if (cur_car != next_car) {
      if (cur_car == cCartoon_dash) {
        c2 = c1;
      } else if (next_car == cCartoon_dash) {
        c1 = c2;
      } else if (cur_car == cCartoon_loop) {
        c2 = c1;
      } else if (next_car == cCartoon_loop) {
        c1 = c2;
      }
    } else if (n_p == 0 && contigFlag &&
               (cur_car == cCartoon_dash || cur_car == cCartoon_loop)) {
      // beginning of loop or dash segment
      c1 = c2;
    }
  }
  *c1a = c1;
  *c2a = c2;
}

/**
 * Get cartoon quality setting, adapt to number of atoms if -1
 */
/* mMol: the atoms the scene draws as molecules, for the automatic quality
   tiers: every enabled molecule, times its states when all_states shows them
   together. Tiering on one coordinate set gave each of 60 split assembly
   copies (or 60 all_states states) full detail -- 3x the geometry of the same
   atoms in one object, and a GPU-bound viewport. */
static int SceneMoleculeAtoms(PyMOLGlobals* G, CoordSet* cs) {
  long long total = 0;
  if (auto objs = ExecutiveFindObjectsByType(G, cObjectMolecule)) {
    for (size_t i = 0, n = VLAGetSize(objs); i < n; ++i) {
      auto obj = static_cast<ObjectMolecule*>(objs[i]);
      if (!obj->Enabled && obj != cs->Obj) continue;
      bool all =
          SettingGet_b(G, obj->Setting.get(), nullptr, cSetting_all_states);
      total += (long long)obj->NAtom * (all ? std::max(obj->NCSet, 1) : 1);
    }
    VLAFreeP(objs);
  }
  return (int)std::min<long long>(std::max<long long>(total, cs->NIndex),
                                  1 << 30);
}

static int GetCartoonQuality(CoordSet* cs, int setting, int v1, int v2, int v3,
                             int v4, int min_ = 3) {
  int quality = SettingGet<int>(cs->G, cs->Setting.get(),
                                cs->Obj->Setting.get(), setting);

  if (quality == -1) {
    int natom = SceneMoleculeAtoms(cs->G, cs);
    quality = (natom < 100000)   ? v1
              : (natom < 500000) ? v2
              : (natom < 999999) ? v3
                                 : v4;
  } else if (quality < min_) {
    quality = min_;
  }

  return quality;
}

static CGO* GenerateRepCartoonCGO(
    CoordSet* cs, ObjectMolecule* obj, nuc_acid_data* ndata,
    short use_cylinders_for_strands, float* pv, int nAt, float* tv, float* pvo,
    float* dl, const CCInOut* car, int* seg, int* at, int* nuc_flag,
    float* putty_vals, float alpha,
    const std::unordered_map<int, CyclicSegmentInfo>* cyclicSegs = nullptr) {
  PyMOLGlobals* G = cs->G;
  int ok = true;
  CGO* cgo;
  int contigFlag, contFlag, extrudeFlag, n_p;
  CExtrude* ex = nullptr;
  float dev;
  unsigned int* vi;
  int atom_index1, atom_index2;
  float *v, *v1, *v2, *vo;
  float *d, *vc = nullptr, *vn;
  float* valpha = nullptr;
  int* atp;
  int c1, c2;
  int a;
  int sampling;
  float* sampling_tmp;
  int* segptr;
  const CCInOut* cc;
  int cur_car;
  float loop_radius;
  int nucleic_color = 0;
  float throw_;
  float power_a = 5;
  float power_b = 5;
  int refine;
  float tube_radius;
  float putty_radius;
  int cartoon_debug, cylindrical_helices, cartoon_color, highlight_color,
      discrete_colors, loop_quality, oval_quality, tube_quality, putty_quality;
  float length, width;
  float oval_width, oval_length;
  float dumbbell_radius, dumbbell_width, dumbbell_length;
  float ring_width;

  auto EXTRUDE_TRUNCATE = [&ex, &n_p, &v, &vc, &valpha, &vn, &vi]() {
    ExtrudeTruncate(ex, 0);
    n_p = 0;
    v = ex->p;
    vc = ex->c;
    valpha = ex->alpha;
    vn = ex->n;
    vi = ex->i;
  };

  cartoon_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_cartoon_color);
  ring_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_cartoon_ring_width);
  if (ring_width < 0.0F) {
    ring_width = fabs(SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_stick_radius)) *
                 0.5F;
  }
  length = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                        cSetting_cartoon_rect_length);
  width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                       cSetting_cartoon_rect_width);
  oval_length = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_oval_length);
  dumbbell_length = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_cartoon_dumbbell_length);
  oval_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_cartoon_oval_width);
  dumbbell_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_cartoon_dumbbell_width);
  dumbbell_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_cartoon_dumbbell_radius);
  if (dumbbell_radius < 0.01F) dumbbell_radius = 0.01F;

  auto tube_cap = static_cast<cCylCap>(SettingGet<int>(
      G, cs->Setting.get(), obj->Setting.get(), cSetting_cartoon_tube_cap));
  auto loop_cap = static_cast<cCylCap>(SettingGet<int>(
      G, cs->Setting.get(), obj->Setting.get(), cSetting_cartoon_loop_cap));

  tube_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_tube_radius);
  if (tube_radius < 0.01F) tube_radius = 0.01F;
  putty_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_cartoon_putty_radius);
  /* WLD removed: if(putty_radius<0.01F) putty_radius=0.01F; --
     should not constrain what is effectively a scale factor */

  tube_quality =
      GetCartoonQuality(cs, cSetting_cartoon_tube_quality, 9, 7, 6, 5);
  oval_quality =
      GetCartoonQuality(cs, cSetting_cartoon_oval_quality, 10, 8, 7, 6);
  putty_quality =
      GetCartoonQuality(cs, cSetting_cartoon_putty_quality, 11, 9, 7, 5);
  loop_quality =
      GetCartoonQuality(cs, cSetting_cartoon_loop_quality, 6, 6, 5, 4);
  sampling = GetCartoonQuality(cs, cSetting_cartoon_sampling, 7, 5, 3, 2, 1);
  // cartoon_spline pairs with denser sampling: the spline path has no
  // per-residue kinks to hide, but shadow-casting silhouettes (Metal RT) still
  // show the ring polyline at the classic 7/residue. Double the AUTO tiers
  // only; an explicit cartoon_sampling is respected as-is.
  const bool spline_mode = SettingGet_b(
      G, cs->Setting.get(), obj->Setting.get(), cSetting_cartoon_spline);
  if (spline_mode && SettingGet<int>(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_cartoon_sampling) < 0) {
    sampling *= 2;
    // Small structures can afford ChimeraX's full 20 rings/residue outright:
    // cast-shadow silhouettes of a 14-ring strand still show the ring polyline
    // at high zoom.
    if (SceneMoleculeAtoms(G, cs) < 100000 && sampling < 20) sampling = 20;
  }

  PRINTFB(G, FB_RepCartoon, FB_Blather)
  " RepCartoon: Use settings tube_quality=%d oval_quality=%d putty_quality=%d "
  "loop_quality=%d sampling=%d\n",
      tube_quality, oval_quality, putty_quality, loop_quality,
      sampling ENDFB(G);

  if (SettingGetGlobal_i(G, cSetting_ray_trace_mode) > 0)
    if (loop_quality < 12) loop_quality *= 2;
  refine = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                        cSetting_cartoon_refine);

  // cartoon_spline: one natural cubic spline per segment instead of the
  // per-interval blend below (see CartoonGenerateSampleSpline).
  CartoonSplineCoefs spline_coefs;
  if (spline_mode) {
    spline_coefs = CartoonBuildSplineCoefs(nAt, seg, pv, cyclicSegs);
  }
  power_a = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_cartoon_power);
  power_b = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_cartoon_power_b);
  throw_ = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                        cSetting_cartoon_throw);
  loop_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_loop_radius);
  if (loop_radius < 0.01F) loop_radius = 0.01F;
  discrete_colors = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_cartoon_discrete_colors);
  nucleic_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_cartoon_nucleic_acid_color);
  if (nucleic_color == -1) nucleic_color = cartoon_color;
  highlight_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_cartoon_highlight_color);

  cylindrical_helices = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_cartoon_cylindrical_helices);
  int const sampling_cylindrical_helices = sampling / 8 + 1;

  sampling_tmp = pymol::malloc<float>(sampling * 3);
  cartoon_debug = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_cartoon_debug);

  cgo = CGONew(G);
  if (alpha != 1.0F) CGOAlpha(cgo, alpha);
  /* debugging output */
  if (cartoon_debug > 0.5 && cartoon_debug < 2.5) {
    GenerateRepCartoonDrawDebugLineAlongPath(cgo, nAt, pv);
  }
  PRINTFD(G, FB_RepCartoon)
  " RepCartoon-Debug: creating 3D scaffold...\n" ENDFD;

  /* okay, we now have enough info to generate smooth interpolations */
  if (nAt > 1) {
    ex = ExtrudeNew(G);
    CHECKOK(ok, ex);
    if (ok) {
      // Spline mode shades twisted strand/arrow quads with true surface normals
      // (see ExtrudeTwistNormals); the classic path keeps its face normals.
      ex->twist_normals = spline_mode ? 1 : 0;
      ok &=
          ExtrudeAllocPointsNormalsColors(ex, cs->NIndex * (3 * sampling + 3));
    }
  }
  /* process cylindrical helices first */
  if (ok && (nAt > 1) &&
      cylindrical_helices == CARTOON_CYLINDRICAL_HELICES_STRAIGHT) {
    ok = GenerateRepCartoonProcessCylindricalHelices(
        G, obj, cs, cgo, ex, nAt, seg, pv, tv, pvo, car, at, dl, cartoon_color,
        discrete_colors, loop_radius, alpha);
  }
  if (ok && nAt > 1) {
    EXTRUDE_TRUNCATE();
    v1 = pv; /* points */
    v2 = tv; /* tangents */
    vo = pvo;
    d = dl;
    segptr = seg;
    cc = car;
    atp = at; /* cs index pointer */
    a = 0;
    contFlag = true;
    cur_car = cCartoon_skip;
    extrudeFlag = false;
    contigFlag = false;

    // Pre-compute which atoms at the start of a cyclic segment will be consumed
    // by the closing edge extension. The closing edge extrusion continues past
    // firstIdx so that the tube is one continuous piece with no seam.
    std::unordered_set<int> cyclicSkipAtoms;
    if (cyclicSegs) {
      for (auto& [segId, cinfo] : *cyclicSegs) {
        int closureCar = car[cinfo.lastIdx].getCCIn();
        // Check if the first atoms match the closure cartoon type
        for (int i = cinfo.firstIdx; i < cinfo.lastIdx; i++) {
          if (seg[i] != segId) break;
          if (car[i].getCCIn() != closureCar) break;
          cyclicSkipAtoms.insert(i);
        }
      }
    }

    const auto helix_radius =
        SettingGet<float>(G, cs->Setting.get(), obj->Setting.get(),
                          cSetting_cartoon_helix_radius);

    while (contFlag) {
      // Skip atoms that will be consumed by cyclic closure extension
      if (!cyclicSkipAtoms.empty() && cyclicSkipAtoms.count(a)) {
        v1 += 3;
        v2 += 3;
        vo += 3;
        d++;
        atp += 1;
        segptr++;
        cc++;
        a++;
        if (a == nAt) {
          contFlag = false;
        }
        continue;
      }

      if (CheckExtrudeContigFlags(nAt, n_p, a, &cur_car, cc, segptr,
                                  &contigFlag, &extrudeFlag)) {
        EXTRUDE_TRUNCATE();
      }

      if (ok && !extrudeFlag) {
        if ((a < (nAt - 1)) &&
            (*segptr == *(segptr + 1))) { /* working in the same segment... */
          AtomInfoType *ai1, *ai2;
          atom_index1 = cs->IdxToAtm[*atp];
          atom_index2 = cs->IdxToAtm[*(atp + 1)];
          ai1 = obj->AtomInfo + atom_index1;
          ai2 = obj->AtomInfo + atom_index2;

          float alpha1 = alpha;
          float alpha2 = alpha;

          ComputeCartoonAtomColors(G, obj, cs, nuc_flag, atom_index1,
                                   atom_index2, &c1, &c2, atp, cc, cur_car,
                                   cartoon_color, alpha1, alpha2, nucleic_color,
                                   discrete_colors, n_p, contigFlag);
          dev = throw_ * (*d);

          auto const cur_sampling = (cur_car == cCartoon_cylinder)
                                        ? sampling_cylindrical_helices
                                        : sampling;

          int const interval =
              (int)((v1 - pv) / 3);  // atom index of this interval's start
          if (!spline_coefs.valid.empty() && spline_coefs.valid[interval]) {
            CartoonGenerateSampleSpline(
                G, cur_sampling, &n_p, spline_coefs.coef.data() + 12 * interval,
                vo, c1, c2, alpha1, alpha2, ai1->masked ? -1 : atom_index1,
                ai2->masked ? -1 : atom_index2, &vc, &valpha, &vi, &v, &vn);
          } else {
            CartoonGenerateSample(
                G, cur_sampling, &n_p, dev, vo, v1, v2, c1, c2, alpha1, alpha2,
                ai1->masked ? -1 : atom_index1, ai2->masked ? -1 : atom_index2,
                power_a, power_b, &vc, &valpha, &vi, &v, &vn);

            /* now do a smoothing pass along orientation
               vector to smooth helices, etc... */
            CartoonGenerateRefine(refine, cur_sampling, v, vn, vo,
                                  sampling_tmp);
          }
        }
        v1 += 3;
        v2 += 3;
        vo += 3;
        d++;
        atp += 1;
        segptr++;
        cc++;
      }

      a++;
      if (a == nAt) {  // if at end, don't continue and extrude if needed
        contFlag = false;
        if (n_p) extrudeFlag = true;
      }
      if (ok && extrudeFlag) {
        contigFlag = true;
        if ((a < nAt) && extrudeFlag) {
          if (*(segptr - 1) != *(segptr)) contigFlag = false;
        }

        // Detect if this extrusion closes a cyclic segment.
        // Only trigger when:
        // 1. We've reached the last atom of the cyclic segment
        // 2. The first atom's cartoon type matches the current extrusion type
        //    (otherwise we'd stretch a tube/loop into an arrow/sheet)
        bool isCyclicExtrusion = false;
        if (ok && cyclicSegs && n_p > 0) {
          int lastSegId = *(segptr - 1);
          auto cit = cyclicSegs->find(lastSegId);
          if (cit != cyclicSegs->end() && (a - 1) == cit->second.lastIdx) {
            auto const& cinfo = cit->second;
            int firstCarType = car[cinfo.firstIdx].getCCIn();
            if (firstCarType == cur_car) {
              isCyclicExtrusion = true;

              // Generate closing edge samples (last atom -> first atom)
              float close_v1[6], close_v2[6], close_vo[6];
              copy3f(pv + cinfo.lastIdx * 3, close_v1);
              copy3f(pv + cinfo.firstIdx * 3, close_v1 + 3);
              copy3f(tv + cinfo.lastIdx * 3, close_v2);
              copy3f(tv + cinfo.firstIdx * 3, close_v2 + 3);
              copy3f(pvo + cinfo.lastIdx * 3, close_vo);
              copy3f(pvo + cinfo.firstIdx * 3, close_vo + 3);

              int closeAtom1 = cs->IdxToAtm[at[cinfo.lastIdx]];
              int closeAtom2 = cs->IdxToAtm[at[cinfo.firstIdx]];
              auto* ai1 = obj->AtomInfo + closeAtom1;
              auto* ai2 = obj->AtomInfo + closeAtom2;

              float alpha1 = alpha, alpha2 = alpha;
              ComputeCartoonAtomColors(
                  G, obj, cs, nuc_flag, closeAtom1, closeAtom2, &c1, &c2,
                  &at[cinfo.lastIdx], &car[cinfo.lastIdx], cur_car,
                  cartoon_color, alpha1, alpha2, nucleic_color, discrete_colors,
                  n_p, contigFlag);

              dev = throw_ * dl[cinfo.lastIdx];

              auto const cur_sampling = (cur_car == cCartoon_cylinder)
                                            ? sampling_cylindrical_helices
                                            : sampling;

              CartoonGenerateSample(G, cur_sampling, &n_p, dev, close_vo,
                                    close_v1, close_v2, c1, c2, alpha1, alpha2,
                                    ai1->masked ? -1 : closeAtom1,
                                    ai2->masked ? -1 : closeAtom2, power_a,
                                    power_b, &vc, &valpha, &vi, &v, &vn);

              CartoonGenerateRefine(refine, cur_sampling, v, vn, close_vo,
                                    sampling_tmp);

              // Continue past firstIdx through skipped atoms so this extrusion
              // covers the full loop section with no seam at the joint
              for (int ci = cinfo.firstIdx;
                   ci < cinfo.lastIdx && seg[ci] == lastSegId &&
                   car[ci].getCCIn() == cur_car && seg[ci] == seg[ci + 1];
                   ci++) {
                float cont_v1[6], cont_v2[6], cont_vo[6];
                copy3f(pv + ci * 3, cont_v1);
                copy3f(pv + (ci + 1) * 3, cont_v1 + 3);
                copy3f(tv + ci * 3, cont_v2);
                copy3f(tv + (ci + 1) * 3, cont_v2 + 3);
                copy3f(pvo + ci * 3, cont_vo);
                copy3f(pvo + (ci + 1) * 3, cont_vo + 3);

                int contAtom1 = cs->IdxToAtm[at[ci]];
                int contAtom2 = cs->IdxToAtm[at[ci + 1]];
                auto* cai1 = obj->AtomInfo + contAtom1;
                auto* cai2 = obj->AtomInfo + contAtom2;

                float ca1 = alpha, ca2 = alpha;
                ComputeCartoonAtomColors(
                    G, obj, cs, nuc_flag, contAtom1, contAtom2, &c1, &c2,
                    &at[ci], &car[ci], cur_car, cartoon_color, ca1, ca2,
                    nucleic_color, discrete_colors, n_p, contigFlag);

                dev = throw_ * dl[ci];
                CartoonGenerateSample(G, cur_sampling, &n_p, dev, cont_vo,
                                      cont_v1, cont_v2, c1, c2, ca1, ca2,
                                      cai1->masked ? -1 : contAtom1,
                                      cai2->masked ? -1 : contAtom2, power_a,
                                      power_b, &vc, &valpha, &vi, &v, &vn);
                CartoonGenerateRefine(refine, cur_sampling, v, vn, cont_vo,
                                      sampling_tmp);
              }
            }  // if firstCarType matches
          }
        }

        if (ok && (cur_car != cCartoon_skip) &&
            (cur_car != cCartoon_skip_helix)) {
          if ((cartoon_debug > 0.5) && (cartoon_debug < 2.5)) {
            ok = GenerateRepCartoonDrawDebugNormals(cgo, ex, n_p);
          }

          if (ok) {
            ExtrudeTruncate(ex, n_p);
            ok &= ExtrudeComputeTangents(ex);
          }
          if (ok) {
            /* set up shape */
            switch (cur_car) {
              case cCartoon_tube:
                ok = CartoonExtrudeTube(
                    use_cylinders_for_strands, ex, cgo, tube_radius,
                    tube_quality, isCyclicExtrusion ? cCylCap::Flat : tube_cap);
                break;
              case cCartoon_putty:
                ok = CartoonExtrudePutty(G, obj, cs, cgo, ex, putty_quality,
                                         putty_radius, putty_vals, sampling);
                if (!ok) contFlag = false;
                break;
              case cCartoon_loop:
                ok = CartoonExtrudeCircle(
                    ex, cgo, use_cylinders_for_strands, loop_quality,
                    loop_radius, isCyclicExtrusion ? cCylCap::None : loop_cap);
                break;
              case cCartoon_dash:
                ok = CartoonExtrudeCircle(ex, cgo, use_cylinders_for_strands,
                                          loop_quality, loop_radius, loop_cap,
                                          2);
                break;
              case cCartoon_rect:
                ok = CartoonExtrudeRect(G, ex, cgo, width, length,
                                        highlight_color);
                break;
              case cCartoon_oval:
                ok = CartoonExtrudeOval(G, ex, cgo, use_cylinders_for_strands,
                                        oval_quality, oval_width, oval_length,
                                        highlight_color);
                break;
              case cCartoon_arrow:
                ok = CartoonExtrudeArrow(G, ex, cgo, sampling, width, length,
                                         highlight_color);
                break;
              case cCartoon_dumbbell:
                ok = CartoonExtrudeDumbbell(
                    G, ex, cgo, sampling, dumbbell_width, dumbbell_length,
                    highlight_color, loop_quality, dumbbell_radius,
                    use_cylinders_for_strands);
                break;
              case cCartoon_cylinder:
                CartoonExtrudeCurvedCylindricalHelix(
                    ex, cgo, loop_quality * 2, helix_radius,
                    sampling_cylindrical_helices);
                break;
            }
            if (!ok) contFlag = false;
          }
        }
        a--; /* undo above... */
        extrudeFlag = false;
        if (ok) {
          EXTRUDE_TRUNCATE();  // doesn't include vi = ex->i, not used?
        }
        /*
          ExtrudeTruncate(ex, 0);
        n_p = 0;
        v = ex->p;
        vc = ex->c;
        vn = ex->n;*/
      }
    }
  }

  if (ok && nAt > 1) {
    if ((cartoon_debug > 0.5) && (cartoon_debug < 2.5)) {
      ok = GenerateRepCartoonDrawDebugOrient(cgo, nAt, pv, pvo, tv);
    }
  }
  if (ex) {
    ExtrudeFree(ex);
  }
  /* draw the rings */
  if (ok && ndata->ring_anchor && ndata->n_ring) {
    ok = GenerateRepCartoonDrawRings(G, ndata, obj, cs, cgo, ring_width,
                                     cartoon_color, alpha);
  }
  if (ok) CGOStop(cgo);

  FreeP(sampling_tmp);

  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

bool RepCartoon::sameVis() const {
  if (!LastVisib) return false;

  for (int idx = 0; idx < cs->NIndex; idx++) {
    const auto* ai = cs->getAtomInfo(idx);
    if (LastVisib[idx] != GET_BIT(ai->visRep, cRepCartoon)) {
      return false;
    }
  }

  return true;
}

void RepCartoon::invalidate(cRepInv_t level) {
  if (level >= cRepInvColor) {
    FreeP(LastVisib);
  }
  Rep::invalidate(level);
}

/**
 * compute and fill "ndata" with:
 * - cartoon trace and segments
 * - cartoon types and secondary structure
 * - rings
 */
static void RepCartoonGeneratePASS1(PyMOLGlobals* G, RepCartoon* I,
                                    ObjectMolecule* obj, CoordSet* cs,
                                    nuc_acid_data* ndata) {
  int st, nd;
  int a, a1, a3, a4 = 0;
  char* lv = I->LastVisib;
  int trace, trace_mode;
  AtomInfoType *ai, *last_ai = nullptr;
  int cartoon_side_chain_helper;
  int cylindrical_helices;
  int fancy_helices;
  int fancy_sheets;
  int parity = 1;
  float *v_c, *v_n, *v_o;
  int cur_car;
  nuc_acid_cap leading_O5p(G, ndata, cs, 3);
  nuc_acid_cap trailing_O3p(G, ndata, cs, 2);

  // settings
  fancy_sheets = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_cartoon_fancy_sheets);
  fancy_helices = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_cartoon_fancy_helices);
  cylindrical_helices = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_cartoon_cylindrical_helices);
  cartoon_side_chain_helper =
      SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                   cSetting_cartoon_side_chain_helper);
  int trace_ostate = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_cartoon_trace_atoms);
  trace_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_trace_atoms_mode);
  auto gap_cutoff = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_cartoon_gap_cutoff);

  /**
   * JJ: I think this for loop is a bit ill conceived. Nuc_acid_data
   * skips over consideration of invisible cartoons, which has historically
   * worked since it's compensated for later, but this approach doesn't consider
   * the difference between removed vs hidden residues which can have different
   * representations. For now, use this flag to track hidden (but valid atoms)
   * cartoons, which should be skipped and not dashed. Removed residues should
   * be dashed. See PYMOL-3168
   */
  bool invisFound = false;
  // iterate over (sorted) atoms
  for (CoordSetAtomIterator iter(cs); iter.next();) {
    ai = iter.getAtomInfo();

    // cartoon rep for this atom?
    if (!(*(lv++) = GET_BIT(ai->visRep, cRepCartoon))) {
      invisFound = true;
      continue;
    }

    const char* ai_name = LexStr(G, ai->name);

    // atom indices
    a = iter.getIdx();
    a1 = iter.getAtm();

    // rings
    if (ndata->ring_anchor && ai->protons != cAN_H &&
        (ndata->ring_finder_eff >= 3     /* all 5-7 atom rings */
         || (ndata->ring_finder_eff <= 2 /*  C4-containing rings */
             && (WordMatchExact(G, "C4", ai_name, 1))) ||
         (ndata->ring_finder_eff == 1 &&
          (WordMatchExact(G, "C4*", ai_name, 1) ||
           WordMatchExact(G, "C4'", ai_name, 1))))) {
      VLACheck(ndata->ring_anchor, int, ndata->n_ring);
      ndata->ring_anchor[ndata->n_ring] = a1;
      ndata->n_ring++;
    }

    // handle alternative conformations
    if (ai->alt[0]) {
      if (!ndata->alt) {
        ndata->alt = ai->alt[0];
      } else if (ai->alt[0] != ndata->alt) {
        if (ai->alt[0] > ndata->alt &&
            (!ndata->next_alt || ai->alt[0] < ndata->next_alt))
          ndata->next_alt = ai->alt[0];
        continue;
      }
    }

    // atom level setting
    trace = AtomSettingGetWD(G, ai, cSetting_cartoon_trace_atoms, trace_ostate);

    // CA or cartoon_trace_atoms
    if (trace ||
        (ai->protons == cAN_C && WordMatchExact(G, "CA", ai_name, true) &&
         !AtomInfoSameResidueP(G, last_ai, ai))) {
      PRINTFD(G, FB_RepCartoon)
      " RepCartoon: found CA in %d%c; a2 %d\n", ai->resv, ai->getInscode(true),
          ndata->a2 ENDFD;

      // 3' nucleic acid cap
      if (trailing_O3p.cap()) ndata->a2 = -1;

      bool is_gap = false;

      // auto-detect CA-only models
      if (!ai->bonded) trace = true;

      // check for gap
      if (ndata->a2 >= 0) {
        if (!trace) {
          // CA->N->C->CA = 3 bonds
          if (!ObjectMoleculeCheckBondSep(obj, a1, ndata->a2, 3)) is_gap = true;
        } else {
          if (!AtomInfoSequential(G, obj->AtomInfo + ndata->a2, ai, trace_mode))
            is_gap = true;
        }
      }

      cur_car = ai->cartoon;

      if (is_gap) {
        int delta = ai->resv - last_ai->resv;
        if (delta < 1 || delta > gap_cutoff ||
            !AtomInfoSameChainP(G, ai, last_ai)) {
          ndata->a2 = -1;
        } else {
          auto gapRep = invisFound ? cCartoon_skip : cCartoon_dash;
          (ndata->cc - 1)->setCCOut(gapRep);
          invisFound = false;
        }
      }

      last_ai = ai;

      // start new segment?
      if (ndata->a2 < 0) ndata->nSeg++;

      (*ndata->fp) = ai->flags; /* store atom flags */

      // side_chain_helper
      if ((ai->visRep & (cRepLineBit | cRepCylBit | cRepSphereBit)) &&
          AtomSettingGetWD(G, ai, cSetting_cartoon_side_chain_helper,
                           cartoon_side_chain_helper)) {
        (*ndata->fp) |= cAtomFlag_no_smooth;
      }

      // secondary structure type
      switch (ai->ssType[0]) {
        case 'H':
        case 'h':
          if (cur_car == cCartoon_auto) {
            if (cylindrical_helices == CARTOON_CYLINDRICAL_HELICES_STRAIGHT) {
              cur_car = cCartoon_skip_helix;
            } else if (cylindrical_helices ==
                       CARTOON_CYLINDRICAL_HELICES_CURVED) {
              cur_car = cCartoon_cylinder;
            } else if (fancy_helices) {
              cur_car = cCartoon_dumbbell;
            } else {
              cur_car = cCartoon_oval;
            }
          }
          (*ndata->ss) = ss_t::HELIX;
          parity = 0;
          break;
        case 'S':
        case 's':
          if (cur_car == cCartoon_auto) {
            cur_car = fancy_sheets ? cCartoon_arrow : cCartoon_rect;
          }
          (*ndata->ss) = ss_t::SHEET;
          parity = !parity;
          break;
        default: /* 'L', 'T', 0, etc. */
          if (cur_car == cCartoon_auto) {
            cur_car = cCartoon_loop;
          }
          parity = 0;
          (*ndata->ss) = ss_t::NONE;
          break;
      }

      if (cur_car == cCartoon_putty) ndata->putty_flag = true;

      // coordinates
      copy3f(cs->coordPtr(a), ndata->vptr);

      *((ndata->cc)++) = cur_car;
      ndata->a2 = a1;
      ndata->vptr += 3;
      ndata->ss++;
      ndata->fp++;
      *(ndata->sptr++) = ndata->nSeg;
      ndata->nAt++;
      *(ndata->iptr++) = a;

      if (trace) {
        if (a1 == 0 || a1 + 1 == obj->NAtom ||
            (a3 = cs->atmToIdx(a1 - 1)) == -1 ||
            (a4 = cs->atmToIdx(a1 + 1)) == -1) {
          zero3f(ndata->voptr);
        } else {
          float t0[3], t1[3];
          subtract3f(cs->coordPtr(a), cs->coordPtr(a3), t0);
          subtract3f(cs->coordPtr(a), cs->coordPtr(a4), t1);
          add3f(t0, t1, ndata->voptr);
          normalize3f(ndata->voptr);
        }
        ndata->voptr += 3;
        continue;
      }

      // pointers to C+N+O coordinates
      v_c = nullptr;
      v_n = nullptr;
      v_o = nullptr;

      // get start (st) and end (nd) indices of residue atoms
      AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, a1, &st, &nd);

      // find C+N+O coordinates
      for (a3 = st; a3 <= nd; a3++) {
        a4 = cs->atmToIdx(a3);

        if (a4 == -1) continue;

        const char* a3name = LexStr(G, obj->AtomInfo[a3].name);

        if (WordMatchExact(G, "C", a3name, true)) {
          v_c = cs->coordPtr(a4);
        } else if (WordMatchExact(G, "N", a3name, true)) {
          v_n = cs->coordPtr(a4);
        } else if (WordMatchExact(G, "O", a3name, true)) {
          v_o = cs->coordPtr(a4);
        }
      }

      // orientation vector
      if (!(v_c && v_n && v_o)) {
        zero3f(ndata->voptr);
      } else {
        float t0[3], t1[3];
        subtract3f(v_n, v_c, t0); /* t0 = N<---C */
        normalize3f(t0);
        subtract3f(v_n, v_o, t1); /* t1 = N<---O */
        normalize3f(t1);
        cross_product3f(t0, t1, ndata->voptr);
        normalize3f(ndata->voptr);
        if (parity) {
          invert3f(ndata->voptr);
        }
      }
      ndata->voptr += 3;

    } else if (!AtomInfoSameResidueP(G, last_ai, ai) &&
               (ndata->na_mode != 1
                    ? AtomInfoIsNucBackboneTrace(ai_name, ai->protons)
                    :
                    // C3* C3' atom
                    (ai->protons == cAN_C &&
                     (WordMatchExact(G, NUCLEIC_NORMAL1, ai_name, 1) ||
                      WordMatchExact(G, NUCLEIC_NORMAL2, ai_name, 1))))) {
      // check for gap
      if (ndata->a2 >= 0 &&
          // six bonds between phosphates
          !ObjectMoleculeCheckBondSep(obj, a1, ndata->a2, 6)) {
        /*  3' cap */
        trailing_O3p.cap();
        ndata->a2 = -1;
      }

      last_ai = ai;
      trailing_O3p.set(-1, -1, nullptr);

      /*  5' cap */
      if (leading_O5p.active() && ndata->a2 < 0 &&
          !AtomInfoSameResidueP(G, ai, leading_O5p.ai) &&
          ObjectMoleculeCheckBondSep(obj, a1, leading_O5p.atm, 5)) {
        leading_O5p.cap();
      }

      leading_O5p.set(-1, -1, nullptr);

      /* this is the main nucleic acid cartoon section... */
      nuc_acid(G, ndata, a, a1, ai, cs, obj, true);

    } else if (
        // P -> O3* bond
        AtomInfoSameResidueP(G, last_ai, ai) && ndata->a2 >= 0 && last_ai &&
        cAN_P == last_ai->protons && cAN_O == ai->protons &&
        trailing_O3p.enabled &&
        (WordMatchExact(G, "O3'", ai_name, 1) ||
         WordMatchExact(G, "O3*", ai_name, 1)) &&
        ObjectMoleculeCheckBondSep(obj, a1, ndata->a2, 5)) {
      // remember trailing O3*
      trailing_O3p.set(a, a1, ai);

    } else if (
        // O5* atom
        ai->protons == cAN_O && leading_O5p.enabled &&
        (WordMatchExact(G, "O5'", ai_name, 1) ||
         WordMatchExact(G, "O5*", ai_name, 1))) {
      leading_O5p.set(a, a1, ai);
    }
  }

  /* BEGIN 3' cap */
  if (trailing_O3p.cap()) {
    ndata->a2 = -1;
  }
}

Rep* RepCartoonNew(CoordSet* cs, int state) {
  PyMOLGlobals* G = cs->G;
  ObjectMolecule* obj;
  int *i, *sptr, *at, *seg, nAt;
  CCInOut *car, *cc;
  float* pv = nullptr;
  float *pvo = nullptr, *pva = nullptr;
  float* dv = nullptr;
  float* nv = nullptr;
  float* tv = nullptr;
  float* tmp = nullptr;
  float* dl = nullptr;

  int ladder_mode;
  int round_helices;
  int na_mode;
  int* flag_tmp;
  float putty_vals[4] = {
      10.0F, 0.0F, FLT_MAX,
      -FLT_MAX};  // putty_mean, putty_stdev, putty_min, putty_max
  int* ring_anchor = nullptr;
  int* nuc_flag = nullptr;
  float alpha;
  int ok = true;
  nuc_acid_data ndata;
  short use_shaders = SettingGetGlobal_b(G, cSetting_use_shaders);
  short na_strands_as_cylinders =
      use_shaders &&
      (SettingGetGlobal_i(G, cSetting_cartoon_nucleic_acid_as_cylinders) & 2) &&
      SettingGetGlobal_b(G, cSetting_render_as_cylinders);

  // skip if not visible
  if (!cs->hasRep(cRepCartoonBit)) return nullptr;

  /* THIS IS BY FAR THE WORST ROUTINE IN PYMOL!
   * DEVELOP ON IT ONLY AT EXTREME RISK TO YOUR MENTAL HEALTH */

  auto I = new RepCartoon(cs, state);

  PRINTFD(G, FB_RepCartoon)
  " RepCartoonNew-Debug: entered.\n" ENDFD;

  obj = cs->Obj;

  alpha = 1.0F - SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_cartoon_transparency);
  round_helices = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_cartoon_round_helices);
  na_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_cartoon_nucleic_acid_mode);
  ladder_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_ladder_mode);

  /* find all of the CA points */

  auto const nAtIndex = obj->NAtom;  // was NAtIndex
  at = pymol::malloc<int>(nAtIndex); /* cs index pointers */
  pv = pymol::malloc<float>(nAtIndex * 3);
  tmp = pymol::malloc<float>(nAtIndex * 3);
  pvo = pymol::malloc<float>(nAtIndex * 3); /* orientation vector */
  pva = pymol::malloc<float>(
      nAtIndex * 6); /* alternative orientation vectors, two per atom */
  seg = pymol::malloc<int>(nAtIndex);
  car = pymol::calloc<CCInOut>(nAtIndex); /* cartoon type for each atom */
  auto sstype = pymol::malloc<ss_t>(nAtIndex);
  flag_tmp = pymol::calloc<int>(nAtIndex);
  nuc_flag = pymol::calloc<int>(nAtIndex);

  I->LastVisib = pymol::calloc<char>(nAtIndex);

  auto cartoon_all_alt = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                      cSetting_cartoon_all_alt);

  ndata.next_alt = 0;

  do {
    ndata.alt = ndata.next_alt;
    ndata.next_alt = 0;
    memset(car, 0, sizeof(*car) * nAtIndex);
    memset(sstype, 0, sizeof(*sstype) * nAtIndex);
    memset(flag_tmp, 0, sizeof(*flag_tmp) * nAtIndex);
    memset(nuc_flag, 0, sizeof(*nuc_flag) * nAtIndex);

    i = at;
    sptr = seg;
    cc = car;
    auto* ss = sstype;
    nAt = 0;

    ndata.na_mode = na_mode;
    ndata.nuc_flag = nuc_flag;
    ndata.a2 = -1;
    ndata.nSeg = 0;
    ndata.v_o_last = nullptr;
    ndata.sptr = sptr;
    ndata.iptr = i;
    ndata.cc = cc;
    ndata.nAt = nAt;
    ndata.ss = ss;
    ndata.putty_flag = false;
    ndata.fp = flag_tmp;
    ndata.vptr = pv;
    ndata.voptr = pvo;
    ndata.ring_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_cartoon_ring_mode);
    ndata.ring_finder = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_cartoon_ring_finder);
    ndata.ring_finder_eff = ndata.ring_finder;
    if ((!ndata.ring_mode) || (ndata.ring_finder == 2))
      ndata.ring_finder_eff = 1;
    if (ndata.ring_mode || ladder_mode) {
      ring_anchor = VLAlloc(int, nAtIndex / 10 + 1);
    }
    ndata.ring_anchor = ring_anchor;
    ndata.n_ring = 0;

    RepCartoonGeneratePASS1(G, I, obj, cs, &ndata);
    nAt = ndata.nAt;
    if (nAt && ndata.putty_flag) {
      RepCartoonComputePuttyValues(obj, putty_vals);
    }

    // Detect cyclic segments (e.g. cyclic peptides like 2NB5)
    auto cyclicSegs = RepCartoonDetectCyclicSegments(obj, cs, nAt, seg, at);

    PRINTFD(G, FB_RepCartoon)
    " RepCartoon-Debug: path outlined, interpolating... nAt=%d\n", nAt ENDFD;

    if (nAt) {
      dv = pymol::malloc<float>(
          nAt * 3); /* differences between next and current 3f */
      nv = pymol::malloc<float>(nAt * 3); /* normal */
      dl = pymol::malloc<float>(
          nAt); /* length (i.e., normal * length = difference) */
      if (cyclicSegs.empty()) {
        RepCartoonComputeDifferencesAndNormals(G, nAt, seg, pv, dv, nv, dl,
                                               true);
      } else {
        RepCartoonComputeDifferencesAndNormalsCyclic(G, nAt, seg, pv, dv, nv,
                                                     dl, true, cyclicSegs);
      }

      /* compute tangents */
      tv = pymol::malloc<float>(nAt * 3 + 6);
      if (cyclicSegs.empty()) {
        RepCartoonComputeTangents(nAt, seg, nv, tv);
      } else {
        RepCartoonComputeTangentsCyclic(nAt, seg, nv, tv, cyclicSegs);
      }

      PRINTFD(G, FB_RepCartoon)
      " RepCartoon-Debug: generating coordinate systems...\n" ENDFD;

      if (round_helices) {
        ndata.voptr = pvo;
        RepCartoonComputeRoundHelices(&ndata, nAt, seg, sstype, tv, pv);
      }

      RepCartoonRefineNormals(G, I, obj, cs, &ndata, nAt, seg, tv, pvo, pva,
                              sstype, nv);

      {
        int smooth_loops =
            SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_cartoon_smooth_loops);
        bool cartoon_flat_sheets =
            SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_cartoon_flat_sheets);
        if (smooth_loops || cartoon_flat_sheets) {
          if (cartoon_flat_sheets) {
            RepCartoonFlattenSheets(G, obj, cs, &ndata, nAt, seg, car, pv, pvo,
                                    sstype, tv, tmp, flag_tmp);
          }
          if (smooth_loops) {
            RepCartoonSmoothLoops(G, obj, cs, &ndata, nAt, seg, pv, sstype, pvo,
                                  tv, tmp, flag_tmp);
          }
          /* recompute differences and normals */
          if (cyclicSegs.empty()) {
            RepCartoonComputeDifferencesAndNormals(G, nAt, seg, pv, dv, nv, dl,
                                                   true);
          } else {
            RepCartoonComputeDifferencesAndNormalsCyclic(
                G, nAt, seg, pv, dv, nv, dl, true, cyclicSegs);
          }
          /* recompute tangents */
          if (cyclicSegs.empty()) {
            RepCartoonComputeTangents(nAt, seg, nv, tv);
          } else {
            RepCartoonComputeTangentsCyclic(nAt, seg, nv, tv, cyclicSegs);
          }
          if (cartoon_flat_sheets) {
            RepCartoonFlattenSheetsRefineTips(G, obj, cs, nAt, seg, sstype, tv);
          }
        }
      }
    }

    auto* cyclicSegsPtr = cyclicSegs.empty() ? nullptr : &cyclicSegs;
    CGO* preshadercgo = GenerateRepCartoonCGO(
        cs, obj, &ndata, na_strands_as_cylinders, pv, nAt, tv, pvo, dl, car,
        seg, at, nuc_flag, putty_vals, alpha, cyclicSegsPtr);

    if (preshadercgo && preshadercgo->has_begin_end) {
      CGOCombineBeginEnd(&preshadercgo);
    }

    if (I->preshader) {
      I->preshader->free_append(preshadercgo);
    } else {
      I->preshader = preshadercgo;
    }
  } while (ndata.next_alt && cartoon_all_alt);

  CHECKOK(ok, I->preshader);

  ok &= !G->Interrupt;
  if (!ok || !CGOHasOperations(I->preshader)) {
    /* cannot generate RepCartoon */
    delete I;
    I = nullptr;
  }
  FreeP(dv);
  FreeP(dl);
  FreeP(tv);
  FreeP(nv);
  FreeP(at);
  FreeP(seg);
  FreeP(pv);
  FreeP(pvo);
  FreeP(pva);
  FreeP(car);
  FreeP(tmp);
  FreeP(sstype);
  FreeP(flag_tmp);
  FreeP(nuc_flag);
  VLAFreeP(ndata.ring_anchor);
  return (Rep*)I;
}
