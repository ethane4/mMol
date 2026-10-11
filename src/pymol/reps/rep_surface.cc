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
#include "pymol/base/busy.h"

#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/python/p.h"
#include "pymol/reps/rep_surface.h"
#include "pymol/scene/scene.h"
#include "pymol/select/selector.h"

#include "pymol/reps/rep_surface_internal.h"

static void RepSurfaceSmoothEdges(RepSurface* I);

void setShaderCGO(RepSurface* I, CGO* cgo) {
  if (I->shaderCGO != I->pickingCGO) {
    CGOFree(I->shaderCGO);
  }
  I->shaderCGO = cgo;
}

static void setPickingCGO(RepSurface* I, CGO* cgo) {
  if (I->shaderCGO != I->pickingCGO) {
    CGOFree(I->pickingCGO);
  }
  I->pickingCGO = cgo;
}

RepSurface::~RepSurface() {
  setPickingCGO(this, nullptr);
  setShaderCGO(this, nullptr);
}

void SolventDotFree(SolventDot* I) {
  if (I) {
    VLAFreeP(I->dot);
    VLAFreeP(I->dotNormal);
    VLAFreeP(I->dotCode);
  }
  DeleteP(I);
}

bool visibility_test(bool proximityFlag,
                     const int* vi,  // visibility array
                     const int* t)   // indices
{
  if (proximityFlag) return (vi[t[0]] || vi[t[1]] || vi[t[2]]);
  return (vi[t[0]] && vi[t[1]] && vi[t[2]]);
}

static int check_and_add(int* cache, int spacing, int t0, int t1) {
  int* rec;
  int cnt;
  t0++;
  t1++;

  rec = cache + spacing * t0;
  cnt = spacing;
  while (cnt > 0) {
    if (*rec == t1) return 1;
    if (!*rec) {
      *rec = t1;
      break;
    }
    rec++;
    cnt--;
  }
  rec = cache + spacing * t1;
  cnt = spacing;
  while (cnt > 0) {
    if (*rec == t0) return 1;
    if (!*rec) {
      *rec = t0;
      break;
    }
    rec++;
    cnt--;
  }
  return 0;
}

#define CLAMP_VALUE(v) ((v > 1.f) ? 1.f : (v < 0.f) ? 0.f : v)

static int AtomInfoIsMasked(ObjectMolecule* obj, int atm) {
  AtomInfoType* ait;
  if (atm < 0) return cPickableNoPick;
  ait = &obj->AtomInfo[atm];
  return (ait->masked ? cPickableNoPick : cPickableAtom);
}

static int RepSurfaceCGOGenerate(RepSurface* I, RenderInfo* info) {
  PyMOLGlobals* G = I->G;
  float* v = I->V.data();
  float* vn = I->VN.data();
  float* vc = I->VC.data();
  float* va = I->VA.data();
  int* t = I->T.data();
  int* s = I->S.data();
  int c = I->N;
  int* vi = I->Vis.data();
  int* at = I->AT.data();
  int ok = true;
  float alpha;
  int t_mode;
  CGO* convertcgo = nullptr;

  auto* const cs = I->cs;
  auto* const obj = I->cs->Obj;

  bool pick_surface = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_pick_surface);
  short dot_as_spheres = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                      cSetting_dot_as_spheres);
  alpha = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                       cSetting_transparency);
  alpha = 1.0F - alpha;
  if (fabs(alpha - 1.0) < R_SMALL4) alpha = 1.0F;

  setPickingCGO(I, nullptr);
  setShaderCGO(I, CGONew(G));

  if (!I->shaderCGO) return false;

  I->shaderCGO->use_shader = true;
  I->dot_as_spheres = dot_as_spheres;

  if (alpha < 1) {
    I->setHasTransparency();
  }

  if (I->Type == SurfaceType::DotDefault) {
    /* no triangle information, so we're rendering dots only */
    int normals = SettingGet<int>(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_dot_normals);
    if (!normals) {
      CGOResetNormal(I->shaderCGO, true);
    }
    if ((alpha != 1.0)) {
      CGOAlpha(I->shaderCGO, alpha);
    }
    if (dot_as_spheres) {
      if (c) {
        ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
        if (ok && I->oneColorFlag) {
          ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
        }
        while (ok && c--) {
          if (*vi) {
            if (!I->oneColorFlag) {
              ok &= CGOColorv(I->shaderCGO, vc);
            }
            if (ok && normals) ok &= CGONormalv(I->shaderCGO, vn);
            if (ok && pick_surface)
              ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked(obj, *at));

            if (ok) ok &= CGOSphere(I->shaderCGO, v, 1.f);
          }
          vi++;
          vc += 3;
          vn += 3;
          v += 3;
          at++;
        }
      }
    } else {
      ok &= CGODotwidth(I->shaderCGO,
                        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_dot_width));
      if (ok && c) {
        ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
        ok &= CGOBegin(I->shaderCGO, GL_POINTS);
        if (ok && I->oneColorFlag) {
          ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
        }

        while (ok && c--) {
          if (*vi) {
            if (!I->oneColorFlag) {
              ok &= CGOColorv(I->shaderCGO, vc);
            }
            if (normals) {
              CGONormalv(I->shaderCGO, vn);
            }
            if (ok && pick_surface)
              ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked(obj, *at));
            if (ok) ok &= CGOVertexv(I->shaderCGO, v);
          }
          vi++;
          vc += 3;
          vn += 3;
          v += 3;
          at++;
        }
        if (ok) ok &= CGOEnd(I->shaderCGO);
      }
    }
  } else if (I->Type == SurfaceType::Triangle) {
    int normals = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_mesh_normals);
    if (ok && !normals) {
      ok &= CGOResetNormal(I->shaderCGO, true);
    }
    if (ok) {
      float line_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                      cSetting_mesh_width);
      line_width = SceneGetDynamicLineWidth(info, line_width);

      ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_MESH);

      c = I->NT;
      if (ok && c) {
        if (I->oneColorFlag) {
          ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
          while (ok && c--) {
            if (visibility_test(I->proximity, vi, t)) {
              if (normals) {
                int idx;
                ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
                idx = (*(t + 2));
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

                idx = (*t);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                if (ok) ok &= CGOEnd(I->shaderCGO);
              } else {
                int idx;
                ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);

                idx = (*(t + 2));
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                idx = *t;
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = *t;
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = *t;
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                if (ok) ok &= CGOEnd(I->shaderCGO);
              }
            } else
              t += 3;
          }
        } else { /* not oneColorFlag */
          while (ok && c--) {
            if (visibility_test(I->proximity, vi, t)) {
              if (normals) {
                int idx;
                ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);

                idx = (*(t + 2));
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                if (ok) ok &= CGOEnd(I->shaderCGO);
              } else {
                int idx;
                ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);

                idx = (*(t + 2));
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                idx = (*t);
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                     AtomInfoIsMasked(obj, I->AT[idx]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
                t++;
                if (ok) ok &= CGOEnd(I->shaderCGO);
              }
            } else
              t += 3;
          }
        }
      }
    }
  } else {
    /* we're rendering triangles */
    if ((alpha != 1.0) || va) {
      t_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_transparency_mode);

      if (info && info->alpha_cgo) {
        t_mode = 0;
      }

      if (t_mode) {
        float **t_buf = nullptr, **tb;
        float *z_value = nullptr, *zv;
        int* ix = nullptr;
        int n_tri = 0;
        float sum[3];
        float matrix[16];

        glGetFloatv(GL_MODELVIEW_MATRIX, matrix);

        if (I->oneColorFlag) {
          t_buf = pymol::malloc<float*>(I->NT * 6);
        } else {
          t_buf = pymol::malloc<float*>(I->NT * 12);
        }
        CHECKOK(ok, t_buf);
        if (ok) {
          z_value = pymol::malloc<float>(I->NT);
          CHECKOK(ok, z_value);
        }
        if (ok) {
          ix = pymol::malloc<int>(I->NT);
          CHECKOK(ok, ix);
        }
        zv = z_value;
        tb = t_buf;
        c = I->NT;
        if (ok) {
          if (I->oneColorFlag) {
            while (c--) {
              if (visibility_test(I->proximity, vi, t)) {
                *(tb++) = vn + (*t) * 3;
                *(tb++) = v + (*t) * 3;
                *(tb++) = vn + (*(t + 1)) * 3;
                *(tb++) = v + (*(t + 1)) * 3;
                *(tb++) = vn + (*(t + 2)) * 3;
                *(tb++) = v + (*(t + 2)) * 3;

                add3f(tb[-1], tb[-3], sum);
                add3f(sum, tb[-5], sum);

                *(zv++) = matrix[2] * sum[0] + matrix[6] * sum[1] +
                          matrix[10] * sum[2];
                n_tri++;
              }
              t += 3;
            }
          } else {
            while (c--) {
              if (visibility_test(I->proximity, vi, t)) {
                if (va)
                  *(tb++) = va + (*t);
                else
                  *(tb++) = &alpha;
                *(tb++) = vc + (*t) * 3;
                *(tb++) = vn + (*t) * 3;
                *(tb++) = v + (*t) * 3;

                if (va)
                  *(tb++) = va + (*(t + 1));
                else
                  *(tb++) = &alpha;
                *(tb++) = vc + (*(t + 1)) * 3;
                *(tb++) = vn + (*(t + 1)) * 3;
                *(tb++) = v + (*(t + 1)) * 3;

                if (va)
                  *(tb++) = va + (*(t + 2));
                else
                  *(tb++) = &alpha;
                *(tb++) = vc + (*(t + 2)) * 3;
                *(tb++) = vn + (*(t + 2)) * 3;
                *(tb++) = v + (*(t + 2)) * 3;

                add3f(tb[-1], tb[-5], sum);
                add3f(sum, tb[-9], sum);
                *(zv++) = matrix[2] * sum[0] + matrix[6] * sum[1] +
                          matrix[10] * sum[2];
                n_tri++;
              }
              t += 3;
            }
          }
        }
        switch (t_mode) {
          case 1:
            ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, true);
            /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
            break;
          default:
            ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, false);
            /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
            break;
        }
        c = n_tri;
        if (I->oneColorFlag) {
          float col[3];
          ColorGetEncoded(G, I->oneColor, col);
          if (ok) {
            CGOAlpha(I->shaderCGO, alpha);
            if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
            if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
            for (c = 0; ok && c < n_tri; c++) {
              int idx;
              tb = t_buf + 6 * c;
              //		    tb = t_buf + 6 * ix[c];
              if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
              idx = ((*tb - v) / 3);
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                   AtomInfoIsMasked(obj, I->AT[idx]));
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
              }
              if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
              if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
              idx = ((*tb - v) / 3);
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                   AtomInfoIsMasked(obj, I->AT[idx]));
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
              }
              if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
              if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
              idx = ((*tb - v) / 3);
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                   AtomInfoIsMasked(obj, I->AT[idx]));
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
              }
              if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
            }
            if (ok) ok &= CGOEnd(I->shaderCGO);
          }
        } else { /* else I->oneColorFlag */
          if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
          for (c = 0; ok && c < n_tri; c++) {
            float *vv, *v_alpha;
            int idx;
            tb = t_buf + 12 * c; /* need to update index every frame */
            //		    tb = t_buf + 12 * ix[c];
            v_alpha = *(tb++);
            vv = *(tb++);
            ok &= CGOAlpha(I->shaderCGO, *v_alpha);
            if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
            if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
            idx = ((*tb - v) / 3);
            if (ok && pick_surface)
              ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                 AtomInfoIsMasked(obj, I->AT[idx]));
            if (ok && !I->VAO.empty()) {
              ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
            }
            if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));

            v_alpha = *(tb++);
            vv = *(tb++);
            if (ok) ok &= CGOAlpha(I->shaderCGO, *v_alpha);
            if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
            if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
            idx = ((*tb - v) / 3);
            if (ok && pick_surface)
              ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                 AtomInfoIsMasked(obj, I->AT[idx]));
            if (ok && !I->VAO.empty()) {
              ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
            }
            if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));

            v_alpha = *(tb++);
            vv = *(tb++);
            if (ok) ok &= CGOAlpha(I->shaderCGO, *v_alpha);
            if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
            if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
            idx = ((*tb - v) / 3);
            if (ok && pick_surface)
              ok &= CGOPickColor(I->shaderCGO, I->AT[idx],
                                 AtomInfoIsMasked(obj, I->AT[idx]));
            if (ok && !I->VAO.empty()) {
              ok &= CGOAccessibility(I->shaderCGO, I->VAO[idx]);
            }
            if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
          }
          if (ok) ok &= CGOEnd(I->shaderCGO);
        }
        FreeP(ix);
        FreeP(z_value);
        FreeP(t_buf);
      } else if (ok) {
        if (info->alpha_cgo) { /* global transparency sort */
          if (I->allVisibleFlag) {
            if (I->oneColorFlag) {
              float col[3];
              ColorGetEncoded(G, I->oneColor, col);

              c = *(s++);
              while (c) {
                int parity = 0;
                s += 2;
                while (ok && c--) {
                  ok &= CGOAlphaTriangle(info->alpha_cgo, v + s[-2] * 3,
                                         v + s[-1] * 3, v + (*s) * 3,
                                         vn + s[-2] * 3, vn + s[-1] * 3,
                                         vn + (*s) * 3, col, col, col, alpha,
                                         alpha, alpha, parity);
                  s++;
                  parity = !parity;
                }
                c = *(s++);
              }
            } else {
              c = *(s++);
              while (ok && c) {
                float *col0, *col1, *col2;
                int parity = 0;
                float alpha0, alpha1, alpha2;
                col0 = vc + (*s) * 3;
                if (va) {
                  alpha0 = va[(*s)];
                } else {
                  alpha0 = alpha;
                }
                s++;
                col1 = vc + (*s) * 3;
                if (va) {
                  alpha1 = va[(*s)];
                } else {
                  alpha1 = alpha;
                }
                s++;
                while (ok && c--) {
                  col2 = vc + (*s) * 3;
                  if (va) {
                    alpha2 = va[(*s)];
                  } else {
                    alpha2 = alpha;
                  }
                  ok &= CGOAlphaTriangle(info->alpha_cgo, v + s[-2] * 3,
                                         v + s[-1] * 3, v + (*s) * 3,
                                         vn + s[-2] * 3, vn + s[-1] * 3,
                                         vn + (*s) * 3, col0, col1, col2,
                                         alpha0, alpha1, alpha2, parity);
                  alpha0 = alpha1;
                  alpha1 = alpha2;
                  col0 = col1;
                  col1 = col2;
                  s++;
                  parity = !parity;
                }
                c = *(s++);
              }
            }
          } else if (ok) { /* subset s */
            c = I->NT;
            if (c) {
              if (I->oneColorFlag) {
                float color[3];
                ColorGetEncoded(G, I->oneColor, color);
                while (ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {
                    ok &= CGOAlphaTriangle(info->alpha_cgo, v + t[0] * 3,
                                           v + t[1] * 3, v + t[2] * 3,
                                           vn + t[0] * 3, vn + t[1] * 3,
                                           vn + t[2] * 3, color, color, color,
                                           alpha, alpha, alpha, 0);
                  }
                  t += 3;
                }
              } else {
                while (ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {
                    if (va) {
                      ok &= CGOAlphaTriangle(
                          info->alpha_cgo, v + t[0] * 3, v + t[1] * 3,
                          v + t[2] * 3, vn + t[0] * 3, vn + t[1] * 3,
                          vn + t[2] * 3, vc + t[0] * 3, vc + t[1] * 3,
                          vc + t[2] * 3, va[t[0]], va[t[1]], va[t[2]], 0);
                    } else {
                      ok &= CGOAlphaTriangle(
                          info->alpha_cgo, v + t[0] * 3, v + t[1] * 3,
                          v + t[2] * 3, vn + t[0] * 3, vn + t[1] * 3,
                          vn + t[2] * 3, vc + t[0] * 3, vc + t[1] * 3,
                          vc + t[2] * 3, alpha, alpha, alpha, 0);
                    }
                  }
                  t += 3;
                }
              }
              if (ok) CGOEnd(info->alpha_cgo);
            }
          }
        } else if (ok) {
          /* fast and ugly */
          /*          glCullFace(GL_BACK);
             glEnable(GL_CULL_FACE);
             glDepthMask(GL_FALSE); */
          if (I->allVisibleFlag) {
            if (I->oneColorFlag) {
              float col[3];
              ColorGetEncoded(G, I->oneColor, col);

              if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
              if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
              c = *(s++);
              while (ok && c) {
                if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
                while (ok && c--) {
                  ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                  if (ok && !I->VAO.empty()) {
                    ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                  }
                  if (ok && pick_surface)
                    ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                       AtomInfoIsMasked(obj, I->AT[*s]));
                  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                  s++;
                }
                if (ok) ok &= CGOEnd(I->shaderCGO);
                c = *(s++);
              }
            } else { /* I->oneColorFlag */
              c = *(s++);
              while (ok && c) {
                float* col;
                if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
                col = vc + (*s) * 3;
                if (va) {
                  if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
                  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                } else {
                  if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
                  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                }
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
                col = vc + (*s) * 3;
                if (va) {
                  if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
                  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                } else {
                  if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
                  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                }
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
                while (ok && c--) {
                  col = vc + (*s) * 3;
                  if (va) {
                    ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
                    if (ok)
                      ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                  } else {
                    ok &= CGOAlpha(I->shaderCGO, alpha);
                    if (ok)
                      ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
                  }
                  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                  if (ok && !I->VAO.empty()) {
                    ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                  }
                  if (ok && pick_surface)
                    ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                       AtomInfoIsMasked(obj, I->AT[*s]));
                  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                  s++;
                }
                if (ok) ok &= CGOEnd(I->shaderCGO);
                c = *(s++);
              }
            }
          } else { /* subset s */
            c = I->NT;
            if (ok && c) {
              if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
              if (I->oneColorFlag) {
                float color[3];
                ColorGetEncoded(G, I->oneColor, color);
                if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
                if (ok)
                  ok &= CGOColor(I->shaderCGO, color[0], color[1], color[2]);
                while (ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {
                    CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                    if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                    if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                  } else
                    t += 3;
                }
              } else {
                float* col;
                while (ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {
                    col = vc + (*t) * 3;
                    if (va) {
                      ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
                    } else {
                      ok &= CGOAlpha(I->shaderCGO, alpha);
                    }
                    if (ok) ok &= CGOColorv(I->shaderCGO, col);
                    if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                    col = vc + (*t) * 3;
                    if (va) {
                      if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
                    } else {
                      if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
                    }
                    if (ok) ok &= CGOColorv(I->shaderCGO, col);
                    if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                    col = vc + (*t) * 3;
                    if (va) {
                      if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
                    } else {
                      if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
                    }
                    if (ok) ok &= CGOColorv(I->shaderCGO, col);
                    if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                    if (ok && !I->VAO.empty()) {
                      ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                    }
                    if (ok && pick_surface)
                      ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                         AtomInfoIsMasked(obj, I->AT[*t]));
                    if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                    t++;
                  } else
                    t += 3;
                }
              }
              if (ok) ok &= CGOEnd(I->shaderCGO);
            }
          }
        }
        /*          glDisable(GL_CULL_FACE);
           glDepthMask(GL_TRUE); */
      }
    } else if (ok) { /* opaque */
      if (I->allVisibleFlag) {
        if (I->oneColorFlag) {
          if (ok) {
            CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
            c = *(s++);
            while (ok && c) {
              if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
              if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                   AtomInfoIsMasked(obj, I->AT[*s]));
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
              }
              if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
              s++;
              if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                   AtomInfoIsMasked(obj, I->AT[*s]));
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
              }
              if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
              s++;
              while (ok && c--) {
                ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
              }
              if (ok) ok &= CGOEnd(I->shaderCGO);
              c = *(s++);
            }
          }
        } else { /* not one color */
          {
            c = *(s++);
            while (ok && c) {
              ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
              if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
              if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
              }
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                   AtomInfoIsMasked(obj, I->AT[*s]));
              if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
              s++;
              if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
              if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
              if (ok && !I->VAO.empty()) {
                ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
              }
              if (ok && pick_surface)
                ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                   AtomInfoIsMasked(obj, I->AT[*s]));
              if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
              s++;
              while (ok && c--) {
                ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*s]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*s],
                                     AtomInfoIsMasked(obj, I->AT[*s]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
                s++;
              }
              if (ok) ok &= CGOEnd(I->shaderCGO);
              c = *(s++);
            }
          }
        } /* one color */
      } else if (ok) { /* subsets */
        c = I->NT;
        if (ok && c) {
          ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
          if (I->oneColorFlag) {
            if (ok) ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
            while (ok && c--) {
              if (visibility_test(I->proximity, vi, t)) {
                ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
              } else
                t += 3;
            }
          } else {
            while (ok && c--) {
              if (visibility_test(I->proximity, vi, t)) {
                ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
                if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
                if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
                if (ok && !I->VAO.empty()) {
                  ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
                }
                if (ok && pick_surface)
                  ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                                     AtomInfoIsMasked(obj, I->AT[*t]));
                if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
                t++;
              } else
                t += 3;
            }
          }
          if (ok) ok &= CGOEnd(I->shaderCGO);
        }
      }
      /*          if (use_shader) {
        CShaderPrg_Disable(shaderPrg);
        }*/
    }
  }
  if (ok && SettingGetGlobal_i(G, cSetting_surface_debug)) {
    t = I->T.data();
    c = I->NT;
    if (c) {
      ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
      while (ok && c--) {
        if (I->allVisibleFlag || visibility_test(I->proximity, vi, t)) {
          ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && !I->VAO.empty()) {
            ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
          }
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
          if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && !I->VAO.empty()) {
            ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
          }
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
          if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && !I->VAO.empty()) {
            ok &= CGOAccessibility(I->shaderCGO, I->VAO[*t]);
          }
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
        } else {
          t += 3;
        }
      }
      if (ok) ok &= CGOEnd(I->shaderCGO);
    }
    t = I->T.data();
    c = I->NT;

    if (ok && c) {
      ok &= CGOColor(I->shaderCGO, 0.0, 1.0, 0.0);
      if (ok) ok &= CGODotwidth(I->shaderCGO, 1.0F);
      while (ok && c--) {
        ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
        if (I->allVisibleFlag || visibility_test(I->proximity, vi, t)) {
          if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
          if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
          if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
          if (ok && pick_surface)
            ok &= CGOPickColor(I->shaderCGO, I->AT[*t],
                               AtomInfoIsMasked(obj, I->AT[*t]));
          if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
          t++;
        } else {
          t += 3;
        }
        if (ok) ok &= CGOEnd(I->shaderCGO);
      }
    }

    c = I->N;
    if (ok && c) {
      ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
      if (ok) ok &= CGOResetNormal(I->shaderCGO, true);
      if (ok) ok &= CGOBegin(I->shaderCGO, GL_LINES);
      while (ok && c--) {
        ok &= CGOVertexv(I->shaderCGO, v);
        if (ok)
          ok &= CGOVertex(I->shaderCGO, v[0] + vn[0] / 2, v[1] + vn[1] / 2,
                          v[2] + vn[2] / 2);
        v += 3;
        vn += 3;
      }
      if (ok) ok &= CGOEnd(I->shaderCGO);
    }
  }

  if (ok) ok &= CGOStop(I->shaderCGO);
  if (I->Type != SurfaceType::Triangle) {
    CGOCombineBeginEnd(&I->shaderCGO);
  }
  if (I->Type == SurfaceType::DotDefault) {
    if (dot_as_spheres) {
      CGO* tmpCGO = CGONew(G);
      CHECKOK(ok, tmpCGO);
      ok &= CGOEnable(tmpCGO, GL_SPHERE_SHADER);
      ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);  // TODO: this needs normals in
                                                 // sphere shader PYMOL-1870
      ok &= CGOSpecial(tmpCGO, DOTSIZE_WITH_SPHERESCALE);
      convertcgo = CGOOptimizeSpheresToVBONonIndexedNoShader(
          I->shaderCGO,
          CGO_BOUNDING_BOX_SZ + fsizeof<cgo::draw::sphere_buffers>() + 2);
      ok &= CGOAppendNoStop(tmpCGO, convertcgo);
      CGOFreeWithoutVBOs(convertcgo);
      ok &= CGODisable(tmpCGO, GL_SPHERE_SHADER);
      CGOStop(tmpCGO);
      convertcgo = tmpCGO;
    } else {
      convertcgo = CGOOptimizeToVBONotIndexedNoShader(I->shaderCGO);
    }
    CHECKOK(ok, convertcgo);
    if (ok) convertcgo->use_shader = true;
  } else if (I->Type == SurfaceType::Triangle) {
    CGO* tmpCGO = CGONew(G);
    CHECKOK(ok, tmpCGO);
    auto convertcgo2 = CGOConvertLinesToShaderCylinders(I->shaderCGO, 0);
    CHECKOK(ok, convertcgo2);
    CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
    CGOSpecial(tmpCGO, MESH_WIDTH_FOR_SURFACES);
    convertcgo = CGOConvertShaderCylindersToCylinderShader(convertcgo2, tmpCGO);
    CGOAppendNoStop(tmpCGO, convertcgo);
    CGOFreeWithoutVBOs(convertcgo);
    CGODisable(tmpCGO, GL_CYLINDER_SHADER);
    CGOStop(tmpCGO);
    convertcgo = tmpCGO;
    convertcgo->use_shader = true;
    CGOFree(convertcgo2);
  } else {
    if ((alpha != 1.0) || va) {  // semi-transparent
      if (ok)
        convertcgo = CGOOptimizeToVBOIndexedWithColorEmbedTransparentInfo(
            I->shaderCGO, 0, 0, 0);
      CHECKOK(ok, convertcgo);
    } else {
      if (ok) convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0, false);
      CHECKOK(ok, convertcgo);
    }
    if (ok) convertcgo->use_shader = true;
    {
      CGO* tmpCGO = nullptr;
      tmpCGO = CGONew(G);
      CGOEnable(tmpCGO, GL_SURFACE_SHADER);
      // CGOEnable(tmpCGO, CGO_GL_LIGHTING); // do we need this?
      CGOSpecial(tmpCGO, SET_SURFACE_UNIFORMS);
      CGOAppendNoStop(tmpCGO, convertcgo);
      CGODisable(tmpCGO, GL_SURFACE_SHADER);
      CGOStop(tmpCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = tmpCGO;
      convertcgo->use_shader = true;
    }
  }

  if (ok && I->Type == SurfaceType::DotDefault && !dot_as_spheres) {
    setShaderCGO(I, CGONew(G));
    CHECKOK(ok, I->shaderCGO);
    if (ok) {
      I->shaderCGO->use_shader = true;
      if (ok) ok &= CGOResetNormal(I->shaderCGO, true);
      if (ok) ok &= CGOEnable(I->shaderCGO, GL_SURFACE_SHADER);
      if (ok) ok &= CGOSpecial(I->shaderCGO, POINTSIZE_DYNAMIC_DOT_WIDTH);
      if (ok) CGOAppendNoStop(I->shaderCGO, convertcgo);
      if (ok) ok &= CGODisable(I->shaderCGO, GL_SURFACE_SHADER);
      if (ok) ok &= CGOStop(I->shaderCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = nullptr;
    }
  } else {
    setShaderCGO(I, convertcgo);
    convertcgo = nullptr;
  }

  if (I->shaderCGO) {
    I->shaderCGO->no_pick = !pick_surface;
    setPickingCGO(I, I->shaderCGO);
  }

  return ok;
}

void RepSurface::render(RenderInfo* info) {
  auto* I = this;
  CRay* ray = info->ray;
  auto pick = info->pick;
  float* v = I->V.data();
  float* vn = I->VN.data();
  float* vc = I->VC.data();
  float* va = I->VA.data();
  int* rc = I->RC.data();
  int* t = I->T.data();
  int* s = I->S.data();
  int c = I->N;
  int* vi = I->Vis.data();
  int ok = true;
  float ambient_occlusion_scale = 0.f;
  int ambient_occlusion_mode =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                   cSetting_ambient_occlusion_mode);
  int ambient_occlusion_mode_div_4 = 0;
  if (ambient_occlusion_mode) {
    ambient_occlusion_scale =
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_ambient_occlusion_scale);
    ambient_occlusion_mode_div_4 = ambient_occlusion_mode / 4;
  }

  if ((I->Type != SurfaceType::DotDefault) && (!s)) {
    return;
  }

  auto alpha = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_transparency);
  alpha = 1.0F - alpha;
  if (fabs(alpha - 1.0) < R_SMALL4) alpha = 1.0F;
  if (ray) {
    ray->transparentf(1.0F - alpha);
    if (I->Type == SurfaceType::DotDefault) {
      /* dot surface */
      float radius;
      radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_dot_radius);
      if (radius == 0.0F) {
        radius = ray->PixelRadius *
                 SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_dot_width) /
                 1.4142F;
      }

      if (I->oneColorFlag) {
        float col[3];
        ColorGetEncoded(G, I->oneColor, col);
        ray->color3fv(col);
      }

      if (c)
        while (ok && c--) {
          if (*vi) {
            if (!I->oneColorFlag) {
              ray->color3fv(vc);
            }
            ok &= ray->sphere3fv(v, radius);
          }
          vi++;
          vc += 3;
          v += 3;
        }
    } else if ((I->Type == SurfaceType::Solid) ||
               (I->Type == SurfaceType::SolidDeterministic) ||
               (I->Type == SurfaceType::SolidWeightings) ||
               (I->Type == SurfaceType::SolidScribed)) { /* solid surface */
      c = I->NT;

      if (I->oneColorFlag) {
        float col[3], col1[3], col2[3], col3[3];
        ColorGetEncoded(G, I->oneColor, col);
        while (ok && c--) {
          if (visibility_test(I->proximity, vi, t)) {
            copy3f(col, col1);
            copy3f(col, col2);
            copy3f(col, col3);
            if (!I->VAO.empty()) {
              float ao1, ao2, ao3;
              switch (ambient_occlusion_mode_div_4) {
                case 1:
                  ao1 = 1.f - ambient_occlusion_scale * (I->VAO[*t]);
                  ao2 = 1.f - ambient_occlusion_scale * (I->VAO[*(t + 1)]);
                  ao3 = 1.f - ambient_occlusion_scale * (I->VAO[*(t + 2)]);
                  break;
                case 2:
                  ao1 =
                      cos(.5f * PI *
                          CLAMP_VALUE(ambient_occlusion_scale * (I->VAO[*t])));
                  ao2 = cos(.5f * PI *
                            CLAMP_VALUE(ambient_occlusion_scale *
                                        (I->VAO[*(t + 1)])));
                  ao3 = cos(.5f * PI *
                            CLAMP_VALUE(ambient_occlusion_scale *
                                        (I->VAO[*(t + 2)])));
                  break;
                default:
                  ao1 = CLAMP_VALUE(
                      1.f /
                      (1.f +
                       exp(.5f *
                           ((ambient_occlusion_scale * (I->VAO[*t])) - 10.f))));
                  ao2 = CLAMP_VALUE(1.f /
                                    (1.f + exp(.5f * ((ambient_occlusion_scale *
                                                       (I->VAO[*(t + 1)])) -
                                                      10.f))));
                  ao3 = CLAMP_VALUE(1.f /
                                    (1.f + exp(.5f * ((ambient_occlusion_scale *
                                                       (I->VAO[*(t + 2)])) -
                                                      10.f))));
              }
              mult3f(col1, ao1, col1);
              mult3f(col2, ao2, col2);
              mult3f(col3, ao3, col3);
            }
            ok &= ray->triangle3fv(v + (*t) * 3, v + (*(t + 1)) * 3,
                                   v + (*(t + 2)) * 3, vn + (*t) * 3,
                                   vn + (*(t + 1)) * 3, vn + (*(t + 2)) * 3,
                                   col1, col2, col3);
          }
          t += 3;
        }
      } else {
        while (ok && c--) {
          int ttA = *t, ttB = *(t + 1), ttC = *(t + 2);
          if (visibility_test(I->proximity, vi, t)) {
            int ttA3 = ttA * 3, ttB3 = ttB * 3, ttC3 = ttC * 3;
            float cA[3], cB[3], cC[3];
            copy3f(vc + ttA3, cA);
            copy3f(vc + ttB3, cB);
            copy3f(vc + ttC3, cC);
            //            register float *cA = vc + ttA3, *cB = vc + ttB3, *cC =
            //            vc + ttC3;
            if (rc) {
              if (rc[ttA] < -1) ColorGetEncoded(G, rc[ttA], cA);
              //                ColorGetEncoded(G, rc[ttA], (cA = colA));
              if (rc[ttB] < -1) ColorGetEncoded(G, rc[ttB], cB);
              //                ColorGetEncoded(G, rc[ttB], (cB = colB));
              if (rc[ttC] < -1) ColorGetEncoded(G, rc[ttC], cC);
              //                ColorGetEncoded(G, rc[ttC], (cC = colC));
            }
            if ((*(vi + ttA)) || (*(vi + ttB)) || (*(vi + ttC))) {
              if (!I->VAO.empty()) {
                float ao1, ao2, ao3;
                switch (ambient_occlusion_mode_div_4) {
                  case 1:
                    ao1 = 1.f - ambient_occlusion_scale * (I->VAO[*t]);
                    ao2 = 1.f - ambient_occlusion_scale * (I->VAO[*(t + 1)]);
                    ao3 = 1.f - ambient_occlusion_scale * (I->VAO[*(t + 2)]);
                    break;
                  case 2:
                    ao1 = cos(
                        .5f * PI *
                        CLAMP_VALUE(ambient_occlusion_scale * (I->VAO[*t])));
                    ao2 = cos(.5f * PI *
                              CLAMP_VALUE(ambient_occlusion_scale *
                                          (I->VAO[*(t + 1)])));
                    ao3 = cos(.5f * PI *
                              CLAMP_VALUE(ambient_occlusion_scale *
                                          (I->VAO[*(t + 2)])));
                    break;
                  default:
                    ao1 = CLAMP_VALUE(
                        1.f / (1.f + exp(.5f * ((ambient_occlusion_scale *
                                                 (I->VAO[*t])) -
                                                10.f))));
                    ao2 = CLAMP_VALUE(
                        1.f / (1.f + exp(.5f * ((ambient_occlusion_scale *
                                                 (I->VAO[*(t + 1)])) -
                                                10.f))));
                    ao3 = CLAMP_VALUE(
                        1.f / (1.f + exp(.5f * ((ambient_occlusion_scale *
                                                 (I->VAO[*(t + 2)])) -
                                                10.f))));
                }
                mult3f(cA, ao1, cA);
                mult3f(cB, ao2, cB);
                mult3f(cC, ao3, cC);
              }
              if (va) {
                ok &= ray->triangleTrans3fv(v + ttA3, v + ttB3, v + ttC3,
                                            vn + ttA3, vn + ttB3, vn + ttC3, cA,
                                            cB, cC, 1.0F - va[ttA],
                                            1.0F - va[ttB], 1.0F - va[ttC]);
              } else {
                ok &= ray->triangle3fv(v + ttA3, v + ttB3, v + ttC3, vn + ttA3,
                                       vn + ttB3, vn + ttC3, cA, cB, cC);
              }
            }
          }
          t += 3;
        }
      }
    } else if (I->Type == SurfaceType::Triangle) {
      float radius;
      int t0, t1, t2;
      int spacing = 10;
      int* cache = pymol::calloc<int>(spacing * (I->N + 1));
      CHECKOK(ok, cache);

      radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_mesh_radius);

      if (ok && radius == 0.0F) {
        float line_width = SettingGet_f(
            G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_width);
        line_width = SceneGetDynamicLineWidth(info, line_width);

        radius = ray->PixelRadius * line_width / 2.0F;
      }
      if (ok) {
        float col[3], *col0, *col1, *col2;
        c = I->NT;
        if (I->oneColorFlag) {
          ColorGetEncoded(G, I->oneColor, col);
          col0 = col1 = col2 = col;
        }
        while (ok && c--) {
          t0 = (*t);
          t1 = (*(t + 1));
          t2 = (*(t + 2));
          if (visibility_test(I->proximity, vi, t)) {
            if (!I->oneColorFlag) {
              col0 = vc + t0 * 3;
              col1 = vc + t1 * 3;
              col2 = vc + t2 * 3;
            }
            if (!check_and_add(cache, spacing, t0, t1))
              ok &= ray->sausage3fv(v + t0 * 3, v + t1 * 3, radius, col0, col1);
            if (!check_and_add(cache, spacing, t1, t2))
              ok &= ray->sausage3fv(v + t1 * 3, v + t2 * 3, radius, col1, col2);
            if (!check_and_add(cache, spacing, t2, t0))
              ok &= ray->sausage3fv(v + t2 * 3, v + t0 * 3, radius, col2, col0);
          }
          t += 3;
        }
        FreeP(cache);
      }
    }
    if (ok) {
      ray->transparentf(0.0);
    } else {
      /* If not ok, then Clear Entire RepSurface, not just the ray object */
    }
  } else if (G->HaveGUI && G->ValidContext) {
    /* Not ray tracing, but rendering */
    if (pick) {
      // Don't render transparent unpickable surfaces. Do render solid but
      // unpickable surfaces to write the depth buffer and prevent
      // through-picking.
      if (I->pickingCGO && !(I->pickingCGO->no_pick && alpha < 1.F)) {
        CGORenderPicking(I->pickingCGO, info, &I->context, cs->Setting.get(),
                         obj->Setting.get());
      }
    } else {
      {
        bool dot_as_spheres = SettingGet_b(
            G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_as_spheres);

        if (!I->shaderCGO || CGOCheckWhetherToFree(G, I->shaderCGO) ||
            (I->Type == SurfaceType::DotDefault &&
             I->dot_as_spheres != dot_as_spheres)) {
          ok &= RepSurfaceCGOGenerate(I, info);
        }

        if (ok && I->shaderCGO) {
          const float* color = ColorGet(G, obj->Color);
          CGORender(I->shaderCGO, color, nullptr, nullptr, info, I);
          return;
        }

        PRINTFB(G, FB_RepSurface, FB_Errors)
        " RepSurfaceCGOGenerate failed\n" ENDFB(G);
      }

      setShaderCGO(I, nullptr);
    }
  }
  if (!ok) {
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepSurface] = false;
  }
}

bool RepSurface::sameVis() const {
  for (int idx = 0; idx < cs->NIndex; ++idx) {
    auto* ai = cs->getAtomInfo(idx);
    if (LastVisib[idx] != GET_BIT(ai->visRep, cRepSurface)) {
      return false;
    }
  }

  return true;
}

void RepSurface::invalidate(cRepInv_t level) {
  Rep::invalidate(level);

  if (level >= cRepInvColor) ColorInvalidated = true;
}

static void RepSurfaceSetSettings(PyMOLGlobals* G, CoordSet* cs,
                                  ObjectMolecule* obj, int surface_quality,
                                  SurfaceType surface_type, float* point_sep,
                                  int* sphere_idx, int* solv_sph_idx,
                                  int* circumscribe) {
  if (surface_quality >= 4) { /* totally impractical */
    *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_surface_best) /
                 4.f;
    *sphere_idx = 4;
    *solv_sph_idx = 4;
    if (*circumscribe < 0) *circumscribe = 91;
  } else {
    switch (surface_quality) {
      case 3: /* nearly impractical */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_best) /
                     3.f;
        *sphere_idx = 4;
        *solv_sph_idx = 3;
        if (*circumscribe < 0) *circumscribe = 71;
        break;
      case 2:
        /* nearly perfect */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_best) /
                     2.f;
        *sphere_idx = 3;
        *solv_sph_idx = 3;
        if (*circumscribe < 0) *circumscribe = 41;
        break;
      case 1:
        /* good */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_best);
        *sphere_idx = 2;
        *solv_sph_idx = 3;
        if ((*circumscribe < 0) &&
            (surface_type == SurfaceType::SolidEmpirical))
          *circumscribe = 40;
        break;
      case 0:
        /* 0 - normal */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_normal);
        *sphere_idx = 1;
        *solv_sph_idx = 2;
        if ((*circumscribe < 0) &&
            (surface_type == SurfaceType::SolidEmpirical))
          *circumscribe = 30;
        break;
      case -1:
        /* -1 */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_poor);
        *sphere_idx = 1;
        *solv_sph_idx = 2;
        if ((*circumscribe < 0) &&
            (surface_type == SurfaceType::SolidEmpirical))
          *circumscribe = 10;
        break;
      case -2:
        /* -2 god awful */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_poor) *
                     1.5F;
        *sphere_idx = 1;
        *solv_sph_idx = 1;
        break;
      case -3:
        /* -3 miserable */
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_miserable);
        *sphere_idx = 1;
        *solv_sph_idx = 1;
        break;
      default:
        *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_miserable) *
                     1.18F;
        *sphere_idx = 0;
        *solv_sph_idx = 1;
    }
  }
  /* Fixed problem with surface holes when surface_quality>2, it seems like
     circumscribe can only be used with surface_solvent */
  if ((*circumscribe < 0) ||
      (!SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_surface_solvent)))
    *circumscribe = 0;
}

static int RepSurfacePrepareSurfaceJob(
    PyMOLGlobals* G, SurfaceJob* surf_job, RepSurface* I, CoordSet* cs,
    ObjectMolecule* obj, std::vector<SurfaceJobAtomInfo>& atom_info,
    float* carve_vla, int n_present, std::vector<int>& present_vla,
    int optimize, int sphere_idx, int solv_sph_idx, SurfaceType surface_type,
    int circumscribe, float probe_radius, float point_sep, float carve_cutoff) {
  int ok = true;
  surf_job->maxVdw = I->max_vdw;
  surf_job->allVisibleFlag = I->allVisibleFlag;

  surf_job->atomInfo = std::move(atom_info);

  surf_job->nPresent = n_present;
  if (!present_vla.empty() && optimize) {
    /* implies that n_present < cs->NIndex, so eliminate
       irrelevant atoms & coordinates if we are optimizing subsets */
    surf_job->coord.resize(n_present * 3);
    CHECKOK(ok, !surf_job->coord.empty());
    if (ok) {
      int* p = present_vla.data();
      SurfaceJobAtomInfo* ai_src = surf_job->atomInfo.data();
      SurfaceJobAtomInfo* ai_dst = surf_job->atomInfo.data();
      const float* v_src = cs->Coord.data();
      float* v_dst = surf_job->coord.data();
      int a;

      for (a = 0; a < cs->NIndex; a++) {
        if (*(p++)) {
          copy3f(v_src, v_dst);
          v_dst += 3;
          if (ai_dst != ai_src) *ai_dst = *ai_src;
          ai_dst++;
        }
        v_src += 3;
        ai_src++;
      }
    }
    surf_job->atomInfo.resize(n_present);
    CHECKOK(ok, !surf_job->atomInfo.empty());
  } else {
    surf_job->presentVla = std::move(present_vla);
    surf_job->coord.resize(cs->NIndex * 3);
    CHECKOK(ok, !surf_job->coord.empty());
    if (ok)
      std::copy_n(cs->Coord.data(), 3 * cs->NIndex, surf_job->coord.data());
  }
  if (ok) {
    surf_job->sphereIndex = sphere_idx;
    surf_job->solventSphereIndex = solv_sph_idx;

    surf_job->surfaceType = surface_type;
    surf_job->circumscribe = circumscribe;
    surf_job->probeRadius = probe_radius;
    surf_job->pointSep = point_sep;

    surf_job->trimCutoff = SettingGet_f(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_trim_cutoff);
    surf_job->trimFactor = SettingGet_f(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_trim_factor);

    surf_job->cavityMode = SettingGet_i(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_cavity_mode);
    surf_job->cavityRadius =
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_surface_cavity_radius);
    surf_job->cavityCutoff =
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_surface_cavity_cutoff);
    if (carve_vla) {
      surf_job->carveVla = std::vector<float>(VLAGetSize(carve_vla));
      std::copy_n(carve_vla, surf_job->carveVla.size(),
                  surf_job->carveVla.data());
    }
    surf_job->carveCutoff = carve_cutoff;
    surf_job->surfaceMode = SettingGet_i(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_mode);
    surf_job->surfaceSolvent = SettingGet_b(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_solvent);
    surf_job->cavityCull = SettingGet_i(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_cavity_cull);
  }
  return ok;
}

static void RepSurfaceConvertSurfaceJobToPyObject(
    PyMOLGlobals* G, SurfaceJob* surf_job, CoordSet* cs, ObjectMolecule* obj,
    PyObject** entry, PyObject** input, PyObject** output, int* found) {
  int cache_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_cache_mode);

  if (cache_mode > 0) {
    int blocked = PAutoBlock(G);
    *input = SurfaceJobInputAsTuple(G, surf_job);

    if (PCacheGet(G, output, entry, *input) == OV_STATUS_YES) {
      if (OV_OK(SurfaceJobResultFromTuple(G, surf_job, *output))) {
        *found = true;
        PXDecRef(*input);
        *input = nullptr;
        PXDecRef(*entry);
        *entry = nullptr;
      }
      PXDecRef(*output);
      *output = nullptr;
    }
    if (PyErr_Occurred()) PyErr_Print();
    PAutoUnblock(G, blocked);
  }
}

static void RepSurfaceFindAllPresentAtoms(ObjectMolecule* obj, CoordSet* cs,
                                          int* present_vla, int inclH,
                                          int cullByFlag) {
  int* ap = present_vla;
  const int* idx_to_atm = cs->IdxToAtm.data();
  const AtomInfoType* obj_AtomInfo = obj->AtomInfo.data();
  int a, cs_NIndex = cs->NIndex;
  for (a = 0; a < cs_NIndex; a++) {
    const AtomInfoType* ai1 = obj_AtomInfo + *(idx_to_atm++);
    if ((ai1->visRep & cRepSurfaceBit) && (inclH || (!ai1->isHydrogen())) &&
        ((!cullByFlag) ||
         (!(ai1->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
      *ap = 2;
    else
      *ap = 0;
    ap++;
  }
}

static int RepSurfaceAddNearByAtomsIfNotSurfaced(PyMOLGlobals* G, MapType* map,
                                                 ObjectMolecule* obj,
                                                 CoordSet* cs, int* present_vla,
                                                 int inclH, int cullByFlag,
                                                 float probe_radius,
                                                 int optimize) {
  int ok = true;
  float probe_radiusX2 = probe_radius * 2;
  int a;
  for (a = 0; ok && a < cs->NIndex; a++) {
    if (!present_vla[a]) {
      AtomInfoType* ai1 = obj->AtomInfo + cs->IdxToAtm[a];
      if ((inclH || (!ai1->isHydrogen())) &&
          ((!cullByFlag) || !(ai1->flags & cAtomFlag_ignore))) {
        const float* v0 = cs->coordPtr(a);
        int i = *(MapLocusEStart(map, v0));
        if (optimize) {
          if (i && !map->EList.empty()) {
            int j = map->EList[i++];
            while (j >= 0) {
              if (present_vla[j] > 1) {
                AtomInfoType* ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                if (within3f(cs->coordPtr(j), v0,
                             ai1->vdw + ai2->vdw + probe_radiusX2)) {
                  present_vla[a] = 1;
                  break;
                }
              }
              j = map->EList[i++];
            }
          }
        } else
          present_vla[a] = 1;
      }
    }
    ok &= !G->Interrupt;
  }
  return ok;
}

static int RepSurfaceRemoveAtomsNotWithinCutoff(PyMOLGlobals* G,
                                                MapType* carve_map,
                                                float* carve_vla, CoordSet* cs,
                                                int* present_vla,
                                                float carve_cutoff) {
  int ok = true;
  int a;
  for (a = 0; ok && a < cs->NIndex; a++) {
    int include_flag = false;
    if (carve_map) {
      const float* v0 = cs->coordPtr(a);
      int i = *(MapLocusEStart(carve_map, v0));
      if (i && !carve_map->EList.empty()) {
        int j = carve_map->EList[i++];
        while (j >= 0) {
          if (within3f(carve_vla + 3 * j, v0, carve_cutoff)) {
            include_flag = true;
            break;
          }
          j = carve_map->EList[i++];
        }
      }
    }
    if (!include_flag) present_vla[a] = 0;
    ok &= !G->Interrupt;
  }
  return ok;
}

Rep* RepSurfaceNew(CoordSet* cs, int state) {
  // The app's progress bar follows this build's time-weighted steps (about
  // 3% solvent dots, 62% solvent pass, 5% clean-up, 27% triangulation) rather
  // than OrthoBusyFast's even fifths.
  struct ExplicitProgress {
    ExplicitProgress() {
      g_mmol_busy_explicit = true;
      BusyProgress(0.0f);
    }
    ~ExplicitProgress() {
      BusyProgress(1.0f);
      g_mmol_busy_explicit = false;
    }
  } explicit_progress;
  int ok = true;
  PyMOLGlobals* G = cs->G;
  ObjectMolecule* obj = cs->Obj;

  int surface_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_surface_mode);
  int cullByFlag = (surface_mode == cRepSurface_by_flags);
  int inclH = !((surface_mode == cRepSurface_heavy_atoms) ||
                (surface_mode == cRepSurface_vis_heavy_only));
  int inclInvis = !((surface_mode == cRepSurface_vis_only) ||
                    (surface_mode == cRepSurface_vis_heavy_only));
  int visFlag = false;

  if (GET_BIT(obj->RepVisCache, cRepSurface)) {
    const int* idx_to_atm = cs->IdxToAtm.data();
    const AtomInfoType* obj_AtomInfo = obj->AtomInfo.data();
    int a, cs_NIndex = cs->NIndex;

    // If all atoms have "flag ignore", then don't cull by flags
    if (cullByFlag) {
      bool all_ignore = true;

      for (int idx = 0; idx < cs_NIndex; ++idx) {
        if (!(obj_AtomInfo[idx_to_atm[idx]].flags & cAtomFlag_ignore)) {
          all_ignore = false;
          break;
        }
      }

      if (all_ignore) {
        cullByFlag = false;
        surface_mode = cRepSurface_all;
      }
    }

    for (a = 0; a < cs_NIndex; a++) {
      const AtomInfoType* ai1 = obj_AtomInfo + *(idx_to_atm++);
      if ((ai1->visRep & cRepSurfaceBit) && (inclH || (!ai1->isHydrogen())) &&
          ((!cullByFlag) ||
           (!(ai1->flags & (cAtomFlag_exfoliate | cAtomFlag_ignore))))) {
        visFlag = true;
        break;
      }
    }
  }
  if (!visFlag) {
    return (nullptr); /* skip if no thing visible */
  }

  auto I = new RepSurface(cs, state);
  I->surface_mode = surface_mode;

  {
    int surface_flag = false;
    auto surface_type = SettingGet<SurfaceType>(
        G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_type);
    int surface_quality = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                       cSetting_surface_quality);
    float probe_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                      cSetting_solvent_radius);
    int optimize = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_surface_optimize_subsets);
    int circumscribe = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                    cSetting_surface_circumscribe);
    int sphere_idx = 0, solv_sph_idx = 0;
    MapType* map = nullptr;
    float point_sep;
    int n_present = 0;
    int carve_state = 0;
    int carve_flag = false;
    float carve_cutoff;
    const char* carve_selection = nullptr;
    float* carve_vla = nullptr;
    MapType* carve_map = nullptr;
    bool smooth_edges = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_surface_smooth_edges);

    I->Type = surface_type;

    I->max_vdw = ObjectMoleculeGetMaxVDW(obj);

    RepSurfaceSetSettings(G, cs, obj, surface_quality, surface_type, &point_sep,
                          &sphere_idx, &solv_sph_idx, &circumscribe);

    I->allVisibleFlag = true;
    obj = cs->Obj;

    /* don't waist time computing a Surface unless we need it!! */
    {
      const int* idx_to_atm = cs->IdxToAtm.data();
      const AtomInfoType* obj_AtomInfo = obj->AtomInfo.data();
      int a, cs_NIndex = cs->NIndex;
      for (a = 0; a < cs_NIndex; a++) {
        const AtomInfoType* ai1 = obj_AtomInfo + *(idx_to_atm++);
        if ((ai1->visRep & cRepSurfaceBit) && (inclH || (!ai1->isHydrogen())) &&
            ((!cullByFlag) ||
             (!(ai1->flags & (cAtomFlag_exfoliate | cAtomFlag_ignore)))))
          surface_flag = true;
        else {
          I->allVisibleFlag = false;
          if (surface_flag) break;
        }
      }
    }
    if (surface_flag) {
      std::vector<SurfaceJobAtomInfo> atom_info(cs->NIndex);
      CHECKOK(ok, !atom_info.empty());
      if (ok && !atom_info.empty()) {
        const AtomInfoType *i_ai, *obj_atom_info = obj->AtomInfo.data();
        const int* idx_to_atm = cs->IdxToAtm.data();
        int n_index = cs->NIndex;
        SurfaceJobAtomInfo* i_atom_info = atom_info.data();
        int i;
        /* fill in surfacing flags into SurfaceJobAtomInfo array */
        for (i = 0; i < n_index; i++) {
          i_ai = obj_atom_info + idx_to_atm[i];
          i_atom_info->flags =
              i_ai->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate);
          i_atom_info->vdw = i_ai->vdw;
          i_atom_info++;
        }
      }

      OrthoBusyFast(G, 0, 1);

      if (ok) {
        n_present = cs->NIndex;
        carve_selection = SettingGet_s(G, cs->Setting.get(), obj->Setting.get(),
                                       cSetting_surface_carve_selection);
        carve_cutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                    cSetting_surface_carve_cutoff);
        if ((!carve_selection) || (!carve_selection[0])) carve_cutoff = 0.0F;
        if (carve_cutoff > 0.0F) {
          carve_state = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_surface_carve_state) -
                        1;
          carve_cutoff += 2 * I->max_vdw + probe_radius;

          if (carve_selection)
            carve_map = SelectorGetSpacialMapFromSeleCoord(
                G, SelectorIndexByName(G, carve_selection), carve_state,
                carve_cutoff, &carve_vla);
          if (carve_map) ok &= MapSetupExpress(carve_map);
          carve_flag = true;
          I->allVisibleFlag = false;
        }
      }
      std::vector<int> present_vla;
      if (ok && !I->allVisibleFlag) {
        /* optimize the space over which we calculate a surface */
        /* first find out which atoms are actually to be surfaced */
        present_vla.resize(cs->NIndex);
        CHECKOK(ok, !present_vla.empty());
        if (ok) {
          RepSurfaceFindAllPresentAtoms(obj, cs, present_vla.data(), inclH,
                                        cullByFlag);
        }

        if (ok)
          map = new MapType(G, 2 * I->max_vdw + probe_radius, cs->Coord,
                            cs->NIndex, nullptr, present_vla.data());
        CHECKOK(ok, map);
        if (ok) ok &= MapSetupExpress(map);

        if (ok && inclInvis) {
          /* then add in the nearby atoms which are not surfaced and not ignored
           */
          ok &= RepSurfaceAddNearByAtomsIfNotSurfaced(
              G, map, obj, cs, present_vla.data(), inclH, cullByFlag,
              probe_radius, optimize);
        }

        if (ok && carve_flag && (!optimize)) {
          /* and optimize for carved region */
          ok &= RepSurfaceRemoveAtomsNotWithinCutoff(
              G, carve_map, carve_vla, cs, present_vla.data(), carve_cutoff);
        }
        MapFree(map);
        map = nullptr;

        /* now count how many atoms we actually need to think about */
        n_present = 0;
        if (ok) {
          int a;
          for (a = 0; a < cs->NIndex; a++) {
            if (present_vla[a]) {
              n_present++;
            }
          }
        }
      }

      if (ok) {
        auto surf_jobU = std::make_unique<SurfaceJob>();
        auto surf_job = surf_jobU.get();
        CHECKOK(ok, surf_job);
        ok &= RepSurfacePrepareSurfaceJob(
            G, surf_job, I, cs, obj, atom_info, carve_vla, n_present,
            present_vla, optimize, sphere_idx, solv_sph_idx, surface_type,
            circumscribe, probe_radius, point_sep, carve_cutoff);

        ok &= !G->Interrupt;

        if (ok) {
          int found = false;
          PyObject* entry = nullptr;
          PyObject* output = nullptr;
          PyObject* input = nullptr;
          int cache_mode = SettingGet_i(
              G, cs->Setting.get(), obj->Setting.get(), cSetting_cache_mode);
          RepSurfaceConvertSurfaceJobToPyObject(G, surf_job, cs, obj, &entry,
                                                &input, &output, &found);
          if (ok && !found) {
            ok &= SurfaceJobRun(G, surf_job);

            if (cache_mode > 1) {
              int blocked = PAutoBlock(G);
              output = SurfaceJobResultAsTuple(G, surf_job);
              PCacheSet(G, entry, output);
              PXDecRef(entry);
              entry = nullptr;
              PXDecRef(output);
              output = nullptr;
              PXDecRef(input);
              input = nullptr;
              PAutoUnblock(G, blocked);
            }
          }
          if (entry || input || output) {
            int blocked = PAutoBlock(G);
            PXDecRef(entry);
            PXDecRef(input);
            PXDecRef(output);
            PAutoUnblock(G, blocked);
          }
        }
        /* surf_job must be valid at this point */
        if (ok) {
          I->N = std::move(surf_job->N);
          I->V = std::move(surf_job->V);
          I->VN = std::move(surf_job->VN);
          I->NT = surf_job->NT;
          I->T = std::move(surf_job->T);
          I->S = std::move(surf_job->S);
        }
      }

      ok &= !G->Interrupt;
      if (ok) I->recolor();

      if (ok && smooth_edges) RepSurfaceSmoothEdges(I);
    }
    if (carve_map) MapFree(carve_map);
    VLAFreeP(carve_vla);
    OrthoBusyFast(G, 4, 4);
  }

  if (!ok) {
    delete I;
    I = nullptr;
  }
  return (Rep*)I;
}

void RepSurfaceSmoothEdges(RepSurface* I) {
  if (I->allVisibleFlag) return;

  std::vector<std::vector<int>> edges(I->N);
  // Find all edges.
  // Uses the following edge structure:
  // edge[v0] -> { v1, t0, t1 }
  // where (v0, v1) defines an edge and (t0, t1)
  // are triangles the edge belongs to.
  // The edge is unique when t0 == t1.
  int* t = I->T.data();
  int* vi = I->Vis.data();
  for (auto i = 0; i < I->NT; i++) {
    if (visibility_test(I->proximity, vi, t)) {
      int v0 = t[0];
      int v1 = t[1];
      int v2 = t[2];
      if (v0 < v1) {
        edges[v0].push_back(v1);
        edges[v0].push_back(i);
        edges[v0].push_back(i);
      } else {
        edges[v1].push_back(v0);
        edges[v1].push_back(i);
        edges[v1].push_back(i);
      }
      if (v1 < v2) {
        edges[v1].push_back(v2);
        edges[v1].push_back(i);
        edges[v1].push_back(i);
      } else {
        edges[v2].push_back(v1);
        edges[v2].push_back(i);
        edges[v2].push_back(i);
      }
      if (v2 < v0) {
        edges[v2].push_back(v0);
        edges[v2].push_back(i);
        edges[v2].push_back(i);
      } else {
        edges[v0].push_back(v2);
        edges[v0].push_back(i);
        edges[v0].push_back(i);
      }
    }
    t += 3;
  }
  // Assign triangles to every edge
  for (auto i = 0; i < I->N; i++) {
    for (auto j = 0; j < edges[i].size(); j += 3) {
      for (auto k = 0; k < j; k += 3) {
        if ((j != k) && (edges[i][j] == edges[i][k])) {
          edges[i][j + 2] = edges[i][k + 1];
          edges[i][k + 2] = edges[i][j + 1];
        }
      }
    }
  }
  std::vector<int> unique_edges;
  // Build a list of unique edges
  for (auto i = 0; i < I->N; i++) {
    for (auto j = 0; j < edges[i].size(); j += 3) {
      if (edges[i][j + 1] == edges[i][j + 2]) {
        unique_edges.push_back(i);
        unique_edges.push_back(edges[i][j]);
      }
    }
  }
  if (unique_edges.size() > 0) {
    auto new_vertices = I->V;
    // Average positions of vertices shared by consecutive unique edges
    for (auto i = 0; i < unique_edges.size(); i += 2) {
      for (auto j = 0; j < i; j += 2) {
        auto v0 = -1, v1 = -1, v2 = -1;
        if (unique_edges[i] == unique_edges[j]) {
          v0 = unique_edges[i];
          v1 = unique_edges[i + 1];
          v2 = unique_edges[j + 1];
        } else if (unique_edges[i + 1] == unique_edges[j]) {
          v0 = unique_edges[i + 1];
          v1 = unique_edges[i];
          v2 = unique_edges[j + 1];
        } else if (unique_edges[i] == unique_edges[j + 1]) {
          v0 = unique_edges[i];
          v1 = unique_edges[i + 1];
          v2 = unique_edges[j];
        } else if (unique_edges[i + 1] == unique_edges[j + 1]) {
          v0 = unique_edges[i + 1];
          v1 = unique_edges[i];
          v2 = unique_edges[j];
        }
        if (v0 >= 0) {
          for (auto k = 0; k < 3; k++) {
            new_vertices[3 * v0 + k] =
                (1. / 3.) *
                (I->V[3 * v0 + k] + I->V[3 * v1 + k] + I->V[3 * v2 + k]);
          }
        }
      }
    }
    I->V = std::move(new_vertices);
  }
}
