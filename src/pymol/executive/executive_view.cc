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

// The camera: extents, zooming, orienting, centering and curves.

#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/geometry/matrix.h"
#include "pymol/objects/object_curve.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/geometry/ttt.h"
#include "pymol/executive/executive_internal.h"
pymol::Result<> ExecutiveOrient(PyMOLGlobals* G, const char* sele, int state,
                                float animate, int complete, float buffer,
                                int quiet) {
  double egval[3], egvali[3];
  double evect[3][3];
  float m[4][4], mt[4][4];
  float t[3];
  const float _0 = 0.0F;
  int a, b;

  double mi[16];

  SelectorTmp tmpsele(G, sele);
  sele = tmpsele.getName();

  if (!ExecutiveGetMoment(G, sele, mi, state)) {
    return {};
  }

  if (!MatrixEigensolveC33d(G, mi, egval, egvali, (double*)(void*)evect)) {
    normalize3d(evect[0]);
    normalize3d(evect[1]);
    normalize3d(evect[2]);

    for (a = 0; a < 3; a++) {
      for (b = 0; b < 3; b++) {
        m[a][b] = (float)evect[b][a]; /* fill columns */
      }
    }

    for (a = 0; a < 3; a++) { /* expand to 4x4 */
      m[3][a] = 0;
      m[a][3] = 0;
    }
    m[3][3] = 1.0;

    normalize3f(m[0]); /* cross normalization (probably unnec.)  */
    normalize3f(m[1]);
    normalize3f(m[2]);

    for (a = 0; a < 3; a++) /* convert to row-major */
      for (b = 0; b < 3; b++) mt[a][b] = m[b][a];

    cross_product3f(mt[0], mt[1], t); /* insure right-handed matrix */
    if (dot_product3f(t, mt[2]) < 0.0) {
      mt[2][0] = -mt[2][0];
      mt[2][1] = -mt[2][1];
      mt[2][2] = -mt[2][2];
    }

    for (a = 0; a < 3; a++) /* convert back to column major */
      for (b = 0; b < 3; b++) m[a][b] = mt[b][a];

    if (animate < 0.0F) {
      if (SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }
    if (animate != 0.0F) ScenePrimeAnimation(G);

    {
      float old_mat[16];
      float new_mat[16];
      float x, y, z;
      copy44f(SceneGetMatrix(G), old_mat);

      SceneSetMatrix(G, m[0]); /* load matrix */

      /* there must  be a more elegant to get the PC on X and the SC
       * on Y then what is shown below, but I couldn't get it to work.
       * I tried swapping the eigen-columns around but either that is
       * a bogus approach (?) or my code was buggy.  Hence the following...*/

      if ((egval[0] < egval[2]) && (egval[2] < egval[1])) { /* X < Z < Y */
        SceneRotate(G, 90, 1, 0, 0);                        /*1<-->2 */
      } else if ((egval[1] < egval[0]) &&
                 (egval[0] < egval[2])) { /* Y < X < Z */
        SceneRotate(G, 90, 0, 0, 1);      /*0<-->1 */
      } else if ((egval[1] < egval[2]) &&
                 (egval[2] < egval[0])) { /* Y < Z < X */
        SceneRotate(G, 90, 0, 1, 0);      /*1<-->2 */
        SceneRotate(G, 90, 0, 0, 1);      /*0<-->1 */
      } else if ((egval[2] < egval[1]) &&
                 (egval[1] < egval[0])) { /* Z < Y < X */
        SceneRotate(G, 90, 0, 1, 0);      /*0<-->2 */
      } else if ((egval[2] < egval[0]) &&
                 (egval[0] < egval[1])) { /* Z < X < Y */
        SceneRotate(G, 90, 0, 1, 0);      /*0<-->2 */
        SceneRotate(G, 90, 1, 0, 0);      /*0<-->1 */
      }

      /* now choose orientation that has the least perturbation from the
       * starting matrix */

      copy44f(SceneGetMatrix(G), new_mat);

      x = old_mat[0] * new_mat[0] + old_mat[4] * new_mat[4] +
          old_mat[8] * new_mat[8];
      y = old_mat[1] * new_mat[1] + old_mat[5] * new_mat[5] +
          old_mat[9] * new_mat[9];
      z = old_mat[2] * new_mat[2] + old_mat[6] * new_mat[6] +
          old_mat[10] * new_mat[10];

      if ((x > _0) && (y < _0) && (z < _0)) {
        SceneRotate(G, 180, 1, 0, 0);
      } else if ((x < _0) && (y > _0) && (z < _0)) {
        SceneRotate(G, 180, 0, 1, 0);
      } else if ((x < _0) && (y < _0) && (z > _0)) {
        SceneRotate(G, 180, 0, 0, 1);
      }
    }

    /* X < Y < Z  - do nothing - that's what we want */

    ExecutiveWindowZoom(G, sele, buffer, state, complete, false, quiet);
    if (animate != 0.0F) SceneLoadAnimation(G, animate, 0);
  }
  return {};
}
static int count_objects(PyMOLGlobals* G, int public_only) {
  int count = 0;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      if (!public_only)
        count++;
      else if (rec->obj->Name[0] != '_')
        count++;
    }
  }
  return count;
}
int ExecutiveGetCameraExtent(PyMOLGlobals* G, const char* name, float* mn,
                             float* mx, int transformed, int state) {
  int sele;
  ObjectMoleculeOpRec op;
  int flag = false;

  if ((state == -2) || (state == -3)) /* TO DO: support per-object states */
    state = SceneGetState(G);

  PRINTFD(G, FB_Executive)
  " %s: name %s state %d\n", __func__, name, state ENDFD;

  sele = SelectorIndexByName(G, name);

  if (sele >= 0) {
    ObjectMoleculeOpRecInit(&op);
    if (state < 0) {
      op.code = OMOP_CameraMinMax;
    } else {
      op.code = OMOP_CSetCameraMinMax;
      op.cs1 = state;
    }
    op.v1[0] = FLT_MAX;
    op.v1[1] = FLT_MAX;
    op.v1[2] = FLT_MAX;
    op.v2[0] = -FLT_MAX;
    op.v2[1] = -FLT_MAX;
    op.v2[2] = -FLT_MAX;
    op.i1 = 0;
    op.i2 = transformed;
    op.mat1 = SceneGetMatrix(G);

    ExecutiveObjMolSeleOp(G, sele, &op);

    PRINTFD(G, FB_Executive)
    " %s: minmax over %d vertices\n", __func__, op.i1 ENDFD;
    if (op.i1) flag = true;
  }
  copy3f(op.v1, mn);
  copy3f(op.v2, mx);

  PRINTFD(G, FB_Executive)
  " %s: returning %d\n", __func__, flag ENDFD;

  return (flag);
}
int ExecutiveGetExtent(PyMOLGlobals* G, const char* name, float* mn, float* mx,
                       int transformed, int state, int weighted) {
  int sele;
  ObjectMoleculeOpRec op, op2;
  CExecutive* I = G->Executive;
  pymol::CObject* obj;
  int result = false;
  float f1, f2, fmx;
  int a;

  if (WordMatchExact(G, cKeywordCenter, name, true)) {
    SceneGetCenter(G, mn);
    copy3f(mn, mx);
    return 1;
  }
  if (WordMatchExact(G, cKeywordOrigin, name, true)) {
    SceneOriginGet(G, mn);
    copy3f(mn, mx);
    return 1;
  }

  PRINTFD(G, FB_Executive)
  " %s: name %s state %d\n", __func__, name, state ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

  if ((state == -2) ||
      (state == -3)) { /* we want the currently displayed state */
    state = SceneGetState(G);
    op.include_static_singletons =
        true; /* make sure we get the static singletons too */
    op2.include_static_singletons = true;
  }

  op2.i1 = 0;
  op2.v1[0] = -1.0;
  op2.v1[1] = -1.0;
  op2.v1[2] = -1.0;
  op2.v2[0] = 1.0;
  op2.v2[1] = 1.0;
  op2.v2[2] = 1.0;

  {
    auto matched_recs = ExecutiveGetSpecRecsFromPattern(G, name);
    int have_atoms_flag = false;
    int have_extent_flag = false;

    /* first, compute atomic extents */

    if (weighted) {
      op2.i1 = 0;

      op2.v1[0] = 0.0F;
      op2.v1[1] = 0.0F;
      op2.v1[2] = 0.0F;

      op.i1 = 0;

      op.v1[0] = FLT_MAX;
      op.v1[1] = FLT_MAX;
      op.v1[2] = FLT_MAX;

      op.v2[0] = -FLT_MAX;
      op.v2[1] = -FLT_MAX;
      op.v2[2] = -FLT_MAX;
    }

    /* first, handle molecular objects */

    for (auto& recref : matched_recs) {
      auto* rec = &recref;
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectMolecule &&
              rec->obj->type != cObjectAlignment)
            break;
        case cExecSelection:
        case cExecAll:
          if (rec->type == cExecAll)
            sele = SelectorIndexByName(G, cKeywordAll);
          else
            sele = SelectorIndexByName(G, rec->name);
          if (sele >= 0) {
            if (state < 0) {
              op.code = OMOP_MNMX;
            } else {
              op.code = OMOP_CSetMinMax;
              op.cs1 = state;
            }
            op.i2 = transformed;
            ExecutiveObjMolSeleOp(G, sele, &op);
            if (op.i1) {
              have_atoms_flag = true;
            }
            PRINTFD(G, FB_Executive)
            " %s: minmax over %d vertices\n", __func__, op.i1 ENDFD;
          }

          if (weighted) {
            if (state < 0)
              op2.code = OMOP_SUMC;
            else {
              op2.code = OMOP_CSetSumVertices;
              op2.cs1 = state;
            }
            op2.i2 = transformed;
            ExecutiveObjMolSeleOp(G, sele, &op2);
          }
          break;
      }
    }
    if (have_atoms_flag) have_extent_flag = true;

    /* now handle nonmolecular objects */

    for (auto& recref : matched_recs) {
      auto* rec = &recref;
      switch (rec->type) {
        case cExecAll:
          rec = nullptr;
          while (ListIterate(I->Spec, rec, next)) {
            if (rec->type == cExecObject) {
              obj = rec->obj;
              if (!obj->ExtentFlag) {
                switch (obj->type) {
                  case cObjectMap:
                  case cObjectMesh:
                  case cObjectSurface:
                    if (!rec->obj->ExtentFlag) {
                      /* allow object to update extents, if necessary */
                      rec->obj->update();
                    }
                }
              }
              if (obj->ExtentFlag) switch (obj->type) {
                  case cObjectMolecule:
                    break;
                    /* intentional fall-through */
                  default:
                    if (!have_extent_flag) {
                      copy3f(obj->ExtentMin, op.v1);
                      copy3f(obj->ExtentMax, op.v2);
                      have_extent_flag = true;
                    } else {
                      min3f(obj->ExtentMin, op.v1, op.v1);
                      max3f(obj->ExtentMax, op.v2, op.v2);
                    }
                    break;
                }
            }
          }
          break;
        case cExecObject:
          obj = rec->obj;
          if (!obj->ExtentFlag) {
            switch (obj->type) {
              case cObjectMap:
              case cObjectMesh:
              case cObjectSurface:
                if (!rec->obj->ExtentFlag) {
                  /* allow object to update extents, if necessary */
                  rec->obj->update();
                }
            }
          }
          if (obj->ExtentFlag) switch (obj->type) {
              case cObjectMolecule: /* will have been handled above... */
                break;
              default:
                if (!have_extent_flag) {
                  copy3f(obj->ExtentMin, op.v1);
                  copy3f(obj->ExtentMax, op.v2);
                  have_extent_flag = true;
                } else {
                  min3f(obj->ExtentMin, op.v1, op.v1);
                  max3f(obj->ExtentMax, op.v2, op.v2);
                }
                break;
            }
          break;
      }
    }

    if (have_atoms_flag && weighted) {
      if (op2.i1) {
        op2.v1[0] /= op2.i1; /* compute average */
        op2.v1[1] /= op2.i1;
        op2.v1[2] /= op2.i1;

        for (a = 0; a < 3; a++) { /* this puts origin at the weighted center */
          f1 = op2.v1[a] - op.v1[a];
          f2 = op.v2[a] - op2.v1[a];
          if (f1 > f2)
            fmx = f1;
          else
            fmx = f2;
          op.v1[a] = op2.v1[a] - fmx;
          op.v2[a] = op2.v1[a] + fmx;
        }
      }
    }

    if (have_extent_flag) {
      copy3f(op.v1, mn);
      copy3f(op.v2, mx);
    } else {
      zero3f(mn);
      zero3f(mx);
    }

    result = have_extent_flag;
  }

  PRINTFD(G, FB_Executive)
  " %s: returning %d\n", __func__, result ENDFD;

  return result;
}
static int ExecutiveGetMaxDistance(PyMOLGlobals* G, const char* name,
                                   float* pos, float* dev, int transformed,
                                   int state) {
  int sele;
  ObjectMoleculeOpRec op, op2;
  CExecutive* I = G->Executive;
  pymol::CObject* obj;
  int flag = false;
  float f1, fmx = 0.0F;

  if ((state == -2) || (state == -3)) /* TO DO: support per-object states */
    state = SceneGetState(G);

  PRINTFD(G, FB_Executive)
  " %s: name %s state %d\n", __func__, name, state ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

  {
    auto matched_recs = ExecutiveGetSpecRecsFromPattern(G, name);

    op2.i1 = 0;
    op2.v1[0] = -1.0;
    op2.v1[1] = -1.0;
    op2.v1[2] = -1.0;
    op2.v2[0] = 1.0;
    op2.v2[1] = 1.0;
    op2.v2[2] = 1.0;

    {
      /* first handle molecular objects */

      for (auto& recref : matched_recs) {
        auto* rec = &recref;
        switch (rec->type) {
          case cExecObject:
          case cExecSelection:
          case cExecAll:
            if (rec->type == cExecAll)
              sele = SelectorIndexByName(G, cKeywordAll);
            else
              sele = SelectorIndexByName(G, rec->name);
            if (sele >= 0) {
              if (state < 0) {
                op.code = OMOP_MaxDistToPt;
              } else {
                op.code = OMOP_CSetMaxDistToPt;
                op.cs1 = state;
              }
              op.v1[0] = pos[0];
              op.v1[1] = pos[1];
              op.v1[2] = pos[2];
              op.i1 = 0;
              op.f1 = 0.0F;
              op.i2 = transformed;
              ExecutiveObjMolSeleOp(G, sele, &op);
              fmx = op.f1;

              if (op.i1) flag = true;
            }
            break;
        }
      }
    }

    {
      /* now handle nonmolecular objects */

      for (auto& recref : matched_recs) {
        auto* rec = &recref;
        switch (rec->type) {
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                obj = rec->obj;
                if (obj->ExtentFlag) {
                  switch (obj->type) {
                    case cObjectMolecule:
                      break;
                    default:
                      if (obj->ExtentFlag) {
                        f1 = (float)diff3f(obj->ExtentMin, pos);
                        if (fmx < f1) fmx = f1;
                        f1 = (float)diff3f(obj->ExtentMax, pos);
                        if (fmx < f1) fmx = f1;
                        flag = true;
                        break;
                      }
                  }
                }
              }
            }
            break;
          case cExecObject:
            obj = rec->obj;
            switch (rec->obj->type) {
              case cObjectMolecule:
                break;
              default:
                if (obj->ExtentFlag) {
                  f1 = (float)diff3f(obj->ExtentMin, pos);
                  if (fmx < f1) fmx = f1;
                  f1 = (float)diff3f(obj->ExtentMax, pos);
                  if (fmx < f1) fmx = f1;
                  flag = true;
                }
                break;
            }
        }
      }
    }
  }
  *dev = fmx;
  return (flag);
}
pymol::Result<> ExecutiveWindowZoom(PyMOLGlobals* G, const char* name,
                                    float buffer, int state, int inclusive,
                                    float animate, int quiet) {
  float center[3], radius;
  float mn[3], mx[3], df[3];
  int sele0;

  PRINTFD(G, FB_Executive)
  " ExecutiveWindowZoom-DEBUG: entered\n" ENDFD;
  if (ExecutiveGetExtent(G, name, mn, mx, true, state, true)) {
    if (buffer != 0.0F) {
      mx[0] += buffer;
      mx[1] += buffer;
      mx[2] += buffer;
      mn[0] -= buffer;
      mn[1] -= buffer;
      mn[2] -= buffer;
    }
    subtract3f(mx, mn, df);
    average3f(mn, mx, center);
    if (inclusive) {
      if (!ExecutiveGetMaxDistance(G, name, center, &radius, true, state))
        radius = 0.0;
      radius += buffer;
    } else {
      radius = df[0];
      if (radius < df[1]) radius = df[1];
      if (radius < df[2]) radius = df[2];
      radius = radius / 2.0F;
    }
    if (radius < MAX_VDW) radius = MAX_VDW;
    PRINTFD(G, FB_Executive)
    " %s: zooming with radius %8.3f...state %d\n", __func__, radius,
        state ENDFD;
    PRINTFD(G, FB_Executive)
    " %s: on center %8.3f %8.3f %8.3f...\n", __func__, center[0], center[1],
        center[2] ENDFD;
    if (animate < 0.0F) {
      if (SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }
    if (animate != 0.0F) ScenePrimeAnimation(G);
    SceneOriginSet(G, center, false);
    SceneWindowSphere(G, center, radius);
    if (animate != 0.0F)
      SceneLoadAnimation(G, animate, 0);
    else
      SceneAbortAnimation(G);
    SceneInvalidate(G);
  } else {
    sele0 = SelectorIndexByName(G, name);
    if (sele0 > 0) { /* any valid selection except "all" */
      /* no longer an error to zoom on an empty selection -- just has no effect
       */
      if (!quiet) {
        PRINTFB(G, FB_Executive, FB_Warnings)
        "ExecutiveWindowZoom-Warning: selection doesn't specify any "
        "coordinates.\n" ENDFB(G);
      }
    } else if (ExecutiveValidName(G, name)) {
      PRINTFD(G, FB_Executive)
      " ExecutiveWindowZoom-DEBUG: name valid, but no extents -- using default "
      "view\n" ENDFD;
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      return pymol::make_error(__func__, "selection or object unknown.");
    }
  }
  return {};
}
pymol::Result<> ExecutiveCenter(PyMOLGlobals* G, const char* name, int state,
                                int origin, float animate, float* pos,
                                int quiet) {
  float center[3];
  float mn[3], mx[3];
  int sele0;
  int have_center = false;

  if (name && ExecutiveGetExtent(G, name, mn, mx, true, state, true)) {
    average3f(mn, mx, center);
    have_center = true;
    PRINTFD(G, FB_Executive)
    " %s: centering state %d\n", __func__, state ENDFD;
    PRINTFD(G, FB_Executive)
    " %s: on center %8.3f %8.3f %8.3f...\n", __func__, center[0], center[1],
        center[2] ENDFD;
  } else if (pos) {
    have_center = true;
    copy3f(pos, center);
  }
  if (have_center) {
    if (animate < 0.0F) {
      if (SettingGetGlobal_b(G, cSetting_animation))
        animate = SettingGetGlobal_f(G, cSetting_animation_duration);
      else
        animate = 0.0F;
    }

    if (animate != 0.0F) ScenePrimeAnimation(G);
    if (origin) SceneOriginSet(G, center, false);
    SceneRelocate(G, center);
    SceneInvalidate(G);
    if (animate != 0.0F) SceneLoadAnimation(G, animate, 0);
  } else {
    sele0 = SelectorIndexByName(G, name);
    if (sele0 >= 0) { /* any valid selection except "all" */
      if (!quiet) {
        /* no longer an error to center on an empty selection -- just have no
         * effect */
        PRINTFB(G, FB_Executive, FB_Warnings)
        "ExecutiveCenter-Warning: selection doesn't specify any "
        "coordinates.\n" ENDFB(G);
      }
    } else if (ExecutiveValidName(G, name)) {
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      return pymol::make_error("Selection or object unknown.");
    }
  }
  return {};
}
pymol::Result<> ExecutiveOrigin(PyMOLGlobals* G, const char* sele, int preserve,
                                const char* oname, const float* pos,
                                int state) {
  float center[3];
  float mn[3], mx[3];
  pymol::CObject* obj = nullptr;
  int have_center = false;
  if (oname && oname[0]) {
    obj = ExecutiveFindObjectByName(G, oname);
    if (!obj) return pymol::make_error("Object ", oname, " not found.");
  }
  if (sele && sele[0]) {
    auto s1 = SelectorTmp2::make(G, sele);
    auto has_extent =
        ExecutiveGetExtent(G, s1->getName(), mn, mx, true, state, true);
    if (!has_extent) {
      return pymol::make_error("Could not determine extent of selection.");
    }
    average3f(mn, mx, center);
    have_center = true;
  } else if (pos) {
    copy3f(pos, center);
    have_center = true;
  }

  if (!have_center) {
    return pymol::make_error("Center could not be determined.");
  }

  if (obj) {
    ObjectSetTTTOrigin(obj, center);
    PRINTFB(G, FB_Executive, FB_Blather)
    " %s: origin for %s set to %8.3f %8.3f %8.3f\n", __func__, oname, center[0],
        center[1], center[2] ENDFB(G);
  } else {
    PRINTFB(G, FB_Executive, FB_Blather)
    " %s: scene origin set to %8.3f %8.3f %8.3f\n", __func__, center[0],
        center[1], center[2] ENDFB(G);
    SceneOriginSet(G, center, preserve);
  }
  SceneInvalidate(G);
  return {};
}
int ExecutiveGetMoment(PyMOLGlobals* G, const char* name, double* mi,
                       int state) {
  int sele;
  ObjectMoleculeOpRec op;
  int a, b;
  int c = 0;

  if ((state == -2) || (state == -3)) /* TO DO: support per-object states */
    state = SceneGetState(G);

  sele = SelectorIndexByName(G, name);
  if (sele >= 0) {
    ObjectMoleculeOpRecInit(&op);
    if (state < 0) {
      op.code = OMOP_SUMC;
    } else {
      op.code = OMOP_CSetSumVertices;
      op.cs1 = state;
    }

    op.v1[0] = 0.0;
    op.v1[1] = 0.0;
    op.v1[2] = 0.0;
    op.i1 = 0;
    op.i2 = 0; /* untransformed...is this right? */

    ExecutiveObjMolSeleOp(G, sele, &op);

    if (op.i1) { /* any vertices? */
      c += op.i1;
      scale3f(op.v1, 1.0F / op.i1, op.v1); /* compute raw average */
      if (state < 0) {
        op.code = OMOP_MOME;
      } else {
        op.code = OMOP_CSetMoment;
        op.cs1 = state;
      }
      for (a = 0; a < 3; a++)
        for (b = 0; b < 3; b++) op.d[a][b] = 0.0;
      ExecutiveObjMolSeleOp(G, sele, &op);
      {
        double* p = mi;
        for (a = 0; a < 3; a++)
          for (b = 0; b < 3; b++) *(p++) = op.d[a][b];
      }
    }
  } else {
    identity33d(mi);
  }

  return (c);
}

/**
 * The `zoom` argument takes the following values:
 *
 *   -1 = use `auto_zoom` setting
 *    0 = do nothing
 *    1 = zoom `obj` if `is_new` is ture, otherwise do nothing
 *    2 = zoom `obj`
 *    3 = zoom current state of `obj`
 *    4 = zoom all
 *    5 = zoom `obj` if it is the only object
 *
 */
void ExecutiveDoZoom(PyMOLGlobals* G, pymol::CObject* obj, int is_new, int zoom,
                     int quiet) {
  if (zoom) {
    if (zoom < 0) {
      zoom = SettingGetGlobal_i(G, cSetting_auto_zoom);
      if (zoom < 0) {
        zoom = 1;
      }
    }
    switch (zoom) {
      case 1: /* zoom when new */
        if (is_new) ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);
        /* (all states) */
        break;
      case 2: /* zoom always */
        ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);
        /* (all states) */
        break;
      case 3: /* always zoom current state */
        ExecutiveWindowZoom(G, obj->Name, 0.0,
                            ObjectGetCurrentState(obj, false), 0, 0, quiet);
        break;
      case 4: /* zoom all objects */
        ExecutiveWindowZoom(G, cKeywordAll, 0.0, -1, 0, 0, quiet);
        break;
      case 5: /* zoom first object only */
        if (count_objects(G, true) == 1)
          ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0,
                              quiet); /* (all states) */
        break;
    }
  }
}

pymol::Result<> ExecutiveCurveNew(PyMOLGlobals* G,
                                  pymol::zstring_view curveName,
                                  pymol::zstring_view curvetype) {
  auto oldObj = ExecutiveFindObject<ObjectCurve>(G, curveName);
  if (oldObj) {
    return pymol::make_error("Curve of this name already exists.");
  }
  auto obj = std::make_unique<ObjectCurve>(G);
  obj->setName(curveName);
  ExecutiveManageObject(G, obj.release(), false, true);
  return {};
}

static pymol::Result<> ExecutiveMoveObjectOnCurve(PyMOLGlobals* G,
                                                  pymol::CObject& mobileObj,
                                                  ObjectCurve& curveObj,
                                                  float t) {
  auto pos = curveObj.getPosition(t);
  if (!mobileObj.TTTFlag) {
    initializeTTT44f(mobileObj.TTT);
    mobileObj.TTTFlag = true;
  }
  auto ttt = pymol::TTT::from_pymol_2_legacy(mobileObj.TTT);
  ttt.setTranslation(pos);
  std::copy_n(glm::value_ptr(pymol::TTT::as_pymol_2_legacy(ttt)), 16,
              mobileObj.TTT);
  return {};
}

/**
 * @return the center of an object's extent
 * @param objName the name of object to face
 */
static pymol::Result<glm::vec3> ExecutiveGetObjectExtentCenter(
    PyMOLGlobals* G, pymol::zstring_view name) {
  glm::vec3 min;
  glm::vec3 max;
  bool asTransformed = true;
  bool quiet = true;
  auto ext = ExecutiveGetExtent(G, name.c_str(), glm::value_ptr(min),
                                glm::value_ptr(max), asTransformed,
                                cStateCurrent, quiet);
  if (!ext) {
    return pymol::make_error("Couldn't get extent of: ", name.c_str());
  }
  return (min + max) / 2.0f;
}

/**
 * Moves the camera along a curve
 * @param tarObj the curve to move the camera along
 * @param t the progress along the curve
 */
static pymol::Result<> ExecutiveMoveCameraOnCurve(PyMOLGlobals* G,
                                                  ObjectCurve& curveObj,
                                                  float t) {
  auto& cam = G->Scene->m_view;
  auto pos = curveObj.getPosition(t);
  auto dir = glm::normalize(curveObj.getNormalizedDirection(t));
  auto tfx = glm::lookAt(pos, pos + dir, glm::vec3(0.0f, 1.0f, 0.0f));
  cam.setView(SceneView::FromWorldHomogeneous(tfx, cam.getView()));

  return {};
}

pymol::Result<> ExecutiveMoveOnCurve(PyMOLGlobals* G,
                                     pymol::zstring_view mobileObjectName,
                                     pymol::zstring_view curveObjectName,
                                     float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  auto curveObj = ExecutiveFindObject<ObjectCurve>(G, curveObjectName);
  if (!curveObj) {
    return pymol::make_error("Curve object not found.");
  }
  if (mobileObjectName == "_Camera") {
    return ExecutiveMoveCameraOnCurve(G, *curveObj, t);
  }
  auto mobileObj = ExecutiveFindObjectByName(G, mobileObjectName);
  if (!mobileObj) {
    return pymol::make_error("Mobile object not found.");
  }
  if (curveObj == mobileObj) {
    return pymol::make_error("Mobile and curve cannot be the same.");
  }
  return ExecutiveMoveObjectOnCurve(G, *mobileObj, *curveObj, t);
}

/**
 * Directs the view toward an object
 * @param tarObj the object to face
 */
static pymol::Result<> ExecutiveCameraLookAt(PyMOLGlobals* G,
                                             pymol::CObject& tarObj) {
  auto& cam = G->Scene->m_view;

  auto target = ExecutiveGetObjectExtentCenter(G, tarObj.Name);
  if (!target) {
    return target.error();
  }

  auto tfx = glm::lookAt(cam.worldPos(), *target, glm::vec3(0.0f, 1.0f, 0.0f));
  cam.setView(SceneView::FromWorldHomogeneous(tfx, cam.getView()));
  return {};
}

static pymol::Result<> ExecutiveObjectLookAt(PyMOLGlobals* G,
                                             pymol::CObject& targetObj,
                                             pymol::CObject& mobileObj) {
  return {};
}

pymol::Result<> ExecutiveLookAt(PyMOLGlobals* G,
                                pymol::zstring_view targetObjectName,
                                pymol::zstring_view mobileObjectName) {
  auto targetObj = ExecutiveFindObjectByName(G, targetObjectName);
  if (!targetObj) {
    return pymol::make_error("Target object not found.");
  }
  if (mobileObjectName == "_Camera") {
    return ExecutiveCameraLookAt(G, *targetObj);
  }
  auto mobileObj = ExecutiveFindObjectByName(G, mobileObjectName);
  if (!mobileObj) {
    return pymol::make_error("Mobile object not found.");
  }
  return ExecutiveObjectLookAt(G, *targetObj, *mobileObj);
}
