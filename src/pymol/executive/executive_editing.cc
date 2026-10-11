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

// Editing structures: bonds, hydrogens, atom flags, sculpting and measurements.

#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list.h"
#include "pymol/base/list_macros.h"
#include "pymol/objects/object_dist.h"
#include "pymol/objects/object_molecule3.h"
#include "pymol/objects/object_volume.h"
#include "pymol/reps/rep_dot.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

#define cIndicateSele "indicate"

pymol::Result<> ExecutiveFixChemistry(PyMOLGlobals* G, const char* s1,
                                      const char* s2, int invalidate,
                                      int quiet) {
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;

  assert((sele1 >= 0) && (sele2 >= 0));
  {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject)
        if (rec->obj->type == cObjectMolecule) {
          ObjectMoleculeFixChemistry((ObjectMolecule*)rec->obj, sele1, sele2,
                                     invalidate);
        }
    }
  }
  return {};
}

pymol::Result<> ExecutiveAssignSS(PyMOLGlobals* G, const char* target,
                                  int state, const char* context, int preserve,
                                  ObjectMolecule* single_object, int quiet) {
  if (!target[0]) {
    return pymol::make_error("selection must not be empty");
  }
  pymol::Result<SelectorTmp> targetTmp;
  pymol::Result<SelectorTmp> contextTmp;
  int sele0 = -1;
  int sele1 = -1;
  sele0 = SelectorIndexByName(G, target);
  if (sele0 < 0) {
    targetTmp = SelectorTmp::make(G, target);
    p_return_if_error(targetTmp);
    sele0 = targetTmp->getIndex();
    assert(sele0 >= 0);
  }
  if (!context || !context[0]) {
    sele1 = sele0;
  } else {
    contextTmp = SelectorTmp::make(G, context);
    p_return_if_error(contextTmp);
    sele1 = contextTmp->getIndex();
    assert(sele1 >= 0);
  }
  SelectorAssignSS(G, sele0, sele1, state, preserve, single_object, quiet);
  return {};
}

/**
 * Performs a window average of coordinate states (trajectories).
 *
 * @param pbc Consider periodic boundary conditions
 */
pymol::Result<> ExecutiveSmooth(PyMOLGlobals* G, const char* selection,
                                int cycles, int window, int first, int last,
                                int ends, int quiet, float dist_cutoff,
                                bool pbc) {
  SETUP_SELE(selection, tmpsele1, sele);
  const char* name = tmpsele1->getName();

  ObjectMoleculeOpRec op;
  int range, offset;
  int end_skip = 0;
  bool loop = false;

  PRINTFD(G, FB_Executive)
  " %s: entered %s,%d,%d,%d,%d,%d\n", __func__, name, cycles, first, last,
      window, ends ENDFD;

  /* count the number of states over which to smooth */
  int max_state = ExecutiveCountStates(G, name) - 1;
  if (last < 0) last = max_state;
  if (first < 0) first = 0;
  if (last < first) {
    std::swap(last, first);
  }
  if (last > max_state) last = max_state;

  int const n_state = last - first + 1;

  int const backward = window / 2;
  int const forward = window / 2;

  switch (ends) {
    case 0:
      end_skip = 1;
      break;
    case 1:
      end_skip = 0;
      break;
    case 2:
      end_skip = backward;
      break;
    case 3: /* cyclic averaging */
      end_skip = 0;
      loop = true;
      break;
    default:
      end_skip = 0;
      break;
  }

  if (ends) {
    range = (last - first) + 1;
    offset = 0;
  } else {
    range = (last - end_skip) - (first + end_skip) + 1;
    offset = end_skip;
  }

  PRINTFD(G, FB_Executive)
  " %s: first %d last %d n_state %d backward %d forward %d range %d\n",
      __func__, first, last, n_state, backward, forward, range ENDFD;

  auto const window_abs = std::abs(window);
  if (window_abs < 2) {
    return pymol::make_error("window must be at least size 2");
  }

  if (n_state < window_abs) {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " %s: Window size (%d) larger than number of frames (%d)\n", __func__,
          window_abs, n_state ENDFB(G);
    }
    return {};
  }

  float const dist_cutoff_sq =
      dist_cutoff > 0 ? (dist_cutoff * dist_cutoff) : -1;

  /* determine storage req */
  ObjectMoleculeOpRecInit(&op);
  op.code = OMOP_CountAtoms;
  op.i1 = 0;
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto const n_atom = size_t(op.i1);

  // TODO Could we call SelectorCountAtoms instead?
  assert(n_atom == SelectorCountAtoms(G, sele, cStateAll));

  if (!n_atom) {
    return {};
  }

  auto pbc_obj = pbc ? SelectorGetFirstObjectMolecule(G, sele) : nullptr;
  if (pbc_obj) {
    ObjectMoleculePBCUnwrap(*pbc_obj);
  }

  auto const n_index = n_atom * n_state;
  auto coord0 = std::vector<float>(3 * n_index);
  auto flag0 = std::vector<int>(n_index);
  auto flag1 = std::vector<int>(n_index);

  /* get the data */
  if (!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
    " Smooth: copying coordinates to temporary arrays..\n" ENDFB(G);
  }
  op.code = OMOP_CSetIdxGetAndFlag;
  op.i1 = n_atom;
  op.i2 = 0;
  op.cs1 = first;
  op.cs2 = last;
  op.vv1 = coord0.data();
  op.ii1 = flag0.data();
  op.nvv1 = 0;
  ExecutiveObjMolSeleOp(G, sele, &op);

  PRINTFD(G, FB_Executive)
  " %s: got %d %d\n", __func__, op.i2, op.nvv1 ENDFD;

  // make a copy of the input coordinates
  auto coord1 = coord0;

  for (int a = 0; a < cycles; ++a) {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
      " Smooth: smoothing (pass %d)...\n", a + 1 ENDFB(G);
    }

    // use coordiantes from previous pass as input
    std::swap(coord0, coord1);

    // views for friendly indexing
    auto const* coord0x3 = pymol::reshape<3>(coord0.data());
    auto* const coord1x3 = pymol::reshape<3>(coord1.data());

#pragma omp parallel for
    for (int b = 0; b < range; ++b) {
      int const st_b = b + offset;

      if (st_b < end_skip || st_b >= (n_state - end_skip)) {
        continue;
      }

      for (int c = 0; c < n_atom; ++c) {
        auto const index_c = n_atom * st_b + c;

        if (!flag0[index_c]) {
          continue;
        }

        float sum[3] = {};
        int cnt = 0;
        int st_prev = 0;
        for (int d = -backward; d <= forward; ++d) {
          int st = st_b + d;
          if (loop) {
            if (st < 0) {
              st = n_state + st;
            } else if (st >= n_state) {
              st = st - n_state;
            }
          } else {
            if (st < 0) {
              st = 0;
            } else if (st >= n_state) {
              st = n_state - 1;
            }
          }

          auto const index = n_atom * st + c;
          if (!flag0[index]) {
            continue;
          }

          float const* v0 = coord0x3[index];

          if (dist_cutoff_sq > 0 && cnt) {
            auto const* v0_prev = coord0x3[n_atom * st_prev + c];
            if (diffsq3f(v0, v0_prev) > dist_cutoff_sq) {
              if (d <= 0) {
                scale3f(v0, cnt, sum);
              } else {
                for (; d <= forward; ++d) {
                  add3f(sum, v0_prev, sum);
                  ++cnt;
                }
                break;
              }
            }
          }

          /* atom's avg position */
          add3f(sum, v0, sum);
          ++cnt;
          st_prev = st;
        }
        if (cnt) {
          flag1[index_c] = true;
          scale3f(sum, 1.0F / cnt, coord1x3[index_c]);
        }
      }
    }
  }

  if (!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
    " Smooth: updating coordinates...\n" ENDFB(G);
  }

  /* set the new coordinates */

  op.code = OMOP_CSetIdxSetFlagged;
  op.i1 = n_atom;
  op.i2 = 0;
  if (ends) {
    op.cs1 = first;
    op.cs2 = last;
    op.vv1 = coord1.data();
    op.ii1 = flag1.data();
  } else {
    op.cs1 = first + end_skip;
    op.cs2 = last - end_skip;
    op.vv1 = coord1.data() + (end_skip * n_atom * 3);
    op.ii1 = flag1.data() + (end_skip * n_atom);
  }
  op.nvv1 = 0;

  ExecutiveObjMolSeleOp(G, sele, &op);
  PRINTFD(G, FB_Executive)
  " %s: put %d %d\n", __func__, op.i2, op.nvv1 ENDFD;

  if (pbc_obj) {
    ObjectMoleculePBCWrap(*pbc_obj);
  }

  return {};
}
int ExecutiveSculptIterateAll(PyMOLGlobals* G) {
  int active = false;
  float center_array[8] = {0.0F, 0.0F, 0.0F, 0.0F};
  float* center = center_array;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* objMol;
  CGOReset(G->DebugCGO);

  if (SettingGetGlobal_b(G, cSetting_sculpting)) {
    if (!SettingGetGlobal_b(G, cSetting_sculpt_auto_center)) center = nullptr;

    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule*)rec->obj;
          if (SettingGet_b(G, nullptr, objMol->Setting.get(),
                           cSetting_sculpting)) {
            constexpr int state = -2;  // current state
            ObjectMoleculeSculptIterate(
                objMol, state,
                SettingGet_i(G, nullptr, objMol->Setting.get(),
                             cSetting_sculpting_cycles),
                center);
            active = true;
          }
        }
      }
    }
    if (center && (center[3] > 1.0F)) {
      float pos[3];
      SceneGetCenter(G, pos);
      center[3] = 1.0F / center[3];
      scale3f(center, center[3], center);
      center[7] = 1.0F / center[7];
      scale3f(center + 4, center[7], center + 4);
      subtract3f(center, center + 4, center);
      add3f(pos, center, center);
      ExecutiveCenter(G, nullptr, -1, true, false, center, true);
    }
  }
  if (active) {
    EditorInvalidateShaderCGO(G);
  }
  return (active);
}
float ExecutiveSculptIterate(PyMOLGlobals* G, const char* name, int state,
                             int n_cycle) {
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* objMol;
  float total_strain = 0.0F;

  if (WordMatchExact(G, name, cKeywordAll, true)) {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule*)rec->obj;
          total_strain +=
              ObjectMoleculeSculptIterate(objMol, state, n_cycle, nullptr);
        }
      }
    }
  } else if (!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s not found.\n", name ENDFB(G);
  } else if (obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
  } else {
    total_strain = ObjectMoleculeSculptIterate((ObjectMolecule*)obj, state,
                                               n_cycle, nullptr);
  }
  return (total_strain);
}
int ExecutiveSculptActivate(PyMOLGlobals* G, const char* name, int state,
                            int match_state, int match_by_segment) {
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  SpecRec* rec = nullptr;
  ObjectMolecule* objMol;
  CExecutive* I = G->Executive;
  int ok = true;
  if (state < 0) state = SceneGetState(G);

  if (WordMatchExact(G, name, cKeywordAll, true)) {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptImprint(objMol, state, match_state,
                                      match_by_segment);
        }
      }
    }
  } else if (!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else if (obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectMoleculeSculptImprint((ObjectMolecule*)obj, state, match_state,
                                match_by_segment);
  }
  return (ok);
}
int ExecutiveSculptDeactivate(PyMOLGlobals* G, const char* name) {
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  SpecRec* rec = nullptr;
  ObjectMolecule* objMol;
  CExecutive* I = G->Executive;

  int ok = true;

  if (WordMatchExact(G, name, cKeywordAll, true)) {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          objMol = (ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptClear(objMol);
        }
      }
    }
  } else if (!obj) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else if (obj->type != cObjectMolecule) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Executive-Error: object %s is not a molecular object.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectMoleculeSculptClear((ObjectMolecule*)obj);
  }
  return (ok);
}
pymol::Result<> ExecutiveSetGeometry(PyMOLGlobals* G, const char* s1, int geom,
                                     int valence) {
  SETUP_SELE_DEFAULT(1);

  int count = 0;
  auto recs = pymol::make_list_adapter(G->Executive->Spec);
  for (auto& rec : recs) {
    auto obj = ExecutiveIsObjectType(rec, cObjectMolecule) ? rec.obj : nullptr;
    if (obj) {
      count += ObjectMoleculeSetGeometry(G, static_cast<ObjectMolecule*>(obj),
                                         sele1, geom, valence);
    }
  }
  if (count == 0) {
    return pymol::make_error("Empty selection.");
  }
  return {};
}

int ExecutivePhiPsi(PyMOLGlobals* G, const char* s1, ObjectMolecule*** objVLA,
                    int** iVLA, float** phiVLA, float** psiVLA, int state) {
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  int result = false;
  ObjectMoleculeOpRec op1;
  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;
    op1.i2 = state;
    op1.obj1VLA = VLAlloc(ObjectMolecule*, 1000);
    op1.i1VLA = VLAlloc(int, 1000);
    op1.f1VLA = VLAlloc(float, 1000);
    op1.f2VLA = VLAlloc(float, 1000);
    op1.code = OMOP_PhiPsi;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    result = op1.i1;
    VLASize(op1.i1VLA, int, op1.i1);
    VLASize(op1.obj1VLA, ObjectMolecule*, op1.i1);
    VLASize(op1.f1VLA, float, op1.i1);
    VLASize(op1.f2VLA, float, op1.i1);
    *iVLA = op1.i1VLA;
    *objVLA = op1.obj1VLA;
    *phiVLA = op1.f1VLA;
    *psiVLA = op1.f2VLA;
  } else {
    *objVLA = nullptr;
    *iVLA = nullptr;
    *phiVLA = nullptr;
    *psiVLA = nullptr;
  }
  return (result);
}
/**
 * Macro with error handling to get a single atom vertex and associated temp
 * selection. Assigns `vertexN`, `tmpseleN`, and `resultN`
 */
#define SETUP_SINGLE_ATOM_VERTEX(N)                                    \
  auto tmpsele##N = SelectorTmp::make(G, s##N);                        \
  p_return_if_error_prefixed(tmpsele##N, "Selection " #N ": ");        \
  auto result##N =                                                     \
      SelectorGetSingleAtomVertex(G, (tmpsele##N)->getIndex(), state); \
  p_return_if_error_prefixed(result##N, "Selection " #N ": ");         \
  const float* vertex##N = (result##N).result().data()

/**
 * Get the distance between two atoms.
 * Implementation of `cmd.get_distance()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param state object state
 * @return distance in Angstrom
 */
pymol::Result<float> ExecutiveGetDistance(PyMOLGlobals* G, const char* s1,
                                          const char* s2, int state) {
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  return static_cast<float>(diff3f(vertex1, vertex2));
}
/**
 * Get the angle between tree atoms s0-s1-s2.
 * Implementation of `cmd.get_angle()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param s2 single-atom selection expression
 * @param state object state
 * @return angle in degrees
 */
pymol::Result<float> ExecutiveGetAngle(PyMOLGlobals* G, const char* s1,
                                       const char* s2, const char* s3,
                                       int state) {
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);

  float d1[3], d2[3];
  subtract3f(vertex1, vertex2, d1);
  subtract3f(vertex3, vertex2, d2);
  return rad_to_deg(get_angle3f(d1, d2));
}
/**
 * Get the dihedral angle between four atoms.
 * Implementation of `cmd.get_dihedral()`
 *
 * @param s0 single-atom selection expression
 * @param s1 single-atom selection expression
 * @param s2 single-atom selection expression
 * @param s3 single-atom selection expression
 * @param state object state
 * @return dihedral angle in degrees
 */
pymol::Result<float> ExecutiveGetDihe(PyMOLGlobals* G, const char* s1,
                                      const char* s2, const char* s3,
                                      const char* s4, int state) {
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);
  SETUP_SINGLE_ATOM_VERTEX(4);
  return rad_to_deg(get_dihedral3f(vertex1, vertex2, vertex3, vertex4));
}
pymol::Result<> ExecutiveSetDihe(PyMOLGlobals* G, const char* s1,
                                 const char* s2, const char* s3, const char* s4,
                                 float value, int state, int quiet) {
  SETUP_SINGLE_ATOM_VERTEX(1);
  SETUP_SINGLE_ATOM_VERTEX(2);
  SETUP_SINGLE_ATOM_VERTEX(3);
  SETUP_SINGLE_ATOM_VERTEX(4);
  float current =
      rad_to_deg(get_dihedral3f(vertex1, vertex2, vertex3, vertex4));
  float change = value - current;
  auto save_state = SceneGetState(G);
  SceneSetFrame(G, -1, state); /* KLUDGE ALERT!
                                * necessary because the editor
                                * can only work on the current state...this
                                * needs to be changed.*/
  EditorSelect(G, tmpsele3->getName(), tmpsele2->getName(), "", "", false, true,
               true);
  EditorTorsion(G, change);
  SceneSetFrame(G, -1, save_state);
  if (!quiet) {
    PRINTFB(G, FB_Editor, FB_Actions)
    " SetDihedral: adjusted to %5.3f\n", value ENDFB(G);
  }
  return {};
}
pymol::Result<float> ExecutiveGetArea(PyMOLGlobals* G, const char* sele,
                                      int state, bool load_b) {
  SETUP_SELE(sele, tmpsele0, sele0);

  auto obj0 = SelectorGetSingleObjectMolecule(G, sele0);
  if (!obj0) {
    if (SelectorCountAtoms(G, sele0, state) > 0)
      return pymol::make_error("Selection must be within a single object");
    return 0.f;
  }

  auto cs = obj0->getCoordSet(state);
  if (!cs) return pymol::make_error("Invalid state");

  auto rep = (RepDot*)RepDotDoNew(cs, cRepDotAreaType, state);
  if (!rep) return pymol::make_error("Can't get dot representation.");

  if (load_b) {
    /* zero out B-values within selection */
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_SetB;
    op.f1 = 0.0;
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele0, &op);
  }

  float const* const area = rep->A;
  int const* const ati = rep->Atom;

  AtomInfoType* ai = nullptr;
  int known_member = -1;
  bool is_member = false;
  float result = 0.f;

  for (int a = 0; a < rep->N; ++a) {
    if (known_member != ati[a]) {
      known_member = ati[a];
      ai = obj0->AtomInfo + known_member;
      is_member = SelectorIsMember(G, ai->selEntry, sele0);
    }

    if (is_member) {
      result += area[a];
      if (load_b) ai->b += area[a];
    }
  }

  delete rep;
  return result;
}
pymol::Result<> ExecutiveRenameObjectAtoms(PyMOLGlobals* G,
                                           const char* selection, int force,
                                           int quiet) {
  SETUP_SELE(selection, tmpsele1, sele);

  {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_RenameAtoms;
    op.i1 = 0;
    op.i2 = force;
    ExecutiveObjMolSeleOp(G, sele, &op);

    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Actions)
      " Rename: renamed %d atoms.\n", op.i1 ENDFB(G);
    }
  }
  return {};
}
pymol::Result<> ExecutiveFuse(PyMOLGlobals* G, const char* s0, const char* s1,
                              int mode, int recolor, int move_flag) {
  int i0 = -1;
  int i1 = -1;
  SelectorID_t sele2 = cSelectionInvalid;
  ObjectMolecule *obj0, *obj1;
  ObjectMoleculeOpRec op;
#define tmp_fuse_sele "tmp_fuse_sele"

  SETUP_SELE_DEFAULT_PREFIXED(0, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);

  {
    {
      EditorInactivate(G);
      obj0 = SelectorGetSingleObjectMolecule(G, sele0);
      obj1 = SelectorGetSingleObjectMolecule(G, sele1);
      if (obj0) i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);
      if (obj1) i1 = ObjectMoleculeGetAtomIndex(obj1, sele1);
      if (obj0 && obj1 && (i0 >= 0) && (i1 >= 0) && (obj0 != obj1)) {
        ObjectMoleculeVerifyChemistry(obj0, -1);
        ObjectMoleculeVerifyChemistry(obj1, -1);

        SelectorCreate(G, tmp_fuse_sele, nullptr, obj0, 1, nullptr);
        sele2 = SelectorIndexByName(G, tmp_fuse_sele);
        if (mode) {
          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_PrepareFromTemplate;
          op.ai = obj1->AtomInfo + i1;
          op.i1 = mode;
          op.i2 = 0;
          op.i3 = recolor;
          if (recolor) op.i4 = obj1->Color;
          ExecutiveObjMolSeleOp(G, sele2, &op);
        }
        SelectorDelete(G, tmp_fuse_sele);

        p_return_if_error(
            ObjectMoleculeFuse(obj1, i1, obj0, i0, mode != 3, move_flag));
      }
    }
  }
  return {};
}
pymol::Result<> ExecutiveSpheroid(PyMOLGlobals* G, const char* name,
                                  int average) { /* EXPERIMENTAL */
  CExecutive* I = G->Executive;
  pymol::CObject* os = nullptr;

  if (strlen(name)) {
    os = ExecutiveFindObjectByName(G, name);
    if (!os)
      return pymol::make_error("Object not found.");
    else if (os->type != cObjectMolecule) {
      return pymol::make_error("Bad object type.");
    }
  }

  if (os || (!strlen(name))) { /* sort one or all */
    SpecRec* rec = nullptr;
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject)
        if (rec->obj->type == cObjectMolecule)
          if ((!os) || (rec->obj == os)) {
            auto obj = (ObjectMolecule*)rec->obj;
            ObjectMoleculeCreateSpheroid(obj, average);
            obj->invalidate(cRepAll, cRepInvRep, -1);
          }
    }
    SceneChanged(G);
  }
  return {};
}
pymol::Result<> ExecutiveSort(PyMOLGlobals* G, const char* name) {
  CExecutive* I = G->Executive;
  ObjectMolecule* obj;
  SpecRec* rec = nullptr;
  ObjectMoleculeOpRec op;
  int sele;
  int ok = true;
  if ((!name) || (!name[0])) name = cKeywordAll;

  {
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int changed = false;
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if ((rec->type == cExecObject) &&
                  (rec->obj->type == cObjectMolecule)) {
                obj = (ObjectMolecule*)rec->obj;
                if (ok) ok &= ObjectMoleculeSort(obj);
                if (ok) {
                  changed = true;
                  sele = SelectorIndexByName(G, rec->name);
                  if (sele >= 0) {
                    ObjectMoleculeOpRecInit(&op);
                    op.code = OMOP_INVA;
                    op.i1 = cRepCartoonBit | cRepRibbonBit;
                    op.i2 = cRepInvRep;
                    ExecutiveObjMolSeleOp(G, sele, &op);
                  }
                }
              }
            }
            break;
          case cExecSelection:
            sele = SelectorIndexByName(G, rec->name);
            if (sele >= 0) {
              op.code = OMOP_Sort;
              ExecutiveObjMolSeleOp(G, sele, &op);
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_INVA;
              op.i1 = cRepCartoonBit | cRepRibbonBit;
              op.i2 = cRepInvRep;
              ExecutiveObjMolSeleOp(G, sele, &op);
              ObjectMoleculeOpRecInit(&op);
            }
            break;
          case cExecObject:
            if (rec->obj->type == cObjectMolecule) {
              obj = (ObjectMolecule*)rec->obj;
              if (ok) ok &= ObjectMoleculeSort(obj);
              changed = true;
              sele = SelectorIndexByName(G, rec->name);
              if (sele >= 0) {
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_INVA;
                op.i1 = cRepCartoonBit | cRepRibbonBit;
                op.i2 = cRepInvRep;
                ExecutiveObjMolSeleOp(G, sele, &op);
              }
            }
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
    if (changed) SceneChanged(G);
  }
  return {};
}
pymol::Result<> ExecutiveRemoveAtoms(PyMOLGlobals* G, const char* s1,
                                     int quiet) {
  SETUP_SELE_DEFAULT(1);

  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  ObjectMoleculeOpRec op;

  {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_Remove;
          op.i1 = 0;
          obj = (ObjectMolecule*)rec->obj;
          ObjectMoleculeVerifyChemistry(obj,
                                        -1); /* remember chemistry for later */
          ObjectMoleculeSeleOp(obj, sele1, &op);
          if (op.i1) {
            if (!quiet) {
              PRINTFD(G, FB_Editor)
              " ExecutiveRemove-Debug: purging %i of %i atoms in %s\n", op.i1,
                  obj->NAtom, obj->Name ENDFD;
            }
            ObjectMoleculePurge(obj);
            if (!quiet) {
              PRINTFB(G, FB_Editor, FB_Actions)
              " Remove: eliminated %d atoms in model \"%s\".\n", op.i1,
                  obj->Name ENDFB(G);
            }
          }
        }
      }
    }
  }
  EditorRemoveStale(G);
  return {};
}
pymol::Result<> ExecutiveAddHydrogens(PyMOLGlobals* G, const char* s1,
                                      int quiet, int state, bool legacy) {
  ObjectMoleculeOpRec op;

  if (legacy) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    " %s-Warning: legacy mode was removed\n", __FUNCTION__ ENDFB(G);
  }

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_AddHydrogens;
    op.i1 = state;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
  return {};
}
void ExecutiveFixHydrogens(PyMOLGlobals* G, const char* s1, int quiet) {
  int sele1;
  ObjectMoleculeOpRec op;

  sele1 = SelectorIndexByName(G, s1);
  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_FixHydrogens;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
}
pymol::Result<> ExecutiveFlag(PyMOLGlobals* G, int flag, const char* sele,
                              int action, int quiet) {
  enum {
    cFlagActionReset = 0,
    cFlagActionSet = 1,
    cFlagActionClear = 2,
  };

  if (flag < 0 || flag > 31) {
    return pymol::make_error("flag ", flag, " out of range [0, 31]");
  }

  ObjectMoleculeOpRec op;
  ObjectMoleculeOpRecInit(&op);
  switch (action) {
    case cFlagActionReset:
      op.code = OMOP_Flag;
      break;
    case cFlagActionSet:
      op.code = OMOP_FlagSet;
      break;
    case cFlagActionClear:
      op.code = OMOP_FlagClear;
      break;
    default:
      return pymol::make_error("invalid action ", action);
  }

  op.i1 = 1u << flag;
  op.i2 = ~op.i1;
  op.i3 = 0;
  op.i4 = 0;

  if ((op.i1 & cAtomFlag_exfoliate) && action != cFlagActionClear) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    "The 'exfoliate' flag is deprecated. "
    "Use 'hide surface, (%s)' instead.\n",
        sele ENDFB(G);
  }

  SETUP_SELE(sele, s1, sele1);

  ExecutiveObjMolSeleOp(G, sele1, &op);
  if (Feedback(G, FB_Executive, FB_Actions)) {
    if (!quiet) {
      switch (action) {
        case cFlagActionReset:
          if (op.i3) {
            PRINTF " Flag: flag %d is set in %d of %d atoms.\n", flag, op.i3,
                op.i4 ENDF(G);
          } else {
            PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF(G);
          }
          break;
        case cFlagActionSet:
          PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
        case cFlagActionClear:
          PRINTF " Flag: flag %d cleared on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
      }
    }
  }
  if (SettingGetGlobal_b(G, cSetting_auto_indicate_flags)) {
    auto buffer = pymol::string_format("(flag %d)", flag);
    SelectorCreate(G, cIndicateSele, buffer.c_str(), nullptr, true, nullptr);
    ExecutiveSetObjVisib(G, cIndicateSele, true, false);
    SceneInvalidate(G);
  }
  return {};
}
pymol::Result<> ExecutiveProtect(PyMOLGlobals* G, const char* s1, int mode,
                                 int quiet) {
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Protect;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    if (!quiet) {
      if (Feedback(G, FB_Executive, FB_Actions)) {
        if (op.i2) {
          if (mode) {
            PRINTF " Protect: %d atoms protected from movement.\n",
                op.i2 ENDF(G);
          } else {
            PRINTF " Protect: %d atoms deprotected.\n", op.i2 ENDF(G);
          }
        }
      }
    }
  }
  return {};
}
pymol::Result<> ExecutiveMask(PyMOLGlobals* G, const char* s1, int mode,
                              int quiet) {
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Mask;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    if (!quiet) {
      if (Feedback(G, FB_Executive, FB_Actions)) {
        if (op.i2) {
          if (mode) {
            PRINTF " Mask: %d atoms masked (cannot be picked or selected).\n",
                op.i2 ENDF(G);
          } else {
            PRINTF " Mask: %d atoms unmasked.\n", op.i2 ENDF(G);
          }
        }
      }
    }
    op.code = OMOP_INVA; /* need to invalidate all pickable representations */
    op.i1 = cRepsAtomMask;
    op.i2 = cRepInvPick;
    ExecutiveObjMolSeleOp(G, sele1, &op);
  }
  return {};
}
/**
 * flag > 0:  Set stereo_mode and turn stereo on
 * flag = 0:  Turn off stereo
 * flag = -1: Swap eyes (stereo_shift *= -1)
 * flag = -2: Turn on stereo with current stereo_mode
 * flag = -3: Turn on chromadepth and turn off stereo
 */
pymol::Result<> ExecutiveStereo(PyMOLGlobals* G, int flag) {
  switch (flag) {
    case -3:
      SettingSet(G, cSetting_chromadepth, 1);
      SceneSetStereo(G, 0);
      break;
    case -1:
      SettingSetGlobal_f(G, cSetting_stereo_shift,
                         -SettingGetGlobal_f(G, cSetting_stereo_shift));
      break;
    default:
      SettingSet(G, cSetting_chromadepth, 0);

      if (flag == cStereo_quadbuffer && !G->StereoCapable) {
        return pymol::make_error(
            "no 'quadbuffer' support detected (force with 'pymol -S')");
      }

      if (flag == cStereo_openvr) {
        return pymol::make_error(
            "'openvr' stereo mode not available in this build");
      }

      if (flag > 0) {
        SettingSet(G, cSetting_stereo_mode, flag);
      }

      SceneSetStereo(G, flag != 0);
  }

  SceneDirty(G);
  return {};
}
pymol::Result<> ExecutiveRevalence(PyMOLGlobals* G, const char* s1,
                                   const char* s2, const char* src,
                                   int target_state, int source_state,
                                   int reset, int quiet) {
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  {
    if (src && src[0]) {
      SETUP_SELE(src, tmpsele3, sele3);

      {
        ObjectMolecule* obj3 = SelectorGetSingleObjectMolecule(G, sele3);
        if (!obj3) {
          return pymol::make_error(
              "Revalence can only source a single object at a time.");
        } else {
          ObjectMoleculeOpRec op;

          ObjectMoleculeOpRecInit(&op);
          op.code = OMOP_RevalenceFromSource;
          op.i1 = sele1;
          op.i2 = sele2;
          op.i3 = target_state;
          op.obj3 = obj3;
          op.i4 = sele3;
          op.i5 = source_state;
          op.i6 = quiet;

          ExecutiveObjMolSeleOp(G, sele1, &op);

          /*
             if(ObjectMoleculeXferValences(obj1,sele1,sele2,target_state,obj3,sele3,source_state,quiet))
             { ObjectMoleculeVerifyChemistry(obj1,target_state);
             ObjectMoleculeInvalidate(obj1,cRepAll,cRepInvBonds,target_state);
             }
           */
        }
      }
    } else { /* guess valences */
      ObjectMoleculeOpRec op;

      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_RevalenceByGuessing;
      op.i1 = sele1;
      op.i2 = sele2;
      op.i3 = target_state;
      op.i4 = reset;
      op.i6 = quiet;

      ExecutiveObjMolSeleOp(G, sele1, &op);
    }
  }
  return {};
}
pymol::Result<> ExecutiveBond(PyMOLGlobals* G, const char* s1, const char* s2,
                              int order, int mode, int quiet,
                              pymol::zstring_view symop) {
  int cnt;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  int flag = false;

  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, cSelectionInvalid);

  {
    ObjectMolecule* obj1 = SelectorGetSingleObjectMolecule(G, sele1);
    ObjectMolecule* obj2 = SelectorGetSingleObjectMolecule(G, sele2);
    if ((!obj1) || (!obj2) || (obj1 != obj2)) {
      if ((!quiet) && (mode == 1)) {
        PRINTFB(G, FB_Editor, FB_Warnings)
        "Editor-Warning: bonds cannot be created between objects, only "
        "within.\n" ENDFB(G);
      }
    }
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          switch (mode) {
            case 1: /* add */
              cnt = ObjectMoleculeAddBond((ObjectMolecule*)rec->obj, sele1,
                                          sele2, order, symop);
              if (cnt) {
                if (!quiet) {
                  PRINTFB(G, FB_Editor, FB_Actions)
                  " Bond: %d bonds added to model \"%s\".\n", cnt,
                      rec->obj->Name ENDFB(G);
                  flag = true;
                }
              }
              break;
            case 2: /* adjust */
              cnt = ObjectMoleculeAdjustBonds((ObjectMolecule*)rec->obj, sele1,
                                              sele2, 1, order, symop);
              if (cnt) {
                if (!quiet) {
                  PRINTFB(G, FB_Editor, FB_Actions)
                  " Valence: %d bond valences adjusted in model \"%s\".\n", cnt,
                      rec->obj->Name ENDFB(G);
                  flag = true;
                }
              }
              break;
            case 0: /* remove */
            default:
              cnt = ObjectMoleculeRemoveBonds((ObjectMolecule*)rec->obj, sele1,
                                              sele2);
              if (cnt) {
                if (!quiet) {
                  PRINTFB(G, FB_Editor, FB_Actions)
                  " Unbond: %d bonds removed from model \"%s\".\n", cnt,
                      rec->obj->Name ENDFB(G);
                }
                flag = true;
              }
          }
        }
      }
    }
    if (!flag) {
      if (!quiet) {
        switch (mode) {
          case 1:
            PRINTFB(G, FB_Editor, FB_Warnings)
            "Bond-Warning: no bonds added." ENDFB(G);
            break;
          case 2:
            PRINTFB(G, FB_Editor, FB_Warnings)
            "Valence-Warning: no bond valences changed." ENDFB(G);
            break;
          case 0:
          default:
            PRINTFB(G, FB_Editor, FB_Warnings)
            "Unbond-Warning: no bonds removed." ENDFB(G);
            break;
        }
      }
    }
  }
  return {};
}

/**
 * brief Returns an error if quiet flag is not set
 * @param quiet quiet flag
 * @param msgs error messages
 * @return pymol::make_error with conditionally filled message
 */

template <typename... Ts>
static pymol::Error ErrorMsgIfQuiet(bool quiet, Ts&&... msgs) {
  if (!quiet) {
    return pymol::make_error(std::forward<Ts>(msgs)...);
  } else {
    return pymol::make_error("");
  }
}
pymol::Result<float> ExecutiveAngle(PyMOLGlobals* G, const char* nam,
                                    const char* s1, const char* s2,
                                    const char* s3, int mode, int labels,
                                    int reset, int zoom, int quiet, int state,
                                    int state1, int state2, int state3) {
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);
  SETUP_SELE_DEFAULT_PREFIXED(3, sele2);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  obj =
      ObjectDistNewFromAngleSele(G, obj, sele1, sele2, sele3, mode, labels,
                                 &result, reset, state, state1, state2, state3);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
    ObjectSetName(obj, nam);
    ExecutiveManageObject(G, obj, zoom, quiet);
    if (!labels) ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return rad_to_deg(result);
}
pymol::Result<float> ExecutiveDihedral(PyMOLGlobals* G, const char* nam,
                                       const char* s1, const char* s2,
                                       const char* s3, const char* s4, int mode,
                                       int labels, int reset, int zoom,
                                       int quiet, int state) {
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);
  SETUP_SELE_DEFAULT_PREFIXED(3, sele2);
  SETUP_SELE_DEFAULT_PREFIXED(4, sele3);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  obj = ObjectDistNewFromDihedralSele(G, obj, sele1, sele2, sele3, sele4, mode,
                                      labels, &result, reset, state);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
    ObjectSetName(obj, nam);
    ExecutiveManageObject(G, obj, zoom, quiet);
    if (!labels) ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return rad_to_deg(result);
}

/**
 * Create a distance measurement object
 *
 * @param nam: name of measurement object to create or add to
 * @param s1: selection expression
 * @param s2: selection expression or "same" keyword (shortcut for s1 = s2)
 * @param mode: 0 (any), 1 (bonds), 2 (hbonds), 3 (distance_exclusion),
 * 4 (centroids), 5 (pi-*), 6 (pi-pi), 7 (pi-cat), 8 (vdw-dist-ratio)
 * @param cutoff: Distance cutoff in Angstrom, or vdw-distance ratio cutoff
 * (mode 8)
 * @param labels: Boolean flag to show distance labels
 * @param reset: Boolean flag to clear a distance object if it already exists
 * @return average of measured distances in Angstrom
 */
pymol::Result<float> ExecutiveDistance(PyMOLGlobals* G, const char* nam,
                                       const char* s1, const char* s2, int mode,
                                       float cutoff, int labels, int quiet,
                                       int reset, int state, int zoom,
                                       int state1, int state2) {
  if (strcmp(s1, s2) == 0) {
    s2 = cKeywordSame;
  }

  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  auto obj = ExecutiveFindOrDeleteObject<ObjectDist>(G, nam);
  auto need_manage = !obj;
  float result = -1.0F;

  /* create a new distance from the two selections */
  obj = ObjectDistNewFromSele(G, obj, sele1, sele2, mode, cutoff, labels, reset,
                              &result, state, state1, state2);

  // ObjectDistNewFrom... always succeeds
  assert(obj);

  if (need_manage) {
    switch (mode) {
      case 6:  // pi-pi
        SettingSet(cSetting_dash_color, "0x06c5ff" /* light blue */, obj);
        break;
      case 7:  // pi-cat
        SettingSet(cSetting_dash_color, "0x468c00" /* dark green */, obj);
        break;
      case 8:  // "clashes"
        SettingSet(cSetting_dash_color, "0xff8800" /* light red */, obj);
        break;
      case 9:  // "halogen-bonds"
        SettingSet(cSetting_dash_color, "0xaa00ff" /* dark magenta */, obj);
        break;
      case 10:  // "salted-bridge"
        SettingSet(cSetting_dash_color, "0xff2eff" /* light magenta */, obj);
        break;
    }

    ObjectSetName(obj, nam);
    ExecutiveManageObject(G, obj, zoom, quiet);
    if (!labels) ExecutiveSetRepVisib(G, nam, cRepLabel, 0);
  }

  return result;
}

int ExecutiveAssignAtomTypes(PyMOLGlobals* G, const char* s1, int format,
                             int state, int quiet) {
  int result = 0;
  int sele1;

  sele1 = SelectorIndexByName(G, s1);
  if (state < 0) state = 0;
  {
    if (sele1 >= 0) {
      result = SelectorAssignAtomTypes(G, sele1, state, quiet, format);
    }
  }
  return (result);
}

/**
 * Discard all bonds and do distance based bonding.
 * Implementation of `cmd.rebond()`
 *
 * @param oname object name
 * @param state object state, negative values fall back to current state
 */
pymol::Result<> ExecutiveRebond(PyMOLGlobals* G, const char* oname, int state,
                                bool pbc) {
  auto obj = ExecutiveFindObjectMoleculeByName(G, oname);
  if (!obj) {
    return pymol::make_error("cannot find object");
  }

  auto cs = obj->getCoordSet(state);
  if (!cs) {
    return pymol::make_error("no such state");
  }

  ObjectMoleculeRemoveBonds(obj, 0, 0);

  // Cases where we want to discretely rebond that isn't a pbc case?
  if (obj->DiscreteFlag && pbc) {
    ObjectMoleculeConnectDiscrete(obj, true, 3, pbc);
  } else {
    ObjectMoleculeConnect(obj, cs, true, 3, pbc);
  }

  obj->invalidate(cRepAll, cRepInvAll, -1);

  return {};
}

pymol::Result<> ExecutiveAddBondByIndices(PyMOLGlobals* G,
                                          pymol::zstring_view oname,
                                          unsigned int atm1, unsigned int atm2,
                                          int order) {
  auto obj = ExecutiveFindObject<ObjectMolecule>(G, oname.c_str());
  if (!obj) {
    return pymol::make_error("Cannot find object ", oname);
  }
  return ObjectMoleculeAddBondByIndices(obj, atm1, atm2, order);
}
