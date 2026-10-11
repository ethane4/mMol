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

// Structure alignment and fitting: align, fit, RMS, CE-align and US-align.

#include "pymol/select/atom_iterators.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/geometry/matrix.h"
#include "pymol/objects/object_alignment.h"
#include "pymol/objects/object_molecule3.h"
#include "pymol/ui/seeker.h"
#include "pymol/geometry/ttt.h"
#include "pymol/align/ce_types.h"
#include "pymol/executive/executive_internal.h"

int ExecutiveAlign(PyMOLGlobals* G, const char* s1, const char* s2,
                   const char* mat_file, float gap, float extend, int max_gap,
                   int max_skip, float cutoff, int cycles, int quiet,
                   const char* oname, int state1, int state2,
                   ExecutiveRMSInfo* rms_info, int transform, int reset,
                   float seq_wt, float radius, float scale, float base,
                   float coord_wt, float expect, int window, float ante) {
  int sele1 = SelectorIndexByName(G, s1);
  int sele2 = SelectorIndexByName(G, s2);
  int* vla1 = nullptr;
  int* vla2 = nullptr;
  int na, nb;
  int c;
  int ok = true;
  int use_sequence = (mat_file && mat_file[0] && (seq_wt != 0.0F));
  int use_structure =
      (seq_wt >= 0.0F); /* negative seq_wt means sequence only! */
  ObjectMolecule* mobile_obj = nullptr;
  CMatch* match = nullptr;

  if (!use_structure) window = 0;

  if ((scale == 0.0F) && (seq_wt == 0.0F) && (ante < 0.0F) && window)
    ante = window;

  if (ante < 0.0F) ante = 0.0F;

  if ((sele1 >= 0)) {
    mobile_obj = SelectorGetSingleObjectMolecule(G, sele1);
    if (!mobile_obj) {
      ok = false;
      PRINTFB(G, FB_Executive, FB_Errors)
      " ExecutiveAlign: mobile selection must derive from one object "
      "only.\n" ENDFB(G);
    }
  }
  if (ok && (sele1 >= 0) && (sele2 >= 0) && rms_info) {
    vla1 = SelectorGetResidueVLA(G, sele1, use_structure, nullptr);
    vla2 = SelectorGetResidueVLA(G, sele2, use_structure, mobile_obj);
    if (vla1 && vla2) {
      na = VLAGetSize(vla1) / 3;
      nb = VLAGetSize(vla2) / 3;
      if (na && nb) {
        match = MatchNew(G, na, nb, window);
        if (match) {
          if (use_sequence) {
            if (ok) ok = MatchResidueToCode(match, vla1, na);
            if (ok) ok = MatchResidueToCode(match, vla2, nb);
            if (ok) ok = MatchMatrixFromFile(match, mat_file, quiet);
            if (ok) ok = MatchPreScore(match, vla1, na, vla2, nb, quiet);
          }
          if (use_structure) {
            /* avoid degenerate alignments */
            ok = ((na > 1) && (nb > 1) && ok);
            if (ok) {
              ok = SelectorResidueVLAsTo3DMatchScores(
                  G, match, vla1, na, state1, vla2, nb, state2, seq_wt, radius,
                  scale, base, coord_wt, expect);
            } else {
              PRINTFB(G, FB_Executive, FB_Errors)
              " ExecutiveAlign: No alignment found.\n" ENDFB(G);
            }
          }
          if (ok)
            ok = MatchAlign(match, gap, extend, max_gap, max_skip, quiet,
                            window, ante);
          if (ok) {
            rms_info->raw_alignment_score = match->score;
            rms_info->n_residues_aligned = match->n_pair;
            if (match->pair) {
              c = SelectorCreateAlignments(G, match->pair, sele1, vla1, sele2,
                                           vla2, "_align1", "_align2", false,
                                           false);

              if (c) {
                int mode = 2;
                if (!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                  " %s: %d atoms aligned.\n", __func__, c ENDFB(G);
                }
                if (oname && oname[0] && reset) ExecutiveDelete(G, oname);
                if (!transform) mode = 1;
                ok = ExecutiveRMS(G, "_align1", "_align2", mode, cutoff, cycles,
                                  quiet, oname, state1, state2, false, 0,
                                  rms_info);
              } else {
                if (!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                  " ExecutiveAlign-Error: atomic alignment failed (mismatched "
                  "identifiers?).\n" ENDFB(G);
                }
                ok = false;
              }
            }
          }
          MatchFree(match);
        }
      } else {
        ok = false;
        PRINTFB(G, FB_Executive, FB_Errors)
        " ExecutiveAlign: invalid selections for alignment.\n" ENDFB(G);
      }
    }
    ExecutiveUpdateCoordDepends(
        G, mobile_obj);  // Updates dynamic_measures - see PYMOL-3090
  }

  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return ok;
}

/**
 * Implementation of `cmd.find_pairs()`
 *
 * @param s1 atom selection expression
 * @param s2 atom selection expression
 * @param state1 object state
 * @param state2 object state
 * @param mode 0: find all pairs, 1: find hydrogen bonds
 * @param cutoff maximum distance in Angstrom
 * @param h_angle hydrogen bond angle cutoff for mode=1
 * @param[out] indexVLA pairs of atom indices (0-based)
 * @param[out] objVLA pairs of objects
 * @return number of pairs
 */
pymol::Result<int> ExecutivePairIndices(PyMOLGlobals* G, const char* s1,
                                        const char* s2, int state1, int state2,
                                        int mode, float cutoff, float h_angle,
                                        int** indexVLA,
                                        ObjectMolecule*** objVLA) {
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(2, sele1);

  return SelectorGetPairIndices(G, sele1, state1, sele2, state2, mode, cutoff,
                                h_angle, indexVLA, objVLA);
}

typedef struct {
  int priority;
  float vertex[3];
  AtomInfoType* ai;
} FitVertexRec;

static int fVertexOrdered(const FitVertexRec* array, int l, int r) {
  return (array[l].priority <= array[r].priority);
}

static int fAtomOrdered(PyMOLGlobals* G, const AtomInfoType* const* array,
                        int l, int r) {
  return (AtomInfoCompare(G, array[l], array[r]));
}

static int fAtomIDOrdered(const AtomInfoType* const* array, int l, int r) {
  return (array[l]->id <= array[r]->id);
}

static int fAtomRankOrdered(const AtomInfoType* const* array, int l, int r) {
  return (array[l]->rank <= array[r]->rank);
}

static int fAtomTemp1Ordered(const AtomInfoType* const* array, int l, int r) {
  return (array[l]->temp1 <= array[r]->temp1);
}

static void PackSortedIndices(int n, int* x, int rec_size, void* data) {
  int a;
  for (a = 0; a < n; a++) {
    if (a != x[a]) {
      memcpy(((char*)data) + (a * rec_size), ((char*)data) + (x[a] * rec_size),
             rec_size);
    }
  }
}
int ExecutiveRMS(PyMOLGlobals* G, const char* s1, const char* s2, int mode,
                 float refine, int max_cyc, int quiet, const char* oname,
                 int state1, int state2, int ordered_selections, int matchmaker,
                 ExecutiveRMSInfo* rms_info) {
  /* mode 0: measure rms without superposition
     mode 1: measure rms with trial superposition
     mode 2: measure rms with actual superposition */

  int sele1, sele2;
  float rms = -1.0;
  int a, b;
  float inv, *f, *f1, *f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int* flag;
  int ok = true;
  int repeat;
  float v1[3], *v2;
  ObjectAlignment* align_to_update = nullptr;

  bool ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  bool ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);

  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if (matrix_mode < 0) matrix_mode = 0; /* for now */

  if (matchmaker == -1) {
    /* matchmaker -1 is the same as matchmaker 0 except that the
       selections are not pre-matched prior to calling of this routine */
    matchmaker = 0;
  }

  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);

  /* this function operates on stored coordinates -- thus transformation
     matrices will need to be applied to the resulting atoms */

  // get coordinates
  {
    auto sele = sele1;
    auto state = state1;
    auto op = &op1;

    // for both selections
    do {
      ObjectMoleculeOpRecInit(op);

      if (sele >= 0) {
        if (state < 0) {
          op->code = OMOP_AVRT;
        } else {
          op->code = OMOP_StateVRT;
          op->i1 = state;
        }

        op->nvv1 = 0;  // length of vc1 (number of atoms with coordinates)
        op->vc1 = VLACalloc(int, 1000);    // number of states per atom
        op->vv1 = VLACalloc(float, 1000);  // coordinates (sum over states)

        if (mode == 0)
          op->i2 = true; /* if measuring current coordinates, then get global
                            txfd values */

        if (matchmaker || (oname && oname[0]))
          op->ai1VLA = VLACalloc(AtomInfoType*, 1000);

        if (ordered_selections)
          op->vp1 = VLAlloc(
              int, 1000);  // selection member "priority"? (MemberType::tag)

        ExecutiveObjMolSeleOp(G, sele, op);

        for (a = 0; a < op->nvv1; a++) {
          inv = (float)op->vc1[a]; /* average over coordinate sets */
          if (inv) {
            f = op->vv1 + (a * 3);
            scale3f(f, 1.F / inv, f);
          }
        }
      }

      // second iteration done
      if (sele == sele2) break;

      sele = sele2;
      state = state2;
      op = &op2;

    } while (true);
  }

  if (op1.vv1 && op2.vv1) {
    if (op1.nvv1 && op2.nvv1) {
      ObjectMolecule* mobile_obj = nullptr;

      int n_pair = 0;

      if (!(mobile_obj = SelectorGetSingleObjectMolecule(G, sele1))) {
        if (mode != 2) {
          PRINTFB(G, FB_Executive, FB_Warnings)
          "Executive-Warning: Mobile selection spans more than one "
          "object.\n" ENDFB(G);
        } else {
          PRINTFB(G, FB_Executive, FB_Errors)
          "Executive-Error: Mobile selection spans more than one object. "
          "Aborting.\n" ENDFB(G);
          ok = false;
        }
      }

      if (ok && op1.nvv1 && op2.nvv1 && (matchmaker > 0)) {
        /* matchmaker 0 is the default... by internal atom ordering only */
        int* idx1 = pymol::malloc<int>(op1.nvv1);
        int* idx2 = pymol::malloc<int>(op2.nvv1);
        int sort_flag = false;
        if (!(idx1 && idx2))
          ok = false;
        else {
          switch (matchmaker) {
            case 1: /* by atom info-based ordering */
              UtilSortIndexGlobals(G, op1.nvv1, op1.ai1VLA, idx1,
                                   (UtilOrderFnGlobals*)fAtomOrdered);
              UtilSortIndexGlobals(G, op2.nvv1, op2.ai1VLA, idx2,
                                   (UtilOrderFnGlobals*)fAtomOrdered);
              sort_flag = true;
              break;
            case 2: /* by matching atom identifiers */
              UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1,
                            (UtilOrderFn*)fAtomIDOrdered);
              UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2,
                            (UtilOrderFn*)fAtomIDOrdered);
              sort_flag = true;
              break;
            case 3: /* by matching atom ranks */
              UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1,
                            (UtilOrderFn*)fAtomRankOrdered);
              UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2,
                            (UtilOrderFn*)fAtomRankOrdered);
              sort_flag = true;
              break;
            case 4: /* by internal atom indexes (stored in temp1 kludge field)
                     */
              UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1,
                            (UtilOrderFn*)fAtomTemp1Ordered);
              UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2,
                            (UtilOrderFn*)fAtomTemp1Ordered);
              sort_flag = true;
              break;
          }
          if (sort_flag) {
            /* GOD this is SO ugly! */

            if (op1.vv1) {
              float* tmp = VLAlloc(float, op1.nvv1 * 3);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, 3 * sizeof(float),
                                       op1.vv1, tmp);
                VLAFreeP(op1.vv1);
                op1.vv1 = tmp;
              }
            }
            if (op1.vc1) {
              int* tmp = VLAlloc(int, op1.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vc1,
                                       tmp);
                VLAFreeP(op1.vc1);
                op1.vc1 = tmp;
              }
            }
            if (op1.vp1) {
              int* tmp = VLAlloc(int, op1.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vp1,
                                       tmp);
                VLAFreeP(op1.vp1);
                op1.vp1 = tmp;
              }
            }
            if (op1.ai1VLA) {
              AtomInfoType** tmp = VLACalloc(AtomInfoType*, op1.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(AtomInfoType*),
                                       op1.ai1VLA, tmp);
                VLAFreeP(op1.ai1VLA);
                op1.ai1VLA = tmp;
              }
            }

            if (op2.vv1) {
              float* tmp = VLAlloc(float, op2.nvv1 * 3);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, 3 * sizeof(float),
                                       op2.vv1, tmp);
                VLAFreeP(op2.vv1);
                op2.vv1 = tmp;
              }
            }
            if (op2.vc1) {
              int* tmp = VLAlloc(int, op2.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vc1,
                                       tmp);
                VLAFreeP(op2.vc1);
                op2.vc1 = tmp;
              }
            }
            if (op2.vp1) {
              int* tmp = VLAlloc(int, op2.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vp1,
                                       tmp);
                VLAFreeP(op2.vp1);
                op2.vp1 = tmp;
              }
            }
            if (op2.ai1VLA) {
              AtomInfoType** tmp = VLACalloc(AtomInfoType*, op2.nvv1);
              if (!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(AtomInfoType*),
                                       op2.ai1VLA, tmp);
                VLAFreeP(op2.ai1VLA);
                op2.ai1VLA = tmp;
              }
            }
          }
        }

        if (matchmaker != 0) {
          int n1 = 0, n2 = 0, c1 = 0, c2 = 0;
          int cmp;

          while ((n1 < op1.nvv1) && (n2 < op2.nvv1)) {
            cmp = 0;
            switch (matchmaker) {
              case 1: /* insure that AtomInfoType matches */
                if (AtomInfoMatch(G, op1.ai1VLA[n1], op2.ai1VLA[n2],
                                  ignore_case, ignore_case_chain))
                  cmp = 0;
                else
                  cmp = AtomInfoCompare(G, op1.ai1VLA[n1], op2.ai1VLA[n2]);
                printf("%d-%d %d-%d: %d\n", c1, n1, c2, n2, cmp);
                break;
              case 2: /* ID */
              case 3: /* rank */
              {
                int val1;
                int val2;

                switch (matchmaker) {
                  case 2: /* ID */
                    val1 = op1.ai1VLA[n1]->id;
                    val2 = op2.ai1VLA[n2]->id;
                    break;
                  case 3: /* rank */
                    val1 = op1.ai1VLA[n1]->rank;
                    val2 = op2.ai1VLA[n2]->rank;
                    break;
                  case 4: /* index (via temp1) */
                    val1 = op1.ai1VLA[n1]->temp1;
                    val2 = op2.ai1VLA[n2]->temp1;
                    break;
                  default:
                    val1 = 0;
                    val2 = 0;
                    break;
                }
                if (val1 == val2)
                  cmp = 0;
                else if (val1 < val2)
                  cmp = -1;
                else
                  cmp = 1;
              } break;
            }
            if (!cmp) { /* match found */
              idx1[c1++] = n1++;
              idx2[c2++] = n2++;
              n_pair++;
            } else if (cmp < 0) { /* op1 below op2 */
              n1++;
            } else { /* op2 below op1 */
              n2++;
            }
          }

          if (n_pair) {
            if (op1.vv1)
              PackSortedIndices(n_pair, idx1, 3 * sizeof(float), op1.vv1);
            if (op1.vc1) PackSortedIndices(n_pair, idx1, sizeof(int), op1.vc1);
            if (op1.vp1) PackSortedIndices(n_pair, idx1, sizeof(int), op1.vp1);
            if (op1.ai1VLA)
              PackSortedIndices(n_pair, idx1, sizeof(AtomInfoType*),
                                op1.ai1VLA);

            if (op2.vv1)
              PackSortedIndices(n_pair, idx2, 3 * sizeof(float), op2.vv1);
            if (op2.vc1) PackSortedIndices(n_pair, idx2, sizeof(int), op2.vc1);
            if (op2.vp1) PackSortedIndices(n_pair, idx2, sizeof(int), op2.vp1);
            if (op2.ai1VLA)
              PackSortedIndices(n_pair, idx2, sizeof(AtomInfoType*),
                                op2.ai1VLA);
          }
        }
        FreeP(idx1);
        FreeP(idx2);
      } else if (op1.nvv1 != op2.nvv1) {
        sprintf(buffer, "Atom counts between selections don't match (%d vs %d)",
                op1.nvv1, op2.nvv1);
        ErrMessage(G, "ExecutiveRMS", buffer);
        n_pair = 0;
        ok = false;
      } else {
        n_pair = 0;
        for (a = 0; a < op1.nvv1; ++a) {   // for atoms in selection
          if (op1.vc1[a] && op2.vc1[a]) {  // check state counts
            if (n_pair < a) {              // copy over if necessary
              copy3(op1.vv1 + 3 * a, op1.vv1 + 3 * n_pair);
              copy3(op2.vv1 + 3 * a, op2.vv1 + 3 * n_pair);
              if (op1.ai1VLA) op1.ai1VLA[n_pair] = op1.ai1VLA[a];
              if (op2.ai1VLA) op2.ai1VLA[n_pair] = op2.ai1VLA[a];
              if (op1.vp1) op1.vp1[n_pair] = op1.vp1[a];
              if (op2.vp1) op2.vp1[n_pair] = op2.vp1[a];
              op1.vc1[n_pair] = op1.vc1[a];
              op2.vc1[n_pair] = op2.vc1[a];
            }
            ++n_pair;
          }
        }
      }

      if (n_pair) {
        /* okay -- we're on track to do an alignment */

        if (ordered_selections && op1.vp1 && op2.vp1) {
          /* if we expected ordered selections and have priorities,
             then we may need to sort vertices */

          int sort_flag1 = false, sort_flag2 = false;
          int well_defined1 = true, well_defined2 = true;

          for (a = 0; a < (n_pair - 1); a++) {
            /*          printf("op1 vertex %d priority %d\n",a,op1.vp1[a]);
               printf("op2 vertex %d priority %d\n",a,op2.vp1[a]); */

            if (op1.vp1[a] > op1.vp1[a + 1])
              sort_flag1 = true;
            else if (op1.vp1[a] == op1.vp1[a + 1])
              well_defined1 = false;
            if (op2.vp1[a] > op2.vp1[a + 1])
              sort_flag2 = true;
            else if (op2.vp1[a] == op2.vp1[a + 1])
              well_defined2 = false;
          }

          if (sort_flag1 || sort_flag2) {
            if (!(well_defined1 || well_defined2)) {
              PRINTFB(G, FB_Executive, FB_Warnings)
              "Executive-Warning: Ordering requested but not well "
              "defined.\n" ENDFB(G);
            } else {
              FitVertexRec* vert = pymol::malloc<FitVertexRec>(n_pair);

              if (sort_flag1) {
                float *src, *dst;
                src = op1.vv1;
                for (a = 0; a < n_pair; a++) {
                  vert[a].priority = op1.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn*)fVertexOrdered);
                dst = op1.vv1;
                for (a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              if (sort_flag2) {
                float *src, *dst;
                src = op2.vv1;
                for (a = 0; a < n_pair; a++) {
                  vert[a].priority = op2.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn*)fVertexOrdered);
                dst = op2.vv1;
                for (a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              FreeP(vert);
            }
          }
        }

        if (rms_info) {
          rms_info->initial_n_atom = n_pair;
          rms_info->n_cycles_run = 0;
          rms_info->final_n_atom = n_pair; /* in case there is no refinement */
        }

        if (mode != 0) {
          rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, nullptr, op2.ttt);
          if (rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
          repeat = true;
          b = 0;
          while (repeat) {
            repeat = false;
            b++;
            if (b > max_cyc) break;
            if ((refine > R_SMALL4) && (rms > R_SMALL4)) {
              int n_next = n_pair;
              AtomInfoType **ai1, **ai2;

              flag = pymol::malloc<int>(n_pair);

              if (flag) {
                for (a = 0; a < n_pair; a++) {
                  MatrixTransformTTTfN3f(1, v1, op2.ttt, op1.vv1 + (a * 3));
                  v2 = op2.vv1 + (a * 3);
                  if ((diff3f(v1, v2) / rms) > refine) {
                    flag[a] = false;
                    repeat = true;
                  } else
                    flag[a] = true;
                }
                f1 = op1.vv1;
                f2 = op2.vv1;
                ai1 = op1.ai1VLA;
                ai2 = op2.ai1VLA;
                for (a = 0; a < n_pair; a++) {
                  if (!flag[a]) {
                    n_next--;
                  } else {
                    copy3f(op1.vv1 + (3 * a), f1);
                    copy3f(op2.vv1 + (3 * a), f2);
                    f1 += 3;
                    f2 += 3;
                    if (ai1 && ai2) { /* make sure we keep track of which atoms
                                         are aligned */
                      *(ai1++) = op1.ai1VLA[a];
                      *(ai2++) = op2.ai1VLA[a];
                    }
                  }
                }
                if (!quiet && (n_next != n_pair)) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                  " %s: %d atoms rejected during cycle %d (RMSD=%0.2f).\n",
                      __func__, n_pair - n_next, b, rms ENDFB(G);
                }
                n_pair = n_next;
                FreeP(flag);
                if (n_pair) {
                  rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, nullptr,
                                         op2.ttt);
                  if (rms_info) {
                    rms_info->n_cycles_run = b;
                    rms_info->final_n_atom = n_pair;
                    rms_info->final_rms = rms;
                  }
                } else
                  break;
              }
            }
          }
        } else { /* mode == 0 -- simple RMS, with no coordinate movement */
          rms = MatrixGetRMS(G, n_pair, op1.vv1, op2.vv1, nullptr);
          if (rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
        }
      }
      if (!n_pair) {
        PRINTFB(G, FB_Executive, FB_Results)
        " Executive: Error -- no atoms left after refinement!\n" ENDFB(G);
        ok = false;
      }

      if (ok) {
        if (!quiet) {
          PRINTFB(G, FB_Executive, FB_Results)
          " Executive: RMSD = %8.3f (%d to %d atoms)\n", rms, n_pair,
              n_pair ENDFB(G);
        }
        if (oname && oname[0]) {
          int align_state = state2;
          ObjectMolecule* trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

          if (align_state < 0) {
            align_state = SceneGetState(G);
          }

          /* we're going to create/update an alignment object */

          {
            /* Get unique ids and construct the alignment vla */
            pymol::vla<int> align_vla(n_pair * 3);

            {
              int* id_p = align_vla.data();
              int i;
              for (i = 0; i < n_pair; i++) {
                id_p[0] = AtomInfoCheckUniqueID(G, op2.ai1VLA[i]); /* target */
                id_p[1] = AtomInfoCheckUniqueID(G, op1.ai1VLA[i]);
                id_p[2] = 0;
                id_p += 3;
              }
              VLASize(align_vla, int, n_pair * 3);
            }
            {
              ObjectAlignment* obj = nullptr;

              /* does object already exist? */
              {
                pymol::CObject* execObj = ExecutiveFindObjectByName(G, oname);
                if (execObj && (execObj->type != cObjectAlignment))
                  ExecutiveDelete(G, oname);
                else
                  obj = (ObjectAlignment*)execObj;
              }
              obj = ObjectAlignmentDefine(G, obj, align_vla, align_state, true,
                                          trg_obj, mobile_obj);
              obj->Color = ColorGetIndex(G, "yellow");
              ObjectSetName(obj, oname);
              ExecutiveManageObject(G, obj, 0, quiet);
              align_to_update = obj;
              SceneInvalidate(G);
            }
          }
        }
        if (ok && mode == 2) {
          if (matrix_mode > 0) {
            ObjectMolecule *src_obj, *trg_obj;
            src_obj = SelectorGetFirstObjectMolecule(
                G, sele1); /* get at least one object */
            trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

            /* first we need to make sure that the object being moved
               matches the target with respect to both the TTT and the
               object's state matrix (if any) */

            if (src_obj && trg_obj) {
              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 1,
                                  1, /* TTT mode */
                                  state2, state1, false, 0, quiet);

              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 2,
                                  2, /* Object state mode */
                                  state2, state1, false, 0, quiet);

              switch (matrix_mode) {
                case 1: /* TTTs */
                  ExecutiveCombineObjectTTT(G, src_obj->Name, op2.ttt, true,
                                            -1);
                  break;
                case 2: {
                  double homo[16], *src_homo;
                  convertTTTfR44d(op2.ttt, homo);
                  if (ExecutiveGetObjectMatrix(G, src_obj->Name, state1,
                                               &src_homo, false)) {
                    left_multiply44d44d(src_homo, homo);
                    ExecutiveSetObjectMatrix(G, src_obj->Name, state1, homo);
                  }
                } break;
              }
              /* next we need to update the object's TTT matrix to reflect
                 the transformation */
            }
          } else { /* matrix_mode is zero -- legacy behavior */
            /* this will transform the actual coordinates */
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(G, sele1, &op2);
          }
        }
      }
    } else {
      ErrMessage(G, __func__, "No atoms selected.");
      ok = false;
    }
  }

  if (align_to_update) {
    align_to_update->update();
  }

  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  VLAFreeP(op1.vp1);
  VLAFreeP(op2.vp1);
  VLAFreeP(op1.ai1VLA);
  VLAFreeP(op2.ai1VLA);
  return (ok);
}
/**
 * Fit states or calculate ensemble RMSD
 *
 * @param s1 atom selection expression
 * @param target reference state
 * @param mode 2=intra_fit, 1=intra_rms, 0=intra_rms_cur
 * @param mix intra_fit only, average the prior target coordinates
 * @param pbc Consider periodic boundary conditions
 */
pymol::Result<pymol::vla<float>> ExecutiveRMSStates(PyMOLGlobals* G,
                                                    const char* s1, int target,
                                                    int mode, int quiet,
                                                    int mix, bool pbc) {
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float* result = nullptr;
  int ok = true;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.vv1 = nullptr;
  op2.vv1 = nullptr;

  ObjectMolecule* obj = SelectorGetSingleObjectMolecule(G, sele1);

  if (!obj) {
    if (mode != 2) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      "Executive-Warning: Mobile selection spans more than one object.\n" ENDFB(
          G);
    } else {
      return pymol::make_error("Mobile selection spans more than one object.");
    }
  }

  if (target == cStateCurrent) {
    target = obj ? obj->getCurrentState() : SceneGetState(G);
  }

  if (target < 0) {
    target = 0;
  }

  if (mode != 2) {
    pbc = false;
  }

  if (ok && sele1 >= 0) {
    op1.code = OMOP_SVRT;
    op1.nvv1 = 0;
    op1.i1 = target;
    op1.vv1 = (float*)VLAMalloc(1000, sizeof(float), 5, 0);
    op1.i1VLA = VLAlloc(int, 1000);
    ExecutiveObjMolSeleOp(G, sele1, &op1);

    if (pbc) {
      ObjectMoleculePBCUnwrap(*obj);
    }

    op2.vv2 = op1.vv1;
    op2.nvv2 = op1.nvv1;
    op2.i1VLA = op1.i1VLA;
    op2.i2 = target;
    op2.i1 = mode;
    op2.i3 = mix;
    op2.f1VLA = VLAlloc(float, 10);
    VLASize(op2.f1VLA, float, 0); /* failsafe */
    op2.vv1 = (float*)VLAMalloc(1000, sizeof(float), 5, 0);
    op2.code = OMOP_SFIT;
    op2.nvv1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    result = op2.f1VLA;
    VLAFreeP(op1.i1VLA);
    VLAFreeP(op2.vv1);

    if (pbc) {
      float center[3];
      pymol::meanNx3(op1.vv1, op1.nvv1, center);
      ObjectMoleculePBCWrap(*obj, center);
    }

    VLAFreeP(op1.vv1);

    if (mode == 2) {
      ExecutiveUpdateCoordDepends(G, obj);
    }
  }
  return pymol::vla_take_ownership(result);
}
float ExecutiveRMSPairs(PyMOLGlobals* G, const std::vector<SelectorTmp>& sele,
                        int mode, bool quiet) {
  int sele1, sele2;
  int a, c;
  float rms = -1.0, inv, *f;
  OrthoLineType buffer;

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType combi, s1;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.nvv1 = 0;
  op1.vc1 = (int*)VLAMalloc(1000, sizeof(int), 5, 1);
  op1.vv1 = (float*)VLAMalloc(1000, sizeof(float), 5, 1); /* auto-zero */
  op1.code = OMOP_AVRT;

  op2.nvv1 = 0;
  op2.vc1 = (int*)VLAMalloc(1000, sizeof(int), 5, 1);
  op2.vv1 = (float*)VLAMalloc(1000, sizeof(float), 5, 1); /* auto-zero */
  op2.code = OMOP_AVRT;

  strcpy(combi, "(");
  c = 0;
  auto pairs = sele.size() / 2;
  for (a = 0; a < pairs; a++) {
    sele1 = sele[c].getIndex();
    if (sele1 >= 0) ExecutiveObjMolSeleOp(G, sele1, &op1);
    strcat(combi, sele[c].getName());
    if (a < (pairs - 1)) strcat(combi, " or ");
    c++;
    sele2 = sele[c].getIndex();
    if (sele2 >= 0) ExecutiveObjMolSeleOp(G, sele2, &op2);
    c++;
  }
  strcat(combi, ")");
  for (a = 0; a < op1.nvv1; a++) {
    inv = (float)op1.vc1[a];
    if (inv) {
      f = op1.vv1 + (a * 3);
      inv = 1.0F / inv;
      *(f++) *= inv;
      *(f++) *= inv;
      *(f++) *= inv;
    }
  }
  for (a = 0; a < op2.nvv1; a++) {
    inv = (float)op2.vc1[a];
    if (inv) {
      f = op2.vv1 + (a * 3);
      inv = 1.0F / inv;
      *(f++) *= inv;
      *(f++) *= inv;
      *(f++) *= inv;
    }
  }
  if (op1.vv1 && op2.vv1) {
    if (op1.nvv1 != op2.nvv1) {
      sprintf(buffer,
              "Atom counts between selection sets don't match (%d != %d).",
              op1.nvv1, op2.nvv1);
      ErrMessage(G, __func__, buffer);
    } else if (op1.nvv1) {
      if (mode != 0)
        rms = MatrixFitRMSTTTf(G, op1.nvv1, op1.vv1, op2.vv1, nullptr, op2.ttt);
      else
        rms = MatrixGetRMS(G, op1.nvv1, op1.vv1, op2.vv1, nullptr);

      if (!quiet) PRINTFB(G, FB_Executive, FB_Results)
      " %s: RMSD = %8.3f (%d to %d atoms)\n", __func__, rms, op1.nvv1,
          op2.nvv1 ENDFB(G);

      op2.code = OMOP_TTTF;
      SelectorGetTmp(G, combi, s1);
      sele1 = SelectorIndexByName(G, s1);
      ExecutiveObjMolSeleOp(G, sele1, &op2);
      SelectorFreeTmp(G, s1);
    } else {
      ErrMessage(G, __func__, "No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return (rms);
}

PyObject* ExecutiveCEAlign(PyMOLGlobals* G, PyObject* listA, PyObject* listB,
                           int lenA, int lenB, float d0, float d1,
                           int windowSize, int gapMax) {
  int i = 0;
  int smaller;
  double **dmA, **dmB, **S;
  int bufferSize;
  pcePoint coordsA, coordsB;
  pathCache paths = nullptr;
  PyObject* result;

  smaller = lenA < lenB ? lenA : lenB;

  /* get the coodinates from the Python objects */
  coordsA = (pcePoint)getCoords(listA, lenA);
  coordsB = (pcePoint)getCoords(listB, lenB);

  /* calculate the distance matrix for each protein */
  dmA = (double**)calcDM(coordsA, lenA);
  dmB = (double**)calcDM(coordsB, lenB);

  /* calculate the CE Similarity matrix */
  S = (double**)calcS(dmA, dmB, lenA, lenB, windowSize);

  /* find the best path through the CE Sim. matrix */
  bufferSize = 0;

  /* the following line HANGS PyMOL */
  paths = (pathCache)findPath(S, dmA, dmB, lenA, lenB, d0, d1, windowSize,
                              gapMax, &bufferSize);

  /* Get the optimal superposition here... */
  result = (PyObject*)findBest(coordsA, coordsB, paths, bufferSize, smaller,
                               windowSize);

  /* release memory */
  free(coordsA);
  free(coordsB);
  for (i = 0; i < bufferSize; ++i) free(paths[i]);
  free(paths);

  /* distance matrices	 */
  for (i = 0; i < lenA; i++) free(dmA[i]);
  free(dmA);

  for (i = 0; i < lenB; i++) free(dmB[i]);
  free(dmB);

  /* similarity matrix */
  for (i = 0; i < lenA; i++) free(S[i]);
  free(S);

  return (PyObject*)result;
}

pymol::Result<ExecutiveRMSInfo> ExecutiveFit(
    PyMOLGlobals* G, pymol::zstring_view str1, pymol::zstring_view str2,
    int mode, int cutoff, int cycles, int quiet, pymol::zstring_view object,
    int state1, int state2, int matchmaker) {
  SelectorTmp s1(G, str1.c_str());
  SelectorTmp s2(G, str2.c_str());
  ExecutiveRMSInfo rmsinfo;
  ExecutiveRMS(G, s1.getName(), s2.getName(), mode, cutoff, cycles, quiet,
               object.c_str(), state1, state2, false, matchmaker, &rmsinfo);
  return rmsinfo;
}

/**
 * Run TM-align on two selections and return results.
 *
 * @param mobile_sele mobile selection (will be transformed)
 * @param target_sele target selection (stays fixed)
 * @param mobile_state state of mobile selection (0-based)
 * @param target_state state of target selection (0-based)
 * @param quiet suppress output
 * @param transform apply superposition transform
 * @param oname name for alignment object (empty = don't create)
 * @param fast use fast mode (fewer iterations)
 */
pymol::Result<pymol::usalign::TMAlignResult> ExecutiveUSalign(
    PyMOLGlobals* G, const char* mobile_sele, const char* target_sele,
    int mobile_state, int target_state, int quiet, int transform,
    const char* oname, int fast) {
  // Resolve selections
  auto sele_mobile = SelectorIndexByName(G, mobile_sele);
  if (sele_mobile < 0)
    return pymol::make_error("Invalid mobile selection: ", mobile_sele);

  auto sele_target = SelectorIndexByName(G, target_sele);
  if (sele_target < 0)
    return pymol::make_error("Invalid target selection: ", target_sele);

  // Extract CA coordinates and sequences
  struct ResidueInfo {
    glm::dvec3 coord;
    char seq_char;
    AtomInfoType* ai;
  };

  auto extract_ca = [&](SelectorID_t sele,
                        int state) -> std::vector<ResidueInfo> {
    std::vector<ResidueInfo> residues;
    SeleCoordIterator iter(G, sele, state);
    while (iter.next()) {
      auto* ai = iter.getAtomInfo();
      if (ai->flags & cAtomFlag_guide) {
        float* c = iter.getCoord();
        ResidueInfo ri;
        ri.coord = glm::dvec3(c[0], c[1], c[2]);
        ri.seq_char = SeekerGetAbbr(G, LexStr(G, ai->resn), 'O', 'X');
        ri.ai = ai;
        residues.push_back(ri);
      }
    }
    return residues;
  };

  auto mobile_res = extract_ca(sele_mobile, mobile_state);
  auto target_res = extract_ca(sele_target, target_state);

  if (mobile_res.size() < 3) {
    return pymol::make_error("Mobile selection has fewer than 3 guide atoms (",
                             mobile_res.size(), ")");
  }
  if (target_res.size() < 3) {
    return pymol::make_error("Target selection has fewer than 3 guide atoms (",
                             target_res.size(), ")");
  }

  // Build coordinate vectors and sequences
  std::vector<glm::dvec3> mobile_ca, target_ca;
  std::string mobile_seq, target_seq;
  mobile_ca.reserve(mobile_res.size());
  target_ca.reserve(target_res.size());

  for (const auto& r : mobile_res) {
    mobile_ca.push_back(r.coord);
    mobile_seq.push_back(r.seq_char);
  }
  for (const auto& r : target_res) {
    target_ca.push_back(r.coord);
    target_seq.push_back(r.seq_char);
  }

  // Run TM-align
  auto result = pymol::usalign::TMalign(target_ca, mobile_ca, target_seq,
                                        mobile_seq, fast != 0);

  if (result.aligned_length < 1) {
    return pymol::make_error("TM-align failed to find any alignment");
  }

  // Print results
  if (!quiet) {
    PRINTFB(G, FB_Executive, FB_Results)
    " USalign: TM-score= %6.4f (normalized by target, N=%d, d0=%.2f)\n",
        result.tm_score_target, static_cast<int>(target_ca.size()),
        result.d0_target ENDFB(G);
    PRINTFB(G, FB_Executive, FB_Results)
    " USalign: TM-score= %6.4f (normalized by mobile, N=%d, d0=%.2f)\n",
        result.tm_score_mobile, static_cast<int>(mobile_ca.size()),
        result.d0_mobile ENDFB(G);
    PRINTFB(G, FB_Executive, FB_Results)
    " USalign: Aligned length= %d, RMSD= %5.2f, Seq_ID=n_identical/n_aligned= "
    "%4.3f\n",
        result.aligned_length, result.rmsd, result.seq_identity ENDFB(G);
  }

  // Apply transform to mobile object
  if (transform) {
    // Convert double-precision Superposition to float TTT
    // USalign convention: y_aligned = R * x + t
    // where x = mobile coords, y = target coords
    // The rotation R and translation t transform mobile -> target space

    const auto& sup = result.transform;

    // Build a legacy-style 16-float TTT matrix
    // TTT format: [R00 R01 R02 pre_x] [R10 R11 R12 pre_y]
    //             [R20 R21 R22 pre_z] [tx  ty  tz  1]
    // where pre is the pre-translation (origin), and t is post-translation
    // For a simple rotation+translation (no origin): pre=0, R=rotation,
    // t=translation

    glm::mat3 rot_f(sup.rotation);
    glm::quat q = glm::quat_cast(rot_f);
    glm::vec3 t(sup.translation);

    // Create TTT: pretranslate=0, rotate=q, posttranslate=t
    pymol::TTT ttt(glm::vec3(0.0f), q, t);

    // Convert to legacy float[16] format for ExecuteCombineObjectTTT
    auto legacy = pymol::TTT::as_pymol_2_legacy(ttt);
    float tttf[16];
    std::memcpy(tttf, glm::value_ptr(legacy), 16 * sizeof(float));

    // Follow the same pattern as ExecutiveAlign:
    // 1. Copy target's TTT and state matrix to mobile (reset to same frame)
    // 2. Combine the alignment transform (reverse_order=true)
    // Note: Only the first object in the mobile selection is transformed,
    // matching ExecutiveAlign behavior for multi-object selections.
    ObjectMolecule* mobile_obj = SelectorGetFirstObjectMolecule(G, sele_mobile);
    ObjectMolecule* target_obj =
        SelectorGetSingleObjectMolecule(G, sele_target);
    if (mobile_obj && target_obj) {
      ExecutiveMatrixCopy(G, target_obj->Name, mobile_obj->Name, 1, 1,
                          target_state, mobile_state, false, 0, quiet);
      ExecutiveMatrixCopy(G, target_obj->Name, mobile_obj->Name, 2, 2,
                          target_state, mobile_state, false, 0, quiet);
      ExecutiveCombineObjectTTT(G, mobile_obj->Name, tttf, true, -1);
    }
  }

  // Create alignment object
  if (oname && oname[0]) {
    int align_state = target_state;
    if (align_state < 0) {
      align_state = SceneGetState(G);
    }

    ObjectMolecule* trg_obj = SelectorGetSingleObjectMolecule(G, sele_target);
    ObjectMolecule* mob_obj = SelectorGetFirstObjectMolecule(G, sele_mobile);

    if (trg_obj && mob_obj) {
      int n_pair = result.aligned_length;
      pymol::vla<int> align_vla(n_pair * 3);
      int* id_p = align_vla.data();

      for (int k = 0; k < n_pair; k++) {
        int mi = result.mobile_indices[k];
        int ti = result.target_indices[k];
        if (mi < static_cast<int>(mobile_res.size()) &&
            ti < static_cast<int>(target_res.size())) {
          id_p[0] = AtomInfoCheckUniqueID(G, target_res[ti].ai);
          id_p[1] = AtomInfoCheckUniqueID(G, mobile_res[mi].ai);
          id_p[2] = 0;
          id_p += 3;
        }
      }

      ObjectAlignment* obj = nullptr;
      {
        pymol::CObject* execObj = ExecutiveFindObjectByName(G, oname);
        if (execObj && execObj->type != cObjectAlignment) {
          ExecutiveDelete(G, oname);
        } else {
          obj = dynamic_cast<ObjectAlignment*>(execObj);
        }
      }
      obj = ObjectAlignmentDefine(G, obj, align_vla, align_state, true, trg_obj,
                                  mob_obj);
      obj->Color = ColorGetIndex(G, "yellow");
      ObjectSetName(obj, oname);
      ExecutiveManageObject(G, obj, 0, quiet);
      SceneInvalidate(G);
    }
  }

  return result;
}
