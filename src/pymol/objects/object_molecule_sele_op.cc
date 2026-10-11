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

// Operations on the atoms in a selection (ObjectMoleculeSeleOp).

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/geometry/matrix.h"
#include "pymol/python/p.h"
#include "pymol/executive/executive.h"
#include "pymol/chem/hydrogen_adder.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_molecule_internal.h"
int ObjectMoleculeTransformSelection(ObjectMolecule* I, int state, int sele,
                                     const float* matrix, int log,
                                     const char* sname, int homogenous,
                                     int global) {
  /* called from "translate [5,5,5], objSele" */
  /* if sele == -1, then the whole object state is transformed */
  PyMOLGlobals* G = I->G;
  int a, s;
  int flag = false;
  CoordSet* cs;
  const AtomInfoType* ai;
  int logging;
  int all_states = false, inp_state;
  int ok = true;
  float homo_matrix[16], tmp_matrix[16];
  const float* input_matrix = matrix;

  inp_state = state;
  if (state == -2) state = ObjectGetCurrentState(I, false);
  if (state < 0) {
    all_states = true;
    state = -1;
  }
  PRINTFD(G, FB_ObjectMolecule)
  "ObjMolTransSele-Debug: state %d\n", state ENDFD;
  while (1) {
    if (all_states) {
      state++;
      if (state >= I->NCSet) break;
    }
    if (state < I->NCSet) {
      cs = I->CSet[state];
      if (cs) {
        int use_matrices =
            SettingGet_i(G, I->Setting.get(), nullptr, cSetting_matrix_mode);
        if (use_matrices < 0) use_matrices = 0;

        if (global && !homogenous) { /* convert matrix to homogenous */
          convertTTTfR44f(matrix, homo_matrix);
          matrix = homo_matrix;
          input_matrix = homo_matrix;
          homogenous = true;
        }

        if (global && ((use_matrices && !cs->Matrix.empty()) || I->TTTFlag)) {
          /* if input coordinates are in the global system,
             they may need to be converted to local coordinates */

          matrix = input_matrix;

          /* global to object */

          if (I->TTTFlag) {
            float ttt[16];
            if (matrix != tmp_matrix) {
              copy44f(matrix, tmp_matrix);
            }
            convertTTTfR44f(I->TTT, ttt);
            {
              float ttt_inv[16];
              invert_special44f44f(ttt, ttt_inv);
              left_multiply44f44f(ttt_inv, tmp_matrix);
              right_multiply44f44f(tmp_matrix, ttt);
            }
            matrix = tmp_matrix;
          }

          /* object to state */

          if (use_matrices) {
            if (!cs->Matrix.empty()) {
              double tmp_double[16], tmp_inv[16];
              copy44f44d(matrix, tmp_double);
              invert_special44d44d(cs->Matrix.data(), tmp_inv);
              left_multiply44d44d(tmp_inv, tmp_double);
              right_multiply44d44d(tmp_double, cs->Matrix.data());
              copy44d44f(tmp_double, tmp_matrix);
              matrix = tmp_matrix;
            }
          }
        }

        if (sele >= 0) { /* transforming select atoms */
          ai = I->AtomInfo;
          for (a = 0; a < I->NAtom; a++) {
            s = ai->selEntry;
            if (!(ai->protekted == cAtomProtected_explicit))
              if (SelectorIsMember(G, s, sele)) {
                if (homogenous)
                  CoordSetTransformAtomR44f(cs, a, matrix);
                else
                  CoordSetTransformAtomTTTf(cs, a, matrix);
                flag = true;
              }
            ai++;
          }
        } else { /* transforming whole coordinate set */
          if (!use_matrices) {
            ai = I->AtomInfo;
            for (a = 0; a < I->NAtom; a++) {
              if (!(ai->protekted == cAtomProtected_explicit)) {
                if (homogenous)
                  CoordSetTransformAtomR44f(cs, a, matrix);
                else
                  CoordSetTransformAtomTTTf(cs, a, matrix);
              }
              ai++;
            }
            flag = true;
            CoordSetRecordTxfApplied(cs, matrix, homogenous);
          } else {
            ObjectMoleculeTransformState44f(I, state, matrix, false, homogenous,
                                            false);
          }
        }
        if (flag) {
          cs->invalidateRep(cRepAll, cRepInvCoord);
          ExecutiveUpdateCoordDepends(G, I);
        }
      }
    }
    if (!all_states) break;
  }

  if (log) {
    OrthoLineType line;
    WordType sele_str = ",'";
    logging = SettingGetGlobal_i(G, cSetting_logging);
    if (sele >= 0) {
      strcat(sele_str, sname);
    }
    strcat(sele_str, "'");
    switch (logging) {
      case cPLog_pml:
        sprintf(line,
                "_ cmd.transform_object('%s',[\\\n_ "
                "%15.9f,%15.9f,%15.9f,%15.9f,\\\n_ "
                "%15.9f,%15.9f,%15.9f,%15.9f,\\\n_ "
                "%15.9f,%15.9f,%15.9f,%15.9f,\\\n_ "
                "%15.9f,%15.9f,%15.9f,%15.9f\\\n_     ],%d,%d%s,%d)\n",
                I->Name, matrix[0], matrix[1], matrix[2], matrix[3], matrix[4],
                matrix[5], matrix[6], matrix[7], matrix[8], matrix[9],
                matrix[10], matrix[11], matrix[12], matrix[13], matrix[14],
                matrix[15], inp_state + 1, 0, sele_str, homogenous);
        PLog(G, line, cPLog_no_flush);
        break;
      case cPLog_pym:

        sprintf(line,
                "cmd.transform_object('%s',[\n%15.9f,%15.9f,%15.9f,%15.9f,\n%"
                "15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f,\n%"
                "15.9f,%15.9f,%15.9f,%15.9f\n],%d,%d%s,%d)\n",
                I->Name, matrix[0], matrix[1], matrix[2], matrix[3], matrix[4],
                matrix[5], matrix[6], matrix[7], matrix[8], matrix[9],
                matrix[10], matrix[11], matrix[12], matrix[13], matrix[14],
                matrix[15], inp_state + 1, 0, sele_str, homogenous);
        PLog(G, line, cPLog_no_flush);
        break;
      default:
        break;
    }
  }
  return (ok);
}
bool ObjectMoleculeSeleOp(ObjectMolecule* I, int sele,
                          ObjectMoleculeOpRec* op) {
  float* coord;
  int a, b, s;
  int c, t_i;
  int a1 = 0, ind;
  float rms;
  float v1[3], v2, *vv1, *vv2, *vt, *vt1, *vt2;
  int hit_flag = false;
  int ok = true;
  int cnt;
  int skip_flag;
  int match_flag = false;
  int offset;
  int priority;
  int use_matrices = false;
  CoordSet* cs;
  AtomInfoType *ai, *ai0;
  PyMOLGlobals* G = I->G;
  PyObject* expr_co = nullptr;
  int compileType = Py_single_input;
  PRINTFD(G, FB_ObjectMolecule)
  " %s-DEBUG: sele %d op->code %d\n", __func__, sele, op->code ENDFD;
  if (sele >= 0) {
    const char* errstr = "Alter";
    /* always run on entry */
    switch (op->code) {
      case OMOP_LABL:
        errstr = "Label";
        if (op->i2 != cExecutiveLabelEvalOn) {
          break;
        }
        compileType = Py_eval_input;
      case OMOP_ALTR:
      case OMOP_AlterState:
        // assume blocked interpreter

        if (op->s1 && op->s1[0]) {
          expr_co = Py_CompileString(op->s1, "", compileType);
          if (expr_co == nullptr) {
            ok = ErrMessage(G, errstr, "failed to compile expression");
          }
        }
        /* PBlockAndUnlockAPI() is not safe.
         * what if "v" is invalidated by another thread? */
        break;
    }
    switch (op->code) {
      case OMOP_ReferenceStore:
      case OMOP_ReferenceRecall:
      case OMOP_ReferenceValidate:
      case OMOP_ReferenceSwap:
        for (b = 0; b < I->NCSet; b++) {
          if (I->CSet[b]) {
            if ((b == op->i1) || (op->i1 < 0)) {
              int inv_flag = false;
              cs = I->CSet[b];
              if (cs && CoordSetValidateRefPos(cs)) {
                for (a = 0; a < I->NAtom; a++) {
                  s = I->AtomInfo[a].selEntry;
                  if (SelectorIsMember(G, s, sele)) {
                    ind = cs->atmToIdx(a);
                    if (ind >= 0) {
                      float* v = cs->coordPtr(ind);
                      RefPosType* rp = cs->RefPos + ind;
                      switch (op->code) {
                        case OMOP_ReferenceStore:
                          copy3f(v, rp->coord);
                          rp->specified = true;
                          break;
                        case OMOP_ReferenceRecall:
                          if (rp->specified) {
                            copy3f(rp->coord, v);
                            inv_flag = true;
                          }
                          break;
                        case OMOP_ReferenceValidate:
                          if (!rp->specified) {
                            copy3f(v, rp->coord);
                            rp->specified = true;
                          }
                          break;
                        case OMOP_ReferenceSwap:
                          if (rp->specified) {
                            copy3f(rp->coord, v1);
                            copy3f(v, rp->coord);
                            copy3f(v1, v);
                            inv_flag = true;
                          }
                          break;
                      }
                      op->i2++;
                    }
                  }
                }
              }
              if (inv_flag && cs) {
                cs->invalidateRep(cRepAll, cRepInvRep);
              }
            }
          }
        }
        break;
      case OMOP_AddHydrogens:
        if (ok) {
          ok &= ObjectMoleculeAddSeleHydrogensRefactored(I, sele, op->i1);
        }
        break;
      case OMOP_FixHydrogens:
        if (ok) ok &= ObjectMoleculeFixSeleHydrogens(I, sele, -1); /* state? */
        break;
      case OMOP_RevalenceFromSource:
      case OMOP_RevalenceByGuessing:
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          if (SelectorIsMember(G, ai->selEntry, sele)) {
            hit_flag = true;
            break;
          }
          ai++;
        }
        break;
      case OMOP_PrepareFromTemplate:
        ai0 = op->ai; /* template atom */
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            ai = I->AtomInfo + a;
            if (op->i1 != 3) {
              ai->hetatm = ai0->hetatm;
              ai->flags = ai0->flags;
              LexAssign(G, ai->chain, ai0->chain);
              strcpy(ai->alt, ai0->alt);
              LexAssign(G, ai->segi, ai0->segi);
            }
            if (op->i1 == 1) { /* mode 1, merge residue information */
              ai->inscode = ai0->inscode;
              ai->resv = ai0->resv;
              LexAssign(G, ai->resn, ai0->resn);
            }
            AtomInfoAssignColors(G, ai);
            if (op->i3) {
              if ((ai->elem[0] == ai0->elem[0]) &&
                  (ai->elem[1] == ai0->elem[1]))
                ai->color = ai0->color;
              else if (ai->protons == cAN_C) {
                ai->color = op->i4;
              }
            }
            ai->visRep = ai0->visRep;
            ai->id = -1;
            ai->rank = -1;
            op->i2++;
          }
        }
        break;
      case OMOP_Sort:
        for (a = 0; ok && a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            if (ok) ok &= ObjectMoleculeSort(I);
            break;
          }
        }
        break;
      case OMOP_SetAtomicSetting:
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          if (SelectorIsMember(G, ai->selEntry, sele)) {
            int uid = AtomInfoCheckUniqueID(G, ai);
            ai->has_setting = true;
            if (SettingUniqueSetTypedValue(G, uid, op->i1, op->i2, op->ii1))
              op->i4++;
          }
          ai++;
        }
        break;
      case OMOP_Pop:
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            if (SelectorMoveMember(G, s, sele, op->i1)) {
              op->i2--;
              op->i3++;
            }
            if (!op->i2) break;
          }
        }
        break;
      case OMOP_AVRT:     /* average vertex coordinate */
      case OMOP_StateVRT: /* state vertex coordinate */
      {
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int b_end = I->NCSet;
        if (op->code == OMOP_StateVRT && op->i1 < b_end) {
          b_end = op->i1 + 1;
        }
        if (op_i2) {
          use_matrices = SettingGet_i(I->G, I->Setting.get(), nullptr,
                                      cSetting_matrix_mode);
          if (use_matrices < 0) use_matrices = 0;
        }
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (!(priority = SelectorIsMember(G, s, sele))) continue;

          cnt = 0;

          // all states for AVRT, one state for StateVRT (don't use
          // StateIterator which depends on settings)
          for (b = op->i1; b < b_end; ++b) {
            if (!(cs = I->CSet[b])) continue;

            if ((a1 = cs->atmToIdx(a)) == -1) continue;

            if (!cnt) {
              VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
            }
            cnt++;
            vv2 = cs->coordPtr(a1);

            if (op_i2) { /* do we want transformed coordinates? */
              if (use_matrices) {
                if (!cs->Matrix.empty()) { /* state transformation */
                  transform44d3f(cs->Matrix.data(), vv2, v1);
                  vv2 = v1;
                }
              }
              if (obj_TTTFlag) {
                transformTTT44f3f(I->TTT, vv2, v1);
                vv2 = v1;
              }
            }

            // sum up coordinates over states
            vv1 = op->vv1 + (op->nvv1 * 3);
            add3f(vv1, vv2, vv1);
          }

          // number of summed up coordinates (states) for this atom
          VLACheck(op->vc1, int, op->nvv1);
          op->vc1[op->nvv1] = cnt;

          // ordered_selections
          if (op->vp1) {
            VLACheck(op->vp1, int, op->nvv1);
            op->vp1[op->nvv1] = priority;
          }

          // atom pointer VLA
          if (op->ai1VLA) {
            VLACheck(op->ai1VLA, AtomInfoType*, op->nvv1);
            op->ai1VLA[op->nvv1] = I->AtomInfo + a;
            I->AtomInfo[a].temp1 = a;
            /* KLUDGE ALERT!!! storing atom index in the temp1 field... */
          }

          // number of atoms in selection (incl. the ones with no coordinates)
          op->nvv1++;
        }
      } break;
      case OMOP_SFIT: /* state fitting within a single object */
        vt = pymol::malloc<float>(
            3 * op->nvv2); /* temporary (matching) target vertex pointers */
        cnt = 0;
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            cnt++;
            break;
          }
        }
        if (cnt) { /* only perform action for selected object */

          for (b = 0; b < I->NCSet; b++) {
            rms = -1.0;
            vt1 = vt; /* reset target vertex pointers */
            vt2 = op->vv2;
            t_i = 0; /* original target vertex index */
            if (I->CSet[b] && (b != op->i2)) {
              op->nvv1 = 0;
              for (a = 0; a < I->NAtom; a++) {
                s = I->AtomInfo[a].selEntry;
                if (SelectorIsMember(G, s, sele)) {
                  a1 = I->CSet[b]->atmToIdx(a);
                  if (a1 >= 0) {
                    match_flag = false;
                    while (t_i < op->nvv2) {
                      if (op->i1VLA[t_i] == a) { /* same atom? */
                        match_flag = true;
                        break;
                      }
                      if (op->i1VLA[t_i] < a) { /* catch up? */
                        t_i++;
                        vt2 += 3;
                      } else
                        break;
                    }
                    if (match_flag) {
                      VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                      vv2 = I->CSet[b]->coordPtr(a1);
                      vv1 = op->vv1 + (op->nvv1 * 3);
                      *(vv1++) = *(vv2++);
                      *(vv1++) = *(vv2++);
                      *(vv1++) = *(vv2++);
                      *(vt1++) = *(vt2);
                      *(vt1++) = *(vt2 + 1);
                      *(vt1++) = *(vt2 + 2);
                      op->nvv1++;
                    }
                  }
                }
              }
              if (op->nvv1 != op->nvv2) {
                PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: Missing atoms in state %d (%d instead of "
                "%d).\n",
                    b + 1, op->nvv1, op->nvv2 ENDFB(G);
              }
              if (op->nvv1) {
                if (op->i1 != 0) /* fitting flag */
                  rms = MatrixFitRMSTTTf(G, op->nvv1, op->vv1, vt, nullptr,
                                         op->ttt);
                else
                  rms = MatrixGetRMS(G, op->nvv1, op->vv1, vt, nullptr);
                if (op->i1 == 2) {
                  ObjectMoleculeTransformTTTf(I, op->ttt, b);

                  if (op->i3) {
                    const float divisor = (float)op->i3;
                    const float premult = (float)op->i3 - 1.0F;

                    /* mix flag is set, so average the prior target
                       coordinates with these coordinates */

                    vt2 = op->vv2;
                    t_i = 0; /* original target vertex index */
                    for (a = 0; a < I->NAtom; a++) {
                      s = I->AtomInfo[a].selEntry;
                      if (SelectorIsMember(G, s, sele)) {
                        a1 = I->CSet[b]->atmToIdx(a);
                        if (a1 >= 0) {
                          match_flag = false;
                          while (t_i < op->nvv2) {
                            if (op->i1VLA[t_i] == a) { /* same atom? */
                              match_flag = true;
                              break;
                            }
                            if (op->i1VLA[t_i] < a) { /* catch up? */
                              t_i++;
                              vt2 += 3;
                            } else
                              break;
                          }
                          if (match_flag) {
                            vv2 = I->CSet[b]->coordPtr(a1);
                            *(vt2) = ((premult * (*vt2)) + *(vv2++)) / divisor;
                            *(vt2 + 1) =
                                ((premult * (*(vt2 + 1))) + *(vv2++)) / divisor;
                            *(vt2 + 2) =
                                ((premult * (*(vt2 + 2))) + *(vv2++)) / divisor;
                          }
                        }
                      }
                    }
                  }
                }
              } else {
                PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: No matches found for state %d.\n",
                    b + 1 ENDFB(G);
              }
            }
            VLACheck(op->f1VLA, float, b);
            op->f1VLA[b] = rms;
          }
          VLASize(op->f1VLA, float,
                  I->NCSet); /* NOTE this action is object-specific! */
        }
        FreeP(vt);
        break;
      case OMOP_OnOff:
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            hit_flag = true;
            break;
          }
        }
        break;
      case OMOP_SaveUndo: /* save undo */
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            hit_flag = true;
            break;
          }
        }
        break;
      case OMOP_IdentifyObjects: /* identify atoms */
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            VLACheck(op->i1VLA, int, op->i1);
            op->i1VLA[op->i1] = I->AtomInfo[a].id;
            if (op->obj1VLA != nullptr) {
              VLACheck(op->obj1VLA, ObjectMolecule*, op->i1);
              op->obj1VLA[op->i1] = I;
            }
            op->i1++;
          }
        }
        break;
      case OMOP_Index: /* identify atoms */
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            VLACheck(op->i1VLA, int, op->i1);
            op->i1VLA[op->i1] = a; /* NOTE: need to incr by 1 before python */
            VLACheck(op->obj1VLA, ObjectMolecule*, op->i1);
            op->obj1VLA[op->i1] = I;
            op->i1++;
          }
        }
        break;
      case OMOP_GetObjects: /* identify atoms */
        for (a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if (SelectorIsMember(G, s, sele)) {
            VLACheck(op->obj1VLA, ObjectMolecule*, op->i1);
            op->obj1VLA[op->i1] = I;
            op->i1++;
            break;
          }
        }
        break;
      case OMOP_CountAtoms: /* count atoms in object, in selection */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) op->i1++;
          ai++;
        }
        break;
      case OMOP_PhiPsi:
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            VLACheck(op->i1VLA, int, op->i1);
            op->i1VLA[op->i1] = a;
            VLACheck(op->obj1VLA, ObjectMolecule*, op->i1);
            op->obj1VLA[op->i1] = I;
            VLACheck(op->f1VLA, float, op->i1);
            VLACheck(op->f2VLA, float, op->i1);
            if (ObjectMoleculeGetPhiPsi(I, a, op->f1VLA + op->i1,
                                        op->f2VLA + op->i1, op->i2))
              op->i1++;
          }
          ai++;
        }
        break;
      case OMOP_Cartoon: /* adjust cartoon type */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            if (ai->cartoon != op->i1) {
              op->i3++;
            }
            ai->cartoon = op->i1;
            op->i2++;
          }
          ai++;
        }
        break;
      case OMOP_Protect: /* protect atoms from movement */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            ai->protekted = op->i1;
            op->i2++;
          }
          ai++;
        }
        break;
      case OMOP_Mask: /* protect atoms from selection */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            ai->masked = op->i1;
            op->i2++;
          }
          ai++;
        }
        break;
      case OMOP_SetB: /* set B-value */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            ai->b = op->f1;
            op->i2++;
          }
          ai++;
        }
        break;
      case OMOP_Remove: /* flag atoms for deletion */
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          ai->deleteFlag = false;
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            ai->deleteFlag = true;
            op->i1++;
          }
          ai++;
        }
        break;
      case OMOP_GetChains:
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            // pointer hack
            ((std::set<lexidx_t>*)(void*)op->ii1)->insert(ai->chain);
            op->i1++;
          }
          ai++;
        }
        break;

      case OMOP_Spectrum:
        ai = I->AtomInfo.data();
        ai0 = nullptr;
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            skip_flag = false;
            if (op->i4 && ai0) /* byres and we've done a residue */
              if (AtomInfoSameResidue(G, ai, ai0)) skip_flag = true;
            if (!skip_flag) {
              c = (int)(0.49999 + op->i1 * (op->ff1[op->i3] - op->f1) / op->f2);
              if (c < 0) c = 0;
              if (c >= op->i1) c = op->i1 - 1;
              ai->color = op->ii1[c];

              /*               printf("%8.3 %8.3\n",ai->partial_charge, */
              if (op->i4) { /* byres */
                offset = -1;
                while ((a + offset) >= 0) {
                  ai0 = I->AtomInfo + a + offset;
                  if (AtomInfoSameResidue(G, ai, ai0)) {
                    ai0->color = op->ii1[c];
                    hit_flag = true;
                  } else
                    break;
                  offset--;
                }
                offset = 1;
                while ((a + offset) < I->NAtom) {
                  ai0 = I->AtomInfo + a + offset;
                  if (AtomInfoSameResidue(G, ai, ai0)) {
                    ai0->color = op->ii1[c];
                    hit_flag = true;
                  } else
                    break;
                  offset++;
                }
              }
              ai0 = ai;
            }
            op->i3++;
          }
          ai++;
        }
        break;

      case OMOP_SingleStateVertices: /* same as OMOP_VERT for a single state */
        ai = I->AtomInfo.data();
        if (op->cs1 < I->NCSet) {
          if (I->CSet[op->cs1]) {
            b = op->cs1;
            for (a = 0; a < I->NAtom; a++) {
              s = ai->selEntry;
              if (SelectorIsMember(G, s, sele)) {
                op->i1++;
                a1 = I->CSet[b]->atmToIdx(a);
                if (a1 >= 0) {
                  VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                  vv2 = I->CSet[b]->coordPtr(a1);
                  vv1 = op->vv1 + (op->nvv1 * 3);
                  *(vv1++) = *(vv2++);
                  *(vv1++) = *(vv2++);
                  *(vv1++) = *(vv2++);
                  op->nvv1++;
                }
              }
              ai++;
            }
          }
        }
        break;
      case OMOP_CSetIdxGetAndFlag:
        ai = I->AtomInfo.data();
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            for (b = op->cs1; b <= op->cs2; b++) {
              offset = b - op->cs1;
              if (b < I->NCSet) {
                if (I->CSet[b]) {
                  a1 = I->CSet[b]->atmToIdx(a);
                  if (a1 >= 0) {
                    op->ii1[op->i1 * offset + op->i2] = 1; /* presence flag */
                    vv1 = op->vv1 + 3 * (op->i1 * offset +
                                         op->i2); /* atom-based offset */
                    vv2 = I->CSet[b]->coordPtr(a1);
                    *(vv1++) = *(vv2++);
                    *(vv1++) = *(vv2++);
                    *(vv1++) = *(vv2++);
                    op->nvv1++;
                  }
                }
              }
            }
            op->i2++; /* atom index field for atoms within selection... */
          }
          ai++;
        }
        break;
      case OMOP_CSetIdxSetFlagged:
        ai = I->AtomInfo.data();
        hit_flag = false;
        for (a = 0; a < I->NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            for (b = op->cs1; b <= op->cs2; b++) {
              offset = b - op->cs1;
              if (b < I->NCSet) {
                if (I->CSet[b]) {
                  a1 = I->CSet[b]->atmToIdx(a);
                  if (a1 >= 0) {
                    if (op->ii1[op->i1 * offset + op->i2]) { /* copy flag */
                      vv1 = op->vv1 + 3 * (op->i1 * offset +
                                           op->i2); /* atom-based offset */
                      vv2 = I->CSet[b]->coordPtr(a1);
                      *(vv2++) = *(vv1++);
                      *(vv2++) = *(vv1++);
                      *(vv2++) = *(vv1++);
                      op->nvv1++;
                      hit_flag = true;
                    }
                  }
                }
              }
            }
            op->i2++; /* atom index field for atoms within selection... */
          }
          ai++;
        }
        break;
      case OMOP_SUMC: /* performance optimized to speed center & zoom actions */
      {
        /* given a selection, sum up all the coordinates (for centering) */
        float* op_v1 = op->v1;
        int op_i1 = op->i1;
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int i_NCSet = I->NCSet;
        int i_NAtom = I->NAtom;
        int i_DiscreteFlag = I->DiscreteFlag;
        CoordSet* const* i_CSet = I->CSet.data();
        if (op_i2) {
          use_matrices = SettingGet_i(I->G, I->Setting.get(), nullptr,
                                      cSetting_matrix_mode);
          if (use_matrices < 0) use_matrices = 0;
        }
        ai = I->AtomInfo.data();
        for (a = 0; a < i_NAtom; a++) {
          s = ai->selEntry;
          /* for each atom, if this current atom is in the selection */
          if (SelectorIsMember(G, s, sele)) {
            /* loop over all atoms; regardless of state */
            for (b = 0; b < i_NCSet; b++) {
              if (i_DiscreteFlag) {
                if ((cs = I->DiscreteCSet[a])) a1 = I->DiscreteAtmToIdx[a];
              } else {
                if ((cs = i_CSet[b])) a1 = cs->AtmToIdx[a];
              }
              /* if valid coordinate set and atom info for this atom */
              if (cs && (a1 >= 0)) {
                coord = cs->coordPtr(a1);
                if (op_i2) { /* do we want transformed coordinates? */
                  if (use_matrices) {
                    if (!cs->Matrix.empty()) { /* state transformation */
                      transform44d3f(cs->Matrix.data(), coord, v1);
                      coord = v1;
                    }
                  }
                  if (obj_TTTFlag) {
                    transformTTT44f3f(I->TTT, coord, v1);
                    coord = v1;
                  }
                }
                /* op_v1 += coord */
                add3f(op_v1, coord, op_v1);
                /* count += 1 */
                op_i1++;
              }
              if (i_DiscreteFlag) break;
            }
          }
          /* next atom */
          ai++;
        }
        op->i1 = op_i1;
      } break;
      case OMOP_MNMX: /* performance optimized to speed center & zoom actions */
      {
        float* op_v1 = op->v1;
        float* op_v2 = op->v2;
        int op_i1 = op->i1;
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int i_NCSet = I->NCSet;
        int i_NAtom = I->NAtom;
        int i_DiscreteFlag = I->DiscreteFlag;
        CoordSet* const* i_CSet = I->CSet.data();
        if (op_i2) {
          use_matrices = SettingGet_i(I->G, I->Setting.get(), nullptr,
                                      cSetting_matrix_mode);
          if (use_matrices < 0) use_matrices = 0;
        }
        ai = I->AtomInfo.data();
        for (a = 0; a < i_NAtom; a++) {
          s = ai->selEntry;
          if (SelectorIsMember(G, s, sele)) {
            for (b = 0; b < i_NCSet; b++) {
              if (i_DiscreteFlag) {
                if ((cs = I->DiscreteCSet[a])) a1 = I->DiscreteAtmToIdx[a];
              } else {
                if ((cs = i_CSet[b])) a1 = cs->AtmToIdx[a];
              }
              if (cs && (a1 >= 0)) {
                coord = cs->coordPtr(a1);
                if (op_i2) { /* do we want transformed coordinates? */
                  if (use_matrices) {
                    if (!cs->Matrix.empty()) { /* state transformation */
                      transform44d3f(cs->Matrix.data(), coord, v1);
                      coord = v1;
                    }
                  }
                  if (obj_TTTFlag) {
                    transformTTT44f3f(I->TTT, coord, v1);
                    coord = v1;
                  }
                }
                if (op_i1) {
                  if (op_v1[0] > coord[0]) op_v1[0] = coord[0];
                  if (op_v1[1] > coord[1]) op_v1[1] = coord[1];
                  if (op_v1[2] > coord[2]) op_v1[2] = coord[2];
                  if (op_v2[0] < coord[0]) op_v2[0] = coord[0];
                  if (op_v2[1] < coord[1]) op_v2[1] = coord[1];
                  if (op_v2[2] < coord[2]) op_v2[2] = coord[2];
                } else {
                  op_v1[0] = coord[0];
                  op_v1[1] = coord[1];
                  op_v1[2] = coord[2];
                  op_v2[0] = coord[0];
                  op_v2[1] = coord[1];
                  op_v2[2] = coord[2];
                }
                op_i1++;
              }
              if (i_DiscreteFlag) break;
            }
          }
          ai++;
        }
        op->i1 = op_i1;
      } break;
      default: {
        int inv_flag;

        switch (op->code) {
          case OMOP_INVA:
            /* set up an important optimization... */
            for (b = 0; b < I->NCSet; b++) {
              cs = I->CSet[b];
              if (cs) cs->objMolOpInvalidated = false;
            }
            break;
        }
        use_matrices =
            SettingGet_i(I->G, I->Setting.get(), nullptr, cSetting_matrix_mode);
        if (use_matrices < 0) use_matrices = 0;
        ai = I->AtomInfo.data();

        for (a = 0; a < I->NAtom; a++) {
          switch (op->code) {
            case OMOP_Flag:
              ai->flags &= op->i2; /* clear flag using mask */
              op->i4++;
              /* no break here is intentional!  */
            case OMOP_FlagSet:
            case OMOP_FlagClear:
            case OMOP_COLR: /* normal atom based loops */
            case OMOP_VISI:
            case OMOP_CheckVis:
            case OMOP_TTTF:
            case OMOP_ALTR:
            case OMOP_LABL:
            case OMOP_AlterState:
              s = ai->selEntry;
              if (SelectorIsMember(G, s, sele)) {
                switch (op->code) {
                  case OMOP_Flag:
                  case OMOP_FlagSet:
                    ai->flags |= op->i1; /* set flag */
                    op->i3++;
                    break;
                  case OMOP_FlagClear:
                    ai->flags &= op->i2; /* clear flag */
                    op->i3++;
                    break;
                  case OMOP_VISI:
                    switch (op->i2) {
                      case cVis_HIDE:
                        ai->visRep &= ~(op->i1);
                        I->visRep &= ~(op->i1);  // cell
                        break;
                      case cVis_SHOW:
                        ai->visRep |= op->i1 & cRepsAtomMask;
                        I->visRep |= op->i1 & cRepsObjectMask;  // cell
                        break;
                      case cVis_AS:
                        ai->visRep = op->i1 & cRepsAtomMask;
                        I->visRep = op->i1 & cRepsObjectMask;  // cell
                        break;
                    }
                    break;
                  case OMOP_CheckVis:
                    if ((ai->visRep & op->i1)) {
                      op->i2 = true;
                    }
                    break;
                  case OMOP_COLR:
                    if (op->i1 == cColorAtomic)
                      ai->color = AtomInfoGetColor(G, ai);
                    else
                      ai->color = op->i1;
                    hit_flag = true;
                    op->i2++;
                    break;
                  case OMOP_TTTF:
                    hit_flag = true;
                    break;
                  case OMOP_LABL:
                    if (ok) {
                      if (!op->s1[0]) {
                        if (ai->label) {
                          op->i1--; /* negative if unlabelling */
                          LexAssign(G, ai->label, 0);
                        }
                        ai->visRep &= ~cRepLabelBit;
                        hit_flag = true;
                      } else {
                        switch (op->i2) {
                          case cExecutiveLabelEvalOn: {
                            /* python label expression evaluation */
                            CoordSet* cs = nullptr;
                            if (I->DiscreteFlag && I->DiscreteCSet) {
                              cs = I->DiscreteCSet[a];
                            } else if (I->NCSet == 1) {
                              cs = I->CSet[0];
                            }
                            if (PLabelAtom(I->G, I, cs, expr_co, a)) {
                              if (ai->label) {
                                op->i1++; /* only if the string has been set,
                                             report labelled */
                              }
                              ai->visRep |= cRepLabelBit;
                              hit_flag = true;
                            } else {
                              ok = false;
                            }
                          } break;
                          case cExecutiveLabelEvalAlt: {
                            if (PLabelAtomAlt(I->G, &I->AtomInfo[a], I->Name,
                                              op->s1, a)) {
                              if (ai->label) {
                                op->i1++; /* only if the string has been set,
                                             report labelled */
                              }
                              ai->visRep |= cRepLabelBit;
                              hit_flag = true;
                            } else {
                              ok = false;
                            }
                          } break;
                          case cExecutiveLabelEvalOff: {
                            /* simple string label text */
                            AtomInfoType* ai = I->AtomInfo + a;
                            LexDec(G, ai->label);
                            ai->label = LexIdx(G, op->s1);
                          } break;
                        }
                      }
                    }
                    break;
                  case OMOP_ALTR:
                    if (ok) {
                      CoordSet* cs = nullptr;
                      if (I->DiscreteFlag && I->DiscreteCSet) {
                        cs = I->DiscreteCSet[a];
                      } else if (I->NCSet == 1) {
                        cs = I->CSet[0];
                      }
                      if (PAlterAtom(I->G, I, cs, expr_co, op->i2, a,
                                     op->py_ob1))
                        op->i1++;
                      else
                        ok = false;
                    }
                    break;
                  case OMOP_AlterState:
                    if (ok) {
                      if (op->i2 < I->NCSet) {
                        cs = I->CSet[op->i2];
                        if (cs) {
                          a1 = cs->atmToIdx(a);
                          if (a1 >= 0) {
                            if (PAlterAtomState(I->G, expr_co, op->i3, I, cs, a,
                                                a1, op->i2, op->py_ob1)) {
                              op->i1++;
                              hit_flag = true;
                            } else
                              ok = false;
                          }
                        }
                      }
                    }
                    break;
                }
                break;
              }
              break;

              /* coord-set based properties, iterating only a single coordinate
               * set */
            case OMOP_CSetMinMax:
            case OMOP_CSetCameraMinMax:
            case OMOP_CSetMaxDistToPt:
            case OMOP_CSetSumSqDistToPt:
            case OMOP_CSetSumVertices:
            case OMOP_CSetMoment:
              cs = nullptr;
              if ((op->cs1 >= 0) && (op->cs1 < I->NCSet)) {
                cs = I->CSet[op->cs1];
              } else if (op->include_static_singletons) {
                if ((I->NCSet == 1) &&
                    (SettingGet_b(G, nullptr, I->Setting.get(),
                                  cSetting_static_singletons))) {
                  cs = I->CSet[0]; /*treat static singletons as present in each
                                      state */
                }
              }

              if (cs) {
                s = ai->selEntry;
                if (SelectorIsMember(G, s, sele)) {
                  switch (op->code) {
                    case OMOP_CSetSumVertices:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        coord = cs->coordPtr(a1);
                        if (op->i2) { /* do we want transformed coordinates? */
                          if (use_matrices) {
                            if (!cs->Matrix
                                     .empty()) { /* state transformation */
                              transform44d3f(cs->Matrix.data(), coord, v1);
                              coord = v1;
                            }
                          }
                          if (I->TTTFlag) {
                            transformTTT44f3f(I->TTT, coord, v1);
                            coord = v1;
                          }
                        }
                        add3f(op->v1, coord, op->v1);
                        op->i1++;
                      }
                      break;
                    case OMOP_CSetMinMax:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        coord = cs->coordPtr(a1);
                        if (op->i2) { /* do we want transformed coordinates? */
                          if (use_matrices) {
                            if (!cs->Matrix
                                     .empty()) { /* state transformation */
                              transform44d3f(cs->Matrix.data(), coord, v1);
                              coord = v1;
                            }
                          }
                          if (I->TTTFlag) {
                            transformTTT44f3f(I->TTT, coord, v1);
                            coord = v1;
                          }
                        }
                        if (op->i1) {
                          for (c = 0; c < 3; c++) {
                            if (*(op->v1 + c) > *(coord + c))
                              *(op->v1 + c) = *(coord + c);
                            if (*(op->v2 + c) < *(coord + c))
                              *(op->v2 + c) = *(coord + c);
                          }
                        } else {
                          for (c = 0; c < 3; c++) {
                            *(op->v1 + c) = *(coord + c);
                            *(op->v2 + c) = *(coord + c);
                          }
                        }
                        op->i1++;
                      }
                      break;
                    case OMOP_CSetCameraMinMax:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        coord = cs->coordPtr(a1);
                        if (op->i2) { /* do we want transformed coordinates? */
                          if (use_matrices) {
                            if (!cs->Matrix
                                     .empty()) { /* state transformation */
                              transform44d3f(cs->Matrix.data(), coord, v1);
                              coord = v1;
                            }
                          }
                          if (I->TTTFlag) {
                            transformTTT44f3f(I->TTT, coord, v1);
                            coord = v1;
                          }
                        }
                        MatrixTransformC44fAs33f3f(op->mat1, coord, v1);
                        /* convert to view-space */
                        coord = v1;
                        if (op->i1) {
                          for (c = 0; c < 3; c++) {
                            if (*(op->v1 + c) > *(coord + c))
                              *(op->v1 + c) = *(coord + c);
                            if (*(op->v2 + c) < *(coord + c))
                              *(op->v2 + c) = *(coord + c);
                          }
                        } else {
                          for (c = 0; c < 3; c++) {
                            *(op->v1 + c) = *(coord + c);
                            *(op->v2 + c) = *(coord + c);
                          }
                        }
                        op->i1++;
                      }
                      break;
                    case OMOP_CSetSumSqDistToPt:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        float dist;
                        coord = cs->coordPtr(a1);
                        if (op->i2) { /* do we want transformed coordinates? */
                          if (use_matrices) {
                            if (!cs->Matrix
                                     .empty()) { /* state transformation */
                              transform44d3f(cs->Matrix.data(), coord, v1);
                              coord = v1;
                            }
                          }
                          if (I->TTTFlag) {
                            transformTTT44f3f(I->TTT, coord, v1);
                            coord = v1;
                          }
                        }
                        dist = (float)diff3f(op->v1, coord);
                        op->d1 += dist * dist;
                        op->i1++;
                      }
                      break;
                    case OMOP_CSetMaxDistToPt:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        float dist;
                        coord = cs->coordPtr(a1);
                        if (op->i2) { /* do we want transformed coordinates? */
                          if (use_matrices) {
                            if (!cs->Matrix
                                     .empty()) { /* state transformation */
                              transform44d3f(cs->Matrix.data(), coord, v1);
                              coord = v1;
                            }
                          }
                          if (I->TTTFlag) {
                            transformTTT44f3f(I->TTT, coord, v1);
                            coord = v1;
                          }
                        }
                        dist = (float)diff3f(op->v1, coord);
                        if (dist > op->f1) op->f1 = dist;
                        op->i1++;
                      }
                      break;
                    case OMOP_CSetMoment:
                      a1 = cs->atmToIdx(a);
                      if (a1 >= 0) {
                        subtract3f(cs->coordPtr(a1), op->v1, v1);
                        v2 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
                        op->d[0][0] += v2 - v1[0] * v1[0];
                        op->d[0][1] += -v1[0] * v1[1];
                        op->d[0][2] += -v1[0] * v1[2];
                        op->d[1][0] += -v1[1] * v1[0];
                        op->d[1][1] += v2 - v1[1] * v1[1];
                        op->d[1][2] += -v1[1] * v1[2];
                        op->d[2][0] += -v1[2] * v1[0];
                        op->d[2][1] += -v1[2] * v1[1];
                        op->d[2][2] += v2 - v1[2] * v1[2];
                      }
                      break;
                  }
                }
              }
              break;
            default:
              /* coord-set based properties, iterating as all coordsets within
               * atoms */
              for (b = 0; b < I->NCSet; b++) {
                if (I->DiscreteFlag) {
                  cs = I->DiscreteCSet[a];
                } else {
                  cs = I->CSet[b];
                }
                if (cs) {
                  s = ai->selEntry;
                  inv_flag = false;
                  if (SelectorIsMember(G, s, sele)) {
                    switch (op->code) {
                      case OMOP_CameraMinMax:
                        a1 = cs->atmToIdx(a);
                        if (a1 >= 0) {
                          coord = cs->coordPtr(a1);
                          if (op->i2) { /* do we want transformed coordinates?
                                         */
                            if (use_matrices) {
                              if (!cs->Matrix
                                       .empty()) { /* state transformation */
                                transform44d3f(cs->Matrix.data(), coord, v1);
                                coord = v1;
                              }
                            }
                            if (I->TTTFlag) {
                              transformTTT44f3f(I->TTT, coord, v1);
                              coord = v1;
                            }
                          }
                          MatrixTransformC44fAs33f3f(op->mat1, coord, v1);
                          /* convert to view-space */
                          coord = v1;
                          if (op->i1) {
                            for (c = 0; c < 3; c++) {
                              if (*(op->v1 + c) > *(coord + c))
                                *(op->v1 + c) = *(coord + c);
                              if (*(op->v2 + c) < *(coord + c))
                                *(op->v2 + c) = *(coord + c);
                            }
                          } else {
                            for (c = 0; c < 3; c++) {
                              *(op->v1 + c) = *(coord + c);
                              *(op->v2 + c) = *(coord + c);
                            }
                          }
                          op->i1++;
                        }
                        break;
                      case OMOP_MaxDistToPt:
                        a1 = cs->atmToIdx(a);
                        if (a1 >= 0) {
                          float dist;
                          coord = cs->coordPtr(a1);
                          if (op->i2) { /* do we want transformed coordinates?
                                         */
                            if (use_matrices) {
                              if (!cs->Matrix
                                       .empty()) { /* state transformation */
                                transform44d3f(cs->Matrix.data(), coord, v1);
                                coord = v1;
                              }
                            }
                            if (I->TTTFlag) {
                              transformTTT44f3f(I->TTT, coord, v1);
                              coord = v1;
                            }
                          }
                          dist = (float)diff3f(coord, op->v1);
                          if (dist > op->f1) op->f1 = dist;
                          op->i1++;
                        }
                        break;
                      case OMOP_INVA:
                        if (!cs->objMolOpInvalidated) {
                          /* performance optimization: avoid repeatedly
                             calling invalidate on the same coord sets */
                          a1 = cs->atmToIdx(a);
                          if (a1 >=
                              0) /* selection touches this coordinate set */
                            inv_flag = true; /* so set the invalidation flag */
                        }
                        break;
                      case OMOP_VERT:
                        /* get the atom index whether it's discrete or not */
                        a1 = cs->atmToIdx(a);
                        if (a1 >= 0) {
                          /* if a1 is a valid atom index, then copy it's xyz
                           * coordinates into vv1; increment the counter, nvv1
                           */
                          VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                          vv2 = cs->coordPtr(a1);
                          vv1 = op->vv1 + (op->nvv1 * 3);
                          *(vv1++) = *(vv2++);
                          *(vv1++) = *(vv2++);
                          *(vv1++) = *(vv2++);
                          op->nvv1++;
                        }
                        break;
                      case OMOP_SVRT:
                        /* gives us only vertices for a specific coordinate set
                         */
                        if (b == op->i1) {
                          a1 = cs->atmToIdx(a);
                          if (a1 >= 0) {
                            VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                            VLACheck(op->i1VLA, int, op->nvv1);
                            op->i1VLA[op->nvv1] =
                                a; /* save atom index for later comparisons */
                            vv2 = cs->coordPtr(a1);
                            vv1 = op->vv1 + (op->nvv1 * 3);
                            *(vv1++) = *(vv2++);
                            *(vv1++) = *(vv2++);
                            *(vv1++) = *(vv2++);
                            op->nvv1++;
                          }
                        }
                        break;
                      case OMOP_MOME:
                        /* Moment of inertia tensor - unweighted - assumes v1 is
                         * center of molecule */
                        a1 = cs->atmToIdx(a);
                        if (a1 >= 0) {
                          subtract3f(cs->coordPtr(a1), op->v1, v1);
                          v2 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
                          op->d[0][0] += v2 - v1[0] * v1[0];
                          op->d[0][1] += -v1[0] * v1[1];
                          op->d[0][2] += -v1[0] * v1[2];
                          op->d[1][0] += -v1[1] * v1[0];
                          op->d[1][1] += v2 - v1[1] * v1[1];
                          op->d[1][2] += -v1[1] * v1[2];
                          op->d[2][0] += -v1[2] * v1[0];
                          op->d[2][1] += -v1[2] * v1[1];
                          op->d[2][2] += v2 - v1[2] * v1[2];
                        }
                        break;
                    }
                  }
                  switch (op->code) {
                      /* full coord-set based */
                    case OMOP_INVA:
                      /* shouldn't this be calling the object invalidation
                       * routine instead? */
                      if (inv_flag) {
                        cs->objMolOpInvalidated = true;
                        for (auto d = cRep_t(0); d < cRepCnt; ++d) {
                          if ((1 << d) & op->i1) {
                            cs->invalidateRep(d, cRepInv_t(op->i2));
                          }
                        }
                      }
                      break;
                  }
                  if (I->DiscreteFlag) {
                    /* don't iterate every coordinate set for discrete objects!
                     */
                    break;
                  }
                }
              } /* end coordset section */
              break;
          }
          ai++;
        }
      } /* case: default */
      break;
    }
    if (hit_flag) {
      switch (op->code) {
        case OMOP_COLR:
          ExecutiveUpdateColorDepends(I->G, I);
          break;
        case OMOP_TTTF:
          ObjectMoleculeTransformTTTf(I, op->ttt, -1);
          break;
        case OMOP_LABL:
          I->invalidate(cRepLabel, cRepInvText, -1);
          break;
        case OMOP_AlterState: /* overly coarse - doing all states, could do just
                                 1 */
          if (!op->i3) {      /* not read_only? */
            I->invalidate(cRepAll, cRepInvRep, -1);
            SceneChanged(G);
          }
          break;
        case OMOP_CSetIdxSetFlagged:
          I->invalidate(cRepAll, cRepInvRep, -1);
          SceneChanged(G);
          break;
        case OMOP_SaveUndo:
          op->i2 = true;
          ObjectMoleculeSaveUndo(I, op->i1, false);
          break;
        case OMOP_OnOff:
          ExecutiveSetObjVisib(G, I->Name, op->i1, false);
          break;
        case OMOP_RevalenceFromSource:
          if (ObjectMoleculeXferValences(I, op->i1, op->i2, op->i3, op->obj3,
                                         op->i4, op->i5, op->i6)) {
            ObjectMoleculeVerifyChemistry(I, op->i3);
            I->invalidate(cRepAll, cRepInvBonds, op->i3);
          }
          break;
        case OMOP_RevalenceByGuessing: {
          int* flag1 = pymol::calloc<int>(I->NAtom);
          int* flag2 = pymol::calloc<int>(I->NAtom);
          if (flag1 && flag2) {
            int a;
            int* f1 = flag1;
            int* f2 = flag2;
            const AtomInfoType* ai = I->AtomInfo.data();
            for (a = 0; a < I->NAtom; a++) {
              *(f1++) = SelectorIsMember(G, ai->selEntry, op->i1);
              *(f2++) = SelectorIsMember(G, ai->selEntry, op->i2);
              ai++;
            }
            {
              int target_state = op->i3;
              if (target_state < 0) target_state = 0; /* TO DO */
              ObjectMoleculeGuessValences(I, target_state, flag1, flag2,
                                          op->i4);
              ObjectMoleculeVerifyChemistry(I, target_state);
              I->invalidate(cRepAll, cRepInvBonds, target_state);
            }
            FreeP(flag1);
            FreeP(flag2);
          }
        } break;
      }
    }

    /* always run on exit... */
    switch (op->code) {
      case OMOP_LABL:
        if (op->i2 != cExecutiveLabelEvalOn) {
          break;
        }
      case OMOP_ALTR:
      case OMOP_AlterState:
        Py_XDECREF(expr_co);
        break;
    }
    /* */
  }

  return ok;
}
