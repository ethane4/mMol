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

// Object matrices and transformations.

#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

static int ExecutiveSetObjectMatrix2(PyMOLGlobals* G, pymol::CObject* obj,
                                     int state, double* matrix);
static int ExecutiveGetObjectMatrix2(PyMOLGlobals* G, pymol::CObject* obj,
                                     int state, double** matrix, int incl_ttt);
static pymol::Result<> ExecutiveTransformObjectSelection2(
    PyMOLGlobals* G, pymol::CObject* obj, int state, const char* s1, int log,
    const float* matrix, int homogenous, int global);

int ExecutiveMatrixCopy2(PyMOLGlobals* G, pymol::CObject* source_obj,
                         pymol::CObject* target_obj, int source_mode,
                         int target_mode, int source_state, int target_state,
                         int target_undo, int log, int quiet) {
  /*  mode 0: raw coordinates, as per the txf history
     mode 1: object TTT matrix
     mode 2: state matrix */

  int ok = true;
  int copy_ttt_too = false;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if (matrix_mode < 0) matrix_mode = 0; /* for now */

  if ((source_mode < 0) && (target_mode < 0)) {
    copy_ttt_too = true;
  }
  if (source_mode < 0) source_mode = matrix_mode;
  if (target_mode < 0) target_mode = matrix_mode;

  switch (source_mode) {
    case 0: /* txf history is the source matrix */
    {
      double* history = nullptr;
      int found = ExecutiveGetObjectMatrix2(G, source_obj, source_state,
                                            &history, false);
      if (found) {
        switch (target_mode) {
          case 0: /* apply changes to coordinates in the target object */
          {
            double temp_inverse[16];
            if (target_undo) {
              double* target_history = nullptr;
              int target_found = ExecutiveGetObjectMatrix2(
                  G, source_obj, target_state, &target_history, false);
              if (target_found && target_history) {
                invert_special44d44d(target_history, temp_inverse);
                if (history) {
                  right_multiply44d44d(temp_inverse, history);
                  history = temp_inverse;
                } else {
                  history = temp_inverse;
                }
              }
              {
                float historyf[16];
                if (history) {
                  convert44d44f(history, historyf);
                } else {
                  identity44f(historyf);
                }
                ExecutiveTransformObjectSelection2(G, target_obj, target_state,
                                                   "", log, historyf, true,
                                                   false);
              }
            }
            if (copy_ttt_too) {
              const float* tttf;
              int found = ObjectGetTTT(source_obj, &tttf, -1);
              if (found) {
                ObjectSetTTT(target_obj, tttf, -1, -1);
                target_obj->invalidate(cRepNone, cRepInvExtents, -1);
              }
            }
          } break;
          case 1: /* applying changes to the object's TTT matrix */
            if (history) {
              float tttf[16];
              convertR44dTTTf(history, tttf);
              ObjectSetTTT(target_obj, tttf, -1, -1);
            } else {
              ObjectSetTTT(target_obj, nullptr, -1, -1);
            }
            target_obj->invalidate(cRepNone, cRepInvExtents, -1);
            break;
          case 2: /* applying changes to the state matrix */
            ok =
                ExecutiveSetObjectMatrix2(G, target_obj, target_state, history);
            break;
        }
        break;
      }
    } break;
    case 1: /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */

      const float* tttf;
      int found = ObjectGetTTT(source_obj, &tttf, -1);
      if (found) {
        switch (target_mode) {
          case 0: /* coordinates & history unsupported.. */
            /* should complain */
            break;
          case 1: /* TTT */
            ObjectSetTTT(target_obj, tttf, -1, -1);
            target_obj->invalidate(cRepNone, cRepInvExtents, -1);
            break;
          case 2: /* State */
            if (tttf) {
              double homo[16];
              convertTTTfR44d(tttf, homo);
              ok = ExecutiveSetObjectMatrix2(G, target_obj, -1, homo);
            } else {
              ok = ExecutiveSetObjectMatrix2(G, target_obj, -1, nullptr);
            }
            break;
        }
      }
    } break;
    case 2: /* from the state matrix */
    {
      double* homo;
      int found =
          ExecutiveGetObjectMatrix2(G, source_obj, source_state, &homo, false);
      if (found) {
        switch (target_mode) {
          case 0: /* coordinates & history */
            /* TODO */
            break;
          case 1: /* TTT */
            if (homo) {
              float tttf[16];
              convertR44dTTTf(homo, tttf);
              ObjectSetTTT(target_obj, tttf, -1, -1);
              target_obj->invalidate(cRepNone, cRepInvExtents, -1);
            } else {
              ObjectSetTTT(target_obj, nullptr, -1, -1);
              target_obj->invalidate(cRepNone, cRepInvExtents, -1);
            }
            break;
          case 2: /* State */
            ok = ExecutiveSetObjectMatrix2(G, target_obj, target_state, homo);
            if (copy_ttt_too) {
              const float* tttf;
              int found = ObjectGetTTT(source_obj, &tttf, -1);
              if (found) {
                ObjectSetTTT(target_obj, tttf, -1, -1);
                target_obj->invalidate(cRepNone, cRepInvExtents, -1);
              }
            }
            break;
        }
      }
    } break;
  }
  SceneInvalidate(G);
  return ok;
}

int ExecutiveMatrixCopy(PyMOLGlobals* G, const char* source_name,
                        const char* target_name, int source_mode,
                        int target_mode, int source_state, int target_state,
                        int target_undo, int log, int quiet) {
  /*  mode 0: raw coordinates, as per the txf history
     mode 1: object TTT matrix
     mode 2: state matrix
     mode 3 (source only): camera matrix transformation */

  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  SpecRec* src_rec = nullptr;
  int ok = true;
  int copy_ttt_too = false;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if (matrix_mode < 0) matrix_mode = 0; /* for now */

  if ((source_mode < 0) && (target_mode < 0)) {
    copy_ttt_too = true;
  }

  if (source_mode < 0) source_mode = matrix_mode;
  if (target_mode < 0) target_mode = matrix_mode;
  if (source_name[0] == 0) {
    source_mode = 3;
    target_undo = 0;
  } else
    src_rec = ExecutiveFindSpec(G, source_name);

  if (source_mode != 3 && !src_rec) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    " %s-Warning: Can't find source object '%s'.\n", __FUNCTION__,
        source_name ENDFB(G);
  }

  int list_id = ExecutiveGetNamesListFromPattern(G, target_name, true,
                                                 cExecExpandKeepGroups);

  if (!list_id) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    " %s-Warning: No match for target '%s'.\n", __FUNCTION__,
        target_name ENDFB(G);
  }

  switch (source_mode) {
    case 0: /* txf history is the source matrix */
    {
      double* src_history = nullptr;
      int found = ExecutiveGetObjectMatrix(G, source_name, source_state,
                                           &src_history, false);
      if (found) {
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec* rec;

        while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                         (TrackerRef**)(void*)&rec)) {
          if (rec && (rec != src_rec)) {
            switch (rec->type) {
              case cExecObject:

                switch (target_mode) {
                  case 0: /* apply changes to coordinates in the target object
                           */
                  {
                    /* Use a per-target copy of the source history so that
                       modifications for one target don't corrupt the matrix
                       for subsequent targets (avoids pointer aliasing bug) */
                    double history_buf[16];
                    double* history = nullptr;
                    if (src_history) {
                      copy44d(src_history, history_buf);
                      history = history_buf;
                    }
                    double temp_inverse[16];
                    if (target_undo) {
                      double* target_history = nullptr;
                      int target_found = ExecutiveGetObjectMatrix(
                          G, rec->name, target_state, &target_history, false);
                      if (target_found && target_history) {
                        invert_special44d44d(target_history, temp_inverse);
                        if (history) {
                          right_multiply44d44d(temp_inverse, history);
                          history = temp_inverse;
                        } else {
                          history = temp_inverse;
                        }
                      }
                    }
                    {
                      float historyf[16];
                      if (history) {
                        convert44d44f(history, historyf);
                      } else {
                        identity44f(historyf);
                      }
                      ExecutiveTransformObjectSelection(G, rec->name,
                                                        target_state, "", log,
                                                        historyf, true, false);
                    }
                    if (copy_ttt_too) {
                      const float* tttf;
                      int found = ExecutiveGetObjectTTT(G, source_name, &tttf,
                                                        -1, quiet);
                      if (found) {
                        ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet,
                                              -1);
                      }
                    }
                  } break;
                  case 1: /* applying changes to the object's TTT matrix */
                    if (src_history) {
                      float tttf[16];
                      convertR44dTTTf(src_history, tttf);
                      ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                    } else {
                      ExecutiveSetObjectTTT(G, rec->name, nullptr, -1, quiet,
                                            -1);
                    }
                    /* to do: logging, return values, etc. */
                    break;
                  case 2: /* applying changes to the state matrix */
                    ok = ExecutiveSetObjectMatrix(G, rec->name, target_state,
                                                  src_history);
                    break;
                }
                break;
            }
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    } break;
    case 1: /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */

      const float* tttf;
      int found = ExecutiveGetObjectTTT(G, source_name, &tttf, -1, quiet);
      if (found) {
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec* rec;

        while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                         (TrackerRef**)(void*)&rec)) {
          if (rec && (rec != src_rec)) {
            switch (rec->type) {
              case cExecObject:

                switch (target_mode) {
                  case 0: /* coordinates & history unsupported.. */
                    /* should complain */
                    break;
                  case 1: /* TTT */
                    ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                    break;
                  case 2: /* State */
                    if (tttf) {
                      double homo[16];
                      convertTTTfR44d(tttf, homo);
                      ok = ExecutiveSetObjectMatrix(G, rec->name, -1, homo);
                    } else {
                      ok = ExecutiveSetObjectMatrix(G, rec->name, -1, nullptr);
                    }
                    break;
                }
            }
          }
        }

        TrackerDelIter(I_Tracker, iter_id);
      }
    } break;
    case 2: /* from the state matrix */
    {
      double* homo;
      int found =
          ExecutiveGetObjectMatrix(G, source_name, source_state, &homo, false);
      if (found) {
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec* rec;

        while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                         (TrackerRef**)(void*)&rec)) {
          if (rec && (rec != src_rec)) {
            switch (rec->type) {
              case cExecObject:

                switch (target_mode) {
                  case 0: /* coordinates & history */
                    /* TODO */
                    break;
                  case 1: /* TTT */
                    if (homo) {
                      float tttf[16];
                      convertR44dTTTf(homo, tttf);
                      ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                    } else {
                      ExecutiveSetObjectTTT(G, rec->name, nullptr, -1, quiet,
                                            -1);
                    }
                    break;
                  case 2: /* State */
                    ok = ExecutiveSetObjectMatrix(G, rec->name, target_state,
                                                  homo);
                    if (copy_ttt_too) {
                      const float* tttf;
                      int found = ExecutiveGetObjectTTT(G, source_name, &tttf,
                                                        -1, quiet);
                      if (found) {
                        ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet,
                                              -1);
                      }
                    }
                    break;
                }
            }
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    } break;
    case 3: /* camera */
    {
      SceneViewType view;
      double homo[16], *history;
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec* rec;
      SceneGetView(G, view);
      homo[0] = view[0];
      homo[1] = view[4];
      homo[2] = view[8];
      homo[3] = -(view[0] * view[19] + view[4] * view[20] + view[8] * view[21]);
      homo[4] = view[1];
      homo[5] = view[5];
      homo[6] = view[9];
      homo[7] = -(view[1] * view[19] + view[5] * view[20] + view[9] * view[21]);
      homo[8] = view[2];
      homo[9] = view[6];
      homo[10] = view[10];
      homo[11] =
          -(view[2] * view[19] + view[6] * view[20] + view[10] * view[21]);
      homo[12] = 0.0;
      homo[13] = 0.0;
      homo[14] = 0.0;
      homo[15] = 1.0;
      history = homo;
      while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                       (TrackerRef**)(void*)&rec)) {
        if (rec && (rec != src_rec)) {
          switch (rec->type) {
            case cExecObject:

              switch (target_mode) {
                case 0: /* apply changes to coordinates in the target object */
                {
                  double temp_inverse[16];
                  if (target_undo) {
                    double* target_history = nullptr;
                    int target_found = ExecutiveGetObjectMatrix(
                        G, rec->name, target_state, &target_history, false);
                    if (target_found && target_history) {
                      invert_special44d44d(target_history, temp_inverse);
                      if (history) {
                        right_multiply44d44d(temp_inverse, history);
                        history = temp_inverse;
                      } else {
                        history = temp_inverse;
                      }
                    }
                  }
                  {
                    float historyf[16];
                    if (history) {
                      convert44d44f(history, historyf);
                    } else {
                      identity44f(historyf);
                    }
                    ExecutiveTransformObjectSelection(G, rec->name,
                                                      target_state, "", log,
                                                      historyf, true, false);
                  }
                } break;
                case 1: /* applying changes to the object's TTT matrix */
                  if (history) {
                    float tttf[16];
                    convertR44dTTTf(history, tttf);
                    ExecutiveSetObjectTTT(G, rec->name, tttf, -1, quiet, -1);
                  } else {
                    ExecutiveSetObjectTTT(G, rec->name, nullptr, -1, quiet, -1);
                  }
                  /* to do: logging, return values, etc. */
                  break;
                case 2: /* applying changes to the state matrix */
                  ok = ExecutiveSetObjectMatrix(G, rec->name, target_state,
                                                history);
                  break;
              }
              break;
          }
        }
        TrackerDelIter(I_Tracker, iter_id);
      }
    } break;
  }

  TrackerDelList(I_Tracker, list_id);

  SceneInvalidate(G);
  return ok;
}

pymol::Result<> ExecutiveResetMatrix(PyMOLGlobals* G, const char* name,
                                     int mode, int state, int log, int quiet) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;
  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if (matrix_mode < 0) matrix_mode = 0; /* for now */

  if (mode < 0) mode = matrix_mode;

  bool obj_found = false;
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec && (rec->type == cExecObject)) {
      /*  pymol::CObject *obj = ExecutiveFindObjectByName(G,name); */
      pymol::CObject* obj = rec->obj;
      if (obj) {
        obj_found = true;
        switch (obj->type) {
          case cObjectMolecule:
            switch (mode) {
              case 0: /* transformations already applied to the coordinates */
                for (StateIterator iter(rec->obj, state); iter.next();) {
                  double* history = nullptr;
                  bool found = ExecutiveGetObjectMatrix2(
                      G, rec->obj, iter.state, &history, false);
                  if (found && history) {
                    double temp_inverse[16];
                    float historyf[16];
                    invert_special44d44d(history, temp_inverse);
                    convert44d44f(temp_inverse, historyf);
                    ExecutiveTransformObjectSelection2(G, rec->obj, iter.state,
                                                       "", log, historyf, true,
                                                       false);
                  }
                }
                break;
              case 1: /* operate on the TTT display matrix */
                ObjectResetTTT(
                    obj, SettingGetGlobal_b(G, cSetting_movie_auto_store));
                obj->invalidate(cRepNone, cRepInvExtents, -1);

                break;
              case 2: /* applying changes to the state matrix */
              {
                double ident[16];
                identity44d(ident);
                ExecutiveSetObjectMatrix(G, rec->name, state, ident);
              } break;
            }
            break;
          default:
            if (auto* objstate = obj->getObjectState(state)) {
              ObjectStateResetMatrix(objstate);
              obj->invalidate(cRepNone, cRepInvExtents, state);
            }
            break;
        }
      }
    }
  }
  if (!obj_found) {
    return pymol::make_error("No object found");
  }
  return {};
}

static double ret_mat[16]; /* UGH ..not thread-safe */

static int ExecutiveGetObjectMatrix2(PyMOLGlobals* G, pymol::CObject* obj,
                                     int state, double** matrix, int incl_ttt) {
  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */

  int ok = false;
  if (state < 0) {
    /* to do -- TTT only */
  } else {
    if (auto* objstate = obj->getObjectState(state)) {
      *matrix = ObjectStateGetMatrix(objstate);
      ok = true;
    }

    if (ok && incl_ttt) {
      const float* ttt;
      double tttd[16];
      if (ObjectGetTTT(obj, &ttt, -1)) {
        convertTTTfR44d(ttt, tttd);
        if (*matrix) {
          copy44d(*matrix, ret_mat);
        } else {
          identity44d(ret_mat);
        }
        left_multiply44d44d(tttd, ret_mat);
        *matrix = ret_mat;
      }
    }
  }
  return ok;
}

int ExecutiveGetObjectMatrix(PyMOLGlobals* G, const char* name, int state,
                             double** matrix, int incl_ttt) {
  int ok = false;
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    return ExecutiveGetObjectMatrix2(G, obj, state, matrix, incl_ttt);
  }
  return ok;
}

static int ExecutiveSetObjectMatrix2(PyMOLGlobals* G, pymol::CObject* obj,
                                     int state, double* matrix) {
  /* -1 for the TTT matrix, 0 or greater for the state matrix */

  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */
  int ok = false;
  if (state < 0) {
  } else {
    if (auto* objstate = obj->getObjectState(state)) {
      ObjectStateSetMatrix(objstate, matrix);
      ok = true;
    }
  }
  return ok;
}

int ExecutiveSetObjectMatrix(PyMOLGlobals* G, const char* name, int state,
                             double* matrix) {
  int ok = false;
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    return ExecutiveSetObjectMatrix2(G, obj, state, matrix);
  }
  return ok;
}

pymol::Result<> ExecutiveTranslateAtom(PyMOLGlobals* G, const char* sele,
                                       const float* v, int state, int mode,
                                       int log) {
  SETUP_SELE(sele, tmpsele0, sele0);

  {
    auto obj0 = SelectorGetSingleObjectMolecule(G, sele0);
    if (!obj0) {
      return pymol::make_error("Selection isn't a single atom.");
    } else {
      auto i0 = ObjectMoleculeGetAtomIndex(obj0, sele0);
      if (i0 < 0) {
        return pymol::make_error("Selection isn't a single atom.");
      } else {
        ObjectMoleculeMoveAtom(obj0, state, i0, v, mode, log);
      }
    }
  }
  return {};
}

/**
 * Performs a TTT-altering function whose transformation is optionally stored as
 * a keyframe
 * @param name name of object or pattern
 * @param store to store transformation as keyframe
 * @param func function to call (ex. ExecutiveSetObjectTTT)
 * @param args arguments to func
 */

template <typename Func, typename... Args>
void ExecutiveObjectFuncTTT(PyMOLGlobals* G, pymol::zstring_view name,
                            int store, Func func, Args... args) {
  auto I = G->Executive;
  if (name.empty() || name == cKeywordAll || name == cKeywordSame) {
    for (auto& rec : pymol::make_list_adapter(I->Spec)) {
      switch (rec.type) {
        case cExecObject: {
          pymol::CObject* obj = rec.obj;
          if (ObjectGetSpecLevel(rec.obj, 0) >= 0 || name == cKeywordAll) {
            func(obj, args...);
            obj->invalidate(cRepNone, cRepInvExtents, -1);
          }
        } break;
      }
    }
    if (store && SettingGet<bool>(G, cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  } else { /* pattern */
    for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, name)) {
      switch (rec.type) {
        case cExecObject: {
          pymol::CObject* obj = rec.obj;
          func(obj, args...);
          obj->invalidate(cRepNone, cRepInvExtents, -1);
        } break;
      }
    }
    if (store && SettingGet<bool>(G, cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  }
  SceneInvalidate(G);
}
pymol::Result<> ExecutiveCombineObjectTTT(PyMOLGlobals* G,
                                          pymol::zstring_view name,
                                          const float* ttt, int reverse_order,
                                          int store) {
  ExecutiveObjectFuncTTT(G, name, store, ObjectCombineTTT, ttt, reverse_order,
                         store);
  return {};
}

pymol::Result<> ExecutiveTranslateObjectTTT(PyMOLGlobals* G,
                                            pymol::zstring_view name,
                                            const float* trans, int store,
                                            int quiet) {
  ExecutiveObjectFuncTTT(G, name, store, ObjectTranslateTTT, trans, store);
  return {};
}

pymol::Result<> ExecutiveSetObjectTTT(PyMOLGlobals* G, pymol::zstring_view name,
                                      const float* ttt, int state, int quiet,
                                      int store) {
  ExecutiveObjectFuncTTT(G, name, store, ObjectSetTTT, ttt, state, store);
  return {};
}

int ExecutiveGetObjectTTT(PyMOLGlobals* G, const char* name, const float** ttt,
                          int state, int quiet) {
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  int ok = true;

  if (!obj) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    "Error: object %s not found.\n", name ENDFB(G);
    ok = false;
  } else {
    ObjectGetTTT(obj, ttt, state);
  }
  return (ok);
}

pymol::Result<> ExecutiveTransformSelection(PyMOLGlobals* G, int state,
                                            const char* s1, int log,
                                            const float* ttt, int homogenous) {
  SETUP_SELE_DEFAULT(1);

  auto vla = pymol::vla_take_ownership(SelectorGetObjectMoleculeVLA(G, sele1));
  if (!vla) return pymol::make_error("Could not find selection");
  for (auto obj : vla) {
    ObjectMoleculeTransformSelection(obj, state, sele1, ttt, log,
                                     tmpsele1->getName(), homogenous, true);
  }
  SceneInvalidate(G);
  return {};
}

pymol::Result<> ExecutiveTransformObjectSelection2(PyMOLGlobals* G,
                                                   pymol::CObject* obj,
                                                   int state, const char* s1,
                                                   int log, const float* matrix,
                                                   int homogenous, int global) {
  int ok = true;

  switch (obj->type) {
    case cObjectMolecule: {
      int sele = -1;
      ObjectMolecule* objMol = (ObjectMolecule*)obj;

      if (s1 && s1[0]) {
        sele = SelectorIndexByName(G, s1);
        if (sele < 0) ok = false;
      }
      if (!ok) {
        return pymol::make_error("Selection object ", s1, " not found.");
      } else {
        ObjectMoleculeTransformSelection(objMol, state, sele, matrix, log, s1,
                                         homogenous, global);
      }
      EditorDihedralInvalid(G, objMol);
      SceneInvalidate(G);
    } break;
    default:
      if (auto* objstate = obj->getObjectState(state)) {
        double matrixd[116];
        if (homogenous) {
          convert44f44d(matrix, matrixd);
        } else {
          convertTTTfR44d(matrix, matrixd);
        }
        ObjectStateTransformMatrix(objstate, matrixd);
        obj->invalidate(cRepNone, cRepInvExtents, state);
      }
      break;
  }
  return {};
}

pymol::Result<> ExecutiveTransformObjectSelection(PyMOLGlobals* G,
                                                  const char* name, int state,
                                                  const char* s1, int log,
                                                  const float* matrix,
                                                  int homogenous, int global) {
  pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    return ExecutiveTransformObjectSelection2(G, obj, state, s1, log, matrix,
                                              homogenous, global);
  }
  return {};
}
