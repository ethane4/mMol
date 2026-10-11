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

// Showing, hiding and coloring representations.

#include "pymol/select/atom_iterators.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/seeker.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

static void ExecutiveSetAllRepVisMask(PyMOLGlobals* G, int repmask, int state);

static int* getRepArrayFromBitmask(int visRep);

PyObject* ExecutiveGetVisAsPyDict(PyMOLGlobals* G) {
  assert(PyGILState_Check());

  PyObject *result = nullptr, *list;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  result = PyDict_New();
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->name[0] != '_') {
      list = PyList_New(4);
      PyList_SetItem(list, 0, PyInt_FromLong(rec->visible));

      PyList_SetItem(list, 1, PyList_New(0));

      if (rec->type != cExecObject) {
        PyList_SetItem(list, 2, PConvAutoNone(Py_None));
        PyList_SetItem(list, 3, PConvAutoNone(Py_None));
      } else {
        auto vla = getRepArrayFromBitmask(rec->obj->visRep);
        PyList_SetItem(list, 2, PConvIntVLAToPyList(vla));
        VLAFreeP(vla);

        PyList_SetItem(list, 3, PyInt_FromLong(rec->obj->Color));
      }

      PyDict_SetItemString(result, rec->name, list);
      Py_DECREF(list);
    }
  }
  return (result);
}

static int* getRepArrayFromBitmask(int visRep) {
  int n_vis = 0;
  int* RepVis = VLACalloc(int, cRepCnt);

  for (int a = 0; a < cRepCnt; a++)
    if (GET_BIT(visRep, a)) RepVis[n_vis++] = a;

  VLASize(RepVis, int, n_vis);
  return RepVis;
}

int ExecutiveSetVisFromPyDict(PyMOLGlobals* G, PyObject* dict) {
  assert(PyGILState_Check());

  int ok = true;
  WordType name;
  PyObject *key, *list, *col;
  PyObject* vis_list = nullptr;
  Py_ssize_t pos = 0;
  SpecRec *rec, *grec, **recstack = nullptr;
  int n_vis;
  int rep;
  int a;
  int ll = 0;
  if (ok) ok = (dict != nullptr);
  if (ok) ok = PyDict_Check(dict);
  if (ok) {
    SceneObjectDel(G, nullptr, true); /* remove all objects from scene */
    ExecutiveInvalidateSceneMembers(G);

    // stack for putative visible records
    recstack = pymol::calloc<SpecRec*>(PyDict_Size(dict) + 1);

    while (PyDict_Next(dict, &pos, &key, &list)) {
      if (!PConvPyStrToStr(key, name, sizeof(WordType))) {
        ok = false;
      } else {
        rec = ExecutiveFindSpec(G, name);
        if (rec) {
          if (ok) ok = (list != nullptr);
          if (ok) ok = PyList_Check(list);
          if (ok) ll = PyList_Size(list);
          if (ok) ok = (ll >= 2);
          if (ok)
            ok = PConvPyObjectToInt(PyList_GetItem(list, 0), &rec->visible);

          /* before version 1.8 item 1 was rec reps (repOn) */

          if (ok && (rec->type == cExecObject)) { /* object properties */

            if (ll > 2) { /* object visibility */
              vis_list = PyList_GetItem(list, 2);
              if (ok) ok = (vis_list != nullptr);
              if (ok) {
                if (PyList_Check(vis_list)) {
                  n_vis = PyList_Size(vis_list);
                  rec->obj->visRep = 0;
                  for (a = 0; a < n_vis; a++) {
                    if (PConvPyObjectToInt(PyList_GetItem(vis_list, a), &rep)) {
                      if ((rep >= 0) && (rep < cRepCnt))
                        SET_BIT(rec->obj->visRep, rep);
                    }
                  }
                } else if (PyInt_Check(vis_list)) {
                  PConvPyObjectToInt(vis_list, &rec->obj->visRep);
                }
              }
            }
            if (ll > 3) { /* object color */
              col = PyList_GetItem(list, 3);
              if (ok) ok = (col != nullptr);
              if (ok)
                if (PyInt_Check(col)) {
                  ok = PConvPyObjectToInt(col, &rec->obj->Color);
                  rec->obj->invalidate(cRepAll, cRepInvColor, -1);
                }
            }
          }
          if (rec->visible && (rec->type == cExecObject)) {
            (*(++recstack)) = rec;
          }
        }
      }
    }

    // add visible objects to scene
    for (; (rec = *recstack); recstack--) {
      // check visibility of all parent groups
      for (grec = rec; grec->visible && (grec = grec->group););
      if (!grec) {
        // ok, no invisible parent found
        rec->in_scene = SceneObjectAdd(G, rec->obj);
        ExecutiveInvalidateSceneMembers(G);
      }
    }
    mfree(recstack);
  }
  return ok;
}

pymol::Result<std::pair<float, float>> ExecutiveSpectrum(
    PyMOLGlobals* G, pymol::zstring_view s_s1, pymol::zstring_view s_expr,
    float min, float max, int first, int last, pymol::zstring_view s_prefix,
    int digits, int byres, int quiet) {
  std::pair<float, float> ret;
  int n_color, n_atom;
  ObjectMoleculeOpRec op;
  WordType buffer;
  std::vector<int> color_index;
  std::vector<float> value;
  int a, b;
  char pat[] = "%0Xd";
  int pref_len;
  char* at;
  float range;

  auto s1 = s_s1.c_str();
  auto prefix = s_prefix.c_str();
  auto expr = s_expr.c_str();

  auto tmpsele1 = SelectorTmp::make(G, s1);
  p_return_if_error(tmpsele1);
  int sele1 = tmpsele1->getIndex();

  if (sele1 >= 0) {
    if (digits > 9) digits = 9;
    pat[2] = ('0' + digits);
    UtilNCopy(buffer, prefix, sizeof(WordType) - digits);

    pref_len = strlen(prefix);
    at = buffer + pref_len;

    n_color = abs(first - last) + 1;
    if (n_color) {
      color_index.resize(n_color);
      for (a = 0; a < n_color; a++) {
        b = first + ((last - first) * a) / (n_color - 1);
        sprintf(at, pat, b);
        color_index[a] = ColorGetIndex(G, buffer);
      }

      // set up iterator
      SeleAtomIterator iter(G, sele1);
      SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

      // count atoms
      for (n_atom = 0; iter.next();) {
        ++n_atom;
      }

      if (n_atom) {
        value.resize(n_atom);

        if (WordMatchExact(G, "count", expr, true)) {
          for (a = 0; a < n_atom; a++) {
            value[a] = (float)a + 1;
          }
        } else {
          if (WordMatchExact(G, "pc", expr, true)) {
            expr = "partial_charge";
          } else if (WordMatchExact(G, "resi", expr, true)) {
            expr = "resv";
          }

          // look up expression definition
          auto ap = PyMOL_GetAtomPropertyInfo(G->PyMOL, expr);
          if (!ap) {
            return pymol::make_error("Unknown expression: ", expr);
          }

          // for enumerated values
          std::map<size_t, unsigned int> enumerated_values;
          union {
            size_t value_e;
            char value_s[sizeof(size_t)];
          };

          for (a = 0, iter.reset(); iter.next(); ++a) {
            const auto ai = iter.getAtomInfo();
            const auto raw_ptr = reinterpret_cast<const char*>(ai) + ap->offset;

            // numeric values
            switch (ap->Ptype) {
              case cPType_float:
                value[a] = *reinterpret_cast<const float*>(raw_ptr);
                continue;
              case cPType_int:
              case cPType_int_custom_type:
                value[a] = *reinterpret_cast<const int*>(raw_ptr);
                continue;
              case cPType_uint32:
                value[a] = *reinterpret_cast<const uint32_t*>(raw_ptr);
                continue;
              case cPType_schar:
                value[a] = *reinterpret_cast<const signed char*>(raw_ptr);
                continue;
              case cPType_char_as_type:
                value[a] = ai->hetatm;
                continue;
              case cPType_index:
                value[a] = iter.getAtm() + 1.f;
                continue;
            }

            // enumerated values
            switch (ap->Ptype) {
              case cPType_int_as_string:
                value_e =
                    LexNumeric(*reinterpret_cast<const lexidx_t*>(raw_ptr));
                break;
              case cPType_string:
                // works for small strings
                strncpy(value_s, raw_ptr, sizeof(value_e));
                break;
              case cPType_model:
                value_e = (size_t)iter.obj;
                break;
              default:
                return pymol::make_error("Unsupported Ptype for expr: ", expr);
            }

            // lookup or insert value
            auto& e = enumerated_values[value_e];
            if (e == 0) e = enumerated_values.size();
            value[a] = e - 1.f;
          }

          if (!quiet && !enumerated_values.empty()) {
            PRINTFB(G, FB_Executive, FB_Actions)
            " Spectrum: Expression is non-numeric, enumerating values\n" ENDFB(
                G);
          }
        }

        if (max < min) {
          max = value[0];
          min = value[0];
          for (a = 1; a < n_atom; a++) {
            if (value[a] < min) min = value[a];
            if (value[a] > max) max = value[a];
          }
        }
        range = max - min;

        if (!quiet) {
          PRINTFB(G, FB_Executive, FB_Actions)
          " Spectrum: range (%8.5f to %8.5f).\n", min, max ENDFB(G);
        }
        if (range == 0.0F) range = 1.0F;
        ret.first = min;
        ret.second = max;

        op.code = OMOP_Spectrum;
        op.i1 = n_color - 1;
        op.i2 = n_atom;
        op.i3 = 0;
        op.i4 = byres;
        op.ii1 = color_index.data();
        op.ff1 = value.data();
        op.f1 = min;
        op.f2 = range;

        ExecutiveObjMolSeleOp(G, sele1, &op);

        op.code = OMOP_INVA;
        op.i1 = cRepBitmask;
        op.i2 = cRepInvColor;
        ExecutiveObjMolSeleOp(G, sele1, &op);
      }
    }
  }
  return ret;
}

pymol::Result<int> ExecutiveCartoon(PyMOLGlobals* G, int type, const char* s1) {
  SETUP_SELE_DEFAULT(1);

  ObjectMoleculeOpRec op1;

  ObjectMoleculeOpRecInit(&op1);
  op1.i2 = 0;
  {
    op1.code = OMOP_Cartoon;
    op1.i1 = type;
    op1.i2 = 0;
    op1.i3 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    if (op1.i3 > 0) {
      op1.code = OMOP_INVA;
      op1.i1 = cRepCartoonBit;
      op1.i2 = cRepInvRep;
      ExecutiveObjMolSeleOp(G, sele1, &op1);
    }
  }
  return (op1.i2);
}
pymol::Result<> ExecutiveColorFromSele(PyMOLGlobals* G, const char* sele,
                                       const char* color, int flags,
                                       int quiet) {
  auto s1 = SelectorTmp2::make(G, sele);
  p_return_if_error(s1);
  return ExecutiveColor(G, s1->getName(), color, flags, quiet);
}

pymol::Result<> ExecutiveColor(PyMOLGlobals* G, const char* name,
                               const char* color, int flags, int quiet) {
  /* flags:
     0x1 -- ignore or suppress selection name matches
   */

  CExecutive* I = G->Executive;
  int col_ind;
  int ok = false;
  col_ind = ColorGetIndex(G, color);
  if ((!name) || (!name[0])) name = cKeywordAll;
  if (col_ind == -1) {
    return pymol::make_error("Unknown color.");
  } else {
    CTracker* I_Tracker = I->Tracker;
    SpecRec* rec = nullptr;
    int n_atm = 0;
    int n_obj = 0;

    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecSelection:
          case cExecObject:
          case cExecAll:
            if ((rec->type == cExecSelection) || /* coloring a selection */
                (rec->type == cExecAll) ||       /* coloring all */
                ((rec->type == cExecObject) &&   /* coloring object and its
                                                    backing selection */
                 (rec->obj->type == cObjectMolecule))) {
              if (!(flags & 0x1)) {
                int sele = SelectorIndexByName(G, rec->name);
                ObjectMoleculeOpRec op;
                if (sele >= 0) {
                  ok = true;
                  ObjectMoleculeOpRecInit(&op);
                  op.code = OMOP_COLR;
                  op.i1 = col_ind;
                  op.i2 = n_atm;
                  ExecutiveObjMolSeleOp(G, sele, &op);
                  n_atm = op.i2;
                  op.code = OMOP_INVA;
                  op.i1 = cRepBitmask;
                  op.i2 = cRepInvColor;
                  ExecutiveObjMolSeleOp(G, sele, &op);
                }
              }
            }
            break;
        }

        switch (rec->type) { /* sets object color */
          case cExecObject:
            rec->obj->Color = col_ind;
            rec->obj->invalidate(cRepAll, cRepInvColor, -1);
            n_obj++;
            ok = true;
            SceneInvalidate(G);
            break;
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                rec->obj->Color = col_ind;
                rec->obj->invalidate(cRepAll, cRepInvColor, -1);
                n_obj++;
                ok = true;
                SceneInvalidate(G);
              }
            }
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    if (n_obj || n_atm) {
      char atms[] = "s";
      char objs[] = "s";
      if (n_obj < 2) objs[0] = 0;
      if (n_atm < 2) atms[0] = 0;
      if (!quiet) {
        if (n_obj && n_atm) {
          PRINTFB(G, FB_Executive, FB_Actions)
          " Executive: Colored %d atom%s and %d object%s.\n", n_atm, atms,
              n_obj, objs ENDFB(G);
        } else if (n_obj) {
          PRINTFB(G, FB_Executive, FB_Actions)
          " Executive: Colored %d object%s.\n", n_obj, objs ENDFB(G);
        } else {
          PRINTFB(G, FB_Executive, FB_Actions)
          " Executive: Colored %d atom%s.\n", n_atm, atms ENDFB(G);
        }
      }
    }
  }
  return {};
}
pymol::Result<bool> ExecutiveSetObjVisib(PyMOLGlobals* G,
                                         pymol::zstring_view name, int onoff,
                                         int parents) {
  CExecutive* I = G->Executive;
  bool changed = false;
  PRINTFD(G, FB_Executive)
  " ExecutiveSetObjVisib: entered.\n" ENDFD;
  {
    CTracker* I_Tracker = I->Tracker;
    SpecRec* rec;
    int list_id = ExecutiveGetNamesListFromPattern(G, name.data(), true, false);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        if (!changed && rec->visible != onoff) {
          changed = true;
        }
        switch (rec->type) {
          case cExecAll: {
            bool const suppress_hidden =
                SettingGet<bool>(G, cSetting_suppress_hidden);
            bool const hide_underscore =
                SettingGet<bool>(G, cSetting_hide_underscore_names);
            if (suppress_hidden && hide_underscore)
              ExecutiveUpdateGroups(G, false);
            SpecRec* tRec = nullptr;
            while (ListIterate(I->Spec, tRec, next)) {
              if (onoff != tRec->visible) {
                if (tRec->type == cExecObject) {
                  if (tRec->visible) {
                    tRec->in_scene = SceneObjectDel(G, tRec->obj, true);
                    ExecutiveInvalidateSceneMembers(G);
                    tRec->visible = !tRec->visible;
                    ReportEnabledChange(G, rec);
                  } else {
                    if (!(suppress_hidden && tRec->isHidden(hide_underscore))) {
                      tRec->in_scene = SceneObjectAdd(G, tRec->obj);
                      ExecutiveInvalidateSceneMembers(G);
                      tRec->visible = !tRec->visible;
                      ReportEnabledChange(G, rec);
                    }
                  }
                } else if ((tRec->type != cExecSelection) ||
                           (!onoff)) { /* hide all selections, but show all */
                  tRec->visible = !tRec->visible;
                  ReportEnabledChange(G, rec);
                }
              }
            }
          } break;
          case cExecObject:
            /*
               if(rec->visible!=onoff) {
               if(rec->visible) {
               rec->in_scene = SceneObjectDel(G,rec->obj);
               ExecutiveInvalidateSceneMembers(G);
               } else {
               rec->in_scene = SceneObjectAdd(G,rec->obj);
               ExecutiveInvalidateSceneMembers(G);
               }
               rec->visible=!rec->visible;
               }
             */
            if (onoff) { /* enable */
              ExecutiveSpecEnable(G, rec, parents, false);
            } else { /* disable */
              if (rec->visible) {
                if (rec->in_scene)
                  rec->in_scene = SceneObjectDel(G, rec->obj, true);
                rec->visible = false;
                ExecutiveInvalidateSceneMembers(G);
                ReportEnabledChange(G, rec);
              }
              SceneChanged(G);
            }
            break;
          case cExecSelection:
            if (rec->visible != onoff) {
              int previousVisible = rec->visible;
              rec->visible = !rec->visible;
              if (rec->visible)
                if (SettingGetGlobal_b(G, cSetting_active_selections)) {
                  ExecutiveHideSelections(G);
                  rec->visible = true;
                }
              SceneInvalidate(G);
              SeqDirty(G);
              if (previousVisible != rec->visible) {
                ReportEnabledChange(G, rec);
              }
            }
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  PRINTFD(G, FB_Executive)
  " ExecutiveSetObjVisib: leaving...\n" ENDFD;
  return changed;
}
static void fInvalidateRepMask(pymol::CObject* obj, cRepBitmask_t repmask,
                               int state = -1) {
  for (auto a = cRep_t(0); a < cRepCnt; ++a) {
    if ((1 << a) & repmask) obj->invalidate(a, cRepInvVisib, state);
  }
}
pymol::Result<> ExecutiveToggleRepVisib(PyMOLGlobals* G, const char* name,
                                        int rep) {
  int sele = -1;
  SpecRec* tRec;
  ObjectMoleculeOpRec op;
  OrthoLineType tmpname;

  PRINTFD(G, FB_Executive)
  " ExecutiveToggleRepVisib: entered.\n" ENDFD;

  tRec = ExecutiveFindSpec(G, name);

  if (rep == -2) {
    // special case: toggle object visibility (should that be in this function?)
    if (tRec) {
      ExecutiveSetObjVisib(G, name, !tRec->visible, 0);
    } else {
      return pymol::make_error(name, " not found.");
    }
  } else if (tRec && tRec->type == cExecObject &&
             tRec->obj->type != cObjectMolecule) {
    // non-atom object
    tRec->obj->visRep ^= rep;
    fInvalidateRepMask(tRec->obj, rep, 0);
    SceneChanged(G);
  } else if (SelectorGetTmp(G, name, tmpname) >= 0) {
    // atom selection
    sele = SelectorIndexByName(G, tmpname);
    if (sele >= 0) {
      ObjectMoleculeOpRecInit(&op);

      op.code = OMOP_CheckVis;
      op.i1 = rep;
      op.i2 = false;
      ExecutiveObjMolSeleOp(G, sele, &op);
      op.i2 = !op.i2;

      op.code = OMOP_VISI;
      op.i1 = rep;
      ExecutiveObjMolSeleOp(G, sele, &op);
      op.code = OMOP_INVA;
      op.i2 = cRepInvVisib;
      ExecutiveObjMolSeleOp(G, sele, &op);
    }
    SelectorFreeTmp(G, tmpname);
  }
  PRINTFD(G, FB_Executive)
  " ExecutiveToggleRepVisib: leaving...\n" ENDFD;
  return {};
}
pymol::Result<> ExecutiveSetRepVisib(PyMOLGlobals* G, pymol::zstring_view name,
                                     int rep, int state) {
  int repmask = (rep == cRepAll) ? cRepBitmask : (1 << rep);
  return ExecutiveSetRepVisMask(G, name, repmask, state);
}

pymol::Result<> ExecutiveSetRepVisMask(PyMOLGlobals* G,
                                       pymol::zstring_view name, int repmask,
                                       int state) {
  PRINTFD(G, FB_Executive)
  " ExecutiveSetRepVisib: entered.\n" ENDFD;

  {
    CExecutive* I = G->Executive;
    CTracker* I_Tracker = I->Tracker;
    SpecRec* rec = nullptr;
    int list_id = ExecutiveGetNamesListFromPattern(G, name.data(), true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        /* per-atom */

        switch (rec->type) {
          case cExecObject:
            if (rec->obj->type != cObjectMolecule &&
                rec->obj->type != cObjectAlignment)
              break;
          case cExecSelection: {
            int sele = SelectorIndexByName(G, rec->name);
            if (sele >= 0) {
              ObjectMoleculeOpRec op;
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_VISI;
              op.i1 = repmask;
              op.i2 = state;
              ExecutiveObjMolSeleOp(G, sele, &op);
              op.code = OMOP_INVA;
              if (state == cVis_AS) op.i1 = cRepBitmask;
              op.i2 = cRepInvVisib;
              ExecutiveObjMolSeleOp(G, sele, &op);
            }
          } break;
        }

        /* per-object/name */

        switch (rec->type) {
          case cExecObject:
            ObjectSetRepVisMask(rec->obj, repmask, state);
            fInvalidateRepMask(rec->obj, repmask, 0);
            SceneChanged(G);
            break;
          case cExecAll:
            ExecutiveSetAllRepVisMask(G, repmask, state);
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  PRINTFD(G, FB_Executive)
  " ExecutiveSetRepVisib: leaving...\n" ENDFD;
  return {};
}
pymol::Result<> ExecutiveSetOnOffBySele(PyMOLGlobals* G,
                                        pymol::zstring_view sname, int onoff) {
  SelectorTmp2 tmp{G, sname.data()};
  const char* name = tmp.getName();

  auto tRec = ExecutiveFindSpec(G, name);
  if (!tRec && sname == cKeywordAll) {
    ExecutiveSetObjVisib(G, name, onoff, false);
  }
  if (tRec) {
    int sele = tmp.getIndex();
    if (sele >= 0) {
      ObjectMoleculeOpRec op;
      ObjectMoleculeOpRecInit(&op);

      op.code = OMOP_OnOff;
      op.i1 = onoff;
      ExecutiveObjMolSeleOp(G, sele, &op);
    }
  }

  return {};
}
/**
 * @param repmask rep bit mask
 * @param state 0 (hide), 1 (show), 2 (as)
 */
static void ExecutiveSetAllRepVisMask(PyMOLGlobals* G, int repmask, int state) {
  ObjectMoleculeOpRec op;
  ObjectMolecule* obj;
  int sele;

  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  PRINTFD(G, FB_Executive)
  " ExecutiveSetAllRepVisib: entered.\n" ENDFD;
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      if (rec->type == cExecObject) {
        switch (rec->obj->type) {
          case cObjectMolecule:
            obj = (ObjectMolecule*)rec->obj;
            sele = SelectorIndexByName(G, obj->Name);
            ObjectMoleculeOpRecInit(&op);

            op.code = OMOP_VISI;
            op.i1 = repmask;
            op.i2 = state;
            ObjectMoleculeSeleOp(obj, sele, &op);
            op.code = OMOP_INVA;
            if (state == cVis_AS) op.i1 = cRepBitmask;
            op.i2 = cRepInvVisib;
            ObjectMoleculeSeleOp(obj, sele, &op);
            break;
          default:
            ObjectSetRepVisMask(rec->obj, repmask, state);
            fInvalidateRepMask(rec->obj, repmask, -1);
            SceneInvalidate(G);
            break;
        }
      }
    }
  }
  PRINTFD(G, FB_Executive)
  " ExecutiveSetAllRepVisib: leaving...\n" ENDFD;
}
pymol::Result<> ExecutiveInvalidateRep(PyMOLGlobals* G, const char* str1,
                                       cRep_t rep, cRepInv_t level) {
  CExecutive* I = G->Executive;
  ObjectMoleculeOpRec op;
  SpecRec* rec = nullptr;
  const char* name = "";

  SelectorTmp2 s1;
  if (str1 && !WordMatchExact(G, str1, "all", true)) {
    s1 = SelectorTmp2(G, str1);
    name = s1.getName();
  } else {
    name = str1;
  }

  if ((!name) || (!name[0])) name = cKeywordAll;
  {
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);

    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecSelection:
          case cExecObject: {
            int sele = SelectorIndexByName(G, rec->name);
            if (sele >= 0) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_INVA;
              op.i1 = (rep == cRepAll) ? cRepBitmask : (1 << rep);
              op.i2 = level;
              ExecutiveObjMolSeleOp(G, sele, &op);
            } else {
              rec->obj->invalidate(rep, level, -1);
            }
          } break;
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                rec->obj->invalidate(rep, level, -1);
              }
            }
            SceneInvalidate(G);
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return {};
}
