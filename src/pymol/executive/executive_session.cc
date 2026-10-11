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

// Saving and restoring the object list in sessions.

#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/scene/movie.h"
#include "pymol/objects/object_alignment.h"
#include "pymol/objects/object_cgo.h"
#include "pymol/objects/object_callback.h"
#include "pymol/objects/object_curve.h"
#include "pymol/objects/object_dist.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/objects/object_group.h"
#include "pymol/objects/object_mesh.h"
#include "pymol/objects/object_slice.h"
#include "pymol/objects/object_surface.h"
#include "pymol/objects/object_volume.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/executive/setting_effects.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/app/main.h"
#include "pymol/scene/movie_scene.h"
#include "pymol/executive/executive_internal.h"

static int ExecutiveSetNamedEntries(PyMOLGlobals* G, PyObject* names,
                                    int version, int part_rest, int part_sess) {
  CExecutive* I = G->Executive;
  int ok = true;
  int skip = false;
  int a = 0, l = 0, ll = 0;
  PyObject *cur, *el;
  SpecRec* rec = nullptr;
  int extra_int;
  int incomplete = false;
  ObjectNameType new_name;

  if (ok) ok = (names != nullptr);
  if (ok) ok = PyList_Check(names);
  if (ok) l = PyList_Size(names);

  while (ok && (a < l)) {
    cur = PyList_GetItem(names, a);
    if (cur != Py_None) { /* skip over None w/o aborting */
      skip = false;
      rec = nullptr;
      ListElemCalloc(G, rec, SpecRec);
      rec->next = nullptr;
      rec->name[0] = 0;
      if (ok) ok = PyList_Check(cur);
      if (ok) ll = PyList_Size(cur);
      if (ok)
        ok = PConvPyStrToStr(PyList_GetItem(cur, 0), rec->name,
                             sizeof(WordType));
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(cur, 1), &rec->type);
      if (ok)
        ok = CPythonVal_PConvPyIntToInt_From_List(G, cur, 2, &rec->visible);

      /* before version 1.8 item 3 was rec reps (repOn) */

      if (ok) ok = PConvPyIntToInt(PyList_GetItem(cur, 4), &extra_int);
      switch (rec->type) {
        case cExecObject:
          if (!ok) break;

          el = PyList_GetItem(cur, 5);

          switch (extra_int) {
            case cObjectMolecule:
              ok = ObjectMoleculeNewFromPyList(
                  G, el, (ObjectMolecule**)(void*)&rec->obj);
              break;
            case cObjectMeasurement:
              ok = ObjectDistNewFromPyList(G, el,
                                           (ObjectDist**)(void*)&rec->obj);
              break;
            case cObjectMap:
              ok = ObjectMapNewFromPyList(G, el, (ObjectMap**)(void*)&rec->obj);
              break;
            case cObjectMesh:
              ok = ObjectMeshNewFromPyList(G, el,
                                           (ObjectMesh**)(void*)&rec->obj);
              break;
            case cObjectSlice:
              ok = ObjectSliceNewFromPyList(G, el,
                                            (ObjectSlice**)(void*)&rec->obj);
              break;
            case cObjectSurface:
              ok = ObjectSurfaceNewFromPyList(
                  G, el, (ObjectSurface**)(void*)&rec->obj);
              break;
            case cObjectCGO:
              ok = ObjectCGONewFromPyList(G, el, (ObjectCGO**)(void*)&rec->obj,
                                          version);
              break;
            case cObjectGadget:
              ok = ObjectGadgetNewFromPyList(
                  G, el, (ObjectGadget**)(void*)&rec->obj, version);
              break;
            case cObjectAlignment:
              ok = ObjectAlignmentNewFromPyList(
                  G, el, (ObjectAlignment**)(void*)&rec->obj, version);
              break;
            case cObjectGroup:
              if (part_rest) {
                // if group already exists, do not create new one
                pymol::CObject* obj = ExecutiveFindObjectByName(G, rec->name);
                if (obj && obj->type == cObjectGroup) {
                  skip = 1;
                  break;
                }
              }
              ok = ObjectGroupNewFromPyList(
                  G, el, (ObjectGroup**)(void*)&rec->obj, version);
              break;
            case cObjectVolume:
              ok = ObjectVolumeNewFromPyList(G, el,
                                             (ObjectVolume**)(void*)&rec->obj);
              break;
            case cObjectCallback:
              // skip dummy entries from old sessions and failed-to-pickle
              // sessions
              skip = !ObjectCallbackNewFromPyList(
                  G, el, (ObjectCallback**)(void*)&rec->obj);
              break;
            case cObjectCurve: {
              rec->obj = new ObjectCurve(G, el);
            } break;
            default:
              PRINTFB(G, FB_Executive, FB_Errors)
              " Executive: skipping unrecognized object \"%s\" of type %d.\n",
                  rec->name, extra_int ENDFB(G);
              skip = true;
              break;
          }

          CPythonVal_Free(el);

          break;
        case cExecSelection:  // on the first pass, just create an entry in the
                              // rec list
          rec->sele_color = extra_int;
          if (part_rest || part_sess) {  // don't attempt to restore selections
                                         // with partial sessions
            skip = true;
          }
          break;
      }

      if (ll > 6) {
        if (ok) {
          ok = PConvPyStrToStr(PyList_GetItem(cur, 6), rec->group_name,
                               sizeof(WordType));
        }
      }

      if (PyErr_Occurred()) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetNamedEntries-Error: after object \"%s\".\n",
            rec->name ENDFB(G);
        PyErr_Print();
      }

      if (ok && !skip) {
        if (part_rest && ExecutiveProcessObjectName(G, rec->name, new_name)) {
          // rename duplicates
          strcpy(rec->obj->Name, new_name);
          strcpy(rec->name, new_name);
        }

        // replace existing object (unless auto_rename_duplicate_objects=1)
        if (ExecutiveValidName(G, rec->name)) {
          ExecutiveDelete(G, rec->name);
        }

        switch (rec->type) {
          case cExecObject:
            if (rec->visible) {
              rec->in_scene = SceneObjectAdd(G, rec->obj);
              ExecutiveInvalidateSceneMembers(G);
            }
            ExecutiveUpdateObjectSelection(G, rec->obj);
            break;
        }

        rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef*)rec);
        TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);

        switch (rec->type) {
          case cExecObject:
            TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id, 1);
            break;
          case cExecSelection:
            TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id, 1);
            break;
        }
        ListAppend(I->Spec, rec, next, SpecRec);
        ExecutiveAddKey(I, rec);
        ExecutiveInvalidateGroups(G, false);
        ExecutiveInvalidatePanelList(G);
      } else {
        ListElemFree(rec);
      }
    }
    a++;
    if (!ok) {
      incomplete = true;
      ok = true;
    }
  }
  return (!incomplete);
}

static int ExecutiveSetSelectionsFromPyList(PyMOLGlobals* G, PyObject* names) {
  /* must already have objects loaded at this point... */

  int ok = true;
  int a = 0, l = 0;
  PyObject* cur;
  SpecRec* rec = nullptr;
  int extra;
  int incomplete = false;

  if (ok) ok = (names != nullptr);
  if (ok) ok = PyList_Check(names);
  if (ok) l = PyList_Size(names);
  while (ok && (a < l)) {
    cur = PyList_GetItem(names, a);
    if (cur != Py_None) { /* skip over None w/o aborting */
      rec = nullptr;
      ListElemCalloc(G, rec, SpecRec);
      rec->next = nullptr;

      if (ok) ok = PyList_Check(cur);
      if (ok)
        ok = PConvPyStrToStr(PyList_GetItem(cur, 0), rec->name,
                             sizeof(WordType));
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(cur, 1), &rec->type);
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(cur, 2), &rec->visible);
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(cur, 4), &extra);
      switch (rec->type) {
        case cExecSelection:
          ok = SelectorFromPyList(G, rec->name, PyList_GetItem(cur, 5));
          break;
      }
      ListElemFree(rec);
    }
    a++;
    if (!ok) {
      incomplete = true;
      ok = true;
    }
  }
  return (!incomplete);
}

static PyObject* ExecutiveGetExecObjectAsPyList(PyMOLGlobals* G, SpecRec* rec) {
  PyObject* result = nullptr;
  int recobjtype = rec->obj->type;
  switch (recobjtype) {
    case cObjectMesh: { /* If a mesh no longer has its dependent map, then it
                           gets saved as a CGO */
      int allMapsExist = ObjectMeshAllMapsInStatesExist((ObjectMesh*)rec->obj);
      if (!allMapsExist) {
        recobjtype = cObjectCGO;
      }
    }
  }
  result = PyList_New(7);
  PyList_SetItem(result, 0, PyString_FromString(rec->obj->Name));
  PyList_SetItem(result, 1, PyInt_FromLong(cExecObject));
  PyList_SetItem(result, 2, PyInt_FromLong(rec->visible));
  /* before version 1.8 item 3 was rec reps (repOn) */
  PyList_SetItem(result, 3, PConvAutoNone(nullptr));
  PyList_SetItem(result, 4, PyInt_FromLong(recobjtype));
  switch (rec->obj->type) {
    case cObjectGadget:
      PyList_SetItem(result, 5, ObjectGadgetAsPyList((ObjectGadget*)rec->obj));
      break;
    case cObjectMolecule:
      PyList_SetItem(result, 5,
                     ObjectMoleculeAsPyList((ObjectMolecule*)rec->obj));
      break;
    case cObjectMeasurement:
      PyList_SetItem(result, 5, ObjectDistAsPyList((ObjectDist*)rec->obj));
      break;
    case cObjectMap:
      PyList_SetItem(result, 5, ObjectMapAsPyList((ObjectMap*)rec->obj));
      break;
    case cObjectMesh:
      PyList_SetItem(result, 5, ObjectMeshAsPyList((ObjectMesh*)rec->obj));
      break;
    case cObjectSlice:
      PyList_SetItem(result, 5, ObjectSliceAsPyList((ObjectSlice*)rec->obj));
      break;
    case cObjectSurface:
      PyList_SetItem(result, 5,
                     ObjectSurfaceAsPyList((ObjectSurface*)rec->obj));
      break;
    case cObjectCGO:
      PyList_SetItem(result, 5, ObjectCGOAsPyList((ObjectCGO*)rec->obj));
      break;
    case cObjectAlignment:
      PyList_SetItem(result, 5,
                     ObjectAlignmentAsPyList((ObjectAlignment*)rec->obj));
      break;
    case cObjectGroup:
      PyList_SetItem(result, 5, ObjectGroupAsPyList((ObjectGroup*)rec->obj));
      break;
    case cObjectVolume:
      PyList_SetItem(result, 5, ObjectVolumeAsPyList((ObjectVolume*)rec->obj));
      break;
    case cObjectCallback:
      PyList_SetItem(result, 5,
                     ObjectCallbackAsPyList((ObjectCallback*)rec->obj));
      break;
    case cObjectCurve:
      PyList_SetItem(result, 5,
                     static_cast<ObjectCurve*>(rec->obj)->asPyList());
      break;
    default:
      PyList_SetItem(result, 5, PConvAutoNone(nullptr));
      break;
  }
  PyList_SetItem(result, 6, PyString_FromString(rec->group_name));

  return (result);
}

static PyObject* ExecutiveGetExecSeleAsPyList(PyMOLGlobals* G, SpecRec* rec) {
  PyObject* result = nullptr;
  int sele;

  sele = SelectorIndexByName(G, rec->name);
  if (sele >= 0) {
    result = PyList_New(7);
    PyList_SetItem(result, 0, PyString_FromString(rec->name));
    PyList_SetItem(result, 1, PyInt_FromLong(cExecSelection));
    PyList_SetItem(result, 2, PyInt_FromLong(rec->visible));
    /* before version 1.8 item 3 was rec reps (repOn) */
    PyList_SetItem(result, 3, PConvAutoNone(nullptr));
    PyList_SetItem(result, 4, PyInt_FromLong(-1));
    PyList_SetItem(result, 5, SelectorAsPyList(G, sele));
    PyList_SetItem(result, 6, PyString_FromString(rec->group_name));
  }
  return (PConvAutoNone(result));
}

static PyObject* ExecutiveGetNamedEntries(PyMOLGlobals* G, int list_id,
                                          int partial) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  PyObject* result = nullptr;
  int count = 0, total_count = 0;
  int iter_id = 0;
  SpecRec *rec = nullptr, *list_rec = nullptr;

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  if (list_id) {
    total_count = TrackerGetNCandForList(I_Tracker, list_id);
    iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  } else {
    total_count = ExecutiveCountNames(G);
  }
  result = PyList_New(total_count);

  /* critical reliance on short-circuit behavior */

  while ((iter_id && TrackerIterNextCandInList(
                         I_Tracker, iter_id, (TrackerRef**)(void*)&list_rec)) ||
         ((!iter_id) && ListIterate(I->Spec, rec, next))) {
    if (list_id) rec = list_rec;
    if (count >= total_count) break;
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          PyList_SetItem(result, count, ExecutiveGetExecObjectAsPyList(G, rec));
          break;
        case cExecSelection:
          if (!partial) {
            PyList_SetItem(result, count, ExecutiveGetExecSeleAsPyList(G, rec));
          } else {
            /* cannot currently save selections in partial sessions */
            PyList_SetItem(result, count, PConvAutoNone(nullptr));
          }
          break;
        default:
          PyList_SetItem(result, count, PConvAutoNone(nullptr));
          break;
      }
    } else {
      PyList_SetItem(result, count, PConvAutoNone(nullptr));
    }
    count++;
  }

  while (
      count <
      total_count) { /* insure that all members of outgoing list are defined */
    PyList_SetItem(result, count, PConvAutoNone(nullptr));
    count++;
  }

  if (iter_id) {
    TrackerDelIter(I_Tracker, iter_id);
  }
  return (PConvAutoNone(result));
}

int ExecutiveGetSession(PyMOLGlobals* G, PyObject* dict, const char* names,
                        int partial, int quiet) {
  assert(PyGILState_Check());

  int list_id = 0;
  SceneViewType sv;
  PyObject* tmp;

  if (names && names[0]) {
    list_id =
        ExecutiveGetNamesListFromPattern(G, names, true, cExecExpandKeepGroups);
  }

  tmp = MovieScenesAsPyList(G);
  PyDict_SetItemString(dict, "moviescenes", tmp);
  Py_XDECREF(tmp);

  tmp = PyInt_FromLong(_PyMOL_VERSION_int);
  PyDict_SetItemString(dict, "version", tmp);
  Py_XDECREF(tmp);

  tmp = ExecutiveGetNamedEntries(G, list_id, partial);
  PyDict_SetItemString(dict, "names", tmp);
  Py_XDECREF(tmp);

  tmp = ColorAsPyList(G);
  PyDict_SetItemString(dict, "colors", tmp);
  Py_XDECREF(tmp);

  tmp = ColorExtAsPyList(G);
  PyDict_SetItemString(dict, "color_ext", tmp);
  Py_XDECREF(tmp);

  tmp = SettingUniqueAsPyList(G);
  PyDict_SetItemString(dict, "unique_settings", tmp);
  Py_XDECREF(tmp);

  if (partial) { /* mark this as a partial session */

    PyDict_SetItemString(dict, "partial", PConvAutoNone(Py_None));

  } else {
    /* none of the following information is saved in partial sessions */

    tmp = SelectorSecretsAsPyList(G);
    PyDict_SetItemString(dict, "selector_secrets", tmp);
    Py_XDECREF(tmp);

    tmp = SettingGetGlobalsAsPyList(G);
    PyDict_SetItemString(dict, "settings", tmp);
    Py_XDECREF(tmp);

    SceneGetView(G, sv);
    tmp = PConvFloatArrayToPyList(sv, cSceneViewSize);
    PyDict_SetItemString(dict, "view", tmp);
    Py_XDECREF(tmp);

    tmp = MovieAsPyList(G);
    PyDict_SetItemString(dict, "movie", tmp);
    Py_XDECREF(tmp);

    tmp = EditorAsPyList(G);
    PyDict_SetItemString(dict, "editor", tmp);
    Py_XDECREF(tmp);

    tmp = MainAsPyList(G);
    PyDict_SetItemString(dict, "main", tmp);
    Py_XDECREF(tmp);
  }

  return true;
}

static void ExecutiveMigrateSession(PyMOLGlobals* G, int session_version) {
  if (session_version < 1700) {
    if (SettingGetGlobal_i(G, cSetting_seq_view_label_color) == 0 /* white */) {
      SettingSetGlobal_i(G, cSetting_seq_view_label_color, cColorFront);
    }
  }
  if (session_version < 100) {
    /* migrate lighting model */
    SettingSetGlobal_f(G, cSetting_direct,
                       1.8 * SettingGetGlobal_f(G, cSetting_direct));
    SettingSetGlobal_f(G, cSetting_reflect,
                       0.5 * SettingGetGlobal_f(G, cSetting_reflect));
    SettingSetGlobal_f(G, cSetting_ambient,
                       1.166 * SettingGetGlobal_f(G, cSetting_ambient));
    SettingSetGlobal_f(G, cSetting_gamma,
                       0.769 * SettingGetGlobal_f(G, cSetting_gamma));

    /* try best to meet existing expectations with existing sessions */
    SettingSetGlobal_f(G, cSetting_ray_legacy_lighting, 1.0F);

    /* force use of movie_delay in preference to movie_fps */

    SettingSetGlobal_f(G, cSetting_movie_fps, 0.0F);

    /* and labels */
    SettingSetGlobal_i(G, cSetting_label_digits, 2);
  }
  if (session_version < 99) {
    SettingSetGlobal_f(G, cSetting_cartoon_ladder_mode, 0);
    SettingSetGlobal_f(G, cSetting_cartoon_tube_cap, 0);
    SettingSetGlobal_f(G, cSetting_cartoon_nucleic_acid_mode, 1);
    {
      float old_sulfur[3] = {1.0, 0.5, 0.0};
      ColorDef(G, "sulfur", old_sulfur, 0, true);
    }
  }
  if (session_version < 98) {
    /* produce expected rendering quality & performance with old sessions */
    SettingSetGlobal_b(G, cSetting_ray_orthoscopic, 1);
  }
  if (session_version < 96) {
    SettingSetGlobal_f(G, cSetting_ray_transparency_contrast, 1.0F);
  }
  if (session_version < 95) {
    {
      /* adjust fog to reflect current importance of seeing to the Z-slab center
       * w/o fog */

      float fog_start = SettingGetGlobal_f(G, cSetting_fog_start);
      float ray_trace_fog_start =
          SettingGetGlobal_f(G, cSetting_ray_trace_fog_start);
      if ((fog_start == 0.40F) || (fog_start == 0.35F) ||
          (fog_start == 0.30F)) {
        SettingSetGlobal_f(G, cSetting_fog_start, 0.45F);
      }
      if ((ray_trace_fog_start == 0.45F) || (ray_trace_fog_start == 0.40F) ||
          (ray_trace_fog_start == 0.35F)) {
        SettingSetGlobal_f(G, cSetting_ray_trace_fog_start, 0.50F);
      }
    }

    { /* adjust GUI width */

      int gui_width = SettingGetGlobal_i(G, cSetting_internal_gui_width);

      if (gui_width == 160) {
        SettingSetGlobal_i(G, cSetting_internal_gui_width, 220);
      }
    }

    { /* enable antialiasing */

      int antialias = SettingGetGlobal_i(G, cSetting_antialias);

      if (antialias == 0) {
        SettingSetGlobal_i(G, cSetting_antialias, 1);
      }
    }
  }
}

int ExecutiveSetSession(PyMOLGlobals* G, PyObject* session, int partial_restore,
                        int quiet) {
  assert(PyGILState_Check());

  int ok = true;
  int incomplete = false;
  PyObject* tmp;
  SceneViewType sv;
  int version = -1, version_full;
  int migrate_sessions = SettingGetGlobal_b(G, cSetting_session_migration);
  char active[WordLength] = "";
  int have_active = false;
  int partial_session = false;

  if (!partial_restore) { /* if user has requested partial restore */
    ExecutiveDelete(G, "all");
    ColorReset(G);
  }

  if (!session || !PyDict_Check(session)) {
    PRINTFB(G, FB_Executive, FB_Errors)
    "Error: not a dict\n" ENDFB(G);
    return 0;
  }

  if (ok) { /* if session is partial, then don't error about missing stuff */
    tmp = PyDict_GetItemString(session, "partial");
    if (tmp) {
      partial_session = true;
    }
  }

  if (partial_restore) {
    G->SettingUnique->old2new =
        std::make_unique<std::unordered_map<int, int>>();
  }

  if (ok) {
    tmp = PyDict_GetItemString(session, "version");
    if (tmp) {
      ok = PConvPyIntToInt(tmp, &version);
      if (ok) {
        version_full = version;
        while (version_full <
               210) /* any version less than 2.1 (account for next major version
                       2) should be 4 digits, otherwise 3 */
          version_full *= 10;
        float version_float =
            version_full / (version >= 1000000 ? 1000000.f : 1000.f);
        if (version > _PyMOL_VERSION_int) {
          if (!quiet) {
            PRINTFB(G, FB_Executive, FB_Errors)
            "Warning: This session was created with a newer version of PyMOL "
            "(%1.6f).\n",
                version_float ENDFB(G);
            if (SettingGetGlobal_i(G, cSetting_session_version_check)) {
              PRINTFB(G, FB_Executive, FB_Errors)
              "Error: Please update first -- see http://www.pymol.org\n" ENDFB(
                  G);
              ok = false;
            } else {
              PRINTFB(G, FB_Executive, FB_Errors)
              "Warning: Some content may not load completely.\n" ENDFB(G);
            }
          }
        } else {
          if (!quiet) {
            PRINTFB(G, FB_Executive, FB_Details)
            " Executive: Loading version %1.6f session...\n",
                version_float ENDFB(G);
          }
        }
      }
    }
  }
  if (ok) {
    tmp = PyDict_GetItemString(session, "eval_nag");
    if (tmp) {
      ok = PyString_Check(tmp);
      if (ok) {
        const char* st = PyString_AsString(tmp);
        if (st) {
          if (Feedback(G, FB_Nag, FB_Warnings)) {
            OrthoAddOutput(G, st);
          }
        }
      }
    }
  }

  if (ok) { /* colors must be restored before settings and objects */
    tmp = PyDict_GetItemString(session, "colors");
    if (tmp) {
      ok = ColorFromPyList(G, tmp, partial_restore);
    }

    if (tmp || (!partial_restore)) { /* ignore missing if partial restore */

      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after colors.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }

  if (ok) {
    tmp = PyDict_GetItemString(session, "color_ext");
    if (tmp) {
      ok = ColorExtFromPyList(G, tmp, partial_restore);
    }

    if (tmp || (!partial_session)) { /* ignore missing if partial restore */
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after color_ext.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }
  if (ok) {
    tmp = PyDict_GetItemString(session, "settings");
    if (tmp && !partial_restore) {
      SettingInitGlobal(G, false, false, /* use_defaults */ true);
      SettingSetGlobalsFromPyList(G, tmp);
      ColorUpdateFrontFromSettings(G);
    }

    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after settings.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }
  if (ok) {
    tmp = PyDict_GetItemString(session, "unique_settings");
    if (tmp) {
      ok = SettingUniqueFromPyList(G, tmp, partial_restore);
    }

    if (tmp || (!partial_session)) { /* ignore missing if partial restore */
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after settings.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }
  if (ok) {
    tmp = PyDict_GetItemString(session, "names");
    if (tmp) {
      if (ok)
        ok = ExecutiveSetNamedEntries(G, tmp, version, partial_restore,
                                      partial_session);
      if (!(partial_restore || partial_session)) {
        if (ok) ok = ExecutiveSetSelectionsFromPyList(G, tmp);
        if (ok)
          have_active = ExecutiveGetActiveSeleName(G, active, false, false);
      }
    }
    if (PyErr_Occurred()) {
      PyErr_Print();
      ok = false;
    }
    if (!ok) {
      PRINTFB(G, FB_Executive, FB_Errors)
      "ExectiveSetSession-Error: after names.\n" ENDFB(G);
    }
    if (!ok) {
      incomplete = true;
      ok = true; /* keep trying...don't give up */
    }
  }
  if (ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "selector_secrets");
    if (tmp) {
      if (ok) ok = SelectorSecretsFromPyList(G, tmp);
    }

    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after selector secrets.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }
  if (ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "view");
    if (tmp) {
      ok = PConvPyListToFloatArrayInPlace(tmp, sv, cSceneViewSize);
    }
    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */
      if (ok) SceneSetView(G, sv, true, 0, 0);
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after view.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }

  if (ok && !partial_restore) {
    tmp = CPythonVal_PyDict_GetItemString(G, session, "moviescenes");
    if (tmp) {
      MovieScenesFromPyList(G, tmp);
      CPythonVal_Free(tmp);
    }
  }

  if (ok && !(partial_restore)) {
    int warning;
    tmp = PyDict_GetItemString(session, "movie");
    if (tmp) {
      ok = MovieFromPyList(G, tmp, &warning);
    }
    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */

      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after movie.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }

  if (ok && !(partial_restore)) {
    tmp = PyDict_GetItemString(session, "editor");
    if (tmp) {
      ok = EditorFromPyList(G, tmp);
    }
    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */

      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after editor.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }

  if (ok) {
    PParse(G, "cmd.mouse(quiet=1)");
  }
  // Do not load viewport size when we have a GUI
  if (ok) {
    tmp = PyDict_GetItemString(session, "main");
    if (tmp) {
      if (!G->HaveGUI &&
          /* PYMOL-775 added suspend_updates check, but does it make sense? */
          !SettingGetGlobal_b(G, cSetting_suspend_updates) &&
          !partial_restore) {
        ok = MainFromPyList(G, tmp);
      } else if (!quiet) {
        int viewport[2];
        PConvPyListToIntArrayInPlace(tmp, viewport, 2);
        PRINTFB(G, FB_Executive, FB_Actions)
        " Session was saved with: viewport %d, %d\n", viewport[0],
            viewport[1] ENDFB(G);
      }
    }
    if (tmp || (!(partial_restore |
                  partial_session))) { /* ignore missing if partial restore */
      if (PyErr_Occurred()) {
        PyErr_Print();
        ok = false;
      }
      if (!ok) {
        PRINTFB(G, FB_Executive, FB_Errors)
        "ExectiveSetSession-Error: after main.\n" ENDFB(G);
      }
      if (!ok) {
        incomplete = true;
        ok = true; /* keep trying...don't give up */
      }
    }
  }

  if (ok && migrate_sessions) { /* migrate sessions */
    tmp = PyDict_GetItemString(session, "version");
    if (tmp) {
      ok = PConvPyIntToInt(tmp, &version);
      if (ok) {
        ExecutiveMigrateSession(G, version);
      }
    }
  }
  if (ok) {
    if (have_active) ExecutiveSetObjVisib(G, active, true, false);
  }

  G->SettingUnique->old2new.reset();

  if (incomplete) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    "ExectiveSetSession-Warning: restore may be incomplete.\n" ENDFB(G);
  }
  OrthoBackgroundTextureNeedsUpdate(G);
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
  OrthoInvalidateDoDraw(G);
  SceneChanged(G);

  G->Color->HaveOldSessionColors = false;
  G->Color->HaveOldSessionExtColors = false;

  return (ok);
}
