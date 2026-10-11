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

// Getting, setting and unsetting settings on objects, atoms and bonds.

#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/objects/object_volume.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/executive/setting_effects.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"
/**
 * @pre GIL
 */
PyObject* ExecutiveGetSettingOfType(PyMOLGlobals* G, int index,
                                    const char* object, int state, int type) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  pymol::CObject* obj = nullptr;
  CSetting *set_ptr1 = nullptr, *set_ptr2 = nullptr;
  pymol::copyable_ptr<CSetting>* handle = nullptr;

  if (object && object[0]) {
    obj = ExecutiveFindObjectByName(G, object);
    if (!obj)
      return PyErr_Format(P_CmdException, "object \"%s\" not found", object);
    handle = obj->getSettingHandle(-1);
    if (handle) set_ptr1 = handle->get();
    if (state >= 0) {
      handle = obj->getSettingHandle(state);
      if (handle)
        set_ptr2 = handle->get();
      else {
        return PyErr_Format(P_CmdException, "object \"%s\" lacks state %d",
                            object, state + 1);
      }
    }
  }
  {
    switch (type) {
      case cSetting_boolean: {
        int value = SettingGet_b(G, set_ptr2, set_ptr1, index);
        result = PyBool_FromLong(value);
      } break;
      case cSetting_int: {
        int value = SettingGet_i(G, set_ptr2, set_ptr1, index);
        result = Py_BuildValue("i", value);
      } break;
      case cSetting_float: {
        float value = SettingGet_f(G, set_ptr2, set_ptr1, index);
        result = PyFloat_FromDouble(pymol::pretty_f2d(value));
      } break;
      case cSetting_float3: {
        const float* value = SettingGet_3fv(G, set_ptr2, set_ptr1, index);
        if (value) {
          result = Py_BuildValue("fff", pymol::pretty_f2d(value[0]),
                                 pymol::pretty_f2d(value[1]),
                                 pymol::pretty_f2d(value[2]));
        } else {
          PyErr_SetNone(PyExc_ValueError);
        }
      } break;
      case cSetting_color: {
        int value = SettingGet_color(G, set_ptr2, set_ptr1, index);
        result = Py_BuildValue("i", value);
      } break;
      case cSetting_string: {
        OrthoLineType buffer = "";
        result = Py_BuildValue(
            "s", SettingGetTextPtr(G, set_ptr2, set_ptr1, index, buffer));
      } break;
      case cSetting_tuple:
        result = SettingGetTuple(G, set_ptr2, set_ptr1, index);
        break;
      default:
        PyErr_Format(PyExc_ValueError, "invalid setting type %d", type);
        break;
    }
  }
  return (result);
}
/**
 * @pre GIL
 */
PyObject* ExecutiveGetBondSetting(PyMOLGlobals* G, int index, char* s1,
                                  const char* s2, int state, int quiet,
                                  int updates) {
  assert(PyGILState_Check());

  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int sele1, sele2;
  SettingName name;
  PyObject* result = PyList_New(0);

  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);

  if ((sele1 >= 0) && (sele2 >= 0)) {
    while ((ListIterate(I->Spec, rec, next))) {
      if ((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
        obj = (ObjectMolecule*)rec->obj;
        {
          int a, nBond = obj->NBond;
          int nSet = 0;
          const BondType* bi = obj->Bond.data();
          const AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();

          PyObject* pyObjList = nullptr;
          PyObject* pyBondList = nullptr;

          for (a = 0; a < nBond; a++) {
            ai1 = ai + bi->index[0];
            ai2 = ai + bi->index[1];
            if ((SelectorIsMember(G, ai1->selEntry, sele1) &&
                 SelectorIsMember(G, ai2->selEntry, sele2)) ||
                (SelectorIsMember(G, ai2->selEntry, sele1) &&
                 SelectorIsMember(G, ai1->selEntry, sele2))) {
              PyObject* pyBondInfo = PyList_New(3);
              PyObject* bond_setting_value = nullptr;
              if (!pyObjList) {
                pyObjList = PyList_New(2);
                pyBondList = PyList_New(0);
                PyList_SetItem(pyObjList, 0, PyString_FromString(obj->Name));
                PyList_SetItem(pyObjList, 1, pyBondList);
                PyList_Append(result, pyObjList);
                Py_DECREF(pyObjList);
              }
              PyList_SetItem(pyBondInfo, 0,
                             PyInt_FromLong((long)bi->index[0] + 1));
              PyList_SetItem(pyBondInfo, 1,
                             PyInt_FromLong((long)bi->index[1] + 1));
              if (bi->has_setting) {
                bond_setting_value =
                    SettingUniqueGetPyObject(G, bi->unique_id, index);
              }
              PyList_SetItem(pyBondInfo, 2, PConvAutoNone(bond_setting_value));
              PyList_Append(pyBondList, pyBondInfo);
              Py_DECREF(pyBondInfo);
              nSet++;
            }
            bi++;
          }
          if (nSet && !quiet) {
            SettingGetName(G, index, name);
            PRINTF
            " Getting: %s for %d bonds in object \"%s\".\n", name, nSet,
                obj->Name ENDF(G);
          }
        }
      }
    }
  }
  return result;
}
int ExecutiveSetBondSetting(PyMOLGlobals* G, int index, PyObject* tuple,
                            const char* s1, const char* s2, int state,
                            int quiet, int updates) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int sele1, sele2;
  SettingName name = "";
  int unblock;
  int ok = true;
  int side_effects = false;
  union {
    int value_ptr[1];
    float float_storage[1];
  };
  int value_type = 0;

  PRINTFD(G, FB_Executive)
  " %s: entered. '%s' '%s'\n", __func__, s1, s2 ENDFD;
  unblock = PAutoBlock(G);
  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);
  if ((sele1 >= 0) && (sele2 >= 0)) {
    int have_value = false;
    int type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
    PyObject* value = PyTuple_GetItem(tuple, 1);
    if (value) {
      switch (type) {
        case cSetting_boolean:
          *(value_ptr) = PyInt_AsLong(value);
          value_type = cSetting_boolean;
          have_value = true;
          break;
        case cSetting_int:
          *(value_ptr) = PyInt_AsLong(value);
          value_type = cSetting_int;
          have_value = true;
          break;
        case cSetting_float:
          float_storage[0] = PyFloat_AsDouble(value);
          value_type = cSetting_float;
          have_value = true;
          break;
        case cSetting_color: {
          int color_index = ColorGetIndex(G, PyString_AsString(value));
          if ((color_index < 0) && (color_index > cColorExtCutoff))
            color_index = 0;
          *(value_ptr) = color_index;
          value_type = cSetting_color;
          have_value = true;
        } break;
      }
      if (have_value) {
        rec = nullptr;
        while ((ListIterate(I->Spec, rec, next))) {
          if ((rec->type == cExecObject) &&
              (rec->obj->type == cObjectMolecule)) {
            obj = (ObjectMolecule*)rec->obj;
            {
              int a, nBond = obj->NBond;
              int nSet = 0;
              BondType* bi = obj->Bond.data();
              const AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();
              for (a = 0; a < nBond; a++) {
                ai1 = ai + bi->index[0];
                ai2 = ai + bi->index[1];
                if ((SelectorIsMember(G, ai1->selEntry, sele1) &&
                     SelectorIsMember(G, ai2->selEntry, sele2)) ||
                    (SelectorIsMember(G, ai2->selEntry, sele1) &&
                     SelectorIsMember(G, ai1->selEntry, sele2))) {
                  int uid = AtomInfoCheckUniqueBondID(G, bi);
                  bi->has_setting = true;
                  SettingUniqueSetTypedValue(G, uid, index, value_type,
                                             value_ptr);
                  if (updates) side_effects = true;
                  nSet++;
                }
                bi++;
              }
              if (nSet && !quiet) {
                SettingGetName(G, index, name);
                PRINTF
                " Setting: %s set for %d bonds in object \"%s\".\n", name, nSet,
                    obj->Name ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  if (side_effects) {
    SettingGenerateSideEffects(G, index, s1, state,
                               quiet); /* not strictly correct */
    /*    SettingGenerateSideEffects(G,index,s2,state); */
  }

  if (!SettingLevelCheck(G, index, cSettingLevel_bond)) {
    if (!name[0]) SettingGetName(G, index, name);

    PRINTFB(G, FB_Setting, FB_Warnings)
    " Setting-Warning: '%s' is not a bond-level setting\n", name ENDFB(G);
  }

  PAutoUnblock(G, unblock);
  return (ok);
}
/**
 * @return Always true
 */
int ExecutiveUnsetBondSetting(PyMOLGlobals* G, int index, const char* s1,
                              const char* s2, int state, int quiet,
                              int updates) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  SettingName name;
  /*  int unblock; */
  int ok = true;
  int side_effects = false;
  int sele1, sele2;
  PRINTFD(G, FB_Executive)
  " %s: entered. sele '%s' '%s'\n", __func__, s1, s2 ENDFD;
  /* unblock = PAutoBlock(G); */
  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);
  if ((sele1 >= 0) && (sele2 >= 0)) {
    rec = nullptr;
    while ((ListIterate(I->Spec, rec, next))) {
      if ((rec->type == cExecObject) && (rec->obj->type == cObjectMolecule)) {
        obj = (ObjectMolecule*)rec->obj;
        {
          int nSet = 0;
          BondType* bi = obj->Bond.data();
          BondType* bi_end = bi + obj->NBond;
          AtomInfoType *ai1, *ai2, *ai = obj->AtomInfo.data();
          for (; bi != bi_end; ++bi) {
            if (!bi->has_setting) continue;
            ai1 = ai + bi->index[0];
            ai2 = ai + bi->index[1];
            if ((SelectorIsMember(G, ai1->selEntry, sele1) &&
                 SelectorIsMember(G, ai2->selEntry, sele2)) ||
                (SelectorIsMember(G, ai2->selEntry, sele1) &&
                 SelectorIsMember(G, ai1->selEntry, sele2))) {
              int uid = AtomInfoCheckUniqueBondID(G, bi);
              if (!SettingUniqueUnset(G, uid, index)) continue;
              if (updates) side_effects = true;
              nSet++;
            }
          }
          if (nSet && !quiet) {
            SettingGetName(G, index, name);
            PRINTF
            " Setting: %s unset for %d bonds in object \"%s\".\n", name, nSet,
                rec->obj->Name ENDF(G);
          }
        }
      }
    }
  }
  if (side_effects) {
    SettingGenerateSideEffects(G, index, s1, state, quiet);
    /*    SettingGenerateSideEffects(G,index,s2,state); */
  }
  /* PAutoUnblock(G, unblock); */
  return (ok);
}
pymol::Result<> ExecutiveSetSetting(PyMOLGlobals* G, int index, PyObject* tuple,
                                    pymol::zstring_view preSele, int state,
                                    int quiet, int updates) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name = "";
  int nObj = 0;
  const char* sele = preSele.c_str();
  int ok = true;

  pymol::Result<SelectorTmp2> s1;

  if (!preSele.empty()) {
    s1 = SelectorTmp2::make(G, preSele.c_str());
    p_return_if_error(s1);
    sele = s1->getName();
  }

  PRINTFD(G, FB_Executive)
  " %s: entered. sele \"%s\" updates=%d index=%d\n", __func__, sele, updates,
      index ENDFD;

  if (!quiet) {
    SettingGetName(G, index, name);
  }

  pymol::pautoblock unblock(G);

  if ((!sele) || (sele[0] == 0)) { /* global setting */
    ok = SettingSetFromTuple(G, nullptr, index, tuple);
    if (ok) {
      if (!quiet) {
        if (Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, nullptr, nullptr, index, value);
          PRINTF " Setting: %s set to %s.\n", name, value ENDF(G);
        }
      }
      if (updates) {
        SettingGenerateSideEffects(G, index, nullptr, state, quiet);
      }
    }
  } else {
    unsigned char levelmask = 0;
    int side_effects = false;

    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecAll:
            levelmask |= SettingLevelInfo[state < 0 ? cSettingLevel_object
                                                    : cSettingLevel_ostate]
                             .mask;
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                {
                  handle = rec->obj->getSettingHandle(state);
                  if (handle) {
                    SettingCheckHandle(G, *handle);
                    ok = SettingSetFromTuple(G, handle->get(), index, tuple);
                    if (updates) side_effects = true;
                    nObj++;
                  }
                }
              }
            }
            if (Feedback(G, FB_Setting, FB_Actions)) {
              if (nObj && handle) {
                SettingGetTextValue(G, handle->get(), nullptr, index, value);
                if (!quiet) {
                  if (state < 0) {
                    PRINTF
                    " Setting: %s set to %s in %d objects.\n", name, value,
                        nObj ENDF(G);
                  } else {
                    PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n", name,
                        value, nObj, state + 1 ENDF(G);
                  }
                }
              }
            }
            break;
          case cExecSelection:
            if (SettingLevelCheckMask(
                    G, index, SettingLevelInfo[cSettingLevel_bond].mask)) {
              // handle bond-level settings (PYMOL-2726)
              ok = ExecutiveSetBondSetting(G, index, tuple, sele, sele, state,
                                           quiet, false);
              if (updates) side_effects = true;
              sele1 = -1;
            } else {
              levelmask |= SettingLevelInfo[cSettingLevel_atom].mask;
              sele1 = SelectorIndexByName(G, rec->name);
            }

            if (sele1 >= 0) {
              int have_atomic_value = false;
              int type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
              PyObject* value = PyTuple_GetItem(tuple, 1);
              if (value) {
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_SetAtomicSetting;
                op.i1 = index;
                op.ii1 = &op.i3;
                switch (type) {
                  case cSetting_boolean:
                    *(op.ii1) = PyInt_AsLong(value);
                    op.i2 = cSetting_boolean;
                    have_atomic_value = true;
                    break;
                  case cSetting_int:
                    *(op.ii1) = PyInt_AsLong(value);
                    op.i2 = cSetting_int;
                    have_atomic_value = true;
                    break;
                  case cSetting_float:
                    *(float*)op.ii1 = (float)PyFloat_AsDouble(value);
                    op.i2 = cSetting_float;
                    have_atomic_value = true;
                    break;
                  case cSetting_float3: {
                    PConvPyListOrTupleToFloatArrayInPlace(value, op.ttt, 3);
                    op.mat1 = op.ttt;  // for passing (float**)
                    op.ii1 = (int*)&op.mat1;
                    op.i2 = cSetting_float3;
                    have_atomic_value = true;
                  } break;
                  case cSetting_color: {
                    int color_index =
                        ColorGetIndex(G, PyString_AsString(value));
                    if ((color_index < 0) && (color_index > cColorExtCutoff)) {
                      switch (color_index) {
                        case cColorAtomic:
                          color_index = -1;
                          break;
                        case cColorFront:
                        case cColorBack:
                        case cColorDefault:
                          break;
                        default:
                          color_index = 0;
                          break;
                      }
                    }
                    *(op.ii1) = color_index;
                    op.i2 = cSetting_color;
                    have_atomic_value = true;
                  } break;
                }
                if (have_atomic_value) {
                  rec = nullptr;
                  while ((ListIterate(I->Spec, rec, next))) {
                    if ((rec->type == cExecObject) &&
                        (rec->obj->type == cObjectMolecule)) {
                      obj = (ObjectMolecule*)rec->obj;
                      op.i4 = 0;
                      ObjectMoleculeSeleOp(obj, sele1, &op);
                      if (op.i4) {
                        if (updates) side_effects = true;
                        if (!quiet) {
                          PRINTF
                          " Setting: %s set for %d atoms in object \"%s\".\n",
                              name, op.i4, rec->obj->Name ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
            break;
          case cExecObject:
            levelmask |= SettingLevelInfo[state < 0 ? cSettingLevel_object
                                                    : cSettingLevel_ostate]
                             .mask;
            {
              handle = rec->obj->getSettingHandle(state);
              if (handle) {
                SettingCheckHandle(G, *handle);
                ok = SettingSetFromTuple(G, handle->get(), index, tuple);
                if (ok) {
                  if (updates) side_effects = true;
                  if (!quiet) {
                    if (state < 0) { /* object-specific */
                      if (Feedback(G, FB_Setting, FB_Actions)) {
                        SettingGetTextValue(G, handle->get(), nullptr, index,
                                            value);
                        PRINTF
                        " Setting: %s set to %s in object \"%s\".\n", name,
                            value, rec->obj->Name ENDF(G);
                      }
                    } else { /* state-specific */
                      if (Feedback(G, FB_Setting, FB_Actions)) {
                        SettingGetTextValue(G, handle->get(), nullptr, index,
                                            value);
                        PRINTF
                        " Setting: %s set to %s in object \"%s\", state %d.\n",
                            name, value, rec->obj->Name, state + 1 ENDF(G);
                      }
                    }
                  }
                }
              }
            }
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    if (side_effects) SettingGenerateSideEffects(G, index, sele, state, quiet);

    if (!SettingLevelCheckMask(G, index, levelmask)) {
      if (!name[0]) SettingGetName(G, index, name);

      PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: '%s' is a %s-level setting\n", name,
          SettingLevelGetName(index) ENDFB(G);
    }
  }

  if (!ok) {
    return pymol::make_error("Error");
  }
  return {};
}

int ExecutiveSetSettingFromString(PyMOLGlobals* G, int index, const char* value,
                                  const char* sele, int state, int quiet,
                                  int updates) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value2;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name;
  int nObj = 0;
  int ok = true;

  PRINTFD(G, FB_Executive)
  " %s: entered. sele \"%s\"\n", __func__, sele ENDFD;
  if (sele[0] == 0) { /* global setting */
    ok = SettingSetFromString(G, nullptr, index, value);
    if (ok) {
      if (!quiet) {
        if (Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, nullptr, nullptr, index, value2);
          SettingGetName(G, index, name);
          PRINTF " Setting: %s set to %s.\n", name, value2 ENDF(G);
        }
      }
      if (updates) SettingGenerateSideEffects(G, index, sele, state, quiet);
    }
  } else {
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                {
                  handle = rec->obj->getSettingHandle(state);
                  if (handle) {
                    SettingCheckHandle(G, *handle);
                    ok = SettingSetFromString(G, handle->get(), index, value);
                    if (updates)
                      SettingGenerateSideEffects(G, index, rec->name, state,
                                                 quiet);
                    nObj++;
                  }
                }
              }
            }
            if (Feedback(G, FB_Setting, FB_Actions)) {
              if (nObj && handle) {
                SettingGetTextValue(G, handle->get(), nullptr, index, value2);
                SettingGetName(G, index, name);
                if (!quiet) {
                  if (state < 0) {
                    PRINTF
                    " Setting: %s set to %s in %d objects.\n", name, value2,
                        nObj ENDF(G);
                  } else {
                    PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n", name,
                        value2, nObj, state + 1 ENDF(G);
                  }
                }
              }
            }
            break;
          case cExecSelection:
            /* this code has not yet been tested... */

            sele1 = SelectorIndexByName(G, rec->name);
            if (sele1 >= 0) {
              int type;
              int value_store;
              if (SettingStringToTypedValue(G, index, value, &type,
                                            &value_store)) {
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_SetAtomicSetting;
                op.i1 = index;
                op.i2 = type;
                op.ii1 = &value_store;
                rec = nullptr;
                while ((ListIterate(I->Spec, rec, next))) {
                  if ((rec->type == cExecObject) &&
                      (rec->obj->type == cObjectMolecule)) {
                    obj = (ObjectMolecule*)rec->obj;
                    op.i4 = 0;
                    ObjectMoleculeSeleOp(obj, sele1, &op);
                    if (op.i4) {
                      if (updates)
                        SettingGenerateSideEffects(G, index, rec->name, state,
                                                   quiet);
                      if (!quiet) {
                        SettingGetName(G, index, name);
                        PRINTF
                        " Setting: %s set for %d atoms in object \"%s\".\n",
                            name, op.i4, rec->obj->Name ENDF(G);
                      }
                    }
                  }
                }
              }
            }
            break;
          case cExecObject: {
            handle = rec->obj->getSettingHandle(state);
            if (handle) {
              SettingCheckHandle(G, *handle);
              ok = SettingSetFromString(G, handle->get(), index, value);
              if (ok) {
                if (updates)
                  SettingGenerateSideEffects(G, index, sele, state, quiet);
                if (!quiet) {
                  if (state < 0) { /* object-specific */
                    if (Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), nullptr, index,
                                          value2);
                      SettingGetName(G, index, name);
                      PRINTF
                      " Setting: %s set to %s in object \"%s\".\n", name,
                          value2, rec->obj->Name ENDF(G);
                    }
                  } else { /* state-specific */
                    if (Feedback(G, FB_Setting, FB_Actions)) {
                      SettingGetTextValue(G, handle->get(), nullptr, index,
                                          value2);
                      SettingGetName(G, index, name);
                      PRINTF
                      " Setting: %s set to %s in object \"%s\", state %d.\n",
                          name, value2, rec->obj->Name, state + 1 ENDF(G);
                    }
                  }
                }
              }
            }
          } break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return (ok);
}

int ExecutiveSetObjSettingFromString(PyMOLGlobals* G, int index,
                                     const char* value, pymol::CObject* obj,
                                     int state, int quiet, int updates) {
  OrthoLineType value2;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  SettingName name;
  int ok = true;

  PRINTFD(G, FB_Executive)
  " ExecutiveSetObjSettingFromString: entered \n" ENDFD;
  if (!obj) { /* global */
    ok = SettingSetFromString(G, nullptr, index, value);
    if (ok) {
      if (!quiet) {
        if (Feedback(G, FB_Setting, FB_Actions)) {
          SettingGetTextValue(G, nullptr, nullptr, index, value2);
          SettingGetName(G, index, name);
          PRINTF " Setting: %s set to %s.\n", name, value2 ENDF(G);
        }
      }
      if (updates)
        SettingGenerateSideEffects(G, index, obj->Name, state, quiet);
    }
  } else { /* based on a single object */
    {
      handle = obj->getSettingHandle(state);
      if (handle) {
        SettingCheckHandle(G, *handle);
        ok = SettingSetFromString(G, handle->get(), index, value);
        if (ok) {
          if (updates)
            SettingGenerateSideEffects(G, index, obj->Name, state, quiet);
          if (!quiet) {
            if (state < 0) { /* object-specific */
              if (Feedback(G, FB_Setting, FB_Actions)) {
                SettingGetTextValue(G, handle->get(), nullptr, index, value2);
                SettingGetName(G, index, name);
                PRINTF
                " Setting: %s set to %s in object \"%s\".\n", name, value2,
                    obj->Name ENDF(G);
              }
            } else { /* state-specific */
              if (Feedback(G, FB_Setting, FB_Actions)) {
                SettingGetTextValue(G, handle->get(), nullptr, index, value2);
                SettingGetName(G, index, name);
                PRINTF
                " Setting: %s set to %s in object \"%s\", state %d.\n", name,
                    value2, obj->Name, state + 1 ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  return (ok);
}
/**
 * Restore setting to default value.
 *
 * The effective value will be the value at the next higher level (e.g. unset
 * state-level -> use object-level) or the value stored in `G->Default` for
 * global settings.
 *
 * @return Always true
 */
pymol::Result<> ExecutiveUnsetSetting(PyMOLGlobals* G, int index,
                                      pymol::zstring_view preSele, int state,
                                      int quiet, int updates) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int sele1;
  ObjectMoleculeOpRec op;
  pymol::copyable_ptr<CSetting>* handle = nullptr;
  const char* name = SettingGetName(index);
  int nObj = 0;
  int ok = true;
  const char* sele = preSele.c_str();

  pymol::Result<SelectorTmp2> s1;

  if (!preSele.empty()) {
    s1 = SelectorTmp2::make(G, preSele.c_str());
    p_return_if_error(s1);
    sele = s1->getName();
  }

  if (sele[0] == 0) {
    SettingRestoreDefault(G->Setting, index, G->Default);
    if (!quiet && Feedback(G, FB_Executive, FB_Actions)) {
      OrthoLineType value = "";
      SettingGetTextValue(G, nullptr, nullptr, index, value);
      PRINTF " Setting: %s restored to default (%s)\n", name, value ENDF(G);
    }
  } else {
    // Undefine per-object, per-state, or per-atom settings.
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, sele, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecAll:
            rec = nullptr;
            while (ListIterate(I->Spec, rec, next)) {
              if (rec->type == cExecObject) {
                {
                  handle = rec->obj->getSettingHandle(state);
                  if (handle && *handle && SettingUnset(handle->get(), index)) {
                    nObj++;
                  }
                }
              }
            }
            if (Feedback(G, FB_Setting, FB_Actions)) {
              if (nObj && handle) {
                if (!quiet) {
                  if (state < 0) {
                    PRINTF " Setting: %s unset in %d objects.\n", name,
                        nObj ENDF(G);
                  } else {
                    PRINTF
                    " Setting: %s unset in %d objects, state %d.\n", name, nObj,
                        state + 1 ENDF(G);
                  }
                }
              }
            }
            break;
          case cExecSelection:
            if (SettingLevelCheckMask(
                    G, index, SettingLevelInfo[cSettingLevel_bond].mask)) {
              // handle bond-level settings (PYMOL-2726)
              ok = ExecutiveUnsetBondSetting(G, index, sele, sele, state, quiet,
                                             false);
              sele1 = -1;
            } else {
              sele1 = SelectorIndexByName(G, rec->name);
            }

            if (sele1 >= 0) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_SetAtomicSetting;
              op.i1 = index;
              op.i2 = cSetting_blank;
              op.ii1 = nullptr;

              rec = nullptr;
              while ((ListIterate(I->Spec, rec, next))) {
                if ((rec->type == cExecObject) &&
                    (rec->obj->type == cObjectMolecule)) {
                  obj = (ObjectMolecule*)rec->obj;
                  op.i4 = 0;
                  ObjectMoleculeSeleOp(obj, sele1, &op);
                  if (op.i4) {
                    if (!quiet) {
                      PRINTF
                      " Setting: %s unset for %d atoms in object \"%s\".\n",
                          name, op.i4, rec->obj->Name ENDF(G);
                    }
                  }
                }
              }
            }
            break;
          case cExecObject: {
            handle = rec->obj->getSettingHandle(state);
            if (handle && *handle && SettingUnset(handle->get(), index)) {
              if (!quiet) {
                if (state < 0) { /* object-specific */
                  if (Feedback(G, FB_Setting, FB_Actions)) {
                    PRINTF
                    " Setting: %s unset in object \"%s\".\n", name,
                        rec->obj->Name ENDF(G);
                  }
                } else { /* state-specific */
                  if (Feedback(G, FB_Setting, FB_Actions)) {
                    PRINTF
                    " Setting: %s unset in object \"%s\", state %d.\n", name,
                        rec->obj->Name, state + 1 ENDF(G);
                  }
                }
              }
            }
          } break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  if (updates) SettingGenerateSideEffects(G, index, sele, state, quiet);
  if (!ok) {
    return pymol::make_error("Error");
  }
  return {};
}
