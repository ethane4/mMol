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

// Maps, isosurfaces, isomeshes, volumes, ramps and slices.

#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list_macros.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/objects/object_mesh.h"
#include "pymol/objects/object_slice.h"
#include "pymol/objects/object_surface.h"
#include "pymol/objects/object_volume.h"
#include "pymol/base/parse.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

/**
 * Does the object and state validation which is common to isosurface, isomesh,
 * and volume creation.
 */
template <typename ObjectT>
static pymol::Result<> EtcHelper(PyMOLGlobals* G,  //
                                 const char* obj_name, int& obj_state,
                                 ObjectT*& origObj,  //
                                 const char* map_name, int& map_state,
                                 ObjectMap*& mapObj,  //
                                 int& multi) {
  if (obj_state < -3) {
    return pymol::make_error("Invalid state ", obj_state + 1);
  }

  if (map_state < -3) {
    return pymol::make_error("Invalid source_state ", map_state + 1);
  }

  mapObj = ExecutiveFindObject<ObjectMap>(G, map_name);
  if (!mapObj) {
    return pymol::make_error("Map object \"", map_name, "\" not found");
  }

  // Fail if object with same name but different type exists
  auto anyOrigObj = ExecutiveFindObjectByName(G, obj_name);
  origObj = dynamic_cast<ObjectT*>(anyOrigObj);
  if (anyOrigObj && !origObj) {
    return pymol::make_error(
        "Object name ", obj_name,
        " already exists as another type. Cannot overwrite.");
  }

  switch (obj_state) {
    case -1:
      // all states
      obj_state = 0;
      map_state = -1;
      break;
    case -2:
      // current state
      obj_state = SceneGetState(G);
      if (map_state < 0) map_state = obj_state;
      break;
    case -3:
      // append mode
      obj_state = origObj ? origObj->getNFrame() : 0;
      if (map_state < 0) map_state = obj_state;
      break;
  }

  switch (map_state) {
    case -1:
      // all states
      map_state = 0;
      multi = true;
      break;
    case -2:
      // current state
      map_state = SceneGetState(G);
      break;
    case -3:
      // append mode
      map_state = mapObj->getNFrame() - 1;
      break;
  }

  return {};
}

pymol::Result<> ExecutiveIsosurfaceEtc(PyMOLGlobals* G, const char* surf_name,
                                       const char* map_name, float lvl,
                                       const char* sele, float fbuf, int state,
                                       float carve, int map_state, int side,
                                       int quiet, int surf_mode) {
  int c;
  ObjectSurface *obj = nullptr, *origObj = nullptr;
  ObjectMap* mapObj;
  float mn[3] = {0, 0, 0};
  float mx[3] = {15, 15, 15};
  pymol::vla<float> vert_vla;
  ObjectMapState* ms;
  int multi = false;

  auto res0 = EtcHelper(G, surf_name, state, origObj, map_name, map_state,
                        mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    while (1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if (ms) {
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
          case 0: /* using map to get extents */
            for (c = 0; c < 3; c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3 * 7 + c];
            }
            if (!ms->Matrix.empty()) {
              transform44d3f(ms->Matrix.data(), mn, mn);
              transform44d3f(ms->Matrix.data(), mx, mx);
              {
                float tmp;
                int a;
                for (a = 0; a < 3; a++)
                  if (mn[a] > mx[a]) {
                    tmp = mn[a];
                    mn[a] = mx[a];
                    mx[a] = tmp;
                  }
              }
            }
            carve = 0.0F;
            break;
          case 1: /* using selection to get extents */
          {
            auto tmpsele = SelectorTmp2::make(G, sele);
            p_return_if_error(tmpsele);
            auto s1 = tmpsele->getName();
            ExecutiveGetExtent(G, s1, mn, mx, false, -1,
                               false); /* TODO state */
            if (carve != 0.0F) {
              vert_vla = pymol::vla_take_ownership(
                  ExecutiveGetVertexVLA(G, s1, state));
              if (fbuf <= R_SMALL4) fbuf = fabs(carve);
            }
          }
            for (c = 0; c < 3; c++) {
              mn[c] -= fbuf;
              mx[c] += fbuf;
            }
            break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
        " Isosurface: buffer %8.3f carve %8.3f\n", fbuf, carve ENDFB(G);
        obj = ObjectSurfaceFromBox(G, origObj, mapObj, map_state, state, mn, mx,
                                   lvl, static_cast<cIsosurfaceMode>(surf_mode),
                                   carve, std::move(vert_vla),
                                   static_cast<cIsosurfaceSide>(side), quiet);
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);

        if (!origObj) {
          ObjectSetName(obj, surf_name);
          ExecutiveManageObject(G, obj, -1, quiet);
        }
        if (SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if (obj) ObjectGotoState(obj, state);
        if (!quiet) {
          PRINTFB(G, FB_ObjectSurface, FB_Actions)
          " Isosurface: created \"%s\", setting level to %5.3f\n", surf_name,
              lvl ENDFB(G);
        }
      } else if (!multi) {
        return pymol::make_error("state ", map_state + 1,
                                 " not present in map \"", map_name, "\"");
      }
      if (multi) {
        origObj = obj;
        map_state++;
        state++;
        if (map_state >= mapObj->State.size()) break;
      } else {
        break;
      }
    }
  }
  return {};
}

pymol::Result<> ExecutiveIsomeshEtc(PyMOLGlobals* G, const char* mesh_name,
                                    const char* map_name, float lvl,
                                    const char* sele, float fbuf, int state,
                                    float carve, int map_state, int quiet,
                                    int mesh_mode_, float alt_lvl) {
  auto const mesh_mode = static_cast<cIsomeshMode>(mesh_mode_);
  ObjectMesh *obj = nullptr, *origObj = nullptr;
  ObjectMap* mapObj;
  float mn[3] = {0, 0, 0};
  float mx[3] = {15, 15, 15};
  float* vert_vla = nullptr;
  int multi = false;
  ObjectMapState* ms;
  ObjectMolecule* sele_obj = nullptr;
  CSymmetry* symm;

  auto res0 = EtcHelper(G, mesh_name, state, origObj, map_name, map_state,
                        mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    while (1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if (ms) {
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
          case 0: /* do the whole map */
          {
            int c;
            for (c = 0; c < 3; c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3 * 7 + c];
            }
          }
            if (!ms->Matrix.empty()) {
              transform44d3f(ms->Matrix.data(), mn, mn);
              transform44d3f(ms->Matrix.data(), mx, mx);
              {
                float tmp;
                int a;
                for (a = 0; a < 3; a++)
                  if (mn[a] > mx[a]) {
                    tmp = mn[a];
                    mn[a] = mx[a];
                    mx[a] = tmp;
                  }
              }
            }
            carve = -0.0; /* impossible */
            break;
          case 1: /* just do area around selection */
                  /* determine the selected object */
            {
              auto tmpsele = SelectorTmp2::make(G, sele);
              p_return_if_error(tmpsele);
              int sele1 = tmpsele->getIndex();
              if (sele1 >= 0)
                sele_obj = SelectorGetSingleObjectMolecule(G, sele1);
              auto s1 = tmpsele->getName();
              ExecutiveGetExtent(G, s1, mn, mx, false, -1,
                                 false); /* TODO state */
              if (carve != 0.0) {
                vert_vla = ExecutiveGetVertexVLA(G, s1, state);
                if (fbuf <= R_SMALL4) fbuf = fabs(carve);
              }
            }
            {
              int c;
              for (c = 0; c < 3; c++) {
                mn[c] -= fbuf;
                mx[c] += fbuf;
              }
            }
            break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
        " Isomesh: buffer %8.3f carve %8.3f \n", fbuf, carve ENDFB(G);

        symm = nullptr;
        if (sele_obj && ObjectMapValidXtal(mapObj, state)) {
          if (SettingGet_b(G, nullptr, sele_obj->Setting.get(),
                           cSetting_map_auto_expand_sym) &&
              (sele_obj->Symmetry)) {
            // legacy default: take symmetry from molecular object
            symm = sele_obj->Symmetry.get();
          } else if (SettingGet_b(G, nullptr, mapObj->Setting.get(),
                                  cSetting_map_auto_expand_sym)) {
            // fallback: take symmetry from map state
            symm = ms->Symmetry.get();
          }
        }

        if (symm) {
          obj = ObjectMeshFromXtalSym(G, origObj, mapObj, symm, map_state,
                                      state, mn, mx, lvl, mesh_mode, carve,
                                      vert_vla, alt_lvl, quiet);
        } else {
          obj = nullptr;
        }
        if (!obj) {
          obj = ObjectMeshFromBox(G, origObj, mapObj, map_state, state, mn, mx,
                                  lvl, mesh_mode, carve, vert_vla, alt_lvl,
                                  quiet);
        }
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);

        if (!origObj) {
          ObjectSetName(obj, mesh_name);
          ExecutiveManageObject(G, obj, false, quiet);
        }

        if (SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if (obj) ObjectGotoState(obj, state);
        if (!quiet) {
          if (mesh_mode != cIsomeshMode::gradient) {
            PRINTFB(G, FB_ObjectMesh, FB_Actions)
            " Isomesh: created \"%s\", setting level to %5.3f\n", mesh_name,
                lvl ENDFB(G);
          } else {
            PRINTFB(G, FB_ObjectMesh, FB_Actions)
            " Gradient: created \"%s\"\n", mesh_name ENDFB(G);
          }
        }
      } else if (!multi) {
        return pymol::make_error("state ", map_state + 1,
                                 " not present in map \"", map_name, "\"");
      }
      if (multi) {
        origObj = obj;
        map_state++;
        state++;
        if (map_state >= mapObj->State.size()) break;
      } else {
        break;
      }
    }
  }
  return {};
}

pymol::Result<> ExecutiveVolume(PyMOLGlobals* G, const char* volume_name,
                                const char* map_name, float lvl,
                                const char* sele, float fbuf, int state,
                                float carve, int map_state, int quiet) {
  ObjectVolume *obj = nullptr, *origObj = nullptr;
  ObjectMap* mapObj;
  float mn[3] = {0, 0, 0};
  float mx[3] = {15, 15, 15};
  float* vert_vla = nullptr;
  int multi = false;
  ObjectMapState* ms;
  ObjectMolecule* sele_obj = nullptr;
  CSymmetry* symm;

  /* Goal: make a new volume map from the object or selection
   *
   * If the user specifies an VOLUME OBJECT NAME that already exists, then
   * check to make sure it's a volume.  If it is, keep it, else abort
   *
   * If the user specifies a MAP OBJECT NAME that already exists, then
   * make sure it's really a map.  Keep it if it is (to append the state)
   * otherwise, delete it.
   *
   * Given that the MAP is valid, find the current states for the MAP object
   * and the scene/protein.
   *
   * If the user is loading a map with N-states, into an N-state object, we need
   * to loop over each state.  So,
   *
   * FOR EACH state in the MAP:
   *   - determine the extents of the MAP for THIS STATE
   *
   *   - create the VOLUME object from the MAP/XTAL SYMM
   *   - if VOLUME CREATION didn't work, then try creating from a BOX
   *
   *   - transform (translate/rotate) the VOLUME vis so that it
   *     matches the map's ObjectMatrix.
   *
   *   - if this is a new VOLUME, (eg., origObj==nullptr) then set the internal
   *     PyMOL state to manage the object and add the name to the object list
   *
   *   - update the map state for the next loop and re-iterate
   */

  auto res0 = EtcHelper(G, volume_name, state, origObj, map_name, map_state,
                        mapObj, multi);
  if (!res0) {
    return res0.error();
  }

  {
    /* do for each state */
    while (1) {
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if (ms) {
        /* determine extents */
        int box_mode = (sele && sele[0]) ? 1 : 0;
        switch (box_mode) {
          case 0: /* do the whole map */
          {
            int c;
            for (c = 0; c < 3; c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3 * 7 + c];
            }
          }
            if (!ms->Matrix.empty()) {
              transform44d3f(ms->Matrix.data(), mn, mn);
              transform44d3f(ms->Matrix.data(), mx, mx);
              {
                float tmp;
                int a;
                for (a = 0; a < 3; a++)
                  if (mn[a] > mx[a]) {
                    tmp = mn[a];
                    mn[a] = mx[a];
                    mx[a] = tmp;
                  }
              }
            }
            carve = -0.0; /* impossible */
            break;
          case 1: /* just do area around selection */
          {
            auto tmpsele = SelectorTmp2::make(G, sele);
            p_return_if_error(tmpsele);
            int sele1 = tmpsele->getIndex();
            if (sele1 >= 0)
              sele_obj = SelectorGetSingleObjectMolecule(G, sele1);
            auto const s1 = tmpsele->getName();
            ExecutiveGetExtent(G, s1, mn, mx, false, -1,
                               false); /* TODO state */
            if (carve != 0.0) {
              vert_vla = ExecutiveGetVertexVLA(G, s1, state);
              if (fbuf <= R_SMALL4) fbuf = fabs(carve);
            }
          }
            {
              int c;
              for (c = 0; c < 3; c++) {
                mn[c] -= fbuf;
                mx[c] += fbuf;
              }
            }
            break;
        }
        PRINTFB(G, FB_CCmd, FB_Blather)
        " Volume: buffer %8.3f carve %8.3f \n", fbuf, carve ENDFB(G);

        symm = nullptr;
        if (sele_obj && ObjectMapValidXtal(mapObj, state)) {
          if (SettingGet_b(G, nullptr, sele_obj->Setting.get(),
                           cSetting_map_auto_expand_sym) &&
              (sele_obj->Symmetry)) {
            // legacy default: take symmetry from molecular object
            symm = sele_obj->Symmetry.get();
          } else if (SettingGet_b(G, nullptr, mapObj->Setting.get(),
                                  cSetting_map_auto_expand_sym)) {
            // fallback: take symmetry from map state
            symm = ms->Symmetry.get();
          }
        }

        if (symm) {
          obj = ObjectVolumeFromXtalSym(G, origObj, mapObj, symm, map_state,
                                        state, mn, mx, lvl, box_mode, carve,
                                        vert_vla, quiet);
        } else {
          obj = nullptr;
        }

        if (!obj) {
          obj = ObjectVolumeFromBox(G, origObj, mapObj, map_state, state, mn,
                                    mx, lvl, box_mode, carve, vert_vla, quiet);
        }
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, mapObj, obj, 1, 1, -1, -1, false, 0, quiet);
        /* set the object name
         * manage the object in the UI */
        if (!origObj) {
          ObjectSetName(obj, volume_name);
          ExecutiveManageObject(G, obj, false, quiet);
        }

        if (SettingGetGlobal_b(G, cSetting_isomesh_auto_state))
          if (obj) ObjectGotoState(obj, state);
        if (!quiet) {
          PRINTFB(G, FB_ObjectVolume, FB_Actions)
          " Volume: created \"%s\"\n", volume_name ENDFB(G);
        }
      } else if (!multi) {
        return pymol::make_error("state ", map_state + 1,
                                 " not present in map \"", map_name, "\"");
      }
      if (multi) {
        origObj = obj;
        map_state++;
        state++;
        if (map_state >= mapObj->State.size()) break;
      } else {
        break;
      }
    }
  }
  return {};
}

void ExecutiveInvalidateMapDependents(PyMOLGlobals* G, const char* map_name,
                                      const char* new_name) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      switch (rec->obj->type) {
        case cObjectMesh:
          ObjectMeshInvalidateMapName((ObjectMesh*)rec->obj, map_name,
                                      new_name);
          break;
        case cObjectSurface:
          ObjectSurfaceInvalidateMapName((ObjectSurface*)rec->obj, map_name,
                                         new_name);
          break;
        case cObjectVolume:
          ObjectVolumeInvalidateMapName((ObjectVolume*)rec->obj, map_name,
                                        new_name);
          break;
      }
    }
  }
  SceneInvalidate(G);
}

/**
 * returns a pointer to the data in a volume or map object
 */
CField* ExecutiveGetVolumeField(PyMOLGlobals* G, const char* objName,
                                int state) {
  ObjectMapState* oms;
  pymol::CObject* obj;

  obj = ExecutiveFindObjectByName(G, objName);
  ok_assert(1, obj);

  switch (obj->type) {
    case cObjectVolume:
      return ObjectVolumeGetField((ObjectVolume*)obj);
    case cObjectMap:
      oms = ObjectMapGetState((ObjectMap*)obj, state);
      ok_assert(1, oms && oms->Field);
      return oms->Field->data.get();
  }

ok_except1:
  return nullptr;
}

/**
 * returns allocated memory
 */
pymol::Result<std::vector<float>> ExecutiveGetHistogram(PyMOLGlobals* G,
                                                        const char* objName,
                                                        int n_points,
                                                        float min_val,
                                                        float max_val) {
  pymol::CObject* obj;
  ObjectMapState* oms = nullptr;

  obj = ExecutiveFindObjectByName(G, objName);

  if (!obj) {
    return pymol::make_error("could not find object ", objName);
  }

  switch (obj->type) {
    case cObjectMap:
      oms = ObjectMapGetState((ObjectMap*)obj, 0);
      break;
    case cObjectVolume:
      oms = ObjectVolumeGetMapState((ObjectVolume*)obj);
      break;
    default:
      return pymol::make_error("object type must be map or volume");
  }

  if (oms) {
    auto hist = std::vector<float>(n_points + 4);
    float range = SettingGet_f(G, obj->Setting.get(), nullptr,
                               cSetting_volume_data_range);
    ObjectMapStateGetHistogram(G, oms, n_points, range, hist.data(), min_val,
                               max_val);
    return hist;
  }

  return pymol::make_error("failed to get map state");
}

PyObject* ExecutiveGetVolumeRamp(PyMOLGlobals* G, const char* objName,
                                 int state) {
  PyObject* result = nullptr;

  PRINTFD(G, FB_Executive) "Executive-GetVolumeRamp Entered.\n" ENDFD;

  if (auto obj = ExecutiveFindObject<ObjectVolume>(G, objName)) {
    result = ObjectVolumeGetRamp(obj, state);
  }

  PRINTFD(G, FB_Executive) "Executive-GetVolumeRamp Exited.\n" ENDFD;

  return result;
}

pymol::Result<> ExecutiveSetVolumeRamp(PyMOLGlobals* G, const char* objName,
                                       std::vector<float> ramp_list,
                                       int state) {
  auto obj = ExecutiveFindObject<ObjectVolume>(G, objName);
  if (obj) {
    return ObjectVolumeSetRamp(obj, std::move(ramp_list), state);
  }

  return pymol::make_error("Object ", objName, " not found");
}

pymol::Result<> ExecutiveIsolevel(PyMOLGlobals* G, const char* name,
                                  float level, int state, int quiet) {
  auto obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    switch (obj->type) {
      case cObjectMesh:
        ObjectMeshSetLevel((ObjectMesh*)obj, level, state, quiet);
        SceneChanged(G);
        return {};
      case cObjectSurface:
        ObjectSurfaceSetLevel((ObjectSurface*)obj, level, state, quiet);
        SceneChanged(G);
        return {};
      default:
        return pymol::make_error("Object ", name, " is of wrong type.");
    }
  }
  return pymol::make_error("Object not found");
}

pymol::Result<float> ExecutiveGetIsolevel(PyMOLGlobals* G, const char* name,
                                          int state) {
  auto obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    switch (obj->type) {
      case cObjectMesh:
        return ObjectMeshGetLevel((ObjectMesh*)obj, state);
      case cObjectSurface:
        return ObjectSurfaceGetLevel((ObjectSurface*)obj, state);
      default:
        return pymol::make_error("Object ", name, " is of wrong type.");
    }
  }
  return pymol::make_error("Object not found");
}

pymol::Result<> ExecutiveRampNew(PyMOLGlobals* G, const char* name,
                                 const char* src_name, pymol::vla<float> range,
                                 pymol::vla<float> color, int src_state,
                                 const char* sele, float beyond, float within,
                                 float sigma, int zero, int calc_mode,
                                 int quiet) {
  ObjectGadgetRamp* obj = nullptr;
  ObjectGadgetRamp* origRamp = nullptr;
  pymol::CObject* src_obj = nullptr;
  pymol::CObject* origObj = ExecutiveFindObjectByName(G, name);
  float* vert_vla = nullptr;
  int rampType = -1;

  if (origObj && origObj->type == cObjectGadget &&
      ((ObjectGadget*)origObj)->GadgetType == cGadgetRamp) {
    origRamp = (ObjectGadgetRamp*)origObj;
    rampType = origRamp->RampType;
  } else if (!range || !(color || calc_mode)) {
    return pymol::make_error("Missing 'range' or 'color' to create new ramp.");
  }

  if (src_name && src_name[0]) {
    if (WordMatchExact(G, src_name, cKeywordNone, true)) {
      rampType = cRampNone;
    } else {
      src_obj = ExecutiveFindObjectByName(G, src_name);
      if (src_obj) {
        switch (src_obj->type) {
          case cObjectMap:
            rampType = cRampMap;
            break;
          case cObjectMolecule:
            rampType = cRampMol;
            break;
          default:
            pymol::make_error(src_name, " is not a map or molecule.");
        }
      } else {
        return pymol::make_error(src_name, " not found.");
      }
    }
  }

  switch (rampType) {
    case cRampMap:
      /* mapping this ramp from a selection */
      if (sele && sele[0]) {
        auto tmpsele = SelectorTmp::make(G, sele);
        p_return_if_error(tmpsele);
        sele = tmpsele->getName();
        assert(sele[0]);

        vert_vla = ExecutiveGetVertexVLA(G, sele, src_state);
      }
      obj = ObjectGadgetRampMapNewAsDefined(
          G, origRamp, (ObjectMap*)src_obj, std::move(range), std::move(color),
          src_state, vert_vla, beyond, within, sigma, zero, calc_mode);
      VLAFreeP(vert_vla);
      break;
    case cRampNone:
    case cRampMol:
      obj = ObjectGadgetRampMolNewAsDefined(
          G, origRamp, (ObjectMolecule*)src_obj, std::move(range),
          std::move(color), src_state, calc_mode);
      break;
    default:
      return pymol::make_error("Missing 'name' to create new ramp.");
  }

  if (!obj) return pymol::make_error("Object not found");

  if (obj != origRamp) {
    ExecutiveDelete(G, name);
    ObjectSetName(obj, name);
    ColorRegisterExt(G, obj->Name, obj);
    ExecutiveManageObject(G, obj, false, quiet);
  }

  ExecutiveInvalidateRep(G, cKeywordAll, cRepAll,
                         cRepInvColor); /* recolor everything */
  return {};
}
#define cMapOperatorMinimum 0
#define cMapOperatorMaximum 1
#define cMapOperatorSum 2
#define cMapOperatorAverage 3
#define cMapOperatorDifference 4
#define cMapOperatorCopy 5
#define cMapOperatorUnique 6

pymol::Result<> ExecutiveMapSet(PyMOLGlobals* G, const char* name,
                                int operator_, const char* operands,
                                int target_state, int source_state, int zoom,
                                int quiet) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  int isNew = false;
  ObjectMap* target = ExecutiveFindObjectMapByName(G, name);
  ObjectMap* first_operand = nullptr;
  int src_state_start = 0, src_state_stop = 0;
  int list_id = ExecutiveGetNamesListFromPattern(G, operands, true, true);

  if (target_state < 0) /* if we're targeting all states, 0 is the offset */
    target_state = 0;

  /* first, figure out what the range of input states is */

  if (source_state < 0) { /* all source states */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int max_n_state = 0;
    SpecRec* rec;
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecObject:
            if (rec->obj->type == cObjectMap) {
              ObjectMap* obj = (ObjectMap*)rec->obj;
              if (obj->State.size() > max_n_state)
                max_n_state = obj->State.size(); /* count states */
            }
        }
      }
    }
    TrackerDelIter(I_Tracker, iter_id);
    src_state_start = 0;
    src_state_stop = max_n_state;
  } else {
    src_state_start = source_state;
    src_state_stop = source_state + 1;
  }

  {
    /* next, find the first operand */

    OrthoLineType first_op_st;
    ParseWordCopy(first_op_st, operands,
                  sizeof(OrthoLineType) - 1); /* copy the first operand */
    {
      int sub_list_id =
          ExecutiveGetNamesListFromPattern(G, first_op_st, true, true);
      int sub_iter_id = TrackerNewIter(I_Tracker, 0, sub_list_id);
      SpecRec* rec;

      while (TrackerIterNextCandInList(I_Tracker, sub_iter_id,
                                       (TrackerRef**)(void*)&rec)) {
        if (rec) {
          switch (rec->type) {
            case cExecObject:
              if (rec->obj->type == cObjectMap) {
                ObjectMap* obj = (ObjectMap*)rec->obj;
                first_operand = obj;
              }
              break;
          }
        }
        if (first_operand) break;
      }
      TrackerDelList(I_Tracker, sub_list_id);
      TrackerDelIter(I_Tracker, sub_iter_id);
    }
  }

  {
    /* okay, next thing we need to worry about is where we're putting the data.

       Case 1. If the map already exists, then we'll use the existing map points
       for storing the result.

       Case 2. If the operation implies a copy of existing map geometry, then
       we'll create that copy first before performing the calulation.

       Case 3. If the operation implies a new map geometry, then we need to
       compute that geometry and create the map.

     */

    if (!target) { /* target map doesn't exist... */
      int need_union_geometry = false;
      int need_first_geometry = false;
      switch (operator_) {
        case cMapOperatorSum:
        case cMapOperatorAverage:
        case cMapOperatorMinimum:
        case cMapOperatorMaximum:
        case cMapOperatorDifference:
          need_union_geometry = true;
          break;
        case cMapOperatorUnique:
        case cMapOperatorCopy:
          need_first_geometry = true;
          break;
      }

      if (need_union_geometry) {
        int src_state, trg_state;
        ObjectMapDesc desc;
        target = new ObjectMap(G);

        ObjectSetName(target, name);
        isNew = true;

        for (src_state = src_state_start; src_state < src_state_stop;
             src_state++) {
          trg_state = src_state + target_state;
          desc.mode = cObjectMap_OrthoMinMaxGrid; /* Orthorhombic: min, max,
                                                     spacing,
                                                     centered over range  */
          desc.init_mode = 0;                     /* zeros */
          desc.Grid[0] = 1.0F;
          desc.Grid[1] = 1.0F;
          desc.Grid[2] = 1.0F;

          {
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            float grid_sum[3] = {0.0F, 0.0F, 0.0F};
            int grid_count = 0;
            int first_extent = true;

            SpecRec* rec;
            while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                             (TrackerRef**)(void*)&rec)) {
              if (rec) {
                /* compute an average grid and get the effective "union" extent
                 */
                switch (rec->type) {
                  case cExecObject:
                    if (rec->obj->type == cObjectMap) {
                      ObjectMap* obj = (ObjectMap*)rec->obj;
                      ObjectMapState* ms = &obj->State[src_state];
                      if (src_state < obj->State.size()) {
                        if (ms->Active) {
                          if (first_extent) {
                            copy3f(ms->ExtentMin, desc.MinCorner);
                            copy3f(ms->ExtentMax, desc.MaxCorner);
                            first_extent = false;
                          } else {
                            int b;
                            for (b = 0; b < 3; b++) {
                              if (ms->ExtentMin[b] < desc.MinCorner[b])
                                desc.MinCorner[b] = ms->ExtentMin[b];
                              if (ms->ExtentMax[b] > desc.MaxCorner[b])
                                desc.MaxCorner[b] = ms->ExtentMax[b];
                            }
                          }
                          if (!ObjectMapStateValidXtal(ms)) {
                            /* other general-purpose maps currently handled */
                            int b;
                            for (b = 0; b < 3; b++) {
                              grid_sum[b] += ms->Grid[b];
                            }
                            grid_count++;
                          }
                        }
                      }
                    }
                }
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
            if (grid_count) {
              int b;
              for (b = 0; b < 3; b++) {
                desc.Grid[b] = grid_sum[b] / grid_count;
              }
            }
            if (!first_extent) {
              add3f(desc.Grid, desc.MaxCorner, desc.MaxCorner);
              subtract3f(desc.MinCorner, desc.Grid, desc.MinCorner);
              ObjectMapNewStateFromDesc(G, target, &desc, trg_state, quiet);
              target->State[trg_state].Active = true;
            }
          }
        }
        /* need union geometry */
      } else if (need_first_geometry) {
        if (first_operand) {
          if (ObjectMapNewCopy(G, first_operand, &target, source_state,
                               target_state)) {
            if (target) {
              ObjectSetName(target, name);
              isNew = true;
            }
          }
        }
      }
    }
  }

  if (!target) {
    TrackerDelList(I_Tracker, list_id);
    return pymol::make_error("Cannot find or construct target map.");
  }

  /* now do the actual operation */

  int src_state;
  for (src_state = src_state_start; src_state < src_state_stop; src_state++) {
    int trg_state = src_state + target_state;
    ObjectMapState* ms;
    VecCheckEmplace(target->State, trg_state, G);

    ms = &target->State[target_state];
    if (ms->Active) {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      int n_pnt =
          (ms->Field->points->size() / ms->Field->points->base_size) / 3;
      float* pnt = (float*)ms->Field->points->data.data();
      float* r_value = pymol::malloc<float>(n_pnt);
      float* l_value = pymol::calloc<float>(n_pnt);
      int* present = pymol::calloc<int>(n_pnt);
      std::vector<std::uint8_t> inside(n_pnt);
      SpecRec* rec;

      while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                       (TrackerRef**)(void*)&rec)) {
        if (rec) {
          if (rec->type == cExecObject) {
            if (rec->obj->type == cObjectMap) {
              ObjectMap* obj = (ObjectMap*)rec->obj;
              if (ObjectMapInterpolate(obj, src_state, pnt, r_value,
                                       inside.data(), n_pnt)) {
                int a;
                float* rv = r_value;
                float* lv = l_value;
                auto* flg = inside.data();
                int* pre = present;

                switch (operator_) {
                  case cMapOperatorCopy:
                    for (a = 0; a < n_pnt; a++) {
                      if (flg) {
                        *lv = *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                    break;
                  case cMapOperatorMinimum:
                    for (a = 0; a < n_pnt; a++) {
                      if (flg) {
                        if (*pre) {
                          if (*lv > *rv) *lv = *rv;
                        } else { /* first map */
                          *pre = 1;
                          *lv = *rv;
                        }
                      }
                      rv++;
                      lv++;
                      flg++;
                      pre++;
                    }
                    break;
                  case cMapOperatorMaximum:
                    for (a = 0; a < n_pnt; a++) {
                      if (flg) {
                        if (*pre) {
                          if (*lv < *rv) *lv = *rv;
                        } else { /* first map */
                          *pre = 1;
                          *lv = *rv;
                        }
                      }
                      rv++;
                      lv++;
                      flg++;
                      pre++;
                    }
                    break;
                  case cMapOperatorSum:
                    for (a = 0; a < n_pnt; a++) {
                      if (flg) {
                        *lv += *rv;
                      }
                      rv++;
                      lv++;
                      flg++;
                    }
                    break;
                  case cMapOperatorAverage:
                    for (a = 0; a < n_pnt; a++) {
                      if (flg) {
                        *lv += *rv;
                      }
                      (*pre)++;
                      rv++;
                      lv++;
                      flg++;
                      pre++;
                    }
                    break;
                  case cMapOperatorDifference:
                    if (obj != first_operand) {
                      for (a = 0; a < n_pnt; a++) {
                        if (flg) {
                          *lv -= *rv;
                        }
                        rv++;
                        lv++;
                        flg++;
                      }
                    } else {
                      for (a = 0; a < n_pnt; a++) {
                        if (flg) {
                          *lv += *rv;
                        }
                        rv++;
                        lv++;
                        flg++;
                      }
                    }
                    break;
                  case cMapOperatorUnique:
                    if (obj != first_operand) {
                      for (a = 0; a < n_pnt; a++) {
                        if (flg) {
                          *lv -= *rv;
                        }
                        rv++;
                        lv++;
                        flg++;
                      }
                    } else {
                      for (a = 0; a < n_pnt; a++) {
                        if (flg) {
                          *lv += *rv;
                        }
                        rv++;
                        lv++;
                        flg++;
                      }
                    }

                    break;
                }
              }
            }
          }
        }
      }

      {
        int a;
        float* lv = l_value;
        int* pre = present;

        switch (operator_) {
          case cMapOperatorUnique:
            lv = l_value;
            for (a = 0; a < n_pnt; a++) {
              if (*lv < 0.0F) *lv = 0.0F;
              lv++;
            }
            break;
          case cMapOperatorAverage:
            lv = l_value;
            pre = present;
            for (a = 0; a < n_pnt; a++) {
              if (*pre) *lv /= *pre;
              lv++;
              pre++;
            }
        }
      }

      /* copy after calculation so that operand can include target */

      memcpy(ms->Field->data->data.data(), l_value, n_pnt * sizeof(float));

      FreeP(present);
      FreeP(l_value);
      FreeP(r_value);
      TrackerDelIter(I_Tracker, iter_id);
    }
  }

  /* and finally, update */

  if (target) {
    ObjectMapUpdateExtents(target);
    if (isNew) {
      ExecutiveManageObject(G, target, -1, quiet);
    } else {
      ExecutiveDoZoom(G, target, false, zoom, true);
    }
    SceneChanged(G);
  }
  TrackerDelList(I_Tracker, list_id);

  return {};
}

const char* ExecutiveMapGenerate(PyMOLGlobals* G, const char* name,
                                 const char* reflection_file,
                                 const char* tempFile, char* amplitudes,
                                 const char* phases, const char* weights,
                                 double reso_low, double reso_high,
                                 const char* space_group, double cell[6],
                                 int quiet, int zoom) {
  /* FIXME: this should be returned in memory!! */
  /* In the meantime, use mkstemp and load in Python */
  int ok;

  ok = 0;

  if (weights && (!strncmp(weights, "None", 4))) weights = nullptr;

  /* printf("Passing to primex driver: space_group=%s, cell=[%f %f %f %f %f
   * %f], reso_high=%f, rseo_low=%f, refl_file=%s, ampl=%s, phases=%s,
   * weights=%s, map_file=%s", space_group, cell[0], cell[1], cell[2],
   * cell[3], cell[4], cell[5], reso_high, reso_low, reflection_file,
   * amplitudes, phases, weights, tempFile); */

  PRINTFB(G, FB_Executive, FB_Errors)
  " Error: MTZ map loading not supported in this PyMOL build.\n" ENDFB(G);

  if (!ok)
    return nullptr;
  else
    return (const char*)tempFile;
}

pymol::Result<> ExecutiveMapNew(PyMOLGlobals* G, const char* name, int type,
                                float grid_spacing, const char* s0,
                                float buffer, const float* minCorner,
                                const float* maxCorner, int state,
                                int have_corners, int quiet, int zoom,
                                int normalize, float clamp_floor,
                                float clamp_ceiling, float resolution) {
  pymol::CObject* origObj = nullptr;
  ObjectMap* objMap;
  ObjectMapState* ms = nullptr;
  int a;
  float v[3];
  ObjectMapDesc _md, *md;

  SETUP_SELE_DEFAULT(0);
  const char* sele = tmpsele0->getName();

  int isNew = true;
  int n_state;
  int valid_extent = false;
  int st;
  int st_once_flag = true;
  int n_st;
  int extent_state;
  int clamp_flag = (clamp_floor <= clamp_ceiling);

  md = &_md;

  float grid[3] = {grid_spacing, grid_spacing, grid_spacing};

  if ((state == -2) || (state == -3)) /* TO DO: support per-object states */
    state = SceneGetState(G);

  /* remove object if it already exists */

  origObj = ExecutiveFindObjectByName(G, name);

  if (origObj) {
    if (origObj->type != cObjectMap) {
      ExecutiveDelete(G, origObj->Name);
    } else {
      isNew = false;
    }
  }

  n_st = ExecutiveCountStates(G, nullptr);

  for (st = 0; st < n_st; st++) {
    if (state == -1)
      st_once_flag = false; /* each state, separate map, separate extent */
    if (!st_once_flag) state = st;
    extent_state = state;
    if (state <= -3) extent_state = -1;
    if (strlen(sele) && (!have_corners)) {
      valid_extent = ExecutiveGetExtent(G, sele, md->MinCorner, md->MaxCorner,
                                        true, extent_state, false);
      /* TODO restrict to state */
    } else {
      valid_extent = 1;
      copy3f(minCorner, md->MinCorner);
      copy3f(maxCorner, md->MaxCorner);
    }
    copy3f(grid, md->Grid);

    subtract3f(md->MaxCorner, md->MinCorner, v);
    for (a = 0; a < 3; a++) {
      if (v[a] < 0.0) std::swap(md->MaxCorner[a], md->MinCorner[a]);
    };
    subtract3f(md->MaxCorner, md->MinCorner, v);

    if (buffer < -0.5F) {  // buffer == -1
      buffer = SettingGetGlobal_f(G, cSetting_gaussian_resolution);
    }

    if (buffer > 0.0F) {
      for (a = 0; a < 3; a++) {
        md->MinCorner[a] -= buffer;
        md->MaxCorner[a] += buffer;
      }
    }
    md->mode = cObjectMap_OrthoMinMaxGrid;
    md->init_mode = -1; /* no initialization */

    /* validate grid */
    for (a = 0; a < 3; a++)
      if (md->Grid[a] <= R_SMALL8) md->Grid[a] = R_SMALL8;

    if (isNew)
      objMap = new ObjectMap(G);
    else
      objMap = (ObjectMap*)origObj;
    if (objMap) {
      int once_flag = true;
      n_state = SelectorCountStates(G, sele0);
      if (valid_extent)
        for (a = 0; a < n_state; a++) {
          if (state == -5)
            once_flag = false; /* api: state=-4 = each state, separate map,
                                  shared extent */
          if (state == -4)
            state = -1; /* api: state=-3 all states, but one map */
          if (!once_flag) state = a;
          ms = ObjectMapNewStateFromDesc(G, objMap, md, state, quiet);
          if (!ms) {
            return pymol::make_error("Invalid state ", state);
          }

          switch (type) {
            case 0: /* vdw */
              SelectorMapMaskVDW(G, sele0, ms, 0.0F, state);
              break;
            case 1: /* coulomb */
              SelectorMapCoulomb(G, sele0, ms, 0.0F, state, false, false, 1.0F);
              break;
            case 2: /* gaussian */
              SelectorMapGaussian(G, sele0, ms, 0.0F, state, normalize, false,
                                  quiet, resolution);
              break;
            case 3: /* coulomb_neutral */
              SelectorMapCoulomb(G, sele0, ms, 0.0F, state, true, false, 1.0F);
              break;
            case 4: /* coulomb_local */
              SelectorMapCoulomb(G, sele0, ms,
                                 SettingGetGlobal_f(G, cSetting_coulomb_cutoff),
                                 state, false, true, 2.0F);
              break;
            case 5: /* gaussian_max */
              SelectorMapGaussian(G, sele0, ms, 0.0F, state, normalize, true,
                                  quiet, resolution);
              break;
          }
          if (!ms->Active)
            ObjectMapStatePurge(G, ms);
          else if (clamp_flag) {
            ObjectMapStateClamp(ms, clamp_floor, clamp_ceiling);
          }

          if (once_flag) break;
        }

      ObjectSetName(objMap, name);
      ObjectMapUpdateExtents(objMap);
      if (isNew) {
        ExecutiveManageObject(G, objMap, -1, quiet);
      } else {
        ExecutiveDoZoom(G, objMap, false, zoom, true);
      }
      isNew = false;
      origObj = objMap;
    }
    SceneChanged(G);

    if (st_once_flag) break;
  }
  return {};
}
int ExecutiveMapSetBorder(PyMOLGlobals* G, const char* name, float level,
                          int state) {
  CExecutive* I = G->Executive;
  int result = true;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;

  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type == cObjectMap) {
            ObjectMap* obj = (ObjectMap*)rec->obj;
            result = ObjectMapSetBorder(obj, level, state);

            if (result) {
              ExecutiveInvalidateMapDependents(G, obj->Name);
            }
          }
          break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

pymol::Result<> ExecutiveMapDouble(PyMOLGlobals* G, const char* name,
                                   int state) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;

  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type == cObjectMap) {
            ObjectMap* obj = (ObjectMap*)rec->obj;
            auto result = ObjectMapDouble(obj, state);
            if (result) {
              ExecutiveInvalidateMapDependents(G, obj->Name);
            } else {
              return result;
            }
            if (result && rec->visible) SceneChanged(G);
          }
          break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return {};
}

pymol::Result<> ExecutiveMapHalve(PyMOLGlobals* G, const char* name, int state,
                                  int smooth) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;

  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type == cObjectMap) {
            ObjectMap* obj = (ObjectMap*)rec->obj;
            auto result = ObjectMapHalve(obj, state, smooth);
            if (result) {
              ExecutiveInvalidateMapDependents(G, obj->Name);
            } else {
              return result;
            }
            if (result && rec->visible) SceneChanged(G);
          }
          break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return {};
}

pymol::Result<> ExecutiveMapTrim(PyMOLGlobals* G, const char* name,
                                 const char* sele, float buffer, int map_state,
                                 int sele_state, int quiet) {
  auto s1 = SelectorTmp2::make(G, sele);
  CExecutive* I = G->Executive;
  float mn[3], mx[3];
  if (ExecutiveGetExtent(G, s1->getName(), mn, mx, true, sele_state, false)) {
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    SpecRec* rec;

    {
      int a;
      float t;
      for (a = 0; a < 3; a++) {
        mn[a] -= buffer;
        mx[a] += buffer;
        if (mn[a] > mx[a]) {
          t = mn[a];
          mn[a] = mx[a];
          mx[a] = t;
        }
      }
    }
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecObject:
            if (rec->obj->type == cObjectMap) {
              ObjectMap* obj = (ObjectMap*)rec->obj;
              auto result = ObjectMapTrim(obj, map_state, mn, mx, quiet);
              if (result) {
                ExecutiveInvalidateMapDependents(G, obj->Name);
              } else {
                return result;
              }
              if (result && rec->visible) SceneChanged(G);
            }
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return {};
}

pymol::Result<> ExecutiveSliceNew(PyMOLGlobals* G, const char* slice_name,
                                  const char* map_name, int state,
                                  int map_state) {
  int multi = false;
  pymol::CObject *obj = nullptr, *mObj, *origObj;
  ObjectMap* mapObj;
  ObjectMapState* ms;

  origObj = ExecutiveFindObjectByName(G, slice_name);
  if (origObj) {
    if (origObj->type != cObjectSlice) {
      return pymol::make_error("Object ", slice_name,
                               " is not an ObjectSlice.");
    }
  }

  mObj = ExecutiveFindObjectByName(G, map_name);
  if (mObj) {
    if (mObj->type != cObjectMap) mObj = nullptr;
  }
  if (mObj) {
    mapObj = (ObjectMap*)mObj;
    if (state == -1) {
      multi = true;
      state = 0;
      map_state = 0;
    } else if (state == -2) {
      state = SceneGetState(G);
      if (map_state < 0) map_state = state;
    } else if (state == -3) { /* append mode */
      state = 0;
      if (origObj) state = origObj->getNFrame();
    } else {
      if (map_state == -1) {
        map_state = 0;
        multi = true;
      } else {
        multi = false;
      }
    }
    while (1) {
      if (map_state == -2) map_state = SceneGetState(G);
      if (map_state == -3) map_state = mapObj->getNFrame() - 1;
      ms = ObjectMapStateGetActive(mapObj, map_state);
      if (ms) {
        obj = ObjectSliceFromMap(G, (ObjectSlice*)origObj, mapObj, state,
                                 map_state);

        if (!origObj) {
          ObjectSetName(obj, slice_name);
          ExecutiveManageObject(G, obj, -1, false);
        }
        PRINTFB(G, FB_ObjectMesh, FB_Actions)
        " SliceMap: created \"%s\".\n", slice_name ENDFB(G);

      } else if (!multi) {
        return pymol::make_error("State ", map_state + 1,
                                 " not present in map ", map_name);
      }
      if (multi) {
        origObj = obj;
        map_state++;
        state++;
        if (map_state >= mapObj->State.size()) break;
      } else {
        break;
      }
    }
  } else {
    return pymol::make_error("Map or brick object ", map_name, " not found.");
  }
  return {};
}
