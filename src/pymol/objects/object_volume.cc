
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

#include <algorithm>

#include "pymol/platform/os_python.h"

#include "pymol/platform/os_gl.h"
#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"

#include <math.h>

#include "pymol/base/base.h"
#include "pymol/color/color.h"
#include "pymol/base/err.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/fields/field.h"
#include "pymol/fields/isosurf.h"
#include "pymol/geometry/map.h"
#include "pymol/geometry/matrix.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/objects/object_volume.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/parse.h"
#include "pymol/scene/scene.h"
#include "pymol/settings/setting.h"
#include "pymol/graphics/shader_mgr.h"
#include "pymol/geometry/vector.h"
#include "pymol/app/main.h"

static void ObjectVolumeRecomputeExtent(ObjectVolume* I);

static ObjectVolumeState* ObjectVolumeGetActiveState(ObjectVolume* I) {
  if (!I) {
    return nullptr;
  }

  for (auto& state : I->State) {
    if (state.Active) {
      return &state;
    }
  }
  return nullptr;
}

CObjectState* ObjectVolume::_getObjectState(int state) { return &State[state]; }

static ObjectMapState* ObjectVolumeStateGetMapState(ObjectVolumeState* vs) {
  ObjectMap* map = nullptr;

  PyMOLGlobals* G = vs->G;

  map = ExecutiveFindObjectMapByName(G, vs->MapName);
  if (!map) {
    PRINTFB(G, FB_ObjectVolume, FB_Errors)
    "ObjectVolume-Error: map '%s' has been deleted.\n", vs->MapName ENDFB(G);
    return nullptr;
  }

  return ObjectMapGetState(map, vs->MapState);
}

ObjectMapState* ObjectVolumeGetMapState(ObjectVolume* I) {
  ObjectVolumeState* ovs = ObjectVolumeGetActiveState(I);
  if (ovs) return ObjectVolumeStateGetMapState(ovs);
  return nullptr;
}

static PyObject* ObjectVolumeStateAsPyList(ObjectVolumeState* I) {
  PyObject* result = nullptr;
  result = PyList_New(19);
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  PyList_SetItem(result, 1, PyString_FromString(I->MapName));
  PyList_SetItem(result, 2, PyInt_FromLong(I->MapState));
  PyList_SetItem(result, 3,
                 PConvAutoNone(nullptr) /* CrystalAsPyList(&I->Crystal) */);
  PyList_SetItem(result, 4, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(
      result, 7,
      PConvAutoNone(nullptr) /* PConvIntArrayToPyList(I->Range, 6) */);
  PyList_SetItem(result, 8, PyFloat_FromDouble(0.0 /* I->Level */));
  PyList_SetItem(result, 9, PyFloat_FromDouble(0.0 /* I->Radius */));
  PyList_SetItem(
      result, 10,
      PyInt_FromLong(/* I->CarveFlag */ I->AtomVertex.data() != nullptr));
  PyList_SetItem(result, 11, PyFloat_FromDouble(I->CarveBuffer));
  PyList_SetItem(result, 12,
                 I->AtomVertex ? PConvFloatVLAToPyList(I->AtomVertex)
                               : PConvAutoNone(nullptr));
  PyList_SetItem(result, 13, PyInt_FromLong(0 /* I->VolumeMode */));
  PyList_SetItem(result, 14, PyFloat_FromDouble(0.0 /* I->AltLevel */));
  PyList_SetItem(result, 15, PyInt_FromLong(1 /* I->quiet */));
  if (I->Field) {
    PyList_SetItem(result, 16, IsosurfAsPyList(I->G, I->Field.get()));
  } else {
    PyList_SetItem(result, 16, PConvAutoNone(nullptr));
  }
  PyList_SetItem(result, 17, PyInt_FromLong(I->RampSize()));
  if (!I->Ramp.empty()) {
    PyList_SetItem(result, 18, PConvToPyObject(I->Ramp));
  } else {
    PyList_SetItem(result, 18, PConvAutoNone(nullptr));
  }
  return (PConvAutoNone(result));
}

static PyObject* ObjectVolumeAllStatesAsPyList(ObjectVolume* I) {
  auto result = PyList_New(I->State.size());
  for (int a = 0; a < I->State.size(); a++) {
    if (I->State[a].Active) {
      PyList_SetItem(result, a, ObjectVolumeStateAsPyList(&I->State[a]));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(nullptr));
    }
  }
  return (PConvAutoNone(result));
}

static int ObjectVolumeStateFromPyList(PyMOLGlobals* G, ObjectVolumeState* I,
                                       PyObject* list) {
  int ok = true;
  int ll = 0;
  PyObject* tmp;
  if (ok) ok = (list != nullptr);
  if (ok) {
    if (!PyList_Check(list))
      I->Active = false;
    else {
      *I = ObjectVolumeState(G);
      if (ok) ok = (list != nullptr);
      if (ok) ok = PyList_Check(list);
      if (ok) ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
      if (ok)
        ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->MapName, WordLength);
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->MapState);
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->ExtentFlag);
      if (ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 5),
                                            I->ExtentMin, 3);
      if (ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6),
                                            I->ExtentMax, 3);
      if (ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 11), &I->CarveBuffer);
      if (ok) {
        tmp = PyList_GetItem(list, 12);
        if (tmp == Py_None)
          I->AtomVertex = nullptr;
        else
          ok = PConvPyListToFloatVLA(tmp, &I->AtomVertex);
      }
      if (ok) {
        I->RefreshFlag = true;
        I->ResurfaceFlag = true;
      }
      if (ok && (ll > 16)) {
        tmp = PyList_GetItem(list, 16);
        if (CPythonVal_IsNone(tmp))
          I->Field.reset();
        else {
          I->Field.reset(IsosurfNewFromPyList(G, tmp));
          ok = I->Field != nullptr;
        }
        CPythonVal_Free(tmp);
      }
      if (ok && (ll > 18)) {
        tmp = PyList_GetItem(list, 18);
        if (!CPythonVal_IsNone(tmp)) ok = PConvFromPyObject(G, tmp, I->Ramp);
        CPythonVal_Free(tmp);
      }
    }
  }
  return (ok);
}

static int ObjectVolumeAllStatesFromPyList(ObjectVolume* I, PyObject* list) {
  int ok = true;
  VecCheckEmplace(I->State, I->State.size(), I->G);
  if (ok) ok = PyList_Check(list);
  if (ok) {
    for (int a = 0; a < I->State.size(); a++) {
      CPythonVal* val = CPythonVal_PyList_GetItem(I->G, list, a);
      ok = ObjectVolumeStateFromPyList(I->G, &I->State[a], val);
      CPythonVal_Free(val);
      if (!ok) break;
    }
  }
  return (ok);
}

int ObjectVolumeNewFromPyList(PyMOLGlobals* G, PyObject* list,
                              ObjectVolume** result) {
  int ok = true;
  ObjectVolume* I = nullptr;
  (*result) = nullptr;
  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  I = new ObjectVolume(G);
  if (ok) ok = (I != nullptr);
  if (ok) {
    auto* val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if (ok) {
    CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 2);
    ok = ObjectVolumeAllStatesFromPyList(I, val);
    CPythonVal_Free(val);
  }
  if (ok) {
    (*result) = I;
    ObjectVolumeRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

PyObject* ObjectVolumeAsPyList(ObjectVolume* I) {
  auto result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->State.size()));
  PyList_SetItem(result, 2, ObjectVolumeAllStatesAsPyList(I));
  return (PConvAutoNone(result));
}

int ObjectVolumeInvalidateMapName(ObjectVolume* I, const char* name,
                                  const char* new_name) {
  int result = false;
  for (int a = 0; a < I->State.size(); a++) {
    auto vs = &I->State[a];
    if (vs->Active) {
      if (strcmp(vs->MapName, name) == 0) {
        if (new_name) strcpy(vs->MapName, new_name);
        I->invalidate(cRepAll, cRepInvAll, a);
        result = true;
      }
    }
  }
  return result;
}

void ObjectVolume::invalidate(cRep_t rep, cRepInv_t level, int state) {
  auto I = this;
  int a;
  int once_flag = true;
  if (level >= cRepInvExtents) {
    I->ExtentFlag = false;
  }

  PRINTFB(I->G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeInvalidate-Msg: %zu states.\n", I->State.size() ENDFB(I->G);

  if ((rep == cRepVolume) || (rep == cRepAll) || (rep == cRepExtent)) {
    for (a = 0; a < I->State.size(); a++) {
      if (state < 0) once_flag = false;
      if (!once_flag) state = a;
      if (level == cRepInvColor || level == cRepInvAll) {
        I->State[state].RecolorFlag = true;
      }
      if (level != cRepInvColor) {
        I->State[state].ResurfaceFlag = true;
        I->State[state].RefreshFlag = true;
      }
      SceneChanged(I->G);
      if (once_flag) break;
    }
  }
}

/**
 * Get the field either from the associated map, or from vs->Field in case
 * this is a reduced or symmetry expanded volume.
 */
static CField* ObjectVolumeStateGetField(ObjectVolumeState* vs) {
  if (!vs) return nullptr;
  if (vs->Field) return vs->Field->data.get();
  return ObjectVolumeStateGetMapState(vs)->Field->data.get();
}

CField* ObjectVolumeGetField(ObjectVolume* I) {
  return ObjectVolumeStateGetField(ObjectVolumeGetActiveState(I));
}

/**
 * Get a 4x4 (incl. translation) FracToReal from corner array
 */
static void get44FracToRealFromCorner(const float* corner, float* frac2real) {
  float tmp[16];
  identity44f(tmp);
  subtract3f(corner + 3, corner, tmp);
  subtract3f(corner + 6, corner, tmp + 4);
  subtract3f(corner + 12, corner, tmp + 8);
  copy3f(corner, tmp + 12);
  transpose44f44f(tmp, frac2real);
}

void ObjectVolume::update() {
  auto I = this;
  ObjectMapState* oms = nullptr;
  float carve_buffer;
  int avoid_flag = false;
  int flag;
  int h, k, l;
  int i, j;
  float range;
  MapType* voxelmap; /* this has nothing to do with isosurfaces... */

  for (int a = 0; a < I->State.size(); a++) {
    auto vs = &I->State[a];
    if (!vs || !vs->Active) continue;

    PRINTFD(G, FB_ObjectVolume)
    "ObjectVolumeUpdate: state=%d, refresh=%d, resurface=%d.\n", a,
        vs->RefreshFlag, vs->ResurfaceFlag ENDFD;

    oms = ObjectVolumeStateGetMapState(vs);

    if (!oms) {
      vs->ResurfaceFlag = false;
      continue;
    }

    if (vs->RefreshFlag || vs->ResurfaceFlag) {
      if (!oms->Matrix.empty()) {
        ObjectStateSetMatrix(vs, oms->Matrix.data());
      } else if (!vs->Matrix.empty()) {
        ObjectStateResetMatrix(vs);
      }

      // data min/max/mean/stdev
      range = SettingGet_f(I->G, I->Setting.get(), nullptr,
                           cSetting_volume_data_range);
      ObjectMapStateGetHistogram(I->G, oms, 0, range, vs->min_max_mean_stdev,
                                 0.f, 0.f);
    }

    // handle legacy or default color ramp
    if (vs->Ramp.empty() || (vs->RampSize() && vs->Ramp[0] == 0.f &&
                             vs->Ramp[5 * (vs->RampSize() - 1)] == 359.f)) {
      if (!vs->Ramp.empty()) {
        // legacy color ramp (0..359)
        range = vs->min_max_mean_stdev[1] - vs->min_max_mean_stdev[0];
        PRINTFB(G, FB_ObjectVolume, FB_Warnings)
        " ObjectVolumeUpdate: detected legacy color ramp\n" ENDFB(G);
        for (i = 0; i < vs->Ramp.size(); i += 5) {
          vs->Ramp[i] = vs->Ramp[i] / 359.f * range + vs->min_max_mean_stdev[0];
        }
      } else {
        // default color ramp (1.0 sigma peak)
        if (vs->Ramp.empty()) {
          vs->Ramp = std::vector<float>(
              {vs->min_max_mean_stdev[2] + 0.7f * vs->min_max_mean_stdev[3],
               0.f, 0.f, 1.f, 0.0f,
               vs->min_max_mean_stdev[2] + 1.0f * vs->min_max_mean_stdev[3],
               0.f, 1.f, 1.f, 0.2f,
               vs->min_max_mean_stdev[2] + 1.3f * vs->min_max_mean_stdev[3],
               0.f, 0.f, 1.f, 0.0f});
          vs->RecolorFlag = true;
        }
      }
    }

    if ((I->visRep & cRepVolumeBit) && vs->ResurfaceFlag) {
      Isofield* field = nullptr;
      vs->ResurfaceFlag = false;
      if (vs->Field) {
        field = vs->Field.get();
      } else if (oms->Field) {
        field = oms->Field.get();
      } else {
        field = nullptr;
      }

      if (field) {
        // get bounds and dimension data from field
        copy3(field->data->dim.data(), vs->dim);
        IsofieldGetCorners(G, field, vs->Corner);

        // transform corners by state matrix
        if (!vs->Matrix.empty()) {
          for (i = 0; i < 8; i++)
            transform44d3f(vs->Matrix.data(), vs->Corner + 3 * i,
                           vs->Corner + 3 * i);
        }
      }

      if (/* CarveFlag */ vs->AtomVertex) {
        carve_buffer = vs->CarveBuffer;
        if (vs->CarveBuffer < 0.0F) {
          avoid_flag = true;
          carve_buffer = -carve_buffer;
        }

        // cull my friend, cull */
        voxelmap = new MapType(I->G, -carve_buffer, vs->AtomVertex,
                               vs->AtomVertex.size() / 3, nullptr);
        if (voxelmap) {
          int x, y, z;
          int dx, dy, dz;
          float vv[3];
          float frac2real[16];

          MapSetupExpress(voxelmap);

          dx = vs->dim[0];
          dy = vs->dim[1];
          dz = vs->dim[2];

          get44FracToRealFromCorner(vs->Corner, frac2real);

          // initialize carve mask
          vs->carvemask = pymol::make_copyable<CField>(
              G, (int*)vs->dim, 3, sizeof(GLubyte), cFieldOther);

          // loop over voxels
          for (z = 0; z < dz; z++) {
            for (y = 0; y < dy; y++) {
              for (x = 0; x < dx; x++) {
                float frac[3] = {(x + .5f) / dx, (y + .5f) / dy,
                                 (z + .5f) / dz};

                transform44f3f(frac2real, frac, vv);
                flag = avoid_flag;

                // loop over close atoms
                MapLocus(voxelmap, vv, &h, &k, &l);
                for (i = *(MapEStart(voxelmap, h, k, l));
                     i && (j = voxelmap->EList[i]) >= 0; i++) {
                  if (within3f(vs->AtomVertex.data() + 3 * j, vv,
                               carve_buffer)) {
                    flag = !flag;
                    break;
                  }
                }

                // 0xFF (masked) or 0 (not masked), will be 1.0 or 0.0 in shader
                vs->carvemask->get<GLubyte>(x, y, z) = flag ? 0x0 : 0xFF;
              }
            }
          }
          MapFree(voxelmap);
        }
      }
    }
    vs->isUpdated = true;
    SceneInvalidate(I->G);
  }
  if (!I->ExtentFlag) {
    ObjectVolumeRecomputeExtent(I);
    if (I->ExtentFlag) SceneInvalidate(I->G);
  }
}

// Volumes were drawn only by the GL volume shader, which no longer exists.
void ObjectVolume::render(RenderInfo* info) {}

int ObjectVolume::getNFrame() const { return State.size(); }
ObjectVolume::ObjectVolume(PyMOLGlobals* G) : pymol::CObject(G) {
  type = cObjectVolume;
}
ObjectVolumeState::ObjectVolumeState(PyMOLGlobals* G)
    : CObjectState(G), Crystal(G), Active(true) {}
ObjectVolume* ObjectVolumeFromXtalSym(PyMOLGlobals* G, ObjectVolume* obj,
                                      ObjectMap* map, CSymmetry* sym,
                                      int map_state, int state, float* mn,
                                      float* mx, float level, int box_mode,
                                      float carve, float* vert_vla, int quiet) {
  ObjectVolume* I;
  ObjectVolumeState* vs;
  ObjectMapState* oms;
  int created = !obj;

  if (created) {
    I = new ObjectVolume(G);
  } else {
    I = obj;
  }
  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, G);
  }
  vs = &I->State[state];

  strcpy(vs->MapName, map->Name);
  vs->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);

  if (oms) {
    copy3f(mn, vs->ExtentMin); /* this is not exactly correct...should actually
                                  take vertex points from range */
    copy3f(mx, vs->ExtentMax);

    if (!oms->Matrix.empty()) {
      ObjectStateSetMatrix(vs, oms->Matrix.data());
    } else if (!vs->Matrix.empty()) {
      ObjectStateResetMatrix(vs);
    }

    {
      float *min_ext, *max_ext;
      float tmp_min[3], tmp_max[3];
      if (MatrixInvTransformExtentsR44d3f(vs->Matrix.data(), vs->ExtentMin,
                                          vs->ExtentMax, tmp_min, tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = vs->ExtentMin;
        max_ext = vs->ExtentMax;
      }

      if (sym && box_mode) {
        int eff_range[6];

        IsosurfGetRange(G, oms->Field.get(), &oms->Symmetry->Crystal, min_ext,
                        max_ext, eff_range, false);

        {
          int fdim[3];
          int expand_result;
          /* need to generate symmetry-expanded temporary map */

          fdim[0] = eff_range[3] - eff_range[0];
          fdim[1] = eff_range[4] - eff_range[1];
          fdim[2] = eff_range[5] - eff_range[2];
          vs->Field = pymol::make_copyable<Isofield>(I->G, fdim);

          expand_result =
              IsosurfExpand(oms->Field.get(), vs->Field.get(),
                            &oms->Symmetry->Crystal, sym, eff_range);

          if (expand_result == 0) {
            if (!quiet) {
              PRINTFB(G, FB_ObjectVolume, FB_Warnings)
              " ObjectVolume-Warning: no symmetry expanded map points "
              "found.\n" ENDFB(G);
            }
          } else {
            if (!quiet) {
              PRINTFB(G, FB_ObjectVolume, FB_Warnings)
              " ObjectVolume-Warning: not all symmetry expanded points covered "
              "by map.\n" ENDFB(G);
            }
          }
        }
      }
    }
    vs->ExtentFlag = true;
  }

  vs->CarveBuffer = carve;
  vs->AtomVertex = pymol::vla_take_ownership(vert_vla);

  I->ExtentFlag = false;

  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}
ObjectVolume* ObjectVolumeFromBox(PyMOLGlobals* G, ObjectVolume* obj,
                                  ObjectMap* map, int map_state, int state,
                                  float* mn, float* mx, float level,
                                  int meshMode, float carve, float* vert_vla,
                                  int quiet) {
  return ObjectVolumeFromXtalSym(G, obj, map, nullptr, map_state, state, mn, mx,
                                 level, meshMode, carve, vert_vla, quiet);
}
void ObjectVolumeRecomputeExtent(ObjectVolume* I) {
  int extent_flag = false;

  for (int a = 0; a < I->State.size(); a++) {
    auto vs = &I->State[a];
    if (vs->Active) {
      if (vs->ExtentFlag) {
        if (!extent_flag) {
          extent_flag = true;
          copy3f(vs->ExtentMax, I->ExtentMax);
          copy3f(vs->ExtentMin, I->ExtentMin);
        } else {
          max3f(vs->ExtentMax, I->ExtentMax, I->ExtentMax);
          min3f(vs->ExtentMin, I->ExtentMin, I->ExtentMin);
        }
      }
    }
  }

  I->ExtentFlag = extent_flag;

  if (I->TTTFlag && I->ExtentFlag) {
    const float* ttt;
    double tttd[16];
    if (ObjectGetTTT(I, &ttt, -1)) {
      convertTTTfR44d(ttt, tttd);
      MatrixTransformExtentsR44d3f(tttd, I->ExtentMin, I->ExtentMax,
                                   I->ExtentMin, I->ExtentMax);
    }
  }
}
PyObject* ObjectVolumeGetRamp(ObjectVolume* I, int state) {
  /* TODO: Allow for multi-state maps? */
  if (!I) {
    return (PConvAutoNone(nullptr));
  }
  if (auto ovs = static_cast<ObjectVolumeState*>(I->getObjectState(state))) {
    if (!ovs->isUpdated) I->update();

    return PConvFloatArrayToPyList(ovs->Ramp.data(), ovs->Ramp.size());
  }
  return (PConvAutoNone(nullptr));
}
pymol::Result<> ObjectVolumeSetRamp(ObjectVolume* I,
                                    std::vector<float>&& ramp_list, int state) {
  /* TODO: Allow for multi-state maps? */
  if (!I) {
    return pymol::make_error("Invalid volume.");
  }
  if (ramp_list.empty()) {
    return pymol::make_error("Empty ramp.");
  }
  auto ovs = static_cast<ObjectVolumeState*>(I->getObjectState(state));
  if (!ovs) {
    return pymol::make_error("Invalid volume state.");
  }

  ovs->Ramp = std::move(ramp_list);
  ovs->RecolorFlag = true;

  SceneChanged(I->G);

  return {};
}

pymol::CObject* ObjectVolume::clone() const { return new ObjectVolume(*this); }
