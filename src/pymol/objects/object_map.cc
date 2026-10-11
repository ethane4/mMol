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

#include "pymol/base/err.h"
#include "pymol/objects/object_map.h"
#include "pymol/color/color.h"
#include "pymol/scene/scene.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/geometry/matrix.h"
#include "pymol/base/feedback.h"

/* MapState::ValidXtal -- determines whether the MapState's Xtal type passed in
 * is valid PARAMS ms, MapState RETURNS true/false
 */
int ObjectMapStateValidXtal(ObjectMapState* ms) {
  if (ms && ms->Active) {
    switch (ms->MapSource) {
        /* these are the only formats which are known to contain xtal
           information */
      case cMapSourceCrystallographic:
      case cMapSourceCCP4:
      case cMapSourceBRIX:
      case cMapSourceGRD:
        return true;
        break;
    }
  }
  return false;
}

/* Map::ValidXtal -- detemines whether the Map is valid
 * PARAMS
 *  I, Map
 *  state, map's state
 * RETURNS
 *  true/false
 */
int ObjectMapValidXtal(ObjectMap* I, int state) {
  if ((state >= 0) && (state < I->State.size())) {
    ObjectMapState* ms = &I->State[state];
    return ObjectMapStateValidXtal(ms);
  }
  return false;
}

void ObjectMapRegeneratePoints(ObjectMap* om) {
  int i;
  for (i = 0; i < om->State.size(); i++) {
    ObjectMapStateRegeneratePoints(&om->State[i]);
  }
}

void ObjectMapStateRegeneratePoints(ObjectMapState* ms) {
  int a, b, c, e;
  float v[3], vr[3];

  if (ObjectMapStateValidXtal(ms)) {
    for (c = 0; c < ms->FDim[2]; c++) {
      v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
      for (b = 0; b < ms->FDim[1]; b++) {
        v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
        for (a = 0; a < ms->FDim[0]; a++) {
          v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          for (e = 0; e < 3; e++) F4(ms->Field->points, a, b, c, e) = vr[e];
        }
      }
    }
  } else {
    for (c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for (b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for (a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for (e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }
  }
}

static PyObject* ObjectMapStateAsPyList(ObjectMapState* I) {
  PyObject* result = nullptr;

  result = PyList_New(16);
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  if (I->Symmetry) {
    PyList_SetItem(result, 1, SymmetryAsPyList(I->Symmetry.get()));
  } else {
    PyList_SetItem(result, 1, PConvAutoNone(Py_None));
  }
  if (!I->Origin.empty()) {
    PyList_SetItem(result, 2, PConvFloatArrayToPyList(I->Origin.data(), 3));
  } else {
    PyList_SetItem(result, 2, PConvAutoNone(Py_None));
  }
  if (!I->Range.empty()) {
    PyList_SetItem(result, 3, PConvFloatArrayToPyList(I->Range.data(), 3));
  } else {
    PyList_SetItem(result, 3, PConvAutoNone(Py_None));
  }
  if (!I->Dim.empty()) {
    PyList_SetItem(result, 4, PConvIntArrayToPyList(I->Dim.data(), 3));
  } else {
    PyList_SetItem(result, 4, PConvAutoNone(Py_None));
  }
  if (!I->Grid.empty()) {
    PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->Grid.data(), 3));
  } else {
    PyList_SetItem(result, 5, PConvAutoNone(Py_None));
  }
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->Corner, 24));
  PyList_SetItem(result, 7, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 8, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 9, PyInt_FromLong(I->MapSource));

  PyList_SetItem(result, 10, PConvIntArrayToPyList(I->Div, 3));
  PyList_SetItem(result, 11, PConvIntArrayToPyList(I->Min, 3));
  PyList_SetItem(result, 12, PConvIntArrayToPyList(I->Max, 3));
  PyList_SetItem(result, 13, PConvIntArrayToPyList(I->FDim, 4));

  PyList_SetItem(result, 14, IsosurfAsPyList(I->G, I->Field.get()));
  PyList_SetItem(result, 15, ObjectStateAsPyList(I));
  return (PConvAutoNone(result));
}

static PyObject* ObjectMapAllStatesAsPyList(ObjectMap* I) {
  PyObject* result = nullptr;
  int a;
  result = PyList_New(I->State.size());
  for (a = 0; a < I->State.size(); a++) {
    if (I->State[a].Active) {
      PyList_SetItem(result, a, ObjectMapStateAsPyList(&I->State[a]));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(nullptr));
    }
  }
  return (PConvAutoNone(result));
}

static int ObjectMapStateCopy(const ObjectMapState* src, ObjectMapState* I) {
  int ok = true;
  if (ok) {
    I->Active = src->Active;
    if (I->Active) {
      I->Symmetry = src->Symmetry;
      I->Origin = src->Origin;
      I->Range = src->Range;
      I->Grid = src->Grid;
      I->Dim = src->Dim;
      std::copy_n(src->Corner, 24, I->Corner);

      copy3f(src->ExtentMin, I->ExtentMin);
      copy3f(src->ExtentMax, I->ExtentMax);

      I->MapSource = src->MapSource;

      copy3f(src->Div, I->Div);
      copy3f(src->Min, I->Min);
      copy3f(src->Max, I->Max);
      copy3f(src->FDim, I->FDim);

      I->Field = src->Field;
      if (ok) ObjectMapStateRegeneratePoints(I);
    }
  }
  return (ok);
}

ObjectMapState::ObjectMapState(const ObjectMapState& src) : CObjectState(src) {
  ObjectMapStateCopy(&src, this);
}

ObjectMapState& ObjectMapState::operator=(const ObjectMapState& src) {
  CObjectState::operator=(src);
  ObjectMapStateCopy(&src, this);
  return *this;
}

static int ObjectMapStateFromPyList(PyMOLGlobals* G, ObjectMapState* I,
                                    PyObject* list) {
  int ok = true;
  int ll = 0;
  PyObject* tmp;
  if (ok) ok = (list != nullptr);
  if (ok) {
    if (!PyList_Check(list))
      I->Active = false;
    else {
      if (ok) ok = PyList_Check(list);
      if (ok) ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */

      if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
      if (ok) {
        tmp = PyList_GetItem(list, 1);
        if (tmp == Py_None)
          I->Symmetry = nullptr;
        else {
          I->Symmetry.reset(SymmetryNewFromPyList(G, tmp));
          ok = I->Symmetry != nullptr;
        }
        CPythonVal_Free(tmp);
      }
      if (ok) {
        tmp = PyList_GetItem(list, 2);
        if (CPythonVal_IsNone(tmp))
          I->Origin.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Origin);
        CPythonVal_Free(tmp);
      }
      if (ok) {
        tmp = PyList_GetItem(list, 3);
        if (CPythonVal_IsNone(tmp))
          I->Range.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Range);
        CPythonVal_Free(tmp);
      }
      if (ok) {
        tmp = PyList_GetItem(list, 4);
        if (CPythonVal_IsNone(tmp))
          I->Dim.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Dim);
        CPythonVal_Free(tmp);
      }
      if (ok) {
        tmp = PyList_GetItem(list, 5);
        if (CPythonVal_IsNone(tmp))
          I->Grid.clear();
        else
          ok = PConvFromPyObject(G, tmp, I->Grid);
        CPythonVal_Free(tmp);
      }
      if (ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->Corner,
                                            24);
      if (ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 7),
                                            I->ExtentMin, 3);
      if (ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 8),
                                            I->ExtentMax, 3);
      if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 9), &I->MapSource);
      if (ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 10), I->Div, 3);
      if (ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 11), I->Min, 3);
      if (ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 12), I->Max, 3);
      if (ok)
        ok = CPythonVal_PConvPyListToIntArrayInPlace_From_List(G, list, 13,
                                                               I->FDim, 4);
      if (ok) {
        tmp = CPythonVal_PyList_GetItem(G, list, 14);
        I->Field.reset(IsosurfNewFromPyList(G, tmp));
        ok = I->Field != nullptr;
        CPythonVal_Free(tmp);
      }
      if (ok && (ll > 15)) {
        tmp = CPythonVal_PyList_GetItem(G, list, 15);
        ok = ObjectStateFromPyList(G, tmp, I);
        CPythonVal_Free(tmp);
      }
      if (ok) ObjectMapStateRegeneratePoints(I);
    }
  }
  return (ok);
}

static int ObjectMapAllStatesFromPyList(ObjectMap* I, PyObject* list) {
  int ok = true;
  int a;
  if (ok) ok = PyList_Check(list);
  if (ok) {
    I->State.resize(PyList_Size(list), ObjectMapState(I->G));
    for (a = 0; a < I->State.size(); a++) {
      CPythonVal* val = CPythonVal_PyList_GetItem(I->G, list, a);
      ok = ObjectMapStateFromPyList(I->G, &I->State[a], val);
      CPythonVal_Free(val);
      if (!ok) break;
    }
  }
  return (ok);
}

PyObject* ObjectMapAsPyList(ObjectMap* I) {
  PyObject* result = nullptr;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->State.size()));
  PyList_SetItem(result, 2, ObjectMapAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}

int ObjectMapNewFromPyList(PyMOLGlobals* G, PyObject* list,
                           ObjectMap** result) {
  int ok = true;
  ObjectMap* I = nullptr;
  (*result) = nullptr;

  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  I = new ObjectMap(G);
  if (ok) ok = (I != nullptr);

  if (ok) {
    auto* val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if (ok) {
    CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 2);
    ok = ObjectMapAllStatesFromPyList(I, val);
    CPythonVal_Free(val);
  }
  if (ok) {
    (*result) = I;
    ObjectMapUpdateExtents(I);
  } else {
    /* cleanup? */
  }

  return (ok);
}

int ObjectMapNewCopy(PyMOLGlobals* G, const ObjectMap* src, ObjectMap** result,
                     int source_state, int target_state) {
  int ok = true;
  ObjectMap* I = nullptr;
  I = new ObjectMap(G);
  if (ok) ok = (I != nullptr);
  if (ok) ok = ObjectCopyHeader(I, src);
  if (ok) {
    if (source_state == -1) { /* all states */
      int state;
      VecCheckEmplace(I->State, I->State.size(), I->G);
      for (state = 0; state < src->State.size(); state++) {
        I->State[state] = src->State[state];
      }
    } else {
      if (target_state < 0) target_state = 0;
      if (source_state < 0) source_state = 0;
      VecCheckEmplace(I->State, target_state, G);
      if (source_state < src->State.size()) {
        I->State[target_state] = src->State[source_state];
      } else {
        ok = false;
        /* to do */
      }
    }
  }
  if (ok) *result = I;
  return ok;
}

ObjectMapState* ObjectMapStatePrime(ObjectMap* I, int state) {
  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  return &I->State[state];
}

void ObjectMapUpdateExtents(ObjectMap* I) {
  int a;
  float *min_ext, *max_ext;
  float tr_min[3], tr_max[3];
  I->ExtentFlag = false;

  for (a = 0; a < I->State.size(); a++) {
    ObjectMapState* ms = &I->State[a];
    if (ms->Active) {
      if (!ms->Matrix.empty()) {
        transform44d3f(ms->Matrix.data(), ms->ExtentMin, tr_min);
        transform44d3f(ms->Matrix.data(), ms->ExtentMax, tr_max);
        {
          float tmp;
          int a;
          for (a = 0; a < 3; a++)
            if (tr_min[a] > tr_max[a]) {
              tmp = tr_min[a];
              tr_min[a] = tr_max[a];
              tr_max[a] = tmp;
            }
        }
        min_ext = tr_min;
        max_ext = tr_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }

      if (!I->ExtentFlag) {
        copy3f(min_ext, I->ExtentMin);
        copy3f(max_ext, I->ExtentMax);
        I->ExtentFlag = true;
      } else {
        min3f(min_ext, I->ExtentMin, I->ExtentMin);
        max3f(max_ext, I->ExtentMax, I->ExtentMax);
      }
    }
  }

  if (I->TTTFlag && I->ExtentFlag) {
    const float* ttt;
    double tttd[16];
    if (ObjectGetTTT(I, &ttt, -1)) {
      convertTTTfR44d(ttt, tttd);
      MatrixTransformExtentsR44d3f(tttd, I->ExtentMin, I->ExtentMax,
                                   I->ExtentMin, I->ExtentMax);
    }
  }

  PRINTFD(I->G, FB_ObjectMap)
  " ObjectMapUpdateExtents-DEBUG: ExtentFlag %d\n", I->ExtentFlag ENDFD;
}

void ObjectMapStatePurge(PyMOLGlobals* G, ObjectMapState* I) {
  ObjectStatePurge(I);
  I->Field = nullptr;
  I->Origin.clear();
  I->Dim.clear();
  I->Range.clear();
  I->Grid.clear();
  I->shaderCGO = nullptr;

  I->Symmetry = nullptr;

  I->Active = false;
}

void ObjectMap::update() {
  auto I = this;
  if (!I->ExtentFlag) {
    ObjectMapUpdateExtents(I);
    if (I->ExtentFlag) SceneInvalidate(I->G);
  }
}

void ObjectMap::invalidate(cRep_t rep, cRepInv_t level, int state) {
  auto I = this;
  if (level >= cRepInvExtents) {
    I->ExtentFlag = false;
  }
  if ((rep < 0) || (rep == cRepDot)) {
    int a;
    for (a = 0; a < I->State.size(); a++) {
      if (I->State[a].Active) I->State[a].have_range = false;
      I->State[a].shaderCGO = nullptr;
    }
  }
  SceneInvalidate(I->G);
}

/* Has no prototype */
static CGO* ObjectMapCGOGenerate(PyMOLGlobals* G, float* corner) {
  CGO cgo_value(G);
  CGO* const shaderCGO = &cgo_value;
  CGOBegin(shaderCGO, GL_LINES);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 1);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 2);

  CGOVertexv(shaderCGO, corner + 3 * 2);
  CGOVertexv(shaderCGO, corner + 3 * 3);

  CGOVertexv(shaderCGO, corner + 3 * 1);
  CGOVertexv(shaderCGO, corner + 3 * 3);

  CGOVertexv(shaderCGO, corner + 3 * 0);
  CGOVertexv(shaderCGO, corner + 3 * 4);

  CGOVertexv(shaderCGO, corner + 3 * 1);
  CGOVertexv(shaderCGO, corner + 3 * 5);

  CGOVertexv(shaderCGO, corner + 3 * 2);
  CGOVertexv(shaderCGO, corner + 3 * 6);

  CGOVertexv(shaderCGO, corner + 3 * 3);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOVertexv(shaderCGO, corner + 3 * 4);
  CGOVertexv(shaderCGO, corner + 3 * 5);

  CGOVertexv(shaderCGO, corner + 3 * 4);
  CGOVertexv(shaderCGO, corner + 3 * 6);

  CGOVertexv(shaderCGO, corner + 3 * 6);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOVertexv(shaderCGO, corner + 3 * 5);
  CGOVertexv(shaderCGO, corner + 3 * 7);

  CGOEnd(shaderCGO);

  CGOStop(shaderCGO);

  return CGOOptimizeToVBONotIndexed(shaderCGO, 0, false);
}

void ObjectMap::render(RenderInfo* info) {
  auto I = this;
  int state = info->state;
  CRay* ray = info->ray;
  auto pick = info->pick;
  const RenderPass pass = info->pass;
  ObjectMapState* ms = nullptr;

  if (pass != RenderPass::Antialias) return;

  for (StateIterator iter(G, I->Setting.get(), state, I->State.size());
       iter.next();) {
    state = iter.state;
    if (I->State[state].Active) ms = &I->State[state];

    if (ms) {
      float* corner = ms->Corner;
      float tr_corner[24];
      ObjectPrepareContext(I, info);

      if (!ms->Matrix.empty()) { /* transform the corners before drawing */
        int a;
        for (a = 0; a < 8; a++) {
          transform44d3f(ms->Matrix.data(), corner + 3 * a, tr_corner + 3 * a);
        }
        corner = tr_corner;
      }

      if ((I->visRep & cRepExtentBit)) {
        if (ray) {
          const float* vc;
          float radius = ray->PixelRadius / 1.4142F;
          vc = ColorGet(G, I->Color);
          ray->color3fv(vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 1, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 2, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 2, corner + 3 * 3, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 1, corner + 3 * 3, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 0, corner + 3 * 4, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 1, corner + 3 * 5, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 2, corner + 3 * 6, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 3, corner + 3 * 7, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 4, corner + 3 * 5, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 4, corner + 3 * 6, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 6, corner + 3 * 7, radius, vc, vc);
          ray->sausage3fv(corner + 3 * 5, corner + 3 * 7, radius, vc, vc);
        } else if (G->HaveGUI && G->ValidContext) {
          if (pick) {
          } else {
            // shader
            if (!ms->shaderCGO) {
              ms->shaderCGO.reset(ObjectMapCGOGenerate(G, corner));
            }

            if (ms->shaderCGO) {
              CShaderPrg* shaderPrg =
                  G->ShaderMgr->Enable_DefaultShader(info->pass);
              if (shaderPrg) {
                shaderPrg->SetLightingEnabled(0);

                CGORender(ms->shaderCGO.get(), ColorGet(G, I->Color), nullptr,
                          nullptr, info, nullptr);
                shaderPrg->Disable();
              }
            }
          }
        }
      }

      if ((I->visRep & cRepDotBit)) {
        if (!ms->have_range) {
          double sum = 0.0, sumsq = 0.0;
          CField* data = ms->Field->data.get();
          int cnt = data->dim[0] * data->dim[1] * data->dim[2];
          float* raw_data = (float*)data->data.data();
          int a;
          for (a = 0; a < cnt; a++) {
            double f_val = *(raw_data++);
            sum += f_val;
            sumsq += (f_val * f_val);
          }
          if (cnt) {
            float mean, stdev;
            mean = (float)(sum / cnt);
            stdev = (float)sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
            ms->high_cutoff = mean + stdev;
            ms->low_cutoff = mean - stdev;
            ms->have_range = true;
          }
        }
        if (ms->have_range &&
            SettingGet_b(G, nullptr, I->Setting.get(), cSetting_dot_normals)) {
          IsofieldComputeGradients(G, ms->Field.get());
        }
        if (ms->have_range) {
          int a;
          CField* data = ms->Field->data.get();
          int cnt = data->dim[0] * data->dim[1] * data->dim[2];
          CField* points = ms->Field->points.get();
          CField* gradients = nullptr;

          if (SettingGet_b(G, nullptr, I->Setting.get(),
                           cSetting_dot_normals)) {
            gradients = ms->Field->gradients.get();
          }
          if (data && points) {
            float* raw_data = (float*)data->data.data();
            float raw_point[3], *raw_point_ptr = (float*)points->data.data();

#define RAW_POINT_TRANSFORM(ptr, v3f)              \
  {                                                \
    if (!ms->Matrix.empty())                       \
      transform44d3f(ms->Matrix.data(), ptr, v3f); \
    else                                           \
      copy3f(ptr, v3f);                            \
    ptr += 3;                                      \
  }

            float* raw_gradient = nullptr;
            float high_cut = ms->high_cutoff, low_cut = ms->low_cutoff;
            float width =
                SettingGet_f(G, nullptr, I->Setting.get(), cSetting_dot_width);

            if (ray) {
              float radius = ray->PixelRadius * width / 1.4142F;
              float vc[3];
              int color = I->Color;
              int ramped = ColorCheckRamped(G, I->Color);

              {
                const float* tmp = ColorGet(G, I->Color);
                copy3f(tmp, vc);
              }

              for (a = 0; a < cnt; a++) {
                float f_val = *(raw_data++);
                RAW_POINT_TRANSFORM(raw_point_ptr, raw_point);
                if ((f_val >= high_cut) || (f_val <= low_cut)) {
                  if (ramped) {
                    ColorGetRamped(G, color, raw_point, vc, state);
                    ray->color3fv(vc);
                  }
                  ray->sphere3fv(raw_point, radius);
                }
              }
            } else if (G->HaveGUI && G->ValidContext) {
              if (pick) {
              }
            }
          }
        }
      }
    }
  }
}

ObjectMapState::ObjectMapState(PyMOLGlobals* G) : CObjectState(G) {
  auto I = this;
  ObjectMapStatePurge(G, I);
  ObjectStateInit(G, I);
  I->Symmetry.reset(new CSymmetry(G));
  I->Field = nullptr;
  I->Origin.clear();
  I->Dim.clear();
  I->Range.clear();
  I->Grid.clear();
  I->MapSource = cMapSourceUndefined;
  I->have_range = false;
}

int ObjectMap::getNFrame() const { return State.size(); }
CObjectState* ObjectMap::_getObjectState(int state) {
  if (!State[state].Active) return nullptr;
  return &State[state];
}
CSymmetry const* ObjectMap::getSymmetry(int state) const {
  auto* ms = ObjectMapGetState(const_cast<ObjectMap*>(this), state);

  if (ms) {
    return ms->Symmetry.get();
  }

  return nullptr;
}
bool ObjectMap::setSymmetry(CSymmetry const& symmetry, int state) {
  bool success = false;

  for (StateIterator iter(G, Setting.get(), state, State.size());
       iter.next();) {
    auto& oms = State[iter.state];
    if (oms.Active) {
      oms.Symmetry.reset(new CSymmetry(symmetry));
      success = true;
    }
  }

  if (success) {
    ObjectMapRegeneratePoints(this);
  }

  return success;
}
ObjectMap::ObjectMap(PyMOLGlobals* G) : pymol::CObject(G) {
  auto I = this;
  I->type = cObjectMap;

  I->visRep = cRepExtentBit;
}
ObjectMapState* ObjectMapNewStateFromDesc(PyMOLGlobals* G, ObjectMap* I,
                                          ObjectMapDesc* inp_md, int state,
                                          int quiet) {
  int ok = true;
  float v[3];
  int a, b, c, d;
  float* fp;
  ObjectMapState* ms = nullptr;
  ObjectMapDesc _md, *md;
  ms = ObjectMapStatePrime(I, state);

  md = &_md;
  *(md) = *(inp_md);

  if (I) {
    ms->Origin = std::vector<float>(3);
    ms->Range = std::vector<float>(3);
    ms->Grid = std::vector<float>(3);
    ms->MapSource = cMapSourceDesc;
  }
  switch (md->mode) {
    case cObjectMap_OrthoMinMaxGrid: /* Orthorhombic: min, max, spacing,
                                        centered over range  */

      subtract3f(md->MaxCorner, md->MinCorner, v);
      for (a = 0; a < 3; a++) {
        if (v[a] < 0.0) std::swap(md->MaxCorner[a], md->MinCorner[a]);
      };
      subtract3f(md->MaxCorner, md->MinCorner, v);
      for (a = 0; a < 3; a++) {
        md->Dim[a] = (int)(v[a] / md->Grid[a]);
        if (md->Dim[a] < 1) md->Dim[a] = 1;
        if ((md->Dim[a] * md->Grid[a]) < v[a]) md->Dim[a]++;
      }

      PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMap: Dim %d %d %d\n", md->Dim[0], md->Dim[1],
          md->Dim[2] ENDFB(I->G);

      average3f(md->MaxCorner, md->MinCorner, v);
      for (a = 0; a < 3; a++) {
        md->MinCorner[a] = v[a] - 0.5F * md->Dim[a] * md->Grid[a];
      }

      if (Feedback(I->G, FB_ObjectMap, FB_Blather)) {
        dump3f(md->MinCorner, " ObjectMap: MinCorner:");
        dump3f(md->MaxCorner, " ObjectMap: MaxCorner:");
        dump3f(md->Grid, " ObjectMap: Grid:");
      }

      /* now populate the map data structure */

      copy3f(md->MinCorner, ms->Origin.data());
      copy3f(md->Grid, ms->Grid.data());
      for (a = 0; a < 3; a++) ms->Range[a] = md->Grid[a] * (md->Dim[a] - 1);

      /* these maps start at zero */
      for (a = 0; a < 3; a++) {
        ms->Min[a] = 0;
        ms->Max[a] = md->Dim[a] - 1;
        ms->Div[a] = md->Dim[a] - 1;
      }

      /* define corners */

      for (a = 0; a < 8; a++) copy3f(ms->Origin.data(), ms->Corner + 3 * a);

      d = 0;
      for (c = 0; c < 2; c++) {
        {
          v[2] = (c ? ms->Range[2] : 0.0F);
          for (b = 0; b < 2; b++) {
            v[1] = (b ? ms->Range[1] : 0.0F);
            for (a = 0; a < 2; a++) {
              v[0] = (a ? ms->Range[0] : 0.0F);
              add3f(v, ms->Corner + 3 * d, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }
      for (a = 0; a < 3; a++) ms->FDim[a] = ms->Max[a] - ms->Min[a] + 1;
      ms->FDim[3] = 3;

      ms->Field.reset(new Isofield(I->G, ms->FDim));
      if (!ms->Field)
        ok = false;
      else {
        for (a = 0; a < md->Dim[0]; a++) {
          v[0] = md->MinCorner[0] + a * md->Grid[0];
          for (b = 0; b < md->Dim[1]; b++) {
            v[1] = md->MinCorner[1] + b * md->Grid[1];
            for (c = 0; c < md->Dim[2]; c++) {
              v[2] = md->MinCorner[2] + c * md->Grid[2];
              fp = F4Ptr(ms->Field->points, a, b, c, 0);
              copy3f(v, fp);
            }
          }
        }
      }
      break;
    default:
      ok = false;
  }
  if (ok) {
    switch (md->init_mode) {
      case 0:
        for (a = 0; a < md->Dim[0]; a++) {
          for (b = 0; b < md->Dim[1]; b++) {
            for (c = 0; c < md->Dim[2]; c++) {
              F3(ms->Field->data, a, b, c) = 0.0F;
            }
          }
        }
        break;
      case 1:
        for (a = 0; a < md->Dim[0]; a++) {
          for (b = 0; b < md->Dim[1]; b++) {
            for (c = 0; c < md->Dim[2]; c++) {
              F3(ms->Field->data, a, b, c) = 1.0F;
            }
          }
        }
        break;
      case -2: /* for testing... */
        for (a = 0; a < md->Dim[0]; a++) {
          for (b = 0; b < md->Dim[1]; b++) {
            for (c = 0; c < md->Dim[2]; c++) {
              F3(ms->Field->data, a, b, c) =
                  (float)sqrt1d(a * a + b * b + c * c);
            }
          }
        }
        break;
    }
  }

  if (ok) {
    copy3f(ms->Origin.data(), ms->ExtentMin);
    copy3f(ms->Origin.data(), ms->ExtentMax);
    add3f(ms->Range.data(), ms->ExtentMax, ms->ExtentMax);
    ObjectMapUpdateExtents(I);
  }
  if (!ok) {
    ErrMessage(I->G, "ObjectMap", "Unable to create map");
    DeleteP(I);
  } else {
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Actions)
      " ObjectMap: Map created.\n" ENDFB(I->G);
    }
  }

  return (ms);
}

void ObjectMapDump(const ObjectMap* om, const char* fname, int state,
                   int quiet) {
  auto* oms = om->getObjectMapState(state);

  if (!oms) {
    ErrMessage(om->G, __func__, "state out of range");
    return;
  }

  auto file = fopen(fname, "wb");
  if (!file) {
    ErrMessage(om->G, __func__, "can't open file for writing");
    return;
  }

  auto* field = oms->Field.get();

  for (int xi = 0; xi < field->dimensions[0]; xi++) {
    for (int yi = 0; yi < field->dimensions[1]; yi++) {
      for (int zi = 0; zi < field->dimensions[2]; zi++) {
        float x = field->points->get<float>(xi, yi, zi, 0);
        float y = field->points->get<float>(xi, yi, zi, 1);
        float z = field->points->get<float>(xi, yi, zi, 2);

        switch (field->data->type) {
          case cFieldFloat: {
            float value = field->data->get<float>(xi, yi, zi);
            fprintf(file, "%10.4f%10.4f%10.4f%10.4f\n", x, y, z, value);
            break;
          }
          case cFieldInt: {
            int value = field->data->get<int>(xi, yi, zi);
            fprintf(file, "%10.4f%10.4f%10.4f%10d\n", x, y, z, value);
            break;
          }
          default:
            ErrMessage(om->G, __func__, "unknown field type");
            fclose(file);
            return;
        }
      }
    }
  }
  fclose(file);
  if (!quiet) {
    PRINTFB(om->G, FB_ObjectMap, FB_Actions)
    " ObjectMapDump: %s written to %s\n", om->Name, fname ENDFB(om->G);
  }
}

pymol::CObject* ObjectMap::clone() const { return new ObjectMap(*this); }
