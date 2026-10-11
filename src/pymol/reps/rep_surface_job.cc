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

// The surface job: computing a molecular surface from atoms, off the main
// object.

#include "pymol/base/busy.h"
#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/reps/rep_surface.h"
#include "pymol/scene/scene.h"
#include "pymol/geometry/triangle.h"
#include "pymol/reps/rep_surface_internal.h"

static void SurfaceJobPurgeResult(PyMOLGlobals* G, SurfaceJob* I) {
  I->N = 0;
  I->NT = 0;
  I->V.clear();
  I->VN.clear();
  I->T.clear();
  I->S.clear();
}

static PyObject* SurfaceJobAtomInfoAsPyTuple(
    std::vector<SurfaceJobAtomInfo>& atom_info) {
  PyObject* result = nullptr;
  if (!atom_info.empty()) {
    auto size = 2 * atom_info.size() + 1;
    result = PyTuple_New(size);
    if (result) {
      PyTuple_SetItem(result, 0, PyInt_FromLong(2)); /* width of array */
      auto ai = atom_info.data();
      for (std::size_t i = 1; i < size; i += 2) {
        PyTuple_SetItem(result, i, PyFloat_FromDouble(ai->vdw));
        PyTuple_SetItem(result, i + 1, PyInt_FromLong(ai->flags));
        ai++;
      }
    }
  }
  return (PConvAutoNone(result));
}

PyObject* SurfaceJobInputAsTuple(PyMOLGlobals* G, SurfaceJob* I) {
  PyObject* result = PyTuple_New(24);
  if (result) {
    PyTuple_SetItem(result, 0, PyString_FromString("SurfaceJob"));
    PyTuple_SetItem(result, 1, PyInt_FromLong(1)); /* version */
    PyTuple_SetItem(result, 2, PConvToPyObject(I->coord));
    PyTuple_SetItem(result, 3, SurfaceJobAtomInfoAsPyTuple(I->atomInfo));
    PyTuple_SetItem(result, 4, PyFloat_FromDouble(I->maxVdw));
    PyTuple_SetItem(result, 5, PyInt_FromLong(I->allVisibleFlag));
    PyTuple_SetItem(result, 6, PyInt_FromLong(I->nPresent));
    PyTuple_SetItem(result, 7, PConvToPyObject(I->presentVla));
    PyTuple_SetItem(result, 8, PyInt_FromLong(I->solventSphereIndex));
    PyTuple_SetItem(result, 9, PyInt_FromLong(I->sphereIndex));
    PyTuple_SetItem(result, 10,
                    PyInt_FromLong(static_cast<int>(I->surfaceType)));
    PyTuple_SetItem(result, 11, PyInt_FromLong(I->circumscribe));
    PyTuple_SetItem(result, 12, PyFloat_FromDouble(I->probeRadius));
    PyTuple_SetItem(result, 13, PyFloat_FromDouble(I->carveCutoff));
    PyTuple_SetItem(result, 14, PConvToPyObject(I->carveVla));
    PyTuple_SetItem(result, 15, PyInt_FromLong(I->surfaceMode));
    PyTuple_SetItem(result, 16, PyInt_FromLong(I->surfaceSolvent));
    PyTuple_SetItem(result, 17, PyInt_FromLong(I->cavityCull));
    PyTuple_SetItem(result, 18, PyFloat_FromDouble(I->pointSep));
    PyTuple_SetItem(result, 19, PyFloat_FromDouble(I->trimCutoff));
    PyTuple_SetItem(result, 20, PyFloat_FromDouble(I->trimFactor));

    PyTuple_SetItem(result, 21, PyInt_FromLong(I->cavityMode));
    PyTuple_SetItem(result, 22, PyFloat_FromDouble(I->cavityRadius));
    PyTuple_SetItem(result, 23, PyFloat_FromDouble(I->cavityCutoff));
  }
  return result;
}

PyObject* SurfaceJobResultAsTuple(PyMOLGlobals* G, SurfaceJob* I) {
  PyObject* result = PyTuple_New(6);
  if (result) {
    PyTuple_SetItem(result, 0, PyInt_FromLong(I->N));
    PyTuple_SetItem(result, 1, PConvToPyObject(I->V));  // WAS TO TUPLE!
    PyTuple_SetItem(result, 2, PConvToPyObject(I->VN));
    PyTuple_SetItem(result, 3, PyInt_FromLong(I->NT));
    PyTuple_SetItem(result, 4, PConvToPyObject(I->T));
    PyTuple_SetItem(result, 5, PConvToPyObject(I->S));
  }
  return result;
}

ov_status SurfaceJobResultFromTuple(PyMOLGlobals* G, SurfaceJob* I,
                                    PyObject* tuple) {
  ov_status status = OV_STATUS_FAILURE;
  SurfaceJobPurgeResult(G, I);
  if (tuple && PyTuple_Check(tuple)) {
    ov_size size = PyTuple_Size(tuple);
    if (size >= 6) {
      status = OV_STATUS_SUCCESS;

      I->N = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
      if (OV_OK(status))
        status = PConvFromPyObject(G, PyTuple_GetItem(tuple, 1), I->V);
      if (OV_OK(status))
        status = PConvFromPyObject(G, PyTuple_GetItem(tuple, 2), I->VN);
      I->NT = PyInt_AsLong(PyTuple_GetItem(tuple, 3));
      if (OV_OK(status))
        status = PConvFromPyObject(G, PyTuple_GetItem(tuple, 4), I->T);
      if (OV_OK(status))
        status = PConvFromPyObject(G, PyTuple_GetItem(tuple, 5), I->S);
    }
    if (OV_ERR(status)) SurfaceJobPurgeResult(G, I);
  }
  return status;
}

static int SurfaceJobEliminateCloseDotsType3orMore(PyMOLGlobals* G,
                                                   SurfaceJob* I,
                                                   int* repeat_flag,
                                                   int* dot_flag) {
  int ok = true;
  int jj;
  float dist;
  float nearest;
  float point_sep = I->pointSep;
  float min_sep2 = point_sep * point_sep;
  float diff[3];
  {
    int a;
    for (a = 0; a < I->N; a++) dot_flag[a] = 1;
  }
  {
    MapType* map =
        new MapType(G, point_sep + 0.05F, I->V.data(), I->N, nullptr);
    int a;
    float* v = I->V.data();
    float* vn = I->VN.data();
    float min_dot = 0.1F;
    CHECKOK(ok, map);
    if (ok) ok &= MapSetupExpress(map);
    for (a = 0; ok && a < I->N; a++) {
      if (dot_flag[a]) {
        int i = *(MapLocusEStart(map, v));
        if (i && !map->EList.empty()) {
          int j = map->EList[i++];
          jj = I->N;
          nearest = point_sep + 1.0F;
          while (j >= 0) {
            if (j > a) {
              if (dot_flag[j]) {
                if (dot_product3f(I->VN.data() + (3 * j), vn) > min_dot) {
                  if (within3fret(I->V.data() + (3 * j), v, point_sep, min_sep2,
                                  diff, &dist)) {
                    *repeat_flag = true;
                    if (dist < nearest) {
                      /* try to be as determinstic as possible
                         in terms of how we collapse points */
                      jj = j;
                      nearest = dist;
                    } else if ((j < jj) && (fabs(dist - nearest) < R_SMALL4)) {
                      jj = j;
                      nearest = dist;
                    }
                  }
                }
              }
            }
            j = map->EList[i++];
          }
          if (jj < I->N) {
            dot_flag[jj] = 0;
            add3f(vn, I->VN.data() + (3 * jj), vn);
            average3f(I->V.data() + (3 * jj), v, v);
            *repeat_flag = true;
          }
        }
      }
      v += 3;
      vn += 3;
      ok &= !G->Interrupt;
    }
    MapFree(map);
  }
  return ok;
}

static int SurfaceJobEliminateCloseDotsTypeLessThan3(PyMOLGlobals* G,
                                                     SurfaceJob* I,
                                                     int* repeat_flag,
                                                     int* dot_flag) {
  int ok = true;
  int a;
  float point_sep = I->pointSep;
  MapType* map = new MapType(G, -point_sep, I->V.data(), I->N, nullptr);
  float* v = I->V.data();
  float* vn = I->VN.data();

  CHECKOK(ok, map);
  if (ok) {
    for (a = 0; a < I->N; a++) dot_flag[a] = 1;
    ok &= MapSetupExpress(map);
  }
  for (a = 0; ok && a < I->N; a++) {
    if (dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      if (i && !map->EList.empty()) {
        int j = map->EList[i++];
        while (j >= 0) {
          if (j != a) {
            if (dot_flag[j]) {
              if (within3f(I->V.data() + (3 * j), v, point_sep)) {
                dot_flag[j] = 0;
                add3f(vn, I->VN.data() + (3 * j), vn);
                average3f(I->V.data() + (3 * j), v, v);
                *repeat_flag = true;
              }
            }
          }
          j = map->EList[i++];
        }
      }
    }
    v += 3;
    vn += 3;
    ok &= !G->Interrupt;
  }
  MapFree(map);
  return ok;
}

static void SurfaceJobEliminateFromVArrays(PyMOLGlobals* G, SurfaceJob* I,
                                           int* dot_flag, short normalize) {
  float* v0 = I->V.data();
  float* vn0 = I->VN.data();
  int* p = dot_flag;
  int c = I->N;
  int a;
  float* v = I->V.data();
  float* vn = I->VN.data();
  I->N = 0;
  for (a = 0; a < c; a++) {
    if (*(p++)) {
      *(v0++) = *(v++);
      *(v0++) = *(v++);
      *(v0++) = *(v++);
      if (normalize) normalize3f(vn);
      *(vn0++) = *(vn++);
      *(vn0++) = *(vn++);
      *(vn0++) = *(vn++);
      I->N++;
    } else {
      v += 3;
      vn += 3;
    }
  }
}

static int SurfaceJobEliminateCloseDots(PyMOLGlobals* G, SurfaceJob* I) {
  int ok = true;
  if (I->N) {
    int repeat_flag = true;
    int* dot_flag = pymol::malloc<int>(I->N);
    CHECKOK(ok, dot_flag);
    while (ok && repeat_flag) {
      repeat_flag = false;
      if (static_cast<int>(I->surfaceType) >= 3) {
        ok = SurfaceJobEliminateCloseDotsType3orMore(G, I, &repeat_flag,
                                                     dot_flag);
      } else { /* surface types < 3 */
        ok = SurfaceJobEliminateCloseDotsTypeLessThan3(G, I, &repeat_flag,
                                                       dot_flag);
      }
      if (ok) {
        SurfaceJobEliminateFromVArrays(G, I, dot_flag, true);
      }
      ok &= !G->Interrupt;
    }
    FreeP(dot_flag);
  }
  return ok;
}

/* For each vertex, lookup all vertices within the neighborhood, and sum the
   dot_product of the normals. If the average of the dot_products of the normals
   is less than the trim_cutoff, then the middle vertex is eliminated. */
static int SurfaceJobEliminateTroublesomeVerticesMark(
    PyMOLGlobals* G, SurfaceJob* I, int* repeat_flag, MapType* map,
    int* dot_flag, float neighborhood, float trim_cutoff) {
  int ok = true;
  int a;
  float* v = I->V.data();
  float* vn = I->VN.data();
  for (a = 0; ok && a < I->N; a++) {
    if (dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      if (i && !map->EList.empty()) {
        int j = map->EList[i++];
        int n_nbr = 0;
        float dot_sum = 0.0F;
        while (j >= 0) {
          if (j != a) {
            if (dot_flag[j]) {
              float* v0 = I->V.data() + 3 * j;
              if (within3f(v0, v, neighborhood)) {
                float* n0 = I->VN.data() + 3 * j;
                dot_sum += dot_product3f(n0, vn);
                n_nbr++;
              }
            }
          }
          j = map->EList[i++];
        }
        if (n_nbr) {
          dot_sum /= n_nbr;
          if (dot_sum < trim_cutoff) {
            dot_flag[a] = false;
            *repeat_flag = true;
          }
        }
      }
    }
    v += 3;
    vn += 3;
    ok &= !G->Interrupt;
  }
  return ok;
}

static int SurfaceJobEliminateTroublesomeVertices(PyMOLGlobals* G,
                                                  SurfaceJob* I) {
  int ok = true;
  if ((I->surfaceType != SurfaceType::SolidDeterministic) && I->N &&
      (I->trimCutoff > 0.0F) && (I->trimFactor > 0.0F)) {
    float trim_cutoff = I->trimCutoff;
    float trim_factor = I->trimFactor;
    int repeat_flag = true;
    float point_sep = I->pointSep;
    float neighborhood = trim_factor * point_sep;
    int* dot_flag = pymol::malloc<int>(I->N);
    CHECKOK(ok, dot_flag);
    if (ok &&
        I->surfaceType == SurfaceType::SolidEmpirical) { /* emprical tweaks */
      trim_factor *= 2.5;
      trim_cutoff *= 1.5;
    }
    while (ok && repeat_flag) {
      MapType* map = new MapType(G, neighborhood, I->V.data(), I->N, nullptr);
      CHECKOK(ok, map);
      if (ok) {
        int a;
        for (a = 0; a < I->N; a++) dot_flag[a] = 1;
        ok &= MapSetupExpress(map);
      }
      repeat_flag = false;
      ok &= SurfaceJobEliminateTroublesomeVerticesMark(
          G, I, &repeat_flag, map, dot_flag, neighborhood, trim_cutoff);
      if (ok) {
        SurfaceJobEliminateFromVArrays(G, I, dot_flag, true);
      }
      MapFree(map);
      ok &= !G->Interrupt;
    }
    FreeP(dot_flag);
  }
  return ok;
}

static int SurfaceJobAtomProximityCleanupPass(PyMOLGlobals* G, SurfaceJob* I,
                                              int* dot_flag, int* present_vla,
                                              float probe_radius) {
  int ok = true;
  float cutoff = 0.5 * probe_radius;
  float* I_coord = I->coord.data();
  SurfaceJobAtomInfo* I_atom_info = I->atomInfo.data();
  int n_index = I->atomInfo.size();
  MapType* map = new MapType(G, I->maxVdw + probe_radius, I_coord, n_index,
                             nullptr, present_vla);
  int a;
  float* v;
  CHECKOK(ok, map);
  if (ok) ok &= MapSetupExpress(map);
  v = I->V.data();
  for (a = 0; ok && a < I->N; a++) {
    int i = *(MapLocusEStart(map, v));
    if (i && !map->EList.empty()) {
      int j = map->EList[i++];
      while (j >= 0) {
        SurfaceJobAtomInfo* atom_info = I_atom_info + j;
        if ((!present_vla) || present_vla[j]) {
          if (within3f(I_coord + 3 * j, v, atom_info->vdw + cutoff)) {
            dot_flag[a] = true;
          }
        }
        j = map->EList[i++];
      }
    }
    v += 3;
    ok &= !G->Interrupt;
  }
  MapFree(map);
  return ok;
}

static int SurfaceJobRefineCopyNewPoints(SurfaceJob* I, float* new_dot,
                                         int n_new) {
  int ok = true;
  float* n1 = new_dot + 3;
  float* v1 = new_dot;
  float *v, *vn;
  I->V.resize(3 * (I->N + n_new));
  CHECKOK(ok, !I->V.empty());
  if (ok) {
    I->VN.resize(3 * (I->N + n_new));
  }
  CHECKOK(ok, !I->VN.empty());
  if (ok) {
    v = I->V.data() + 3 * I->N;
    vn = I->VN.data() + 3 * I->N;
    I->N += n_new;
  }
  while (ok && n_new--) {
    copy3f(v1, v);
    copy3f(n1, vn);
    v += 3;
    vn += 3;
    v1 += 6;
    n1 += 6;
  }
  return ok;
}
static int SurfaceJobRefineAddNewVerticesCheckPoint(
    SurfaceJob* I, MapType* map, int* n_new, float** new_dot, int j, float* v,
    float* vn, float map_cutoff, float neighborhood, float insert_cutoff) {
  int ok = true;
  float* v0 = I->V.data() + 3 * j;
  if (within3f(v0, v, map_cutoff)) {
    int add_new = false;
    float* n0 = I->VN.data() + 3 * j;
    VLACheck(*new_dot, float, (*n_new) * 6 + 5);
    CHECKOK(ok, *new_dot);
    if (ok) {
      float* v1 = (*new_dot) + (*n_new) * 6;
      average3f(v, v0, v1);
      if ((dot_product3f(n0, vn) <
           0.666 /* dot_cutoff, was hardcoded as a variable */) &&
          (within3f(v0, v, neighborhood))) {
        // if the normals are further than dot_cutoff apart
        // and the related points are close to each other, than add new
        add_new = true;
      } else {
        /* if points are too far apart, insert a new one, i.e.,
           search for any point within insert_cutoff, if not, add */
        int ii = *(MapLocusEStart(map, v1));
        if (ii) {
          int found = false;
          int jj = map->EList[ii++];
          while (jj >= 0) {
            if (jj != j) {
              float* vv0 = I->V.data() + 3 * jj;
              if (within3f(vv0, v1, insert_cutoff)) {
                found = true;
                break;
              }
            }
            jj = map->EList[ii++];
          }
          if (!found) add_new = true;
        }
      }
      if (add_new) {
        /* highly divergent, add dot in-between v and v0
           (averaged, v1 set above, compute the normal (averaged below))
           and n_new incremented */
        float* n1 = v1 + 3;
        (*n_new)++;
        average3f(vn, n0, n1);
        normalize3f(n1);
      }
    }
  }
  return ok;
}

static int SurfaceJobRefineAddNewVertices(PyMOLGlobals* G, SurfaceJob* I) {
  int ok = true;
  float point_sep = I->pointSep;
  int n_new = 0;
  float neighborhood =
      2.6 * point_sep; /* these constants need more tuning... */
  float insert_cutoff = 1.1 * point_sep;
  float map_cutoff = neighborhood;
  float* new_dot = VLAlloc(float, 1000);
  float *v, *vn;
  if (map_cutoff <
      (2.9 * point_sep)) { /* these constants need more tuning... */
    map_cutoff = 2.9 * point_sep;
  }
  {
    MapType* map = nullptr;
    int a;
    map = new MapType(G, map_cutoff, I->V.data(), I->N, nullptr);
    CHECKOK(ok, map);
    if (ok) ok &= MapSetupExpress(map);
    v = I->V.data();
    vn = I->VN.data();
    for (a = 0; ok && a < I->N; a++) {
      int i = *(MapLocusEStart(map, v));
      if (i && !map->EList.empty()) {
        int j = map->EList[i++];
        while (ok && j >= 0) {
          if (j > a) {
            SurfaceJobRefineAddNewVerticesCheckPoint(
                I, map, &n_new, &new_dot, j, v, vn, map_cutoff, neighborhood,
                insert_cutoff);
          }
          j = map->EList[i++];
          ok &= !G->Interrupt;
        }
      }
      v += 3;
      vn += 3;
      ok &= !G->Interrupt;
    }
    MapFree(map);
  }
  if (ok && n_new) {
    ok = SurfaceJobRefineCopyNewPoints(I, new_dot, n_new);
  }
  VLAFreeP(new_dot);
  return ok;
}

static void SurfaceJobCheckInteriorSolventSurface(MapType* solv_map, float* v,
                                                  SolventDot* sol_dot,
                                                  float probe_rad_less,
                                                  float dist2, int a,
                                                  int* flag) {
  int ii;
  ii = *(MapLocusEStart(solv_map, v));
  if (ii && !solv_map->EList.empty()) {
    float* i_dot = sol_dot->dot;
    float dist = probe_rad_less;
    int* elist_ii = solv_map->EList.data() + ii;
    float v_0 = v[0];
    int jj_next, jj = *(elist_ii++);
    float v_1 = v[1];
    float* v1 = i_dot + 3 * jj;
    float v_2 = v[2];
    while (jj >= 0) {
      /* huge bottleneck -- optimized for superscaler processors */
      float dx = v1[0], dy, dz;
      jj_next = *(elist_ii++);
      dx -= v_0;
      if (jj != a) {
        dx = (dx < 0.0F) ? -dx : dx;
        dy = v1[1] - v_1;
        if (!(dx > dist)) {
          dy = (dy < 0.0F) ? -dy : dy;
          dz = v1[2] - v_2;
          if (!(dy > dist)) {
            dx = dx * dx;
            dz = (dz < 0.0F) ? -dz : dz;
            dy = dy * dy;
            if (!(dz > dist)) {
              dx = dx + dy;
              dz = dz * dz;
              if (!(dx > dist2))
                if ((dx + dz) <= dist2) {
                  *flag = false;
                  break;
                }
            }
          }
        }
      }
      v1 = i_dot + 3 * jj_next;
      jj = jj_next;
    }
  }
}

static void SurfaceJobCheckPresentAndWithin(MapType* map, SurfaceJob* I,
                                            int* present_vla, float* v,
                                            float probe_rad_more, int* flag) {
  SurfaceJobAtomInfo* I_atom_info = I->atomInfo.data();
  float* I_coord = I->coord.data();

  int i = *(MapLocusEStart(map, v));
  if (i && !map->EList.empty()) {
    int j = map->EList[i++];
    while (j >= 0) {
      SurfaceJobAtomInfo* atom_info = I_atom_info + j;
      if ((!present_vla) || present_vla[j]) {
        if (within3f(I_coord + 3 * j, v, atom_info->vdw + probe_rad_more)) {
          *flag = false;
          break;
        }
      }
      j = map->EList[i++];
    }
  }
}

static void SurfaceJobSetProbeRadius(SurfaceType surface_type, float point_sep,
                                     float* probe_radius, float* probe_rad_more,
                                     float* probe_rad_less,
                                     float* probe_rad_less2) {
  float solv_tole = point_sep * 0.04F;

  if (*probe_radius < (2.5F * point_sep)) { /* minimum probe radius allowed */
    *probe_radius = 2.5F * point_sep;
  }
  *probe_rad_more = *probe_radius * (1.0F + solv_tole);
  switch (surface_type) {
    case SurfaceType::Solid:
    case SurfaceType::SolidDeterministic:
    case SurfaceType::SolidWeightings:
    case SurfaceType::SolidScribed:
    case SurfaceType::SolidEmpirical:
      *probe_rad_less = *probe_radius;
      break;
    default:
      *probe_rad_less = *probe_radius * (1.0F - solv_tole);
      break;
  }
  *probe_rad_less2 = (*probe_rad_less) * (*probe_rad_less);
}

int SurfaceJobRun(PyMOLGlobals* G, SurfaceJob* I) {
  int ok = true;
  int MaxN;
  int n_present = I->nPresent;
  SphereRec* sp = G->Sphere->Sphere[I->sphereIndex];
  SphereRec* ssp = G->Sphere->Sphere[I->solventSphereIndex];

  SurfaceJobPurgeResult(G, I);

  {
    /* compute limiting storage requirements */
    int tmp = n_present;
    if (tmp < 1) tmp = 1;
    if (sp->nDot < ssp->nDot)
      MaxN = tmp * ssp->nDot;
    else
      MaxN = tmp * sp->nDot;
  }
  MaxN = n_present;
  //  MaxN = 0;
  I->V.resize((MaxN + 1) * 3);
  CHECKOK(ok, !I->V.empty());
  if (ok) {
    I->VN.resize((MaxN + 1) * 3);
  }
  CHECKOK(ok, !I->VN.empty());

  ok &= !G->Interrupt;

  if (!(!I->V.empty() && !I->VN.empty()) ||
      (!ok)) { /* bail out point -- try to reduce crashes */
    PRINTFB(G, FB_RepSurface, FB_Errors)
    "Error-RepSurface: insufficient memory to calculate surface at this "
    "quality.\n" ENDFB(G);
    I->V.clear();
    I->VN.clear();
    ok = false;
  } else {
    SolventDot* sol_dot = nullptr;
    float* v = I->V.data();
    float* vn = I->VN.data();
    MapType* carve_map = nullptr;
    float probe_radius = I->probeRadius;
    int circumscribe = I->circumscribe;
    SurfaceType surface_type = I->surfaceType;
    float point_sep = I->pointSep;
    int* present_vla = I->presentVla.data();

    I->N = 0;

    const bool surf_timing = getenv("PYMOL_SURFACE_TIMING") != nullptr;
    double t_soldot0 = surf_timing ? UtilGetSeconds(G) : 0.0;

    sol_dot = SolventDotNew(G, I->coord.data(), I->atomInfo, probe_radius, ssp,
                            present_vla, circumscribe, I->surfaceMode,
                            I->surfaceSolvent, I->cavityCull, I->allVisibleFlag,
                            I->maxVdw, I->cavityMode, I->cavityRadius,
                            I->cavityCutoff);
    CHECKOK(ok, sol_dot);
    BusyProgress(0.03f);
    if (surf_timing && sol_dot) {
      fprintf(stderr,
              "SURFACE_TIMING: stage=solvent_dot nAtom=%d nDot=%d ms=%.3f\n",
              (int)I->atomInfo.size(), sol_dot->nDot,
              (UtilGetSeconds(G) - t_soldot0) * 1000.0);
      fflush(stderr);
    }
    ok &= !G->Interrupt;
    if (ok) {
      if (!I->surfaceSolvent) {
        float probe_rad_more, probe_rad_less, probe_rad_less2;

        SurfaceJobSetProbeRadius(surface_type, point_sep, &probe_radius,
                                 &probe_rad_more, &probe_rad_less,
                                 &probe_rad_less2);

        if (surface_type == SurfaceType::SolidScribed ||
            surface_type == SurfaceType::SolidEmpirical) {
          /* effectively double-weights atom points */
          if (sol_dot->nDot) {
            if (sol_dot->nDot > MaxN) {
              MaxN = sol_dot->nDot;
              VecCheck(I->V, 3 * (MaxN + 1));
              VecCheck(I->VN, 3 * (MaxN + 1));
              v = I->V.data();
              vn = I->VN.data();
            }
            int a;
            float* v0 = sol_dot->dot;
            float* n0 = sol_dot->dotNormal;
            for (a = 0; a < sol_dot->nDot; a++) {
              scale3f(n0, -probe_radius, v);
              add3f(v0, v, v);
              copy3f(n0, vn);
              v += 3;
              vn += 3;
              n0 += 3;
              v0 += 3;
              I->N++;
            }
          }
        }
        ok &= !G->Interrupt;
        if (ok) {
          MapType* solv_map{};
          auto map = new MapType(G, I->maxVdw + probe_rad_more, I->coord.data(),
                                 I->coord.size() / 3, nullptr, nullptr);
          CHECKOK(ok, map);
          if (ok)
            solv_map = new MapType(G, probe_rad_less, sol_dot->dot,
                                   sol_dot->nDot, nullptr);
          CHECKOK(ok, solv_map);
          if (ok) {
            ok &= MapSetupExpress(solv_map);
            if (ok) ok &= MapSetupExpress(map);
            ok &= !G->Interrupt;
            ok &= !map->EList.empty() && !solv_map->EList.empty();
            if (sol_dot->nDot && ok) {
              Vector3f* dot = pymol::malloc<Vector3f>(sp->nDot);
              float *v0, *n0;
              CHECKOK(ok, dot);
              if (ok) {
                int b;
                for (b = 0; b < sp->nDot; b++) {
                  scale3f(sp->dot[b], probe_radius, dot[b]);
                }
              }
              v0 = sol_dot->dot;
              n0 = sol_dot->dotNormal;
              if (ok) {
                int a, b;
                int sp_nDot = sp->nDot;
                for (a = 0; ok && a < sol_dot->nDot; a++) {
                  if (sol_dot->dotCode[a] ||
                      (surface_type != SurfaceType::SolidEmpirical)) {
                    OrthoBusyFast(G, a + sol_dot->nDot * 2,
                                  sol_dot->nDot * 5); /* 2/5 to 3/5 */
                    BusyProgress(0.03f + 0.62f * a / sol_dot->nDot);
                    for (b = 0; ok && b < sp_nDot; b++) {
                      float* dot_b = dot[b];
                      v[0] = v0[0] + dot_b[0];
                      v[1] = v0[1] + dot_b[1];
                      v[2] = v0[2] + dot_b[2];
                      {
                        int flag = true;
                        SurfaceJobCheckInteriorSolventSurface(
                            solv_map, v, sol_dot, probe_rad_less,
                            probe_rad_less2, a, &flag);
                        /* at this point, we have points on the interior of the
                           solvent surface, so now we need to further trim that
                           surface to cover atoms that are present */
                        if (flag) {
                          SurfaceJobCheckPresentAndWithin(
                              map, I, present_vla, v, probe_rad_more, &flag);
                          if (!flag) { /* compute the normals */
                            vn[0] = -sp->dot[b][0];
                            vn[1] = -sp->dot[b][1];
                            vn[2] = -sp->dot[b][2];
                            I->N++;
                            VecCheck(I->V, 3 * (I->N + 1));
                            VecCheck(I->VN, 3 * (I->N + 1));
                            CHECKOK(ok, !I->V.empty());
                            CHECKOK(ok, !I->VN.empty());
                            v = I->V.data() + I->N * 3;
                            vn = I->VN.data() + I->N * 3;
                          }
                        }
                      }
                      ok &= !G->Interrupt;
                    }
                  }
                  v0 += 3;
                  n0 += 3;
                  ok &= !G->Interrupt;
                }
              }
              FreeP(dot);
            }
          }
          MapFree(solv_map);
          MapFree(map);
        }
      } else {
        float* v0 = sol_dot->dot;
        float* n0 = sol_dot->dotNormal;
        int a;
        circumscribe = 0;
        if (sol_dot->nDot) {
          VecCheck(I->V, 3 * (I->N + sol_dot->nDot));
          VecCheck(I->VN, 3 * (I->N + sol_dot->nDot));
          v = I->V.data();
          vn = I->VN.data();
          for (a = 0; a < sol_dot->nDot; a++) {
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            I->N++;
          }
        }
      }
    }
    SolventDotFree(sol_dot);
    sol_dot = nullptr;
    ok &= !G->Interrupt;
    if (ok) {
      BusyProgress(0.65f);
      int refine, ref_count = 1;

      if ((surface_type == SurfaceType::Solid) && (circumscribe)) {
        ref_count = 2; /* these constants need more tuning... */
      }
      for (refine = 0; ok && refine < ref_count; refine++) {
        /* add new vertices in regions where curvature is very high
           or where there are gaps with no points */
        if (I->N && (surface_type == SurfaceType::Solid) && (circumscribe)) {
          ok = SurfaceJobRefineAddNewVertices(G, I);
        }

        if (ok && I->N && (surface_type == SurfaceType::Solid) &&
            (circumscribe)) {
          /* combine scribing with an atom proximity cleanup pass */
          int* dot_flag = pymol::calloc<int>(I->N);
          CHECKOK(ok, dot_flag);
          ok &= SurfaceJobAtomProximityCleanupPass(G, I, dot_flag, present_vla,
                                                   probe_radius);
          /* purge unused dots */
          if (ok)
            SurfaceJobEliminateFromVArrays(G, I, dot_flag,
                                           false);  // not normalize?
          FreeP(dot_flag);
        }

        /* now, eliminate dots that are too close to each other */
        ok &= SurfaceJobEliminateCloseDots(G, I);

        /* now eliminate troublesome vertices in regions of extremely high
         * curvature */
        ok &= SurfaceJobEliminateTroublesomeVertices(G, I);
        ok &= !G->Interrupt;
      }
    }

    if (ok && I->N && !I->V.empty() && !I->VN.empty()) {
      // VLASizeForSure(I->V, float, 3 * I->N);
      I->V.resize(3 * I->N);
      CHECKOK(ok, !I->V.empty());
      if (ok) {
        // VLASizeForSure(I->VN, float, 3 * I->N);
        I->VN.resize(3 * I->N);
      }
      CHECKOK(ok, !I->VN.empty());
    }

    PRINTFB(G, FB_RepSurface, FB_Blather)
    " RepSurface: %i surface points.\n", I->N ENDFB(G);

    ok &= !G->Interrupt;

    OrthoBusyFast(G, 3, 5);
    BusyProgress(0.70f);
    if (I->N) {
      if (ok && surface_type != SurfaceType::DotDefault) {
        float cutoff = point_sep * 5.0F;
        if ((cutoff > probe_radius) && (!I->surfaceSolvent))
          cutoff = probe_radius;
        double t_tri0 = surf_timing ? UtilGetSeconds(G) : 0.0;
        I->T =
            TrianglePointsToSurface(G, I->V.data(), I->VN.data(), I->N, cutoff,
                                    &I->NT, I->S, nullptr, I->cavityMode);
        CHECKOK(ok, !I->T.empty());
        if (surf_timing) {
          fprintf(
              stderr,
              "SURFACE_TIMING: stage=triangulate nVert=%d nTri=%d ms=%.3f\n",
              I->N, I->NT, (UtilGetSeconds(G) - t_tri0) * 1000.0);
          fflush(stderr);
        }
        BusyProgress(0.97f);
        PRINTFB(G, FB_RepSurface, FB_Blather)
        " RepSurface: %i triangles.\n", I->NT ENDFB(G);
      }
    } else {
      if (ok) I->V.resize(1);
      CHECKOK(ok, !I->V.empty());
      if (ok) I->VN.resize(1);
      CHECKOK(ok, !I->VN.empty());
    }
    if (carve_map) MapFree(carve_map);
  }
  return ok;
}
