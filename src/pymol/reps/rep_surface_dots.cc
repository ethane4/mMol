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

// Placing solvent dots around atoms for the surface job.

#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/reps/rep_surface.h"
#include "pymol/scene/scene.h"
#include <omp.h>
#include "pymol/reps/rep_surface_internal.h"

static int SolventDotFilterOutSameXYZ(PyMOLGlobals* G, MapType* map,
                                      SurfaceJobAtomInfo* atom_info,
                                      SurfaceJobAtomInfo* a_atom_info,
                                      float* coord, int a, int* present,
                                      int* skip_flag) {
  int ok = true;
  float* v0 = coord + 3 * a;
  int i = *(MapLocusEStart(map, v0));
  if (i && !map->EList.empty()) {
    int j = map->EList[i++];
    while (ok && j >= 0) {
      SurfaceJobAtomInfo* j_atom_info = atom_info + j;
      if (j > a) /* only check if this is atom trails */
        if ((!present) || present[j]) {
          if (j_atom_info->vdw == a_atom_info->vdw) { /* handle singularities */
            float* v1 = coord + 3 * j;
            if ((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2])) {
              // is this necessary? if some different atoms have exact same
              // x/y/z
              *skip_flag = true;
            }
          }
        }
      j = map->EList[i++];
      ok &= !G->Interrupt;
    }
  }
  return ok;
}

static int SolventDotGetDotsAroundVertexInSphere(
    PyMOLGlobals* G, SolventDot* I, MapType* map, SurfaceJobAtomInfo* atom_info,
    SurfaceJobAtomInfo* a_atom_info, float* coord, int a, int* present,
    SphereRec* sp, float radius, int* dotCnt, int stopDot, float* dotPtr,
    float* dotNormal, int* nDot) {
  float vdw = a_atom_info->vdw + radius;
  float* v0 = coord + 3 * a;
  int b, ok = true;
  float* v = dotPtr + (*nDot) * 3;
  float* n = nullptr;
  Vector3f* sp_dot = sp->dot;
  if (dotNormal) {
    n = dotNormal + (*nDot) * 3;
  }
  for (b = 0; ok && b < sp->nDot; b++) {
    float* sp_dot_b = (float*)(sp_dot + b);
    int i;
    int flag = true;
    if (n) {
      n[0] = sp_dot_b[0];
      n[1] = sp_dot_b[1];
      n[2] = sp_dot_b[2];
    }
    v[0] = v0[0] + vdw * sp_dot_b[0];
    v[1] = v0[1] + vdw * sp_dot_b[1];
    v[2] = v0[2] + vdw * sp_dot_b[2];
    i = *(MapLocusEStart(map, v));
    if (i) {
      int j = map->EList[i++];
      while (ok && j >= 0) {
        SurfaceJobAtomInfo* j_atom_info = atom_info + j;
        if ((!present) || present[j]) {
          if (j != a) {
            int skip_flag = false;
            if (j_atom_info->vdw ==
                a_atom_info->vdw) { /* handle singularities */
              float* v1 = coord + 3 * j;
              if ((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2])) {
                skip_flag = true;
              }
            }
            if (!skip_flag)
              if (within3f(coord + 3 * j, v, j_atom_info->vdw + radius)) {
                flag = false;
                break;
              }
          }
        }
        j = map->EList[i++];
        ok &= !G->Interrupt;
      }
    }
    if (ok && flag && (*dotCnt < stopDot)) {
      (*dotCnt)++;
      v += 3;
      if (n) n += 3;
      (*nDot)++;
    }
  }
  return ok;
}

static int SolventDotCircumscribeAroundVertex(
    PyMOLGlobals* G, SolventDot* I, MapType* map, float* vdw, float dist,
    float* v0, float* v2, int circumscribe, SurfaceJobAtomInfo* atom_info,
    SurfaceJobAtomInfo* a_atom_info, SurfaceJobAtomInfo* jj_atom_info,
    int* present, int a, int jj, float* coord, float probe_radius, int* dotCnt,
    int stopDot) {
  int ok = true;
  float vz[3], vx[3], vy[3], vp[3];
  float tri_a = vdw[0], tri_b = vdw[2], tri_c = dist;
  float tri_s = (tri_a + tri_b + tri_c) * 0.5F;
  float area = (float)sqrt1f(tri_s * (tri_s - tri_a) * (tri_s - tri_b) *
                             (tri_s - tri_c));
  float radius = (2 * area) / dist;
  float adj = (float)sqrt1f(vdw[1] - radius * radius);
  int b;
  float* v = I->dot + I->nDot * 3;
  float* n = I->dotNormal + I->nDot * 3;

  subtract3f(v2, v0, vz);
  get_system1f3f(vz, vx, vy);

  copy3f(vz, vp);
  scale3f(vp, adj, vp);
  add3f(v0, vp, vp);

  for (b = 0; ok && b <= circumscribe; b++) {
    float xcos = (float)cos((b * 2 * cPI) / circumscribe);
    float ysin = (float)sin((b * 2 * cPI) / circumscribe);
    float xcosr = xcos * radius;
    float ysinr = ysin * radius;
    int flag = true, i;
    v[0] = vp[0] + vx[0] * xcosr + vy[0] * ysinr;
    v[1] = vp[1] + vx[1] * xcosr + vy[1] * ysinr;
    v[2] = vp[2] + vx[2] * xcosr + vy[2] * ysinr;

    i = *(MapLocusEStart(map, v));
    if (i && !map->EList.empty()) {
      int j = map->EList[i++];
      while (ok && j >= 0) {
        SurfaceJobAtomInfo* j_atom_info = atom_info + j;
        if ((!present) || present[j])
          if ((j != a) && (j != jj)) {
            int skip_flag = false;
            if (a_atom_info->vdw ==
                j_atom_info->vdw) { /* handle singularities */
              float* v1 = coord + 3 * j;
              if ((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                skip_flag = true;
            }
            if (jj_atom_info->vdw ==
                j_atom_info->vdw) { /* handle singularities */
              float* v1 = coord + 3 * j;
              if ((v2[0] == v1[0]) && (v2[1] == v1[1]) && (v2[2] == v1[2]))
                skip_flag = true;
            }
            if (!skip_flag)
              if (within3f(coord + 3 * j, v, j_atom_info->vdw + probe_radius)) {
                flag = false;
                break;
              }
          }
        j = map->EList[i++];
        ok &= !G->Interrupt;
      }
    }
    if (ok && flag && (*dotCnt < stopDot)) {
      float vt0[3], vt2[3];
      subtract3f(v0, v, vt0);
      subtract3f(v2, v, vt2);
      normalize3f(vt0);
      normalize3f(vt2);
      add3f(vt0, vt2, n);
      invert3f(n);
      normalize3f(n);
      /*
        n[0] = vx[0] * xcos + vy[0] * ysin;
        n[1] = vx[1] * xcos + vy[1] * ysin;
        n[2] = vx[2] * xcos + vy[2] * ysin;
      */
      I->dotCode[I->nDot] = 1; /* mark as exempt */
      (*dotCnt)++;
      v += 3;
      n += 3;
      I->nDot++;
    }
  }
  return ok;
}

static int SolventDotMarkDotsWithinCutoff(PyMOLGlobals* G, SolventDot* I,
                                          MapType* map, float* I_dot, int nDot,
                                          float* cavityDot, int* dot_flag,
                                          float cutoff) {
  int ok = true, *p = dot_flag;
  float* v = I->dot;
  int a;
  for (a = 0; ok && a < I->nDot; a++) {
    int i = *(MapLocusEStart(map, v));
    if (i && !map->EList.empty()) {
      int j = map->EList[i++];
      while (j >= 0) {
        if (within3f(cavityDot + (3 * j), v, cutoff)) {
          *p = true;
          break;
        }
        j = map->EList[i++];
      }
    }
    v += 3;
    p++;
    if (G->Interrupt) {
      ok = false;
    }
  }
  return ok;
}

static void SolventDotSlideDotsAndInfo(PyMOLGlobals* G, SolventDot* I,
                                       int* dot_flag, int flag_value) {
  float* v0 = I->dot;
  float* n0 = I->dotNormal;
  int* dc0 = I->dotCode;
  int* p = dot_flag;
  int c = I->nDot;
  float* n = n0;
  float* v = v0;
  int* dc = dc0;
  int a;
  I->nDot = 0;
  for (a = 0; a < c; a++) {
    if (flag_value ? *(p++) : !*(p++)) {
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(dc0++) = *(dc++);
      I->nDot++;
    } else {
      v += 3;
      n += 3;
    }
  }
  PRINTFD(G, FB_RepSurface)
  " SolventDotNew-DEBUG: %d->%d\n", c, I->nDot ENDFD;
}

static int SolventDotMarkDotsWithinProbeRadius(PyMOLGlobals* G, SolventDot* I,
                                               MapType* map, int cavity_cull,
                                               float probe_radius_plus,
                                               int* dot_flag, int* flag) {
  int ok = true;
  int* p = dot_flag;
  float* v = I->dot;
  int a;
  for (a = 0; ok && a < I->nDot; a++) {
    if (!dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      int cnt = 0;
      if (i && !map->EList.empty()) {
        int j = map->EList[i++];
        while (j >= 0) {
          if (j != a) {
            if (within3f(I->dot + (3 * j), v, probe_radius_plus)) {
              if (dot_flag[j]) {
                *p = true;
                (*flag) = true;
                break;
              }
              cnt++;
              if (cnt > cavity_cull) {
                *p = true;
                (*flag) = true;
                break;
              }
            }
          }
          j = map->EList[i++];
        }
      }
    }
    v += 3;
    p++;
    ok &= !G->Interrupt;
  }
  ok &= !G->Interrupt;
  return ok;
}

/* ===================================================================
 * Prototype: parallelized solvent-dot sampling (env-flagged)
 *
 * The per-atom solvent-dot sampling loop is embarrassingly parallel: every
 * atom samples a fixed set of sphere points and tests each against a
 * read-only neighbor map.  The only serial dependency in the original code
 * is that accepted dots are appended into shared output arrays in atom
 * order.  The parallel path below preserves that exact ordering (each thread
 * owns a contiguous atom range and writes to a thread-local buffer; buffers
 * are merged in thread/atom order), so it produces bit-identical output to
 * the serial path and therefore an identical downstream triangulation.
 *
 * Controlled at runtime (no rebuild needed to A/B):
 *   PYMOL_SURFACE_PARALLEL = unset|0|off  -> serial reference (default)
 *                          = 1|on         -> parallel
 *                          = verify        -> run both, compare, use serial
 *   PYMOL_SURFACE_THREADS  = N            -> override thread count
 * =================================================================== */

/* Serial reference implementation of the per-atom solvent-dot sampling loop,
 * factored out of SolventDotNew so the verify path can reuse it. Writes
 * accepted dots/normals into the caller-provided arrays. */
static int SolventDotComputeMainLoopSerial(
    PyMOLGlobals* G, SolventDot* I, MapType* map,
    std::vector<SurfaceJobAtomInfo>& atom_info, float* coord, int n_coord,
    int* present, SphereRec* sp, float probe_radius, int* dotCnt, int stopDot,
    float* dotOut, float* normalOut, int* nDotOut) {
  int ok = true;
  SurfaceJobAtomInfo* a_atom_info = atom_info.data();
  for (int a = 0; ok && a < n_coord; a++) {
    OrthoBusyFast(G, a, n_coord * 5);
    if ((!present) || (present[a])) {
      int skip_flag = false;
      ok = SolventDotFilterOutSameXYZ(G, map, atom_info.data(), a_atom_info,
                                      coord, a, present, &skip_flag);
      if (ok && !skip_flag) {
        ok = SolventDotGetDotsAroundVertexInSphere(
            G, I, map, atom_info.data(), a_atom_info, coord, a, present, sp,
            probe_radius, dotCnt, stopDot, dotOut, normalOut, nDotOut);
      }
    }
    a_atom_info++;
  }
  return ok;
}

/* Parallel implementation. Each thread processes a contiguous range of atoms
 * into thread-local buffers, which are then concatenated in thread order.
 * With explicit contiguous ranges and dynamic threads disabled, the merged
 * order is exactly the serial atom order. */
static int SolventDotComputeMainLoopParallel(
    PyMOLGlobals* G, SolventDot* I, MapType* map,
    std::vector<SurfaceJobAtomInfo>& atom_info, float* coord, int n_coord,
    int* present, SphereRec* sp, float probe_radius, int* dotCnt, int stopDot,
    float* dotOut, float* normalOut, int* nDotOut, int nthreads) {
  if (nthreads < 1) nthreads = 1;
  if (nthreads > n_coord) nthreads = (n_coord > 0) ? n_coord : 1;

  std::vector<std::vector<float>> tDot(nthreads);
  std::vector<std::vector<float>> tNorm(nthreads);

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    int nt = omp_get_num_threads();
    int chunk = (n_coord + nt - 1) / nt;
    int lo = tid * chunk;
    int hi = lo + chunk;
    if (hi > n_coord) hi = n_coord;

    std::vector<float>& vd = tDot[tid];
    std::vector<float>& vn = tNorm[tid];
    /* per-atom scratch, reused across this thread's atoms */
    std::vector<float> sd(sp->nDot * 3);
    std::vector<float> sn(sp->nDot * 3);
    SurfaceJobAtomInfo* ainfo = atom_info.data();

    for (int a = lo; a < hi; a++) {
      if (G->Interrupt) continue;
      if ((!present) || (present[a])) {
        SurfaceJobAtomInfo* a_atom_info = ainfo + a;
        int skip_flag = false;
        SolventDotFilterOutSameXYZ(G, map, atom_info.data(), a_atom_info, coord,
                                   a, present, &skip_flag);
        if (!skip_flag) {
          int localN = 0, localCnt = 0;
          /* identical per-atom computation as the serial path */
          SolventDotGetDotsAroundVertexInSphere(
              G, I, map, atom_info.data(), a_atom_info, coord, a, present, sp,
              probe_radius, &localCnt, sp->nDot, sd.data(), sn.data(), &localN);
          if (localN > 0) {
            vd.insert(vd.end(), sd.begin(), sd.begin() + localN * 3);
            vn.insert(vn.end(), sn.begin(), sn.begin() + localN * 3);
          }
        }
      }
    }
  }

  /* Merge thread-local buffers in atom order, reproducing the serial append
   * order and the global stopDot cap. */
  for (int t = 0; t < nthreads; t++) {
    int navail = (int)(tDot[t].size() / 3);
    int room = stopDot - *dotCnt;
    int take = (navail < room) ? navail : room;
    if (take > 0) {
      memcpy(dotOut + (*nDotOut) * 3, tDot[t].data(),
             (size_t)take * 3 * sizeof(float));
      memcpy(normalOut + (*nDotOut) * 3, tNorm[t].data(),
             (size_t)take * 3 * sizeof(float));
      *nDotOut += take;
      *dotCnt += take;
    }
  }
  return true;
}

/* Dispatch the solvent-dot main loop based on the runtime flag. */
static int SolventDotMainLoop(PyMOLGlobals* G, SolventDot* I, MapType* map,
                              std::vector<SurfaceJobAtomInfo>& atom_info,
                              float* coord, int n_coord, int* present,
                              SphereRec* sp, float probe_radius, int* dotCnt,
                              int stopDot) {
  int mode = 0; /* 0 = serial, 1 = parallel, 2 = verify */
  /* Parallel by default when OpenMP is compiled in. PYMOL_SURFACE_PARALLEL
   * overrides: "0"/"off" forces serial, "verify" runs both and compares. */
  mode = 1;
  int nthreads = omp_get_max_threads();
  {
    const char* e = getenv("PYMOL_SURFACE_PARALLEL");
    if (e && *e) {
      if (!strcmp(e, "verify"))
        mode = 2;
      else if (!strcmp(e, "0") || !strcmp(e, "off"))
        mode = 0;
      else
        mode = 1;
    }
    const char* tn = getenv("PYMOL_SURFACE_THREADS");
    if (tn && atoi(tn) > 0) nthreads = atoi(tn);
  }

  const bool timing = getenv("PYMOL_SURFACE_TIMING") != nullptr;
  double t0 = timing ? UtilGetSeconds(G) : 0.0;

  if (mode == 0) {
    int rc = SolventDotComputeMainLoopSerial(
        G, I, map, atom_info, coord, n_coord, present, sp, probe_radius, dotCnt,
        stopDot, I->dot, I->dotNormal, &I->nDot);
    if (timing) {
      fprintf(stderr,
              "SURFACE_TIMING: stage=dot_loop mode=serial threads=1 "
              "ms=%.3f\n",
              (UtilGetSeconds(G) - t0) * 1000.0);
      fflush(stderr);
    }
    return rc;
  }
  if (mode == 1) {
    int rc = SolventDotComputeMainLoopParallel(
        G, I, map, atom_info, coord, n_coord, present, sp, probe_radius, dotCnt,
        stopDot, I->dot, I->dotNormal, &I->nDot, nthreads);
    if (timing) {
      fprintf(stderr,
              "SURFACE_TIMING: stage=dot_loop mode=parallel "
              "threads=%d ms=%.3f\n",
              nthreads, (UtilGetSeconds(G) - t0) * 1000.0);
      fflush(stderr);
    }
    return rc;
  }

  /* mode == 2: verify — run both into separate buffers and compare. The
   * serial result is used downstream so behavior is unchanged in verify. */
  {
    float* parDot = VLAlloc(float, (stopDot + 1) * 3);
    float* parNorm = VLAlloc(float, (stopDot + 1) * 3);
    int parN = 0, parCnt = 0;
    int ok = true;

    ok &= SolventDotComputeMainLoopParallel(
        G, I, map, atom_info, coord, n_coord, present, sp, probe_radius,
        &parCnt, stopDot, parDot, parNorm, &parN, nthreads);
    ok &= SolventDotComputeMainLoopSerial(
        G, I, map, atom_info, coord, n_coord, present, sp, probe_radius, dotCnt,
        stopDot, I->dot, I->dotNormal, &I->nDot);

    int nmism = 0;
    double maxd = 0.0;
    if (parN != I->nDot) {
      fprintf(
          stderr,
          "SURFACE_VERIFY: COUNT MISMATCH serial=%d parallel=%d threads=%d\n",
          I->nDot, parN, nthreads);
    } else {
      for (int k = 0; k < I->nDot * 3; k++) {
        double d = fabs((double)I->dot[k] - (double)parDot[k]);
        double dn = fabs((double)I->dotNormal[k] - (double)parNorm[k]);
        if (d > maxd) maxd = d;
        if (dn > maxd) maxd = dn;
        if (I->dot[k] != parDot[k] || I->dotNormal[k] != parNorm[k]) nmism++;
      }
      fprintf(stderr,
              "SURFACE_VERIFY: nDot=%d threads=%d mismatched_floats=%d "
              "max_abs_diff=%g %s\n",
              I->nDot, nthreads, nmism, maxd,
              (nmism == 0) ? "IDENTICAL" : "DIFFER");
    }
    fflush(stderr);
    VLAFreeP(parDot);
    VLAFreeP(parNorm);
    return ok;
  }
  return true;
}

SolventDot* SolventDotNew(PyMOLGlobals* G, float* coord,
                          std::vector<SurfaceJobAtomInfo>& atom_info,
                          float probe_radius, SphereRec* sp, int* present,
                          int circumscribe, int surface_mode,
                          int surface_solvent, int cavity_cull,
                          int all_visible_flag, float max_vdw, int cavity_mode,
                          float cavity_radius, float cavity_cutoff) {
  int ok = true;
  int stopDot;
  int n_coord = atom_info.size();
  auto I = new SolventDot();
  CHECKOK(ok, I);
  /*  printf("%p %p %p %f %p %p %p %d %d %d %d %d %f\n",
     G,
     coord,
     atom_info,
     probe_radius,sp,
     extent,present,
     circumscribe,  surface_mode,
     surface_solvent,  cavity_cull,
     all_visible_flag, max_vdw);
   */

  stopDot = n_coord * sp->nDot + 2 * circumscribe;
  if (ok) I->dot = VLAlloc(float, (stopDot + 1) * 3);
  CHECKOK(ok, I->dot);
  if (ok) {
    I->dotNormal = VLAlloc(float, (stopDot + 1) * 3);
    CHECKOK(ok, I->dotNormal);
  }
  if (ok) {
    I->dotCode = VLACalloc(int, stopDot + 1);
    CHECKOK(ok, I->dotCode);
  }
  I->nDot = 0;
  if (ok) {
    int dotCnt = 0;
    MapType* map = new MapType(G, max_vdw + probe_radius, coord, n_coord,
                               nullptr, present);
    CHECKOK(ok, map);
    ok &= !G->Interrupt;
    if (map && ok) {
      ok &= MapSetupExpress(map);
      if (ok) {
        ok = SolventDotMainLoop(G, I, map, atom_info, coord, n_coord, present,
                                sp, probe_radius, &dotCnt, stopDot);
      }

      /* for each pair of proximal atoms, circumscribe a circle for their
       * intersection */
      if (ok) {
        MapType* map2 = nullptr;
        if (circumscribe && (!surface_solvent)) {
          map2 = new MapType(G, 2 * (max_vdw + probe_radius), coord, n_coord,
                             nullptr, present);
          CHECKOK(ok, map2);
        }
        ok &= !G->Interrupt;
        if (ok && map2) {
          int a;
          int skip_flag;

          SurfaceJobAtomInfo* a_atom_info = atom_info.data();
          ok &= MapSetupExpress(map2);
          for (a = 0; ok && a < n_coord; a++) {
            if ((!present) || present[a]) {
              float* v0 = coord + 3 * a;
              skip_flag = false;

              ok = SolventDotFilterOutSameXYZ(G, map2, atom_info.data(),
                                              a_atom_info, coord, a, present,
                                              &skip_flag);
              if (ok && !skip_flag) {
                int ii = *(MapLocusEStart(map2, v0));
                if (ii) {
                  int jj = map2->EList[ii++];
                  float vdw[3];
                  vdw[0] = a_atom_info->vdw + probe_radius;
                  vdw[1] = vdw[0] * vdw[0];
                  while (ok && jj >= 0) {
                    auto* jj_atom_info = &atom_info[jj];
                    float dist;
                    if (jj > a) /* only check if this is atom trails */
                      if ((!present) || present[jj]) {
                        float* v2 = coord + 3 * jj;
                        vdw[2] = jj_atom_info->vdw + probe_radius;
                        dist = (float)diff3f(v0, v2);
                        if ((dist > R_SMALL4) && (dist < (vdw[0] + vdw[2]))) {
                          ok = SolventDotCircumscribeAroundVertex(
                              G, I, map, vdw, dist, v0, v2, circumscribe,
                              atom_info.data(), a_atom_info, jj_atom_info,
                              present, a, jj, coord, probe_radius, &dotCnt,
                              stopDot);
                        }
                      }
                    jj = map2->EList[ii++];
                  }
                }
              }
            }
            a_atom_info++;
            ok &= !G->Interrupt;
          }
        }
        MapFree(map2);
      }
    }
    MapFree(map);
  }

  if (ok && cavity_mode) {
    int nCavityDot = 0;
    int dotCnt = 0;
    float* cavityDot = VLAlloc(float, (stopDot + 1) * 3);
    CHECKOK(ok, cavityDot);
    if (cavity_radius < 0.0F) {
      cavity_radius = -probe_radius * cavity_radius;
    }
    if (cavity_cutoff < 0.0F) {
      cavity_cutoff = cavity_radius - cavity_cutoff * probe_radius;
    }
    {
      MapType* map = new MapType(G, max_vdw + cavity_radius, coord, n_coord,
                                 nullptr, present);
      CHECKOK(ok, map);
      if (G->Interrupt) ok = false;
      if (ok && map) {
        ok &= MapSetupExpress(map);
        if (ok) {
          int a;
          int skip_flag;
          SurfaceJobAtomInfo* a_atom_info = atom_info.data();
          for (a = 0; a < n_coord; a++) {
            if ((!present) || (present[a])) {
              skip_flag = false;
              ok = SolventDotFilterOutSameXYZ(G, map, atom_info.data(),
                                              a_atom_info, coord, a, present,
                                              &skip_flag);
              if (ok && !skip_flag) {
                ok = SolventDotGetDotsAroundVertexInSphere(
                    G, I, map, atom_info.data(), a_atom_info, coord, a, present,
                    sp, cavity_radius, &dotCnt, stopDot, cavityDot, nullptr,
                    &nCavityDot);
              }
            }
            a_atom_info++;
          }
        }
      }
      MapFree(map);
    }
    {
      int* dot_flag = pymol::calloc<int>(I->nDot);
      ErrChkPtr(G, dot_flag);
      {
        MapType* map =
            new MapType(G, cavity_cutoff, cavityDot, nCavityDot, nullptr);
        if (map) {
          MapSetupExpress(map);
          ok = SolventDotMarkDotsWithinCutoff(
              G, I, map, I->dot, I->nDot, cavityDot, dot_flag, cavity_cutoff);
        }
        MapFree(map);
      }
      SolventDotSlideDotsAndInfo(G, I, dot_flag, false);
      FreeP(dot_flag);
    }
    VLAFreeP(cavityDot);
  }

  if (ok && (cavity_mode != 1) && (cavity_cull > 0) && (probe_radius > 0.75F) &&
      (!surface_solvent)) {
    int* dot_flag = pymol::calloc<int>(I->nDot);
    float probe_radius_plus;
    probe_radius_plus = probe_radius * 1.5F;

    ErrChkPtr(G, dot_flag);

    {
      MapType* map =
          new MapType(G, probe_radius_plus, I->dot, I->nDot, nullptr);
      if (map) {
        int flag = true;
        MapSetupExpress(map);
        while (ok && flag) {
          flag = false;
          ok = SolventDotMarkDotsWithinProbeRadius(
              G, I, map, cavity_cull, probe_radius_plus, dot_flag, &flag);
        }
      }
      MapFree(map);
    }

    if (ok) {
      SolventDotSlideDotsAndInfo(G, I, dot_flag, true);
    }

    FreeP(dot_flag);
  }
  if (!ok) {
    SolventDotFree(I);
    I = nullptr;
  }
  return I;
}
