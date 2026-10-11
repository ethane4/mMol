
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
-* Cameron Mura
-*
-*
Z* -------------------------------------------------------------------
*/

// The cartoon guide path: putty values, normals, tangents, helices, sheets and
// loops.

#include "pymol/reps/rep_cartoon.h"
#include "pymol/base/feedback.h"
#include "pymol/reps/rep_cartoon_internal.h"

void RepCartoonComputePuttyValues(ObjectMolecule* obj, float* putty_vals) {
  double sum = 0.0, sumsq = 0.0;
  float value;
  int cnt = 0;
  AtomInfoType* ai;
  int a;
  for (a = 0; a < obj->NAtom; a++) {
    ai = obj->AtomInfo + a;

    if (ai->visRep & cRepCartoonBit) {
      value = ai->b;
      sum += value;
      sumsq += (value * value);
      if (value < putty_vals[2]) putty_vals[2] = value;
      if (value > putty_vals[3]) putty_vals[3] = value;
      cnt++;
    }
  }

  if (cnt) {
    putty_vals[0] = (float)(sum / cnt);
    putty_vals[1] = (float)sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
  } else {
    /* aren't these assignments unnecessary? */
    putty_vals[0] = 10.0F;
    putty_vals[1] = 10.0F;
    putty_vals[2] = 0.0F;
    putty_vals[3] = 10.0F;
  }
}

/**
 * @brief Detect cyclic segments after PASS1 trace collection.
 *
 * Scans segment boundaries and uses ObjectMoleculeCheckBondSep to determine
 * if the last CA in a segment is connected back to the first CA through
 * 3 bonds (CA-C-N-CA), indicating a cyclic peptide backbone.
 *
 * @param obj    The molecule object (used for bond connectivity check)
 * @param cs     The coordinate set
 * @param nAt    Number of atoms in the trace arrays
 * @param seg    Segment ID array (per trace atom)
 * @param at     CoordSet atom index array (per trace atom)
 * @return Map from segment ID to CyclicSegmentInfo for each cyclic segment
 */
std::unordered_map<int, CyclicSegmentInfo> RepCartoonDetectCyclicSegments(
    ObjectMolecule* obj, CoordSet* cs, int nAt, const int* seg, const int* at) {
  std::unordered_map<int, CyclicSegmentInfo> result;
  if (nAt < 3) return result;

  // Find first and last index for each segment
  std::unordered_map<int, CyclicSegmentInfo> segBounds;
  for (int a = 0; a < nAt; a++) {
    int segId = seg[a];
    auto it = segBounds.find(segId);
    if (it == segBounds.end()) {
      segBounds[segId] = {a, a};
    } else {
      it->second.lastIdx = a;
    }
  }

  // Check each segment for cyclic connectivity (CA-C-N-CA = 3 bonds)
  for (auto& [segId, info] : segBounds) {
    int segLen = info.lastIdx - info.firstIdx + 1;
    if (segLen < 3) continue;

    int firstAtom = cs->IdxToAtm[at[info.firstIdx]];
    int lastAtom = cs->IdxToAtm[at[info.lastIdx]];

    if (ObjectMoleculeCheckBondSep(obj, firstAtom, lastAtom, 3)) {
      assert(info.firstIdx >= 0 && info.firstIdx < nAt);
      assert(info.lastIdx >= 0 && info.lastIdx < nAt);
      assert(info.firstIdx < info.lastIdx);
      result[segId] = info;
    }
  }

  return result;
}

void RepCartoonComputeDifferencesAndNormals(PyMOLGlobals* G, int nAt, int* seg,
                                            float* pv, float* dv, float* nv,
                                            float* dl, int quiet) {
  float *v1, *v2, *vptr, *d;
  int *sptr, a;
  /* compute differences and normals */
  sptr = seg;
  vptr = pv;
  v1 = dv;
  v2 = nv;
  d = dl;
  for (a = 0; a < (nAt - 1); a++) {
    if (!quiet) {
      PRINTFD(G, FB_RepCartoon)
      " RepCartoon: seg %d *s %d , *(s+1) %d\n", a, *sptr, *(sptr + 1) ENDFD;
    }
    if (*sptr == *(sptr + 1)) {
      subtract3f(vptr + 3, vptr, v1);
      *d = (float)length3f(v1);
      if (*d > R_SMALL4) {
        float d_1;
        d_1 = 1.0F / (*d);
        scale3f(v1, d_1, v2);
      } else if (a) { /* if zero, copy previous */
        copy3f(v2 - 3, v2);
      } else {
        zero3f(v2);
      }
    } else {
      zero3f(v2);
    }
    d++;
    vptr += 3;
    v1 += 3;
    v2 += 3;
    sptr++;
  }
}

void RepCartoonComputeTangents(int nAt, int* seg, float* nv, float* tv) {
  float *vptr, *v1;
  int *sptr, a;
  sptr = seg;
  vptr = nv;
  v1 = tv;

  *(v1++) = *(vptr++); /* first segment */
  *(v1++) = *(vptr++);
  *(v1++) = *(vptr++);
  sptr++;

  for (a = 1; a < (nAt - 1); a++) {
    if ((*sptr == *(sptr - 1)) && (*sptr == *(sptr + 1))) {
      add3f(vptr, (vptr - 3),
            v1); /* tangent vectors are head-to-tail sums within a segment */
      normalize3f(v1);
    } else if (*sptr == *(sptr - 1)) {
      *(v1) = *(vptr - 3); /* end a segment */
      *(v1 + 1) = *(vptr - 2);
      *(v1 + 2) = *(vptr - 1);
    } else if (*sptr == *(sptr + 1)) {
      *(v1) = *(vptr); /* new segment */
      *(v1 + 1) = *(vptr + 1);
      *(v1 + 2) = *(vptr + 2);
    }
    vptr += 3;
    v1 += 3;
    sptr++;
  }
  *(v1++) = *(vptr - 3); /* last segment */
  *(v1++) = *(vptr - 2);
  *(v1++) = *(vptr - 1);
}

/**
 * @brief Compute difference vectors, normals, and distances with cyclic
 * wrap-around.
 *
 * Like RepCartoonComputeDifferencesAndNormals, but for the last atom in a
 * cyclic segment, computes the difference vector from lastIdx to firstIdx
 * (wrap-around) instead of leaving it zeroed. This provides correct spline
 * tangent data at the cyclic closure point.
 *
 * @param G          PyMOL globals
 * @param nAt        Number of trace atoms
 * @param seg        Segment ID array
 * @param pv         Positions array (3f per atom)
 * @param dv         [out] Difference vectors (3f per atom)
 * @param nv         [out] Normalized direction vectors (3f per atom)
 * @param dl         [out] Distance between consecutive atoms (1f per atom)
 * @param quiet      If false, print debug feedback
 * @param cyclicSegs Map of cyclic segment info (from
 * RepCartoonDetectCyclicSegments)
 */
void RepCartoonComputeDifferencesAndNormalsCyclic(
    PyMOLGlobals* G, int nAt, int* seg, float* pv, float* dv, float* nv,
    float* dl, int quiet,
    const std::unordered_map<int, CyclicSegmentInfo>& cyclicSegs) {
  float *v1, *v2, *vptr, *d;
  int *sptr, a;
  sptr = seg;
  vptr = pv;
  v1 = dv;
  v2 = nv;
  d = dl;
  for (a = 0; a < (nAt - 1); a++) {
    if (*sptr == *(sptr + 1)) {
      subtract3f(vptr + 3, vptr, v1);
      *d = (float)length3f(v1);
      if (*d > R_SMALL4) {
        float d_1 = 1.0F / (*d);
        scale3f(v1, d_1, v2);
      } else if (a) {
        copy3f(v2 - 3, v2);
      } else {
        zero3f(v2);
      }
    } else {
      // Segment boundary — check if last atom of a cyclic segment
      auto it = cyclicSegs.find(*sptr);
      if (it != cyclicSegs.end() && a == it->second.lastIdx) {
        // Wrap: difference from last to first atom of segment
        float* firstPos = pv + it->second.firstIdx * 3;
        subtract3f(firstPos, vptr, v1);
        *d = (float)length3f(v1);
        if (*d > R_SMALL4) {
          float d_1 = 1.0F / (*d);
          scale3f(v1, d_1, v2);
        } else {
          copy3f(v2 - 3, v2);
        }
      } else {
        zero3f(v2);
      }
    }
    d++;
    vptr += 3;
    v1 += 3;
    v2 += 3;
    sptr++;
  }

  // Handle last atom (a = nAt - 1) for cyclic segments
  {
    auto it = cyclicSegs.find(*sptr);
    if (it != cyclicSegs.end() && (nAt - 1) == it->second.lastIdx) {
      float* firstPos = pv + it->second.firstIdx * 3;
      subtract3f(firstPos, vptr, v1);
      *d = (float)length3f(v1);
      if (*d > R_SMALL4) {
        float d_1 = 1.0F / (*d);
        scale3f(v1, d_1, v2);
      } else {
        copy3f(v2 - 3, v2);
      }
    }
  }
}

/**
 * @brief Compute tangent vectors with cyclic wrap-around at segment boundaries.
 *
 * Like RepCartoonComputeTangents, but at the first and last atoms of cyclic
 * segments, averages the tangent with the wrap-around neighbor's direction
 * vector instead of simply copying the endpoint direction.
 *
 * @param nAt        Number of trace atoms
 * @param seg        Segment ID array
 * @param nv         Normalized direction vectors (input, from differences pass)
 * @param tv         [out] Tangent vectors (3f per atom)
 * @param cyclicSegs Map of cyclic segment info
 */
void RepCartoonComputeTangentsCyclic(
    int nAt, int* seg, float* nv, float* tv,
    const std::unordered_map<int, CyclicSegmentInfo>& cyclicSegs) {
  float *vptr, *v1;
  int *sptr, a;
  sptr = seg;
  vptr = nv;
  v1 = tv;

  // First atom
  {
    auto it = cyclicSegs.find(*sptr);
    if (it != cyclicSegs.end() && it->second.firstIdx == 0) {
      // Cyclic: average with the last normal of this segment
      float* lastNv = nv + it->second.lastIdx * 3;
      add3f(vptr, lastNv, v1);
      normalize3f(v1);
      v1 += 3;
      vptr += 3;
      sptr++;
    } else {
      *(v1++) = *(vptr++);
      *(v1++) = *(vptr++);
      *(v1++) = *(vptr++);
      sptr++;
    }
  }

  for (a = 1; a < (nAt - 1); a++) {
    if ((*sptr == *(sptr - 1)) && (*sptr == *(sptr + 1))) {
      add3f(vptr, (vptr - 3), v1);
      normalize3f(v1);
    } else if (*sptr == *(sptr - 1)) {
      // End of segment
      auto it = cyclicSegs.find(*sptr);
      if (it != cyclicSegs.end() && a == it->second.lastIdx) {
        // Cyclic: average current with wrap-around normal
        add3f(vptr - 3, nv + it->second.lastIdx * 3, v1);
        normalize3f(v1);
      } else {
        copy3f(vptr - 3, v1);
      }
    } else if (*sptr == *(sptr + 1)) {
      // Start of segment
      auto it = cyclicSegs.find(*sptr);
      if (it != cyclicSegs.end() && a == it->second.firstIdx) {
        float* lastNv = nv + it->second.lastIdx * 3;
        add3f(vptr, lastNv, v1);
        normalize3f(v1);
      } else {
        copy3f(vptr, v1);
      }
    }
    vptr += 3;
    v1 += 3;
    sptr++;
  }

  // Last atom
  {
    auto it = cyclicSegs.find(*sptr);
    if (it != cyclicSegs.end() && (nAt - 1) == it->second.lastIdx) {
      add3f(vptr - 3, nv + it->second.lastIdx * 3, v1);
      normalize3f(v1);
    } else {
      *(v1++) = *(vptr - 3);
      *(v1++) = *(vptr - 2);
      *(v1++) = *(vptr - 1);
    }
  }
}

void RepCartoonComputeRoundHelices(nuc_acid_data* ndata, const int nAt,
                                   const int* sptr, const ss_t* ss,
                                   const float* v0, const float* vptr) {
  const float *v1 = nullptr, *v2 = nullptr, *v3 = nullptr, *v4 = nullptr,
              *v5 = nullptr;
  int last, a;
  float t0[3], t1[3], t2[3];
  last = 0;
  if (nAt > 1) {
    for (a = 0; a < nAt; a++) {
      if (a) {
        if (*sptr !=
            *(sptr - 1)) { /* contiguous helices in disconnected segments */
          v1 = nullptr;
          v2 = nullptr;
          v3 = nullptr;
          v4 = nullptr;
          v5 = nullptr;
          last = 0;
        }
      }
      v5 = v4;
      v4 = v3;
      v3 = v2;
      v2 = v1;
      if (*ss == ss_t::HELIX)
        v1 = vptr;
      else { /* early termination ? */
        if (last < 2) {
          zero3f(t0);
          if (v2 && v3) {
            subtract3f(v2, vptr, t0);
            normalize3f(t0);
            subtract3f(v3, v2, t1);
            normalize3f(t1);
            add3f(t1, t0, t0);
            if (v4) {
              subtract3f(v4, v3, t1);
              normalize3f(t1);
              add3f(t1, t0, t0);
            }
            if (v5) {
              subtract3f(v5, v4, t1);
              normalize3f(t1);
              add3f(t1, t0, t0);
            }
            normalize3f(t0);
            cross_product3f(t0, v0 - 3, ndata->voptr - 3);
            normalize3f(ndata->voptr - 3);
            cross_product3f(t0, v0 - 6, ndata->voptr - 6);
            normalize3f(ndata->voptr - 6);
            if (v4) {
              cross_product3f(t0, v0 - 9, ndata->voptr - 9);
              normalize3f(ndata->voptr - 9);
            }
            if (v5) {
              cross_product3f(t0, v0 - 12, ndata->voptr - 12);
              normalize3f(ndata->voptr - 12);
            }
            if (v4 && v5) {
              /* now make sure there's no goofy flip on the end...
                 of a short, tight helix */
              if (dot_product3f(ndata->voptr - 9, ndata->voptr - 12) < -0.8F)
                invert3f(ndata->voptr - 12);
            }
          }
        }
        v1 = nullptr;
        v2 = nullptr;
        v3 = nullptr;
        v4 = nullptr;
        v5 = nullptr;
        last = 0;
      }
      if (v1 && v2 && v3 && v4) {
        add3f(v1, v4, t0);
        add3f(v2, v3, t1);
        scale3f(t0, 0.2130F, t0);
        scale3f(t1, 0.2870F, t1);

        add3f(t0, t1, t0);
        if (last) { /* 5th CA or later... */
          subtract3f(t2, t0, t1);
          normalize3f(t1);
          cross_product3f(t1, v0, ndata->voptr);
          normalize3f(ndata->voptr);
          cross_product3f(t1, v0 - 3, ndata->voptr - 3);
          normalize3f(ndata->voptr - 3);
          cross_product3f(t1, v0 - 6, ndata->voptr - 6);
          normalize3f(ndata->voptr - 6);
          if (last == 1) { /* 5th */
            cross_product3f(t1, v0 - 9, ndata->voptr - 9);
            normalize3f(ndata->voptr - 9);
            cross_product3f(t1, v0 - 12, ndata->voptr - 12);
            normalize3f(ndata->voptr - 12);
          }
        }
        last++;
        copy3f(t0, t2);
      }
      vptr += 3;
      ss++;
      ndata->voptr += 3;
      v0 += 3;
      sptr++;
    }
  }
}

void RepCartoonRefineNormals(PyMOLGlobals* G, RepCartoon* I,
                             ObjectMolecule* obj, CoordSet* cs,
                             nuc_acid_data* ndata, const int nAt,
                             const int* seg, const float* tv, float* pvo,
                             float* pva, const ss_t* ss, const float* nv) {
  int refine_normals = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                    cSetting_cartoon_refine_normals);
  if (refine_normals < 0) {
    if (obj->NCSet > 1) {
      int i, n_set = 0;
      for (i = 0; i < obj->NCSet; i++)
        if (obj->CSet[i]) {
          n_set++;
          if (n_set > 1) refine_normals = 0;
          /* default behavior is to not refine normals for multi-state objects
           */
        }
    }
  }

  if (refine_normals) {
    /* first, make sure orientiation vectors are orthogonal to the tangent */
    float *va, max_dot;
    int a;
    float t0[3], t1[3], t2[3], t3[3], o0[12], o1[12];
    float dp;
    const float *v0, *v1 = tv + 3;
    ndata->voptr = pvo + 3;
    const int* sptr = seg + 1;
    for (a = 1; a < (nAt - 1); a++) {
      if ((*sptr == *(sptr - 1)) && (*sptr == *(sptr + 1))) {
        /* only operate on vectors within the cartoon itself --
           not the end vectors */
        remove_component3f(ndata->voptr, v1, t0);
        normalize23f(t0, ndata->voptr);
        /* go on to next vertex */
      }
      v1 += 3;
      ndata->voptr += 3;
      sptr++;
    }
    /* now generate alternative inverted orientation vectors */
    va = pva;
    ndata->voptr = pvo;
    for (a = 0; a < nAt; a++) {
      /* original */
      copy3f(ndata->voptr, va);
      va += 3;
      /* inverse */
      copy3f(ndata->voptr, va);
      if (*ss != ss_t::HELIX) {
        invert3f(va);
        /* for helix, don't allow inversion of normals, since that
           would confuse the inside & outside of the helix  */
      }
      va += 3;
      /* go on to next vertex */
      ndata->voptr += 3;
      ss++;
    }
    /* now iterate forward through pairs */
    ndata->voptr = pvo + 3;
    va = pva + 6;
    auto vptr = nv + 3; /* normals in direction of chain */
    sptr = seg + 1;
    for (a = 1; a < (nAt - 1); a++) {
      if ((*sptr == *(sptr + 1)) &&
          (*sptr == *(sptr - 1))) { /* only operate within a segment */
        remove_component3f(ndata->voptr - 3, vptr - 3,
                           o0); /* previous orientation vector */
        normalize3f(o0);        /* is now perp to chain direction */
        v1 = va; /* candidate orientation vectors for current CA */
        remove_component3f(
            v1, vptr - 3,
            o1); /* removes chain direction from the two candidates */
        remove_component3f(v1 + 3, vptr - 3, o1 + 3);
        normalize3f(o1);
        normalize3f(o1 + 3);
        max_dot = dot_product3f(o0, o1);
        v0 = v1;
        dp = dot_product3f(o0, o1 + 3);
        if (dp > max_dot) {
          v0 = v1 + 3;
          max_dot = dp;
        }
        copy3f(v0,
               ndata->voptr); /* updates atom with optimal orientation vector */
      }
      ndata->voptr += 3;
      va += 6;   /* candidate orientation vectors */
      vptr += 3; /* normal */
      sptr++;
    }
    /* now soften up the kinks */
    v1 = tv + 3;
    va = pva + 6;
    ndata->voptr = pvo + 3;
    sptr = seg + 1;
    for (a = 1; a < (nAt - 1); a++) {
      if ((*sptr == *(sptr - 1)) && (*sptr == *(sptr + 1))) {
        dp = (dot_product3f(ndata->voptr, ndata->voptr + 3) *
              dot_product3f(ndata->voptr, ndata->voptr - 3));
        if (dp < -0.10F) { /* threshold value -- could be a setting */
          add3f(ndata->voptr + 3, ndata->voptr - 3, t0);
          scale3f(ndata->voptr, 0.001, t1);
          add3f(t1, t0, t0);
          remove_component3f(t0, v1, t0);
          normalize3f(t0);
          if (dot_product3f(ndata->voptr, t0) < 0.0F) {
            subtract3f(ndata->voptr, t0, t2);
          } else {
            add3f(ndata->voptr, t0, t2);
          }
          normalize3f(t2);
          dp = 2 * (-0.10F - dp);
          if (dp > 1.0F) dp = 1.0F;
          mix3f(ndata->voptr, t2, dp, t3);
          copy3f(t3, va); /* store modified vector */
          invert3f3f(va, va + 3);
        } else {
          copy3f(ndata->voptr, va); /* keep as is */
        }
      }
      v1 += 3;
      ndata->voptr += 3;
      va += 6;
      sptr++;
    }
    /* now update */
    va = pva + 6;
    ndata->voptr = pvo + 3;
    sptr = seg + 1;
    for (a = 1; a < (nAt - 1); a++) {
      if ((*sptr == *(sptr - 1)) && (*sptr == *(sptr + 1))) {
        copy3f(va, ndata->voptr);
      }
      ndata->voptr += 3;
      va += 6;
      sptr++;
    }
  }
}

void RepCartoonFlattenSheets(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                             nuc_acid_data* ndata, int nAt, const int* sptr,
                             const CCInOut* cc, float* pv, float* pvo,
                             const ss_t* ss, const float* v0, float* tmp,
                             const int* flag_tmp) {
  int last, first, cur_car, end_flag, a, b, c, e, f;
  int flat_cycles;
  float t0[3];
  flat_cycles = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_flat_cycles);

  last = 0;
  first = -1;
  cur_car = *cc;
  end_flag = false;
  if (nAt > 1) {
    for (a = 0; a < nAt; a++) {
      if (a) {
        if (*sptr != *(sptr - 1)) {
          end_flag = true;
        } else if (*ss != ss_t::SHEET) {
          end_flag = true;
        }
        if (a == (nAt - 1)) {
          end_flag = 1;
        }
      }
      if (end_flag) {
        if ((cur_car != cCartoon_loop) && (cur_car != cCartoon_tube)) {
          f = 1;
          for (c = 0; c < flat_cycles; c++) {
            for (b = first + f; b <= last - f; b++) { /* iterative averaging */
              zero3f(t0);
              for (e = -f; e <= f; e++) {
                add3f(pv + 3 * (b + e), t0, t0);
              }
              scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
            }
            for (b = first + f; b <= last - f; b++) {
              if (!((*(flag_tmp + b) & cAtomFlag_no_smooth))) {
                copy3f(tmp + b * 3, pv + b * 3);
              }
            }
            for (b = first + f; b <= last - f; b++) {
              zero3f(t0);
              for (e = -f; e <= f; e++) {
                add3f(pvo + 3 * (b + e), t0, t0);
              }
              scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
            }
            for (b = first + f; b <= last - f; b++) {
              copy3f(tmp + b * 3, pvo + b * 3);
              /*                  normalize3f(pvo+b*3); */
            }
            for (b = first + f; b <= last - f; b++) {
              subtract3f(pv + (b + 1) * 3, pv + (b - 1) * 3, tmp + b * 3);
              normalize3f(tmp + b * 3);
              remove_component3f(pvo + b * 3, tmp + b * 3, pvo + b * 3);
              normalize3f(pvo + b * 3);
            }
          }
        }
        first = -1;
        last = -1;
        end_flag = false;
      }
      if (*ss == ss_t::SHEET) {
        if (first < 0) first = a;
        cur_car = *cc;
        last = a;
      }
      ss++;
      v0 += 3;
      sptr++;
      cc++;
    }
  }
}

void RepCartoonFlattenSheetsRefineTips(PyMOLGlobals* G, ObjectMolecule* obj,
                                       CoordSet* cs, int nAt, int* seg,
                                       const ss_t* ss, float* tv) {
  int *sptr, a;
  float* v2;
  float refine_tips;
  float t0[3];
  refine_tips = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_refine_tips);
  sptr = seg + 1;
  ++ss;
  v2 = tv + 3; /* normal */
  for (a = 1; a < (nAt - 1); a++) {
    if ((*ss == ss_t::SHEET) && (*sptr == *(sptr + 1)) &&
        (*sptr == *(sptr - 1))) { /* sheet in same segment */
      if ((*ss == *(ss + 1)) && (*ss != *(ss - 1))) { /* start, bias forwards */
        scale3f(v2 + 3, refine_tips, t0);
        add3f(t0, v2, v2);
        normalize3f(v2);
      } else if ((*ss != *(ss + 1)) &&
                 (*ss == *(ss - 1))) { /* end, bias backwards */
        scale3f(v2 - 3, refine_tips, t0);
        add3f(t0, v2, v2);
        normalize3f(v2);
      }
    }
    v2 += 3;
    sptr++;
    ss++;
  }
}

void RepCartoonSmoothLoops(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                           nuc_acid_data* ndata, const int nAt, const int* seg,
                           float* pv, const ss_t* ss, float* pvo,
                           const float* /* tv */, float* tmp,
                           const int* flag_tmp) {
  int last, first, end_flag, a, b, c, e, f;
  float t0[3];
  int smooth_first, smooth_last, smooth_cycles;
  smooth_first = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_cartoon_smooth_first);
  smooth_last = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_smooth_last);
  smooth_cycles = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_cartoon_smooth_cycles);

  auto* sptr = seg;
  last = 0;
  first = -1;
  end_flag = false;
  if (nAt > 1) {
    for (a = 0; a < nAt; a++) {
      if (a) {
        if (*sptr != *(sptr - 1)) {
          end_flag = true;
        } else if (*ss != ss_t::NONE) {
          end_flag = true;
        }
        if (a == (nAt - 1)) end_flag = 1;
      }
      if (end_flag) {
        if (a)
          if (first > 0) /* 011130 WLD */
            if (*(seg + first) == *(seg + first - 1)) first--;

        if (last > 0)
          if (*sptr == *(sptr - 1))
            if (last < (nAt - 1)) last++;

        for (f = smooth_first; f <= smooth_last; f++) {
          for (c = 0; c < smooth_cycles; c++) {
            for (b = first + f; b <= last - f; b++) { /* iterative averaging */
              zero3f(t0);
              for (e = -f; e <= f; e++) {
                add3f(pv + 3 * (b + e), t0, t0);
              }
              scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
            }
            for (b = first + f; b <= last - f; b++) {
              if (!(*(flag_tmp + b) & cAtomFlag_no_smooth)) {
                copy3f(tmp + b * 3, pv + b * 3);
              }
            }
            for (b = first + f; b <= last - f; b++) {
              zero3f(t0);
              for (e = -f; e <= f; e++) {
                add3f(pvo + 3 * (b + e), t0, t0);
              }
              scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
            }
            for (b = first + f; b <= last - f; b++) {
              copy3f(tmp + b * 3, pvo + b * 3);
              normalize3f(pvo + b * 3);
            }
          }
        }
        first = -1;
        last = -1;
        end_flag = false;
      }
      if (*ss == ss_t::NONE) {
        if (first < 0) first = a;
        last = a;
      }
      ss++;
      sptr++;
    }
  }
}
