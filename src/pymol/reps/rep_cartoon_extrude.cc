
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

// Sampling the cartoon spline and extruding its profiles (tube, oval,
// rectangle, arrow and others).

#include "pymol/base/err.h"
#include "pymol/reps/rep_cartoon.h"
#include "pymol/reps/rep_cartoon_internal.h"

void CartoonGenerateSample(PyMOLGlobals* G, int sampling, int* n_p, float dev,
                           float* vo, float* v1, float* v2, int c1, int c2,
                           float alpha1, float alpha2, int atom_index1,
                           int atom_index2, float power_a, float power_b,
                           float** vc_p, float** valpha_p, unsigned int** vi_p,
                           float** v_p, float** vn_p) {
  int b, i0;
  float f0, f1, f2, f3, f4;
  float a0;
  const float* v0;
  unsigned int* vi = *vi_p;
  float* valpha = *valpha_p;
  float *vc = *vc_p, *v = *v_p, *vn = *vn_p;
  for (b = 0; b < sampling; b++) { /* needs optimization */
    if (*n_p == 0) {
      /* provide starting point on first point in segment only... */
      f0 = ((float)b) / sampling; /* fraction of completion */
      if (f0 <= 0.5) {
        v0 = ColorGet(G, c1);
        i0 = atom_index1;
        a0 = alpha1;
      } else {
        v0 = ColorGet(G, c2);
        i0 = atom_index2;
        a0 = alpha2;
      }
      f0 = smooth(f0,
                  power_a); /* bias sampling towards the center of the curve */
      /* store colors and alpha*/
      *(vc++) = *(v0++);
      *(vc++) = *(v0++);
      *(vc++) = *(v0++);
      *(valpha++) = a0;
      *(vi++) = i0;
      /* start of line/cylinder */
      f1 = 1.0F - f0;
      f2 = smooth(f0, power_b);
      f3 = smooth(f1, power_b);
      f4 = dev * f2 * f3; /* displacement magnitude */

      *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);
      *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);
      *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);
      vn += 9;
      copy3f(vo, vn - 6); /* starter... */
      (*n_p)++;
    }
    f0 = ((float)b + 1) / sampling;
    if (f0 <= 0.5) {
      v0 = ColorGet(G, c1);
      i0 = atom_index1;
      a0 = alpha1;
    } else {
      v0 = ColorGet(G, c2);
      i0 = atom_index2;
      a0 = alpha2;
    }
    f0 =
        smooth(f0, power_a); /* bias sampling towards the center of the curve */

    /* store colors and alpha*/
    *(vc++) = *(v0++);
    *(vc++) = *(v0++);
    *(vc++) = *(v0++);
    *(valpha++) = a0;
    *(vi++) = i0;

    /* end of line/cylinder */
    f1 = 1.0F - f0;
    f2 = smooth(f0, power_b);
    f3 = smooth(f1, power_b);
    f4 = dev * f2 * f3; /* displacement magnitude */

    *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);
    *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);
    *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);
    /*                remove_component3f(vo,v2,o0);
                      remove_component3f(vo+3,v2,o0+3); */
    vn += 3;
    *(vn++) = f1 * (vo[0] * f2) + f0 * (vo[3] * f3);
    *(vn++) = f1 * (vo[1] * f2) + f0 * (vo[4] * f3);
    *(vn++) = f1 * (vo[2] * f2) + f0 * (vo[5] * f3);
    vn += 3;

    if (b == sampling - 1) copy3f(vo + 3, vn - 6); /* starter... */
    (*n_p)++;
  }
  (*vc_p) = vc;
  (*valpha_p) = valpha;
  (*vi_p) = vi;
  (*v_p) = v;
  (*vn_p) = vn;
}

void CartoonGenerateRefine(int refine, int sampling, float* v, float* vn,
                           float* vo, float* sampling_tmp) {
  int b, c;
  float t0[3], t1[3];
  float *p0, *p1, *p2, *p3;
  float f0, f1, f2, f3;
  c = refine;
  cross_product3f(vn + 3 - (sampling * 9), vn + 3 - 9,
                  t0);  // why is this here?
  cross_product3f(vo, vo + 3, t0);
  if ((sampling > 1) && length3f(t0) > R_SMALL4) {
    normalize3f(t0);
    while (c--) {
      p0 = v - (sampling * 3) - 3;
      p1 = v - (sampling * 3);
      p2 = v - (sampling * 3) + 3;
      for (b = 0; b < (sampling - 1); b++) {
        f0 = dot_product3f(t0, p0);
        f1 = dot_product3f(t0, p1);
        f2 = dot_product3f(t0, p2);

        f3 = (f2 + f0) / 2.0F;
        scale3f(t0, f3 - f1, t1);
        p3 = sampling_tmp + b * 3;
        add3f(t1, p1, p3);

        p0 = p1;
        p1 = p2;
        p2 += 3;
      }
      p1 = v - (sampling * 3);
      for (b = 0; b < (sampling - 1); b++) {
        p3 = sampling_tmp + b * 3;
        copy3f(p3, p1);
        p1 += 3;
      }
    }
  }
}

/**
 * Natural cubic spline coefficients for every contiguous cartoon segment in
 * pv (parameter s in [0,1] per interval).  Cyclic segments are left invalid
 * so that they fall back to the classic sampler (their closing edge is not
 * part of the open spline).
 */
CartoonSplineCoefs CartoonBuildSplineCoefs(
    int nAt, const int* seg, const float* pv,
    const std::unordered_map<int, CyclicSegmentInfo>* cyclicSegs) {
  CartoonSplineCoefs out;
  if (nAt < 2) return out;
  out.coef.assign(nAt * 12, 0.0f);
  out.valid.assign(nAt, 0);
  int start = 0;
  while (start < nAt) {
    int end = start;
    while (end + 1 < nAt && seg[end + 1] == seg[start]) end++;
    int const n = end - start + 1;
    bool const cyclic = cyclicSegs && cyclicSegs->count(seg[start]);
    if (n >= 2 && !cyclic) {
      // second derivatives M (natural ends: M[0] = M[n-1] = 0) via the
      // tridiagonal system  M[i-1] + 4 M[i] + M[i+1] = 6 (y[i+1] - 2 y[i] +
      // y[i-1])
      std::vector<double> M(3 * n, 0.0);
      if (n > 2) {
        int const k = n - 2;
        std::vector<double> diag(k), rhs(k);
        for (int ax = 0; ax < 3; ax++) {
          for (int i = 0; i < k; i++) {
            const float* y = pv + 3 * (start + i + 1) + ax;
            diag[i] = 4.0;
            rhs[i] = 6.0 * ((double)y[3] - 2.0 * (double)y[0] + (double)y[-3]);
          }
          for (int i = 1; i < k; i++) {
            double const w = 1.0 / diag[i - 1];
            diag[i] -= w;
            rhs[i] -= w * rhs[i - 1];
          }
          M[3 * (k) + ax] = rhs[k - 1] / diag[k - 1];  // point index k == n-2
          for (int i = k - 2; i >= 0; i--) {
            M[3 * (i + 1) + ax] = (rhs[i] - M[3 * (i + 2) + ax]) / diag[i];
          }
        }
      }
      for (int i = 0; i < n - 1; i++) {
        int const a = start + i;
        float* c = out.coef.data() + 12 * a;
        for (int ax = 0; ax < 3; ax++) {
          double const y0 = pv[3 * a + ax], y1 = pv[3 * (a + 1) + ax];
          double const m0 = M[3 * i + ax], m1 = M[3 * (i + 1) + ax];
          c[4 * ax + 0] = (float)y0;
          c[4 * ax + 1] = (float)((y1 - y0) - (2.0 * m0 + m1) / 6.0);
          c[4 * ax + 2] = (float)(m0 / 2.0);
          c[4 * ax + 3] = (float)((m1 - m0) / 6.0);
        }
        out.valid[a] = 1;
      }
    }
    start = end + 1;
  }
  return out;
}

/// Rodrigues rotation of v about unit axis k by angle ang (in place).
static inline void CartoonRotateAbout(float* v, const float* k, float ang) {
  float const c = cosf(ang), s = sinf(ang);
  float kxv[3], out[3];
  cross_product3f(k, v, kxv);
  float const kd = dot_product3f(k, v) * (1.0F - c);
  for (int i = 0; i < 3; i++) out[i] = v[i] * c + kxv[i] * s + k[i] * kd;
  copy3f(out, v);
}

/// Sigmoid twist profile (ChimeraX smooth_twist), rescaled to hit 0 and 1
/// exactly.
static inline float CartoonTwistProfile(float f) {
  auto const sig = [](float x) {
    return 1.0F / (1.0F + expf(-8.0F * (x - 0.5F)));
  };
  return (sig(f) - sig(0.0F)) / (sig(1.0F) - sig(0.0F));
}

/**
 * cartoon_spline replacement for CartoonGenerateSample for one interval.
 * Output layout is identical: positions in v, orientation vector in the middle
 * triple of each 9-float block of vn (the tangent triple is filled later by
 * ExtrudeComputeTangents), plus colors, alpha and pick ids per point.
 */
void CartoonGenerateSampleSpline(PyMOLGlobals* G, int sampling, int* n_p,
                                 const float* coef, const float* vo, int c1,
                                 int c2, float alpha1, float alpha2,
                                 int atom_index1, int atom_index2, float** vc_p,
                                 float** valpha_p, unsigned int** vi_p,
                                 float** v_p, float** vn_p) {
  int const np = sampling + 1;
  std::vector<float> pos(3 * np), tan(3 * np), nrm(3 * np);

  // positions and unit tangents on the cubic
  for (int b = 0; b < np; b++) {
    double const s = (double)b / sampling;
    float* p = pos.data() + 3 * b;
    float* t = tan.data() + 3 * b;
    for (int ax = 0; ax < 3; ax++) {
      const float* c = coef + 4 * ax;
      p[ax] = (float)(c[0] + s * (c[1] + s * (c[2] + s * c[3])));
      t[ax] = (float)(c[1] + s * (2.0 * c[2] + 3.0 * s * c[3]));
    }
    if (length3f(t) > R_SMALL8) {
      normalize3f(t);
    } else if (b > 0) {
      copy3f(t - 3, t);
    } else {
      // degenerate start tangent: use the chord
      subtract3f(pos.data() + 3 * (np - 1), p, t);
      normalize3f(t);
    }
  }

  // parallel transport of the first orientation vector along the interval
  float n[3];
  remove_component3f(vo, tan.data(), n);
  if (length3f(n) > R_SMALL8) {
    normalize3f(n);
  } else {
    float dummy[3];
    get_system1f3f(tan.data(), n, dummy);
  }
  copy3f(n, nrm.data());
  for (int b = 1; b < np; b++) {
    const float* t0 = tan.data() + 3 * (b - 1);
    const float* t1 = tan.data() + 3 * b;
    float axis[3];
    cross_product3f(t0, t1, axis);
    float const al = length3f(axis);
    if (al > R_SMALL8) {
      scale3f(axis, 1.0F / al, axis);
      CartoonRotateAbout(n, axis, atan2f(al, dot_product3f(t0, t1)));
    }
    remove_component3f(n, t1, n);
    if (length3f(n) > R_SMALL8) normalize3f(n);
    copy3f(n, nrm.data() + 3 * b);
  }

  // twist so the end lands on the next residue's orientation vector
  const float* tEnd = tan.data() + 3 * sampling;
  float tgt[3];
  remove_component3f(vo + 3, tEnd, tgt);
  if (length3f(tgt) > R_SMALL8) {
    normalize3f(tgt);
    float cx[3];
    cross_product3f(nrm.data() + 3 * sampling, tgt, cx);
    float const ang = atan2f(dot_product3f(tEnd, cx),
                             dot_product3f(nrm.data() + 3 * sampling, tgt));
    for (int b = 1; b < sampling; b++) {
      CartoonRotateAbout(nrm.data() + 3 * b, tan.data() + 3 * b,
                         ang * CartoonTwistProfile((float)b / sampling));
    }
    copy3f(tgt, nrm.data() + 3 * sampling);
  }

  // emit points, mirroring CartoonGenerateSample's bookkeeping
  unsigned int* vi = *vi_p;
  float* valpha = *valpha_p;
  float *vc = *vc_p, *v = *v_p, *vn = *vn_p;
  auto emit = [&](int b, int color, int idx, float a0) {
    const float* col = ColorGet(G, color);
    *(vc++) = col[0];
    *(vc++) = col[1];
    *(vc++) = col[2];
    *(valpha++) = a0;
    *(vi++) = idx;
    copy3f(pos.data() + 3 * b, v);
    v += 3;
    vn += 3;
    copy3f(nrm.data() + 3 * b, vn);
    vn += 6;
    (*n_p)++;
  };
  for (int b = 0; b < sampling; b++) {
    if (*n_p == 0) {
      // first point of the extrusion only; later intervals share it
      emit(0, c1, atom_index1, alpha1);
    }
    float const f0 = ((float)b + 1) / sampling;
    if (f0 <= 0.5F) {
      emit(b + 1, c1, atom_index1, alpha1);
    } else {
      emit(b + 1, c2, atom_index2, alpha2);
    }
  }
  *vc_p = vc;
  *valpha_p = valpha;
  *vi_p = vi;
  *v_p = v;
  *vn_p = vn;
}

int CartoonExtrudeTube(short use_cylinders_for_strands, CExtrude* ex, CGO* cgo,
                       float tube_radius, int tube_quality, cCylCap tube_cap) {
  int ok = true;
  if (use_cylinders_for_strands) {
    ok &= ExtrudeCylindersToCGO(ex, cgo, tube_radius);
  } else {
    ok &= ExtrudeCircle(ex, tube_quality, tube_radius);
    if (ok) ExtrudeBuildNormals1f(ex);
    if (ok)
      ok &= ExtrudeCGOSurfaceTube(ex, cgo, tube_cap, nullptr,
                                  use_cylinders_for_strands);
  }
  return ok;
}

int CartoonExtrudePutty(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                        CGO* cgo, CExtrude* ex, int putty_quality,
                        float putty_radius, float* putty_vals, int sampling) {
  int ok;
  ok = ExtrudeCircle(ex, putty_quality, putty_radius);
  if (ok) ExtrudeBuildNormals1f(ex);
  if (ok)
    ok &= ExtrudeComputePuttyScaleFactors(
        ex, obj,
        SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_cartoon_putty_transform),
        putty_vals[0], putty_vals[1], putty_vals[2], putty_vals[3],
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_cartoon_putty_scale_power),
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_cartoon_putty_range),
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_cartoon_putty_scale_min),
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                     cSetting_cartoon_putty_scale_max),
        sampling / 2);
  if (ok) ok &= ExtrudeCGOSurfaceVariableTube(ex, cgo, cCylCap::Flat);
  return ok;
}

int CartoonExtrudeCircle(CExtrude* ex, CGO* cgo,
                         short use_cylinders_for_strands, int loop_quality,
                         float loop_radius, cCylCap loop_cap, int dash) {
  int ok;
  ok = ExtrudeCircle(ex, loop_quality, loop_radius);
  if (ok) ExtrudeBuildNormals1f(ex);
  if (ok)
    ok &= ExtrudeCGOSurfaceTube(ex, cgo, loop_cap, nullptr,
                                use_cylinders_for_strands, dash);
  return ok;
}

/**
 * Modify `ex` to approximate a curved helix axis and render it as a tube.
 *
 * @param[in,out] ex Helix trace which will be converted to helix center trace
 * @param[in,out] cgo CGO to render to
 * @param quality Number of circular sampling points
 * @param radius Cylindrical helix radius
 * @param sampling Samples per residue
 */
void CartoonExtrudeCurvedCylindricalHelix(CExtrude* ex, CGO* cgo, int quality,
                                          float radius, int sampling) {
  ExtrudeShiftToAxis(ex, radius, sampling);
  ExtrudeCircle(ex, quality, radius);
  ExtrudeCGOSurfaceTube(ex, cgo, cCylCap::Flat, nullptr, false);
}

int CartoonExtrudeRect(PyMOLGlobals* G, CExtrude* ex, CGO* cgo, float width,
                       float length, int highlight_color) {
  int ok;
  if (highlight_color < 0) {
    ok = ExtrudeRectangle(ex, width, length, 0);
    if (ok) ExtrudeBuildNormals2f(ex);
    if (ok) ok &= ExtrudeCGOSurfacePolygon(ex, cgo, cCylCap::Flat, nullptr);
  } else {
    ok = ExtrudeRectangle(ex, width, length, 1);
    if (ok) ExtrudeBuildNormals2f(ex);
    if (ok) ok &= ExtrudeCGOSurfacePolygon(ex, cgo, cCylCap::None, nullptr);
    if (ok) {
      ok &= ExtrudeRectangle(ex, width, length, 2);
      if (ok) ExtrudeBuildNormals2f(ex);
      if (ok)
        ok &= ExtrudeCGOSurfacePolygon(ex, cgo, cCylCap::Flat,
                                       ColorGet(G, highlight_color));
    }
  }
  return ok;
}

int CartoonExtrudeOval(PyMOLGlobals* G, CExtrude* ex, CGO* cgo,
                       short use_cylinders_for_strands, int oval_quality,
                       float oval_width, float oval_length,
                       int highlight_color) {
  int ok;
  ok = ExtrudeOval(ex, oval_quality, oval_width, oval_length);
  if (ok) ExtrudeBuildNormals2f(ex);
  if (ok) {
    if (highlight_color < 0)
      ok &= ExtrudeCGOSurfaceTube(ex, cgo, cCylCap::Flat, nullptr,
                                  use_cylinders_for_strands);
    else
      ok &= ExtrudeCGOSurfaceTube(ex, cgo, cCylCap::Flat,
                                  ColorGet(G, highlight_color),
                                  use_cylinders_for_strands);
  }
  return ok;
}

int CartoonExtrudeArrow(PyMOLGlobals* G, CExtrude* ex, CGO* cgo, int sampling,
                        float width, float length, int highlight_color) {
  int ok;
  ok = ExtrudeRectangle(ex, width, length,
                        0);  // ex->Ns is always 8 for ExtrudeCGOSurfaceStrand
  if (ok) ExtrudeBuildNormals2f(ex);
  if (ok) {
    if (highlight_color < 0)
      ok &= ExtrudeCGOSurfaceStrand(ex, cgo, sampling, nullptr);
    else
      ok &= ExtrudeCGOSurfaceStrand(ex, cgo, sampling,
                                    ColorGet(G, highlight_color));
  }
  /* for PLY files
     ExtrudeCircle(ex,loop_quality,loop_radius);
     ExtrudeBuildNormals1f(ex);
     ExtrudeCGOSurfaceTube(ex,cgo,loop_cap,NULL);
  */
  return ok;
}

int CartoonExtrudeDumbbell(PyMOLGlobals* G, CExtrude* ex, CGO* cgo,
                           int sampling, float dumbbell_width,
                           float dumbbell_length, int highlight_color,
                           int loop_quality, float dumbbell_radius,
                           short use_cylinders_for_strands) {
  int ok;
  CExtrude* ex1 = nullptr;
  if (highlight_color < 0) {
    ok = ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 0);
    if (ok) ExtrudeBuildNormals2f(ex);
    if (ok) ok &= ExtrudeCGOSurfacePolygonTaper(ex, cgo, sampling, nullptr);
  } else {
    ok = ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 1);
    if (ok) ExtrudeBuildNormals2f(ex);
    if (ok) ok &= ExtrudeCGOSurfacePolygonTaper(ex, cgo, sampling, nullptr);
    if (ok) ok &= ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 2);
    if (ok) ExtrudeBuildNormals2f(ex);
    if (ok)
      ok &= ExtrudeCGOSurfacePolygonTaper(ex, cgo, sampling,
                                          ColorGet(G, highlight_color));
  }
  /*
    ExtrudeCGOSurfacePolygonX(ex,cgo,1); */

  if (ok) {
    ex1 = ExtrudeCopyPointsNormalsColors(ex);
    CHECKOK(ok, ex1);
    if (ok) ExtrudeDumbbellEdge(ex1, sampling, -1, dumbbell_length);
    if (ok) ok &= ExtrudeComputeTangents(ex1);
  }
  if (ok) ok &= ExtrudeCircle(ex1, loop_quality, dumbbell_radius);
  if (ok) ExtrudeBuildNormals1f(ex1);

  if (ok)
    ok &= ExtrudeCGOSurfaceTube(ex1, cgo, cCylCap::Flat, nullptr,
                                use_cylinders_for_strands);
  if (ok) {
    ExtrudeFree(ex1);
    ex1 = ExtrudeCopyPointsNormalsColors(ex);
    CHECKOK(ok, ex1);
    if (ok) ExtrudeDumbbellEdge(ex1, sampling, 1, dumbbell_length);
    if (ok) ok &= ExtrudeComputeTangents(ex1);
    if (ok) ok &= ExtrudeCircle(ex1, loop_quality, dumbbell_radius);
    if (ok) ExtrudeBuildNormals1f(ex1);
    if (ok)
      ok &= ExtrudeCGOSurfaceTube(ex1, cgo, cCylCap::Flat, nullptr,
                                  use_cylinders_for_strands);
  }
  if (ex1) ExtrudeFree(ex1);
  return ok;
}
