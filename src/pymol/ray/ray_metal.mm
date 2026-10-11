/*
 * RayMetal — see RayMetal.h. The kernel below is a line-by-line port of
 * BasisHitOrthoscopic, ZLineToSphere, ZLineToSphereCapped, ZLineClipPoint
 * (src/pymol/ray/basis.cc) and the per-pixel body of RayTraceThread
 * (src/pymol/ray/ray.cc), keeping their order of operations so the GPU image
 * matches the CPU one to within float rounding.
 */

#include "pymol/base/busy.h"
#include "pymol/ray/ray_metal.h"
#include "pymol/ui/ortho.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <atomic>
#include <memory>
#include <chrono>
#include <cstring>
#include <mutex>

#if __has_feature(objc_arc)
#define RM_RELEASE(x) ((void)0)
#else
#define RM_RELEASE(x) [(x) release]
#endif

static const char* kRayMetalSource = R"MSL(
#include <metal_stdlib>
using namespace metal;

struct Prim {
  float c1[3], c2[3], c3[3], tr[3];
  float r1, l1, trans;
  int vert, type, cap1, cap2, cull, no_lighting;
};

struct Params {
  int width, height, x_start, x_stop, y_start, y_stop;
  float invWdthRange, invHgtRange, vol0, vol2, border_offset;
  float front, back, fudge0, fudge1;
  int n_vert, n_eElem;
  int dim0, dim1, dim2, d1d2;
  int iMin0, iMin1, iMin2, iMax0, iMax1, iMax2;
  float mapMin0, mapMin1, mapMin2, recipDiv;
  float ambient, lreflect, direct, direct_shade, legacy, legacy_1m, ray_scatter;
  float settingPower, settingReflectPower, settingSpecPower, settingSpecReflect,
      settingSpecDirect, settingSpecDirectPower;
  int spec_count, n_basis, two_sided_lighting, has_bp2;
  float lightN[30], specN[30];
  float project_triangle;
  int fogFlag, fogRangeFlag;
  float fog, fog_start, inv1minusFogStart, invFrontMinusBack;
  int orig_opaque_back;
  uint fore_mask;
  float bkrd0, bkrd1, bkrd2;
  int trans_cont_flag, backface_cull, trans_mode, max_pass;
  float trans_cont, inv_trans_cont, trans_oblique, oblique_power;
  float ray_trans_spec, trans_spec_cut, trans_spec_scale;
  float trans_cutoff, persist_cutoff;
  int blend_colors;
  float red_blend, green_blend, blue_blend;
  int write_depth;
  uint background;
  int ox0, oy0, ox1, oy1;
  int collect;
};

#define cPrimSphere 1
#define cPrimCylinder 2
#define cPrimTriangle 3
#define cPrimSausage 4
#define CapNone 0
#define CapFlat 1
#define CapRound 2
#define MapBorder 2
constant float kR_SMALL4 = 0.0001f;
constant float FLTMAX = 3.402823466e+38f;

struct Bufs {
  device const Prim* prim;
  device const float* V;     // Vertex
  device const float* N;     // Normal
  device const int* V2N;
  device const int* V2P;
  device const int* EHead;
  device const int* EList;
  device const int* EMask;
  device const float* VMax;  // per voxel list: the farthest dist it can hit
  device const uint4* T;     // test records, read as integers (their int
                             // fields are denormal as floats and would flush)
};

static inline float sqrt1f(float v) { return (v > 0) ? precise::sqrt(v) : 0; }
/* correctly rounded, as on the CPU (plain / may be off by an ulp) */
static inline float pdiv(float a, float b) { return precise::divide(a, b); }
static inline float ppow(float a, float b) { return precise::pow(a, b); }

static inline void normalize3(thread float* v) {
  float vlen = sqrt1f((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
  if (vlen > 1e-8f) {
    float s = 1 / vlen;
    v[0] = v[0] * s; v[1] = v[1] * s; v[2] = v[2] * s;
  } else {
    v[0] = v[1] = v[2] = 0;
  }
}

static inline float dot3(thread const float* a, device const float* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline float dot3c(thread const float* a, constant float* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline float dot3cd(constant float* a, device const float* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* ZLineClipPoint: base[2] is 0 for orthoscopic rays */
static inline float ZLineClipPoint(float b0, float b1, thread const float* point,
    thread float& alongNormalSq, float cutoff) {
  float result = FLTMAX;
  float hyp0 = point[0] - b0;
  if (fabs(hyp0) > cutoff) return result;
  float hyp1 = point[1] - b1;
  if (fabs(hyp1) > cutoff) return result;
  float hyp2 = point[2] - 0.0f;
  if (hyp2 < 0.0f) {
    alongNormalSq = (hyp2 * hyp2);
    result = (hyp0 * hyp0) + (hyp1 * hyp1);
  }
  return result;
}

static int ZLineToSphere(float b0, float b1, thread const float* point,
    thread const float* dir, float radius, float maxial, thread float* sphere,
    thread float& asum, thread const float* pre) {
  const float _0 = 0.0f;
  float point0 = point[0], point1 = point[1], point2 = point[2];
  float perpAxis0 = pre[0], perpAxis1 = pre[1];
  float intra0 = point0 - b0;
  float intra1 = point1 - b1;
  const float dir0 = dir[0], dir1 = dir[1], dir2 = dir[2];
  float perpDist = intra0 * perpAxis0 + intra1 * perpAxis1;
  if ((perpDist < _0 ? -perpDist : perpDist) > radius) return 0;
  float dangle = -dir2;
  float ab_dangle = fabs(dangle);
  if (ab_dangle > (1.0f - kR_SMALL4)) {
    if (dangle > _0) {
      sphere[0] = point0; sphere[1] = point1; sphere[2] = point2;
    } else {
      sphere[0] = dir0 * maxial + point0;
      sphere[1] = dir1 * maxial + point1;
      sphere[2] = dir2 * maxial + point2;
    }
    return 1;
  }
  float tan_acos_dangle;
  if (ab_dangle > kR_SMALL4)
    tan_acos_dangle = sqrt1f(1.0f - dangle * dangle) / dangle;
  else
    tan_acos_dangle = FLTMAX;
  float intra_p2 = point2 - 0.0f;
  float dot = intra0 * perpAxis0 + intra1 * perpAxis1;
  float intra_p0 = intra0 - perpAxis0 * dot;
  float intra_p1 = intra1 - perpAxis1 * dot;
  dot = intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2;
  float vradial0 = intra_p0 - dir0 * dot;
  float vradial1 = intra_p1 - dir1 * dot;
  float vradial2 = intra_p2 - dir2 * dot;
  float radialsq = ((vradial0 * vradial0) + (vradial1 * vradial1) + (vradial2 * vradial2));
  float axial_perp;
  if (ab_dangle < kR_SMALL4)
    axial_perp = _0;
  else
    axial_perp = sqrt1f(radialsq) / tan_acos_dangle;
  float axial = sqrt1f(((intra_p0 * intra_p0) + (intra_p1 * intra_p1) + (intra_p2 * intra_p2)) - radialsq);
  if ((intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2) >= _0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;
  float radial = radius * radius - perpDist * perpDist;
  radial = sqrt1f(radial);
  float axial_sum;
  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;
  if (axial_sum < _0)
    axial_sum = _0;
  else if (axial_sum > maxial)
    axial_sum = maxial;
  sphere[0] = dir0 * axial_sum + point0;
  sphere[1] = dir1 * axial_sum + point1;
  sphere[2] = dir2 * axial_sum + point2;
  asum = axial_sum;
  return 1;
}

static int ZLineToSphereCapped(float b0, float b1, thread const float* point,
    thread const float* dir, float radius, float maxial, thread float* sphere,
    thread float& asum, int cap1, int cap2, thread const float* pre) {
  float base[3] = {b0, b1, 0.0f};
  float perpAxis[3], intra_p[3], intra[3], vradial[3], diff[3], fpoint[3], proj[3];
  perpAxis[0] = pre[0];
  perpAxis[1] = pre[1];
  intra[0] = point[0] - base[0];
  intra[1] = point[1] - base[1];
  float perpDist = intra[0] * perpAxis[0] + intra[1] * perpAxis[1];
  if (fabs(perpDist) > radius) return 0;
  perpAxis[2] = 0.0f;
  intra[2] = point[2] - base[2];
  float dangle = -dir[2];
  float ab_dangle = fabs(dangle);
  if (ab_dangle > (1 - kR_SMALL4)) {
    vradial[0] = point[0] - base[0];
    vradial[1] = point[1] - base[1];
    vradial[2] = 0.0f;
    float radial = sqrt1f(vradial[0] * vradial[0] + vradial[1] * vradial[1] + vradial[2] * vradial[2]);
    if (radial > radius) return 0;
    if (dangle > 0.0f) {
      switch (cap1) {
      case CapFlat:
        sphere[0] = base[0]; sphere[1] = base[1]; sphere[2] = point[2] - radius;
        break;
      case CapRound:
        sphere[0] = point[0]; sphere[1] = point[1]; sphere[2] = point[2];
      }
      return 1;
    } else {
      switch (cap1) {
      case CapFlat:
        sphere[0] = base[0]; sphere[1] = base[1];
        sphere[2] = dir[2] * maxial + point[2] - radius;
        break;
      case CapRound:
        sphere[0] = dir[0] * maxial + point[0];
        sphere[1] = dir[1] * maxial + point[1];
        sphere[2] = dir[2] * maxial + point[2];
      }
      return 1;
    }
  }
  float tan_acos_dangle = sqrt1f(1 - dangle * dangle) / dangle;
  {
    float d = intra[0] * perpAxis[0] + intra[1] * perpAxis[1] + intra[2] * perpAxis[2];
    intra_p[0] = intra[0] - perpAxis[0] * d;
    intra_p[1] = intra[1] - perpAxis[1] * d;
    intra_p[2] = intra[2] - perpAxis[2] * d;
  }
  {
    float d = intra_p[0] * dir[0] + intra_p[1] * dir[1] + intra_p[2] * dir[2];
    vradial[0] = intra_p[0] - dir[0] * d;
    vradial[1] = intra_p[1] - dir[1] * d;
    vradial[2] = intra_p[2] - dir[2] * d;
  }
  float radialsq = (vradial[0] * vradial[0]) + (vradial[1] * vradial[1]) + (vradial[2] * vradial[2]);
  float axial_perp;
  if (ab_dangle < kR_SMALL4)
    axial_perp = 0;
  else
    axial_perp = sqrt1f(radialsq) / tan_acos_dangle;
  float axial = ((intra_p[0] * intra_p[0]) + (intra_p[1] * intra_p[1]) + (intra_p[2] * intra_p[2])) - radialsq;
  axial = sqrt1f(axial);
  if ((intra_p[0] * dir[0] + intra_p[1] * dir[1] + intra_p[2] * dir[2]) >= 0.0f)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;
  float radial = radius * radius - perpDist * perpDist;
  radial = sqrt1f(radial);
  float axial_sum;
  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;
  if (axial_sum < 0) {
    switch (cap1) {
    case CapFlat: {
      diff[0] = point[0] - base[0]; diff[1] = point[1] - base[1]; diff[2] = point[2] - base[2];
      float d = diff[0] * dir[0] + diff[1] * dir[1] + diff[2] * dir[2];
      proj[0] = dir[0] * d; proj[1] = dir[1] * d; proj[2] = dir[2] * d;
      float len_proj = sqrt1f(proj[0] * proj[0] + proj[1] * proj[1] + proj[2] * proj[2]);
      dangle = -proj[2] / len_proj;
      if (fabs(dangle) < kR_SMALL4) return 0;
      sphere[0] = base[0]; sphere[1] = base[1];
      sphere[2] = base[2] - len_proj / dangle;
      {
        float dx = sphere[0] - point[0], dy = sphere[1] - point[1], dz = sphere[2] - point[2];
        if (sqrt1f(dz * dz + dy * dy + dx * dx) > radius) return 0;
      }
      sphere[0] += dir[0] * radius; sphere[1] += dir[1] * radius; sphere[2] += dir[2] * radius;
      asum = 0;
      break;
    }
    case CapRound:
      axial_sum = 0;
      sphere[0] = point[0]; sphere[1] = point[1]; sphere[2] = point[2];
      asum = axial_sum;
      break;
    default:
      return 0;
    }
  } else if (axial_sum > maxial) {
    switch (cap2) {
    case CapFlat: {
      fpoint[0] = dir[0] * maxial; fpoint[1] = dir[1] * maxial; fpoint[2] = dir[2] * maxial;
      fpoint[0] = fpoint[0] + point[0]; fpoint[1] = fpoint[1] + point[1]; fpoint[2] = fpoint[2] + point[2];
      diff[0] = fpoint[0] - base[0]; diff[1] = fpoint[1] - base[1]; diff[2] = fpoint[2] - base[2];
      float d = diff[0] * dir[0] + diff[1] * dir[1] + diff[2] * dir[2];
      proj[0] = dir[0] * d; proj[1] = dir[1] * d; proj[2] = dir[2] * d;
      float len_proj = sqrt1f(proj[0] * proj[0] + proj[1] * proj[1] + proj[2] * proj[2]);
      dangle = -proj[2] / len_proj;
      if (fabs(dangle) < kR_SMALL4) return 0;
      sphere[0] = base[0]; sphere[1] = base[1];
      sphere[2] = base[2] - len_proj / dangle;
      {
        float dx = sphere[0] - fpoint[0], dy = sphere[1] - fpoint[1], dz = sphere[2] - fpoint[2];
        if (sqrt1f(dz * dz + dy * dy + dx * dx) > radius) return 0;
      }
      sphere[0] -= dir[0] * radius; sphere[1] -= dir[1] * radius; sphere[2] -= dir[2] * radius;
      asum = maxial;
      break;
    }
    case CapRound:
      axial_sum = maxial;
      sphere[0] = dir[0] * axial_sum + point[0];
      sphere[1] = dir[1] * axial_sum + point[1];
      sphere[2] = dir[2] * axial_sum + point[2];
      asum = axial_sum;
      break;
    default:
      return 0;
    }
  } else {
    sphere[0] = dir[0] * axial_sum + point[0];
    sphere[1] = dir[1] * axial_sum + point[1];
    sphere[2] = dir[2] * axial_sum + point[2];
    asum = axial_sum;
  }
  return 1;
}

struct HitResult {
  int prim;
  float tri1, tri2, dist;
  float sphere[3];
};

/* One list entry's intersection test, as in BasisHitOrthoscopic. Each entry
   is a 64-byte record (uint4 x4) of everything its test reads,
     [0] point xyz (triangle: its first vertex), prim | type << 27 | cull << 31
     [1] radius, radius^2, l1, trans   (triangle: Precomp[6], Precomp[7], -, trans)
     [2] cylinder axis xyz, caps       (triangle: Precomp[0], [1], [3], -)
     [3] Precomp[0], [1], -, vert      (triangle: Precomp[4], [2], [5], vert)
   True for a hit at `dist`; tri1, tri2 and sph are PyMOL's per-ray scratch,
   carried from test to test exactly as there. */
struct Entry {
  int v2p, type, vert;
  float trans, l1;
};

static bool TestEntry(constant Params& P, device const uint4* rec, float b0, float b1,
    thread Entry& E, thread float& dist, thread float& tri1, thread float& tri2,
    thread float* sph) {
  uint4 ua = rec[0], ub = rec[1], uc = rec[2], ud = rec[3];
  E.v2p = (int) (ua.w & 0x07FFFFFFu);
  E.type = (int) ((ua.w >> 27) & 0xF);
  E.vert = (int) ud.w;
  E.trans = as_type<float>(ub.w);
  E.l1 = as_type<float>(ub.z);
  float R = as_type<float>(ub.x), R2 = as_type<float>(ub.y), oppSq;
  float pt[3] = {as_type<float>(ua.x), as_type<float>(ua.y), as_type<float>(ua.z)};
  float axis[3] = {as_type<float>(uc.x), as_type<float>(uc.y), as_type<float>(uc.z)};
  float pre[2] = {as_type<float>(ud.x), as_type<float>(ud.y)};
  switch (E.type) {
  case cPrimTriangle:
    if (!(ua.w >> 31) && R != 0) {     /* R = Precomp[6], R2 = Precomp[7] */
      float tvec0 = b0 - pt[0];
      float tvec1 = b1 - pt[1];
      tri1 = (tvec0 * pre[0] - tvec1 * axis[2]) * R2;
      tri2 = -(tvec0 * axis[1] - tvec1 * axis[0]) * R2;
      if (!((tri1 < P.fudge0) || (tri2 < P.fudge0) ||
            (tri1 > P.fudge1) || ((tri1 + tri2) > P.fudge1))) {
        dist = (0.0f - (tri1 * pre[1]) - (tri2 * as_type<float>(ud.z)) - pt[2]);
        return true;
      }
    }
    return false;
  case cPrimSphere:
    oppSq = ZLineClipPoint(b0, b1, pt, dist, R);
    if (oppSq <= R2) {
      dist = (sqrt1f(dist) - sqrt1f((R2 - oppSq)));
      return true;
    }
    return false;
  case cPrimCylinder:
  case cPrimSausage:
    if (E.type == cPrimCylinder
            ? ZLineToSphereCapped(b0, b1, pt, axis, R, E.l1, sph, tri1,
                  (int) (uc.w & 0xF), (int) ((uc.w >> 4) & 0xF), pre)
            : ZLineToSphere(b0, b1, pt, axis, R, E.l1, sph, tri1, pre)) {
      oppSq = ZLineClipPoint(b0, b1, sph, dist, R);
      if (oppSq <= R2) {
        dist = (sqrt1f(dist) - sqrt1f((R2 - oppSq)));
        return true;
      }
    }
    return false;
  }
  return false;
}

/* PyMOL's acceptance of a hit (front, back, transparency, sausage carry-through) */
static inline bool Accept(thread const Entry& E, float dist, float r_dist, float front,
    float back, float excl_trans) {
  return (dist < r_dist) && (dist >= front) && (dist <= back) && (E.trans != 1.0f) &&
         !(E.type == cPrimSausage && excl_trans != 0.0f && (E.trans > 0.0f) && (dist < excl_trans));
}

/* the ray's grid column and its first voxel (MapInsideXY); false: no column */
static bool Column(constant Params& P, thread const Bufs& B, float b0, float b1,
    thread int& xx, thread int& c) {
  float iDiv = P.recipDiv;
  int a = (int) ((b0 - P.mapMin0) * iDiv) + MapBorder;
  int b = (int) ((b1 - P.mapMin1) * iDiv) + MapBorder;
  c = (int) ((0.0f - P.mapMin2) * iDiv) + MapBorder + 1;
  if (a < P.iMin0) { if ((P.iMin0 - a) > 1) return false; else a = P.iMin0; }
  else if (a > P.iMax0) { if ((a - P.iMax0) > 1) return false; else a = P.iMax0; }
  if (b < P.iMin1) { if ((P.iMin1 - b) > 1) return false; else b = P.iMin1; }
  else if (b > P.iMax1) { if ((b - P.iMax1) > 1) return false; else b = P.iMax1; }
  if (!B.EMask[P.dim1 * a + b]) return false;
  c = clamp(c, P.iMin2, P.iMax2);
  xx = (a * P.d1d2) + (b * P.dim2) + c;
  return true;
}

/* the voxel holding the point at `dist` along the ray (MapLocus) */
static inline int VoxelAt(constant Params& P, float dist) {
  return clamp((int) ((0.0f - dist - P.mapMin2) * P.recipDiv) + MapBorder, P.iMin2, P.iMax2);
}

static void SetHit(thread const Bufs& B, int minIndex, float r_dist, float r_tri1, float r_tri2,
    thread const float* r_sph, thread HitResult& H) {
  H.prim = B.V2P[minIndex];
  bool sphere = B.prim[H.prim].type == cPrimSphere;
  for (int k = 0; k < 3; k++) H.sphere[k] = sphere ? B.V[minIndex * 3 + k] : r_sph[k];
  H.tri1 = r_tri1;
  H.tri2 = r_tri2;
  H.dist = r_dist;
}

/* BasisHitOrthoscopic over the column-deduplicated lists (see the host) */
static int BasisHitOrthoscopic(constant Params& P, thread const Bufs& B,
    float b0, float b1, float front, float back, float excl_trans,
    int except1, int except2, thread HitResult& H) {
  int xx, c;
  if (!Column(P, B, b0, b1, xx, c)) return -1;
  int minIndex = -1;
  float dist = 0, tri1 = 0, tri2 = 0, sph[3] = {0, 0, 0};
  float r_tri1 = 0, r_tri2 = 0, r_dist = FLTMAX, r_sph[3] = {0, 0, 0};
  if (except1 >= 0) except1 = B.V2P[except1];
  if (except2 >= 0) except2 = B.V2P[except2];
  for (; c >= MapBorder; c--, xx--) {
    int h = B.EHead[xx];
    /* a list whose primitives all lie in front of `front` can't be hit */
    if ((h > 0) && (h < P.n_eElem) && B.VMax[xx] >= front) {
      for (int ip = h, i = B.EList[ip]; (i >= 0) && (i < P.n_vert); i = B.EList[++ip]) {
        uint v2p = B.T[i * 4].w & 0x07FFFFFFu;
        if ((int) v2p == except1 || (int) v2p == except2) continue;
        Entry E;
        if (TestEntry(P, B.T + i * 4, b0, b1, E, dist, tri1, tri2, sph) &&
            Accept(E, dist, r_dist, front, back, excl_trans)) {
          minIndex = E.vert;
          r_dist = dist;
          if (E.type == cPrimTriangle) { r_tri1 = tri1; r_tri2 = tri2; }
          else if (E.type != cPrimSphere) {
            if (E.l1 > kR_SMALL4) r_tri1 = tri1 / E.l1;
            for (int k = 0; k < 3; k++) r_sph[k] = sph[k];
          }
        }
      }
    }
    /* stop once the hit lies in front of the next voxel */
    if (minIndex > -1 && VoxelAt(P, r_dist) > c) break;
  }
  if (minIndex > -1) SetHit(B, minIndex, r_dist, r_tri1, r_tri2, r_sph, H);
  return minIndex;
}

/* Transparency passes. Each pass of PyMOL walks the column again for the
   nearest hit behind the last one. A test's outcome doesn't depend on the
   pass, so one walk can record every hit behind `front` (in test order, with
   its voxel) and each pass picks from those, applying the pass's front,
   exclusions and PyMOL's per-voxel stop rule. (Not used when a test can
   read another's leftovers — end-on or zero-length cylinders — see host.) */
#define MAX_HITS 64
struct Hit {
  float dist, tri1, tri2, sph[3];
  int vert, v2p, type, voxel;
  float trans;
};

static int CollectHits(constant Params& P, thread const Bufs& B, float b0, float b1,
    float front, float back, thread Hit* hits) {
  int xx, c, n = 0;
  if (!Column(P, B, b0, b1, xx, c)) return 0;
  float dist = 0, tri1 = 0, tri2 = 0, sph[3] = {0, 0, 0};
  for (; c >= MapBorder; c--, xx--) {
    int h = B.EHead[xx];
    if (!((h > 0) && (h < P.n_eElem) && B.VMax[xx] >= front)) continue;
    for (int ip = h, i = B.EList[ip]; (i >= 0) && (i < P.n_vert); i = B.EList[++ip]) {
      Entry E;
      if (TestEntry(P, B.T + i * 4, b0, b1, E, dist, tri1, tri2, sph) &&
          dist >= front && dist <= back && E.trans != 1.0f) {
        if (n == MAX_HITS) return -2;
        thread Hit& t = hits[n++];
        t.dist = dist;
        t.tri1 = (E.type == cPrimTriangle) ? tri1 : tri1 / E.l1;
        t.tri2 = tri2;
        for (int k = 0; k < 3; k++) t.sph[k] = sph[k];
        t.vert = E.vert; t.v2p = E.v2p; t.type = E.type; t.voxel = c; t.trans = E.trans;
      }
    }
  }
  return n;
}

static int PickHit(constant Params& P, thread const Bufs& B, thread const Hit* hits, int n,
    float front, float back, float excl_trans, int except1, int except2, thread HitResult& H) {
  if (except1 >= 0) except1 = B.V2P[except1];
  if (except2 >= 0) except2 = B.V2P[except2];
  int best = -1;
  float r_dist = FLTMAX;
  for (int j = 0; j < n; j++) {
    thread const Hit& t = hits[j];
    /* PyMOL stops after the first voxel the best hit lies in front of */
    if (best >= 0 && t.voxel != hits[j - 1].voxel && VoxelAt(P, r_dist) > t.voxel + 1) break;
    Entry E = {t.v2p, t.type, t.vert, t.trans, 0};
    if (t.v2p != except1 && t.v2p != except2 && Accept(E, t.dist, r_dist, front, back, excl_trans)) {
      best = j;
      r_dist = t.dist;
    }
  }
  if (best < 0) return -1;
  SetHit(B, hits[best].vert, r_dist, hits[best].tri1, hits[best].tri2, hits[best].sph, H);
  return hits[best].vert;
}

static inline uint tobyte(float f) {
  uint c = (uint) max(f, 0.0f);
  return c > 255 ? 255 : c;
}

kernel void raytrace(constant Params& P [[buffer(0)]],
    device const Prim* prim [[buffer(1)]],
    device const float* V [[buffer(2)]],
    device const float* N [[buffer(3)]],
    device const int* V2N [[buffer(4)]],
    device const int* V2P [[buffer(5)]],
    device const int* EHead [[buffer(6)]],
    device const int* EList [[buffer(7)]],
    device const int* EMask [[buffer(8)]],
    device uint* image [[buffer(9)]],
    device float* depth [[buffer(10)]],
    device const uint4* T [[buffer(11)]],
    device const float* VMax [[buffer(12)]],
    constant int& row0 [[buffer(13)]],
    uint2 gid [[thread_position_in_grid]]) {
  int x = P.ox0 + (int) gid.x;
  int y = row0 + (int) gid.y;
  if (x >= P.ox1 || y >= P.oy1) return;
  /* image/depth hold the region; outside the traced rectangle it is the
     background and zero depth, as the host left it */
  uint idx = (uint) ((y - P.oy0) * (P.ox1 - P.ox0) + (x - P.ox0));
  if (x < P.x_start || x >= P.x_stop || y < P.y_start || y >= P.y_stop) {
    image[idx] = P.background;
    if (P.write_depth) depth[idx] = 0.0f;
    return;
  }

  Bufs B = {prim, V, N, V2N, V2P, EHead, EList, EMask, VMax, T};
  const float _0 = 0.0f, _1 = 1.0f, _p5 = 0.5f, _2 = 2.0f, _255 = 255.0f, _p499 = 0.499f;
  uint pix = P.background;
  float dep = 0.0f;
  uint last_pixel = 0;
  float bkrd[3] = {P.bkrd0, P.bkrd1, P.bkrd2};

  float b0 = (((x + 0.5f + P.border_offset)) * P.invWdthRange) + P.vol0;
  float b1 = ((y + 0.5f + P.border_offset) * P.invHgtRange) + P.vol2;

  int exclude1 = -1, exclude2 = -1;
  float persist = _1, first_excess = _0, excl_trans = _0;
  int pass = 0;
  float new_front = P.front;
  float ffact = _0, ffact1m;
  float trans = _0;

  Hit hits[MAX_HITS];
  int n_hits = -1;  /* -1: not collected yet; -2: too many, walk every pass */
  while ((persist > 0.0001f) && (pass <= P.max_pass)) {
    bool fogFlagTmp = P.fogFlag;
    bool pixel_flag = false;
    float fc[4] = {0, 0, 0, 0};
    float sn[3] = {0, 0, 0}, impact[3] = {0, 0, 0};
    HitResult H;
    int i;
    if (pass && P.collect && n_hits == -1)  /* first transparency pass: record */
      n_hits = max(CollectHits(P, B, b0, b1, new_front, P.back, hits), -2);
    if (pass && n_hits >= 0)
      i = PickHit(P, B, hits, n_hits, new_front, P.back, excl_trans, exclude1, exclude2, H);
    else
      i = BasisHitOrthoscopic(P, B, b0, b1, new_front, P.back, excl_trans, exclude1, exclude2, H);
    int ptype = 0;
    float pr1 = 0;
    int no_lighting = 0;

    if ((i >= 0) && (pass < P.max_pass)) {
      device const Prim& prm = prim[H.prim];
      ptype = prm.type;
      pr1 = prm.r1;
      no_lighting = prm.no_lighting;
      int n_basis_tmp = prm.no_lighting ? 0 : P.n_basis;
      pixel_flag = true;
      if (((trans = prm.trans) != _0) && P.trans_cont_flag)
        trans = ppow(trans, P.inv_trans_cont);
      new_front = H.dist;
      float flat_dotgle, rdotgle;

      if (prm.type == cPrimTriangle) {
        impact[0] = b0; impact[1] = b1; impact[2] = 0.0f - H.dist;
        device const float* n0 = N + (3 * V2N[i]) + 3;
        float w2 = 1.0f - (H.tri1 + H.tri2);
        fc[0] = (prm.c2[0] * H.tri1) + (prm.c3[0] * H.tri2) + (prm.c1[0] * w2);
        fc[1] = (prm.c2[1] * H.tri1) + (prm.c3[1] * H.tri2) + (prm.c1[1] * w2);
        fc[2] = (prm.c2[2] * H.tri1) + (prm.c3[2] * H.tri2) + (prm.c1[2] * w2);
        trans = (prm.tr[1] * H.tri1) + (prm.tr[2] * H.tri2) + (prm.tr[0] * w2);
        float vt1[3];
        sn[0] = n0[3] * H.tri1; sn[1] = n0[4] * H.tri1; sn[2] = n0[5] * H.tri1;
        vt1[0] = n0[6] * H.tri2; vt1[1] = n0[7] * H.tri2; vt1[2] = n0[8] * H.tri2;
        sn[0] = vt1[0] + sn[0]; sn[1] = vt1[1] + sn[1]; sn[2] = vt1[2] + sn[2];
        vt1[0] = n0[0] * w2; vt1[1] = n0[1] * w2; vt1[2] = n0[2] * w2;
        sn[0] = vt1[0] + sn[0]; sn[1] = vt1[1] + sn[1]; sn[2] = vt1[2] + sn[2];
        normalize3(sn);
        trans = ppow(trans, P.inv_trans_cont);

        if (P.has_bp2) { /* RayProjectTriangle */
          constant float* light = P.lightN;
          device const float* v0 = V + i * 3;
          int cnt = 0;
          if (dot3cd(light, n0 - 3) >= _0) cnt++;
          else if (dot3cd(light, n0) >= _0) cnt++;
          else if (dot3cd(light, n0 + 3) >= _0) cnt++;
          else if (dot3cd(light, n0 + 6) >= _0) cnt++;
          if (cnt) {
            float d1[3], d2[3], d3[3], p1[3], p2[3], p3[3];
            float pw2 = 1.0f - (H.tri1 + H.tri2);
            for (int k = 0; k < 3; k++) {
              d1[k] = v0[k] - impact[k];
              d2[k] = v0[3 + k] - impact[k];
              d3[k] = v0[6 + k] - impact[k];
            }
            float t = dot3(d1, n0);
            p1[0] = n0[0] * t; p1[1] = n0[1] * t; p1[2] = n0[2] * t;
            t = dot3(d2, n0 + 3);
            p2[0] = n0[3] * t; p2[1] = n0[4] * t; p2[2] = n0[5] * t;
            t = dot3(d3, n0 + 6);
            p3[0] = n0[6] * t; p3[1] = n0[7] * t; p3[2] = n0[8] * t;
            for (int k = 0; k < 3; k++) {
              d1[k] = p1[k] * pw2;
              d2[k] = p2[k] * H.tri1;
              d3[k] = p3[k] * H.tri2;
            }
            for (int k = 0; k < 3; k++) d2[k] = d1[k] + d2[k];
            for (int k = 0; k < 3; k++) d3[k] = d2[k] + d3[k];
            for (int k = 0; k < 3; k++) d3[k] = d3[k] * P.project_triangle;
            if ((sn[0] * d3[0] + sn[1] * d3[1] + sn[2] * d3[2]) >= _0)
              for (int k = 0; k < 3; k++) impact[k] = d3[k] + impact[k];
          }
        }
        rdotgle = -sn[2];
        flat_dotgle = N[3 * V2N[i] + 2];
      } else {
        impact[0] = b0; impact[1] = b1; impact[2] = 0.0f - H.dist;
        sn[0] = impact[0] - H.sphere[0];
        sn[1] = impact[1] - H.sphere[1];
        sn[2] = impact[2] - H.sphere[2];
        normalize3(sn);
        rdotgle = -sn[2];
        flat_dotgle = sn[2];
        if (prm.type == cPrimCylinder || prm.type == cPrimSausage) {
          float ft = H.tri1;
          fc[0] = (prm.c1[0] * (_1 - ft)) + (prm.c2[0] * ft);
          fc[1] = (prm.c1[1] * (_1 - ft)) + (prm.c2[1] * ft);
          fc[2] = (prm.c1[2] * (_1 - ft)) + (prm.c2[2] * ft);
        } else {
          fc[0] = prm.c1[0]; fc[1] = prm.c1[1]; fc[2] = prm.c1[2];
        }
      }

      if ((P.trans_oblique != _0) && (trans != _0)) {
        if ((sn[2] > _0) || P.two_sided_lighting) {
          float oblique_factor = sn[2];
          if (oblique_factor < _0) oblique_factor = -oblique_factor;
          if (oblique_factor != _1) {
            if (oblique_factor > _p5)
              oblique_factor = (_p5 + _p5 * (_1 - ppow((_1 - oblique_factor) * _2, P.oblique_power)));
            else
              oblique_factor = (_p5 * ppow(oblique_factor * _2, P.oblique_power));
          }
          trans *= (P.trans_oblique * oblique_factor) + (1.0f - P.trans_oblique);
          if (trans < 0.06f) trans = 0.06f;
        }
      }

      float dotgle = -rdotgle;
      if (flat_dotgle < _0) {
        if (P.two_sided_lighting) {
          dotgle = -dotgle;
          sn[0] = -sn[0]; sn[1] = -sn[1]; sn[2] = -sn[2];
        }
      }

      float direct_cmp;
      {
        float pow_dotgle, pow_surfnormal2;
        if (P.settingPower != _1) {
          pow_dotgle = ppow(dotgle, P.settingPower);
          pow_surfnormal2 = ppow(sn[2], P.settingPower);
        } else {
          pow_dotgle = dotgle;
          pow_surfnormal2 = sn[2];
        }
        direct_cmp = P.legacy_1m * pow_surfnormal2 + P.legacy * ((dotgle + pow_dotgle) * _p5);
      }

      float reflect_cmp = _0, excess;
      if (P.settingSpecDirect != _0) {
        if (sn[2] > _0)
          excess = ppow(sn[2], P.settingSpecDirectPower) * P.settingSpecDirect;
        else
          excess = _0;
      } else {
        excess = _0;
      }

      float lit = _1;
      if (n_basis_tmp < 3) {
        reflect_cmp = direct_cmp;
      } else {
        for (int bc = 2; bc < P.n_basis; bc++) {
          lit = _1;
          constant float* L = P.lightN + (bc - 2) * 3;
          constant float* S = P.specN + (bc - 2) * 3;
          {
            float pow_dotgle;
            dotgle = -dot3c(sn, L);
            if (dotgle < _0) dotgle = _0;
            if (P.settingReflectPower != _1)
              pow_dotgle = ppow(dotgle, P.settingReflectPower);
            else
              pow_dotgle = dotgle;
            reflect_cmp += P.legacy_1m * (lit * pow_dotgle) + P.legacy * (lit * (dotgle + pow_dotgle) * _p5);
          }
          if (bc < (P.spec_count + 2)) {
            if (P.ray_scatter != _0) excess += P.ray_scatter * dotgle;
            dotgle = -dot3c(sn, S);
            if (dotgle < _0) dotgle = _0;
            excess += ppow(dotgle, P.settingSpecPower) * P.settingSpecReflect * lit;
          }
        }
      }

      float bright = P.ambient + (((_1 - P.direct_shade) + P.direct_shade * lit) * P.direct * direct_cmp +
                                  P.lreflect * reflect_cmp * (P.legacy_1m + P.legacy * direct_cmp));
      if (excess > _1) excess = _1;
      if (bright > _1) bright = _1;
      else if (bright < _0) bright = _0;

      fc[0] = (bright * fc[0] + excess);
      fc[1] = (bright * fc[1] + excess);
      fc[2] = (bright * fc[2] + excess);

      if (n_basis_tmp && fogFlagTmp) {
        ffact = (P.front - H.dist) * P.invFrontMinusBack;
        if (P.fogRangeFlag) ffact = (ffact - P.fog_start) * P.inv1minusFogStart;
        ffact *= P.fog;
        if (ffact < _0) ffact = _0;
        if (ffact > _1) ffact = _1;
        ffact1m = _1 - ffact;
        if (P.orig_opaque_back) {
          fc[0] = ffact * bkrd[0] + fc[0] * ffact1m;
          fc[1] = ffact * bkrd[1] + fc[1] * ffact1m;
          fc[2] = ffact * bkrd[2] + fc[2] * ffact1m;
          fc[3] = _1;
        } else {
          fc[3] = ffact1m * (_1 - trans);
        }
        if (!pass) {
          if (trans < P.trans_spec_cut)
            first_excess = excess * ffact1m * P.ray_trans_spec;
          else
            first_excess = excess * ffact1m * P.ray_trans_spec * P.trans_spec_scale * (_1 - trans);
        } else {
          fc[0] += first_excess; fc[1] += first_excess; fc[2] += first_excess;
        }
      } else {
        if (!pass) {
          if (trans < P.trans_spec_cut)
            first_excess = excess * P.ray_trans_spec;
          else
            first_excess = excess * P.ray_trans_spec * P.trans_spec_scale * (_1 - trans);
        } else {
          fc[0] += first_excess; fc[1] += first_excess; fc[2] += first_excess;
        }
        if (P.orig_opaque_back)
          fc[3] = _1;
        else
          fc[3] = _1 - trans;
      }
    } else if (pass) {
      i = -1;
      fc[0] = first_excess + bkrd[0];
      fc[1] = first_excess + bkrd[1];
      fc[2] = first_excess + bkrd[2];
      fc[3] = P.orig_opaque_back ? _1 : _0;
      ffact = 1.0f;
      ffact1m = 0.0f;
      pixel_flag = true;
      if (P.trans_cont_flag) persist = ppow(persist, P.trans_cont);
    }

    if (pixel_flag) {
      uint cc0 = tobyte(fc[0] * _255), cc1 = tobyte(fc[1] * _255), cc2 = tobyte(fc[2] * _255);
      if (P.orig_opaque_back) {
        pix = P.fore_mask | (cc2 << 16) | (cc1 << 8) | cc0;
      } else {
        uint cc3 = tobyte(fc[3] * _255);
        pix = (cc3 << 24) | (cc2 << 16) | (cc1 << 8) | cc0;
      }
    }

    if (pass && persist < 1.f) {
      float mix_in;
      if (i >= 0) {
        if (fogFlagTmp) {
          if (P.trans_cont_flag && (ffact > _p5))
            mix_in = 2 * (persist * (_1 - ffact) + (ppow(persist, P.trans_cont) * (ffact - _p5))) * (_1 - trans * ffact);
          else
            mix_in = persist * (_1 - trans * ffact);
        } else {
          mix_in = persist * (_1 - trans);
        }
      } else {
        mix_in = persist;
      }
      float persist_inv = _1 - mix_in;
      if (!P.orig_opaque_back) {
        if (i < 0) {
          fc[0] = (float) (0xFF & (last_pixel >> 24));
          fc[1] = (float) (0xFF & (last_pixel >> 16));
          fc[2] = (float) (0xFF & (last_pixel >> 8));
          fc[3] = (float) (0xFF & (last_pixel));
          if (P.trans_cont_flag) {
            float m = _1 - (float) (0xFF & (last_pixel >> 24)) / _255;
            m = _1 - ppow(m, P.trans_cont);
            fc[0] = m * _255 + _p499;
          }
        } else {
          fc[0] = (0xFF & (pix >> 24)) * mix_in + (0xFF & (last_pixel >> 24)) * persist_inv;
          fc[1] = (0xFF & (pix >> 16)) * mix_in + (0xFF & (last_pixel >> 16)) * persist_inv;
          fc[2] = (0xFF & (pix >> 8)) * mix_in + (0xFF & (last_pixel >> 8)) * persist_inv;
          fc[3] = (0xFF & (pix)) * mix_in + (0xFF & (last_pixel)) * persist_inv;
          {
            float o1 = (float) (0xFF & (last_pixel >> 24)) / _255;
            float o2 = (float) (0xFF & (pix >> 24)) / _255;
            if (o1 < o2) { float m = o1; o1 = o2; o2 = m; }
            float m = o1 + (1.0f - o1) * o2;
            fc[0] = m * _255 + _p499;
          }
        }
      } else {
        fc[0] = (0xFF & (pix >> 24)) * mix_in + (0xFF & (last_pixel >> 24)) * persist_inv;
        fc[1] = (0xFF & (pix >> 16)) * mix_in + (0xFF & (last_pixel >> 16)) * persist_inv;
        fc[2] = (0xFF & (pix >> 8)) * mix_in + (0xFF & (last_pixel >> 8)) * persist_inv;
        fc[3] = (0xFF & (pix)) * mix_in + (0xFF & (last_pixel)) * persist_inv;
      }
      pix = (tobyte(fc[0]) << 24) | (tobyte(fc[1]) << 16) | (tobyte(fc[2]) << 8) | tobyte(fc[3]);
    }

    if (P.write_depth && (i >= 0) && (trans < P.trans_cutoff) && (persist > P.persist_cutoff))
      dep = (P.front + impact[2]);

    if (i >= 0) {
      if (ptype == cPrimSausage)
        excl_trans = new_front + (2 * sn[2] * pr1);
      if ((!P.backface_cull) && (P.trans_mode != 2)) {
        persist = persist * trans;
      } else {
        if (!no_lighting && (persist < 0.9999f) && (trans > 0.05f))
          pix = last_pixel;
        else
          persist = persist * trans;
      }
    }

    if (i < 0) break;
    last_pixel = pix;
    exclude2 = exclude1;
    exclude1 = i;
    pass++;
  }

  if (P.blend_colors) {
    float fc0, fc1, fc2;
    uint cc3 = (0xFF & (pix >> 24));
    fc2 = (float) (0xFF & (pix >> 16));
    fc1 = (float) (0xFF & (pix >> 8));
    fc0 = (float) (0xFF & (pix));
    float red_part = P.red_blend * fc0;
    float green_part = P.green_blend * fc1;
    float blue_part = P.blue_blend * fc2;
    float red_min = (green_part > blue_part) ? green_part : blue_part;
    float green_min = (red_part > blue_part) ? red_part : blue_part;
    float blue_min = (green_part > red_part) ? green_part : red_part;
    if (fc0 < red_min) fc0 = red_min;
    if (fc1 < green_min) fc1 = green_min;
    if (fc2 < blue_min) fc2 = blue_min;
    pix = (cc3 << 24) | (tobyte(fc2) << 16) | (tobyte(fc1) << 8) | tobyte(fc0);
  }
  image[idx] = pix;
  if (P.write_depth) depth[idx] = dep;
}
/* --- RayRender's ray_trace_mode outline pass, line for line ------------ */

struct OParams {
  int width, height, ox0, oy0, ox1, oy1;
  int mode, fogFlag, fogRangeFlag, opaque;
  uint trace_word, background;
  float fog, fog_start, inv1minusFogStart, invFrontMinusBack;
  float slope_f, depth_f, disco_f, disco_raw;
};

/* a region buffer's value at absolute (x, y); outside the region it is 0,
   as on the host */
static inline float at(constant OParams& O, device const float* b, int x, int y, int k = 0, int n = 1) {
  if (x < O.ox0 || x >= O.ox1 || y < O.oy0 || y >= O.oy1) return 0.0f;
  return b[((y - O.oy0) * (O.ox1 - O.ox0) + (x - O.ox0)) * n + k];
}

kernel void outline_grad(constant OParams& O [[buffer(0)]], device const float* depth [[buffer(1)]],
    device float* delta [[buffer(2)]], uint2 gid [[thread_position_in_grid]]) {
  int x = O.ox0 + (int) gid.x, y = O.oy0 + (int) gid.y;
  if (x >= O.ox1 || y >= O.oy1) return;
  float dzdx = 0.0F, dzdy = 0.0F, dd, d = at(O, depth, x, y);
  int xc = 0, yc = 0;
  if (x) { dd = d - at(O, depth, x - 1, y); dzdx += dd; xc++; }
  if (x < (O.width - 1)) {
    dd = at(O, depth, x + 1, y) - d;
    if ((!xc) || (fabs(dd) > fabs(dzdx))) dzdx = dd;
  }
  if (y) { dd = d - at(O, depth, x, y - 1); dzdy += dd; yc++; }
  if (y < (O.height - 1)) {
    dd = at(O, depth, x, y + 1) - d;
    if ((!yc) || (fabs(dd) > fabs(dzdy))) dzdy = dd;
  }
  device float* q = delta + ((y - O.oy0) * (O.ox1 - O.ox0) + (x - O.ox0)) * 3;
  q[0] = dzdx;
  q[1] = dzdy;
  q[2] = sqrt1f(fma(dzdx, dzdx, dzdy * dzdy));  /* fused as the host compiler does */
}

/* one depth fill pass: give depthless pixels their neighbours' mean; the
   rest of dst is left as it was (the host ping-pongs the same way) */
kernel void outline_fill(constant OParams& O [[buffer(0)]], device const float* src [[buffer(1)]],
    device float* dst [[buffer(2)]], uint2 gid [[thread_position_in_grid]]) {
  int x = O.ox0 + (int) gid.x, y = O.oy0 + (int) gid.y;
  if (x >= O.ox1 || y >= O.oy1) return;
  const float R_SMALL4 = 0.0001F;
  float p = at(O, src, x, y);
  device float* q = dst + (y - O.oy0) * (O.ox1 - O.ox0) + (x - O.ox0);
  if (fabs(p) < R_SMALL4) {
    float dep = 0.0F, n;
    int cnt = 0;
    if (x && fabs(n = at(O, src, x - 1, y)) > R_SMALL4) { dep += n; cnt++; }
    if (x < (O.width - 1) && fabs(n = at(O, src, x + 1, y)) > R_SMALL4) { dep += n; cnt++; }
    if (y && fabs(n = at(O, src, x, y - 1)) > R_SMALL4) { dep += n; cnt++; }
    if (y < (O.height - 1) && fabs(n = at(O, src, x, y + 1)) > R_SMALL4) { dep += n; cnt++; }
    if (cnt) {
      dep = pdiv(dep, (float) cnt);
      *q = dep;
    }
  } else {
    *q = p;
  }
}

kernel void outline_edges(constant OParams& O [[buffer(0)]], device const float* delta [[buffer(1)]],
    device const float* depth [[buffer(2)]], device uint* image [[buffer(3)]],
    uint2 gid [[thread_position_in_grid]]) {
  int x = O.ox0 + (int) gid.x, y = O.oy0 + (int) gid.y;
  if (x >= O.ox1 || y >= O.oy1) return;
  const float _0 = 0.0F, _1 = 1.0F, R_SMALL4 = 0.0001F;
  const float _8 = 0.08F, _4 = 0.4F, _25 = 0.25F, _m25 = -0.25F;
  /* as in RayRender: from the factor before the gain is applied */
  float disco_f_625 = O.disco_raw * 0.625F, disco_f_5 = O.disco_raw * 0.5F, disco_f_45 = O.disco_raw * 0.45F;
  float max_slope = 0.0F, max_depth = 0.0F, min_dot = 1.0F, max_dz = 0.0F, max_pz = 0.0F;
  float dx = at(O, delta, x, y, 0, 3), dy = at(O, delta, x, y, 1, 3), dz = at(O, delta, x, y, 2, 3);
  float px = 0.0F, py = 0.0F, pz = 0.0F, ddx, ddy, diff, dot;
  /* the CPU loop carries px/py/pz across pixels; only x == 0 reads it:
     case 5 of the row above's last pixel */
  if (x == 0 && y && O.width > 1) {
    px = at(O, delta, O.width - 2, y, 0, 3); py = at(O, delta, O.width - 2, y, 1, 3); pz = at(O, delta, O.width - 2, y, 2, 3);
  }
  for (int i = 0; i < 8; i++) {
    int nx = x, ny = y;
    bool ok = false;
    switch (i) {
    case 0: ok = x; nx = x - 1; break;
    case 1: ok = x < (O.width - 1); nx = x + 1; break;
    case 2: ok = y; ny = y - 1; break;
    case 3: ok = y < (O.height - 1); ny = y + 1; break;
    case 4: ok = x && y; nx = x - 1; ny = y - 1; break;
    case 5: ok = x && (y < (O.height - 1)); nx = x - 1; ny = y + 1; break;
    case 6: ok = y && (x < (O.width - 1)); nx = x + 1; ny = y - 1; break;
    case 7: ok = (y < (O.height - 1)) && (x < (O.width - 1)); nx = x + 1; ny = y + 1; break;
    }
    if (ok) {
      px = at(O, delta, nx, ny, 0, 3); py = at(O, delta, nx, ny, 1, 3); pz = at(O, delta, nx, ny, 2, 3);
    }
    ddx = dx - px;
    ddy = dy - py;
    diff = fma(ddx, ddx, ddy * ddy);  /* the host's fused multiply-adds */
    if (max_depth < diff) max_depth = diff;
    if ((dz > R_SMALL4) && (pz > R_SMALL4)) {
      dot = fma(pdiv(dx, dz), pdiv(px, pz), pdiv(dy, dz) * pdiv(py, pz));
      if (dot < min_dot) { min_dot = dot; max_dz = dz; max_pz = pz; }
    }
    diff = fabs(dz - pz);
    if (diff > max_slope) max_slope = diff;
  }
  device uint* q = image + (y - O.oy0) * (O.ox1 - O.ox0) + (x - O.ox0);
  if ((max_slope > (O.slope_f)) || (max_depth > (O.depth_f)) ||
      ((min_dot < _8) && ((max_dz > O.disco_f) || (max_pz > O.disco_f))) ||
      ((min_dot < _4) && ((max_dz > disco_f_625) || (max_pz > disco_f_625))) ||
      ((min_dot < _25) && ((max_dz > disco_f_5) || (max_pz > disco_f_5))) ||
      ((min_dot < _m25) && ((max_dz > disco_f_45) && (max_pz > disco_f_45)))) {
    if (O.fogFlag) {
      float ffact = at(O, depth, x, y) * O.invFrontMinusBack, ffact1m, fc[4];
      if (O.fogRangeFlag) ffact = (ffact - O.fog_start) * O.inv1minusFogStart;
      ffact *= O.fog;
      if (ffact < _0) ffact = _0;
      if (ffact > _1) ffact = _1;
      ffact1m = _1 - ffact;
      uint bg = O.background, tw = O.trace_word;
      if (O.opaque) {
        fc[0] = fma((float) (0xFF & (bg >> 24)), ffact, (0xFF & (tw >> 24)) * ffact1m);
        fc[1] = fma((float) (0xFF & (bg >> 16)), ffact, (0xFF & (tw >> 16)) * ffact1m);
        fc[2] = fma((float) (0xFF & (bg >> 8)), ffact, (0xFF & (tw >> 8)) * ffact1m);
        fc[3] = fma((float) (0xFF & (bg)), ffact, (0xFF & (tw)) * ffact1m);
      } else {
        fc[1] = (0xFF & (tw >> 16));
        fc[2] = (0xFF & (tw >> 8));
        fc[0] = fma((float) (0xFF & (bg >> 24)), ffact, (0xFF & (tw >> 24)) * ffact1m);
        fc[3] = (0xFF & (tw));
      }
      *q = (tobyte(fc[0]) << 24) | (tobyte(fc[1]) << 16) | (tobyte(fc[2]) << 8) | tobyte(fc[3]);
    } else {
      *q = O.trace_word;
    }
  } else if (O.mode == 2) {
    *q = O.background;
  } else if (O.mode == 3) {
    *q = (*q & 0xC0C0C0C0);
    *q = *q | ((*q) >> 2) | ((*q) >> 4) | ((*q) >> 6);
  }
}
/* RayAntiThread's 3x case: a 5x5 alpha-weighted filter (weights 1/5/8, sum
   64) per output pixel; with no alpha, a plain weighted RGB mean whose last
   row has only four taps, as in the original */
kernel void antialias3(constant OParams& O [[buffer(0)]], device const uint* image [[buffer(1)]],
    device uint* out [[buffer(2)]], uint2 gid [[thread_position_in_grid]]) {
  int ow = O.width / 3 - 2, oh = O.height / 3 - 2;
  if ((int) gid.x >= ow || (int) gid.y >= oh) return;
  const uint w[5][5] = {{1, 1, 1, 1, 1}, {1, 5, 5, 5, 1}, {1, 5, 8, 5, 1}, {1, 5, 5, 5, 1}, {1, 1, 1, 1, 1}};
  uint c1 = 0, c2 = 0, c3 = 0, c4 = 0, px[5][5];
  for (int j = 0; j < 5; j++)
    for (int i = 0; i < 5; i++) {
      int x = (int) gid.x * 3 + i, y = (int) gid.y * 3 + j;
      px[j][i] = (x < O.ox0 || x >= O.ox1 || y < O.oy0 || y >= O.oy1) ? O.background
                 : image[(y - O.oy0) * (O.ox1 - O.ox0) + (x - O.ox0)];
      uint a = (px[j][i] >> 24) * w[j][i];
      c4 += a;
      c1 += (px[j][i] & 0xFF) * a;
      c2 += ((px[j][i] >> 8) & 0xFF) * a;
      c3 += ((px[j][i] >> 16) & 0xFF) * a;
    }
  if (c4) {
    c1 /= c4; c2 /= c4; c3 /= c4;
  } else {
    c1 = c2 = c3 = 0;
    for (int j = 0; j < 5; j++)
      for (int i = 0; i < (j == 4 ? 4 : 5); i++) {
        c1 += w[j][i] * (px[j][i] & 0xFF);
        c2 += w[j][i] * ((px[j][i] >> 8) & 0xFF);
        c3 += w[j][i] * ((px[j][i] >> 16) & 0xFF);
      }
    c1 >>= 6; c2 >>= 6; c3 >>= 6;
  }
  out[gid.y * ow + gid.x] = (c1 & 0xFF) | ((c2 & 0xFF) << 8) | ((c3 & 0xFF) << 16) | (((c4 >> 6) & 0xFF) << 24);
}
)MSL";

/* --- host ---------------------------------------------------------------- */

static int s_enabled = 1; /* `ray` uses the GPU whenever it can */
static std::string s_report = "not run", s_stages;

void RayMetalSetEnabled(int on) { s_enabled = on ? 1 : 0; }
int RayMetalEnabled() { return s_enabled; }
void RayMetalSetLastReport(const char* msg) { s_report = msg ? msg : ""; }
const char* RayMetalLastReport() { return s_report.c_str(); }
void RayMetalSetStages(const char* msg) { s_stages = msg ? msg : ""; }
const char* RayMetalStages() { return s_stages.c_str(); }

static id<MTLDevice> s_device = nil;
static id<MTLCommandQueue> s_queue = nil;
static id<MTLComputePipelineState> s_pipeline = nil, s_grad = nil, s_fill = nil, s_edges = nil,
                                   s_aa = nil;
static std::string s_init_err;

static bool RayMetalInit(std::string& err) {
  static std::once_flag once;
  std::call_once(once, [] {
    @autoreleasepool {
      s_device = MTLCreateSystemDefaultDevice();
      if (!s_device) {
        s_init_err = "no Metal device";
        return;
      }
      s_queue = [s_device newCommandQueue];
      MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
      if (@available(macOS 15.0, *))
        opts.mathMode = MTLMathModeSafe;
      else
        opts.fastMathEnabled = NO;
      NSError* e = nil;
      id<MTLLibrary> lib = [s_device newLibraryWithSource:@(kRayMetalSource) options:opts error:&e];
      RM_RELEASE(opts);
      auto make = [&](NSString* name) -> id<MTLComputePipelineState> {
        id<MTLFunction> fn = [lib newFunctionWithName:name];
        id<MTLComputePipelineState> ps =
            fn ? [s_device newComputePipelineStateWithFunction:fn error:&e] : nil;
        RM_RELEASE(fn);
        return ps;
      };
      if (lib) {
        s_grad = make(@"outline_grad");
        s_fill = make(@"outline_fill");
        s_edges = make(@"outline_edges");
        s_aa = make(@"antialias3");
        s_pipeline = (s_grad && s_fill && s_edges && s_aa) ? make(@"raytrace") : nil;
        RM_RELEASE(lib);
      }
      if (!s_pipeline)
        s_init_err =
            std::string("kernel build failed: ") + (e ? e.localizedDescription.UTF8String : "?");
    }
  });
  if (!s_pipeline) err = s_init_err;
  return s_pipeline != nil;
}

static id<MTLBuffer> RMBuf(const void* p, size_t bytes) {
  static const uint32_t zeros[4] = {0, 0, 0, 0};
  if (!p || !bytes) p = zeros, bytes = sizeof(zeros);
  return [s_device newBufferWithBytes:p length:bytes options:MTLResourceStorageModeShared];
}

static void ParallelFor(size_t n, void (^fn)(size_t)) {
  dispatch_apply(n, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), fn);
}

/* Each grid column's voxel lists for the kernel. An orthoscopic ray starts in
   the same voxel layer whatever its pixel, so all rays through column (a, b)
   walk the same voxels in the same order, and PyMOL's per-ray MapCache ("test
   each primitive once") skips the same entries for all of them: keep only an
   entry's first appearance in its column, in PyMOL's order. Also drop entries
   whose (x, y) footprint misses the column, as no ray there can hit them —
   unless an end-on cylinder exists, whose test reuses the previous test's
   tri1/sphere (ZLineToSphere[Capped]), so removing tests could change it.
   A lite map has no express lists: its voxel list is MapSetupExpressXY's, the
   3x3 (x, y) neighbours' grid chains in that layer, remapped to Vertex. */
static bool BuildColumnLists(const RayMetalScene& s, std::vector<int32_t>& head,
                             std::vector<int32_t>& list, std::vector<int32_t>& mask,
                             std::vector<float>& vmax) {
  const RayMetalParams& p = s.p;
  const int d0 = p.dim0, d1 = p.dim1, d2 = p.dim2, nv = p.n_vert, np = s.n_prims;
  int c0 = (int)((0.0f - p.mapMin2) * p.recipDiv) + 2 + 1;
  c0 = std::min(std::max(c0, p.iMin2), p.iMax2);
  head.assign((size_t)d0 * d1 * d2, 0);
  vmax.assign((size_t)d0 * d1 * d2, 0.f);
  mask.assign((size_t)d0 * d1, 0);

  /* per prim (x0, x1, y0, y1) footprint and the farthest dist it can hit */
  std::vector<float> box((size_t)np * 5);
  float* BX = box.data();
  __block bool endOn = false, zeroLen = false;
  const float fud = std::fabs(p.fudge1 - 1.0f), inf = 1e30f;
  ParallelFor((np + 4095) / 4096, ^(size_t ck) {
    for (size_t k = ck * 4096; k < std::min<size_t>(np, ck * 4096 + 4096); k++) {
      const RayMetalPrim& pr = s.prims[k];
      const float* v = s.vertex + (size_t)pr.vert * 3;
      float* b = BX + k * 5;
      float x0 = -inf, x1 = inf, y0 = -inf, y1 = inf, z0 = -inf, pad = pr.r1;
      if (pr.type == 1) {
        x0 = x1 = v[0];
        y0 = y1 = v[1];
        z0 = v[2];
      } else if (pr.type == 3) {
        x0 = std::min({v[0], v[3], v[6]});
        x1 = std::max({v[0], v[3], v[6]});
        y0 = std::min({v[1], v[4], v[7]});
        y1 = std::max({v[1], v[4], v[7]});
        z0 = std::min({v[2], v[5], v[8]});
        float z1 = std::max({v[2], v[5], v[8]});
        /* the barycentric test accepts up to ray_triangle_fudge outside */
        pad = 2.f * fud * ((x1 - x0) + (y1 - y0) + (z1 - z0)) + 1e-3f;
      } else if (pr.type == 2 || pr.type == 4) {
        const float* d = s.normal + s.vert2normal[pr.vert] * 3;
        if (std::fabs(d[2]) > 1.0f - 2e-4f) endOn = true;
        if (pr.l1 <= 1e-4f) zeroLen = true;
        float ex = v[0] + d[0] * pr.l1, ey = v[1] + d[1] * pr.l1;
        x0 = std::min(v[0], ex);
        x1 = std::max(v[0], ex);
        y0 = std::min(v[1], ey);
        y1 = std::max(v[1], ey);
        z0 = std::min(v[2], v[2] + d[2] * pr.l1);
      } else {
        pad = 0;
      }
      b[0] = x0 - pad;
      b[1] = x1 + pad;
      b[2] = y0 - pad;
      b[3] = y1 + pad;
      /* rays run along -z from z = 0: no hit lies beyond the farthest point
         (padded for rounding) */
      b[4] = -z0 + pad + 1e-3f + 1e-5f * std::fabs(z0);
    }
  });
  const bool filter = !endOn;
  /* a ray's column is (int) ((x - Min) * recipDiv) + 2, clamped by one voxel
     at the grid edges: pad the column's bounds generously */
  const float vox = 1.0f / p.recipDiv;
  auto bounds = [&](int a, float mn, int imin, int imax, float& lo, float& hi) {
    lo = mn + (a - 2.5f) * vox - ((a <= imin + 1 || a <= 3) ? 2 * vox : 0);
    hi = mn + (a - 0.5f) * vox + ((a >= imax - 1) ? 2 * vox : 0);
  };

  std::vector<std::vector<int32_t>> rows(d0);
  /* blocks capture C++ objects as const copies: use raw pointers */
  std::vector<int32_t>* RW = rows.data();
  int32_t *HD = head.data(), *MK = mask.data();
  float* VM = vmax.data();
  const int nw = (int)[[NSProcessInfo processInfo] activeProcessorCount];
  ParallelFor(nw, ^(size_t w) {
    std::vector<int32_t> stamp(std::max(np, 1), -1);
    for (int a = (int)w + p.iMin0; a <= p.iMax0; a += nw) {
      std::vector<int32_t>& L = RW[a];
      float xl, xh, yl, yh;
      bounds(a, p.mapMin0, p.iMin0, p.iMax0, xl, xh);
      for (int b = p.iMin1; b <= p.iMax1; b++) {
        if (!s.lite && !s.emask[d1 * a + b]) continue;
        bounds(b, p.mapMin1, p.iMin1, p.iMax1, yl, yh);
        const int32_t col = a * d1 + b;
        float far = -inf;
        auto emit = [&](int i) {
          int pr = s.vert2prim[i];
          if (stamp[pr] == col) return;
          stamp[pr] = col;
          const float* bb = BX + (size_t)pr * 5;
          if (filter && (bb[1] < xl || bb[0] > xh || bb[3] < yl || bb[2] > yh)) return;
          L.push_back(i);
          far = filter ? std::max(far, bb[4]) : inf;
        };
        for (int c = c0; c >= 2; c--) {
          size_t idx = ((size_t)a * d1 + b) * d2 + c, start = L.size();
          far = -inf;
          if (s.lite) {
            for (int dd = a - 1; dd <= a + 1; dd++)
              for (int ee = b - 1; ee <= b + 1; ee++)
                for (int i = s.mhead[((size_t)dd * d1 + ee) * d2 + c]; i >= 0; i = s.mlink[i])
                  emit(i >= nv ? s.literef[i] : i);
          } else {
            int h = s.ehead[idx];
            if (h > 0 && h < p.n_eElem)
              for (int ip = h; s.elist[ip] >= 0 && s.elist[ip] < nv; ip++) emit(s.elist[ip]);
          }
          if (L.size() > start) {
            L.push_back(-1);
            HD[idx] = (int32_t)start + 1; /* row-local + 1 for now */
            VM[idx] = far;
            MK[(size_t)d1 * a + b] = 1;
          }
        }
      }
    }
  });
  /* concatenate the rows behind a dummy entry 0, rebasing the heads */
  std::vector<size_t> off(d0 + 1, 1);
  for (int a = 0; a < d0; a++) off[a + 1] = off[a] + rows[a].size();
  list.assign(off[d0] + 1, -1);
  int32_t* LS = list.data();
  const size_t* OFF = off.data();
  ParallelFor(d0, ^(size_t a) {
    std::copy(RW[a].begin(), RW[a].end(), LS + OFF[a]);
    for (size_t k = a * d1 * d2; k < (a + 1) * d1 * d2; k++)
      if (HD[k]) HD[k] += (int32_t)OFF[a] - 1;
  });
  /* collecting hits for transparency passes (kernel) is exact when no test
     reads another's leftovers */
  return filter && !zeroLen;
}

/* The kernel's 64-byte test record for each map vertex (layout above
   BasisHitOrthoscopic), written as raw words: the int fields as ints. */
static void BuildRecords(const RayMetalScene& s, uint32_t* R) {
  auto fb = [](float f) {
    uint32_t u;
    memcpy(&u, &f, 4);
    return u;
  };
  const int nv = s.p.n_vert;
  ParallelFor((nv + 16383) / 16384, ^(size_t ck) {
    for (int i = (int)ck * 16384; i < std::min(nv, (int)ck * 16384 + 16384); i++) {
      uint32_t* r = R + (size_t)i * 16;
      memset(r, 0, 64);
      const RayMetalPrim& pr = s.prims[s.vert2prim[i]];
      r[3] = (uint32_t)s.vert2prim[i] | ((uint32_t)(pr.type & 0xF) << 27) |
             ((uint32_t)(pr.cull ? 1 : 0) << 31);
      r[7] = fb(pr.trans);
      r[15] = (uint32_t)pr.vert;
      if (pr.type == 3) { /* triangle */
        const float* pre = s.precomp + s.vert2normal[i] * 3;
        const float* v0 = s.vertex + (size_t)pr.vert * 3;
        r[0] = fb(v0[0]);
        r[1] = fb(v0[1]);
        r[2] = fb(v0[2]);
        r[4] = fb(pre[6]);
        r[5] = fb(pre[7]);
        r[8] = fb(pre[0]);
        r[9] = fb(pre[1]);
        r[10] = fb(pre[3]);
        r[12] = fb(pre[4]);
        r[13] = fb(pre[2]);
        r[14] = fb(pre[5]);
      } else {
        const float* v = s.vertex + (size_t)i * 3;
        r[0] = fb(v[0]);
        r[1] = fb(v[1]);
        r[2] = fb(v[2]);
        r[4] = fb(s.radius[i]);
        r[5] = fb(s.radius2[i]);
        r[6] = fb(pr.l1);
        if (pr.type != 1) { /* cylinder / sausage; caps masked: triangles leave them unset */
          const float* n = s.normal + s.vert2normal[i] * 3;
          const float* pre = s.precomp + s.vert2normal[i] * 3;
          r[8] = fb(n[0]);
          r[9] = fb(n[1]);
          r[10] = fb(n[2]);
          r[11] = ((uint32_t)pr.cap1 & 0xF) | (((uint32_t)pr.cap2 & 0xF) << 4);
          r[12] = fb(pre[0]);
          r[13] = fb(pre[1]);
        }
      }
    }
  });
}

bool RayMetalTrace(const RayMetalScene& s, unsigned int* image, float* depth,
                   const RayMetalOutline* outline, unsigned int* aa3, std::string& err,
                   double* gpu_ms) {
  if (!RayMetalInit(err)) return false;
  if (s.n_prims >= (1 << 27)) {
    err = "too many primitives";
    return false;
  }
  const RayMetalParams& p = s.p;
  int w = p.ox1 - p.ox0, h = p.oy1 - p.oy0; /* the region the GPU writes */
  if (w <= 0 || h <= 0 || p.x_stop <= p.x_start || p.y_stop <= p.y_start) {
    *gpu_ms = 0;
    return !outline;
  }

  @autoreleasepool {
    std::vector<int32_t> head, list, mask;
    std::vector<float> vmax;
    RayMetalParams pk = p;
    pk.collect = BuildColumnLists(s, head, list, mask, vmax);
    pk.n_eElem = (int32_t)list.size();
    size_t npix = (size_t)w * h;
    id<MTLBuffer> rec = [s_device newBufferWithLength:std::max<size_t>(64, (size_t)p.n_vert * 64)
                                              options:MTLResourceStorageModeShared];
    BuildRecords(s, (uint32_t*)rec.contents);
    id<MTLBuffer> bufs[13] = {nil,
                              RMBuf(s.prims, sizeof(RayMetalPrim) * s.n_prims),
                              RMBuf(s.vertex, sizeof(float) * s.n_vertex_f),
                              RMBuf(s.normal, sizeof(float) * s.n_normal_f),
                              RMBuf(s.vert2normal, sizeof(int32_t) * s.n_vert2normal),
                              RMBuf(s.vert2prim, sizeof(int32_t) * s.n_vert2prim),
                              RMBuf(head.data(), sizeof(int32_t) * head.size()),
                              RMBuf(list.data(), sizeof(int32_t) * list.size()),
                              RMBuf(mask.data(), sizeof(int32_t) * mask.size()),
                              [s_device newBufferWithLength:sizeof(uint32_t) * npix
                                                    options:MTLResourceStorageModeShared],
                              depth ? [s_device newBufferWithLength:sizeof(float) * npix
                                                            options:MTLResourceStorageModeShared]
                                    : RMBuf(nullptr, 0),
                              rec,
                              RMBuf(vmax.data(), sizeof(float) * vmax.size())};

    auto t0 = std::chrono::steady_clock::now();
    /* bands of rows, each its own command buffer, so no single GPU
       submission runs long enough to trip the display watchdog */
    NSMutableArray* cbs = [NSMutableArray array];
    const int nbands = (p.oy1 - p.oy0 + 127) / 128;
    auto bandsDone = std::make_shared<std::atomic<int>>(0);
    for (int r = p.oy0; r < p.oy1; r += 128) {
      int32_t row0 = r;
      id<MTLCommandBuffer> cb = [s_queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
      [enc setComputePipelineState:s_pipeline];
      [enc setBytes:&pk length:sizeof(pk) atIndex:0];
      for (int k = 1; k < 13; k++) [enc setBuffer:bufs[k] offset:0 atIndex:k];
      [enc setBytes:&row0 length:sizeof(row0) atIndex:13];
      /* 8x8 tiles: a SIMD-group is an 8x4 block of neighbouring rays */
      [enc dispatchThreads:MTLSizeMake(w, std::min(128, p.oy1 - r), 1)
          threadsPerThreadgroup:MTLSizeMake(8, 8, 1)];
      [enc endEncoding];
      /* the trace is most of a large render: advance the progress bar per band */
      [cb addCompletedHandler:^(id<MTLCommandBuffer>) {
        RayProgress(0.18f + 0.78f * (bandsDone->fetch_add(1) + 1) / nbands);
      }];
      [cb commit];
      [cbs addObject:cb];
    }
    id<MTLBuffer> delta = nil, fill = nil, aaOut = nil;
    const size_t aaN = (size_t)(p.width / 3 - 2) * (p.height / 3 - 2);
    if (!outline) aa3 = nullptr;
    if (outline) {
      /* RayRender's outline pass: gradients of the traced depth, then (with
         fog) three depth-fill passes ping-ponging as on the host — the second
         writes over the traced depth — then the edges */
      struct {
        int32_t width, height, ox0, oy0, ox1, oy1, mode, fogFlag, fogRangeFlag, opaque;
        uint32_t trace_word, background;
        float fog, fog_start, inv1minusFogStart, invFrontMinusBack, slope_f, depth_f, disco_f,
            disco_raw;
      } O = {p.width,
             p.height,
             p.ox0,
             p.oy0,
             p.ox1,
             p.oy1,
             outline->mode,
             outline->fogFlag,
             outline->fogRangeFlag,
             outline->opaque,
             outline->trace_word,
             p.background,
             outline->fog,
             outline->fog_start,
             outline->inv1minusFogStart,
             outline->invFrontMinusBack,
             outline->slope_f,
             outline->depth_f,
             outline->disco_f,
             outline->disco_raw};
      delta = [s_device newBufferWithLength:sizeof(float) * 3 * npix
                                    options:MTLResourceStorageModePrivate];
      fill = [s_device newBufferWithLength:sizeof(float) * npix
                                   options:MTLResourceStorageModePrivate];
      id<MTLCommandBuffer> cb = [s_queue commandBuffer];
      id<MTLBlitCommandEncoder> blit = [cb blitCommandEncoder];
      [blit fillBuffer:fill range:NSMakeRange(0, fill.length) value:0];
      [blit endEncoding];
      /* one encoder per pass: each reads the previous one's output */
      auto pass = [&](id<MTLComputePipelineState> ps, std::initializer_list<id<MTLBuffer>> b) {
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:ps];
        [enc setBytes:&O length:sizeof(O) atIndex:0];
        int k = 1;
        for (id<MTLBuffer> x : b) [enc setBuffer:x offset:0 atIndex:k++];
        [enc dispatchThreads:MTLSizeMake(w, h, 1) threadsPerThreadgroup:MTLSizeMake(8, 8, 1)];
        [enc endEncoding];
      };
      id<MTLBuffer> dep = bufs[10], img = bufs[9];
      pass(s_grad, {dep, delta});
      if (outline->fogFlag) {
        pass(s_fill, {dep, fill});
        pass(s_fill, {fill, dep});
        pass(s_fill, {dep, fill});
      }
      pass(s_edges, {delta, outline->fogFlag ? fill : dep, img});
      if (aa3) { /* and RayRender's 3x antialias: only that image comes back */
        aaOut = [s_device newBufferWithLength:sizeof(uint32_t) * aaN
                                      options:MTLResourceStorageModeShared];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:s_aa];
        [enc setBytes:&O length:sizeof(O) atIndex:0];
        [enc setBuffer:img offset:0 atIndex:1];
        [enc setBuffer:aaOut offset:0 atIndex:2];
        [enc dispatchThreads:MTLSizeMake(p.width / 3 - 2, p.height / 3 - 2, 1)
            threadsPerThreadgroup:MTLSizeMake(8, 8, 1)];
        [enc endEncoding];
      }
      [cb commit];
      [cbs addObject:cb];
    }
    bool ok = true;
    for (id<MTLCommandBuffer> cb in cbs) {
      [cb waitUntilCompleted];
      if (ok && cb.status != MTLCommandBufferStatusCompleted) {
        ok = false;
        err = std::string("GPU error: ") +
              (cb.error ? cb.error.localizedDescription.UTF8String : "?");
      }
    }
    *gpu_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    if (ok && aa3) {
      memcpy(aa3, aaOut.contents, sizeof(uint32_t) * aaN);
    } else if (ok) { /* the region back into the full-size image (and depth,
                        if the host still needs it) */
      const uint32_t* gi = (const uint32_t*)bufs[9].contents;
      const float* gd = (const float*)bufs[10].contents;
      for (int y = 0; y < h; y++) {
        size_t o = (size_t)(p.oy0 + y) * p.width + p.ox0;
        memcpy(image + o, gi + (size_t)y * w, sizeof(uint32_t) * w);
        if (depth && !outline) memcpy(depth + o, gd + (size_t)y * w, sizeof(float) * w);
      }
    }
    RM_RELEASE(delta);
    RM_RELEASE(fill);
    RM_RELEASE(aaOut);
    for (int k = 1; k < 13; k++) RM_RELEASE(bufs[k]);
    return ok;
  }
}
