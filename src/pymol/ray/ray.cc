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
-*   Larry Coopet (various optimizations)
-*   Chris Want (RayRenderVRML2, via the public domain )
-*
Z* -------------------------------------------------------------------
*/
#include "pymol/base/busy.h"

#include "pymol/base/err.h"
#include "pymol/ray/ray.h"
#include "pymol/color/color.h"
#include "pymol/geometry/matrix.h"
#include "pymol/glyph/character.h"
#include "pymol/scene/scene.h"
#include "pymol/image/my_png.h"
#include "pymol/base/feedback.h"

#include <thread>
#include <dispatch/dispatch.h>

/* Run fn(y) for y in [0, n) across cores (rows are independent). */
static void RayParallelRows(int n, const std::function<void(int)>& fn) {
  if (n >= 64) {
    struct Ctx {
      const std::function<void(int)>* fn;
      int n, chunk;
    } ctx{&fn, n, 16};
    dispatch_apply_f((n + ctx.chunk - 1) / ctx.chunk,
                     dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                     &ctx, [](void* c, size_t k) {
                       auto* x = (Ctx*)c;
                       int y0 = (int)k * x->chunk,
                           y1 = std::min(x->n, y0 + x->chunk);
                       for (int y = y0; y < y1; y++) (*x->fn)(y);
                     });
    return;
  }
  for (int y = 0; y < n; y++) fn(y);
}
#include "pymol/ray/ray_metal.h"
#include "pymol/ray/ray_internal.h"

#ifndef RAY_SMALL
#endif
typedef float float4[4];

static void RayRelease(CRay* I);

static void RayApplyMatrix33(unsigned int n, float3* q, const float m[16],
                             float3* p);
void RaySetContext(CRay* I, pymol::RenderContext context) {
  I->context = context;
}

float RayGetScreenVertexScale(CRay* I, float* v1) {
  /* what size should a screen pixel be at the coordinate provided? */

  float vt[3];
  float ratio;
  RayApplyMatrix33(1, (float3*)vt, I->ModelView, (float3*)v1);

  if (I->Ortho) {
    ratio = 2 * (float)(fabs(I->Pos.z) * tan((I->Fov / 2.0) * cPI / 180.0)) /
            (I->Height);
  } else {
    float front_size = 2 * I->Volume[4] *
                       ((float)tan((I->Fov / 2.0F) * PI / 180.0F)) /
                       (I->Height);
    ratio = fabs(front_size * (-vt[2] / I->Volume[4]));
  }
  return ratio;
}

int RayGetNPrimitives(CRay* I) { return (I->NPrimitive); }
#define MIN(x, y) ((x < y) ? x : y)
#define MAX(x, y) ((x > y) ? x : y)
int RayExpandPrimitives(CRay* I) {
  int a;
  float *v0, *v1, *n0, *n1;
  CBasis* basis;
  int nVert, nNorm;
  float voxel_floor;
  int ok = true;

  nVert = 0;
  nNorm = 0;
  for (a = 0; a < I->NPrimitive; a++) {
    switch (I->Primitive[a].type) {
      case cPrimSphere:
        nVert++;
        break;
      case cPrimEllipsoid:
        nVert++;
        nNorm += 3;
        break;
      case cPrimCone:
      case cPrimCylinder:
      case cPrimSausage:
        nVert++;
        nNorm++;
        break;
      case cPrimTriangle:
      case cPrimCharacter:
        nVert += 3;
        nNorm += 4;
        break;
    }
  }

  basis = I->Basis;

  VLASize2<float>(basis->Vertex, 3 * nVert);
  VLASize2<float>(basis->Radius, nVert);
  VLASize2<float>(basis->Radius2, nVert);
  VLASize2<int>(basis->Vert2Normal, nVert);
  VLASize2<float>(basis->Normal, 3 * nNorm);

  I->Vert2Prim.resize(nVert);

  voxel_floor = I->PixelRadius / 2.0F;

  basis->MaxRadius = 0.0F;
  basis->MinVoxel = 0.0F;
  basis->NVertex = nVert;
  basis->NNormal = nNorm;

  nVert = 0;
  nNorm = 0;
  v0 = basis->Vertex;
  n0 = basis->Normal;
  ok &= !I->G->Interrupt;
  for (a = 0; ok && a < I->NPrimitive; a++) {
    switch (I->Primitive[a].type) {
      case cPrimTriangle:
      case cPrimCharacter:

        I->Primitive[a].vert = nVert;
        I->Vert2Prim[nVert] = a;
        I->Vert2Prim[nVert + 1] = a;
        I->Vert2Prim[nVert + 2] = a;
        basis->Radius[nVert] = I->Primitive[a].r1;
        basis->Radius2[nVert] =
            I->Primitive[a].r1 * I->Primitive[a].r1; /*necessary?? */
        /*              if(basis->Radius[nVert]>basis->MinVoxel)
           basis->MinVoxel=basis->Radius[nVert]; */
        if (basis->MinVoxel < voxel_floor) basis->MinVoxel = voxel_floor;
        basis->Vert2Normal[nVert] = nNorm;
        basis->Vert2Normal[nVert + 1] = nNorm;
        basis->Vert2Normal[nVert + 2] = nNorm;
        n1 = I->Primitive[a].n0;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        n1 = I->Primitive[a].n1;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        n1 = I->Primitive[a].n2;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        n1 = I->Primitive[a].n3;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        nNorm += 4;
        v1 = I->Primitive[a].v1;
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        v1 = I->Primitive[a].v2;
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        v1 = I->Primitive[a].v3;
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        nVert += 3;
        break;
      case cPrimSphere:
        I->Primitive[a].vert = nVert;
        I->Vert2Prim[nVert] = a;
        v1 = I->Primitive[a].v1;
        basis->Radius[nVert] = I->Primitive[a].r1;
        basis->Radius2[nVert] =
            I->Primitive[a].r1 * I->Primitive[a].r1; /*precompute */
        if (basis->Radius[nVert] > basis->MaxRadius)
          basis->MaxRadius = basis->Radius[nVert];
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        nVert++;
        break;
      case cPrimEllipsoid:
        I->Primitive[a].vert = nVert;
        I->Vert2Prim[nVert] = a;
        v1 = I->Primitive[a].v1;
        basis->Radius[nVert] = I->Primitive[a].r1;
        basis->Radius2[nVert] =
            I->Primitive[a].r1 * I->Primitive[a].r1; /*precompute */
        if (basis->Radius[nVert] > basis->MaxRadius)
          basis->MaxRadius = basis->Radius[nVert];
        basis->Vert2Normal[nVert] = nNorm;
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        nVert++;
        n1 = I->Primitive[a].n1;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        n1 = I->Primitive[a].n2;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        n1 = I->Primitive[a].n3;
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        (*n0++) = (*n1++);
        nNorm += 3;
        break;
      case cPrimCone:
      case cPrimCylinder:
      case cPrimSausage:
        I->Primitive[a].vert = nVert;
        I->Vert2Prim[nVert] = a;
        basis->Radius[nVert] = I->Primitive[a].r1;
        basis->Radius2[nVert] =
            I->Primitive[a].r1 * I->Primitive[a].r1; /*precompute */
        if (basis->MinVoxel < voxel_floor) basis->MinVoxel = voxel_floor;
        subtract3f(I->Primitive[a].v2, I->Primitive[a].v1, n0);
        I->Primitive[a].l1 = (float)length3f(n0);
        normalize3f(n0);
        n0 += 3;
        basis->Vert2Normal[nVert] = nNorm;
        nNorm++;
        v1 = I->Primitive[a].v1;
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        (*v0++) = (*v1++);
        nVert++;
        break;
    }
    ok &= !I->G->Interrupt;
  }
  if (nVert > basis->NVertex) {
    fprintf(stderr, "Error: basis->NVertex exceeded\n");
  }
  PRINTFB(I->G, FB_Ray, FB_Blather)
  " Ray: minvoxel  %8.3f\n Ray: NPrimit  %d nvert %d\n", basis->MinVoxel,
      I->NPrimitive, nVert ENDFB(I->G);
  return ok;
}
void RayComputeBox(CRay* I) {
#define minmax(v, r)          \
  {                           \
    xp = v[0] + r;            \
    xm = v[0] - r;            \
    yp = v[1] + r;            \
    ym = v[1] - r;            \
    zp = v[2] + r;            \
    zm = v[2] - r;            \
    if (xmin > xm) xmin = xm; \
    if (xmax < xp) xmax = xp; \
    if (ymin > ym) ymin = ym; \
    if (ymax < yp) ymax = yp; \
    if (zmin > zm) zmin = zm; \
    if (zmax < zp) zmax = zp; \
  }

  CPrimitive* prm;
  CBasis* basis1;

  float xmin = 0.0F, ymin = 0.0F, xmax = 0.0F, ymax = 0.0F, zmin = 0.0F,
        zmax = 0.0F;
  float xp, xm, yp, ym, zp, zm;

  float *v, r;
  float vt[3];
  const float _0 = 0.0F;
  int a;

  basis1 = I->Basis + 1;
  if (basis1->NVertex) {
    xmin = xmax = basis1->Vertex[0];
    ymin = ymax = basis1->Vertex[1];
    zmin = zmax = basis1->Vertex[2];

    for (a = 0; a < I->NPrimitive; a++) {
      prm = I->Primitive + a;

      switch (prm->type) {
        case cPrimTriangle:
        case cPrimCharacter:

          r = _0;
          v = basis1->Vertex + prm->vert * 3;
          minmax(v, r);
          v = basis1->Vertex + prm->vert * 3 + 3;
          minmax(v, r);
          v = basis1->Vertex + prm->vert * 3 + 6;
          minmax(v, r);
          break;
        case cPrimSphere:
        case cPrimEllipsoid:
          r = prm->r1;
          v = basis1->Vertex + prm->vert * 3;
          minmax(v, r);
          break;
        case cPrimCone:
        case cPrimCylinder:
        case cPrimSausage:
          r = prm->r1;
          v = basis1->Vertex + prm->vert * 3;
          minmax(v, r);
          v = basis1->Normal + basis1->Vert2Normal[prm->vert] * 3;
          scale3f(v, prm->l1, vt);
          v = basis1->Vertex + prm->vert * 3;
          add3f(v, vt, vt);
          minmax(vt, r);
          break;
      } /* end of switch */
    }
  }
  // without the R_SMALL4 padding, this caused ray tracing failure
  // on OS X (X11) with a flat ObjectCGO (zmin == zmax).
  I->min_box[0] = xmin - R_SMALL4;
  I->min_box[1] = ymin - R_SMALL4;
  I->min_box[2] = zmin - R_SMALL4;
  I->max_box[0] = xmax + R_SMALL4;
  I->max_box[1] = ymax + R_SMALL4;
  I->max_box[2] = zmax + R_SMALL4;
}

int RayTransformFirst(CRay* I, int perspective, int identity) {
  CBasis *basis0, *basis1;
  CPrimitive* prm;
  int a;
  float* v0;
  int backface_cull;
  int ok = true;
  int two_sided_lighting =
      SettingGetGlobal_b(I->G, cSetting_two_sided_lighting);

  if (two_sided_lighting < 0) {
    if (SettingGetGlobal_i(I->G, cSetting_surface_cavity_mode))
      two_sided_lighting = true;
    else
      two_sided_lighting = false;
  }

  backface_cull = SettingGetGlobal_b(I->G, cSetting_backface_cull);

  if (two_sided_lighting ||
      (SettingGetGlobal_i(I->G, cSetting_transparency_mode) == 1) ||
      (SettingGetGlobal_i(I->G, cSetting_ray_interior_color) != -1) ||
      I->CheckInterior)
    backface_cull = 0;

  basis0 = I->Basis;
  basis1 = I->Basis + 1;

  if (ok) {
    VLASize2<float>(basis1->Vertex, 3 * basis0->NVertex);
    CHECKOK(ok, basis1->Vertex);
  }
  if (ok) {
    VLASize2<float>(basis1->Normal, 3 * basis0->NNormal);
    CHECKOK(ok, basis1->Normal);
  }
  if (ok) {
    VLASize2<float>(basis1->Precomp, 3 * basis0->NNormal);
    CHECKOK(ok, basis1->Precomp);
  }
  if (ok) {
    VLASize2<int>(basis1->Vert2Normal, basis0->NVertex);
    CHECKOK(ok, basis1->Vert2Normal);
  }
  if (ok) {
    VLASize2<float>(basis1->Radius, basis0->NVertex);
    CHECKOK(ok, basis1->Radius);
  }
  if (ok) {
    VLASize2<float>(basis1->Radius2, basis0->NVertex);
    CHECKOK(ok, basis1->Radius2);
  }
  ok &= !I->G->Interrupt;
  if (ok) {
    if (identity) {
      UtilCopyMem(basis1->Vertex, basis0->Vertex,
                  basis0->NVertex * sizeof(float) * 3);
    } else {
      RayApplyMatrix33(basis0->NVertex, (float3*)basis1->Vertex, I->ModelView,
                       (float3*)basis0->Vertex);
    }
  }
  ok &= !I->G->Interrupt;

  if (ok) {
    memcpy(basis1->Radius, basis0->Radius, basis0->NVertex * sizeof(float));
    memcpy(basis1->Radius2, basis0->Radius2, basis0->NVertex * sizeof(float));
    memcpy(basis1->Vert2Normal, basis0->Vert2Normal,
           basis0->NVertex * sizeof(int));
  }
  ok &= !I->G->Interrupt;
  if (ok) {
    basis1->MaxRadius = basis0->MaxRadius;
    basis1->MinVoxel = basis0->MinVoxel;
    basis1->NVertex = basis0->NVertex;
  }
  ok &= !I->G->Interrupt;
  if (ok) {
    if (identity) {
      UtilCopyMem(basis1->Normal, basis0->Normal,
                  basis0->NNormal * sizeof(float) * 3);
    } else {
      RayTransformNormals33(basis0->NNormal, (float3*)basis1->Normal,
                            I->ModelView, (float3*)basis0->Normal);
    }
    basis1->NNormal = basis0->NNormal;
  }
  ok &= !I->G->Interrupt;
  if (perspective) {
    for (a = 0; ok && a < I->NPrimitive; a++) {
      prm = I->Primitive + a;

      prm = I->Primitive + a;
      switch (prm->type) {
        case cPrimTriangle:
        case cPrimCharacter:
          BasisTrianglePrecomputePerspective(
              basis1->Vertex + prm->vert * 3,
              basis1->Vertex + prm->vert * 3 + 3,
              basis1->Vertex + prm->vert * 3 + 6,
              basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
          break;
      }
      ok &= !I->G->Interrupt;
    }
  } else {
    for (a = 0; ok && a < I->NPrimitive; a++) {
      prm = I->Primitive + a;
      switch (prm->type) {
        case cPrimTriangle:
        case cPrimCharacter:
          BasisTrianglePrecompute(
              basis1->Vertex + prm->vert * 3,
              basis1->Vertex + prm->vert * 3 + 3,
              basis1->Vertex + prm->vert * 3 + 6,
              basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
          v0 = basis1->Normal + (basis1->Vert2Normal[prm->vert] * 3 + 3);
          prm->cull = (!identity) && backface_cull &&
                      ((v0[2] < 0.0F) && (v0[5] < 0.0F) && (v0[8] < 0.0F));
          break;
        case cPrimCone:
        case cPrimSausage:
        case cPrimCylinder:
          BasisCylinderSausagePrecompute(
              basis1->Normal + basis1->Vert2Normal[prm->vert] * 3,
              basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
          break;
      }
      ok &= !I->G->Interrupt;
    }
  }
  return ok;
}
static int RayTransformBasis(CRay* I, CBasis* basis1) {
  CBasis* basis0;
  int a;
  float *v0, *v1;
  CPrimitive* prm;
  int ok = true;

  basis0 = I->Basis + 1;

  VLASize2<float>(basis1->Vertex, 3 * basis0->NVertex);
  CHECKOK(ok, basis1->Vertex);
  if (ok) VLASize2<float>(basis1->Normal, 3 * basis0->NNormal);
  CHECKOK(ok, basis1->Normal);
  if (ok) VLASize2<float>(basis1->Precomp, 3 * basis0->NNormal);
  CHECKOK(ok, basis1->Precomp);
  if (ok) VLASize2<int>(basis1->Vert2Normal, basis0->NVertex);
  CHECKOK(ok, basis1->Vert2Normal);
  if (ok) VLASize2<float>(basis1->Radius, basis0->NVertex);
  CHECKOK(ok, basis1->Radius);
  if (ok) VLASize2<float>(basis1->Radius2, basis0->NVertex);
  CHECKOK(ok, basis1->Radius2);
  v0 = basis0->Vertex;
  v1 = basis1->Vertex;
  for (a = 0; ok && a < basis0->NVertex; a++) {
    matrix_transform33f3f(basis1->Matrix, v0, v1);
    v0 += 3;
    v1 += 3;
    basis1->Radius[a] = basis0->Radius[a];
    basis1->Radius2[a] = basis0->Radius2[a];
    basis1->Vert2Normal[a] = basis0->Vert2Normal[a];
  }
  if (ok) {
    v0 = basis0->Normal;
    v1 = basis1->Normal;
  }
  for (a = 0; ok && a < basis0->NNormal; a++) {
    matrix_transform33f3f(basis1->Matrix, v0, v1);
    v0 += 3;
    v1 += 3;
  }
  if (ok) {
    basis1->MaxRadius = basis0->MaxRadius;
    basis1->MinVoxel = basis0->MinVoxel;
    basis1->NVertex = basis0->NVertex;
    basis1->NNormal = basis0->NNormal;
  }
  for (a = 0; ok && a < I->NPrimitive; a++) {
    prm = I->Primitive + a;
    switch (prm->type) {
      case cPrimTriangle:
      case cPrimCharacter:

        BasisTrianglePrecompute(
            basis1->Vertex + prm->vert * 3, basis1->Vertex + prm->vert * 3 + 3,
            basis1->Vertex + prm->vert * 3 + 6,
            basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
        break;
      case cPrimCone:
      case cPrimSausage:
      case cPrimCylinder:
        BasisCylinderSausagePrecompute(
            basis1->Normal + basis1->Vert2Normal[prm->vert] * 3,
            basis1->Precomp + basis1->Vert2Normal[prm->vert] * 3);
        break;
    }
  }
  return ok;
}
void RayRenderTest(CRay* I, int width, int height, float front, float back,
                   float fov) {
  PRINTFB(I->G, FB_Ray, FB_Details)
  " RayRenderTest: obtained %i graphics primitives.\n",
      I->NPrimitive ENDFB(I->G);
}

#define noIDTF_COLOR
void RayProjectTriangle(CRay* I, RayInfo* r, float* light, float* v0, float* n0,
                        float scale) {
  float w2;
  float d1[3], d2[3], d3[3];
  float p1[3], p2[3], p3[3];
  int c = 0;
  const float _0 = 0.0F;
  float* impact = r->impact;

  if (dot_product3f(light, n0 - 3) >= _0)
    c++;
  else if (dot_product3f(light, n0) >= _0)
    c++;
  else if (dot_product3f(light, n0 + 3) >= _0)
    c++;
  else if (dot_product3f(light, n0 + 6) >= _0)
    c++;

  if (c) {
    w2 = 1.0F - (r->tri1 + r->tri2);

    subtract3f(v0, impact, d1);
    subtract3f(v0 + 3, impact, d2);
    subtract3f(v0 + 6, impact, d3);
    project3f(d1, n0, p1);
    project3f(d2, n0 + 3, p2);
    project3f(d3, n0 + 6, p3);
    scale3f(p1, w2, d1);
    scale3f(p2, r->tri1, d2);
    scale3f(p3, r->tri2, d3);
    add3f(d1, d2, d2);
    add3f(d2, d3, d3);
    scale3f(d3, scale, d3);
    if (dot_product3f(r->surfnormal, d3) >= _0) add3f(d3, impact, impact);
  }
}

/* This is both an antialias and a slight blur */

/* for whatever reason, greatly GCC perfers a linear sequence of
   accumulates over a single large expression -- the difference is
   huge: over 10% !!! */

#define combine4by4(var, src, mask)   \
  {                                   \
    var = ((src)[0] & mask);          \
    var += ((src)[1] & mask);         \
    var += ((src)[2] & mask);         \
    var += ((src)[3] & mask);         \
    var += ((src)[4] & mask);         \
    var += (((src)[5] & mask) * 13);  \
    var += (((src)[6] & mask) * 13);  \
    var += ((src)[7] & mask);         \
    var += ((src)[8] & mask);         \
    var += (((src)[9] & mask) * 13);  \
    var += (((src)[10] & mask) * 13); \
    var += ((src)[11] & mask);        \
    var += ((src)[12] & mask);        \
    var += ((src)[13] & mask);        \
    var += ((src)[14] & mask);        \
    var += ((src)[15] & mask);        \
    var = (var >> 6) & mask;          \
  }

#define combine5by5(var, src, mask)  \
  {                                  \
    var = ((src)[0] & mask);         \
    var += ((src)[1] & mask);        \
    var += ((src)[2] & mask);        \
    var += ((src)[3] & mask);        \
    var += ((src)[4] & mask);        \
    var += ((src)[5] & mask);        \
    var += (((src)[6] & mask) * 5);  \
    var += (((src)[7] & mask) * 5);  \
    var += (((src)[8] & mask) * 5);  \
    var += ((src)[9] & mask);        \
    var += ((src)[10] & mask);       \
    var += (((src)[11] & mask) * 5); \
    var += (((src)[12] & mask) * 8); \
    var += (((src)[13] & mask) * 5); \
    var += ((src)[14] & mask);       \
    var += ((src)[15] & mask);       \
    var += (((src)[16] & mask) * 5); \
    var += (((src)[17] & mask) * 5); \
    var += (((src)[18] & mask) * 5); \
    var += ((src)[19] & mask);       \
    var += ((src)[20] & mask);       \
    var += ((src)[21] & mask);       \
    var += ((src)[22] & mask);       \
    var += ((src)[23] & mask);       \
    var += ((src)[24] & mask);       \
    var = (var >> 6) & mask;         \
  }

#define combine6by6(var, src, mask)  \
  {                                  \
    var = ((src)[0] & mask);         \
    var += ((src)[1] & mask);        \
    var += ((src)[2] & mask);        \
    var += ((src)[3] & mask);        \
    var += ((src)[4] & mask);        \
    var += ((src)[5] & mask);        \
    var += ((src)[6] & mask);        \
    var += (((src)[7] & mask) * 5);  \
    var += (((src)[8] & mask) * 7);  \
    var += (((src)[9] & mask) * 7);  \
    var += (((src)[10] & mask) * 5); \
    var += ((src)[11] & mask);       \
    var += ((src)[12] & mask);       \
    var += (((src)[13] & mask) * 7); \
    var += (((src)[14] & mask) * 8); \
    var += (((src)[15] & mask) * 8); \
    var += (((src)[16] & mask) * 7); \
    var += ((src)[17] & mask);       \
    var += ((src)[18] & mask);       \
    var += (((src)[19] & mask) * 7); \
    var += (((src)[20] & mask) * 8); \
    var += (((src)[21] & mask) * 8); \
    var += (((src)[22] & mask) * 7); \
    var += ((src)[23] & mask);       \
    var += ((src)[24] & mask);       \
    var += (((src)[25] & mask) * 5); \
    var += (((src)[26] & mask) * 7); \
    var += (((src)[27] & mask) * 7); \
    var += (((src)[28] & mask) * 5); \
    var += ((src)[29] & mask);       \
    var += ((src)[30] & mask);       \
    var += ((src)[31] & mask);       \
    var += ((src)[32] & mask);       \
    var += ((src)[33] & mask);       \
    var += ((src)[34] & mask);       \
    var += ((src)[35] & mask);       \
    var = (var >> 7) & mask;         \
  }

#define m00FF 0x00FF
#define mFF00 0xFF00
#define mFFFF 0xFFFF

float* rayDepthPixels = nullptr;
int rayVolume = 0, rayWidth = 0, rayHeight = 0;
/*
 * RayTraceThread on the GPU (src/pymol/ray/ray_metal.mm): the same
 * per-pixel setup, flattened for the kernel. Returns false — and the caller
 * runs RayTraceThread — when the scene uses anything the kernel doesn't port
 * or the GPU can't run it.
 */
/* Why the GPU tracer can't take this scene (nullptr: it can). Settings and
   primitives only — decided before the voxel map is built. */
static const char* RayMetalUnsupported(CRay* I, int perspective) {
  PyMOLGlobals* G = I->G;
  if (perspective) return "perspective (set orthoscopic on)";
  if (SettingGetGlobal_i(G, cSetting_ray_shadows) && I->NBasis > 2)
    return "ray_shadows";
  if (SettingGetGlobal_i(G, cSetting_ray_interior_color) != -1 ||
      I->CheckInterior)
    return "interior colour";
  if (SettingGetGlobal_i(G, cSetting_chromadepth) > 0) return "chromadepth";
  if (I->BigEndian) return "unsupported mode";
  if (!I->NPrimitive) return "empty scene";
  if (I->NBasis - 2 > 10) return "too many lights";
  for (int a = 0; a < I->NPrimitive; a++) {
    const CPrimitive* p = I->Primitive + a;
    switch (p->type) {
      case cPrimSphere:
      case cPrimCylinder:
      case cPrimTriangle:
      case cPrimSausage:
        break;
      case cPrimCharacter:
        return "labels";
      default:
        return "ellipsoids/cones";
    }
    if (p->wobble) return "ray_texture";
    /* only the colour slots the type shades with: spheres leave c2/c3 unset */
    if (p->ramped || p->c1[0] <= cColorExtCutoff ||
        (p->type != cPrimSphere && p->c2[0] <= cColorExtCutoff) ||
        (p->type == cPrimTriangle && p->c3[0] <= cColorExtCutoff))
      return "ramped colours";
  }
  return nullptr;
}

static bool RayTraceMetal(CRayThreadInfo* T, unsigned int trace_word,
                          int ray_trace_mode, bool* outlined,
                          unsigned int* aa_dst, bool* antialiased) {
  CRay* I = T->ray;
  PyMOLGlobals* G = I->G;
  char msg[256];
  auto fallback = [&](const char* why) {
    snprintf(msg, sizeof(msg), "cpu: %s", why);
    RayMetalSetLastReport(msg);
    return false;
  };

  if (const char* why = RayMetalUnsupported(I, T->perspective))
    return fallback(why);
  if (T->bkrd_data || T->bkrd_is_gradient)
    return fallback("gradient/image background");
  if (T->edging) return fallback("unsupported mode");
  CBasis* bp1 = I->Basis + 1;
  if (!bp1->Map) return fallback("empty scene");

  std::vector<RayMetalPrim> prims(I->NPrimitive);
  for (int a = 0; a < I->NPrimitive; a++) {
    const CPrimitive* p = I->Primitive + a;
    RayMetalPrim& g = prims[a];
    copy3f(p->c1, g.c1);
    copy3f(p->c2, g.c2);
    copy3f(p->c3, g.c3);
    copy3f(p->tr, g.tr);
    g.r1 = p->r1;
    g.l1 = p->l1;
    g.trans = p->trans;
    g.vert = p->vert;
    g.type = p->type;
    g.cap1 = (int)p->cap1;
    g.cap2 = (int)p->cap2;
    g.cull = p->cull;
    g.no_lighting = p->no_lighting;
  }

  RayMetalParams P;
  memset(&P, 0, sizeof(P));
  const float _0 = 0.0F, _1 = 1.0F;

  /* --- the scalar setup of RayTraceThread, in the same order --- */
  float ray_scatter = SettingGetGlobal_f(G, cSetting_ray_scatter);
  float legacy = SettingGetGlobal_f(G, cSetting_ray_legacy_lighting);
  int spec_count = SettingGetGlobal_i(G, cSetting_spec_count);
  {
    float fudge = SettingGetGlobal_f(G, cSetting_ray_triangle_fudge);
    P.fudge0 = 0.0F - fudge;
    P.fudge1 = 1.0F + fudge;
  }
  if (spec_count < 0) spec_count = SettingGetGlobal_i(G, cSetting_light_count);
  int trans_mode = SettingGetGlobal_i(G, cSetting_transparency_mode);
  int backface_cull = SettingGetGlobal_i(G, cSetting_backface_cull);
  int opaque_back = SettingGetGlobal_i(G, cSetting_ray_opaque_background);
  if (opaque_back < 0)
    opaque_back = SettingGetGlobal_i(G, cSetting_opaque_background);
  int two_sided_lighting = SettingGetGlobal_i(G, cSetting_two_sided_lighting);
  if (two_sided_lighting < 0)
    two_sided_lighting =
        SettingGetGlobal_i(G, cSetting_surface_cavity_mode) ? true : false;
  float trans_cont = SettingGetGlobal_f(G, cSetting_ray_transparency_contrast);
  if (trans_mode == 1) two_sided_lighting = true;
  P.trans_cont = trans_cont;
  P.inv_trans_cont = 1.0F;
  if (trans_cont > 1.0F) {
    P.trans_cont_flag = true;
    P.inv_trans_cont = 1.0F / trans_cont;
  }
  float ambient = T->ambient, lreflect, direct;
  {
    float reflect_scale = SceneGetReflectScaleValue(G, 10);
    lreflect =
        reflect_scale * (SettingGetGlobal_f(G, cSetting_reflect) - ray_scatter);
    if (lreflect < _0) lreflect = _0;
    ray_scatter = ray_scatter * reflect_scale;
  }
  direct = SettingGetGlobal_f(G, cSetting_direct);
  ambient *= (1.0F - legacy) + (legacy * (0.22F / 0.12F));
  lreflect *= (1.0F - legacy) + (legacy * (0.72F / 0.45F));
  direct *= (1.0F - legacy) + (legacy * (0.24F / 0.45F));
  P.ambient = ambient;
  P.lreflect = lreflect;
  P.direct = direct;
  P.ray_scatter = ray_scatter;
  P.legacy = legacy;
  P.legacy_1m = _1 - legacy;
  P.direct_shade = SettingGetGlobal_f(G, cSetting_ray_direct_shade);
  float trans_spec_cut =
      SettingGetGlobal_f(G, cSetting_ray_transparency_spec_cut);
  P.trans_spec_cut = trans_spec_cut;
  P.trans_spec_scale = (trans_spec_cut < _1) ? _1 / (_1 - trans_spec_cut) : _0;
  P.blend_colors = SettingGetGlobal_i(G, cSetting_ray_blend_colors);
  if (P.blend_colors) {
    P.red_blend = SettingGetGlobal_f(G, cSetting_ray_blend_red);
    P.green_blend = SettingGetGlobal_f(G, cSetting_ray_blend_green);
    P.blue_blend = SettingGetGlobal_f(G, cSetting_ray_blend_blue);
  }
  P.max_pass = SettingGetGlobal_i(G, cSetting_ray_max_passes);
  P.ray_trans_spec = SettingGetGlobal_f(G, cSetting_ray_transparency_specular);
  P.trans_oblique = SettingGetGlobal_f(G, cSetting_ray_transparency_oblique);
  P.oblique_power =
      SettingGetGlobal_f(G, cSetting_ray_transparency_oblique_power);
  P.trans_cutoff = SettingGetGlobal_f(G, cSetting_ray_trace_trans_cutoff);
  P.persist_cutoff = SettingGetGlobal_f(G, cSetting_ray_trace_persist_cutoff);
  P.trans_mode = trans_mode;
  P.settingPower = SettingGetGlobal_f(G, cSetting_power);
  P.settingReflectPower = SettingGetGlobal_f(G, cSetting_reflect_power);
  SceneGetAdjustedLightValues(G, &P.settingSpecReflect, &P.settingSpecPower,
                              &P.settingSpecDirect, &P.settingSpecDirectPower,
                              10);
  if (two_sided_lighting || (trans_mode == 1)) backface_cull = 0;
  P.backface_cull = backface_cull;
  P.two_sided_lighting = two_sided_lighting;
  P.spec_count = spec_count;

  P.inv1minusFogStart = _1;
  float fog = SettingGetGlobal_f(G, cSetting_ray_trace_fog);
  if (fog < 0.0F)
    fog = SettingGetGlobal_b(G, cSetting_depth_cue)
              ? SettingGetGlobal_f(G, cSetting_fog)
              : _0;
  if (fog != _0) {
    P.fogFlag = true;
    float fog_start = SettingGetGlobal_f(G, cSetting_ray_trace_fog_start);
    if (fog_start < 0.0F) fog_start = SettingGetGlobal_f(G, cSetting_fog_start);
    if (fog_start > 1.0F) fog_start = 1.0F;
    if (fabs(fog_start) > R_SMALL4) {
      P.fogRangeFlag = true;
      if (fabs(fog_start - 1.0F) < R_SMALL4) P.fogFlag = false;
    }
    P.inv1minusFogStart = _1 / (_1 - fog_start);
    P.fog_start = fog_start;
  }
  P.fog = fog;

  float invHgt, invWdth;
  if (T->border) {
    invHgt = _1 / (float)(T->height - (3.0F + T->border));
    invWdth = _1 / (float)(T->width - (3.0F + T->border));
  } else {
    invHgt = _1 / (float)(T->height);
    invWdth = _1 / (float)(T->width);
  }
  P.invWdthRange = invWdth * I->Range[0];
  P.invHgtRange = invHgt * I->Range[1];
  P.vol0 = I->Volume[0];
  P.vol2 = I->Volume[2];
  P.invFrontMinusBack = _1 / (T->front - T->back);
  P.border_offset = T->border ? (-1.50F + T->border / 2.0F) : 0.0F;
  P.project_triangle = SettingGetGlobal_f(G, cSetting_ray_improve_shadows);

  P.n_basis = I->NBasis;
  P.has_bp2 = I->NBasis > 2;
  for (int bc = 2; bc < I->NBasis; bc++) {
    copy3f(I->Basis[bc].LightNormal, P.lightN + (bc - 2) * 3);
    copy3f(I->Basis[bc].SpecNormal, P.specN + (bc - 2) * 3);
  }

  P.width = T->width;
  P.height = T->height;
  P.x_start = T->x_start;
  P.x_stop = T->x_stop;
  P.y_start = T->y_start;
  P.y_stop = T->y_stop;
  P.front = T->front;
  P.back = T->back;
  P.orig_opaque_back = opaque_back;
  P.fore_mask = T->fore_mask;
  P.bkrd0 = T->bkrd_top[0];
  P.bkrd1 = T->bkrd_top[1];
  P.bkrd2 = T->bkrd_top[2];
  P.write_depth = T->depth != nullptr;
  P.background = T->background;

  /* RayRender's outline pass runs on the GPU too (not for ray_volume, which
     reads the depth back): the GPU then writes the area it can change */
  RayMetalOutline O;
  *outlined =
      ray_trace_mode && T->depth && !SettingGetGlobal_b(G, cSetting_ray_volume);
  P.ox0 = T->x_start;
  P.ox1 = T->x_stop;
  P.oy0 = T->y_start;
  P.oy1 = T->y_stop;
  if (*outlined) {
    const int m = 6; /* see the margin in RayRender */
    bool rect = ray_trace_mode == 1;
    P.ox0 = rect ? std::max(0, T->x_start - m) : 0;
    P.ox1 = rect ? std::min(T->width, T->x_stop + m) : T->width;
    P.oy0 = rect ? std::max(0, T->y_start - m) : 0;
    P.oy1 = rect ? std::min(T->height, T->y_stop + m) : T->height;
    float gain =
        I->PixelRadius /
        (SettingGetGlobal_f(G, cSetting_ray_trace_gain) * I->Magnified);
    gain /= T->border + 1; /* RayRender's antialias */
    O.mode = ray_trace_mode;
    O.trace_word = trace_word;
    O.opaque = opaque_back;
    O.fogFlag = P.fogFlag;
    O.fogRangeFlag = P.fogRangeFlag;
    O.fog = std::min(P.fog, 1.0F);
    O.fog_start = P.fog_start;
    O.inv1minusFogStart = P.inv1minusFogStart;
    O.invFrontMinusBack = P.invFrontMinusBack;
    O.slope_f = SettingGetGlobal_f(G, cSetting_ray_trace_slope_factor) * gain;
    O.depth_f = SettingGetGlobal_f(G, cSetting_ray_trace_depth_factor) * gain;
    O.disco_raw = SettingGetGlobal_f(G, cSetting_ray_trace_disco_factor);
    O.disco_f = O.disco_raw * gain;
  }

  MapType* map = bp1->Map;
  P.n_vert = bp1->NVertex;
  P.n_eElem = map->size();
  P.dim0 = map->Dim[0];
  P.dim1 = map->Dim[1];
  P.dim2 = map->Dim[2];
  P.d1d2 = map->D1D2;
  P.iMin0 = map->iMin[0];
  P.iMin1 = map->iMin[1];
  P.iMin2 = map->iMin[2];
  P.iMax0 = map->iMax[0];
  P.iMax1 = map->iMax[1];
  P.iMax2 = map->iMax[2];
  P.mapMin0 = map->Min[0];
  P.mapMin1 = map->Min[1];
  P.mapMin2 = map->Min[2];
  P.recipDiv = map->recipDiv;

  RayMetalScene S;
  S.p = P;
  S.prims = prims.data();
  S.n_prims = I->NPrimitive;
  S.vertex = bp1->Vertex;
  S.n_vertex_f = VLAGetSize(bp1->Vertex);
  S.normal = bp1->Normal;
  S.n_normal_f = VLAGetSize(bp1->Normal);
  S.precomp = bp1->Precomp;
  S.radius = bp1->Radius;
  S.radius2 = bp1->Radius2;
  S.vert2normal = bp1->Vert2Normal;
  S.n_vert2normal = VLAGetSize(bp1->Vert2Normal);
  S.vert2prim = I->Vert2Prim.data();
  S.n_vert2prim = I->Vert2Prim.size();
  S.ehead = map->EHead.data();
  S.elist = map->EList.data();
  S.emask = map->EMask.data();
  S.lite = bp1->Lite;
  S.mhead = map->Head.data();
  S.mlink = map->Link.data();
  S.literef = bp1->LiteRef;

  std::string err;
  double ms = 0.0;
  /* and the 3x antialias (antialias 2 with outlines), straight into aa_dst */
  unsigned int* aa3 =
      (*outlined && aa_dst && T->border == 2) ? aa_dst : nullptr;
  *antialiased = false;
  if (!RayMetalTrace(S, T->image, T->depth, *outlined ? &O : nullptr, aa3, err,
                     &ms))
    return fallback(err.c_str());
  snprintf(msg, sizeof(msg), "gpu: %dx%d, %d primitives, %.1f ms", T->width,
           T->height, I->NPrimitive, ms);
  RayMetalSetLastReport(msg);
  *antialiased = aa3 != nullptr;
  return true;
}

void RayRender(CRay* I, unsigned int* image, double timing, float angle,
               int antialias, unsigned int* return_bg) {
  /* the traced rectangle: depth is zero outside it, so the outline pass
     needs only this area plus a margin (set once the rectangle is known) */
  int trace_x0 = 0, trace_x1 = -1, trace_y0 = 0, trace_y1 = -1;
  bool outlined_on_gpu = false; /* RayTraceMetal did the outline pass */
  bool aa_on_gpu = false;       /* ... and the antialias, into image_copy */
  /* per-stage wall times (ms), reported by _mmol_render.gpu_ray_stages() */
  double rm_t[8] = {0};
  int rm_n = 0;
  auto rm_stamp = [&]() {
    if (rm_n < 8) rm_t[rm_n++] = UtilGetSeconds(I->G);
  };
  rm_stamp();
#define RM_STAMP() rm_stamp()
  int a, x, y;
  unsigned int* image_copy = nullptr;
  unsigned int back_mask, fore_mask = 0, trace_word = 0;
  unsigned int background;
  size_t buffer_size;
  int orig_opaque_back = 0, opaque_back = 0;
  int n_hit = 0;
  const float* bkrd_ptr;
  float bkrd_top[3], bkrd_bottom[3];
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  double now;
  int shadows;
  int n_thread;
  int mag = 1;
  int oversample_cutoff;
  int perspective = SettingGetGlobal_i(I->G, cSetting_ray_orthoscopic);
  int n_light = SettingGetGlobal_i(I->G, cSetting_light_count);
  float ambient;
  float* depth = nullptr;
  float front = I->Volume[4];
  float back = I->Volume[5];
  float fov = I->Fov;
  float* pos = glm::value_ptr(I->Pos);
  size_t width = I->Width;
  size_t height = I->Height;
  int ray_trace_mode;
  const float _0 = 0.0F, _p499 = 0.499F;
  int volume;
  const char* bg_image_filename;
  int ok = true;

  if (n_light > 10) n_light = 10;

  if (perspective < 0) perspective = SettingGetGlobal_b(I->G, cSetting_ortho);
  perspective = !perspective;

  VLASize2<CPrimitive>(I->Primitive, I->NPrimitive);

  n_thread = SettingGetGlobal_i(I->G, cSetting_max_threads);
  if (n_thread < 1) n_thread = 1;
  if (n_thread > PYMOL_MAX_THREADS) n_thread = PYMOL_MAX_THREADS;
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if (opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

  orig_opaque_back = opaque_back;
  ray_trace_mode = SettingGetGlobal_i(I->G, cSetting_ray_trace_mode);

  shadows = SettingGetGlobal_i(I->G, cSetting_ray_shadows);

  oversample_cutoff = SettingGetGlobal_i(I->G, cSetting_ray_oversample_cutoff);

  if (antialias < 0) {
    antialias = SettingGetGlobal_i(I->G, cSetting_antialias);
  }

  if (ray_trace_mode && (antialias == 1))
    antialias = 2;
  else if (ray_trace_mode && antialias)
    antialias++;

  if (antialias < 0) antialias = 0;
  if (antialias > 4) antialias = 4;

  if ((!antialias) || ray_trace_mode) oversample_cutoff = 0;

  mag = antialias;
  if (mag < 1) mag = 1;

  if (antialias > 1) {
    width = (width + 2) * mag;
    height = (height + 2) * mag;
    image_copy = image;
    buffer_size = width * height;
    image = pymol::malloc<unsigned int>(buffer_size);
    ErrChkPtr(I->G, image);
  } else {
    buffer_size = width * height;
  }
  if (ray_trace_mode) {
    depth = pymol::calloc<float>(width * height);
  } else if (oversample_cutoff) {
    depth = pymol::calloc<float>(width * height);
  }
  ambient = SettingGetGlobal_f(I->G, cSetting_ambient);

  bkrd_is_gradient = SettingGetGlobal_b(I->G, cSetting_bg_gradient);
  bg_image_filename =
      SettingGet_s(I->G, nullptr, nullptr, cSetting_bg_image_filename);
  I->bkgrd_data = OrthoBackgroundDataGet(*I->G->Ortho);
  if (!I->bkgrd_data && bg_image_filename && bg_image_filename[0]) {
    I->bkgrd_data = MyPNGRead(bg_image_filename);
    bkrd_is_gradient = 0;
  }
  if (I->bkgrd_data) {
    opaque_back = 1;
  }

  if (!opaque_back) {
    bkrd_is_gradient = 0;
  }
  if (bkrd_is_gradient) {
    bkrd_ptr = ColorGet(
        I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb_top));
    copy3f(bkrd_ptr, bkrd_top);
    bkrd_ptr = ColorGet(
        I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb_bottom));
    copy3f(bkrd_ptr, bkrd_bottom);
  } else {
    bkrd_ptr = ColorGet(
        I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
    copy3f(bkrd_ptr, bkrd_top);
    copy3f(bkrd_ptr, bkrd_bottom);
  }
  { /* adjust bkrd and trace to offset the effect of gamma correction */
    float gamma = SettingGetGlobal_f(I->G, cSetting_gamma);
    float inp;
    float sig;
    inp = (bkrd_top[0] + bkrd_top[1] + bkrd_top[2]) / 3.0F;
    if (inp < R_SMALL4)
      sig = 1.0F;
    else
      sig = (float)(pow(inp, gamma)) / inp;
    bkrd_top[0] *= sig;
    bkrd_top[1] *= sig;
    bkrd_top[2] *= sig;
    if (bkrd_top[0] > 1.0F) bkrd_top[0] = 1.0F;
    if (bkrd_top[1] > 1.0F) bkrd_top[1] = 1.0F;
    if (bkrd_top[2] > 1.0F) bkrd_top[2] = 1.0F;

    if (bkrd_is_gradient) {
      float inp;
      float sig;
      inp = (bkrd_bottom[0] + bkrd_bottom[1] + bkrd_bottom[2]) / 3.0F;
      if (inp < R_SMALL4)
        sig = 1.0F;
      else
        sig = (float)(pow(inp, gamma)) / inp;
      bkrd_bottom[0] *= sig;
      bkrd_bottom[1] *= sig;
      bkrd_bottom[2] *= sig;
      if (bkrd_bottom[0] > 1.0F) bkrd_bottom[0] = 1.0F;
      if (bkrd_bottom[1] > 1.0F) bkrd_bottom[1] = 1.0F;
      if (bkrd_bottom[2] > 1.0F) bkrd_bottom[2] = 1.0F;
    } else {
      copy3f(bkrd_top, bkrd_bottom);
    }
    if (ray_trace_mode) {
      float inp;
      float sig;
      int trace_color = SettingGetGlobal_color(I->G, cSetting_ray_trace_color);
      float trgb[3];
      const float* trgb_v = ColorGet(I->G, trace_color);

      copy3f(trgb_v, trgb);
      inp = (trgb[0] + trgb[1] + trgb[2]) / 3.0F;
      if (inp < R_SMALL4)
        sig = 1.0F;
      else
        sig = (float)(pow(inp, gamma)) / inp;
      trgb[0] *= sig;
      trgb[1] *= sig;
      trgb[2] *= sig;
      if (trgb[0] > 1.0F) trgb[0] = 1.0F;
      if (trgb[1] > 1.0F) trgb[1] = 1.0F;
      if (trgb[2] > 1.0F) trgb[2] = 1.0F;

      if (I->BigEndian) {
        trace_word = ((0xFF & ((unsigned int)(trgb[0] * 255 + _p499))) << 24) |
                     ((0xFF & ((unsigned int)(trgb[1] * 255 + _p499))) << 16) |
                     ((0xFF & ((unsigned int)(trgb[2] * 255 + _p499))) << 8) |
                     0xFF;
      } else {
        trace_word = 0xFF000000 |
                     ((0xFF & ((unsigned int)(trgb[2] * 255 + _p499))) << 16) |
                     ((0xFF & ((unsigned int)(trgb[1] * 255 + _p499))) << 8) |
                     ((0xFF & ((unsigned int)(trgb[0] * 255 + _p499))));
      }
    }
  }
  if (opaque_back) {
    if (I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
    fore_mask = back_mask;
  } else {
    back_mask = 0x00000000;
  }
  if (!bkrd_is_gradient) {
    if (I->BigEndian) {
      background =
          back_mask |
          ((0xFF & ((unsigned int)(bkrd_top[0] * 255 + _p499))) << 24) |
          ((0xFF & ((unsigned int)(bkrd_top[1] * 255 + _p499))) << 16) |
          ((0xFF & ((unsigned int)(bkrd_top[2] * 255 + _p499))) << 8);
    } else {
      background =
          back_mask |
          ((0xFF & ((unsigned int)(bkrd_top[2] * 255 + _p499))) << 16) |
          ((0xFF & ((unsigned int)(bkrd_top[1] * 255 + _p499))) << 8) |
          ((0xFF & ((unsigned int)(bkrd_top[0] * 255 + _p499))));
    }
  } else {
    background = back_mask;
  }
  OrthoBusyFast(I->G, 2, 20);

  if (!bkrd_is_gradient) {
    PRINTFB(I->G, FB_Ray, FB_Blather)
    " RayNew: Background = %x %d %d %d\n", background, (int)(bkrd_top[0] * 255),
        (int)(bkrd_top[1] * 255), (int)(bkrd_top[2] * 255) ENDFB(I->G);
  }
  if (return_bg) *return_bg = background;

  if (!I->NPrimitive) { /* nothing to render! */
    if (I->bkgrd_data) {
      fill_background_image(I, image, width, height,
                            width * (unsigned int)height);
    } else if (bkrd_is_gradient) {
      fill_gradient(I, opaque_back, image, bkrd_top, bkrd_bottom, width, height,
                    width * (unsigned int)height);
    } else {
      fill(image, background, width * (unsigned int)height);
    }
  } else {
    if (I->PrimSizeCnt) {
      float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_camera);
      I->PrimSize = I->PrimSize / (I->PrimSizeCnt * factor);
      /*      printf("avg dist %8.7f\n",I->PrimSize); */
    } else {
      I->PrimSize = 0.0F;
    }
    ok &= !I->G->Interrupt;
    if (ok) ok &= RayExpandPrimitives(I);
    if (ok) ok &= RayTransformFirst(I, perspective, false);
    RM_STAMP(); /* 1: expand */
    RayProgress(0.1f);

    OrthoBusyFast(I->G, 3, 20);

    now = UtilGetSeconds(I->G) - timing;

    PRINTFB(I->G, FB_Ray, FB_Blather)
    " Ray: processed %i graphics primitives in %4.2f sec.\n", I->NPrimitive,
        now ENDFB(I->G);

    if (ok) { /* light sources */
      int bc;
      I->NBasis = n_light + 1;
      if (I->NBasis > MAX_BASIS) I->NBasis = MAX_BASIS;
      if (I->NBasis < 2) I->NBasis = 2;
      for (bc = 2; ok && bc < I->NBasis; bc++) {
        ok &= BasisInit(I->G, I->Basis + bc);
      }
      for (bc = 2; ok && bc < I->NBasis; bc++) {
        { /* setup light & rotate if necessary  */
          float light[3];
          const float* lightv;
          switch (bc) {
            default:
            case 2:
              lightv = SettingGetfv(I->G, cSetting_light);
              break;
            case 3:
              lightv = SettingGetfv(I->G, cSetting_light2);
              break;
            case 4:
              lightv = SettingGetfv(I->G, cSetting_light3);
              break;
            case 5:
              lightv = SettingGetfv(I->G, cSetting_light4);
              break;
            case 6:
              lightv = SettingGetfv(I->G, cSetting_light5);
              break;
            case 7:
              lightv = SettingGetfv(I->G, cSetting_light6);
              break;
            case 8:
              lightv = SettingGetfv(I->G, cSetting_light7);
              break;
            case 9:
              lightv = SettingGetfv(I->G, cSetting_light8);
              break;
            case 10:
              lightv = SettingGetfv(I->G, cSetting_light9);
              break;
          }
          copy3f(lightv, light);
          normalize3f(light);

          if (angle) {
            float temp[16];
            identity44f(temp);
            MatrixRotateC44f(temp, (float)-PI * angle / 180, 0.0F, 1.0F, 0.0F);
            MatrixTransformC44fAs33f3f(temp, light, light);
          }

          I->Basis[bc].LightNormal[0] = light[0];
          I->Basis[bc].LightNormal[1] = light[1];
          I->Basis[bc].LightNormal[2] = light[2];
          normalize3f(I->Basis[bc].LightNormal);

          {
            float spec_vector[3];
            copy3f(I->Basis[bc].LightNormal, spec_vector);
            spec_vector[2]--;
            normalize3f(spec_vector);
            copy3f(spec_vector, I->Basis[bc].SpecNormal);
          }
        }

        if (ok && shadows) { /* don't waste time on shadows unless needed */
          BasisSetupMatrix(I->Basis + bc);
          ok &= RayTransformBasis(I, I->Basis + bc);
        }
        ok &= !I->G->Interrupt;
      }
    }

    OrthoBusyFast(I->G, 4, 20);
    if (shadows && (n_thread > 1)) { /* parallel execution */

      CRayHashThreadInfo* thread_info =
          pymol::calloc<CRayHashThreadInfo>(I->NBasis);

      /* rendering map */
      int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();

      thread_info[0].basis = I->Basis + 1;
      thread_info[0].vert2prim = vert2prim_ptr;
      thread_info[0].prim = I->Primitive;
      thread_info[0].n_prim = I->NPrimitive;
      thread_info[0].clipBox = I->Volume;
      thread_info[0].phase = 0;
      thread_info[0].perspective = perspective;
      thread_info[0].front = front;

      thread_info[0].image = image;
      thread_info[0].bkrd_is_gradient = bkrd_is_gradient;
      thread_info[0].width = width;
      thread_info[0].height = height;
      if (I->bkgrd_data) {
        thread_info[0].background = background;
        thread_info[0].opaque_back = opaque_back;
      } else if (bkrd_is_gradient) {
        thread_info[0].bkrd_top = bkrd_top;
        thread_info[0].bkrd_bottom = bkrd_bottom;
        thread_info[0].opaque_back = opaque_back;
      } else {
        thread_info[0].background = background;
      }
      thread_info[0].bytes = width * (unsigned int)height;
      thread_info[0].ray = I; /* for compute box */
      thread_info[0].size_hint = I->PrimSize;
      /* shadow map */

      {
        int bc;
        float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_shadow);
        for (bc = 2; bc < I->NBasis; bc++) {
          int* vert2prim_ptr =
              I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();
          thread_info[bc - 1].basis = I->Basis + bc;
          thread_info[bc - 1].vert2prim = vert2prim_ptr;
          thread_info[bc - 1].prim = I->Primitive;
          thread_info[bc - 1].n_prim = I->NPrimitive;
          thread_info[bc - 1].clipBox = nullptr;
          thread_info[bc - 1].phase = bc - 1;
          thread_info[bc - 1].perspective = false;
          thread_info[bc - 1].front = _0;
          /* allowing these maps to be more fine helps performance */
          thread_info[bc - 1].size_hint = I->PrimSize * factor;
        }
      }

      /* NOTE that we're not limiting the number of threads in this phase
         under the assumption that it will usually just be a few threads */
      RayHashSpawn(thread_info, n_thread, I->NBasis - 1);

      FreeP(thread_info);
    } else if (ok) {
      int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();
      /* the GPU tracer builds its own (column) lists from the voxel grid,
         so skip PyMOL's express lists when it is going to run */
      g_basis_lite_request = RayMetalEnabled() && !bkrd_is_gradient &&
                             !I->bkgrd_data &&
                             !RayMetalUnsupported(I, perspective);
      ok &=
          BasisMakeMap(I->Basis + 1, vert2prim_ptr, I->Primitive, I->NPrimitive,
                       I->Volume, perspective, front, I->PrimSize);
      g_basis_lite_request = 0;
      if (ok && shadows) {
        int bc;
        float factor = SettingGetGlobal_f(I->G, cSetting_ray_hint_shadow);
        for (bc = 2; ok && bc < I->NBasis; bc++) {
          ok &= BasisMakeMap(I->Basis + bc, vert2prim_ptr, I->Primitive,
                             I->NPrimitive, nullptr, false, _0,
                             I->PrimSize * factor);
        }
      }

      /* serial tasks which RayHashThread does in parallel mode using the first
       * thread */

      if (ok) {
        if (I->bkgrd_data) {
          fill_background_image(I, image, width, height,
                                width * (unsigned int)height);
        } else if (bkrd_is_gradient) {
          fill_gradient(I, opaque_back, image, bkrd_top, bkrd_bottom, width,
                        height, width * (unsigned int)height);
        } else {
          fill(image, background, width * (unsigned int)height);
        }
        RayComputeBox(I);
      }
    }
    OrthoBusyFast(I->G, 5, 20);
    RM_STAMP(); /* 2: map + fill */
    RayProgress(0.18f);
    now = UtilGetSeconds(I->G) - timing;

    if (ok) {
      if (shadows) {
        PRINTFB(I->G, FB_Ray, FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], [%4.2f:%dx%dx%d], %4.2f sec.\n",
            I->Basis[1].Map->Div, I->Basis[1].Map->Dim[0],
            I->Basis[1].Map->Dim[1], I->Basis[1].Map->Dim[2],
            I->Basis[2].Map->Div, I->Basis[2].Map->Dim[0],
            I->Basis[2].Map->Dim[2], I->Basis[2].Map->Dim[2], now ENDFB(I->G);
      } else {
        PRINTFB(I->G, FB_Ray, FB_Blather)
        " Ray: voxels: [%4.2f:%dx%dx%d], %4.2f sec.\n", I->Basis[1].Map->Div,
            I->Basis[1].Map->Dim[0], I->Basis[1].Map->Dim[1],
            I->Basis[1].Map->Dim[2], now ENDFB(I->G);
      }
    }
    /* IMAGING */

    if (ok) {
      /* now spawn threads as needed */
      CRayThreadInfo* rt = pymol::calloc<CRayThreadInfo>(n_thread);

      int x_start = 0, y_start = 0;
      int x_stop = 0, y_stop = 0;
      float x_test = _0, y_test = _0;
      int x_pixel, y_pixel;

      if (perspective) {
        int c;

        if (I->min_box[2] > -front) I->min_box[2] = -front;
        if (I->max_box[2] > -front) I->max_box[2] = -front;

        for (c = 0; c < 4; c++) {
          switch (c) {
            case 0:
              x_test = -I->min_box[0] / I->min_box[2];
              y_test = -I->min_box[1] / I->min_box[2];
              break;
            case 1:
              x_test = -I->min_box[0] / I->max_box[2];
              y_test = -I->min_box[1] / I->max_box[2];
              break;
            case 2:
              x_test = -I->max_box[0] / I->min_box[2];
              y_test = -I->max_box[1] / I->min_box[2];
              break;
            case 3:
              x_test = -I->max_box[0] / I->max_box[2];
              y_test = -I->max_box[1] / I->max_box[2];
              break;
          }

          /* project onto back to get the effective range */

          x_pixel = (int)(width * (((x_test * I->Volume[5]) - I->Volume[0]) /
                                   I->Range[0]));
          y_pixel = (int)(height * (((y_test * I->Volume[5]) - I->Volume[2]) /
                                    I->Range[1]));

          if (!c) {
            x_start = x_pixel;
            x_stop = x_pixel;
            y_start = y_pixel;
            y_stop = y_pixel;
          } else {
            if (x_start > x_pixel) x_start = x_pixel;
            if (x_stop < x_pixel) x_stop = x_pixel;
            if (y_start > y_pixel) y_start = y_pixel;
            if (y_stop < y_pixel) y_stop = y_pixel;
          }
        }
        x_start -= 2;
        x_stop += 2;
        y_start -= 2;
        y_stop += 2;

        /*
           x_start = 0;
           y_start = 0;
           x_stop = width;
           y_stop = height; */

      } else {
        x_start =
            (int)((width * (I->min_box[0] - I->Volume[0])) / I->Range[0]) - 2;
        x_stop =
            (int)((width * (I->max_box[0] - I->Volume[0])) / I->Range[0]) + 2;

        y_stop =
            (int)((height * (I->max_box[1] - I->Volume[2])) / I->Range[1]) + 2;
        y_start =
            (int)((height * (I->min_box[1] - I->Volume[2])) / I->Range[1]) - 2;
      }
      if (x_start < 0) x_start = 0;
      if (y_start < 0) y_start = 0;
      if (x_stop > width) x_stop = width;
      if (y_stop > height) y_stop = height;
      trace_x0 = x_start;
      trace_x1 = x_stop;
      trace_y0 = y_start;
      trace_y1 = y_stop;

      for (a = 0; a < n_thread; a++) {
        rt[a].ray = I;
        rt[a].width = width;
        rt[a].height = height;
        rt[a].x_start = x_start;
        rt[a].x_stop = x_stop;
        rt[a].y_start = y_start;
        rt[a].y_stop = y_stop;
        rt[a].image = image;
        rt[a].border = mag - 1;
        rt[a].front = front;
        rt[a].back = back;
        rt[a].fore_mask = fore_mask;
        rt[a].bkrd_is_gradient = bkrd_is_gradient;
        if (bkrd_is_gradient) {
          rt[a].bkrd_top = bkrd_top;
          rt[a].bkrd_bottom = bkrd_bottom;
        } else {
          rt[a].bkrd_top = bkrd_top;
          rt[a].bkrd_bottom = bkrd_top;
        }
        rt[a].ambient = ambient;
        rt[a].background = background;
        rt[a].phase = a;
        rt[a].n_thread = n_thread;
        rt[a].edging = nullptr;
        rt[a].edging_cutoff =
            oversample_cutoff; /* info needed for busy indicator */
        rt[a].perspective = perspective;
        rt[a].fov = fov;
        rt[a].pos[2] = pos[2];
        rt[a].depth = depth;
        if (I->bkgrd_data) {
          rt[a].bgWidth = I->bkgrd_data->getWidth();
          rt[a].bgHeight = I->bkgrd_data->getHeight();
        }
        rt[a].bkrd_data = I->bkgrd_data ? I->bkgrd_data->bits() : nullptr;
      }

      bool traced_on_gpu =
          RayMetalEnabled() &&
          RayTraceMetal(rt, trace_word, ray_trace_mode, &outlined_on_gpu,
                        antialias > 1 ? image_copy : nullptr, &aa_on_gpu);
      outlined_on_gpu &= traced_on_gpu;
      aa_on_gpu &= traced_on_gpu;
      if (!RayMetalEnabled()) RayMetalSetLastReport("cpu: GPU tracer off");
      /* the CPU tracer (fallback, or the edge-oversampling pass) needs the
         express lists a lite map skipped */
      if ((!traced_on_gpu || oversample_cutoff) && I->Basis[1].Lite)
        ok &= BasisLiteFinish(I->Basis + 1);
      if (!traced_on_gpu)
        if (n_thread > 1)
          RayTraceSpawn(rt, n_thread);
        else
          RayTraceThread(rt);
      RM_STAMP(); /* 3: trace */

      if (oversample_cutoff) { /* perform edge oversampling, if requested */
        unsigned int* edging;

        edging = pymol::malloc<unsigned int>(buffer_size);

        memcpy(edging, image, buffer_size * sizeof(unsigned int));

        for (a = 0; a < n_thread; a++) {
          rt[a].edging = edging;
        }

        if (n_thread > 1)
          RayTraceSpawn(rt, n_thread);
        else
          RayTraceThread(rt);

        FreeP(edging);
      }
      FreeP(rt);
    }
  }

  RM_STAMP(); /* 4: edge oversampling */
  if (ok && depth && ray_trace_mode && !outlined_on_gpu) {
    /* zeroed: outside the processed area it must read as "no gradient" */
    float* delta = pymol::calloc<float>(3 * width * height);
    int x, y;
    /* rows/columns that can change: the traced rectangle plus a margin for
       the gradient (1), the three fog-fill passes (3) and the outline
       neighbourhood (1); everywhere else depth and gradient are zero and
       the pixel is left as is. The whole image if nothing was traced. */
    const int m = 6;
    /* modes 2 and 3 rewrite every non-edge pixel: whole image */
    const bool rect =
        ray_trace_mode == 1 && trace_x1 > trace_x0 && trace_y1 > trace_y0;
    const int ox0 = rect ? std::max(0, trace_x0 - m) : 0;
    const int ox1 = rect ? std::min((int)width, trace_x1 + m) : (int)width;
    const int oy0 = rect ? std::max(0, trace_y0 - m) : 0;
    const int oy1 = rect ? std::min((int)height, trace_y1 + m) : (int)height;
    ErrChkPtr(I->G, delta);
    if (ok) {
      RayParallelRows(oy1 - oy0, [&](int yr) {
        const int y = oy0 + yr;
        int x, xc, yc;
        float d, dzdx, dzdy, *p, *q, dd;
        p = depth + (size_t)y * width + ox0;
        q = delta + (size_t)3 * ((size_t)y * width + ox0);
        for (x = ox0; x < ox1; x++) {
          dzdx = 0.0F;
          dzdy = 0.0F;
          xc = 0;
          yc = 0;
          d = *p;
          if (x) {
            dd = d - p[-1];
            dzdx += dd;
            xc++;
          }
          if (x < (width - 1)) {
            dd = p[1] - d;
            if ((!xc) || (fabs(dd) > fabs(dzdx))) dzdx = dd;
            xc = 1;
          }
          if (y) {
            dd = d - p[-width];
            dzdy += dd;
            yc++;
          }
          if (y < (height - 1)) {
            dd = p[width] - d;
            if ((!yc) || (fabs(dd) > fabs(dzdy))) dzdy = dd;
            yc = 1;
          }
          p++;
          *(q++) = dzdx;
          *(q++) = dzdy;
          /*
             if(((x == y )||(x==width/2)||(y==height/2))&&((dzdx!=0.0F) ||
             (dzdy!=0.0F))) { printf("%5d %5d : %8.3f %8.3f\n",y,x,dzdx,dzdy);
             } */

          *(q++) = sqrt1f(dzdx * dzdx + dzdy * dzdy);
        }
      });
    }

    if (ok) {
      int i;
      {
        const float _1 = 1.0F;

        float invFrontMinusBack = _1 / (front - back);
        float inv1minusFogStart = _1;
        int fogFlag = false;
        float fog_start = 0.0F;
        int fogRangeFlag = false;
        float fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
        if (fog < 0.0F) {
          if (SettingGetGlobal_b(I->G, cSetting_depth_cue)) {
            fog = SettingGetGlobal_f(I->G, cSetting_fog);
          } else
            fog = _0;
        }

        if (fog != _0) {
          if (fog > 1.0F) fog = 1.0F;
          fogFlag = true;
          fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
          if (fog_start < 0.0F)
            fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
          if (fog_start > 1.0F) fog_start = 1.0F;
          if (fabs(fog_start) > R_SMALL4) {
            fogRangeFlag = true;
            if (fabs(fog_start - 1.0F) < R_SMALL4) /* prevent div/0 */
              fogFlag = false;
          }
          inv1minusFogStart = _1 / (_1 - fog_start);
        }

        if (fogFlag) { /* make sure we have depth values at every potentially
                          drawn pixel */
          /* zeroed: pixels with no depth and no neighbour depth are not
             written below, and an uninitialised buffer made the fogged
             outline at those pixels depend on leftover heap contents */
          float* tmp = pymol::calloc<float>(width * height);
          float dep;
          float *p, *q;
          int cnt;
          for (i = 0; i < 3; i++) { /* three passes required */
            {
              float *const src = depth, *const dst = tmp;
              RayParallelRows(oy1 - oy0, [&, src, dst](int yr) {
                const int y = oy0 + yr;
                float dep;
                int x, cnt;
                float *p = src + (size_t)y * width + ox0,
                      *q = dst + (size_t)y * width + ox0;
                for (x = ox0; x < ox1; x++) {
                  if (fabs(*p) < R_SMALL4) {
                    dep = 0.0F;
                    cnt = 0;
                    if (x) {
                      if (fabs(p[-1]) > R_SMALL4) {
                        dep += p[-1];
                        cnt++;
                      }
                    }
                    if (x < (width - 1)) {
                      if (fabs(p[1]) > R_SMALL4) {
                        dep += p[1];
                        cnt++;
                      }
                    }
                    if (y) {
                      if (fabs(p[-width]) > R_SMALL4) {
                        dep += p[-width];
                        cnt++;
                      }
                    }
                    if (y < (height - 1)) {
                      if (fabs(p[width]) > R_SMALL4) {
                        dep += p[width];
                        cnt++;
                      }
                    }
                    if (cnt) {
                      dep /= cnt;
                      *q = dep;
                    }
                  } else {
                    *q = *p;
                  }
                  p++;
                  q++;
                }
              });
            }
            p = tmp;
            tmp = depth;
            depth = p;
          }
          FreeP(tmp);
        }
        {
          unsigned int* q = image;
          float* p = delta;
          int width3 = width * 3;

          float slope_f =
              SettingGetGlobal_f(I->G, cSetting_ray_trace_slope_factor);
          float depth_f =
              SettingGetGlobal_f(I->G, cSetting_ray_trace_depth_factor);
          float disco_f =
              SettingGetGlobal_f(I->G, cSetting_ray_trace_disco_factor);
          float diff, max_depth;
          float dot, min_dot, max_slope, max_dz, max_pz;
          float dx, dy, dz, px = 0.0F, py = 0.0F, pz = 0.0F, ddx, ddy;
          const float _8 = 0.08F;
          const float _4 = 0.4F;
          const float _25 = 0.25F;
          const float _m25 = -0.25F;
          float disco_f_625 = disco_f * 0.625F;
          float disco_f_5 = disco_f * 0.5F;
          float disco_f_45 = disco_f * 0.45F;

          {
            float gain = I->PixelRadius /
                         (SettingGetGlobal_f(I->G, cSetting_ray_trace_gain) *
                          I->Magnified);
            if (antialias) gain /= antialias;
            slope_f *= gain;
            depth_f *= gain;
            disco_f *= gain;
          }
          RayParallelRows(oy1 - oy0, [&](int yr) {
            /* rows run in parallel; each starts from what the serial loop
               carried into it (px/py/pz: case 5 of the row above's last
               pixel — only read at x == 0; elsewhere case 0 overwrites it),
               so the result is unchanged */
            const int y = oy0 + yr;
            unsigned int* q = image + (size_t)y * width + ox0;
            float* p = delta + (size_t)3 * ((size_t)y * width + ox0);
            unsigned int background_row = background;
            float diff, max_depth, dot, min_dot, max_slope, max_dz, max_pz;
            float dx, dy, dz, px = 0.0F, py = 0.0F, pz = 0.0F, ddx, ddy;
            int x, i;
            if (y && width > 1) {
              const float* c = delta + ((size_t)y * width + width - 2) * 3;
              px = c[0];
              py = c[1];
              pz = c[2];
            }
            float bkrd[3], perc = 0.;
            if (bkrd_is_gradient) {
              /* for RayRender, y is from bottom to top */
              perc = y / (float)height;
              bkrd[0] = bkrd_bottom[0] + perc * (bkrd_top[0] - bkrd_bottom[0]);
              bkrd[1] = bkrd_bottom[1] + perc * (bkrd_top[1] - bkrd_bottom[1]);
              bkrd[2] = bkrd_bottom[2] + perc * (bkrd_top[2] - bkrd_bottom[2]);
              if (I->BigEndian) {
                background_row =
                    back_mask |
                    ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))) << 24) |
                    ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 16) |
                    ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 8);
              } else {
                background_row =
                    back_mask |
                    ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 16) |
                    ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 8) |
                    ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))));
              }
            }
            for (x = ox0; x < ox1; x++) {
              max_slope = 0.0F;
              max_depth = 0.0F;
              min_dot = 1.0F;
              max_dz = 0.0F;
              max_pz = 0.0F;
              dx = p[0];
              dy = p[1];
              dz = p[2];
              for (i = 0; i < 8; i++) {
                switch (i) {
                  case 0:
                    if (x) {
                      px = p[-3];
                      py = p[-2];
                      pz = p[-1];
                    }
                    break;
                  case 1:
                    if (x < (width - 1)) {
                      px = p[3];
                      py = p[4];
                      pz = p[5];
                    }
                    break;
                  case 2:
                    if (y) {
                      px = p[-width3];
                      py = p[-width3 + 1];
                      pz = p[-width3 + 2];
                    }
                    break;
                  case 3:
                    if (y < (height - 1)) {
                      px = p[width3];
                      py = p[width3 + 1];
                      pz = p[width3 + 2];
                    }
                    break;
                  case 4:
                    if (x && y) {
                      px = p[-width3 - 3];
                      py = p[-width3 - 2];
                      pz = p[-width3 - 1];
                    }
                    break;
                  case 5:
                    if (x && (y < (height - 1))) {
                      px = p[width3 - 3];
                      py = p[width3 - 2];
                      pz = p[width3 - 1];
                    }
                    break;
                  case 6:
                    if (y && (x < (width - 1))) {
                      px = p[-width3 + 3];
                      py = p[-width3 + 4];
                      pz = p[-width3 + 5];
                    }
                    break;
                  case 7:
                    if ((y < (height - 1)) && (x < (width - 1))) {
                      px = p[width3 + 3];
                      py = p[width3 + 4];
                      pz = p[width3 + 5];
                    }
                    break;
                }
                ddx = dx - px;
                ddy = dy - py;
                diff = ddx * ddx + ddy * ddy;
                if (max_depth < diff) max_depth = diff;
                if ((dz > R_SMALL4) && (pz > R_SMALL4)) {
                  dot = (dx / dz) * (px / pz) + (dy / dz) * (py / pz);
                  if (dot < min_dot) {
                    min_dot = dot;
                    max_dz = dz;
                    max_pz = pz;
                  }
                }
                /*                if(dz>max_dz) max_dz = dz;
                   if(pz>max_pz) max_pz = pz; */
                diff = fabs(dz - pz);
                if (diff > max_slope) max_slope = diff;
              }
              if ((max_slope > (slope_f)) /* depth */
                  || (max_depth > (depth_f))
                  /* slope */
                  /* gradient discontinuities -- could probably use more
                     tuning... */
                  || ((min_dot < _8) &&
                      ((max_dz > disco_f) || (max_pz > disco_f))) ||
                  ((min_dot < _4) &&
                   ((max_dz > disco_f_625) || (max_pz > disco_f_625))) ||
                  ((min_dot < _25) &&
                   ((max_dz > disco_f_5) || (max_pz > disco_f_5))) ||
                  ((min_dot < _m25) &&
                   ((max_dz > disco_f_45) && (max_pz > disco_f_45)))) {
                if (fogFlag) {
                  float ffact = depth[q - image] * invFrontMinusBack;
                  float ffact1m;
                  float fc[4];
                  unsigned int cc0, cc1, cc2, cc3;

                  if (fogRangeFlag)
                    ffact = (ffact - fog_start) * inv1minusFogStart;

                  ffact *= fog;

                  if (ffact < _0) ffact = _0;
                  if (ffact > _1) ffact = _1;

                  ffact1m = _1 - ffact;

                  if (orig_opaque_back) {
                    fc[0] = (0xFF & (background_row >> 24)) * ffact +
                            (0xFF & (trace_word >> 24)) * ffact1m;
                    fc[1] = (0xFF & (background_row >> 16)) * ffact +
                            (0xFF & (trace_word >> 16)) * ffact1m;
                    fc[2] = (0xFF & (background_row >> 8)) * ffact +
                            (0xFF & (trace_word >> 8)) * ffact1m;
                    fc[3] = (0xFF & (background_row)) * ffact +
                            (0xFF & (trace_word)) * ffact1m;
                  } else { /* if non-opaque background_row, then use alpha to
                              blend */
                    fc[1] = (0xFF & (trace_word >> 16));
                    fc[2] = (0xFF & (trace_word >> 8));
                    if (I->BigEndian) {
                      fc[0] = (0xFF & (trace_word >> 24));
                      fc[3] = (0xFF & (background_row)) * ffact +
                              (0xFF & (trace_word)) * ffact1m;
                    } else {
                      fc[0] = (0xFF & (background_row >> 24)) * ffact +
                              (0xFF & (trace_word >> 24)) * ffact1m;
                      fc[3] = (0xFF & (trace_word));
                    }
                  }
                  cc0 = (uint)(fc[0]);
                  cc1 = (uint)(fc[1]);
                  cc2 = (uint)(fc[2]);
                  cc3 = (uint)(fc[3]);

                  if (cc0 > 255) cc0 = 255;
                  if (cc1 > 255) cc1 = 255;
                  if (cc2 > 255) cc2 = 255;
                  if (cc3 > 255) cc3 = 255;

                  *q = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
                } else {
                  *q = trace_word;
                }
              } else if (ray_trace_mode == 2) { /* only draw edge */
                *q = background_row;
              } else if (ray_trace_mode == 3) { /* quantize */
                *q = (*q & 0xC0C0C0C0);
                *q = *q | ((*q) >> 2) | ((*q) >> 4) | ((*q) >> 6);
              }
              p += 3;
              q++;
            }
          });
        }
      }
    }
    FreeP(delta);
  }

  RM_STAMP(); /* 5: outline */
  RayProgress(0.97f);
  if (ok && aa_on_gpu) {
    FreeP(image);
    image = image_copy;
  } else if (ok && antialias > 1) {
    /* now spawn threads as needed */
    /* each scanline is independent, so any thread count gives the same
       image: use every core (GCD) rather than Python-spawned threads */
    int n_thread = std::max(1, (int)std::thread::hardware_concurrency());
    CRayAntiThreadInfo* rt = pymol::calloc<CRayAntiThreadInfo>(n_thread);

    for (a = 0; a < n_thread; a++) {
      rt[a].width = width;
      rt[a].height = height;
      rt[a].image = image;
      rt[a].image_copy = image_copy;
      rt[a].phase = a;
      rt[a].mag = mag; /* fold magnification */
      rt[a].n_thread = n_thread;
      rt[a].ray = I;
      rt[a].no_busy = 1;
    }

    dispatch_apply_f(
        n_thread, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), rt,
        [](void* c, size_t k) { RayAntiThread(((CRayAntiThreadInfo*)c) + k); });
    FreeP(rt);
    FreeP(image);
    image = image_copy;
  }

  PRINTFD(I->G, FB_Ray)
  " RayRender: n_hit %d\n", n_hit ENDFD;

  if (ok) {
    /* EXPERIMENTAL RAY-VOLUME CODE */
    volume = SettingGetGlobal_b(I->G, cSetting_ray_volume);

    if (volume) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          float dd = depth[x + width * y];
          if (dd == 0.0) dd = -back;
          depth[x + width * y] = -dd / (back - front) + 0.1;
        }
      }
      if (rayDepthPixels) FreeP(rayDepthPixels);
      rayDepthPixels = depth;
      rayWidth = width;
      rayHeight = height;
      rayVolume = 3;
    } else
      FreeP(depth);
  }
  I->bkgrd_data = nullptr;
  RM_STAMP(); /* 6: antialias + rest */
  if (rm_n == 7) {
    char msg[256];
    snprintf(msg, sizeof(msg),
             "expand %.0f, map %.0f, trace %.0f, edge %.0f, outline %.0f, aa "
             "%.0f ms (%dx%d)",
             1e3 * (rm_t[1] - rm_t[0]), 1e3 * (rm_t[2] - rm_t[1]),
             1e3 * (rm_t[3] - rm_t[2]), 1e3 * (rm_t[4] - rm_t[3]),
             1e3 * (rm_t[5] - rm_t[4]), 1e3 * (rm_t[6] - rm_t[5]), width,
             height);
    RayMetalSetStages(msg);
  }
}

void RayRenderColorTable(CRay* I, int width, int height, int* image) {
  int x, y;
  unsigned int r = 0, g = 0, b = 0;
  unsigned int *pixel, mask, *p;

  if (I->BigEndian)
    mask = 0x000000FF;
  else
    mask = 0xFF000000;

  p = (unsigned int*)image;
  for (x = 0; x < width; x++)
    for (y = 0; y < height; y++) *(p++) = mask;

  if ((width >= 512) && (height >= 512)) {
    for (y = 0; y < 512; y++)
      for (x = 0; x < 512; x++) {
        pixel = (unsigned int*)(image + ((width)*y) + x);
        if (I->BigEndian) {
          *(pixel) = mask | (r << 24) | (g << 16) | (b << 8);
        } else {
          *(pixel) = mask | (b << 16) | (g << 8) | r;
        }
        b = b + 4;
        if (!(0xFF & b)) {
          b = 0;
          g = g + 4;
          if (!(0xFF & g)) {
            g = 0;
            r = r + 4;
          }
        }
      }
  }
}
CRay* RayNew(PyMOLGlobals* G, int antialias) {
  unsigned int test;
  unsigned char* testPtr;
  int a;

  auto I = new CRay();
  I->G = G;
  test = 0xFF000000;
  testPtr = (unsigned char*)&test;
  I->BigEndian = (*testPtr) & 0x01;
  I->Trans = 0.0F;
  I->Wobble = 0;
  I->TTTFlag = false;
  zero3f(I->WobbleParam);
  PRINTFB(I->G, FB_Ray, FB_Blather)
  " RayNew: BigEndian = %d\n", I->BigEndian ENDFB(I->G);

  I->Basis = pymol::malloc<CBasis>(12);
  BasisInit(I->G, I->Basis);
  BasisInit(I->G, I->Basis + 1);
  I->NBasis = 2;
  I->Primitive = nullptr;
  I->NPrimitive = 0;
  I->CheckInterior = false;
  if (antialias < 0) antialias = SettingGetGlobal_i(I->G, cSetting_antialias);
  I->Sampling = antialias;
  if (I->Sampling < 2) /* always supersample fonts by at least 2X */
    I->Sampling = 2;
  for (a = 0; a < 256; a++) {
    I->Random[a] = (float)((rand() / (1.0 + RAND_MAX)) - 0.5);
  }

  I->Wobble = SettingGet_i(I->G, nullptr, nullptr, cSetting_ray_texture);
  {
    auto v = SettingGet<const float*>(I->G, cSetting_ray_texture_settings);
    int color = SettingGetGlobal_color(I->G, cSetting_ray_interior_color);
    copy3f(v, I->WobbleParam);
    v = ColorGet(I->G, color);
    copy3f(v, I->IntColor);
  }

  return (I);
}
void RayPrepare(CRay* I, float v0, float v1, float v2, float v3, float v4,
                float v5, float fov, glm::vec3 pos, float* mat,
                const glm::mat4& rotMat, float aspRat, int width, int height,
                float pixel_scale, int ortho, float pixel_ratio,
                float front_back_ratio, float magnified)

/*prepare for vertex calls */
{
  int a;
  if (!I->Primitive) I->Primitive = VLAlloc(CPrimitive, 10000);
  I->Volume[0] = v0;
  I->Volume[1] = v1;
  I->Volume[2] = v2;
  I->Volume[3] = v3;
  I->Volume[4] = v4;
  I->Volume[5] = v5;
  I->Range[0] = I->Volume[1] - I->Volume[0];
  I->Range[1] = I->Volume[3] - I->Volume[2];
  I->Range[2] = I->Volume[5] - I->Volume[4];
  I->AspRatio = aspRat;
  I->Width = width;
  I->Height = height;
  CharacterSetRetention(I->G, true);

  if (mat)
    for (a = 0; a < 16; a++) I->ModelView[a] = mat[a];
  else {
    identity44f(I->ModelView);
  }
  identity44f(I->ProMatrix);
  if (ortho) {
    I->ProMatrix[0] = 2.f / I->Range[0];
    I->ProMatrix[5] = 2.f / I->Range[1];
    I->ProMatrix[10] = -2.f / I->Range[2];
    I->ProMatrix[12] = -(I->Volume[0] + I->Volume[1]) / I->Range[0];
    I->ProMatrix[13] = -(I->Volume[2] + I->Volume[3]) / I->Range[1];
    I->ProMatrix[14] = -(I->Volume[4] + I->Volume[5]) / I->Range[2];
  } else {
    I->ProMatrix[0] = I->Volume[4] / (front_back_ratio * I->Volume[1]);
    I->ProMatrix[5] = I->Volume[4] / (front_back_ratio * I->Volume[3]);
    I->ProMatrix[10] = -(I->Volume[5] + I->Volume[4]) / I->Range[2];
    I->ProMatrix[11] = -1.f;
    I->ProMatrix[14] = (-2.f * I->Volume[5] * I->Volume[4]) / I->Range[2];
    I->ProMatrix[15] = 0.f;
  }
  I->Rotation = rotMat;
  I->Ortho = ortho;
  if (ortho) {
    I->PixelRadius = (((float)I->Range[0]) / width) * pixel_scale;
  } else {
    I->PixelRadius = (((float)I->Range[0]) / width) * pixel_scale * pixel_ratio;
  }
  I->PixelRatio = pixel_ratio;
  I->Magnified = magnified;
  I->FrontBackRatio = front_back_ratio;
  I->PrimSizeCnt = 0;
  I->PrimSize = 0.0;
  I->Fov = fov;
  I->Pos = pos;
}
void RaySetTTT(CRay* I, int flag, float* ttt) {
  I->TTTFlag = flag;
  if (I->TTTFlag) {
    I->TTT = glm::make_mat4(ttt);
  }
}

void RayGetTTT(CRay* I, float* ttt) {
  if (!I->TTTFlag) {
    identity44f(ttt);
  } else {
    std::copy_n(glm::value_ptr(I->TTT), 16, ttt);
  }
}
void RayRelease(CRay* I) {
  int a;

  for (a = 0; a < I->NBasis; a++) {
    BasisFinish(&I->Basis[a]);
  }
  I->NBasis = 0;
  VLAFreeP(I->Primitive);
}
void RayFree(CRay* I) {
  RayRelease(I);
  CharacterSetRetention(I->G, false);
  FreeP(I->Basis);
  DeleteP(I);
}
void RayPushTTT(CRay* I) {
  if (I->TTTFlag) {
    I->TTTStack.push_back(I->TTT);
  }
}

void RayPopTTT(CRay* I) {
  if (I->TTTStack.empty()) {
    I->TTTFlag = false;
  } else {
    I->TTT = I->TTTStack.back();
    I->TTTStack.pop_back();
    I->TTTFlag = true;
  }
}

static void RayApplyMatrix33(unsigned int n, float3* q, const float m[16],
                             float3* p) {
  {
    unsigned int i;
    float m0 = m[0], m4 = m[4], m8 = m[8], m12 = m[12];
    float m1 = m[1], m5 = m[5], m9 = m[9], m13 = m[13];
    float m2 = m[2], m6 = m[6], m10 = m[10], m14 = m[14];
    for (i = 0; i < n; i++) {
      float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
      q[i][0] = m0 * p0 + m4 * p1 + m8 * p2 + m12;
      q[i][1] = m1 * p0 + m5 * p1 + m9 * p2 + m13;
      q[i][2] = m2 * p0 + m6 * p1 + m10 * p2 + m14;
    }
  }
}

void RayApplyMatrixInverse33(unsigned int n, float3* q, const float m[16],
                             float3* p) {
  {
    unsigned int i;
    float m0 = m[0], m4 = m[4], m8 = m[8], m12 = m[12];
    float m1 = m[1], m5 = m[5], m9 = m[9], m13 = m[13];
    float m2 = m[2], m6 = m[6], m10 = m[10], m14 = m[14];
    for (i = 0; i < n; i++) {
      float p0 = p[i][0] - m12, p1 = p[i][1] - m13, p2 = p[i][2] - m14;
      q[i][0] = m0 * p0 + m1 * p1 + m2 * p2;
      q[i][1] = m4 * p0 + m5 * p1 + m6 * p2;
      q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
    }
  }
}

void RayTransformNormals33(unsigned int n, float3* q, const float m[16],
                           float3* p) {
  unsigned int i;
  float m0 = m[0], m4 = m[4], m8 = m[8];
  float m1 = m[1], m5 = m[5], m9 = m[9];
  float m2 = m[2], m6 = m[6], m10 = m[10];
  for (i = 0; i < n; i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m4 * p1 + m8 * p2;
    q[i][1] = m1 * p0 + m5 * p1 + m9 * p2;
    q[i][2] = m2 * p0 + m6 * p1 + m10 * p2;
  }
  for (i = 0; i < n;
       i++) { /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}

void RayTransformInverseNormals33(unsigned int n, float3* q, const float m[16],
                                  float3* p) {
  unsigned int i;
  float m0 = m[0], m4 = m[4], m8 = m[8];
  float m1 = m[1], m5 = m[5], m9 = m[9];
  float m2 = m[2], m6 = m[6], m10 = m[10];
  for (i = 0; i < n; i++) {
    float p0 = p[i][0], p1 = p[i][1], p2 = p[i][2];
    q[i][0] = m0 * p0 + m1 * p1 + m2 * p2;
    q[i][1] = m4 * p0 + m5 * p1 + m6 * p2;
    q[i][2] = m8 * p0 + m9 * p1 + m10 * p2;
  }
  for (i = 0; i < n;
       i++) { /* renormalize - can we do this to the matrix instead? */
    normalize3f(q[i]);
  }
}
void RayGetScreenVertex(CRay* I, float* v, float* res) {
  MatrixTransformC44f4f(I->ModelView, v, res);
  normalize4f(res);
}

void RayAdjustZtoScreenZ(CRay* ray, float* pos, float zarg) {
  PyMOLGlobals* G = ray->G;
  float BackSafe = ray->Volume[5], FrontSafe = ray->Volume[4];
  float clipRange = (BackSafe - FrontSafe);
  float z = (zarg + 1.f) / 2.f;
  float zInPreProj = -(z * clipRange + FrontSafe);
  float pos4[4], tpos[4], npos[4];
  float InvModMatrix[16];
  copy3f(pos, pos4);
  pos4[3] = 1.f;
  MatrixTransformC44f4f(ray->ModelView, pos4, tpos);
  normalize4f(tpos);
  /* NEED TO ACCOUNT FOR ORTHO */
  if (SettingGetGlobal_b(G, cSetting_ortho)) {
    npos[0] = tpos[0];
    npos[1] = tpos[1];
  } else {
    npos[0] = zInPreProj * tpos[0] / tpos[2];
    npos[1] = zInPreProj * tpos[1] / tpos[2];
  }
  npos[2] = zInPreProj;
  npos[3] = 1.f;
  MatrixInvertC44f(ray->ModelView, InvModMatrix);
  MatrixTransformC44f4f(InvModMatrix, npos, npos);
  normalize4f(npos);
  copy3f(npos, pos);
}

static float RayGetRawDepth(CRay* ray, float* pos) {
  float pos4[4], tpos[4];
  copy3f(pos, pos4);
  pos4[3] = 1.f;
  MatrixTransformC44f4f(ray->ModelView, pos4, tpos);
  normalize4f(tpos);
  return -tpos[2];
}

void RayAdjustZtoScreenZofPoint(CRay* ray, float* pos, float* zpoint) {
  float BackSafe = ray->Volume[5], FrontSafe = ray->Volume[4];
  float clipRange = (BackSafe - FrontSafe);
  float zInClipping = RayGetRawDepth(ray, zpoint);
  float z = ((2.f * (zInClipping - FrontSafe) / clipRange) - 1.f);
  RayAdjustZtoScreenZ(ray, pos, z);
}

void RaySetPointToWorldScreenRelative(CRay* ray, float* pos, float* screenPt) {
  float npos[4];
  float PmvMatrix[16], InvPmvMatrix[16];
  int width = ray->Width, height = ray->Height;
  multiply44f44f44f(ray->ModelView, RayGetProMatrix(ray), PmvMatrix);

  npos[0] = (floor(screenPt[0] * width)) / width;
  npos[1] = (floor(screenPt[1] * height)) / height;
  npos[2] = 0.f;
  npos[3] = 1.f;
  MatrixInvertC44f(PmvMatrix, InvPmvMatrix);
  MatrixTransformC44f4f(InvPmvMatrix, npos, npos);
  normalize4f(npos);
  RayAdjustZtoScreenZ(ray, npos, screenPt[2]);
  copy3f(npos, pos);
}

float RayGetScaledAllAxesAtPoint(CRay* I, float* pt, float* xn, float* yn,
                                 float* zn) {
  float xn0[3] = {1.0F, 0.0F, 0.0F};
  float yn0[3] = {0.0F, 1.0F, 0.0F};
  float zn0[3] = {0.0F, 0.0F, 1.0F};
  float v_scale;

  v_scale = RayGetScreenVertexScale(I, pt) / I->Sampling;

  RayApplyMatrixInverse33(1, (float3*)xn0, glm::value_ptr(I->Rotation),
                          (float3*)xn0);
  RayApplyMatrixInverse33(1, (float3*)yn0, glm::value_ptr(I->Rotation),
                          (float3*)yn0);
  RayApplyMatrixInverse33(1, (float3*)zn0, glm::value_ptr(I->Rotation),
                          (float3*)zn0);

  scale3f(xn0, v_scale, xn);
  scale3f(yn0, v_scale, yn);
  scale3f(zn0, v_scale, zn);
  return v_scale;
}
float* RayGetProMatrix(CRay* I) { return I->ProMatrix; }
