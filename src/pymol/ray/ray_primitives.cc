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

// Adding spheres, cylinders, cones, triangles, labels and other primitives to a
// ray scene.

#include "pymol/base/err.h"
#include "pymol/ray/ray.h"
#include "pymol/glyph/character.h"
#include "pymol/text/text.h"
#include "pymol/graphics/cgo.h"
#include "pymol/ray/ray_internal.h"

static void RayApplyContextToVertex(CRay* I, float* v) {
  switch (I->context) {
    case pymol::RenderContext::UnitWindow: {
      float tw;
      float th;

      if (I->AspRatio > 1.0F) {
        tw = I->AspRatio;
        th = 1.0F;
      } else {
        th = 1.0F / I->AspRatio;
        tw = 1.0F;
      }

      if (!SettingGetGlobal_b(I->G, cSetting_ortho)) {
        float scale = v[2] + 0.5F;
        scale = I->FrontBackRatio * scale + 1.0F - scale;

        /* z-coodinate is easy... */

        v[2] = v[2] * I->Range[2] - (I->Volume[4] + I->Volume[5]) / 2.0F;
        v[0] -= 0.5F;
        v[1] -= 0.5F;
        v[0] = scale * v[0] * I->Range[0] / tw +
               (I->Volume[0] + I->Volume[1]) / 2.0F;
        v[1] = scale * v[1] * I->Range[1] / th +
               (I->Volume[2] + I->Volume[3]) / 2.0F;

        RayApplyMatrixInverse33(1, (float3*)v, I->ModelView, (float3*)v);
      } else {
        v[0] += (tw - 1.0F) / 2;
        v[1] += (th - 1.0F) / 2;
        v[0] = v[0] * (I->Range[0] / tw) + I->Volume[0];
        v[1] = v[1] * (I->Range[1] / th) + I->Volume[2];
        v[2] = v[2] * I->Range[2] - (I->Volume[4] + I->Volume[5]) / 2.0F;
        RayApplyMatrixInverse33(1, (float3*)v, I->ModelView, (float3*)v);
      }

    } break;
  }
}

static void RayApplyContextToNormal(CRay* I, float* v) {
  switch (I->context) {
    case pymol::RenderContext::UnitWindow:
      RayTransformInverseNormals33(1, (float3*)v, I->ModelView, (float3*)v);
      break;
  }
}
void CRay::wobble(int mode, const float* v) {
  CRay* I = this;
  I->Wobble = mode;
  if (v) copy3f(v, I->WobbleParam);
}
void CRay::transparentf(float v) {
  CRay* I = this;
  if (v > 1.0F) v = 1.0F;
  if (v < 0.0F) v = 0.0F;
  I->Trans = v;
}

void CRay::interiorColor3fv(const float* v, int passive) {
  CRay* I = this;
  I->IntColor[0] = (*v++);
  I->IntColor[1] = (*v++);
  I->IntColor[2] = (*v++);
  if (!passive) I->CheckInterior = true;
}
void CRay::color3fv(const float* v) {
  CRay* I = this;
  I->CurColor[0] = (*v++);
  I->CurColor[1] = (*v++);
  I->CurColor[2] = (*v++);
}
int CRay::sphere3fv(const float* v, float r) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimSphere;
  p->r1 = r;
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = (I->CurColor[0] < 0.0F);
  p->no_lighting = 0;

  I->PrimSize += 2 * r;
  I->PrimSizeCnt++;

  /*
     copy3f(I->WobbleParam,p->wobble_param); */
  vv = p->v1;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->c1;
  v = I->CurColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->ic;
  v = I->IntColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  if (I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
  }

  RayApplyContextToVertex(I, p->v1);

  I->NPrimitive++;
  return true;
}

void RayGetScaledAxes(CRay* I, float* xn, float* yn) {
  float* v;
  float vt[3];
  float xn0[3] = {1.0F, 0.0F, 0.0F};
  float yn0[3] = {0.0F, 1.0F, 0.0F};
  float v_scale;

  v = TextGetPos(I->G);

  if (I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), v, vt);
  } else {
    copy3f(v, vt);
  }

  v_scale = RayGetScreenVertexScale(I, vt) / I->Sampling;

  RayApplyMatrixInverse33(1, (float3*)xn0, glm::value_ptr(I->Rotation),
                          (float3*)xn0);
  RayApplyMatrixInverse33(1, (float3*)yn0, glm::value_ptr(I->Rotation),
                          (float3*)yn0);

  scale3f(xn0, v_scale, xn);
  scale3f(yn0, v_scale, yn);
}
int CRay::character(int char_id) {
  CRay* I = this;
  CPrimitive* p;
  float* v;
  float vt[3];
  float* vv;
  float width, height;
  float v_scale;
  int ok = true;

  v = TextGetPos(I->G);
  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive + 1);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCharacter;
  p->trans = I->Trans;
  p->char_id = char_id;
  p->wobble = I->Wobble;
  p->ramped = 0;
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = v[0];
  (*vv++) = v[1];
  (*vv++) = v[2];

  if (I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
  }
  /* what's the width of 1 screen window pixel at this point in space? */

  v_scale = RayGetScreenVertexScale(I, p->v1) / I->Sampling;

  RayApplyContextToVertex(I, p->v1);

  {
    float xn[3] = {1.0F, 0.0F, 0.0F};
    float yn[3] = {0.0F, 1.0F, 0.0F};
    float zn[3] = {0.0F, 0.0F, 1.0F};
    float sc[3];
    float scale;
    float xorig, yorig, advance;
    int width_i, height_i;
    CPrimitive* pp = p + 1;

    RayApplyMatrixInverse33(1, (float3*)xn, glm::value_ptr(I->Rotation),
                            (float3*)xn);
    RayApplyMatrixInverse33(1, (float3*)yn, glm::value_ptr(I->Rotation),
                            (float3*)yn);
    RayApplyMatrixInverse33(1, (float3*)zn, glm::value_ptr(I->Rotation),
                            (float3*)zn);

    CharacterGetGeometry(I->G, char_id, &width_i, &height_i, &xorig, &yorig,
                         &advance);
    width = (float)width_i;
    height = (float)height_i;

    scale = v_scale * advance;
    scale3f(xn, scale, vt); /* advance raster position in 3-space */
    add3f(v, vt, vt);
    TextSetPos(I->G, vt);

    /* position the pixmap relative to raster position */

    /*    scale = ((-xorig)-0.5F)*I->PixelRadius; */
    scale = ((-xorig) - 0.0F) * v_scale;
    scale3f(xn, scale, sc);
    add3f(sc, p->v1, p->v1);

    scale = ((-yorig) - 0.0F) * v_scale;
    scale3f(yn, scale, sc);
    add3f(sc, p->v1, p->v1);

    scale = v_scale * width;
    scale3f(xn, scale, xn);
    scale = v_scale * height;
    scale3f(yn, scale, yn);

    copy3f(zn, p->n0);
    copy3f(zn, p->n1);
    copy3f(zn, p->n2);
    copy3f(zn, p->n3);

    *(pp) = (*p);

    /* define coordinates of first triangle */

    add3f(p->v1, xn, p->v2);
    add3f(p->v1, yn, p->v3);

    I->PrimSize += 2 * (diff3f(p->v1, p->v2) + diff3f(p->v1, p->v3) +
                        diff3f(p->v2, p->v3));
    I->PrimSizeCnt += 6;

    /* encode characters coordinates in the colors  */

    zero3f(p->c1);
    set3f(p->c2, width, 0.0F, 0.0F);
    set3f(p->c3, 0.0F, height, 0.0F);

    /* define coordinates of second triangle */

    add3f(yn, xn, pp->v1);
    add3f(p->v1, pp->v1, pp->v1);
    add3f(p->v1, yn, pp->v2);
    add3f(p->v1, xn, pp->v3);

    {
      float *v, *vv;
      vv = p->ic;
      v = I->IntColor;
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      vv = pp->ic;
      v = I->IntColor;
      (*vv++) = (*v++);
      (*vv++) = (*v++);
      (*vv++) = (*v++);
    }

    /* encode integral character coordinates into the vertex colors  */

    set3f(pp->c1, width, height, 0.0F);
    set3f(pp->c2, 0.0F, height, 0.0F);
    set3f(pp->c3, width, 0.0F, 0.0F);
  }

  I->NPrimitive += 2;
  return true;
}
int CRay::cylinder3fv(const cgo::draw::cylinder& cyl) {
  return cylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                     cyl.color2, 1.0f - Trans, 1.0f - Trans);
}
int CRay::cylinder3fv(const cgo::draw::cylinder& cyl, const float alpha1,
                      const float alpha2) {
  return cylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                     cyl.color2, alpha1, alpha2);
}

int CRay::cylinder3fv(const float* v1, const float* v2, float r,
                      const float* c1, const float* c2, const float alpha1,
                      const float alpha2) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1 = r;
  p->cap1 = cCylCapFlat;
  p->cap2 = cCylCapFlat;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if (I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);

  // FIXME: alpha1 is not used
  p->trans = 1.0 - alpha2;
  {
    float* v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}
int CRay::customCylinder3fv(const cgo::draw::custom_cylinder& cyl,
                            const float alpha1, const float alpha2) {
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                           cyl.color2, cyl.get_cap1(), cyl.get_cap2(), alpha1,
                           alpha2);
}

int CRay::customCylinder3fv(const cgo::draw::custom_cylinder& cyl) {
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                           cyl.color2, cyl.get_cap1(), cyl.get_cap2(),
                           1.f - Trans, 1.f - Trans);
}

int CRay::customCylinderAlpha3fv(const cgo::draw::custom_cylinder_alpha& cyl) {
  return customCylinder3fv(cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                           cyl.color2, cyl.get_cap1(), cyl.get_cap2(),
                           cyl.color1[3], cyl.color2[3]);
}

int CRay::customCylinder3fv(const float* v1, const float* v2, float r,
                            const float* c1, const float* c2,
                            const cCylCap cap1, const cCylCap cap2) {
  return customCylinder3fv(v1, v2, r, c1, c2, cap1, cap2, 1.f - Trans,
                           1.f - Trans);
}

int CRay::customCylinder3fv(const float* v1, const float* v2, float r,
                            const float* c1, const float* c2,
                            const cCylCap cap1, const cCylCap cap2,
                            const float alpha1, const float alpha2) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCylinder;
  p->r1 = r;
  p->cap1 = cap1;
  p->cap2 = cap2;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if (I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;

  // FIXME: alpha1 is not used
  p->trans = 1.0f - alpha2;

  {
    float* v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}

int CRay::cone3fv(const float* v1, const float* v2, float r1, float r2,
                  const float* c1, const float* c2, cCylCap cap1,
                  cCylCap cap2) {
  CRay* I = this;
  CPrimitive* p;
  float r_max = (r1 > r2) ? r1 : r2;
  float* vv;
  int ok = true;

  if (r2 > r1) { /* make sure r1 is always larger */
    std::swap(r1, r2);
    std::swap(c1, c2);
    std::swap(v1, v2);
    std::swap(cap1, cap2);
  }

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimCone;
  p->r1 = r1;
  p->r2 = r2;
  p->trans = I->Trans;
  p->cap1 = cap1;

  if (cap2 >= cCylCapFlat) cap2 = cCylCapFlat;
  if (cap1 >= cCylCapFlat) cap1 = cCylCapFlat;

  p->cap2 = cap2;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r_max;
  I->PrimSizeCnt++;

  if (I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;
  {
    float* v;
    vv = p->ic;
    v = I->IntColor;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}
int CRay::sausage3fv(const float* v1, const float* v2, float r, const float* c1,
                     const float* c2) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimSausage;
  p->r1 = r;
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);

  I->PrimSize += diff3f(p->v1, p->v2) + 2 * r;
  I->PrimSizeCnt++;

  if (I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->ic;

  {
    float* v = I->IntColor;
    vv = p->ic;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  I->NPrimitive++;
  return true;
}

int CRay::ellipsoid3fv(const float* v, float r, const float* n1,
                       const float* n2, const float* n3) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;

  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimEllipsoid;
  p->r1 = r; /* maximum extent */
  p->trans = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = (I->CurColor[0] < 0.0F);
  p->no_lighting = 0;

  I->PrimSize += 2 * r;
  I->PrimSizeCnt++;

  vv = p->n0; /* storing lengths of the direction vectors in n0 */

  (*vv++) = length3f(n1);
  (*vv++) = length3f(n2);
  (*vv++) = length3f(n3);

  /* normalize the ellipsoid axes */

  vv = p->n1;
  if (p->n0[0] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[0];
    (*vv++) = (*n1++) * factor;
    (*vv++) = (*n1++) * factor;
    (*vv++) = (*n1++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->n2;
  if (p->n0[1] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[1];
    (*vv++) = (*n2++) * factor;
    (*vv++) = (*n2++) * factor;
    (*vv++) = (*n2++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->n3;
  if (p->n0[2] > R_SMALL8) {
    float factor;
    factor = 1.0F / p->n0[2];
    (*vv++) = (*n3++) * factor;
    (*vv++) = (*n3++) * factor;
    (*vv++) = (*n3++) * factor;
  } else {
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
    (*vv++) = 0.0F;
  }

  vv = p->v1;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->c1;
  v = I->CurColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  vv = p->ic;
  v = I->IntColor;
  (*vv++) = (*v++);
  (*vv++) = (*v++);
  (*vv++) = (*v++);

  if (I->TTTFlag) {
    p->r1 *= length3f(glm::value_ptr(I->TTT));
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n1, p->n1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n2, p->n2);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n3, p->n3);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToNormal(I, p->n1);
  RayApplyContextToNormal(I, p->n2);
  RayApplyContextToNormal(I, p->n3);

  I->NPrimitive++;
  return true;
}

int CRay::setLastToNoLighting(char no_lighting) {
  CRay* I = this;
  CPrimitive* p;
  if (!I->NPrimitive) return false;
  p = I->Primitive + I->NPrimitive - 1;
  p->no_lighting = no_lighting;
  return true;
}
int CRay::triangle3fv(const float* v1, const float* v2, const float* v3,
                      const float* n1, const float* n2, const float* n3,
                      const float* c1, const float* c2, const float* c3) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  float* vv;
  float n0[3] = {0.f, 0.f, 1.f}, nx[3], s1[3], s2[3], s3[3];
  float l1, l2, l3;
  short normals_exist = n1 && n2 && n3;

  /*  dump3f(v1," v1");
     dump3f(v2," v2");
     dump3f(v3," v3");
     dump3f(n1," n1");
     dump3f(n2," n2");
     dump3f(n3," n3");
     dump3f(c1," c1");
     dump3f(c2," c2");
     dump3f(c3," c3"); */
  VLACheck2<CPrimitive>(I->Primitive, I->NPrimitive);
  CHECKOK(ok, I->Primitive);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive;

  p->type = cPrimTriangle;
  p->trans = I->Trans;
  p->tr[0] = I->Trans;
  p->tr[1] = I->Trans;
  p->tr[2] = I->Trans;
  p->wobble = I->Wobble;
  p->ramped = ((c1[0] < 0.0F) || (c2[0] < 0.0F) || (c3[0] < 0.0F));
  p->no_lighting = 0;
  /*
     copy3f(I->WobbleParam,p->wobble_param); */

  /* determine exact triangle normal */
  if (normals_exist) {
    add3f(n1, n2, nx);
    add3f(n3, nx, nx);
  }
  subtract3f(v1, v2, s1);
  subtract3f(v3, v2, s2);
  subtract3f(v1, v3, s3);
  cross_product3f(s1, s2, n0);
  if (normals_exist) {
    if ((fabs(n0[0]) < RAY_SMALL) && (fabs(n0[1]) < RAY_SMALL) &&
        (fabs(n0[2]) < RAY_SMALL)) {
      copy3f(nx, n0);
    } /* fall-back */
    else if (dot_product3f(n0, nx) < 0)
      invert3f(n0);
  }
  normalize3f(n0);

  vv = p->n0;
  (*vv++) = n0[0];
  (*vv++) = n0[1];
  (*vv++) = n0[2];

  /* determine maximum distance from vertex to point */
  l1 = (float)length3f(s1);
  l2 = (float)length3f(s2);
  l3 = (float)length3f(s3);
  if (l2 > l1) {
    if (l3 > l2)
      l1 = l3;
    else
      l1 = l2;
  }
  /* store cutoff distance */

  p->r1 = l1 * 0.6F;

  /*  if(l1>20) {
     printf("%8.3f\n",l1);
     printf("%8.3f %8.3f %8.3f\n",s1[0],s1[1],s1[2]);
     printf("%8.3f %8.3f %8.3f\n",s2[0],s2[1],s2[2]);
     printf("%8.3f %8.3f %8.3f\n",s3[0],s3[1],s3[2]);
     } */

  vv = p->v1;
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  (*vv++) = (*v1++);
  vv = p->v2;
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  (*vv++) = (*v2++);
  vv = p->v3;
  (*vv++) = (*v3++);
  (*vv++) = (*v3++);
  (*vv++) = (*v3++);

  I->PrimSize +=
      diff3f(p->v1, p->v2) + diff3f(p->v1, p->v3) + diff3f(p->v2, p->v3);
  I->PrimSizeCnt += 3;

  vv = p->c1;
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  (*vv++) = (*c1++);
  vv = p->c2;
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  (*vv++) = (*c2++);
  vv = p->c3;
  (*vv++) = (*c3++);
  (*vv++) = (*c3++);
  (*vv++) = (*c3++);

  {
    float* v = I->IntColor;
    vv = p->ic;
    (*vv++) = (*v++);
    (*vv++) = (*v++);
    (*vv++) = (*v++);
  }

  if (normals_exist) {
    vv = p->n1;
    (*vv++) = (*n1++);
    (*vv++) = (*n1++);
    (*vv++) = (*n1++);
    vv = p->n2;
    (*vv++) = (*n2++);
    (*vv++) = (*n2++);
    (*vv++) = (*n2++);
    vv = p->n3;
    (*vv++) = (*n3++);
    (*vv++) = (*n3++);
    (*vv++) = (*n3++);
  } else {
    vv = p->n1;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
    vv = p->n2;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
    vv = p->n3;
    (*vv++) = n0[0];
    (*vv++) = n0[1];
    (*vv++) = n0[2];
  }

  if (I->TTTFlag) {
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v1, p->v1);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v2, p->v2);
    transformTTT44f3f(glm::value_ptr(I->TTT), p->v3, p->v3);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n0, p->n0);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n1, p->n1);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n2, p->n2);
    transform_normalTTT44f3f(glm::value_ptr(I->TTT), p->n3, p->n3);
  }

  RayApplyContextToVertex(I, p->v1);
  RayApplyContextToVertex(I, p->v2);
  RayApplyContextToVertex(I, p->v3);
  RayApplyContextToNormal(I, p->n0);
  RayApplyContextToNormal(I, p->n1);
  RayApplyContextToNormal(I, p->n2);
  RayApplyContextToNormal(I, p->n3);

  I->NPrimitive++;
  return true;
}

int CRay::triangleTrans3fv(const float* v1, const float* v2, const float* v3,
                           const float* n1, const float* n2, const float* n3,
                           const float* c1, const float* c2, const float* c3,
                           float t1, float t2, float t3) {
  CRay* I = this;
  CPrimitive* p;
  int ok = true;
  ok = I->triangle3fv(v1, v2, v3, n1, n2, n3, c1, c2, c3);
  if (!ok) return false;
  p = I->Primitive + I->NPrimitive - 1;

  p->tr[0] = t1;
  p->tr[1] = t2;
  p->tr[2] = t3;
  p->trans = (t1 + t2 + t3) / 3.0F;
  return true;
}
