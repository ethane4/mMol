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
-*
-*
Z* -------------------------------------------------------------------
*/

// Intersecting rays with spheres, cylinders, cones and ellipsoids.

#include "pymol/base/base.h"
#include "pymol/ray/basis.h"
#include "pymol/ray/basis_internal.h"

static const float kR_SMALL5 = 0.0001F;
int ZLineToSphere(float* base, float* point, float* dir, float radius,
                  float maxial, float* sphere, float* asum, float* pre) {
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float intra_p0, intra_p1, intra_p2;
  float radial, axial, axial_sum, dangle, ab_dangle, axial_perp;
  float radialsq, tan_acos_dangle;
  float vradial0, vradial1, vradial2;
  const float _0 = 0.0f;
  float point0 = point[0];
  float point1 = point[1];
  float point2 = point[2];
  float perpAxis0 =
      pre[0]; /* was cross_product(MinusZ,dir,perpAxis),normalize */
  float perpAxis1 = pre[1];
  float intra0 = point0 - base[0];
  float intra1 = point1 - base[1];
  const float dir0 = dir[0];
  const float dir1 = dir[1];
  const float dir2 = dir[2];
  float dot, perpDist;

  /* the perpAxis defines a perp-plane which includes the cyl-axis */

  /* get minimum distance between the lines */

  perpDist = intra0 * perpAxis0 + intra1 * perpAxis1;
  /* was dot_product3f(intra,perpAxis); */

  if ((perpDist < _0 ? -perpDist : perpDist) > radius) return 0;

  /*if(fabs(perpDist) > radius)   return 0; */

  dangle = -dir2; /* was dot(MinusZ,dir) */
  ab_dangle = (float)fabs(dangle);
  if (ab_dangle > (1.0f - kR_SMALL4)) {
    if (dangle > _0) {
      sphere[0] = point0;
      sphere[1] = point1;
      sphere[2] = point2;
    } else {
      sphere[0] = dir0 * maxial + point0;
      sphere[1] = dir1 * maxial + point1;
      sphere[2] = dir2 * maxial + point2;
    }
    return (1);
  }

  if (ab_dangle > kR_SMALL4)
    tan_acos_dangle = (float)(sqrt1d(1.0 - dangle * dangle) / dangle);
  else
    tan_acos_dangle = FLT_MAX;

  /* now we need to define the triangle in the perp-plane
     to figure out where the projected line intersection point is */

  intra_p2 = point2 - base[2];

  /* first, compute radial distance in the perp-plane between the two starting
   * points */

  dot = intra0 * perpAxis0 + intra1 * perpAxis1;
  intra_p0 = intra0 - perpAxis0 * dot;
  intra_p1 = intra1 - perpAxis1 * dot;

  dot = intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2;
  vradial0 = intra_p0 - dir0 * dot;
  vradial1 = intra_p1 - dir1 * dot;
  vradial2 = intra_p2 - dir2 * dot;

  radialsq =
      ((vradial0 * vradial0) + (vradial1 * vradial1) + (vradial2 * vradial2));

  /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

  if (ab_dangle < kR_SMALL4)
    axial_perp = _0;
  else
    axial_perp = (float)(sqrt1f(radialsq) / tan_acos_dangle);

  axial = (float)sqrt1f(
      ((intra_p0 * intra_p0) + (intra_p1 * intra_p1) + (intra_p2 * intra_p2)) -
      radialsq);

  if ((intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2) >= _0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;

  /* now we have to think about where the vector will actually strike the
   * cylinder */

  /* by definition, it must be perpdist away from the perp-plane becuase the
     perp-plane is parallel to the line, so we can compute the radial component
     to this point */

  radial = radius * radius - perpDist * perpDist;
  radial = (float)sqrt1f(radial);

  /* now the trick is figuring out how to adjust the axial distance to get the
     actual position along the cyl line which will give us a representative
     sphere */

  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;
  /*
     printf("radial2 %8.3f \n",radial); */

  if (axial_sum < _0)
    axial_sum = _0;
  else if (axial_sum > maxial)
    axial_sum = maxial;

  sphere[0] = dir0 * axial_sum + point0;
  sphere[1] = dir1 * axial_sum + point1;
  sphere[2] = dir2 * axial_sum + point2;

  *asum = axial_sum;
  /*  printf("==>%8.3f sphere %8.3f %8.3f
   * %8.3f\n",base[1],sphere[1],axial_perp,axial); */
  return (1);
}

int LineToSphere(float* base, float* ray, float* point, float* dir,
                 float radius, float maxial, float* sphere, float* asum) {
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float perpDist, radial, axial, axial_sum, dangle, ab_dangle, axial_perp;
  float radialsq, tan_acos_dangle;
  float perpAxis0, perpAxis1, perpAxis2;
  float intra0, intra1, intra2;
  float intra_p0, intra_p1, intra_p2;
  float vradial0, vradial1, vradial2;
  float dir0 = dir[0], dir1 = dir[1], dir2 = dir[2];
  float ray0 = ray[0], ray1 = ray[1], ray2 = ray[2];
  float point0 = point[0], point1 = point[1], point2 = point[2];
  float dot;

  const float _0 = 0.0F;
  const float _1 = 1.0F;

  /*    cross_product3f(ray,dir,perpAxis); */

  perpAxis0 = ray1 * dir2 - ray2 * dir1;
  perpAxis1 = ray2 * dir0 - ray0 * dir2;
  perpAxis2 = ray0 * dir1 - ray1 * dir0;

  /* subtract3f(point,base,intra); */
  intra0 = point0 - base[0];

  /* normalize3f(perpAxis) */
  {
    float len =
        (float)sqrt1d((perpAxis0 * perpAxis0) + (perpAxis1 * perpAxis1) +
                      (perpAxis2 * perpAxis2));
    intra1 = point1 - base[1];
    /* subtract3f(point,base,intra); */
    if (len > R_SMALL8) {
      len = _1 / len;
      intra2 = point2 - base[2];
      perpAxis0 *= len;
      perpAxis1 *= len;
      perpAxis2 *= len;
    } else {
      intra2 = point2 - base[2];
    }
  }
  /* the perpAxis defines a perp-plane which includes the cyl-axis */

  /* get minimum distance between the lines */

  /* dot_product3f(intra, perpAxis); */
  perpDist = intra0 * perpAxis0;
  perpDist += intra1 * perpAxis1 + intra2 * perpAxis2;

  /*if(fabs(perpDist) > radius)   return 0; */

  dangle = ray0 * dir0;

  if ((perpDist < _0 ? -perpDist : perpDist) > radius) return 0;

  /* dangle  = dot_product3f(ray, dir); */
  dangle += ray1 * dir1 + ray2 * dir2;

  ab_dangle = (float)fabs(dangle);

  if (ab_dangle > (1.0f - kR_SMALL4)) {
    if (dangle > _0) {
      sphere[0] = point0;
      sphere[1] = point1;
      sphere[2] = point2;
    } else {
      sphere[0] = dir0 * maxial + point0;
      sphere[1] = dir1 * maxial + point1;
      sphere[2] = dir2 * maxial + point2;
    }
    return (1);
  }

  if (ab_dangle > kR_SMALL4)
    tan_acos_dangle = (float)(sqrt1d(1.0 - dangle * dangle) / dangle);
  else
    tan_acos_dangle = FLT_MAX;

  /* now we need to define the triangle in the perp-plane
     to figure out where the projected line intersection point is */

  /* first, compute radial distance in the perp-plane between the two starting
   * points */

  /* dot = dot_product3f(intra,perpAxis); */
  dot = intra0 * perpAxis0 + intra1 * perpAxis1 + intra2 * perpAxis2;

  intra_p0 = intra0 - perpAxis0 * dot;
  intra_p1 = intra1 - perpAxis1 * dot;
  intra_p2 = intra2 - perpAxis2 * dot;

  /* dot = dot_product3f(intra_p, dir); */
  dot = intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2;

  vradial0 = intra_p0 - dir0 * dot;
  vradial1 = intra_p1 - dir1 * dot;
  vradial2 = intra_p2 - dir2 * dot;

  radialsq =
      ((vradial0 * vradial0) + (vradial1 * vradial1) + (vradial2 * vradial2));

  /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

  if (ab_dangle < kR_SMALL4)
    axial_perp = _0;
  else
    axial_perp = (float)(sqrt1f(radialsq) / tan_acos_dangle);

  axial = (float)sqrt1f(
      ((intra_p0 * intra_p0) + (intra_p1 * intra_p1) + (intra_p2 * intra_p2)) -
      radialsq);

  if ((intra_p0 * dir0 + intra_p1 * dir1 + intra_p2 * dir2) >= _0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;

  /* now we have to think about where the vector will actually strike the
   * cylinder */

  /* by definition, it must be perpdist away from the perp-plane becuase the
     perp-plane is parallel to the line, so we can compute the radial component
     to this point */

  radial = radius * radius - perpDist * perpDist;
  radial = (float)sqrt1f(radial);

  /* now the trick is figuring out how to adjust the axial distance to get the
     actual position along the cyl line which will give us a representative
     sphere */

  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;
  /*
     printf("radial2 %8.3f \n",radial); */

  if (axial_sum < _0)
    axial_sum = _0;
  else if (axial_sum > maxial)
    axial_sum = maxial;

  sphere[0] = dir0 * axial_sum + point0;
  sphere[1] = dir1 * axial_sum + point1;
  sphere[2] = dir2 * axial_sum + point2;

  *asum = axial_sum;
  /*  printf("==>%8.3f sphere %8.3f %8.3f
   * %8.3f\n",base[1],sphere[1],axial_perp,axial); */
  return (1);
}

int FrontToInteriorSphere(float* front, float* point, float* dir, float radius,
                          float radius2, float maxial) {
  float intra_p[3];
  float axial;
  float intra[3], axis[3];
  float sphere[3];

  subtract3f(point, front, intra);
  remove_component3f(intra, dir, intra_p);
  add3f(front, intra_p, intra_p);
  subtract3f(point, intra_p, axis);
  axial = -dot_product3f(axis, dir);

  if (axial < 0.0F)
    axial = 0.0F;
  else if (axial > maxial)
    axial = maxial;

  sphere[0] = axial * dir[0] + point[0];
  sphere[1] = axial * dir[1] + point[1];
  sphere[2] = axial * dir[2] + point[2];

  return (diffsq3f(sphere, front) <= radius2);
}
int ZLineToSphereCapped(float* base, float* point, float* dir, float radius,
                        float maxial, float* sphere, float* asum, cCylCap cap1,
                        cCylCap cap2, float* pre) {
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float perpAxis[3], intra_p[3];
  float perpDist, radial, axial, axial_sum, dangle, ab_dangle, axial_perp;
  float radialsq, tan_acos_dangle;
  float len_proj;
  float intra[3], vradial[3];
  float diff[3], fpoint[3];
  float proj[3];

  perpAxis[0] = pre[0]; /* was cross_product(MinusZ,dir,perpAxis),normalize */
  perpAxis[1] = pre[1];

  /* the perpAxis defines a perp-plane which includes the cyl-axis */

  intra[0] = point[0] - base[0];
  intra[1] = point[1] - base[1];

  /* get minimum distance between the lines */

  perpDist = intra[0] * perpAxis[0] + intra[1] * perpAxis[1];
  /* was dot_product3f(intra,perpAxis); */

  if (fabs(perpDist) > radius) {
    return (0);
  }

  perpAxis[2] = 0.0;
  intra[2] = point[2] - base[2];

  dangle = -dir[2]; /* was dot(MinusZ,dir) */
  ab_dangle = (float)fabs(dangle);
  if (ab_dangle > (1 - kR_SMALL4)) { /* vector inline with light ray... */
    vradial[0] = point[0] - base[0];
    vradial[1] = point[1] - base[1];
    vradial[2] = 0.0F;
    radial = (float)length3f(vradial);
    if (radial > radius) return 0;
    if (dangle > 0.0) {
      switch (cap1) {
        case cCylCapFlat:
          sphere[0] = base[0];
          sphere[1] = base[1];
          sphere[2] = point[2] - radius;
          break;
        case cCylCapRound:
          sphere[0] = point[0];
          sphere[1] = point[1];
          sphere[2] = point[2];
      }
      return (1);
    } else {
      switch (cap1) {
        case cCylCapFlat:
          sphere[0] = base[0];
          sphere[1] = base[1];
          sphere[2] = dir[2] * maxial + point[2] - radius;
          break;
        case cCylCapRound:
          sphere[0] = dir[0] * maxial + point[0];
          sphere[1] = dir[1] * maxial + point[1];
          sphere[2] = dir[2] * maxial + point[2];
      }
      return (1);
    }
  }

  /*tan_acos_dangle = tan(acos(dangle)); */
  tan_acos_dangle = (float)sqrt1f(1 - dangle * dangle) / dangle;

  /*
     printf("perpDist %8.3f\n",perpDist);
     printf("dir %8.3f %8.3f %8.3f\n",dir[0],dir[1],dir[2]);
     printf("base %8.3f %8.3f %8.3f\n",base[0],base[1],base[2]);
     printf("point %8.3f %8.3f %8.3f\n",point[0],point[1],point[2]);
     printf("intra %8.3f %8.3f %8.3f\n",intra[0],intra[1],intra[2]);
     printf("perpAxis %8.3f %8.3f %8.3f\n",perpAxis[0],perpAxis[1],perpAxis[2]);
   */

  /* now we need to define the triangle in the perp-plane
     to figure out where the projected line intersection point is */

  /* first, compute radial distance in the perp-plane between the two starting
   * points */

  remove_component3f(intra, perpAxis, intra_p);
  remove_component3f(intra_p, dir, vradial);

  radialsq = lengthsq3f(vradial);

  /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

  if (ab_dangle < kR_SMALL4)
    axial_perp = 0;
  else
    axial_perp = (float)sqrt1f(radialsq) / tan_acos_dangle;

  axial = (float)lengthsq3f(intra_p) - radialsq;
  axial = (float)sqrt1f(axial);

  /*
     printf("radial %8.3f\n",radial);
     printf("vradial %8.3f %8.3f %8.3f\n",vradial[0],vradial[1],vradial[2]);
     printf("radial %8.3f\n",radial);
     printf("dangle %8.3f \n",dangle);
     printf("axial_perp %8.3f \n",axial_perp);
     printf("axial1 %8.3f \n",axial);
     printf("%8.3f\n",dot_product3f(intra_p,dir));
   */

  if (dot_product3f(intra_p, dir) >= 0.0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;

  /*
     printf("axial2 %8.3f\n",axial);
   */

  /* now we have to think about where the vector will actually strike the
   * cylinder */

  /* by definition, it must be perpdist away from the perp-plane becuase the
     perp-plane is parallel to the line, so we can compute the radial component
     to this point */

  radial = radius * radius - perpDist * perpDist;
  radial = (float)sqrt1f(radial);

  /* now the trick is figuring out how to adjust the axial distance to get the
     actual position along the cyl line which will give us a representative
     sphere */

  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;

  /*    printf("ab_dangle %8.3f \n",ab_dangle);

     printf("axial_sum %8.3f \n",axial_sum);
   */
  if (axial_sum < 0) {
    switch (cap1) {
      case cCylCapFlat:
        subtract3f(point, base, diff);
        project3f(diff, dir, proj);
        len_proj = (float)length3f(proj);
        dangle = -proj[2] / len_proj;
        if (fabs(dangle) < kR_SMALL4) return 0;
        sphere[0] = base[0];
        sphere[1] = base[1];
        sphere[2] = base[2] - len_proj / dangle;
        if (diff3f(sphere, point) > radius) return 0;
        sphere[0] += dir[0] * radius;
        sphere[1] += dir[1] * radius;
        sphere[2] += dir[2] * radius;
        *asum = 0;
        break;
      case cCylCapRound:
        axial_sum = 0;
        sphere[0] = point[0];
        sphere[1] = point[1];
        sphere[2] = point[2];
        *asum = axial_sum;
        break;
      case cCylCapNone:
      default:
        return 0;
        break;
    }
  } else if (axial_sum > maxial) {
    switch (cap2) {
      case cCylCapFlat:
        scale3f(dir, maxial, fpoint);
        add3f(fpoint, point, fpoint);
        subtract3f(fpoint, base, diff);
        project3f(diff, dir, proj);
        len_proj = (float)length3f(proj);
        dangle = -proj[2] / len_proj;
        if (fabs(dangle) < kR_SMALL4) return 0;
        sphere[0] = base[0];
        sphere[1] = base[1];
        sphere[2] = base[2] - len_proj / dangle;
        if (diff3f(sphere, fpoint) > radius) return 0;
        sphere[0] -= dir[0] * radius;
        sphere[1] -= dir[1] * radius;
        sphere[2] -= dir[2] * radius;
        *asum = maxial;
        break;
      case cCylCapRound:
        axial_sum = maxial;
        sphere[0] = dir[0] * axial_sum + point[0];
        sphere[1] = dir[1] * axial_sum + point[1];
        sphere[2] = dir[2] * axial_sum + point[2];
        *asum = axial_sum;
        break;
      case cCylCapNone:
      default:
        return 0;
        break;
    }
  } else {
    sphere[0] = dir[0] * axial_sum + point[0];
    sphere[1] = dir[1] * axial_sum + point[1];
    sphere[2] = dir[2] * axial_sum + point[2];

    *asum = axial_sum;

    /*  printf("==>%8.3f sphere %8.3f %8.3f
     * %8.3f\n",base[1],sphere[1],axial_perp,axial); */
  }
  return (1);
}

int LineToSphereCapped(float* base, float* ray, float* point, float* dir,
                       float radius, float maxial, float* sphere, float* asum,
                       cCylCap cap1, cCylCap cap2) {
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float perpAxis[3], intra_p[3];
  float perpDist, radial, axial, axial_sum, dangle, ab_dangle, axial_perp;
  float radialsq, tan_acos_dangle;
  float len_proj;
  float intra[3], vradial[3];
  float diff[3], fpoint[3];
  float proj[3];

  subtract3f(point, base, intra);

  cross_product3f(ray, dir, perpAxis);

  normalize3f(perpAxis);

  /* the perpAxis defines a perp-plane which includes the cyl-axis */

  /* get minimum distance between the lines */

  perpDist = dot_product3f(intra, perpAxis);

  /* was dot_product3f(intra,perpAxis); */

  if (fabs(perpDist) > radius) {
    return (0);
  }

  dangle = dot_product3f(ray, dir);
  ab_dangle = (float)fabs(dangle);

  if (ab_dangle > (1 - kR_SMALL4)) { /* vector inline with light ray... */
    vradial[0] = point[0] - base[0];
    vradial[1] = point[1] - base[1];
    vradial[2] = point[2] - base[2];
    radial = (float)length3f(vradial);
    if (radial > radius) return 0;
    if (dangle > 0.0) {
      switch (cap1) {
        case cCylCapFlat:
          sphere[0] = dir[0] * radius + point[0];
          sphere[1] = dir[1] * radius + point[1];
          sphere[2] = dir[2] * radius + point[2];
          break;
        case cCylCapRound:
          sphere[0] = point[0];
          sphere[1] = point[1];
          sphere[2] = point[2];
      }
      return (1);
    } else {
      switch (cap1) {
        case cCylCapFlat:
          maxial -= radius;
          break;
      }
      sphere[0] = dir[0] * maxial + point[0];
      sphere[1] = dir[1] * maxial + point[1];
      sphere[2] = dir[2] * maxial + point[2];
      return (1);
    }
  }

  /*tan_acos_dangle = tan(acos(dangle)); */
  tan_acos_dangle = (float)sqrt1f(1 - dangle * dangle) / dangle;

  /*
     printf("perpDist %8.3f\n",perpDist);
     printf("dir %8.3f %8.3f %8.3f\n",dir[0],dir[1],dir[2]);
     printf("base %8.3f %8.3f %8.3f\n",base[0],base[1],base[2]);
     printf("point %8.3f %8.3f %8.3f\n",point[0],point[1],point[2]);
     printf("intra %8.3f %8.3f %8.3f\n",intra[0],intra[1],intra[2]);
     printf("perpAxis %8.3f %8.3f %8.3f\n",perpAxis[0],perpAxis[1],perpAxis[2]);
   */

  /* now we need to define the triangle in the perp-plane
     to figure out where the projected line intersection point is */

  /* first, compute radial distance in the perp-plane between the two starting
   * points */

  remove_component3f(intra, perpAxis, intra_p);
  remove_component3f(intra_p, dir, vradial);

  radialsq = lengthsq3f(vradial);

  /* now figure out the axial distance along the cyl-line that will give us
     the point of closest approach */

  if (ab_dangle < kR_SMALL4)
    axial_perp = 0;
  else
    axial_perp = (float)sqrt1f(radialsq) / tan_acos_dangle;

  axial = (float)lengthsq3f(intra_p) - radialsq;
  axial = (float)sqrt1f(axial);

  /*
     printf("radial %8.3f\n",radial);
     printf("vradial %8.3f %8.3f %8.3f\n",vradial[0],vradial[1],vradial[2]);
     printf("radial %8.3f\n",radial);
     printf("dangle %8.3f \n",dangle);
     printf("axial_perp %8.3f \n",axial_perp);
     printf("axial1 %8.3f \n",axial);
     printf("%8.3f\n",dot_product3f(intra_p,dir));
   */

  if (dot_product3f(intra_p, dir) >= 0.0)
    axial = axial_perp - axial;
  else
    axial = axial_perp + axial;

  /*
     printf("axial2 %8.3f\n",axial);
   */

  /* now we have to think about where the vector will actually strike the
   * cylinder */

  /* by definition, it must be perpdist away from the perp-plane becuase the
     perp-plane is parallel to the line, so we can compute the radial component
     to this point */

  radial = radius * radius - perpDist * perpDist;
  radial = (float)sqrt1f(radial);

  /* now the trick is figuring out how to adjust the axial distance to get the
     actual position along the cyl line which will give us a representative
     sphere */

  if (ab_dangle > kR_SMALL4)
    axial_sum = axial - radial / tan_acos_dangle;
  else
    axial_sum = axial;

  /*    printf("ab_dangle %8.3f \n",ab_dangle);

     printf("axial_sum %8.3f \n",axial_sum);
   */
  if (axial_sum < 0) {
    switch (cap1) {
      case cCylCapFlat:
        subtract3f(point, base, diff);
        project3f(diff, dir, proj);
        len_proj = (float)length3f(proj);
        dangle = dot_product3f(proj, ray) / len_proj;
        if (fabs(dangle) < kR_SMALL4) return 0;
        len_proj /= dangle;
        sphere[0] = base[0] + ray[0] * len_proj;
        sphere[1] = base[1] + ray[1] * len_proj;
        sphere[2] = base[2] + ray[2] * len_proj;
        if (diff3f(sphere, point) > radius) return 0;
        sphere[0] += dir[0] * radius;
        sphere[1] += dir[1] * radius;
        sphere[2] += dir[2] * radius;
        *asum = 0;
        break;
      case cCylCapRound:
        axial_sum = 0;
        sphere[0] = point[0];
        sphere[1] = point[1];
        sphere[2] = point[2];
        *asum = axial_sum;
        break;
      case cCylCapNone:
      default:
        return 0;
        break;
    }
  } else if (axial_sum > maxial) {
    switch (cap2) {
      case cCylCapFlat:
        scale3f(dir, maxial, fpoint);
        add3f(fpoint, point, fpoint);
        subtract3f(fpoint, base, diff);
        project3f(diff, dir, proj);
        len_proj = (float)length3f(proj);
        dangle = dot_product3f(proj, ray) / len_proj;
        if (fabs(dangle) < kR_SMALL4) return 0;
        len_proj /= dangle;
        sphere[0] = base[0] + ray[0] * len_proj;
        sphere[1] = base[1] + ray[1] * len_proj;
        sphere[2] = base[2] + ray[2] * len_proj;
        if (diff3f(sphere, fpoint) > radius) return 0;
        sphere[0] -= dir[0] * radius;
        sphere[1] -= dir[1] * radius;
        sphere[2] -= dir[2] * radius;
        *asum = maxial;
        break;
      case cCylCapRound:
        axial_sum = maxial;
        sphere[0] = dir[0] * axial_sum + point[0];
        sphere[1] = dir[1] * axial_sum + point[1];
        sphere[2] = dir[2] * axial_sum + point[2];
        *asum = axial_sum;
        break;
      case cCylCapNone:
      default:
        return 0;
        break;
    }
  } else {
    sphere[0] = dir[0] * axial_sum + point[0];
    sphere[1] = dir[1] * axial_sum + point[1];
    sphere[2] = dir[2] * axial_sum + point[2];

    *asum = axial_sum;

    /*  printf("==>%8.3f sphere %8.3f %8.3f
     * %8.3f\n",base[1],sphere[1],axial_perp,axial); */
  }
  return (1);
}

int ConeLineToSphereCapped(float* base, float* ray, float* point, float* dir,
                           float radius, float small_radius, float maxial,
                           float* sphere, float* asum, float* sph_rad,
                           float* sph_rad_sq, cCylCap cap1, cCylCap cap2) {
  /* Strategy - find an imaginary sphere that lies at the correct point on
     the line segment, then treat as a sphere reflection */

  float axial_sum, dangle, ab_dangle;
  float len_proj;
  float diff[3], fpoint[3];
  float proj[3];

  {
    float perp_axis[3];
    float intra[3];
    float perp_dist;

    subtract3f(point, base, intra);
    cross_product3f(ray, dir, perp_axis);

    normalize3f(perp_axis);

    /* the perp_axis defines a perp-plane which includes the cyl-axis */

    /* get minimum distance between the lines */

    perp_dist = fabs(dot_product3f(intra, perp_axis));

    if (perp_dist > radius) {
      /* the infinite ray and the cone direction lines don't pass close
         enough to intersect within the bounding cylinder */
      return 0;
    }
  }

  dangle = dot_product3f(ray, dir);
  ab_dangle = (float)fabs(dangle);

  /* set up the cone */

  {
    double spread = (radius - small_radius) / maxial;
    float orig_axial_len = radius / spread;
    float base2orig_radial[3];
    float base2orig_normal[3];
    float base2orig_radial_len, base2orig_radial_len_sq;
    float base2orig_len_sq;
    float base2orig_axial_len, base2orig_spread;
    float orig[3];
    float base2orig[3];
    float near_p[3];

    int base_inside_cone = false;
    float shift1, shift2; /* this is what we are solving for --
                             the distance along the cone axis to the point of
                             intersection */

    scale3f(dir, orig_axial_len, orig);
    add3f(point, orig, orig);
    subtract3f(orig, base, base2orig);

    remove_component3f(base2orig, dir, base2orig_radial);

    base2orig_radial_len_sq = lengthsq3f(base2orig_radial);

    base2orig_len_sq = lengthsq3f(base2orig);

    base2orig_axial_len = sqrt1f(base2orig_len_sq - base2orig_radial_len_sq);

    base2orig_radial_len = sqrt1f(base2orig_radial_len_sq);

    base2orig_spread = base2orig_radial_len / base2orig_axial_len;

    base_inside_cone = (base2orig_spread < spread);

    normalize23f(base2orig, base2orig_normal);

    if (ab_dangle > kR_SMALL4) {
      float ray_extend = base2orig_axial_len / dangle;
      if (dot_product3f(base2orig_normal, dir) < 0.0) ray_extend = -ray_extend;

      scale3f(ray, ray_extend, near_p);
      add3f(base, near_p, near_p);

      /* Now we punt entirely and throw the solution of this quadratic
         relationship over to Mathematica.  Surely this calculation
         could be significantly optimized... */

      {
        double partA, partB, partC;

        double dir0 = dir[0], dir1 = dir[1], dir2 = dir[2];
        double ray0 = ray[0], ray1 = ray[1], ray2 = ray[2];
        double dir0Sq = dir0 * dir0, dir1Sq = dir1 * dir1, dir2Sq = dir2 * dir2;
        double ray0Sq = ray0 * ray0, ray1Sq = ray1 * ray1, ray2Sq = ray2 * ray2;

        double cone0 = orig[0], cone1 = orig[1], cone2 = orig[2];
        double near0 = near_p[0], near1 = near_p[1], near2 = near_p[2];
        double cone0Sq = cone0 * cone0, cone1Sq = cone1 * cone1,
               cone2Sq = cone2 * cone2;
        double near0Sq = near0 * near0, near1Sq = near1 * near1,
               near2Sq = near2 * near2;

        double dAngle = ray0 * dir0 + ray1 * dir1 + ray2 * dir2;

        double spreadSq = spread * spread;
        double dAngleSq = dAngle * dAngle;

        partB =
            dAngleSq *
            (4 * pow(cone0 * dAngle * dir0 + cone1 * dAngle * dir1 +
                         cone2 * dAngle * dir2 - dAngle * dir0 * near0 -
                         dAngle * dir1 * near1 - dAngle * dir2 * near2 -
                         cone0 * ray0 + near0 * ray0 - cone1 * ray1 +
                         near1 * ray1 - cone2 * ray2 + near2 * ray2,
                     2.0) -
             4 *
                 (cone0Sq + cone1Sq + cone2Sq - 2 * cone0 * near0 + near0Sq -
                  2 * cone1 * near1 + near1Sq - 2 * cone2 * near2 + near2Sq) *
                 (ray0Sq + ray1Sq -
                  2 * dAngle * (dir0 * ray0 + dir1 * ray1 + dir2 * ray2) +
                  ray2Sq + dAngleSq * (dir0Sq + dir1Sq + dir2Sq - spreadSq)));

        if (partB < 0.0) { /* negative? then there are NO real solutions */
          return 0;
        } else {
          double partBroot = sqrt(partB);

          partA = -(cone0 * dAngleSq * dir0) - cone1 * dAngleSq * dir1 -
                  cone2 * dAngleSq * dir2 + dAngleSq * dir0 * near0 +
                  dAngleSq * dir1 * near1 + dAngleSq * dir2 * near2 +
                  cone0 * dAngle * ray0 - dAngle * near0 * ray0 +
                  cone1 * dAngle * ray1 - dAngle * near1 * ray1 +
                  cone2 * dAngle * ray2 - dAngle * near2 * ray2;

          partC = (ray0Sq + ray1Sq -
                   2 * dAngle * (dir0 * ray0 + dir1 * ray1 + dir2 * ray2) +
                   ray2Sq + dAngleSq * (dir0Sq + dir1Sq + dir2Sq - spreadSq));

          shift1 = (float)((partA + partBroot * 0.5) / partC);
          shift2 = (float)((partA - partBroot * 0.5) / partC);
        }
      }

      {
        float axial_sum1 = orig_axial_len + shift1;
        float axial_sum2 = orig_axial_len + shift2;

        if (dangle > 0.0F) { /* cone is narrowing in parallel with ray */
          if (shift1 < shift2) {
            axial_sum = axial_sum1;
          } else {
            axial_sum = axial_sum2;
          }
          if ((axial_sum < 0.0F) ||
              (base_inside_cone && (axial_sum < orig_axial_len))) {
            switch (cap1) {
              case cCylCapFlat:
                subtract3f(point, base, diff);
                project3f(diff, dir, proj);
                len_proj = (float)length3f(proj);
                dangle = dot_product3f(proj, ray) / len_proj;
                if (fabs(dangle) < kR_SMALL5) return 0;
                len_proj /= dangle;
                sphere[0] = base[0] + ray[0] * len_proj;
                sphere[1] = base[1] + ray[1] * len_proj;
                sphere[2] = base[2] + ray[2] * len_proj;
                if (diff3f(sphere, point) > radius) return 0;
                sphere[0] += dir[0] * radius;
                sphere[1] += dir[1] * radius;
                sphere[2] += dir[2] * radius;
                *sph_rad = radius;
                *sph_rad_sq = radius * radius;
                *asum = 0;
                return 1;
                break;
              case cCylCapNone:
              default:
                return 0;
                break;
            }
          } else if (axial_sum > maxial) {
            return 0;
          }
        } else { /* cone is narrowing against ray */
          if (shift1 < shift2) {
            axial_sum = axial_sum2;
            if (axial_sum > orig_axial_len) {
              axial_sum = axial_sum1;
            }
          } else {
            axial_sum = axial_sum1;
            if (axial_sum > orig_axial_len) {
              axial_sum = axial_sum2;
            }
          }
          if (axial_sum < 0.0F) {
            return 0;
          } else if (axial_sum >= maxial) {
            switch (cap2) {
              case cCylCapFlat:
                scale3f(dir, maxial, fpoint);
                add3f(fpoint, point, fpoint);
                subtract3f(fpoint, base, diff);
                project3f(diff, dir, proj);
                len_proj = (float)length3f(proj);
                dangle = dot_product3f(proj, ray) / len_proj;
                if (fabs(dangle) < kR_SMALL5) return 0;
                len_proj /= dangle;
                sphere[0] = base[0] + ray[0] * len_proj;
                sphere[1] = base[1] + ray[1] * len_proj;
                sphere[2] = base[2] + ray[2] * len_proj;
                if (diff3f(sphere, fpoint) > small_radius)
                  /* need to handle this case */
                  return 0;
                sphere[0] -= dir[0] * small_radius;
                sphere[1] -= dir[1] * small_radius;
                sphere[2] -= dir[2] * small_radius;
                *sph_rad = small_radius;
                *sph_rad_sq = small_radius * small_radius;
                *asum = maxial;
                return 1;
                break;
              case cCylCapNone:
              default:
                return 0;
                break;
            }
          }
        }
      }
    } else {
      axial_sum = orig_axial_len - base2orig_axial_len;
      if ((axial_sum < 0.0F) || (axial_sum > maxial)) return 0;
    }

    /* normal hit in mid-section of the cone */

    {
      float radius_at_hit = radius - spread * axial_sum;

      float adjustment = radius_at_hit * spread;

      *asum = axial_sum; /* color blend based on actual hit location */

      axial_sum -= adjustment;

      /*        printf(" %8.3f ",axial_sum); */

      sphere[0] = dir[0] * axial_sum + point[0];
      sphere[1] = dir[1] * axial_sum + point[1];
      sphere[2] = dir[2] * axial_sum + point[2];

      /* and provide virtual sphere radius info */

      *sph_rad_sq =
          (float)((radius_at_hit * radius_at_hit) + (adjustment * adjustment));
      *sph_rad = (float)sqrt(*sph_rad_sq);
      /*        printf("%8.3f %8.3f ",*sph_rad, tan_acos_dangle); */

      /*         dump3f(sphere,"sph"); */
      return 1;
    }
  }
  return 0;
}

int FrontToInteriorSphereCapped(float* front, float* point, float* dir,
                                float radius, float radius2, float maxial,
                                cCylCap cap1, cCylCap cap2) {
  float intra_p[3];
  float axial;
  float intra[3], axis[3];
  float sphere[3];

  subtract3f(point, front, intra);
  remove_component3f(intra, dir, intra_p);
  add3f(front, intra_p, intra_p);
  subtract3f(point, intra_p, axis);
  axial = -dot_product3f(axis, dir);

  if (axial < 0.0F)
    return 0;
  else if (axial > maxial)
    return 0;

  sphere[0] = axial * dir[0] + point[0];
  sphere[1] = axial * dir[1] + point[1];
  sphere[2] = axial * dir[2] + point[2];

  return (diffsq3f(sphere, front) < radius2);
}
float ZLineClipPoint(float* base, float* point, float* alongNormalSq,
                     float cutoff) {
  float hyp0, hyp1, hyp2;
  float result = FLT_MAX;

  /* this routine determines whether or not a vector starting at "base"
     heading in the direction "normal" intersects a sphere located at "point".

     It returns how far along the vector the intersection with the plane is or
     MAXFLOAT if there isn't a relevant intersection

     NOTE: this routine has been optimized for normals along Z
     Optimizes-out vectors that are more than "cutoff" from "point" in x,y plane
   */

  hyp0 = point[0] - base[0];
  if (fabs(hyp0) > cutoff) return result;

  hyp1 = point[1] - base[1];
  if (fabs(hyp1) > cutoff) return result;

  hyp2 = point[2] - base[2];

  if (hyp2 < 0.0) {
    (*alongNormalSq) = (hyp2 * hyp2);
    result = (hyp0 * hyp0) + (hyp1 * hyp1);
  }
  return result;
}

float ZLineClipPointNoZCheck(float* base, float* point, float* alongNormalSq,
                             float cutoff) {
  float hyp0, hyp1, hyp2;
  float result = FLT_MAX;

  /* this routine determines whether or not a vector starting at "base"
     heading in the direction "normal" intersects a sphere located at "point".

     It returns how far along the vector the intersection with the plane is or
     MAXFLOAT if there isn't a relevant intersection

     NOTE: this routine has been optimized for normals along Z
     Optimizes-out vectors that are more than "cutoff" from "point" in x,y plane
   */

  hyp0 = point[0] - base[0];
  if (fabs(hyp0) > cutoff) return result;

  hyp1 = point[1] - base[1];
  if (fabs(hyp1) > cutoff) return result;

  hyp2 = point[2] - base[2];

  (*alongNormalSq) = (hyp2 * hyp2);
  result = (hyp0 * hyp0) + (hyp1 * hyp1);

  return result;
}

int LineClipPoint(float* base, float* ray, float* point, float* dist,
                  float cutoff, float cutoff2) {
  float hyp0, hyp1, hyp2;
  float opp0, opp1, opp2;
  float adj0, adj1, adj2;
  float ray0, ray1, ray2;
  float proj;
  double dcutoff = (double)cutoff;
  float opp_len_sq;

  /* this routine determines whether or not a vector starting at "base"
     heading in the direction "ray" intersects a sphere located at "point".

     It returns how far along the vector the intersection with the plane is or
     MAXFLOAT if there isn't a relevant intersection

     NOTE: this routine has been optimized for normals along Z
     Optimizes-out vectors that are more than "cutoff" from "point"
   */

  /* compute the hypo */

  hyp2 = point[2] - base[2];
  hyp1 = point[1] - base[1];
  hyp0 = point[0] - base[0];

  ray0 = ray[0];
  ray1 = ray[1];
  ray2 = ray[2];

  /* compute the adjacent edge (dot-projection) */

  proj = (ray0 * hyp0) + (ray1 * hyp1) + (ray2 * hyp2);

  adj0 = ray0 * proj;
  adj1 = ray1 * proj;
  opp0 = hyp0 - adj0;
  adj2 = ray2 * proj;
  if (fabs(opp0) > dcutoff) return 0;
  opp1 = hyp1 - adj1;
  opp2 = hyp2 - adj2;
  if (fabs(opp1) > dcutoff) return 0;
  if (fabs(opp2) > dcutoff) return 0;

  opp_len_sq = (opp0 * opp0) + (opp1 * opp1) + (opp2 * opp2);

  if (opp_len_sq <= cutoff2) {
    *dist = proj - (float)sqrt1f(cutoff2 - opp_len_sq);
    return 1;
  }

  return 0;
}

int LineClipEllipsoidPoint(float* base, float* ray, float* point, float* dist,
                           float cutoff, float cutoff2, float* scale, float* n1,
                           float* n2, float* n3) {
  float point_to_base[3];
  float scaled_base[3];
  float scaled_ray[3];
  float d1, d2, d3, s1, s2, s3;
  float comp1[3], comp2[3], comp3[3];

  subtract3f(base, point, point_to_base);

  /* project difference vector onto ellipsoid axes */

  d1 = dot_product3f(point_to_base, n1);
  d2 = dot_product3f(point_to_base, n2);
  d3 = dot_product3f(point_to_base, n3);

  s1 = d1 / scale[0];
  s2 = d2 / scale[1];
  s3 = d3 / scale[2];

  scale3f(n1, s1, comp1);
  scale3f(n2, s2, comp2);
  scale3f(n3, s3, comp3);

  copy3f(point, scaled_base);
  add3f(comp1, scaled_base, scaled_base);
  add3f(comp2, scaled_base, scaled_base);
  add3f(comp3, scaled_base, scaled_base);

  d1 = dot_product3f(ray, n1);
  d2 = dot_product3f(ray, n2);
  d3 = dot_product3f(ray, n3);

  s1 = d1 / scale[0];
  s2 = d2 / scale[1];
  s3 = d3 / scale[2];

  scale3f(n1, s1, comp1);
  scale3f(n2, s2, comp2);
  scale3f(n3, s3, comp3);

  copy3f(comp1, scaled_ray);
  add3f(comp2, scaled_ray, scaled_ray);
  add3f(comp3, scaled_ray, scaled_ray);

  d1 = length3f(scaled_ray);

  normalize3f(scaled_ray);

  /*  dump3f(ray,"ray");
     dump3f(scaled_ray,"scaled_ray");
   */

  {
    float hyp0, hyp1, hyp2;
    float opp0, opp1, opp2;
    float adj0, adj1, adj2;
    float ray0, ray1, ray2;
    float proj;
    double dcutoff = (double)cutoff;
    float opp_len_sq;

    /* this routine determines whether or not a vector starting at "base"
       heading in the direction "ray" intersects a sphere located at "point".

       It returns how far along the vector the intersection with the plane is or
       MAXFLOAT if there isn't a relevant intersection

     */

    /* compute the hypo */

    hyp2 = point[2] - scaled_base[2];
    hyp1 = point[1] - scaled_base[1];
    hyp0 = point[0] - scaled_base[0];

    ray0 = scaled_ray[0];
    ray1 = scaled_ray[1];
    ray2 = scaled_ray[2];

    /* compute the adjacent edge (dot-projection) */

    proj = (ray0 * hyp0) + (ray1 * hyp1) + (ray2 * hyp2);

    adj0 = ray0 * proj;
    adj1 = ray1 * proj;
    opp0 = hyp0 - adj0;
    adj2 = ray2 * proj;
    if (fabs(opp0) > dcutoff) return 0;
    opp1 = hyp1 - adj1;
    opp2 = hyp2 - adj2;
    if (fabs(opp1) > dcutoff) return 0;
    if (fabs(opp2) > dcutoff) return 0;

    opp_len_sq = (opp0 * opp0) + (opp1 * opp1) + (opp2 * opp2);

    if (opp_len_sq <= cutoff2) { /* line hits the virtual sphere */

      float scaled_dist = proj - (float)sqrt1f(cutoff2 - opp_len_sq);

      *(dist) = scaled_dist / d1;

      return 1;
    }

    return 0;
  }
}
