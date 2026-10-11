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

// Helpers shared by the basis*.cc files; not part of the CBasis API.

#pragma once

#include "pymol/ray/basis.h"

static const float kR_SMALL4 = 0.0001F;

#define EPSILON 0.000001F

// Defined in basis_intersect.cc.
int ZLineToSphere(float* base, float* point, float* dir, float radius,
                  float maxial, float* sphere, float* asum, float* pre);
int LineToSphere(float* base, float* ray, float* point, float* dir,
                 float radius, float maxial, float* sphere, float* asum);
int FrontToInteriorSphere(float* front, float* point, float* dir, float radius,
                          float radius2, float maxial);
int ZLineToSphereCapped(float* base, float* point, float* dir, float radius,
                        float maxial, float* sphere, float* asum, cCylCap cap1,
                        cCylCap cap2, float* pre);
int LineToSphereCapped(float* base, float* ray, float* point, float* dir,
                       float radius, float maxial, float* sphere, float* asum,
                       cCylCap cap1, cCylCap cap2);
int ConeLineToSphereCapped(float* base, float* ray, float* point, float* dir,
                           float radius, float small_radius, float maxial,
                           float* sphere, float* asum, float* sph_rad,
                           float* sph_rad_sq, cCylCap cap1, cCylCap cap2);
int FrontToInteriorSphereCapped(float* front, float* point, float* dir,
                                float radius, float radius2, float maxial,
                                cCylCap cap1, cCylCap cap2);
float ZLineClipPoint(float* base, float* point, float* alongNormalSq,
                     float cutoff);
float ZLineClipPointNoZCheck(float* base, float* point, float* alongNormalSq,
                             float cutoff);
int LineClipPoint(float* base, float* ray, float* point, float* dist,
                  float cutoff, float cutoff2);
int LineClipEllipsoidPoint(float* base, float* ray, float* point, float* dist,
                           float cutoff, float cutoff2, float* scale, float* n1,
                           float* n2, float* n3);
