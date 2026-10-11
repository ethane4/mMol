
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

// Helpers shared by the ray*.cc files; not part of the ray API.

#pragma once

#include "pymol/ray/ray.h"

#define SettingGetfv SettingGetGlobal_3fv
#define RAY_SMALL 0.00001

/* BASES
   0 contains untransformed vertices (vector size = 3)
        1 contains transformed vertices (vector size = 3)
        2 contains transformed vertices for shadowing
*/

/* note: the following value must be at least one greater than the max
   number of lights */
#define MAX_BASIS 12

typedef float float3[3];

struct _CRayThreadInfo {
  CRay* ray;
  int width, height;
  unsigned int* image;
  float front, back;
  unsigned int fore_mask;
  float *bkrd_top, *bkrd_bottom;
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  float ambient;
  unsigned int background;
  int border;
  int phase, n_thread;
  int x_start, x_stop;
  int y_start, y_stop;
  unsigned int* edging;
  unsigned int edging_cutoff;
  int perspective;
  float fov, pos[3];
  float* depth;

  int bgWidth, bgHeight;
  void* bkrd_data; /* used for image-based background */
};

struct _CRayHashThreadInfo {
  CBasis* basis;
  int* vert2prim;
  CPrimitive* prim;
  int n_prim;
  float* clipBox;
  unsigned int* image;
  unsigned int background;
  size_t bytes;
  int perspective;
  float front;
  int phase;
  float size_hint;
  CRay* ray;
  float *bkrd_top, *bkrd_bottom;
  short bkrd_is_gradient; /* if not gradient, use bkrd_top as bkrd */
  int width, height;
  int opaque_back;
};

struct _CRayAntiThreadInfo {
  unsigned int* image;
  unsigned int* image_copy;
  unsigned int width, height;
  int mag;
  int phase, n_thread;
  CRay* ray;
  int no_busy; /* off Python's threads: don't touch the busy indicator */
};

// Defined in ray.cc.
void RayApplyMatrixInverse33(unsigned int n, float3* q, const float m[16],
                             float3* p);
void RayTransformNormals33(unsigned int n, float3* q, const float m[16],
                           float3* p);
void RayTransformInverseNormals33(unsigned int n, float3* q, const float m[16],
                                  float3* p);
void RayProjectTriangle(CRay* I, RayInfo* r, float* light, float* v0, float* n0,
                        float scale);

// Defined in ray_trace.cc.
void fill(unsigned int* buffer, unsigned int value, size_t cnt);
void fill_background_image(CRay* I, unsigned int* buffer, int width, int height,
                           unsigned int cnt);
void fill_gradient(CRay* I, int opaque_back, unsigned int* buffer,
                   float* bkrd_bottom, float* bkrd_top, int width, int height,
                   unsigned int cnt);
void RayHashSpawn(CRayHashThreadInfo* Thread, int n_thread, int n_total);
void RayTraceSpawn(CRayThreadInfo* Thread, int n_thread);
