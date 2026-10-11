/*
 * PyMOL COLLADA export
 *
 * (c) 2014 Jared Sampson
 * (c) 2014 Schrodinger, Inc.
 */

#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"

#include "pymol/base/feedback.h"
#include "pymol/ray/ray.h"

/* Generates COLLADA output and appends it to `vla_ptr`. */
void RayRenderCOLLADA(CRay* I, int width, int height, char** vla_ptr,
                      float front, float back, float fov) {
  PRINTFB(I->G, FB_Ray, FB_Errors)
  " ColladaRender-Error: No libxml2 support, can't output COLLADA "
  "file.\n" ENDFB(I->G);
}
