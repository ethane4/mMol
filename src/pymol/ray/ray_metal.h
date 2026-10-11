/*
 * RayMetal — PyMOL's `ray` tracer (RayTraceThread + BasisHitOrthoscopic) run
 * as a Metal compute kernel, one thread per pixel.
 *
 * Ray.cpp builds the primitives, bases and voxel map exactly as for the CPU
 * tracer, then (when enabled and the scene uses only what the kernel ports)
 * hands the flattened data here instead of running RayTraceThread. Everything
 * after it — edge oversampling, antialiasing, the ray_trace_mode outline —
 * stays on the CPU, so the output goes through the same post-processing.
 *
 * Ported: orthoscopic rays; spheres, cylinders, sausages, triangles; the full
 * lighting model (ambient / direct / reflect / legacy / power / specular,
 * spec_count, ray_scatter, two-sided lighting); fog; the multi-pass
 * transparency loop; ray_blend_colors; flat backgrounds.
 * Not ported (the CPU tracer runs instead): perspective, ray_shadows,
 * labels, ellipsoids, cones, textures (wobble), ramped colours, interior
 * colour, chromadepth, gradient / image backgrounds.
 */
#pragma once

#include <cstdint>
#include <string>

/* Plain-scalar layouts, mirrored field for field by the MSL source. */
struct RayMetalPrim {
  float c1[3], c2[3], c3[3], tr[3];
  float r1, l1, trans;
  int32_t vert, type, cap1, cap2, cull, no_lighting;
};

struct RayMetalParams {
  int32_t width, height, x_start, x_stop, y_start, y_stop;
  float invWdthRange, invHgtRange, vol0, vol2, border_offset;
  float front, back, fudge0, fudge1;

  int32_t n_vert, n_eElem;
  int32_t dim0, dim1, dim2, d1d2;
  int32_t iMin0, iMin1, iMin2, iMax0, iMax1, iMax2;
  float mapMin0, mapMin1, mapMin2, recipDiv;

  float ambient, lreflect, direct, direct_shade, legacy, legacy_1m, ray_scatter;
  float settingPower, settingReflectPower, settingSpecPower, settingSpecReflect,
      settingSpecDirect, settingSpecDirectPower;
  int32_t spec_count, n_basis, two_sided_lighting, has_bp2;
  float lightN[30], specN[30]; /* bases 2..n_basis-1, index (bc-2)*3 */
  float project_triangle;

  int32_t fogFlag, fogRangeFlag;
  float fog, fog_start, inv1minusFogStart, invFrontMinusBack;

  int32_t orig_opaque_back;
  uint32_t fore_mask;
  float bkrd0, bkrd1, bkrd2;

  int32_t trans_cont_flag, backface_cull, trans_mode, max_pass;
  float trans_cont, inv_trans_cont, trans_oblique, oblique_power;
  float ray_trans_spec, trans_spec_cut, trans_spec_scale;
  float trans_cutoff, persist_cutoff;

  int32_t blend_colors;
  float red_blend, green_blend, blue_blend;

  int32_t write_depth;
  uint32_t background;        /* what the host filled the image with */
  int32_t ox0, oy0, ox1, oy1; /* region the GPU writes: the traced rectangle,
                                 plus the outline's margin when it runs here */
  int32_t collect;            /* set by RayMetalTrace */
};

struct RayMetalScene {
  RayMetalParams p;
  const RayMetalPrim* prims;
  int n_prims;
  const float *vertex, *normal, *precomp, *radius, *radius2;
  size_t n_vertex_f, n_normal_f;
  const int32_t *vert2normal, *vert2prim;
  size_t n_vert2normal, n_vert2prim;
  const int32_t *ehead, *elist, *emask; /* PyMOL's express lists */
  /* lite map (no express lists): grid chains over the remapped vertices */
  int lite;
  const int32_t *mhead, *mlink, *literef;
};

/* RayRender's ray_trace_mode outline pass, run on the GPU after the trace */
struct RayMetalOutline {
  int32_t mode, fogFlag, fogRangeFlag, opaque;
  uint32_t trace_word;
  float fog, fog_start, inv1minusFogStart, invFrontMinusBack;
  float slope_f, depth_f, disco_f;
  float disco_raw; /* RayRender scales disco_f by the gain only after deriving
                      its 0.625/0.5/0.45 thresholds from it */
};

/* Trace `scene` into image (width*height RGBA words, prefilled with the
   background) and depth (may be null); with `outline`, run the outline pass
   too, and with `aa3` also the 3x antialias, writing only the final image to
   aa3. Returns false with `err` set if the GPU could not run it; the caller
   then uses the CPU tracer. */
bool RayMetalTrace(const RayMetalScene& scene, unsigned int* image,
                   float* depth, const RayMetalOutline* outline,
                   unsigned int* aa3, std::string& err, double* gpu_ms);

/* GPU tracer switch (on by default) and last-run report (read by _mmol_render;
   C linkage so the Swift bridge can declare them without the core's headers).
 */
extern "C" {
void RayMetalSetEnabled(int on);
int RayMetalEnabled();
void RayMetalSetLastReport(const char* msg);
const char* RayMetalLastReport();
void RayMetalSetStages(const char* msg);
const char* RayMetalStages();
}
