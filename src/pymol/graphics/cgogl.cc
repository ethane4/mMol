#include <cassert>
#include "pymol/graphics/cgogl.h"

#include "pymol/graphics/cgo.h"
#include "pymol/graphics/cgo_renderer.h"
#include "pymol/color/color.h"
#include "pymol/objects/coord_set.h"
#include "pymol/base/feedback.h"
#include "pymol/graphics/gl/gl_vertex_buffer.h"
#include "pymol/graphics/graphics_util.h"
#include "pymol/base/pymol_globals.h"
#include "pymol/graphics/renderer.h"
#include "pymol/scene/scene.h"
#include "pymol/scene/scene_def.h"
#include "pymol/graphics/shader_mgr.h"
#include "pymol/geometry/sphere.h"
#include "pymol/platform/os_gl.h"
#include "pymol/platform/os_gl_cgo.h"

#include "pymol/graphics/vertex_format.h"

#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>

// Helper: extract attribute offsets from a VBO's CPU descriptor and draw
// via the Renderer. Returns true if handled (Metal path), false if not.
// Per-object interior-cap flag for a SURFACE mesh (closed solid). Cartoons
// share the indexed/non-indexed mesh paths but are not closed solids, so gate
// on the surface-shader flag captured at shader-enable. Resolve
// ray_interior_color for the interior cap and push it to the renderer. A real
// color (>=0) overrides the cap color; -1 (default) keeps the per- primitive
// default (atom/bond darkened for impostors, gray for surface).
static void metalApplyInteriorCapColor(CCGORenderer* I) {
  auto* G = I->G;
  CSetting* s1 = (I->rep && I->rep->cs) ? I->rep->cs->Setting.get() : nullptr;
  CSetting* s2 = (I->rep && I->rep->obj) ? I->rep->obj->Setting.get() : nullptr;
  int ci = SettingGet_color(G, s1, s2, cSetting_ray_interior_color);
  if (ci >= 0) {
    const float* c = ColorGet(G, ci);
    G->Renderer->setInteriorCapColor(c[0], c[1], c[2], true);
  } else {
    G->Renderer->setInteriorCapColor(0.32f, 0.32f, 0.36f, false);
  }
}

static int metalSurfaceInteriorCap(CCGORenderer* I) {
  if (!I || !I->metalIsSurfaceShader) return 0;
  CSetting* s1 = (I->rep && I->rep->cs) ? I->rep->cs->Setting.get() : nullptr;
  CSetting* s2 = (I->rep && I->rep->obj) ? I->rep->obj->Setting.get() : nullptr;
  if (!SettingGet_b(I->G, s1, s2, cSetting_metal_interior_cap)) return 0;
  metalApplyInteriorCapColor(I);
  return 1;
}

// Per-representation clipping: program the renderer's per-rep clip planes for
// the lit VBO draw that follows. Only the surface rep carries clip fractions
// (surface_clip_front/back); every other lit draw (cartoon/lines) must reset to
// "disabled" so a surface clip can't leak onto it (the renderer member
// persists).
//
// Reference is the surface object's center of mass (bounding-sphere center) and
// radius, NOT the camera slab: each fraction is normalized 0..1 across the
// molecule's depth, where 0 = no clip and 1 = clip that side all the way to the
// COM. Equal front/back give a slab symmetric about the COM, and the feel is
// the same on any structure (no camera-slab dead-zone). The COM depth is
// projected through the live modelview, so the slab stays centered as the
// camera moves.
static void metalApplyRepClip(CCGORenderer* I) {
  auto* G = I->G;
  if (!G->Renderer) return;
  // Gate on the ACTUAL rep type, not the metalIsSurfaceShader flag (that flag
  // is set on GL_SURFACE_SHADER enable; it must not leak the clip onto a later
  // rep).
  pymol::CObject* obj = I->rep ? I->rep->obj : nullptr;
  CoordSet* cs = I->rep ? I->rep->cs : nullptr;
  // Per-rep screen-space SSAO exemption (#79): cartoon/ribbon are excluded from
  // the SSAO crease/contour darkening (which paints spurious lines on their
  // silhouettes/self-folds) unless metal_ssao_cartoon opts them back in. Armed
  // on every path (before the surface early-returns) so a prior draw can't
  // leak in.
  {
    int repType = I->rep ? I->rep->type() : -1;
    bool aoExempt = (repType == cRepCartoon || repType == cRepRibbon) &&
                    !SettingGetGlobal_b(G, cSetting_metal_ssao_cartoon);
    G->Renderer->setRepScreenAO(aoExempt);
  }
  if (I->rep && I->rep->type() == cRepSurface && obj && cs &&
      cs->getNIndex() > 0) {
    CSetting* s1 = cs->Setting.get();
    CSetting* s2 = obj->Setting.get();
    // Outer-contour outline: arm capture for this surface draw when enabled.
    if (SettingGet_b(G, s1, s2, cSetting_surface_contour)) {
      float width = SettingGet_f(G, s1, s2, cSetting_surface_contour_width);
      int ci = SettingGet_color(G, s1, s2, cSetting_surface_contour_color);
      if (ci < 0) ci = SettingGet_color(G, s1, s2, cSetting_surface_color);
      float rgba[4] = {0.0f, 0.0f, 0.0f, 1.0f};  // default black, opaque
      if (ci >= 0) {
        const float* c = ColorGet(G, ci);
        rgba[0] = c[0];
        rgba[1] = c[1];
        rgba[2] = c[2];
      }
      if (!SettingGet_b(G, s1, s2, cSetting_surface_contour_opaque)) {
        float tr = SettingGet_f(G, s1, s2, cSetting_transparency);
        rgba[3] = 1.0f - (tr < 0.0f ? 0.0f : (tr > 1.0f ? 1.0f : tr));
      }
      G->Renderer->setRepContour(true, rgba, width);
    } else {
      G->Renderer->setRepContour(false, nullptr, 0.0f);
    }
    float fFrac = SettingGet_f(G, s1, s2, cSetting_surface_clip_front);
    float bFrac = SettingGet_f(G, s1, s2, cSetting_surface_clip_back);
    if (fFrac != 0.0f || bFrac != 0.0f) {
      // Project every atom to eye-space DEPTH and take the actual [min,max]
      // range. Using the true depth range (not a bounding-box estimate) matters
      // because `orient` puts the molecule's shortest axis toward the camera,
      // and a box projection over-estimates that depth — making the 0..1
      // fraction clip far too little until near 1. The modelview here is the
      // same matrix the shader uses, so the planes line up with the rasterized
      // eyeDist (verified).
      const glm::mat4& mv = SceneGetModelViewMatrix(G);
      int n = cs->getNIndex();
      float dMin = 1e30f, dMax = -1e30f;
      for (int i = 0; i < n; i++) {
        const float* c = cs->coordPtr(i);
        // eye.z = mv * c (column-major); depth from camera = -eye.z (positive).
        float depth =
            -(mv[0][2] * c[0] + mv[1][2] * c[1] + mv[2][2] * c[2] + mv[3][2]);
        dMin = std::min(dMin, depth);
        dMax = std::max(dMax, depth);
      }
      // Pad so the surface shell beyond the atoms stays enclosed (fraction 0 =
      // no clip on that side). Center is the COM along the view axis.
      float pad = SettingGet_f(G, s1, s2, cSetting_solvent_radius) + 1.0f;
      float center = 0.5f * (dMin + dMax);
      float halfD = 0.5f * (dMax - dMin) + pad;
      float ff = fFrac < 0.0f ? 0.0f : (fFrac > 1.0f ? 1.0f : fFrac);
      float bb = bFrac < 0.0f ? 0.0f : (bFrac > 1.0f ? 1.0f : bFrac);
      float front =
          center - halfD * (1.0f - ff);           // ff=0 -> near face (no clip)
      float back = center + halfD * (1.0f - bb);  // bb=0 -> far face (no clip)
      G->Renderer->setRepClip(front, back);
      return;
    }
    // Surface, no per-rep clip: keep the contour state set above, slab
    // disabled.
    G->Renderer->setRepClip(-1.0f, 1e6f);
    return;
  }
  G->Renderer->setRepContour(false, nullptr, 0.0f);  // non-surface: no contour
  G->Renderer->setRepClip(-1.0f, 1e6f);  // disabled: global slab only
}

static bool drawVBOViaMetal(CCGORenderer* I, VertexBufferGL* vbo,
                            pymol::PrimitiveType mode, int vertexCount) {
  auto* renderer = I ? I->G->Renderer : nullptr;
  if (!renderer || !vbo || !vbo->hasCPUData()) return false;

  const auto& desc = vbo->getDesc();
  size_t stride = vbo->cpuStride();
  int posOffset = -1, normalOffset = -1, colorOffset = -1;
  int colorType = 0;  // 0 = UByte4Norm, 1 = Float4

  for (const auto& d : desc.descs) {
    if (d.attr_name == "a_Vertex") {
      posOffset = static_cast<int>(d.offset);
    } else if (d.attr_name == "a_Normal") {
      // RendererMetal binds attribute 1 as Float3 unconditionally.
      // CGOPackedNormals (CGO.cpp) never packs normals on the Metal path, so
      // this holds by construction.
      assert(d.m_format == VertexFormat::Float3);
      normalOffset = static_cast<int>(d.offset);
    } else if (d.attr_name == "a_Color") {
      colorOffset = static_cast<int>(d.offset);
      if (d.m_format == VertexFormat::Float4) {
        colorType = 1;
      }
    }
  }

  metalApplyRepClip(I);
  renderer->drawVBO(mode, vertexCount, vbo->cpuData(), vbo->cpuDataSize(),
                    stride, posOffset, normalOffset, colorOffset, colorType,
                    metalSurfaceInteriorCap(I));
  return true;
}

static bool drawVBOIndexedViaMetal(CCGORenderer* I, VertexBufferGL* vbo,
                                   IndexBufferGL* ibo,
                                   pymol::PrimitiveType mode, int indexCount) {
  auto* renderer = I ? I->G->Renderer : nullptr;
  if (!renderer || !vbo || !vbo->hasCPUData() || !ibo || !ibo->hasCPUData())
    return false;

  const auto& desc = vbo->getDesc();
  size_t stride = vbo->cpuStride();
  int posOffset = -1, normalOffset = -1, colorOffset = -1;
  int colorType = 0;

  for (const auto& d : desc.descs) {
    if (d.attr_name == "a_Vertex") {
      posOffset = static_cast<int>(d.offset);
    } else if (d.attr_name == "a_Normal") {
      // RendererMetal binds attribute 1 as Float3 unconditionally.
      // CGOPackedNormals (CGO.cpp) never packs normals on the Metal path, so
      // this holds by construction.
      assert(d.m_format == VertexFormat::Float3);
      normalOffset = static_cast<int>(d.offset);
    } else if (d.attr_name == "a_Color") {
      colorOffset = static_cast<int>(d.offset);
      if (d.m_format == VertexFormat::Float4) {
        colorType = 1;
      }
    }
  }

  metalApplyRepClip(I);
  renderer->drawVBOIndexed(mode, indexCount, vbo->cpuData(), vbo->cpuDataSize(),
                           stride, posOffset, normalOffset, colorOffset,
                           colorType, ibo->cpuData(), ibo->cpuDataSize(),
                           metalSurfaceInteriorCap(I));
  return true;
}

// Metal path for analytic sphere impostors. The GL path below uses the GL
// "sphere" shader (CShaderPrg), absent in NO_OPENGL builds; instead route the
// impostor VBO to Renderer::drawSphereImpostors (per-pixel ray-sphere casting).
static void drawSphereImpostorsViaMetal(CCGORenderer* I,
                                        const cgo::draw::sphere_buffers* sp) {
  auto* G = I->G;
  auto* vbo = G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
  if (!vbo || !vbo->hasCPUData()) return;

  pymol::Renderer::SphereImpostorDrawCall call;
  call.data = vbo->cpuData();
  call.dataSize = vbo->cpuDataSize();
  call.stride = vbo->cpuStride();
  call.sphereCount = sp->num_spheres;
  for (const auto& d : vbo->getDesc().descs) {
    int off = static_cast<int>(d.offset);
    if (d.attr_name == "a_vertex_radius")
      call.posRadiusOff = off;
    else if (d.attr_name == "a_Color")
      call.colorOff = off;
    else if (d.attr_name == "a_rightUpFlags") {
      call.rightUpOff = off;
      call.rightUpIsFloat = (d.m_format == VertexFormat::Float) ? 1 : 0;
    }
  }
  // Radius is already baked into a_vertex_radius.w by the rep, so no extra
  // scale.
  call.sphereSizeScale = 1.0f;
  call.ortho = SettingGetGlobal_b(G, cSetting_ortho) ? 1 : 0;
  {
    // Per-object: fill the slab cross-section of clipped spheres with a solid
    // interior color instead of seeing through them (metal_interior_cap).
    CSetting* s1 = (I->rep && I->rep->cs) ? I->rep->cs->Setting.get() : nullptr;
    CSetting* s2 =
        (I->rep && I->rep->obj) ? I->rep->obj->Setting.get() : nullptr;
    call.interiorCap =
        SettingGet_b(G, s1, s2, cSetting_metal_interior_cap) ? 1 : 0;
    if (call.interiorCap) metalApplyInteriorCapColor(I);
  }

  G->Renderer->drawSphereImpostors(call);
}

// True if this VBO carries the cylinder-impostor attributes (sticks/bonds are
// emitted as a generic CGO_DRAW_CUSTOM op whose VBO has these names).
static bool vboLooksLikeCylinders(VertexBufferGL* vbo) {
  if (!vbo) return false;
  bool v1 = false, v2 = false, rad = false;
  for (const auto& d : vbo->getDesc().descs) {
    if (d.attr_name == "attr_vertex1")
      v1 = true;
    else if (d.attr_name == "attr_vertex2")
      v2 = true;
    else if (d.attr_name == "attr_radius")
      rad = true;
  }
  return v1 && v2 && rad;
}

// Metal path for analytic cylinder impostors. The GL "cylinder" shader is
// absent in NO_OPENGL builds, so route the impostor VBO + index buffer to
// Renderer::drawCylinderImpostors (per-pixel ray-cylinder). The cylinder rep
// emits its VBO via a generic CGO_DRAW_CUSTOM op, so this is called from
// CGO_gl_draw_custom once the VBO is recognized as cylinders.
static void drawCylinderImpostorsViaMetal(CCGORenderer* I, VertexBufferGL* vbo,
                                          IndexBufferGL* ibo, int indexCount) {
  auto* G = I->G;
  if (!vbo || !vbo->hasCPUData() || !ibo || !ibo->hasCPUData()) return;

  pymol::Renderer::CylinderImpostorDrawCall call;
  call.vdata = vbo->cpuData();
  call.vdataSize = vbo->cpuDataSize();
  call.stride = vbo->cpuStride();
  call.idata = ibo->cpuData();
  call.idataSize = ibo->cpuDataSize();
  call.indexCount =
      (indexCount > 0)
          ? indexCount
          : static_cast<int>(ibo->cpuDataSize() / sizeof(VertexIndex_t));
  call.cylinderCount = call.indexCount / 36;  // 36 indices per cylinder box
  for (const auto& d : vbo->getDesc().descs) {
    int off = static_cast<int>(d.offset);
    if (d.attr_name == "attr_vertex1")
      call.v1Off = off;
    else if (d.attr_name == "attr_vertex2")
      call.v2Off = off;
    else if (d.attr_name == "a_Color") {
      call.colorOff = off;
      call.colorIsFloat = (d.m_format == VertexFormat::Float4) ? 1 : 0;
    } else if (d.attr_name == "a_Color2")
      call.color2Off = off;
    else if (d.attr_name == "attr_radius")
      call.radiusOff = off;
    else if (d.attr_name == "attr_flags") {
      call.flagsOff = off;
      call.flagsIsFloat = (d.m_format == VertexFormat::Float) ? 1 : 0;
    } else if (d.attr_name == "a_cap")
      call.capOff = off;
  }
  // Sticks bake their radius into attr_radius and leave uniRadius 0; dashes
  // (distances/dihedrals/angles) supply the real radius via the
  // CYLINDER_WIDTH_FOR_DISTANCES special op, captured into metalCylUniRadius.
  call.uniRadius = I->metalCylUniRadius;
  // a_cap: per-vertex when the CGO baked the cap/interp bits per cylinder (the
  // stick rep always does — RepCylBond varies caps per bond). Only when the CGO
  // emitted ONE constant for every cylinder does it arrive via
  // CGO_VERTEX_ATTRIBUTE_1F, captured into metalCylCapConst. Before GitHub
  // issue #441, capOff didn't exist and every cylinder used the captured
  // constant, so one cgo.CYLINDER (the Move gizmo) permanently re-capped every
  // stick.
  call.capConst = I->metalCylCapConst;
  call.noFlatCaps = 1;  // round caps (matches Get_CylinderShader GL default)
  call.ortho = SettingGetGlobal_b(G, cSetting_ortho) ? 1 : 0;
  {
    // Per-object slab interior cap (metal_interior_cap) — same as spheres.
    CSetting* s1 = (I->rep && I->rep->cs) ? I->rep->cs->Setting.get() : nullptr;
    CSetting* s2 =
        (I->rep && I->rep->obj) ? I->rep->obj->Setting.get() : nullptr;
    call.interiorCap =
        SettingGet_b(G, s1, s2, cSetting_metal_interior_cap) ? 1 : 0;
    if (call.interiorCap) metalApplyInteriorCapColor(I);
  }

  G->Renderer->drawCylinderImpostors(call);
}

// Label/text Metal helpers (defined later, after Texture.h is included).
static bool vboLooksLikeLabels(VertexBufferGL* vbo);
static void drawLabelsViaMetal(CCGORenderer* I, size_t vboid, int vertexCount);
static void drawSphereImpostorsViaMetal(CCGORenderer* I,
                                        const cgo::draw::sphere_buffers* sp);
static bool vboLooksLikeCylinders(VertexBufferGL* vbo);
static void drawCylinderImpostorsViaMetal(CCGORenderer* I, VertexBufferGL* vbo,
                                          IndexBufferGL* ibo, int indexCount);

/* ======== GL Rendering ======== */

static pymol::PrimitiveType glModeToPrimitive(int mode) {
  switch (mode) {
    case GL_POINTS:
      return pymol::PrimitiveType::Points;
    case GL_LINES:
      return pymol::PrimitiveType::Lines;
    case GL_LINE_STRIP:
      return pymol::PrimitiveType::LineStrip;
    case GL_TRIANGLES:
      return pymol::PrimitiveType::Triangles;
    case GL_TRIANGLE_STRIP:
      return pymol::PrimitiveType::TriangleStrip;
    case GL_TRIANGLE_FAN:
      return pymol::PrimitiveType::TriangleFan;
    default:
      return pymol::PrimitiveType::Triangles;
  }
}

static int CGO_gl_begin_WARNING_CALLED = false,
           CGO_gl_end_WARNING_CALLED = false,
           CGO_gl_vertex_WARNING_CALLED = false;

const float g_ones4f[4] = {1.f, 1.f, 1.f, 1.f};

static int CGOConvertDebugMode(int debug, int modeArg) {
  int mode = modeArg;
  if (debug == 1) {
    switch (mode) {
      case GL_TRIANGLES:
        mode = GL_LINES;
        break;
      case GL_TRIANGLE_STRIP:
        mode = GL_LINE_STRIP;
        break;
      case GL_TRIANGLE_FAN:
        mode = GL_LINES;
        break;
    }
  } else {
    mode = GL_POINTS;
  }
  return mode;
}

static void CGO_gl_begin(CCGORenderer* I, CGO_op_data pc) {
  if (!CGO_gl_begin_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings)
    " CGO_gl_begin() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_begin_WARNING_CALLED = true;
  }
}
static void CGO_gl_end(CCGORenderer* I, CGO_op_data) {
  if (!CGO_gl_end_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings)
    " CGO_gl_end() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_end_WARNING_CALLED = true;
  }
}
static void CGO_gl_vertex(CCGORenderer* I, CGO_op_data v) {
  if (!CGO_gl_vertex_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings)
    " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_vertex_WARNING_CALLED = true;
  }
}

static void CGO_gl_vertex_cross(CCGORenderer* I, CGO_op_data v) {
  if (!CGO_gl_vertex_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings)
    " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_vertex_WARNING_CALLED = true;
  }
}

static void CGO_gl_line(CCGORenderer* I, CGO_op_data v) {}

static void CGO_gl_splitline(CCGORenderer* I, CGO_op_data v) {}

static void CGO_gl_normal(CCGORenderer* I, CGO_op_data varg) {
  const float* v = *varg;
  if (I->use_shader) {
    glVertexAttrib3fv(VERTEX_NORMAL, v);
  } else {
  }
}

static void CGO_gl_draw_buffers_indexed(CCGORenderer* I, CGO_op_data pc) {
  auto sp = reinterpret_cast<const cgo::draw::buffers_indexed*>(*pc);

  // Metal path: draw indexed using retained CPU data through the Renderer
  if (I->G->Renderer) {
    auto* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
    auto* ibo = I->G->ShaderMgr->getGPUBuffer<IndexBufferGL>(sp->iboid);
    drawVBOIndexedViaMetal(I, vbo, ibo, glModeToPrimitive(sp->mode),
                           sp->nindices);
  }
}

static void CGO_gl_draw_buffers_not_indexed(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::buffers_not_indexed* sp =
      reinterpret_cast<decltype(sp)>(*pc);

  // Metal path: draw using retained CPU data through the Renderer
  if (I->G->Renderer) {
    auto* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
    drawVBOViaMetal(I, vbo, glModeToPrimitive(sp->mode), sp->nverts);
  }
}

static void CGO_gl_draw_custom(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::custom* sp = reinterpret_cast<decltype(sp)>(*pc);
  if (I->G->Renderer) {
    auto* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
    // Text/label CGOs are emitted as custom draws but need the screen-aligned
    // glyph-atlas shader rather than the geometry VBO shader.
    if (vboLooksLikeLabels(vbo)) {
      drawLabelsViaMetal(I, sp->vboid, sp->nverts);
      return;
    }
    if (vboLooksLikeCylinders(vbo)) {
      auto* ibo = sp->iboid
                      ? I->G->ShaderMgr->getGPUBuffer<IndexBufferGL>(sp->iboid)
                      : nullptr;
      drawCylinderImpostorsViaMetal(I, vbo, ibo, sp->nindices);
      return;
    }
    if (sp->iboid) {
      auto* ibo = I->G->ShaderMgr->getGPUBuffer<IndexBufferGL>(sp->iboid);
      drawVBOIndexedViaMetal(I, vbo, ibo, glModeToPrimitive(sp->mode),
                             sp->nindices);
    } else {
      drawVBOViaMetal(I, vbo, glModeToPrimitive(sp->mode), sp->nverts);
    }
  }
}

static void CGO_gl_draw_sphere_buffers(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::sphere_buffers* sp = reinterpret_cast<decltype(sp)>(*pc);
  if (I->G->Renderer) {
    drawSphereImpostorsViaMetal(I, sp);
  }
}

static void CGO_gl_draw_bezier_buffers(CCGORenderer* I, CGO_op_data cgo_data) {
  const auto bezier =
      reinterpret_cast<const cgo::draw::bezier_buffers*>(*cgo_data);
  if (I->G->Renderer) {
    // Metal: GPU-tessellate the cubic Bezier control points into a smooth tube
    // (the GL "bezier" tessellation shader is absent in NO_OPENGL builds).
    auto* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBufferGL>(bezier->vboid);
    if (vbo && vbo->hasCPUData()) {
      float radius =
          SettingGet_f(I->G, I->set1, I->set2, cSetting_cartoon_tube_radius);
      if (radius <= 0.f) radius = 0.5f;
      const float* col = I->color ? I->color : g_ones4f;
      I->G->Renderer->drawBezierTubes(vbo->cpuData(), vbo->cpuDataSize(),
                                      radius, col[0], col[1], col[2]);
    }
  }
}

static void CGO_gl_draw_cylinder_buffers(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::cylinder_buffers* sp = reinterpret_cast<decltype(sp)>(*pc);
  if (I->G->Renderer) {
    // Legacy cylinder_buffers op. The current stick rep instead emits a generic
    // CGO_DRAW_CUSTOM op (handled in CGO_gl_draw_custom), but route this here
    // too for completeness.
    auto* vbo = I->G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
    auto* ibo = I->G->ShaderMgr->getGPUBuffer<IndexBufferGL>(sp->iboid);
    drawCylinderImpostorsViaMetal(I, vbo, ibo, 0);
  }
}
#include "pymol/graphics/texture.h"

// True if this VBO is a text/label VBO (built by CGOConvertToLabelShader),
// identified by its screen-aligned glyph attributes.
static bool vboLooksLikeLabels(VertexBufferGL* vbo) {
  if (!vbo) return false;
  bool hasWorld = false, hasTex = false;
  for (const auto& d : vbo->getDesc().descs) {
    if (d.attr_name == "attr_worldpos")
      hasWorld = true;
    else if (d.attr_name == "attr_texcoords")
      hasTex = true;
  }
  return hasWorld && hasTex;
}

// Metal path for text/label quads. The GL path uses the "label" GL shader
// (CShaderPrg), which doesn't exist in NO_OPENGL builds; instead route the
// interleaved label VBO + glyph atlas to Renderer::drawLabels.
static void drawLabelsViaMetal(CCGORenderer* I, size_t vboid, int vertexCount) {
  auto* G = I->G;
  if (!G->Renderer) return;

  auto* vbo = G->ShaderMgr->getGPUBuffer<VertexBufferGL>(vboid);
  if (!vbo || !vbo->hasCPUData()) return;

  pymol::Renderer::LabelDrawCall call;
  call.data = vbo->cpuData();
  call.dataSize = vbo->cpuDataSize();
  call.stride = vbo->cpuStride();
  call.vertexCount = vertexCount;

  for (const auto& d : vbo->getDesc().descs) {
    int off = static_cast<int>(d.offset);
    if (d.attr_name == "attr_worldpos")
      call.worldPosOff = off;
    else if (d.attr_name == "attr_targetpos")
      call.targetPosOff = off;
    else if (d.attr_name == "attr_screenoffset")
      call.screenOff = off;
    else if (d.attr_name == "attr_texcoords")
      call.texOff = off;
    else if (d.attr_name == "attr_screenworldoffset")
      call.screenWorldOff = off;
    else if (d.attr_name == "attr_relative_mode")
      call.relModeOff = off;
  }

  auto extent = SceneGetExtentStereo(G);
  call.screenW = extent.width;
  call.screenH = extent.height;
  call.screenOriginVertexScale = SceneGetScreenVertexScale(G, nullptr) / 2.f;
  float front = SceneGetCurrentFrontSafe(G);
  float back = SceneGetCurrentBackSafe(G);
  call.front = front;
  call.clipRange = back - front;

  // Mirror SET_LABEL_SCALE_UNIFORMS: scale-by-vertex when label_size < 0.
  CSetting *s1 = nullptr, *s2 = nullptr;
  if (I->rep && I->rep->cs) s1 = I->rep->cs->Setting.get();
  if (I->rep && I->rep->obj) s2 = I->rep->obj->Setting.get();
  float label_size = SettingGet_f(G, s1, s2, cSetting_label_size);
  call.scaleByVertexScale = label_size < 0.f ? 1.f : 0.f;
  float tfs = I->info ? static_cast<float>(I->info->texture_font_size) : 0.f;
  call.labelTextureSize = label_size < 0.f ? (-2.f * tfs / label_size) : 1.f;

  const unsigned char* px = nullptr;
  int aw = 0, ah = 0;
  unsigned long long gen = 0;
  if (!TextureGetTextTextureMetalData(G, &px, &aw, &ah, &gen)) return;
  call.atlasPixels = px;
  call.atlasW = aw;
  call.atlasH = ah;
  call.atlasGen = gen;

  G->Renderer->drawLabels(call);
}

static void CGO_gl_draw_labels(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::labels* sp = reinterpret_cast<decltype(sp)>(*pc);

  if (I->G->Renderer) {
    drawLabelsViaMetal(I, sp->vboid, sp->ntextures * 6);
  }
}

// Metal path for label connectors (background box, background outline and the
// connector line). The GL path below relies on connector.gs to amplify one
// vertex per connector into all of that geometry; Metal has no geometry
// shaders, so hand the un-expanded per-connector records to
// Renderer::drawConnectors, which does the amplification in its vertex shader.
static void drawConnectorsViaMetal(CCGORenderer* I,
                                   const cgo::draw::connectors* sp) {
  auto* G = I->G;
  auto* vbo = G->ShaderMgr->getGPUBuffer<VertexBufferGL>(sp->vboid);
  if (!vbo || !vbo->hasCPUData()) return;

  pymol::Renderer::ConnectorDrawCall call;
  call.data = vbo->cpuData();
  call.dataSize = vbo->cpuDataSize();
  call.stride = vbo->cpuStride();
  call.connectorCount = sp->nconnectors;

  for (const auto& d : vbo->getDesc().descs) {
    int off = static_cast<int>(d.offset);
    if (d.attr_name == "a_target_pt3d")
      call.targetOff = off;
    else if (d.attr_name == "a_center_pt3d")
      call.centerOff = off;
    else if (d.attr_name == "a_indentFactor")
      call.indentOff = off;
    else if (d.attr_name == "a_screenWorldOffset")
      call.screenWorldOff = off;
    else if (d.attr_name == "a_textSize")
      call.textSizeOff = off;
    else if (d.attr_name == "a_Color")
      call.colorOff = off;
    else if (d.attr_name == "a_relative_mode")
      call.relModeOff = off;
    else if (d.attr_name == "a_draw_flags")
      call.drawFlagsOff = off;
    else if (d.attr_name == "a_bkgrd_color")
      call.bgColorOff = off;
    else if (d.attr_name == "a_rel_ext_length")
      call.relExtLenOff = off;
    else if (d.attr_name == "a_con_width")
      call.conWidthOff = off;
  }

  auto extent = SceneGetExtentStereo(G);
  call.screenW = extent.width;
  call.screenH = extent.height;
  call.screenOriginVertexScale = SceneGetScreenVertexScale(G, nullptr) / 2.f;
  float front = SceneGetCurrentFrontSafe(G);
  float back = SceneGetCurrentBackSafe(G);
  call.front = front;
  call.clipRange = back - front;
  call.lineSmooth = SettingGetGlobal_b(G, cSetting_line_smooth) ? 1 : 0;

  // Mirrors the "textureToLabelSize" uniform of the GL connector shader.
  CSetting *s1 = nullptr, *s2 = nullptr;
  if (I->rep && I->rep->cs) s1 = I->rep->cs->Setting.get();
  if (I->rep && I->rep->obj) s2 = I->rep->obj->Setting.get();
  float label_size = SettingGet_f(G, s1, s2, cSetting_label_size);
  if (label_size < 0.f && I->info) {
    float v_scale = SceneGetScreenVertexScale(G, nullptr);
    call.textureToLabelSize =
        v_scale * static_cast<float>(I->info->texture_font_size) / label_size;
  } else {
    call.textureToLabelSize = 1.f;
  }

  G->Renderer->drawConnectors(call);
}

static void CGO_gl_draw_connectors(CCGORenderer* I, CGO_op_data pc) {
  const cgo::draw::connectors* sp = reinterpret_cast<decltype(sp)>(*pc);

  if (I->G->Renderer) {
    drawConnectorsViaMetal(I, sp);
  }
}

static void CGO_gl_linewidth(CCGORenderer* I, CGO_op_data pc) {
  glLineWidth(**pc);
}

/* CGO_gl_special - this is the implementation function for
   CGOSpecial/CGO_SPECIAL. Only the ops that feed the Metal renderer do
   anything; the rest set GL shader uniforms and are no-ops.
 */
static void CGO_gl_special(CCGORenderer* I, CGO_op_data pc) {
  int mode = CGO_get_int(*pc);

  CSetting *csSetting = nullptr, *objSetting = nullptr;
  if (I->rep && I->rep->cs) {
    csSetting = I->rep->cs->Setting.get();
  }
  if (I->rep && I->rep->obj) {
    objSetting = I->rep->obj->Setting.get();
  }
  switch (mode) {
    case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON:
      // Per-rep ribbon_width, as the raw width (issue #71; see
      // LINEWIDTH_FOR_LINES).
      if (I->G->Renderer) {
        I->G->Renderer->lineWidth(
            SettingGet_f(I->G, csSetting, objSetting, cSetting_ribbon_width));
      }
      break;
    case LINEWIDTH_DYNAMIC_MESH: {
      float line_width;
      if (I->rep) {
        line_width =
            SettingGet_f(I->G, I->rep->cs->Setting.get(),
                         I->rep->obj->Setting.get(), cSetting_mesh_width);
      } else {
        line_width = SettingGet_f(I->G, nullptr, nullptr, cSetting_mesh_width);
      }
      // Pre-dynamic-scale value matches the renderer's _lineWidth units.
      if (I->G->Renderer) {
        I->G->Renderer->lineWidth(line_width);
      }
    } break;
    case CYLINDER_WIDTH_FOR_DISTANCES:
      if (I->G->Renderer) {
        // RepDistDash (and RepDihedral/RepAngle) bake a 1.0 placeholder into
        // attr_radius and rely on this op to supply the real dash radius, so
        // capture it for the subsequent cylinder-impostor draw.
        float radius =
            SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_radius);
        if (radius == 0.0F) {
          float dash_size =
              SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_width);
          float line_width = SceneGetDynamicLineWidth(I->info, dash_size);
          I->metalCylUniRadius = SceneGetLineWidthForCylindersStatic(
              I->G, I->info, line_width, dash_size);
        } else {
          I->metalCylUniRadius = radius;
        }
      }
      break;
    case LINEWIDTH_FOR_LINES:
      // Per-rep/object width, as the RAW width: drawLinesAA reads _lineWidth
      // in these raw units (issue #71).
      if (I->G->Renderer) {
        I->G->Renderer->lineWidth(
            SettingGet_f(I->G, csSetting, objSetting, cSetting_line_width));
      }
      break;
    case LINEWIDTH_DYNAMIC_WITH_SCALE:
    case POINTSIZE_DYNAMIC_DOT_WIDTH:
    case LINEWIDTH_DYNAMIC_WITH_SCALE_DASH:
    case CYLINDERWIDTH_DYNAMIC_MESH:
    case LINEWIDTH_WITH_SCALE:
    case DOTSIZE_WITH_SPHERESCALE:
    case MESH_WIDTH_FOR_SURFACES:
    case CYLINDER_WIDTH_FOR_RIBBONS:
    case DOT_WIDTH_FOR_DOTS:
    case DOT_WIDTH_FOR_DOT_SPHERES:
    case CYLINDER_WIDTH_FOR_NONBONDED:
    case CYLINDER_WIDTH_FOR_REPWIRE:
    case CYLINDER_WIDTH_FOR_REPWIRE_VARWIDTH:
    case ENABLE_BACK_FACES_IF_NOT_TWO_SIDED:
    case DISABLE_BACK_FACES_IF_NOT_TWO_SIDED:
    case SET_SURFACE_UNIFORMS:
    case SET_ALIGNMENT_UNIFORMS_ATTRIBS:
    case SET_LABEL_SCALE_UNIFORMS:
      break;
    default:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGO_gl_special(): bad mode=%d\n", mode ENDFB(I->G);
  }
}

/* CGO_gl_special_with_arg - this is the implementation function for
   CGOSpecialWithArg/CGO_SPECIAL_WITH_ARG.  Each op has its own implementation.
 */
static void CGO_gl_special_with_arg(CCGORenderer* I, CGO_op_data pc) {}

static void CGO_gl_dotwidth(CCGORenderer* I, CGO_op_data pc) {}

/* CGO_gl_enable/CGO_gl_disable - only the Metal per-rep state the shader
   enable/disable ops carry is kept; the GL state and shader programs are
   gone.
 */
static void CGO_gl_enable(CCGORenderer* I, CGO_op_data pc) {
  if (!I->use_shader || !I->G->Renderer) return;
  switch (CGO_get_int(*pc)) {
    case GL_DEFAULT_SHADER:
    case GL_LINE_SHADER:
      I->metalIsSurfaceShader = false;
      break;
    case GL_SURFACE_SHADER:
      // Metal interior-cap: mark that following indexed-mesh draws are a
      // (closed) surface, so they may be stencil-capped when clipped.
      I->metalIsSurfaceShader = true;
      break;
    case GL_CYLINDER_SHADER:
      // uni_radius and the a_cap constant are per-rep state. Reset both at
      // shader-enable so sticks (which emit neither
      // CYLINDER_WIDTH_FOR_DISTANCES nor a constant a_cap) start from the
      // defaults; a following special op / attribute op sets them before that
      // rep's draw. Without the cap reset, one CGO's constant leaked into
      // every later rep for the life of the process (GitHub issue #441).
      I->metalCylUniRadius = 0.0f;
      I->metalCylCapConst = 15.0f;  // cCylShaderBothCapsRound
      break;
  }
}

static void CGO_gl_disable(CCGORenderer* I, CGO_op_data pc) {
  // Reset the surface-shader flag on disable so it can't stay sticky and leak
  // surface-only behavior (per-rep clip, metal_interior_cap) onto a rep drawn
  // later in the frame.
  if (I->use_shader && I->G->Renderer &&
      CGO_get_int(*pc) == GL_SURFACE_SHADER) {
    I->metalIsSurfaceShader = false;
  }
}

static void CGO_gl_alpha(CCGORenderer* I, CGO_op_data pc) { I->alpha = **pc; }

static void CGO_gl_reset_normal(CCGORenderer* I, CGO_op_data pc) {
  SceneResetNormalUseShader(I->G, CGO_get_int(*pc), I->use_shader);
}

static void CGO_gl_null(CCGORenderer* I, CGO_op_data pc) {}

static void CGO_gl_error(CCGORenderer* I, CGO_op_data pc) {
  PRINTFB(I->G, FB_CGO, FB_Warnings)
  " CGO_gl_error() is not suppose to be called op=%d\n",
      CGO_get_int((*pc) - 1) ENDFB(I->G);
}

static void CGO_gl_vertex_attribute_1f(CCGORenderer* I, CGO_op_data varg) {
  auto vertex_attr =
      reinterpret_cast<const cgo::draw::vertex_attribute_1f*>(*varg);
  const char* name =
      I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx);
  if (I->G->Renderer) {
    // Metal path: no GL shader to set a generic attribute on. The cylinder
    // impostor reads `a_cap` as a uniform; capture the constant here so the
    // subsequent draw_cylinder_buffers op can pass it along.
    if (name && strcmp(name, "a_cap") == 0)
      I->metalCylCapConst = vertex_attr->value;
  }
}

/* dispatch table for OpenGL */

CGO_op_fn CGO_gl[] = {
    CGO_gl_null,   /* 0x00 */
    CGO_gl_null,   /* 0x01 */
    CGO_gl_begin,  /* 0x02 */
    CGO_gl_end,    /* 0x03 */
    CGO_gl_vertex, /* 0x04 */
    CGO_gl_normal, /* 0x05 */
    CGO_gl_null,   /* 0x06 */
    CGO_gl_null,   /* 0x07 */
    CGO_gl_null,   /* 0x08 */
    CGO_gl_null,   /* 0x09 */

    CGO_gl_linewidth, /* 0x0A */
    CGO_gl_null,      /* 0x0B */
    CGO_gl_enable,    /* 0x0C */
    CGO_gl_disable,   /* 0x0D */
    CGO_gl_null,      /* 0x0E */
    CGO_gl_null,      /* 0x0F */

    CGO_gl_dotwidth, /* 0X10 */
    CGO_gl_null,     /* 0x11 */
    CGO_gl_null,     /* 0x12 */
    CGO_gl_null,     /* 0X13 */

    CGO_gl_null, /* 0X14 */
    CGO_gl_null, /* 0x15 */
    CGO_gl_null, /* 0x16 */
    CGO_gl_null, /* 0X17 */

    CGO_gl_null,                     /* 0X18 */
    CGO_gl_alpha,                    /* 0x19 */
    CGO_gl_null,                     /* 0x1A */
    CGO_gl_null,                     /* 0X1B */
    CGO_gl_null,                     /* 0x1C DrawArrays() */
    CGO_gl_null,                     /* 0x1D */
    CGO_gl_reset_normal,             /* 0x1E */
    CGO_gl_null,                     /* pick color  0X1F */
    CGO_gl_null,                     /* 0x20 draw buffers REMOVED */
    CGO_gl_draw_buffers_indexed,     /* 0x21 draw buffers indexed */
    CGO_gl_null,                     /* 0x22 bounding box */
    CGO_gl_draw_buffers_not_indexed, /* 0x23 draw buffers not indexed */
    CGO_gl_special,                  /* 0x24 special */
    CGO_gl_draw_cylinder_buffers,    /* 0x25 draw GLSL cylinders */
    CGO_gl_null,                     /* 0x26 shader cylinder */
    CGO_gl_null,                     /* 0x27 shader cylinder with 2nd color */
    CGO_gl_draw_sphere_buffers,      /* 0x28 draw sphere buffers */
    CGO_gl_null,  /* 0x29 accessibility used for ambient occlusion */
    CGO_gl_error, /* 0x2A draw texture */
    CGO_gl_null,  /* 0x2B draw textures */
    CGO_gl_null,  /* 0x2C draw screen textures and
                                                  polygons */
    CGO_gl_error,
    CGO_gl_error,
    CGO_gl_draw_labels,
    CGO_gl_error,
    CGO_gl_draw_connectors,
    CGO_gl_null,
    CGO_gl_null,
    CGO_gl_special_with_arg,
    CGO_gl_line,
    CGO_gl_splitline,
    CGO_gl_draw_custom,
    CGO_gl_null,
    CGO_gl_null,
    CGO_gl_vertex_attribute_1f,
    CGO_gl_null,
    CGO_gl_null,
    CGO_gl_vertex,
    CGO_gl_null,          // interpolated
    CGO_gl_vertex_cross,  // CGO_VERTEX_CROSS
    CGO_gl_null,
    CGO_gl_null,  // custom cylinder alpha
    CGO_gl_null,  // bezier
    CGO_gl_draw_bezier_buffers,
    CGO_gl_error};

void CGORenderGL(CGO* I, const float* color, CSetting* set1, CSetting* set2,
                 RenderInfo* info, Rep* rep)
/* this should be as fast as you can make it...

 * the ASM loop is about 2X long as raw looped GL calls,

 * but hopefully superscaler processors won't care */
{
  PyMOLGlobals* G = I->G;

  const float zee[] = {0.f, 0.f, 1.f};
  const float color_tmp[] = {1.f, 1.f, 1.f};

  if (I->render_alpha) {
    // for now, the render_alpha_only flag calls CGOSetZVector/CGORenderGLAlpha
    auto modMatrix = SceneGetModelViewMatrixPtr(G);
    CGOSetZVector(I, modMatrix[2], modMatrix[6], modMatrix[10]);
    CGORenderAlpha(I, info, 1);
    if (I->render_alpha == 1)  // right now, render_alpha 1: renders alpha only,
                               // 2: renders both alpha and rest
      return;
  }

  if (!G->ValidContext) {
    return;
  }

  if (!I->c) {
    return;
  }

  {
    CCGORenderer* R = G->CGORenderer;
    R->info = info;
    R->use_shader = I->use_shader;
    R->debug = I->debug;
    R->sphere_quality = I->sphere_quality;
    R->rep = rep;
    R->color = color;
    R->alpha = 1.0F - SettingGet_f(G, set1, set2, cSetting_cgo_transparency);
    R->set1 = set1;
    R->set2 = set2;
    // normals should be initialized to the view vector
    // (changed BB 9/14 from SceneResetNormalUseShader(), to
    // CScene->LinesNormal, which was arbitrary, I believe)
    SceneResetNormalToViewVector(I->G, I->use_shader);

    if (!color) {
      color = color_tmp;
    }

    if (!(info && info->alpha_cgo)) {
      // Regular CGO dispatch table rendering
      for (auto it = I->begin(); !it.is_stop(); ++it) {
        const auto op = it.op_code();
        assert(op < CGO_sz_size());
        CGO_OP_DATA_CONST float* const pc = it.data();

        CGO_gl[op](R, &pc);
      }
      return;
    }

    /* we're sorting transparent triangles globally */
    {
      {
        int mode = -1;
        int vc = 0;
        // triangle normals
        const float *n0 = nullptr, *n1 = nullptr, *n2 = nullptr;
        // triangle vertices
        const float *v0 = nullptr, *v1 = nullptr, *v2 = nullptr;
        // triangle colors
        const float *c0 = color, *c1 = nullptr, *c2 = nullptr;

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          const auto op = it.op_code();
          assert(op < CGO_sz_size());
          CGO_OP_DATA_CONST float* const pc = it.data();

          if ((R->alpha != 1.f)) {
            switch (op) { /* transparency */
              case CGO_BEGIN:
                mode = CGO_get_int(pc);
                CGO_gl_begin(R, &pc);
                vc = 0;
                n0 = zee;
                break;
              case CGO_END:
                CGO_gl_end(R, &pc);
                mode = -1;
                break;
              case CGO_NORMAL:
                switch (mode) {
                  case GL_TRIANGLES:
                  case GL_TRIANGLE_STRIP:
                  case GL_TRIANGLE_FAN:
                    n0 = pc;
                    break;
                  default:
                    CGO_gl_normal(R, &pc);
                }
                break;
              case CGO_COLOR:
                c0 = pc;
                break;
              case CGO_TRIANGLE:
                CGOAlphaTriangle(info->alpha_cgo, pc, pc + 3, pc + 6, pc + 9,
                                 pc + 12, pc + 15, pc + 18, pc + 21, pc + 24,
                                 R->alpha, R->alpha, R->alpha, false);
                break;
              case CGO_DRAW_ARRAYS: {
                const cgo::draw::arrays* sp =
                    reinterpret_cast<decltype(sp)>(pc);
                int mode = sp->mode, arrays = sp->arraybits,
                    nverts = sp->nverts;
                float *vertexVals = 0, *nxtVals = 0, *colorVals = 0,
                      *normalVals;
                float *vertexVals_tmp = 0, *colorVals_tmp = 0,
                      *normalVals_tmp = 0;
                int step;
                short nxtn = 3;
                nxtVals = vertexVals = vertexVals_tmp = sp->floatdata;
                if (arrays & CGO_NORMAL_ARRAY) {
                  nxtVals = normalVals = normalVals_tmp =
                      vertexVals + (nxtn * nverts);
                }
                if (arrays & CGO_COLOR_ARRAY) {
                  nxtVals = colorVals = colorVals_tmp =
                      nxtVals + (nxtn * nverts);
                  nxtn = 4;
                }
                switch (mode) {
                  case GL_TRIANGLES: {
                    for (step = 0; step < nverts; step += 3) {
                      if (colorVals_tmp) {
                        c0 = colorVals_tmp;
                        c1 = colorVals_tmp + 4;
                        c2 = colorVals_tmp + 8;
                      } else {
                        c1 = c2 = c0;
                      }
                      if (normalVals_tmp) {
                        n0 = normalVals_tmp;
                        n1 = normalVals_tmp + 3;
                        n2 = normalVals_tmp + 6;
                      } else {
                        n1 = n2 = n0;
                      }
                      CGOAlphaTriangle(info->alpha_cgo, vertexVals_tmp,
                                       vertexVals_tmp + 3, vertexVals_tmp + 6,
                                       n0, n1, n2, c0, c1, c2, R->alpha,
                                       R->alpha, R->alpha, false);
                      vertexVals_tmp += 9;
                      if (normalVals_tmp) {
                        normalVals_tmp += 9;
                      }
                      if (colorVals_tmp) {
                        colorVals_tmp += 12;
                      }
                    }
                  } break;
                  case GL_TRIANGLE_STRIP: {
                    if (colorVals_tmp) {
                      c1 = colorVals_tmp;
                      c2 = colorVals_tmp + 4;
                      colorVals_tmp += 8;
                    } else {
                      c1 = c2 = c0;
                    }
                    if (normalVals_tmp) {
                      n1 = normalVals_tmp;
                      n2 = normalVals_tmp + 3;
                      normalVals_tmp += 6;
                    } else {
                      n1 = n2 = n0;
                    }
                    vertexVals_tmp += 6;
                    for (step = 2; step < nverts; step++) {
                      if (colorVals_tmp) {
                        c0 = c1;
                        c1 = c2;
                        c2 = colorVals_tmp;
                      }
                      if (normalVals_tmp) {
                        n0 = n1;
                        n1 = n2;
                        n2 = normalVals_tmp;
                      }
                      CGOAlphaTriangle(info->alpha_cgo, vertexVals_tmp - 6,
                                       vertexVals_tmp - 3, vertexVals_tmp, n0,
                                       n1, n2, c0, c1, c2, R->alpha, R->alpha,
                                       R->alpha, false);
                      vertexVals_tmp += 3;
                      if (normalVals_tmp) {
                        normalVals_tmp += 3;
                      }
                      if (colorVals_tmp) {
                        colorVals_tmp += 4;
                      }
                    }
                  } break;
                  case GL_TRIANGLE_FAN: {
                    float* firstVertex = vertexVals_tmp;
                    if (colorVals_tmp) {
                      c0 = colorVals_tmp;
                      c2 = colorVals_tmp + 4;
                      colorVals_tmp += 8;
                    } else {
                      c1 = c2 = c0;
                    }
                    if (normalVals_tmp) {
                      n0 = normalVals_tmp;
                      n2 = normalVals_tmp + 3;
                      normalVals_tmp += 6;
                    }
                    vertexVals_tmp += 6;
                    for (step = 2; step < nverts; step++) {
                      if (colorVals_tmp) {
                        c1 = c2;
                        c2 = colorVals_tmp;
                      }
                      if (normalVals_tmp) {
                        n1 = n2;
                        n2 = normalVals_tmp;
                      }
                      CGOAlphaTriangle(info->alpha_cgo, firstVertex,
                                       vertexVals_tmp - 3, vertexVals_tmp, n0,
                                       n1, n2, c0, c1, c2, R->alpha, R->alpha,
                                       R->alpha, false);
                      vertexVals_tmp += 3;
                      if (normalVals_tmp) {
                        normalVals_tmp += 3;
                      }
                      if (colorVals_tmp) {
                        colorVals_tmp += 4;
                      }
                    }
                  } break;
                }
              } break;
              case CGO_VERTEX:
                v0 = pc;
                switch (mode) {
                  case GL_TRIANGLES:
                    if (3 * ((vc + 1) / 3) == vc + 1) {
                      CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2,
                                       c0, c1, c2, R->alpha, R->alpha, R->alpha,
                                       true);
                    }
                    v2 = v1;
                    c2 = c1;
                    n2 = n1;
                    v1 = v0;
                    c1 = c0;
                    n1 = n0;
                    vc++;
                    break;
                  case GL_TRIANGLE_STRIP:
                    if (vc > 1) {
                      CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2,
                                       c0, c1, c2, R->alpha, R->alpha, R->alpha,
                                       !(vc & 0x1));
                    }
                    v2 = v1;
                    c2 = c1;
                    n2 = n1;
                    v1 = v0;
                    c1 = c0;
                    n1 = n0;
                    vc++;
                    break;
                  case GL_TRIANGLE_FAN:
                    if (vc > 1) {
                      CGOAlphaTriangle(info->alpha_cgo, v0, v1, v2, n0, n1, n2,
                                       c0, c1, c2, R->alpha, R->alpha, R->alpha,
                                       false);
                    } else if (!vc) {
                      n2 = n0;
                      v2 = v0;
                      c2 = c0;
                    }
                    v1 = v0;
                    c1 = c0;
                    n1 = n0;
                    vc++;
                    break;
                  default:
                    CGO_gl_vertex(R, &pc);
                    break;
                }
                break;
              default:
                CGO_gl[op](R, &pc);
                break;
            }
          } else { /* opaque */
            switch (op) {
              case CGO_COLOR:
                /* Since CGO operations are done in sequence, alpha could happen
                   after color is set.  In this case, we still need to keep
                   track of the color in case there is a transparent object */
                c0 = pc;
                break;
              default:
                break;
            }
            CGO_gl[op](R, &pc);
          }
        }
      }
    }
  }
}

void CGORenderGLAlpha(CGO* I, RenderInfo* info, bool calcDepth) {
  PyMOLGlobals* G = I->G;
  if (G->ValidContext && I->c) {
    int mode = GL_TRIANGLES;
    if (I->debug) {
      mode = CGOConvertDebugMode(I->debug, GL_TRIANGLES);
    }
    /* 1. transform and measure range (if not already known)
       2. bin into linked lists based on Z-centers
       3. render by layer */

    if (I->z_flag) {
      if (!I->i_start) {
        I->i_size = 256;
        I->i_start = pymol::calloc<int>(I->i_size);
      } else {
        UtilZeroMem(I->i_start, sizeof(int) * I->i_size);
      }
      {
        const int i_size = I->i_size;
        const float* base = I->op;
        int* start = I->i_start;
        int delta = 1, ntris = 0;
        /* bin the triangles */
        if (calcDepth) {
          for (auto it = I->begin(); !it.is_stop(); ++it) {
            if (it.op_code() == CGO_ALPHA_TRIANGLE) {
              float* const pc = it.data();
              const float z = dot_product3f(pc + 1, I->z_vector);
              if (z > I->z_max) I->z_max = z;
              if (z < I->z_min) I->z_min = z;
              pc[4] = z;
              ntris++;
            }
          }
        }

        const float range_factor = (0.9999F * i_size) / (I->z_max - I->z_min);

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          if (it.op_code() == CGO_ALPHA_TRIANGLE) {
            float* const pc = it.data();
            assert(base < pc && pc < I->op + I->c);
            auto i =
                std::clamp<int>((pc[4] - I->z_min) * range_factor, 0, i_size);
            CGO_put_int(pc, start[i]);
            start[i] = (pc - base); /* NOTE: will always be > 0 since we have
                                       CGO_read_int'd */
          }
        }

        // for single-layer transparency, render front-to-back
        if (SettingGetGlobal_i(G, cSetting_transparency_mode) == 2) {
          delta = -1;
          start += (i_size - 1);
        }

        /* now render by bin */
      }
    } else {
    }
  }
}
