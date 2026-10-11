
/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

// Helpers shared by the cgo*.cc files; not part of the CGO API.

#pragma once

#include "pymol/graphics/cgo.h"
#include "pymol/graphics/gl/gl_vertex_buffer.h"

#define VERTEX_NORMAL_SIZE 3

constexpr unsigned VERTEX_PICKCOLOR_RGBA_SIZE = 1;  // 4 unsigned bytes

constexpr unsigned VERTEX_PICKCOLOR_INDEX_SIZE = 2;  // index + bond

constexpr unsigned VERTEX_PICKCOLOR_SIZE = VERTEX_PICKCOLOR_RGBA_SIZE +  //
                                           VERTEX_PICKCOLOR_INDEX_SIZE;

template <typename T>
inline T CLAMPVALUE(T val, T minimum, T maximum) {
  return (val < minimum) ? minimum : (val > maximum) ? maximum : val;
}

inline cCylCap cap1_from_cyl_shader_bits(const unsigned char bits) {
  return (bits & 1) ? ((bits & cCylShaderCap1RoundBit) ? cCylCap::Round
                                                       : cCylCap::Flat)
                    : cCylCap::None;
}

inline cCylCap cap2_from_cyl_shader_bits(const unsigned char bits) {
  return (bits & 2) ? ((bits & cCylShaderCap2RoundBit) ? cCylCap::Round
                                                       : cCylCap::Flat)
                    : cCylCap::None;
}

#define WARN_UNEXPECTED_OPERATION(G, op) \
  PRINTFB(G, FB_CGO, FB_Warnings)        \
  " %s-Warning: unexpected op=0x%x (line %d)\n", __func__, op, __LINE__ ENDFB(G)

// like g_return_val_if_fail from glib
#define RETURN_VAL_IF_FAIL(expr, val) \
  {                                   \
    if (!(expr)) return (val);        \
  }

// Defined in cgo_simplify.cc.
int CGOSimpleCylinder(CGO* I, const float* v1, const float* v2,
                      const float tube_size, const float* c1, const float* c2,
                      float a1, const float a2, const bool interp,
                      const cCylCap cap1, const cCylCap cap2,
                      const Pickable* pickcolor2 = nullptr,
                      const bool stick_round_nub = false);
int CGOSimpleCone(CGO* I, const float* v1, const float* v2, float r1, float r2,
                  const float* c1, const float* c2, cCylCap cap1, cCylCap cap2);
int CGOSimpleCone(CGO* I, const cgo::draw::cone& cone);

template <typename CylinderT>
int CGOSimpleCylinder(CGO* I, const CylinderT& cyl, const float a1,
                      const float a2, const bool interp, const cCylCap cap1,
                      const cCylCap cap2, const Pickable* pickcolor2 = nullptr,
                      const bool stick_round_nub = false) {
  return CGOSimpleCylinder(I, cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
                           cyl.color2, a1, a2, interp, cap1, cap2, pickcolor2,
                           stick_round_nub);
}

// Defined in cgo.cc.
bool CGOPackedNormals(PyMOLGlobals* G);

// Defined in cgo_optimize.cc.
VertexBuffer& CreateVertexBuffer(
    PyMOLGlobals* G, VertexBufferLayout layout = VertexBufferLayout::Separate,
    MemoryUsageProperty memProperty = MemoryUsageProperty::GpuOnly);
void SetVertexBufferData(PyMOLGlobals* G, VertexBuffer& vertexBuffer,
                         BufferDataDesc&& dataDesc,
                         pymol::span<const std::byte> dataOverride = {});
void SetVertexBufferDataSized(PyMOLGlobals* G, VertexBuffer& vertexBuffer,
                              BufferDataDesc&& dataDesc, std::size_t bufferSize,
                              std::byte* dataOverride = nullptr);
IndexBuffer& CreateIndexBuffer(PyMOLGlobals* G);
void SetIndexBufferData(PyMOLGlobals* G, IndexBuffer& indexBuffer,
                        pymol::span<const std::uint32_t> data);

// Defined in cgo_simplify.cc.
int CGOQuadricToEllipsoid(const float* v, float r, const float* q, float* r_el,
                          float* n0, float* n1, float* n2);
