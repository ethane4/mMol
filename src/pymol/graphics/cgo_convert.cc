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

// Converting CGOs for the shaders: lines, cylinders, crosses, labels and
// points.

#include "pymol/graphics/cgo.h"
#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/scene/scene.h"
#include "pymol/graphics/graphics_util.h"
#include <iostream>
#include "pymol/graphics/cgo_internal.h"

inline unsigned char cyl_shader_bits_from_caps(const cCylCap cap1,
                                               const cCylCap cap2) {
  return ((cap1 == cCylCap::Flat)
              ? cCylShaderCap1Flat
              : ((cap1 == cCylCap::Round) ? cCylShaderCap1Round : 0)) |
         ((cap2 == cCylCap::Flat)
              ? cCylShaderCap2Flat
              : ((cap2 == cCylCap::Round) ? cCylShaderCap2Round : 0));
}

CGO* CGOConvertToLabelShader(const CGO* I, CGO* addTo) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  AttribDataOp world_pos_op = {{CGO_DRAW_LABEL, 1, FLOAT3_TO_FLOAT3,
                                offsetof(cgo::draw::label, world_pos), 0}};
  AttribDataOp screen_offset_op = {
      {CGO_DRAW_LABEL, 2, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::label, screen_world_offset), 0}};
  AttribDataOp screen_min_op = {{CGO_DRAW_LABEL, 3, FLOAT3_TO_FLOAT3,
                                 offsetof(cgo::draw::label, screen_min), 0}};
  AttribDataOp screen_max_op = {{CGO_DRAW_LABEL, 4, FLOAT3_TO_FLOAT3,
                                 offsetof(cgo::draw::label, screen_max), 0}};
  AttribDataOp text_extent_op = {{CGO_DRAW_LABEL, 5, FLOAT2_TO_FLOAT2,
                                  offsetof(cgo::draw::label, text_extent), 0}};
  AttribDataOp relative_mode_op = {{CGO_DRAW_LABEL, 6, FLOAT_TO_FLOAT,
                                    offsetof(cgo::draw::label, relative_mode),
                                    0}};
  AttribDataOp target_pos_op = {{CGO_DRAW_LABEL, 7, FLOAT3_TO_FLOAT3,
                                 offsetof(cgo::draw::label, target_pos), 6}};

  AttribDataDesc attrDesc = {
      {"attr_worldpos", VertexFormat::Float3, world_pos_op},
      {"attr_targetpos", VertexFormat::Float3, target_pos_op},
      {"attr_screenoffset", VertexFormat::Float3, screen_min_op},
      {"attr_texcoords", VertexFormat::Float2, text_extent_op},
      {"attr_screenworldoffset", VertexFormat::Float3, screen_offset_op},
      {"attr_relative_mode", VertexFormat::Float, relative_mode_op}};

  auto ComputeScreenValues = [](void* varData, const float* pc,
                                void* screenData, int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label*>(pc);
    const auto& smin = sp->screen_min;
    const auto& smax = sp->screen_max;
    float* v = reinterpret_cast<float*>(varData);
    switch (idx) {
      case 0:
        v[0] = smin[0];
        v[1] = smin[1];
        v[2] = smin[2];
        break;
      case 1:
        v[0] = smin[0];
        v[1] = smax[1];
        v[2] = smin[2];
        break;
      case 2:
        v[0] = smax[0];
        v[1] = smin[1];
        v[2] = smin[2];
        break;
      case 3:
        v[0] = smin[0];
        v[1] = smax[1];
        v[2] = smin[2];
        break;
      case 4:
        v[0] = smax[0];
        v[1] = smin[1];
        v[2] = smin[2];
        break;
      case 5:
        v[0] = smax[0];
        v[1] = smax[1];
        v[2] = smin[2];
        break;
    };
  };

  auto ComputeTexCoords = [](void* varData, const float* pc, void* discard,
                             int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label*>(pc);
    float* v = reinterpret_cast<float*>(varData);
    const auto& te = sp->text_extent;
    static struct {
      int x, y;
    } const idxs[6] = {{0, 1}, {0, 3}, {2, 1}, {0, 3}, {2, 1}, {2, 3}};
    v[0] = te[idxs[idx].x];
    v[1] = te[idxs[idx].y];
  };

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {ComputeScreenValues, nullptr, "attr_screenoffset"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {ComputeTexCoords, nullptr, "attr_texcoords"});

  uchar pickdata[4] = {0, 0, 0, 0};
  addTo->add<cgo::draw::vertex_attribute_4ub>(
      G->ShaderMgr->GetAttributeUID("attr_pickcolor"), pickdata);

  AttribDataOp pickOp = {{CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"attr_pickcolor", VertexFormat::UByte4Norm, pickOp}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved, true);
}

CGO* CGOExpandDrawTextures(const CGO* I, int est) {
  CGO* cgo = CGONew(I->G);
  int ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeTextures() CGO_DRAW_BUFFERS_INDEXED or "
        "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_TEXTURE: {
        float screenMin[3], screenMax[3], textExtent[4];
        float alpha = cgo->alpha;
        CGOAlpha(cgo, 0.f);
        CGOColor(cgo, 0.f, 0.f, 0.f);
        copy3f(pc + 3, screenMin);
        copy3f(pc + 6, screenMax);
        copy4f(pc + 9, textExtent);
        CGOBegin(cgo, GL_TRIANGLES);
        CGOTexCoord2f(cgo, textExtent[0], textExtent[1]);
        CGOVertexv(cgo, screenMin);
        CGOTexCoord2f(cgo, textExtent[0], textExtent[3]);
        CGOVertex(cgo, screenMin[0], screenMax[1], screenMin[2]);
        CGOTexCoord2f(cgo, textExtent[2], textExtent[1]);
        CGOVertex(cgo, screenMax[0], screenMin[1], screenMin[2]);
        CGOTexCoord2f(cgo, textExtent[0], textExtent[3]);
        CGOVertex(cgo, screenMin[0], screenMax[1], screenMin[2]);
        CGOTexCoord2f(cgo, textExtent[2], textExtent[1]);
        CGOVertex(cgo, screenMax[0], screenMin[1], screenMin[2]);
        CGOTexCoord2f(cgo, textExtent[2], textExtent[3]);
        CGOVertex(cgo, screenMax[0], screenMax[1], screenMin[2]);
        CGOEnd(cgo);
        CGOAlpha(cgo, alpha);
      } break;
      default:
        cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  CGOStop(cgo);
  return cgo;
}

CGO* CGOConvertLinesToShaderCylinders(const CGO* I, int est) {
  int tot_nverts = 0, tot_ncyls = 0;

  CGO* cgo = CGONewSized(I->G, I->c + est);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
      case CGO_DRAW_ARRAYS: {
        const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        float* vals =
            cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
        int nvals = sp->narrays * sp->nverts;
        memcpy(vals, sp->floatdata, nvals);
      } break;
      case CGO_END:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOConvertLinesToShaderCylinders: CGO_END encountered without "
        "CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);
        break;
      case CGO_VERTEX:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOConvertLinesToShaderCylinders: CGO_VERTEX encountered without "
        "CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);
        break;
      case CGO_BEGIN: {
        const float *last_vertex = nullptr, *last_color = nullptr,
                    *current_color = nullptr, *color = nullptr;
        unsigned int last_pick_color_idx = 0;
        int last_pick_color_bnd = cPickableNoPick;
        int nverts = 0, err = 0;
        int mode = CGO_get_int(pc);

        for (++it; !err && it != CGO_END; ++it) {
          auto pc = it.data();
          int op = it.op_code();

          switch (op) {
            case CGO_VERTEX:
              if (last_vertex) {
                switch (mode) {
                  case GL_LINES:
                  case GL_LINE_STRIP: {
                    float axis[3];
                    bool pick_color_diff = false;
                    axis[0] = pc[0] - last_vertex[0];
                    axis[1] = pc[1] - last_vertex[1];
                    axis[2] = pc[2] - last_vertex[2];
                    pick_color_diff =
                        (cgo->current_pick_color_index != last_pick_color_idx ||
                         cgo->current_pick_color_bond != last_pick_color_bnd);
                    if (last_color && current_color &&
                        (!equal3f(last_color, current_color) ||
                         pick_color_diff)) {
                      CGOColorv(cgo, last_color);
                      if (pick_color_diff) {
                        Pickable pickcolor2 = {cgo->current_pick_color_index,
                                               cgo->current_pick_color_bond};
                        CGOPickColor(cgo, last_pick_color_idx,
                                     last_pick_color_bnd);
                        cgo->add<cgo::draw::shadercylinder2ndcolor>(
                            cgo, last_vertex, axis, 1.f,
                            cCylShaderBothCapsRound, current_color,
                            &pickcolor2);
                        CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
                      } else {
                        cgo->add<cgo::draw::shadercylinder2ndcolor>(
                            cgo, last_vertex, axis, 1.f,
                            cCylShaderBothCapsRound, current_color);
                      }
                      CGOColorv(cgo, current_color);
                    } else {
                      cgo->add<cgo::draw::shadercylinder>(
                          last_vertex, axis, 1.f, cCylShaderBothCapsRound);
                    }
                    last_vertex = pc;
                    last_pick_color_idx = cgo->current_pick_color_index;
                    last_pick_color_bnd = cgo->current_pick_color_bond;
                    tot_ncyls++;
                  }
                    if (mode == GL_LINES) {
                      last_vertex = nullptr;
                      last_color = nullptr;
                    }
                }
              } else {
                last_vertex = pc;
                current_color = color;
                last_pick_color_idx = cgo->current_pick_color_index;
                last_pick_color_bnd = cgo->current_pick_color_bond;
              }
              nverts++;
              break;
            case CGO_LINE: {
              float axis[3];
              auto line = reinterpret_cast<const cgo::draw::line*>(pc);
              subtract3f(line->vertex2, line->vertex1, axis);
              cgo->add<cgo::draw::shadercylinder>(line->vertex1, axis, 1.f,
                                                  cCylShaderBothCapsRound);
              tot_ncyls++;
            } break;
            case CGO_SPLITLINE: {
              float axis[3];
              auto splitline =
                  reinterpret_cast<const cgo::draw::splitline*>(pc);
              Pickable pickcolor2 = {splitline->index, splitline->bond};
              float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
                                CONVERT_COLOR_VALUE(splitline->color2[1]),
                                CONVERT_COLOR_VALUE(splitline->color2[2])};
              unsigned char flags = splitline->flags;
              subtract3f(splitline->vertex2, splitline->vertex1, axis);
              if ((flags & cgo::draw::splitline::equal_colors) &&
                  (flags & cgo::draw::splitline::no_split_for_pick)) {
                cgo->add<cgo::draw::shadercylinder>(
                    splitline->vertex1, axis, 1., cCylShaderBothCapsRound);
              } else {
                int cap = cCylShaderBothCapsRound;
                if (flags & splitline->flags &
                    cgo::draw::splitline::interpolation) {
                  cap |= cCylShaderInterpColor;
                }
                cgo->add<cgo::draw::shadercylinder2ndcolor>(
                    cgo, splitline->vertex1, axis, 1., cap, color2,
                    &pickcolor2);
                last_pick_color_idx = splitline->index;
                last_pick_color_bnd = splitline->bond;
              }
              tot_ncyls++;
            } break;
            case CGO_COLOR:
              if (op == CGO_COLOR) {
                last_color = current_color;
                current_color = pc;
                color = pc;
              }
            case CGO_PICK_COLOR:
              if (op == CGO_PICK_COLOR) {
                cgo->current_pick_color_index = CGO_get_uint(pc);
                cgo->current_pick_color_bond = CGO_get_int(pc + 1);
              }
            default:
              cgo->add_to_cgo(op, pc);
          }
        }

        tot_nverts += nverts;
      } break;
      default:
        cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
  }
  if (tot_ncyls) {
    return (cgo);
  } else {
    CGOFree(cgo);
    return nullptr;
  }
}

static void trilinesBufferAddVertex(
    float*& buffer,
    const float* v1,     // vertex
    const float* v2,     // vertex other end of line
    const float* color,  // RGB color
    float alpha,         // alpha
    signed char uv)      // uv
{
  // vertex
  (*buffer++) = v1[0];
  (*buffer++) = v1[1];
  (*buffer++) = v1[2];

  // othervertex
  (*buffer++) = v2[0];
  (*buffer++) = v2[1];
  (*buffer++) = v2[2];

  (*buffer++) = (float)uv;

  // RGBA
  unsigned char* byte_view = (unsigned char*)(buffer++);
  (*byte_view++) = CLIP_COLOR_VALUE(color[0]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[1]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[2]);
  (*byte_view++) = CLIP_COLOR_VALUE(alpha);
}

static void trilinesBufferAddVertices(
    float*& buffer,
    const float* v1,     // vertex
    const float* v2,     // vertex other end of line
    const float* color,  // RGB color
    float alpha)         // alpha
{
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1);  //-1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3);  // 1, 1);
  // Vertex 2
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 0);  //-1, -1);
  // Vertex 4
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3);  // 1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 2);  // 1, -1);
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1);  //-1, 1);
}

static void CGOTrilines_GetCurrentColor(const float*& current_color,
                                        const float* colorv,
                                        const float* last_color,
                                        const float* cc) {
  if (!current_color) {
    if (colorv) {
      current_color = colorv;
    } else if (last_color) {
      current_color = last_color;
    } else {
      current_color = cc;
    }
  }
}

/**
 * Changes
 *
 *     CGO_ENABLE <frommode>
 *
 * to
 *
 *     CGO_ENABLE <tomode>
 *
 * in place.
 */
void CGOChangeShadersTo(CGO* I, int frommode, int tomode) {
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_ENABLE) {
      auto eo = it.cast<cgo::draw::enable>();
      if (eo->mode == frommode) {
        eo->mode = tomode;
      }
    }
  }
}

CGO* CGOConvertSpheresToPoints(const CGO* I) {
  CGO* cgo;

  int ok = true;
  cgo = CGONew(I->G);
  CHECKOK(ok, cgo);
  CGOBegin(cgo, GL_POINTS);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        CGOPickColor(cgo, cgo->current_pick_color_index,
                     cgo->current_pick_color_bond);
        break;
      case CGO_SHADER_CYLINDER:
      case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
      case CGO_CYLINDER:
      case CGO_CONE:
      case CGO_SAUSAGE:
      case CGO_CUSTOM_CYLINDER:
      case CGO_CUSTOM_CYLINDER_ALPHA:
      case CGO_END:
      case CGO_VERTEX:
      case CGO_BEGIN:
      case CGO_ELLIPSOID:
      case CGO_QUADRIC:
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
      case CGO_DRAW_SPHERE_BUFFERS:
      case CGO_DRAW_CYLINDER_BUFFERS:
      case CGO_DRAW_LABELS:
        break;
      case CGO_SPHERE:
        CGOVertexv(cgo, pc);
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
      default:
        cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  CGOEnd(cgo);
  if (ok) {
    ok &= CGOStop(cgo);
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

// Will create an interleaved VBO with { vertex, otherVertex, uv, and texcoord
// info } Currently, this function does not support/parse lines inside
// CGODrawArrays, i.e., a CGO that CGOCombineBeginEnd was used on
CGO* CGOConvertLinesToTrilines(const CGO* I, bool addshaders) {
  static std::set<int> lineops = {CGO_VERTEX, CGO_LINE, CGO_SPLITLINE};
  auto G = I->G;
  const int nLines = CGOCountNumberOfOperationsOfTypeN(I, lineops) + 1;

  if (nLines == 0) {
    return nullptr;
  }

  int line_counter = 0;
  GLuint glbuff = 0;
  const float* colorv = nullptr;
  unsigned int buff_size = nLines * 6 * (8 * sizeof(float));

  // VBO memory -- number of lines x 4 vertices per line x (vertex + otherVertex
  // + normal + texCoord + color)
  std::vector<float> buffer_start(buff_size);
  float* buffer = buffer_start.data();

  std::unique_ptr<CGO> cgo(new CGO(G));

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_DRAW_ARRAYS: {
        auto sp = it.cast<cgo::draw::arrays>();
        float* vals =
            cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
        int nvals = sp->narrays * sp->nverts;
        memcpy(vals, sp->floatdata, nvals);
      } break;
      case CGO_END:
        WARN_UNEXPECTED_OPERATION(G, op);
        return nullptr;
      case CGO_BEGIN: {
        const float *last_vertex = nullptr, *last_color = nullptr,
                    *current_color = nullptr, *color = nullptr;
        const int mode = it.cast<cgo::draw::begin>()->mode;

        for (++it;; ++it) {
          if (it.is_stop()) {
            WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
            return nullptr;
          }

          const auto op = it.op_code();
          if (op == CGO_END) {
            break;
          }

          const auto pc = it.data();

          switch (op) {
            case CGO_VERTEX:
              if (last_vertex) {
                switch (mode) {
                  case GL_LINES:
                  case GL_LINE_STRIP: {
                    float cc[3] = {1, 1, 1};
                    float alpha = cgo->alpha;
                    CGOTrilines_GetCurrentColor(current_color, colorv,
                                                last_color, cc);
                    trilinesBufferAddVertices(buffer, pc, last_vertex,
                                              current_color, alpha);
                    line_counter++;
                    last_vertex = pc;
                  }
                    if (mode == GL_LINES) {
                      last_vertex = nullptr;
                      last_color = nullptr;
                    }
                }
              } else {
                last_vertex = pc;
                current_color = color;
              }
              break;
            case CGO_LINE: {
              auto line = it.cast<cgo::draw::line>();
              float cc[3] = {1, 1, 1};
              float alpha = cgo->alpha;
              CGOTrilines_GetCurrentColor(current_color, colorv, last_color,
                                          cc);
              trilinesBufferAddVertices(buffer, line->vertex1, line->vertex2,
                                        current_color, alpha);
              line_counter++;
            } break;
            case CGO_SPLITLINE: {
              auto splitline = it.cast<cgo::draw::splitline>();
              float cc[3] = {1, 1, 1};
              float alpha = cgo->alpha;
              float mid[3];
              float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
                                CONVERT_COLOR_VALUE(splitline->color2[1]),
                                CONVERT_COLOR_VALUE(splitline->color2[2])};
              add3f(splitline->vertex1, splitline->vertex2, mid);
              mult3f(mid, .5f, mid);
              CGOTrilines_GetCurrentColor(current_color, colorv, last_color,
                                          cc);
              trilinesBufferAddVertices(buffer, splitline->vertex1, mid,
                                        current_color, alpha);
              trilinesBufferAddVertices(buffer, mid, splitline->vertex2, color2,
                                        alpha);
              line_counter += 2;
            } break;
            case CGO_COLOR:
              last_color = current_color;
              current_color = pc;
              color = pc;
              break;
          }
        }
      } break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_COLOR:
        colorv = pc;
        break;
    }
  }

  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
  }

  {
    /* TODO: Rmove GL calls and use Generic Vertex Buffers */
    glGenBuffers(1, &glbuff);
    glBindBuffer(GL_ARRAY_BUFFER, glbuff);
    glBufferData(GL_ARRAY_BUFFER, line_counter * 6 * 8 * sizeof(float),
                 buffer_start.data(), GL_STATIC_DRAW);
    CheckGLErrorOK(G,
                   "ERROR: CGOConvertLinesToTriangleStrips() glBindBuffer "
                   "returns err=%d\n");
    if (addshaders) cgo->add<cgo::draw::enable>(GL_TRILINES_SHADER);
    cgo->add<cgo::draw::trilines>(line_counter * 6, glbuff);
    cgo->has_draw_buffers = true;
    if (addshaders) cgo->add<cgo::draw::disable>(GL_TRILINES_SHADER);
    CGOStop(cgo.get());
  }

  return cgo.release();
}

/**
 * copies data for a particular attribute operation into the array used to load
 * the VBO. this takes into account whether it is interleaved or not.
 *
 * isInterleaved    : whether the VBO is interleaved
 * nvert            : which vertex in the VBO
 * attribOp         : the attribute op
 * vertexDataSize   : total vertex data size in VBO (for interleaved)
 * dataPtrs         : all data pointers for attributes (for interleaved, they
 * are all the pointer to the one array) attrOffset       : offsets of the
 * attributes (for interleaved) pcarg            : pc pointer to CGO operation
 * data pick_data        : pointer to pick data for current vertex (writes to if
 * pick data) has_pick_colorBS : keeps track of which pick attributes have been
 * set
 *
 */
static void copyAttributeForOp(bool isInterleaved, int& nvert,
                               AttribOp* attribOp, int vertexDataSize,
                               std::vector<void*>& dataPtrs,
                               std::vector<int>& attrOffset, const float* pcarg,
                               float* pick_data, int& has_pick_colorBS,
                               int pstride) {
  auto attrDesc = attribOp->desc;
  int ord = attrDesc->order;
  int copyord = -1;
  void* dataPtr = dataPtrs[ord];
  unsigned char* pc = ((unsigned char*)pcarg) + attribOp->offset;
  if (isInterleaved) {
    dataPtr =
        (unsigned char*)dataPtr + nvert * vertexDataSize + attrOffset[ord];
    if (attribOp->copyAttribDesc) {
      copyord = attribOp->copyAttribDesc->order;
      pc = ((unsigned char*)dataPtrs[ord]) + nvert * vertexDataSize +
           attrOffset[copyord];
    }
  } else {
    auto sz = GetSizeOfVertexFormat(attrDesc->m_format);
    dataPtr = (unsigned char*)dataPtr + nvert * sz;
    if (attribOp->copyAttribDesc) {
      copyord = attribOp->copyAttribDesc->order;
      auto copysz = GetSizeOfVertexFormat(attribOp->copyAttribDesc->m_format);
      pc = (unsigned char*)dataPtr + nvert * copysz;
    }
  }
  switch (attribOp->conv_type) {
    case NO_COPY:
      break;
    case FLOAT_TO_FLOAT:
      *((float*)dataPtr) = *((float*)pc);
      break;
    case FLOAT2_TO_FLOAT2:
      *((float*)dataPtr) = *((float*)pc);
      *((float*)dataPtr + 1) = *((float*)pc + 1);
      break;
    case FLOAT3_TO_FLOAT3:
      copy3f((float*)pc, (float*)dataPtr);
      break;
    case FLOAT4_TO_FLOAT4:
      *((float*)dataPtr) = *((float*)pc);
      *((float*)dataPtr + 1) = *((float*)pc + 1);
      *((float*)dataPtr + 2) = *((float*)pc + 2);
      *((float*)dataPtr + 3) = *((float*)pc + 3);
      break;
    case FLOAT3_TO_UB3: {
      auto dataPtrUB = (unsigned char*)dataPtr;
      float* pcf = (float*)pc;
      dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
      dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
      dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
    } break;
    case FLOAT1_TO_UB_4TH: {
      auto dataPtrUB = (unsigned char*)dataPtr;
      float* pcf = (float*)pc;
      dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[0]);
    } break;
    case UB3_TO_UB3: {
      auto dataPtrUB = (unsigned char*)dataPtr;
      auto pcUB = (unsigned char*)pc;
      dataPtrUB[0] = pcUB[0];
      dataPtrUB[1] = pcUB[1];
      dataPtrUB[2] = pcUB[2];
      break;
    }
    case UINT_INT_TO_PICK_DATA: {
      float* pcf = (float*)pc;
      unsigned int index = CGO_get_uint(pcf);
      int bond = CGO_get_int(pcf + 1);
      CGO_put_uint(ord * 2 + pick_data, index);
      CGO_put_int(ord * 2 + pick_data + 1, bond);
      has_pick_colorBS |= (1 << ord);
    } break;
    case FLOAT4_TO_UB4: {
      auto dataPtrUB = (unsigned char*)dataPtr;
      float* pcf = (float*)pc;
      dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
      dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
      dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
      dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[3]);
    } break;
    case CYL_CAP_TO_CAP: {
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = *pc;
    } break;
    case CYL_CAPS_ARE_ROUND: {
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = cCylShaderBothCapsRound | cCylShaderInterpColor;
    } break;
    case CYL_CAPS_ARE_FLAT: {
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = cCylShaderBothCapsFlat | cCylShaderInterpColor;
    } break;
    case CYL_CAPS_ARE_CUSTOM: {
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      float* pcf = (float*)pc;
      const cCylCap cap1 = static_cast<cCylCap>(int(pcf[0]));
      const cCylCap cap2 = static_cast<cCylCap>(int(pcf[1]));
      dataPtrUB[0] =
          cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
    } break;
    case UB1_TO_INTERP: {
      bool interp = (pc[0] & cgo::draw::splitline::interpolation);
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = interp ? 1 : 0;
    } break;
    case UB1_INTERP_TO_CAP: {
      bool interp = (pc[0] & cgo::draw::splitline::interpolation);
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] =
          (cCylShaderBothCapsRound | (interp ? cCylShaderInterpColor : 0));
    } break;
    case FLOAT1_TO_INTERP: {
      float interp = *((float*)pc);
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = (interp > .5f) ? 1 : 0;
    } break;
    case FLOAT1_INTERP_TO_CAP: {
      float interp = *((float*)pc);
      unsigned char* dataPtrUB = (unsigned char*)dataPtr;
      dataPtrUB[0] = (cCylShaderBothCapsRound |
                      ((interp > .5f) ? cCylShaderInterpColor : 0));
    } break;
    case UB4_TO_UB4: {
      auto dataPtrUB = (unsigned char*)dataPtr;
      auto pcUB = (unsigned char*)pc;
      dataPtrUB[0] = pcUB[0];
      dataPtrUB[1] = pcUB[1];
      dataPtrUB[2] = pcUB[2];
      dataPtrUB[3] = pcUB[3];
      break;
    }
    case PICK_DATA_TO_PICK_DATA: {
      float* pcf;
      if (copyord < 0) {
        pcf = (float*)pc;
      } else {
        pcf = (copyord * 2 + pick_data);
        if (nvert) {
          pcf -= pstride;
        }
      }
      unsigned int index = CGO_get_uint(pcf);
      int bond = CGO_get_int(pcf + 1);
      CGO_put_uint(ord * 2 + pick_data, index);
      CGO_put_int(ord * 2 + pick_data + 1, bond);
      has_pick_colorBS |= (1 << ord);
      break;
    }
  }
}

/**
 * copies data for a particular attribute into the array used to load the VBO.
 * this takes into account whether it is interleaved or not, and if an attribute
 * has repeat values
 *
 * isInterleaved  : whether the VBO is interleaved
 * nvert          : which vertex in the VBO
 * attribDesc     : the attribute description
 * vertexDataSize : total vertex data size in VBO (for interleaved)
 * dataPtrs       : all data pointers for attributes (for interleaved, they are
 * all the pointer to the one array) attrOffset     : offsets of the attributes
 * (for interleaved)
 *
 */
static void copyAttributeForVertex(bool isInterleaved, int& nvert,
                                   AttribDesc& attribDesc,
                                   const int vertexDataSize,
                                   std::vector<void*>& dataPtrs,
                                   std::vector<int>& attrOffset) {
  int ord = attribDesc.order;
  void* dataPtr = dataPtrs[ord];
  unsigned char* pc = nullptr;
  auto attrSize = GetSizeOfVertexFormat(attribDesc.m_format);
  if (isInterleaved) {
    dataPtr =
        (unsigned char*)dataPtr + nvert * vertexDataSize + attrOffset[ord];
    pc = (unsigned char*)dataPtr - vertexDataSize;
  } else {
    dataPtr = (unsigned char*)dataPtr + nvert * attrSize;
    pc = (unsigned char*)dataPtr - attrSize;
  }
  if (attribDesc.repeat_value && attribDesc.repeat_value_length) {
    int pos = (nvert % attribDesc.repeat_value_length);
    pc = attribDesc.repeat_value + pos * attrSize;
    memcpy(dataPtr, pc, attrSize);
  } else {
    memcpy(dataPtr, pc, attrSize);
  }
}
/**
 * check all attributes (pick and non-pick) to see if they are specified in the
 * CGO (I) also checks to see if any picking is specified and sets has_picking
 * argument if any of the attributes are not specified and if a default_value is
 * set (in the AttribDesc) then the associated vertex_attribute CGO operation is
 * inserted into the cgo that is passed in.
 *
 * I:           primitive CGO that is processed
 * attrData:    definition of attributes
 * pickData:    definition of pick attributes
 * cgo:         new cgo that could have vertex_attribute CGO operations added
 * has_picking: if there are any operations in the CGO (I) that specifies
 * different values for picking. if has_picking is set, then a VBO for picking
 * is generated in CGOConvertToShader()
 *
 */
static void CheckAttributesForUsage(const CGO* I, AttribDataDesc& attrData,
                                    AttribDataDesc& pickData, CGO* cgo,
                                    bool& has_picking) {
  size_t attrIdx = 0;

  // need to check attributes:
  //  - remove any that are not needed
  //  - add glVertexAttrib for those that are removed
  std::map<int, int>
      opToAttrUsed;  // bitmask for each op to which attributes are set
  for (auto& attrDesc : attrData) {
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt != attrOps->end();
         ++attrOpIt) {
      auto attrOp = &(*attrOpIt);
      if (opToAttrUsed.find(attrOp->op) == opToAttrUsed.end())
        opToAttrUsed[attrOp->op] = 1 << attrDesc.order;
      else
        opToAttrUsed[attrOp->op] |= 1 << attrDesc.order;
    }
  }
  // add picking ops (1 << attrIdx) i.e., any pick op is the last bit
  int pidx = 0;
  for (auto pickDataIt = pickData.begin(); pickDataIt != pickData.end();
       ++pickDataIt) {
    auto pickDesc = &(*pickDataIt);
    auto pickOps = &pickDesc->attrOps;
    pickDesc->order = pidx++;
    for (auto pickOpIt = pickOps->begin(); pickOpIt != pickOps->end();
         ++pickOpIt) {
      auto pickOp = &(*pickOpIt);
      pickOp->desc = pickDesc;
      if (opToAttrUsed.find(pickOp->op) == opToAttrUsed.end())
        opToAttrUsed[pickOp->op] = 1 << attrIdx;
      else
        opToAttrUsed[pickOp->op] |= 1 << attrIdx;
    }
  }
  size_t totAttrIdx = (1 << (attrIdx + 1)) - 1;
  size_t allAttrIdxUsed = 0;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    if (opToAttrUsed.find(op) != opToAttrUsed.end()) {
      int attrUsed = opToAttrUsed[op];
      if (attrUsed & (1 << attrIdx)) {  // if picking, need to check values, and
                                        // take them out if cPickableNoPick
        switch (op) {
          case cgo::draw::shadercylinder2ndcolor::op_code:
            if (reinterpret_cast<const cgo::draw::shadercylinder2ndcolor*>(pc)
                    ->pick_color_bond == cPickableNoPick)
              attrUsed ^= (1 << attrIdx);
            break;
          case cgo::draw::splitline::op_code:
            if (reinterpret_cast<const cgo::draw::splitline*>(pc)->bond ==
                cPickableNoPick)
              attrUsed ^= (1 << attrIdx);
        }
      }
      allAttrIdxUsed |= attrUsed;
      if (allAttrIdxUsed == totAttrIdx) break;
    }
  }
  has_picking =
      allAttrIdxUsed & (1 << attrIdx);  // has_picking if the last bit is set

  if (allAttrIdxUsed != totAttrIdx) {
    // go through any attributes that aren't used:
    //   - add associated vertex_attribute type (if default_value is set)
    //   - remove attribute from attrData description so that it isn't included
    //   in VBO
    AttribDataDesc attrDataNew;
    for (auto idx = 0; idx < attrIdx; ++idx) {
      if (!attrData[idx].repeat_value && !(allAttrIdxUsed & (1 << idx))) {
        // attribute not used, need to create glVertexAttrib
        if (attrData[idx].default_value) {
          // need to add glVertexAttrib CGO OP
          int attr_lookup_idx =
              I->G->ShaderMgr->GetAttributeUID(attrData[idx].attr_name);
          switch (GetVertexFormatBaseType(attrData[idx].m_format)) {
            case VertexFormatBaseType::Float:
              switch (VertexFormatToGLSize(attrData[idx].m_format)) {
                case 1:
                  cgo->add<cgo::draw::vertex_attribute_1f>(
                      attr_lookup_idx, *(float*)attrData[idx].default_value);
                  break;
                case 3:
                  cgo->add<cgo::draw::vertex_attribute_3f>(
                      attr_lookup_idx, attrData[idx].default_value);
                  break;
                default:
                  std::cerr << "\tBAD SIZE: attrData[idx].m_format="
                            << static_cast<int>(attrData[idx].m_format)
                            << std::endl;
              }
              break;
            case VertexFormatBaseType::UByte:
              switch (VertexFormatToGLSize(attrData[idx].m_format)) {
                case 1: {
                  float val;
                  unsigned char valuc = *attrData[idx].default_value;
                  if (VertexFormatIsNormalized(attrData[idx].m_format)) {
                    val = CLAMPVALUE(valuc / 255.f, 0.f, 1.f);
                  } else {
                    val = (float)valuc;
                  }
                  cgo->add<cgo::draw::vertex_attribute_1f>(attr_lookup_idx,
                                                           val);
                } break;
                case 4:
                  cgo->add<cgo::draw::vertex_attribute_4ub>(
                      attr_lookup_idx, attrData[idx].default_value);
                  break;
                default:
                  std::cerr << "\tNOT IMPLEMENTED: attrData[idx].m_format="
                            << static_cast<int>(attrData[idx].m_format)
                            << std::endl;
              }
          }
        }
      } else {
        attrDataNew.push_back(attrData[idx]);
      }
    }
    attrData.swap(attrDataNew);  // only keep attributes that are used
  }
}

/**
 * Populates two structures and sets vertsperpickinfo
 * opToCntPer           : CGO op to how many vertices are generated for each op.
 * opToOrderedAttribOps : CGO op to an ordered map of AttribOps that define how
 * we operate on the attribute arrays for each CGO op. attrData             :
 * definition of attributes pickData             : definition of pick attributes
 * vertsperpickinfo     : number of vertices per pick (only used when
 *
 */
static void PopulateOpsIntoStructuresForConversion(
    std::map<int, int>& opToCntPer,
    std::map<int, std::map<int, AttribOp*>>& opToOrderedAttribOps,
    AttribDataDesc& attrData, AttribDataDesc& pickData, int& vertsperpickinfo,
    const bool has_picking) {
  size_t attrIdx = 0;
  for (auto& attrDesc : attrData) {
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt != attrOps->end();
         ++attrOpIt) {
      auto attrOp = &(*attrOpIt);
      attrOp->desc = &attrDesc;
      if (attrOp->copyFromAttr >= 0) {
        attrOp->copyAttribDesc = &attrData[attrOp->copyFromAttr];
      }
      if (attrOp->incr_vertices > 0) {
        if (!vertsperpickinfo) {
          vertsperpickinfo = attrOp->incr_vertices;
        } else {
          if (attrOp->incr_vertices != vertsperpickinfo) {
            std::cerr << "WARNING: attrOp->incr_vertices set to multiple "
                         "values, vertsperpickinfo="
                      << vertsperpickinfo
                      << " attrOp->incr_vertices=" << attrOp->incr_vertices
                      << " : picking might get confused" << std::endl;
          }
        }
        if (opToCntPer.find(attrOp->op) == opToCntPer.end())
          opToCntPer[attrOp->op] = attrOp->incr_vertices;
        else
          opToCntPer[attrOp->op] += attrOp->incr_vertices;
      }
      if (opToOrderedAttribOps.find(attrOp->op) == opToOrderedAttribOps.end())
        opToOrderedAttribOps[attrOp->op] = std::map<int, AttribOp*>({});
      opToOrderedAttribOps[attrOp->op][attrOp->order] = attrOp;
    }
  }
  if (has_picking) {
    for (auto pickDataIt = pickData.begin(); pickDataIt != pickData.end();
         ++pickDataIt) {
      auto pickDesc = &(*pickDataIt);
      auto pickOps = &pickDesc->attrOps;
      for (auto pickOpIt = pickOps->begin(); pickOpIt != pickOps->end();
           ++pickOpIt) {
        auto pickOp = &(*pickOpIt);
        if (pickOp->copyFromAttr >= 0) {
          pickOp->copyAttribDesc = &pickData[pickOp->copyFromAttr];
        }
        if (opToOrderedAttribOps.find(pickOp->op) == opToOrderedAttribOps.end())
          opToOrderedAttribOps[pickOp->op] = std::map<int, AttribOp*>({});
        opToOrderedAttribOps[pickOp->op][pickOp->order] = pickOp;
      }
    }
  }
}

/**
 * converts a "primitive" CGO into a CGO that renders a custom operation
 *
 * I:                   primitive CGO that is processed
 * attrData:            definition of attributes that are accumulated and put
 * into the VBO pickData:            definition of pick attributes that
 * accumulate pick data and put into the interleaved picking VBO mode: which
 * openGL mode to use for rendering (e.g., GL_POINTS, GL_TRIANGLES,
 * GL_LINE_STRIPS, etc.) layout:              SEPARATE, SEQUENTIAL, or
 * INTERLEAVED : how the VBO is layed out in memory (default: INTERLEAVED)
 * check_attr_for_data: if true, this function checks whether all attributes are
 * used, and if any are not specified, and default values are defined in
 * AttribDesc, then glVertexAttrib CGO ops are created on the returned CGO
 * (default: true) idx_array:           if specified, this array is used to
 * specify indices for each fragment in an Indexed Buffer and glDrawElements is
 * used to render (instead of glDrawArrays) nvertsperfrag:       in conjunction
 * with idx_array, how many vertices for each set of indices (default: 0) both
 * idx_array and nvertsperfrag are used to specify vertices and geometry for
 * each fragment, such as the box for cylinders nfragspergroup:      Currently,
 * this represents the number of fragments that are inside of a group.  For
 * example, crosses as cylinders (CGOConvertCrossesToCylinderShader) have 36
 * indices per fragment (nvertsperfrag=36) and 3 fragments per group
 * (nfragspergroup=3) (i.e., one for each line). note: idx_array, nvertsperfrag
 * and nfragspergroup should probably be moved to AttrOp in the future to
 * support different types of fragments.
 *
 * returns a CGO that consists of vertex_attrib_* (e.g., glVertexAttrib) and a
 * custom CGO operation that calls either glDrawArrays or glDrawElements. It
 * also supports picking if specified in the pickData.
 *
 */
CGO* CGOConvertToShader(const CGO* I, AttribDataDesc& attrData,
                        AttribDataDesc& pickData, int mode,
                        const VertexBufferLayout layout,
                        bool check_attr_for_data, int* idx_array,
                        int nvertsperfrag, int nfragspergroup) {
  CGO* cgo;
  int ok = true;
  bool isInterleaved = (layout == VertexBufferLayout::Interleaved);
  bool has_picking = true;
  std::map<std::string, AttribDesc*> attrToDesc;

  cgo = CGONew(I->G);
  cgo->use_shader = true;

  if (check_attr_for_data) {
    CheckAttributesForUsage(I, attrData, pickData, cgo, has_picking);
  } else {
    // if attributes aren't checked, still need to set pick order and desc
    // pointer
    int pidx = 0;
    for (auto& pickDesc : pickData) {
      pickDesc.order = pidx++;
      for (auto& pickOp : pickDesc.attrOps) {
        pickOp.desc = &pickDesc;
      }
    }
  }
  std::map<int, int> opToCntPer;
  std::map<int, std::map<int, AttribOp*>> opToOrderedAttribOps;

  int vertsperpickinfo = 0;
  PopulateOpsIntoStructuresForConversion(opToCntPer, opToOrderedAttribOps,
                                         attrData, pickData, vertsperpickinfo,
                                         has_picking);

  // Populate these variables used for accumulating and setting VBO data arrays
  int vertexDataSize = 0;
  std::vector<int> attrSizes;
  std::vector<int> attrOffset;  // for interleaved
  int curoffset = 0;            // for interleaved
  size_t attrIdx = 0;
  for (auto& attrDesc : attrData) {
    attrDesc.order = attrIdx++;
    auto attrSize = GetSizeOfVertexFormat(attrDesc.m_format);
    attrSizes.push_back(attrSize);
    vertexDataSize += attrSize;

    if (isInterleaved) {
      attrOffset.push_back(curoffset);
      curoffset += attrSize;
    } else {
      attrOffset.push_back(0);  // no offset when not interleaved
    }
    attrToDesc[attrDesc.attr_name] = &attrDesc;
  }

  // populate funcData.attrib
  for (auto& attrDesc : attrData)
    for (auto& attrop : attrDesc.attrOps)
      for (auto& funcData : attrop.funcDataConversions) {
        if (attrToDesc.find(funcData.attribName) != attrToDesc.end()) {
          funcData.attrib = attrToDesc[funcData.attribName];
        }
      }

  // Since some cards require word-aligned strides (e.g., ATI)
  // we need to make sure it is word-aligned.  If it isn't, and you
  // want to save memory, maybe SEQUENTIAL is a better option
  if (vertexDataSize % 4) {
    vertexDataSize += 4 - (vertexDataSize % 4);
  }

  int ntotalverts = CGOCountNumberOfOperationsOfTypeN(I, opToCntPer);

  // PYMOL-2668
  if (ntotalverts == 0) {
    CGOStop(cgo);
    return cgo;
  }

  size_t pickvbohash = 0;
  int pickDataSize = pickData.size();
  int pstride = 2 * pickDataSize;
  // Generate VBOs: both for vertex data (vbo) and picking color data (pickvbo)
  // if necessary
  // - pickvbo is interleaved so that for multiple channels, pick data for each
  // vertex
  //   is contiguous
  auto& vbo = CreateVertexBuffer(I->G, layout);
  VertexBuffer* pickvbo = nullptr;
  if (pickDataSize) {
    pickvbo = &CreateVertexBuffer(I->G, VertexBufferLayout::Sequential,
                                  MemoryUsageProperty::CpuToGpu);
    pickvbohash = pickvbo->get_hash_id();
  }

  // adding SPECIAL OPERATIONS (for now) before adding custom OP
  // - for now, this is the only operation that needs to be passed to the new
  // CGO
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    switch (op) {
      case CGO_SPECIAL:
        cgo->add_to_cgo(op, pc);
        break;
    }
  }

  // defines how many passes we have
  // We *could* set npickcolattr=1 if we had access to PickColorConverter here
  // and getTotalBits() == 32
  const int npickcolattr = SHADER_PICKING_PASSES_MAX;
  /* for picking, we need to mask the attributes and bind the buffer in the the
   * picking VBO */
  int pl = 0;
  int npickattr = pickData.size();
  for (auto& pickDesc : pickData) {
    cgo->add<cgo::draw::mask_attribute_if_picking>(
        I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name),
        vbo.get_hash_id());
    if (has_picking) {
      cgo->add<cgo::draw::bind_vbo_for_picking>(pickvbohash, pl++, npickattr);
    } else {
      /* if no picking, should render black */
      static unsigned char zerocolor[]{0, 0, 0, 0};
      cgo->add<cgo::draw::vertex_attribute_4ub_if_picking>(
          I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name), zerocolor);
    }
  }
  size_t iboid = 0;
  int num_total_indexes = 0;
  if (nvertsperfrag) {
    int nfrags = nfragspergroup * ntotalverts / vertsperpickinfo;
    int nvertsperindivfrag = vertsperpickinfo / nfragspergroup;
    num_total_indexes = nfrags * nvertsperfrag;
    std::vector<VertexIndex_t> vertexIndices(num_total_indexes);
    int idxpl = 0;
    // using vertsperpickinfo as verts per frag
    for (int cnt = 0, vpl = 0; cnt < nfrags; ++cnt) {
      for (int idx_array_pl = 0; idx_array_pl < nvertsperfrag; ++idx_array_pl) {
        vertexIndices[idxpl] = idx_array[idx_array_pl] + vpl;
        idxpl++;
      }
      vpl += nvertsperindivfrag;
    }
    auto& ibo = CreateIndexBuffer(I->G);
    SetIndexBufferData(I->G, ibo, vertexIndices);
    iboid = ibo.get_hash_id();
  }

  // pick_data is interleaved if more than one attribute
  float* pick_data = cgo->add<cgo::draw::custom>(
      mode, ntotalverts, vbo.get_hash_id(), pickvbohash, vertsperpickinfo,
      pickDataSize, iboid, num_total_indexes);
  std::vector<unsigned char> allData(ntotalverts * vertexDataSize);
  std::vector<void*> dataPtrs;
  std::vector<int> repeat_attr_idx;
  int allAttrBS = 0;

  // Initialize first entry in array(s) with default values and populate
  // dataPtrs
  if (isInterleaved) {
    int pl = 0;
    auto attrDataIt = attrData.begin();
    auto attrOffsetIt = attrOffset.begin();
    for (; attrDataIt != attrData.end() && attrOffsetIt != attrOffset.end();
         ++attrDataIt, ++attrOffsetIt) {
      auto attrDesc = &(*attrDataIt);
      if (attrDesc->repeat_value) {
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc->order);
      }
      auto attrOffset = *attrOffsetIt;
      unsigned char* first_value = nullptr;
      first_value =
          (attrDesc->default_value
               ? attrDesc->default_value
               : (attrDesc->repeat_value ? attrDesc->repeat_value : nullptr));
      if (first_value) {
        auto attrSize = GetSizeOfVertexFormat(attrDesc->m_format);
        memcpy(allData.data() + attrOffset, first_value, attrSize);
      }
      dataPtrs.push_back(allData.data());
      ++pl;
    }
  } else {
    auto* curAllDataPtr = allData.data();
    int pl = 0;
    for (auto& attrDesc : attrData) {
      if (attrDesc.repeat_value) {
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc.order);
      }
      dataPtrs.push_back(curAllDataPtr);
      unsigned char* first_value = nullptr;
      first_value =
          (attrDesc.default_value
               ? attrDesc.default_value
               : (attrDesc.repeat_value ? attrDesc.repeat_value : nullptr));
      if (first_value) {
        memcpy(curAllDataPtr, first_value, attrSizes[pl]);
      }
      curAllDataPtr = curAllDataPtr + ntotalverts * attrSizes[pl];
      ++pl;
    }
  }

  int nvert = 0;
  int attrBS = 0;
  int allPickAttrBS = (1 << pickData.size()) - 1;
  int has_pick_colorBS = allPickAttrBS;
  for (int pi = 0; pi < pickDataSize; ++pi) {
    CGO_put_uint(2 * pi + pick_data, 0);
    CGO_put_int(2 * pi + pick_data + 1, cPickableNoPick);
  }
  bool cont =
      true;  // need to break while statement as well if past all vertices

  // This is the loop that goes through the CGO and accumulates all of the data
  // for both the rendering and picking VBOs.
  // - For each OP, go through the list of Attribute OPs
  // - Each Attribute OP:
  //    - copies attribute data from the OP
  //    - can generate vertices (incr_vertices)
  // - Attributes are kept track of (attrBS) and when vertices are generated,
  //   the attributes that have not been newly written for the current vertex
  //   are copied from the previous one.
  for (auto it = I->begin(); cont && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    if (opToOrderedAttribOps.find(op) != opToOrderedAttribOps.end()) {
      std::map<int, AttribOp*>* attribOpsInOrder = &opToOrderedAttribOps[op];
      for (auto attribOpIt : *attribOpsInOrder) {
        AttribOp* attribOp = attribOpIt.second;
        int ord = attribOp->desc->order;
        cont = nvert < ntotalverts;
        if (!cont) break;
        copyAttributeForOp(isInterleaved, nvert, attribOp, vertexDataSize,
                           dataPtrs, attrOffset, pc, pick_data,
                           has_pick_colorBS, pstride);
        if (ord >= 0)  // picking is negative, has_pick_colorBS is used instead
          attrBS |= (1 << ord);
        else
          std::cout << "   ord=%d\n" << ord << std::endl;
        if (attribOp->incr_vertices) {
          if (has_pick_colorBS != allPickAttrBS) {
            // copy pick colors that haven't been set from previous vertex
            for (int pi = 0; pi < pickDataSize; ++pi) {
              if (has_pick_colorBS ^ (1 << pi)) {
                CGO_put_uint(pick_data, CGO_get_uint(pick_data - pstride));
                CGO_put_int(pick_data + 1,
                            CGO_get_int(pick_data - pstride + 1));
              }
              pick_data += 2;
            }
          } else {
            pick_data += pstride;
          }
          has_pick_colorBS = 0;
          if (!nvert && attrBS != allAttrBS) {
            // for the first vertex, all attributes should be set
            for (auto idx = 0; idx < attrData.size(); ++idx) {
              if (!(attrBS & (1 << idx))) {
                if (!attrData[idx].default_value) {
                  std::cerr << "WARNING: attribute #" << idx << " ("
                            << attrData[idx].attr_name
                            << ") not set for first"
                               " vertex and does not have default value"
                            << std::endl;
                }
              }
            }
          }
          if (nvert && attrBS != allAttrBS) {
            // for each vertex that hasn't been written for the current vertex,
            // copy it from the previous vertex
            for (auto idx = 0; idx < attrData.size(); ++idx) {
              if (!(attrBS & (1 << idx))) {
                copyAttributeForVertex(isInterleaved, nvert, attrData[idx],
                                       vertexDataSize, dataPtrs, attrOffset);
              }
            }
          }
          attrBS = 0;

          // creating new vertices, all attribute data should be copied into new
          // vertex.
          for (int nxt = 0; nxt < attribOp->incr_vertices; ++nxt) {
            /* for now, always copy values from previous */
            ++nvert;
            if (nvert < ntotalverts) {
              {
                // last should not need to copy into next, since we call
                // copyAttributeForVertex (above) for all attributes that
                // haven't been set
                // - note: it might be faster (especially for interleaved) to
                // copy into the next
                //         vertex, then the above copyAttributeForVertex() would
                //         not be needed
                if (isInterleaved) {
                  auto* dest = allData.data() + vertexDataSize * nvert;
                  memcpy(dest, dest - vertexDataSize, vertexDataSize);
                } else {
                  auto dataPtrIt = dataPtrs.begin();
                  auto attrDataIt = attrData.begin();
                  for (; attrDataIt != attrData.end() &&
                         dataPtrIt != dataPtrs.end();
                       ++attrDataIt, ++dataPtrIt) {
                    auto attrDesc = &(*attrDataIt);
                    auto dataPtr = *dataPtrIt;
                    auto attrSize = GetSizeOfVertexFormat(attrDesc->m_format);
                    void* dest = ((unsigned char*)dataPtr) + attrSize * nvert;
                    memcpy((unsigned char*)dest,
                           ((unsigned char*)dest) - attrSize, attrSize);
                  }
                }
              }
              // always copy repeat attributes
              if (!repeat_attr_idx.empty()) {
                for (auto ridx = repeat_attr_idx.begin();
                     ridx != repeat_attr_idx.end(); ++ridx) {
                  copyAttributeForVertex(isInterleaved, nvert, attrData[*ridx],
                                         vertexDataSize, dataPtrs, attrOffset);
                }
              }
            }
            if (!attribOp->funcDataConversions.empty()) {
              // for all attributes, call the funcDataConversion() if defined
              int nvert_m_1 = nvert - 1;
              for (auto funcData : attribOp->funcDataConversions) {
                auto funcAttrib = funcData.attrib;
                auto order = funcAttrib->order;
                if (isInterleaved) {
                  unsigned char* dest = allData.data() +
                                        vertexDataSize * nvert_m_1 +
                                        attrOffset[order];
                  funcData.funcDataConversion(dest, pc,
                                              funcData.funcDataGlobalArg, nxt);
                } else {
                  auto dataPtr = dataPtrs[order];
                  void* dest =
                      ((unsigned char*)dataPtr) + attrSizes[order] * nvert_m_1;
                  funcData.funcDataConversion(dest, pc,
                                              funcData.funcDataGlobalArg, nxt);
                }
              }
            }
          }
        }
      }
    }
  }

  /* Generate Pick Buffers with all pick attributes (if necessary) */
  if (pickvbo) {
    BufferDataDesc pickBufferData;
    for (int i = 0; i < npickcolattr; i++) {
      for (auto& pickDesc : pickData) {
        auto pickSize = GetSizeOfVertexFormat(pickDesc.m_format);
        pickBufferData.descs.push_back(BufferDesc(
            pickDesc.attr_name, pickDesc.m_format, pickSize * nvert));
      }
    }
    SetVertexBufferData(I->G, *pickvbo, std::move(pickBufferData));
  }

  /* Generate VBO Buffers with all pick attributes based on the VertexBuffer
   * type SEPARATE/SEQUENTIAL/INTERLEAVED*/
  BufferDataDesc bufferData;
  switch (layout) {
    case VertexBufferLayout::Separate:
    case VertexBufferLayout::Sequential: {
      auto attrDataIt = attrData.begin();
      auto dataPtrIt = dataPtrs.begin();
      auto attrSizeIt = attrSizes.begin();
      for (; attrDataIt != attrData.end() && dataPtrIt != dataPtrs.end() &&
             attrSizeIt != attrSizes.end();
           ++attrDataIt, ++dataPtrIt, ++attrSizeIt) {
        auto attrDesc = &(*attrDataIt);
        auto dataPtr = *dataPtrIt;
        auto attrSize = *attrSizeIt;
        bufferData.descs.push_back(BufferDesc(attrDesc->attr_name,
                                              attrDesc->m_format,
                                              nvert * attrSize, dataPtr));
      }
      SetVertexBufferData(I->G, vbo, std::move(bufferData));
      break;
    } break;
    case VertexBufferLayout::Interleaved: {
      auto attrDataIt = attrData.begin();
      auto attrOffsetIt = attrOffset.begin();
      for (; attrDataIt != attrData.end() && attrOffsetIt != attrOffset.end();
           ++attrDataIt, ++attrOffsetIt) {
        auto attrDesc = &(*attrDataIt);
        auto offset = *attrOffsetIt;
        bufferData.descs.push_back(BufferDesc{attrDesc->attr_name,
                                              attrDesc->m_format, 0, nullptr,
                                              (std::uint32_t)offset});
      }
      bufferData.stride = static_cast<std::size_t>(vertexDataSize);
      SetVertexBufferDataSized(I->G, vbo, std::move(bufferData),
                               nvert * vertexDataSize,
                               reinterpret_cast<std::byte*>(allData.data()));
      break;
    }
  }

  CGOStop(cgo);
  return cgo;
}

// CGOCheckSplitLineInterpolationIsSame:
//   - returns true if always the same
//   - returns false if not always the same
bool CGOCheckSplitLineInterpolationIsSame(const CGO* I, bool& interp_value) {
  bool interp_value_first = false;
  bool interp_value_is_set = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
      case cgo::draw::splitline::op_code:
        interp_value = (it.cast<cgo::draw::splitline>()->flags &
                        cgo::draw::splitline::interpolation);
        break;
      case CGO_INTERPOLATED:
        interp_value = it.cast<float>()[0] > 0.5f;
        break;
      default:
        continue;
    }
    if (!interp_value_is_set) {
      interp_value_first = interp_value;
      interp_value_is_set = true;
    } else if (interp_value != interp_value_first) {
      return false;
    }
  }
  return true;
}

// CGOCheckShaderCylinderCapInfoIsSame:
//   - returns true if always the same
//   - returns false if not always the same
static bool CGOCheckShaderCylinderCapInfoIsSame(const CGO* I,
                                                unsigned char& cap_value) {
  unsigned char cap_value_first = 0;
  bool cap_value_first_is_set = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
      case cgo::draw::shadercylinder::op_code:
        cap_value = it.cast<cgo::draw::shadercylinder>()->cap;
        break;
      case cgo::draw::shadercylinder2ndcolor::op_code:
        cap_value = it.cast<cgo::draw::shadercylinder2ndcolor>()->cap;
        break;
      case cgo::draw::sausage::op_code:
        cap_value = cCylShaderBothCapsRound | cCylShaderInterpColor;
        break;
      case cgo::draw::cylinder::op_code:
        cap_value = cCylShaderBothCapsFlat | cCylShaderInterpColor;
        break;
      case cgo::draw::custom_cylinder::op_code: {
        auto cc = it.cast<cgo::draw::custom_cylinder>();
        const cCylCap cap1 = cc->get_cap1();
        const cCylCap cap2 = cc->get_cap2();
        cap_value =
            cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
      } break;
      case cgo::draw::custom_cylinder_alpha::op_code: {
        auto cc = it.cast<cgo::draw::custom_cylinder_alpha>();
        const cCylCap cap1 = cc->get_cap1();
        const cCylCap cap2 = cc->get_cap2();
        cap_value =
            cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
      } break;
      default:
        continue;
    }

    if (!cap_value_first_is_set) {
      cap_value_first = cap_value;
      cap_value_first_is_set = true;
    } else if (cap_value != cap_value_first) {
      return false;
    }
  }

  return true;
}

CGO* CGOConvertToTrilinesShader(const CGO* I, CGO* addTo, bool add_color) {
  PyMOLGlobals* G = I->G;

  AttribDataOp vertexOps = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 0},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex1), 0}};
  AttribDataOp vertexOtherOps = {
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 6},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex2), 6}};
  AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 6, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp color2Ops = {
      {CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 7, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::splitline, index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::splitline, index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  AttribDataDesc attrDesc = {
      {"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_OtherVertex", VertexFormat::Float3, vertexOtherOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"a_UV", VertexFormat::UByte}};

  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
    attrDesc[2].default_value = default_color;
    attrDesc[3].default_value = default_color;
  }
  AttribDesc* uvdesc = &attrDesc[attrDesc.size() - 1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = {1, 3, 0, 3, 2, 1};
  uvdesc->repeat_value = uv_bits;

  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
                               offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 2);  // a_Color
    attrDesc.erase(attrDesc.begin() + 2);  // a_Color2
  }

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved);
}

CGO* CGOConvertToLinesShader(const CGO* I, CGO* addTo, bool add_color) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertexOps = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 1},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex1), 1},
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 1},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex2), 1}};
  AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::splitline, index), 0}};

  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
                             {"a_Color", VertexFormat::UByte4Norm, colorOps}};
  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
    attrDesc[1].default_value = default_color;
  }
  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
                               offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 1);  // a_Color
  }

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_LINES,
                            VertexBufferLayout::Interleaved);
}

CGO* CGOConvertLinesToCylinderShader(const CGO* I, CGO* addTo, bool add_color) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  AttribDataOp vertex1Ops = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 0},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex1), 0}};
  AttribDataOp vertex2Ops = {
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 8},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::splitline, vertex2), 8}};
  static AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 6, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  static AttribDataOp color2Ops = {
      {CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};

  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
                             {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
                             {"a_Color", VertexFormat::UByte4Norm, colorOps},
                             {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
                             {"attr_radius", VertexFormat::Float}};

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {0, 4, 6, 2,
                                      1, 5, 7, 3};  // right(4)/up(2)/out(1)

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
    fdesc = &attrDesc[2];
    fdesc->default_value = default_color;
    fdesc = &attrDesc[3];
    fdesc->default_value = default_color;
  }

  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = {// box indices
                         0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0,
                         0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_cap"),
        (cCylShaderBothCapsRound | (interp_value ? cCylShaderInterpColor : 0)));
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_INTERP_TO_CAP,
                               offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_cap attribute
    attrDesc.push_back({"a_cap", VertexFormat::UByte, interpOps});
  }

  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 2);  // attr_colors
    attrDesc.erase(attrDesc.begin() + 2);  // attr_colors2
  }

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 7, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::splitline, index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::splitline, index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved, true,
                            box_indices_ptr, 36);
}

struct CrossSizeData {
  float cross_size;
  bool forward;
  CrossSizeData(float _cross_size, bool _forward)
      : cross_size(_cross_size), forward(_forward) {}
};

static void CrossVertexConversion(void* varData, const float* pc,
                                  void* crossData, int idx) {
  CrossSizeData* csd = (CrossSizeData*)crossData;
  int idxpl = idx / 8;  // X Y or Z
  float* varDataF = ((float*)varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToCylinderShader(const CGO* I, CGO* addTo,
                                       float cross_size_arg) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertex1Ops = {{CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 0}};
  AttribDataOp vertex2Ops = {
      {CGO_VERTEX_CROSS, 2, FLOAT3_TO_FLOAT3, 0, 3 * 8, 0}};

  static AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
                                  {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  static AttribDataOp color2Ops = {{CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
                                   {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0}};

  CrossSizeData crossData[] = {{cross_size_arg, false}, {cross_size_arg, true}};
  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
                             {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
                             {"a_Color", VertexFormat::UByte4Norm, colorOps},
                             {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
                             {"attr_radius", VertexFormat::Float}};

  attrDesc.reserve(10);
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversion, &crossData[0], "attr_vertex1"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversion, &crossData[1], "attr_vertex2"});

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {0, 4, 6, 2,
                                      1, 5, 7, 3};  // right(4)/up(2)/out(1)

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = {
      255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = {// box indices
                         0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0,
                         0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  addTo->add<cgo::draw::vertex_attribute_1f>(
      G->ShaderMgr->GetAttributeUID("a_cap"), cCylShaderBothCapsRound);

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved, true,
                            box_indices_ptr, 36, 3);
}

struct CrossSizeDataLines {
  float cross_size;
  CrossSizeDataLines(float _cross_size) : cross_size(_cross_size) {}
};

static void CrossVertexConversionLines(void* varData, const float* pc,
                                       void* crossData, int idx) {
  CrossSizeDataLines* csd = (CrossSizeDataLines*)crossData;
  int idxpl = idx / 2;  // X Y or Z
  bool forward = idx % 2;
  float* varDataF = ((float*)varData);
  varDataF[idxpl] += (forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToLinesShader(const CGO* I, CGO* addTo,
                                    float cross_size_arg) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertexOps = {
      {CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 6}};  // 6 vertices for a cross
  AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
                           {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};

  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
                             {"a_Color", VertexFormat::UByte4Norm, colorOps}};
  unsigned char default_color[] = {
      255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
  attrDesc[1].default_value = default_color;

  CrossSizeDataLines crossData = {cross_size_arg};
  attrDesc[0].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionLines, &crossData, "a_Vertex"});

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
                               offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_LINES,
                            VertexBufferLayout::Interleaved);
}

static void CrossVertexConversionTrilines(void* varData, const float* pc,
                                          void* crossData, int idx) {
  CrossSizeData* csd = (CrossSizeData*)crossData;
  int idxpl = idx / 6;  // X Y or Z
  float* varDataF = ((float*)varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToTrilinesShader(const CGO* I, CGO* addTo,
                                       float cross_size_arg) {
  PyMOLGlobals* G = I->G;

  AttribDataOp vertexOps = {{CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 0}};
  AttribDataOp vertexOtherOps = {
      {CGO_VERTEX_CROSS, 2, FLOAT3_TO_FLOAT3, 0, 6 * 3}};
  AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
                           {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp color2Ops = {{CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
                            {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  AttribDataDesc attrDesc = {
      {"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_OtherVertex", VertexFormat::Float3, vertexOtherOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"a_UV", VertexFormat::UByte}};

  CrossSizeData crossData[] = {{cross_size_arg, false}, {cross_size_arg, true}};

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionTrilines, &crossData[0], "a_Vertex"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionTrilines, &crossData[1], "a_OtherVertex"});

  unsigned char default_color[] = {
      255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
  attrDesc[2].default_value = default_color;
  attrDesc[3].default_value = default_color;

  AttribDesc* uvdesc = &attrDesc[attrDesc.size() - 1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = {1, 3, 0, 3, 2, 1};
  uvdesc->repeat_value = uv_bits;

  addTo->add<cgo::draw::vertex_attribute_1f>(
      G->ShaderMgr->GetAttributeUID("a_interpolate"), 0.f);

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved);
}

static void SetVertexFromOriginAxisForCylinder(void* varData, const float* pc,
                                               void* np, int idx) {
  float* varDataF = ((float*)varData);
  add3f(pc, pc + 3, varDataF);  // adding origin and axis for both
                                // shadercylinder and shadercylinder2ndcolor
}

/**
 * converts all cylinders in the input CGO to a CGO custom operation, which
 * includes picking information (if it exists)
 *
 * I     - input CGO (includes cylinders)
 * addTo - CGO that vertex_attribute operations are added to (if needed), in
 * this case for caps if values for all cylinders are the same
 *
 */
CGO* CGOConvertShaderCylindersToCylinderShader(const CGO* I, CGO* addTo) {
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  // TODO: NEED TO ADD: CGO_CUSTOM_CYLINDER and CGO_CYLINDER
  AttribDataOp vertex1Ops = {
      {CGO_SHADER_CYLINDER, 1, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::shadercylinder, origin), 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 1, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::shadercylinder2ndcolor, origin), 0},
      {CGO_SAUSAGE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::sausage, vertex1),
       0},
      {CGO_CYLINDER, 1, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::cylinder, vertex1), 0},
      {CGO_CUSTOM_CYLINDER, 1, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::custom_cylinder, vertex1), 0},
      {CGO_CUSTOM_CYLINDER_ALPHA, 1, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::custom_cylinder_alpha, vertex1), 0}};
  AttribDataOp vertex2Ops = {
      {CGO_SHADER_CYLINDER, 5, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::shadercylinder, axis), 8},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 6, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::shadercylinder2ndcolor, axis), 8},
      {CGO_SAUSAGE, 6, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::sausage, vertex2),
       8},
      {CGO_CYLINDER, 6, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::cylinder, vertex2), 8},
      {CGO_CUSTOM_CYLINDER, 6, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::custom_cylinder, vertex2), 8},
      {CGO_CUSTOM_CYLINDER_ALPHA, 6, FLOAT3_TO_FLOAT3,
       offsetof(cgo::draw::custom_cylinder_alpha, vertex2), 8}};
  static AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SAUSAGE, 4, FLOAT3_TO_UB3, offsetof(cgo::draw::sausage, color1)},
      {CGO_CYLINDER, 4, FLOAT3_TO_UB3, offsetof(cgo::draw::cylinder, color1)},
      {CGO_CUSTOM_CYLINDER, 4, FLOAT3_TO_UB3,
       offsetof(cgo::draw::custom_cylinder, color1)},
      {CGO_CUSTOM_CYLINDER_ALPHA, 4, FLOAT4_TO_UB4,
       offsetof(cgo::draw::custom_cylinder_alpha, color1)}};
  static AttribDataOp color2Ops = {
      {CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 2, FLOAT3_TO_UB3,
       offsetof(cgo::draw::shadercylinder2ndcolor, color2)},
      {CGO_SAUSAGE, 5, FLOAT3_TO_UB3, offsetof(cgo::draw::sausage, color2)},
      {CGO_CYLINDER, 5, FLOAT3_TO_UB3, offsetof(cgo::draw::cylinder, color2)},
      {CGO_CUSTOM_CYLINDER, 5, FLOAT3_TO_UB3,
       offsetof(cgo::draw::custom_cylinder, color2)},
      {CGO_CUSTOM_CYLINDER_ALPHA, 5, FLOAT4_TO_UB4,
       offsetof(cgo::draw::custom_cylinder_alpha, color2)}};
  AttribDataOp radiusOps = {
      {CGO_SHADER_CYLINDER, 2, FLOAT_TO_FLOAT,
       offsetof(cgo::draw::shadercylinder, tube_size), 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 3, FLOAT_TO_FLOAT,
       offsetof(cgo::draw::shadercylinder2ndcolor, tube_size), 0},
      {CGO_SAUSAGE, 3, FLOAT_TO_FLOAT, offsetof(cgo::draw::sausage, radius), 0},
      {CGO_CYLINDER, 3, FLOAT_TO_FLOAT, offsetof(cgo::draw::cylinder, radius),
       0},
      {CGO_CUSTOM_CYLINDER, 3, FLOAT_TO_FLOAT,
       offsetof(cgo::draw::custom_cylinder, radius), 0},
      {CGO_CUSTOM_CYLINDER_ALPHA, 3, FLOAT_TO_FLOAT,
       offsetof(cgo::draw::custom_cylinder_alpha, radius), 0}};

  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
                             {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
                             {"a_Color", VertexFormat::UByte4Norm, colorOps},
                             {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
                             {"attr_radius", VertexFormat::Float, radiusOps}};

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {0, 4, 6, 2,
                                      1, 5, 7, 3};  // right(4)/up(2)/out(1)

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {SetVertexFromOriginAxisForCylinder, nullptr, "attr_vertex2"});
  attrDesc[1].attrOps[1].funcDataConversions.push_back(
      {SetVertexFromOriginAxisForCylinder, nullptr, "attr_vertex2"});

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = {
      255, 255, 255, 255};  // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = {// box indices
                         0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0,
                         0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  bool interp_same;
  unsigned char cap_value = 0;
  if ((interp_same = CGOCheckShaderCylinderCapInfoIsSame(I, cap_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_cap"), cap_value);
  } else {
    AttribDataOp interpOps = {
        {CGO_SHADER_CYLINDER, 3, CYL_CAP_TO_CAP,
         offsetof(cgo::draw::shadercylinder, cap), 0},
        {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 4, CYL_CAP_TO_CAP,
         offsetof(cgo::draw::shadercylinder2ndcolor, cap), 0},
        {CGO_SAUSAGE, 2, CYL_CAPS_ARE_ROUND, 0, 0},
        {CGO_CYLINDER, 2, CYL_CAPS_ARE_FLAT, 0, 0},
        {CGO_CUSTOM_CYLINDER, 2, CYL_CAPS_ARE_CUSTOM,
         offsetof(cgo::draw::custom_cylinder, cap1), 0},
        {CGO_CUSTOM_CYLINDER_ALPHA, 2, CYL_CAPS_ARE_CUSTOM,
         offsetof(cgo::draw::custom_cylinder_alpha, cap1), 0},
    };
    attrDesc.push_back({"a_cap", VertexFormat::UByte, interpOps});
  }

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 8, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 5, UINT_INT_TO_PICK_DATA,
       offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
                            VertexBufferLayout::Interleaved, true,
                            box_indices_ptr, 36);
}
