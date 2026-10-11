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

// Packing CGOs into vertex and index buffers for the GPU.

#include "pymol/graphics/cgo.h"
#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/scene/scene.h"
#include "pymol/graphics/graphics_util.h"
#include "pymol/graphics/cgo_internal.h"

#define VAR_FOR_NORMAL pl
#define VAR_FOR_NORMAL_CNT_PLUS

#define VALUES_PER_IMPOSTER_SPACE_COORD 1
constexpr unsigned VERTEX_ACCESSIBILITY_SIZE = 1;

static void CGOAddVertexToDrawArrays(CGO* cgo, int pl, int plc, int pla,
                                     const float* vertex, short notHaveValue,
                                     float* vertexVals, float* normalVals,
                                     float* colorVals, float* pickColorVals,
                                     float* accessibilityVals) {
  float* tmp_ptr;
  if (notHaveValue & CGO_NORMAL_ARRAY) {
    if (pl) {
      tmp_ptr = &normalVals[pl - 3];
      copy3f(tmp_ptr, &normalVals[pl]);
    } else {
      copy3f(cgo->normal, &normalVals[pl]);
    }
  }
  if (notHaveValue & CGO_COLOR_ARRAY) {
    if (plc) {
      tmp_ptr = &colorVals[plc - 4];
      copy4f(tmp_ptr, &colorVals[plc]);
    } else {
      copy3f(&colorVals[plc], cgo->color);
      colorVals[plc + 3] = cgo->alpha;
    }
  }
  if (pickColorVals) {
    CGO_put_uint(pickColorVals + pla * 2, cgo->current_pick_color_index);
    CGO_put_int(pickColorVals + pla * 2 + 1, cgo->current_pick_color_bond);
  }
  if (accessibilityVals) {
    accessibilityVals[pla] = cgo->current_accessibility;
  }
  copy3f(vertex, &vertexVals[pl]);
}

bool CGOCombineBeginEnd(CGO** I, bool do_not_split_lines) {
  CGO* cgo = CGOCombineBeginEnd(*I, 0, do_not_split_lines);
  CGOFree(*I);
  *I = cgo;
  return (cgo != nullptr);
}

/**
 * Converts Begin/End blocks into CGO_DRAW_ARRAYS operations.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param do_not_split_lines ???
 *
 * @return New converted CGO with has_begin_end=false
 */
CGO* CGOCombineBeginEnd(const CGO* I, int est, bool do_not_split_lines) {
  CGO* cgo;

  int ok = true;
  if (!I) return nullptr;
  cgo = CGONewSized(I->G, 0);
  ok &= cgo ? true : false;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
      case CGO_END:
      case CGO_VERTEX:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOCombineBeginEnd: op=0x%02x encountered without CGO_BEGIN\n",
            op ENDFB(I->G);
        break;
      case CGO_BEGIN: {
        float firstColor[3], firstAlpha;
        char hasFirstColor = 0, hasFirstAlpha = 0;
        int nverts = 0, damode = CGO_VERTEX_ARRAY, err = 0;

        // read int argument of the BEGIN operation
        int mode = CGO_get_int(it.data());
        ++it;

        // we want to iterate twice over the BEGIN/END block
        auto it2 = it;

        // first iteration over BEGIN/END block (consumes 'it')
        for (; !err && it != CGO_END; ++it) {
          auto pc = it.data();
          switch (it.op_code()) {
            case CGO_DRAW_ARRAYS:
            case CGO_STOP:
              PRINTFB(I->G, FB_CGO, FB_Errors)
              " CGO-Error: CGOCombineBeginEnd: invalid op=0x%02x inside "
              "BEGIN/END\n",
                  it.op_code() ENDFB(I->G);
              err = true;
              continue;
            case CGO_NORMAL:
              damode |= CGO_NORMAL_ARRAY;
              break;
            case CGO_COLOR:
              if (!nverts) {
                hasFirstColor = 1;
                copy3f(pc, firstColor);
              } else {
                hasFirstColor = 0;
                damode |= CGO_COLOR_ARRAY;
              }
              copy3f(pc, cgo->color);
              break;
            case CGO_PICK_COLOR:
              damode |= CGO_PICK_COLOR_ARRAY;
              break;
            case CGO_ACCESSIBILITY:
              damode |= CGO_ACCESSIBILITY_ARRAY;
              break;
            case CGO_VERTEX:
              nverts++;
              break;
            case CGO_ALPHA:
              cgo->alpha = *pc;
              if (!nverts) {
                hasFirstAlpha = 1;
                firstAlpha = cgo->alpha;
              } else {
                hasFirstAlpha = 0;
                damode |= CGO_COLOR_ARRAY;
              }
              break;
            case CGO_LINE:
              nverts += 2;
              break;
            case CGO_SPLITLINE: {
              auto splitline =
                  reinterpret_cast<const cgo::draw::splitline*>(pc);
              if (do_not_split_lines ||
                  (splitline->flags & cgo::draw::splitline::equal_colors)) {
                nverts += 2;
              } else {
                nverts += 4;
              }
            } break;
          }
        }
        if (nverts > 0 && !err) {
          int pl = 0, plc = 0, pla = 0;
          float *vertexVals = nullptr, *normalVals = nullptr,
                *colorVals = nullptr, *pickColorVals = nullptr,
                *accessibilityVals = nullptr;

          if (hasFirstAlpha || hasFirstColor) {
            if (hasFirstAlpha) {
              CGOAlpha(cgo, firstAlpha);
            }
            if (hasFirstColor) {
              CGOColorv(cgo, firstColor);
            }
          }

          float* nxtVals = cgo->add<cgo::draw::arrays>(mode, damode, nverts);
          ok &= nxtVals ? true : false;
          if (!ok) continue;

          assert(damode & CGO_VERTEX_ARRAY);
          vertexVals = nxtVals;
          nxtVals += VERTEX_POS_SIZE * nverts;

          if (damode & CGO_NORMAL_ARRAY) {
            normalVals = nxtVals;
            nxtVals += nverts * VERTEX_NORMAL_SIZE;
          }

          if (damode & CGO_COLOR_ARRAY) {
            colorVals = nxtVals;
            nxtVals += nverts * VERTEX_COLOR_SIZE;
          }

          if (damode & CGO_PICK_COLOR_ARRAY) {
            pickColorVals = nxtVals + VERTEX_PICKCOLOR_RGBA_SIZE * nverts;
            nxtVals += nverts * VERTEX_PICKCOLOR_SIZE;
          }

          if (damode & CGO_ACCESSIBILITY_ARRAY) {
            accessibilityVals = nxtVals;
            nxtVals += nverts * VERTEX_ACCESSIBILITY_SIZE;
          }

          auto notHaveValue = damode;

          // second iteration (with copy of iterator, doesn't consume 'it')
          for (; ok && it2 != CGO_END; ++it2) {
            auto pc = it2.data();
            switch (it2.op_code()) {
              case CGO_NORMAL:
                copy3f(pc, &normalVals[pl]);
                notHaveValue &= ~CGO_NORMAL_ARRAY;
                break;
              case CGO_COLOR:
                if (colorVals) {
                  copy3f(pc, &colorVals[plc]);
                  colorVals[plc + 3] = cgo->alpha;
                  notHaveValue &= ~CGO_COLOR_ARRAY;
                }
                copy3f(pc, cgo->color);
                break;
              case CGO_PICK_COLOR:
                cgo->current_pick_color_index = CGO_get_uint(pc);
                cgo->current_pick_color_bond = CGO_get_int(pc + 1);
                notHaveValue &= ~CGO_PICK_COLOR_ARRAY;
                break;
              case CGO_ACCESSIBILITY:
                cgo->current_accessibility = pc[0];
                break;
              case CGO_SPLITLINE: {
                auto splitline =
                    reinterpret_cast<const cgo::draw::splitline*>(pc);
                float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
                                  CONVERT_COLOR_VALUE(splitline->color2[1]),
                                  CONVERT_COLOR_VALUE(splitline->color2[2])};
                if (do_not_split_lines ||
                    (splitline->flags & cgo::draw::splitline::equal_colors)) {
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla,
                                           splitline->vertex1, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  pl += 3;
                  plc += 4;
                  pla++;
                  notHaveValue = damode;

                  if (!(splitline->flags &
                        cgo::draw::splitline::equal_colors)) {
                    if (colorVals) {
                      copy3f(color2, &colorVals[plc]);
                      colorVals[plc + 3] = cgo->alpha;
                      notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
                    }
                    copy3f(color2, cgo->color);
                  }
                  if (pickColorVals) {
                    cgo->current_pick_color_index = splitline->index;
                    cgo->current_pick_color_bond = splitline->bond;
                    notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
                  }
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla,
                                           splitline->vertex2, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  pl += 3;
                  plc += 4;
                  pla++;
                  notHaveValue = damode;
                } else {
                  float mid[3];
                  add3f(splitline->vertex1, splitline->vertex2, mid);
                  mult3f(mid, .5f, mid);
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla,
                                           splitline->vertex1, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  notHaveValue = damode;
                  pl += 3;
                  plc += 4;
                  pla++;
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  pl += 3;
                  plc += 4;
                  pla++;
                  if (colorVals) {
                    copy3f(color2, &colorVals[plc]);
                    colorVals[plc + 3] = cgo->alpha;
                    notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
                  }
                  copy3f(color2, cgo->color);
                  if (pickColorVals) {
                    cgo->current_pick_color_index = splitline->index;
                    cgo->current_pick_color_bond = splitline->bond;
                    notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
                  }
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  notHaveValue = damode;
                  pl += 3;
                  plc += 4;
                  pla++;
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla,
                                           splitline->vertex2, notHaveValue,
                                           vertexVals, normalVals, colorVals,
                                           pickColorVals, accessibilityVals);
                  pl += 3;
                  plc += 4;
                  pla++;
                  notHaveValue = damode;
                }
              } break;
              case CGO_LINE: {
                auto line = reinterpret_cast<const cgo::draw::line*>(pc);
                CGOAddVertexToDrawArrays(
                    cgo, pl, plc, pla, line->vertex1, notHaveValue, vertexVals,
                    normalVals, colorVals, pickColorVals, accessibilityVals);
                pl += 3;
                plc += 4;
                pla++;
                notHaveValue = damode;
                CGOAddVertexToDrawArrays(
                    cgo, pl, plc, pla, line->vertex2, notHaveValue, vertexVals,
                    normalVals, colorVals, pickColorVals, accessibilityVals);
                pl += 3;
                plc += 4;
                pla++;
              } break;
              case CGO_VERTEX:
                CGOAddVertexToDrawArrays(cgo, pl, plc, pla, pc, notHaveValue,
                                         vertexVals, normalVals, colorVals,
                                         pickColorVals, accessibilityVals);
                pl += 3;
                plc += 4;
                pla++;
                notHaveValue = damode;
                break;
              case CGO_ALPHA:
                // in case we're before CGO_COLOR
                cgo->alpha = *pc;
                if (colorVals) {
                  // in case we're after CGO_COLOR
                  colorVals[plc + 3] = *pc;
                }
                break;
            }
          }
        }
      } break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        cgo->add_to_cgo(op, pc);
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
      default:
        cgo->add_to_cgo(op, pc);
    }
  }
  if (ok) {
    ok &= CGOStop(cgo);
    if (ok) {
      cgo->use_shader = I->use_shader;
      if (cgo->use_shader) {
        cgo->cgo_shader_ub_color =
            SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
        cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
      }
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

#define set_min_max(mn, mx, pt)               \
  {                                           \
    if (mn[0] > *pt) mn[0] = *pt;             \
    if (mn[1] > *(pt + 1)) mn[1] = *(pt + 1); \
    if (mn[2] > *(pt + 2)) mn[2] = *(pt + 2); \
    if (mx[0] < *pt) mx[0] = *pt;             \
    if (mx[1] < *(pt + 1)) mx[1] = *(pt + 1); \
    if (mx[2] < *(pt + 2)) mx[2] = *(pt + 2); \
  }

struct CGOCount {
  int num_total_vertices{};
  int num_total_indexes{};
  int num_total_vertices_lines{};
  int num_total_indexes_lines{};
  int num_total_vertices_points{};
};

static CGOCount CGOCountNumVertices(const CGO* I);

static CGOCount CGOCountNumVertices(const CGO* I) {
  CGOCount count{};
  short err = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    err = 0;
    switch (op) {
      case CGO_DRAW_ARRAYS: {
        const auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false, shouldCompressLines = false,
              shouldCompressPoints = false;
        switch (sp->mode) {
          case GL_TRIANGLE_FAN:
          case GL_TRIANGLE_STRIP:
          case GL_TRIANGLES:
            shouldCompress = true;
            break;
          case GL_LINES:
          case GL_LINE_STRIP:
          case GL_LINE_LOOP:
            shouldCompressLines = true;
            break;
          case GL_POINTS:
            shouldCompressPoints = true;
            break;
          default:
            break;
        }
        if (!shouldCompress && !shouldCompressLines && !shouldCompressPoints) {
          // skipped
        } else if (shouldCompressLines) {
          count.num_total_vertices_lines += sp->nverts;
          switch (sp->mode) {
            case GL_LINE_LOOP:
              count.num_total_indexes_lines += 2 * sp->nverts;
              break;
            case GL_LINE_STRIP:
              count.num_total_indexes_lines += 2 * (sp->nverts - 1);
              break;
            case GL_LINES:
              count.num_total_indexes_lines += sp->nverts;
              break;
          }
        } else if (shouldCompress) {
          count.num_total_vertices += sp->nverts;
          switch (sp->mode) {
            case GL_TRIANGLE_FAN:
              count.num_total_indexes += 3 * (sp->nverts - 2);
              break;
            case GL_TRIANGLE_STRIP:
              count.num_total_indexes += 3 * (sp->nverts - 2);
              break;
            case GL_TRIANGLES:
              count.num_total_indexes += sp->nverts;
              break;
          }
        } else if (shouldCompressPoints) {
          count.num_total_vertices_points += sp->nverts;
        }
      } break;
      case CGO_END:
        if (!err) {
          PRINTFB(I->G, FB_CGO, FB_Warnings)
          " CGOCountNumVertices: CGO_END encountered, should call "
          "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
          err = true;
        }
      case CGO_VERTEX:
        if (!err) {
          PRINTFB(I->G, FB_CGO, FB_Warnings)
          " CGOCountNumVertices: CGO_VERTEX encountered, should call "
          "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
          err = true;
        }
      case CGO_BEGIN:
        if (!err) {
          PRINTFB(I->G, FB_CGO, FB_Warnings)
          " CGOCountNumVertices: CGO_BEGIN encountered, should call "
          "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
          err = true;
        }
      default:
        break;
    }
  }
  return count;
}

/**
 * @param[in] I Primitive CGO
 * @return number of vertices, or -1 on error
 */
static int CGOCountNumVerticesForScreen(const CGO* I) {
  auto G = I->G;
  constexpr int MODE_INVALID = -1;
  int begin_mode = MODE_INVALID;
  int num_total_indexes = 0;
  int nverts = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    // sanity check: non-primitive operations
    switch (op) {
      case CGO_DRAW_ARRAYS:
        WARN_UNEXPECTED_OPERATION(G, op);
        return -1;
    }

    if (begin_mode == MODE_INVALID) {
      switch (op) {
        case CGO_BEGIN:
          begin_mode = it.cast<cgo::draw::begin>()->mode;
          break;
        case CGO_VERTEX:
        case CGO_END:
          WARN_UNEXPECTED_OPERATION(G, op);
          return -1;
      }
    } else {
      switch (op) {
        case CGO_BEGIN:
          WARN_UNEXPECTED_OPERATION(G, op);
          return -1;
        case CGO_VERTEX:
          ++nverts;
          break;
        case CGO_END:
          switch (begin_mode) {
            case GL_TRIANGLE_FAN:
            case GL_TRIANGLE_STRIP:
              num_total_indexes += 3 * (nverts - 2);
              break;
            case GL_TRIANGLES:
              num_total_indexes += nverts;
              break;
            default:
              // not implemented
              assert(false);
          }
          begin_mode = MODE_INVALID;
          nverts = 0;
          break;
      }
    }
  }

  return num_total_indexes;
}

static void SetVertexValuesForVBO(
    PyMOLGlobals* G, CGO* cgo, int pl, int plc, int cnt, int incr,
    const float* vertexValsDA, const float* normalValsDA,
    const float* colorValsDA, const float* pickColorValsDA, float* vertexVals,
    uchar* normalValsC, float* normalVals, uchar* colorValsUC, float* colorVals,
    float* pickColorVals, float* accessibilityVals = nullptr,
    const float* accessibilityValsDA = nullptr) {
  int pl2 = pl + 1, pl3 = pl + 2;
  int pln1 = VAR_FOR_NORMAL, pln2 = VAR_FOR_NORMAL + 1,
      pln3 = VAR_FOR_NORMAL + 2;
  int plc2 = plc + 1, plc3 = plc + 2, plc4 = plc + 3;
  int c, c2, c3;
  int cc, cc2, cc3, cc4;
  int pcc = incr * 2, pcco = cnt * 2;
  c = cnt * 3;
  c2 = c + 1;
  c3 = c + 2;
  cc = cnt * 4;
  cc2 = cc + 1;
  cc3 = cc + 2;
  cc4 = cc + 3;
  vertexVals[pl] = vertexValsDA[c];
  vertexVals[pl2] = vertexValsDA[c2];
  vertexVals[pl3] = vertexValsDA[c3];

  if (normalValsC) {
    if (normalValsDA) {
      normalValsC[pln1] = CLIP_NORMAL_VALUE(normalValsDA[c]);
      normalValsC[pln2] = CLIP_NORMAL_VALUE(normalValsDA[c2]);
      normalValsC[pln3] = CLIP_NORMAL_VALUE(normalValsDA[c3]);
    } else {
      normalValsC[pln1] = CLIP_NORMAL_VALUE(cgo->normal[0]);
      normalValsC[pln2] = CLIP_NORMAL_VALUE(cgo->normal[1]);
      normalValsC[pln3] = CLIP_NORMAL_VALUE(cgo->normal[2]);
    }

#ifdef ALIGN_VBOS_TO_4_BYTE_ARRAYS
    normalValsC[pln3 + 1] = 127;
#endif
  } else {
    if (normalValsDA) {
      normalVals[pln1] = normalValsDA[c];
      normalVals[pln2] = normalValsDA[c2];
      normalVals[pln3] = normalValsDA[c3];
    } else {
      normalVals[pln1] = cgo->normal[0];
      normalVals[pln2] = cgo->normal[1];
      normalVals[pln3] = cgo->normal[2];
    }
  }

  if (colorValsUC) {
    if (colorValsDA) {
      colorValsUC[plc] = CLIP_COLOR_VALUE(colorValsDA[cc]);
      colorValsUC[plc2] = CLIP_COLOR_VALUE(colorValsDA[cc2]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(colorValsDA[cc3]);
      colorValsUC[plc4] = CLIP_COLOR_VALUE(colorValsDA[cc4]);
    } else {
      colorValsUC[plc] = CLIP_COLOR_VALUE(cgo->color[0]);
      colorValsUC[plc2] = CLIP_COLOR_VALUE(cgo->color[1]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(cgo->color[2]);
      colorValsUC[plc4] = CLIP_COLOR_VALUE(cgo->alpha);
    }
  } else {
    if (colorValsDA) {
      colorVals[plc] = colorValsDA[cc];
      colorVals[plc2] = colorValsDA[cc2];
      colorVals[plc3] = colorValsDA[cc3];
      colorVals[plc4] = colorValsDA[cc4];
    } else {
      colorVals[plc] = cgo->color[0];
      colorVals[plc2] = cgo->color[1];
      colorVals[plc3] = cgo->color[2];
      colorVals[plc4] = cgo->alpha;
    }
  }
  if (pickColorValsDA) {
    cgo->current_pick_color_index = CGO_get_uint(pickColorValsDA + pcco);
    cgo->current_pick_color_bond = CGO_get_int(pickColorValsDA + pcco + 1);
  }
  CGO_put_uint(pickColorVals + pcc, cgo->current_pick_color_index);
  CGO_put_int(pickColorVals + pcc + 1, cgo->current_pick_color_bond);
  if (accessibilityValsDA) {
    accessibilityVals[pl / 3] = accessibilityValsDA[cnt];
  }
}

struct NormalColorFormatSize {
  VertexFormat normalFormat{};
  VertexFormat colorFormat{};
  std::size_t normalSize{};
  std::size_t colorSize{};
};

static NormalColorFormatSize GetNormalColorFormatSize(PyMOLGlobals* G) {
  NormalColorFormatSize fmt{};

  fmt.normalFormat =
      VERTEX_NORMAL_SIZE == 3 ? VertexFormat::Float3 : VertexFormat::Float4;
  if (CGOPackedNormals(G)) {
    fmt.normalFormat = VERTEX_NORMAL_SIZE == 3 ? VertexFormat::Byte3Norm
                                               : VertexFormat::Byte4Norm;
  }
  fmt.normalSize = GetSizeOfVertexFormat(fmt.normalFormat);

  fmt.colorFormat = VertexFormat::Float4;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    fmt.colorFormat = VertexFormat::UByte4Norm;
  }
  fmt.colorSize = GetSizeOfVertexFormat(fmt.colorFormat);

  return fmt;
}

// TODO: The following functions are temporary until OpenGL is a backend

/**
 * @brief Create a vertex buffer
 * @param layout Vertex buffer layout
 * @param memProperty Usage property for Vertex buffer
 * @return Reference to created Vertex Buffer (no mem allocated)
 */
VertexBuffer& CreateVertexBuffer(PyMOLGlobals* G, VertexBufferLayout layout,
                                 MemoryUsageProperty memProperty) {
  return *G->ShaderMgr->newGPUBuffer<VertexBufferGL>(layout, memProperty);
}

/**
 * @brief Allocates GPU memory for vertex buffer and sets data
 * @param dataDesc Buffer data description
 * @param dataOverride Optional data override
 */
void SetVertexBufferData(PyMOLGlobals* G, VertexBuffer& vertexBuffer,
                         BufferDataDesc&& dataDesc,
                         pymol::span<const std::byte> dataOverride) {
  auto* vbo = static_cast<VertexBufferGL*>(&vertexBuffer);
  vbo->bufferData(std::move(dataDesc));
}

/**
 * @brief Allocates GPU memory for vertex buffer and sets data
 * @param vertexBuffer Vertex buffer
 * @param dataDesc Buffer data description
 * @param bufferSize Size of the buffer
 * @param dataOverride Optional data override
 */
void SetVertexBufferDataSized(PyMOLGlobals* G, VertexBuffer& vertexBuffer,
                              BufferDataDesc&& dataDesc, std::size_t bufferSize,
                              std::byte* dataOverride) {
  auto* vbo = static_cast<VertexBufferGL*>(&vertexBuffer);
  vbo->bufferData(std::move(dataDesc), dataOverride, bufferSize);
}

/**
 * @brief Create an index buffer
 * @return Reference to created Index Buffer (no mem allocated)
 */
IndexBuffer& CreateIndexBuffer(PyMOLGlobals* G) {
  return *G->ShaderMgr->newGPUBuffer<IndexBufferGL>();
}

/**
 * @brief Allocates GPU memory for index buffer and sets data
 * @param indexBuffer Index buffer
 * @param data Index data
 */
void SetIndexBufferData(PyMOLGlobals* G, IndexBuffer& indexBuffer,
                        pymol::span<const std::uint32_t> data) {
  auto* ibo = static_cast<IndexBufferGL*>(&indexBuffer);
  ibo->copyFrom(data);
}

static int OptimizePointsToVBO(const CGO* I, CGO* cgo,
                               int num_total_vertices_points, float* min,
                               float* max, short* has_draw_buffer,
                               bool addshaders) {
  auto G = I->G;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  bool has_normals = false, has_colors = false;
  int ok = true;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  mul += CGOPackedNormals(G) ? 1 : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(num_total_vertices_points) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * num_total_vertices_points;
  nxtn = 3;
  if (CGOPackedNormals(I->G)) {
    normalValsC = (uchar*)normalVals;
    nxtn = 1;
  }
  auto* colorVals = normalVals + nxtn * num_total_vertices_points;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*)colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * num_total_vertices_points);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_BOUNDING_BOX:
      case CGO_DRAW_SPHERE_BUFFERS:
      case CGO_DRAW_LABELS:
      case CGO_DRAW_TEXTURES:
      case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      case CGO_DRAW_CYLINDER_BUFFERS:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
      case CGO_DRAW_BUFFERS_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Errors)
        "ERROR: OptimizePointsToVBO used with unsupported CGO ops" ENDFB(I->G);
        return 0;
        break;
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        has_normals = true;
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        has_colors = true;
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_ARRAYS: {
        const auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
        switch (sp->mode) {
          case GL_POINTS:
            shouldCompress = true;
          default:
            break;
        }
        /*	TODO : DO WE NEED TO COMPENSATE FOR THIS? if (!has_normals &&
          arrays & CGO_NORMAL_ARRAY){ arrays = arrays ^ CGO_NORMAL_ARRAY;
          narrays
          -= 1;
          }*/
        if (shouldCompress) {
          int cnt, nxtn = VERTEX_POS_SIZE, incr = 0;
          float *vertexValsDA = nullptr, *nxtVals = nullptr,
                *colorValsDA = nullptr, *normalValsDA = nullptr;
          float *pickColorValsDA = nullptr, *pickColorValsTMP;

          nxtVals = vertexValsDA = sp->floatdata;

          for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
            set_min_max(min, max, &vertexValsDA[cnt]);
          }
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            has_normals = true;
            nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
          }
          if (sp->arraybits & CGO_COLOR_ARRAY) {
            has_colors = true;
            nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
            nxtn = VERTEX_COLOR_SIZE;
          }
          if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
            nxtVals = nxtVals + (nxtn * sp->nverts);
            pickColorValsDA = nxtVals + sp->nverts;
            nxtn = VERTEX_PICKCOLOR_SIZE;
          }
          pickColorValsTMP = pickColorVals + (idxpl * 2);
          switch (sp->mode) {
            case GL_POINTS:
              for (cnt = 0; cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(
                    I->G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP);
                idxpl++;
                pl += 3;
                plc += 4;
              }
              break;
          }
          vpl += sp->nverts;
        }
      } break;
    }
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    short arrays = CGO_VERTEX_ARRAY | CGO_PICK_COLOR_ARRAY;
    auto fmt = GetNormalColorFormatSize(I->G);

    auto& vbo = CreateVertexBuffer(I->G);

    BufferDataDesc bufData{
        {{"a_Vertex", VertexFormat::Float3,
          sizeof(float) * num_total_vertices_points * 3, vertexVals}}};

    if (has_normals) {
      bufData.descs.push_back({"a_Normal", fmt.normalFormat,
                               num_total_vertices_points * fmt.normalSize,
                               normalVals});
    }
    if (has_colors) {
      bufData.descs.push_back({"a_Color", fmt.colorFormat,
                               num_total_vertices_points * fmt.colorSize,
                               colorVals});
    }
    SetVertexBufferData(I->G, vbo, std::move(bufData));

    if (ok && has_colors) {
      arrays |= CGO_COLOR_ARRAY;
    }

    auto vboid = vbo.get_hash_id();

    if (ok) {
      float* newPickColorVals;
      if (addshaders) CGOEnable(cgo, GL_DEFAULT_SHADER);
      newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(
          GL_POINTS, arrays, num_total_vertices_points, vboid);
      CHECKOK(ok, newPickColorVals);
      if (ok && addshaders) ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      if (!newPickColorVals) I->G->ShaderMgr->freeGPUBuffer(vboid);
      if (ok)
        memcpy(newPickColorVals + num_total_vertices_points, pickColorVals,
               num_total_vertices_points * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(vboid);
    }
  }
  return ok;
  /* END GL_POINTS */
  //    printf("num_total_vertices_points=%d\n", num_total_vertices_points);
}

static void FixPickColorsForLine(float* pick1, float* pick2) {
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  if (p1 != p2 || b1 != b2) {
    // if the pick colors are different, then pick the first one
    CGO_put_uint(pick1, p2);
    CGO_put_int(pick1 + 1, b2);
  }
}

static void FixPickColorsForTriangle(float* pick1, float* pick2, float* pick3) {
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  unsigned int p3 = CGO_get_uint(pick3);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  int b3 = CGO_get_int(pick3 + 1);
  if (p1 != p2 || p1 != p3 || p2 != p3 || b1 != b2 || b1 != b3 || b2 != b3) {
    // right now, if the pick colors are different, then pick majority,
    // otherwise, pick first one
    if (p1 == p2 && b1 == b2) {
      CGO_put_uint(pick3, p1);
      CGO_put_int(pick3 + 1, b1);
    } else if (p1 == p3 && b1 == b3) {
      CGO_put_uint(pick2, p1);
      CGO_put_int(pick2 + 1, b1);
    } else if (p2 == p3 && b2 == b3) {
      CGO_put_uint(pick1, p2);
      CGO_put_int(pick1 + 1, b2);
    } else {
      CGO_put_uint(pick2, p1);
      CGO_put_int(pick2 + 1, b1);
      CGO_put_uint(pick3, p1);
      CGO_put_int(pick3 + 1, b1);
    }
  }
}

static bool CGOProcessCGOtoArrays(const CGO* I, CGO* cgo, CGO* addtocgo,
                                  float* min, float* max,
                                  int* ambient_occlusion, float* vertexVals,
                                  float* normalVals, uchar* normalValsC,
                                  float* colorVals, uchar* colorValsUC,
                                  float* pickColorVals,
                                  float* accessibilityVals, bool& has_normals,
                                  bool& has_colors, bool& has_accessibility) {
  auto G = I->G;
  int idxpl = 0;
  int pl = 0, plc = 0, vpl = 0;
  int ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_BOUNDING_BOX:
      case CGO_DRAW_SPHERE_BUFFERS:
      case CGO_DRAW_LABELS:
      case CGO_DRAW_TEXTURES:
      case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      case CGO_DRAW_CYLINDER_BUFFERS:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
      case CGO_DRAW_BUFFERS_INDEXED: {
        const float* newpc = pc;
        if (addtocgo) addtocgo->add_to_cgo(op, newpc);
      } break;
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        has_normals = true;
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        has_colors = true;
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_ACCESSIBILITY:
        cgo->current_accessibility = pc[0];
        has_accessibility = true;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_ARRAYS: {
        const auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
        switch (sp->mode) {
          case GL_TRIANGLE_FAN:
          case GL_TRIANGLE_STRIP:
          case GL_TRIANGLES:
            shouldCompress = true;
          default:
            break;
        }
        if (shouldCompress) {
          int cnt, incr = 0;
          const float* nxtVals = sp->floatdata;
          const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                      *normalValsDA = nullptr, *pickColorValsDA = nullptr,
                      *accessibilityValsDA = nullptr;

          assert(sp->arraybits & CGO_VERTEX_ARRAY);
          vertexValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_POS_SIZE;

          for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
            set_min_max(min, max, &vertexValsDA[cnt]);
          }
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            has_normals = true;
            normalValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
          }

          if (sp->arraybits & CGO_COLOR_ARRAY) {
            has_colors = true;
            colorValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
          }
          if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
            nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
            pickColorValsDA = nxtVals;
            nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
          }
          float* pickColorValsTMP = pickColorVals + (idxpl * 2);
          if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY) {
            if (!(*ambient_occlusion) && incr) {
              for (cnt = 0; cnt < incr; cnt++) {
                /* if ambient_occlusion, need to fill in the array */
                accessibilityVals[cnt] = 1.f;
              }
            }
            (*ambient_occlusion) = 1;
            accessibilityValsDA = nxtVals;
            nxtVals += VERTEX_ACCESSIBILITY_SIZE * sp->nverts;
            has_accessibility = true;
          } else {
            if (*ambient_occlusion) {
              for (cnt = incr; cnt < incr + sp->nverts; cnt++) {
                /* if ambient_occlusion, need to fill in the array */
                accessibilityVals[cnt] = 1.f;
              }
            }
          }
          switch (sp->mode) {
            case GL_TRIANGLES:
              for (cnt = 0; ok && cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP,
                    accessibilityVals, accessibilityValsDA);
                if (incr && (incr % 3) == 0) {
                  FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                                           pickColorValsTMP + (incr - 2) * 2,
                                           pickColorValsTMP + (incr - 1) * 2);
                }
                idxpl++;
                pl += 3;
                plc += 4;
                ok &= !G->Interrupt;
              }
              break;
            case GL_TRIANGLE_STRIP: {
              short flip = 0;
              for (cnt = 2; ok && cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt - (flip ? 0 : 2), incr++, vertexValsDA,
                    normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                    normalValsC, normalVals, colorValsUC, colorVals,
                    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt - 1, incr++, vertexValsDA,
                    normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                    normalValsC, normalVals, colorValsUC, colorVals,
                    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt - (flip ? 2 : 0), incr++, vertexValsDA,
                    normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                    normalValsC, normalVals, colorValsUC, colorVals,
                    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
                FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                                         pickColorValsTMP + (incr - 2) * 2,
                                         pickColorValsTMP + (incr - 1) * 2);
                idxpl++;
                pl += 3;
                plc += 4;
                ok &= !G->Interrupt;
                flip = !flip;
              }
            } break;
            case GL_TRIANGLE_FAN:
              for (cnt = 2; ok && cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, 0, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP,
                    accessibilityVals, accessibilityValsDA);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt - 1, incr++, vertexValsDA,
                    normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                    normalValsC, normalVals, colorValsUC, colorVals,
                    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP,
                    accessibilityVals, accessibilityValsDA);
                FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                                         pickColorValsTMP + (incr - 2) * 2,
                                         pickColorValsTMP + (incr - 1) * 2);
                idxpl++;
                pl += 3;
                plc += 4;
                ok &= !G->Interrupt;
              }
              break;
          }
          vpl += sp->nverts;
        }
      } break;
    }
    ok &= !G->Interrupt;
  }
  ok &= !G->Interrupt;
  return ok;
}

static bool CGOProcessScreenCGOtoArrays(PyMOLGlobals* G, CGO* cgo,
                                        float* vertexVals, float* texcoordVals,
                                        float* colorVals, uchar* colorValsUC) {
  int pl = 0;
  cgo->alpha = 1.f;

  for (auto it = cgo->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_BOUNDING_BOX:
      case CGO_DRAW_SPHERE_BUFFERS:
      case CGO_DRAW_LABELS:
      case CGO_DRAW_TEXTURES:
      case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      case CGO_DRAW_CYLINDER_BUFFERS:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_ARRAYS:
      case CGO_ACCESSIBILITY:
        WARN_UNEXPECTED_OPERATION(G, op);
        return false;
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        break;
      case CGO_TEX_COORD:
        cgo->texture[0] = *pc;
        cgo->texture[1] = *(pc + 1);
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_BEGIN: {
        int mode = it.cast<cgo::draw::begin>()->mode;
        int nverts = 0, ipl = 0;
        cgo->texture[0] = cgo->texture[1] = 0.f;

        for (++it;; ++it) {
          if (it.is_stop()) {
            WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
            return false;
          }

          const auto op = it.op_code();
          if (op == CGO_END) {
            break;
          }

          const auto pc = it.data();

          switch (op) {
            case CGO_TEX_COORD:
              cgo->texture[0] = *pc;
              cgo->texture[1] = *(pc + 1);
              break;
            case CGO_COLOR:
              cgo->color[0] = *pc;
              cgo->color[1] = *(pc + 1);
              cgo->color[2] = *(pc + 2);
              break;
            case CGO_ALPHA:
              cgo->alpha = *pc;
              break;
            case CGO_VERTEX: {
              switch (mode) {
                case GL_TRIANGLES: {
                  int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
                  vertexVals[vpl] = *pc;
                  vertexVals[vpl + 1] = *(pc + 1);
                  vertexVals[vpl + 2] = *(pc + 2);
                  texcoordVals[tpl] = cgo->texture[0];
                  texcoordVals[tpl + 1] = cgo->texture[1];
                  if (colorValsUC) {
                    colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
                    colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
                    colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
                    colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
                  } else {
                    colorVals[cpl] = cgo->color[0];
                    colorVals[cpl + 1] = cgo->color[1];
                    colorVals[cpl + 2] = cgo->color[2];
                    colorVals[cpl + 3] = cgo->alpha;
                  }
                  pl++;
                } break;
                case GL_TRIANGLE_STRIP: {
                  int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
                  if (ipl < 3) {
                    vertexVals[vpl] = *pc;
                    vertexVals[vpl + 1] = *(pc + 1);
                    vertexVals[vpl + 2] = *(pc + 2);
                    texcoordVals[tpl] = cgo->texture[0];
                    texcoordVals[tpl + 1] = cgo->texture[1];
                    if (colorValsUC) {
                      colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
                      colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
                      colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
                      colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
                    } else {
                      colorVals[cpl] = cgo->color[0];
                      colorVals[cpl + 1] = cgo->color[1];
                      colorVals[cpl + 2] = cgo->color[2];
                      colorVals[cpl + 3] = cgo->alpha;
                    }
                    pl++;
                    ipl++;
                  } else {
                    vertexVals[vpl] = vertexVals[vpl - 6];
                    vertexVals[vpl + 1] = vertexVals[vpl - 5];
                    vertexVals[vpl + 2] = vertexVals[vpl - 4];
                    texcoordVals[tpl] = texcoordVals[tpl - 4];
                    texcoordVals[tpl + 1] = texcoordVals[tpl - 3];
                    if (colorValsUC) {
                      colorValsUC[cpl] = colorValsUC[cpl - 8];
                      colorValsUC[cpl + 1] = colorValsUC[cpl - 7];
                      colorValsUC[cpl + 2] = colorValsUC[cpl - 6];
                      colorValsUC[cpl + 3] = colorValsUC[cpl - 5];
                    } else {
                      colorVals[cpl] = colorVals[cpl - 8];
                      colorVals[cpl + 1] = colorVals[cpl - 7];
                      colorVals[cpl + 2] = colorVals[cpl - 6];
                      colorVals[cpl + 3] = colorVals[cpl - 5];
                    }
                    pl++;
                    vpl += 3;
                    tpl += 2;
                    cpl += 4;
                    ipl++;
                    vertexVals[vpl] = vertexVals[vpl - 6];
                    vertexVals[vpl + 1] = vertexVals[vpl - 5];
                    vertexVals[vpl + 2] = vertexVals[vpl - 4];
                    texcoordVals[tpl] = texcoordVals[tpl - 4];
                    texcoordVals[tpl + 1] = texcoordVals[tpl - 3];
                    if (colorValsUC) {
                      colorValsUC[cpl] = colorValsUC[cpl - 8];
                      colorValsUC[cpl + 1] = colorValsUC[cpl - 7];
                      colorValsUC[cpl + 2] = colorValsUC[cpl - 6];
                      colorValsUC[cpl + 3] = colorValsUC[cpl - 5];
                    } else {
                      colorVals[cpl] = colorVals[cpl - 8];
                      colorVals[cpl + 1] = colorVals[cpl - 7];
                      colorVals[cpl + 2] = colorVals[cpl - 6];
                      colorVals[cpl + 3] = colorVals[cpl - 5];
                    }
                    pl++;
                    vpl += 3;
                    tpl += 2;
                    cpl += 4;
                    ipl++;
                    vertexVals[vpl] = *pc;
                    vertexVals[vpl + 1] = *(pc + 1);
                    vertexVals[vpl + 2] = *(pc + 2);
                    texcoordVals[tpl] = cgo->texture[0];
                    texcoordVals[tpl + 1] = cgo->texture[1];
                    if (colorValsUC) {
                      colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
                      colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
                      colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
                      colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
                    } else {
                      colorVals[cpl] = cgo->color[0];
                      colorVals[cpl + 1] = cgo->color[1];
                      colorVals[cpl + 2] = cgo->color[2];
                      colorVals[cpl + 3] = cgo->alpha;
                    }
                    pl++;
                    ipl++;
                  }
                } break;
                case GL_TRIANGLE_FAN:
                default:
                  printf(
                      "CGOProcessScreenCGOtoArrays: WARNING: mode=%d not "
                      "implemented "
                      "yet GL_LINES=%d GL_LINE_STRIP=%d GL_LINE_LOOP=%d\n",
                      mode, GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP);
                  break;
              }
              nverts++;
            }
          }
        }
      } break;
    }
  }
  return true;
}

static bool OptimizeVertsToVBONotIndexed(const CGO* I, CGO* cgo,
                                         const CGOCount& count, float* min,
                                         float* max, short* has_draw_buffer,
                                         bool addshaders) {
  auto G = I->G;
  bool ok = true;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul =
      VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
  mul += CGOPackedNormals(G) ? 1 : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_indexes) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_indexes;
  unsigned nxtn = VERTEX_NORMAL_SIZE;
  if (CGOPackedNormals(G)) {
    normalValsC = (uchar*)normalVals;
    nxtn = 1;
  }
  auto* colorVals = normalVals + nxtn * count.num_total_indexes;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*)colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_indexes);
  nxtn = 3;
  auto* accessibilityVals = pickColorVals + nxtn * count.num_total_indexes;

  bool has_normals = false, has_colors = false, has_accessibility = false;
  int ambient_occlusion{};
  ok = CGOProcessCGOtoArrays(I, cgo, cgo, min, max, &ambient_occlusion,
                             vertexVals, normalVals, normalValsC, colorVals,
                             colorValsUC, pickColorVals, accessibilityVals,
                             has_normals, has_colors, has_accessibility);
  if (!ok) {
    if (!G->Interrupt) PRINTFB(G, FB_CGO, FB_Errors)
    "ERROR: CGOProcessCGOtoArrays() could not allocate enough "
    "memory\n" ENDFB(G);
    return false;
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(G);

    auto& vbo = CreateVertexBuffer(G, VertexBufferLayout::Sequential);
    BufferDataDesc bufData{
        {{"a_Vertex", VertexFormat::Float3,
          sizeof(float) * count.num_total_indexes * 3, vertexVals}}};
    if (has_normals) {
      bufData.descs.push_back({"a_Normal", fmt.normalFormat,
                               count.num_total_indexes * fmt.normalSize,
                               normalVals});
    }
    if (has_colors) {
      bufData.descs.push_back({"a_Color", fmt.colorFormat,
                               count.num_total_indexes * fmt.colorSize,
                               colorVals});
    }
    if (has_accessibility) {
      bufData.descs.push_back({"a_Accessibility", VertexFormat::Float,
                               sizeof(float) * count.num_total_indexes,
                               accessibilityVals});
    }
    SetVertexBufferData(I->G, vbo, std::move(bufData));

    auto vboid = vbo.get_hash_id();
    // picking VBO: generate a buffer twice the size needed, for each picking
    // pass
    auto& pickvbo = CreateVertexBuffer(G, VertexBufferLayout::Sequential,
                                       MemoryUsageProperty::CpuToGpu);
    BufferDataDesc pickDataDesc{
        {//
         BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes},
         BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes}}};
    SetVertexBufferData(G, pickvbo, std::move(pickDataDesc));
    auto pickvboid = pickvbo.get_hash_id();

    if (ok) {
      float* newPickColorVals;
      int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
                   CGO_PICK_COLOR_ARRAY;
      if (ambient_occlusion) {
        arrays |= CGO_ACCESSIBILITY_ARRAY;
      }
      if (addshaders) CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
      newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(
          GL_TRIANGLES, arrays, count.num_total_indexes, vboid, pickvboid);
      if (ok && addshaders) ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      CHECKOK(ok, newPickColorVals);
      if (!newPickColorVals) {
        G->ShaderMgr->freeGPUBuffer(pickvboid);
        G->ShaderMgr->freeGPUBuffer(vboid);
      }
      if (!ok) {
        PRINTFB(G, FB_CGO, FB_Errors)
        "CGOOptimizeToVBONotIndexed: ERROR: "
        "CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(
            G);
        return false;
      }
      memcpy(newPickColorVals + count.num_total_indexes, pickColorVals,
             count.num_total_indexes * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      G->ShaderMgr->freeGPUBuffer(vboid);
      G->ShaderMgr->freeGPUBuffer(pickvboid);
    }
  }
  return ok;
}

static bool OptimizeLinesToVBONotIndexed(const CGO* I, CGO* cgo,
                                         const CGOCount& count, float* min,
                                         float* max, short* has_draw_buffer,
                                         bool addshaders) {
  auto G = I->G;
  bool ok = true;
  bool has_color = false, has_normals = false;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  mul += CGOPackedNormals(G) ? 1 : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_indexes_lines) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_indexes_lines;
  nxtn = 3;
  if (CGOPackedNormals(G)) {
    normalValsC = (uchar*)normalVals;
    nxtn = 1;
  }

  auto* colorVals = normalVals + nxtn * count.num_total_indexes_lines;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*)colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_indexes_lines);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_SPECIAL:
      case CGO_RESET_NORMAL: {
        const float* newpc = pc;
        cgo->add_to_cgo(op, newpc);
      } break;
      case CGO_NORMAL:
        has_normals = true;
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        break;
      case CGO_COLOR:
        has_color = true;
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        break;
      case CGO_ACCESSIBILITY:
        cgo->current_accessibility = pc[0];
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_ARRAYS: {
        auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
        switch (sp->mode) {
          case GL_LINE_LOOP:
          case GL_LINE_STRIP:
          case GL_LINES:
            shouldCompress = true;
          default:
            break;
        }

        if (shouldCompress) {
          int cnt, incr = 0;
          const float* nxtVals = sp->floatdata;
          const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                      *normalValsDA = nullptr, *pickColorValsDA = nullptr;

          assert(sp->arraybits & CGO_VERTEX_ARRAY);
          vertexValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_POS_SIZE;

          for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
            set_min_max(min, max, &vertexValsDA[cnt]);
          }
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            has_normals = true;
            normalValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
          }

          if (sp->arraybits & CGO_COLOR_ARRAY) {
            has_color = true;
            colorValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
          }
          if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
            nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
            pickColorValsDA = nxtVals;
            nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
          }
          float* pickColorValsTMP = pickColorVals + (idxpl * 2);
          switch (sp->mode) {
            case GL_LINES:
              for (cnt = 0; cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP);
                if (incr && (incr % 2) == 0) {
                  FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                                       pickColorValsTMP + (incr - 1) * 2);
                }
                idxpl++;
                pl += 3;
                plc += 4;
              }
              break;
            case GL_LINE_STRIP:
              for (cnt = 1; cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                                      vertexValsDA, normalValsDA, colorValsDA,
                                      pickColorValsDA, vertexVals, normalValsC,
                                      normalVals, colorValsUC, colorVals,
                                      pickColorValsTMP);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP);
                FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                                     pickColorValsTMP + (incr - 1) * 2);
                idxpl++;
                pl += 3;
                plc += 4;
              }
              break;
            case GL_LINE_LOOP:
              for (cnt = 1; cnt < sp->nverts; cnt++) {
                SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                                      vertexValsDA, normalValsDA, colorValsDA,
                                      pickColorValsDA, vertexVals, normalValsC,
                                      normalVals, colorValsUC, colorVals,
                                      pickColorValsTMP);
                idxpl++;
                pl += 3;
                plc += 4;
                SetVertexValuesForVBO(
                    G, cgo, pl, plc, cnt, incr++, vertexValsDA, normalValsDA,
                    colorValsDA, pickColorValsDA, vertexVals, normalValsC,
                    normalVals, colorValsUC, colorVals, pickColorValsTMP);
                FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                                     pickColorValsTMP + (incr - 1) * 2);
                idxpl++;
                pl += 3;
                plc += 4;
              }
              SetVertexValuesForVBO(G, cgo, pl, plc, 0, incr++, vertexValsDA,
                                    normalValsDA, colorValsDA, pickColorValsDA,
                                    vertexVals, normalValsC, normalVals,
                                    colorValsUC, colorVals, pickColorValsTMP);
              idxpl++;
              pl += 3;
              plc += 4;
              SetVertexValuesForVBO(G, cgo, pl, plc, sp->nverts - 1, incr++,
                                    vertexValsDA, normalValsDA, colorValsDA,
                                    pickColorValsDA, vertexVals, normalValsC,
                                    normalVals, colorValsUC, colorVals,
                                    pickColorValsTMP);
              FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                                   pickColorValsTMP + (incr - 1) * 2);
              idxpl++;
              pl += 3;
              plc += 4;
              break;
          }

          //	  pl += 3 * nverts;
          //	  plc += 4 * nverts;
          vpl += sp->nverts;
        }
      } break;
    }
  }
  auto fmt = GetNormalColorFormatSize(G);

  auto& vbo = CreateVertexBuffer(G, VertexBufferLayout::Sequential);
  BufferDataDesc bufData{
      {{"a_Vertex", VertexFormat::Float3,
        sizeof(float) * count.num_total_indexes_lines * 3, vertexVals}}};

  if (has_normals) {
    bufData.descs.push_back({"a_Normal", fmt.normalFormat,
                             count.num_total_indexes_lines * fmt.normalSize,
                             normalVals});
  }
  if (has_color) {
    bufData.descs.push_back({"a_Color", fmt.colorFormat,
                             count.num_total_indexes_lines * fmt.colorSize,
                             colorVals});
  }
  SetVertexBufferData(I->G, vbo, std::move(bufData));
  auto vboid = vbo.get_hash_id();

  // picking VBO: generate a buffer twice the size needed, for each picking
  // pass
  auto& pickvbo = CreateVertexBuffer(G, VertexBufferLayout::Sequential,
                                     MemoryUsageProperty::CpuToGpu);
  BufferDataDesc pickDataDesc{
      {//
       BufferDesc("a_Color", VertexFormat::UByte4Norm,
                  sizeof(float) * count.num_total_indexes_lines),
       BufferDesc("a_Color", VertexFormat::UByte4Norm,
                  sizeof(float) * count.num_total_indexes_lines)}};
  SetVertexBufferData(G, pickvbo, std::move(pickDataDesc));
  auto pickvboid = pickvbo.get_hash_id();

  if (ok) {
    float* newPickColorVals;
    if (addshaders) CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
    CGODisable(cgo, GL_SHADER_LIGHTING);
    newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(
        GL_LINES,
        CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
            CGO_PICK_COLOR_ARRAY,
        count.num_total_indexes_lines, vboid, pickvboid);
    if (ok && addshaders) ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
    CHECKOK(ok, newPickColorVals);
    if (!ok) {
      PRINTFB(G, FB_CGO, FB_Errors)
      "CGOOptimizeToVBONotIndexed: ERROR: "
      "CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(G);
      if (!newPickColorVals) {
        G->ShaderMgr->freeGPUBuffer(pickvboid);
        G->ShaderMgr->freeGPUBuffer(vboid);
      }
      return false;
    }
    memcpy(newPickColorVals + count.num_total_indexes_lines, pickColorVals,
           count.num_total_indexes_lines * 2 * sizeof(float));
    *has_draw_buffer = true;
  } else {
    G->ShaderMgr->freeGPUBuffer(pickvboid);
    G->ShaderMgr->freeGPUBuffer(vboid);
  }
  return ok;
}

bool CGOOptimizeToVBONotIndexed(CGO** I) {
  CGO* cgo = CGOOptimizeToVBONotIndexed(*I, 0, true);
  CGOFree(*I);
  *I = cgo;
  return (cgo != nullptr);
}

/**
 * Converts primitive operations into VBO operations.
 *
 * If the input CGO has sortable alpha triangles, you should use
 * CGOOptimizeToVBOIndexed() instead.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param addshaders Add Enable/Disable shader operations
 *
 * @return New converted CGO with use_shader=true
 */
CGO* CGOOptimizeToVBONotIndexed(const CGO* I, int est, bool addshaders) {
  auto G = I->G;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  short has_draw_buffer = false;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX},
        max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int ok = true;
  auto cgo = CGONew(G);

  auto count = CGOCountNumVertices(I);
  if (count.num_total_vertices_points > 0) {
    if (!OptimizePointsToVBO(I, cgo, count.num_total_vertices_points, min, max,
                             &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (count.num_total_indexes > 0) {
    if (!OptimizeVertsToVBONotIndexed(I, cgo, count, min, max, &has_draw_buffer,
                                      addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok && count.num_total_indexes_lines > 0) {
    if (!OptimizeLinesToVBONotIndexed(I, cgo, count, min, max, &has_draw_buffer,
                                      addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }

  if (ok &&
      (count.num_total_vertices > 0 || count.num_total_vertices_lines > 0 ||
       count.num_total_vertices_points > 0)) {
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok) ok &= CGOStop(cgo);

  if (!ok) {
    CGOFree(cgo);
    return nullptr;
  }

  cgo->has_draw_buffers |= has_draw_buffer;
  cgo->use_shader = true;
  cgo->cgo_shader_ub_color = SettingGet<int>(G, cSetting_cgo_shader_ub_color);
  cgo->cgo_shader_ub_normal = CGOPackedNormals(G);
  return (cgo);
}

static bool OptimizeVertsToVBOIndexed(const CGO* I, CGO* cgo,
                                      const CGOCount& count, float* min,
                                      float* max, short* has_draw_buffer,
                                      bool addshaders,
                                      bool embedTransparencyInfo) {
  auto G = I->G;
  bool ok = true;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  short ambient_occlusion = 0;
  float* sumarray = nullptr;
  int n_data = 0;

  if (embedTransparencyInfo) {
    int n_tri = count.num_total_indexes / 3;
    int bytes_to_allocate =
        2 * count.num_total_indexes *
            sizeof(VertexIndex_t) +  // vertexIndicesOriginal, vertexIndices
        3 * count.num_total_indexes * sizeof(float) +  // 3 * for sum
        n_tri * sizeof(float) +
        2 * n_tri * sizeof(int) +
        256 * sizeof(int);  // z_value (float * n_tri), ix (n_tri * int),
                            // sort_mem ((n_tri + 256) * int)
    // round to 4 byte words for the length of the CGO
    n_data = bytes_to_allocate / 4 + (((bytes_to_allocate % 4) == 0) ? 0 : 1);
  }
  std::vector<VertexIndex_t> vertexIndices(count.num_total_indexes);

  unsigned mul =
      VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
  mul += CGOPackedNormals(G) ? 1 : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_vertices) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_vertices;
  nxtn = 3;
  if (CGOPackedNormals(I->G)) {
    normalValsC = (uchar*)normalVals;
    nxtn = 1;
  }

  auto* colorVals = normalVals + nxtn * count.num_total_vertices;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*)colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_vertices);
  auto* accessibilityVals = pickColorVals + 3 * count.num_total_vertices;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_ACCESSIBILITY:
        cgo->current_accessibility = *pc;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_ARRAYS: {
        auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
        switch (sp->mode) {
          case GL_TRIANGLE_FAN:
          case GL_TRIANGLE_STRIP:
          case GL_TRIANGLES:
            shouldCompress = true;
          default:
            break;
        }
        if (shouldCompress) {
          int cnt, nxtn = 3;
          float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0,
                *normalValsDA, *accessibilityValsDA;
          float *pickColorValsDA, *pickColorValsTMP, *accessibilityValsTMP;

          nxtVals = vertexValsDA = sp->floatdata;
          for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
            set_min_max(min, max, &vertexValsDA[cnt]);
          }
          for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
            vertexVals[pl + cnt] = vertexValsDA[cnt];
          }
          if (CGOPackedNormals(I->G)) {
            if (sp->arraybits & CGO_NORMAL_ARRAY) {
              nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                    CLIP_NORMAL_VALUE(normalValsDA[cnt]);
              }
            } else {
              uchar norm[3] = {CLIP_NORMAL_VALUE(cgo->normal[0]),
                               CLIP_NORMAL_VALUE(cgo->normal[1]),
                               CLIP_NORMAL_VALUE(cgo->normal[2])};
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                    norm[cnt % 3];
              }
            }
          } else {
            if (sp->arraybits & CGO_NORMAL_ARRAY) {
              nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalVals[pl + cnt] = normalValsDA[cnt];
              }
            } else {
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalVals[pl + cnt] = cgo->normal[cnt % 3];
              }
            }
          }
          nxtn = 3;
          if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
            if (sp->arraybits & CGO_COLOR_ARRAY) {
              nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 4; cnt += 4) {
                colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
                colorValsUC[plc + cnt + 1] =
                    CLIP_COLOR_VALUE(colorValsDA[cnt + 1]);
                colorValsUC[plc + cnt + 2] =
                    CLIP_COLOR_VALUE(colorValsDA[cnt + 2]);
                colorValsUC[plc + cnt + 3] =
                    CLIP_COLOR_VALUE(colorValsDA[cnt + 3]);
              }
              nxtn = VERTEX_COLOR_SIZE;
            } else {
              uchar col[4] = {CLIP_COLOR_VALUE(cgo->color[0]),
                              CLIP_COLOR_VALUE(cgo->color[1]),
                              CLIP_COLOR_VALUE(cgo->color[2]),
                              CLIP_COLOR_VALUE(cgo->alpha)};
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorValsUC[plc + cnt] = col[cnt % 4];
              }
            }
          } else {
            if (sp->arraybits & CGO_COLOR_ARRAY) {
              nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 4; cnt += 4) {
                colorVals[plc + cnt] = colorValsDA[cnt];
                colorVals[plc + cnt + 1] = colorValsDA[cnt + 1];
                colorVals[plc + cnt + 2] = colorValsDA[cnt + 2];
                colorVals[plc + cnt + 3] = colorValsDA[cnt + 3];
              }
              nxtn = VERTEX_COLOR_SIZE;
            } else {
              float col[4] = {cgo->color[0], cgo->color[1], cgo->color[2],
                              cgo->alpha};
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorVals[plc + cnt] = col[cnt % 4];
              }
            }
          }
          if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
            nxtVals = nxtVals + (nxtn * sp->nverts);
            pickColorValsDA = nxtVals + sp->nverts;
            pickColorValsTMP = pickColorVals + (vpl * 2);
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
              CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
            }
            nxtn = VERTEX_PICKCOLOR_SIZE;
          } else {
            pickColorValsTMP = pickColorVals + (vpl * 2);
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
              CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
            }
          }
          if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY) {
            if (!ambient_occlusion) {
              for (cnt = 0; cnt < vpl; cnt++) {
                accessibilityVals[cnt] = 1.f;
              }
            }
            ambient_occlusion = 1;
            nxtVals = nxtVals + (nxtn * sp->nverts);
            accessibilityValsDA = nxtVals;
            accessibilityValsTMP = accessibilityVals + vpl;
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              accessibilityValsTMP[cnt] = accessibilityValsDA[cnt];
            }
          } else {
            if (ambient_occlusion) {
              accessibilityValsTMP = accessibilityVals + vpl;
              for (cnt = 0; cnt < sp->nverts; cnt++) {
                accessibilityValsTMP[cnt] = 1.f;
              }
            }
          }
          switch (sp->mode) {
            case GL_TRIANGLES:
              for (cnt = 0; cnt < sp->nverts; cnt++) {
                vertexIndices[idxpl++] = vpl + cnt;
              }
              break;
            case GL_TRIANGLE_STRIP: {
              short flip = 0;
              for (cnt = 2; cnt < sp->nverts; cnt++) {
                vertexIndices[idxpl++] = vpl + cnt - (flip ? 0 : 2);
                vertexIndices[idxpl++] = vpl + cnt - 1;
                vertexIndices[idxpl++] = vpl + cnt - (flip ? 2 : 0);
                flip = !flip;
              }
            } break;
            case GL_TRIANGLE_FAN:
              for (cnt = 2; cnt < sp->nverts; cnt++) {
                vertexIndices[idxpl++] = vpl;
                vertexIndices[idxpl++] = vpl + cnt - 1;
                vertexIndices[idxpl++] = vpl + cnt;
              }
              break;
          }
          pl += 3 * sp->nverts;
          plc += 4 * sp->nverts;
          vpl += sp->nverts;
        }
      } break;
      default:
        break;
    }
    ok &= !I->G->Interrupt;
  }
  if (sumarray) {
    for (idxpl = 0; idxpl < count.num_total_indexes; idxpl += 3) {
      add3f(&vertexVals[3 * vertexIndices[idxpl]],
            &vertexVals[3 * vertexIndices[idxpl + 1]], sumarray);
      add3f(&vertexVals[3 * vertexIndices[idxpl + 2]], sumarray, sumarray);
      sumarray += 3;
    }
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(I->G);

    auto& vbo = CreateVertexBuffer(I->G);
    BufferDataDesc dataDesc{
        {//
         BufferDesc{"a_Vertex", VertexFormat::Float3,
                    sizeof(float) * count.num_total_vertices * 3, vertexVals},
         BufferDesc{"a_Normal", fmt.normalFormat,
                    count.num_total_vertices * fmt.normalSize, normalVals},
         BufferDesc{"a_Color", fmt.colorFormat,
                    count.num_total_vertices * fmt.colorSize, colorVals},
         BufferDesc{"a_Accessibility", VertexFormat::Float,
                    sizeof(float) * count.num_total_vertices,
                    accessibilityVals}}};
    SetVertexBufferData(I->G, vbo, std::move(dataDesc));

    auto& ibo = CreateIndexBuffer(I->G);
    SetIndexBufferData(I->G, ibo, vertexIndices);

    auto vboid = vbo.get_hash_id();
    auto iboid = ibo.get_hash_id();

    auto& pickvbo = CreateVertexBuffer(I->G, VertexBufferLayout::Sequential,
                                       MemoryUsageProperty::CpuToGpu);
    BufferDataDesc pickDataDesc{
        {//
         BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes},
         BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes}}};
    SetVertexBufferData(I->G, pickvbo, std::move(pickDataDesc));
    auto pickvboid = pickvbo.get_hash_id();

    if (ok) {
      float* newPickColorVals;
      int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
                   CGO_PICK_COLOR_ARRAY;
      if (ambient_occlusion) {
        arrays |= CGO_ACCESSIBILITY_ARRAY;
      }
      if (addshaders) CGOEnable(cgo, GL_DEFAULT_SHADER);
      newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(
          GL_TRIANGLES, arrays, count.num_total_indexes,
          count.num_total_vertices, vboid, iboid, n_data, pickvboid);
      if (embedTransparencyInfo) {
        int n_tri = count.num_total_indexes / 3;
        float* sumarray;
        float* sum = sumarray = newPickColorVals + count.num_total_vertices * 3;
        float* z_value = sum + (count.num_total_indexes * 3);
        int* ix = (int*)z_value + n_tri;
        int* sort_mem = ix + n_tri;
        auto vertexIndicesOriginalTI = (VertexIndex_t*)(sort_mem + n_tri + 256);

        for (idxpl = 0; idxpl < count.num_total_indexes; idxpl += 3) {
          add3f(&vertexVals[3 * vertexIndices[idxpl]],
                &vertexVals[3 * vertexIndices[idxpl + 1]], sumarray);
          add3f(&vertexVals[3 * vertexIndices[idxpl + 2]], sumarray, sumarray);
          sumarray += 3;
        }
        memcpy(vertexIndicesOriginalTI, vertexIndices.data(),
               sizeof(VertexIndex_t) * count.num_total_indexes);
      }

      if (addshaders && ok) ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      CHECKOK(ok, newPickColorVals);
      if (!newPickColorVals) {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
      if (ok)
        memcpy(newPickColorVals + count.num_total_vertices, pickColorVals,
               count.num_total_vertices * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      I->G->ShaderMgr->freeGPUBuffer(vboid);
      I->G->ShaderMgr->freeGPUBuffer(iboid);
    }
  }
  return ok;
}

static bool OptimizeLinesToVBOIndexed(const CGO* I, CGO* cgo,
                                      const CGOCount& count, float* min,
                                      float* max, short* has_draw_buffer,
                                      bool addshaders) {
  auto G = I->G;
  bool ok = true;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, sz;
  bool hasNormals = 0;

  hasNormals = !CGOHasAnyLineVerticesWithoutNormals(I);
  std::vector<VertexIndex_t> vertexIndexes(count.num_total_indexes_lines);

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  if (hasNormals) {
    mul += CGOPackedNormals(G) ? 1 : VERTEX_NORMAL_SIZE;
  }
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_vertices_lines) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* nxtVals = vertexVals + VERTEX_POS_SIZE * count.num_total_vertices_lines;
  float* normalVals = nullptr;
  if (hasNormals) {
    normalVals = nxtVals;
    if (CGOPackedNormals(I->G)) {
      normalValsC = (uchar*)normalVals;
      sz = 1;
    } else {
      sz = VERTEX_NORMAL_SIZE;
    }
    nxtVals += sz * count.num_total_vertices_lines;
  }

  auto* colorVals = nxtVals;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*)colorVals;
    sz = 1;
  } else {
    sz = VERTEX_COLOR_SIZE;
  }
  nxtVals += sz * count.num_total_vertices_lines;

  auto* pickColorVals = nxtVals;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        break;
      case CGO_ACCESSIBILITY:
        cgo->current_accessibility = *pc;
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_ARRAYS: {
        auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
        switch (sp->mode) {
          case GL_LINE_LOOP:
          case GL_LINE_STRIP:
          case GL_LINES:
            shouldCompress = true;
          default:
            break;
        }
        if (shouldCompress) {
          int cnt, nxtn = 3;
          float *vertexValsDA = 0, *nxtVals2 = 0, *colorValsDA = 0,
                *normalValsDA = 0;
          float *pickColorValsDA = 0, *pickColorValsTMP;

          nxtVals2 = vertexValsDA = sp->floatdata;
          for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
            set_min_max(min, max, &vertexValsDA[cnt]);
          }
          for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
            vertexVals[pl + cnt] = vertexValsDA[cnt];
          }
          if (normalVals) {
            if (sp->arraybits & CGO_NORMAL_ARRAY) {
              if (CGOPackedNormals(I->G)) {
                nxtVals2 = normalValsDA = nxtVals2 + (nxtn * sp->nverts);
                for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                  normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                      CLIP_NORMAL_VALUE(normalValsDA[cnt]);
                }
              } else {
                nxtVals2 = normalValsDA = nxtVals2 + (nxtn * sp->nverts);
                for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                  normalVals[VAR_FOR_NORMAL + cnt] = normalValsDA[cnt];
                }
              }
            }
          }
          if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
            if (sp->arraybits & CGO_COLOR_ARRAY) {
              nxtVals2 = colorValsDA = nxtVals2 + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
              }
              nxtn = 4;
            } else {
              uchar col[4] = {CLIP_COLOR_VALUE(cgo->color[0]),
                              CLIP_COLOR_VALUE(cgo->color[1]),
                              CLIP_COLOR_VALUE(cgo->color[2]),
                              CLIP_COLOR_VALUE(cgo->alpha)};
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorValsUC[plc + cnt] = col[cnt % 4];
              }
            }
          } else {
            if (sp->arraybits & CGO_COLOR_ARRAY) {
              nxtVals2 = colorValsDA = nxtVals2 + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorVals[plc + cnt] = colorValsDA[cnt];
              }
              nxtn = 4;
            } else {
              float col[4] = {cgo->color[0], cgo->color[1], cgo->color[2],
                              cgo->alpha};
              for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
                colorVals[plc + cnt] = col[cnt % 4];
              }
            }
          }
          if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
            nxtVals2 = nxtVals2 + (nxtn * sp->nverts);
            pickColorValsDA = nxtVals2 + sp->nverts;
            pickColorValsTMP = pickColorVals + (vpl * 2);
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
              CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
            }
            nxtn = VERTEX_PICKCOLOR_SIZE;
          } else {
            pickColorValsTMP = pickColorVals + (vpl * 2);
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
              CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
            }
          }
          if (idxpl + sp->nverts > count.num_total_indexes_lines) {
            PRINTFB(I->G, FB_CGO, FB_Errors)
            "ERROR: CGOOptimizeToVBOIndexed() num_total_indexes_lines=%d "
            "mode=%d nverts=%d idxpl=%d\n",
                count.num_total_indexes_lines, sp->mode, sp->nverts,
                idxpl ENDFB(I->G);
          }
          switch (sp->mode) {
            case GL_LINES:
              for (cnt = 0; cnt < sp->nverts; cnt++) {
                vertexIndexes[idxpl++] = vpl + cnt;
              }
              break;
            case GL_LINE_STRIP:
              for (cnt = 1; cnt < sp->nverts; cnt++) {
                vertexIndexes[idxpl++] = vpl + cnt - 1;
                vertexIndexes[idxpl++] = vpl + cnt;
              }
              break;
            case GL_LINE_LOOP:
              for (cnt = 1; cnt < sp->nverts; cnt++) {
                vertexIndexes[idxpl++] = vpl + cnt - 1;
                vertexIndexes[idxpl++] = vpl + cnt;
              }
              vertexIndexes[idxpl++] = vpl;
              vertexIndexes[idxpl++] = vpl + sp->nverts - 1;
              break;
          }

          pl += 3 * sp->nverts;
          plc += 4 * sp->nverts;
          vpl += sp->nverts;
        }
      } break;
      case CGO_SPECIAL:
        CGOSpecial(cgo, CGO_get_int(pc));
    }
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(I->G);

    auto& vbo = CreateVertexBuffer(I->G);
    BufferDataDesc dataDesc{{
        //
        BufferDesc{"a_Vertex", VertexFormat::Float3,
                   sizeof(float) * count.num_total_vertices_lines * 3,
                   vertexVals},
        BufferDesc{"a_Normal", fmt.normalFormat,
                   count.num_total_vertices_lines * fmt.normalSize, normalVals},
        BufferDesc{"a_Color", fmt.colorFormat,
                   count.num_total_vertices_lines * fmt.colorSize, colorVals},
    }};
    SetVertexBufferData(I->G, vbo, std::move(dataDesc));

    auto& ibo = CreateIndexBuffer(I->G);
    SetIndexBufferData(I->G, ibo, vertexIndexes);

    auto vboid = vbo.get_hash_id();
    auto iboid = ibo.get_hash_id();

    auto& pickvbo = CreateVertexBuffer(I->G, VertexBufferLayout::Sequential,
                                       MemoryUsageProperty::CpuToGpu);
    BufferDataDesc pickDataDesc{
        {//
         BufferDesc("a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes),
         BufferDesc("a_Color", VertexFormat::UByte4Norm,
                    sizeof(float) * count.num_total_indexes)}};
    SetVertexBufferData(I->G, pickvbo, std::move(pickDataDesc));
    auto pickvboid = pickvbo.get_hash_id();

    if (ok) {
      if (addshaders) {
        CGOEnable(cgo, GL_DEFAULT_SHADER);
        // TODO: Check if this next line is supposed to be in this
        // if-statement. VBONotIndexed has it otherwise.
        CGODisable(cgo, GL_SHADER_LIGHTING);
      }
      auto newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(
          GL_LINES,
          CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
              CGO_PICK_COLOR_ARRAY,
          count.num_total_indexes_lines, count.num_total_vertices_lines, vboid,
          iboid, 0, pickvboid);
      CHECKOK(ok, newPickColorVals);
      if (addshaders && ok) ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      if (!newPickColorVals) {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
      if (ok)
        memcpy(newPickColorVals + count.num_total_vertices_lines, pickColorVals,
               count.num_total_vertices_lines * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      I->G->ShaderMgr->freeGPUBuffer(vboid);
      I->G->ShaderMgr->freeGPUBuffer(iboid);
    }
  }
  return ok;
}

/**
 * Similar to CGOOptimizeToVBONotIndexed(), but works with sortable alpha
 * triangles.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param color Initial color to use until a CGO_COLOR operation is observed
 * @param addshaders Add Enable/Disable shader operations
 * @param embedTransparencyInfo ???
 *
 * @return New converted CGO with use_shader=true
 */
CGO* CGOOptimizeToVBOIndexed(const CGO* I, int est, const float* color,
                             bool addshaders, bool embedTransparencyInfo) {
  CGO* cgo;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  short has_draw_buffer = false;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX},
        max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int ok = true;

  auto count = CGOCountNumVertices(I);

  cgo = CGONewSized(I->G, I->c + est);
  CHECKOK(ok, cgo);
  if (ok) {
    if (color) {
      cgo->color[0] = color[0];
      cgo->color[1] = color[1];
      cgo->color[2] = color[2];
      cgo->alpha = color[3];
    } else {
      cgo->color[0] = 1.f;
      cgo->color[1] = 1.f;
      cgo->color[2] = 1.f;
      cgo->alpha = 1.f;
    }
  }

  if (count.num_total_vertices_points > 0) {
    /* This does not need to be indexed (for now) */
    if (!OptimizePointsToVBO(I, cgo, count.num_total_vertices_points, min, max,
                             &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }

  if (count.num_total_vertices > 0) {
    if (!OptimizeVertsToVBOIndexed(I, cgo, count, min, max, &has_draw_buffer,
                                   addshaders, embedTransparencyInfo)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok && count.num_total_vertices_lines > 0) {
    if (!OptimizeLinesToVBOIndexed(I, cgo, count, min, max, &has_draw_buffer,
                                   addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok &&
      (count.num_total_vertices > 0 || count.num_total_vertices_lines > 0)) {
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok) ok &= CGOStop(cgo);
  if (ok) {
    if (has_draw_buffer) {
      cgo->has_draw_buffers = true;
    }
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader) {
      cgo->cgo_shader_ub_color =
          SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
      cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

struct OptimizeSphereData {
  std::vector<float> vert;
  std::vector<unsigned char> color;
  std::vector<int> pickColor;
  std::vector<unsigned char> rightUpFlagsUB;
  std::vector<float> rightUpFlags;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int total_vert = 0;
  int total_spheres = 0;
};

/**
 * @param[in] I input CGO
 * @param[out] cgo output CGO
 * @param[in] num_total_spheres Number of spheres to optimize
 * @param[out] leftOverCGO CGO data not relevant to spheres.
 * @pre cgo must not be nullptr
 *
 * Note: To use leftOverCGO, it must already be preallocated.
 */

static OptimizeSphereData GetOptimizeSphereData(const CGO* I, CGO*& cgo,
                                                int num_total_spheres,
                                                CGO* leftOverCGO) {
  OptimizeSphereData sphereData;
  const int tot = VerticesPerSphere() * 4 * num_total_spheres;
  int rightup_flags[6] = {0, 1, 3, 3, 2, 0};

  auto cgo_shader_ub_flags =
      SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);

  sphereData.vert.resize(tot);
  auto vertVals = sphereData.vert.data();

  sphereData.color.resize(tot);
  auto colorValsUB = sphereData.color.data();

  unsigned char* rightUpFlagValsUB = nullptr;
  float* rightUpFlagVals = nullptr;
  if (cgo_shader_ub_flags) {
    sphereData.rightUpFlagsUB.resize(VALUES_PER_IMPOSTER_SPACE_COORD *
                                     VerticesPerSphere() * num_total_spheres);
    rightUpFlagValsUB = sphereData.rightUpFlagsUB.data();
  } else {
    sphereData.rightUpFlags.resize(VALUES_PER_IMPOSTER_SPACE_COORD *
                                   VerticesPerSphere() * num_total_spheres);
    rightUpFlagVals = sphereData.rightUpFlags.data();
  }

  bool has_picking = CGOHasOperationsOfType(I, CGO_PICK_COLOR);
  int* pickcolorVals = nullptr;
  if (has_picking) {
    // atom/bond info for picking, 2 ints for each sphere
    sphereData.pickColor.resize(num_total_spheres * 2 * 4);
    pickcolorVals = sphereData.pickColor.data();
  }

  cgo->alpha = 1.f;
  float min_alpha = 1.f;
  bool copyNormalToLeftOver = false;
  bool copyColorToLeftOver = false;
  bool copyPickColorToLeftOver = false;
  bool copyAlphaToLeftOver = false;
  bool ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_NORMAL:
        cgo->normal[0] = *pc;
        cgo->normal[1] = *(pc + 1);
        cgo->normal[2] = *(pc + 2);
        copyNormalToLeftOver = true;
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
        copyColorToLeftOver = true;
        break;
      case CGO_ALPHA:
        cgo->alpha = *pc;
        if (cgo->alpha < min_alpha) min_alpha = cgo->alpha;
        copyAlphaToLeftOver = true;
        break;
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        copyPickColorToLeftOver = true;
        break;
      case CGO_SPHERE:
        for (int vv = 0; vv < VerticesPerSphere();
             vv++) {  // generate eight vertices of a bounding box for each
                      // cylinder
          vertVals[0] = *(pc);
          vertVals[1] = *(pc + 1);
          vertVals[2] = *(pc + 2);
          vertVals[3] = *(pc + 3);
          set_min_max(sphereData.min, sphereData.max, vertVals);
          if (cgo_shader_ub_flags) {
            rightUpFlagValsUB[0] = rightup_flags[vv];
            rightUpFlagValsUB++;
          } else {
            rightUpFlagVals[0] = rightup_flags[vv];
            rightUpFlagVals++;
          }
          colorValsUB[0] = CLIP_COLOR_VALUE(cgo->color[0]);
          colorValsUB[1] = CLIP_COLOR_VALUE(cgo->color[1]);
          colorValsUB[2] = CLIP_COLOR_VALUE(cgo->color[2]);
          colorValsUB[3] = CLIP_COLOR_VALUE(cgo->alpha);
          colorValsUB += 4;
          vertVals += 4;
          sphereData.total_vert++;
        }
        if (has_picking) {
          *(pickcolorVals++) = cgo->current_pick_color_index;
          *(pickcolorVals++) = cgo->current_pick_color_bond;
        }
        sphereData.total_spheres++;
        break;
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeSpheresToVBONonIndexed() CGO_DRAW_BUFFERS_INDEXED "
        "or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeCylindersToVBO() "
        "CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS encountered op=0x%X\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_LABELS:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_LABELS encountered "
        "op=0x%X\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_TEXTURES:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_TEXTURES encountered "
        "op=0x%X\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_ARRAYS:
      default:
        if (!leftOverCGO) break;
        if (copyAlphaToLeftOver) {
          copyAlphaToLeftOver = false;
          CGOAlpha(leftOverCGO, cgo->alpha);
        }
        if (copyColorToLeftOver) {
          copyColorToLeftOver = false;
          CGOColor(leftOverCGO, cgo->color[0], cgo->color[1], cgo->color[2]);
        }
        if (copyNormalToLeftOver) {
          CGONormalv(leftOverCGO, cgo->normal);
        }
        if (copyPickColorToLeftOver) {
          copyPickColorToLeftOver = false;
          CGOPickColor(leftOverCGO, cgo->current_pick_color_index,
                       cgo->current_pick_color_bond);
        }
        leftOverCGO->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  return sphereData;
}

/**
 * Creates and fills OpenGL buffers with optimized sphere Data
 * @param[in] I input CGO
 * @param[out] cgo output optimized CGO
 * @param[in] addshaders adds shader call to CGO
 * @param[in] num_total_spheres number of total spheres to optimize
 * @param[in] sphereData Optimized Sphere Data
 * @pre cgo must not be nullptr
 */

static bool PopulateBufferOptimizedSphereData(
    const CGO* I, CGO* cgo, bool addshaders, int num_total_spheres,
    const OptimizeSphereData& sphereData) {
  bool ok = true;
  auto rtp = VertexFormat::Float;
  short rsz = sizeof(float);
  const void* radiusptr = sphereData.rightUpFlags.data();
  auto cgo_shader_ub_flags =
      SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);
  if (cgo_shader_ub_flags) {
    rtp = VertexFormat::UByte;
    rsz = sizeof(std::uint8_t);
    radiusptr = sphereData.rightUpFlagsUB.data();
  }

  auto& vbo = CreateVertexBuffer(I->G);
  BufferDataDesc vertDescs{{
      //
      BufferDesc("a_vertex_radius", VertexFormat::Float4,
                 sizeof(float) * sphereData.total_vert * 4,
                 sphereData.vert.data()),
      BufferDesc("a_Color", VertexFormat::UByte4Norm,
                 sizeof(float) * sphereData.total_vert,
                 sphereData.color.data()),
      BufferDesc("a_rightUpFlags", rtp,
                 rsz * sphereData.total_vert * VALUES_PER_IMPOSTER_SPACE_COORD,
                 radiusptr),
  }};
  SetVertexBufferData(I->G, vbo, std::move(vertDescs));
  auto vboid = vbo.get_hash_id();

  auto& pickvbo = CreateVertexBuffer(I->G, VertexBufferLayout::Sequential,
                                     MemoryUsageProperty::CpuToGpu);
  BufferDataDesc pickDataDesc{
      {//
       BufferDesc("a_Color", VertexFormat::UByte4Norm, 0),
       BufferDesc("a_Color", VertexFormat::UByte4Norm,
                  sizeof(float) * sphereData.total_vert)}};
  auto size = sizeof(float) * sphereData.total_vert * 2;
  SetVertexBufferDataSized(I->G, pickvbo, std::move(pickDataDesc), size);

  auto pickvboid = pickvbo.get_hash_id();

  cgo->has_draw_buffers = true;
  cgo->has_draw_sphere_buffers = true;

  auto freebuffers = [vboid, pickvboid, I]() {
    I->G->ShaderMgr->freeGPUBuffer(vboid);
    I->G->ShaderMgr->freeGPUBuffer(pickvboid);
  };
  if (ok) {
    if (addshaders) CGOEnable(cgo, GL_SPHERE_SHADER);
    auto pickcolor_data = (int*)cgo->add<cgo::draw::sphere_buffers>(
        sphereData.total_spheres, (cgo_shader_ub_flags ? 3 : 1), vboid,
        pickvboid);  // always cgo_shader_ub_color
    CHECKOK(ok, pickcolor_data);
    if (ok && !sphereData.pickColor.empty()) {
      memcpy(pickcolor_data, sphereData.pickColor.data(),
             num_total_spheres * 2 * 4);
    }
    if (ok && addshaders) ok &= CGODisable(cgo, GL_SPHERE_SHADER);
    if (!ok) {
      freebuffers();
    }
  } else {
    freebuffers();
  }
  return ok;
}

/**
 * Optimizes Sphere CGO and populates data into GPU buffers
 * @param[in] I input CGO
 * @param[in] est estimated size of output CGO
 * @param[in] addshaders adds shader call to CGO
 * @param[out] leftOverCGO CGO data not relevant to spheres.
 * @return optimized CGO
 *
 * Note: If estimated size is unknown beforehand, provide 0 to est.
 * Note: To use leftOverCGO, it must already be preallocated.
 */

CGO* CGOOptimizeSpheresToVBONonIndexed(const CGO* I, int est, bool addshaders,
                                       CGO* leftOverCGO) {
  bool ok = true;
  int num_total_spheres = CGOCountNumberOfOperationsOfType(I, CGO_SPHERE);

  if (num_total_spheres <= 0) {
    return nullptr;
  }
  auto cgo = CGONewSized(I->G, I->c + est);
  auto sphereData =
      GetOptimizeSphereData(I, cgo, num_total_spheres, leftOverCGO);

  if (ok && sphereData.total_spheres > 0) {
    ok = PopulateBufferOptimizedSphereData(I, cgo, addshaders,
                                           num_total_spheres, sphereData);
  }

  if (ok && num_total_spheres > 0) {
    ok &= CGOBoundingBox(cgo, sphereData.min, sphereData.max);
  }

  if (ok) ok &= CGOStop(cgo);

  if (ok) {
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader) {
      cgo->cgo_shader_ub_color = true;
      cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

CGO* CGOOptimizeBezier(const CGO* I) {
  auto cgo = std::make_unique<CGO>(I->G);
  int num_splines = CGOCountNumberOfOperationsOfType(I, CGO_BEZIER);
  auto& vbo = CreateVertexBuffer(I->G);
  std::vector<float> vertData;
  vertData.reserve(num_splines * CGO_BEZIER_SZ);
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_BEZIER: {
        vertData.resize(vertData.size() + CGO_BEZIER_SZ);
        std::copy_n(pc, CGO_BEZIER_SZ, vertData.end() - CGO_BEZIER_SZ);
      } break;
    }
  }

  std::size_t numDimensions = 3;
  // Upload ALL bezier patches (4 control points each), not just the first —
  // the original hardcoded one patch (4 verts), truncating multi-segment
  // splines. The Metal tube path draws every patch; the GL path still draws
  // the first via glDrawArrays(GL_PATCHES, 0, 4).
  std::size_t numVerts = 4 * static_cast<std::size_t>(num_splines);
  BufferDataDesc dataDesc{{
      BufferDesc("position", VertexFormat::Float3,
                 sizeof(float) * numVerts * numDimensions, vertData.data()),
  }};
  SetVertexBufferData(I->G, vbo, std::move(dataDesc));
  auto vboid = vbo.get_hash_id();

  CGOEnable(cgo.get(), GL_BEZIER_SHADER);
  cgo->add<cgo::draw::bezier_buffers>(vboid);
  CGODisable(cgo.get(), GL_BEZIER_SHADER);
  cgo->use_shader = true;
  return cgo.release();
}

CGO* CGOOptimizeTextures(const CGO* I, int est) {
  CGO* cgo = nullptr;
  int ok = true;
  auto num_total_textures =
      CGOCountNumberOfOperationsOfType(I, CGO_DRAW_TEXTURE);
  //  printf("CGOOptimizeTextures: num_total_textures=%d\n",
  //  num_total_textures);
  if (num_total_textures) {
    int place3 = 0, place2 = 0;
    std::vector<float> worldPos(num_total_textures * 18);
    std::vector<float> screenValues(num_total_textures * 18);
    std::vector<float> textExtents(num_total_textures * 12);
    std::vector<float> pickColorValsVec(num_total_textures *
                                        12); /* pick index and bond */
    auto* pickColorVals = pickColorValsVec.data();

    cgo = CGONewSized(I->G, 0);

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto op = it.op_code();
      const auto pc = it.data();

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
          copy3f(pc, &worldPos[place3]);
          copy3f(pc, &worldPos[place3 + 3]);
          copy3f(pc, &worldPos[place3 + 6]);
          copy3f(pc, &worldPos[place3 + 9]);
          copy3f(pc, &worldPos[place3 + 12]);
          copy3f(pc, &worldPos[place3 + 15]);
          copy3f(pc + 3, screenMin);
          copy3f(pc + 6, screenMax);
          copy4f(pc + 9, textExtent);
          copy3f(screenMin, &screenValues[place3]);
          copy3f(screenMin, &screenValues[place3 + 3]);
          copy3f(screenMin, &screenValues[place3 + 6]);
          copy3f(screenMin, &screenValues[place3 + 9]);
          copy3f(screenMin, &screenValues[place3 + 12]);
          copy3f(screenMax, &screenValues[place3 + 15]);
          screenValues[place3 + 4] = screenMax[1];
          screenValues[place3 + 6] = screenMax[0];
          screenValues[place3 + 10] = screenMax[1];
          screenValues[place3 + 12] = screenMax[0];
          screenValues[place3 + 17] = screenMin[2];
          place3 += 18;
          CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[0];
          textExtents[place2++] = textExtent[1];
          CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[0];
          textExtents[place2++] = textExtent[3];
          CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[2];
          textExtents[place2++] = textExtent[1];
          CGO_put_int(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_uint(pickColorVals + place2 + 1,
                       cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[0];
          textExtents[place2++] = textExtent[3];
          CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[2];
          textExtents[place2++] = textExtent[1];
          CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
          CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
          textExtents[place2++] = textExtent[2];
          textExtents[place2++] = textExtent[3];
        } break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      auto& vbo = CreateVertexBuffer(I->G, VertexBufferLayout::Sequential);
      BufferDataDesc dataDesc{
          {BufferDesc("attr_worldpos", VertexFormat::Float3,
                      sizeof(float) * num_total_textures * 18, worldPos.data()),
           BufferDesc("attr_screenoffset", VertexFormat::Float3,
                      sizeof(float) * num_total_textures * 18,
                      screenValues.data()),
           BufferDesc("attr_texcoords", VertexFormat::Float3,
                      sizeof(float) * num_total_textures * 18,
                      textExtents.data())}};
      SetVertexBufferData(I->G, vbo, std::move(dataDesc));
      auto vboid = vbo.get_hash_id();

      if (ok) {
        float* pickArray =
            cgo->add<cgo::draw::textures>(num_total_textures, vboid);
        CHECKOK(ok, pickArray);
        if (!pickArray) I->G->ShaderMgr->freeGPUBuffer(vboid);
        if (ok)
          memcpy(pickArray + num_total_textures * 6, pickColorVals,
                 num_total_textures * 12 * sizeof(float));
        if (ok) ok &= CGOStop(cgo);
      } else {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
      }
      if (!ok) {
        CGOFree(cgo);
      }
    }
  }
  return cgo;
}

CGO* CGOOptimizeConnectors(const CGO* I, int est) {
  CGO* cgo = nullptr;
  int num_total_connectors;
  int ok = true;
  // The Metal renderer expands connectors in its vertex shader, so it wants
  // the same compact one-record-per-connector layout the GL geometry-shader
  // path uses (and it never reads a_isCenterPt).
  int use_geometry_shaders =
      I->G->Renderer || SettingGetGlobal_b(I->G, cSetting_use_geometry_shaders);
  int factor = (use_geometry_shaders ? 1 : 4);
  num_total_connectors =
      CGOCountNumberOfOperationsOfType(I, CGO_DRAW_CONNECTOR);

  if (num_total_connectors) {
    uchar* isCenterPt = nullptr;
    int place3 = 0, place2 = 0, place = 0;
    std::vector<float> targetPt3d(num_total_connectors * 20 * factor,
                                  0); /* too much, relativeMode only needs 1
                                      byte per vertex, instead of 1 float */
    auto* labelCenterPt3d =
        targetPt3d.data() + (num_total_connectors * 3 * factor);
    auto* indentFactor = labelCenterPt3d + (num_total_connectors * 3 * factor);
    auto* screenWorldOffset =
        indentFactor + (num_total_connectors * 2 * factor);
    auto* connectorColor =
        screenWorldOffset + (num_total_connectors * 3 * factor);
    auto* textSize = connectorColor + (num_total_connectors * factor);
    auto* relativeMode =
        (uchar*)(textSize + (num_total_connectors * 2 * factor));
    auto* drawBkgrd = (uchar*)(relativeMode + (num_total_connectors * factor));
    auto* bkgrdColor = (float*)(drawBkgrd + (num_total_connectors * factor));
    auto* relExtLength = (float*)(bkgrdColor + (num_total_connectors * factor));
    auto* connectorWidth =
        (float*)(relExtLength + (num_total_connectors * factor));
    if (!use_geometry_shaders)
      isCenterPt = (uchar*)(connectorWidth + (num_total_connectors * factor));
    else
      isCenterPt = nullptr;
    cgo = CGONewSized(I->G, 0);

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto op = it.op_code();
      const auto pc = it.data();

      switch (op) {
        case CGO_PICK_COLOR:
          cgo->current_pick_color_index = CGO_get_uint(pc);
          cgo->current_pick_color_bond = CGO_get_int(pc + 1);
          break;
        case CGO_DRAW_BUFFERS_INDEXED:
        case CGO_DRAW_BUFFERS_NOT_INDEXED:
          PRINTFB(I->G, FB_CGO, FB_Warnings)
          "WARNING: CGOOptimizeConnectors() CGO_DRAW_BUFFERS_INDEXED or "
          "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
              op ENDFB(I->G);
          break;
        case CGO_DRAW_CONNECTOR: {
          uchar* uc;
          int f;
          if (!use_geometry_shaders) {
            isCenterPt[place] = 0;
            isCenterPt[place + 1] = 2;
            isCenterPt[place + 2] = 2;
            isCenterPt[place + 3] = 1;
          }
          copy3f(pc, &targetPt3d[place3]);
          copy3f(pc + 3, &labelCenterPt3d[place3]);
          copy2f(pc + 6, &indentFactor[place2]);
          copy3f(pc + 9, &screenWorldOffset[place3]);
          copy2f(pc + 12, &textSize[place2]);
          relativeMode[place] = (uchar)(int)pc[17];
          drawBkgrd[place] = (uchar)(int)pc[18];
          uc = (uchar*)&bkgrdColor[place];
          uc[0] = CLIP_COLOR_VALUE(*(pc + 19));
          uc[1] = CLIP_COLOR_VALUE(*(pc + 20));
          uc[2] = CLIP_COLOR_VALUE(*(pc + 21));
          uc[3] = CLIP_COLOR_VALUE(*(pc + 22));
          uc = (uchar*)&connectorColor[place];
          uc[0] = CLIP_COLOR_VALUE(*(pc + 14));
          uc[1] = CLIP_COLOR_VALUE(*(pc + 15));
          uc[2] = CLIP_COLOR_VALUE(*(pc + 16));
          uc[3] = 255;
          relExtLength[place] = *(pc + 8);
          connectorWidth[place] = *(pc + 23);
          place3 += 3;
          place2 += 2;
          place += 1;
          for (f = 1; f < factor; f++) {
            copy3f(pc, &targetPt3d[place3]);
            copy3f(pc + 3, &labelCenterPt3d[place3]);
            copy2f(pc + 6, &indentFactor[place2]);
            copy3f(pc + 9, &screenWorldOffset[place3]);
            copy2f(pc + 12, &textSize[place2]);
            relativeMode[place] = (uchar)(int)pc[17];
            drawBkgrd[place] = (uchar)(int)pc[18];
            relExtLength[place] = *(pc + 8);
            connectorWidth[place] = *(pc + 23);
            uc = (uchar*)&bkgrdColor[place];
            uc[0] = uc[-4];
            uc[1] = uc[-3];
            uc[2] = uc[-2];
            uc[3] = uc[-1];
            uc = (uchar*)&connectorColor[place];
            uc[0] = uc[-4];
            uc[1] = uc[-3];
            uc[2] = uc[-2];
            uc[3] = uc[-1];
            place3 += 3;
            place2 += 2;
            place += 1;
          }
        } break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      const size_t quant = factor * num_total_connectors;
      auto& vbo = CreateVertexBuffer(I->G);
      BufferDataDesc bufferDataDesc{
          {//
           BufferDesc("a_target_pt3d", VertexFormat::Float3,
                      sizeof(float) * 3 * quant, targetPt3d.data()),
           BufferDesc("a_center_pt3d", VertexFormat::Float3,
                      sizeof(float) * 3 * quant, labelCenterPt3d),
           BufferDesc("a_indentFactor", VertexFormat::Float2,
                      sizeof(float) * 2 * quant, indentFactor),
           BufferDesc("a_screenWorldOffset", VertexFormat::Float3,
                      sizeof(float) * 3 * quant, screenWorldOffset),
           BufferDesc("a_textSize", VertexFormat::Float2,
                      sizeof(float) * 2 * quant, textSize),
           BufferDesc("a_Color", VertexFormat::UByte4Norm,
                      sizeof(float) * quant, connectorColor),
           BufferDesc("a_relative_mode", VertexFormat::UByte,
                      sizeof(std::uint8_t) * quant, relativeMode),
           BufferDesc("a_draw_flags", VertexFormat::UByte,
                      sizeof(std::uint8_t) * quant, drawBkgrd),
           BufferDesc("a_bkgrd_color", VertexFormat::UByte4Norm,
                      sizeof(float) * quant, bkgrdColor),
           BufferDesc("a_rel_ext_length", VertexFormat::Float,
                      sizeof(float) * quant, relExtLength),
           BufferDesc("a_con_width", VertexFormat::Float, sizeof(float) * quant,
                      connectorWidth),
           BufferDesc("a_isCenterPt", VertexFormat::UByte,
                      sizeof(std::uint8_t) * quant, isCenterPt)}};
      SetVertexBufferData(I->G, vbo, std::move(bufferDataDesc));
      auto vboid = vbo.get_hash_id();
      if (ok) {
        cgo->add<cgo::draw::connectors>(num_total_connectors, vboid);
        if (ok) ok &= CGOStop(cgo);
      }
      if (!ok) {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        CGOFree(cgo);
      }
    }
  }
  CheckGLErrorOK(I->G, "ERROR: CGOOptimizeConnectors() end returns err=%d\n");
  return cgo;
}

CGO* CGOOptimizeScreenTexturesAndPolygons(CGO* I, int est) {
  auto G = I->G;
  int ok = true;

  int num_total_indices = CGOCountNumVerticesForScreen(I);
  if (num_total_indices <= 0) return nullptr;

  std::unique_ptr<CGO> cgo_managed(new CGO(G));
  auto* const cgo = cgo_managed.get();

  {
    float *colorVals = 0, *texcoordVals;
    int tot, nxtn;
    uchar* colorValsUC = 0;
    CGOAlpha(cgo, 1.f);
    cgo->alpha = 1.f;
    cgo->color[0] = 1.f;
    cgo->color[1] = 1.f;
    cgo->color[2] = 1.f;

    {
      int mul =
          6;  // 3 - screenoffset/vertex, 2 - texture coordinates, 1 - color
      /*
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
        mul++;
      } else {
        mul += 4;
        }*/
      tot = num_total_indices * mul;
    }

    auto vertexVals_managed = std::vector<float>(tot);
    float* vertexVals = vertexVals_managed.data();
    texcoordVals = vertexVals + 3 * num_total_indices;
    nxtn = 2;
    colorVals = texcoordVals + nxtn * num_total_indices;
    colorValsUC = (uchar*)colorVals;
    nxtn = 1;
    ok = CGOProcessScreenCGOtoArrays(G, I, vertexVals, texcoordVals, colorVals,
                                     colorValsUC);
    RETURN_VAL_IF_FAIL(ok && !G->Interrupt, nullptr);
    if (ok) {
      auto& vbo = CreateVertexBuffer(I->G);
      BufferDataDesc dataDesc{
          {BufferDesc("attr_screenoffset", VertexFormat::Float3,
                      sizeof(float) * num_total_indices * 3, vertexVals),
           BufferDesc("attr_texcoords", VertexFormat::Float2,
                      sizeof(float) * num_total_indices * 2, texcoordVals),
           BufferDesc("attr_backgroundcolor", VertexFormat::UByte4Norm,
                      sizeof(std::uint8_t) * num_total_indices * 4,
                      colorValsUC)}};
      SetVertexBufferData(I->G, vbo, std::move(dataDesc));
      auto vboid = vbo.get_hash_id();
      if (ok) {
        CGOEnable(cgo, GL_SCREEN_SHADER);
        cgo->add<cgo::draw::screen_textures>(num_total_indices, vboid);
        if (ok) ok &= CGODisable(cgo, GL_SCREEN_SHADER);
        RETURN_VAL_IF_FAIL(ok, nullptr);
      } else {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
      }
    }
    cgo->use_shader = true;
  }

  return cgo_managed.release();
}
