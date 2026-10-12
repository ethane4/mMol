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

#include "pymol/graphics/cgo.h"
#include "pymol/graphics/cgogl.h"
#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/print_utils.h"
#include "pymol/scene/scene.h"
#include "pymol/geometry/sphere.h"
#include "pymol/geometry/triangle.h"
#include "pymol/text/v_font.h"

#include "pymol/graphics/cgo_internal.h"

/* Whether CGO VBOs pack normals into bytes (cgo_shader_ub_normal). A GL-path
 * option ONLY: the Metal renderer (G->Renderer set) binds a_Normal as Float3
 * for the raster and the RT per-vertex normals alike, and CGOGL asserts it.
 * Every place that SIZES, FILLS or DECLARES the normal slot must ask this one
 * function -- asking the setting directly in any of them makes the declared
 * format and the data disagree. */
bool CGOPackedNormals(PyMOLGlobals* G) {
  return !G->Renderer && SettingGet<bool>(G, cSetting_cgo_shader_ub_normal);
}

static void set_current_pick_color(CGO* cgo, unsigned int idx, int bnd) {
  if (cgo) {
    cgo->current_pick_color_index = idx;
    cgo->current_pick_color_bond = bnd;
  }
}

/**
 * Assign a new pick color index for {context, index, bond} (only in the
 * first pass) and convert that to an RGBA color for the current picking pass.
 *
 * @param[out] cgo Set the "current pick color" in this CGO
 * @param[in,out] pickmgr Pick color manager
 * @param[out] color RGBA pick color to use in this pass
 * @param[in] context Object state identifier
 * @param[in] index Primary index
 * @param[in] bond Secondary index
 */
void AssignNewPickColor(CGO* cgo, PickColorManager* pickmgr,
                        unsigned char* color, const PickContext* context,
                        unsigned int index, int bond) {
  set_current_pick_color(cgo, index, bond);
  pickmgr->colorNext(color, context, index, bond);
}

std::size_t CGO_sz[] = {
    CGO_NULL_SZ,                            //
    CGO_NULL_SZ,                            //
    CGO_BEGIN_SZ,                           //
    CGO_END_SZ,                             //
    CGO_VERTEX_SZ,                          //
    CGO_NORMAL_SZ,                          //
    fsizeof<cgo::draw::color>(),            //
    fsizeof<cgo::draw::sphere>(),           //
    CGO_TRIANGLE_SZ,                        //
    fsizeof<cgo::draw::cylinder>(),         //
    CGO_LINEWIDTH_SZ,                       //
    CGO_WIDTHSCALE_SZ,                      //
    CGO_ENABLE_SZ,                          //
    CGO_DISABLE_SZ,                         //
    fsizeof<cgo::draw::sausage>(),          //
    fsizeof<cgo::draw::custom_cylinder>(),  //
    CGO_DOTWIDTH_SZ,                        //
    CGO_ALPHA_TRIANGLE_SZ,                  //
    CGO_ELLIPSOID_SZ,                       //
    CGO_FONT_SZ,                            //
    CGO_FONT_SCALE_SZ,                      //
    CGO_FONT_VERTEX_SZ,                     //
    CGO_FONT_AXES_SZ,                       //
    CGO_CHAR_SZ,                            //
    CGO_INDENT_SZ,                          //
    CGO_ALPHA_SZ,                           //
    CGO_QUADRIC_SZ,                         //
    CGO_CONE_SZ,                            //
    fsizeof<cgo::draw::arrays>(),           //
    CGO_NULL_SZ,
    CGO_RESET_NORMAL_SZ,  //
    CGO_PICK_COLOR_SZ,    //
    CGO_NULL_SZ,          // CGO_DRAW_BUFFERS_SZ no longer used
    fsizeof<cgo::draw::buffers_indexed>(),                  //
    CGO_BOUNDING_BOX_SZ,                                    //
    fsizeof<cgo::draw::buffers_not_indexed>(),              //
    CGO_SPECIAL_SZ,                                         //
    fsizeof<cgo::draw::cylinder_buffers>(),                 //
    fsizeof<cgo::draw::shadercylinder>(),                   //
    fsizeof<cgo::draw::shadercylinder2ndcolor>(),           //
    fsizeof<cgo::draw::sphere_buffers>(),                   //
    CGO_ACCESSIBILITY_SZ,                                   //
    CGO_DRAW_TEXTURE_SZ,                                    //
    fsizeof<cgo::draw::textures>(),                         //
    fsizeof<cgo::draw::screen_textures>(),                  //
    CGO_TEX_COORD_SZ,                                       //
    fsizeof<cgo::draw::label>(),                            //
    fsizeof<cgo::draw::labels>(),                           //
    CGO_DRAW_CONNECTOR_SZ,                                  //
    fsizeof<cgo::draw::connectors>(),                       //
    CGO_DRAW_TRILINES_SZ,                                   //
    CGO_UNIFORM3F_SZ,                                       //
    CGO_SPECIAL_WITH_ARG_SZ,                                //
    fsizeof<cgo::draw::line>(),                             //
    fsizeof<cgo::draw::splitline>(),                        //
    fsizeof<cgo::draw::custom>(),                           //
    fsizeof<cgo::draw::vertex_attribute_3f>(),              //
    fsizeof<cgo::draw::vertex_attribute_4ub>(),             //
    fsizeof<cgo::draw::vertex_attribute_1f>(),              //
    fsizeof<cgo::draw::mask_attribute_if_picking>(),        //
    fsizeof<cgo::draw::bind_vbo_for_picking>(),             //
    CGO_VERTEX_BEGIN_LINE_STRIP_SZ,                         //
    CGO_INTERPOLATED_SZ,                                    //
    CGO_VERTEX_CROSS_SZ,                                    //
    fsizeof<cgo::draw::vertex_attribute_4ub_if_picking>(),  //
    fsizeof<cgo::draw::custom_cylinder_alpha>(),            //
    CGO_BEZIER_SZ,                                          //
    fsizeof<cgo::draw::bezier_buffers>(),                   //
    CGO_NULL_SZ                                             //
};

/**
 * Get the number of elements in `CGO_sz`
 */
size_t CGO_sz_size() { return sizeof(CGO_sz) / sizeof(*CGO_sz); }

static float* CGO_add(CGO* I, unsigned c);
static float* CGO_size(CGO* I, int sz);

/**
 * Inverse function of CGOArrayFromPyListInPlace
 *
 * I: (input) Primitive CGO (may contain CGO_DRAW_ARRAYS)
 *
 * Return: All-float Python list primitive CGO
 */
static PyObject* CGOArrayAsPyList(const CGO* I) {
  std::vector<float> flat;
  flat.reserve(I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto op = it.op_code();
    auto pc = it.data();
    auto sz = CGO_sz[op];

    flat.push_back(op);

    switch (op) {
      case CGO_BEGIN:
      case CGO_ENABLE:
      case CGO_DISABLE:
      case CGO_SPECIAL:
        // first member int
        flat.push_back(*reinterpret_cast<const int*>(pc));
        ++pc;
        --sz;
        break;
      case CGO_PICK_COLOR:
        assert(sz == 2);
        flat.push_back(*reinterpret_cast<const int*>(pc++));  // atom index
        flat.push_back(*reinterpret_cast<const int*>(pc++));  // bond index
        sz -= 2;
        break;
      case CGO_DRAW_ARRAYS: {
        auto sp = reinterpret_cast<const cgo::draw::arrays*>(pc);
        flat.push_back(sp->mode);
        flat.push_back(sp->arraybits);
        flat.push_back(sp->narrays);  // (redundant)
        flat.push_back(sp->nverts);
        pc = sp->get_data();
        sz = sp->get_data_length();
      }
    }

    // float members
    for (; sz; --sz) {
      flat.push_back(*(pc++));
    }
  }

  return PConvToPyObject(flat);
}

PyObject* CGOAsPyList(CGO* I) {
  PyObject* result;
  result = PyList_New(2);
  PyObject* list = CGOArrayAsPyList(I);
  PyList_SetItem(result, 0, PyInt_FromLong(PyList_Size(list)));
  PyList_SetItem(result, 1, list);
  return (result);
}

static float CPythonVal_PyFloat_AsDouble_From_List(void* G, PyObject* list,
                                                   size_t i) {
  float out;
  PConvPyFloatToFloat(PyList_GetItem(list, i), &out);
  return out;
}

/**
 * Inverse function of CGOArrayAsPyList
 *
 * list: (input) All-float Python list primitive CGO (may contain
 * CGO_DRAW_ARRAYS) I: (output) empty CGO
 */
static int CGOArrayFromPyListInPlace(PyObject* list, CGO* I) {
  // sanity check
  if (!list || !PyList_Check(list)) return false;

  auto G = I->G;

#define GET_FLOAT(i) \
  ((float)CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))
#define GET_INT(i) ((int)CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))

  for (int i = 0, l = PyList_Size(list); i < l;) {
    unsigned op = GET_INT(i++);
    ok_assert(1, op < CGO_sz_size());
    int sz = CGO_sz[op];
    float* fdata = I->add_to_buffer(sz + 1);
    CGO_write_int(fdata, op);

    switch (op) {
      case CGO_STOP:
        // don't increment size for null terminator
        I->c -= 1;
        return true;
      case CGO_BEGIN:
        I->has_begin_end = true;
      case CGO_ENABLE:
      case CGO_DISABLE:
      case CGO_SPECIAL:
        // first member int
        ok_assert(1, i < l);
        CGO_write_int(fdata, GET_INT(i++));
        sz--;
        break;
      case CGO_PICK_COLOR:
        ok_assert(1, i + 1 < l);
        CGO_write_int(fdata, GET_INT(i++));  // atom index
        CGO_write_int(fdata, GET_INT(i++));  // bond index
        sz -= 2;
        break;
      case CGO_DRAW_ARRAYS: {
        // has abstract superclass, need to be constructed!
        ok_assert(1, i + 3 < l);
        auto sp = new (fdata)
            cgo::draw::arrays(GET_INT(i), GET_INT(i + 1), GET_INT(i + 3));

        // sanity check
        int narrays_check = GET_INT(i + 2);
        if (sp->narrays != narrays_check) {
          PRINTFB(I->G, FB_CGO, FB_Warnings)
          " CGO-Warning: narrays mismatch: %d != %d\n", sp->narrays,
              narrays_check ENDFB(I->G);
        }

        // data
        sz = sp->get_data_length();
        sp->floatdata = fdata = I->allocate_in_data_heap(sz);

        i += 4;
      } break;
    }

    // float members
    for (; sz; --sz) {
      ok_assert(1, i < l);
      *(fdata++) = GET_FLOAT(i++);
    }
  }

#undef GET_FLOAT
#undef GET_INT

  return true;

ok_except1:
  PRINTFB(G, FB_CGO, FB_Errors) " %s-Error: Corrupt data\n", __func__ ENDFB(G);
  return false;
}

CGO* CGONewFromPyList(PyMOLGlobals* G, PyObject* list, int version,
                      bool shouldCombine) {
  int ok = true;
  auto I = CGONew(G);
  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if ((version > 0) && (version <= 86)) {
    if (ok) ok = PConvFromPyListItem(G, list, 0, I->c);
    if (ok) VLACheck(I->op, float, I->c);
    if (ok)
      ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 1), I->op, I->c);
  } else {
    if (ok) ok = CGOArrayFromPyListInPlace(PyList_GetItem(list, 1), I);
  }
  if (!ok) {
    CGOFree(I);
  }
  {
    CGO* cgo = nullptr;
    if (shouldCombine && I && I->has_begin_end) {
      cgo = CGOCombineBeginEnd(I, 0);
      CGOFree(I);
    } else {
      cgo = I;
    }
    return cgo;
  }
}

CGO::CGO(PyMOLGlobals* G, int size) : G(G) {
  op = VLACalloc(float, size + 32);
  cgo_shader_ub_color = SettingGet<bool>(G, cSetting_cgo_shader_ub_color);
  cgo_shader_ub_normal = CGOPackedNormals(G);
}

void CGOSetUseShader(CGO* I, int use_shader) {
  I->use_shader = use_shader;
  if (use_shader) {
    I->cgo_shader_ub_color =
        SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color);
    I->cgo_shader_ub_normal = CGOPackedNormals(I->G);
  } else {
    I->cgo_shader_ub_color = 0;
    I->cgo_shader_ub_normal = 0;
  }
}
void CGOReset(CGO* I) {
  I->c = 0;
  I->z_flag = false;
  I->alpha = 1.f;
  I->has_begin_end = false;
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f;
  I->normal[1] = 0.f;
  I->normal[2] = 1.f;
  I->color[0] = 0.f;
  I->color[1] = 0.f;
  I->color[2] = 1.f;
  I->pickColor[0] = 0;
  I->pickColor[1] = 0;
  I->pickColor[2] = 0;
  I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
}

void CGOFree(CGO*& I, bool withVBOs) {
  if (I) {
    if (!withVBOs) {
      I->has_draw_buffers = false;
    }
    DeleteP(I);
  }
}

CGO::~CGO() {
  if (has_draw_buffers) {
    CGOFreeVBOs(this);
  }
  FreeP(i_start);
  VLAFreeP(op);
}

static float* CGO_add(CGO* I, unsigned c) {
  float* at;
  VLACheck(I->op, float, I->c + c);
  if (!I->op) {
    return nullptr;
  }
  at = I->op + I->c;
  I->c += c;
  return (at);
}

static float* CGO_size(CGO* I, int sz) {
  float* at;
  VLASize(I->op, float, sz);
  if (!I->op) {
    return nullptr;
  }
  at = I->op + I->c;
  I->c = sz;
  return (at);
}

/*===== Object Creation Routines =======*/

int CGOFromFloatArray(CGO* I, const float* src, int len) {
  int iarg;
  int ok;
  int all_ok = true;
  int bad_entry = 0;
  int sz;
  int a;
  int cc = 0;
  float val;
  float *pc, *save_pc, *tf;
  VLACheck(I->op, float, I->c + len + 32);
  save_pc = I->op + I->c;
  while (len-- > 0) {
    cc++;
    const auto op = static_cast<unsigned>(*(src++));

    if (op >= CGO_sz_size()) {
      bad_entry = cc;
      break;
    }

    sz = CGO_sz[op];
    if (len < sz) break; /* discard short instruction */
    len -= sz;
    pc = save_pc;
    CGO_write_int(pc, op);
    ok = true;
    for (a = 0; a < sz; a++) {
      cc++;
      val = *(src++);
      if (std::abs(val) <= R_SMALL8) {
        val = 0;
      }
      if ((FLT_MAX - val) > 0.0F) { /* make sure we have a real float */
        *(pc++) = val;
      } else {
        *(pc++) = 0.0;
        ok = false;
      }
    }
    if (ok) {
      switch (op) {
        case CGO_END:
        case CGO_VERTEX:
        case CGO_BEGIN:
          I->has_begin_end = true;
      }
      switch (op) { /* now convert any instructions with int arguments */
        case CGO_BEGIN:
        case CGO_ENABLE:
        case CGO_DISABLE:
        case CGO_SPECIAL:
          tf = save_pc + 1;
          iarg = (int)*(tf);
          CGO_write_int(tf, iarg);
          break;
        case CGO_PICK_COLOR:
          tf = save_pc + 1;
          CGO_write_int(tf, int(save_pc[1]));  // atom index
          CGO_write_int(tf, int(save_pc[2]));  // bond index
          break;
      }
      save_pc = pc;
      I->c += sz + 1;
    } else { /* discard illegal instructions */
      if (all_ok) bad_entry = cc;
      all_ok = false;
    }
  }
  return (bad_entry);
}

int CGOBegin(CGO* I, int mode) {
  float* pc = CGO_add(I, CGO_BEGIN_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_BEGIN);
  CGO_write_int(pc, mode);
  I->has_begin_end = true;
  I->texture[0] = 0.f;
  I->texture[1] = 0.f;
  return true;
}

int CGOEnd(CGO* I) {
  float* pc = CGO_add(I, CGO_END_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_END);
  I->has_begin_end = true;
  return true;
}

int CGOEnable(CGO* I, int mode) {
  float* pc = CGO_add(I, CGO_ENABLE_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_ENABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGODisable(CGO* I, int mode) {
  float* pc = CGO_add(I, CGO_DISABLE_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_DISABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGOLinewidth(CGO* I, float v) {
  float* pc = CGO_add(I, CGO_LINEWIDTH_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_LINEWIDTH);
  *(pc++) = v;
  return true;
}

/**
 * implements special-case operations inside a CGO
 *
 * v: lookup value defined for each special operation (see CGO.h)
 */
int CGOSpecial(CGO* I, int v) {
  float* pc = CGO_add(I, CGO_SPECIAL_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_SPECIAL);
  CGO_write_int(pc, v);
  return true;
}

/**
 * implements special-case operations with an argument
 * inside a CGO
 *
 * v: lookup value defined for each special operation (see CGO.h)
 * argval : argument value
 */
int CGOSpecialWithArg(CGO* I, int v, float argval) {
  float* pc = CGO_add(I, CGO_SPECIAL_WITH_ARG_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_SPECIAL_WITH_ARG);
  CGO_write_int(pc, v);
  *pc = argval;
  return true;
}

int CGODotwidth(CGO* I, float v) {
  float* pc = CGO_add(I, CGO_DOTWIDTH_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_DOTWIDTH);
  *(pc++) = v;
  return true;
}

/* CGOUniform3f - specifies a 3f uniform variable and
   its value.  This function returns the offset of where
   these values are stored inside the CGO float array
   so that they can be accessed and changed from outside
   the CGO.

 */
int CGOUniform3f(CGO* I, int uniform_id, const float* value) {
  float* pc = CGO_add(I, CGO_UNIFORM3F_SZ + 1);
  if (!pc) return 0;
  CGO_write_int(pc, CGO_UNIFORM3F);
  CGO_write_int(pc, uniform_id);
  copy3f(value, pc);
  return pc - I->op;
}

int CGOBoundingBox(CGO* I, const float* min, const float* max) {
  float* pc = CGO_add(I, CGO_BOUNDING_BOX_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_BOUNDING_BOX);
  *(pc++) = *(min);
  *(pc++) = *(min + 1);
  *(pc++) = *(min + 2);
  *(pc++) = *(max);
  *(pc++) = *(max + 1);
  *(pc++) = *(max + 2);
  return true;
}

int CGOAccessibility(CGO* I, float a) {
  float* pc = CGO_add(I, CGO_ACCESSIBILITY_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_ACCESSIBILITY);
  *(pc++) = a;
  return true;
}

int CGODrawTexture(CGO* I, float* worldPos, float* screenMin, float* screenMax,
                   float* textExtent) {
  float* pc = CGO_add(I, CGO_DRAW_TEXTURE_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_DRAW_TEXTURE);
  *(pc++) = worldPos[0];
  *(pc++) = worldPos[1];
  *(pc++) = worldPos[2];
  *(pc++) = screenMin[0];
  *(pc++) = screenMin[1];
  *(pc++) = screenMin[2];
  *(pc++) = screenMax[0];
  *(pc++) = screenMax[1];
  *(pc++) = screenMax[2];
  *(pc++) = textExtent[0];
  *(pc++) = textExtent[1];
  *(pc++) = textExtent[2];
  *(pc++) = textExtent[3];
  return true;
}

int CGODrawConnector(CGO* I, float* targetPt3d, float* labelCenterPt3d,
                     float text_width, float text_height, float* indentFactor,
                     float* screenWorldOffset, float* connectorColor,
                     short relativeMode, int draw_flags, float bkgrd_transp,
                     float* bkgrd_color, float rel_ext_length,
                     float connectorWidth) {
  float* pc = CGO_add(I, CGO_DRAW_CONNECTOR_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_DRAW_CONNECTOR);
  *(pc++) = targetPt3d[0];
  *(pc++) = targetPt3d[1];
  *(pc++) = targetPt3d[2];
  *(pc++) = labelCenterPt3d[0];
  *(pc++) = labelCenterPt3d[1];
  *(pc++) = labelCenterPt3d[2];
  *(pc++) = indentFactor[0];
  *(pc++) = indentFactor[1];
  *(pc++) = rel_ext_length; /* place for ext_length relative to height (i.e.,
                               text_height which is total height */
  *(pc++) = screenWorldOffset[0];
  *(pc++) = screenWorldOffset[1];
  *(pc++) = screenWorldOffset[2];
  *(pc++) = text_width;
  *(pc++) = text_height;
  *(pc++) = connectorColor[0];
  *(pc++) = connectorColor[1];
  *(pc++) = connectorColor[2];
  *(pc++) = (float)relativeMode;
  *(pc++) = (float)draw_flags;
  *(pc++) = bkgrd_color[0];
  *(pc++) = bkgrd_color[1];
  *(pc++) = bkgrd_color[2];
  *(pc++) = bkgrd_transp;
  *(pc++) = connectorWidth;  // place for label_connector_width
  return true;
}

int CGOPickColor(CGO* I, unsigned int index, int bond) {
  // check if uchar is -1 since extrude does this for masked atoms
  if (index == (unsigned int)-1) {
    bond = cPickableNoPick;
  }
  if (I->current_pick_color_index == index &&
      I->current_pick_color_bond == bond)
    return true;
  float* pc = CGO_add(I, CGO_PICK_COLOR_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_PICK_COLOR);
  CGO_write_uint(pc, index);
  CGO_write_int(pc, bond);
  I->current_pick_color_index = index;
  I->current_pick_color_bond = bond;
  return true;
}

int CGOAlpha(CGO* I, float alpha) {
  float* pc = CGO_add(I, CGO_ALPHA_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_ALPHA);
  *(pc++) = alpha;
  I->alpha = alpha;
  return true;
}

int CGOSphere(CGO* I, const float* v1, float r) {
  float* pc = CGO_add(I, CGO_SPHERE_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_SPHERE);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
  return true;
}

int CGOEllipsoid(CGO* I, const float* v1, float r, const float* n1,
                 const float* n2, const float* n3) {
  float* pc = CGO_add(I, CGO_ELLIPSOID_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_ELLIPSOID);

  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
  *(pc++) = *(n1++);
  *(pc++) = *(n1++);
  *(pc++) = *(n1++);
  *(pc++) = *(n2++);
  *(pc++) = *(n2++);
  *(pc++) = *(n2++);
  *(pc++) = *(n3++);
  *(pc++) = *(n3++);
  *(pc++) = *(n3++);
  return true;
}

void CGOSetZVector(CGO* I, float z0, float z1, float z2) {
  I->z_flag = true;
  I->z_vector[0] = z0;
  I->z_vector[1] = z1;
  I->z_vector[2] = z2;
  I->z_min = FLT_MAX;
  I->z_max = -FLT_MAX;
}

const static float one_third = 1.0F / 3.0F;
const static float _0 = 0.0F;

int CGOAlphaTriangle(CGO* I, const float* v1, const float* v2, const float* v3,
                     const float* n1, const float* n2, const float* n3,
                     const float* c1, const float* c2, const float* c3,
                     float a1, float a2, float a3, int reverse) {
  if (v1 && v2 && v3) {
    float* pc = CGO_add(I, CGO_ALPHA_TRIANGLE_SZ + 1);
    float z = _0;
    if (!pc) return false;
    CGO_write_int(pc, CGO_ALPHA_TRIANGLE);
    CGO_write_int(pc, 0);  // this is the place for the next triangle in the bin
    *(pc++) = (v1[0] + v2[0] + v3[0]) * one_third;
    *(pc++) = (v1[1] + v2[1] + v3[1]) * one_third;
    *(pc++) = (v1[2] + v2[2] + v3[2]) * one_third;
    if (I->z_flag) {
      float* zv = I->z_vector;
      z = pc[-3] * zv[0] + pc[-2] * zv[1] + pc[-1] * zv[2];
      if (z > I->z_max) I->z_max = z;
      if (z < I->z_min) I->z_min = z;
    }
    *(pc++) = z;

    if (reverse) {
      *(pc++) = *(v2++); /* vertices @ +5 */
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
    } else {
      *(pc++) = *(v1++); /* vertices @ +5 */
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
    }

    *(pc++) = *(v3++);
    *(pc++) = *(v3++);
    *(pc++) = *(v3++);

    if (reverse) {
      *(pc++) = *(n2++); /* normals @ +14 */
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
    } else {
      *(pc++) = *(n1++); /* normals @ +14 */
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
    }
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);

    if (reverse) {
      *(pc++) = *(c2++); /* colors @ +23 */
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = a2;
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = a1;
    } else {
      *(pc++) = *(c1++); /* colors @ +23 */
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = a1;
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = a2;
    }
    *(pc++) = *(c3++);
    *(pc++) = *(c3++);
    *(pc++) = *(c3++);
    *(pc++) = a3;
  }
  return true;
}

int CGOVertex(CGO* I, float v1, float v2, float v3) {
  float* pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  return true;
}

int CGOVertexv(CGO* I, const float* v) {
  float* pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOVertexCrossv(CGO* I, const float* v) {
  float* pc = CGO_add(I, CGO_VERTEX_CROSS_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_VERTEX_CROSS);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOColor(CGO* I, float v1, float v2, float v3) {
  float* pc = CGO_add(I, CGO_COLOR_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_COLOR);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  I->color[0] = v1;
  I->color[1] = v2;
  I->color[2] = v3;
  return true;
}

int CGOColorv(CGO* I, const float* v) { return CGOColor(I, v[0], v[1], v[2]); }

int CGOTexCoord2f(CGO* I, float v1, float v2) {
  float* pc = CGO_add(I, CGO_TEX_COORD_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_TEX_COORD);
  *(pc++) = v1;
  *(pc++) = v2;
  I->texture[0] = v1;
  I->texture[1] = v2;
  return true;
}

int CGONormal(CGO* I, float v1, float v2, float v3) {
  float* pc = CGO_add(I, CGO_NORMAL_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  I->normal[0] = v1;
  I->normal[1] = v2;
  I->normal[2] = v3;
  return true;
}

int CGOResetNormal(CGO* I, int mode) {
  float* pc = CGO_add(I, CGO_RESET_NORMAL_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_RESET_NORMAL);
  CGO_write_int(pc, mode);
  SceneGetResetNormal(I->G, I->normal, mode);
  return true;
}

int CGOFontScale(CGO* I, float v1, float v2) {
  float* pc = CGO_add(I, CGO_FONT_SCALE_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_FONT_SCALE);
  *(pc++) = v1;
  *(pc++) = v2;
  return true;
}

int CGONormalv(CGO* I, const float* v) {
  float* pc = CGO_add(I, CGO_NORMAL_SZ + 1);
  if (!pc) return false;
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

/**
 * Add a null terminator to the CGO buffer, but don't increment
 * the size variable (CGO::c).
 */
int CGOStop(CGO* I) {
#define CGO_STOP_ZEROS 1

  float* pc = CGO_size(I, I->c + CGO_STOP_ZEROS);
  if (!pc) return false;
  UtilZeroMem(pc, sizeof(float) * CGO_STOP_ZEROS);
  I->c -= CGO_STOP_ZEROS;
  return true;
}

int CGOCheckComplex(CGO* I) {
  int fc = 0;
  const SphereRec* sp = I->G->Sphere->Sphere[1];

  /* stick_quality needs to match *every* CGO? */
  auto nEdge = SettingGet<int>(I->G, cSetting_stick_quality);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
      case CGO_CYLINDER:
      case CGO_CONE:
      case CGO_SAUSAGE:
      case CGO_CUSTOM_CYLINDER:
      case CGO_CUSTOM_CYLINDER_ALPHA:
        fc += 3 * (3 + (nEdge + 1) * 9) + 9;
        break;
      case CGO_ELLIPSOID:
      case CGO_QUADRIC:
      case CGO_SPHERE:
        fc += (sp->NVertTot * 6) + (sp->NStrip * 3) + 3;
        break;
      case CGO_DRAW_ARRAYS: {
        cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        fc += sp->nverts;
      } break;
      case CGO_DRAW_BUFFERS_INDEXED: {
        cgo::draw::buffers_indexed* sp = reinterpret_cast<decltype(sp)>(pc);
        switch (sp->mode) {
          case GL_TRIANGLES:
            fc += sp->nindices / 3;
            break;
          case GL_LINES:
            fc += sp->nindices / 2;
            break;
        }
      } break;
      case CGO_DRAW_BUFFERS_NOT_INDEXED: {
        cgo::draw::buffers_not_indexed* sp = reinterpret_cast<decltype(sp)>(pc);
        switch (sp->mode) {
          case GL_TRIANGLES:
            fc += sp->nverts / 3;
            break;
          case GL_LINES:
            fc += sp->nverts / 2;
            break;
        }
      } break;
      case CGO_DRAW_SPHERE_BUFFERS: {
        cgo::draw::sphere_buffers* sp = reinterpret_cast<decltype(sp)>(pc);
        fc += sp->num_spheres * VerticesPerSphere();
      } break;
      case CGO_DRAW_CYLINDER_BUFFERS: {
        cgo::draw::cylinder_buffers* sp = reinterpret_cast<decltype(sp)>(pc);
        fc += sp->num_cyl * NumVerticesPerCylinder();
      } break;
    }
  }
  return (fc);
}

int CGOPreloadFonts(CGO* I) {
  int ok = true;
  int font_seen = false;
  int font_id;

  auto blocked = PAutoBlock(I->G);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
      case CGO_FONT:
        ok = ok && (VFontLoad(I->G, 1.0, 1, 1, true));
        font_seen = true;
        break;
      case CGO_CHAR:
        if (!font_seen) {
          font_id = VFontLoad(I->G, 1.0, 1, 1, true);
          ok = ok && font_id;
          font_seen = true;
        }
        break;
    }
  }
  if (blocked) PUnblock(I->G);
  return (ok);
}

int CGOCheckForText(CGO* I) {
  int fc = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
      case CGO_FONT:
      case CGO_FONT_AXES:
      case CGO_FONT_SCALE:
        fc++;
        break;
      case CGO_INDENT:
      case CGO_FONT_VERTEX:
        fc++;
        break;
      case CGO_CHAR:
        fc += 3 + 2 * 3 * 10; /* est 10 lines per char */
        break;
    }
  }
  PRINTFD(I->G, FB_CGO)
  " CGOCheckForText-Debug: %d\n", fc ENDFD;

  return (fc);
}

CGO* CGODrawText(const CGO* I, int est,
                 float* camera) { /* assumes blocked intepreter */
  CGO* cgo;
  int font_id = 0;
  char text[2] = " ";
  float pos[] = {0.0F, 0.0F, 0.0F};
  float axes[] = {1.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1.0F};
  float scale[2] = {1.0, 1.0};

  cgo = CGONewSized(I->G, I->c + est);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_FONT:
        break;
      case CGO_FONT_AXES:
        break;
      case CGO_FONT_SCALE:
        scale[0] = pc[0];
        scale[1] = pc[1];
        break;
      case CGO_FONT_VERTEX:
        copy3f(pc, pos);
        break;
      case CGO_INDENT:
        text[0] = (unsigned char)*pc;
        VFontIndent(I->G, font_id, text, pos, scale, axes, pc[1]);
        break;
      case CGO_CHAR:
        if (!font_id) {
          font_id = VFontLoad(I->G, 1.0, 1, 1, false);
        }
        text[0] = (unsigned char)*pc;
        VFontWriteToCGO(I->G, font_id, cgo, text, pos, scale, axes, cgo->color);
        break;
      case CGO_COLOR:
        cgo->color[0] = *pc;
        cgo->color[1] = *(pc + 1);
        cgo->color[2] = *(pc + 2);
      default:
        cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  if (cgo && cgo->has_begin_end) {
    /* this is mainly for VFontWriteToCGO() that still creates CGOBegin/CGOEnd
     */
    if (cgo && cgo->has_begin_end) {
      CGO* convertcgo = nullptr;
      convertcgo = CGOCombineBeginEnd(cgo, 0);
      CGOFree(cgo);
      cgo = convertcgo;
    }
  }
  return (cgo);
}

/**
 * Release all shader resources from this CGO
 */
void CGOFreeVBOs(CGO* I) {
  constexpr bool freevbos = true;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
      case CGO_DRAW_TRILINES: {
        unsigned buf = it.cast<cgo::draw::trilines>()->buffer;
        if (freevbos) I->G->ShaderMgr->AddVBOToFree(buf);
      } break;
      case CGO_DRAW_CUSTOM: {
        auto sp = it.cast<cgo::draw::custom>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
        }
      } break;
      case CGO_DRAW_SPHERE_BUFFERS: {
        auto sp = it.cast<cgo::draw::sphere_buffers>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
        }
      } break;
      case CGO_DRAW_LABELS: {
        auto sp = it.cast<cgo::draw::labels>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
        }
      } break;
      case CGO_DRAW_TEXTURES: {
        auto sp = it.cast<cgo::draw::textures>();
        if (freevbos) I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
      } break;
      case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS: {
        auto sp = it.cast<cgo::draw::screen_textures>();
        if (freevbos) I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
      } break;
      case CGO_DRAW_CYLINDER_BUFFERS: {
        auto sp = it.cast<cgo::draw::cylinder_buffers>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
        }
      } break;
      case CGO_DRAW_BUFFERS_NOT_INDEXED: {
        auto sp = it.cast<cgo::draw::buffers_not_indexed>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
          I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
        }
      } break;
      case CGO_DRAW_BUFFERS_INDEXED: {
        auto sp = it.cast<cgo::draw::buffers_indexed>();
        if (freevbos) {
          I->G->ShaderMgr->freeGPUBuffers(
              {sp->vboid, sp->iboid, sp->pickvboid});
        }
      } break;
      case CGO_DRAW_CONNECTORS: {
        auto sp = it.cast<cgo::draw::connectors>();
        if (freevbos) I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
      } break;
    }
  }
}

/* ======== Raytrace Renderer ======== */

int CGOGetExtent(const CGO* I, float* mn, float* mx) {
  int result = false;

#define check_extent(v, r)                                        \
  {                                                               \
    if (!result) {                                                \
      mn[0] = ((*(v)) - r);                                       \
      mx[0] = ((*(v)) + r);                                       \
      mn[1] = ((*(v + 1)) - r);                                   \
      mx[1] = ((*(v + 1)) + r);                                   \
      mn[2] = ((*(v + 2)) - r);                                   \
      mx[2] = ((*(v + 2)) + r);                                   \
      result = true;                                              \
    } else {                                                      \
      if (mn[0] > ((*(v)) - r)) mn[0] = ((*(v)) - r);             \
      if (mx[0] < ((*(v)) + r)) mx[0] = ((*(v)) + r);             \
      if (mn[1] > ((*((v) + 1)) - r)) mn[1] = ((*((v) + 1)) - r); \
      if (mx[1] < ((*((v) + 1)) + r)) mx[1] = ((*((v) + 1)) + r); \
      if (mn[2] > ((*((v) + 2)) - r)) mn[2] = ((*((v) + 2)) - r); \
      if (mx[2] < ((*((v) + 2)) + r)) mx[2] = ((*((v) + 2)) + r); \
    }                                                             \
  }

#define check_extent4(v, r)                                       \
  {                                                               \
    if (!result) {                                                \
      mn[0] = ((*(v)) - r);                                       \
      mx[0] = ((*(v)) + r);                                       \
      mn[1] = ((*(v + 1)) - r);                                   \
      mx[1] = ((*(v + 1)) + r);                                   \
      mn[2] = ((*(v + 2)) - r);                                   \
      mx[2] = ((*(v + 2)) + r);                                   \
      mn[3] = ((*(v + 3)) - r);                                   \
      mx[3] = ((*(v + 3)) + r);                                   \
      result = true;                                              \
    } else {                                                      \
      if (mn[0] > ((*(v)) - r)) mn[0] = ((*(v)) - r);             \
      if (mx[0] < ((*(v)) + r)) mx[0] = ((*(v)) + r);             \
      if (mn[1] > ((*((v) + 1)) - r)) mn[1] = ((*((v) + 1)) - r); \
      if (mx[1] < ((*((v) + 1)) + r)) mx[1] = ((*((v) + 1)) + r); \
      if (mn[2] > ((*((v) + 2)) - r)) mn[2] = ((*((v) + 2)) - r); \
      if (mx[2] < ((*((v) + 2)) + r)) mx[2] = ((*((v) + 2)) + r); \
      if (mn[3] > ((*((v) + 3)) - r)) mn[3] = ((*((v) + 3)) - r); \
      if (mx[3] < ((*((v) + 3)) + r)) mx[3] = ((*((v) + 3)) + r); \
    }                                                             \
  }

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_VERTEX:
        check_extent(pc, 0);
        break;
      case CGO_SPHERE:
      case CGO_ELLIPSOID:
        check_extent(pc, *(pc + 3));
        break;
      case CGO_CYLINDER:
      case CGO_CONE:
      case CGO_SAUSAGE:
      case CGO_CUSTOM_CYLINDER:
      case CGO_CUSTOM_CYLINDER_ALPHA:
        check_extent(pc, *(pc + 6));
        check_extent(pc + 3, *(pc + 6));
        break;
      case CGO_TRIANGLE:
        check_extent(pc, 0);
        check_extent(pc + 3, 0);
        check_extent(pc + 6, 0);
        break;
      case CGO_DRAW_ARRAYS: {
        const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        const float* pct = sp->floatdata;
        int pl;

        if (sp->arraybits & CGO_VERTEX_ARRAY) {
          for (pl = 0; pl < sp->nverts; pl++) {
            check_extent(pct, 0);
            pct += VERTEX_POS_SIZE;
          }
        }
      } break;
      case CGO_BOUNDING_BOX: {
        if (!result) {
          mn[0] = (*pc);
          mn[1] = *(pc + 1);
          mn[2] = *(pc + 2);
          mx[0] = *(pc + 3);
          mx[1] = *(pc + 4);
          mx[2] = *(pc + 5);
          result = true;
        } else {
          if (mn[0] > *pc) mn[0] = (*pc);
          if (mn[1] > *(pc + 1)) mn[1] = *(pc + 1);
          if (mn[2] > *(pc + 2)) mn[2] = *(pc + 2);
          if (mx[0] < *(pc + 3)) mx[0] = *(pc + 3);
          if (mx[1] < *(pc + 4)) mx[1] = *(pc + 4);
          if (mx[2] < *(pc + 5)) mx[2] = *(pc + 5);
        }
      }
    }
  }
  return (result);
}

int CGOHasNormals(const CGO* I) {
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
      case CGO_NORMAL:
      case CGO_SPHERE:
      case CGO_ELLIPSOID:
      case CGO_CYLINDER:
      case CGO_CONE:
      case CGO_SAUSAGE:
      case CGO_CUSTOM_CYLINDER:
      case CGO_CUSTOM_CYLINDER_ALPHA:
        return true;
      case CGO_DRAW_ARRAYS:
        if (it.cast<cgo::draw::arrays>()->arraybits & CGO_NORMAL_ARRAY) {
          return true;
        }
        break;
    }
  }
  return false;
}

bool CGO::append(const CGO& source_, bool stopAtEnd) {
  const CGO* const source = &source_;
  int ok = 1;

  for (auto it = source->begin(); !it.is_stop(); ++it) {
    add_to_cgo(it.op_code(), it.data());
  }

  if (stopAtEnd) ok &= CGOStop(this);
  has_draw_buffers |= source->has_draw_buffers;
  has_draw_cylinder_buffers |= source->has_draw_cylinder_buffers;
  return ok;
}

/**
 * Appends `src` to the end of this CGO. Takes ownership of data
 * (incl. VBOs) and leaves `src` as a valid but empty CGO.
 */
void CGO::move_append(CGO&& src_) {
  CGO* const src = &src_;
  if (!src->c) return;

  // copy buffer
  VLACheck(op, float, c + src->c);
  memcpy(op + c, src->op, src->c * sizeof(float));

  // update sizes
  c += src->c;
  src->c = 0;

  // null terminators (CGO_STOP)
  *(op + c) = 0;
  *(src->op) = 0;

  // move heap data
  for (auto& ref : src->_data_heap) {
    _data_heap.emplace_back(std::move(ref));
  }
  src->_data_heap.clear();

  // copy boolean flags
  has_draw_buffers |= src->has_draw_buffers;
  has_draw_cylinder_buffers |= src->has_draw_cylinder_buffers;
  has_draw_sphere_buffers |= src->has_draw_sphere_buffers;
  has_begin_end |= src->has_begin_end;
  use_shader |= src->use_shader;
  render_alpha |= src->render_alpha;
  src->has_draw_buffers = false;
}

/**
 * Appends `src` to the end of this CGO and then free's `src`
 * and sets the pointer to nullptr.
 */
void CGO::free_append(CGO*& src) {
  free_append(std::move(src));
  assert(src == nullptr);
}
void CGO::free_append(CGO*&& src) {
  if (!src) return;
  move_append(std::move(*src));
  DeleteP(src);
}

int CGOAppend(CGO* dest, const CGO* source, bool stopAtEnd) {
  return dest->append(*source, stopAtEnd);
}

int CGOCountNumberOfOperationsOfType(const CGO* I, int optype) {
  std::set<int> ops = {optype};
  return CGOCountNumberOfOperationsOfTypeN(I, ops);
}

int CGOCountNumberOfOperationsOfTypeN(const CGO* I,
                                      const std::set<int>& optype) {
  int numops = 0;
  for (auto cgoit = I->begin(); !cgoit.is_stop(); ++cgoit) {
    if (optype.count(cgoit.op_code())) numops++;
  }
  return (numops);
}

int CGOCountNumberOfOperationsOfTypeN(const CGO* I,
                                      const std::map<int, int>& optype) {
  int numops = 0;
  for (auto cgoit = I->begin(); !cgoit.is_stop(); ++cgoit) {
    auto it = optype.find(cgoit.op_code());
    if (it != optype.end()) numops += it->second;
  }
  return (numops);
}

bool CGOHasOperationsOfType(const CGO* I, int optype) {
  std::set<int> ops = {optype};
  return CGOHasOperationsOfTypeN(I, ops);
}

bool CGOHasOperations(const CGO* I) { return !I->begin().is_stop(); }

bool CGOHasOperationsOfTypeN(const CGO* I, const std::set<int>& optype) {
  if (!I->op) return false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (optype.count(it.op_code())) return 1;
  }
  return (0);
}

static bool CGOFilterOutOperationsOfTypeN(const CGO* I, CGO* cgo,
                                          const std::set<int>& optype) {
  if (!I->op) return false;

  bool ret = false;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto op = it.op_code();
    if (optype.find(op) == optype.end()) {
      auto pc = it.data();
      cgo->add_to_cgo(op, pc);
    } else {
      ret = true;  // returns if filtered anything
    }
  }
  return ret;
}

bool CGOFilterOutCylinderOperationsInto(const CGO* I, CGO* cgo) {
  static std::set<int> optypes = {
      CGO_SHADER_CYLINDER, CGO_SHADER_CYLINDER_WITH_2ND_COLOR,
      CGO_SAUSAGE,         CGO_CYLINDER,
      CGO_CUSTOM_CYLINDER, CGO_CUSTOM_CYLINDER_ALPHA};
  return CGOFilterOutOperationsOfTypeN(I, cgo, optypes);
}

bool CGOFilterOutBezierOperationsInto(const CGO* I, CGO* cgo) {
  static std::set<int> optypes = {CGO_BEZIER};
  return CGOFilterOutOperationsOfTypeN(I, cgo, optypes);
}

bool CGOHasCylinderOperations(const CGO* I) {
  static std::set<int> optypes = {
      CGO_SHADER_CYLINDER, CGO_SHADER_CYLINDER_WITH_2ND_COLOR,
      CGO_SAUSAGE,         CGO_CYLINDER,
      CGO_CUSTOM_CYLINDER, CGO_CUSTOM_CYLINDER_ALPHA};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOHasSphereOperations(const CGO* I) {
  static std::set<int> optypes = {CGO_SPHERE};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOHasBezierOperations(const CGO* I) {
  static std::set<int> optypes = {CGO_BEZIER};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOCheckWhetherToFree(PyMOLGlobals* G, CGO* I) {
  if (I->use_shader) {
    if (I->cgo_shader_ub_color !=
            SettingGetGlobal_i(G, cSetting_cgo_shader_ub_color) ||
        I->cgo_shader_ub_normal != CGOPackedNormals(G)) {
      return true;
    }
  }
  return false;
}

CGO* CGOColorByRamp(PyMOLGlobals* G, const CGO* I, ObjectGadgetRamp* ramp,
                    int state, CSetting* set1) {
  if (!I) {
    return nullptr;
  }

  auto cgo = CGONewSized(G, 0);
  if (!cgo) {
    return nullptr;
  }

  int ok = true;
  float white[3] = {1.f, 1.f, 1.f};
  float probe_radius = SettingGet_f(G, set1, nullptr, cSetting_solvent_radius);
  float v_above[3], n0[3] = {0.f, 0.f, 0.f};
  int ramp_above =
      SettingGet_i(G, set1, nullptr, cSetting_surface_ramp_above_mode) == 1;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    bool skipCopy = false;
    switch (op) {
      case CGO_DRAW_ARRAYS: {
        const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        float* vals =
            cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
        int nvals = sp->narrays * sp->nverts;
        ok &= vals ? true : false;
        if (ok) memcpy(vals, sp->floatdata, nvals);
        skipCopy = true;
      } break;
      case CGO_NORMAL:
        copy3f(pc, n0);
        break;
      case CGO_VERTEX: {
        float color[3];
        copy3f(white, color);
        if (ramp_above) {
          copy3f(n0, v_above);
          scale3f(v_above, probe_radius, v_above);
          add3f(pc, v_above, v_above);
        } else {
          copy3f(pc, v_above);
        }
        if (ObjectGadgetRampInterVertex(ramp, v_above, color, state)) {
          CGOColorv(cgo, color);
        } else {
          CGOColorv(cgo, white);
        }
      } break;
    }
    if (!skipCopy) {
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

/**
 * FIXME: This function always returns true for `checkOpaque=ture`
 */
int CGOHasTransparency(const CGO* I, bool checkTransp, bool checkOpaque) {
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_ALPHA) {
      const auto pc = it.data();
      if (checkTransp && *pc < 1.f) return 1;
      if (checkOpaque && *pc == 1.f) return 1;
    }
  }

  return checkOpaque;
}

CGO* CGOConvertTrianglesToAlpha(const CGO* I) {
  int tot_nverts = 0;

  CGO* cgo = CGONewSized(I->G, I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
      case CGO_DRAW_ARRAYS: {
        const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        const int mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
        const float* nxtVals = sp->floatdata;

        assert(arrays & CGO_VERTEX_ARRAY);
        const float* vertexValsDA = nxtVals;
        nxtVals += VERTEX_POS_SIZE * nverts;

        const float* normalValsDA = nullptr;
        if (arrays & CGO_NORMAL_ARRAY) {
          normalValsDA = nxtVals;
          nxtVals += VERTEX_NORMAL_SIZE * nverts;
        }

        const float* colorValsDA = nullptr;
        if (arrays & CGO_COLOR_ARRAY) {
          colorValsDA = nxtVals;
          nxtVals += VERTEX_COLOR_SIZE * nverts;
        }

        const float* const vertexVals0 = vertexValsDA;
        const float* const normalVals0 = normalValsDA;
        const float* const colorVals0 = colorValsDA;

        switch (mode) {
          case GL_TRIANGLES: {
            for (int cnt = 0; cnt < nverts; cnt += 3) {
              if (colorVals0) {
                CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA + 3,
                                 vertexValsDA + 6, normalValsDA,
                                 normalValsDA + 3, normalValsDA + 6,
                                 colorValsDA, colorValsDA + 4, colorValsDA + 8,
                                 *(colorValsDA + 3), *(colorValsDA + 7),
                                 *(colorValsDA + 11), 0);
              } else {
                CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA + 3,
                                 vertexValsDA + 6, normalValsDA,
                                 normalValsDA + 3, normalValsDA + 6, cgo->color,
                                 cgo->color, cgo->color, cgo->alpha, cgo->alpha,
                                 cgo->alpha, 0);
              }
              vertexValsDA += 9;
              normalValsDA += 9;
              if (colorVals0) colorValsDA += 12;
            }
          }
            tot_nverts += nverts;
            break;
          case GL_TRIANGLE_STRIP: {
            short flip = 0;
            vertexValsDA += 6;
            normalValsDA += 6;
            if (colorVals0) colorValsDA += 8;

            for (int cnt = 2; cnt < nverts; cnt++) {
              if (colorVals0) {
                CGOAlphaTriangle(cgo, vertexValsDA - 6, vertexValsDA - 3,
                                 vertexValsDA, normalValsDA - 6,
                                 normalValsDA - 3, normalValsDA,
                                 colorValsDA - 8, colorValsDA - 4, colorValsDA,
                                 *(colorValsDA - 5), *(colorValsDA - 1),
                                 *(colorValsDA + 3), flip);
              } else {
                CGOAlphaTriangle(cgo, vertexValsDA - 6, vertexValsDA - 3,
                                 vertexValsDA, normalValsDA - 6,
                                 normalValsDA - 3, normalValsDA, cgo->color,
                                 cgo->color, cgo->color, cgo->alpha, cgo->alpha,
                                 cgo->alpha, flip);
              }
              vertexValsDA += 3;
              normalValsDA += 3;
              if (colorVals0) colorValsDA += 4;
              flip = !flip;
            }
          }
            tot_nverts += nverts;
            break;
          case GL_TRIANGLE_FAN: {
            vertexValsDA += 6;
            normalValsDA += 6;
            if (colorVals0) colorValsDA += 8;
            for (int cnt = 2; cnt < nverts; cnt++) {
              if (colorVals0) {
                CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA - 3,
                                 vertexValsDA, normalVals0, normalValsDA - 3,
                                 normalValsDA, colorVals0, colorValsDA - 4,
                                 colorValsDA, *(colorVals0 + 3),
                                 *(colorValsDA - 1), *(colorValsDA + 3), 0);
              } else {
                CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA - 3,
                                 vertexValsDA, normalVals0, normalValsDA - 3,
                                 normalValsDA, cgo->color, cgo->color,
                                 cgo->color, cgo->alpha, cgo->alpha, cgo->alpha,
                                 0);
              }
              vertexValsDA += 3;
              normalValsDA += 3;
              if (colorVals0) colorValsDA += 4;
            }
          }
            tot_nverts += nverts;
            break;
        }
      } break;
      case CGO_END:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOConvertTrianglesToAlpha: CGO_END encountered without CGO_BEGIN "
        "but "
        "skipped for OpenGLES\n" ENDFB(I->G);
        break;
      case CGO_VERTEX:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOConvertTrianglesToAlpha: CGO_VERTEX encountered without CGO_BEGIN "
        "but skipped for OpenGLES\n" ENDFB(I->G);
        break;
      case CGO_BEGIN: {
        float vertices[3][3], colors[4][3], normals[4][3], alpha[4];
        short verticespl = 2, colorspl = 2, normalspl = 2, alphapl = 2;
        short hasShifted = 0;
        int nverts = 0, err = 0;
        int mode = CGO_get_int(pc);
        short mode_is_triangles = 0, flip = 0, mode_is_fan = 0;
        copy3f(cgo->color, colors[3]);
        copy3f(cgo->normal, normals[3]);
        alpha[3] = cgo->alpha;
        switch (mode) {
          case GL_TRIANGLE_FAN:
            mode_is_fan = 1;
          case GL_TRIANGLES:
          case GL_TRIANGLE_STRIP:
            mode_is_triangles = 1;
        }
        if (!mode_is_triangles) {
          CGOBegin(cgo, mode);
        }

        for (++it; !err && it != CGO_END; ++it) {
          auto pc = it.data();
          int op = it.op_code();
          short add_to_cgo = 1;
          switch (op) {
            case CGO_VERTEX:
              if (mode_is_triangles) {
                if (!(hasShifted & 1)) {  // colors
                  if (colorspl >= 0) {
                    copy3f(colors[colorspl + 1], colors[colorspl]);
                    colorspl--;
                  } else {
                    if (!mode_is_fan) copy3f(colors[1], colors[2]);
                    copy3f(colors[0], colors[1]);
                  }
                }
                if (!(hasShifted & 2)) {  // normals
                  if (normalspl >= 0) {
                    copy3f(normals[normalspl + 1], normals[normalspl]);
                    normalspl--;
                  } else {
                    if (!mode_is_fan) copy3f(normals[1], normals[2]);
                    copy3f(normals[0], normals[1]);
                  }
                }
                if (!(hasShifted & 4)) {  // alphas
                  if (alphapl >= 0) {
                    alpha[alphapl] = alpha[alphapl + 1];
                    alphapl--;
                  } else {
                    if (!mode_is_fan) alpha[2] = alpha[1];
                    alpha[1] = alpha[0];
                  }
                }
                if (verticespl >= 0) {
                  copy3f(pc, vertices[verticespl]);
                  verticespl--;
                } else {
                  if (!mode_is_fan) copy3f(vertices[1], vertices[2]);
                  copy3f(vertices[0], vertices[1]);
                  copy3f(pc, vertices[0]);
                }

                nverts++;
                switch (mode) {
                  case GL_TRIANGLES:
                    if (!(nverts % 3)) {
                      CGOAlphaTriangle(cgo, vertices[2], vertices[1],
                                       vertices[0], normals[2], normals[1],
                                       normals[0], colors[2], colors[1],
                                       colors[0], alpha[2], alpha[1], alpha[0],
                                       0);
                    }
                    break;
                  case GL_TRIANGLE_STRIP:
                    if (verticespl < 0) {
                      int off0, off2;
                      if (flip) {
                        off0 = 0;
                        off2 = 2;
                      } else {
                        off0 = 2;
                        off2 = 0;
                      }
                      flip = !flip;
                      CGOAlphaTriangle(cgo, vertices[off0], vertices[1],
                                       vertices[off2], normals[off0],
                                       normals[1], normals[off2], colors[off0],
                                       colors[1], colors[off2], alpha[off0],
                                       alpha[1], alpha[off2], 0);
                    }
                    break;
                  case GL_TRIANGLE_FAN:
                    if (verticespl < 0) {
                      CGOAlphaTriangle(cgo, vertices[2], vertices[1],
                                       vertices[0], normals[2], normals[1],
                                       normals[0], colors[2], colors[1],
                                       colors[0], alpha[2], alpha[1], alpha[0],
                                       0);
                    }
                }
                add_to_cgo = !mode_is_triangles;
                hasShifted = 0;
              } else {
                add_to_cgo = 1;
              }
            case CGO_COLOR:
              if (op == CGO_COLOR) {
                add_to_cgo = !mode_is_triangles;
                if (mode_is_triangles) {
                  if (colorspl >= 0) {
                    copy3f(pc, colors[colorspl]);
                    colorspl--;
                  } else {
                    if (!mode_is_fan) copy3f(colors[1], colors[2]);
                    copy3f(colors[0], colors[1]);
                    copy3f(pc, colors[0]);
                  }
                  hasShifted |= 1;
                }
              }
            case CGO_NORMAL:
              if (op == CGO_NORMAL) {
                add_to_cgo = !mode_is_triangles;
                if (mode_is_triangles) {
                  if (normalspl >= 0) {
                    copy3f(pc, normals[normalspl]);
                    normalspl--;
                  } else {
                    if (!mode_is_fan) copy3f(normals[1], normals[2]);
                    copy3f(normals[0], normals[1]);
                    copy3f(pc, normals[0]);
                  }
                  hasShifted |= 2;
                }
              }
            case CGO_ALPHA:
              if (op == CGO_ALPHA) {
                add_to_cgo = !mode_is_triangles;
                if (mode_is_triangles) {
                  if (alphapl >= 0)
                    alpha[alphapl--] = *pc;
                  else {
                    if (!mode_is_fan) alpha[2] = alpha[1];
                    alpha[1] = alpha[0];
                    alpha[0] = *pc;
                  }
                  hasShifted |= 4;
                }
              }
            default:
              if (add_to_cgo) {
                cgo->add_to_cgo(op, pc);
              }
          }
        }

        if (!mode_is_triangles) {
          CGOEnd(cgo);
        }

        tot_nverts += nverts;
      } break;
      case CGO_COLOR:
        if (op == CGO_COLOR) {
          copy3f(pc, cgo->color);
        }
      case CGO_NORMAL:
        if (op == CGO_NORMAL) {
          copy3f(pc, cgo->normal);
        }
      case CGO_ALPHA:
        if (op == CGO_ALPHA) {
          cgo->alpha = *pc;
        }
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
  if (tot_nverts) {
    return (cgo);
  } else {
    CGOFree(cgo);
    return nullptr;
  }
}

/**
 * Converts TRIANGLE(_STRIP|_FAN) to TRIANGLES and generates
 * normals for all triangles. Discards any existing normals for
 * triangles.
 *
 * I: primitive CGO
 * return: new primitive CGO with normals on triangles
 */
CGO* CGOGenerateNormalsForTriangles(const CGO* I) {
  auto G = I->G;
  auto cgo = CGONewSized(G, I->c);

  float vertices[3][3];
  float current_color[3] = {0.f, 0.f, 0.f}, colors[3][3];
  float current_normal[3];
  float current_alpha = 0, alphas[3];

  bool has_alpha = false;
  bool has_color = false;

  int mode = 0;
  bool inside_begin_triangles = false;
  int current_i = 0;
  int vertex_count = 0;
  bool flip = false;
  bool emit;

  const int indices_regular[] = {0, 1, 2};
  const int indices_flipped[] = {0, 2, 1};

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    auto op = it.op_code();

    if (op == CGO_BEGIN) {
      mode = *reinterpret_cast<const int*>(pc);

      switch (mode) {
        case GL_TRIANGLE_STRIP:
        case GL_TRIANGLE_FAN:
        case GL_TRIANGLES:
          current_i = 0;
          vertex_count = 0;
          flip = false;
          inside_begin_triangles = true;

          CGOBegin(cgo, GL_TRIANGLES);
          continue;  // for-loop, no add_to_cgo
      }

      inside_begin_triangles = false;
    } else if (op == CGO_END) {
      inside_begin_triangles = false;
    }

    if (!inside_begin_triangles) {
      cgo->add_to_cgo(op, pc);
      continue;
    }

    // handle operations inside BEGIN/END TRIANGLE(S|_STRIP|_FAN)
    switch (op) {
      case CGO_VERTEX:
        copy3(reinterpret_cast<const float*>(pc), vertices[current_i]);
        copy3(current_color, colors[current_i]);
        alphas[current_i] = current_alpha;

        ++vertex_count;

        switch (mode) {
          case GL_TRIANGLE_STRIP:
            current_i = vertex_count % 3;
            emit = (vertex_count > 2);
            break;
          case GL_TRIANGLE_FAN:
            current_i = ((vertex_count + 1) % 2) + 1;
            emit = (vertex_count > 2);
            break;
          default:
            current_i = vertex_count % 3;
            emit = (current_i == 0);
        }

        if (emit) {
          auto* indices = flip ? indices_flipped : indices_regular;

          if (mode != GL_TRIANGLES) {
            flip = !flip;
          }

          CalculateTriangleNormal(vertices[0], vertices[indices[1]],
                                  vertices[indices[2]], current_normal);

          CGONormalv(cgo, current_normal);

          for (int j = 0; j < 3; ++j) {
            int k = indices[j];
            if (has_color) CGOColorv(cgo, colors[k]);
            if (has_alpha) CGOAlpha(cgo, alphas[k]);
            CGOVertexv(cgo, vertices[k]);
          }
        }

        break;
      case CGO_COLOR:
        copy3(reinterpret_cast<const float*>(pc), current_color);
        has_color = true;
        break;
      case CGO_ALPHA:
        current_alpha = *reinterpret_cast<const float*>(pc);
        has_alpha = true;
        break;
      case CGO_NORMAL:
        // discard, we will generate new normals
        break;
      default:
        PRINTFB(G, FB_CGO, FB_Warnings)
        " CGO-Warning: CGOGenerateNormalsForTriangles: unhandled op=0x%02x "
        "inside BEGIN/END\n",
            op ENDFB(G);
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
  return (cgo);
}

CGO* CGOTurnLightingOnLinesOff(const CGO* I, bool use_shader) {
  bool cur_mode_is_lines = false;
  auto cgo = CGONewSized(I->G, I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_DRAW_ARRAYS: {
        const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        float* vals;
        int nvals = sp->narrays * sp->nverts;
        switch (sp->mode) {
          case GL_LINES:
          case GL_LINE_STRIP:
            CGODisable(cgo, CGO_GL_LIGHTING);
            cur_mode_is_lines = true;
        }
        vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
        memcpy(vals, sp->floatdata, nvals);
        if (cur_mode_is_lines) {
          CGOEnable(cgo, CGO_GL_LIGHTING);
          cur_mode_is_lines = false;
        }
      } break;
      case CGO_DRAW_BUFFERS_INDEXED: {
        const cgo::draw::buffers_indexed* sp =
            reinterpret_cast<decltype(sp)>(pc);
        int mode = sp->mode, mode_is_lines = 0;
        switch (mode) {
          case GL_LINES:
          case GL_LINE_STRIP:
            mode_is_lines = true;
        }
        if (mode_is_lines) {
          CGODisable(cgo, CGO_GL_LIGHTING);
        }
        cgo->copy_op_from<cgo::draw::buffers_indexed>(pc);
        if (mode_is_lines) {
          CGOEnable(cgo, CGO_GL_LIGHTING);
        }
      } break;
      case CGO_DRAW_BUFFERS_NOT_INDEXED: {
        const cgo::draw::buffers_not_indexed* sp =
            reinterpret_cast<decltype(sp)>(pc);
        int mode = sp->mode, mode_is_lines = 0;
        switch (mode) {
          case GL_LINES:
          case GL_LINE_STRIP:
            mode_is_lines = true;
        }
        if (mode_is_lines) {
          CGODisable(cgo, CGO_GL_LIGHTING);
        }
        cgo->copy_op_from<cgo::draw::buffers_not_indexed>(pc);
        if (mode_is_lines) {
          CGOEnable(cgo, CGO_GL_LIGHTING);
        }
      } break;
      case CGO_END: {
        CGOEnd(cgo);
        if (cur_mode_is_lines) {
          CGOEnable(cgo, CGO_GL_LIGHTING);
          cur_mode_is_lines = 0;
        }
      } break;
      case CGO_BEGIN: {
        int mode = CGO_get_int(pc);
        switch (mode) {
          case GL_LINES:
          case GL_LINE_STRIP:
            CGODisable(cgo, CGO_GL_LIGHTING);
            cur_mode_is_lines = true;
            break;
          default:
            if (!use_shader) {  // no shaders, not lines, turn lighting on
              CGOEnable(cgo, CGO_GL_LIGHTING);
            }
        }
        CGOBegin(cgo, mode);
      } break;
      default:
        cgo->add_to_cgo(op, pc);
    }
  }
  cgo->use_shader = use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = CGOPackedNormals(cgo->G);
  }
  return (cgo);
}

bool CGOHasAnyTriangleVerticesWithoutNormals(const CGO* I,
                                             bool checkTriangles) {
  bool inside = false;
  bool hasNormal = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
      case CGO_BEGIN:
        switch (CGO_get_int(pc)) {
          case GL_TRIANGLE_FAN:
          case GL_TRIANGLES:
          case GL_TRIANGLE_STRIP:
            if (checkTriangles) inside = 1;
            break;
          case GL_LINE_STRIP:
          case GL_LINES:
            if (!checkTriangles) inside = 1;
            break;
        }
        break;
      case CGO_END:
        inside = 0;
        break;
      case CGO_NORMAL:
        hasNormal = 1;
        break;
      case CGO_VERTEX:
        if (inside && !hasNormal) return 1;
        break;
      case CGO_DRAW_ARRAYS: {
        const auto sp = it.cast<cgo::draw::arrays>();
        switch (sp->mode) {
          case GL_TRIANGLE_FAN:
          case GL_TRIANGLES:
          case GL_TRIANGLE_STRIP:
            if (checkTriangles) {
              if (!(sp->arraybits & CGO_NORMAL_ARRAY)) {
                return 1;
              }
            }
            break;
          case GL_LINE_STRIP:
          case GL_LINES:
            if (!checkTriangles) {
              if (!(sp->arraybits & CGO_NORMAL_ARRAY)) {
                return 1;
              }
            }
            break;
        }
      } break;
    }
  }
  return 0;
}

void CGO::add_to_cgo(int op, const float* pc) {
  switch (op) {
    case CGO_STOP:
      // only append to buffer, don't increment size
      CGOStop(this);
      break;
    case CGO_DRAW_ARRAYS:
      copy_op_from<cgo::draw::arrays>(pc);
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
      copy_op_from<cgo::draw::buffers_indexed>(pc);
      break;
    case CGO_DRAW_TEXTURES:
      copy_op_from<cgo::draw::textures>(pc);
      break;
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      copy_op_from<cgo::draw::screen_textures>(pc);
      break;
    case CGO_DRAW_LABELS:
      copy_op_from<cgo::draw::labels>(pc);
      break;
    case CGO_DRAW_CONNECTORS:
      copy_op_from<cgo::draw::connectors>(pc);
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      copy_op_from<cgo::draw::buffers_not_indexed>(pc);
      break;
    case CGO_DRAW_SPHERE_BUFFERS:
      copy_op_from<cgo::draw::sphere_buffers>(pc);
      break;
    case CGO_DRAW_CYLINDER_BUFFERS:
      copy_op_from<cgo::draw::cylinder_buffers>(pc);
      break;
    case CGO_DRAW_CUSTOM:
      copy_op_from<cgo::draw::custom>(pc);
      break;
    default:
      int sz = CGO_sz[op];
      std::copy_n(pc - 1, sz + 1, add_to_buffer(sz + 1));
  };
}

cgo::draw::shadercylinder2ndcolor::shadercylinder2ndcolor(
    CGO* I, const float* _origin, const float* _axis, const float _tube_size,
    int _cap, const float* _color2, Pickable* pickcolor2, const float _alpha)
    : tube_size(_tube_size), alpha(_alpha) {
  copy3f(_origin, origin);
  copy3f(_axis, axis);
  cap = _cap;
  copy3f(_color2, color2);
  if (pickcolor2) {
    I->current_pick_color_index = pick_color_index = pickcolor2->index;
    I->current_pick_color_bond = pick_color_bond = pickcolor2->bond;
  } else {
    pick_color_index = I->current_pick_color_index;
    pick_color_bond = I->current_pick_color_bond;
  }
};

/**
 * CGO iterator increment
 */
CGO::const_iterator& CGO::const_iterator::operator++() {
  const unsigned op = op_code();

  // Corrupted OP codes should never make it into a CGO, so we don't want to
  // handle this gracefully. Only import of CGOs (from PSEs, from Python, etc.)
  // should handle invalid codes gracefully.
  assert(op < CGO_sz_size());

  m_pc += CGO_sz[op] + 1;
  return *this;
}

void CGORender(CGO* I, const float* color, CSetting* set1, CSetting* set2,
               RenderInfo* info, Rep* rep) {
  CGORenderGL(I, color, set1, set2, info, rep);
}

void CGORenderAlpha(CGO* I, RenderInfo* info, bool calcDepth) {
  CGORenderGLAlpha(I, info, calcDepth);
}
