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

// Rendering CGOs with the ray tracer.

#include "pymol/graphics/cgo.h"
#include "pymol/graphics/cgo_renderer.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/scene/scene.h"
#include "pymol/graphics/cgo_internal.h"

static int CGORenderQuadricRay(CRay* ray, float* v, float r, float* q) {
  float r_el, n0[3], n1[3], n2[3];
  int ok = true;
  if (CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ok &= ray->ellipsoid3fv(v, r_el, n0, n1, n2);
  return ok;
}

/* ======== Raytrace Renderer ======== */

int CGORenderRay(CGO* I, CRay* ray, RenderInfo* info, const float* color,
                 ObjectGadgetRamp* ramp, CSetting* set1, CSetting* set2) {
  int vc = 0;
  float linewidth = 1.0F;
  float widthscale = 0.15F;
  float lineradius, dotradius, dotwidth;
  float white[] = {1.0, 1.0, 1.0};
  float zee[] = {0.0, 0.0, 1.0};
  int ok = true;
  const float *n0 = nullptr, *n1 = nullptr, *n2 = nullptr, *v0 = nullptr,
              *v1 = nullptr, *v2 = nullptr, *c0 = nullptr, *c1 = nullptr,
              *c2 = nullptr;
  float rampc0[3], rampc1[3], rampc2[3];
  int mode = -1;
  /* workaround; multi-state ray-trace bug */
  if (!I) {
    assert("TODO investigate" && false);
    return 0; /* not sure if it should return 0 or 1, 0 - fails, but is it a
                 memory issue? might not be since the arg is nullptr */
  }

  I->G->CGORenderer->alpha =
      1.0F - SettingGet_f(I->G, set1, set2, cSetting_cgo_transparency);

  widthscale = SettingGet_f(I->G, set1, set2, cSetting_cgo_ray_width_scale);

  /*  printf("debug %8.9f\n",SceneGetScreenVertexScale(I->G,zee)); */
  linewidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_width);
  if (linewidth < 0.0F) linewidth = 1.0F;
  lineradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_radius);
  dotwidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_width);
  dotradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_radius);
  if (lineradius < 0.0F) lineradius = linewidth * ray->PixelRadius / 2.0F;
  if (dotradius < 0.0F) dotradius = dotwidth * ray->PixelRadius / 2.0F;
  if (widthscale < 0.0F) widthscale = ray->PixelRadius / 2.0F;
  if (color)
    c0 = color;
  else
    c0 = white;
  ray->transparentf(1.0F - I->G->CGORenderer->alpha);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
      case CGO_BEGIN:
        mode = CGO_get_int(pc);
        vc = 0;
        n0 = zee;
        break;
      case CGO_END:
        switch (mode) {
          case GL_LINE_LOOP:
            if (vc > 1) ok &= ray->sausage3fv(v0, v2, lineradius, c0, c2);
            break;
        }
        mode = -1;
        break;
      case CGO_WIDTHSCALE:
        widthscale = *pc;
        lineradius = widthscale * linewidth;
        dotradius = widthscale * dotwidth;
        break;
      case CGO_DOTWIDTH:
        dotwidth = *pc;
        dotradius = widthscale * dotwidth;
        break;
      case CGO_LINEWIDTH:
        linewidth = *pc;
        lineradius = widthscale * linewidth;
        break;
      case CGO_NORMAL:
        n0 = pc;
        break;
      case CGO_SPECIAL_WITH_ARG: {
        float argval = *(pc + 1);
        switch (*(int*)pc) {
          case LINEWIDTH_FOR_LINES:
            linewidth = argval;
            lineradius = widthscale * linewidth;
        }
      } break;
      case CGO_SPECIAL: {
        switch ((*(int*)pc)) {
          case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON: {
            float radius =
                SettingGet_f(I->G, set1, set2, cSetting_ribbon_radius);
            if (radius == 0.0F) {
              float ribbon_width =
                  SettingGet_f(I->G, set1, set2, cSetting_ribbon_width);
              float line_width = SceneGetDynamicLineWidth(info, ribbon_width);
              SceneGetDynamicLineWidth(info, line_width);
              radius = ray->PixelRadius * line_width / 2.0F;
            }
            lineradius = radius;
          } break;
          case LINEWIDTH_FOR_LINES: {
            float radius = SettingGet_f(I->G, set1, set2, cSetting_line_radius);
            if (radius <= 0.0F) {
              float line_width =
                  SettingGet_f(I->G, set1, set2, cSetting_line_width);
              line_width = SceneGetDynamicLineWidth(info, line_width);
              radius = ray->PixelRadius * line_width / 2.0F;
            }
            lineradius = radius;
          } break;
          case CYLINDER_WIDTH_FOR_NONBONDED:
          case LINEWIDTH_WITH_SCALE: {
            float line_width =
                SettingGet_f(I->G, set1, set2, cSetting_line_width);
            line_width = SceneGetDynamicLineWidth(info, line_width);
            lineradius = widthscale * line_width / 2.f;
          } break;
        }
      } break;
      case CGO_COLOR:
        c0 = pc;
        ray->color3fv(c0);
        break;
      case CGO_ALPHA:
        I->G->CGORenderer->alpha = *pc;
        ray->transparentf(1.0F - *pc);
        break;
      case CGO_LINE: {
        auto line = reinterpret_cast<cgo::draw::line*>(pc);
        ok &= ray->sausage3fv(line->vertex1, line->vertex2, lineradius, c0, c0);
      } break;
      case CGO_SPLITLINE: {
        auto splitline = reinterpret_cast<cgo::draw::splitline*>(pc);
        float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
                          CONVERT_COLOR_VALUE(splitline->color2[1]),
                          CONVERT_COLOR_VALUE(splitline->color2[2])};
        if (splitline->flags & cgo::draw::splitline::interpolation) {
          ok &= ray->sausage3fv(splitline->vertex1, splitline->vertex2,
                                lineradius, c0, color2);
        } else {
          float mid[3];
          add3f(splitline->vertex1, splitline->vertex2, mid);
          mult3f(mid, .5f, mid);
          ok &= ray->customCylinder3fv(splitline->vertex1, mid, lineradius, c0,
                                       c0, cCylCap::Round, cCylCap::None);
          ok &= ray->customCylinder3fv(mid, splitline->vertex2, lineradius,
                                       color2, color2, cCylCap::None,
                                       cCylCap::Round);
        }
      } break;
      case CGO_VERTEX_CROSS: {
        float pt1[3], pt2[3];
        float nonbonded_size =
            SettingGet_f(I->G, set1, set2, cSetting_nonbonded_size);
        copy3f(pc, pt1);
        copy3f(pc, pt2);
        pt1[0] -= nonbonded_size;
        pt2[0] += nonbonded_size;
        ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);

        copy3f(pc, pt1);
        copy3f(pc, pt2);
        pt1[1] -= nonbonded_size;
        pt2[1] += nonbonded_size;
        ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);

        copy3f(pc, pt1);
        copy3f(pc, pt2);
        pt1[2] -= nonbonded_size;
        pt2[2] += nonbonded_size;
        ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);
      } break;
      case CGO_VERTEX_BEGIN_LINE_STRIP:
      case CGO_VERTEX:
        v0 = pc;

        if (ramp) {
          if (!ObjectGadgetRampInterVertex(ramp, v0, rampc0, -1)) {
            copy3f(white, rampc0);
          }
          c0 = rampc0;
        }
        switch (mode) {
          case GL_POINTS:
            ok &= ray->sphere3fv(v0, dotradius);
            break;
          case GL_LINES:
            if (vc & 0x1) ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
            v1 = v0;
            if (!ramp) {
              c1 = c0;
            }
            break;
          case GL_LINE_STRIP:
            if (vc) {
              ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
            }
            v1 = v0;
            if (!ramp) {
              c1 = c0;
            }
            break;
          case GL_LINE_LOOP:
            if (vc)
              ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
            else {
              v2 = v0;
              c2 = c0;
            }
            v1 = v0;
            if (!ramp) c1 = c0;
            break;
          case GL_TRIANGLES:
            if (((vc + 1) % 3) == 0)
              ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
            v2 = v1;
            n2 = n1;
            v1 = v0;
            n1 = n0;
            if (!ramp) {
              c2 = c1;
              c1 = c0;
            }
            break;
          case GL_TRIANGLE_STRIP:
            if (vc > 1)
              ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
            v2 = v1;
            n2 = n1;
            v1 = v0;
            n1 = n0;
            if (!ramp) {
              c2 = c1;
              c1 = c0;
            }
            break;
          case GL_TRIANGLE_FAN:
            if (vc > 1)
              ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
            else if (!vc) {
              n2 = n0;
              v2 = v0;
              if (!ramp) c2 = c0;
            }
            v1 = v0;
            n1 = n0;
            if (!ramp) c1 = c0;
            break;
        }
        if (ramp) {
          switch (mode) {
            case GL_TRIANGLES:
            case GL_TRIANGLE_STRIP:
            case GL_TRIANGLE_FAN:
              copy3f(rampc1, rampc2);
              c2 = rampc2;
            case GL_LINES:
            case GL_LINE_STRIP:
            case GL_LINE_LOOP:
              copy3f(rampc0, rampc1);
              c1 = rampc1;
              break;
          }
        }
        vc++;
        break;
      case CGO_SPHERE:
        ray->color3fv(c0);
        ok &= ray->sphere3fv(pc, *(pc + 3));
        break;
      case CGO_ELLIPSOID:
        ray->color3fv(c0);
        ok &= ray->ellipsoid3fv(pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
        break;
      case CGO_QUADRIC:
        ray->color3fv(c0);
        ok &= CGORenderQuadricRay(ray, pc, *(pc + 3), pc + 4);
        break;
      case CGO_CONE:
        ok &= ray->cone3fv(pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
                           static_cast<cCylCap>(int(pc[14])),
                           static_cast<cCylCap>(int(pc[15])));
        break;
      case CGO_CUSTOM_CYLINDER: {
        auto cyl = reinterpret_cast<cgo::draw::custom_cylinder*>(pc);
        ok &= ray->customCylinder3fv(*cyl);
      } break;
      case CGO_CUSTOM_CYLINDER_ALPHA: {
        auto cyl = reinterpret_cast<cgo::draw::custom_cylinder_alpha*>(pc);
        ok &= ray->customCylinderAlpha3fv(*cyl);
      } break;
      case CGO_SHADER_CYLINDER: {
        float p2[3];
        int cap = CGO_get_int(pc + 7);
        const cCylCap cap1 = cap1_from_cyl_shader_bits(cap);
        const cCylCap cap2 = cap2_from_cyl_shader_bits(cap);
        add3f(pc, pc + 3, p2);
        ok &= ray->customCylinder3fv(pc, p2, *(pc + 6), ray->CurColor,
                                     ray->CurColor, cap1, cap2);
      } break;
      case CGO_SHADER_CYLINDER_WITH_2ND_COLOR: {
        auto cyl = reinterpret_cast<cgo::draw::shadercylinder2ndcolor*>(pc);
        float v1[3];
        int cap = cyl->cap;
        const cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        const cCylCap bcap = cap2_from_cyl_shader_bits(cap);
        int colorinterp = cap & cCylShaderInterpColor;
        const float* color1 = c0;
        const float* color2 = cyl->color2;
        add3f(cyl->origin, cyl->axis, v1);
        float alpha1 = I->G->CGORenderer->alpha;
        float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : alpha1;
        if (colorinterp || equal3f(color1, color2)) {
          ok &= ray->customCylinder3fv(pc, v1, cyl->tube_size, color1, color2,
                                       fcap, bcap, alpha1, alpha2);
        } else {
          float mid[3];
          mult3f(cyl->axis, .5f, mid);
          add3f(cyl->origin, mid, mid);

          ray->color3fv(c0);
          ok &= ray->customCylinder3fv(cyl->origin, mid, cyl->tube_size, color1,
                                       color1, fcap, cCylCap::None, alpha1,
                                       alpha2);
          ray->color3fv(cyl->color2);
          ok &= ray->customCylinder3fv(mid, v1, cyl->tube_size, color2, color2,
                                       cCylCap::None, bcap, alpha1, alpha2);
        }
      } break;
      case CGO_CYLINDER: {
        auto* cyl = reinterpret_cast<cgo::draw::cylinder*>(pc);
        ok &= ray->cylinder3fv(*cyl);
      } break;
      case CGO_SAUSAGE:
        ok &= ray->sausage3fv(pc, pc + 3, *(pc + 6), pc + 7, pc + 10);
        break;
      case CGO_TRIANGLE:
        ok &= ray->triangle3fv(pc, pc + 3, pc + 6, pc + 9, pc + 12, pc + 15,
                               pc + 18, pc + 21, pc + 24);
        break;
      case CGO_DRAW_ARRAYS: {
        cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
        int const mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
        float* vertexVals = sp->floatdata;
        float *normalVals = 0, *colorVals = 0;
        int offset = 0;
        if (arrays & CGO_VERTEX_ARRAY) {
          vertexVals = sp->floatdata;
          offset += nverts * VERTEX_POS_SIZE;
        }
        if (arrays & CGO_NORMAL_ARRAY) {
          normalVals = sp->floatdata + offset;
          offset += nverts * VERTEX_NORMAL_SIZE;
        }
        if (arrays & CGO_COLOR_ARRAY) {
          colorVals = sp->floatdata + offset;
          offset += nverts * VERTEX_COLOR_SIZE;
        }
        vc = 0;
        for (int v = 0, pl = 0, plc = 0; ok && v < nverts;
             v++, pl += 3, plc += 4) {
          if (normalVals) {
            n0 = &normalVals[pl];
          }
          if (colorVals) {
            c0 = &colorVals[plc];
            ray->color3fv(c0);
            ray->transparentf(1.0f - c0[3]);
          }
          if (vertexVals) {
            v0 = &vertexVals[pl];
          }
          switch (mode) {
            case GL_POINTS:
              ok &= ray->sphere3fv(v0, dotradius);
              break;
            case GL_LINES:
              if (vc & 0x1) ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
              v1 = v0;
              c1 = c0;
              break;
            case GL_LINE_STRIP:
              if (vc) ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
              v1 = v0;
              c1 = c0;
              break;
            case GL_LINE_LOOP:
              if (vc)
                ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
              else {
                v2 = v0;
                c2 = c0;
              }
              v1 = v0;
              c1 = c0;
              break;
            case GL_TRIANGLES:
              if (((vc + 1) % 3) == 0)
                ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
              v2 = v1;
              c2 = c1;
              n2 = n1;
              v1 = v0;
              c1 = c0;
              n1 = n0;
              break;
            case GL_TRIANGLE_STRIP:
              if (vc > 1)
                ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
              v2 = v1;
              c2 = c1;
              n2 = n1;
              v1 = v0;
              c1 = c0;
              n1 = n0;
              break;
            case GL_TRIANGLE_FAN:
              if (vc > 1)
                ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
              else if (!vc) {
                n2 = n0;
                v2 = v0;
                c2 = c0;
              }
              v1 = v0;
              c1 = c0;
              n1 = n0;
              break;
          }
          vc++;
        }
      } break;
      default:
        break;
    }
  }

  if (ok) ray->transparentf(0.0F);
  return ok;
}
