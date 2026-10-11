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

// Exporting the ray scene as G3d, VRML, IDTF, OBJ/MTL and POV-Ray.

#include "pymol/ray/ray.h"
#include "pymol/geometry/triangle.h"
#include "pymol/color/color.h"
#include "pymol/geometry/matrix.h"
#include "pymol/scene/scene.h"
#include "pymol/base/feedback.h"
#include "pymol/ray/ray_internal.h"
G3dPrimitive* RayRenderG3d(CRay* I, int width, int height, float front,
                           float back, float fov, int quiet) {
  /* generate a rendering stream for Miguel's G3d java rendering engine */

  float scale_x, scale_y;
  int shift_x, shift_y;
  float* d;
  CBasis* base;
  CPrimitive* prim;
  float* vert;
  float vert2[3];
  int a;
  G3dPrimitive *jprim = VLAlloc(G3dPrimitive, 10000), *jp;
  int n_jp = 0;

#define convert_r(r) 2 * (int)(r * scale_x);
#define convert_x(x) shift_x + (int)(x * scale_x);
#define convert_y(y) height - (shift_y + (int)(y * scale_y));
#define convert_z(z) -(int)((z + front) * scale_x);
#define convert_col(c)                                                         \
  (0xFF000000 | (((int)(c[0] * 255.0)) << 16) | (((int)(c[1] * 255.0)) << 8) | \
   (((int)(c[2] * 255.0))))

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, false);

  if (!quiet) {
    PRINTFB(I->G, FB_Ray, FB_Details)
    " RayRenderG3d: processed %i graphics primitives.\n",
        I->NPrimitive ENDFB(I->G);
  }
  base = I->Basis + 1;

  /* always orthoscopic */

  /* front should give a zero Z,
     -I->Range[0] should be off the right hand size
     I->Range[1] should be off the top */
  scale_x = width / I->Range[0];
  scale_y = height / I->Range[1];
  shift_x = width / 2;
  shift_y = height / 2;

  for (a = 0; a < I->NPrimitive; a++) {
    prim = I->Primitive + a;
    vert = base->Vertex + 3 * (prim->vert);
    switch (prim->type) {
      case cPrimSphere:
        VLACheck(jprim, G3dPrimitive, n_jp);
        jp = jprim + n_jp;
        jp->op = 1;
        jp->r = convert_r(prim->r1);
        jp->x1 = convert_x(vert[0]);
        jp->y1 = convert_y(vert[1]);
        jp->z1 = convert_z(vert[2]);
        jp->c = convert_col(prim->c1);
        n_jp++;
        break;
      case cPrimSausage:
        VLACheck(jprim, G3dPrimitive, n_jp);
        d = base->Normal + 3 * base->Vert2Normal[prim->vert];
        scale3f(d, prim->l1, vert2);
        add3f(vert, vert2, vert2);

        jp = jprim + n_jp;
        jp->op = 3;
        jp->r = convert_r(prim->r1);
        jp->x1 = convert_x(vert[0]);
        jp->y1 = convert_y(vert[1]);
        jp->z1 = convert_z(vert[2]);
        jp->x2 = convert_x(vert2[0]);
        jp->y2 = convert_y(vert2[1]);
        jp->z2 = convert_z(vert2[2]);
        jp->c = convert_col(prim->c1);
        n_jp++;
        break;
      case cPrimTriangle:
        VLACheck(jprim, G3dPrimitive, n_jp);
        jp = jprim + n_jp;
        jp->op = 2;
        jp->x1 = convert_x(vert[0]);
        jp->y1 = convert_y(vert[1]);
        jp->z1 = convert_z(vert[2]);
        jp->x2 = convert_x(vert[3]);
        jp->y2 = convert_y(vert[4]);
        jp->z2 = convert_z(vert[5]);
        jp->x3 = convert_x(vert[6]);
        jp->y3 = convert_y(vert[7]);
        jp->z3 = convert_z(vert[8]);
        jp->c = convert_col(prim->c1);
        n_jp++;
        break;
    }
  }
  VLASize(jprim, G3dPrimitive, n_jp);
  return jprim;
}

void RayRenderVRML1(CRay* I, int width, int height, char** vla_ptr, float front,
                    float back, float fov, float angle, float z_corr) {
  char* vla = *vla_ptr;
  ov_size cc = 0; /* character count */
  OrthoLineType buffer;

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, false);

  strcpy(buffer, "#VRML V1.0 ascii\n\n");
  UtilConcatVLA(&vla, &cc, buffer);

  UtilConcatVLA(&vla, &cc, "MaterialBinding { value OVERALL }\n");

  sprintf(buffer,
          "Material {\n ambientColor 0 0 0\n diffuseColor 1 1 1\n "
          "specularColor 1 1 1\nshininess 0.2\n}\n");
  UtilConcatVLA(&vla, &cc, buffer);

  {
    int a;
    CPrimitive* prim;
    float* vert;
    CBasis* base = I->Basis + 1;

    UtilConcatVLA(&vla, &cc, "Separator {\n");

    UtilConcatVLA(&vla, &cc, "MatrixTransform {\n");
    UtilConcatVLA(&vla, &cc, "matrix 1.0 0.0 0.0 0.0\n");
    UtilConcatVLA(&vla, &cc, "       0.0 1.0 0.0 0.0\n");
    UtilConcatVLA(&vla, &cc, "       0.0 0.0 1.0 0.0\n");
    sprintf(buffer, "    %8.6f %8.6f %8.6f 1.0\n",
            (I->Volume[0] + I->Volume[1]) / 2,
            (I->Volume[2] + I->Volume[3]) / 2, 0.0F);
    UtilConcatVLA(&vla, &cc, buffer);
    UtilConcatVLA(&vla, &cc, "}\n");

    for (a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);
      switch (prim->type) {
        case cPrimSphere:
          sprintf(buffer, "Material {\ndiffuseColor %6.4f %6.4f %6.4f\n}\n\n",
                  prim->c1[0], prim->c1[1], prim->c1[2]);
          UtilConcatVLA(&vla, &cc, buffer);
          UtilConcatVLA(&vla, &cc, "Separator {\n");
          sprintf(buffer,
                  "Transform {\ntranslation %8.6f %8.6f %8.6f\nscaleFactor "
                  "%8.6f %8.6f %8.6f\n}\n",
                  vert[0], vert[1], vert[2] - z_corr, prim->r1, prim->r1,
                  prim->r1);
          UtilConcatVLA(&vla, &cc, buffer);
          sprintf(buffer, "Sphere {}\n");
          UtilConcatVLA(&vla, &cc, buffer);
          UtilConcatVLA(&vla, &cc, "}\n\n");
          break;
        case cPrimCylinder:
        case cPrimSausage:
        case cPrimCone:
          break;
        case cPrimTriangle:
          break;
      }
    }

    UtilConcatVLA(&vla, &cc, "}\n");
  }

  *vla_ptr = vla;
}

int TriangleReverse(CPrimitive* p) {
  float s1[3], s2[3], n0[3];

  subtract3f(p->v1, p->v2, s1);
  subtract3f(p->v3, p->v2, s2);
  cross_product3f(s1, s2, n0);

  if (dot_product3f(p->n0, n0) < 0.0F)
    return 0;
  else
    return 1;
}

void RayRenderVRML2(CRay* I, int width, int height, char** vla_ptr, float front,
                    float back, float fov, float angle, float z_corr) {
  /*

     From: pymol-users-admin@lists.sourceforge.net on behalf of Chris Want
     Sent: Tuesday, February 07, 2006 1:47 PM
     To: pymol-users@lists.sourceforge.net
     Subject: [PyMOL] VRML patch

     Hi Warren,

     I took your advice and modified the RayRenderVRML2() function to
     support triangles. I also threw out the sphere code that was already
     there and rewrote it (the code there was for VRML1, not
     VRML2). While I was at it, I also implemented export for cylinders
     and sausages.

     The code in the attached patch (diff-ed against cvs, and tested with
     two VRML2 readers) can be regarded as being in the public domain.

     Regards,
     Chris

     cwant_at_ualberta.ca

   */

  char* vla = *vla_ptr;
  ov_size cc = 0; /* character count */
  OrthoLineType buffer;
  float mid[3]; /*, wid[3]; */
  float h_fov = cPI * (fov * width) / (180 * height);
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);
  RayComputeBox(I);
  /*
     mid[0] = (I->max_box[0] + I->min_box[0]) / 2.0;
     mid[1] = (I->max_box[1] + I->min_box[1]) / 2.0;
     mid[2] = (I->max_box[2] + I->min_box[2]) / 2.0;
     wid[0]  = (I->max_box[0] - I->min_box[0]);
     wid[1]  = (I->max_box[1] - I->min_box[1]);
     wid[2]  = (I->max_box[2] - I->min_box[2]);
   */
  auto Pos = glm::value_ptr(I->Pos);
  copy3f(Pos, mid);
  UtilConcatVLA(&vla, &cc,
                "#VRML V2.0 utf8\n" /* WLD: most VRML2 readers req. utf8 */
                "\n");
  if (!identity) {
    sprintf(buffer,
            "Viewpoint {\n"
            " position 0 0 %6.8f\n"
            " orientation 1 0 0 0\n"
            " description \"Z view\"\n"
            " fieldOfView %8.6f\n" /* WLD: use correct FOV */
            "}\n"
            /* WLD: only write the viewpoint which matches PyMOL
               "Viewpoint {\n"
               " position %6.8f 0 0\n"
               " orientation 0 1 0 1.570796\n"
               " description \"X view\"\n"
               "}\n"
               "Viewpoint {\n"
               " position 0 %6.8f 0\n"
               " orientation 0 -0.707106 -0.7071061 3.141592\n"
               " description \"Y view\"\n"
               "}\n" */
            ,
            -z_corr, /* *0.96646  for some reason, PyMOL and C4D cameras differ
                        by about 3.5% ... */
            h_fov
            /*(wid[2] + wid[1]),
               (wid[0] + wid[1]),
               (wid[1] + wid[2]) */
    );
    UtilConcatVLA(&vla, &cc, buffer);
  }
  if (!identity) {
    float light[3];
    auto lightv = SettingGet<const float*>(I->G, cSetting_light);
    copy3f(lightv, light);
    normalize3f(light);
    sprintf(buffer,
            "DirectionalLight {\n"
            " direction %8.6f %8.6f %8.3f\n"
            "}\n",
            light[0], light[1], light[2]);
    UtilConcatVLA(&vla, &cc, buffer);
  }
  UtilConcatVLA(&vla, &cc,
                "NavigationInfo {\n"
                " headlight TRUE\n"
                " type \"EXAMINE\"\n"
                "}\n");
  {
    int a, b;
    CPrimitive* prim;
    float* vert;
    int mesh_obj = false, mesh_start = 0;

    CBasis* base = I->Basis + 1;

    for (a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);

      if (prim->type == cPrimTriangle) {
        if (!mesh_obj) {
          /* start mesh */
          mesh_start = a;
          UtilConcatVLA(&vla, &cc,
                        "Shape {\n"
                        " appearance Appearance {\n"
                        "  material Material { diffuseColor 1.0 1.0 1.0 }\n"
                        " }\n"
                        " geometry IndexedFaceSet {\n"
                        "  coord Coordinate {\n"
                        "   point [\n");
          mesh_obj = true;
        }
      } else if (mesh_obj) {
        CPrimitive* cprim;
        int tri = 0;
        /* output connectivity */
        UtilConcatVLA(&vla, &cc,
                      "   ]\n"
                      "  }\n"
                      "  coordIndex [\n");
        for (b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          if (TriangleReverse(cprim))
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
          else
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
          UtilConcatVLA(&vla, &cc, buffer);
          tri += 3;
        }

        /* output vertex colors */
        UtilConcatVLA(&vla, &cc,
                      "  ]\n"
                      "  colorPerVertex TRUE\n"
                      "  color Color {\n"
                      "   color [\n");
        for (b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          sprintf(buffer,
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n",
                  cprim->c1[0], cprim->c1[1], cprim->c1[2], cprim->c2[0],
                  cprim->c2[1], cprim->c2[2], cprim->c3[0], cprim->c3[1],
                  cprim->c3[2]);
          UtilConcatVLA(&vla, &cc, buffer);
        }

        /* output vertex normals */
        UtilConcatVLA(&vla, &cc,
                      "  ] } \n"
                      "  normalPerVertex TRUE\n"
                      "  normal Normal {\n"
                      "   vector [\n");
        for (b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          {
            float* norm = base->Normal + 3 * base->Vert2Normal[cprim->vert];
            sprintf(buffer,
                    "%6.4f %6.4f %6.4f,\n"
                    "%6.4f %6.4f %6.4f,\n"
                    "%6.4f %6.4f %6.4f,\n",
                    norm[3], norm[4], norm[5],    /* transformed cprim->n1 */
                    norm[6], norm[7], norm[8],    /* transformed cprim->n2 */
                    norm[9], norm[10], norm[11]); /* transformed cprim->n3 */
            UtilConcatVLA(&vla, &cc, buffer);
          }
        }
        UtilConcatVLA(&vla, &cc,
                      "  ] }\n"
                      "  normalIndex [ \n");
        tri = 0;
        for (b = mesh_start; b < a; b++) {
          cprim = I->Primitive + b;
          if (TriangleReverse(cprim))
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
          else
            sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
          UtilConcatVLA(&vla, &cc, buffer);
          tri += 3;
        }

        /* close mesh */
        UtilConcatVLA(&vla, &cc,
                      " ] \n"
                      " }\n"
                      "}\n");
        mesh_obj = false;
      }

      switch (prim->type) {
        case cPrimSphere:
          sprintf(buffer,
                  "Transform {\n"
                  " translation %8.6f %8.6f %8.6f\n"
                  " children Shape {\n"
                  "  geometry Sphere { radius %8.6f }\n"
                  "  appearance Appearance {\n"
                  "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                  "                       specularColor 0.8 0.8 0.8 \n"
                  "                       shininess 0.8 }\n"
                  "  }\n"
                  " }\n"
                  "}\n",
                  vert[0] - mid[0], vert[1] - mid[1], vert[2] - mid[2],
                  prim->r1, prim->c1[0], prim->c1[1], prim->c1[2]);
          UtilConcatVLA(&vla, &cc, buffer);
          break;
        case cPrimCone:
          /* TO DO */
          break;
        case cPrimCylinder:
        case cPrimSausage: {
          float *d, vert2[3], axis[3], angle;
          OrthoLineType geometry;
          /* find the axis and angle that will rotate the y axis onto
           * the direction of the length of the cylinder
           */
          d = base->Normal + 3 * base->Vert2Normal[prim->vert];
          if ((d[0] * d[0] + d[2] * d[2]) < 0.000001) {
            /* parallel with y */
            axis[0] = 1.0;
            axis[1] = 0.0;
            axis[2] = 0.0;
            if (d[1] > 0) {
              angle = 0.0;
            } else {
              angle = cPI;
            }
          } else {
            axis[0] = d[2];
            axis[1] = 0.0;
            axis[2] = -d[0];
            normalize3f(axis);
            angle = d[1];
            if (angle > 1.0)
              angle = 1.0;
            else if (angle < -1.0)
              angle = -1.0;
            angle = acos(angle);
          }
          /* vrml cylinders have origin in middle, not tip, that is why we
           * use prim->l1/2
           */
          scale3f(d, prim->l1 / 2, vert2);
          add3f(vert, vert2, vert2);
          if (prim->type == cPrimSausage) {
            OrthoLineType geom_add;
            sprintf(geometry,
                    "  Shape {\n"
                    "   geometry Cylinder {\n"
                    "    radius %8.6f\n"
                    "    height %8.6f\n"
                    "    bottom FALSE\n"
                    "    top    FALSE\n"
                    "   }\n"
                    "   appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "   }\n",
                    prim->r1, prim->l1, (prim->c1[0] + prim->c2[0]) / 2,
                    (prim->c1[1] + prim->c2[1]) / 2,
                    (prim->c1[2] + prim->c2[2]) / 2);
            /* WLD: format string split to comply with ISO C89 standards */
            sprintf(geom_add,
                    "  }\n"
                    "  Transform {\n"
                    "   translation 0.0 %8.6f 0.0\n"
                    "   children Shape {\n"
                    "    geometry Sphere { radius %8.6f }\n"
                    "    appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "    }\n"
                    "   }\n"
                    "  }\n",
                    prim->l1 / 2, prim->r1, prim->c1[0], prim->c1[1],
                    prim->c1[2]);
            strcat(geometry, geom_add);
            /* WLD: format string split to comply with ISO C89 standards */
            sprintf(geom_add,
                    "  Transform {\n"
                    "   translation 0.0 %8.6f 0.0\n"
                    "   children Shape {\n"
                    "    geometry Sphere { radius %8.6f }\n"
                    "    appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "    }\n"
                    "   }\n"
                    "  }\n",
                    -prim->l1 / 2, prim->r1, prim->c2[0], prim->c2[1],
                    prim->c2[2]);
            strcat(geometry, geom_add);
          } else {
            sprintf(geometry,
                    "  Shape {\n"
                    "   geometry Cylinder {\n"
                    "    radius %8.6f\n"
                    "    height %8.6f\n"
                    "   }\n"
                    "   appearance Appearance {\n"
                    "   material Material { diffuseColor %6.4f %6.4f %6.4f \n"
                    "                       specularColor 0.8 0.8 0.8 \n"
                    "                       shininess 0.8 }\n"
                    "   }\n"
                    "  }\n",
                    prim->r1, prim->l1, (prim->c1[0] + prim->c2[0]) / 2,
                    (prim->c1[1] + prim->c2[1]) / 2,
                    (prim->c1[2] + prim->c2[2]) / 2);
          }
          sprintf(buffer,
                  "Transform {\n"
                  " translation %8.6f %8.6f %8.6f\n"
                  " rotation %8.6f %8.6f %8.6f %8.6f\n"
                  " children [\n"
                  "%s"
                  " ]\n"
                  "}\n",
                  vert2[0] - mid[0], vert2[1] - mid[1], vert2[2] - mid[2],
                  axis[0], axis[1], axis[2], angle, geometry);
          UtilConcatVLA(&vla, &cc, buffer);
        } break;
        case cPrimTriangle:
          /* output coords. connectivity and vertex colors handled above/below
           */
          sprintf(buffer,
                  "%8.6f %8.6f %8.6f,\n"
                  "%8.6f %8.6f %8.6f,\n"
                  "%8.6f %8.6f %8.6f,\n",
                  vert[0] - mid[0], vert[1] - mid[1], vert[2] - mid[2],
                  vert[3] - mid[0], vert[4] - mid[1], vert[5] - mid[2],
                  vert[6] - mid[0], vert[7] - mid[1], vert[8] - mid[2]);
          UtilConcatVLA(&vla, &cc, buffer);
          break;
      }
    }

    if (mesh_obj) {
      CBasis* base = I->Basis + 1;
      CPrimitive* cprim;
      int tri = 0;
      /* output connectivity */
      UtilConcatVLA(&vla, &cc,
                    "   ]\n"
                    "  }\n"
                    "  coordIndex [\n");
      for (b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        if (TriangleReverse(cprim))
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
        else
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
        UtilConcatVLA(&vla, &cc, buffer);
        tri += 3;
      }

      /* output vertex colors */
      UtilConcatVLA(&vla, &cc,
                    "  ]\n"
                    "  colorPerVertex TRUE\n"
                    "  color Color {\n"
                    "   color [\n");
      for (b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        sprintf(buffer,
                "%6.4f %6.4f %6.4f,\n"
                "%6.4f %6.4f %6.4f,\n"
                "%6.4f %6.4f %6.4f,\n",
                cprim->c1[0], cprim->c1[1], cprim->c1[2], cprim->c2[0],
                cprim->c2[1], cprim->c2[2], cprim->c3[0], cprim->c3[1],
                cprim->c3[2]);
        UtilConcatVLA(&vla, &cc, buffer);
      }

      /* output vertex normals */
      UtilConcatVLA(&vla, &cc,
                    "  ] } \n"
                    "  normalPerVertex TRUE\n"
                    "  normal Normal {\n"
                    "   vector [\n");
      for (b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        {
          float* norm = base->Normal + 3 * base->Vert2Normal[cprim->vert];
          sprintf(buffer,
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n"
                  "%6.4f %6.4f %6.4f,\n",
                  norm[3], norm[4], norm[5],    /* transformed cprim->n1 */
                  norm[6], norm[7], norm[8],    /* transformed cprim->n2 */
                  norm[9], norm[10], norm[11]); /* transformed cprim->n3 */
          UtilConcatVLA(&vla, &cc, buffer);
        }
      }
      UtilConcatVLA(&vla, &cc,
                    "  ] }\n"
                    "  normalIndex [ \n");
      tri = 0;
      for (b = mesh_start; b < a; b++) {
        cprim = I->Primitive + b;
        if (TriangleReverse(cprim))
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 2, tri + 1);
        else
          sprintf(buffer, "%d %d %d -1,\n", tri, tri + 1, tri + 2);
        UtilConcatVLA(&vla, &cc, buffer);
        tri += 3;
      }

      /* close mesh */
      UtilConcatVLA(&vla, &cc,
                    " ] \n"
                    " }\n"
                    "}\n");
      mesh_obj = false;
    }
  }

  *vla_ptr = vla;
}

/* simple write-once/read-many hash for float-3/4 vectors  */

#define VECTOR_HASH_MASK 0xFFFF

typedef struct {
  float key[4];
  int value;
  int next; /* 1-based offsets, 0 = terminal */
} VectorHashElem;

typedef struct {
  int first[VECTOR_HASH_MASK + 1];
  VectorHashElem* elem;
  int size;
} VectorHash;

static void VectorHash_Free(VectorHash* I) {
  if (I) {
    VLAFreeP(I->elem);
  }
  FreeP(I);
}

static VectorHash* VectorHash_New(void) {
  VectorHash* I = pymol::calloc<VectorHash>(1);
  if (I) {
    I->elem = VLACalloc(VectorHashElem, 100);
    if (!I->elem) {
      VectorHash_Free(I);
      I = nullptr;
    }
  }
  return I;
}

static int VectorHash_GetOrSetKeyValue(VectorHash* I, float* key, float* alpha,
                                       int* value) {
  unsigned int hash;
  /* returns non-zero if the entry is new */
  {
    unsigned int a, b, c;

    a = ((unsigned int*)key)[0];
    b = ((unsigned int*)key)[1];
    c = ((unsigned int*)key)[2];

    /* Robert Jenkin's 96 to 32 bit hash (public domain)
       this is probably way overkill */

    a = a - b;
    a = a - c;
    a = a ^ (c >> 13);
    b = b - c;
    b = b - a;
    b = b ^ (a << 8);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 13);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 12);
    b = b - c;
    b = b - a;
    b = b ^ (a << 16);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 5);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 3);
    b = b - c;
    b = b - a;
    b = b ^ (a << 10);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 15);

    /* mix in the fourth key (if present) */

    if (alpha) {
      unsigned int d = *((unsigned int*)alpha);
      c = c + d;
    }

    /* fold those 32 bits to 16 */

    c ^= (c >> 16);

    /* apply mask */

    hash = c & VECTOR_HASH_MASK;
  }
  {
    int offset = I->first[hash];
    while (offset) {
      VectorHashElem* elem = I->elem + offset;
      float* v = elem->key;
      if ((v[0] == key[0]) && (v[1] == key[1]) && (v[2] == key[2])) {
        if ((!alpha) || (*alpha == v[3])) {
          *value = elem->value; /* matched, so return value */
          return 0;             /* key/value exists */
        }
      }
      offset = elem->next;
    }
    /* not matched -- add new key/value  */
    if (VLACheck(I->elem, VectorHashElem, ++(I->size))) {
      VectorHashElem* elem = I->elem + I->size;
      elem->next = I->first[hash];
      I->first[hash] = I->size;
      copy3f(key, elem->key);
      if (alpha) elem->key[3] = *alpha;
      elem->value = *value;
      return 1; /* inform caller */
    } else {
      I->size--;
      return -1;
    }
  }
}

static void unique_vector_add(VectorHash* vh, float* vector,
                              float* vector_array, int* vector_count,
                              int* index_array, int* index_count) {
  int index = *vector_count;
  switch (VectorHash_GetOrSetKeyValue(vh, vector, nullptr, &index)) {
    case 1: {
      float* vector_slot = vector_array + 3 * (*vector_count);
      copy3f(vector, vector_slot);
      (*vector_count)++;
    }
      /* INTENTIONAL omission of break */
    case 0:
      index_array[(*index_count)++] = index;
      break;
  }
}

static void unique_color_add(VectorHash* vh, float* vector, float* vector_array,
                             int* vector_count, int* index_array,
                             int* index_count, float alpha) {
  int index = *vector_count;

  switch (VectorHash_GetOrSetKeyValue(vh, vector, &alpha, &index)) {
    case 1: {
      float* vector_slot =
          vector_array + 4 * (*vector_count); /* NOTE 4x spacing */
      copy3f(vector, vector_slot);
      vector_slot[3] = alpha;
      (*vector_count)++;
    }
      /* INTENTIONAL omission of break */
    case 0:
      index_array[(*index_count)++] = index;
      break;
  }
}

typedef struct {
  int face_count;
  int position_count;
  int normal_count;
  int* face_position_list;
  int* face_normal_list;
  int* face_shading_list;
  float* model_position_list;
  float* model_normal_list;
  VectorHash* position_hash;
  VectorHash* normal_hash;
  int color_count;
  int* face_color_list;
  float* model_diffuse_color_list;
  VectorHash* color_hash;
} IdtfResourceMesh;

typedef struct {
  float* color_list;
  int color_count;
  VectorHash* color_hash;
} IdtfMaterial;

static ov_size idtf_dump_file_header(char** vla, ov_size cnt) {
  UtilConcatVLA(vla, &cnt, "FILE_FORMAT \"IDTF\"\nFORMAT_VERSION 100\n\n");

  UtilConcatVLA(vla, &cnt, "NODE \"VIEW\" {\n");
  UtilConcatVLA(vla, &cnt, "\tNODE_NAME \"DefaultView\"\n");
  UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n\t\t\tPARENT_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t1.000000 0.000000 0.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 1.000000 0.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 1.000000 0.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 0.000000 1.0\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_NAME \"SceneViewResource\"\n");
  UtilConcatVLA(vla, &cnt, "\tVIEW_DATA {\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_TYPE \"PERSPECTIVE\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_PROJECTION 34.515877\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  UtilConcatVLA(vla, &cnt, "NODE \"LIGHT\"\n");
  UtilConcatVLA(vla, &cnt, "{\n");
  UtilConcatVLA(vla, &cnt, "\tNODE_NAME \"Omni01\"\n");
  UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t1.000000 0.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 1.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t0.000000 0.000000 1.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\t50.000000 -50.00000 50.000000 1.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_NAME \"DefaultPointLight\"\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");
  return cnt;
}

static ov_size idtf_dump_model_nodes(char** vla, ov_size cnt,
                                     IdtfResourceMesh* mesh_vla, int n_mesh) {
  int a;
  IdtfResourceMesh* mesh = mesh_vla;
  for (a = 0; a < n_mesh; a++) {
    OrthoLineType buffer;

    UtilConcatVLA(vla, &cnt, "NODE \"MODEL\" {\n");

    sprintf(buffer, "\tNODE_NAME \"Mesh%d\"\n", a);
    UtilConcatVLA(vla, &cnt, buffer);

    UtilConcatVLA(vla, &cnt, "\tPARENT_LIST {\n");
    UtilConcatVLA(vla, &cnt, "\t\tPARENT_COUNT 1\n");
    UtilConcatVLA(vla, &cnt, "\t\tPARENT 0 {\n");
    UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_NAME \"<NULL>\"\n");
    UtilConcatVLA(vla, &cnt, "\t\t\tPARENT_TM {\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t1.000000 0.000000 0.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 1.000000 0.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 0.000000 1.000000 0.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t0.000000 0.000000 0.000000 1.0\n");
    UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
    UtilConcatVLA(vla, &cnt, "\t\t}\n");
    UtilConcatVLA(vla, &cnt, "\t}\n");

    sprintf(buffer, "\tRESOURCE_NAME \"Mesh%d\"\n", a);
    UtilConcatVLA(vla, &cnt, buffer);

    UtilConcatVLA(vla, &cnt, "}\n\n");

    mesh++;
  }
  return cnt;
}

static ov_size idtf_dump_resource_header(char** vla, ov_size cnt) {
  UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"VIEW\" {\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\tRESOURCE_NAME \"SceneViewResource\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_PASS_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\t\tVIEW_ROOT_NODE_LIST {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\tROOT_NODE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t\tROOT_NODE_NAME \"<NULL>\"\n");
  UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t\t}\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"LIGHT\" {\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE_COUNT 1\n");
  UtilConcatVLA(vla, &cnt, "\tRESOURCE 0 {\n");
  UtilConcatVLA(vla, &cnt, "\t\tRESOURCE_NAME \"DefaultPointLight\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_TYPE \"POINT\"\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_COLOR 1.000000 1.000000 1.000000\n");
  UtilConcatVLA(vla, &cnt,
                "\t\tLIGHT_ATTENUATION 1.000000 0.000000 0.000000\n");
  UtilConcatVLA(vla, &cnt, "\t\tLIGHT_INTENSITY 1.000000\n");
  UtilConcatVLA(vla, &cnt, "\t}\n");
  UtilConcatVLA(vla, &cnt, "}\n\n");

  return cnt;
}

static ov_size idtf_dump_resources(char** vla, ov_size cnt,
                                   IdtfResourceMesh* mesh_vla, int n_mesh,
                                   IdtfMaterial* material) {
  {
    OrthoLineType buffer;
    int n_color = material->color_count;

    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"SHADER\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_color);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int c;
      for (c = 0; c < n_color; c++) {
        sprintf(buffer, "\tRESOURCE %d {\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tRESOURCE_NAME \"Shader%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tSHADER_MATERIAL_NAME \"Material%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t\tSHADER_ACTIVE_TEXTURE_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t}\n");
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }

  {
    OrthoLineType buffer;
    int n_color = material->color_count;

    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"MATERIAL\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_color);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int c;
      float* fp = material->color_list;

      for (c = 0; c < n_color; c++) {
        sprintf(buffer, "\tRESOURCE %d {\n", c);
        UtilConcatVLA(vla, &cnt, buffer);
        sprintf(buffer, "\t\tRESOURCE_NAME \"Material%06d\"\n", c);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tMATERIAL_AMBIENT %0.6f %0.6f %0.6f\n", fp[0] * 0,
                fp[1] * 0, fp[2] * 0);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tMATERIAL_DIFFUSE %0.6f %0.6f %0.6f\n", fp[0],
                fp[1], fp[2]);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt,
                      "\t\tMATERIAL_SPECULAR 0.750000 0.750000 0.750000\n");

        sprintf(buffer, "\t\tMATERIAL_EMISSIVE %0.6f %0.6f %0.6f\n",
                fp[0] * 0.13, fp[1] * 0.13, fp[2] * 0.13);
        UtilConcatVLA(vla, &cnt, buffer);
        UtilConcatVLA(vla, &cnt, "\t\tMATERIAL_REFLECTIVITY 0.40000\n");

        sprintf(buffer, "\t\tMATERIAL_OPACITY %0.6f\n", fp[3]);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t}\n");

        fp += 4;
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }

  {
    OrthoLineType buffer;
    UtilConcatVLA(vla, &cnt, "RESOURCE_LIST \"MODEL\" {\n");

    sprintf(buffer, "\tRESOURCE_COUNT %d\n", n_mesh);
    UtilConcatVLA(vla, &cnt, buffer);

    {
      int a;
      IdtfResourceMesh* mesh = mesh_vla;
      for (a = 0; a < n_mesh; a++) {
        sprintf(buffer, "\tRESOURCE %d {\n", a);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\tRESOURCE_NAME \"Mesh%d\"\n", a);
        UtilConcatVLA(vla, &cnt, buffer);

        UtilConcatVLA(vla, &cnt, "\t\tMODEL_TYPE \"MESH\"\n");
        UtilConcatVLA(vla, &cnt, "\t\tMESH {\n");

        sprintf(buffer, "\t\t\tFACE_COUNT %d\n", mesh->face_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_POSITION_COUNT %d\n",
                mesh->position_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_NORMAL_COUNT %d\n", mesh->normal_count);
        UtilConcatVLA(vla, &cnt, buffer);
#ifdef IDTF_COLOR
        sprintf(buffer, "\t\t\tMODEL_DIFFUSE_COLOR_COUNT %d\n",
                mesh->color_count);
        UtilConcatVLA(vla, &cnt, buffer);

        sprintf(buffer, "\t\t\tMODEL_SPECULAR_COLOR_COUNT %d\n",
                mesh->color_count);
        UtilConcatVLA(vla, &cnt, buffer);
#else
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_DIFFUSE_COLOR_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SPECULAR_COLOR_COUNT 0\n");
#endif
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_TEXTURE_COORD_COUNT 0\n");
        UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_BONE_COUNT 0\n");

        {
          int n_color = material->color_count;

          sprintf(buffer, "\t\t\tMODEL_SHADING_COUNT %d\n", n_color);
          UtilConcatVLA(vla, &cnt, buffer);

          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SHADING_DESCRIPTION_LIST {\n");

          {
            int c;
            for (c = 0; c < n_color; c++) {
              sprintf(buffer, "\t\t\t\tSHADING_DESCRIPTION %d {\n", c);
              UtilConcatVLA(vla, &cnt, buffer);

              UtilConcatVLA(vla, &cnt, "\t\t\t\tTEXTURE_LAYER_COUNT 0\n");

              sprintf(buffer, "\t\t\t\tSHADER_ID %d\n", c + 1);
              UtilConcatVLA(vla, &cnt, buffer);

              UtilConcatVLA(vla, &cnt, "\t\t\t\t}\n");
            }
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int* ip = mesh->face_position_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_POSITION_LIST {\n");

          for (b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int* ip = mesh->face_normal_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_NORMAL_LIST {\n");

          for (b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          int* ip = mesh->face_shading_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_SHADING_LIST {\n");

          for (b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d\n", ip[0]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip++;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

#ifdef IDTF_COLOR
        {
          int b;
          int* ip = mesh->face_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_DIFFUSE_COLOR_LIST {\n");

          for (b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
        {
          int b;
          int* ip = mesh->face_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMESH_FACE_SPECULAR_COLOR_LIST {\n");

          for (b = 0; b < mesh->face_count; b++) {
            sprintf(buffer, "\t\t\t%d %d %d\n", ip[0], ip[1], ip[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            ip += 3;
          }
          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#endif

        {
          int b;
          float* fp = mesh->model_position_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_POSITION_LIST {\n");

          for (b = 0; b < mesh->position_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 3;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }

        {
          int b;
          float* fp = mesh->model_normal_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_NORMAL_LIST {\n");

          for (b = 0; b < mesh->normal_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f\n", fp[0], fp[1], fp[2]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 3;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#ifdef IDTF_COLOR
        {
          int b;
          float* fp = mesh->model_diffuse_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_DIFFUSE_COLOR_LIST {\n");

          for (b = 0; b < mesh->color_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f %1.6f\n", fp[0], fp[1],
                    fp[2], fp[3]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 4;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
        {
          int b;
          float* fp = mesh->model_diffuse_color_list;
          UtilConcatVLA(vla, &cnt, "\t\t\tMODEL_SPECULAR_COLOR_LIST {\n");

          for (b = 0; b < mesh->color_count; b++) {
            sprintf(buffer, "\t\t\t\t%1.6f %1.6f %1.6f %1.6f\n", fp[0], fp[1],
                    fp[2], fp[3]);
            UtilConcatVLA(vla, &cnt, buffer);
            fp += 4;
          }

          UtilConcatVLA(vla, &cnt, "\t\t\t}\n");
        }
#endif

        UtilConcatVLA(vla, &cnt, "\t\t}\n");
        UtilConcatVLA(vla, &cnt, "\t}\n");

        mesh++;
      }
    }
    UtilConcatVLA(vla, &cnt, "}\n\n");
  }
  return cnt;
}
void RayRenderIDTF(CRay* I, char** node_vla, char** rsrc_vla) {
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  {
    CBasis* base = I->Basis + 1;
    CPrimitive* prim = I->Primitive;
    int mesh_cnt = 0;
    IdtfResourceMesh* mesh_vla = VLACalloc(IdtfResourceMesh, 1);
    if (mesh_vla) {
      IdtfResourceMesh* mesh = nullptr;
      int a;

      for (a = 0; a < I->NPrimitive; a++) {
        switch (prim->type) {
          case cPrimTriangle:
          case cPrimSphere:
            if (!mesh) {
              /* create a new triangle mesh */
              if (VLACheck(mesh_vla, IdtfResourceMesh, mesh_cnt)) {
                mesh = mesh_vla + mesh_cnt;
                if ((mesh->face_position_list = VLACalloc(int, 3)) &&
                    (mesh->face_normal_list = VLACalloc(int, 3)) &&
                    (mesh->face_shading_list =
                         VLACalloc(int, 1)) && /* defaults to zero */

                    (mesh->model_position_list = VLAlloc(float, 3)) &&
                    (mesh->face_color_list = VLACalloc(int, 3)) &&
                    (mesh->model_diffuse_color_list = VLAlloc(float, 4)) &&
                    ((mesh->color_hash = VectorHash_New())) &&
                    (mesh->model_normal_list = VLAlloc(float, 3)) &&
                    ((mesh->position_hash = VectorHash_New())) &&
                    ((mesh->normal_hash = VectorHash_New()))) {
                  mesh_cnt++;
                } else {
                  mesh = nullptr;
                }
              }
            }
            break;
          default: /* close/terminate mesh */
            if (mesh) {
              mesh = nullptr;
            }
            break;
        }

        switch (prim->type) {
          case cPrimTriangle:
            if (mesh) {
              if (VLACheck(mesh->face_position_list, int,
                           mesh->face_count * 3 + 2) &&
                  VLACheck(mesh->face_normal_list, int,
                           mesh->face_count * 3 + 2) &&
                  VLACheck(mesh->face_shading_list, int, mesh->face_count) &&
                  VLACheck(mesh->model_position_list, float,
                           (mesh->position_count + 3) * 3) &&
                  VLACheck(mesh->face_color_list, int,
                           mesh->face_count * 3 + 2) &&
                  VLACheck(mesh->model_diffuse_color_list, float,
                           (mesh->color_count + 3) * 4) &&
                  VLACheck(mesh->model_normal_list, float,
                           (mesh->normal_count + 3) * 3)) {
                float* vert = base->Vertex + 3 * (prim->vert);
                float* norm =
                    base->Normal + 3 * base->Vert2Normal[prim->vert] + 3;
                int reverse = TriangleReverse(prim);
                int face_position_count = mesh->face_count * 3;
                int face_normal_count = face_position_count;
                int face_color_count = face_position_count;

                unique_vector_add(
                    mesh->position_hash, vert, mesh->model_position_list,
                    &mesh->position_count, mesh->face_position_list,
                    &face_position_count);
                unique_vector_add(mesh->normal_hash, norm,
                                  mesh->model_normal_list, &mesh->normal_count,
                                  mesh->face_normal_list, &face_normal_count);
                unique_color_add(mesh->normal_hash, prim->c1,
                                 mesh->model_diffuse_color_list,
                                 &mesh->color_count, mesh->face_color_list,
                                 &face_color_count, 1.0F - prim->trans);
                if (reverse) {
                  vert += 6;
                  norm += 6;
                  unique_vector_add(
                      mesh->position_hash, vert, mesh->model_position_list,
                      &mesh->position_count, mesh->face_position_list,
                      &face_position_count);
                  unique_vector_add(mesh->normal_hash, norm,
                                    mesh->model_normal_list,
                                    &mesh->normal_count, mesh->face_normal_list,
                                    &face_normal_count);
                  unique_color_add(mesh->normal_hash, prim->c3,
                                   mesh->model_diffuse_color_list,
                                   &mesh->color_count, mesh->face_color_list,
                                   &face_color_count, 1.0F - prim->trans);
                  vert -= 3;
                  norm -= 3;
                  unique_vector_add(
                      mesh->position_hash, vert, mesh->model_position_list,
                      &mesh->position_count, mesh->face_position_list,
                      &face_position_count);
                  unique_vector_add(mesh->normal_hash, norm,
                                    mesh->model_normal_list,
                                    &mesh->normal_count, mesh->face_normal_list,
                                    &face_normal_count);
                  unique_color_add(mesh->normal_hash, prim->c2,
                                   mesh->model_diffuse_color_list,
                                   &mesh->color_count, mesh->face_color_list,
                                   &face_color_count, 1.0F - prim->trans);
                } else {
                  vert += 3;
                  norm += 3;
                  unique_vector_add(
                      mesh->position_hash, vert, mesh->model_position_list,
                      &mesh->position_count, mesh->face_position_list,
                      &face_position_count);
                  unique_vector_add(mesh->normal_hash, norm,
                                    mesh->model_normal_list,
                                    &mesh->normal_count, mesh->face_normal_list,
                                    &face_normal_count);
                  unique_color_add(mesh->normal_hash, prim->c2,
                                   mesh->model_diffuse_color_list,
                                   &mesh->color_count, mesh->face_color_list,
                                   &face_color_count, 1.0F - prim->trans);
                  vert += 3;
                  norm += 3;
                  unique_vector_add(
                      mesh->position_hash, vert, mesh->model_position_list,
                      &mesh->position_count, mesh->face_position_list,
                      &face_position_count);
                  unique_vector_add(mesh->normal_hash, norm,
                                    mesh->model_normal_list,
                                    &mesh->normal_count, mesh->face_normal_list,
                                    &face_normal_count);
                  unique_color_add(mesh->normal_hash, prim->c3,
                                   mesh->model_diffuse_color_list,
                                   &mesh->color_count, mesh->face_color_list,
                                   &face_color_count, 1.0F - prim->trans);
                }
                mesh->face_count++;
              }
            }
            break;
          case cPrimSphere:

            break;
          case cPrimCone:
            break;
          case cPrimCylinder:
          case cPrimSausage:
            break;
        }

        /* looping through each primitive */
        prim++;
      }

      /* now we need to consolidate materials for each color
         combination (creating averages as necessary) and update
         mesh->face_shading_list appropriately for each face */

      {
        IdtfMaterial* material = pymol::calloc<IdtfMaterial>(1);
        if (material && (material->color_list = VLAlloc(float, 4)) &&
            (material->color_hash = VectorHash_New())) {
          const float one_third = (1 / 3.0F);
          int a;
          IdtfResourceMesh* mesh = mesh_vla;

          for (a = 0; a < mesh_cnt; a++) {
            int* ip = mesh->face_color_list;
            int shading_count = 0;
            int b;
            for (b = 0; b < mesh->face_count; b++) {
              float* fp0 = mesh->model_diffuse_color_list + 4 * ip[0];
              float* fp1 = mesh->model_diffuse_color_list + 4 * ip[1];
              float* fp2 = mesh->model_diffuse_color_list + 4 * ip[2];
              float avg_rgba[4];

              avg_rgba[0] = (fp0[0] + fp1[0] + fp2[0]) * one_third;
              avg_rgba[1] = (fp0[1] + fp1[1] + fp2[1]) * one_third;
              avg_rgba[2] = (fp0[2] + fp1[2] + fp2[2]) * one_third;
              avg_rgba[3] = (fp0[3] + fp1[3] + fp2[3]) * one_third;

              if (VLACheck(material->color_list, float,
                           material->color_count * 4 + 3)) {
                unique_color_add(material->color_hash, avg_rgba,
                                 material->color_list, &material->color_count,
                                 mesh->face_shading_list, &shading_count,
                                 avg_rgba[3]);
              }
              ip += 3;
            }
            mesh++;
          }

          {
            int cnt = 0;
            cnt = idtf_dump_file_header(node_vla, cnt);
            cnt = idtf_dump_model_nodes(node_vla, cnt, mesh_vla, mesh_cnt);
            VLASize((*node_vla), char, cnt + 1);
          }
          {
            int cnt = 0;
            cnt = idtf_dump_resource_header(rsrc_vla, cnt);
            cnt = idtf_dump_resources(rsrc_vla, cnt, mesh_vla, mesh_cnt,
                                      material);
            VLASize((*rsrc_vla), char, cnt + 1);
          }

          VLAFreeP(material->color_list);
          VectorHash_Free(material->color_hash);
        }
        FreeP(material);
      }

      {
        /* refactor into a delete method */
        IdtfResourceMesh* mesh = mesh_vla;
        int i;
        for (i = 0; i < mesh_cnt; i++) {
          VLAFreeP(mesh->face_position_list);
          VLAFreeP(mesh->face_normal_list);
          VLAFreeP(mesh->face_shading_list);
          VLAFreeP(mesh->face_color_list);
          VLAFreeP(mesh->model_position_list);
          VLAFreeP(mesh->model_normal_list);
          VLAFreeP(mesh->model_diffuse_color_list);
          VectorHash_Free(mesh->color_hash);
          VectorHash_Free(mesh->position_hash);
          VectorHash_Free(mesh->normal_hash);
          mesh++;
        }
        VLAFreeP(mesh_vla);
      }
    }
  }
}
void RayRenderObjMtl(CRay* I, int width, int height, char** objVLA_ptr,
                     char** mtlVLA_ptr, float front, float back, float fov,
                     float angle, float z_corr) {
  char* objVLA = *objVLA_ptr;
  char* mtlVLA = *mtlVLA_ptr;
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  ov_size oc = 0;    /* obj character count */
  /*  int mc = 0; */ /* mtl character count */

  OrthoLineType buffer;

  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  {
    int a;
    CPrimitive* prim;
    float *vert, *norm;
    int vc = 0;
    int nc = 0;
    CBasis* base = I->Basis + 1;

    for (a = 0; a < I->NPrimitive; a++) {
      prim = I->Primitive + a;
      vert = base->Vertex + 3 * (prim->vert);
      norm = base->Normal + 3 * base->Vert2Normal[prim->vert] + 3;
      switch (prim->type) {
        case cPrimTriangle:
          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1],
                  vert[2] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[3], vert[4],
                  vert[5] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[6], vert[7],
                  vert[8] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[0], norm[1], norm[2]);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[3], norm[4], norm[5]);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "vn %8.6f %8.6f %8.6f\n", norm[6], norm[7], norm[8]);
          UtilConcatVLA(&objVLA, &oc, buffer);
          if (TriangleReverse(prim)) {
            sprintf(buffer, "f %d//%d %d//%d %d//%d\n", vc + 1, nc + 1, vc + 3,
                    nc + 3, vc + 2, nc + 2);
          } else {
            sprintf(buffer, "f %d//%d %d//%d %d//%d\n", vc + 1, nc + 1, vc + 2,
                    nc + 2, vc + 3, nc + 3);
          }
          UtilConcatVLA(&objVLA, &oc, buffer);
          nc += 3;
          vc += 3;

          /*
             prim->c1[0],prim->c1[1],prim->c1[2])
             prim->c2[0],prim->c2[1],prim->c2[2],
             prim->c3[0],prim->c3[1],prim->c3[2]
             UtilConcatVLA(&vla,&oc,buffer);
             UtilConcatVLA(&vla,&oc,buffer);
           */

          break;
        case cPrimSphere:
          /*      sprintf(buffer,"v %8.6f %8.6f %8.6f\np %d\n",
             vert[0],vert[1],vert[2]-z_corr,vc+1);
             UtilConcatVLA(&objVLA,&oc,buffer); */

          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1],
                  vert[2] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1],
                  vert[2] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "v %8.6f %8.6f %8.6f\n", vert[0], vert[1],
                  vert[2] - z_corr);
          UtilConcatVLA(&objVLA, &oc, buffer);
          sprintf(buffer, "f %d %d %d\n", vc + 1, vc + 2, vc + 3);
          UtilConcatVLA(&objVLA, &oc, buffer);
          vc += 3;
          break;
      }
    }
  }

  *objVLA_ptr = objVLA;
  *mtlVLA_ptr = mtlVLA;
}
void RayRenderPOV(CRay* I, int width, int height, char** headerVLA_ptr,
                  char** charVLA_ptr, float front, float back, float fov,
                  float angle, int antialias) {
  float fog;
  const float* bkrd;
  float fog_start = 0.0F;
  float gamma;
  float* d;
  CBasis* base;
  CPrimitive* prim;
  OrthoLineType buffer;
  float *vert, *norm;
  float vert2[3];
  ov_size cc, hc;
  int a;
  int smooth_color_triangle;
  int mesh_obj = false;
  char *charVLA, *headerVLA;
  char transmit[64];
  float light[3];
  const float* lightv;
  float spec_power = SettingGetGlobal_f(I->G, cSetting_spec_power);
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  if (spec_power < 0.0F) {
    spec_power = SettingGetGlobal_f(I->G, cSetting_shininess);
  }
  spec_power /= 4.0F;

  charVLA = *charVLA_ptr;
  headerVLA = *headerVLA_ptr;
  smooth_color_triangle =
      SettingGetGlobal_b(I->G, cSetting_smooth_color_triangle);
  PRINTFB(I->G, FB_Ray, FB_Blather)
  " RayRenderPOV: w %d h %d f %8.3f b %8.3f\n", width, height, front,
      back ENDFB(I->G);
  if (Feedback(I->G, FB_Ray, FB_Blather)) {
    dump3f(I->Volume, " RayRenderPOV: vol");
    dump3f(I->Volume + 3, " RayRenderPOV: vol");
  }
  cc = 0;
  hc = 0;
  gamma = SettingGetGlobal_f(I->G, cSetting_gamma);
  if (gamma > R_SMALL4)
    gamma = 1.0F / gamma;
  else
    gamma = 1.0F;

  lightv = SettingGetfv(I->G, cSetting_light);
  copy3f(lightv, light);

  fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
  if (fog < 0.0F) fog = SettingGetGlobal_b(I->G, cSetting_depth_cue);
  if (fog != 0.0F) {
    if (fog > 1.0F) fog = 1.0F;
    fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
    if (fog_start < 0.0F)
      fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
    if (fog_start > 1.0F) fog_start = 1.0F;
  }

  /* SETUP */

  if (antialias < 0) antialias = SettingGetGlobal_i(I->G, cSetting_antialias);

  bkrd =
      ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);

  PRINTFB(I->G, FB_Ray, FB_Details)
  " RayRenderPovRay: processed %i graphics primitives.\n",
      I->NPrimitive ENDFB(I->G);
  base = I->Basis + 1;

  {
    int ortho = SettingGetGlobal_b(I->G, cSetting_ortho);
    if (!identity) {
      if (!ortho) {
        sprintf(buffer,
                "camera {direction<0.0,0.0,%8.3f>\n location <0.0 , 0.0 , "
                "0.0>\n right %12.10f*x up y \n }\n",
                -57.3F * cos(fov * cPI / (180 * 2.4)) /
                    fov, /* by trial and error */
                I->Range[0] / I->Range[1]);
      } else {
        sprintf(buffer,
                "camera {orthographic location <0.0 , 0.0 , %12.10f>\nlook_at  "
                "<0.0 , 0.0 , -1.0> right %12.10f*x up %12.10f*y}\n",
                front, -I->Range[0], I->Range[1]);
      }
    } else {
      float look[3];
      float loc[3] = {0.0F, 0.0F, 0.0F};
      SceneViewType view;

      SceneGetCenter(I->G, look);
      SceneGetView(I->G, view);
      loc[2] = -view[18];
      MatrixInvTransformC44fAs33f3f(view, loc, loc);
      add3f(look, loc, loc);

      if (!ortho) {
        sprintf(
            buffer,
            "camera {angle %12.10f "
            "sky<%12.10f,%12.10f,%12.10f>\nlocation<%12.10f,%12.10f,%12.10f>"
            "\nlook_at<%12.10f,%12.10f,%12.10f> right %12.10f*x up y }\n",
            fov * I->Range[0] / I->Range[1], view[1], view[5], view[9], loc[0],
            loc[1], loc[2], look[0], look[1], look[2],
            -I->Range[0] / I->Range[1]);
      } else {
        sprintf(buffer,
                "camera {orthographic "
                "sky<%12.10f,%12.10f,%12.10f>\nlocation<%12.10f,%12.10f,%12."
                "10f>\nlook_at<%12.10f,%12.10f,%12.10f> right %12.10f*x up "
                "%12.10f*y}\n",
                view[1], view[5], view[9], loc[0], loc[1], loc[2], look[0],
                look[1], look[2], -I->Range[0], I->Range[1]);
      }
    }
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  {
    float ambient = SettingGetGlobal_f(I->G, cSetting_ambient) +
                    SettingGetGlobal_f(I->G, cSetting_direct);
    float reflect = SettingGetGlobal_f(I->G, cSetting_reflect);

    if (ambient > 0.5) ambient = 0.5;

    reflect = 1.2F - 1.5F * ambient;

    sprintf(buffer,
            "#default { finish{phong %8.3f ambient %8.3f diffuse %8.3f "
            "phong_size %8.6f}}\n",
            SettingGetGlobal_f(I->G, cSetting_spec_reflect), ambient, reflect,
            spec_power);
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  if (!identity) {
    if (angle) {
      float temp[16];
      identity44f(temp);
      MatrixRotateC44f(temp, (float)-PI * angle / 180, 0.0F, 1.0F, 0.0F);
      MatrixTransformC44fAs33f3f(temp, light, light);
    }
  }

  {
    float lite[3];
    if (!identity) {
      lite[0] = -light[0] * 10000.0F;
      lite[1] = -light[1] * 10000.0F;
      lite[2] = -light[2] * 10000.0F - front;
    } else {
      /* this is not correct unless the camera is in the default orientation */
      float look[3];
      SceneViewType view;

      SceneGetCenter(I->G, look);
      SceneGetView(I->G, view);

      lite[0] = -light[0] * 10000.0F;
      lite[1] = -light[1] * 10000.0F;
      lite[2] = -light[2] * 10000.0F;
      MatrixInvTransformC44fAs33f3f(view, lite, lite);
      add3f(look, lite, lite);
    }
    sprintf(buffer, "light_source{<%6.4f,%6.4f,%6.4f>  rgb<1.0,1.0,1.0>}\n",
            lite[0], lite[1], lite[2]);
    UtilConcatVLA(&headerVLA, &hc, buffer);
  }

  if (!identity) {
    int opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
    if (opaque_back < 0)
      opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

    if (opaque_back) { /* drop a plane into the background for the background
                          color */
      sprintf(buffer,
              "plane{z , %6.4f \n pigment{color rgb<%6.4f,%6.4f,%6.4f>}\n "
              "finish{phong 0 specular 0 diffuse 0 ambient 1.0}}\n",
              -back, bkrd[0], bkrd[1], bkrd[2]);
      UtilConcatVLA(&headerVLA, &hc, buffer);
    }
  }

  for (a = 0; a < I->NPrimitive; a++) {
    auto cap1 = cCylCap::Round;
    auto cap2 = cCylCap::Round;
    prim = I->Primitive + a;
    vert = base->Vertex + 3 * (prim->vert);
    if (prim->type == cPrimTriangle) {
      if (smooth_color_triangle)
        if (!mesh_obj) {
          sprintf(buffer, "mesh {\n");
          UtilConcatVLA(&charVLA, &cc, buffer);
          mesh_obj = true;
        }
    } else if (mesh_obj) {
      sprintf(buffer, " pigment{color rgb <1,1,1>}}");
      UtilConcatVLA(&charVLA, &cc, buffer);
      mesh_obj = false;
    }
    switch (prim->type) {
      case cPrimSphere:
        sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n", vert[0],
                vert[1], vert[2], prim->r1);
        UtilConcatVLA(&charVLA, &cc, buffer);
        sprintf(buffer, "pigment{color rgb<%6.4f,%6.4f,%6.4f>}}\n", prim->c1[0],
                prim->c1[1], prim->c1[2]);
        UtilConcatVLA(&charVLA, &cc, buffer);
        break;
      case cPrimCylinder:
        cap1 = prim->cap1;
        cap2 = prim->cap2;
      case cPrimSausage:
        d = base->Normal + 3 * base->Vert2Normal[prim->vert];
        scale3f(d, prim->l1, vert2);
        add3f(vert, vert2, vert2);
        sprintf(buffer,
                "cylinder{<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>"
                ",\n %12.10f\n",
                vert[0], vert[1], vert[2], vert2[0], vert2[1], vert2[2],
                prim->r1);
        UtilConcatVLA(&charVLA, &cc, buffer);

        if (cap1 != cCylCapFlat && cap2 != cCylCapFlat) {
          UtilConcatVLA(&charVLA, &cc, "open\n");
        }

        sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
                (prim->c1[0] + prim->c2[0]) / 2,
                (prim->c1[1] + prim->c2[1]) / 2,
                (prim->c1[2] + prim->c2[2]) / 2);
        UtilConcatVLA(&charVLA, &cc, buffer);

        if (cap1 == cCylCapRound) {
          sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
                  vert[0], vert[1], vert[2], prim->r1);
          UtilConcatVLA(&charVLA, &cc, buffer);
          sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
                  prim->c1[0], prim->c1[1], prim->c1[2]);
          UtilConcatVLA(&charVLA, &cc, buffer);
        }

        if (cap2 == cCylCapRound) {
          sprintf(buffer, "sphere{<%12.10f,%12.10f,%12.10f>, %12.10f\n",
                  vert2[0], vert2[1], vert2[2], prim->r1);
          UtilConcatVLA(&charVLA, &cc, buffer);
          sprintf(buffer, "pigment{color rgb<%6.4f1,%6.4f,%6.4f>}}\n",
                  prim->c2[0], prim->c2[1], prim->c2[2]);
          UtilConcatVLA(&charVLA, &cc, buffer);
        }

        break;
      case cPrimTriangle:
        norm = base->Normal + 3 * base->Vert2Normal[prim->vert] +
               3; /* first normal is the average */

        if (!TriangleDegenerate(vert, norm, vert + 3, norm + 3, vert + 6,
                                norm + 6)) {
          if (smooth_color_triangle) {
            sprintf(buffer,
                    "smooth_color_triangle{<%12.10f,%12.10f,%12.10f>,\n<%12."
                    "10f,%12.10f,%12.10f>,\n<%6.4f1,%6.4f,%6.4f>,\n<%12.10f,%"
                    "12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%6.4f1,%6."
                    "4f,%6.4f>,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%"
                    "12.10f>,\n<%6.4f1,%6.4f,%6.4f> }\n",
                    vert[0], vert[1], vert[2], norm[0], norm[1], norm[2],
                    prim->c1[0], prim->c1[1], prim->c1[2], vert[3], vert[4],
                    vert[5], norm[3], norm[4], norm[5], prim->c2[0],
                    prim->c2[1], prim->c2[2], vert[6], vert[7], vert[8],
                    norm[6], norm[7], norm[8], prim->c3[0], prim->c3[1],
                    prim->c3[2]);
            UtilConcatVLA(&charVLA, &cc, buffer);
          } else {
            /* nowadays we use mesh2 to generate smooth_color_triangles */

            UtilConcatVLA(&charVLA, &cc, "mesh2 { ");
            sprintf(buffer,
                    "vertex_vectors { 3, "
                    "<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,\n<%"
                    "12.10f,%12.10f,%12.10f>}\n normal_vectors { "
                    "3,\n<%12.10f,%12.10f,%12.10f>,\n<%12.10f,%12.10f,%12.10f>,"
                    "\n<%12.10f,%12.10f,%12.10f>}\n",
                    vert[0], vert[1], vert[2], vert[3], vert[4], vert[5],
                    vert[6], vert[7], vert[8], norm[0], norm[1], norm[2],
                    norm[3], norm[4], norm[5], norm[6], norm[7], norm[8]);
            UtilConcatVLA(&charVLA, &cc, buffer);

            if (prim->trans > R_SMALL4)
              sprintf(transmit, "transmit %4.6f", prim->trans);
            else
              transmit[0] = 0;

            sprintf(buffer, "texture_list { 3, ");
            UtilConcatVLA(&charVLA, &cc, buffer);

            sprintf(buffer,
                    "texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                    prim->c1[0], prim->c1[1], prim->c1[2], transmit);
            UtilConcatVLA(&charVLA, &cc, buffer);

            sprintf(buffer,
                    ",texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}}\n",
                    prim->c2[0], prim->c2[1], prim->c2[2], transmit);
            UtilConcatVLA(&charVLA, &cc, buffer);

            sprintf(buffer,
                    ",texture { pigment{color rgb<%6.4f1,%6.4f,%6.4f> %s}} }\n",
                    prim->c3[0], prim->c3[1], prim->c3[2], transmit);
            UtilConcatVLA(&charVLA, &cc, buffer);

            sprintf(buffer, "face_indices { 1, <0,1,2>, 0, 1, 2 } }\n");
            UtilConcatVLA(&charVLA, &cc, buffer);
          }
        }
        break;
    }
  }

  if (mesh_obj) {
    sprintf(buffer, " pigment{color rgb <1,1,1>}}");
    UtilConcatVLA(&charVLA, &cc, buffer);
    mesh_obj = false;
  }
  *charVLA_ptr = charVLA;
  *headerVLA_ptr = headerVLA;
}
