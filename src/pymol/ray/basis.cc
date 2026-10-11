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
-*
-*
Z* -------------------------------------------------------------------
*/

#include "pymol/platform/memory_debug.h"
#include "pymol/base/base.h"
#include "pymol/ray/basis.h"
#include "pymol/base/err.h"
#include "pymol/base/feedback.h"
#include "pymol/ray/basis_internal.h"
void BasisSetupMatrix(CBasis* I) {
  float oldZ[3] = {0.0, 0.0, 1.0};
  float newY[3];
  float dotgle, angle;

  cross_product3f(oldZ, I->LightNormal, newY);

  dotgle = dot_product3f(oldZ, I->LightNormal);

  if ((1.0 - fabs(dotgle)) < kR_SMALL4) {
    dotgle = (float)(dotgle / fabs(dotgle));
    newY[0] = 0.0;
    newY[1] = 1.0;
    newY[2] = 0.0;
  }

  normalize3f(newY);

  angle = (float)(-acos(dotgle));

  /* now all we gotta do is effect a rotation about the new Y axis to line up
   * new Z with Z */

  rotation_to_matrix33f(newY, angle, I->Matrix);

  /*
     printf("%8.3f %8.3f %8.3f
     %8.3f\n",angle*180.0/cPI,newY[0],newY[1],newY[2]);

     matrix_transform33f3f(I->Matrix,newY,test);
     printf("   %8.3f %8.3f %8.3f\n",test[0],test[1],test[2]);

     printf("   %8.3f %8.3f
     %8.3f\n",I->LightNormal[0],I->LightNormal[1],I->LightNormal[2]);
     matrix_transform33f3f(I->Matrix,I->LightNormal,test);
     printf("   %8.3f %8.3f %8.3f\n",test[0],test[1],test[2]);

     printf(">%8.3f %8.3f
     %8.3f\n",I->Matrix[0][0],I->Matrix[0][1],I->Matrix[0][2]); printf(">%8.3f
     %8.3f %8.3f\n",I->Matrix[1][0],I->Matrix[1][1],I->Matrix[1][2]);
     printf(">%8.3f %8.3f
     %8.3f\n",I->Matrix[2][0],I->Matrix[2][1],I->Matrix[2][2]);
   */
}
void BasisGetTriangleFlatDotgle(CBasis* I, RayInfo* r, int i) {
  float* n0 = I->Normal + (3 * I->Vert2Normal[i]);
  r->flat_dotgle = n0[2];
}

void BasisGetTriangleFlatDotglePerspective(CBasis* I, RayInfo* r, int i) {
  float* n0 = I->Normal + (3 * I->Vert2Normal[i]);
  r->flat_dotgle = -dot_product3f(r->dir, n0);
}
void BasisGetEllipsoidNormal(CBasis* I, RayInfo* r, int i, int perspective) {
  if (perspective) {
    r->impact[0] = r->base[0] + r->dir[0] * r->dist;
    r->impact[1] = r->base[1] + r->dir[1] * r->dist;
    r->impact[2] = r->base[2] + r->dir[2] * r->dist;
  } else {
    r->impact[0] = r->base[0];
    r->impact[1] = r->base[1];
    r->impact[2] = r->base[2] - r->dist;
  }

  {
    float* n1 = I->Normal + (3 * I->Vert2Normal[i]);
    float *n2 = n1 + 3, *n3 = n1 + 6;
    float* scale = r->prim->n0;
    float d1, d2, d3, s1, s2, s3;
    float comp1[3], comp2[3], comp3[3];
    float direct[3], surfnormal[3];

    direct[0] = r->impact[0] - r->sphere[0];
    direct[1] = r->impact[1] - r->sphere[1];
    direct[2] = r->impact[2] - r->sphere[2];

    normalize3f(direct);

    d1 = dot_product3f(direct, n1);
    d2 = dot_product3f(direct, n2);
    d3 = dot_product3f(direct, n3);

    if (scale[0] > R_SMALL8) {
      s1 = d1 / (scale[0] * scale[0]);
    } else {
      s1 = 0.0F;
    }
    if (scale[1] > R_SMALL8) {
      s2 = d2 / (scale[1] * scale[1]);
    } else {
      s2 = 0.0F;
    }
    if (scale[2] > R_SMALL8) {
      s3 = d3 / (scale[2] * scale[2]);
    } else {
      s3 = 0.0F;
    }

    scale3f(n1, s1, comp1);
    scale3f(n2, s2, comp2);
    scale3f(n3, s3, comp3);

    copy3f(comp1, surfnormal);
    add3f(comp2, surfnormal, surfnormal);
    add3f(comp3, surfnormal, surfnormal);

    normalize23f(surfnormal, r->surfnormal);
  }
}

void BasisGetTriangleNormal(CBasis* I, RayInfo* r, int i, float* fc,
                            int perspective) {
  float *n0, w2, fc0, fc1, fc2;
  float vt1[3];
  CPrimitive* lprim = r->prim;

  if (perspective) {
    r->impact[0] = r->base[0] + r->dir[0] * r->dist;
    r->impact[1] = r->base[1] + r->dir[1] * r->dist;
    r->impact[2] = r->base[2] + r->dir[2] * r->dist;
  } else {
    r->impact[0] = r->base[0];
    r->impact[1] = r->base[1];
    r->impact[2] = r->base[2] - r->dist;
  }

  n0 = I->Normal + (3 * I->Vert2Normal[i]) + 3; /* skip triangle normal */
  w2 = 1.0F - (r->tri1 + r->tri2);
  /*  printf("%8.3f %8.3f\n",r->tri[1],r->tri[2]); */

  fc0 =
      (lprim->c2[0] * r->tri1) + (lprim->c3[0] * r->tri2) + (lprim->c1[0] * w2);
  fc1 =
      (lprim->c2[1] * r->tri1) + (lprim->c3[1] * r->tri2) + (lprim->c1[1] * w2);
  fc2 =
      (lprim->c2[2] * r->tri1) + (lprim->c3[2] * r->tri2) + (lprim->c1[2] * w2);

  r->trans =
      (lprim->tr[1] * r->tri1) + (lprim->tr[2] * r->tri2) + (lprim->tr[0] * w2);

  scale3f(n0 + 3, r->tri1, r->surfnormal);
  scale3f(n0 + 6, r->tri2, vt1);
  add3f(vt1, r->surfnormal, r->surfnormal);

  scale3f(n0, w2, vt1);
  add3f(vt1, r->surfnormal, r->surfnormal);

  normalize3f(r->surfnormal);

  fc[0] = fc0;
  fc[1] = fc1;
  fc[2] = fc2;
}
int g_basis_lite_request = 0;
static void BasisRemapExpressXY(MapType* map, int* tempRef, int i_nVertex);
int BasisMakeMap(CBasis* I, int* vert2prim, CPrimitive* prim, int n_prim,
                 float* volume, int perspective, float front, float size_hint) {
  float* v;
  float ll;
  CPrimitive* prm;
  int i;
  int* tempRef = nullptr;
  int n = 0, h, q, x, y, z, j, k, l, e;
  int extra_vert = 0;
  float p[3], dd[3], *d1, *d2, vd[3], cx[3], cy[3];
  float* tempVertex = nullptr;
  float xs, ys;
  int remapMode = true; /* remap mode means that some objects will span more
                         * than one voxel, so we have to worry about populating
                         * those voxels and also about eliminating duplicates
                         * when traversing the neighbor lists */
  float min[3], max[3], extent[6];
  float sep;
  float diagonal[3];
  float l1, l2;
  float bh, ch;
  int n_voxel;
  int ok = true;
  const float _0 = 0.0;
  const float _p5 = 0.5;

  FreeP(I->LiteRef);
  I->Lite = 0;

  PRINTFD(I->G, FB_Ray)
  " BasisMakeMap: I->NVertex %d [(%8.3f, %8.3f, %8.3f),...]\n", I->NVertex,
      I->Vertex[0], I->Vertex[1], I->Vertex[2] ENDFD;

  sep = I->MinVoxel;
  if (sep == _0) {
    remapMode = false;
    sep = I->MaxRadius; /* also will imply no remapping of vertices */
  }
  /* we need to get a sense of the actual size in order to avoid sep being too
   * small */

  v = I->Vertex;

  min[0] = max[0] = v[0];
  min[1] = max[1] = v[1];
  min[2] = max[2] = v[2];

  v += 3;

  {
    int a;
    for (a = 1; a < I->NVertex; a++) {
      if (min[0] > v[0]) min[0] = v[0];
      if (max[0] < v[0]) max[0] = v[0];

      if (min[1] > v[1]) min[1] = v[1];
      if (max[1] < v[1]) max[1] = v[1];

      if (min[2] > v[2]) min[2] = v[2];
      if (max[2] < v[2]) max[2] = v[2];

      v += 3;
    }
  }

  if (volume) {
    /*
       if(min[0] > volume[0])      min[0]   = volume[0];
       if(max[0] < volume[1])      max[0]   = volume[1];
       if(min[1] > volume[2])      min[1]   = volume[2];
       if(max[1] < volume[3])      max[1]   = volume[3];
       if(min[2] > (-volume[5]))   min[2]   = (-volume[5]);
       if(max[2] < (-volume[4]))   max[2]   = (-volume[4]);
     */

    if (min[0] < volume[0]) min[0] = volume[0];
    if (max[0] > volume[1]) max[0] = volume[1];
    if (min[1] < volume[2]) min[1] = volume[2];
    if (max[1] > volume[3]) max[1] = volume[3];
    if (min[2] > (-volume[5])) min[2] = (-volume[5]);
    if (max[2] < (-volume[4])) max[2] = (-volume[4]);

    PRINTFB(I->G, FB_Ray, FB_Debugging)
    " BasisMakeMap: (%8.3f,%8.3f),(%8.3f,%8.3f),(%8.3f,%8.3f)\n", volume[0],
        volume[1], volume[2], volume[3], volume[4], volume[5] ENDFB(I->G);
  }

  /* don't break up space unnecessarily if we only have a few vertices... */

  if (I->NVertex) {
    l1 = (float)fabs(max[0] - min[0]);
    l2 = (float)fabs(max[1] - min[1]);

    if (l1 < l2) l1 = l2;

    l2 = (float)fabs(max[2] - min[2]);

    if (l1 < l2) l1 = l2;

    if (l1 < kR_SMALL4) l1 = 100.0;

    if (I->NVertex < (l1 / sep)) sep = (l1 / I->NVertex);
  }

  if (Feedback(I->G, FB_Ray, FB_Debugging)) {
    dump3f(min, " BasisMakeMap: min");
    dump3f(max, " BasisMakeMap: max");
    dump3f(I->Vertex, " BasisMakeMap: I->Vertex");
    fflush(stdout);
  }

  sep = MapGetSeparation(I->G, sep, max, min,
                         diagonal); /* this needs to be a minimum
                                     * estimate of the actual value */
  /* here we have to carry out a complicated work-around in order to
   * efficiently encode our primitives into the map in a way that doesn't
   * require expanding the map cutoff to the size of the largest object*/
  if (remapMode) {
    int a, b, c;

    if (sep < size_hint) /* this keeps us from wasting time & memory on
                            unnecessary subdivision */
      sep = size_hint;

    {
      int* vert2prim_a = vert2prim;
      for (a = 0; a < I->NVertex; a++) {
        prm = prim + *(vert2prim_a++);

        switch (prm->type) {
          case cPrimTriangle:
          case cPrimCharacter:
            if (a == prm->vert) { /* only do this calculation for one of the
                                     three vertices */
              l1 = (float)length3f(I->Precomp + I->Vert2Normal[a] * 3);
              l2 = (float)length3f(I->Precomp + I->Vert2Normal[a] * 3 + 3);
              if ((l1 >= sep) || (l2 >= sep)) {
                b = (int)ceil(l1 / sep) + 1;
                c = (int)ceil(l2 / sep) + 1;
                extra_vert += 4 * b * c;
              }
            }
            break;

          case cPrimCone:
          case cPrimCylinder:
          case cPrimSausage:
            if ((prm->l1 + 2 * prm->r1) >= sep) {
              q = ((int)(2 * (floor(prm->r1 / sep) + 1))) + 1;
              q = q * q * ((int)ceil((prm->l1 + 2 * prm->r1) / sep) + 2);
              extra_vert += q;
            }
            break;
          case cPrimEllipsoid:
          case cPrimSphere:
            if (prm->r1 >= sep) {
              b = (int)(2 * floor(prm->r1 / sep) + 1);
              extra_vert += (b * b * b);
            }
            break;
        }
      } /* for */
    }

    extra_vert += I->NVertex;
    if (ok) tempVertex = pymol::malloc<float>(extra_vert * 3);
    CHECKOK(ok, tempVertex);
    if (ok) tempRef = pymol::malloc<int>(extra_vert);
    CHECKOK(ok, tempRef);

    ErrChkPtr(I->G, tempVertex); /* can happen if extra vert is unreasonable */
    ErrChkPtr(I->G, tempRef);

    /* lower indexes->flags, top is ref->lower index */
    ok &= !I->G->Interrupt;

    if (ok) {
      float *vv, *d;
      int* vert2prim_a = vert2prim;

      n = I->NVertex;

      v = tempVertex;
      vv = I->Vertex;

      memcpy(v, vv, n * sizeof(float) * 3);
      vv += n * 3;
      v += n * 3;

      for (a = 0; ok && a < I->NVertex; a++) {
        prm = prim + *(vert2prim_a++);

        switch (prm->type) {
          case cPrimTriangle:
          case cPrimCharacter:
            if (a == prm->vert) {
              /* only do this calculation for one of the three vertices */
              d1 = I->Precomp + I->Vert2Normal[a] * 3;
              d2 = I->Precomp + I->Vert2Normal[a] * 3 + 3;
              vv = I->Vertex + a * 3;
              l1 = (float)length3f(d1);
              l2 = (float)length3f(d2);
              if ((l1 >= sep) || (l2 >= sep)) {
                b = (int)floor(l1 / sep) + 1;
                c = (int)floor(l2 / sep) + 1;
                extra_vert += b * c;
                bh = (float)(b / 2) + 1;
                ch = (float)(c / 2) + 1;

                for (x = 0; x < bh; x++) {
                  const float xb = (float)x / b;

                  for (y = 0; y < ch; y++) {
                    const float yc = (float)y / c;

                    *(v++) = vv[0] + (d1[0] * xb) + (d2[0] * yc);
                    *(v++) = vv[1] + (d1[1] * xb) + (d2[1] * yc);
                    *(v++) = vv[2] + (d1[2] * xb) + (d2[2] * yc);

                    tempRef[n++] = a;
                  }
                }

                for (x = 0; x < bh; x++) {
                  const float xb = (float)x / b;

                  for (y = 0; y < ch; y++) {
                    const float yc = (float)y / c;

                    if ((xb + yc) < _p5) {
                      *(v++) = vv[0] + d1[0] * (_p5 + xb) + (d2[0] * yc);
                      *(v++) = vv[1] + d1[1] * (_p5 + xb) + (d2[1] * yc);
                      *(v++) = vv[2] + d1[2] * (_p5 + xb) + (d2[2] * yc);

                      tempRef[n++] = a;

                      *(v++) = vv[0] + (d1[0] * xb) + d2[0] * (_p5 + yc);
                      *(v++) = vv[1] + (d1[1] * xb) + d2[1] * (_p5 + yc);
                      *(v++) = vv[2] + (d1[2] * xb) + d2[2] * (_p5 + yc);

                      tempRef[n++] = a;
                    }
                  }
                }
              }
            } /* if */
            break;
          case cPrimCone:
          case cPrimCylinder:
          case cPrimSausage:

            ll = prm->l1 + 2 * prm->r1;

            if (ll >= sep) {
              d = I->Normal + 3 * I->Vert2Normal[a];
              vv = I->Vertex + a * 3;
              get_system1f3f(d, cx,
                             cy); /* creates an orthogonal system about d */

              p[0] = vv[0] - d[0] * prm->r1;
              p[1] = vv[1] - d[1] * prm->r1;
              p[2] = vv[2] - d[2] * prm->r1;
              dd[0] = d[0] * sep;
              dd[1] = d[1] * sep;
              dd[2] = d[2] * sep;

              q = (int)floor(prm->r1 / sep) + 1;

              while (1) {
                vd[0] = (p[0] += dd[0]);
                vd[1] = (p[1] += dd[1]);
                vd[2] = (p[2] += dd[2]);

                for (x = -q; x <= q; x++) {
                  for (y = -q; y <= q; y++) {
                    xs = x * sep;
                    ys = y * sep;
                    *(v++) = vd[0] + xs * cx[0] + ys * cy[0];
                    *(v++) = vd[1] + xs * cx[1] + ys * cy[1];
                    *(v++) = vd[2] + xs * cx[2] + ys * cy[2];

                    tempRef[n++] = a;
                  }
                }

                if (ll <= _0) break;
                ll -= sep;
              }
            }
            break;
          case cPrimEllipsoid:
          case cPrimSphere:
            if (prm->r1 >= sep) {
              q = (int)floor(prm->r1 / sep);
              vv = I->Vertex + a * 3;

              for (x = -q; x <= q; x++) {
                for (y = -q; y <= q; y++) {
                  for (z = -q; z <= q; z++) {
                    *(v++) = vv[0] + x * sep;
                    *(v++) = vv[1] + y * sep;
                    *(v++) = vv[2] + z * sep;

                    tempRef[n++] = a;
                  }
                }
              }
            }
            break;
        } /* end of switch */
        ok &= !I->G->Interrupt;
      }
    }
    if (n > extra_vert) {
      printf("BasisMakeMap: %d>%d\n", n, extra_vert);
      ErrFatal(I->G, "BasisMakeMap",
               "used too many extra vertices (this indicates a bug)...\n");
    }
    PRINTFB(I->G, FB_Ray, FB_Blather)
    " BasisMakeMap: %d total vertices\n", n ENDFB(I->G);

    if (ok) {
      if (volume) {
        v = tempVertex;

        min[0] = max[0] = v[0];
        min[1] = max[1] = v[1];
        min[2] = max[2] = v[2];

        v += 3;

        if (Feedback(I->G, FB_Ray, FB_Debugging)) {
          dump3f(min, " BasisMakeMap: remapped min");
          dump3f(max, " BasisMakeMap: remapped max");
          fflush(stdout);
        }

        for (a = 1; a < n; a++) {
          if (min[0] > v[0]) min[0] = v[0];
          if (max[0] < v[0]) max[0] = v[0];

          if (min[1] > v[1]) min[1] = v[1];
          if (max[1] < v[1]) max[1] = v[1];

          if (min[2] > v[2]) min[2] = v[2];
          if (max[2] < v[2]) max[2] = v[2];
          v += 3;
        }

        if (Feedback(I->G, FB_Ray, FB_Debugging)) {
          dump3f(min, " BasisMakeMap: remapped min");
          dump3f(max, " BasisMakeMap: remapped max");
          fflush(stdout);
        }
        if (min[0] < volume[0]) min[0] = volume[0];
        if (max[0] > volume[1]) max[0] = volume[1];
        if (min[1] < volume[2]) min[1] = volume[2];
        if (max[1] > volume[3]) max[1] = volume[3];

        if (min[2] < (-volume[5])) min[2] = (-volume[5]);
        if (max[2] > (-volume[4])) max[2] = (-volume[4]);

        extent[0] = min[0];
        extent[1] = max[0];
        extent[2] = min[1];
        extent[3] = max[1];
        extent[4] = min[2];
        extent[5] = max[2];
        PRINTFB(I->G, FB_Ray, FB_Blather)
        " BasisMakeMap: Extent [%8.2f %8.2f] [%8.2f %8.2f] [%8.2f %8.2f]\n",
            extent[0], extent[1], extent[2], extent[3], extent[4],
            extent[5] ENDFB(I->G);
        I->Map = new MapType(I->G, -sep, tempVertex, n, extent);
        CHECKOK(ok, I->Map);
      } else {
        I->Map = new MapType(I->G, sep, tempVertex, n, nullptr);
        CHECKOK(ok, I->Map);
      }
    }
    if (ok) n_voxel = I->Map->Dim[0] * I->Map->Dim[1] * I->Map->Dim[2];

    if (!ok) {
    } else if (perspective) {
      /* this is a new optimization which prevents primitives
         contained entirely within a single voxel from spilling over
         into neighboring voxels for performance of intersection
         checks (NOTE: this currently only helps with triangles, and
         is only useful in optimizing between Z layers).  The main
         impact of this optimization is to reduce the size of the
         express list (I->Map->Elist) array by about 3-fold */

      int* prm_spanner = nullptr;
      int* spanner = nullptr;
      float* v = tempVertex;
      int j, k, l, jj, kk, ll;
      int nVertex = I->NVertex;
      int prm_index;
      MapType* map = I->Map;

      prm_spanner = pymol::calloc<int>(n_prim);
      CHECKOK(ok, prm_spanner);
      if (ok) spanner = pymol::calloc<int>(n);
      CHECKOK(ok, spanner);

      /* figure out which primitives span more than one voxel */
      for (a = 0; ok && a < n; a++) {
        if (a < nVertex)
          prm_index = vert2prim[a];
        else
          prm_index = vert2prim[tempRef[a]];
        if (!prm_spanner[prm_index]) {
          prm = prim + prm_index;
          {
            float* vv = prm->vert * 3 + I->Vertex;
            switch (prm->type) {
              case cPrimTriangle:
              case cPrimCharacter:
                MapLocus(map, v, &j, &k, &l);
                MapLocus(map, vv, &jj, &kk, &ll);
                if ((j != jj) || (k != kk) || (l != ll)) {
                  prm_spanner[prm_index] = 1;
                }
                break;
              default: /* currently we aren't optimizing other primitives */
                prm_spanner[prm_index] = 1;
                break;
            }
          }
        }
        v += 3;
        ok &= !I->G->Interrupt;
      }
      /* now flag associated vertices */
      for (a = 0; ok && a < n; a++) {
        if (a < nVertex)
          prm_index = vert2prim[a];
        else
          prm_index = vert2prim[tempRef[a]];
        spanner[a] = prm_spanner[prm_index];
        ok &= !I->G->Interrupt;
      }
      /* and do the optimized expansion */
      ok &= MapSetupExpressPerp(I->Map, tempVertex, front, n, true, spanner);
      FreeP(spanner);
      FreeP(prm_spanner);
    } else if (n_voxel < (3 * n) && g_basis_lite_request) {
      /* lite: the GPU tracer builds its lists from the grid itself */
      FreeP(I->LiteRef);
      I->Lite = 1;
      I->LiteRef = tempRef;
      I->LiteN = n;
      tempRef = nullptr;
    } else if (n_voxel < (3 * n)) {
      ok &= MapSetupExpressXY(I->Map, n, true);
    } else {
      ok &= MapSetupExpressXYVert(I->Map, tempVertex, n, true);
    }

    if (ok && !I->Lite) {
      MapType* map = I->Map;
      int *sp, *ip, *ip0, ii;
      int* elist = map->EList.data();
      int* ehead = map->EHead.data();
      int newelem = 0, neelem = -map->size();
      int i_nVertex = I->NVertex;
      const int iMin0 = map->iMin[0];
      const int iMin1 = map->iMin[1];
      const int iMin2 = map->iMin[2];
      const int iMax0 = map->iMax[0];
      const int iMax1 = map->iMax[1];
      const int iMax2 = map->iMax[2];

      /* now do a filter-reassignment pass to remap fake vertices
         to the original line vertex while deleting duplicate entries */

      memset(tempRef, 0, sizeof(int) * i_nVertex);

      if (n_voxel < (3 * n)) { /* faster to traverse the entire map */
        BasisRemapExpressXY(map, tempRef, i_nVertex);
      } else {
        int* start;
        int jm1, jp1, km1, kp1, lm1, lp1;
        MapType* map = I->Map;

        v = tempVertex;

        for (e = 0; e < n; e++) {
          MapLocus(map, v, &j, &k, &l);
          jm1 = j - 1;
          jp1 = j + 1;
          km1 = k - 1;
          kp1 = k + 1;
          lm1 = l - 1;
          lp1 = l + 1;

          if (perspective) { /* for normal maps */

            int* iPtr1 = map->EHead.data() + ((j - 1) * map->D1D2) +
                         ((k - 1) * map->Dim[2]) + (l - 1);
            for (a = jm1; a <= jp1; a++) {
              int* iPtr2 = iPtr1;
              if ((a >= iMin0) && (a <= iMax0)) {
                for (b = km1; b <= kp1; b++) {
                  int* iPtr3 = iPtr2;
                  if ((b >= iMin1) && (b <= iMax1)) {
                    for (c = lm1; c <= lp1; c++) {
                      if ((c >= iMin2) && (c <= iMax2)) {
                        start = iPtr3;
                        /*start   = MapEStart(map,a,b,c); */
                        h = *start;
                        if ((h < 0) && (h >= neelem)) {
                          sp = elist - h;
                          (*start) = -h; /* no repeat visits */
                          i = *(sp++);
                          ip = ip0 = (sp - 1);

                          while (i >= 0) {
                            int ii = *(sp++);
                            if (i >= i_nVertex) i = tempRef[i];
                            if (!tempRef[i]) { /*eliminate duplicates */
                              tempRef[i] = 1;
                              *(ip++) = i;
                            }
                            i = ii;
                          }

                          *(ip++) = -1; /* terminate list */
                          newelem += (ip - ip0);

                          /* now reset flags efficiently */
                          i = *(ip0++);

                          while (i >= 0) { /* unroll */
                            ii = *(ip0++);
                            tempRef[i] = 0;
                            if ((i = ii) >= 0) {
                              ii = *(ip0++);
                              tempRef[i] = 0;
                              if ((i = ii) >= 0) {
                                ii = *(ip0++);
                                tempRef[i] = 0;
                                i = ii;
                              }
                            }
                          }
                        } /* h > 0 */
                      }
                      iPtr3++;
                    }
                  }
                  iPtr2 += map->Dim[2];
                }
              }
              iPtr1 += map->D1D2;
            }
          } else { /* for XY maps... */

            c = l;
            {
              int* iPtr1 =
                  ehead + ((j - 1) * map->D1D2) + ((k - 1) * map->Dim[2]) + c;
              if ((c >= iMin2) && (c <= iMax2)) {
                for (a = jm1; a <= jp1; a++) {
                  int* iPtr2 = iPtr1;
                  if ((a >= iMin0) && (a <= iMax0)) {
                    for (b = km1; b <= kp1; b++) {
                      if ((b >= iMin1) && (b <= iMax1)) {
                        start = iPtr2;
                        /*start   = MapEStart(map,a,b,c); */
                        h = *start;
                        if ((h < 0) && (h >= neelem)) {
                          sp = elist - h;
                          (*start) = -h; /* no repeat visits */
                          i = *(sp++);
                          ip = ip0 = sp - 1;

                          while (i >= 0) {
                            int ii = *(sp++);

                            if (i >= i_nVertex) i = tempRef[i];

                            if (!tempRef[i]) { /*eliminate duplicates */
                              tempRef[i] = 1;
                              *(ip++) = i;
                            }
                            i = ii;
                          }

                          *(ip++) = -1; /* terminate list */
                          newelem += (ip - ip0);

                          /* now reset flags efficiently */
                          i = *(ip0++);
                          while (i >= 0) { /* unroll */
                            ii = *(ip0++);
                            tempRef[i] = 0;
                            if ((i = ii) >= 0) {
                              ii = *(ip0++);
                              tempRef[i] = 0;
                              if ((i = ii) >= 0) {
                                ii = *(ip0++);
                                tempRef[i] = 0;
                                i = ii;
                              }
                            }
                          }
                        } /* h > 0 */
                      }
                      iPtr2 += map->Dim[2];
                    } /* for b */
                  }
                  iPtr1 += map->D1D2;
                } /* for a */
              }
            }
          }
          v += 3; /* happens for EVERY e! */
        } /* for e */
      }
    }

    FreeP(tempVertex);
    FreeP(tempRef);

  } else {
    /* simple sphere mode */
    I->Map = new MapType(I->G, -sep, I->Vertex, I->NVertex, nullptr);
    CHECKOK(ok, I->Map);
    if (ok) {
      if (perspective) {
        ok &= MapSetupExpressPerp(I->Map, I->Vertex, front, I->NVertex, false,
                                  nullptr);
      } else {
        ok &= MapSetupExpressXYVert(I->Map, I->Vertex, I->NVertex, false);
      }
    }
  }
  return ok;
}
int BasisInit(PyMOLGlobals* G, CBasis* I) {
  int ok = true;
  I->G = G;
  I->Radius = nullptr;
  I->Radius2 = nullptr;
  I->Normal = nullptr;
  I->Vert2Normal = nullptr;
  I->Precomp = nullptr;
  I->Vertex = VLAlloc(float, 1);
  CHECKOK(ok, I->Vertex);
  if (ok) I->Radius = VLAlloc(float, 1);
  CHECKOK(ok, I->Radius);
  if (ok) I->Radius2 = VLAlloc(float, 1);
  CHECKOK(ok, I->Radius2);
  if (ok) I->Normal = VLAlloc(float, 1);
  CHECKOK(ok, I->Normal);
  if (ok) I->Vert2Normal = VLAlloc(int, 1);
  CHECKOK(ok, I->Vert2Normal);
  if (ok) I->Precomp = VLAlloc(float, 1);
  CHECKOK(ok, I->Precomp);
  I->Map = nullptr;
  I->NVertex = 0;
  I->NNormal = 0;
  I->Lite = 0;
  I->LiteRef = nullptr;
  I->LiteN = 0;
  return ok;
}
void BasisFinish(CBasis* I) {
  if (I->Map) {
    MapFree(I->Map);
    I->Map = nullptr;
  }
  VLAFreeP(I->Radius2);
  VLAFreeP(I->Radius);
  VLAFreeP(I->Vertex);
  VLAFreeP(I->Vert2Normal);
  VLAFreeP(I->Normal);
  VLAFreeP(I->Precomp);
  I->Vertex = nullptr;
  FreeP(I->LiteRef);
  I->Lite = 0;
}

/* the path-1 (ExpressXY) express lists + remap pass of BasisMakeMap, for a
   lite map: remap fake vertices to the original vertex, drop duplicates */
static void BasisRemapExpressXY(MapType* map, int* tempRef, int i_nVertex) {
  int a, b, c, h, i;
  int *sp, *ip, *ip0, ii;
  int* elist = map->EList.data();
  int neelem = -map->size();
  memset(tempRef, 0, sizeof(int) * i_nVertex);
  for (a = map->iMin[0]; a <= map->iMax[0]; a++) {
    for (b = map->iMin[1]; b <= map->iMax[1]; b++) {
      for (c = map->iMin[2]; c <= map->iMax[2]; c++) {
        int* start = MapEStart(map, a, b, c);
        h = *start;
        if ((h < 0) && (h >= neelem)) {
          sp = elist - h;
          *(start) = -h;
          i = *(sp++);
          ip = ip0 = (sp - 1);
          while (i >= 0) {
            ii = *(sp++);
            if (i >= i_nVertex) i = tempRef[i];
            if (!tempRef[i]) {
              tempRef[i] = 1;
              *(ip++) = i;
            }
            i = ii;
          }
          *(ip++) = -1;
          i = *(ip0++);
          while (i >= 0) {
            ii = *(ip0++);
            tempRef[i] = 0;
            i = ii;
          }
        }
      }
    }
  }
}

int BasisLiteFinish(CBasis* I) {
  if (!I->Lite) return true;
  int ok = MapSetupExpressXY(I->Map, I->LiteN, true);
  if (ok) BasisRemapExpressXY(I->Map, I->LiteRef, I->NVertex);
  FreeP(I->LiteRef);
  I->Lite = 0;
  return ok;
}
void BasisTrianglePrecompute(float* v0, float* v1, float* v2, float* pre) {
  float det;

  subtract3f(v1, v0, pre);
  subtract3f(v2, v0, pre + 3);

  det = pre[0] * pre[4] - pre[1] * pre[3];

  if (fabs(det) < EPSILON)
    *(pre + 6) = 0.0F;
  else {
    *(pre + 6) = 1.0F;
    *(pre + 7) = 1.0F / det;
  }
}

void BasisTrianglePrecomputePerspective(float* v0, float* v1, float* v2,
                                        float* pre) {
  subtract3f(v1, v0, pre);
  subtract3f(v2, v0, pre + 3);
}

void BasisCylinderSausagePrecompute(float* dir, float* pre) {
  float ln = (float)(1.0f / sqrt1f(dir[1] * dir[1] + dir[0] * dir[0]));
  pre[0] = dir[1] * ln;
  pre[1] = -dir[0] * ln;
}
