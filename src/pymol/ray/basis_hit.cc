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

// Finding the primitive a ray hits first, in perspective, orthoscopic and
// shadow modes.

#include "pymol/base/base.h"
#include "pymol/ray/basis.h"
#include "pymol/glyph/character.h"
#include "pymol/ray/basis_internal.h"

int BasisHitPerspective(BasisCallRec* BC) {
  CBasis* BI = BC->Basis;
  MapType* map = BI->Map;
  int iMin0 = map->iMin[0];
  int iMin1 = map->iMin[1];
  int iMin2 = map->iMin[2];
  int iMax0 = map->iMax[0];
  int iMax1 = map->iMax[1];
  int iMax2 = map->iMax[2];
  int a, b, c;

  float iDiv = map->recipDiv;
  float base0, base1, base2;

  float min0 = map->Min[0] * iDiv;
  float min1 = map->Min[1] * iDiv;
  float min2 = map->Min[2] * iDiv;

  int new_ray = !BC->pass;
  RayInfo* r = BC->rr;

  auto* cache = &BC->cache;
  auto* cache_cache = cache->Cache.data();
  auto* cache_CacheLink = cache->CacheLink.data();

  CPrimitive* r_prim = nullptr;

  if (new_ray) { /* see if we can eliminate this ray right away using the mask
                  */

    base0 = (r->base[0] * iDiv) - min0;
    base1 = (r->base[1] * iDiv) - min1;

    a = (int)base0;
    b = (int)base1;
    a += MapBorder;
    b += MapBorder;
    if (a < iMin0)
      a = iMin0;
    else if (a > iMax0)
      a = iMax0;
    if (b < iMin1)
      b = iMin1;
    else if (b > iMax1)
      b = iMax1;

    if (!*(map->EMask.data() + a * map->Dim[1] + b)) return -1;
  }

  {
    int last_a = -1, last_b = -1, last_c = -1;
    int allow_break;
    int minIndex = -1;

    float step0, step1, step2;
    float back_dist = BC->back_dist;

    const float _0 = 0.0F, _1 = 1.0F;
    float r_tri1 = _0, r_tri2 = _0, r_dist,
          dist; /* zero inits to suppress compiler warnings */
    float r_sphere0 = _0, r_sphere1 = _0, r_sphere2 = _0;
    int h, *ip;
    int excl_trans_flag;
    int *elist, local_iflag = false;
    int terminal = -1;
    int* ehead = map->EHead.data();
    int d1d2 = map->D1D2;
    int d2 = map->Dim[2];
    const int* vert2prim = BC->vert2prim;
    const float excl_trans = BC->excl_trans;
    const float BasisFudge0 = BC->fudge0;
    const float BasisFudge1 = BC->fudge1;
    int v2p;
    int i, ii;
    int n_vert = BI->NVertex, n_eElem = map->size();
    int except1 = BC->except1;
    int except2 = BC->except2;
    int check_interior_flag = BC->check_interior && !BC->pass;
    float sph[3], vt[3], tri1 = _0, tri2;
    CPrimitive* BC_prim = BC->prim;
    int* BI_Vert2Normal = BI->Vert2Normal;
    float* BI_Vertex = BI->Vertex;
    float* BI_Precomp = BI->Precomp;
    float* BI_Normal = BI->Normal;
    float* BI_Radius = BI->Radius;
    float* BI_Radius2 = BI->Radius2;
    copy3f(r->base, vt);

    elist = map->EList.data();

    r_dist = FLT_MAX;

    excl_trans_flag = (excl_trans != _0);

    if (except1 >= 0) except1 = vert2prim[except1];
    if (except2 >= 0) except2 = vert2prim[except2];

    MapCacheReset(*cache);

    { /* take steps with a Z-size equil to the grid spacing */
      float div = iDiv * (-MapGetDiv(BI->Map) / r->dir[2]);
      step0 = r->dir[0] * div;
      step1 = r->dir[1] * div;
      step2 = r->dir[2] * div;
    }

    base0 = (r->skip[0] * iDiv) - min0;
    base1 = (r->skip[1] * iDiv) - min1;
    base2 = (r->skip[2] * iDiv) - min2;

    allow_break = false;
    while (1) {
      int inside_code;
      int clamped;

      a = ((int)base0);
      b = ((int)base1);
      c = ((int)base2);

      inside_code = 1;
      clamped = false;

      a += MapBorder;
      b += MapBorder;
      c += MapBorder;
#define EDGE_ALLOWANCE 1

      if (a < iMin0) {
        if (((iMin0 - a) > EDGE_ALLOWANCE) && allow_break)
          break;
        else {
          a = iMin0;
          clamped = true;
        }
      } else if (a > iMax0) {
        if (((a - iMax0) > EDGE_ALLOWANCE) && allow_break)
          break;
        else {
          a = iMax0;
          clamped = true;
        }
      }
      if (b < iMin1) {
        if (((iMin1 - b) > EDGE_ALLOWANCE) && allow_break)
          break;
        else {
          b = iMin1;
          clamped = true;
        }
      } else if (b > iMax1) {
        if (((b - iMax1) > EDGE_ALLOWANCE) && allow_break)
          break;
        else {
          b = iMax1;
          clamped = true;
        }
      }
      if (c < iMin2) {
        if ((iMin2 - c) > EDGE_ALLOWANCE)
          break;
        else {
          c = iMin2;
          clamped = true;
        }
      } else if (c > iMax2) {
        if ((c - iMax2) > EDGE_ALLOWANCE)
          inside_code = 0;
        else {
          c = iMax2;
          clamped = true;
        }
      }
      if (inside_code && (((a != last_a) || (b != last_b) || (c != last_c)))) {
        int new_min_index;
        h = *(ehead + (d1d2 * a) + (d2 * b) + c);

        new_min_index = -1;

        if (!clamped) /* don't discard a ray until it has hit the objective at
                         least once */
          allow_break = true;

        if ((terminal > 0) && (last_c != c)) {
          if (!terminal--) break;
        }
        if ((h > 0) && (h < n_eElem)) {
          int do_loop;

          ip = elist + h;
          last_a = a;
          i = *(ip++);
          last_b = b;
          do_loop = ((i >= 0) && (i < n_vert));
          last_c = c;

          while (do_loop) { /* n_vert checking is a bug workaround */
            CPrimitive* prm;
            v2p = vert2prim[i];
            ii = *(ip++);
            prm = BC_prim + v2p;
            do_loop = ((ii >= 0) && (ii < n_vert));
            /*            if((v2p != except1) && (v2p != except2) &&
             * (!cache->cached(v2p))) { */
            if ((v2p != except1) && (v2p != except2) && (!cache_cache[v2p])) {
              int prm_type = prm->type;

              /*cache->cache(v2p); */
              cache_cache[v2p] = 1;
              cache_CacheLink[v2p] = cache->CacheStart;
              cache->CacheStart = v2p;

              switch (prm_type) {
                case cPrimTriangle:
                case cPrimCharacter: {
                  float* dir = r->dir;
                  float* d10 = BI_Precomp + BI_Vert2Normal[i] * 3;
                  float* d20 = d10 + 3;
                  float* v0;
                  float det, inv_det;
                  float pvec0, pvec1, pvec2;
                  float dir0 = dir[0], dir1 = dir[1], dir2 = dir[2];
                  float d20_0 = d20[0], d20_1 = d20[1], d20_2 = d20[2];
                  float d10_0 = d10[0], d10_1 = d10[1], d10_2 = d10[2];

                  /* cross_product3f(dir, d20, pvec); */

                  pvec0 = dir1 * d20_2 - dir2 * d20_1;
                  pvec1 = dir2 * d20_0 - dir0 * d20_2;
                  pvec2 = dir0 * d20_1 - dir1 * d20_0;

                  /* det = dot_product3f(pvec, d10); */

                  det = pvec0 * d10_0 + pvec1 * d10_1 + pvec2 * d10_2;

                  v0 = BI_Vertex + prm->vert * 3;
                  if ((det >= EPSILON) || (det <= -EPSILON)) {
                    float tvec0, tvec1, tvec2;
                    float qvec0, qvec1, qvec2;

                    inv_det = _1 / det;

                    /* subtract3f(vt,v0,tvec); */

                    tvec0 = vt[0] - v0[0];
                    tvec1 = vt[1] - v0[1];
                    tvec2 = vt[2] - v0[2];

                    /* dot_product3f(tvec,pvec) * inv_det; */
                    tri1 = (tvec0 * pvec0 + tvec1 * pvec1 + tvec2 * pvec2) *
                           inv_det;

                    /* cross_product3f(tvec,d10,qvec); */

                    qvec0 = tvec1 * d10_2 - tvec2 * d10_1;
                    qvec1 = tvec2 * d10_0 - tvec0 * d10_2;

                    if ((tri1 >= BasisFudge0) && (tri1 <= BasisFudge1)) {
                      qvec2 = tvec0 * d10_1 - tvec1 * d10_0;

                      /* dot_product3f(dir, qvec) * inv_det; */
                      tri2 = (dir0 * qvec0 + dir1 * qvec1 + dir2 * qvec2) *
                             inv_det;

                      /* dot_product3f(d20, qvec) * inv_det; */
                      dist = (d20_0 * qvec0 + d20_1 * qvec1 + d20_2 * qvec2) *
                             inv_det;

                      if ((tri2 >= BasisFudge0) && (tri2 <= BasisFudge1) &&
                          ((tri1 + tri2) <= BasisFudge1)) {
                        if ((dist < r_dist) && (dist >= _0) &&
                            (dist <= back_dist) && (prm->trans != _1)) {
                          new_min_index = prm->vert;
                          r_tri1 = tri1;
                          r_tri2 = tri2;
                          r_dist = dist;
                        }
                      }
                    }
                  }
                } break;
                case cPrimSphere: {
                  if (LineClipPoint(r->base, r->dir, BI_Vertex + i * 3, &dist,
                                    BI_Radius[i], BI_Radius2[i])) {
                    if ((dist < r_dist) && (prm->trans != _1)) {
                      if ((dist >= _0) && (dist <= back_dist)) {
                        new_min_index = prm->vert;
                        r_dist = dist;
                      } else if (check_interior_flag && (dist <= back_dist)) {
                        if (diffsq3f(vt, BI_Vertex + i * 3) < BI_Radius2[i]) {
                          local_iflag = true;
                          r_prim = prm;
                          r_dist = _0;
                          new_min_index = prm->vert;
                        }
                      }
                    }
                  }
                } break;
                case cPrimEllipsoid: {
                  if (LineClipPoint(r->base, r->dir, BI_Vertex + i * 3, &dist,
                                    BI_Radius[i], BI_Radius2[i])) {
                    if ((dist < r_dist) && (prm->trans != _1)) {
                      float* n1 = BI_Normal + BI_Vert2Normal[i] * 3;
                      if (LineClipEllipsoidPoint(r->base, r->dir,
                                                 BI_Vertex + i * 3, &dist,
                                                 BI_Radius[i], BI_Radius2[i],
                                                 prm->n0, n1, n1 + 3, n1 + 6)) {
                        if (dist < r_dist) {
                          if ((dist >= _0) && (dist <= back_dist)) {
                            new_min_index = prm->vert;
                            r_dist = dist;
                          }
                        }
                      }
                    }
                  }
                } break;

                case cPrimCylinder:
                  if (LineToSphereCapped(r->base, r->dir, BI_Vertex + i * 3,
                                         BI_Normal + BI_Vert2Normal[i] * 3,
                                         BI_Radius[i], prm->l1, sph, &tri1,
                                         prm->cap1, prm->cap2)) {
                    if (LineClipPoint(r->base, r->dir, sph, &dist, BI_Radius[i],
                                      BI_Radius2[i])) {
                      if ((dist < r_dist) && (prm->trans != _1)) {
                        if ((dist >= _0) && (dist <= back_dist)) {
                          if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                          r_sphere0 = sph[0];
                          r_sphere1 = sph[1];
                          r_sphere2 = sph[2];
                          new_min_index = prm->vert;
                          r_dist = dist;
                        } else if (check_interior_flag && (dist <= back_dist)) {
                          if (FrontToInteriorSphereCapped(
                                  vt, BI_Vertex + i * 3,
                                  BI_Normal + BI_Vert2Normal[i] * 3,
                                  BI_Radius[i], BI_Radius2[i], prm->l1,
                                  prm->cap1, prm->cap2)) {
                            local_iflag = true;
                            r_prim = prm;
                            r_dist = _0;

                            new_min_index = prm->vert;
                          }
                        }
                      }
                    }
                  }
                  break;
                case cPrimCone: {
                  float sph_rad, sph_rad_sq;
                  if (ConeLineToSphereCapped(r->base, r->dir, BI_Vertex + i * 3,
                                             BI_Normal + BI_Vert2Normal[i] * 3,
                                             BI_Radius[i], prm->r2, prm->l1,
                                             sph, &tri1, &sph_rad, &sph_rad_sq,
                                             prm->cap1, prm->cap2)) {
                    if (LineClipPoint(r->base, r->dir, sph, &dist, sph_rad,
                                      sph_rad_sq)) {
                      if ((dist < r_dist) && (prm->trans != _1)) {
                        if ((dist >= _0) && (dist <= back_dist)) {
                          if (prm->l1 > kR_SMALL4)
                            r_tri1 = tri1 / prm->l1; /* color blending */
                          r_sphere0 = sph[0];
                          r_sphere1 = sph[1];
                          r_sphere2 = sph[2];
                          new_min_index = prm->vert;
                          r_dist = dist;
                        } else if (check_interior_flag && (dist <= back_dist)) {
                          if (FrontToInteriorSphereCapped(
                                  vt, BI_Vertex + i * 3,
                                  BI_Normal + BI_Vert2Normal[i] * 3,
                                  BI_Radius[i], BI_Radius2[i], prm->l1,
                                  prm->cap1, prm->cap2)) {
                            local_iflag = true;
                            r_prim = prm;
                            r_dist = _0;
                            new_min_index = prm->vert;
                          }
                        }
                      }
                    }
                  }
                } break;
                case cPrimSausage:
                  if (LineToSphere(r->base, r->dir, BI_Vertex + i * 3,
                                   BI_Normal + BI_Vert2Normal[i] * 3,
                                   BI_Radius[i], prm->l1, sph, &tri1)) {
                    if (LineClipPoint(r->base, r->dir, sph, &dist, BI_Radius[i],
                                      BI_Radius2[i])) {
                      int tmp_flag = false;
                      if ((dist < r_dist) && (prm->trans != _1)) {
                        if ((dist >= _0) && (dist <= back_dist)) {
                          tmp_flag = true;
                          if (excl_trans_flag) {
                            if ((prm->trans > _0) && (dist < excl_trans))
                              tmp_flag = false;
                          }
                          if (tmp_flag) {
                            if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                            r_sphere0 = sph[0];
                            r_sphere1 = sph[1];
                            r_sphere2 = sph[2];
                            new_min_index = prm->vert;
                            r_dist = dist;
                          }
                        } else if (check_interior_flag && (dist <= back_dist)) {
                          if (FrontToInteriorSphere(
                                  vt, BI_Vertex + i * 3,
                                  BI_Normal + BI_Vert2Normal[i] * 3,
                                  BI_Radius[i], BI_Radius2[i], prm->l1)) {
                            local_iflag = true;
                            r_prim = prm;
                            r_dist = _0;
                            new_min_index = prm->vert;
                          }
                        }
                      }
                    }
                  }
                  break;
              } /* end of switch */
            }
            /* end of if */
            i = ii;

          } /* end of while */

          if (local_iflag) {
            r->prim = r_prim;
            r->dist = r_dist;

            break;
          }

          if (new_min_index > -1) {
            minIndex = new_min_index;

            r_prim = BC_prim + vert2prim[minIndex];

            if ((r_prim->type == cPrimSphere) ||
                (r_prim->type == cPrimEllipsoid)) {
              const float* vv = BI->Vertex + minIndex * 3;
              r_sphere0 = vv[0];
              r_sphere1 = vv[1];
              r_sphere2 = vv[2];
            }

            BC->interior_flag = local_iflag;
            r->tri1 = r_tri1;
            r->tri2 = r_tri2;
            r->prim = r_prim;
            r->dist = r_dist;
            r->sphere[0] = r_sphere0;
            r->sphere[1] = r_sphere1;
            r->sphere[2] = r_sphere2;
          }
        } /* if -- h valid */
      }
      /* end of if */
      if (minIndex > -1) {
        if (terminal < 0) terminal = EDGE_ALLOWANCE + 1;
      }

      base0 += step0;
      base1 += step1;
      base2 += step2;
      /* advance through the map one block at a time -- note that this is a
       * crappy way to walk through the map... */
    }

    BC->interior_flag = local_iflag;
    return (minIndex);
  }
}

int BasisHitOrthoscopic(BasisCallRec* BC) {
  const float _0 = 0.0F, _1 = 1.0F;
  float oppSq, dist = _0, sph[3], vt[3], tri1, tri2;
  int a, b, c, h, *ip;
  int excl_trans_flag;
  int check_interior_flag;
  int *elist, local_iflag = false;
  float minusZ[3] = {0.0F, 0.0F, -1.0F};

  CBasis* BI = BC->Basis;
  RayInfo* r = BC->rr;

  if (MapInsideXY(BI->Map, r->base, &a, &b, &c)) {
    int minIndex = -1;
    int v2p;
    int i, ii;
    int* xxtmp;
    int do_loop;
    int except1 = BC->except1;
    int except2 = BC->except2;
    int n_vert = BI->NVertex, n_eElem = BI->Map->size();
    const int* vert2prim = BC->vert2prim;
    const float front = BC->front;
    const float back = BC->back;
    const float excl_trans = BC->excl_trans;
    const float BasisFudge0 = BC->fudge0;
    const float BasisFudge1 = BC->fudge1;

    auto* cache = &BC->cache;

    float r_tri1 = _0, r_tri2 = _0,
          r_dist = _0; /* zero inits to suppress compiler warnings */
    float r_sphere0 = _0, r_sphere1 = _0, r_sphere2 = _0;
    CPrimitive* r_prim = nullptr;

    check_interior_flag = BC->check_interior && (!BC->pass);

    /* assumption: always heading in the negative Z direction with our vector...
     */
    vt[0] = r->base[0];
    vt[1] = r->base[1];
    vt[2] = r->base[2] - front;

    if (except1 >= 0) except1 = vert2prim[except1];
    if (except2 >= 0) except2 = vert2prim[except2];

    excl_trans_flag = (excl_trans != _0);

    r_dist = FLT_MAX;

    xxtmp =
        BI->Map->EHead.data() + (a * BI->Map->D1D2) + (b * BI->Map->Dim[2]) + c;

    MapCacheReset(*cache);

    elist = BI->Map->EList.data();

    while (c >= MapBorder) {
      h = *xxtmp;
      if ((h > 0) && (h < n_eElem)) {
        ip = elist + h;
        i = *(ip++);
        do_loop = ((i >= 0) && (i < n_vert));
        while (do_loop) {
          ii = *(ip++);
          v2p = vert2prim[i];
          do_loop = ((ii >= 0) && (ii < n_vert));

          if ((v2p != except1) && (v2p != except2) && (!cache->cached(v2p))) {
            CPrimitive* prm = BC->prim + v2p;
            cache->cache(v2p);

            switch (prm->type) {
              case cPrimTriangle:
              case cPrimCharacter:
                if (!prm->cull) {
                  float* pre = BI->Precomp + BI->Vert2Normal[i] * 3;

                  if (pre[6]) {
                    float* vert0 = BI->Vertex + prm->vert * 3;

                    float tvec0 = vt[0] - vert0[0];
                    float tvec1 = vt[1] - vert0[1];

                    tri1 = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                    tri2 = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];

                    if (!((tri1 < BasisFudge0) || (tri2 < BasisFudge0) ||
                          (tri1 > BasisFudge1) ||
                          ((tri1 + tri2) > BasisFudge1))) {
                      dist = (r->base[2] - (tri1 * pre[2]) - (tri2 * pre[5]) -
                              vert0[2]);

                      if ((dist < r_dist) && (dist >= front) &&
                          (dist <= back) && (prm->trans != _1)) {
                        minIndex = prm->vert;
                        r_tri1 = tri1;
                        r_tri2 = tri2;
                        r_dist = dist;
                      }
                    }
                  }
                }
                break;

              case cPrimSphere:
                oppSq = ZLineClipPoint(r->base, BI->Vertex + i * 3, &dist,
                                       BI->Radius[i]);
                if (oppSq <= BI->Radius2[i]) {
                  dist =
                      (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i] - oppSq)));

                  if ((dist < r_dist) && (prm->trans != _1)) {
                    if ((dist >= front) && (dist <= back)) {
                      minIndex = prm->vert;
                      r_dist = dist;
                    } else if (check_interior_flag) {
                      if (diffsq3f(vt, BI->Vertex + i * 3) < BI->Radius2[i]) {
                        local_iflag = true;
                        r_prim = prm;
                        r_dist = front;
                        minIndex = prm->vert;
                      }
                    }
                  }
                }
                break;
              case cPrimEllipsoid:
                oppSq = ZLineClipPoint(r->base, BI->Vertex + i * 3, &dist,
                                       BI->Radius[i]);
                if (oppSq <= BI->Radius2[i]) {
                  dist =
                      (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i] - oppSq)));

                  if ((dist < r_dist) && (prm->trans != _1)) {
                    float* n1 = BI->Normal + BI->Vert2Normal[i] * 3;
                    if (LineClipEllipsoidPoint(r->base, minusZ,
                                               BI->Vertex + i * 3, &dist,
                                               BI->Radius[i], BI->Radius2[i],
                                               prm->n0, n1, n1 + 3, n1 + 6)) {
                      if (dist < r_dist) {
                        if ((dist >= _0) && (dist <= back)) {
                          minIndex = prm->vert;
                          r_dist = dist;
                        }
                      }
                    }
                  }
                }
                break;

              case cPrimCylinder:
                if (ZLineToSphereCapped(r->base, BI->Vertex + i * 3,
                                        BI->Normal + BI->Vert2Normal[i] * 3,
                                        BI->Radius[i], prm->l1, sph, &tri1,
                                        prm->cap1, prm->cap2,
                                        BI->Precomp + BI->Vert2Normal[i] * 3)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, BI->Radius[i]);
                  if (oppSq <= BI->Radius2[i]) {
                    dist = (float)(sqrt1f(dist) -
                                   sqrt1f((BI->Radius2[i] - oppSq)));

                    if ((dist < r_dist) && (prm->trans != _1)) {
                      if ((dist >= front) && (dist <= back)) {
                        if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                        r_sphere0 = sph[0];
                        r_sphere1 = sph[1];
                        r_sphere2 = sph[2];
                        minIndex = prm->vert;
                        r_dist = dist;
                      } else if (check_interior_flag) {
                        if (FrontToInteriorSphereCapped(
                                vt, BI->Vertex + i * 3,
                                BI->Normal + BI->Vert2Normal[i] * 3,
                                BI->Radius[i], BI->Radius2[i], prm->l1,
                                prm->cap1, prm->cap2)) {
                          local_iflag = true;
                          r_prim = prm;
                          r_dist = front;
                          minIndex = prm->vert;
                        }
                      }
                    }
                  }
                }
                break;
              case cPrimCone: {
                float sph_rad, sph_rad_sq;
                if (ConeLineToSphereCapped(r->base, minusZ, BI->Vertex + i * 3,
                                           BI->Normal + BI->Vert2Normal[i] * 3,
                                           BI->Radius[i], prm->r2, prm->l1, sph,
                                           &tri1, &sph_rad, &sph_rad_sq,
                                           prm->cap1, prm->cap2)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, sph_rad);
                  if (oppSq <= sph_rad_sq) {
                    dist = (float)(sqrt1f(dist) - sqrt1f((sph_rad_sq - oppSq)));

                    if ((dist < r_dist) && (prm->trans != _1)) {
                      if ((dist >= front) && (dist <= back)) {
                        if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                        r_sphere0 = sph[0];
                        r_sphere1 = sph[1];
                        r_sphere2 = sph[2];
                        minIndex = prm->vert;
                        r_dist = dist;
                      } else if (check_interior_flag) {
                        if (FrontToInteriorSphereCapped(
                                vt, BI->Vertex + i * 3,
                                BI->Normal + BI->Vert2Normal[i] * 3, sph_rad,
                                sph_rad_sq, prm->l1, prm->cap1, prm->cap2)) {
                          local_iflag = true;
                          r_prim = prm;
                          r_dist = front;
                          minIndex = prm->vert;
                        }
                      }
                    }
                  }
                }
              } break;
              case cPrimSausage:
                if (ZLineToSphere(r->base, BI->Vertex + i * 3,
                                  BI->Normal + BI->Vert2Normal[i] * 3,
                                  BI->Radius[i], prm->l1, sph, &tri1,
                                  BI->Precomp + BI->Vert2Normal[i] * 3)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, BI->Radius[i]);
                  if (oppSq <= BI->Radius2[i]) {
                    int tmp_flag = false;

                    dist = (float)(sqrt1f(dist) -
                                   sqrt1f((BI->Radius2[i] - oppSq)));
                    if ((dist < r_dist) && (prm->trans != _1)) {
                      if ((dist >= front) && (dist <= back)) {
                        tmp_flag = true;
                        if (excl_trans_flag) {
                          if ((prm->trans > _0) && (dist < excl_trans))
                            tmp_flag = false;
                        }
                        if (tmp_flag) {
                          if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                          r_sphere0 = sph[0];
                          r_sphere1 = sph[1];
                          r_sphere2 = sph[2];
                          minIndex = prm->vert;
                          r_dist = dist;
                        }
                      } else if (check_interior_flag) {
                        if (FrontToInteriorSphere(
                                vt, BI->Vertex + i * 3,
                                BI->Normal + BI->Vert2Normal[i] * 3,
                                BI->Radius[i], BI->Radius2[i], prm->l1)) {
                          local_iflag = true;
                          r_prim = prm;
                          r_dist = front;
                          minIndex = prm->vert;
                        }
                      }
                    }
                  }
                }
                break;
            } /* end of switch */
          }
          /* end of if */
          i = ii;

        } /* end of while */
      }

      /* and of course stop when we hit the edge of the map */

      if (local_iflag) break;

      /* we've processed all primitives associated with this voxel,
         so if an intersection has been found which occurs in front of
         the next voxel, then we can stop */

      if (minIndex > -1) {
        int aa, bb, cc;

        vt[2] = r->base[2] - r_dist;
        MapLocus(BI->Map, vt, &aa, &bb, &cc);
        if (cc > c)
          break;
        else
          vt[2] = r->base[2] - front;
      }

      c--;
      xxtmp--;

    } /* end of while */

    if (minIndex > -1) {
      r_prim = BC->prim + vert2prim[minIndex];

      if ((r_prim->type == cPrimSphere) || (r_prim->type == cPrimEllipsoid)) {
        const float* vv = BI->Vertex + minIndex * 3;
        r_sphere0 = vv[0];
        r_sphere1 = vv[1];
        r_sphere2 = vv[2];
      }
    }

    BC->interior_flag = local_iflag;
    r->tri1 = r_tri1;
    r->tri2 = r_tri2;
    r->prim = r_prim;
    r->dist = r_dist;
    r->sphere[0] = r_sphere0;
    r->sphere[1] = r_sphere1;
    r->sphere[2] = r_sphere2;
    return (minIndex);
  } /* end of if */
  BC->interior_flag = local_iflag;
  return (-1);
}

int BasisHitShadow(BasisCallRec* BC) {
  const float _0 = 0.0F;
  const float _1 = 1.0F;
  float oppSq, dist = _0, tri1, tri2;
  float sph[3], vt[3];
  int h, *ip;
  int a, b, c;
  int *elist, local_iflag = false;
  float minusZ[3] = {0.0F, 0.0F, -1.0F};
  /* local copies (eliminate these extra copies later on) */

  CBasis* BI = BC->Basis;
  RayInfo* r = BC->rr;

  if (MapInsideXY(BI->Map, r->base, &a, &b, &c)) {
    int minIndex = -1;
    int v2p;
    int i, ii;
    int* xxtmp;

    int n_vert = BI->NVertex, n_eElem = BI->Map->size();
    int except1 = BC->except1;
    int except2 = BC->except2;
    const int* vert2prim = BC->vert2prim;
    const int trans_shadows = BC->trans_shadows;
    const int nearest_shadow = BC->nearest_shadow;
    const float BasisFudge0 = BC->fudge0;
    const float BasisFudge1 = BC->fudge1;
    const int label_shadow_mode = BC->label_shadow_mode;
    auto* cache = &BC->cache;
    auto* cache_cache = cache->Cache.data();
    auto* cache_CacheLink = cache->CacheLink.data();
    CPrimitive* BC_prim = BC->prim;

    float r_tri1 = _0, r_tri2 = _0,
          r_dist; /* zero inits to suppress compiler warnings */
    float r_sphere0 = _0, r_sphere1 = _0, r_sphere2 = _0;
    float r_trans = _0;
    CPrimitive* r_prim = nullptr;

    /* assumption: always heading in the negative Z direction with our vector...
     */
    vt[0] = r->base[0];
    vt[1] = r->base[1];

    if (except1 >= 0) except1 = vert2prim[except1];
    if (except2 >= 0) except2 = vert2prim[except2];

    r_trans = _1;
    r_dist = FLT_MAX;

    xxtmp =
        BI->Map->EHead.data() + (a * BI->Map->D1D2) + (b * BI->Map->Dim[2]) + c;

    MapCacheReset(*cache);

    elist = BI->Map->EList.data();

    while (c >= MapBorder) {
      h = *xxtmp;
      if ((h > 0) && (h < n_eElem)) {
        int do_loop;
        ip = elist + h;
        i = *(ip++);
        do_loop = ((i >= 0) && (i < n_vert));
        while (do_loop) {
          ii = *(ip++);
          v2p = vert2prim[i];
          do_loop = ((ii >= 0) && (ii < n_vert));
          if ((v2p != except1) && (v2p != except2) && !cache->cached(v2p)) {
            CPrimitive* prm = BC_prim + v2p;
            int prm_type;

            /*cache->cache(v2p); */
            cache_cache[v2p] = 1;
            prm_type = prm->type;
            cache_CacheLink[v2p] = cache->CacheStart;
            cache->CacheStart = v2p;

            switch (prm_type) {
              case cPrimCharacter: /* will need special handling for character
                                      shadows */
                if (label_shadow_mode & 0x2) { /* if labels case shadows... */
                  float* pre = BI->Precomp + BI->Vert2Normal[i] * 3;

                  if (pre[6]) {
                    float* vert0 = BI->Vertex + prm->vert * 3;

                    float tvec0 = vt[0] - vert0[0];
                    float tvec1 = vt[1] - vert0[1];

                    tri1 = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                    tri2 = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];

                    if (!((tri1 < BasisFudge0) || (tri2 < BasisFudge0) ||
                          (tri1 > BasisFudge1) ||
                          ((tri1 + tri2) > BasisFudge1))) {
                      dist = (r->base[2] - (tri1 * pre[2]) - (tri2 * pre[5]) -
                              vert0[2]);

                      {
                        float fc[3];
                        float trans;

                        r->tri1 = tri1;
                        r->tri2 = tri2;
                        r->dist = dist;
                        r->prim = prm;

                        {
                          float w2;
                          w2 = _1 - (r->tri1 + r->tri2);

                          fc[0] = (prm->c2[0] * r->tri1) +
                                  (prm->c3[0] * r->tri2) + (prm->c1[0] * w2);
                          fc[1] = (prm->c2[1] * r->tri1) +
                                  (prm->c3[1] * r->tri2) + (prm->c1[1] * w2);
                          fc[2] = (prm->c2[2] * r->tri1) +
                                  (prm->c3[2] * r->tri2) + (prm->c1[2] * w2);
                        }

                        trans = CharacterInterpolate(BI->G, prm->char_id, fc);

                        if (trans == _0) { /* opaque? return immed. */
                          if (dist > -kR_SMALL4) {
                            if (nearest_shadow) {
                              if (dist < r_dist) {
                                minIndex = prm->vert;
                                r_tri1 = tri1;
                                r_tri2 = tri2;
                                r_dist = dist;
                                r_trans = (r->trans = trans);
                              }
                            } else {
                              r->prim = prm;
                              r->trans = _0;
                              r->dist = dist;
                              return (1);
                            }
                          }
                        } else if (trans_shadows) {
                          if ((dist > -kR_SMALL4) &&
                              ((r_trans > trans) ||
                               (nearest_shadow && (dist < r_dist) &&
                                (r_trans >= trans)))) {
                            minIndex = prm->vert;
                            r_tri1 = tri1;
                            r_tri2 = tri2;
                            r_dist = dist;
                            r_trans = (r->trans = trans);
                          }
                        }
                      }
                    }
                  }
                }
                break;

              case cPrimTriangle: {
                float* pre = BI->Precomp + BI->Vert2Normal[i] * 3;

                if (pre[6]) {
                  float* vert0 = BI->Vertex + prm->vert * 3;

                  float tvec0 = vt[0] - vert0[0];
                  float tvec1 = vt[1] - vert0[1];

                  tri1 = (tvec0 * pre[4] - tvec1 * pre[3]) * pre[7];
                  tri2 = -(tvec0 * pre[1] - tvec1 * pre[0]) * pre[7];
                  if (!((tri1 < BasisFudge0) || (tri2 < BasisFudge0) ||
                        (tri1 > BasisFudge1) ||
                        ((tri1 + tri2) > BasisFudge1))) {
                    float* tr = prm->tr;
                    float trans = _0;

                    dist = (r->base[2] - (tri1 * pre[2]) - (tri2 * pre[5]) -
                            vert0[2]);

                    if (prm->trans != _0) {
                      trans = (tr[1] * tri1) + (tr[2] * tri2) +
                              (tr[0] * (_1 - (tri1 + tri2)));
                    }

                    if (trans == _0) {
                      if (dist > -kR_SMALL4) {
                        if (nearest_shadow) { /* do we need the nearest shadow?
                                               */
                          if (dist < r_dist) {
                            minIndex = prm->vert;
                            r_tri1 = tri1;
                            r_tri2 = tri2;
                            r_dist = dist;
                            r_trans = (r->trans = trans);
                          }
                        } else {
                          r->prim = prm;
                          r->trans = _0;
                          r->dist = dist;
                          return (1);
                        }
                      }
                    } else if (trans_shadows) {
                      if ((dist > -kR_SMALL4) &&
                          ((r_trans > trans) ||
                           (nearest_shadow && (dist < r_dist) &&
                            (r_trans >= trans)))) {
                        minIndex = prm->vert;
                        r_tri1 = tri1;
                        r_tri2 = tri2;
                        r_dist = dist;
                        r_trans = (r->trans = trans);
                      }
                    }
                  }
                }
              } break;

              case cPrimSphere:

                oppSq = ZLineClipPoint(r->base, BI->Vertex + i * 3, &dist,
                                       BI->Radius[i]);
                if (oppSq <= BI->Radius2[i]) {
                  dist =
                      (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i] - oppSq)));

                  if (prm->trans == _0) {
                    if (dist > -kR_SMALL4) {
                      if (nearest_shadow) {
                        if (dist < r_dist) {
                          minIndex = prm->vert;
                          r_dist = dist;
                          r_trans = (r->trans = prm->trans);
                        }
                      } else {
                        r->prim = prm;
                        r->trans = prm->trans;
                        r->dist = dist;
                        return (1);
                      }
                    }
                  } else if (trans_shadows) {
                    if ((dist > -kR_SMALL4) &&
                        ((r_trans > prm->trans) ||
                         (nearest_shadow && (dist < r_dist) &&
                          (r_trans >= prm->trans)))) {
                      minIndex = prm->vert;
                      r_dist = dist;
                      r_trans = (r->trans = prm->trans);
                    }
                  }
                }
                break;

              case cPrimEllipsoid:

                oppSq = ZLineClipPointNoZCheck(r->base, BI->Vertex + i * 3,
                                               &dist, BI->Radius[i]);
                if (oppSq <= BI->Radius2[i]) {
                  dist =
                      (float)(sqrt1f(dist) - sqrt1f((BI->Radius2[i] - oppSq)));

                  if ((dist < r_dist) || (trans_shadows && (r_trans != _0))) {
                    float* n1 = BI->Normal + BI->Vert2Normal[i] * 3;
                    if (LineClipEllipsoidPoint(r->base, minusZ,
                                               BI->Vertex + i * 3, &dist,
                                               BI->Radius[i], BI->Radius2[i],
                                               prm->n0, n1, n1 + 3, n1 + 6)) {
                      if (prm->trans == _0) {
                        if (dist > -kR_SMALL4) {
                          if (nearest_shadow) {
                            if (dist < r_dist) {
                              minIndex = prm->vert;
                              r_dist = dist;
                              r_trans = (r->trans = prm->trans);
                            }
                          } else {
                            r->prim = prm;
                            r->trans = prm->trans;
                            r->dist = dist;
                            return (1);
                          }
                        }
                      } else if (trans_shadows) {
                        if ((dist > -kR_SMALL4) &&
                            ((r_trans > prm->trans) ||
                             (nearest_shadow && (dist < r_dist) &&
                              (r_trans >= prm->trans)))) {
                          minIndex = prm->vert;
                          r_dist = dist;
                          r_trans = (r->trans = prm->trans);
                        }
                      }
                    }
                  }
                }
                break;
              case cPrimCone: {
                float sph_rad, sph_rad_sq;
                if (ConeLineToSphereCapped(r->base, minusZ, BI->Vertex + i * 3,
                                           BI->Normal + BI->Vert2Normal[i] * 3,
                                           BI->Radius[i], prm->r2, prm->l1, sph,
                                           &tri1, &sph_rad, &sph_rad_sq,
                                           cCylCap::Flat, cCylCap::Flat)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, sph_rad);
                  if (oppSq <= sph_rad_sq) {
                    dist = (float)(sqrt1f(dist) - sqrt1f((sph_rad_sq - oppSq)));

                    if (prm->trans == _0) {
                      if (dist > -kR_SMALL4) {
                        if (nearest_shadow) {
                          if (dist < r_dist) {
                            if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;
                            r_sphere0 = sph[0];
                            r_sphere1 = sph[1];
                            r_sphere2 = sph[2];
                            minIndex = prm->vert;
                            r->trans = prm->trans;
                            r_dist = dist;
                            r_trans = (r->trans = prm->trans);
                          }
                        } else {
                          r->prim = prm;
                          r->trans = prm->trans;
                          r->dist = dist;
                          return (1);
                        }
                      }
                    } else if (trans_shadows) {
                      if ((dist > -kR_SMALL4) &&
                          ((r_trans > prm->trans) ||
                           (nearest_shadow && (dist < r_dist) &&
                            (r_trans >= prm->trans)))) {
                        if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;
                        r_sphere0 = sph[0];
                        r_sphere1 = sph[1];
                        r_sphere2 = sph[2];
                        minIndex = prm->vert;
                        r->trans = prm->trans;
                        r_dist = dist;
                        r_trans = (r->trans = prm->trans);
                      }
                    }
                  }
                }
              } break;
              case cPrimCylinder:
                if (ZLineToSphereCapped(r->base, BI->Vertex + i * 3,
                                        BI->Normal + BI->Vert2Normal[i] * 3,
                                        BI->Radius[i], prm->l1, sph, &tri1,
                                        prm->cap1, prm->cap2,
                                        BI->Precomp + BI->Vert2Normal[i] * 3)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, BI->Radius[i]);
                  if (oppSq <= BI->Radius2[i]) {
                    dist = (float)(sqrt1f(dist) -
                                   sqrt1f((BI->Radius2[i] - oppSq)));

                    if (prm->trans == _0) {
                      if (dist > -kR_SMALL4) {
                        if (nearest_shadow) {
                          if (dist < r_dist) {
                            if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;
                            r_sphere0 = sph[0];
                            r_sphere1 = sph[1];
                            r_sphere2 = sph[2];
                            minIndex = prm->vert;
                            r->trans = prm->trans;
                            r_dist = dist;
                            r_trans = (r->trans = prm->trans);
                          }
                        } else {
                          r->prim = prm;
                          r->trans = prm->trans;
                          r->dist = dist;
                          return (1);
                        }
                      }
                    } else if (trans_shadows) {
                      if ((dist > -kR_SMALL4) &&
                          ((r_trans > prm->trans) ||
                           (nearest_shadow && (dist < r_dist) &&
                            (r_trans >= prm->trans)))) {
                        if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;
                        r_sphere0 = sph[0];
                        r_sphere1 = sph[1];
                        r_sphere2 = sph[2];
                        minIndex = prm->vert;
                        r->trans = prm->trans;
                        r_dist = dist;
                        r_trans = (r->trans = prm->trans);
                      }
                    }
                  }
                }
                break;

              case cPrimSausage:
                if (ZLineToSphere(r->base, BI->Vertex + i * 3,
                                  BI->Normal + BI->Vert2Normal[i] * 3,
                                  BI->Radius[i], prm->l1, sph, &tri1,
                                  BI->Precomp + BI->Vert2Normal[i] * 3)) {
                  oppSq = ZLineClipPoint(r->base, sph, &dist, BI->Radius[i]);
                  if (oppSq <= BI->Radius2[i]) {
                    dist = (float)(sqrt1f(dist) -
                                   sqrt1f((BI->Radius2[i] - oppSq)));

                    if (prm->trans == _0) {
                      if (dist > -kR_SMALL4) {
                        if (nearest_shadow) {
                          if (dist < r_dist) {
                            if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;
                            r_sphere0 = sph[0];
                            r_sphere1 = sph[1];
                            r_sphere2 = sph[2];
                            minIndex = prm->vert;
                            r_dist = dist;
                            r_trans = (r->trans = prm->trans);
                          }
                        } else {
                          r->prim = prm;
                          r->trans = prm->trans;
                          r->dist = dist;
                          return (1);
                        }
                      }
                    } else if (trans_shadows) {
                      if ((dist > -kR_SMALL4) &&
                          ((r_trans > prm->trans) ||
                           (nearest_shadow && (dist < r_dist) &&
                            (r_trans >= prm->trans)))) {
                        if (prm->l1 > kR_SMALL4) r_tri1 = tri1 / prm->l1;

                        r_sphere0 = sph[0];
                        r_sphere1 = sph[1];
                        r_sphere2 = sph[2];
                        minIndex = prm->vert;
                        r_dist = dist;
                        r_trans = (r->trans = prm->trans);
                      }
                    }
                  }
                }
                break;
            } /* end of switch */
          }
          /* end of if */
          i = ii;
        } /* end of while */
      }

      /* and of course stop when we hit the edge of the map */

      if (local_iflag) break;

      /* we've processed all primitives associated with this voxel,
         so if an intersection has been found which occurs in front of
         the next voxel, then we can stop */

      /* this optimization invalid for transparent surfaces

         if( minIndex > -1 )
         {
         int   aa,bb,cc;

         vt[2]   = r->base[2] - r_dist;
         MapLocus(BI->Map,vt,&aa,&bb,&cc);
         if(cc > c)
         break;
         }
       */

      c--;
      xxtmp--;

    } /* end of while */

    if (minIndex > -1) {
      r_prim = BC->prim + vert2prim[minIndex];

      if ((r_prim->type == cPrimSphere) || (r_prim->type == cPrimEllipsoid)) {
        const float* vv = BI->Vertex + minIndex * 3;
        r_sphere0 = vv[0];
        r_sphere1 = vv[1];
        r_sphere2 = vv[2];
      }
    }

    BC->interior_flag = local_iflag;
    r->tri1 = r_tri1;
    r->tri2 = r_tri2;
    r->prim = r_prim;
    r->dist = r_dist;
    r->trans = r_trans;
    r->sphere[0] = r_sphere0;
    r->sphere[1] = r_sphere1;
    r->sphere[2] = r_sphere2;
    return (minIndex);
  } /* end of if */
  BC->interior_flag = local_iflag;
  return (-1);
}
