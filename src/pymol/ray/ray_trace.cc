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

// The CPU ray tracer: tracing, shading and antialiasing threads.

#include "pymol/base/busy.h"
#include "pymol/ray/ray.h"
#include "pymol/python/p.h"
#include "pymol/glyph/character.h"
#include "pymol/scene/scene.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/feedback.h"
#include "pymol/ray/ray_internal.h"
static void RayGetSphereNormal(CRay* I, RayInfo* r) {
  r->impact[0] = r->base[0];
  r->impact[1] = r->base[1];
  r->impact[2] = r->base[2] - r->dist;

  r->surfnormal[0] = r->impact[0] - r->sphere[0];
  r->surfnormal[1] = r->impact[1] - r->sphere[1];
  r->surfnormal[2] = r->impact[2] - r->sphere[2];

  normalize3f(r->surfnormal);
}

static void RayGetSphereNormalPerspective(CRay* I, RayInfo* r) {
  r->impact[0] = r->base[0] + r->dist * r->dir[0];
  r->impact[1] = r->base[1] + r->dist * r->dir[1];
  r->impact[2] = r->base[2] + r->dist * r->dir[2];

  r->surfnormal[0] = r->impact[0] - r->sphere[0];
  r->surfnormal[1] = r->impact[1] - r->sphere[1];
  r->surfnormal[2] = r->impact[2] - r->sphere[2];

  normalize3f(r->surfnormal);
}

void fill(unsigned int* buffer, unsigned int value, size_t cnt) {
  // std::fill_n(buffer, cnt, value);
  for (auto it_end = buffer + cnt; buffer != it_end;) {
    *(buffer++) = value;
  }
}

static void add4ucf(unsigned char* ucval, float* fval, float mulv) {
  fval[0] += ucval[0] * mulv;
  fval[1] += ucval[1] * mulv;
  fval[2] += ucval[2] * mulv;
  fval[3] += ucval[3] * mulv;
}

static void copy4fuc(float* fval, unsigned char* ucval) {
  ucval[0] = (0xFF & (int)pymol_roundf(fval[0]));
  ucval[1] = (0xFF & (int)pymol_roundf(fval[1]));
  ucval[2] = (0xFF & (int)pymol_roundf(fval[2]));
  ucval[3] = (0xFF & (int)pymol_roundf(fval[3]));
}

static float fmodpos(float a, float b) {
  float ret = fmod(a, b);
  if (ret < 0.f) {
    ret = fmod(-ret, b);
    ret = fmod(b - ret, b);
  }
  return ret;
}

static void compute_background_for_pixel(
    unsigned char* bkrd, short isOutsideInY, int bg_image_mode,
    const float* bg_image_tilesize, const float* bg_rgb, int bg_image_linear,
    unsigned char* bkgrd_data, int bkgrd_width, int bkgrd_height, float w,
    float wr, float hl, float hpixelx, short blend_bg) {
  short isBackground = isOutsideInY;
  float wl;
  switch (bg_image_mode) {
    case 1:  // isCentered
    {
      float tmpx = floor(w - hpixelx);
      isBackground |= (tmpx < 0.f || tmpx > (float)bkgrd_width);
      wl = fmodpos(tmpx, (float)bkgrd_width);
    } break;
    case 2:  // isTiled
      wl = bkgrd_width *
           (fmodpos((float)w, bg_image_tilesize[0]) / bg_image_tilesize[0]);
      break;
    case 3:  // isCenteredRepeated
      wl = fmodpos(floor(w - hpixelx), (float)bkgrd_width);
      break;
    default:
      wl = w * wr;
  }
  if (isBackground) {
    copy3f(bg_rgb, bkrd);
    bkrd[3] = 1.f;
  } else if (bg_image_linear) {
    int wlf = floor(wl), hlf = floor(hl);
    float xp = (wl - wlf) - .5f, yp = (hl - hlf) - .5f;
    float xpa = fabs(xp), ypa = fabs(yp);
    float xpam = 1.f - xpa, ypam = 1.f - ypa;
    float valx[4] = {0.f, 0.f, 0.f, 0.f}, valy[4] = {0.f, 0.f, 0.f, 0.f},
          val[4] = {0.f, 0.f, 0.f, 0.f};
    int xd = (xp > 0.f) ? 1 : -1, yd = (yp > 0.f) ? 1 : -1;
    int wlfn = (wlf + xd + bkgrd_width) % bkgrd_width,
        hlfn = (hlf + yd + bkgrd_height) % bkgrd_height;

    add4ucf(&bkgrd_data[(bkgrd_width * hlf + wlf) * 4], valx, xpam);
    add4ucf(&bkgrd_data[(bkgrd_width * hlf + wlfn) * 4], valx, xpa);

    add4ucf(&bkgrd_data[(bkgrd_width * hlfn + wlf) * 4], valy, xpam);
    add4ucf(&bkgrd_data[(bkgrd_width * hlfn + wlfn) * 4], valy, xpa);

    mult4f(valx, ypam, val);
    mult4f(valy, ypa, valy);
    add4f(valy, val, val);
    copy4fuc(val, bkrd);
  } else {
    copy4f(&bkgrd_data[(bkgrd_width * (int)floor(hl) + (int)floor(wl)) * 4],
           bkrd);
  }
  // alpha blending with bg_rgb if needed
  if (blend_bg && bkrd[3] < 255) {
    float alpha = (255.f - bkrd[3]) / 255.f;
    float tmpadd1[3], tmpadd2[3];
    tmpadd1[0] = bkrd[0];
    tmpadd1[1] = bkrd[1];
    tmpadd1[2] = bkrd[2];
    mult3f(tmpadd1, 1.f - alpha, tmpadd1);
    mult3f(bg_rgb, alpha, tmpadd2);
    add3f(tmpadd1, tmpadd2, tmpadd2);
    bkrd[0] = 0xFF & (int)pymol_roundf(tmpadd2[0]);
    bkrd[1] = 0xFF & (int)pymol_roundf(tmpadd2[1]);
    bkrd[2] = 0xFF & (int)pymol_roundf(tmpadd2[2]);
    bkrd[3] = 255;
  }
}

void fill_background_image(CRay* I, unsigned int* buffer, int width, int height,
                           unsigned int cnt) {
  int bkgrd_width = I->bkgrd_data->getWidth(),
      bkgrd_height = I->bkgrd_data->getHeight();
  unsigned char* bkgrd_data = I->bkgrd_data->bits();
  int bg_image_mode = SettingGetGlobal_i(I->G, cSetting_bg_image_mode);
  int bg_image_linear = SettingGetGlobal_b(I->G, cSetting_bg_image_linear);
  float bg_rgb[3];
  const float* tmpf;
  int w, h;
  unsigned int value;
  float wr = bkgrd_width / (float)width, hr = bkgrd_height / (float)height;
  float hl;
  unsigned int back_mask;
  unsigned char bkrd[4];
  float hpixelx = floor(width / 2.f) - floor(bkgrd_width / 2.f),
        hpixely = floor(height / 2.f) - floor(bkgrd_height / 2.f);
  auto bg_image_tilesize =
      SettingGet<const float*>(I->G, cSetting_bg_image_tilesize);
  short isOutsideInY = 0;
  int opaque_back;
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if (opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);

  tmpf =
      ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  mult3f(tmpf, 255.f, bg_rgb);

  if (opaque_back) {
    if (I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
  } else {
    back_mask = 0x00000000;
  }
  // tiled or stretched
  for (h = 0; h < height; h++) {
    switch (bg_image_mode) {
      case 1:  // isCentered
      {
        float tmpy = floor(h - hpixely);
        isOutsideInY = (tmpy < 0.f || tmpy > (float)bkgrd_height);
        hl = fmodpos(tmpy, (float)bkgrd_height);
      } break;
      case 2:  // isTiled
        hl = bkgrd_height *
             (fmodpos((float)h, bg_image_tilesize[1]) / bg_image_tilesize[1]);
        break;
      case 3:  // isCenteredRepeated
        hl = fmodpos(floor(h - hpixely), (float)bkgrd_height);
        break;
      default:
        hl = h * hr;
        break;
    }
    for (w = 0; w < width; w++) {
      compute_background_for_pixel(bkrd, isOutsideInY, bg_image_mode,
                                   bg_image_tilesize, bg_rgb, bg_image_linear,
                                   bkgrd_data, bkgrd_width, bkgrd_height, w, wr,
                                   hl, hpixelx, opaque_back);
      if (I->BigEndian) {
        value = ((0xFF & bkrd[0]) << 24) | ((0xFF & bkrd[1]) << 16) |
                ((0xFF & bkrd[2]) << 8) | (0xFF & bkrd[3]);
      } else {
        value = ((0xFF & bkrd[3]) << 24) | ((0xFF & bkrd[2]) << 16) |
                ((0xFF & bkrd[1]) << 8) | (0xFF & bkrd[0]);
      }
      if (opaque_back) {
        value = back_mask | value;
      }
      *(buffer++) = value;
    }
  }
}

void fill_gradient(CRay* I, int opaque_back, unsigned int* buffer,
                   float* bkrd_bottom, float* bkrd_top, int width, int height,
                   unsigned int cnt) {
  const float _p499 = 0.499F;
  int w, h;
  unsigned int value;
  float bkrd[3], perc;
  unsigned int back_mask;

  if (opaque_back) {
    if (I->BigEndian)
      back_mask = 0x000000FF;
    else
      back_mask = 0xFF000000;
  } else {
    back_mask = 0x00000000;
  }
  for (h = 0; h < height; h++) {
    /* for fill_gradient, y is from top to bottom */
    perc = h / (float)height;
    bkrd[0] = bkrd_top[0] + perc * (bkrd_bottom[0] - bkrd_top[0]);
    bkrd[1] = bkrd_top[1] + perc * (bkrd_bottom[1] - bkrd_top[1]);
    bkrd[2] = bkrd_top[2] + perc * (bkrd_bottom[2] - bkrd_top[2]);
    if (I->BigEndian) {
      value = back_mask |
              ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))) << 24) |
              ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 16) |
              ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 8);
    } else {
      value = back_mask |
              ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 16) |
              ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 8) |
              ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))));
    }
    for (w = 0; w < width; w++) {
      *(buffer++) = value;
    }
  }
}
static void RayReflectAndTexture(CRay* I, RayInfo* r, int perspective) {
  if (r->prim->wobble) switch (r->prim->wobble) {
      case 1:
        scatter3f(r->surfnormal, I->WobbleParam[0]);
        break;
      case 2:
        wiggle3f(r->surfnormal, r->impact, I->WobbleParam);
        break;
      case 3: {
        float3 v;
        float3 n;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = (float)cos((v[0] + v[1] + v[2]) * I->WobbleParam[1]);
        n[1] = (float)cos((v[0] - v[1] + v[2]) * I->WobbleParam[1]);
        n[2] = (float)cos((v[0] + v[1] - v[2]) * I->WobbleParam[1]);
        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, I->WobbleParam[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      }
      case 4: {
        float3 v;
        float3 n;
        float* tp = I->WobbleParam;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = I->Random[0xFF & (int)((cos((v[0]) * tp[1]) * 256 * tp[2]))];
        n[1] =
            I->Random[0xFF & (int)((cos((v[1]) * tp[1]) * 256 * tp[2] + 96))];
        n[2] =
            I->Random[0xFF & (int)((cos((v[2]) * tp[1]) * 256 * tp[2] + 148))];
        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, tp[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      } break;
      case 5: {
        float3 v;
        float3 n;
        float* tp = I->WobbleParam;
        copy3f(r->impact, v);
        RayApplyMatrixInverse33(1, &v, I->ModelView, &v);
        n[0] = I->Random[0xFF & (int)((v[0] * tp[1]) + 0)] +
               I->Random[0xFF & (int)((v[1] * tp[1]) + 20)] +
               I->Random[0xFF & (int)((v[2] * tp[1]) + 40)];
        n[1] = I->Random[0xFF & (int)((-v[0] * tp[1]) + 90)] +
               I->Random[0xFF & (int)((v[1] * tp[1]) + 100)] +
               I->Random[0xFF & (int)((-v[2] * tp[1]) + 120)];
        n[2] = I->Random[0xFF & (int)((v[0] * tp[1]) + 200)] +
               I->Random[0xFF & (int)((-v[1] * tp[1]) + 70)] +
               I->Random[0xFF & (int)((v[2] * tp[1]) + 30)];

        n[0] += I->Random[0xFF & ((int)((v[0] - v[1]) * tp[1]) + 0)] +
                I->Random[0xFF & ((int)((v[1] - v[2]) * tp[1]) + 20)] +
                I->Random[0xFF & ((int)((v[2] - v[0]) * tp[1]) + 40)];
        n[1] += I->Random[0xFF & ((int)((v[0] + v[1]) * tp[1]) + 10)] +
                I->Random[0xFF & ((int)((v[1] + v[2]) * tp[1]) + 90)] +
                I->Random[0xFF & ((int)((v[2] + v[0]) * tp[1]) + 30)];
        n[2] += I->Random[0xFF & ((int)((-v[0] + v[1]) * tp[1]) + 220)] +
                I->Random[0xFF & ((int)((-v[1] + v[2]) * tp[1]) + 20)] +
                I->Random[0xFF & ((int)((-v[2] + v[0]) * tp[1]) + 50)];

        n[0] += I->Random[0xFF & ((int)((v[0] + v[1] + v[2]) * tp[1]) + 5)] +
                I->Random[0xFF & ((int)((v[0] + v[1] + v[2]) * tp[1]) + 25)] +
                I->Random[0xFF & ((int)((v[0] + v[1] + v[2]) * tp[1]) + 46)];
        n[1] += I->Random[0xFF & ((int)((-v[0] - v[1] + v[2]) * tp[1]) + 90)] +
                I->Random[0xFF & ((int)((-v[0] - v[1] + v[2]) * tp[1]) + 45)] +
                I->Random[0xFF & ((int)((-v[0] - v[1] + v[2]) * tp[1]) + 176)];
        n[2] += I->Random[0xFF & ((int)((v[0] + v[1] - v[2]) * tp[1]) + 192)] +
                I->Random[0xFF & ((int)((v[0] + v[1] - v[2]) * tp[1]) + 223)] +
                I->Random[0xFF & ((int)((v[0] + v[1] - v[2]) * tp[1]) + 250)];

        RayTransformNormals33(1, &n, I->ModelView, &n);
        scale3f(n, tp[0], n);
        add3f(n, r->surfnormal, r->surfnormal);
        normalize3f(r->surfnormal);
      } break;
    }
  if (perspective) {
    r->dotgle = dot_product3f(r->dir, r->surfnormal);
    r->flat_dotgle = -r->dotgle;

    r->reflect[0] = r->dir[0] - (2 * r->dotgle * r->surfnormal[0]);
    r->reflect[1] = r->dir[1] - (2 * r->dotgle * r->surfnormal[1]);
    r->reflect[2] = r->dir[2] - (2 * r->dotgle * r->surfnormal[2]);
  } else {
    r->dotgle = -r->surfnormal[2];
    r->flat_dotgle = r->surfnormal[2];

    r->reflect[0] = -(2 * r->dotgle * r->surfnormal[0]);
    r->reflect[1] = -(2 * r->dotgle * r->surfnormal[1]);
    r->reflect[2] = -1.0F - (2 * r->dotgle * r->surfnormal[2]);
  }
}

void RayHashSpawn(CRayHashThreadInfo* Thread, int n_thread, int n_total) {
  int blocked;
  PyObject* info_list;
  int a, c, n = 0;
  CRay* I = Thread->ray;
  PyMOLGlobals* G = I->G;

  blocked = PAutoBlock(G);

  PRINTFB(I->G, FB_Ray, FB_Blather)
  " Ray: filling voxels with %d threads...\n", n_thread ENDFB(I->G);
  while (n < n_total) {
    c = n;
    info_list = PyList_New(n_thread);
    for (a = 0; a < n_thread; a++) {
      if ((c + a) < n_total) {
        PyList_SetItem(info_list, a,
                       PyCapsule_New(Thread + c + a, nullptr, nullptr));
      } else {
        PyList_SetItem(info_list, a, PConvAutoNone(nullptr));
      }
      n++;
    }
    PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_ray_hash_spawn", "OO",
                                 info_list, G->P_inst->cmd));
    Py_DECREF(info_list);
  }
  PAutoUnblock(G, blocked);
}

int RayHashThread(CRayHashThreadInfo* T) {
  BasisMakeMap(T->basis, T->vert2prim, T->prim, T->n_prim, T->clipBox,
               T->perspective, T->front, T->size_hint);

  /* utilize a little extra wasted CPU time in thread 0 which computes the
   * smaller map... */
  if (!T->phase) {
    if (T->ray->bkgrd_data) {
      fill_background_image(T->ray, T->image, T->width, T->height,
                            T->width * (unsigned int)T->height);
    } else if (T->bkrd_is_gradient) {
      fill_gradient(T->ray, T->opaque_back, T->image, T->bkrd_top,
                    T->bkrd_bottom, T->width, T->height,
                    T->width * (unsigned int)T->height);
    } else {
      fill(T->image, T->background, T->bytes);
    }
    RayComputeBox(T->ray);
  }
  return 1;
}

void RayTraceSpawn(CRayThreadInfo* Thread, int n_thread) {
  int blocked;
  PyObject* info_list;
  int a;
  CRay* I = Thread->ray;
  PyMOLGlobals* G = I->G;

  blocked = PAutoBlock(G);

  PRINTFB(I->G, FB_Ray, FB_Blather)
  " Ray: rendering with %d threads...\n", n_thread ENDFB(I->G);
  info_list = PyList_New(n_thread);
  for (a = 0; a < n_thread; a++) {
    PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
  }
  PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_ray_spawn", "OO", info_list,
                               G->P_inst->cmd));
  Py_DECREF(info_list);
  PAutoUnblock(G, blocked);
}

static int find_edge(
    unsigned int* ptr, float* depth, unsigned int width, int threshold,
    int back) { /* can only be called for a pixel NOT on the edge */
  {             /* color testing */
    int compare0, compare1, compare2, compare3, compare4, compare5, compare6,
        compare7, compare8;
    {
      int back_test, back_two = false;
      compare0 = (signed int)*(ptr);
      compare1 = (signed int)*(ptr - 1);
      back_test = (compare0 == back);
      compare2 = (signed int)*(ptr + 1);
      back_two = back_two || ((compare1 == back) == back_test);
      compare3 = (signed int)*(ptr - width);
      back_two = back_two || ((compare2 == back) == back_test);
      compare4 = (signed int)*(ptr + width);
      back_two = back_two || ((compare3 == back) == back_test);
      compare5 = (signed int)*(ptr - width - 1);
      back_two = back_two || ((compare4 == back) == back_test);
      compare6 = (signed int)*(ptr + width - 1);
      back_two = back_two || ((compare5 == back) == back_test);
      compare7 = (signed int)*(ptr - width + 1);
      back_two = back_two || ((compare6 == back) == back_test);
      compare8 = (signed int)*(ptr + width + 1);
      back_two = back_two || ((compare7 == back) == back_test);
      if (back_two)
        threshold =
            (threshold >> 1);  // halve threshold for pixels that hit background
    }

    {
      int current;
      unsigned int shift = 0;
      int sum1 = 0, sum2 = 3, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0,
          sum8 = 0;
      int a;
      for (a = 0; a < 4; a++) {
        current = ((compare0 >> shift) & 0xFF);
        sum1 += abs(current - ((compare1 >> shift) & 0xFF));
        sum2 += abs(current - ((compare2 >> shift) & 0xFF));
        if (sum1 >= threshold) return 1;
        sum3 += abs(current - ((compare3 >> shift) & 0xFF));
        if (sum2 >= threshold) return 1;
        sum4 += abs(current - ((compare4 >> shift) & 0xFF));
        if (sum3 >= threshold) return 1;
        sum5 += abs(current - ((compare5 >> shift) & 0xFF));
        if (sum4 >= threshold) return 1;
        sum6 += abs(current - ((compare6 >> shift) & 0xFF));
        if (sum5 >= threshold) return 1;
        sum7 += abs(current - ((compare7 >> shift) & 0xFF));
        if (sum6 >= threshold) return 1;
        sum8 += abs(current - ((compare8 >> shift) & 0xFF));
        if (sum7 >= threshold) return 1;
        if (sum8 >= threshold) return 1;
        shift += 8;
      }
    }
  }

  if (depth) { /* depth testing */
    float compare0, compare1, compare2, compare3, compare4, compare5, compare6,
        compare7, compare8;
    float dcutoff = threshold / 128.0F;

    compare1 = *(depth - 1);
    compare0 = *(depth);
    compare2 = *(depth + 1);
    if (fabs(compare0 - compare1) > dcutoff) return 1;
    compare5 = *(depth - width - 1);
    if (fabs(compare0 - compare2) > dcutoff) return 1;
    compare3 = *(depth - width);
    if (fabs(compare0 - compare5) > dcutoff) return 1;
    compare7 = *(depth - width + 1);
    if (fabs(compare0 - compare3) > dcutoff) return 1;
    compare6 = *(depth + width - 1);
    if (fabs(compare0 - compare7) > dcutoff) return 1;
    compare4 = *(depth + width);
    if (fabs(compare0 - compare6) > dcutoff) return 1;
    compare8 = *(depth + width + 1);
    if (fabs(compare0 - compare4) > dcutoff) return 1;
    if (fabs(compare0 - compare8) > dcutoff) return 1;
    /*    if(fabs(compare0-compare1)>0.001F)
       printf("%8.3f \n",compare0-compare1);
       if(fabs(compare0-compare2)>0.001F)
       printf("%8.3f \n",compare0-compare2);
       if(fabs(compare0-compare3)>0.001F)
       printf("%8.3f \n",compare0-compare3);
       if(fabs(compare0-compare4)>0.001F)
       printf("%8.3f \n",compare0-compare4);
       if(fabs(compare0-compare5)>0.001F)
       printf("%8.3f \n",compare0-compare5);
       if(fabs(compare0-compare6)>0.001F)
       printf("%8.3f \n",compare0-compare6);
       if(fabs(compare0-compare7)>0.001F)
       printf("%8.3f \n",compare0-compare7);
       if(fabs(compare0-compare8)>0.001F)
       printf("%8.3f \n",compare0-compare8); */
  }
  return 0;
}

static void RayPrimGetColorRamped(PyMOLGlobals* G, float* matrix, RayInfo* r,
                                  float* fc) {
  float fc1[3], fc2[3], fc3[3];
  float *c1, *c2, *c3, w2;
  float back_pact[3];
  const float _0 = 0.0F, _1 = 1.0F, _01 = 0.1F;
  CPrimitive* lprim = r->prim;
  inverse_transformC44f3f(matrix, r->impact, back_pact);

  switch (lprim->type) {
    case cPrimTriangle:
      w2 = 1.0F - (r->tri1 + r->tri2);
      c1 = lprim->c1;
      if (c1[0] <= _0) {
        ColorGetRamped(G, (int)(c1[0] - _01), back_pact, fc1, -1);
        c1 = fc1;
      }
      c2 = lprim->c2;
      if (c2[0] <= _0) {
        ColorGetRamped(G, (int)(c2[0] - _01), back_pact, fc2, -1);
        c2 = fc2;
      }
      c3 = lprim->c3;
      if (c3[0] <= _0) {
        ColorGetRamped(G, (int)(c3[0] - _01), back_pact, fc3, -1);
        c3 = fc3;
      }
      fc[0] = (c2[0] * r->tri1) + (c3[0] * r->tri2) + (c1[0] * w2);
      fc[1] = (c2[1] * r->tri1) + (c3[1] * r->tri2) + (c1[1] * w2);
      fc[2] = (c2[2] * r->tri1) + (c3[2] * r->tri2) + (c1[2] * w2);
      break;
    case cPrimSphere:
      c1 = lprim->c1;
      if (c1[0] <= _0) {
        ColorGetRamped(G, (int)(c1[0] - _01), back_pact, fc1, -1);
        c1 = fc1;
      }
      copy3f(c1, fc);
      break;
    case cPrimEllipsoid:
      /* TO DO */
      break;
    case cPrimCone:
    case cPrimCylinder:
    case cPrimSausage:
      w2 = r->tri1;
      c1 = lprim->c1;
      if (c1[0] <= _0) {
        ColorGetRamped(G, (int)(c1[0] - _01), back_pact, fc1, -1);
        c1 = fc1;
      }
      c2 = lprim->c2;
      if (c2[0] <= _0) {
        ColorGetRamped(G, (int)(c2[0] - _01), back_pact, fc2, -1);
        c2 = fc2;
      }
      /* TODO : here is where we would set the correct color if we implemented
         non-linear (smoothstep) half bonds. We probably need to add another
         parameter/variable to the primitive to tell this function what
         interpolation to do */
      fc[0] = (c1[0] * (_1 - w2)) + (c2[0] * w2);
      fc[1] = (c1[1] * (_1 - w2)) + (c2[1] * w2);
      fc[2] = (c1[2] * (_1 - w2)) + (c2[2] * w2);
      break;
    default:
      fc[0] = _1;
      fc[1] = _1;
      fc[2] = _1;
      break;
  }
}

int RayTraceThread(CRayThreadInfo* T) {
  CRay* I = T->ray;
  int x, y, yy;
  float excess = 0.0F;
  float dotgle;
  float bright, direct_cmp, reflect_cmp, fc[4];
  float ambient, direct, lreflect, ft, ffact = 0.0F, ffact1m;
  unsigned int cc0, cc1, cc2, cc3;
  int i;
  RayInfo r1, r2;
  int fogFlag = false;
  int fogRangeFlag = false;
  int opaque_back = 0, orig_opaque_back = 0;
  int n_hit = 0;
  int two_sided_lighting;
  float fog;
  float inter[3] = {0.0F, 0.0F, 0.0F};
  float fog_start = 0.0F;
  /*    float gamma,inp,sig=1.0F; */
  float persist, persist_inv;
  float new_front;
  int pass;
  unsigned int last_pixel = 0, *pixel;
  int exclude1, exclude2;
  float lit;
  int backface_cull;
  float project_triangle;
  float excl_trans;
  int shadows;
  int trans_shadows;
  int trans_mode;
  float first_excess;
  int pixel_flag;
  float ray_trans_spec, ray_lab_spec;
  float shadow_fudge;
  int label_shadow_mode;
  int interior_color;
  int interior_flag;
  int interior_shadows;
  int interior_wobble;
  int interior_mode;
  float interior_reflect;
  int wobble_save;
  float settingPower, settingReflectPower, settingSpecPower, settingSpecReflect,
      settingSpecDirect;
  float settingSpecDirectPower;
  float invHgt, invFrontMinusBack, inv1minusFogStart, invWdth, invHgtRange;
  float invWdthRange, vol0;
  float vol2;
  CBasis *bp1, *bp2;
  int render_height;
  int offset = 0;
  BasisCallRec BasisCall[MAX_BASIS];
  float border_offset;
  int edge_sampling = false;
  unsigned int edge_avg[4] = {0, 0, 0, 0};
  unsigned int edge_alpha_avg[4] = {0, 0, 0, 0};
  int edge_cnt = 0;
  float edge_base[2] = {0.0F, 0.0F};
  float interior_normal[3] = {0.0F, 0.0F, 0.0F};
  float edge_width = 0.35356F;
  float edge_height = 0.35356F;
  float trans_spec_cut, trans_spec_scale, trans_oblique, oblique_power;
  float direct_shade;
  float red_blend = 0.0F;
  float blue_blend = 0.0F;
  float green_blend = 0.0F;
  float trans_cont;
  float pixel_base[3];
  float inv_trans_cont = 1.0F;
  float trans_cutoff, persist_cutoff;
  int trans_cont_flag = false;
  int blend_colors;
  int max_pass;
  float BasisFudge0, BasisFudge1;
  int perspective = T->perspective;
  float eye[3] = {0.0F, 0.0F, 0.0F};
  float start[3] = {0.0F, 0.0F, 0.0F};
  float nudge[3] = {0.0F, 0.0F, 0.0F};
  float back_pact[3];
  float* depth = T->depth;
  float ray_scatter = SettingGetGlobal_f(I->G, cSetting_ray_scatter);
  const float shadow_decay =
      SettingGetGlobal_f(I->G, cSetting_ray_shadow_decay_factor);
  const float shadow_range =
      SettingGetGlobal_f(I->G, cSetting_ray_shadow_decay_range);
  const int clip_shadows = SettingGetGlobal_b(I->G, cSetting_ray_clip_shadows);
  const int spec_local = SettingGetGlobal_i(I->G, cSetting_ray_spec_local);
  float legacy = SettingGetGlobal_f(I->G, cSetting_ray_legacy_lighting);
  int spec_count = SettingGetGlobal_i(I->G, cSetting_spec_count);
  const float _0 = 0.0F;
  const float _1 = 1.0F;
  const float _p5 = 0.5F;
  const float _2 = 2.0F;
  const float _255 = 255.0F;
  const float _p499 = 0.499F;
  const float _persistLimit = 0.0001F;
  float legacy_1m = _1 - legacy;
  int n_basis = I->NBasis;
  int bg_image_mode = SettingGetGlobal_i(I->G, cSetting_bg_image_mode);
  int bg_image_linear = SettingGetGlobal_b(I->G, cSetting_bg_image_linear);
  auto bg_image_tilesize =
      SettingGet<const float*>(I->G, cSetting_bg_image_tilesize);
  float hpixelx = floor(T->width / 2.f) - floor(T->bgWidth / 2.f),
        hpixely = floor(T->height / 2.f) - floor(T->bgHeight / 2.f);
  float hl;
  float wr = T->bgWidth / (float)T->width, hr = T->bgHeight / (float)T->height;
  float bg_rgb[3];
  const float* tmpf;
  int chromadepth;

  tmpf =
      ColorGet(I->G, SettingGet_color(I->G, nullptr, nullptr, cSetting_bg_rgb));
  mult3f(tmpf, 255.f, bg_rgb);

  /*   MemoryDebugDump();
     printf("%d\n",sizeof(CPrimitive));
   */

  {
    float fudge = SettingGetGlobal_f(I->G, cSetting_ray_triangle_fudge);

    BasisFudge0 = 0.0F - fudge;
    BasisFudge1 = 1.0F + fudge;
  }
  if (spec_count < 0) {
    spec_count = SettingGetGlobal_i(I->G, cSetting_light_count);
  }
  /* SETUP */

  /*  if(T->n_thread>1)
     printf(" Ray: Thread %d: Spawned.\n",T->phase+1);
   */

  interior_shadows = SettingGetGlobal_i(I->G, cSetting_ray_interior_shadows);
  interior_wobble = SettingGetGlobal_i(I->G, cSetting_ray_interior_texture);
  interior_color = SettingGetGlobal_i(I->G, cSetting_ray_interior_color);
  interior_reflect =
      1.0F - SettingGetGlobal_f(I->G, cSetting_ray_interior_reflect);
  interior_mode = SettingGetGlobal_i(I->G, cSetting_ray_interior_mode);
  label_shadow_mode = SettingGetGlobal_i(I->G, cSetting_label_shadow_mode);
  project_triangle = SettingGetGlobal_f(I->G, cSetting_ray_improve_shadows);
  shadows = SettingGetGlobal_i(I->G, cSetting_ray_shadows);
  trans_shadows = SettingGetGlobal_i(I->G, cSetting_ray_transparency_shadows);

  backface_cull = SettingGetGlobal_i(I->G, cSetting_backface_cull);
  opaque_back = SettingGetGlobal_i(I->G, cSetting_ray_opaque_background);
  if (opaque_back < 0)
    opaque_back = SettingGetGlobal_i(I->G, cSetting_opaque_background);
  orig_opaque_back = opaque_back;

  if (T->bkrd_data) {
    opaque_back = 1;
  }

  two_sided_lighting = SettingGetGlobal_i(I->G, cSetting_two_sided_lighting);
  if (two_sided_lighting < 0) {
    if (SettingGetGlobal_i(I->G, cSetting_surface_cavity_mode))
      two_sided_lighting = true;
    else
      two_sided_lighting = false;
  }

  ray_trans_spec = SettingGetGlobal_f(I->G, cSetting_ray_transparency_specular);
  ray_lab_spec = SettingGetGlobal_f(I->G, cSetting_ray_label_specular);
  trans_cont = SettingGetGlobal_f(I->G, cSetting_ray_transparency_contrast);
  trans_mode = SettingGetGlobal_i(I->G, cSetting_transparency_mode);
  trans_oblique = SettingGetGlobal_f(I->G, cSetting_ray_transparency_oblique);
  oblique_power =
      SettingGetGlobal_f(I->G, cSetting_ray_transparency_oblique_power);
  trans_cutoff = SettingGetGlobal_f(I->G, cSetting_ray_trace_trans_cutoff);
  persist_cutoff = SettingGetGlobal_f(I->G, cSetting_ray_trace_persist_cutoff);
  chromadepth = SettingGetGlobal_i(I->G, cSetting_chromadepth);

  if (trans_mode == 1) two_sided_lighting = true;
  if (trans_cont > 1.0F) {
    trans_cont_flag = true;
    inv_trans_cont = 1.0F / trans_cont;
  }
  ambient = T->ambient;
  /* divide up the reflected light component over all lights */
  {
    float reflect_scale = SceneGetReflectScaleValue(I->G, 10);
    lreflect = reflect_scale *
               (SettingGetGlobal_f(I->G, cSetting_reflect) - ray_scatter);
    if (lreflect < _0) lreflect = _0;
    ray_scatter = ray_scatter * reflect_scale;
  }
  direct = SettingGetGlobal_f(I->G, cSetting_direct);

  /* apply legacy adjustments */

  ambient *= (1.0F - legacy) + (legacy * (0.22F / 0.12F));
  lreflect *= (1.0F - legacy) + (legacy * (0.72F / 0.45F));
  direct *= (1.0F - legacy) + (legacy * (0.24F / 0.45F));

  direct_shade = SettingGetGlobal_f(I->G, cSetting_ray_direct_shade);
  trans_spec_cut = SettingGetGlobal_f(I->G, cSetting_ray_transparency_spec_cut);
  blend_colors = SettingGetGlobal_i(I->G, cSetting_ray_blend_colors);
  max_pass = SettingGetGlobal_i(I->G, cSetting_ray_max_passes);
  if (blend_colors) {
    red_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_red);
    green_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_green);
    blue_blend = SettingGetGlobal_f(I->G, cSetting_ray_blend_blue);
  }

  if (trans_spec_cut < _1)
    trans_spec_scale = _1 / (_1 - trans_spec_cut);
  else
    trans_spec_scale = _0;

  /* COOP */
  settingPower = SettingGetGlobal_f(I->G, cSetting_power);
  settingReflectPower = SettingGetGlobal_f(I->G, cSetting_reflect_power);

  SceneGetAdjustedLightValues(I->G, &settingSpecReflect, &settingSpecPower,
                              &settingSpecDirect, &settingSpecDirectPower, 10);

  if ((interior_color != -1) || (two_sided_lighting) || (trans_mode == 1) ||
      I->CheckInterior)
    backface_cull = 0;

  shadow_fudge = SettingGetGlobal_f(I->G, cSetting_ray_shadow_fudge);

  inv1minusFogStart = _1;

  fog = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog);
  if (fog < 0.0F) {
    if (SettingGetGlobal_b(I->G, cSetting_depth_cue)) {
      fog = SettingGetGlobal_f(I->G, cSetting_fog);
    } else
      fog = _0;
  }

  if (fog != _0) {
    fogFlag = true;
    fog_start = SettingGetGlobal_f(I->G, cSetting_ray_trace_fog_start);
    if (fog_start < 0.0F)
      fog_start = SettingGetGlobal_f(I->G, cSetting_fog_start);
    if (fog_start > 1.0F) fog_start = 1.0F;
    if (abs(fog_start) > R_SMALL4) {
      fogRangeFlag = true;
      if (fabs(fog_start - 1.0F) < R_SMALL4) /* prevent div/0 */
        fogFlag = false;
    }
    inv1minusFogStart = _1 / (_1 - fog_start);
  }

  /* ray-trace */

  if (T->border) {
    invHgt = _1 / (float)(T->height - (3.0F + T->border));
    invWdth = _1 / (float)(T->width - (3.0F + T->border));
  } else {
    invHgt = _1 / (float)(T->height);
    invWdth = _1 / (float)(T->width);
  }

  if (perspective) {
    float height_range, width_range;

    /* subpixel-offsets for antialiasing naturally correspond to
       effective pixel sizes at the front of the visible slab... */

    height_range = (T->front) * 2 * ((float)tan((T->fov / 2.0F) * PI / 180.0F));
    width_range = height_range * (I->Range[0] / I->Range[1]);
    invWdthRange = invWdth * width_range;
    invHgtRange = invHgt * height_range;
    vol0 = eye[0] - width_range / 2.0F;
    vol2 = eye[1] - height_range / 2.0F;
  } else {
    invWdthRange = invWdth * I->Range[0];
    invHgtRange = invHgt * I->Range[1];
    vol0 = I->Volume[0];
    vol2 = I->Volume[2];
  }
  invFrontMinusBack = _1 / (T->front - T->back);

  edge_width *= invWdthRange;
  edge_height *= invHgtRange;

  bp1 = I->Basis + 1;
  if (I->NBasis > 2)
    bp2 = I->Basis + 2;
  else
    bp2 = nullptr;

  render_height = T->y_stop - T->y_start;

  if (render_height) {
    offset = (T->phase * render_height / T->n_thread);
    offset = offset - (offset % T->n_thread) + T->phase;
  }
  if ((interior_color != -1) || I->CheckInterior) {
    if (interior_color != -1) ColorGetEncoded(I->G, interior_color, inter);
    if (bp2) {
      interior_normal[0] = interior_reflect * bp2->LightNormal[0];
      interior_normal[1] = interior_reflect * bp2->LightNormal[1];
      interior_normal[2] = 1.0F + interior_reflect * bp2->LightNormal[2];
    } else {
      interior_normal[0] = 0.0;
      interior_normal[1] = 0.0;
      interior_normal[2] = 1.0F;
    }
    normalize3f(interior_normal);
  }

  r1.base[2] = _0;

  int* vert2prim_ptr = I->Vert2Prim.empty() ? nullptr : I->Vert2Prim.data();

  BasisCall[0].Basis = I->Basis + 1;
  BasisCall[0].rr = &r1;
  BasisCall[0].vert2prim = vert2prim_ptr;
  BasisCall[0].prim = I->Primitive;
  BasisCall[0].shadow = false;
  BasisCall[0].back = T->back;
  BasisCall[0].trans_shadows = trans_shadows;
  BasisCall[0].nearest_shadow = (shadow_decay != _0);
  BasisCall[0].check_interior = ((interior_color != -1) || I->CheckInterior);
  BasisCall[0].fudge0 = BasisFudge0;
  BasisCall[0].fudge1 = BasisFudge1;

  BasisCall[0].cache = MapCacheType(*I->Basis[1].Map);

  if (shadows && (n_basis > 2)) {
    int bc;
    for (bc = 2; bc < n_basis; bc++) {
      BasisCall[bc].Basis = I->Basis + bc;
      BasisCall[bc].rr = &r2;
      BasisCall[bc].vert2prim = vert2prim_ptr;
      BasisCall[bc].prim = I->Primitive;
      BasisCall[bc].shadow = true;
      BasisCall[bc].front = _0;
      BasisCall[bc].back = _0;
      BasisCall[bc].excl_trans = _0;
      BasisCall[bc].trans_shadows = trans_shadows;
      BasisCall[bc].nearest_shadow = (shadow_decay != _0) || (clip_shadows);
      BasisCall[bc].check_interior = false;
      BasisCall[bc].fudge0 = BasisFudge0;
      BasisCall[bc].fudge1 = BasisFudge1;
      BasisCall[bc].label_shadow_mode = label_shadow_mode;
      BasisCall[bc].cache = MapCacheType(*I->Basis[bc].Map);
    }
  }

  if (T->border) {
    border_offset = -1.50F + T->border / 2.0F;
  } else {
    border_offset = 0.0F;
  }

  unsigned int back_mask = 0x00000000;
  if (T->bkrd_is_gradient) {
    if (opaque_back) {
      if (I->BigEndian)
        back_mask = 0x000000FF;
      else
        back_mask = 0xFF000000;
    }
  }
  for (yy = T->y_start; (yy < T->y_stop); yy++) {
    float perc, bkrd[4] = {0.f, 0.f, 0.f, 1.f};
    unsigned int bkrd_value = 0;
    short isOutsideInY = 0;

    if (I->G->Interrupt) break;

    y = T->y_start +
        ((yy - T->y_start) + offset) %
            (render_height); /* make sure threads write to different pages */
    if (T->bkrd_data) {
      switch (bg_image_mode) {
        case 1:  // isCentered
        {
          float tmpy = floor(y - hpixely);
          isOutsideInY = (tmpy < 0.f || tmpy > (float)T->bgHeight);
          hl = fmodpos(tmpy, (float)T->bgHeight);
        } break;
        case 2:  // isTiled
          hl = (float)T->bgHeight *
               (fmodpos((float)y, bg_image_tilesize[1]) / bg_image_tilesize[1]);
          break;
        case 3:  // isCenteredRepeated
          hl = fmodpos(floor(y - hpixely), (float)T->bgHeight);
          break;
        default:
          hl = y * hr;
          break;
      }
    } else if (T->bkrd_is_gradient) {
      /* for RayTraceThread, y is from bottom to top */
      perc = y / (float)T->height;
      bkrd[0] = T->bkrd_bottom[0] + perc * (T->bkrd_top[0] - T->bkrd_bottom[0]);
      bkrd[1] = T->bkrd_bottom[1] + perc * (T->bkrd_top[1] - T->bkrd_bottom[1]);
      bkrd[2] = T->bkrd_bottom[2] + perc * (T->bkrd_top[2] - T->bkrd_bottom[2]);
      bkrd[3] = 1.f;
      if (T->ray->BigEndian) {
        bkrd_value = back_mask |
                     ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))) << 24) |
                     ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 16) |
                     ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 8);
      } else {
        bkrd_value = back_mask |
                     ((0xFF & ((unsigned int)(bkrd[2] * 255 + _p499))) << 16) |
                     ((0xFF & ((unsigned int)(bkrd[1] * 255 + _p499))) << 8) |
                     ((0xFF & ((unsigned int)(bkrd[0] * 255 + _p499))));
      }
    } else {
      bkrd_value = T->background;
      bkrd[0] = T->bkrd_top[0];
      bkrd[1] = T->bkrd_top[1];
      bkrd[2] = T->bkrd_top[2];
      if (orig_opaque_back) {
        bkrd[3] = 1.f;
      } else {
        bkrd[3] = 0.f;
      }
    }
    if ((!T->phase) && !(yy & 0xF)) { /* don't slow down rendering too much */
      RayProgress(0.18f + 0.78f * (T->edging ? 0.9f : 1.0f) * y /
                              std::max(T->height, 1));
      if (T->edging_cutoff) {
        if (T->edging) {
          OrthoBusyFast(I->G, (int)(2.5F * T->height / 3 + 0.5F * y),
                        4 * T->height / 3);
        } else {
          OrthoBusyFast(I->G, (int)(T->height / 3 + 0.5F * y),
                        4 * T->height / 3);
        }
      } else {
        OrthoBusyFast(I->G, T->height / 3 + y, 4 * T->height / 3);
      }
    }
    pixel = T->image + (T->width * y) + T->x_start;

    if ((y % T->n_thread) == T->phase) { /* this is my scan line */
      pixel_base[1] = ((y + 0.5F + border_offset) * invHgtRange) + vol2;

      for (x = T->x_start; (x < T->x_stop); x++) {
        if (T->bkrd_data) {
          // Need to compute background for every pixel if image-based
          unsigned char bkrd_uc[4];
          compute_background_for_pixel(
              bkrd_uc, isOutsideInY, bg_image_mode, bg_image_tilesize, bg_rgb,
              bg_image_linear, (unsigned char*)T->bkrd_data, T->bgWidth,
              T->bgHeight, x, wr, hl, hpixelx, orig_opaque_back);
          bkrd[0] = bkrd_uc[0] / 255.f;
          bkrd[1] = bkrd_uc[1] / 255.f;
          bkrd[2] = bkrd_uc[2] / 255.f;
          bkrd[3] = bkrd_uc[3] / 255.f;
          if (T->ray->BigEndian) {
            bkrd_value =
                back_mask |
                ((0xFF & ((unsigned int)(bkrd_uc[0] * 255 + _p499))) << 24) |
                ((0xFF & ((unsigned int)(bkrd_uc[1] * 255 + _p499))) << 16) |
                ((0xFF & ((unsigned int)(bkrd_uc[2] * 255 + _p499))) << 8);
          } else {
            bkrd_value =
                back_mask |
                ((0xFF & ((unsigned int)(bkrd_uc[2] * 255 + _p499))) << 16) |
                ((0xFF & ((unsigned int)(bkrd_uc[1] * 255 + _p499))) << 8) |
                ((0xFF & ((unsigned int)(bkrd_uc[0] * 255 + _p499))));
          }
        }

        pixel_base[0] = (((x + 0.5F + border_offset)) * invWdthRange) + vol0;

        while (1) {
          if (T->edging) {
            if (!edge_sampling) {
              if (x && y && (x < (T->width - 1)) &&
                  (y < (T->height - 1))) { /* not on the edge... */
                if (find_edge(T->edging + (pixel - T->image),
                              depth + (pixel - T->image), T->width,
                              T->edging_cutoff, bkrd_value)) {
                  unsigned char* pixel_c = (unsigned char*)pixel;
                  unsigned int c1, c2, c3, c4;
                  edge_cnt = 1;
                  edge_sampling = true;

                  edge_avg[0] = (c1 = pixel_c[0]);
                  edge_avg[1] = (c2 = pixel_c[1]);
                  edge_avg[2] = (c3 = pixel_c[2]);
                  edge_avg[3] = (c4 = pixel_c[3]);

                  edge_alpha_avg[0] = c1 * c4;
                  edge_alpha_avg[1] = c2 * c4;
                  edge_alpha_avg[2] = c3 * c4;
                  edge_alpha_avg[3] = c4;

                  edge_base[0] = pixel_base[0];
                  edge_base[1] = pixel_base[1];
                }
              }
            }
            if (edge_sampling) {
              if (edge_cnt == 5) {
                /* done with edging, so store averaged value */

                unsigned char* pixel_c = (unsigned char*)pixel;
                unsigned int c1, c2, c3, c4;

                edge_sampling = false;
                /* done with edging, so store averaged value */

                if (edge_alpha_avg[3]) {
                  c4 = edge_alpha_avg[3];
                  c1 = edge_alpha_avg[0] / c4;
                  c2 = edge_alpha_avg[1] / c4;
                  c3 = edge_alpha_avg[2] / c4;
                  c4 /= edge_cnt;
                } else {
                  c1 = edge_avg[0] / edge_cnt;
                  c2 = edge_avg[1] / edge_cnt;
                  c3 = edge_avg[2] / edge_cnt;
                  c4 = edge_avg[3] / edge_cnt;
                }
                pixel_c[0] = c1;
                pixel_c[1] = c2;
                pixel_c[2] = c3;
                pixel_c[3] = c4;

                /* Need to account for alpha, if exists */
                if (orig_opaque_back && c4 < 255) {
                  float bkpart = (255 - c4) / 255.f;
                  float ppart = 1.f - bkpart;
                  pixel_c[0] =
                      pymol_roundf(bkpart * bkrd[0] * 255 + ppart * pixel_c[0]);
                  pixel_c[1] =
                      pymol_roundf(bkpart * bkrd[1] * 255 + ppart * pixel_c[1]);
                  pixel_c[2] =
                      pymol_roundf(bkpart * bkrd[2] * 255 + ppart * pixel_c[2]);
                  pixel_c[3] = 255;
                }

                /* restore X,Y coordinates */
                r1.base[0] = pixel_base[0];
                r1.base[1] = pixel_base[1];

              } else {
                *pixel = 0;
                //                *pixel = bkrd_value;
                switch (edge_cnt) {
                  case 1:
                    r1.base[0] = edge_base[0] + edge_width;
                    r1.base[1] = edge_base[1] + edge_height;
                    break;
                  case 2:
                    r1.base[0] = edge_base[0] + edge_width;
                    r1.base[1] = edge_base[1] - edge_height;
                    break;
                  case 3:
                    r1.base[0] = edge_base[0] - edge_width;
                    r1.base[1] = edge_base[1] + edge_height;
                    break;
                  case 4:
                    r1.base[0] = edge_base[0] - edge_width;
                    r1.base[1] = edge_base[1] - edge_height;
                    break;
                }
              }
            }
            if (!edge_sampling) /* not oversampling this edge or already done...
                                 */
              break;
          } else {
            r1.base[0] = pixel_base[0];
            r1.base[1] = pixel_base[1];
          }

          exclude1 = -1;
          exclude2 = -1;
          persist = _1;
          first_excess = _0;
          excl_trans = _0;
          pass = 0;
          new_front = T->front;

          if (perspective) {
            r1.base[2] = -T->front;
            r1.dir[0] = (r1.base[0] - eye[0]);
            r1.dir[1] = (r1.base[1] - eye[1]);
            r1.dir[2] = (r1.base[2] - eye[2]);
            if (BasisCall[0].check_interior) {
              start[0] = r1.base[0];
              start[1] = r1.base[1];
              start[2] = r1.base[2];
            }
            normalize3f(r1.dir);
            {
              float scale = I->max_box[2] / r1.base[2];

              r1.skip[0] = r1.base[0] * scale;
              r1.skip[1] = r1.base[1] * scale;
              r1.skip[2] = I->max_box[2];
            }
          }
          while ((persist > _persistLimit) && (pass <= max_pass)) {
            char fogFlagTmp = fogFlag;
            pixel_flag = false;
            BasisCall[0].except1 = exclude1;
            BasisCall[0].except2 = exclude2;
            BasisCall[0].front = new_front;
            BasisCall[0].excl_trans = excl_trans;
            BasisCall[0].interior_flag = false;
            BasisCall[0].pass = pass;

            if (perspective) {
              if (pass) {
                add3f(nudge, r1.base, r1.base);
                copy3f(r1.base, r1.skip);
              }
              BasisCall[0].back_dist = -(T->back + r1.base[2]) / r1.dir[2];
              i = BasisHitPerspective(&BasisCall[0]);
            } else {
              i = BasisHitOrthoscopic(&BasisCall[0]);
            }
            interior_flag = BasisCall[0].interior_flag && (!pass);

            if (((i >= 0) || interior_flag) && (pass < max_pass)) {
              int n_basis_tmp = r1.prim->no_lighting ? 0 : n_basis;
              pixel_flag = true;
              n_hit++;
              if (((r1.trans = r1.prim->trans) != _0) && trans_cont_flag) {
                r1.trans = (float)pow(r1.trans, inv_trans_cont);
              }
              if (interior_flag) {
                copy3f(interior_normal, r1.surfnormal);
                if (perspective) {
                  copy3f(start, r1.impact);
                  r1.dist = _0;
                } else {
                  copy3f(r1.base, r1.impact);
                  r1.dist = T->front;
                  r1.impact[2] -= T->front;
                }

                if (interior_wobble >= 0) {
                  wobble_save =
                      r1.prim->wobble; /* This is a no-no for multithreading! */
                  r1.prim->wobble = interior_wobble;

                  RayReflectAndTexture(I, &r1, perspective);

                  r1.prim->wobble = wobble_save;
                } else
                  RayReflectAndTexture(I, &r1, perspective);

                dotgle = -r1.dotgle;
                if (interior_color >= 0 || interior_color <= cColorExtCutoff) {
                  copy3f(inter, fc);
                } else if (interior_color == cColorAtomic) {
                  // Affects spheres+cylinders, but not surfaces+cartoon.
                  copy3f(r1.prim->c1, fc);
                } else {
                  copy3f(r1.prim->ic, fc);
                }
              } else {
                if (!perspective) new_front = r1.dist;

                switch (r1.prim->type) {
                  case cPrimTriangle:

                    BasisGetTriangleNormal(bp1, &r1, i, fc, perspective);
                    r1.trans = (float)pow(r1.trans, inv_trans_cont);

                    if (r1.prim->ramped) {
                      RayPrimGetColorRamped(I->G, I->ModelView, &r1, fc);
                    }
                    if (bp2) {
                      RayProjectTriangle(
                          I, &r1, bp2->LightNormal, bp1->Vertex + i * 3,
                          bp1->Normal + bp1->Vert2Normal[i] * 3 + 3,
                          project_triangle);
                    }

                    RayReflectAndTexture(I, &r1, perspective);
                    if (perspective) {
                      BasisGetTriangleFlatDotglePerspective(bp1, &r1, i);
                    } else {
                      BasisGetTriangleFlatDotgle(bp1, &r1, i);
                    }
                    break;
                  case cPrimCharacter:
                    BasisGetTriangleNormal(bp1, &r1, i, fc, perspective);

                    r1.trans = CharacterInterpolate(I->G, r1.prim->char_id, fc);
                    fogFlagTmp = false;
                    RayReflectAndTexture(I, &r1, perspective);
                    BasisGetTriangleFlatDotgle(bp1, &r1, i);
                    break;

                  case cPrimEllipsoid:

                    BasisGetEllipsoidNormal(bp1, &r1, i, perspective);
                    RayReflectAndTexture(I, &r1, perspective);

                    fc[0] = r1.prim->c1[0];
                    fc[1] = r1.prim->c1[1];
                    fc[2] = r1.prim->c1[2];
                    break;

                  default: /* sphere, cylinder, sausage, etc. */

                    /* must be a sphere (effectively speaking) */

                    if (perspective) {
                      RayGetSphereNormalPerspective(I, &r1);
                    } else {
                      RayGetSphereNormal(I, &r1);
                    }

                    RayReflectAndTexture(I, &r1, perspective);

                    if (r1.prim->ramped) {
                      RayPrimGetColorRamped(I->G, I->ModelView, &r1, fc);
                    } else {
                      switch (r1.prim->type) {
                        case cPrimCylinder:
                        case cPrimSausage:
                        case cPrimCone:
                          ft = r1.tri1;
                          fc[0] = (r1.prim->c1[0] * (_1 - ft)) +
                                  (r1.prim->c2[0] * ft);
                          fc[1] = (r1.prim->c1[1] * (_1 - ft)) +
                                  (r1.prim->c2[1] * ft);
                          fc[2] = (r1.prim->c1[2] * (_1 - ft)) +
                                  (r1.prim->c2[2] * ft);
                          break;
                        default:
                          fc[0] = r1.prim->c1[0];
                          fc[1] = r1.prim->c1[1];
                          fc[2] = r1.prim->c1[2];
                          break;
                      }
                    }
                    break;
                }

                if ((trans_oblique != _0) && (r1.trans != _0)) {
                  if ((r1.surfnormal[2] > _0) || two_sided_lighting) {
                    float oblique_factor = r1.surfnormal[2];
                    if (oblique_factor < _0) oblique_factor = -oblique_factor;
                    if (oblique_factor != _1) {
                      if (oblique_factor > _p5) {
                        oblique_factor =
                            (float)(_p5 +
                                    _p5 * (_1 - pow((_1 - oblique_factor) * _2,
                                                    oblique_power)));
                      } else {
                        oblique_factor = (float)(_p5 * pow(oblique_factor * _2,
                                                           oblique_power));
                      }
                    }
                    r1.trans *= (trans_oblique * oblique_factor) +
                                (1.0F - trans_oblique);
                    if (r1.trans < 0.06F)
                      r1.trans =
                          0.06F; /* don't allow transparent to become opaque */
                  }
                }

                dotgle = -r1.dotgle;

                if (r1.flat_dotgle < _0) {
                  if ((!two_sided_lighting) && (BasisCall[0].check_interior) &&
                      (interior_mode != 2)) {
                    interior_flag = true;
                    copy3f(interior_normal, r1.surfnormal);
                    if (perspective) {
                      copy3f(start, r1.impact);
                      r1.dist = _0;
                    } else {
                      copy3f(r1.base, r1.impact);
                      r1.impact[2] -= T->front;
                      r1.dist = T->front;
                    }

                    if (interior_wobble >= 0) {
                      wobble_save = r1.prim->wobble;
                      r1.prim->wobble = interior_wobble;
                      RayReflectAndTexture(I, &r1, perspective);
                      r1.prim->wobble = wobble_save;
                    } else {
                      RayReflectAndTexture(I, &r1, perspective);
                    }

                    dotgle = -r1.dotgle;
                    if (interior_color >= 0 ||
                        interior_color <= cColorExtCutoff) {
                      copy3f(inter, fc);
                    } else if (interior_color == cColorAtomic) {
                      // Affects surfaces+cartoon, but not spheres+cylinders.
                      copy3f(r1.prim->c1, fc);
                    } else {
                      copy3f(r1.prim->ic, fc);
                    }
                  }
                }

                if ((r1.flat_dotgle < _0) && (!interior_flag)) {
                  if (two_sided_lighting) {
                    dotgle = -dotgle;
                    invert3f(r1.surfnormal);
                  }
                }
              }

              {
                double pow_dotgle;
                float pow_surfnormal2;

                if (settingPower != _1) {
                  pow_dotgle = pow(dotgle, settingPower);
                  pow_surfnormal2 = (float)pow(r1.surfnormal[2], settingPower);
                } else {
                  pow_dotgle = dotgle;
                  pow_surfnormal2 = r1.surfnormal[2];
                }
                direct_cmp = legacy_1m * pow_surfnormal2 + /* new model */
                             legacy * ((float)(dotgle + pow_dotgle) *
                                       _p5); /* legacy model */
              }

              reflect_cmp = _0;
              if (settingSpecDirect != _0) {
                if (r1.surfnormal[2] > _0) {
                  excess =
                      (float)(pow(r1.surfnormal[2], settingSpecDirectPower) *
                              settingSpecDirect);
                } else {
                  excess = _0;
                }
              } else {
                excess = _0;
              }

              lit = _1;
              if (n_basis_tmp < 3) {
                reflect_cmp = direct_cmp;
              } else {
                int bc;
                CBasis* bp;
                for (bc = 2; bc < n_basis; bc++) {
                  lit = _1;
                  bp = I->Basis + bc;

                  if (shadows && ((!interior_flag) || (interior_shadows)) &&
                      ((r1.prim->type != cPrimCharacter) ||
                       (label_shadow_mode & 0x1))) {
                    matrix_transform33f3f(bp->Matrix, r1.impact, r2.base);
                    r2.base[2] -= shadow_fudge;
                    BasisCall[bc].except2 = -1;
                    BasisCall[bc].except1 =
                        i; /* exclude current prim from shadow comp */
                    if (BasisHitShadow(&BasisCall[bc]) > -1) {
                      if ((!clip_shadows) || (bp->LightNormal[2] >= _0) ||
                          ((T->front + r1.impact[2] -
                            (r2.dist * bp->LightNormal[2])) < _0)) {
                        lit = (float)pow(r2.trans, _p5);
                        if ((shadow_decay != _0) && (r2.dist > shadow_range)) {
                          if (shadow_decay > 0) {
                            lit += ((_1 - lit) *
                                    (_1 - _1 / exp((r2.dist - shadow_range) *
                                                   shadow_decay)));
                          } else {
                            lit += ((_1 - lit) *
                                    (_1 - _1 / pow(r2.dist / shadow_range,
                                                   -shadow_decay)));
                          }
                        }
                      }
                    }
                  }

                  if (lit > _0) {
                    {
                      double pow_dotgle;

                      dotgle = -dot_product3f(r1.surfnormal, bp->LightNormal);
                      if (dotgle < _0) dotgle = _0;

                      if (settingReflectPower != _1)
                        pow_dotgle = pow(dotgle, settingReflectPower);
                      else
                        pow_dotgle = dotgle;

                      reflect_cmp +=
                          legacy_1m *
                              ((float)(lit * pow_dotgle)) + /* new model */
                          legacy * ((float)(lit * (dotgle + pow_dotgle) *
                                            _p5)); /* legacy model */
                    }

                    if (bc < (spec_count + 2)) {
                      if (ray_scatter != _0) /* scattered specular light */
                        excess += ray_scatter * dotgle;

                      if (spec_local && perspective) {
                        /* slower, C4D-like local specular */
                        float tmp[3];

                        add3f(r1.surfnormal, r1.surfnormal, tmp);
                        add3f(tmp, bp->LightNormal, tmp);
                        normalize3f(tmp);
                        dotgle = -(dot_product3f(r1.dir, tmp)) * 1.004F;
                        if (dotgle > _1)
                          dotgle = _1;
                        else if (dotgle < _0)
                          dotgle = _0;
                        dotgle = pow(dotgle, 0.29F);
                      } else {
                        dotgle = -dot_product3f(
                            r1.surfnormal,
                            bp->SpecNormal); /* fast OpenGL-like global specular
                                              */
                      }
                      if (dotgle < _0) dotgle = _0;
                      if (r1.prim->type != cPrimCharacter) {
                        excess += (float)(pow(dotgle, settingSpecPower) *
                                          settingSpecReflect * lit);
                      } else {
                        excess +=
                            (float)(pow(dotgle, settingSpecPower) *
                                    settingSpecReflect * lit * ray_lab_spec);
                      }
                    }
                  }
                }
                if (chromadepth > 0) {
                  float d = -(r1.impact[2] + T->front) / (T->back - T->front);
                  const float margin = 1. / 4.;
                  float h = d * (1. + 2. * margin) - margin;
                  if (h < 0.0) h = 0.0;
                  if (h > 1.0) h = 1.0;
                  float hue6 = 6. * 240. / 360. * h;
                  float rgb[3];
                  rgb[0] = fabs(hue6 - 3.) - 1.;
                  rgb[1] = 2. - fabs(hue6 - 2.);
                  rgb[2] = 2. - fabs(hue6 - 4.);
                  clamp3f(rgb);
                  if (chromadepth == 2) {
                    float lum = 0.30 * fc[0] + 0.59 * fc[1] + 0.11 * fc[2];
                    rgb[0] *= lum;
                    rgb[1] *= lum;
                    rgb[2] *= lum;
                  }
                  copy3(rgb, fc);
                }
              }

              if (fc[0] <= ((float)cColorExtCutoff)) { /* ramped color */
                inverse_transformC44f3f(I->ModelView, r1.impact, back_pact);
                ColorGetRamped(I->G, (int)(fc[0] - 0.1F), back_pact, fc, -1);
              }

              bright =
                  ambient +
                  (((_1 - direct_shade) + direct_shade * lit) * direct *
                       direct_cmp +
                   lreflect * reflect_cmp *
                       (legacy_1m + legacy * direct_cmp)); /* blend legacy */
              if (excess > _1) excess = _1;
              if (bright > _1)
                bright = _1;
              else if (bright < _0)
                bright = _0;

              /*                      bright *= (_1-excess); */

              fc[0] = (bright * fc[0] + excess);
              fc[1] = (bright * fc[1] + excess);
              fc[2] = (bright * fc[2] + excess);

              if (n_basis_tmp && fogFlagTmp) {
                if (perspective) {
                  ffact = (T->front + r1.impact[2]) * invFrontMinusBack;
                } else {
                  ffact = (T->front - r1.dist) * invFrontMinusBack;
                }
                if (fogRangeFlag)
                  ffact = (ffact - fog_start) * inv1minusFogStart;

                ffact *= fog;

                if (ffact < _0) ffact = _0;
                if (ffact > _1) ffact = _1;

                ffact1m = _1 - ffact;

                if (orig_opaque_back) {
                  fc[0] = ffact * bkrd[0] + fc[0] * ffact1m;
                  fc[1] = ffact * bkrd[1] + fc[1] * ffact1m;
                  fc[2] = ffact * bkrd[2] + fc[2] * ffact1m;
                  fc[3] = _1;
                } else {
                  fc[3] = ffact1m * (_1 - r1.trans);
                }

                if (!pass) {
                  if (r1.trans < trans_spec_cut) {
                    first_excess = excess * ffact1m * ray_trans_spec;
                  } else {
                    first_excess = excess * ffact1m * ray_trans_spec *
                                   trans_spec_scale * (_1 - r1.trans);
                  }
                } else {
                  fc[0] += first_excess; /* dubious? */
                  fc[1] += first_excess;
                  fc[2] += first_excess;
                }
              } else {
                if (!pass) {
                  if (r1.trans < trans_spec_cut) {
                    first_excess = excess * ray_trans_spec;
                  } else {
                    first_excess = excess * ray_trans_spec * trans_spec_scale *
                                   (_1 - r1.trans);
                  }
                } else {
                  fc[0] += first_excess;
                  fc[1] += first_excess;
                  fc[2] += first_excess;
                }
                if (orig_opaque_back) {
                  fc[3] = _1;
                } else {
                  fc[3] = _1 - r1.trans;
                }
              }
            } else if (pass) {
              /* hit nothing, and we're on on second or greater pass,
                 or we're on the last pass of a dead-end loop */
              i = -1;

              fc[0] = first_excess + bkrd[0];
              fc[1] = first_excess + bkrd[1];
              fc[2] = first_excess + bkrd[2];

              if (orig_opaque_back) {
                fc[3] = _1;
              } else {
                fc[3] = _0;
                //		fc[3] = bkrd[3];
              }

              ffact = 1.0F;
              ffact1m = 0.0F;

              pixel_flag = true;
              if (trans_cont_flag) persist = (float)pow(persist, trans_cont);
            }

            if (pixel_flag) {
              /*
                 inp    = (fc[0]+fc[1]+fc[2]) * _inv3;
                 if(inp < R_SMALL4)
                 sig = _1;
                 else
                 sig = (float)(pow(inp,gamma) / inp);

                 cc0 = (uint)(sig * fc[0] * _255);
                 cc1 = (uint)(sig * fc[1] * _255);
                 cc2 = (uint)(sig * fc[2] * _255);
               */

              cc0 = (uint)(fc[0] * _255);
              cc1 = (uint)(fc[1] * _255);
              cc2 = (uint)(fc[2] * _255);

              if (cc0 > 255) cc0 = 255;
              if (cc1 > 255) cc1 = 255;
              if (cc2 > 255) cc2 = 255;

              if (orig_opaque_back) {
                if (I->BigEndian)
                  *pixel =
                      T->fore_mask | (cc0 << 24) | (cc1 << 16) | (cc2 << 8);
                else
                  *pixel = T->fore_mask | (cc2 << 16) | (cc1 << 8) | cc0;
              } else {
                /* use alpha channel for fog with transparent backgrounds */
                cc3 = (uint)(fc[3] * _255);
                if (cc3 > 255) cc3 = 255;

                if (I->BigEndian)
                  *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
                else
                  *pixel = (cc3 << 24) | (cc2 << 16) | (cc1 << 8) | cc0;
              }
            }
            if (pass && persist < 1.f) { /* average all four channels */
              float mix_in;
              if (i >= 0) {
                if (fogFlagTmp) {
                  if (trans_cont_flag && (ffact > _p5)) {
                    mix_in =
                        2 *
                        (persist * (_1 - ffact) +
                         ((float)pow(persist, trans_cont) * (ffact - _p5))) *
                        (_1 - r1.trans * ffact);
                  } else {
                    mix_in = persist * (_1 - r1.trans * ffact);
                  }
                } else {
                  mix_in = persist * (_1 - r1.trans);
                }
              } else {
                mix_in = persist;
              }

              persist_inv = _1 - mix_in;

              if (!orig_opaque_back) {
                if (i < 0) { /* hit nothing -- so don't blend */
                  fc[0] = (float)(0xFF & (last_pixel >> 24));
                  fc[1] = (float)(0xFF & (last_pixel >> 16));
                  fc[2] = (float)(0xFF & (last_pixel >> 8));
                  fc[3] = (float)(0xFF & (last_pixel));
                  if (trans_cont_flag) { /* unless we are increasing contrast */
                    float m;
                    if (I->BigEndian) {
                      m = _1 - (float)(0xFF & (last_pixel)) / _255;
                    } else {
                      m = _1 - (float)(0xFF & (last_pixel >> 24)) / _255;
                    }
                    m = _1 - (float)pow(m, trans_cont);
                    if (I->BigEndian) {
                      fc[3] = m * _255 + _p499;
                    } else {
                      fc[0] = m * _255 + _p499;
                    }
                  }
                } else { /* hit something -- so keep blend and compute
                            cumulative alpha */

                  fc[0] = (0xFF & ((*pixel) >> 24)) * mix_in +
                          (0xFF & (last_pixel >> 24)) * persist_inv;
                  fc[1] = (0xFF & ((*pixel) >> 16)) * mix_in +
                          (0xFF & (last_pixel >> 16)) * persist_inv;
                  fc[2] = (0xFF & ((*pixel) >> 8)) * mix_in +
                          (0xFF & (last_pixel >> 8)) * persist_inv;
                  fc[3] = (0xFF & ((*pixel))) * mix_in +
                          (0xFF & (last_pixel)) * persist_inv;

                  if (i >= 0) { /* make sure opaque objects get opaque alpha */
                    float o1, o2;
                    float m;

                    if (I->BigEndian) {
                      o1 = (float)(0xFF & (last_pixel)) / _255;
                      o2 = (float)(0xFF & (*pixel)) / _255;
                    } else {
                      o1 = (float)(0xFF & (last_pixel >> 24)) / _255;
                      o2 = (float)(0xFF & ((*pixel) >> 24)) / _255;
                    }

                    if (o1 < o2) { /* make sure o1 is largest opacity */
                      m = o1;
                      o1 = o2;
                      o2 = m;
                    }
                    m = o1 + (1.0F - o1) * o2;
                    if (I->BigEndian) {
                      fc[3] = m * _255 + _p499;
                    } else {
                      fc[0] = m * _255 + _p499;
                    }
                  }
                }
              } else { /* opaque background, so just blend */
                fc[0] = (0xFF & ((*pixel) >> 24)) * mix_in +
                        (0xFF & (last_pixel >> 24)) * persist_inv;
                fc[1] = (0xFF & ((*pixel) >> 16)) * mix_in +
                        (0xFF & (last_pixel >> 16)) * persist_inv;
                fc[2] = (0xFF & ((*pixel) >> 8)) * mix_in +
                        (0xFF & (last_pixel >> 8)) * persist_inv;
                fc[3] = (0xFF & ((*pixel))) * mix_in +
                        (0xFF & (last_pixel)) * persist_inv;
              }

              cc0 = (uint)(fc[0]);
              cc1 = (uint)(fc[1]);
              cc2 = (uint)(fc[2]);
              cc3 = (uint)(fc[3]);

              if (cc0 > 255) cc0 = 255;
              if (cc1 > 255) cc1 = 255;
              if (cc2 > 255) cc2 = 255;
              if (cc3 > 255) cc3 = 255;

              *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
            }

            if (depth && (i >= 0) && (r1.trans < trans_cutoff) &&
                (persist > persist_cutoff)) {
              depth[pixel - T->image] = (T->front + r1.impact[2]);
            }

            if (i >= 0) {
              if (r1.prim->type ==
                  cPrimSausage) { /* carry ray through the stick */
                if (perspective)
                  excl_trans = (2 * r1.surfnormal[2] * r1.prim->r1 / r1.dir[2]);
                else
                  excl_trans = new_front + (2 * r1.surfnormal[2] * r1.prim->r1);
              }

              if ((!backface_cull) && (trans_mode != 2))
                persist = persist * r1.trans;
              else {
                if (!r1.prim->no_lighting && (persist < 0.9999F) &&
                    (r1.trans > 0.05F)) {
                  /* don't combine transparent surfaces */
                  *pixel = last_pixel;
                } else {
                  persist = persist * r1.trans;
                }
              }
            }

            if (i < 0) { /* nothing hit */
              break;
            } else {
              if (perspective) {
                if (r1.prim->type != cPrimCharacter) {
                  float extend = r1.dist + 0.00001F;
                  scale3f(r1.dir, extend, nudge);
                } else {
                  float extend = r1.dist;
                  scale3f(r1.dir, extend, nudge);
                }
              }
              last_pixel = *pixel;
              exclude2 = exclude1;
              exclude1 = i;
              pass++;
            }

          } /* end of ray while */

          if (blend_colors) {
            float red_min = _0;
            float green_min = _0;
            float blue_min = _0;
            float red_part;
            float green_part;
            float blue_part;

            if (I->BigEndian) {
              fc[0] = (float)(0xFF & (*pixel >> 24));
              fc[1] = (float)(0xFF & (*pixel >> 16));
              fc[2] = (float)(0xFF & (*pixel >> 8));
              cc3 = (0xFF & (*pixel));
            } else {
              cc3 = (0xFF & (*pixel >> 24));
              fc[2] = (float)(0xFF & (*pixel >> 16));
              fc[1] = (float)(0xFF & (*pixel >> 8));
              fc[0] = (float)(0xFF & (*pixel));
            }

            red_part = red_blend * fc[0];
            green_part = green_blend * fc[1];
            blue_part = blue_blend * fc[2];

            red_min = (green_part > blue_part) ? green_part : blue_part;
            green_min = (red_part > blue_part) ? red_part : blue_part;
            blue_min = (green_part > red_part) ? green_part : red_part;

            if (fc[0] < red_min) fc[0] = red_min;
            if (fc[1] < green_min) fc[1] = green_min;
            if (fc[2] < blue_min) fc[2] = blue_min;

            cc0 = (uint)(fc[0]);
            cc1 = (uint)(fc[1]);
            cc2 = (uint)(fc[2]);

            if (cc0 > 255) cc0 = 255;
            if (cc1 > 255) cc1 = 255;
            if (cc2 > 255) cc2 = 255;

            if (I->BigEndian)
              *pixel = (cc0 << 24) | (cc1 << 16) | (cc2 << 8) | cc3;
            else
              *pixel = (cc3 << 24) | (cc2 << 16) | (cc1 << 8) | cc0;
          }

          {
            /* If final pixel has alpha, then the background should be
               blended into it */
            unsigned char* pixel_c = (unsigned char*)pixel;
            if (pixel_c[3] < 255) {
              float pixa = (pixel_c[3] / 255.f);
              float pixam1 = _1 - pixa;
              float pixam1255 = pixam1 * 255.f;
              if (T->bkrd_data || T->bkrd_is_gradient) {
                pixel_c[0] = (unsigned char)pymol_roundf(pixel_c[0] * pixa +
                                                         bkrd[0] * pixam1255);
                pixel_c[1] = (unsigned char)pymol_roundf(pixel_c[1] * pixa +
                                                         bkrd[1] * pixam1255);
                pixel_c[2] = (unsigned char)pymol_roundf(pixel_c[2] * pixa +
                                                         bkrd[2] * pixam1255);
                pixel_c[3] = (unsigned char)pymol_roundf(pixel_c[3] * pixa +
                                                         bkrd[3] * pixam1255);
              }
            }
          }

          if (!T->edging) break;
          /* if here, then we're edging...
             so accumulate averages */
          {
            unsigned char* pixel_c = (unsigned char*)pixel;
            unsigned int c1, c2, c3, c4;

            edge_avg[0] += (c1 = pixel_c[0]);
            edge_avg[1] += (c2 = pixel_c[1]);
            edge_avg[2] += (c3 = pixel_c[2]);
            edge_avg[3] += (c4 = pixel_c[3]);

            edge_alpha_avg[0] += c1 * c4;
            edge_alpha_avg[1] += c2 * c4;
            edge_alpha_avg[2] += c3 * c4;
            edge_alpha_avg[3] += c4;

            edge_cnt++;
          }

        } /* end of edging while */
        pixel++;
      } /* end of for */
    }
    /* end of if */
  } /* end of for */
  /*  if(T->n_thread>1)
     printf(" Ray: Thread %d: Complete.\n",T->phase+1); */
  return (n_hit);
}

int RayAntiThread(CRayAntiThreadInfo* T) {
  int src_row_pixels;

  unsigned int* pSrc;
  unsigned int* pDst;
  /*   unsigned int m00FF=0x00FF,mFF00=0xFF00,mFFFF=0xFFFF; */
  int width;
  int height;
  int x, y, yy;
  unsigned int* p;
  int offset = 0;
  CRay* I = T->ray;

  if (!T->no_busy) OrthoBusyFast(I->G, 9, 10);
  width = (T->width / T->mag) - 2;
  height = (T->height / T->mag) - 2;

  src_row_pixels = T->width;

  offset = (T->phase * height) / T->n_thread;
  offset = offset - (offset % T->n_thread) + T->phase;

  for (yy = 0; yy < height; yy++) {
    y = (yy + offset) % height; /* make sure threads write to different pages */

    if ((y % T->n_thread) == T->phase) { /* this is my scan line */
      unsigned long c1, c2, c3, c4, a;
      unsigned char* c;

      pSrc = T->image + src_row_pixels * (y * T->mag);
      pDst = T->image_copy + width * y;
      switch (T->mag) {
        case 2: {
          for (x = 0; x < width; x++) {
            c = (unsigned char*)(p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 13);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if (c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else { /* compute straight RGB average */

              c = (unsigned char*)(p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += 13 * c[0];
              c2 += 13 * c[1];
              c3 += 13 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 6;
              c2 = c2 >> 6;
              c3 = c3 >> 6;
            }

            c = (unsigned char*)(pDst++);

            *(c++) = (unsigned char)c1;
            *(c++) = (unsigned char)c2;
            *(c++) = (unsigned char)c3;
            *(c++) = (unsigned char)(c4 >> 6);
          }
        } break;
        case 3: {
          for (x = 0; x < width; x++) {
            c = (unsigned char*)(p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if (c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else { /* compute straight RGB average */

              c = (unsigned char*)(p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 6;
              c2 = c2 >> 6;
              c3 = c3 >> 6;
            }

            c = (unsigned char*)(pDst++);

            *(c++) = (unsigned char)c1;
            *(c++) = (unsigned char)c2;
            *(c++) = (unsigned char)c3;
            *(c++) = (unsigned char)(c4 >> 6);
          }
        } break;
        case 4: {
          for (x = 0; x < width; x++) {
            c = (unsigned char*)(p = pSrc + (x * T->mag));
            c1 = c2 = c3 = c4 = a = 0;

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 8);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 7);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3] * 5);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            c = (unsigned char*)(p += src_row_pixels);

            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;
            c4 += (a = c[3]);
            c1 += c[0] * a;
            c2 += c[1] * a;
            c3 += c[2] * a;
            c += 4;

            if (c4) {
              c1 /= c4;
              c2 /= c4;
              c3 /= c4;
            } else { /* compute straight RGB average */

              c = (unsigned char*)(p = pSrc + (x * T->mag));
              c1 = c2 = c3 = 0;

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 8 * c[0];
              c2 += 8 * c[1];
              c3 += 8 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 7 * c[0];
              c2 += 7 * c[1];
              c3 += 7 * c[2];
              c += 4;
              c1 += 5 * c[0];
              c2 += 5 * c[1];
              c3 += 5 * c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c = (unsigned char*)(p += src_row_pixels);

              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;
              c1 += c[0];
              c2 += c[1];
              c3 += c[2];
              c += 4;

              c1 = c1 >> 7;
              c2 = c2 >> 7;
              c3 = c3 >> 7;
            }

            c = (unsigned char*)(pDst++);

            *(c++) = (unsigned char)c1;
            *(c++) = (unsigned char)c2;
            *(c++) = (unsigned char)c3;
            *(c++) = (unsigned char)(c4 >> 7);
          }
        } break;
      }
    }
  }
  return 1;
}
