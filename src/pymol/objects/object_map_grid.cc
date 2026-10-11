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
-*
-*
-*
Z* -------------------------------------------------------------------
*/

// Operations on map grids: statistics, interpolation, trimming, resampling and
// borders.

#include "pymol/objects/object_map.h"
#include "pymol/base/feedback.h"

/* MapState::GetExcludedStats --
 * PARARMS
 *   G, usual PyMOL globals
 *   ms, MapState
 *   vert_vla, variable length array of vertices
 *   beyond, radius of exclusion
 *   within, radius of inlcusion
 *   level, map level
 *
 * RETURNS
 *  number of exluded pts?
 */
int ObjectMapStateGetExcludedStats(PyMOLGlobals* G, ObjectMapState* ms,
                                   float* vert_vla, float beyond, float within,
                                   float* level) {
  double sum = 0.0, sumsq = 0.0;
  float mean, stdev;
  int cnt = 0;
  int list_size;
  float cutoff = beyond;
  MapType* voxelmap = nullptr;

  /* size of the VLA */
  if (vert_vla) {
    list_size = VLAGetSize(vert_vla) / 3;
  } else {
    list_size = 0;
  }
  if (cutoff < within) cutoff = within;

  /* make a new map from the VLA .............. */
  if (list_size)
    voxelmap = new MapType(G, -cutoff, vert_vla, list_size, nullptr);

  if (voxelmap || (!list_size)) {
    int a, b, c;
    int h, k, l, i, j;
    const int* fdim = ms->FDim;
    int within_flag, within_default = false;
    int beyond_flag;

    const Isofield* field = ms->Field.get();
    if (list_size) MapSetupExpress(voxelmap);

    within_flag = true;
    beyond_flag = true;

    if (within < R_SMALL4) within_default = true;
    for (c = 0; c < fdim[2]; c++) {
      for (b = 0; b < fdim[1]; b++) {
        for (a = 0; a < fdim[0]; a++) {
          if (list_size) {
            within_flag = within_default;
            beyond_flag = true;

            const float* v = F4Ptr(field->points, a, b, c, 0);

            MapLocus(voxelmap, v, &h, &k, &l);
            i = *(MapEStart(voxelmap, h, k, l));
            if (i) {
              j = voxelmap->EList[i++];
              while (j >= 0) {
                if (!within_flag) {
                  if (within3f(vert_vla + 3 * j, v, within)) {
                    within_flag = true;
                  }
                }
                if (within3f(vert_vla + 3 * j, v, beyond)) {
                  beyond_flag = false;
                  break;
                }
                j = voxelmap->EList[i++];
              }
            }
          }

          if (within_flag &&
              beyond_flag) { /* point isn't too close to any vertex */
            const float f_val = F3(field->data, a, b, c);
            sum += f_val;
            sumsq += (f_val * f_val);
            cnt++;
          }
        }
      }
    }
    if (voxelmap) MapFree(voxelmap);
  }
  if (cnt) {
    mean = (float)(sum / cnt);
    stdev = (float)sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
    level[1] = mean;
    level[0] = mean - stdev;
    level[2] = mean + stdev;
  }
  return cnt;
}

/* MapState::ObjectMapStateGetDataRange -- get min/max from the scalar field
 * PARAMS
 * G
 *    usual PyMOLGlobals
 * ms
 *    I
 * RETURNS
 * npts, but stores range in min/max
 */
int ObjectMapStateGetDataRange(PyMOLGlobals* G, ObjectMapState* ms, float* min,
                               float* max) {
  float max_val = 0.0F, min_val = 0.0F;
  CField* data = ms->Field->data.get();
  int cnt = data->dim[0] * data->dim[1] * data->dim[2];
  float* raw_data = (float*)data->data.data();
  if (cnt) {
    int a;
    min_val = (max_val = *(raw_data++));
    for (a = 1; a < cnt; a++) {
      double f_val = *(raw_data++);
      if (min_val > f_val) min_val = f_val;
      if (max_val < f_val) max_val = f_val;
    }
  }
  *min = min_val;
  *max = max_val;
  return cnt;
}

/* MapState::ObjectMapStateGetHistogram -- compute a map histogram
 *
 * INPUT PARAMS
 *
 * h_points - number of histogram points
 * limit - limit the data to (mean - limit * stdev, mean + limit * stdev)
 *         if limit <= 0, don't trim the histogram
 * min_arg, max_arg - limit the data, has priority over "limit" param;
 *         Unused if min_arg==max_arg
 *
 * OUTPUT PARAMS
 *
 * histogram - output histogram buffer. first four values are
 *             minimum, maximum, mean and stdev, followed by n_points
 *             of non-normalized histogram counts.
 */
int ObjectMapStateGetHistogram(PyMOLGlobals* G, ObjectMapState* ms,
                               int n_points, float limit, float* histogram,
                               float min_arg, float max_arg) {
  float max_val = 0.0f, min_val = 0.0f;
  float sum = 0.0f, sumsq = 0.0f;
  float min_his, max_his, irange, mean, stdev;
  int pos;
  CField* data = ms->Field->data.get();
  int cnt = data->dim[0] * data->dim[1] * data->dim[2];
  float* raw_data = (float*)data->data.data();
  if (cnt) {
    int a;

    // compute min/max/mean/stdev
    sum = min_val = (max_val = *(raw_data++));
    sumsq = sum * sum;
    for (a = 1; a < cnt; a++) {
      double f_val = *(raw_data++);
      if (min_val > f_val) min_val = f_val;
      if (max_val < f_val) max_val = f_val;
      sum += f_val;
      sumsq += f_val * f_val;
    }
    mean = (float)(sum / cnt);
    stdev = (float)sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));

    // adjust min/max to limit
    if (min_arg != max_arg) {
      min_his = min_arg;
      max_his = max_arg;
    } else if (limit > 0.0F) {
      min_his = mean - limit * stdev;
      if (min_his < min_val) min_his = min_val;
      max_his = mean + limit * stdev;
      if (max_his > max_val) max_his = max_val;
    } else {
      min_his = min_val;
      max_his = max_val;
    }

    // Compute the histogram
    if (n_points > 0) {
      irange = (float)(n_points - 1) / (max_his - min_his);
      for (a = 0; a < n_points; a++) histogram[a + 4] = 0.0f;
      raw_data = (float*)data->data.data();
      for (a = 0; a < cnt; a++) {
        double f_val = *(raw_data++);
        pos = (int)(irange * (f_val - min_his));
        if (pos >= 0 && pos < n_points) {
          histogram[pos + 4] += 1.0;
        }
      }
    }
    histogram[0] = min_his;
    histogram[1] = max_his;
    histogram[2] = mean;
    histogram[3] = stdev;
  } else {
    histogram[0] = 0.0;
    histogram[1] = 1.0;
    histogram[2] = 1.0;
    histogram[3] = 1.0;
  }
  return cnt;
}

/**
 * @see ObjectMapStateInterpolate for parameter description
 * @param state Object state (can be -2 for current state)
 * @return True if the map state exists and `result` has been written to
 */
int ObjectMapInterpolate(ObjectMap* I, int state, const float* array,
                         float* result, std::uint8_t* flag, int n) {
  int ok = false;
  float txf_buffer[3];
  float* txf = txf_buffer;

  ObjectMapState* ms = ObjectMapGetState(I, state);

  if (ms) {
    const double* matrix = ObjectStateGetInvMatrix(ms);

    if (matrix) {
      /* we have to back-transform points */
      if (n > 1) {
        txf = pymol::malloc<float>(3 * n);
      }

      const float* src = array;
      float* dst = txf;
      array = txf;

      for (int nn = n; nn--; src += 3, dst += 3) {
        transform44d3f(matrix, src, dst);
      }
    }

    ok = ObjectMapStateInterpolate(ms, array, result, flag, n);
    ok = true;
  }

  if (txf != txf_buffer) FreeP(txf);

  return (ok);
}

static void ObjectMapStateTrim(PyMOLGlobals* G, ObjectMapState* ms, float* mn,
                               float* mx, int quiet) {
  int div[3];
  int min[3];
  int fdim[4];
  int new_min[3], new_max[3], new_fdim[3];
  int a, b, c, d, e, f;
  float *vt, *vr;
  float v[3];
  float grid[3];
  Isofield* field;
  float orig_size = 1.0F;
  float new_size = 1.0F;

  if (ObjectMapStateValidXtal(ms)) {
    float tst[3], frac_tst[3];
    float frac_mn[3];
    float frac_mx[3];
    int hit_flag = false;

    /* compute the limiting box extents in fractional space */

    for (a = 0; a < 8; a++) {
      tst[0] = (a & 0x1) ? mn[0] : mx[0];
      tst[1] = (a & 0x2) ? mn[1] : mx[1];
      tst[2] = (a & 0x4) ? mn[2] : mx[2];
      transform33f3f(ms->Symmetry->Crystal.realToFrac(), tst, frac_tst);
      if (!a) {
        copy3f(frac_tst, frac_mn);
        copy3f(frac_tst, frac_mx);
      } else {
        for (b = 0; b < 3; b++) {
          frac_mn[b] = (frac_mn[b] > frac_tst[b]) ? frac_tst[b] : frac_mn[b];
          frac_mx[b] = (frac_mx[b] < frac_tst[b]) ? frac_tst[b] : frac_mx[b];
        }
      }
    }
    for (a = 0; a < 3; a++) {
      div[a] = ms->Div[a];
      min[a] = ms->Min[a];
      fdim[a] = ms->FDim[a];
    }
    fdim[3] = 3;

    {
      int first_flag[3] = {false, false, false};
      for (d = 0; d < 3; d++) {
        int tst_min, tst_max;
        float v_min, v_max;
        for (c = 0; c < (fdim[d] - 1); c++) {
          tst_min = c + min[d];
          tst_max = tst_min + 1;
          v_min = tst_min / ((float)div[d]);
          v_max = tst_max / ((float)div[d]);
          if (((v_min >= frac_mn[d]) && (v_min <= frac_mx[d])) ||
              ((v_max >= frac_mn[d]) && (v_max <= frac_mx[d]))) {
            if (!first_flag[d]) {
              first_flag[d] = true;
              new_min[d] = tst_min;
              new_max[d] = tst_max;
            } else {
              new_min[d] = (new_min[d] > tst_min) ? tst_min : new_min[d];
              new_max[d] = (new_max[d] < tst_max) ? tst_max : new_max[d];
            }
          }
        }
        new_fdim[d] = (new_max[d] - new_min[d]) + 1;
      }
      hit_flag = first_flag[0] && first_flag[1] & first_flag[2];
    }

    if (hit_flag)
      hit_flag = (new_fdim[0] != fdim[0]) || (new_fdim[1] != fdim[1]) ||
                 (new_fdim[2] != fdim[2]);

    if (hit_flag) {
      orig_size = fdim[0] * fdim[1] * fdim[2];
      new_size = new_fdim[0] * new_fdim[1] * new_fdim[2];

      field = new Isofield(G, new_fdim);
      field->save_points = ms->Field->save_points;

      for (c = 0; c < new_fdim[2]; c++) {
        f = c + (new_min[2] - min[2]);
        for (b = 0; b < new_fdim[1]; b++) {
          e = b + (new_min[1] - min[1]);
          for (a = 0; a < new_fdim[0]; a++) {
            d = a + (new_min[0] - min[0]);
            vt = F4Ptr(field->points, a, b, c, 0);
            vr = F4Ptr(ms->Field->points, d, e, f, 0);
            copy3f(vr, vt);
            F3(field->data, a, b, c) = F3(ms->Field->data, d, e, f);
          }
        }
      }
      for (a = 0; a < 3; a++) {
        ms->Min[a] = new_min[a];
        ms->Max[a] = new_max[a];
        ms->FDim[a] = new_fdim[a];
      }
      ms->Field.reset(field);

      /* compute new extents */
      v[2] = (ms->Min[2]) / ((float)ms->Div[2]);
      v[1] = (ms->Min[1]) / ((float)ms->Div[1]);
      v[0] = (ms->Min[0]) / ((float)ms->Div[0]);

      auto const& cryst = ms->Symmetry->Crystal;

      transform33f3f(cryst.fracToReal(), v, ms->ExtentMin);

      v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
      v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
      v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);

      transform33f3f(cryst.fracToReal(), v, ms->ExtentMax);

      /* new corner */
      {
        float vv[3];
        d = 0;
        for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
          for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
            for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
              transform33f3f(cryst.fracToReal(), v, vv);
              copy3f(vv, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }
    }
  } else { /* not a crystal map */
    int hit_flag = false;
    float* origin = ms->Origin.data();

    for (a = 0; a < 3; a++) {
      min[a] = ms->Min[a];
      grid[a] = ms->Grid[a];
      fdim[a] = ms->FDim[a];
    }
    fdim[3] = 3;

    {
      int first_flag[3] = {false, false, false};
      for (d = 0; d < 3; d++) {
        int tst_min, tst_max;
        float v_min, v_max;
        for (c = 0; c < (fdim[d] - 1); c++) {
          tst_min = c + min[d];
          tst_max = tst_min + 1;
          v_min = origin[d] + grid[d] * (tst_min);
          v_max = origin[d] + grid[d] * (tst_max);
          if (((v_min >= mn[d]) && (v_min <= mx[d])) ||
              ((v_max >= mn[d]) && (v_max <= mx[d]))) {
            if (!first_flag[d]) {
              first_flag[d] = true;
              hit_flag = true;
              new_min[d] = tst_min;
              new_max[d] = tst_max;
            } else {
              new_min[d] = (new_min[d] > tst_min) ? tst_min : new_min[d];
              new_max[d] = (new_max[d] < tst_max) ? tst_max : new_max[d];
            }
          }
        }
        new_fdim[d] = (new_max[d] - new_min[d]) + 1;
      }
      hit_flag = first_flag[0] && first_flag[1] & first_flag[2];
    }

    if (hit_flag)
      hit_flag = (new_fdim[0] != fdim[0]) || (new_fdim[1] != fdim[1]) ||
                 (new_fdim[2] != fdim[2]);

    if (hit_flag) {
      orig_size = fdim[0] * fdim[1] * fdim[2];
      new_size = new_fdim[0] * new_fdim[1] * new_fdim[2];

      field = new Isofield(G, new_fdim);
      field->save_points = ms->Field->save_points;

      for (c = 0; c < new_fdim[2]; c++) {
        f = c + (new_min[2] - min[2]);
        for (b = 0; b < new_fdim[1]; b++) {
          e = b + (new_min[1] - min[1]);
          for (a = 0; a < new_fdim[0]; a++) {
            d = a + (new_min[0] - min[0]);
            vt = F4Ptr(field->points, a, b, c, 0);
            vr = F4Ptr(ms->Field->points, d, e, f, 0);
            copy3f(vr, vt);
            F3(field->data, a, b, c) = F3(ms->Field->data, d, e, f);
          }
        }
      }
      for (a = 0; a < 3; a++) {
        ms->Min[a] = new_min[a];
        ms->Max[a] = new_max[a];
        ms->FDim[a] = new_fdim[a];
        if (!ms->Dim.empty()) ms->Dim[a] = new_fdim[a];
      }

      ms->Field.reset(field);

      for (e = 0; e < 3; e++) {
        ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
        ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
      }

      d = 0;
      for (c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
        v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);

        for (b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
          v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

          for (a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
            v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
            copy3f(v, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  }
  if (!quiet) {
    PRINTFB(G, FB_ObjectMap, FB_Actions)
    " ObjectMap: Map volume reduced by %2.0f%%.\n",
        (100 * (orig_size - new_size)) / orig_size ENDFB(G);
  }
}

static void ObjectMapStateDouble(PyMOLGlobals* G, ObjectMapState* ms) {
  int div[3];
  int min[3];
  int max[3];
  int fdim[4];
  int a, b, c;
  float v[3], vr[3];
  float* vt;
  float x, y, z;
  float grid[3];

  Isofield* field;

  if (ObjectMapStateValidXtal(ms)) {
    for (a = 0; a < 3; a++) {
      div[a] = ms->Div[a] * 2;
      min[a] = ms->Min[a] * 2;
      max[a] = ms->Max[a] * 2;
      fdim[a] = ms->FDim[a] * 2 - 1;
    }
    fdim[3] = 3;
    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;
    for (c = 0; c < fdim[2]; c++) {
      v[2] = (c + min[2]) / ((float)div[2]);
      z = (c & 0x1) ? 0.5F : 0.0F;
      for (b = 0; b < fdim[1]; b++) {
        v[1] = (b + min[1]) / ((float)div[1]);
        y = (b & 0x1) ? 0.5F : 0.0F;
        for (a = 0; a < fdim[0]; a++) {
          v[0] = (a + min[0]) / ((float)div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          x = (a & 0x1) ? 0.5F : 0.0F;
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(vr, vt);
          if ((a & 0x1) || (b & 0x1) || (c & 0x1)) {
            F3(field->data, a, b, c) = FieldInterpolatef(
                ms->Field->data.get(), a / 2, b / 2, c / 2, x, y, z);
          } else {
            F3(field->data, a, b, c) = F3(ms->Field->data, a / 2, b / 2, c / 2);
          }
        }
      }
    }
    for (a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      ms->Div[a] = div[a];
    }

    ms->Field.reset(field);
  } else {
    for (a = 0; a < 3; a++) {
      grid[a] = ms->Grid[a] / 2.0F;
      min[a] = ms->Min[a] * 2;
      max[a] = ms->Max[a] * 2;
      fdim[a] = ms->FDim[a] * 2 - 1;
    }
    fdim[3] = 3;

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    for (c = 0; c < fdim[2]; c++) {
      v[2] = ms->Origin[2] + grid[2] * (c + min[2]);
      z = (c & 0x1) ? 0.5F : 0.0F;
      for (b = 0; b < fdim[1]; b++) {
        v[1] = ms->Origin[1] + grid[1] * (b + min[1]);
        y = (b & 0x1) ? 0.5F : 0.0F;
        for (a = 0; a < fdim[0]; a++) {
          v[0] = ms->Origin[0] + grid[0] * (a + min[0]);
          x = (a & 0x1) ? 0.5F : 0.0F;
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(v, vt);
          if ((a & 0x1) || (b & 0x1) || (c & 0x1)) {
            F3(field->data, a, b, c) = FieldInterpolatef(
                ms->Field->data.get(), a / 2, b / 2, c / 2, x, y, z);
          } else {
            F3(field->data, a, b, c) = F3(ms->Field->data, a / 2, b / 2, c / 2);
          }
        }
      }
    }
    for (a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      if (!ms->Dim.empty()) ms->Dim[a] = fdim[a];
      if (!ms->Grid.empty()) ms->Grid[a] = grid[a];
    }
    ms->Field.reset(field);
  }
}

static void ObjectMapStateHalve(PyMOLGlobals* G, ObjectMapState* ms,
                                int smooth) {
  int div[3];
  int min[3];
  int max[3];
  int fdim[4];
  int a, b, c;
  float v[3], vr[3];
  float* vt;
  float x, y, z;
  float grid[3];

  Isofield* field;

  if (ObjectMapStateValidXtal(ms)) {
    int *old_div, *old_min, *old_max;
    int a_2, b_2, c_2;

    for (a = 0; a < 3; a++) {
      div[a] = ms->Div[a] / 2;
      /* note that when Div is odd, a one-cell extrapolation will
         occur in order to preserve map size and get us onto a even
         number of sampling intervals for the cell */

      min[a] = ms->Min[a] / 2;
      max[a] = ms->Max[a] / 2;
      while ((min[a] * 2) < ms->Min[a]) min[a]++;
      while ((max[a] * 2) > ms->Max[a]) max[a]--;

      fdim[a] = (max[a] - min[a]) + 1;
    }
    fdim[3] = 3;
    old_div = ms->Div;
    old_min = ms->Min;
    old_max = ms->Max;

    if (smooth) FieldSmooth3f(ms->Field->data.get());

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    /*
       printf("old_div %d %d %d\n",old_div[0],old_div[1],old_div[2]);
       printf("old_min %d %d %d\n",old_min[0],old_min[1],old_min[2]);
       printf("min %d %d %d\n",min[0],min[1],min[2]);
       printf("%d %d %d\n",ms->FDim[0],ms->FDim[1],ms->FDim[2]);
     */

    for (c = 0; c < fdim[2]; c++) {
      v[2] = (c + min[2]) / ((float)div[2]);
      c_2 = (c + min[2]) * 2 - old_min[2];
      z = (v[2] - ((c_2 + old_min[2]) / (float)old_div[2])) * old_div[2];
      if (c_2 >= old_max[2]) {
        c_2 = old_max[2] - 1;
        z = (v[2] - ((c_2 + old_min[2]) / (float)old_div[2])) * old_div[2];
      }
      for (b = 0; b < fdim[1]; b++) {
        v[1] = (b + min[1]) / ((float)div[1]);
        b_2 = (b + min[1]) * 2 - old_min[1];
        y = (v[1] - ((b_2 + old_min[1]) / (float)old_div[1])) * old_div[1];
        if (b_2 >= old_max[1]) {
          b_2 = old_max[1] - 1;
          y = (v[1] - ((b_2 + old_min[1]) / (float)old_div[1])) * old_div[1];
        }
        for (a = 0; a < fdim[0]; a++) {
          v[0] = (a + min[0]) / ((float)div[0]);
          a_2 = (a + min[0]) * 2 - old_min[0];
          x = (v[0] - ((a_2 + old_min[0]) / (float)old_div[0])) * old_div[0];
          if (a_2 >= old_max[0]) {
            a_2 = old_max[0] - 1;
            x = (v[0] - ((a_2 + old_min[0]) / (float)old_div[0])) * old_div[0];
          }
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(vr, vt);
          F3(field->data, a, b, c) =
              FieldInterpolatef(ms->Field->data.get(), a_2, b_2, c_2, x, y, z);
        }
      }
    }
    for (a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      ms->Div[a] = div[a];
    }

    ms->Field.reset(field);

    /* compute new extents */
    v[2] = (ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

    /* new corner */
    {
      float vv[3];
      int d = 0;

      for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
        v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
        for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
          v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
          for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
            v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
            transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vv);
            copy3f(vv, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  } else {
    for (a = 0; a < 3; a++) {
      grid[a] = ms->Grid[a] * 2.0F;
      min[a] = ms->Min[a] / 2;
      max[a] = ms->Max[a] / 2;
      fdim[a] = (ms->FDim[a] + 1) / 2;
    }
    fdim[3] = 3;

    field = new Isofield(G, fdim);
    field->save_points = ms->Field->save_points;

    for (c = 0; c < fdim[2]; c++) {
      v[2] = ms->Origin[2] + grid[2] * (c + min[2]);
      for (b = 0; b < fdim[1]; b++) {
        v[1] = ms->Origin[1] + grid[1] * (b + min[1]);
        for (a = 0; a < fdim[0]; a++) {
          v[0] = ms->Origin[0] + grid[0] * (a + min[0]);
          vt = F4Ptr(field->points, a, b, c, 0);
          copy3f(v, vt);
          F3(field->data, a, b, c) = F3(ms->Field->data, a * 2, b * 2, c * 2);
        }
      }
    }
    for (a = 0; a < 3; a++) {
      ms->Min[a] = min[a];
      ms->Max[a] = max[a];
      ms->FDim[a] = fdim[a];
      if (!ms->Dim.empty()) ms->Dim[a] = fdim[a];
      if (!ms->Grid.empty()) ms->Grid[a] = grid[a];
    }
    ms->Field.reset(field);
  }
}

pymol::Result<> ObjectMapTrim(ObjectMap* I, int state, float* mn, float* mx,
                              int quiet) {
  int update = false;

  /* TO DO: convert mn and mx into map local coordinates if map itself is
   * transformed...  */

  if (state < 0) {
    for (auto& ms : I->State) {
      if (ms.Active) {
        ObjectMapStateTrim(I->G, &ms, mn, mx, quiet);
        update = true;
      }
    }
  } else if ((state >= 0) && (state < I->State.size()) &&
             (I->State[state].Active)) {
    ObjectMapStateTrim(I->G, &I->State[state], mn, mx, quiet);
  } else {
    return pymol::make_error("Invalid state.");
  }
  if (update) ObjectMapUpdateExtents(I);
  return {};
}

pymol::Result<> ObjectMapDouble(ObjectMap* I, int state) {
  if (state < 0) {
    for (auto& state : I->State) {
      if (state.Active) ObjectMapStateDouble(I->G, &state);
    }
  } else if ((state >= 0) && (state < I->State.size()) &&
             (I->State[state].Active)) {
    ObjectMapStateDouble(I->G, &I->State[state]);
  } else {
    return pymol::make_error("Invalidate state.");
  }
  return {};
}

pymol::Result<> ObjectMapHalve(ObjectMap* I, int state, int smooth) {
  if (state < 0) {
    for (auto& state : I->State) {
      if (state.Active) ObjectMapStateHalve(I->G, &state, smooth);
    }
  } else if ((state >= 0) && (state < I->State.size()) &&
             (I->State[state].Active)) {
    ObjectMapStateHalve(I->G, &I->State[state], smooth);
  } else {
    return pymol::make_error("Invalidate state.");
  }
  ObjectMapUpdateExtents(I);

  return {};
}

int ObjectMapStateContainsPoint(ObjectMapState* ms, float* point) {
  int result = false;
  float x, y, z;
  int x_floor, y_floor, z_floor;
  int x_ceil, y_ceil, z_ceil;

  if (ObjectMapStateValidXtal(ms)) {
    float frac[3];

    transform33f3f(ms->Symmetry->Crystal.realToFrac(), point, frac);

    x = (ms->Div[0] * frac[0]);
    y = (ms->Div[1] * frac[1]);
    z = (ms->Div[2] * frac[2]);
    x_floor = floor(x);
    x_ceil = ceil(x);
    y_floor = floor(y);
    y_ceil = ceil(y);
    z_floor = floor(z);
    z_ceil = ceil(z);

    if ((x_floor >= ms->Min[0]) && (x_ceil <= ms->Max[0]) &&
        (y_floor >= ms->Min[1]) && (y_ceil <= ms->Max[1]) &&
        (z_floor >= ms->Min[2]) && (z_ceil <= ms->Max[2]))
      result = true;
  } else {
    x = (point[0] - ms->Origin[0]) / ms->Grid[0];
    y = (point[1] - ms->Origin[1]) / ms->Grid[1];
    z = (point[2] - ms->Origin[2]) / ms->Grid[2];

    x_floor = floor(x);
    x_ceil = ceil(x);
    y_floor = floor(y);
    y_ceil = ceil(y);
    z_floor = floor(z);
    z_ceil = ceil(z);

    if ((x_floor >= ms->Min[0]) && (x_ceil <= ms->Max[0]) &&
        (y_floor >= ms->Min[1]) && (y_ceil <= ms->Max[1]) &&
        (z_floor >= ms->Min[2]) && (z_ceil <= ms->Max[2]))
      result = true;

    if ((x >= ms->Min[0]) && (x <= ms->Max[0]) && (y >= ms->Min[1]) &&
        (y <= ms->Max[1]) && (z >= ms->Min[2]) && (z <= ms->Max[2]))
      result = true;
  }
  return (result);
}

/**
 * @param array Coordinate array of length `3*n`
 * @param[out] result Array of length `n` which will be populated with map
 * values at coorinate positions (linear interpolated between grid points)
 * @param[out] flag Array of length `n` which will be populated with booleans
 * indicating if points were within map bounds (optional, can be nullptr)
 * @return False if any coordinate was out of bounds
 */
int ObjectMapStateInterpolate(ObjectMapState* ms, const float* array,
                              float* result, std::uint8_t* flag, int n) {
  int ok = true;
  const float* inp;
  int a, b, c;
  float x, y, z;
  inp = array;

  if (ObjectMapStateValidXtal(ms)) {
    float frac[3];

    while (n--) {
      /* get the fractional coordinate */
      transform33f3f(ms->Symmetry->Crystal.realToFrac(), inp, frac);

      inp += 3;

      /* compute the effective lattice offset as a function of cell spacing */

      x = (ms->Div[0] * frac[0]);
      y = (ms->Div[1] * frac[1]);
      z = (ms->Div[2] * frac[2]);

      /* now separate the integral and fractional parts for interpolation */

      a = (int)floor(x + R_SMALL8);
      b = (int)floor(y + R_SMALL8);
      c = (int)floor(z + R_SMALL8);
      x -= a;
      y -= b;
      z -= c;

      if (flag) *flag = 1;
      /* wrap into the map */

      if (a < ms->Min[0]) {
        if (x < 0.99F) {
          ok = false;
          if (flag) *flag = 0;
        }
        x = 0.0F;
        a = ms->Min[0];
      } else if (a >= ms->FDim[0] + ms->Min[0] - 1) {
        if (x > 0.01F) {
          ok = false;
          if (flag) *flag = 0;
        }
        x = 0.0F;
        a = ms->FDim[0] + ms->Min[0] - 1;
      }

      if (b < ms->Min[1]) {
        if (y < 0.99F) {
          ok = false;
          if (flag) *flag = 0;
        }
        y = 0.0F;
        b = ms->Min[1];
      } else if (b >= ms->FDim[1] + ms->Min[1] - 1) {
        if (y > 0.01F) {
          ok = false;
          if (flag) *flag = 0;
        }
        y = 0.0F;
        b = ms->FDim[1] + ms->Min[1] - 1;
      }

      if (c < ms->Min[2]) {
        if (z < 0.99F) {
          ok = false;
          if (flag) *flag = 0;
        }
        z = 0.0F;
        c = ms->Min[2];
      } else if (c >= ms->FDim[2] + ms->Min[2] - 1) {
        if (z > 0.01) {
          ok = false;
          if (flag) *flag = 0;
        }
        z = 0.0F;
        c = ms->FDim[2] + ms->Min[2] - 1;
      }
      /*      printf("%d %d %d %8.3f %8.3f %8.3f\n",a,b,c,x,y,z); */
      *(result++) = FieldInterpolatef(ms->Field->data.get(), a - ms->Min[0],
                                      b - ms->Min[1], c - ms->Min[2], x, y, z);
      if (flag) flag++;
    }
  } else {
    while (n--) {
      x = (inp[0] - ms->Origin[0]) / ms->Grid[0];
      y = (inp[1] - ms->Origin[1]) / ms->Grid[1];
      z = (inp[2] - ms->Origin[2]) / ms->Grid[2];
      inp += 3;

      a = (int)floor(x + R_SMALL8);
      b = (int)floor(y + R_SMALL8);
      c = (int)floor(z + R_SMALL8);
      x -= a;
      y -= b;
      z -= c;

      if (flag) *flag = 1;
      if (a < ms->Min[0]) {
        x = 0.0F;
        a = ms->Min[0];
        ok = false;
        if (flag) *flag = 0;
      } else if (a >= ms->Max[0]) {
        x = 1.0F;
        a = ms->Max[0] - 1;
        ok = false;
        if (flag) *flag = 0;
      }

      if (b < ms->Min[1]) {
        y = 0.0F;
        b = ms->Min[1];
        ok = false;
        if (flag) *flag = 0;
      } else if (b >= ms->Max[1]) {
        y = 1.0F;
        b = ms->Max[1] - 1;
        ok = false;
        if (flag) *flag = 0;
      }

      if (c < ms->Min[2]) {
        z = 0.0F;
        c = ms->Min[2];
        ok = false;
        if (flag) *flag = 0;
      } else if (c >= ms->Max[2]) {
        z = 1.0F;
        c = ms->Max[2] - 1;
        ok = false;
        if (flag) *flag = 0;
      }
      /*      printf("%d %d %d %8.3f %8.3f %8.3f\n",a,b,c,x,y,z); */
      *(result++) = FieldInterpolatef(ms->Field->data.get(), a - ms->Min[0],
                                      b - ms->Min[1], c - ms->Min[2], x, y, z);
      if (flag) flag++;
    }
  }
  return (ok);
}

void ObjectMapStateClamp(ObjectMapState* I, float clamp_floor,
                         float clamp_ceiling) {
  int a, b, c;
  float* fp;

  for (a = 0; a < I->FDim[0]; a++)
    for (b = 0; b < I->FDim[1]; b++)
      for (c = 0; c < I->FDim[2]; c++) {
        fp = F3Ptr(I->Field->data, a, b, c);
        if (*fp < clamp_floor)
          *fp = clamp_floor;
        else if (*fp > clamp_ceiling)
          *fp = clamp_ceiling;
      }
}

int ObjectMapStateSetBorder(ObjectMapState* I, float level) {
  int result = true;
  int a, b, c;

  c = I->FDim[2] - 1;
  for (a = 0; a < I->FDim[0]; a++)
    for (b = 0; b < I->FDim[1]; b++) {
      F3(I->Field->data, a, b, 0) = level;
      F3(I->Field->data, a, b, c) = level;
    }

  a = I->FDim[0] - 1;
  for (b = 0; b < I->FDim[1]; b++)
    for (c = 0; c < I->FDim[2]; c++) {
      F3(I->Field->data, 0, b, c) = level;
      F3(I->Field->data, a, b, c) = level;
    }

  b = I->FDim[1] - 1;
  for (a = 0; a < I->FDim[0]; a++)
    for (c = 0; c < I->FDim[2]; c++) {
      F3(I->Field->data, a, 0, c) = level;
      F3(I->Field->data, a, b, c) = level;
    }
  return (result);
}
int ObjectMapSetBorder(ObjectMap* I, float level, int state) {
  for (StateIterator iter(I, state); iter.next();) {
    auto* oms = &I->State[iter.state];
    if (oms->Active) {
      if (!ObjectMapStateSetBorder(oms, level)) return false;
    }
  }
  return true;
}
