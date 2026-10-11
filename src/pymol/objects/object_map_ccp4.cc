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

// Reading and writing CCP4 maps.

#include "pymol/base/err.h"
#include "pymol/objects/object_map.h"
#include "pymol/geometry/matrix.h"
#include "pymol/base/file.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"

#define n_space_group_numbers 231
static const char* space_group_numbers[] = {
    "",           "P 1",        "P -1",       "P 2",        "P 21",
    "C 2",        "P M",        "P C",        "C M",        "C C",
    "P 2/M",      "P 21/M",     "C 2/M",      "P 2/C",      "P 21/C",
    "C 2/C",      "P 2 2 2",    "P 2 2 21",   "P 21 21 2",  "P 21 21 21",
    "C 2 2 21",   "C 2 2 2",    "F 2 2 2",    "I 2 2 2",    "I 21 21 21",
    "P M M 2",    "P M C 21",   "P C C 2",    "P M A 2",    "P C A 21",
    "P N C 2",    "P M N 21",   "P B A 2",    "P N A 21",   "P N N 2",
    "C M M 2",    "C M C 21",   "C C C 2",    "A M M 2",    "A B M 2",
    "A M A 2",    "A B A 2",    "F M M 2",    "F D D 2",    "I M M 2",
    "I B A 2",    "I M A 2",    "P M M M",    "P N N N",    "P C C M",
    "P B A N",    "P M M A",    "P N N A",    "P M N A",    "P C C A",
    "P B A M",    "P C C N",    "P B C M",    "P N N M",    "P M M N",
    "P B C N",    "P B C A",    "P N M A",    "C M C M",    "C M C A",
    "C M M M",    "C C C M",    "C M M A",    "C C C A",    "F M M M",
    "F D D D",    "I M M M",    "I B A M",    "I B C A",    "I M M A",
    "P 4",        "P 41",       "P 42",       "P 43",       "I 4",
    "I 41",       "P -4",       "I -4",       "P 4/M",      "P 42/M",
    "P 4/N",      "P 42/N",     "I 4/M",      "I 41/A",     "P 4 2 2",
    "P 4 21 2",   "P 41 2 2",   "P 41 21 2",  "P 42 2 2",   "P 42 21 2",
    "P 43 2 2",   "P 43 21 2",  "I 4 2 2",    "I 41 2 2",   "P 4 M M",
    "P 4 B M",    "P 42 C M",   "P 42 N M",   "P 4 C C",    "P 4 N C",
    "P 42 M C",   "P 42 B C",   "I 4 M M",    "I 4 C M",    "I 41 M D",
    "I 41 C D",   "P -4 2 M",   "P -4 2 C",   "P -4 21 M",  "P -4 21 C",
    "P -4 M 2",   "P -4 C 2",   "P -4 B 2",   "P -4 N 2",   "I -4 M 2",
    "I -4 C 2",   "I -4 2 M",   "I -4 2 D",   "P 4/M M M",  "P 4/M C C",
    "P 4/N B M",  "P 4/N N C",  "P 4/M B M",  "P 4/M N C",  "P 4/N M M",
    "P 4/N C C",  "P 42/M M C", "P 42/M C M", "P 42/N B C", "P 42/N N M",
    "P 42/M B C", "P 42/M N M", "P 42/N M C", "P 42/N C M", "I 4/M M M",
    "I 4/M C M",  "I 41/A M D", "I 41/A C D", "P 3",        "P 31",
    "P 32",       "R 3",        "P -3",       "R -3",       "P 3 1 2",
    "P 3 2 1",    "P 31 1 2",   "P 31 2 1",   "P 32 1 2",   "P 32 2 1",
    "R 3 2",      "P 3 M 1",    "P 3 1 M",    "P 3 C 1",    "P 3 1 C",
    "R 3 M",      "R 3 C",      "P -3 1 M",   "P -3 1 C",   "P -3 M 1",
    "P -3 C 1",   "R -3 M",     "R -3 C",     "P 6",        "P 61",
    "P 65",       "P 62",       "P 64",       "P 63",       "P -6",
    "P 6/M",      "P 63/M",     "P 6 2 2",    "P 61 2 2",   "P 65 2 2",
    "P 62 2 2",   "P 64 2 2",   "P 63 2 2",   "P 6 M M",    "P 6 C C",
    "P 63 C M",   "P 63 M C",   "P -6 M 2",   "P -6 C 2",   "P -6 2 M",
    "P -6 2 C",   "P 6/M M M",  "P 6/M C C",  "P 63/M C M", "P 63/M M C",
    "P 2 3",      "F 2 3",      "I 2 3",      "P 21 3",     "I 21 3",
    "P M -3",     "P N -3",     "F M -3",     "F D -3",     "I M -3",
    "P A -3",     "I A -3",     "P 4 3 2",    "P 42 3 2",   "F 4 3 2",
    "F 41 3 2",   "I 4 3 2",    "P 43 3 2",   "P 41 3 2",   "I 41 3 2",
    "P -4 3 M",   "F -4 3 M",   "I -4 3 M",   "P -4 3 N",   "F -4 3 C",
    "I -4 3 D",   "P M -3 M",   "P N -3 N",   "P M -3 N",   "P N -3 M",
    "F M -3 M",   "F M -3 C",   "F D -3 M",   "F D -3 C",   "I M -3 M",
    "I A -3 D",
};
static float ccp4_next_value(char** pp, int mode) {
  char* p = *pp;
  switch (mode) {
    case 0:
      *pp += 1;
      return (float)*((int8_t*)p);
    case 1:
      *pp += 2;
      return (float)*((int16_t*)p);
    case 2:
      *pp += 4;
      return *((float*)p);
  }
  printf("ERROR unsupported mode\n");
  return 0.f;
}

/* swaps n*width bytes in memory starting at p */
static void swap_endian(char* p, int n, int width) {
  char tmp, *q, *pstop = p + (n - 1) * width + 1;
  int w2 = width / 2, wm1 = width - 1;
  for (; p < pstop; p += w2) {
    for (q = p + wm1; p < q; p++, q--) {
      tmp = *p;
      *p = *q;
      *q = tmp;
    }
  }
}

static bool validateCCP4LoadType(int& format) {
  switch (format) {
    case cLoadTypeCCP4Map:
      format = cLoadTypeCCP4Str;
    case cLoadTypeCCP4Str:
    case cLoadTypeCCP4Unspecified:
    case cLoadTypeMRC:
      return true;
    case cLoadTypeCCP4UnspecifiedStr:
      format = cLoadTypeCCP4Unspecified;
      return true;
    case cLoadTypeMRCStr:
      format = cLoadTypeMRC;
      return true;
    default:
      return false;
  }
}

static int ObjectMapCCP4StrToMap(ObjectMap* I, char* CCP4Str, int bytes,
                                 int state, int quiet, int format) {
  auto G = I->G;
  char* p;
  int* i;
  size_t bytes_per_pt;
  char* q;
  float dens;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  int little_endian = 1, map_endian;
  /* CCP4 named from their docs */
  int nc, nr, ns;
  int map_mode;
  int ncstart, nrstart, nsstart;
  int nx, ny, nz;
  float xlen, ylen, zlen, alpha, beta, gamma;
  int ispg;  // space group number
  int sym_skip;
  int mapc, mapr, maps;
  int cc[3];
  int n_pts;
  double sum, sumsq;
  float mean, stdev;
  int normalize;
  ObjectMapState* ms;
  int expectation;

  if (!validateCCP4LoadType(format)) {
    ErrMessage(G, __func__, "wrong format");
    return false;
  }

  if (bytes < 256 * sizeof(int)) {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
    " ObjectMapCCP4: Map appears to be truncated -- aborting." ENDFB(I->G);
    return (0);
  }

  /* state check */
  if (state < 0) state = I->State.size();
  /* alloc/init a new MapState */
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_ccp4_maps);

  p = CCP4Str;
  little_endian = *((char*)&little_endian);
  map_endian = (*p || *(p + 1));  // NOTE: this assumes 0x0 < NC < 0x10000

  if (little_endian != map_endian) {
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Blather)
      " ObjectMapCCP4: Map appears to be reverse endian, swapping...\n" ENDFB(
          I->G);
    }
    swap_endian(p, 256, sizeof(int));
  }

  i = (int*)p;
  nc = *(i++); /* columns */
  nr = *(i++); /* rows */
  ns = *(i++); /* sections */
  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: NC %d   NR %d   NS %d\n", nc, nr, ns ENDFB(I->G);
  }
  map_mode = *(i++); /* mode */

  switch (map_mode) {
    case 0:
      bytes_per_pt = 1;  // int8_t
      break;
    case 1:
      bytes_per_pt = 2;  // int16_t
      break;
    case 2:
      bytes_per_pt = 4;  // float
      break;
    default:
      PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      "ObjectMapCCP4-ERR: Only map mode 0-2 currently supported (this map is "
      "mode %d)",
          map_mode ENDFB(I->G);
      return (0);
  }

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: Map is mode %d.\n", map_mode ENDFB(I->G);
  }
  ncstart = *(i++);
  nrstart = *(i++);
  nsstart = *(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: NCSTART %d   NRSTART %d   NSSTART  %d\n", ncstart, nrstart,
        nsstart ENDFB(I->G);
  }

  nx = *(i++);
  ny = *(i++);
  nz = *(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: NX %d   NY %d   NZ  %d \n", nx, ny, nz ENDFB(I->G);
  }

  xlen = *(float*)(i++);
  ylen = *(float*)(i++);
  zlen = *(float*)(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: X %8.3f   Y %8.3f  Z  %8.3f \n", xlen, ylen,
        zlen ENDFB(I->G);
  }

  alpha = *(float*)(i++);
  beta = *(float*)(i++);
  gamma = *(float*)(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: alpha %8.3f   beta %8.3f  gamma %8.3f \n", alpha, beta,
        gamma ENDFB(I->G);
  }

  mapc = *(i++);
  mapr = *(i++);
  maps = *(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: MAPC %d   MAPR %d  MAPS  %d \n", mapc, mapr,
        maps ENDFB(I->G);
  }

  // AMIN, AMAX, AMEAN
  mind = *(float*)(i++);
  maxd = *(float*)(i++);
  mean = *(float*)(i++);

  ispg = *(i++);
  sym_skip = *(i++);

  // CCP4 Format:
  // 25-37 skew transformation
  // 38-52 future use "Storage space sometimes used by specific programs
  // (default = 0)"

  // MRC 2000 Format:
  // 25-49 "Extra, user defined storage space (default = 0)"
  // 50-52 ORIGIN

  {
    int skew = *(i++);

    if (format != cLoadTypeMRC && skew) {
      double matrix[16];
      auto i_float = (const float*)i;

      // SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc)
      copy33f44d(i_float, matrix);
      xx_matrix_invert(matrix, matrix, 4);

      // SKWTRN          Skew translation t
      matrix[3] = i_float[9];
      matrix[7] = i_float[10];
      matrix[11] = i_float[11];

      // Xo(map) = S * (Xo(atoms) - t)
      ObjectStateSetMatrix(ms, matrix);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ObjectMapCCP4: Applied skew transformation\n" ENDFB(I->G);
    }
  }

  // XORIGIN, YORIGIN, ZORIGIN (50 - 52)
  // TODO See "Origin Conventions" in http://situs.biomachina.org/fmap.pdf
  float* mrc2000origin = (float*)(i + 49 - 25);
  if (format != cLoadTypeCCP4Str && lengthsq3f(mrc2000origin) > R_SMALL4) {
    double matrix[] = {
        1., 0., 0., mrc2000origin[0], 0., 1., 0., mrc2000origin[1],
        0., 0., 1., mrc2000origin[2], 0., 0., 0., 1.};

    ObjectStateSetMatrix(ms, matrix);

    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ObjectMapCCP4: Applied MRC 2000 ORIGIN %.2f %.2f %.2f\n",
          mrc2000origin[0], mrc2000origin[1], mrc2000origin[2] ENDFB(I->G);

      if (ncstart || nrstart || nsstart) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapCCP4: Ignoring N*START\n" ENDFB(I->G);
      }
    }

    ncstart = 0;
    nrstart = 0;
    nsstart = 0;
  }

  i += 54 - 25;
  stdev = *(float*)(i++);

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: AMIN %f AMAX %f AMEAN %f ARMS %f\n", mind, maxd, mean,
        stdev ENDFB(I->G);
  }

  n_pts = nc * ns * nr;

  /* at least one EM map encountered lacks NZ, so we'll try to guess it */

  if ((nx == ny) && (!nz)) {
    nz = nx;

    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
      " ObjectMapCCP4: NZ value is zero, but NX = NY, so we'll guess NZ = NX = "
      "NY.\n" ENDFB(I->G);
    }
  }

  expectation = sym_skip + sizeof(int) * 256 + bytes_per_pt * n_pts;

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " ObjectMapCCP4: sym_skip %d bytes %d expectation %d\n", sym_skip, bytes,
        expectation ENDFB(I->G);
  }

  if (bytes < expectation) {
    if (bytes == (expectation - sym_skip)) {
      /* accomodate bogus CCP4 map files with bad symmetry length information */
      if (!quiet) {
        PRINTFB(I->G, FB_ObjectMap, FB_Blather)
        " ObjectMapCCP4: Map has invalid symmetry length -- working "
        "around.\n" ENDFB(I->G);
      }

      expectation -= sym_skip;
      sym_skip = 0;

    } else {
      PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting.\n" ENDFB(I->G);
      return (0);
    }
  }

  q = p + (sizeof(int) * 256) + sym_skip;

  if (little_endian != map_endian && bytes_per_pt > 1) {
    swap_endian(q, n_pts, bytes_per_pt);
  }

  // with normalize == 2, use mean and stdev from file header
  if (normalize == 1 && n_pts > 1) {
    c = n_pts;
    sum = 0.0;
    sumsq = 0.0;
    while (c--) {
      dens = ccp4_next_value(&q, map_mode);
      sumsq += dens * dens;
      sum += dens;
    }
    mean = (float)(sum / n_pts);
    stdev = (float)sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));
    if (stdev < 0.000001) stdev = 1.0;
  }

  q = p + (sizeof(int) * 256) + sym_skip;
  mapc--; /* convert to C indexing... */
  mapr--;
  maps--;

  ms->Div[0] = nx;
  ms->Div[1] = ny;
  ms->Div[2] = nz;

  ms->FDim[mapc] = nc;
  ms->Min[mapc] = ncstart;
  ms->Max[mapc] = nc + ncstart - 1;

  ms->FDim[mapr] = nr;
  ms->Min[mapr] = nrstart;
  ms->Max[mapr] = nr + nrstart - 1;

  ms->FDim[maps] = ns;
  ms->Min[maps] = nsstart;
  ms->Max[maps] = ns + nsstart - 1;

  if (!quiet) {
    if (Feedback(I->G, FB_ObjectMap, FB_Blather)) {
      dump3i(ms->Div, " ObjectMapCCP4: ms->Div");
      dump3i(ms->Min, " ObjectMapCCP4: ms->Min");
      dump3i(ms->Max, " ObjectMapCCP4: ms->Max");
      dump3i(ms->FDim, " ObjectMapCCP4: ms->FDim");
    }
  }

  ms->Symmetry->Crystal.setDims(xlen, ylen, zlen);
  ms->Symmetry->Crystal.setAngles(alpha, beta, gamma);

  if (ispg < n_space_group_numbers) {
    ms->Symmetry->setSpaceGroup(space_group_numbers[ispg]);
  }

  /* -- JV; for vol
  ms->Mean = mean;
  ms->SD = stdev;
  ms->MaxValue = maxd;
  ms->MinValue = mind;
  */

  ms->FDim[3] = 3;
  if (!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
    ok = false;
  else {
    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceCCP4;
    ms->Field->save_points = false;

    maxd = -FLT_MAX;
    mind = FLT_MAX;

    for (cc[maps] = 0; cc[maps] < ms->FDim[maps]; cc[maps]++) {
      v[maps] = (cc[maps] + ms->Min[maps]) / ((float)ms->Div[maps]);

      for (cc[mapr] = 0; cc[mapr] < ms->FDim[mapr]; cc[mapr]++) {
        v[mapr] = (cc[mapr] + ms->Min[mapr]) / ((float)ms->Div[mapr]);

        for (cc[mapc] = 0; cc[mapc] < ms->FDim[mapc]; cc[mapc]++) {
          v[mapc] = (cc[mapc] + ms->Min[mapc]) / ((float)ms->Div[mapc]);

          dens = ccp4_next_value(&q, map_mode);

          if (normalize) dens = (dens - mean) / stdev;
          F3(ms->Field->data, cc[0], cc[1], cc[2]) = dens;
          if (maxd < dens) maxd = dens;
          if (mind > dens) mind = dens;
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          for (e = 0; e < 3; e++)
            F4(ms->Field->points, cc[0], cc[1], cc[2], e) = vr[e];
        }
      }
    }
  }
  if (ok) {
    d = 0;
    for (c = 0; c < 2; c++) {
      v[2] = (c * (ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
      for (b = 0; b < 2; b++) {
        v[1] = (b * (ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
        for (a = 0; a < 2; a++) {
          v[0] = (a * (ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          copy3f(vr, ms->Corner + 3 * d);
          d++;
        }
      }
    }
  }

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " ObjectMapCCP4: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1],
        ms->FDim[2] ENDFB(I->G);
  }

  if (!quiet) {
    if (n_pts > 1) {
      if (normalize) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapCCP4: Normalizing with mean = %8.6f and stdev = %8.6f.\n",
            mean, stdev ENDFB(I->G);
      } else {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " ObjectMapCCP4: Map will not be normalized.\n ObjectMapCCP4: Current "
        "mean = %8.6f and stdev = %8.6f.\n",
            mean, stdev ENDFB(I->G);
      }
    }
  }

  if (ok) {
    v[2] = (ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);
  }
  if (!ok) {
    ErrMessage(I->G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
      " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }

  return (ok);
}
ObjectMapState* getObjectMapState(PyMOLGlobals* G, const char* name,
                                  int state) {
  auto* om = ExecutiveFindObject<ObjectMap>(G, name);
  return om ? om->getObjectMapState(state) : nullptr;
}
/**
 * Export to CCP4 or MRC map format.
 *
 * CCP4:
 * - can have SKWMAT/SKWTRN, but almost no other progam supports them (Maestro
 * does, to some extent)
 * - can't have ORIGIN
 *
 * MRC:
 * - can have ORIGIN, most programs support it
 * - can't have SKWMAT/SKWTRN
 * - can't have N*START (or shouldn't, because of weird ORIGIN conventions)
 * - can't have non-orthogonal axes (or shouldn't, because others might ignore
 * it)
 */
std::vector<char> ObjectMapStateToCCP4Str(const ObjectMapState* ms, int quiet,
                                          int format) {
  std::vector<char> buffer;

  if (!ms || !ms->Active) return buffer;  // empty

  auto G = ms->G;
  auto field = ms->Field->data;

  if (field->type != cFieldFloat || field->base_size != 4) {
    PRINTFB(G, FB_ObjectMap, FB_Errors)
    " MapStateToCCP4-Error: Unsupported field type\n" ENDFB(G);
    return buffer;  // empty
  }

  if (!validateCCP4LoadType(format)) {
    ErrMessage(G, __func__, "wrong format");
    return {};
  }

  buffer.resize(1024 + field->size(), 0);
  auto buffer_s = reinterpret_cast<char*>(&buffer.front());
  auto buffer_i = reinterpret_cast<int32_t*>(&buffer.front());
  auto buffer_f = reinterpret_cast<float*>(&buffer.front());

  buffer_i[0] = ms->FDim[2];  // NC / NX
  buffer_i[1] = ms->FDim[1];  // NR / NY
  buffer_i[2] = ms->FDim[0];  // NS / NZ
  buffer_i[3] = 2;            // MODE (2 = 32-bit reals)
  buffer_i[4] = ms->Min[2];   // NCSTART / NXSTART
  buffer_i[5] = ms->Min[1];   // NRSTART / NYSTART
  buffer_i[6] = ms->Min[0];   // NSSTART / NZSTART
  buffer_i[7] = ms->Div[0];   // NX / MX
  buffer_i[8] = ms->Div[1];   // NY / MY
  buffer_i[9] = ms->Div[2];   // NZ / MZ

  if (ms->Div[0] == 0) {
    // fallback
    buffer_i[7] = ms->FDim[0] - 1;  // NX / MX
    buffer_i[8] = ms->FDim[1] - 1;  // NY / MY
    buffer_i[9] = ms->FDim[2] - 1;  // NZ / MZ
  }

  bool cell_fallback = true;

  // Cell
  if (ms->Symmetry) {
    auto const* dim = ms->Symmetry->Crystal.dims();
    auto const* angle = ms->Symmetry->Crystal.angles();
    buffer_f[10] = dim[0];    // X length / CELL A
    buffer_f[11] = dim[1];    // Y length
    buffer_f[12] = dim[2];    // Z length
    buffer_f[13] = angle[0];  // Alpha    / CELL B
    buffer_f[14] = angle[1];  // Beta
    buffer_f[15] = angle[2];  // Gamma

    // check for 1x1x1 dummy cell
    cell_fallback = fabs(lengthsq3f(dim) - 3.f) < R_SMALL4;
  }

  if (cell_fallback) {
    // fallback: orthogonal extents box
    subtract3f(ms->ExtentMax, ms->ExtentMin, buffer_f + 10);  // CELL A
    set3f(buffer_f + 13, 90.f, 90.f, 90.f);                   // CELL B
  }

  buffer_i[16] = 3;  // MAPC
  buffer_i[17] = 2;  // MAPR
  buffer_i[18] = 1;  // MAPS

  // dummy values
  buffer_f[19] = -5.f;  // AMIN
  buffer_f[20] = 5.f;   // AMAX
  buffer_f[21] = 0.f;   // AMEAN

  // Space group number
  if (ms->Symmetry) {
    for (int ispg = 0; ispg < n_space_group_numbers; ++ispg) {
      if (strcmp(ms->Symmetry->spaceGroup(), space_group_numbers[ispg]) == 0) {
        buffer_i[22] = ispg;  // ISPG
        break;
      }
    }
  }

  buffer_i[23] = 0;  // NSYMBT

  float* b_skwmat = buffer_f + 25;
  float* b_skwtrn = buffer_f + 34;
  float* b_origin = buffer_f + 49;

  // MRC typically only supports orthogonal axes
  const float _f3_90[3]{90.f, 90.f, 90.f};
  bool mrc_possible =
      format != cLoadTypeCCP4Str && equal3f(buffer_f + 13 /* CELL B */, _f3_90);

  if (format == cLoadTypeMRC && !mrc_possible) {
    PRINTFB(G, FB_ObjectMap, FB_Warnings)
    " %s-Warning: MRC expects orthogonal axes\n", __func__ ENDFB(G);
  }

  // skew transformation
  if (!ms->Matrix.empty()) {
    double m[16];
    copy44d(ms->Matrix.data(), m);

    // Skew translation t
    set3f(buffer_f + 34, m[3], m[7], m[11]);  // SKWTRN

    // Skew matrix S (in order S11, S12, S13, S21 etc)
    m[3] = m[7] = m[11] = 0.;
    xx_matrix_invert(m, m, 4);
    copy44d33f(m, buffer_f + 25);  // SKWMAT

    if (mrc_possible) {
      mrc_possible = is_identityf(3, b_skwmat);

      if (format == cLoadTypeMRC && !mrc_possible) {
        PRINTFB(G, FB_ObjectMap, FB_Warnings)
        " %s-Warning: MRC expects orthonormal map\n", __func__ ENDFB(G);
      }
    }

    if (mrc_possible) {
      // translation only -> use origin instead of skew matrix
      copy3(b_skwtrn, b_origin);  // MRC ORIGIN
      zero3(b_skwtrn);
    } else if (format != cLoadTypeMRC) {
      buffer_i[24] = 1;  // LSKFLG
    }
  }

  // origin (stored with skew transformation)
  if (!ms->Origin.empty() && lengthsq3f(ms->Origin.data()) > R_SMALL4) {
    if (!mrc_possible && !buffer_i[24] /* LSKFLG */) {
      identity33f(b_skwmat);
      buffer_i[24] = 1;  // LSKFLG
    }

    if (buffer_i[24] /* LSKFLG */) {
      add3f(ms->Origin.data(), b_skwtrn, b_skwtrn);  // add to SKWTRN
    } else {
      add3f(ms->Origin.data(), b_origin, b_origin);  // add to MRC ORIGIN
    }
  }

  // Don't write ORIGIN and N*START. The convention is to ignore N*START if
  // ORIGIN != (0,0,0), likely because IMOD uses ORIGIN but ignores N*START.
  if (lengthsq3f(b_origin) > R_SMALL4 || format == cLoadTypeMRC) {
    // data origin in fractional coordinates
    float n_start_frac[3] = {
        float(buffer_i[6] /* NSSTART */) / (buffer_i[7] /* NX */),
        float(buffer_i[5] /* NRSTART */) / (buffer_i[8] /* NY */),
        float(buffer_i[4] /* NCSTART */) / (buffer_i[9] /* NZ */)};

    // data origin in Angstrom
    float n_start_real[3];
    transform33f3f(ms->Symmetry->Crystal.fracToReal(), n_start_frac,
                   n_start_real);

    // add to MRC origin
    add3f(n_start_real, b_origin, b_origin);

    // get rid of N*START
    buffer_i[4] = 0;
    buffer_i[5] = 0;
    buffer_i[6] = 0;
  }

  if (format == cLoadTypeCCP4Unspecified) {
    format = mrc_possible ? cLoadTypeMRC : cLoadTypeCCP4Str;
  }

  // Character string 'MAP ' to identify file type
  memcpy(buffer_s + 52 * 4, "MAP ", 4);  // MAP

  // Machine stamp indicating the machine type which wrote file
  if (1 == *(int32_t*)"\1\0\0\0") {
    buffer_i[53] = 0x00004144;  // MACHST (little endian)
  } else {
    buffer_i[53] = 0x11110000;  // MACHST (big endian)
  }

  buffer_f[54] = 1.f;  // ARMS
  buffer_i[55] = 1;    // NLABL
  sprintf(buffer_s + 56 * 4, "PyMOL %s format=%s", _PyMOL_VERSION,
          (format == cLoadTypeMRC) ? "MRC" : "CCP4");  // 57-256 LABEL(20,10)

  // Map data array follows
  memcpy(buffer_s + 1024, field->data.data(), field->size());

  return buffer;
}
static ObjectMap* ObjectMapReadCCP4Str(PyMOLGlobals* G, ObjectMap* I,
                                       char* XPLORStr, int bytes, int state,
                                       int quiet, int format) {
  int ok = true;
  int isNew = true;

  if (!I)
    isNew = true;
  else
    isNew = false;
  if (ok) {
    if (isNew) {
      I = (ObjectMap*)new ObjectMap(G);
      isNew = true;
    } else {
      isNew = false;
    }
    ObjectMapCCP4StrToMap(I, XPLORStr, bytes, state, quiet, format);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}
ObjectMap* ObjectMapLoadCCP4(PyMOLGlobals* G, ObjectMap* obj, const char* fname,
                             int state, int is_string, int bytes, int quiet,
                             int format) {
  ObjectMap* I = nullptr;
  char* buffer;
  long size;

  if (!is_string) {
    if (!quiet) PRINTFB(G, FB_ObjectMap, FB_Actions)
    " ObjectMapLoadCCP4File: Loading from '%s'.\n", fname ENDFB(G);

    buffer = FileGetContents(fname, &size);

    if (!buffer) ErrMessage(G, "ObjectMapLoadCCP4File", "Unable to open file!");
  } else {
    buffer = (char*)fname;
    size = (long)bytes;
  }

  if (buffer) {
    I = ObjectMapReadCCP4Str(G, obj, buffer, size, state, quiet, format);

    if (!is_string) mfree(buffer);

    if (!quiet) {
      if (state < 0) state = I->State.size() - 1;
      if (state < I->State.size()) {
        ObjectMapState* ms;
        ms = &I->State[state];
        if (ms->Active) {
          CrystalDump(&ms->Symmetry->Crystal);
        }
      }
    }
  }
  return (I);
}
