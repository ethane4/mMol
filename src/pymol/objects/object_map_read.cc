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

// Reading PHI, XPLOR, FLD, BRIX, GRD, DX and ACNT maps, and ChemPy and NumPy
// arrays.

#include "pymol/base/err.h"
#include "pymol/objects/object_map.h"
#include "pymol/base/parse.h"
#include "pymol/scene/scene.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/file.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"

static int ObjectMapNumPyArrayToMapState(PyMOLGlobals* G, ObjectMapState* I,
                                         PyObject* ary, int quiet);
static int ObjectMapPHIStrToMap(ObjectMap* I, char* PHIStr, int bytes,
                                int state, int quiet) {
  char* p;
  float dens, dens_rev;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 0;
  int map_dim;
  int map_bytes;

  ObjectMapState* ms;

  char cc[MAXLINELEN];
  char* rev;

  little_endian = *((char*)&little_endian);

  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = PHIStr;

  if (*p) /* use FORMATTED IO record to determine map endiness */
    map_endian = 1;
  else
    map_endian = 0;

  p += 4;

  ParseNCopy(cc, p, 20);
  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 20;
  p += 4;

  p += 4;
  ParseNCopy(cc, p, 10);
  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 10;
  ParseNCopy(cc, p, 60);
  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 60;
  p += 4;

  rev = (char*)&dens_rev;

  if (little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] =
        p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
  }

  map_bytes = *((int*)rev);

  map_dim = (int)(pow((map_bytes / 4.0), 1 / 3.0) + 0.5);

  if ((4 * map_dim * map_dim * map_dim) != map_bytes) /* consistency check */
    map_dim = 65;

  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " PHIStrToMap: Map Size %d x %d x %d\n", map_dim, map_dim,
        map_dim ENDFB(I->G);
  }
  p += 4;

  ms->FDim[0] = map_dim;
  ms->FDim[1] = map_dim;
  ms->FDim[2] = map_dim;
  ms->FDim[3] = 3;

  ms->Div[0] = (map_dim - 1) / 2;
  ms->Div[1] = (map_dim - 1) / 2;
  ms->Div[2] = (map_dim - 1) / 2;
  ms->Min[0] = -ms->Div[0];
  ms->Min[1] = -ms->Div[1];
  ms->Min[2] = -ms->Div[2];
  ms->Max[0] = ms->Div[0];
  ms->Max[1] = ms->Div[1];
  ms->Max[2] = ms->Div[2];

  ms->Field.reset(new Isofield(I->G, ms->FDim));
  ms->MapSource = cMapSourceGeneralPurpose;
  ms->Field->save_points = false;

  for (c = 0; c < ms->FDim[2];
       c++) { /* z y x ordering into c b a  so that x = a, etc. */
    for (b = 0; b < ms->FDim[1]; b++) {
      for (a = 0; a < ms->FDim[0]; a++) {
        if (little_endian != map_endian) {
          rev[0] = p[3];
          rev[1] = p[2];
          rev[2] = p[1];
          rev[3] = p[0];
        } else {
          rev[0] = p[0]; /* gotta go char by char because of memory alignment
                            issues ... */
          rev[1] = p[1];
          rev[2] = p[2];
          rev[3] = p[3];
        }
        dens = *((float*)rev);
        F3(ms->Field->data, a, b, c) = dens;
        if (maxd < dens) maxd = dens;
        if (mind > dens) mind = dens;
        p += 4;
      }
    }
  }
  p += 4;

  p += 4;
  ParseNCopy(cc, p, 16);
  if (!quiet) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " PHIStrToMap: %s\n", cc ENDFB(I->G);
  }
  p += 16;
  p += 4;

  ms->Grid = std::vector<float>(3);
  p += 4;
  if (little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] =
        p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
  }
  ms->Grid[0] = 1.0F / (*((float*)rev));
  ms->Grid[1] = ms->Grid[0];
  ms->Grid[2] = ms->Grid[0];
  p += 4;

  ms->Origin = std::vector<float>(3);
  if (little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] =
        p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
    ;
  }
  ms->Origin[0] = *((float*)rev);
  p += 4;

  if (little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] =
        p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
    ;
  }
  ms->Origin[1] = *((float*)rev);
  p += 4;
  if (little_endian != map_endian) {
    rev[0] = p[3];
    rev[1] = p[2];
    rev[2] = p[1];
    rev[3] = p[0];
  } else {
    rev[0] =
        p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1] = p[1];
    rev[2] = p[2];
    rev[3] = p[3];
    ;
  }
  ms->Origin[2] = *((float*)rev);
  p += 4;

  p += 4;

  if (ok) {
    for (e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }
  }

  for (c = 0; c < ms->FDim[2]; c++) {
    v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
    for (b = 0; b < ms->FDim[1]; b++) {
      v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
      for (a = 0; a < ms->FDim[0]; a++) {
        v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
        for (e = 0; e < 3; e++) {
          F4(ms->Field->points, a, b, c, e) = v[e];
        }
      }
    }
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

  /* interpolation test code
     {
     float test[3];
     float result;
     float cmp;

     for(c=0;c<ms->FDim[0]-1;c++) {
     for(b=0;b<ms->FDim[1]-1;b++) {
     for(a=0;a<ms->FDim[2]-1;a++) {
     for(e=0;e<3;e++) {
     test[e] = (F4(ms->Field->points,a,b,c,e)+
     F4(ms->Field->points,a,b,c,e))/2.0;
     }
     ObjectMapStateInterpolate(ms,test,&result,1);
     cmp = (F3(ms->Field->data,a,b,c)+
     F3(ms->Field->data,a,b,c))/2.0;
     if(fabs(cmp-result)>0.001) {
     printf("%d %d %d\n",a,b,c);
     printf("%8.5f %8.5f\n",
     cmp,
     result);
     }
     }
     }
     }
     }
   */

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
static int ObjectMapXPLORStrToMap(ObjectMap* I, char* XPLORStr, int state,
                                  int quiet) {
  char* p;
  int a, b, c, d, e;
  float v[3], vr[3], dens, maxd, mind;
  char cc[MAXLINELEN];
  int n;
  int ok = true;
  ObjectMapState* ms;

  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }

  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = XPLORStr;

  while (*p) {
    p = ParseNCopy(cc, p, 8);
    if (!*cc)
      p = ParseNextLine(p);
    else if (sscanf(cc, "%i", &n) == 1) {
      p = ParseWordCopy(cc, p, MAXLINELEN);
      if (strstr(cc, "!NTITLE") || (!*cc)) {
        p = ParseNextLine(p);
        while (n--) {
          p = ParseNextLine(p);
        }
      } else if (strstr(cc, "REMARKS")) {
        p = ParseNextLine(p);
      } else {
        break;
      }
    }
  }
  if (*p) { /* n contains first dimension */
    ms->Div[0] = n;
    if (sscanf(cc, "%i", &ms->Min[0]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Max[0]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Div[1]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Min[1]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Max[1]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Div[2]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Min[2]) != 1) ok = false;
    p = ParseNCopy(cc, p, 8);
    if (sscanf(cc, "%i", &ms->Max[2]) != 1) ok = false;
    p = ParseNextLine(p);
    p = ParseNCopy(cc, p, 12);

    float cellparams[3];

    if (sscanf(cc, "%f", &cellparams[0]) != 1) ok = false;
    p = ParseNCopy(cc, p, 12);
    if (sscanf(cc, "%f", &cellparams[1]) != 1) ok = false;
    p = ParseNCopy(cc, p, 12);
    if (sscanf(cc, "%f", &cellparams[2]) != 1) ok = false;

    ms->Symmetry->Crystal.setDims(cellparams);

    p = ParseNCopy(cc, p, 12);
    if (sscanf(cc, "%f", &cellparams[0]) != 1) ok = false;
    p = ParseNCopy(cc, p, 12);
    if (sscanf(cc, "%f", &cellparams[1]) != 1) ok = false;
    p = ParseNCopy(cc, p, 12);
    if (sscanf(cc, "%f", &cellparams[2]) != 1) ok = false;

    ms->Symmetry->Crystal.setAngles(cellparams);

    p = ParseNextLine(p);
    p = ParseNCopy(cc, p, 3);
    if (strcmp(cc, "ZYX")) ok = false;
    p = ParseNextLine(p);

  } else {
    ok = false;
  }
  if (ok) {
    ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
    ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
    ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
    ms->FDim[3] = 3;
    if (!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(I->G, ms->FDim));
      ms->MapSource = cMapSourceCrystallographic;
      ms->Field->save_points = false;
      for (c = 0; c < ms->FDim[2]; c++) {
        v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
        p = ParseNextLine(p);
        for (b = 0; b < ms->FDim[1]; b++) {
          v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
          for (a = 0; a < ms->FDim[0]; a++) {
            v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
            p = ParseNCopy(cc, p, 12);
            if (!cc[0]) {
              p = ParseNextLine(p);
              p = ParseNCopy(cc, p, 12);
            }
            if (sscanf(cc, "%f", &dens) != 1) {
              ok = false;
            } else {
              F3(ms->Field->data, a, b, c) = dens;
              if (maxd < dens) maxd = dens;
              if (mind > dens) mind = dens;
            }
            transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
            for (e = 0; e < 3; e++) {
              F4(ms->Field->points, a, b, c, e) = vr[e];
            }
          }
        }
        p = ParseNextLine(p);
      }
      if (ok) {
        d = 0;
        for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
          for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
            for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
              transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
              copy3f(vr, ms->Corner + 3 * d);
              d++;
            }
          }
        }
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
      " ObjectMap: Map read.  Range = %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  }

  return (ok);
}
static int ObjectMapFLDStrToMap(ObjectMap* I, char* PHIStr, int bytes,
                                int state, int quiet) {
  char* p;
  float dens, dens_rev;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 1;
  char cc[MAXLINELEN];
  char* rev;
  int ndim = 0;
  int veclen = 0;
  int nspace = 0;
  ObjectMapState* ms;
  int got_ndim = false;
  int got_dim1 = false;
  int got_dim2 = false;
  int got_dim3 = false;
  int got_data = false;
  int got_veclen = false;
  int got_min_ext = false;
  int got_max_ext = false;
  int got_field = false;
  int got_nspace = false;

  little_endian = *((char*)&little_endian);
  map_endian = little_endian;

  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = PHIStr;

  while ((*p) && (!(got_ndim && got_dim1 && got_dim2 && got_dim3 &&
                    got_veclen && got_min_ext && got_max_ext && got_field &&
                    got_nspace && got_data))) {
    if (!got_ndim) {
      ParseWordCopy(cc, p, 4);
      if (!strcmp(cc, "ndim")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ndim) == 1) got_ndim = true;
      }
    }

    if (!got_dim1) {
      ParseWordCopy(cc, p, 4);
      if (!strcmp(cc, "dim1")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ms->FDim[0]) == 1) got_dim1 = true;
      }
    }

    if (!got_dim2) {
      ParseWordCopy(cc, p, 4);
      if (!strcmp(cc, "dim2")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ms->FDim[1]) == 1) got_dim2 = true;
      }
    }

    if (!got_dim3) {
      ParseWordCopy(cc, p, 4);
      if (!strcmp(cc, "dim3")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ms->FDim[2]) == 1) got_dim3 = true;
      }
    }

    if (!got_nspace) {
      ParseWordCopy(cc, p, 6);
      if (!strcmp(cc, "nspace")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &nspace) == 1) got_nspace = true;
      }
    }

    if (!got_veclen) {
      ParseWordCopy(cc, p, 6);
      if (!strcmp(cc, "veclen")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &veclen) == 1) got_veclen = true;
      }
    }

    if (!got_data) {
      ParseWordCopy(cc, p, 4);
      if (!strcmp(cc, "data")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (!strcmp(cc, "float")) got_data = true;
      }
    }

    if (!got_min_ext) {
      ParseWordCopy(cc, p, 7);
      if (!strcmp(cc, "min_ext")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%f", &ms->ExtentMin[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if (sscanf(cc, "%f", &ms->ExtentMin[1]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%f", &ms->ExtentMin[2]) == 1) {
              got_min_ext = true;
            }
          }
        }
      }
    }

    if (!got_max_ext) {
      ParseWordCopy(cc, p, 7);
      if (!strcmp(cc, "max_ext")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%f", &ms->ExtentMax[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if (sscanf(cc, "%f", &ms->ExtentMax[1]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%f", &ms->ExtentMax[2]) == 1) {
              got_max_ext = true;
            }
          }
        }
      }
    }

    if (!got_field) {
      ParseWordCopy(cc, p, 5);
      if (!strcmp(cc, "field")) {
        p = ParseSkipEquals(p);
        p = ParseWordCopy(cc, p, 50);
        if (!strcmp(cc, "uniform")) got_field = true;
      }
    }
    p = ParseNextLine(p);
  }

  if (got_ndim && got_dim1 && got_dim2 && got_dim3 && got_veclen &&
      got_min_ext && got_max_ext && got_field && got_nspace && got_data) {
    int pass = 0;

    ms->Origin = std::vector<float>(3);
    ms->Range = std::vector<float>(3);
    ms->Grid = std::vector<float>(3);

    copy3f(ms->ExtentMin, ms->Origin.data());
    subtract3f(ms->ExtentMax, ms->ExtentMin, ms->Range.data());
    ms->FDim[3] = 3;

    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " FLDStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1],
        ms->FDim[2] ENDFB(I->G);

    for (a = 0; a < 3; a++) {
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
      ms->Div[a] = ms->FDim[a] - 1;

      if (ms->FDim[a])
        ms->Grid[a] = ms->Range[a] / (ms->Max[a]);
      else
        ms->Grid[a] = 0.0F;
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceFLD;
    ms->Field->save_points = false;

    while (1) {
      maxd = -FLT_MAX;
      mind = FLT_MAX;
      p = PHIStr;

      while (*p) { /* ^L^L sentinel */
        if ((p[0] == 12) && (p[1] == 12)) {
          p += 2;
          break;
        }
        p++;
      }

      rev = (char*)&dens_rev;
      for (c = 0; c < ms->FDim[2];
           c++) { /* z y x ordering into c b a  so that x = a, etc. */
        for (b = 0; b < ms->FDim[1]; b++) {
          for (a = 0; a < ms->FDim[0]; a++) {
            if (little_endian != map_endian) {
              rev[0] = p[3];
              rev[1] = p[2];
              rev[2] = p[1];
              rev[3] = p[0];
            } else {
              rev[0] = p[0]; /* gotta go char by char because of memory
                                alignment issues ... */
              rev[1] = p[1];
              rev[2] = p[2];
              rev[3] = p[3];
            }
            dens = *((float*)rev);
            F3(ms->Field->data, a, b, c) = dens;
            if (maxd < dens) maxd = dens;
            if (mind > dens) mind = dens;
            p += 4;
          }
        }
      }

      /* There's no way to determine the original handedness of input
         field files.  So instead, we simplymake an educated guess about
         whether we're byte-swapped based on the range of the density
         values obtained. */

      if (((maxd / FLT_MAX) > 0.1F) && ((mind / (-FLT_MAX)) > 0.1F)) {
        if (pass == 0) {
          map_endian = (!map_endian); /* okay, try again swapped */
        } else if (pass == 1) {
          /* didn't help, so resort to original order */
          map_endian = (!map_endian);
        } else {
          break;
        }
      } else {
        break;
      }
      pass++;
    }

    for (c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for (b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for (a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for (e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
    }

    d = 0;
    for (c = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      if (!c)
        v[2] = ms->ExtentMin[2];
      else
        v[2] = ms->ExtentMax[2];

      for (b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);

        if (!b)
          v[1] = ms->ExtentMin[1];
        else
          v[1] = ms->ExtentMax[1];

        for (a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);

          if (!a)
            v[0] = ms->ExtentMin[0];
          else
            v[0] = ms->ExtentMax[0];

          copy3f(v, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
      " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
    " Error: unable to read FLD file.\n" ENDFB(I->G);
    /*  printf(" got_ndim %d\n",got_ndim);
       printf(" got_dim1 %d\n",got_dim1);
       printf(" got_dim2 %d\n",got_dim2);
       printf(" got_dim3 %d\n",got_dim3);
       printf(" got_veclen %d\n",got_veclen);
       printf(" got_min_ext %d\n",got_min_ext);
       printf(" got_max_ext %d\n",got_max_ext);
       printf(" got_field %d\n",got_field);
       printf(" got_nspace %d\n",got_nspace);
       printf(" got_data %d\n",got_data);
     */
  }

  return (ok);
}
static int ObjectMapBRIXStrToMap(ObjectMap* I, char* BRIXStr, int bytes,
                                 int state, int quiet) {
  char *p, *pp;
  float dens;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  char cc[MAXLINELEN];
  ObjectMapState* ms;
  int got_origin = false;
  int got_extent = false;
  int got_grid = false;
  int got_cell = false;
  int got_prod = false;
  int got_plus = false;
  int got_sigma = false;
  int got_smiley = false;
  float prod = 1.0F;
  float plus = 0.0F;
  float sigma = 0.0F;
  float calc_sigma = 0.0F;
  float calc_mean = 0.0F;
  int normalize;
  int swap_bytes;

  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_o_maps);
  swap_bytes = SettingGetGlobal_b(I->G, cSetting_swap_dsn6_bytes);
  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];

  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = BRIXStr;

  ParseNCopy(cc, p, 3);
  got_smiley = (strcmp(cc, ":-)") == 0);

  if (got_smiley) {
    /* ASCII BRIX format header */

    while ((*p) && (!(got_origin && got_extent && got_grid && got_cell &&
                      got_prod && got_plus && got_sigma))) {
      if (!got_origin) {
        pp = ParseWordCopy(cc, p, 6);
        if (WordMatch(I->G, "origin", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%d", &ms->Min[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%d", &ms->Min[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if (sscanf(cc, "%d", &ms->Min[2]) == 1) {
                got_origin = true;
              }
            }
          }
        }
      }

      if (!got_extent) {
        pp = ParseWordCopy(cc, p, 6);
        if (WordMatch(I->G, "extent", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%d", &ms->Max[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%d", &ms->Max[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if (sscanf(cc, "%d", &ms->Max[2]) == 1) {
                got_extent = true;
                ms->Max[0] += ms->Min[0] - 1;
                ms->Max[1] += ms->Min[1] - 1;
                ms->Max[2] += ms->Min[2] - 1;
              }
            }
          }
        }
      }

      if (!got_grid) {
        pp = ParseWordCopy(cc, p, 4);
        if (WordMatch(I->G, "grid", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%d", &ms->Div[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%d", &ms->Div[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if (sscanf(cc, "%d", &ms->Div[2]) == 1) {
                got_grid = true;
              }
            }
          }
        }
      }

      if (!got_cell) {
        pp = ParseWordCopy(cc, p, 4);
        if (WordMatch(I->G, "cell", cc, true) < 0) {
          p = pp;
          got_cell = true;
          float cellparams[6];

          for (int i = 0; i != 6; ++i) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%f", cellparams + i) != 1) {
              got_cell = false;
              break;
            }
          }

          if (got_cell) {
            ms->Symmetry->Crystal.setDims(cellparams);
            ms->Symmetry->Crystal.setAngles(cellparams + 3);
          }
        }
      }

      if (!got_plus) {
        pp = ParseWordCopy(cc, p, 4);
        if (WordMatch(I->G, "plus", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%f", &plus) == 1) {
            got_plus = true;
          }
        }
      }

      if (!got_prod) {
        pp = ParseWordCopy(cc, p, 4);
        if (WordMatch(I->G, "prod", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%f", &prod) == 1) {
            got_prod = true;
          }
        }
      }

      if (!got_sigma) {
        pp = ParseWordCopy(cc, p, 5);
        if (WordMatch(I->G, "sigma", cc, true) < 0) {
          p = ParseWordCopy(cc, pp, 50);
          if (sscanf(cc, "%f", &sigma) == 1) {
            got_sigma = true;
          }
        }
      }

      p++;
    }
  } else {
    /* Binary DSN6 format */

    /* first, determine whether or not we need to byte-swap the header... */

    int passed_endian_check = false;
    short int* shint_ptr;
    float scale_factor;
    char tmp_char;
    int a;
    int start_swap_at = 0;
    int stop_swap_at = bytes;

    shint_ptr = (short int*)(p + 18 * 2);
    if (*shint_ptr == 100) {
      passed_endian_check = true;
      start_swap_at = 512; /* don't byte-swap header */
    }

    if (!swap_bytes) {
      stop_swap_at = 512;
    }
    /* for DSN6, always swap the bricks */

    p = BRIXStr + start_swap_at;
    for (a = start_swap_at; a < stop_swap_at; a += 2) {
      tmp_char = *p;
      *p = *(p + 1);
      *(p + 1) = tmp_char;
      p += 2;
    }

    p = BRIXStr;

    if (*shint_ptr == 100) {
      passed_endian_check = true;
    }

    if (!passed_endian_check) {
      PRINTFB(I->G, FB_ObjectMap, FB_Errors)
      " Error: This looks like a DSN6 map file, but I can't match "
      "endianness.\n" ENDFB(I->G);
    } else {
      shint_ptr = (short int*)p;

      ms->Min[0] = *(shint_ptr++);
      ms->Min[1] = *(shint_ptr++);
      ms->Min[2] = *(shint_ptr++);
      got_origin = true;

      ms->Max[0] = *(shint_ptr++);
      ms->Max[1] = *(shint_ptr++);
      ms->Max[2] = *(shint_ptr++);
      got_extent = true;
      ms->Max[0] += ms->Min[0] - 1;
      ms->Max[1] += ms->Min[1] - 1;
      ms->Max[2] += ms->Min[2] - 1;

      ms->Div[0] = *(shint_ptr++);
      ms->Div[1] = *(shint_ptr++);
      ms->Div[2] = *(shint_ptr++);
      got_grid = true;

      float cellparams[6] = {
          float(*(shint_ptr++)), float(*(shint_ptr++)), float(*(shint_ptr++)),
          float(*(shint_ptr++)), float(*(shint_ptr++)), float(*(shint_ptr++)),
      };
      got_cell = true;

      prod = (float)(*(shint_ptr++)) / 100.0F;
      got_prod = true;

      plus = (float)(*(shint_ptr++));
      got_plus = true;

      scale_factor = (float)(*(shint_ptr++));

      for (a = 0; a < 6; ++a) {
        cellparams[a] /= scale_factor;
      }
      ms->Symmetry->Crystal.setDims(cellparams);
      ms->Symmetry->Crystal.setAngles(cellparams + 3);
    }
  }

  if (got_origin && got_extent && got_grid && got_cell && got_plus &&
      got_prod) {
    PRINTFB(I->G, FB_ObjectMap, FB_Blather)
    " BRIXStrToMap: Prod = %8.3f, Plus = %8.3f\n", prod, plus ENDFB(I->G);

    ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
    ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
    ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
    ms->FDim[3] = 3;
    if (!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(I->G, ms->FDim));
      ms->MapSource = cMapSourceBRIX;
      ms->Field->save_points = false;

      {
        int block_size = 8;
        int xblocks = ((ms->FDim[0] - 1) / block_size) + 1;
        int yblocks = ((ms->FDim[1] - 1) / block_size) + 1;
        int zblocks = ((ms->FDim[2] - 1) / block_size) + 1;
        int cc, bb, aa, ia, ib, ic, xa, xb, xc;
        double sum = 0.0;
        double sumsq = 0.0;
        int n_pts = 0;

        p = BRIXStr + 512;
        for (cc = 0; cc < zblocks; cc++) {
          for (bb = 0; bb < yblocks; bb++) {
            for (aa = 0; aa < xblocks; aa++) {
              for (c = 0; c < block_size; c++) {
                xc = c + cc * block_size;
                ic = xc + ms->Min[2];
                v[2] = ic / ((float)ms->Div[2]);

                for (b = 0; b < block_size; b++) {
                  xb = b + bb * block_size;
                  ib = xb + ms->Min[1];
                  v[1] = ib / ((float)ms->Div[1]);

                  for (a = 0; a < block_size; a++) {
                    xa = a + aa * block_size;
                    ia = xa + ms->Min[0];
                    v[0] = ia / ((float)ms->Div[0]);

                    dens = (((float)(*((unsigned char*)(p++)))) - plus) / prod;
                    if ((ia <= ms->Max[0]) && (ib <= ms->Max[1]) &&
                        (ic <= ms->Max[2])) {
                      F3(ms->Field->data, xa, xb, xc) = dens;
                      sumsq += dens * dens;
                      sum += dens;
                      n_pts++;
                      if (maxd < dens) maxd = dens;
                      if (mind > dens) mind = dens;
                      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                      for (e = 0; e < 3; e++) {
                        F4(ms->Field->points, xa, xb, xc, e) = vr[e];
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if (n_pts > 1) {
          calc_mean = (float)(sum / n_pts);
          calc_sigma =
              (float)sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));
        }
      }

      if (ok) {
        d = 0;
        for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
          v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
          for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
            v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
            for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
              v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
              transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
              copy3f(vr, ms->Corner + 3 * d);
              d++;
            }
          }
        }
      }

      if (ok && normalize) {
        if (calc_sigma > R_SMALL8) {
          for (c = 0; c < ms->FDim[2]; c++) {
            for (b = 0; b < ms->FDim[1]; b++) {
              for (a = 0; a < ms->FDim[0]; a++) {
                F3(ms->Field->data, a, b, c) =
                    (F3(ms->Field->data, a, b, c) - calc_mean) / calc_sigma;
              }
            }
          }
        }
      }

      v[2] = (ms->Min[2]) / ((float)ms->Div[2]);
      v[1] = (ms->Min[1]) / ((float)ms->Div[1]);
      v[0] = (ms->Min[0]) / ((float)ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

      v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
      v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
      v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);

      transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " BRIXStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1],
          ms->FDim[2] ENDFB(I->G);

      if (got_sigma) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " BRIXStrToMap: Reported Sigma = %8.3f\n", sigma ENDFB(I->G);
      }

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " BRIXStrToMap: Range = %5.6f to %5.6f\n", mind, maxd ENDFB(I->G);

      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " BRIXStrToMap: Calculated Mean = %8.3f, Sigma = %8.3f\n", calc_mean,
          calc_sigma ENDFB(I->G);

      if (normalize) {
        PRINTFB(I->G, FB_ObjectMap, FB_Details)
        " BRIXStrToMap: Normalizing...\n" ENDFB(I->G);
      }

      ms->Active = true;
      ObjectMapUpdateExtents(I);

      if (!quiet) {
        CrystalDump(&ms->Symmetry->Crystal);
      }
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
    " Error: unable to read BRIX/DSN6 file.\n" ENDFB(I->G);
  }

  return (ok);
}
static int ObjectMapGRDStrToMap(ObjectMap* I, char* GRDStr, int bytes,
                                int state, int quiet) {
  /* NOTE: binary GRD reader not yet validated */

  char* p;
  float dens;
  float* f = nullptr;
  int a, b, c, d, e;
  float v[3], vr[3], maxd, mind;
  int ok = true;
  char cc[MAXLINELEN];
  ObjectMapState* ms;
  int got_cell = false;
  int got_brick = false;
  int fast_axis = 1; /* assumes fast X */

  int got_ranges = false;
  int normalize = false;
  float mean, stdev;
  double sum = 0.0;
  double sumsq = 0.0;
  int n_pts = 0;
  int ascii = true;
  int little_endian = 1, map_endian = 0;
  union int_or_char_array {
    int block_len;
    char rev[4];
  };
  union int_or_char_array rev_union;
  rev_union.block_len = 0;

  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }
  ms = &I->State[state];
  normalize = SettingGetGlobal_b(I->G, cSetting_normalize_grd_maps);
  maxd = -FLT_MAX;
  mind = FLT_MAX;

  p = GRDStr;

  if (bytes > 4) {
    if ((!p[0]) || (!p[1]) || (!p[2]) || (!p[3])) {
      /* if zero appears in the first four bytes, then this is a binary map */
      ascii = false;

      little_endian = *((char*)&little_endian);
      map_endian = (*p || *(p + 1));
    }
  }
  /* print map title */

  if (ascii) {
    p = ParseNCopy(cc, p, 100);
  } else {
    /*

    according to one site on the internet...GRD binary formats look like this:

            write(14) title
            write(14) scale,oldmid
            write(14) ivary,nbyte,intdat,extent,extent
            write(14) extent,xang,yang,zang,xstart
            write(14) xend,ystart,yend,zstart,zend
            write(14) intx,inty,intz

            do 100 i = 1,intz+1
                do 100 j = 1,inty+1
    100            write(14)(phimap(k,j,i),k=1,intx+1)

        where: scale is the inverse grid spacing,
        oldmid are x,y,z coordinates of the grid center,
        intx,inty,intz = grid points per side - 1
        ivary = 0
        nbyte = 4
        intdat = 0
        xang,yang,zang = 90
        extent is the absolute value of the x,y, or z coordinate of the grid
        corners with the largest absolute value
        xstart,xend,ystart, yend,zstart,zend are the fractional limits (-1 to 1)
        of the molecule in each direction.
        title is a descriptive header for the molecule

    However, that doesn't make sense...in the files I see, scale and oldmid
    don't seem to exist!

    So here's another claim which I find more credible...

    character*132 toplbl !ascii header
    integer*4 ivary !0 => x index varys most rapidly
    integer*4 nbyte !=4, # of bytes in data
    integer*4 inddat !=0, floating point data
    real*4 xang,yang,zang !=90,90,90 unit cell angles
    integer*4 intx,inty,intz !=igrid-1, # of intervals/grid side
    real*4 extent !maximum extent of grid
    real*4 xstart,xend !beginning, end of grid sides
    real*4 ystart,yend !in fractional
    real*4 zstart,zend !units of extent
    write(14)toplbl write(14)ivary, nbyte, intdat, extent, extent, extent, xang,
    yang, zang, xstart, xend, ystart, yend, zstart, zend, intx, inty, intz do k
    = 1,igrid do j = 1,igrid write(14)(phimap(i,j,k),i=1,igrid) end do end d

    */

    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
      " ObjectMapGRD-Warning: Binary GRD reader not yet validated.\n" ENDFB(
          I->G);
    }

    if (little_endian != map_endian) {
      rev_union.rev[0] = p[3];
      rev_union.rev[1] = p[2];
      rev_union.rev[2] = p[1];
      rev_union.rev[3] = p[0];
    } else {
      rev_union.rev[0] = p[0]; /* gotta go char by char because of memory
                                  alignment issues ... */
      rev_union.rev[1] = p[1];
      rev_union.rev[2] = p[2];
      rev_union.rev[3] = p[3];
    }
    ParseNCopy(cc, p + 4, rev_union.block_len);
    /* now flip file to correct endianess */
    if (little_endian != map_endian) {
      char* c = p;
      int cnt = bytes >> 2;
      unsigned char c0, c1, c2, c3;
      while (cnt) {
        c0 = c[0];
        c1 = c[1];
        c2 = c[2];
        c3 = c[3];
        c[0] = c3;
        c[1] = c2;
        c[2] = c1;
        c[3] = c0;
        c += 4;
        cnt--;
      }
    }
    p += 4 + rev_union.block_len;
    f = (float*)p;
  }

  PRINTFB(I->G, FB_ObjectMap, FB_Details)
  " ObjectMap: %s\n", cc ENDFB(I->G);

  if (ascii) p = ParseNextLine(p);

  /* skip format -- we're reading float regardless... */
  if (ascii)
    p = ParseNextLine(p);
  else {
    {
      int block_len_check;
      int ivary;
      int nbyte;
      int intdat;

      block_len_check = *((int*)(f++));

      if (rev_union.block_len != block_len_check) {
        PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapGRD-Warning: block length not matched -- not a true GRD "
        "binary?\n" ENDFB(I->G);
      }

      rev_union.block_len = *((int*)(f++));
      ivary = *((int*)(f++));
      nbyte = *((int*)(f++));
      intdat = *((int*)(f++));

      if ((ivary) || (nbyte != 4) || (intdat)) {
        if (!quiet) {
          PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
          " ObjectMapGRD-Warning: funky ivary, nbyte, intdat -- not a true GRD "
          "binary?\n" ENDFB(I->G);
        }
      }
    }
  }

  /* read unit cell */

  if (ok) {
    float cellparams[6];

    if (ascii) {
      p = ParseWordCopy(cc, p, 50);
      if (sscanf(cc, "%f", &cellparams[0]) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%f", &cellparams[1]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if (sscanf(cc, "%f", &cellparams[2]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%f", &cellparams[3]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if (sscanf(cc, "%f", &cellparams[4]) == 1) {
                p = ParseWordCopy(cc, p, 50);
                if (sscanf(cc, "%f", &cellparams[5]) == 1) {
                  got_cell = true;
                }
              }
            }
          }
        }
      }
      ok = got_cell;
    } else {
      cellparams[0] = *(f++); /* x-extent */
      cellparams[1] = *(f++); /* y-extent */
      cellparams[2] = *(f++); /* z-extent */
      cellparams[3] = *(f++); /* xang */
      cellparams[4] = *(f++); /* yang */
      cellparams[5] = *(f++); /* zang */
      got_cell = 1;
    }

    ms->Symmetry->Crystal.setDims(cellparams);
    ms->Symmetry->Crystal.setAngles(cellparams + 3);
  }

  if (ascii) p = ParseNextLine(p);

  if (ascii) {
    /* read brick dimensions */

    if (ok) {
      p = ParseWordCopy(cc, p, 50);
      if (sscanf(cc, "%d", &ms->FDim[0]) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ms->FDim[1]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if (sscanf(cc, "%d", &ms->FDim[2]) == 1) {
            got_brick = true;
            ms->FDim[0]++;
            ms->FDim[1]++;
            ms->FDim[2]++; /* stupid 0-based counters */
          }
        }
      }
      ok = got_brick;
    }

    p = ParseNextLine(p);

    /* read ranges */

    if (ok) {
      p = ParseWordCopy(cc, p, 50);
      if (sscanf(cc, "%d", &fast_axis) == 1) {
        p = ParseWordCopy(cc, p, 50);
        if (sscanf(cc, "%d", &ms->Min[0]) == 1) {
          p = ParseWordCopy(cc, p, 50);
          if (sscanf(cc, "%d", &ms->Div[0]) == 1) {
            p = ParseWordCopy(cc, p, 50);
            if (sscanf(cc, "%d", &ms->Min[1]) == 1) {
              p = ParseWordCopy(cc, p, 50);
              if (sscanf(cc, "%d", &ms->Div[1]) == 1) {
                p = ParseWordCopy(cc, p, 50);
                if (sscanf(cc, "%d", &ms->Min[2]) == 1) {
                  p = ParseWordCopy(cc, p, 50);
                  if (sscanf(cc, "%d", &ms->Div[2]) == 1) {
                    got_ranges = true;
                  }
                }
              }
            }
          }
        }
      }
      ok = got_ranges;
    }
  } else {
    float fmin[3];
    float fmax[3];

    fmin[0] = *(f++);
    fmax[0] = *(f++);
    fmin[1] = *(f++);
    fmax[1] = *(f++);
    fmin[2] = *(f++);
    fmax[2] = *(f++);

    dump3f(fmin, "fmin");
    dump3f(fmax, "fmax");
    ms->FDim[0] = *((int*)f++) + 1;
    ms->FDim[1] = *((int*)f++) + 1;
    ms->FDim[2] = *((int*)f++) + 1;

    ms->Div[0] = pymol_roundf((ms->FDim[0] - 1) / (fmax[0] - fmin[0]));
    ms->Div[1] = pymol_roundf((ms->FDim[1] - 1) / (fmax[1] - fmin[1]));
    ms->Div[2] = pymol_roundf((ms->FDim[2] - 1) / (fmax[2] - fmin[2]));

    ms->Min[0] = pymol_roundf(ms->Div[0] * fmin[0]);
    ms->Min[1] = pymol_roundf(ms->Div[1] * fmin[1]);
    ms->Min[2] = pymol_roundf(ms->Div[2] * fmin[2]);

    ms->Max[0] = ms->Min[0] + ms->FDim[0] - 1;
    ms->Max[1] = ms->Min[1] + ms->FDim[1] - 1;
    ms->Max[2] = ms->Min[2] + ms->FDim[2] - 1;
    got_ranges = true;

    {
      int block_len_check;

      block_len_check = *((int*)(f++));
      if (rev_union.block_len != block_len_check) {
        PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapGRD-Warning: block length not matched -- not a true GRD "
        "binary?\n" ENDFB(I->G);
      }
    }
  }

  if (ok) {
    if (ascii) {
      ms->Div[0] = ms->Div[0] - ms->Min[0];
      ms->Div[1] = ms->Div[1] - ms->Min[1];
      ms->Div[2] = ms->Div[2] - ms->Min[2];

      ms->Max[0] = ms->Min[0] + ms->FDim[0] - 1;
      ms->Max[1] = ms->Min[1] + ms->FDim[1] - 1;
      ms->Max[2] = ms->Min[2] + ms->FDim[2] - 1;
    }

    ms->FDim[3] = 3;

    if (Feedback(I->G, FB_ObjectMap, FB_Blather)) {
      dump3i(ms->Div, "ms->Div");
      dump3i(ms->Min, "ms->Min");
      dump3i(ms->Max, "ms->Max");
      dump3i(ms->FDim, "ms->FDim");
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceGRD;
    ms->Field->save_points = false;

    switch (fast_axis) {
      case 3: /* Fast Y - BROKEN! */
        PRINTFB(I->G, FB_ObjectMap, FB_Warnings)
        " ObjectMapGRD-Warning: fast_axis %d unsupported!\n",
            fast_axis ENDFB(I->G);
        /* intentional fall though... */
      case 1: /* Fast X */
      default:
        for (c = 0; c < ms->FDim[2]; c++) {
          v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
          for (b = 0; b < ms->FDim[1]; b++) {
            if (!ascii) f++; /* skip block delimiter */
            v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
            for (a = 0; a < ms->FDim[0]; a++) {
              v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
              if (ascii) {
                p = ParseNextLine(p);
                p = ParseNCopy(cc, p, 24);
                if (sscanf(cc, "%f", &dens) != 1) {
                  ok = false;
                }
              } else {
                dens = *(f++);
              }
              if (ok) {
                sumsq += dens * dens;
                sum += dens;
                n_pts++;
                F3(ms->Field->data, a, b, c) = dens;
                if (maxd < dens) maxd = dens;
                if (mind > dens) mind = dens;
              }
              transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
              for (e = 0; e < 3; e++) {
                F4(ms->Field->points, a, b, c, e) = vr[e];
              }
            }
            if (!ascii) f++; /* skip fortran block delimiter */
          }
        }
        break;
    }
  }

  if (n_pts > 1) {
    mean = (float)(sum / n_pts);
    stdev = (float)sqrt1d((sumsq - (sum * sum / n_pts)) / (n_pts - 1));

    if (normalize) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ObjectMapGRDStrToMap: Normalizing: mean = %8.6f and stdev = %8.6f.\n",
          mean, stdev ENDFB(I->G);
      if (stdev < R_SMALL8) stdev = 1.0F;
      for (c = 0; c < ms->FDim[2]; c++)
        for (b = 0; b < ms->FDim[1]; b++) {
          for (a = 0; a < ms->FDim[0]; a++) {
            dens = F3(ms->Field->data, a, b, c);
            F3(ms->Field->data, a, b, c) = (dens - mean) / stdev;
          }
        }
    } else {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ObjectMapGRDStrToMap: Mean = %8.6f and stdev = %8.6f.\n", mean,
          stdev ENDFB(I->G);
    }
  }

  if (ok) {
    d = 0;
    for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
      v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
      for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
        v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
        for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
          v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
          transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
          copy3f(vr, ms->Corner + 3 * d);
          d++;
        }
      }
    }

    v[2] = (ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = (ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = (ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMin);

    v[2] = ((ms->FDim[2] - 1) + ms->Min[2]) / ((float)ms->Div[2]);
    v[1] = ((ms->FDim[1] - 1) + ms->Min[1]) / ((float)ms->Div[1]);
    v[0] = ((ms->FDim[0] - 1) + ms->Min[0]) / ((float)ms->Div[0]);

    transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, ms->ExtentMax);

    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " GRDXStrToMap: Map Size %d x %d x %d\n", ms->FDim[0], ms->FDim[1],
        ms->FDim[2] ENDFB(I->G);

    ms->Active = true;
    ObjectMapUpdateExtents(I);
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectMap, FB_Results)
      " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(I->G);

      CrystalDump(&ms->Symmetry->Crystal);
    }
  } else {
    PRINTFB(I->G, FB_ObjectMap, FB_Errors)
    " Error: unable to read GRD file.\n" ENDFB(I->G);
  }

  return (ok);
}
static ObjectMap* ObjectMapReadXPLORStr(PyMOLGlobals* G, ObjectMap* I,
                                        char* XPLORStr, int state, int quiet) {
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
    ObjectMapXPLORStrToMap(I, XPLORStr, state, quiet);

    SceneChanged(I->G);
    SceneCountFrames(I->G);
  }
  return (I);
}
static ObjectMap* ObjectMapReadFLDStr(PyMOLGlobals* G, ObjectMap* I,
                                      char* MapStr, int bytes, int state,
                                      int quiet) {
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
    ObjectMapFLDStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}
static ObjectMap* ObjectMapReadBRIXStr(PyMOLGlobals* G, ObjectMap* I,
                                       char* MapStr, int bytes, int state,
                                       int quiet) {
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
    ObjectMapBRIXStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}
static ObjectMap* ObjectMapReadGRDStr(PyMOLGlobals* G, ObjectMap* I,
                                      char* MapStr, int bytes, int state,
                                      int quiet) {
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
    ObjectMapGRDStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}
static ObjectMap* ObjectMapReadPHIStr(PyMOLGlobals* G, ObjectMap* I,
                                      char* MapStr, int bytes, int state,
                                      int quiet) {
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
    ObjectMapPHIStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}
ObjectMap* ObjectMapLoadPHI(PyMOLGlobals* G, ObjectMap* obj, const char* fname,
                            int state, int is_string, int bytes, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  if (!is_string) {
    if (!quiet) PRINTFB(G, FB_ObjectMap, FB_Actions)
    " ObjectMapLoadPHIFile: Loading from '%s'.\n", fname ENDFB(G);

    buffer = FileGetContents(fname, &size);

    if (!buffer) ErrMessage(G, "ObjectMapLoadPHIFile", "Unable to open file!");
  } else {
    buffer = (char*)fname;
    size = (long)bytes;
  }

  if (buffer) {
    I = ObjectMapReadPHIStr(G, obj, buffer, size, state, quiet);

    if (!is_string) mfree(buffer);
  }
  return (I);
}
static int is_number(char* p) {
  int result = (*p != 0);
  char c;
  if (result)
    while ((c = *(p++))) {
      if (!((c == '.') || (c == '-') || (c == '+') || (c == 'e') ||
            (c == 'E') || ((c >= '0') && (c <= '9'))))
        result = false;
      break;
    }
  return result;
}

static pymol::Result<std::unique_ptr<ObjectMapState>> ObjectMapDXStrToMap(
    PyMOLGlobals* G, const char* DXStr, int bytes, bool quiet) {
  if (DXStr[0] == '\x1f' && DXStr[1] == '\x8b') {
    return pymol::make_error("gzipped data not supported");
  }

  enum {
    ASCII_NUMERIC,
    BINARY_DOUBLE,
    BINARY_FLOAT,
  } data_type = ASCII_NUMERIC;

  int n_items = 0;

  const char *p, *pp;
  float dens;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  /* DX named from their docs */

  int stage = 0;

  char cc[MAXLINELEN];

  auto ms = std::make_unique<ObjectMapState>(G);

  ms->Origin = std::vector<float>(3);
  ms->Grid = std::vector<float>(3);

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = DXStr;

  /* get the dimensions */

  ms->FDim[3] = 3;

  while ((*p) && stage == 0) {
    pp = p;
    p = ParseNCopy(cc, p, 35);
    if ((strcmp(cc, "object 1 class gridpositions counts") == 0) ||
        is_number(cc)) {
      if (is_number(cc)) p = pp;
      p = ParseWordCopy(cc, p, 10);
      if (sscanf(cc, "%d", &ms->FDim[0]) == 1) {
        p = ParseWordCopy(cc, p, 10);
        if (sscanf(cc, "%d", &ms->FDim[1]) == 1) {
          p = ParseWordCopy(cc, p, 10);
          if (sscanf(cc, "%d", &ms->FDim[2]) == 1) {
            stage = 1;
          }
        }
      }
    }
    p = ParseNextLine(p);
  }

  if (stage != 1) {
    return pymol::make_error(
        "missing 'object 1 class gridpositions count' line");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
    " DXStrToMap: Dimensions: %d %d %d\n", ms->FDim[0], ms->FDim[1],
        ms->FDim[2] ENDFB(G);
  }

  /* get the origin */

  while ((*p) && stage == 1) {
    pp = p;
    p = ParseNCopy(cc, p, 6);
    if ((strcmp(cc, "origin") == 0) || is_number(cc)) {
      if (is_number(cc)) p = pp;
      p = ParseWordCopy(cc, p, 20);
      if (sscanf(cc, "%f", &ms->Origin[0]) == 1) {
        p = ParseWordCopy(cc, p, 20);
        if (sscanf(cc, "%f", &ms->Origin[1]) == 1) {
          p = ParseWordCopy(cc, p, 20);
          if (sscanf(cc, "%f", &ms->Origin[2]) == 1) {
            stage = 2;
          }
        }
      }
    }
    p = ParseNextLine(p);
  }

  if (stage != 2) {
    return pymol::make_error("missing 'origin'");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
    " DXStrToMap: Origin %8.3f %8.3f %8.3f\n", ms->Origin[0], ms->Origin[1],
        ms->Origin[2] ENDFB(G);
  }

  float delta[9];
  int delta_i = 0;

  while ((*p) && stage == 2) {
    pp = p;
    p = ParseNCopy(cc, p, 5);

    if (strcmp(cc, "delta") != 0) {
      if (is_number(cc)) {
        p = pp;
      } else {
        p = ParseNextLine(p);
        continue;
      }
    }

    if (3 != sscanf(p, " %f %f %f", delta + delta_i, delta + delta_i + 1,
                    delta + delta_i + 2)) {
      return pymol::make_error("expected 3 floats");
    }

    p = ParseNextLine(p);
    delta_i += 3;

    if (delta_i == 9) {
      stage = 3;

      if (is_diagonalf(3, delta)) {
        ms->Grid[0] = delta[0];
        ms->Grid[1] = delta[4];
        ms->Grid[2] = delta[8];
      } else {
        if (ms->Matrix.empty()) ms->Matrix = std::vector<double>(16);

        copy33f44d(delta, ms->Matrix.data());
        ms->Matrix[3] = ms->Origin[0];
        ms->Matrix[7] = ms->Origin[1];
        ms->Matrix[11] = ms->Origin[2];

        ones3f(ms->Grid.data());
        zero3f(ms->Origin.data());
      }
    }
  }

  if (stage != 3) {
    return pymol::make_error("missing 'delta'");
  }

  {
    PRINTFB(G, FB_ObjectMap, FB_Details)
    " DXStrToMap: Grid %8.3f %8.3f %8.3f\n", ms->Grid[0], ms->Grid[1],
        ms->Grid[2] ENDFB(G);
  }

  while ((*p) && stage == 3) {
    p = ParseNCopy(cc, p, 6);
    if (strcmp(cc, "object") == 0) {
      p = ParseWordCopy(cc, p, 20);
      char cc_type[8] = "";
      char cc_binary[8] = "";
      if (2 <= sscanf(p, " class array type %7s rank %*s items %s %7s", cc_type,
                      cc, cc_binary)) {
        if (sscanf(cc, "%d", &n_items) == 1) {
          if (n_items == ms->FDim[0] * ms->FDim[1] * ms->FDim[2]) stage = 4;
        }

        if (strcmp(cc_binary, "binary") == 0) {
          if (strcmp(cc_type, "double") == 0) {
            data_type = BINARY_DOUBLE;
          } else if (strcmp(cc_type, "float") == 0) {
            data_type = BINARY_FLOAT;
          } else {
            return pymol::make_error(
                pymol::string_format("type '%s' not supported", cc_type));
          }
        }
      }
    } else if (is_number(cc)) {
      n_items = ms->FDim[0] * ms->FDim[1] * ms->FDim[2];
      stage = 4;
      break;
    }
    p = ParseNextLine(p);
  }

  if (stage != 4) {
    return pymol::make_error("missing 'object . class array'");
  }

  PRINTFB(G, FB_ObjectMap, FB_Details)
  " DXStrToMap: %d data points.\n", n_items ENDFB(G);

  unsigned malformed_ascii_floats_count = 0;

  {
    ms->Field.reset(new Isofield(G, ms->FDim));
    ms->MapSource = cMapSourceGeneralPurpose;
    ms->Field->save_points = false;

    for (a = 0; a < 3; a++) {
      ms->Div[a] = ms->FDim[a] - 1;
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
    }

    for (a = 0; a < ms->FDim[0]; a++) {
      for (b = 0; b < ms->FDim[1]; b++) {
        for (c = 0; c < ms->FDim[2]; c++) {
          switch (data_type) {
            case BINARY_DOUBLE:
              dens = *reinterpret_cast<double const*>(p);
              p += sizeof(double);
              break;
            case BINARY_FLOAT:
              dens = *reinterpret_cast<float const*>(p);
              p += sizeof(float);
              break;
            case ASCII_NUMERIC:
              p = ParseWordCopy(cc, p, 20);
              if (!cc[0]) {
                p = ParseNextLine(p);
                p = ParseWordCopy(cc, p, 20);
              }
              if (sscanf(cc, "%f", &dens) != 1) {
                ++malformed_ascii_floats_count;
                continue;
              }
              break;
            default:
              assert(false);
          }

          if (maxd < dens) maxd = dens;
          if (mind > dens) mind = dens;
          F3(ms->Field->data, a, b, c) = dens;
        }
      }
    }

    for (e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }

    for (c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for (b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for (a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for (e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
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

    if (!quiet) {
      PRINTFB(G, FB_ObjectMap, FB_Results)
      " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(G);
    }
  }

  if (malformed_ascii_floats_count > 0) {
    PRINTFB(G, FB_ObjectMap, FB_Warnings)
    " %s-Warning: Failed to parse %u values\n", __func__,
        malformed_ascii_floats_count ENDFB(G);
  }

  ms->Active = true;
  return std::move(ms);
}
ObjectMap* ObjectMapReadDXStr(PyMOLGlobals* G, ObjectMap* I, const char* MapStr,
                              int bytes, int state, bool quiet) {
  auto mapstateresult = ObjectMapDXStrToMap(G, MapStr, bytes, quiet);
  if (!mapstateresult) {
    ErrMessage(G, __func__, mapstateresult.error().what().c_str());
    return nullptr;
  }

  assert(mapstateresult.result());
  assert(mapstateresult.result()->Active);

  if (!I) {
    I = new ObjectMap(G);
  }

  if (state < 0) {
    state = I->State.size();
  }

  VecCheckEmplace(I->State, state, G);

  I->State[state] = std::move(*mapstateresult.result());

  ObjectMapUpdateExtents(I);
  SceneChanged(G);
  SceneCountFrames(G);

  return (I);
}
ObjectMap* ObjectMapLoadDXFile(PyMOLGlobals* G, ObjectMap* obj,
                               const char* fname, int state, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  buffer = FileGetContents(fname, &size);

  if (!buffer) {
    ErrMessage(G, "ObjectMapLoadDXFile", "Unable to open file!");
    PRINTFB(G, FB_ObjectMap, FB_Errors)
    "ObjectMapLoadDXFile: Does '%s' exist?\n", fname ENDFB(G);
  } else {
    if (Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadDXFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadDXStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);
}
static int ObjectMapACNTStrToMap(ObjectMap* I, char* ACNTStr, int bytes,
                                 int state, int quiet) {
  int n_items = 0;

  char* p;
  float dens;
  int a, b, c, d, e;
  float v[3], maxd, mind;
  int ok = true;

  int stage = 0;

  ObjectMapState* ms;

  char cc[MAXLINELEN];

  if (state < 0) state = I->State.size();
  if (I->State.size() <= state) {
    VecCheckEmplace(I->State, state, I->G);
  }

  ms = &I->State[state];

  ms->Origin = std::vector<float>(3);
  ms->Grid = std::vector<float>(3);

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  p = ACNTStr;

  /* skip header */

  p = ParseNextLine(p);

  /* get the origin, spacing, and dimensions */

  ms->FDim[3] = 3;

  /* read the map info */

  {
    int i;
    for (i = 0; i < 3; i++) {
      int ii = (4 - i) % 3; /* y, x , z */
      p = ParseWordCopy(cc, p, MAXLINELEN);
      if (sscanf(cc, "%f", &ms->Origin[ii]) == 1) {
        p = ParseWordCopy(cc, p, MAXLINELEN);
        if (sscanf(cc, "%f", &ms->Grid[ii]) == 1) {
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if (sscanf(cc, "%d", &ms->FDim[ii]) == 1) {
            p = ParseNextLine(p);
            stage++;
          }
        }
      }
    }
  }

  /* skip the angles (ignored -- should probably check instead */

  p = ParseNextLine(p);

  /* read the data block */

  if (ok && (stage == 3)) {
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " ACNTStrToMap: Dimensions: %d %d %d\n", ms->FDim[0], ms->FDim[1],
        ms->FDim[2] ENDFB(I->G);
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " ACNTStrToMap: Origin %8.3f %8.3f %8.3f\n", ms->Origin[0], ms->Origin[1],
        ms->Origin[2] ENDFB(I->G);
    PRINTFB(I->G, FB_ObjectMap, FB_Details)
    " ACNTStrToMap: Grid %8.3f %8.3f %8.3f\n", ms->Grid[0], ms->Grid[1],
        ms->Grid[2] ENDFB(I->G);

    n_items = ms->FDim[0] * ms->FDim[1] * ms->FDim[2];

    if (ok && (stage == 1)) {
      PRINTFB(I->G, FB_ObjectMap, FB_Details)
      " ACNTStrToMap: %d data points.\n", n_items ENDFB(I->G);
    }

    ms->Field.reset(new Isofield(I->G, ms->FDim));
    ms->MapSource = cMapSourceGeneralPurpose;
    ms->Field->save_points = false;

    for (a = 0; a < 3; a++) {
      ms->Div[a] = ms->FDim[a] - 1;
      ms->Min[a] = 0;
      ms->Max[a] = ms->FDim[a] - 1;
    }

    for (c = 0; c < ms->FDim[2]; c++) {
      for (a = 0; a < ms->FDim[0]; a++) {
        for (b = 0; b < ms->FDim[1]; b++) {
          p = ParseWordCopy(cc, p, MAXLINELEN);
          p = ParseNextLine(p);
          if (sscanf(cc, "%f", &dens) == 1) {
            if (maxd < dens) maxd = dens;
            if (mind > dens) mind = dens;
            F3(ms->Field->data, a, b, c) = dens;
          } else {
            ok = false;
          }
        }
      }
    }

    for (e = 0; e < 3; e++) {
      ms->ExtentMin[e] = ms->Origin[e] + ms->Grid[e] * ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e] + ms->Grid[e] * ms->Max[e];
    }

    for (c = 0; c < ms->FDim[2]; c++) {
      v[2] = ms->Origin[2] + ms->Grid[2] * (c + ms->Min[2]);
      for (b = 0; b < ms->FDim[1]; b++) {
        v[1] = ms->Origin[1] + ms->Grid[1] * (b + ms->Min[1]);
        for (a = 0; a < ms->FDim[0]; a++) {
          v[0] = ms->Origin[0] + ms->Grid[0] * (a + ms->Min[0]);
          for (e = 0; e < 3; e++) {
            F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
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

    if (ok) stage = 4;
  }

  if (stage != 4) ok = false;

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

static ObjectMap* ObjectMapReadACNTStr(PyMOLGlobals* G, ObjectMap* I,
                                       char* MapStr, int bytes, int state,
                                       int quiet) {
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
    ObjectMapACNTStrToMap(I, MapStr, bytes, state, quiet);
    SceneChanged(G);
    SceneCountFrames(G);
  }
  return (I);
}

ObjectMap* ObjectMapLoadACNTFile(PyMOLGlobals* G, ObjectMap* obj,
                                 const char* fname, int state, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  buffer = FileGetContents(fname, &size);

  if (!buffer) {
    ErrMessage(G, "ObjectMapLoadACNTFile", "Unable to open file!");
    PRINTFB(G, FB_ObjectMap, FB_Errors)
    "ObjectMapLoadACNTFile: Does '%s' exist?\n", fname ENDFB(G);
  } else {
    if (Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadACNTFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadACNTStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);
}
ObjectMap* ObjectMapLoadFLDFile(PyMOLGlobals* G, ObjectMap* obj,
                                const char* fname, int state, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  buffer = FileGetContents(fname, &size);

  if (!buffer) {
    ErrMessage(G, "ObjectMapLoadFLDFile", "Unable to open file!");
  } else {
    if (Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadFLDFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadFLDStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);
}
ObjectMap* ObjectMapLoadBRIXFile(PyMOLGlobals* G, ObjectMap* obj,
                                 const char* fname, int state, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  buffer = FileGetContents(fname, &size);

  if (!buffer) {
    ErrMessage(G, "ObjectMapLoadBRIXFile", "Unable to open file!");
  } else {
    if (Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadBRIXFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadBRIXStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);
}
ObjectMap* ObjectMapLoadGRDFile(PyMOLGlobals* G, ObjectMap* obj,
                                const char* fname, int state, int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  buffer = FileGetContents(fname, &size);

  if (!buffer) {
    ErrMessage(G, "ObjectMapLoadGRDFile", "Unable to open file!");
  } else {
    if (Feedback(G, FB_ObjectMap, FB_Actions)) {
      printf(" ObjectMapLoadGRDFile: Loading from '%s'.\n", fname);
    }

    I = ObjectMapReadGRDStr(G, obj, buffer, size, state, quiet);

    mfree(buffer);
  }
  return (I);
}
ObjectMap* ObjectMapLoadXPLOR(PyMOLGlobals* G, ObjectMap* obj,
                              const char* fname, int state, int is_file,
                              int quiet) {
  ObjectMap* I = nullptr;
  long size;
  char* buffer;

  if (is_file) {
    buffer = FileGetContents(fname, &size);

    if (!buffer) ErrMessage(G, "ObjectMapLoadXPLOR", "Unable to open file!");
  } else {
    buffer = (char*)fname;
  }

  if (buffer) {
    if ((!quiet) && (Feedback(G, FB_ObjectMap, FB_Actions))) {
      if (is_file) {
        printf(" ObjectMapLoadXPLOR: Loading from '%s'.\n", fname);
      } else {
        printf(" ObjectMapLoadXPLOR: Loading...\n");
      }
    }

    I = ObjectMapReadXPLORStr(G, obj, buffer, state, quiet);

    if (is_file) mfree(buffer);
    if (!quiet) {
      if (Feedback(G, FB_ObjectMap, FB_Details)) {
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
  }
  return (I);
}
static int ObjectMapNumPyArrayToMapState(PyMOLGlobals* G, ObjectMapState* ms,
                                         PyObject* ary, int quiet) {
  int a, b, c, d, e;
  float v[3], dens, maxd, mind;
  int ok = true;
  void* ptr;

  maxd = -FLT_MAX;
  mind = FLT_MAX;
  if (ok) {
    ms->FDim[0] = ms->Dim[0];
    ms->FDim[1] = ms->Dim[1];
    ms->FDim[2] = ms->Dim[2];
    ms->FDim[3] = 3;

    if (!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
      ok = false;
    else {
      ms->Field.reset(new Isofield(G, ms->FDim));
      for (c = 0; c < ms->FDim[2]; c++) {
        v[2] = ms->Origin[2] + ms->Grid[2] * c;
        for (b = 0; b < ms->FDim[1]; b++) {
          v[1] = ms->Origin[1] + ms->Grid[1] * b;
          for (a = 0; a < ms->FDim[0]; a++) {
            v[0] = ms->Origin[0] + ms->Grid[0] * a;
            dens = 0.0;
            F3(ms->Field->data, a, b, c) = dens;
            if (maxd < dens) maxd = dens;
            if (mind > dens) mind = dens;
            for (e = 0; e < 3; e++) F4(ms->Field->points, a, b, c, e) = v[e];
          }
        }
      }
      d = 0;
      for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
        v[2] = ms->Origin[2] + ms->Grid[2] * c;
        for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
          v[1] = ms->Origin[1] + ms->Grid[1] * b;
          for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
            v[0] = ms->Origin[0] + ms->Grid[0] * a;
            copy3f(v, ms->Corner + 3 * d);
            d++;
          }
        }
      }
    }
  }
  if (ok) {
    copy3f(ms->Origin.data(), ms->ExtentMin);
    copy3f(ms->Origin.data(), ms->ExtentMax);
    add3f(ms->Range.data(), ms->ExtentMax, ms->ExtentMax);
  }
  if (!ok) {
    ErrMessage(G, "ObjectMap", "Error reading map");
  } else {
    ms->Active = true;
    if (!quiet) {
      PRINTFB(G, FB_ObjectMap, FB_Results)
      " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind, maxd ENDFB(G);
    }
  }
  return (ok);
}
ObjectMap* ObjectMapLoadChemPyBrick(PyMOLGlobals* G, ObjectMap* I,
                                    PyObject* Map, int state, int discrete,
                                    int quiet) {
  ErrMessage(G, "ObjectMap",
             "missing numpy support, required for chempy.brick");
  return nullptr;

  int ok = true;
  int isNew = true;
  PyObject* tmp;
  ObjectMapState* ms;

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

    if (state < 0) state = I->State.size();
    if (I->State.size() <= state) {
      VecCheckEmplace(I->State, state, I->G);
    }
    ms = &I->State[state];

    if (PyObject_HasAttrString(Map, "origin") &&
        PyObject_HasAttrString(Map, "dim") &&
        PyObject_HasAttrString(Map, "range") &&
        PyObject_HasAttrString(Map, "grid") &&
        PyObject_HasAttrString(Map, "lvl")) {
      tmp = PyObject_GetAttrString(Map, "origin");
      if (tmp) {
        PConvFromPyObject(G, tmp, ms->Origin);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick origin.");
      tmp = PyObject_GetAttrString(Map, "dim");
      if (tmp) {
        PConvFromPyObject(G, tmp, ms->Dim);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick dimension.");
      tmp = PyObject_GetAttrString(Map, "range");
      if (tmp) {
        PConvFromPyObject(G, tmp, ms->Range);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick range.");
      tmp = PyObject_GetAttrString(Map, "grid");
      if (tmp) {
        PConvFromPyObject(G, tmp, ms->Grid);
        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick grid.");
      tmp = PyObject_GetAttrString(Map, "lvl");
      if (tmp) {
        ObjectMapNumPyArrayToMapState(G, ms, tmp, quiet);

        Py_DECREF(tmp);
      } else
        ok = ErrMessage(G, "ObjectMap", "missing brick density.");

    } else
      ok = ErrMessage(G, "ObjectMap", "missing any brick attribute.");

    SceneChanged(G);
    SceneCountFrames(G);
    if (ok) {
      int a;
      for (a = 0; a < 3; a++) {
        ms->Min[a] = 0;
        ms->Max[a] = ms->Dim[a] - 1;
      }
      ms->Active = true;
      ms->MapSource = cMapSourceChempyBrick;
      ObjectMapUpdateExtents(I);
    }
  }
  return (I);
}
ObjectMap* ObjectMapLoadChemPyMap(PyMOLGlobals* G, ObjectMap* I, PyObject* Map,
                                  int state, int discrete, int quiet) {
  int ok = true;
  int isNew = true;
  float* cobj;
  WordType format;
  float v[3], vr[3], dens, maxd, mind;
  int a, b, c, d, e;
  ObjectMapState* ms;

  maxd = -FLT_MAX;
  mind = FLT_MAX;

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

    if (state < 0) state = I->State.size();
    if (I->State.size() <= state) {
      VecCheckEmplace(I->State, state, I->G);
    }
    ms = &I->State[state];

    float dim[3], angle[3];

    if (!PConvAttrToStrMaxLen(Map, "format", format, sizeof(WordType) - 1))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'format' parameter.");
    else if (!PConvAttrToFloatArrayInPlace(Map, "cell_dim", dim, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_dim' parameter.");
    else if (!PConvAttrToFloatArrayInPlace(Map, "cell_ang", angle, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_ang' parameter.");
    else if (!PConvAttrToIntArrayInPlace(Map, "cell_div", ms->Div, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'cell_div' parameter.");
    else if (!PConvAttrToIntArrayInPlace(Map, "first", ms->Min, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'first' parameter.");
    else if (!PConvAttrToIntArrayInPlace(Map, "last", ms->Max, 3))
      ok = ErrMessage(G, "LoadChemPyMap", "bad 'last' parameter.");

    ms->Symmetry->Crystal.setDims(dim);
    ms->Symmetry->Crystal.setAngles(angle);

    if (ok) {
      if (strcmp(format, "CObjectZYXfloat") == 0) {
        ok = PConvAttrToPtr(Map, "c_object", (void**)(void*)&cobj);
        if (!ok) ErrMessage(G, "LoadChemPyMap", "pymol::CObject unreadable.");
      } else {
        ok = ErrMessage(G, "LoadChemPyMap", "unsupported format.");
      }
    }
    /* good to go */

    if (ok) {
      if (strcmp(format, "CObjectZYXfloat") == 0) {
        ms->FDim[0] = ms->Max[0] - ms->Min[0] + 1;
        ms->FDim[1] = ms->Max[1] - ms->Min[1] + 1;
        ms->FDim[2] = ms->Max[2] - ms->Min[2] + 1;
        if (Feedback(G, FB_ObjectMap, FB_Actions)) {
          printf(" LoadChemPyMap: CObjectZYXdouble %dx%dx%d\n", ms->FDim[0],
                 ms->FDim[1], ms->FDim[2]);
        }
        ms->FDim[3] = 3;
        if (!(ms->FDim[0] && ms->FDim[1] && ms->FDim[2]))
          ok = false;
        else {
          ms->Field.reset(new Isofield(G, ms->FDim));
          for (c = 0; c < ms->FDim[2]; c++) {
            v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
            for (b = 0; b < ms->FDim[1]; b++) {
              v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
              for (a = 0; a < ms->FDim[0]; a++) {
                v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);

                dens = *(cobj++);

                F3(ms->Field->data, a, b, c) = dens;
                if (maxd < dens) maxd = dens;
                if (mind > dens) mind = dens;
                transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                for (e = 0; e < 3; e++)
                  F4(ms->Field->points, a, b, c, e) = vr[e];
              }
            }
          }

          if (ok) {
            d = 0;
            for (c = 0; c < ms->FDim[2]; c += (ms->FDim[2] - 1)) {
              v[2] = (c + ms->Min[2]) / ((float)ms->Div[2]);
              for (b = 0; b < ms->FDim[1]; b += (ms->FDim[1] - 1)) {
                v[1] = (b + ms->Min[1]) / ((float)ms->Div[1]);
                for (a = 0; a < ms->FDim[0]; a += (ms->FDim[0] - 1)) {
                  v[0] = (a + ms->Min[0]) / ((float)ms->Div[0]);
                  transform33f3f(ms->Symmetry->Crystal.fracToReal(), v, vr);
                  copy3f(vr, ms->Corner + 3 * d);
                  d++;
                }
              }
            }
          }
        }
      }
    }

    if (ok) {
      CrystalDump(&ms->Symmetry->Crystal);

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
      ErrMessage(G, "ObjectMap", "Error reading map");
    } else {
      ms->Active = true;
      ObjectMapUpdateExtents(I);
      if (!quiet) {
        PRINTFB(I->G, FB_ObjectMap, FB_Results)
        " ObjectMap: Map read.  Range: %5.3f to %5.3f\n", mind,
            maxd ENDFB(I->G);
      }
    }

    if (ok) {
      SceneChanged(G);
      SceneCountFrames(G);
    }
  }
  return (I);
}
