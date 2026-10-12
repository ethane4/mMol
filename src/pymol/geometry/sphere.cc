

/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC.
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

#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"
#include "pymol/platform/os_gl.h"

#include "pymol/base/base.h"
#include "pymol/geometry/sphere.h"
#include "pymol/geometry/vector.h"
#include "pymol/base/err.h"

#include "pymol/platform/memory_debug.h"

static const int mesh[30][2] = {
    {0, 3},  {0, 4},  {0, 5}, {0, 8}, {0, 11}, {1, 2},  {1, 6},  {1, 7},
    {1, 8},  {1, 11}, {2, 6}, {2, 7}, {2, 9},  {2, 10}, {3, 4},  {3, 5},
    {3, 9},  {3, 10}, {4, 7}, {4, 8}, {4, 9},  {5, 6},  {5, 10}, {5, 11},
    {6, 10}, {6, 11}, {7, 8}, {7, 9}, {8, 11}, {9, 10}};

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

/**
 * The precomputed icosphere meshes (levels 0-4), from
 * $PYMOL_DATA/sphere_meshes.bin: each table as an int32 count then its values,
 * in the order below. Loaded once and shared by every PyMOL instance.
 */
namespace {
struct SphereTables {
  std::vector<int32_t> NTri, NStrip, NVertTot, nDot, dot_start;
  std::vector<float> dot, area;
  std::vector<int32_t> StripLen_start, StripLen, Sequence_start, Sequence;
  std::vector<int32_t> Tri_start, Tri;

  template <typename T>
  static bool read(std::ifstream& in, std::vector<T>& v) {
    int32_t n = 0;
    if (!in.read(reinterpret_cast<char*>(&n), sizeof(n)) || n < 0) return false;
    v.resize(n);
    return bool(in.read(reinterpret_cast<char*>(v.data()), n * sizeof(T)));
  }

  bool load() {
    const char* data_dir = getenv("PYMOL_DATA");
    if (!data_dir) return false;
    std::ifstream in(std::string(data_dir) + "/sphere_meshes.bin",
                     std::ios::binary);
    return read(in, NTri) && read(in, NStrip) && read(in, NVertTot) &&
           read(in, nDot) && read(in, dot_start) && read(in, dot) &&
           read(in, area) && read(in, StripLen_start) && read(in, StripLen) &&
           read(in, Sequence_start) && read(in, Sequence) &&
           read(in, Tri_start) && read(in, Tri) &&
           NTri.size() == NUMBER_OF_SPHERE_LEVELS;
  }
};

const SphereTables& sphereTables() {
  static SphereTables tables;
  static bool loaded = tables.load();
  if (!loaded) {
    fprintf(stderr, "PyMOL: can't read $PYMOL_DATA/sphere_meshes.bin\n");
    abort();
  }
  return tables;
}
}  // namespace

void SphereInit(PyMOLGlobals* G) {
  CSphere* I = (G->Sphere = pymol::calloc<CSphere>(1));

  auto& T = sphereTables();
  I->Array = pymol::malloc<SphereRec>(NUMBER_OF_SPHERE_LEVELS);

  for (int i = 0; i < NUMBER_OF_SPHERE_LEVELS; i++) {
    // The records point into the shared tables, which are never modified.
    auto& rec = I->Array[i];
    rec.area = const_cast<float*>(&T.area[T.dot_start[i]]);
    rec.dot = reinterpret_cast<Vector3f*>(
        const_cast<float*>(&T.dot[3 * T.dot_start[i]]));
    rec.StripLen = const_cast<int*>(&T.StripLen[T.StripLen_start[i]]);
    rec.Sequence = const_cast<int*>(&T.Sequence[T.Sequence_start[i]]);
    rec.NStrip = T.NStrip[i];
    rec.NVertTot = T.NVertTot[i];
    rec.nDot = T.nDot[i];
    rec.Tri = const_cast<int*>(&T.Tri[T.Tri_start[i]]);
    rec.NTri = T.NTri[i];

    if (i) {
      rec.Mesh = nullptr;
      rec.NMesh = 0;
    } else {
      rec.Mesh = (int*)(void*)mesh;
      rec.NMesh = 30;
    }
    I->Sphere[i] = &rec;
  }
}

void SphereFree(PyMOLGlobals* G) {
  CSphere* I = G->Sphere;

  FreeP(I->Array);
  FreeP(I);
}

SphereRec* GetSpheroidSphereRec(PyMOLGlobals* G) {
  return G->Sphere->Sphere[2];
}
