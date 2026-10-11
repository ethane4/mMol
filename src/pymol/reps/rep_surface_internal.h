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

// Helpers shared by the rep_surface*.cc files; not part of the RepSurface API.

#pragma once

#include "pymol/graphics/cgo.h"
#include "pymol/reps/rep_surface.h"
#include "pymol/geometry/sphere.h"

struct RepSurface;
struct SolventDot;
struct SurfaceJobAtomInfo;
struct SurfaceJob;

// Defined in rep_surface_dots.cc.
SolventDot* SolventDotNew(PyMOLGlobals* G, float* coord,
                          std::vector<SurfaceJobAtomInfo>& atom_info,
                          float probe_radius, SphereRec* sp, int* present,
                          int circumscribe, int surface_mode,
                          int surface_solvent, int cavity_cull,
                          int all_visible_flag, float max_vdw, int cavity_mode,
                          float cavity_radius, float cavity_cutoff);

// Defined in rep_surface.cc.
void setShaderCGO(RepSurface* I, CGO* cgo);
void SolventDotFree(SolventDot* I);
bool visibility_test(bool proximityFlag,
                     const int* vi,  // visibility array
                     const int* t);  // indices

// Defined in rep_surface_job.cc.
PyObject* SurfaceJobInputAsTuple(PyMOLGlobals* G, SurfaceJob* I);
PyObject* SurfaceJobResultAsTuple(PyMOLGlobals* G, SurfaceJob* I);
ov_status SurfaceJobResultFromTuple(PyMOLGlobals* G, SurfaceJob* I,
                                    PyObject* tuple);
int SurfaceJobRun(PyMOLGlobals* G, SurfaceJob* I);

enum class SurfaceType {
  Solid = 0,
  DotDefault = 1,
  Triangle = 2,
  SolidDeterministic = 3,
  SolidWeightings = 4,
  SolidScribed = 5,
  SolidEmpirical = 6
};

struct RepSurface : Rep {
  using Rep::Rep;

  ~RepSurface() override;

  cRep_t type() const override { return cRepSurface; }
  void render(RenderInfo* info) override;
  void invalidate(cRepInv_t level) override;
  Rep* recolor() override;
  bool sameVis() const override;
  bool sameColor() const override;

  int N = 0;
  int NT = 0;              //!< number of triangles (?)
  bool proximity = false;  //!< cSetting_surface_proximity
  std::vector<float> V;    //!< Verts
  std::vector<float> VN;   //!< Normals
  std::vector<float> VC;   //!< Colors
  std::vector<float> VA;   //!< Alpha
  std::vector<float> VAO;  //!< Ambient Occlusion per vertex
  std::vector<int> RC;
  std::vector<int> Vis;
  std::vector<int> T;   //!< vertices (of triangles?)
  std::vector<int> S;   //!< strips
  std::vector<int> AT;  //!< closest atom for vertices
  bool oneColorFlag = false;
  int oneColor;
  bool allVisibleFlag = false;
  char* LastVisib = nullptr;
  int* LastColor = nullptr;
  bool ColorInvalidated = false;
  SurfaceType Type;
  float max_vdw;

  /* These variables are for using the shader.  All of them */
  /* are allocated/set when generate_shader_cgo to minimize */
  /* allocation during the rendering loop. */
  CGO* shaderCGO = nullptr;
  CGO* pickingCGO = nullptr;
  bool dot_as_spheres = false;

  int surface_mode = cRepSurface_by_flags;
};

struct SolventDot {
  int nDot{};
  float* dot{};
  float* dotNormal{};
  int* dotCode{};
};

struct SurfaceJobAtomInfo {
  float vdw;
  int flags;
};

struct SurfaceJob {
  /* input */
  std::vector<float> coord;
  std::vector<SurfaceJobAtomInfo> atomInfo{};

  float maxVdw{};
  int allVisibleFlag{};

  int nPresent{};
  std::vector<int> presentVla{};

  int solventSphereIndex{};
  int sphereIndex{};

  SurfaceType surfaceType{};
  int circumscribe{};
  float probeRadius{};
  float carveCutoff{};
  std::vector<float> carveVla{};

  int surfaceMode{};
  int surfaceSolvent{};
  int cavityCull{};
  float pointSep{};
  float trimCutoff{};
  float trimFactor{};

  int cavityMode{};
  float cavityRadius{};
  float cavityCutoff{};

  /* results */
  std::vector<float> V{};
  std::vector<float> VN{};
  int N{};
  std::vector<int> T{};
  std::vector<int> S{};
  int NT{};
};
