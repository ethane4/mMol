
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
-* Cameron Mura
-*
-*
Z* -------------------------------------------------------------------
*/

// Helpers shared by the rep_cartoon*.cc files; not part of the RepCartoon API.

#pragma once

#include "pymol/graphics/cgo.h"
#include "pymol/reps/extrude.h"
#include "pymol/objects/coord_set.h"

class CCInOut;
struct RepCartoon;
class nuc_acid_cap;

#define NUCLEIC_NORMAL1 "C3*"
#define NUCLEIC_NORMAL2 "C3'"

enum class ss_t {
  NONE = 0,
  HELIX = 1,
  SHEET = 2,
  NUCLEIC = 3,
};

/**
 * @brief Stores the trace array index range for a cyclic cartoon segment.
 *
 * Used to identify segments where the last residue's CA bonds back to the
 * first residue's CA (e.g. cyclic peptides with a covale bond in _struct_conn).
 */
struct CyclicSegmentInfo {
  int firstIdx;  ///< Index of the first atom in the segment's trace arrays
  int lastIdx;   ///< Index of the last atom in the segment's trace arrays
};

typedef struct nuc_acid_data {
  int na_mode;
  int* nuc_flag;          // whether atom is part of nucleotide
  int a2;                 // defaults to -1
  int nSeg;               // defaults to 0
  const float* v_o_last;  // defaults to nullptr
  int* sptr;
  int* iptr;
  CCInOut* cc;
  int nAt;
  ss_t* ss;
  int putty_flag;
  int* fp;
  float *vptr, *voptr;
  int* ring_anchor;
  int ring_mode, ring_finder, ring_finder_eff;
  int n_ring;
  char alt;
  char next_alt;
} nuc_acid_data;

/* ------------------------------------------------------------------------
 * cartoon_spline: ChimeraX-style path tessellation.
 *
 * Instead of blending each CA->CA interval separately (CartoonGenerateSample,
 * with cartoon_throw/cartoon_power and a CartoonGenerateRefine pass), fit ONE
 * natural cubic spline through the control points of every cartoon segment
 * and place the `sampling` points of each interval on it.  The orientation
 * vector is carried along the interval by parallel transport and twisted with
 * a sigmoid profile so it lands exactly on the next residue's orientation
 * vector (RepCartoonRefineNormals / RepCartoonFlattenSheets have already made
 * consecutive orientation vectors sign-consistent, so no flip test is needed).
 * Everything downstream (ExtrudeComputeTangents, ExtrudeBuildNormals2f, the
 * per-type extruders including the strand arrow) is unchanged.
 * ------------------------------------------------------------------------ */

struct CartoonSplineCoefs {
  std::vector<float> coef;  // 12 floats per atom index a: [axis][c0 c1 c2 c3]
  std::vector<char> valid;  // coef valid for interval a -> a+1
};

// Defined in rep_cartoon_nucleic.cc.
void nuc_acid(PyMOLGlobals* G, nuc_acid_data* ndata, int a, int a1,
              const AtomInfoType* ai, const CoordSet* cs,
              const ObjectMolecule* obj, int set_flags);
int GenerateRepCartoonDrawRings(PyMOLGlobals* G, nuc_acid_data* ndata,
                                ObjectMolecule* obj, CoordSet* cs, CGO* cgo,
                                float ring_width, int cartoon_color,
                                float alpha);

// Defined in rep_cartoon_extrude.cc.
void CartoonGenerateSample(PyMOLGlobals* G, int sampling, int* n_p, float dev,
                           float* vo, float* v1, float* v2, int c1, int c2,
                           float alpha1, float alpha2, int atom_index1,
                           int atom_index2, float power_a, float power_b,
                           float** vc_p, float** valpha_p, unsigned int** vi_p,
                           float** v_p, float** vn_p);
void CartoonGenerateRefine(int refine, int sampling, float* v, float* vn,
                           float* vo, float* sampling_tmp);
CartoonSplineCoefs CartoonBuildSplineCoefs(
    int nAt, const int* seg, const float* pv,
    const std::unordered_map<int, CyclicSegmentInfo>* cyclicSegs);
void CartoonGenerateSampleSpline(PyMOLGlobals* G, int sampling, int* n_p,
                                 const float* coef, const float* vo, int c1,
                                 int c2, float alpha1, float alpha2,
                                 int atom_index1, int atom_index2, float** vc_p,
                                 float** valpha_p, unsigned int** vi_p,
                                 float** v_p, float** vn_p);
int CartoonExtrudeTube(short use_cylinders_for_strands, CExtrude* ex, CGO* cgo,
                       float tube_radius, int tube_quality, cCylCap tube_cap);
int CartoonExtrudePutty(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                        CGO* cgo, CExtrude* ex, int putty_quality,
                        float putty_radius, float* putty_vals, int sampling);
int CartoonExtrudeCircle(CExtrude* ex, CGO* cgo,
                         short use_cylinders_for_strands, int loop_quality,
                         float loop_radius, cCylCap loop_cap, int dash = 0);
void CartoonExtrudeCurvedCylindricalHelix(CExtrude* ex, CGO* cgo, int quality,
                                          float radius, int sampling);
int CartoonExtrudeRect(PyMOLGlobals* G, CExtrude* ex, CGO* cgo, float width,
                       float length, int highlight_color);
int CartoonExtrudeOval(PyMOLGlobals* G, CExtrude* ex, CGO* cgo,
                       short use_cylinders_for_strands, int oval_quality,
                       float oval_width, float oval_length,
                       int highlight_color);
int CartoonExtrudeArrow(PyMOLGlobals* G, CExtrude* ex, CGO* cgo, int sampling,
                        float width, float length, int highlight_color);
int CartoonExtrudeDumbbell(PyMOLGlobals* G, CExtrude* ex, CGO* cgo,
                           int sampling, float dumbbell_width,
                           float dumbbell_length, int highlight_color,
                           int loop_quality, float dumbbell_radius,
                           short use_cylinders_for_strands);

// Defined in rep_cartoon_path.cc.
void RepCartoonComputePuttyValues(ObjectMolecule* obj, float* putty_vals);
std::unordered_map<int, CyclicSegmentInfo> RepCartoonDetectCyclicSegments(
    ObjectMolecule* obj, CoordSet* cs, int nAt, const int* seg, const int* at);
void RepCartoonComputeDifferencesAndNormals(PyMOLGlobals* G, int nAt, int* seg,
                                            float* pv, float* dv, float* nv,
                                            float* dl, int quiet);
void RepCartoonComputeTangents(int nAt, int* seg, float* nv, float* tv);
void RepCartoonComputeDifferencesAndNormalsCyclic(
    PyMOLGlobals* G, int nAt, int* seg, float* pv, float* dv, float* nv,
    float* dl, int quiet,
    const std::unordered_map<int, CyclicSegmentInfo>& cyclicSegs);
void RepCartoonComputeTangentsCyclic(
    int nAt, int* seg, float* nv, float* tv,
    const std::unordered_map<int, CyclicSegmentInfo>& cyclicSegs);
void RepCartoonComputeRoundHelices(nuc_acid_data* ndata, const int nAt,
                                   const int* sptr, const ss_t* ss,
                                   const float* v0, const float* vptr);
void RepCartoonRefineNormals(PyMOLGlobals* G, RepCartoon* I,
                             ObjectMolecule* obj, CoordSet* cs,
                             nuc_acid_data* ndata, const int nAt,
                             const int* seg, const float* tv, float* pvo,
                             float* pva, const ss_t* ss, const float* nv);
void RepCartoonFlattenSheets(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                             nuc_acid_data* ndata, int nAt, const int* sptr,
                             const CCInOut* cc, float* pv, float* pvo,
                             const ss_t* ss, const float* v0, float* tmp,
                             const int* flag_tmp);
void RepCartoonFlattenSheetsRefineTips(PyMOLGlobals* G, ObjectMolecule* obj,
                                       CoordSet* cs, int nAt, int* seg,
                                       const ss_t* ss, float* tv);
void RepCartoonSmoothLoops(PyMOLGlobals* G, ObjectMolecule* obj, CoordSet* cs,
                           nuc_acid_data* ndata, const int nAt, const int* seg,
                           float* pv, const ss_t* ss, float* pvo,
                           const float* /* tv */, float* tmp,
                           const int* flag_tmp);

class CCInOut {
  signed char cc_in{cCartoon_auto /* 0 */};
  signed char cc_out{0};

 public:
  signed char getCCIn() const { return cc_in; }
  signed char getCCOut() const { return cc_out ? cc_out : getCCIn(); }

  void setCCOut(int c) { cc_out = c; }

  // legacy
  void operator=(int c) { cc_in = c; }
  operator int() const { return cc_in; }
  bool operator==(const CCInOut& other) const { return cc_in == other.cc_in; }
};

struct RepCartoon : Rep {
  using Rep::Rep;

  ~RepCartoon() override;

  cRep_t type() const override { return cRepCartoon; }
  void render(RenderInfo* info) override;
  void invalidate(cRepInv_t level) override;
  bool sameVis() const override;

  CGO* ray = nullptr;
  CGO* std = nullptr;
  CGO* preshader = nullptr;

  /**
   * Free the preshader CGO or move to another owner.
   * @post preshader == nullptr
   */
  void disposePreshaderCGO() {
    if (!ray) {
      std::swap(ray, preshader);
    } else {
      CGOFree(preshader);
    }
  }

  char* LastVisib = nullptr;
};

/**
 * nucleic acid cap
 */
class nuc_acid_cap {
  PyMOLGlobals* G;
  nuc_acid_data* ndata;
  CoordSet* cs;
  int idx;

 public:
  int atm;
  AtomInfoType* ai;
  bool enabled;

  // constructor
  nuc_acid_cap(PyMOLGlobals* G_, nuc_acid_data* ndata_, CoordSet* cs_, int mode)
      : G(G_), ndata(ndata_), cs(cs_), idx(0), atm(0), ai(nullptr) {
    enabled = (ndata->na_mode == 4 || ndata->na_mode == mode);
  }

  // set atom indices
  void set(int idx_, int atm_, AtomInfoType* ai_) {
    idx = idx_;
    atm = atm_;
    ai = ai_;
  }

  // check if this cap is enabled and has atom information
  bool active() { return (ai && enabled); }

  // make the cap, if active, and clear atom indices
  bool cap() {
    if (!active()) return false;

    nuc_acid(G, ndata, idx, atm, ai, cs, cs->Obj, false);
    set(-1, -1, nullptr);
    return true;
  }
};
