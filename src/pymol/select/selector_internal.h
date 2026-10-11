/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2004 by Warren Lyford Delano of DeLano Scientific.
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

// Helpers shared by the selector*.cc files; not part of the selector API.

#pragma once

#include "pymol/base/err.h"
#include "pymol/base/word.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/select/selector_def.h"

using SelectorInfoIter_t = decltype(CSelectorManager::Info)::iterator;

/// Helper type which replaces `int*` return types
typedef std::unique_ptr<int[]> sele_array_t;

inline void sele_array_calloc(sele_array_t& sele, size_t count) {
  sele.reset(new int[count]());
}

struct EvalElem {
  int level, imp_op_level;
  int type; /* 0 = value 1 = operation 2 = pre-operation */
  unsigned int code;
  std::string m_text;
  sele_array_t sele;

  // Helpers for refactoring `sele` type
  int* sele_data() { return sele.get(); }
  void sele_free() { sele.reset(); }
  void sele_calloc(size_t count) { sele_array_calloc(sele, count); }

  // TODO replace with pymol::Error handling
  void sele_check_ok(int& ok) { CHECKOK(ok, sele_data()); }
  void sele_err_chk_ptr(PyMOLGlobals* G) { ErrChkPtr(G, sele_data()); }

  /// read-only access to text
  const char* text() const { return m_text.c_str(); }
};

// Defined in selector_eval.cc.
pymol::Result<sele_array_t> SelectorSelect(PyMOLGlobals* G, const char* sele,
                                           int state, SelectorID_t domain,
                                           int quiet);
int SelectorModulate1(PyMOLGlobals* G, EvalElem* base, int state);
int SelectorSelect0(PyMOLGlobals* G, EvalElem* base);
int SelectorSelect2(PyMOLGlobals* G, EvalElem* base, int state);
int SelectorLogic1(PyMOLGlobals* G, EvalElem* base, int state);
int SelectorLogic2(PyMOLGlobals* G, EvalElem* base);
int SelectorOperator22(PyMOLGlobals* G, EvalElem* base, int state);

// Defined in selector_parse.cc.
pymol::Result<sele_array_t> SelectorEvaluate(PyMOLGlobals* G,
                                             std::vector<std::string>& word,
                                             int state, int quiet);
std::vector<std::string> SelectorParse(PyMOLGlobals* G, const char* s);

// Defined in selector.cc.
int SelectorEmbedSelection(PyMOLGlobals* G, const int* atom,
                           pymol::zstring_view name, ObjectMolecule* obj,
                           int no_dummies, int exec_manage);
int* SelectorGetIndexVLA(PyMOLGlobals* G, SelectorID_t sele);
void SelectorUpdateTableSingleObject(PyMOLGlobals* G, ObjectMolecule* obj,
                                     int req_state, bool no_dummies = false);

#define STYP_OPR1 1
#define STYP_OPR2 2
#define STYP_SEL0 3
#define STYP_SEL1 4
#define STYP_SEL2 5
#define STYP_LIST 6
#define STYP_PRP1 7
#define STYP_SEL3 8
#define STYP_OP22 9 /* sele oper arg1 arg2 sele */

/*                  code   |   type    | priority */
#define SELE_NOT1 (0x0100 | STYP_OPR1 | 0x70)
#define SELE_BYR1 (0x0200 | STYP_OPR1 | 0x20)
#define SELE_AND2 (0x0300 | STYP_OPR2 | 0x60)
#define SELE_OR_2 (0x0400 | STYP_OPR2 | 0x40)
#define SELE_IN_2 (0x0500 | STYP_OPR2 | 0x40)
#define SELE_ALLz (0x0600 | STYP_SEL0 | 0x90)
#define SELE_NONz (0x0700 | STYP_SEL0 | 0x90)
#define SELE_HETz (0x0800 | STYP_SEL0 | 0x80)
#define SELE_HYDz (0x0900 | STYP_SEL0 | 0x90)
#define SELE_VISz (0x0A00 | STYP_SEL0 | 0x90)
#define SELE_ARD_ (0x0B00 | STYP_PRP1 | 0x30)
#define SELE_EXP_ (0x0C00 | STYP_PRP1 | 0x30)
#define SELE_NAMs (0x0D00 | STYP_SEL1 | 0x80)
#define SELE_ELEs (0x0E00 | STYP_SEL1 | 0x80)
#define SELE_RSIs (0x0F00 | STYP_SEL1 | 0x80)
#define SELE_CHNs (0x1000 | STYP_SEL1 | 0x80)
#define SELE_SEGs (0x1100 | STYP_SEL1 | 0x80)
#define SELE_MODs (0x1200 | STYP_SEL1 | 0x80)
#define SELE_IDXs (0x1300 | STYP_SEL1 | 0x80)
#define SELE_RSNs (0x1400 | STYP_SEL1 | 0x80)
#define SELE_SELs (0x1500 | STYP_SEL1 | 0x80)
#define SELE_BVLx (0x1600 | STYP_SEL2 | 0x80)
#define SELE_ALTs (0x1700 | STYP_SEL1 | 0x80)
#define SELE_FLGs (0x1800 | STYP_SEL1 | 0x80)
#define SELE_GAP_ (0x1900 | STYP_PRP1 | 0x80)
#define SELE_TTYs (0x1A00 | STYP_SEL1 | 0x80)
#define SELE_NTYs (0x1B00 | STYP_SEL1 | 0x80)
#define SELE_PCHx (0x1C00 | STYP_SEL2 | 0x80)
#define SELE_FCHx (0x1D00 | STYP_SEL2 | 0x80)
#define SELE_ID_s (0x1E00 | STYP_SEL1 | 0x80)
#define SELE_BNDz (0x1F00 | STYP_SEL0 | 0x80)
#define SELE_LIK2 (0x2000 | STYP_OPR2 | 0x40)
#define SELE_NGH1 (0x2100 | STYP_OPR1 | 0x20)
#define SELE_QVLx (0x2200 | STYP_SEL2 | 0x80)
#define SELE_BYO1 (0x2300 | STYP_OPR1 | 0x20)
#define SELE_SSTs (0x2400 | STYP_SEL1 | 0x80)
#define SELE_STAs (0x2500 | STYP_SEL1 | 0x80)
#define SELE_PREz (0x2500 | STYP_SEL0 | 0x80)
#define SELE_WIT_ (0x2600 | STYP_OP22 | 0x30)
#define SELE_ORIz (0x2700 | STYP_SEL0 | 0x90)
#define SELE_CENz (0x2800 | STYP_SEL0 | 0x90)
#define SELE_ENAz (0x2900 | STYP_SEL0 | 0x90)
#define SELE_REPs (0x2A00 | STYP_SEL1 | 0x80)
#define SELE_COLs (0x2B00 | STYP_SEL1 | 0x80)
#define SELE_HBDs (0x2C00 | STYP_SEL0 | 0x80)
#define SELE_HBAs (0x2D00 | STYP_SEL0 | 0x80)
#define SELE_BYC1 (0x2E00 | STYP_OPR1 | 0x20)
#define SELE_BYS1 (0x2F00 | STYP_OPR1 | 0x20)
#define SELE_BYM1 (0x3000 | STYP_OPR1 | 0x20)
#define SELE_BYF1 (0x3100 | STYP_OPR1 | 0x20)
#define SELE_EXT_ (0x3200 | STYP_PRP1 | 0x30)
#define SELE_BON1 (0x3300 | STYP_OPR1 | 0x50)
#define SELE_FST1 (0x3400 | STYP_OPR1 | 0x30)
#define SELE_CAS1 (0x3500 | STYP_OPR1 | 0x30)
#define SELE_BEY_ (0x3600 | STYP_OP22 | 0x30)
#define SELE_POLz (0x3700 | STYP_SEL0 | 0x90)
#define SELE_SOLz (0x3800 | STYP_SEL0 | 0x90)
#define SELE_ORGz (0x3900 | STYP_SEL0 | 0x90)
#define SELE_INOz (0x3A00 | STYP_SEL0 | 0x90)
#define SELE_GIDz (0x3B00 | STYP_SEL0 | 0x90)
#define SELE_RNKs (0x3C00 | STYP_SEL1 | 0x80)
#define SELE_PEPs (0x3D00 | STYP_SEL1 | 0x80)
#define SELE_ACCz (0x3E00 | STYP_SEL0 | 0x90)
#define SELE_DONz (0x3F00 | STYP_SEL0 | 0x90)
#define SELE_LST1 (0x4000 | STYP_OPR1 | 0x30)
#define SELE_NTO_ (0x4100 | STYP_OP22 | 0x30)
#define SELE_CCLs (0x4200 | STYP_SEL1 | 0x80)
#define SELE_RCLs (0x4300 | STYP_SEL1 | 0x80)
#define SELE_PTDz (0x4400 | STYP_SEL0 | 0x90)
#define SELE_MSKz (0x4500 | STYP_SEL0 | 0x90)
#define SELE_IOR2 (0x4600 | STYP_OPR2 | 0x10)
#define SELE_FXDz (0x4700 | STYP_SEL0 | 0x90)
#define SELE_RSTz (0x4800 | STYP_SEL0 | 0x90)
#define SELE_ANT2 (0x4900 | STYP_OPR2 | 0x60)
#define SELE_BYX1 (0x4A00 | STYP_OPR1 | 0x20)
#define SELE_STRO (0x4B00 | STYP_SEL1 | 0x80)
#define SELE_METz (0x4C00 | STYP_SEL0 | 0x90)
#define SELE_BB_z (0x4D00 | STYP_SEL0 | 0x90)
#define SELE_SC_z (0x4E00 | STYP_SEL0 | 0x90)
#define SELE_PROP (0x4F00 | STYP_SEL3 | 0x80)
#define SELE_XVLx (0x5000 | STYP_SEL2 | 0x80)
#define SELE_YVLx (0x5100 | STYP_SEL2 | 0x80)
#define SELE_ZVLx (0x5200 | STYP_SEL2 | 0x80)
#define SELE_CUST (0x5300 | STYP_SEL1 | 0x80)
#define SELE_RING (0x5400 | STYP_OPR1 | 0x20)
#define SELE_LABs (0x5500 | STYP_SEL1 | 0x80)
#define SELE_PROz (0x5600 | STYP_SEL0 | 0x90)
#define SELE_NUCz (0x5700 | STYP_SEL0 | 0x90)
#define SELE_DESz (0x5800 | STYP_SEL0 | 0x90)

// Defined in selector.cc.
void SelectorManagerInsertMember(CSelectorManager& self, AtomInfoType& ai,
                                 int sele, int tag = 1);
int SelectorGetObjAtmOffset(CSelector* I, const ObjectMolecule* obj,
                            int offset);
extern WordKeyValue Keyword[];

// Defined in selector.cc.
SelectorInfoIter_t SelectGetInfoIter(PyMOLGlobals* G, const char* name,
                                     ov_size minMatch, int ignCase);
void SelectorDeleteSeleAtIter(PyMOLGlobals* G, SelectorInfoIter_t it);

// Defined in selector_eval.cc.
pymol::Result<> SelectorSelect1(PyMOLGlobals* G, EvalElem* base, int quiet);
pymol::Result<> SelectorSelect3(PyMOLGlobals* G, EvalElem* base, int state);
