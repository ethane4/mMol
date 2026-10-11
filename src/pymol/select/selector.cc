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

#include "pymol/base/feedback.h"

#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/ui/seq.h"
#include "pymol/executive/editor.h"
#include "pymol/chem/mol2_typing.h"

#include "pymol/base/parse.h"

#include "pymol/base/util2.h"

#include "pymol/objects/property.h"

#include "pymol/select/selector_internal.h"

/**
 * Prefix for temporary selections.
 *
 * Must be a valid name for the selection language, so it can't contain any
 * characters which are selection operators (e.g. "_!" is not valid).
 *
 * Ideally a prefix which can't be used by the user. This is currently
 * only true if the `validate_object_names` setting is off (default).
 */
#define cSelectorTmpPrefix "_#"

typedef struct {
  int depth1;
  int depth2;
  int depth3;
  int depth4;
  int sum;
  int frag;
} WalkDepthRec;
static void SelectorPurgeMembers(PyMOLGlobals* G, SelectorID_t sele);
static int* SelectorGetIndexVLAImpl(PyMOLGlobals* G, CSelector* I, int sele);
static void SelectorClean(PyMOLGlobals* G);

static sele_array_t SelectorGetSeleArrayForAtomIndices(CSelector* I,
                                                       ObjectMolecule* obj,
                                                       const int* idx,
                                                       int n_idx,
                                                       bool numbered_tags);
bool SelectorAtomIterator::next() {
  if ((++a) >= selector->Table.size()) return false;

  auto& table_a = selector->Table[a];

  atm = table_a.atom;
  obj = selector->Obj[table_a.model];

  return true;
}

/**
 * Add atom `ai` to selection `sele`
 */
void SelectorManagerInsertMember(CSelectorManager& self, AtomInfoType& ai,
                                 int sele, int tag) {
  int m;
  if (self.FreeMember > 0) {
    m = self.FreeMember;
    self.FreeMember = self.Member[m].next;
  } else {
    m = self.Member.size();
    self.Member.emplace_back();
  }
  self.Member[m].selection = sele;
  self.Member[m].tag = tag;
  self.Member[m].next = ai.selEntry;
  ai.selEntry = m;
}
static void SelectorGetUniqueTmpName(PyMOLGlobals* G, char* out) {
  sprintf(out, "%s%d", cSelectorTmpPrefix, G->SelectorMgr->NSelection);
}

static bool SelectorIsTmp(pymol::zstring_view name) {
  assert(name);
  return name.starts_with(cSelectorTmpPrefix);
}
int SelectorGetObjAtmOffset(CSelector* I, const ObjectMolecule* obj,
                            int offset) {
  if (I->SeleBaseOffsetsValid) {
    return obj->SeleBase + offset;
  } else {
    ov_diff stop_below = obj->SeleBase;
    ov_diff stop_above = I->Table.size() - 1;
    int result = stop_below;
    int step = offset;
    int cur;
    int proposed;
    int prior1 = -1, prior2 = -1;

    /* non-linear hunt to find atom */

    result = stop_below;
    cur = I->Table[result].atom;
    while (step > 1) {
      if (cur < offset) {
        stop_below = result + 1;
        while (step > 1) {
          proposed = result + step;
          if (proposed <= stop_above) {
            if (I->Obj[I->Table[proposed].model] == obj) {
              if (proposed == prior1) {
                proposed--;
                step--; /* guarantee progress (avoid flip flop) */
              }
              result = prior1 = proposed;
              break;
            } else if (stop_above > proposed) {
              stop_above = proposed - 1;
            }
          }
          step = (step >> 1);
        }
      } else if (cur > offset) {
        stop_above = result - 1;
        while (step > 1) {
          proposed = result - step;
          if (proposed >= stop_below) {
            if (I->Obj[I->Table[proposed].model] == obj) {
              if (proposed == prior2) {
                proposed++;
                step--; /* guarantee progress (avoid flip flop) */
              }
              result = prior2 = proposed;
              break;
            }
          }
          step = (step >> 1);
        }
      } else
        return result;
      cur = I->Table[result].atom;
      if (cur == offset) return result;
    }

    {
      /* failsafe / linear search */
      int dir = 1;
      if (cur > offset) dir = -1;
      while (1) { /* TODO: optimize this search algorithm! */
        if (cur == offset) return result;
        if (dir > 0) {
          if (result >= stop_above) break;
          result++;
        } else {
          if (result <= stop_below) break;
          result--;
        }
        if (I->Obj[I->Table[result].model] != obj) break;
        cur = I->Table[result].atom;
      }
    }
  }
  return -1;
}

#define SEL_PREMAX 0x8

WordKeyValue Keyword[] = {
    {"not", SELE_NOT1},
    {"!", SELE_NOT1},

    {"neighbor", SELE_NGH1},
    {"nbr;", SELE_NGH1}, /* deprecated */
    {"nbr.", SELE_NGH1},

    {"byfragment", SELE_BYF1},
    {"byfrag", SELE_BYF1},
    {"bf.", SELE_BYF1},

    {"byresidue", SELE_BYR1},
    {"byresi", SELE_BYR1}, /* unofficial */
    {"byres", SELE_BYR1},
    {"br;", SELE_BYR1}, /* deprecated */
    {"br.", SELE_BYR1},
    {"b;", SELE_BYR1}, /* deprecated */

    {"bychain", SELE_BYC1},
    {"bc.", SELE_BYC1},

    {"byobject", SELE_BYO1},
    {"byobj", SELE_BYO1},
    {"bo;", SELE_BYO1}, /* deprecated */
    {"bo.", SELE_BYO1},

    {"bound_to", SELE_BON1},
    {"bto.", SELE_BON1},

    {"bymolecule", SELE_BYM1},
    {"bymol", SELE_BYM1},
    {"bm.", SELE_BYM1},

    {"bysegment", SELE_BYS1},
    {"byseg", SELE_BYS1},
    {"bysegi", SELE_BYS1}, /* unofficial */
    {"bs.", SELE_BYS1},

    {"bycalpha", SELE_CAS1},
    {"bca.", SELE_CAS1},

    {"first", SELE_FST1},
    {"last", SELE_LST1},

    {"and", SELE_AND2},
    {"&", SELE_AND2},
    {"or", SELE_OR_2},
    {"+", SELE_OR_2}, /* added to mitigate damage caused by the obj1+obj2 parser
                         bug */
    {"-", SELE_ANT2}, /* added to provide natural complement to the above: an
                         AND NOT or SUBTRACT operation */
    {"|", SELE_OR_2},
    {"in", SELE_IN_2},

    {"like", SELE_LIK2},
    {"l;", SELE_LIK2},
    {"l.", SELE_LIK2},

    {cKeywordAll, SELE_ALLz}, /* 0 parameter */
    {"*", SELE_ALLz},         /* 0 parameter */

    {cKeywordNone, SELE_NONz}, /* 0 parameter */
    {"hetatm", SELE_HETz},     /* 0 parameter */
    {"het", SELE_HETz},        /* 0 parameter */

    {"hydrogens", SELE_HYDz}, /* 0 parameter */
    {"hydro", SELE_HYDz},     /* 0 parameter */
    {"h;", SELE_HYDz},        /* deprecated */
    {"h.", SELE_HYDz},        /* 0 parameter */

    {"hba.", SELE_HBAs},
    {"hbd.", SELE_HBDs},

    {"visible", SELE_VISz}, /* 0 parameter */
    {"v;", SELE_VISz},      /* 0 parameter */
    {"v.", SELE_VISz},      /* 0 parameter */

    {"around", SELE_ARD_}, /* 1 parameter */
    {"a;", SELE_ARD_},     /* deprecated */
    {"a.", SELE_ARD_},     /* 1 parameter */

    {"expand", SELE_EXP_}, /* 1 parameter */
    {"x;", SELE_EXP_},     /* 1 parameter */
    {"x.", SELE_EXP_},     /* 1 parameter */

    {"extend", SELE_EXT_}, /* 1 parameter */
    {"xt.", SELE_EXT_},    /* 1 parameter */

    {"name", SELE_NAMs},
    {"n;", SELE_NAMs}, /* deprecated */
    {"n.", SELE_NAMs},

    {"symbol", SELE_ELEs},
    {"element", SELE_ELEs},
    {"elem", SELE_ELEs},
    {"e;", SELE_ELEs}, /* deprecated */
    {"e.", SELE_ELEs},

    {"enabled", SELE_ENAz},

    {"residue", SELE_RSIs},
    {"resi", SELE_RSIs},
    {"resident", SELE_RSIs},
    {"resid", SELE_RSIs},
    {"i;", SELE_RSIs}, /* deprecated */
    {"i.", SELE_RSIs},

    {"rep", SELE_REPs},

    {"color", SELE_COLs},
    {"cartoon_color", SELE_CCLs},
    {"ribbon_color", SELE_RCLs},

    {"altloc", SELE_ALTs},
    {"alt", SELE_ALTs},

    {"flag", SELE_FLGs},
    {"f;", SELE_FLGs}, /* deprecated */
    {"f.", SELE_FLGs},

    {"gap", SELE_GAP_},

    {"partial_charge", SELE_PCHx},
    {"pc;", SELE_PCHx}, /* deprecated */
    {"pc.", SELE_PCHx},

    {"masked", SELE_MSKz},
    {"msk.", SELE_MSKz},

    {"protected", SELE_PTDz},

    {"formal_charge", SELE_FCHx},
    {"fc;", SELE_FCHx}, /* deprecated */
    {"fc.", SELE_FCHx},

    {"numeric_type", SELE_NTYs},
    {"nt;", SELE_NTYs}, /* deprecated */
    {"nt.", SELE_NTYs},

    {"text_type", SELE_TTYs},
    {"custom", SELE_CUST},
    {"tt;", SELE_TTYs}, /* deprecated */
    {"tt.", SELE_TTYs},

    {"chain", SELE_CHNs},
    {"c;", SELE_CHNs}, /* deprecated */
    {"c.", SELE_CHNs},

    {cKeywordCenter, SELE_CENz},
    {"bonded", SELE_BNDz},

    {"segment", SELE_SEGs},
    {"segid", SELE_SEGs},
    {"segi", SELE_SEGs},
    {"s;", SELE_SEGs}, /* deprecated */
    {"s.", SELE_SEGs},

    {"ss", SELE_SSTs},

    {"state", SELE_STAs},

    {"object", SELE_MODs},
    {"o.", SELE_MODs},

    {cKeywordOrigin, SELE_ORIz},

    {"model", SELE_MODs},
    {"m;", SELE_MODs}, /* deprecated */
    {"m.", SELE_MODs},

    {"index", SELE_IDXs},
    {"idx.", SELE_IDXs},

    {"id", SELE_ID_s},
    {"ID", SELE_ID_s},
    {"rank", SELE_RNKs},

    {"within", SELE_WIT_},
    {"w.", SELE_WIT_},

    {"near_to", SELE_NTO_},
    {"nto.", SELE_NTO_},

    {"beyond", SELE_BEY_},
    {"be.", SELE_BEY_},

    {"donors", SELE_DONz},
    {"don.", SELE_DONz},

    {"acceptors", SELE_ACCz},
    {"acc.", SELE_ACCz},

    {"delocalized", SELE_DESz},
    {"deloc.", SELE_DESz},

    {"pepseq", SELE_PEPs},
    {"ps.", SELE_PEPs},

    /*
       {  "nucseq",  SELE_NUCs },
       {  "ns.",      SELE_NUCs },
      */

    {"fixed", SELE_FXDz},
    {"fxd.", SELE_FXDz},

    {"restrained", SELE_RSTz},
    {"rst.", SELE_RSTz},

    {"polymer", SELE_POLz},
    {"pol.", SELE_POLz},

    {"polymer.protein", SELE_PROz},
    {"polymer.nucleic", SELE_NUCz},

    {"organic", SELE_ORGz},
    {"org.", SELE_ORGz},

    {"inorganic", SELE_INOz},
    {"ino.", SELE_INOz},

    {"solvent", SELE_SOLz},
    {"sol.", SELE_SOLz},

    {"guide", SELE_GIDz},

    {"present", SELE_PREz},
    {"pr.", SELE_PREz},

    {"resname", SELE_RSNs},
    {"resn", SELE_RSNs},
    {"r;", SELE_RSNs}, /* deprecated */
    {"r.", SELE_RSNs},

    {"%", SELE_SELs},
    {"b", SELE_BVLx}, /* 2 operand selection operator */
    {"q", SELE_QVLx}, /* 2 operand selection operator */

    {"stereo", SELE_STRO},

    {"bycell", SELE_BYX1},

    {"metals", SELE_METz}, /* 0 parameter */

    {"backbone", SELE_BB_z},
    {"bb.", SELE_BB_z},

    {"sidechain", SELE_SC_z},
    {"sc.", SELE_SC_z},

    {"p.", SELE_PROP},

    {"x", SELE_XVLx},
    {"y", SELE_YVLx},
    {"z", SELE_ZVLx},

    {"byring", SELE_RING},
    {"label", SELE_LABs},

    {"", 0}};

int SelectorRenameObjectAtoms(PyMOLGlobals* G, ObjectMolecule* obj,
                              SelectorID_t sele, bool force,
                              bool update_table) {
  int result = 0;
  int obj_nAtom = obj->NAtom;

  if (update_table) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  }
  if (obj_nAtom) {
    int* flag = pymol::calloc<int>(obj_nAtom);
    if (!flag) {
      result = -1;
    } else {
      const AtomInfoType* ai = obj->AtomInfo.data();
      int a;
      for (a = 0; a < obj_nAtom; a++) {
        if (SelectorIsMember(G, ai->selEntry, sele)) {
          flag[a] = true;
          result = true;
        }
        ai++;
      }
      if (!result && !force) {
        // nothing selected, no need to continue
        return 0;
      }
      result = ObjectMoleculeRenameAtoms(obj, flag, force);
    }
    FreeP(flag);
  }
  return result;
}

int SelectorClassifyAtoms(PyMOLGlobals* G, int sele, int preserve,
                          ObjectMolecule* only_object) {
  CSelector* I = G->Selector;
  ObjectMolecule *obj, *obj0, *obj1 = nullptr;
  int a, aa, at, a0, a1;
  AtomInfoType *ai, *last_ai = nullptr, *ai0, *ai1;
  unsigned int mask;
  int n_dummies = 0;

  int auto_show_classified =
      SettingGetGlobal_i(G, cSetting_auto_show_classified);
  int auto_show_mask = (auto_show_classified != 2) ? 0 : cRepBitmask;

  int visRep_organic = cRepCylBit | cRepNonbondedSphereBit;
  int visRep_inorganic = cRepSphereBit;
  int visRep_polymer = cRepCartoonBit;

  const lexborrow_t lex_pseudo = LexBorrow(G, "pseudo");

  // detect large systems
  if (auto_show_classified == -1 && only_object &&
      only_object->NAtom *
              (only_object->DiscreteFlag ? 1 : only_object->NCSet) >
          5e5) {
    auto_show_classified = 3;
  }

  if (auto_show_classified == 3) {
    visRep_organic = visRep_inorganic = cRepLineBit | cRepNonbondedBit;
    visRep_polymer = cRepRibbonBit;
  }

  if (only_object) {
    SelectorUpdateTableSingleObject(G, only_object,
                                    cSelectorUpdateTableAllStates, true);
    n_dummies = 0;
  } else {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
    n_dummies = cNDummyAtoms;
  }
  a = 0;
  while (a < I->Table.size()) {
    obj = I->Obj[I->Table[a].model];
    at = I->Table[a].atom;
    ai = obj->AtomInfo + at;

    if (SelectorIsMember(G, ai->selEntry, sele) &&
        ((!AtomInfoSameResidueP(G, ai, last_ai)))) {
      AtomInfoType* guide_atom = nullptr;

      /* delimit residue */

      a0 = a - 1;
      while (a0 >= n_dummies) {
        obj0 = I->Obj[I->Table[a0].model];
        if (obj0 != obj) break;
        ai0 = obj0->AtomInfo + I->Table[a0].atom;
        if (!AtomInfoSameResidue(G, ai0, ai)) break;
        a0--;
      }

      a1 = a + 1;
      while (a1 < I->Table.size()) {
        obj1 = I->Obj[I->Table[a1].model];
        if (obj1 != obj) break;
        ai1 = obj1->AtomInfo + I->Table[a1].atom;
        if (!AtomInfoSameResidue(G, ai1, ai)) break;
        a1++;
      }

      a0++;
      a1--;

      mask = 0;
      if (!ai->hetatm && AtomInfoKnownProteinResName(LexStr(G, ai->resn)))
        mask = cAtomFlag_polymer | cAtomFlag_protein;
      else if ((!ai->hetatm || AtomInfoKnownPNAResName(LexStr(G, ai->resn))) &&
               AtomInfoKnownNucleicResName(LexStr(G, ai->resn)))
        mask = cAtomFlag_polymer | cAtomFlag_nucleic;
      else if (AtomInfoKnownWaterResName(G, LexStr(G, ai->resn)))
        mask = cAtomFlag_solvent;
      else {
        /* does this residue have a canonical atoms? */

        bool found_only_h = true;
        int found_ca = false;
        int found_n = false;
        int found_c = false;
        int found_o = false;
        int found_oh2 = false;
        int found_carbon = false;
        int found_cn_bond = false;
        int found_nc_bond = false;
        int found_o3_bond = false;
        int found_o3star = false;
        int found_c3star = false;
        int found_c4star = false;
        int found_c5star = false;
        int found_o5star = false;
        int found_p_bond = false;

        ai0 = obj->AtomInfo + I->Table[a0].atom;
        for (aa = a0; aa <= a1; aa++) {
          if (ai0->protons == cAN_C) {
            const char* name = LexStr(G, ai0->name);
            found_carbon = true;
            switch (name[0]) {
              case 'C':
                switch (name[1]) {
                  case 0:
                    found_c = true;
                    found_cn_bond = ObjectMoleculeIsAtomBondedToName(
                        obj, I->Table[aa].atom, "N", 0);
                    break;
                  case 'A':
                    switch (name[2]) {
                      case 0:
                        found_ca = true;
                        guide_atom = ai0;
                        break;
                    }
                  case '3':
                    switch (name[2]) {
                      case '*':
                      case '\'':
                        guide_atom = ai0;
                        found_c3star = true;
                        break;
                    }
                    break;
                  case '4':
                    switch (name[2]) {
                      case '*':
                      case '\'':
                        found_c4star = true;
                        break;
                    }
                    break;
                  case '5':
                    switch (name[2]) {
                      case '*':
                      case '\'':
                        found_c5star = true;
                        break;
                    }
                    break;
                }
            }
          } else if (ai0->protons == cAN_N) {
            const char* name = LexStr(G, ai0->name);
            switch (name[0]) {
              case 'N':
                switch (name[1]) {
                  case 0:
                    found_n = true;
                    found_nc_bond = ObjectMoleculeIsAtomBondedToName(
                        obj, I->Table[aa].atom, "C", 0);
                    break;
                }
            }
          } else if (ai0->protons == cAN_O) {
            const char* name = LexStr(G, ai0->name);
            switch (name[0]) {
              case 'O':
                switch (name[1]) {
                  case 0:
                    found_o = true;
                    break;
                  case 'H':
                    switch (name[2]) {
                      case '2':
                        found_oh2 = true;
                        break;
                    }
                  case '3':
                    switch (name[2]) {
                      case '*':
                      case '\'':
                        found_o3star = true;
                        found_o3_bond = ObjectMoleculeIsAtomBondedToName(
                            obj, I->Table[aa].atom, "P", 0);
                        break;
                    }
                    break;
                  case '5':
                    switch (name[2]) {
                      case '*':
                      case '\'':
                        found_o5star = true;
                        break;
                    }
                    break;
                }
            }
          } else if (ai0->protons == cAN_P) {
            const char* name = LexStr(G, ai0->name);
            switch (name[0]) {
              case 'P':
                switch (name[1]) {
                  case 0:
                    found_p_bond = (ObjectMoleculeIsAtomBondedToName(
                                        obj, I->Table[aa].atom, "O3*", 0) ||
                                    ObjectMoleculeIsAtomBondedToName(
                                        obj, I->Table[aa].atom, "O3'", 0));
                    break;
                }
            }
          }

          if (!ai0->isHydrogen()) {
            found_only_h = false;
          }

          ai0++;
        }

        if (found_ca && found_n && found_c && found_o &&
            (found_cn_bond || found_nc_bond)) {
          mask = cAtomFlag_polymer | cAtomFlag_protein;
        } else if (found_o3star && found_c3star && found_c4star &&
                   found_c5star && found_o5star &&
                   (found_o3_bond || found_p_bond)) {
          mask = cAtomFlag_polymer | cAtomFlag_nucleic;
        } else if (found_carbon)
          mask = cAtomFlag_organic;
        else if ((found_o || found_oh2) && (a1 == a0))
          mask = cAtomFlag_solvent;
        else if (!found_only_h) {
          // exclude hydrogens, they get misclassified as
          // 'inorganic' if they are not sorted
          mask = cAtomFlag_inorganic;
        }
      }

      /* mark which atoms we can write to */

      ai0 = obj->AtomInfo + I->Table[a0].atom;
      if (preserve) {
        printf("NOT IMPLEMENTED\n");
      } else {
        auto visRep_polymer_obj = visRep_polymer;
        if (obj->NAtom < 50) {
          // prevent single residue objects from disappearing
          visRep_polymer_obj |= visRep_organic;
        }

        for (aa = a0; aa <= a1; aa++) {
          if (SelectorIsMember(G, ai0->selEntry, sele)) {
            // apply styles if atom was unclassified
            if (auto_show_classified && !(ai0->flags & cAtomFlag_class)) {
              if (mask & cAtomFlag_organic) {
                ai0->visRep = (ai0->visRep & auto_show_mask) | visRep_organic;
              } else if (mask & cAtomFlag_inorganic) {
                ai0->visRep = (ai0->visRep & auto_show_mask) | visRep_inorganic;
              } else if (mask & cAtomFlag_polymer) {
                ai0->visRep =
                    (ai0->visRep & auto_show_mask) | visRep_polymer_obj;
              }

              // hide Desmond virtual sites and off-centered partial charges
              if (ai0->name == lex_pseudo) {
                ai0->visRep = 0;
              }
            }
            ai0->flags = (ai0->flags & cAtomFlag_class_mask) | mask;
          }
          ai0++;
        }
      }

      if ((mask & cAtomFlag_polymer)) {
        AtomInfoType* guide_atom_c3 = nullptr;
        ai0 = obj->AtomInfo + I->Table[a0].atom;
        for (aa = a0; !guide_atom && aa <= a1; aa++) {
          if (ai0->protons == cAN_C) {
            const char* name = LexStr(G, ai0->name);
            switch (name[0]) {
              case 'C':
                switch (name[1]) {
                  case 'A':
                    switch (name[2]) {
                      case 0:
                        guide_atom = ai0;
                        break;
                    }
                    break;
                  case '4':
                    switch (
                        name[2]) { /* use C4* as guide atom for nucleic acids */
                      case '*':
                      case '\'':
                        guide_atom = ai0;
                        break;
                    }
                    break;
                  case '3':
                    if ((mask & cAtomFlag_nucleic) && !guide_atom_c3) {
                      switch (name[2]) { /* C3' as fallback guide for nucleic
                                            acids without C4' (PNA) */
                        case '*':
                        case '\'':
                          guide_atom_c3 = ai0;
                          break;
                      }
                    }
                    break;
                }
            }
          }
          ai0++;
        }
        if (!guide_atom && guide_atom_c3) {
          guide_atom = guide_atom_c3;
        }
      }

      if (guide_atom) guide_atom->flags |= cAtomFlag_guide;

      if (a1 > (a + 1)) a = a1;
    }
    a++;
  }
  return true;
}

MapType* SelectorGetSpacialMapFromSeleCoord(PyMOLGlobals* G, int sele,
                                            int state, float cutoff,
                                            float** coord_vla) {
  int* index_vla = nullptr;
  float* coord = nullptr;
  int n, nc = 0;
  MapType* result = nullptr;
  if (sele < 0)
    return nullptr;
  else {
    auto ptr = std::make_unique<CSelector>(G, G->SelectorMgr);
    CSelector mapSele(G, G->SelectorMgr);
    auto I = &mapSele;
    SelectorUpdateTableImpl(G, I, state, -1);
    index_vla = SelectorGetIndexVLAImpl(G, I, sele);

    if (index_vla) {
      n = VLAGetSize(index_vla);
      if (n) coord = VLAlloc(float, n * 3);
      if (coord) {
        int i, a;
        int st, sta;
        ObjectMolecule* obj;
        CoordSet* cs;
        int at;
        int idx;
        for (i = 0; i < n; i++) {
          a = index_vla[i];

          obj = I->Obj[I->Table[a].model];
          at = +I->Table[a].atom;
          for (st = 0; st < I->NCSet; st++) {
            if ((state < 0) || (st == state)) {
              sta = st;
              if (sta < obj->NCSet)
                cs = obj->CSet[sta];
              else
                cs = nullptr;
              if (cs) {
                idx = cs->atmToIdx(at);
              } else {
                idx = -1;
              }
              if (idx >= 0) {
                VLACheck(coord, float, nc * 3 + 2);
                const float* src = cs->coordPtr(idx);
                float* dst = coord + 3 * nc;
                copy3f(src, dst);
                nc++;
              }
            }
          }
        }
        if (nc) {
          result = new MapType(G, cutoff, coord, nc, nullptr);
        }
      }
    }
  }
  VLAFreeP(index_vla);
  if (coord) VLASize(coord, float, nc * 3);
  *(coord_vla) = coord;
  return (result);
}

SelectorInfoIter_t SelectGetInfoIter(PyMOLGlobals* G, const char* name,
                                     ov_size minMatch, int ignCase) {
  auto& Info = G->SelectorMgr->Info;
  auto end_offset = Info.end();

  while (name[0] == '?') name++;

  for (auto offset = Info.begin(); offset != end_offset; ++offset) {
    if (offset->name == name) {
      return offset;
    }
  }
  /* not found, so try partial/ignored-case match */
  int best_match = -1;
  auto best_offset = end_offset;

  for (auto offset = Info.begin(); offset != end_offset; ++offset) {
    int wm = WordMatch(G, name, offset->name.c_str(), ignCase);
    if (wm < 0) { /* exact match is always good */
      best_offset = offset;
      best_match = wm;
      break;
    }
    if (wm > 0) {
      if (best_match < wm) {
        best_match = wm;
        best_offset = offset;
      } else if (best_match == wm) { /* uh oh -- ambiguous match */
        best_offset = end_offset;
      }
    }
  }
  if ((best_match < 0) || (best_match > (int)minMatch)) return best_offset;

  return end_offset;
}

/**
 * Only called once: When doing `cmd.delete("all")`
 */
void SelectorDefragment(PyMOLGlobals* G) {
  CSelector* S = G->Selector;
  auto I = S->mgr;
  /* restore new member ordering so that CPU can continue to get good cache hit
   */

  int n_free = 0;
  auto m = I->FreeMember;
  while (m) {
    n_free++;
    m = I->Member[m].next;
  }
  if (n_free) {
    std::vector<int> list(n_free);
    auto l = list.data();
    auto m = I->FreeMember;
    while (m) {
      *(l++) = m;
      m = I->Member[m].next;
    }
    std::sort(list.begin(), list.end());
    auto NMember = int(I->Member.size()) - 1;
    while (n_free > 5000) { /* compact inactive members when possible */
      if (list[n_free - 1] == NMember) {
        NMember--;
        n_free--;
      } else
        break;
    }
    for (int a = 0; a < (n_free - 1); a++) {
      I->Member[list[a]].next = list[a + 1];
    }
    I->Member[list[n_free - 1]].next = 0;
    I->FreeMember = list[0];
    I->Member.resize(NMember + 1);
  }
}

void SelectorDeleteSeleAtIter(PyMOLGlobals* G, SelectorInfoIter_t it) {
  SelectorPurgeMembers(G, it->ID);
  G->SelectorMgr->Info.erase(it);
}

const char* SelectorGetNameFromIndex(PyMOLGlobals* G, SelectorID_t index) {
  auto I = G->SelectorMgr;
  for (int a = 1; a < I->Info.size(); a++) {
    if (I->Info[a].ID == index) {
      return I->Info[a].name.c_str();
    }
  }
  return nullptr;
}
int SelectorCountStates(PyMOLGlobals* G, int sele) {
  CSelector* I = G->Selector;
  int a;
  int result = 0;
  int n_frame;
  int at1;
  ObjectMolecule* last = nullptr;
  ObjectMolecule* obj;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];
    if (obj != last) {
      at1 = I->Table[a].atom;
      if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
        {
          n_frame = obj->getNFrame();
          if (result < n_frame) result = n_frame;
        }
        last = obj;
      }
    }
  }
  return (result);
}
int SelectorCheckIntersection(PyMOLGlobals* G, int sele1, int sele2) {
  CSelector* I = G->Selector;

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  for (int a = cNDummyAtoms; a < I->Table.size(); a++) {
    auto obj = I->Obj[I->Table[a].model];
    auto at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele1) &&
        SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele2))
      return 1;
  }
  return 0;
}
int SelectorCountAtoms(PyMOLGlobals* G, int sele, int state) {
  CSelector* I = G->Selector;
  int result = 0;

  SelectorUpdateTable(G, state, -1);
  for (int a = cNDummyAtoms; a < I->Table.size(); a++) {
    auto obj = I->Obj[I->Table[a].model];
    auto at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
      result++;
    }
  }
  return (result);
}
int* SelectorGetResidueVLA(PyMOLGlobals* G, SelectorID_t sele, int ca_only,
                           ObjectMolecule* exclude) {
  /* returns a VLA containing atom indices followed by residue integers
     (residue names packed as characters into integers)
     The indices are the first and last residue in the selection...
   */
  CSelector* I = G->Selector;
  int *result = nullptr, *r;
  AtomInfoType *ai1 = nullptr, *ai2;

  /* update the selector's table */
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  /* make room for model#, at#, rcode, per atom */
  result = VLAlloc(int, I->Table.size() * 3);

  r = result;
  PRINTFD(G, FB_Selector)
  " %s-DEBUG: entry, sele = %d\n", __func__, sele ENDFD;

  for (SeleAtomIterator iter(G, sele); iter.next();) {
    if (iter.obj == exclude) continue;

    ai2 = iter.getAtomInfo();

    if (ca_only) {
      if (!(ai2->flags & cAtomFlag_guide)) continue;
    } else if (ai1 && AtomInfoSameResidue(G, ai1, ai2)) {
      continue;
    }

    *(r++) = I->Table[iter.a].model;
    *(r++) = I->Table[iter.a].atom;

    const char* resn = LexStr(G, ai2->resn);
    *(r) = (resn[0] << (8 * 2));
    if (resn[0] && resn[1]) {
      *(r) |= (resn[1] << (8 * 1));
      *(r) |= (resn[2] << (8 * 0));
    }
    r++;

    ai1 = ai2;
  }
  if (result) {
    VLASize(result, int, (ov_size)(r - result));
  }
  PRINTFD(G, FB_Selector)
  " %s-DEBUG: exit, result = %p, size = %d\n", __func__, (void*)result,
      (unsigned int)VLAGetSize(result) ENDFD;

  return (result);
}
/* bad to make this non-static? */
/* static int *SelectorGetIndexVLA(PyMOLGlobals * G, int sele) */
int* SelectorGetIndexVLA(PyMOLGlobals* G, SelectorID_t sele) {
  return (SelectorGetIndexVLAImpl(G, G->Selector, sele));
}
static int* SelectorGetIndexVLAImpl(
    PyMOLGlobals* G, CSelector* I,
    SelectorID_t sele) { /* assumes updated tables */
  int a, c = 0;
  int* result = nullptr;
  ObjectMolecule* obj;
  int at1;

  result = VLAlloc(int, (I->Table.size() / 10) + 1);
  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];
    at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
      VLACheck(result, int, c);
      result[c++] = a;
    }
  }
  VLASize(result, int, c);
  return (result);
}
void SelectorUpdateObjectSele(PyMOLGlobals* G, ObjectMolecule* obj) {
  if (obj->Name[0]) {
    SelectorCreate(G, obj->Name, nullptr, obj, true, nullptr);
    /* create a selection with same name */
    if (SettingGetGlobal_b(G, cSetting_auto_classify_atoms)) {
      SelectorClassifyAtoms(G, 0, false, obj);

      // for file formats other than PDB
      if (obj->need_hetatm_classification) {
        for (auto ai = obj->AtomInfo.data(), ai_end = ai + obj->NAtom;
             ai != ai_end; ++ai) {
          if (!(ai->flags & cAtomFlag_polymer)) {
            ai->hetatm = true;
            ai->flags |= cAtomFlag_ignore;
          }
        }
        obj->need_hetatm_classification = false;
      }
    }
  }
}
void SelectorLogSele(PyMOLGlobals* G, const char* name) {
  CSelector* I = G->Selector;
  int a;
  std::string line, buf1;
  int cnt = -1;
  int first = 1;
  int append = 0;
  ObjectMolecule* obj;
  int at1;
  int sele;
  int logging;
  int robust;
  logging = SettingGetGlobal_i(G, cSetting_logging);
  robust = SettingGetGlobal_b(G, cSetting_robust_logs);
  if (logging) {
    sele = SelectorIndexByName(G, name);
    if (sele >= 0) {
      SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        obj = I->Obj[I->Table[a].model];
        at1 = I->Table[a].atom;
        if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
          if (cnt < 0) {
            if (first) {
              switch (logging) {
                case cPLog_pml:
                  line = pymol::string_format("_ cmd.select(\"%s\",\"(", name);
                  break;
                case cPLog_pym:
                  line = pymol::string_format("cmd.select(\"%s\",\"(", name);
                  break;
              }
              append = 0;
              cnt = 0;
              first = 0;
            } else {
              switch (logging) {
                case cPLog_pml:
                  line = pymol::string_format("_ cmd.select(\"%s\",\"(%s", name,
                                              name);
                  break;
                case cPLog_pym:
                  line = pymol::string_format("cmd.select(\"%s\",\"(%s", name,
                                              name);
                  break;
              }
              append = 1;
              cnt = 0;
            }
          }
          if (append) line += "|";
          if (robust)
            buf1 = ObjectMoleculeGetAtomSeleFast(obj, at1);
          else
            buf1 = pymol::string_format("%s`%d", obj->Name, at1 + 1);
          line += buf1;
          append = 1;
          cnt++;
          if (line.size() > (sizeof(OrthoLineType) / 2)) {
            line += ")\")\n";
            PLog(G, line.c_str(), cPLog_no_flush);
            cnt = -1;
          }
        }
      }
      if (cnt > 0) {
        line += ")\")\n";
        PLog(G, line.c_str(), cPLog_no_flush);
        PLogFlush(G);
      }
    }
  }
}
/**
 * This is the most heavily called routine in interactive PyMOL
 *
 * @param s    AtomInfoType.selEntry
 * @param sele selection index or 0 for "all"
 */
int SelectorIsMember(PyMOLGlobals* G, SelectorMemberOffset_t s,
                     SelectorID_t sele) {
  if (sele > 1) {
    const MemberType *mem, *member = G->SelectorMgr->Member.data();
    for (; s; s = mem->next) {
      mem = member + s;
      if (mem->selection == sele) return mem->tag;
    }
  } else if (!sele)
    return true; /* "all" is selection number 0, unordered */
  return false;
}
bool SelectorMoveMember(PyMOLGlobals* G, SelectorMemberOffset_t s,
                        SelectorID_t sele_old, SelectorID_t sele_new) {
  auto I = G->SelectorMgr;
  int result = false;
  while (s) {
    if (I->Member[s].selection == sele_old) {
      I->Member[s].selection = sele_new;
      result = true;
    }
    s = I->Member[s].next;
  }
  return result;
}
ObjectMolecule* SelectorGetFastSingleObjectMolecule(PyMOLGlobals* G,
                                                    SelectorID_t sele) {
  auto I = G->SelectorMgr;
  ObjectMolecule* result = nullptr;
  auto it = std::find_if(
      I->Info.begin(), I->Info.end(),
      [sele](const SelectionInfoRec& rec) { return rec.ID == sele; });
  if (it != I->Info.end()) {
    auto& info = *it;
    if (info.justOneObject()) {
      if (ExecutiveValidateObjectPtr(G, info.theOneObject, cObjectMolecule))
        result = info.theOneObject;
    } else {
      result = SelectorGetSingleObjectMolecule(
          G, sele); /* fallback onto slow approach */
    }
  }
  return (result);
}
ObjectMolecule* SelectorGetFastSingleAtomObjectIndex(PyMOLGlobals* G,
                                                     SelectorID_t sele,
                                                     int* index) {
  auto I = G->SelectorMgr;
  ObjectMolecule* result = nullptr;
  auto it = std::find_if(
      I->Info.begin(), I->Info.end(),
      [sele](const SelectionInfoRec& rec) { return rec.ID == sele; });
  if (it != I->Info.end()) {
    auto& info = *it;
    if (info.justOneAtom()) {
      ObjectMolecule* obj = info.theOneObject;
      int at = info.theOneAtom;
      if (ExecutiveValidateObjectPtr(G, obj, cObjectMolecule)) {
        if ((at < obj->NAtom) &&
            SelectorIsMember(G, obj->AtomInfo[at].selEntry, sele)) {
          *index = at;
          return obj;
        }
      }
    }
    /* fallback onto slow approach */
    {
      auto res = SelectorGetSingleAtomObjectIndex(G, sele);
      if (res) {
        std::tie(result, *index) = res.result();
      }
    }
  }
  return (result);
}

/*========================================================================
 * SelectorGetSingleObjectMolecule -- get a ptr to the molecule indiecated
 *    by the selection parameter
 * PARAMS
 *   (int) selection #
 * RETURNS
 *   (ptr) pts to the ObjectMolecule or nullptr if not found
 */
ObjectMolecule* SelectorGetSingleObjectMolecule(PyMOLGlobals* G,
                                                SelectorID_t sele) {
  /* slow way */

  int a;
  ObjectMolecule* result = nullptr;
  ObjectMolecule* obj;
  CSelector* I = G->Selector;
  int at1;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];
    at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
      if (result) {
        if (result != obj) {
          result = nullptr;
          break;
        }
      } else {
        result = obj;
      }
    }
  }
  return (result);
}
ObjectMolecule* SelectorGetFirstObjectMolecule(PyMOLGlobals* G,
                                               SelectorID_t sele) {
  /* slow way */

  int a;
  ObjectMolecule* result = nullptr;
  ObjectMolecule* obj;
  CSelector* I = G->Selector;
  int at1;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];
    at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
      result = obj;
      break;
    }
  }
  return (result);
}
ObjectMolecule** SelectorGetObjectMoleculeVLA(PyMOLGlobals* G,
                                              SelectorID_t sele) {
  int a;
  ObjectMolecule* last = nullptr;
  ObjectMolecule *obj, **result = nullptr;
  CSelector* I = G->Selector;
  int at1;
  int n = 0;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  result = VLAlloc(ObjectMolecule*, 10);
  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];
    at1 = I->Table[a].atom;
    if (SelectorIsMember(G, obj->AtomInfo[at1].selEntry, sele)) {
      if (obj != last) {
        VLACheck(result, ObjectMolecule*, n);
        result[n] = obj;
        last = obj;
        n++;
      }
    }
  }
  VLASize(result, ObjectMolecule*, n);
  return (result);
}
pymol::Result<std::pair<ObjectMolecule*, int>> SelectorGetSingleAtomObjectIndex(
    PyMOLGlobals* G, SelectorID_t sele) {
  /* slow way */

  bool found_it = false;
  void* iterator = nullptr;
  ObjectMolecule* obj = nullptr;
  std::pair<ObjectMolecule*, int> result;

  while (ExecutiveIterateObjectMolecule(G, &obj, &iterator)) {
    const AtomInfoType* ai = obj->AtomInfo.data();
    for (int a = 0; a < obj->NAtom; a++) {
      int s = (ai++)->selEntry;
      if (SelectorIsMember(G, s, sele)) {
        if (found_it) {
          return pymol::make_error(
              "More than one atom found"); /* ADD'L EXIT POINT */
        } else {
          result = std::make_pair(obj, a);
          found_it = true;
        }
      }
    }
  }
  if (found_it) {
    return result;
  } else {
    return pymol::make_error("Not found");
  }
}
pymol::Result<std::array<float, 3>> SelectorGetSingleAtomVertex(PyMOLGlobals* G,
                                                                int sele,
                                                                int state) {
  auto atom_index_result = SelectorGetSingleAtomObjectIndex(G, sele);
  p_return_if_error(atom_index_result);
  {
    auto obj_idx = atom_index_result.result();
    std::array<float, 3> v;
    auto found_it = ObjectMoleculeGetAtomTxfVertex(obj_idx.first, state,
                                                   obj_idx.second, v.data());
    if (found_it) {
      return v;
    } else {
      return pymol::make_error("Invalid Atom");
    }
  }
}
void SelectorDeletePrefixSet(PyMOLGlobals* G, const char* pref) {
  auto I = G->SelectorMgr;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

  while (1) {
    auto it = SelectGetInfoIter(G, pref, strlen(pref), ignore_case);
    if (it == I->Info.end()) {
      break;
    }

    // important to use a copy, otherwise you'll delete all objects
    auto name_copy = it->name;
    ExecutiveDelete(G, name_copy.c_str());
  }
}
#define MAX_DEPTH 1000

int SelectorCheckNeighbors(PyMOLGlobals* G, int maxDist, ObjectMolecule* obj,
                           int at1, int at2, int* zero, int* scratch) {
  int stkDepth = 0;
  int si = 0;
  int stk[MAX_DEPTH];

  zero[at1] = 0;
  scratch[si++] = at1;
  stk[stkDepth] = at1;
  stkDepth++;

  while (stkDepth) { /* this will explore a tree */
    stkDepth--;
    auto const a = stk[stkDepth];
    auto const dist = zero[a] + 1;

    for (auto const& neighbor : AtomNeighbors(obj, a)) {
      auto const a1 = neighbor.atm;
      if (a1 == at2) {
        while (si--) {
          zero[scratch[si]] = 0;
        }
        /* EXIT POINT 1 */
        return 1;
      }
      if ((!zero[a1]) && (stkDepth < MAX_DEPTH) && (dist < maxDist)) {
        zero[a1] = dist;
        scratch[si++] = a1;
        stk[stkDepth] = a1;
        stkDepth++;
      }
    }
  }
  while (si--) {
    zero[scratch[si]] = 0;
  }
  /* EXIT POINT 2 */
  return 0;
}
static int SelectorWalkTree(PyMOLGlobals* G, int* atom, int* comp, int* toDo,
                            int** stk, int stkDepth, ObjectMolecule* obj,
                            int sele1, int sele2, int sele3, int sele4) {
  int c = 0;

  while (stkDepth) { /* this will explore a tree, stopping at protected atoms */
    stkDepth--;
    auto const a = (*stk)[stkDepth];
    toDo[a] = 0;
    auto const* ai = obj->AtomInfo.data() + a;
    auto const s = ai->selEntry;
    bool const seleFlag =
        SelectorIsMember(G, s, sele1) || SelectorIsMember(G, s, sele2) ||
        SelectorIsMember(G, s, sele3) || SelectorIsMember(G, s, sele4);
    if (!seleFlag) {
      if (!(ai->protekted ==
            cAtomProtected_explicit)) { /* if not explicitly protected... */
        atom[a] = 1;                    /* mark this atom into the selection */
        comp[a] = 1;
      }

      // add neighbors onto the stack
      for (auto const& neighbor : AtomNeighbors(obj, a)) {
        if (toDo[neighbor.atm]) {
          VLACheck((*stk), int, stkDepth);
          (*stk)[stkDepth] = neighbor.atm;
          stkDepth++;
        }
      }
      c++;
    }
  }
  return (c);
}
static int SelectorWalkTreeDepth(PyMOLGlobals* G, int* atom, int* comp,
                                 int* toDo, int** stk, int stkDepth,
                                 ObjectMolecule* obj, int sele1, int sele2,
                                 int sele3, int sele4, int** extraStk,
                                 WalkDepthRec* wd) {
  int s;
  int c = 0;
  int seleFlag;
  int depth;
  AtomInfoType* ai;

  wd->depth1 = -1;
  wd->depth2 = -1;
  wd->depth3 = -1;
  wd->depth4 = -1;
  VLACheck(*extraStk, int, stkDepth);
  UtilZeroMem(*extraStk, sizeof(int) * stkDepth);

  while (stkDepth) { /* this will explore a tree, stopping at protected atoms */
    stkDepth--;
    auto const a = (*stk)[stkDepth];
    depth = ((*extraStk)[stkDepth] + 1);
    seleFlag = false;
    ai = obj->AtomInfo + a;
    s = ai->selEntry;

    /* record how many cycles it take to reach each & any picked atoms */

    seleFlag = false;
    if (SelectorIsMember(G, s, sele1)) {
      if (((wd->depth1 < 0) || (wd->depth1 > depth))) {
        wd->depth1 = depth;
      }
      seleFlag = true;
    }
    if (SelectorIsMember(G, s, sele2)) {
      if (((wd->depth2 < 0) || (wd->depth2 > depth))) {
        wd->depth2 = depth;
      }
      seleFlag = true;
    }
    if (SelectorIsMember(G, s, sele3)) {
      if (((wd->depth3 < 0) || (wd->depth3 > depth))) {
        wd->depth3 = depth;
      }
      seleFlag = true;
    }
    if (SelectorIsMember(G, s, sele4)) {
      if (((wd->depth4 < 0) || (wd->depth4 > depth))) {
        wd->depth4 = depth;
      }
      seleFlag = true;
    }

    if (!seleFlag) {
      toDo[a] = 0;
      if (!(ai->protekted ==
            cAtomProtected_explicit)) { /* if not explicitly protected... */
        atom[a] = 1;                    /* mark this atom into the selection */
        comp[a] = 1;
      }

      /* add neighbors onto the stack */
      for (auto const& neighbor : AtomNeighbors(obj, a)) {
        if (toDo[neighbor.atm]) {
          VLACheck((*stk), int, stkDepth);
          (*stk)[stkDepth] = neighbor.atm;
          VLACheck((*extraStk), int, stkDepth);
          (*extraStk)[stkDepth] = depth;
          stkDepth++;
        }
      }
      c++;
    }
  }
  return (c);
}
int SelectorIsAtomBondedToSele(PyMOLGlobals* G, ObjectMolecule* obj,
                               SelectorID_t sele1atom,  //
                               SelectorID_t sele2) {
  auto const atm = ObjectMoleculeGetAtomIndex(obj, sele1atom);
  return (atm >= 0) && ObjectMoleculeIsAtomBondedToSele(obj, atm, sele2);
}

static void update_min_walk_depth(WalkDepthRec* minWD, int frag,
                                  WalkDepthRec* wd, int sele1, int sele2,
                                  int sele3, int sele4) {
  /* first, does this fragment even qualify ? */
  int qualifies = true;
  int cnt = 0;
  wd->sum = 0;
  if (sele1 >= 0) {
    if (wd->depth1 < 0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth1;
      cnt++;
    }
  }
  if (sele2 >= 0) {
    if (wd->depth2 < 0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth2;
      cnt++;
    }
  }
  if (sele3 >= 0) {
    if (wd->depth3 < 0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth3;
      cnt++;
    }
  }
  if (sele4 >= 0) {
    if (wd->depth4 < 0) {
      qualifies = false;
    } else {
      wd->sum += wd->depth4;
      cnt++;
    }
  }
  if (qualifies && (cnt > 1)) {
    /* is it better than the current min? */

    if ((!minWD->frag) || (wd->sum < minWD->sum)) {
      (*minWD) = (*wd);
      minWD->frag = frag;
    }
  }
}
int SelectorSubdivide(PyMOLGlobals* G,     //
                      const char* pref,    //
                      SelectorID_t sele1,  //
                      SelectorID_t sele2,  //
                      SelectorID_t sele3,  //
                      SelectorID_t sele4,  //
                      const char* fragPref, const char* compName,
                      int* bondMode) {
  CSelector* I = G->Selector;
  int a0 = 0, a1 = 0, a2;
  int* atom = nullptr;
  int* toDo = nullptr;
  int* comp = nullptr;
  int* pkset = nullptr;
  int set_cnt = 0;
  int nFrag = 0;
  int stkDepth;
  int c;
  int cycFlag = false;
  std::string name, link_sele;
  ObjectMolecule *obj1 = nullptr, *obj2 = nullptr, *obj3 = nullptr,
                 *obj4 = nullptr;
  int index1 = 0, index2 = 0, index3 = 0, index4 = 0;

  /* this is seriously getting out of hand -- need to switch over to arrays soon
   */

  int *atom1_base = nullptr, *atom2_base = nullptr, *atom3_base = nullptr,
      *atom4_base = nullptr;
  int *toDo1_base = nullptr, *toDo2_base = nullptr, *toDo3_base = nullptr,
      *toDo4_base = nullptr;
  int *comp1_base = nullptr, *comp2_base = nullptr, *comp3_base = nullptr,
      *comp4_base = nullptr;
  int *pkset1_base = nullptr, *pkset2_base = nullptr, *pkset3_base = nullptr,
      *pkset4_base = nullptr;

  PRINTFD(G, FB_Selector)
  " SelectorSubdivideObject: entered...\n" ENDFD;
  SelectorDeletePrefixSet(G, pref);
  SelectorDeletePrefixSet(G, fragPref);
  ExecutiveDelete(G, cEditorLink);
  ExecutiveDelete(G, cEditorSet);
  /* delete any existing matches */

  obj1 = SelectorGetFastSingleAtomObjectIndex(G, sele1, &index1);
  obj2 = SelectorGetFastSingleAtomObjectIndex(G, sele2, &index2);
  obj3 = SelectorGetFastSingleAtomObjectIndex(G, sele3, &index3);
  obj4 = SelectorGetFastSingleAtomObjectIndex(G, sele4, &index4);

  if (obj1 || obj2 || obj3 || obj4) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

    comp = pymol::calloc<int>(I->Table.size());
    atom = pymol::malloc<int>(I->Table.size());
    toDo = pymol::malloc<int>(I->Table.size());
    pkset = pymol::calloc<int>(I->Table.size());

    /* NOTE: SeleBase only safe with cSelectorUpdateTableAllStates!  */

    if (obj1) {
      atom1_base = atom + obj1->SeleBase;
      toDo1_base = toDo + obj1->SeleBase;
      comp1_base = comp + obj1->SeleBase;
      pkset1_base = pkset + obj1->SeleBase;
    }

    if (obj2) {
      atom2_base = atom + obj2->SeleBase;
      toDo2_base = toDo + obj2->SeleBase;
      comp2_base = comp + obj2->SeleBase;
      pkset2_base = pkset + obj2->SeleBase;
    }

    if (obj3) {
      atom3_base = atom + obj3->SeleBase;
      toDo3_base = toDo + obj3->SeleBase;
      comp3_base = comp + obj3->SeleBase;
      pkset3_base = pkset + obj3->SeleBase;
    }

    if (obj4) {
      atom4_base = atom + obj4->SeleBase;
      toDo4_base = toDo + obj4->SeleBase;
      comp4_base = comp + obj4->SeleBase;
      pkset4_base = pkset + obj4->SeleBase;
    }

    auto stk = pymol::vla<int>(100);

    {
      int a;
      int* p1;
      p1 = toDo;
      for (a = 0; a < I->Table.size(); a++) *(p1++) = true;
    }

    if (*bondMode) {
      /* verify bond mode, or clear the flag */

      *bondMode = false;

      if ((sele1 >= 0) && (sele2 >= 0) && (sele3 < 0) && (sele4 < 0) &&
          (obj1 == obj2)) {
        /* two selections only, in same object... */

        a0 = index1;
        a1 = index2;

        if ((a0 >= 0) && (a1 >= 0)) {
          for (auto const& neighbor : AtomNeighbors(obj1, a0)) {
            if (neighbor.atm == a1) {
              *bondMode = true;
              break;
            }
          }
        }
      }
    }

    /* ===== BOND MODE ===== (sele0 and sele1 only) */

    if (*bondMode) {
      if (obj1 == obj2) { /* just to be safe */

        pkset1_base[a0] = 1;
        pkset1_base[a1] = 1;
        SelectorEmbedSelection(G, pkset, cEditorBond, nullptr, false, -1);

        a0 = index1;
        if (a0 >= 0) {
          stkDepth = 0;
          for (auto const& neighbor : AtomNeighbors(obj1, a0)) {
            if (toDo1_base[neighbor.atm]) {
              *(stk.check(stkDepth++)) = neighbor.atm;
            }
          }
          UtilZeroMem(atom, sizeof(int) * I->Table.size());
          atom1_base[a0] = 1; /* create selection for this atom alone as
                                 fragment base atom */
          comp1_base[a0] = 1;
          name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
          SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
          c = SelectorWalkTree(G, atom1_base, comp1_base, toDo1_base, &stk,
                               stkDepth, obj1, sele1, sele2, -1, -1) +
              1;
          name = pymol::string_format("%s%1d", pref, nFrag + 1);

          /* check for cyclic situation */
          cycFlag = false;
          a2 = index2;
          if (a2 >= 0) {
            for (auto const& neighbor : AtomNeighbors(obj1, a2)) {
              auto const a1 = neighbor.atm;
              if (a1 != a0 && !toDo1_base[a1]) {
                cycFlag = true; /* we have a cycle... */
                break;
              }
            }
          }
          if (cycFlag) { /* cyclic situation is a bit complex... */

            a0 = index2;
            if (a0 >= 0) {
              stkDepth = 0;
              for (auto const& neighbor : AtomNeighbors(obj1, a0)) {
                if (toDo1_base[neighbor.atm]) {
                  *(stk.check(stkDepth++)) = neighbor.atm;
                }
              }
              atom1_base[a0] = 1;
              comp1_base[a0] = 1;
              c = SelectorWalkTree(G, atom1_base, comp1_base, toDo1_base, &stk,
                                   stkDepth, obj1, sele1, sele2, -1, -1) +
                  1;
            }
          }
          SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
          nFrag++;
        }

        if (!cycFlag) {
          a0 = index2;
          if (a0 >= 0) {
            stkDepth = 0;
            for (auto const& neighbor : AtomNeighbors(obj1, a0)) {
              if (toDo1_base[neighbor.atm]) {
                *(stk.check(stkDepth++)) = neighbor.atm;
              }
            }

            UtilZeroMem(atom, sizeof(int) * I->Table.size());
            atom1_base[a0] = 1; /* create selection for this atom alone as
                                   fragment base atom */
            comp1_base[a0] = 1;
            name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
            SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
            c = SelectorWalkTree(G, atom1_base, comp1_base, toDo1_base, &stk,
                                 stkDepth, obj1, sele1, sele2, -1, -1) +
                1;
            name = pymol::string_format("%s%1d", pref, nFrag + 1);
            SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
            nFrag++;
          }
        }
      }
    } else {
      /* ===== WALK MODE ===== (any combination of sele0, sele1, sele2, sele3 */

      int* extraStk = VLAlloc(int, 50);
      WalkDepthRec curWalk, minWalk;
      minWalk.sum = 0;
      minWalk.frag = 0;

      if (obj1) {
        a0 = index1;
        if (a0 >= 0) {
          pkset1_base[a0] = 1;
          set_cnt++;
          comp1_base[a0] = 1;
          stkDepth = 0;
          for (auto const& neighbor : AtomNeighbors(obj1, a0)) {
            auto const a1 = neighbor.atm;
            if (toDo1_base[a1]) {
              stkDepth = 1;
              stk[0] = a1;
              UtilZeroMem(atom, sizeof(int) * I->Table.size());
              atom1_base[a1] = 1; /* create selection for this atom alone as
                                     fragment base atom */
              comp1_base[a1] = 1;
              name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
              SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
              atom1_base[a1] = 0;
              c = SelectorWalkTreeDepth(G, atom1_base, comp1_base, toDo1_base,
                                        &stk, stkDepth, obj1, sele1, sele2,
                                        sele3, sele4, &extraStk, &curWalk);
              if (c) {
                nFrag++;
                name = pymol::string_format("%s%1d", pref, nFrag);
                SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
                update_min_walk_depth(&minWalk, nFrag, &curWalk, sele1, sele2,
                                      sele3, sele4);
              }
            }
          }
        }
      }

      if (obj2) {
        a0 = index2;
        if (a0 >= 0) {
          pkset2_base[a0] = 1;
          set_cnt++;
          comp2_base[a0] = 1;
          stkDepth = 0;
          for (auto const& neighbor : AtomNeighbors(obj2, a0)) {
            auto const a1 = neighbor.atm;
            if (toDo2_base[a1]) {
              stkDepth = 1;
              stk[0] = a1;
              UtilZeroMem(atom, sizeof(int) * I->Table.size());
              atom2_base[a1] = 1; /* create selection for this atom alone as
                                     fragment base atom */
              comp2_base[a1] = 1;
              name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
              SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
              atom2_base[a1] = 0;
              c = SelectorWalkTreeDepth(G, atom2_base, comp2_base, toDo2_base,
                                        &stk, stkDepth, obj2, sele1, sele2,
                                        sele3, sele4, &extraStk, &curWalk);
              if (c) {
                nFrag++;
                name = pymol::string_format("%s%1d", pref, nFrag);
                SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
                update_min_walk_depth(&minWalk, nFrag, &curWalk, sele1, sele2,
                                      sele3, sele4);
              }
            }
          }
        }
      }

      if (obj3) {
        a0 = index3;
        if (a0 >= 0) {
          pkset3_base[a0] = 1;
          set_cnt++;
          comp3_base[a0] = 1;
          stkDepth = 0;
          for (auto const& neighbor : AtomNeighbors(obj3, a0)) {
            auto const a1 = neighbor.atm;
            if (toDo3_base[a1]) {
              stkDepth = 1;
              stk[0] = a1;
              UtilZeroMem(atom, sizeof(int) * I->Table.size());
              atom3_base[a1] = 1; /* create selection for this atom alone as
                                     fragment base atom */
              comp3_base[a1] = 1;
              name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
              SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
              atom3_base[a1] = 0;
              c = SelectorWalkTreeDepth(G, atom3_base, comp3_base, toDo3_base,
                                        &stk, stkDepth, obj3, sele1, sele2,
                                        sele3, sele4, &extraStk, &curWalk);
              if (c) {
                nFrag++;
                name = pymol::string_format("%s%1d", pref, nFrag);
                SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
                update_min_walk_depth(&minWalk, nFrag, &curWalk, sele1, sele2,
                                      sele3, sele4);
              }
            }
          }
        }
      }

      if (obj4) {
        a0 = index4;
        if (a0 >= 0) {
          pkset4_base[a0] = 1;
          set_cnt++;
          comp4_base[a0] = 1;
          stkDepth = 0;
          for (auto const& neighbor : AtomNeighbors(obj4, a0)) {
            auto const a1 = neighbor.atm;
            if (toDo4_base[a1]) {
              stkDepth = 1;
              stk[0] = a1;
              UtilZeroMem(atom, sizeof(int) * I->Table.size());
              atom4_base[a1] = 1; /* create selection for this atom alone as
                                     fragment base atom */
              comp4_base[a1] = 1;
              name = pymol::string_format("%s%1d", fragPref, nFrag + 1);
              SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
              atom4_base[a1] = 0;
              c = SelectorWalkTreeDepth(G, atom4_base, comp4_base, toDo4_base,
                                        &stk, stkDepth, obj4, sele1, sele2,
                                        sele3, sele4, &extraStk, &curWalk);
              if (c) {
                nFrag++;
                name = pymol::string_format("%s%1d", pref, nFrag);
                SelectorEmbedSelection(G, atom, name, nullptr, false, -1);
                update_min_walk_depth(&minWalk, nFrag, &curWalk, sele1, sele2,
                                      sele3, sele4);
              }
            }
          }
        }
      }

      if (minWalk.frag) { /* create the linking selection if one exists */
        link_sele = pymol::string_format("%s%d|?pk1|?pk2|?pk3|?pk4", pref,
                                         minWalk.frag);
      }
      VLAFreeP(extraStk);
    }

    if (set_cnt > 1) {
      SelectorEmbedSelection(G, pkset, cEditorSet, nullptr, false, -1);
    }

    if (nFrag) {
      SelectorEmbedSelection(G, comp, compName, nullptr, false, -1);
    }

    if (!link_sele.empty())
      SelectorCreate(G, cEditorLink, link_sele.c_str(), nullptr, true, nullptr);

    FreeP(toDo);
    FreeP(atom);
    FreeP(comp);
    FreeP(pkset);
    SelectorClean(G);
  }
  PRINTFD(G, FB_Selector)
  " SelectorSubdivideObject: leaving...nFrag %d\n", nFrag ENDFD;

  return (nFrag);
}
int SelectorGetSeleNCSet(PyMOLGlobals* G, SelectorID_t sele) {
  CSelector* I = G->Selector;

  int a, s, at = 0;
  ObjectMolecule *obj, *last_obj = nullptr;
  int result = 0;

  if ((obj = SelectorGetFastSingleAtomObjectIndex(G, sele, &at))) {
    int a = obj->NCSet;
    CoordSet* cs;
    int idx;

    while (a--) {
      cs = obj->CSet[a];
      idx = cs->atmToIdx(at);
      if (idx >= 0) {
        result = a + 1;
        break;
      }
    }
  } else {
    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      obj = I->Obj[I->Table[a].model];
      if (obj != last_obj) {
        at = I->Table[a].atom;
        s = obj->AtomInfo[at].selEntry;
        if (SelectorIsMember(G, s, sele)) {
          if (result < obj->NCSet) {
            result = obj->NCSet;
            last_obj = obj;
          }
        }
      }
    }
  }
  return (result);
}

std::vector<int> SelectorGetInterstateVector(
    PyMOLGlobals* G, int sele1, int state1, int sele2, int state2,
    float cutoff) { /* Assumes valid tables */
  const size_t table_size = G->Selector->Table.size();
  auto coords_flat = std::vector<float>(3 * table_size);
  auto* coords = pymol::reshape<3>(coords_flat.data());

  // number of atoms in `sele1`
  int n1 = 0;

  // mask on selected atoms in `sele1`
  auto flags = std::vector<MapFlag_t>(table_size);

  // copy coordinates of selection 1
  for (SeleCoordIterator iter(G, sele1, state1, false); iter.next();) {
    copy3(iter.getCoord(), coords[iter.a]);
    flags[iter.a] = true;
    n1++;
  }

  if (n1 == 0) {
    // no atoms in `sele1`
    return {};
  }

  std::unique_ptr<MapType> map(new MapType(G, -cutoff, pymol::flatten(coords),
                                           table_size, nullptr, flags.data()));

  if (!map) {
    PRINTFB(G, FB_Selector, FB_Errors)
    " Selector-Error: unexpected map allocation failure\n" ENDFB(G);
    return {};
  }

  std::vector<int> out;

  for (SeleCoordIterator iter(G, sele2, state2, false); iter.next();) {
    const float* v2 = iter.getCoord();
    for (const auto a1 : MapEIter(*map, v2)) {
      if (within3f(coords[a1], v2, cutoff)) {
        out.push_back(a1);
        out.push_back(iter.a);
      }
    }
  }

  return out;
}
int SelectorAssignAtomTypes(PyMOLGlobals* G, int sele, int state, int quiet,
                            int format) {
  if (format != 1) {
    PRINTFB(G, FB_Selector, FB_Errors)
    " Error: assign_atom_types only supports format='mol2'\n" ENDFB(G);
    return 0;
  }

  SelectorUpdateTable(G, state, -1);

  ObjectMolecule* obj = nullptr;

  for (SeleAtomIterator iter(G, sele); iter.next();) {
    if (obj != iter.obj) {
      obj = iter.obj;
      ObjectMoleculeVerifyChemistry(obj, state);
    }

    LexAssign(G, iter.getAtomInfo()->textType, getMOL2Type(obj, iter.getAtm()));
  }
  return 1;
}
int SelectorSetName(PyMOLGlobals* G, const char* new_name,
                    const char* old_name) {
  auto I = G->SelectorMgr;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

  auto it = SelectGetInfoIter(G, old_name, 1, ignore_case);
  if (it != I->Info.end()) {
    it->name = new_name;
    return true;
  } else {
    return false;
  }
}

/*========================================================================
 * SelectorIndexByName -- fetch the global selector's ID for sname
 * PARAMS
 *  (string) sname, object name
 * RETURNS
 *   (int) index #, or -1 if not found
 */
SelectorID_t SelectorIndexByName(PyMOLGlobals* G, const char* sname,
                                 int ignore_case) {
  auto I = G->SelectorMgr;

  if (sname) {
    if (ignore_case < 0)
      ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

    while ((sname[0] == '%') || (sname[0] == '?')) sname++;

    auto it = SelectGetInfoIter(G, sname, 1, ignore_case);
    if (it == I->Info.end()) return cSelectionInvalid;

    if (sname[0] != '_') { /* don't do checking on internal selections */
      const char* best;
      best = ExecutiveFindBestNameMatch(
          G, sname); /* suppress spurious matches
                        of selections with non-selections */
      if (best != sname && best != it->name) return cSelectionInvalid;
    }

    return it->ID;
  }
  return cSelectionInvalid;
}
static void SelectorPurgeMembers(PyMOLGlobals* G, SelectorID_t sele) {
  auto I = G->SelectorMgr;
  void* iterator = nullptr;
  ObjectMolecule* obj = nullptr;
  short changed = 0;

  if (!I->Member.empty()) {
    while (ExecutiveIterateObjectMolecule(G, &obj, &iterator)) {
      if (obj->type == cObjectMolecule) {
        AtomInfoType* ai = obj->AtomInfo.data();
        int a, n_atom = obj->NAtom;
        for (a = 0; a < n_atom; a++) {
          int s = (ai++)->selEntry;
          int l = -1;
          while (s) {
            auto& i_member_s = I->Member[s];
            int nxt = i_member_s.next;
            if (i_member_s.selection == sele) {
              if (l > 0)
                I->Member[l].next = i_member_s.next;
              else
                ai[-1].selEntry = i_member_s.next;
              changed = 1;
              i_member_s.next = I->FreeMember;
              I->FreeMember = s;
            }
            l = s;
            s = nxt;
          }
        }
      }
    }
  }
  if (changed) {
    // not sure if this is needed since its in SelectorClean()
    ExecutiveInvalidateSelectionIndicatorsCGO(G);
  }
}
int SelectorPurgeObjectMembers(PyMOLGlobals* G, ObjectMolecule* obj) {
  bool changed = false;

  auto I = G->SelectorMgr;
  if (!I->Member.empty()) {
    for (int a = 0; a < obj->NAtom; a++) {
      auto s = obj->AtomInfo[a].selEntry;
      while (s) {
        auto nxt = I->Member[s].next;
        I->Member[s].next = I->FreeMember;
        I->FreeMember = s;
        s = nxt;
      }
      obj->AtomInfo[a].selEntry = 0;
      changed = true;
    }
  }
  if (changed) {
    // not sure if this is needed since its in SelectorClean()
    ExecutiveInvalidateSelectionIndicatorsCGO(G);
  }

  return 1;
}
void SelectorDelete(PyMOLGlobals* G, const char* sele)

/* should (only) be called by Executive or by Selector, unless the

   named selection has never been registered with the executive

   (i.e., temporary on-the-fly selection) */
{
  auto& Info = G->SelectorMgr->Info;
  auto it = SelectGetInfoIter(G, sele, 999,
                              SettingGetGlobal_b(G, cSetting_ignore_case));

  // Does it exist?
  if (it == Info.end()) return;

  // Never delete the "all" selection
  if (it->ID == cSelectionAll) return;

  assert(!SelectorIsTmp(sele) ||
         sele == pymol::string_format("%s%d", cSelectorTmpPrefix, it->ID));

  // get rid of existing selection
  SelectorDeleteSeleAtIter(G, it);
}
/**
 * If `input` is already a name of an object or a valid position keyword
 * (center, origin, all, ...), then simply copy it to `store`. Otherwise
 * process the selection expression and create a temporary named selection.
 *
 * Unlike `SelectorGetTmp`, this will accept names of maps, groups, etc.
 */
int SelectorGetTmp2(PyMOLGlobals* G, const char* input, char* store,
                    bool quiet) {
  auto res = SelectorGetTmp2Result(G, input, store, quiet);
  if (res) {
    return res.result();
  }
  PRINTFB(G, FB_Selector, FB_Errors)
  " Selector-Error: %s\n", res.error().what().c_str() ENDFB(G);
  return -1;
}

pymol::Result<int> SelectorGetTmp2Result(PyMOLGlobals* G, const char* input,
                                         char* store, bool quiet) {
  /* ASSUMES that store is at least as big as an OrthoLineType */
  auto I = G->SelectorMgr;
  PRINTFD(G, FB_Selector)
  " %s-Debug: entered with \"%s\".\n", __func__, input ENDFD;

  store[0] = 0;

  /* skip trivial cases */

  if (input[0] && !((input[0] == '\'') && (input[1] == '\'') && (!input[2]))) {
    /* OKAY, this routine is in flux.  Eventually this routine will...

       (1) fully parse the input recognizing selection keywords, nested
       parens, quotes, escaped strings, etc.

       (2) replace selection blocks with temporary selection names

       (3) return a space-separated list of names for processing

       However, right now, this routine simply handles two cases.

       A. where the input is a selection, in which case store is set to a
       temporary selection name

       B. where the input is simply a list of space-separated name patterns,
       in which case store is simply passed along as a copy of the input

     */

    // make selection if "input" doesn't fit into "store"
    int is_selection = strlen(input) >= OrthoLineLength;

    // can't pass through temp selections (current SelectorFreeTmp limitation)
    if (!is_selection) {
      is_selection = SelectorIsTmp(input);
    }

    const char* p = input;
    OrthoLineType word;

    if (!is_selection)
      while (*p) {
        /* copy first word/token of p into "word", remainder of string in p */
        p = ParseWord(word, p, sizeof(OrthoLineType));
        /* see a paren? then this must be a selection */

        if (word[0] == '(') {
          is_selection = true;
          break;
        }

        if (strchr(word, '/')) {
          is_selection = true;
          break;
        }

        /* encounterd a selection keyword? then this must be a selection */

        {
          auto it = I->Key.find(word);
          if (it != I->Key.end()) {
            if (it->second != SELE_ALLz && it->second != SELE_ORIz &&
                it->second != SELE_CENz) {
              is_selection = true;
              break;
            }
          }
        }

        if (!ExecutiveValidName(G, word)) { /* don't recognize the name? */
          if (!ExecutiveValidNamePattern(
                  G, word)) {    /* don't recognize this as a pattern? */
            is_selection = true; /* must be a selection */
            break;
          }
        }
      }
    if (is_selection) { /* incur the computational expense of
                           parsing the input as an atom selection */
      SelectorGetUniqueTmpName(G, store);
      auto res = SelectorCreate(G, store, input, nullptr, quiet, nullptr);
      if (!res) {
        store[0] = 0;
      }
      return res;
    } else { /* otherwise, just parse the input as a space-separated list of
                names */
      /* not a selection */
      strcpy(store, input);
    }
  }
  return 0;
}
/**
 * Like SelectorGetTmp2, but doesn't accept names from any non-molecular
 * entities like groups or map objects (those will be processed as selection
 * expressions).
 */
int SelectorGetTmp(PyMOLGlobals* G, const char* input, char* store,
                   bool quiet) {
  auto res = SelectorGetTmpResult(G, input, store, quiet);
  if (res) {
    return res.result();
  }
  PRINTFB(G, FB_Selector, FB_Errors)
  " Selector-Error: %s\n", res.error().what().c_str() ENDFB(G);
  return -1;
}

pymol::Result<int> SelectorGetTmpResult(PyMOLGlobals* G, const char* input,
                                        char* store, bool quiet) {
  store[0] = 0;

  // trivial (but valid) case: empty selection string
  if (!input[0]) return 0;

  // if object molecule or named selection, then don't create a temp selection
  if (ExecutiveIsMoleculeOrSelection(G, input) && !SelectorIsTmp(input)) {
    strcpy(store, input);
    return 0;
  }

  // evaluate expression and create a temp selection
  SelectorGetUniqueTmpName(G, store);
  auto res = SelectorCreate(G, store, input, nullptr, quiet, nullptr);

  if (!res) {
    store[0] = 0;
  }

  return res;
}
void SelectorFreeTmp(PyMOLGlobals* G,
                     const char* name) { /* remove temporary selections */
  if (name && SelectorIsTmp(name)) {
    ExecutiveDelete(G, name);
  }
}
int SelectorEmbedSelection(PyMOLGlobals* G, const int* atom,
                           pymol::zstring_view name, ObjectMolecule* obj,
                           int no_dummies, int exec_managed) {
  /* either atom or obj should be nullptr, not both and not neither */

  CSelector* I = G->Selector;
  auto IM = I->mgr;
  int tag;
  int newFlag = true;
  int a, sele;
  int c = 0;
  int start = 0;
  int singleAtomFlag = true;
  int singleObjectFlag = true;
  ObjectMolecule *singleObject = nullptr, *selObj;
  int singleAtom = -1;
  int index;
  AtomInfoType* ai;

  if (exec_managed < 0) {
    if (atom) /* automatic behavior: manage selections defined via atom masks */
      exec_managed = true;
    else
      exec_managed = false;
  }

  // already exist?
  auto it = SelectGetInfoIter(G, name.c_str(), 999,
                              SettingGetGlobal_b(G, cSetting_ignore_case));
  if (it != IM->Info.end()) {
    assert(!SelectorIsTmp(name));

    // don't allow redefinition of "all"
    if (it->ID == cSelectionAll) return 0;

    // get rid of existing selection
    SelectorDeleteSeleAtIter(G, it);
    newFlag = false;
  }

  sele = IM->NSelection++;
  IM->Info.emplace_back(SelectionInfoRec(sele, name.c_str()));

  assert(!SelectorIsTmp(name) ||
         name == pymol::string_format("%s%d", cSelectorTmpPrefix,
                                      IM->Info.back().ID));

  if (no_dummies) {
    start = 0;
  } else {
    start = cNDummyAtoms;
  }
  for (a = start; a < I->Table.size(); a++) {
    tag = false;
    /* set tag based on passed in atom list or on global atom table */
    if (atom) {
      if (atom[a]) tag = atom[a];
    } else {
      if (I->Obj[I->Table[a].model] == obj) tag = 1;
    }
    if (tag) {
      /* if this this atom is tagged, grab its object
       * index, and info record */
      selObj = I->Obj[I->Table[a].model];
      index = I->Table[a].atom;
      ai = selObj->AtomInfo + index;

      /* update whether or not this is a selection w/only one object */
      if (singleObjectFlag) {
        if (singleObject) {
          if (selObj != singleObject) {
            singleObjectFlag = false;
          }
        } else {
          singleObject = selObj;
        }
      }
      /* update whether or not this is a selection w/only one atom */
      if (singleAtomFlag) {
        if (singleAtom >= 0) {
          if (index != singleAtom) {
            singleAtomFlag = false;
          }
        } else {
          singleAtom = index;
        }
      }

      /* store this is the Selectors->Member table, so make sure there's room */
      c++;
      /* at runtime, selections can now have transient ordering --
         but these are not yet persistent through session saves & restores */
      SelectorManagerInsertMember(*IM, *ai, sele, tag);
    }
  }

  /* after scanning, update whether or not we touched multiple objects/atoms */
  if (c) {
    auto& info = IM->Info.back();
    if (singleObjectFlag) {
      info.theOneObject = singleObject;
      if (singleAtomFlag) {
        assert(singleAtom >= 0);
        info.theOneAtom = singleAtom;
      }
    }
  }

  if (exec_managed) {
    if (newFlag) ExecutiveManageSelection(G, name.c_str());
  }
  PRINTFD(G, FB_Selector)
  " Selector: Embedded %s, %d atoms.\n", name.c_str(), c ENDFD;
  return (c);
}
static sele_array_t SelectorApplyMultipick(PyMOLGlobals* G, Multipick* mp) {
  CSelector* I = G->Selector;
  sele_array_t result;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  sele_array_calloc(result, I->Table.size());
  for (const auto& p : mp->picked) {
    assert(p.context.object->type == cObjectMolecule);
    auto obj = static_cast<const ObjectMolecule*>(p.context.object);
    /* NOTE: SeleBase only safe with cSelectorUpdateTableAllStates!  */
    result[obj->SeleBase + p.src.index] = true;
  }
  return (result);
}
static sele_array_t SelectorSelectFromTagDict(
    PyMOLGlobals* G, const std::unordered_map<int, int>& id2tag) {
  CSelector* I = G->Selector;
  sele_array_t result{};

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates,
                      -1); /* for now, update the entire table */
  {
    sele_array_calloc(result, I->Table.size());
    if (result) {
      for (int a = cNDummyAtoms; a < I->Table.size(); a++) {
        auto& table_a = I->Table[a];
        auto ai = I->Obj[table_a.model]->AtomInfo + table_a.atom;
        if (ai->unique_id) {
          auto it = id2tag.find(ai->unique_id);
          if (it != id2tag.end()) {
            result[a] = it->second;
          }
        }
      }
    }
  }
  return (result);
}
static SelectorCreateResult_t _SelectorCreate(
    PyMOLGlobals* G, pymol::zstring_view sname, const char* sele,
    ObjectMolecule** obj, int quiet, Multipick* mp, CSeqRow* rowVLA, int nRow,
    int** obj_idx, int* n_idx, int n_obj,
    const std::unordered_map<int, int>* id2tag, int executive_manage, int state,
    SelectorID_t domain) {
  sele_array_t atom{};
  std::string name;
  int c = 0;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  ObjectMolecule* embed_obj = nullptr;

  if (!SelectorIsTmp(sname)) {
    if (sname.starts_with('%')) {
      sname.remove_prefix(1);
    }

    if (!WordMatchExact(G, cKeywordAll, sname.c_str(), ignore_case)) {
      name = UtilCleanStdStr(sname.c_str());
    }

    if (name.empty()) {
      assert(executive_manage);
      return pymol::make_error("Invalid selection name '", sname.c_str(), "'");
    }

    sname = name;
  }

  {
    if (sele) {
      auto res = SelectorSelect(G, sele, state, domain, quiet);
      p_return_if_error(res);
      atom = std::move(res.result());
    } else if (id2tag) {
      atom = SelectorSelectFromTagDict(G, *id2tag);
    } else if (obj && obj[0]) { /* optimized full-object selection */
      assert(n_obj <= 0);
      embed_obj = *obj;

      SelectorUpdateTableSingleObject(G, embed_obj,
                                      cSelectorUpdateTableAllStates);

      if (obj_idx && n_idx) {
        atom = SelectorGetSeleArrayForAtomIndices(
            G->Selector, embed_obj, *obj_idx, *n_idx, (n_obj == 0));
      }
    } else if (mp) {
      atom = SelectorApplyMultipick(G, mp);
    } else {
      return pymol::make_error(__func__, " insufficient arguments");
    }
  }

  c = SelectorEmbedSelection(G, atom.get(), sname, embed_obj, false,
                             executive_manage);
  SelectorClean(G);
  /* ignore reporting on quiet */
  if (!quiet) {
    /* ignore reporting on internal/private names */
    if (!sname.starts_with('_')) {
      PRINTFB(G, FB_Selector, FB_Actions)
      " Selector: selection \"%s\" defined with %d atoms.\n", sname.c_str(),
          c ENDFB(G);
    }
  }

  PyMOL_NeedRedisplay(G->PyMOL);
  return (c);
}

SelectorCreateResult_t SelectorCreateFromTagDict(
    PyMOLGlobals* G, const char* sname,
    const std::unordered_map<int, int>& id2tag, int exec_managed) {
  return _SelectorCreate(G, sname, nullptr, nullptr, true, nullptr, nullptr, 0,
                         nullptr, nullptr, 0, &id2tag, exec_managed, -1, -1);
}

SelectorCreateResult_t SelectorCreateEmpty(PyMOLGlobals* G, const char* name,
                                           int exec_managed) {
  return _SelectorCreate(G, name, "none", nullptr, 1, nullptr, nullptr, 0,
                         nullptr, 0, 0, nullptr, exec_managed, -1, -1);
}

SelectorCreateResult_t SelectorCreateFromObjectIndices(PyMOLGlobals* G,
                                                       const char* sname,
                                                       ObjectMolecule* obj,
                                                       int* idx, int n_idx) {
  return _SelectorCreate(G, sname, nullptr, &obj, true, nullptr, nullptr, 0,
                         &idx, &n_idx, -1, nullptr, -1, -1, -1);
  /* n_obj = -1 disables numbered tags */
}

SelectorCreateResult_t SelectorCreateOrderedFromObjectIndices(
    PyMOLGlobals* G, const char* sname, ObjectMolecule* obj, int* idx,
    int n_idx) {
  return _SelectorCreate(G, sname, nullptr, &obj, true, nullptr, nullptr, 0,
                         &idx, &n_idx, 0, nullptr, -1, -1, -1);
  /* assigned numbered tags */
}

SelectorCreateResult_t SelectorCreate(PyMOLGlobals* G, const char* sname,
                                      const char* sele, ObjectMolecule* obj,
                                      int quiet, Multipick* mp) {
  return _SelectorCreate(G, sname, sele, &obj, quiet, mp, nullptr, 0, nullptr,
                         0, 0, nullptr, -1, -1, -1);
}

SelectorCreateResult_t SelectorCreateWithStateDomain(
    PyMOLGlobals* G, const char* sname, const char* sele, ObjectMolecule* obj,
    int quiet, Multipick* mp, int state, const char* domain) {
  SelectorID_t domain_sele = cSelectionInvalid;
  ObjectNameType valid_name;

  UtilNCopy(valid_name, sname, sizeof(valid_name));
  if (SettingGetGlobal_b(G, cSetting_validate_object_names)) {
    ObjectMakeValidName(G, valid_name);
    sname = valid_name;
  }

  if (domain && domain[0]) {
    if (!WordMatchExact(G, cKeywordAll, domain, true)) { /* allow domain=all */
      domain_sele = SelectorIndexByName(G, domain);
      if (domain_sele < 0) {
        PRINTFB(G, FB_Selector, FB_Errors)
        "Selector-Error: Invalid domain selection name \"%s\".\n",
            domain ENDFB(G);
        return -1;
      }
    }
  }
  return _SelectorCreate(G, sname, sele, &obj, quiet, mp, nullptr, 0, nullptr,
                         0, 0, nullptr, -1, state, domain_sele);
}
CSelector::~CSelector() { ExecutiveInvalidateSelectionIndicatorsCGO(G); }

static void SelectorClean(PyMOLGlobals* G) {
  auto I = G->Selector;
  I->Table.clear();
  I->Obj.clear();
}
void SelectorUpdateTableSingleObject(PyMOLGlobals* G, ObjectMolecule* obj,
                                     int req_state, bool no_dummies) {
  int state = req_state;
  CSelector* I = G->Selector;

  PRINTFD(G, FB_Selector)
  "SelectorUpdateTableSingleObject-Debug: entered for %s...\n", obj->Name ENDFD;

  SelectorClean(G);

  switch (req_state) {
    case cSelectorUpdateTableAllStates:
      state = req_state;
      break;
    case cSelectorUpdateTableEffectiveStates:
      state = obj->getCurrentState();
      break;
    case cSelectorUpdateTableCurrentState:
      state = SceneGetState(G);
      break;
    default:
      if (req_state < 0) state = cSelectorUpdateTableAllStates; /* fail safe */
      break;
  }

  switch (req_state) {
    case cSelectorUpdateTableAllStates:
      I->SeleBaseOffsetsValid =
          true; /* all states -> all atoms -> offsets valid */
      break;
    default:
      I->SeleBaseOffsetsValid = false; /* not including all atoms, so atom-based
                                          offsets are invalid */
      break;
  }

  int modelCnt = cNDummyModels;
  int c = cNDummyAtoms;

  if (no_dummies) {
    modelCnt = 0;
    c = 0;
  }

  I->NCSet = obj->NCSet;
  I->Table = std::vector<TableRec>(c + obj->NAtom);
  I->Obj = std::vector<ObjectMolecule*>(modelCnt + 1, nullptr);
  I->Obj[modelCnt] = obj;

  obj->SeleBase = c;

  if (state < 0) {
    for (int atm = 0; atm < obj->NAtom; ++atm) {
      I->Table[c].model = modelCnt;
      I->Table[c].atom = atm;
      c++;
    }
  } else if (state < obj->NCSet) {
    const CoordSet* cs = obj->CSet[state];
    if (cs) {
      for (int atm = 0; atm < obj->NAtom; ++atm) {
        if (cs->atmToIdx(atm) >= 0) {
          I->Table[c].model = modelCnt;
          I->Table[c].atom = atm;
          c++;
        }
      }
    }
    I->Table.resize(c);
  }

  assert(c == I->Table.size());
}

/**
 * @param idx List of atom indices
 * @param n_idx Size of atom indices list, or -1 for -1 terminated list
 * @param numbered_tags Create ordered selection
 * @pre SelectorUpdateTableSingleObject(cSelectorUpdateTableAllStates,
 * no_dummies=false) was called
 */
static sele_array_t SelectorGetSeleArrayForAtomIndices(CSelector* I,
                                                       ObjectMolecule* obj,
                                                       const int* idx,
                                                       int n_idx,
                                                       bool numbered_tags) {
  assert(I->Obj.size() == cNDummyModels + 1);
  assert(I->Table.size() == cNDummyAtoms + obj->NAtom);

  sele_array_t result;
  sele_array_calloc(result, I->Table.size());

  if (n_idx == -1) {
    // find end of -1 terminated list, used by SeekerBuildSeleFromAtomList
    for (n_idx = 0; idx[n_idx] != -1;) {
      ++n_idx;
    }
  }

  assert(n_idx >= 0);

  int tag = numbered_tags ? SELECTOR_BASE_TAG : 1;

  for (int i = 0; i < n_idx; ++i) {
    int const atm = idx[i];

    if (atm >= 0 && atm < obj->NAtom) {
      // create an ordered selection based on the input order of the atom
      // indices
      result[obj->SeleBase + atm] = tag;
    }

    if (numbered_tags) {
      ++tag;
    }
  }

  return result;
}
int SelectorUpdateTable(PyMOLGlobals* G, int req_state, SelectorID_t domain) {
  return (SelectorUpdateTableImpl(G, G->Selector, req_state, domain));
}

int SelectorUpdateTableImpl(PyMOLGlobals* G, CSelector* I, int req_state,
                            SelectorID_t domain) {
  int a = 0;
  ov_size c = 0;
  int modelCnt;
  int state = req_state;
  void* iterator = nullptr;
  ObjectMolecule* obj = nullptr;

  /* Origin and Center are dummy objects */
  if (!I->Origin)
    I->Origin.reset(ObjectMoleculeDummyNew(G, cObjectMoleculeDummyOrigin));

  if (!I->Center)
    I->Center.reset(ObjectMoleculeDummyNew(G, cObjectMoleculeDummyCenter));

  SelectorClean(G);
  I->NCSet = 0;

  /* take a summary of PyMOL's current state; foreach molecular object
   * sum up the number of atoms, count how many models, states, etc... */
  modelCnt = cNDummyModels;
  c = cNDummyAtoms;
  while (ExecutiveIterateObjectMolecule(G, &obj, &iterator)) {
    c += obj->NAtom;
    if (I->NCSet < obj->NCSet) I->NCSet = obj->NCSet;
    modelCnt++;
  }
  /* allocate space for each atom, in the record table */
  I->Table = std::vector<TableRec>(c);
  I->Obj = std::vector<ObjectMolecule*>(modelCnt, nullptr);

  switch (req_state) {
    case cSelectorUpdateTableAllStates:
      I->SeleBaseOffsetsValid =
          true; /* all states -> all atoms -> offsets valid */
      break;
    default:
      I->SeleBaseOffsetsValid = false; /* not including all atoms, so atom-based
                                          offsets are invalid */
      break;
  }

  c = 0;
  modelCnt = 0;

  /* update the origin and center dummies */
  obj = I->Origin.get();
  if (obj) {
    I->Obj[modelCnt] = I->Origin.get();
    obj->SeleBase = c; /* make note of where this object starts */
    for (a = 0; a < obj->NAtom; a++) {
      I->Table[c].model = modelCnt;
      I->Table[c].atom = a;
      c++;
    }
    modelCnt++;
  }

  obj = I->Center.get();
  if (obj) {
    I->Obj[modelCnt] = I->Center.get();
    obj->SeleBase = c; /* make note of where this object starts */
    for (a = 0; a < obj->NAtom; a++) {
      I->Table[c].model = modelCnt;
      I->Table[c].atom = a;
      c++;
    }
    modelCnt++;
  }

  while (ExecutiveIterateObjectMolecule(G, &obj, &iterator)) {
    int skip_flag = false;
    if (req_state < 0) {
      switch (req_state) {
        case cSelectorUpdateTableAllStates:
          state = -1; /* all states */
          /* proceed... */
          break;
        case cSelectorUpdateTableCurrentState:
          state = SettingGetGlobal_i(G, cSetting_state) - 1;
          break;
        case cSelectorUpdateTableEffectiveStates:
          state = obj->getCurrentState();
          break;
        default: /* unknown input -- fail safe (all states) */
          state = -1;
          break;
      }
    } else {
      if (state >= obj->NCSet)
        skip_flag = true;
      else if (!obj->CSet[state])
        skip_flag = true;
    }

    if (!skip_flag) {
      /* fill in the table */
      I->Obj[modelCnt] = obj;
      {
        int n_atom = obj->NAtom;
        auto rec = I->Table.data() + c;
        TableRec* start_rec = rec;
        if (state < 0) {    /* all states */
          if (domain < 0) { /* domain=all */
            for (a = 0; a < n_atom; a++) {
              rec->model = modelCnt;
              rec->atom = a;
              rec++;
            }
          } else {
            const AtomInfoType* ai = obj->AtomInfo.data();
            int included_one = false;
            int excluded_one = false;
            for (a = 0; a < n_atom; a++) {
              if (SelectorIsMember(G, ai->selEntry, domain)) {
                rec->model = modelCnt;
                rec->atom = a;
                rec++;
                included_one = true;
              } else {
                excluded_one = true;
              }
              ai++;
            }
            if (included_one && excluded_one)
              I->SeleBaseOffsetsValid = false; /* partial objects in domain, so
                                                  base offsets are invalid */
          }
        } else { /* specific states */
          CoordSet* cs;
          int idx;
          if (domain < 0) {
            for (a = 0; a < n_atom; a++) {
              /* does coordinate exist for this atom in the requested state? */
              if (state < obj->NCSet)
                cs = obj->CSet[state];
              else
                cs = nullptr;
              if (cs) {
                idx = cs->atmToIdx(a);
                if (idx >= 0) {
                  rec->model = modelCnt;
                  rec->atom = a;
                  rec++;
                }
              }
            }
          } else {
            const AtomInfoType* ai = obj->AtomInfo.data();
            for (a = 0; a < n_atom; a++) {
              /* does coordinate exist for this atom in the requested state? */
              if (state < obj->NCSet)
                cs = obj->CSet[state];
              else
                cs = nullptr;
              if (cs) {
                idx = cs->atmToIdx(a);
                if (idx >= 0) {
                  if (SelectorIsMember(G, ai->selEntry, domain)) {
                    rec->model = modelCnt;
                    rec->atom = a;
                    rec++;
                  }
                }
              }
              ai++;
            }
          }
        }
        if (rec != start_rec) { /* skip excluded models */
          modelCnt++;
          obj->SeleBase = c; /* make note of where this object starts */
          c += (rec - start_rec);
        } else {
          obj->SeleBase = 0;
        }
      }
    }
  }
  I->Obj.resize(modelCnt);
  I->Table.resize(c);
  /* printf("selector update table state=%d, natom=%d\n",req_state,c); */
  return (true);
}
void SelectorMemoryDump(PyMOLGlobals* G) {
  auto I = G->SelectorMgr;
  printf(" SelectorMemory: NSelection %d\n", I->NSelection);
  printf(" SelectorMemory: NActive %zu\n", I->Info.size());
  printf(" SelectorMemory: NMember %d\n", int(I->Member.size()) - 1);
}

CSelectorManager::CSelectorManager() {
  auto I = this;

  // indices are >0 by convention
  Member.resize(1);

  /* create placeholder "all" selection, which is selection 0
     and "none" selection, which is selection 1 */
  I->Info.emplace_back(I->NSelection++, cKeywordAll);
  I->Info.emplace_back(I->NSelection++, cKeywordNone);

  assert(I->Info[0].ID == cSelectionAll);
  assert(I->Info[1].ID == cSelectionNone);

  for (auto kw : Keyword) {
    if (!kw.word[0]) {
      break;
    }
    I->Key[kw.word] = kw.value;
  }
}

void SelectorReinit(PyMOLGlobals* G) {
  SelectorClean(G);
  *G->SelectorMgr = CSelectorManager();
}
CSelector::CSelector(PyMOLGlobals* G, CSelectorManager* mgr) : G(G), mgr(mgr) {}
/*

example selections
cas
backbone
(model 1 and backbone)
(name ca)
(resi 123)
(resi 200:400 and chain A)
(model a)

 */

/* In order for selections to be robust during atom insertions
   deletions, they are stored not as lists of selected atoms, but
   rather in the inverse - as atoms with selection membership
   information */

/* Each atom points into the selection heap, a series of linked
   entries where each member is selection index */

/* Definition of the selection language:

        <sele> = [(] [<not>] <prop> <val-range> [<qual> [<qual-range>] ] [)] {
   <SEL1> <sele> }

        Example selections:

   name ca
        ( name ca )
        name ca around 5 {all atoms within 5 angstroms of any ca atom) }
        ( resi 10 )
        resi 10:50
        resi 10A:50A
        resi 10,50
        chain A
        segi A
        model 1 and segi B
        model a and segi B
        not name ca

*/

/* Selection processing, left to right by default, but with parenthesis for
   grouping stack based functional language processing each stack level has a
   full selection matrix logical binary operators (or,and) and negation
   SEL1ation solely on these matrices.
*/

/*

(not (name ca or name c) and (name s around 5 and (name c) )) around 6

0:
1: not
2: name
2: ca

0:

1: not
2:<name ca>
not name ca around 5

force compute

0:
1: not

attrib b < 0

*/

int SelectorSetAtomPropertyForSelection(PyMOLGlobals* G, int sele,
                                        const char* propname, CPythonVal* value,
                                        const PropertyType& proptype, int state,
                                        int quiet) {
  CSelector* I = G->Selector;
  int nAtom = 0;
  int ok = true;

  SelectorUpdateTable(G, state, -1);

  if (ok) {
    int a;
    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      int at = I->Table[a].atom;
      ObjectMolecule* obj = I->Obj[I->Table[a].model];
      int s = obj->AtomInfo[at].selEntry;
      I->Table[a].index = 0;
      if (SelectorIsMember(G, s, sele)) {
        nAtom++;
        PropertySet(G, &obj->AtomInfo[at], propname, value, proptype);
      }
    }
  }
  return nAtom;
}
