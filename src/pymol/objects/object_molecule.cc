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

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/geometry/matrix.h"
#include "pymol/python/p.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"

#include "pymol/objects/property.h"

#define cMaxNegResi 100

#include "pymol/objects/object_molecule_internal.h"

/**
 * Order function for sorting bonds by atom indices (ignores other properties
 * like bond order).
 *
 * Pre-condition: index pairs are in order (index[0] < index[1])
 */
static int BondTypeInOrder(PyMOLGlobals* G, const BondType* bonds, int i1,
                           int i2) {
  auto lhs = bonds[i1].index;
  auto rhs = bonds[i2].index;
  return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

/**
 * Remove duplicated bonds. Disregards if two atoms would be bonded by two
 * bonds of different bond order (only one will be kept).
 */
static void ObjectMoleculeRemoveDuplicateBonds(PyMOLGlobals* G,
                                               ObjectMolecule* I) {
  // make sure index pairs are in order
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    if (bond.index[0] > bond.index[1] && !bond.symop_2) {
      std::swap(bond.index[0], bond.index[1]);
    }
  }

  // get sorted indices
  int* sorted = pymol::malloc<int>(I->NBond);
  UtilSortIndexGlobals(G, I->NBond, I->Bond.data(), sorted,
                       (UtilOrderFnGlobals*)BondTypeInOrder);

  // purge duplicated bonds and mark for removal
  for (int i = 0, j = -1; i < I->NBond; ++i) {
    int b = sorted[i];
    auto ptr = I->Bond + sorted[i];

    bool equal = false;
    if (j != -1) {
      const auto& other = I->Bond[j];
      if (ptr->index[0] == other.index[0] && ptr->index[1] == other.index[1])
        equal = true;
    }

    if (!equal) {
      j = b;
    } else {
      AtomInfoPurgeBond(G, ptr);

      // mark for removal
      ptr->index[0] = ptr->index[1] = 0;
    }
  }

  FreeP(sorted);

  // remove purged bonds from array
  int j = 0;
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    if (bond.index[0] != bond.index[1]) {
      if (j != i) std::swap(I->Bond[j], bond);
      ++j;
    }
  }
  I->NBond = j;
  VLASize(I->Bond, BondType, I->NBond);
}

/**
 * Convert a discrete object to non-discrete. Will merge atoms from states
 * by identifiers. Atom level properties like color, b, q, etc. are lost
 * on the merged atoms.
 */
static bool ObjectMoleculeSetNotDiscrete(PyMOLGlobals* G, ObjectMolecule* I) {
  if (!I->DiscreteFlag) return true;

  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  I->DiscreteFlag = false;

  for (int atm = 0; atm < I->NAtom; ++atm) {
    I->AtomInfo[atm].discrete_state = 0;
  }

  int* outdex;
  int* aindex = AtomInfoGetSortedIndex(G, I, I->AtomInfo, I->NAtom, &outdex);

  // merge duplicate atoms
  for (int i = 0, j = -1; i < I->NAtom; ++i) {
    int atm = aindex[i];
    auto ai = I->AtomInfo + atm;

    if (j == -1 || !AtomInfoMatch(G, ai, I->AtomInfo + j, false, false)) {
      j = atm;
    } else {
      ai->deleteFlag = true;
    }

    // remapping of merged indices
    outdex[atm] = j;
  }

  // point coordsets to merged atoms
  for (int state = 0; state < I->NCSet; ++state) {
    auto cs = I->CSet[state];
    if (!cs) continue;

    for (int idx = 0; idx < cs->NIndex; ++idx) {
      int atm = cs->IdxToAtm[idx];
      cs->IdxToAtm[idx] = outdex[atm];
    }
  }

  // point bonds to merged atoms
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    bond.index[0] = outdex[bond.index[0]];
    bond.index[1] = outdex[bond.index[1]];
  }

  AtomInfoFreeSortedIndexes(G, &aindex, &outdex);

  ObjectMoleculeRemoveDuplicateBonds(G, I);

  I->updateAtmToIdx();
  ObjectMoleculePurge(I);

  return true;
}

/**
 * Make a non-discrete object discrete, or vice versa.
 */
int ObjectMoleculeSetDiscrete(PyMOLGlobals* G, ObjectMolecule* I,
                              int discrete) {
  int state, idx, ao, an, ao1, ao2, an1, an2;
  int maxnatom, natom = I->NAtom, nbond = I->NBond;
  int* aostate2an = nullptr;
  char* bondseen = nullptr;
  BondType* bond;
  CoordSet* cs;

  if (!discrete) {
    if (!I->DiscreteFlag) return true;

    return ObjectMoleculeSetNotDiscrete(G, I);
  }

  if (I->DiscreteFlag) return true;

  // upper bound for number of discrete atoms
  maxnatom = I->NAtom * I->NCSet;

  // mapping (for bonds): atom_old -> atom_new
  ok_assert(1, aostate2an = pymol::malloc<int>(I->NAtom));
  ok_assert(1, bondseen = pymol::calloc<char>(I->NBond));

  // discrete setup
  I->DiscreteFlag = discrete;
  ok_assert(1, I->DiscreteAtmToIdx = pymol::vla<int>(maxnatom));
  ok_assert(1, I->DiscreteCSet = pymol::vla<CoordSet*>(maxnatom));

  // for all coordinate sets
  for (state = 0; state < I->NCSet; state++) {
    cs = I->CSet[state];

    if (!cs) continue;

    // init the atom_old -> atom_new array
    for (ao = 0; ao < I->NAtom; ao++) aostate2an[ao] = -1;

    // for all atoms in coordinate set
    for (idx = 0; idx < cs->NIndex; idx++) {
      ao = an = cs->IdxToAtm[idx];

      if (I->DiscreteCSet[ao]) {
        // seen before, have to copy
        an = natom++;

        VLACheck(I->AtomInfo, AtomInfoType, an);
        ok_assert(1, I->AtomInfo);

        AtomInfoCopy(G, I->AtomInfo + ao, I->AtomInfo + an);
        cs->IdxToAtm[idx] = an;
      }

      I->AtomInfo[an].discrete_state = state + 1;  // 1-based :-(
      I->DiscreteCSet[an] = cs;
      I->DiscreteAtmToIdx[an] = idx;

      // mapping (for bonds): (atom_old, state) -> atom_new
      aostate2an[ao] = an;
    }

    // unused in discrete coordsets
    cs->AtmToIdx.clear();

    // for all old bonds
    for (idx = 0; idx < I->NBond; idx++) {
      bond = I->Bond + idx;

      // old atom indices
      ao1 = bond->index[0];
      ao2 = bond->index[1];

      // new atom indices
      an1 = aostate2an[ao1];
      an2 = aostate2an[ao2];

      // do both atoms exist in this coord set?
      if (an1 == -1 || an2 == -1) continue;

      if (bondseen[idx]) {
        // seen before, have to copy
        VLACheck(I->Bond, BondType, nbond);
        ok_assert(1, I->Bond);

        bond = I->Bond + (nbond++);
        AtomInfoBondCopy(G, I->Bond + idx, bond);
      } else {
        bondseen[idx] = true;
      }

      bond->index[0] = an1;
      bond->index[1] = an2;
    }
  }

  // not needed anymore
  mfree(aostate2an);
  mfree(bondseen);

  // update N* count fields
  I->NAtom = natom;
  I->NBond = nbond;

  // trim VLAs memory to actual size
  if (I->NBond) VLASize(I->Bond, BondType, I->NBond);
  if (I->NAtom) VLASize(I->AtomInfo, AtomInfoType, I->NAtom);

  I->setNDiscrete(I->NAtom);

  I->invalidate(cRepAll, cRepInvAll, -1);

  return true;

ok_except1:
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
  " %s: memory allocation failed\n", __func__ ENDFB(G);
  return false;
}

int ObjectMoleculeCheckFullStateSelection(ObjectMolecule* I, int sele,
                                          int state) {
  PyMOLGlobals* G = I->G;
  int result = false;
  if ((state >= 0) && (state < I->NCSet)) {
    const AtomInfoType* ai = I->AtomInfo.data();
    CoordSet* cs = I->CSet[state];
    if (cs) {
      int a;
      int at;
      result = true;
      for (a = 0; a < cs->NIndex; a++) {
        at = cs->IdxToAtm[a];
        if (!SelectorIsMember(G, ai[at].selEntry, sele)) {
          result = false;
          break;
        }
      }
    }
  }
  return result;
}

/* PARAMS
 *   ch -- ptr to empty str of length 'len'
 *  len -- str len
 * RETURNS
 *   ch, with the caption in place
 * NOTES
 *   User owns the buffer so must clean up after it
 */
char* ObjectMolecule::getCaption(char* ch, int len) const {
  auto I = this;
  int objState;
  int n = 0;
  int show_state = 0;
  int show_as_fraction = 0;
  const char* frozen_str = "";

  int state = ObjectGetCurrentState(I, false);
  int counter_mode = SettingGet_i(I->G, I->Setting.get(), nullptr,
                                  cSetting_state_counter_mode);
  int frozen =
      SettingGetIfDefined_i(I->G, I->Setting.get(), cSetting_state, &objState);

  /* if frozen print (blue) STATE / NSTATES
   * if not frozen, print STATE/NSTATES
   * if beyond NSTATES, print * /NSTATES.
   */

  if (frozen) { /* frozen color */
    frozen_str = "\\789";
  } else {
    if (I->DiscreteFlag) { /* discrete states */
      frozen_str = "\\993";
    } else { /* normal case */
      frozen_str = "";
    }
  }

  switch (counter_mode) {
    case 0: /* off */
      show_state = show_as_fraction = 0;
      break;
    case 2: /* just state */
      show_state = 1;
      show_as_fraction = 0;
      break;
    case -1: /* fraction, full-on */
    case 1:
    default:
      show_state = show_as_fraction = 1;
      break;
  }

  /* bail on null string or no room */
  if (!ch || len == 0) return nullptr;

  ch[0] = 0;

  /* if the state is valid, setup the label */
  if (state >= 0) {
    if (state < I->NCSet) {
      const CoordSet* cs = I->CSet[state];
      if (cs) {
        if (show_state) {
          if (show_as_fraction) {
            if (strlen(cs->Name)) { /* NAME */
              n = snprintf(ch, len, "%s %s%d/%d", cs->Name, frozen_str,
                           state + 1, I->NCSet);
            } else { /* no name */
              n = snprintf(ch, len, "%s%d/%d", frozen_str, state + 1, I->NCSet);
            }
          } else { /* not fraction */
            if (strlen(cs->Name)) {
              n = snprintf(ch, len, "%s %s%d", cs->Name, frozen_str, state + 1);
            } else { /* no name */
              n = snprintf(ch, len, "%s%d", frozen_str, state + 1);
            }
          }
        } else { /* no state */
          n = snprintf(ch, len, "%s", cs->Name);
        }
      } else { /* no coord set, can be an N-state object missing a CSet */
      }
    } else { /* state > NCSet, out of range due to other object or global
                setting */
      if (show_state) {
        if (show_as_fraction) {
          n = snprintf(ch, len, "%s--/%d", frozen_str, I->NCSet);
        } else { /* no fraction */
          n = snprintf(ch, len, "%s--", frozen_str);
        }
      }
    }
  } else if (state == -1) {
    // all states
    n = snprintf(ch, len, "%s*/%d", frozen_str, I->NCSet);
  } else {
    /* blank out the title if outside the valid # of states */
  }

  if (n > len) return nullptr;

  return ch;
}

#define MAX_BOND_DIST 50

void ObjectMoleculeTransformState44f(ObjectMolecule* I, int state,
                                     const float* matrix, int log_trans,
                                     int homogenous, int transformed) {
  int a;
  float tmp_matrix[16];
  CoordSet* cs;
  int use_matrices =
      SettingGet_i(I->G, I->Setting.get(), nullptr, cSetting_matrix_mode);
  if (use_matrices < 0) use_matrices = 0;
  if (!use_matrices) {
    ObjectMoleculeTransformSelection(I, state, -1, matrix, log_trans, I->Name,
                                     homogenous, true);
  } else {
    double dbl_matrix[16];
    if (state == -2) state = ObjectGetCurrentState(I, false);
    /* ensure homogenous matrix to preserve programmer sanity */
    if (!homogenous) {
      convertTTTfR44d(matrix, dbl_matrix);
      copy44d44f(dbl_matrix, tmp_matrix);
      matrix = tmp_matrix;
    } else {
      copy44f44d(matrix, dbl_matrix);
    }

    if (state < 0) { /* all states */
      for (a = 0; a < I->NCSet; a++) {
        cs = I->CSet[a];
        if (cs) ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
      }
    } else if (state < I->NCSet) { /* single state */
      cs = I->CSet[state];
      if (cs) ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
    } else if (I->NCSet == 1) { /* static singleton state */
      cs = I->CSet[0];
      if (cs && SettingGet_b(I->G, I->Setting.get(), nullptr,
                             cSetting_static_singletons)) {
        ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
      }
    }
  }
}

void ObjectMoleculeOpRecInit(ObjectMoleculeOpRec* op) {
  UtilZeroMem((char*)op, sizeof(ObjectMoleculeOpRec));
}

bool ObjectMoleculeGetNeighborVector(ObjectMolecule* I, int atom, int state,
                                     float* v_result) {
  float v_atom[3] = {0.0, 0.0, 0.0};

  auto const* cs = I->getCoordSet(state);
  if (cs == nullptr) {
    return false;
  }

  if (!CoordSetGetAtomVertex(cs, atom, v_atom)) {
    return false;
  }

  for (auto const& neighbor : AtomNeighbors(I, atom)) {
    auto const a2 = neighbor.atm;

    // ignore hydrogens
    if (I->AtomInfo[a2].protons != cAN_H &&
        CoordSetGetAtomVertex(cs, a2, v_result)) {
      return true;
    }
  }

  return false;
}

/**
 * Returns the most proton-rich element with the lowest priority value (OG1
 * before OG2, CG, HB1)
 */
int ObjectMoleculeGetTopNeighbor(PyMOLGlobals* G, ObjectMolecule* I, int start,
                                 int excluded) {
  int highest_at = -1, highest_prot = 0, lowest_pri = 9999;

  for (auto const& neighbor : AtomNeighbors(I, start)) {
    auto const at = neighbor.atm;
    auto const* ai = I->AtomInfo.data() + at;
    if ((highest_at < 0) && (at != excluded)) {
      highest_prot = ai->protons;
      lowest_pri = ai->priority;
      highest_at = at;
    } else if (((ai->protons > highest_prot) ||
                ((ai->protons == highest_prot) &&
                 (ai->priority < lowest_pri))) &&
               (at != excluded)) {
      highest_prot = ai->protons;
      highest_at = at;
      lowest_pri = ai->priority;
    }
  }
  return highest_at;
}

void ObjectMoleculeSculptClear(ObjectMolecule* I) {
  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: entered.\n", __func__ ENDFD;
  if (I->Sculpt) DeleteP(I->Sculpt);
}

void ObjectMoleculeSculptImprint(ObjectMolecule* I, int state, int match_state,
                                 int match_by_segment) {
  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: entered.\n", __func__ ENDFD;

  if (!I->Sculpt) I->Sculpt = new CSculpt(I->G);
  SculptMeasureObject(I->Sculpt, I, state, match_state, match_by_segment);
}

float ObjectMoleculeSculptIterate(ObjectMolecule* I, int state, int n_cycle,
                                  float* center) {
  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: entered.\n", __func__ ENDFD;
  if (I->Sculpt) {
    return SculptIterateObject(I->Sculpt, I, state, n_cycle, center);
  } else
    return 0.0F;
}

/**
 * - Assigns new id to all atoms with AtomInfoType.id == -1
 * - Assigns ObjectMolecule.AtomCounter if -1
 *
 * Cost: O(NAtom + NBond)
 */
void ObjectMoleculeUpdateIDNumbers(ObjectMolecule* I) {
  int a;
  int max;
  AtomInfoType* ai;

  if (I->AtomCounter < 0) {
    max = -1;
    ai = I->AtomInfo.data();
    for (a = 0; a < I->NAtom; a++) {
      if (ai->id > max) max = ai->id;
      ai++;
    }
    I->AtomCounter = max + 1;
  }
  ai = I->AtomInfo.data();
  for (a = 0; a < I->NAtom; a++) {
    if (ai->id < 0) ai->id = I->AtomCounter++;
    ai++;
  }
}
int ObjectMoleculeGetPhiPsi(ObjectMolecule* I, int ca, float* phi, float* psi,
                            int state) {
  int np = -1;
  int cm = -1;
  int c = -1;
  int n = -1;
  int result = false;
  const AtomInfoType* ai;
  float v_ca[3];
  float v_n[3];
  float v_c[3];
  float v_cm[3];
  float v_np[3];
  auto G = I->G;

  ai = I->AtomInfo;

  if (ai[ca].name == G->lex_const.CA) {
    /* find C */
    for (auto const& neighbor : AtomNeighbors(I, ca)) {
      if (ai[neighbor.atm].name == G->lex_const.C) {
        c = neighbor.atm;
        break;
      }
    }

    /* find N */
    for (auto const& neighbor : AtomNeighbors(I, ca)) {
      if (ai[neighbor.atm].name == G->lex_const.N) {
        n = neighbor.atm;
        break;
      }
    }

    /* find NP */
    if (c >= 0) {
      for (auto const& neighbor : AtomNeighbors(I, c)) {
        if (ai[neighbor.atm].name == G->lex_const.N) {
          np = neighbor.atm;
          break;
        }
      }
    }

    /* find CM */
    if (n >= 0) {
      for (auto const& neighbor : AtomNeighbors(I, n)) {
        if (ai[neighbor.atm].name == G->lex_const.C) {
          cm = neighbor.atm;
          break;
        }
      }
    }
    if ((ca >= 0) && (np >= 0) && (c >= 0) && (n >= 0) && (cm >= 0)) {
      if (ObjectMoleculeGetAtomVertex(I, state, ca, v_ca) &&
          ObjectMoleculeGetAtomVertex(I, state, n, v_n) &&
          ObjectMoleculeGetAtomVertex(I, state, c, v_c) &&
          ObjectMoleculeGetAtomVertex(I, state, cm, v_cm) &&
          ObjectMoleculeGetAtomVertex(I, state, np, v_np)) {
        (*phi) = rad_to_deg(get_dihedral3f(v_c, v_ca, v_n, v_cm));
        (*psi) = rad_to_deg(get_dihedral3f(v_np, v_c, v_ca, v_n));
        result = true;
      }
    }
  }
  return (result);
}
int ObjectMoleculeCheckBondSep(ObjectMolecule* I, int a0, int a1, int dist) {
  int result = false;
  int n0;
  int stack[MAX_BOND_DIST + 1];
  int history[MAX_BOND_DIST + 1];
  int depth = 0;
  int distinct;
  int a;
  if (dist > MAX_BOND_DIST) return false;

  auto* const Neighbor = I->getNeighborArray();

  depth = 1;
  history[depth] = a0;
  stack[depth] = Neighbor[a0] + 1; /* go to first neighbor */
  while (depth) {                  /* keep going until we've traversed tree */
    while (Neighbor[stack[depth]] >=
           0) {                    /* end of branches? go back up one bond */
      n0 = Neighbor[stack[depth]]; /* get current neighbor index */
      stack[depth] += 2;           /* set up next neighbor */
      distinct = true; /* check to see if current candidate is distinct from
                          ancestors */
      for (a = 1; a < depth; a++) {
        if (history[a] == n0) distinct = false;
      }
      if (distinct) {
        if (depth < dist) { /* are not yet at the proper distance? */
          if (distinct) {
            depth++;
            stack[depth] = Neighbor[n0] + 1; /* then keep moving outward */
            history[depth] = n0;
          }
        } else if (n0 == a1) /* otherwise, see if we have a match */
          result = true;
      }
    }
    depth--;
  }
  return result;
}
void ObjectGotoState(pymol::CObject* I, int state) {
  auto nstates = I->getNFrame();
  if (nstates > 1 || !SettingGet<bool>(I->G, cSetting_static_singletons)) {
    if (state > nstates) state = nstates - 1;
    if (state < 0) state = nstates - 1;
    SceneSetFrame(I->G, 0, state);
  }
}
CObjectState* ObjectMolecule::_getObjectState(int state) { return CSet[state]; }
pymol::copyable_ptr<CSetting>* ObjectMolecule::getSettingHandle(int state) {
  auto I = this;

  if (state < -1) {
    state = I->getCurrentState();
  }

  if (state < 0) {
    return (&I->Setting);
  } else if (state < I->NCSet) {
    if (I->CSet[state]) {
      return (&I->CSet[state]->Setting);
    } else {
      return (nullptr);
    }
  } else {
    return (nullptr);
  }
}
CSymmetry const* ObjectMolecule::getSymmetry(int state) const {
  auto cs = getCoordSet(state);
  if (cs && cs->Symmetry) {
    return cs->Symmetry.get();
  }

  return Symmetry.get();
}

/**
 * @return False if state does not exist
 */
bool ObjectMolecule::setSymmetry(CSymmetry const& symmetry, int state) {
  bool const all_states = (state == cSelectorUpdateTableAllStates);
  bool success = false;

  if (all_states) {
    Symmetry.reset(new CSymmetry(symmetry));
    success = true;
  }

  for (StateIterator iter(G, Setting.get(), state, NCSet); iter.next();) {
    auto cs = CSet[iter.state];
    if (cs) {
      cs->Symmetry.reset(all_states ? nullptr : new CSymmetry(symmetry));
      cs->invalidateRep(cRepCell, cRepInvRep);
      success = true;
    }
  }

  return success;
}
pymol::Result<> ObjectMoleculeSetStateTitle(ObjectMolecule* I, int state,
                                            const char* text) {
  auto cs = I->getCoordSet(state);
  if (!cs) {
    return pymol::make_error("Invalid state ", state + 1);
  }
  cs->setTitle(text);
  return {};
}
const char* ObjectMoleculeGetStateTitle(const ObjectMolecule* I, int state) {
  auto cs = I->getCoordSet(state);
  if (!cs) {
    PRINTFB(I->G, FB_ObjectMolecule, FB_Errors)
    "Error: invalid state %d\n", state + 1 ENDFB(I->G);
    return nullptr;
  }
  return cs->Name;
}
int ObjectMoleculeGetAtomIndex(const ObjectMolecule* I, SelectorID_t sele) {
  if (sele < 0) return (-1);

  for (int a = 0; a < I->NAtom; a++) {
    auto const s = I->AtomInfo[a].selEntry;
    if (SelectorIsMember(I->G, s, sele)) return (a);
  }
  return (-1);
}
/**
 * Update all AtomInfoType.bonded
 *
 * Cost: O(NAtom + NBond)
 */
void ObjectMoleculeUpdateNonbonded(ObjectMolecule* I) {
  int a;
  const BondType* b;
  AtomInfoType* ai;
  int nAtom = I->NAtom;
  int nBond = I->NBond;

  ai = I->AtomInfo.data();

  for (a = 0; a < nAtom; a++) (ai++)->bonded = false;

  b = I->Bond;
  ai = I->AtomInfo.data();
  for (a = 0; a < nBond; a++) {
    ai[b->index[0]].bonded = true;
    ai[b->index[1]].bonded = true;
    b++;
  }
}

/**
 * Get the raw ObjectMolecule::Neighbor storage structure, which is generated
 * from the ObjectMolecule::Bond array. This structure is cached, to clear the
 * cache, call ObjectMolecule::invalidate(level=cRepInvBonds).
 *
 * Direct usage of this array should be avoided, use the AtomNeighbors() class
 * instead which is a higher-level and more readable interface to the same
 * data.
 *
 * Changed in PyMOL 2.1.1: Ignore zero-order bonds (PYMOL-3025)
 *
 * @return Pointer to cached array, or nullptr if memory allocation failed
 *
 * @verbatim
     0       list offset for atom 0 = n
     1       list offset for atom 1 = n + m + 1
     ...
     n-1     list offset for atom n-1

     n       count for atom 0
     n+1     neighbor of atom 0
     n+2     bond index
     n+3     neighbor of atom 0
     n+4     bond index
     ...
     n+m     -1 terminator for atom 0

     n+m+1   count for atom 1
     n+m+2   neighbor of atom 1
     n+m+3   bond index
     n+m+4   neighbor of atom 1
     n+m+5   bond index
     etc.

     NOTE: all atoms have an offset and a terminator whether or not they have
 any bonds:

     Here's an example of an ALA
     [ offet for atom1, offset for atom2, ..., offset for atom9, \
     {10, 16, 26, 32, 36, 46, 50, 54, 58, 62, \
     (atom1)  nbonds=2, neighbor index1=5, bond index=1, neighbor index2=1, bond
 index=0, null terminator=-1 \ 2, 5, 1, 1, 0, -1,
     (atom2) nbonds=4, nbr index1=6, bond id=4; nbr idx2=4, bond idx=3; nbr
 index3=2, bond idx=2; nbr index=0, bond index=0, null=-1 \ 4, 6, 4, 4, 3, 2, 2,
 0, 0, -1, (atom3) ... 2, 3, 5, 1, 2, -1, 1, 2, 5, -1, 4, 9, 8, 8, 7, 7, 6, 1,
 3, -1, 1, 0, 1, -1, 1, 1, 4, -1, 1, 4, 6, -1, 1, 4, 7, -1, 1, 4, 8, -1 }
   @endverbatim

   */
int const* ObjectMolecule::getNeighborArray() const {
  auto const I = this;

  /* If no neighbors have been calculated, fill in the table */
  if (!I->Neighbor) {
    /* Create/check the VLA */
    auto const size = (I->NAtom * 3) + (I->NBond * 4);

    const_cast<ObjectMolecule*>(this)->Neighbor.reset(new int[size]);
    auto* const neighbor = this->Neighbor.get();

    p_return_val_if_fail(neighbor, nullptr);

    /* initialize; zero out neighbors */
    std::fill_n(neighbor, I->NAtom, 0);

    /* count neighbors for each atom */
    const auto* bnd = I->Bond.data();
    for (int b = 0; b < I->NBond; b++) {
      if (bnd->order && !bnd->hasSymOp()) {
        neighbor[bnd->index[0]]++;
        neighbor[bnd->index[1]]++;
      }
      bnd++;
    }

    /* set up offsets and list terminators */
    auto c = I->NAtom;
    for (int a = 0; a < I->NAtom; a++) {
      auto const d = neighbor[a]; /* get number of neighbors */
      neighbor[c] = d;            /* store neighbor count */
      neighbor[a] =
          c + d + d +
          1; /* set initial position to end of list, we'll fill backwards */
      neighbor[neighbor[a]] = -1; /* store terminator */
      c += d + d + 2;
    }

    /* now load neighbors in a sequential list for each atom (reverse order) */
    bnd = I->Bond.data();
    for (int b = 0; b < I->NBond; b++, ++bnd) {
      auto const l0 = bnd->index[0];
      auto const l1 = bnd->index[1];

      if (!bnd->order || bnd->hasSymOp()) continue;

      neighbor[l0]--;
      neighbor[neighbor[l0]] = b; /* store bond indices (for I->Bond) */
      neighbor[l0]--;
      neighbor[neighbor[l0]] =
          l1; /* store neighbor references (I->AtomInfo, etc.) */

      neighbor[l1]--;
      neighbor[neighbor[l1]] = b; /* store bond indices (for I->Bond) */
      neighbor[l1]--;
      neighbor[neighbor[l1]] =
          l0; /* store neighbor references (I->AtomInfo, etc.) */
    }

    // adjust down to point to the count, not the first entry
    for (int a = 0; a < I->NAtom; a++) {
      if (neighbor[a] >= 0) --neighbor[a];
    }
  }
  return this->Neighbor.get();
}
typedef int CompareFn(PyMOLGlobals*, const AtomInfoType*, const AtomInfoType*);
int ObjectMoleculeMerge(ObjectMolecule* I, pymol::vla<AtomInfoType>&& ai,
                        CoordSet* cs, int bondSearchFlag, int aic_mask,
                        int invalidate) {
  PyMOLGlobals* G = I->G;
  int *index, *outdex;
  int a, b, lb = 0, ac;
  int c, nb, a1, a2;
  int found;
  int nAt, nBd, nBond;
  int expansionFlag = false;
  int ok = true;

  const auto n_index = cs->getNIndex();
  const auto oldNAtom = I->NAtom;
  const auto oldNBond = I->NBond;

  /* first, sort the coodinate set */

  index = AtomInfoGetSortedIndex(G, I, ai, n_index, &outdex);
  CHECKOK(ok, index);
  if (!ok) return false;

  auto ai2 = pymol::vla<AtomInfoType>(n_index);
  if (ok) {
    for (a = 0; a < n_index; a++)
      ai2[a] = std::move(
          ai[index[a]]); /* creates a sorted list of atom info records */
  }
  ai = std::move(ai2);

  /* now, match it up with the current object's atomic information */

  if (ok) {
    for (a = 0; a < n_index; a++) {
      index[a] = -1;
    }
  }

  if (ok) {
    const auto n_atom = I->NAtom;
    const AtomInfoType* const atInfo = I->AtomInfo.data();
    CompareFn* fCompare;

    if (SettingGetGlobal_b(G, cSetting_pdb_hetatm_sort)) {
      fCompare = AtomInfoCompareIgnoreRank;
    } else {
      fCompare = AtomInfoCompareIgnoreRankHet;
    }

    c = 0;
    b = 0;
    lb = 0;

    for (a = 0; a < n_index; a++) {
      int reverse = false;
      const auto* const ai_a = ai + a;
      found = false;
      if (!I->DiscreteFlag) { /* don't even try matching for discrete objects */
        while (b < n_atom) {
          ac = (fCompare(G, ai_a, atInfo + b));
          if (!ac) {
            found = true;
            break;
          } else if (ac < 0) { /* atom is before current, so try going the other
                                  way */
            reverse = true;
            break;
          } else if (!b) {
            int idx;
            ac = (fCompare(G, ai_a, atInfo + n_atom - 1));
            if (ac > 0) { /* atom is after all atoms in list, so don't even
                             bother searching */
              break;
            }
            /* next, try to jump to an appropriate position in the list */
            idx = ((a * n_atom) / n_index) - 1;
            if ((idx > 0) && (b != idx) && (idx < n_atom)) {
              ac = (fCompare(G, ai_a, atInfo + idx));
              if (ac > 0) b = idx;
            }
          }
          lb = b; /* last b before atom */
          b++;
        }
        if (reverse && !found) { /* searching going backwards... */
          while (b >= 0) {
            ac = (fCompare(G, ai_a, atInfo + b));
            if (!ac) {
              found = true;
              break;
            } else if (ac > 0) { /* atom after current -- no good */
              break;
            }
            lb = b;
            b--;
          }
          if (b < 0) b = 0;
        }
      }
      if (found) {
        index[a] = b; /* store real atom index b for a in index[a] */
        b++;
      } else {
        index[a] = I->NAtom + c; /* otherwise, this is a new atom */
        c++;
        b = lb;
      }
    }
  }

  /* first, reassign atom info for matched atoms */

  /* allocate additional space */
  if (ok) {
    if (c) {
      expansionFlag = true;
      nAt = I->NAtom + c;
    } else {
      nAt = I->NAtom;
    }
    if (expansionFlag) {
      VLACheck(I->AtomInfo, AtomInfoType, nAt);
      CHECKOK(ok, I->AtomInfo);
    }
  }

  if (ok) {
    for (a = 0; a < n_index; a++) { /* a is in original file space */
      a1 = outdex[cs->IdxToAtm[a]]; /* a1 is in sorted atom info space */
      a2 = index[a1];
      cs->IdxToAtm[a] = a2; /* a2 is in object space */
      if (a2 < oldNAtom)
        AtomInfoCombine(G, I->AtomInfo + a2, std::move(ai[a1]), aic_mask);
      else
        I->AtomInfo[a2] = std::move(ai[a1]);
    }
  }

  if (ok && I->DiscreteFlag) {
    ok = I->setNDiscrete(nAt);
  }

  I->NAtom = nAt;

  if (ok) {
    cs->updateNonDiscreteAtmToIdx(nAt);
  }

  VLAFreeP(
      ai); /* note that we're trusting AtomInfoCombine to have
              already purged all atom-dependent storage (such as strings) */

  AtomInfoFreeSortedIndexes(G, &index, &outdex);

  /* now find and integrate and any new bonds */
  if (ok && expansionFlag) { /* expansion flag means we have introduced at least
                                1 new atom */
    pymol::vla<BondType> bond;
    ok &= ObjectMoleculeConnect(I, nBond, bond, cs, bondSearchFlag, -1);
    if (nBond) {
      index = pymol::malloc<int>(nBond);
      CHECKOK(ok, index);
      c = 0;
      b = 0;
      nb = 0;
      if (ok) {
        for (a = 0; a < nBond; a++) { /* iterate over new bonds */
          found = false;
          if (!I->DiscreteFlag) { /* don't even try matching for discrete
                                     objects */
            b = nb;               /* pick up where we left off */
            while (b < I->NBond) {
              ac = BondCompare(bond + a, I->Bond + b);
              if (!ac) { /* zero is a match */
                found = true;
                break;
              } else if (ac < 0) { /* gone past position of this bond */
                break;
              } else if (!b) {
                int idx;
                ac = BondCompare(bond + a, I->Bond + I->NBond - 1);

                if (ac > 0) /* bond is after all bonds in list, so don't even
                               bother searching */
                  break;
                /* next, try to jump to an appropriate position in the list */
                idx = ((a * I->NBond) / nBond) - 1;
                if ((idx > 0) && (b != idx) && (idx < I->NBond)) {
                  ac = BondCompare(bond + a, I->Bond + idx);
                  if (ac > 0) b = idx;
                }
              }
              b++; /* no match yet, keep looking */
            }
          }
          if (found) {
            index[a] = b; /* existing bond... */
            nb = b + 1;
          } else { /* this is a new bond, save index and increment */
            index[a] = I->NBond + c;
            c++;
          }
        }
      }
      /* first, reassign atom info for matched atoms */
      if (c) {
        /* allocate additional space */
        nBd = I->NBond + c;

        if (!I->Bond) I->Bond = pymol::vla<BondType>(1);
        CHECKOK(ok, I->Bond);
        if (ok) VLACheck(I->Bond, BondType, nBd);
        CHECKOK(ok, I->Bond);
        for (a = 0; a < nBond; a++) { /* copy the new bonds */
          a2 = index[a];
          if (a2 >= I->NBond) {
            I->Bond[a2] = bond[a];
          }
        }
        I->NBond = nBd;
      }
      FreeP(index);
    }
  }
  if (invalidate) {
    if (oldNAtom) {
      if (oldNAtom == I->NAtom) {
        if (oldNBond != I->NBond) {
          I->invalidate(cRepAll, cRepInvBonds, -1);
        }
      } else {
        I->invalidate(cRepAll, cRepInvAtoms, -1);
      }
    }
  }
  return ok;
}
/**
 * @param state Object state or -2 for current state
 * @return nullptr if there is no CoordSet for the given state
 */
CoordSet* ObjectMolecule::getCoordSet(int state) {
  return static_cast<CoordSet*>(getObjectState(state));
}
const CoordSet* ObjectMolecule::getCoordSet(int state) const {
  return static_cast<const CoordSet*>(getObjectState(state));
}
void ObjectMoleculeTransformTTTf(ObjectMolecule* I, float* ttt, int frame) {
  int b;
  CoordSet* cs;
  for (b = 0; b < I->NCSet; b++) {
    if ((frame < 0) || (frame == b)) {
      cs = I->CSet[b];
      if (cs) {
        cs->invalidateRep(cRepAll, cRepInvCoord);
        MatrixTransformTTTfN3f(cs->NIndex, cs->Coord.data(), ttt,
                               cs->Coord.data());
        CoordSetRecordTxfApplied(cs, ttt, false);
      }
    }
  }
}
void ObjectMoleculeGetAtomSele(const ObjectMolecule* I, int index,
                               char* buffer) {
  PyMOLGlobals* G = I->G;
  assert(index < I->NAtom);
  auto* ai = I->AtomInfo + index;
  char inscode_str[2] = {ai->inscode, '\0'};

  snprintf(buffer, OrthoLineLength, "/%s/%s/%s/%s`%d%s/%s`%s", I->Name,
           LexStr(G, ai->segi), LexStr(G, ai->chain), LexStr(G, ai->resn),
           ai->resv, inscode_str, LexStr(G, ai->name), ai->alt);
}

/**
 * Get Atom selection string
 * @param I target Object Molecule
 * @param index atom index of I
 */

std::string ObjectMoleculeGetAtomSele(const ObjectMolecule* I, int index) {
  PyMOLGlobals* G = I->G;
  assert(index < I->NAtom);
  auto* ai = I->AtomInfo + index;
  char inscode_str[2] = {ai->inscode, '\0'};
  std::string buffer = pymol::string_format(
      "/%s/%s/%s/%s`%d%s/%s`%s", I->Name, LexStr(G, ai->segi),
      LexStr(G, ai->chain), LexStr(G, ai->resn), ai->resv, inscode_str,
      LexStr(G, ai->name), ai->alt);
  return buffer;
}
/**
 * Get Atom selection string
 * @param I target Object Molecule
 * @param index atom index of I
 */

std::string ObjectMolecule::describeElement(int index) const {
  auto I = this;
  auto buffer = ObjectMoleculeGetAtomSele(I, index);
  if (!I->AtomInfo[index].alt[0]) {
    // don't include the trailing backtick
    buffer.pop_back();
  }
  return buffer;
}
std::string ObjectMoleculeGetAtomSeleLog(const ObjectMolecule* I, int index,
                                         int quote) {
  OrthoLineType buffer;
  ObjectMoleculeGetAtomSeleLog(I, index, buffer, quote);
  return buffer;
}

void ObjectMoleculeGetAtomSeleLog(const ObjectMolecule* I, int index,
                                  char* buffer, int quote) {
  char* p = quote ? buffer + 1 : buffer;

  if (SettingGetGlobal_b(I->G, cSetting_robust_logs)) {
    ObjectMoleculeGetAtomSele(I, index, p);
  } else {
    sprintf(p, "(%s`%d)", I->Name, index + 1);
  }

  if (quote) {
    int len = strlen(p);
    buffer[0] = buffer[len + 1] = '"';
    buffer[len + 2] = 0;
  }
}

std::string ObjectMoleculeGetAtomSeleFast(const ObjectMolecule* I, int index) {
  auto G = I->G;
  auto* ai = I->AtomInfo + index;
  char inscodestr[2] = {ai->inscode, 0};
  return pymol::string_format("(/'%s'/'%s'/'%s'/'%s'`%d%s/'%s'`'%s')", I->Name,
                              LexStr(G, ai->segi), LexStr(G, ai->chain),
                              LexStr(G, ai->resn), ai->resv, inscodestr,
                              ai->name, ai->alt);
}
int ObjectMolecule::getNFrame() const { return NCSet; }

struct _CCoordSetUpdateThreadInfo {
  CoordSet* cs;
  int a;
};

void CoordSetUpdateThread(CCoordSetUpdateThreadInfo* T) {
  if (T->cs) {
    T->cs->update(T->a);
  }
}

static void ObjMolCoordSetUpdateSpawn(PyMOLGlobals* G,
                                      CCoordSetUpdateThreadInfo* Thread,
                                      int n_thread, int n_total) {
  if (n_total == 1) {
    CoordSetUpdateThread(Thread);
  } else if (n_total) {
    int blocked;
    PyObject* info_list;
    int a, n = 0;
    blocked = PAutoBlock(G);

    PRINTFB(G, FB_Scene, FB_Blather)
    " Scene: updating coordinate sets with %d threads...\n", n_thread ENDFB(G);
    info_list = PyList_New(n_total);
    for (a = 0; a < n_total; a++) {
      PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
      n++;
    }
    PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_coordset_update_spawn", "Oi",
                                 info_list, n_thread));
    Py_DECREF(info_list);
    PAutoUnblock(G, blocked);
  }
}
void ObjectMolecule::update() {
  auto I = this;
  int a; /*, ok; */

  OrthoBusyPrime(G);
  /* if the cached representation is invalid, reset state */
  if (!I->RepVisCacheValid) {
    /* note which representations are active */
    /* for each atom in each coordset, blank out the representation cache */
    if (I->NCSet > 1) {
      const AtomInfoType* ai = I->AtomInfo.data();
      I->RepVisCache = 0;
      for (a = 0; a < I->NAtom; a++) {
        I->RepVisCache |= ai->visRep;
        ai++;
      }
    } else {
      I->RepVisCache = cRepBitmask; /* if only one coordinate set, then
                                     * there's no benefit to pre-filtering
                                     * the representations... */
    }
    I->RepVisCacheValid = true;
  }
  {
    /* determine the start/stop states */
    int start = 0;
    int stop = I->NCSet;
    /* set start and stop given an object */
    ObjectAdjustStateRebuildRange(I, &start, &stop);
    if ((I->NCSet == 1) && (SettingGet_b(G, I->Setting.get(), nullptr,
                                         cSetting_static_singletons))) {
      start = 0;
      stop = 1;
    }
    if (stop > I->NCSet) stop = I->NCSet;

    /* single and multithreaded coord set updates */
    {
      int n_thread = SettingGetGlobal_i(G, cSetting_max_threads);
      int multithread = SettingGetGlobal_i(G, cSetting_async_builds);

      if (multithread && (n_thread) && (stop - start) > 1) {
        int cnt = 0;

        /* must precalculate to avoid race-condition since this isn't
           mutexed yet and neighbors are needed by cartoons */
        this->getNeighborArray();

        for (a = start; a < stop; a++)
          if ((a < I->NCSet) && I->CSet[a]) cnt++;
        {
          CCoordSetUpdateThreadInfo* thread_info =
              pymol::malloc<CCoordSetUpdateThreadInfo>(cnt);
          if (thread_info) {
            cnt = 0;
            for (a = start; a < stop; a++) {
              if ((a < I->NCSet) && I->CSet[a]) {
                thread_info[cnt].cs = I->CSet[a];
                thread_info[cnt].a = a;
                cnt++;
              }
            }
            ObjMolCoordSetUpdateSpawn(G, thread_info, n_thread, cnt);
            FreeP(thread_info);
          }
        }

      } else { /* single thread */
        for (a = start; a < stop; a++) {
          if ((a < I->NCSet) && I->CSet[a] && (!G->Interrupt)) {
            /* status bar */
            OrthoBusySlow(G, a, I->NCSet);
            PRINTFB(G, FB_ObjectMolecule, FB_Blather)
            " ObjectMolecule-DEBUG: updating representations for state %d of "
            "\"%s\".\n",
                a + 1, I->Name ENDFB(G);
            I->CSet[a]->update(a);
          }
        }
      }
    }
  } /* end block */

  PRINTFD(G, FB_ObjectMolecule)
  " ObjectMolecule: updates complete for object %s.\n", I->Name ENDFD;
}
void ObjectMolecule::invalidate(cRep_t rep, cRepInv_t level, int state) {
  auto I = this;
  int a;
  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: entered. rep: %d level: %d\n", __func__, rep, level ENDFD;

  auto const level_actual = level;

  // Remove the "purge" bit
  level = static_cast<decltype(level)>(level & ~cRepInvPurgeMask);

  if (level >= cRepInvVisib) {
    I->RepVisCacheValid = false;
  }

  if (level >= cRepInvBondsNoNonbonded) {
    if (level < cRepInvBonds) {
      level = cRepInvBonds;
    } else {
      ObjectMoleculeUpdateNonbonded(I);
    }
  }

  if (level >= cRepInvBonds) {
    this->Neighbor.reset();
    if (I->Sculpt) {
      DeleteP(I->Sculpt);
    }
    if (level >= cRepInvAtoms) {
      SelectorUpdateObjectSele(I->G, I);
    }
  }
  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: invalidating representations...\n", __func__ ENDFD;

  if (level >= cRepInvColor) {
    /* after label, this gets called, so we shouldn't invalidate types b/c
       PYMOL-317 not sure if that is the exact level we should cut this off at
       1/28/13 BB: changed from >cRepInvText to >=cRepInvColor for invalidating
       surface colors during gadget invalidation */

    int start = 0;
    int stop = I->NCSet;

    if (state >= 0) {
      start = state;
      stop = state + 1;
    }
    if (stop > I->NCSet) stop = I->NCSet;
    for (a = start; a < stop; a++) {
      CoordSet* cset = 0;
      cset = I->CSet[a];
      if (cset) {
        cset->invalidateRep(rep, level_actual);
      }
    }
  }

  PRINTFD(I->G, FB_ObjectMolecule)
  " %s: leaving...\n", __func__ ENDFD;
}

int ObjectMoleculeMoveAtom(ObjectMolecule* I, int state, int index,
                           const float* v, int mode, int log) {
  int result = 0;
  PyMOLGlobals* G = I->G;
  CoordSet* cs;
  if (!(I->AtomInfo[index].protekted == cAtomProtected_explicit)) {
    if (state < 0) state = 0;
    if (I->NCSet == 1) state = 0;
    state = state % I->NCSet;
    if ((!I->CSet[state]) &&
        (SettingGet_b(G, I->Setting.get(), nullptr, cSetting_all_states)))
      state = 0;
    cs = I->CSet[state];
    if (cs) {
      result = CoordSetMoveAtom(I->CSet[state], index, v, mode);
      cs->invalidateRep(cRepAll, cRepInvCoord);
      ExecutiveUpdateCoordDepends(G, I);
    }
  }
  if (log) {
    OrthoLineType line, buffer;
    if (SettingGetGlobal_i(G, cSetting_logging)) {
      ObjectMoleculeGetAtomSele(I, index, buffer);
      sprintf(line,
              "cmd.translate_atom(\"%s\",%15.9f,%15.9f,%15.9f,%d,%d,%d)\n",
              buffer, v[0], v[1], v[2], state + 1, mode, 0);
      PLog(G, line, cPLog_no_flush);
    }
  }
  return (result);
}
int ObjectMoleculeMoveAtomLabel(ObjectMolecule* I, int state, int index,
                                float* v, int log, float* diff) {
  int result = 0;
  CoordSet* cs;
  if (!(I->AtomInfo[index].protekted == cAtomProtected_explicit)) {
    if (state < 0) state = 0;
    if (I->NCSet == 1) state = 0;
    state = state % I->NCSet;
    if ((!I->CSet[state]) &&
        (SettingGet_b(I->G, I->Setting.get(), nullptr, cSetting_all_states)))
      state = 0;
    cs = I->CSet[state];
    if (cs) {
      result = CoordSetMoveAtomLabel(I->CSet[state], index, v, diff);
      cs->invalidateRep(cRepLabel, cRepInvCoord);
    }
  }
  return (result);
}
int ObjectMoleculeInitBondPath(ObjectMolecule* I, ObjectMoleculeBPRec* bp) {
  int a;
  bp->dist = pymol::malloc<int>(I->NAtom);
  bp->list = pymol::malloc<int>(I->NAtom);
  for (a = 0; a < I->NAtom; a++) bp->dist[a] = -1;
  bp->n_atom = 0;
  return 1;
}
int ObjectMoleculePurgeBondPath(ObjectMolecule* I, ObjectMoleculeBPRec* bp) {
  FreeP(bp->dist);
  FreeP(bp->list);
  return 1;
}
int ObjectMoleculeGetBondPaths(ObjectMolecule* I, int atom, int max,
                               ObjectMoleculeBPRec* bp) {
  /* returns list of bond counts from atom to all others
     dist and list must be vla array pointers or nullptr */

  int b_cnt = 0;

  /* reinitialize dist array (if we've done at least one pass) */

  for (int a = 0; a < bp->n_atom; a++) bp->dist[bp->list[a]] = -1;

  bp->n_atom = 0;
  bp->dist[atom] = 0;
  bp->list[bp->n_atom] = atom;
  bp->n_atom++;

  int cur = 0;
  while (1) {
    b_cnt++;
    if (b_cnt > max) break;

    unsigned n_cur = bp->n_atom - cur;

    /* iterate through all current atoms */

    if (!n_cur) break;
    while (n_cur--) {
      int const a1 = bp->list[cur++];
      for (auto const& neighbor : AtomNeighbors(I, a1)) {
        auto const a2 = neighbor.atm;
        if (bp->dist[a2] < 0) { /* for each atom not yet sampled... */
          bp->dist[a2] = b_cnt;
          bp->list[bp->n_atom] = a2;
          bp->n_atom++;
        }
      }
    }
  }
  return (bp->n_atom);
}
int*** ObjectMoleculeGetBondPrint(ObjectMolecule* I, int max_bond, int max_type,
                                  int* dim) {
  int a, b, i, c;
  int at1, at2;
  int*** result = nullptr;
  ObjectMoleculeBPRec bp;

  dim[0] = max_type + 1;
  dim[1] = max_type + 1;
  dim[2] = max_bond + 1;

  result = (int***)UtilArrayCalloc((unsigned int*)dim, 3, sizeof(int));

  ObjectMoleculeInitBondPath(I, &bp);
  for (a = 0; a < I->NAtom; a++) {
    at1 = I->AtomInfo[a].customType;
    if ((at1 >= 0) && (at1 <= max_type)) {
      ObjectMoleculeGetBondPaths(I, a, max_bond, &bp);
      for (b = 0; b < bp.n_atom; b++) {
        i = bp.list[b];
        at2 = I->AtomInfo[i].customType;
        if ((at2 >= 0) && (at2 <= max_type)) {
          c = bp.dist[i];
          result[at1][at2][c]++;
        }
      }
    }
  }
  ObjectMoleculePurgeBondPath(I, &bp);
  return (result);
}
int ObjectMoleculeGetAtomVertex(const ObjectMolecule* I, int state, int index,
                                float* v) {
  int result = 0;
  if (state < 0)
    state = SettingGet_i(I->G, nullptr, I->Setting.get(), cSetting_state) - 1;
  if (state < 0) state = SceneGetState(I->G);
  if (I->NCSet == 1)
    state = 0; /* static singletons always active here it seems */
  state = state % I->NCSet;
  if ((!I->CSet[state]) &&
      (SettingGet_b(I->G, I->Setting.get(), nullptr, cSetting_all_states)))
    state = 0;
  if (I->CSet[state]) result = CoordSetGetAtomVertex(I->CSet[state], index, v);

  return (result);
}
int ObjectMoleculeGetAtomTxfVertex(const ObjectMolecule* I, int state,
                                   int index, float* v) {
  int result = 0;
  const CoordSet* cs = nullptr;
  if (I->DiscreteFlag) {
    cs = I->DiscreteCSet[index];
  }
  if (state < 0) {
    state = SettingGet_i(I->G, nullptr, I->Setting.get(), cSetting_state) - 1;
  }
  if (state < 0) state = SceneGetState(I->G);
  if (I->NCSet == 1)
    state = 0; /* static singletons always active here it seems */
  state = state % I->NCSet;
  {
    if (!cs) cs = I->CSet[state];
    if ((!cs) &&
        (SettingGet_b(I->G, I->Setting.get(), nullptr, cSetting_all_states))) {
      state = 0;
      cs = I->CSet[state];
    }
    if (cs) {
      result = CoordSetGetAtomTxfVertex(cs, index, v);
    }
  }
  return (result);
}
int ObjectMoleculeSetAtomVertex(ObjectMolecule* I, int state, int index,
                                float* v) {
  int result = 0;
  if (state < 0)
    state = SettingGet_i(I->G, nullptr, I->Setting.get(), cSetting_state) - 1;
  if (state < 0) state = SceneGetState(I->G);
  if (I->NCSet == 1) state = 0;
  state = state % I->NCSet;
  if ((!I->CSet[state]) &&
      (SettingGet_b(I->G, I->Setting.get(), nullptr, cSetting_all_states)))
    state = 0;
  if (I->CSet[state]) result = CoordSetSetAtomVertex(I->CSet[state], index, v);
  return (result);
}
void ObjectMolecule::render(RenderInfo* info) {
  auto I = this;
  int state = info->state;
  const RenderPass pass = info->pass;
  CoordSet* cs;
  int pop_matrix = false;
  int use_matrices =
      SettingGet_i(I->G, I->Setting.get(), nullptr, cSetting_matrix_mode);
  if (use_matrices < 0) use_matrices = 0;
  PRINTFD(I->G, FB_ObjectMolecule)
  " ObjectMolecule: rendering %s pass %d...\n", I->Name,
      static_cast<int>(pass) ENDFD;

  ObjectPrepareContext(I, info);

  for (StateIterator iter(G, I->Setting.get(), state, I->NCSet); iter.next();) {
    cs = I->CSet[iter.state];
    if (cs) {
      if (use_matrices) pop_matrix = ObjectStatePushAndApplyMatrix(cs, info);
      cs->render(info);
      if (pop_matrix) ObjectStatePopMatrix(cs, info);
    }
  }
  PRINTFD(I->G, FB_ObjectMolecule)
  " ObjectMolecule: rendering complete for object %s.\n", I->Name ENDFD;
}
void ObjectMoleculeDummyUpdate(ObjectMolecule* I, int mode) {
  switch (mode) {
    case cObjectMoleculeDummyOrigin:
      SceneOriginGet(I->G, I->CSet[0]->Coord.data());
      break;
    case cObjectMoleculeDummyCenter:
      SceneGetCenter(I->G, I->CSet[0]->Coord.data());
      break;
  }
}
ObjectMolecule* ObjectMoleculeDummyNew(PyMOLGlobals* G, int type) {
  auto I = new ObjectMolecule(G, false);
  I->AtomInfo.resize(I->NAtom = 1);

  auto* cset = new CoordSet(G);
  cset->Obj = I;
  cset->setNIndex(1);
  cset->enumIndices();

  I->CSet.resize(I->NCSet = 1);
  I->CSet[0] = cset;

  return I;
}
ObjectMolecule::ObjectMolecule(PyMOLGlobals* G, int discreteFlag)
    : pymol::CObject(G) {
  auto I = this;
  int a;
  I->type = cObjectMolecule;
  I->CSet = pymol::vla<CoordSet*>(10); /* auto-zero */
  I->DiscreteFlag = discreteFlag;
  if (I->DiscreteFlag) { /* discrete objects don't share atoms between states */
    I->DiscreteAtmToIdx = pymol::vla<int>(0);
    I->DiscreteCSet = pymol::vla<CoordSet*>(0);
  } else {
    I->DiscreteAtmToIdx = nullptr;
    I->DiscreteCSet = nullptr;
  }
  I->AtomInfo = pymol::vla<AtomInfoType>(10);
  for (a = 0; a <= cUndoMask; a++) {
    I->UndoCoord[a] = nullptr;
    I->UndoState[a] = -1;
  }
  I->UndoIter = 0;
}
void ObjectMoleculeCopyNoAlloc(const ObjectMolecule* obj, ObjectMolecule* I) {
  PyMOLGlobals* G = const_cast<PyMOLGlobals*>(obj->G);

  int a;
  BondType* i0;
  const BondType* i1;
  (*I) = (*obj);
  I->Sculpt = nullptr;
  I->Setting.reset(SettingCopyAll(G, obj->Setting.get(), nullptr));

  I->ViewElem = nullptr;
  I->gridSlotSelIndicatorsCGO = nullptr;

  for (a = 0; a <= cUndoMask; a++) I->UndoCoord[a] = nullptr;
  I->CSet = pymol::vla<CoordSet*>(I->NCSet); /* auto-zero */
  for (a = 0; a < I->NCSet; a++) {
    I->CSet[a] = CoordSetCopy(obj->CSet[a]);
    if (I->CSet[a]) I->CSet[a]->Obj = I;
  }

  if (obj->CSTmpl) I->CSTmpl = CoordSetCopy(obj->CSTmpl);

  if (obj->DiscreteFlag) {
    int sz = VLAGetSize(obj->DiscreteAtmToIdx);
    I->DiscreteAtmToIdx = VLACopy2(obj->DiscreteAtmToIdx);
    I->DiscreteCSet = pymol::vla<CoordSet*>(sz);
    I->updateAtmToIdx();
  }
  I->Bond = pymol::vla<BondType>(I->NBond);
  i0 = I->Bond.data();
  i1 = obj->Bond.data();
  for (a = 0; a < I->NBond; a++) {
    AtomInfoBondCopy(G, i1++, i0++);
  }

  // unfortunately, the AtomInfoType is not copy-constructable yet
  // assert copy done correct
  if (I->AtomInfo.size() != obj->AtomInfo.size()) {
    throw "AtomInfo copy failed";
  }
  AtomInfoType* a0 = I->AtomInfo.data();
  const AtomInfoType* a1 = obj->AtomInfo.data();
  memset(a0, 0, sizeof(AtomInfoType) * I->NAtom);
  for (a = 0; a < I->NAtom; a++) AtomInfoCopy(G, a1++, a0++);
}

ObjectMolecule* ObjectMoleculeCopy(const ObjectMolecule* obj) {
  PyMOLGlobals* G = const_cast<PyMOLGlobals*>(obj->G);
  auto I = new ObjectMolecule(G, obj->DiscreteFlag);
  ObjectMoleculeCopyNoAlloc(obj, I);
  return (I);
}
/**
 * Set the order of coordinate sets with an index array
 */
int ObjectMoleculeSetStateOrder(ObjectMolecule* I, int* order, int len) {
  int a;
  CoordSet** csets = VLAlloc(CoordSet*, I->NCSet);

  ok_assert(1, len == I->NCSet);

  // invalidate
  I->invalidate(cRepAll, cRepInvAll, -1);

  // new coord set array
  for (a = 0; a < I->NCSet; a++) {
    int i = order[a];
    ok_assert(1, 0 <= i && i < I->NCSet);
    csets[a] = I->CSet[i];
  }

  VLAFreeP(I->CSet);
  I->CSet = pymol::vla_take_ownership(csets);

  return true;
ok_except1:
  ErrMessage(I->G, "ObjectMoleculeSetStateOrder", "failed");
  VLAFreeP(csets);
  return false;
}
pymol::Result<> ObjectMoleculeDeleteStates(ObjectMolecule* I,
                                           const std::vector<int>& states) {
  auto& CSet = I->CSet;

  // first pass, check for invalid states
  for (auto state : states) {
    if (state < 0 || state >= I->NCSet) {
      auto msg = pymol::string_format("Invalid state index: %d", state);
      I->G->Feedback->addColored(msg.c_str(), FB_Errors);
      return {};
    }
  }

  // second pass, delete states
  for (auto it = states.rbegin(); it != states.rend(); ++it) {
    int state = *it;
    DeleteP(CSet[state]);
    CSet.erase(state, 1);
  }
  I->NCSet -= states.size();
  CSet.resize(I->NCSet);

  // Update Rep States
  for (int i = 0; i < I->NCSet; ++i) {
    if (CSet[i]) {
      auto& cset = *CSet[i];
      for (int r = 0; r < cRepCnt; ++r) {
        if (cset.Rep[r] && cset.Rep[r]->context.state != 0) {
          cset.Rep[r]->context.state = i;
        }
      }
    }
  }

  return {};
}
ObjectMolecule::~ObjectMolecule() {
  auto I = this;
  int a;
  SelectorPurgeObjectMembers(I->G, I);
  for (a = 0; a < I->NCSet; a++) {
    if (I->CSet[a]) {
      delete I->CSet[a];
      I->CSet[a] = nullptr;
    }
  }
  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  VLAFreeP(I->CSet);
  I->m_ciffile.reset();  // free data

  {
    int nAtom = I->NAtom;
    AtomInfoType* ai = I->AtomInfo.data();

    for (a = 0; a < nAtom; a++) {
      AtomInfoPurge(I->G, ai);
      ai++;
    }
    VLAFreeP(I->AtomInfo);
  }
  {
    int nBond = I->NBond;
    BondType* bi = I->Bond.data();

    for (a = 0; a < nBond; a++) {
      AtomInfoPurgeBond(I->G, bi);
      bi++;
    }
    VLAFreeP(I->Bond);
  }
  for (a = 0; a <= cUndoMask; a++) FreeP(I->UndoCoord[a]);
  if (I->Sculpt) DeleteP(I->Sculpt);
  delete I->CSTmpl;
}

CPythonVal* ObjectMoleculeGetProperty(ObjectMolecule* I, const char* propname,
                                      int state, int quiet) {
  PyMOLGlobals* G = I->G;
  for (StateIterator iter(G, I->Setting.get(), state, I->NCSet); iter.next();) {
    if (I->CSet[iter.state])
      return PropertyGetPyObject(G, I->CSet[iter.state], propname);
  }
  return nullptr;
}

int ObjectMoleculeSetProperty(ObjectMolecule* I, const char* propname,
                              CPythonVal* value, const PropertyType& proptype,
                              int state, int quiet) {
  int ok = true;
  PyMOLGlobals* G = I->G;
  for (StateIterator iter(G, I->Setting.get(), state, I->NCSet); iter.next();) {
    if (I->CSet[iter.state])
      ok &= PropertySet(G, I->CSet[iter.state], propname, value, proptype);
  }
  return (ok);
}

void AtomInfoSettingGenerateSideEffects(PyMOLGlobals* G, ObjectMolecule* obj,
                                        int index, int id) {
  switch (index) {
    case cSetting_label_position:
    case cSetting_label_placement_offset:
    case cSetting_label_screen_point:
    case cSetting_label_relative_mode:
      obj->invalidate(cRepLabel, cRepInvCoord, -1);
  }
}

static int AtomInfoInOrder(PyMOLGlobals* G, const AtomInfoType* atom, int atom1,
                           int atom2) {
  return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrderIgnoreHet(PyMOLGlobals* G, const AtomInfoType* atom,
                                    int atom1, int atom2) {
  return (AtomInfoCompareIgnoreHet(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrigOrder(PyMOLGlobals* G, const AtomInfoType* atom,
                               int atom1, int atom2) {
  if (atom[atom1].rank == atom[atom2].rank)
    return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
  return (atom[atom1].rank < atom[atom2].rank);
}

int* AtomInfoGetSortedIndex(PyMOLGlobals* G, const ObjectMolecule* obj,
                            const AtomInfoType* rec, int n, int** outdex) {
  int* index;
  int a;
  const CSetting* setting = nullptr;

  ok_assert(1, index = pymol::malloc<int>(n + 1));
  ok_assert(1, (*outdex) = pymol::malloc<int>(n + 1));

  if (obj && obj->DiscreteFlag) {
    for (a = 0; a < n; a++) index[a] = a;
  } else {
    if (obj) setting = obj->Setting.get();

    UtilSortIndexGlobals(
        G, n, rec, index,
        (UtilOrderFnGlobals*)(SettingGet_b(G, setting, nullptr,
                                           cSetting_retain_order)
                                  ? AtomInfoInOrigOrder
                              : SettingGet_b(G, setting, nullptr,
                                             cSetting_pdb_hetatm_sort)
                                  ? AtomInfoInOrder
                                  : AtomInfoInOrderIgnoreHet));
  }

  for (a = 0; a < n; a++) (*outdex)[index[a]] = a;

  return index;

ok_except1:
  FreeP(index);
  return nullptr;
}

/**
 * Set the size of the DiscreteAtmToIdx and DiscreteCSet VLAs and pad
 * them with -1/NULL if necessary.
 */
bool ObjectMolecule::setNDiscrete(int natom) {
  int n = VLAGetSize(DiscreteAtmToIdx);

  if (n == natom) return true;

  VLASize(DiscreteAtmToIdx, int, natom);
  VLASize(DiscreteCSet, CoordSet*, natom);

  if (!DiscreteAtmToIdx || !DiscreteCSet) return false;

  for (int i = n; i < natom; ++i) {
    DiscreteAtmToIdx[i] = -1;
    DiscreteCSet[i] = nullptr;
  }

  return true;
}

/**
 * Update the AtmToIdx or DiscreteAtmToIdx/DiscreteCSet VLAs from the
 * IdxToAtm arrays
 */
bool ObjectMolecule::updateAtmToIdx() {
  if (DiscreteFlag) {
    ok_assert(1, setNDiscrete(NAtom));
  }

  for (int i = -1; i < NCSet; ++i) {
    CoordSet* cset = (i < 0) ? CSTmpl : CSet[i];

    if (!cset) continue;

    if (!DiscreteFlag) {
      cset->updateNonDiscreteAtmToIdx(NAtom);
    } else {
      for (int idx = 0; idx < cset->NIndex; ++idx) {
        int atm = cset->IdxToAtm[idx];
        assert(atm < NAtom);
        DiscreteAtmToIdx[atm] = idx;
        DiscreteCSet[atm] = cset;
        AtomInfo[atm].discrete_state = i + 1;
      }
    }
  }

  return true;
ok_except1:
  return false;
}

/**
 * Check if this atom has coordinates in any state
 */
bool ObjectMolecule::atomHasAnyCoordinates(size_t atm) const {
  for (size_t i = 0; i < NCSet; ++i) {
    auto cset = CSet[i];
    if (cset && cset->atmToIdx(atm) != -1) {
      return true;
    }
  }

  return false;
}

pymol::CObject* ObjectMolecule::clone() const {
  return ObjectMoleculeCopy(this);
}

AtomNeighbors::AtomNeighbors(const ObjectMolecule* I, int atm) {
  auto* Neighbor = I->getNeighborArray();

  m_neighbor = Neighbor + Neighbor[atm];
}
