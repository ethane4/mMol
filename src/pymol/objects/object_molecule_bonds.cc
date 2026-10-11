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

// Bonding atoms by distance, testing bonds, and finding hydrogen bonds.

#include "pymol/objects/object_molecule.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util.h"
#include "pymol/objects/object_molecule_internal.h"

/* Check if atom is bonded to an atom with given name
 *
 * param same_res:
 *    0 = must not be in same residue
 *    1 = must be in same residue
 *   -1 = don't check residue
 */
int ObjectMoleculeIsAtomBondedToName(ObjectMolecule* obj, int a0,
                                     const char* name, int same_res) {
  PyMOLGlobals* G = obj->G;
  AtomInfoType const* const ai0 = obj->AtomInfo.data() + a0;

  if (a0 >= 0) {
    for (auto const& neighbor : AtomNeighbors(obj, a0)) {
      auto const* const ai2 = obj->AtomInfo.data() + neighbor.atm;
      if (WordMatchExact(G, LexStr(G, ai2->name), name, true) &&
          (same_res < 0 || (same_res == AtomInfoSameResidue(G, ai0, ai2))))
        return true;
    }
  }
  return false;
}

int ObjectMoleculeAreAtomsBonded2(ObjectMolecule* obj0, int a0,
                                  ObjectMolecule* obj1, int a1) {
  if (obj0 == obj1 && a0 >= 0) {
    assert(a1 >= 0);
    for (auto const& neighbor : AtomNeighbors(obj0, a0)) {
      if (a1 == neighbor.atm) return true;
    }
  }
  return false;
}

bool ObjectMoleculeIsAtomBondedToSele(const ObjectMolecule* I, int atm,
                                      SelectorID_t sele) {
  if (atm < I->NAtom) {
    for (auto const& neighbor : AtomNeighbors(I, atm)) {
      if (SelectorIsMember(I->G, I->AtomInfo[neighbor.atm].selEntry, sele)) {
        return true;
      }
    }
  }
  return false;
}
void ObjectMoleculeInitHBondCriteria(PyMOLGlobals* G, HBondCriteria* hbc) {
  hbc->maxAngle = SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_max_angle);
  hbc->maxDistAtMaxAngle =
      SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_cutoff_edge);
  hbc->maxDistAtZero =
      SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_cutoff_center);
  hbc->power_a = SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_power_a);
  hbc->power_b = SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_power_b);
  hbc->cone_dangle = (float)cos(
      PI * 0.5 * SettingGet_f(G, nullptr, nullptr, cSetting_h_bond_cone) /
      180.0F);
  if (hbc->maxDistAtMaxAngle != 0.0F) {
    hbc->factor_a = (float)(0.5 / pow(hbc->maxAngle, hbc->power_a));
    hbc->factor_b = (float)(0.5 / pow(hbc->maxAngle, hbc->power_b));
  }
}

static int ObjectMoleculeTestHBond(float* donToAcc, float* donToH,
                                   float* hToAcc, float* accPlane,
                                   HBondCriteria* hbc) {
  float nDonToAcc[3], nDonToH[3], nAccPlane[3], nHToAcc[3];
  double angle;
  double cutoff;
  double curve;
  double dist;
  float dangle;

  /* A ~~ H - D */

  normalize23f(donToAcc, nDonToAcc);
  normalize23f(hToAcc, nHToAcc);
  if (accPlane) { /* remember, plane need not exist if it's water... */
    normalize23f(accPlane, nAccPlane);
    if (dot_product3f(nHToAcc, nAccPlane) >
        (-hbc->cone_dangle)) /* don't allow H behind Acceptor plane */
      return 0;
  }

  normalize23f(donToH, nDonToH);
  normalize23f(donToAcc, nDonToAcc);

  dangle = dot_product3f(nDonToH, nDonToAcc);
  if ((dangle < 1.0F) && (dangle > 0.0F))
    angle = 180.0 * acos((double)dangle) / PI;
  else if (dangle > 0.0F)
    angle = 0.0;
  else
    angle = 90.0;

  if (angle > hbc->maxAngle) return 0;

  /* interpolate cutoff based on ADH angle */

  if (hbc->maxDistAtMaxAngle != 0.0F) {
    curve = (pow(angle, (double)hbc->power_a) * hbc->factor_a +
             pow(angle, (double)hbc->power_b) * hbc->factor_b);

    cutoff =
        (hbc->maxDistAtMaxAngle * curve) + (hbc->maxDistAtZero * (1.0 - curve));
  } else {
    cutoff = hbc->maxDistAtZero;
  }

  /*
     printf("angle %8.3f curve %8.3f %8.3f %8.3f %8.3f\n",angle,
     curve,cutoff,hbc->maxDistAtMaxAngle,hbc->maxDistAtZero);
   */

  dist = length3f(donToAcc);

  if (dist > cutoff)
    return 0;
  else
    return 1;
}
static int ObjectMoleculeFindBestDonorH(ObjectMolecule* I, int atom, int state,
                                        float* dir, float* best,
                                        AtomInfoType** h_real) {
  int result = 0;
  CoordSet* cs;
  float cand[3], cand_dir[3];
  float best_dot = 0.0F, cand_dot;

  if ((state >= 0) && (state < I->NCSet) && (cs = I->CSet[state]) &&
      (atom < I->NAtom)) {
    auto idx = cs->atmToIdx(atom);

    if (idx >= 0) {
      const float* orig = cs->coordPtr(idx);

      /*  do we need to add any new hydrogens? */

      auto const neighbors = AtomNeighbors(I, atom);
      int const nn = neighbors.size();

      if ((nn < I->AtomInfo[atom].valence) ||
          I->AtomInfo[atom].hb_donor) { /* is there an implicit hydrogen? */
        if (CoordSetFindOpenValenceVector(cs, atom, best, dir)) {
          result = true;
          best_dot = dot_product3f(best, dir);
          add3f(orig, best, best);
          if (h_real) *h_real = nullptr;
        }
      }
      /* iterate through real hydrogens looking for best match
         with desired direction */

      /* look for an attached non-hydrogen as a base */
      for (int i = 0; i < nn; ++i) {
        int const a1 = neighbors[i].atm;
        if (I->AtomInfo[a1].protons == 1) {                      /* hydrogen */
          if (ObjectMoleculeGetAtomVertex(I, state, a1, cand)) { /* present */

            subtract3f(cand, orig, cand_dir);
            normalize3f(cand_dir);
            cand_dot = dot_product3f(cand_dir, dir);
            if (result) { /* improved */
              if ((best_dot < cand_dot) || (h_real && !*h_real)) {
                best_dot = cand_dot;
                copy3f(cand, best);
                if (h_real) *h_real = I->AtomInfo + a1;
              }
            } else { /* first */
              result = true;
              copy3f(cand, best);
              best_dot = cand_dot;
              if (h_real) *h_real = I->AtomInfo + a1;
            }
          }
        }
      }
    }
  }
  return result;
}
int ObjectMoleculeGetCheckHBond(AtomInfoType** h_real, float* h_crd_ret,
                                ObjectMolecule* don_obj, int don_atom,
                                int don_state, ObjectMolecule* acc_obj,
                                int acc_atom, int acc_state,
                                HBondCriteria* hbc) {
  int result = 0;
  const CoordSet *csD, *csA;
  float donToAcc[3];
  float donToH[3];
  float bestH[3];
  float hToAcc[3];
  float accPlane[3], *vAccPlane = nullptr;

  /* first, check for existence of coordinate sets */

  if ((don_state >= 0) && (don_state < don_obj->NCSet) &&
      (csD = don_obj->CSet[don_state]) && (acc_state >= 0) &&
      (acc_state < acc_obj->NCSet) && (csA = acc_obj->CSet[acc_state]) &&
      (don_atom < don_obj->NAtom) && (acc_atom < acc_obj->NAtom)) {
    /* now check for coordinates of these actual atoms */

    auto idxD = csD->atmToIdx(don_atom);
    auto idxA = csA->atmToIdx(acc_atom);

    if ((idxA >= 0) && (idxD >= 0)) {
      /* now get local geometries, including
         real or virtual hydrogen atom positions */

      const float* vDon = csD->coordPtr(idxD);
      const float* vAcc = csA->coordPtr(idxA);

      subtract3f(vAcc, vDon, donToAcc);

      if (ObjectMoleculeFindBestDonorH(don_obj, don_atom, don_state, donToAcc,
                                       bestH, h_real)) {
        subtract3f(bestH, vDon, donToH);
        subtract3f(vAcc, bestH, hToAcc);

        if (ObjectMoleculeGetAvgHBondVector(acc_obj, acc_atom, acc_state,
                                            accPlane, hToAcc) > 0.1) {
          vAccPlane = &accPlane[0];
        }
        result =
            ObjectMoleculeTestHBond(donToAcc, donToH, hToAcc, vAccPlane, hbc);
        if (result && h_crd_ret && h_real && *h_real) copy3f(bestH, h_crd_ret);
      }
    }
  }

  return (result);
}
static float connect_cutoff_adjustment(const AtomInfoType* ai1,
                                       const AtomInfoType* ai2) {
  if (ai1->isHydrogen() || ai2->isHydrogen()) return -0.2f;

  if (ai1->protons == cAN_S || ai2->protons == cAN_S) return 0.2f;

  return 0.f;
}

/**
 * True if two atoms should be bonded
 */
static bool is_distance_bonded(PyMOLGlobals* G, const CoordSet* cs,
                               const AtomInfoType* ai1, const AtomInfoType* ai2,
                               const float* v1, const float* v2, float cutoff,
                               int connect_mode, int discrete_chains,
                               bool connect_bonded, bool unbond_cations) {
  auto dst = (float)diff3f(v1, v2);

  if (dst < R_SMALL4) return false;

  dst -= (ai1->vdw + ai2->vdw) / 2;

  cutoff += connect_cutoff_adjustment(ai1, ai2);

  if (dst > cutoff) return false;

  if (ai1->isHydrogen() && ai2->isHydrogen()) return false;

  if (discrete_chains > 0 && ai1->chain != ai2->chain) return false;

  if (!connect_bonded && ai1->bonded && ai2->bonded) return false;

  bool water_flag = (AtomInfoKnownWaterResName(G, LexStr(G, ai1->resn)) ||
                     AtomInfoKnownWaterResName(G, LexStr(G, ai2->resn)));

  if (connect_mode != 3 &&
      cs->TmpBond && /* connectivity information present in file */
      ai1->hetatm && ai2->hetatm && !water_flag &&
      !(AtomInfoKnownPolymerResName(LexStr(G, ai1->resn)) &&
        AtomInfoKnownPolymerResName(LexStr(G, ai2->resn))))
    return false;

  // don't connect water atoms in different residues
  if (water_flag && !AtomInfoSameResidue(G, ai1, ai2)) return false;

  // don't connect atoms with different, non-NULL alternate conformations
  if (ai1->alt[0] != ai2->alt[0] && ai1->alt[0] && ai2->alt[0]) return false;

  // if either is a cation, unbond is user wants
  if (unbond_cations &&
      (AtomInfoIsFreeCation(G, ai1) || AtomInfoIsFreeCation(G, ai2)))
    return false;

  return true;
}

/**
 * Do bonding of atoms in `I`, using distances and/or temporary bonds in `cs`.
 *
 * Incorporates `cs->TmpBond` unless `connect_mode` is 2.
 * Incorporates `cs->TmpLinkBond`.
 *
 * @param I Molecule to modify
 * @param cs Coordinates and temporary bonds to consider
 * @param bondSearchMode If false and `connect_mode` != 2, do not search for new
 * bonds (only use TmpBond/TmpLinkBond).
 * @param connectModeOverride Overrides `connect_mode` setting if not -1
 * @param pbc Use periodic boundary conditions (find symop bonds)
 *
 * `connect_mode` options:
 * 0 = distance-based (excluding HETATM to HETATM) and CONECT records (default)
 * 1 = CONECT records
 * 2 = distance-based, ignores CONECT records
 * 3 = distance-based (including HETATM to HETATM) and CONECT records
 * 4 = same as `connect_mode` = 0 (special meaning during mmCIF loading)
 */
bool ObjectMoleculeConnect(ObjectMolecule* I, CoordSet* cs, bool bondSearchMode,
                           int connectModeOverride, bool pbc) {
  return ObjectMoleculeConnect(I, I->NBond, I->Bond, cs, bondSearchMode,
                               connectModeOverride, pbc);
}
bool ObjectMoleculeConnect(ObjectMolecule* I, int& nBond,
                           pymol::vla<BondType>& bondvla, struct CoordSet* cs,
                           int bondSearchMode, int connectModeOverride,
                           bool pbc) {
  PyMOLGlobals* G = I->G;
  AtomInfoType* const ai = I->AtomInfo.data();
  auto discrete_chains = SettingGet<int>(G, cSetting_pdb_discrete_chains);
  auto const connect_bonded = SettingGet<bool>(G, cSetting_connect_bonded);
  auto const unbond_cations = SettingGet<int>(G, cSetting_pdb_unbond_cations);
  auto const cutoff_v = SettingGet<float>(G, cSetting_connect_cutoff);
  auto const max_cutoff = cutoff_v + 0.2F;  ///< Sulfur cutoff
  auto connect_mode = (connectModeOverride >= 0)
                          ? connectModeOverride
                          : SettingGet<int>(G, cSetting_connect_mode);

  if (connect_mode == 2) {
    // Force use of distance-based connectivity, ignoring that
    // provided with file.
    bondSearchMode = true;
    cs->NTmpBond = 0;
    VLAFreeP(cs->TmpBond);
  } else if (connect_mode == 4) {
    // mmCIF specific, fall back to default to get any bonds for PDB, XYZ, etc.
    connect_mode = 0;
  }

  nBond = 0;
  auto const maxBond = cs->NIndex * 8;
  // Number of bonds is typically close to number of atoms
  bondvla.reserve(cs->NIndex * 1.2);
  p_return_val_if_fail(bondvla, false);  // memory error

  bool repeat = false;
  switch (connect_mode) {
    case 0: /* distance-based and explicit (not HETATM to HETATM) */
    case 2: /* distance-based only */
    case 3: /* distance-based and explicit (even HETATM to HETATM) */
      repeat = bondSearchMode && cs->NIndex > 0;
  }

  // Distance-based bond location
  while (repeat) {
    repeat = false;
    nBond = 0;

    // For monitoring excessive numbers of bonds
    int violations = 0;
    auto const max_violations = cs->NIndex >> 3;  // 12%
    auto const cnt = std::make_unique<signed char[]>(size_t(cs->NIndex));
    p_return_val_if_fail(cnt, false);  // memory error

    /* initialize each atom's (max) expected valence */
    for (unsigned i = 0; i < cs->NIndex; ++i) {
      auto valcnt = AtomInfoGetExpectedValence(G, ai + cs->IdxToAtm[i]);
      cnt[i] = (valcnt < 0) ? 6 : valcnt;
    }

    // Search for symop bonds (periodic boundary conditions)?
    int offset_begin = 0, offset_end = 1, symmat_end = 1;
    if (pbc && cs->getSymmetry() &&  //
        !cs->getSymmetry()->Crystal.isSuspicious()) {
      offset_begin = -1;
      offset_end = 2;
      symmat_end = cs->getSymmetry()->getNSymMat();
      assert(symmat_end > 0);
    }

    /* make a map of the local neighborhood in space */
    auto const map = std::unique_ptr<MapType>(
        new MapType(G, (max_cutoff + MAX_VDW) * (offset_begin ? -1 : 1),  //
                    cs->Coord, cs->NIndex));
    p_return_val_if_fail(map, false);  // memory error
    MapSetupExpress(map.get());  // Don't let MapEIter call this in omp parallel

    /// Return false on error
    auto const find_bonds_for_atom = [&](unsigned i, float const* v1,
                                         pymol::SymOp const& symop) -> bool {
      auto const a1 = cs->IdxToAtm[i];
      auto* const ai1 = ai + a1;

      for (const auto j : MapEIter(*map, v1)) {
        if (i <= j && !symop) continue;

        /* position in space for atom 2 */
        auto const* const v2 = cs->coordPtr(j);
        auto const a2 = cs->IdxToAtm[j];
        auto* const ai2 = ai + a2;

        if (!is_distance_bonded(G, cs, ai1, ai2, v1, v2, cutoff_v, connect_mode,
                                discrete_chains, connect_bonded,
                                unbond_cations))
          continue;

        /* we have a bond, now process it */

        int order = 1;
        if (!ai1->hetatm || ai1->resn == G->lex_const.MSE) {
          if (AtomInfoSameResidue(I->G, ai1, ai2)) {
            /* hookup standard disconnected PDB residue */
            assign_pdb_known_residue(G, ai1, ai2, &order);
          }
        }

#pragma omp critical
        {
          auto const bnd = bondvla.check(nBond++);
          BondTypeInit2(bnd, a2, a1,
                        -order /* store tentative valence as negative */);
          bnd->symop_2 = symop;

          /* if we allow bonds between chains and it screws up
           * the bonding, disallow inter-chain bonds */
          if (discrete_chains < 0) {
            /* decrement free valences, since we have a bond */
            if (--cnt[i] == -2) violations++;
            if (--cnt[j] == -2) violations++;

            if (violations > max_violations) {
              PRINTFB(G, FB_ObjectMolecule, FB_Blather)
              " %s: Assuming chains are discrete...\n", __func__ ENDFB(G);

              discrete_chains = 1;
              repeat = 1;
            }
          }
        }

        if (repeat) {
          return false;
        }
      }

      return true;
    };

    bool break_all = false;

    // Do bond search in parallel for every atom
#pragma omp parallel for
    for (int i = 0; i < cs->NIndex; ++i) {
      float _v1_buf[3];
      pymol::SymOp symop{};
      for (symop.x = offset_begin; symop.x < offset_end; ++symop.x) {
        for (symop.y = offset_begin; symop.y < offset_end; ++symop.y) {
          for (symop.z = offset_begin; symop.z < offset_end; ++symop.z) {
            for (symop.index = 0; symop.index != symmat_end; ++symop.index) {
              auto const* const v1 = cs->coordPtrSym(i, symop, _v1_buf);
              assert(v1);

              if (break_all || !find_bonds_for_atom(i, v1, symop) ||
                  nBond > maxBond) {
                goto break_all_find_bonds_for_atom;
              }
            }
          }
        }
      }

      continue;
    break_all_find_bonds_for_atom:
      // can't "break" inside an omp parallel block
      break_all = true;
    }

    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " %s: Found %d bonds.\n", __func__, nBond ENDFB(G);
  }

  /* if we have explicit connectivity, determine if we need to set
   * check_conect_all */
  if (cs->NTmpBond && cs->TmpBond) {
    bool check_conect_all = false;
    bool pdb_conect_all = false;

    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " %s: incorporating explicit bonds. %d %d\n", __func__, nBond,
        cs->NTmpBond ENDFB(G);
    if ((nBond == 0) && (cs->NTmpBond > 0) && bondSearchMode &&
        (connect_mode == 0) && cs->NIndex) {
      /* if no bonds were found, and we have explicit connectivity,
       * try to determine if we need to set pdb_conect_mode */
      for (unsigned i = 0; i < cs->NIndex; ++i) {
        auto const& ai1 = ai[cs->IdxToAtm[i]];
        if (ai1.bonded && !ai1.hetatm) {
          /* apparent PDB ATOM record with explicit bonding... */
          check_conect_all = true;
          break;
        }
      }
    }

    bondvla.check(nBond + cs->NTmpBond);
    p_return_val_if_fail(bondvla, false);  // memory error

    auto* ii1 = bondvla.data() + nBond;
    auto* ii2 = cs->TmpBond.data();
    int n_atom = I->NAtom;
    for (unsigned a = 0; a < cs->NTmpBond; ++a) {
      auto const a1 = cs->IdxToAtm[ii2->index[0]];
      auto const a2 = cs->IdxToAtm[ii2->index[1]];
      if ((a1 >= 0) && (a2 >= 0) && (a1 < n_atom) && (a2 < n_atom)) {
        if (check_conect_all) {
          if ((!ai[a1].hetatm) && (!ai[a2].hetatm)) {
            /* found bond between non-HETATMs -- so tell PyMOL to CONECT all
             * ATOMs when writing out a PDB file */
            pdb_conect_all = true;
          }
        }
        ai[a1].bonded = true;
        ai[a2].bonded = true;
        *ii1 = std::move(*ii2);
        ii1->index[0] = a1;
        ii1->index[1] = a2;
        ii1++;
        ii2++;
        nBond++;
      }
    }

    VLAFreeP(cs->TmpBond);
    cs->NTmpBond = 0;

    if (pdb_conect_all) {
      int dummy;
      if (!SettingGetIfDefined_b(G, I->Setting.get(), cSetting_pdb_conect_all,
                                 &dummy)) {
        {
          auto handle = I->getSettingHandle(-1);
          if (handle) {
            SettingCheckHandle(G, *handle);
            SettingSet_b(handle->get(), cSetting_pdb_conect_all, true);
          }
        }
      }
    }
  }

  /* Link b/t ligand and residue? */
  if (cs->NTmpLinkBond && cs->TmpLinkBond) {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    "%s: incorporating linkage bonds. %d %d\n", __func__, nBond,
        cs->NTmpLinkBond ENDFB(G);

    bondvla.check(nBond + cs->NTmpLinkBond);
    p_return_val_if_fail(bondvla, false);  // memory error

    auto* ii1 = bondvla.data() + nBond;
    auto* ii2 = cs->TmpLinkBond.data();
    for (unsigned a = 0; a < cs->NTmpLinkBond; ++a) {
      auto const a1 = ii2->index[0];               /* first atom is in object */
      auto const a2 = cs->IdxToAtm[ii2->index[1]]; /* second is in the cset */
      ai[a1].bonded = true;
      ai[a2].bonded = true;
      (*ii1) = std::move(*ii2);
      ii1->index[0] = a1;
      ii1->index[1] = a2;
      ii1++;
      ii2++;
    }
    nBond += cs->NTmpLinkBond;
    VLAFreeP(cs->TmpLinkBond);
    cs->NTmpLinkBond = 0;
  }

  // Eliminate duplicates
  // TODO do we expect any?
  // TODO should we also check with swapped indices?
  // TODO why not for discrete objects?
  if (nBond > 1 && !I->DiscreteFlag) {
    PRINTFD(G, FB_ObjectMolecule)
    " %s: elminating duplicates with %d bonds...\n", __func__, nBond ENDFD;

    UtilSortInPlace(G, bondvla.data(), nBond, sizeof(BondType),
                    (UtilOrderFn*)BondInOrder);
    auto* ii1 = bondvla.data();
    auto* ii2 = bondvla.data() + 1;
    for (int a = nBond; --a;) {
      if (BondCompare(ii2, ii1) != 0) {
        if (++ii1 != ii2) {
          *ii1 = std::move(*ii2);
        }
      } else if (ii2->order > 0 && ii1->order < 0) {
        // use most certain valence
        ii1->order = ii2->order;
      }
      ii2++;
    }
    nBond = ii1 - bondvla.data() + 1;
  }

  bondvla.resize(nBond);

  // restore bond order positivity
  for (auto bnd = bondvla.begin(), bnd_end = bnd + nBond; bnd != bnd_end;
       ++bnd) {
    if (bnd->order < 0) bnd->order = -bnd->order;
  }

  PRINTFD(G, FB_ObjectMolecule)
  " %s: leaving with %d bonds...\n", __func__, nBond ENDFD;
  return true;
}

void ObjectMoleculeConnectDiscrete(ObjectMolecule* I, int searchFlag,
                                   int connectModeOverride, bool pbc) {
  for (int i = 0; i < I->NCSet; i++) {
    if (!I->CSet[i]) {
      continue;
    }

    int nbond = 0;
    pymol::vla<BondType> bond;

    ObjectMoleculeConnect(I, nbond, bond, I->CSet[i], searchFlag,
                          connectModeOverride, pbc);

    if (!bond) {
      continue;
    }

    if (!I->Bond) {
      I->Bond = std::move(bond);
    } else {
      I->Bond.check(I->NBond + nbond - 1);
      std::copy_n(bond.data(), nbond, I->Bond.data() + I->NBond);
    }

    I->NBond += nbond;
  }
}
