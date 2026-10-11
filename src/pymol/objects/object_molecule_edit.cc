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

// Editing molecules: adding and removing atoms and bonds, fusing, hydrogens and
// undo.

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/executive/executive.h"
#include "pymol/geometry/sphere.h"
#include "pymol/executive/editor.h"
#include "pymol/chem/hydrogen_adder.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"
#include "pymol/objects/object_molecule_internal.h"
int ObjectMoleculeFixSeleHydrogens(ObjectMolecule* I, int sele, int state) {
  int a;
  int seleFlag = false;
  const AtomInfoType* ai0;
  int ok = true;

  ai0 = I->AtomInfo;
  for (a = 0; a < I->NAtom; a++) {
    if (SelectorIsMember(I->G, ai0->selEntry, sele)) {
      seleFlag = true;
      break;
    }
    ai0++;
  }
  if (seleFlag) {
    seleFlag = false;
    if (!ObjectMoleculeVerifyChemistry(I, state)) {
      ErrMessage(I->G, " AddHydrogens",
                 "missing chemical geometry information.");
    } else {
      ai0 = I->AtomInfo;
      for (a = 0; a < I->NAtom; a++) {
        if (!ai0->isHydrogen()) { /* only do heavies */
          if (SelectorIsMember(I->G, ai0->selEntry, sele)) {
            for (StateIterator iter(I->G, I->Setting.get(), state, I->NCSet);
                 iter.next();) {
              auto cs = I->CSet[iter.state];
              if (!cs) continue;

              seleFlag |=
                  ObjectMoleculeSetMissingNeighborCoords(I, cs, a, true);
            }
          }
        }
        ai0++;
      }
    }
    if (seleFlag) I->invalidate(cRepAll, cRepInvAll, -1);
  }
  return ok;
}
int ObjectMoleculeRenameAtoms(ObjectMolecule* I, int* flag, int force) {
  PyMOLGlobals* G = I->G;
  AtomInfoType* ai;
  int a;
  int result;
  if (force) {
    ai = I->AtomInfo.data();
    if (!flag) {
      for (a = 0; a < I->NAtom; a++) {
        LexAssign(G, ai->name, 0);
        ai++;
      }
    } else {
      for (a = 0; a < I->NAtom; a++) {
        if (flag[a]) LexAssign(G, ai->name, 0);
        ai++;
      }
    }
  }
  result = AtomInfoUniquefyNames(I->G, nullptr, 0, I->AtomInfo.data(), flag,
                                 I->NAtom);
  return result;
}
/**
 * Transform coordinates of `cs` in place and then append it to `tcs`.
 *
 * @param coord_orig Original cs coordinates
 * @param at0 Anchor atom index
 * @param mat1_inv Source coordinate system
 * @param valence_vec Valence vector. If nullptr, don't move. If null-vector,
 * don't merge.
 */
static bool AddCoordinateIntoCoordSet(CoordSet* tcs, CoordSet* cs,
                                      const float* coord_orig, int at0,
                                      const float* mat1_inv,
                                      const float* valence_vec) {
  if (!tcs) {
    return true;
  }

  if (valence_vec) {
    if (lengthsq3f(valence_vec) < 0.1) {
      // null-vector
      return true;
    }

    // anchor
    auto idx0 = tcs->atmToIdx(at0);
    assert(idx0 >= 0);
    const float* va0 = tcs->coordPtr(idx0);

    // Target coordinate system
    float matT[16], mat[16];
    identity44f(matT);
    copy3f(valence_vec, matT);
    get_system1f3f(matT, matT + 4, matT + 8);
    copy3(va0, matT + 12);
    transpose44f44f(matT, mat);

    // combine with source coordinate system
    right_multiply44f44f(mat, mat1_inv);

    // transform coordset
    for (int idx = 0; idx < cs->NIndex; idx++) {
      transform44f3f(mat, &coord_orig[3 * idx], cs->coordPtr(idx));
    }
  }

  return CoordSetMerge(tcs->Obj, tcs, cs);
}

/**
 * @param at0 Atom index
 * @param index0 If different from at0: Index of atom to replace (e.g. Hydrogen)
 * @param[out] out (3x1) Normalized open valence vector
 * @return True if valid open valence vector found
 */
static bool GetTargetValenceVector(const CoordSet* tcs, int at0, int index0,
                                   float* out) {
  if (!tcs) {
    return false;
  }

  if (at0 != index0) {
    // anchor
    auto ca0 = tcs->atmToIdx(at0);
    if (ca0 < 0) {
      return false;
    }

    // hydrogen (or other replaced atom)
    auto ch0 = tcs->atmToIdx(index0);
    if (ch0 < 0) {
      return false;
    }

    const float* vh0 = tcs->coordPtr(ch0);
    const float* va0 = tcs->coordPtr(ca0);
    subtract3f(vh0, va0, out);
    return true;
  }

  return CoordSetFindOpenValenceVector(tcs, at0, out);
}
/**
 * Merge src into I.
 *
 * Optionally creates a bond between index0 and index1, or the atoms bound to
 * index0 and index1 if they are "single geometry" atoms like hydrogens.
 *
 * @param index0 Atom index in I
 * @param index1 Atom index in src
 * @param create_bond If false, then don't create a bond, just combine the
 * objects. Implies move_flag=false.
 * @param move_flag If true, then transform the source coordinates to form a
 * reasonable bond geometry.
 */
pymol::Result<> ObjectMoleculeFuse(ObjectMolecule* I, int const index0,
                                   const ObjectMolecule* src, int const index1,
                                   bool create_bond, bool move_flag) {
  constexpr int state1 = 0;
  constexpr bool edit = true;

  auto* G = I->G;

  if (!create_bond) {
    move_flag = false;
  }

  auto const* scs = src->getCoordSet(state1);
  if (!scs) {
    return pymol::make_error("no source coordset");
  }

  auto* ai0 = I->AtomInfo.data();
  auto const* ai1 = src->AtomInfo.data();

  /**
   * If we want a bond and atm has cAtomInfoSingle geometry, then return the
   * atom index of its (only) neighbor. Otherwise, return atm.
   */
  auto const get_anchor_atm = [create_bond](const ObjectMolecule* obj,
                                            int atm) {
    if (create_bond) {
      if (obj->AtomInfo[atm].geom == cAtomInfoSingle) {
        auto neighbors = AtomNeighbors(obj, atm);
        if (neighbors.size() == 1) {
          atm = neighbors[0].atm;
        }
      }
    }
    return atm;
  };

  int const at0 = get_anchor_atm(I, index0);
  int const at1 = get_anchor_atm(src, index1);

  assert(!(at0 < 0 || at1 < 0));

  auto const anch1 = scs->atmToIdx(at1);
  if (anch1 < 0) {
    return pymol::make_error("no coordinate for source anchor atom");
  }

  /* copy atoms and atom info into a 1:1 direct mapping */
  auto cs = std::make_unique<CoordSet>(G);
  cs->setNIndex(scs->NIndex);
  cs->enumIndices();

  auto nai = pymol::vla<AtomInfoType>(scs->NIndex);

  for (int idx = 0; idx < scs->NIndex; ++idx) {
    copy3f(scs->coordPtr(idx), cs->coordPtr(idx));
    auto atm1 = scs->IdxToAtm[idx];

    AtomInfoCopy(G, ai1 + atm1, nai + idx);

    if (edit) {
      if (atm1 == at1) {
        nai[idx].temp1 = 2; /* clear marks */
      } else {
        nai[idx].temp1 = 0; /* clear marks */
      }
    }
  }

  /* copy internal bond information */
  cs->TmpBond.reserve(src->NBond);
  cs->NTmpBond = 0;
  for (int b = 0; b != src->NBond; ++b) {
    auto const& b1 = src->Bond[b];
    auto a0 = scs->atmToIdx(b1.index[0]);
    auto a1 = scs->atmToIdx(b1.index[1]);
    if (a0 >= 0 && a1 >= 0) {
      auto& b0 = cs->TmpBond[cs->NTmpBond++];
      b0 = b1;
      b0.index[0] = a0;
      b0.index[1] = a1;
    }
  }

  // source coordinate system
  float mat1_inv[16];

  if (create_bond) {
    const float* va1 = scs->coordPtr(anch1);

    if (at0 != index0) {
      ai0[index0].deleteFlag = true;
    }

    identity44f(mat1_inv);

    if (at1 != index1) {
      auto const hydr1 = scs->atmToIdx(index1);
      if (hydr1 == -1) {
        return pymol::make_error("no source attachment vector found");
      }
      nai[hydr1].deleteFlag = true;
      const float* vh1 = scs->coordPtr(hydr1);
      subtract3f(va1, vh1, mat1_inv);
    } else {
      auto found = CoordSetFindOpenValenceVector(scs, at1, mat1_inv);
      if (!found) {
        return pymol::make_error("no source attachment vector found");
      }
      scale3f(mat1_inv, -1.0F, mat1_inv);
    }

    // bond length
    float const d = AtomInfoGetBondLength(G, ai0 + at0, ai1 + at1);

    // rotation matrix
    get_system1f3f(mat1_inv, mat1_inv + 4, mat1_inv + 8);

    // translation vector
    float mat1_trans[16];
    identity44f(mat1_trans);
    mat1_trans[3] = (mat1_inv[0] * d - va1[0]);
    mat1_trans[7] = (mat1_inv[1] * d - va1[1]);
    mat1_trans[11] = (mat1_inv[2] * d - va1[2]);

    right_multiply44f44f(mat1_inv, mat1_trans);

    /* set up the linking bond */
    cs->TmpLinkBond.resize(1);
    cs->NTmpLinkBond = 1;
    auto* bond = cs->TmpLinkBond.data();
    BondTypeInit2(bond, at0, anch1);
  }

  AtomInfoUniquefyNames(I, nai.data(), cs->NIndex);

  /* set up tags which will enable use to continue editing bond */

  if (edit) {
    for (int atm = 0; atm < I->NAtom; ++atm) {
      ai0[atm].temp1 = 0;
    }
    ai0[at0].temp1 = 1;
  }

  int const state0 = I->DiscreteFlag ? (ai0[at0].discrete_state - 1)
                                     : cSelectorUpdateTableAllStates;
  assert(!I->DiscreteFlag || I->DiscreteCSet[at0] == I->CSet[state0]);

  // Get target valence vectors for all relevant coordinate sets. Do this before
  // merging atoms to not leave the object molecule in an half-done state in
  // case we can't find valence vectors.
  std::vector<std::array<float, 3>> target_vectors;
  if (move_flag) {
    bool found_any_target_vector = false;
    target_vectors.resize(I->NCSet);

    for (StateIterator iter(G, nullptr, state0, I->NCSet); iter.next();) {
      float* vec = target_vectors[iter.state].data();
      if (GetTargetValenceVector(I->CSet[iter.state], at0, index0, vec)) {
        found_any_target_vector = true;
      } else {
        zero3(vec);
      }
    }

    if (!found_any_target_vector) {
      return pymol::make_error("no target attachment vector found");
    }
  }

  // will free nai, cs->TmpBond and cs->TmpLinkBond
  // invalidates the ai0 pointer!
  bool ok = ObjectMoleculeMerge(I, std::move(nai), cs.get(), false,
                                cAIC_AllMask, true) &&
            ObjectMoleculeExtendIndices(I, state0);
  p_return_val_if_fail(ok, pymol::Error::MEMORY);

  // Get untransformed copy of source coordinates
  pymol::vla<float> coord_orig;
  if (move_flag) {
    coord_orig = cs->Coord;
    p_return_val_if_fail(coord_orig, pymol::Error::MEMORY);
  }

  // Add coordinates into the coordinate set(s)
  for (StateIterator iter(G, nullptr, state0, I->NCSet); iter.next();) {
    const float* const vec =
        move_flag ? target_vectors[iter.state].data() : nullptr;
    ok = AddCoordinateIntoCoordSet(I->CSet[iter.state], cs.get(),
                                   coord_orig.data(), at0, mat1_inv, vec);
    p_return_val_if_fail(ok, pymol::Error::MEMORY);
  }

  ok = ObjectMoleculeSort(I);
  p_return_val_if_fail(ok, pymol::Error::MEMORY);

  ObjectMoleculeUpdateIDNumbers(I);

  if (at0 != index0 || at1 != index1) {
    ObjectMoleculePurge(I);
  }

  // edit the resulting bond
  if (edit) {
    int atm_pk1 = -1;
    int atm_pk2 = -1;
    for (int atm = 0; atm < I->NAtom; ++atm) {
      if (I->AtomInfo[atm].temp1 == 1) {
        atm_pk2 = atm;
      } else if (I->AtomInfo[atm].temp1 == 2) {
        atm_pk1 = atm;
      }
    }
    if (atm_pk2 >= 0 && atm_pk1 >= 0) {
      auto sele1 = pymol::string_format("%s`%d", I->Name, atm_pk1 + 1);
      auto sele2 = pymol::string_format("%s`%d", I->Name, atm_pk2 + 1);
      EditorSelect(G, sele1.data(), sele2.data(), "", "", false, true, true);
    }
  }

  return {};
}
int ObjectMoleculeAttach(ObjectMolecule* I, int index,
                         pymol::vla<AtomInfoType>&& nai) {
  int a;
  AtomInfoType* ai;
  BondType* bond;
  float v[3], v0[3], d;
  CoordSet* cs = nullptr;
  int ok = false;

  ai = I->AtomInfo + index;

  ok_assert(1, cs = CoordSetNew(I->G));
  ok_assert(1, cs->Coord = pymol::vla<float>(3));

  cs->NIndex = 1;
  ok_assert(1, cs->TmpLinkBond = pymol::vla<BondType>(1));

  cs->NTmpLinkBond = 1;
  bond = cs->TmpLinkBond.data();
  BondTypeInit2(bond, index, 0);
  cs->enumIndices();

  ok_assert(1, ObjectMoleculePrepareAtom(I, index, nai.data()));
  d = AtomInfoGetBondLength(I->G, ai, nai);

  ok_assert(1, ObjectMoleculeMerge(I, std::move(nai), cs, false, cAIC_AllMask,
                                   true));  // will free nai and cs->TmpLinkBond
  ok_assert(1, ObjectMoleculeExtendIndices(I, -1));

  for (a = 0; a < I->NCSet; a++) { /* add atom to each coordinate set */
    if (const auto* cs_a = I->CSet[a]) {
      CoordSetGetAtomVertex(cs_a, index, v0);
      CoordSetFindOpenValenceVector(cs_a, index, v);
      scale3f(v, d, v);
      add3f(v0, v, cs->Coord.data());
      ok_assert(1, CoordSetMerge(I, I->CSet[a], cs));
    }
  }

  ok_assert(1, ObjectMoleculeSort(I));
  ObjectMoleculeUpdateIDNumbers(I);

  ok = true;
ok_except1:
  delete cs;
  return ok;
}
int ObjectMoleculeFillOpenValences(ObjectMolecule* I, int index) {
  int a;
  AtomInfoType* ai;
  int result = 0;
  int flag = true;
  float v[3], v0[3], d;
  CoordSet* cs = nullptr;
  int ok = true;

  if ((index >= 0) && (index <= I->NAtom)) {
    while (ok) {
      ai = I->AtomInfo + index;
      auto const nn = AtomNeighbors(I, index).size();

      assert(flag);

      if ((nn >= ai->valence) || (!flag)) break;
      flag = false;

      if (ok) cs = CoordSetNew(I->G);
      CHECKOK(ok, cs);
      if (ok) {
        cs->Coord = pymol::vla<float>(3);
        CHECKOK(ok, cs->Coord);

        cs->NIndex = 1;
        if (ok) cs->TmpLinkBond = pymol::vla<BondType>(1);
        CHECKOK(ok, cs->TmpLinkBond);
        if (ok) {
          cs->NTmpLinkBond = 1;
          auto* bond = cs->TmpLinkBond.data();
          BondTypeInit2(bond, index, 0);
        }
      }

      if (ok) cs->enumIndices();
      auto atInfo = pymol::vla<AtomInfoType>(1);
      AtomInfoType* nai = atInfo.data();
      if (ok) {
        UtilNCopy(nai->elem, "H", 2);
        nai->geom = cAtomInfoSingle;
        nai->valence = 1;
        ok &= ObjectMoleculePrepareAtom(I, index, nai);
        d = AtomInfoGetBondLength(I->G, ai, nai);
        if (ok)
          ok &=
              ObjectMoleculeMerge(I, std::move(atInfo), cs, false, cAIC_AllMask,
                                  true); /* will free nai and cs->TmpLinkBond */
      }
      if (ok) ok &= ObjectMoleculeExtendIndices(I, -1);
      for (a = 0; ok && a < I->NCSet;
           a++) { /* add atom to each coordinate set */
        if (auto* cs_a = I->CSet[a]) {
          CoordSetGetAtomVertex(cs_a, index, v0);
          CoordSetFindOpenValenceVector(cs_a, index, v);
          scale3f(v, d, v);
          add3f(v0, v, cs->Coord.data());
          ok &= CoordSetMerge(I, cs_a, cs);
        }
      }
      delete cs;
      result++;
      flag = true;
    }
  }
  ObjectMoleculeUpdateIDNumbers(I);
  return (result);
}

#define MaxOcc 100
/**
 * @param[out] normal (3x1) Normal vector of planar geometry at atom `atm`
 *
 * @return True if the atom has planar geometry and the normal could be computed
 *
 * @pre neighbors are defined
 */
static bool get_planer_normal(const CoordSet* cs, int const atm,
                              float* normal) {
  int found = false;
  int nOcc = 0;
  float occ[MaxOcc * 3];
  float v0[3], v1[3], v2[3], n0[3];

  if (CoordSetGetAtomVertex(cs, atm, v0)) {
    const auto* I = cs->Obj;

    // look for an attached non-hydrogen as a base
    for (auto const& neighbor : AtomNeighbors(I, atm)) {
      if (CoordSetGetAtomVertex(cs, neighbor.atm, v1)) {
        subtract3f(v1, v0, n0);
        normalize3f(n0); /* n0's point away from center atom */
        copy3f(n0, occ + 3 * nOcc);
        nOcc++;
        if (nOcc == MaxOcc) /* safety valve */
          break;
      }
    }
    switch (I->AtomInfo[atm].geom) {
      case cAtomInfoPlanar:
        if (nOcc > 1) {
          cross_product3f(occ, occ + 3, normal);
          if (nOcc > 2) {
            cross_product3f(occ, occ + 6, v2);
            if (dot_product3f(normal, v2) < 0) {
              subtract3f(normal, v2, normal);
            } else {
              add3f(normal, v2, normal);
            }
            cross_product3f(occ + 3, occ + 6, v2);
            if (dot_product3f(normal, v2) < 0) {
              subtract3f(normal, v2, normal);
            } else {
              add3f(normal, v2, normal);
            }
          }
          normalize3f(normal);
          found = true;
        }
        break;
    }
  }
  return found;
}
/**
 * @param atm Atom index
 * @param[out] v (3x1) Normalized open valence vector
 * @param seek (3x1) Primer for the direction, or nullptr for random priming.
 * @param ignore_index Atom index of neighbor to ignore
 * @return True if valid open valence vector found
 */
bool CoordSetFindOpenValenceVector(const CoordSet* cs, int atm, float* v,
                                   const float* seek, int ignore_index) {
  int nOcc = 0;
  float occ[MaxOcc * 3];
  int last_occ = -1;
  float v0[3], v1[3], n0[3] = {0.0F, 0.0F, 0.0F}, t[3];
  int result = false;
  float y[3], z[3];

  /* default is +X */
  v[0] = 1.0;
  v[1] = 0.0;
  v[2] = 0.0;

  if (cs) {
    const auto* I = cs->Obj;

    if (atm >= 0 && atm <= I->NAtom) {
      if (CoordSetGetAtomVertex(cs, atm, v0)) {
        const auto* ai = I->AtomInfo.data() + atm;

        // look for an attached non-hydrogen as a base
        for (auto const& neighbor : AtomNeighbors(I, atm)) {
          auto const a1 = neighbor.atm;
          if (a1 != ignore_index) {
            if (CoordSetGetAtomVertex(cs, a1, v1)) {
              last_occ = a1;
              subtract3f(v1, v0, n0);
              normalize3f(n0); /* n0's point away from center atom */
              copy3f(n0, occ + 3 * nOcc);
              nOcc++;
              if (nOcc == MaxOcc) /* safety valve */
                break;
            }
          }
        }
        if ((!nOcc) || (nOcc > 4) || (ai->geom == cAtomInfoNone)) {
          if (!seek)
            get_random3f(v);
          else
            copy3f(seek, v);
          result = true;
        } else {
          switch (nOcc) {
            case 1: /* only one current occupied position */
              switch (ai->geom) {
                case cAtomInfoTetrahedral:
                  if (!seek) {
                    get_system1f3f(occ, y, z);
                    scale3f(occ, -0.334F, v);
                    scale3f(z, 0.943F, t);
                    add3f(t, v, v);
                  } else { /* point hydrogen towards sought vector */
                    copy3f(seek, z);
                    get_system2f3f(occ, z, y);
                    scale3f(occ, -0.334F, v);
                    scale3f(z, 0.943F, t);
                    add3f(t, v, v);
                  }
                  result = true;
                  break;
                case cAtomInfoPlanar: {
                  if (!seek) {
                    if (last_occ >= 0 && get_planer_normal(cs, last_occ, n0)) {
                      copy3f(n0, y);
                      get_system2f3f(occ, y, z);
                    } else {
                      get_system1f3f(occ, y, z);
                    }
                    scale3f(occ, -0.500F, v);
                    scale3f(z, 0.866F, t);
                    add3f(t, v, v);
                  } else {
                    copy3f(seek, z);
                    get_system2f3f(occ, z, y);
                    scale3f(occ, -0.500F, v);
                    scale3f(z, 0.866F, t);
                    add3f(t, v, v);
                  }
                  result = true;
                } break;
                case cAtomInfoLinear:
                  scale3f(occ, -1.0F, v);
                  result = true;
                  break;
                default:
                  if (!seek)
                    get_random3f(v);
                  else
                    copy3f(seek, v);
                  result = false;
                  break;
              }
              break;
            case 2: /* only two current occupied positions */
              switch (ai->geom) {
                case cAtomInfoTetrahedral:
                  add3f(occ, occ + 3, t);
                  get_system2f3f(t, occ, z);
                  scale3f(t, -1.0F, v);
                  if (seek) {
                    if (dot_product3f(z, seek) < 0.0F) {
                      invert3f(z);
                    }
                  }
                  scale3f(z, 1.41F, t);
                  add3f(t, v, v);
                  result = true;
                  break;
                case cAtomInfoPlanar:
                  add3f(occ, occ + 3, t);
                  scale3f(t, -1.0F, v);
                  result = true;
                  break;
                default:
                  if (!seek) {
                    add3f(occ, occ + 3, t);
                    scale3f(t, -1.0F, v);
                    if (length3f(t) < 0.1) get_random3f(v);
                  } else
                    copy3f(seek, v);
                  /* hypervalent */
                  result = false;
                  break;
              }
              break;
            case 3: /* only three current occupied positions */
              switch (ai->geom) {
                case cAtomInfoTetrahedral:
                  add3f(occ, occ + 3, t);
                  add3f(occ + 6, t, t);
                  scale3f(t, -1.0F, v);
                  result = true;
                  break;
                default:
                  if (!seek) {
                    add3f(occ, occ + 3, t);
                    add3f(occ + 6, t, t);
                    scale3f(t, -1.0F, v);
                    if (length3f(t) < 0.1) get_random3f(v);
                  } else
                    copy3f(seek, v);
                  /* hypervalent */
                  result = false;
                  break;
              }
              break;
            case 4:
              if (!seek)
                get_random3f(v);
              else
                copy3f(seek, v);
              /* hypervalent */
              result = false;
              break;
          }
        }
      }
    }
  }
  normalize3f(v);
  return (result);
#undef MaxOcc
}
void ObjectMoleculeCreateSpheroid(ObjectMolecule* I, int average) {
  CoordSet* cs;
  int a, b, c, a0;
  SphereRec* sp;
  float *v, *v0, *s, *f, ang, min_dist, *max_sq;
  int* i;
  float* center = nullptr;
  float d0[3], n0[3], d1[3], d2[3];
  float p0[3], p1[3], p2[3];
  int t0, t1, t2, bt0, bt1, bt2;
  float dp, l, *fsum = nullptr;
  float spheroid_smooth;
  float spheroid_fill;
  float spheroid_ratio = 0.1F;    /* minimum ratio of width over length */
  float spheroid_minimum = 0.02F; /* minimum size - to insure valid normals */
  int row, *count = nullptr, base;
  int nRow;
  int first = 0;
  int last = 0;
  int current;
  int cscount;
  int n_state = 0;
  sp = GetSpheroidSphereRec(I->G);

  nRow = I->NAtom * sp->nDot;

  center = pymol::malloc<float>(I->NAtom * 3);
  count = pymol::malloc<int>(I->NAtom);
  fsum = pymol::malloc<float>(nRow);
  max_sq = pymol::malloc<float>(I->NAtom);

  spheroid_smooth = SettingGetGlobal_f(I->G, cSetting_spheroid_smooth);
  spheroid_fill = SettingGetGlobal_f(I->G, cSetting_spheroid_fill);
  /* first compute average coordinate */

  if (average < 1) average = I->NCSet;
  current = 0;
  cscount = 0;
  while (current < I->NCSet) {
    if (I->CSet[current]) {
      if (!cscount) first = current;
      cscount++;
      last = current + 1;
    }

    if (cscount == average || current == I->NCSet - 1) {
      PRINTFB(I->G, FB_ObjectMolecule, FB_Details)
      " ObjectMolecule: computing spheroid from states %d to %d.\n", first + 1,
          last ENDFB(I->G);

      auto spheroid = std::vector<float>(nRow);

      v = center;
      i = count;
      for (a = 0; a < I->NAtom; a++) {
        *(v++) = 0.0;
        *(v++) = 0.0;
        *(v++) = 0.0;
        *(i++) = 0;
      }

      for (b = first; b < last; b++) {
        cs = I->CSet[b];
        if (cs) {
          v = cs->Coord.data();
          for (a = 0; a < cs->NIndex; a++) {
            a0 = cs->IdxToAtm[a];
            v0 = center + 3 * a0;
            add3f(v, v0, v0);
            (*(count + a0))++;
            v += 3;
          }
        }
      }

      i = count;
      v = center;
      for (a = 0; a < I->NAtom; a++)
        if (*i) {
          (*(v++)) /= (*i);
          (*(v++)) /= (*i);
          (*(v++)) /= (*i++);
        } else {
          v += 3;
          i++;
        }

      /* now go through and compute radial distances */

      f = fsum;
      s = spheroid.data();
      for (a = 0; a < nRow; a++) {
        *(f++) = 0.0;
        *(s++) = 0.0;
      }

      v = max_sq;
      for (a = 0; a < I->NAtom; a++) *(v++) = 0.0;

      for (b = first; b < last; b++) {
        cs = I->CSet[b];
        if (cs) {
          v = cs->Coord.data();
          for (a = 0; a < cs->NIndex; a++) {
            a0 = cs->IdxToAtm[a];
            base = (a0 * sp->nDot);
            v0 = center + (3 * a0);
            subtract3f(v, v0, d0); /* subtract from average */
            l = lengthsq3f(d0);
            if (l > max_sq[a0]) max_sq[a0] = l;
            if (l > 0.0) {
              float isq = (float)(1.0 / sqrt1d(l));
              scale3f(d0, isq, n0);
              for (c = 0; c < sp->nDot; c++) { /* average over spokes */
                dp = dot_product3f(sp->dot[c], n0);
                row = base + c;
                if (dp >= 0.0) {
                  ang = (float)((acos(dp) / spheroid_smooth) * (cPI / 2.0));
                  if (ang > spheroid_fill) ang = spheroid_fill;
                  /* take envelop to zero over that angle */
                  if (ang <= (cPI / 2.0)) {
                    dp = (float)cos(ang);
                    fsum[row] += dp * dp;
                    spheroid[row] += l * dp * dp * dp;
                  }
                }
              }
            }
            v += 3;
          }
        }
      }

      f = fsum;
      s = spheroid.data();
      for (a = 0; a < I->NAtom; a++) {
        min_dist = (float)(spheroid_ratio * sqrt(max_sq[a]));
        if (min_dist < spheroid_minimum) min_dist = spheroid_minimum;
        for (b = 0; b < sp->nDot; b++) {
          if (*f > R_SMALL4) {
            (*s) =
                (float)(sqrt1d((*s) / (*(f++)))); /* we put the "rm" in "rms" */
          } else {
            f++;
          }
          if (*s < min_dist) *s = min_dist;
          s++;
        }
      }

      /* set frame 0 coordinates to the average */

      cs = I->CSet[first];
      if (cs) {
        v = cs->Coord.data();
        for (a = 0; a < cs->NIndex; a++) {
          a0 = cs->IdxToAtm[a];
          v0 = center + 3 * a0;
          copy3f(v0, v);
          v += 3;
        }
      }

      /* now compute surface normals */

      auto norm = std::vector<float>(nRow * 3);
      for (a = 0; a < nRow; a++) {
        zero3f(norm.data() + a * 3);
      }
      for (a = 0; a < I->NAtom; a++) {
        base = a * sp->nDot;
        for (b = 0; b < sp->NTri; b++) {
          t0 = sp->Tri[b * 3];
          t1 = sp->Tri[b * 3 + 1];
          t2 = sp->Tri[b * 3 + 2];
          bt0 = base + t0;
          bt1 = base + t1;
          bt2 = base + t2;
          copy3f(sp->dot[t0], p0);
          copy3f(sp->dot[t1], p1);
          copy3f(sp->dot[t2], p2);
          /*      scale3f(sp->dot[t0].v,spheroid[bt0],p0);
             scale3f(sp->dot[t1].v,spheroid[bt1],p1);
             scale3f(sp->dot[t2].v,spheroid[bt2],p2); */
          subtract3f(p1, p0, d1);
          subtract3f(p2, p0, d2);
          cross_product3f(d1, d2, n0);
          normalize3f(n0);
          v = norm.data() + bt0 * 3;
          add3f(n0, v, v);
          v = norm.data() + bt1 * 3;
          add3f(n0, v, v);
          v = norm.data() + bt2 * 3;
          add3f(n0, v, v);
        }
      }

      f = norm.data();
      for (a = 0; a < I->NAtom; a++) {
        base = a * sp->nDot;
        for (b = 0; b < sp->nDot; b++) {
          normalize3f(f);
          f += 3;
        }
      }

      if (I->CSet[first]) {
        I->CSet[first]->Spheroid = std::move(spheroid);
        I->CSet[first]->SpheroidNormal = std::move(norm);
      }

      for (b = first + 1; b < last; b++) {
        delete I->CSet[b];
        I->CSet[b] = nullptr;
      }

      if (n_state != first) {
        I->CSet[n_state] = I->CSet[first];
        I->CSet[first] = nullptr;
      }
      n_state++;

      cscount = 0;
    }
    current++;
  }
  I->NCSet = n_state;
  FreeP(center);
  FreeP(count);
  FreeP(fsum);
  FreeP(max_sq);

  I->invalidate(cRepSphere, cRepInvProp, -1);
}
void ObjectMoleculeReplaceAtom(ObjectMolecule* I, int index,
                               AtomInfoType&& ai) {
  if ((index >= 0) && (index <= I->NAtom)) {
    I->AtomInfo[index] = std::move(ai);
    I->invalidate(cRepAll, cRepInvAtoms, -1);
    /* could we put in a refinement step here? */
  }
}
int ObjectMoleculePrepareAtom(ObjectMolecule* I, int index, AtomInfoType* ai,
                              bool uniquefy) {
  /* match existing properties of the old atom */
  AtomInfoType* ai0;
  int ok = true;

  if ((index >= 0) && (index <= I->NAtom)) {
    ai0 = I->AtomInfo + index;
    ai->resv = ai0->resv;
    ai->hetatm = ai0->hetatm;
    ai->flags = ai0->flags;

    if (!ai->geom) ai->geom = ai0->geom;

    ai->discrete_state = ai0->discrete_state;
    ai->q = ai0->q;
    ai->b = ai0->b;
    strcpy(ai->alt, ai0->alt);
    ai->inscode = ai0->inscode;
    LexAssign(I->G, ai->segi, ai0->segi);
    LexAssign(I->G, ai->chain, ai0->chain);
    LexAssign(I->G, ai->resn, ai0->resn);
    ai->visRep = ai0->visRep;
    ai->id = -1;
    ai->rank = -1;

    AtomInfoAssignParameters(I->G, ai);

    if (uniquefy) {
      AtomInfoUniquefyNames(I->G, I->AtomInfo, I->NAtom, ai, nullptr, 1);
    }

    if ((ai->elem[0] == ai0->elem[0]) && (ai->elem[1] == ai0->elem[1]))
      ai->color = ai0->color;
    else if ((ai->elem[0] == 'C') && (ai->elem[1] == 0)) {
      int found = false;
      for (auto const& neighbor : AtomNeighbors(I, index)) {
        AtomInfoType const* ai1 = I->AtomInfo.data() + neighbor.atm;
        if (ai1->protons == cAN_C) {
          ai->color = ai1->color;
          found = true;
          break;
        }
      }
      if (ok && !found) {
        /* if no carbon nearby, then color according to the object color */
        ai->color = I->Color;
      }
    } else {
      AtomInfoAssignColors(I->G, ai);
    }
  }
  return ok;
}
int ObjectMoleculePreposReplAtom(ObjectMolecule* I, int index,
                                 AtomInfoType* ai) {
  float v0[3], v1[3], v[3];
  float d0[3], d, n0[3];
  int cnt;
  float t[3], sum[3];
  int a;
  int ncycle;
  int ok = true;
  if (ok) {
    for (a = 0; a < I->NCSet; a++) {
      if (I->CSet[a]) {
        if (ObjectMoleculeGetAtomVertex(I, a, index, v0)) {
          copy3f(v0, v); /* default is direct superposition */
          ncycle = -1;
          while (ncycle) {
            cnt = 0;
            zero3f(sum);
            /* look for an attached non-hydrogen as a base */
            for (auto const& neighbor : AtomNeighbors(I, index)) {
              auto const a1 = neighbor.atm;
              auto const* ai1 = I->AtomInfo.data() + a1;
              if (ai1->protons != 1)
                if (ObjectMoleculeGetAtomVertex(I, a, a1, v1)) {
                  d = AtomInfoGetBondLength(I->G, ai, ai1);
                  subtract3f(v0, v1, n0);
                  normalize3f(n0);
                  scale3f(n0, d, d0);
                  add3f(d0, v1, t);
                  add3f(t, sum, sum);
                  cnt++;
                }
            }
            if (cnt) {
              scale3f(sum, 1.0F / cnt, sum);
              copy3f(sum, v0);
              if ((cnt > 1) && (ncycle < 0)) ncycle = 5;
            }
            ncycle = abs(ncycle) - 1;
          }
          if (cnt) copy3f(sum, v);
          ObjectMoleculeSetAtomVertex(I, a, index, v);
        }
      }
    }
  }
  return ok;
}
void ObjectMoleculeSaveUndo(ObjectMolecule* I, int state, int log) {
  CoordSet* cs;
  PyMOLGlobals* G = I->G;
  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter] = -1;
  if (state < 0) state = 0;
  if (I->NCSet == 1) state = 0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if (cs) {
    I->UndoCoord[I->UndoIter] = pymol::malloc<float>(cs->NIndex * 3);
    memcpy(I->UndoCoord[I->UndoIter], cs->Coord,
           sizeof(float) * cs->NIndex * 3);
    I->UndoState[I->UndoIter] = state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }
  I->UndoIter = cUndoMask & (I->UndoIter + 1);
  ExecutiveSetLastObjectEdited(G, I);
  if (log) {
    OrthoLineType line;
    if (SettingGetGlobal_i(I->G, cSetting_logging)) {
      sprintf(line, "cmd.push_undo(\"%s\",%d)\n", I->Name, state + 1);
      PLog(G, line, cPLog_no_flush);
    }
  }
}
void ObjectMoleculeUndo(ObjectMolecule* I, int dir) {
  CoordSet* cs;
  int state;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter] = -1;
  state = SceneGetState(I->G);
  if (state < 0) state = 0;
  if (I->NCSet == 1) state = 0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if (cs) {
    I->UndoCoord[I->UndoIter] = pymol::malloc<float>(cs->NIndex * 3);
    memcpy(I->UndoCoord[I->UndoIter], cs->Coord,
           sizeof(float) * cs->NIndex * 3);
    I->UndoState[I->UndoIter] = state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }

  I->UndoIter = cUndoMask & (I->UndoIter + dir);
  if (!I->UndoCoord[I->UndoIter]) I->UndoIter = cUndoMask & (I->UndoIter - dir);

  if (I->UndoState[I->UndoIter] >= 0) {
    state = I->UndoState[I->UndoIter];
    if (state < 0) state = 0;

    if (I->NCSet == 1) state = 0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if (cs) {
      if (cs->NIndex == I->UndoNIndex[I->UndoIter]) {
        memcpy(cs->Coord.data(), I->UndoCoord[I->UndoIter],
               sizeof(float) * cs->NIndex * 3);
        I->UndoState[I->UndoIter] = -1;
        FreeP(I->UndoCoord[I->UndoIter]);
        cs->invalidateRep(cRepAll, cRepInvCoord);
        SceneChanged(I->G);
      }
    }
  }
}

int ObjectMoleculeAddBond(ObjectMolecule* I, int sele0, int sele1, int order,
                          pymol::zstring_view symop) {
  int a1, a2;
  AtomInfoType *ai1, *ai2;
  int s1, s2;
  int c = 0;

  /* TO DO: optimize for performance -- we shouldn't be doing full
     table scans */

  ai1 = I->AtomInfo.data();
  for (a1 = 0; a1 < I->NAtom; a1++) {
    s1 = ai1->selEntry;
    if (SelectorIsMember(I->G, s1, sele0)) {
      ai2 = I->AtomInfo.data();
      for (a2 = 0; a2 < I->NAtom; a2++) {
        s2 = ai2->selEntry;
        if (SelectorIsMember(I->G, s2, sele1)) {
          if (!I->Bond) {
            I->Bond = pymol::vla<BondType>(1);
          }
          if (I->Bond) {
            BondType* bnd = I->Bond.check(I->NBond);
            BondTypeInit2(bnd, a1, a2, order);

            assert(!bnd->symop_2);
            if (!symop.empty()) {
              bnd->symop_2.reset(symop.c_str());
            }

            I->NBond++;
            c++;
            I->AtomInfo[a1].chemFlag = false;
            I->AtomInfo[a2].chemFlag = false;
            I->AtomInfo[a1].bonded = true;
            I->AtomInfo[a2].bonded = true;
          }
        }
        ai2++;
      }
    }
    ai1++;
  }
  if (c) {
    I->invalidate(cRepAll, cRepInvBondsNoNonbonded, -1);
  }
  return (c);
}
pymol::Result<> ObjectMoleculeAddBondByIndices(ObjectMolecule* I, unsigned atm1,
                                               unsigned atm2, int order) {
  if (atm1 >= I->NAtom || atm2 >= I->NAtom) {
    return pymol::make_error("atom index out of bounds");
  }

  if (!I->Bond) {
    I->Bond = pymol::vla<BondType>(1);
  } else {
    I->Bond.check(I->NBond);
  }

  if (!I->Bond) {
    return pymol::Error::MEMORY;
  }

  auto& bnd = I->Bond[I->NBond++];
  BondTypeInit2(&bnd, atm1, atm2, order);

  I->AtomInfo[atm1].chemFlag = false;
  I->AtomInfo[atm2].chemFlag = false;
  I->AtomInfo[atm1].bonded = true;
  I->AtomInfo[atm2].bonded = true;

  I->invalidate(cRepAll, cRepInvBondsNoNonbonded, -1);

  return {};
}
int ObjectMoleculeAdjustBonds(ObjectMolecule* I, int sele0, int sele1, int mode,
                              int order, pymol::zstring_view symop) {
  auto const G = I->G;
  int a0, a1;
  int cnt = 0;
  BondType* b0;
  int both;
  int a;

  if (I->Bond) {
    b0 = I->Bond.data();
    for (a = 0; a < I->NBond; a++) {
      a0 = b0->index[0];
      a1 = b0->index[1];

      both = 0;

      if (SelectorIsMember(G, I->AtomInfo[a0].selEntry, sele0) &&
          SelectorIsMember(G, I->AtomInfo[a1].selEntry, sele1)) {
        both = 2;
      } else if (  //
          SelectorIsMember(G, I->AtomInfo[a1].selEntry, sele0) &&
          SelectorIsMember(G, I->AtomInfo[a0].selEntry, sele1)) {
        std::swap(a0, a1);
        both = 2;
      }

      if (both == 2) {
        cnt++;
        switch (mode) {
          case 0: /* cycle */
            switch (SettingGet_i(I->G, I->Setting.get(), nullptr,
                                 cSetting_editor_bond_cycle_mode)) {
              case 1: /* 1 arom 2 3 */
                switch (b0->order) {
                  case 1:
                    b0->order = 4;
                    break;
                  case 4:
                    b0->order = 2;
                    break;
                  case 2:
                    b0->order = 3;
                    break;
                  default:
                    b0->order = 1;
                    break;
                }
                break;
              case 2: /* 1 2 3 arom */
                b0->order++;
                if (b0->order > 4) b0->order = 1;
                break;
              default: /* old way -> 1 2 3 */
                b0->order++;
                if (b0->order > 3) b0->order = 1;
                break;
            }
            I->AtomInfo[a0].chemFlag = false;
            I->AtomInfo[a1].chemFlag = false;
            break;
          case 1: /* set */
            b0->order = order;
            I->AtomInfo[a0].chemFlag = false;
            I->AtomInfo[a1].chemFlag = false;
            break;
        }

        if (!symop.empty()) {
          b0->symop_2.reset(symop.c_str());
        }
      }
      b0++;
    }
    if (cnt) {
      I->invalidate(cRepLine, cRepInvBonds, -1);
      I->invalidate(cRepCyl, cRepInvBonds, -1);
      I->invalidate(cRepNonbonded, cRepInvBonds, -1);
      I->invalidate(cRepNonbondedSphere, cRepInvBonds, -1);
      I->invalidate(cRepRibbon, cRepInvBonds, -1);
      I->invalidate(cRepCartoon, cRepInvBonds, -1);
    }
  }

  return (cnt);
}
int ObjectMoleculeRemoveBonds(ObjectMolecule* I, int sele0, int sele1) {
  int a0, a1;
  int offset = 0;
  BondType *b0, *b1;
  int both;
  int s;
  int a;

  if (I->Bond) {
    offset = 0;
    b0 = I->Bond.data();
    b1 = I->Bond.data();
    for (a = 0; a < I->NBond; a++) {
      a0 = b0->index[0];
      a1 = b0->index[1];

      both = 0;
      s = I->AtomInfo[a0].selEntry;
      if (SelectorIsMember(I->G, s, sele0)) both++;
      s = I->AtomInfo[a1].selEntry;
      if (SelectorIsMember(I->G, s, sele1)) both++;
      if (both < 2) { /* reverse combo */
        both = 0;
        s = I->AtomInfo[a1].selEntry;
        if (SelectorIsMember(I->G, s, sele0)) both++;
        s = I->AtomInfo[a0].selEntry;
        if (SelectorIsMember(I->G, s, sele1)) both++;
      }

      if (both == 2) {
        AtomInfoPurgeBond(I->G, b0);
        offset--;
        b0++;
        I->AtomInfo[a0].chemFlag = false;
        I->AtomInfo[a1].chemFlag = false;
      } else if (offset) {
        *(b1++) = *(b0++); /* copy bond info */
      } else {
        *(b1++) = *(b0++); /* copy bond info */
      }
    }
    if (offset) {
      I->NBond += offset;
      VLASize(I->Bond, BondType, I->NBond);
      I->invalidate(cRepLine, cRepInvBonds, -1);
      I->invalidate(cRepCyl, cRepInvBonds, -1);
      I->invalidate(cRepNonbonded, cRepInvBonds, -1);
      I->invalidate(cRepNonbondedSphere, cRepInvBonds, -1);
      I->invalidate(cRepRibbon, cRepInvBonds, -1);
      I->invalidate(cRepCartoon, cRepInvBonds, -1);
    }
  }

  return (-offset);
}
void ObjectMoleculePurge(ObjectMolecule* I) {
  PyMOLGlobals* G = I->G;
  int offset = 0;

  // remove the object selection and free up any selection entries.
  // note that we don't delete atom selection members -- those may be needed in
  // the new object.
  SelectorDelete(G, I->Name);

  // old-to-new mapping
  auto oldToNew = std::vector<int>(size_t(I->NAtom), -1);
  for (int atm = 0; atm < I->NAtom; ++atm) {
    auto& ai0 = I->AtomInfo[atm];
    if (ai0.deleteFlag) {
      AtomInfoPurge(G, &ai0);
      offset--;
      assert(oldToNew[atm] == -1);
    } else {
      if (offset) {
        I->AtomInfo[atm + offset] = std::move(ai0);
      }
      oldToNew[atm] = atm + offset;
    }
  }
  if (offset) {
    I->NAtom += offset;
    I->AtomInfo.resize(I->NAtom);
    for (int a = 0; a < I->NCSet; ++a) {
      if (auto cs = I->CSet[a]) CoordSetAdjustAtmIdx(cs, oldToNew.data());
    }
    if (I->CSTmpl) {
      CoordSetAdjustAtmIdx(I->CSTmpl, oldToNew.data());
    }
  }

  I->updateAtmToIdx();

  // step 4, bonds

  offset = 0;
  auto b0 = I->Bond.data();
  auto b1 = I->Bond.data();
  for (int b = 0; b < I->NBond; ++b) {
    auto a0 = b0->index[0];
    auto a1 = b0->index[1];
    if (a0 < 0 || a1 < 0 || (oldToNew[a0] < 0) || (oldToNew[a1] < 0)) {
      /* deleting bond */
      AtomInfoPurgeBond(I->G, b0);
      offset--;
      b0++;
    } else {
      if (offset) {
        *b1 = *b0;
      }
      b1->index[0] = oldToNew[a0]; /* copy bond info */
      b1->index[1] = oldToNew[a1];
      b0++;
      b1++;
    }
  }
  if (offset) {
    I->NBond += offset;
    VLASize(I->Bond, BondType, I->NBond);
  }

  I->invalidate(cRepAll, cRepInvAtoms, -1);
}
