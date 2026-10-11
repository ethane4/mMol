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

// Distances, angles, dihedrals and VDW overlaps between selections.

#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/select/selector_internal.h"

int SelectorVdwFit(PyMOLGlobals* G, int sele1, int state1, int sele2,
                   int state2, float buffer, int quiet) {
  CSelector* I = G->Selector;
  float sumVDW = 0.0, dist;
  int a1, a2;
  AtomInfoType *ai1, *ai2;
  int at1, at2;
  CoordSet *cs1, *cs2;
  ObjectMolecule *obj1, *obj2;
  int idx1, idx2;
  int a;

  if (state1 < 0) state1 = 0;
  if (state2 < 0) state2 = 0;

  if (state1 != state2) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  auto vla = SelectorGetInterstateVector(G, sele1, state1, sele2, state2,
                                         2 * MAX_VDW + buffer);
  const int c = vla.size() / 2;

  if (c) {
    auto adj = std::vector<float>(vla.size());

    for (a = 0; a < c; a++) {
      a1 = vla[a * 2];
      a2 = vla[a * 2 + 1];

      at1 = I->Table[a1].atom;
      at2 = I->Table[a2].atom;

      obj1 = I->Obj[I->Table[a1].model];
      obj2 = I->Obj[I->Table[a2].model];

      if ((state1 < obj1->NCSet) && (state2 < obj2->NCSet)) {
        cs1 = obj1->CSet[state1];
        cs2 = obj2->CSet[state2];
        if (cs1 && cs2) { /* should always be true */

          ai1 = obj1->AtomInfo + at1;
          ai2 = obj2->AtomInfo + at2;

          idx1 = cs1->atmToIdx(at1);
          idx2 = cs2->atmToIdx(at2);

          sumVDW = ai1->vdw + ai2->vdw;
          dist = (float)diff3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2));

          if (dist < (sumVDW + buffer)) {
            float shift = (dist - (sumVDW + buffer)) / 2.0F;
            adj[2 * a] = ai1->vdw + shift;
            adj[2 * a + 1] = ai2->vdw + shift;
          } else {
            adj[2 * a] = ai1->vdw;
            adj[2 * a + 1] = ai2->vdw;
          }
        }
      }
    }

    for (a = 0; a < c; a++) {
      a1 = vla[a * 2];
      a2 = vla[a * 2 + 1];

      at1 = I->Table[a1].atom;
      at2 = I->Table[a2].atom;

      obj1 = I->Obj[I->Table[a1].model];
      obj2 = I->Obj[I->Table[a2].model];

      if ((state1 < obj1->NCSet) && (state2 < obj2->NCSet)) {
        cs1 = obj1->CSet[state1];
        cs2 = obj2->CSet[state2];
        if (cs1 && cs2) { /* should always be true */

          ai1 = obj1->AtomInfo + at1;
          ai2 = obj2->AtomInfo + at2;

          if (adj[2 * a] < ai1->vdw) {
            ai1->vdw = adj[2 * a];
          }

          if (adj[2 * a + 1] < ai2->vdw) {
            ai2->vdw = adj[2 * a + 1];
          }
        }
      }
    }
  }

  return true;
}
float SelectorSumVDWOverlap(PyMOLGlobals* G, int sele1, int state1, int sele2,
                            int state2, float adjust) {
  CSelector* I = G->Selector;
  float result = 0.0;
  float sumVDW = 0.0, dist;
  int a1, a2;
  AtomInfoType *ai1, *ai2;
  int at1, at2;
  CoordSet *cs1, *cs2;
  ObjectMolecule *obj1, *obj2;
  int idx1, idx2;
  int a;

  if (state1 < 0) state1 = 0;
  if (state2 < 0) state2 = 0;

  if (state1 != state2) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  auto vla = SelectorGetInterstateVector(G, sele1, state1, sele2, state2,
                                         2 * MAX_VDW + adjust);
  const int c = vla.size() / 2;

  for (a = 0; a < c; a++) {
    a1 = vla[a * 2];
    a2 = vla[a * 2 + 1];

    at1 = I->Table[a1].atom;
    at2 = I->Table[a2].atom;

    obj1 = I->Obj[I->Table[a1].model];
    obj2 = I->Obj[I->Table[a2].model];

    if ((state1 < obj1->NCSet) && (state2 < obj2->NCSet)) {
      cs1 = obj1->CSet[state1];
      cs2 = obj2->CSet[state2];
      if (cs1 && cs2) { /* should always be true */

        ai1 = obj1->AtomInfo + at1;
        ai2 = obj2->AtomInfo + at2;

        idx1 = cs1->atmToIdx(at1);
        idx2 = cs2->atmToIdx(at2);

        sumVDW = ai1->vdw + ai2->vdw + adjust;
        dist = (float)diff3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2));

        if (dist < sumVDW) {
          result += ((sumVDW - dist) / 2.0F);
        }
      }
    }
  }
  return (result);
}
DistSet* SelectorGetDistSet(PyMOLGlobals* G, DistSet* ds, int sele1, int state1,
                            int sele2, int state2, int mode, float cutoff,
                            float* result) {
  CSelector* I = G->Selector;
  std::vector<int> vla;
  int c;
  float dist;
  int a1, a2;
  AtomInfoType *ai1, *ai2;
  int at, at1, at2;
  CoordSet *cs1, *cs2;
  ObjectMolecule *obj, *obj1, *obj2, *lastObj;
  int idx1, idx2;
  int a;
  int nv = 0;
  float *vv0, *vv1;
  float dist_sum = 0.0;
  int dist_cnt = 0;
  int s;
  int a_keeper = false;
  int *zero = nullptr, *scratch = nullptr;
  std::vector<bool> coverage;
  HBondCriteria hbcRec, *hbc;
  int exclusion = 0;
  int bonds_only = 0;
  int from_proton = SettingGetGlobal_b(G, cSetting_h_bond_from_proton);
  AtomInfoType* h_ai;

  bool cutoff_is_ratio_distance_to_vdW = false;

  /* if we're creating hydrogen bonds, then set some distance cutoffs */
  switch (mode) {
    case 1:
      bonds_only = 1;
      break;
    case 2:
      exclusion = SettingGetGlobal_i(G, cSetting_h_bond_exclusion);
      break;
    case 8:
      cutoff_is_ratio_distance_to_vdW = true;
      mode = 3;
      // no break, continue with case 3
    case 3:
      exclusion = SettingGetGlobal_i(G, cSetting_distance_exclusion);
      break;
  }

  hbc = &hbcRec;
  *result = 0.0;
  /* if the dist set exists, get info from it, otherwise get a new one */
  if (!ds) {
    ds = DistSetNew(G);
  } else {
    nv = ds->NIndex; /* number of vertices */
  }

  auto& vv = ds->Coord;
  vv.reserve(10);

  /* update states: if the two are the same, update that one state, else update
   * all states */
  if ((state1 < 0) || (state2 < 0) || (state1 != state2)) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  /* find and prepare (neighbortables) in any participating Molecular objects */
  if ((mode == 1) || (mode == 2) ||
      (mode == 3)) { /* fill in all the neighbor tables */
    int max_n_atom = I->Table.size();
    lastObj = nullptr;
    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      /* foreach atom in the session, get its identifier and ObjectMolecule to
       * which it belongs */
      at = I->Table[a].atom; /* grab the atom ID from the Selectors->Table */
      obj = I->Obj[I->Table[a].model]; /* -- JV -- quick way to get an object
                                          from an atom */
      s = obj->AtomInfo[at]
              .selEntry; /* grab the selection entry# from this Atoms Info */
      if (obj != lastObj) {
        if (max_n_atom < obj->NAtom) max_n_atom = obj->NAtom;
        /* if the current atom is in sele1 or sele2 then update it's object's
         * neighbor table */
        if (SelectorIsMember(G, s, sele1) || SelectorIsMember(G, s, sele2)) {
          /* if hbonds (so, more than just distance) */
          if (mode == 2) ObjectMoleculeVerifyChemistry(obj, -1);
          lastObj = obj;
        }
      }
    }
    /* prepare these for the next round */
    zero = pymol::calloc<int>(max_n_atom);
    scratch = pymol::malloc<int>(max_n_atom);
  }

  /* if we're hydrogen bonding, setup the cutoff */
  if (mode == 2) {
    ObjectMoleculeInitHBondCriteria(G, hbc);
    if (cutoff < 0.0F) {
      cutoff = hbc->maxDistAtMaxAngle;
      if (cutoff < hbc->maxDistAtZero) {
        cutoff = hbc->maxDistAtZero;
      }
    }
  }
  if (cutoff < 0) cutoff = 1000.0;

  if (mode == 4) {
    // centroid distance
    float centroid1[3], centroid2[3];
    ObjectMoleculeOpRec op;

    // get centroid 1
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_CSetSumVertices;
    op.cs1 = state1;
    ExecutiveObjMolSeleOp(G, sele1, &op);

    if (op.i1 > 0) {
      scale3f(op.v1, 1.f / op.i1, centroid1);

      // get centroid 2
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CSetSumVertices;
      op.cs1 = state2;
      ExecutiveObjMolSeleOp(G, sele2, &op);

      if (op.i1 > 0) {
        scale3f(op.v1, 1.f / op.i1, centroid2);

        // store positions in measurement object
        VLACheck(vv, float, (nv * 3) + 6);
        vv0 = vv + (nv * 3);
        nv += 2;
        copy3f(centroid1, vv0);
        copy3f(centroid2, vv0 + 3);

        // for the return value
        dist_cnt = 1;
        dist_sum = diff3f(centroid1, centroid2);
      }
    }

    // skip searching for pairwise distances
    c = 0;
  } else {
    /* coverage determines if a given atom appears in sel1 and sel2 */
    coverage.resize(I->Table.size());

    for (SelectorAtomIterator iter(I); iter.next();) {
      s = iter.getAtomInfo()->selEntry;
      if (SelectorIsMember(G, s, sele1) && SelectorIsMember(G, s, sele2))
        coverage[iter.a] = true;
    }

    float cutoff_map = cutoff;
    if (cutoff_is_ratio_distance_to_vdW) {
      constexpr float vdw_upper_bounds = 3.f;
      cutoff_map *= 2 * vdw_upper_bounds;
    }

    /* this creates an interleaved list of ints for mapping ids to states within
     * a given neighborhood */
    vla = SelectorGetInterstateVector(G, sele1, state1, sele2, state2,
                                      cutoff_map);
    c = vla.size() / 2;
  }

  /* for each state */
  for (a = 0; a < c; a++) {
    /* get the interstate atom identifier for the two atoms to distance */
    a1 = vla[a * 2];
    a2 = vla[a * 2 + 1];

    /* check their coverage to avoid duplicates */
    if (a1 < a2 || (a1 != a2 && !(coverage[a1] && coverage[a2])) ||
        (state1 != state2)) { /* eliminate reverse duplicates */
      /* get the object-local atom ID */
      at1 = I->Table[a1].atom;
      at2 = I->Table[a2].atom;
      /* get the object for this global atom ID */
      obj1 = I->Obj[I->Table[a1].model];
      obj2 = I->Obj[I->Table[a2].model];

      /* the states are valid for these two atoms */
      if ((state1 < obj1->NCSet) && (state2 < obj2->NCSet)) {
        /* get the coordinate sets for both atoms */
        cs1 = obj1->CSet[state1];
        cs2 = obj2->CSet[state2];
        if (cs1 && cs2) {
          /* for bonding */
          float* don_vv = nullptr;
          float* acc_vv = nullptr;

          /* grab the appropriate atom information for this object-local atom */
          ai1 = obj1->AtomInfo + at1;
          ai2 = obj2->AtomInfo + at2;

          idx1 = cs1->atmToIdx(at1);
          idx2 = cs2->atmToIdx(at2);

          if ((idx1 >= 0) && (idx2 >= 0)) {
            /* actual distance calculation from ptA to ptB */
            dist = (float)diff3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2));

            if (cutoff_is_ratio_distance_to_vdW) {
              dist /= ai1->vdw + ai2->vdw;
            }

            /* if we pass the boding cutoff */
            if (dist < cutoff) {
              float h_crd[3];
              h_ai = nullptr;

              a_keeper = true;
              if (exclusion && (obj1 == obj2)) {
                a_keeper = !SelectorCheckNeighbors(G, exclusion, obj1, at1, at2,
                                                   zero, scratch);
              } else if (bonds_only) {
                a_keeper =
                    SelectorCheckNeighbors(G, 1, obj1, at1, at2, zero, scratch);
              }
              if (a_keeper && (mode == 2)) {
                // We need to check:
                // situation 1:                      ai1 == donor and ai2 ==
                // acceptor situation 2 (if not situation 1): ai1 == acceptor
                // and ai2 == donor

                a_keeper = false;

                if (ai1->hb_donor && ai2->hb_acceptor) {
                  /* proton comes from ai1 */
                  a_keeper = ObjectMoleculeGetCheckHBond(
                      &h_ai, h_crd, obj1, at1, state1, obj2, at2, state2, hbc);
                  if (a_keeper) {
                    if (h_ai && from_proton) {
                      don_vv = h_crd;
                      ai1 = h_ai;
                    } else {
                      don_vv = cs1->coordPtr(idx1);
                    }
                    acc_vv = cs2->coordPtr(idx2);
                  }
                }

                // Check situation 2 only if no H-bond is found beforehand
                if ((!a_keeper) && (ai1->hb_acceptor && ai2->hb_donor)) {
                  /* proton comes from ai2 */
                  a_keeper = ObjectMoleculeGetCheckHBond(
                      &h_ai, h_crd, obj2, at2, state2, obj1, at1, state1, hbc);

                  if (a_keeper) {
                    if (h_ai && from_proton) {
                      don_vv = h_crd;
                      ai2 = h_ai;
                    } else {
                      don_vv = cs2->coordPtr(idx2);
                    }
                    acc_vv = cs1->coordPtr(idx1);
                  }
                }
              }

              if ((sele1 == sele2) && (at1 > at2)) a_keeper = false;

              if (a_keeper) {
                /* Insert DistInfo records for updating distances */
                /* Init/Add the elem to the DistInfo list */
                ds->MeasureInfo.emplace_front();
                auto* atom1Info = &ds->MeasureInfo.front();

                // TH
                atom1Info->id[0] = AtomInfoCheckUniqueID(G, ai1);
                atom1Info->id[1] = AtomInfoCheckUniqueID(G, ai2);

                atom1Info->offset = nv; /* offset into this DSet's Coord */
                atom1Info->state[0] = state1; /* state1 of sel1 */
                atom1Info->state[1] = state2;
                atom1Info->measureType = cRepDash; /* DISTANCE-dash */

                /* we have a distance we want to keep */
                dist_cnt++;
                dist_sum += dist;
                /* see if vv has room at another 6 floats */
                VLACheck(vv, float, (nv * 3) + 6);
                vv0 = vv + (nv * 3);

                if ((mode == 2) && (don_vv) && (acc_vv)) {
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(don_vv++);
                  *(vv0++) = *(acc_vv++);
                  *(vv0++) = *(acc_vv++);
                  *(vv0++) = *(acc_vv++);
                } else {
                  vv1 = cs1->coordPtr(idx1);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  vv1 = cs2->coordPtr(idx2);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                  *(vv0++) = *(vv1++);
                }

                nv += 2;
              }
            }
          }
        }
      }
    }
  }
  if (dist_cnt) (*result) = dist_sum / dist_cnt;
  FreeP(zero);
  FreeP(scratch);
  if (vv) VLASize(vv, float, (nv + 1) * 3);
  ds->NIndex = nv;
  ds->Coord = vv;
  return (ds);
}

DistSet* SelectorGetAngleSet(PyMOLGlobals* G, DistSet* ds, int sele1,
                             int state1, int sele2, int state2, int sele3,
                             int state3, int mode, float* angle_sum,
                             int* angle_cnt) {
  CSelector* I = G->Selector;
  int nv = 0;
  std::vector<bool> coverage;

  if (!ds) {
    ds = DistSetNew(G);
  } else {
    nv = ds->NAngleIndex;
  }

  auto& vv = ds->AngleCoord;
  vv.reserve(10);

  if ((state1 < 0) || (state2 < 0) || (state3 < 0) || (state1 != state2) ||
      (state1 != state3)) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  /* which atoms are involved? */

  {
    int a, s, at;
    ObjectMolecule* obj;

    coverage.resize(I->Table.size());
    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      at = I->Table[a].atom;
      obj = I->Obj[I->Table[a].model];
      s = obj->AtomInfo[at].selEntry;
      if (SelectorIsMember(G, s, sele1) && SelectorIsMember(G, s, sele3))
        coverage[a] = true;
    }
  }

  {
    int a, s, at;
    ObjectMolecule* obj;
    int* list1 = VLAlloc(int, 1000);
    int* list2 = VLAlloc(int, 1000);
    int* list3 = VLAlloc(int, 1000);
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int bonded12, bonded23;

    /* now generate three lists of atoms, one for each selection set */

    if (list1 && list2 && list3) {
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        at = I->Table[a].atom;
        obj = I->Obj[I->Table[a].model];
        s = obj->AtomInfo[at].selEntry;
        if (SelectorIsMember(G, s, sele1)) {
          VLACheck(list1, int, n1);
          list1[n1++] = a;
        }
        if (SelectorIsMember(G, s, sele2)) {
          VLACheck(list2, int, n2);
          list2[n2++] = a;
        }
        if (SelectorIsMember(G, s, sele3)) {
          VLACheck(list3, int, n3);
          list3[n3++] = a;
        }
      }

      /* for each set of 3 atoms in each selection... */

      {
        int i1, i2, i3;
        int a1, a2, a3;
        int at1, at2, at3;

        /*        AtomInfoType *ai1,*ai2,ai3; */
        CoordSet *cs1, *cs2, *cs3;
        ObjectMolecule *obj1, *obj2, *obj3;

        int idx1, idx2, idx3;
        float angle;
        float d1[3], d2[3];
        float *v1, *v2, *v3, *vv0;

        for (i1 = 0; i1 < n1; i1++) {
          a1 = list1[i1];
          at1 = I->Table[a1].atom;
          obj1 = I->Obj[I->Table[a1].model];

          if (state1 < obj1->NCSet) {
            cs1 = obj1->CSet[state1];

            if (cs1) {
              idx1 = cs1->atmToIdx(at1);

              if (idx1 >= 0) {
                for (i2 = 0; i2 < n2; i2++) {
                  a2 = list2[i2];
                  at2 = I->Table[a2].atom;
                  obj2 = I->Obj[I->Table[a2].model];

                  if (state2 < obj2->NCSet) {
                    cs2 = obj2->CSet[state2];

                    if (cs2) {
                      idx2 = cs2->atmToIdx(at2);

                      if (idx2 >= 0) {
                        /* neighbor table like BPRec */
                        bonded12 =
                            ObjectMoleculeAreAtomsBonded2(obj1, at1, obj2, at2);

                        for (i3 = 0; i3 < n3; i3++) {
                          a3 = list3[i3];

                          if ((a1 != a2 || state1 != state2) &&
                              (a2 != a3 || state2 != state3) &&
                              (a1 != a3 || state1 != state3)) {
                            if (!(coverage[a1] && coverage[a3]) || (a1 < a3) ||
                                (state1 !=
                                 state3)) { /* eliminate alternate-order
                                               duplicates */

                              at3 = I->Table[a3].atom;
                              obj3 = I->Obj[I->Table[a3].model];

                              if (state3 < obj3->NCSet) {
                                cs3 = obj3->CSet[state3];

                                if (cs3) {
                                  idx3 = cs3->atmToIdx(at3);

                                  if (idx3 >= 0) {
                                    bonded23 = ObjectMoleculeAreAtomsBonded2(
                                        obj2, at2, obj3, at3);

                                    if (!mode || ((mode == 1) &&
                                                  (bonded12 && bonded23))) {
                                      /* store the 3 coordinates */

                                      v1 = cs1->coordPtr(idx1);
                                      v2 = cs2->coordPtr(idx2);
                                      v3 = cs3->coordPtr(idx3);

                                      subtract3f(v1, v2, d1);
                                      subtract3f(v3, v2, d2);

                                      /* Insert DistInfo records for updating
                                       * distances */
                                      /* Init/Add the elem to the DistInfo list
                                       */
                                      ds->MeasureInfo.emplace_front();
                                      auto* atom1Info =
                                          &ds->MeasureInfo.front();

                                      // TH
                                      atom1Info->id[0] = AtomInfoCheckUniqueID(
                                          G, obj1->AtomInfo + at1);
                                      atom1Info->id[1] = AtomInfoCheckUniqueID(
                                          G, obj2->AtomInfo + at2);
                                      atom1Info->id[2] = AtomInfoCheckUniqueID(
                                          G, obj3->AtomInfo + at3);

                                      atom1Info->offset =
                                          nv; /* offset into this DSet's Coord
                                               */
                                      atom1Info->state[0] =
                                          state1; /* state1 of sel1 */
                                      atom1Info->state[1] = state2;
                                      atom1Info->state[2] = state3;
                                      atom1Info->measureType = cRepAngle;

                                      angle = get_angle3f(d1, d2);

                                      (*angle_sum) += angle;
                                      (*angle_cnt)++;

                                      VLACheck(vv, float, (nv * 3) + 14);
                                      vv0 = vv + (nv * 3);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v1++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v2++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = *(v3++);
                                      *(vv0++) = (float)!bonded12;
                                      /* show line 1 flag */
                                      *(vv0++) = (float)!bonded23;
                                      *(vv0++) = 0.0F;
                                      *(vv0++) =
                                          0.0F; /* label x relative to v2 */
                                      *(vv0++) =
                                          0.0F; /* label y relative to v2 */
                                      *(vv0++) =
                                          0.0F; /* label z relative to v2 */
                                      nv += 5;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    VLAFreeP(list1);
    VLAFreeP(list2);
    VLAFreeP(list3);
  }

  if (vv) VLASize(vv, float, (nv + 1) * 3);
  ds->NAngleIndex = nv;
  ds->AngleCoord = vv;
  return (ds);
}

DistSet* SelectorGetDihedralSet(PyMOLGlobals* G, DistSet* ds, int sele1,
                                int state1, int sele2, int state2, int sele3,
                                int state3, int sele4, int state4, int mode,
                                float* angle_sum, int* angle_cnt) {
  CSelector* I = G->Selector;
  int nv = 0;
  std::vector<bool> coverage14;
  std::vector<bool> coverage23;
  ObjectMolecule* just_one_object = nullptr;
  int just_one_atom[4] = {-1, -1, -1, -1};

  if (!ds) {
    ds = DistSetNew(G);
  } else {
    nv = ds->NDihedralIndex;
  }

  auto& vv = ds->DihedralCoord;
  vv.reserve(10);

  if ((state1 < 0) || (state2 < 0) || (state3 < 0) || (state4 < 0) ||
      (state1 != state2) || (state1 != state3) || (state1 != state4)) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }

  /* which atoms are involved? */

  {
    int a, s, at;
    ObjectMolecule* obj;

    coverage14.resize(I->Table.size());
    coverage23.resize(I->Table.size());
    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      bool coverage1 = false;
      bool coverage2 = false;
      at = I->Table[a].atom;
      obj = I->Obj[I->Table[a].model];
      if (!a) just_one_object = obj;
      s = obj->AtomInfo[at].selEntry;
      if (SelectorIsMember(G, s, sele1)) {
        if (obj != just_one_object)
          just_one_object = nullptr;
        else if (just_one_atom[0] == -1)
          just_one_atom[0] = a;
        else
          just_one_atom[0] = -2;
        coverage1 = true;
      }
      if (SelectorIsMember(G, s, sele2)) {
        if (obj != just_one_object)
          just_one_object = nullptr;
        else if (just_one_atom[1] == -1)
          just_one_atom[1] = a;
        else
          just_one_atom[1] = -2;
        coverage2 = true;
      }
      if (SelectorIsMember(G, s, sele3)) {
        if (obj != just_one_object)
          just_one_object = nullptr;
        else if (just_one_atom[2] == -1)
          just_one_atom[2] = a;
        else
          just_one_atom[2] = -2;
        coverage23[a] = coverage2;
      }
      if (SelectorIsMember(G, s, sele4)) {
        if (obj != just_one_object)
          just_one_object = nullptr;
        else if (just_one_atom[3] == -1)
          just_one_atom[3] = a;
        else
          just_one_atom[3] = -2;
        coverage14[a] = coverage1;
      }
    }
  }

  {
    int a, s, at;
    ObjectMolecule* obj;
    int* list1 = VLAlloc(int, 1000);
    int* list2 = VLAlloc(int, 1000);
    int* list3 = VLAlloc(int, 1000);
    int* list4 = VLAlloc(int, 1000);
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    int bonded12, bonded23, bonded34;

    /* now generate three lists of atoms, one for each selection set */

    if (list1 && list2 && list3 && list4) {
      if (just_one_object && (just_one_atom[0] >= 0) &&
          (just_one_atom[1] >= 0) && (just_one_atom[2] >= 0) &&
          (just_one_atom[3] >= 0)) {
        /* optimal case */

        list1[0] = just_one_atom[0];
        list2[0] = just_one_atom[1];
        list3[0] = just_one_atom[2];
        list4[0] = just_one_atom[3];

        n1 = n2 = n3 = n4 = 1;

      } else {
        for (a = cNDummyAtoms; a < I->Table.size(); a++) {
          at = I->Table[a].atom;
          obj = I->Obj[I->Table[a].model];
          s = obj->AtomInfo[at].selEntry;
          if (SelectorIsMember(G, s, sele1)) {
            VLACheck(list1, int, n1);
            list1[n1++] = a;
          }
          if (SelectorIsMember(G, s, sele2)) {
            VLACheck(list2, int, n2);
            list2[n2++] = a;
          }
          if (SelectorIsMember(G, s, sele3)) {
            VLACheck(list3, int, n3);
            list3[n3++] = a;
          }
          if (SelectorIsMember(G, s, sele4)) {
            VLACheck(list4, int, n4);
            list4[n4++] = a;
          }
        }
      }

      /* for each set of 3 atoms in each selection... */

      {
        int i1, i2, i3, i4;
        int a1, a2, a3, a4;
        int at1, at2, at3, at4;

        /*        AtomInfoType *ai1,*ai2,ai3; */
        CoordSet *cs1, *cs2, *cs3, *cs4;
        ObjectMolecule *obj1, *obj2, *obj3, *obj4;

        int idx1, idx2, idx3, idx4;
        float angle;
        float *v1, *v2, *v3, *v4, *vv0;

        for (i1 = 0; i1 < n1; i1++) {
          a1 = list1[i1];
          at1 = I->Table[a1].atom;
          obj1 = I->Obj[I->Table[a1].model];
          if (state1 < obj1->NCSet) {
            cs1 = obj1->CSet[state1];

            if (cs1) {
              idx1 = cs1->atmToIdx(at1);

              if (idx1 >= 0) {
                for (i2 = 0; i2 < n2; i2++) {
                  a2 = list2[i2];
                  at2 = I->Table[a2].atom;
                  obj2 = I->Obj[I->Table[a2].model];

                  if (state2 < obj2->NCSet) {
                    cs2 = obj2->CSet[state2];

                    if (cs2) {
                      idx2 = cs2->atmToIdx(at2);

                      if (idx2 >= 0) {
                        bonded12 =
                            ObjectMoleculeAreAtomsBonded2(obj1, at1, obj2, at2);

                        if (!mode || ((mode == 1) && bonded12))
                          for (i3 = 0; i3 < n3; i3++) {
                            a3 = list3[i3];
                            at3 = I->Table[a3].atom;
                            obj3 = I->Obj[I->Table[a3].model];

                            if (state3 < obj3->NCSet) {
                              cs3 = obj3->CSet[state3];

                              if (cs3) {
                                idx3 = cs3->atmToIdx(at3);

                                if (idx3 >= 0) {
                                  bonded23 = ObjectMoleculeAreAtomsBonded2(
                                      obj2, at2, obj3, at3);
                                  if (!mode || ((mode == 1) && bonded23))
                                    for (i4 = 0; i4 < n4; i4++) {
                                      a4 = list4[i4];

                                      if ((a1 != a2) && (a1 != a3) &&
                                          (a1 != a4) && (a2 != a3) &&
                                          (a2 != a4) && (a3 != a4)) {
                                        if (!(coverage14[a1] &&
                                              coverage14[a4] &&
                                              coverage23[a2] &&
                                              coverage23[a3]) ||
                                            (a1 < a4)) {
                                          /* eliminate alternate-order
                                           * duplicates */

                                          at4 = I->Table[a4].atom;
                                          obj4 = I->Obj[I->Table[a4].model];

                                          if (state4 < obj4->NCSet) {
                                            cs4 = obj4->CSet[state4];

                                            if (cs4) {
                                              idx4 = cs3->atmToIdx(at4);

                                              if (idx4 >= 0) {
                                                bonded34 =
                                                    ObjectMoleculeAreAtomsBonded2(
                                                        obj3, at3, obj4, at4);

                                                if (!mode ||
                                                    ((mode == 1) && bonded34)) {
                                                  /* store the 3 coordinates */

                                                  v1 = cs1->coordPtr(idx1);
                                                  v2 = cs2->coordPtr(idx2);
                                                  v3 = cs3->coordPtr(idx3);
                                                  v4 = cs4->coordPtr(idx4);

                                                  /* Insert DistInfo records for
                                                   * updating distances */
                                                  /* Init/Add the elem to the
                                                   * DistInfo list */
                                                  ds->MeasureInfo
                                                      .emplace_front();
                                                  auto* atom1Info =
                                                      &ds->MeasureInfo.front();

                                                  // TH
                                                  atom1Info->id[0] =
                                                      AtomInfoCheckUniqueID(
                                                          G,
                                                          obj1->AtomInfo + at1);
                                                  atom1Info->id[1] =
                                                      AtomInfoCheckUniqueID(
                                                          G,
                                                          obj2->AtomInfo + at2);
                                                  atom1Info->id[2] =
                                                      AtomInfoCheckUniqueID(
                                                          G,
                                                          obj3->AtomInfo + at3);
                                                  atom1Info->id[3] =
                                                      AtomInfoCheckUniqueID(
                                                          G,
                                                          obj4->AtomInfo + at4);

                                                  atom1Info->offset =
                                                      nv; /* offset into this
                                                             DSet's Coord */

                                                  atom1Info->state[0] =
                                                      state1; /* state1 of sel1
                                                               */
                                                  atom1Info->state[1] = state2;
                                                  atom1Info->state[2] = state3;
                                                  atom1Info->state[3] = state4;

                                                  atom1Info->measureType =
                                                      cRepDihedral;

                                                  angle = get_dihedral3f(
                                                      v1, v2, v3, v4);

                                                  (*angle_sum) += angle;
                                                  (*angle_cnt)++;

                                                  VLACheck(vv, float,
                                                           (nv * 3) + 17);
                                                  vv0 = vv + (nv * 3);
                                                  ObjectMoleculeGetAtomTxfVertex(
                                                      obj1, state1, at1, vv0);
                                                  ObjectMoleculeGetAtomTxfVertex(
                                                      obj2, state2, at2,
                                                      vv0 + 3);
                                                  ObjectMoleculeGetAtomTxfVertex(
                                                      obj3, state3, at3,
                                                      vv0 + 6);
                                                  ObjectMoleculeGetAtomTxfVertex(
                                                      obj4, state4, at4,
                                                      vv0 + 9);
                                                  vv0 += 12;
                                                  *(vv0++) = (float)!bonded12;
                                                  *(vv0++) = (float)!bonded23;
                                                  *(vv0++) = (float)!bonded34;
                                                  *(vv0++) =
                                                      0.0F; /* label x relative
                                                               to v2+v3/2 */
                                                  *(vv0++) =
                                                      0.0F; /* label y relative
                                                               to v2+v3/2 */
                                                  *(vv0++) =
                                                      0.0F; /* label z relative
                                                               to v2+v3/2 */
                                                  nv += 6;
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                }
                              }
                            }
                          }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    VLAFreeP(list1);
    VLAFreeP(list2);
    VLAFreeP(list3);
    VLAFreeP(list4);
  }

  if (vv) VLASize(vv, float, (nv + 1) * 3);
  ds->NDihedralIndex = nv;
  ds->DihedralCoord = vv;
  return (ds);
}
