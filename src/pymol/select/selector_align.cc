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

// Pairing the atoms of two selections for alignment.

#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/select/selector_internal.h"

#define cINTER_ENTRIES 11

int SelectorResidueVLAsTo3DMatchScores(PyMOLGlobals* G, CMatch* match,
                                       int* vla1, int n1, int state1, int* vla2,
                                       int n2, int state2, float seq_wt,
                                       float radius, float scale, float base,
                                       float coord_wt, float rms_exp) {
  CSelector* I = G->Selector;
  int a, b, *vla;
  int n_max = (n1 > n2) ? n1 : n2;
  float* inter1 = pymol::calloc<float>(cINTER_ENTRIES * n1);
  float* inter2 = pymol::calloc<float>(cINTER_ENTRIES * n2);
  float* v_ca = pymol::calloc<float>(3 * n_max);
  if (inter1 && inter2 && v_ca) {
    int pass;

    for (pass = 0; pass < 2; pass++) {
      ObjectMolecule* obj;
      const CoordSet* cs;
      const int* neighbor = nullptr;
      const AtomInfoType* atomInfo = nullptr;
      const ObjectMolecule* last_obj = nullptr;
      float** dist_mat;
      float* inter;
      int state;
      int n;
      if (!pass) {
        vla = vla1;
        state = state1;
        inter = inter1;
        n = n1;
        dist_mat = match->da;
      } else {
        vla = vla2;
        state = state2;
        inter = inter2;
        n = n2;
        dist_mat = match->db;
      }

      if (state < 0) state = 0;
      for (a = 0; a < n; a++) {
        int at_ca1;
        float* vv_ca = v_ca + a * 3;

        obj = I->Obj[vla[0]];
        at_ca1 = vla[1];
        if (obj != last_obj) {
          last_obj = obj;
          neighbor = obj->getNeighborArray();
          atomInfo = obj->AtomInfo;
        }

        if (state < obj->NCSet)
          cs = obj->CSet[state];
        else
          cs = nullptr;
        if (cs && neighbor && atomInfo) {
          int idx_ca1 = cs->atmToIdx(at_ca1);

          if (idx_ca1 >= 0) {
            int mem0, mem1, mem2, mem3, mem4;
            int nbr0, nbr1, nbr2, nbr3;
            const float* v_ca1 = cs->coordPtr(idx_ca1);
            int idx_cb1 = -1;
            int cnt = 0;

            copy3f(v_ca1, vv_ca);
            copy3f(v_ca1, inter + 8);

            /* find attached CB */

            mem0 = at_ca1;
            nbr0 = neighbor[mem0] + 1;
            while ((mem1 = neighbor[nbr0]) >= 0) {
              if ((atomInfo[mem1].protons == cAN_C) &&
                  (atomInfo[mem1].name == G->lex_const.CB)) {
                idx_cb1 = cs->atmToIdx(mem1);
                break;
              }
              nbr0 += 2;
            }

            /* find remote CA, CB */

            if (idx_cb1 >= 0) {
              const float* v_cb1 = cs->coordPtr(idx_cb1);

              mem0 = at_ca1;
              nbr0 = neighbor[mem0] + 1;
              while ((mem1 = neighbor[nbr0]) >= 0) {
                nbr1 = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr1]) >= 0) {
                  if (mem2 != mem0) {
                    int idx_ca2 = -1;

                    nbr2 = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr2]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0)) {
                        if ((atomInfo[mem3].protons == cAN_C) &&
                            (atomInfo[mem3].name == G->lex_const.CA)) {
                          idx_ca2 = cs->atmToIdx(mem3);
                          break;
                        }
                      }
                      nbr2 += 2;
                    }
                    if (idx_ca2 >= 0) {
                      const float* v_ca2 = cs->coordPtr(idx_ca2);

                      nbr2 = neighbor[mem2] + 1;
                      while ((mem3 = neighbor[nbr2]) >= 0) {
                        if ((mem3 != mem1) && (mem3 != mem0)) {
                          int idx_cb2 = -1;
                          nbr3 = neighbor[mem3] + 1;
                          while ((mem4 = neighbor[nbr3]) >= 0) {
                            if ((mem4 != mem2) && (mem4 != mem1) &&
                                (mem4 != mem0)) {
                              if ((atomInfo[mem4].protons == cAN_C) &&
                                  (atomInfo[mem4].name == G->lex_const.CB)) {
                                idx_cb2 = cs->atmToIdx(mem4);
                                break;
                              }
                            }
                            nbr3 += 2;
                          }

                          if (idx_cb2 >= 0) {
                            const float* v_cb2 = nullptr;
                            v_cb2 = cs->coordPtr(idx_cb2);
                            {
                              float angle =
                                  get_dihedral3f(v_cb1, v_ca1, v_ca2, v_cb2);
                              if (idx_cb1 < idx_cb2) {
                                inter[0] = (float)cos(angle);
                                inter[1] = (float)sin(angle);
                              } else {
                                inter[2] = (float)cos(angle);
                                inter[3] = (float)sin(angle);
                              }
                            }
                            cnt++;
                          }
                        }
                        nbr2 += 2;
                      }
                    }
                  }
                  nbr1 += 2;
                }
                nbr0 += 2;
              }
            }
          }
        }
        vla += 3;
        inter += cINTER_ENTRIES;
      }
      if (dist_mat) {
        for (a = 0; a < n; a++) { /* optimize this later */
          float* vv_ca = v_ca + a * 3;
          for (b = 0; b < n; b++) {
            float* vv_cb = v_ca + b * 3;
            float diff = (float)diff3f(vv_ca, vv_cb);
            dist_mat[a][b] = diff;
            dist_mat[b][a] = diff;
          }
        }
      }
      {
        std::unique_ptr<MapType> map(new MapType(G, radius, v_ca, n, nullptr));
        if (!pass) {
          inter = inter1;
        } else {
          inter = inter2;
        }
        if (map) {
          for (a = 0; a < n; a++) {
            float* v_ca1 = v_ca + 3 * a;
            float* i_ca1 = inter + cINTER_ENTRIES * a;
            for (const auto b : MapEIter(*map, v_ca1)) {
              float* v_ca2 = v_ca + 3 * b;
              if (a != b) {
                if (within3f(v_ca1, v_ca2, radius)) {
                  float* i_ca2 = inter + cINTER_ENTRIES * b;
                  i_ca1[4] += i_ca2[0]; /* add dihedral vectors head-to-tail */
                  i_ca1[5] += i_ca2[1];
                  i_ca1[6] += i_ca2[2];
                  i_ca1[7] += i_ca2[3];
                }
              }
            }
          }
          for (a = 0; a < n; a++) {
            float nf = (float)sqrt(inter[4] * inter[4] + inter[5] * inter[5]);
            if (nf > 0.0001F) {
              inter[4] = inter[4] / nf;
              inter[5] = inter[5] / nf;
            }
            nf = (float)sqrt(inter[6] * inter[6] + inter[7] * inter[7]);
            if (nf > 0.0001F) {
              inter[6] = inter[6] / nf;
              inter[7] = inter[7] / nf;
            }
            inter += cINTER_ENTRIES;
          }
        }
      }
    }
    {
      const float _0F = 0.0F;

      if ((scale != 0.0F) || (seq_wt != 0.0F)) {
        for (a = 0; a < n1; a++) {
          float* i1 = inter1 + cINTER_ENTRIES * a;
          for (b = 0; b < n2; b++) {
            float* i2 = inter2 + cINTER_ENTRIES * b;
            float sm[cINTER_ENTRIES], comp1, comp2, comp3 = 1.0F;
            float score;
            int c;
            for (c = 0; c < (cINTER_ENTRIES - 1); c += 2) {
              if (((i1[c] == _0F) && (i1[c + 1] == _0F)) ||
                  ((i2[c] == _0F) && (i2[c + 1] == _0F))) {
                /* handle glycine case */
                sm[c] = 1.0F;
                sm[c + 1] = 1.0F;
              } else {
                sm[c] = i1[c] + i2[c];
                sm[c + 1] = i1[c + 1] + i2[c + 1];
              }
            }
            comp1 = (float)((sqrt(sm[0] * sm[0] + sm[1] * sm[1]) +
                             sqrt(sm[2] * sm[2] + sm[3] * sm[3])) *
                            0.25);
            comp2 = (float)((sqrt(sm[4] * sm[4] + sm[5] * sm[5]) +
                             sqrt(sm[6] * sm[6] + sm[7] * sm[7])) *
                            0.25);
            score = scale * (comp1 * comp2 - base);
            if (coord_wt != 0.0) {
              float diff = (float)diff3f(i1 + 8, i2 + 8);
              comp3 = (float)-log(diff / rms_exp);
              score = (1 - coord_wt) * score + coord_wt * comp3 * scale;
            }
            match->mat[a][b] = seq_wt * match->mat[a][b] + score;
          }
        }
      }
    }
  }
  FreeP(inter1);
  FreeP(inter2);
  FreeP(v_ca);
  return 1;
}
int SelectorGetPairIndices(PyMOLGlobals* G, int sele1, int state1, int sele2,
                           int state2, int mode, float cutoff, float h_angle,
                           int** indexVLA, ObjectMolecule*** objVLA) {
  CSelector* I = G->Selector;
  float dist;
  int a1, a2;
  int at1, at2;
  CoordSet *cs1, *cs2;
  ObjectMolecule *obj1, *obj2;
  int idx1, idx2;
  int a;
  int dist_cnt = 0;
  float dir[3];
  float v1[3], v2[3];
  int flag;
  float angle_cutoff = 0.0;

  if (mode == 1) {
    angle_cutoff = (float)cos(PI * h_angle / 180.0);
  }

  if (state1 < 0) state1 = 0;
  if (state2 < 0) state2 = 0;

  if (state1 != state2) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, state1, -1);
  }
  if (cutoff < 0) cutoff = 1000.0;

  auto vla =
      SelectorGetInterstateVector(G, sele1, state1, sele2, state2, cutoff);
  const int c = vla.size() / 2;

  (*indexVLA) = VLAlloc(int, 1000);
  (*objVLA) = VLAlloc(ObjectMolecule*, 1000);

  for (a = 0; a < c; a++) {
    a1 = vla[a * 2];
    a2 = vla[a * 2 + 1];

    if (a1 != a2) {
      at1 = I->Table[a1].atom;
      at2 = I->Table[a2].atom;

      obj1 = I->Obj[I->Table[a1].model];
      obj2 = I->Obj[I->Table[a2].model];

      if (state1 < obj1->NCSet && state2 < obj2->NCSet) {
        cs1 = obj1->CSet[state1];
        cs2 = obj2->CSet[state2];
        if (cs1 && cs2) {
          idx1 = cs1->atmToIdx(at1);
          idx2 = cs2->atmToIdx(at2);

          if ((idx1 >= 0) && (idx2 >= 0)) {
            subtract3f(cs1->coordPtr(idx1), cs2->coordPtr(idx2), dir);
            dist = (float)length3f(dir);
            if (dist > R_SMALL4) {
              float dist_1 = 1.0F / dist;
              scale3f(dir, dist_1, dir);
            }
            if (dist < cutoff) {
              if (mode == 1) { /* coarse hydrogen bonding assessment */
                flag = false;
                if (ObjectMoleculeGetAvgHBondVector(obj1, at1, state1, v1,
                                                    nullptr) > 0.3)
                  if (dot_product3f(v1, dir) < -angle_cutoff) flag = true;
                if (ObjectMoleculeGetAvgHBondVector(obj2, at2, state2, v2,
                                                    nullptr) > 0.3)
                  if (dot_product3f(v2, dir) > angle_cutoff) flag = true;
              } else
                flag = true;

              if (flag) {
                VLACheck((*objVLA), ObjectMolecule*, dist_cnt + 1);
                VLACheck((*indexVLA), int, dist_cnt + 1);
                (*objVLA)[dist_cnt] = obj1;
                (*indexVLA)[dist_cnt] = at1;
                dist_cnt++;
                (*objVLA)[dist_cnt] = obj2;
                (*indexVLA)[dist_cnt] = at2;
                dist_cnt++;
              }
            }
          }
        }
      }
    }
  }

  VLASize((*objVLA), ObjectMolecule*, dist_cnt);
  VLASize((*indexVLA), int, dist_cnt);
  dist_cnt = dist_cnt / 2;
  return (dist_cnt);
}
int SelectorCreateAlignments(PyMOLGlobals* G, int* pair, int sele1, int* vla1,
                             int sele2, int* vla2, const char* name1,
                             const char* name2, int identical,
                             int atomic_input) {
  CSelector* I = G->Selector;
  int *flag1 = nullptr, *flag2 = nullptr;
  int* p;
  int i, np;
  int cnt;
  int mod1, mod2;                        /* model indexes */
  int at1, at2, at1a, at2a;              /* atoms indexes */
  int vi1, vi2;                          /* vla indexes */
  int index1, index2;                    /* indices in the selection array */
  AtomInfoType *ai1, *ai2, *ai1a, *ai2a; /* atom information pointers */
  ObjectMolecule *obj1, *obj2;
  int cmp;
  PRINTFD(G, FB_Selector)
  " %s-DEBUG: entry.\n", __func__ ENDFD cnt = 0;
  /* number of pairs of atoms */
  np = VLAGetSize(pair) / 2;
  if (np) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates,
                        -1); /* unnecessary? */
    /* flags initialized to false */
    flag1 = pymol::calloc<int>(I->Table.size());
    flag2 = pymol::calloc<int>(I->Table.size());

    /* we need to create two selection arrays: for the matched
     * atoms in the original selections */
    p = pair;
    for (i = 0; i < np;
         i++) { /* iterate through all pairs of matched residues */
      vi1 = *(p++);
      vi2 = *(p++);

      /* find positions in the selection arrays */
      /* SelectorGetResidueVLA returns a VLA with 3 entries per atom,
       * index0=model, index1=atom#, hence: *3, and *3+1 */
      mod1 = vla1[vi1 * 3];
      at1 = vla1[vi1 * 3 + 1];

      mod2 = vla2[vi2 * 3];
      at2 = vla2[vi2 * 3 + 1];

      PRINTFD(G, FB_Selector)
      " S.C.A.-DEBUG: mod1 %d at1 %d mod2 %d at2 %d\n", mod1, at1, mod2,
          at2 ENDFD obj1 = I->Obj[mod1];
      obj2 = I->Obj[mod2];

      ai1 = obj1->AtomInfo + at1;
      ai2 = obj2->AtomInfo + at2;
      at1a = at1;
      at2a = at2;
      ai1a = ai1;
      ai2a = ai2;

      if (atomic_input) {
        index1 = SelectorGetObjAtmOffset(I, obj1, at1a);
        index2 = SelectorGetObjAtmOffset(I, obj2, at2a);
        flag1[index1] = true;
        flag2[index2] = true;
        cnt++;
      } else {
        /* search back to first atom in residue */
        while (at1a > 0 && AtomInfoSameResidue(G, ai1a, ai1a - 1)) {
          ai1a--;
          at1a--;
        }
        while (at2a > 0 && AtomInfoSameResidue(G, ai2a, ai2a - 1)) {
          ai2a--;
          at2a--;
        }

        while (1) { /* match up all matching atom names in each residue */
          cmp = AtomInfoNameOrder(G, ai1a, ai2a);
          if (cmp == 0) { /* atoms match */
            index1 = SelectorGetObjAtmOffset(I, obj1, at1a);
            index2 = SelectorGetObjAtmOffset(I, obj2, at2a);

            PRINTFD(G, FB_Selector)
            " S.C.A.-DEBUG: compare %s %s %d\n", LexStr(G, ai1a->name),
                LexStr(G, ai2a->name),
                cmp ENDFD PRINTFD(G,
                                  FB_Selector) " S.C.A.-DEBUG: entry %d %d\n",
                ai1a -> selEntry,
                ai2a->selEntry ENDFD if ((index1 >= 0) && (index2 >= 0)) {
              if (SelectorIsMember(G, ai1a->selEntry, sele1) &&
                  SelectorIsMember(G, ai2a->selEntry, sele2)) {
                if ((!identical) || (ai1a->resn == ai2a->resn)) {
                  flag1[index1] = true;
                  flag2[index2] = true;
                  cnt++;
                }
              }
            }
            at1a++;
            at2a++;
          } else if (cmp < 0) { /* 1 is before 2 */
            at1a++;
          } else if (cmp > 0) { /* 1 is after 2 */
            at2a++;
          }
          if (at1a >= obj1->NAtom) break;
          if (at2a >= obj2->NAtom) break;
          ai1a = obj1->AtomInfo + at1a;
          ai2a = obj2->AtomInfo + at2a;
          /* make sure we're still in the same residue */
          if (!AtomInfoSameResidue(G, ai1a, ai1)) break;
          if (!AtomInfoSameResidue(G, ai2a, ai2)) break;
        }
      }
    }
    if (cnt) {
      SelectorEmbedSelection(G, flag1, name1, nullptr, false, -1);
      SelectorEmbedSelection(G, flag2, name2, nullptr, false, -1);
    }
    FreeP(flag1);
    FreeP(flag2);
  }
  PRINTFD(G, FB_Selector)
  " %s-DEBUG: exit, cnt = %d.\n", __func__, cnt ENDFD return cnt;
}
