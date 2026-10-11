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

// Neighbor and nearest-atom queries, pseudoatoms, sorting and geometry.

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/base/feedback.h"
#include "pymol/executive/executive.h"

#include "pymol/objects/object_molecule_internal.h"

#define cMaxOther 6

typedef struct {
  int n_cyclic_arom, cyclic_arom[cMaxOther];
  int n_arom, arom[cMaxOther];
  int n_high_val, high_val[cMaxOther];
  int n_cyclic, cyclic[cMaxOther];
  int n_planer, planer[cMaxOther];
  int n_rest, rest[cMaxOther];
  int score;
} OtherRec;

static int populate_other(OtherRec* other, int at, const AtomInfoType* ai,
                          const BondType* bd, const ObjectMolecule* obj) {
  int five_cycle = false;
  int six_cycle = false;

  {
    // don't get bogged down with structures that have unreasonable connectivity
    const int ESCAPE_MAX = 500;
    int escape_count = ESCAPE_MAX;

    auto const mem0_atm = bd->index[0];
    auto const mem1_atm = bd->index[1];

    for (auto const& mem2 : AtomNeighbors(obj, mem1_atm)) {
      if (mem2.atm == mem0_atm) {
        continue;
      }

      for (auto const& mem3 : AtomNeighbors(obj, mem2.atm)) {
        if (mem3.atm == mem1_atm) {
          continue;
        }

        for (auto const& mem4 : AtomNeighbors(obj, mem3.atm)) {
          if (mem4.atm == mem2.atm || mem4.atm == mem1_atm ||
              mem4.atm == mem0_atm) {
            continue;
          }

          for (auto const& mem5 : AtomNeighbors(obj, mem4.atm)) {
            if (!(escape_count--)) {
              goto escape;
            }

            if (mem5.atm == mem3.atm || mem5.atm == mem2.atm ||
                mem5.atm == mem1_atm) {
              continue;
            }

            if (mem5.atm == mem0_atm) { /* five-cycle */
              five_cycle = true;

              if (six_cycle) {
                goto escape;
              }
            }

            for (auto const& mem6 : AtomNeighbors(obj, mem5.atm)) {
              if (mem6.atm == mem4.atm || mem6.atm == mem3.atm ||
                  mem6.atm == mem2.atm || mem6.atm == mem1_atm) {
                continue;
              }

              if (mem6.atm == mem0_atm) { /* six-cycle */
                six_cycle = true;

                if (five_cycle) {
                  goto escape;
                }

                break;
              }
            }
          }
        }
      }
    }
  }
escape:

  if (bd->order == 4) { /* aromatic */
    if ((five_cycle || six_cycle) && (other->n_cyclic_arom < cMaxOther)) {
      other->cyclic_arom[other->n_cyclic_arom++] = at;
      if (five_cycle && six_cycle)
        other->score += 34;
      else if (five_cycle)
        other->score += 33;
      else
        other->score += 32;
      return 1;
    } else if (other->n_arom < cMaxOther) {
      other->arom[other->n_arom++] = at;
      other->score += 64;
      return 1;
    }
  }
  if (bd->order > 1) {
    if (other->n_high_val < cMaxOther) {
      other->high_val[other->n_high_val++] = at;
      other->score += 16;
      return 1;
    }
  }
  if (five_cycle || six_cycle) {
    if (other->n_cyclic < cMaxOther) {
      other->cyclic[other->n_cyclic++] = at;
      other->score += 8;
      return 1;
    }
  }
  if (ai->geom == cAtomInfoPlanar) {
    if (other->n_planer < cMaxOther) {
      other->planer[other->n_planer++] = at;
      other->score += 4;
      return 1;
    }
  }
  if (other->n_rest < cMaxOther) {
    other->rest[other->n_rest++] = at;
    other->score += 1;
    return 1;
  }
  return 0;
}

static int append_index(int* result, int offset, int a1, int a2, int score,
                        int ar_count) {
  int c;
  c = result[a1];
  while (c < offset) {
    if (result[c] == a2) { /* already entered */
      if (result[c + 1] < score) {
        result[c + 1] = score;
        result[c + 2] = ar_count;
      }
      return offset;
    }
    c += 3;
  }
  result[offset++] = a2;
  result[offset++] = score;
  result[offset++] = ar_count;
  return offset;
}

int ObjectMoleculeAddPseudoatom(ObjectMolecule* I, int sele_index,
                                const char* name, const char* resn,
                                const char* resi, const char* chain,
                                const char* segi, const char* elem, float vdw,
                                int hetatm, float b, float q, const char* label,
                                const float* pos, int color, int state,
                                int mode, int quiet) {
  PyMOLGlobals* G = I->G;
  int start_state = 0, stop_state = 0;
  int extant_only = false;
  int ai_merged = false;
  float pos_array[3] = {0.0F, 0.0F, 0.0F};
  int ok = true;

  pymol::vla<AtomInfoType> atInfo(1);
  AtomInfoType* ai = atInfo.data();

  // FIXME this should be cStateCurrent
  if (state == cStateAll) {
    state = I->getCurrentState();
  }

  if (state >= 0) { /* specific state */
    start_state = state;
    stop_state = state + 1;
  } else { /* all states */
    if (sele_index >= 0) {
      start_state = 0;
      stop_state = SelectorCountStates(G, sele_index);

      // Here, -3 does not mean cSelectorUpdateTableEffectiveStates. It means
      // all states present in `sele_index` AND `I`. State != -3 means all
      // states in `sele_index`, so it would create new states in `I` if they
      // don't exist yet.
      if (state == -3) extant_only = true;
    } else {
      start_state = 0;
      stop_state = ExecutiveCountStates(G, cKeywordAll);
      if (stop_state < 1) stop_state = 1;
    }
  }
  {
    /* match existing properties of the old atom */
    ai->setResi(resi);
    ai->hetatm = hetatm;
    ai->geom = cAtomInfoNone;
    ai->q = q;
    ai->b = b;
    ai->chain = LexIdx(G, chain);
    ai->segi = LexIdx(G, segi);
    ai->resn = LexIdx(G, resn);
    ai->name = LexIdx(G, name);
    strcpy(ai->elem, elem);
    ai->id = -1;
    ai->rank = -1;
    if (vdw >= 0.0F)
      ai->vdw = vdw;
    else
      ai->vdw = 1.0F;
    if (label[0]) {
      ai->label = LexIdx(G, label);
      ai->visRep = cRepLabelBit;
    } else {
      ai->visRep = RepGetAutoShowMask(G);
    }

    ai->flags |= cAtomFlag_inorganic;  // suppress auto_show_classified

    if (color < 0) {
      AtomInfoAssignColors(I->G, ai);
      if ((ai->elem[0] == 'C') && (ai->elem[1] == 0))
        /* carbons are always colored according to the object color */
        ai->color = I->Color;
    } else {
      ai->color = color;
    }
    AtomInfoAssignParameters(I->G, ai);
    AtomInfoUniquefyNames(I->G, I->AtomInfo, I->NAtom, ai, nullptr, 1);
    if (!quiet) {
      PRINTFB(G, FB_ObjectMolecule, FB_Actions)
      " ObjMol: created %s/%s/%s/%s`%d%c/%s\n", I->Name, LexStr(G, ai->segi),
          LexStr(G, ai->chain), LexStr(G, ai->resn), ai->resv,
          ai->getInscode(true), LexStr(G, ai->name) ENDFB(G);
    }
  }

  CoordSet* cset = CoordSetNew(G);
  cset->NIndex = 1;
  cset->Coord = pymol::vla<float>(3);
  cset->Obj = I;
  cset->enumIndices();

  for (state = start_state; state < stop_state; state++) {
    if ((extant_only && (state < I->NCSet) && I->CSet[state]) || !extant_only) {
      if (sele_index >= 0) {
        ObjectMoleculeOpRec op;
        ObjectMoleculeOpRecInit(&op);
        op.code = OMOP_CSetSumVertices;
        op.cs1 = state;

        ExecutiveObjMolSeleOp(I->G, sele_index, &op);

        if (op.i1) {
          float factor = 1.0F / op.i1;
          scale3f(op.v1, factor, pos_array);
          pos = pos_array;

          if (vdw < 0.0F) {
            switch (mode) {
              case 1:
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_CSetMaxDistToPt;
                copy3f(pos_array, op.v1);
                op.cs1 = state;
                ExecutiveObjMolSeleOp(I->G, sele_index, &op);
                vdw = op.f1;
                break;
              case 2:
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_CSetSumSqDistToPt;
                copy3f(pos_array, op.v1);
                op.cs1 = state;
                ExecutiveObjMolSeleOp(I->G, sele_index, &op);
                vdw = sqrt1f(op.d1 / op.i1);
                break;
              case 0:
              default:
                vdw = 0.5F;
                break;
            }
            if (vdw >= 0.0F)
              ai->vdw =
                  vdw; /* NOTE: only uses vdw from first state selection... */
          }
        } else {
          pos = nullptr; /* skip this state */
        }
      } else if (!pos) {
        SceneGetCenter(I->G, pos_array);
        pos = pos_array;
      }

      if (pos) { /* only add coordinate to state if we have position for it */

        float* coord = cset->Coord.data();

        copy3f(pos, coord);

        if (!ai_merged) {
          if (ok)
            ok &= ObjectMoleculeMerge(I, std::move(atInfo), cset, false,
                                      cAIC_AllMask,
                                      true); /* NOTE: will release atInfo */
          if (ok) ok &= ObjectMoleculeExtendIndices(I, -1);
          ai_merged = true;
        }
        if (state >= I->NCSet) {
          VLACheck(I->CSet, CoordSet*, state);
          I->NCSet = state + 1;
        }
        if (!I->CSet[state]) {
          /* new coordinate set */
          I->CSet[state] = CoordSetCopy(cset);
        } else {
          /* merge coordinate set */
          if (ok) ok &= CoordSetMerge(I, I->CSet[state], cset);
        }
      }
    }
  }

  delete cset;

  if (ai_merged) {
    if (ok) ok &= ObjectMoleculeSort(I);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
    I->invalidate(cRepAll, cRepInvAtoms, -1);
  } else {
    VLAFreeP(atInfo);
  }
  return (ok);
}

int* ObjectMoleculeGetPrioritizedOtherIndexList(ObjectMolecule* I,
                                                CoordSet* cs) {
  int a, b;
  int b1, b2, a1, a2, a3;
  OtherRec* o;
  OtherRec* other = pymol::calloc<OtherRec>(cs->NIndex);
  int* result = nullptr;
  int offset;
  int n_alloc = 0;
  const BondType* bd;
  int ok = true;

  CHECKOK(ok, other);

  bd = I->Bond;
  for (a = 0; ok && a < I->NBond; a++) {
    b1 = bd->index[0];
    b2 = bd->index[1];
    a1 = cs->atmToIdx(b1);
    a2 = cs->atmToIdx(b2);
    if ((a1 >= 0) && (a2 >= 0)) {
      n_alloc += populate_other(other + a1, a2, I->AtomInfo + b2, bd, I);
      n_alloc += populate_other(other + a2, a1, I->AtomInfo + b1, bd, I);
    }
    bd++;
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    n_alloc = 3 * (n_alloc + cs->NIndex);
    o = other;
    result = pymol::malloc<int>(n_alloc);
    CHECKOK(ok, result);
  }
  if (ok) {
    for (a = 0; a < cs->NIndex; a++) {
      result[a] = -1;
    }
    offset = cs->NIndex;
    bd = I->Bond;
  }
  for (a = 0; ok && a < I->NBond; a++) {
    b1 = bd->index[0];
    b2 = bd->index[1];
    a1 = cs->atmToIdx(b1);
    a2 = cs->atmToIdx(b2);
    if ((a1 >= 0) && (a2 >= 0)) {
      if (result[a1] < 0) {
        o = other + a1;
        result[a1] = offset;
        for (b = 0; b < o->n_cyclic_arom; b++) {
          a3 = o->cyclic_arom[b];
          offset =
              append_index(result, offset, a1, a3, 128 + other[a3].score, 1);
        }
        for (b = 0; b < o->n_arom; b++) {
          a3 = o->arom[b];
          offset =
              append_index(result, offset, a1, a3, 64 + other[a3].score, 1);
        }
        for (b = 0; b < o->n_high_val; b++) {
          a3 = o->high_val[b];
          offset =
              append_index(result, offset, a1, a3, 16 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_cyclic; b++) {
          a3 = o->cyclic[b];
          offset = append_index(result, offset, a1, a3, 8 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_planer; b++) {
          a3 = o->planer[b];
          offset = append_index(result, offset, a1, a3, 2 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_rest; b++) {
          a3 = o->rest[b];
          offset = append_index(result, offset, a1, a3, 1 + other[a3].score, 0);
        }
        result[offset++] = -1;
      }

      if (result[a2] < 0) {
        o = other + a2;
        result[a2] = offset;
        for (b = 0; b < o->n_cyclic_arom; b++) {
          a3 = o->cyclic_arom[b];
          offset =
              append_index(result, offset, a2, a3, 128 + other[a3].score, 1);
        }
        for (b = 0; b < o->n_arom; b++) {
          a3 = o->arom[b];
          offset =
              append_index(result, offset, a2, a3, 64 + other[a3].score, 1);
        }
        for (b = 0; b < o->n_high_val; b++) {
          a3 = o->high_val[b];
          offset =
              append_index(result, offset, a2, a3, 16 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_cyclic; b++) {
          a3 = o->cyclic[b];
          offset = append_index(result, offset, a2, a3, 8 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_planer; b++) {
          a3 = o->planer[b];
          offset = append_index(result, offset, a2, a3, 2 + other[a3].score, 0);
        }
        for (b = 0; b < o->n_rest; b++) {
          a3 = o->rest[b];
          offset = append_index(result, offset, a2, a3, 1 + other[a3].score, 0);
        }
        result[offset++] = -1;
      }
    }
    bd++;
    ok &= !I->G->Interrupt;
  }
  FreeP(other);
  return result;
}

int ObjectMoleculeGetNearestBlendedColor(ObjectMolecule* I, const float* point,
                                         float cutoff, int state, float* dist,
                                         float* color, int sub_vdw) {
  int result = -1;
  float tot_weight = 0.0F;
  float cutoff2 = cutoff * cutoff;
  float nearest = -1.0F;

  color[0] = 0.0F;
  color[1] = 0.0F;
  color[2] = 0.0F;

  assert(state != -1 /* all states */);
  auto* cs = I->getCoordSet(state);

  {
    if (cs) {
      MapType* map;
      CoordSetUpdateCoord2IdxMap(cs, cutoff);
      if (sub_vdw) {
        cutoff -= MAX_VDW;
        cutoff2 = cutoff * cutoff;
      }
      nearest = cutoff2;
      if ((map = cs->Coord2Idx)) {
        int a, b, c, d, e, f, j;
        float test;
        float* v;
        MapLocus(map, point, &a, &b, &c);
        for (d = a - 1; d <= a + 1; d++)
          for (e = b - 1; e <= b + 1; e++)
            for (f = c - 1; f <= c + 1; f++) {
              j = *(MapFirst(map, d, e, f));
              while (j >= 0) {
                v = cs->coordPtr(j);
                test = diffsq3f(v, point);
                if (sub_vdw) {
                  test = sqrt1f(test);
                  test -= I->AtomInfo[cs->IdxToAtm[j]].vdw;
                  if (test < 0.0F) test = 0.0F;
                  test = test * test;
                }
                if (test < cutoff2) {
                  float weight = cutoff - sqrt1f(test);
                  const float* at_col =
                      ColorGet(I->G, I->AtomInfo[cs->IdxToAtm[j]].color);
                  color[0] += at_col[0] * weight;
                  color[1] += at_col[1] * weight;
                  color[2] += at_col[2] * weight;
                  tot_weight += weight;
                }
                if (test <= nearest) {
                  result = j;
                  nearest = test;
                }
                j = MapNext(map, j);
              }
            }
      } else {
        int j;
        float test;
        const float* v = cs->Coord.data();
        for (j = 0; j < cs->NIndex; j++) {
          test = diffsq3f(v, point);
          if (sub_vdw) {
            test = sqrt1f(test);
            test -= I->AtomInfo[cs->IdxToAtm[j]].vdw;
            if (test < 0.0F) test = 0.0F;
            test = test * test;
          }
          if (test < cutoff2) {
            float weight = cutoff - sqrt1f(test);
            const float* at_col =
                ColorGet(I->G, I->AtomInfo[cs->IdxToAtm[j]].color);
            color[0] += at_col[0] * weight;
            color[1] += at_col[1] * weight;
            color[2] += at_col[2] * weight;
            tot_weight += weight;
          }
          if (test <= nearest) {
            result = j;
            nearest = test;
          }
          v += 3;
        }
      }
      if (result >= 0) result = cs->IdxToAtm[result];
    }
  }
  if (dist) {
    if (result >= 0) {
      *dist = sqrt1f(nearest);
      if (tot_weight > 0.0F) {
        color[0] /= tot_weight;
        color[1] /= tot_weight;
        color[2] /= tot_weight;
      }
    } else {
      *dist = -1.0F;
    }
  }
  return result;
}

int ObjectMoleculeGetNearestAtomIndex(ObjectMolecule* I, const float* point,
                                      float cutoff, int state, float* dist) {
  int result = -1;
  float nearest = -1.0F;

  assert(state != -1 /* all states */);
  auto* cs = I->getCoordSet(state);

  {
    if (cs) {
      MapType* map;
      CoordSetUpdateCoord2IdxMap(cs, cutoff);
      nearest = cutoff * cutoff;
      if ((map = cs->Coord2Idx)) {
        int a, b, c, d, e, f, j;
        float test;
        float* v;
        MapLocus(map, point, &a, &b, &c);
        for (d = a - 1; d <= a + 1; d++)
          for (e = b - 1; e <= b + 1; e++)
            for (f = c - 1; f <= c + 1; f++) {
              j = *(MapFirst(map, d, e, f));
              while (j >= 0) {
                v = cs->coordPtr(j);
                test = diffsq3f(v, point);
                if (test <= nearest) {
                  result = j;
                  nearest = test;
                }
                j = MapNext(map, j);
              }
            }
      } else {
        int j;
        float test;
        const float* v = cs->Coord.data();
        for (j = 0; j < cs->NIndex; j++) {
          test = diffsq3f(v, point);
          if (test <= nearest) {
            result = j;
            nearest = test;
          }
          v += 3;
        }
      }
      if (result >= 0) result = cs->IdxToAtm[result];
    }
  }
  if (dist) {
    if (result >= 0) {
      *dist = sqrt1f(nearest);
    } else {
      *dist = -1.0F;
    }
  }
  return result;
}

int ObjectMoleculeGetPrioritizedOther(const int* other, int a1, int a2,
                                      int* double_sided) {
  int a3 = -1;
  int lvl = -1, ck, ck_lvl;
  int offset;
  int ar_count = 0;

  a3 = -1;
  lvl = -1;
  if (a1 >= 0) {
    offset = other[a1];
    if (offset >= 0) {
      while (1) {
        ck = other[offset];
        if (ck != a2) {
          if (ck >= 0) {
            ck_lvl = other[offset + 1];
            if (ck_lvl > lvl) {
              a3 = ck;
              lvl = ck_lvl;
            }
            ar_count += other[offset + 2];
          } else
            break;
        }
        offset += 3;
      }
    }
  }
  if (a2 >= 0) {
    offset = other[a2];
    if (offset >= 0) {
      while (1) {
        ck = other[offset];
        if (ck != a1) {
          if (ck >= 0) {
            ck_lvl = other[offset + 1];
            if (ck_lvl > lvl) {
              a3 = ck;
              lvl = ck_lvl;
            }
            ar_count += other[offset + 2];
          } else
            break;
        }
        offset += 3;
      }
    }
  }

  if (double_sided) {
    if (ar_count == 4)
      *double_sided = true;
    else
      *double_sided = false;
  }
  return a3;
}

void ObjectMoleculeFixChemistry(ObjectMolecule* I, int sele1, int sele2,
                                int invalidate) {
  int b;
  int flag = false;
  int s1, s2;
  AtomInfoType *ai1, *ai2;
  int order;
  BondType* bond = I->Bond.data();
  for (b = 0; b < I->NBond; b++) {
    flag = false;
    ai1 = I->AtomInfo + bond->index[0];
    ai2 = I->AtomInfo + bond->index[1];
    s1 = ai1->selEntry;
    s2 = ai2->selEntry;

    if ((SelectorIsMember(I->G, s1, sele1) &&
         SelectorIsMember(I->G, s2, sele2)) ||
        (SelectorIsMember(I->G, s2, sele1) &&
         SelectorIsMember(I->G, s1, sele2))) {
      order = -1;
      if (strlen(LexStr(I->G, ai1->resn)) <
          4) { /* Standard disconnected PDB residue */
        if (AtomInfoSameResidue(I->G, ai1, ai2)) {
          assign_pdb_known_residue(I->G, ai1, ai2, &order);
        }
      }
      if (order > 0) {
        bond->order = order;
        ai1->chemFlag = false;
        ai2->chemFlag = false;
        flag = true;
      } else if (invalidate) {
        ai1->chemFlag = false;
        ai2->chemFlag = false;
        flag = true;
      }
    }
    bond++;
  }
  if (flag) {
    I->invalidate(cRepAll, cRepInvAll, -1);
    SceneChanged(I->G);
  }
}
float ObjectMoleculeGetMaxVDW(ObjectMolecule* I) {
  float max_vdw = 0.0F;
  int a;
  const AtomInfoType* ai;
  if (I->NAtom) {
    ai = I->AtomInfo;
    for (a = 0; a < I->NAtom; a++) {
      if (max_vdw < ai->vdw) max_vdw = ai->vdw;
      ai++;
    }
  }
  return (max_vdw);
}
/**
 * Sort the ObjectMolecule::AtomInfo and ObjectMolecule::Bond arrays and adjust
 * IdxToAtm/AtmToIdx in all coordiante sets.
 *
 * This function has no effect on discrete objects.
 */
int ObjectMoleculeSort(ObjectMolecule* I) { /* sorts atoms and bonds */
  int* index;
  int* outdex = nullptr;
  int a, b;
  int ok = true;
  if (!I->DiscreteFlag) { /* currently, discrete objects are never sorted */
    int already_in_order = true;
    int i_NAtom = I->NAtom;
    index = AtomInfoGetSortedIndex(I->G, I, I->AtomInfo, i_NAtom, &outdex);
    CHECKOK(ok, index);
    if (ok) {
      for (a = 0; a < i_NAtom; a++) {
        if (index[a] != a) {
          already_in_order = false;
          break;
        }
      }
    }
    if (ok && !already_in_order) { /* if we aren't already in perfect order */

      for (a = 0; a < I->NBond; a++) { /* bonds */
        I->Bond[a].index[0] = outdex[I->Bond[a].index[0]];
        I->Bond[a].index[1] = outdex[I->Bond[a].index[1]];
      }

      for (a = -1; a < I->NCSet; a++) { /* coordinate set mapping */
        auto* cs = (a < 0) ? I->CSTmpl : I->CSet[a];

        if (cs) {
          int cs_NIndex = cs->NIndex;
          int* cs_IdxToAtm = cs->IdxToAtm.data();
          for (b = 0; b < cs_NIndex; b++)
            cs_IdxToAtm[b] = outdex[cs_IdxToAtm[b]];
        }
      }

      I->updateAtmToIdx();

      ExecutiveUniqueIDAtomDictInvalidate(I->G);

      pymol::vla<AtomInfoType> atInfo(i_NAtom);
      CHECKOK(ok, atInfo);
      if (ok) {
        for (a = 0; a < i_NAtom; a++)
          atInfo[a] = std::move(I->AtomInfo[index[a]]);
        I->AtomInfo = std::move(atInfo);
      }
    }
    AtomInfoFreeSortedIndexes(I->G, &index, &outdex);
    if (ok) {
      UtilSortInPlace(I->G, I->Bond.data(), I->NBond, sizeof(BondType),
                      (UtilOrderFn*)BondInOrder);
      /* sort...important! */
      I->invalidate(cRepAll, cRepInvAtoms, -1); /* important */
    }
  }
  return ok;
}

int ObjectMoleculeSetGeometry(PyMOLGlobals* G, ObjectMolecule* I, int sele,
                              int geom, int valence) {
  int count = 0;
  for (int a = 0; a < I->NAtom; ++a) {
    auto s = I->AtomInfo[a].selEntry;

    if (SelectorIsMember(G, s, sele)) {
      auto& ai = I->AtomInfo[a];
      ai.geom = geom;
      ai.valence = valence;
      ai.chemFlag = true;
      count++;
      break;
    }
  }
  return count;
}
