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

// Inferring chemistry: valences, bond orders, geometry and hydrogen bonds.

#include "pymol/objects/object_molecule.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_molecule_internal.h"

/* find sets of atoms with identical skeletons */

struct match_info {
  AtomInfoType* ai_a;
  AtomInfoType* ai_b;
  BondType* bi_a;
  BondType* bi_b;
  const int *nbr_a, *nbr_b;
  int* matched;

  // values: -1 .. 2
  using mark_type = signed char;

  std::vector<mark_type> atom_mark_a;
  std::vector<mark_type> atom_mark_b;
  std::vector<mark_type> bond_mark_a;
  std::vector<mark_type> bond_mark_b;
};

#define recmat3(x, y, z)                                  \
  (recursive_match(u_a[0], u_b[x], x_a[0], x_b[x], mi) && \
   recursive_match(u_a[1], u_b[y], x_a[1], x_b[y], mi) && \
   recursive_match(u_a[2], u_b[z], x_a[2], x_b[z], mi))

#define recmat4(w, x, y, z)                               \
  (recursive_match(u_a[0], u_b[w], x_a[0], x_b[w], mi) && \
   recursive_match(u_a[1], u_b[x], x_a[1], x_b[x], mi) && \
   recursive_match(u_a[2], u_b[y], x_a[2], x_b[y], mi) && \
   recursive_match(u_a[3], u_b[z], x_a[3], x_b[z], mi))

static void undo_match(int* start, match_info* mi) {
  /* remove the matched atom set */
  int* match = mi->matched;
  while (match > start) {
    int a = match[-4];
    int b = match[-3];
    int aa = match[-2];
    int bb = match[-1];
    mi->atom_mark_a[a] = false;
    mi->atom_mark_b[b] = false;
    mi->bond_mark_a[aa] = false;
    mi->bond_mark_b[bb] = false;
    match -= 4;
  }
  mi->matched = start;
}

static int recursive_match(int a, int b, int aa, int bb, match_info* mi) {
  if (mi->atom_mark_a[a] && mi->atom_mark_b[b]) {
    if (aa >= 0 && bb >= 0 && !(mi->bond_mark_a[aa] || mi->bond_mark_b[bb])) {
      /* note bond */
      mi->atom_mark_a[a] = true;
      mi->atom_mark_b[b] = true;
      mi->bond_mark_a[aa] = true;
      mi->bond_mark_b[bb] = true;
      mi->matched[0] = a; /* atoms */
      mi->matched[1] = b;
      mi->matched[2] = aa; /* bonds */
      mi->matched[3] = bb;
      mi->matched += 4;
    }
    return true;
  } else if (mi->atom_mark_a[a] != mi->atom_mark_b[b])
    return false;
  else if (mi->ai_a[a].protons != mi->ai_b[b].protons)
    return false;
  else {
    int ni_a = mi->nbr_a[a];
    int ni_b = mi->nbr_b[b];
    int num_a = mi->nbr_a[ni_a++];
    int num_b = mi->nbr_b[ni_b++];

    if ((num_a != num_b) ||
        (num_a > 4)) /* must have the same number of neighbors (four or less) */
      return false;
    else {
      int match_flag = false;
      int u_a[4], u_b[4], x_a[4], x_b[4];
      int n_u_a = 0;
      int n_u_b = 0;
      int* before = mi->matched;
      mi->atom_mark_a[a] = true;
      mi->atom_mark_b[b] = true;
      mi->bond_mark_a[aa] = true;
      mi->bond_mark_b[bb] = true;
      mi->matched[0] = a; /* atoms */
      mi->matched[1] = b;
      mi->matched[2] = aa; /* bonds */
      mi->matched[3] = bb;
      mi->matched += 4;

      /* collect unvisited bonds */

      while (num_a--) {
        int i_a = mi->nbr_a[ni_a + 1];
        if (!mi->bond_mark_a[i_a]) {    /* bond not yet visited */
          u_a[n_u_a] = mi->nbr_a[ni_a]; /* atom index */
          x_a[n_u_a++] = i_a;
        }
        ni_a += 2;
      }

      while (num_b--) {
        int i_b = mi->nbr_b[ni_b + 1];
        if (!mi->bond_mark_b[i_b]) {
          u_b[n_u_b] = mi->nbr_b[ni_b];
          x_b[n_u_b++] = i_b;
        }
        ni_b += 2;
      }
      /* NOTE: implicitly relying upon C short-circuit evaluation */
      if (n_u_a == n_u_b) {
        if (!n_u_a)
          match_flag = true;
        else if (n_u_a == 1) {
          match_flag = recursive_match(u_a[0], u_b[0], x_a[0], x_b[0], mi);
        } else if (n_u_a == 2) {
          match_flag = (recursive_match(u_a[0], u_b[0], x_a[0], x_b[0], mi) &&
                        recursive_match(u_a[1], u_b[1], x_a[1], x_b[1], mi)) ||
                       (recursive_match(u_a[0], u_b[1], x_a[0], x_b[1], mi) &&
                        recursive_match(u_a[1], u_b[0], x_a[1], x_b[0], mi));
        } else if (n_u_a == 3) {
          match_flag = recmat3(0, 1, 2) || recmat3(0, 2, 1) ||
                       recmat3(1, 0, 2) || recmat3(1, 2, 0) ||
                       recmat3(2, 0, 1) || recmat3(2, 1, 0);
        } else if (n_u_a == 4) {
          match_flag = recmat4(0, 1, 2, 3) || recmat4(0, 2, 1, 3) ||
                       recmat4(1, 0, 2, 3) || recmat4(1, 2, 0, 3) ||
                       recmat4(2, 0, 1, 3) || recmat4(2, 1, 0, 3);
          if (!match_flag) {
            match_flag = recmat4(0, 1, 3, 2) || recmat4(0, 2, 3, 1) ||
                         recmat4(1, 0, 3, 2) || recmat4(1, 2, 3, 0) ||
                         recmat4(2, 0, 3, 1) || recmat4(2, 1, 3, 0);
            if (!match_flag) {
              match_flag = recmat4(0, 3, 1, 2) || recmat4(0, 3, 2, 1) ||
                           recmat4(1, 3, 0, 2) || recmat4(1, 3, 2, 0) ||
                           recmat4(2, 3, 0, 1) || recmat4(2, 3, 1, 0);
              if (!match_flag) {
                match_flag = recmat4(3, 0, 1, 2) || recmat4(3, 0, 2, 1) ||
                             recmat4(3, 1, 0, 2) || recmat4(3, 1, 2, 0) ||
                             recmat4(3, 2, 0, 1) || recmat4(3, 2, 1, 0);
              }
            }
          }
        }
      }
      if (!match_flag) { /* clean the match tree */
        undo_match(before, mi);
      }
      return match_flag;
    }
  }
}

int ObjectMoleculeXferValences(ObjectMolecule* Ia, int sele1, int sele2,
                               int target_state, ObjectMolecule* Ib, int sele3,
                               int source_state, int quiet) {
  int* matched = nullptr;
  int match_found = false;
  PyMOLGlobals* G = Ia->G;
  if (Ia == Ib) return false;

  {
    int max_match = Ia->NAtom + Ia->NBond;
    if (max_match < (Ib->NAtom + Ib->NBond))
      max_match = (Ib->NAtom + Ib->NBond);
    matched = pymol::calloc<int>(max_match * 4);
  }

  {
    int a, b;
    AtomInfoType* ai_a = Ia->AtomInfo.data();
    AtomInfoType* ai_b = Ib->AtomInfo.data();
    BondType* bi_a = Ia->Bond.data();
    BondType* bi_b = Ib->Bond.data();

    match_info mi;

    mi.atom_mark_a.resize(Ia->NAtom);
    mi.atom_mark_b.resize(Ib->NAtom);
    mi.bond_mark_a.resize(Ia->NBond);
    mi.bond_mark_b.resize(Ib->NBond);

    // confirm zero-initialization
    assert(std::none_of(mi.atom_mark_a.begin(), mi.atom_mark_a.end(),
                        [](bool m) { return m; }));

    mi.ai_a = ai_a;
    mi.ai_b = ai_b;
    mi.bi_a = bi_a;
    mi.bi_b = bi_b;
    mi.nbr_a = Ia->getNeighborArray();
    mi.nbr_b = Ib->getNeighborArray();
    mi.matched = matched;
    for (a = 0; a < Ia->NAtom; a++) {
      if (!mi.atom_mark_a[a]) {
        int a_entry = ai_a[a].selEntry;
        if (SelectorIsMember(G, a_entry, sele1) ||
            SelectorIsMember(G, a_entry, sele2)) {
          for (b = 0; b < Ib->NAtom; b++) {
            if (SelectorIsMember(G, ai_b[b].selEntry, sele3)) {
              if (recursive_match(a, b, -1, -1, &mi)) {
                int* match = mi.matched;
                match_found = true;
                /* now graft bond orders for bonds in matching atoms */

                while (match > matched) {
                  int at_a = match[-4];
                  int at_b = match[-3];
                  int bd_a = match[-2];
                  int bd_b = match[-1];

                  if ((bd_a >= 0) && (bd_a >= 0)) {
                    int a1 = bi_a[bd_a].index[0];
                    int a2 = bi_a[bd_a].index[1];

                    int a1_entry = ai_a[a1].selEntry;
                    int a2_entry = ai_a[a2].selEntry;

                    if ((SelectorIsMember(G, a1_entry, sele1) &&
                         SelectorIsMember(G, a2_entry, sele2)) ||
                        (SelectorIsMember(G, a2_entry, sele1) &&
                         SelectorIsMember(G, a1_entry, sele2))) {
                      /* only update bonds which actually sit between the two
                       * selections */
                      bi_a[bd_a].order = bi_b[bd_b].order;
                      ai_a[at_a].chemFlag = false;
                    }
                  }
                  mi.atom_mark_b[at_b] =
                      false; /* release source for future matching */
                  if (bd_b >= 0) mi.bond_mark_b[bd_b] = false;
                  match -= 4;
                }
                break;
              }
            }
          }
        }
      }
    }
  }
  FreeP(matched);
  return match_found;
}
int ObjectMoleculeVerifyChemistry(ObjectMolecule* I, int state) {
  int result = false;
  const AtomInfoType* ai;
  int a;
  int flag;

  if (state < 0) {
    /* use the first defined state */
    for (a = 0; a < I->NCSet; a++) {
      if (I->CSet[a]) {
        state = a;
        break;
      }
    }
  }
  ai = I->AtomInfo;
  flag = true;
  for (a = 0; a < I->NAtom; a++) {
    if (!ai->chemFlag) {
      flag = false;
    }
    ai++;
  }
  if ((!flag) && (state >= 0) && (state < I->NCSet)) {
    if (I->CSet[state]) {
      ObjectMoleculeInferChemFromBonds(I, state);
      ObjectMoleculeInferChemFromNeighGeom(I, state);
      ObjectMoleculeInferHBondFromChem(I);
      /*      ObjectMoleculeInferChemForProtein(I,0); */
    }
    flag = true;
    ai = I->AtomInfo;
    for (a = 0; a < I->NAtom; a++) {
      if (!ai->chemFlag) {
        flag = false;
        break;
      }
      ai++;
    }
  }
  if (flag) result = true;
  return (result);
}
/**
 * Determines hybridization from coordinates in those few cases where it is
 * unambiguous.
 */
int ObjectMoleculeGetAtomGeometry(const ObjectMolecule* I, int state, int at) {
  int result = -1;
  float v0[3], v1[3], v2[3], v3[3];
  float d1[3], d2[3], d3[3];
  float cp1[3], cp2[3], cp3[3];
  float avg;
  float dp;
  auto const neighbors = AtomNeighbors(I, at);
  auto const nn = neighbors.size();
  if (nn == 4)
    result = cAtomInfoTetrahedral;
  else if (nn == 3) {
    /* check cross products */
    ObjectMoleculeGetAtomVertex(I, state, at, v0);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[0].atm, v1);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[1].atm, v2);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[2].atm, v3);
    subtract3f(v1, v0, d1);
    subtract3f(v2, v0, d2);
    subtract3f(v3, v0, d3);
    cross_product3f(d1, d2, cp1);
    cross_product3f(d2, d3, cp2);
    cross_product3f(d3, d1, cp3);
    normalize3f(cp1);
    normalize3f(cp2);
    normalize3f(cp3);
    avg = (dot_product3f(cp1, cp2) + dot_product3f(cp2, cp3) +
           dot_product3f(cp3, cp1)) /
          3.0F;
    if (avg > 0.75)
      result = cAtomInfoPlanar;
    else
      result = cAtomInfoTetrahedral;
  } else if (nn == 2) {
    ObjectMoleculeGetAtomVertex(I, state, at, v0);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[0].atm, v1);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[1].atm, v2);
    subtract3f(v1, v0, d1);
    subtract3f(v2, v0, d2);
    normalize3f(d1);
    normalize3f(d2);
    dp = dot_product3f(d1, d2);
    if (dp < -0.75) result = cAtomInfoLinear;
  }
  return (result);
}
static float compute_avg_center_dot_cross_fn(ObjectMolecule* I, CoordSet* cs,
                                             int n_atom, int* atix) {
  float result = 0.0F;
  float* v[9];
  int missing_flag = false;
  {
    int a, i;
    for (i = 0; i < n_atom; i++) {
      int a1 = atix[i];
      a = cs->atmToIdx(a1);
      if (a < 0) {
        missing_flag = true;
        break;
      } else {
        v[i] = cs->coordPtr(a);
      }
    }
  }
  if (!missing_flag) {
    int i;
    float d10[3], d20[3], cp[8][3];
    float avg = 0.0;

    v[n_atom] = v[1];
    for (i = 1; i < n_atom; i++) {
      subtract3f(v[i], v[0], d10);
      subtract3f(v[i + 1], v[0], d20);
      normalize3f(d10);
      normalize3f(d20);
      cross_product3f(d10, d20, cp[i]);
      normalize3f(cp[i]);
      if (i > 1) {
        if (dot_product3f(cp[i - 1], cp[i]) < 0.0) invert3f(cp[i]);
      }
    }
    copy3f(cp[1], cp[n_atom]);
    for (i = 1; i < n_atom; i++) {
      avg += dot_product3f(cp[i], cp[i + 1]);
    }
    result = avg / (n_atom - 1);
  }
  return result;
}

static int verify_planer_bonds(const ObjectMolecule* I, const CoordSet* cs,
                               int n_atom, const int* atix, const int* neighbor,
                               const float* dir, float cutoff) {
  int a, i;
  for (i = 0; i < n_atom; i++) {
    int a1 = atix[i];
    a = cs->atmToIdx(a1);
    if (a >= 0) {
      const float* v0 = cs->coordPtr(a);
      int n = neighbor[a1] + 1;
      int a2;
      while ((a2 = neighbor[n]) >= 0) {
        n += 2;
        a = cs->atmToIdx(a2);
        if (a >= 0) {
          const float* v1 = cs->coordPtr(a);
          float d10[] = {0.f, 0.f, 0.f};
          subtract3f(v1, v0, d10);
          normalize3f(d10);
          {
            float dot = fabs(dot_product3f(d10, dir));
            if (dot > cutoff) {
              switch (I->AtomInfo[a1].protons) {
                case cAN_C:
                case cAN_N:
                case cAN_O:
                case cAN_S:
                  switch (I->AtomInfo[a2].protons) {
                    case cAN_C:
                    case cAN_N:
                    case cAN_O:
                    case cAN_S:
                      return 0;
                      break;
                  }
                  break;
              }
            }
          }
        }
      }
    }
  }
  return 1;
}

static float compute_avg_ring_dot_cross_fn(ObjectMolecule* I, CoordSet* cs,
                                           int n_atom, int* atix, float* dir) {
  float result = 0.0F;
  float* v[9];
  int missing_flag = false;
  {
    int a, i;
    for (i = 0; i < n_atom; i++) {
      int a1 = atix[i];
      a = cs->atmToIdx(a1);
      if (a < 0) {
        missing_flag = true;
        break;
      } else {
        v[i] = cs->coordPtr(a);
      }
    }
  }
  if (!missing_flag) {
    int i;
    float d10[3], d21[3], cp[8][3];
    float avg = 0.0;

    v[n_atom] = v[0];
    v[n_atom + 1] = v[1];
    for (i = 0; i < n_atom; i++) {
      subtract3f(v[i + 0], v[i + 1], d10);
      subtract3f(v[i + 2], v[i + 1], d21);
      normalize3f(d10);
      normalize3f(d21);
      cross_product3f(d10, d21, cp[i]);
      normalize3f(cp[i]);
      if (i) {
        if (dot_product3f(cp[i - 1], cp[i]) < 0.0) invert3f(cp[i]);
      }
      add3f(cp[i], dir, dir);
    }
    copy3f(cp[0], cp[n_atom]);
    for (i = 0; i < n_atom; i++) {
      avg += dot_product3f(cp[i], cp[i + 1]);
    }
    result = avg / n_atom;
  }
  normalize3f(dir);
  return result;
}

typedef struct {
  int cyclic, planer, aromatic;
} ObservedInfo;

void ObjectMoleculeGuessValences(ObjectMolecule* I, int state, int* flag1,
                                 int* flag2, int reset) {
  /* this a hacked 80% solution ...it will get things wrong, but it is
     better than nothing! */

  const float planer_cutoff = 0.96F;
  CoordSet* cs = nullptr;
  ObservedInfo* obs_atom = nullptr;
  ObservedInfo* obs_bond = nullptr;
  int* flag = nullptr;
  int warning1 = 0, warning2 = 0;

  /* WORKAROUND of a possible -funroll-loops inlining optimizer bug in gcc 3.3.3
   */

  float (*compute_avg_center_dot_cross)(ObjectMolecule*, CoordSet*, int, int*);
  float (*compute_avg_ring_dot_cross)(ObjectMolecule*, CoordSet*, int, int*,
                                      float*);

  compute_avg_center_dot_cross = compute_avg_center_dot_cross_fn;
  compute_avg_ring_dot_cross = compute_avg_ring_dot_cross_fn;

  /* end WORKAROUND */

  if ((state >= 0) && (state < I->NCSet)) {
    cs = I->CSet[state];
  }
  if (cs) {
    obs_atom = pymol::calloc<ObservedInfo>(I->NAtom);
    obs_bond = pymol::calloc<ObservedInfo>(I->NBond);
  }
  flag = pymol::calloc<int>(I->NAtom);
  if (flag) {
    if (!flag1) {
      int a, *flag_a = flag;
      const AtomInfoType* ai = I->AtomInfo.data();
      /* default behavior: only reset hetatm valences */
      for (a = 0; a < I->NAtom; a++) {
        *(flag_a++) = (ai++)->hetatm;
      }
    } else if (flag1 && flag2) {
      int a, *flag_a = flag, *flag1_a = flag1, *flag2_a = flag2;
      for (a = 0; a < I->NAtom; a++) {
        *(flag_a++) = (*(flag1_a++) || *(flag2_a++));
      }
    } else if (flag1) {
      int a, *flag_a = flag, *flag1_a = flag1;
      for (a = 0; a < I->NAtom; a++) {
        *(flag_a++) = *(flag1_a++);
      }
    }
  }
  if (!flag1) flag1 = flag;
  if (!flag2) flag2 = flag;
  if (reset) {
    /* reset chemistry information and bond orders for selected atoms */
    {
      int a;
      AtomInfoType* ai = I->AtomInfo.data();
      for (a = 0; a < I->NAtom; a++) {
        if (flag[a]) ai->chemFlag = 0;
        ai++;
      }
    }
    {
      int b;
      BondType* bi = I->Bond.data();
      for (b = 0; b < I->NBond; b++) {
        int at0 = bi->index[0];
        int at1 = bi->index[1];
        if ((flag1[at0] && flag2[at1]) || (flag1[at1] && flag2[at0]))
          bi->order = 1;
        bi++;
      }
    }
  }
  if (cs && obs_bond && obs_atom && flag && flag1 && flag2) {
    int a;
    int const* const neighbor = I->getNeighborArray();
    AtomInfoType* atomInfo = I->AtomInfo.data();
    BondType* bondInfo = I->Bond.data();

    for (a = 0; a < I->NAtom; a++) {
      AtomInfoType* ai = atomInfo + a;
      if ((!ai->chemFlag) && (flag[a])) {
        /*  determine whether or not atom participates in a planer system with 5
         * or 6 atoms */

        {
          int mem[9];
          int nbr[7];
          const int ESCAPE_MAX = 500;

          const int* atmToIdx = I->DiscreteFlag ? nullptr : cs->AtmToIdx.data();

          int escape_count =
              ESCAPE_MAX; /* don't get bogged down with structures
                         that have unreasonable connectivity */
          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while (((mem[1] = neighbor[nbr[0]]) >= 0) &&
                 ((!atmToIdx) || (atmToIdx[mem[0]] >= 0))) {
            nbr[1] = neighbor[mem[1]] + 1;
            while (((mem[2] = neighbor[nbr[1]]) >= 0) &&
                   ((!atmToIdx) || (atmToIdx[mem[1]] >= 0))) {
              if (mem[2] != mem[0]) {
                nbr[2] = neighbor[mem[2]] + 1;
                while (((mem[3] = neighbor[nbr[2]]) >= 0) &&
                       ((!atmToIdx) || (atmToIdx[mem[2]] >= 0))) {
                  if (mem[3] != mem[1]) {
                    nbr[3] = neighbor[mem[3]] + 1;
                    while (((mem[4] = neighbor[nbr[3]]) >= 0) &&
                           ((!atmToIdx) || (atmToIdx[mem[3]] >= 0))) {
                      if ((mem[4] != mem[2]) && (mem[4] != mem[1]) &&
                          (mem[4] != mem[0])) {
                        nbr[4] = neighbor[mem[4]] + 1;
                        while (((mem[5] = neighbor[nbr[4]]) >= 0) &&
                               ((!atmToIdx) || (atmToIdx[mem[4]] >= 0))) {
                          if (!(escape_count--)) goto escape;
                          if ((mem[5] != mem[3]) && (mem[5] != mem[2]) &&
                              (mem[5] != mem[1])) {
                            if (mem[5] == mem[0]) { /* five-cycle */
                              int i;
                              float dir[] = {0.f, 0.f, 0.f};
                              float avg_dot_cross = compute_avg_ring_dot_cross(
                                  I, cs, 5, mem, dir);
                              for (i = 0; i < 5; i++) {
                                obs_atom[mem[i]].cyclic = true;
                                obs_bond[neighbor[nbr[0] + 1]].cyclic = true;
                              }
                              if (avg_dot_cross > planer_cutoff) {
                                if (verify_planer_bonds(I, cs, 5, mem, neighbor,
                                                        dir, 0.35F)) {
                                  for (i = 0; i < 5; i++) {
                                    obs_atom[mem[i]].planer = true;
                                    obs_bond[neighbor[nbr[i] + 1]].planer =
                                        true;
                                  }
                                }
                              }
                            }

                            nbr[5] = neighbor[mem[5]] + 1;
                            while (((mem[6] = neighbor[nbr[5]]) >= 0) &&
                                   ((!atmToIdx) || (atmToIdx[mem[5]] >= 0))) {
                              if ((mem[6] != mem[4]) && (mem[6] != mem[3]) &&
                                  (mem[6] != mem[2]) && (mem[6] != mem[1])) {
                                if (mem[6] == mem[0]) { /* six-cycle */
                                  int i;
                                  float dir[3] = {0.f, 0.f, 0.f};
                                  float avg_dot_cross =
                                      compute_avg_ring_dot_cross(I, cs, 6, mem,
                                                                 dir);
                                  for (i = 0; i < 6; i++) {
                                    obs_atom[mem[i]].cyclic = true;
                                    obs_bond[neighbor[nbr[i] + 1]].cyclic =
                                        true;
                                  }
                                  if (avg_dot_cross > planer_cutoff) {
                                    if (verify_planer_bonds(I, cs, 6, mem,
                                                            neighbor, dir,
                                                            0.35F)) {
                                      for (i = 0; i < 6; i++) {
                                        obs_atom[mem[i]].planer = true;
                                        obs_bond[neighbor[nbr[i] + 1]].planer =
                                            true;
                                      }
                                    }
                                  }
                                }
                              }
                              nbr[5] += 2;
                            }
                          }
                          nbr[4] += 2;
                        }
                      }
                      nbr[3] += 2;
                    }
                  }
                  nbr[2] += 2;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }

        escape:
          escape_count = ESCAPE_MAX; /* don't get bogged down with structures
                                        that have unreasonable connectivity */
          warning1 = 1;
        }
      }
    }

    for (a = 0; a < I->NAtom; a++) {
      AtomInfoType* ai = atomInfo + a;
      if ((!ai->chemFlag) && flag[a]) {
        ObservedInfo* ob_at = obs_atom + a;

        {
          switch (ai->protons) {
            case cAN_P:
            case cAN_S:
            case cAN_N:
            case cAN_C: {
              int atm[5];
              int bnd[5];
              int n = neighbor[a];
              int nn = neighbor[n++];
              if (nn > 4) nn = 4;
              atm[0] = a;
              {
                int i;
                for (i = 1; i <= nn; i++) {
                  atm[i] = neighbor[n];
                  bnd[i] = neighbor[n + 1];
                  n += 2;
                }
              }
              {
                int o1_at = -1, o2_at = -1, o3_at = -1, o4_at = -1;
                int o1_bd = 0, o2_bd = 0, o3_bd = 0, o4_bd = 0;
                float o1_len = 0.0F, o2_len = 0.0F;
                float n1_v[3] = {0.0F, 0.0F, 0.0F};
                int n1_at = -1, n2_at = -1, n3_at = -1;
                int n1_bd = 0, n2_bd = 0, n3_bd = 0;
                float n1_len = 0.0F, n2_len = 0.0F, n3_len = 0.0F;
                int c1_at = -1, c2_at = -1;
                float c1_v[3] = {0.0F, 0.0F, 0.0F};
                float* v0 = nullptr;

                {
                  int idx0 = cs->atmToIdx(a);
                  if (idx0 >= 0) v0 = cs->coordPtr(idx0);
                }
                {
                  int i;
                  for (i = 1; i <= nn; i++) {
                    float* v1 = nullptr;

                    {
                      int idx1 = cs->atmToIdx(atm[i]);
                      if (idx1 >= 0) v1 = cs->coordPtr(idx1);
                    }
                    if (v0 && v1) {
                      float diff[3];
                      subtract3f(v1, v0, diff);
                      {
                        float bond_len = length3f(diff);

                        switch (atomInfo[atm[i]].protons) {
                          case cAN_C:
                            if (c1_at < 0) {
                              c1_at = atm[i];
                              copy3f(diff, c1_v);
                            } else if (c2_at < 0) {
                              c2_at = atm[i];
                            }
                            break;
                          case cAN_O:
                            if (o1_at < 0) {
                              o1_at = atm[i];
                              o1_len = bond_len;
                              o1_bd = bnd[i];
                            } else if (o2_at < 0) {
                              o2_at = atm[i];
                              o2_len = bond_len;
                              o2_bd = bnd[i];
                            } else if (o3_at < 0) {
                              o3_at = atm[i];
                              o3_bd = bnd[i];
                            } else if (o4_at < 0) {
                              o4_at = atm[i];
                              o4_bd = bnd[i];
                            }
                            break;
                          case cAN_N:
                            if (n1_at < 0) {
                              copy3f(diff, n1_v);
                              n1_at = atm[i];
                              n1_len = bond_len;
                              n1_bd = bnd[i];
                            } else if (n2_at < 0) {
                              n2_at = atm[i];
                              n2_len = bond_len;
                              n2_bd = bnd[i];
                            } else if (n3_at < 0) {
                              n3_at = atm[i];
                              n3_len = bond_len;
                              n3_bd = bnd[i];
                            }
                            break;
                        }
                      }
                    }
                  }
                }
                {
                  float avg_dot_cross = 0.0F;

                  switch (ai->protons) {
                    case cAN_C: /* planer carbons */
                      if (nn == 3) {
                        avg_dot_cross =
                            compute_avg_center_dot_cross(I, cs, 4, atm);

                        if (avg_dot_cross > planer_cutoff) {
                          if ((n1_at >= 0) && (o1_at >= 0) && (o2_at < 0)) {
                            /* simple amide? */
                            if (o1_len < 1.38F) {
                              if (neighbor[neighbor[o1_at]] == 1)
                                if ((flag1[a] && flag2[o1_at]) ||
                                    (flag2[a] && flag1[o1_at]))
                                  bondInfo[o1_bd].order = 2;
                            }
                          } else if ((n1_at >= 0) && (o1_at >= 0) &&
                                     (o2_at >= 0)) {
                            /* carbamyl */
                            if ((o1_len < 1.38F) &&
                                (neighbor[neighbor[o1_at]] == 1) &&
                                (o2_len < 1.38F) &&
                                (neighbor[neighbor[o2_at]] == 1)) {
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 4;
                              if ((flag1[a] && flag2[o2_at]) ||
                                  (flag2[a] && flag1[o2_at]))
                                bondInfo[o2_bd].order = 4;
                            } else if ((o1_len < 1.38F) &&
                                       (neighbor[neighbor[o1_at]] == 1)) {
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 2;
                            } else if ((o2_len < 1.38F) &&
                                       (neighbor[neighbor[o2_at]] == 1))
                              if ((flag1[a] && flag2[o2_at]) ||
                                  (flag2[a] && flag1[o2_at]))
                                bondInfo[o2_bd].order = 2;
                          } else if ((n1_at < 0) && (o1_at >= 0) &&
                                     (o2_at < 0)) {
                            /* ketone */
                            if ((o1_len < 1.31F) &&
                                (neighbor[neighbor[o1_at]] == 1)) {
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 2;
                            }
                          } else if ((o1_at >= 0) && (o2_at >= 0) &&
                                     (n1_at < 0)) {
                            /* simple carboxylate? */
                            if ((o1_len < 1.38F) && (o2_len < 1.38F) &&
                                (neighbor[neighbor[o1_at]] == 1) &&
                                (neighbor[neighbor[o2_at]] == 1)) {
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 4;
                              if ((flag1[a] && flag2[o2_at]) ||
                                  (flag2[a] && flag1[o2_at]))
                                bondInfo[o2_bd].order = 4;
                            } else if ((o1_len < 1.38F) &&
                                       (neighbor[neighbor[o1_at]] ==
                                        1)) { /* esters */
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 2;
                            } else if ((o2_len < 1.38F) &&
                                       (neighbor[neighbor[o2_at]] == 1)) {
                              if ((flag1[a] && flag2[o2_at]) ||
                                  (flag2[a] && flag1[o2_at]))
                                bondInfo[o2_bd].order = 2;
                            }
                          } else if ((n1_at >= 0) && (n2_at >= 0) &&
                                     (n3_at < 0) && (c1_at >= 0) &&
                                     (n1_len < 1.43F) && (n2_len < 1.43F) &&
                                     obs_atom[c1_at].planer &&
                                     obs_atom[c1_at].cyclic &&
                                     (!ob_at->cyclic) &&
                                     (!obs_atom[n1_at].cyclic) &&
                                     (!obs_atom[n2_at].cyclic)) {
                            if ((flag1[a] && flag2[n1_at]) ||
                                (flag2[a] && flag1[n1_at]))
                              bondInfo[n1_bd].order = 4;
                            atomInfo[n1_at].valence = 3;
                            atomInfo[n1_at].geom = cAtomInfoPlanar;
                            atomInfo[n1_at].chemFlag = 2;
                            if ((flag1[a] && flag2[n2_at]) ||
                                (flag2[a] && flag1[n2_at]))
                              bondInfo[n2_bd].order = 4;
                            atomInfo[n2_at].valence = 3;
                            atomInfo[n2_at].geom = cAtomInfoPlanar;
                            atomInfo[n2_at].chemFlag = 2;
                          } else if ((n1_at >= 0) && (n2_at >= 0) &&
                                     (n3_at >= 0)) {
                            /* guanido with no hydrogens */
                            if ((n1_len < 1.44F) && (n2_len < 1.44F) &&
                                (n3_len < 1.44F)) {
                              if ((neighbor[neighbor[n1_at]] == 1) &&
                                  (neighbor[neighbor[n2_at]] == 1) &&
                                  (neighbor[neighbor[n3_at]] >= 2)) {
                                if ((flag1[a] && flag2[n1_at]) ||
                                    (flag2[a] && flag1[n1_at]))
                                  bondInfo[n1_bd].order = 4;
                                atomInfo[n1_at].valence = 3;
                                atomInfo[n1_at].geom = cAtomInfoPlanar;
                                atomInfo[n1_at].chemFlag = 2;
                                if ((flag1[a] && flag2[n2_at]) ||
                                    (flag2[a] && flag1[n2_at]))
                                  bondInfo[n2_bd].order = 4;
                                atomInfo[n2_at].valence = 3;
                                atomInfo[n2_at].geom = cAtomInfoPlanar;
                                atomInfo[n2_at].chemFlag = 2;
                              } else if ((neighbor[neighbor[n1_at]] == 1) &&
                                         (neighbor[neighbor[n2_at]] >= 2) &&
                                         (neighbor[neighbor[n3_at]] == 1)) {
                                if ((flag1[a] && flag2[n1_at]) ||
                                    (flag2[a] && flag1[n1_at]))
                                  bondInfo[n1_bd].order = 4;
                                atomInfo[n1_at].valence = 3;
                                atomInfo[n1_at].geom = cAtomInfoPlanar;
                                atomInfo[n1_at].chemFlag = 2;
                                if ((flag1[a] && flag2[n3_at]) ||
                                    (flag2[a] && flag1[n3_at]))
                                  bondInfo[n3_bd].order = 4;
                                atomInfo[n3_at].valence = 3;
                                atomInfo[n3_at].geom = cAtomInfoPlanar;
                                atomInfo[n3_at].chemFlag = 2;
                              } else if ((neighbor[neighbor[n1_at]] >= 2) &&
                                         (neighbor[neighbor[n2_at]] == 1) &&
                                         (neighbor[neighbor[n3_at]] == 1)) {
                                if ((flag1[a] && flag2[n2_at]) ||
                                    (flag2[a] && flag1[n2_at]))
                                  bondInfo[n2_bd].order = 4;
                                atomInfo[n2_at].valence = 3;
                                atomInfo[n2_at].geom = cAtomInfoPlanar;
                                atomInfo[n2_at].chemFlag = 2;
                                if ((flag1[a] && flag2[n3_at]) ||
                                    (flag2[a] && flag1[n3_at]))
                                  bondInfo[n3_bd].order = 4;
                                atomInfo[n3_at].valence = 3;
                                atomInfo[n3_at].geom = cAtomInfoPlanar;
                                atomInfo[n3_at].chemFlag = 2;
                              }
                            }
                          }
                        }
                      }
                      /* any carbon */

                      /* handle imines and nitriles */

                      if ((nn >= 2) && (nn <= 3) && (n1_at >= 0) &&
                          (o1_at < 0) && (n2_at < 0) && (n1_len < 1.36F) &&
                          (!ob_at->cyclic) && (!obs_bond[n1_bd].cyclic) &&
                          (!obs_atom[n1_at].planer) &&
                          ((nn == 2) ||
                           ((nn == 3) && (avg_dot_cross > planer_cutoff)))) {
                        float n1_dot_cross = 1.0F;
                        int n2 = neighbor[n1_at];
                        int nn2 = neighbor[n2++];

                        { /* check nitrogen planarity */
                          int atm2[5];
                          if (nn2 > 2) {
                            nn2 = 3;
                            atm2[0] = n1_at;
                            {
                              int i2;
                              for (i2 = 1; i2 <= nn2; i2++) {
                                atm2[i2] = neighbor[n2];
                                n2 += 2;
                              }
                            }
                            n1_dot_cross =
                                compute_avg_center_dot_cross(I, cs, 4, atm2);
                          }
                        }
                        if (n1_dot_cross > planer_cutoff) {
                          if ((flag1[a] && flag2[n1_at]) ||
                              (flag2[a] && flag1[n1_at]))
                            bondInfo[n1_bd].order = 2;
                          if ((n1_len < 1.24F) && (c1_at >= 0) && (nn2 == 1)) {
                            normalize3f(n1_v);
                            normalize3f(c1_v);
                            if (dot_product3f(n1_v, c1_v) < -0.9) {
                              if ((flag1[a] && flag2[n1_at]) ||
                                  (flag2[a] && flag1[n1_at]))
                                bondInfo[n1_bd].order = 3;
                            }
                          }
                        }
                      }
                      break;
                    case cAN_N:
                      if (nn == 3) {
                        avg_dot_cross =
                            compute_avg_center_dot_cross(I, cs, 4, atm);

                        if ((avg_dot_cross > planer_cutoff)) {
                          if ((o1_at >= 0) && (o2_at >= 0) && (o3_at < 0)) {
                            /* nitro */
                            if (neighbor[neighbor[o1_at]] == 1) {
                              if ((flag1[a] && flag2[o1_at]) ||
                                  (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 4;
                            }
                            if (neighbor[neighbor[o2_at]] == 1) {
                              if ((flag1[a] && flag2[o2_at]) ||
                                  (flag2[a] && flag1[o2_at]))
                                bondInfo[o2_bd].order = 4;
                            }
                          }
                        }
                      }
                      break;
                    case cAN_S:
                    case cAN_P:
                      if ((o1_at >= 0) && (o2_at >= 0) && (o3_at >= 0) &&
                          (o4_at >= 0)) {
                        /* sulfate, phosphate */
                        int o1 = -1, o2 = -1, o3 = -1;
                        int a1 = 0;
                        if (neighbor[neighbor[o1_at]] == 1) {
                          o1 = o1_bd;
                          a1 = o1_at;
                        }
                        if (neighbor[neighbor[o2_at]] == 1) {
                          if (o1 < 0) {
                            o1 = o2_bd;
                            a1 = o2_at;
                          } else if (o2 < 0) {
                            o2 = o2_bd;
                          }
                        }
                        if (neighbor[neighbor[o3_at]] == 1) {
                          if (o1 < 0) {
                            o1 = o3_bd;
                            a1 = o3_at;
                          } else if (o2 < 0)
                            o2 = o3_bd;
                          else if (o3 < 0)
                            o3 = o3_bd;
                        }
                        if (neighbor[neighbor[o4_at]] == 1) {
                          if (o1 < 0) {
                            o1 = o4_bd;
                            a1 = o4_at;
                          } else if (o2 < 0)
                            o2 = o4_bd;
                          else if (o3 < 0)
                            o3 = o4_bd;
                        }
                        if (o2 >= 0) {
                          if ((flag1[a] && flag2[a1]) ||
                              (flag2[a] && flag1[a1]))
                            bondInfo[o1].order = 2;
                          if (o2 == o2_bd) {
                            atomInfo[o2_at].formalCharge = -1;
                          } else if (o2 == o3_bd) {
                            atomInfo[o3_at].formalCharge = -1;
                          } else if (o2 == o4_bd) {
                            atomInfo[o4_at].formalCharge = -1;
                          }
                        }
                      } else if ((o1_at >= 0) && (o2_at >= 0) && (o3_at >= 0) &&
                                 (o4_at < 0)) {
                        /* sulfonamide */
                        int o1 = -1, o2 = -1;
                        int a1 = 0, a2 = 0;
                        if (neighbor[neighbor[o1_at]] == 1) {
                          o1 = o1_bd;
                          a1 = o1_at;
                        }
                        if (neighbor[neighbor[o2_at]] == 1) {
                          if (o1 < 0) {
                            o1 = o2_bd;
                            a1 = o2_at;
                          } else if (o2 < 0) {
                            o2 = o2_bd;
                            a2 = o2_at;
                          }
                        }
                        if (neighbor[neighbor[o3_at]] == 1) {
                          if (o1 < 0) {
                            o1 = o3_bd;
                            a1 = o3_at;
                          } else if (o2 < 0) {
                            o2 = o3_bd;
                            a2 = o3_at;
                          }
                        }
                        if (o1 >= 0) {
                          if ((flag1[a] && flag2[a1]) ||
                              (flag2[a] && flag1[a1]))
                            bondInfo[o1].order = 2;
                        }
                        if (o2 >= 0) {
                          if ((flag1[a] && flag2[a2]) ||
                              (flag2[a] && flag1[a2]))
                            bondInfo[o2].order = 2;
                        }
                      } else if ((o1_at >= 0) && (o2_at >= 0) && (o3_at < 0)) {
                        /* sulphone */
                        if (neighbor[neighbor[o1_at]] == 1) {
                          if ((flag1[a] && flag2[o1_at]) ||
                              (flag2[a] && flag1[o1_at]))
                            bondInfo[o1_bd].order = 2;
                        }
                        if (neighbor[neighbor[o2_at]] == 1) {
                          if ((flag1[a] && flag2[o2_at]) ||
                              (flag2[a] && flag1[o2_at]))
                            bondInfo[o2_bd].order = 2;
                        }
                      }
                      break;
                  }
                }
              }
            }
          }
          /* this sets aromatic bonds for cyclic planer systems */

          if (ob_at->cyclic && ob_at->planer) {
            switch (ai->protons) {
              case cAN_C:
              case cAN_N:
              case cAN_O:
              case cAN_S: {
                int n, a0, b0;
                n = neighbor[a] + 1;
                while (1) {
                  a0 = neighbor[n];
                  if (a0 < 0) break;
                  b0 = neighbor[n + 1];
                  n += 2;
                  if (obs_atom[a0].cyclic && obs_atom[a0].planer &&
                      obs_bond[b0].cyclic) {
                    obs_atom[a0].aromatic = true;
                    switch (I->AtomInfo[a0].protons) {
                      case cAN_C:
                      case cAN_N:
                      case cAN_O:
                      case cAN_S:
                        if ((flag1[a] && flag2[a0]) || (flag2[a] && flag1[a0]))
                          I->Bond[b0].order = 4;
                        break;
                    }
                  }
                }
              } break;
            }
          }
        }
      }
    }

    /* now try to address some simple cases with aromatic nitrogens */
    for (a = 0; a < I->NAtom; a++) {
      AtomInfoType* ai = atomInfo + a;
      if ((!ai->chemFlag) && (ai->protons == cAN_N) &&
          (ai->formalCharge == 0) && flag[a] && obs_atom[a].cyclic &&
          obs_atom[a].aromatic) {
        int n = neighbor[a];
        int nn = neighbor[n++];

        if (nn == 2) { /* only two explicit neighbors */

          int mem[9];
          int nbr[7];
          const int ESCAPE_MAX = 500;

          const int* atmToIdx = I->DiscreteFlag ? nullptr : cs->AtmToIdx.data();

          int escape_count =
              ESCAPE_MAX; /* don't get bogged down with structures
                         that have unreasonable connectivity */
          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while (((mem[1] = neighbor[nbr[0]]) >= 0) &&
                 ((!atmToIdx) || (atmToIdx[mem[0]] >= 0))) {
            nbr[1] = neighbor[mem[1]] + 1;
            while (((mem[2] = neighbor[nbr[1]]) >= 0) &&
                   ((!atmToIdx) || (atmToIdx[mem[1]] >= 0))) {
              if (mem[2] != mem[0]) {
                nbr[2] = neighbor[mem[2]] + 1;
                while (((mem[3] = neighbor[nbr[2]]) >= 0) &&
                       ((!atmToIdx) || (atmToIdx[mem[2]] >= 0))) {
                  if (mem[3] != mem[1]) {
                    nbr[3] = neighbor[mem[3]] + 1;
                    while (((mem[4] = neighbor[nbr[3]]) >= 0) &&
                           ((!atmToIdx) || (atmToIdx[mem[3]] >= 0))) {
                      if ((mem[4] != mem[2]) && (mem[4] != mem[1]) &&
                          (mem[4] != mem[0])) {
                        nbr[4] = neighbor[mem[4]] + 1;
                        while (((mem[5] = neighbor[nbr[4]]) >= 0) &&
                               ((!atmToIdx) || (atmToIdx[mem[4]] >= 0))) {
                          if (!(escape_count--))
                            goto escape2; /* BUG FIX: need a new escape2,
                                             instead of mistakenly going back to
                                             the first escape, which is in the
                                             loop above */
                          if ((mem[5] != mem[3]) && (mem[5] != mem[2]) &&
                              (mem[5] != mem[1])) {
                            if (mem[5] == mem[0] && (!ai->chemFlag)) {
                              /* unassigned aromatic nitrogen-containing
                               * five-cycle */

                              /* c1ccnc1 becomes c1ccn[H]c1 */

                              if ((atomInfo[mem[1]].protons == cAN_C) &&
                                  (atomInfo[mem[2]].protons == cAN_C) &&
                                  (atomInfo[mem[3]].protons == cAN_C) &&
                                  (atomInfo[mem[4]].protons == cAN_C) &&
                                  obs_atom[mem[1]].aromatic &&
                                  obs_atom[mem[2]].aromatic &&
                                  obs_atom[mem[3]].aromatic &&
                                  obs_atom[mem[4]].aromatic) {
                                ai->valence = 3;
                                ai->chemFlag = 2;
                                ai->geom = cAtomInfoPlanar;
                              }

                              /* c1ncnc1 becomes c1n[H]cnc1 */

                              if ((atomInfo[mem[1]].protons == cAN_C) &&
                                  (atomInfo[mem[2]].protons == cAN_N) &&
                                  (atomInfo[mem[2]].formalCharge == 0) &&
                                  (!atomInfo[mem[2]].chemFlag) &&
                                  (atomInfo[mem[3]].protons == cAN_C) &&
                                  (atomInfo[mem[4]].protons == cAN_C) &&
                                  obs_atom[mem[1]].aromatic &&
                                  obs_atom[mem[2]].aromatic &&
                                  obs_atom[mem[3]].aromatic &&
                                  obs_atom[mem[4]].aromatic) {
                                int n2 = neighbor[mem[2]];
                                int nn2 = neighbor[n2++];
                                if (nn2 ==
                                    2) { /* second nitrogen also ambiguous */
                                  ai->valence = 3;
                                  ai->chemFlag = 2;
                                  ai->geom = cAtomInfoPlanar;
                                }
                              }

                              /* c1cnnc1 becomes c1cn[H]nc1 */

                              if ((atomInfo[mem[1]].protons == cAN_N) &&
                                  (atomInfo[mem[1]].formalCharge == 0) &&
                                  (!atomInfo[mem[1]].chemFlag) &&
                                  (atomInfo[mem[2]].protons == cAN_C) &&
                                  (atomInfo[mem[3]].protons == cAN_C) &&
                                  (atomInfo[mem[4]].protons == cAN_C) &&
                                  obs_atom[mem[1]].aromatic &&
                                  obs_atom[mem[2]].aromatic &&
                                  obs_atom[mem[3]].aromatic &&
                                  obs_atom[mem[4]].aromatic) {
                                int n2 = neighbor[mem[1]];
                                int nn2 = neighbor[n2++];
                                if (nn2 ==
                                    2) { /* second nitrogen also ambiguous */
                                  ai->valence = 3;
                                  ai->chemFlag = 2;
                                  ai->geom = cAtomInfoPlanar;
                                }
                              }
                            }
                          }
                          nbr[4] += 2;
                        }
                      }
                      nbr[3] += 2;
                    }
                  }
                  nbr[2] += 2;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }
        escape2: /* BUG FIX: Need separate escape for this loop */
          escape_count = ESCAPE_MAX; /* don't get bogged down with structures
                                        that have unreasonable connectivity */
          warning2 = 1;
        }
      }
    }
  }
  if (warning1 || warning2) {
    PRINTFB(I->G, FB_ObjectMolecule, FB_Blather)
    " %s(%d,%d): Unreasonable connectivity in heteroatom,\n  unsuccessful in "
    "guessing valences.\n",
        __func__, warning1, warning2 ENDFB(I->G);
  }
  FreeP(obs_bond);
  FreeP(obs_atom);
  FreeP(flag);
}
void ObjectMoleculeInferChemForProtein(ObjectMolecule* I, int state) {
  /* Infers chemical relations for a molecules under protein assumptions.
   *
   * NOTE: this routine needs an all-atom model (with hydrogens!)
   * and it will make mistakes on non-protein atoms (if they haven't
   * already been assigned)
   */

  int changedFlag = true;

  /* first, try to find all amids and acids */
  while (changedFlag) {
    changedFlag = false;
    for (int a = 0; a < I->NAtom; a++) {
      auto const* const ai = I->AtomInfo.data() + a;
      if (ai->chemFlag) {
        if (ai->geom == cAtomInfoPlanar)
          if (ai->protons == cAN_C) {
            auto const neighbors = AtomNeighbors(I, a);
            if (neighbors.size() > 1) {
              AtomInfoType* ai1 = nullptr;
              for (auto const& neighbor : neighbors) {
                auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
                if ((ai0->protons == cAN_O) && (!ai0->chemFlag)) {
                  ai1 = ai0; /* found candidate carbonyl */
                  break;
                }
              }
              if (ai1) {
                for (auto const& neighbor : neighbors) {
                  auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
                  if (ai0 != ai1) {
                    if (ai0->protons == cAN_O) {
                      if (!ai0->chemFlag) {
                        ai0->chemFlag = true; /* acid */
                        ai0->geom = cAtomInfoPlanar;
                        ai0->valence = 1;
                        ai1->chemFlag = true;
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
                      }
                    } else if (ai0->protons == cAN_N) {
                      if (!ai0->chemFlag) {
                        ai0->chemFlag = true; /* amide N */
                        ai0->geom = cAtomInfoPlanar;
                        ai0->valence = 3;
                        ai1->chemFlag = true; /* amide =O */
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
                      } else if (ai0->geom == cAtomInfoPlanar) {
                        ai1->chemFlag = true; /* amide =O */
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
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
  /* then handle aldehydes and amines (partial amides - both missing a valence)
   */

  changedFlag = true;
  while (changedFlag) {
    changedFlag = false;
    for (int a = 0; a < I->NAtom; a++) {
      auto* const ai = I->AtomInfo.data() + a;
      if (!ai->chemFlag) {
        if (ai->protons == cAN_C) {
          auto const neighbors = AtomNeighbors(I, a);
          if (neighbors.size() > 1) {
            for (auto const& neighbor : neighbors) {
              auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
              if ((ai0->protons == cAN_O) && (!ai0->chemFlag)) { /* =O */
                ai->chemFlag = true;
                ai->geom = cAtomInfoPlanar;
                ai->valence = 1;
                ai0->chemFlag = true;
                ai0->geom = cAtomInfoPlanar;
                ai0->valence = 3;
                changedFlag = true;
                break;
              }
            }
          }
        } else if (ai->protons == cAN_N) {
          if ((!ai->chemFlag) || ai->geom != cAtomInfoLinear) {
            if (ai->formalCharge == 0) {
              ai->chemFlag = true;
              ai->geom = cAtomInfoPlanar;
              ai->valence = 3;
            }
          }
        }
      }
    }
  }
}
/**
 * Assigns:
 * - geom
 * - valence
 * - chemFlag
 */
void ObjectMoleculeInferChemFromNeighGeom(ObjectMolecule* I, int state) {
  /* infers chemical relations from neighbors and geometry
   * NOTE: very limited in scope */

  int a;
  int changedFlag = true;
  int geom;
  int carbonVal[10];

  AtomInfoType *ai, *ai2;

  carbonVal[cAtomInfoTetrahedral] = 4;
  carbonVal[cAtomInfoPlanar] = 3;
  carbonVal[cAtomInfoLinear] = 2;

  while (changedFlag) {
    changedFlag = false;
    for (a = 0; a < I->NAtom; a++) {
      ai = I->AtomInfo + a;
      if (!ai->chemFlag) {
        geom = ObjectMoleculeGetAtomGeometry(I, state, a);
        switch (ai->protons) {
          case cAN_K:
            ai->chemFlag = 1;
            ai->geom = cAtomInfoNone;
            ai->valence = 0;
            break;
          case cAN_H:
          case cAN_F:
          case cAN_I:
          case cAN_Br:
            ai->chemFlag = 1;
            ai->geom = cAtomInfoSingle;
            ai->valence = 1;
            break;
          case cAN_O: {
            auto const neighbors = AtomNeighbors(I, a);
            auto const nn = neighbors.size();
            if (nn != 1) { /* water, hydroxy, ether */
              ai->chemFlag = 1;
              ai->geom = cAtomInfoTetrahedral;
              ai->valence = 2;
            } else { /* hydroxy or carbonyl? check carbon geometry */
              ai2 = I->AtomInfo.data() + neighbors[0].atm;
              if (ai2->chemFlag) {
                if ((ai2->geom == cAtomInfoTetrahedral) ||
                    (ai2->geom == cAtomInfoLinear)) {
                  ai->chemFlag = 1; /* hydroxy */
                  ai->geom = cAtomInfoTetrahedral;
                  ai->valence = 2;
                }
              }
            }
            break;
          }
          case cAN_C:
            if (geom >= 0) {
              ai->geom = geom;
              ai->valence = carbonVal[geom];
              ai->chemFlag = true;
            } else {
              auto const neighbors = AtomNeighbors(I, a);
              auto const nn = neighbors.size();
              if (nn == 1) { /* only one neighbor */
                ai2 = I->AtomInfo.data() + neighbors[0].atm;
                if (ai2->chemFlag && (ai2->geom == cAtomInfoTetrahedral)) {
                  ai->chemFlag =
                      true; /* singleton carbon bonded to tetC must be tetC */
                  ai->geom = cAtomInfoTetrahedral;
                  ai->valence = 4;
                }
              }
            }
            break;
          case cAN_N:
            if (geom == cAtomInfoPlanar) {
              ai->chemFlag = true;
              ai->geom = cAtomInfoPlanar;
              ai->valence = 3;
            } else if (geom == cAtomInfoTetrahedral) {
              ai->chemFlag = true;
              ai->geom = cAtomInfoTetrahedral;
              ai->valence = 4;
            }
            break;
          case cAN_S: {
            auto const nn = AtomNeighbors(I, a).size();
            if (nn == 4) { /* sulfone */
              ai->chemFlag = true;
              ai->geom = cAtomInfoTetrahedral;
              ai->valence = 4;
            } else if (nn == 3) { /* suloxide */
              ai->chemFlag = true;
              ai->geom = cAtomInfoTetrahedral;
              ai->valence = 3;
            } else if (nn == 2) { /* thioether */
              ai->chemFlag = true;
              ai->geom = cAtomInfoTetrahedral;
              ai->valence = 2;
            }
            break;
          }
          case cAN_Cl:
            ai->chemFlag = 1;
            if (ai->formalCharge == 0) {
              ai->geom = cAtomInfoSingle;
              ai->valence = 1;
            } else {
              ai->geom = cAtomInfoNone;
              ai->valence = 0;
            }
            break;
        }
        if (ai->chemFlag) changedFlag = true;
      }
    }
  }
}
/**
 * Assigns, based on valence, geom, formalCharge, and bonded atoms:
 * - hb_donor
 * - hb_acceptor
 */
void ObjectMoleculeInferHBondFromChem(ObjectMolecule* I) {
  int a;
  AtomInfoType* ai;
  int a1;
  int has_hydro;
  /* initialize accumulators on uncategorized atoms */

  const lexborrow_t lex_pseudo = LexBorrow(I->G, "pseudo");

  ai = I->AtomInfo.data();
  for (a = 0; a < I->NAtom; a++) {
    auto const neighbors = AtomNeighbors(I, a);
    int nn = neighbors.size();
    ai->hb_donor = false;
    ai->hb_acceptor = false;

    has_hydro = (nn < ai->valence); /* implicit hydrogens? */

    if (!has_hydro) {
      /* explicit hydrogens? */
      switch (ai->protons) {
        case cAN_N:
        case cAN_O:
          for (auto const& neighbor : neighbors) {
            a1 = neighbor.atm;
            if (I->AtomInfo[a1].protons == 1) {
              has_hydro = true;
              break;
            }
            if (I->AtomInfo[a1].name == lex_pseudo && --nn < ai->valence) {
              has_hydro = true;
              break;
            }
          }
          break;
      }
    }

    switch (ai->protons) {
        /* cat-ions, lewis acids etc. */
      case cAN_Fe:
      case cAN_Ca:
      case cAN_Cu:
      case cAN_K:
      case cAN_Na:
      case cAN_Mg:
      case cAN_Zn:
      case cAN_Hg:
      case cAN_Sr:
      case cAN_Ba:
        ai->hb_donor = true;
        break;
      case cAN_N:
        if (has_hydro)
          ai->hb_donor = true;
        else {
          int delocalized = false;
          int has_double_bond = false;
          int neighbor_has_double_bond = false;
          int mem[3];
          int nbr[3];
          int const* neighbor = I->getNeighborArray();

          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while (((mem[1] = neighbor[nbr[0]]) >= 0)) {
            int b_order = I->Bond[neighbor[nbr[0] + 1]].order;
            if (b_order > 1) { /* any double/triple/aromatic bonds? */
              delocalized = true;
            }
            if (b_order == 2) {
              has_double_bond = true;
            }
            nbr[1] = neighbor[mem[1]] + 1;
            while (((mem[2] = neighbor[nbr[1]]) >= 0)) {
              if (mem[2] != mem[0]) {
                int b_order2 = I->Bond[neighbor[nbr[1] + 1]].order;
                if (b_order2 == 2) {
                  neighbor_has_double_bond = true;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }
          if ((ai->formalCharge <= 0) && delocalized && (nn < 3)) {
            /* delocalized nitrogen can likely serve as an acceptor */
            ai->hb_acceptor = true;
          }
          if (delocalized && (neighbor_has_double_bond) && (!has_double_bond) &&
              (ai->geom == cAtomInfoPlanar) && (nn == 2) &&
              (ai->formalCharge >= 0)) {
            /* there's a fair chance of a resonance structure with this nitrogen
             * as a donor */
            ai->hb_donor = true;
          }
          if ((ai->geom != cAtomInfoPlanar) && (nn == 3) &&
              (ai->formalCharge >= 0) && (!delocalized)) {
            /* tertiary amine case -- assume potential donor */
            ai->hb_donor = true;
          }
        }
        break;
      case cAN_O:
        if (ai->formalCharge <= 0) ai->hb_acceptor = true;
        if (has_hydro)
          ai->hb_donor = true;
        else {
          int has_double_bond = false;
          int neighbor_has_aromatic_bond = false;
          int mem[3];
          int nbr[3];
          auto* const neighbor = I->getNeighborArray();

          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while (((mem[1] = neighbor[nbr[0]]) >= 0)) {
            int b_order = I->Bond[neighbor[nbr[0] + 1]].order;
            if (b_order == 2) {
              has_double_bond = true;
            }
            nbr[1] = neighbor[mem[1]] + 1;
            while (((mem[2] = neighbor[nbr[1]]) >= 0)) {
              if (mem[2] != mem[0]) {
                int b_order2 = I->Bond[neighbor[nbr[1] + 1]].order;
                if (b_order2 == 4) {
                  neighbor_has_aromatic_bond = true;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }

          if (has_double_bond && neighbor_has_aromatic_bond &&
              (ai->formalCharge >= 0)) {
            /* allow for phenolic resonance structures (and the like) */
            ai->hb_donor = true;
          }
        }
        break;
    }
    ai++;
  }
}
/**
 * Assigns:
 * - geom
 * - valence
 * - chemFlag
 */
void ObjectMoleculeInferChemFromBonds(ObjectMolecule* I, int state) {
  int a, b;
  const BondType* b0;
  AtomInfoType *ai, *ai0, *ai1 = nullptr;
  int a0, a1;
  int expect, order;
  int changedFlag;
  /* initialize accumulators on uncategorized atoms */

  ai = I->AtomInfo.data();
  for (a = 0; a < I->NAtom; a++) {
    if (!ai->chemFlag) {
      ai->geom = 0;
      ai->valence = 0;
    }
    ai++;
  }

  // Ignore "pseudo" atoms (Desmond virtual sites)
  const lexborrow_t lex_pseudo = LexBorrow(I->G, "pseudo");
  std::vector<unsigned char> pseudo_neighbor_count(
      lex_pseudo == LEX_BORROW_NOTFOUND ? 0 : I->NAtom);

  /* find maximum bond order for each atom */

  b0 = I->Bond;
  for (b = 0; b < I->NBond; b++) {
    a0 = b0->index[0];
    a1 = b0->index[1];
    ai0 = I->AtomInfo + a0;
    ai1 = I->AtomInfo + a1;
    order = b0->order;
    b0++;

    // count "pseudo" neighbors
    if (!pseudo_neighbor_count.empty()) {
      if (ai0->name == lex_pseudo) {
        pseudo_neighbor_count[a1] += 1;
      }
      if (ai1->name == lex_pseudo) {
        pseudo_neighbor_count[a0] += 1;
      }
    }

    if (!ai0->chemFlag) {
      if (order > ai0->geom) ai0->geom = order;
      ai0->valence += order;
    }
    if (!ai1->chemFlag) {
      if (order > ai1->geom) ai1->geom = order;
      ai1->valence += order;
    }
    if (order == 3) {
      /* override existing chemistry * this is a temp fix to a pressing
         problem... we need to rethink the chemisty assignment ordering (should
         bond information come first? */
      ai0->geom = cAtomInfoLinear;
      ai1->geom = cAtomInfoLinear;
      if (ai0->chemFlag != 2) {
        switch (ai0->protons) {
          case cAN_C:
            ai0->valence = 2;
            break;
          default:
            ai0->valence = 1;
        }
        ai0->chemFlag = true;
      }
      if (ai1->chemFlag != 2) {
        switch (ai1->protons) {
          case cAN_C:
            ai1->valence = 2;
            break;
          default:
            ai1->valence = 1;
        }
        ai1->chemFlag = true;
      }
    } else if (order == 4) {
      ai0->geom = cAtomInfoPlanar;
      ai1->geom = cAtomInfoPlanar;
      if (ai0->chemFlag != 2) {
        switch (ai0->protons) {
          case cAN_O:
            ai0->valence = 1;
            break;
          case cAN_N:
            ai0->valence = (ai0->formalCharge == 1)
                               ? 3
                               : 2;  // TODO wrong for neutral -NH-
            break;
          case cAN_C:
            ai0->valence = 3;
            break;
          case cAN_S:
            ai0->valence = 2;
            break;
          default:
            ai0->valence = 4;
        }
        ai0->chemFlag = true;
      }
      if (ai1->chemFlag != 2) {
        switch (ai1->protons) {
          case cAN_O:
            ai1->valence = 1;
            break;
          case cAN_N:
            ai1->valence = (ai1->formalCharge == 1)
                               ? 3
                               : 2;  // TODO wrong for neutral -NH-
            break;
          case cAN_C:
            ai1->valence = 3;
            break;
          default:
            ai1->valence = 1;
        }
        ai1->chemFlag = true;
      }
    }
  }

  /* now set up valences and geometries */

  ai = I->AtomInfo.data();
  for (a = 0; a < I->NAtom; a++) {
    if (!ai->chemFlag) {
      expect = AtomInfoGetExpectedValence(I->G, ai);
      int nn = AtomNeighbors(I, a).size();

      // don't count "pseudo" atom neighbors
      if (!pseudo_neighbor_count.empty()) {
        nn -= pseudo_neighbor_count[a];
      }

      if (ai->geom == 3) {
        ai->geom = cAtomInfoLinear;
        switch (ai->protons) {
          case cAN_C:
            ai->valence = 2;
            break;
          default:
            ai->valence = 1;
        }
        ai->chemFlag = true;
      } else {
        if (expect < 0) expect = -expect; /* for now, just ignore this issue */
        /*      printf("%d %d %d %d\n",ai->geom,ai->valence,nn,expect); */
        if (ai->valence == expect) { /* sum of bond orders equals valence */
          ai->chemFlag = true;
          ai->valence = nn;
          switch (ai->geom) { /* max bond order observed */
            case 0:
              ai->geom = cAtomInfoNone;
              break;
            case 2:
              ai->geom = cAtomInfoPlanar;
              break;
            case 3:
              ai->geom = cAtomInfoLinear;
              break;
            default:
              if (expect == 1)
                ai->geom = cAtomInfoSingle;
              else
                ai->geom = cAtomInfoTetrahedral;
              break;
          }
        } else if (ai->valence < expect) { /* missing a bond */
          ai->chemFlag = true;
          ai->valence = nn + (expect - ai->valence);
          switch (ai->geom) {
            case 2:
              ai->geom = cAtomInfoPlanar;
              break;
            case 3:
              ai->geom = cAtomInfoLinear;
              break;
            default:
              if (expect == 1)
                ai->geom = cAtomInfoSingle;
              else
                ai->geom = cAtomInfoTetrahedral;
              break;
          }
        } else if (ai->valence > expect) {
          ai->chemFlag = true;
          ai->valence = nn;
          switch (ai->geom) {
            case 2:
              ai->geom = cAtomInfoPlanar;
              break;
            case 3:
              ai->geom = cAtomInfoLinear;
              break;
            default:
              if (expect == 1)
                ai->geom = cAtomInfoSingle;
              else
                ai->geom = cAtomInfoTetrahedral;
              break;
          }
          if (nn > 3) ai->geom = cAtomInfoTetrahedral;
        }
      }
    }
    ai++;
  }

  /* now go through and make sure conjugated amines are planer */
  changedFlag = true;
  while (changedFlag) {
    changedFlag = false;
    ai = I->AtomInfo.data();
    for (a = 0; a < I->NAtom; a++) {
      if (ai->chemFlag) {
        if (ai->protons == cAN_N)
          if (ai->formalCharge < 1)
            if (ai->geom == cAtomInfoTetrahedral) {
              /* search for uncharged tetrahedral nitrogen */
              for (auto const& neighbor : AtomNeighbors(I, a)) {
                auto const* ai0 = I->AtomInfo.data() + neighbor.atm;
                if ((ai0->chemFlag) && (ai0->geom == cAtomInfoPlanar) &&
                    ((ai0->protons == cAN_C) || (ai0->protons == cAN_N))) {
                  ai->geom =
                      cAtomInfoPlanar; /* found probable delocalization */
                  if (ai->formalCharge == 0)
                    ai->valence = 3; /* just in case... */
                  changedFlag = true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }

  /* now go through and make sure conjugated anions are planer */
  changedFlag = true;
  while (changedFlag) {
    changedFlag = false;
    ai = I->AtomInfo.data();
    for (a = 0; a < I->NAtom; a++) {
      if (ai->chemFlag) {
        if (ai->protons == cAN_O)
          if (ai->formalCharge == -1)
            if ((ai->geom == cAtomInfoTetrahedral) ||
                (ai->geom == cAtomInfoSingle)) {
              /* search for anionic tetrahedral oxygen */
              for (auto const& neighbor : AtomNeighbors(I, a)) {
                auto const* ai0 = I->AtomInfo + neighbor.atm;
                if ((ai0->chemFlag) && (ai0->geom == cAtomInfoPlanar) &&
                    ((ai0->protons == cAN_C) || (ai0->protons == cAN_N))) {
                  ai->geom =
                      cAtomInfoPlanar; /* found probable delocalization */
                  changedFlag = true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }
}
float ObjectMoleculeGetAvgHBondVector(ObjectMolecule* I, int atom, int state,
                                      float* v, float* incoming)

/* computes average / optima hydrogen bonding vector for an atom */
{
  float result = 0.0;
  int vec_cnt = 0;
  float v_atom[3], v_neigh[3], v_diff[3], v_acc[3] = {0.0, 0.0, 0.0};
  int sp2_flag = false;
  int order;

  auto const* cs = I->getCoordSet(state);
  if (cs) {
    if (CoordSetGetAtomVertex(cs, atom,
                              v_atom)) { /* atom exists in this C-set */
      for (auto const& neighbor : AtomNeighbors(I, atom)) {
        auto const a2 = neighbor.atm;
        order = I->Bond[neighbor.bond].order;
        if ((order == 2) || (order == 4)) {
          sp2_flag = true;
        }

        const auto& ai = I->AtomInfo[atom];

        bool isH = I->AtomInfo[a2].protons == cAN_H;

        if (!isH ||
            (ai.protons == cAN_O && isH)) { /* ignore hydrogens unless water*/
          if (CoordSetGetAtomVertex(cs, a2, v_neigh)) {
            subtract3f(v_atom, v_neigh, v_diff);
            normalize3f(v_diff);
            add3f(v_diff, v_acc, v_acc);
            vec_cnt++;
          }
        }
      }
      if (vec_cnt) {
        result = (float)length3f(v_acc);
        result = result / vec_cnt;
        normalize23f(v_acc, v);
      } else {
        copy3f(v_acc, v);
      }

      if (incoming && (vec_cnt == 1) &&
          (fabs(dot_product3f(v, incoming)) < 0.99F)) {
        /* if we know where the donor is, and the acceptor can
           potentially rotate the lone pair, then we should optimally
           orient the acceptor, if possible */
        AtomInfoType* ai = I->AtomInfo + atom;
        float v_perp[3];
        float v_tmp1[3], v_tmp2[3];
        if (((ai->protons == cAN_O) && (!sp2_flag)) || /* C-O-H */
            ((ai->protons == cAN_N) && (sp2_flag))) {  /* C=N-H */

          remove_component3f(incoming, v, v_perp);
          normalize3f(v_perp);
          scale3f(v, 0.333644F, v_tmp1);
          scale3f(v_perp, 0.942699F, v_tmp2);
          add3f(v_tmp1, v_tmp2, v_tmp2);
          subtract3f(v, v_tmp2, v);
          normalize3f(v);
        }
      }
    }
  }
  return (result);
}
