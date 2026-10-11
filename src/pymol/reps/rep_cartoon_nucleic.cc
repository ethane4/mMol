
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

// Cartoons of nucleic acids: backbones, caps and base rings.

#include "pymol/reps/rep_cartoon.h"
#include "pymol/reps/rep_cartoon_internal.h"

#define ESCAPE_MAX 500

#define NUCLEIC_NORMAL0 "C2"

struct BondPathStep {
  const char* name;
  const char* name_alt;  // e.g. star variant, or nullptr
  int protons;
};

/**
 * @brief follows a bond path from `start` through atoms matching the given
 * names. Each step finds a bonded neighbor matching name/name_alt that hasn't
 * been visited.
 * @param obj of atoms whose bonds are followed
 * @param start index of starting atom
 * @param steps sequence of Bond paths to navigate
 * @param marked per-atom flags for atoms already claimed by a ring; skipped
 * during traversal
 * @return the atom index at the end of the path, or -1 if not found.
 */
template <size_t N>
static int FollowBondPath(const ObjectMolecule* obj, int start,
                          const std::array<BondPathStep, N>& steps,
                          const int* marked) {
  auto* atomInfo = obj->AtomInfo.data();
  int prev = -1;
  int cur = start;

  for (const auto& step : steps) {
    int next = -1;
    for (const auto& nbr : AtomNeighbors(obj, cur)) {
      if (nbr.atm != prev && atomInfo[nbr.atm].protons == step.protons &&
          (!marked || !marked[nbr.atm])) {
        auto aname = LexStr(obj->G, atomInfo[nbr.atm].name);
        if (WordMatchExact(obj->G, step.name, aname, 1) ||
            (step.name_alt &&
             WordMatchExact(obj->G, step.name_alt, aname, 1))) {
          next = nbr.atm;
          break;
        }
      }
    }
    if (next < 0) {
      return -1;
    }
    prev = cur;
    cur = next;
  }
  return cur;
}

#define MAX_RING_ATOM 10

/**
 * Return true if a connector between the two atoms should be drawn.
 *
 * TODO: should RepSphere and side_chain_helper check really be part of this?
 * TODO: should line|cyl|sphere also be checked for at2?
 */
static bool ring_connector_visible(PyMOLGlobals* G, const AtomInfoType* ai1,
                                   const AtomInfoType* ai2, bool sc_helper) {
  if (!(ai1->visRep & ai2->visRep & cRepCartoonBit)) return false;

  if (!(ai1->visRep & (cRepLineBit | cRepCylBit | cRepSphereBit))) return true;

  return !(
      AtomSettingGetWD(G, ai1, cSetting_cartoon_side_chain_helper, sc_helper) ||
      AtomSettingGetWD(G, ai2, cSetting_cartoon_side_chain_helper, sc_helper));
}

/**
 * depth-first neighbor search for nucleic acid atom
 *
 * @param nuc_flag    NAtom-length boolean nucleic acid flags array
 * @param obj         Molecule
 * @param atm         atom index
 * @param max_depth   maximum search depth
 * @param seen        set of visited atom indices (read/write)
 *
 * @return True if any nucleic acid atom found within max_depth radius
 */
static bool has_nuc_neighbor(const int* nuc_flag, const ObjectMolecule* obj,
                             const int atm, const int max_depth,
                             std::set<int>& seen) {
  for (auto const& neighbor : AtomNeighbors(obj, atm)) {
    auto const atm_neighbor = neighbor.atm;

    if (nuc_flag[atm_neighbor]) return true;

    if (seen.count(atm_neighbor)) continue;

    seen.insert(atm_neighbor);

    if (max_depth > 1 &&
        has_nuc_neighbor(nuc_flag, obj, atm_neighbor, max_depth - 1, seen))
      return true;
  }

  return false;
}

/* atix must contain n_atom + 1 elements, with the first atom repeated at the
 * end */

static void do_ring(PyMOLGlobals* G, nuc_acid_data* ndata, int n_atom,
                    int* atix, ObjectMolecule* obj, CoordSet* cs,
                    float ring_width, CGO* cgo, int ring_color,
                    float ladder_radius, int ladder_color, int ladder_mode,
                    int sc_helper, float ring_alpha, float alpha, int* marked,
                    float* moved, float ring_radius) {
  const float* v_i[MAX_RING_ATOM];
  const float* col[MAX_RING_ATOM];
  float n_up[MAX_RING_ATOM][3];
  float n_dn[MAX_RING_ATOM][3];
  const AtomInfoType* ai_i[MAX_RING_ATOM];
  int have_all = true;
  int all_marked = true;
  const AtomInfoType* ai;
  int have_C4 = -1;
  int have_C4_prime = -1;
  int have_C_number = -1;
  int nf = false;
  int ring_mode = ndata->ring_mode;
  int finder = ndata->ring_finder;
  const auto& nuc_flag = ndata->nuc_flag;

  /* first, make sure all atoms have known coordinates */
  {
    int a, i;
    for (i = 0; i <= n_atom; i++) {
      int a1 = atix[i];
      int have_atom = false;
      if (nuc_flag[a1]) nf = true;
      a = cs->atmToIdx(a1);
      if (a >= 0) {
        ai = obj->AtomInfo + a1;
        if (ai->visRep & cRepCartoonBit) {
          ai_i[i] = ai;

          // take atom level settings from any ring atom (effectifly from the
          // last one with settings)
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ring_mode,
                                  &ring_mode);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ring_color,
                                  &ring_color);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ring_radius,
                                  &ring_radius);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ring_width,
                                  &ring_width);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ladder_mode,
                                  &ladder_mode);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ladder_color,
                                  &ladder_color);
          AtomSettingGetIfDefined(G, ai, cSetting_cartoon_ladder_radius,
                                  &ladder_radius);

          col[i] = ColorGet(G, ai->color);
          v_i[i] = cs->coordPtr(a);
          have_atom = true;
          const char* ai_name = LexStr(G, ai->name);
          if (WordMatchExact(G, "C4", ai_name, 1))
            have_C4 = a1;
          else if (WordMatchExact(G, "C4'", ai_name, 1) ||
                   WordMatchExact(G, "C4*", ai_name, 1))
            have_C4_prime = a1;
          if (((ai_name[0] == 'C') || (ai_name[0] == 'c')) &&
              isdigit(ai_name[1]))
            have_C_number = a1;
        }
        if (!marked[a1]) all_marked = false;
      }
      if (!have_atom) {
        have_all = false;
        break;
      }
    }
  }

  if (!ring_mode) finder = 0;

  ring_width *= 0.5F;

  if (n_atom && have_all && (!all_marked)) {
    if (ladder_mode) {
      int i;
      int a1;
      const AtomInfoType* ai2;
      const AtomInfoType* atomInfo = obj->AtomInfo.data();
      int mem0, mem1, mem2, mem3, mem4, mem5, mem6, mem7;
      auto* const neighbor = obj->getNeighborArray();
      int nbr[7];
      int sugar_at = -1, base_at = -1;
      int phos3_at = -1, phos5_at = -1;
      int o3_at = -1, o5_at = -1, c1_at = -1;
      int purine_flag = false;
      int c1 = -1, c1_linked = -1;
      int c2 = -1, c2_linked = -1;
      int c3 = -1, c3_linked = -1;
      int c4 = -1, c4_linked = -1;
      int c5 = -1, c5_linked = -1;

      if (n_atom == 5) { /* going from the sugar to the base */
        for (i = 0; i < n_atom; i++) {
          a1 = atix[i];
          if (!nf) nf = nuc_flag[a1];
          ai = atomInfo + a1;
          if ((ai->protons == cAN_C) && (!marked[a1]) &&
              (WordMatchExact(G, "C3*", LexStr(G, ai->name), 1) ||
               WordMatchExact(G, "C3'", LexStr(G, ai->name), 1))) {
            sugar_at = a1;
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while ((mem1 = neighbor[nbr[0]]) >= 0) {
              if ((atomInfo[mem1].protons == cAN_O) && (!marked[mem1])) {
                ai = atomInfo + mem1;
                if (WordMatchExact(G, "O3*", LexStr(G, ai->name), 1) ||
                    WordMatchExact(G, "O3'", LexStr(G, ai->name), 1))
                  o3_at = mem1;
                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (!marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_P)) {
                    phos3_at = mem2;
                  }
                  nbr[1] += 2;
                }
              }

              if ((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
                ai2 = atomInfo + mem1;
                if (WordMatchExact(G, NUCLEIC_NORMAL1, LexStr(G, ai2->name),
                                   1) ||
                    WordMatchExact(G, NUCLEIC_NORMAL2, LexStr(G, ai2->name), 1))
                  sugar_at = mem1;

                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (!marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_C)) {
                    ai = atomInfo + mem2;
                    if (WordMatchExact(G, "C1*", LexStr(G, ai->name), 1) ||
                        WordMatchExact(G, "C1'", LexStr(G, ai->name), 1))
                      c1_at = mem2;
                    nbr[2] = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr[2]]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0)) {
                        if ((atomInfo[mem3].protons == cAN_O) &&
                            (!marked[mem3])) {
                          ai = atomInfo + mem3;
                          if (WordMatchExact(G, "O5*", LexStr(G, ai->name),
                                             1) ||
                              WordMatchExact(G, "O5'", LexStr(G, ai->name), 1))
                            o5_at = mem3;

                          nbr[3] = neighbor[mem3] + 1;
                          while ((mem4 = neighbor[nbr[3]]) >= 0) {
                            if ((mem4 != mem2) && (!marked[mem4]) &&
                                (atomInfo[mem4].protons == cAN_P)) {
                              phos5_at = mem4;
                            }
                            nbr[3] += 2;
                          }
                        }
                        if (atomInfo[mem3].protons == cAN_N) {
                          if (ring_mode) {
                            ai2 = atomInfo + mem3;
                            if ((!marked[mem3]) &&
                                (WordMatchExact(G, "N1", LexStr(G, ai2->name),
                                                1) ||
                                 WordMatchExact(G, "N9", LexStr(G, ai2->name),
                                                1))) {
                              base_at = mem3;
                              if (ring_mode != 3) {
                                ladder_radius = ring_width * 1.5;
                              } else {
                                ladder_radius = ring_width * 3;
                              }
                            }
                          } else {
                            nbr[3] = neighbor[mem3] + 1;
                            while ((mem4 = neighbor[nbr[3]]) >= 0) {
                              if ((mem4 != mem2) && (mem4 != mem1) &&
                                  (mem4 != mem0) &&
                                  (atomInfo[mem4].protons == cAN_C)) {
                                nbr[4] = neighbor[mem4] + 1;
                                while ((mem5 = neighbor[nbr[4]]) >= 0) {
                                  if ((mem5 != mem3) && (mem5 != mem2) &&
                                      (mem5 != mem1) && (mem5 != mem0) &&
                                      (atomInfo[mem5].protons == cAN_N) &&
                                      (marked[mem5])) {
                                    /* must be in a mapped ring */
                                    /* clear flag here */
                                    purine_flag = false;

                                    nbr[5] = neighbor[mem5] + 1;
                                    while ((mem6 = neighbor[nbr[5]]) >= 0) {
                                      if ((mem6 != mem4) && (mem6 != mem3) &&
                                          (mem6 != mem2) && (mem6 != mem1) &&
                                          (mem6 != mem0) &&
                                          (atomInfo[mem6].protons == cAN_C) &&
                                          (marked[mem6])) {
                                        nbr[6] = neighbor[mem6] + 1;
                                        while ((mem7 = neighbor[nbr[6]]) >= 0) {
                                          ai2 = atomInfo + mem7;
                                          if ((mem7 != mem5) &&
                                              (mem7 != mem4) &&
                                              (mem7 != mem3) &&
                                              (mem7 != mem2) &&
                                              (mem7 != mem2) &&
                                              (mem7 != mem1) &&
                                              (mem7 != mem0) &&
                                              (ai2->protons == cAN_N) &&
                                              (marked[mem7])) {
                                            if (WordMatchExact(
                                                    G, "N1",
                                                    LexStr(G, ai2->name), 1)) {
                                              /* and set flag */
                                              base_at = mem7;
                                              purine_flag = true;
                                            }
                                          }
                                          nbr[6] += 2;
                                        }
                                      }
                                      nbr[5] += 2;
                                    }

                                    if (!purine_flag) {
                                      ai2 = atomInfo + mem5;
                                      if (marked[mem5] &&
                                          WordMatchExact(G, "N3",
                                                         LexStr(G, ai2->name),
                                                         1)) {
                                        base_at = mem5;
                                      }
                                    }
                                  }
                                  nbr[4] += 2;
                                }
                              }
                              nbr[3] += 2;
                            }
                          }
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              }
              nbr[0] += 2;
            }
          }
        }
      } else if ((!ring_mode) &&
                 (n_atom == 6)) { /* going from the base to the sugar */

        for (i = 0; i < n_atom; i++) {
          a1 = atix[i];
          if (!nf) nf = nuc_flag[a1];
          ai = atomInfo + a1;

          /* base-hunting */

          if ((ai->protons == cAN_N) && (!marked[a1]) &&
              (WordMatchExact(G, "N1", LexStr(G, ai->name), 1) ||
               WordMatchExact(G, "N3", LexStr(G, ai->name), 1))) {
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while ((mem1 = neighbor[nbr[0]]) >= 0) {
              if ((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (!marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_N)) {
                    nbr[2] = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr[2]]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0) &&
                          (atomInfo[mem3].protons == cAN_C)) {
                        nbr[3] = neighbor[mem3] + 1;
                        while ((mem4 = neighbor[nbr[3]]) >= 0) {
                          if ((mem4 != mem2) && (mem4 != mem1) &&
                              (mem4 != mem0)) {
                            if ((atomInfo[mem4].protons == cAN_N) ||
                                WordMatchExact(G, "C5",
                                               LexStr(G, atomInfo[mem4].name),
                                               1)) { /* purine case */
                              nbr[4] = neighbor[mem4] + 1;
                              while ((mem5 = neighbor[nbr[4]]) >= 0) {
                                if ((mem5 != mem3) && (mem5 != mem2) &&
                                    (mem5 != mem1) && (mem5 != mem0) &&
                                    (marked[mem5]) &&
                                    (atomInfo[mem5].protons == cAN_C)) {
                                  nbr[5] = neighbor[mem5] + 1;
                                  while ((mem6 = neighbor[nbr[5]]) >= 0) {
                                    if ((mem6 != mem4) && (mem6 != mem3) &&
                                        (mem6 != mem2) && (mem6 != mem1) &&
                                        (mem6 != mem0) &&
                                        ((atomInfo[mem6].protons == cAN_C) ||
                                         (atomInfo[mem6].protons == cAN_O)) &&
                                        (marked[mem6])) {
                                      nbr[6] = neighbor[mem6] + 1;
                                      while ((mem7 = neighbor[nbr[6]]) >= 0) {
                                        ai2 = atomInfo + mem7;
                                        if ((mem7 != mem5) && (mem7 != mem4) &&
                                            (mem7 != mem3) && (mem7 != mem2) &&
                                            (mem7 != mem2) && (mem7 != mem1) &&
                                            (mem7 != mem0) &&
                                            (ai2->protons == cAN_C) &&
                                            (marked[mem7])) {
                                          if (WordMatchExact(
                                                  G, NUCLEIC_NORMAL1,
                                                  LexStr(G, ai2->name), 1) ||
                                              WordMatchExact(
                                                  G, NUCLEIC_NORMAL2,
                                                  LexStr(G, ai2->name), 1)) {
                                            base_at = a1;
                                            sugar_at = mem7;
                                          }
                                        }
                                        nbr[6] += 2;
                                      }
                                    }
                                    nbr[5] += 2;
                                  }
                                }
                                nbr[4] += 2;
                              }
                            } else if (((atomInfo[mem4].protons == cAN_C) ||
                                        (atomInfo[mem4].protons == cAN_O)) &&
                                       (marked[mem4])) {
                              /* pyrimidine case */

                              nbr[4] = neighbor[mem4] + 1;
                              while ((mem5 = neighbor[nbr[4]]) >= 0) {
                                ai2 = atomInfo + mem5;
                                if ((mem5 != mem3) && (mem5 != mem2) &&
                                    (mem5 != mem1) && (mem5 != mem0) &&
                                    (ai2->protons == cAN_C) && (marked[mem5])) {
                                  if (WordMatchExact(G, NUCLEIC_NORMAL1,
                                                     LexStr(G, ai2->name), 1) ||
                                      WordMatchExact(G, NUCLEIC_NORMAL2,
                                                     LexStr(G, ai2->name), 1)) {
                                    base_at = a1;
                                    sugar_at = mem5;
                                  }
                                }
                                nbr[4] += 2;
                              }
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
              }
              nbr[0] += 2;
            }
          }
        }
      }
      /* PNA: find base_at and sugar_at via N1/N9 -> C8' -> C7' -> N4' path */
      if (sugar_at < 0 && base_at < 0 && (n_atom == 5 || n_atom == 6)) {
        static const std::array<BondPathStep, 3> pna_path = {{
            {"C8'", "C8*", cAN_C},
            {"C7'", "C7*", cAN_C},
            {"N4'", "N4*", cAN_N},
        }};
        for (i = 0; i < n_atom; i++) {
          a1 = atix[i];
          ai = atomInfo + a1;
          if (ai->protons == cAN_N && !marked[a1] &&
              (WordMatchExact(G, "N1", LexStr(G, ai->name), 1) ||
               WordMatchExact(G, "N9", LexStr(G, ai->name), 1))) {
            int n4_at = FollowBondPath(obj, a1, pna_path, marked);
            if (n4_at >= 0) {
              base_at = a1;
              sugar_at = n4_at;
              if (!nf) nf = nuc_flag[a1];
              break;
            }
          }
        }
      }

      if (n_atom == 6) {
        for (i = 0; i < n_atom; i++) {
          a1 = atix[i];
          ai = atomInfo + a1;

          /* glycosylation hunting */

          if ((ai->protons == cAN_C) &&
              (!marked[a1])) { /* unmarked C in ring */
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while ((mem1 = neighbor[nbr[0]]) >= 0) {
              if ((atomInfo[mem1].protons == cAN_C) && /* exocyclic C */
                  (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (!marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_O)) {
                    /* exocyclic O */
                    nbr[2] = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr[2]]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0) && (marked[mem3]) &&
                          (atomInfo[mem3].protons == cAN_C)) {
                        /* cyclic C */
                        if (WordMatchExact(G, "C5", LexStr(G, ai->name), 1) &&
                            WordMatchExact(G, "C6",
                                           LexStr(G, atomInfo[mem1].name), 1)) {
                          c5_linked = mem3;
                          c5 = mem0;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              } else if ((atomInfo[mem1].protons == cAN_O) && /* exocyclic O */
                         (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_C)) {
                    /* cyclic C */
                    const char* ai_name = LexStr(G, ai->name);
                    if (WordMatchExact(G, "C1", ai_name, 1)) {
                      c1_linked = mem2;
                      c1 = mem0;
                    } else if (WordMatchExact(G, "C2", ai_name, 1)) {
                      c2_linked = mem2;
                      c2 = mem0;
                    } else if (WordMatchExact(G, "C3", ai_name, 1)) {
                      c3_linked = mem2;
                      c3 = mem0;
                    } else if (WordMatchExact(G, "C4", ai_name, 1)) {
                      c4_linked = mem2;
                      c4 = mem0;
                    }
                  } else if ((mem2 != mem0) && (!marked[mem2]) &&
                             (atomInfo[mem2].protons == cAN_C)) {
                    /* exocyclic C */
                    nbr[2] = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr[2]]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0) && (marked[mem3]) &&
                          (atomInfo[mem3].protons == cAN_C)) {
                        /* cyclic C */
                        if (WordMatchExact(G, "C5",
                                           LexStr(G, atomInfo[mem3].name), 1) &&
                            WordMatchExact(G, "C6",
                                           LexStr(G, atomInfo[mem2].name), 1)) {
                          c5 = mem0;
                          c5_linked = mem3;
                        }
                      } else if ((mem3 != mem1) && (mem3 != mem0) &&
                                 (!marked[mem3]) &&
                                 (atomInfo[mem3].protons == cAN_C)) {
                        /* exocyclic */
                        if (WordMatchExact(G, "C1", LexStr(G, ai->name), 1) &&
                            WordMatchExact(G, "CA",
                                           LexStr(G, atomInfo[mem3].name), 1)) {
                          c1 = mem0;
                          c1_linked = mem3;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              } else if ((atomInfo[mem1].protons == cAN_N) && /* exocyclic N */
                         (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while ((mem2 = neighbor[nbr[1]]) >= 0) {
                  if ((mem2 != mem0) && (!marked[mem2]) &&
                      (atomInfo[mem2].protons == cAN_C)) {
                    /* exocyclic C */
                    nbr[2] = neighbor[mem2] + 1;
                    while ((mem3 = neighbor[nbr[2]]) >= 0) {
                      if ((mem3 != mem1) && (mem3 != mem0) && (!marked[mem3]) &&
                          (atomInfo[mem3].protons == cAN_C)) {
                        /* exocyclic */
                        nbr[3] = neighbor[mem3] + 1;
                        while ((mem4 = neighbor[nbr[3]]) >= 0) {
                          if ((mem4 != mem2) && (mem3 != mem1) &&
                              (!marked[mem4]) &&
                              (atomInfo[mem4].protons == cAN_C)) {
                            /* exocyclic */
                            if (WordMatchExact(G, "C1", LexStr(G, ai->name),
                                               1) &&
                                WordMatchExact(G, "CA",
                                               LexStr(G, atomInfo[mem4].name),
                                               1)) {
                              c1 = mem0;
                              c1_linked = mem4;
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
              }
              nbr[0] += 2;
            }
          }
        }
      }

      if (sugar_at >= 0) { /* still need to find the phosphates... */
        int c3_index = -1;
        ai = atomInfo + sugar_at;
        if (WordMatchExact(G, "C3*", LexStr(G, ai->name), 1) ||
            WordMatchExact(G, "C3'", LexStr(G, ai->name), 1)) {
          c3_index = sugar_at;
        } else {
          mem0 = sugar_at;
          nbr[0] = neighbor[mem0] + 1;
          while ((mem1 = neighbor[nbr[0]]) >= 0) {
            if ((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
              ai = atomInfo + mem1;
              if (!(WordMatchExact(G, "C3*", LexStr(G, ai->name), 1) ||
                    WordMatchExact(G, "C3'", LexStr(G, ai->name), 1))) {
                c3_index = mem1;
              }
            }
            nbr[0] += 2;
          }
        }

        if (c3_index >= 0) { /* now we know where we are... */

          mem0 = c3_index;
          nbr[0] = neighbor[mem0] + 1;
          while ((mem1 = neighbor[nbr[0]]) >= 0) {
            if ((atomInfo[mem1].protons == cAN_O) && (!marked[mem1])) {
              ai = atomInfo + mem1;
              if (WordMatchExact(G, "O3*", LexStr(G, ai->name), 1) ||
                  WordMatchExact(G, "O3'", LexStr(G, ai->name), 1))
                o3_at = mem1;
              nbr[1] = neighbor[mem1] + 1;
              while ((mem2 = neighbor[nbr[1]]) >= 0) {
                if ((mem2 != mem0) && (!marked[mem2]) &&
                    (atomInfo[mem2].protons == cAN_P)) {
                  phos3_at = mem2;
                }
                nbr[1] += 2;
              }
            }
            if ((atomInfo[mem1].protons == cAN_C) && (marked[mem1])) {
              nbr[1] = neighbor[mem1] + 1;
              while ((mem2 = neighbor[nbr[1]]) >= 0) {
                if ((mem2 != mem0) && (!marked[mem2]) &&
                    (atomInfo[mem2].protons == cAN_C)) {
                  ai = atomInfo + mem2;
                  if (WordMatchExact(G, "C1*", LexStr(G, ai->name), 1) ||
                      WordMatchExact(G, "C1'", LexStr(G, ai->name), 1))
                    c1_at = mem2;

                  nbr[2] = neighbor[mem2] + 1;
                  while ((mem3 = neighbor[nbr[2]]) >= 0) {
                    if ((mem3 != mem1) && (mem3 != mem0) &&
                        (atomInfo[mem3].protons == cAN_O) && (!marked[mem3])) {
                      ai = atomInfo + mem3;
                      if (WordMatchExact(G, "O5*", LexStr(G, ai->name), 1) ||
                          WordMatchExact(G, "O5'", LexStr(G, ai->name), 1))
                        o5_at = mem3;
                      nbr[3] = neighbor[mem3] + 1;
                      while ((mem4 = neighbor[nbr[3]]) >= 0) {
                        if ((mem4 != mem2) && (!marked[mem4]) &&
                            (atomInfo[mem4].protons == cAN_P)) {
                          phos5_at = mem4;
                        }
                        nbr[3] += 2;
                      }
                    }
                    nbr[2] += 2;
                  }
                }
                nbr[1] += 2;
              }
            }
            nbr[0] += 2;
          }
        }
      }
      /* glycosylation connectors */
      if (1) {
        if ((ring_mode) && (finder >= 3)) {
          if (finder >= 3) {
            int b;
            float glyco_radius;
            switch (ring_mode) {
              case 3:
              case 4:
                glyco_radius = ring_width * 3;
                break;
              case 5:
                glyco_radius = ladder_radius;
                break;
              default:
                glyco_radius = ring_width * 1.5;
                break;
            }

            if ((alpha != 1.0F) || (ring_alpha != alpha)) CGOAlpha(cgo, alpha);

            for (b = 0; b < 5; b++) {
              int g1 = -1, g2 = -1;
              switch (b) {
                case 0:
                  g1 = c1;
                  g2 = c1_linked;
                  break;
                case 1:
                  g1 = c2;
                  g2 = c2_linked;
                  break;
                case 2:
                  g1 = c3;
                  g2 = c3_linked;
                  break;
                case 3:
                  g1 = c4;
                  g2 = c4_linked;
                  break;
                case 4:
                  g1 = c5;
                  g2 = c5_linked;
                  break;
              }

              if ((g1 >= 0) && (g2 >= 0)) {
                const AtomInfoType* g1_ai = atomInfo + g1;
                const AtomInfoType* g2_ai = atomInfo + g2;

                if (ring_connector_visible(G, g1_ai, g2_ai, sc_helper)) {
                  float avg[3];
                  int g1_x = cs->atmToIdx(g1);
                  int g2_x = cs->atmToIdx(g2);
                  const float* g1p = cs->coordPtr(g1_x);
                  const float* g2p = cs->coordPtr(g2_x);

                  if (!((!((ring_mode == 0) || (ring_mode == 4) ||
                           (ring_mode == 5))) ||
                        (!marked[g2]))) {
                    g2p = moved + 3 * g2;
                  }
                  if ((ring_mode == 0) || (ring_mode == 4) ||
                      (ring_mode == 5)) {
                    /* ring center */
                    int i;
                    /* compute average coordinate and mark atoms so that ring is
                     * only drawn once */
                    zero3f(avg);
                    for (i = 0; i < n_atom; i++) {
                      add3f(avg, v_i[i], avg);
                    }
                    scale3f(avg, 1.0F / n_atom, avg);
                    g1p = avg;
                  }

                  {
                    const float *color1, *color2;
                    if (ladder_color >= 0) {
                      color1 = color2 = ColorGet(G, ladder_color);
                    } else {
                      color1 = ColorGet(G, g1_ai->color);
                      color2 = ColorGet(G, g2_ai->color);
                    }
                    CGOPickColor(
                        cgo, g1,
                        g1_ai->masked ? cPickableNoPick : cPickableAtom);
                    Pickable pickcolor2 = {
                        static_cast<unsigned int>(g2),
                        g2_ai->masked ? cPickableNoPick : cPickableAtom};
                    float axis[3];
                    subtract3f(g2p, g1p, axis);
                    CGOColorv(cgo, color1);
                    float ladder_alpha =
                        1.0f - AtomSettingGetWD(G, ai_i[i],
                                                cSetting_cartoon_transparency,
                                                1.0f - alpha);
                    cgo->add<cgo::draw::shadercylinder2ndcolor>(
                        cgo, g1p, axis, glyco_radius, 0x1f, color2, &pickcolor2,
                        ladder_alpha);
                  }
                }
              }
            }
          }
        }
      }

      /* see if any of the neighbors are confirmed nucleic acids... */
      if (sugar_at >= 0) {
        if (!nf) {
          auto seen = std::set<int>();
          nf = has_nuc_neighbor(nuc_flag, obj, sugar_at, 5, seen);
        }

        if (nf) {
          if ((ring_mode) && ((finder == 1) || (finder >= 3))) {
            if ((c1_at >= 0) && (base_at >= 0)) {
              int save_at = sugar_at;
              sugar_at = c1_at;
              {
                const AtomInfoType* sug_ai = atomInfo + sugar_at;
                const AtomInfoType* bas_ai = atomInfo + base_at;

                if (ring_connector_visible(G, bas_ai, sug_ai, sc_helper)) {
                  int sug = cs->atmToIdx(sugar_at);
                  int bas = cs->atmToIdx(base_at);

                  if ((sug >= 0) && (bas >= 0)) {
                    {
                      const float *color1, *color2;
                      if (ladder_color >= 0) {
                        color1 = color2 = ColorGet(G, ladder_color);
                      } else {
                        color1 = ColorGet(G, sug_ai->color);
                        color2 = ColorGet(G, bas_ai->color);
                      }
                      CGOPickColor(
                          cgo, sugar_at,
                          sug_ai->masked ? cPickableNoPick : cPickableAtom);
                      Pickable pickcolor2 = {
                          static_cast<unsigned int>(base_at),
                          bas_ai->masked ? cPickableNoPick : cPickableAtom};
                      float axis[3];
                      subtract3f(cs->coordPtr(bas), cs->coordPtr(sug), axis);
                      CGOColorv(cgo, color1);
                      float ladder_alpha =
                          1.0f - AtomSettingGetWD(G, ai_i[i],
                                                  cSetting_cartoon_transparency,
                                                  1.0f - alpha);
                      cgo->add<cgo::draw::shadercylinder2ndcolor>(
                          cgo, cs->coordPtr(sug), axis, ladder_radius, 0x1f,
                          color2, &pickcolor2, ladder_alpha);
                    }
                  }
                }
              }
              base_at = save_at;
              sugar_at = save_at;
            }
          }
          if ((base_at >= 0) && (sugar_at >= 0)) {
            const AtomInfoType* sug_ai = atomInfo + sugar_at;
            const AtomInfoType* bas_ai = atomInfo + base_at;

            if (ring_connector_visible(G, bas_ai, sug_ai, sc_helper)) {
              int sug = cs->atmToIdx(sugar_at);
              int bas = cs->atmToIdx(base_at);

              if ((sug >= 0) && (bas >= 0)) {
                float tmp[3], outer[3];
                const float* v_outer = cs->coordPtr(sug);

                if ((o3_at >= 0) && (phos3_at < 0)) phos3_at = o3_at;
                if ((o5_at >= 0) && (phos5_at < 0)) phos5_at = o5_at;
                if ((ndata->na_mode != 1) && (phos3_at >= 0) &&
                    (phos5_at >= 0)) {
                  int p3 = cs->atmToIdx(phos3_at);
                  int p5 = cs->atmToIdx(phos5_at);
                  if ((p3 >= 0) && (p5 >= 0)) {
                    if (ring_mode) {
                      scale3f(cs->coordPtr(p5), 0.333333F, outer);
                      scale3f(cs->coordPtr(p3), 0.666667F, tmp);
                    } else {
                      scale3f(cs->coordPtr(p3), 0.5F, outer);
                      scale3f(cs->coordPtr(p5), 0.5F, tmp);
                    }
                    add3f(tmp, outer, outer);
                    v_outer = outer;
                  }
                }
                {
                  const float *color1, *color2;
                  if (ladder_color >= 0) {
                    color1 = color2 = ColorGet(G, ladder_color);
                  } else {
                    color1 = ColorGet(G, sug_ai->color);
                    color2 = ColorGet(G, bas_ai->color);
                  }

                  CGOPickColor(
                      cgo, sugar_at,
                      sug_ai->masked ? cPickableNoPick : cPickableAtom);
                  Pickable pickcolor2 = {
                      static_cast<unsigned int>(base_at),
                      bas_ai->masked ? cPickableNoPick : cPickableAtom};
                  float axis[3];
                  subtract3f(cs->coordPtr(bas), v_outer, axis);
                  CGOColorv(cgo, color1);
                  float ladder_alpha =
                      1.0f - AtomSettingGetWD(G, sug_ai,
                                              cSetting_cartoon_transparency,
                                              1.0f - alpha);
                  cgo->add<cgo::draw::shadercylinder2ndcolor>(
                      cgo, v_outer, axis, ladder_radius, 0x1f, color2,
                      &pickcolor2, ladder_alpha);
                }
              }
            }
          }
        }
      }
    }
    if ((!ring_mode) || (finder == 2)) {
      if (ladder_mode) { /* mark sure all rings traversed are mark */
        int i;
        for (i = 0; i <= n_atom; i++) {
          marked[atix[i]] = true;
        }
      }
    }
    if ((!nf) && ((have_C4_prime >= 0) || (have_C4 >= 0))) {
      int mem0;
      /* see if any of the neighbors are confirmed nucleic acids... */
      if (have_C4_prime >= 0)
        mem0 = have_C4_prime;
      else if (have_C4 >= 0)
        mem0 = have_C4;
      else
        mem0 = -1;
      if (mem0 >= 1) {
        auto seen = std::set<int>();
        nf = has_nuc_neighbor(nuc_flag, obj, mem0, 9, seen);
      }
    }
    if (n_atom) { /* store center of ring */

      float avg[3];
      float avg_col[3];
      int i;
      float up[3], upi[3];
      float vc0[3], vc1[3];
      const float* color = nullptr;
      /* compute average coordinate and mark atoms so that ring is only drawn
       * once */
      zero3f(avg);
      zero3f(avg_col);

      for (i = 0; i < n_atom; i++) {
        add3f(avg, v_i[i], avg);
        add3f(avg_col, col[i], avg_col);
        marked[atix[i]] = true;
      }

      scale3f(avg, 1.0F / n_atom, avg);
      scale3f(avg_col, 1.0F / n_atom, avg_col);

      for (i = 0; i < n_atom; i++) {
        float* v_moved = moved + atix[i] * 3;
        copy3f(avg, v_moved); /* store ring center for later use */
      }

      if ((nf || (!ladder_mode) || (finder >= 3)) && ring_mode &&
          (((finder == 1) && ((have_C4 >= 0) || (have_C4_prime >= 0))) ||
           ((finder == 2) && ((have_C4 >= 0))) ||
           ((finder == 3) && ((have_C_number >= 0))) || ((finder == 4)))) {
        auto atom_alpha =
            1.0f - AtomSettingGetWD(G, ai_i[i], cSetting_cartoon_transparency,
                                    1.0f - ring_alpha);
        if ((alpha != 1.0F) || (ring_alpha != alpha) || atom_alpha != 1.0) {
          if (atom_alpha != ring_alpha) {
            CGOAlpha(cgo, atom_alpha);
          } else {
            CGOAlpha(cgo, ring_alpha);
          }
        }

        if (ring_color >= 0) {
          color = ColorGet(G, ring_color);
        } else {
          color = avg_col;
        }

        CGOColorv(cgo, color);

        if ((ring_mode == 4) || (ring_mode == 5)) { /* spherical ring  */
          float radius = ring_radius;
          if (radius < 0.0F) {
            radius = 0.0F;
            if (ring_mode == 4) {
              for (i = 0; i < n_atom; i++) {
                float dist = diff3f(avg, v_i[i]);
                if (radius < dist) radius = dist;
              }
            } else {
              radius = ladder_radius;
            }
          }

          if (n_atom) {
            CGOColorv(cgo, avg_col);
            CGOPickColor(cgo, atix[0],
                         ai_i[0]->masked ? cPickableNoPick : cPickableAtom);
            CGOSphere(cgo, avg, radius);
          }

        } else {
          /* clear the normals */

          for (i = 0; i <= n_atom; i++) {
            zero3f(n_up[i]);
            zero3f(n_dn[i]);
          }

          /* compute average normals */

          {
            float acc[3];
            int ii;

            zero3f(acc);
            for (i = 0; i < n_atom; i++) {
              ii = i + 1;
              subtract3f(v_i[i], avg, vc0);
              subtract3f(v_i[ii], avg, vc1);
              cross_product3f(vc0, vc1, up);
              add3f(up, n_up[i], n_up[i]);
              add3f(up, n_up[ii], n_dn[ii]);
              if (!i) {
                add3f(up, n_up[n_atom], n_up[n_atom]);
              } else if (ii == n_atom) {
                add3f(up, n_up[0], n_up[0]);
              }
              add3f(up, acc, acc);
            }
            normalize3f(up);
            scale3f(up, -1.0F, upi);
          }

          for (i = 0; i <= n_atom; i++) {
            normalize3f(n_up[i]);
            scale3f(n_up[i], -1.0F, n_dn[i]);
          }

          {
            int ii;
            float mid[3];
            float up_add[3];
            float ct[3], cb[3];
            float v0t[3], v0b[3];
            float v1t[3], v1b[3];
            float out[3];

            CGOBegin(cgo, GL_TRIANGLES);
            for (i = 0; i < n_atom; i++) {
              ii = i + 1;
              average3f(v_i[ii], v_i[i], mid);

              subtract3f(mid, avg, out); /* compute outward-facing normal */
              normalize3f(out);

              scale3f(up, ring_width, up_add);

              add3f(avg, up_add, ct);
              subtract3f(avg, up_add, cb);

              add3f(v_i[i], up_add, v0t);
              subtract3f(v_i[i], up_add, v0b);

              add3f(v_i[ii], up_add, v1t);
              subtract3f(v_i[ii], up_add, v1b);

              CGONormalv(cgo, up);
              if (ring_color < 0) CGOColorv(cgo, color);
              CGOPickColor(cgo, atix[i],
                           ai_i[i]->masked ? cPickableNoPick : cPickableAtom);
              CGOVertexv(cgo, ct);
              CGONormalv(cgo, n_up[i]);
              if (ring_color < 0) CGOColorv(cgo, col[i]);
              //              CGOPickColor(cgo, atix[i], cPickableAtom);
              CGOVertexv(cgo, v0t);
              CGONormalv(cgo, n_up[ii]);
              if (ring_color < 0) CGOColorv(cgo, col[ii]);
              CGOPickColor(cgo, atix[ii],
                           ai_i[ii]->masked ? cPickableNoPick : cPickableAtom);
              //              CGOPickColor(cgo, atix[ii], cPickableAtom);
              CGOVertexv(cgo, v1t);

              if (ring_mode > 1) {
                CGONormalv(cgo, out);

                if (ring_color < 0) CGOColorv(cgo, col[i]);
                //                CGOPickColor(cgo, atix[i], cPickableAtom);
                CGOVertexv(cgo, v0t);
                CGOVertexv(cgo, v0b);
                if (ring_color < 0) CGOColorv(cgo, col[ii]);
                CGOVertexv(cgo, v1t);
                CGOVertexv(cgo, v1t);
                if (ring_color < 0) CGOColorv(cgo, col[i]);
                CGOVertexv(cgo, v0b);
                if (ring_color < 0) CGOColorv(cgo, col[ii]);
                CGOVertexv(cgo, v1b);
              }

              CGONormalv(cgo, upi);
              if (ring_color < 0) CGOColorv(cgo, color);
              CGOVertexv(cgo, cb);
              CGONormalv(cgo, n_dn[ii]);
              if (ring_color < 0) CGOColorv(cgo, col[ii]);
              CGOVertexv(cgo, v1b);
              CGONormalv(cgo, n_dn[i]);
              if (ring_color < 0) CGOColorv(cgo, col[i]);
              CGOVertexv(cgo, v0b);
            }
            CGOEnd(cgo);

            if ((alpha != 1.0F) || (ring_alpha != alpha)) CGOAlpha(cgo, alpha);

            {
              float ring_width_for_mode = ring_width;
              if (ring_mode == 3) {
                ring_width_for_mode = 3.f * ring_width;
              }
              for (i = 0; i < n_atom; i++) {
                ii = i + 1;
                {
                  const float *color1, *color2;
                  if (ring_color < 0) {
                    color1 = col[i];
                    color2 = col[ii];
                  } else {
                    color1 = color2 = color;
                  }
                  CGOPickColor(
                      cgo, atix[i],
                      ai_i[i]->masked ? cPickableNoPick : cPickableAtom);
                  Pickable pickcolor2 = {
                      static_cast<unsigned int>(atix[ii]),
                      ai_i[ii]->masked ? cPickableNoPick : cPickableAtom};
                  float axis[3];
                  subtract3f(v_i[ii], v_i[i], axis);
                  CGOColorv(cgo, color1);
                  float ladder_alpha =
                      1.0f - AtomSettingGetWD(G, ai_i[i],
                                              cSetting_cartoon_transparency,
                                              1.0f - alpha);
                  cgo->add<cgo::draw::shadercylinder2ndcolor>(
                      cgo, v_i[i], axis, ring_width_for_mode, 0x1f, color2,
                      &pickcolor2, ladder_alpha);
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
 * for nucleic acid polymers, fill "ndata" with:
 * - cartoon type and secondary structure
 * - check if single nucleotide or actual polymer
 */
void nuc_acid(PyMOLGlobals* G, nuc_acid_data* ndata, int a, int a1,
              const AtomInfoType* ai, const CoordSet* cs,
              const ObjectMolecule* obj, int set_flags) {
  int a3, a4, st, nd;
  const float *v_o, *v_c;
  const float* v1;
  int cur_car;
  const auto& nuc_flag = ndata->nuc_flag;

  if (ndata->a2 < 0) {
    ndata->nSeg++;
    ndata->v_o_last = nullptr;
  }
  *(ndata->sptr++) = ndata->nSeg;
  *(ndata->iptr++) = a;
  cur_car = ai->cartoon;
  if (cur_car == cCartoon_auto) cur_car = cCartoon_tube;
  (*ndata->ss) = ss_t::NUCLEIC; /* DNA/RNA */

  if (cur_car == cCartoon_putty) ndata->putty_flag = true;

  *(ndata->cc++) = cur_car;
  v1 = cs->coordPtr(a);
  copy3f(v1, ndata->vptr);
  ndata->vptr += 3;

  if (ndata->a2 >= 0) {
    if (set_flags) {
      if (AtomInfoIsNucBackboneTrace(LexStr(G, obj->AtomInfo[ndata->a2].name),
                                     obj->AtomInfo[ndata->a2].protons) &&
          (!nuc_flag[ndata->a2])) {
        int* nf = nullptr;
        AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, ndata->a2, &st,
                                   &nd);

        nf = nuc_flag + st;
        for (a3 = st; a3 <= nd; a3++) {
          *(nf++) = true;
        }
      }
    } else if ((ndata->na_mode >= 2) &&
               (!nuc_flag[ndata->a2])) { /* just a single nucleotide -- skip */
      cur_car = cCartoon_skip;
      *(ndata->cc - 2) = cCartoon_skip;
      *(ndata->cc - 1) = cCartoon_skip;
    }
  }

  ndata->a2 = a1;

  ndata->ss++;

  v_c = nullptr;
  v_o = nullptr;

  AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, a1, &st, &nd);

  {
    int* nf = nullptr;
    if (set_flags && ndata->v_o_last) nf = nuc_flag + st;
    for (a3 = st; a3 <= nd; a3++) {
      if (nf)
        *(nf++) =
            true; /* mark this residue as being part of a nucleic acid chain */
      a4 = cs->atmToIdx(a3);
      if (a4 >= 0) {
        if (ndata->na_mode == 1) {
          if (WordMatchExact(G, NUCLEIC_NORMAL1,
                             LexStr(G, obj->AtomInfo[a3].name), 1) ||
              WordMatchExact(G, NUCLEIC_NORMAL2,
                             LexStr(G, obj->AtomInfo[a3].name), 1)) {
            v_c = cs->coordPtr(a4);
          }
        } else if (a3 == a1) {
          v_c = cs->coordPtr(a4);
        }
        if (WordMatchExact(G, NUCLEIC_NORMAL0,
                           LexStr(G, obj->AtomInfo[a3].name), 1)) {
          v_o = cs->coordPtr(a4);
        }
      }
    }
  }
  if (!(v_c && v_o)) {
    zero3f(ndata->voptr);
    ndata->v_o_last = nullptr;
  } else {
    if (ndata->v_o_last) {
      float t0[3];
      add3f(v_o, ndata->v_o_last, t0);
      add3f(ndata->v_o_last, t0, t0);
      scale3f(t0, 0.333333F, t0);
      subtract3f(v_c, t0, ndata->voptr);
    } else {
      subtract3f(v_c, v_o, ndata->voptr);
    }
    ndata->v_o_last = v_o;
    normalize3f(ndata->voptr);
  }
  ndata->voptr += 3;
  ndata->nAt++;
  return;
}

int GenerateRepCartoonDrawRings(PyMOLGlobals* G, nuc_acid_data* ndata,
                                ObjectMolecule* obj, CoordSet* cs, CGO* cgo,
                                float ring_width, int cartoon_color,
                                float alpha) {
  int ring_i;
  int mem[8];
  int nbr[7];
  int* marked = pymol::calloc<int>(obj->NAtom);
  float* moved = pymol::calloc<float>(obj->NAtom * 3);
  int ring_color;
  int ok = true;
  int escape_count;
  int ladder_mode, ladder_color;
  float ladder_radius, ring_radius;
  int cartoon_side_chain_helper;
  float ring_alpha;
  ring_alpha = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_cartoon_ring_transparency);
  if (ring_alpha < 0.0F)
    ring_alpha = alpha;
  else
    ring_alpha = 1.0F - ring_alpha;
  cartoon_side_chain_helper =
      SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                   cSetting_cartoon_side_chain_helper);
  ladder_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_ladder_mode);
  ladder_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_cartoon_ladder_radius);
  ladder_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_cartoon_ladder_color);
  ring_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                             cSetting_cartoon_ring_radius);
  if (ladder_color == -1) ladder_color = cartoon_color;
  ring_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_cartoon_ring_color);
  if (ring_color == -1) ring_color = cartoon_color;

  const int* atmToIdx = obj->DiscreteFlag ? nullptr : cs->AtmToIdx.data();
  const auto* const neighbor = obj->getNeighborArray();

  escape_count = ESCAPE_MAX; /* don't get bogged down with structures
                                that have unreasonable connectivity */
  for (ring_i = 0; ok && ring_i < ndata->n_ring; ring_i++) {
    mem[0] = ndata->ring_anchor[ring_i];
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
                        /*    printf(" 5: %s(%d) %s(%d) %s(%d) %s(%d) %s(%d)\n",
                              obj->AtomInfo[mem[0]].name,mem[0],
                              obj->AtomInfo[mem[1]].name,mem[1],
                              obj->AtomInfo[mem[2]].name,mem[2],
                              obj->AtomInfo[mem[3]].name,mem[3],
                              obj->AtomInfo[mem[4]].name,mem[4]); */
                        do_ring(G, ndata, 5, mem, obj, cs, ring_width, cgo,
                                ring_color, ladder_radius, ladder_color,
                                ladder_mode, cartoon_side_chain_helper,
                                ring_alpha, alpha, marked, moved, ring_radius);
                      }
                      nbr[5] = neighbor[mem[5]] + 1;
                      while (((mem[6] = neighbor[nbr[5]]) >= 0) &&
                             ((!atmToIdx) || (atmToIdx[mem[5]] >= 0))) {
                        if ((mem[6] != mem[4]) && (mem[6] != mem[3]) &&
                            (mem[6] != mem[2]) && (mem[6] != mem[1])) {
                          if (mem[6] == mem[0]) { /* six-cycle */
                            /* printf(" 6: %s %s %s %s %s %s\n",
                               obj->AtomInfo[mem[0]].name,
                               obj->AtomInfo[mem[1]].name,
                               obj->AtomInfo[mem[2]].name,
                               obj->AtomInfo[mem[3]].name,
                               obj->AtomInfo[mem[4]].name,
                               obj->AtomInfo[mem[5]].name
                               ); */
                            do_ring(G, ndata, 6, mem, obj, cs, ring_width, cgo,
                                    ring_color, ladder_radius, ladder_color,
                                    ladder_mode, cartoon_side_chain_helper,
                                    ring_alpha, alpha, marked, moved,
                                    ring_radius);
                          }
                          nbr[6] = neighbor[mem[6]] + 1;
                          while (((mem[7] = neighbor[nbr[6]]) >= 0) &&
                                 ((!atmToIdx) || (atmToIdx[mem[6]] >= 0))) {
                            if ((mem[7] != mem[5]) && (mem[7] != mem[4]) &&
                                (mem[7] != mem[3]) && (mem[7] != mem[2]) &&
                                (mem[7] != mem[1])) {
                              if (mem[7] == mem[0]) {
                                do_ring(G, ndata, 7, mem, obj, cs, ring_width,
                                        cgo, ring_color, ladder_radius,
                                        ladder_color, ladder_mode,
                                        cartoon_side_chain_helper, ring_alpha,
                                        alpha, marked, moved, ring_radius);
                              }
                            }
                            nbr[6] += 2;
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
  }
  FreeP(marked);
  FreeP(moved);
  return ok;
}
