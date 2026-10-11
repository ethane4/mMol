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

// Assigning secondary structure (dss).

#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/select/selector_internal.h"

#define cSSMaxHBond 6

#define cSSHelix3HBond 0x0001
#define cSSHelix4HBond 0x0002
#define cSSHelix5HBond 0x0004
#define cSSGotPhiPsi 0x0008
#define cSSPhiPsiHelix 0x0010
#define cSSPhiPsiNotHelix 0x0020
#define cSSPhiPsiStrand 0x0040
#define cSSPhiPsiNotStrand 0x0080
#define cSSAntiStrandSingleHB 0x0100
#define cSSAntiStrandDoubleHB 0x0200
#define cSSAntiStrandBuldgeHB 0x0400
#define cSSAntiStrandSkip 0x0800
#define cSSParaStrandSingleHB 0x1000
#define cSSParaStrandDoubleHB 0x2000
#define cSSParaStrandSkip 0x4000

#define cSSBreakSize 5

typedef struct {
  int real;
  int ca, n, c, o; /* indices in selection-table space */
  float phi, psi;
  char ss, ss_save;
  int flags;
  int n_acc, n_don;
  int acc[cSSMaxHBond]; /* interactions where this residue is an acceptor */
  int don[cSSMaxHBond]; /* interactions where this residue is a donor */
  ObjectMolecule* obj;
  int preserve;
  int present;
} SSResi;

int SelectorAssignSS(PyMOLGlobals* G, int target, int present, int state_value,
                     int preserve, ObjectMolecule* single_object, int quiet) {
  /* PyMOL's secondary structure assignment algorithm:

     General principal -- if it looks like a duck, then it's a duck:

     I. Helices
     - must have reasonably helical geometry within the helical span
     - near-ideal geometry guarantees helix assignment
     - a continuous ladder stre i+3, i+4, or i+5 hydrogen bonding
     with permissible geometry can reinforce marginal cases
     - a minimum helix is three residues with i+3 H-bond

     II. Sheets
     - Hydrogen bonding ladders are the primary guide
     - Out-of-the envelope
     - 1-residue gaps in sheets are filled unless there
     is a turn.
   */

  CSelector* I = G->Selector;
  SSResi* res;
  int n_res = 0;
  int state_start, state_stop, state;
  int consensus = true;
  int first_last_only = false;
  int first_pass = true;

  if (!single_object) {
    if (state_value < 0) {
      switch (state_value) {
        case cSelectorUpdateTableCurrentState:
        case cSelectorUpdateTableEffectiveStates:
          SelectorUpdateTable(G, state_value, -1);
          break;
        default:
          SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
          break;
      }
    } else {
      SelectorUpdateTable(G, state_value, -1);
    }
  } else {
    SelectorUpdateTableSingleObject(G, single_object, state_value);
  }

  res = VLACalloc(SSResi, 1000);

  if (state_value < 0) {
    if (state_value == -4) consensus = false;
    if (state_value == -5) first_last_only = true;
    state_start = 0;
    state_stop = SelectorGetSeleNCSet(G, target);

    if (state_value == cStateCurrent) {
      StateIterator iter(G, nullptr, state_value, state_stop);
      if (iter.next()) {
        state_start = iter.state;
        state_stop = iter.state + 1;
      }
    }
  } else {
    state_start = state_value;
    state_stop = state_value + 1;
  }
  for (state = state_start; state < state_stop; state++) {
    int a;
    ObjectMolecule* obj;
    int aa, a0, a1, at, idx;
    AtomInfoType *ai, *ai0, *ai1;
    CoordSet* cs;
    ObjectMolecule* last_obj = nullptr;
    /* first, we need to count the number of residues under consideration */

    if (first_pass) {
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        obj = I->Obj[I->Table[a].model];
        at = +I->Table[a].atom;
        ai = obj->AtomInfo + at;

        /* see if CA coordinates exist... */

        if (SelectorIsMember(G, ai->selEntry, present)) {
          if ((ai->protons == cAN_C) &&
              (WordMatchExact(G, G->lex_const.CA, ai->name, true))) {
            if (last_obj != obj) {
              ObjectMoleculeVerifyChemistry(obj, state_value);
              last_obj = obj;
            }
            /* delimit residue */

            a0 = a - 1;
            while (a0 >= cNDummyAtoms) {
              ai0 = I->Obj[I->Table[a0].model]->AtomInfo + I->Table[a0].atom;
              if (!AtomInfoSameResidue(G, ai0, ai)) break;
              a0--;
            }

            a1 = a + 1;
            while (a1 < I->Table.size()) {
              ai1 = I->Obj[I->Table[a1].model]->AtomInfo + I->Table[a1].atom;
              if (!AtomInfoSameResidue(G, ai1, ai)) break;
              a1++;
            }

            {
              int found_N = 0;
              int found_O = 0;
              int found_C = 0;

              /* locate key atoms */

              for (aa = a0 + 1; aa < a1; aa++) {
                ai = I->Obj[I->Table[aa].model]->AtomInfo + I->Table[aa].atom;
                if ((ai->protons == cAN_C) &&
                    (WordMatchExact(G, G->lex_const.C, ai->name, true))) {
                  found_C = aa;
                }
                if ((ai->protons == cAN_N) &&
                    (WordMatchExact(G, G->lex_const.N, ai->name, true))) {
                  found_N = aa;
                }
                if ((ai->protons == cAN_O) &&
                    (WordMatchExact(G, G->lex_const.O, ai->name, true))) {
                  found_O = aa;
                }
              }

              if ((found_C) && (found_N) && (found_O)) {
                VLACheck(res, SSResi, n_res);
                res[n_res].n = found_N;
                res[n_res].o = found_O;
                res[n_res].c = found_C;
                res[n_res].ca = a;
                res[n_res].obj = I->Obj[I->Table[a].model];
                res[n_res].real = true;

                n_res++;

              } else {
                if (!quiet) {
                  PRINTFB(G, FB_Selector, FB_Warnings)
                  " AssignSS-Warning: Ignoring incomplete residue "
                  "/%s/%s/%s/%d%c ...\n",
                      obj->Name, LexStr(G, ai->segi), LexStr(G, ai->chain),
                      ai->resv, ai->getInscode(true) ENDFB(G);
                }
              }
            }
          }
        }
      } /* count pass */

      if (preserve) { /* if we're in preserve mode, then mark which objects
                         don't get changed */
        int a, b;
        char ss;
        ObjectMolecule* p_obj = nullptr;
        SSResi *r, *r2;
        for (a = 0; a < n_res; a++) {
          r = res + a;
          if (r->real) {
            if (p_obj != r->obj) {
              ss = r->obj->AtomInfo[I->Table[r->ca].atom].ssType[0];
              if ((ss == 'S') || (ss == 'H') || (ss == 's') || (ss == 'h')) {
                p_obj = r->obj;

                b = a;
                while (b >= 0) {
                  r2 = res + b;
                  if (p_obj == r2->obj) r2->preserve = true;
                  b--;
                }
                b = a + 1;
                while (b < n_res) {
                  r2 = res + b;
                  if (p_obj == r2->obj) r2->preserve = true;
                  b++;
                }
              }
            }
          }
        }
      }
      /*  printf("n_res %d\n",n_res); */

      /* now, let's repack res. into discrete chunks so that we can do easy gap
       * & ladder analysis */

      {
        SSResi* res2;
        int a;
        int n_res2 = 0;
        int add_break;
        int at_ca0, at_ca1;

        res2 = VLACalloc(SSResi, n_res * 2);

        for (a = 0; a < n_res; a++) {
          add_break = false;

          if (!a) {
            add_break = true;
          } else if (res[a].obj != res[a - 1].obj) {
            add_break = true;
          } else if (res[a].obj) {
            at_ca0 = I->Table[res[a].ca].atom;
            at_ca1 = I->Table[res[a - 1].ca].atom;
            if (!ObjectMoleculeCheckBondSep(res[a].obj, at_ca0, at_ca1,
                                            3)) { /* CA->N->C->CA = 3 bonds */
              add_break = true;
            }
          }

          if (add_break) {
            n_res2 += cSSBreakSize;
          }

          VLACheck(res2, SSResi, n_res2);
          res2[n_res2] = res[a];
          n_res2++;
        }

        n_res2 += cSSBreakSize;
        VLACheck(res2, SSResi, n_res2);

        VLAFreeP(res);
        res = res2;
        n_res = n_res2;
      }
      first_pass = false;
    }

    /* okay, the rest of this loop runs for each coordinate set */

    {
      int b;
      for (a = 0; a < n_res; a++) {
        res[a].present = res[a].real;

        if (res[a].present) {
          obj = res[a].obj;
          if (state < obj->NCSet)
            cs = obj->CSet[state];
          else
            cs = nullptr;
          for (b = 0; b < 4; b++) {
            if (cs) {
              switch (b) {
                case 0:
                  at = I->Table[res[a].n].atom;
                  break;
                case 1:
                  at = I->Table[res[a].o].atom;
                  break;
                case 2:
                  at = I->Table[res[a].c].atom;
                  break;
                default:
                case 3:
                  at = I->Table[res[a].ca].atom;
                  break;
              }
              idx = cs->atmToIdx(at);
            } else
              idx = -1;
            if (idx < 0) {
              res[a].present = false;
            }
          }
        }
      }
    }

    /* next, we need to record hydrogen bonding relationships */

    {
      float *v0, *v1;
      int n1;
      int at;

      int a, aa;
      int a0, a1;   /* SS res space */
      int at0, at1; /* object-atom space */
      int exclude;

      ObjectMolecule *obj0, *obj1;

      CoordSet* cs;
      float cutoff;
      HBondCriteria hbcRec, *hbc;
      int *zero = nullptr, *scratch = nullptr;

      {
        int max_n_atom = I->Table.size();
        ObjectMolecule* lastObj = nullptr;
        for (a = cNDummyAtoms; a < I->Table.size(); a++) {
          ObjectMolecule* obj = I->Obj[I->Table[a].model];
          if (obj != lastObj) {
            if (max_n_atom < obj->NAtom) max_n_atom = obj->NAtom;
            lastObj = obj;
          }
        }
        zero = pymol::calloc<int>(max_n_atom);
        scratch = pymol::malloc<int>(max_n_atom);
      }

      for (a = 0; a < n_res; a++) {
        res[a].n_acc = 0;
        res[a].n_don = 0;
      }
      hbc = &hbcRec;
      ObjectMoleculeInitHBondCriteria(G, hbc);

      /* use parameters which reflect the spirit of Kabsch and Sander
         ( i.e. long hydrogen-bonds/polar electrostatic interactions ) */

      hbc->maxAngle = 63.0F;
      hbc->maxDistAtMaxAngle = 3.2F;
      hbc->maxDistAtZero = 4.0F;
      hbc->power_a = 1.6F;
      hbc->power_b = 5.0F;
      hbc->cone_dangle = 0.0F; /* 180 deg. */
      if (hbc->maxDistAtMaxAngle != 0.0F) {
        hbc->factor_a = 0.5F / (float)pow(hbc->maxAngle, hbc->power_a);
        hbc->factor_b = 0.5F / (float)pow(hbc->maxAngle, hbc->power_b);
      }

      cutoff = hbc->maxDistAtMaxAngle;
      if (cutoff < hbc->maxDistAtZero) {
        cutoff = hbc->maxDistAtZero;
      }

      n1 = 0;

      const size_t table_size = I->Table.size();
      auto coords_flat = std::vector<float>(3 * table_size);
      auto* coords = pymol::reshape<3>(coords_flat.data());
      auto Flag1 = std::vector<MapFlag_t>(table_size, 0);
      auto Flag2 = std::vector<int>(table_size, 0);

      for (a = 0; a < n_res; a++) {
        if (res[a].present) {
          obj0 = res[a].obj;

          if (obj0) {
            /* map will contain the h-bond backbone nitrogens */

            aa = res[a].n;
            at = I->Table[aa].atom;
            Flag2[aa] = a; /* so we can find the atom again... */

            cs = obj0->getCoordSet(state);
            if (!cs) continue;

            if (CoordSetGetAtomVertex(cs, at, coords[aa])) {
              Flag1[aa] = true;
              n1++;
            }

            /* also copy O coordinates for usage below */

            aa = res[a].o;
            at = I->Table[aa].atom;
            CoordSetGetAtomVertex(cs, at, coords[aa]);
          }
        }
      }

      if (n1) {
        short too_many_atoms = false;
        std::unique_ptr<MapType> map(
            new MapType(G, -cutoff, pymol::flatten(coords), table_size, nullptr,
                        Flag1.data()));
        if (map) {
          for (a0 = 0; a0 < n_res; a0++) {
            if (res[a0].obj) {
              /* now iterate through carbonyls */
              obj0 = res[a0].obj;
              const auto as0 = res[a0].o;
              at0 = I->Table[as0].atom;

              v0 = coords[as0];

              int nat = 0;
              for (const auto as1 : MapEIter(*map, v0)) {
                v1 = coords[as1];

                if (within3f(v0, v1, cutoff)) {
                  obj1 = I->Obj[I->Table[as1].model];
                  at1 = I->Table[as1].atom;

                  if (obj0 ==
                      obj1) { /* don't count hbonds between adjacent residues */
                    exclude = SelectorCheckNeighbors(G, 5, obj0, at0, at1, zero,
                                                     scratch);
                  } else {
                    exclude = false;
                  }

                  /*                      if(!exclude) {
                     printf("at1 %s %s vs at0 %s %s\n",
                     obj1->AtomInfo[at1].resi,
                     obj1->AtomInfo[at1].name,
                     obj0->AtomInfo[at0].resi,
                     obj0->AtomInfo[at0].name
                     );
                     }
                   */
                  if ((!exclude) && ObjectMoleculeGetCheckHBond(
                                        NULL, nullptr, obj1, /* donor first */
                                        at1, state, obj0,    /* then acceptor */
                                        at0, state, hbc)) {
                    /*                        printf(" found hbond between
                       acceptor resi %s and donor resi %s\n",
                       res[a0].obj->AtomInfo[at0].resi,
                       res[I->Flag2[as1]].obj->AtomInfo[I->Table[as1].atom].resi);
                     */

                    a1 = Flag2[as1]; /* index in SS n_res space */

                    /* store acceptor link */

                    n1 = res[a0].n_acc;
                    if (n1 < (cSSMaxHBond - 1)) {
                      res[a0].acc[n1] = a1;
                      res[a0].n_acc = n1 + 1;
                    }

                    /* store donor link */

                    n1 = res[a1].n_don;
                    if (n1 < (cSSMaxHBond - 1)) {
                      res[a1].don[n1] = a0;
                      res[a1].n_don = n1 + 1;
                    }
                  }
                }
                nat++;
              }
              if (nat > 1000) {  // if map returns more than 1000 atoms within
                                 // 4, should be a dss error
                too_many_atoms = true;
                break;
              }
            }
          }
        }
        if (too_many_atoms) {
          PRINTFB(G, FB_Selector, FB_Errors)
          " %s: ERROR: Unreasonable number of neighbors for dss, cannot assign "
          "secondary structure.\n",
              __func__ ENDFB(G);
        }
      }
      FreeP(zero);
      FreeP(scratch);
    }

    { /* compute phi, psi's */

      SSResi* r;
      int a;

      float helix_psi_delta, helix_phi_delta;
      float strand_psi_delta, strand_phi_delta;

      float helix_psi_target =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_psi_target);
      float helix_psi_include =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_psi_include);
      float helix_psi_exclude =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_psi_exclude);

      float helix_phi_target =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_phi_target);
      float helix_phi_include =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_phi_include);
      float helix_phi_exclude =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_helix_phi_exclude);

      float strand_psi_target =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_psi_target);
      float strand_psi_include =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_psi_include);
      float strand_psi_exclude =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_psi_exclude);

      float strand_phi_target =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_phi_target);
      float strand_phi_include =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_phi_include);
      float strand_phi_exclude =
          SettingGet_f(G, nullptr, nullptr, cSetting_ss_strand_phi_exclude);

      for (a = 0; a < n_res; a++) {
        r = res + a;
        if (r->real && ((r - 1)->real)) {
          r->flags = 0;

          if (ObjectMoleculeGetPhiPsi(r->obj, I->Table[r->ca].atom, &r->phi,
                                      &r->psi, state)) {
            r->flags |= cSSGotPhiPsi;

            helix_psi_delta = (float)fabs(r->psi - helix_psi_target);
            strand_psi_delta = (float)fabs(r->psi - strand_psi_target);
            helix_phi_delta = (float)fabs(r->phi - helix_phi_target);
            strand_phi_delta = (float)fabs(r->phi - strand_phi_target);

            if (helix_psi_delta > 180.0F)
              helix_psi_delta = 360.0F - helix_psi_delta;
            if (strand_psi_delta > 180.0F)
              strand_psi_delta = 360.0F - strand_psi_delta;
            if (helix_phi_delta > 180.0F)
              helix_phi_delta = 360.0F - helix_phi_delta;
            if (strand_phi_delta > 180.0F)
              strand_phi_delta = 360.0F - strand_phi_delta;

            /* printf("helix %d strand %d\n",helix_delta,strand_delta); */

            if ((helix_psi_delta > helix_psi_exclude) ||
                (helix_phi_delta > helix_phi_exclude)) {
              r->flags |= cSSPhiPsiNotHelix;
            } else if ((helix_psi_delta < helix_psi_include) &&
                       (helix_phi_delta < helix_phi_include)) {
              r->flags |= cSSPhiPsiHelix;
            }

            if ((strand_psi_delta > strand_psi_exclude) ||
                (strand_phi_delta > strand_phi_exclude)) {
              r->flags |= cSSPhiPsiNotStrand;
            } else if ((strand_psi_delta < strand_psi_include) &&
                       (strand_phi_delta < strand_phi_include)) {
              r->flags |= cSSPhiPsiStrand;
            }
          }
        }
      }
    }

    /* by default, tentatively assign everything as loop */

    {
      int a;
      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        if (res[a].present) res[a].ss = 'L';
      }
    }

    {
      SSResi *r, *r2;
      int a, b, c;

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;
        if (r->real) {
          /* look for tell-tale i+3,4,5 hydrogen bonds for helix  */

          /* is residue an acceptor for i+3,4,5 residue? */
          for (b = 0; b < r->n_acc; b++) {
            r->flags |= ((r->acc[b] == (a + 3)) ? cSSHelix3HBond : 0) |
                        ((r->acc[b] == (a + 4)) ? cSSHelix4HBond : 0) |
                        ((r->acc[b] == (a + 5)) ? cSSHelix5HBond : 0);
          }

          /* is residue a donor for i-3,4,5 residue */
          for (b = 0; b < r->n_don; b++) {
            r->flags |= ((r->don[b] == (a - 3)) ? cSSHelix3HBond : 0) |
                        ((r->don[b] == (a - 4)) ? cSSHelix4HBond : 0) |
                        ((r->don[b] == (a - 5)) ? cSSHelix5HBond : 0);
          }

          /*        if(r->flags & (cSSHelix3HBond | cSSHelix4HBond |
             cSSHelix5HBond)) { printf("HelixHB %s \n",
             r->obj->AtomInfo[I->Table[r->ca].atom].resi);
             }
           */

          /* look for double h-bonded antiparallel beta sheet pairs:
           *
           *  \ /\ /
           *   N  C
           *   #  O
           *   O  #
           *   C  N
           *  / \/ \
           *
           */

          for (b = 0; b < r->n_acc; b++) { /* iterate through acceptors */
            r2 = (res + r->acc[b]);
            if (r2->real) {
              for (c = 0; c < r2->n_acc; c++) {
                if (r2->acc[c] == a) { /* found a pair */
                  r->flags |= cSSAntiStrandDoubleHB;
                  r2->flags |= cSSAntiStrandDoubleHB;

                  /*                printf("anti double %s to %s\n",
                     r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                     r2->obj->AtomInfo[I->Table[r2->ca].atom].resi); */
                }
              }
            }
          }

          /* look for antiparallel beta buldges
           *
           *     CCNC
           *  \ / O  \ /
           *   N      C
           *   #      O
           *    O    #
           *     C  N
           *    / \/ \
           *
           */

          for (b = 0; b < r->n_acc; b++) { /* iterate through acceptors */
            r2 = (res + r->acc[b]) + 1;    /* go forward 1 */
            if (r2->real) {
              for (c = 0; c < r2->n_acc; c++) {
                if (r2->acc[c] == a) { /* found a buldge */
                  r->flags |= cSSAntiStrandDoubleHB;
                  r2->flags |= cSSAntiStrandBuldgeHB;
                  (r2 - 1)->flags |= cSSAntiStrandBuldgeHB;

                  /*                printf("anti BULDGE %s to %s %s\n",
                     r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                     r2->obj->AtomInfo[I->Table[r2->ca].atom].resi,
                     r2->obj->AtomInfo[I->Table[(r2-1)->ca].atom].resi); */
                }
              }
            }
          }

          /* look for antiparallel beta sheet ladders (single or double)
           *
           *        O
           *     N  C
           *  \ / \/ \ /
           *   C      N
           *   O      #
           *   #      O
           *   N      C
           *  / \ /\ / \
           *     C  N
           *     O
           */

          if ((r + 1)->real && (r + 2)->real) {
            for (b = 0; b < r->n_acc; b++) { /* iterate through acceptors */
              r2 = (res + r->acc[b]) - 2;    /* go back 2 */
              if (r2->real) {
                for (c = 0; c < r2->n_acc; c++) {
                  if (r2->acc[c] == a + 2) { /* found a ladder */

                    (r)->flags |= cSSAntiStrandSingleHB;
                    (r + 1)->flags |= cSSAntiStrandSkip;
                    (r + 2)->flags |= cSSAntiStrandSingleHB;

                    (r2)->flags |= cSSAntiStrandSingleHB;
                    (r2 + 1)->flags |= cSSAntiStrandSkip;
                    (r2 + 2)->flags |= cSSAntiStrandSingleHB;

                    /*                  printf("anti ladder %s %s to %s %s\n",
                       r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                       r->obj->AtomInfo[I->Table[(r+2)->ca].atom].resi,
                       r2->obj->AtomInfo[I->Table[r2->ca].atom].resi,
                       r2->obj->AtomInfo[I->Table[(r2+2)->ca].atom].resi); */
                  }
                }
              }
            }
          }

          /* look for parallel beta sheet ladders
           *

           *    \ /\ /
           *     C  N
           *    O    #
           *   #      O
           *   N      C
           *  / \ /\ / \
           *     C  N
           *     O
           */

          if ((r + 1)->real && (r + 2)->real) {
            for (b = 0; b < r->n_acc; b++) { /* iterate through acceptors */
              r2 = (res + r->acc[b]);
              if (r2->real) {
                for (c = 0; c < r2->n_acc; c++) {
                  if (r2->acc[c] == a + 2) { /* found a ladder */

                    (r)->flags |= cSSParaStrandSingleHB;
                    (r + 1)->flags |= cSSParaStrandSkip;
                    (r + 2)->flags |= cSSParaStrandSingleHB;

                    (r2)->flags |= cSSParaStrandDoubleHB;

                    /*                                    printf("parallel
                       ladder %s %s to %s \n",
                       r->obj->AtomInfo[I->Table[r->ca].atom].resi,
                       r->obj->AtomInfo[I->Table[(r+2)->ca].atom].resi,
                       r2->obj->AtomInfo[I->Table[r2->ca].atom].resi); */
                  }
                }
              }
            }
          }
        }
      }
    }

    {
      int a;
      SSResi* r;
      /* convert flags to assignments */

      /* HELICES FIRST */

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;

        if (r->real) {
          /* clean internal helical residues are easy to find using H-bonds */

          if (((r - 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r + 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond))) {
            if (!(r->flags & (cSSPhiPsiNotHelix))) {
              r->ss = 'H';
            }
          }

          /*
             if(((r-1)->flags & (cSSHelix3HBond )) &&
             ((r  )->flags & (cSSHelix3HBond )) &&
             ((r+1)->flags & (cSSHelix3HBond ))) {
             if(!(r->flags & (cSSPhiPsiNotHelix))) {
             r->ss = 'H';
             }
             }

             if(((r-1)->flags & (cSSHelix4HBond)) &&
             ((r  )->flags & (cSSHelix4HBond)) &&
             ((r+1)->flags & (cSSHelix4HBond))) {
             if(!(r->flags & (cSSPhiPsiNotHelix))) {
             r->ss = 'H';
             }
             }

             if(((r-1)->flags & (cSSHelix5HBond)) &&
             ((r  )->flags & (cSSHelix5HBond)) &&
             ((r+1)->flags & (cSSHelix5HBond))) {
             if(!(r->flags & (cSSPhiPsiNotHelix))) {
             r->ss = 'H';
             }
             }
           */
        }
      }

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;

        if (r->real) {
          /* occasionally they'll be one whacked out residue missing h-bonds...
             in an otherwise good segment */

          if (((r - 2)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r - 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r - 1)->flags & (cSSPhiPsiHelix)) &&
              ((r)->flags & (cSSPhiPsiHelix)) &&
              ((r + 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r + 1)->flags & (cSSPhiPsiHelix)) &&
              ((r + 2)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond))) {
            r->ss = 'h';
          }
        }
      }

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;
        if (r->real) {
          if (r->ss == 'h') {
            r->flags |= (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond);
            r->ss = 'H';
          }
        }
      }

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;

        if (r->real) {
          /* deciding where the helix ends is trickier -- here we use helix
           * geometry */

          if (((r)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r)->flags & (cSSPhiPsiHelix)) &&
              ((r + 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r + 1)->flags & (cSSPhiPsiHelix)) &&
              ((r + 2)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r + 2)->flags & (cSSPhiPsiHelix)) && ((r + 1)->ss == 'H')) {
            r->ss = 'H';
          }

          if (((r)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r)->flags & (cSSPhiPsiHelix)) &&
              ((r - 1)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r - 1)->flags & (cSSPhiPsiHelix)) &&
              ((r - 2)->flags &
               (cSSHelix3HBond | cSSHelix4HBond | cSSHelix5HBond)) &&
              ((r - 2)->flags & (cSSPhiPsiHelix)) && ((r - 1)->ss == 'H')) {
            r->ss = 'H';
          }
        }
      }

      /* THEN SHEETS/STRANDS */

      for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
        r = res + a;
        if (r->real) {
          /* Antiparallel Sheets */

          if (((r)->flags & (cSSAntiStrandDoubleHB)) &&
              (!((r->flags & (cSSPhiPsiNotStrand))))) {
            (r)->ss = 'S';
          }

          if (((r)->flags &
               (cSSAntiStrandBuldgeHB)) && /* no strand geometry filtering for
                                              buldges.. */
              ((r + 1)->flags & (cSSAntiStrandBuldgeHB))) {
            (r)->ss = 'S';
            (r + 1)->ss = 'S';
          }

          if (((r - 1)->flags & (cSSAntiStrandDoubleHB)) &&
              ((r)->flags & (cSSAntiStrandSkip)) &&
              (!(((r)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r + 1)->flags &
               (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB))) {
            (r)->ss = 'S';
          }

          if (((r - 1)->flags &
               (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
              ((r)->flags & (cSSAntiStrandSkip)) &&
              (!(((r)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r + 1)->flags & (cSSAntiStrandDoubleHB))) {
            (r)->ss = 'S';
          }

          /* include open "ladders" if PHIPSI geometry supports assignment */

          if (((r - 1)->flags &
               (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
              ((r - 1)->flags & (cSSPhiPsiStrand)) &&
              (!(((r - 1)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r)->flags & (cSSPhiPsiStrand)) &&
              (!(((r - 1)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r + 1)->flags &
               (cSSAntiStrandSingleHB | cSSAntiStrandDoubleHB)) &&
              ((r + 1)->flags & (cSSPhiPsiStrand))) {
            (r - 1)->ss = 'S';
            (r)->ss = 'S';
            (r + 1)->ss = 'S';
          }

          /* Parallel Sheets */

          if (((r)->flags & (cSSParaStrandDoubleHB)) &&
              (!(((r)->flags & (cSSPhiPsiNotStrand))))) {
            (r)->ss = 'S';
          }

          if (((r - 1)->flags & (cSSParaStrandDoubleHB)) &&
              ((r)->flags & (cSSParaStrandSkip)) &&
              (!(((r)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r + 1)->flags &
               (cSSParaStrandSingleHB | cSSParaStrandDoubleHB))) {
            (r)->ss = 'S';
          }

          if (((r - 1)->flags &
               (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
              ((r)->flags & (cSSParaStrandSkip)) &&
              (!(((r)->flags & (cSSPhiPsiNotStrand)))) &&
              ((r + 1)->flags & (cSSParaStrandDoubleHB))) {
            (r)->ss = 'S';
          }

          /* include open "ladders" if PHIPSI geometry supports assignment */

          if (((r - 1)->flags &
               (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
              ((r - 1)->flags & (cSSPhiPsiStrand)) &&
              ((r)->flags & (cSSParaStrandSkip)) &&
              ((r)->flags & (cSSPhiPsiStrand)) &&
              ((r + 1)->flags &
               (cSSParaStrandSingleHB | cSSParaStrandDoubleHB)) &&
              ((r + 1)->flags & (cSSPhiPsiStrand))) {
            (r - 1)->ss = 'S';
            (r)->ss = 'S';
            (r + 1)->ss = 'S';
          }
        }
      }
    }

    {
      int a, b;
      SSResi *r, *r2;
      int repeat = true;
      int found;

      while (repeat) {
        repeat = false;

        for (a = cSSBreakSize; a < (n_res - cSSBreakSize); a++) {
          r = res + a;
          if (r->real) {
            /* make sure we don't have any 2-residue segments */

            if ((r->ss == 'S') && ((r + 1)->ss == 'S') &&
                (((r - 1)->ss != 'S') && ((r + 2)->ss != 'S'))) {
              r->ss = 'L';
              (r + 1)->ss = 'L';
              repeat = true;
            }
            if ((r->ss == 'H') && ((r + 1)->ss == 'H') &&
                (((r - 1)->ss != 'H') && ((r + 2)->ss != 'H'))) {
              r->ss = 'L';
              (r + 1)->ss = 'L';
              repeat = true;
            }

            /* make sure we don't have any 1-residue segments */

            if ((r->ss == 'S') &&
                (((r - 1)->ss != 'S') && ((r + 1)->ss != 'S'))) {
              r->ss = 'L';
              repeat = true;
            }
            if ((r->ss == 'H') &&
                (((r - 1)->ss != 'H') && ((r + 1)->ss != 'H'))) {
              r->ss = 'L';
              repeat = true;
            }

            /* double-check to make sure every terminal strand residue
               that should have a partner has one */

            if ((r->ss == 'S') &&
                (((r - 1)->ss != 'S') || ((r + 1)->ss != 'S'))) {
              found = false;

              for (b = 0; b < r->n_acc; b++) {
                r2 = res + r->acc[b];
                if (r2->ss == r->ss) {
                  found = true;
                  break;
                }
              }

              if (!found) {
                for (b = 0; b < r->n_don; b++) {
                  r2 = res + r->don[b];
                  if (r2->ss == r->ss) {
                    found = true;
                    break;
                  }
                }
              }

              if (!found) {
                /* allow these strand "skip" residues to persist if a neighbor
                 * has hydrogen bonds */
                if (r->flags & (cSSAntiStrandSkip | cSSParaStrandSkip)) {
                  if ((r + 1)->ss == r->ss)
                    for (b = 0; b < (r + 1)->n_acc; b++) {
                      r2 = res + (r + 1)->acc[b];
                      if (r2->ss == r->ss) {
                        found = true;
                        break;
                      }
                    }

                  if (!found) {
                    if ((r - 1)->ss == r->ss) {
                      for (b = 0; b < (r - 1)->n_don; b++) {
                        r2 = res + (r - 1)->don[b];
                        if (r2->ss == r->ss) {
                          found = true;
                          break;
                        }
                      }
                    }
                  }
                }
              }

              if (!found) {
                r->ss = 'L';
                repeat = true;
              }
            }
          }
        }
      }
    }

    {
      int a;
      for (a = 0; a < n_res;
           a++) { /* now apply consensus or union behavior, if appropriate */
        if (res[a].present) {
          if (res[a].ss_save) {
            if (res[a].ss != res[a].ss_save) {
              if (consensus) {
                res[a].ss = res[a].ss_save = 'L';
              } else if (res[a].ss == 'L')
                res[a].ss = res[a].ss_save;
            }
          }
          res[a].ss_save = res[a].ss;
        }
      }
    }

    {
      int a, aa;
      ObjectMolecule *obj = nullptr, *last_obj = nullptr;
      AtomInfoType* ai;
      int changed_flag = false;

      for (a = 0; a < n_res; a++) {
        if (res[a].present && (!res[a].preserve)) {
          aa = res[a].ca;
          obj = I->Obj[I->Table[aa].model];

          if (obj != last_obj) {
            if (changed_flag && last_obj) {
              last_obj->invalidate(cRepCartoon, cRepInvRep, -1);
              SceneChanged(G);
              changed_flag = false;
            }
            last_obj = obj;
          }
          ai = obj->AtomInfo + I->Table[aa].atom;

          if (SelectorIsMember(G, ai->selEntry, target)) {
            ai->ssType[0] = res[a].ss;
            ai->cartoon = 0; /* switch back to auto */
            ai->ssType[1] = 0;
            changed_flag = true;
          }
        }
      }

      if (changed_flag && last_obj) {
        last_obj->invalidate(cRepCartoon, cRepInvRep, -1);
        SceneChanged(G);
        changed_flag = false;
      }
    }
    if (first_last_only && (state == state_start)) state = state_stop - 2;
  }

  VLAFreeP(res);
  return 1;
}
