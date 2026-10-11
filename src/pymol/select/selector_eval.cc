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

// Evaluating selection operators on the atom table.

#include "pymol/chem/ring_finder.h"
#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/executive/editor.h"
#include "pymol/ui/seeker.h"
#include "pymol/base/util2.h"
#include "pymol/objects/property.h"
#include "pymol/select/selector_internal.h"

#define cDummyOrigin 0
#define cDummyCenter 1

static WordKeyValue rep_names[] = {
    {"spheres", cRepSphereBit},
    {"sticks", cRepCylBit},
    {"surface", cRepSurfaceBit},
    {"labels", cRepLabelBit},
    {"nb_spheres", cRepNonbondedSphereBit},
    {"cartoon", cRepCartoonBit},
    {"ribbon", cRepRibbonBit},
    {"lines", cRepLineBit},
    {"dots", cRepDotBit},
    {"mesh", cRepMeshBit},
    {"nonbonded", cRepNonbondedBit},
    {"ellipsoid", cRepEllipsoidBit},
    {
        "",
    },
};

static const char* backbone_names[] = {
    // protein
    "CA", "C", "O", "N", "OXT", "H",
    // nucleic acid
    "P", "OP1", "OP2", "OP3", "C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'",
    "C5'", "O5'", "H1'", "H3'", "H4'", "H2'", "H2''", "H12'", "H22'", "H5'",
    "H5''", "H15'", "H25'", "HO2'", "HO3'", "HO5'", ""};

#define SCMP_GTHN 0x01
#define SCMP_LTHN 0x02
#define SCMP_RANG 0x03
#define SCMP_EQAL 0x04
#define SCMP_GEQL 0x05
#define SCMP_LEQL 0x06

static WordKeyValue AtOper[] = {
    {">", SCMP_GTHN},  {"<", SCMP_LTHN},  {"in", SCMP_RANG}, {"=", SCMP_EQAL},
    {"==", SCMP_EQAL}, {">=", SCMP_GEQL}, {"<=", SCMP_LEQL}, {"", 0},
};

static short fcmp(float a, float b, int oper) {
  switch (oper) {
    case SCMP_GTHN:
      return (a > b);
    case SCMP_LTHN:
      return (a < b);
    case SCMP_EQAL:
      return fabs(a - b) < R_SMALL4;
    case SCMP_GEQL:
      return (a >= b);
    case SCMP_LEQL:
      return (a <= b);
  }
  printf("ERROR: invalid operator %d\n", oper);
  return false;
}
static int SelectorGetArrayNCSet(PyMOLGlobals* G, const sele_array_t& uptr,
                                 int no_dummies) {
  const int* array = uptr.get();
  CSelector* I = G->Selector;
  int a;
  ObjectMolecule* obj;
  int result = 0;
  int start = 0;
  if (no_dummies) start = cNDummyAtoms;
  for (a = start; a < I->Table.size(); a++) {
    if (*(array++)) {
      if (a >= cNDummyAtoms) {
        obj = I->Obj[I->Table[a].model];
        if (result < obj->NCSet) result = obj->NCSet;
      } else {
        if (result < 1) result = 1; /* selected dummy has at least one CSet */
      }
    }
  }
  return (result);
}
pymol::Result<sele_array_t> SelectorSelect(PyMOLGlobals* G, const char* sele,
                                           int state, SelectorID_t domain,
                                           int quiet) {
  SelectorUpdateTable(G, state, domain);
  auto parsed = SelectorParse(G, sele);
  if (!parsed.empty()) {
    return SelectorEvaluate(G, parsed, state, quiet);
  }
  return {};
}
int SelectorModulate1(PyMOLGlobals* G, EvalElem* base, int state) {
  CSelector* I = G->Selector;
  int a, d, e;
  int c = 0;
  float dist;
  int nbond;
  float* v2;
  CoordSet* cs;
  int ok = true;
  int nCSet;
  int n1, at, idx;
  ObjectMolecule* obj;

  if (state < 0) {
    switch (state) {
      case cSelectorUpdateTableCurrentState:
      case cSelectorUpdateTableEffectiveStates:
        state = SceneGetState(G);
        break;
    }
  }

  base[1].sele = std::move(base[0].sele); /* base1 has the mask */
  base->sele_calloc(I->Table.size());
  base->sele_check_ok(ok);
  if (!ok) return false;
  switch (base[1].code) {
    case SELE_ARD_:
    case SELE_EXP_:
      if (!sscanf(base[2].text(), "%f", &dist))
        ok = ErrMessage(G, "Selector", "Invalid distance.");
      if (ok) {
        for (d = 0; d < I->NCSet; d++) {
          if ((state < 0) || (d == state)) {
            n1 = 0;

            const size_t table_size = I->Table.size();
            auto coords_flat = std::vector<float>(table_size * 3);
            auto* coords = pymol::reshape<3>(coords_flat.data());
            auto Flag1 = std::vector<MapFlag_t>(table_size, 0);

            // Potential atoms to be selected (exclude dummies)
            for (a = cNDummyAtoms; a < I->Table.size(); a++) {
              at = I->Table[a].atom;
              obj = I->Obj[I->Table[a].model];
              if (d < obj->NCSet)
                cs = obj->CSet[d];
              else
                cs = nullptr;
              if (cs) {
                if (CoordSetGetAtomVertex(cs, at, coords[a])) {
                  Flag1[a] = true;
                  n1++;
                }
              }
            }
            if (n1) {
              std::unique_ptr<MapType> map(
                  new MapType(G, -dist, pymol::flatten(coords), table_size,
                              nullptr, Flag1.data()));
              CHECKOK(ok, map);
              if (ok) {
                nCSet = SelectorGetArrayNCSet(G, base[1].sele, false);
                for (e = 0; ok && e < nCSet; e++) {
                  if ((state < 0) || (e == state)) {
                    // Input selection (include dummies)
                    for (a = 0; ok && a < I->Table.size(); a++) {
                      if (base[1].sele[a]) {
                        at = I->Table[a].atom;
                        obj = I->Obj[I->Table[a].model];
                        if (e < obj->NCSet)
                          cs = obj->CSet[e];
                        else
                          cs = nullptr;
                        if (cs) {
                          idx = cs->atmToIdx(at);
                          if (idx >= 0) {
                            v2 = cs->coordPtr(idx);
                            for (const auto j : MapEIter(*map, v2, false)) {
                              if (!base[0].sele[j] &&
                                  (!base[1].sele[j] ||
                                   base[1].code == SELE_EXP_)) {
                                /*exclude current selection */
                                if (within3f(coords[j], v2, dist))
                                  base[0].sele[j] = true;
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
      break;

    case SELE_EXT_:
      if (sscanf(base[2].text(), "%d", &nbond) != 1)
        ok = ErrMessage(G, "Selector", "Invalid bond count.");
      if (ok) {
        std::copy_n(base[1].sele_data(), I->Table.size(), base[0].sele_data());
        while ((nbond--) > 0) {
          std::swap(base[1].sele, base[0].sele);
          for (unsigned a = cNDummyAtoms; a < I->Table.size(); a++) {
            if (base[1].sele[a]) {
              auto const* lastObj = I->Obj[I->Table[a].model];
              for (auto const& neighbor :
                   AtomNeighbors(lastObj, I->Table[a].atom)) {
                auto const a2 =
                    SelectorGetObjAtmOffset(I, lastObj, neighbor.atm);
                assert(a2 >= 0);
                if (a2 >= 0) {
                  base[0].sele[a2] = 1;
                }
              }
            }
          }
        }
        base[1].sele_free();
      }
      break;

    case SELE_GAP_:
      if (!sscanf(base[2].text(), "%f", &dist))
        ok = ErrMessage(G, "Selector", "Invalid distance.");
      if (ok) {
        for (a = 0; a < I->Table.size(); a++) {
          obj = I->Obj[I->Table[a].model];
          at = I->Table[a].atom;
          I->Table[a].f1 = obj->AtomInfo[at].vdw;
          base[0].sele[a] = true; /* start selected, subtract off */
          c = I->Table.size();
        }
        for (d = 0; d < I->NCSet; d++) {
          if ((state < 0) || (d == state)) {
            n1 = 0;

            auto Flag1 = std::vector<MapFlag_t>(I->Table.size(), 0);
            auto Vertex = std::vector<float>(I->Table.size() * 3, 0.0f);

            for (a = 0; a < I->Table.size(); a++) {
              obj = I->Obj[I->Table[a].model];
              at = I->Table[a].atom;
              if (d < obj->NCSet)
                cs = obj->CSet[d];
              else
                cs = nullptr;
              if (cs) {
                if (CoordSetGetAtomVertex(cs, at, Vertex.data() + 3 * a)) {
                  Flag1[a] = true;
                  n1++;
                }
              }
            }
            if (n1) {
              std::unique_ptr<MapType> map(
                  new MapType(G, -(dist + 2 * MAX_VDW), Vertex.data(),
                              I->Table.size(), nullptr, Flag1.data()));
              CHECKOK(ok, map);
              if (ok) {
                nCSet = SelectorGetArrayNCSet(G, base[1].sele, false);
                for (e = 0; ok && e < nCSet; e++) {
                  if ((state < 0) || (e == state)) {
                    for (a = 0; ok && a < I->Table.size(); a++) {
                      if (base[1].sele[a]) {
                        at = I->Table[a].atom;
                        obj = I->Obj[I->Table[a].model];
                        if (e < obj->NCSet)
                          cs = obj->CSet[e];
                        else
                          cs = nullptr;
                        if (cs) {
                          idx = cs->atmToIdx(at);

                          if (idx >= 0) {
                            v2 = cs->coordPtr(idx);
                            for (const auto j : MapEIter(*map, v2, false)) {
                              if ((base[0].sele[j]) &&
                                  (!base[1].sele[j])) { /*exclude current
                                                           selection */
                                if (within3f(
                                        Vertex.data() + 3 * j, v2,
                                        dist + /* eliminate atoms w/o gap */
                                            I->Table[a].f1 + I->Table[j].f1)) {
                                  base[0].sele[j] = false;
                                  c--;
                                }
                              } else if (base[1].sele[j]) {
                                base[0].sele[j] = false;
                                c--;
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
      break;
  }
  base[1].sele_free();
  if (Feedback(G, FB_Selector, FB_Debugging)) {
    c = 0;
    for (a = cNDummyAtoms; a < I->Table.size(); a++)
      if (base[0].sele[a]) c++;
    fprintf(stderr, "SelectorModulate0: %d atoms selected.\n", c);
  }
  return (ok);
}
int SelectorSelect0(PyMOLGlobals* G, EvalElem* passed_base) {
  CSelector* I = G->Selector;
  int a, b, flag;
  EvalElem* base = passed_base;
  int c = 0;
  ObjectMolecule *obj, *cur_obj = nullptr;
  CoordSet* cs;

  base->type = STYP_LIST;
  base->sele_calloc(I->Table.size());
  base->sele_err_chk_ptr(G);

  switch (base->code) {
    case SELE_HBAs:
    case SELE_HBDs:
    case SELE_DONz:
    case SELE_ACCz:
    case SELE_DESz: {
      /* first, verify chemistry for all atoms... */
      ObjectMolecule *lastObj = nullptr, *obj;
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        obj = I->Obj[I->Table[a].model];
        if (obj != lastObj) {
          ObjectMoleculeVerifyChemistry(obj, -1);
          lastObj = obj;
        }
      }
    }
      switch (base->code) {
        case SELE_HBAs:
        case SELE_ACCz:
          for (a = cNDummyAtoms; a < I->Table.size(); a++)
            base[0].sele[a] = I->Obj[I->Table[a].model]
                                  ->AtomInfo[I->Table[a].atom]
                                  .hb_acceptor;
          break;
        case SELE_HBDs:
        case SELE_DONz:
          for (a = cNDummyAtoms; a < I->Table.size(); a++)
            base[0].sele[a] =
                I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hb_donor;
          break;
        case SELE_DESz:
          for (a = cNDummyAtoms; a < I->Table.size(); a++) {
            const auto* m = I->Obj[I->Table[a].model];
            const auto at_i = I->Table[a].atom;
            const auto deloc = (float)getExplicitDegree(m, at_i) /
                               (float)getExplicitValence(m, at_i);
            base[0].sele[a] = floor(deloc) != deloc;
          }
          break;
      }
      break;
    case SELE_NONz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++) base[0].sele[a] = false;
      break;
    case SELE_BNDz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].bonded;
      break;
    case SELE_HETz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].hetatm;
      break;
    case SELE_HYDz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].isHydrogen();
      break;
    case SELE_METz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].isMetal();
      }
      break;
    case SELE_BB_z:
    case SELE_SC_z: {
      AtomInfoType* ai;
      flag = (base->code == SELE_BB_z);
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        ai = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
        if (!(ai->flags & cAtomFlag_polymer)) {
          base[0].sele[a] = 0;
          continue;
        }
        base[0].sele[a] = !flag;
        for (b = 0; backbone_names[b][0]; b++) {
          if (!(strcmp(LexStr(G, ai->name), backbone_names[b]))) {
            base[0].sele[a] = flag;
            break;
          }
        }
      }
    } break;
    case SELE_FXDz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_fix;
      break;
    case SELE_RSTz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_restrain;
      break;
    case SELE_POLz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_polymer;
      break;
    case SELE_PROz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_protein;
      break;
    case SELE_NUCz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_nucleic;
      break;
    case SELE_SOLz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_solvent;
      break;
    case SELE_PTDz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].protekted !=
            cAtomProtected_off;
      break;
    case SELE_MSKz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].masked;
      break;
    case SELE_ORGz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_organic;
      break;
    case SELE_INOz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            cAtomFlag_inorganic;
      break;
    case SELE_GIDz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++)
        base[0].sele[a] =
            bool(I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
                 cAtomFlag_guide);
      break;

    case SELE_PREz:
      cs = nullptr;
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        base[0].sele[a] = false;
        obj = I->Obj[I->Table[a].model];
        if (obj != cur_obj) { /* different object */
          cs = obj->getCoordSet(cSelectorUpdateTableCurrentState);
          cur_obj = obj;
        }
        if (cs) {
          if (cs->atmToIdx(I->Table[a].atom) >= 0) {
            base[0].sele[a] = true;
            c++;
          }
        }
      }
      break;
    case SELE_ALLz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        base[0].sele[a] = true;
        c++;
      }
      break;
    case SELE_ORIz:
      for (a = 0; a < I->Table.size(); a++) {
        base[0].sele[a] = false;
        c++;
      }
      if (I->Origin)
        ObjectMoleculeDummyUpdate(I->Origin.get(), cObjectMoleculeDummyOrigin);
      base[0].sele[cDummyOrigin] = true;
      break;
    case SELE_CENz:
      for (a = 0; a < I->Table.size(); a++) {
        base[0].sele[a] = false;
        c++;
      }
      if (I->Center)
        ObjectMoleculeDummyUpdate(I->Center.get(), cObjectMoleculeDummyCenter);
      base[0].sele[cDummyCenter] = true;
      break;
    case SELE_VISz: {
      ObjectMolecule* last_obj = nullptr;
      AtomInfoType* ai;
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        flag = false;
        obj = I->Obj[I->Table[a].model];
        if (obj->Enabled) {
          ai = obj->AtomInfo + I->Table[a].atom;

          if (last_obj != obj) {
            ObjectMoleculeVerifyChemistry(obj, -1);
            last_obj = obj;
          }

          flag = ai->isVisible();
        }
        base[0].sele[a] = flag;
        if (flag) c++;
      }
    } break;
    case SELE_ENAz:
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        flag = (I->Obj[I->Table[a].model]->Enabled);
        base[0].sele[a] = flag;
        if (flag) c++;
      }
      break;
  }
  PRINTFD(G, FB_Selector)
  " %s: %d atoms selected.\n", __func__, c ENDFD;

  return (1);
}
pymol::Result<> SelectorSelect1(PyMOLGlobals* G, EvalElem* base, int quiet) {
  CSelector* I = G->Selector;
  auto IM = I->mgr;
  CWordMatcher* matcher = nullptr;
  int a, b, c = 0, hit_flag;
  ObjectMolecule *obj, *last_obj;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  int ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);
  int I_NAtom = I->Table.size();
  int* base_0_sele_a;

  int model, s, col_idx;
  int flag;
  int index, state;
  int rep_mask;
  const char* wildcard = SettingGetGlobal_s(G, cSetting_wildcard);

  ObjectMolecule* cur_obj = nullptr;
  CoordSet* cs = nullptr;

  base->type = STYP_LIST;
  base->sele_calloc(I_NAtom); /* starting with zeros */
  base->sele_err_chk_ptr(G);
  switch (base->code) {
    case SELE_PEPs:
      if (base[1].text()[0]) {
        AtomInfoType *last_ai0 = nullptr, *ai0;
        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          ai0 = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
          if (!AtomInfoSameResidueP(G, ai0,
                                    last_ai0)) { /* new starting residue */
            int match_found = false;
            const char* ch = base[1].text(); /* sequence argument */
            AtomInfoType *ai1, *last_ai1 = nullptr;
            for (b = a; b < I_NAtom; b++) {
              ai1 = I->Obj[I->Table[b].model]->AtomInfo + I->Table[b].atom;
              if (!AtomInfoSameResidueP(G, ai1, last_ai1)) {
                if (*ch != '-') { /* if not skipping this residue */
                  if (!((*ch == '+') ||
                        (SeekerGetAbbr(G, LexStr(G, ai1->resn), 'O', 0) ==
                         *ch))) { /* if a mismatch */
                    break;
                  }
                }
                ch++;
                if (!*ch) { /* end of sequence pattern */
                  match_found = true;
                  break;
                }
                last_ai1 = ai1;
              }
            }
            if (match_found) {
              const char* ch = base[1].text(); /* sequence argument */
              AtomInfoType *ai1, *last_ai1 = nullptr, *ai2;
              for (b = a; b < I_NAtom; b++) {
                ai1 = I->Obj[I->Table[b].model]->AtomInfo + I->Table[b].atom;
                if (!AtomInfoSameResidueP(G, ai1, last_ai1)) {
                  if (*ch != '-') { /* if not skipping this residue */
                    if ((*ch == '+') ||
                        (SeekerGetAbbr(G, LexStr(G, ai1->resn), 'O', 0) ==
                         *ch)) { /* if matched */
                      int d;
                      for (d = b; d < I_NAtom; d++) {
                        ai2 = I->Obj[I->Table[d].model]->AtomInfo +
                              I->Table[d].atom; /* complete residue */
                        if (AtomInfoSameResidue(G, ai1, ai2)) {
                          c++;
                          base[0].sele[d] = true;
                        }
                      }
                    }
                  }
                  ch++;
                  if (!*ch) { /* end of sequence pattern */
                    break;
                  }
                  last_ai1 = ai1;
                }
              }
            }
          }
        }
      }
      break;
    case SELE_IDXs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigInteger(&options);

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a =
                   WordMatcherMatchInteger(matcher, table_a.atom + 1)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }

    } break;
    case SELE_ID_s: {
      CWordMatchOptions options;

      WordMatchOptionsConfigInteger(&options);

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchInteger(
                   matcher, I->Obj[table_a.model]->AtomInfo[table_a.atom].id)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_RNKs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigInteger(&options);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchInteger(
                   matcher,
                   I->Obj[table_a.model]->AtomInfo[table_a.atom].rank)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_NAMs: {
      CWordMatchOptions options;
      const char* atom_name_wildcard =
          SettingGetGlobal_s(G, cSetting_atom_name_wildcard);

      if (!atom_name_wildcard[0]) atom_name_wildcard = wildcard;

      WordMatchOptionsConfigAlphaList(&options, atom_name_wildcard[0],
                                      ignore_case);

      matcher = WordMatcherNew(G, base[1].text(), &options, false);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];
      last_obj = nullptr;
      for (a = cNDummyAtoms; a < I_NAtom; a++) {
        auto& table_a = I->Table[a];
        obj = I->Obj[table_a.model];
        if (obj != last_obj) {
          /* allow objects to have their own atom_name_wildcards...this is a
             tricky workaround for handling nucleic acid structures that use "*"
             in atom names */

          const char* atom_name_wildcard = SettingGet_s(
              G, obj->Setting.get(), nullptr, cSetting_atom_name_wildcard);

          if (!atom_name_wildcard[0]) atom_name_wildcard = wildcard;

          if (options.wildcard != atom_name_wildcard[0]) {
            options.wildcard = atom_name_wildcard[0];
            if (matcher) WordMatcherFree(matcher);
            matcher = WordMatcherNew(G, base[1].text(), &options, false);
            if (!matcher)
              WordPrimeCommaMatch(
                  G, &base[1].m_text[0] /* replace '+' with ',' */);
          }
          last_obj = obj;
        }

        const char* name = LexStr(G, obj->AtomInfo[table_a.atom].name);
        if (matcher)
          hit_flag = WordMatcherMatchAlpha(matcher, name);
        else
          hit_flag =
              (WordMatchCommaExact(G, base[1].text(), name, ignore_case) < 0);

        if ((*base_0_sele_a = hit_flag)) c++;
        base_0_sele_a++;
      }
      if (matcher) WordMatcherFree(matcher);
    } break;
    case SELE_TTYs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          auto ai = I->Obj[table_a.model]->AtomInfo + table_a.atom;
          if ((*base_0_sele_a =
                   WordMatcherMatchAlpha(matcher, LexStr(G, ai->textType))))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_ELEs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchAlpha(
                   matcher,
                   I->Obj[table_a.model]->AtomInfo[table_a.atom].elem)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_STRO: {
      CWordMatchOptions options;
      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];
        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          const char* mmstereotype = AtomInfoGetStereoAsStr(
              I->Obj[table_a.model]->AtomInfo + table_a.atom);
          if ((*base_0_sele_a = WordMatcherMatchAlpha(matcher, mmstereotype)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_REPs:
      rep_mask = 0;
      WordPrimeCommaMatch(G, &base[1].m_text[0] /* replace '+' with ',' */);
      for (a = 0; rep_names[a].word[0]; a++) {
        if (WordMatchComma(G, base[1].text(), rep_names[a].word, ignore_case) <
            0)
          rep_mask |= rep_names[a].value;
      }
      for (SelectorAtomIterator iter(I); iter.next();) {
        if (iter.getAtomInfo()->visRep & rep_mask) {
          base[0].sele[iter.a] = true;
          c++;
        } else {
          base[0].sele[iter.a] = false;
        }
      }
      break;
    case SELE_COLs:
      col_idx = ColorGetIndex(G, base[1].text());
      for (a = cNDummyAtoms; a < I_NAtom; a++) {
        base[0].sele[a] = false;
        if (I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].color ==
            col_idx) {
          base[0].sele[a] = true;
          c++;
        }
      }
      break;
    case SELE_CCLs:
    case SELE_RCLs:
      // setting index
      index = (base->code == SELE_CCLs) ? cSetting_cartoon_color
                                        : cSetting_ribbon_color;
      col_idx = ColorGetIndex(G, base[1].text());
      for (a = cNDummyAtoms; a < I_NAtom; a++) {
        base[0].sele[a] = false;
        {
          AtomInfoType* ai =
              I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
          int value;
          if (AtomSettingGetIfDefined(G, ai, index, &value)) {
            if (value == col_idx) {
              base[0].sele[a] = true;
              c++;
            }
          }
        }
      }
      break;
    case SELE_CHNs:
    case SELE_SEGs:
    case SELE_CUST:
    case SELE_LABs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case_chain);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      int offset = 0;
      switch (base->code) {
        case SELE_CHNs:
          offset = offsetof(AtomInfoType, chain);
          break;
        case SELE_SEGs:
          offset = offsetof(AtomInfoType, segi);
          break;
        case SELE_CUST:
          offset = offsetof(AtomInfoType, custom);
          break;
        case SELE_LABs:
          offset = offsetof(AtomInfoType, label);
          break;
        default:
          printf("coding error: missing case\n");
      }

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchAlpha(
                   matcher,
                   LexStr(G, *reinterpret_cast<decltype(AtomInfoType::chain)*>(
                                 ((char*)(I->Obj[table_a.model]->AtomInfo +
                                          table_a.atom)) +
                                 offset)))))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_SSTs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchAlpha(
                   matcher,
                   I->Obj[table_a.model]->AtomInfo[table_a.atom].ssType)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_STAs:
      sscanf(base[1].text(), "%d", &state);
      state = state - 1;
      obj = nullptr;

      if (state < 0 && state != cSelectorUpdateTableCurrentState) {
        return pymol::make_error("state ", state + 1,
                                 " unsupported (must be -1 (current) or >=1)");
      } else {
        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          base[0].sele[a] = false;
          obj = I->Obj[I->Table[a].model];
          if (obj != cur_obj) { /* different object */
            cs = obj->getCoordSet(state);
            cur_obj = obj;
          }
          if (cs) {
            if (cs->atmToIdx(I->Table[a].atom) >= 0) {
              base[0].sele[a] = true;
              c++;
            }
          }
        }
      }
      break;
    case SELE_ALTs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchAlpha(
                   matcher, I->Obj[table_a.model]->AtomInfo[table_a.atom].alt)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_FLGs:
      sscanf(base[1].text(), "%d", &flag);
      flag = (1 << flag);
      for (a = cNDummyAtoms; a < I_NAtom; a++) {
        if (I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].flags &
            flag) {
          base[0].sele[a] = true;
          c++;
        } else
          base[0].sele[a] = false;
      }
      break;
    case SELE_NTYs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigInteger(&options);

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          if ((*base_0_sele_a = WordMatcherMatchInteger(
                   matcher,
                   I->Obj[table_a.model]->AtomInfo[table_a.atom].customType)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_RSIs: {
      CWordMatchOptions options;
      AtomInfoType* ai;

      WordMatchOptionsConfigMixed(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          ai = I->Obj[table_a.model]->AtomInfo + table_a.atom;
          char resi[8];
          AtomResiFromResv(resi, sizeof(resi), ai);
          if ((*base_0_sele_a = WordMatcherMatchMixed(matcher, resi, ai->resv)))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_RSNs: {
      CWordMatchOptions options;

      WordMatchOptionsConfigAlphaList(&options, wildcard[0], ignore_case);

      base_0_sele_a = &base[0].sele[cNDummyAtoms];

      if ((matcher = WordMatcherNew(G, base[1].text(), &options, true))) {
        base_0_sele_a = &base[0].sele[cNDummyAtoms];

        for (a = cNDummyAtoms; a < I_NAtom; a++) {
          auto& table_a = I->Table[a];
          auto& resn = I->Obj[table_a.model]->AtomInfo[table_a.atom].resn;
          if ((*base_0_sele_a =
                   WordMatcherMatchAlpha(matcher, LexStr(G, resn))))
            c++;
          base_0_sele_a++;
        }
        WordMatcherFree(matcher);
      }
    } break;
    case SELE_SELs: {
      const char* word = base[1].text();
      WordType activeselename = "";
      int enabled_only = false;
      CWordMatchOptions options;

      if (word[0] == '?') {
        word++;
        if (word[0] == '?') {
          ExecutiveGetActiveSeleName(G, activeselename, false, false);
          enabled_only = true;
          word++;
        }
      }
      WordMatchOptionsConfigAlpha(&options, wildcard[0], ignore_case);

      if ((matcher = WordMatcherNew(G, word, &options, false))) {
        for (a = 0; a < I_NAtom;
             a++) /* zero out first before iterating through selections */
          base[0].sele[a] = false;

        for (const auto& rec : I->mgr->Info) {
          if (rec.name.empty()) {
            // TODO Can this happen? Why?
            PRINTFB(G, FB_Selector, FB_Warnings)
            " Selector-Unexpected: Empty selection name (ID:%d)\n",
                rec.ID ENDFB(G);
            break;
          }
          if (WordMatcherMatchAlpha(matcher, rec.name.c_str())) {
            if (!enabled_only || activeselename == rec.name) {
              for (a = cNDummyAtoms; a < I_NAtom; a++) {
                s = I->Obj[I->Table[a].model]
                        ->AtomInfo[I->Table[a].atom]
                        .selEntry;
                while (s) {
                  if (I->mgr->Member[s].selection == rec.ID) {
                    if (!base[0].sele[a]) {
                      base[0].sele[a] = I->mgr->Member[s].tag;
                      c++;
                    }
                  }
                  s = IM->Member[s].next;
                }
              }
            }
          }
        }
        WordMatcherFree(matcher);

        /* must also allow for group name pattern matches */

        {
          int group_list_id;
          if ((group_list_id =
                   ExecutiveGetExpandedGroupListFromPattern(G, word))) {
            int last_was_member = false;
            last_obj = nullptr;
            for (a = cNDummyAtoms; a < I_NAtom; a++) {
              if (last_obj != I->Obj[I->Table[a].model]) {
                last_obj = I->Obj[I->Table[a].model];
                last_was_member =
                    ExecutiveCheckGroupMembership(G, group_list_id, last_obj);
              }
              if (last_was_member && !base[0].sele[a]) {
                base[0].sele[a] = true;
                c++;
              }
            }
          }
          ExecutiveFreeGroupList(G, group_list_id);
        }

      } else if (!enabled_only ||
                 WordMatchExact(G, activeselename, word, ignore_case)) {
        auto it = SelectGetInfoIter(G, word, 1, ignore_case);
        if (it != IM->Info.end()) {
          for (a = cNDummyAtoms; a < I_NAtom; a++) {
            base[0].sele[a] = false;
            s = I->Obj[I->Table[a].model]->AtomInfo[I->Table[a].atom].selEntry;
            while (s) {
              if (IM->Member[s].selection == it->ID) {
                base[0].sele[a] = IM->Member[s].tag;
                c++;
              }
              s = IM->Member[s].next;
            }
          }
        } else {
          int group_list_id;
          if ((group_list_id = ExecutiveGetExpandedGroupList(G, word))) {
            int last_was_member = false;
            last_obj = nullptr;
            for (a = 0; a < I_NAtom;
                 a++) /* zero out first before iterating through selections */
              base[0].sele[a] = false;
            for (a = cNDummyAtoms; a < I_NAtom; a++) {
              if (last_obj != I->Obj[I->Table[a].model]) {
                last_obj = I->Obj[I->Table[a].model];
                last_was_member =
                    ExecutiveCheckGroupMembership(G, group_list_id, last_obj);
              }
              if ((base[0].sele[a] = last_was_member)) c++;
            }
            ExecutiveFreeGroupList(G, group_list_id);
          } else if (base[1].m_text[0] == '?') { /* undefined ?sele allowed */
            for (a = cNDummyAtoms; a < I_NAtom; a++) base[0].sele[a] = false;
          } else {
            return pymol::make_error("Invalid selection name \"", word, "\".");
          }
        }
      }
    } break;
    case SELE_MODs:

      /* need to change this to handle wildcarded model names */

      /* first, trim off and record the atom index if one exists */

      index = -1;
      auto pos = base[1].m_text.find('`');
      if (pos != std::string::npos) {
        const char* np = base[1].text() + pos;
        if (sscanf(np + 1, "%d", &index) != 1)
          index = -1;
        else
          index--;
        base[1].m_text.resize(pos);
      }
      model = 0;

      {
        CWordMatchOptions options;
        WordMatchOptionsConfigAlpha(&options, wildcard[0], ignore_case);

        if ((matcher = WordMatcherNew(G, base[1].text(), &options, false))) {
          int obj_matches = false;

          for (a = 0; a < I_NAtom;
               a++) /* zero out first before iterating through selections */
            base[0].sele[a] = false;

          base_0_sele_a = &base[0].sele[cNDummyAtoms];
          last_obj = nullptr;
          for (a = cNDummyAtoms; a < I_NAtom; a++) {
            auto& table_a = I->Table[a];
            obj = I->Obj[table_a.model];
            if (obj != last_obj) {
              obj_matches =
                  WordMatcherMatchAlpha(matcher, I->Obj[table_a.model]->Name);
              last_obj = obj;
            }
            if (obj_matches) {
              if ((index < 0) || (table_a.atom == index)) {
                *base_0_sele_a = true;
                c++;
              }
            }
            base_0_sele_a++;
          }
          WordMatcherFree(matcher);
        } else {
          obj = (ObjectMolecule*)ExecutiveFindObjectByName(G, base[1].text());
          if (obj) {
            for (a = cNDummyModels; a < I->Obj.size(); a++)
              if (I->Obj[a] == obj) {
                model = a + 1;
                break;
              }
          }
          if (!model)
            if (sscanf(base[1].text(), "%i", &model) == 1) {
              if (model <= 0)
                model = 0;
              else if (model > I->Obj.size())
                model = 0;
              else if (!I->Obj[model])
                model = 0;
            }
          if (model) {
            model--;
            if (index >= 0) {
              for (a = cNDummyAtoms; a < I_NAtom; a++) {
                if (I->Table[a].model == model)
                  if (I->Table[a].atom == index) {
                    base[0].sele[a] = true;
                    c++;
                  } else {
                    base[0].sele[a] = false;
                  }
                else
                  base[0].sele[a] = false;
              }
            } else {
              for (a = cNDummyAtoms; a < I_NAtom; a++) {
                if (I->Table[a].model == model) {
                  base[0].sele[a] = true;
                  c++;
                } else
                  base[0].sele[a] = false;
              }
            }
          } else {
            return pymol::make_error("invalid model \"", base[1].text(), "\"");
          }
        }
      }
      break;
  }
  PRINTFD(G, FB_Selector)
  " %s:  %d atoms selected.\n", __func__, c ENDFD;
  return {};
}
int SelectorSelect2(PyMOLGlobals* G, EvalElem* base, int state) {
  int a;
  int c = 0;
  int ok = true;
  int oper;
  float comp1;
  int exact;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

  AtomInfoType* at1;
  CSelector* I = G->Selector;
  base->type = STYP_LIST;
  base->sele_calloc(I->Table.size());
  base->sele_err_chk_ptr(G);
  switch (base->code) {
    case SELE_XVLx:
    case SELE_YVLx:
    case SELE_ZVLx:
      oper = WordKey(G, AtOper, base[1].text(), 4, ignore_case, &exact);
      switch (oper) {
        case SCMP_GTHN:
        case SCMP_LTHN:
        case SCMP_EQAL:
        case SCMP_GEQL:
        case SCMP_LEQL:
          if (sscanf(base[2].text(), "%f", &comp1) != 1)
            ok = ErrMessage(G, "Selector", "Invalid Number");
          break;
        default:
          ok = ErrMessage(G, "Selector", "Invalid Operator.");
          break;
      }

      if (ok) {
        ObjectMolecule* obj;
        CoordSet* cs;
        int at, idx, s, s0 = 0, sN = I->NCSet;

        if (state != cStateAll) {
          s0 = (state < cStateAll) ? SceneGetState(G) : state;
          sN = s0 + 1;
        }

        for (s = s0; s < sN; s++) {
          for (a = cNDummyAtoms; a < I->Table.size(); a++) {
            if (base[0].sele[a]) continue;

            obj = I->Obj[I->Table[a].model];
            if (s >= obj->NCSet) continue;

            at = I->Table[a].atom;
            cs = obj->CSet[s];
            idx = cs->atmToIdx(at);
            if (idx < 0) continue;

            idx *= 3;
            switch (base->code) {
              case SELE_ZVLx:
                idx++;
              case SELE_YVLx:
                idx++;
              case SELE_XVLx:
                break;
            }

            if (fcmp(cs->Coord[idx], comp1, oper)) {
              base[0].sele[a] = true;
            } else {
              base[0].sele[a] = false;
            }
          }
        }
      }
      break;
    case SELE_PCHx:
    case SELE_FCHx:
    case SELE_BVLx:
    case SELE_QVLx:
      oper = WordKey(G, AtOper, base[1].text(), 4, ignore_case, &exact);
      if (!oper) {
        ok = ErrMessage(G, "Selector", "Invalid Operator.");
        break;
      }
      switch (oper) {
        case SCMP_GTHN:
        case SCMP_LTHN:
        case SCMP_EQAL:
        case SCMP_GEQL:
        case SCMP_LEQL:
          if (sscanf(base[2].text(), "%f", &comp1) != 1)
            ok = ErrMessage(G, "Selector", "Invalid Number");
          break;
        default:
          ok = ErrMessage(G, "Selector", "Invalid Operator.");
          break;
      }

      if (ok) {
        for (a = cNDummyAtoms; a < I->Table.size(); a++) {
          at1 = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
          float atomValue;
          switch (base->code) {
            case SELE_BVLx:
              atomValue = at1->b;
              break;
            case SELE_PCHx:
              atomValue = at1->partialCharge;
              break;
            case SELE_FCHx:
              atomValue = at1->formalCharge;
              break;
            case SELE_QVLx:
              atomValue = at1->q;
              break;
            default: {
              ErrMessage(G, "Selector", "Invalid Operand.");
              return false;
            }
          }

          if (fcmp(atomValue, comp1, oper)) {
            base[0].sele[a] = true;
            c++;
          } else {
            base[0].sele[a] = false;
          }
        }
      }
      break;
  }

  PRINTFD(G, FB_Selector)
  " %s: %d atoms selected.\n", __func__, c ENDFD;
  return (ok);
}
pymol::Result<> SelectorSelect3(PyMOLGlobals* G, EvalElem* base, int state) {
  int a;
  int oper;
  float comp1, fval;
  int exact;
  AtomInfoType* at1;
  CSelector* I = G->Selector;

  base->type = STYP_LIST;
  base->sele_calloc(I->Table.size());
  base->sele_err_chk_ptr(G);
  switch (base->code) {
    case SELE_PROP:
      oper = WordKey(G, AtOper, base[2].text(), 4, 0, &exact);
      switch (oper) {
        case SCMP_GTHN:
        case SCMP_LTHN:
        case SCMP_EQAL:
          if (sscanf(base[3].text(), "%f", &comp1) != 1) {
            return pymol::make_error("Invalid Number: ", base[3].text());
          }
          for (a = cNDummyAtoms; a < I->Table.size(); a++) {
            at1 = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
            if (at1->prop_id) {
              fval = PropertyGetAsFloat(G, at1->prop_id, base[1].text());
              base[0].sele[a] = fcmp(fval, comp1, oper);
            }
          }
          break;
        case SCMP_RANG: {
          char buffer[64];
          CWordMatchOptions options;
          CWordMatcher* matcher;
          WordMatchOptionsConfigAlphaList(&options, '*', 0);
          if (!(matcher = WordMatcherNew(G, base[3].text(), &options, true))) {
            assert(false);  // should never fail
            return pymol::Error();
          }
          for (a = cNDummyAtoms; a < I->Table.size(); a++) {
            at1 = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
            if (at1->prop_id) {
              auto sval =
                  PropertyGetAsString(G, at1->prop_id, base[1].text(), buffer);
              base[0].sele[a] = sval && WordMatcherMatchAlpha(matcher, sval);
            }
          }
          WordMatcherFree(matcher);
        } break;
        default:
          return pymol::make_error("Invalid Operator: ", base[2].text());
      }
      break;
    default:
      assert(false);
  }
  return {};
}
/**
 * Ring finder subroutine
 * Modifies base[0].sele
 */
class SelectorRingFinder : public AbstractRingFinder {
  CSelector* m_selector;
  EvalElem* m_base;

 protected:
  void onRingFound(ObjectMolecule* obj, const int* indices,
                   size_t size) override {
    for (size_t i = 0; i < size; ++i) {
      int offset = SelectorGetObjAtmOffset(m_selector, obj, indices[i]);
      if (offset >= 0) m_base->sele[offset] = 1;
    }
  }

 public:
  SelectorRingFinder(CSelector* selector, EvalElem* base, int maxringsize = 7)
      : AbstractRingFinder(maxringsize), m_selector(selector), m_base(base) {}
};
int SelectorLogic1(PyMOLGlobals* G, EvalElem* inp_base, int state) {
  /* some cases in this function still need to be optimized
     for performance (see BYR1 for example) */

  CSelector* I = G->Selector;
  int a, b, tag;
  int c = 0;
  int flag;
  EvalElem* base = inp_base;
  AtomInfoType *at1, *at2;
  int n_atom = I->Table.size();
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  ObjectMolecule* lastObj = nullptr;

  base[0].sele = std::move(base[1].sele);
  base[0].type = STYP_LIST;
  switch (base->code) {
    case SELE_NOT1: {
      int* base_0_sele_a;

      base_0_sele_a = base[0].sele_data();
      for (a = 0; a < n_atom; a++) {
        if ((*base_0_sele_a = !*base_0_sele_a)) c++;
        base_0_sele_a++;
      }
    } break;
    case SELE_RING: {
      std::vector<bool> selemask(base[0].sele_data(),
                                 base[0].sele_data() + n_atom);
      SelectorRingFinder ringfinder(I, base);

      std::fill_n(base[0].sele_data(), n_atom, 0);

      for (SelectorAtomIterator iter(I); iter.next();) {
        if (selemask[iter.a]) ringfinder.apply(iter.obj, iter.getAtm());
      }
    } break;
    case SELE_NGH1:
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);

      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = base[1].sele[a])) {
          auto const* lastObj = I->Obj[table_a.model];
          for (auto const& neighbor : AtomNeighbors(lastObj, table_a.atom)) {
            auto const a2 = SelectorGetObjAtmOffset(I, lastObj, neighbor.atm);
            if (a2 >= 0) {
              if (!base[1].sele[a2]) base[0].sele[a2] = tag;
            }
          }
        }
      }
      base[1].sele_free();
      break;
    case SELE_BON1:
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = base[1].sele[a])) {
          auto const* lastObj = I->Obj[table_a.model];
          for (auto const& neighbor : AtomNeighbors(lastObj, table_a.atom)) {
            auto const a2 = SelectorGetObjAtmOffset(I, lastObj, neighbor.atm);
            assert(a2 >= 0);
            if (a2 >= 0 && !base[0].sele[a2]) {
              base[0].sele[a2] = 1;
            }
          }
        }
      }
      base[1].sele_free();
      break;
    case SELE_BYO1:
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);
      for (a = cNDummyAtoms; a < n_atom; a++) {
        if (base[1].sele[a]) {
          if (I->Obj[I->Table[a].model] != lastObj) {
            lastObj = I->Obj[I->Table[a].model];
            b = a;
            while (b >= 0) {
              if (I->Obj[I->Table[b].model] != lastObj) break;
              base[0].sele[b] = 1;
              b--;
            }
            b = a + 1;
            while (b < n_atom) {
              if (I->Obj[I->Table[b].model] != lastObj) break;
              base[0].sele[b] = 1;
              b++;
            }
          }
        }
      }
      base[1].sele_free();
      break;
    case SELE_BYR1: /* ASSUMES atoms are sorted & grouped by residue */
    case SELE_CAS1: {
      int* base_0_sele = base[0].sele_data();
      int break_atom = -1;
      int last_tag = 0;
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = base_0_sele[a]) &&
            ((a >= break_atom) || (base_0_sele[a] != last_tag))) {
          at1 = &I->Obj[table_a.model]->AtomInfo[table_a.atom];
          b = a - 1;
          while (b >= 0) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (AtomInfoSameResidue(G, at1, at2)) {
                  base_0_sele[b] = tag;
                  c++;
                  flag = 1;
                }
              }
              if (!flag) {
                break;
              }
            }
            b--;
          }
          b = a + 1;
          while (b < n_atom) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (AtomInfoSameResidue(G, at1, at2)) {
                  base_0_sele[b] = tag;
                  c++;
                  flag = 1;
                }
              }
              if (!flag) {
                break_atom = b - 1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
      }
      if (base->code == SELE_CAS1) {
        c = 0;
        for (a = cNDummyAtoms; a < n_atom; a++) {
          auto& table_a = I->Table[a];
          if (base_0_sele[a]) {
            base_0_sele[a] = false;

            if (I->Obj[table_a.model]->AtomInfo[table_a.atom].protons == cAN_C)
              if (WordMatchExact(
                      G, G->lex_const.CA,
                      I->Obj[table_a.model]->AtomInfo[table_a.atom].name,
                      ignore_case)) {
                base_0_sele[a] = true;
                c++;
              }
          }
        }
      }
    } break;
    case SELE_BYC1: /* ASSUMES atoms are sorted & grouped by chain */
    {
      int* base_0_sele = base[0].sele_data();
      int break_atom_high = -1;
      int break_atom_low = 0;
      int last_tag = 0;
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = base_0_sele[a]) &&
            ((a >= break_atom_high) || (base_0_sele[a] != last_tag))) {
          if (tag != last_tag) break_atom_low = 0;
          at1 = &I->Obj[table_a.model]->AtomInfo[table_a.atom];
          b = a - 1;
          while (b >= break_atom_low) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (at1->chain == at2->chain)
                  if (at1->segi == at2->segi) {
                    base_0_sele[b] = tag;
                    c++;
                    flag = 1;
                  }
              }
              if (!flag) {
                break_atom_low = b + 1;
                break;
              }
            }
            b--;
          }
          if (b < 0) break_atom_low = 0;
          b = a + 1;
          while (b < n_atom) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (at1->chain == at2->chain)
                  if (at1->segi == at2->segi) {
                    base_0_sele[b] = tag;
                    c++;
                    flag = 1;
                  }
              }
              if (!flag) {
                break_atom_high = b - 1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
      }
    } break;
    case SELE_BYS1: /* ASSUMES atoms are sorted & grouped by segi */
    {
      int* base_0_sele = base[0].sele_data();
      int break_atom_high = -1;
      int break_atom_low = 0;
      int last_tag = 0;
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = base_0_sele[a]) &&
            ((a >= break_atom_high) || (base_0_sele[a] != last_tag))) {
          if (tag != last_tag) break_atom_low = 0;
          at1 = &I->Obj[table_a.model]->AtomInfo[table_a.atom];
          b = a - 1;
          while (b >= break_atom_low) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (at1->segi == at2->segi) {
                  base_0_sele[b] = tag;
                  c++;
                  flag = 1;
                }
              }
              if (!flag) {
                break_atom_low = b + 1;
                break;
              }
            }
            b--;
          }
          b = a + 1;
          while (b < n_atom) {
            if (!base_0_sele[b]) {
              flag = false;
              if (table_a.model == I->Table[b].model) {
                at2 = &I->Obj[I->Table[b].model]->AtomInfo[I->Table[b].atom];
                if (at1->segi == at2->segi) {
                  base_0_sele[b] = tag;
                  c++;
                  flag = 1;
                }
              }
              if (!flag) {
                break_atom_high = b - 1;
                last_tag = tag;
                break;
              }
            }
            b++;
          }
        }
      }
    } break;
    case SELE_BYF1: /* first, identify all atom by fragment selection */
    {
      /* NOTE: this algorithm looks incompatible with selection
         tags...need to do some more thinking & work... */

      int n_frag = EditorGetNFrag(G);

      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);

      if (n_frag) {
        int a, f, at, s;
        int* fsele;
        ObjectMolecule* obj;

        fsele = pymol::malloc<int>(n_frag + 1);

        for (f = 0; f < n_frag; f++) {
          auto name = pymol::string_format("%s%1d", cEditorFragPref, f + 1);
          fsele[f] = SelectorIndexByName(G, name.c_str());
        }

        /* mark atoms by fragment */
        for (a = 0; a < n_atom; a++) {
          at = I->Table[a].atom;
          obj = I->Obj[I->Table[a].model];
          s = obj->AtomInfo[at].selEntry;
          for (f = 0; f < n_frag; f++) {
            if (SelectorIsMember(G, s, fsele[f])) {
              base[0].sele[a] = f + 1;
            }
          }
        }

        /* mark fragments we keep */
        for (f = 0; f <= n_frag; f++) {
          fsele[f] = 0;
        }
        for (a = 0; a < n_atom; a++) {
          int f = base[0].sele[a];
          if (base[1].sele[a] && f) fsele[f] = 1;
        }

        /* now set flags */
        for (a = 0; a < n_atom; a++) {
          c += (base[0].sele[a] = fsele[base[0].sele[a]]);
        }

        FreeP(fsele);
      }
      base[1].sele_free();
    } break;
    case SELE_BYM1: {
      int c = 0;
      int a, at, aa;
      int* stk;
      int stkDepth = 0;
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);

      stk = VLAlloc(int, 50);

      for (a = 0; a < n_atom; a++) {
        if ((tag = base[1].sele[a]) && (!base[0].sele[a])) {
          VLACheck(stk, int, stkDepth);
          stk[stkDepth] = a;
          stkDepth++;

          auto const obj = I->Obj[I->Table[a].model];

          while (stkDepth) { /* this will explore a tree */
            stkDepth--;
            a = stk[stkDepth];
            base[0].sele[a] = tag;
            c++;
            at = I->Table[a].atom; /* start walk from this location */

            /* add neighbors onto the stack */
            for (auto const& neighbor : AtomNeighbors(obj, at)) {
              if (obj->Bond[neighbor.bond].order > 0) {
                if ((aa = SelectorGetObjAtmOffset(I, obj, neighbor.atm)) >= 0) {
                  if (!base[0].sele[aa]) {
                    VLACheck(stk, int, stkDepth);
                    stk[stkDepth] = aa; /* add index in selector space */
                    stkDepth++;
                  }
                }
              }
            }
          }
        }
      }
      base[1].sele_free();
      VLAFreeP(stk);
    } break;
    case SELE_BYX1: /* by cell */
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);
      {
        ObjectMolecule* obj;
        CoordSet* cs;
        int d, n1, at;
        for (d = 0; d < I->NCSet; d++) {
          if ((state < 0) || (d == state)) {
            n1 = 0;

            auto Flag1 = std::vector<MapFlag_t>(I->Table.size(), 0);
            auto Vertex = std::vector<float>(I->Table.size() * 3, 0.0f);

            for (a = 0; a < I->Table.size(); a++) {
              at = I->Table[a].atom;
              obj = I->Obj[I->Table[a].model];
              if (d < obj->NCSet)
                cs = obj->CSet[d];
              else
                cs = nullptr;
              if (cs) {
                const auto* sym = cs->getSymmetry();
                if (sym) {
                  int idx;
                  idx = cs->atmToIdx(at);
                  if (idx >= 0) {
                    transform33f3f(sym->Crystal.realToFrac(), cs->coordPtr(idx),
                                   Vertex.data() + 3 * a);
                    Flag1[a] = true;
                    n1++;
                  }
                }
              }
            }
            if (n1) {
              std::unique_ptr<MapType> map(new MapType(G, -1.1, Vertex.data(),
                                                       I->Table.size(), nullptr,
                                                       Flag1.data()));
              if (map) {
                int e, nCSet;
                nCSet = SelectorGetArrayNCSet(G, base[1].sele, false);
                for (e = 0; e < nCSet; e++) {
                  if ((state < 0) || (e == state)) {
                    for (a = 0; a < I->Table.size(); a++) {
                      if (base[1].sele[a]) {
                        at = I->Table[a].atom;
                        obj = I->Obj[I->Table[a].model];
                        if (e < obj->NCSet)
                          cs = obj->CSet[e];
                        else
                          cs = nullptr;
                        if (cs) {
                          const auto* sym = cs->getSymmetry();
                          if (sym) {
                            int idx;
                            idx = cs->atmToIdx(at);
                            if (idx >= 0) {
                              float probe[3], probe_i[3];

                              transform33f3f(sym->Crystal.realToFrac(),
                                             cs->coordPtr(idx), probe);
                              probe_i[0] = (int)floor(probe[0]);
                              probe_i[1] = (int)floor(probe[1]);
                              probe_i[2] = (int)floor(probe[2]);

                              for (const auto j :
                                   MapEIter(*map, probe, false)) {
                                if (!base[0].sele[j]) {
                                  float* tst = Vertex.data() + 3 * j;
                                  base[0].sele[j] =
                                      ((probe_i[0] == (int)floor(tst[0])) &&
                                       (probe_i[1] == (int)floor(tst[1])) &&
                                       (probe_i[2] == (int)floor(tst[2])));
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
      base[1].sele_free();
      break;
    case SELE_FST1:
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);
      for (a = cNDummyAtoms; a < n_atom; a++) {
        if (base[1].sele[a]) {
          base[0].sele[a] = base[1].sele[a]; /* preserve tag */
          break;
        }
      }
      base[1].sele_free();
      break;
    case SELE_LST1: {
      int last = -1;
      base[1].sele = std::move(base[0].sele);
      base[0].sele_calloc(n_atom);
      for (a = cNDummyAtoms; a < n_atom; a++) {
        if (base[1].sele[a]) {
          last = a;
        }
      }
      if (last >= 0) base[0].sele[last] = base[1].sele[last]; /* preserve tag */
    }
      base[1].sele_free();
      break;
  }
  PRINTFD(G, FB_Selector)
  " %s: %d atoms selected.\n", __func__, c ENDFD;
  return (1);
}
int SelectorLogic2(PyMOLGlobals* G, EvalElem* base) {
  CSelector* I = G->Selector;
  int a, b, tag;
  int c = 0;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  int ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);
  int *base_0_sele_a, *base_2_sele_a;
  int n_atom = I->Table.size();

  AtomInfoType *at1, *at2;

  switch (base[1].code) {
    case SELE_OR_2:
    case SELE_IOR2: {
      base_0_sele_a = base[0].sele_data();
      base_2_sele_a = base[2].sele_data();

      for (a = 0; a < n_atom; a++) {
        if (((*base_0_sele_a) =
                 (((*base_0_sele_a) > (*base_2_sele_a)) ? (*base_0_sele_a)
                                                        : (*base_2_sele_a)))) {
          /* use higher tag */
          c++;
        }
        base_0_sele_a++;
        base_2_sele_a++;
      }
    } break;
    case SELE_AND2:

      base_0_sele_a = base[0].sele_data();
      base_2_sele_a = base[2].sele_data();

      for (a = 0; a < n_atom; a++) {
        if ((*base_0_sele_a) && (*base_2_sele_a)) {
          (*base_0_sele_a) =
              (((*base_0_sele_a) > (*base_2_sele_a)) ? (*base_0_sele_a)
                                                     : (*base_2_sele_a));
          /* use higher tag */
          c++;
        } else {
          (*base_0_sele_a) = 0;
        }
        base_0_sele_a++;
        base_2_sele_a++;
      }
      break;
    case SELE_ANT2:
      base_0_sele_a = base[0].sele_data();
      base_2_sele_a = base[2].sele_data();

      for (a = 0; a < n_atom; a++) {
        if ((*base_0_sele_a) && !(*base_2_sele_a)) {
          c++;
        } else {
          (*base_0_sele_a) = 0;
        }
        base_0_sele_a++;
        base_2_sele_a++;
      }
      break;
    case SELE_IN_2: {
      int* base_2_sele_b;
      base_0_sele_a = &base[0].sele[cNDummyAtoms];
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = *base_0_sele_a)) {
          at1 = &I->Obj[table_a.model]->AtomInfo[table_a.atom];
          *base_0_sele_a = 0;
          base_2_sele_b = &base[2].sele[cNDummyAtoms];
          for (b = cNDummyAtoms; b < n_atom; b++) {
            auto& table_b = I->Table[b];
            if (*base_2_sele_b) {
              at2 = &I->Obj[table_b.model]->AtomInfo[table_b.atom];
              if (at1->resv == at2->resv)
                if (WordMatchExact(G, at1->chain, at2->chain,
                                   ignore_case_chain))
                  if (WordMatchExact(G, at1->name, at2->name, ignore_case))
                    if (WordMatchExact(G, at1->inscode, at2->inscode,
                                       ignore_case))
                      if (WordMatchExact(G, at1->resn, at2->resn, ignore_case))
                        if (WordMatchExact(G, at1->segi, at2->segi,
                                           ignore_case_chain)) {
                          *base_0_sele_a = tag;
                          break;
                        }
            }
            base_2_sele_b++;
          }
        }
        if (*(base_0_sele_a++)) c++;
      }
    } break;
    case SELE_LIK2: {
      int* base_2_sele_b;
      base_0_sele_a = &base[0].sele[cNDummyAtoms];
      for (a = cNDummyAtoms; a < n_atom; a++) {
        auto& table_a = I->Table[a];
        if ((tag = *base_0_sele_a)) {
          at1 = &I->Obj[table_a.model]->AtomInfo[table_a.atom];
          *base_0_sele_a = 0;
          base_2_sele_b = &base[2].sele[cNDummyAtoms];
          for (b = cNDummyAtoms; b < n_atom; b++) {
            auto& table_b = I->Table[b];
            if (*base_2_sele_b) {
              at2 = &I->Obj[table_b.model]->AtomInfo[table_b.atom];
              if (at1->resv == at2->resv)
                if (WordMatchExact(G, at1->name, at2->name, ignore_case))
                  if (WordMatchExact(G, at1->inscode, at2->inscode,
                                     ignore_case)) {
                    *base_0_sele_a = tag;
                    break;
                  }
            }
            base_2_sele_b++;
          }
        }
        if (*(base_0_sele_a++)) c++;
      }
    } break;
  }
  base[2].sele_free();
  PRINTFD(G, FB_Selector)
  " %s: %d atoms selected.\n", __func__, c ENDFD;
  return (1);
}
int SelectorOperator22(PyMOLGlobals* G, EvalElem* base, int state) {
  int c = 0;
  int a, d, e;
  CSelector* I = G->Selector;
  ObjectMolecule* obj;

  float dist;
  CoordSet* cs;
  int ok = true;
  int nCSet;
  int n1, at, idx;
  int code = base[1].code;

  if (state < 0) {
    switch (state) {
      case cSelectorUpdateTableCurrentState:
      case cSelectorUpdateTableEffectiveStates:
        state = SceneGetState(G);
        break;
    }
  }

  switch (code) {
    case SELE_WIT_:
    case SELE_BEY_:
    case SELE_NTO_:
      if (!sscanf(base[2].text(), "%f", &dist))
        ok = ErrMessage(G, "Selector", "Invalid distance.");
      if (ok) {
        if (dist < 0.0) dist = 0.0;

        const size_t table_size = I->Table.size();
        auto coords_flat = std::vector<float>(table_size * 3);
        auto* coords = pymol::reshape<3>(coords_flat.data());

        /* copy starting mask */
        const auto Flag2 = std::move(base[0].sele);
        base[0].sele_calloc(table_size);

        for (d = 0; d < I->NCSet; d++) {
          if ((state < 0) || (d == state)) {
            n1 = 0;
            auto Flag1 = std::vector<MapFlag_t>(table_size);
            for (a = 0; a < table_size; a++) {
              at = I->Table[a].atom;
              obj = I->Obj[I->Table[a].model];
              if (d < obj->NCSet)
                cs = obj->CSet[d];
              else
                cs = nullptr;
              if (cs) {
                if (CoordSetGetAtomVertex(cs, at, coords[a])) {
                  Flag1[a] = true;
                  n1++;
                }
              }
            }
            if (n1) {
              std::unique_ptr<MapType> map(
                  new MapType(G, -dist, pymol::flatten(coords), table_size,
                              nullptr, Flag1.data()));
              CHECKOK(ok, map);
              if (ok) {
                nCSet = SelectorGetArrayNCSet(G, base[4].sele, false);
                for (e = 0; ok && e < nCSet; e++) {
                  if ((state < 0) || (e == state)) {
                    for (a = 0; a < I->Table.size(); a++) {
                      if (base[4].sele[a]) {
                        at = I->Table[a].atom;
                        obj = I->Obj[I->Table[a].model];
                        if (e < obj->NCSet)
                          cs = obj->CSet[e];
                        else
                          cs = nullptr;
                        if (cs) {
                          idx = cs->atmToIdx(at);
                          if (idx >= 0) {
                            const float* v2 = cs->coordPtr(idx);
                            for (const auto j : MapEIter(*map, v2, false)) {
                              if (!base[0].sele[j] && Flag2[j] &&
                                  within3f(coords[j], v2, dist) &&
                                  (code != SELE_NTO_ || !base[4].sele[j])) {
                                base[0].sele[j] = true;
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
        if (code == SELE_BEY_) {
          for (a = 0; a < I->Table.size(); a++) {
            if (Flag2[a]) base[0].sele[a] = !base[0].sele[a];
          }
        }
        for (a = cNDummyAtoms; a < I->Table.size(); a++)
          if (base[0].sele[a]) c++;
      }
      break;
  }
  base[4].sele_free();
  PRINTFD(G, FB_Selector)
  " %s: %d atoms selected.\n", __func__, c ENDFD;
  return (1);
}
