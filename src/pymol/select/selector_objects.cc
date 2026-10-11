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

// Creating molecules from selections and copying coordinates in and out.

#include "pymol/base/feedback.h"
#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/select/selector_internal.h"
static bool SelectorIsSelectionDiscrete(PyMOLGlobals* G, SelectorID_t sele,
                                        bool update_table) {
  CSelector* I = G->Selector;

  if (update_table) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  }

  for (int a = cNDummyAtoms; a < I->Table.size(); a++) {
    auto& table_a = I->Table[a];
    auto obj = I->Obj[table_a.model];
    auto ai = obj->AtomInfo + table_a.atom;
    if (SelectorIsMember(G, ai->selEntry, sele)) {
      if (obj->DiscreteFlag) {
        return true;
      }
    }
  }
  return false;
}
/**
 * Get selection coordinates as Nx3 numpy array. Equivalent to
 *
 *     PyMOL> coords = []
 *     PyMOL> cmd.iterate_state(state, sele, 'coords.append([x,y,z])')
 *     PyMOL> coords = numpy.array(coords)
 */
PyObject* SelectorGetCoordsAsNumPy(PyMOLGlobals* G, int sele, int state) {
  printf("No numpy support\n");
  return nullptr;
}
/**
 * Load coordinates from a Nx3 sequence into the given selection.
 * Most efficiant with numpy arrays. Equivalent to
 *
 * PyMOL> coords = iter(coords)
 * PyMOL> cmd.alter_state(state, sele, '(x,y,z) = coords.next()')
 */
pymol::Result<> SelectorLoadCoords(PyMOLGlobals* G, PyObject* coords, int sele,
                                   int state) {
  double matrix[16];
  double* matrix_ptr = nullptr;
  float v_xyz[3];
  int a, b, nAtom = 0, itemsize;
  SeleCoordIterator iter(G, sele, state);
  CoordSet* mat_cs = nullptr;
  PyObject *v, *w;
  bool is_np_array = false;
  void* ptr;

  if (!PySequence_Check(coords)) {
    return pymol::make_error("Passed argument is not a sequence");
  }

  // atom count in selection
  while (iter.next()) nAtom++;

  // sequence length must match atom count
  if (nAtom != PySequence_Size(coords)) {
    return pymol::make_error("Atom count mismatch");
  }

  // detect numpy arrays, allows faster data access (see below)

  for (a = 0, iter.reset(); iter.next(); a++) {
    // get xyz from python
    if (is_np_array) {
      // fast implementation for numpy arrays only
    } else {
      // general implementation for any 2d sequence
      v = PySequence_ITEM(coords, a);

      // get xyz from python sequence item
      for (b = 0; b < 3; b++) {
        if (!(w = PySequence_GetItem(v, b))) break;

        v_xyz[b] = (float)PyFloat_AsDouble(w);
        Py_DECREF(w);
      }

      Py_DECREF(v);
    }

    if (PyErr_Occurred()) {
      return pymol::make_error("Load Coords error occurred.");
    }

    // coord set specific stuff
    if (mat_cs != iter.cs) {
      // update matrix
      matrix_ptr = ObjectGetTotalMatrix(iter.obj, state, false, matrix)
                       ? matrix
                       : nullptr;
      mat_cs = iter.cs;

      // invalidate reps
      iter.cs->invalidateRep(cRepAll, cRepInvRep);
    }

    // handle matrix
    if (matrix_ptr) {
      inverse_transform44d3f(matrix_ptr, v_xyz, v_xyz);
    }

    // copy coordinates
    copy3f(v_xyz, iter.getCoord());
  }

  return {};
}
pymol::Result<> SelectorUpdateCmd(PyMOLGlobals* G,     //
                                  SelectorID_t sele0,  //
                                  SelectorID_t sele1,  //
                                  int sta0, int sta1, int matchmaker,
                                  int quiet) {
  CSelector* I = G->Selector;
  int a, b;
  int at0 = 0, at1;
  int c0 = 0, c1 = 0;
  int i0 = 0, i1;
  ObjectMolecule *obj0 = nullptr, *obj1;
  CoordSet* cs0;
  const CoordSet* cs1;
  int matched_flag;
  int b_start;
  int ccc = 0;

  bool ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  bool ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);

  PRINTFD(G, FB_Selector)
  " %s-Debug: entered sta0 %d sta1 %d", __func__, sta0, sta1 ENDFD;

  // either both or none must be "all states"
  if (sta0 != sta1) {
    if (sta0 == cSelectorUpdateTableAllStates) {
      sta0 = sta1;
    } else if (sta1 == cSelectorUpdateTableAllStates) {
      sta1 = sta0;
    }
  }

  if ((sta0 < 0) || (sta1 < 0) || (sta0 != sta1)) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  } else {
    SelectorUpdateTable(G, sta0, -1);
  }

  auto vla0 = pymol::vla_take_ownership(SelectorGetIndexVLA(G, sele0));
  auto vla1 = pymol::vla_take_ownership(SelectorGetIndexVLA(G, sele1));

  if (vla0 && vla1) {
    c0 = VLAGetSize(vla0);
    c1 = VLAGetSize(vla1);
  }

  if (c0 < 1 || c1 < 1)
    return pymol::make_error("No coordinates updated.");
  else {
    b = 0;
    for (a = 0; a < c1; a++) { /* iterate over source atoms */
      /* NOTE, this algorithm is N^2 and slow in the worst case...
         however the best case (N) is quite common, especially when merging
         files written out of PyMOL */

      i1 = vla1[a];
      at1 = I->Table[i1].atom;
      obj1 = I->Obj[I->Table[i1].model];
      matched_flag = false;

      switch (matchmaker) {
        case 0:
          /* simply assume that atoms are stored in PyMOL in the identical
           * order, one for one */
          if (b < c0) {
            i0 = vla0[b];
            at0 = I->Table[i0].atom;
            obj0 = I->Obj[I->Table[i0].model];
            b++;
            matched_flag = true;
          }
          break;
        case 1:
          /* match each pair based on atom info */
          b_start = b;
          matched_flag = false;
          while (1) {
            i0 = vla0[b];
            at0 = I->Table[i0].atom;
            obj0 = I->Obj[I->Table[i0].model];
            if (obj0 != obj1) {
              if (AtomInfoMatch(G, obj1->AtomInfo + at1, obj0->AtomInfo + at0,
                                ignore_case, ignore_case_chain)) {
                matched_flag = true;
                break;
              }
            } else if (at0 == at1) {
              matched_flag = true;
              break;
            }
            b++;
            if (b >= c0) b = 0;
            if (b == b_start) break;
          }
          break;
        case 2: /* match based on ID */
        {
          int target = obj1->AtomInfo[at1].id;
          b_start = b;
          matched_flag = false;
          while (1) {
            i0 = vla0[b];
            at0 = I->Table[i0].atom;
            obj0 = I->Obj[I->Table[i0].model];
            if (obj0 != obj1) {
              if (obj0->AtomInfo[at0].id == target) {
                matched_flag = true;
                break;
              }
            } else if (at0 == at1) {
              matched_flag = true;
              break;
            }
            b++;
            if (b >= c0) b = 0;
            if (b == b_start) break;
          }
        } break;
        case 3: /* match based on rank */
        {
          int target = obj1->AtomInfo[at1].rank;
          b_start = b;
          matched_flag = false;
          while (1) {
            i0 = vla0[b];
            at0 = I->Table[i0].atom;
            obj0 = I->Obj[I->Table[i0].model];
            if (obj0 != obj1) {
              if (obj0->AtomInfo[at0].rank == target) {
                matched_flag = true;
                break;
              }
            } else if (at0 == at1) {
              matched_flag = true;
            }
            b++;
            if (b >= c0) b = 0;
            if (b == b_start) break;
          }
        } break;
        case 4: /* match based on index */
        {
          b_start = b;
          matched_flag = false;
          while (1) {
            i0 = vla0[b];
            at0 = I->Table[i0].atom;
            obj0 = I->Obj[I->Table[i0].model];
            if (obj0 != obj1) {
              if (at0 == at1) {
                matched_flag = true;
                break;
              }
            } else if (at0 == at1) {
              matched_flag = true;
              break;
            }
            b++;
            if (b >= c0) b = 0;
            if (b == b_start) break;
          }
        } break;
      }

      if (matched_flag) { /* atom matched, so copy coordinates */
        ccc++;

        StateIterator iter0(G, obj0->Setting.get(), sta0, obj0->NCSet);
        StateIterator iter1(G, obj1->Setting.get(), sta1, obj1->NCSet);

        while (iter0.next() && iter1.next()) {
          cs0 = obj0->CSet[iter0.state];
          cs1 = obj1->CSet[iter1.state];
          if (cs1 && cs0) {
            int idx0 = cs0->atmToIdx(at0);
            int idx1 = cs1->atmToIdx(at1);
            if (idx0 >= 0 && idx1 >= 0) {
              copy3f(cs1->coordPtr(idx1), cs0->coordPtr(idx0));
            }
          }
        }
      }
    }
    obj0 = nullptr;

    {
      ObjectMolecule** objs = SelectorGetObjectMoleculeVLA(G, sele0);
      int sz = VLAGetSize(objs);
      for (b = 0; b < sz; b++) {
        objs[b]->invalidate(cRepAll, cRepInvCoord, -1);
        ExecutiveUpdateCoordDepends(G, objs[b]);
      }
      VLAFree(objs);
    }
    SceneChanged(G);
    if (!quiet) {
      PRINTFB(G, FB_Selector, FB_Actions)
      " Update: coordinates updated for %d atoms.\n", ccc ENDFB(G);
    }
  }
  return {};
}
int SelectorCreateObjectMolecule(PyMOLGlobals* G, SelectorID_t sele,
                                 const char* name, int target, int source,
                                 int discrete, int zoom, int quiet,
                                 int singletons, int copy_properties) {
  CSelector* I = G->Selector;
  int ok = true;
  int a, b, a2, b1, b2, c, d, s, at;
  const BondType* ii1;
  int nBond = 0;
  int nCSet, nAtom, ts;
  int isNew;
  CoordSet* cs = nullptr;
  CoordSet *cs1, *cs2;
  ObjectMolecule* obj;
  pymol::CObject* ob;
  ObjectMolecule* targ = nullptr;
  ObjectMolecule* info_src = nullptr;
  int static_singletons = SettingGetGlobal_b(G, cSetting_static_singletons);

  if (singletons < 0) singletons = static_singletons;

  ob = ExecutiveFindObjectByName(G, name);
  if (ob)
    if (ob->type == cObjectMolecule) targ = (ObjectMolecule*)ob;

  SelectorUpdateTable(G, source, -1);

  if (!targ) {
    isNew = true;
    if (discrete < 0) discrete = SelectorIsSelectionDiscrete(G, sele, false);
    targ = new ObjectMolecule(G, discrete);
    targ->Bond = pymol::vla<BondType>(1);
    {
      /* copy object color of previous object (if any) */
      ObjectMolecule* singleObj = nullptr;
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        at = I->Table[a].atom;
        I->Table[a].index = -1;
        obj = I->Obj[I->Table[a].model];
        s = obj->AtomInfo[at].selEntry;
        if (SelectorIsMember(G, s, sele)) {
          if (!singleObj)
            singleObj = obj;
          else if (singleObj && (obj != singleObj)) {
            singleObj = nullptr;
            break;
          }
        }
      }
      if (singleObj) targ->Color = singleObj->Color;
      /* should also consider copying lots of other stuff from the source object
       * ... */
    }
  } else {
    isNew = false;
  }

  std::function<void(int)> const body = [&](int const source) {
    c = 0;

    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      at = I->Table[a].atom;
      I->Table[a].index = -1;
      obj = I->Obj[I->Table[a].model];
      s = obj->AtomInfo[at].selEntry;
      if (SelectorIsMember(G, s, sele)) {
        I->Table[a].index = c; /* Mark records as to which atom they'll be */
        c++;
        if (!info_src) info_src = obj;
      }
    }
    if (isNew && info_src) { /* copy symmetry information, etc. */
      if (targ->Symmetry == nullptr && info_src->Symmetry != nullptr) {
        targ->Symmetry.reset(new CSymmetry(*info_src->Symmetry));
      }
    }

    if (info_src && source == cSelectorUpdateTableAllStates &&
        targ->DiscreteFlag && !info_src->DiscreteFlag) {
      for (int state = 0; state < info_src->getNFrame(); ++state) {
        body(state);
      }
      return;
    }

    nAtom = c;

    nBond = 0;
    auto bond = pymol::vla<BondType>(nAtom * 4);
    for (a = cNDummyModels; a < I->Obj.size();
         a++) { /* find bonds wholly contained in the selection */
      obj = I->Obj[a];
      ii1 = obj->Bond;
      for (b = 0; b < obj->NBond; b++) {
        b1 = SelectorGetObjAtmOffset(I, obj, ii1->index[0]);
        b2 = SelectorGetObjAtmOffset(I, obj, ii1->index[1]);
        if ((b1 >= 0) && (b2 >= 0)) {
          if ((I->Table[b1].index >= 0) && (I->Table[b2].index >= 0)) {
            BondType* dst_bond = bond.check(nBond);
            {
              AtomInfoBondCopy(G, ii1, dst_bond);
              dst_bond->index[0] =
                  I->Table[b1].index; /* store what will be the new index */
              dst_bond->index[1] = I->Table[b2].index;
              /*            printf("Selector-DEBUG %d
               * %d\n",dst_bond->index[0],dst_bond->index[1]); */
              nBond++;
            }
          }
        }
        ii1++;
      }
    }

    pymol::vla<AtomInfoType> atInfo(nAtom);
    /* copy the atom info records and create new zero-based IDs */
    c = 0;
    {
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        if (I->Table[a].index >= 0) {
          obj = I->Obj[I->Table[a].model];
          at = I->Table[a].atom;
          VLACheck(atInfo, AtomInfoType, c);
          AtomInfoCopy(G, obj->AtomInfo + at, atInfo + c);
          c++;
        }
      }
    }

    cs = CoordSetNew(G); /* set up a dummy coordinate set for the merge xref */
    cs->NIndex = nAtom;
    cs->enumIndices();
    cs->TmpBond = std::move(bond); /* load up the bonds */
    cs->NTmpBond = nBond;

    /*  printf("Selector-DEBUG nAtom %d\n",nAtom); */
    ObjectMoleculeMerge(targ, std::move(atInfo), cs, false, cAIC_AllMask,
                        true); /* will free atInfo */
    /* cs->IdxToAtm will now have the reverse mapping from the new subset
       to the new merged molecule */

    ObjectMoleculeExtendIndices(targ, -1);
    ObjectMoleculeUpdateIDNumbers(targ);
    ObjectMoleculeUpdateNonbonded(targ);

    if (!isNew) { /* recreate selection table */

      SelectorUpdateTable(G, source, -1);
    }

    /* get maximum state index for the selection...note that
       we'll be creating states from 1 up to the maximum required to
       capture atoms in the selection
     */

    nCSet = 0;

    {
      c = 0;
      for (a = cNDummyAtoms; a < I->Table.size(); a++) {
        at = I->Table[a].atom;
        I->Table[a].index = -1;
        obj = I->Obj[I->Table[a].model];
        s = obj->AtomInfo[at].selEntry;
        if (SelectorIsMember(G, s, sele)) {
          I->Table[a].index = c; /* Mark records  */
          if (nCSet < obj->NCSet) nCSet = obj->NCSet;
          c++;
        }
      }
    }

    if (c != nAtom) ErrFatal(G, "SelectorCreate", "inconsistent selection.");
    /* cs->IdxToAtm now has the relevant indexes for the coordinate transfer */

    for (StateIterator iter(G, nullptr, source, nCSet); iter.next();) {
      d = iter.state;
      {
        cs2 = CoordSetNew(G);
        c = 0;
        cs2->setNIndex(nAtom);
        for (a = cNDummyAtoms; a < I->Table.size();
             a++) /* any selected atoms in this state? */
          if (I->Table[a].index >= 0) {
            at = I->Table[a].atom;
            obj = I->Obj[I->Table[a].model];
            cs1 = nullptr;
            if (d < obj->NCSet) {
              cs1 = obj->CSet[d];
            } else if (singletons && (obj->NCSet == 1)) {
              cs1 = obj->CSet[0];
            }
            if (cs1) {
              if ((!cs2->Name[0]) &&
                  (cs1->Name[0])) /* copy the molecule name (if any) */
                strcpy(cs2->Name, cs1->Name);

              if (CoordSetGetAtomVertex(cs1, at, cs2->coordPtr(c))) {
                a2 = cs->IdxToAtm[I->Table[a].index]; /* actual merged atom
                                                         index */
                cs2->IdxToAtm[c] = a2;
                c++;
              }
            }
          }
        cs2->setNIndex(c);
        if (target >= 0) {
          ts = target++;
        } else {
          ts = d;
        }
        VLACheck(targ->CSet, CoordSet*, ts);
        if (targ->NCSet <= ts) targ->NCSet = ts + 1;
        delete targ->CSet[ts];
        targ->CSet[ts] = cs2;
        cs2->Obj = targ;
      }
    }
    delete cs;
  };  // end body

  body(source);

  targ->updateAtmToIdx();

  SceneCountFrames(G);
  if (!quiet) {
    PRINTFB(G, FB_Selector, FB_Details)
    " Selector: found %d atoms.\n", nAtom ENDFB(G);
  }
  if (ok) ok &= ObjectMoleculeSort(targ);
  if (isNew) {
    ObjectSetName(targ, name);
    ExecutiveManageObject(G, targ, zoom, quiet);
  } else {
    ExecutiveUpdateObjectSelection(G, targ);
  }
  SceneChanged(G);
  return ok;
}
