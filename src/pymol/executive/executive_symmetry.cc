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

// Crystal symmetry and symmetry expansion.

#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list.h"
#include "pymol/base/list_macros.h"
#include "pymol/objects/object_volume.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

/**
 * Get the list of objects which match the pattern, or all objects
 * if pattern equals 'all'.
 * @param str pattern string provided by user
 * @return List of borrowed object pointers
 */
static std::vector<pymol::CObject*> ExecutiveGetObjectsFromPattern(
    PyMOLGlobals* G, pymol::zstring_view pattern) {
  std::vector<pymol::CObject*> objects;

  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, pattern)) {
    switch (rec.type) {
      case cExecObject:
        objects.push_back(rec.obj);
        break;
      case cExecAll:
        for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
          if (rec.type == cExecObject) {
            objects.push_back(rec.obj);
          }
        }
    }
  }

  return objects;
}

/**
 * @param sele object name (or single-object atom selection expression)
 * @param state object state (only for maps, ignored for molecules)
 * @return true if symmetry is defined
 */
pymol::Result<bool> ExecutiveGetSymmetry(PyMOLGlobals* G, const char* sele,
                                         int state, float* a, float* b,
                                         float* c, float* alpha, float* beta,
                                         float* gamma, char* sgroup) {
  // try `sele` as an object name
  auto obj = ExecutiveFindObjectByName(G, sele);

  // odd feature: accept atom selections to select single object
  if (!obj) {
    auto tmpsele1 = SelectorTmp::make(G, sele);
    p_return_if_error(tmpsele1);
    obj = SelectorGetSingleObjectMolecule(G, tmpsele1->getIndex());
    if (!obj) {
      return pymol::make_error("selection must refer to exactly one object");
    }
  }

  CSymmetry const* symm = obj->getSymmetry(state);

  if (symm) {
    *a = symm->Crystal.dims()[0];
    *b = symm->Crystal.dims()[1];
    *c = symm->Crystal.dims()[2];
    *alpha = symm->Crystal.angles()[0];
    *beta = symm->Crystal.angles()[1];
    *gamma = symm->Crystal.angles()[2];
    UtilNCopy(sgroup, symm->spaceGroup(), sizeof(WordType));
    return true;
  }

  return false;
}

/**
 * Set symmetry for one or more objects.
 * @param names Object name pattern
 * @param state Object state (supports all=-1 and current=-2)
 * @return True if symmetry could be set for at least one object
 */
static bool ExecutiveSetSymmetry(PyMOLGlobals* G, pymol::zstring_view names,
                                 int const state, CSymmetry const& symmetry,
                                 bool const quiet = false) {
  auto objects = ExecutiveGetObjectsFromPattern(G, names);
  bool success = false;

  for (auto* obj : objects) {
    if (!obj->setSymmetry(symmetry, state)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " %s-Warning: Cannot set symmetry for '%s' (type %s)\n", __func__,
          obj->Name, typeid(obj).name() ENDFB(G);
      continue;
    }

    success = true;

    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Details)
      " %s-Details: Updated symmetry for '%s'\n", __func__, obj->Name ENDFB(G);
    }
  }

  return success;
}

pymol::Result<> ExecutiveSetSymmetry(PyMOLGlobals* G, const char* sele,
                                     int state, float a, float b, float c,
                                     float alpha, float beta, float gamma,
                                     const char* sgroup, int quiet) {
  /* create a new symmetry object for copying */
  CSymmetry symmetry(G);
  symmetry.Crystal.setDims(a, b, c);
  symmetry.Crystal.setAngles(alpha, beta, gamma);
  symmetry.setSpaceGroup(sgroup);

  if (!ExecutiveSetSymmetry(G, sele, state, symmetry, quiet)) {
    return pymol::make_error("no object selected");
  }

  return {};
}

pymol::Result<> ExecutiveSymmetryCopy(PyMOLGlobals* G, const char* source_name,
                                      const char* target_name, int source_state,
                                      int target_state, int quiet) {
  /* Copy the symmetry info from source to target; currently maps can have
   * multiple states for symmetry, but ObjectMolecule cannot */

  auto const* source_obj = ExecutiveFindObjectByName(G, source_name);
  if (!source_obj) {
    return pymol::make_error("source object not found");
  }

  auto const* source_symm = source_obj->getSymmetry(source_state);
  if (!source_symm) {
    return pymol::make_error(pymol::string_format(
        "no symmetry in object '%s' state %d", source_name, source_state));
  }

  if (!ExecutiveSetSymmetry(G, target_name, target_state, *source_symm,
                            quiet)) {
    return pymol::make_error("target object not found");
  }

  return {};
}

/**
 * Makes a custom symmetry operation label.
 */
static std::string make_symexp_segi_label(int a, int x, int y, int z) {
  SegIdent seg = "_000";

  if (a > 61) {
    // beyond what can be encoded with a single alphanumeric
    // character (PYMOL-2475)
  } else if (a > 35) {
    seg[0] = 'a' + (a - 36);
  } else if (a > 25) {
    seg[0] = '0' + (a - 26);
  } else {
    seg[0] = 'A' + a;
  }

  if (x > 0) {
    seg[1] = 'A' + x - 1;
  } else if (x < 0) {
    seg[1] = 'Z' + x + 1;
  }

  if (y > 0) {
    seg[2] = 'A' + y - 1;
  } else if (y < 0) {
    seg[2] = 'Z' + y + 1;
  }

  if (z > 0) {
    seg[3] = 'A' + z - 1;
  } else if (z < 0) {
    seg[3] = 'Z' + z + 1;
  }

  return seg;
}
/**
 * Generates symmetry mates of object @a name as new objects, with all states.
 *
 * Only includes mates within @a cutoff of selection @a s1, and only in a range
 * of +/-1 unit cells. The cutoff applies across states, so a symmetry atom will
 * be included if it's within any atom of @a s1 in any state.
 *
 * @param name Prefix for new objects
 * @param oname Object to replicate
 * @param s1 Atom selection (for cutoff)
 * @param cutoff Distance cutoff from selection
 * @param segi (bool) If true, write symmetry operation code to segment
 * identifier
 */
void ExecutiveSymExp(PyMOLGlobals* G, const char* name, const char* oname,
                     const char* s1, float cutoff, int segi,
                     int quiet) { /* TODO state */
  SelectorTmp tmpsele1(G, s1);
  auto sele = tmpsele1.getIndex();

  /* object to expand */
  auto const* const obj = ExecutiveFindObject<ObjectMolecule>(G, oname);

  if (!(obj && sele)) {
    ErrMessage(G, __func__, "Invalid object");
    return;
  }

  // Get symmetry from first coordset
  const CSymmetry* sym = nullptr;
  for (int b = 0; b < obj->NCSet && !sym; ++b) {
    if (auto const* cs = obj->CSet[b]) {
      sym = cs->getSymmetry();
    }
  }

  if (!sym) {
    ErrMessage(G, __func__, "No symmetry loaded!");
    return;
  }

  int const nsymmat = sym->getNSymMat();

  if (nsymmat < 1) {
    ErrMessage(G, __func__, "unknown space group");
    return;
  }

  if (!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
    " ExecutiveSymExp: Generating symmetry mates...\n" ENDFB(G);
  }

  bool const matrix_mode = SettingGet<int>(G, cSetting_matrix_mode) > 0;

  // Get the center of mass for this object/selection (over all states!)
  ObjectMoleculeOpRec op;
  ObjectMoleculeOpRecInit(&op);
  op.code = OMOP_SUMC;
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto tc_vec = glm::make_vec3(op.v1) / float(std::max(1, op.i1));
  auto tc = glm::value_ptr(tc_vec);

  /* Transformation: RealToFrac (3x3) * tc (3x1) = (3x1) */
  transform33f3f(sym->Crystal.realToFrac(), tc, tc);

  /* 2.  Copy the coordinates for the atoms in this selection into op */
  op.code = OMOP_VERT;
  op.nvv1 = 0;
  op.vv1 = VLAlloc(float, 10000);
  ExecutiveObjMolSeleOp(G, sele, &op);
  auto const vv1 = pymol::vla_take_ownership(op.vv1);

  /* op.nvv1 is the number of atom coordinates we copied in the previous step */
  if (!op.nvv1) {
    ErrMessage(G, __func__, "No atoms indicated!");
    return;
  }

  auto map = std::unique_ptr<MapType>(
      new MapType(G, -cutoff, vv1.data(), op.nvv1, nullptr));
  if (!map) {
    ErrMessage(G, __func__, "map is nullptr");
    return;
  }

  std::vector<glm::vec3> cs_centers(obj->NCSet);
  for (int b = 0; b < obj->NCSet; ++b) {
    if (auto const* cs = obj->CSet[b]) {
      CoordSetGetAverage(cs, glm::value_ptr(cs_centers[b]));
    }
  }

  /* go out no more than one lattice step in each direction: -1, 0, +1 */
  for (int x = -1; x < 2; ++x) {
    for (int y = -1; y < 2; ++y) {
      for (int z = -1; z < 2; ++z) {
        for (int a = 0; a < nsymmat; a++) {
          /* make a copy of the original */
          auto new_obj = ObjectMoleculeCopy(obj);
          bool keepFlag = false;

          for (int b = 0; b < new_obj->NCSet; ++b) {
            auto* cs = new_obj->CSet[b];
            if (!cs) {
              continue;
            }

            float mat[16];
            copy33f44f(sym->Crystal.realToFrac(), mat);
            left_multiply44f44f(sym->getSymMat(a), mat);

            float ts[3];
            transform44f3f(mat, glm::value_ptr(cs_centers[b]), ts);

            /* compute the effective translation resulting
               from application of the symmetry operator so
               that we can shift it into the cell of the
               target selection */
            for (int c = 0; c < 3; c++) {
              ts[c] = std::round(tc[c] - ts[c]);
            }
            float const shift[3] = {ts[0] + x, ts[1] + y, ts[2] + z};
            float m[16];
            identity44f(m);
            m[3] = shift[0];
            m[7] = shift[1];
            m[11] = shift[2];

            left_multiply44f44f(const_cast<float const*>(m), mat);
            copy33f44f(sym->Crystal.fracToReal(), m);
            left_multiply44f44f(const_cast<float const*>(m), mat);

            if (is_identityf(4, mat)) {
              continue;
            }

            double mat_d[16];
            copy44f44d(mat, mat_d);
            ObjectStateLeftCombineMatrixR44d(cs, mat_d);

            if (!matrix_mode) {
              CoordSetTransform44f(cs, mat);
            }

            if (keepFlag) {
              continue;
            }

            /* for each coordinate in this coordinate set */
            for (unsigned idx = 0; idx < cs->NIndex; ++idx) {
              const auto* v2 = cs->coordPtr(idx);

              if (matrix_mode) {
                transform44f3f(mat, v2, ts);
                v2 = ts;
              }

              if (MapAnyWithin(*map, vv1.data(), v2, cutoff)) {
                keepFlag = true;
                break;
              }
            }

            // CIF-style symop label (e.g. "1_555")
            cs->setTitle(pymol::string_format("%d_%.0f%.0f%.0f", a + 1,
                                              shift[0] + 5, shift[1] + 5,
                                              shift[2] + 5));
          }

          if (!keepFlag) {
            DeleteP(new_obj);
            continue;
          }

          if (segi) {
            auto seg = make_symexp_segi_label(a, x, y, z);
            lexidx_t segi = LexIdx(G, seg.c_str());
            for (unsigned atm = 0; atm < new_obj->NAtom; ++atm) {
              LexAssign(G, new_obj->AtomInfo[atm].segi, segi);
            }
            LexDec(G, segi);
          }

          // manage the new object

          auto const new_name =
              pymol::string_format("%s%02d%02d%02d%02d", name, a, x, y, z);

          PRINTFB(G, FB_Executive, FB_Blather)
          "Making new object: %s from name=%s, a=%d, x=%d, y=%d, z=%d\n",
              new_name.c_str(), name, a, x, y, z ENDFB(G);

          ObjectSetName(new_obj, new_name.c_str());
          ExecutiveDelete(G, new_obj->Name);
          ExecutiveManageObject(G, new_obj, false, quiet);
        }
      }
    }
  }
}
