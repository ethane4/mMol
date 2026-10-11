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

// Reading molecules from PDB, MOL, SDF, MOL2, XYZ, MacroModel and ChemPy
// models.

#include "pymol/base/parse.h"
#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/python/p.h"
#include "pymol/executive/executive.h"
#include "pymol/io/mol_v3000.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/property.h"
#include "pymol/objects/object_molecule_internal.h"

#define nskip ParseNSkip
static CoordSet* ObjectMoleculeXYZStr2CoordSet(PyMOLGlobals* G,
                                               const char* buffer,
                                               AtomInfoType** atInfoPtr,
                                               const char** restart) {
  const char *p, *p_store;
  int nAtom;
  int a, c;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType *atInfo = nullptr, *ai;
  char cc[MAXLINELEN];
  int atomCount;
  BondType* bond = nullptr;
  int nBond = 0;
  int b1, b2;
  WordType tmp_name;
  int auto_show = RepGetAutoShowMask(G);
  int tinker_xyz = true;
  int valid_atom;
  int have_n_atom = false;
  BondType* ii;

  p = buffer;
  nAtom = 0;
  atInfo = *atInfoPtr;

  p_store = p;
  p = ncopy(cc, p, MAXLINELEN - 1);
  if (sscanf(cc, "%d", &nAtom) != 1) {
    nAtom = 0;
    tinker_xyz = false;
    p = p_store;
  } else {
    have_n_atom = true;
    p = nskip(p, 2);
    p = ncopy(tmp_name, p, sizeof(WordType) - 1);
    p = nextline(p);
  }

  if (tinker_xyz && nAtom) { /* test Tinker XYZ formatting assumption */
    const char* pp = p;
    int dummy_int;
    float dummy_float;
    AtomName dummy_name;

    pp = ncopy(cc, pp, 6);
    if (!sscanf(cc, "%d", &dummy_int)) tinker_xyz = false; /* id */
    pp = nskip(pp, 2);
    pp = ncopy(cc, pp, 3);
    if (sscanf(cc, "%s", dummy_name) != 1) tinker_xyz = false; /* name */
    pp = ncopy(cc, pp, 12);
    if (sscanf(cc, "%f", &dummy_float) != 1) tinker_xyz = false; /* x */
    pp = ncopy(cc, pp, 12);
    if (sscanf(cc, "%f", &dummy_float) != 1) tinker_xyz = false; /* y */
    pp = ncopy(cc, pp, 12);
    if (sscanf(cc, "%f", &dummy_float) != 1) tinker_xyz = false; /* z */
    pp = ncopy(cc, pp, 6);
    if (sscanf(cc, "%d", &dummy_int) != 1)
      tinker_xyz = false; /* numeric type */
  }

  if (!tinker_xyz) {
    const char* pp = p;
    int have_atom_line = true;
    float dummy_float;
    AtomName dummy_name;

    pp = ParseWordCopy(cc, pp, sizeof(AtomName) - 1);
    if (sscanf(cc, "%s", dummy_name) != 1) have_atom_line = false; /* name */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if (sscanf(cc, "%f", &dummy_float) != 1) have_atom_line = false; /* x */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if (sscanf(cc, "%f", &dummy_float) != 1) have_atom_line = false; /* y */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if (sscanf(cc, "%f", &dummy_float) != 1) have_atom_line = false; /* z */
    if (!have_atom_line) { /* copy the comment line into the title field */
      p = ncopy(tmp_name, p, sizeof(WordType) - 1);
      p = nextline(p);
    }
  }

  if (nAtom) {
    coord = VLAlloc(float, 3 * nAtom);
    if (atInfo) VLACheck(atInfo, AtomInfoType, nAtom);
  } else {
    coord = VLAlloc(float, 3);
  }

  if (tinker_xyz) {
    nBond = 0;
    bond = VLACalloc(BondType, 6 * nAtom); /* is this a safe assumption? */
    ii = bond;
  }

  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
  " ObjectMoleculeReadXYZ: Found %i atoms...\n", nAtom ENDFB(G);

  a = 0;
  atomCount = 0;
  while (*p) {
    VLACheck(atInfo, AtomInfoType, atomCount);
    ai = atInfo + atomCount;

    if (!tinker_xyz) {
      valid_atom = true;

      p = ParseWordCopy(cc, p, MAXLINELEN - 1);
      UtilCleanStr(cc);
      if (!cc[0]) valid_atom = false;
      if (valid_atom) {
        strncpy(ai->elem, cc, cElemNameLen);
        ai->name = LexIdx(G, cc);

        ai->rank = atomCount;
        ai->id = atomCount + 1;

        VLACheck(coord, float, a * 3 + 2);
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if (sscanf(cc, "%f", coord + a) != 1) valid_atom = false;
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if (sscanf(cc, "%f", coord + a + 1) != 1) valid_atom = false;
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if (sscanf(cc, "%f", coord + a + 2) != 1) valid_atom = false;

        ai->resn = LexIdx(G, "UNK");

        ai->alt[0] = 0;
        ai->chain = 0;
        ai->resv = atomCount + 1;

        ai->q = 1.0;
        ai->b = 0.0;

        ai->segi = 0;

        ai->visRep = auto_show;

        /* in the absense of external tinker information, assume hetatm */

        ai->hetatm = 1;

        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }
    } else { /* tinker XYZ */

      valid_atom = true;

      p = ncopy(cc, p, 6);
      if (!sscanf(cc, "%d", &ai->id)) break;
      ai->rank = atomCount;

      p = nskip(p, 2); /* to 12 */
      p = ntrim(cc, p, 3);
      ai->name = LexIdx(G, cc);

      ai->resn = LexIdx(G, "UNK");

      ai->alt[0] = 0;
      ai->chain = 0;

      ai->resv = atomCount + 1;

      valid_atom = true;

      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + a);
      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + (a + 1));
      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + (a + 2));

      ai->q = 1.0;
      ai->b = 0.0;

      ai->segi = 0;
      ai->elem[0] = 0; /* let atom info guess/infer atom type */

      ai->visRep = auto_show;

      p = ncopy(cc, p, 6);
      sscanf(cc, "%d", &ai->customType);

      /* in the absense of external tinker information, assume hetatm */

      ai->hetatm = 1;

      AtomInfoAssignParameters(G, ai);
      AtomInfoAssignColors(G, ai);

      b1 = atomCount;
      for (c = 0; c < 6; c++) {
        p = ncopy(cc, p, 6);
        if (!cc[0]) break;
        if (!sscanf(cc, "%d", &b2)) break;
        if (b1 < (b2 - 1)) {
          VLACheck(bond, BondType, nBond);
          ii = bond + nBond;
          nBond++;
          ii->index[0] = b1;
          ii->index[1] = b2 - 1;
          ii->order = 1; /* missing bond order information */
          ii++;
        }
      }
    }

    if (valid_atom) {
      PRINTFD(G, FB_ObjectMolecule)
      " ObjectMolecule-DEBUG: %s %s %d %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n",
          LexStr(G, ai->name), LexStr(G, ai->resn), ai->resv,
          LexStr(G, ai->chain), *(coord + a), *(coord + a + 1),
          *(coord + a + 2), ai->b, ai->q, LexStr(G, ai->segi) ENDFD;

      a += 3;
      atomCount++;
    }
    p = nextline(p);
    if (have_n_atom && (atomCount >= nAtom)) {
      int dummy;
      ncopy(cc, p, MAXLINELEN - 1);
      if (sscanf(cc, "%d", &dummy) == 1) *restart = p;
      break;
    }
  }

  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
  " XYZStr2CoordSet: Read %d bonds.\n", nBond ENDFB(G);

  if (!tinker_xyz) nAtom = atomCount; /* use number of atoms actually read */

  cset = CoordSetNew(G);
  cset->NIndex = nAtom;
  cset->Coord = pymol::vla_take_ownership(coord);
  cset->TmpBond = pymol::vla_take_ownership(bond);
  cset->NTmpBond = nBond;
  strcpy(cset->Name, tmp_name);
  if (atInfoPtr) *atInfoPtr = atInfo;
  return (cset);
}
static void copy_and_trim_repr(char* dest, PyObject* obj) {
  PyObject* repr = PyObject_Repr(obj);
  const char* src = PyString_AS_STRING(obj);
  int srclen;
  if (src[0] == '\'') src++;
  srclen = strlen(src);
  if (src[srclen - 1] == '\'') srclen -= 1;
  strncpy(dest, src, srclen);
  dest[srclen] = 0;
  Py_DECREF(repr);
}

static CoordSet* ObjectMoleculeChemPyModel2CoordSet(PyMOLGlobals* G,
                                                    PyObject* model,
                                                    AtomInfoType** atInfoPtr) {
  int nAtom, nBond;
  int a, c;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType *atInfo = nullptr, *ai;
  float* f;
  BondType *ii, *bond = nullptr;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  int hetatm;
  int ignore_ids;
  PyObject* atomList = nullptr;
  PyObject* bondList = nullptr;
  PyObject* atom = nullptr;
  PyObject* bnd = nullptr;
  PyObject* index = nullptr;
  PyObject* crd = nullptr;
  PyObject* tmp = nullptr;

  ignore_ids = !SettingGetGlobal_b(G, cSetting_preserve_chempy_ids);

  nAtom = 0;
  nBond = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;

  atomList = PyObject_GetAttrString(model, "atom");
  if (atomList && PyList_Check(atomList))
    nAtom = PyList_Size(atomList);
  else
    ok = ErrMessage(G, __func__, "can't get atom list");

  if (ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if (atInfo) VLACheck(atInfo, AtomInfoType, nAtom);
  }

  if (ok) {
    f = coord;
    for (a = 0; a < nAtom; a++) {
      atom = PyList_GetItem(atomList, a);
      if (!atom) ok = ErrMessage(G, __func__, "can't get atom");
      crd = PyObject_GetAttrString(atom, "coord");
      if (!crd)
        ok = ErrMessage(G, __func__, "can't get coordinates");
      else {
        for (c = 0; c < 3; c++) {
          tmp = PySequence_GetItem(crd, c);
          if (tmp) ok = PConvPyObjectToFloat(tmp, f++);
          Py_XDECREF(tmp);
          if (!ok) {
            ErrMessage(G, __func__, "can't read coordinates");
            break;
          }
        }
      }
      Py_XDECREF(crd);

      ai = atInfo + a;
      ai->id = a; /* chempy models are zero-based */
      if (!ignore_ids) {
        if (ok) { /* get chempy atom id if extant */
          if (PTruthCallStr(atom, "has", "id")) {
            tmp = PyObject_GetAttrString(atom, "id");
            if (tmp) ok = PConvPyObjectToInt(tmp, &ai->id);
            if (!ok) ErrMessage(G, __func__, "can't read atom identifier");
            Py_XDECREF(tmp);
          } else {
            ai->id = -1;
          }
        }
      }
      ai->rank = a; /* ranks are always zero-based */

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "name");
        if (tmp) {
          WordType tmp_word;

          ok = PConvPyObjectToStrMaxClean(tmp, tmp_word, sizeof(WordType) - 1);
          AtomInfoCleanAtomName(tmp_word);
          ai->name = LexIdx(G, tmp_word);
        }
        if (!ok) ErrMessage(G, __func__, "can't read name");
        Py_XDECREF(tmp);
      }

      if (ok) {
        ai->textType = 0;
        if (PTruthCallStr(atom, "has", "text_type")) {
          tmp = PyObject_GetAttrString(atom, "text_type");
          if (tmp) {
            OrthoLineType temp;
            ok = PConvPyObjectToStrMaxClean(tmp, temp,
                                            sizeof(OrthoLineType) - 1);
            ai->textType = LexIdx(G, temp);
          }
          if (!ok) ErrMessage(G, __func__, "can't read text_type");
          Py_XDECREF(tmp);
        }
      }

      if (ok) {
        ai->custom = 0;
        if (PTruthCallStr(atom, "has", "custom")) {
          tmp = PyObject_GetAttrString(atom, "custom");
          if (tmp) {
            OrthoLineType temp;
            ok = PConvPyObjectToStrMaxClean(tmp, temp,
                                            sizeof(OrthoLineType) - 1);
            ai->custom = LexIdx(G, temp);
          }
          if (!ok) ErrMessage(G, __func__, "can't read custom");
          Py_XDECREF(tmp);
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "vdw")) {
          tmp = PyObject_GetAttrString(atom, "vdw");
          if (tmp) ok = PConvPyObjectToFloat(tmp, &ai->vdw);
          if (!ok) ErrMessage(G, __func__, "can't read vdw radius");
          Py_XDECREF(tmp);
        } else {
          ai->vdw = 0.0f;
        }
      }
      if (ok) {
        if (PTruthCallStr(atom, "has", "bohr")) {
          tmp = PyObject_GetAttrString(atom, "bohr");
          if (tmp) ok = PConvPyObjectToFloat(tmp, &ai->elec_radius);
          if (!ok) ErrMessage(G, __func__, "can't read elec. radius");
          Py_XDECREF(tmp);
        } else {
          ai->elec_radius = 0.0F;
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "stereo")) {
          tmp = PyObject_GetAttrString(atom, "stereo");
          if (tmp) {
            char tmp_char;
            ok = PConvPyObjectToChar(tmp, &tmp_char);
            ai->stereo = tmp_char;
          }
          if (!ok) ErrMessage(G, __func__, "can't read stereo");
          Py_XDECREF(tmp);
        } else {
          ai->stereo = 0;
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "numeric_type")) {
          tmp = PyObject_GetAttrString(atom, "numeric_type");
          if (tmp) ok = PConvPyObjectToInt(tmp, &ai->customType);
          if (!ok) ErrMessage(G, __func__, "can't read numeric_type");
          Py_XDECREF(tmp);
        } else {
          ai->customType = cAtomInfoNoType;
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "formal_charge")) {
          tmp = PyObject_GetAttrString(atom, "formal_charge");
          if (tmp) ok = PConvPyObjectToChar(tmp, (char*)&ai->formalCharge);
          if (!ok) ErrMessage(G, __func__, "can't read formal_charge");
          Py_XDECREF(tmp);
        } else {
          ai->formalCharge = 0;
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "partial_charge")) {
          tmp = PyObject_GetAttrString(atom, "partial_charge");
          if (tmp) ok = PConvPyObjectToFloat(tmp, &ai->partialCharge);
          if (!ok) ErrMessage(G, __func__, "can't read partial_charge");
          Py_XDECREF(tmp);
        } else {
          ai->partialCharge = 0.0;
        }
      }

      if (ok) {
        if (PTruthCallStr(atom, "has", "flags")) {
          tmp = PyObject_GetAttrString(atom, "flags");
          if (tmp) ok = PConvPyObjectToInt(tmp, (int*)&ai->flags);
          if (!ok) ErrMessage(G, __func__, "can't read flags");
          Py_XDECREF(tmp);
        } else {
          ai->flags = 0;
        }
      }

      if (ok) {
        char buf[8] = "";
        tmp = PyObject_GetAttrString(atom, "resn");
        if (tmp) ok = PConvPyObjectToStrMaxClean(tmp, buf, sizeof(buf) - 1);
        ai->resn = LexIdx(G, buf);
        if (!ok) ErrMessage(G, __func__, "can't read resn");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "ins_code");
        if (tmp) {
          ResIdent tmp_ins_code;
          ok = PConvPyObjectToStrMaxClean(tmp, tmp_ins_code,
                                          sizeof(ResIdent) - 1);
          if (!ok)
            ErrMessage(G, __func__, "can't read ins_code");
          else if (tmp_ins_code[0] != '?') {
            ai->setInscode(tmp_ins_code[0]);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok) {
        if (PTruthCallStr(atom, "has", "resi_number")) {
          tmp = PyObject_GetAttrString(atom, "resi_number");
          if (tmp) ok = PConvPyObjectToInt(tmp, &ai->resv);
          if (!ok) ErrMessage(G, __func__, "can't read resi_number");
          Py_XDECREF(tmp);
        }
      }

      if (ok) {
        OrthoLineType temp = "";
        tmp = PyObject_GetAttrString(atom, "segi");
        if (tmp)
          PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
        ai->segi = LexIdx(G, temp);
        if (!ok) ErrMessage(G, __func__, "can't read segi");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "b");
        if (tmp) ok = PConvPyObjectToFloat(tmp, &ai->b);
        if (!ok) ErrMessage(G, __func__, "can't read b value");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "u");
        if (tmp) {
          ok = PConvPyObjectToFloat(tmp, &ai->b);
          if (!ok)
            ErrMessage(G, __func__, "can't read u value");
          else
            ai->b *= (8 * cPI * cPI); /* B-factor = 8 pi^2 U-factor */
        } else {
          assert(PyErr_Occurred() &&
                 PyErr_ExceptionMatches(PyExc_AttributeError));
          PyErr_Clear();
        }
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "u_aniso");
        if (tmp) {
          float u[6];
          if (PConvPyListToFloatArrayInPlace(tmp, u, 6)) {
            // only allocate if not all zero
            if (u[0] || u[1] || u[2] || u[3] || u[4] || u[5])
              std::copy_n(u, 6, ai->get_anisou());
          }
          Py_DECREF(tmp);
        } else {
          assert(PyErr_Occurred() &&
                 PyErr_ExceptionMatches(PyExc_AttributeError));
          PyErr_Clear();
        }
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "q");
        if (tmp) ok = PConvPyObjectToFloat(tmp, &ai->q);
        if (!ok) ErrMessage(G, __func__, "can't read occupancy");
        Py_XDECREF(tmp);
      }

      if (ok) {
        OrthoLineType temp = "";
        tmp = PyObject_GetAttrString(atom, "chain");
        if (tmp)
          PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
        ai->chain = LexIdx(G, temp);
        if (!ok) ErrMessage(G, __func__, "can't read chain");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "hetatm");
        if (tmp) ok = PConvPyObjectToInt(tmp, &hetatm);
        if (!ok)
          ErrMessage(G, __func__, "can't read hetatm");
        else {
          ai->hetatm = hetatm;
          if (!PTruthCallStr(atom, "has", "flags")) {
            if (ai->hetatm) ai->flags = cAtomFlag_ignore;
          }
        }
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "alt");
        if (tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->alt, sizeof(Chain) - 1);
        if (!ok) ErrMessage(G, __func__, "can't read alternate conformation");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "symbol");
        if (tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->elem, sizeof(ElemName) - 1);
        if (!ok) ErrMessage(G, __func__, "can't read symbol");
        Py_XDECREF(tmp);
      }

      if (ok) {
        tmp = PyObject_GetAttrString(atom, "ss");
        if (tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->ssType, sizeof(SSType) - 1);
        if (!ok) ErrMessage(G, __func__, "can't read secondary structure");
        Py_XDECREF(tmp);
      }

      if (ok && PyObject_HasAttrString(atom, "label")) {
        tmp = PyObject_GetAttrString(atom, "label");
        if (tmp) {
          if (PyString_Check(tmp))
            ai->label = LexIdx(G, PyString_AsString(tmp));
        }
        Py_XDECREF(tmp);
      }

      if (ok && PyObject_HasAttrString(atom, "sphere_scale")) {
        tmp = PyObject_GetAttrString(atom, "sphere_scale");
        if (tmp) {
          float value;
          if (PConvPyFloatToFloat(tmp, &value)) {
            SettingSet(G, cSetting_sphere_scale, value, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok && PyObject_HasAttrString(atom, "cartoon_color")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_color");
        if (tmp) {
          int color_index;
          ok = PConvPyObjectToInt(tmp, &color_index);
          if (!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            SettingSet(G, cSetting_cartoon_color, color_index, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok && PyObject_HasAttrString(atom, "cartoon_transparency")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_transparency");
        if (tmp) {
          float alpha_val;
          ok = PConvPyObjectToFloat(tmp, &alpha_val);
          if (!ok)
            ErrMessage(G, __FUNCTION__, "bad alpha value");
          else {
            SettingSet(G, cSetting_cartoon_transparency, alpha_val, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok && PyObject_HasAttrString(atom, "cartoon_trgb")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_trgb");
        if (tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int*)&trgb);
          if (!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_cartoon_color, ColorGetIndex(G, color_name),
                       ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok && PyObject_HasAttrString(atom, "label_trgb")) {
        tmp = PyObject_GetAttrString(atom, "label_trgb");
        if (tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int*)&trgb);
          if (!ok)
            ErrMessage(G, __func__, "bad label color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_label_color, ColorGetIndex(G, color_name),
                       ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if (ok && PyObject_HasAttrString(atom, "ribbon_color")) {
        tmp = PyObject_GetAttrString(atom, "ribbon_color");
        if (tmp) {
          int color_index;
          ok = PConvPyObjectToInt(tmp, &color_index);
          if (!ok)
            ErrMessage(G, __func__, "bad ribbon color info");
          else {
            SettingSet(G, cSetting_ribbon_color, color_index, ai);
          }
        }
        Py_XDECREF(tmp);
      }

      if (ok && PyObject_HasAttrString(atom, "ribbon_trgb")) {
        tmp = PyObject_GetAttrString(atom, "ribbon_trgb");
        if (tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int*)&trgb);
          if (!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_ribbon_color, ColorGetIndex(G, color_name),
                       ai);
          }
        }
        Py_XDECREF(tmp);
      }

      atInfo[a].visRep = auto_show;

      if (ok && PyObject_HasAttrString(atom, "visible")) {
        unsigned int vis;
        tmp = PyObject_GetAttrString(atom, "visible");
        if (tmp) {
          ok = PConvPyObjectToInt(tmp, (signed int*)&vis);
          if (!ok)
            ErrMessage(G, __func__, "bad visibility info");
          else {
            atInfo[a].visRep = vis;

            if (!(ai->flags & cAtomFlag_class)) {
              ai->flags |=
                  cAtomFlag_inorganic;  // suppress auto_show_classified
            }
          }
        }
        Py_XDECREF(tmp);
      }

      if (ok && atInfo) {
        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }

      if (ok && PyObject_HasAttrString(atom, "trgb")) {
        /* were we given a Transparency-Red-Blue-Green value? */
        unsigned int trgb;
        tmp = PyObject_GetAttrString(atom, "trgb");
        if (tmp) {
          ok = PConvPyObjectToInt(tmp, (signed int*)&trgb);
          if (!ok)
            ErrMessage(G, __func__, "bad color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            atInfo[a].color = ColorGetIndex(G, color_name);
          }
          Py_DECREF(tmp);
        }
      }

      if (!ok) break;

      assert(!PyErr_Occurred());
    }
  }

  if (nAtom) {
    bondList = PyObject_GetAttrString(model, "bond");
    if (bondList && PyList_Check(bondList))
      nBond = PyList_Size(bondList);
    else
      ok = ErrMessage(G, __func__, "can't get bond list");

    if (ok) {
      bond = VLACalloc(BondType, nBond);
      ii = bond;
      for (a = 0; a < nBond; a++) {
        bnd = PyList_GetItem(bondList, a);
        if (!bnd) ok = ErrMessage(G, __func__, "can't get bond");
        index = PyObject_GetAttrString(bnd, "index");
        if (!index)
          ok = ErrMessage(G, __func__, "can't get bond indices");
        else {
          for (c = 0; c < 2; c++) {
            tmp = PyList_GetItem(index, c);
            if (tmp) ok = PConvPyObjectToInt(tmp, &ii->index[c]);
            if (!ok) {
              ErrMessage(G, __func__, "can't read coordinates");
              break;
            }
          }
        }
        if (ok) {
          int order = 0;
          tmp = PyObject_GetAttrString(bnd, "order");
          if (tmp) ok = PConvPyObjectToInt(tmp, &order);
          if (!ok) ErrMessage(G, __func__, "can't read bond order");
          Py_XDECREF(tmp);
          ii->order = order;
        }

        auto const set_symop = [&](const char* key) {
          if (ok && PyObject_HasAttrString(bnd, key)) {
            ok = false;
            if (auto tmp = PyObject_GetAttrString(bnd, key)) {
              if (auto const* code = PyUnicode_AsUTF8(tmp)) {
                ii->symop_2.reset(code);
                ok = true;
              }
              Py_DECREF(tmp);
            }
          }
        };

        set_symop("symmetry_2");

        if (ok && PyObject_HasAttrString(bnd, "stick_radius")) {
          tmp = PyObject_GetAttrString(bnd, "stick_radius");
          if (tmp) {
            float value;
            if (PConvPyFloatToFloat(tmp, &value)) {
              SettingSet(G, cSetting_stick_radius, value, ii);
            }
          }
          Py_XDECREF(tmp);
        }
        if (ok && PyObject_HasAttrString(bnd, "valence")) {
          tmp = PyObject_GetAttrString(bnd, "valence");
          if (tmp) {
            int value;
            if (PConvPyIntToInt(tmp, &value)) {
              SettingSet(G, cSetting_valence, value, ii);
            }
          }
          Py_XDECREF(tmp);
        }

        Py_XDECREF(index);
        ii++;
      }
    }
  }

  Py_XDECREF(atomList);
  Py_XDECREF(bondList);

  if (ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);

    {
      PyObject *propList = nullptr, *propName, *propValue, *propKeyValue;
      OrthoLineType nameString;
      int nProps = 0;

      // need to check since fragments and other molecules use unpickle
      // and pml files might not have molecule_properties
      if (PyObject_HasAttrString(model, "molecule_properties")) {
        propList = PyObject_GetAttrString(model, "molecule_properties");
        if (propList && PyList_Check(propList))
          nProps = PyList_Size(propList);
        else
          ok = ErrMessage(G, __func__, "can't get molecule properties");
      }
      if (nProps) {
        for (a = 0; a < nProps; a++) {
          propKeyValue = PyList_GetItem(propList, a);
          propName = PyList_GetItem(propKeyValue, 0);
          propValue = PyList_GetItem(propKeyValue, 1);
          copy_and_trim_repr(nameString, propName);
          PropertySet(G, cset, nameString, propValue);
        }
      }
      Py_XDECREF(propList);
    }
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if (atInfoPtr) *atInfoPtr = atInfo;

  if (PyErr_Occurred()) PyErr_Print();
  return (cset);
}
ObjectMolecule* ObjectMoleculeLoadChemPyModel(PyMOLGlobals* G,
                                              ObjectMolecule* I,
                                              PyObject* model, int frame,
                                              int discrete) {
  CoordSet* cset = nullptr;
  auto atInfo = pymol::vla<AtomInfoType>(10);
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;
  int fractional = false;
  int connect_mode = -1;
  int auto_bond = false;
  PyObject *tmp, *mol;

  if (!I)
    isNew = true;
  else
    isNew = false;

  if (ok) {
    if (isNew) {
      I = (ObjectMolecule*)new ObjectMolecule(G, discrete);
      std::swap(atInfo, I->AtomInfo);
      isNew = true;
    } else {
      if (discrete > -1) ObjectMoleculeSetDiscrete(G, I, discrete);
    }

    if (isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }
    cset = ObjectMoleculeChemPyModel2CoordSet(G, model, &atInfo);

    if (!cset)
      ok = false;
    else {
      mol = PyObject_GetAttrString(model, "molecule");
      if (mol) {
        if (PyObject_HasAttrString(mol, "title")) {
          tmp = PyObject_GetAttrString(mol, "title");
          if (tmp) {
            UtilNCopy(cset->Name, PyString_AsString(tmp), sizeof(WordType));
            Py_DECREF(tmp);
            if (!strcmp(cset->Name, "untitled")) /* ignore untitled */
              cset->Name[0] = 0;
          }
        }
        Py_DECREF(mol);
      } else if (PyErr_Occurred() == PyExc_AttributeError) {
        PyErr_Clear();
      }
      if (PyObject_HasAttrString(model, "spheroid") &&
          PyObject_HasAttrString(model, "spheroid_normals")) {
        tmp = PyObject_GetAttrString(model, "spheroid");
        if (tmp) {
          PConvFromPyObject(G, tmp, cset->Spheroid);
          Py_DECREF(tmp);
        }
        tmp = PyObject_GetAttrString(model, "spheroid_normals");
        if (tmp) {
          PConvFromPyObject(G, tmp, cset->SpheroidNormal);
          Py_DECREF(tmp);
        }
      }
      if (PyObject_HasAttrString(model, "spacegroup") &&
          PyObject_HasAttrString(model, "cell")) {
        CSymmetry* symmetry = new CSymmetry(G);
        if (symmetry) {
          tmp = PyObject_GetAttrString(model, "spacegroup");
          if (tmp) {
            const char* tmp_str = nullptr;
            if (PConvPyStrToStrPtr(tmp, &tmp_str)) {
              symmetry->setSpaceGroup(tmp_str);
            }
            Py_DECREF(tmp);
          }
          tmp = PyObject_GetAttrString(model, "cell");
          if (tmp) {
            float cell[6];
            if (PConvPyListToFloatArrayInPlace(tmp, cell, 6)) {
              symmetry->Crystal.setDims(cell);
              symmetry->Crystal.setAngles(cell + 3);
            }
            Py_DECREF(tmp);
          }
          cset->Symmetry = std::unique_ptr<CSymmetry>(symmetry);
        }
      }
      if (PyObject_HasAttrString(model, "fractional")) {
        tmp = PyObject_GetAttrString(model, "fractional");
        if (tmp) {
          int tmp_int = 0;
          if (PConvPyIntToInt(tmp, &tmp_int)) {
            fractional = tmp_int;
          }
          Py_DECREF(tmp);
        }
      }
      if (PyObject_HasAttrString(model, "connect_mode")) {
        tmp = PyObject_GetAttrString(model, "connect_mode");
        if (tmp) {
          int tmp_int = 0;
          if (PConvPyIntToInt(tmp, &tmp_int)) {
            auto_bond = true;
            connect_mode = tmp_int;
          }
          Py_DECREF(tmp);
        }
      }
      nAtom = cset->NIndex;
    }
  }
  /* include coordinate set */
  if (ok) {
    if (frame < 0) frame = I->NCSet;

    if (I->DiscreteFlag && atInfo) {
      unsigned int a;
      int fp1 = frame + 1;
      AtomInfoType* ai = atInfo.data();
      for (a = 0; a < nAtom; a++) {
        (ai++)->discrete_state = fp1;
      }
    }

    cset->Obj = I;
    cset->enumIndices();
    cset->invalidateRep(cRepAll, cRepInvRep);
    if (isNew) {
      std::swap(I->AtomInfo, atInfo);
    } else {
      ObjectMoleculeMerge(I, std::move(atInfo), cset, false, cAIC_AllMask,
                          true);
    }
    if (isNew) I->NAtom = nAtom;
    VLACheck(I->CSet, CoordSet*, frame);
    if (I->NCSet <= frame) I->NCSet = frame + 1;
    delete I->CSet[frame];
    I->CSet[frame] = cset;
    if (fractional && cset->Symmetry) {
      CoordSetFracToReal(cset, &cset->Symmetry->Crystal);
    }
    if (ok && isNew)
      ok &= ObjectMoleculeConnect(I, cset, auto_bond, connect_mode);
    if (cset->Symmetry && (!I->Symmetry)) {
      I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
    }
    SceneCountFrames(G);
    if (ok) ok &= ObjectMoleculeExtendIndices(I, frame);
    if (ok) ok &= ObjectMoleculeSort(I);
    if (ok) {
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
    }
  }
  return (I);
}
/**
 * Update coordinates of an exisiting coordset or insert/append a new
 * coordset.
 *
 * coords:      flat coordinate array of length NIndex * 3
 * coords_len:  must be NIndex * 3
 * frame:       coordset index
 */
ObjectMolecule* ObjectMoleculeLoadCoords(PyMOLGlobals* G, ObjectMolecule* I,
                                         const float* coords, int coords_len,
                                         int frame) {
  CoordSet* cset = nullptr;
  int a;
  bool is_new = false;

  if (frame < 0) {
    frame = I->NCSet;
  } else if (frame < I->NCSet) {
    cset = I->CSet[frame];
  }

  if (!cset) {
    // template coordinate set, if available
    cset = I->CSTmpl;

    // find any coordinate set
    for (a = 0; !cset && a < I->NCSet; ++a) cset = I->CSet[a];
    ok_assert(1, cset);
    cset = CoordSetCopy(cset);
    is_new = true;
  }

  // check atom count
  if (coords_len != cset->NIndex * 3) {
    ErrMessage(G, "LoadCoords", "atom count mismatch");
    ok_raise(1);
  }

  // copy coordinates
  for (a = 0; a < coords_len; ++a) {
    cset->Coord[a] = coords[a];
  }

  cset->invalidateRep(cRepAll, cRepInvRep);

  // include coordinate set
  if (is_new) {
    VLACheck(I->CSet, CoordSet*, frame);
    if (I->NCSet <= frame) I->NCSet = frame + 1;
    I->CSet[frame] = cset;
    SceneCountFrames(G);
  }

  // success
  return (I);

  // error handling
ok_except1:
  if (is_new && cset) delete cset;
  ErrMessage(G, "LoadCoords", "failed");
  return nullptr;
}
/* see above... but look up object by name
 */
ObjectMolecule* ObjectMoleculeLoadCoords(PyMOLGlobals* G, const char* name,
                                         const float* coords, int coords_len,
                                         int frame) {
  pymol::CObject* cobj = ExecutiveFindObjectByName(G, name);

  if (!cobj || cobj->type != cObjectMolecule) {
    ErrMessage(G, "LoadCoords", "named object molecule not found.");
    return nullptr;
  }

  return ObjectMoleculeLoadCoords(G, (ObjectMolecule*)cobj, coords, coords_len,
                                  frame);
}
/* see above... but take coordinates from Python list
 *
 * coords: 2d Python float sequence with shape (NIndex, 3)
 */
ObjectMolecule* ObjectMoleculeLoadCoords(PyMOLGlobals* G, ObjectMolecule* I,
                                         PyObject* coords, int frame) {
  CoordSet* cset = nullptr;
  int a, b, l;
  PyObject *v, *w;
  float* f;
  bool is_new = false;

  if (!PySequence_Check(coords)) {
    ErrMessage(G, "LoadCoords", "passed argument is not a sequence");
    ok_raise(1);
  }

  if (frame < 0) {
    frame = I->NCSet;
  } else if (frame < I->NCSet) {
    cset = I->CSet[frame];
  }

  if (!cset) {
    // template coordinate set, if available
    cset = I->CSTmpl;

    // find any coordinate set
    for (a = 0; !cset && a < I->NCSet; ++a) cset = I->CSet[a];
    ok_assert(1, cset);
    cset = CoordSetCopy(cset);
    is_new = true;
  }

  // check atom count
  l = PySequence_Size(coords);
  if (l != cset->NIndex) {
    ErrMessage(G, "LoadCoords", "atom count mismatch");
    ok_raise(1);
  }

  // copy coordinates
  f = cset->Coord.data();
  for (a = 0; a < l; a++) {
    v = PySequence_ITEM(coords, a);

    for (b = 0; b < 3; b++) {
      if (!(w = PySequence_GetItem(v, b))) break;

      f[a * 3 + b] = (float)PyFloat_AsDouble(w);
      Py_DECREF(w);
    }

    Py_DECREF(v);
    ok_assert(2, !PyErr_Occurred());
  }

  cset->invalidateRep(cRepAll, cRepInvRep);

  // include coordinate set
  if (is_new) {
    VLACheck(I->CSet, CoordSet*, frame);
    if (I->NCSet <= frame) I->NCSet = frame + 1;
    I->CSet[frame] = cset;
    SceneCountFrames(G);
  }

  // success
  return (I);

  // error handling
ok_except2:
  PyErr_Print();
ok_except1:
  if (is_new && cset) delete cset;
  ErrMessage(G, "LoadCoords", "failed");
  return nullptr;
}
int ObjectMoleculeExtendIndices(ObjectMolecule* I, int state) {
  int a;
  CoordSet* cs;

  if (I->DiscreteFlag && (state >= 0)) {
    /* if object is discrete, then we don't need to extend each state,
       just the current one (which updates object DiscreteAtmToIdx) */
    cs = I->CSTmpl;
    ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    if (state < I->NCSet) {
      cs = I->CSet[state];
      ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    }
  } else { /* do all states */
    for (a = -1; a < I->NCSet; a++) {
      cs = (a < 0) ? I->CSTmpl : I->CSet[a];
      ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    }
  }
  return true;
ok_except1:
  return false;
}
static char* sdl_mol_prop_name(const char* prop, short* isprop) {
  int i, j, len;
  char* name = nullptr;

  /* for now, need to handle all special characters */
  if (strchr(prop, '/') || strchr(prop, '\\')) {
    if (isprop) *isprop = false;
    name = (char*)strdup(prop);
    return name;
  }
  /* Rule 1: If <prop>, return substring between angle brackets */
  for (i = 0; i < strlen(prop); i++) {
    if (prop[i] == '<') {
      for (len = 0, j = i + 1; j < strlen(prop); j++) {
        if (prop[j] == '>') {
          len = j - i - 1; /* change to j-i+1  to keep <> */
          name = (char*)malloc(len + 1);
          strncpy(name, prop + i + 1, len); /* change to prop+i to keep <> */
          name[len] = '\0';
          if (isprop) *isprop = true;
          return name;
        }
      }
    }
  }

  /* Rule 2: Else if prop contains DT\d+, return that. */
  for (i = 0; i < strlen(prop); i++) {
    if (prop[i] == 'D' && prop[i + 1] == 'T' && isdigit(prop[i + 2])) {
      for (len = 0, j = i + 3; j < strlen(prop); j++) {
        if (!isdigit(prop[j])) {
          len = j - i;
          name = (char*)malloc(len + 1);
          strncpy(name, prop + i, len);
          name[len] = '\0';
          if (isprop) *isprop = true;
          return name;
        }
      }
    }
  }

  /* Rule 3: Else, strip [ribs]_*_ prefix itself. */
  for (i = 0; i < strlen(prop); i++) {
    if (prop[i] == '_') {
      for (j = i + 1; j < strlen(prop); j++) {
        if (prop[j] == '_') {
          name = (char*)malloc(strlen(prop + j + 1) + 1);
          strcpy(name, prop + j + 1);
          if (isprop) *isprop = true;
          return name;
        }
      }
    }
  }

  /* Rule 4: else return prop itself */
  if (isprop) *isprop = false;
  name = (char*)strdup(prop);
  return name;
}

static short is_sdl_mol_prop_name(const char* prop) {
  short isprop;
  char* cctmp = pymol::malloc<char>(strlen(prop) + 1);
  strcpy(cctmp, prop);
  free(sdl_mol_prop_name(cctmp, &isprop));
  FreeP(cctmp);
  return isprop;
}
static CoordSet* ObjectMoleculeMOLStr2CoordSet(
    PyMOLGlobals* G, const char* buffer, AtomInfoType** atInfoPtr,
    const char** restart, short loadpropertiesall, OVLexicon* loadproplex) {
  const char* p;
  int nAtom, nBond;
  int a, cnt, atm, chg;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType* atInfo = nullptr;
  char cc[MAXLINELEN], cc1[MAXLINELEN], resn[MAXLINELEN] = "UNK";
  float* f;
  BondType* ii;
  BondType* bond = nullptr;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  WordType nameTmp;

  p = buffer;
  nAtom = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;

  /*  p=ParseWordCopy(nameTmp,p,sizeof(WordType)-1); */
  p = ParseNCopy(nameTmp, p, sizeof(WordType) - 1);
  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
  " ObjMolMOLStr2CoordSet: title '%s'\n", nameTmp ENDFB(G) p = nextline(p);
  p = nextline(p);
  p = nextline(p);

  if (ok) {
    p = ncopy(cc, p, 3);
    if (sscanf(cc, "%d", &nAtom) != 1)
      ok = ErrMessage(G, "ReadMOLFile", "bad atom count");
  }

  if (ok) {
    p = ncopy(cc, p, 3);
    if (sscanf(cc, "%d", &nBond) != 1)
      ok = ErrMessage(G, "ReadMOLFile", "bad bond count");
  }

  if (ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if (atInfo) VLACheck(atInfo, AtomInfoType, nAtom);
  }

  p = nextline(p);

  /* read coordinates and atom names */

  if (ok) {
    f = coord;
    for (a = 0; a < nAtom; a++) {
      if (ok) {
        p = ncopy(cc, p, 10);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if (ok) {
        p = ncopy(cc, p, 10);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if (ok) {
        p = ncopy(cc, p, 10);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if (ok) {
        p = nskip(p, 1);
        p = ntrim(cc, p, 3);
        strncpy(atInfo[a].elem, cc, cElemNameLen);
        atInfo[a].name = LexIdx(G, cc);
        atInfo[a].visRep = auto_show;
      }
      if (ok) {
        int tmp_int;
        p = nskip(p, 2);
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%hhi", &atInfo[a].formalCharge) == 1) {
          if (atInfo[a].formalCharge) {
            atInfo[a].formalCharge = 4 - atInfo[a].formalCharge;
          }
        }
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%d", &tmp_int) != 1)
          atInfo[a].stereo = 0;
        else
          atInfo[a].stereo = tmp_int;
      }
      if (ok && atInfo) {
        atInfo[a].id = a + 1;
        atInfo[a].rank = a;
        atInfo[a].resn = LexIdx(G, resn);
        atInfo[a].hetatm = true;
        AtomInfoAssignParameters(G, atInfo + a);
        AtomInfoAssignColors(G, atInfo + a);
        atInfo[a].alt[0] = 0;
        atInfo[a].segi = 0;
        atInfo[a].inscode = 0;
      }
      p = nextline(p);
      if (!ok) break;
    }
  }
  if (ok) {
    bond = VLACalloc(BondType, nBond);
    ii = bond;
    for (a = 0; a < nBond; a++) {
      if (ok) {
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%d", &ii->index[0]) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond atom");
      }

      if (ok) {
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%d", &ii->index[1]) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond atom");
      }

      if (ok) {
        int order = 0;
        p = ncopy(cc, p, 3);
        if (sscanf(cc, "%d", &order) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond order");
        ii->order = order;
      }
      if (ok) {
        p = ncopy(cc, p, 3);
        // stereo
      }
      ii++;
      if (!ok) break;
      p = nextline(p);
    }
    ii = bond;
    for (a = 0; a < nBond; a++) {
      ii->index[0]--; /* adjust bond indexs down one */
      ii->index[1]--;
      ii++;
    }
  }
  while (*p) { /* read M  CHG records */
    auto p_line = p;
    p = ncopy(cc, p, 6);
    if (cc[5] == ' ') cc[5] = 0;
    if ((!strcmp(cc, "M  END")) || (!strcmp(cc, "M END"))) {
      /* denotes end of MOL block */
      break;
    }
    if ((!strcmp(cc, "M  CHG")) || (!strcmp(cc, "M CHG"))) {
      p = ncopy(cc, p, 3);
      if (sscanf(cc, "%d", &cnt) == 1) {
        while (cnt--) {
          p = ncopy(cc, p, 4);
          p = ncopy(cc1, p, 4);
          if (!((*cc) || (*cc1))) break;
          if ((sscanf(cc, "%d", &atm) == 1) && (sscanf(cc1, "%d", &chg) == 1)) {
            atm--;
            if ((atm >= 0) && (atm < nAtom)) atInfo[atm].formalCharge = chg;
          }
        }
      }
    } else if (!strcmp(cc, "M  V30")) {
      p = MOLV3000Parse(G, p_line, atInfo, bond, coord, nAtom, nBond);
      if (!p) {
        ok = false;
        break;
      }
      continue;
    }
    p = nextline(p);
  }
  if (ok) {
    (*restart) = p;
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, nameTmp);
    p = nextline(p);
    // Need to read in the properties
    {
      OVreturn_word ovresult;
      const char* ptmp = p;
      const char* np = ptmp;
      char *cc = VLACalloc(char, 1024), *cctmp = VLACalloc(char, 1024);
      int ccp = 0;
      char* prop_name = nullptr;
      short is_prop_name, line = 0;
      int nchs;
      cc[0] = 0;
      ncopy(cctmp, np, 4);
      cctmp[4] = 0;
      while (strcmp(cctmp, "$$$$")) {
        ptmp = np;
        np = nextline(np);
        if (np == ptmp) break;
        nchs = (np - ptmp - 1);
        if (nchs) {
          VLACheck(cctmp, char, nchs + 1);
          ncopy(cctmp, ptmp, nchs);
          cctmp[nchs] = 0;
          is_prop_name = is_sdl_mol_prop_name(cctmp);
          if (is_prop_name) {
            // is_prop_name
            if (prop_name) {
              short load = loadpropertiesall;
              if (!load && loadproplex) {
                ovresult = OVLexicon_BorrowFromCString(loadproplex, prop_name);
                load = (ovresult.status == OVstatus_SUCCESS);
              }
              if (load) {
                PropertyCheckUniqueID(G, cset);
                PropertySetromString(G, cset->prop_id, prop_name, cc);
                ccp = 0;
                cc[0] = 0;
                line = 0;
              }
              free(prop_name);
            }
            prop_name = sdl_mol_prop_name(cctmp, nullptr);
          } else {
            if (prop_name) {
              int sz = strlen(cctmp);
              if (line == 1) {
                VLACheck(cc, char, ccp + 2);
                cc[ccp++] = '\n';
              }
              VLACheck(cc, char, ccp + sz + 2);
              ncopy(&cc[ccp], cctmp, sz);
              ccp += sz;
              VLACheck(cc, char, ccp + 2);
              if (line) {
                cc[ccp++] = '\n';
              }
              cc[ccp] = 0;
              line++;
            }
          }
        }
        ncopy(cctmp, np, 4);
      }
      if (prop_name) {
        short load = loadpropertiesall;
        if (!load && loadproplex) {
          ovresult = OVLexicon_BorrowFromCString(loadproplex, prop_name);
          load = (ovresult.status == OVstatus_SUCCESS);
        }
        if (load) {
          PropertyCheckUniqueID(G, cset);
          PropertySetromString(G, cset->prop_id, prop_name, cc);
          ccp = 0;
          cc[0] = 0;
          line = 0;
        }
        free(prop_name);
      }
      VLAFreeP(cc);
      VLAFreeP(cctmp);
    }
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
    (*restart) = nullptr;
  }
  if (atInfoPtr) *atInfoPtr = atInfo;
  return (cset);
}
static CoordSet* ObjectMoleculeSDF2Str2CoordSet(
    PyMOLGlobals* G, const char* buffer, AtomInfoType** atInfoPtr,
    const char** next_mol, short loadpropertiesall, OVLexicon* loadproplex) {
  char cc[MAXLINELEN];
  const char* p;
  CoordSet* result = nullptr;
  result = ObjectMoleculeMOLStr2CoordSet(G, buffer, atInfoPtr, next_mol,
                                         loadpropertiesall, loadproplex);
  p = *next_mol;
  if (p) {
    while (*p) { /* we simply need to skip until we've read past the end of the
                    SDF record */
      p = ncopy(cc, p, 4);
      p = nextline(p);
      if (!strcmp(cc, "$$$$")) { /* SDF record separator... */
        break;
      }
    }
    if (!*p) p = nullptr;
  }
  *next_mol = p;
  return result;
}

/**
 * Get the sum of bond orders for every atom. For this purpose, aromatic
 * bonds will be considered either single or double bond, based on
 * the entire system (heuristic method, room for improvement).
 */
static std::vector<signed char> get_bond_order_sums(ObjectMolecule* obj) {
  std::vector<signed char> valences(obj->NAtom);
  std::vector<signed char> freevalences(obj->NAtom);
  std::vector<signed char> orders(obj->NBond);

  // bond order sums as if all aromatic bonds were single bonds
  for (int b = 0; b < obj->NBond; ++b) {
    auto bond = obj->Bond + b;
    auto order = bond->order;
    orders[b] = order;
    if (order > 3) order = 1;
    valences[bond->index[0]] += order;
    valences[bond->index[1]] += order;
  }

  // determine free valences as if all aromatic bonds were single bonds
  for (int atm = 0; atm < obj->NAtom; ++atm) {
    int tmp = 0;
    switch (obj->AtomInfo[atm].protons) {
      case cAN_C:
        tmp = 4;
        break;
      case cAN_N:
      case cAN_P:
        tmp = 5;
        break;
      case cAN_O:
      case cAN_S:
        tmp = 6;
        break;
    }
    if (tmp) {
      tmp -= valences[atm];
      if (tmp) {
        freevalences[atm] = (tmp - 1) % 2 + 1;
      }
    }
  }

  // (This is a heuristic, gets most of the ZINC dataset correct)
  // Do two passes over all aromatic bonds.
  // 1st pass: assign double bonds between atoms with `freevalence == 1`
  // 2nd pass: assign double bonds between atoms with `freevalence >= 1`
  for (int secondpass = 0; secondpass < 2; ++secondpass) {
    for (int b = 0; b < obj->NBond; ++b) {
      if (orders[b] == 4) {
        auto atm1 = obj->Bond[b].index[0];
        auto atm2 = obj->Bond[b].index[1];
        if (secondpass ? (freevalences[atm1] >= 1 && freevalences[atm2] >= 1)
                       : (freevalences[atm1] == 1 && freevalences[atm2] == 1)) {
          freevalences[atm1] = 0;
          freevalences[atm2] = 0;
          valences[atm1] += 1;
          valences[atm2] += 1;
          orders[b] = 2;
        }
      }
    }
  }

  return valences;
}

/**
 * For each atom, set `formalCharge` based on mol2 `textType`
 *
 * See also: getMOL2Type() in src/pymol/chem/mol2_typing.cc
 */
static void ObjectMoleculeMOL2SetFormalCharges(PyMOLGlobals* G,
                                               ObjectMolecule* obj) {
  // check if structure has explicit hydrogens
  bool has_hydrogens = false;
  for (int at = 0; at < obj->NAtom; ++at) {
    auto ai = obj->AtomInfo + at;
    if (ai->isHydrogen()) {
      has_hydrogens = true;
      break;
    }
  }

  if (!has_hydrogens) {
    // could eventually do PDB nomenclature charge assignment here...
    return;
  }

  // Unfortunately, aromatic bonds (type 4) can't be used to determine
  // formal charges. The following uses a heuristic to guess bond orders
  // for aromatic bonds.
  auto valences = get_bond_order_sums(obj);

  // (period currently incompatible with G->lex_const)
  auto lex_N_4 = LexBorrow(G, "N.4");

  for (int at = 0; at < obj->NAtom; ++at) {
    int fcharge = 0;
    auto ai = obj->AtomInfo + at;

    if (ai->protons == cAN_N) {
      if (ai->textType == lex_N_4) {
        fcharge = 1;
      } else if (valences[at] == 2) {
        fcharge = -1;
      } else if (valences[at] == 4) {
        fcharge = 1;
      }
    } else if (ai->protons == cAN_O) {
      if (valences[at] == 1) {
        fcharge = -1;
      }
    }

    ai->formalCharge = fcharge;
  }
}
static CoordSet* ObjectMoleculeMOL2Str2CoordSet(PyMOLGlobals* G,
                                                const char* buffer,
                                                AtomInfoType** atInfoPtr,
                                                const char** next_mol) {
  const char* p;
  int nAtom, nBond, nSubst, nFeat, nSet;
  int a;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType* atInfo = nullptr;
  char cc[MAXLINELEN];
  float* f;
  const char* last_p;
  BondType* ii;
  BondType* bond = nullptr;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  int have_molecule = false;
  WordType nameTmp;
  bool has_non_hetatm = false;

  // Maestro writes <0> as subst_name for waters
  const lexidx_t empty_subst_name = LexIdx(G, "<0>");

  p = buffer;
  nAtom = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;

  while ((*p) && ok) {
    last_p = p;
    /* top level -- looking for an RTI, assuming p points to the beginning of a
     * line */
    p = ParseWordCopy(cc, p, MAXLINELEN);
    if (cc[0] == '@') {
      if (WordMatchExact(G, cc, "@<TRIPOS>MOLECULE", true) ||
          WordMatchExact(G, cc, "@MOLECULE", true)) {
        if (have_molecule) {
          *next_mol = last_p;
          break; /* next record of multi-mol2 */
        }
        p = ParseNextLine(p);
        p = ParseNTrim(nameTmp, p, sizeof(WordType) - 1); /* get mol name */
        p = ParseNextLine(p);
        p = ParseWordCopy(cc, p, MAXLINELEN);
        if (sscanf(cc, "%d", &nAtom) != 1)
          ok = ErrMessage(G, "ReadMOL2File", "bad atom count");
        else {
          coord = VLAlloc(float, 3 * nAtom);
          if (atInfo) VLACheck(atInfo, AtomInfoType, nAtom);

          p = ParseWordCopy(cc, p, MAXLINELEN);
          if (sscanf(cc, "%d", &nBond) != 1) {
            nBond = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if (sscanf(cc, "%d", &nSubst) != 1) {
            nSubst = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if (sscanf(cc, "%d", &nFeat) != 1) {
            nFeat = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if (sscanf(cc, "%d", &nSet) != 1) {
            nSet = 0;
          }
        }
        p = ParseNextLine(p);
        p = ParseNextLine(p);
        have_molecule = true;
      } else if (WordMatchExact(G, cc, "@<TRIPOS>ATOM", true) ||
                 WordMatchExact(G, cc, "@ATOM", true)) {
        if (!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@ATOM before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        f = coord;
        for (a = 0; a < nAtom; a++) {
          AtomInfoType* ai = atInfo + a;

          if (ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if (sscanf(cc, "%d", &ai->id) != 1) {
              ok = ErrMessage(G, "ReadMOL2File", "bad atom id");
            }
          }
          if (ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            cc[cAtomNameLen] = 0;
            UtilCleanStr(cc);
            if (cc[0])
              LexAssign(G, ai->name, cc);
            else
              ok = ErrMessage(G, "ReadMOL2File", "bad atom name");
          }
          if (ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if (sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad x coordinate");
          }
          if (ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if (sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad y coordinate");
          }
          if (ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if (sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad z coordinate");
          }
          if (ok) {
            OrthoLineType temp;
            p = ParseWordCopy(cc, p, OrthoLineLength - 1);
            if (sscanf(cc, "%s", temp) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad atom type");
            else { /* convert atom type to elem symbol */
              char* tt = temp;
              char* el = ai->elem;
              int elc = 0;
              while (*tt && ((*tt) != '.')) {
                *(el++) = *(tt++);
                elc++;
                if (elc > cElemNameLen) break;
              }
              *el = 0;
              if (el[2]) el[0] = 0;

              ai->textType = LexIdx(G, temp);
            }
          }
          if (ok) {
            char cc_resi[8];
            p = ParseWordCopy(cc_resi, p, sizeof(cc_resi) - 1);
            if (cc_resi[0]) { /* subst_id is residue identifier */
              ai->setResi(cc_resi);
              p = ParseWordCopy(cc, p, MAXLINELEN);
              if (cc[0]) {
                if (nSubst == 0) {
                  // without substructure information (e.g. OpenBabel exported
                  // MOL2): if subst_name includes the subst_id (e.g. 5 ALA5)
                  // then strip the number
                  int len_resi = strlen(cc_resi);
                  int len_resn = strlen(cc);
                  if (len_resn > len_resi) {
                    if (strcmp(cc_resi, cc + len_resn - len_resi) == 0) {
                      cc[len_resn - len_resi] = 0;
                    }
                  }
                }

                LexAssign(G, ai->resn, cc);

                p = ParseWordCopy(cc, p, MAXLINELEN);
                if (cc[0]) {
                  if (sscanf(cc, "%f", &ai->partialCharge) != 1)
                    ok = ErrMessage(G, "ReadMOL2File", "bad atom charge");
                }

                // status_bit
                p = ParseWordCopy(cc, p, MAXLINELEN);
                if (cc[0]) {
                  // for water, substitute resn "<0>" with "HOH" for
                  // correct classification as "solvent"
                  if (ai->resn == empty_subst_name && strstr(cc, "WATER")) {
                    LexAssign(G, ai->resn, G->lex_const.HOH);
                  }
                }
              }
            }
          }
          p = ParseNextLine(p);

          ai->visRep = auto_show;

          ai->id = a + 1;
          ai->rank = a;
          ai->hetatm = true;
          AtomInfoAssignParameters(G, atInfo + a);
          AtomInfoAssignColors(G, atInfo + a);
          ai->alt[0] = 0;
          ai->segi = 0;
        }
      } else if (WordMatchExact(G, cc, "@<TRIPOS>SET", true) ||
                 WordMatchExact(G, cc, "@SET", true)) {
        if (!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@SET before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        if (ok) {
          for (a = 0; a < nSet; a++) {
            if (ok) {
              char ss = 0;
              p = ParseWordCopy(cc, p, 50);
              cc[5] = 0;
              if (!strcmp("HELIX", cc)) ss = 'H';
              if (!strcmp("SHEET", cc)) ss = 'S';
              p = ParseWordCopy(cc, p, 50);
              if (!strcmp("STATIC", cc)) {
                int nMember;

                p = ParseNextLine(p);
                p = ParseWordCopy(cc, p, 20);
                if (sscanf(cc, "%d", &nMember) != 1)
                  ok = ErrMessage(G, "ReadMOL2File", "bad member count");
                else {
                  while (ok && (nMember > 0)) {
                    p = ParseWordCopy(cc, p, 20);
                    if ((!cc[0]) && (!*p)) {
                      ok = false;
                      break;
                    }
                    if (cc[0] != '\\') {
                      int atom_id;
                      if (sscanf(cc, "%d", &atom_id) != 1) {
                        ok = ErrMessage(G, "ReadMOL2File", "bad member");
                      } else {
                        atom_id--;
                        if ((atom_id >= 0) && (atom_id < nAtom)) {
                          atInfo[atom_id].ssType[0] = ss;
                          atInfo[atom_id].ssType[1] = 0;
                        }
                        nMember--;
                      }
                    } else {
                      p = ParseNextLine(p);
                    }
                  }
                }
                p = ParseNextLine(p);
              } else {
                p = ParseNextLine(p);
                p = ParseNextLine(p);
              }
            }
          }
        }
      } else if (WordMatchExact(G, cc, "@<TRIPOS>BOND", true) ||
                 WordMatchExact(G, cc, "@BOND", true)) {
        if (!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@BOND before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);

        if (ok) {
          bond = VLACalloc(BondType, nBond);
          ii = bond;
          for (a = 0; a < nBond; a++) {
            if (ok) {
              p = ParseWordCopy(cc, p, 20);
            }

            if (ok) {
              p = ParseWordCopy(cc, p, 20);
              if (sscanf(cc, "%d", &ii->index[0]) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad source atom id");
            }

            if (ok) {
              p = ParseWordCopy(cc, p, 20);
              if (sscanf(cc, "%d", &ii->index[1]) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad target atom id");
            }

            if (ok) {
              p = ParseWordCopy(cc, p, 20);
              if (!cc[1]) {
                switch (cc[0]) {
                  case '1':
                    ii->order = 1;
                    break;
                  case '2':
                    ii->order = 2;
                    break;
                  case '3':
                    ii->order = 3;
                    break;
                  case '4':
                    ii->order = 4;
                    break;
                }
              } else if (WordMatchExact(G, "ar", cc, true)) {
                ii->order = 4;
              } else if (WordMatchExact(G, "am", cc, true)) {
                ii->order = 1;
              } else if (WordMatchExact(G, "un", cc, true)) {
                ii->order = 1;
              } else if (WordMatchExact(G, "nc", cc, true)) {
                ii->order = 0; /* is this legal in PyMOL? */
              } else if (WordMatchExact(G, "du", cc, true)) {
                ii->order = 1;
              } else
                ok = ErrMessage(G, "ReadMOL2File", "bad bond type");
            }
            ii++;
            if (!ok) break;
            p = ParseNextLine(p);
          }
          ii = bond;
          for (a = 0; a < nBond; a++) {
            ii->index[0]--; /* adjust bond indexs down one */
            ii->index[1]--;
            ii++;
          }
        }
      } else if (WordMatchExact(G, cc, "@<TRIPOS>SUBSTRUCTURE", true) ||
                 WordMatchExact(G, cc, "@SUBSTRUCTURE", true)) {
        if (!have_molecule) {
          ok =
              ErrMessage(G, "ReadMOL2File", "@SUBSTSRUCTURE before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        if (ok) {
          WordType subst_name;
          SegIdent segment; /* what MOL2 calls chain */
          lexidx_t chain = 0;
          WordType subst_type;
          ResIdent resi;
          ResName resn;
          int id, dict_type, root, resv;
          int end_line, seg_flag, subst_flag, resi_flag;
          int chain_flag, resn_flag;
          std::unordered_map<int, std::vector<int>> resv_map;

          {
            int b;
            AtomInfoType* ai = atInfo;
            for (b = 0; b < nAtom; b++) {
              resv_map[ai->resv].push_back(b);
              ai++;
            }
          }

          for (a = 0; a < nSubst; a++) {
            bool hetatm = true;
            segment[0] = 0;
            subst_name[0] = 0;
            LexAssign(G, chain, 0);
            resn[0] = 0;
            resi[0] = 0;
            end_line = false;
            seg_flag = false;
            subst_flag = false;

            resi_flag = false;
            chain_flag = false;
            resn_flag = false;

            if (ok) {
              const char* save_p = p;
              p = ParseWordCopy(cc, p, 20);
              if (sscanf(cc, "%d", &id) != 1) {
                if (cc[0] != '@')
                  ok = ErrMessage(G, "ReadMOL2File", "bad substructure id");
                else {
                  p = save_p;
                  break; /* missing substructure... */
                }
              }
            }
            if (ok) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              if (sscanf(cc, "%s", subst_name) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad substructure name");
            }
            if (ok) {
              p = ParseWordCopy(cc, p, 20);
              if (sscanf(cc, "%d", &root) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad target root atom id");
              else
                root--; /* convert 1-based to 0-based */
            }

            /* optional data */
            if (ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              if (sscanf(cc, "%s", subst_type) != 1) {
                end_line = true;
                subst_type[0] = 0;
              } else {
                subst_flag = 1;
              }
            }

            if (ok && (!end_line)) {
              p = ParseWordCopy(cc, p, 20);
              if (sscanf(cc, "%d", &dict_type) != 1) end_line = true;
            }

            if (ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              cc[cSegiLen] = 0;
              if (sscanf(cc, "%s", segment) != 1) {
                end_line = true;
                segment[0] = 0;
              } else if (strcmp(segment, "****")) {
                seg_flag = true;
                if (!segment[1]) { /* if segment is single letter, then also
                                      assign to chain field */
                  LexAssign(G, chain, segment);
                  chain_flag = true;
                }
              }
            }

            if (ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              cc[cResnLen] = 0;
              if (!strcmp(cc, "****")) { /* whoops -- no residue name... */
                char* pp;
                UtilNCopy(resn, subst_name, sizeof(ResName));
                pp = resn;
                /* any numbers? */
                while (*pp) {
                  if (((*pp) >= '0') && ((*pp) <= '9')) break;
                  pp++;
                }
                /* trim them off */
                while (pp >= resn) {
                  if (((*pp) >= '0') && ((*pp) <= '9'))
                    *pp = 0;
                  else
                    break;
                  pp--;
                }

                if (resn[0]) resn_flag = true;
              } else {
                if (sscanf(cc, "%s", resn) != 1) {
                  end_line = true;
                  resn[0] = 0;
                } else {
                  resn_flag = true;
                }
              }
            }

            if (ok && (root >= 0) && (root < nAtom)) {
              if ((strlen(subst_name) == 4) &&
                  ((subst_name[0] >= '1') && (subst_name[0] <= '9')) &&
                  ((((subst_name[1] >= 'A') && (subst_name[1] <= 'Z')) ||
                    ((subst_name[1] >= 'a') && (subst_name[1] <= 'z'))) ||
                   (((subst_name[2] >= 'A') && (subst_name[2] <= 'Z')) ||
                    ((subst_name[2] >= 'a') && (subst_name[2] <= 'z'))) ||
                   (((subst_name[3] >= 'A') && (subst_name[3] <= 'Z')) ||
                    ((subst_name[3] >= 'a') && (subst_name[3] <= 'z'))))) {
                /* if subst_name looks like a PDB code, then it isn't a residue
                   identifier so do nothing... */
              } else {
                if (!subst_flag) {
                  /* Merck stuffs the chain ID and the residue ID in the
                     subst_name field, so handle that case first */

                  /* get the residue value */

                  ParseIntCopy(cc, subst_name, 20);
                  if (sscanf(cc, "%d", &resv) == 1) {
                    /* we have the resv, so now get the chain if there is one */
                    char* pp = subst_name;
                    if (!((pp[0] >= '0') && (pp[0] <= '9'))) {
                      char tmp[2] = {pp[0], 0};
                      LexAssign(G, chain, tmp);
                      chain_flag = true;
                      pp++;
                    }
                    UtilNCopy(
                        resi, pp,
                        cResiLen); /* now get the textual residue identifier */

                    if (resi[0]) {
                      resi_flag = true;
                    }
                  }
                } else {
                  /* normal mode: assume that the substructure ID is a vanilla
                   * residue identifier */

                  char* pp = subst_name;

                  if (resn_flag) { /* if we have a resn, then remove it from the
                                      subst_name */
                    char* qq = resn;
                    while ((*pp) && (*qq)) {
                      if (*pp == *qq) {
                        pp++;
                        qq++;
                      } else
                        break;
                    }
                  }

                  ParseIntCopy(cc, pp, 20);
                  if (sscanf(cc, "%d", &resv) == 1) {
                    UtilNCopy(
                        resi, pp,
                        cResiLen); /* now get the textual residue identifier */
                    if (resi[0]) {
                      resi_flag = true;
                    }
                  }

                  if (strcmp(subst_type, "RESIDUE") == 0) {
                    hetatm = false;
                    has_non_hetatm = true;
                  }
                }
              }

              /* for now, assume that atom ids are 1-based and sequential */

              if (ok) {
                if (resi_flag || chain_flag || resn_flag || seg_flag) {
                  auto it = resv_map.find(id);
                  if (it != resv_map.end()) {
                    for (auto b : it->second) {
                      assert(b >= 0 && b < nAtom);
                      {
                        AtomInfoType* ai = atInfo + b;
                        if (resi_flag) ai->setResi(resi);
                        if (chain_flag) {
                          LexAssign(G, ai->chain, chain);
                        }
                        if (resn_flag) LexAssign(G, ai->resn, resn);
                        if (seg_flag) LexAssign(G, ai->segi, segment);
                        ai->hetatm = hetatm;
                      }
                    }
                  }
                }
              }
            }
            if (!ok) break;
            p = ParseNextLine(p);
          }
          LexDec(G, chain);
        }
      } else
        p = ParseNextLine(p);
    } else
      p = ParseNextLine(p);
  }

  LexDec(G, empty_subst_name);

  if (has_non_hetatm) {
    // contains "residues" (assume polymer) so ignore hetatm for surfacing
    for (auto ai = atInfo, ai_end = atInfo + nAtom; ai != ai_end; ++ai) {
      if (ai->hetatm) {
        ai->flags |= cAtomFlag_ignore;
      }
    }
  }

  if (ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, nameTmp);
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if (atInfoPtr) *atInfoPtr = atInfo;
  return (cset);
}

/**
 * Read one molecule in MOL2, MOL, SDF, MMD or XYZ format from the string
 * pointed to by (*next_entry). All these formats (except MOL) support multiple
 * concatenated entries in one file. If multiplex=1, then read only one
 * molecule (one state) and set the next_entry pointer to the beginning of the
 * next entry. Otherwise, read a multi-state molecule.
 */
ObjectMolecule* ObjectMoleculeReadStr(PyMOLGlobals* G, ObjectMolecule* I,
                                      const char** next_entry,
                                      cLoadType_t content_format, int frame,
                                      int discrete, int quiet, int multiplex,
                                      char* new_name, short loadpropertiesall,
                                      OVLexicon* loadproplex) {
  int ok = true;
  CoordSet* cset = nullptr;
  pymol::vla<AtomInfoType> atInfo;
  int isNew;
  int nAtom;
  const char *restart = nullptr, *start = *next_entry;
  int repeatFlag = true;
  int successCnt = 0;
  char tmpName[WordLength];
  int deferred_tasks = false;
  int skip_out;
  int connect = false;
  int set_formal_charges = false;
  *next_entry = nullptr;
  int aic_mask = cAIC_MOLMask;

  while (repeatFlag) {
    repeatFlag = false;
    skip_out = false;

    if (!I)
      isNew = true;
    else
      isNew = false;

    if (isNew) {
      I = (ObjectMolecule*)new ObjectMolecule(G, (discrete > 0));
      std::swap(atInfo, I->AtomInfo);
    } else {
      atInfo = pymol::vla<AtomInfoType>(10);
    }

    if (isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }

    restart = nullptr;
    switch (content_format) {
      case cLoadTypeMOL2:
      case cLoadTypeMOL2Str:
        cset = ObjectMoleculeMOL2Str2CoordSet(G, start, &atInfo, &restart);
        if (cset) {
          set_formal_charges = true;
        }
        break;
      case cLoadTypeMOL:
      case cLoadTypeMOLStr:
        cset = ObjectMoleculeMOLStr2CoordSet(G, start, &atInfo, &restart,
                                             loadpropertiesall, loadproplex);
        restart = nullptr;
        break;
      case cLoadTypeSDF2:
      case cLoadTypeSDF2Str:
        cset = ObjectMoleculeSDF2Str2CoordSet(G, start, &atInfo, &restart,
                                              loadpropertiesall, loadproplex);
        break;
      case cLoadTypeXYZ:
      case cLoadTypeXYZStr:
        cset = ObjectMoleculeXYZStr2CoordSet(G, start, &atInfo, &restart);
        if (!cset->TmpBond) connect = true;
        break;
      case cLoadTypeMMD:
      case cLoadTypeMMDStr:
        cset = ObjectMoleculeMMDStr2CoordSet(G, start, &atInfo, &restart);
        aic_mask = cAIC_MMDMask;
        break;
    }

    if (!cset) {
      if (!isNew) VLAFreeP(atInfo);
      if (!successCnt) {
        if (isNew) std::swap(I->AtomInfo, atInfo);
        DeleteP(I);
        ok = false;
      } else {
        skip_out = true;
      }
    }

    if (ok && !skip_out) {
      if ((discrete > 0 &&
           !I->DiscreteFlag) ||  // should set to discrete if explicitly defined
          ((discrete < 0) && (restart) && isNew && (multiplex <= 0))) {
        /* if default discrete behavior and new object, with
           multi-coordinate set file, and not multiplex, then set
           discrete */

        ObjectMoleculeSetDiscrete(G, I, true);
      }

      if (frame < 0) frame = I->NCSet;
      if (I->NCSet <= frame) I->NCSet = frame + 1;
      VLACheck(I->CSet, CoordSet*, frame);

      nAtom = cset->NIndex;

      if (I->DiscreteFlag && atInfo) {
        int a;
        int fp1 = frame + 1;
        AtomInfoType* ai = atInfo.data();
        for (a = 0; a < nAtom; a++) {
          (ai++)->discrete_state = fp1;
        }
      }

      if (multiplex > 0) UtilNCopy(tmpName, cset->Name, WordLength);

      cset->Obj = I;
      cset->enumIndices();
      cset->invalidateRep(cRepAll, cRepInvRep);
      if (isNew) {
        std::swap(I->AtomInfo, atInfo);
      } else {
        ObjectMoleculeMerge(I, std::move(atInfo), cset, false, aic_mask, false);
        /* NOTE: will release atInfo */
      }

      if (isNew) I->NAtom = nAtom;
      if (frame < 0) frame = I->NCSet;
      VLACheck(I->CSet, CoordSet*, frame);
      if (I->NCSet <= frame) I->NCSet = frame + 1;
      delete I->CSet[frame];
      I->CSet[frame] = cset;

      if (ok && isNew) ok &= ObjectMoleculeConnect(I, cset, connect);

      if (ok) {
        if (I->DiscreteFlag) {
          ok &= ObjectMoleculeExtendIndices(I, frame);
        }
      }

      deferred_tasks = true;
      successCnt++;
      if (!quiet) {
        if (successCnt > 1) {
          if (successCnt == 2) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
            " %s: read through molecule %d.\n", __func__, 1 ENDFB(G);
          }
          if (UtilShouldWePrintQuantity(successCnt)) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
            " %s: read through molecule %d.\n", __func__, successCnt ENDFB(G);
          }
        }
      }
      if (multiplex > 0) {
        UtilNCopy(new_name, tmpName, WordLength);
        if (restart) {
          *next_entry = restart;
        }
      } else if (restart) {
        repeatFlag = true;
        start = restart;
        frame = frame + 1;
      }
    }
  }
  if (deferred_tasks &&
      I) { /* defer time-consuming tasks until all states have been loaded */
    if (!I->DiscreteFlag) {
      ObjectMoleculeExtendIndices(I, cStateAll);
    }
    ObjectMoleculeSort(I);
    if (set_formal_charges) {
      ObjectMoleculeMOL2SetFormalCharges(G, I);
    }
    SceneCountFrames(G);
    I->invalidate(cRepAll, cRepInvAll, -1);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return (I);
}
ObjectMolecule* ObjectMoleculeReadPDBStr(PyMOLGlobals* G, ObjectMolecule* I,
                                         const char* PDBStr, int state,
                                         int discrete, char* pdb_name,
                                         const char** next_pdb,
                                         PDBInfoRec* pdb_info, int quiet,
                                         int* model_number) {
  CoordSet* cset = nullptr;
  pymol::vla<AtomInfoType> atInfo;
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;
  const char *start, *restart = nullptr;
  int repeatFlag = true;
  int successCnt = 0;
  unsigned int aic_mask = cAIC_PDBMask;

  SegIdent segi_override = ""; /* saved segi for corrupted NMR pdb files */

  start = PDBStr;
  while (repeatFlag) {
    repeatFlag = false;

    if (!I)
      isNew = true;
    else
      isNew = false;

    if (ok) {
      if (isNew) {
        I = (ObjectMolecule*)new ObjectMolecule(G, discrete);
        CHECKOK(ok, I);
        if (ok) std::swap(atInfo, I->AtomInfo);
        isNew = true;
      } else {
        atInfo = pymol::vla<AtomInfoType>(10);
        CHECKOK(ok, atInfo);
        isNew = false;
      }
      if (ok && isNew) {
        I->Color = AtomInfoUpdateAutoColor(G);

        if (pdb_info->variant == PDB_VARIANT_VDB ||
            pdb_info->variant == PDB_VARIANT_PQR) {
          // pqr files have no chain identifier by default
          // vdb files have same chain identifiers in all symmetry copies
          SettingSet(cSetting_retain_order, 1, I);
        }
      }
      if (ok)
        cset = ObjectMoleculePDBStr2CoordSet(G, start, &atInfo, &restart,
                                             segi_override, pdb_name, next_pdb,
                                             pdb_info, quiet, model_number);
      CHECKOK(ok, cset);
      if (ok) {
        nAtom = cset->NIndex;
      }
    }
    if (ok && pdb_name && (*next_pdb)) {
      /* problematic scenario? */
    }

    /* include coordinate set */
    if (ok) {
      if (I->DiscreteFlag && atInfo) {
        unsigned int a;
        int fp1 = state + 1;
        AtomInfoType* ai = atInfo.data();
        for (a = 0; a < nAtom; a++) {
          (ai++)->discrete_state = fp1;
        }
      }

      cset->Obj = I;
      cset->enumIndices();
      cset->invalidateRep(cRepAll, cRepInvRep);
      if (isNew) {
        std::swap(I->AtomInfo, atInfo);
      } else {
        ok &= ObjectMoleculeMerge(I, std::move(atInfo), cset, true, aic_mask,
                                  true);
        /* NOTE: will release atInfo */
      }
      if (isNew) I->NAtom = nAtom;
      if (state < 0) state = I->NCSet;
      if (*model_number > 0) {
        if (SettingGetGlobal_b(G, cSetting_pdb_honor_model_number))
          state = *model_number - 1;
      }
      VLACheck(I->CSet, CoordSet*, state);
      CHECKOK(ok, I->CSet);
      if (ok) {
        if (I->NCSet <= state) I->NCSet = state + 1;
        delete I->CSet[state];
        I->CSet[state] = cset;
      }
      if (ok && isNew) ok &= ObjectMoleculeConnect(I, cset);
      if (ok && cset->Symmetry) {
        I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
      }
      if (I->Symmetry) {
        /* check scale records */
        if (pdb_info && pdb_info->scale.flag[0] && pdb_info->scale.flag[1] &&
            pdb_info->scale.flag[2]) {
          float* sca = pdb_info->scale.matrix;
          sca[15] = 1.0F;

          CoordSetInsureOrthogonal(G, cset, sca, &I->Symmetry->Crystal, quiet);
        }
      }
      SceneCountFrames(G);
      if (ok) ok &= ObjectMoleculeExtendIndices(I, state);
      if (ok) ok &= ObjectMoleculeSort(I);
      if (ok) {
        ObjectMoleculeUpdateIDNumbers(I);
        ObjectMoleculeUpdateNonbonded(I);
        ObjectMoleculeAutoDisableAtomNameWildcard(I);
      }
      if (SettingGetGlobal_b(G, cSetting_pdb_hetatm_guess_valences)) {
        ObjectMoleculeGuessValences(I, state, nullptr, nullptr, false);
      }

      successCnt++;
      if (!quiet) {
        if (successCnt > 1) {
          if (successCnt == 2) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
            " %s: read MODEL %d\n", __func__, 1 ENDFB(G);
          }
          PRINTFB(G, FB_ObjectMolecule, FB_Actions)
          " %s: read MODEL %d\n", __func__, successCnt ENDFB(G);
        }
      }
    }
    if (restart) {
      repeatFlag = true;
      start = restart;
      state = state + 1;
    }
  }
  if (!ok && isNew) {
    DeleteP(I);
  }
  return (I);
}
CoordSet* ObjectMoleculeMMDStr2CoordSet(PyMOLGlobals* G, const char* buffer,
                                        AtomInfoType** atInfoPtr,
                                        const char** restart) {
  const char* p;
  int nAtom, nBond;
  int a, c, bPart, bOrder;
  float* coord = nullptr;
  CoordSet* cset = nullptr;
  AtomInfoType *atInfo = nullptr, *ai;
  char cc[MAXLINELEN];
  WordType title;

  float* f;
  BondType *ii, *bond = nullptr;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);

  p = buffer;
  nAtom = 0;
  if (atInfoPtr) atInfo = *atInfoPtr;

  if (ok) {
    p = ncopy(cc, p, 6);  // could this be too small for large molecules with
                          // #atoms > 999999 ?
    if (sscanf(cc, "%d", &nAtom) != 1)
      ok = ErrMessage(G, "ReadMMDFile", "bad atom count");
  }

  if (ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if (atInfo) VLACheck(atInfo, AtomInfoType, nAtom);
  }

  if (!atInfo) {
    ErrFatal(G, "MMDStr2CoordSet", "need atom information record!");
    /* failsafe for old version.. */
  }

  nBond = 0;
  if (ok) {
    bond = VLACalloc(BondType, 6 * nAtom);
  }

  p = ParseWordCopy(title, p, sizeof(WordType) - 1);
  UtilCleanStr(title);

  p = nextline(p);

  /* read coordinates and atom names */

  if (ok) {
    f = coord;
    ii = bond;
    for (a = 0; a < nAtom; a++) {
      ai = atInfo + a;

      ai->id = a + 1;
      ai->rank = a;
      if (ok) {
        p = ncopy(cc, p, 4);
        if (sscanf(cc, "%d", &ai->customType) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad atom type");
      }
      if (ok) {
        if (ai->customType <= 14)
          strcpy(ai->elem, "C");
        else if (ai->customType <= 23)
          strcpy(ai->elem, "O");
        else if (ai->customType <= 40)
          strcpy(ai->elem, "N");
        else if (ai->customType <= 48)
          strcpy(ai->elem, "H");
        else if (ai->customType <= 52)
          strcpy(ai->elem, "S");
        else if (ai->customType <= 53)
          strcpy(ai->elem, "P");
        else if (ai->customType <= 55)
          strcpy(ai->elem, "B");
        else if (ai->customType <= 56)
          strcpy(ai->elem, "F");
        else if (ai->customType <= 57)
          strcpy(ai->elem, "Cl");
        else if (ai->customType <= 58)
          strcpy(ai->elem, "Br");
        else if (ai->customType <= 59)
          strcpy(ai->elem, "I");
        else if (ai->customType <= 60)
          strcpy(ai->elem, "Si");
        else if (ai->customType <= 61)
          strcpy(ai->elem, "Du");
        else if (ai->customType <= 62)
          strcpy(ai->elem, "Z0");
        else if (ai->customType <= 63)
          strcpy(ai->elem, "Lp");
        else
          ai->elem[0] = 0;
        /* else strcpy(ai->elem,"?"); WLD 090305 -- guess instead. */
      }
      for (c = 0; c < 6; c++) {
        if (ok) {
          p = ncopy(cc, p, 8);
          if (sscanf(cc, "%d%d", &bPart, &bOrder) != 2)
            ok = ErrMessage(G, "ReadMMDFile", "bad bond record");
          else {
            if (bPart && bOrder && (a < (bPart - 1))) {
              nBond++;
              ii->index[0] = a;
              ii->index[1] = bPart - 1;
              ii->order = bOrder;
              ii++;
            }
          }
        }
      }
      if (ok) {
        p = ncopy(cc, p, 12);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if (ok) {
        p = ncopy(cc, p, 12);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if (ok) {
        p = ncopy(cc, p, 12);
        if (sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if (ok) {
        p = nskip(p, 1);
        p = ncopy(cc, p, 5);
        ai->setResi(cc);

        // single-letter code
        p = nskip(p, 1);

        // chain
        p = ncopy(cc, p, 1);
        LexAssign(G, ai->chain, cc);
      }
      if (ok) {
        p = nskip(p, 4);
        p = ncopy(cc, p, 9);
        if (sscanf(cc, "%f", &ai->partialCharge) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad charge");
      }
      if (ok) {
        p = nskip(p, 10);
        p = ncopy(cc, p, 3);
        UtilCleanStr(cc);
        LexAssign(G, ai->resn, cc);
        ai->hetatm = true;
      }

      ai->segi = 0;
      ai->alt[0] = 0;

      if (ok) {
        p = nskip(p, 2);
        p = ntrim(cc, p, 4);
        if (!cc[0]) {
          sprintf(cc, "%s%02d", ai->elem, a + 1);
        }
        ai->name = LexIdx(G, cc);

        ai->visRep = auto_show;
      }
      if (ok) {
        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }
      if (!ok) break;
      p = nextline(p);
    }
  }
  if (ok) VLASize(bond, BondType, nBond);
  if (ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, title);
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if (atInfoPtr) *atInfoPtr = atInfo;
  *restart = *p ? p : nullptr;
  return (cset);
}
