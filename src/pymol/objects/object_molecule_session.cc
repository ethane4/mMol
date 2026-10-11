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

// Converting molecules, atoms, bonds and coordinate sets to and from session
// lists.

#include "pymol/base/err.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/base/feedback.h"
#include "pymol/python/sym_op_pconv.h"
#include "pymol/objects/atom_info_history.h"
#include "pymol/objects/bond_type_history.h"
#include "pymol/objects/property.h"
#include "pymol/objects/object_molecule_internal.h"

static PyObject* ObjectMoleculeCSetAsPyList(ObjectMolecule* I) {
  PyObject* result = nullptr;
  int a;
  result = PyList_New(I->NCSet);
  for (a = 0; a < I->NCSet; a++) {
    if (I->CSet[a]) {
      PyList_SetItem(result, a, CoordSetAsPyList(I->CSet[a]));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(Py_None));
    }
  }
  return (PConvAutoNone(result));
}

/*static PyObject *ObjectMoleculeDiscreteCSetAsPyList(ObjectMolecule *I)
  {
  PyObject *result = nullptr;
  return(PConvAutoNone(result));
  }*/
static int ObjectMoleculeCSetFromPyList(ObjectMolecule* I, PyObject* list) {
  int ok = true;
  int a;
  if (ok) ok = PyList_Check(list);
  if (ok) {
    VLACheck(I->CSet, CoordSet*, I->NCSet);
    for (a = 0; a < I->NCSet; a++) {
      if (ok)
        ok = CoordSetFromPyList(I->G, PyList_GetItem(list, a), &I->CSet[a]);
      PRINTFB(I->G, FB_ObjectMolecule, FB_Debugging)
      " %s: ok %d after CoordSet %d\n", __func__, ok, a ENDFB(I->G);

      if (ok)
        if (I->CSet[a]) /* WLD 030205 */
          I->CSet[a]->Obj = I;
    }
  }
  return (ok);
}

static PyObject* ObjectMoleculeBondAsPyList(ObjectMolecule* I) {
  PyObject* result = nullptr;
  PyObject* bond_list;
  const BondType* bond;
  int a;

  PyMOLGlobals* G = I->G;
  int pse_export_version =
      SettingGetGlobal_f(I->G, cSetting_pse_export_version) * 1000;

  if (SettingGetGlobal_b(G, cSetting_pse_binary_dump) &&
      (!pse_export_version || pse_export_version >= 1765)) {
    /* For the pse_binary_dump, save entire Bond array to a binary string array
     */

    // supported versions
    auto version = (!pse_export_version || pse_export_version >= 1810) ? 181
                   : (pse_export_version >= 1770)                      ? 177
                                                                       : 176;

    auto blob = Copy_To_BondType_Version(version, I->Bond.data(), I->NBond);
    auto blobsize = VLAGetByteSize(blob);

    result = PyList_New(2);
    PyList_SetItem(result, 0, PyInt_FromLong(version));
    PyList_SetItem(result, 1,
                   PyBytes_FromStringAndSize(
                       reinterpret_cast<const char*>(blob), blobsize));

    VLAFreeP(blob);

    return result;
  }
  result = PyList_New(I->NBond);
  bond = I->Bond;
  for (a = 0; a < I->NBond; a++) {
    size_t const list_size = bond->hasSymOp() ? 8 : 7;
    bond_list = PyList_New(list_size);
    PyList_SetItem(bond_list, 0, PyInt_FromLong(bond->index[0]));
    PyList_SetItem(bond_list, 1, PyInt_FromLong(bond->index[1]));
    PyList_SetItem(bond_list, 2, PyInt_FromLong(bond->order));
    PyList_SetItem(bond_list, 3, PyInt_FromLong(-1));  // id
    PyList_SetItem(bond_list, 4, PyInt_FromLong(0));   // stereo
    PyList_SetItem(bond_list, 5, PyInt_FromLong(bond->unique_id));
    PyList_SetItem(bond_list, 6, PyInt_FromLong(bond->has_setting));
    if (list_size > 7) {
      PyList_SetItem(bond_list, 7, PConvToPyObject(bond->symop_2));
    }
    PyList_SetItem(result, a, bond_list);
    bond++;
  }

  return (PConvAutoNone(result));
}

static int ObjectMoleculeBondFromPyList(ObjectMolecule* I, PyObject* list) {
  PyMOLGlobals* G = I->G;
  int ok = true;
  int a;
  int stereo, ll = 0;
  PyObject* bond_list = nullptr;
  BondType* bond;

  if (ok) ok = PyList_Check(list);
  if (ok) ll = PyList_Size(list);

  bool pse_binary_dump = false;

  if (ll >= 2) {
    // checking if from pse_binary_dump
    // pse_binary_dump saves 2 values: bondInfo_version, BondType binary
    CPythonVal* val1 = CPythonVal_PyList_GetItem(G, list, 1);
    pse_binary_dump = PyBytes_Check(val1);
    CPythonVal_Free(val1);
  }
  if (pse_binary_dump) {
    CPythonVal* verobj = CPythonVal_PyList_GetItem(G, list, 0);
    int bondInfo_version;
    ok = PConvPyIntToInt(verobj, &bondInfo_version);

    CPythonVal* strobj = CPythonVal_PyList_GetItem(G, list, 1);
    auto strval = PyBytes_AsSomeString(strobj);

    if (ok) ok = bool((I->Bond = pymol::vla<BondType>(I->NBond)));

    Copy_Into_BondType_From_Version(strval.data(), bondInfo_version,
                                    I->Bond.data(), I->NBond);

    CPythonVal_Free(verobj);
    CPythonVal_Free(strobj);
  } else {
    if (ok) ok = bool((I->Bond = pymol::vla<BondType>(I->NBond)));
    bond = I->Bond.data();
    for (a = 0; a < I->NBond; a++) {
      bond_list = nullptr;
      if (ok) bond_list = PyList_GetItem(list, a);
      if (ok) ok = PyList_Check(bond_list);
      if (ok) ll = PyList_Size(bond_list);
      if (ok)
        ok = PConvPyIntToInt(PyList_GetItem(bond_list, 0), &bond->index[0]);
      if (ok)
        ok = PConvPyIntToInt(PyList_GetItem(bond_list, 1), &bond->index[1]);
      if (ok)
        if ((ok = CPythonVal_PConvPyIntToInt_From_List(I->G, bond_list, 2,
                                                       &stereo)))
          bond->order = stereo;
      if (ok && (ll > 5)) {
        int has_setting;
        if (ok)
          ok = PConvPyIntToInt(PyList_GetItem(bond_list, 5), &bond->unique_id);
        if (ok)
          ok = PConvPyIntToInt(PyList_GetItem(bond_list, 6), &has_setting);
        if (ok) bond->has_setting = (short int)has_setting;
        if (ok && bond->unique_id) { /* reserve existing IDs */
          bond->unique_id =
              SettingUniqueConvertOldSessionID(G, bond->unique_id);
        }
      }
      if (ll > 7) {
        PConvFromPyListItem(G, bond_list, 7, bond->symop_2);
      }
      bond++;
      CPythonVal_Free(bond_list);
    }
  }
  PRINTFB(G, FB_ObjectMolecule, FB_Debugging)
  " %s: ok %d after restore\n", __func__, ok ENDFB(G);

  return (ok);
}

/**
 * Extract an atom property "column" as a Python list.
 *
 * As an optimization, trailing None values are removed, so the returned list
 * may be shorter than the atom array or even be empty.
 *
 * @param mol Object molecule which provides the atom array
 * @param func Function which extracts a property from an atom
 * @return List of properties
 */
static PyObject* AtomColumnAsPyList(
    const ObjectMolecule& mol,
    std::function<PyObject*(const AtomInfoType&)> func) {
  auto list = PyList_New(mol.NAtom);
  PyObject* prev = Py_None;
  int pruned_size = 0;

  for (int i = 0; i < mol.NAtom; ++i) {
    PyObject* value = func(mol.AtomInfo[i]);

    // Simple optimization: Repeated property lists will reference the same
    // Python object
    if (prev != value && PyObject_RichCompareBool(prev, value, Py_EQ)) {
      Py_INCREF(prev);
      Py_DECREF(value);
      value = prev;
    } else {
      prev = value;
    }

    if (value != Py_None) {
      pruned_size = i + 1;
    }

    PyList_SetItem(list, i, value);
  }

  assert(pruned_size <= mol.NAtom);

  // Simple optimization: Prune "None" tail
  PyList_SetSlice(list, pruned_size, mol.NAtom, nullptr);

  assert(PyList_Size(list) == pruned_size);

  return list;
}

/**
 * Restore an atom property "column" from a Python list.
 * @param mol Object molecule to update atoms in
 * @param list Property list (may be shorter than atom array)
 * @param func Function which sets a property for an atom
 */
static void AtomColumnFromPyList(
    ObjectMolecule& mol, PyObject* list,
    std::function<void(AtomInfoType&, PyObject*)> func) {
  if (!list || !PyList_Check(list)) {
    return;
  }

  int size = PyList_Size(list);
  assert(size <= mol.NAtom);

  for (int i = 0; i < size; ++i) {
    func(mol.AtomInfo[i], PyList_GetItem(list, i));
  }
}

static PyObject* ObjectMoleculeAtomAsPyList(ObjectMolecule* I) {
  PyMOLGlobals* G = I->G;
  PyObject* result = nullptr;
  const AtomInfoType* ai;
  int a;
  int pse_export_version =
      SettingGetGlobal_f(I->G, cSetting_pse_export_version) * 1000;

  if (SettingGetGlobal_b(G, cSetting_pse_binary_dump) &&
      (!pse_export_version || pse_export_version >= 1765)) {
    /* For the pse_binary_dump, record all strings in lex and
       write them into separate binary string
     */
    AtomInfoTypeConverter converter(G, I->NAtom);
    std::set<lexidx_t> lexIDs;
    int totalstlen = 0;
    ai = I->AtomInfo.data();
    for (a = 0; a < I->NAtom; a++) {
      if (ai->textType) lexIDs.insert(ai->textType);
      if (ai->chain) lexIDs.insert(ai->chain);
      if (ai->label) lexIDs.insert(ai->label);
      if (ai->custom) lexIDs.insert(ai->custom);
      if (ai->segi) lexIDs.insert(ai->segi);
      if (ai->resn) lexIDs.insert(ai->resn);
      if (ai->name) lexIDs.insert(ai->name);
      ++ai;
    }
    for (const auto& lexID :
         lexIDs) {  // need to calculate totalstlen so we can allocate
      const char* lexstr = LexStr(G, lexID);
      int lexlen = strlen(lexstr);
      totalstlen += lexlen + 1;
    }
    int strinfolen = totalstlen + sizeof(int) * (lexIDs.size() + 1);
    void* strinfo = pymol::malloc<unsigned char>(strinfolen);
    int* strval = (int*)strinfo;
    *(strval++) =
        lexIDs.size();  // first write number of strings into binary data string
    char* strpl = (char*)((char*)strinfo + (1 + lexIDs.size()) * sizeof(int));
    /* write map of lex ids and strings into binary data string as an array of
       ids and null-terminated strings */
    for (const auto& lexID : lexIDs) {
      *(strval++) = converter.to_lexidx_int(lexID);
      const char* strptr = LexStr(G, lexID);
      strcpy(strpl, strptr);
      strpl += strlen(strptr) + 1;
    }

    auto version = AtomInfoVERSION;
    if (pse_export_version && pse_export_version < 1810) {
      if (pse_export_version < 1770) {
        version = 176;
      } else {
        version = 177;
      }
    }

    auto blob = converter.allocCopy(version, I->AtomInfo);
    auto blobsize = VLAGetByteSize(blob);

    // PyMOL versions up to 2.3.5 can only restore list size 3
    size_t result_size = 3;

    // Atom properties (not binary)
    PyObject* prop_list = nullptr;
    if (pse_export_version > 2399) {
      prop_list = AtomColumnAsPyList(*I, [G](const AtomInfoType& atom) {
        return PConvAutoNone(
            atom.prop_id ? PropertyAsPyList(G, atom.prop_id, true) : nullptr);
      });

      result_size = 4;
    }

    result = PyList_New(result_size);
    PyList_SetItem(result, 0, PyInt_FromLong(version));
    PyList_SetItem(result, 1,
                   PyBytes_FromStringAndSize(
                       reinterpret_cast<const char*>(blob), blobsize));
    PyList_SetItem(result, 2,
                   PyBytes_FromStringAndSize(
                       reinterpret_cast<const char*>(strinfo), strinfolen));

    if (result_size > 3) {
      PyList_SetItem(result, 3, PConvAutoNone(prop_list));
    }

    VLAFreeP(blob);
    FreeP(strinfo);
    return result;
  }
  result = PyList_New(I->NAtom);
  ai = I->AtomInfo;
  for (a = 0; a < I->NAtom; a++) {
    PyList_SetItem(result, a, AtomInfoAsPyList(I->G, ai));
    ai++;
  }
  return (PConvAutoNone(result));
}

static int ObjectMoleculeAtomFromPyList(ObjectMolecule* I, PyObject* list) {
  PyMOLGlobals* G = I->G;
  int ok = true;
  int a, ll = 0;
  AtomInfoType* ai;

  if (ok) ok = PyList_Check(list);
  if (ok) ll = PyList_Size(list);

  bool pse_binary_dump = false;

  if (ll >= 3) {
    // checking if from pse_binary_dump
    // pse_binary_dump saves 3 values: atomInfo_version, AtomInfo binary, and
    // strings array
    CPythonVal* val1 = CPythonVal_PyList_GetItem(G, list, 1);
    CPythonVal* val2 = CPythonVal_PyList_GetItem(G, list, 2);
    pse_binary_dump = PyBytes_Check(val1) && PyBytes_Check(val2);
    CPythonVal_Free(val1);
    CPythonVal_Free(val2);
  }
  if (pse_binary_dump) {
    CPythonVal* verobj = CPythonVal_PyList_GetItem(G, list, 0);
    int atomInfo_version;
    ok = PConvPyIntToInt(verobj, &atomInfo_version);

    CPythonVal* strlookupobj = CPythonVal_PyList_GetItem(G, list, 2);
    auto strval_1 = PyBytes_AsSomeString(strlookupobj);
    int* strval = (int*)strval_1.data();

    AtomInfoTypeConverter converter(G, I->NAtom);

    auto& oldIDtoLexID = converter.lexidxmap;
    int nstrings = *(strval++);
    char* strpl = (char*)(strval + nstrings);
    int strcnt = nstrings;
    int stlen;
    // populate oldIDtoLexID with nstrings from binary string data (3rd entry in
    // list)
    while (strcnt) {
      lexidx_t idx =
          LexIdx(G, strpl);  // increments ref count, need to take into account
      int oldidx = *(strval++);
      oldIDtoLexID[oldidx] = idx;
      stlen = strlen(strpl);
      strpl += stlen + 1;
      strcnt--;
    }

    CPythonVal* strobj = CPythonVal_PyList_GetItem(G, list, 1);
    auto strval_2 = PyBytes_AsSomeString(strobj);

    VLACheck(I->AtomInfo, AtomInfoType, I->NAtom + 1);
    converter.copy(I->AtomInfo.data(), strval_2.data(), atomInfo_version);

    // go through AtomInfo array, swap new strings, convert colors, convert
    // settings (everything that AtomInfoFromPyList does except set properties,
    // which are currently
    //  not saved for pse_binary_dump)
    AtomInfoType* ai = I->AtomInfo.data();
    for (a = 0; a < I->NAtom; ++a, ++ai) {
      ai->color = ColorConvertOldSessionIndex(G, ai->color);
      if (ai->unique_id) {
        ai->unique_id = SettingUniqueConvertOldSessionID(G, ai->unique_id);
      }
    }
    // need to decrement since we call LexIdx() above on each
    for (auto it = oldIDtoLexID.begin(); it != oldIDtoLexID.end(); ++it) {
      LexDec(G, it->second);
    }
    CPythonVal_Free(verobj);
    CPythonVal_Free(strobj);
    CPythonVal_Free(strlookupobj);

    if (ll > 3) {
      // Restore atom properties
      AtomColumnFromPyList(*I, PyList_GetItem(list, 3),  //
                           [G](AtomInfoType& atom, PyObject* value) {
                             assert(atom.prop_id == 0);
                             atom.prop_id = PropertyFromPyList(G, value);
                           });
    }

  } else {
    // The old slow way of loading in AtomInfo, using python lists
    if (ok) VLACheck(I->AtomInfo, AtomInfoType, I->NAtom + 1);
    CHECKOK(ok, I->AtomInfo);
    ai = I->AtomInfo.data();
    for (a = 0; ok && a < I->NAtom; a++) {
      PyObject* val = PyList_GetItem(list, a);
      ok &= AtomInfoFromPyList(I->G, ai, val);
      ai++;
    }
  }
  PRINTFB(I->G, FB_ObjectMolecule, FB_Debugging)
  " %s: ok %d \n", __func__, ok ENDFB(I->G);
  return (ok);
}

int ObjectMoleculeNewFromPyList(PyMOLGlobals* G, PyObject* list,
                                ObjectMolecule** result) {
  int ok = true;
  ObjectMolecule* I = nullptr;
  int discrete_flag = 0;
  (*result) = nullptr;

  if (ok) ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 8), &discrete_flag);

  if (ok) I = new ObjectMolecule(G, discrete_flag);
  CHECKOK(ok, I);

  if (ok) {
    PyObject* val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->NCSet);
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->NBond);
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 3), &I->NAtom);
  if (ok) ok = ObjectMoleculeCSetFromPyList(I, PyList_GetItem(list, 4));
  if (ok) {
    ok = CoordSetFromPyList(G, PyList_GetItem(list, 5), &I->CSTmpl);

    if (I->CSTmpl) I->CSTmpl->Obj = I;
  }
  if (ok) {
    CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 6);
    ok = ObjectMoleculeBondFromPyList(I, val);
    CPythonVal_Free(val);
  }
  if (!ok && I) I->NBond = 0;
  if (ok) {
    CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 7);
    ok = ObjectMoleculeAtomFromPyList(I, val);
    CPythonVal_Free(val);
  }
  if (!ok && I) I->NAtom = 0;
  if (ok) {
    CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 10);
    I->Symmetry.reset(SymmetryNewFromPyList(G, val));
    CPythonVal_Free(val);
  }
  /* 11 was CurCSet */
  /* 12 was BondCounter */
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->AtomCounter);

  I->updateAtmToIdx();

  if (ok) I->invalidate(cRepAll, cRepInvAll, -1);
  if (ok)
    (*result) = I;
  else {
    /* cleanup */
    if (I) DeleteP(I);
    (*result) = nullptr;
  }
  return (ok);
}
PyObject* ObjectMoleculeAsPyList(ObjectMolecule* I) {
  PyObject* result = nullptr;

  /* first, dump the atoms */

  result = PyList_New(16);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NCSet));
  PyList_SetItem(result, 2, PyInt_FromLong(I->NBond));
  PyList_SetItem(result, 3, PyInt_FromLong(I->NAtom));
  PyList_SetItem(result, 4, ObjectMoleculeCSetAsPyList(I));
  PyList_SetItem(result, 5, CoordSetAsPyList(I->CSTmpl));
  PyList_SetItem(result, 6, ObjectMoleculeBondAsPyList(I));
  PyList_SetItem(result, 7, ObjectMoleculeAtomAsPyList(I));
  PyList_SetItem(result, 8, PyInt_FromLong(I->DiscreteFlag));
  PyList_SetItem(
      result, 9,
      PyInt_FromLong(I->DiscreteFlag ? I->NAtom : 0 /* NDiscrete */));
  PyList_SetItem(result, 10, SymmetryAsPyList(I->Symmetry.get()));
  PyList_SetItem(result, 11, PyInt_FromLong(0 /* CurCSet */));
  PyList_SetItem(result, 12, PyInt_FromLong(-1 /* BondCounter */));
  PyList_SetItem(result, 13, PyInt_FromLong(I->AtomCounter));

  float pse_export_version =
      SettingGetGlobal_f(I->G, cSetting_pse_export_version);

  if (I->DiscreteFlag &&
      (pse_export_version ||
       !SettingGetGlobal_b(I->G, cSetting_pse_binary_dump)) &&
      pse_export_version < 1.7699) {
    int* dcs;
    int a;
    CoordSet* cs;

    /* get coordinate set indices */

    for (a = 0; a < I->NCSet; a++) {
      cs = I->CSet[a];
      if (cs) {
        cs->tmp_index = a;
      }
    }

    dcs = pymol::malloc<int>(I->NAtom);

    for (a = 0; a < I->NAtom; a++) {
      cs = I->DiscreteCSet[a];
      if (cs)
        dcs[a] = cs->tmp_index;
      else
        dcs[a] = -1;
    }

    PyList_SetItem(result, 14,
                   PConvIntArrayToPyList(I->DiscreteAtmToIdx, I->NAtom));
    PyList_SetItem(result, 15, PConvIntArrayToPyList(dcs, I->NAtom));
    FreeP(dcs);
  } else {
    PyList_SetItem(result, 14, PConvAutoNone(nullptr));
    PyList_SetItem(result, 15, PConvAutoNone(nullptr));
  }

  return (PConvAutoNone(result));
}
