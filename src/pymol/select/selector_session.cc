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

// Selections and colorections in sessions.

#include "pymol/select/selector.h"
#include "pymol/executive/executive.h"
#include "pymol/base/util2.h"
#include "pymol/select/selector_internal.h"

/* special selections, unknown to executive */
#define cSelectorSecretsPrefix "_!"
#define cColorectionFormat "_!c_%s_%d"

typedef struct {
  int color;
  SelectorID_t sele;
} ColorectionRec;

static void SelectorDeleteIndex(PyMOLGlobals* G, SelectorID_t index) {
  auto I = G->SelectorMgr;
  auto it = std::find_if(
      I->Info.begin() + 1, I->Info.end(),
      [index](const SelectionInfoRec& rec) { return rec.ID == index; });
  if (it != I->Info.end()) {
    SelectorDeleteSeleAtIter(G, it);
  }
}

PyObject* SelectorColorectionGet(PyMOLGlobals* G, const char* prefix) {
  auto I = G->Selector;
  auto IM = G->SelectorMgr;
  PyObject* result = nullptr;
  int n_used = 0;
  ColorectionRec *used = nullptr, tmp;
  ov_size a, b;
  int found;
  int color;
  AtomInfoType* ai;
  used = VLAlloc(ColorectionRec, 1000);

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    ai = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
    color = ai->color;
    found = false;
    for (b = 0; b < n_used; b++) {
      if (used[b].color == color) {
        tmp = used[0]; /* optimize to minimize N^2 effects */
        used[0] = used[b];
        used[b] = tmp;
        found = true;
        break;
      }
    }
    if (!found) {
      VLACheck(used, ColorectionRec, n_used);
      used[n_used] = used[0];
      used[0].color = color;
      n_used++;
    }
  }
  for (a = 0; a < n_used; a++) {
    /* create selections */

    SelectorID_t sele = IM->NSelection++;
    used[a].sele = sele;
    IM->Info.emplace_back(SelectionInfoRec(
        sele, pymol::string_format(cColorectionFormat, prefix, used[a].color)));
  }

  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    ai = I->Obj[I->Table[a].model]->AtomInfo + I->Table[a].atom;
    color = ai->color;
    for (b = 0; b < n_used; b++) {
      if (used[b].color == color) {
        tmp = used[0]; /* optimize to minimize N^2 effects */
        used[0] = used[b];
        used[b] = tmp;

        /* add selection onto atom */
        SelectorManagerInsertMember(*IM, *ai, used[0].sele);
        break;
      }
    }
  }

  VLASize(used, ColorectionRec, n_used * 2);
  result = PConvIntVLAToPyList((int*)used);
  VLAFreeP(used);
  return (result);
}

int SelectorColorectionApply(PyMOLGlobals* G, PyObject* list,
                             const char* prefix) {
  CSelector* I = G->Selector;
  int ok = true;
  ColorectionRec* used = nullptr;
  ov_size n_used = 0;
  int a, b;
  AtomInfoType* ai;
  ObjectMolecule *obj, *last = nullptr;

  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  if (ok) n_used = PyList_Size(list) / 2;
  if (ok) ok = ((used = VLAlloc(ColorectionRec, n_used)) != nullptr);
  if (ok) ok = PConvPyListToIntArrayInPlace(list, (int*)used, n_used * 2);
  if (ok) {
    SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

    for (b = 0; b < n_used; b++) { /* update selection indices */
      auto name =
          pymol::string_format(cColorectionFormat, prefix, used[b].color);
      used[b].sele = SelectorIndexByName(G, name.c_str());
    }

    for (a = cNDummyAtoms; a < I->Table.size(); a++) {
      obj = I->Obj[I->Table[a].model];
      ai = obj->AtomInfo + I->Table[a].atom;

      for (b = 0; b < n_used; b++) {
        if (SelectorIsMember(G, ai->selEntry, used[b].sele)) {
          ai->color = used[b].color;
          if (obj != last) {
            obj->invalidate(cRepAll, cRepInvColor, -1);
            last = obj;
          }
          break;
        }
      }
    }
  }
  VLAFreeP(used);
  return (ok);
}

int SelectorColorectionFree(PyMOLGlobals* G, PyObject* list,
                            const char* prefix) {
  int ok = true;
  ColorectionRec* used = nullptr;
  ov_size n_used = 0;
  ov_size b;

  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  if (ok) n_used = PyList_Size(list) / 2;
  if (ok) ok = ((used = VLAlloc(ColorectionRec, n_used)) != nullptr);
  if (ok) ok = PConvPyListToIntArrayInPlace(list, (int*)used, n_used * 2);
  if (ok) {
    for (b = 0; b < n_used; b++) { /* update selection indices */
      auto name =
          pymol::string_format(cColorectionFormat, prefix, used[b].color);
      used[b].sele = SelectorIndexByName(G, name.c_str());
    }

    for (b = 0; b < n_used; b++) {
      SelectorDeleteIndex(G, used[b].sele);
    }
  }
  VLAFreeP(used);
  return (ok);
}

PyObject* SelectorSecretsAsPyList(PyMOLGlobals* G) {
  auto I = G->SelectorMgr;

  auto n_secret = std::count_if(I->Info.begin(), I->Info.end(),  //
                                [](const SelectionInfoRec& rec) {
                                  return pymol::starts_with(
                                      rec.name, cSelectorSecretsPrefix);
                                });
  auto result = PyList_New(n_secret);
  n_secret = 0;
  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
  for (int a = 0; a < I->Info.size(); a++) {
    if (pymol::starts_with(I->Info[a].name, cSelectorSecretsPrefix)) {
      auto list = PyList_New(2);
      PyList_SetItem(list, 0, PyString_FromString(I->Info[a].name.c_str()));
      PyList_SetItem(list, 1, SelectorAsPyList(G, I->Info[a].ID));
      PyList_SetItem(result, n_secret, list);
      n_secret++;
    }
  }
  return (result);
}

int SelectorSecretsFromPyList(PyMOLGlobals* G, PyObject* list) {
  int ok = true;
  ov_size n_secret = 0;
  ov_size a;
  PyObject* entry = nullptr;
  std::string name;
  ov_size ll = 0;
  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  if (ok) n_secret = PyList_Size(list);
  if (ok) {
    for (a = 0; a < n_secret; a++) {
      if (ok) entry = PyList_GetItem(list, a);
      if (ok) ok = (entry != nullptr);
      if (ok) ok = PyList_Check(entry);
      if (ok) ll = PyList_Size(entry);
      if (ok & (ll > 1)) {
        if (ok) {
          ok = PConvFromPyListItem(G, entry, 0, name);
        }
        if (ok) {
          CPythonVal* val = CPythonVal_PyList_GetItem(G, entry, 1);
          ok = SelectorFromPyList(G, name.c_str(), val);
          CPythonVal_Free(val);
        }
      }
      if (!ok) break;
    }
  }
  return (ok);
}

typedef struct {
  int atom;
  int tag;
} SelAtomTag;

PyObject* SelectorAsPyList(
    PyMOLGlobals* G,
    SelectorID_t sele1) { /* assumes SelectorUpdateTable has been called */
  CSelector* I = G->Selector;
  int a, b;
  int at;
  int s;
  SelAtomTag** vla_list = nullptr;
  int n_obj = 0;
  int n_idx = 0;
  int cur = -1;
  ObjectMolecule** obj_list = nullptr;
  ObjectMolecule *obj, *cur_obj = nullptr;
  PyObject* result = nullptr;
  PyObject* obj_pyobj;
  PyObject* idx_pyobj;
  PyObject* tag_pyobj;

  vla_list = VLACalloc(SelAtomTag*, 10);
  obj_list = VLAlloc(ObjectMolecule*, 10);

  n_idx = 0;
  for (a = cNDummyAtoms; a < I->Table.size(); a++) {
    int tag;
    at = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];
    s = obj->AtomInfo[at].selEntry;
    if ((tag = SelectorIsMember(G, s, sele1))) {
      if (cur_obj != obj) {
        if (n_idx) {
          VLASize(vla_list[cur], SelAtomTag, n_idx);
        }
        cur++;
        VLACheck(vla_list, SelAtomTag*, n_obj);
        vla_list[cur] = VLAlloc(SelAtomTag, 1000);
        VLACheck(obj_list, ObjectMolecule*, n_obj);
        obj_list[cur] = obj;
        cur_obj = obj;
        n_obj++;
        n_idx = 0;
      }
      VLACheck(vla_list[cur], SelAtomTag, n_idx);
      vla_list[cur][n_idx].atom = at;
      vla_list[cur][n_idx].tag = tag;
      n_idx++;
    }
  }
  if (cur_obj) {
    if (n_idx) {
      VLASize(vla_list[cur], SelAtomTag, n_idx);
    }
  }
  if (n_obj) {
    result = PyList_New(n_obj);
    for (a = 0; a < n_obj; a++) {
      obj_pyobj = PyList_New(3);
      n_idx = VLAGetSize(vla_list[a]);
      idx_pyobj = PyList_New(n_idx);
      tag_pyobj = PyList_New(n_idx);
      for (b = 0; b < n_idx; b++) {
        PyList_SetItem(idx_pyobj, b, PyInt_FromLong(vla_list[a][b].atom));
        PyList_SetItem(tag_pyobj, b, PyInt_FromLong(vla_list[a][b].tag));
      }
      VLAFreeP(vla_list[a]);
      PyList_SetItem(obj_pyobj, 0, PyString_FromString(obj_list[a]->Name));
      PyList_SetItem(obj_pyobj, 1, idx_pyobj);
      PyList_SetItem(obj_pyobj, 2, tag_pyobj);
      PyList_SetItem(result, a, obj_pyobj);
    }
  } else {
    result = PyList_New(0);
  }
  VLAFreeP(vla_list);
  VLAFreeP(obj_list);
  return (result);
}

int SelectorFromPyList(PyMOLGlobals* G, const char* name, PyObject* list) {
  int ok = true;
  auto I = G->SelectorMgr;
  ov_size a, b;
  ov_size ll;
  PyObject* obj_list = nullptr;
  PyObject *idx_list = nullptr, *tag_list;
  ov_size n_obj = 0, n_idx = 0;
  int idx, tag;
  const char* oname;
  ObjectMolecule* obj;
  int singleAtomFlag = true;
  int singleObjectFlag = true;
  ObjectMolecule* singleObject = nullptr;
  int singleAtom = -1;

  if (ok) ok = PyList_Check(list);
  if (ok) n_obj = PyList_Size(list);

  /* get rid of existing selection */
  SelectorDelete(G, name);

  int sele = I->NSelection++;
  I->Info.emplace_back(SelectionInfoRec(sele, name));
  if (ok) {
    for (a = 0; a < n_obj; a++) {
      ll = 0;
      if (ok) obj_list = PyList_GetItem(list, a);
      if (ok) ok = PyList_Check(obj_list);
      if (ok) ll = PyList_Size(obj_list);
      if (ok) ok = PConvPyStrToStrPtr(PyList_GetItem(obj_list, 0), &oname);
      obj = nullptr;
      if (ok) obj = ExecutiveFindObjectMoleculeByName(G, oname);
      if (ok && obj) {
        if (ok) idx_list = PyList_GetItem(obj_list, 1);
        if (ll > 2)
          tag_list = PyList_GetItem(obj_list, 2);
        else
          tag_list = nullptr;
        if (ok) ok = PyList_Check(idx_list);
        if (ok) n_idx = PyList_Size(idx_list);
        for (b = 0; b < n_idx; b++) {
          if (ok) ok = PConvPyIntToInt(PyList_GetItem(idx_list, b), &idx);
          if (tag_list)
            PConvPyIntToInt(PyList_GetItem(tag_list, b), &tag);
          else
            tag = 1;
          if (ok && (idx < obj->NAtom)) {
            SelectorManagerInsertMember(*I, obj->AtomInfo[idx], sele, tag);

            /* take note of selections which are one atom/one object */
            if (singleObjectFlag) {
              if (singleObject) {
                if (obj != singleObject) {
                  singleObjectFlag = false;
                }
              } else {
                singleObject = obj;
              }
            }

            if (singleAtomFlag) {
              if (singleAtom >= 0) {
                if (idx != singleAtom) {
                  singleAtomFlag = false;
                }
              } else {
                singleAtom = idx;
              }
            }
          }
        }
      }
    }
    { /* make note of single atom/object selections */
      auto& info = I->Info.back();
      if (singleObjectFlag && singleObject) {
        info.theOneObject = singleObject;
        if (singleAtomFlag && (singleAtom >= 0)) {
          info.theOneAtom = singleAtom;
        }
      }
    }
  }
  return (ok);
}
