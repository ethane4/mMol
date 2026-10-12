/**
 * Session (de)serialization of the fields ObjectMeshState and
 * ObjectSurfaceState share: list items 0-12. Each caller handles the rest.
 */

#pragma once

#include "pymol/base/word.h"
#include "pymol/platform/os_python.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/symmetry/crystal.h"

template <typename State>
void ObjectIsoStateHeadAsPyList(const State* I, PyObject* result) {
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  PyList_SetItem(result, 1, PyString_FromString(I->MapName));
  PyList_SetItem(result, 2, PyInt_FromLong(I->MapState));
  PyList_SetItem(result, 3, CrystalAsPyList(&I->Crystal));
  PyList_SetItem(result, 4, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 7, PConvIntArrayToPyList(I->Range, 6));
  PyList_SetItem(result, 8, PyFloat_FromDouble(I->Level));
  PyList_SetItem(result, 9, PyFloat_FromDouble(I->Radius));
  PyList_SetItem(result, 10, PyInt_FromLong(I->CarveFlag));
  PyList_SetItem(result, 11, PyFloat_FromDouble(I->CarveBuffer));
  if (I->CarveFlag && I->AtomVertex) {
    PyList_SetItem(result, 12, PConvFloatVLAToPyList(I->AtomVertex));
  } else {
    PyList_SetItem(result, 12, PConvAutoNone(nullptr));
  }
}

template <typename State>
int ObjectIsoStateHeadFromPyList(State* I, PyObject* list) {
  int ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
  if (ok) ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->MapName, WordLength);
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->MapState);
  if (ok) ok = CrystalFromPyList(&I->Crystal, PyList_GetItem(list, 3));
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->ExtentFlag);
  if (ok)
    ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 5), I->ExtentMin,
                                        3);
  if (ok)
    ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->ExtentMax,
                                        3);
  if (ok)
    ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 7), I->Range, 6);
  if (ok) ok = PConvPyFloatToFloat(PyList_GetItem(list, 8), &I->Level);
  if (ok) ok = PConvPyFloatToFloat(PyList_GetItem(list, 9), &I->Radius);
  if (ok) ok = PConvPyIntToInt(PyList_GetItem(list, 10), &I->CarveFlag);
  if (ok) ok = PConvPyFloatToFloat(PyList_GetItem(list, 11), &I->CarveBuffer);
  if (ok) {
    PyObject* tmp = PyList_GetItem(list, 12);
    if (tmp == Py_None)
      I->AtomVertex = nullptr;
    else
      ok = PConvPyListToFloatVLA(tmp, &I->AtomVertex);
  }
  return ok;
}
