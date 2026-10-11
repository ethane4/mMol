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

// Helpers shared by the executive_*.cc files; not part of the executive API.

#pragma once

#include "pymol/executive/executive.h"

/**
 * Macro with error handling to get a selection ID from a non-empty expression.
 * @param expression Selection expression
 * @param[out] tmpsele Name of the pymol::Result<SelectorTmp> instance
 * @param[out] seleID Name of the SelectorID_t instance
 */
#define SETUP_SELE(expression, tmpsele, seleID)                                \
  auto tmpsele = SelectorTmp::make(G, expression);                             \
  p_return_if_error(tmpsele);                                                  \
  SelectorID_t seleID = tmpsele->getIndex();                                   \
  if (seleID < 0) {                                                            \
    return pymol::make_error("This should not happen - PyMOL may have a bug"); \
  }                                                                            \
  assert(seleID != cSelectionInvalid)

/**
 * Macro with error handling to get a selection ID from a non-empty expression
 * `sN`. Assigns `tmpseleN` and `seleN`.
 */
#define SETUP_SELE_DEFAULT(N) SETUP_SELE(s##N, tmpsele##N, sele##N)

/**
 * Macro with error handling to get a selection ID from a non-empty expression
 * `sN`. Assigns `tmpseleN` and `seleN`.
 * Adds "Selection N: " prefix to error messages.
 * @param same_value If `sN` is "same" then use this value for `seleN`
 */
#define SETUP_SELE_DEFAULT_PREFIXED(N, same_value)                \
  pymol::Result<SelectorTmp> tmpsele##N;                          \
  SelectorID_t sele##N = same_value;                              \
  if (!WordMatchExact(G, s##N, cKeywordSame, true)) {             \
    tmpsele##N = SelectorTmp::make(G, s##N);                      \
    p_return_if_error_prefixed(tmpsele##N, "Selection " #N ": "); \
    sele##N = (tmpsele##N)->getIndex();                           \
  }                                                               \
  if (sele##N == cSelectionInvalid) {                             \
    return pymol::make_error("Invalid selection " #N);            \
  }

/// Enables `rec` and, if `parents`, the groups that contain it; logs the
/// command if `log` and logging is on.
void ExecutiveSpecEnable(PyMOLGlobals* G, SpecRec* rec, int parents, int log);

/**
 * Find object of given type, or delete object if it exists but has the wrong
 * type.
 * @param name Object name
 * @return nullptr if object can't be found or has the wrong type
 */
template <typename ObjectT>
ObjectT* ExecutiveFindOrDeleteObject(PyMOLGlobals* G,
                                     pymol::zstring_view name) {
  auto anyObj = ExecutiveFindObjectByName(G, name.c_str());
  auto obj = dynamic_cast<ObjectT*>(anyObj);
  if (anyObj && !obj) {
    // incompatible object with the same name
    ExecutiveDelete(G, name.c_str());
  }
  return obj;
}

#define cExecExpandKeepGroups 2

/// Redraws after `rec` was enabled or disabled.
void ReportEnabledChange(PyMOLGlobals* G, SpecRec* rec);

/// Indexes `rec` by name in the executive's lexicon; false on failure.
int ExecutiveAddKey(CExecutive* I, SpecRec* rec);

/// Number of entries (objects and selections) in the object list.
int ExecutiveCountNames(PyMOLGlobals* G);

/// Assigns grid-mode slots to the objects, unless they are current and
/// `force` is false.
void ExecutiveUpdateGridSlots(PyMOLGlobals* G, int force);
