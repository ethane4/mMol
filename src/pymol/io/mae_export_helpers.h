/*
 * MAE format export helper functions
 *
 * (c) 2018 Schrodinger, Inc.
 */

#pragma once

#include <string>

struct PyMOLGlobals;
class SeleCoordIterator;
struct AtomInfoType;

namespace pymol {
struct CObject;
}

int MaeExportGetAtomStyle(PyMOLGlobals* G, const SeleCoordIterator& iter);

int MaeExportGetBondStyle(const AtomInfoType* ai1, const AtomInfoType* ai2);

int MaeExportGetRibbonStyle(const AtomInfoType* ai);

void MaeExportGetRibbonColor(PyMOLGlobals* G, const SeleCoordIterator& iter,
                             char* ribbon_color_rgb);

std::string MaeExportGetLabelUserText(PyMOLGlobals* G, const AtomInfoType* ai);

std::string MaeExportGetSubGroupId(PyMOLGlobals* G, const pymol::CObject* obj);

std::string MaeExportStrRepr(const char* text);
