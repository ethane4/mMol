/*
 * Molecular file formats export
 *
 * (c) 2016 Schrodinger, Inc.
 */

#pragma once

#include "pymol/platform/os_python.h"
#include "pymol/base/vla.h"

#include "pymol/base/pymol_globals.h"

pymol::vla<char> MoleculeExporterGetStr(
    PyMOLGlobals* G, const char* format, const char* sele = "all",
    int state = cStateCurrent,  // current (-1 in Python API)
    const char* ref_object = "", int ref_state = cStateAll, int multi = -1,
    bool quiet = true);

PyObject* MoleculeExporterGetPyBonds(PyMOLGlobals* G, const char* selection,
                                     int state);
