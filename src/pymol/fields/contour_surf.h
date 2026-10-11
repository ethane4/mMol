/**
 * @file
 * Isosurface with VTKm
 *
 * (c) 2020 Schrodinger, Inc.
 */

#pragma once

#include "pymol/base/pymol_enums.h"
#include "pymol/base/vla.h"

struct PyMOLGlobals;
struct Isofield;
class CarveHelper;

int ContourSurfVolume(PyMOLGlobals* G, Isofield* field, float level,
                      pymol::vla<int>& num,     //
                      pymol::vla<float>& vert,  //
                      const int* range,         //
                      cIsosurfaceMode mode,     //
                      const CarveHelper*,       //
                      cIsosurfaceSide side);
