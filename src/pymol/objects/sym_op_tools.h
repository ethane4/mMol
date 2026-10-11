/**
 * @file Tools for manipulating SymOp bonds
 *
 * (c) Schrodinger, Inc.
 */

#pragma once

#include "pymol/base/pymol_globals.h"
#include "pymol/symmetry/sym_op.h"

#include <vector>

struct CoordSet;

namespace pymol {
std::vector<std::pair<SymOp, float>> find_bond_symops(CoordSet const&, unsigned,
                                                      unsigned,
                                                      float tolerance = 1.05);
}  // namespace pymol
