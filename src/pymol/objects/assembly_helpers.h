/*
 * Biological Assembly Helpers
 *
 * (c) 2017 Schrodinger, Inc.
 */

#pragma once

#include <set>

#include "pymol/objects/atom_info.h"
#include "pymol/objects/coord_set.h"
#include "pymol/base/lex.h"
#include "pymol/objects/object_molecule.h"

CoordSet* CoordSetCopyFilterChains(const CoordSet* other,
                                   const AtomInfoType* atInfo,
                                   const std::set<lexborrow_t>& chains_set);

void ObjectMoleculeSetAssemblyCSets(ObjectMolecule* I,
                                    CoordSet** assembly_csets);
