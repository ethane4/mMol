/*
 * This file contains source code for the PyMOL computer program
 * Copyright (c) Schrodinger, LLC.
 */

#pragma once

#include "pymol/base/pymol_globals.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/objects/coord_set.h"
#include "pymol/objects/atom_info.h"

void SideChainHelperMarkNonCartoonBonded(bool* marked,
                                         const ObjectMolecule* obj,
                                         const CoordSet* cs,
                                         bool cartoon_side_chain_helper,
                                         bool ribbon_side_chain_helper);

bool SideChainHelperFilterBond(PyMOLGlobals* G, const bool* marked,
                               const AtomInfoType* ati1,
                               const AtomInfoType* ati2, int b1, int b2,
                               int na_mode, int* c1, int* c2);
