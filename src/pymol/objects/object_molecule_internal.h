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

// Helpers shared by the object_molecule*.cc files; not part of the
// ObjectMolecule API.

#pragma once

#include "pymol/select/selector.h"

// Short names for the line parsers in base/parse.h.
#define ntrim ParseNTrim
#define nextline ParseNextLine
#define ncopy ParseNCopy

// Defined in object_molecule_read.cc.
CoordSet* ObjectMoleculeMMDStr2CoordSet(PyMOLGlobals* G, const char* buffer,
                                        AtomInfoType** atInfoPtr,
                                        const char** restart);

// Defined in object_molecule.cc.
void ObjectMoleculeTransformTTTf(ObjectMolecule* I, float* ttt, int state);

// Defined in object_molecule_chem.cc.
int ObjectMoleculeGetAtomGeometry(const ObjectMolecule* I, int state, int at);
void ObjectMoleculeInferHBondFromChem(ObjectMolecule* I);

// Defined in object_molecule_edit.cc.
int ObjectMoleculeFixSeleHydrogens(ObjectMolecule* I, int sele, int state);

// Defined in object_molecule_pdb.cc.
void assign_pdb_known_residue(PyMOLGlobals* G, AtomInfoType* ai1,
                              AtomInfoType* ai2, int* bond_order);
