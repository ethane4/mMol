#pragma once

#include "pymol/objects/object_molecule.h"

void ObjectMoleculePBCUnwrap(ObjectMolecule&, bool bymol = true);
void ObjectMoleculePBCWrap(ObjectMolecule&, float const* center = nullptr);
