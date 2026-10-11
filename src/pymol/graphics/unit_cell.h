// Drawing a crystal's unit cell.

#pragma once

class CGO;
struct CCrystal;

/// The unit cell's 12 edges as a line CGO, or null if it couldn't be built.
CGO* CrystalGetUnitCellCGO(const CCrystal* I);
