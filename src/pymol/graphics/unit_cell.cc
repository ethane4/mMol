#include "pymol/graphics/unit_cell.h"

#include "pymol/base/pymol_globals.h"
#include "pymol/geometry/matrix.h"
#include "pymol/graphics/cgo.h"
#include "pymol/settings/setting.h"
#include "pymol/symmetry/crystal.h"

static float unitCellVertices[][3] = {
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},   // bottom 4 vertices
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};  // top 4 vertices
static float const unitCellVerticesCentered[][3] = {
    {-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
    {-0.5, -0.5, 0.5},  {0.5, -0.5, 0.5},  {0.5, 0.5, 0.5},  {-0.5, 0.5, 0.5},
};

static int unitCellLineIndices[] = {
    0, 1, 1, 2, 2, 3, 3, 0,   // bottom 4 lines
    4, 5, 5, 6, 6, 7, 7, 4,   // top 4 lines
    0, 4, 1, 5, 2, 6, 3, 7};  // 4 connector lines

CGO* CrystalGetUnitCellCGO(const CCrystal* I) {
  PyMOLGlobals* G = I->G;
  float v[3];
  CGO* cgo = nullptr;

  auto const ucv = SettingGet<bool>(G, cSetting_cell_centered)
                       ? unitCellVerticesCentered
                       : unitCellVertices;

  if (I) {
    cgo = CGONew(G);
    CGODisable(cgo, GL_LIGHTING);

    float* vertexVals = CGODrawArrays(cgo, GL_LINES, CGO_VERTEX_ARRAY, 24);
    for (int pl = 0; pl < 24; pl++) {
      transform33f3f(I->fracToReal(), ucv[unitCellLineIndices[pl]], v);
      copy3f(v, &vertexVals[pl * 3]);
    }

    CGOEnable(cgo, GL_LIGHTING);
    CGOStop(cgo);
  }
  return (cgo);
}
