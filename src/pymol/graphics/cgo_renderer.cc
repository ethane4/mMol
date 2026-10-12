#include "pymol/graphics/cgo_renderer.h"

#include "pymol/base/base.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/base/pymol_globals.h"
#include "pymol/reps/rep.h"
#include "pymol/settings/setting.h"

struct RenderInfo;
struct Rep;
struct PyMOLGlobals;
struct CSetting;

bool CGORendererInit(PyMOLGlobals* G) {
  G->CGORenderer = new CCGORenderer();
  if (!G->CGORenderer) {
    return false;
  }
  G->CGORenderer->G = G;
  G->CGORenderer->alpha = 1.0F;
  return true;
}

void CGORendererFree(PyMOLGlobals* G) { DeleteP(G->CGORenderer); }
