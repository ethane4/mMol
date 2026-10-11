#include "pymol/base/busy.h"

extern "C" volatile float g_mmol_ray_progress = -1.0f;
extern "C" volatile float g_mmol_busy_progress = 0.0f;
bool g_mmol_busy_explicit = false;

static BusyFastHandler s_fast_handler = nullptr;

void BusySetFastHandler(BusyFastHandler handler) { s_fast_handler = handler; }

void BusyFast(PyMOLGlobals* G, int progress, int total) {
  if (s_fast_handler) s_fast_handler(G, progress, total);
}
