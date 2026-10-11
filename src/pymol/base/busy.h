// Progress reporting from long operations.
//
// Low-level code (triangulation, ray tracing, surface building) reports
// progress here; the GUI that displays it registers a handler at startup, so
// that code doesn't depend on the GUI. The mMol app reads the two progress
// values from another thread to draw its progress bar.

#pragma once

struct PyMOLGlobals;

/* The running `ray`'s progress, 0..1, or < 0 when none runs. Written by the
   ray (and GPU completion handlers), read by the app's progress bar from
   another thread, hence volatile. */
extern "C" volatile float g_mmol_ray_progress;
inline void RayProgress(float f) {
  if (g_mmol_ray_progress >= 0 && f > g_mmol_ray_progress)
    g_mmol_ray_progress = f;
}

/* Progress of the rep being built (surface, mesh, dots), 0..1, for the same
   progress bar. BusyFast reports PyMOL's coarse fractions; a build that
   reports its own time-weighted steps (RepSurface) sets g_mmol_busy_explicit
   so those don't overwrite it. */
extern "C" volatile float g_mmol_busy_progress;
extern bool g_mmol_busy_explicit;
inline void BusyProgress(float f) { g_mmol_busy_progress = f; }

/// Handler that displays fine-grained progress (OrthoBusyFast).
using BusyFastHandler = void (*)(PyMOLGlobals* G, int progress, int total);

/// Registers the handler BusyFast forwards to.
void BusySetFastHandler(BusyFastHandler handler);

/// Reports `progress` of `total` steps to the registered handler, if any.
void BusyFast(PyMOLGlobals* G, int progress, int total);
