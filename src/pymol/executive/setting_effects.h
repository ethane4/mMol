// Side effects of changing a setting: invalidating reps, redrawing,
// updating the scene and the GUI. Kept out of the settings store, which the
// low-level components use and which must not depend on them.

#pragma once

#include <vector>

struct PyMOLGlobals;

/// Settings changed since the last call, for `name` (empty: global) and
/// `state`, as setting indices.
std::vector<int> SettingGetUpdateList(PyMOLGlobals* G, const char* name = "",
                                      int state = 0);

/// Sets the front (contrast) color from the background settings.
void ColorUpdateFrontFromSettings(PyMOLGlobals* G);

/// Applies the side effects of setting `index` changing on `sele` in `state`.
void SettingGenerateSideEffects(PyMOLGlobals* G, int index, const char* sele,
                                int state, int quiet);
