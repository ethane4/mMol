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
-* sc
-*
-*
Z* -------------------------------------------------------------------
*/

// Helpers shared by the scene*.cc files; not part of the scene API.

#pragma once

#include "pymol/base/util.h"

/// Smallest distance between the front and back clipping planes.
#define cSliceMin 1.0F

/// Releases the cached scene image (see scene.cc).
void ScenePurgeImage(PyMOLGlobals* G);

/// Recomputes the inverse of the camera rotation matrix.
void SceneUpdateInvMatrix(PyMOLGlobals* G);

/// Advances the camera rock or sweep animation by the last frame's time.
void SceneUpdateCameraRock(PyMOLGlobals* G, int dirty);
