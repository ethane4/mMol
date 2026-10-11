#pragma once

#include "pymol/platform/os_gl.h"

#include <string_view>

struct PyMOLGlobals;

static inline bool CheckGLErrorOK(PyMOLGlobals*, std::string_view) {
  return true;
}

size_t gl_sizeof(GLenum type);
