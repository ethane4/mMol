#include "pymol/graphics/graphics_util.h"

#include "pymol/base/pymol_globals.h"
#include "pymol/base/feedback.h"

#include <stdlib.h>

#include <iostream>
// UTIL
// Prints a backtrace during runtime of the last ^ stack frames
void print_trace();
void print_trace() {
  printf("Use debugger with `b %s` to get a backtrace\n", __func__);
}

// Returns the size in bytes of the opengl type
size_t gl_sizeof(GLenum type) {
  size_t size = 1;
  switch (type) {
    case GL_UNSIGNED_BYTE:
    case GL_BYTE:
      size = 1;
      break;
    case GL_UNSIGNED_SHORT:
    case GL_SHORT:
      size = 2;
      break;
    case GL_UNSIGNED_INT:
    case GL_INT:
    case GL_FLOAT:
      size = 4;
      break;
    default:
      printf("Unsupported GL Type!");
      break;
  }
  return size;
}
