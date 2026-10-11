#pragma once

#include "pymol/platform/os_gl.h"

#define GLEW_EXT_gpu_shader4 false
#define GLEW_EXT_geometry_shader4 false
#define GLEW_ARB_tessellation_shader false
#define GLEW_ARB_gpu_shader5 false
#define GLEW_EXT_draw_buffers2 false

#define TM3_IS_ONEBUF !GLEW_EXT_draw_buffers2
