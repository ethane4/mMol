/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright Schrodinger, LLC.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

#include "pymol/platform/os_gl.h"
#include "pymol/platform/os_python.h"
#include <string.h>
#include <iostream>
#include "pymol/graphics/shader_mgr.h"
#include "pymol/base/list_macros.h"
#include "pymol/base/pymol_options.h"
#include "pymol/base/feedback.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/settings/setting.h"
#include "pymol/scene/scene.h"
#include "pymol/color/color.h"
#include "pymol/geometry/vector.h"
#include "pymol/base/util.h"
#include "pymol/base/util2.h"
#include "pymol/graphics/texture.h"
#include "pymol/base/file_stream.h"
#include "pymol/geometry/matrix.h"
#include "pymol/base/parse.h"
#include "pymol/graphics/renderer.h"

#include "pymol/graphics/cgo.h"

#include <algorithm>
#include <sstream>
#include <stack>
#include <vector>
#include <functional>
/* ============================================================================
 * ShaderMgrInit is called from PyMOL.c during start up; it just allocates
 * the global ShaderMgr
 */
bool ShaderMgrInit(PyMOLGlobals* G) {
  G->ShaderMgr = new CShaderMgr(G);

  if (!G->ShaderMgr) return false;

  return true;
}

/* ============================================================================
 * CShaderMgr -- Simple Shader Manager class
 */
CShaderMgr::CShaderMgr(PyMOLGlobals* G_) {
  G = G_;
  vbos_to_free.reserve(256);
}

CShaderMgr::~CShaderMgr() {
  freeGPUBuffer(offscreen_ortho_rt);
  freeGPUBuffer(offscreen_sized_image_rt);

  FreeAllVBOs();
}

/**
 * glDeleteBuffers for vbos_to_free
 */
void CShaderMgr::FreeAllVBOs() {
  freeAllGPUBuffers();

  LOCK_GUARD_MUTEX(lock, vbos_to_free_mutex);

  if (vbos_to_free.empty()) return;

  glDeleteBuffers(vbos_to_free.size(), &vbos_to_free[0]);

  vbos_to_free.clear();
}

/**
 * thread-safe deferred glDeleteBuffers(1, &vboid)
 */
void CShaderMgr::AddVBOToFree(GLuint vboid) {
  LOCK_GUARD_MUTEX(lock, vbos_to_free_mutex);
  vbos_to_free.push_back(vboid);
}

/*
 * Lighting setting indices are not contiguous, so we need a mapping array
 */
int light_setting_indices[] = {
    cSetting_light,  cSetting_light2, cSetting_light3,
    cSetting_light4, cSetting_light5, cSetting_light6,
    cSetting_light7, cSetting_light8, cSetting_light9};

GLfloat* CShaderMgr::GetLineWidthRange() { return line_width_range; }

void CShaderMgr::bindGPUBuffer(size_t hashid) {
  auto search = _gpu_object_map.find(hashid);
  if (search != _gpu_object_map.end()) search->second->bind();
}

void CShaderMgr::freeGPUBuffer(size_t hashid) {
  if (!hashid) return;
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  _gpu_objects_to_free_vector.push_back(hashid);
}

void CShaderMgr::freeGPUBuffers(std::vector<size_t>&& hashids) {
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  _gpu_objects_to_free_vector.insert(_gpu_objects_to_free_vector.end(),
                                     hashids.begin(), hashids.end());
}

void CShaderMgr::freeGPUBuffers(size_t* arr, size_t len) {
  for (unsigned int i = 0; i < len; ++i) freeGPUBuffer(arr[i]);
}

void CShaderMgr::freeAllGPUBuffers() {
  LOCK_GUARD_MUTEX(lock, gpu_objects_to_free_mutex);
  for (auto hashid : _gpu_objects_to_free_vector) {
    auto search = _gpu_object_map.find(hashid);
    if (search != _gpu_object_map.end()) {
      if (search->second) {
        // Non-GL renderers cache their GPU upload keyed on the address of this
        // buffer's CPU copy. Drop that entry before the address is freed —
        // the replacement buffer for a rebuilt rep very often lands on the
        // same block, and would otherwise be drawn with the old contents.
        if (G->Renderer) {
          if (const void* cpu = search->second->cpuDataPtr())
            G->Renderer->invalidateVBOCacheEntry(cpu);
        }
        delete search->second;
      }
      _gpu_object_map.erase(search);
    }
  }
  _gpu_objects_to_free_vector.clear();
}

int CShaderMgr::GetAttributeUID(const char* name) {
  auto uloc = attribute_uids_by_name.find(name);
  if (uloc != attribute_uids_by_name.end()) return uloc->second;

  int uid = attribute_uids_by_name.size() + 1;
  attribute_uids_by_name[name] = uid;
  attribute_uids[uid] = name;
  return uid;
}
const char* CShaderMgr::GetAttributeName(int uid) {
  auto uloc = attribute_uids.find(uid);
  if (uloc == attribute_uids.end()) return nullptr;
  return attribute_uids[uid].c_str();
}

GLFramebufferConfig CShaderMgr::bindOffscreenOrtho(const Extent2D& extent,
                                                   bool clear) {
  using namespace tex;
  RenderTargetGL::shape_type req_size(extent.width, extent.height);
  if (!offscreen_ortho_rt) {
    auto rt = newGPUBuffer<RenderTargetGL>(req_size);
    rt->layout({{4, rt_layout_t::UBYTE}});
    offscreen_ortho_rt = rt->get_hash_id();
  }

  auto rt = getGPUBuffer<RenderTargetGL>(offscreen_ortho_rt);
  if (rt->size() != req_size) {
    freeGPUBuffer(offscreen_ortho_rt);
    rt = newGPUBuffer<RenderTargetGL>(req_size);
    rt->layout({{4, rt_layout_t::UBYTE}});
    offscreen_ortho_rt = rt->get_hash_id();
  }
  rt->bind(clear);
  return GLFramebufferConfig{
      static_cast<std::uint32_t>(offscreen_ortho_rt),  //
      GL_COLOR_ATTACHMENT0                             //
  };
}

GLFramebufferConfig CShaderMgr::bindOffscreenSizedImage(const Extent2D& extent,
                                                        bool clear) {
  using namespace tex;
  RenderTargetGL::shape_type req_size(extent.width, extent.height);
  if (!offscreen_sized_image_rt) {
    auto rt = newGPUBuffer<RenderTargetGL>(req_size);
    rt->layout({{4, rt_layout_t::UBYTE}});
    offscreen_sized_image_rt = rt->get_hash_id();
  }
  auto rt = getGPUBuffer<RenderTargetGL>(offscreen_sized_image_rt);
  if (rt->size() != req_size) {
    freeGPUBuffer(offscreen_sized_image_rt);
    rt = newGPUBuffer<RenderTargetGL>(req_size);
    rt->layout({{4, rt_layout_t::UBYTE}});
    offscreen_sized_image_rt = rt->get_hash_id();
  }
  rt->bind(clear);
  return GLFramebufferConfig{
      static_cast<std::uint32_t>(offscreen_sized_image_rt),  //
      GL_COLOR_ATTACHMENT0                                   //
  };
}

void CShaderMgr::setDrawBuffer(GLenum mode) {
  if (mode == GL_BACK) {
    mode = G->ShaderMgr->defaultBackbuffer.drawBuffer;
  }

  if (!hasFrameBufferBinding() &&
      (mode != G->ShaderMgr->currentFBConfig.drawBuffer) && G->HaveGUI &&
      G->ValidContext) {
    G->ShaderMgr->currentFBConfig.drawBuffer = mode;
  }
}

void CShaderMgr::setDrawBuffer(GLFramebufferConfig config) {
  if (config.drawBuffer == GL_BACK) {
    config = defaultBackbuffer;
  }
  if (glBindFramebuffer != nullptr &&
      config.framebuffer == defaultBackbuffer.framebuffer) {
    glBindFramebuffer(GL_FRAMEBUFFER, config.framebuffer);
    setDrawBuffer(config.drawBuffer);
    return;
  }

  if (auto rt = getGPUBuffer<RenderTargetGL>(config.framebuffer)) {
    rt->bind(false);
  }
}

std::vector<unsigned char> CShaderMgr::readPixelsFrom(
    PyMOLGlobals* G, const Rect2D& rect, const GLFramebufferConfig& srcConfig) {
  constexpr std::size_t pixelSizeBytes = 4;
  std::vector<unsigned char> dstPixels(
      rect.extent.width * rect.extent.height * pixelSizeBytes, 0);
  int prevReadFBO;
  int prevDrawFBO;
  int prevReadBuffer;
  glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prevReadFBO);
  glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &prevDrawFBO);
  glGetIntegerv(GL_READ_BUFFER, &prevReadBuffer);

  if (srcConfig.framebuffer == OpenGLDefaultFramebufferID) {
    glBindFramebuffer(GL_READ_FRAMEBUFFER, OpenGLDefaultFramebufferID);
    glReadBuffer(srcConfig.drawBuffer);
  } else {
    if (auto rt = getGPUBuffer<RenderTargetGL>(srcConfig.framebuffer)) {
      // TODO: bindOnlyAsRead
      rt->fbo()->bind();
    }
  }
  PyMOLReadPixels(rect.offset.x, rect.offset.y, rect.extent.width,
                  rect.extent.height, GL_RGBA, GL_UNSIGNED_BYTE,
                  dstPixels.data());

  // Restore State
  glBindFramebuffer(GL_READ_FRAMEBUFFER, prevReadFBO);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, prevDrawFBO);
  glReadBuffer(prevReadBuffer);
  return dstPixels;
}

void CShaderMgr::drawPixelsTo(PyMOLGlobals* G, const Rect2D& rect,
                              const std::byte* srcPixels,
                              const GLFramebufferConfig& dstConfig) {
  int prevReadFBO;
  int prevDrawFBO;
  int prevDrawBuffer;
  glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prevReadFBO);
  glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &prevDrawFBO);
  glGetIntegerv(GL_READ_BUFFER, &prevDrawBuffer);

  if (dstConfig.framebuffer == OpenGLDefaultFramebufferID) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, dstConfig.framebuffer);
  } else {
    if (auto rt = getGPUBuffer<RenderTargetGL>(dstConfig.framebuffer)) {
      // TODO: bindOnlyAsDraw?
      rt->fbo()->bind();
    }
  }
  glDrawBuffer(dstConfig.drawBuffer);
  PyMOLDrawPixels(rect.extent.width, rect.extent.height, GL_RGBA,
                  GL_UNSIGNED_BYTE, srcPixels);
  // Restore State
  glBindFramebuffer(GL_READ_FRAMEBUFFER, prevReadFBO);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, prevDrawFBO);
  glReadBuffer(prevDrawBuffer);
}
