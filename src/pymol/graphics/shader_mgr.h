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
#ifndef _SHADER_MGR_H
#define _SHADER_MGR_H

#include "pymol/platform/os_gl.h"
#include "pymol/base/pymol_globals.h"
#include "pymol/reps/rep.h"
#include "pymol/graphics/generic_buffer.h"
#include "pymol/scene/scene_def.h"
#include <map>
#include <set>
#include <string>
#include <unordered_map>

#define GET_FRAGDEPTH_SUPPORT() 1

#include <mutex>
#define LOCK_GUARD_MUTEX(name, var) std::lock_guard<std::mutex> name(var)

/* ============================================================================
 * CShaderMgr class -- simple ShaderMgr for PyMOL
 * ============================================================================*/
class CShaderMgr {
 public:
  CShaderMgr(PyMOLGlobals* G);
  ~CShaderMgr();

  GLfloat* GetLineWidthRange();

  template <typename T, typename... TArgs>
  T* newGPUBuffer(TArgs&&... args) {
    std::hash<GPUBuffer*> hash_fn;
    T* buffer = new T(std::forward<TArgs>(args)...);
    auto bufptr = dynamic_cast<GPUBuffer*>(buffer);
    const size_t hashid = hash_fn(bufptr);
    buffer->set_hash_id(hashid);
    _gpu_object_map[hashid] = buffer;
    return buffer;
  }

  template <typename T>
  T* getGPUBuffer(size_t hashid) {
    return dynamic_cast<T*>(getGPUBuffer(hashid));
  }

  GPUBuffer* getGPUBuffer(size_t hashid) {
    auto search = _gpu_object_map.find(hashid);
    if (search != _gpu_object_map.end())
      return search->second;
    else
      return nullptr;
  }

  void bindGPUBuffer(size_t hashid);
  void freeGPUBuffer(size_t handle);
  void freeGPUBuffers(std::vector<size_t>&& handles);
  void freeGPUBuffers(size_t* arr, size_t len);

  void AddVBOToFree(GLuint vboid);
  void FreeAllVBOs();

  inline static std::uint32_t OpenGLDefaultFramebufferID = 0;

 private:
  void freeAllGPUBuffers();

 public:
  PyMOLGlobals* G;

  // for deleting buffers in the correct thread
  std::vector<GLuint> vbos_to_free;
  std::mutex vbos_to_free_mutex;
  std::mutex gpu_objects_to_free_mutex;

 private:
  std::unordered_map<size_t, GPUBuffer*> _gpu_object_map;
  std::vector<size_t> _gpu_objects_to_free_vector;

 public:
  std::map<int, std::string> attribute_uids;
  std::map<const std::string, int> attribute_uids_by_name;
  int GetAttributeUID(const char* name);
  const char* GetAttributeName(int);

  GLfloat line_width_range[2];

  // for glDrawBuffer (e.g. GL_BACK, unless we're using QOpenGLWidget)
  GLFramebufferConfig defaultBackbuffer{};  // Default framebuffer configuration
  GLFramebufferConfig currentFBConfig{};    // Current framebuffer configuration
  GLFramebufferConfig topLevelConfig{};     // Final render target in draw loop
                                            // (currently only known by Ortho)

 public:
  std::size_t offscreen_ortho_rt{};
  std::size_t offscreen_sized_image_rt{};

  /**
   * @brief Binds an offscreen render target for UI overlay rendering
   * @param extent extent of the offscreen render target
   * @param clear whether to clear the offscreen render target
   * @return framebuffer configuration of bound ortho offscreen render target
   */
  GLFramebufferConfig bindOffscreenOrtho(const Extent2D& extent,
                                         bool clear = true);

  /**
   * @brief Binds an offscreen render target for sized rendering
   * @param extent extent of the offscreen render target
   * @param clear whether to clear the offscreen render target
   * @return framebuffer configuration of bound offscreen render target
   */
  GLFramebufferConfig bindOffscreenSizedImage(const Extent2D& extent,
                                              bool clear = true);

  /**
   * @brief Reads pixel data from framebuffer
   * @param rect the rectangle to read from
   * @param srcConfig the framebuffer configuration to read from
   * @return pixel data in RGBA format
   */
  std::vector<unsigned char> readPixelsFrom(
      PyMOLGlobals* G, const Rect2D& rect,
      const GLFramebufferConfig& srcConfig);

  /**
   * @brief Draws pixel data to framebuffer
   * @param rect the rectangle to draw to
   * @param srcPixels the pixel data (RGBA) to draw
   * @param dstConfig the framebuffer configuration to draw to
   */
  void drawPixelsTo(PyMOLGlobals* G, const Rect2D& rect,
                    const std::byte* srcPixels,
                    const GLFramebufferConfig& dstConfig);

  /**
   * Sets the current draw buffer for rendering
   * @param drawBuffer Buffer to draw to (e.g. GL_BACK)
   * @note: This currently assumes the default framebuffer. Prefer
   * GLFramebufferConfig coverload.
   */
  void setDrawBuffer(GLenum drawBuffer);

  /**
   * Sets the current framebuffer or draw buffer for rendering
   * @param drawBuffer Framebuffer to draw to
   */
  void setDrawBuffer(GLFramebufferConfig config);
};

bool ShaderMgrInit(PyMOLGlobals* G);

#define VERTEX_POS_SIZE 3
#define VERTEX_COLOR_SIZE 4

#define VERTEX_POS 0
#define VERTEX_NORMAL 1
#define VERTEX_COLOR 2

#define CYLINDER_VERTEX1 0
#define CYLINDER_VERTEX2 1
#define CYLINDER_COLOR 2
#define CYLINDER_COLOR2 3
#define CYLINDER_RADIUS 4
#define CYLINDER_CAP 5

#define RAMP_OFFSETPT 0

#endif
