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
#define MAX_LOG_LEN 1024

#include <algorithm>
#include <sstream>
#include <stack>
#include <vector>
#include <functional>
/*  Texture Usage:

    0 - for not-PURE_OPENGL_ES_2: ObjectVolume: volumeTex
        for WEBGL: SCHRODINGER logo
    1 - for not-PURE_OPENGL_ES_2: Volume: either colorTex1D or colorTex2D
        for _PYMOL_PRECOMPUTED_LIGHTING: Lighting Texture
   (ShaderMgr->lightingTexture) 2 - FXAA - color_texture SMAA1 - colorTex SMAA3
   - colorTex 3 - SMAA3 - blendTex Label/Indicator Shader : textureMap (both
   PURE_OPENGL_ES_2 and non-PURE_OPENGL_ES_2 4 - Background Texture:
   bgTextureMap 5 - OIT - 2nd pass : accumTex Volumes - carvemask 6 - OIT - 2nd
   pass : revealageTex SMAA2 - edgesTex 7 - SMAA2 - areaTex OIT Copy - colorTex
    8 - SMAA2 - searchTex

 */

#define MASK_SHADERS_PRESENT_GEOMETRY 0x2;

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
  current_shader = nullptr;
  shaders_present = 0;
  stereo_flag = 0;
  stereo_blend = 0;
  print_warnings = 1;

  lightingTexture = 0;
  is_picking = 0;

  vbos_to_free.reserve(256);
}

CShaderMgr::~CShaderMgr() {
  for (auto& prog : programs) {
    delete prog.second;
  }
  programs.clear();

  freeGPUBuffer(offscreen_rt);
  freeGPUBuffer(offscreen_ortho_rt);
  freeGPUBuffer(offscreen_sized_image_rt);

  FreeAllVBOs();
}

int CShaderMgr::RemoveShaderPrg(const std::string& name) {
  if (programs.find(name) != programs.end()) {
    delete programs[name];
  }
  return 1;
}

/**
 * Lookup a shader program by name and set it as the `current_shader` of the
 * shader manager. If `pass` is provided and is less than zero, and we are
 * in transparency_mode 3, then look up the NO_ORDER_TRANSP derivative.h
 */
CShaderPrg* CShaderMgr::GetShaderPrg(std::string name, short set_current_shader,
                                     RenderPass pass) {
  if (pass == RenderPass::Transparent &&
      SettingGetGlobal_i(G, cSetting_transparency_mode) == 3) {
    name += "_t";
  }

  auto it = programs.find(name);
  if (it == programs.end()) return nullptr;

  if (set_current_shader) current_shader = it->second;

  return it->second;
}

int CShaderMgr::ShaderPrgExists(const char* name) {
  return (programs.find(name) != programs.end());
}

int CShaderMgr::ShadersPresent() { return shaders_present; }

int CShaderMgr::GeometryShadersPresent() {
  return shaders_present & MASK_SHADERS_PRESENT_GEOMETRY;
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

CShaderPrg* CShaderMgr::Enable_DefaultShaderWithSettings(const CSetting* set1,
                                                         const CSetting* set2,
                                                         RenderPass pass) {
  CShaderPrg* shaderPrg = Get_DefaultShader(pass);
  return Setup_DefaultShader(shaderPrg, set1, set2);
}

CShaderPrg* CShaderMgr::Enable_DefaultShader(RenderPass pass) {
  CShaderPrg* shaderPrg = Get_DefaultShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg* CShaderMgr::Enable_LineShader(RenderPass pass) {
  CShaderPrg* shaderPrg = Get_LineShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg* CShaderMgr::Enable_SurfaceShader(RenderPass pass) {
  CShaderPrg* shaderPrg = Get_SurfaceShader(pass);
  return Setup_DefaultShader(shaderPrg, nullptr, nullptr);
}

CShaderPrg* CShaderMgr::Enable_ConnectorShader(RenderPass pass) {
  CShaderPrg* shaderPrg = Get_ConnectorShader(pass);
  if (!shaderPrg) return nullptr;
  shaderPrg = Setup_DefaultShader(shaderPrg, nullptr, nullptr);
  shaderPrg->SetLightingEnabled(0);
  {
    float front, back;
    front = SceneGetCurrentFrontSafe(G);
    back = SceneGetCurrentBackSafe(G);
    shaderPrg->Set1f("front", front);
    shaderPrg->Set1f("clipRange", back - front);
  }

  auto extent = SceneGetExtentStereo(G);
  shaderPrg->Set2f("screenSize", extent.width, extent.height);

  {
    float v_scale = SceneGetScreenVertexScale(G, nullptr);
    shaderPrg->Set1f("screenOriginVertexScale", v_scale / 2.f);
  }

  return shaderPrg;
}

CShaderPrg* CShaderMgr::Setup_DefaultShader(CShaderPrg* shaderPrg,
                                            const CSetting* set1,
                                            const CSetting* set2) {
  if (!shaderPrg) {
    current_shader = nullptr;
    return shaderPrg;
  }
  shaderPrg->Enable();
  shaderPrg->SetBgUniforms();
  shaderPrg->Set_Stereo_And_AnaglyphMode();

  bool two_sided_lighting_enabled =
      SceneGetTwoSidedLightingSettings(G, set1, set2);

  shaderPrg->SetLightingEnabled(1);  // lighting on by default
  shaderPrg->Set1i("two_sided_lighting_enabled", two_sided_lighting_enabled);
  shaderPrg->Set1f("ambient_occlusion_scale", 0.f);
  shaderPrg->Set1i("accessibility_mode",
                   SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) / 4);
  shaderPrg->Set1f(
      "accessibility_mode_on",
      SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) ? 1.f : 0.f);

  // interior_color
  {
    int interior_color =
        SettingGet_i(G, set1, set2, cSetting_ray_interior_color);
    if (interior_color == cColorDefault || two_sided_lighting_enabled) {
      shaderPrg->Set1i("use_interior_color", 0);
    } else {
      float inter[] = {0.f, 0.f, 0.f};
      ColorGetEncoded(G, interior_color, inter);
      shaderPrg->Set1i("use_interior_color", 1);
      shaderPrg->Set4f("interior_color", inter[0], inter[1], inter[2], 1.f);
    }
  }

  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();
  return (shaderPrg);
}
CShaderPrg* CShaderMgr::Enable_CylinderShader(RenderPass pass) {
  return Enable_CylinderShader("cylinder", pass);
}

CShaderPrg* CShaderMgr::Enable_CylinderShader(const char* shader_name,
                                              RenderPass pass) {
  CShaderPrg* shaderPrg;

  auto extent = SceneGetExtentStereo(G);
  shaderPrg = GetShaderPrg(shader_name, 1, pass);
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();

  shaderPrg->SetLightingEnabled(1);  // lighting on by default

  shaderPrg->Set1f("uni_radius", 0.f);

  shaderPrg->Set_Stereo_And_AnaglyphMode();

  shaderPrg->Set1f("inv_height", 1.0f / extent.height);
  shaderPrg->Set1i("no_flat_caps", 1);
  {
    float smooth_half_bonds =
        (SettingGetGlobal_i(G, cSetting_smooth_half_bonds)) ? .2f : 0.f;
    shaderPrg->Set1f("half_bond", smooth_half_bonds);
  }
  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();

  shaderPrg->SetBgUniforms();

  // always enable backface culling for cylinders
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  return shaderPrg;
}

CShaderPrg* CShaderMgr::Get_DefaultSphereShader(RenderPass pass) {
  return GetShaderPrg("sphere", 1, pass);
}

CShaderPrg* CShaderMgr::Enable_DefaultSphereShader(RenderPass pass) {
  CShaderPrg* shaderPrg = Get_DefaultSphereShader(pass);
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();
  shaderPrg->SetLightingEnabled(1);
  shaderPrg->Set1f("sphere_size_scale", 1.f);

  shaderPrg->Set_Stereo_And_AnaglyphMode();

  shaderPrg->Set_Specular_Values();
  shaderPrg->Set_Matrices();

  shaderPrg->SetBgUniforms();

  return (shaderPrg);
}

CShaderPrg* CShaderMgr::Get_ConnectorShader(RenderPass pass) {
  return GetShaderPrg("connector", 1, pass);
}

CShaderPrg* CShaderMgr::Get_DefaultShader(RenderPass pass) {
  return GetShaderPrg("default", 1, pass);
}

CShaderPrg* CShaderMgr::Get_LineShader(RenderPass pass) {
  return GetShaderPrg("line", 1, pass);
}

CShaderPrg* CShaderMgr::Get_SurfaceShader(RenderPass pass) {
  return GetShaderPrg("surface", 1, pass);
}

CShaderPrg* CShaderMgr::Get_CylinderShader(RenderPass pass,
                                           short set_current_shader) {
  return GetShaderPrg("cylinder", set_current_shader, pass);
}

CShaderPrg* CShaderMgr::Get_Current_Shader() { return current_shader; }

CShaderPrg* CShaderMgr::Get_BackgroundShader() { return GetShaderPrg("bg"); }

CShaderPrg* CShaderMgr::Enable_BackgroundShader() {
  CShaderPrg* shaderPrg = Get_BackgroundShader();
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  glDisable(GL_DEPTH_TEST);
  shaderPrg->SetBgUniforms();

  return shaderPrg;
}

CShaderPrg* CShaderMgr::Enable_TriLinesShader() {
  CShaderPrg* shaderPrg = GetShaderPrg("trilines");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();
  shaderPrg->SetBgUniforms();
  shaderPrg->Set_Stereo_And_AnaglyphMode();
  shaderPrg->Set_Matrices();
  auto extent = SceneGetExtentStereo(G);
  shaderPrg->Set2f("inv_dimensions", 1.f / extent.width, 1.f / extent.height);
  return shaderPrg;
}

CShaderPrg* CShaderMgr::Enable_OITShader() {
  CShaderPrg* shaderPrg = GetShaderPrg("oit");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  constexpr GLuint accumTexUnit = 5;
  constexpr GLuint revealageTexUnit = 6;
  oit_pp->activateRTAsTexture(OIT_PostProcess::OITRT::ACCUM, accumTexUnit);
  oit_pp->activateRTAsTexture(OIT_PostProcess::OITRT::REVEALAGE,
                              revealageTexUnit);
  shaderPrg->Set1i("accumTex", accumTexUnit);
  shaderPrg->Set1i("revealageTex", revealageTexUnit);
  shaderPrg->Set1f("isRight", stereo_flag > 0 ? 1. : 0);
  glEnable(GL_BLEND);
  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE,
                      GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_DEPTH_TEST);
  return shaderPrg;
}

CShaderPrg* CShaderMgr::Enable_OITCopyShader() {
  CShaderPrg* shaderPrg = GetShaderPrg("copy");
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  constexpr GLuint colorTexUnit = 7;
  activateOffscreenTexture(colorTexUnit);
  shaderPrg->Set1i("colorTex", colorTexUnit);
  if (G->ShaderMgr->stereo_blend) {
    // for full-screen stereo
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
  } else {
    glDisable(GL_BLEND);
  }
  glDisable(GL_DEPTH_TEST);
  return shaderPrg;
}

CShaderPrg* CShaderMgr::Enable_LabelShader(RenderPass pass) {
  CShaderPrg* shaderPrg;
  shaderPrg = Get_LabelShader(pass);
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();
  return Setup_LabelShader(shaderPrg);
}

CShaderPrg* CShaderMgr::Enable_ScreenShader() {
  CShaderPrg* shaderPrg;
  shaderPrg = Get_ScreenShader();
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();

  int ortho_width, ortho_height;
  std::tie(ortho_width, ortho_height) = OrthoGetSize(*G->Ortho);
  shaderPrg->Set2f("t2PixelSize", 2.f / ortho_width, 2.f / ortho_height);

  return Setup_LabelShader(shaderPrg);
}

CShaderPrg* CShaderMgr::Enable_RampShader() {
  CShaderPrg* shaderPrg;
  shaderPrg = Get_RampShader();
  if (!shaderPrg) return nullptr;
  shaderPrg->Enable();
  return Setup_LabelShader(shaderPrg);
}

CShaderPrg* CShaderMgr::Setup_LabelShader(CShaderPrg* shaderPrg) {
  shaderPrg->Set_Matrices();

  glActiveTexture(GL_TEXTURE3);
  TextureBindTexture(G);
  if (!(shaderPrg->uniform_set & 8)) {
    shaderPrg->uniform_set |= 8;
    shaderPrg->Set1i("textureMap", 3);
  }

  auto extent = SceneGetExtentStereo(G);

  shaderPrg->Set2f("screenSize", extent.width, extent.height);

  shaderPrg->SetBgUniforms();

  {
    float v_scale = SceneGetScreenVertexScale(G, nullptr);
    shaderPrg->Set1f("screenOriginVertexScale", v_scale / 2.f);
  }
  {
    float front, back;
    front = SceneGetCurrentFrontSafe(G);
    back = SceneGetCurrentBackSafe(G);
    shaderPrg->Set1f("front", front);
    shaderPrg->Set1f("clipRange", back - front);
  }

  return shaderPrg;
}

CShaderPrg* CShaderMgr::Get_LabelShader(RenderPass pass) {
  return GetShaderPrg("label", 1, pass);
}

CShaderPrg* CShaderMgr::Get_ScreenShader() {
  if (is_picking) return nullptr;
  return GetShaderPrg("screen");
}

CShaderPrg* CShaderMgr::Get_RampShader() { return GetShaderPrg("ramp"); }

CShaderPrg* CShaderMgr::Get_IndicatorShader() {
  return GetShaderPrg("indicator");
}

CShaderPrg* CShaderMgr::Enable_IndicatorShader() {
  CShaderPrg* shaderPrg = Get_IndicatorShader();
  if (!shaderPrg) return shaderPrg;
  shaderPrg->Enable();

  shaderPrg->Set_Stereo_And_AnaglyphMode();
  shaderPrg->Set_Matrices();

  glActiveTexture(GL_TEXTURE3);
  TextureBindTexture(G);
  if (!(shaderPrg->uniform_set & 8)) {
    shaderPrg->Set1i("textureMap", 3);
    shaderPrg->uniform_set |= 8;
  }
  shaderPrg->SetMat4fc("g_ModelViewMatrix", SceneGetModelViewMatrixPtr(G));
  shaderPrg->SetMat4fc("g_ProjectionMatrix", SceneGetProjectionMatrixPtr(G));

  return (shaderPrg);
}

CShaderPrg* CShaderMgr::Enable_BezierShader() {
  auto shaderPrg = Get_BezierShader();
  if (!shaderPrg) {
    return nullptr;
  }
  shaderPrg->Enable();

  glPatchParameteri(GL_PATCH_VERTICES, 4);

  const float segmentCount = 256.0f;
  const float stripCount = 1.0f;

  shaderPrg->Set1f("segmentCount", segmentCount);
  shaderPrg->Set1f("stripCount", stripCount);

  shaderPrg->SetMat4fc("g_ModelViewMatrix", SceneGetModelViewMatrixPtr(G));
  shaderPrg->SetMat4fc("g_ProjectionMatrix", SceneGetProjectionMatrixPtr(G));
  return shaderPrg;
}

CShaderPrg* CShaderMgr::Get_BezierShader() { return GetShaderPrg("bezier"); }

void CShaderMgr::ResetUniformSet() {
  for (auto& prog : programs) {
    prog.second->uniform_set = 0;
  }
}

void CShaderMgr::SetIsPicking(int is_picking) { this->is_picking = is_picking; }

#define LIGHTINGTEXTUREWIDTH 64

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

// SceneRenderBindToOffscreen
void CShaderMgr::bindOffscreen(int width, int height, GridInfo* grid) {
  using namespace tex;
  RenderTargetGL::shape_type req_size(width, height);
  RenderTargetGL* rt = nullptr;

  // Doesn't exist, create
  if (!offscreen_rt) {
    CGOFree(G->Scene->offscreenCGO);
    rt = newGPUBuffer<RenderTargetGL>(req_size);
    rt->layout({{4, rt_layout_t::UBYTE}});
    offscreen_rt = rt->get_hash_id();
  } else {
    rt = getGPUBuffer<RenderTargetGL>(offscreen_rt);

    // resize
    if (req_size != rt->size()) {
      rt->resize(req_size);
    }
  }

  if (rt) rt->bind(!stereo_blend);
  glEnable(GL_BLEND);

  SceneInitializeViewport(G, true);
  if (grid->active) {
    grid->cur_view.offset = Offset2D{};
    grid->cur_view.extent.width = static_cast<std::uint32_t>(req_size.x);
    grid->cur_view.extent.height = static_cast<std::uint32_t>(req_size.y);
  }
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

void CShaderMgr::activateOffscreenTexture(GLuint textureIdx) {
  glActiveTexture(GL_TEXTURE0 + textureIdx);
  auto t = getGPUBuffer<RenderTargetGL>(offscreen_rt);
  if (t->_textures[0]) t->_textures[0]->bind();
}

void CShaderMgr::Disable_Current_Shader() {
  if (current_shader) {
    current_shader->Disable();
  }
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
