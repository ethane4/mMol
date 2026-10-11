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

#include "pymol/graphics/immediate_helper.h"

#include "pymol/geometry/matrix.h"
#include "pymol/scene/scene.h"
#include "pymol/ray/scene_ray.h"
#include "pymol/scene/scene_mouse.h"
#include "pymol/scene/scene_picking.h"
#include "pymol/ui/control.h"
#include "pymol/scene/movie.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/ui/wizard.h"
#include "pymol/ui/seq.h"
#include "pymol/text/text.h"
#include "pymol/base/pymol_options.h"
#include "pymol/base/feedback.h"
#include "pymol/graphics/gfx_manager.h"

// #define _OPENVR_STEREO_DEBUG_VIEWS

#include "pymol/scene/scene_internal.h"

static void glReadBufferError(PyMOLGlobals* G, GLenum b, GLenum e) {
  PRINTFB(G, FB_OpenGL, FB_Warnings)
  " WARNING: glReadBuffer caused GL error 0x%04x\n", e ENDFB(G);
}
// TH 2013-11-01: glReadBuffer fails in JyMOL when picking, OSX 10.9, Intel
// Graphics for minor cases (i.e., png is called) this might get called outside
// of the main thread (from ExecutiveDrawNow()) in this case, just don't call
// glReadBuffer for now, it should be ok because i believe it is used to figure
// out the size of the 3D window (using SceneImagePrepareImpl) in situations
// where the size gets changed.
#define glReadBuffer(b)                                 \
  {                                                     \
    int e;                                              \
    if (PIsGlutThread()) glReadBuffer(b);               \
    if ((e = glGetError())) glReadBufferError(G, b, e); \
  }

#define SceneLineHeight 127
#define SceneTopMargin 0
#define SceneBottomMargin 3
#define SceneLeftMargin 3

/* Shared with ShaderMgr */
/** Coefficients from: http://3dtv.at/Knowhow/AnaglyphComparison_en.aspx */
/** Optimize the look and feel of anaglyph 3D */
/* the last mode is the 3x3 identity */
// matrices are column major
float anaglyphL_constants[6][9] = {
    {0.299, 0.000, 0.000, 0.587, 0.000, 0.000, 0.114, 0.000, 0.000},
    {0.299, 0.000, 0.000, 0.587, 0.000, 0.000, 0.114, 0.000, 0.000},
    {1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
    {0.299, 0.000, 0.000, 0.587, 0.000, 0.000, 0.114, 0.000, 0.000},
    {0.000, 0.000, 0.000, 0.700, 0.000, 0.000, 0.300, 0.000, 0.000},
    {1.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000}};

float anaglyphR_constants[6][9] = {
    {0.000, 0.000, 0.299, 0.000, 0.000, 0.587, 0.000, 0.000, 0.114},
    {0.000, 0.299, 0.299, 0.000, 0.587, 0.587, 0.000, 0.114, 0.114},
    {0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000},
    {0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000},
    {0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000},
    {1.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000}};

#define F2UI(a) (unsigned int)((a) * 255.0)

/* allow up to 10 seconds at 30 FPS */

/* EXPERIMENTAL VOLUME RAYTRACING DATA */
extern float* rayDepthPixels;
extern int rayVolume, rayWidth, rayHeight;

static void SceneRestartPerfTimer(PyMOLGlobals* G);

int SceneViewEqual(SceneViewType left, SceneViewType right) {
  int i;
  for (i = 0; i < cSceneViewSize; i++) {
    if (fabs(left[i] - right[i]) > R_SMALL4) return false;
  }
  return true;
}

Rect2D GridSetRayViewport(GridInfo& I, int slot) {
  Rect2D view{};
  if (slot)
    I.slot = slot + I.first_slot - 1;
  else
    I.slot = slot;
  /* if we are in grid mode, then prepare the grid slot viewport */
  if (slot < 0) {
    return I.cur_view;
  } else if (slot == 0) {
    view.offset = Offset2D{};
    view.extent.width =
        static_cast<std::uint32_t>(I.cur_view.extent.width / I.n_col);
    view.extent.height =
        static_cast<std::uint32_t>(I.cur_view.extent.height / I.n_row);
    if (I.n_col < I.n_row) {
      view.extent.width *= I.n_col;
      view.extent.height *= I.n_col;
    } else {
      view.extent.width *= I.n_row;
      view.extent.height *= I.n_row;
    }
    view.offset.x +=
        I.cur_view.offset.x + (I.cur_view.extent.width - view.extent.width) / 2;
    view.offset.y += I.cur_view.offset.y;
  } else {
    int abs_grid_slot = slot - I.first_slot;
    int grid_col = abs_grid_slot % I.n_col;
    int grid_row = (abs_grid_slot / I.n_col);
    view.offset.x = static_cast<std::int32_t>(
        (grid_col * I.cur_view.extent.width) / I.n_col);
    view.offset.y = static_cast<std::int32_t>(
        I.cur_view.extent.height -
        ((grid_row + 1) * I.cur_view.extent.height) / I.n_row);
    view.extent.width =
        ((grid_col + 1) * I.cur_view.extent.width) / I.n_col - view.offset.x;
    view.extent.height = static_cast<std::uint32_t>(
        (I.cur_view.extent.height -
         ((grid_row)*I.cur_view.extent.height) / I.n_row) -
        view.offset.y);
    view.offset.x += I.cur_view.offset.x;
    view.offset.y += I.cur_view.offset.y;
  }
  return view;
}

void GridUpdate(GridInfo* I, float asp_ratio, GridMode mode, int size) {
  if (mode != GridMode::NoGrid) {
    I->size = size;
    I->mode = mode;
    {
      int n_row = 1;
      int n_col = 1;
      int r_size = size;
      while ((n_row * n_col) < r_size) {
        float asp1 = asp_ratio * (n_row + 1.0) / n_col;
        float asp2 = asp_ratio * (n_row) / (n_col + 1.0);
        if (asp1 < 1.0F) asp1 = 1.0 / asp1;
        if (asp2 < 1.0F) asp2 = 1.0 / asp2;
        if (fabs(asp1) > fabs(asp2))
          n_col++;
        else
          n_row++;
      }

      // the above algorithm can generate a 3x2 grid for size=4, but we want a
      // 2x2 in that case.
      while ((n_col - 1) * n_row >= size && size) {
        n_col -= 1;
      }
      while ((n_row - 1) * n_col >= size && size) {
        n_row -= 1;
      }

      I->n_row = n_row;
      I->n_col = n_col;
    }
    if (I->size > 1) {
      I->active = true;
      I->asp_adjust = (float)I->n_row / I->n_col;
      I->first_slot = 1;
      I->last_slot = I->size;
    } else {
      I->active = false;
    }
  } else {
    I->active = false;
  }
}

void SceneInvalidateStencil(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->StencilValid = false;
}

int SceneGetGridSize(PyMOLGlobals* G, GridMode grid_mode) {
  CScene* I = G->Scene;
  int size = 0;

  switch (grid_mode) {
    case GridMode::ByObject:
      if (I->m_slots.empty()) {
        I->m_slots.push_back(0);
      } else {
        std::fill(I->m_slots.begin(), I->m_slots.end(), 0);
      }
      {
        int max_slot = 0;
        for (auto& obj : I->Obj) {
          if (auto slot = obj->grid_slot) {
            max_slot = std::max(slot, max_slot);
            if (slot > 0) {
              VecCheck(I->m_slots, slot);
              I->m_slots[slot] = 1;
            }
          }
        }
        for (int slot = 0; slot <= max_slot; slot++) {
          if (I->m_slots[slot]) I->m_slots[slot] = ++size;
        }
      }
      break;
    case GridMode::ByObjectStates:
    case GridMode::ByObjectByState:
      if (!I->m_slots.empty()) {
        I->m_slots.clear();
      }
      {
        int max_slot = 0;
        for (auto& obj : I->Obj) {
          auto slot = obj->getNFrame();
          if (grid_mode == GridMode::ByObjectByState) {
            obj->grid_slot = max_slot;  // slot offset for 1st state
            max_slot += slot;
          } else if (max_slot < slot) {
            max_slot = slot;
          }
        }
        size = max_slot;
      }
      break;
  }
  auto grid_max = SettingGet<int>(G, cSetting_grid_max);
  if (grid_max >= 0) size = std::min(size, grid_max);
  return size;
}

int SceneMustDrawBoth(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (G->StereoCapable &&
          ((I->StereoMode == 1) ||
           SettingGetGlobal_b(G, cSetting_stereo_double_pump_mono)));
}

static int SceneDeferClickWhen(Block* block, int button, int x, int y,
                               double when, int mod);

int stereo_via_adjacent_array(int stereo_mode) {
  switch (stereo_mode) {
    case cStereo_crosseye:
    case cStereo_walleye:
    case cStereo_sidebyside:
      return true;
  }
  return false;
}

int StereoIsAdjacent(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return stereo_via_adjacent_array(I->StereoMode);
}

void SceneAbortAnimation(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  if (I->cur_ani_elem < I->n_ani_elem) { /* allow user to override animation */
    I->cur_ani_elem = I->n_ani_elem;
  }
}

void ScenePrimeAnimation(PyMOLGlobals* G) {
  if (G->HaveGUI) {
    CScene* I = G->Scene;
    UtilZeroMem(I->ani_elem, sizeof(CViewElem));
    SceneToViewElem(G, I->ani_elem, nullptr);
    I->ani_elem[0].specification_level = 2;
    I->n_ani_elem = 0;
  }
}

static float SceneGetFPS(PyMOLGlobals* G) {
  float fps = SettingGetGlobal_f(G, cSetting_movie_fps);
  float minTime;
  if (fps <= 0.0F) {
    if (fps < 0.0)
      minTime = 0.0; /* negative fps means full speed */
    else             /* 0 fps means use movie_delay instead */
      minTime = SettingGetGlobal_f(G, cSetting_movie_delay) / 1000.0;
    if (minTime >= 0.0F)
      fps = 1.0F / minTime;
    else
      fps = 1000.0F;
  }
  return fps;
}

/**
 * Release `G->Scene->Image` and clear related flags (`CopyType`).
 *
 * Invalidates the Ortho CGO and requests a redisplay. It's not entirely clear
 * why this is done, a code comment says "need to invalidate since text could be
 * shown".
 */
void ScenePurgeImage(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->CopyType = false;
  I->Image = nullptr;

  // TODO does this belong here?
  if (true /* !noinvalid */)
    OrthoInvalidateDoDraw(
        G);  // right now, need to invalidate since text could be shown
}

void SceneInvalidateCopy(PyMOLGlobals* G, int free_buffer) {
  CScene* I = G->Scene;
  if (I) {
    if (free_buffer) {
      ScenePurgeImage(G);
    } else {
      I->Image = nullptr;
    }
    if (I->CopyType)
      OrthoInvalidateDoDraw(
          G);  // right now, need to invalidate since text could be shown
    I->CopyType = false;
  }
}

void SceneInvalidate(PyMOLGlobals* G) {
  SceneInvalidateCopy(G, false);
  SceneDirty(G);
  PyMOL_NeedRedisplay(G->PyMOL);
}

void SceneLoadAnimation(PyMOLGlobals* G, double duration, int hand) {
  if (G->HaveGUI) {
    double now;
    int target = (int)(duration * 30);
    CScene* I = G->Scene;
    if (target < 1) target = 1;
    if (target > MAX_ANI_ELEM) target = MAX_ANI_ELEM;
    UtilZeroMem(I->ani_elem + 1, sizeof(CViewElem) * target);
    SceneToViewElem(G, I->ani_elem + target, nullptr);
    I->ani_elem[target].specification_level = 2;
    now = UtilGetSeconds(G);
    I->ani_elem[0].timing_flag = true;
    I->ani_elem[0].timing = now + 0.01;
    I->ani_elem[target].timing_flag = true;
    I->ani_elem[target].timing = now + duration;
    ViewElemInterpolate(G, I->ani_elem, I->ani_elem + target, 2.0F, 1.0F, true,
                        0.0F, hand, 0.0F);
    SceneFromViewElem(G, I->ani_elem, true);
    I->cur_ani_elem = 0;
    I->n_ani_elem = target;
    I->AnimationStartTime = UtilGetSeconds(G);
    I->AnimationStartFlag = true;
    I->AnimationStartFrame = SceneGetFrame(G);
    I->AnimationLagTime = 0.0;
  }
}

/**
 * Set FrontSafe and BackSafe
 */
void UpdateFrontBackSafe(CScene* I) {
  float front = I->m_view.m_clip().m_front;
  float back = I->m_view.m_clip().m_back;

  // minimum slab
  if (back - front < cSliceMin) {
    float avg = (back + front) / 2.0;
    back = avg + cSliceMin / 2.0;
    front = avg - cSliceMin / 2.0;
  }

  // minimum front
  if (front < cSliceMin) {
    front = cSliceMin;

    // minimum slab
    if (back < front + cSliceMin) back = front + cSliceMin;
  }

  I->m_view.m_clipSafe().m_front = front;
  I->m_view.m_clipSafe().m_back = back;
}

#define SELE_MODE_MAX 7

static const char SelModeKW[][20] = {
    "", "byresi", "bychain", "bysegi", "byobject", "bymol", "bca.",
};

void SceneUpdateInvMatrix(PyMOLGlobals* G) {
  // TODO: Test glm::inverse(I->m_view.rotMatrix())
  CScene* I = G->Scene;
  const auto rm = glm::value_ptr(I->m_view.rotMatrix());
  float* im = I->InvMatrix;
  im[0] = rm[0];
  im[1] = rm[4];
  im[2] = rm[8];
  im[3] = 0.0F;
  im[4] = rm[1];
  im[5] = rm[5];
  im[6] = rm[9];
  im[7] = 0.0F;
  im[8] = rm[2];
  im[9] = rm[6];
  im[10] = rm[10];
  im[11] = 0.0F;
  im[12] = 0.0F;
  im[13] = 0.0F;
  im[14] = 0.0F;
  im[15] = 1.0F;
}

void SceneUpdateStereo(PyMOLGlobals* G) {
  SceneSetStereo(G, SettingGetGlobal_b(G, cSetting_stereo));
  PyMOL_NeedRedisplay(G->PyMOL);
}

/**
 * Get the selection operator for the `mouse_selection_mode` setting.
 * Example: `mouse_selection_mode=2` -> "bychain"
 */
const char* SceneGetSeleModeKeyword(PyMOLGlobals* G) {
  int sel_mode = SettingGetGlobal_i(G, cSetting_mouse_selection_mode);
  if ((sel_mode >= 0) && (sel_mode < SELE_MODE_MAX))
    return (char*)SelModeKW[sel_mode];
  return (char*)SelModeKW[0];
}

SceneUnitContext ScenePrepareUnitContext(const Extent2D& extent) {
  SceneUnitContext context{};
  float tw = 1.0F;
  float th = 1.0F;
  float aspRat = extent.height != 0
                     ? (extent.width / static_cast<float>(extent.height))
                     : 1.0f;

  if (aspRat > 1.0F) {
    tw = aspRat;
  } else {
    th = 1.0F / aspRat;
  }

  context.unit_left = (1.0F - tw) / 2;
  context.unit_right = (tw + 1.0F) / 2;
  context.unit_top = (1.0F - th) / 2;
  context.unit_bottom = (th + 1.0F) / 2;
  context.unit_front = -0.5F;
  context.unit_back = 0.5F;
  /*
     printf(
     "ScenePrepareUnitContext:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
     context->unit_left,
     context->unit_right,
     context->unit_top,
     context->unit_bottom,
     context->unit_front,
     context->unit_back);
   */
  return context;
}

/**
 * cmd.get_viewport()
 */
void SceneGetWidthHeight(PyMOLGlobals* G, int* width, int* height) {
  CScene* I = G->Scene;
  *width = I->Width;
  *height = I->Height;
}

Extent2D SceneGetExtent(PyMOLGlobals* G) {
  return Extent2D{static_cast<std::uint32_t>(G->Scene->Width),
                  static_cast<std::uint32_t>(G->Scene->Height)};
}

void SceneSetExtent(PyMOLGlobals* G, const Extent2D& extent) {
  G->Scene->Width = extent.width;
  G->Scene->Height = extent.height;
}

Rect2D SceneGetRect(PyMOLGlobals* G) {
  auto I = G->Scene;
  return Rect2D{Offset2D{I->rect.left, I->rect.bottom}, SceneGetExtent(G)};
}

float SceneGetAspectRatio(PyMOLGlobals* G) {
  auto extent = SceneGetExtent(G);
  return static_cast<float>(extent.width) / static_cast<float>(extent.height);
}

/**
 * Get the actual current (sub-)viewport size, considering grid mode and
 * side-by-side stereo
 */
Extent2D SceneGetExtentStereo(PyMOLGlobals* G) {
  auto I = G->Scene;

  if (I->grid.active) {
    // TODO: this considers "draw W, H" (PYMOL-2775)
    return I->grid.cur_viewport_size;
  }

  Extent2D extent{static_cast<std::uint32_t>(I->Width),
                  static_cast<std::uint32_t>(I->Height)};
  // TODO: this does NOT consider "draw W, H" (PYMOL-2775)
  if (stereo_via_adjacent_array(I->StereoMode)) extent.width /= 2.f;
  return extent;
}

int SceneGetStereo(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (I->StereoMode);
}

void SceneGetCardInfo(PyMOLGlobals* G, char** vendor, char** renderer,
                      char** version) {
  CScene* I = G->Scene;
  (*vendor) = I->vendor;
  (*renderer) = I->renderer;
  (*version) = I->version;
}

void SceneSuppressMovieFrame(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->MovieFrameFlag = false;
}

/**
 * Get center of screen in world coordinates
 *
 * @param[out] pos 3f output vector
 */
void SceneGetCenter(PyMOLGlobals* G, float* pos) {
  CScene* I = G->Scene;
  auto camPos = I->m_view.pos();
  const auto& ori = I->m_view.origin();

  MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()),
                             glm::value_ptr(ori), pos);

  pos[0] -= camPos.x;
  pos[1] -= camPos.y;

  MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), pos,
                                pos);
}
int SceneGetNFrame(PyMOLGlobals* G, int* has_movie) {
  CScene* I = G->Scene;
  if (has_movie) *has_movie = I->HasMovie;
  return (I->NFrame);
}
void SceneDontCopyNext(PyMOLGlobals* G)

/* disables automatic copying of the image for the next rendering run */
{
  CScene* I = G->Scene;
  I->CopyNextFlag = false;
}
void SceneUpdateStereoMode(PyMOLGlobals* G) {
  if (G->Scene->StereoMode) {
    SceneSetStereo(G, true);
  }
}
void SceneSetStereo(PyMOLGlobals* G, bool flag) {
  CScene* I = G->Scene;
  int cur_stereo_mode = I->StereoMode;

  if (flag) {
    I->StereoMode = SettingGetGlobal_i(G, cSetting_stereo_mode);
  } else {
    I->StereoMode = 0;
  }

  SettingSetGlobal_b(G, cSetting_stereo, flag);

  if (cur_stereo_mode != I->StereoMode) {
    if (cur_stereo_mode == cStereo_geowall ||
        I->StereoMode == cStereo_geowall) {
      OrthoReshape(G, G->Option->winX, G->Option->winY, true);
    }

    SceneInvalidateStencil(G);
    SceneInvalidate(G);
  }
}
void SceneSetMatrix(PyMOLGlobals* G, float* m) {
  CScene* I = G->Scene;
  I->m_view.setRotMatrix(glm::make_mat4(m));
  SceneUpdateInvMatrix(G);
}
void SceneGetViewNormal(PyMOLGlobals* G, float* v) {
  CScene* I = G->Scene;
  copy3f(I->ViewNormal, v);
}
int SceneGetState(PyMOLGlobals* G) {
  return (SettingGetGlobal_i(G, cSetting_state) - 1);
}
float* SceneGetMatrix(PyMOLGlobals* G) {
  return const_cast<float*>(glm::value_ptr(G->Scene->m_view.rotMatrix()));
}

float* SceneGetPmvMatrix(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  multiply44f44f44f(SceneGetModelViewMatrixPtr(G),
                    SceneGetProjectionMatrixPtr(G), I->PmvMatrix);
  return (I->PmvMatrix);
}
GLFramebufferConfig SceneDrawBothGetConfig(PyMOLGlobals* G) {
  GLFramebufferConfig config;
  config.framebuffer = G->ShaderMgr->defaultBackbuffer.framebuffer;
  config.drawBuffer = SceneMustDrawBoth(G)
                          ? GL_BACK_LEFT
                          : G->ShaderMgr->defaultBackbuffer.drawBuffer;
  return config;
}
int SceneGetFrame(PyMOLGlobals* G) {
  if (MovieDefined(G))
    return (SettingGetGlobal_i(G, cSetting_frame) - 1);
  else
    return (SettingGetGlobal_i(G, cSetting_state) - 1);
}
/**
 * Returns the number of movie frames, or the number of states if no movie
 * is defined.
 */
int SceneCountFrames(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  int mov_len = MovieGetLength(G);
  I->HasMovie = (mov_len != 0);
  if (mov_len > 0) {
    I->NFrame = mov_len;
  } else {
    if (mov_len < 0) {
      // allows you to see cached movie even w/o object
      I->NFrame = -mov_len;
    } else {
      I->NFrame = 0;
    }

    for (const auto& obj : I->Obj) {
      int n = obj->getNFrame();
      if (n > I->NFrame) I->NFrame = n;
    }
  }
  PRINTFD(G, FB_Scene)
  " %s: leaving... I->NFrame %d\n", __func__, I->NFrame ENDFD return I->NFrame;
}
void SceneSetFrame(PyMOLGlobals* G, int mode, int frame) {
  CScene* I = G->Scene;
  int newFrame;
  int newState = 0;
  int movieCommand = false;
  int suppress = false;

  newFrame = SettingGetGlobal_i(G, cSetting_frame) - 1;
  PRINTFD(G, FB_Scene)
  " %s: entered.\n", __func__ ENDFD;
  switch (mode) {
    case -1: /* movie/frame override - go to this state absolutely! */
      newState = frame;
      break;
    case 0: /* absolute frame */
      newFrame = frame;
      break;
    case 1: /* relative frame */
      newFrame += frame;
      break;
    case 2: /* end */
      newFrame = I->NFrame - 1;
      break;
    case 3: /* middle with automatic movie command */
      newFrame = I->NFrame / 2;
      movieCommand = true;
      break;
    case 4: /* absolute with automatic movie command */
      newFrame = frame;
      movieCommand = true;
      break;
    case 5: /* relative with automatic movie command */
      newFrame += frame;
      movieCommand = true;
      break;
    case 6: /* end with automatic movie command */
      newFrame = I->NFrame - 1;
      movieCommand = true;
      break;
    case 7: /* absolute with forced movie command */
      newFrame = frame;
      movieCommand = true;
      break;
    case 8: /* relative with forced movie command */
      newFrame += frame;
      movieCommand = true;
      break;
    case 9: /* end with forced movie command */
      newFrame = I->NFrame - 1;
      movieCommand = true;
      break;
    case 10: /* seek forward to current scene (if present) */
    {
      frame = MovieSeekScene(G, true);
      if (frame >= 0) {
        newFrame = frame;
        movieCommand = true;
      } else {
        suppress = true;
      }
    } break;
  }
  if (!suppress) {
    SceneCountFrames(G);
    if (mode >= 0) {
      if (newFrame >= I->NFrame) newFrame = I->NFrame - 1;
      if (newFrame < 0) newFrame = 0;
      newState = MovieFrameToIndex(G, newFrame);
      if (newFrame == 0) {
        if (MovieMatrix(G, cMovieMatrixRecall)) {
          SceneAbortAnimation(G); /* if we have a programmed initial
                                     orientation, don't allow animation
                                     to override it */
        }
      }
      SettingSetGlobal_i(G, cSetting_frame, newFrame + 1);
      SettingSetGlobal_i(G, cSetting_state, newState + 1);
      ExecutiveInvalidateSelectionIndicatorsCGO(G);
      SceneInvalidatePicking(G);
      if (movieCommand) {
        MovieDoFrameCommand(G, newFrame);
        MovieFlushCommands(G);
      }
      if (SettingGetGlobal_b(G, cSetting_cache_frames))
        I->MovieFrameFlag = true;
    } else {
      SettingSetGlobal_i(G, cSetting_frame, newFrame + 1);
      SettingSetGlobal_i(G, cSetting_state, newState + 1);
      ExecutiveInvalidateSelectionIndicatorsCGO(G);
      SceneInvalidatePicking(G);
    }
    MovieSetScrollBarFrame(G, newFrame);
    SeqChanged(G);  // SceneInvalidate(G);
  }
  PRINTFD(G, FB_Scene)
  " %s: leaving...\n", __func__ ENDFD;
  OrthoInvalidateDoDraw(G);
}
void SceneDirty(PyMOLGlobals* G)

/* This means that the current image on the screen (and/or in the buffer)

   needs to be updated */
{
  CScene* I = G->Scene;

  PRINTFD(G, FB_Scene)
  " %s: called.\n", __func__ ENDFD;

  if (I) {
    if (!I->DirtyFlag) {
      I->DirtyFlag = true;
      /* SceneInvalidateCopy(G,false); */
      OrthoDirty(G);
    }
  }
}
void SceneChanged(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->ChangedFlag = true;
  SceneInvalidateCopy(G, false);
  SceneDirty(G);
  SeqChanged(G);
  PyMOL_NeedRedisplay(G->PyMOL);
}
Block* SceneGetBlock(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (I);
}
void SceneIdle(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  double renderTime;
  double minTime;
  int frameFlag = false;

  if (I->PossibleSingleClick == 2) {
    double now = UtilGetSeconds(G);
    double single_click_delay = I->SingleClickDelay;
    double diff = now - I->LastReleaseTime;
    if (diff > single_click_delay) {
      /* post a single click processing event */
      SceneDeferClickWhen(I, I->LastButton + P_GLUT_SINGLE_LEFT, I->LastWinX,
                          I->LastWinY, I->LastClickTime,
                          I->LastMod); /* push a click onto the queue */

      I->PossibleSingleClick = 0;
      OrthoDirty(G); /* force an update */
    }
  }
  if (!OrthoDeferredWaiting(G)) {
    if (MoviePlaying(G)) {
      renderTime = UtilGetSeconds(G) - I->LastFrameTime;
      {
        float fps = SettingGetGlobal_f(G, cSetting_movie_fps);
        if (fps <= 0.0F) {
          if (fps < 0.0)
            minTime = 0.0; /* negative fps means full speed */
          else             /* 0 fps means use movie_delay instead */
            minTime = SettingGetGlobal_f(G, cSetting_movie_delay) / 1000.0;
          if (minTime >= 0)
            fps = 1.0 / minTime;
          else
            fps = 1000.0F;
        } else {
          minTime = 1.0 / fps;
        }
        if (renderTime >= (minTime - I->LastFrameAdjust)) {
          float adjust = (renderTime - minTime);
          if ((fabs(adjust) < minTime) &&
              (fabs(I->LastFrameAdjust) < minTime)) {
            float new_adjust = (renderTime - minTime) + I->LastFrameAdjust;
            I->LastFrameAdjust =
                (new_adjust + fps * I->LastFrameAdjust) / (1 + fps);
          } else {
            I->LastFrameAdjust = 0.0F;
          }
          frameFlag = true;
        }
      }
    } else if (ControlRocking(G)) {
      renderTime = -I->LastSweepTime + UtilGetSeconds(G);
      minTime = SettingGetGlobal_f(G, cSetting_rock_delay) / 1000.0;
      if (renderTime >= minTime) {
        I->RenderTime = renderTime;
        SceneUpdateCameraRock(G, true);
      }
    }

    if (MoviePlaying(G) && frameFlag) {
      I->LastFrameTime = UtilGetSeconds(G);
      if ((SettingGetGlobal_i(G, cSetting_frame) - 1) == (I->NFrame - 1)) {
        if (SettingGetGlobal_b(G, cSetting_movie_loop)) {
          SceneSetFrame(G, 7, 0);
        } else
          MoviePlay(G, cMovieStop);
      } else {
        SceneSetFrame(G, 5, 1);
      }
      PyMOL_NeedRedisplay(G->PyMOL);
    }
  }
}
int SceneObjectAdd(PyMOLGlobals* G, pymol::CObject* obj) {
  CScene* I = G->Scene;
  obj->Enabled = true;
  I->Obj.push_back(obj);
  if (obj->type == cObjectGadget) {
    I->GadgetObjs.push_back(obj);
  } else {
    I->NonGadgetObjs.push_back(obj);
  }
  SceneCountFrames(G);
  SceneChanged(G);
  SceneInvalidatePicking(G);  // PYMOL-2793
  return 1;
}
int SceneObjectIsActive(PyMOLGlobals* G, pymol::CObject* obj) {
  int result = false;
  CScene* I = G->Scene;
  if (find(I->Obj.begin(), I->Obj.end(), obj) != I->Obj.end()) result = true;
  return result;
}

int SceneObjectDel(PyMOLGlobals* G, pymol::CObject* obj, int allow_purge) {
  CScene* I = G->Scene;
  int defer_builds_mode = SettingGetGlobal_i(G, cSetting_defer_builds_mode);

  if (!obj) { /* deletes all members */
    if (allow_purge && (defer_builds_mode >= 3)) {
      for (auto& obj : I->Obj) {
        /* purge graphics representation when no longer used */
        obj->invalidate(cRepAll, cRepInvPurge, -1);
      }
    }
    I->Obj.clear();
    I->GadgetObjs.clear();
    I->NonGadgetObjs.clear();
  } else {
    auto& obj_list =
        (obj->type == cObjectGadget) ? I->GadgetObjs : I->NonGadgetObjs;
    auto itg = find(obj_list.begin(), obj_list.end(), obj);
    if (itg != obj_list.end()) obj_list.erase(itg);

    auto it = find(I->Obj.begin(), I->Obj.end(), obj);
    if (it != I->Obj.end()) {
      if (allow_purge && (defer_builds_mode >= 3)) {
        /* purge graphics representation when no longer used */
        (*it)->invalidate(cRepAll, cRepInvPurge, -1);
      }
      obj->Enabled = false;
      I->Obj.erase(it);
    }
  }
  SceneCountFrames(G);
  SceneInvalidate(G);
  SceneInvalidatePicking(G);
  return 0;
}

bool SceneObjectRemove(PyMOLGlobals* G, pymol::CObject* obj) {
  if (obj == nullptr) {
    return true;
  }
  CScene* I = G->Scene;
  auto& obj_list =
      (obj->type == cObjectGadget) ? I->GadgetObjs : I->NonGadgetObjs;
  auto it = std::find(obj_list.begin(), obj_list.end(), obj);
  if (it == obj_list.end()) {
    return false;
  }
  obj_list.erase(it, obj_list.end());
  return true;
}

/*static unsigned int byte_max(unsigned int value)
{
  return (value>0xFF) ? 0xFF : value;
}
*/

#define SceneClickMargin DIP2PIXEL(2)
#define SceneTopMargin 0
#define SceneToggleMargin DIP2PIXEL(2)
#define SceneRightMargin 0
#define SceneToggleWidth DIP2PIXEL(17)
#define SceneToggleSize DIP2PIXEL(16)
#define SceneToggleTextShift DIP2PIXEL(4)
#define SceneTextLeftMargin DIP2PIXEL(1)
static void draw_button(int x2, int y2, int z, int w, int h, float* light,
                        float* dark, float* inside, CGO* orthoCGO) {
  if (orthoCGO) {
    CGOColorv(orthoCGO, light);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2, y2, z);
    CGOVertex(orthoCGO, x2, y2 + h, z);
    CGOVertex(orthoCGO, x2 + w, y2, z);
    CGOVertex(orthoCGO, x2 + w, y2 + h, z);
    CGOEnd(orthoCGO);
  } else {
    ImmBatch batch;
    batch.begin(GL_POLYGON);
    batch.color3fv(light);
    batch.vertex3i(x2, y2, z);
    batch.vertex3i(x2, y2 + h, z);
    batch.vertex3i(x2 + w, y2 + h, z);
    batch.vertex3i(x2 + w, y2, z);
    batch.end();
  }

  if (orthoCGO) {
    CGOColorv(orthoCGO, dark);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2 + 1, y2, z);
    CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, z);
    CGOVertex(orthoCGO, x2 + w, y2, z);
    CGOVertex(orthoCGO, x2 + w, y2 + h - 1, z);
    CGOEnd(orthoCGO);
  } else {
    ImmBatch batch;
    batch.begin(GL_POLYGON);
    batch.color3fv(dark);
    batch.vertex3i(x2 + 1, y2, z);
    batch.vertex3i(x2 + 1, y2 + h - 1, z);
    batch.vertex3i(x2 + w, y2 + h - 1, z);
    batch.vertex3i(x2 + w, y2, z);
    batch.end();
  }

  if (inside) {
    if (orthoCGO) {
      CGOColorv(orthoCGO, inside);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, z);
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, z);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, z);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, z);
      CGOEnd(orthoCGO);
    } else {
      ImmBatch batch;
      batch.begin(GL_POLYGON);
      batch.color3fv(inside);
      batch.vertex3i(x2 + 1, y2 + 1, z);
      batch.vertex3i(x2 + 1, y2 + h - 1, z);
      batch.vertex3i(x2 + w - 1, y2 + h - 1, z);
      batch.vertex3i(x2 + w - 1, y2 + 1, z);
      batch.end();
    }
  } else { /* rainbow */
    if (orthoCGO) {
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOColor(orthoCGO, 0.1F, 1.0F, 0.1F);  // green
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, z);
      CGOColor(orthoCGO, 1.0F, 1.0F, 0.1F);  // yellow
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, z);
      CGOColor(orthoCGO, 1.f, 0.1f, 0.1f);  // red
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, z);
      CGOColor(orthoCGO, 0.1F, 0.1F, 1.0F);  // blue
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, z);
      CGOEnd(orthoCGO);
    } else {
      ImmBatch batch;
      batch.begin(GL_POLYGON);
      batch.color3f(1.0F, 0.1F, 0.1F);
      batch.vertex3i(x2 + 1, y2 + 1, z);
      batch.color3f(0.1F, 1.0F, 0.1F);
      batch.vertex3i(x2 + 1, y2 + h - 1, z);
      batch.color3f(1.0F, 1.0F, 0.1F);
      batch.vertex3i(x2 + w - 1, y2 + h - 1, z);
      batch.color3f(0.1F, 0.1F, 1.0F);
      batch.vertex3i(x2 + w - 1, y2 + 1, z);
      batch.end();
    }
  }
}

/**
 * Update the G->Scene->SceneVLA names array which is used for scene buttons
 *
 * @param list List of scene names
 */
void SceneSetNames(PyMOLGlobals* G, const std::vector<std::string>& list) {
  CScene* I = G->Scene;
  I->SceneVec.clear();
  I->SceneVec.reserve(list.size());

  for (const auto& name : list) {
    I->SceneVec.emplace_back(name, false);
  }

  OrthoDirty(G);
}
static void SceneDrawButtons(Block* block, int draw_for_real, CGO* orthoCGO) {
  PyMOLGlobals* G = block->m_G;
  CScene* I = G->Scene;
  int x, y, xx, x2;
  float enabledColor[3] = {0.5F, 0.5F, 0.5F};
  float pressedColor[3] = {0.7F, 0.7F, 0.7F};
  float disabledColor[3] = {0.25F, 0.25F, 0.25F};
  float lightEdge[3] = {0.6F, 0.6F, 0.6F};
  float darkEdge[3] = {0.35F, 0.35F, 0.35F};
  int charWidth = DIP2PIXEL(8);
  int n_ent;
  int n_disp;
  int skip = 0;
  int row = -1;
  int lineHeight =
      DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
  int text_lift = (lineHeight / 2) - DIP2PIXEL(5);
  int op_cnt = 1;

  if (((G->HaveGUI && G->ValidContext) || (!draw_for_real)) &&
      ((block->rect.right - block->rect.left) > 6) && (!I->SceneVec.empty())) {
    int max_char;
    int nChar;
    I->ButtonsShown = true;

    /* do we have enough structures to warrant a scroll bar? */
    n_ent = I->SceneVec.size();

    n_disp =
        (((I->rect.top - I->rect.bottom) - (SceneTopMargin)) / lineHeight) - 1;
    if (n_disp < 1) n_disp = 1;

    {
      for (auto& elem : I->SceneVec) elem.drawn = false;
    }
    if (n_ent > n_disp) {
      int bar_maxed = I->m_ScrollBar.isMaxed();
      if (!I->ScrollBarActive) {
        I->m_ScrollBar.setLimits(n_ent, n_disp);
        if (bar_maxed) {
          I->m_ScrollBar.maxOut();
          I->NSkip = static_cast<int>(I->m_ScrollBar.getValue());
        } else {
          I->m_ScrollBar.setValue(0);
          I->NSkip = 0;
        }
      } else {
        I->m_ScrollBar.setLimits(n_ent, n_disp);
        if (bar_maxed) I->m_ScrollBar.maxOut();
        I->NSkip = static_cast<int>(I->m_ScrollBar.getValue());
      }
      I->ScrollBarActive = 1;

    } else {
      I->ScrollBarActive = 0;
      I->NSkip = 0;
    }

    max_char = (((I->rect.right - I->rect.left) -
                 (SceneTextLeftMargin + SceneRightMargin + 4)) -
                (op_cnt * SceneToggleWidth));
    if (I->ScrollBarActive) {
      max_char -= (SceneScrollBarMargin + SceneScrollBarWidth);
    }
    max_char /= charWidth;

    if (I->ScrollBarActive) {
      I->m_ScrollBar.setBox(
          I->rect.top - SceneScrollBarMargin,
          I->rect.left + SceneScrollBarMargin, I->rect.bottom + 2,
          I->rect.left + SceneScrollBarMargin + SceneScrollBarWidth);
      if (draw_for_real) I->m_ScrollBar.draw(orthoCGO);
    }

    skip = I->NSkip;
    x = I->rect.left + SceneTextLeftMargin;

    /*    y = ((I->rect.top-lineHeight)-SceneTopMargin)-lineHeight; */

    {
      int n_vis = n_disp;
      if (n_ent < n_vis) n_vis = n_ent;
      y = (I->rect.bottom + SceneBottomMargin) + (n_vis - 1) * lineHeight;
    }

    /*    xx = I->rect.right-SceneRightMargin-SceneToggleWidth*(cRepCnt+op_cnt);
     */
    xx = I->rect.right - SceneRightMargin - SceneToggleWidth * (op_cnt);

    if (I->ScrollBarActive) {
      x += SceneScrollBarWidth + SceneScrollBarMargin;
    }
    {
      int i;

      for (i = 0; i < n_ent; i++) {
        if (skip) {
          skip--;
        } else {
          row++;
          x2 = xx;
          nChar = max_char;

          if ((x - SceneToggleMargin) - (xx - SceneToggleMargin) > -10) {
            x2 = x + 10;
          }
          {
            float toggleColor[3] = {0.5F, 0.5F, 1.0F};

            if (draw_for_real) {
              glColor3fv(toggleColor);

              TextSetColor(G, I->TextColor);
              TextSetPos2i(G, x + DIP2PIXEL(2), y + text_lift);
            }
            {
              const char* cur_name =
                  SettingGetGlobal_s(G, cSetting_scene_current_name);
              auto elem = &I->SceneVec[i];
              int item = I->NSkip + row;
              auto c = elem->name.c_str();
              int len = static_cast<int>(elem->name.size());

              x2 = xx;
              if (len > max_char) len = max_char;
              x2 = x + len * charWidth + DIP2PIXEL(6);

              /* store rectangles for finding clicks */

              elem->drawn = true;

              elem->rect = pymol::Rect<int>(x, x2, y, y + lineHeight);

              if (I->ButtonMargin < x2) I->ButtonMargin = x2;

              if (draw_for_real) {
                if ((item == I->Pressed) && (item == I->Over)) {
                  draw_button(x, y, 0, (x2 - x) - 1, (lineHeight - 1),
                              lightEdge, darkEdge, pressedColor, orthoCGO);
                } else if (cur_name && elem->name == cur_name) {
                  draw_button(x, y, 0, (x2 - x) - 1, (lineHeight - 1),
                              lightEdge, darkEdge, enabledColor, orthoCGO);
                } else {
                  draw_button(x, y, 0, (x2 - x) - 1, (lineHeight - 1),
                              lightEdge, darkEdge, disabledColor, orthoCGO);
                }

                TextSetColor(G, I->TextColor);

                if (c) {
                  while (*c) {
                    if ((nChar--) > 0)
                      TextDrawChar(G, *(c++), orthoCGO);
                    else
                      break;
                  }
                }
              }
            }
          }
          y -= lineHeight;
          if (y < (I->rect.bottom)) break;
        }
      }
    }
    I->HowFarDown = y;
    I->ButtonsValid = true;
  }
}

void CScene::draw(
    CGO* orthoCGO) /* returns true if scene was drawn (using a cached image) */
{
  PyMOLGlobals* G = m_G;
  CScene* I = G->Scene;
  int drawn = false;

  if (G->HaveGUI && G->ValidContext) {
    I->ButtonsShown = false;

    drawn = SceneDrawImageOverlay(G, 0, orthoCGO);

    if (SettingGetGlobal_b(G, cSetting_scene_buttons)) {
      SceneDrawButtons(this, true, orthoCGO);
    } else {
      I->ButtonMargin = 0;
    }
  }
  if (drawn)
    OrthoDrawWizardPrompt(
        G, orthoCGO); /* ugly hack necessitated because wizard
                              prompt is overwritten when image is drawn */
}

float ScenePushRasterMatrix(PyMOLGlobals* G, float* v) {
  float scale = SceneGetExactScreenVertexScale(G, v);
  return scale;
}

void ScenePopRasterMatrix(PyMOLGlobals* G) { /* TODO */ }

/**
 * Compose the ModelViewMatrix from Pos, RotMatrix and Origin
 * See also: CScene.ModMatrix (queried from OpenGL)
 *
 * @param[out] modelView 4x4 matrix
 */
static void SceneComposeModelViewMatrix(CScene* I, float* modelView) {
  identity44f(modelView);
  const auto& pos = I->m_view.pos();
  const auto& ori = I->m_view.origin();
  MatrixTranslateC44f(modelView, pos.x, pos.y, pos.z);
  MatrixMultiplyC44f(glm::value_ptr(I->m_view.rotMatrix()), modelView);
  MatrixTranslateC44f(modelView, -ori.x, -ori.y, -ori.z);
}
void SceneGetEyeNormal(PyMOLGlobals* G, float* v1, float* normal) {
  CScene* I = G->Scene;
  float p1[4], p2[4];
  float modelView[16];

  SceneComposeModelViewMatrix(I, modelView);

  copy3f(v1, p1);
  p1[3] = 1.0;
  MatrixTransformC44f4f(modelView, p1, p2); /* modelview transformation */
  copy3f(p2, p1);
  normalize3f(p1);
  MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), p1, p2);
  invert3f3f(p2, normal);
}

/**
 * Get the depth (camera space Z) of v1
 *
 * @param v1 point (3f) in world space or nullptr (= origin)
 */
float SceneGetRawDepth(PyMOLGlobals* G, const float* v1) {
  CScene* I = G->Scene;
  float vt[3];
  float modelView[16];

  if (!v1 || SettingGetGlobal_i(G, cSetting_ortho)) return -I->m_view.pos().z;

  SceneComposeModelViewMatrix(I, modelView);

  MatrixTransformC44f3f(modelView, v1, vt);
  return -vt[2];
}

/**
 * Get the depth (camera space Z) of v1 in normalized clip space
 * from 0.0 (near) to 1.0 (far)
 *
 * @param v1 point (3f) in world space or nullptr (= origin)
 */
float SceneGetDepth(PyMOLGlobals* G, const float* v1) {
  CScene* I = G->Scene;
  float rawDepth = SceneGetRawDepth(G, v1);
  return ((rawDepth - I->m_view.m_clipSafe().m_front) /
          (I->m_view.m_clipSafe().m_back - I->m_view.m_clipSafe().m_front));
}
/**
 * Get the angstrom per pixel factor at v1. If v1 is nullptr, return the
 * factor at the origin, but clamped to an empirical positive value.
 *
 * @param v1 point (3f) in world space or nullptr (= origin)
 */
float SceneGetScreenVertexScale(PyMOLGlobals* G, const float* v1)

/* does not require OpenGL-provided matrices */
{
  float depth = SceneGetRawDepth(G, v1);
  float ratio = depth * GetFovWidth(G) / G->Scene->Height;

  if (!v1 && ratio < R_SMALL4)
    // origin depth, return a safe clipped value (origin must not be
    // behind or very close in front of the camera)
    ratio = R_SMALL4;

  return ratio;
}

int CScene::click(int button, int x, int y,
                  int mod)  // Originally SceneDeferClick!!
{
  return SceneDeferClickWhen(this, button, x, y, UtilGetSeconds(m_G), mod);
}

static int SceneDeferClickWhen(Block* block, int button, int x, int y,
                               double when, int mod) {
  PyMOLGlobals* G = block->m_G;
  std::function<void()> deferred = [=]() {
    SceneClick(block, button, x, y, mod, when);
  };
  OrthoDefer(G, std::move(deferred));
  return 1;
}

int CScene::drag(int x, int y, int mod)  // Originally SceneDeferDrag
{
  PyMOLGlobals* G = m_G;
  auto when = UtilGetSeconds(G);
  std::function<void()> deferred = [=]() { SceneDrag(this, x, y, mod, when); };
  OrthoDefer(G, std::move(deferred));
  return 1;
}

// static int SceneDeferredRelease(DeferredMouse * dm)
//{
//   SceneRelease(dm->block, dm->button, dm->x, dm->y, dm->mod, dm->when);
//   return 1;
// }

int CScene::release(int button, int x, int y,
                    int mod)  // Originally SceneDeferRelease
{
  PyMOLGlobals* G = m_G;
  auto when = UtilGetSeconds(G);
  std::function<void()> deferred = [=]() {
    SceneRelease(this, button, x, y, mod, when);
  };

  OrthoDefer(G, std::move(deferred));
  return 1;
}
void SceneFree(PyMOLGlobals* G) {
  CScene* I = G->Scene;

  CGOFree(I->AlphaCGO);
  CGOFree(I->offscreenCGO);
  CGOFree(I->offscreenOIT_CGO);
  CGOFree(I->offscreenOIT_CGO_copy);
  I->m_slots.clear();
  I->Obj.clear();
  I->GadgetObjs.clear();
  I->NonGadgetObjs.clear();

  ScenePurgeImage(G);
  CGOFree(G->DebugCGO);
  delete G->Scene;
}
void SceneResetMatrix(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->m_view.setRotMatrix(glm::mat4(1.0f));
  SceneUpdateInvMatrix(G);
}
void SceneSetDefaultView(PyMOLGlobals* G) {
  CScene* I = G->Scene;

  I->m_view.setRotMatrix(glm::mat4(1.0f));
  SceneUpdateInvMatrix(G);

  I->ViewNormal[0] = 0.0F;
  I->ViewNormal[1] = 0.0F;
  I->ViewNormal[2] = 1.0F;

  I->m_view.setPos(0.0f, 0.0f, -50.0f);
  I->m_view.setOrigin(0.0f, 0.0f, 0.0f);

  I->m_view.m_clip().m_front = 40.0F;
  I->m_view.m_clip().m_back = 100.0F;
  UpdateFrontBackSafe(I);

  I->Scale = 1.0F;
}

int SceneReinitialize(PyMOLGlobals* G) {
  int ok = true;
  SceneSetDefaultView(G);
  SceneCountFrames(G);
  SceneSetFrame(G, 0, 0);
  SceneInvalidate(G);
  G->Scene->SceneVec.clear();
  return (ok);
}
int SceneInit(PyMOLGlobals* G) {
  CScene* I = nullptr;
  I = (G->Scene = new CScene(G));
  if (I) {
    assert(!I->RovingDirtyFlag);
    assert(I->DirtyFlag);

    /* all defaults to zero, so only initialize non-zero elements */
    G->DebugCGO = CGONew(G);

    I->LastClickTime = UtilGetSeconds(G);

    SceneSetDefaultView(G);

    I->active = true;

    OrthoAttach(G, I, cOrthoScene);

    I->LastRender = UtilGetSeconds(G);
    I->LastFrameTime = UtilGetSeconds(G);

    I->LastSweepTime = UtilGetSeconds(G);

    SceneRestartFrameTimer(G);
    SceneRestartPerfTimer(G);

    /* scene list */

    I->Pressed = -1;
    I->Over = -1;

    return 1;
  } else
    return 0;
}
void CScene::reshape(int width, int height) {
  PyMOLGlobals* G = m_G;
  CScene* I = G->Scene;

  if (I->margin.right) {
    width -= I->margin.right;
    if (width < 1) width = 1;
  }

  if (I->margin.top) {
    height -= I->margin.top;
  }

  I->Width = width;

  I->Height = height;

  I->rect.top = I->Height;
  I->rect.left = 0;
  I->rect.bottom = 0;
  I->rect.right = I->Width;

  if (I->margin.bottom) {
    height -= I->margin.bottom;
    if (height < 1) height = 1;
    I->Height = height;
    I->rect.bottom = I->rect.top - I->Height;
  }
  SceneDirty(G);

  if (I->CopyType && (!I->CopyForced)) {
    SceneInvalidateCopy(G, false);
  }
  /*MovieClearImages(G); */
  MovieSetSize(G, I->Width, I->Height);
  SceneInvalidateStencil(G);
}
void SceneResetNormal(PyMOLGlobals* G, int lines) {
  CScene* I = G->Scene;
  if (G->HaveGUI && G->ValidContext) {
    if (lines)
      glNormal3fv(I->LinesNormal);
    else
      glNormal3fv(I->ViewNormal);
  }
}

void SceneResetNormalCGO(PyMOLGlobals* G, CGO* cgo, int lines) {
  CScene* I = G->Scene;
  if (G->HaveGUI && G->ValidContext) {
    if (lines)
      CGONormalv(cgo, I->LinesNormal);
    else
      CGONormalv(cgo, I->ViewNormal);
  }
}

void SceneResetNormalToViewVector(PyMOLGlobals* G, short use_shader) {
  auto modMatrix = SceneGetModelViewMatrixPtr(G);
  if (G->HaveGUI && G->ValidContext) {
    glVertexAttrib3f(VERTEX_NORMAL, modMatrix[2], modMatrix[6], modMatrix[10]);
  }
}

void SceneResetNormalUseShader(PyMOLGlobals* G, int lines, short use_shader) {
  CScene* I = G->Scene;
  if (G->HaveGUI && G->ValidContext) {
    if (lines)
      glVertexAttrib3fv(VERTEX_NORMAL, I->LinesNormal);
    else
      glVertexAttrib3fv(VERTEX_NORMAL, I->ViewNormal);
  }
}

void SceneGetResetNormal(PyMOLGlobals* G, float* normal, int lines) {
  CScene* I = G->Scene;
  float* norm;
  if (G->HaveGUI && G->ValidContext) {
    if (lines)
      norm = I->LinesNormal;
    else
      norm = I->ViewNormal;
    normal[0] = norm[0];
    normal[1] = norm[1];
    normal[2] = norm[2];
  }
}

void SceneUpdateAnimation(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  int rockFlag = false;
  int dirty = false;
  int movie_rock = SettingGetGlobal_i(G, cSetting_movie_rock);

  if (movie_rock < 0) movie_rock = ControlRocking(G);

  if (MoviePlaying(G) && movie_rock) {
    if (MovieGetRealtime(G) &&
        !SettingGetGlobal_b(G, cSetting_movie_animate_by_frame)) {
      I->RenderTime = UtilGetSeconds(G) - I->LastSweepTime;
      rockFlag = true;
      dirty = true; /* force a subsequent update */
    } else {
      float fps = SceneGetFPS(G); /* guaranteed to be >= 0.0F */
      if (fps > 0.0F) {
        int rock_frame = SceneGetFrame(G);
        if (rock_frame != I->RockFrame) {
          I->RockFrame = rock_frame;
          rockFlag = true;
          I->RenderTime = 1.0 / fps;
        }
      } else {
        I->RenderTime = UtilGetSeconds(G) - I->LastSweepTime;
        rockFlag = true;
      }
    }
  } else
    dirty = true;

  if (I->cur_ani_elem < I->n_ani_elem) { /* play motion animation */
    double now;

    int cur = I->cur_ani_elem;

    if (I->AnimationStartFlag) {
      /* allow animation timing to lag since it may take a few seconds
         to get here given geometry updates, etc. */

      I->AnimationLagTime = UtilGetSeconds(G) - I->AnimationStartTime;
      I->AnimationStartFlag = false;
    }

    if ((!MoviePlaying(G)) ||
        ((MovieGetRealtime(G) &&
          !SettingGetGlobal_b(G, cSetting_movie_animate_by_frame)))) {
      now = UtilGetSeconds(G) - I->AnimationLagTime;
    } else {
      float fps = SceneGetFPS(G); /* guaranteed to be >= 0.0F */
      int frame = SceneGetFrame(G);
      int n_frame = 0;

      cur = 0; /* allow backwards interpolation */
      if (frame >= I->AnimationStartFrame) {
        n_frame = frame - I->AnimationStartFrame;
      } else {
        n_frame = frame + (I->NFrame - I->AnimationStartFrame);
      }
      now = I->AnimationStartTime + n_frame / fps;
    }

    while (I->ani_elem[cur].timing < now) {
      cur++;
      if (cur >= I->n_ani_elem) {
        cur = I->n_ani_elem;
        break;
      }
    }
    I->cur_ani_elem = cur;
    SceneFromViewElem(G, I->ani_elem + cur, dirty);
    OrthoDirty(G);
  }
  if (rockFlag && (I->RenderTime != 0.0)) {
    SceneUpdateCameraRock(G, dirty);
  }
}

int SceneGetDrawFlag(GridInfo* grid, int* slot_vla, int slot) {
  int draw_flag = false;
  if (grid && grid->active) {
    switch (grid->mode) {
      case GridMode::ByObject: /* assigned grid slots (usually by group) */
      {
        if (((slot < 0) && grid->slot) || ((slot == 0) && (grid->slot == 0)) ||
            (slot_vla && (slot >= 0 && slot_vla[slot] == grid->slot))) {
          draw_flag = true;
        }
      } break;
      case GridMode::ByObjectStates:
      case GridMode::ByObjectByState:
        draw_flag = true;
        break;
    }
  } else {
    draw_flag = true;
  }
  return draw_flag;
}

int SceneGetDrawFlagGrid(PyMOLGlobals* G, GridInfo* grid, int slot) {
  CScene* I = G->Scene;
  return SceneGetDrawFlag(grid, I->m_slots.data(), slot);
}
void SceneCopy(PyMOLGlobals* G, GLFramebufferConfig config, int force,
               int entire_window) {
  CScene* I = G->Scene;

  if (force || (!(I->StereoMode ||
                  SettingGetGlobal_b(G, cSetting_stereo_double_pump_mono) ||
                  I->ButtonsShown))) {
    /* no copies while in stereo mode */
    if (force || ((!I->DirtyFlag) && (!I->CopyType))) {
      Rect2D rect;
      if (entire_window) {
        rect = OrthoGetRect(G);
      } else {
        rect = Rect2D{{}, SceneGetExtent(G)};
      }
      ScenePurgeImage(G);
      auto imgData = G->ShaderMgr->readPixelsFrom(G, rect, config);
      if (!imgData.empty()) {
        I->Image = std::make_shared<pymol::Image>(rect.extent.width,
                                                  rect.extent.height);
        I->Image->setVecData(std::move(imgData));
      }
      I->CopyType = true;
      I->Image->m_needs_alpha_reset = true;
      I->CopyForced = force;
    }
  }
}

struct CObjectUpdateThreadInfo {
  pymol::CObject* obj;
};

void SceneObjectUpdateThread(CObjectUpdateThreadInfo* T) {
  if (T->obj) {
    T->obj->update();
  }
}

static void SceneObjectUpdateSpawn(PyMOLGlobals* G,
                                   CObjectUpdateThreadInfo* Thread,
                                   int n_thread, int n_total) {
  if (n_total == 1) {
    SceneObjectUpdateThread(Thread);
  } else if (n_total) {
    int blocked;
    PyObject* info_list;
    int a, n = 0;
    blocked = PAutoBlock(G);

    PRINTFB(G, FB_Scene, FB_Blather)
    " Scene: updating objects with %d threads...\n", n_thread ENDFB(G);
    info_list = PyList_New(n_total);
    for (a = 0; a < n_total; a++) {
      PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
      n++;
    }
    PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->cmd, "_object_update_spawn", "Oi",
                                 info_list, n_thread));
    Py_DECREF(info_list);
    PAutoUnblock(G, blocked);
  }
}

static void SceneStencilCheck(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  if (I->StereoMode == cStereo_stencil_by_row) {
    int bottom = 0;

#ifndef _PYMOL_PRETEND_GLUT
    if (G->Main)
      bottom = p_glutGet(P_GLUT_WINDOW_Y) + p_glutGet(P_GLUT_WINDOW_HEIGHT);
#endif

    int parity = bottom & 0x1;
    if (parity != I->StencilParity) {
      I->StencilValid = false;
      I->StencilParity = parity;
      SceneDirty(G);
    }
  }
}
void SceneUpdate(PyMOLGlobals* G, int force) {
  CScene* I = G->Scene;

  int cur_state = SettingGetGlobal_i(G, cSetting_state) - 1;
  int defer_builds_mode = SettingGetGlobal_i(G, cSetting_defer_builds_mode);

  PRINTFD(G, FB_Scene)
  " SceneUpdate: entered.\n" ENDFD;

  OrthoBusyPrime(G);
  WizardDoPosition(G, false);
  WizardDoView(G, false);
  EditorUpdate(G);
  SceneStencilCheck(G);

  if (defer_builds_mode == 0) {
    if (SettingGetGlobal_i(G, cSetting_draw_mode) == -2) {
      defer_builds_mode = 1;
    }
  }

  if (force || I->ChangedFlag ||
      ((cur_state != I->LastStateBuilt) && (defer_builds_mode > 0))) {
    SceneCountFrames(G);

    if (force || (defer_builds_mode != 5)) { /* mode 5 == immediate mode */

      PyMOL_SetBusy(G->PyMOL,
                    true); /*  race condition -- may need to be fixed */

      /* update all gadgets first (single-threaded since they're thread-unsafe)
       */
      for (auto& GadgetObj : I->GadgetObjs) {
        GadgetObj->update();
      }

      {
        int n_thread = SettingGetGlobal_i(G, cSetting_max_threads);
        int multithread = SettingGetGlobal_i(G, cSetting_async_builds);
        if (multithread && (n_thread > 1)) {
          int min_start = -1;
          int max_stop = -1;
          int n_obj = 0;
          for (auto& obj : I->Obj) {
            int start = 0;
            n_obj++;
            int stop = obj->getNFrame();
            /* set start/stop to define the range for this object
             * depending upon various build settings */
            ObjectAdjustStateRebuildRange(obj, &start, &stop);
            if (min_start < 0) {
              min_start = start;
              max_stop = stop;
            } else {
              if (min_start > start) min_start = start;
              if (max_stop < stop) max_stop = stop;
            }
          }

          int n_frame = max_stop - min_start;

          if (n_frame > n_thread) {
            n_thread = 1;
            /* prevent n_thread * n_thread -- only multithread within
               individual object states (typically more balanced) */
          } else if (n_frame > 1) {
            n_thread = n_thread / n_frame;
          }

          if (n_thread < 1) n_thread = 1;
        }

        /* Note: we might want to optimize this by doing multi-threaded updates
           for all objects. */
        if (multithread && (n_thread > 1)) {
          /* multi-threaded geometry update */
          int cnt = I->NonGadgetObjs.size();

          if (cnt) {
            CObjectUpdateThreadInfo* thread_info =
                pymol::malloc<CObjectUpdateThreadInfo>(cnt);
            if (thread_info) {
              cnt = 0;
              for (auto& NonGadgetObj : I->NonGadgetObjs) {
                thread_info[cnt++].obj = NonGadgetObj;
              }
              SceneObjectUpdateSpawn(G, thread_info, n_thread, cnt);
              FreeP(thread_info);
            }
          }
        } else
          /* single-threaded update */
          for (auto& obj : I->Obj) {
            obj->update();
          }
      }
      PyMOL_SetBusy(G->PyMOL,
                    false); /*  race condition -- may need to be fixed */
    } else { /* defer builds mode == 5 -- for now, only update non-molecular
                objects */
      /* single-threaded update */
      for (auto& obj : I->Obj) {
        if (obj->type != cObjectMolecule) {
          obj->update();
        }
      }
    }

    I->ChangedFlag = false;

    if ((defer_builds_mode >= 2) && (force || (defer_builds_mode != 5)) &&
        (cur_state != I->LastStateBuilt)) {
      /* purge graphics representation when no longer used */
      if (I->LastStateBuilt >= 0) {
        for (auto it = I->Obj.begin(); it != I->Obj.end(); ++it) {
          if ((*it)->type != cObjectMolecule || force ||
              defer_builds_mode != 5) {
            int static_singletons = SettingGet_b(
                G, (*it)->Setting.get(), nullptr, cSetting_static_singletons);
            int async_builds = SettingGet_b(G, (*it)->Setting.get(), nullptr,
                                            cSetting_async_builds);
            int max_threads = SettingGet_i(G, (*it)->Setting.get(), nullptr,
                                           cSetting_max_threads);
            int nFrame = 0;
            nFrame = (*it)->getNFrame();
            if ((nFrame > 1) || (!static_singletons)) {
              int start = I->LastStateBuilt;
              int stop = start + 1;
              int ste;
              if (async_builds && (max_threads > 1)) {
                if ((start / max_threads) == (cur_state / max_threads)) {
                  stop = start; /* don't purge current batch */
                } else {
                  int base = start / max_threads; /* now purge previous batch */
                  start = base * max_threads;
                  stop = (base + 1) * max_threads;
                }
              }
              for (ste = start; ste < stop; ste++) {
                (*it)->invalidate(cRepAll, cRepInvPurge, ste);
              }
            }
          }
        }
      }
    }
    I->LastStateBuilt = cur_state;
    WizardDoScene(G);
    if (!MovieDefined(G)) {
      if (SettingGetGlobal_i(G, cSetting_frame) != (cur_state + 1))
        SettingSetGlobal_i(G, cSetting_frame, (cur_state + 1));
    }
  }

  PRINTFD(G, FB_Scene)
  " %s: leaving...\n", __func__ ENDFD;
}

float SceneGetSpecularValue(PyMOLGlobals* G, float spec, int limit) {
  int n_light = SettingGetGlobal_i(G, cSetting_spec_count);
  if (n_light < 0) n_light = SettingGetGlobal_i(G, cSetting_light_count);
  if (n_light > limit) n_light = limit;
  if (n_light > 2) {
    spec = spec / pow(n_light - 1, 0.6F);
  }
  return (spec > 1.F) ? 1.F : (spec < 0.F) ? 0.F : spec;
}

float SceneGetReflectScaleValue(PyMOLGlobals* G, int limit) {
  float result = 1.0F;
  int n_light = SettingGetGlobal_i(G, cSetting_light_count);
  if (n_light > limit) n_light = limit;
  if (n_light > 1) {
    float tmp[3];
    float sum = 0.0F;
    for (int i = 0; i < n_light - 1; ++i) {
      copy3f(SettingGetGlobal_3fv(G, light_setting_indices[i]), tmp);
      normalize3f(tmp);
      sum += 1.f - tmp[2];
    }
    sum *= 0.5;
    return result / sum;
  }
  return result;
}

/**
 * Get specular and shininess, adjusted to the number of lights.
 *
 * @param[out] ptr_spec              specular for lights 2-N
 * @param[out] ptr_spec_power        shininess for lights 2-N
 * @param[out] ptr_spec_direct       specular for light 1
 * @param[out] ptr_spec_direct_power shininess for light 1
 * @param limit number of lights (e.g. `light_count` setting)
 */
void SceneGetAdjustedLightValues(PyMOLGlobals* G, float* ptr_spec,
                                 float* ptr_spec_power, float* ptr_spec_direct,
                                 float* ptr_spec_direct_power, int limit) {
  float specular = SettingGetGlobal_f(G, cSetting_specular);
  if (specular == 1.0F)
    specular = SettingGetGlobal_f(G, cSetting_specular_intensity);
  if (specular < R_SMALL4) specular = 0.0F;

  float spec_power = SettingGetGlobal_f(G, cSetting_spec_power);
  if (spec_power < 0.0F) spec_power = SettingGetGlobal_f(G, cSetting_shininess);

  float spec_reflect = SettingGetGlobal_f(G, cSetting_spec_reflect);
  if (spec_reflect < 0.0F) spec_reflect = specular;

  float spec_direct = SettingGetGlobal_f(G, cSetting_spec_direct);
  if (spec_direct < 0.0F) spec_direct = specular;

  float spec_direct_power = SettingGetGlobal_f(G, cSetting_spec_direct_power);
  if (spec_direct_power < 0.0F) spec_direct_power = spec_power;

  *ptr_spec = SceneGetSpecularValue(G, spec_reflect, limit);
  *ptr_spec_power = spec_power;
  *ptr_spec_direct = spec_direct > 1.F ? 1.F : spec_direct;
  *ptr_spec_direct_power = spec_direct_power;
}

/*
 * Shader attribute names
 */
#define TEMPLATE(i) "g_LightSource[" #i "].position"
const char* lightsource_position_names[] = {
    TEMPLATE(0), TEMPLATE(1), TEMPLATE(2), TEMPLATE(3), TEMPLATE(4),
    TEMPLATE(5), TEMPLATE(6), TEMPLATE(7), TEMPLATE(8), TEMPLATE(9)};
#undef TEMPLATE

#define TEMPLATE(i) "g_LightSource[" #i "].diffuse"
const char* lightsource_diffuse_names[] = {
    TEMPLATE(0), TEMPLATE(1), TEMPLATE(2), TEMPLATE(3), TEMPLATE(4),
    TEMPLATE(5), TEMPLATE(6), TEMPLATE(7), TEMPLATE(8), TEMPLATE(9)};
#undef TEMPLATE

/**
 * Sets up lighting for immediate mode if shaderPrg=nullptr, otherwise
 * sets lighting uniforms for the given shader program.
 *
 * Supports up to light_count=8
 */
void SceneProgramLighting(PyMOLGlobals* G, CShaderPrg* shaderPrg) {
  /* load up the light positions relative to the camera while
     MODELVIEW still has the identity */
  int n_light = glm::clamp(SettingGetGlobal_i(G, cSetting_light_count), 0, 8);
  int spec_count = SettingGetGlobal_i(G, cSetting_spec_count);
  float direct = SettingGetGlobal_f(G, cSetting_direct);
  float reflect = SettingGetGlobal_f(G, cSetting_reflect) *
                  SceneGetReflectScaleValue(G, n_light);
  float spec[4];
  float diff[4];
  const float zero[4] = {0.0F, 0.0F, 0.0F, 1.0F};
  float vv[4] = {0.F, 0.F, 1.F, 0.F};  // position
  float spec_value, shine, spec_direct, spec_direct_power;

  SceneGetAdjustedLightValues(G, &spec_value, &shine, &spec_direct,
                              &spec_direct_power, n_light);

  if (n_light < 2) {
    direct += reflect;
    if (direct > 1.0F) direct = 1.0F;
  }

  if (spec_count < 0) {
    spec_count = n_light;
  }

  // light 0

  white4f(diff, SettingGetGlobal_f(G, cSetting_ambient));

  {
    shaderPrg->Set4fv("g_LightModel.ambient", diff);
    white4f(diff, (direct > R_SMALL4) ? direct : 0.f);
    shaderPrg->Set4fv(lightsource_diffuse_names[0], diff);
    shaderPrg->Set4fv(lightsource_position_names[0], vv);
  }

  // light 1-N

  white4f(spec, spec_value);
  white4f(diff, reflect);

  for (int i = 1; i < n_light; ++i) {
    // normalized/inverted light direction
    copy3f(SettingGetGlobal_3fv(G, light_setting_indices[i - 1]), vv);
    normalize3f(vv);
    invert3f(vv);

    {
      shaderPrg->Set4fv(lightsource_position_names[i], vv);
      shaderPrg->Set4fv(lightsource_diffuse_names[i], diff);
    }
  }
}

/**
 * Set up the Scene Fog* member variables and immediate mode fog (incl.
 * gl_Fog struct for non-ES2 shaders)
 */
int SceneSetFog(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  int fog_active = false;
  float fog_density = SettingGetGlobal_f(G, cSetting_fog);
  I->FogStart =
      (I->m_view.m_clipSafe().m_back - I->m_view.m_clipSafe().m_front) *
          SettingGetGlobal_f(G, cSetting_fog_start) +
      I->m_view.m_clipSafe().m_front;
  if ((fog_density > R_SMALL8) && (fog_density != 1.0F)) {
    I->FogEnd = I->FogStart +
                (I->m_view.m_clipSafe().m_back - I->FogStart) / fog_density;
  } else {
    I->FogEnd = I->m_view.m_clipSafe().m_back;
  }

  if (SettingGetGlobal_b(G, cSetting_depth_cue) && fog_density != 0.0F) {
    fog_active = true;
  }

  return fog_active;
}

/**
 * Set the g_Fog_* uniforms for ES2 shaders
 */
void SceneSetFogUniforms(PyMOLGlobals* G, CShaderPrg* shaderPrg) {
  CScene* I = G->Scene;
  if (shaderPrg) {
    float fogScale = 1.0f / (I->FogEnd - I->FogStart);
    shaderPrg->Set1f("g_Fog_end", I->FogEnd);
    shaderPrg->Set1f("g_Fog_scale", fogScale);
  }
}

void SceneSetupGLPicking(PyMOLGlobals* G) {
  /* picking mode: we want flat, unshaded, unblended, unsmooth colors */

  glDisable(GL_FOG);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_DITHER);
  glDisable(GL_BLEND);
  glDisable(GL_POLYGON_SMOOTH);
  if (G->Option->multisample) glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
  glShadeModel(GL_FLAT);
}
void SceneRestartFrameTimer(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->LastFrameTime = UtilGetSeconds(G);
}

static void SceneRestartPerfTimer(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->LastRender = UtilGetSeconds(G);
  I->RenderTime = 0.0;
}

void SceneRestartSweepTimer(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->LastSweep = 0.0F; /* continue to defer rocking until this is done */
  I->LastSweepX = 0.0F;
  I->LastSweepY = 0.0F;
  I->SweepTime = 0.0;
  I->LastSweepTime = UtilGetSeconds(G);
  SceneRestartPerfTimer(G);
}
void ScenePrepareMatrix(PyMOLGlobals* G, int mode, int stereo_mode /* = 0 */) {
  CScene* I = G->Scene;

  float stAng, stShift;
  const auto& pos = I->m_view.pos();
  const auto& ori = I->m_view.origin();

  {
    if (!mode) {
      SceneComposeModelViewMatrix(I, SceneGetModelViewMatrixPtr(G));
    } else {
      /* stereo */
      float tmpMatrix[16];
      stAng = SettingGetGlobal_f(G, cSetting_stereo_angle);  // * cPI / 180.f;
      stShift = SettingGetGlobal_f(G, cSetting_stereo_shift);

      stShift = (float)(stShift * fabs(pos.z) / 100.0);
      stAng = (float)(-stAng * atan(stShift / fabs(pos.z)) / 2.f);

      if (mode == 2) { /* left hand */
        stAng = -stAng;
        stShift = -stShift;
      }

      PRINTFD(G, FB_Scene)
      " StereoMatrix-Debug: mode %d stAng %8.3f stShift %8.3f \n", mode, stAng,
          stShift ENDFD;
      identity44f(tmpMatrix);
      I->modelViewMatrix = glm::mat4(1.0f);
      MatrixRotateC44f(SceneGetModelViewMatrixPtr(G), stAng, 0.f, 1.f, 0.f);
      MatrixTranslateC44f(tmpMatrix, pos.x + stShift, pos.y, pos.z);
      MatrixMultiplyC44f(tmpMatrix, SceneGetModelViewMatrixPtr(G));
      MatrixMultiplyC44f(glm::value_ptr(I->m_view.rotMatrix()),
                         SceneGetModelViewMatrixPtr(G));
      MatrixTranslateC44f(SceneGetModelViewMatrixPtr(G), -ori.x, -ori.y,
                          -ori.z);
    }
  }
}

int SceneGetTwoSidedLightingSettings(PyMOLGlobals* G, const CSetting* set1,
                                     const CSetting* set2) {
  int two_sided_lighting =
      SettingGet_b(G, set1, set2, cSetting_two_sided_lighting);
  if (two_sided_lighting < 0) {
    if (SettingGet_i(G, set1, set2, cSetting_surface_cavity_mode))
      two_sided_lighting = true;
    else
      two_sided_lighting = false;
  }
  two_sided_lighting =
      two_sided_lighting ||
      (SettingGet_i(G, set1, set2, cSetting_transparency_mode) == 1);
  return two_sided_lighting;
}

float SceneGetDynamicLineWidth(RenderInfo* info, float line_width) {
  if (info && info->dynamic_width) {
    float factor;
    if (info->vertex_scale > R_SMALL4) {
      factor = info->dynamic_width_factor / info->vertex_scale;
      if (factor > info->dynamic_width_max) factor = info->dynamic_width_max;
      if (factor < info->dynamic_width_min) {
        factor = info->dynamic_width_min;
      }
    } else {
      factor = info->dynamic_width_max;
    }
    return factor * line_width;
  }
  return line_width;
}

float SceneGetLineWidthForCylinders(PyMOLGlobals* G, RenderInfo* info,
                                    float line_width_arg) {
  float line_width = SceneGetDynamicLineWidth(info, line_width_arg);

  float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);

  if (pixel_scale_value < 0) pixel_scale_value = 1.0F;
  /* the radius of the cylinders is the vertex_scale * ray_pixel_scale */
  /* this turns out to be exactly right, but changes if the scene or user
     moves */
  return info->vertex_scale * pixel_scale_value * line_width / 2.f;
}

// line width arg has been processed as dynamic already
float SceneGetLineWidthForCylindersStatic(PyMOLGlobals* G, RenderInfo* info,
                                          float dynamic_line_width_arg,
                                          float line_width_arg) {
  float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
  if (pixel_scale_value < 0) pixel_scale_value = 1.0F;
  /* the radius of the cylinders is the vertex_scale * ray_pixel_scale */
  /* this turns out to be exactly right, but changes if the scene or user
     moves */
  if (SceneGetStereo(G) == cStereo_openvr)
    // note: that is reversion of magic dynamic line modifications in PYMOL
    return pixel_scale_value * 0.07f * line_width_arg / 2.0f;

  return info->vertex_scale * pixel_scale_value * dynamic_line_width_arg / 2.f;
}

void ScenePushModelViewMatrix(PyMOLGlobals* G) {
  auto I = G->Scene;
  I->m_ModelViewMatrixStack.push_back(I->modelViewMatrix);
}

void ScenePopModelViewMatrix(PyMOLGlobals* G, bool immediate) {
  CScene* I = G->Scene;
  auto& stack = I->m_ModelViewMatrixStack;

  if (stack.empty()) {
    printf("ERROR: depth == 0\n");
    return;
  }

  I->modelViewMatrix = stack.back();
  stack.pop_back();

  // Metal backend: re-sync the renderer's cached modelview to the restored
  // matrix, so a per-object TTT applied by ObjectPrepareContext (which pushes
  // camera*TTT into the renderer, see PyMOLObject.cpp) does not leak into the
  // next object. The glLoadMatrixf below is a no-op on Metal.
  if (G->Renderer) {
    G->Renderer->matrixMode(0);  // modelview
    G->Renderer->loadMatrixf(SceneGetModelViewMatrixPtr(G));
  }
}

glm::mat4& SceneGetModelViewMatrix(PyMOLGlobals* G) {
  return G->Scene->modelViewMatrix;
}

float* SceneGetModelViewMatrixPtr(PyMOLGlobals* G) {
  auto& mat = SceneGetModelViewMatrix(G);
  return glm::value_ptr(mat);
}

glm::mat4& SceneGetProjectionMatrix(PyMOLGlobals* G) {
  return G->Scene->projectionMatrix;
}

float* SceneGetProjectionMatrixPtr(PyMOLGlobals* G) {
  auto& mat = SceneGetProjectionMatrix(G);
  return glm::value_ptr(mat);
}

void SceneGLClear(PyMOLGlobals* G, GLbitfield mask) { glClear(mask); }

int SceneIncrementTextureRefreshes(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return ++(I->n_texture_refreshes);
}

void SceneResetTextureRefreshes(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->n_texture_refreshes = 0;
}

void SceneGetScaledAxes(PyMOLGlobals* G, pymol::CObject* obj, float* xn,
                        float* yn) {
  CScene* I = G->Scene;
  float* v;
  float vt[3];
  float xn0[3] = {1.0F, 0.0F, 0.0F};
  float yn0[3] = {0.0F, 1.0F, 0.0F};
  float v_scale;

  v = TextGetPos(G);

  if (obj->TTTFlag) {
    transformTTT44f3f(obj->TTT, v, vt);
  } else {
    copy3f(v, vt);
  }

  v_scale = SceneGetScreenVertexScale(G, vt);

  MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), xn0,
                                xn0);
  MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), yn0,
                                yn0);
  scale3f(xn0, v_scale, xn);
  scale3f(yn0, v_scale, yn);
}

int SceneGetCopyType(PyMOLGlobals* G) { return G->Scene->CopyType; }

void SceneAdjustZtoScreenZ(PyMOLGlobals* G, float* pos, float zarg) {
  CScene* I = G->Scene;
  float clipRange =
      (I->m_view.m_clipSafe().m_back - I->m_view.m_clipSafe().m_front);
  float z = (zarg + 1.f) / 2.f;
  float zInPreProj = -(z * clipRange + I->m_view.m_clipSafe().m_front);
  float pos4[4], tpos[4], npos[4];
  float InvModMatrix[16];
  copy3f(pos, pos4);
  pos4[3] = 1.f;
  MatrixTransformC44f4f(SceneGetModelViewMatrixPtr(G), pos4, tpos);
  normalize4f(tpos);
  /* NEED TO ACCOUNT FOR ORTHO */
  if (SettingGetGlobal_b(G, cSetting_ortho)) {
    npos[0] = tpos[0];
    npos[1] = tpos[1];
  } else {
    npos[0] = zInPreProj * tpos[0] / tpos[2];
    npos[1] = zInPreProj * tpos[1] / tpos[2];
  }
  npos[2] = zInPreProj;
  npos[3] = 1.f;

  MatrixInvertC44f(SceneGetModelViewMatrixPtr(G), InvModMatrix);
  MatrixTransformC44f4f(InvModMatrix, npos, npos);
  normalize4f(npos);
  copy3f(npos, pos);
}

/* this function takes a screen point, where z is normalized between the
   clipping planes, and converts it to the world coordinates */
void SceneSetPointToWorldScreenRelative(PyMOLGlobals* G, float* pos,
                                        float* screenPt) {
  float npos[4];
  float InvPmvMatrix[16];
  auto extent = SceneGetExtentStereo(G);
  npos[0] = (.5f + floor(screenPt[0] * extent.width)) /
            extent.width;  // add .5, in middle of pixels?
  npos[1] = (.5f + floor(screenPt[1] * extent.height)) /
            extent.height;  // add .5, in middle of pixels?
  npos[2] = 0.f;
  npos[3] = 1.f;
  MatrixInvertC44f(SceneGetPmvMatrix(G), InvPmvMatrix);
  MatrixTransformC44f4f(InvPmvMatrix, npos, npos);
  normalize4f(npos);
  SceneAdjustZtoScreenZ(G, npos, screenPt[2]);
  copy3f(npos, pos);
}

float SceneGetCurrentBackSafe(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (I->m_view.m_clipSafe().m_back);
}
float SceneGetCurrentFrontSafe(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (I->m_view.m_clipSafe().m_front);
}

/**
 * Get the field-of-view width at a depth of 1.0
 */
float GetFovWidth(PyMOLGlobals* G) {
  float fov = SettingGetGlobal_f(G, cSetting_field_of_view);
  return 2.f * tanf(fov * PI / 360.f);
}

void SceneInvalidatePicking(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  I->pickmgr.invalidate();
}

SceneElem::SceneElem(std::string name_, bool drawn_)
    : name(std::move(name_)), drawn(drawn_) {}

void SceneSetViewport(PyMOLGlobals* G, const Rect2D& rect) {
  switch (G->GFXMgr->backend()) {
    case GFXAPIBackend::OPENGL:
      glViewport(rect.offset.x, rect.offset.y, rect.extent.width,
                 rect.extent.height);
    default:
      break;
  }
}

void SceneSetViewport(PyMOLGlobals* G, int x, int y, int width, int height) {
  assert(width >= 0 && height >= 0);
  SceneSetViewport(G, Rect2D{x, y, static_cast<std::uint32_t>(width),
                             static_cast<std::uint32_t>(height)});
}

Rect2D SceneGetViewport(PyMOLGlobals* G) {
  Rect2D viewport{};
  int viewBuffer[4];
  glGetIntegerv(GL_VIEWPORT, (GLint*)(void*)viewBuffer);
  viewport.offset = Offset2D{static_cast<std::int32_t>(viewBuffer[0]),
                             static_cast<std::int32_t>(viewBuffer[1])};
  viewport.extent = Extent2D{static_cast<std::uint32_t>(viewBuffer[2]),
                             static_cast<std::uint32_t>(viewBuffer[3])};
  return viewport;
}
