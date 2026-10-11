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

// Moving the camera: views, translation, rotation, zoom, clipping and the
// origin.

#include "pymol/geometry/matrix.h"
#include "pymol/scene/scene.h"
#include "pymol/ui/but_mode.h"
#include "pymol/scene/movie.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/scene/scene_internal.h"

#define SceneRotateWithDirty SceneRotate
static void SceneClipSetWithDirty(PyMOLGlobals* G, float front, float back,
                                  int dirty);

void SceneToViewElem(PyMOLGlobals* G, CViewElem* elem, const char* scene_name) {
  double* dp;
  CScene* I = G->Scene;
  const auto& pos = I->m_view.pos();
  const auto& ori = I->m_view.origin();

  float dY = 0, dZ = 0;
  float fov = SettingGetGlobal_f(G, cSetting_field_of_view);
  float scale = 1.0f / I->Scale;

  /* copy rotation matrix */
  elem->matrix_flag = true;
  dp = elem->matrix;
  auto fp = glm::value_ptr(I->m_view.rotMatrix());
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);

  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);

  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);
  *(dp++) = (double)*(fp++);

  *(dp++) = 0.0F;
  *(dp++) = 0.0F;
  *(dp++) = 0.0F;
  *(dp++) = 1.0F;

  /* copy position */
  elem->pre_flag = true;
  dp = elem->pre;
  *(dp++) = (double)pos.x * scale;
  *(dp++) = (double)pos.y * scale - dY;
  *(dp++) = (double)pos.z * scale - dZ;

  /* copy origin (negative) */
  elem->post_flag = true;
  dp = elem->post;
  *(dp++) = (double)-ori.x;
  *(dp++) = (double)-ori.y;
  *(dp++) = (double)-ori.z;

  elem->clip_flag = true;
  elem->front = I->m_view.m_clip().m_front * scale + dZ;
  elem->back = I->m_view.m_clip().m_back * scale + dZ;

  elem->ortho_flag = true;
  elem->ortho = SettingGetGlobal_b(G, cSetting_ortho) ? fov : -fov;

  {
    if (elem->scene_flag && elem->scene_name) {
      OVLexicon_DecRef(G->Lexicon, elem->scene_name);
      elem->scene_name = 0;
      elem->scene_flag = 0;
    }
  }
  {
    if (!scene_name)
      scene_name = SettingGetGlobal_s(G, cSetting_scene_current_name);
    if (scene_name && scene_name[0]) {
      OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon, scene_name);
      if (OVreturn_IS_OK(result)) {
        elem->scene_name = result.word;
        elem->scene_flag = true;
      }
    }
  }
}

void SceneFromViewElem(PyMOLGlobals* G, CViewElem* elem, int dirty) {
  CScene* I = G->Scene;
  float* fp;
  double* dp;
  int changed_flag = false;

  float dY = 0, dZ = 0;
  float fov = elem->ortho;
  float scale = I->Scale;
  auto pos = I->m_view.pos();
  auto ori = I->m_view.origin();
  auto rot = I->m_view.rotMatrix();

  if (elem->matrix_flag) {
    rot = glm::make_mat4(elem->matrix);
    changed_flag = true;
    SceneUpdateInvMatrix(G);
  }

  if (elem->pre_flag) {
    dp = elem->pre;
    fp = glm::value_ptr(pos);
    *(fp++) = (float)*(dp++) * scale;
    *(fp++) = (float)*(dp++) * scale - dY;
    *(fp++) = (float)*(dp++) * scale - dZ;
    changed_flag = true;
  }

  if (elem->post_flag) {
    dp = elem->post;
    fp = glm::value_ptr(ori);
    *(fp++) = (float)(-*(dp++));
    *(fp++) = (float)(-*(dp++));
    *(fp++) = (float)(-*(dp++));
    changed_flag = true;
  }

  if (elem->clip_flag) {
    SceneClipSetWithDirty(G, elem->front * scale + dZ, elem->back * scale + dZ,
                          dirty);
  }

  if (elem->ortho_flag) {
    if (fov < 0.0F) {
      SettingSetGlobal_b(G, cSetting_ortho, 0);
      if (fov < -(1.0F - R_SMALL4)) {
        SettingSetGlobal_f(G, cSetting_field_of_view, -fov);
      }
    } else {
      SettingSetGlobal_b(G, cSetting_ortho, (fov > 0.5F));
      if (fov > (1.0F + R_SMALL4)) {
        SettingSetGlobal_f(G, cSetting_field_of_view, fov);
      }
    }
  }
  if (elem->state_flag && !MovieDefined(G)) {
    SettingSetGlobal_i(G, cSetting_state, (elem->state) + 1);
  }
  if (changed_flag) {
    SceneRestartSweepTimer(G);
    I->RockFrame = 0;
    SceneRovingDirty(G);
    I->m_view.setPos(pos);
    I->m_view.setOrigin(ori);
    I->m_view.setRotMatrix(rot);
  }
}
/**
   Get information required to define the geometry
   of a particular view, for shipping to and from python
   as a list of floats
   @verbatim
   0-15 = 4x4 rotation matrix
   16-18 = position
   19-21 = origin
   22    = front plane
   23    = rear plane
   24    = orthoscopic flag
   @endverbatim
   @param[out] view buffer to fill
*/
void SceneGetView(PyMOLGlobals* G, SceneViewType view) {
  float* p;
  CScene* I = G->Scene;
  p = view;

  float dY = 0, dZ = 0;
  float fov = SettingGetGlobal_f(G, cSetting_field_of_view);
  float scale = 1.0f / I->Scale;
  const auto& pos = I->m_view.pos();
  const auto& ori = I->m_view.origin();

  std::copy_n(glm::value_ptr(I->m_view.rotMatrix()), 16, p);
  p += 16;
  *(p++) = pos.x * scale;
  *(p++) = pos.y * scale - dY;
  *(p++) = pos.z * scale - dZ;
  *(p++) = ori.x;
  *(p++) = ori.y;
  *(p++) = ori.z;
  *(p++) = I->m_view.m_clip().m_front * scale + dZ;
  *(p++) = I->m_view.m_clip().m_back * scale + dZ;
  *(p++) = SettingGetGlobal_b(G, cSetting_ortho) ? fov : -fov;
}
void SceneSetView(PyMOLGlobals* G, const SceneViewType view, int quiet,
                  float animate, int hand) {
  const float* p;
  CScene* I = G->Scene;

  if (animate < 0.0F) {
    if (SettingGetGlobal_b(G, cSetting_animation))
      animate = SettingGetGlobal_f(G, cSetting_animation_duration);
    else
      animate = 0.0F;
  }
  if (animate != 0.0F)
    ScenePrimeAnimation(G);
  else {
    SceneAbortAnimation(G);
  }

  float dY = 0, dZ = 0;
  float fov = view[24];
  float scale = I->Scale;

  p = view;
  I->m_view.setRotMatrix(glm::make_mat4(p));
  p += 16;
  SceneUpdateInvMatrix(G);
  auto newX = *(p++) * scale;
  auto newY = *(p++) * scale - dY;
  auto newZ = *(p++) * scale - dZ;
  I->m_view.setPos(newX, newY, newZ);
  auto newOriX = *(p++);
  auto newOriY = *(p++);
  auto newOriZ = *(p++);
  I->m_view.setOrigin(newOriX, newOriY, newOriZ);

  I->LastSweep = 0.0F;
  I->LastSweepX = 0.0F;
  I->LastSweepY = 0.0F;
  I->SweepTime = 0.0;
  I->RockFrame = 0;

  SceneClipSet(G, p[0] * scale + dZ, p[1] * scale + dZ);

  p += 2;
  if (fov < 0.0F) {
    SettingSetGlobal_b(G, cSetting_ortho, 0);
    if (fov < -(1.0F - R_SMALL4)) {
      SettingSetGlobal_f(G, cSetting_field_of_view, -fov);
    }
  } else {
    SettingSetGlobal_b(G, cSetting_ortho, (fov > 0.5F));
    if (fov > (1.0F + R_SMALL4)) {
      SettingSetGlobal_f(G, cSetting_field_of_view, fov);
    }
  }
  if (!quiet) {
    PRINTFB(G, FB_Scene, FB_Actions)
    " Scene: view updated.\n" ENDFB(G);
  }

  if (animate != 0.0F) SceneLoadAnimation(G, animate, hand);

  SceneRovingDirty(G);
}
void SceneTranslate(PyMOLGlobals* G, float x, float y, float z) {
  CScene* I = G->Scene;
  I->m_view.translate(x, y, z);
  SceneClipSet(G, I->m_view.m_clip().m_front - z,
               I->m_view.m_clip().m_back - z);
}

void SceneTranslateScaled(PyMOLGlobals* G, float x, float y, float z,
                          int sdof_mode) {
  CScene* I = G->Scene;
  int invalidate = false;

  switch (sdof_mode) {
    case SDOF_NORMAL_MODE:
      if ((x != 0.0F) || (y != 0.0F)) {
        float vScale = SceneGetExactScreenVertexScale(G, nullptr);
        float factor = vScale * (I->Height + I->Width) / 2;
        I->m_view.translate(x * factor, y * factor, 0.0f);
        invalidate = true;
      }
      if (z != 0.0F) {
        float factor =
            ((I->m_view.m_clipSafe().m_front + I->m_view.m_clipSafe().m_back) /
             2); /* average distance within visible space */
        if (factor > 0.0F) {
          factor *= z;
          I->m_view.translate(0.0f, 0.0f, factor);
          I->m_view.m_clip().m_front -= factor;
          I->m_view.m_clip().m_back -= factor;
          UpdateFrontBackSafe(I);
        }
        invalidate = true;
      }
      break;
    case SDOF_CLIP_MODE:
      if ((x != 0.0F) || (y != 0.0F)) {
        float vScale = SceneGetExactScreenVertexScale(G, nullptr);
        float factor = vScale * (I->Height + I->Width) / 2;
        I->m_view.translate(x * factor, y * factor, 0.0f);
        invalidate = true;
      }
      if (z != 0.0F) {
        float factor =
            ((I->m_view.m_clipSafe().m_front + I->m_view.m_clipSafe().m_back) /
             2); /* average distance within visible space */
        if (factor > 0.0F) {
          factor *= z;
          {
            float old_front = I->m_view.m_clip().m_front;
            float old_back = I->m_view.m_clip().m_back;
            float old_origin = -I->m_view.pos().z;
            SceneClip(G, SceneClipMode::Linear, factor, nullptr, 0);
            SceneDoRoving(G, old_front, old_back, old_origin, true, true);
          }
          invalidate = true;
        }
      }
      break;
    case SDOF_DRAG_MODE: {
      float v2[3];
      float scale = SettingGetGlobal_f(G, cSetting_sdof_drag_scale);

      {
        /* when dragging, we treat all axes proportionately */
        float vScale = SceneGetExactScreenVertexScale(G, nullptr);
        float factor = vScale * (I->Height + I->Width) / 2;
        x *= factor;
        y *= factor;
        z *= factor;
      }

      v2[0] = x * scale;
      v2[1] = y * scale;
      v2[2] = z * scale;

      /* transform into model coodinate space */
      MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v2,
                                    v2);

      EditorDrag(G, nullptr, -1, cButModeMovDrag,
                 SettingGetGlobal_i(G, cSetting_state) - 1, nullptr, v2,
                 nullptr);
    } break;
  }
  if (invalidate) {
    SceneInvalidate(G);
    if (SettingGetGlobal_b(G, cSetting_roving_origin)) {
      float v2[3];
      SceneGetCenter(G, v2); /* gets position of center of screen */
      SceneOriginSet(G, v2, true);
    }
    if (SettingGetGlobal_b(G, cSetting_roving_detail)) {
      SceneRovingPostpone(G);
    }
  }
}

void SceneRotateScaled(PyMOLGlobals* G, float rx, float ry, float rz,
                       int sdof_mode) {
  CScene* I = G->Scene;
  int invalidate = false;
  float axis[3];
  switch (sdof_mode) {
    case SDOF_NORMAL_MODE:
      axis[0] = rx;
      axis[1] = ry;
      axis[2] = rz;
      {
        float angle = length3f(axis);
        normalize3f(axis);
        SceneRotate(G, 60 * angle, axis[0], axis[1], axis[2]);
      }
      break;
    case SDOF_CLIP_MODE:
      if ((fabs(rz) > fabs(rx)) || (fabs(rz) > fabs(rx))) {
        rx = 0.0;
        ry = 0.0;
      } else {
        rz = 0.0;
      }
      axis[0] = rx;
      axis[1] = ry;
      axis[2] = 0.0;
      {
        float angle = length3f(axis);
        normalize3f(axis);
        SceneRotate(G, 60 * angle, axis[0], axis[1], axis[2]);
      }
      if (axis[2] != rz) {
        SceneClip(G, SceneClipMode::Scaling, 1.0F + rz, nullptr, 0);
      }
      break;
    case SDOF_DRAG_MODE: {
      float scale = SettingGetGlobal_f(G, cSetting_sdof_drag_scale);
      float v1[3], v2[3];
      axis[0] = rx;
      axis[1] = ry;
      axis[2] = rz;

      EditorReadyDrag(G, SettingGetGlobal_i(G, cSetting_state) - 1);

      {
        float angle = length3f(axis);
        normalize3f(axis);

        v1[0] = cPI * (60 * angle / 180.0F) * scale;

        /* transform into model coodinate space */
        MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()),
                                      axis, v2);

        EditorDrag(G, nullptr, -1, cButModeRotDrag,
                   SettingGetGlobal_i(G, cSetting_state) - 1, v1, v2, nullptr);
        invalidate = true;
      }
    } break;
  }
  if (invalidate) {
    SceneInvalidate(G);
  }
}
static void SceneClipSetWithDirty(PyMOLGlobals* G, float front, float back,
                                  int dirty) {
  CScene* I = G->Scene;

  // minimum slab
  float minSlab = cSliceMin * I->Scale;
  if (back - front < minSlab) {
    float avg = (back + front) / 2.0;
    back = avg + minSlab / 2.0;
    front = avg - minSlab / 2.0;
  }

  I->m_view.m_clip().m_front = front;
  I->m_view.m_clip().m_back = back;

  UpdateFrontBackSafe(I);

  if (dirty)
    SceneInvalidate(G);
  else
    SceneInvalidateCopy(G, false);
}
void SceneClipSet(PyMOLGlobals* G, float front, float back) {
  SceneClipSetWithDirty(G, front, back, true);
}
static SceneClipMode SceneClipGetEnum(pymol::zstring_view mode) {
  static const std::unordered_map<pymol::zstring_view, SceneClipMode> modes{
      {"near", SceneClipMode::Near},
      {"far", SceneClipMode::Far},
      {"move", SceneClipMode::Move},
      {"slab", SceneClipMode::Slab},
      {"atoms", SceneClipMode::Atoms},
      {"scaling", SceneClipMode::Scaling},
      {"linear", SceneClipMode::Linear},
      {"near_set", SceneClipMode::Near_Set},
      {"far_set", SceneClipMode::Far_Set},
  };

  auto it = modes.find(mode);
  return (it == modes.end()) ? SceneClipMode::Invalid : it->second;
}

pymol::Result<> SceneClipFromMode(PyMOLGlobals* G, pymol::zstring_view mode,
                                  float movement, pymol::zstring_view sele,
                                  int state) {
  auto plane = SceneClipGetEnum(mode);
  if (plane == SceneClipMode::Invalid) {
    return pymol::make_error("invalid clip mode");
  }
  SceneClip(G, plane, movement, sele.c_str(), state);
  return {};
}
void SceneClip(PyMOLGlobals* G, SceneClipMode mode, float movement,
               const char* sele, int state) { /* 0=front, 1=back */
  CScene* I = G->Scene;
  float avg;
  float mn[3], mx[3], cent[3], v0[3], offset[3], origin[3];
  const auto& pos = I->m_view.pos();
  switch (mode) {
    case SceneClipMode::Near:
      SceneClipSet(G, I->m_view.m_clip().m_front - movement,
                   I->m_view.m_clip().m_back);
      break;
    case SceneClipMode::Far:
      SceneClipSet(G, I->m_view.m_clip().m_front,
                   I->m_view.m_clip().m_back - movement);
      break;
    case SceneClipMode::Move:
      SceneClipSet(G, I->m_view.m_clip().m_front - movement,
                   I->m_view.m_clip().m_back - movement);
      break;
    case SceneClipMode::Slab:
      if (sele[0]) {
        if (!ExecutiveGetExtent(G, sele, mn, mx, true, state, false))
          sele = nullptr;
        else {
          average3f(mn, mx, cent); /* get center of selection */
          subtract3f(cent, glm::value_ptr(I->m_view.origin()),
                     v0); /* how far from origin? */
          MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v0,
                                     offset); /* convert to view-space */
        }
      } else {
        sele = nullptr;
      }
      avg = (I->m_view.m_clip().m_front + I->m_view.m_clip().m_back) / 2.0F;
      movement /= 2.0F;
      if (sele) {
        avg = -pos.z - offset[2];
      }
      SceneClipSet(G, avg - movement, avg + movement);
      break;
    case SceneClipMode::Atoms:
      if (!sele)
        sele = cKeywordAll;
      else if (!sele[0]) {
        sele = cKeywordAll;
      }
      if (WordMatchExact(G, sele, cKeywordCenter, true)) {
        MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()),
                                   glm::value_ptr(I->m_view.origin()),
                                   origin); /* convert to view-space */
        SceneClipSet(G, origin[2] - movement, origin[2] + movement);
      } else if (WordMatchExact(G, sele, cKeywordOrigin, true)) {
        SceneClipSet(G, -pos.z - movement, -pos.z + movement);
      } else {
        if (!ExecutiveGetCameraExtent(G, sele, mn, mx, true, state))
          sele = nullptr;
        if (sele) {
          if (sele[0]) {
            average3f(mn, mx, cent); /* get center of selection */
            MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()),
                                       glm::value_ptr(I->m_view.origin()),
                                       origin); /* convert to view-space */
            subtract3f(mx, origin, mx);         /* how far from origin? */
            subtract3f(mn, origin, mn);         /* how far from origin? */
            SceneClipSet(G, -pos.z - mx[2] - movement,
                         -pos.z - mn[2] + movement);
          } else {
            sele = nullptr;
          }
        }
      }
      break;
    case SceneClipMode::Scaling: {
      double avg = (I->m_view.m_clip().m_front / 2.0) +
                   (I->m_view.m_clip().m_back / 2.0);
      double width_half = I->m_view.m_clip().m_back - avg;
      double new_w_half = std::min(
          movement * width_half,
          width_half + 1000.0);  // prevent exploding of clipping planes

      SceneClipSet(G, avg - new_w_half, avg + new_w_half);
    } break;
    case SceneClipMode::Proportional: {
      float shift =
          (I->m_view.m_clip().m_front - I->m_view.m_clip().m_back) * movement;
      SceneClipSet(G, I->m_view.m_clip().m_front + shift,
                   I->m_view.m_clip().m_back + shift);
    } break;
    case SceneClipMode::Linear: {
      SceneClipSet(G, I->m_view.m_clip().m_front + movement,
                   I->m_view.m_clip().m_back + movement);
    } break;
    case SceneClipMode::Near_Set: {
      SceneClipSet(G, movement, I->m_view.m_clip().m_back);
    } break;
    case SceneClipMode::Far_Set: {
      SceneClipSet(G, I->m_view.m_clip().m_front, movement);
    } break;
  }
}

/**
 * Retrieves clipping plane distances from camera
 * @return near and far clipping plane distances
 * @TODO: Needs camera param index when multiple cameras are implemented
 * @TODO: Return error when camera idx found.
 */
pymol::Result<std::pair<float, float>> SceneGetClip(PyMOLGlobals* G) {
  auto clip = G->Scene->getSceneView().m_clip;
  return std::make_pair(clip.m_front, clip.m_back);
}
void SceneUpdateCameraRock(PyMOLGlobals* G, int dirty) {
  CScene* I = G->Scene;
  float ang_cur, disp, diff;
  float sweep_angle = SettingGetGlobal_f(G, cSetting_sweep_angle);
  float sweep_speed = SettingGetGlobal_f(G, cSetting_sweep_speed);
  float sweep_phase = SettingGetGlobal_f(G, cSetting_sweep_phase);
  int sweep_mode = SettingGetGlobal_i(G, cSetting_sweep_mode);
  float shift = (float)(PI / 2.0F);

  I->SweepTime += I->RenderTime;
  I->LastSweepTime = UtilGetSeconds(G);

  switch (sweep_mode) {
    case 0:
    case 1:
    case 2:
      if (sweep_angle <= 0.0F) {
        diff =
            (float)((PI / 180.0F) * I->RenderTime * 10 * sweep_speed / 0.75F);
      } else {
        ang_cur = (float)(I->SweepTime * sweep_speed) + sweep_phase;
        disp = (float)(sweep_angle * (PI / 180.0F) * sin(ang_cur) / 2);
        diff = (float)(disp - I->LastSweep);
        I->LastSweep = disp;
      }
      switch (sweep_mode) {
        case 0:
          SceneRotateWithDirty(G, (float)(180 * diff / PI), 0.0F, 1.0F, 0.0F,
                               dirty);
          break;
        case 1:
          SceneRotateWithDirty(G, (float)(180 * diff / PI), 1.0F, 0.0F, 0.0F,
                               dirty);
          break;
        case 2: /* z-rotation...useless! */
          SceneRotateWithDirty(G, (float)(180 * diff / PI), 0.0F, 0.0F, 1.0F,
                               dirty);
          break;
      }
      break;
    case 3: /* nutate */
      SceneRotateWithDirty(G, (float)(-I->LastSweepY), 0.0F, 1.0F, 0.0F, dirty);
      SceneRotateWithDirty(G, (float)(-I->LastSweepX), 1.0F, 0.0F, 0.0F, dirty);
      ang_cur = (float)(I->SweepTime * sweep_speed) + sweep_phase;

      I->LastSweepX = (float)(sweep_angle * sin(ang_cur) / 2);
      I->LastSweepY = (float)(sweep_angle * sin(ang_cur + shift) / 2);

      if (I->SweepTime * sweep_speed < PI) {
        float factor = (float)((I->SweepTime * sweep_speed) / PI);
        I->LastSweepX *= factor;
        I->LastSweepY *= factor;
      }
      SceneRotateWithDirty(G, (float)I->LastSweepX, 1.0F, 0.0F, 0.0F, dirty);
      SceneRotateWithDirty(G, (float)I->LastSweepY, 0.0F, 1.0F, 0.0F, dirty);
      break;
  }
}
/**
 * Zoom to location and radius
 */
void SceneWindowSphere(PyMOLGlobals* G, const float* location, float radius) {
  CScene* I = G->Scene;
  float v0[3];
  auto pos = I->m_view.pos();

  if (I->StereoMode == cStereo_openvr) {
    I->Scale = 1.0f / radius;
    radius = 1.0;
  } else {
    I->Scale = 1.0;
  }
  float dist = 2.f * radius / GetFovWidth(G);

  /* find where this point is in relationship to the origin */
  subtract3f(glm::value_ptr(I->m_view.origin()), location, v0);

  MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v0,
                             glm::value_ptr(pos)); /* convert to view-space */

  if (I->Height > I->Width && I->Height && I->Width)
    dist *= (float)I->Height / (float)I->Width;

  pos.z -= dist;
  I->m_view.m_clip().m_front = (-pos.z - radius * 1.2F);
  I->m_view.m_clip().m_back = (-pos.z + radius * 1.2F);
  UpdateFrontBackSafe(I);
  SceneRovingDirty(G);

  I->m_view.setPos(pos);
}
void SceneRelocate(PyMOLGlobals* G, const float* location) {
  CScene* I = G->Scene;
  float v0[3];
  auto pos = I->m_view.pos();

  auto slab_width = I->m_view.m_clip().m_back - I->m_view.m_clip().m_front;

  /* find out how far camera was from previous origin */
  auto dist = pos.z;

  // stay in front of camera, empirical value to show at least 1 bond
  if (dist > -5.f && I->StereoMode != cStereo_openvr) dist = -5.f;

  /* find where this point is in relationship to the origin */
  subtract3f(glm::value_ptr(I->m_view.origin()), location, v0);

  /*  printf("%8.3f %8.3f
   * %8.3f\n",I->m_view.m_clip().m_front,pos.z,I->m_view.m_clip().m_back); */

  auto pos_ptr = glm::value_ptr(pos);
  MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v0,
                             pos_ptr); /* convert to view-space */

  pos.z = dist;
  if (I->StereoMode == cStereo_openvr) {
    pos += glm::vec3(0.0f, 1.0f, 0.0f);
  }
  I->m_view.m_clip().m_front = (-pos.z - (slab_width * 0.50F));
  I->m_view.m_clip().m_back = (-pos.z + (slab_width * 0.50F));
  I->m_view.setPos(pos);
  UpdateFrontBackSafe(I);
  SceneRovingDirty(G);
}
/**
 * Get the origin of rotation in model space
 * cmd.get_view()[12:15]
 *
 * @param[out] origin
 */
void SceneOriginGet(PyMOLGlobals* G, float* origin) {
  CScene* I = G->Scene;
  copy3f(glm::value_ptr(I->m_view.origin()), origin);
}
/**
 * Set the origin of rotation in model space
 * (`cmd.get_view()[12:15]`)
 *
 * @param origin New origin
 * @param preserve preserve current viewing location
 */
void SceneOriginSet(PyMOLGlobals* G, const float* origin, int preserve) {
  CScene* I = G->Scene;
  float v0[3];
  glm::vec3 v1;

  if (preserve) { /* preserve current viewing location */
    subtract3f(origin, glm::value_ptr(I->m_view.origin()),
               v0); /* model-space translation */
    MatrixTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v0,
                               glm::value_ptr(v1)); /* convert to view-space */
    I->m_view.translate(v1); /* offset view to compensate */
  }
  I->m_view.setOrigin(origin[0], origin[1], origin[2]); /* move origin */
  SceneInvalidate(G);
}
/**
 * Update the scene rotation matrix (m_view.m_rotMatrix)
 *
 * @param angle Angle in degrees
 * @param x,y,z Axis
 * @param dirty Call SceneInvalidate()
 */
void SceneRotate(PyMOLGlobals* G, float angle, float x, float y, float z,
                 bool dirty) {
  CScene* I = G->Scene;
  float temp[16];
  angle = (float)(-PI * angle / 180.0);
  identity44f(temp);
  MatrixRotateC44f(temp, angle, x, y, z);
  MatrixMultiplyC44f(glm::value_ptr(I->m_view.rotMatrix()), temp);
  I->m_view.setRotMatrix(glm::make_mat4(temp));
  SceneUpdateInvMatrix(G);
  if (dirty) {
    SceneInvalidate(G);
  }
}

void SceneRotateAxis(PyMOLGlobals* G, float angle, char axis) {
  switch (axis) {
    case 'x':
      SceneRotate(G, angle, 1.0f, 0.0f, 0.0f);
      break;
    case 'y':
      SceneRotate(G, angle, 0.0f, 1.0f, 0.0f);
      break;
    case 'z':
      SceneRotate(G, angle, 0.0f, 0.0f, 1.0f);
      break;
  }
}
