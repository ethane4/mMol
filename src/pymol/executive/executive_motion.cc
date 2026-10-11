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
   -*
   -*
   -*
   Z* -------------------------------------------------------------------
*/

// Camera and object motions for movies.

#include "pymol/executive/executive.h"
#include "pymol/base/list_macros.h"
#include "pymol/ui/menu.h"
#include "pymol/scene/movie.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/executive/executive_internal.h"

int ExecutiveGroupMotionModify(PyMOLGlobals* G, pymol::CObject* group,
                               ViewElemAction action, int index, int count,
                               int target, int freeze) {
  CExecutive* I = G->Executive;
  int result = true;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G, group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectGroup) {
            ObjectMotionModify(rec->obj, action, index, count, target, freeze,
                               true);
          }
          break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupMotion(PyMOLGlobals* G, pymol::CObject* group,
                         MViewAction action, int first, int last, float power,
                         float bias, int simple, float linear, int wrap,
                         int hand, int window, int cycles, int state,
                         int quiet) {
  CExecutive* I = G->Executive;
  int result = true;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G, group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectGroup) {
            ObjectMotion(rec->obj, action, first, last, power, bias, simple,
                         linear, wrap, hand, window, cycles, state, quiet);
          }
          break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupCombineTTT(PyMOLGlobals* G, pymol::CObject* group,
                             const float* ttt, int reverse_order, int store) {
  CExecutive* I = G->Executive;
  int result = true;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G, group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectGroup) {
            ObjectCombineTTT(rec->obj, ttt, reverse_order, store);
          }
          break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveGroupTranslateTTT(PyMOLGlobals* G, pymol::CObject* group,
                               const float* v, int store) {
  CExecutive* I = G->Executive;
  int result = true;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetExpandedGroupList(G, group->Name);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec* rec;
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecObject:
          if (rec->obj->type != cObjectGroup) {
            ObjectTranslateTTT(rec->obj, v, store);
          }
          break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveMotionView(PyMOLGlobals* G, MViewAction action, int first,
                        int last, float power, float bias, int simple,
                        float linear, const char* name, int wrap, int hand,
                        int window, int cycles, const char* scene_name,
                        float scene_cut, int state, int quiet, int autogen) {
  int ok = true;

  CExecutive* I = G->Executive;

  if (wrap < 0) {
    wrap = SettingGetGlobal_b(G, cSetting_movie_loop);
  }

  if ((!name) || (!name[0]) || (!strcmp(name, cKeywordNone)) ||
      (!strcmp(name, cKeywordAll)) || (!strcmp(name, cKeywordSame))) {
    if (autogen) {
      /* if autogenerated, then use current settings */
      power = SettingGetGlobal_f(G, cSetting_motion_power);
      bias = SettingGetGlobal_f(G, cSetting_motion_bias);
      linear = SettingGetGlobal_f(G, cSetting_motion_linear);
      hand = SettingGetGlobal_i(G, cSetting_motion_hand);
    }
    /* camera */

    ok = MovieView(G, action, first, last, power, bias, true, linear, wrap,
                   hand, window, cycles, scene_name, scene_cut, state, quiet);

    if (name && name[0] && strcmp(name, cKeywordNone)) {
      /* also do all other objects */
      SpecRec* rec = nullptr;
      while (ListIterate(I->Spec, rec, next)) {
        switch (rec->type) {
          case cExecObject:
            if (autogen) {
              power = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                   cSetting_motion_power);
              bias = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                  cSetting_motion_bias);
              simple = SettingGet_i(G, nullptr, rec->obj->Setting.get(),
                                    cSetting_motion_simple);
              linear = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                    cSetting_motion_linear);
              hand = SettingGet_i(G, nullptr, rec->obj->Setting.get(),
                                  cSetting_motion_hand);
            }
            if ((ObjectGetSpecLevel(rec->obj, 0) >= 0) ||
                (!strcmp(name, cKeywordAll))) {
              ok = ObjectMotion(rec->obj, action, first, last, power, bias,
                                simple < 0 ? 0 : 1, linear, wrap, hand, window,
                                cycles, state, quiet);
            }
            break;
        }
      }
    }
  } else { /* pattern */
    CTracker* I_Tracker = I->Tracker;
    SpecRec* rec = nullptr;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecObject:
            if (autogen) {
              power = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                   cSetting_motion_power);
              bias = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                  cSetting_motion_bias);
              simple = SettingGet_i(G, nullptr, rec->obj->Setting.get(),
                                    cSetting_motion_simple);
              linear = SettingGet_f(G, nullptr, rec->obj->Setting.get(),
                                    cSetting_motion_linear);
              hand = SettingGet_i(G, nullptr, rec->obj->Setting.get(),
                                  cSetting_motion_hand);
            }

            ok = ObjectMotion(rec->obj, action, first, last, power, bias,
                              simple < 0 ? 0 : simple, linear, wrap, hand,
                              window, cycles, state, quiet);
            break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    // fix for PYMOL-1465
    OrthoReshape(G, -1, -1, false);
  }
  ExecutiveCountMotions(G);
  return ok;
}

pymol::Result<> ExecutiveMotionViewModify(PyMOLGlobals* G,
                                          ViewElemAction action, int index,
                                          int count, int target,
                                          const char* name, int freeze,
                                          int quiet) {
  CExecutive* I = G->Executive;
  if ((!name) || (!name[0]) || (!strcmp(name, cKeywordNone)) ||
      (!strcmp(name, cKeywordSame)) || (!strcmp(name, cKeywordAll))) {
    /* camera */
    if (MovieGetSpecLevel(G, 0) >= 0) {
      MovieViewModify(G, action, index, count, target, true, true);
    }
    if ((!name) || strcmp(name, cKeywordNone)) {
      /* also do all other objects */
      SpecRec* rec = nullptr;
      while (ListIterate(I->Spec, rec, next)) {
        switch (rec->type) {
          case cExecObject:
            if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
              ObjectMotionModify(rec->obj, action, index, count, target, true,
                                 true);
            }
            break;
        }
      }
      ExecutiveMotionTrim(G);
    } else {
      ExecutiveMotionExtend(G, true);
    }

    if ((!freeze) && SettingGetGlobal_i(G, cSetting_movie_auto_interpolate)) {
      ExecutiveMotionReinterpolate(G);
    }
  } else { /* pattern */
    CTracker* I_Tracker = I->Tracker;
    SpecRec* rec = nullptr;
    int list_id = ExecutiveGetNamesListFromPattern(G, name, true, true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec) {
        switch (rec->type) {
          case cExecObject:
            if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
              /* only modify objects with motion matrices */
              ObjectMotionModify(rec->obj, action, index, count, target, freeze,
                                 false);
            }
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  ExecutiveCountMotions(G);
  SceneCountFrames(G);
  return {};
}

void ExecutiveMotionReinterpolate(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) {
    switch (rec->type) {
      case cExecAll:
        if (MovieGetSpecLevel(G, 0) >= 0) {
          MovieViewReinterpolate(G);
        }
        break;
      case cExecObject:
        if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
          ObjectMotionReinterpolate(rec->obj);
        }
        break;
    }
  }
}

void ExecutiveMotionTrim(PyMOLGlobals* G) {
  int n_frame = MovieGetLength(G);
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) {
    switch (rec->type) {
      case cExecObject:
        if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
          ObjectMotionTrim(rec->obj, n_frame);
        }
        break;
    }
  }
}

void ExecutiveMotionExtend(PyMOLGlobals* G, int freeze) {
  int n_frame = 0;
  int max_length = 0;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  if (MovieGetSpecLevel(G, -1) > 0) n_frame = MovieGetLength(G);
  while (ListIterate(I->Spec, rec, next)) {
    switch (rec->type) {
      case cExecObject:
        if (ObjectGetSpecLevel(rec->obj, -1) > 0) {
          int length = ObjectMotionGetLength(rec->obj);
          if (max_length < length) max_length = length;
        }
        break;
    }
  }
  if (max_length) {
    if (n_frame < max_length) MovieViewTrim(G, max_length);
    while (ListIterate(I->Spec, rec, next)) {
      switch (rec->type) {
        case cExecObject:
          if (ObjectGetSpecLevel(rec->obj, -1) > 0) {
            ObjectMotionTrim(rec->obj, max_length);
          }
          break;
      }
    }
  }
  if ((!freeze) && SettingGetGlobal_i(G, cSetting_movie_auto_interpolate)) {
    ExecutiveMotionReinterpolate(G);
  }
}

int ExecutiveCountMotions(PyMOLGlobals* G) {
  int count = 0;
  CExecutive* I = G->Executive;
  if (MovieGetLength(G)) {
    SpecRec* rec = nullptr;
    while (ListIterate(I->Spec, rec, next)) {
      switch (rec->type) {
        case cExecAll:
          if (MovieGetSpecLevel(G, 0) >= 0) count++;
          break;
        case cExecObject:
          if (ObjectGetSpecLevel(rec->obj, 0) >= 0) count++;
          break;
      }
    }
  }

  if (count < 1 && SceneGetNFrame(G) > 1) count = 1;

  if (count != I->LastMotionCount) {
    if (SettingGetGlobal_i(G, cSetting_movie_panel)) {
      OrthoDoViewportWhenReleased(G);
    }
    I->LastMotionCount = count;
  }

  return (count);
}

void ExecutiveMotionDraw(PyMOLGlobals* G, BlockRect* rect, int expected,
                         CGO* orthoCGO) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  int frames = MovieGetLength(G);
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  while (ListIterate(I->Spec, rec, next)) {
    switch (rec->type) {
      case cExecAll:
        if (MovieGetSpecLevel(G, 0) >= 0) {
          int presentation = SettingGetGlobal_b(G, cSetting_presentation);
          if (presentation) {
            expected = 1;
          }
          draw_rect.top = rect->top - (height * count) / expected;
          draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
          MovieDrawViewElem(G, &draw_rect, frames, orthoCGO);
          count++;
          if (presentation) {
            goto done;
          }
        }
        break;
      case cExecObject:
        if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
          draw_rect.top = rect->top - (height * count) / expected;
          draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
          ObjectDrawViewElem(rec->obj, &draw_rect, frames, orthoCGO);
          count++;
        }
        break;
    }
  }
done:
  return;
}

void ExecutiveMotionMenuActivate(PyMOLGlobals* G, BlockRect* rect, int expected,
                                 int passive, int x, int y, int same) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  if (same) {
    if (MovieGetSpecLevel(G, 0) >= 0) {
      int n_frame = MovieGetLength(G);
      int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
      WordType frame_str = "0";
      if ((frame >= 0) && (frame < n_frame)) {
        sprintf(frame_str, "%d", frame + 1);
      }
      MenuActivate2Arg(G, x, y, x, y, passive, "obj_motion", cKeywordSame,
                       frame_str);
    }
  } else {
    while (ListIterate(I->Spec, rec, next)) {
      switch (rec->type) {
        case cExecAll:
          if (MovieGetSpecLevel(G, 0) >= 0) {
            draw_rect.top = rect->top - (height * count) / expected;
            draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
            if ((y > draw_rect.bottom) && (y < draw_rect.top)) {
              int n_frame = MovieGetLength(G);
              int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
              WordType frame_str = "0";
              if ((frame >= 0) && (frame < n_frame)) {
                sprintf(frame_str, "%d", frame + 1);
              }
              MenuActivate1Arg(G, x, y, x, y, passive, "camera_motion",
                               frame_str);
              goto done;
            }
            count++;
          }
          break;
        case cExecObject:
          if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
            draw_rect.top = rect->top - (height * count) / expected;
            draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
            if ((y > draw_rect.bottom) && (y < draw_rect.top)) {
              int n_frame = MovieGetLength(G);
              int frame = MovieXtoFrame(G, &draw_rect, n_frame, x, false);
              WordType frame_str = "0";
              if ((frame >= 0) && (frame < n_frame)) {
                sprintf(frame_str, "%d", frame + 1);
              }
              MenuActivate2Arg(G, x, y, x, y, passive, "obj_motion",
                               rec->obj->Name, frame_str);
              goto done;
            }
            count++;
          }
          break;
      }
    }
  }
done:
  return;
}

void ExecutiveMotionClick(PyMOLGlobals* G, BlockRect* rect, int mode,
                          int expected, int x, int y, int nearest) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  BlockRect draw_rect = *rect;
  int count = 0;
  int height = rect->top - rect->bottom;
  while (ListIterate(I->Spec, rec, next)) {
    switch (rec->type) {
      case cExecAll:
        if (MovieGetSpecLevel(G, 0) >= 0) {
          draw_rect.top = rect->top - (height * count) / expected;
          draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
          if ((y >= draw_rect.bottom) && (y <= draw_rect.top)) {
            MoviePrepareDrag(G, &draw_rect, NULL, mode, x, y, nearest);
            goto done;
          }
          count++;
        }
        break;
      case cExecObject:
        if (ObjectGetSpecLevel(rec->obj, 0) >= 0) {
          MoviePrepareDrag(G, rect, NULL, mode, x, y, nearest);
          draw_rect.top = rect->top - (height * count) / expected;
          draw_rect.bottom = rect->top - (height * (count + 1)) / expected;
          if ((y >= draw_rect.bottom) && (y <= draw_rect.top)) {
            MoviePrepareDrag(G, &draw_rect, rec->obj, mode, x, y, nearest);
            goto done;
          }
          count++;
        }
        break;
    }
  }
done:
  return;
}
