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

// Roving detail: updating the representations near the camera center as it
// moves.

#include "pymol/geometry/matrix.h"
#include "pymol/scene/scene.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/scene/scene_internal.h"

void SceneRovingPostpone(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  float delay;
  if (SettingGetGlobal_b(G, cSetting_roving_detail)) {
    delay = SettingGetGlobal_f(G, cSetting_roving_delay);
    if (delay < 0.0F) {
      I->RovingLastUpdate = UtilGetSeconds(G); /* put off delay */
    }
  }
}

void SceneRovingDirty(PyMOLGlobals* G) {
  CScene* I = G->Scene;

  if (SettingGetGlobal_b(G, cSetting_roving_detail)) {
    SceneRovingPostpone(G);
    I->RovingDirtyFlag = true;
  }
}
void SceneDoRoving(PyMOLGlobals* G, float old_front, float old_back,
                   float old_origin, int adjust_flag, int zoom_flag) {
  EditorFavorOrigin(G, nullptr);
  if (SettingGetGlobal_b(G, cSetting_roving_origin)) {
    CScene* I = G->Scene;
    float delta_front, delta_back;
    float front_weight, back_weight, slab_width;
    float z_buffer = 3.0;
    float old_pos2 = 0.0F;
    float v2[3];

    z_buffer = SettingGetGlobal_f(G, cSetting_roving_origin_z_cushion);

    delta_front = I->m_view.m_clip().m_front - old_front;
    delta_back = I->m_view.m_clip().m_back - old_back;

    zero3f(v2);

    slab_width = I->m_view.m_clip().m_back - I->m_view.m_clip().m_front;

    /* first, check to make sure that the origin isn't too close to either plane
     */
    if ((z_buffer * 2) > slab_width) z_buffer = slab_width * 0.5F;

    if (old_origin < (I->m_view.m_clip().m_front +
                      z_buffer)) { /* old origin behind front plane */
      front_weight = 1.0F;
      delta_front = (I->m_view.m_clip().m_front + z_buffer) -
                    old_origin; /* move origin into allowed regioin */
    } else if (old_origin > (I->m_view.m_clip().m_back -
                             z_buffer)) { /* old origin was behind back plane */
      front_weight = 0.0F;
      delta_back = (I->m_view.m_clip().m_back - z_buffer) - old_origin;

    } else if (slab_width >= R_SMALL4) { /* otherwise, if slab exists */
      front_weight = (old_back - old_origin) /
                     slab_width; /* weight based on relative proximity */
    } else {
      front_weight = 0.5F;
    }

    back_weight = 1.0F - front_weight;

    if ((front_weight > 0.2) &&
        (back_weight > 0.2)) { /* origin not near edge */
      if (delta_front * delta_back >
          0.0F) { /* planes moving in same direction */
        if (fabs(delta_front) >
            fabs(delta_back)) { /* so stick with whichever moves less */
          v2[2] = delta_back;
        } else {
          v2[2] = delta_front;
        }
      } else {
        /* planes moving in opposite directions (increasing slab size) */
        /* don't move origin */
      }
    } else { /* origin is near edge -- move origin with plane having highest
                weight */
      if (front_weight < back_weight) {
        v2[2] = delta_back;
      } else {
        v2[2] = delta_front;
      }
    }

    old_pos2 = I->m_view.pos().z;

    MatrixInvTransformC44fAs33f3f(glm::value_ptr(I->m_view.rotMatrix()), v2,
                                  v2); /* transform offset into realspace */
    subtract3f(glm::value_ptr(I->m_view.origin()), v2,
               v2);              /* calculate new origin location */
    SceneOriginSet(G, v2, true); /* move origin, preserving camera location */

    if (SettingGetGlobal_b(G, cSetting_ortho) || zoom_flag) {
      /* we're orthoscopic, so we don't want the effective field of view
         to change.  Thus, we have to hold Pos[2] constant, and instead
         move the planes.
       */
      float delta = old_pos2 - I->m_view.pos().z;
      I->m_view.translate(0, 0, delta);
      SceneClipSet(G, I->m_view.m_clip().m_front - delta,
                   I->m_view.m_clip().m_back - delta);
    }
    slab_width = I->m_view.m_clip().m_back - I->m_view.m_clip().m_front;

    /* first, check to make sure that the origin isn't too close to either plane
     */
    if ((z_buffer * 2) > slab_width) z_buffer = slab_width * 0.5F;
  }
  if ((adjust_flag) && SettingGetGlobal_b(G, cSetting_roving_detail)) {
    SceneRovingPostpone(G);
  }
  if (SettingGetGlobal_b(G, cSetting_roving_detail)) {
    SceneRovingDirty(G);
  }
}

void SceneRovingChanged(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  SceneRovingDirty(G);
  I->RovingCleanupFlag = true;
}

static void SceneRovingCleanup(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  const char* s;
  char buffer[OrthoLineLength];

  I->RovingCleanupFlag = false;

  s = SettingGet_s(G, nullptr, nullptr, cSetting_roving_selection);

  sprintf(buffer, "cmd.hide('lines','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('sticks','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('spheres','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('ribbon','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('cartoon','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('labels','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('nonbonded','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
  sprintf(buffer, "cmd.hide('nb_spheres','''%s''')", s);
  PParse(G, buffer);
  PFlush(G);
}

void SceneRovingUpdate(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  char buffer[OrthoLineLength];
  float sticks, lines, spheres, labels, ribbon, cartoon;
  float polar_contacts, polar_cutoff, nonbonded, nb_spheres;
  char byres[10] = "byres";
  char not_[4] = "not";
  char empty[1] = "";
  char* p1;
  char* p2;
  const char* s;
  int refresh_flag = false;
  const char* name;
  float level;
  float isosurface, isomesh;

  if (I->RovingDirtyFlag &&
      ((UtilGetSeconds(G) - I->RovingLastUpdate) >
       fabs(SettingGetGlobal_f(G, cSetting_roving_delay)))) {
    if (I->RovingCleanupFlag) SceneRovingCleanup(G);

    s = SettingGet_s(G, nullptr, nullptr, cSetting_roving_selection);
    sticks = SettingGetGlobal_f(G, cSetting_roving_sticks);
    lines = SettingGetGlobal_f(G, cSetting_roving_lines);
    labels = SettingGetGlobal_f(G, cSetting_roving_labels);
    spheres = SettingGetGlobal_f(G, cSetting_roving_spheres);
    ribbon = SettingGetGlobal_f(G, cSetting_roving_ribbon);
    cartoon = SettingGetGlobal_f(G, cSetting_roving_cartoon);
    polar_contacts = SettingGetGlobal_f(G, cSetting_roving_polar_contacts);
    polar_cutoff = SettingGetGlobal_f(G, cSetting_roving_polar_cutoff);
    nonbonded = SettingGetGlobal_f(G, cSetting_roving_nonbonded);
    nb_spheres = SettingGetGlobal_f(G, cSetting_roving_nb_spheres);

    isomesh = SettingGetGlobal_f(G, cSetting_roving_isomesh);
    isosurface = SettingGetGlobal_f(G, cSetting_roving_isosurface);

    if (SettingGetGlobal_b(G, cSetting_roving_byres))
      p2 = byres;
    else
      p2 = empty;

    if (sticks != 0.0F) {
      if (sticks < 0.0F) {
        p1 = not_;
        sticks = (float)fabs(sticks);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('sticks','''%s''');cmd.show('sticks','%s & enabled & "
              "%s %s (center expand %1.3f)')",
              s, s, p1, p2, sticks);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (lines != 0.0F) {
      if (lines < 0.0F) {
        p1 = not_;
        lines = (float)fabs(lines);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('lines','''%s''');cmd.show('lines','%s & enabled & %s "
              "%s (center expand %1.3f)')",
              s, s, p1, p2, lines);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (labels != 0.0F) {
      if (labels < 0.0F) {
        p1 = not_;
        labels = (float)fabs(labels);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('labels','''%s''');cmd.show('labels','%s & enabled & "
              "%s %s (center expand %1.3f)')",
              s, s, p1, p2, labels);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (spheres != 0.0F) {
      if (spheres < 0.0F) {
        p1 = not_;
        spheres = (float)fabs(spheres);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('spheres','''%s''');cmd.show('spheres','%s & enabled & "
              "%s %s (center expand %1.3f)')",
              s, s, p1, p2, spheres);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (cartoon != 0.0F) {
      if (cartoon < 0.0F) {
        p1 = not_;
        cartoon = (float)fabs(cartoon);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('cartoon','''%s''');cmd.show('cartoon','%s & enabled & "
              "%s %s (center expand %1.3f)')",
              s, s, p1, p2, cartoon);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (ribbon != 0.0F) {
      if (ribbon < 0.0F) {
        p1 = not_;
        ribbon = (float)fabs(ribbon);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('ribbon','''%s''');cmd.show('ribbon','%s & enabled & "
              "%s %s (center expand %1.3f)')",
              s, s, p1, p2, ribbon);
      PParse(G, buffer);
      PFlush(G);

      refresh_flag = true;
    }

    if (polar_contacts != 0.0F) {
      int label_flag = 0;
      if (polar_contacts < 0.0F) {
        p1 = not_;
        polar_contacts = (float)fabs(polar_contacts);
      } else {
        p1 = empty;
      }
      if (polar_cutoff < 0.0F) {
        label_flag = true;
        polar_cutoff = (float)fabs(polar_cutoff);
      }
      sprintf(buffer,
              "cmd.delete('rov_pc');cmd.dist('rov_pc','%s & enabled & %s %s "
              "(center expand %1.3f)','same',%1.4f,mode=2,label=%d,quiet=2)",
              s, p1, p2, polar_contacts, polar_cutoff, label_flag);
      PParse(G, buffer);
      PFlush(G);

      refresh_flag = true;
    }

    if (nonbonded != 0.0F) {
      if (nonbonded < 0.0F) {
        p1 = not_;
        nonbonded = (float)fabs(nonbonded);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('nonbonded','''%s''');cmd.show('nonbonded','%s & "
              "enabled & %s %s (center expand %1.3f)')",
              s, s, p1, p2, nonbonded);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (nb_spheres != 0.0F) {
      if (nb_spheres < 0.0F) {
        p1 = not_;
        nb_spheres = (float)fabs(nb_spheres);
      } else {
        p1 = empty;
      }
      sprintf(buffer,
              "cmd.hide('nb_spheres','''%s''');cmd.show('nb_spheres','%s & "
              "enabled & %s %s (center expand %1.3f)')",
              s, s, p1, p2, nb_spheres);
      PParse(G, buffer);
      PFlush(G);
      refresh_flag = true;
    }

    if (isomesh != 0.0F) {
      int auto_save;

      auto_save = SettingGetGlobal_i(G, cSetting_auto_zoom);
      SettingSetGlobal_i(G, cSetting_auto_zoom, 0);

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map1_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map1_level);
            sprintf(buffer, "cmd.isomesh('rov_m1','%s',%8.6f,'center',%1.3f)",
                    name, level, isomesh);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map2_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map2_level);
            sprintf(buffer, "cmd.isomesh('rov_m2','%s',%8.6f,'center',%1.3f)",
                    name, level, isomesh);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map3_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map3_level);
            sprintf(buffer, "cmd.isomesh('rov_m3','%s',%8.6f,'center',%1.3f)",
                    name, level, isomesh);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }
      SettingSetGlobal_i(G, cSetting_auto_zoom, auto_save);
    }

    if (isosurface != 0.0F) {
      int auto_save;

      auto_save = SettingGetGlobal_i(G, cSetting_auto_zoom);
      SettingSetGlobal_i(G, cSetting_auto_zoom, 0);

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map1_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map1_level);
            sprintf(buffer,
                    "cmd.isosurface('rov_s1','%s',%8.6f,'center',%1.3f)", name,
                    level, isosurface);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map2_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map2_level);
            sprintf(buffer,
                    "cmd.isosurface('rov_s2','%s',%8.6f,'center',%1.3f)", name,
                    level, isosurface);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }

      name = SettingGet_s(G, nullptr, nullptr, cSetting_roving_map3_name);
      if (name)
        if (name[0])
          if (ExecutiveFindObjectByName(G, name)) {
            level = SettingGetGlobal_f(G, cSetting_roving_map3_level);
            sprintf(buffer,
                    "cmd.isosurface('rov_s3','%s',%8.6f,'center',%1.3f)", name,
                    level, isosurface);
            PParse(G, buffer);
            PFlush(G);
            refresh_flag = true;
          }
      SettingSetGlobal_i(G, cSetting_auto_zoom, auto_save);
    }

    if (refresh_flag) {
      PParse(G, "cmd.refresh()");
      PFlush(G);
    }

    I->RovingLastUpdate = UtilGetSeconds(G);
    I->RovingDirtyFlag = false;
  }
}
int SceneRovingCheckDirty(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  return (I->RovingDirtyFlag);
}
