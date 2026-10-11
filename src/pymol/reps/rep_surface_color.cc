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

// Coloring surfaces by their atoms, ramps and transparency.

#include "pymol/base/feedback.h"
#include "pymol/reps/rep_surface.h"
#include "pymol/scene/scene.h"
#include "pymol/select/selector.h"
#include "pymol/reps/rep_surface_internal.h"

bool RepSurface::sameColor() const {
  if (ColorInvalidated) return false;

  const auto* lc = LastColor;
  for (int idx = 0; idx < cs->NIndex; idx++) {
    auto* ai = cs->getAtomInfo(idx);
    if (ai->visRep & cRepSurfaceBit) {
      if (*(lc++) != ai->color) {
        return false;
      }
    }
  }

  return true;
}

Rep* RepSurface::recolor() {
  auto const I = this;
  assert(cs == this->cs);

  // Reset the transparency flag before recomputing per-vertex alpha below. A
  // `transparency` setting change invalidates at cRepInvColor, which re-runs
  // recolor() in place on this same Rep (see Rep::update); the per-vertex loop
  // only ever calls setHasTransparency() (=true) when at_transp > 0 and never
  // clears it. Without this reset the flag stays stuck true after transparency
  // returns to 0, so the now-opaque surface keeps being routed into the
  // transparent pass and renders semi-transparent. Mirrors RepCartoon, which
  // passes the computed boolean to setHasTransparency().
  I->setHasTransparency(false);

  MapType *map = nullptr, *ambient_occlusion_map = nullptr;
  int a, i0, i, j, c1;
  float *v0, *vc, *va;
  const float* c0;
  float* n0;
  int* lc;
  char* lv;
  int first_color;
  float *v_pos, v_above[3];
  int ramp_above;
  ObjectMolecule* obj;
  float probe_radius;
  float dist;
  float cutoff;
  int inclH;
  int inclInvis;
  int cullByFlag = false;
  int surface_mode, ambient_occlusion_mode;
  int surface_color;
  int* present = nullptr;
  int* rc;
  int ramped_flag = false;

  int carve_state = 0;
  int carve_flag = false;
  float carve_cutoff;
  float carve_normal_cutoff;
  int carve_normal_flag;
  const char* carve_selection = nullptr;
  float* carve_vla = nullptr;
  MapType* carve_map = nullptr;

  int clear_state = 0;
  int clear_flag = false;
  float clear_cutoff;
  const char* clear_selection = nullptr;
  float* clear_vla = nullptr;

  int state = getState();

  float transp;
  int variable_alpha = false;

  MapType* clear_map = nullptr;

  AtomInfoType *ai2 = nullptr, *ai1;

  obj = cs->Obj;
  ambient_occlusion_mode =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                   cSetting_ambient_occlusion_mode);
  surface_mode = I->surface_mode;
  ramp_above = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                            cSetting_surface_ramp_above_mode);
  surface_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                   cSetting_surface_color);
  cullByFlag = (surface_mode == cRepSurface_by_flags);
  inclH = !((surface_mode == cRepSurface_heavy_atoms) ||
            (surface_mode == cRepSurface_vis_heavy_only));
  inclInvis = !((surface_mode == cRepSurface_vis_only) ||
                (surface_mode == cRepSurface_vis_heavy_only));
  probe_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_solvent_radius);
  I->proximity = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_surface_proximity);
  carve_cutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_surface_carve_cutoff);
  clear_cutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_surface_clear_cutoff);
  transp = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                        cSetting_transparency);
  carve_normal_cutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_surface_carve_normal_cutoff);
  carve_normal_flag = carve_normal_cutoff > (-1.0F);

  cutoff = I->max_vdw + 2 * probe_radius;

  if (!I->LastVisib) I->LastVisib = pymol::malloc<char>(cs->NIndex);
  if (!I->LastColor) I->LastColor = pymol::malloc<int>(cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  for (a = 0; a < cs->NIndex; a++) {
    ai2 = cs->getAtomInfo(a);
    *(lv++) = GET_BIT(ai2->visRep, cRepSurface);
    *(lc++) = ai2->color;
  }

  if (I->N) {
    if (carve_cutoff > 0.0F) {
      carve_state = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_surface_carve_state) -
                    1;
      carve_selection = SettingGet_s(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_surface_carve_selection);
      if (carve_selection)
        carve_map = SelectorGetSpacialMapFromSeleCoord(
            G, SelectorIndexByName(G, carve_selection), carve_state,
            carve_cutoff, &carve_vla);
      if (carve_map) MapSetupExpress(carve_map);
      carve_flag = true;
    }

    if (clear_cutoff > 0.0F) {
      clear_state = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                                 cSetting_surface_clear_state) -
                    1;
      clear_selection = SettingGet_s(G, cs->Setting.get(), obj->Setting.get(),
                                     cSetting_surface_clear_selection);
      if (clear_selection)
        clear_map = SelectorGetSpacialMapFromSeleCoord(
            G, SelectorIndexByName(G, clear_selection), clear_state,
            clear_cutoff, &clear_vla);
      if (clear_map) MapSetupExpress(clear_map);
      clear_flag = true;
    }

    if (I->VC.empty()) I->VC.resize(3 * I->N);
    vc = I->VC.data();
    if (I->VA.empty()) I->VA.resize(I->N);
    va = I->VA.data();
    if (I->RC.empty()) I->RC.resize(I->N);
    rc = I->RC.data();
    if (I->Vis.empty()) I->Vis.resize(I->N);
    if (ColorCheckRamped(G, surface_color)) {
      I->oneColorFlag = false;
    } else {
      I->oneColorFlag = true;
    }
    first_color = -1;

    present = pymol::malloc<int>(cs->NIndex);
    {
      int* ap = present;
      for (a = 0; a < cs->NIndex; a++) {
        ai1 = obj->AtomInfo + cs->IdxToAtm[a];
        if ((ai1->visRep & cRepSurfaceBit) && (inclH || (!ai1->isHydrogen())) &&
            ((!cullByFlag) ||
             (!(ai1->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
          *ap = 2;
        else
          *ap = 0;
        ap++;
      }
    }

    if (inclInvis) {
      float probe_radiusX2 = probe_radius * 2;
      map = new MapType(G, 2 * I->max_vdw + probe_radius, cs->Coord, cs->NIndex,
                        nullptr, present);
      MapSetupExpress(map);
      /* add in nearby invisibles */
      for (a = 0; a < cs->NIndex; a++) {
        if (!present[a]) {
          ai1 = obj->AtomInfo + cs->IdxToAtm[a];
          if ((!cullByFlag) || !(ai1->flags & cAtomFlag_ignore)) {
            v0 = cs->coordPtr(a);
            i = *(MapLocusEStart(map, v0));
            if (i && !map->EList.empty()) {
              j = map->EList[i++];
              while (j >= 0) {
                if (present[j] > 1) {
                  ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                  if (within3f(cs->coordPtr(j), v0,
                               ai1->vdw + ai2->vdw + probe_radiusX2)) {
                    present[a] = 1;
                    break;
                  }
                }
                j = map->EList[i++];
              }
            }
          }
        }
      }
      MapFree(map);
      map = nullptr;
    }

    /**
     * Update the ambient occlusion accessibility array (VAO)
     */
    auto update_VAO = [&]() {
      if (!ambient_occlusion_mode) {
        I->VAO.clear();
        return;
      }

      float maxDist = 0.f, maxDistA = 0.f;
      int level_min = 64, level_max = 0;
      double start_time, cur_time;
      short vertex_map = 0; /* vertex or atom map */
      float map_cutoff = cutoff;

      switch (ambient_occlusion_mode % 4) {
        case 1:
        case 3:
          vertex_map =
              0; /* ambient_occlusion_mode - 1 atoms in map (default), 2
                    - vertices in map */
          break;
        case 2:
          vertex_map = 1;
      }

      I->VAO.resize(I->N);
      start_time = UtilGetSeconds(G);

      if (vertex_map) {
        ambient_occlusion_map =
            new MapType(G, map_cutoff, I->V.data(), I->N, nullptr, nullptr);
      } else {
        ambient_occlusion_map =
            new MapType(G, map_cutoff, cs->Coord, cs->NIndex, nullptr, present);
      }
      MapSetupExpress(ambient_occlusion_map);

      if (ambient_occlusion_mode == 3) {
        /* per atom */
        float* VAO = pymol::malloc<float>(cs->NIndex);
        short* nVAO = pymol::malloc<short>(cs->NIndex);
        memset(VAO, 0, sizeof(float) * cs->NIndex);
        memset(nVAO, 0, sizeof(short) * cs->NIndex);

        for (a = 0; a < I->N; a++) {
          int nbits = 0;
          short level1, level2, has;
          unsigned long bits = 0L, bit;
          float d[3], *vn0, v0mod[3];
          int closeA = -1;
          float closeDist = FLT_MAX;
          has = 0;

          v0 = I->V.data() + 3 * a;
          vn0 = I->VN.data() + 3 * a;
          mult3f(vn0, .01f, v0mod);
          add3f(v0, v0mod, v0mod);

          i = *(MapLocusEStart(ambient_occlusion_map, v0));
          if (i && !map->EList.empty()) {
            j = ambient_occlusion_map->EList[i++];
            while (j >= 0) {
              subtract3f(cs->coordPtr(j), v0, d);
              dist = (float)length3f(d);
              if (dist < closeDist) {
                closeA = j;
                closeDist = dist;
              }
              j = ambient_occlusion_map->EList[i++];
            }
          }
          if (closeA >= 0) {
            if (nVAO[closeA]) {
              I->VAO[a] = VAO[closeA];
            } else {
              v0 = cs->coordPtr(closeA);
              i = *(MapLocusEStart(ambient_occlusion_map, v0));
              if (i) {
                j = ambient_occlusion_map->EList[i++];
                while (j >= 0) {
                  if (closeA == j) {
                    j = ambient_occlusion_map->EList[i++];
                    continue;
                  }
                  subtract3f(cs->coordPtr(j), v0, d);
                  dist = (float)length3f(d);
                  if (dist > 12.f) {
                    j = ambient_occlusion_map->EList[i++];
                    continue;
                  }
                  has = 1;
                  level1 = (d[2] < 0.f) ? 4 : 0;
                  level1 |= (d[1] < 0.f) ? 2 : 0;
                  level1 |= (d[0] < 0.f) ? 1 : 0;
                  d[0] = fabs(d[0]);
                  d[1] = fabs(d[1]);
                  d[2] = fabs(d[2]);
                  level2 = (d[0] <= d[1]) ? 4 : 0;
                  level2 |= (d[1] <= d[2]) ? 2 : 0;
                  level2 |= (d[0] <= d[2]) ? 1 : 0;

                  bit = level1 * 8 + level2;
                  bits |= (1L << bit);
                  j = ambient_occlusion_map->EList[i++];
                }
              }
              if (has) {
                nbits = countBits(bits);
                if (nbits > level_max) level_max = nbits;
                if (nbits < level_min) level_min = nbits;
                I->VAO[a] = nbits;
              } else {
                level_min = 0;
                I->VAO[a] = 0.f;
              }
              VAO[closeA] = I->VAO[a];
              nVAO[closeA] = 1;
            }
          }
        }
        FreeP(VAO);
        FreeP(nVAO);
      } else {
        float natomsL = 0;
        for (a = 0; a < I->N; a++) {
          int natoms = 0, nbits = 0;
          short level1, level2, has;
          unsigned long bits = 0L, bit;
          float d[3], *vn0, pt[3], v0mod[3];

          if (a % 1000 == 0) {
            PRINTFB(G, FB_RepSurface, FB_Debugging)
            "RepSurfaceColor():  Ambient Occlusion computing mode=%d "
            "#vertices=%d done=%d\n",
                ambient_occlusion_mode, I->N, a ENDFB(G);
          }
          v0 = I->V.data() + 3 * a;
          vn0 = I->VN.data() + 3 * a;
          mult3f(vn0, .01f, v0mod);
          add3f(v0, v0mod, v0mod);
          i = *(MapLocusEStart(ambient_occlusion_map, v0));
          if (i && !ambient_occlusion_map->EList.empty()) {
            j = ambient_occlusion_map->EList[i++];
            maxDistA = 0.f;
            has = 0;
            while (j >= 0) {
              natomsL++;
              if (vertex_map && a == j) {
                j = ambient_occlusion_map->EList[i++];
                continue;
              }
              if (vertex_map) {
                copy3f(I->V.data() + j * 3, pt);
                subtract3f(I->V.data() + j * 3, v0mod, d);
              } else {
                copy3f(cs->coordPtr(j), pt);
                subtract3f(cs->coordPtr(j), v0mod, d);
              }
              dist = (float)length3f(d);
              normalize3f(d);
              if (get_angle3f(vn0, d) >= (PI / 2.f)) {
                j = ambient_occlusion_map->EList[i++];
                continue;
              }
              if (dist <= .0001f || dist > 12.f) {
                j = ambient_occlusion_map->EList[i++];
                continue;
              }
              has = 1;
              (dist > maxDistA) ? maxDistA = dist : 0;
              level1 = (d[2] < 0.f) ? 4 : 0;
              level1 |= (d[1] < 0.f) ? 2 : 0;
              level1 |= (d[0] < 0.f) ? 1 : 0;
              d[0] = fabs(d[0]);
              d[1] = fabs(d[1]);
              d[2] = fabs(d[2]);
              level2 = (d[0] <= d[1]) ? 4 : 0;
              level2 |= (d[1] <= d[2]) ? 2 : 0;
              level2 |= (d[0] <= d[2]) ? 1 : 0;

              bit = level1 * 8 + level2;
              bits |= (1L << bit);
              j = ambient_occlusion_map->EList[i++];
              natoms++;
            }
            if (has) {
              nbits = countBits(bits);
              if (nbits > level_max) level_max = nbits;
              if (nbits < level_min) level_min = nbits;
              (maxDistA > maxDist) ? maxDist = maxDistA : 0;
              I->VAO[a] = nbits;
            } else {
              level_min = 0;
              I->VAO[a] = 0.f;
            }
            maxDistA = 0.f;
          }
        }
        PRINTFB(G, FB_RepSurface, FB_Debugging)
        "RepSurfaceColor():  #vertices=%d Ambient Occlusion average #atoms "
        "looked at per vertex = %f\n",
            I->N, (natomsL / I->N) ENDFB(G);
      }
      {
        int ambient_occlusion_smooth =
            SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_ambient_occlusion_smooth);
        int surface_quality = SettingGet_i(
            G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_quality);
        float min_max_diff = level_max - level_min;
        if (surface_quality > 0) ambient_occlusion_smooth *= surface_quality;
        /* Now we should set VAO from min/max */
        if (min_max_diff) {
          for (a = 0; a < I->N; a++) {
            I->VAO[a] = ((I->VAO[a] - level_min) / min_max_diff);
          }
        } else {
          for (a = 0; a < I->N; a++) {
            I->VAO[a] = 0.f;
          }
        }

        /* SMOOTH Accessibility VALUES */
        if (ambient_occlusion_smooth && !I->T.empty()) {
          int i, j, pt1, pt2, pt3;
          float ave;
          float* tmpVAO = pymol::malloc<float>(I->N);
          int *nVAO = pymol::malloc<int>(I->N), c, *t;

          for (j = 0; j < ambient_occlusion_smooth; j++) {
            memset(nVAO, 0, sizeof(int) * I->N);
            memset(tmpVAO, 0, sizeof(float) * I->N);

            t = I->T.data();
            c = I->NT;
            while (c--) {
              if (I->allVisibleFlag ||
                  visibility_test(I->proximity, I->Vis.data(), t)) {
                pt1 = *t;
                pt2 = *(t + 1);
                pt3 = *(t + 2);
                nVAO[pt1] += 1;
                nVAO[pt2] += 1;
                nVAO[pt3] += 1;

                ave = ave3(I->VAO[pt1], I->VAO[pt2], I->VAO[pt3]);

                tmpVAO[pt1] += ave;
                tmpVAO[pt2] += ave;
                tmpVAO[pt3] += ave;
              }
              t += 3;
            }

            for (i = 0; i < I->N; i++) {
              if (nVAO[i]) {
                /* only update if added and greater than original */
                //		if ((tmpVAO[i]/(float)nVAO[i]) > I->VAO[i])
                I->VAO[i] = tmpVAO[i] / (float)nVAO[i];
              }
            }
          }
          FreeP(tmpVAO);
          FreeP(nVAO);
        }
      }
      MapFree(ambient_occlusion_map);
      cur_time = UtilGetSeconds(G);

      PRINTFB(G, FB_RepSurface, FB_Debugging)
      "RepSurfaceColor():  Ambient Occlusion computed #atoms=%d #vertices=%d "
      "time=%lf seconds\n",
          cs->NIndex, I->N, (cur_time - start_time) ENDFB(G);
      ambient_occlusion_map = nullptr;
    };  // update_VAO

    /* now, assign colors to each point */
    map = new MapType(G, cutoff, cs->Coord, cs->NIndex, nullptr, present);
    if (map) {
      short color_smoothing =
          SettingGetGlobal_i(G, cSetting_surface_color_smoothing);
      float color_smoothing_threshold =
          SettingGetGlobal_f(G, cSetting_surface_color_smoothing_threshold);
      int atm, ok = true;
      MapSetupExpress(map);
      ok &= !G->Interrupt;
      if (ok && I->AT.empty()) I->AT.resize(I->N);
      for (a = 0; ok && a < I->N; a++) {
        float at_transp = transp;

        AtomInfoType* ai0 = nullptr;
        float minDist = FLT_MAX, minDist2 = FLT_MAX, distDiff = FLT_MAX;
        int pi = -1, catm = -1; /* variables for color smoothing */
        AtomInfoType *pai = nullptr,
                     *pai2 = nullptr; /* variables for color smoothing */
        c1 = 1;
        i0 = -1;
        v0 = I->V.data() + 3 * a;
        n0 = I->VN.data() + 3 * a;
        auto vi = &I->Vis[a];
        /* colors */
        i = *(MapLocusEStart(map, v0));
        if (i && !map->EList.empty()) {
          j = map->EList[i++];
          while (j >= 0) {
            atm = cs->IdxToAtm[j];
            ai2 = obj->AtomInfo + atm;
            if ((inclH || (!ai2->isHydrogen())) &&
                ((!cullByFlag) || (!(ai2->flags & cAtomFlag_ignore)))) {
              dist = (float)diff3f(v0, cs->coordPtr(j)) - ai2->vdw;
              if (color_smoothing) {
                if (dist < minDist) {
                  /* switching closest to 2nd closest */
                  pai2 = pai;
                  minDist2 = minDist;
                  pi = j;
                  catm = atm;
                  pai = ai2;
                  minDist = dist;
                } else if (dist < minDist2) {
                  /* just setting second closest */
                  pai2 = ai2;
                  minDist2 = dist;
                }
              } else if (dist < minDist) {
                i0 = j;
                catm = atm;
                ai0 = ai2;
                minDist = dist;
              }
            }
            j = map->EList[i++];
          }
        }
        I->AT[a] = catm;
        if (color_smoothing) {
          i0 = pi;
          ai0 = pai;
          /* TODO: should find point closest to v0 between
                   atoms points (cs->coordPtr(pi)) and (cs->coordPtr(pi2))
             (including vdw), then set this distance to the distance between the
             vertex v0 and that point. We might want to use the normal to
             compute this distance.
           */
          distDiff = fabs(minDist2 - minDist);
        }
        if (i0 >= 0) {
          int at_surface_color =
              AtomSettingGetWD(G, ai0, cSetting_surface_color, surface_color);
          at_transp = AtomSettingGetWD(G, ai0, cSetting_transparency, transp);

          if (at_surface_color != -1) {
            c1 = at_surface_color;
            distDiff = FLT_MAX;
          } else {
            c1 = ai0->color;
          }

          if (I->oneColorFlag) {
            if (first_color >= 0) {
              if (first_color != c1) I->oneColorFlag = false;
            } else
              first_color = c1;
          }
          if (I->allVisibleFlag)
            *vi = 1;
          else {
            ai2 = obj->AtomInfo + cs->IdxToAtm[i0];
            if ((ai2->visRep & cRepSurfaceBit) &&
                (inclH || (!ai2->isHydrogen())) &&
                ((!cullByFlag) ||
                 (!(ai2->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
              *vi = 1;
            else
              *vi = 0;
          }
        } else {
          *vi = 0;
        }
        if (carve_flag && (*vi)) { /* is point visible, and are we carving? */
          *vi = 0;

          if (carve_map) {
            i = *(MapLocusEStart(carve_map, v0));
            if (i && !carve_map->EList.empty()) {
              j = carve_map->EList[i++];
              while (j >= 0) {
                float* v_targ = carve_vla + 3 * j;
                if (within3f(v_targ, v0, carve_cutoff)) {
                  if (!carve_normal_flag) {
                    *vi = 1;
                    break;
                  } else {
                    float v_to[3];
                    subtract3f(v_targ, v0, v_to);
                    if (dot_product3f(v_to, n0) >= carve_normal_cutoff) {
                      *vi = 1;
                      break;
                    }
                  }
                }
                j = carve_map->EList[i++];
              }
            }
          }
        }
        if (clear_flag && (*vi)) { /* is point visible, and are we clearing? */
          if (clear_map) {
            i = *(MapLocusEStart(clear_map, v0));
            if (i && !clear_map->EList.empty()) {
              j = clear_map->EList[i++];
              while (j >= 0) {
                if (within3f(clear_vla + 3 * j, v0, clear_cutoff)) {
                  *vi = 0;
                  break;
                }
                j = clear_map->EList[i++];
              }
            }
          }
        }
        /*
           if(ColorCheckRamped(G,surface_color)) {
           c1 = surface_color;
           }
         */

        if (ColorCheckRamped(G, c1)) {
          I->oneColorFlag = false;
          switch (ramp_above) {
            case 1:
              copy3f(n0, v_above);
              scale3f(v_above, probe_radius, v_above);
              add3f(v0, v_above, v_above);
              v_pos = v_above;
              rc[0] = -1;
              break;
            default:
              v_pos = v0;
              rc[0] = c1;
              ramped_flag = true;
              break;
          }
          ColorGetRamped(G, c1, v_pos, vc, state);
          vc += 3;
          rc++;
        } else {
          if (color_smoothing && distDiff < color_smoothing_threshold && pai2) {
            const float* c2;
            float weight, weight2;
            if (color_smoothing == 1) {
              weight =
                  1.f + sin(.5f * PI * (distDiff / color_smoothing_threshold));
            } else {
              weight = 1.f + (distDiff / color_smoothing_threshold);
            }
            weight2 = 2.f - weight;
            c0 = ColorGet(G, c1);
            c2 = ColorGet(G, pai2->color);
            *(rc++) = c1;
            *(vc++) = ((weight * (*(c0++))) + (weight2 * (*(c2++)))) / 2.f;
            *(vc++) = ((weight * (*(c0++))) + (weight2 * (*(c2++)))) / 2.f;
            *(vc++) = ((weight * (*(c0++))) + (weight2 * (*(c2++)))) / 2.f;
          } else {
            c0 = ColorGet(G, c1);
            *(rc++) = c1;
            *(vc++) = *(c0++);
            *(vc++) = *(c0++);
            *(vc++) = *(c0++);
          }
        }
        if (at_transp != transp) variable_alpha = true;
        *(va++) = 1.0F - at_transp;

        if (at_transp > 0) {
          I->setHasTransparency();
        }

        if (!*vi) I->allVisibleFlag = false;
      }
      MapFree(map);
    }
    if (variable_alpha) I->oneColorFlag = false;

    if (I->oneColorFlag) {
      I->oneColor = first_color;
    }

    // ambient occlusion
    update_VAO();
  }
  /*
     if(surface_color>=0) {
     I->oneColorFlag=true;
     I->oneColor=surface_color;
     }
   */

  if (G->HaveGUI) {
    setShaderCGO(I, nullptr);
  }

  if (!I->VA.empty() && (!variable_alpha)) {
    I->VA.clear();
  }

  if (carve_map) MapFree(carve_map);
  VLAFreeP(carve_vla);
  if (clear_map) MapFree(clear_map);
  VLAFreeP(clear_vla);
  if ((!ramped_flag) || (!SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                       cSetting_ray_color_ramps))) {
    I->RC.clear();
  }
  I->ColorInvalidated = false;
  FreeP(present);

  return this;
}
