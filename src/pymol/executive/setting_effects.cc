#include "pymol/executive/setting_effects.h"

#include "pymol/platform/os_python.h"
#include "pymol/platform/os_proprietary.h"
#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"
#include "pymol/base/base.h"
#include "pymol/base/err.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/ui/ortho.h"
#include "pymol/settings/setting.h"
#include "pymol/color/color.h"
#include "pymol/scene/scene.h"
#include "pymol/ui/but_mode.h"
#include "pymol/graphics/cgo.h"
#include "pymol/executive/executive.h"
#include "pymol/executive/editor.h"
#include "pymol/python/p.h"
#include "pymol/base/util.h"
#include "pymol/app/main.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/wizard.h"
#include "pymol/ui/seq.h"
#include "pymol/base/pymol_options.h"
#include "pymol/graphics/shader_mgr.h"
#include "pymol/geometry/sphere.h"
#include "pymol/select/selector.h"
#include "pymol/base/parse.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"
#include <thread>

void ColorUpdateFrontFromSettings(PyMOLGlobals* G) {
  int bg_gradient = SettingGet_b(G, nullptr, nullptr, cSetting_bg_gradient);
  const char* bg_image_filename =
      SettingGet_s(G, nullptr, nullptr, cSetting_bg_image_filename);
  short bg_image = bg_image_filename && bg_image_filename[0];

  if (!bg_gradient) {
    if (!bg_image && !OrthoBackgroundDataIsSet(*G->Ortho)) {
      const float* v =
          ColorGet(G, SettingGet_color(G, nullptr, nullptr, cSetting_bg_rgb));
      ColorUpdateFront(G, v);
    } else {
      float v[] = {0.f, 0.f, 0.f};
      ColorUpdateFront(G, v);
    }
  } else {
    // Copied before the second lookup (#542): for two 24-bit RGB colours
    // ColorGet returns the same scratch buffer twice, so the second call
    // overwrote the first and the "average" was just the top colour.
    float vv[3], v[3];
    copy3f(ColorGet(G, SettingGet_color(G, nullptr, nullptr,
                                        cSetting_bg_rgb_bottom)),
           v);
    const float* vb =
        ColorGet(G, SettingGet_color(G, nullptr, nullptr, cSetting_bg_rgb_top));
    average3f(v, vb, vv);
    ColorUpdateFront(G, vv);
  }
}

static int SettingCheckUseShaders(CSetting* I, int quiet) {
  PyMOLGlobals* G = I->G;
  if (SettingGetGlobal_i(G, cSetting_use_shaders)) {
    if (G->ShaderMgr->IsConfigured() && !G->ShaderMgr->ShadersPresent()) {
      SettingSet_b(I, cSetting_use_shaders, 0);
      if (!quiet) {
        PRINTFB(G, FB_Setting, FB_Warnings)
        "Setting-Error: use_shaders cannot be set when Shaders are not "
        "available, setting use_shaders back to false\n" ENDFB(G);
      }
      return 1;
    }
  }
  return 0;
}
/**
 * Get the indices of all settings that have changed since last calling
 * this function. Resets the "changed" flag.
 *
 * @param name object name or nullptr/"" for global settings
 * @param state object state
 */
std::vector<int> SettingGetUpdateList(PyMOLGlobals* G, const char* name,
                                      int state) {
  CSetting* I = G->Setting;
  pymol::copyable_ptr<CSetting>* handle;
  std::vector<int> result;

  if (name && name[0]) {
    // object-state settings

    pymol::CObject* obj = ExecutiveFindObjectByName(G, name);

    if (!obj || !(handle = obj->getSettingHandle(state)) ||
        !(I = handle->get()))
      // not found -> empty list
      return result;
  }

  for (int a = 0; a < cSetting_INIT; ++a) {
    if (I->info[a].changed) {
      I->info[a].changed = false;
      result.push_back(a);
    }
  }
  return (result);
}
void SettingGenerateSideEffects(PyMOLGlobals* G, int index, const char* sele,
                                int state, int quiet) {
  const char* inv_sele = (sele && sele[0]) ? sele : cKeywordAll;
  auto& rec = SettingGetInfo(index);

  if (rec.level == cSettingLevel_unused) {
    const char* name = rec.name;

    if (!quiet && name && name[0]) {
      PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: '%s' is no longer used\n", name ENDFB(G);
    }

    return;
  }

  // range check for int (global only)
  if (rec.type == cSetting_int && rec.hasMinMax() && !(sele && sele[0])) {
    int value = SettingGetGlobal_i(G, index);
    bool clamped = true;

    if (value < rec.value.i[1]) {
      value = rec.value.i[1];
    } else if (value > rec.value.i[2]) {
      value = rec.value.i[2];
    } else {
      clamped = false;
    }

    if (clamped) {
      PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: %s range = [%d,%d]; setting to %d.\n", rec.name,
          rec.value.i[1], rec.value.i[2], value ENDFB(G);
      SettingSetGlobal_i(G, index, value);
    }
  }

  switch (index) {
    case cSetting_stereo:
      SceneUpdateStereo(G);
      break;
    case cSetting_pick_surface:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        SceneInvalidatePicking(G);  // right now, when pick_surface is off,
                                    // wipes each CGO's pickColor array
      }
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
      break;
    case cSetting_pickable:
      ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvAll);
      SceneChanged(G);
      break;
    case cSetting_grid_mode:
      if (!SettingGetGlobal_i(G, cSetting_grid_mode))
        G->ShaderMgr->ResetUniformSet();
      ExecutiveUpdateSceneMembers(G);
    case cSetting_grid_slot:
      ExecutiveInvalidateGroups(G, false);
      ExecutiveUpdateSceneMembers(G);
      SceneChanged(G);
      break;
    case cSetting_grid_max:
      SceneChanged(G);
      break;
    case cSetting_defer_builds_mode:
      ExecutiveRebuildAll(G);
      break;
    case cSetting_seq_view:
    case cSetting_seq_view_label_spacing:
    case cSetting_seq_view_label_mode:
    case cSetting_seq_view_label_start:
    case cSetting_seq_view_format:
    case cSetting_seq_view_color:
    case cSetting_seq_view_unaligned_mode:
    case cSetting_seq_view_gap_mode:
      SeqChanged(G);
      break;
    case cSetting_seq_view_fill_color:
    case cSetting_seq_view_fill_char:
    case cSetting_seq_view_label_color:
      OrthoDirty(G);
      break;
    case cSetting_show_frame_rate:
      OrthoDirty(G);
      break;
    case cSetting_group_full_member_names:
    case cSetting_group_arrow_prefix:
      OrthoDirty(G);
      break;
    case cSetting_static_singletons:
      SeqChanged(G);
      break;
    case cSetting_seq_view_location:
      PParse(G, "cmd.viewport(-1,-1)");
      SeqChanged(G);
      break;
    case cSetting_seq_view_overlay:
      PParse(G, "cmd.viewport(-1,-1)");
      break;
    case cSetting_stereo_mode:
    case cSetting_anaglyph_mode:
      SceneUpdateStereoMode(G);
      OrthoInvalidateDoDraw(G);
      OrthoDirty(G);
      PyMOL_NeedRedisplay(G->PyMOL);
      break;
    case cSetting_precomputed_lighting:
    case cSetting_light_count:
    case cSetting_spec_count:
    case cSetting_dot_lighting:
    case cSetting_mesh_lighting:
    case cSetting_cgo_lighting:
    case cSetting_field_of_view:
    case cSetting_fog_start:
    case cSetting_two_sided_lighting:
    case cSetting_transparency_global_sort:
    case cSetting_dot_normals:
    case cSetting_mesh_normals:
      SceneInvalidate(G);
      break;
    case cSetting_spec_power:
      if (!quiet) {
        PRINTFB(G, FB_Setting, FB_Debugging)
        "Setting-Details: spec_power is depreciated in PyMOL 1.5.  This option "
        "will not work in future versions. Please set shininess to set the "
        "specular exponent for movable light sources.\n" ENDFB(G);
      }
      SceneInvalidate(G);
      break;
    case cSetting_light:
    case cSetting_light2:
    case cSetting_light3:
    case cSetting_light4:
    case cSetting_light5:
    case cSetting_light6:
    case cSetting_light7:
    case cSetting_reflect:
    case cSetting_direct:
    case cSetting_ambient:
    case cSetting_specular:
    case cSetting_specular_intensity:
    case cSetting_shininess:
    case cSetting_spec_reflect:
    case cSetting_spec_direct:
    case cSetting_spec_direct_power:
    case cSetting_power:
    case cSetting_reflect_power:
      SceneInvalidate(G);
      break;
    case cSetting_use_display_lists:
    case cSetting_simplify_display_lists:
    case cSetting_excl_display_lists_shaders:
      if (!quiet) {
        PRINTFB(G, FB_Setting, FB_Debugging)
        "Setting-Details: display lists were depreciated in PyMOL 1.7.x.  The "
        "settings use_display_lists, simplify_display_lists, and "
        "excl_display_lists_shaders no longer work.\n" ENDFB(G);
      }
    case cSetting_use_geometry_shaders:
      if (SettingGetGlobal_i(G, cSetting_use_geometry_shaders) &&
          G->ShaderMgr->IsConfigured() &&
          !G->ShaderMgr->GeometryShadersPresent()) {
        SettingSet_b(G->Setting, cSetting_use_geometry_shaders, 0);
        if (!quiet) {
          PRINTFB(G, FB_Setting, FB_Warnings)
          "Setting-Error: geometry shaders not available\n" ENDFB(G);
        }
        return;
      }

      ExecutiveInvalidateRep(G, inv_sele, cRepLabel, cRepInvRep);
      {
        // check if lines need to be invalidated
        bool line_as_cylinders =
            SettingGetGlobal_b(G, cSetting_use_shaders) &&
            SettingGetGlobal_b(G, cSetting_render_as_cylinders) &&
            SettingGetGlobal_b(G, cSetting_line_as_cylinders);
        if (!line_as_cylinders && !SettingGetGlobal_b(G, cSetting_trilines)) {
          ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
        }
      }
      break;
    case cSetting_shaders_from_disk:
      SceneInvalidate(G);
      break;
    case cSetting_use_shaders: {
      short changed = 0;
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        if (SettingCheckUseShaders(G->Setting, quiet)) {
          return;
        }
      }
      SceneInvalidate(G);
      if (SettingGetGlobal_b(G, cSetting_sphere_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_ribbon_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_nonbonded_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader)) {
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_dash_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepAngle, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepDihedral, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepDash, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_line_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_cartoon_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_cgo_use_shader) ||
          SettingGetGlobal_b(G, cSetting_stick_as_cylinders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_stick_use_shader) ||
          SettingGetGlobal_b(G, cSetting_stick_as_cylinders) ||
          (SettingGetGlobal_b(G, cSetting_stick_ball) &&
           SettingGetGlobal_b(G, cSetting_valence))) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
        changed = 1;
      }
      if (SettingGetGlobal_b(G, cSetting_surface_use_shader) ||
          SettingGetGlobal_b(G, cSetting_dot_use_shader) ||
          SettingGetGlobal_b(G, cSetting_mesh_use_shader)) {
        changed = 1;
      }
      if (changed) {
        SceneChanged(G);
      }
    } break;
    case cSetting_stereo_shift:
    case cSetting_stereo_angle:
    case cSetting_stereo_dynamic_strength:
      SceneInvalidate(G);
      break;
    case cSetting_scene_buttons:
    case cSetting_scene_buttons_mode:
      SceneInvalidate(G);
      OrthoInvalidateDoDraw(G);
      break;
    case cSetting_dash_round_ends:
    case cSetting_dash_color:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_dash_use_shader)) {
        ExecutiveInvalidateRep(G, "all", cRepDash, cRepInvRep);
      }
      SceneInvalidate(G);
      break;
    case cSetting_angle_color:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_dash_use_shader)) {
        ExecutiveInvalidateRep(G, "all", cRepAngle, cRepInvRep);
      }
      SceneInvalidate(G);
      break;
    case cSetting_dihedral_color:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_dash_use_shader)) {
        ExecutiveInvalidateRep(G, "all", cRepDihedral, cRepInvRep);
      }
      SceneInvalidate(G);
      break;
    case cSetting_mouse_selection_mode:
      OrthoDirty(G);
      break;
    case cSetting_internal_gui_control_size:
      WizardRefresh(G);
      OrthoDirty(G);
      break;
    case cSetting_hide_underscore_names:
      ExecutiveInvalidateGroups(G, false);
      OrthoDirty(G);
      break;
    case cSetting_gradient_spacing:
    case cSetting_gradient_max_length:
    case cSetting_gradient_min_length:
    case cSetting_gradient_normal_min_dot:
    case cSetting_gradient_step_size:
    case cSetting_gradient_min_slope:
    case cSetting_gradient_symmetry:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_min_mesh_spacing:
    case cSetting_mesh_grid_max:
    case cSetting_mesh_cutoff:
    case cSetting_mesh_carve_state:
    case cSetting_mesh_carve_cutoff:
    case cSetting_mesh_carve_selection:
    case cSetting_mesh_clear_state:
    case cSetting_mesh_clear_cutoff:
    case cSetting_mesh_clear_selection:
    case cSetting_mesh_mode:
    case cSetting_mesh_type:
    case cSetting_mesh_solvent:
    case cSetting_mesh_quality:
    case cSetting_mesh_skip:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_valence:
    case cSetting_valence_mode:
    case cSetting_valence_size:
    case cSetting_valence_zero_mode:
    case cSetting_valence_zero_scale:
    case cSetting_half_bonds:
    case cSetting_line_stick_helper:
    case cSetting_hide_long_bonds:
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_stick_transparency:
      /* line_stick_helper reads this too, and the line rep was never
         invalidated for it -- a pre-existing hole. */
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_stick_debug:
    case cSetting_stick_round_nub:
    case cSetting_stick_as_cylinders:
    case cSetting_stick_good_geometry:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_line_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_ribbon_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_dot_as_spheres:
      SceneInvalidate(G);
      SceneChanged(G);
      break;
    case cSetting_dot_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_nonbonded_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_nb_spheres_size:
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_nb_spheres_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        SceneChanged(G);
      }
      break;
    case cSetting_render_as_cylinders:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
    case cSetting_mesh_as_cylinders:
    case cSetting_line_as_cylinders:
    case cSetting_ribbon_as_cylinders:
    case cSetting_dash_as_cylinders:
      if (index == cSetting_dash_as_cylinders) {
        ExecutiveInvalidateRep(G, inv_sele, cRepAngle, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepDihedral, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepDash, cRepInvRep);
      }
    case cSetting_nonbonded_as_cylinders:
    case cSetting_alignment_as_cylinders:
    case cSetting_cartoon_nucleic_acid_as_cylinders:
      if (SettingGetGlobal_b(G, cSetting_render_as_cylinders)) {
        if (G->ShaderMgr->shaders_present &&
            !G->ShaderMgr->ShaderPrgExists("cylinder")) {
          SettingSet_b(G->Setting, cSetting_render_as_cylinders, 0);
          if (!quiet) {
            PRINTFB(G, FB_Setting, FB_Warnings)
            "Setting-Error: render_as_cylinders cannot be set when the "
            "Cylinder Shader is not available, setting render_as_cylinder back "
            "to false\n" ENDFB(G);
          }
          return;
        }
        switch (index) {
          case cSetting_render_as_cylinders:
          case cSetting_cartoon_nucleic_acid_as_cylinders:
          case cSetting_use_shaders:
            if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
                SettingGetGlobal_b(G, cSetting_render_as_cylinders)) {
              ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
            }
        }
        SceneChanged(G);
      }
      break;
    case cSetting_mesh_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        SceneChanged(G);
      }
      break;
    case cSetting_slice_height_scale:
    case cSetting_slice_height_map:
    case cSetting_slice_grid:
    case cSetting_slice_dynamic_grid:
    case cSetting_slice_dynamic_grid_resolution:
      ExecutiveInvalidateRep(G, inv_sele, cRepSlice, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_label_font_id:
    case cSetting_label_size:
      ExecutiveInvalidateRep(G, inv_sele, cRepLabel, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_retain_order:
    case cSetting_pdb_hetatm_sort:
    case cSetting_pdb_insertions_go_first:
      ExecutiveSort(G, inv_sele);
      break;
    case cSetting_roving_lines:
    case cSetting_roving_sticks:
    case cSetting_roving_spheres:
    case cSetting_roving_labels:
    case cSetting_roving_selection:
    case cSetting_roving_ribbon:
    case cSetting_roving_cartoon:
    case cSetting_roving_polar_contacts:
    case cSetting_roving_polar_cutoff:
    case cSetting_roving_nonbonded:
    case cSetting_roving_nb_spheres:
    case cSetting_roving_map1_level:
    case cSetting_roving_map2_level:
    case cSetting_roving_map3_level:
    case cSetting_roving_map1_name:
    case cSetting_roving_map2_name:
    case cSetting_roving_map3_name:
    case cSetting_roving_isosurface:
    case cSetting_roving_isomesh:
      SceneRovingChanged(G);
      break;
    case cSetting_roving_byres:
    case cSetting_roving_detail:
      SceneRovingDirty(G);
      break;
    case cSetting_dash_transparency:
    case cSetting_dash_length:
    case cSetting_dash_gap:
    case cSetting_dash_radius:
    case cSetting_dash_width:
    case cSetting_angle_size:
    case cSetting_label_digits:
    case cSetting_label_distance_digits:
    case cSetting_label_angle_digits:
    case cSetting_label_dihedral_digits:
    case cSetting_angle_label_position:
    case cSetting_dihedral_size:
    case cSetting_dihedral_label_position:
      ExecutiveRebuildAllObjectDist(G);
      SceneChanged(G);
      break;
    case cSetting_button_mode:
      EditorMouseInvalid(G);
      OrthoDirty(G);
      break;
    case cSetting_stick_radius:
    case cSetting_stick_h_scale:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon,
                             cRepInvRep); /* base width */
      SceneChanged(G);
      break;
    case cSetting_nb_spheres_quality: {
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
      SceneChanged(G);
    }
    case cSetting_cgo_sphere_quality:
    case cSetting_cgo_debug: {
      ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
      SceneChanged(G);
      break;
    }
    case cSetting_stick_quality:
    case cSetting_stick_ball:
    case cSetting_stick_nub:
    case cSetting_stick_ball_ratio:
    case cSetting_stick_ball_color:
    case cSetting_stick_fixed_radius:
    case cSetting_stick_valence_scale:
    case cSetting_stick_overlap:
    case cSetting_stick_color:
    case cSetting_stick_use_shader:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_clamp_colors:
    case cSetting_ramp_blend_nearby_colors:
      ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_label_color:
    case cSetting_label_outline_color:
    case cSetting_label_position:
      ExecutiveRebuildAllObjectDist(G);
      ExecutiveInvalidateRep(G, inv_sele, cRepLabel, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_cartoon_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_ribbon_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_cgo_line_width:
    case cSetting_line_width: {
      GLfloat* range = G->ShaderMgr->GetLineWidthRange();
      float line_width = SettingGetGlobal_f(G, index);
      {
        if (line_width <= 0.f) {
          PRINTFB(G, FB_Setting, FB_Warnings)
          " Setting-Warning: %s is set incorrectly (%f), setting to 1\n",
              rec.name, line_width ENDFB(G);
          SettingSetGlobal_f(G, index, 1.f);
        } else if (G->HaveGUI && range[1] > 0.f && line_width > range[1]) {
          PRINTFB(G, FB_Setting, FB_Warnings)
          " Setting-Warning: %s is out of range of the graphics card's "
          "capability (range: %f-%f), lines might not be rendered correctly\n",
              rec.name, range[0], range[1] ENDFB(G);
        }
      }
    }
    case cSetting_line_color:
    case cSetting_line_radius:
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_scenes_changed: {
      int scene_buttons = SettingGetGlobal_i(G, cSetting_scene_buttons);
      if (scene_buttons) OrthoInvalidateDoDraw(G);
    } break;
    case cSetting_button_mode_name: {
      auto internal_gui_mode =
          SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting);
      if (internal_gui_mode != InternalGUIMode::Default) {
        OrthoInvalidateDoDraw(G);
      }
    } break;
    case cSetting_scene_current_name:
      SceneRestartFrameTimer(G);
      break;
    case cSetting_mesh_width:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_ellipsoid_probability:
    case cSetting_ellipsoid_scale:
    case cSetting_ellipsoid_color:
    case cSetting_ellipsoid_transparency:
      ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_ellipsoid_quality:
    case cSetting_cgo_ellipsoid_quality:
      ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_mesh_color:
    case cSetting_mesh_negative_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_ray_color_ramps:
      ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_sphere_mode:
      if (SettingGet<int>(G, cSetting_sphere_mode) > 9) {
        PRINTFB(G, FB_Setting, FB_Warnings)
        " Setting-Warning: sphere_mode > 9 is not supported in Open-Source "
        "version of PyMOL\n" ENDFB(G);
      }
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      SceneInvalidate(G);
      break;
    case cSetting_cell_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepCell, cRepInvColor);
      SceneInvalidate(G);
    case cSetting_cull_spheres:
    case cSetting_sphere_scale:
    case cSetting_sphere_transparency:
    case cSetting_sphere_solvent:
    case cSetting_sphere_point_max_size:
    case cSetting_sphere_point_size:
    case cSetting_sphere_use_shader:
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      SceneInvalidate(G);
      break;
    case cSetting_sphere_quality:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_nonbonded_size:
    case cSetting_nonbonded_transparency:
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbondedSphere, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_mesh_radius:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_ambient_occlusion_scale:
      SceneChanged(G);
      break;
    case cSetting_ambient_occlusion_mode:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
      switch (SettingGetGlobal_i(G, cSetting_ambient_occlusion_mode) % 4) {
        case 1:
          SettingSetGlobal_i(G, cSetting_ambient_occlusion_smooth, 10);
          SettingSetGlobal_f(G, cSetting_ambient_occlusion_scale, 25.f);
          break;
        case 2:
          SettingSetGlobal_i(G, cSetting_ambient_occlusion_smooth, 5);
          SettingSetGlobal_f(G, cSetting_ambient_occlusion_scale, .9f);
          break;
        case 3:
          SettingSetGlobal_i(G, cSetting_ambient_occlusion_smooth, 5);
          SettingSetGlobal_f(G, cSetting_ambient_occlusion_scale, 1.1f);
          break;
      }
      SceneChanged(G);
      break;
    case cSetting_ambient_occlusion_smooth:
    case cSetting_surface_negative_color:
    case cSetting_surface_color:
    case cSetting_surface_ramp_above_mode:
    case cSetting_transparency:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_dot_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_sphere_color:
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvColor);
      SceneChanged(G);
      break;

    case cSetting_surface_quality:
    case cSetting_surface_mode:
    case cSetting_surface_normal:
    case cSetting_surface_type:
    case cSetting_surface_carve_state:
    case cSetting_surface_carve_cutoff:
    case cSetting_surface_carve_selection:
    case cSetting_surface_carve_normal_cutoff:
    case cSetting_surface_clear_state:
    case cSetting_surface_clear_cutoff:
    case cSetting_surface_clear_selection:
    case cSetting_surface_trim_cutoff:
    case cSetting_surface_trim_factor:
    case cSetting_surface_circumscribe:
    case cSetting_surface_solvent:
    case cSetting_surface_proximity:
    case cSetting_surface_cavity_mode:
    case cSetting_surface_cavity_radius:
    case cSetting_surface_cavity_cutoff:
    case cSetting_cavity_cull:
    case cSetting_surface_smooth_edges:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_surface_use_shader:
      SceneChanged(G);
      break;
    case cSetting_isosurface_algorithm:
    case cSetting_surface_negative_visible:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_mesh_negative_visible:
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvAll);
      SceneChanged(G);
      break;
    case cSetting_solvent_radius:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepMesh, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_trace_atoms_mode:
      ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_ribbon_smooth:
    case cSetting_ribbon_power:
    case cSetting_ribbon_power_b:
    case cSetting_ribbon_sampling:
    case cSetting_ribbon_radius:
    case cSetting_ribbon_width:
    case cSetting_ribbon_throw:
    case cSetting_ribbon_trace_atoms:
    case cSetting_ribbon_transparency:
      ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_draw_mode:
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      break;
    case cSetting_cartoon_side_chain_helper:
    case cSetting_cartoon_nucleic_acid_mode:
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_ribbon_side_chain_helper:
    case cSetting_ribbon_nucleic_acid_mode:
      ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepEllipsoid, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_ray_trace_mode: /* affects loop quality */
    case cSetting_cartoon_transparency:
    case cSetting_cartoon_ring_transparency:
    case cSetting_cartoon_trace_atoms:
    case cSetting_cartoon_refine:
    case cSetting_cartoon_nucleic_acid_color:
    case cSetting_cartoon_ring_mode:
    case cSetting_cartoon_ring_finder:
    case cSetting_cartoon_ring_width:
    case cSetting_cartoon_ring_color:
    case cSetting_cartoon_ladder_mode:
    case cSetting_cartoon_ladder_radius:
    case cSetting_cartoon_ladder_color:
    case cSetting_cartoon_sampling:
    case cSetting_cartoon_loop_quality:
    case cSetting_cartoon_loop_radius:
    case cSetting_cartoon_loop_cap:
    case cSetting_cartoon_tube_quality:
    case cSetting_cartoon_tube_radius:
    case cSetting_cartoon_tube_cap:
    case cSetting_cartoon_putty_quality:
    case cSetting_cartoon_putty_radius:
    case cSetting_cartoon_putty_range:
    case cSetting_cartoon_putty_scale_min:
    case cSetting_cartoon_putty_scale_max:
    case cSetting_cartoon_putty_scale_power:
    case cSetting_cartoon_putty_transform:
    case cSetting_cartoon_power:
    case cSetting_cartoon_power_b:
    case cSetting_cartoon_ring_radius:
    case cSetting_cartoon_rect_length:
    case cSetting_cartoon_rect_width:
    case cSetting_cartoon_oval_length:
    case cSetting_cartoon_oval_width:
    case cSetting_cartoon_oval_quality:
    case cSetting_cartoon_round_helices:
    case cSetting_cartoon_flat_sheets:
    case cSetting_cartoon_spline:
    case cSetting_cartoon_refine_normals:
    case cSetting_cartoon_smooth_cylinder_cycles:
    case cSetting_cartoon_smooth_cylinder_window:
    case cSetting_cartoon_smooth_loops:
    case cSetting_cartoon_dumbbell_width:
    case cSetting_cartoon_dumbbell_length:
    case cSetting_cartoon_dumbbell_radius:
    case cSetting_cartoon_fancy_helices:
    case cSetting_cartoon_fancy_sheets:
    case cSetting_cartoon_cylindrical_helices:
    case cSetting_cartoon_refine_tips:
    case cSetting_cartoon_helix_radius:
    case cSetting_cartoon_throw:
    case cSetting_cartoon_debug:
    case cSetting_cartoon_highlight_color:
    case cSetting_cartoon_discrete_colors:
    case cSetting_cartoon_smooth_first:
    case cSetting_cartoon_smooth_last:
    case cSetting_cartoon_smooth_cycles:
    case cSetting_cartoon_flat_cycles:
    case cSetting_cartoon_gap_cutoff:
    case cSetting_cartoon_all_alt:
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_cartoon_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_cgo_use_shader:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_cgo_shader_ub_flags:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_cgo_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_cgo_shader_ub_color:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_cgo_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_cgo_shader_ub_normal:
      if (SettingGetGlobal_b(G, cSetting_use_shaders) &&
          SettingGetGlobal_b(G, cSetting_cgo_use_shader)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvRep);
        ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
        // if spheres rendered with geometry, then normals are used
        // this should really only invalidate spheres that use normals
        // (sphere_mode=0, not sure about other modes) but invalidating all
        // spheres for now
        ExecutiveInvalidateRep(G, inv_sele, cRepSphere, cRepInvRep);
        SceneChanged(G);
      }
      break;
    case cSetting_trilines:
      ExecutiveInvalidateRep(G, inv_sele, cRepLine, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCyl, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepRibbon, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepCartoon, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepDash, cRepInvRep);
      ExecutiveInvalidateRep(G, inv_sele, cRepNonbonded, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_cell_centered:
      ExecutiveInvalidateRep(G, inv_sele, cRepCell, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_dot_width:
    case cSetting_dot_radius:
    case cSetting_dot_density:
    case cSetting_dot_mode:
    case cSetting_dot_solvent:
    case cSetting_dot_hydrogens:
    case cSetting_trim_dots:
      ExecutiveInvalidateRep(G, inv_sele, cRepDot, cRepInvRep);
      SceneChanged(G);
      break;
    case cSetting_bg_gradient:
      ColorUpdateFrontFromSettings(G);
      ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
      OrthoBackgroundTextureNeedsUpdate(G);
      SceneChanged(G);
      break;
    case cSetting_bg_image_mode:
    case cSetting_bg_image_filename:
    case cSetting_bg_image_linear:
    case cSetting_bg_image_tilesize:
      OrthoBackgroundTextureNeedsUpdate(G);
      SceneChanged(G);
      break;
    case cSetting_bg_rgb_top:
    case cSetting_bg_rgb_bottom: {
      /* clamp this value */
      const char* bg_image_filename =
          SettingGet_s(G, nullptr, nullptr, cSetting_bg_image_filename);
      if (!(bg_image_filename && bg_image_filename[0]) &&
          SettingGetGlobal_b(G, cSetting_bg_gradient) &&
          !OrthoBackgroundDataIsSet(*G->Ortho)) {
        ColorUpdateFrontFromSettings(G);
        ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
        OrthoBackgroundTextureNeedsUpdate(G);
      }
    }
      SceneChanged(G);
      break;
    case cSetting_bg_rgb: {
      /* clamp this value */
      const float* v =
          ColorGet(G, SettingGet_color(G, nullptr, nullptr, cSetting_bg_rgb));
      {
        const char* bg_image_filename =
            SettingGet_s(G, nullptr, nullptr, cSetting_bg_image_filename);
        if (!(bg_image_filename && bg_image_filename[0]) &&
            !OrthoBackgroundDataIsSet(*G->Ortho)) {
          ColorUpdateFront(G, v);
          ExecutiveInvalidateRep(G, inv_sele, cRepAll, cRepInvColor);
        }
      }
    }
      SceneChanged(G);
      break;
    case cSetting_selection_width:
    case cSetting_selection_width_scale:
    case cSetting_selection_width_max:
    case cSetting_selection_round_points:
      ExecutiveInvalidateSelectionIndicatorsCGO(G);
    case cSetting_line_smooth:
    case cSetting_ortho:
    case cSetting_chromadepth:
    case cSetting_transparency_mode:
      if (index == cSetting_transparency_mode)
        ExecutiveInvalidateRep(G, inv_sele, cRepCGO, cRepInvAll);
    case cSetting_depth_cue:
    case cSetting_fog:
    case cSetting_ray_transparency_oblique:
    case cSetting_ray_transparency_oblique_power:
      SceneInvalidate(G);
      break;
    case cSetting_sculpting:
      OrthoDirty(G);
      break;
    case cSetting_auto_overlay:
      OrthoRemoveAutoOverlay(G); /* always start clean */
      break;
    case cSetting_overlay:
    case cSetting_overlay_lines:
    case cSetting_text:
      OrthoDirty(G);
      break;
    case cSetting_internal_gui_mode:
    case cSetting_internal_gui_width:
    case cSetting_internal_gui:
    case cSetting_internal_feedback:
    case cSetting_mouse_grid:
    case cSetting_movie_panel_row_height:
    case cSetting_movie_panel:
      if (!SettingGetGlobal_b(G, cSetting_suspend_updates)) {
        OrthoCommandIn(G, "viewport");
      }
      break;
    case cSetting_suspend_updates:
      if (!SettingGetGlobal_b(G, cSetting_suspend_updates)) {
        SceneChanged(G); /* force big update upon resumption */
        OrthoDirty(G);
      }
      break;
    case cSetting_security:
      G->Security = SettingGetGlobal_i(G, cSetting_security);
      break;
    case cSetting_state:
      if (SettingGet<int>(G, index) < 0 /* all */) {
        PRINTFB(G, FB_Setting, FB_Warnings)
        " Setting-Warning: state can't be less than 0.\n" ENDFB(G);
      }
    case cSetting_all_states:
    case cSetting_frame:
      ExecutiveInvalidateSelectionIndicatorsCGO(G);
      SceneInvalidatePicking(G);
      SceneChanged(G);
      break;
    case cSetting_dynamic_width:
    case cSetting_dynamic_width_factor:
    case cSetting_dynamic_width_min:
    case cSetting_dynamic_width_max:
      SceneChanged(G);
      break;
    case cSetting_rock:
    case cSetting_sweep_mode:
    case cSetting_sweep_phase:
    case cSetting_sweep_angle:
    case cSetting_sweep_speed:
      SceneRestartSweepTimer(G);
      break;
    case cSetting_motion_power:
    case cSetting_motion_bias:
    case cSetting_motion_simple:
    case cSetting_motion_linear:
    case cSetting_motion_hand:
    case cSetting_movie_loop:
      if (SettingGetGlobal_b(G, cSetting_movie_auto_interpolate)) {
        ExecutiveMotionReinterpolate(G);
      }
      break;
    case cSetting_volume_bit_depth:
      ExecutiveInvalidateRep(G, inv_sele, cRepVolume, cRepInvAll);
      SceneInvalidate(G);
      break;
    case cSetting_volume_layers:
      ExecutiveInvalidateRep(G, inv_sele, cRepVolume, cRepInvColor);
      SceneInvalidate(G);
      break;
    case cSetting_cgo_transparency:
      SceneInvalidate(G);
      SceneChanged(G);
      break;
    case cSetting_label_connector_mode: {
      int lc_mode = SettingGetGlobal_i(G, cSetting_label_connector_mode);
      if (lc_mode < 0 || lc_mode > 4) {
        if (!quiet) {
          PRINTFB(G, FB_Setting, FB_Warnings)
          "Setting-Warning: label_connector_mode range = [0,4]" ENDFB(G);
        }
      }
    }
    case cSetting_float_labels:
    case cSetting_label_z_target:
    case cSetting_label_connector:
    case cSetting_label_connector_color:
    case cSetting_label_connector_width:
    case cSetting_label_connector_ext_length:
    case cSetting_label_bg_color:
    case cSetting_label_placement_offset:
    case cSetting_label_relative_mode:
    case cSetting_label_screen_point:
    case cSetting_label_multiline_spacing:
    case cSetting_label_multiline_justification:
    case cSetting_label_padding:
    case cSetting_label_bg_transparency:
    case cSetting_label_bg_outline:
    case cSetting_ray_label_connector_flat:
      ExecutiveInvalidateRep(G, inv_sele, cRepLabel, cRepInvAll);
      break;
    case cSetting_surface_color_smoothing:
    case cSetting_surface_color_smoothing_threshold:
      ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
      SceneChanged(G);
      break;
    case cSetting_smooth_half_bonds:
      SceneChanged(G);
      break;
    case cSetting_antialias_shader:
    case cSetting_atom_type_format:
    case cSetting_colored_feedback:
    case cSetting_suspend_undo:
    case cSetting_volume_mode:
      PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: %s is not supported in Open-Source version of PyMOL\n",
          rec.name ENDFB(G);
      break;
      break;
    case cSetting_surface_debug:
      if (SettingGetGlobal_b(G, cSetting_use_shaders)) {
        ExecutiveInvalidateRep(G, inv_sele, cRepSurface, cRepInvColor);
      }
      break;
    case cSetting_pick32bit:
      SceneInvalidatePicking(G);
      break;
    case cSetting_display_scale_factor: {
      int scaleFactor = SettingGetGlobal_i(G, cSetting_display_scale_factor);
      if (scaleFactor > 0) {
        _gScaleFactor = scaleFactor;
        ExecutiveInvalidateRep(G, nullptr, cRepLabel, cRepInvRep);
        OrthoCommandIn(G, "viewport");
      } else {
        SettingSetGlobal_i(G, cSetting_display_scale_factor, 1);
        PRINTFB(G, FB_Setting, FB_Warnings)
        "Setting-Error: Cannot set a display scale factor of 0\n" ENDFB(G);
      }
      break;
    }
    default:
      // A setting with no dedicated side effect above still changes what the
      // next frame looks like (e.g. the metal_* renderer knobs,
      // surface_clip_front / _back): mark the scene dirty and ask for a
      // redisplay so an on-demand draw loop repaints without waiting for a
      // pointer event. Reps are not rebuilt and the ray-traced copy is kept —
      // only a repaint is requested.
      SceneDirty(G);
      PyMOL_NeedRedisplay(G->PyMOL);
      break;
  }
}
