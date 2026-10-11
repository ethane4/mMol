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

// The object panel: drawing, clicks and drags, and the selection indicators.

#include "pymol/graphics/immediate_helper.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list.h"
#include "pymol/base/list_macros.h"
#include "pymol/ui/menu.h"
#include "pymol/scene/movie.h"
#include "pymol/objects/object_group.h"
#include "pymol/ui/seeker.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/text/text.h"
#include "pymol/graphics/texture.h"
#include "pymol/executive/executive_internal.h"

static void ExecutiveInvalidateGridSlots(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  I->ValidGridSlots = false;
}

static void ExecutiveSetGridSlot(SpecRec* rec, int new_grid_slot) {
  if (rec->grid_slot != new_grid_slot) {
    CGOFree(rec->gridSlotSelIndicatorsCGO);
    rec->gridSlotSelIndicatorsCGO = nullptr;
    rec->grid_slot = new_grid_slot;
  }
}
void ExecutiveUpdateGridSlots(PyMOLGlobals* G, int force) {
  CExecutive* I = G->Executive;
  int grid_slot_count = 0;
  int grid_by_group = 1; /* grid slots are inherited this many levels */

  ExecutiveUpdateGroups(G, false);
  if (force || (!I->ValidGridSlots)) {
    CTracker* I_Tracker = I->Tracker;
    I->ValidGridSlots = true;
    {
      SpecRec* rec = nullptr;
      while (ListIterate(I->Spec, rec, next)) {
        ExecutiveSetGridSlot(rec, 0);
        if (rec->type == cExecObject) {
          /* make sure every object (potentially) needing a grid slot gets one
           */
          switch (rec->obj->type) {
            case cObjectMolecule:
            case cObjectMap:
            case cObjectMesh:
            case cObjectMeasurement:
            case cObjectCallback:
            case cObjectCGO:
            case cObjectSurface:
            case cObjectSlice:
            case cObjectGadget:
            case cObjectGroup:
            case cObjectVolume:
              ExecutiveSetGridSlot(rec, ++grid_slot_count);
              break;
          }
        }
      }
    }

    if (grid_by_group) {
      SpecRec *rec = nullptr, *group_rec = nullptr;
      while (ListIterate(I->Spec, rec, next)) {
        OVreturn_word result;
        if (OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(
                                I->Lex, rec->group_name)))) {
          auto keyRes = I->Key.find(result.word);
          if (keyRes != I->Key.end()) {
            if (TrackerGetCandRef(I_Tracker, keyRes->second,
                                  (TrackerRef**)(void*)&group_rec)) {
              int grid_slot_group_depth = grid_by_group;
              {
                SpecRec* check_rec = group_rec;
                while (check_rec && grid_slot_group_depth) {
                  if (grid_slot_group_depth == 1)
                    ExecutiveSetGridSlot(rec, check_rec->grid_slot);
                  if (check_rec == rec) { /* cycle */
                    break;
                  } else {
                    check_rec = check_rec->group;
                    grid_slot_group_depth--;
                  }
                }
              }
            }
          }
        }
      }
    }

    {
      SpecRec* rec = nullptr;
      while (ListIterate(I->Spec, rec, next)) {
        if (rec->type == cExecObject) {
          int obj_slot = SettingGet_i(G, rec->obj->Setting.get(), nullptr,
                                      cSetting_grid_slot);
          if (obj_slot == -1) {
            rec->obj->grid_slot = rec->grid_slot;
          } else
            rec->obj->grid_slot = obj_slot;
        }
      }
    }
  }
}

void ExecutiveInvalidatePanelList(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  I->Panel.clear();
  ExecutiveInvalidateGridSlots(G);
}

/**
 * Add all members which belong to group
 */
static void PanelListGroup(CExecutive* I, SpecRec const* group, unsigned level,
                           bool hide_underscore) {
  for (auto& rec : pymol::make_list_adapter(I->Spec)) {
    if (rec.group != group) {
      continue;
    }

    assert(!rec.in_panel);

    if (rec.isHiddenNotRecursive(hide_underscore)) {
      continue;
    }

    if (!level) {
      assert(!rec.group_name[0]);
      rec.group_name[0] =
          0; /* force open any cycles which have been created... */
    }

    I->Panel.emplace_back(&rec, level);
    rec.in_panel = true;

    if (auto obj_group = dynamic_cast<ObjectGroup*>(rec.obj)) {
      auto& panelitem = I->Panel.back();
      panelitem.is_group = true;
      if (obj_group->OpenOrClosed) {
        panelitem.is_open = true;
        PanelListGroup(I, &rec, level + 1, hide_underscore);
      }
    }
  }
}

static void ExecutiveUpdatePanelList(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
  if (I->Panel.empty()) {
    for (auto& rec : pymol::make_list_adapter(I->Spec)) {
      rec.in_panel = false;
    }

    /* brute-force & inefficient -- need to optimize algorithm */
    PanelListGroup(I, nullptr, 0, hide_underscore);
  }
}

static int get_op_cnt(PyMOLGlobals* G) {
  int result = 5;
  if (!strcmp(SettingGetGlobal_s(G, cSetting_button_mode_name),
              "3-Button Motions"))
    result = 6;
  return result;
}

/**
 * Scroll the i'th match in the object menu panel to the top.
 * Scroll to last match if i < 0 and to first match if i > #-1.
 * Open groups if hit is inside.
 * Highlight the hit (same as mouse click highlight).
 *
 * Returns the number of hits
 */
int ExecutiveScrollTo(PyMOLGlobals* G, const char* name, int i) {
  CExecutive* I = G->Executive;
  int pos = 0, numhits = 0;
  ObjectGroup* group;
  SpecRec *tmp, *spec = nullptr, *first = nullptr;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  int j, lendiff, plen = strlen(name);

  ok_assert(1, I->Spec);

  // i'th substring match, skip the "all" item
  for (tmp = I->Spec->next; tmp; tmp = tmp->next) {
    lendiff = strlen(tmp->name) - plen;
    for (j = 0; j <= lendiff; j++)
      if (WordMatchNoWild(G, name, tmp->name + j, ignore_case)) {
        if (numhits++ == i || i < 0) spec = tmp;
        if (!first) first = tmp;
        break;
      }
    tmp->hilight = 0;
  }

  // if i was out of range
  if (!spec) spec = first;

  ok_assert(1, spec);

  // flash button until panel is clicked for the next time
  spec->hilight = 1;

  // open parent groups
  for (tmp = spec->group; tmp; tmp = tmp->group) {
    if (!(tmp->type == cExecObject && tmp->obj->type == cObjectGroup)) break;
    group = (ObjectGroup*)tmp->obj;
    if (!group->OpenOrClosed) {
      group->OpenOrClosed = 1;
      ExecutiveInvalidatePanelList(G);
    }
  }

  // in case any parent got opened
  ExecutiveUpdatePanelList(G);

  // scroll that record to the top
  for (auto const& panelitem : I->Panel) {
    if (panelitem.spec == spec) {
      I->m_ScrollBar.setValueNoCheck(pos);
      return numhits;
    }
    pos++;
  }

ok_except1:
  return numhits;
}

/* #define ExecLineHeight 18 */
#define ExecClickMargin DIP2PIXEL(2)
#define ExecTopMargin 0
#define ExecToggleMargin DIP2PIXEL(2)
#define ExecLeftMargin DIP2PIXEL(1)
#define ExecRightMargin 0
#define ExecToggleWidth DIP2PIXEL(17)
#define ExecToggleSize DIP2PIXEL(16)
#define ExecToggleTextShift DIP2PIXEL(4)

int ExecutiveSetDrag(PyMOLGlobals* G, const char* name, int quiet, int mode) {
  char drag_name[] = cEditorDrag;
  int set_flag = false;
  int need_sele = true;
  int result = true;
  if (name[0]) {
    pymol::CObject* obj = ExecutiveFindObjectByName(G, name);
    if (obj) {
      EditorSetDrag(G, obj, -1, quiet, SceneGetState(G));
      set_flag = true;
    } else {
      SpecRec* rec = ExecutiveFindSpec(G, name);
      if (rec) {
        if (rec->type == cExecSelection) {
          SelectorCreate(G, drag_name, name, nullptr, true, nullptr);
          need_sele = false;
          {
            int sele = SelectorIndexByName(G, drag_name);
            ObjectMolecule* objMol = SelectorGetSingleObjectMolecule(G, sele);
            if (objMol) {
              if (mode > 0) sele = -1; /* force drag by matrix */
              EditorSetDrag(G, objMol, sele, quiet, SceneGetState(G));
              set_flag = true;
            } else {
              PRINTFB(G, FB_Executive, FB_Errors)
              " Drag-Error: selection spans more than one object.\n" ENDFB(G);
            }
          }
        } else if (rec->type == cExecObject) {
          switch (rec->obj->type) {
            case cObjectGroup:
              PRINTFB(G, FB_Executive, FB_Errors)
              " Drag-Error: cannot drag group objects yet.\n" ENDFB(G);
              break;
          }
          result = false;
        }
      }
    }
    result = set_flag;
    if (!result) {
      EditorInactivate(G);
      PRINTFB(G, FB_Executive, FB_Errors)
      " Drag-Error: invalid or empty selection." ENDFB(G);
    } else if (EditorDraggingObjectMatrix(G)) {
      SelectorCreate(G, drag_name, "none", nullptr, true, nullptr);
    } else if (need_sele && (obj->type == cObjectMolecule) &&
               (!EditorDraggingObjectMatrix(G))) {
      SelectorCreate(G, drag_name, obj->Name, (ObjectMolecule*)(void*)obj, true,
                     nullptr); /* for indication only */
    }
  } else {
    EditorInactivate(G);
  }
  return result;
}

int ExecutivePop(PyMOLGlobals* G, const char* target, const char* source,
                 int quiet) {
  int ok = true;
  int src;
  int result = 0;

  ExecutiveDelete(G, target);
  if (ExecutiveFindObjectMoleculeByName(G, source)) {
    ok = false;
    PRINTFB(G, FB_Executive, FB_Errors)
    " Pop-Error: source selection '%s' can't be an object.\n", source ENDFB(G);

  } else {
    src = SelectorIndexByName(G, source);
    if (src < 0) ok = false;
    if (!ok) {
      PRINTFB(G, FB_Executive, FB_Errors)
      " Pop-Error: invalid source selection name '%s'\n", source ENDFB(G);
    } else {
      ObjectMoleculeOpRec op;

      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_Pop;
      SelectorCreateEmpty(G, target, true);
      op.i1 = SelectorIndexByName(G, target);
      op.i2 = 1;
      op.i3 = 0;
      ExecutiveObjMolSeleOp(G, src, &op);
      result = op.i3;
    }
  }
  if (!result) ExecutiveDelete(G, target);
  if (!ok)
    return -1;
  else
    return result;
}

#define ExecScrollBarMargin DIP2PIXEL(1)
#define ExecScrollBarWidth DIP2PIXEL(13)
void ExecutiveHideSelections(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;

  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecSelection) {
      if (rec->visible) {
        rec->visible = false;
        SceneInvalidate(G);
        SeqDirty(G);
        ReportEnabledChange(G, rec);
      }
    }
  }
}

void ExecutiveInvalidateSelectionIndicators(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
  I->selectorTextureSize = 0;
  // The sole caller is the glyph-atlas wrap in TextureIsCharTextured(), which
  // restarts the packing cursor at the origin. Forget where the indicator used
  // to live too, or ExecutiveRegenerateTextureForSelector() keeps re-filling
  // the stale sub-rect that the glyph packer is now handing out again — the
  // pink indicator and freshly rasterized label glyphs then overwrite each
  // other in the shared texture.
  I->selectorTextureAllocatedSize = 0;
}

void ExecutiveInvalidateSelectionIndicatorsCGO(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  if (I) {
    if (I->selIndicatorsCGO) {
      CGOFree(I->selIndicatorsCGO);
      I->selIndicatorsCGO = 0;
    }
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        CGOFree(rec->gridSlotSelIndicatorsCGO);
      }
    }
  }
}

static void ExecutiveRegenerateTextureForSelector(PyMOLGlobals* G,
                                                  int round_points,
                                                  int* widths_arg) {
  CExecutive* I = G->Executive;
  unsigned char* temp_buffer =
      pymol::malloc<unsigned char>(widths_arg[0] * widths_arg[0] * 4);
  int a, b;
  float mid_point, disty, distx, dist, wminusd;
  float widths[] = {widths_arg[0] / 2.f, widths_arg[1] / 2.f,
                    widths_arg[2] / 2.f};
  unsigned char* q;
  mid_point = ((widths_arg[0] - 1.f) / 2.f);
  TextureInitTextTexture(G);
  if (I->selectorTextureAllocatedSize < widths_arg[0]) {
    // need to re-allocate a new part of the texture for the selection indicator
    TextureGetPlacementForNewSubtexture(G, widths_arg[0], widths_arg[0],
                                        &I->selectorTexturePosX,
                                        &I->selectorTexturePosY);
    I->selectorTextureAllocatedSize = widths_arg[0];
  }
  q = temp_buffer;
  if (round_points) {
    for (b = 0; b < widths_arg[0]; b++) {
      disty = fabs(mid_point - b);
      for (a = 0; a < widths_arg[0]; a++) {
        distx = fabs(mid_point - a);
        dist = sqrt(distx * distx + disty * disty);
        wminusd = widths[0] - dist;
        if (dist < widths[2]) {
          // white
          q[0] = q[1] = q[2] = q[3] = 255;
        } else if (dist < widths[1]) {
          // black
          q[0] = q[1] = q[2] = 0;
          q[3] = 255;
        } else if (fabs(wminusd) < .5f) {
          // color plus transparency
          q[0] = 255;
          q[1] = 51;
          q[2] = 153;
          q[3] = (unsigned char)((wminusd + .5) * 255);
        } else if (dist < widths[0]) {
          // color (pink by default 1.f, .2f, .6f)
          q[0] = 255;
          q[1] = 51;
          q[2] = 153;
          q[3] = 255;
        } else {
          // black
          q[0] = q[1] = q[2] = q[3] = 0;
        }
        //	printf("%3d ", q[1]);
        q += 4;
      }
      //      printf("\n");
    }
  } else {
    for (b = 0; b < widths_arg[0]; b++) {
      dist = disty = fabs(mid_point - b);
      for (a = 0; a < widths_arg[0]; a++) {
        distx = fabs(mid_point - a);
        dist = disty;
        if (distx > dist) {
          dist = distx;
        }
        if (dist < widths[2]) {
          // white
          q[0] = q[1] = q[2] = q[3] = 255;
        } else if (dist < widths[1]) {
          // black
          q[0] = q[1] = q[2] = 0;
          q[3] = 255;
        } else {
          // color (pink by default 1.f, .2f, .6f)
          q[0] = 255;
          q[1] = 51;
          q[2] = 153;
          q[3] = 255;
        }
        //	printf("%3d ", q[1]);
        q += 4;
      }
      //      printf("\n");
    }
  }
  TextureFillNewSubtexture(G, widths_arg[0], widths_arg[0],
                           I->selectorTexturePosX, I->selectorTexturePosY,
                           temp_buffer);
  FreeP(temp_buffer);
}

/**
 * Draws the Selection Indicator CGO
 * @param selIndicatorsCGO CGO to draw
 */
static void ExecutiveRenderIndicatorCGOGLRaster(PyMOLGlobals* G,
                                                CGO* selIndicatorsCGO) {
  CExecutive* I = G->Executive;
  CShaderPrg* shaderPrg;
  float text_texture_dim = TextureGetTextTextureSize(G);
  float textureScale;
  int no_depth = (int)SettingGetGlobal_f(G, cSetting_selection_overlay);
  shaderPrg = G->ShaderMgr->Enable_IndicatorShader();
  if (!shaderPrg) return;
  glEnable(GL_POINT_SPRITE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  shaderPrg->SetLightingEnabled(0);
  shaderPrg->SetAttrib4fLocation("a_Color", 1.f, 1.f, 1.f, 1.f);
  shaderPrg->Set1f("g_pointSize", DIP2PIXEL(I->selectorTextureSize));
  shaderPrg->Set2f("textureLookup", I->selectorTexturePosX / text_texture_dim,
                   I->selectorTexturePosY / text_texture_dim);
  textureScale = I->selectorTextureSize / text_texture_dim;
  shaderPrg->Set2f("textureScale", textureScale, textureScale);
  int v[4];
  glGetIntegerv(GL_VIEWPORT, v);
  shaderPrg->Set4f("viewport", v[0], v[1], v[2], v[3]);
  if (no_depth) glDisable(GL_DEPTH_TEST);
  CGORender(selIndicatorsCGO, nullptr, nullptr, nullptr, nullptr, nullptr);
  if (no_depth) glEnable(GL_DEPTH_TEST);
  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);
  glDisable(GL_POINT_SPRITE);
  shaderPrg->Disable();
}

/**
 * Draws the Selection Indicator CGO for a given pass
 * @param rec SpecRec's whose selection color to use
 * @param pass Pass number
 * @param width Width of the selection indicator
 * @param selIndicatorsCGO CGO to draw into
 */
static void ExecutiveSetupIndicatorPassMultipassImmediate(PyMOLGlobals* G,
                                                          SpecRec* rec,
                                                          int pass, int width,
                                                          CGO* cgo) {}
/**
 * Gets the adjusted selection width
 * @param round_points rounded indicator points
 * @return adjusted selection width
 */
static int ExecutiveGetAdjustedSelectionWidth(PyMOLGlobals* G,
                                              bool round_points) {
  int width;
  float width_scale = SettingGet<float>(G, cSetting_selection_width_scale);
  auto min_width =
      static_cast<int>(SettingGet<float>(G, cSetting_selection_width));
  auto max_width =
      static_cast<int>(SettingGet<float>(G, cSetting_selection_width_max));

  if (width_scale >= 0.0f) {
    width = static_cast<int>(
        width_scale * std::fabs(SettingGetGlobal_f(G, cSetting_stick_radius)) /
        SceneGetScreenVertexScale(G, nullptr));
    width = std::clamp(width, min_width, max_width);
  }
  if (round_points) {
    width = static_cast<int>(width * 1.44f);
  }
  if (width > 6 && width & 0x1) { /* keep it even above 6 */
    width--;
  }
  return width;
}

/**
 * Gets the selection level widths
 * @param width selection width
 * @return level widths
 */
static glm::ivec3 ExecutiveGetSelectionLevelWidths(int width) {
  glm::ivec3 level_widths(width, 0, 0);
  if (width > 2) {
    switch (width) {
      case 1:
      case 2:
      case 3:
        level_widths[1] = 1;
        break;
      case 4:
        level_widths[1] = 2;
        break;
      case 5:
        level_widths[1] = 3;
        break;
      case 6:
      case 7:
      case 8:
      case 9:
        level_widths[1] = 4;
        break;
      default:
        level_widths[1] = 6;
        break;
    }
  }
  if (width > 4) {
    if (width > 5) {
      level_widths[2] = 2;
    } else {
      level_widths[2] = 1;
    }
  }
  return level_widths;
}

/**
 * Renders the selection indicators from the selected targets
 * @param use_shader use shader
 * @param width selection width
 * @param curState current state
 * @param slot grid slot
 * @param grid grid info
 */
static void ExecutiveRenderSelectionsFromTargets(PyMOLGlobals* G,
                                                 bool use_shader, int width,
                                                 int curState, int slot,
                                                 GridInfo* grid) {
  auto I = G->Executive;
  bool vis_only = SettingGet<bool>(G, cSetting_selection_visible_only);

  for (auto& rec : pymol::make_list_adapter(I->Spec)) {
    if (rec.type != cExecSelection) {
      continue;
    }

    if (!rec.visible) {
      continue;
    }
    int enabled = true;
    auto* group_rec = rec.group;
    while (enabled && group_rec) {
      if (!group_rec->visible)
        enabled = false;
      else
        group_rec = group_rec->group;
    }

    if (!enabled) {
      continue;
    }

    int sele = SelectorIndexByName(G, rec.name); /* TODO: speed this up */
    if (sele < 0) {
      continue;
    }
    int totpass = use_shader ? 1 : 3;

    if (use_shader) {  // there is only one pass anyway
      if (!slot) {
        I->selIndicatorsCGO = CGONew(G);
        CGODotwidth(I->selIndicatorsCGO, width);
        CGOBegin(I->selIndicatorsCGO, GL_POINTS);
      }
    } else {
      I->selIndicatorsCGO = CGONew(G);
    }

    for (int pass = 0; pass < totpass; pass++) {
      if (!use_shader) {
        ExecutiveSetupIndicatorPassMultipassImmediate(G, &rec, pass, width,
                                                      I->selIndicatorsCGO);
        CGOBegin(I->selIndicatorsCGO, GL_POINTS);
      }
      for (auto& rec1 : pymol::make_list_adapter(I->Spec)) {
        if (rec1.type != cExecObject) {
          continue;
        }
        if (rec1.obj->type != cObjectMolecule) {
          continue;
        }
        auto obj = static_cast<ObjectMolecule*>(rec1.obj);
        if (SceneGetDrawFlagGrid(G, grid, rec1.obj->grid_slot)) {
          if (!use_shader || !slot) {
            ObjectMoleculeRenderSele(obj, curState, sele, vis_only,
                                     *I->selIndicatorsCGO);
          } else if (!rec1.gridSlotSelIndicatorsCGO) {
            if (!use_shader) {
              ObjectMoleculeRenderSele(obj, curState, sele, vis_only,
                                       *I->selIndicatorsCGO);
            } else {
              auto gssCGO = std::make_unique<CGO>(G);
              CGODotwidth(gssCGO.get(), width);
              CGOBegin(gssCGO.get(), GL_POINTS);
              ObjectMoleculeRenderSele(obj, curState, sele, vis_only, *gssCGO);
              CGOEnd(gssCGO.get());
              CGOStop(gssCGO.get());
              CGO* optimized = CGOOptimizeToVBONotIndexedNoShader(gssCGO.get());
              rec1.gridSlotSelIndicatorsCGO = optimized;
              assert(rec1.gridSlotSelIndicatorsCGO);
              ExecutiveRenderIndicatorCGOGLRaster(
                  G, rec1.gridSlotSelIndicatorsCGO);
            }
          }
        } else {
          CGOFree(rec1.gridSlotSelIndicatorsCGO);
        }
      }
      if (!use_shader) {
        CGOEnd(I->selIndicatorsCGO);
      }
    }  // for pass
    if (use_shader) {
      if (!slot) {
        CGOEnd(I->selIndicatorsCGO);
        CGOStop(I->selIndicatorsCGO);
        CGO* optimized =
            CGOOptimizeToVBONotIndexedNoShader(I->selIndicatorsCGO);
        CGOFree(I->selIndicatorsCGO);
        if (I->selIndicatorsCGO = optimized) {
          assert(I->selIndicatorsCGO->use_shader);
          ExecutiveRenderSelections(G, curState, slot, grid);
        }
        return;
      }
    } else {
      CGOStop(I->selIndicatorsCGO);
      ExecutiveRenderSelections(G, curState, slot, grid);
    }
  }  // for rec
}

// TODO: This should be incorporated into a CGO or GL Pipeline object
class ExecutiveRenderSelectionSetTopLevelGLPipeline {
 public:
  ExecutiveRenderSelectionSetTopLevelGLPipeline(bool round_points,
                                                bool no_depth, bool fog)
      : m_round_points{round_points}, m_no_depth{no_depth}, m_fog{fog} {
    if (no_depth) glDisable(GL_DEPTH_TEST);
  }
  ~ExecutiveRenderSelectionSetTopLevelGLPipeline() {
    if (m_no_depth) glEnable(GL_DEPTH_TEST);
  }

 private:
  bool m_round_points;
  bool m_no_depth;
  bool m_fog;
  GLint alpha_func_before;
  GLfloat alpha_ref_before;
};

// Name of the transient hover-PREVIEW selection (modules/mmol/pick.py
// _PRESELECT). ExecutiveGetSelectionCoords skips it so the preview is never
// drawn in the committed selection color; SceneDrawMetalPreselection gathers it
// separately (via ExecutiveGetNamedSelectionCoords) in its own color/size.
static const char kPreselectionName[] = "_preselect";

// Gather the coordinates of one specific named selection's atoms (same
// grid-cell filter and current-state logic as ExecutiveGetSelectionCoords, but
// restricted to `name`). Used to draw the hover preview separately from the
// committed selection. Appends 3 floats (x,y,z) per matching atom to `coords`.
/* mMol: the selection markers of `sele` in `obj`, where it renders: enabled
   objects only (a selection spanning hidden copies drew floating markers),
   every state when all_states shows them all, with the object's TTT applied. */
static void AppendSelectionCoords(PyMOLGlobals* G, ObjectMolecule* obj,
                                  int sele, std::vector<float>& coords) {
  if (!obj->Enabled) return;
  // members first (one pass), so all_states doesn't redo the test per state
  std::vector<int> members;
  for (int atIdx = 0; atIdx < obj->NAtom; atIdx++)
    if (SelectorIsMember(G, obj->AtomInfo[atIdx].selEntry, sele))
      members.push_back(atIdx);
  if (members.empty()) return;
  const bool every = ObjectGetCurrentState(obj, false) == -1;
  int st0 = std::max(obj->getCurrentState(), 0), st1 = st0 + 1;
  if (every) {
    st0 = 0;
    st1 = obj->NCSet;
  }
  for (int st = st0; st < st1 && st < obj->NCSet; ++st) {
    auto* cs = obj->CSet[st];
    if (!cs) continue;
    for (int atIdx : members) {
      int coordIdx = cs->atmToIdx(atIdx);
      if (coordIdx < 0) continue;
      float v[3];
      copy3f(cs->coordPtr(coordIdx), v);
      if (obj->TTTFlag) transformTTT44f3f(obj->TTT, v, v);
      coords.insert(coords.end(), v, v + 3);
    }
  }
}

void ExecutiveGetNamedSelectionCoords(PyMOLGlobals* G, const char* name,
                                      std::vector<float>& coords) {
  auto I = G->Executive;

  int sele = SelectorIndexByName(G, name);
  if (sele < 0) return;

  for (auto& rec1 : pymol::make_list_adapter(I->Spec)) {
    if (rec1.type != cExecObject || rec1.obj->type != cObjectMolecule) continue;

    // Same grid-cell filter as ExecutiveGetSelectionCoords: when grid is
    // active, only collect coords for objects in the cell being rendered.
    if (!SceneGetDrawFlagGrid(G, &G->Scene->grid, rec1.obj->grid_slot))
      continue;

    AppendSelectionCoords(G, static_cast<ObjectMolecule*>(rec1.obj), sele,
                          coords);
  }
}

void ExecutiveGetSelectionCoords(PyMOLGlobals* G, std::vector<float>& coords) {
  auto I = G->Executive;
  bool vis_only = SettingGet<bool>(G, cSetting_selection_visible_only);

  for (auto& rec : pymol::make_list_adapter(I->Spec)) {
    if (rec.type != cExecSelection || !rec.visible) continue;

    // Skip the transient hover-preview selection — it is drawn separately (in a
    // distinct color/size) by SceneDrawMetalPreselection, so it must never be
    // painted in the committed selection color here.
    if (strcmp(rec.name, kPreselectionName) == 0) continue;

    int sele = SelectorIndexByName(G, rec.name);
    if (sele < 0) continue;

    for (auto& rec1 : pymol::make_list_adapter(I->Spec)) {
      if (rec1.type != cExecObject || rec1.obj->type != cObjectMolecule)
        continue;

      // In grid mode, only collect coords for objects belonging to the cell
      // currently being rendered (G->Scene->grid.slot). When grid is inactive
      // this returns true for every object, so the non-grid path is unchanged.
      if (!SceneGetDrawFlagGrid(G, &G->Scene->grid, rec1.obj->grid_slot))
        continue;

      AppendSelectionCoords(G, static_cast<ObjectMolecule*>(rec1.obj), sele,
                            coords);
    }
  }
}

void ExecutiveRenderSelections(PyMOLGlobals* G, int curState, int slot,
                               GridInfo* grid) {
  CExecutive* I = G->Executive;
  bool use_shader = SettingGet<bool>(G, cSetting_use_shaders);
  bool round_points = SettingGet<bool>(G, cSetting_selection_round_points);
  short inv_indicators = false;

  auto width = ExecutiveGetAdjustedSelectionWidth(G, round_points);

  if (I->selectorTextureSize != width || round_points != I->selectorIsRound) {
    inv_indicators = true;
  } else {
    if (use_shader) {
      if (slot) {
        for (auto& rec : pymol::make_list_adapter(I->Spec)) {
          if (rec.type == cExecObject) {
            if (SceneGetDrawFlagGrid(G, grid, rec.obj->grid_slot)) {
              if (rec.gridSlotSelIndicatorsCGO) {
                ExecutiveRenderIndicatorCGOGLRaster(
                    G, rec.gridSlotSelIndicatorsCGO);
              }
            }
          }
        }
      } else if (I->selIndicatorsCGO) {
        ExecutiveRenderIndicatorCGOGLRaster(G, I->selIndicatorsCGO);
        return;
      }
    } else {
      if (I->selIndicatorsCGO) {
        std::optional<ExecutiveRenderSelectionSetTopLevelGLPipeline>
            immediatePipeline;
        auto no_depth =
            static_cast<int>(SettingGet<float>(G, cSetting_selection_overlay));
        [[maybe_unused]] bool fog = SettingGet<bool>(G, cSetting_depth_cue) &&
                                    SettingGet<float>(G, cSetting_fog);
        immediatePipeline.emplace(round_points, no_depth, fog);
        CGORender(I->selIndicatorsCGO, nullptr, nullptr, nullptr, nullptr,
                  nullptr);
        return;
      }
    }
  }
  if (inv_indicators) {
    CGOFree(I->selIndicatorsCGO);
    if (slot) {
      for (auto& rec : pymol::make_list_adapter(I->Spec)) {
        if (rec.type == cExecObject) {
          if (SceneGetDrawFlagGrid(G, grid, rec.obj->grid_slot)) {
            CGOFree(rec.gridSlotSelIndicatorsCGO);
          }
        }
      }
    }
  }

  auto list_adaptor = pymol::make_list_adapter(I->Spec);
  bool any_active = std::any_of(
      std::begin(list_adaptor), std::end(list_adaptor), [](const auto& rec) {
        return rec.type == cExecSelection && rec.visible;
      });
  if (!any_active) {
    return;
  }

  if (use_shader && inv_indicators) {
    auto level_widths = ExecutiveGetSelectionLevelWidths(width);
    ExecutiveInvalidateSelectionIndicatorsCGO(G);
    ExecutiveRegenerateTextureForSelector(G, round_points,
                                          glm::value_ptr(level_widths));
  }
  I->selectorTextureSize = width;
  I->selectorIsRound = round_points;

  ExecutiveRenderSelectionsFromTargets(G, use_shader, width, curState, slot,
                                       grid);
}
int CExecutive::click(int button, int x, int y, int mod) {
  PyMOLGlobals* G = m_G;
  CExecutive* I = G->Executive;
  int n, a;
  SpecRec* rec = nullptr;
  int t, xx;
  int pass = false;
  int skip;
  int ExecLineHeight =
      DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  int op_cnt = get_op_cnt(G);

  if (y < I->HowFarDown) {
    if (SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) !=
        InternalGUIMode::Default)
      return SceneGetBlock(G)->click(button, x, y, mod);
  }

  switch (button) {
    case P_GLUT_BUTTON_SCROLL_FORWARD:
      I->m_ScrollBar.moveBy(-1);
      return 1;
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      I->m_ScrollBar.moveBy(1);
      return 1;
  }

  n = ((I->rect.top - y) - (ExecTopMargin + ExecClickMargin)) / ExecLineHeight;
  a = n;
  xx = (x - rect.left);
  if (I->ScrollBarActive) {
    if ((x - rect.left) <
        (ExecScrollBarWidth + ExecScrollBarMargin + ExecToggleMargin)) {
      pass = 1;
      I->m_ScrollBar.click(button, x, y, mod);
    }
    xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
  }
  skip = I->NSkip;
  if (!pass) {
    I->RecoverPressed = nullptr;
    /* while(ListIterate(I->Spec,rec,next)) { */
    for (auto& panelitem : I->Panel) {
      auto* const panel = &panelitem;
      rec = panel->spec;

      assert(rec->name[0] != '_' || !hide_underscore);
      {
        if (skip) {
          skip--;
        } else {
          if (!a) {
            t = ((rect.right - ExecRightMargin) - x - 1) / ExecToggleWidth;
            if (t < op_cnt) {
              int my = rect.top - (ExecTopMargin + n * ExecLineHeight) - 3;
              int mx = rect.right - (ExecRightMargin + t * ExecToggleWidth);

              char* namesele =
                  (rec->type == cExecObject) ? rec->obj->Name : rec->name;

              t = (op_cnt - t) - 1;
              switch (t) {
                case 0: /* action */
                  switch (rec->type) {
                    case cExecAll:
                      MenuActivate(G, mx, my, x, y, false, "all_action",
                                   namesele);
                      break;
                    case cExecSelection:
                      MenuActivate(G, mx, my, x, y, false, "sele_action",
                                   namesele);
                      break;
                    case cExecObject:
                      switch (rec->obj->type) {
                        case cObjectGroup:
                          MenuActivate(G, mx, my, x, y, false, "group_action",
                                       namesele);
                          break;
                        case cObjectMolecule:
                          MenuActivate(G, mx, my, x, y, false, "mol_action",
                                       namesele);
                          break;
                        case cObjectMap:
                          MenuActivate(G, mx, my, x, y, false, "map_action",
                                       namesele);
                          break;
                        case cObjectSurface:
                          MenuActivate(G, mx, my, x, y, false, "surface_action",
                                       namesele);
                          break;
                        case cObjectMesh:
                          MenuActivate(G, mx, my, x, y, false, "mesh_action",
                                       namesele);
                          break;
                        case cObjectMeasurement:
                        case cObjectCGO:
                        case cObjectCallback:
                        case cObjectAlignment:
                        case cObjectVolume:
                          MenuActivate(G, mx, my, x, y, false, "simple_action",
                                       namesele);
                          break;
                        case cObjectSlice:
                          MenuActivate(G, mx, my, x, y, false, "slice_action",
                                       namesele);
                          break;
                        case cObjectGadget:
                          MenuActivate(G, mx, my, x, y, false, "ramp_action",
                                       namesele);
                          break;
                      }
                      break;
                  }
                  break;
                case 1:
                  switch (rec->type) {
                    case cExecAll:
                      MenuActivate(G, mx, my, x, y, false, "mol_show",
                                   cKeywordAll);
                      break;
                    case cExecSelection:
                      MenuActivate(G, mx, my, x, y, false, "mol_show",
                                   namesele);
                      break;
                    case cExecObject:
                      switch (rec->obj->type) {
                        case cObjectGroup:
                        case cObjectMolecule:
                          MenuActivate(G, mx, my, x, y, false, "mol_show",
                                       namesele);
                          break;
                        case cObjectCGO:
                        case cObjectAlignment:
                          MenuActivate(G, mx, my, x, y, false, "cgo_show",
                                       namesele);
                          break;
                        case cObjectMeasurement:
                          MenuActivate(G, mx, my, x, y, false,
                                       "measurement_show", namesele);
                          break;
                        case cObjectMap:
                          MenuActivate(G, mx, my, x, y, false, "map_show",
                                       namesele);
                          break;
                        case cObjectMesh:
                          MenuActivate(G, mx, my, x, y, false, "mesh_show",
                                       namesele);
                          break;
                        case cObjectSurface:
                          MenuActivate(G, mx, my, x, y, false, "surface_show",
                                       namesele);
                          break;
                        case cObjectSlice:
                          MenuActivate(G, mx, my, x, y, false, "slice_show",
                                       namesele);
                          break;
                        case cObjectVolume:
                          MenuActivate(G, mx, my, x, y, false, "volume_show",
                                       namesele);
                          break;
                      }
                      break;
                  }
                  break;
                case 2:
                  switch (rec->type) {
                    case cExecAll:
                      MenuActivate(G, mx, my, x, y, false, "mol_hide",
                                   cKeywordAll);
                      break;
                    case cExecSelection:
                      MenuActivate(G, mx, my, x, y, false, "mol_hide",
                                   namesele);
                      break;
                    case cExecObject:
                      switch (rec->obj->type) {
                        case cObjectGroup:
                        case cObjectMolecule:
                          MenuActivate(G, mx, my, x, y, false, "mol_hide",
                                       namesele);
                          break;
                        case cObjectCGO:
                        case cObjectAlignment:
                          MenuActivate(G, mx, my, x, y, false, "cgo_hide",
                                       namesele);
                          break;
                        case cObjectMeasurement:
                          MenuActivate(G, mx, my, x, y, false,
                                       "measurement_hide", namesele);
                          break;
                        case cObjectMap:
                          MenuActivate(G, mx, my, x, y, false, "map_hide",
                                       namesele);
                          break;
                        case cObjectMesh:
                          MenuActivate(G, mx, my, x, y, false, "mesh_hide",
                                       namesele);
                          break;
                        case cObjectSurface:
                          MenuActivate(G, mx, my, x, y, false, "surface_hide",
                                       namesele);
                          break;
                        case cObjectSlice:
                          MenuActivate(G, mx, my, x, y, false, "slice_hide",
                                       namesele);
                          break;
                        case cObjectVolume:
                          MenuActivate(G, mx, my, x, y, false, "volume_hide",
                                       namesele);
                          break;
                      }
                      break;
                  }
                  break;
                case 3:
                  switch (rec->type) {
                    case cExecAll:
                      MenuActivate(G, mx, my, x, y, false, "mol_labels",
                                   "(all)");
                      break;
                    case cExecSelection:
                      MenuActivate(G, mx, my, x, y, false, "mol_labels",
                                   namesele);
                      break;
                    case cExecObject:
                      switch (rec->obj->type) {
                        case cObjectGroup:
                        case cObjectMolecule:
                          MenuActivate(G, mx, my, x, y, false, "mol_labels",
                                       namesele);
                          break;
                        case cObjectMeasurement:
                          break;
                        case cObjectMap:
                        case cObjectSurface:
                        case cObjectMesh:
                        case cObjectSlice:
                          break;
                      }
                      break;
                  }
                  break;
                case 4:
                  switch (rec->type) {
                    case cExecAll:
                    case cExecSelection:
                      MenuActivate(G, mx, my, x, y, false, "mol_color",
                                   namesele);
                      break;
                    case cExecObject:
                      switch (rec->obj->type) {
                        case cObjectGroup:
                        case cObjectMolecule:
                          MenuActivate(G, mx, my, x, y, false, "mol_color",
                                       namesele);
                          break;
                        case cObjectMap:
                        case cObjectCGO:
                        case cObjectAlignment:
                          MenuActivate(G, mx, my, x, y, false, "general_color",
                                       namesele);
                          break;
                        case cObjectMesh:
                          MenuActivate(G, mx, my, x, y, false, "mesh_color",
                                       namesele);
                          break;
                        case cObjectSurface:
                          MenuActivate2Arg(G, mx, my, x, y, false, "mesh_color",
                                           namesele, "surface");
                          break;
                        case cObjectMeasurement:
                          MenuActivate(G, mx, my, x, y, false,
                                       "measurement_color", namesele);
                          break;
                        case cObjectSlice:
                          MenuActivate(G, mx, my, x, y, false, "slice_color",
                                       namesele);
                          break;
                        case cObjectVolume:
                          MenuActivate(G, mx, my, x, y, false, "vol_color",
                                       namesele);
                          break;
                        case cObjectGadget:
                          MenuActivate(G, mx, my, x, y, false, "ramp_color",
                                       namesele);
                          break;
                      }
                      break;
                  }
                  break;
                case 5:
                  switch (rec->type) {
                    case cExecAll:
                      MenuActivate0Arg(G, mx, my, x, y, false, "camera_motion");
                      break;
                    case cExecSelection: /* not relevant at present */
                      break;
                    case cExecObject: /* this kinds of objects can be animated
                                       */
                      switch (rec->obj->type) {
                        case cObjectGroup: /* however, this just groups motions
                                              for a set of objects */
                        case cObjectMolecule:
                        case cObjectMeasurement:
                        case cObjectMap:
                        case cObjectSurface:
                        case cObjectCGO:
                        case cObjectMesh:
                          MenuActivate(G, mx, my, x, y, false, "obj_motion",
                                       namesele);
                          break;
                      }
                      break;
                  }
                  break;
              }
            } else { /* clicked in variable area */

              if (((panel->is_group) &&
                   (((xx)-1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                  ((!panel->is_group) &&
                   (((xx)-1) / DIP2PIXEL(8)) > panel->nest_level)) {
                /* clicked on name */

                rec->hilight = 1;
                switch (button) {
                  case P_GLUT_LEFT_BUTTON:
                    if (I->DragMode != ExecutiveDragMode::Off) {
                      break;
                    }
                    I->Pressed = n;
                    I->OldVisibility = rec->visible;
                    I->Over = n;
                    I->DragMode = ExecutiveDragMode::Visibility;
                    I->ToggleMode = ExecutiveToggleMode::DeferVisibility;
                    I->LastChanged = nullptr;
                    I->LastZoomed = nullptr;
                    if (mod == (cOrthoSHIFT | cOrthoCTRL)) {
                      I->ToggleMode = ExecutiveToggleMode::HoverActivate;
                      if (!rec->visible) {
                        ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility,
                                                   mod, false);
                      }
                      if (rec != I->LastZoomed)
                        ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false,
                                            -1.0F, true);
                      I->LastZoomed = rec;
                      I->LastChanged = rec;
                    } else if (mod & cOrthoSHIFT) {
                      ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod,
                                                 false);
                      I->ToggleMode = ExecutiveToggleMode::ImmediateVisibility;
                    } else if (mod & cOrthoCTRL) {
                      I->ToggleMode = ExecutiveToggleMode::HoverActivate;
                      if (!rec->visible) {
                        ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility,
                                                   mod, false);
                      }
                      I->LastChanged = rec;
                    }
                    I->PressedWhat = 1;
                    I->OverWhat = 1;
                    break;
                  case P_GLUT_MIDDLE_BUTTON:
                    if (I->DragMode != ExecutiveDragMode::Off) {
                      break;
                    }
                    I->Pressed = n;
                    I->OldVisibility = rec->visible;
                    I->Over = n;
                    I->LastOver = I->Over;
                    I->DragMode = ExecutiveDragMode::VisibilityWithCamera;
                    I->ToggleMode = ExecutiveToggleMode::DeferVisibility;
                    I->LastChanged = rec;
                    I->LastZoomed = nullptr;
                    if (mod & cOrthoCTRL) {
                      I->ToggleMode =
                          ExecutiveToggleMode::ZoomActivateDeactivatePrevious;
                      ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false, -1.0F,
                                          true);
                      I->LastZoomed = rec;
                      if (mod & cOrthoSHIFT) { /* exclusive */
                        I->ToggleMode =
                            ExecutiveToggleMode::ZoomExclusiveActivate;
                        ExecutiveSetObjVisib(G, cKeywordAll, false, false);
                        /* need to log this */
                        if (!rec->visible)
                          ExecutiveSpecSetVisibility(G, rec, true, 0, true);
                      }
                    } else {
                      I->ToggleMode =
                          ExecutiveToggleMode::CenterActivateDeactivatePrevious;
                      ExecutiveCenter(G, rec->name, -1, true, -1.0F, nullptr,
                                      true);
                    }
                    if (!rec->visible) {
                      ExecutiveSpecSetVisibility(G, rec, !rec->visible, mod,
                                                 false);
                      I->LastChanged = rec;
                    }
                    I->PressedWhat = 1;
                    I->OverWhat = 1;
                    break;
                  case P_GLUT_RIGHT_BUTTON:
                    if (I->DragMode != ExecutiveDragMode::Off) {
                      break;
                    }
                    I->DragMode = ExecutiveDragMode::Reorder; /* reorder */
                    I->Pressed = n;
                    I->Over = n;
                    I->LastOver = I->Over;
                    I->PressedWhat = 1;
                    I->OverWhat = 1;
                    break;
                }
                OrthoGrab(G, this);
                OrthoDirty(G);
              } else if (panel->is_group) {
                /* clicked on group control */
                rec->hilight = 2;
                I->DragMode = ExecutiveDragMode::Visibility;
                I->Pressed = n;
                I->Over = n;
                I->LastOver = I->Over;
                I->PressedWhat = 2;
                I->OverWhat = 2;

                OrthoGrab(G, this);
                OrthoDirty(G);
              }
            }
          }
          a--;
        }
      }
    }
  }
  PyMOL_NeedRedisplay(G->PyMOL);
  return (1);
}

int CExecutive::release(int button, int x, int y, int mod) {
  PyMOLGlobals* G = m_G;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  int pass = false;
  int skip;
  int xx;
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  if (y < I->HowFarDown) {
    if (SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) !=
        InternalGUIMode::Default)
      return SceneGetBlock(G)->release(button, x, y, mod);
  }

  xx = (x - rect.left);
  if (I->ScrollBarActive) {
    if ((x - rect.left) <
        (ExecScrollBarWidth + ExecScrollBarMargin + ExecToggleMargin)) {
      pass = 1;
      I->m_ScrollBar.release(button, x, y, mod);
      OrthoUngrab(G);
    }
    xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
  }

  skip = I->NSkip;

  if (!pass) {
    I->drag(x, y, mod); /* incorporate final changes in cursor position */
    switch (I->DragMode) {
      case ExecutiveDragMode::Visibility:

        /*while(ListIterate(I->Spec,rec,next)) { */
        for (auto& panelitem : I->Panel) {
          auto* const panel = &panelitem;
          rec = panel->spec;

          assert(rec->name[0] != '_' || !hide_underscore);
          {
            if (skip) {
              skip--;
            } else {
              if ((I->PressedWhat == 1) &&
                  (((panel->is_group) &&
                    ((xx - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                   ((!panel->is_group) &&
                    ((xx - 1) / DIP2PIXEL(8)) > panel->nest_level))) {
                /* over name */
                if (rec->hilight == 1) {
                  if (rec->type == cExecSelection) {
                    ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, 0,
                                               false);
                  } else {
                    ExecutiveSpecSetVisibility(G, rec, !I->OldVisibility, mod,
                                               true);
                  }
                }
              } else if ((I->PressedWhat == 2) && (panel->is_group)) {
                if (rec->hilight == 2) {
                  ObjectGroup* obj = (ObjectGroup*)rec->obj;
                  OrthoLineType buf2;
                  sprintf(buf2, "cmd.group(\"%s\",action='%s')\n",
                          rec->obj->Name, obj->OpenOrClosed ? "close" : "open");
                  PLog(G, buf2, cPLog_no_flush);
                  ExecutiveGroup(G, rec->obj->Name, "", 5, 1);
                }
                /* over group control */
              }
            }
          }
        }
        break;
      case ExecutiveDragMode::Reorder:
        if (I->ReorderFlag) {
          I->ReorderFlag = false;
          PLog(G, I->ReorderLog, cPLog_no_flush);
        }
        break;
    }
  }

  {
    SpecRec* rec = nullptr;
    while (ListIterate(I->Spec, rec, next)) {
      rec->hilight = 0;
    }
  }

  I->Over = -1;
  I->Pressed = -1;
  I->DragMode = ExecutiveDragMode::Off;
  I->PressedWhat = 0;
  OrthoUngrab(G);
  PyMOL_NeedRedisplay(G->PyMOL);
  return (1);
}
int CExecutive::drag(int x, int y, int mod) {
  PyMOLGlobals* G = m_G;
  CExecutive* I = G->Executive;
  int xx, t;
  int ExecLineHeight =
      DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
#ifndef NDEBUG
  int hide_underscore = SettingGetGlobal_b(G, cSetting_hide_underscore_names);
#endif
  int op_cnt = get_op_cnt(G);

  ExecutiveUpdateGroups(G, false);
  ExecutiveUpdatePanelList(G);
  if (y < I->HowFarDown) {
    if (SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) !=
        InternalGUIMode::Default)
      return SceneGetBlock(G)->drag(x, y, mod);
  }

  if (I->DragMode != ExecutiveDragMode::Off) {
    xx = (x - rect.left);
    t = ((rect.right - ExecRightMargin) - x) / ExecToggleWidth;
    if (I->ScrollBarActive) {
      xx -= (ExecScrollBarWidth + ExecScrollBarMargin);
    }

    {
      int row_offset;
      if ((xx >= 0) && (t >= op_cnt)) {
        row_offset = ((rect.top - y) - (ExecTopMargin + ExecClickMargin)) /
                     ExecLineHeight;
        I->Over = row_offset;
      } else {
        I->Over = -1;
        row_offset = -1;
        {
          SpecRec* rec = nullptr;
          while (ListIterate(I->Spec, rec, next)) rec->hilight = 0;
        }
      }

      if (I->RecoverPressed) {
        SpecRec* rec = nullptr;
        int skip = I->NSkip;
        int row = 0;
        for (auto& panelitem : I->Panel) {
          auto* const panel = &panelitem;
          rec = panel->spec;

          assert(rec->name[0] != '_' || !hide_underscore);
          {
            if (skip) {
              skip--;
            } else {
              if (rec == I->RecoverPressed) {
                I->Pressed = row;
                I->RecoverPressed = nullptr;
              }
              row++;
            }
          }
        }
      }

      if (I->PressedWhat == 2) {
        I->OverWhat = 0;
      }
      if (I->Over >= 0) {
        SpecRec* rec = nullptr;
        int skip = I->NSkip;
        int row = 0;
        switch (I->DragMode) {
          case ExecutiveDragMode::Visibility:

            for (auto& panelitem : I->Panel) {
              auto* const panel = &panelitem;
              rec = panel->spec;

              assert(rec->name[0] != '_' || !hide_underscore);
              {
                if (skip) {
                  skip--;
                } else {
                  rec->hilight = 0;
                  if ((I->PressedWhat == 1) && /* name button */
                      (((row >= I->Over) && (row <= I->Pressed)) ||
                       ((row >= I->Pressed) && (row <= I->Over)))) {
                    if (((panel->is_group) &&
                         ((xx - 1) / DIP2PIXEL(8)) > (panel->nest_level + 1)) ||
                        ((!panel->is_group) &&
                         ((xx - 1) / DIP2PIXEL(8)) > panel->nest_level)) {
                      /* dragged over name */

                      I->OverWhat = 1;

                      switch (I->ToggleMode) {
                        case ExecutiveToggleMode::DeferVisibility:
                          if (row || (row == I->Pressed)) rec->hilight = 1;
                          break;
                        case ExecutiveToggleMode::ImmediateVisibility:
                          if (row)
                            ExecutiveSpecSetVisibility(
                                G, rec, !I->OldVisibility, mod, false);
                          break;
                        case ExecutiveToggleMode::HoverActivate:
                          if ((row == I->Over) && row) {
                            if (I->LastChanged != rec) {
                              ExecutiveSpecSetVisibility(G, I->LastChanged,
                                                         false, mod, false);
                            }
                            if (!rec->visible) {
                              ExecutiveSpecSetVisibility(G, rec, true, mod,
                                                         false);
                              I->LastChanged = rec;
                            }
                            if (mod == (cOrthoSHIFT | cOrthoCTRL)) {
                              if (rec != I->LastZoomed)
                                ExecutiveWindowZoom(G, rec->name, 0.0F, -1,
                                                    false, -1.0F, true);
                              I->LastZoomed = rec;
                            }
                          }
                          break;
                      }
                    }
                  } else if ((row == I->Pressed) && (I->PressedWhat == 2)) {
                    if (!((panel->is_group) && ((xx - 1) / DIP2PIXEL(8)) >
                                                   (panel->nest_level + 1))) {
                      /* on group control */

                      I->OverWhat = 2;
                      if (I->PressedWhat == I->OverWhat) {
                        rec->hilight = 2;
                      }
                    } else {
                      I->OverWhat = 0;
                    }
                  }
                  row++;
                }
              }
            }
            break;
          case ExecutiveDragMode::Reorder: /* right buttonBROKEN */
          {
            if ((I->Over != I->Pressed) && (I->LastOver != I->Over)) {
              SpecRec* new_rec = nullptr;
              SpecRec* mov_rec = nullptr;

              for (auto& panelitem : I->Panel) {
                rec = panelitem.spec;
                assert(rec->name[0] != '_' || !hide_underscore);
                {
                  {
                    if (skip) {
                      skip--;
                    } else {
                      if (row == I->Pressed) {
                        mov_rec = rec;
                      }
                      if (row == I->Over) {
                        new_rec = rec;
                      }
                      row++;
                    }
                  }
                }
              }
              {
                int group_flag = false;
                int order_flag = false;
                int location = 0;
                char *first = nullptr, *second = nullptr;
                int is_child = false;

                if (mov_rec && (!new_rec) && (I->Over > I->Pressed) &&
                    mov_rec->group) {
                  first = mov_rec->group->name;
                  second = mov_rec->name;
                  order_flag = true;
                  strcpy(mov_rec->group_name, mov_rec->group->group_name);
                  group_flag = true;
                } else if (mov_rec && new_rec) {
                  if (new_rec->isChildOf(mov_rec)) {
                    /* do nothing when a group is dragged over one of its
                     * members
                     */
                  } else {
                    if (I->Pressed < I->Over) { /* downward */
                      first = new_rec->name;
                      second = mov_rec->name;
                      order_flag = true;
                    } else { /* upward */
                      first = mov_rec->name;
                      second = new_rec->name;
                      order_flag = true;
                      location = -2; /* upper */
                    }

                    if (mov_rec->group ==
                        new_rec->group) { /* reordering within a group level */
                      if ((new_rec->type == cExecObject) &&
                          (new_rec->obj->type == cObjectGroup)) {
                        ObjectGroup* group = (ObjectGroup*)new_rec->obj;
                        if (group->OpenOrClosed && !is_child) {
                          /* put inside an open group */
                          strcpy(mov_rec->group_name, new_rec->name);
                          order_flag = false;
                          group_flag = true;
                        }
                      }
                    } else if ((mov_rec->group != new_rec) &&
                               (new_rec->group != mov_rec)) {
                      if ((new_rec->type == cExecObject) &&
                          (new_rec->obj->type == cObjectGroup)) {
                        ObjectGroup* group = (ObjectGroup*)new_rec->obj;
                        if (group->OpenOrClosed && !is_child) {
                          /* put inside group */
                          strcpy(mov_rec->group_name, new_rec->name);
                        } else {
                          strcpy(mov_rec->group_name, new_rec->group_name);
                        }

                        /* WLD TEST */

                        if (I->Pressed < I->Over) { /* downward */
                          first = new_rec->name;
                          second = mov_rec->name;
                          order_flag = true;
                        } else { /* upward */
                          first = mov_rec->name;
                          second = new_rec->name;
                          order_flag = true;
                        }

                        /* WLD END */

                      } else
                        strcpy(mov_rec->group_name, new_rec->group_name);
                      group_flag = true;
                    }
                  }
                }

                if (group_flag) {
                  OrthoLineType buf;
                  if (mov_rec->group_name[0]) {
                    sprintf(buf, "group %s, %s\n", mov_rec->group_name,
                            mov_rec->name);
                  } else {
                    sprintf(buf, "ungroup %s\n", mov_rec->name);
                  }
                  PLog(G, buf, cPLog_no_flush);
                  ExecutiveInvalidateGroups(G, false);
                  I->RecoverPressed = mov_rec;
                  I->Pressed = 0;
                }
                if (order_flag && first && second) {
                  OrthoLineType order_input;
                  sprintf(order_input, "%s %s", first, second);
                  ExecutiveOrder(G, order_input, false, location);
                  sprintf(I->ReorderLog, "cmd.order(\"%s\",location=\"%s\")\n",
                          order_input, (location == -2) ? "upper" : "current");
                  PLog(G, I->ReorderLog, cPLog_no_flush);
                  I->RecoverPressed = mov_rec;
                  I->Pressed = 0;
                }
              }
            }
          } break;
          case ExecutiveDragMode::VisibilityWithCamera: /* middle button */
            for (auto& panelitem : I->Panel) {
              rec = panelitem.spec;
              assert(rec->name[0] != '_' || !hide_underscore);
              {
                if (skip) {
                  skip--;
                } else {
                  rec->hilight = 0;
                  if (((row >= I->Over) && (row <= I->Pressed)) ||
                      ((row >= I->Pressed) && (row <= I->Over))) {
                    switch (I->ToggleMode) {
                      case ExecutiveToggleMode::
                          CenterActivateDeactivatePrevious: /* center and
                                                               activate, while
                                                               deactivating
                                                               previous */
                        if ((row == I->Over) && row) {
                          if (I->LastChanged != rec) {
                            ExecutiveSpecSetVisibility(G, I->LastChanged, false,
                                                       mod, false);
                            ExecutiveCenter(G, rec->name, -1, true, -1.0F,
                                            nullptr, true);
                            if (!rec->visible)
                              ExecutiveSpecSetVisibility(G, rec, true, mod,
                                                         false);
                            I->LastChanged = rec;
                          }
                          rec->hilight = 0;
                        }
                        break;
                      case ExecutiveToggleMode::
                          ZoomActivateDeactivatePrevious: /* zoom and activate,
                                                             while deactivating
                                                             previous */
                        if ((row == I->Over) && row) {
                          if (I->LastChanged != rec) {
                            ExecutiveSpecSetVisibility(G, I->LastChanged, false,
                                                       mod, false);
                            ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false,
                                                -1.0F, true);
                            if (!rec->visible)
                              ExecutiveSpecSetVisibility(G, rec, true, mod,
                                                         false);
                            I->LastChanged = rec;
                          }
                          rec->hilight = 1;
                        }
                        break;
                      case ExecutiveToggleMode::
                          ZoomExclusiveActivate: /* zoom and make only object
                                                    enabled */
                        if ((row == I->Over) && row) {
                          if (rec != I->LastZoomed) {
                            ExecutiveSpecSetVisibility(G, I->LastZoomed, false,
                                                       mod, false);
                            ExecutiveWindowZoom(G, rec->name, 0.0F, -1, false,
                                                -1.0F, true);
                            I->LastZoomed = rec;
                            ExecutiveSpecSetVisibility(G, rec, true, 0, false);
                          }
                          rec->hilight = 1;
                        }
                    }
                  }
                  row++;
                }
              }
            }
            break;
        }
        I->LastOver = I->Over;

      } else if (I->LastChanged)
        ExecutiveSpecSetVisibility(G, I->LastChanged, false, mod, false);
      OrthoDirty(G);
    }
  }
  return (1);
}

static void draw_button(int x2, int y2, int w, int h, const float* light,
                        const float* dark, const float* inside, CGO* orthoCGO) {
  if (orthoCGO) {
    CGOColorv(orthoCGO, light);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2, y2, 0.f);
    CGOVertex(orthoCGO, x2, y2 + h, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2 + h, 0.f);
    CGOEnd(orthoCGO);
  } else {
    ImmBatch batch;
    batch.begin(GL_POLYGON);
    batch.color3fv(light);
    batch.vertex2i(x2, y2);
    batch.vertex2i(x2, y2 + h);
    batch.vertex2i(x2 + w, y2 + h);
    batch.vertex2i(x2 + w, y2);
    batch.end();
  }

  if (orthoCGO) {
    CGOColorv(orthoCGO, dark);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, x2 + 1, y2, 0.f);
    CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2, 0.f);
    CGOVertex(orthoCGO, x2 + w, y2 + h - 1, 0.f);
    CGOEnd(orthoCGO);
  } else {
    ImmBatch batch;
    batch.begin(GL_POLYGON);
    batch.color3fv(dark);
    batch.vertex2i(x2 + 1, y2);
    batch.vertex2i(x2 + 1, y2 + h - 1);
    batch.vertex2i(x2 + w, y2 + h - 1);
    batch.vertex2i(x2 + w, y2);
    batch.end();
  }

  if (inside) {
    if (orthoCGO) {
      CGOColorv(orthoCGO, inside);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      ImmBatch batch;
      batch.begin(GL_POLYGON);
      batch.color3fv(inside);
      batch.vertex2i(x2 + 1, y2 + 1);
      batch.vertex2i(x2 + 1, y2 + h - 1);
      batch.vertex2i(x2 + w - 1, y2 + h - 1);
      batch.vertex2i(x2 + w - 1, y2 + 1);
      batch.end();
    }
  } else { /* rainbow */
    if (orthoCGO) {
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOColor(orthoCGO, 0.1F, 1.0F, 0.1F);  // green
      CGOVertex(orthoCGO, x2 + 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.0F, 1.0F, 0.1F);  // yellow
      CGOVertex(orthoCGO, x2 + w - 1, y2 + h - 1, 0.f);
      CGOColor(orthoCGO, 1.f, 0.1f, 0.1f);  // red
      CGOVertex(orthoCGO, x2 + 1, y2 + 1, 0.f);
      CGOColor(orthoCGO, 0.1F, 0.1F, 1.0F);  // blue
      CGOVertex(orthoCGO, x2 + w - 1, y2 + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      ImmBatch batch;
      batch.begin(GL_POLYGON);
      batch.color3f(1.0F, 0.1F, 0.1F);
      batch.vertex2i(x2 + 1, y2 + 1);
      batch.color3f(0.1F, 1.0F, 0.1F);
      batch.vertex2i(x2 + 1, y2 + h - 1);
      batch.color3f(1.0F, 1.0F, 0.1F);
      batch.vertex2i(x2 + w - 1, y2 + h - 1);
      batch.color3f(0.1F, 0.1F, 1.0F);
      batch.vertex2i(x2 + w - 1, y2 + 1);
      batch.end();
    }
  }
}

static void draw_button_char(PyMOLGlobals* G, int x2, int y2, char ch,
                             CGO* orthoCGO) {
  TextSetColor3f(G, 0.0F, 0.0F, 0.0F);
  TextSetPos2i(G, x2 + ExecToggleTextShift, y2);
  TextDrawChar(G, ch, orthoCGO);
}

/**
 * Get the text color for the object name.
 * @param obj Object pointer, may be nullptr
 * @param mode 0: default, 1: first carbon atom color, 2: object color
 * @param bg_rgb Background color
 * @param default_rgb Default text color
 * @return pointer to borrowed RGB array
 */
static const float* getNameColor(const pymol::CObject* obj, int mode,
                                 const float* bg_rgb,
                                 const float* default_rgb) {
  enum {
    cNameColorMode_default = 0,
    cNameColorMode_carbon = 1,
    cNameColorMode_object = 2,
  };

  if (mode == cNameColorMode_default || !obj) {
    return default_rgb;
  }

  int color = cColorDefault;

  switch (mode) {
    case cNameColorMode_carbon:
      // First carbon atom color
      // Assuming that there is typically a carbon atom within the first few
      // atoms, this procedure should be cheap (O(1)).
      if (obj->type == cObjectMolecule) {
        auto objmol = static_cast<const ObjectMolecule*>(obj);
        for (auto ai = objmol->AtomInfo.data(), ai_end = ai + objmol->NAtom;
             ai != ai_end; ++ai) {
          if (ai->protons == cAN_C) {
            color = ai->color;
            break;
          }
        }
      }
      break;
    case cNameColorMode_object:
      // Object color
      color = obj->Color;
      break;
  }

  if (color != cColorDefault) {
    const float* rgb = ColorGet(obj->G, color);

    // only use if it's different from the background color
    if (!within3f(bg_rgb, rgb, 0.1f)) {
      return rgb;
    }
  }

  return default_rgb;
}
void CExecutive::draw(CGO* orthoCGO) {
  PyMOLGlobals* G = m_G;
  int x, y, xx, x2, y2;
  WordType ch;
  char* c = nullptr;
  float enabledColor[3] = {0.5F, 0.5F, 0.5F};
  float cloakedColor[3] = {0.35F, 0.35F, 0.35F};
  float pressedColor[3] = {0.7F, 0.7F, 0.7F};
  float disabledColor[3] = {0.25F, 0.25F, 0.25F};
  float lightEdge[3] = {0.6F, 0.6F, 0.6F};
  float darkEdge[3] = {0.35F, 0.35F, 0.35F};
  float captionColor[3] = {0.3F, 0.9F, 0.3F};

  float activeColor[3] = {0.9F, 0.9F, 1.0F};
  float toggleColor3[3] = {0.6F, 0.6F, 0.8F};

  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;
  int n_ent;
  int n_disp;
  int skip = 0;
  int row = -1;
  int ExecLineHeight =
      DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_control_size));
  int text_lift = (ExecLineHeight / 2) - DIP2PIXEL(5);
  int op_cnt = get_op_cnt(G);
  int full_names = SettingGetGlobal_b(G, cSetting_group_full_member_names);
  int arrows = SettingGetGlobal_b(G, cSetting_group_arrow_prefix);
  auto name_color_mode =
      SettingGet<int>(G, cSetting_internal_gui_name_color_mode);

  ExecutiveUpdatePanelList(G);

  /* if we're running with a GUI and have a valid panel */
  if (G->HaveGUI && G->ValidContext && ((rect.right - rect.left) > 6) &&
      !I->Panel.empty()) {
    int max_char;
    int nChar;

    /* count entries
     * do we have enough structures to warrant a scroll bar? */
    n_ent = 0;
    for (auto& panelitem : I->Panel) {
      rec = panelitem.spec;
      assert(rec && (rec->name[0] != '_' ||
                     !SettingGet<bool>(G, cSetting_hide_underscore_names)));
      n_ent++;
    }

    n_disp = ((rect.top - rect.bottom) - (ExecTopMargin)) / ExecLineHeight;
    if (n_disp < 1) n_disp = 1;

    /* we need a scrollbar */
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

    /* determination of longest string based on internal_gui_size, etc... */
    max_char =
        (((rect.right - rect.left) - (ExecLeftMargin + ExecRightMargin + 4)) -
         (op_cnt * ExecToggleWidth));
    if (I->ScrollBarActive) {
      max_char -= (ExecScrollBarMargin + ExecScrollBarWidth);
    }
    max_char /= DIP2PIXEL(8);

    /* fill and outline the entire block */
    if (SettingGet<InternalGUIMode>(cSetting_internal_gui_mode, G->Setting) ==
        InternalGUIMode::Default) {
      if (orthoCGO) CGOColorv(orthoCGO, BackColor);
      fill(orthoCGO);
      drawLeftEdge(orthoCGO);
    }

    /* draw the scroll bar */
    if (I->ScrollBarActive) {
      I->m_ScrollBar.setBox(
          rect.top - ExecScrollBarMargin, rect.left + ExecScrollBarMargin,
          rect.bottom + 2,
          rect.left + ExecScrollBarMargin + ExecScrollBarWidth);
      I->m_ScrollBar.draw(orthoCGO);
    }

    x = rect.left + ExecLeftMargin;
    y = (rect.top - ExecLineHeight) - ExecTopMargin;
    /*    xx = rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+op_cnt); */
    xx = rect.right - ExecRightMargin - ExecToggleWidth * (op_cnt);

    if (I->ScrollBarActive) {
      x += ExecScrollBarWidth + ExecScrollBarMargin;
    }
    skip = I->NSkip;

    /* for each object in the Panel... */
    for (auto& panelitem : I->Panel) {
      auto* const panel = &panelitem;
      rec = panel->spec;
      {
        if (skip) {
          skip--;
        } else {
          /* setup the X,Y offsets for this entry */
          row++;
          x2 = xx;
          y2 = y;
          nChar = max_char;

          if ((x - ExecToggleMargin) - (xx - ExecToggleMargin) > -10) {
            x2 = x + 10;
          }
          /*
           * The ASHLC menus; these access Python so,
           * protect this block from non-python instances
           */
          {
            int a;
            float toggleColor[3] = {0.5F, 0.5F, 1.0F};
            float toggleColor2[3] = {0.4F, 0.4F, 0.6F};
            float toggleDarkEdge[3] = {0.3F, 0.3F, 0.5F};
            float toggleLightEdge[3] = {0.7F, 0.7F, 0.9F};

            glColor3fv(toggleColor);
            for (a = 0; a < op_cnt; a++) {
              switch (a) {
                case 0:
                  /*
                     glColor3fv(toggleColor);
                     glBegin(GL_POLYGON);
                     glVertex2i(x2,y2+(ExecToggleSize)/2);
                     glVertex2i(x2+(ExecToggleSize)/2,y2);
                     glVertex2i(x2+ExecToggleSize,y2+(ExecToggleSize)/2);
                     glVertex2i(x2+(ExecToggleSize)/2,y2+ExecToggleSize);
                     glEnd();
                   */

                  /* the infamous ASHLC! */

                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, toggleColor,
                              orthoCGO);

                  draw_button_char(G, x2, y2 + text_lift, 'A', orthoCGO);
                  break;
                case 1:
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, toggleColor3,
                              orthoCGO);

                  draw_button_char(G, x2, y2 + text_lift, 'S', orthoCGO);
                  break;
                case 2:
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, toggleColor2,
                              orthoCGO);
                  draw_button_char(G, x2, y2 + text_lift, 'H', orthoCGO);
                  break;
                case 3:
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, toggleColor,
                              orthoCGO);
                  draw_button_char(G, x2, y2 + text_lift, 'L', orthoCGO);
                  break;
                case 4:
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, nullptr,
                              orthoCGO);
                  draw_button_char(G, x2, y2 + text_lift, 'C', orthoCGO);
                  break;
                case 5: {
                  float* button_color = toggleColor2;
                  int spec_level = 0;
                  /* choose color / brightness based on specification level */

                  if (rec->type == cExecAll) {
                    spec_level = MovieGetSpecLevel(G, SceneGetFrame(G));
                  } else if (rec->type == cExecObject) {
                    spec_level = ObjectGetSpecLevel(rec->obj, SceneGetFrame(G));
                  }
                  switch (spec_level) {
                    case 1:
                      button_color = toggleColor3;
                      break;
                    case 2:
                      button_color = activeColor;
                      break;
                  }
                  draw_button(x2, y2, ExecToggleSize, (ExecLineHeight - 1),
                              toggleLightEdge, toggleDarkEdge, button_color,
                              orthoCGO);
                }

                  draw_button_char(G, x2, y2 + text_lift, 'M', orthoCGO);
                  break;
              }
              x2 += ExecToggleWidth;
            }
          }
          /* end ASHLC */

          {
            int x3 = x;
            int hidden_prefix = false;

            TextSetColor(G, TextColor);
            TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);

            if ((rec->type == cExecObject) || (rec->type == cExecAll) ||
                (rec->type == cExecSelection)) {
              y2 = y;
              x2 = xx;
              if ((x - ExecToggleMargin) - (xx - ExecToggleMargin) >
                  DIP2PIXEL(-10)) {
                x2 = x + DIP2PIXEL(10);
              }
              x3 += panel->nest_level * DIP2PIXEL(8);
              TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
              nChar -= panel->nest_level;

              const float* but_color = disabledColor;
              {
                int but_width = (x2 - x3) - 1;

                /* drawing a group +/- NAME */
                if (panel->is_group) {
                  const int button_width = DIP2PIXEL(15);
                  if ((rec->hilight == 2) && (I->Over == I->Pressed)) {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1),
                                lightEdge, darkEdge, pressedColor, orthoCGO);
                  } else if (panel->is_open) {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1),
                                lightEdge, darkEdge, disabledColor, orthoCGO);
                  } else {
                    draw_button(x3, y2, button_width, (ExecLineHeight - 1),
                                lightEdge, darkEdge, disabledColor, orthoCGO);
                  }
                  TextSetPos2i(G, x3 + DIP2PIXEL(4), y2 + text_lift);
                  if (panel->is_open) {
                    TextDrawChar(G, '-', orthoCGO);
                  } else {
                    TextDrawChar(G, '+', orthoCGO);
                  }

                  but_width -= DIP2PIXEL(16);
                  x3 += DIP2PIXEL(16);
                  nChar -= 2;

                  TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
                }

                if ((rec->hilight == 1) ||
                    ((row == I->Over) && (I->OverWhat == 1))) {
                  /* button hull */
                  but_color = pressedColor;
                } else if (rec->visible) {
                  int enabled = true;
                  SpecRec* group_rec = rec->group;
                  while (enabled && group_rec) {
                    if (!group_rec->visible)
                      enabled = false;
                    else
                      group_rec = group_rec->group;
                  }

                  if (enabled) {
                    but_color = enabledColor;
                  } else {
                    but_color = cloakedColor;
                  }
                }

                draw_button(x3, y2, but_width, (ExecLineHeight - 1), lightEdge,
                            darkEdge, but_color, orthoCGO);
              }

              TextSetColor(G, getNameColor(rec->obj, name_color_mode, but_color,
                                           TextColor));

              /* object name */
              c = rec->name;

              /* parse out the prefix if group.foo */
              if (!full_names) {
                if (rec->group) { /* if prefix matches group name, then truncate
                                   */
                  char *p = c, *q = rec->group->name;
                  while ((*p == *q) && (*q)) {
                    p++;
                    q++;
                  }
                  if ((*p) && (!*q) && (*p == '.')) {
                    hidden_prefix = true;
                    c = p;
                  }
                }
              }

              /* wrap selection names with "(" and ")" */
              if (rec->type == cExecSelection)
                if ((nChar--) > 0) {
                  TextDrawChar(G, '(', orthoCGO);
                }
            }

            if (c) {
              if (hidden_prefix) {
                /* ^.name */
                if (arrows && ((nChar--) > 0)) {
                  TextDrawChar(G, '^', orthoCGO);
                  TextSetPos2i(G, x3 + DIP2PIXEL(2), y2 + text_lift);
                  TextDrawChar(G, '|', orthoCGO);
                }
              }

              /* draw the object name, char by char */
              while (*c) {
                if ((nChar--) > 0) {
                  TextDrawChar(G, *(c++), orthoCGO);
                } else
                  break;
              }
            }

            /* SELECTIONS: wrap selection names with "(" and ")" */
            if (rec->type == cExecSelection) {
              if ((nChar--) > 0) {
                TextDrawChar(G, ')', orthoCGO);
              }
            }

            /* OBJECTS: output any label captions, like state number, state
             * title */
            if (rec->type == cExecObject) {
              /* get this object's "caption" that goes on its title line,
               * currently, this is "state-title [curState/nState]" */
              c = rec->obj->getCaption(ch, WordLength);
              /* now print the caption */
              if (c && c[0] && nChar > 1) {
                TextSetColor(G, captionColor);
                TextSetPos2i(
                    G, x + DIP2PIXEL(2) + DIP2PIXEL(8) * (max_char - nChar),
                    y2 + text_lift);
                if ((nChar--) > 0) TextDrawChar(G, ' ', orthoCGO);
                while (*c && nChar > 0) {
                  /* allow color encoding for names */
                  if (TextSetColorFromCode(G, c, captionColor)) {
                    c += 4;
                  } else {
                    TextDrawChar(G, *(c++), orthoCGO);
                    --nChar;
                  }
                }
              }
            }
          }

          y -= ExecLineHeight;
          if (y < (rect.bottom)) break;
        }
      }
    }
    I->HowFarDown = y;
  }
}
void CExecutive::reshape(int width, int height) {
  PyMOLGlobals* G = m_G;
  CExecutive* I = G->Executive;

  Block::reshape(width, height);

  I->Width = rect.right - rect.left + 1;
  I->Height = rect.top - rect.bottom + 1;
}
