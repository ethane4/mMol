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
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
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

void CExecutive::reshape(int width, int height) {
  PyMOLGlobals* G = m_G;
  CExecutive* I = G->Executive;

  Block::reshape(width, height);

  I->Width = rect.right - rect.left + 1;
  I->Height = rect.top - rect.bottom + 1;
}
