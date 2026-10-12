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
#include "pymol/base/busy.h"

#include "pymol/ui/control.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/base/list.h"
#include "pymol/base/list_macros.h"
#include "pymol/scene/movie.h"
#include "pymol/image/my_png.h"
#include "pymol/objects/object_alignment.h"
#include "pymol/objects/object_dist.h"
#include "pymol/objects/object_gadget_ramp.h"
#include "pymol/objects/object_group.h"
#include "pymol/objects/object_mesh.h"
#include "pymol/objects/object_surface.h"
#include "pymol/objects/object_volume.h"
#include "pymol/base/pymol_options.h"
#include "pymol/ray/scene_ray.h"
#include "pymol/chem/sculpt_cache.h"
#include "pymol/ui/seeker.h"
#include "pymol/executive/setting_effects.h"
#include "pymol/executive/spec_rec.h"
#include "pymol/ui/wizard.h"

#include "ov/ov_context.h"

#include "pymol/scene/movie_scene.h"

#include "pymol/objects/property.h"

#include "pymol/executive/executive_internal.h"

pymol::TrackerAdapter<SpecRec> ExecutiveGetSpecRecsFromPattern(
    PyMOLGlobals* G, pymol::zstring_view str, bool allow_partial,
    bool expand_groups) {
  return pymol::TrackerAdapter<SpecRec>(
      G->Executive->Tracker, ExecutiveGetNamesListFromPattern(
                                 G, str.c_str(), allow_partial, expand_groups));
}

/**
 * ObjectIterator methods
 */
void ObjectIterator::reset() {
  rec = G->Executive->Spec;

  // DEBUG assume first element is always "all"
  if (rec->type != cExecAll) printf("Error: first SpecRec is not cExecAll\n");
}

bool ObjectIterator::next() {
  if (!rec || !(rec = rec->next)) return false;

  if (rec->type != cExecObject) return next();

  return true;
}

pymol::CObject* ObjectIterator::getObject() { return rec->obj; }

/**
 * True if `rec` and all its parent groups are enabled
 */
static bool SpecRecIsEnabled(const SpecRec* rec) {
  while (rec->visible && (rec = rec->group)) {
  }
  return !rec;
}

bool ExecutiveIsObjectType(const SpecRec& rec, int cObjectType) {
  return rec.type == cExecObject && rec.obj->type == cObjectType;
}
void ReportEnabledChange(PyMOLGlobals* G, SpecRec* rec) {
  OrthoInvalidateDoDraw(G);
  ExecutiveInvalidateSelectionIndicatorsCGO(G);
}

int ExecutiveReference(PyMOLGlobals* G, int action, const char* sele, int state,
                       int quiet) {
  int result = -1;
  int s1 = SelectorIndexByName(G, sele);
  if (s1 >= 0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);

    switch (action) {
      case 1:
        op.code = OMOP_ReferenceStore;
        break;
      case 2:
        op.code = OMOP_ReferenceRecall;
        break;
      case 3:
        op.code = OMOP_ReferenceValidate;
        break;
      case 4:
        op.code = OMOP_ReferenceSwap;
        break;
    }
    op.i1 = state;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G, s1, &op);
    result = op.i2;
  }
  return result;
}

std::string ExecutivePreparePseudoatomName(PyMOLGlobals* G,
                                           pymol::zstring_view object_name) {
  std::string new_object_name;
  if (object_name.empty()) {
    new_object_name = ExecutiveGetUnusedName(G, "pseudo");
  } else {
    ObjectNameType valid_name{};
    assert(object_name.size() < sizeof(ObjectNameType));
    std::copy_n(object_name.c_str(), object_name.size(), valid_name);
    ObjectMakeValidName(G, valid_name);
    new_object_name = valid_name;
  }
  return new_object_name;
}

pymol::Result<> ExecutivePseudoatom(
    PyMOLGlobals* G, pymol::zstring_view object_name_view, const char* sele,
    const char* name, const char* resn, const char* resi, const char* chain,
    const char* segi, const char* elem, float vdw, int hetatm, float b, float q,
    const char* label, const float* pos, int color, int state, int mode,
    int quiet) {
  pymol::Result<SelectorTmp> s1;

  auto object_name = object_name_view.c_str();
  auto obj = ExecutiveFindObject<ObjectMolecule>(G, object_name);

  int is_new = false;
  int sele_index = -1;
  float local_pos[3];

  if (sele && sele[0]) {
    if (WordMatchExact(G, cKeywordCenter, sele, true)) {
      SceneGetCenter(G, local_pos);
      pos = local_pos;
    } else if (WordMatchExact(G, cKeywordOrigin, sele, true)) {
      SceneOriginGet(G, local_pos);
      pos = local_pos;
    } else {
      s1 = SelectorTmp::make(G, sele);
      p_return_if_error(s1);
      sele_index = s1->getIndex();
      assert(sele_index >= 0);
    }
  }
  if (!obj) {
    /* new object */
    is_new = true;
    obj = new ObjectMolecule(G, false);
    ObjectSetName(obj, object_name);
  }

  if (ObjectMoleculeAddPseudoatom(obj, sele_index, name, resn, resi, chain,
                                  segi, elem, vdw, hetatm, b, q, label, pos,
                                  color, state, mode, quiet)) {
    if (is_new) {
      ExecutiveDelete(G, object_name); /* just in case */
      ExecutiveManageObject(G, obj, false, true);
    } else {
      ExecutiveUpdateObjectSelection(G, obj);
    }
  }
  return {};
}

void ExecutiveInvalidateSceneMembers(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  I->ValidSceneMembers = false;
}

void ExecutiveUpdateSceneMembers(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  ExecutiveUpdateGroups(G, false);
  ExecutiveUpdateGridSlots(G, false);
  if (!I->ValidSceneMembers) {
    SpecRec* rec = nullptr;
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        int visible = rec->visible;
        SpecRec* group_rec = rec->group;
        while (visible && group_rec) { /* visibility is a group issue... */
          if (!group_rec->visible)
            visible = false;
          else
            group_rec = group_rec->group;
        }
        if (rec->in_scene && !visible) {
          rec->in_scene = SceneObjectDel(G, rec->obj, true);
        } else if (visible && !rec->in_scene) {
          rec->in_scene = SceneObjectAdd(G, rec->obj);
        }
      }
    }
    I->ValidSceneMembers = true;
  }
}

void ExecutiveInvalidateGroups(PyMOLGlobals* G, bool force) {
  auto I = G->Executive;
  if (!(force || I->ValidGroups)) {
    return;
  }
  for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
    rec.group = nullptr;
    if (ExecutiveIsObjectType(rec, cObjectGroup)) {
      if (rec.group_member_list_id) {
        TrackerDelList(I->Tracker, rec.group_member_list_id);
      }
      rec.group_member_list_id = 0; /* not a list */
    }
  }
  I->ValidGroups = false;
  ExecutiveInvalidateSceneMembers(G);
  ExecutiveInvalidatePanelList(G);
  /* any changes to group structure means that we need to check scene
     members */
}

void ExecutiveUpdateGroups(PyMOLGlobals* G, bool force) {
  CExecutive* I = G->Executive;

  if (force || (!I->ValidGroups)) {
    CTracker* I_Tracker = I->Tracker;

    /* first, get rid of existing group lists */

    if (force || I->ValidGroups) ExecutiveInvalidateGroups(G, true);

    /* create empty lists for each group (also init grid_slot) */

    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      rec.group = nullptr;
      if (ExecutiveIsObjectType(rec, cObjectGroup)) {
        rec.group_member_list_id = TrackerNewList(I_Tracker, nullptr);
      }
    }

    /* iterate through and populate groups lists with their members */

    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      OVreturn_word result;
      if (OVreturn_IS_OK(
              (result = OVLexicon_BorrowFromCString(I->Lex, rec.group_name)))) {
        auto keyRes = I->Key.find(result.word);
        if (keyRes != I->Key.end()) {
          SpecRec* group_rec = nullptr;
          if (TrackerGetCandRef(I_Tracker, keyRes->second,
                                (TrackerRef**)(void*)&group_rec)) {
            int cycle = false;
            { /* don't close infinite loops */
              SpecRec* check_rec = group_rec;
              while (check_rec) {
                if (check_rec == &rec) {
                  cycle = true;
                  break;
                } else {
                  check_rec = check_rec->group;
                }
              }
            }
            if (!cycle) {
              rec.group = group_rec;
              TrackerLink(I_Tracker, rec.cand_id,
                          group_rec->group_member_list_id, 1);
            }
          }
        }
      }
    }

    /* note that it is possible to have infinite loops -- these must be
       allowed for later in the group expansion routine(s) */
    I->ValidGroups = true;
    ExecutiveInvalidatePanelList(G);
  }
}

static int ExecutiveGetObjectParentList(PyMOLGlobals* G, SpecRec* child) {
  int list_id = 0;
  ExecutiveUpdateGroups(G, false);
  {
    CExecutive* I = G->Executive;
    CTracker* I_Tracker = I->Tracker;
    int priority = 1; /* generations removed from child */
    int repeat_flag = true;
    SpecRec* group_rec = nullptr;

    list_id = TrackerNewList(I_Tracker, nullptr);
    while (child && child->group && repeat_flag) {
      OVreturn_word result;
      repeat_flag = false;
      if (OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(
                              I->Lex, child->group_name)))) {
        auto keyRes = I->Key.find(result.word);
        if (keyRes != I->Key.end()) {
          if (TrackerGetCandRef(I_Tracker, keyRes->second,
                                (TrackerRef**)(void*)&group_rec)) {
            if (TrackerLink(I_Tracker, keyRes->second, list_id, priority++)) {
              /* checking this prevents infinite loops */
              if (group_rec->group) {
                repeat_flag = true;
                child = group_rec;
              }
            }
          }
        }
      }
    }
  }
  return list_id;
}

int ExecutiveVdwFit(PyMOLGlobals* G, const char* s1, int state1, const char* s2,
                    int state2, float buffer, int quiet) {
  SelectorTmp tmpsele1(G, s1);
  SelectorTmp tmpsele2(G, s2);
  int sele1 = tmpsele1.getIndex();
  int sele2 = tmpsele2.getIndex();

  int ok = true;

  if ((sele1 >= 0) && (sele2 >= 0)) {
    ok = SelectorVdwFit(G, sele1, state1, sele2, state2, buffer, quiet);
  } else {
    ok = false;
  }
  return ok;
}

int ExecutiveAddKey(CExecutive* I, SpecRec* rec) {
  int ok = false;
  OVreturn_word result;
  if (OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, rec->name)))) {
    I->Key[result.word] = rec->cand_id;
    ok = true;
  }
  return ok;
}

static int ExecutiveDelKey(CExecutive* I, SpecRec* rec) {
  int ok = false;
  OVreturn_word result;
  if (OVreturn_IS_OK(
          (result = OVLexicon_BorrowFromCString(I->Lex, rec->name)))) {
    if (OVreturn_IS_OK(OVLexicon_DecRef(I->Lex, result.word))) {
      auto res = I->Key.find(result.word);
      if (res != I->Key.end()) {
        I->Key.erase(res);
        ok = true;
      }
    }
  }
  return ok;
}

static SpecRec* ExecutiveUnambiguousNameMatch(PyMOLGlobals* G,
                                              const char* name) {
  CExecutive* I = G->Executive;
  SpecRec* result = nullptr;
  SpecRec* rec = nullptr;
  int best = 0;
  int wm;
  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

  while (ListIterate(I->Spec, rec, next)) {
    wm = WordMatch(G, name, rec->name, ignore_case);
    if (wm < 0) { /* exact match, so this is valid */
      result = rec;
      best = wm;
      break;
    } else if ((wm > 0) && (best < wm)) {
      result = rec;
      best = wm;
    } else if ((wm > 0) && (best == wm)) { /* ambiguous match... no good */
      result = nullptr;
    }
  }
  return (result);
}

static SpecRec* ExecutiveAnyCaseNameMatch(PyMOLGlobals* G, const char* name) {
  CExecutive* I = G->Executive;
  SpecRec* result = nullptr;
  SpecRec* rec = nullptr;

  int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  while (ListIterate(I->Spec, rec, next)) {
    if (WordMatchExact(G, name, rec->name, ignore_case)) {
      result = rec;
      break;
    }
  }
  return (result);
}

void ExecutiveUpdateColorDepends(PyMOLGlobals* G, ObjectMolecule* mol) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;

  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      if (rec->obj->type == cObjectGadget) {
        ObjectGadget* gadget = (ObjectGadget*)rec->obj;
        if (gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp* ramp = (ObjectGadgetRamp*)gadget;
          if (ramp->RampType == cRampMol) {
            if (ramp->Mol == mol) {
              ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
              break;
            }
          }
        }
      }
    }
  }
}

void ExecutiveUpdateCoordDepends(
    PyMOLGlobals* G, ObjectMolecule* mol) { /* nasty, ugly, inefficient hack */

  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectGadget* gadget;
  int done_inv_all = false;
  int dynamic_measures = SettingGet_b(G, mol ? mol->Setting.get() : nullptr,
                                      nullptr, cSetting_dynamic_measures);

  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      switch (rec->obj->type) {
        case cObjectGadget:
          if (done_inv_all) break;
          gadget = (ObjectGadget*)rec->obj;
          if (gadget->GadgetType == cGadgetRamp) {
            ObjectGadgetRamp* ramp = (ObjectGadgetRamp*)gadget;
            if (ramp->RampType == cRampMol) {
              if (ramp->Mol == mol) {
                ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
                done_inv_all = true;
              }
            }
          }
          break;
        case cObjectMeasurement:
          if (dynamic_measures)
            ObjectDistMoveWithObject((ObjectDist*)rec->obj, mol);
          break;
        case cObjectAlignment:
          rec->obj->invalidate(cRepAll, cRepInvRep,
                               cSelectorUpdateTableAllStates);
          break;
      }
    }
  }
}

int ExecutiveValidNamePattern(PyMOLGlobals* G, const char* name) {
  int result = false;
  CWordMatcher* matcher;
  CWordMatchOptions options;
  const char* wildcard = SettingGetGlobal_s(G, cSetting_wildcard);

  WordMatchOptionsConfigNameList(&options, *wildcard,
                                 SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if (matcher) { /* this appears to be a pattern */
    result = true;
    WordMatcherFree(matcher);
  } else if (ExecutiveUnambiguousNameMatch(G, name)) {
    /* this does not appear to be a pattern, so it is an unambiguous partial
     * name? */
    result = true;
  }
  return result;
}

#define cExecNoExpand false
#define cExecExpandGroups true

static void ExecutiveExpandGroupsInList(PyMOLGlobals* G, int list_id,
                                        int expand_groups) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  int new_member_added = true;
  SpecRec* rec;
  ExecutiveUpdateGroups(G, false);
  while (new_member_added) { /* keep adding til we can't add no more */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    new_member_added = false;
    if (iter_id) {
      while ((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                  (TrackerRef**)(void*)&rec))) {
        if (rec && (rec->type == cExecObject) && rec->group_member_list_id &&
            (rec->obj->type == cObjectGroup)) {
          int group_iter_id =
              TrackerNewIter(I_Tracker, 0, rec->group_member_list_id);
          int group_cand_id;
          SpecRec* group_rec;
          if (group_iter_id) {
            while ((group_cand_id = TrackerIterNextCandInList(
                        I_Tracker, group_iter_id,
                        (TrackerRef**)(void*)&group_rec))) {
              if (group_rec && group_cand_id) {
                if (TrackerLink(I_Tracker, group_cand_id, list_id, 1))
                  new_member_added = true;
              }
            }
            TrackerDelIter(I_Tracker, group_iter_id);
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
  }
  /* now purge all groups from the expanded list */
  if (expand_groups != cExecExpandKeepGroups) {
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    while ((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                (TrackerRef**)(void*)&rec))) {
      if (rec && (rec->type == cExecObject) &&
          (rec->obj->type == cObjectGroup)) {
        TrackerUnlink(I_Tracker, cand_id, list_id);
      }
    }
  }
}

/* DON'T FORGET TO RELEASE LIST WHEN DONE!!! */
int ExecutiveGetNamesListFromPattern(PyMOLGlobals* G, const char* name,
                                     int allow_partial, int expand_groups) {
  CExecutive* I = G->Executive;
  int result = 0;
  CWordMatcher* matcher;
  CWordMatchOptions options;
  CTracker* I_Tracker = I->Tracker;
  const char* wildcard = SettingGetGlobal_s(G, cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  int group_found = false;
  SpecRec* rec;

  if (!name) return -1;

  // sanity check: name patterns are not object selections, bail if
  // parenthesis or operators in pattern
  if (strchr(name, '(') || strchr(name, ')') || strchr(name, '|')) {
    PRINTFB(G, FB_Executive, FB_Errors)
    " Names-Pattern-Error: Pattern looks like an atom selection"
    " (has parenthesis or operators), this is not supported for"
    " object name patterns.\n" ENDFB(G);
    return -1;
  }

  // special case: allow "not ..."
  bool match_not = false;
  if (WordMatchNoWild(G, "not ", name, false)) {
    match_not = true;
    name += 4;
  } else if (name[0] == '!') {
    match_not = true;
    name += 1;
  }

  // skip whitespace
  while (name[0] == ' ') {
    ++name;
  }

  bool match_enabled = WordMatchExact(G, "enabled", name, false);

  // ignore % and ? prefixes
  while (name[0] && (name[0] == '%' || name[0] == '?')) name++;

  WordMatchOptionsConfigNameList(&options, *wildcard,
                                 SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, /* force= */ match_not);
  if (matcher || match_enabled) {
    if (iter_id) {
      while ((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                  (TrackerRef**)(void*)&rec))) {
        if (rec && !(rec->type == cExecAll)) {
          bool test = match_enabled ? SpecRecIsEnabled(rec)
                                    : WordMatcherMatchAlpha(matcher, rec->name);
          if (test ^ match_not) {
            if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
              group_found = true;
            if (!result) result = TrackerNewList(I_Tracker, nullptr);
            if (result) {
              TrackerLink(I_Tracker, cand_id, result, 1);
            }
          }
        }
      }
    }
  } else if ((rec = ExecutiveFindSpec(G, name))) { /* only one name in list */
    if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
      group_found = true;
    result = TrackerNewList(I_Tracker, nullptr);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  } else if (allow_partial && (rec = ExecutiveUnambiguousNameMatch(G, name))) {
    if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup))
      group_found = true;
    result = TrackerNewList(I_Tracker, nullptr);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  }
  if (matcher) WordMatcherFree(matcher);
  if (iter_id) TrackerDelIter(I->Tracker, iter_id);
  if (group_found && expand_groups) {
    ExecutiveExpandGroupsInList(G, result, expand_groups);
  }
  return result;
}

static int ExecutiveUngroup(PyMOLGlobals* G, const char* members, int quiet) {
  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, members)) {
    rec.group_name[0] = 0;
    rec.group = nullptr;
  }
  ExecutiveInvalidateGroups(G, true);
  return true;
}

/**
 * Removes recs from a group rec
 *
 * @param[in] rec group rec to remove recs from
 * @param[out] discarded discarded recs
 *
 * Note: Caller is responsible for calling ExecutivePurgeSpec on discarded recs
 */

static void ExecutiveGroupPurge(
    PyMOLGlobals* G, SpecRec* rec,
    std::vector<DiscardedRec>* const discarded = nullptr) {
  std::vector<DiscardedRec> recs;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
  SpecRec* rec2 = nullptr;
  while (ListIterate(G->Executive->Spec, rec2, next)) {
    if ((rec2->group == rec) ||
        WordMatchExact(G, rec2->group_name, rec->name, ignore_case)) {
      bool save = discarded != nullptr;
      auto result = ExecutiveDelete(G, rec2->name, save);
      if (discarded && result) {
        discarded->insert(discarded->end(), result->begin(), result->end());
      }
      rec2 = nullptr;
    }
  }
}

int ExecutiveGroup(PyMOLGlobals* G, pymol::zstring_view nameView,
                   pymol::zstring_view membersView, int action, int quiet) {
  auto name = nameView.c_str();
  auto members = membersView.c_str();
  if (action == cExecutiveGroupUngroup) {
    // Up to PyMOL 2.3 the member argument was ignored and 'name' used for
    // ungrouping
    if (name[0]) {
      members = name;
    }
    return ExecutiveUngroup(G, members, quiet);
  }

  int ok = true;
  CExecutive* I = G->Executive;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);

  ObjectNameType valid_name;
  UtilNCopy(valid_name, name, sizeof(ObjectNameType));
  ObjectMakeValidName(G, valid_name);

  pymol::CObject* obj = ExecutiveFindObjectByName(G, valid_name);

  if (obj && (obj->type != cObjectGroup)) {
    PRINTFB(G, FB_Executive, FB_Errors)
    " Group-Error: object '%s' is not a group object.", name ENDFB(G);
    ok = false;
  } else {
    if ((!obj) && (action == cExecutiveGroupAdd)) {
      obj = new ObjectGroup(G);
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, false, true);
      }
    }
  }
  if ((!members[0]) &&
      ((action == cExecutiveGroupOpen) || (action == cExecutiveGroupClose) ||
       (action == cExecutiveGroupToggle) || (action == cExecutiveGroupEmpty) ||
       (action == cExecutiveGroupPurge) || (action == cExecutiveGroupExcise) ||
       (action == cExecutiveGroupRaise))) {
    ExecutiveUpdateGroups(G, false);
    {
      CTracker* I_Tracker = I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G, name, true, false);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec* rec;

      while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                       (TrackerRef**)(void*)&rec)) {
        if (rec) {
          ObjectGroup* objGroup = nullptr;
          if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
            objGroup = (ObjectGroup*)rec->obj;
          }

          switch (action) {
            case cExecutiveGroupOpen:
              if (objGroup) objGroup->OpenOrClosed = 1;
              break;
            case cExecutiveGroupClose:
              if (objGroup) objGroup->OpenOrClosed = 0;
              break;
            case cExecutiveGroupToggle:
              if (objGroup) objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
              break;
            case cExecutiveGroupEmpty:
              if (objGroup) {
                SpecRec* rec2 = nullptr;
                while (ListIterate(I->Spec, rec2, next)) {
                  if ((rec2->group == rec) ||
                      WordMatchExact(G, rec2->group_name, rec->name,
                                     ignore_case)) {
                    rec2->group = nullptr;
                    rec2->group_name[0] = 0;
                  }
                }
              }
              break;
            case cExecutiveGroupPurge:
              if (objGroup) {
                ExecutiveGroupPurge(G, rec);
              }
              break;
            case cExecutiveGroupExcise:
              if (objGroup) {
                if (rec->group_name[0]) {
                  /* cascade group members up to the surrounding group */
                  SpecRec* rec2 = nullptr;
                  while (ListIterate(I->Spec, rec2, next)) {
                    if ((rec2->group == rec) ||
                        WordMatch(G, rec->name, rec2->group_name,
                                  ignore_case)) {
                      strcpy(rec2->group_name, rec->group_name);
                      rec2->group = rec->group;
                    }
                  }
                } else if ((rec->type == cExecObject) &&
                           (rec->obj->type == cObjectGroup)) {
                  /* and/or delete their group membership */
                  SpecRec* rec2 = nullptr;
                  while (ListIterate(I->Spec, rec2, next)) {
                    if ((rec2->group == rec) ||
                        WordMatch(G, rec->name, rec2->group_name,
                                  ignore_case)) {
                      rec2->group_name[0] = 0;
                      rec2->group = nullptr;
                    }
                  }
                }
                ExecutiveDelete(G, rec->name);
              }
              break;
            case cExecutiveGroupRaise:
              if (objGroup) {
                rec->group = nullptr;
                rec->group_name[0] = 0;
              }
              break;
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
      ExecutiveInvalidateGroups(G, true);
    }
  } else {
    if (obj && (obj->type == cObjectGroup)) {
      ObjectGroup* objGroup = (ObjectGroup*)obj;
      switch (action) {
        case cExecutiveGroupOpen:
          objGroup->OpenOrClosed = 1;
          break;
        case cExecutiveGroupClose:
          objGroup->OpenOrClosed = 0;
          break;
        case cExecutiveGroupToggle:
          objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
          break;
      }
      if (members[0] && (action != cExecutiveGroupRemove))
        action = cExecutiveGroupAdd;

      switch (action) {
        case cExecutiveGroupAdd: {
          for (auto& rec :
               ExecutiveGetSpecRecsFromPattern(G, members, true, false)) {
            if (rec.type != cExecObject ||
                (rec.type == cExecObject && rec.obj != obj)) {
              UtilNCopy(rec.group_name, valid_name, sizeof(WordType));
              if (!quiet) {
                PRINTFB(G, FB_Executive, FB_Actions)
                " Executive: adding '%s' to group '%s'.\n", rec.name,
                    rec.group_name ENDFB(G);
              }
            }
          }
        } break;
      }

      ExecutiveInvalidateGroups(G, true);
    }
  }
  return ok;
}

static void ExecutiveGetUniqueIDAtomVLADict(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  I->m_eoo.clear();
  I->m_eoo.reserve(1000);

  std::size_t n_oi = 0;
  auto recs = pymol::make_list_adapter(G->Executive->Spec);
  for (auto& rec : recs) {
    if (rec.type != cExecObject) {
      continue;
    }
    if (rec.obj->type != cObjectMolecule) {
      continue;
    }
    auto obj = static_cast<ObjectMolecule*>(rec.obj);
    auto n_atom = obj->NAtom;
    const AtomInfoType* ai = obj->AtomInfo.data();
    for (int a = 0; a < n_atom; a++) {
      auto id = ai->unique_id;
      if (id != 0) {
        auto o2oIt = I->m_id2eoo.find(id);
        if (o2oIt == I->m_id2eoo.end()) {
          I->m_id2eoo[id] = n_oi++;
          I->m_eoo.push_back(ExecutiveObjectOffset{obj, a});
        }
      }
      ai++;
    }
  }
}

int ExecutiveDrawCmd(PyMOLGlobals* G, int width, int height, int antialias,
                     int entire_window, int quiet) {
  CExecutive* I = G->Executive;
  Extent2D extent{static_cast<std::uint32_t>(width),
                  static_cast<std::uint32_t>(height)};
  if (width == 0 && height == 0) {
    extent = SceneGetExtent(G);
  }
  if (antialias < 0) antialias = SettingGetGlobal_i(G, cSetting_antialias);
  if (entire_window) {
    SceneInvalidateCopy(G, false);
    OrthoDirty(G);
    I->CaptureFlag = true;
  } else {
    if (SettingGetGlobal_i(G, cSetting_draw_mode) == -1) {
      ExecutiveSetSettingFromString(G, cSetting_draw_mode, "-2", "", -1, true,
                                    true);
      SceneUpdate(G, false);
    }
    SceneDeferImage(G, extent, nullptr, antialias, -1.0, cMyPNG_FormatPNG,
                    quiet, nullptr);
  }
  return 1;
}

int ExecutiveCountNames(PyMOLGlobals* G) {
  int count = 0;
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) count++;

  return (count);
}

static int ReorderOrderFn(PyMOLGlobals* G, const SpecRec* const* rec, int l,
                          int r) {
  return (WordCompare(G, rec[l]->name, rec[r]->name, true) <= 0);
}

static int SpecRecListPopulate(SpecRec** list, SpecRec* first,
                               const char* group_name) {
  /* Add items to list such that group items are ordered behind their group
   * records. author: Thomas Holder, 2013-08 runtime: O(N * G) where N is the
   * number of SpecRecs and G the number of groups args: list: pointer to the
   * end of the list first: pointer to G->Executive->Spec (root of linked list)
   *       group_name: if empty string, append ungrouped items, otherwise append
   * group items returns: number of appended items
   */
  SpecRec* rec;
  int a = 0;
  for (rec = first; rec; rec = rec->next) {
    if (!strcmp(group_name, rec->group_name)) {
      list[a++] = rec;
      if (rec->type == cExecObject && rec->obj->type == cObjectGroup)
        a += SpecRecListPopulate(list + a, first, rec->name);
    }
  }
  return a;
}

pymol::Result<> ExecutiveOrder(PyMOLGlobals* G, pymol::zstring_view s1View,
                               int sort, int location) {
  CExecutive* I = G->Executive;
  CTracker* I_Tracker = I->Tracker;
  auto s1 = s1View.c_str();
  CWordList* word_list = WordListNew(G, s1);
  int n_names = ExecutiveCountNames(G);

  if (n_names) {
    SpecRec **list, **subset, **sorted;
    int* index = nullptr;
    int n_sel;
    int source_row = -1;
    int min_row = -1;
    list = pymol::malloc<SpecRec*>(n_names);
    subset = pymol::calloc<SpecRec*>(n_names);
    sorted = pymol::calloc<SpecRec*>(n_names);
    index = pymol::malloc<int>(n_names);
    if (list && subset) {
      /* create an array of current names */
      {
        int a = 0;
        /* copy all names into array */
        /* update Thomas Holder 2013-08: order group members behind group.
         * fixes PYMOL-1382 (eventually).
         * I don't know if this is the best place for the fix, but it's my best
         * guess.
         */
#ifndef NDEBUG
        auto const list_size =
#endif
            SpecRecListPopulate(list, I->Spec, "");
        assert(list_size == n_names);
        /* unlink them */
        for (a = 0; a < n_names; a++) {
          list[a]->next = nullptr;
        }
      }
      /* transfer matching names to the subset array */
      {
        int a;
        int entry;
        int min_entry = word_list->n_word;
        const char* word = nullptr;
        int word_iter = 0;
        while (WordListIterate(G, word_list, &word, &word_iter)) {
          int list_id = ExecutiveGetNamesListFromPattern(G, word, true, false);
          SpecRec* rec = nullptr;
          entry = word_iter - 1;
          for (a = n_names - 1; a > 0; a--) { /* skipping zeroth */
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                             (TrackerRef**)(void*)&rec)) {
              if (rec == list[a]) {
                if ((a < min_row) || (min_row < 0)) min_row = a;
                if (entry <= min_entry) {
                  source_row = a; /* where will new list be inserted... */
                  min_entry = entry;
                }
                /* ensure that each record appears only once */
                rec->next = subset[entry];
                subset[entry] = rec;
                list[a] = nullptr;
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
          }
          TrackerDelList(I_Tracker, list_id);
        }
        if (word_list->n_word &&
            WordMatchExact(G, word_list->start[0], cKeywordAll, true))
          location = -1; /* set to top if "all" is first in list */
      }
      /* expand the selected entries */
      {
        SpecRec *rec, *last;
        int b;
        n_sel = 0;
        for (b = 0; b < word_list->n_word; b++) {
          rec = subset[b];
          while (rec) {
            sorted[n_sel++] = rec;
            last = rec;
            rec = rec->next;
            last->next = nullptr;
          }
        }
      }
      /* sort the selected entries, if requested */
      if (sort) {
        UtilCopyMem(subset, sorted, sizeof(SpecRec*) * n_sel);
        {
          int a;
          UtilSortIndexGlobals(G, n_sel, subset, index,
                               (UtilOrderFnGlobals*)ReorderOrderFn);
          for (a = 0; a < n_sel; a++) {
            sorted[a] = subset[index[a]];
          }
        }
      }
      /* reassemble the list using the new order */
      {
        SpecRec* spec = nullptr;
        SpecRec* last = nullptr;
        int a, b;
        int flag;
        for (a = 0; a < n_names; a++) {
          flag = false;
          if (sorted) { /* not yet added */
            switch (location) {
              case -1: /* top */
                if (a == 1) flag = true;
                break;
              case -2: /* upper */
                if (min_row >= 0) {
                  if (a == min_row) flag = true;
                } else if (!list[a])
                  flag = true;
                break;
              case 0: /* current */
                if (source_row >= 0) {
                  if (a == source_row) flag = true;
                } else if (!list[a])
                  flag = true;
                break;
            }
          }
          if (flag) {
            for (b = 0; b < n_sel; b++) {
              if (sorted[b]) {
                if (last) last->next = sorted[b];
                last = sorted[b];
                if (!spec) spec = last;
              }
            }
            FreeP(sorted);
          }
          if (list[a]) {
            if (last) last->next = list[a];
            last = list[a];
            if (!spec) spec = last;
          }
        }
        if (sorted) { /* still not yet readded? */
          for (b = 0; b < n_sel; b++) {
            if (sorted[b]) {
              if (last) last->next = sorted[b];
              last = sorted[b];
              if (!spec) spec = last;
            }
          }
        }
        I->Spec = spec;
        OrthoDirty(G);
        SeqChanged(G);
      }
      FreeP(index);
      FreeP(sorted);
      FreeP(list);
      FreeP(subset);
    }
    ExecutiveInvalidatePanelList(G);
  }
  WordListFreeP(word_list);
  return {};
}

pymol::vla<ObjectMolecule*> ExecutiveGetObjectMoleculeVLA(PyMOLGlobals* G,
                                                          const char* sele) {
  ObjectMolecule** result = nullptr;
  int s1 = SelectorIndexByName(G, sele);
  if (s1 >= 0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetObjects;
    op.obj1VLA = (ObjectMolecule**)VLAlloc(ObjectMolecule*, 10);
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, s1, &op);
    result = (ObjectMolecule**)op.obj1VLA;
    VLASize(result, ObjectMolecule*, op.i1);
  }
  return pymol::vla_take_ownership(result);
}

/**
 * Return the selector index of the "active" alignment.
 */
int ExecutiveGetActiveAlignmentSele(PyMOLGlobals* G) {
  const char* alignment = ExecutiveGetActiveAlignment(G);
  if (alignment && alignment[0]) {
    return SelectorIndexByName(G, alignment);
  }
  return -1;
}

/**
 * Return the name of the "active" alignment. That is:
 *
 * 1) The "seq_view_alignment" setting if it's set
 * 2) or the name of the first enabled alignment object
 * 3) or nullptr
 */
const char* ExecutiveGetActiveAlignment(PyMOLGlobals* G) {
  const char* alignment = SettingGetGlobal_s(G, cSetting_seq_view_alignment);
  if (alignment && alignment[0]) { /* explicit alignment setting name */
    return alignment;
  } else { /* otherwise, use the first active alignment */
    SpecRec* rec = nullptr;
    CExecutive* I = G->Executive;
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->visible) {
        if (rec->type == cExecObject)
          if (rec->obj->type == cObjectAlignment) {
            return rec->obj->Name;
          }
      }
    }
  }
  return nullptr;
}

int ExecutiveGetActiveSele(PyMOLGlobals* G) {
  ObjectNameType name;
  if (ExecutiveGetActiveSeleName(G, name, false, false))
    return SelectorIndexByName(G, name);
  else
    return -1;
}

int ExecutiveGetActiveSeleName(PyMOLGlobals* G, char* name, int create_new,
                               int log) {
  /* TODO: cache/optimize to avoid table scan */

  int result = false;
  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecSelection) {
      if (rec->visible) {
        strcpy(name, rec->name);
        result = true;
      }
    }
  }
  if ((!result) && create_new) {
    if (SettingGetGlobal_b(G, cSetting_auto_number_selections)) {
      int sel_num = SettingGetGlobal_i(G, cSetting_sel_counter) + 1;

      SettingSetGlobal_i(G, cSetting_sel_counter, sel_num);
      sprintf(name, "sel%02d", sel_num);
      SelectorCreateEmpty(G, name, -1);
      if (log) {
        if (SettingGetGlobal_i(G, cSetting_logging)) {
          OrthoLineType buf2;
          sprintf(buf2, "cmd.select('%s','none')\n", name);
          PLog(G, buf2, cPLog_no_flush);
        }
      }
    } else {
      sprintf(name, "sele");
      SelectorCreateEmpty(G, name, -1);
      if (log) {
        OrthoLineType buf2;
        sprintf(buf2, "cmd.select('%s','none')\n", name);
        PLog(G, buf2, cPLog_no_flush);
      }
    }
  }
  return result;
}

int ExecutiveGetActiveSeleName(PyMOLGlobals* G, std::string& name,
                               int create_new, int log) {
  ObjectNameType name_arg;
  auto result = ExecutiveGetActiveSeleName(G, name_arg, create_new, log);
  name = name_arg;
  return result;
}

pymol::Result<> ExecutiveSetObjectColor(PyMOLGlobals* G, const char* name,
                                        const char* color, int quiet) {
  int col_ind = ColorGetIndex(G, color);
  auto obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    obj->Color = col_ind;
  } else {
    return pymol::make_error("Object ", name, " not found.");
  }
  return {};
}

int ExecutiveGetObjectColorIndex(PyMOLGlobals* G, const char* name) {
  int result = -1;
  pymol::CObject* obj = nullptr;
  obj = ExecutiveFindObjectByName(G, name);
  if (obj) {
    result = obj->Color;
  }
  return (result);
}

pymol::Result<std::array<float, 3>> ExecutiveGetAtomVertex(PyMOLGlobals* G,
                                                           const char* s1,
                                                           int state,
                                                           int index) {
  auto tmpsele1 = SelectorTmp::make(G, s1);
  p_return_if_error(tmpsele1);

  switch (tmpsele1->getAtomCount()) {
    case 0:
      return pymol::make_error("Empty Selection");
    case 1:
      return SelectorGetSingleAtomVertex(G, tmpsele1->getIndex(), state);
  }

  assert(tmpsele1->getAtomCount() > 0);
  return pymol::make_error("More than one atom found");
}

void ExecutiveMakeUnusedName(PyMOLGlobals* G, char* prefix, int length,
                             bool alwaysnumber, int start,
                             const char* pattern) {
  if (!prefix[0]) strcpy(prefix, "obj");

  int prefixlen = strlen(prefix);
  int suffixlen = length - prefixlen;
  char* end = prefix + prefixlen;

  for (int cnt = start; alwaysnumber || ExecutiveValidName(G, prefix); ++cnt) {
    snprintf(end, suffixlen, pattern, cnt);
    alwaysnumber = false;
  }
}

std::string ExecutiveGetUnusedName(PyMOLGlobals* G, const char* prefix,
                                   bool alwaysnumber) {
  OrthoLineType unused_name;
  strcpy(unused_name, prefix);

  ObjectMakeValidName(G, unused_name);

  ExecutiveMakeUnusedName(G, unused_name, OrthoLineLength, alwaysnumber);

  return std::string(unused_name);
}

int ExecutiveProcessObjectName(PyMOLGlobals* G, const char* proposed,
                               char* actual) {
  int result = true;
  UtilNCopy(actual, proposed, sizeof(ObjectNameType));
  if (SettingGetGlobal_b(G, cSetting_validate_object_names))
    ObjectMakeValidName(G, actual);
  if (SettingGetGlobal_b(G, cSetting_auto_rename_duplicate_objects) ||
      !proposed[0]) {
    ExecutiveMakeUnusedName(G, actual, sizeof(ObjectNameType), false, 2, "_%d");
  }
  return result;
}

pymol::Result<std::vector<DiscardedRec>> ExecutiveSetName(
    PyMOLGlobals* G, pymol::zstring_view old_name_view,
    pymol::zstring_view new_name_view, bool save) {
  std::vector<DiscardedRec> discarded;
  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;
  int found = false;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
  auto old_name = old_name_view.c_str();
  auto new_name = new_name_view.c_str();

  ObjectNameType name;
  UtilNCopy(name, new_name, sizeof(ObjectNameType));

  if (ObjectMakeValidName(name)) {
    PRINTFB(G, FB_Executive, FB_Warnings)
    " Warning: Invalid characters in '%s' have been replaced or stripped\n",
        name ENDFB(G);
  }

  if (!name[0]) {
    return pymol::make_error("Blank names not allowed.");
  } else if (WordMatchExact(G, name, cKeywordSame, ignore_case) ||
             SelectorNameIsKeyword(G, name)) {
    return pymol::make_error("Name ", name, " is a selection keyword.");
  }
  if (!WordMatchExact(G, name, old_name, ignore_case)) {
    while (ListIterate(I->Spec, rec, next)) {
      if (found) break;
      switch (rec->type) {
        case cExecObject:
          if (WordMatchExact(G, rec->obj->Name, old_name, ignore_case)) {
            ExecutiveDelKey(I, rec);
            auto r = ExecutiveDelete(G, name, save);
            discarded = std::move(*r);
            assert(r);
            ObjectSetName(rec->obj, name);
            UtilNCopy(rec->name, rec->obj->Name, WordLength);
            ExecutiveAddKey(I, rec);
            if (rec->obj->type == cObjectMolecule) {
              /*
                 SelectorDelete(G,old_name);
                 ExecutiveUpdateObjectSelection(G,rec->obj);
               */
              SelectorSetName(G, name, old_name);
              SceneChanged(G);
              SeqChanged(G);
            }
            if (rec->obj->type == cObjectMap)
              ExecutiveInvalidateMapDependents(G, old_name, name);

            found = true;
          }
          break;
        case cExecSelection:
          if (WordMatchExact(G, rec->name, old_name, ignore_case)) {
            if (SelectorSetName(G, name, old_name)) {
              auto r = ExecutiveDelete(G, name, save); /* just in case */
              discarded = std::move(*r);
              assert(r);
              ExecutiveDelKey(I, rec);
              UtilNCopy(rec->name, name, WordLength);
              ExecutiveAddKey(I, rec);
              found = true;
              OrthoDirty(G);
            }
          }
          break;
      }
    }
    if (!found)
      return pymol::make_error("Could not find object named ", name);
    else {
      rec = nullptr;
      int old_name_len = strlen(old_name);
      int new_name_len = strlen(name);
      ObjectNameType childname;
      UtilNCopy(childname, name, sizeof(ObjectNameType));
      while (ListIterate(I->Spec, rec, next)) {
        if (WordMatchExact(G, rec->group_name, old_name, ignore_case)) {
          UtilNCopy(rec->group_name, name, WordLength);
          // rename group members for group_auto_mode
          if (strncmp(rec->name, old_name, old_name_len) == 0 &&
              rec->name[old_name_len] == '.') {
            UtilNCopy(childname + new_name_len, rec->name + old_name_len,
                      sizeof(ObjectNameType) - new_name_len);
            ExecutiveSetName(G, rec->name, childname);
          }
        }
      }
      ExecutiveInvalidateGroups(G, false);
    }
  }
  ColorRenameExt(G, old_name_view, new_name_view);
  return discarded;
}

static int fStrOrderFn(const char* const* array, int l, int r) {
  return strcmp(array[l], array[r]) < 0;
}

/**
 * Returns an VLA with pointers into G->Lexicon
 */
pymol::Result<std::vector<const char*>> ExecutiveGetChains(PyMOLGlobals* G,
                                                           const char* s1,
                                                           int state) {
  std::set<lexidx_t> chains;
  int c = 0;
  ObjectMoleculeOpRec op;

  SETUP_SELE_DEFAULT(1);

  {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetChains;
    op.ii1 = (int*)(void*)&chains;  // pointer pack
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op);
    std::vector<const char*> result(chains.size());
    for (const auto& chain : chains) {
      result[c++] = LexStr(G, chain);
    }
    // sort the array
    UtilSortInPlace(G, result.data(), chains.size(), sizeof(char*),
                    (UtilOrderFn*)fStrOrderFn);
    return result;
  }
}

int ExecutiveValidateObjectPtr(PyMOLGlobals* G, pymol::CObject* ptr,
                               int object_type) {
  /* this routine needs to be sped up significantly... */

  CExecutive* I = G->Executive;
  int ok = false;
  SpecRec* rec = nullptr;

  while (ListIterate(I->Spec, rec, next)) {
    if (rec->obj == ptr) {
      if (rec->type == cExecObject) {
        if ((!object_type) || (rec->obj->type == object_type)) {
          ok = true;
          break;
        }
      }
    }
  }
  return (ok);
}
int*** ExecutiveGetBondPrint(PyMOLGlobals* G, const char* name, int max_bond,
                             int max_type, int* dim) {
  int*** result = nullptr;
  pymol::CObject* obj;
  ObjectMolecule* objMol;

  obj = ExecutiveFindObjectByName(G, name);
  if (obj->type == cObjectMolecule) {
    objMol = (ObjectMolecule*)obj;
    result = ObjectMoleculeGetBondPrint(objMol, max_bond, max_type, dim);
  }
  return (result);
}

int ExecutiveValidName(PyMOLGlobals* G, const char* name) {
  if (!ExecutiveFindSpec(G, name)) {
    int ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);

    if (!WordMatchExact(G, name, cKeywordAll, ignore_case))
      if (!WordMatchExact(G, name, cKeywordSame, ignore_case))
        if (!WordMatchExact(G, name, cKeywordCenter, ignore_case))
          if (!WordMatchExact(G, name, cKeywordOrigin, ignore_case))
            return false;
  }
  return true;
}
float* ExecutiveGetVertexVLA(PyMOLGlobals* G, const char* s1, int state) {
  /* returns nullptr if none found */

  float* result = nullptr;
  ObjectMoleculeOpRec op1;
  int sele1;
  sele1 = SelectorIndexByName(G, s1);
  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.nvv1 = 0;
    op1.vv1 = VLAlloc(float, 1000);
    if (state >= 0) {
      op1.cs1 = state;
      op1.code = OMOP_SingleStateVertices;
    } else {
      op1.code = OMOP_VERT;
    }
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    VLASize(op1.vv1, float, op1.nvv1 * 3);
    result = op1.vv1;
  }
  return (result);
}
void ExecutiveSetLastObjectEdited(PyMOLGlobals* G, pymol::CObject* o) {
  CExecutive* I = G->Executive;
  I->LastEdited = o;
}
pymol::CObject* ExecutiveGetLastObjectEdited(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  return (I->LastEdited);
}
int ExecutiveSaveUndo(PyMOLGlobals* G, const char* s1, int state) {
  int sele1;
  ObjectMoleculeOpRec op1;

  if (state < 0) state = SceneGetState(G);
  sele1 = SelectorIndexByName(G, s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2 = 0;
  if (sele1 >= 0) {
    op1.code = OMOP_SaveUndo;
    op1.i1 = state;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
  }
  return (op1.i2);
}
pymol::Result<> ExecutiveSetTitle(PyMOLGlobals* G, const char* name, int state,
                                  const char* text) {
  ObjectMolecule* obj;
  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if (!obj) {
    return pymol::make_error("Object ", name, " not found.");
  } else {
    auto res = ObjectMoleculeSetStateTitle(obj, state, text);
    if (!res) {
      return res;
    }
  }
  SceneDirty(G);
  return {};
}
const char* ExecutiveGetTitle(PyMOLGlobals* G, const char* name, int state) {
  ObjectMolecule* obj;
  obj = ExecutiveFindObjectMoleculeByName(G, name);
  if (!obj) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    "Error: object %s not found.\n", name ENDFB(G);
    return nullptr;
  }
  return ObjectMoleculeGetStateTitle(obj, state);
}
/**
 * Implementation of `cmd.get_names()`
 *
 * @param mode see modules/pymol/querying.py
 * @param enabled_only only include enabled names
 * @param s0 atom selection expression, limits results to molecular objects and
 * named selections which have at least one atom in the selection.
 */
pymol::Result<std::vector<const char*>> ExecutiveGetNames(PyMOLGlobals* G,
                                                          int mode,
                                                          int enabled_only,
                                                          const char* s0) {
  enum {
    cGetNames_all = 0,
    cGetNames_objects = 1,
    cGetNames_selections = 2,
    cGetNames_public = 3,
    cGetNames_public_objects = 4,
    cGetNames_public_selections = 5,
    cGetNames_public_nongroup_objects = 6,
    cGetNames_public_group_objects = 7,
    cGetNames_nongroup_objects = 8,
    cGetNames_group_objects = 9,
  };

  bool include_objects =  //
      mode == cGetNames_all || mode == cGetNames_objects ||
      mode == cGetNames_public || mode == cGetNames_public_objects;
  bool include_selections =
      mode == cGetNames_all || mode == cGetNames_selections ||
      mode == cGetNames_public || mode == cGetNames_public_selections;
  bool include_group_objects =
      mode == cGetNames_group_objects || mode == cGetNames_public_group_objects;
  bool include_nongroup_objects =  //
      mode == cGetNames_nongroup_objects ||
      mode == cGetNames_public_nongroup_objects;
  bool public_only =
      mode >= cGetNames_public && mode <= cGetNames_public_group_objects;

  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  std::vector<const char*> result;
  pymol::Result<SelectorTmp> tmpsele0;
  SelectorID_t sele0 = cSelectionInvalid;
  int incl_flag;
  if (s0[0]) {
    tmpsele0 = SelectorTmp::make(G, s0);
    p_return_if_error(tmpsele0);
    sele0 = tmpsele0->getIndex();
    assert(sele0 != cSelectionInvalid);
  }

  while (ListIterate(I->Spec, rec, next)) {
    incl_flag = 0;
    if ((rec->type == cExecObject &&
         (include_objects ||
          (rec->obj->type != cObjectGroup && include_nongroup_objects) ||
          (rec->obj->type == cObjectGroup && include_group_objects))) ||
        (rec->type == cExecSelection && include_selections)) {
      if (!public_only || rec->name[0] != '_') {
        if ((!enabled_only) || (rec->visible)) {
          incl_flag = 0;
          if (sele0 < 0)
            incl_flag = 1;
          else
            switch (rec->type) {
              case cExecObject:
                if (rec->obj->type == cObjectMolecule) {
                  int a;
                  ObjectMolecule* obj_mol = (ObjectMolecule*)rec->obj;
                  const AtomInfoType* ai = obj_mol->AtomInfo.data();
                  for (a = 0; a < obj_mol->NAtom; a++) {
                    if (SelectorIsMember(G, ai->selEntry, sele0)) {
                      incl_flag = 1;
                      break;
                    }
                    ai++;
                  }
                }
                break;
              case cExecSelection:
                if (SelectorCheckIntersection(
                        G, sele0, SelectorIndexByName(G, rec->name))) {
                  incl_flag = 1;
                  break;
                }
                break;
            }
          if (incl_flag) {
            result.push_back(rec->name);
          }
        }
      }
    }
  }
  return (result);
}
/**
 * Return true if `name` is the name of an object molecule or a named selection
 */
bool ExecutiveIsMoleculeOrSelection(PyMOLGlobals* G, const char* name) {
  if (!strcmp(name, cKeywordAll)) return true;

  SpecRec* rec = ExecutiveFindSpec(G, name);
  if (rec && ((rec->type == cExecObject && rec->obj->type == cObjectMolecule) ||
              (rec->type == cExecSelection)))
    return true;

  return false;
}
/**
 * @param name Name of the object or selection
 * @return Type label like "object:molecule" or "selection"
 */
pymol::Result<char const*> ExecutiveGetType(PyMOLGlobals* G, const char* name) {
  auto rec = ExecutiveFindSpec(G, name);

  if (!rec) {
    return pymol::make_error("object not found");
  }

  if (rec->type == cExecSelection) {
    return "selection";
  }

  if (rec->type != cExecObject) {
    return "";
  }

  switch (rec->obj->type) {
    case cObjectMolecule:
      return "object:molecule";
    case cObjectMap:
      return "object:map";
    case cObjectMesh:
      return "object:mesh";
    case cObjectSlice:
      return "object:slice";
    case cObjectSurface:
      return "object:surface";
    case cObjectMeasurement:
      return "object:measurement";
    case cObjectCGO:
      return "object:cgo";
    case cObjectGroup:
      return "object:group";
    case cObjectVolume:
      return "object:volume";
    case cObjectAlignment:
      return "object:alignment";
    case cObjectGadget:
      return "object:ramp";
    default:
      return "object:";
  }
}
pymol::Result<> ExecutiveUpdateCmd(PyMOLGlobals* G, const char* s0,
                                   const char* s1, int sta0, int sta1,
                                   int method, int quiet) {
  SETUP_SELE_DEFAULT_PREFIXED(0, cSelectionInvalid);
  SETUP_SELE_DEFAULT_PREFIXED(1, cSelectionInvalid);

  return SelectorUpdateCmd(G, sele0, sele1, sta0, sta1, method, quiet);
}
void ExecutiveRebuildAll(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  PRINTFD(G, FB_Executive)
  " ExecutiveRebuildAll: entered.\n" ENDFD;
  auto defer_builds_mode = SettingGet<bool>(G, cSetting_defer_builds_mode);
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      auto level = cRepInvAll;
      switch (rec->obj->type) {
        case cObjectMeasurement:
          ObjectDistInvalidateRep((ObjectDist*)rec->obj, cRepAll);
          break;
        case cObjectMolecule:
          level = defer_builds_mode ? cRepInvPurge : cRepInvRep;
        case cObjectSurface:
        case cObjectMesh:
        case cObjectSlice:
        case cObjectAlignment:
        case cObjectCGO:
          rec->obj->invalidate(cRepAll, level, -1);
          break;
      }
    }
  }
  SeqChanged(G);
  SceneChanged(G);
}
void ExecutiveRebuildAllObjectDist(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      if (rec->obj->type == cObjectMeasurement) {
        ObjectDistInvalidateRep((ObjectDist*)rec->obj, cRepAll);
      }
    }
  }
  SceneInvalidate(G);
}
void ExecutiveUndo(PyMOLGlobals* G, int dir) {
  CExecutive* I = G->Executive;
  ObjectMolecule *obj = nullptr, *compObj;
  SpecRec* rec = nullptr;

  auto o = ExecutiveGetLastObjectEdited(G);
  PRINTFB(G, FB_Executive, FB_Debugging)
  " ExecutiveUndo: last object %p\n", (void*)o ENDFB(G);
  if (o)
    if (o->type == cObjectMolecule) obj = (ObjectMolecule*)o;
  /* make sure this is still a real object */
  if (obj) {
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject)
        if (rec->obj->type == cObjectMolecule) {
          compObj = (ObjectMolecule*)rec->obj;
          if (obj == compObj) {
            ObjectMoleculeUndo(obj, dir);
            break;
          }
        }
    }
  }
}
float ExecutiveOverlap(PyMOLGlobals* G, const char* s1, int state1,
                       const char* s2, int state2, float adjust) {
  SelectorTmp tmpsele1(G, s1);
  SelectorTmp tmpsele2(G, s2);
  int sele1 = tmpsele1.getIndex();
  int sele2 = tmpsele2.getIndex();

  float result = 0.0;

  if (state1 < 0) state1 = 0;
  if (state2 < 0) state2 = 0;

  if ((sele1 >= 0) && (sele2 >= 0))
    result = SelectorSumVDWOverlap(G, sele1, state1, sele2, state2, adjust);

  return (result);
}
char* ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals* G, const char* name,
                                    int state, int format, int quiet) {
  char* result = nullptr;
  if ((!name) || (!name[0]) || (strcmp(name, "(all)") == 0)) {
    /* use current alignment as the default */
    name = SettingGetGlobal_s(G, cSetting_seq_view_alignment);
    if (name[0] == 0) {
      SpecRec* rec = nullptr;
      CExecutive* I = G->Executive;
      while (ListIterate(I->Spec, rec, next)) {
        if (rec->visible) {
          if (rec->type == cExecObject)
            if (rec->obj->type == cObjectAlignment) {
              name = rec->obj->Name;
              break;
            }
        }
      }
    }
  }
  if (!name) {
    ErrMessage(G, " Executive", "invalid alignment object name.");
  } else {
    pymol::CObject* obj = ExecutiveFindObjectByName(G, name);

    if (!obj) {
      ErrMessage(G, " Executive", "alignment object not found.");
    } else if (obj->type != cObjectAlignment) {
      ErrMessage(G, " Executive", "invalid object type.");
    } else {
      ObjectAlignmentAsStrVLA(G, (ObjectAlignment*)obj, state, format, &result);
    }
  }
  return (result);
}
pymol::Result<> ExecutiveSeleToObject(PyMOLGlobals* G, const char* name,
                                      const char* s1, int source, int target,
                                      int discrete, int zoom, int quiet,
                                      int singletons, int copy_properties) {
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  int ok = false;
  ObjectNameType valid_name;

  UtilNCopy(valid_name, name, sizeof(valid_name));
  if (SettingGetGlobal_b(G, cSetting_validate_object_names)) {
    ObjectMakeValidName(G, valid_name);
    name = valid_name;
  }
  {
    int exists = (ExecutiveFindObjectMoleculeByName(G, name) != nullptr);

    if (sele1 >= 0) {
      ok = SelectorCreateObjectMolecule(G, sele1, name, target, source,
                                        discrete, false, quiet, singletons,
                                        copy_properties);
      if (ok) {
        int sele2 = SelectorIndexByName(G, name);
        ObjectMolecule *old_obj, *new_obj;
        old_obj = SelectorGetFirstObjectMolecule(
            G, sele1); /* get at least one object */
        new_obj = SelectorGetSingleObjectMolecule(G, sele2);

        /* first we need to make sure that the object being moved
           matches the target with respect to both the TTT and the
           object's state matrix (if any) */

        if (old_obj && new_obj) {
          ExecutiveMatrixCopy(G, old_obj->Name, new_obj->Name, 1,
                              1, /* TTT mode */
                              source, target, false, 0, quiet);

          ExecutiveMatrixCopy(G, old_obj->Name, new_obj->Name, 2,
                              2, /* Object state mode */
                              source, target, false, 0, quiet);

          if (copy_properties) {
            // old_obj to new_obj
            if (old_obj->NCSet == new_obj->NCSet) {
              // only if the objects have the same number of CoordSets (for now)
              int state;
              for (state = 0; state < new_obj->NCSet; state++) {
                CoordSet *new_cs = new_obj->CSet[state],
                         *old_cs = old_obj->CSet[state];
                if (old_cs->prop_id) {
                  PropertyCheckUniqueID(G, new_cs);
                  PropertyCopyProperties(G, old_cs->prop_id, new_cs->prop_id);
                }
              }
            }
          }

          ExecutiveDoZoom(G, new_obj, !exists, zoom, true);
        }
      }
    }
  }
  if (ok) {
    return {};
  } else {
    return pymol::make_error("Failed to Create Object");
  }
}
pymol::Result<> ExecutiveCopy(PyMOLGlobals* G, const char* src, const char* dst,
                              int zoom) {
  const pymol::CObject* os = ExecutiveFindObjectByName(G, src);
  if (!os) {
    return pymol::make_error("Object not found.");
  }

  pymol::CObject* oDst = os->clone();

  if (!oDst) {
    return pymol::make_error("Failed to create copy");
  }

  strcpy(oDst->Name, dst);
  ExecutiveManageObject(G, oDst, zoom, false);

  PRINTFB(G, FB_Executive, FB_Actions)
  " Executive: object %s created.\n", oDst->Name ENDFB(G);
  SceneChanged(G);
  return {};
}

pymol::Result<> ExecutiveMove(PyMOLGlobals* G, pymol::zstring_view axis,
                              float dist) {
  switch (axis[0]) {
    case 'x':
      SceneTranslate(G, dist, 0.0, 0.0);
      break;
    case 'y':
      SceneTranslate(G, 0.0, dist, 0.0);
      break;
    case 'z':
      SceneTranslate(G, 0.0, 0.0, dist);
      break;
    default:
      return pymol::make_error("Axis must be x, y, or z");
  }
  return {};
}
/**
 * Implementation of `cmd.identify()`
 *
 * @param s1 atom selection expression
 * @param mode 0 for index only, 1 for (obj, index)
 * @param[out] indexVLA
 * @param[out] objVLA
 * @return number of atoms or -1 on selection error
 */
int ExecutiveIdentifyObjects(PyMOLGlobals* G, const char* s1, int mode,
                             int** indexVLA, ObjectMolecule*** objVLA) {
  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();
  ObjectMoleculeOpRec op2;
  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code = OMOP_IdentifyObjects;
    if (mode != 0) {
      op2.obj1VLA = VLAlloc(ObjectMolecule*, 1000);
    }
    op2.i1VLA = VLAlloc(int, 1000);
    op2.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    VLASize(op2.i1VLA, int, op2.i1);
    if (mode != 0) {
      VLASize(op2.obj1VLA, ObjectMolecule*, op2.i1);
    }
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } else {
    return -1;
  }
  return (op2.i1);
}
int ExecutiveIndex(PyMOLGlobals* G, const char* s1, int mode, int** indexVLA,
                   ObjectMolecule*** objVLA) {
  ObjectMoleculeOpRec op2;

  SelectorTmp tmpsele1(G, s1);
  int sele1 = tmpsele1.getIndex();

  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code = OMOP_Index;
    op2.obj1VLA = VLAlloc(ObjectMolecule*, 1000);
    op2.i1VLA = VLAlloc(int, 1000);
    op2.i1 = 0;
    ExecutiveObjMolSeleOp(G, sele1, &op2);
    VLASize(op2.i1VLA, int, op2.i1);
    VLASize(op2.obj1VLA, ObjectMolecule*, op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } else {
    return -1;  // invalid selection
  }
  return (op2.i1);
}
void ExecutiveUpdateObjectSelection(PyMOLGlobals* G, pymol::CObject* obj) {
  if (obj->type == cObjectMolecule) {
    SelectorUpdateObjectSele(G, (ObjectMolecule*)obj);
  }
}
/**
 * Reset camera view or object TTT matrix. Stores key frames for modified
 * objects if `movie_auto_store=on`.
 *
 * @param name Empty, "all", "same", or object name pattern
 *
 * - empty name: Reset camera view
 * - "all": Reset TTT matrices and store key frames for all objects
 * - "same": Reset TTT matrices and store key frames for objects which currently
 * have any key frames
 * - pattern: Reset TTT matrices and store key frames for objects which match
 * the pattern
 */
pymol::Result<> ExecutiveReset(PyMOLGlobals* G, pymol::zstring_view name) {
  if (name.empty()) {
    SceneResetMatrix(G);
    ExecutiveWindowZoom(G, cKeywordAll, 0.0, -1, 0, 0,
                        true); /* reset does all states */
    return {};
  }

  bool do_reset_all = name == cKeywordAll;
  auto store = SettingGet<bool>(G, cSetting_movie_auto_store);

  /**
   * @param any_spec_level If false, then filter for objects with spec level >=
   * 0
   */
  auto reset_rec = [&](SpecRec& rec, bool any_spec_level = true) {
    pymol::CObject* obj = rec.obj;
    if (rec.type == cExecObject &&
        (any_spec_level || ObjectGetSpecLevel(obj, 0) >= 0)) {
      ObjectResetTTT(obj, store);
      obj->invalidate(cRepNone, cRepInvExtents, -1);
    }
  };

  if (do_reset_all || name == cKeywordSame) {
    for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
      reset_rec(rec, do_reset_all);
    }
  } else {
    for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, name)) {
      reset_rec(rec);
    }
  }

  if (store && SettingGet<bool>(G, cSetting_movie_auto_interpolate)) {
    ExecutiveMotionReinterpolate(G);
  }

  SceneInvalidate(G);

  return {};
}
void ExecutiveDrawNow(PyMOLGlobals* G, ExecutiveDrawInfo execDrawInfo) {
  CExecutive* I = G->Executive;
  OrthoDrawInfo drawInfo{};
  drawInfo.offscreenRender = execDrawInfo.offscreen;
  drawInfo.clearTarget = execDrawInfo.clearTarget;

  if (PyMOL_GetIdleAndReady(G->PyMOL) &&
      !SettingGetGlobal_b(G, cSetting_suspend_deferred))
    OrthoExecDeferred(G);
  if (!SettingGetGlobal_b(G, cSetting_suspend_updates)) {
    int stereo_mode = SettingGetGlobal_i(G, cSetting_stereo_mode);
    int stereo = SettingGetGlobal_i(G, cSetting_stereo);
    if (G->HaveGUI && G->ValidContext && !G->Renderer) {
      glMatrixMode(GL_MODELVIEW); /* why is this necessary?  is it? */
    }

    ExecutiveUpdateSceneMembers(G);
    SceneUpdate(G, false);
    if (WizardUpdate(G)) SceneUpdate(G, false);
    if (stereo) {
      switch (stereo_mode) {
        case cStereo_geowall: {
          int width = G->Option->winX;
          int height = G->Option->winY;
          SceneSetViewport(G, 0, 0, width / 2, height);
          drawInfo.renderMode = OrthoRenderMode::GeoWallLeft;
          OrthoDoDraw(G, drawInfo);
          drawInfo.renderMode = OrthoRenderMode::GeoWallRight;
          OrthoDoDraw(G, drawInfo);
          SceneSetViewport(G, 0, 0, width, height);
        } break;
        default:
          drawInfo.renderMode = OrthoRenderMode::Main;
          OrthoDoDraw(G, drawInfo);
          break;
      }
    } else {
      drawInfo.renderMode = OrthoRenderMode::Main;
      OrthoDoDraw(G, drawInfo);
    }

    if (G->HaveGUI && G->ValidContext) {
      if (I->CaptureFlag) {
        I->CaptureFlag = false;
        SceneCaptureWindow(G);
      }
    }
    PyMOL_NeedSwap(G->PyMOL);
  }

  //  PRINTFD(G, FB_Executive)
  //    " ExecutiveDrawNow: leaving.\n" ENDFD;
}
int ExecutiveCountStates(PyMOLGlobals* G, const char* s1) {
  CExecutive* I = G->Executive;
  int sele1;
  int result = 0;
  int n_state;
  SpecRec* list_rec = nullptr;
  if ((!s1) || (!s1[0])) s1 = cKeywordAll;
  for (const auto& rec : ExecutiveGetSpecRecsFromPattern(G, s1)) {
    switch (rec.type) {
      case cExecAll:
        while (ListIterate(I->Spec, list_rec, next)) {
          if (list_rec->type == cExecObject) {
            n_state = list_rec->obj->getNFrame();
            if (result < n_state) result = n_state;
          }
        }
        break;
      case cExecSelection:
        sele1 = SelectorIndexByName(G, rec.name);
        if (sele1 >= 0) {
          SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);
          n_state = SelectorGetSeleNCSet(G, sele1);
          if (result < n_state) result = n_state;
        }
        break;
      case cExecObject:
        n_state = rec.obj->getNFrame();
        if (result < n_state) result = n_state;
        break;
    }
  }
  return (result);
}
int ExecutiveRay(PyMOLGlobals* G, int width, int height, int mode, float angle,
                 float shift, int quiet, int defer, int antialias) {
  if ((mode == 0) && G->HaveGUI &&
      SettingGetGlobal_b(G, cSetting_auto_copy_images)) {
    /* force deferred behavior if copying image to clipboard */
    defer = 1;
  }

  ExecutiveUpdateSceneMembers(G);

  g_mmol_ray_progress = 0.0f;
  if (defer && (mode == 0)) {
    SceneDeferRay(G, width, height, mode, angle, shift, quiet, true, antialias);
  } else {
    SceneRay(G, width, height, mode, nullptr, nullptr, angle, shift, quiet,
             nullptr, true, antialias);
  }
  g_mmol_ray_progress = -1.0f;
  return 1;
}
const char* ExecutiveFindBestNameMatch(PyMOLGlobals* G, const char* name) {
  const char* result;
  CExecutive* I = G->Executive;
  SpecRec *rec = nullptr, *best_rec = nullptr;
  int best;
  int wm;
  auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);

  best = 0;
  result = name;

  while (ListIterate(I->Spec, rec, next)) {
    wm = WordMatch(G, name, rec->name, ignore_case);
    if (wm < 0) {
      best_rec = rec;
      best = wm;
      break;
    } else if ((best > 0) && (best < wm)) {
      best_rec = rec;
      best = wm;
    }
  }
  if (best_rec) result = best_rec->name;
  return (result);
}

SpecRec* ExecutiveFindSpec(PyMOLGlobals* G, pymol::zstring_view name_view) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  const char* name = name_view.c_str();
  // ignore % prefix
  if (name[0] && name[0] == '%') name++;
  { /* first, try for perfect, case-specific match */
    OVreturn_word result;
    if (OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, name)))) {
      auto keyRes = I->Key.find(result.word);
      if (keyRes != I->Key.end()) {
        if (!TrackerGetCandRef(I->Tracker, keyRes->second,
                               (TrackerRef**)(void*)&rec)) {
          rec = nullptr;
        }
      }
    }
    if (!rec) { /* otherwise try partial/case-nonspecific match */
      rec = ExecutiveAnyCaseNameMatch(G, name);
    }
  }
  return (rec);
}
bool ExecutiveObjMolSeleOp(PyMOLGlobals* G, int sele, ObjectMoleculeOpRec* op) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  ObjectMolecule* obj = nullptr;
  int update_table = true;

  /* if we're given a valid selection */
  if (sele >= 0) {
    /* iterate over all the objects in the global list */
    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (rec->obj->type == cObjectMolecule) {
          /* if the objects are valid molecules, then perform the operation in
           * op_code */
          obj = (ObjectMolecule*)rec->obj;

          switch (op->code) {
            case OMOP_RenameAtoms: {
              int result =
                  SelectorRenameObjectAtoms(G, obj, sele, op->i2, update_table);
              if (result > 0) op->i1 += result;
              update_table = false;
            } break;
            default:
              /* all other cases, perform the operation on obj */
              if (!ObjectMoleculeSeleOp(obj, sele, op)) {
                return false;
              }
              break;
          }
        }
      }
    }
  }

  return true;
}
/**
 * Full screen state fallback in case we can't get get the state from
 * the window manager.
 */
static bool _is_full_screen = false;
/**
 * Get the fullscreen state from the window manager or return -1 if
 * not available.
 */
bool ExecutiveIsFullScreen(PyMOLGlobals* G) {
  if (!G->HaveGUI || !G->ValidContext) return false;

  int flag = -1;

  PRINTFD(G, FB_Executive)
  " %s: flag=%d fallback=%d.\n", __func__, flag, _is_full_screen ENDFD;

  if (flag > -1) return flag;
  return _is_full_screen;
}
void ExecutiveFullScreen(PyMOLGlobals* G, int flag) {
  if (!G->HaveGUI) return;

  int wm_flag = ExecutiveIsFullScreen(G);

  if (flag < 0) {
    flag = !wm_flag;
  }

  _is_full_screen = (flag != 0);

  PyMOL_NeedReshape(G->PyMOL, flag, 0, 0, 0, 0);
  SceneChanged(G);
}

int ExecutiveCheckGroupMembership(PyMOLGlobals* G, int list_id,
                                  pymol::CObject* obj) {
  CExecutive* I = G->Executive;
  int result = false;
  CTracker* I_Tracker = I->Tracker;
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  if (iter_id) {
    SpecRec* rec = nullptr;
    while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                     (TrackerRef**)(void*)&rec)) {
      if (rec && (rec->type == cExecObject) && (rec->obj == obj)) {
        result = true;
        break;
      }
    }
    TrackerDelIter(I_Tracker, iter_id);
  }
  return result;
}

int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals* G,
                                             const char* name) {
  CExecutive* I = G->Executive;
  int result = 0;
  CWordMatcher* matcher;
  CWordMatchOptions options;
  CTracker* I_Tracker = I->Tracker;
  const char* wildcard = SettingGetGlobal_s(G, cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  SpecRec* rec;

  WordMatchOptionsConfigNameList(&options, *wildcard,
                                 SettingGetGlobal_b(G, cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if (matcher) {
    if (iter_id) {
      while ((cand_id = TrackerIterNextCandInList(I_Tracker, iter_id,
                                                  (TrackerRef**)(void*)&rec))) {
        if (rec && !(rec->type == cExecAll)) {
          if (WordMatcherMatchAlpha(matcher, rec->name)) {
            if ((rec->type == cExecObject) &&
                (rec->obj->type == cObjectGroup)) {
              if (!result) result = TrackerNewList(I_Tracker, nullptr);
              if (result) {
                TrackerLink(I_Tracker, cand_id, result, 1);
              }
            }
          }
        }
      }
    }
  } else if ((rec = ExecutiveFindSpec(G, name))) { /* only one name in list */
    if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
      result = TrackerNewList(I_Tracker, nullptr);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  } else if ((rec = ExecutiveUnambiguousNameMatch(G, name))) {
    if ((rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
      result = TrackerNewList(I_Tracker, nullptr);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  }
  if (matcher) WordMatcherFree(matcher);
  if (iter_id) TrackerDelIter(I->Tracker, iter_id);
  if (result) ExecutiveExpandGroupsInList(G, result, cExecExpandGroups);
  return result;
}

int ExecutiveGetExpandedGroupList(PyMOLGlobals* G, const char* name) {
  CExecutive* I = G->Executive;
  int result = 0;
  int list_id = 0;

  SpecRec* rec = ExecutiveFindSpec(G, name);
  ExecutiveUpdateGroups(G, false);
  if (rec && (rec->type == cExecObject) && (rec->obj->type == cObjectGroup)) {
    list_id = rec->group_member_list_id;
  }
  if (list_id) {
    result = TrackerNewListCopy(I->Tracker, list_id, nullptr);
    ExecutiveExpandGroupsInList(G, result, cExecExpandGroups);
  }
  return result;
}

void ExecutiveFreeGroupList(PyMOLGlobals* G, int list_id) {
  CExecutive* I = G->Executive;
  TrackerDelList(I->Tracker, list_id);
}
pymol::CObject* ExecutiveFindObjectByName(PyMOLGlobals* G,
                                          pymol::zstring_view name) {
  pymol::CObject* obj = nullptr;
  SpecRec* rec = ExecutiveFindSpec(G, name.c_str());
  if (rec && (rec->type == cExecObject)) {
    obj = rec->obj;
  }
  return (obj);
}
/* returns nullptr if none found */
pymol::CObject** ExecutiveFindObjectsByType(PyMOLGlobals* G, int objType) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  int n = 0;
  pymol::CObject** rVal = VLAlloc(pymol::CObject*, 1);

  /* loop over all known objects */
  while (ListIterate(I->Spec, rec, next)) {
    /* make sure it exists and is the right type */
    if (rec->obj && rec->type == cExecObject) {
      if (rec->obj->type == objType) {
        /* this could be optimized */
        VLACheck(rVal, pymol::CObject*, n);
        *(rVal + n) = rec->obj;
        n++;
      }
    }
  }
  VLASize(rVal, pymol::CObject*, n);
  if (n == 0) {
    VLAFree(rVal);
    return nullptr;
  } else
    return rVal;
}
Block* ExecutiveGetBlock(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  return (I);
}

void ExecutivePurgeSpec(PyMOLGlobals* G, SpecRec* rec, bool save) {
  CExecutive* I = G->Executive;

  if (!save) {
    CGOFree(rec->gridSlotSelIndicatorsCGO);
  }

  ExecutiveInvalidateGroups(G, false);
  ExecutiveInvalidatePanelList(G);
  switch (rec->type) {
    case cExecObject:
      if (I->LastEdited == rec->obj) I->LastEdited = nullptr;
      if (rec->obj->type == cObjectMolecule)
        if (EditorIsAnActiveObject(G, (ObjectMolecule*)rec->obj))
          EditorInactivate(G);
      SeqChanged(G);
      if (rec->visible) {
        SceneObjectDel(G, rec->obj, false);
        ExecutiveInvalidateSceneMembers(G);
      }
      ExecutiveDelKey(I, rec);
      SelectorDelete(G, rec->name);
      if (!save) {
        DeleteP(rec->obj);
      }
      TrackerDelCand(I->Tracker, rec->cand_id);
      break;
    case cExecSelection:
      if (rec->visible) {
        SceneInvalidate(G);
        SeqDirty(G);
      }
      ExecutiveDelKey(I, rec);
      SelectorDelete(G, rec->name);
      TrackerDelCand(I->Tracker, rec->cand_id);
      break;
  }
}
pymol::Result<std::vector<DiscardedRec>> ExecutiveDelete(
    PyMOLGlobals* G, pymol::zstring_view nameView, bool save) {
  std::vector<DiscardedRec> discardedRecs;
  CExecutive* I = G->Executive;
  auto name = nameView.c_str();
  std::vector<OrderRec> specPositions;
  if (save) {
    specPositions = ExecutiveGetOrderOf(G, nameView);
  }

  auto getRecPos = [&](SpecRec* rec) {
    auto it =
        std::find_if(specPositions.begin(), specPositions.end(),
                     [&](OrderRec& oRec) { return oRec.name == rec->name; });
    return it == specPositions.end() ? -1 : it->pos;
  };
  auto deleteObjRec = [&](SpecRec* rec) {
    if (save) {
      if (rec->obj->type == cObjectGroup) {
        ExecutiveGroupPurge(G, rec, &discardedRecs);
      }
      ExecutivePurgeSpec(G, rec, save);
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      SceneObjectRemove(G, rec->obj);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      if (rec->obj->type == cObjectGroup) {
        /* remove members of the group */
        ExecutiveGroup(G, rec->name, "", cExecutiveGroupPurge, true);
      }
      ExecutivePurgeSpec(G, rec, save);
      ListDelete(I->Spec, rec, next,
                 SpecRec); /* NOTE: order N in list length! TO FIX */
    }
  };

  auto deleteSeleRec = [&](SpecRec* rec) {
    ExecutivePurgeSpec(G, rec, save);
    if (save) {
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      ListDelete(I->Spec, rec, next,
                 SpecRec); /* NOTE: order N in list length! TO FIX */
    }
  };

  SpecRec* rec = nullptr;
  CTracker* I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, false, false);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                   (TrackerRef**)(void*)&rec)) {
    if (rec) {
      switch (rec->type) {
        case cExecSelection:
          deleteSeleRec(rec);
          break;
        case cExecObject:
          deleteObjRec(rec);
          break;
        case cExecAll:
          rec = nullptr;
          while (ListIterate(I->Spec, rec, next)) {
            switch (rec->type) {
              case cExecAll:
                break;
              case cExecObject:
                deleteObjRec(rec);
                rec = nullptr;
                break;
              case cExecSelection:
                deleteSeleRec(rec);
                rec = nullptr;
                break;
            }
          }
          SelectorDefragment(G);
          break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);

  // fix for PYMOL-757 - Note: This should probably go somewhere else, but we
  // couldn't figure out the correct place so far
  ExecutiveUpdateGroups(G, false);
  return discardedRecs;
}

pymol::Result<> ExecutiveDeleteStates(PyMOLGlobals* G, std::string_view name,
                                      const std::vector<int>& states) {
  for (auto& rec : ExecutiveGetSpecRecsFromPattern(G, name.data())) {
    if (rec.type != cExecObject) {
      continue;
    }
    auto* obj = rec.obj;
    if (obj->type != cObjectMolecule) {
      continue;
    }
    auto* mol = static_cast<ObjectMolecule*>(obj);
    if (mol->DiscreteFlag) {
      G->Feedback->addColored(
          "Error: Cannot delete states from discrete objects.\n", FB_Warnings);
      continue;
    }
    ObjectMoleculeDeleteStates(mol, states);
  }
  SceneChanged(G);
  ExecutiveInvalidatePanelList(G);
  return {};
}

void ExecutiveDump(PyMOLGlobals* G, const char* fname, const char* obj,
                   int state, int quiet) {
  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;

  SceneUpdate(G, false);

  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) {
      if (strcmp(rec->obj->Name, obj) == 0) break;
    }
  }
  if (rec) {
    if (rec->obj->type == cObjectMesh) {
      ObjectMeshDump((ObjectMesh*)rec->obj, fname, state, quiet);
    } else if (rec->obj->type == cObjectSurface) {
      ObjectSurfaceDump((ObjectSurface*)rec->obj, fname, state, quiet);
    } else if (rec->obj->type == cObjectMap) {
      ObjectMapDump((ObjectMap*)rec->obj, fname, state, quiet);
    } else {
      ErrMessage(G, __func__, "Invalid object type for this operation.");
    }
  } else {
    ErrMessage(G, __func__, "Object not found.");
  }
}
void ExecutiveMemoryDump(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  fprintf(stderr, " Executive: %d candidate(s) %d list(s) %d link(s).\n",
          TrackerGetNCand(I->Tracker), TrackerGetNList(I->Tracker),
          TrackerGetNLink(I->Tracker));
}

static void ExecutiveDoAutoGroup(PyMOLGlobals* G, SpecRec* rec) {
  CExecutive* I = G->Executive;
  int auto_mode = SettingGetGlobal_i(G, cSetting_group_auto_mode);
  if (auto_mode && (rec->name[0] != '_')) {
    auto ignore_case = SettingGet<bool>(G, cSetting_ignore_case);
    char* period = rec->name + strlen(rec->name);
    SpecRec* found_group = nullptr;
    WordType seek_group_name;
    UtilNCopy(seek_group_name, rec->name, sizeof(WordType));

    while ((period > rec->name) && (!found_group)) {
      period--;
      if (*period == '.') {
        seek_group_name[period - rec->name] = 0;
        {
          SpecRec* group_rec = nullptr;
          while (ListIterate(I->Spec, group_rec, next)) {
            if ((group_rec->type == cExecObject) &&
                (group_rec->obj->type == cObjectGroup)) {
              if (WordMatchExact(G, group_rec->name, seek_group_name,
                                 ignore_case)) {
                found_group = group_rec;
                strcpy(rec->group_name, seek_group_name);
                break;
              }
            }
          }
        }

        if ((!found_group) && (auto_mode == 2)) {
          auto obj = new ObjectGroup(G);
          if (obj) {
            ObjectSetName(obj, seek_group_name);
            strcpy(rec->group_name, obj->Name);
            ExecutiveManageObject(G, obj, false, true);
            ExecutiveInvalidateGroups(G, false);
            break;
          }
        }
      }
    }
    if (found_group) ExecutiveInvalidateGroups(G, false);
  }
}
/**
 * Manages an object. Adds it to the list of spec records and related trackers,
 * and to the scene for rendering.
 *
 * Also handles:
 * - auto_dss
 * - group_auto_mode
 * - auto_defer_builds
 *
 * If the object is already managed, then only do `auto_dss` and
 * `auto_defer_builds`.
 *
 * If an object with the same name exists, then delete it and re-use the
 * existing spec rec to manage the new object.
 *
 * @param obj Object to manage. Executive takes ownership.
 * @param zoom Zoom the camera, see ExecutiveDoZoom for valid values.
 */
void ExecutiveManageObject(PyMOLGlobals* G, pymol::CObject* obj, int zoom,
                           int quiet) {
  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;
  int exists = false;
  int previousVisible;
  int previousObjType = 0;

  if (SettingGetGlobal_b(G, cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->obj == obj) {
      exists = true;
    }
  }
  if (!exists) {
    if (WordMatchExact(G, cKeywordAll, obj->Name, true)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " Executive: object name \"%s\" is illegal -- renamed to 'all_'.\n",
          obj->Name ENDFB(G);
      strcat(obj->Name, "_"); /* don't allow object named "all" */
    } else if (SelectorNameIsKeyword(G, obj->Name)) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " Executive-Warning: name \"%s\" collides with a selection language "
      "keyword.\n",
          obj->Name ENDFB(G);
    }

    while (ListIterate(I->Spec, rec, next)) {
      if (rec->type == cExecObject) {
        if (strcmp(rec->obj->Name, obj->Name) == 0) break;
      }
    }
    if (rec) { /* another object of this type already exists */
      /* purge it */
      SceneObjectDel(G, rec->obj, false);
      ExecutiveInvalidateSceneMembers(G);
      previousObjType = rec->obj->type;
      DeleteP(rec->obj);
    } else {
      if (!quiet)
        if (obj->Name[0] != '_') { /* suppress internal objects */
          PRINTFB(G, FB_Executive, FB_Actions)
          " Executive: object \"%s\" created.\n", obj->Name ENDFB(G);
        }
    }
    if (!rec) ListElemCalloc(G, rec, SpecRec);

    strcpy(rec->name, obj->Name);
    rec->type = cExecObject;
    rec->obj = obj;
    previousVisible = rec->visible;
    if (previousObjType == rec->obj->type) {
      // skip
    } else if (rec->obj->type == cObjectMap) {
      rec->visible = 0;
    } else {
      rec->visible = 1;
      /*      SceneObjectAdd(G,obj); */
    }
    if (previousVisible != rec->visible) {
      ReportEnabledChange(G, rec);
    }

    // Is this a new rec?
    if (!rec->cand_id) {
      rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef*)(void*)rec);

      TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
      TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id, 1);
      ListAppend(I->Spec, rec, next, SpecRec);
      ExecutiveAddKey(I, rec);
      ExecutiveInvalidatePanelList(G);

      ExecutiveDoAutoGroup(G, rec);
    }

    if (rec->visible) {
      rec->in_scene = SceneObjectAdd(G, obj);
      ExecutiveInvalidateSceneMembers(G);
    }
  }

  ExecutiveUpdateObjectSelection(G, obj);

  if (SettingGetGlobal_b(G, cSetting_auto_dss)) {
    if (obj->type == cObjectMolecule) {
      ObjectMolecule* objMol = (ObjectMolecule*)obj;
      if (objMol->NCSet == 1) {
        ExecutiveAssignSS(G, obj->Name, 0, nullptr, true, objMol, true);
      }
    }
  }

  {
    int n_state = obj->getNFrame();
    int defer_limit = SettingGetGlobal_i(G, cSetting_auto_defer_builds);
    if ((defer_limit >= 0) && (n_state >= defer_limit)) {
      int defer_builds = SettingGetGlobal_b(G, cSetting_defer_builds_mode);
      if (!defer_builds) SettingSetGlobal_b(G, cSetting_defer_builds_mode, 3);
    }
  }

  ExecutiveDoZoom(G, obj, !exists, zoom, true);

  SeqChanged(G);
  OrthoInvalidateDoDraw(G);
}
void ExecutiveManageSelection(PyMOLGlobals* G, const char* name) {
  SpecRec* rec = nullptr;
  CExecutive* I = G->Executive;
  bool hidden_name = name[0] == '_';
  bool hide_all =
      !hidden_name && (SettingGet<bool>(G, cSetting_active_selections) ||
                       SettingGet<bool>(G, cSetting_auto_hide_selections));

  for (auto& rec_ : pymol::make_list_adapter(I->Spec)) {
    if (rec_.type == cExecSelection) {
      if (!rec && pymol::zstring_view(rec_.name) == name) {
        rec = &rec_;
      } else if (hide_all) {
        rec_.setEnabled(G, false);
      }
    }
  }

  if (!rec) {
    ListElemCalloc(G, rec, SpecRec);
    strcpy(rec->name, name);
    rec->type = cExecSelection;
    rec->next = nullptr;
    rec->sele_color = -1;
    assert(!rec->visible);
    rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef*)(void*)rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
    TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id, 1);
    ListAppend(I->Spec, rec, next, SpecRec);
    ExecutiveAddKey(I, rec);
    ExecutiveInvalidatePanelList(G);
  }
  assert(rec);
  if (!hidden_name) {
    if (SettingGetGlobal_b(G, cSetting_auto_show_selections) && !rec->visible) {
      rec->visible = true;
      ReportEnabledChange(G, rec);
    }
  }
  if (rec->visible) SceneInvalidate(G);
  ExecutiveDoAutoGroup(G, rec);
  SeqDirty(G);
}
void ExecutiveSpecEnable(PyMOLGlobals* G, SpecRec* rec, int parents, int log) {
  if (log && SettingGetGlobal_b(G, cSetting_logging)) {
    OrthoLineType buffer = "";
    sprintf(buffer, "cmd.enable('%s',%d)", rec->obj->Name, parents);
    PLog(G, buffer, cPLog_pym);
  }

  if (!rec->visible) {
    rec->visible = true;
    ReportEnabledChange(G, rec);
  }
  if (!rec->in_scene) {
    rec->in_scene = SceneObjectAdd(G, rec->obj);
  }

  if (parents) {
    CExecutive* I = G->Executive;
    CTracker* I_Tracker = I->Tracker;
    int list_id = ExecutiveGetObjectParentList(G, rec);
    if (list_id) {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec* parent_rec = nullptr;

      while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                       (TrackerRef**)(void*)&parent_rec)) {
        if (rec) {
          switch (parent_rec->type) {
            case cExecObject:
              if (!parent_rec->in_scene) {
                parent_rec->in_scene = SceneObjectAdd(G, parent_rec->obj);
              }
              if (!parent_rec->visible) {
                parent_rec->visible = true;
                ReportEnabledChange(G, parent_rec);
              }
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
    TrackerDelList(I_Tracker, list_id);
  }
  ExecutiveInvalidateSceneMembers(G);
}

void ExecutiveSpecSetVisibility(PyMOLGlobals* G, SpecRec* rec, int new_vis,
                                int mod, int parents) {
  std::string buffer;
  int logging = SettingGetGlobal_i(G, cSetting_logging);
  if (rec->type == cExecObject) {
    if (rec->visible && !new_vis) {
      if (logging)
        buffer = pymol::string_format("cmd.disable('%s')", rec->obj->Name);
      SceneObjectDel(G, rec->obj, true);
      ExecutiveInvalidateSceneMembers(G);
      if (rec->visible != new_vis) {
        rec->visible = new_vis;
        ReportEnabledChange(G, rec);
      }
    } else if ((!rec->visible) && new_vis) {
      ExecutiveSpecEnable(G, rec, parents, logging);
      /*
         if(logging)
         sprintf(buffer,"cmd.enable('%s')",rec->obj->Name);
         rec->in_scene = SceneObjectAdd(G,rec->obj);
         rec->visible=new_vis;
         ExecutiveInvalidateSceneMembers(G);
       */
    }
    SceneChanged(G);
    if (logging && buffer[0]) {
      PLog(G, buffer, cPLog_pym);
    }
  } else if (rec->type == cExecAll) {
    if (SettingGetGlobal_i(G, cSetting_logging)) {
      if (rec->visible)
        buffer = "cmd.disable('all')";
      else
        buffer = "cmd.enable('all')";
      PLog(G, buffer, cPLog_pym);
    }
    ExecutiveSetObjVisib(G, cKeywordAll, !rec->visible, false);
  } else if (rec->type == cExecSelection) {
    if (mod & cOrthoCTRL) {
      buffer = pymol::string_format("cmd.enable('%s')", rec->name);
      PLog(G, buffer, cPLog_pym);
      if (!rec->visible) {
        rec->visible = true;
        ReportEnabledChange(G, rec);
      }
    } else {
      if (rec->visible && !new_vis) {
        if (SettingGetGlobal_i(G, cSetting_logging))
          buffer = pymol::string_format("cmd.disable('%s')", rec->name);
      } else if ((!rec->visible) && new_vis) {
        buffer = pymol::string_format("cmd.enable('%s')", rec->name);
      }
      if (new_vis && SettingGetGlobal_b(G, cSetting_active_selections)) {
        ExecutiveHideSelections(G);
      }
      if (SettingGetGlobal_i(G, cSetting_logging)) {
        PLog(G, buffer, cPLog_pym);
      }
      if (rec->visible != new_vis) {
        rec->visible = new_vis;
        ReportEnabledChange(G, rec);
      }
    }
    SceneChanged(G);
  }
}
int ExecutiveIterateObjectMolecule(PyMOLGlobals* G, ObjectMolecule** obj,
                                   void** hidden) {
  SpecRec **rec = (SpecRec**)hidden, *I_Spec = G->Executive->Spec;
  while (ListIterate(I_Spec, (*rec), next)) {
    if (((*rec)->type == cExecObject) && ((*rec)->obj->type == cObjectMolecule))
      break;
  }
  if (*rec)
    (*obj) = (ObjectMolecule*)(*rec)->obj;
  else
    (*obj) = nullptr;
  return ((*rec) != nullptr);
}
pymol::Result<> ExecutiveReinitialize(PyMOLGlobals* G, int what,
                                      pymol::zstring_view inPattern) {
  CExecutive* I = G->Executive;
  int blocked = false;
  /* reinitialize PyMOL */
  const char* pattern = inPattern.data();
  if (what == 2) pattern = nullptr;

  if (pattern && (!pattern[0])) pattern = nullptr;
  if (!pattern) {
    switch (what) {
      case 0: /* everything */
        ExecutiveDelete(G, cKeywordAll);
        ColorReset(G);
        SettingInitGlobal(G, false, false, true);
        ColorUpdateFrontFromSettings(G);
        MovieReset(G);
        EditorInactivate(G);
        ControlRock(G, 0);
        OrthoReshape(G, -1, -1, false);
        MovieScenesInit(G);

        blocked = PAutoBlock(G);
        PRunStringInstance(G, "cmd.view('*','clear')");
        PRunStringInstance(G, "cmd.config_mouse(\"three_button\")");
        WizardSet(G, nullptr, false);
        PAutoUnblock(G, blocked);

        SculptCachePurge(G);
        SceneReinitialize(G);
        SelectorReinit(G);
        SeqChanged(G);
        break;
      case 1: /* settings */
        SettingInitGlobal(G, false, false, true);
        ExecutiveRebuildAll(G);
        break;
      case 2: /* store_defaults */
        SettingStoreDefault(G);
        break;
      case 3: /* original_settings */
        SettingInitGlobal(G, false, false, false);
        ExecutiveRebuildAll(G);
        break;
      case 4: /* purge_defaults */
        SettingPurgeDefault(G);
        break;
        /* reinitialize is called with what + 5 to reset internal_gui if
         * necessary (PYMOL-1425) */
      case 5:
      case 6:
        if (G->Default) {
          SettingSetGlobal_i(
              G, cSetting_internal_gui,
              SettingGet_i(G, G->Default, nullptr, cSetting_internal_gui));
          SettingGenerateSideEffects(G, cSetting_internal_gui, nullptr, -1, 0);
        }
        break;
    }
    SceneUpdateStereo(G);
  } else {
    {
      CTracker* I_Tracker = I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G, pattern, true, true);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec* rec;

      while (TrackerIterNextCandInList(I_Tracker, iter_id,
                                       (TrackerRef**)(void*)&rec)) {
        if (rec) {
          switch (rec->type) {
            case cExecObject:
              switch (what) {
                case 0:
                case 1:
                  if (rec->obj->Setting.get()) {
                    ObjectPurgeSettings(rec->obj);
                    rec->obj->invalidate(cRepAll, cRepInvAll, -1);
                    SceneInvalidate(G);
                    SeqChanged(G);
                  }
                  break;
              }
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
    }

    /* to do */
  }

  return {};
}
int ExecutiveInit(PyMOLGlobals* G) {
  CExecutive* I = nullptr;
  if ((I = (G->Executive = new CExecutive(G)))) {
    SpecRec* rec = nullptr;

    ListInit(I->Spec);
    I->Tracker = TrackerNew(G);
    I->all_names_list_id = TrackerNewList(I->Tracker, nullptr);
    I->all_obj_list_id = TrackerNewList(I->Tracker, nullptr);
    I->all_sel_list_id = TrackerNewList(I->Tracker, nullptr);
    I->active = true;
    OrthoAttach(G, I, cOrthoTool);

    I->Lex = OVLexicon_New(G->Context->heap);

    ColorSetRampFns({
        [](PyMOLGlobals* G, const char* name) {
          return ExecutiveFindObject<ObjectGadgetRamp>(G, name);
        },
        ObjectGadgetRampInterVertex,
        [](PyMOLGlobals* G) {
          ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
          SceneChanged(G);
        },
    });

    /* create "all" entry */

    ListElemCalloc(G, rec, SpecRec);

    strcpy(rec->name, cKeywordAll);
    rec->type = cExecAll;
    rec->visible = true;
    rec->next = nullptr;
    rec->cand_id = TrackerNewCand(I->Tracker, (TrackerRef*)(void*)rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id, 1);
    ListAppend(I->Spec, rec, next, SpecRec);
    ExecutiveAddKey(I, rec);

    return 1;
  } else
    return 0;
}
void ExecutiveFree(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  SpecRec* rec = nullptr;
  CGOFree(I->selIndicatorsCGO);
  while (ListIterate(I->Spec, rec, next)) {
    if (rec->type == cExecObject) DeleteP(rec->obj);
  }
  ListFree(I->Spec, next, SpecRec);
  if (I->Tracker) TrackerFree(I->Tracker);
  OVLexicon_DEL_AUTO_NULL(I->Lex);

  ExecutiveUniqueIDAtomDictInvalidate(G);

  DeleteP(G->Executive);
}

/**
 * Get the coord set for the given object name and state index. If "omp" is
 * not nullptr, then also store a pointer to the object molecule.
 */
CoordSet* ExecutiveGetCoordSet(PyMOLGlobals* G, const char* name, int state,
                               ObjectMolecule** omp) {
  ObjectMolecule* om = nullptr;
  CoordSet* cs = nullptr;

  ok_assert(1, om = ExecutiveFindObject<ObjectMolecule>(G, name));
  ok_assert(1, cs = om->getCoordSet(state));

ok_except1:
  if (omp != nullptr) *omp = om;
  return cs;
}

/**
 * Like ExecutiveGetProperty() but also return the type code.
 *
 * @return A [type, value] pair, or nullptr if there is no such property
 */
CPythonVal* ExecutiveGetPropertyForObject(PyMOLGlobals* G, const char* propname,
                                          const char* objname, int state,
                                          int quiet) {
  CPythonVal* result = nullptr;
  ObjectMolecule* obj = nullptr;
  obj = ExecutiveFindObjectMoleculeByName(G, objname);
  if (obj) {
    result = ObjectMoleculeGetProperty(obj, propname, state, quiet);
  } else {
    PRINTFB(G, FB_Executive, FB_Errors)
    "ExecutiveGetProperty: objname='%s' not found\n", objname ENDFB(G);
  }
  return (result);
}

/**
 * @return The property value, or nullptr if there is no such property
 */
CPythonVal* ExecutiveGetProperty(PyMOLGlobals* G, const char* propname,
                                 const char* objname, int state, int quiet) {
  CPythonVal* result = nullptr;
  auto* withtype =
      ExecutiveGetPropertyForObject(G, propname, objname, state, quiet);
  if (withtype) {
    result = PyList_GetItem(withtype, 1);
    CPythonVal_Free(withtype);
  }
  return result;
}

pymol::Result<> ExecutiveSetPropertyForObject(PyMOLGlobals* G,
                                              const char* propname,
                                              CPythonVal* value,
                                              const char* objname, int state,
                                              int proptype, int quiet) {
  ObjectMolecule* obj = nullptr;
  int ok = true;

  auto vla = ExecutiveGetObjectMoleculeVLA(G, objname);
  if (!vla) return pymol::make_error("No object(s) found");
  if (ok) {
    int nObj = VLAGetSize(vla);
    for (int a = 0; a < nObj; a++) {
      obj = vla[a];
      ok &= ObjectMoleculeSetProperty((ObjectMolecule*)obj, propname, value,
                                      static_cast<PropertyType>(proptype),
                                      state, quiet);
    }
  }
  if (!ok) {
    return pymol::make_error("Error while setting properties...");
  }
  return {};
}

int ExecutiveSetAtomPropertyForSelection(PyMOLGlobals* G, const char* propname,
                                         CPythonVal* value, const char* s1,
                                         int state, int proptype, int quiet) {
  int sele1, result, ok = true;
  sele1 = SelectorIndexByName(G, s1);

  if (state < 0) state = 0;

  {
    int unblock = PAutoBlock(G); /*   PBlock(G);    PBlockAndUnlockAPI(); */
    if (sele1 >= 0) {
      result = SelectorSetAtomPropertyForSelection(
          G, sele1, propname, value, static_cast<PropertyType>(proptype), state,
          quiet);
      if (!quiet && !result) {
        PRINTFB(G, FB_Executive, FB_Warnings)
        "ExecutiveSetAtomPropertyForSelection: no atoms set\n" ENDFB(G);
      }
    }
    if (PyErr_Occurred()) {
      PyErr_Print();
      ok = false;
    }
    PAutoUnblock(G, unblock); /*    PUnblock(G);  PLockAPIAndUnblock(); */
  }
  return (ok);
}

void ExecutiveUniqueIDAtomDictInvalidate(PyMOLGlobals* G) {
  CExecutive* I = G->Executive;
  I->m_id2eoo.clear();
  I->m_eoo.clear();
}

const ExecutiveObjectOffset* ExecutiveUniqueIDAtomDictGet(PyMOLGlobals* G,
                                                          int i) {
  CExecutive* I = G->Executive;

  if (I->m_eoo.empty()) {
    ExecutiveGetUniqueIDAtomVLADict(G);
  }

  auto offsetIt = I->m_id2eoo.find(i);
  return offsetIt == I->m_id2eoo.end() ? nullptr : &I->m_eoo[offsetIt->second];
}

pymol::Result<> ExecutiveSetFeedbackMask(PyMOLGlobals* G, int action,
                                         unsigned int sysmod,
                                         unsigned char mask) {
  switch (action) {
    case 0:
      G->Feedback->setMask(sysmod, mask);
      break;
    case 1:
      G->Feedback->enable(sysmod, mask);
      break;
    case 2:
      G->Feedback->disable(sysmod, mask);
      break;
    case 3:
      G->Feedback->push();
      break;
    case 4:
      G->Feedback->pop();
      break;
  }
  return {};
}

/**
 * Set the `visible` property and emit `ReportEnabledChange` if the property
 * changed. If it doesn't change, do nothing.
 */
void SpecRec::setEnabled(PyMOLGlobals* G, bool enabled) {
  if (enabled != visible) {
    visible = enabled;
    ReportEnabledChange(G, this);
  }
}

std::vector<OrderRec> ExecutiveGetOrderOf(PyMOLGlobals* G,
                                          pymol::zstring_view nameListView) {
  auto I = G->Executive;
  std::vector<OrderRec> recs;
  for (auto& rec :
       ExecutiveGetSpecRecsFromPattern(G, nameListView, true, false)) {
    auto pos = ListGetPosition(I->Spec, &rec);
    recs.emplace_back(rec.name, *pos);
  }
  auto sort_by_pos = [](const OrderRec& a, const OrderRec& b) {
    return a.pos < b.pos;
  };
  std::sort(recs.begin(), recs.end(), sort_by_pos);
  return recs;
}

pymol::Result<> ExecutiveBackgroundColor(PyMOLGlobals* G,
                                         pymol::zstring_view color) {
  SettingSet_color(G->Setting, cSetting_bg_rgb, color.c_str());
  SettingGenerateSideEffects(G, cSetting_bg_rgb, nullptr, -1, 0);
  return {};
}

pymol::Result<> ExecutiveBackgroundColor(PyMOLGlobals* G, int colorIdx) {
  SettingSet_i(G->Setting, cSetting_bg_rgb, colorIdx);
  SettingGenerateSideEffects(G, cSetting_bg_rgb, nullptr, -1, 0);
  return {};
}

/**
 * Updates an object's dependencies.
 * @param obj the object to update
 * @param visited a set of objects that have already been visited
 */
static void ExecutiveUpdateObjectsImpl(
    PyMOLGlobals* G, const pymol::CObject& obj,
    std::unordered_set<const pymol::CObject*>& visited) {
  if (visited.count(&obj)) {
    return;
  }
  auto& deps = G->Executive->m_objDeps[&obj];
  deps.clear();
  visited.insert(&obj);
  auto is_mol_in_dist = [&](const auto& mol, const auto& dist) {
    for (const auto ds : dist->DSet) {
      const auto ds_deps = ds->getDependentObjects();
      if (ds_deps.count(mol)) {
        return true;
      }
    }
    return false;
  };

  // TODO: This only detects mol <--> dist dependencies.
  for (auto& rec : pymol::make_list_adapter(G->Executive->Spec)) {
    if (rec.type != cExecObject || rec.obj == &obj) {
      continue;
    }
    switch (rec.obj->type) {
      case cObjectMolecule: {
        auto mol = static_cast<ObjectMolecule*>(rec.obj);
        if (auto measure = dynamic_cast<const ObjectDist*>(&obj)) {
          if (is_mol_in_dist(mol, measure)) {
            deps.insert(rec.obj);
            ExecutiveUpdateObjectsImpl(G, *rec.obj, visited);
          }
        }
      } break;
      case cObjectMeasurement: {
        auto measure = static_cast<ObjectDist*>(rec.obj);
        if (auto mol = dynamic_cast<const ObjectMolecule*>(&obj)) {
          if (is_mol_in_dist(mol, measure)) {
            deps.insert(rec.obj);
            ExecutiveUpdateObjectsImpl(G, *rec.obj, visited);
          }
        }
      } break;
      default:
        break;
    }
  }

  // Erase self
  deps.erase(&obj);
}
