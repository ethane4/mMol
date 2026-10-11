
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
#include "pymol/platform/os_python.h"

#include "pymol/platform/os_proprietary.h"

#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_std.h"

#include "pymol/base/base.h"
#include "pymol/base/err.h"
#include "pymol/platform/memory_debug.h"
#include "pymol/settings/setting.h"
#include "pymol/settings/color_index.h"
#include "pymol/base/word.h"
#include "pymol/base/util.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/base/pymol_options.h"
#include "pymol/base/parse.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"

#include <algorithm>
#include <thread>

/**
 * Setting level info table
 *
 * Levels are not hierarchical at the atom/bond level, that's why a simple
 * sorted enumeration is not sufficient.
 *
 *     global < object < object-state
 *                       object-state < atom < atom-state
 *                       object-state < bond < bond-state
 */
const SettingLevelInfoType SettingLevelInfo[] = {
    {"unused", 0x00},        // 0b00000000
    {"global", 0x00},        // 0b00000000
    {"object", 0x01},        // 0b00000001
    {"object-state", 0x03},  // 0b00000011
    {"atom", 0x07},          // 0b00000111
    {"atom-state", 0x0F},    // 0b00001111
    {"bond", 0x13},          // 0b00010011
    {"bond-state", 0x33},    // 0b00110011
    {NULL, 0}};

// The following defines the static SettingInfo table
#define SETTINGINFO_IMPLEMENTATION
#include "pymol/settings/setting_info.h"

template <>
const char* _SettingGet<const char*>(int index, const CSetting* I);

// get level name for setting index (for feedback)
const SettingInfoRec& SettingGetInfo(int index) { return SettingInfo[index]; }

static SettingColorFns s_color{};
static PyObject* (*s_color_to_python)(PyMOLGlobals* G, int color) = nullptr;

void SettingSetColorFns(const SettingColorFns& fns) { s_color = fns; }

static SettingUniqueIDFns s_unique_id{};

void SettingSetUniqueIDFns(const SettingUniqueIDFns& fns) { s_unique_id = fns; }

void SettingSetColorToPython(PyObject* (*fn)(PyMOLGlobals* G, int color)) {
  s_color_to_python = fn;
}

const char* SettingLevelGetName(unsigned index) {
  return SettingLevelInfo[SettingInfo[index].level].name;
}

// check if setting index is valid in given level-mask
bool SettingLevelCheckMask(PyMOLGlobals* G, int index, unsigned char mask) {
  unsigned char validmask = SettingLevelInfo[SettingInfo[index].level].mask;
  return (0 == (mask & ~validmask));
}

// check if setting index is valid in given level
bool SettingLevelCheck(PyMOLGlobals* G, int index, unsigned char level) {
  return SettingLevelCheckMask(G, index, SettingLevelInfo[level].mask);
}
static void SettingRecCopy(int const index, SettingRec const& src,
                           SettingRec& dst) {
  switch (SettingInfo[index].type) {
    case cSetting_string:
      dst.set_s(src.str_ ? src.str_->c_str() : nullptr);
      break;
    case cSetting_float3:
      dst.set_3f(src.float3_);
      break;
    default:
      dst.set_i(src.int_);
      break;
  }
  dst.defined = src.defined;
}

CSetting* SettingCopyAll(PyMOLGlobals* G, const CSetting* src, CSetting* dst) {
  if (!src) {
    delete dst;
    return nullptr;
  }

  if (!dst) {
    dst = SettingNew(G);
  }

  for (int index = 0; index < cSetting_INIT; ++index) {
    SettingRecCopy(index, src->info[index], dst->info[index]);
  }

  return dst;
}

void SettingStoreDefault(PyMOLGlobals* G) {
  G->Default = SettingCopyAll(G, G->Setting, G->Default);
}

void SettingPurgeDefault(PyMOLGlobals* G) { DeleteP(G->Default); }

void SettingUniqueDetachChain(PyMOLGlobals* G, int unique_id) {
  CSettingUnique* I = G->SettingUnique;
  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt == I->id2offset.end()) {
    return;
  }
  int offset = offsetIt->second;
  I->id2offset.erase(offsetIt);

  while (offset) {
    auto entry = &I->entry[offset];
    auto next = entry->next;
    entry->next = I->next_free;
    I->next_free = offset;
    offset = next;
  }
}

static void SettingUniqueExpand(PyMOLGlobals* G) {
  CSettingUnique* I = G->SettingUnique;

  if (!I->next_free) {
    int new_n_alloc = (I->entry.size() * 3) / 2;
    auto sizeBefore = I->entry.size();
    I->entry.resize(sizeBefore + new_n_alloc, SettingUniqueEntry{});
    for (int a = sizeBefore; a < new_n_alloc; a++) {
      I->entry[a].next = I->next_free;
      I->next_free = a;
    }
  }
}

static SettingUniqueEntry* SettingFindSettingUniqueEntry(PyMOLGlobals* G,
                                                         int unique_id,
                                                         int setting_id) {
  CSettingUnique* I = G->SettingUnique;
  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt != I->id2offset.end()) {
    SettingUniqueEntry* entry;
    for (int offset = offsetIt->second; offset; offset = entry->next) {
      entry = &I->entry[offset];
      if (entry->setting_id == setting_id) {
        return entry;
      }
    }
  }
  return nullptr;
}

int SettingUniqueCheck(PyMOLGlobals* G, int unique_id, int setting_id) {
  return SettingFindSettingUniqueEntry(G, unique_id, setting_id) != nullptr;
}

/**
 * Return true for convertible types and set int-compatible types to int
 */
inline bool type_upcast(int& type) {
  switch (type) {
    case cSetting_int:
    case cSetting_color:
    case cSetting_boolean:
      type = cSetting_int;
    case cSetting_float:
      return true;
  }
  return false;
}

bool SettingUniqueGetTypedValuePtr(PyMOLGlobals* G, int unique_id,
                                   int setting_id, int setting_type,
                                   void* value) {
  auto entry = SettingFindSettingUniqueEntry(G, unique_id, setting_id);
  if (!entry) return false;

  int type_from = SettingInfo[setting_id].type;

  if (type_from != setting_type) {
    if (!type_upcast(type_from) || !type_upcast(setting_type)) {
      PRINTFB(G, FB_Setting, FB_Errors)
      " Setting-Error: type mismatch\n" ENDFB(G);
      return false;
    }
  }

  if (setting_type == cSetting_float3) {
    *(const float**)value = entry->value.float3_;
  } else if (setting_type == type_from) {
    *(int*)value = entry->value.int_;
  } else if (setting_type == cSetting_int) {
    *(int*)value = (int)entry->value.float_;
  } else {  // setting_type == cSetting_float
    *(float*)value = (float)entry->value.int_;
  }

  return true;
}

/**
 * Warning: Returns colors as (fff) tuple instead of color index
 *
 * @pre GIL
 */
PyObject* SettingUniqueGetPyObject(PyMOLGlobals* G, int unique_id, int index) {
  assert(PyGILState_Check());

  int type = SettingGetType(G, index);

  union {
    int val_i;
    float val_f;
    const float* ptr_3f;
  };

  if (SettingUniqueGetTypedValuePtr(G, unique_id, index, type, &ptr_3f)) {
    switch (type) {
      case cSetting_boolean:
        return CPythonVal_New_Boolean(val_i);
      case cSetting_int:
        return CPythonVal_New_Integer(val_i);
      case cSetting_float:
        return CPythonVal_New_Float(val_f);
      case cSetting_color:
        if (s_color_to_python) return s_color_to_python(G, val_i);
        return CPythonVal_New_Integer(val_i);
      case cSetting_float3: {
        PyObject* result = PyTuple_New(3);
        PyTuple_SET_ITEM(result, 0, CPythonVal_New_Float(ptr_3f[0]));
        PyTuple_SET_ITEM(result, 1, CPythonVal_New_Float(ptr_3f[1]));
        PyTuple_SET_ITEM(result, 2, CPythonVal_New_Float(ptr_3f[2]));
        return result;
      }
    }
  }

  return nullptr;
}

static int SettingUniqueEntry_IsSame(SettingUniqueEntry* entry,
                                     int setting_type, const void* value) {
  if (SettingInfo[entry->setting_id].type != setting_type) {
    return 0;
  }
  if (setting_type == cSetting_float3) {
    const float *v = (float*)value, *ev = entry->value.float3_;
    return v[0] == ev[0] && v[1] == ev[1] && v[2] == ev[2];
  } else {
    return (entry->value.int_ == *(int*)value);
  }
}

static void SettingUniqueEntry_Set(SettingUniqueEntry* entry, int value_type,
                                   const void* value) {
  int setting_type = SettingGetType(entry->setting_id);
  switch (value_type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      if (setting_type == cSetting_float) {
        entry->value.float_ = *(int*)value;
      } else {
        entry->value.int_ = *(int*)value;
      }
      break;
    case cSetting_float:
      if (setting_type != cSetting_float) {
        entry->value.int_ = *(float*)value;
      } else {
        entry->value.float_ = *(float*)value;
      }
      break;
    case cSetting_float3:
      memcpy(entry->value.float3_, *(const float**)value, sizeof(float) * 3);
      break;
    default:
      printf("SettingUniqueEntry_Set-Error: unsupported type %d\n", value_type);
  }
}

/**
 * Return false if setting was not set (nothing changed)
 */
bool SettingUniqueUnset(PyMOLGlobals* G, int unique_id, int setting_id) {
  auto I = G->SettingUnique;
  auto offsetIt = I->id2offset.find(unique_id);

  if (offsetIt != I->id2offset.end()) {
    for (int prev = 0, offset = offsetIt->second; offset;
         prev = offset, offset = I->entry[offset].next) {
      if (I->entry[offset].setting_id != setting_id) continue;

      if (!prev) { /* if first entry in list */
        I->id2offset.erase(offsetIt);
        if (I->entry[offset].next) { /* set new list start */
          I->id2offset[unique_id] = I->entry[offset].next;
        }
      } else { /* otherwise excise from middle or end */
        I->entry[prev].next = I->entry[offset].next;
      }
      I->entry[offset].next = I->next_free;
      I->next_free = offset;

      return true;
    }
  }
  return false;
}

int SettingUniqueSetTypedValue(PyMOLGlobals* G, int unique_id, int setting_id,
                               int setting_type, const void* value)

/* set value to nullptr in order to delete setting */
{
  CSettingUnique* I = G->SettingUnique;
  int isset = false;

  if (!value) {
    return SettingUniqueUnset(G, unique_id, setting_id);
  }

  auto offsetIt = I->id2offset.find(unique_id);
  if (offsetIt != I->id2offset.end()) { /* setting list exists for atom */
    int offset = offsetIt->second;
    int prev = 0;
    int found = false;
    while (offset) {
      auto entry = &I->entry[offset];
      if (entry->setting_id == setting_id) {
        found = true; /* this setting is already defined */
        if (!SettingUniqueEntry_IsSame(entry, setting_type, value)) {
          SettingUniqueEntry_Set(entry, setting_type, value);
          isset = true;
        }
        break;
      }
      prev = offset;
      offset = entry->next;
    }
    if ((!found) &&
        value) { /* setting not found in existing list, so append new value */
      if (!I->next_free) SettingUniqueExpand(G);
      if (I->next_free) {
        offset = I->next_free;
        {
          auto entry = &I->entry[offset];
          I->next_free = entry->next;
          entry->next = 0;

          if (prev) { /* append onto existing list */
            I->entry[prev].next = offset;
            entry->setting_id = setting_id;
            SettingUniqueEntry_Set(entry, setting_type, value);
            isset = true;
          } else {
            I->id2offset[unique_id] = offset;
            /* create new list */
            entry->setting_id = setting_id;
            SettingUniqueEntry_Set(entry, setting_type, value);
            isset = true;
          }
        }
      }
    }
  } else if (value &&
             offsetIt == I->id2offset.end()) { /* new setting list for atom */
    if (!I->next_free) SettingUniqueExpand(G);
    if (I->next_free) {
      int offset = I->next_free;
      auto entry = &I->entry[offset];

      I->id2offset[unique_id] = offset;
      I->next_free = entry->next;
      entry->setting_id = setting_id;
      entry->next = 0;
      SettingUniqueEntry_Set(entry, setting_type, value);
      isset = true;
    }
  } else {
    /* unhandled error */
  }
  return isset;
}

/**
 * @pre GIL
 */
bool SettingUniqueSetPyObject(PyMOLGlobals* G, int unique_id, int index,
                              PyObject* value) {
  assert(PyGILState_Check());

  if (!value) return SettingUniqueUnset(G, unique_id, index);

  int type = SettingGetType(G, index);

  float val_3f[3];
  union {
    int val_i;
    float val_f;
    float* ptr_3f;
  };

  switch (type) {
    case cSetting_boolean:
    case cSetting_int:
      ok_assert(1, PConvPyObjectToInt(value, &val_i));
      break;
    case cSetting_float:
      ok_assert(1, PConvPyObjectToFloat(value, &val_f));
      break;
    case cSetting_color:
      if (!PConvPyIntToInt(value, &val_i)) {
        OrthoLineType sval;
        ok_assert(1, PConvPyStrToStr(value, sval, OrthoLineLength));
        val_i = s_color.getIndex(G, sval);
      }
      break;
    case cSetting_float3:
      if (!PConvPyListOrTupleToFloatArrayInPlace(value, val_3f, 3)) {
        OrthoLineType sval;
        ok_assert(1, PConvPyStrToStr(value, sval, OrthoLineLength) &&
                         sscanf(sval, "%f%f%f", &val_3f[0], &val_3f[1],
                                &val_3f[2]) == 3);
      }
      ptr_3f = val_3f;
      break;
    default:
      PRINTFB(G, FB_Python, FB_Errors)
      " Python-Error: atom-state-level setting unsupported type=%d\n",
          type ENDFB(G);
      return false;
  }

  return SettingUniqueSetTypedValue(G, unique_id, index, type, &val_i);

ok_except1:
  PRINTFB(G, FB_Setting, FB_Errors)
  " Setting-Error: type mismatch\n" ENDFB(G);
  return false;
}

void SettingUniqueResetAll(PyMOLGlobals* G) {
  CSettingUnique* I = G->SettingUnique;
  I->id2offset.clear();
  I->entry.clear();
  I->entry.resize(CSettingUnique::numInitEntries, SettingUniqueEntry{});
  /* note: intentially skip index 0  */
  for (int a = 2; a < I->entry.size(); a++) {
    I->entry[a].next = a - 1; /* 1-based linked list with 0 as sentinel */
  }
  I->next_free = I->entry.size() - 1;
}

int SettingUniqueCopyAll(PyMOLGlobals* G, int src_unique_id,
                         int dst_unique_id) {
  int ok = true;
  CSettingUnique* I = G->SettingUnique;

  auto dstOffsetIt = I->id2offset.find(dst_unique_id);
  if (dstOffsetIt != I->id2offset.end()) { /* setting list exists for atom */
    PRINTFB(G, FB_Setting, FB_Errors)
    " SettingUniqueCopyAll-Bug: merging settings not implemented\n" ENDFB(G);
    ok = false;
  } else if (dstOffsetIt ==
             I->id2offset.end()) { /* new setting list for atom */
    auto srcOffsetIt = I->id2offset.find(src_unique_id);
    if (srcOffsetIt != I->id2offset.end()) {
      int dst_offset = 0;
      for (int src_offset = srcOffsetIt->second; src_offset;
           src_offset = I->entry[src_offset].next) {
        SettingUniqueExpand(G);  // this may reallocate I->entry

        if (!dst_offset) {
          I->id2offset[dst_unique_id] = I->next_free;
        } else {
          I->entry[dst_offset].next = I->next_free;
        }

        dst_offset = I->next_free;
        I->next_free = I->entry[dst_offset].next;
        I->entry[dst_offset] = I->entry[src_offset];
        I->entry[dst_offset].next = 0;
      }
    }
  } else {
    ok = false;
    /* unhandled error */
  }

  return ok;
}

static void SettingUniqueInit(PyMOLGlobals* G) {
  G->SettingUnique = new CSettingUnique();
  SettingUniqueResetAll(G);
}

static void SettingUniqueFree(PyMOLGlobals* G) { DeleteP(G->SettingUnique); }

/**
 * For unique_id remapping during partial session loading
 */
int SettingUniqueConvertOldSessionID(PyMOLGlobals* G, int old_unique_id) {
  CSettingUnique* I = G->SettingUnique;
  int unique_id = old_unique_id;
  if (I->old2new) {
    auto oldIdIt = I->old2new->find(old_unique_id);
    if (oldIdIt != I->old2new->end()) {
      unique_id = oldIdIt->second;
    } else {
      unique_id = s_unique_id.getNew(G);
      I->old2new->emplace(old_unique_id, unique_id);
    }
  } else {
    s_unique_id.reserve(G, unique_id);
  }
  return unique_id;
}

/**
 * Return true if the given setting index should not be stored to PSE.
 *
 * Blacklisted are unused and system-dependent settings.
 */
static bool is_session_blacklisted(int index) {
  if (index >= cSetting_INIT ||
      SettingInfo[index].level == cSettingLevel_unused) {
    return true;
  }

  switch (index) {
    case cSetting_antialias_shader:
    case cSetting_ati_bugs:
    case cSetting_cache_max:
    case cSetting_cgo_shader_ub_color:
    case cSetting_cgo_shader_ub_flags:
    case cSetting_cgo_shader_ub_normal:
    case cSetting_colored_feedback:
    case cSetting_cylinder_shader_ff_workaround:
    case cSetting_defer_updates:
    case cSetting_fast_idle:
    case cSetting_fetch_path:
    case cSetting_internal_feedback:
    case cSetting_internal_gui:
    case cSetting_internal_prompt:
    case cSetting_logging:
    case cSetting_max_threads:
    case cSetting_mouse_grid:
    case cSetting_mouse_scale:
    case cSetting_nb_spheres_use_shader:
    case cSetting_no_idle:
    case cSetting_nvidia_bugs:
    case cSetting_presentation:
    case cSetting_precomputed_lighting:
    case cSetting_render_as_cylinders:
    case cSetting_security:
    case cSetting_session_changed:
    case cSetting_session_file:
    case cSetting_session_migration:
    case cSetting_session_version_check:
    case cSetting_shaders_from_disk:
    case cSetting_show_progress:
    case cSetting_slow_idle:
    case cSetting_stereo:
    case cSetting_stereo_double_pump_mono:
    case cSetting_stereo_mode:
    case cSetting_suspend_deferred:
    case cSetting_suspend_undo:
    case cSetting_suspend_undo_atom_count:
    case cSetting_suspend_updates:
    case cSetting_text:
    case cSetting_trilines:
    case cSetting_use_geometry_shaders:
    case cSetting_use_shaders:
    case cSetting_pick32bit:
    case cSetting_display_scale_factor:
      return true;
  }

  return false;
}

/**
 * @pre GIL
 */
int SettingUniqueFromPyList(PyMOLGlobals* G, PyObject* list,
                            int partial_restore) {
  assert(PyGILState_Check());

  int ok = true;
  if (!partial_restore) {
    SettingUniqueResetAll(G);
  }
  if (list)
    if (PyList_Check(list)) {
      ov_size n_id = PyList_Size(list);
      ov_size a;
      for (a = 0; a < n_id; a++) {
        PyObject* id_list = PyList_GetItem(list, a);
        int unique_id;
        if (ok) ok = PyList_Check(id_list);
        if (ok) ok = (PyList_Size(id_list) > 1);
        if (ok) ok = PConvPyIntToInt(PyList_GetItem(id_list, 0), &unique_id);
        if (ok && partial_restore) {
          unique_id = SettingUniqueConvertOldSessionID(G, unique_id);
        }
        if (ok) {
          ov_size n_set = 0;

          PyObject* setting_list = PyList_GetItem(id_list, 1);
          if (ok) ok = PyList_Check(setting_list);
          if (ok) n_set = PyList_Size(setting_list);
          if (ok) {
            ov_size b;
            for (b = 0; b < n_set; b++) {
              PyObject* entry_list = PyList_GetItem(setting_list, b);
              if (ok) ok = PyList_Check(entry_list);
              if (ok) ok = (PyList_Size(entry_list) > 2);
              if (ok) {
                int setting_id;
                int setting_type;
                float val_3f[3];
                union {
                  int int_;
                  float float_;
                  float* float3_;
                } value_store;
                if (ok)
                  ok = PConvPyIntToInt(PyList_GetItem(entry_list, 0),
                                       &setting_id);
                if (ok)
                  ok = PConvPyIntToInt(PyList_GetItem(entry_list, 1),
                                       &setting_type);
                if (ok) switch (setting_type) {
                    case cSetting_int:
                    case cSetting_color:
                    case cSetting_boolean:
                      ok = PConvPyIntToInt(PyList_GetItem(entry_list, 2),
                                           &value_store.int_);
                      if (setting_type == cSetting_color) {
                        value_store.int_ =
                            s_color.convertOldSessionIndex(G, value_store.int_);
                      }
                      break;
                    case cSetting_float:
                      ok = PConvPyFloatToFloat(PyList_GetItem(entry_list, 2),
                                               &value_store.float_);
                      break;
                    case cSetting_float3: {
                      CPythonVal* el =
                          CPythonVal_PyList_GetItem(G, entry_list, 2);
                      value_store.float3_ = val_3f;
                      ok = PConvPyListToFloatArrayInPlace(
                          el, value_store.float3_, 3);
                      CPythonVal_Free(el);
                    } break;
                  }
                if (ok) {
                  SettingUniqueSetTypedValue(G, unique_id, setting_id,
                                             setting_type, &value_store);
                }
              }
            }
          }
        }
      }
    }
  return ok;
}

/**
 * @pre GIL
 */
PyObject* SettingUniqueAsPyList(PyMOLGlobals* G) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  CSettingUnique* I = G->SettingUnique;
  {
    auto n_entry = I->id2offset.size();
    result = PyList_New(n_entry);
    if (result) {
      n_entry = 0;
      for (auto id2offset : I->id2offset) {
        PyObject* setting_list = nullptr;
        auto unique_id = id2offset.first;
        auto offset = id2offset.second;
        {
          int n_set = 0;

          /* count number of settings for this unique_id */

          auto save_offset = offset;
          while (offset) {
            auto entry = &I->entry[offset];
            n_set++;
            offset = entry->next;
          }

          /* create and insert list for each setting */

          setting_list = PyList_New(n_set);
          n_set = 0;
          offset = save_offset;
          while (offset) {
            PyObject* setting_entry = PyList_New(3);
            auto entry = &I->entry[offset];
            int type = SettingInfo[entry->setting_id].type;
            PyList_SetItem(setting_entry, 0, PyInt_FromLong(entry->setting_id));
            PyList_SetItem(setting_entry, 1, PyInt_FromLong(type));
            switch (type) {
              case cSetting_int:
              case cSetting_color:
              case cSetting_boolean:
                PyList_SetItem(setting_entry, 2,
                               PyInt_FromLong(entry->value.int_));
                break;
              case cSetting_float:
                PyList_SetItem(
                    setting_entry, 2,
                    PyFloat_FromDouble(*(float*)&entry->value.float_));
                break;
              case cSetting_float3:
                PyList_SetItem(
                    setting_entry, 2,
                    PConvFloatArrayToPyList((float*)&entry->value.float3_, 3));
                break;
            }
            PyList_SetItem(setting_list, n_set, setting_entry);
            n_set++;
            offset = entry->next;
          }
        }

        /* add this unique_id set into the overall list */

        {
          PyObject* unique_list = PyList_New(2);
          PyList_SetItem(unique_list, 0, PyInt_FromLong(unique_id));
          PyList_SetItem(unique_list, 1, setting_list);
          PyList_SetItem(result, n_entry, unique_list);
        }
        n_entry++;
      }
    }
  }
  return (PConvAutoNone(result));
}

int SettingSetSmart_i(PyMOLGlobals* G, CSetting* set1, CSetting* set2,
                      int index, int value) {
  int dummy;
  if (set1 && SettingGetIfDefined_i(G, set1, index, &dummy)) {
    return SettingSet_i(set1, index, value);
  }
  if (set2 && SettingGetIfDefined_i(G, set2, index, &dummy)) {
    return SettingSet_i(set2, index, value);
  }
  return SettingSetGlobal_i(G, index, value);
}

/**
 * @pre GIL
 */
int SettingSetGlobalsFromPyList(PyMOLGlobals* G, PyObject* list) {
  assert(PyGILState_Check());

  int ok = true;

  CSetting* I = G->Setting;

  if (list)
    if (PyList_Check(list)) ok = SettingFromPyList(I, list);

  /* restore the following settings  */

  if (G->Option->no_quit) {
    SettingSet_b(I, cSetting_presentation_auto_quit, 0);
  }

  return (ok);
}

/**
 * @pre GIL
 */
PyObject* SettingGetGlobalsAsPyList(PyMOLGlobals* G) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  CSetting* I = G->Setting;
  result = SettingAsPyList(I);
  return (PConvAutoNone(result));
}

/**
 * @pre GIL
 */
static PyObject* get_list(CSetting* I, int index, bool incl_blacklisted) {
  assert(PyGILState_Check());

  PyObject *result = nullptr, *value = nullptr;
  int setting_type = SettingInfo[index].type;

  if (!incl_blacklisted && is_session_blacklisted(index)) {
    return nullptr;
  }

  switch (setting_type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      value = PyInt_FromLong(I->info[index].int_);
      break;
    case cSetting_float:
      value = PyFloat_FromDouble(I->info[index].float_);
      break;
    case cSetting_float3:
      value = PConvFloatArrayToPyList(I->info[index].float3_, 3);
      break;
    case cSetting_string:
      value = PyString_FromString(SettingGet<const char*>(index, I));
      break;
  }

  if (value) {
    result = PyList_New(3);
    PyList_SetItem(result, 0, PyInt_FromLong(index));
    PyList_SetItem(result, 1, PyInt_FromLong(setting_type));
    PyList_SetItem(result, 2, value);
  }

  return result;
}

/**
 * @pre GIL
 */
PyObject* SettingAsPyList(CSetting* I, bool incl_blacklisted) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  int a;

  if (I) {
    std::vector<PyObject*> list;
    list.reserve(cSetting_INIT);

    for (a = 0; a < cSetting_INIT; a++) {
      if (I->info[a].defined) {
        PyObject* item = get_list(I, a, incl_blacklisted);
        if (item != nullptr) {
          list.push_back(item);
        }
      }
    }

    result = PConvToPyObject(list);
  }
  return (PConvAutoNone(result));
}
/**
 * @pre GIL
 */
static int set_list(CSetting* I, PyObject* list) {
  assert(PyGILState_Check());

  int index = -1;
  int setting_type = -1;

  union {
    int val_i;
    float val_f;
    float val_3f[3];
    const char* val_s;
  };

  if (list == nullptr || CPythonVal_IsNone(list)) return true;

  ok_assert(1, PyList_Check(list));
  ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(I->G, list, 0, &index));
  ok_assert(1,
            CPythonVal_PConvPyIntToInt_From_List(I->G, list, 1, &setting_type));

  if (is_session_blacklisted(index)) return true;

  switch (setting_type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(I->G, list, 2, &val_i));
      if (setting_type == cSetting_color)
        val_i = s_color.convertOldSessionIndex(I->G, val_i);
      SettingSet_i(I, index, val_i);
      break;
    case cSetting_float:
      ok_assert(
          1, CPythonVal_PConvPyFloatToFloat_From_List(I->G, list, 2, &val_f));
      SettingSet_f(I, index, val_f);
      break;
    case cSetting_float3:
      ok_assert(1, CPythonVal_PConvPyListToFloatArrayInPlaceAutoZero_From_List(
                       I->G, list, 2, val_3f, 3));
      SettingSet_3fv(I, index, val_3f);
      break;
    case cSetting_string:
      ok_assert(1, val_s = PyString_AsString(PyList_GetItem(list, 2)));
      SettingSet_s(I, index, val_s);
      break;
    default:
      ok_raise(1);
  }

  return true;
ok_except1:
  printf(" set_list-Error: i=%d, t=%d\n", index, setting_type);
  return false;
}
/**
 * Used to set object and object-state level settings from PSEs
 *
 * @pre GIL
 */
CSetting* SettingNewFromPyList(PyMOLGlobals* G, PyObject* list) {
  assert(PyGILState_Check());

  int ok = true;
  ov_size size;
  ov_size a;
  CSetting* I = nullptr;
  if (ok) ok = (list != nullptr);
  if (ok) ok = PyList_Check(list);
  if (ok) {
    I = SettingNew(G);
    size = PyList_Size(list);
    for (a = 0; a < size; a++) {
      if (ok) ok = set_list(I, PyList_GetItem(list, a));
    }
  }
  return (I);
}
/**
 * @pre GIL
 */
int SettingFromPyList(CSetting* I, PyObject* list) {
  assert(PyGILState_Check());

  int ok = true;
  ov_size size;
  ov_size a;

  if (ok) ok = (I != nullptr);
  if (ok) ok = PyList_Check(list);
  if (ok) {
    size = PyList_Size(list);
    for (a = 0; a < size; a++) {
      if (!set_list(I, PyList_GetItem(list, a))) ok = false;
    }
  }
  return (ok);
}
void SettingCheckHandle(PyMOLGlobals* G,
                        pymol::copyable_ptr<CSetting>& handle) {
  if (!handle) handle.reset(SettingNew(G));
}
int SettingGetTextValue(PyMOLGlobals* G, const CSetting* set1,
                        const CSetting* set2, int index, char* buffer) {
  const char* sptr = SettingGetTextPtr(G, set1, set2, index, buffer);
  if (!sptr) return 0;

  if (sptr != buffer) {
    if (strlen(sptr) > OrthoLineLength) {
      PRINTFB(G, FB_Setting, FB_Warnings)
      "Setting-Warning: text longer than OrthoLineLength" ENDFB(G);
    }

    strncpy(buffer, sptr, OrthoLineLength);
  }

  return 1;
}
/**
 * Returns a pointer to the internal string representation if available,
 * or it formats the value into buffer and returns a pointer to buffer.
 */
const char* SettingGetTextPtr(PyMOLGlobals* G, const CSetting* set1,
                              const CSetting* set2, int index, char* buffer) {
  int type;
  const char* sptr = nullptr;
  const float* ptr;
  type = SettingGetType(G, index);
  switch (type) {
    case cSetting_boolean:
      sprintf(buffer, SettingGet_b(G, set1, set2, index) ? "on" : "off");
      break;
    case cSetting_int:
      sprintf(buffer, "%d", SettingGet_i(G, set1, set2, index));
      break;
    case cSetting_float:
      sprintf(buffer, "%1.5f", SettingGet_f(G, set1, set2, index));
      break;
    case cSetting_float3:
      ptr = SettingGet_3fv(G, set1, set2, index);
      sprintf(buffer, "[ %1.5f, %1.5f, %1.5f ]", ptr[0], ptr[1], ptr[2]);
      break;
    case cSetting_color: {
      int color = SettingGet_color(G, set1, set2, index);
      switch (color) {
        case cColorAtomic:
          strcpy(buffer, "atomic");
          break;
        case cColorObject:
          strcpy(buffer, "object");
          break;
        case cColorFront:
          strcpy(buffer, "front");
          break;
        case cColorBack:
          strcpy(buffer, "back");
          break;
        case cColorDefault:
          strcpy(buffer, "default");
          break;
        default:
          sptr = s_color.getName(G, color);
          if (sptr) return sptr;
          strcpy(buffer, "invalid");
      }
    } break;
    case cSetting_string:
      return SettingGet_s(G, set1, set2, index);
    default:
      return nullptr;
  }
  return buffer;
}
/**
 * @pre GIL
 */
int SettingSetFromTuple(PyMOLGlobals* G, CSetting* I, int index,
                        PyObject* tuple) {
  assert(PyGILState_Check());

  PyObject* value;
  int type;
  int ok = true;
  if (!I) I = G->Setting; /* fall back on global settings */

  /* this data structure has been pre-checked at the python level... */

  type = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
  value = PyTuple_GetItem(tuple, 1);
  switch (type) {
    case cSetting_boolean:
    case cSetting_int:
      SettingSet_i(I, index, PyInt_AsLong(value));
      break;
    case cSetting_float:
      SettingSet_f(I, index, (float)PyFloat_AsDouble(value));
      break;
    case cSetting_float3:
      float tmp[3];
      PyArg_ParseTuple(value, "fff", tmp, tmp + 1, tmp + 2);
      SettingSet_3fv(I, index, tmp);
      break;
    case cSetting_color:
      SettingSet_color(I, index, PyString_AsString(value));
      break;
    case cSetting_string:
      SettingSet_s(I, index, PyString_AsString(value));
      break;
    default:
      ok = false;
      break;
  }
  return (ok);
}
int SettingStringToTypedValue(PyMOLGlobals* G, int index, const char* st,
                              int* type, int* value) {
  int ok = true;
  int newvalue;
  float newfvalue;
  /* this data structure has been pre-checked at the python level... */

  *type = SettingGetType(G, index);

  switch (*type) {
    case cSetting_boolean:
      if ((!*st) || (*st == '0') || (*st == 'F') ||
          WordMatchExact(G, st, "on", true) ||
          WordMatchExact(G, st, "false", true)) {
        newvalue = 0;
      } else {
        newvalue = 1;
      }
      if (newvalue != *value) {
        *value = newvalue;
      }
      break;
    case cSetting_int:
      if (sscanf(st, "%d", &newvalue) != 1) {
        ok = false;
      } else if (newvalue != *value) {
        *value = newvalue;
      }
      break;
    case cSetting_float:
      if (sscanf(st, "%f", &newfvalue) != 1) {
        ok = false;
      } else if (newfvalue != *((float*)value)) {
        *(float*)value = newfvalue;
      }
      break;
    case cSetting_color: {
      int color_index = s_color.getIndex(G, st);
      if (*(value) != color_index) {
        *(value) = color_index;
      }
    } break;
    default:
      ok = false;
      break;
  }
  return (ok);
}

int SettingSetFromString(PyMOLGlobals* G, CSetting* I, int index,
                         const char* st) {
  int type;
  int ok = true;
  if (!I) I = G->Setting; /* fall back on global settings */

  /* this data structure has been pre-checked at the python level... */

  type = SettingGetType(G, index);

  switch (type) {
    case cSetting_boolean:
      if ((!*st) || (*st == '0') || (*st == 'F') ||
          WordMatchExact(G, st, "on", true) ||
          WordMatchExact(G, st, "false", true))
        SettingSet_b(I, index, 0);
      else
        SettingSet_b(I, index, 1);
      break;
    case cSetting_int: {
      int tmp;
      if (sscanf(st, "%d", &tmp) == 1)
        SettingSet_i(I, index, tmp);
      else
        ok = false;
    } break;
    case cSetting_float: {
      float tmp;
      if (sscanf(st, "%f", &tmp) == 1)
        SettingSet_f(I, index, tmp);
      else
        ok = false;
    } break;
    case cSetting_float3: {
      float tmp[3];
      if (sscanf(st, "%f%f%f", tmp, tmp + 1, tmp + 2) == 3)
        SettingSet_3fv(I, index, tmp);
      else
        ok = false;
    } break;
    case cSetting_color:
      SettingSet_color(I, index, st);
      break;
    case cSetting_string:
      SettingSet_s(I, index, st);
      break;
    default:
      ok = false;
      break;
  }
  return (ok);
}
/**
 * Warning: Returns colors as (fff) tuple instead of color index
 *
 * @pre GIL
 */
PyObject* SettingGetPyObject(PyMOLGlobals* G, const CSetting* set1,
                             const CSetting* set2, int index) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  const float* ptr;
  int type = SettingGetType(G, index);

  switch (type) {
    case cSetting_boolean:
      result = CPythonVal_New_Boolean(SettingGet_b(G, set1, set2, index));
      break;
    case cSetting_int:
      result = CPythonVal_New_Integer(SettingGet_i(G, set1, set2, index));
      break;
    case cSetting_float:
      result = CPythonVal_New_Float(SettingGet_f(G, set1, set2, index));
      break;
    case cSetting_float3:
      ptr = SettingGet_3fv(G, set1, set2, index);
      result =
          Py_BuildValue("(fff)", pymol::pretty_f2d(ptr[0]),
                        pymol::pretty_f2d(ptr[1]), pymol::pretty_f2d(ptr[2]));
      break;
    case cSetting_color: {
      int retcol = SettingGet_color(G, set1, set2, index);
      if (retcol > 0) {
        const float* col;
        col = s_color.get(G, retcol);
        result = Py_BuildValue("(fff)", col[0], col[1], col[2]);
      }
    } break;
    case cSetting_string:
      result = PyString_FromString(SettingGet_s(G, set1, set2, index));
      break;
  }
  return result;
}

/**
 * @pre GIL
 */
PyObject* SettingGetTuple(PyMOLGlobals* G, const CSetting* set1,
                          const CSetting* set2, int index) {
  assert(PyGILState_Check());

  PyObject* result = nullptr;
  const float* ptr;
  int type = SettingGetType(G, index);

  switch (type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      result = Py_BuildValue("ii", type, SettingGet_i(G, set1, set2, index));
      break;
    case cSetting_float:
      result = Py_BuildValue(
          "if", type, pymol::pretty_f2d(SettingGet_f(G, set1, set2, index)));
      break;
    case cSetting_float3:
      ptr = SettingGet_3fv(G, set1, set2, index);
      result =
          Py_BuildValue("i(fff)", type, pymol::pretty_f2d(ptr[0]),
                        pymol::pretty_f2d(ptr[1]), pymol::pretty_f2d(ptr[2]));
      break;
    case cSetting_string:
      result = Py_BuildValue("is", type, SettingGet_s(G, set1, set2, index));
      break;
    default:
      result = PConvAutoNone(Py_None);
      break;
  }
  return result;
}
CSetting* SettingNew(PyMOLGlobals* G) { return new CSetting(G); }

CSetting::CSetting(const CSetting& other) { *this = other; }

CSetting& CSetting::operator=(const CSetting& other) {
  for (int index = 0; index < cSetting_INIT; ++index) {
    SettingRecCopy(index, other.info[index], this->info[index]);
  }
  return *this;
}
CSetting::~CSetting() {
  // need to free strings
  for (int index = 0; index < cSetting_INIT; ++index) {
    if (SettingInfo[index].type == cSetting_string) {
      info[index].delete_s();
    }
  }
}
CSetting::CSetting(PyMOLGlobals* G) : G(G) {}
/**
 * Restore the default value from `src` or `SettingInfo`
 */
void SettingRestoreDefault(CSetting* I, int index, const CSetting* src) {
  // 1) from stored default if provided
  if (src) {
    SettingRecCopy(index, src->info[index], I->info[index]);
    return;
  }

  // 2) from SettingInfo
  auto& rec = SettingInfo[index];

  switch (rec.type) {
    case cSetting_blank:
      break;
    case cSetting_boolean:
    case cSetting_int:
      I->info[index].set_i(rec.value.i[0]);
      break;
    case cSetting_float:
      I->info[index].set_f(rec.value.f[0]);
      break;
    case cSetting_float3:
      I->info[index].set_3f(rec.value.f);
      break;
    case cSetting_string:
      I->info[index].delete_s();
      break;
    case cSetting_color:
      SettingSet_color(I, index, rec.value.s);
      break;
    default:
      // coding error
      printf(" ERROR: unkown type\n");
  };

  I->info[index].defined = false;
}
int SettingUnset(CSetting* I, int index) {
  if (I) {
    SettingRec& sr = I->info[index];
    if (!sr.defined) {
      return false;
    }
    sr.defined = false;
    sr.changed = true;
  }
  return true;
}
int SettingGetType(int index) { return (SettingInfo[index].type); }
template <>
int _SettingGet<int>(int index, const CSetting* I) {
  PyMOLGlobals* G = I->G;
  int result;
  switch (SettingInfo[index].type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_color:
      result = I->info[index].int_;
      break;
    case cSetting_float:
      result = (int)I->info[index].float_;
      break;
    default:
      PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (int) %d\n", index ENDFB(G);
      result = 0;
      break;
  }
  return (result);
}
template <>
bool _SettingGet<bool>(int index, const CSetting* I) {
  PyMOLGlobals* G = I->G;
  switch (SettingInfo[index].type) {
    case cSetting_boolean:
    case cSetting_int:
    case cSetting_float:
      return I->info[index].int_ != 0;
    default:
      PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (boolean) %d\n", index ENDFB(G);
      return false;
  }
}
template <>
float _SettingGet<float>(int index, const CSetting* I) {
  float result;
  PyMOLGlobals* G = I->G;
  switch (SettingInfo[index].type) {
    case cSetting_color:
      PRINTFB(G, FB_Setting, FB_Warnings)
      " Setting-Warning: type read mismatch (float/color) %d\n", index ENDFB(G);
    case cSetting_boolean:
    case cSetting_int:
      result = (float)I->info[index].int_;
      break;
    case cSetting_float:
      result = I->info[index].float_;
      break;
    default:
      PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (float) %d\n", index ENDFB(G);
      result = 0.0F;
  }
  return (result);
}
template <>
const char* _SettingGet<const char*>(int index, const CSetting* I) {
  const char* result;
  PyMOLGlobals* G = I->G;
  switch (SettingInfo[index].type) {
    case cSetting_string:
      if (I->info[index].str_) {
        result = I->info[index].str_->c_str();
      } else {
        result = SettingInfo[index].value.s;
      }
      break;
    default:
      PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type read mismatch (string) %d\n", index ENDFB(G);
      result = nullptr;
  }
  return (char*)result;
}
template <>
const float* _SettingGet<const float*>(int index, const CSetting* I) {
  if (SettingInfo[index].type != cSetting_float3) {
    PyMOLGlobals* G = I->G;
    PRINTFB(G, FB_Setting, FB_Errors)
    " Setting-Error: type read mismatch (float3) %d\n", index ENDFB(G);
    return nullptr;
  }
  return I->info[index].float3_;
}
int SettingSet_i(CSetting* I, int index, int value) {
  int ok = true;
  if (I) {
    PyMOLGlobals* G = I->G;
    {
      int setting_type = SettingInfo[index].type;
      switch (setting_type) {
        case cSetting_boolean:
        case cSetting_int:
        case cSetting_color:
          I->info[index].set_i(value);
          break;
        case cSetting_float:
          I->info[index].set_f((float)value);
          break;
        default:
          PRINTFB(G, FB_Setting, FB_Errors)
          "Setting-Error: type set mismatch (integer) %d\n", index ENDFB(G);
          ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}
static int SettingSet_color_from_3f(CSetting* I, int index,
                                    const float* vector) {
  float vals[3];
  for (int i = 0; i < 3; ++i) vals[i] = std::clamp(vector[i], 0.f, 1.f);
  return SettingSet_i(I, index, s_color.fromRGB(I->G, vals));
}

int SettingSet_color(CSetting* I, int index, const char* value) {
  int ok = true;
  int color_index;
  if (I) {
    PyMOLGlobals* G = I->G;
    color_index = s_color.getIndex(G, value);
    if ((color_index == -1) &&
        (strcmp(value, "-1") && strcmp(value, "-2") && strcmp(value, "-3") &&
         strcmp(value, "-4") && strcmp(value, "-5") &&
         strcmp(value, "default"))) {
      float vals[3];
      ok = ParseFloat3List(value, vals);
      if (ok) {
        for (auto& v : vals) v = std::clamp(v, 0.f, 1.f);
        color_index = cColor_TRGB_Bits |
                      ((int)(255 * vals[0] + 0.49999F)) << 16 |
                      ((int)(255 * vals[1] + 0.49999F)) << 8 |
                      ((int)(255 * vals[2] + 0.49999F));
      } else {
        PRINTFB(G, FB_Setting, FB_Errors)
        "Setting-Error: unknown color '%s'\n", value ENDFB(G);
      }
    }
    if (ok) {
      SettingSet_i(I, index, color_index);
    }
  }
  return (ok);
}
int SettingSet_f(CSetting* I, int index, float value) {
  int ok = true;
  if (I) {
    PyMOLGlobals* G = I->G;
    {
      int setting_type = SettingInfo[index].type;
      switch (setting_type) {
        case cSetting_boolean:
        case cSetting_int:
        case cSetting_color:
          I->info[index].set_i((int)value);
          break;
        case cSetting_float:
          I->info[index].set_f(value);
          break;
        default:
          PRINTFB(G, FB_Setting, FB_Errors)
          "Setting-Error: type set mismatch (float) %d\n", index ENDFB(G);
          ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}
int SettingSet_s(CSetting* I, int index, const char* value) {
  int ok = true;
  if (I) {
    PyMOLGlobals* G = I->G;
    {
      int setting_type = SettingInfo[index].type;
      switch (setting_type) {
        case cSetting_string:
          I->info[index].set_s(value);
          break;
        case cSetting_color:
          return SettingSet_color(I, index, value);
        default:
          PRINTFB(G, FB_Setting, FB_Errors)
          "Setting-Error: type set mismatch (string) %d\n", index ENDFB(G);
          ok = false;
      }
    }
  } else {
    ok = false;
  }
  return (ok);
}
int SettingSet_3fv(CSetting* I, int index, const float* vector) {
  switch (SettingInfo[index].type) {
    case cSetting_float3:
      I->info[index].set_3f(vector);
      return true;
    case cSetting_color:
      return SettingSet_color_from_3f(I, index, vector);
    default:
      PyMOLGlobals* G = I->G;
      PRINTFB(G, FB_Setting, FB_Errors)
      "Setting-Error: type set mismatch (float3) %d\n", index ENDFB(G);
      return false;
  }
}
int SettingGetIndex(PyMOLGlobals* G, const char* name) {
  if (auto result = get_setting_id(G->PyMOL, name)) {
    return *result;
  }
  return -1;
}
int SettingGetName(PyMOLGlobals* G, int index, SettingName name) {
  UtilNCopy(name, SettingInfo[index].name, sizeof(SettingName));
  return (name[0] != 0);
}
const char* SettingGetName(int index) { return SettingInfo[index].name; }
void SettingFreeGlobal(PyMOLGlobals* G) {
  SettingUniqueFree(G);
  DeleteP(G->Setting);
  DeleteP(G->Default);
}
void SettingInitGlobal(PyMOLGlobals* G, int alloc, int reset_gui,
                       int use_default) {
  CSetting* I = G->Setting;

  /* use function pointers to prevent the compiler from inlining every
     call in this block (a waste of RAM and time) */

  int (*set_i)(CSetting* I, int index, int value) = SettingSet_i;
  int (*set_b)(CSetting* I, int index, int value) = SettingSet_b;

  if (alloc || !I) {
    I = G->Setting = SettingNew(G);
    SettingUniqueInit(G);
  }

  if (G->Default && use_default) {
    SettingCopyAll(G, G->Default, G->Setting);

  } else {
    // copy defaults from SettingInfo table
    for (int index = 0; index < cSetting_INIT; ++index) {
      if (!reset_gui) switch (index) {
          case cSetting_internal_gui_width:
          case cSetting_internal_gui:
            continue;
        }

      SettingRestoreDefault(I, index);
    }

    // open-source has no volume_mode=1
    set_i(I, cSetting_volume_mode, 0);

    // command line arguments overwrites
    set_b(I, cSetting_auto_show_lines, G->Option->sphere_mode < 0);
    set_i(I, cSetting_auto_zoom, G->Option->zoom_mode);
    set_b(I, cSetting_auto_show_nonbonded, G->Option->sphere_mode < 0);
    set_b(I, cSetting_presentation, G->Option->presentation);
    set_i(I, cSetting_defer_builds_mode, G->Option->defer_builds_mode);
    set_b(I, cSetting_presentation_auto_quit, !G->Option->no_quit);
    set_b(I, cSetting_auto_show_spheres, G->Option->sphere_mode >= 0);
    set_i(I, cSetting_internal_feedback, G->Option->internal_feedback);

    if (G->Option->stereo_mode) {
      set_i(I, cSetting_stereo_mode, G->Option->stereo_mode);
    } else if (G->StereoCapable || G->Option->blue_line) {
      set_i(I, cSetting_stereo_mode,
            cStereo_quadbuffer); /* quadbuffer if we can */
    }

    /* In order to get electrostatic potentials in kT from the Coulomb
       equation...

       PyMOL charges: Q, in e
       PyMOL distances: r, in Angstrom
       Coulomb Constant: K = 8.987552e9 ((J*m)/(C^2))
       Angstrom Convertor: 1 A = 1e-10 m
       Coulomb Convertor: 1 e = 1.60217733e-19 C
       Angstrom Convertor: 1 A = 10e-10 m
       Dielectric Constant: D (unitless)

       ePot = (KQ)/(Dr) =

       8.987552e9 J*m     1.6021773e-19 C   1.6021773e-19 C      1 A        Q
       ---------------- * --------------- * --------------- * ---------- * --- =
       C^2              1 e               1 e            1e-10 m     Dr

       2.3070795237e-18 J*A      Q
       = ---------------------- * ---
       e^2               Dr

       Boltzmann Constant: k = 1.380658e-23 (J/K)
       Temperature: 300 Kelvin

       kT = 1.380658e-23 * 300 = 4.141974e-21 J

       2.3070795237e-18 J*A         1 kT             Q
       ePot = --------------------- * ------------------ * ---
       e^2              4.141974e-21 J       Dr

       557.00000 kT*A    Q
       ePot = -------------- * --- which will give kT/e units when applied
       e^2           Dr
     */
  }
}

int SettingCheckFontID(PyMOLGlobals* G, CSetting* set1, CSetting* set2,
                       int font_id) {
  int ret = font_id;
  if (font_id < 5) {  // we are no longer supporting GLUT labels since they are
                      // not resizeable
    PRINTFB(G, FB_Setting, FB_Warnings)
    "RepLabel-Warning: GLUT labels (label_font_id 0-4) are no longer available "
    "for labelling\n    the scene since they are not resizeable "
    "label_font_id=%d setting back to 5 (default) \n",
        font_id ENDFB(G);
    if (SettingGet_i(G, set1, nullptr, cSetting_label_font_id) == font_id &&
        SettingSet_i(set1, cSetting_label_font_id, 5)) {
    } else if (SettingGet_i(G, set2, nullptr, cSetting_label_font_id) ==
                   font_id &&
               SettingSet_i(set2, cSetting_label_font_id, 5)) {
    } else if (SettingGetGlobal_i(G, cSetting_label_font_id) == font_id) {
      SettingSetGlobal_i(G, cSetting_label_font_id, 5);
    };
    ret = 5;
  }
  return ret;
}

/**
 * State index iterator constructor, see Setting.h for documentation.
 *
 * @param set Optional object settings (can be nullptr)
 * @param nstate Maximum number of states
 */
StateIterator::StateIterator(PyMOLGlobals* G, CSetting* set,
                             StateIndex_t state_, int nstate) {
  if (state_ == cStateCurrent) {
    state_ = SettingGet_i(G, set, nullptr, cSetting_state) - 1;
  }

  if (state_ == cStateAll) {
    state = 0;
    end = nstate;
  } else {
    // given state or static singleton
    state = (state_ > 0 && nstate == 1 &&
             SettingGet_b(G, set, nullptr, cSetting_static_singletons))
                ? 0
                : state_;
    end = state + 1;
  }

  if (state < 0) state = 0;

  if (end > nstate) end = nstate;

  state--;
}

/**
 * Helper function to init CPyMOL.Setting, a (name: index) dictionary.
 * Called in PyMOL_InitAPI
 */
bool CPyMOLInitSetting(OVLexicon* Lex, std::unordered_map<int, int>& Setting) {
  for (int index = 0; index < cSetting_INIT; ++index) {
    auto& rec = SettingInfo[index];

    if (rec.level == cSettingLevel_unused) continue;

    OVreturn_word result = OVLexicon_GetFromCString(Lex, rec.name);

    if (!OVreturn_IS_OK(result)) {
      return false;
    }
    Setting[result.word] = index;
  }
  return true;
}

/**
 * Export the settings names to Python a as (name: index) dictionary.
 * Replacement for pymol.settings.SettingIndex
 */
PyObject* SettingGetSettingIndices() {
  PyObject* val;
  PyObject* dict = PyDict_New();

  for (int index = 0; index < cSetting_INIT; ++index) {
    auto& rec = SettingInfo[index];

    if (rec.level == cSettingLevel_unused) continue;

    if ((val = PyInt_FromLong(index))) {
      PyDict_SetItemString(dict, rec.name, val);
      Py_DECREF(val);
    }
  }

  return dict;
}

/**
 * Return a list of all setting indices for the given unique id
 */
PyObject* SettingUniqueGetIndicesAsPyList(PyMOLGlobals* G, int unique_id) {
  CSettingUnique* I = G->SettingUnique;
  PyObject* list = PyList_New(0);

  auto offsetIt = I->id2offset.find(unique_id);
  if (unique_id && offsetIt != I->id2offset.end()) {
    SettingUniqueEntry* entry;
    for (int offset = offsetIt->second; offset; offset = entry->next) {
      entry = &I->entry[offset];
      PyObject* item = PyInt_FromLong(entry->setting_id);
      PyList_Append(list, item);
      Py_DECREF(item);
    }
  }

  return list;
}

/*
 * Getters for templatted programming
 */

const CSetting* _SettingGetFirstDefined(int index, PyMOLGlobals* G,
                                        const CSetting* set1,
                                        const CSetting* set2) {
  if (set1 && set1->info[index].defined) return set1;
  if (set2 && set2->info[index].defined) return set2;
  return G->Setting;
}
