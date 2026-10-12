
/*
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* Copyright (c) Schrodinger, LLC.
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
#include "pymol/platform/os_std.h"
#include "pymol/platform/os_gl.h"

#include <clocale>

#include "pymol/platform/memory_debug.h"

#include "pymol/base/base.h"

#include "ov/ov_context.h"

#include "pymol/base/err.h"
#include "pymol/base/util.h"
#include "pymol/select/selector.h"
#include "pymol/color/color.h"
#include "pymol/ui/ortho.h"
#include "pymol/scene/scene.h"
#include "pymol/objects/pymol_object.h"
#include "pymol/executive/executive.h"
#include "pymol/base/word.h"
#include "pymol/reps/rep_mesh.h"
#include "pymol/objects/object_molecule.h"
#include "pymol/ui/control.h"
#include "pymol/geometry/sphere.h"
#include "pymol/executive/setting_effects.h"
#include "pymol/settings/setting.h"
#include "pymol/ray/ray.h"
#include "pymol/scene/movie.h"
#include "pymol/python/p.h"
#include "pymol/executive/editor.h"
#include "pymol/chem/sculpt_cache.h"
#include "pymol/fields/isosurf.h"
#include "pymol/fields/tetsurf.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/text/v_font.h"
#include "pymol/ui/wizard.h"
#include "pymol/text/text.h"
#include "pymol/glyph/character.h"
#include "pymol/ui/seq.h"
#include "pymol/ui/seeker.h"
#include "pymol/graphics/texture.h"
#include "pymol/app/test_pymol.h"
#include "pymol/glyph/type_face.h"
#include "pymol/io/plug_io_manager.h"
#include "pymol/scene/movie_scene.h"
#include "pymol/base/lex.h"
#include "pymol/select/selector_def.h"
#include "pymol/ui/but_mode.h"
#include "pymol/graphics/cgo_renderer.h"
#include "pymol/graphics/gfx_manager.h"

#include "pymol/objects/property.h"

#include "pymol/app/pymol.h"
#include "pymol/base/pymol_globals.h"
#include "pymol/base/pymol_options.h"
#include "pymol/base/feedback.h"
#include "pymol/graphics/graphics_util.h"
#include "pymol/zstring_view.h"

#include "pymol/graphics/shader_mgr.h"
#include "pymol/base/version.h"

#include <unordered_map>

PyMOLGlobals* SingletonPyMOLGlobals = nullptr;

#define PYMOL_API_LOCK if (!I->ModalDraw) {
#define PYMOL_API_LOCK_MODAL {
#define PYMOL_API_TRYLOCK if (!I->ModalDraw) {
#define PYMOL_API_UNLOCK }
#define PYMOL_API_UNLOCK_NO_FLUSH }
#define IDLE_AND_READY 3

struct CPyMOL {
  PyMOLGlobals* G = nullptr;
  int FakeDragFlag{};
  int RedisplayFlag{};
  int PassiveFlag{};
  int SwapFlag{};
  int BusyFlag{};
  int InterruptFlag{};
  int ReshapeFlag{};
  int ClickReadyFlag{};
  int DrawnFlag{};
  ObjectNameType ClickedObject{};
  int ClickedIndex{}, ClickedButton{}, ClickedModifiers{}, ClickedX{},
      ClickedY{}, ClickedHavePos{}, ClickedPosState{};
  int ClickedBondIndex{};
  float ClickedPos[3]{};
  int ImageRequestedFlag{}, ImageReadyFlag{};
  int DraggedFlag{};
  int Reshape[PYMOL_RESHAPE_SIZE]{};
  int Progress[PYMOL_PROGRESS_SIZE]{};
  int ProgressChanged{};
  int IdleAndReady{};
  bool done_ConfigureShaders{};

  PyMOLModalDrawFn* ModalDraw = nullptr;

  PyMOLSwapBuffersFn* SwapFn = nullptr;

  /* Python stuff */
  int PythonInitStage{};
  /* dynamically mapped string constants */

  OVLexicon* Lex = nullptr;
  std::unordered_map<int, int> Rep;
  ov_word lex_everything{}, lex_sticks{}, lex_spheres{}, lex_surface{};
  ov_word lex_labels{}, lex_nb_spheres{}, lex_cartoon{}, lex_ribbon{};
  ov_word lex_lines{}, lex_mesh{}, lex_dots{}, lex_dashes{}, lex_nonbonded{};
  ov_word lex_cell{}, lex_cgo{}, lex_callback{}, lex_extent{}, lex_slice{};

  std::unordered_map<int, int> Reinit;
  ov_word lex_settings;

  std::unordered_map<int, int> SelectList;
  ov_word lex_index{}, lex_id{}, lex_rank{};

  std::unordered_map<int, int> Setting;

  AtomPropertyInfo AtomPropertyInfos[NUM_ATOM_PROPERTIES]{};
  std::unordered_map<int, int> AtomPropertyLexicon;
  ov_word lex_atom_prop_model{}, lex_atom_prop_index{}, lex_atom_prop_type{},
      lex_atom_prop_name{}, lex_atom_prop_resn{}, lex_atom_prop_resi{},
      lex_atom_prop_resv{}, lex_atom_prop_chain{}, lex_atom_prop_alt{},
      lex_atom_prop_segi{}, lex_atom_prop_elem{}, lex_atom_prop_ss{},
      lex_atom_prop_text_type{}, lex_atom_prop_custom{}, lex_atom_prop_label{},
      lex_atom_prop_numeric_type{}, lex_atom_prop_q{}, lex_atom_prop_b,
      lex_atom_prop_vdw{}, lex_atom_prop_elec_radius{},
      lex_atom_prop_partial_charge{}, lex_atom_prop_formal_charge{},
      lex_atom_prop_stereo{}, lex_atom_prop_cartoon{}, lex_atom_prop_color{},
      lex_atom_prop_ID{}, lex_atom_prop_rank{}, lex_atom_prop_flags{},
      lex_atom_prop_geom{}, lex_atom_prop_valence{},
      lex_atom_prop_explicit_degree{}, lex_atom_prop_explicit_valence{},
      lex_atom_prop_x{}, lex_atom_prop_y, lex_atom_prop_z{},
      lex_atom_prop_settings{}, lex_atom_prop_properties{},
      lex_atom_prop_reps{}, lex_atom_prop_protons{}, lex_atom_prop_oneletter{},
      lex_atom_prop_s{}, lex_atom_prop_p{}, lex_atom_prop_state{};
  /*
    lex_atom_prop_, lex_atom_prop_, lex_atom_prop_,
    lex_atom_prop_, lex_atom_prop_, lex_atom_prop_,*/
};

/* convenience functions -- inline */

inline PyMOLstatus get_status_ok(int ok) {
  if (ok)
    return PyMOLstatus_SUCCESS;
  else
    return PyMOLstatus_FAILURE;
}

inline PyMOLreturn_status return_status_ok(int ok) {
  PyMOLreturn_status result;
  result.status = get_status_ok(ok);
  return result;
}

inline PyMOLreturn_status return_status(int status) {
  PyMOLreturn_status result;
  result.status = status;
  return result;
}

inline PyMOLreturn_float return_result(const pymol::Result<float>& res) {
  PyMOLreturn_float result = {PyMOLstatus_FAILURE};
  if (res) {
    result.status = PyMOLstatus_SUCCESS;
    result.value = res.result();
  }
  return result;
}

static PyMOLreturn_string_array return_result(
    const pymol::Result<std::vector<const char*>>& res) {
  PyMOLreturn_string_array result = {PyMOLstatus_SUCCESS, 0, nullptr};
  if (!res) {
    result.status = PyMOLstatus_FAILURE;
  } else if (!res.result().empty()) {
    const auto& vec = res.result();
    result.size = vec.size();
    result.array = VLAlloc(char*, result.size);

    // allocate space for concatenated string in first array element
    size_t reslen = 0;
    for (const char* s : vec) {
      reslen += strlen(s) + 1;
    }
    result.array[0] = VLAlloc(char, reslen);

    // copy elements
    for (size_t pl = 0, i = 0; i != vec.size(); ++i) {
      result.array[i] = result.array[0] + pl;
      strcpy(result.array[i], vec[i]);
      pl += strlen(vec[i]) + 1;
    }
  }
  return result;
}

pymol::Result<int> get_setting_id(CPyMOL* I, const char* setting) {
  OVreturn_word result;
  if (!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, setting))))
    return pymol::make_error(setting, " not found.");
  auto it = I->Setting.find(result.word);
  if (it == I->Setting.end()) {
    return pymol::make_error(setting, " not found.");
  }
  return it->second;
}

static pymol::Result<int> get_reinit_id(CPyMOL* I, const char* reinit) {
  OVreturn_word result;
  if (!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, reinit))))
    return pymol::make_error(reinit, " not found.");
  auto it = I->Reinit.find(result.word);
  if (it == I->Reinit.end()) {
    return pymol::make_error(reinit, " not found.");
  }
  return it->second;
}

static pymol::Result<int> get_select_list_mode(CPyMOL* I, const char* mode) {
  OVreturn_word result;
  if (!OVreturn_IS_OK((result = OVLexicon_BorrowFromCString(I->Lex, mode))))
    return pymol::make_error(mode, " not found.\n");
  auto it = I->SelectList.find(result.word);
  if (it == I->SelectList.end()) {
    return pymol::make_error(mode, " not found.\n");
  }
  return it->second;
}

extern "C" {

static OVstatus PyMOL_InitAPI(CPyMOL* I) {
  OVContext* C = I->G->Context;
  OVreturn_word result;
  I->Lex = OVLexicon_New(C->heap);
  if (!I->Lex) return_OVstatus_FAILURE;

  /* the following preprocessor macros may require GNU's cpp or VC++
     we'll see... */

#define LEX(ARG)                                                          \
  if (!OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, #ARG)))) \
    return_OVstatus_FAILURE else I->lex_##ARG = result.word;

  /* string constants that are accepted on input */

#define LEX_REP(NAME, CODE) \
  LEX(NAME)                 \
  I->Rep[I->lex_##NAME] = CODE;

  LEX_REP(everything, -1);
  LEX_REP(sticks, 0);
  LEX_REP(spheres, 1);
  LEX_REP(surface, 2);
  LEX_REP(labels, 3);
  LEX_REP(nb_spheres, 4);
  LEX_REP(cartoon, 5);
  LEX_REP(ribbon, 6);
  LEX_REP(lines, 7);
  LEX_REP(mesh, 8);
  LEX_REP(dots, 9);
  LEX_REP(dashes, 10);
  LEX_REP(nonbonded, 11);
  LEX_REP(cell, 12);
  LEX_REP(cgo, 13);
  LEX_REP(callback, 14);
  LEX_REP(extent, 15);
  LEX_REP(slice, 16);

#define LEX_REINIT(NAME, CODE)                                               \
  {                                                                          \
    if (!OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, #NAME)))) \
      return_OVstatus_FAILURE else I->lex_##NAME = result.word;              \
  }                                                                          \
  I->Reinit[I->lex_##NAME] = CODE;

  LEX_REINIT(everything, 0);
  LEX_REINIT(settings, 1);

#define LEX_SELLIST(NAME, CODE)                                              \
  {                                                                          \
    if (!OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, #NAME)))) \
      return_OVstatus_FAILURE else I->lex_##NAME = result.word;              \
  }                                                                          \
  I->SelectList[I->lex_##NAME] = CODE;

  LEX_SELLIST(index, 0);
  LEX_SELLIST(id, 1);
  LEX_SELLIST(rank, 2);

  if (!CPyMOLInitSetting(I->Lex, I->Setting)) return_OVstatus_FAILURE;

#define LEX_ATM_PROP(ARG)                                                 \
  if (!OVreturn_IS_OK((result = OVLexicon_GetFromCString(I->Lex, #ARG)))) \
    return_OVstatus_FAILURE else I->lex_atom_prop_##ARG = result.word;
#define LEX_ATOM_PROP(NAME, CODE, TYPE, OFFSET)           \
  LEX_ATM_PROP(NAME)                                      \
  I->AtomPropertyLexicon[I->lex_atom_prop_##NAME] = CODE; \
  I->AtomPropertyInfos[CODE].id = CODE;                   \
  I->AtomPropertyInfos[CODE].Ptype = TYPE;                \
  I->AtomPropertyInfos[CODE].offset = OFFSET;             \
  I->AtomPropertyInfos[CODE].maxlen = 0;

#define LEX_ATOM_PROP_S(NAME, CODE, TYPE, OFFSET, MAXLEN) \
  LEX_ATM_PROP(NAME)                                      \
  I->AtomPropertyLexicon[I->lex_atom_prop_##NAME] = CODE; \
  I->AtomPropertyInfos[CODE].id = CODE;                   \
  I->AtomPropertyInfos[CODE].Ptype = TYPE;                \
  I->AtomPropertyInfos[CODE].offset = OFFSET;             \
  I->AtomPropertyInfos[CODE].maxlen = MAXLEN;

  /*TEMP*/
  LEX_ATOM_PROP(model, 0, cPType_model, 0);
  LEX_ATOM_PROP(index, 1, cPType_index, 0);
  LEX_ATOM_PROP(type, 2, cPType_char_as_type, 0);
  LEX_ATOM_PROP(name, 3, cPType_int_as_string, offsetof(AtomInfoType, name));
  LEX_ATOM_PROP(resn, 4, cPType_int_as_string, offsetof(AtomInfoType, resn));
  LEX_ATOM_PROP(resi, 5, 0, 0);
  LEX_ATOM_PROP(resv, 6, cPType_int, offsetof(AtomInfoType, resv));
  LEX_ATOM_PROP(chain, 7, cPType_int_as_string, offsetof(AtomInfoType, chain));
  LEX_ATOM_PROP_S(alt, 8, cPType_string, offsetof(AtomInfoType, alt), 1);
  LEX_ATOM_PROP(segi, 9, cPType_int_as_string, offsetof(AtomInfoType, segi));
  LEX_ATOM_PROP_S(elem, 10, cPType_string, offsetof(AtomInfoType, elem),
                  cElemNameLen);
  LEX_ATOM_PROP_S(ss, 11, cPType_string, offsetof(AtomInfoType, ssType), 1);
  LEX_ATOM_PROP(text_type, 12, cPType_int_as_string,
                offsetof(AtomInfoType, textType));
  LEX_ATOM_PROP(custom, 13, cPType_int_as_string,
                offsetof(AtomInfoType, custom));
  LEX_ATOM_PROP(label, 14, cPType_int_as_string, offsetof(AtomInfoType, label));
  LEX_ATOM_PROP(numeric_type, 15, cPType_int_custom_type,
                offsetof(AtomInfoType, customType));
  LEX_ATOM_PROP(q, 16, cPType_float, offsetof(AtomInfoType, q));
  LEX_ATOM_PROP(b, 17, cPType_float, offsetof(AtomInfoType, b));
  LEX_ATOM_PROP(vdw, 18, cPType_float, offsetof(AtomInfoType, vdw));
  LEX_ATOM_PROP(elec_radius, 19, cPType_float,
                offsetof(AtomInfoType, elec_radius));
  LEX_ATOM_PROP(partial_charge, 20, cPType_float,
                offsetof(AtomInfoType, partialCharge));
  LEX_ATOM_PROP(formal_charge, 21, cPType_schar,
                offsetof(AtomInfoType, formalCharge));
  LEX_ATOM_PROP(stereo, 22, 0, 0);
  LEX_ATOM_PROP(cartoon, 23, cPType_schar, offsetof(AtomInfoType, cartoon));
  LEX_ATOM_PROP(color, 24, cPType_int, offsetof(AtomInfoType, color));
  LEX_ATOM_PROP(ID, 25, cPType_int, offsetof(AtomInfoType, id));
  LEX_ATOM_PROP(rank, 26, cPType_int, offsetof(AtomInfoType, rank));
  LEX_ATOM_PROP(flags, 27, cPType_uint32, offsetof(AtomInfoType, flags));
  LEX_ATOM_PROP(geom, 28, cPType_schar, offsetof(AtomInfoType, geom));
  LEX_ATOM_PROP(valence, 29, cPType_schar, offsetof(AtomInfoType, valence));
  LEX_ATOM_PROP(x, 30, cPType_xyz_float, 0);
  LEX_ATOM_PROP(y, 31, cPType_xyz_float, 1);
  LEX_ATOM_PROP(z, 32, cPType_xyz_float, 2);
  LEX_ATOM_PROP(settings, 33, cPType_settings, 0);
  LEX_ATOM_PROP(properties, 34, cPType_properties, 0);
  LEX_ATOM_PROP(s, 35, cPType_settings, 0);
  LEX_ATOM_PROP(p, 36, cPType_properties, 0);
  LEX_ATOM_PROP(state, 37, cPType_state, 0);
  LEX_ATOM_PROP(reps, 38, cPType_int, offsetof(AtomInfoType, visRep));
  LEX_ATOM_PROP(protons, 39, cPType_schar, offsetof(AtomInfoType, protons));
  LEX_ATOM_PROP(oneletter, 40, 0, 0);
  LEX_ATOM_PROP(explicit_degree, ATOM_PROP_EXPLICIT_DEGREE, 0, 0);
  LEX_ATOM_PROP(explicit_valence, ATOM_PROP_EXPLICIT_VALENCE, 0, 0);
  //  LEX_ATOM_PROP(, );

  return_OVstatus_SUCCESS;
}

static OVstatus PyMOL_PurgeAPI(CPyMOL* I) {
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  return_OVstatus_SUCCESS;
}

static int word_count(const char* src) { /* only works for ascii */
  int cnt = 0;
  while ((*src) && ((*src) < 33)) /* skip leading whitespace */
    src++;
  while (*src) {
    if ((*src) > 32) {
      cnt++;
      while ((*src) && ((*src) > 32)) src++;
    }
    while ((*src) && ((*src) < 33)) src++;
  }
  return cnt;
}

static const char* next_word(const char* src, char* dst,
                             int buf_size) { /* only works for ascii */
  while ((*src) && ((*src) < 33))            /* skip leading whitespace */
    src++;
  while (*src) {
    if ((*src) > 32) {
      while ((*src) && ((*src) > 32) && (buf_size > 1)) {
        *(dst++) = *(src++);
        buf_size--;
      }
      break;
    }
  }
  dst[0] = 0;
  return src;
}

/**
 * Supported file formats and their internal codes
 */
struct {
  const char* name;
  cLoadType_t code_buffer;
  cLoadType_t code_filename;
} const ContentTypeTable[] = {
    // molecules
    {"pdb", cLoadTypePDBStr, cLoadTypePDB},
    {"vdb", cLoadTypeVDBStr, cLoadTypeUnknown},
    {"cif", cLoadTypeCIFStr, cLoadTypeCIF},
    {"mae", cLoadTypeMAEStr, cLoadTypeMAE},
    {"sdf", cLoadTypeSDF2Str, cLoadTypeSDF2},
    {"mol", cLoadTypeMOLStr, cLoadTypeMOL},
    {"mol2", cLoadTypeMOL2Str, cLoadTypeMOL2},
    {"xyz", cLoadTypeXYZStr, cLoadTypeXYZ},
    {"pqr", cLoadTypeUnknown, cLoadTypePQR},
    {"macromodel", cLoadTypeMMDStr, cLoadTypeMMD},
    // maps
    {"ccp4", cLoadTypeCCP4Str, cLoadTypeCCP4Map},
    {"mrc", cLoadTypeMRCStr, cLoadTypeMRC},
    {"map", cLoadTypeCCP4UnspecifiedStr, cLoadTypeCCP4Unspecified},
    {"xplor", cLoadTypeXPLORStr, cLoadTypeXPLORMap},
    {"phi", cLoadTypePHIStr, cLoadTypePHIMap},
    {"dx", cLoadTypeDXStr, cLoadTypeDXMap},
    // special
    {"cgo", cLoadTypeCGO, cLoadTypeUnknown},
    {NULL, cLoadTypeUnknown, cLoadTypeUnknown}};

/**
 * Proxy for "ExecutiveLoad" with string "content_format" (and "content_type")
 * argument.
 *
 * content:     Either file name or file contents, depending on "content_type"
 * content_type:        "filename", "string", "raw", or "cgo"
 * content_length:      Length of "content", if it's not a file name or a
 *                      null-terminated string (pass -1).
 * content_format:      The file format, e.g. "pdb", "sdf", "mol2", ...
 * object_name:         New object name. Can be empty if "content_type" is
 *                      "filename".
 */
static PyMOLreturn_status Loader(CPyMOL* I, const char* content,
                                 const char* content_type, int content_length,
                                 const char* content_format,
                                 const char* object_name, int state,
                                 int discrete, int finish, int quiet,
                                 int multiplex, int zoom) {
  PyMOLGlobals* G = I->G;
  bool content_is_filename = false;
  int ok = true;
  WordType obj_name;

  // `content` can be a file name, or the file contents
  if (strcmp(content_type, "filename") == 0) {
    content_is_filename = true;
  } else if (strcmp(content_type, "string") == 0) {
    if (content_length < 0) content_length = strlen(content);
  } else if (strcmp(content_type, "raw") != 0 &&
             strcmp(content_type, "cgo") != 0) {
    PRINTFB(G, FB_Executive, FB_Errors)
    " Error: Unknown content type '%s'\n", content_type ENDFB(G);
    ok = false;
  }

  if (ok) {
    { /* if object_name is blank and content is a filename, then
         compute the object_name from the file prefix */
      if ((!object_name[0]) && content_is_filename) {
        const char *start, *stop;
        stop = start = content + strlen(content) - 1;
        while (start > content) { /* known path separators */
          if ((start[-1] == ':') || (start[-1] == '\'') || (start[-1] == '/'))
            break;
          start--;
        }
        while (stop > start) {
          if (*stop == '.') break;
          stop--;
        }
        if (stop == start) stop = content + strlen(content);
        if ((stop - start) >= sizeof(WordType))
          stop = start + sizeof(WordType) - 1;
        {
          char* q;
          const char* p = start;
          q = obj_name;
          while (p < stop) {
            *(q++) = *(p++);
          }
          *q = 0;
          object_name = obj_name;
        }
      }
    }
    {
      cLoadType_t pymol_content_type = cLoadTypeUnknown;

      for (auto it = ContentTypeTable; it->name; ++it) {
        if (strcmp(it->name, content_format) == 0) {
          pymol_content_type =
              content_is_filename ? it->code_filename : it->code_buffer;
          break;
        }
      }

      if (pymol_content_type == cLoadTypeUnknown) {
        PRINTFB(G, FB_Executive, FB_Errors)
        " Error: Unknown content format '%s' with type '%s'\n", content_format,
            content_type ENDFB(G);
        ok = false;
      }

      if (ok) {
        auto result = ExecutiveLoad(
            I->G, content_is_filename ? content : nullptr,
            content_is_filename ? nullptr : content, content_length,
            pymol_content_type, object_name, state - 1, zoom, discrete, finish,
            multiplex, quiet, nullptr, 0, nullptr);
        ok = static_cast<bool>(result);
      }
    }
  }
  if (ok) PyMOL_NeedRedisplay(I);
  return return_status_ok(ok);
}

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL* I, const char* content,
                                 const char* content_type,
                                 const char* content_format,
                                 const char* object_name, int state,
                                 int discrete, int finish, int quiet,
                                 int multiplex, int zoom) {
  PyMOLreturn_status status = {PyMOLstatus_FAILURE};
  PYMOL_API_LOCK
  status = Loader(I, content, content_type, -1, content_format, object_name,
                  state, discrete, finish, quiet, multiplex, zoom);
  PYMOL_API_UNLOCK return status;
}

static const CPyMOLOptions Defaults = {
    true,  /* pmgui */
    true,  /* internal_gui */
    true,  /* show_splash */
    1,     /* internal_feedback */
    true,  /* security */
    false, /* game mode */
    0,     /* force_stereo */
    640,   /* winX */
    480,   /* winY */
    false, /* blue_line */
    0,     /* winPX */
    175,   /* winPY */
    true,  /* external_gui */
    true,  /* siginthand */
    false, /* reuse helper */
    false, /* auto reinitialize */
    false, /* keep thread alive */
    false, /* quiet */
    false, /* incentive product */
    "",    /* after_load_script */
    0,     /* multisample */
    1,     /* window_visible */
    0,     /* read_stdin */
    0,     /* presentation */
    0,     /* defer builds mode */
    0,     /* full screen mode */
    -1,    /* sphere mode */
    0,     /* stereo capable */
    0,     /* stereo mode */
    -1,    /* zoom mode */
    0,     /* launch_status */
    0,     /* no quit */
    0,     /* gldebug */
    false, /* no openvr stub */
};

CPyMOLOptions* PyMOLOptions_New(void) {
  CPyMOLOptions* result = nullptr;
  result = pymol::calloc<CPyMOLOptions>(1);
  if (result) *result = Defaults;
  return result;
}

void PyMOLOptions_Free(CPyMOLOptions* options) { FreeP(options); }

void PyMOL_ResetProgress(CPyMOL* I) {
  I->ProgressChanged = true;
  UtilZeroMem(I->Progress, sizeof(int) * 6);
}

void PyMOL_SetProgress(CPyMOL* I, int offset, int current, int range) {
  switch (offset) {
    case PYMOL_PROGRESS_SLOW:
    case PYMOL_PROGRESS_MED:
    case PYMOL_PROGRESS_FAST:
      if (current != I->Progress[offset]) {
        I->Progress[offset] = current;
        I->ProgressChanged = true;
      }
      if (range != I->Progress[offset + 1]) {
        I->Progress[offset + 1] = range;
        I->ProgressChanged = true;
      }
  }
}

int PyMOL_GetProgress(CPyMOL* I, int* progress, int reset) {
  int a;
  int result = I->ProgressChanged;
  for (a = 0; a < PYMOL_PROGRESS_SIZE; a++) {
    progress[a] = I->Progress[a];
  }
  if (reset) I->ProgressChanged = false;
  return result;
}

CPyMOL* PyMOL_NewWithOptions(const CPyMOLOptions* option) {
  auto result = new CPyMOL();
  assert(result);

  auto G = pymol::calloc<PyMOLGlobals>(1);
  assert(G);

  result->G = G;
  G->PyMOL = result;

  PyMOL_ResetProgress(result);

  G->Option = pymol::calloc<CPyMOLOptions>(1);
  assert(G->Option);

  if (!option) {
    assert(!Defaults.stereo_capable);
    option = &Defaults;
  }

  *(G->Option) = *option;

  G->Security = option->security;
  G->StereoCapable = option->stereo_capable;

  return result;
}

void PyMOL_Start(CPyMOL* I) {
  PyMOLGlobals* G = I->G;

  // It's possible to change this from Python, functions which rely on
  // C locale should reset it before doing printf, atof, etc.
  std::setlocale(LC_NUMERIC, "C");

  G->Context = OVContext_New();
  G->Lexicon = OVLexicon_New(G->Context->heap);

  if (OVreturn_IS_ERROR(PyMOL_InitAPI(I))) {
    printf("ERROR: PyMOL internal C API initialization failed.\n");
  }

  // global lexicon "constants"
#define LEX_CONSTANTS_IMPL
#include "pymol/base/lex_constants.h"

  G->Feedback = new CFeedback(G, G->Option->quiet);
  WordInit(G);
  UtilInit(G);
  ColorInit(G);
  CGORendererInit(G);
  ShaderMgrInit(G);
  G->GFXMgr = new GFXManager(G->ShaderMgr);
  SettingInitGlobal(G, true, true, false);
  PropertyInit(G);
  SettingSetGlobal_i(G, cSetting_internal_gui, G->Option->internal_gui);
  SettingSetGlobal_i(G, cSetting_internal_feedback,
                     G->Option->internal_feedback);
  TextureInit(G);
  TypeInit(G);
  TextInit(G);
  CharacterInit(G);
  PlugIOManagerInit(G);
  SphereInit(G);
  // OpenVRInit() called in ExecutiveStereo
  OrthoInit(G, G->Option->show_splash);
  SceneInit(G);
  MovieScenesInit(G);
  WizardInit(G); /* must come after ortho & scene */
  G->Movie = new CMovie(G);
  G->SelectorMgr = new CSelectorManager();
  G->Selector = new CSelector(G, G->SelectorMgr);
  SeqInit(G);
  SeekerInit(G);
  ButModeInit(G);
  ControlInit(G);
  AtomInfoInit(G);
  SculptCacheInit(G);
  VFontInit(G);
  ExecutiveInit(G);
  IsosurfInit(G);
  TetsurfInit(G);
  EditorInit(G);

  I->DrawnFlag = false;
  I->RedisplayFlag = true;
  G->Ready = true;
}

void PyMOL_StartWithPython(CPyMOL* I) {
  PyMOL_Start(I);

  /* now locate all the C to Python function hooks and objects we need */

  PInit(I->G, false);

  /* and begin the initialization sequence */

  I->PythonInitStage = 1;
}

/* Set the Python initialization stage.  Stage 1 enables exec_deferred()
   in PyMOL_Idle, which processes queued commands such as config_mouse().
   Used by embedding hosts (e.g. AppKit) that call PyMOL_Start and PInit
   separately instead of using PyMOL_StartWithPython. */
void PyMOL_SetPythonInitStage(CPyMOL* I, int stage) {
  I->PythonInitStage = stage;
}

void PyMOL_Stop(CPyMOL* I) {
  PyMOLGlobals* G = I->G;
  G->Terminating = true;
  TetsurfFree(G);
  IsosurfFree(G);
  WizardFree(G);
  EditorFree(G);
  ExecutiveFree(G);
  VFontFree(G);
  SculptCacheFree(G);
  AtomInfoFree(G);
  ButModeFree(G);
  ControlFree(G);
  SeekerFree(G);
  SeqFree(G);
  DeleteP(G->Selector);
  DeleteP(G->SelectorMgr);
  DeleteP(G->Movie);
  SceneFree(G);
  MovieScenesFree(G);
  OrthoFree(G);
  DeleteP(G->GFXMgr);
  DeleteP(G->ShaderMgr);
  SettingFreeGlobal(G);
  CharacterFree(G);
  TextFree(G);
  TypeFree(G);
  TextureFree(G);
  SphereFree(G);
  PlugIOManagerFree(G);
  PropertyFree(G);
  PFree(G);
  CGORendererFree(G);
  ColorFree(G);
  UtilFree(G);
  WordFree(G);
  DeleteP(G->Feedback);

  PyMOL_PurgeAPI(I);
  /*    printf("%d \n", OVLexicon_GetNActive(G->Lexicon)); */
  OVLexicon_Del(G->Lexicon);
  OVContext_Del(G->Context);
}

void PyMOL_Free(CPyMOL* I) {
#if !defined(_PYMOL_ACTIVEX) && !defined(_MACPYMOL_XCODE)
  PYMOL_API_LOCK
#endif
  /* take PyMOL down gracefully */
  PyMOLOptions_Free(I->G->Option);

  FreeP(I->G->P_inst);
  if (I->G == SingletonPyMOLGlobals) SingletonPyMOLGlobals = nullptr;

  FreeP(I->G);
  FreeP(I);
  return;
#if !defined(_PYMOL_ACTIVEX) && !defined(_MACPYMOL_XCODE)
  PYMOL_API_UNLOCK;
#endif
}

struct PyMOLGlobals* PyMOL_GetGlobals(CPyMOL* I) { return I->G; }

struct PyMOLGlobals** PyMOL_GetGlobalsHandle(CPyMOL* I) { return &(I->G); }

void PyMOL_DrawWithoutLock(CPyMOL* I);

void PyMOL_Draw(CPyMOL* I) {
  PYMOL_API_LOCK_MODAL
  PyMOL_DrawWithoutLock(I);
  PYMOL_API_UNLOCK
}

void PyMOL_DrawWithoutLock(CPyMOL* I) {
  // Metal draws the scene (RendererMetal, from drawInMTKView); this only does
  // the one-time setup.
  if (!I->done_ConfigureShaders) {
    I->done_ConfigureShaders = true;
    I->G->HaveGUI = I->G->Option->pmgui;
    // Shader-style CGOs build the VBOs that Metal renders.
    SettingSetGlobal_b(I->G, cSetting_use_shaders, true);
  }
  I->DrawnFlag = true;
  I->RedisplayFlag = false;
}

void PyMOL_Key(CPyMOL* I, unsigned char k, int x, int y, int modifiers) {
  PYMOL_API_LOCK PyMOLGlobals* G = I->G;
  if (!WizardDoKey(G, k, x, y, modifiers)) OrthoKey(G, k, x, y, modifiers);
  PyMOL_NeedRedisplay(G->PyMOL);
  PYMOL_API_UNLOCK
}

void PyMOL_Special(CPyMOL* I, int k, int x, int y, int modifiers) {
  PYMOL_API_LOCK PyMOLGlobals* G = I->G;

  int grabbed = false;
  char buffer[255];
  (void)buffer;
  if (!grabbed) grabbed = WizardDoSpecial(G, (unsigned char)k, x, y, modifiers);

  switch (k) {
    case P_GLUT_KEY_UP:
    case P_GLUT_KEY_DOWN:
      grabbed = 1;
      OrthoSpecial(G, k, x, y, modifiers);
      break;
    case P_GLUT_KEY_LEFT:
    case P_GLUT_KEY_RIGHT:
      if (OrthoArrowsGrabbed(G)) {
        grabbed = 1;
        OrthoSpecial(G, k, x, y, modifiers);
      }
      break;
  }

  if (!grabbed) {
    sprintf(buffer, "_special %d,%d,%d,%d", k, x, y, modifiers);
    PLog(G, buffer, cPLog_pml);
    PParse(G, buffer);
    PFlush(G);
  }

  PYMOL_API_UNLOCK
}

void PyMOL_Reshape(CPyMOL* I, int width, int height, int force) {
  PYMOL_API_LOCK PyMOLGlobals* G = I->G;

  G->Option->winX = width;
  G->Option->winY = height;

  OrthoReshape(G, width, height, force);
  PYMOL_API_UNLOCK
}

int PyMOL_Idle(CPyMOL* I) {
  int did_work = false;
  PYMOL_API_TRYLOCK PyMOLGlobals* G = I->G;

  I->DraggedFlag = false;
  if (I->IdleAndReady < IDLE_AND_READY) {
    if (I->DrawnFlag) I->IdleAndReady++;
  }
  if (I->FakeDragFlag == 1) {
    I->FakeDragFlag = false;
    OrthoFakeDrag(G);
    did_work = true;
  }

  if (ControlIdling(G)) {
    ExecutiveSculptIterateAll(G);
    ControlSdofIterate(G);
    did_work = true;
  }

  SceneIdle(G);

  if (SceneRovingCheckDirty(G)) {
    SceneRovingUpdate(G);
    did_work = true;
  }
  if (PFlush(G)) {
    did_work = true;
  }

  if (I->PythonInitStage > 0) {
    if (I->PythonInitStage < 2) {
      I->PythonInitStage++;
    } else {
      I->PythonInitStage = -1;
      PBlock(G);

      PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->obj, "adapt_to_hardware", "O",
                                   G->P_inst->obj));

      if (PyErr_Occurred()) PyErr_Print();

      PXDecRef(PYOBJECT_CALLMETHOD(G->P_inst->obj, "exec_deferred", "O",
                                   G->P_inst->obj));

      if (PyErr_Occurred()) PyErr_Print();

      PUnblock(G);
      PFlush(G);
    }
  }

  /* reset the interrupt flag if we're not doing anything */

  if (!(did_work || I->ModalDraw))
    if (PyMOL_GetInterrupt(I, false)) PyMOL_SetInterrupt(I, false);

  PYMOL_API_UNLOCK_NO_FLUSH return (did_work || I->ModalDraw);
}

void PyMOL_NeedFakeDrag(CPyMOL* I) { I->FakeDragFlag = true; }

void PyMOL_NeedRedisplay(CPyMOL* I) { I->RedisplayFlag = true; }

void PyMOL_NeedSwap(CPyMOL* I) { I->SwapFlag = true; }

void PyMOL_NeedReshape(CPyMOL* I, int mode, int x, int y, int width,
                       int height) {
  PyMOLGlobals* G = I->G;
  if (width < 0) {
    if (!G->HaveGUI) return;

    width = SceneGetBlock(G)->getWidth();
    if (SettingGetGlobal_b(G, cSetting_internal_gui))
      width += DIP2PIXEL(SettingGetGlobal_i(G, cSetting_internal_gui_width));
  }

  /* if height is negative, force a reshape based on the current height */

  if (height < 0) {
    int internal_feedback;
    height = SceneGetBlock(G)->getHeight();
    internal_feedback = SettingGetGlobal_i(G, cSetting_internal_feedback);
    if (internal_feedback)
      height +=
          (internal_feedback - 1) * cOrthoLineHeight + cOrthoBottomSceneMargin;
    if (SettingGetGlobal_b(G, cSetting_seq_view) &&
        !SettingGetGlobal_b(G, cSetting_seq_view_overlay)) {
      height += SeqGetHeight(G);
    }
    height += MovieGetPanelHeight(G);
  }

  if (G->HaveGUI) {
    float sf = DIP2PIXEL(1);

    I->Reshape[1] = x / sf;
    I->Reshape[2] = y / sf;
    I->Reshape[3] = width / sf;
    I->Reshape[4] = height / sf;

    I->ReshapeFlag = true;
    I->Reshape[0] = mode;
    PyMOL_NeedRedisplay(I);
  } else {
    /* if no gui, then force immediate reshape */
    PyMOLGlobals* G = I->G;

    G->Option->winX = width;
    G->Option->winY = height;

    OrthoReshape(G, width, height, true);
  }
}

int PyMOL_GetIdleAndReady(CPyMOL* I) {
  return (I->IdleAndReady == IDLE_AND_READY);
}

void PyMOL_SetPassive(CPyMOL* I, int onOff) { I->PassiveFlag = onOff; }

void PyMOL_SetClickReady(CPyMOL* I, const char* name, int index, int button,
                         int mod, int x, int y, const float* pos, int state,
                         int bond) {
  I->ClickReadyFlag = true;
  I->ClickedIndex = index;
  I->ClickedBondIndex = bond;
  I->ClickedButton = button;
  I->ClickedModifiers = mod;
  I->ClickedX = x;
  I->ClickedY = y;
  I->ClickedPosState = state;

  strcpy(I->ClickedObject, name ? name : "");

  if ((I->ClickedHavePos = bool(pos))) {
    copy3f(pos, I->ClickedPos);
  } else {
    zero3f(I->ClickedPos);
  }
}

int PyMOL_GetRedisplay(CPyMOL* I, int reset) {
  int result = false;

  PYMOL_API_TRYLOCK PyMOLGlobals* G = I->G;
  result = I->RedisplayFlag;

  if (result) {
    if (SettingGet_b(G, nullptr, nullptr, cSetting_defer_updates)) {
      result = false;
    } else {
      if (reset) I->RedisplayFlag = false;
    }
  }
  PYMOL_API_UNLOCK_NO_FLUSH return (
      result || I->ModalDraw); /* always true when ModalDraw is set */
}

int PyMOL_GetModalDraw(CPyMOL* I) {
  if (I) return (I->ModalDraw != nullptr);
  return false;
}

void PyMOL_SetModalDraw(CPyMOL* I, PyMOLModalDrawFn* fn) { I->ModalDraw = fn; }

int PyMOL_GetBusy(CPyMOL* I, int reset) { /* lock intentionally omitted */
  int result = I->BusyFlag;
  if (reset) PyMOL_SetBusy(I, false);
  return result;
}

void PyMOL_SetBusy(CPyMOL* I, int value) { /* lock intentionally omitted */
  if (!I->BusyFlag)
    /* if we weren't busy before, then reset the progress indicators */
    PyMOL_ResetProgress(I);

  I->BusyFlag = value;
}

int PyMOL_GetInterrupt(CPyMOL* I, int reset) { /* lock intentionally omitted */
  if (I) {
    int result = I->InterruptFlag;

    if (reset) PyMOL_SetInterrupt(I, false);
    return result;
  } else
    return false;
}

void PyMOL_SetInterrupt(CPyMOL* I, int value) { /* lock intentionally omitted */
  if (I) {
    I->InterruptFlag = value;
    if (I->G) I->G->Interrupt = value;
  }
}

void PyMOL_Drag(CPyMOL* I, int x, int y, int modifiers) {
  PYMOL_API_LOCK OrthoDrag(I->G, x, y, modifiers);
  I->DraggedFlag = true;
  PYMOL_API_UNLOCK
}

/**
 * Mouse button and keyboard press handler
 *
 * button: mouse button or key code
 * state:
 *   -2 = key press with GLUT_KEY_* special code
 *   -1 = key press with ascii code
 *    0 = mouse down
 *    1 = mouse up
 * x, y: mouse pointer position
 * modifiers: SHIFT/CTRL/ALT bitmask
 */
void PyMOL_Button(CPyMOL* I, int button, int state, int x, int y,
                  int modifiers) {
  PYMOL_API_LOCK
  if (state == -1) {
    PyMOL_Key(I, (unsigned char)button, x, y, modifiers);
  } else if (state == -2) {
    PyMOL_Special(I, button, x, y, modifiers);
  } else {
    OrthoButton(I->G, button, state, x, y, modifiers);
  }
  PYMOL_API_UNLOCK
}

void PyMOL_SwapBuffers(CPyMOL* I) {
  if (I->SwapFn && I->G->ValidContext) {
    I->SwapFn();
    I->SwapFlag = false;
  } else {
    I->SwapFlag = true;
  }
}

void PyMOL_PushValidContext(CPyMOL* I) {
  if (I && I->G) I->G->ValidContext++;
}

void PyMOL_PopValidContext(CPyMOL* I) {
  if (I && I->G && (I->G->ValidContext > 0)) I->G->ValidContext--;
}

void PyMOL_SetDefaultMouse(CPyMOL* I) {
  PYMOL_API_LOCK PyMOLGlobals* G = I->G;

  ButModeSet(G, cButModeLeftNone, cButModeRotXYZ);
  ButModeSet(G, cButModeMiddleNone, cButModeTransXY);
  ButModeSet(G, cButModeRightNone, cButModeTransZ);

  ButModeSet(G, cButModeLeftShft, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleShft, cButModePotentialClick);
  ButModeSet(G, cButModeRightShft, cButModeClipNF);

  ButModeSet(G, cButModeLeftCtrl, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleCtrl, cButModePotentialClick);
  ButModeSet(G, cButModeRightCtrl, cButModePotentialClick);

  ButModeSet(G, cButModeLeftCtSh, cButModePotentialClick);
  ButModeSet(G, cButModeMiddleCtSh, cButModePotentialClick);
  ButModeSet(G, cButModeRightCtSh, cButModePotentialClick);

  ButModeSet(G, cButModeWheelNone, cButModeScaleSlab);
  ButModeSet(G, cButModeWheelShft, cButModeMoveSlab);
  ButModeSet(G, cButModeWheelCtrl, cButModeMoveSlabAndZoom);
  ButModeSet(G, cButModeWheelCtSh, cButModeTransZ);

  ButModeSet(G, cButModeMiddleCtSh, cButModeOrigAt); /* SET TWICE?!? */

  ButModeSet(G, cButModeLeftSingle, cButModeSimpleClick);
  ButModeSet(G, cButModeMiddleSingle, cButModeCent);
  ButModeSet(G, cButModeRightSingle, cButModeSimpleClick);

  ButModeSet(G, cButModeLeftDouble, cButModeSimpleClick);
  ButModeSet(G, cButModeRightDouble, cButModeSimpleClick);

  {
    int a;
    for (a = cButModeLeftShftDouble; a <= cButModeRightCtrlAltShftSingle; a++) {
      ButModeSet(G, a, cButModeSimpleClick);
    }
    for (a = cButModeLeftAlt; a <= cButModeRightCtrlAltShft; a++) {
      ButModeSet(G, a, cButModePotentialClick);
    }
  }
  G->Feedback->currentMask(FB_Scene) &=
      ~(FB_Results); /* suppress click messages */
  PYMOL_API_UNLOCK
}

AtomPropertyInfo* PyMOL_GetAtomPropertyInfo(CPyMOL* I,
                                            const char* atompropname) {
  OVreturn_word result;
  if (!OVreturn_IS_OK(
          (result = OVLexicon_BorrowFromCString(I->Lex, atompropname))))
    return nullptr;
  auto it = I->AtomPropertyLexicon.find(result.word);
  if (it == I->AtomPropertyLexicon.end()) {
    return nullptr;
  }
  return &I->AtomPropertyInfos[it->second];
}

}  // extern "C"
