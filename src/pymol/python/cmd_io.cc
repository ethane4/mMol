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

// Python bindings for loading and saving (pymol/importing.py, exporting.py).

#include "pymol/executive/executive.h"
#include "pymol/ray/scene_ray.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/io/plug_io_manager.h"
#include "pymol/io/molecule_exporter.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdFindMolfilePlugin(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* ext = nullptr;
  int mask = 0;
  if (!PyArg_ParseTuple(args, "Os|i", &self, &ext, &mask)) {
    API_HANDLE_ERROR;
  } else {
    API_SETUP_PYMOL_GLOBALS;
    if (G) {
      APIEnter(G);
      const char* plugin = PlugIOManagerFindPluginByExt(G, ext, mask);
      APIExit(G);
      return PyUnicode_FromString(plugin ? plugin : "");
    }
  }
  return APIAutoNone(nullptr);
}

static PyObject* CmdGetCCP4Str(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* name = nullptr;
  int state = 0;
  int quiet = 1;
  int format = cLoadTypeCCP4Unspecified;
  if (!PyArg_ParseTuple(args, "Osii|i", &self, &name, &state, &quiet,
                        &format)) {
    API_HANDLE_ERROR;
  } else {
    API_SETUP_PYMOL_GLOBALS;
    if (G) {
      APIEnterBlocked(G);
      auto v = ObjectMapGetCCP4Str(G, name, state, quiet, format);
      PyObject* result =
          v.empty() ? nullptr : PyBytes_FromStringAndSize(&v.front(), v.size());

      APIExitBlocked(G);
      return APIAutoNone(result);
    }
  }
  return APIAutoNone(nullptr);
}

static PyObject* CmdGetSession(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* dict;
  int partial, quiet;
  const char* names;
  int binary = -1;
  float version = -1.f;

  API_SETUP_ARGS(G, self, args, "OOsii|if", &self, &dict, &names, &partial,
                 &quiet, &binary, &version);
  API_ASSERT(-1 <= binary && binary <= 1);

  APIEnterBlocked(G);

  const auto binary_orig = SettingGet<bool>(G, cSetting_pse_binary_dump);
  if (binary != -1) SettingSet(G, cSetting_pse_binary_dump, bool(binary));

  const auto version_orig = SettingGet<float>(G, cSetting_pse_export_version);
  if (version >= 0.f) SettingSet(G, cSetting_pse_export_version, version);

  ExecutiveGetSession(G, dict, names, partial, quiet);

  SettingSet(G, cSetting_pse_binary_dump, binary_orig);
  SettingSet(G, cSetting_pse_export_version, version_orig);

  APIExitBlocked(G);

  if (PyErr_Occurred()) {
    return nullptr;
  }

  return APISuccess();
}

static PyObject* CmdSetSession(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int quiet, partial;
  PyObject* obj;
  API_SETUP_ARGS(G, self, args, "OOii", &self, &obj, &partial, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  bool ok = ExecutiveSetSession(G, obj, partial, quiet);
  APIExitBlocked(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdLoadColorTable(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* str1;
  float gamma;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osfi", &self, &str1, &gamma, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  bool ok = ColorTableLoad(G, str1, gamma, quiet);
  APIExit(G);
  return APIResultOk(G, ok);
}

static PyObject* CmdGetSeqAlignStr(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  char* seq = nullptr;
  int state;
  int format;
  int quiet;
  PyObject* result = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Osiii", &self, &str1, &state, &format, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    APIEnter(G);
    seq = ExecutiveNameToSeqAlignStrVLA(G, str1, state, format, quiet);
    APIExit(G);
    if (seq) result = Py_BuildValue("s", seq);
    VLAFreeP(seq);
  }
  return (APIAutoNone(result));
}

static PyObject* CmdGetStr(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  pymol::vla<char> vla;
  char* format;
  char* sele;
  int state;
  char* ref;
  int ref_state;
  int quiet;
  int multi;

  API_SETUP_ARGS(G, self, args, "Ossisiii", &self, &format, &sele, &state, &ref,
                 &ref_state, &multi, &quiet);
  APIEnter(G);
  vla = MoleculeExporterGetStr(G, format, sele, state, ref, ref_state, multi,
                               quiet);
  APIExit(G);

  if (vla) {
    result = PyBytes_FromStringAndSize(vla, vla.size());
  }

  return APIAutoNone(result);
}

static PyObject* CmdPNG(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* filename = nullptr;
  int quiet;
  int result = 0;
  int width, height, ray;
  int prior, format;
  float dpi;

  API_SETUP_ARGS(G, self, args, "Oziifiiii", &self, &filename, &width, &height,
                 &dpi, &ray, &quiet, &prior, &format);
  API_ASSERT(APIEnterNotModal(G));

  // if `filename` is None, then return a PNG buffer
  std::vector<unsigned char> pngbuf;
  pymol::null_safe_zstring_view fileview(filename);

  {
    // with prior=1 other arguments (width, height, ray) are ignored

    if (!prior) {
      if (ray || (!G->HaveGUI && (!SceneGetCopyType(G) || width || height))) {
        prior =
            SceneRay(G, width, height,
                     SettingGetGlobal_i(G, cSetting_ray_default_renderer),
                     nullptr, nullptr, 0.0F, 0.0F, quiet, nullptr, true, -1);
      } else if (width || height) {
        Extent2D extent{static_cast<std::uint32_t>(width),
                        static_cast<std::uint32_t>(height)};
        bool withOverlay = true;  // TODO: Retrieve from python
        prior = !OrthoDeferImage(G, extent, fileview.c_str(), -1, dpi, format,
                                 quiet, nullptr, withOverlay);
        result = bool(filename);
      } else if (!SceneGetCopyType(G)) {
        ExecutiveDrawNow(G); /* TODO STATUS */
      }
    }

    if (!result) {
      auto outbuf = !fileview.empty() ? nullptr : &pngbuf;
      if (ScenePNG(G, fileview.c_str(), dpi, quiet, prior, format, outbuf)) {
        /* signal success by returning 1 instead of 0, or -1 for error */
        result = 1;
      }
    }
    APIExit(G);
  }

  if (fileview.empty()) {
    if (pngbuf.empty()) {
      return APIFailure(G, "getting png buffer failed");
    }

    return PyBytes_FromStringAndSize(
        reinterpret_cast<const char*>(pngbuf.data()), pngbuf.size());
  }

  return APIResultCode(result);
}

static PyObject* CmdLoadCoords(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state = 0;
  OrthoLineType s1;
  PyObject* coords = nullptr;

  API_SETUP_ARGS(G, self, args, "OsO|i", &self, &str1, &coords, &state);
  pymol::Result<> result;
  if (!str1[0]) {
    return APIFailure(G, "selection is empty");
  }
  API_ASSERT(APIEnterBlockedNotModal(G));

  if (SelectorGetTmp(G, str1, s1) >= 0) {
    int sele1 = SelectorIndexByName(G, s1);
    if (sele1 >= 0) {
      result = SelectorLoadCoords(G, coords, sele1, state);
    }
    SelectorFreeTmp(G, s1);
  }

  APIExitBlocked(G);

  if (!result && PyErr_Occurred()) {
    return nullptr;
  }

  return APIResult(G, result);
}

static PyObject* CmdLoadCoordSet(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* oname;
  PyObject* model;
  int frame;
  int quiet;

  API_SETUP_ARGS(G, self, args, "OsOii", &self, &oname, &model, &frame, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveLoadCoordset(G, oname, model, frame, static_cast<bool>(quiet));
  APIExit(G);
  return APISuccess();
}

static PyObject* CmdLoad(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *fname, *oname;
  char* object_props = nullptr;
  char* atom_props = nullptr;
  char* plugin = nullptr;
  int frame;
  cLoadType_t type;
  int finish, discrete;
  int quiet;
  int multiplex;
  int zoom;
  Py_ssize_t bytes;
  int mimic;
  const char* contents;

  API_SETUP_ARGS(G, self, args, "Oszz#iiiiiii|zzzi", &self, &oname, &fname,
                 &contents, &bytes, &frame, &type, &finish, &discrete, &quiet,
                 &multiplex, &zoom, &plugin, &object_props, &atom_props,
                 &mimic);
  API_ASSERT(APIEnterNotModal(G));

  auto result = ExecutiveLoad(G, fname, contents, bytes, type, oname, frame,
                              zoom, discrete, finish, multiplex, quiet, plugin,
                              object_props, atom_props, mimic);

  OrthoRestorePrompt(G);
  APIExit(G);

  return APIResult(G, result);
}

static PyObject* CmdLoadTraj(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *fname, *oname;
  int frame, type;
  int interval, average, start, stop, max, image;
  char* str1;
  float shift[3];
  char* plugin = nullptr;
  int quiet = 0; /* TODO */
  API_SETUP_ARGS(G, self, args, "Ossiiiiiiisifffs", &self, &oname, &fname,
                 &frame, &type, &interval, &average, &start, &stop, &max, &str1,
                 &image, &shift[0], &shift[1], &shift[2], &plugin);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveLoadTraj(G, oname, fname, frame, type, interval, average, start,
                        stop, max, str1, image, shift, plugin, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAssignAtomTypes(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* sele;
  int state = -1;
  int ok = false;
  int format;
  int quiet;
  OrthoLineType s1;

  API_SETUP_ARGS(G, self, args, "Osiii", &self, &sele, &format, &state, &quiet);
  {
    APIEnterBlocked(G);
    ok = (SelectorGetTmp(G, sele, s1) >= 0);
    if (ok) {
      /* format : 1: mol/sybyl, 2: macromodel/mmd */
      ok = ExecutiveAssignAtomTypes(G, s1, format, state, quiet);
      SelectorFreeTmp(G, s1);
    }
    APIExitBlocked(G);
  }
  return APIResultOk(G, ok);
}

const PyMethodDef kIoMethods[] = {
    {"assign_atom_types", CmdAssignAtomTypes, METH_VARARGS},
    {"find_molfile_plugin", CmdFindMolfilePlugin, METH_VARARGS},
    {"get_seq_align_str", CmdGetSeqAlignStr, METH_VARARGS},
    {"get_session", CmdGetSession, METH_VARARGS},
    {"get_str", CmdGetStr, METH_VARARGS},
    {"get_ccp4str", CmdGetCCP4Str, METH_VARARGS},
    {"load", CmdLoad, METH_VARARGS},
    {"load_color_table", CmdLoadColorTable, METH_VARARGS},
    {"load_coords", CmdLoadCoords, METH_VARARGS},
    {"load_coordset", CmdLoadCoordSet, METH_VARARGS},
    {"load_traj", CmdLoadTraj, METH_VARARGS},
    {"png", CmdPNG, METH_VARARGS},
    {"set_session", CmdSetSession, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
