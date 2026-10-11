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

// Python bindings for creating and fitting objects and ramps
// (pymol/creating.py, fitting.py, colorramping.py, experimenting.py).

#include "pymol/executive/executive.h"
#include "pymol/executive/executive_python.h"
#include "pymol/python/p.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/app/test_pymol.h"
#include "pymol/base/feedback.h"
#include "pymol/python/cmd_internal.h"

static PyObject* CmdPseudoatom(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *object_name, *sele, *label;
  char *name, *resn, *resi, *chain, *segi, *elem;
  float vdw;
  int hetatm, color;
  float b, q;
  PyObject* pos;
  int state, mode, quiet;

  API_SETUP_ARGS(G, self, args, "OssssssssfiffsOiiii", &self, &object_name,
                 &sele, &name, &resn, &resi, &chain, &segi, &elem, &vdw,
                 &hetatm, &b, &q, &label, &pos, &color, &state, &mode, &quiet);
  float pos_array[3], *pos_ptr = nullptr;
  if (pos && PyTuple_Check(pos) && (PyTuple_Size(pos) == 3))
    if (PyArg_ParseTuple(pos, "fff", pos_array, pos_array + 1, pos_array + 2))
      pos_ptr = pos_array;

  API_ASSERT(APIEnterBlockedNotModal(G));
  auto pseudoatom_name = ExecutivePreparePseudoatomName(G, object_name);
  auto result = ExecutivePseudoatom(G, pseudoatom_name, sele, name, resn, resi,
                                    chain, segi, elem, vdw, hetatm, b, q, label,
                                    pos_ptr, color, state, mode, quiet);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdSetRawAlignment(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* alnname;
  const char* guidename;
  PyObject* raw;
  int state = 0, quiet = 1;

  API_SETUP_ARGS(G, self, args,
                 "sOsii"
                 "O",
                 &alnname, &raw, &guidename, &state, &quiet, &self);
  API_ASSERT(APIEnterNotModal(G));
  auto result =
      ExecutiveSetRawAlignment(G, alnname, raw, guidename, state, quiet);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdGetVolumeRamp(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* result = nullptr;
  int ok = false;
  char* objName;
  int state;
  ok = PyArg_ParseTuple(args, "Osi", &self, &objName, &state);

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterBlockedNotModal(G))) {
    result = ExecutiveGetVolumeRamp(G, objName, state);
    APIExitBlocked(G);
  }

  if (!result) {
    return APIFailure();
  } else
    return result;
}

static PyObject* CmdSetVolumeRamp(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* objName;
  PyObject* ramp_list;
  std::vector<float> float_array;
  int state;

  API_SETUP_ARGS(G, self, args, "OsOi", &self, &objName, &ramp_list, &state);

  if (!PyList_Check(ramp_list) ||
      !PConvFromPyObject(G, ramp_list, float_array)) {
    return APIFailure(G, pymol::make_error("Invalid color array"));
  }

  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result =
      ExecutiveSetVolumeRamp(G, objName, std::move(float_array), state);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdGetBondPrint(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  char* str1;
  int*** array = nullptr;
  PyObject* result = nullptr;
  int int1, int2;
  int dim[3];
  ok = PyArg_ParseTuple(args, "Osii", &self, &str1, &int1, &int2);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    array = ExecutiveGetBondPrint(G, str1, int1, int2, dim);
    APIExit(G);
    if (array) {
      result = PConv3DIntArrayTo3DPyList(array, dim);
      FreeP(array);
    }
  }
  return (APIAutoNone(result));
}

static PyObject* CmdCurveNew(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char* curveName;
  const char* curveType;
  API_SETUP_ARGS(G, self, args, "Oss", &self, &curveName, &curveType);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveCurveNew(G, curveName, curveType);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdRampNew(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  int ok = false;
  char* map;
  int state;
  char* sele;
  float beyond, within;
  float sigma, *range_vla = nullptr;
  float* color_vla = nullptr;
  int zero, quiet, calc_mode = 0;
  PyObject *range, *color;
  API_SETUP_ARGS(G, self, args, "OssOOisfffii", &self, &name, &map, &range,
                 &color, &state, &sele, &beyond, &within, &sigma, &zero,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));
  pymol::Result<> result = pymol::Error();
  ok = G != nullptr;
  if (PyList_Size(range) > 0) ok = PConvPyListToFloatVLA(range, &range_vla);

  if (ok) {
    if (PyList_Check(color)) {
      if (PyList_Size(color) > 0)
        ok = PConvPyList3ToFloatVLA(color, &color_vla);
    } else if (PyInt_Check(color)) {
      ok = PConvPyIntToInt(color, &calc_mode);
    }
  }
  if (ok) {
    result =
        ExecutiveRampNew(G, name, map, pymol::vla_take_ownership(range_vla),
                         pymol::vla_take_ownership(color_vla), state, sele,
                         beyond, within, sigma, zero, calc_mode, quiet);
  }
  APIExit(G);
  return APIResult(G, result);
}

/* Synthesize a map from possibly weighted coefficients */
/* See MapCalculate for the actual calculation */
static PyObject* CmdMapGenerate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;

  char *name, *reflection_file, *tempFile, *amplitudes, *phases, *weights,
      *space_group;
  const char* cResult = nullptr;
  int quiet, zoom;
  double reso_high, reso_low, cell[6];

  ok = PyArg_ParseTuple(args, "Ossssszddsddddddii", &self, &name,
                        &reflection_file, &tempFile, &amplitudes, &phases,
                        &weights, &reso_low, &reso_high, &space_group, &cell[0],
                        &cell[1], &cell[2], &cell[3], &cell[4], &cell[5],
                        &quiet, &zoom);

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }

  if (ok && (ok = APIEnterNotModal(G))) {
    if (ok) {
      PRINTFB(G, FB_CCmd, FB_Blather)
      " Cmd-Update: Start ExecutiveMapGenerate." ENDFB(G);

      cResult = (const char*)ExecutiveMapGenerate(
          G, name, reflection_file, tempFile, amplitudes, phases, weights,
          reso_low, reso_high, space_group, cell, quiet, zoom);

      PRINTFB(G, FB_CCmd, FB_Blather)
      " Cmd-Update: Finished ExecutiveMapGenerate." ENDFB(G);
    }
    APIExit(G);
  }

  return APIAutoNone(Py_BuildValue("s", cResult));
}

static PyObject* CmdMapNew(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* name;
  float minCorner[3], maxCorner[3];
  float grid;
  float buffer, floor, ceiling, resolution;
  int type;
  int state;
  int have_corners;
  int quiet, zoom;
  int normalize;
  char* selection;
  API_SETUP_ARGS(G, self, args, "Osifsf(ffffff)iiiiifff", &self, &name, &type,
                 &grid, &selection, &buffer, &minCorner[0], &minCorner[1],
                 &minCorner[2], &maxCorner[0], &maxCorner[1], &maxCorner[2],
                 &state, &have_corners, &quiet, &zoom, &normalize, &floor,
                 &ceiling, &resolution);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveMapNew(
      G, name, type, grid, selection, buffer, minCorner, maxCorner, state,
      have_corners, quiet, zoom, normalize, floor, ceiling, resolution);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdAlign(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str2, *str3, *mfile, *oname;
  OrthoLineType s2 = "", s3 = "";
  int ok = false;
  int quiet, cycles, max_skip;
  float cutoff, gap, extend, seq;
  int state1, state2;
  int max_gap, transform, reset, window;
  float radius, scale, base, coord, expect, ante;
  ExecutiveRMSInfo rms_info;

  ok = PyArg_ParseTuple(args, "Ossfiffissiiiiiiffffffif", &self, &str2, &str3,
                        &cutoff, &cycles, &gap, &extend, &max_gap, &oname,
                        &mfile, &state1, &state2, &quiet, &max_skip, &transform,
                        &reset, &seq, &radius, &scale, &base, &coord, &expect,
                        &window, &ante);

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    PRINTFD(G, FB_CCmd)
    "CmdAlign-DEBUG %s %s\n", str2, str3 ENDFD;
    if ((ok = APIEnterNotModal(G))) {
      ok = ((SelectorGetTmp(G, str2, s2) >= 0) &&
            (SelectorGetTmp(G, str3, s3) >= 0));
      if (ok) {
        ok = ExecutiveAlign(G, s2, s3, mfile, gap, extend, max_gap, max_skip,
                            cutoff, cycles, quiet, oname, state1, state2,
                            &rms_info, transform, reset, seq, radius, scale,
                            base, coord, expect, window, ante);
      }
      SelectorFreeTmp(G, s2);
      SelectorFreeTmp(G, s3);
      APIExit(G);
    }
  }
  if (ok) {
    return Py_BuildValue("(fiififi)", rms_info.final_rms, rms_info.final_n_atom,
                         rms_info.n_cycles_run, rms_info.initial_rms,
                         rms_info.initial_n_atom, rms_info.raw_alignment_score,
                         rms_info.n_residues_aligned);
  } else {
    return APIFailure();
  }
}

static PyObject* CmdUSalign(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  const char *mobile, *target, *oname;
  int mobile_state, target_state, quiet, transform, fast;
  API_SETUP_ARGS(G, self, args, "Ossiiiisi", &self, &mobile, &target,
                 &mobile_state, &target_state, &quiet, &transform, &oname,
                 &fast);
  API_ASSERT(APIEnterNotModal(G));

  OrthoLineType s1, s2;
  int ok = (SelectorGetTmp(G, mobile, s1) >= 0) &&
           (SelectorGetTmp(G, target, s2) >= 0);

  // Store results in locals — can't call Py_BuildValue until after APIExit
  double tm_target = 0, tm_mobile = 0, rmsd = 0, seq_id = 0;
  int ali_len = 0;
  bool have_result = false;
  std::string err_msg;

  if (ok) {
    auto res = ExecutiveUSalign(G, s1, s2, mobile_state, target_state, quiet,
                                transform, oname, fast);
    if (res) {
      auto& r = res.result();
      tm_target = r.tm_score_target;
      tm_mobile = r.tm_score_mobile;
      rmsd = r.rmsd;
      ali_len = r.aligned_length;
      seq_id = r.seq_identity;
      have_result = true;
    } else {
      err_msg = res.error().what();
    }
  }

  SelectorFreeTmp(G, s1);
  SelectorFreeTmp(G, s2);
  APIExit(G);

  if (have_result) {
    return Py_BuildValue("{s:d,s:d,s:d,s:i,s:d}", "tm_score_target", tm_target,
                         "tm_score_mobile", tm_mobile, "RMSD", rmsd,
                         "alignment_length", ali_len, "seq_identity", seq_id);
  }
  if (!err_msg.empty()) {
    PyErr_SetString(PyExc_RuntimeError, err_msg.c_str());
    return nullptr;
  }
  return APIFailure();
}

static PyObject* CmdDump(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int state, quiet;
  int ok = false;
  ok = PyArg_ParseTuple(args, "Ossii", &self, &str1, &str2, &state, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    ExecutiveDump(G, str1, str2, state, quiet); /* TODO STATUS */
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdIsomesh(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *mesh_name, *map_name, *sele;
  float lvl, fbuf, alt_lvl;
  int mesh_mode;
  int state = -1;
  float carve;
  int map_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssffiifiif", &self, &mesh_name, &map_name,
                 &sele, &fbuf, &lvl, &mesh_mode, &state, &carve, &map_state,
                 &quiet, &alt_lvl);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveIsomeshEtc(G, mesh_name, map_name, lvl, sele, fbuf, state,
                                 carve, map_state, quiet, mesh_mode, alt_lvl);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdSliceNew(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* slice;
  char* map;
  int state, map_state;

  API_SETUP_ARGS(G, self, args, "Ossii", &self, &slice, &map, &state,
                 &map_state);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSliceNew(G, slice, map, state, map_state);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdIsosurface(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *surf_name, *map_name, *sele;
  float lvl, fbuf;
  int surf_mode;
  int state = -1;
  float carve;
  int map_state = 0;
  int side;
  int quiet;

  API_SETUP_ARGS(G, self, args, "Osssffiifiii", &self, &surf_name, &map_name,
                 &sele, &fbuf, &lvl, &surf_mode, &state, &carve, &map_state,
                 &side, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res =
      ExecutiveIsosurfaceEtc(G, surf_name, map_name, lvl, sele, fbuf, state,
                             carve, map_state, side, quiet, surf_mode);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdSymExp(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2, *str3;
  float cutoff;
  pymol::CObject* mObj;
  int segi;
  int quiet;
  /* oper 0 = all, 1 = sele + buffer, 2 = vector */

  int ok = false;
  ok = PyArg_ParseTuple(args, "Osssfii", &self, &str1, &str2, &str3, &cutoff,
                        &segi, &quiet);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    mObj = ExecutiveFindObjectByName(G, str2);
    if (mObj) {
      if (mObj->type != cObjectMolecule) {
        mObj = nullptr;
        ok = false;
      }
    }
    if (mObj) {
      ExecutiveSymExp(G, str1, str2, str3, cutoff, segi,
                      quiet); /* TODO STATUS */
    }
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdCopy(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int zoom;
  API_SETUP_ARGS(G, self, args, "Ossi", &self, &str1, &str2, &zoom);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveCopy(G, str1, str2, zoom);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdMem(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  ok = PyArg_ParseTuple(args, "O", &self);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok) {
    OVHeap_Dump(G->Context->heap, 0);
    SelectorMemoryDump(G);
    ExecutiveMemoryDump(G);
  }
  return APISuccess();
}

static PyObject* CmdCreate(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int target, source, discrete, quiet;
  int singletons;
  int copy_properties = 0;
  int zoom;
  API_SETUP_ARGS(G, self, args, "Ossiiiiiii", &self, &str1, &str2, &source,
                 &target, &discrete, &zoom, &quiet, &singletons,
                 &copy_properties);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveSeleToObject(G, str1, str2, source, target, discrete,
                                   zoom, quiet, singletons, copy_properties);
  APIExit(G);
  return APIResult(G, res);
}

static PyObject* CmdFitPairs(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  PyObject* list;
  int quiet = 0;
  API_SETUP_ARGS(G, self, args, "OOi", &self, &list, &quiet);
  API_ASSERT(APIEnterBlockedNotModal(G));
  auto result = ExecutiveFitPairs(G, list, quiet);
  APIExitBlocked(G);
  return APIResult(G, result);
}

static PyObject* CmdIntraFit(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char* str1;
  int state;
  int mode;
  int quiet;
  int mix;
  int pbc = true;
  API_SETUP_ARGS(G, self, args, "Osiiii|i", &self, &str1, &state, &mode, &quiet,
                 &mix, &pbc);
  API_ASSERT(APIEnterNotModal(G));
  auto fVLA = ExecutiveRMSStates(G, str1, state, mode, quiet, mix, pbc);
  APIExit(G);
  PyObject* result = nullptr;
  if (fVLA) {
    result = PConvFloatVLAToPyList(fVLA.result().data());
  }
  return APIAutoNone(result);
}

static PyObject* CmdFit(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *str1, *str2;
  int mode;
  int quiet;
  float cutoff;
  int state1, state2;
  int matchmaker, cycles;
  char* object;
  API_SETUP_ARGS(G, self, args, "Ossiiiiifis", &self, &str1, &str2, &mode,
                 &state1, &state2, &quiet, &matchmaker, &cutoff, &cycles,
                 &object);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveFit(G, str1, str2, mode, cutoff, cycles, quiet, object,
                             state1, state2, matchmaker);
  APIExit(G);
  if (result) {
    return PConvToPyObject(result.result().final_rms);
  } else {
    return APIFailure(G, result.error());
  }
}

static PyObject* CmdGroup(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *gname, *names;
  int quiet, action;
  int ok = false;
  API_SETUP_ARGS(G, self, args, "Ossii", &self, &gname, &names, &action,
                 &quiet);
  API_ASSERT(APIEnterNotModal(G));
  {
    ok = ExecutiveGroup(G, gname, names, action, quiet);
    APIExit(G);
  }
  return APIResultOk(G, ok);
}

static PyObject* CmdSpheroid(PyObject* self, PyObject* args)

/* EXPERIMENTAL */
{
  PyMOLGlobals* G = nullptr;
  char* name;
  int average;
  API_SETUP_ARGS(G, self, args, "Osi", &self, &name, &average);
  API_ASSERT(APIEnterNotModal(G));
  auto result = ExecutiveSpheroid(G, name, average);
  APIExit(G);
  return APIResult(G, result);
}

static PyObject* CmdTest(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  /* regression tests */

  int ok = false;
  int code;
  int group;

  ok = PyArg_ParseTuple(args, "Oii", &self, &group, &code);
  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }
  if (ok && (ok = APIEnterNotModal(G))) {
    PRINTFB(G, FB_CCmd, FB_Details)
    " Cmd: initiating test %d-%d.\n", group, code ENDFB(G);
    ok = TestPyMOLRun(G, group, code);
    PRINTFB(G, FB_CCmd, FB_Details)
    " Cmd: concluding test %d-%d.\n", group, code ENDFB(G);
    APIExit(G);
  }
  return APIResultOk(ok);
}

static PyObject* CmdIsolevel(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  float level;
  int state;
  char* name;
  int query, quiet;
  API_SETUP_ARGS(G, self, args, "Osfiii", &self, &name, &level, &state, &query,
                 &quiet);
  if (!query) {
    API_ASSERT(APIEnterNotModal(G));
    auto result = ExecutiveIsolevel(G, name, level, state, quiet);
    APIExit(G);
    return APIResult(G, result);
  } else {
    APIEnter(G);
    auto result = ExecutiveGetIsolevel(G, name, state);
    APIExit(G);
    return APIResult(G, result);
  }
}

static PyObject* CmdCEAlign(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  int ok = false;
  int windowSize = 8, gap_max = 30;
  float d0 = 3.0, d1 = 4.0;
  PyObject *listA, *listB, *result;
  Py_ssize_t lenA, lenB;

  /* Unpack the arguments from Python */

  ok = PyArg_ParseTuple(args, "OOO|ffii", &self, &listA, &listB, &d0, &d1,
                        &windowSize, &gap_max);

  /* Handle errors */

  if (ok) {
    API_SETUP_PYMOL_GLOBALS;
    ok = (G != nullptr);
  } else {
    API_HANDLE_ERROR;
  }

  /* Get the list lengths */

  lenA = PyList_Size(listA);
  if (lenA < 1) {
    result = nullptr;
    ok = false;
  }

  if (ok) lenB = PyList_Size(listB);
  if (ok && lenB < 1) {
    result = nullptr;
    ok = false;
  }

  /* Call CEAlign */

  if (ok) {
    APIEnterBlocked(G);
    result = (PyObject*)ExecutiveCEAlign(G, listA, listB, lenA, lenB, d0, d1,
                                         windowSize, gap_max);
    APIExitBlocked(G);
  }
  return result;
}

static PyObject* CmdVolume(PyObject* self, PyObject* args) {
  PyMOLGlobals* G = nullptr;
  char *volume_name, *map_name, *sele;
  float lvl, fbuf;
  int state = -1;
  float carve;
  int map_state;
  int quiet;
  API_SETUP_ARGS(G, self, args, "Osssffifii", &self, &volume_name, &map_name,
                 &sele, &fbuf, &lvl, &state, &carve, &map_state, &quiet);
  API_ASSERT(APIEnterNotModal(G));
  auto res = ExecutiveVolume(G, volume_name, map_name, lvl, sele, fbuf, state,
                             carve, map_state, quiet);
  APIExit(G);
  return APIResult(G, res);
}

const PyMethodDef kCreatingMethods[] = {
    {"align", CmdAlign, METH_VARARGS},
    {"cealign", CmdCEAlign, METH_VARARGS},
    {"copy", CmdCopy, METH_VARARGS},
    {"create", CmdCreate, METH_VARARGS},
    {"curve_new", CmdCurveNew, METH_VARARGS},
    {"dump", CmdDump, METH_VARARGS},
    {"fit", CmdFit, METH_VARARGS},
    {"fit_pairs", CmdFitPairs, METH_VARARGS},
    {"get_bond_print", CmdGetBondPrint, METH_VARARGS},
    {"get_volume_ramp", CmdGetVolumeRamp, METH_VARARGS},
    {"set_volume_ramp", CmdSetVolumeRamp, METH_VARARGS},
    {"group", CmdGroup, METH_VARARGS},
    {"intrafit", CmdIntraFit, METH_VARARGS},
    {"isolevel", CmdIsolevel, METH_VARARGS},
    {"isomesh", CmdIsomesh, METH_VARARGS},
    {"isosurface", CmdIsosurface, METH_VARARGS},
    {"map_generate", CmdMapGenerate, METH_VARARGS},
    {"map_new", CmdMapNew, METH_VARARGS},
    {"mem", CmdMem, METH_VARARGS},
    {"pseudoatom", CmdPseudoatom, METH_VARARGS},
    {"ramp_new", CmdRampNew, METH_VARARGS},
    {"set_raw_alignment", CmdSetRawAlignment, METH_VARARGS},
    {"slice_new", CmdSliceNew, METH_VARARGS},
    {"spheroid", CmdSpheroid, METH_VARARGS},
    {"symexp", CmdSymExp, METH_VARARGS},
    {"test", CmdTest, METH_VARARGS},
    {"volume", CmdVolume, METH_VARARGS},
    {"usalign", CmdUSalign, METH_VARARGS},
    {nullptr, nullptr} /* sentinel */
};
