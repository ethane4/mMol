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

// Selecting, iterating and labeling atoms.

#include "pymol/ui/but_mode.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_volume.h"
#include "pymol/pyconv/pconv.h"
#include "pymol/ui/seeker.h"
#include "pymol/ui/wizard.h"
#include "pymol/executive/executive_internal.h"

// GL rectangle picking is gone, so a box never picks atoms; a "set" box still
// clears the active selection.
pymol::Result<NetSelect> ExecutiveSelectRect(PyMOLGlobals* G, BlockRect* rect,
                                             int mode) {
  if (mode == cButModeSeleSetBox) {
    ObjectNameType name;
    if (ExecutiveGetActiveSeleName(G, name, false,
                                   SettingGetGlobal_i(G, cSetting_logging))) {
      ExecutiveSetObjVisib(G, name, 0, false);
      if (SettingGetGlobal_i(G, cSetting_logging)) {
        auto buf2 = pymol::string_format("cmd.disable('%s')\n", name);
        PLog(G, buf2, cPLog_no_flush);
      }
    }
  }
  return NetSelect{};
}
pymol::Result<> ExecutiveLabel(PyMOLGlobals* G, const char* str1,
                               const char* expr, int quiet, int eval_mode) {
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  SelectorTmp s1(G, str1);
  sele1 = s1.getIndex();
  if (sele1 >= 0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_LABL;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = eval_mode;

    pymol::pautoblock gil(G);

    if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
      return pymol::Error();
    }

    cnt = op1.i1;
    op1.code = OMOP_VISI;
    op1.i1 = cRepLabelBit;
    op1.i2 = cVis_SHOW;
    ExecutiveObjMolSeleOp(G, sele1, &op1);
    op1.code = OMOP_INVA;
    op1.i2 = cRepInvVisib;
    ExecutiveObjMolSeleOp(G, sele1, &op1);

    if (!quiet) {
      {
        const char* unlabelledstr = "";
        if (cnt < 0) { /* if negative, say unlabelled */
          cnt = -cnt;
          unlabelledstr = "un";
        }
        PRINTFB(G, FB_Executive, FB_Actions)
        " Label: %slabelled %i atoms.\n", unlabelledstr, cnt ENDFB(G);
      }
    }
  } else {
    return pymol::make_error("No atoms selected");
  }
  return {};
}
pymol::Result<int> ExecutiveIterate(PyMOLGlobals* G, const char* str1,
                                    const char* expr, int read_only, int quiet,
                                    PyObject* space) {
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRecInit(&op1);
  SelectorTmp tmpsele1(G, str1);
  int sele1 = tmpsele1.getIndex();
  op1.i1 = 0;
  if (sele1 >= 0) {
    op1.code = OMOP_ALTR;
    op1.i1 = 0;
    op1.i2 = read_only;
    op1.s1 = expr;
    op1.py_ob1 = space;

    if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
      return pymol::Error();
    }

    if (!quiet) {
      if (!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
        " Alter: modified %i atoms.\n", op1.i1 ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
        " Iterate: iterated over %i atoms.\n", op1.i1 ENDFB(G);
      }
    }
    if (!read_only) {
      SeqChanged(G);
    }
  } else {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " %s: No atoms selected.\n", __func__ ENDFB(G);
    }
  }
  return (op1.i1);
}

SelectArgs ExecutiveSelectPrepareArgs(PyMOLGlobals* G,
                                      pymol::zstring_view sname,
                                      pymol::zstring_view sele) {
  SelectArgs args;
  args.sname = sname.c_str();
  args.sele = sele.c_str();

  // selection in first argument (cmd.select("expr"))
  if (args.sele.empty()) {
    args.sele = sname.c_str();
    args.sname =
        SettingGet<bool>(G, cSetting_auto_number_selections) ? "" : "sele";
  }

  if (args.sname.empty()) {
    auto sel_num = SettingGet<int>(G, cSetting_sel_counter) + 1;
    SettingSet<int>(G, cSetting_sel_counter, sel_num);
    args.sname =
        pymol::string_format("sel%02u", static_cast<unsigned>(sel_num));
  }
  return args;
}

/**
 * cmd.select() implementation
 *
 * @param name Name of selection to create, or selection expression if `sele` is
 * empty. If `sele` is empty, `name` will become "sele" or "sel01", depending on
 * `auto_number_selections` and `sel_counter`. If `name` is empty, it will
 * become "sel01", depending on `sel_counter` but independent of
 * `auto_number_selections`.
 * @param sele Selection expression
 * @param enable Enable the named selection's SpecRec if 1, disable if 0, keep
 * unchanged (if exists, otherwise enable) if -1.
 * @param merge Discard existing named selection if 0. Merge with existing named
 * selection if 1, merge only if existing one is enabled if 2.
 * @param state Object state for state (coordinate) dependent expressions
 * @param domain Existing selection name, same as selecting `(sele) and
 * (domain)`
 * @return Number of selected atoms
 */
pymol::Result<int> ExecutiveSelect(PyMOLGlobals* G, const SelectArgs& sargs,
                                   int enable, int quiet, int merge, int state,
                                   const char* domain) {
  const char* name = sargs.sname.c_str();
  const char* sele = sargs.sele.c_str();
  // bail if name not available
  if (ExecutiveFindObjectByName(G, name)) {
    return pymol::make_error("name conflicts with an object");
  }

  // merge with existing selection
  std::string selebuf;
  if (merge) {
    if (merge == 2) {
      // merge if exists and active
      selebuf = pymol::join_to_string("(", sele, ") or ??", name);
    } else {
      // merge if exists
      selebuf = pymol::join_to_string("(", sele, ") or ?", name);
    }
    sele = selebuf.c_str();
  }

  auto res = SelectorCreateWithStateDomain(G, name, sele, nullptr, quiet,
                                           nullptr, state, domain);

  p_return_if_error(res);

  if (enable == 1) {
    ExecutiveSetObjVisib(G, name, 1, 0);
  } else if (enable == 0) {
    ExecutiveSetObjVisib(G, name, 0, 0);
  }

  SceneInvalidate(G);
  SeqDirty(G);

  return res.result();
}
pymol::Result<int> ExecutiveSelectList(PyMOLGlobals* G, const char* sele_name,
                                       const char* oname, const int* list,
                                       size_t list_len, int state, int mode,
                                       int quiet) {
  enum {
    MODE_INDEX = 0,
    MODE_ID = 1,
    MODE_RANK = 2,
  };

  auto obj = ExecutiveFindObject<ObjectMolecule>(G, oname);
  if (!obj) {
    return pymol::make_error("object not found");
  }

  std::vector<int> idx_list;
  idx_list.reserve(list_len);

  if (mode == MODE_INDEX) {
    for (size_t i = 0; i != list_len; ++i) {
      // convert 1-based indices to 0-based array offsets
      idx_list.push_back(list[i] - 1);
    }
  } else if (mode == MODE_ID || mode == MODE_RANK) {
    CoordSet const* cs = obj->getCoordSet(state);
    std::set<int> const list_set(list, list + list_len);

    for (int atm = 0; atm < obj->NAtom; ++atm) {
      auto const& ai = obj->AtomInfo[atm];
      int const index = (mode == MODE_ID) ? ai.id : ai.rank;

      if (list_set.count(index)) {
        if (!cs || cs->atmToIdx(atm) >= 0) {
          idx_list.push_back(atm);
        }
      }
    }
  } else {
    return pymol::make_error("invalid mode");
  }

  return SelectorCreateOrderedFromObjectIndices(
      G, sele_name, obj, idx_list.data(), idx_list.size());
}
pymol::Result<int> ExecutiveIterateList(PyMOLGlobals* G, const char* str1,
                                        PyObject* list, int read_only,
                                        int quiet, PyObject* space) {
  assert(PyGILState_Check());

  int ok = true;
  int n_eval = 0;
  SelectorTmp s1(G, str1);
  int sele0 = s1.getIndex();
  PyObject* entry = nullptr;
  ObjectMolecule* obj = nullptr;
  if (sele0 >= 0) obj = SelectorGetSingleObjectMolecule(G, sele0);
  if (obj) {
    int n_atom = obj->NAtom;
    int list_len = 0;
    int a;
    int index = 0;
    const char* expr = nullptr;
    if (ok) ok = PyList_Check(list);
    if (ok) {
      list_len = PyList_Size(list);
      for (a = 0; a < list_len; a++) {
        if (ok) entry = PyList_GetItem(list, a);
        if (ok) ok = PyList_Check(entry);
        if (ok) ok = (PyList_Size(entry) == 2);
        if (ok) ok = PConvPyIntToInt(PyList_GetItem(entry, 0), &index);
        if (ok) ok = PConvPyStrToStrPtr(PyList_GetItem(entry, 1), &expr);
        if (ok) ok = ((index <= n_atom) && (index > 0));
        if (ok) {
          CoordSet* cs = nullptr;
          if (obj->DiscreteFlag && obj->DiscreteCSet) {
            cs = obj->DiscreteCSet[index - 1];
          } else if (obj->NCSet == 1) {
            cs = obj->CSet[0];
          }
          auto expr_co =
              unique_PyObject_ptr(Py_CompileString(expr, "", Py_single_input));
          ok = expr_co && PAlterAtom(G, obj, cs, expr_co.get(), read_only,
                                     index - 1, space);
        }
        if (ok)
          n_eval++;
        else
          break;
      }
    }
  } else {
    return pymol::make_error("Selection cannot span more than one object.");
  }
  if (ok) {
    if (!quiet) {
      if (!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
        " AlterList: modified %i atoms.\n", n_eval ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
        " IterateList: iterated over %i atoms.\n", n_eval ENDFB(G);
      }
    }
    if (!read_only) {
      SeqChanged(G);
    }
  }
  if (!ok)
    return pymol::make_error("An error occurred.");
  else
    return n_eval;
}
pymol::Result<int> ExecutiveIterateState(PyMOLGlobals* G, int state,
                                         const char* str1, const char* expr,
                                         int read_only, int quiet,
                                         PyObject* space) {
  SelectorTmp tmpsele1(G, str1);
  int sele1 = tmpsele1.getIndex();
  if (sele1 >= 0) {
    int start_state = 0, stop_state = 0;
    ObjectMoleculeOpRec op1;

    if (state >= 0) {
      start_state = state;
      stop_state = state + 1;
    } else {
      if ((state == -2) ||
          (state == -3)) { /* current state, TO DO: effective object state */
        state = SceneGetState(G);
        start_state = state;
        stop_state = state + 1;
      } else if (state == -1) { /* all states (for the selection) */
        start_state = 0;
        stop_state = SelectorCountStates(G, sele1);
      }
    }
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;

    for (state = start_state; state < stop_state; state++) {
      op1.code = OMOP_AlterState;
      op1.s1 = expr;
      op1.py_ob1 = space;
      op1.i2 = state;
      op1.i3 = read_only;
      if (!ExecutiveObjMolSeleOp(G, sele1, &op1)) {
        return pymol::Error();
      }
    }
    if (!read_only) {
      // for dynamic_measures
      ExecutiveUpdateCoordDepends(G, nullptr);
      SeqChanged(G);
    }
    if (!quiet) {
      if (!read_only) {
        PRINTFB(G, FB_Executive, FB_Actions)
        " AlterState: modified %i atom coordinate states.\n", op1.i1 ENDFB(G);
      } else {
        PRINTFB(G, FB_Executive, FB_Actions)
        " IterateState: iterated over %i atom coordinate states.\n",
            op1.i1 ENDFB(G);
      }
    }
    return op1.i1;
  } else {
    if (!quiet) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      "ExecutiveIterateState: No atoms selected.\n" ENDFB(G);
    }
    return 0;
  }
}
