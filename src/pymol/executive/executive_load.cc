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

// Loading molecules, maps and trajectories from files and strings.

#include "pymol/io/cif_file.h"
#include "pymol/executive/editor.h"
#include "pymol/executive/executive.h"
#include "pymol/base/feedback.h"
#include "pymol/objects/object_cgo.h"
#include "pymol/objects/object_volume.h"
#include "pymol/io/plug_io_manager.h"
#include "ov/ov_context.h"
#include "pymol/executive/executive_load.h"
#include "pymol/base/file_stream.h"
#include "pymol/executive/executive_internal.h"

/* routines that still need to be updated for Tracker list iteration

Low priority:

ExecutiveSculptIterate
ExecutiveSculptActivate
ExecutiveSculptDeactivate
ExecutiveRenameObjectAtoms
ExecutiveSpheroid

*/

/**
 * @brief Processes CIF string and manages either Map or Molecule object
 * @param I Pointer to a possibly existing object
 * @param object_name Name of the object
 * @param st CIF string
 * @param frame State/Frame number
 * @param discrete Discrete object flag
 * @param quiet Suppress output
 * @param multiplex Multiplex flag
 * @param zoom Zoom flag
 */
static pymol::Result<pymol::CObject*> ExecutiveProcessCif(
    PyMOLGlobals* G, pymol::CObject* I, const char* object_name, const char* st,
    int frame, int discrete, int quiet, int multiplex, int zoom);

/**
 * Load any file type which is implemented in C.
 *
 * fname:       File name, can be empty if content is provided
 * content:     File contents, if not nullptr
 * content_length:      Length of "content", if not nullptr
 * content_format:      File type code
 * object_name:         New object name
 * state:       Object state to start loading new coordsets in
 * zoom:        Zoom on new loaded atoms
 * discrete:    Make discrete states
 * finish:      update object selection & zoom
 * multiplex:   Split new states into objects
 * quiet:       Suppress feedback
 * object_props:        names of object properties to load
 * atom_props:          names of atom properties to load
 */
pymol::Result<> ExecutiveLoad(PyMOLGlobals* G, const char* fname,
                              const char* content, int content_length,
                              cLoadType_t content_format,
                              const char* object_name_proposed, int state,
                              int zoom, int discrete, int finish, int multiplex,
                              int quiet, const char* plugin_arg,
                              const char* object_props, const char* atom_props,
                              bool mimic) {
  auto res = ExecutiveLoadPrepareArgs(
      G, fname, content, content_length, content_format, object_name_proposed,
      state, zoom, discrete, finish, multiplex, quiet, plugin_arg, object_props,
      atom_props, mimic);
  p_return_if_error(res);
  return ExecutiveLoad(G, res.result());
}

/**
 * Prepare arguments for ExecutiveLoad
 */
pymol::Result<ExecutiveLoadArgs> ExecutiveLoadPrepareArgs(
    PyMOLGlobals* G, pymol::null_safe_zstring_view fname, const char* content,
    int content_length, cLoadType_t content_format,
    const char* object_name_proposed, int state, int zoom, int discrete,
    int finish, int multiplex, int quiet, const char* plugin_arg,
    const char* object_props, const char* atom_props, bool mimic) {
  ExecutiveLoadArgs args;
  bool fname_null_ok = false;

  // validate proposed object name
  ObjectNameType object_name = "";
  ExecutiveProcessObjectName(G, object_name_proposed, object_name);
  args.object_name = object_name;

  // Ensure correct float parsing with scanf. It's possible to change this from
  // Python, so don't rely on a persistent global value.
  std::setlocale(LC_NUMERIC, "C");

  if (!object_props)
    object_props = SettingGetGlobal_s(G, cSetting_load_object_props_default);
  if (!atom_props)
    atom_props = SettingGetGlobal_s(G, cSetting_load_atom_props_default);

  // multiplex -2 -> "multiplex" setting
  // multiplex -1 -> file type dependant default
  // multiplex 1 -> split entries (don't try to load into existing object)
  if (multiplex == -2) {
    multiplex = SettingGetGlobal_i(G, cSetting_multiplex);
  }

  switch (content_format) {
    // string loading functions
    case cLoadTypePDBStr:
    case cLoadTypeVDBStr:
    case cLoadTypeCIFStr:
    case cLoadTypeMAEStr:
    case cLoadTypeXPLORStr:
    case cLoadTypeCCP4Str:
    case cLoadTypeCCP4UnspecifiedStr:
    case cLoadTypeMRCStr:
    case cLoadTypePHIStr:
    case cLoadTypeMMDStr:
    case cLoadTypeMOLStr:
    case cLoadTypeMOL2Str:
    case cLoadTypeSDF2Str:
    case cLoadTypeXYZStr:
    case cLoadTypeDXStr:
    case cLoadTypeBCIFStr:
      if (!content) {
        return pymol::make_error("content is nullptr");
      }
      fname_null_ok = true;
      break;
    case cLoadTypePQR:
    case cLoadTypePDBQT:
    case cLoadTypePDB:
    case cLoadTypeCIF:
    case cLoadTypeMAE:
    case cLoadTypeXPLORMap:
    case cLoadTypeCCP4Map:
    case cLoadTypeCCP4Unspecified:
    case cLoadTypeMRC:
    case cLoadTypePHIMap:
    case cLoadTypeMMD:
    case cLoadTypeMOL:
    case cLoadTypeMOL2:
    case cLoadTypeSDF2:
    case cLoadTypeXYZ:
    case cLoadTypeDXMap:
    case cLoadTypeBCIF:
      if (content) {
        fname_null_ok = true;
        break;
      }

      if (fname.empty()) {
        break;
      }

      try {
        args.content = pymol::file_get_contents(fname);
        PRINTFB(G, FB_Executive, FB_Blather)
        " %s: Loading from %s.\n", __func__, fname.c_str() ENDFB(G);
      } catch (...) {
        return pymol::make_error(
            pymol::string_format("Unable to open file '%s'", fname.c_str()));
      }

      break;

    // molfile_plugin based formats
    case cLoadTypeCUBEMap:
      args.plugin = "cube";
      break;
    case cLoadTypeSpider:
      args.plugin = "spider";
      break;
    case cLoadTypeXTC:
      args.plugin = "xtc";
      break;
    case cLoadTypeTRR:
      args.plugin = "trr";
      break;
    case cLoadTypeGRO:
      args.plugin = "gro";
      break;
    case cLoadTypeG96:
      args.plugin = "g96";
      break;
    case cLoadTypeTRJ2:
      args.plugin = "trj";
      break;
    case cLoadTypeDCD:
      args.plugin = "dcd";
      break;
    case cLoadTypeDTR:
      args.plugin = "dtr";
      break;
    case cLoadTypeCMS:
      args.plugin = "mae";
      break;
    default:
      if (plugin_arg) {
        args.plugin = plugin_arg;
        // Hackish way to request a "sub-plugin". For example the "cube" plugin
        // can load maps or molecules:
        // load example.cube, molobj, format=plugin:cube:1
        // load example.cube, mapobj, format=plugin:cube:4
        auto pos = args.plugin.find(':');
        if (pos != std::string::npos) {
          args.plugin_mask = atoi(args.plugin.c_str() + pos + 1);
          args.plugin.resize(pos);
        }
      }
      break;
  }

  if (fname.empty() && !fname_null_ok) {
    return pymol::make_error("This format requires a filename to load");
  }

  if (!args.plugin.empty()) {
    content_format = cLoadTypePlugin;
  }

  if (content) {
    assert(args.content.empty());
    args.content = std::string(content, content_length);
  }

  args.content_format = content_format;
  args.fname = fname.c_str();
  args.state = state;
  args.zoom = zoom;
  args.discrete = discrete;
  args.finish = finish;
  args.multiplex = multiplex;
  args.quiet = quiet;
  args.object_props = object_props;
  args.atom_props = atom_props;
  args.mimic = mimic;
  return args;
}

pymol::Result<> ExecutiveLoad(PyMOLGlobals* G, ExecutiveLoadArgs const& args) {
  pymol::CObject* origObj = nullptr;
  const char* fname = args.fname.c_str();
  const char* content = args.content.data();
  int size = args.content.size();
  const char* object_name = args.object_name.c_str();
  const char* object_props = args.object_props.c_str();
  const char* atom_props = args.atom_props.c_str();
  const char* plugin = args.plugin.c_str();
  auto content_format = args.content_format;
  auto state = args.state;
  auto zoom = args.zoom;
  auto discrete = args.discrete;
  auto finish = args.finish;
  auto multiplex = args.multiplex;
  auto quiet = args.quiet;

  if (multiplex != 1) {
    origObj = ExecutiveGetExistingCompatible(G, object_name, content_format);
  }

  // file type dependent multiplex and discrete default
  if (discrete < 0) {
    if (multiplex == 1) {
      discrete = 0;
    } else {
      switch (content_format) {
        case cLoadTypeMOL2:
        case cLoadTypeMOL2Str:
          discrete = -1;    /* content-dependent behavior... */
        case cLoadTypeSDF2: /* SDF files currently default to discrete */
        case cLoadTypeSDF2Str:
        case cLoadTypeMAE: /* Let MAE handle itself */
        case cLoadTypeMAEStr:
          break;
        default:
          discrete = 0;
          break;
      }
    }
  }

  int ok = true;
  OrthoLineType buf = "";
  int pdb_variant = PDB_VARIANT_DEFAULT;
  pymol::CObject* obj = nullptr;

  // downstream file type reading functions
  switch (content_format) {
    case cLoadTypePQR:
      pdb_variant = PDB_VARIANT_PQR;
    case cLoadTypePDBQT:
      if (content_format == cLoadTypePDBQT) pdb_variant = PDB_VARIANT_PDBQT;
    case cLoadTypeVDBStr:
      if (ok && content_format == cLoadTypeVDBStr)
        pdb_variant = PDB_VARIANT_VDB;
    case cLoadTypePDB:
    case cLoadTypePDBStr:
      ok = ExecutiveProcessPDBFile(G, origObj, fname, content, object_name,
                                   state, discrete, finish, buf, pdb_variant,
                                   quiet, multiplex, zoom);
      break;
    case cLoadTypeCIF:
    case cLoadTypeCIFStr: {
      auto res = ExecutiveProcessCif(G, static_cast<ObjectMolecule*>(origObj),
                                     object_name, content, state, discrete,
                                     quiet, multiplex, zoom);
      p_return_if_error(res);
      obj = res.result();
    } break;
    case cLoadTypeBCIF:
    case cLoadTypeBCIFStr: {
      auto res = ObjectMoleculeReadBCif(
          G, static_cast<ObjectMolecule*>(origObj), content, size, state,
          discrete, quiet, multiplex, zoom);
      p_return_if_error(res);
      obj = res.result();
    } break;
    case cLoadTypeMAE:
    case cLoadTypeMAEStr:
      return pymol::Error::make<pymol::Error::INCENTIVE_ONLY>(
          "'mae' format not supported by this PyMOL build");
      break;
    case cLoadTypeTOP:
      if (origObj) {
        /* always reinitialize topology objects from scratch */
        ExecutiveDelete(G, origObj->Name);
        origObj = nullptr;
      }
      obj = ObjectMoleculeLoadTOPFile(G, nullptr, fname, state, discrete);
      break;
    case cLoadTypeTRJ:
      if (origObj) {
        ObjectMoleculeLoadTRJFile(G, (ObjectMolecule*)origObj, fname, state, 1,
                                  1, 1, -1, -1, nullptr, 1, nullptr, quiet);
      } else {
        return pymol::make_error(
            "must load object topology before loading trajectory!");
      }
      break;
    case cLoadTypeCRD:
      if (origObj) {
        ObjectMoleculeLoadRSTFile(G, (ObjectMolecule*)origObj, fname, state,
                                  quiet, 1);
      } else {
        return pymol::make_error(
            "must load object topology before loading coordinate file!");
      }
      break;
    case cLoadTypeRST:
      if (origObj) {
        ObjectMoleculeLoadRSTFile(G, (ObjectMolecule*)origObj, fname, state,
                                  quiet, 0);
      } else {
        return pymol::make_error(
            "must load object topology before loading restart file!");
      }
      break;
    case cLoadTypePMO:
      return pymol::make_error("PMO format no longer supported.");
    case cLoadTypeDXStr:
      obj = ObjectMapReadDXStr(G, dynamic_cast<ObjectMap*>(origObj), content,
                               size, state, quiet);
      break;
    case cLoadTypeDXMap:
      obj = ObjectMapLoadDXFile(G, (ObjectMap*)origObj, fname, state, quiet);
      break;
    case cLoadTypeFLDMap:
      obj = ObjectMapLoadFLDFile(G, (ObjectMap*)origObj, fname, state, quiet);
      break;
    case cLoadTypeBRIXMap:
      obj = ObjectMapLoadBRIXFile(G, (ObjectMap*)origObj, fname, state, quiet);
      break;
    case cLoadTypeGRDMap:
      obj = ObjectMapLoadGRDFile(G, (ObjectMap*)origObj, fname, state, quiet);
      break;
    case cLoadTypeACNTMap:
      obj = ObjectMapLoadACNTFile(G, (ObjectMap*)origObj, fname, state, quiet);
      break;
    case cLoadTypeXPLORMap:
    case cLoadTypeXPLORStr:
      obj = ObjectMapLoadXPLOR(G, (ObjectMap*)origObj, content, state, false,
                               quiet);
      break;
    case cLoadTypePHIMap:
    case cLoadTypePHIStr:
      obj = ObjectMapLoadPHI(G, (ObjectMap*)origObj, content, state, true, size,
                             quiet);
      break;
    case cLoadTypeCCP4Map:
    case cLoadTypeCCP4Str:
    case cLoadTypeCCP4Unspecified:
    case cLoadTypeCCP4UnspecifiedStr:
    case cLoadTypeMRC:
    case cLoadTypeMRCStr:
      obj = ObjectMapLoadCCP4(G, (ObjectMap*)origObj, content, state, true,
                              size, quiet, content_format);
      break;
    case cLoadTypeCGO:
      obj = ObjectCGOFromFloatArray(G, (ObjectCGO*)origObj, (float*)content,
                                    size, state, quiet);
      break;
    case cLoadTypeMOL:
    case cLoadTypeMOLStr:
    case cLoadTypeMOL2:
    case cLoadTypeMOL2Str:
    case cLoadTypeSDF2:
    case cLoadTypeSDF2Str:
    case cLoadTypeXYZ:
    case cLoadTypeXYZStr:
    case cLoadTypeMMD:
    case cLoadTypeMMDStr: {
      const char* next_entry = content;
      char new_name[WordLength] = "";

      OVLexicon* loadproplex = nullptr;
      bool loadpropertiesall = false;

      if (object_props[0]) {
        if (strcmp(object_props, "*") == 0) {
          loadpropertiesall = true;
        } else {
          auto keys = strsplit(object_props);
          loadproplex = OVLexicon_New(G->Context->heap);
          for (auto& key : keys) {
            OVLexicon_GetFromCString(loadproplex, key.c_str());
          }
        }
      }

      // (some of) these file types support multiple molecules per file,
      // and we support to load them into separate objects (multiplex).
      do {
        obj = ObjectMoleculeReadStr(G, (ObjectMolecule*)origObj, &next_entry,
                                    content_format, state, discrete, quiet,
                                    multiplex, new_name, loadpropertiesall,
                                    loadproplex);

        if (new_name[0]) {
          // multiplexing
          ObjectSetName(obj, new_name);
          ExecutiveDelete(G, obj->Name);  // just in case there is a collision
          ExecutiveManageObject(G, obj, zoom, true);
          new_name[0] = 0;
          obj = nullptr;
        }
      } while (next_entry);

      OVLexicon_Del(loadproplex);
    } break;
    default:
      if (plugin[0]) {
        obj = PlugIOManagerLoad(G, origObj ? &origObj : nullptr, fname, state,
                                quiet, plugin, args.plugin_mask);
      } else {
        return pymol::make_error("Unable to read that file type from C (",
                                 content_format, ", ", plugin, ")");
      }
  }

  if (origObj && obj) {
    if (finish) ExecutiveUpdateObjectSelection(G, origObj);

    if (fname)
      sprintf(buf, " CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
              fname, object_name, state + 1);
  } else if (obj) {
    ObjectSetName(obj, object_name);
    ExecutiveManageObject(G, obj, zoom, true);

    if (fname)
      sprintf(buf, " CmdLoad: \"%s\" loaded as \"%s\".\n", fname, obj->Name);
    else
      sprintf(buf, " CmdLoad: loaded as \"%s\".\n", obj->Name);
  }

  if (!quiet && buf[0]) {
    PRINTFB(G, FB_Executive, FB_Actions)
    "%s", buf ENDFB(G);
  }

  // Reshape movie in case multi-state.
  OrthoReshape(G, -1, -1, false);
  return {};
}

/* ExecutiveGetExistingCompatible
 *
 * PARAMS
 *  oname -- object name
 *  type  -- new object type
 *
 * RETURNS
 *  Base-class object ptr
 *
 * SIDE EFFECTS
 *  If an object with the same name but different type already exists,
 *  then it is deleted.
 */
pymol::CObject* ExecutiveGetExistingCompatible(PyMOLGlobals* G,
                                               const char* oname,
                                               cLoadType_t type) {
  pymol::CObject* origObj = nullptr;
  origObj = ExecutiveFindObjectByName(G, oname);
  /* check for existing object of right type, delete if not */
  if (origObj) {
    int new_type = -1;
    switch (type) {
      case cLoadTypePlugin:
        // let PlugIOManager delete incompatible objects
        return origObj;
      case cLoadTypeChemPyModel:
      case cLoadTypePDB:
      case cLoadTypePDBStr:
      case cLoadTypeVDBStr:
      case cLoadTypeCIF:
      case cLoadTypeCIFStr:
      case cLoadTypeXYZ:
      case cLoadTypeXYZStr:
      case cLoadTypeMOL:
      case cLoadTypeMOLStr:
      case cLoadTypeMMD:
      case cLoadTypeMMDSeparate:
      case cLoadTypeMMDStr:
      case cLoadTypeTOP:
      case cLoadTypeTRJ:
      case cLoadTypeCRD:
      case cLoadTypeRST:
      case cLoadTypeMOL2:
      case cLoadTypeMOL2Str:
      case cLoadTypeSDF2:
      case cLoadTypeSDF2Str:
      case cLoadTypePQR:
      case cLoadTypePDBQT:
      case cLoadTypeXTC:
      case cLoadTypeDTR:
      case cLoadTypeTRR:
      case cLoadTypeGRO:
      case cLoadTypeTRJ2:
      case cLoadTypeG96:
      case cLoadTypeDCD:
        new_type = cObjectMolecule;
        break;
      case cLoadTypeChemPyBrick:
      case cLoadTypeChemPyMap:
      case cLoadTypeXPLORMap:
      case cLoadTypeXPLORStr:
      case cLoadTypeCCP4Map:
      case cLoadTypeCCP4Str:
      case cLoadTypeCCP4Unspecified:
      case cLoadTypeCCP4UnspecifiedStr:
      case cLoadTypeMRC:
      case cLoadTypeMRCStr:
      case cLoadTypeFLDMap:
      case cLoadTypeBRIXMap:
      case cLoadTypeGRDMap:
      case cLoadTypeDXMap:
      case cLoadTypeDXStr:
        new_type = cObjectMap;
        break;
      case cLoadTypeCallback:
        new_type = cObjectCallback;
        break;
      case cLoadTypeCGO:
        new_type = cObjectCGO;
        break;
    }
    if (new_type == -1 || new_type != origObj->type) {
      ExecutiveDelete(G, origObj->Name);
      origObj = nullptr;
    }
  }
  return origObj;
}

int ExecutiveProcessPDBFile(PyMOLGlobals* G, pymol::CObject* origObj,
                            const char* fname, const char* buffer,
                            const char* oname, int frame, int discrete,
                            int finish, OrthoLineType buf, int variant,
                            int quiet, int multiplex, int zoom) {
  int ok = true;
  pymol::CObject* obj;
  char pdb_name[WordLength] = "";
  char cur_name[WordLength] = "";
  const char* next_pdb = nullptr;
  int repeat_flag = true;
  int n_processed = 0;
  PDBInfoRec pdb_info_rec, *pdb_info = nullptr;
  int model_number;
  pymol::CObject* deferred_zoom_obj = nullptr;

  UtilZeroMem(&pdb_info_rec, sizeof(PDBInfoRec));
  pdb_info = &pdb_info_rec;
  pdb_info->multiplex = multiplex;
  pdb_info->variant = variant;

  while (repeat_flag && ok) {
    const char* start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;
    int is_new = false;

    if (next_pdb) {
      start_at = next_pdb;
      is_repeat_pass = true;
    }

    repeat_flag = false;
    next_pdb = nullptr;
    if (!origObj) {
      is_new = true;
      pdb_name[0] = 0;
      model_number = 0;
      obj = ObjectMoleculeReadPDBStr(G, (ObjectMolecule*)origObj, start_at,
                                     eff_frame, discrete, pdb_name, &next_pdb,
                                     pdb_info, quiet, &model_number);

    } else {
      model_number = 0;
      ObjectMoleculeReadPDBStr(G, (ObjectMolecule*)origObj, start_at, eff_frame,
                               discrete, pdb_name, &next_pdb, pdb_info, quiet,
                               &model_number);

      if (finish) {
        ExecutiveUpdateObjectSelection(G, origObj);
        ExecutiveDoZoom(G, origObj, false, zoom, quiet);
      }
      if (eff_frame < 0) eff_frame = ((ObjectMolecule*)origObj)->NCSet - 1;
      if (buf) {
        if (fname)
          sprintf(buf,
                  " CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                  fname, oname, eff_frame + 1);
        else
          sprintf(
              buf,
              " CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
              oname, eff_frame + 1);
      }
      obj = origObj;
    }

    if (obj) {
      if (next_pdb) {
        /* NOTE: if set, assume that multiple PDBs are present in the file */
        repeat_flag = true;
      }
    }

    if (is_new) {
      if (obj) {
        if (next_pdb) {
          if (pdb_name[0] == 0) {
            if (cur_name[0]) {
              sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
            } else {
              sprintf(pdb_name, "%s_%04d", oname, n_processed + 1);
            }
          } else if (multiplex > 0) {
            if (pdb_info->multi_object_status ==
                1) { /* this is a multi-object PDB file */
              strcpy(cur_name, pdb_name);
            } else if (cur_name[0] == 0) {
              strcpy(cur_name, oname);
            }
            if (model_number > 0) {
              sprintf(pdb_name, "%s_%04d", cur_name, model_number);
            } else {
              sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
            }
          }
          ObjectSetName(obj, pdb_name);
          ExecutiveDelete(G, obj->Name); /* just in case */
        } else {
          if (is_repeat_pass) {
            if (pdb_name[0] == 0) {
              if (cur_name[0]) {
                sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
              } else {
                sprintf(pdb_name, "%s_%04d", oname, n_processed + 1);
              }
            } else if (multiplex > 0) {
              if (pdb_info->multi_object_status ==
                  1) { /* this is a multi-object PDB file */
                strcpy(cur_name, pdb_name);
              } else if (cur_name[0] == 0) {
                strcpy(cur_name, oname);
              }
              if (model_number > 0) {
                sprintf(pdb_name, "%s_%04d", cur_name, model_number);
              } else {
                sprintf(pdb_name, "%s_%04d", cur_name, n_processed + 1);
              }
            }
            ObjectSetName(obj, pdb_name);  /* from PDB */
            ExecutiveDelete(G, obj->Name); /* just in case */
          } else {
            ObjectSetName(obj, oname); /* from filename/parameter */
          }
        }

        if (obj) {
          int do_zoom = repeat_flag ? 0 : zoom;
          if (do_zoom != zoom)
            deferred_zoom_obj = obj;
          else
            deferred_zoom_obj = nullptr;
          ExecutiveManageObject(G, obj, do_zoom, true);
          if (eff_frame < 0) eff_frame = ((ObjectMolecule*)obj)->NCSet - 1;
          if (buf) {
            if (n_processed < 1) {
              if (fname)
                sprintf(buf, " CmdLoad: \"%s\" loaded as \"%s\".\n", fname,
                        oname);
              else
                sprintf(buf,
                        " CmdLoad: PDB-string loaded into object \"%s\", state "
                        "%d.\n",
                        oname, eff_frame + 1);
            } else {
              if (fname) {
                sprintf(buf, " CmdLoad: loaded %d objects from \"%s\".\n",
                        n_processed + 1, fname);
              } else {
                sprintf(buf, " CmdLoad: loaded %d objects from string.\n",
                        n_processed + 1);
              }
            }
          }
        }
      }
    }

    if (obj) {
      n_processed++;
    }
  }

  if (deferred_zoom_obj) {
    ExecutiveDoZoom(G, deferred_zoom_obj, true, zoom, true);
  }

  return ok;
}

static pymol::Result<pymol::CObject*> ExecutiveProcessCif(
    PyMOLGlobals* G, pymol::CObject* I, const char* object_name, const char* st,
    int frame, int discrete, int quiet, int multiplex, int zoom) {
  auto cif = std::make_shared<cif_file_with_error_capture>();
  if (!cif->parse_string(st)) {
    return pymol::make_error("Parsing CIF file failed: ", cif->m_error_msg);
  }

  auto has_refln = [](const auto& datablock) {
    return datablock.get_arr("_refln.index_h") != nullptr;
  };

  bool new_obj = I == nullptr;

  pymol::CObject* obj{};
  for (const auto& [code, datablock] : cif->datablocks()) {
    // Map Coeffs / Reflections
    if (has_refln(datablock)) {
      auto map = ObjectMapReadCifStr(G, datablock);
      p_return_if_error(map);
      obj = map.result();
    } else if (auto mol =
                   ObjectMoleculeReadCifData(G, &datablock, discrete, quiet)) {
      obj = mol;
      // we only provide access from the Python API so far
      if (SettingGet<bool>(G, cSetting_cif_keepinmemory)) {
        mol->m_cifdata = &datablock;
        mol->m_ciffile = cif;
      }
    }
    if (!obj) {
      PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
      " mmCIF-Warning: no map coefficients or coordinates found in data_%s\n",
          datablock.code() ENDFB(G);
      continue;
    }

    if (cif->datablocks().size() == 1 || multiplex == 0) {
      return obj;
    }

    // multiplexing
    ObjectSetName(obj, datablock.code());
    ExecutiveDelete(G, obj->Name);
    ExecutiveManageObject(G, obj, zoom, true);
  }  // end for

  return obj;
}

pymol::Result<> ExecutiveLoadCoordset(PyMOLGlobals* G,
                                      pymol::zstring_view oname,
                                      PyObject* model, int frame, bool quiet) {
  auto origObj = ExecutiveFindObjectByName(G, oname.c_str());
  if (!origObj || origObj->type != cObjectMolecule) {
    return pymol::make_error("Invalid object molecule.");
  }

  PBlock(G);
  auto obj =
      ObjectMoleculeLoadCoords(G, (ObjectMolecule*)origObj, model, frame);
  PUnblock(G);
  if (!obj) {
    return pymol::make_error("Load Coordset Error");
  }
  if (frame < 0) frame = obj->NCSet - 1;

  if (!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
    " CmdLoad: Coordinates appended into object \"%s\", state %d.\n",
        oname.c_str(), frame + 1 ENDFB(G);
  }
  return {};
}

pymol::Result<> ExecutiveLoadTraj(PyMOLGlobals* G, pymol::zstring_view oname,
                                  pymol::zstring_view fname, int frame,
                                  int type, int interval, int average,
                                  int start, int stop, int max,
                                  pymol::zstring_view str1, int image,
                                  const float* shift,
                                  pymol::zstring_view plugin, int quiet) {
  auto s1 = SelectorTmp::make(G, str1.c_str());
  p_return_if_error(s1);
  bool ok = true;
  auto origObj = ExecutiveFindObjectByName(G, oname.c_str());
  if (!origObj) {
    return pymol::make_error(
        "Must load object topology before loading trajectory.");
  }
  if (origObj->type != cObjectMolecule) {
    return pymol::make_error("Object '", oname, "' is not a molecular object.");
  }
  if ((type == cLoadTypeTRJ) && (!plugin.empty())) type = cLoadTypeTRJ2;
  switch (type) {
    case cLoadTypeTRJ: /* this is the ascii AMBER trajectory format... */
      PRINTFD(G, FB_CCmd) " ExecutiveLoadTraj-DEBUG: loading TRJ\n" ENDFD;
      ObjectMoleculeLoadTRJFile(G, (ObjectMolecule*)origObj, fname.c_str(),
                                frame, interval, average, start, stop, max,
                                s1->getName(), image, shift, quiet);
      PRINTFB(G, FB_Executive, FB_Actions)
      " ExecutiveLoadTraj: \"%s\" appended into object \"%s\".\n"
      " ExecutiveLoadTraj: %d total states in the object.\n",
          fname.c_str(), oname.c_str(),
          ((ObjectMolecule*)origObj)->NCSet ENDFB(G);
      break;
    default:
      ok = PlugIOManagerLoadTraj(
          G, (ObjectMolecule*)origObj, fname.c_str(), frame, interval, average,
          start, stop, max, s1->getName(), image, shift, quiet, plugin.c_str());
  }
  if (ok) {
    return {};
  } else {
    return pymol::make_error("Could not load trajectory");
  }
}
