/*
 * Read a molecule from CIF
 *
 * (c) 2014 Schrodinger, Inc.
 */

// Bonds for mmCIF structures: the chemical component dictionary, struct_conn
// and geometry.

#include "pymol/python/p.h"
#include "pymol/base/util2.h"
#include "pymol/base/feedback.h"
#include "pymol/io/cif_molecule_reader_internal.h"

/**
 * Lookup two keys in a map, return true if both found and
 * assign output references `value1` and `value2`.
 */
template <typename Map, typename Key, typename T>
inline bool find2(Map& dict, T& value1, const Key& key1, T& value2,
                  const Key& key2) {
  if (!find1(dict, value1, key1)) return false;
  if (!find1(dict, value2, key2)) return false;
  return true;
}

/**
 * Initialize a bond. Only one of symmetry_1 or symmetry_2 must be non-default.
 * If symmetry_2 is default and symmetry_1 is non-default, then swap the
 * indices.
 */
static bool BondTypeInit3(PyMOLGlobals* G, BondType* bond, unsigned i1,
                          unsigned i2, const char* symmetry_1,
                          const char* symmetry_2, int order = 1) {
  auto symop_1 = pymol::SymOp(symmetry_1);
  auto symop_2 = pymol::SymOp(symmetry_2);

  if (symop_1) {
    if (symop_2) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: Bonds with two symmetry operations not supported\n" ENDFB(G);
      return false;
    }

    std::swap(i1, i2);
    std::swap(symop_1, symop_2);
  }

  BondTypeInit2(bond, i1, i2, order);
  bond->symop_2 = symop_2;

  return true;
}

/**
 * Add one bond without checking if it already exists
 */
void ObjectMoleculeAddBond2(ObjectMolecule* I, int i1, int i2, int order) {
  VLACheck(I->Bond, BondType, I->NBond);
  BondTypeInit2(I->Bond + I->NBond, i1, i2, order);
  I->NBond++;
}

/**
 * Get the distance between two atoms in ObjectMolecule
 */
static float GetDistance(ObjectMolecule* I, int i1, int i2) {
  const CoordSet* cset;
  int idx1 = -1, idx2 = -1;

  // find first coordset which contains both atoms
  if (I->DiscreteFlag) {
    cset = I->DiscreteCSet[i1];
    if (cset == I->DiscreteCSet[i2]) {
      idx1 = I->DiscreteAtmToIdx[i1];
      idx2 = I->DiscreteAtmToIdx[i2];
    }
  } else {
    for (int i = 0; i < I->NCSet; ++i) {
      if ((cset = I->CSet[i])) {
        if ((idx1 = cset->AtmToIdx[i1]) != -1 &&
            (idx2 = cset->AtmToIdx[i2]) != -1) {
          break;
        }
      }
    }
  }

  if (idx1 == -1 || idx2 == -1) return 999.f;

  float v[3];
  subtract3f(cset->coordPtr(idx1), cset->coordPtr(idx2), v);
  return length3f(v);
}

/**
 * Bond order string to int
 */
static int bondOrderLookup(const char* order) {
  if (p_strcasestartswith(order, "doub")) return 2;
  if (p_strcasestartswith(order, "trip")) return 3;
  if (p_strcasestartswith(order, "arom")) return 4;
  if (p_strcasestartswith(order, "delo")) return 4;
  // single
  return 1;
}

/**
 * Read bonds from CHEM_COMP_BOND in `bond_dict` dictionary
 */
bool read_chem_comp_bond_dict(const pymol::cif_data* data,
                              bond_dict_t& bond_dict) {
  const cif_array *arr_id_1, *arr_id_2, *arr_order, *arr_comp_id;

  if (!(arr_id_1 = data->get_arr("_chem_comp_bond.atom_id_1")) ||
      !(arr_id_2 = data->get_arr("_chem_comp_bond.atom_id_2")) ||
      !(arr_order = data->get_arr("_chem_comp_bond.value_order")) ||
      !(arr_comp_id = data->get_arr("_chem_comp_bond.comp_id"))) {
    if ((arr_comp_id = data->get_arr("_chem_comp_atom.comp_id"))) {
      // atom(s) but no bonds (e.g. metals)
      bond_dict.set_unknown(arr_comp_id->as_s());
      return true;
    }

    return false;
  }

  const char *name1, *name2, *resn;
  int order_value;
  int nrows = arr_id_1->size();

  for (int i = 0; i < nrows; i++) {
    resn = arr_comp_id->as_s(i);
    name1 = arr_id_1->as_s(i);
    name2 = arr_id_2->as_s(i);

    const char* order = arr_order->as_s(i);
    order_value = bondOrderLookup(order);

    bond_dict[resn].set(name1, name2, order_value);
  }

  // alt_atom_id -> atom_id
  if ((arr_comp_id = data->get_arr("_chem_comp_atom.comp_id")) &&
      (arr_id_1 = data->get_arr("_chem_comp_atom.atom_id")) &&
      (arr_id_2 = data->get_arr("_chem_comp_atom.alt_atom_id"))) {
    nrows = arr_id_1->size();

    // set of all non-alt ids
    std::set<pymol::zstring_view> atom_ids;
    for (int i = 0; i < nrows; ++i) {
      atom_ids.insert(arr_id_1->as_s(i));
    }

    for (int i = 0; i < nrows; ++i) {
      resn = arr_comp_id->as_s(i);
      name1 = arr_id_1->as_s(i);
      name2 = arr_id_2->as_s(i);

      // skip identity mapping
      if (strcmp(name1, name2) == 0) {
        continue;
      }

      // alt id must not also be a non-alt id (PYMOL-3470)
      if (atom_ids.count(name2)) {
        fprintf(
            stderr,
            "Warning: _chem_comp_atom.alt_atom_id %s/%s ignored for bonding\n",
            resn, name2);
        continue;
      }

      bond_dict[resn].add_alt_name(name1, name2);
    }
  }

  return true;
}

/**
 * parse $PYMOL_DATA/chem_comp_bond-top100.cif (subset of components.cif) into
 * a static (global) dictionary.
 */
static bond_dict_t* get_global_components_bond_dict(PyMOLGlobals* G) {
  static bond_dict_t bond_dict;

  if (bond_dict.empty()) {
    const char* pymol_data = getenv("PYMOL_DATA");
    if (!pymol_data || !pymol_data[0]) return nullptr;

    std::string path(pymol_data);
    path.append(PATH_SEP).append("chem_comp_bond-top100.cif");
    cif_file_with_error_capture cif;
    if (!cif.parse_file(path.c_str())) {
      PRINTFB(G, FB_Executive, FB_Warnings)
      " Warning: Loading '%s' failed: %s\n", path.c_str(),
          cif.m_error_msg.c_str() ENDFB(G);
      return nullptr;
    }

    for (const auto& [code, datablock] : cif.datablocks()) {
      read_chem_comp_bond_dict(&datablock, bond_dict);
    }
  }

  return &bond_dict;
}

/**
 * True for N-H1 and N-H3, those are not in the chemical components dictionary.
 */
static bool is_N_H1_or_H3(PyMOLGlobals* G, const AtomInfoType* a1,
                          const AtomInfoType* a2) {
  if (a2->name == G->lex_const.N) {
    a2 = a1;
  } else if (a1->name != G->lex_const.N) {
    return false;
  }

  return (a2->name == G->lex_const.H1 || a2->name == G->lex_const.H3);
}

/**
 * Add bonds for one residue, with atoms spanning from i_start to i_end-1,
 * based on components.cif
 */
static void ConnectComponent(ObjectMolecule* I, int i_start, int i_end,
                             bond_dict_t* bond_dict) {
  if (i_end - i_start < 2) return;

  auto G = I->G;
  const AtomInfoType *a1, *a2, *ai = I->AtomInfo.data();
  int order;

  // get residue bond dictionary
  auto res_dict = bond_dict->get(G, LexStr(G, ai[i_start].resn));
  if (res_dict == nullptr) return;

  // for all pairs of atoms in given set
  for (int i1 = i_start + 1; i1 < i_end; i1++) {
    for (int i2 = i_start; i2 < i1; i2++) {
      a1 = ai + i1;
      a2 = ai + i2;

      // don't connect different alt codes
      if (a1->alt[0] && a2->alt[0] && strcmp(a1->alt, a2->alt) != 0) {
        continue;
      }

      // restart if we hit the next residue in bulk solvent (atoms must
      // not be sorted for this)
      // TODO artoms are sorted at this point
      if (a1->name == a2->name) {
        i_start = i1;
        break;
      }

      // lookup if atoms are bonded
      order = res_dict->get(LexStr(G, a1->name), LexStr(G, a2->name));
      if (order < 0) {
        if (!is_N_H1_or_H3(G, a1, a2) || GetDistance(I, i1, i2) > 1.2) continue;

        order = 1;
      }

      // make bond
      ObjectMoleculeAddBond2(I, i1, i2, order);
    }
  }
}

/**
 * Add intra residue bonds based on components.cif, and common polymer
 * connecting bonds (C->N, O3*->P)
 */
int ObjectMoleculeConnectComponents(ObjectMolecule* I, bond_dict_t* bond_dict) {
  PyMOLGlobals* G = I->G;
  int i_start = 0;
  std::vector<int> i_prev_c[2], i_prev_o3[2];

  const lexborrow_t lex_O3s = LexBorrow(G, "O3*");
  const lexborrow_t lex_O3p = LexBorrow(G, "O3'");

  if (!bond_dict) {
    // read components.cif
    if (!(bond_dict = get_global_components_bond_dict(G))) return false;
  }

  // reserve some memory for new bonds
  I->Bond.reserve(I->NAtom * 4);

  for (int i = 0; i < I->NAtom; ++i) {
    auto const& atom = I->AtomInfo[i];

    // intra-residue
    if (!AtomInfoSameResidue(G, I->AtomInfo + i_start, I->AtomInfo + i)) {
      ConnectComponent(I, i_start, i, bond_dict);
      i_start = i;
      i_prev_c[0] = std::move(i_prev_c[1]);
      i_prev_o3[0] = std::move(i_prev_o3[1]);
      i_prev_c[1].clear();
      i_prev_o3[1].clear();
    }

    // inter-residue polymer bonds
    if (atom.name == G->lex_const.C) {
      i_prev_c[1].push_back(i);
    } else if (atom.name == lex_O3s || atom.name == lex_O3p) {
      i_prev_o3[1].push_back(i);
    } else {
      auto const* i_prev_ptr = (atom.name == G->lex_const.N)   ? i_prev_c
                               : (atom.name == G->lex_const.P) ? i_prev_o3
                                                               : nullptr;

      if (i_prev_ptr && !i_prev_ptr->empty()) {
        for (int i_prev : *i_prev_ptr) {
          bool alt_check = !atom.alt[0] || !I->AtomInfo[i_prev].alt[0] ||
                           atom.alt[0] == I->AtomInfo[i_prev].alt[0];
          if (alt_check && GetDistance(I, i_prev, i) < 1.8) {
            // make bond
            ObjectMoleculeAddBond2(I, i_prev, i, 1);
          }
        }
      }
    }
  }

  // final residue
  ConnectComponent(I, i_start, I->NAtom, bond_dict);

  // clean up
  VLASize(I->Bond, BondType, I->NBond);

  return true;
}

/**
 * Read GEOM_BOND
 *
 * return: BondType VLA
 */
pymol::vla<BondType> read_geom_bond(PyMOLGlobals* G,
                                    const pymol::cif_data* data,
                                    const pymol::vla<AtomInfoType>& atInfo) {
  const cif_array *arr_ID_1, *arr_ID_2;
  if ((arr_ID_1 = data->get_arr("_geom_bond.atom_site_id_1",
                                "_geom_bond_atom_site_label_1")) == nullptr ||
      (arr_ID_2 = data->get_arr("_geom_bond.atom_site_id_2",
                                "_geom_bond_atom_site_label_2")) == nullptr)
    return {};

  const cif_array* arr_symm_1 = data->get_opt("_geom_bond?site_symmetry_1");
  const cif_array* arr_symm_2 = data->get_opt("_geom_bond?site_symmetry_2");

  int nrows = arr_ID_1->size();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  auto bondvla = pymol::vla<BondType>(6 * nAtom);

  // name -> atom index
  std::map<std::string, int> name_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    std::string key(LexStr(G, atInfo[i].name));
    name_dict[key] = i;
  }

  // read table
  for (int i = 0; i < nrows; i++) {
    std::string key1(arr_ID_1->as_s(i));
    std::string key2(arr_ID_2->as_s(i));

    int i1, i2;
    if (find2(name_dict, i1, key1, i2, key2)) {
      auto const bond = bondvla.check(nBond);
      if (BondTypeInit3(G, bond, i1, i2,      //
                        arr_symm_1->as_s(i),  //
                        arr_symm_2->as_s(i))) {
        ++nBond;
      }
    } else {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Detail: _geom_bond name lookup failed: %s %s\n", key1.c_str(),
          key2.c_str() ENDFB(G);
    }
  }

  if (nBond) {
    VLASize(bondvla, BondType, nBond);
  } else {
    VLAFreeP(bondvla);
  }

  return bondvla;
}

/**
 * Read CHEMICAL_CONN_BOND
 *
 * return: BondType VLA
 */
pymol::vla<BondType> read_chemical_conn_bond(PyMOLGlobals* G,
                                             const pymol::cif_data* data) {
  const cif_array *arr_number, *arr_atom_1, *arr_atom_2, *arr_type;

  if ((arr_number = data->get_arr("_atom_site?chemical_conn_number")) ==
          nullptr ||
      (arr_atom_1 = data->get_arr("_chemical_conn_bond?atom_1")) == nullptr ||
      (arr_atom_2 = data->get_arr("_chemical_conn_bond?atom_2")) == nullptr ||
      (arr_type = data->get_arr("_chemical_conn_bond?type")) == nullptr)
    return {};

  int nAtom = arr_number->size();
  int nBond = arr_atom_1->size();

  auto bondvla = pymol::vla<BondType>(nBond);
  auto bond = bondvla.data();

  // chemical_conn_number -> atom index
  std::map<int, int> number_dict;

  // build dictionary
  for (int i = 0; i < nAtom; i++) {
    number_dict[arr_number->as_i(i)] = i;
  }

  // read table
  int i1, i2;
  for (int i = 0; i < nBond; i++) {
    if (find2(number_dict, i1, arr_atom_1->as_i(i), i2, arr_atom_2->as_i(i))) {
      BondTypeInit2(bond++, i1, i2, bondOrderLookup(arr_type->as_s(i)));
    } else {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Detail: _chemical_conn_bond name lookup failed\n" ENDFB(G);
    }
  }

  return bondvla;
}

/**
 * Read bonds from STRUCT_CONN
 *
 * Output:
 *   cset->TmpBond
 *   cset->NTmpBond
 */
bool read_struct_conn_(PyMOLGlobals* G, const pymol::cif_data* data,
                       const pymol::vla<AtomInfoType>& atInfo, CoordSet* cset,
                       CifContentInfo& info) {
  const cif_array* col_type_id = data->get_arr("_struct_conn.conn_type_id");

  if (!col_type_id) return false;

  const cif_array *col_asym_id[2] = {nullptr, nullptr},
                  *col_comp_id[2] = {nullptr, nullptr},
                  *col_seq_id[2] = {nullptr, nullptr},
                  *col_atom_id[2] = {nullptr, nullptr},
                  *col_alt_id[2] = {nullptr, nullptr},
                  *col_ins_code[2] = {nullptr, nullptr},
                  *col_symm[2] = {nullptr, nullptr};

  if (info.use_auth) {
    col_asym_id[0] = data->get_arr("_struct_conn.ptnr1_auth_asym_id");
    col_comp_id[0] = data->get_arr("_struct_conn.ptnr1_auth_comp_id");
    col_seq_id[0] = data->get_arr("_struct_conn.ptnr1_auth_seq_id");
    col_atom_id[0] = data->get_arr("_struct_conn.ptnr1_auth_atom_id");
    col_asym_id[1] = data->get_arr("_struct_conn.ptnr2_auth_asym_id");
    col_comp_id[1] = data->get_arr("_struct_conn.ptnr2_auth_comp_id");
    col_seq_id[1] = data->get_arr("_struct_conn.ptnr2_auth_seq_id");
    col_atom_id[1] = data->get_arr("_struct_conn.ptnr2_auth_atom_id");

    col_alt_id[0] = data->get_arr("_struct_conn.pdbx_ptnr1_auth_alt_id");
    col_alt_id[1] = data->get_arr("_struct_conn.pdbx_ptnr2_auth_alt_id");

    // auth only
    col_ins_code[0] = data->get_arr("_struct_conn.pdbx_ptnr1_pdb_ins_code");
    col_ins_code[1] = data->get_arr("_struct_conn.pdbx_ptnr2_pdb_ins_code");
  }

  // for assembly chain filtering
  const cif_array* col_label_asym_id[2] = {
      data->get_arr("_struct_conn.ptnr1_label_asym_id"),
      data->get_arr("_struct_conn.ptnr2_label_asym_id")};

  if ((!col_asym_id[0] && !(col_asym_id[0] = col_label_asym_id[0])) ||
      (!col_comp_id[0] &&
       !(col_comp_id[0] = data->get_arr("_struct_conn.ptnr1_label_comp_id"))) ||
      (!col_seq_id[0] &&
       !(col_seq_id[0] = data->get_arr("_struct_conn.ptnr1_label_seq_id"))) ||
      (!col_atom_id[0] &&
       !(col_atom_id[0] = data->get_arr("_struct_conn.ptnr1_label_atom_id"))) ||
      (!col_asym_id[1] && !(col_asym_id[1] = col_label_asym_id[1])) ||
      (!col_comp_id[1] &&
       !(col_comp_id[1] = data->get_arr("_struct_conn.ptnr2_label_comp_id"))) ||
      (!col_seq_id[1] &&
       !(col_seq_id[1] = data->get_arr("_struct_conn.ptnr2_label_seq_id"))) ||
      (!col_atom_id[1] &&
       !(col_atom_id[1] = data->get_arr("_struct_conn.ptnr2_label_atom_id"))))
    return false;

  if (!col_alt_id[0])
    col_alt_id[0] = data->get_opt("_struct_conn.pdbx_ptnr1_label_alt_id");
  if (!col_alt_id[1])
    col_alt_id[1] = data->get_opt("_struct_conn.pdbx_ptnr2_label_alt_id");

  col_symm[0] = data->get_opt("_struct_conn.ptnr1_symmetry");
  col_symm[1] = data->get_opt("_struct_conn.ptnr2_symmetry");

  const cif_array* col_order = data->get_opt("_struct_conn.pdbx_value_order");

  int nrows = col_type_id->size();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  cset->TmpBond = pymol::vla<BondType>(6 * nAtom);

  // identifiers -> coord set index
  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    int idx = cset->atmToIdx(i);
    if (idx != -1) name_dict[make_mm_atom_site_label(G, atInfo + i)] = idx;
  }

  for (int i = 0; i < nrows; i++) {
    const char* type_id = col_type_id->as_s(i);
    if (strncasecmp(type_id, "covale", 6) && strcasecmp(type_id, "modres") &&
        strcasecmp(type_id, "disulf"))
      // ignore non-covalent bonds (saltbr, hydrog)
      continue;

    std::string key[2];
    for (int j = 0; j < 2; j++) {
      const char* asym_id = col_asym_id[j]->as_s(i);

      if (col_label_asym_id[j] &&
          info.is_excluded_chain(col_label_asym_id[j]->as_s(i)))
        goto next_row;

      // doen't work with label_seq_id and bulk solvent
      const char* seq_id = col_seq_id[j]->as_s(i);
      if (!seq_id[0]) goto next_row;

      key[j] = make_mm_atom_site_label(
          G, asym_id, col_comp_id[j]->as_s(i), seq_id,
          col_ins_code[j] ? col_ins_code[j]->as_s(i) : "",
          col_atom_id[j]->as_s(i), col_alt_id[j]->as_s(i));
    }

    int i1, i2;
    if (find2(name_dict, i1, key[0], i2, key[1])) {
      // zero-order bond for metal coordination
      int order = strcasecmp(type_id, "metalc") ? 1 : 0;

      if (order) {
        order = bondOrderLookup(col_order->as_s(i));
      }

      auto const bond = cset->TmpBond.check(nBond);
      if (BondTypeInit3(G, bond, i1, i2,       //
                        col_symm[0]->as_s(i),  //
                        col_symm[1]->as_s(i), order)) {
        ++nBond;
      }
    } else {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Detail: _struct_conn name lookup failed: %s %s\n",
          key[0].c_str(), key[1].c_str() ENDFB(G);
    }

    // label to "continue" from inner for-loop
  next_row:;
  }

  if (nBond) {
    VLASize(cset->TmpBond, BondType, nBond);
    cset->NTmpBond = nBond;
  } else {
    VLAFreeP(cset->TmpBond);
  }

  return true;
}

/**
 * Read bonds from CHEM_COMP_BOND
 *
 * return: BondType VLA
 */
pymol::vla<BondType> read_chem_comp_bond(
    PyMOLGlobals* G, const pymol::cif_data* data,
    const pymol::vla<AtomInfoType>& atInfo) {
  const cif_array *col_ID_1, *col_ID_2, *col_comp_id;

  if ((col_ID_1 = data->get_arr("_chem_comp_bond.atom_id_1")) == nullptr ||
      (col_ID_2 = data->get_arr("_chem_comp_bond.atom_id_2")) == nullptr ||
      (col_comp_id = data->get_arr("_chem_comp_bond.comp_id")) == nullptr)
    return {};

  // "_chem_comp_bond.type" seems to be non-standard here. It's found in the
  // wild with values like "double" and "aromatic". mmcif_nmr-star.dic defines
  // it, but with different vocabulary (e.g. "amide", "ether", etc.).

  const cif_array* col_order =
      data->get_opt("_chem_comp_bond.value_order", "_chem_comp_bond.type");

  int nrows = col_ID_1->size();
  int nAtom = VLAGetSize(atInfo);
  int nBond = 0;

  auto bondvla = pymol::vla<BondType>(6 * nAtom);
  auto bond = bondvla.data();

  // name -> atom index
  std::map<std::string, int> name_dict;

  for (int i = 0; i < nAtom; i++) {
    std::string key(LexStr(G, atInfo[i].name));
    name_dict[key] = i;
  }

  for (int i = 0; i < nrows; i++) {
    std::string key1(col_ID_1->as_s(i));
    std::string key2(col_ID_2->as_s(i));
    const char* order = col_order->as_s(i);

    int i1, i2;
    if (find2(name_dict, i1, key1, i2, key2)) {
      int order_value = bondOrderLookup(order);

      nBond++;
      BondTypeInit2(bond++, i1, i2, order_value);

    } else {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Detail: _chem_comp_bond name lookup failed: %s %s\n",
          key1.c_str(), key2.c_str() ENDFB(G);
    }
  }

  if (nBond) {
    VLASize(bondvla, BondType, nBond);
  } else {
    VLAFreeP(bondvla);
  }

  return bondvla;
}

/**
 * Read bonds from _pymol_bond (non-standard extension)
 *
 * return: BondType VLA
 */
pymol::vla<BondType> read_pymol_bond(PyMOLGlobals* G,
                                     const pymol::cif_data* data,
                                     const pymol::vla<AtomInfoType>& atInfo) {
  const cif_array *col_ID_1, *col_ID_2, *col_order;

  if ((col_ID_1 = data->get_arr("_pymol_bond.atom_site_id_1")) == nullptr ||
      (col_ID_2 = data->get_arr("_pymol_bond.atom_site_id_2")) == nullptr ||
      (col_order = data->get_arr("_pymol_bond.order")) == nullptr)
    return {};

  int nrows = col_ID_1->size();
  int nAtom = VLAGetSize(atInfo);

  auto bondvla = pymol::vla<BondType>(nrows);
  auto bond = bondvla.data();

  // ID -> atom index
  std::map<int, int> id_dict;

  for (int atm = 0; atm < nAtom; ++atm) {
    id_dict[atInfo[atm].id] = atm;
  }

  for (int i = 0; i < nrows; i++) {
    auto key1 = col_ID_1->as_i(i);
    auto key2 = col_ID_2->as_i(i);
    auto order_value = col_order->as_i(i);

    int i1, i2;
    if (find2(id_dict, i1, key1, i2, key2)) {
      BondTypeInit2(bond++, i1, i2, order_value);
    } else {
      PRINTFB(G, FB_Executive, FB_Details)
      " Executive-Detail: _pymol_bond name lookup failed: %d %d\n", key1,
          key2 ENDFB(G);
    }
  }

  return bondvla;
}

/**
 * Bond dictionary getter, with on-demand download of residue dictionaries
 */
const bond_dict_t::mapped_type* bond_dict_t::get(PyMOLGlobals* G,
                                                 const char* resn,
                                                 bool try_download) {
  auto key = make_key(resn);
  auto it = m_map.find(key);

  if (it != m_map.end()) return &it->second;

  if (unknown_resn.count(key)) return nullptr;

  if (try_download) {
    pymol::GIL_Ensure gil;
    bool downloaded = false;

    // call into Python
    unique_PyObject_ptr pyfilename(PyObject_CallMethod(
        G->P_inst->cmd, "download_chem_comp", "siO", resn,
        !Feedback(G, FB_Executive, FB_Details), G->P_inst->cmd));

    if (pyfilename) {
      const char* filename = PyString_AsString(pyfilename.get());

      // update
      if ((downloaded = (filename && filename[0]))) {
        cif_file_with_error_capture cif;
        if (!cif.parse_file(filename)) {
          PRINTFB(G, FB_Executive, FB_Warnings)
          " Warning: Loading _chem_comp_bond CIF data for residue '%s' "
          "failed: %s\n",
              resn, cif.m_error_msg.c_str() ENDFB(G);
          return nullptr;
        }

        for (auto& [code, item] : cif.datablocks())
          read_chem_comp_bond_dict(&item, *this);
      }
    }

    if (downloaded) {
      // second attempt to look up, from eventually updated dictionary
      return get(G, resn, false);
    }
  }

  PRINTFB(G, FB_Executive, FB_Warnings)
  " ExecutiveLoad-Warning: No _chem_comp_bond data for residue '%s'\n",
      resn ENDFB(G);

  // don't try downloading again
  unknown_resn.insert(key);

  return nullptr;
}
