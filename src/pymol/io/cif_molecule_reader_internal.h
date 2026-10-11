/*
 * Read a molecule from CIF
 *
 * (c) 2014 Schrodinger, Inc.
 */

// Helpers shared by the mmCIF reader files; not part of the reader API.

#pragma once

#include "pymol/objects/assembly_helpers.h"
#include "pymol/io/cif_file.h"
#include "pymol/io/cif_bond_dict.h"

class seqvec_t;
struct CifContentInfo;

using pymol::cif_array;

// dictionary content types
enum CifDataType {
  CIF_UNKNOWN,
  CIF_CORE,      // small molecule
  CIF_MMCIF,     // macromolecular structure
  CIF_CHEM_COMP  // chemical component
};

// Defined in cif_molecule_reader.cc.
std::string make_mm_atom_site_label(PyMOLGlobals* G, const AtomInfoType* a);
std::string make_mm_atom_site_label(PyMOLGlobals* G, const char* asym_id,
                                    const char* comp_id, const char* seq_id,
                                    const char* ins_code, const char* atom_id,
                                    const char* alt_id);
CSymmetry* read_symmetry(PyMOLGlobals* G, const pymol::cif_data* data);

// Defined in cif_bonds.cc.
void ObjectMoleculeAddBond2(ObjectMolecule* I, int i1, int i2, int order);
bool read_chem_comp_bond_dict(const pymol::cif_data* data,
                              bond_dict_t& bond_dict);
int ObjectMoleculeConnectComponents(ObjectMolecule* I,
                                    bond_dict_t* bond_dict = nullptr);
pymol::vla<BondType> read_geom_bond(PyMOLGlobals* G,
                                    const pymol::cif_data* data,
                                    const pymol::vla<AtomInfoType>& atInfo);
pymol::vla<BondType> read_chemical_conn_bond(PyMOLGlobals* G,
                                             const pymol::cif_data* data);
bool read_struct_conn_(PyMOLGlobals* G, const pymol::cif_data* data,
                       const pymol::vla<AtomInfoType>& atInfo, CoordSet* cset,
                       CifContentInfo& info);
pymol::vla<BondType> read_chem_comp_bond(
    PyMOLGlobals* G, const pymol::cif_data* data,
    const pymol::vla<AtomInfoType>& atInfo);
pymol::vla<BondType> read_pymol_bond(PyMOLGlobals* G,
                                     const pymol::cif_data* data,
                                     const pymol::vla<AtomInfoType>& atInfo);

// simple 1-indexed string storage
class seqvec_t : public std::vector<std::string> {
 public:
  void set(int i, const char* mon_id) {
    if (i < 1) {
      printf("error: i(%d) < 1\n", i);
      return;
    }
    if (i > size()) resize(i);
    (*this)[i - 1] = mon_id;
  }

  const char* get(int i) const {
    if (i < 1 || i > size()) return nullptr;
    return (*this)[i - 1].c_str();
  }
};

// structure to collect information about a data block
struct CifContentInfo {
  PyMOLGlobals* G;
  CifDataType type;
  bool fractional;
  bool use_auth;
  std::set<lexidx_t> chains_filter;
  std::set<std::string> polypeptide_entities;  // entity ids
  std::map<std::string, seqvec_t>
      sequences;  // entity_id -> [resn1, resn2, ...]

  // _pdbx_poly_seq_scheme: label_asym_id -> (label_seq_id -> auth numbering)
  struct SeqSchemeEntry {
    int auth_seq_id;
    char ins_code;
  };
  std::unordered_map<std::string, std::unordered_map<int, SeqSchemeEntry>>
      seq_scheme;

  bool is_excluded_chain(const char* chain) const {
    if (chains_filter.empty()) return false;

    auto borrowed = LexBorrow(G, chain);
    if (borrowed != LEX_BORROW_NOTFOUND) return is_excluded_chain(borrowed);

    return false;
  }

  bool is_excluded_chain(const lexborrow_t& chain) const {
    return (!chains_filter.empty() &&
            chains_filter.count(reinterpret_cast<const lexidx_t&>(chain)) == 0);
  }

  bool is_polypeptide(const char* entity_id) const {
    return polypeptide_entities.count(entity_id);
  }

  CifContentInfo(PyMOLGlobals* G, bool use_auth = true)
      : G(G), type(CIF_UNKNOWN), fractional(false), use_auth(use_auth) {}
};

/**
 * Lookup one key in a map, return true if found and
 * assign output reference `value1`
 */
template <typename Map, typename Key, typename T>
inline bool find1(Map& dict, T& value1, const Key& key1) {
  auto it = dict.find(key1);
  if (it == dict.end()) return false;
  value1 = it->second;
  return true;
}
