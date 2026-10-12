#pragma once

#include <unordered_map>

#include "pymol/ui/ortho.h"
#include "pymol/ui/scroll_bar.h"

class SpecRec;
class CGO;
struct CTracker;
struct OVLexicon;

struct ExecutiveObjectOffset {
  ObjectMolecule* obj;
  int atm;
};

struct PanelRec {
  SpecRec* spec;
  unsigned nest_level;
  bool is_group = false;
  bool is_open = false;

  PanelRec(SpecRec* spec_, unsigned nest_level_)
      : spec(spec_), nest_level(nest_level_) {}
};

struct ListMember {
  int list_id;
  int next;
};

struct CExecutive : public Block {
  SpecRec* Spec{};
  CTracker* Tracker{};
  int Width{}, Height{};
  ScrollBar m_ScrollBar;
  pymol::CObject* LastEdited{nullptr};
  int all_names_list_id{}, all_obj_list_id{}, all_sel_list_id{};
  OVLexicon* Lex{};
  std::unordered_map<ov_word, int> Key;
  bool ValidGroups{false};
  bool ValidSceneMembers{false};
  int ValidGridSlots{};

  std::vector<PanelRec> Panel{};

  int CaptureFlag{};
  int LastMotionCount{};
  CGO* selIndicatorsCGO{nullptr};

  // AtomInfoType::unique_id -> (object, atom-index)
  std::vector<ExecutiveObjectOffset> m_eoo{};  // vector of (object, atom-index)
  std::unordered_map<ov_word, std::size_t>
      m_id2eoo{};  // unique_id -> m_eoo-index
  std::unordered_map<const pymol::CObject*,
                     std::unordered_set<const pymol::CObject*>>
      m_objDeps;

  CExecutive(PyMOLGlobals* G) : Block(G), m_ScrollBar(G, false) {};

  void reshape(int width, int height) override;
};
