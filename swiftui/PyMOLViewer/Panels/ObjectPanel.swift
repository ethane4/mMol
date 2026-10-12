// ObjectPanel.swift — Object/selection list with A/S/H/L/C action buttons
// SwiftUI replacement for modules/pymol/appkit_object_panel.py

import AppKit
import Combine
import SwiftUI

// Per-level indent for nested group members (#255), and the depth at which
// indenting stops. Same width as the disclosure column every row already leads
// with, so a member's chevron lines up under its parent's name. Capped because
// the panel is ~300pt — past three levels the name column would be crushed,
// so deeper nesting stops moving right.
private let kIndentW: CGFloat = 13
private let kMaxIndentDepth: Int = 3

/// Atom predicate behind the S ▸ / H ▸ "side chains" menu items.
///
/// PyMOL's `sidechain` selection drops two backbone atoms the sticks need:
///
/// - **CA**, always — without it the CA–CB bond is missing and every side chain
///   floats detached from the backbone.
/// - **N**, in proline — PRO's ring closes CD back onto the backbone N, so
///   without N the ring renders as an open, dangling zig-zag (#405). The
///   `neighbor` clause restores it there and nowhere else: in every other
///   residue N's heavy neighbours are CA and the preceding C, neither of which
///   is a sidechain atom. Written as a bond query rather than `resn PRO` so HYP
///   and other N-alkylated residues are covered too.
///
/// `and not hydro` inside the neighbour test is load-bearing: on a hydrogenated
/// N-terminus `h_add` produces an extra amide H that PyMOL classifies as
/// *sidechain*, so a bare `neighbor sidechain` would drag residue 1's backbone N
/// into every hydrogenated structure.
///
/// Show and hide must use the identical predicate, or hiding leaves proline's
/// N stick behind.
private let kSidechainPred =
  "(sidechain or name CA or (name N and neighbor (sidechain and not hydro)))"

// MARK: - Representation inspector: polled state models

/// One active representation on an object, with current setting values + color.
struct RepState: Equatable {
  let rep: String  // "cartoon", "surface", …
  var visible: Bool
  var values: [String: Double]  // setting name → current value
  var color: String  // "inherit" or "#rrggbb"
  var effColor: String? = nil  // what it shows: "#rrggbb" (the atoms' color when inherit)
  var settingColors: [String: String] = [:]  // extra color settings → "inherit"/"#rrggbb"
  // Present when this rep has per-atom transparency overriding the object-level
  // slider; carries the transparency setting name and the effective min–max range.
  var atomTransp: AtomTransp? = nil
}

/// Effective per-atom transparency range for a rep whose object-level slider is
/// overridden by atom-level settings (from appkit_inspector's `atom_transp`).
struct AtomTransp: Equatable {
  let setting: String  // e.g. "cartoon_transparency"
  let min: Double
  let max: Double
}

/// Global "Scene" parameters.
struct SceneState: Equatable {
  var values: [String: Double] = [:]  // setting name → value (toggles 0/1)
  var outlineColor: [Double] = [0, 0, 0]  // metal_outline_color r,g,b in 0…1
}

/// Per-object state metadata for the inspector STATE row (multi-state objects).
struct ObjStateMeta: Equatable {
  var state: Int = 1  // effective current state (resolves the frame)
  var overlayAll: Bool = false  // all_states overlay for this object
  // Per-state titles (e.g. compound names from a multi-record SDF). Empty
  // unless at least one state carries a title. Indexed by state-1 (issue #203).
  var titles: [String] = []

  /// Title for a 1-based state, or nil when none/blank.
  func title(forState state: Int) -> String? {
    guard state >= 1, state <= titles.count else { return nil }
    let t = titles[state - 1]
    return t.isEmpty ? nil : t
  }
}

// MARK: - Data Models

/// Represents a PyMOL object or selection in the panel
struct ObjectEntry: Identifiable, Equatable {
  let id: String
  let name: String
  var isEnabled: Bool
  var isSelection: Bool
  var atomCount: Int?
  // Number of coordinate states (NMR models / trajectory frames). >1 surfaces
  // the per-object STATE controls in the inspector. Defaults to 1.
  var stateCount: Int = 1
  // Group tree (#255). `isGroup` marks an `object:group`; `parent` is the name
  // of the group this row lives in (nil = top level). Both default so the
  // existing memberwise initializers keep compiling.
  var isGroup: Bool = false
  var parent: String? = nil
  var isMap: Bool = false

  var displayName: String {
    if isSelection, let count = atomCount {
      return "\(name) (\(count))"
    }
    return name
  }
}

// MARK: - Group tree (#255)

/// One drawn row of the flattened object tree.
private struct TreeRow: Identifiable {
  let entry: ObjectEntry
  let depth: Int
  var id: String { entry.id }
}

/// Three-way visibility state. A group summarises its descendants, so it can be
/// partially on in a way a single object never is.
private enum CheckState { case on, off, mixed }

/// Flatten the object list into the rows the panel draws, parents before children.
///
/// `engine.objects` mirrors PyMOL's own ordering and must NOT be reordered in
/// place — `parseObjectPanelFeedback`'s equality guard and `setObjectEnabled`'s
/// index lookup both assume it does. So nesting is derived here, at render time.
///
/// Children of a *collapsed* group are omitted entirely, which is what makes the
/// disclosure a real collapse rather than a styling choice. Two robustness rules:
/// an entry whose parent is not in the list is drawn at top level, and anything
/// left unreachable (a cycle from a corrupt parent map) is drawn at top level too
/// — an object must never silently disappear from the panel.
private func flattenObjectTree(
  _ objects: [ObjectEntry],
  openGroups: Set<String>
) -> [TreeRow] {
  let known = Set(objects.map { $0.name })
  let byParent = Dictionary(
    grouping: objects.filter {
      if let p = $0.parent { return known.contains(p) } else { return false }
    }, by: { $0.parent! })

  let roots = objects.filter { $0.parent == nil || !known.contains($0.parent!) }

  // Reachability ignores open/closed, so a collapsed child is not mistaken for
  // an orphan and re-emitted at top level.
  var reachable = Set<String>()
  func mark(_ name: String) {
    guard !reachable.contains(name) else { return }
    reachable.insert(name)
    for c in byParent[name] ?? [] { mark(c.name) }
  }
  for r in roots { mark(r.name) }

  var rows: [TreeRow] = []
  func emit(_ entry: ObjectEntry, _ depth: Int) {
    rows.append(TreeRow(entry: entry, depth: depth))
    guard entry.isGroup, openGroups.contains(entry.name) else { return }
    for child in byParent[entry.name] ?? [] { emit(child, depth + 1) }
  }
  for r in roots { emit(r, 0) }
  for entry in objects where !reachable.contains(entry.name) { emit(entry, 0) }
  return rows
}

/// True when every ancestor of `entry` is enabled.
///
/// A member keeps reporting `enabled` while it is hidden by a disabled parent —
/// `disable grp` leaves `get_names('public_objects', enabled_only=1)` returning
/// the members — so a row must not take `isEnabled` at face value or it will show
/// a checked box for something invisible.
private func ancestorsEnabled(_ entry: ObjectEntry, byName: [String: ObjectEntry]) -> Bool {
  var seen = Set<String>()
  var next = entry.parent
  while let name = next, !seen.contains(name), let anc = byName[name] {
    if !anc.isEnabled { return false }
    seen.insert(name)
    next = anc.parent
  }
  return true
}

/// Checkbox state for a row; groups fold in their whole subtree.
private func checkState(
  for entry: ObjectEntry,
  byParent: [String: [ObjectEntry]]
) -> CheckState {
  guard entry.isGroup else { return entry.isEnabled ? .on : .off }
  if !entry.isEnabled { return .off }
  var anyOn = false
  var anyOff = false
  var seen = Set<String>()
  var stack = byParent[entry.name] ?? []
  while let c = stack.popLast() {
    guard !seen.contains(c.name) else { continue }
    seen.insert(c.name)
    if c.isEnabled { anyOn = true } else { anyOff = true }
    stack.append(contentsOf: byParent[c.name] ?? [])
  }
  return (anyOn && anyOff) ? .mixed : (anyOff ? .off : .on)
}

// MARK: - Menu Option Definitions

/// Show/Hide representation options
private let showHideOptions: [(label: String, rep: String?)] = [
  ("everything", "everything"),
  ("---", nil),
  ("lines", "lines"),
  ("nonbonded", "nonbonded"),
  ("---", nil),
  ("sticks", "sticks"),
  ("nb_spheres", "nb_spheres"),
  ("---", nil),
  ("ribbon", "ribbon"),
  ("cartoon", "cartoon"),
  ("labels", "labels"),
  ("cell", "cell"),
  ("dots", "dots"),
  ("spheres", "spheres"),
  ("mesh", "mesh"),
  ("surface", "surface"),
  ("volume", "volume"),
  ("slice", "slice"),
  ("extent", "extent"),
  ("---", nil),
  ("licorice", "licorice"),
  ("wire", "wire"),
  ("dashes", "dashes"),
]

/// Label options
private let labelOptions: [(label: String, expr: String?)] = [
  ("None", ""),
  ("---", nil),
  ("Residues", "resn+resi"),
  ("Chains", "chain"),
  ("Segments", "segi"),
  ("Atoms", "name"),
  ("Elements", "elem"),
]

// MARK: - Action Menu Structure

/// Hierarchical action menu item.
///
/// Internal rather than private so the unit tests can assert which rows a given
/// kind of object gets — menu *composition* is the behaviour in #461/#468, and
/// it is otherwise only reachable by driving SwiftUI.
indirect enum ActionMenuItem {
  case action(label: String, key: String)
  case separator
  case submenu(label: String, children: [ActionMenuItem])
  // Dynamic submenus whose children are the currently-loaded molecule
  // objects / public selections (built at render time from engine.objects),
  // mirroring desktop PyMOL's "to molecule (*/CA)" / "to selection (*/CA)"
  // align entries.
  case alignToMolecule(label: String)
  case alignToSelection(label: String)
  /// Dynamic submenu listing existing groups plus "New Group…" (#255). Built at
  /// render time from engine.objects, like the align cases, because the set of
  /// groups changes with every poll.
  case moveToGroup(label: String)
  /// Dynamic submenu listing "New Object…" plus every loaded molecule object,
  /// mirroring desktop PyMOL's selection-menu "copy to object ▸" (#461). Built
  /// at render time from engine.objects for the same reason as the cases above.
  case copyToObject(label: String)
}

private let baseActionMenuItems: [ActionMenuItem] = [
  .action(label: "Zoom", key: "zoom"),
  .action(label: "Orient", key: "orient"),
  .action(label: "Center", key: "center"),
  .action(label: "Origin", key: "origin"),
  .separator,
  .action(label: "Drag Matrix", key: "drag_matrix"),
  .action(label: "Reset Matrix", key: "reset_matrix"),
  .separator,
  .action(label: "Drag Coordinates", key: "drag_coords"),
  .action(label: "Clean", key: "clean"),
  .separator,
  .submenu(
    label: "Preset",
    children: [
      .action(label: "classified", key: "preset_classified"),
      .separator,
      .action(label: "simple", key: "preset_simple"),
      .action(label: "simple (no solvent)", key: "preset_simple_no_solv"),
      .action(label: "ball and stick", key: "preset_ball_and_stick"),
      .action(label: "b factor putty", key: "preset_b_factor_putty"),
      .action(label: "technical", key: "preset_technical"),
      .action(label: "ligands", key: "preset_ligands"),
      .action(label: "pretty", key: "preset_pretty"),
      .action(label: "pretty (with solvent)", key: "preset_pretty_solv"),
      .action(label: "publication", key: "preset_publication"),
      .action(label: "publication (with solvent)", key: "preset_pub_solv"),
      .separator,
      .action(label: "protein interface", key: "preset_interface"),
      .separator,
      .action(label: "default", key: "preset_default"),
    ]),
  .submenu(
    label: "Find",
    children: [
      .action(label: "polar contacts (within)", key: "find_polar_within"),
      .action(label: "polar contacts (to other)", key: "find_polar_other"),
      .action(label: "polar contacts (any)", key: "find_polar_any"),
      .separator,
      .action(label: "halogen bonds", key: "find_halogen_bond"),
      .action(label: "salt bridges", key: "find_salt_bridge"),
      .separator,
      .action(label: "pi interactions (all)", key: "find_pi_all"),
      .action(label: "pi-pi", key: "find_pi_pi"),
      .action(label: "pi-cation", key: "find_pi_cation"),
    ]),
  .action(label: "Contact Map…", key: "contact_map"),
  .action(label: "Ramachandran Plot…", key: "ramachandran"),
  .submenu(
    label: "Align",
    children: [
      .action(label: "Align Wizard…", key: "align_wizard"),
      .separator,
      .alignToMolecule(label: "to molecule (*/CA)"),
      .alignToSelection(label: "to selection (*/CA)"),
      .separator,
      .action(label: "enabled to this (*/CA)", key: "align_enabled"),
      .action(label: "all to this (*/CA)", key: "align_all"),
      .separator,
      .action(label: "states (*/CA)", key: "align_states_ca"),
      .action(label: "states", key: "align_states"),
      .separator,
      .action(label: "matrix reset", key: "matrix_reset"),
    ]),
  .submenu(
    label: "Generate",
    children: [
      .action(label: "vacuum electrostatics", key: "gen_vacuum_esp"),
      .separator,
      .action(label: "symmetry mates 4 A", key: "gen_symm_4"),
      .action(label: "symmetry mates 8 A", key: "gen_symm_8"),
      .action(label: "symmetry mates 20 A", key: "gen_symm_20"),
      .separator,
      .action(label: "split chains", key: "gen_split_chains"),
      .action(label: "split states", key: "state_split"),
    ]),
  .separator,
  .action(label: "Assign Sec. Struc.", key: "dss"),
  .separator,
  .submenu(
    label: "Hydrogens",
    children: [
      .action(label: "hide", key: "h_hide"),
      .action(label: "show", key: "h_show"),
      .separator,
      .action(label: "add", key: "h_add"),
      .action(label: "add polar", key: "h_add_polar"),
      .separator,
      .action(label: "remove", key: "h_remove"),
      .action(label: "remove nonpolar", key: "h_remove_nonpolar"),
    ]),
  .action(label: "Remove Waters", key: "remove_waters"),
  .separator,
  .submenu(
    label: "State",
    children: [
      .action(label: "freeze", key: "state_freeze"),
      .action(label: "all states", key: "state_all"),
      .action(label: "thaw", key: "state_thaw"),
      .separator,
      .action(label: "split", key: "state_split"),
      .separator,
      .action(label: "load trajectory…", key: "load_traj"),
    ]),
  .submenu(
    label: "Sequence",
    children: [
      .action(label: "include", key: "seq_include"),
      .action(label: "exclude", key: "seq_exclude"),
      .action(label: "default", key: "seq_default"),
    ]),
  .submenu(
    label: "Movement",
    children: [
      .action(label: "protect", key: "movement_protect"),
      .action(label: "deprotect", key: "movement_deprotect"),
    ]),
  .submenu(
    label: "Masking",
    children: [
      .action(label: "mask", key: "masking_mask"),
      .action(label: "unmask", key: "masking_unmask"),
    ]),
  .submenu(
    label: "Compute",
    children: [
      .action(label: "atom count", key: "compute_count"),
      .separator,
      .action(label: "formal charge sum", key: "compute_formal_charge"),
      .action(label: "partial charge sum", key: "compute_partial_charge"),
      .separator,
      .action(label: "molecular surface area", key: "compute_mol_area"),
      .action(label: "solvent accessible area", key: "compute_sasa"),
      .separator,
      .action(label: "mol. weight (explicit)", key: "compute_mass_explicit"),
      .action(label: "mol. weight (with H)", key: "compute_mass_implicit"),
    ]),
  .separator,
  .moveToGroup(label: "Move to Group"),
  .separator,
  .action(label: "Rename", key: "rename"),
  .action(label: "Duplicate", key: "copy"),
  .action(label: "Delete", key: "delete"),
]

/// The object/selection action menu, tailored to the row it's shown on.
///
/// Selections additionally get **Remove Atoms** — PyMOL's selection-menu
/// `remove atoms` (`cmd.remove(sele); cmd.delete(sele)`): it strips the matched
/// atoms *out of their parent object* and then drops the now-empty selection.
/// This is the quick-cleanup path from issue #202 (strip waters, delete an
/// unwanted chain, trim to a binding site) — distinct from "Delete", which only
/// removes the selection marker and leaves every atom in place. Objects don't
/// get it (they keep Delete + Remove Waters), matching desktop PyMOL.
/// Action menu for a GROUP row (#255).
///
/// Deliberately much shorter than `baseActionMenuItems`: most of that menu is
/// atom-level (Preset, Find, Compute, dss) and would be wrong or slow applied to
/// a whole subtree. What remains is camera framing, the visibility verbs, and
/// group lifecycle.
///
/// The two destructive entries are kept apart on purpose. "Ungroup" dissolves the
/// group and leaves every object at top level; "Delete Group + Contents" removes
/// the members too, because `delete <group>` in PyMOL takes the whole subtree with
/// it. Labelling them identically would make an irreversible action look routine.
/// Action menu for a map object (2fo-fc, cryo-EM…): contour it as a mesh or a
/// surface. Levels are in σ for normalized maps (PyMOL normalizes CCP4 on load).
private let mapActionMenuItems: [ActionMenuItem] = [
  .action(label: "Zoom", key: "zoom"),
  .action(label: "Center", key: "center"),
  .separator,
  .submenu(
    label: "Mesh",
    children: ["1.0", "1.5", "2.0", "3.0"].map {
      .action(label: "@ \($0) σ", key: "map_mesh_" + $0)
    } + [
      .separator,
      .action(label: "@ 1.0 σ around sele (2 Å)", key: "map_mesh_sele"),
    ]),
  .submenu(
    label: "Surface",
    children: ["1.0", "1.5", "2.0", "3.0"].map {
      .action(label: "@ \($0) σ", key: "map_surface_" + $0)
    }),
  .separator,
  .action(label: "Rename", key: "rename"),
  .action(label: "Delete", key: "delete"),
]

private let groupActionMenuItems: [ActionMenuItem] = [
  .action(label: "Zoom", key: "zoom"),
  .action(label: "Orient", key: "orient"),
  .action(label: "Center", key: "center"),
  .separator,
  .action(label: "Inspect Group…", key: "group_inspect"),
  .separator,
  .action(label: "Enable All Members", key: "group_enable_members"),
  .action(label: "Disable All Members", key: "group_disable_members"),
  .separator,
  .action(label: "Rename…", key: "rename"),
  .action(label: "Ungroup (keep objects)", key: "group_ungroup"),
  .action(label: "Delete Group + Contents", key: "group_delete"),
]

/// Copy / Extract section, shown on SELECTION rows only (#461, #468).
///
/// Both verbs make a new object out of the selected atoms; they differ in what
/// happens to the source. **Copy to Object** leaves the source intact (upstream
/// `cmd.copy_to` / `cmd.create`), **Extract** moves the atoms out of it
/// (`cmd.extract`). Desktop PyMOL offers both on selections and neither on
/// objects — `modules/pymol/menu.py:1183` (sele_action), :1209 (sele_action2),
/// :1777 and :1836 all sit in the `else:` branch of an `if object:` — so the
/// section is gated the same way here.
///
/// The three Extract variants come from PyMOL's own Extract submenu
/// (`modules/pymol/menu.py:47`): the selection alone, or widened by one bond
/// shell, atom-wise or residue-wise.
private let copyExtractMenuItems: [ActionMenuItem] = [
  .copyToObject(label: "Copy to Object"),
  .submenu(
    label: "Extract",
    children: [
      .action(label: "object", key: "extract_object"),
      .action(label: "extend 1", key: "extract_extend_1"),
      .action(label: "byres extend 1", key: "extract_byres_extend_1"),
    ]),
]

func actionMenuItems(isSelection: Bool, isGroup: Bool = false, isMap: Bool = false)
  -> [ActionMenuItem]
{
  if isGroup { return groupActionMenuItems }
  if isMap { return mapActionMenuItems }
  guard isSelection else { return baseActionMenuItems }
  var items = baseActionMenuItems
  // Place "Remove Atoms" in the cleanup section right after "Remove Waters" —
  // both strip atoms out of the parent object, so they belong together next to
  // Hydrogens. This is distinct from the trailing "Delete", which only drops
  // the selection marker and leaves every atom in place.
  let removeAtoms = ActionMenuItem.action(label: "Remove Atoms", key: "remove_atoms")
  if let watersIdx = items.firstIndex(where: {
    if case .action(_, "remove_waters") = $0 { return true }
    return false
  }) {
    items.insert(removeAtoms, at: watersIdx + 1)
  } else {
    items.append(removeAtoms)
  }
  // Copy / Extract goes in the object-management tail, just above "Move to
  // Group" — next to Rename/Duplicate/Delete, which is where upstream keeps it
  // too ("duplicate, copy to object, extract object" in sele_action). It rides
  // on the separator already sitting before moveToGroup and adds its own after,
  // so the tail reads: … Compute │ Copy/Extract │ Move to Group │ Rename….
  if let groupIdx = items.firstIndex(where: {
    if case .moveToGroup = $0 { return true }
    return false
  }) {
    items.insert(contentsOf: copyExtractMenuItems + [.separator], at: groupIdx)
  } else {
    items.append(contentsOf: [.separator] + copyExtractMenuItems)
  }
  return items
}

// MARK: - Command Dispatch

/// A ▸ State ▸ load trajectory…: append a trajectory's frames as states.
private func loadTrajectory(onto name: String, engine: PyMOLEngine) {
  let panel = NSOpenPanel()
  panel.message = "Load a trajectory onto \(name) (DCD, TRR, XTC, NetCDF, PDB…)"
  guard panel.runModal() == .OK, let url = panel.url else { return }
  engine.runCommand("load_traj \(url.path), \(name)")
}

/// Translates action keys into PyMOL commands and runs them
private func runActionCommand(_ key: String, name: String, engine: PyMOLEngine) {
  // These open a tab / sheet rather than running a command.
  if key == "contact_map" {
    engine.contactMapRequest = name
    return
  }
  if key == "ramachandran" {
    engine.contactMapRequest = "rama:" + name
    return
  }
  if key == "align_wizard" {
    engine.alignMobile = name
    engine.openTool(.align)
    return
  }
  if key == "load_traj" {
    loadTrajectory(onto: name, engine: engine)
    return
  }
  if key == "new_selection" {
    engine.showSelectionBuilder = true
    return
  }
  let n = name  // shorthand
  let cmd: String
  switch key {
  // View / Transform
  case "zoom": cmd = "zoom \(n), animate=-1"
  case "orient": cmd = "orient \(n), animate=-1"
  case "center": cmd = "center \(n), animate=-1"
  case "origin": cmd = "origin \(n)"
  case "drag_matrix": cmd = "drag \(n)"
  case "reset_matrix": cmd = "reset object=\(n)"
  case "drag_coords": cmd = "drag (\(n))"
  case "clean": cmd = "clean \(n)"
  case "dss": cmd = "dss \(n)"
  // Presets: the key's suffix is the pymol.preset function name.
  case let k where k.hasPrefix("preset_"):
    cmd =
      "python\nfrom pymol import preset; preset.\(k.dropFirst(7))('\(n)', _self=cmd)\npython end"
  // Find
  case "find_polar_within":
    cmd =
      "dist \(n)_polar_conts, \(n), \(n), quiet=1, mode=2, label=0, reset=1; enable \(n)_polar_conts"
  case "find_polar_other":
    cmd =
      "dist \(n)_polar_conts, (\(n)), (byobj (\(n))) and (not (\(n))), quiet=1, mode=2, label=0, reset=1; enable \(n)_polar_conts"
  case "find_polar_any":
    cmd =
      "dist \(n)_polar_conts, (\(n)), (not \(n)), quiet=1, mode=2, label=0, reset=1; enable \(n)_polar_conts"
  case "find_halogen_bond": cmd = "distance \(n)_halogen_bond, \(n), same, reset=1, mode=9"
  case "find_salt_bridge": cmd = "distance \(n)_salt_bridge, \(n), same, reset=1, mode=10"
  case "find_pi_all": cmd = "pi_interactions \(n)_pi_interactions, \(n), reset=1"
  case "find_pi_pi": cmd = "distance \(n)_pi_pi, \(n), same, reset=1, mode=6"
  case "find_pi_cation": cmd = "distance \(n)_pi_cation, \(n), same, reset=1, mode=7"
  // Align
  case "align_enabled": cmd = "python\ncmd.util.mass_align('\(n)', 1, _self=cmd)\npython end"
  case "align_all": cmd = "python\ncmd.util.mass_align('\(n)', 0, _self=cmd)\npython end"
  case "align_states_ca": cmd = "intra_fit (\(n)) and name CA"
  case "align_states": cmd = "intra_fit \(n)"
  case "matrix_reset": cmd = "matrix_reset \(n)"
  // Generate
  case "gen_vacuum_esp":
    cmd = "python\ncmd.util.protein_vacuum_esp('\(n)', mode=2, quiet=0, _self=cmd)\npython end"
  case let k where k.hasPrefix("gen_symm_"):
    cmd = "symexp \(n)_, \(n), \(n), cutoff=\(k.dropFirst(9)), segi=1"
  case "gen_split_chains": cmd = "split_chains \(n)"
  case "map_mesh_sele": cmd = "isomesh \(n)_mesh, \(n), 1.0, sele, carve=2"
  case let k where k.hasPrefix("map_mesh_"):
    cmd = "isomesh \(n)_mesh, \(n), \(k.dropFirst(9))"
  case let k where k.hasPrefix("map_surface_"):
    cmd = "isosurface \(n)_surf, \(n), \(k.dropFirst(12))"
  // Hydrogens
  case "h_hide": cmd = "hide everything, (\(n)) and hydro"
  case "h_show": cmd = "show sticks, (\(n)) and hydro"
  case "h_add": cmd = "h_add \(n); sort \(n) extend 1"
  case "h_add_polar": cmd = "h_add \(n) & (don.|acc.); sort \(n) extend 1"
  case "h_remove": cmd = "remove (\(n)) and hydro"
  case "h_remove_nonpolar": cmd = "remove \(n) & hydro & not nbr. (don.|acc.)"
  case "remove_waters": cmd = "remove (solvent and (\(n)))"
  // State
  case "state_freeze": cmd = "python\ncmd.set('state', cmd.get_state(), '\(n)')\npython end"
  case "state_all": cmd = "set state, 0, \(n)"
  case "state_thaw": cmd = "unset all_states, \(n); unset state, \(n)"
  case "state_split": cmd = "split_states \(n)"
  // Sequence
  case "seq_include": cmd = "set seq_view, on, \(n)"
  case "seq_exclude": cmd = "set seq_view, off, \(n)"
  case "seq_default": cmd = "unset seq_view, \(n)"
  // Movement
  case "movement_protect": cmd = "protect \(n)"
  case "movement_deprotect": cmd = "deprotect \(n)"
  // Masking
  case "masking_mask": cmd = "mask \(n)"
  case "masking_unmask": cmd = "unmask \(n)"
  // Compute
  case "compute_count": cmd = "count_atoms \(n)"
  case "compute_formal_charge":
    cmd = "python\ncmd.util.sum_formal_charges('\(n)', quiet=0, _self=cmd)\npython end"
  case "compute_partial_charge":
    cmd = "python\ncmd.util.sum_partial_charges('\(n)', quiet=0, _self=cmd)\npython end"
  case "compute_mol_area":
    cmd = "python\ncmd.util.get_area('\(n)', -1, 0, quiet=0, _self=cmd)\npython end"
  case "compute_sasa": cmd = "python\ncmd.util.get_sasa('\(n)', quiet=0, _self=cmd)\npython end"
  case "compute_mass_explicit":
    cmd = "python\ncmd.util.compute_mass('\(n)', implicit=False, quiet=0, _self=cmd)\npython end"
  case "compute_mass_implicit":
    cmd = "python\ncmd.util.compute_mass('\(n)', implicit=True, quiet=0, _self=cmd)\npython end"
  // Object management
  case "rename":
    // PyMOL's `wizard renaming` has no UI here — request a name-entry modal
    // (presented by ObjectPanel) instead.
    engine.pendingRename = n
    return
  case "copy": cmd = "copy \(n)_copy, \(n)"
  case "delete": cmd = "delete \(n)"
  case "delete_everything": cmd = "delete all"
  // Group lifecycle (#255). `action=excise` dissolves the group and leaves its
  // members at top level; plain `delete <group>` removes the members too, which
  // is why the two are separate menu entries with explicit labels.
  case "group_ungroup": cmd = "group \(n), , action=excise"
  case "group_delete": cmd = "delete \(n)"
  case "group_enable_members": cmd = "enable \(n)"
  case "group_disable_members": cmd = "disable \(n)"
  case "group_inspect":
    // A group's rep card is the union of its members' reps (see #256). The
    // chevron on a group row is spoken for by child expansion, so this is how
    // the card is opened.
    engine.expandedDetail = n
    return
  case "group_new":
    // Create a group holding just this object; PyMOL creates the group on
    // first reference, so no separate create step is needed.
    engine.pendingGroupFor = n
    return
  // Strip the selection's atoms out of their parent object(s), then drop the
  // now-empty selection (PyMOL's selection-menu "remove atoms"). Shown only on
  // selection rows — see actionMenuItems(isSelection:).
  case "remove_atoms": cmd = "remove (\(n)); delete \(n)"
  // Extract the selection into a new auto-named object, moving its atoms out
  // of the parent (#468). Shown only on selection rows, like "Remove Atoms".
  case "extract_object": cmd = extractCommand(sele: n, scope: .object)
  case "extract_extend_1": cmd = extractCommand(sele: n, scope: .extend1)
  case "extract_byres_extend_1": cmd = extractCommand(sele: n, scope: .byresExtend1)
  // Global ("all" row) actions
  case "deselect": cmd = "deselect"
  case "hide_everything": cmd = "hide everything, \(n)"
  case "reset_view": cmd = "reset"
  // Moves the objects, not the camera: every object is shifted so their
  // centers coincide, so switching between unrelated PDBs stops requiring a
  // re-zoom. Non-destructive (object matrix) — `matrix_reset` puts it back.
  // The follow-up zoom frames the stacked result — the camera may be aimed
  // at where things used to be.
  case "center_all": cmd = "center_all\nzoom visible"
  case "center_all_reset":
    cmd =
      "python\nfor _o in cmd.get_names('public_objects'):\n    cmd.matrix_reset(_o, mode=1)\npython end"
  default: return
  }
  engine.runCommand(cmd, naming: name)
}

/// Action ("A") menu for the global "all" row — a focused, scene-wide subset.
/// Reuses the per-object action keys (all valid with name "all"); deliberately
/// omits per-object items (Rename / Duplicate / Delete) and adds global ones
/// (Deselect, Hide everything, Reset camera).
private let allActionMenuItems: [ActionMenuItem] = [
  .action(label: "Zoom", key: "zoom"),
  .action(label: "Orient", key: "orient"),
  .action(label: "Center camera", key: "center"),
  .action(label: "Reset camera", key: "reset_view"),
  .separator,
  // "Center all" moves the objects onto a shared center; "Center camera"
  // above only aims the camera. The labels have to carry that distinction —
  // they sit two rows apart in the same menu.
  .submenu(
    label: "Center all objects",
    children: [
      .action(label: "stack on first object", key: "center_all"),
      .separator,
      .action(label: "undo (reset matrices)", key: "center_all_reset"),
    ]),
  .separator,
  .action(label: "Deselect", key: "deselect"),
  .action(label: "New Selection…", key: "new_selection"),
  .action(label: "Hide everything", key: "hide_everything"),
  .separator,
  .action(label: "Assign Sec. Struc.", key: "dss"),
  .action(label: "Remove Waters", key: "remove_waters"),
  .submenu(
    label: "Hydrogens",
    children: [
      .action(label: "add", key: "h_add"),
      .action(label: "add polar", key: "h_add_polar"),
      .action(label: "remove", key: "h_remove"),
    ]),
  .submenu(
    label: "Find",
    children: [
      .action(label: "polar contacts (any)", key: "find_polar_any"),
      .action(label: "salt bridges", key: "find_salt_bridge"),
    ]),
  .submenu(
    label: "Preset",
    children: [
      .action(label: "pretty", key: "preset_pretty"),
      .action(label: "technical", key: "preset_technical"),
      .action(label: "ball and stick", key: "preset_ball_and_stick"),
      .action(label: "default", key: "preset_default"),
    ]),
  .separator,
  // Desktop PyMOL's all ▸ A ▸ delete everything.
  .action(label: "Delete everything", key: "delete_everything"),
]

// MARK: - Theme

// The classic chrome's neutrals, blended from the panel grey toward the panel
// text so they stay opaque.
enum PanelTheme {
  private static let bg = ClassicUI.chromeBg, fg = ClassicUI.chromeText
  static let rowAltBackground = bg.blended(with: fg, 0.03).color
  static let textColor = fg.color
  static let buttonBackground = bg.blended(with: fg, 0.16).color
  static let accentColor = ClassicUI.accent
  static let headerColor = bg.blended(with: fg, 0.6).color
  static let disabledColor = bg.blended(with: fg, 0.4).color
  // Amber accent for the per-atom transparency badge / detail row — a fixed hue
  // that reads as "transparency".
  static let atomTranspColor = Color(.sRGB, red: 0.90, green: 0.66, blue: 0.30, opacity: 1)
}

// Row sizing: desktop PyMOL's internal GUI is dense
// (ExecToggleSize 16 on an 18px internal_gui_control_size line).
private let kClassicBtnW: CGFloat = 18
private let kClassicBtnH: CGFloat = 19
private let kClassicRowH: CGFloat = 21

/// The face of one A/S/H/L/C menu button: the internal GUI's bevelled toggles
/// (black letter; A/L blue, S lavender, H slate, C four-corner rainbow).
private struct ASHLCFace: View {
  let letter: String

  var body: some View {
    Text(letter)
      .font(ClassicUI.font)
      .foregroundColor(.black)
      .frame(width: kClassicBtnW, height: kClassicBtnH)
      .background(classicFill)
      .contentShape(Rectangle())
  }

  @ViewBuilder private var classicFill: some View {
    if letter == "C" {
      ClassicBevel(
        fill: ClassicRainbowFill(),
        light: ClassicUI.toggleLight, dark: ClassicUI.toggleDark)
    } else {
      ClassicBevel(
        fill: ClassicUI.toggleFill(letter),
        light: ClassicUI.toggleLight, dark: ClassicUI.toggleDark)
    }
  }
}

/// The internal GUI's name button. The bevel's fill is
/// the enable state (0.5 grey on, 0.25 off), which is why classic rows carry no
/// checkbox; clicking the name toggles it, as in desktop PyMOL.
private struct ClassicNameCell: View {
  let text: String
  let enabled: Bool
  var bold: Bool = false
  let toggle: () -> Void

  var body: some View {
    Button(action: toggle) {
      Text(text)
        .font(ClassicUI.font.weight(bold ? .semibold : .regular))
        .foregroundColor(.white)
        .lineLimit(1)
        .truncationMode(.tail)
        .padding(.horizontal, 3)
        .frame(maxWidth: .infinity, maxHeight: kClassicBtnH, alignment: .leading)
        .background(
          ClassicBevel(fill: enabled ? ClassicUI.nameEnabled : ClassicUI.nameDisabled)
        )
        .contentShape(Rectangle())
    }
    .buttonStyle(ClassicPressStyle())
    .help(text)
  }
}

/// Row chrome shared by every A/S/H/L/C row: flat and dense on the panel's grey.
extension View {
  fileprivate func objectRowChrome() -> some View {
    self.padding(.horizontal, 2)
      .frame(height: kClassicRowH)
      .background(ClassicUI.panelBg)
  }
}

// MARK: - ObjectPanel View

struct ObjectPanel: View {
  @EnvironmentObject var engine: PyMOLEngine
  @State private var renameText = ""
  @State private var groupNameText = ""
  @State private var copyToNewText = ""
  // Which group rows are expanded. Deliberately NOT engine.expandedDetail: that
  // is a single-slot accordion for the rep inspector whose didSet re-polls and
  // which ContentView reads to resize the panel, so reusing it would make
  // opening a group collapse the inspector and resize the pane. Held here rather
  // than as @State on the row because rows live in a LazyVStack, which discards
  // offscreen state. PyMOL's own open/closed flag (SpecRec OpenOrClosed) has no
  // Python accessor, so this is the source of truth and we push it down with
  // `group ..., action=open|close` — see the note in toggleGroupOpen.
  @State private var openGroups: Set<String> = []

  // The panel's name-entry modals are split across two computed properties
  // rather than one `.alert` chain: three of them in a single expression put
  // the body past the Swift type-checker's budget ("unable to type-check in
  // reasonable time"), which is the same reason other views in this file are
  // broken up.
  var body: some View {
    panelWithNamingAlerts
      // "Copy to Object ▸ New Object…" (#461) — same request-channel shape
      // as the two below. Unlike Extract, which auto-names the way desktop
      // PyMOL does, the copy asks: a copy is usually made to be kept and
      // compared against its source, and "obj01" tells you nothing then.
      .alert(
        "Copy “\(engine.pendingCopyToNew ?? "")” to a new object",
        isPresented: Binding(
          get: { engine.pendingCopyToNew != nil },
          set: { if !$0 { engine.pendingCopyToNew = nil } })
      ) {
        TextField("Object name", text: $copyToNewText)
        Button("Copy") {
          if let sele = engine.pendingCopyToNew {
            let new = copyToNewText.trimmingCharacters(in: .whitespaces)
            // Guard in the action body, not only via .disabled() on
            // the button: this is the path that actually builds a
            // command string, so it is where a bad name has to stop.
            if canNameNewObject(new, existing: engine.objects.map(\.name)) {
              engine.runCommand(copyToNewObjectCommand(sele: sele, name: new), naming: sele, new)
            }
          }
          engine.pendingCopyToNew = nil
        }
        .disabled(
          !canNameNewObject(
            copyToNewText.trimmingCharacters(in: .whitespaces),
            existing: engine.objects.map(\.name)))
        Button("Cancel", role: .cancel) { engine.pendingCopyToNew = nil }
      } message: {
        Text(copyToNewMessage)
      }
      .onChange(of: engine.pendingCopyToNew) { newValue in
        if newValue != nil {
          copyToNewText = defaultNewObjectName(existing: engine.objects.map(\.name))
        }
      }
  }

  /// The "Copy to Object ▸ New Object…" alert's subtitle — the reassurance when
  /// the name is usable, and the reason the Copy button is greyed when it is
  /// not. Says *why* rather than just refusing, since both failure modes
  /// (illegal character, name already taken) are invisible in the field itself.
  private var copyToNewMessage: String {
    let name = copyToNewText.trimmingCharacters(in: .whitespaces)
    if name.isEmpty { return "Enter a name for the new object." }
    if !isLegalObjectName(name) {
      return "PyMOL names can use letters, digits and + - . ^ _ only."
    }
    if engine.objects.contains(where: { $0.name == name }) {
      return "“\(name)” is already taken. Pick a name no object or selection is using."
    }
    return "The selected atoms stay in their current object too."
  }

  private var panelWithNamingAlerts: some View {
    panelBody
      // Name-entry modal for the action-menu "Rename" (engine.pendingRename).
      .alert(
        "Rename “\(engine.pendingRename ?? "")”",
        isPresented: Binding(
          get: { engine.pendingRename != nil },
          set: { if !$0 { engine.pendingRename = nil } })
      ) {
        TextField("New name", text: $renameText)
        Button("Rename") {
          if let old = engine.pendingRename {
            let new = renameText.trimmingCharacters(in: .whitespaces)
            if !new.isEmpty && new != old { engine.renameObject(old, to: new) }
          }
          engine.pendingRename = nil
        }
        Button("Cancel", role: .cancel) { engine.pendingRename = nil }
      } message: {
        Text("Enter a new name for this object.")
      }
      .onChange(of: engine.pendingRename) { newValue in
        if let n = newValue { renameText = n }  // prefill with current name
      }
      // "New Group…" (#255) — same request-channel shape as the rename alert.
      .alert(
        "New group for “\(engine.pendingGroupFor ?? "")”",
        isPresented: Binding(
          get: { engine.pendingGroupFor != nil },
          set: { if !$0 { engine.pendingGroupFor = nil } })
      ) {
        TextField("Group name", text: $groupNameText)
        Button("Create") {
          if let member = engine.pendingGroupFor {
            let g = groupNameText.trimmingCharacters(in: .whitespaces)
            if !g.isEmpty && g != member {
              // PyMOL creates the group on first reference, so this
              // both creates it and moves the object in.
              engine.runCommand("group \(g), \(member)", naming: member, tokens: [g])
              openGroups.insert(g)  // show the result immediately
            }
          }
          engine.pendingGroupFor = nil
        }
        Button("Cancel", role: .cancel) { engine.pendingGroupFor = nil }
      } message: {
        Text("Enter a name for the new group.")
      }
      .onChange(of: engine.pendingGroupFor) { newValue in
        if newValue != nil { groupNameText = "" }
      }
  }

  // Desktop PyMOL's flat list: "all", the objects, then the selections.
  private var panelBody: some View {
    ScrollView(.vertical) {
      LazyVStack(spacing: 0) {
        let objects = engine.objects.filter { !$0.isSelection }
        let selections = engine.objects.filter { $0.isSelection }

        AllControlsRow()
        // Nested rows come from a flattened tree walk.
        let rows = flattenObjectTree(objects, openGroups: openGroups)
        let byName = Dictionary(
          objects.map { ($0.name, $0) },
          uniquingKeysWith: { a, _ in a })
        let byParent = Dictionary(
          grouping: objects.filter { $0.parent != nil },
          by: { $0.parent! })
        ForEach(rows, id: \.id) { row in
          ObjectCard(
            entry: row.entry,
            depth: row.depth,
            check: checkState(for: row.entry, byParent: byParent),
            dimmed: !ancestorsEnabled(row.entry, byName: byName),
            isOpen: openGroups.contains(row.entry.name),
            onToggleGroup: { toggleGroupOpen(row.entry.name) }
          )
          if row.entry.isGroup, openGroups.contains(row.entry.name),
            (byParent[row.entry.name] ?? []).isEmpty
          {
            emptyHint("Empty group")
              .padding(.leading, kIndentW * CGFloat(min(row.depth + 1, kMaxIndentDepth)))
          }
        }

        ForEach(selections) { obj in
          ObjectRowView(entry: obj)
        }
      }
      .padding(.vertical, 2)
    }
    .background(ClassicUI.panelBg)
    .onAppear {
      refreshObjects()
    }
  }

  /// Expand/collapse a group row, mirroring the state into PyMOL so the two
  /// agree — `action=open|close` is what a `.pse` records, so a session saved
  /// from here reopens with the same groups expanded.
  private func toggleGroupOpen(_ name: String) {
    let opening = !openGroups.contains(name)
    if opening { openGroups.insert(name) } else { openGroups.remove(name) }
    let escaped = name.replacingOccurrences(of: "'", with: "")
    engine.runCommand("group \(escaped), , action=\(opening ? "open" : "close")", naming: escaped)
  }

  private func emptyHint(_ text: String) -> some View {
    Text(text)
      .font(.system(size: 11))
      .foregroundColor(PanelTheme.disabledColor)
      .frame(maxWidth: .infinity, alignment: .leading)
      .padding(.horizontal, 12)
      .padding(.vertical, 8)
  }

  // Pick-granularity menu (mouse_selection_mode): what a viewport tap selects.
  private func refreshObjects() {
    // Same file-based payload as the ~500ms poll (see parseObjectPanelFeedback):
    // poll_panel writes the list to a temp file and prints only the marker, so
    // a large object list can't overflow the feedback-line cap (#231).
    engine.objectPollEpoch = engine.objectListEpoch
    engine.runCommand(
      "python\n"
        + "from mmol import inspector as _ai\n"
        + "_ai.poll_panel()\n"
        + "python end"
    )
  }
}

// MARK: - Selection mode cycling

/// The "what does a click select" modes (atoms / residues / chains / …), in the
/// order Tab / Shift+Tab cycles them (#319, installKeyMonitor in
/// ContentView).
enum SelectionModeMenu {
  static let modes: [(Int, String)] = [
    (0, "Atoms"), (1, "Residues"), (2, "Chains"),
    (3, "Segments"), (4, "Objects"), (5, "Molecules"), (6, "C-α"),
  ]

  /// The mode one entry forward (Tab) or back (Shift+Tab) from `cur`, wrapping
  /// at both ends. The floored modulo also brings an out-of-range setting
  /// value back inside the list rather than going negative.
  static func nextMode(from cur: Int, forward: Bool) -> Int {
    let n = modes.count
    return ((cur + (forward ? 1 : -1)) % n + n) % n
  }

  /// Apply `nextMode` to the live setting — the Tab / Shift+Tab shortcut's action.
  static func cycle(_ engine: PyMOLEngine, forward: Bool) {
    let cur = Int(engine.sceneState.values["mouse_selection_mode"] ?? 1)
    let next = nextMode(from: cur, forward: forward)
    engine.runCommand("set mouse_selection_mode, \(next)")
    // Reflect the new mode immediately instead of waiting for the ~500ms
    // scene poll to report it back, so two quick Tabs advance two modes
    // rather than reading the same stale `cur` twice. The poll still has the
    // last word.
    engine.sceneState.values["mouse_selection_mode"] = Double(next)
  }
}

// MARK: - Object Row

/// A selection row.
private struct ObjectRowView: View {
  let entry: ObjectEntry
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    HStack(spacing: 0) { ObjectRowContent(entry: entry) }
      .objectRowChrome()
      .contentShape(Rectangle())
      .onTapGesture { engine.setObjectEnabled(entry.name, !entry.isEnabled) }
      // Right-click opens the action menu.
      .classicContextMenu { classicActionItems(for: entry, engine: engine) }
  }
}

// MARK: - Global "all" controls row

/// Pinned row above the object list giving the same A/S/H/L/C controls as an
/// object row but acting on the whole scene (selection "all") — mirrors desktop
/// PyMOL's "all" row for quick global Show/Hide/Label/Color and scene actions.
/// Show/Hide/Label/Color reuse the per-object menus with name "all"; the Action
/// (A) menu uses the global subset `allActionMenuItems`.
private struct AllControlsRow: View {
  @EnvironmentObject var engine: PyMOLEngine

  // The objects shown in the panel (selections excluded) — "all current objects".
  private var objects: [ObjectEntry] { engine.objects.filter { !$0.isSelection } }
  private var allEnabled: Bool { !objects.isEmpty && objects.allSatisfy { $0.isEnabled } }

  var body: some View {
    HStack(spacing: 0) {
      ClassicNameCell(text: "all", enabled: allEnabled, toggle: toggleAll)
      ClassicMenuTrigger {
        ASHLCFace(letter: "A")
      } items: {
        classicActionItems(allActionMenuItems, name: "all", engine: engine)
      }
      RepMenuButton(kind: .show, name: "all")
      RepMenuButton(kind: .hide, name: "all")
      RepMenuButton(kind: .label, name: "all")
      ColorMenuButton(name: "all")
    }
    .objectRowChrome()
    .contentShape(Rectangle())
    .onTapGesture { toggleAll() }
    .classicContextMenu { classicActionItems(allActionMenuItems, name: "all", engine: engine) }
  }

  // Enable or disable every object in one shot. Iterating get_names('objects')
  // targets exactly the current objects (not the "all" atom-selection), so it
  // works regardless of how enable/disable interpret the "all" keyword.
  private func toggleAll() {
    let enable = !allEnabled
    // Optimistic-first: flip every object's checkbox immediately (isEnabled
    // is a var → instant re-render) so the whole "all" row updates this frame
    // instead of waiting on the ~500ms pollObjects reconcile.
    for idx in engine.objects.indices where !engine.objects[idx].isSelection {
      engine.objects[idx].isEnabled = enable
    }
    engine.objectListEpoch += 1
    let action = enable ? "enable" : "disable"
    engine.runPython(
      "from pymol import cmd as _c\n"
        + "for _o in (_c.get_names('objects') or []):\n"
        + "    _c.\(action)(_o)\n")
    engine.refetchSequencesSoon()
  }
}

// MARK: - Action (A) Menu

/// Align one object onto a target molecule/selection over polymer CA atoms,
/// creating an alignment object (mirrors desktop PyMOL's align_to_object /
/// align_to_sele in modules/pymol/menu.py).
private func runAlignCommand(mobile: String, target: String, engine: PyMOLEngine) {
  let cmd =
    "python\ncmd.align(\"polymer and name CA and (\(mobile))\", "
    + "\"polymer and name CA and (\(target))\", quiet=0, "
    + "object=\"aln_\(mobile)_to_\(target)\", reset=1)\npython end"
  engine.runCommand(cmd, naming: mobile, target)
}

private struct ActionMenuButton: View {
  let entry: ObjectEntry
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    ClassicMenuTrigger {
      ASHLCFace(letter: "A")
    } items: {
      classicActionItems(for: entry, engine: engine)
    }
  }
}

// MARK: - Show (S) / Hide (H) / Label (L) Menu Buttons

private struct RepMenuButton: View {
  enum Kind: String {
    case show = "S"
    case hide = "H"
    case label = "L"
  }
  let kind: Kind
  let name: String
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    ClassicMenuTrigger {
      ASHLCFace(letter: kind.rawValue)
    } items: {
      switch kind {
      case .show: showItems
      case .hide: hideItems
      case .label: labelItems
      }
    }
  }
}

// MARK: - Color (C) Menu Button

private struct ColorMenuButton: View {
  let name: String
  @EnvironmentObject var engine: PyMOLEngine
  @State private var customColor: Color = .white

  @State private var showCustom = false

  var body: some View {
    ClassicMenuTrigger {
      ASHLCFace(letter: "C")
    } items: {
      classicItems
    }
    // Stable AX hook so UI tests can open a specific row's color menu
    // (the visible label "C" is shared by every row).
    .accessibilityIdentifier("colorMenu.\(name)")
    .popover(isPresented: $showCustom, arrowEdge: .bottom) {
      VStack(spacing: 8) {
        Text("Custom color").font(.system(size: 11, weight: .semibold))
        DebouncedColorPicker(
          get: { customColor },
          apply: {
            customColor = $0
            applyCustomColor($0)
          })
      }
      .padding(12)
    }
  }

  /// The colour menu: schemes, the tiered colour families, Custom….
  private var classicItems: [ClassicMenuItem] {
    var items: [ClassicMenuItem] = colorOptions.map { opt in
      if opt.label == "---" { return .separator }
      let command = opt.command ?? ""
      return .action(opt.label, swatch: opt.swatch) { applyColor(command: command) }
    }
    items.append(.separator)
    for row in PyMOLColorMenu.rows {
      let shades: [ClassicMenuItem] = row.family.shades.map { shade in
        .action(shade.name, swatch: shade.swatch) { applyColor(command: shade.name) }
      }
      let primary = row.primary.map { p in { applyColor(command: p.name) } }
      items.append(.submenu(row.title, swatch: row.swatch, primary: primary, shades))
    }
    items.append(.separator)
    items.append(.action("Custom…") { showCustom = true })
    return items
  }

  private func applyColor(command: String) {
    if command == "by_file" {
      engine.loadHeatmapCSV(for: name)
    } else {
      engine.runCommand(
        colorSchemeCommand(command, sele: name) ?? "color \(command), \(name)", naming: name)
    }
  }

  private func applyCustomColor(_ color: Color) {
    engine.runCommand(
      "set_color mmol_custom, \(rgb01List(color))\ncolor mmol_custom, \(name)", naming: name)
  }
}

// MARK: - Inspector: color helpers

/// The label column of an inspector row.
private func inspectorLabel(_ text: String) -> some View {
  Text(text)
    .font(.system(size: 10))
    .foregroundColor(PanelTheme.textColor)
    .frame(width: 78, alignment: .leading)
}

/// The command applying a C-menu coloring scheme (util.cnc, spectrum, …) to
/// `sele`, or nil when `scheme` is a plain color name.
private func colorSchemeCommand(_ scheme: String, sele: String) -> String? {
  switch scheme {
  case "spectrum": return "spectrum count, selection=\(sele)"
  // menu.py's rainbow: one colour per residue, N → C.
  case "rainbow": return "spectrum count, rainbow, \(sele), byres=1"
  // Color by B-factor: blue (low) → white → red (high), the classic
  // temperature look. Falls back gracefully if b is uniform/zero.
  case "spectrum_b": return "spectrum b, blue_white_red, \(sele)"
  case let s where s.hasPrefix("util."): return "python\ncmd.\(s)('\(sele)')\npython end"
  default: return nil
  }
}

private func colorFromHex(_ hex: String) -> Color? {
  var s = hex
  if s.hasPrefix("#") { s.removeFirst() }
  guard s.count == 6, let v = Int(s, radix: 16) else { return nil }
  return Color(
    .sRGB,
    red: Double((v >> 16) & 0xff) / 255.0,
    green: Double((v >> 8) & 0xff) / 255.0,
    blue: Double(v & 0xff) / 255.0)
}

/// A ColorPicker whose apply action is debounced. SwiftUI fires the binding's
/// `set:` on every tick of the color-wheel drag, which otherwise floods the
/// core with set_color + recolor each tick (review L-58 — a recolor is far
/// heavier than a `set`). We show the dragged color live but coalesce the
/// actual apply to ~150 ms after the last change, mirroring LabeledSlider.
struct DebouncedColorPicker: View {
  let get: () -> Color
  let apply: (Color) -> Void
  @State private var pending: Color? = nil
  @State private var work: DispatchWorkItem? = nil

  var body: some View {
    ClassicColorWell(
      selection: Binding(
        get: { pending ?? get() },
        set: { c in
          pending = c  // live: track the dragged color
          work?.cancel()
          let w = DispatchWorkItem { apply(c) }  // debounced: hit the core
          work = w
          DispatchQueue.main.asyncAfter(deadline: .now() + 0.15, execute: w)
        }))
  }
}

/// SwiftUI Color → sRGB components in 0…1, alpha included.
private func rgba01(_ color: Color) -> (r: Double, g: Double, b: Double, a: Double) {
  let ns = NSColor(color).usingColorSpace(.sRGB) ?? NSColor.white
  return (ns.redComponent, ns.greenComponent, ns.blueComponent, ns.alphaComponent)
}

/// SwiftUI Color → PyMOL set_color list "[r,g,b]" in 0…1. Alpha is dropped:
/// PyMOL colors have none.
private func rgb01List(_ color: Color) -> String {
  let c = rgba01(color)
  return String(format: "[%.3f,%.3f,%.3f]", c.r, c.g, c.b)
}

private let inspectorNamedColors: [(name: String, color: Color)] = [
  ("red", .red), ("green", Color(.sRGB, red: 0, green: 0.9, blue: 0)),
  ("blue", Color(.sRGB, red: 0.1, green: 0.3, blue: 1)),
  ("yellow", .yellow), ("orange", .orange),
  ("magenta", Color(.sRGB, red: 1, green: 0, blue: 1)), ("cyan", .cyan),
  ("grey70", Color(.sRGB, white: 0.7)), ("grey30", Color(.sRGB, white: 0.3)),
  ("white", .white), ("black", .black),
]

// MARK: - Inspector controls

/// Slider + editable numeric field, debounced live updates, exact commit.
private struct LabeledSlider: View {
  let prop: RepProperty
  let value: Double
  let onLive: (Double) -> Void
  let onCommit: (Double) -> Void

  @State private var local: Double = 0
  @State private var text: String = ""
  @State private var editing = false
  @State private var debounce: DispatchWorkItem?
  @State private var lastLiveAt: Date = .distantPast

  var body: some View {
    HStack(spacing: 6) {
      // Continuous (no `step`): a stepped Slider draws busy tick marks on
      // macOS. Values are rounded for display (fmt) and when applied
      // (fmtScene / the mm/zoom mappings), so dropping the step only removes
      // the ticks, not the effective granularity.
      ClassicSlider(
        value: $local, range: prop.min...prop.max,
        onEditingChanged: { began in
          editing = began
          if !began { onCommit(local) }
        })
      TextField("", text: $text)
        .textFieldStyle(.plain)
        .multilineTextAlignment(.trailing)
        .font(.system(size: 10, design: .monospaced))
        .frame(width: 40)
        .foregroundColor(PanelTheme.textColor)
        .onSubmit {
          if let v = Double(text) {
            local = Swift.min(Swift.max(v, prop.min), prop.max)
            onCommit(local)
          }
          text = fmt(local)
        }
    }
    .onAppear {
      local = value
      text = fmt(value)
    }
    .onChange(of: value) { v in
      if !editing {
        local = v
        text = fmt(v)
      }
    }
    .onChange(of: local) { v in
      text = fmt(v)
      if editing { scheduleLive(v) }
    }
  }

  private func fmt(_ v: Double) -> String { String(format: "%.\(prop.decimals)f", v) }
  // Leading-edge throttle (~30 Hz): fire onLive immediately during a
  // continuous drag, rate-limited, with a trailing call so the final value
  // always lands. A plain trailing debounce (cancel + reschedule on every
  // change) starved the update until the finger paused, so the view chased
  // the drag instead of tracking it — the perceived slider lag.
  private func scheduleLive(_ v: Double) {
    let interval = 0.033
    debounce?.cancel()
    let now = Date()
    let wait = interval - now.timeIntervalSince(lastLiveAt)
    if wait <= 0 {
      lastLiveAt = now
      onLive(v)
    } else {
      let w = DispatchWorkItem {
        lastLiveAt = Date()
        onLive(v)
      }
      debounce = w
      DispatchQueue.main.asyncAfter(deadline: .now() + wait, execute: w)
    }
  }
}

private struct SegmentedSetting: View {
  let prop: RepProperty
  let value: Double
  let onSelect: (Double) -> Void

  var body: some View {
    HStack(spacing: 2) {
      ForEach(Array(prop.options.enumerated()), id: \.offset) { _, opt in
        let sel = abs(opt.value - value) < 0.5
        Button(action: { onSelect(opt.value) }) { Text(opt.label).frame(minWidth: 12) }
          .buttonStyle(ClassicButtonStyle(selected: sel, compact: true))
      }
    }
  }
}

private struct ToggleSetting: View {
  let value: Double
  let onToggle: (Bool) -> Void
  // Optimistic local state. Previously the switch derived its position directly
  // from `value` (a ~500ms-lagged poll) on EVERY re-render, so it flickered /
  // snapped back: opening the panel, a rotation-driven refresh, or the gap
  // between a tap and the next poll would re-render with the stale value and
  // flip the switch (and its tint) back. Driving the switch from local @State
  // that changes ONLY on a user flip or a genuine polled-value change makes a
  // re-render with an unchanged value a no-op for the switch.
  @State private var isOn: Bool

  // Seeded from the value (not set in onAppear), so an opening card shows each
  // switch in place instead of animating it from off to on.
  init(value: Double, onToggle: @escaping (Bool) -> Void) {
    self.value = value
    self.onToggle = onToggle
    _isOn = State(initialValue: value > 0.5)
  }

  var body: some View {
    Toggle("", isOn: $isOn)
      .labelsHidden()
      .toggleStyle(ClassicSwitchStyle())
      .onChange(of: value) { v in
        let want = v > 0.5
        if want != isOn { isOn = want }  // adopt real external changes
      }
      .onChange(of: isOn) { on in
        if on != (value > 0.5) { onToggle(on) }  // user flip → push once
      }
  }
}

/// Per-rep color OVERRIDE control (writes `<rep>_color`): Inherit / named / custom RGB.
private struct RepColorControl: View {
  let objName: String
  let rep: String
  let colorSetting: String
  let defaultColor: Int
  let colorState: String  // "inherit" or "#rrggbb"
  var effColor: String? = nil  // the color on screen
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    HStack(spacing: 6) {
      // Coloring schemes are atom-level (the standard PyMOL behavior); this
      // rep is reset to inherit so the scheme shows through on it.
      ClassicMenuTrigger {
        ClassicDropdownFace(text: colorState == "inherit" ? "Inherit" : "Custom", compact: true)
      } items: {
        [
          .action("Inherit", destructive: false) { setOverride("\(defaultColor)") }, .separator,
          .action("by element") { applyScheme("util.cnc") },
          .action("by chain") { applyScheme("util.cbc") },
          .action("by ss") { applyScheme("util.cbss") },
          .action("spectrum") { applyScheme("spectrum") },
          .action("by b-factor") { applyScheme("spectrum_b") }, .separator,
        ]
          + inspectorNamedColors.map { c in .action(c.name, swatch: c.color) { setOverride(c.name) }
          }
      }
      DebouncedColorPicker(
        get: { colorFromHex(effColor ?? colorState) ?? .white },
        apply: { applyCustom($0) })
    }
  }

  private func setOverride(_ c: String) {
    engine.runCommand("set \(colorSetting), \(c), \(objName)", naming: objName)
  }
  // Apply an atom-level coloring scheme, resetting this rep to inherit so the
  // scheme is visible on it (PyMOL has no true per-rep scheme coloring).
  private func applyScheme(_ s: String) {
    setOverride("\(defaultColor)")
    if let c = colorSchemeCommand(s, sele: objName) { engine.runCommand(c, naming: objName) }
  }
  private func applyCustom(_ color: Color) {
    let nm = "tmp_\(sanitizedName(objName))_\(rep)"
    engine.runCommand(
      "set_color \(nm), \(rgb01List(color))\nset \(colorSetting), \(nm), \(objName)",
      naming: objName)
  }
}

/// Color control bound to an arbitrary per-rep color setting (e.g.
/// surface_contour_color). "Inherit" sets -1 (here the surface color); named
/// colors / the custom picker set an explicit color. Mirrors RepColorControl but
/// without the atom-coloring schemes (which don't apply to a contour line).
private struct SettingColorControl: View {
  let objName: String
  let rep: String
  let setting: String
  let colorState: String  // "inherit" or "#rrggbb"
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    HStack(spacing: 6) {
      ClassicMenuTrigger {
        ClassicDropdownFace(text: colorState == "inherit" ? "Inherit" : "Custom", compact: true)
      } items: {
        [.action("Inherit", destructive: false) { setColor("-1") }, .separator]
          + inspectorNamedColors.map { c in .action(c.name, swatch: c.color) { setColor(c.name) } }
      }
      DebouncedColorPicker(
        get: { colorFromHex(colorState) ?? .white },
        apply: { applyCustom($0) })
    }
  }

  private func setColor(_ c: String) {
    engine.runCommand("set \(setting), \(c), \(objName)", naming: objName)
  }
  private func applyCustom(_ color: Color) {
    let nm = "tmp_\(sanitizedName(objName))_\(rep)_ctr"
    engine.runCommand(
      "set_color \(nm), \(rgb01List(color))\nset \(setting), \(nm), \(objName)", naming: objName)
  }
}

// MARK: - Object row content (shared by selection rows and object cards)

private struct ObjectRowContent: View {
  let entry: ObjectEntry
  // Tri-state for group rows (some descendants on). Defaults keep the plain
  // object/selection behaviour for callers that don't compute it.
  var check: CheckState? = nil
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    // Desktop PyMOL brackets selection names, without an atom count:
    // "(sele)", not "(sele (180))".
    ClassicNameCell(
      text: entry.isSelection ? "(\(entry.name))" : entry.displayName,
      enabled: (check ?? (entry.isEnabled ? .on : .off)) != .off,
      bold: entry.isGroup, toggle: toggleEnabled)

    // Action buttons: A S H L C
    ActionMenuButton(entry: entry)
    RepMenuButton(kind: .show, name: entry.name)
    RepMenuButton(kind: .hide, name: entry.name)
    RepMenuButton(kind: .label, name: entry.name)
    ColorMenuButton(name: entry.name)
  }

  private func toggleEnabled() {
    engine.setObjectEnabled(entry.name, !entry.isEnabled)
  }
}

// MARK: - Expandable object card

private struct ObjectCard: View {
  let entry: ObjectEntry
  // Group-tree presentation (#255). All defaulted so the preview/memberwise
  // call sites that predate the tree keep compiling.
  var depth: Int = 0
  var check: CheckState = .on
  var dimmed: Bool = false
  var isOpen: Bool = false
  var onToggleGroup: (() -> Void)? = nil
  @EnvironmentObject var engine: PyMOLEngine
  @State private var selectedRep: String?
  // Live state-slider value while dragging (nil = follow the object poll). Keeps
  // the thumb from snapping back to the last-polled state mid-drag.
  @State private var scrubState: Int? = nil

  private var expanded: Bool { engine.expandedDetail == entry.name }
  private var reps: [RepState] { engine.objectDetails[entry.name] ?? [] }
  // Reps currently SHOWN (have drawn atoms, from the poll).
  private var activeSet: Set<String> { Set(reps.map { $0.rep }) }
  // Reps hidden via the Visible toggle but kept listed as layers.
  private var keptHidden: Set<String> { engine.keptHidden[entry.name] ?? [] }
  // Layers shown in the inspector = shown ∪ kept-hidden, in catalog order.
  /// Every rep the object draws (e.g. cartoon for the protein and sticks for a
  /// ligand), plus layers hidden with the Visible switch, in catalog order.
  private var listedReps: [String] {
    RepCatalog.order.filter(activeSet.union(keptHidden).contains)
  }
  private var currentRep: String? {
    if let s = selectedRep, listedReps.contains(s) { return s }
    return listedReps.first
  }
  private var heatmapSelected: Bool {
    engine.heatmapChipActive.contains(entry.name) && engine.heatmaps[entry.name] != nil
  }

  var body: some View {
    VStack(spacing: 0) {
      HStack(spacing: 0) {
        // Indent nested members. Applied in THIS HStack, before the
        // chevron — ObjectRowContent is a bare column tuple spliced into
        // this stack, so wrapping it to indent would break its alignment
        // with the selection and "all" rows.
        if depth > 0 {
          Spacer().frame(width: kIndentW * CGFloat(min(depth, kMaxIndentDepth)))
        }
        // Groups get +/-, plain objects keep the chevron, because the two
        // disclosures do genuinely DIFFERENT things: +/- expands the
        // group's children, the chevron opens that object's rep inspector.
        // Sharing one glyph read as a single behaviour that sometimes
        // produced a subtree and sometimes a settings card. The group's own
        // rep card is reached from the A menu ("Inspect Group…").
        // CIRCLE variants, not square and not bare: the tri-state checkbox
        // immediately to the right is minus.square.fill, and an expanded
        // group would otherwise put a bare "-" right beside a boxed "[-]".
        // The enclosing shape is what separates them at a glance.
        // A raised square in the hidden-object grey, the shown-object grey
        // while open — a bigger target than the bare chevron was.
        let open = entry.isGroup ? isOpen : expanded
        Button(action: { entry.isGroup ? (onToggleGroup?() ?? ()) : toggleExpand() }) {
          Image(
            systemName: entry.isGroup
              ? (isOpen ? "minus" : "plus")
              : (expanded ? "chevron.down" : "chevron.right")
          )
          .font(.system(size: 8, weight: .bold))
          .foregroundColor(Color(white: 0.95))
          .frame(width: 16, height: 16)
          .background(ClassicBevel(fill: open ? ClassicUI.nameEnabled : ClassicUI.nameDisabled))
          .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .padding(.trailing, 2)
        .accessibilityLabel(
          entry.isGroup
            ? "\(isOpen ? "Collapse" : "Expand") group \(entry.name)"
            : "\(expanded ? "Collapse" : "Expand") \(entry.name)")
        ObjectRowContent(entry: entry, check: check)
      }
      .objectRowChrome()
      // A member hidden by a disabled ancestor is dimmed: it still reports
      // enabled, so without this the box reads as checked while nothing is
      // on screen. Opacity rather than a new colour, because the A/S/H/L/C
      // cells compute their own.
      .opacity(dimmed ? 0.45 : 1.0)
      // Whole HEADER-row tap toggles enable. Applied to the header HStack
      // only (NOT the outer VStack) so the expanded rep-detail below stays
      // interactive. The disclosure chevron Button and the trailing
      // A/S/H/L/C controls consume their own taps (expand-only for the
      // chevron); only the dead Spacer gap falls through to here.
      .contentShape(Rectangle())
      .onTapGesture {
        engine.setObjectEnabled(entry.name, !entry.isEnabled)
      }
      // Right-click opens the action menu.
      .classicContextMenu { classicActionItems(for: entry, engine: engine) }

      // Drawn once its first poll is in, so it never shows an empty card first.
      if expanded, engine.objectDetails[entry.name] != nil {
        VStack(spacing: 3) {
          // Multi-state objects (NMR models / trajectory frames) get a
          // STATE row: pin this object to a state independent of the
          // global timeline, overlay all states, or fit/split them.
          if entry.stateCount > 1 {
            stateRow()
            Divider().background(PanelTheme.disabledColor.opacity(0.3))
          }
          // Object/layer-level coloring (by element/chain/ss/spectrum/
          // named) is the structure row's "C" button — not duplicated
          // here. The per-rep grid below controls per-rep color overrides.
          // Always show the chips bar (it holds the "+" add menu) so a
          // layer can be added even after the last one is deleted.
          RepChips(
            objName: entry.name, listed: listedReps,
            active: activeSet,
            current: heatmapSelected ? kHeatmapChip : currentRep,
            heatmap: engine.heatmaps[entry.name] != nil,
            onSelect: { chip in
              if chip == kHeatmapChip {
                engine.heatmapChipActive.insert(entry.name)
              } else {
                engine.heatmapChipActive.remove(entry.name)
                selectedRep = chip
              }
            })
          if heatmapSelected {
            HeatmapPanel(objName: entry.name)
          } else if let rep = currentRep {
            // Always present (even when hidden): show/hide the layer
            // + delete it. Hiding keeps the layer listed so it can be
            // toggled back on; the X removes it.
            layerHeader(rep)
            // Full per-rep settings only while the layer is shown
            // (a hidden layer reports no state).
            if let spec = RepCatalog.spec(rep),
              let st = reps.first(where: { $0.rep == rep })
            {
              RepPropertyGrid(objName: entry.name, spec: spec, state: st)
            }
          } else {
            Text("No representations shown — tap + to add one.")
              .font(.system(size: 10))
              .foregroundColor(PanelTheme.disabledColor)
              .frame(maxWidth: .infinity, alignment: .leading)
          }
        }
        .padding(.leading, 18)
        .padding(.trailing, 8)
        .padding(.vertical, 4)
        .background(PanelTheme.rowAltBackground.opacity(0.6))
      }
    }
  }

  private func toggleExpand() {
    // Accordion: opening this card closes whatever else was open.
    engine.expandedDetail = expanded ? nil : entry.name
  }

  // Per-object STATE controls for multi-state objects (NMR / trajectory).
  // The slider/steppers PIN this object to a state via `set state, N, obj`
  // (so it stops following the global timeline); "Sync" un-pins it. Distinct
  // from the global "Overlay all states" in the SCENE card.
  @ViewBuilder
  private func stateRow() -> some View {
    let total = max(entry.stateCount, 1)
    let meta = engine.objectMeta[entry.name]
    // Use the object's effective state from the poll; default to 1 (avoid
    // depending on playback.currentFrame so the inspector doesn't re-render
    // on every frame tick during playback).
    // While dragging, show the live scrub value; otherwise follow the poll.
    let cur = min(max(scrubState ?? meta?.state ?? 1, 1), total)
    VStack(spacing: 3) {
      HStack(spacing: 6) {
        inspectorLabel("State")
        Button {
          scrubState = nil
          setState(max(cur - 1, 1))
        } label: {
          Image(systemName: "minus").font(.system(size: 8, weight: .bold))
        }
        .buttonStyle(ClassicButtonStyle(compact: true))
        ClassicSlider(
          value: Binding(
            get: { Double(cur) },
            set: { v in
              let n = min(max(Int(v.rounded()), 1), total)
              if n != cur {
                scrubState = n
                setState(n)
              }  // apply only on integer change
            }),
          range: 1...Double(max(total, 2)), step: 1,
          onEditingChanged: { editing in if !editing { scrubState = nil } })
        Button {
          scrubState = nil
          setState(min(cur + 1, total))
        } label: {
          Image(systemName: "plus").font(.system(size: 8, weight: .bold))
        }
        .buttonStyle(ClassicButtonStyle(compact: true))
        Text("\(cur)/\(total)")
          .font(.system(size: 10, design: .monospaced))
          .foregroundColor(PanelTheme.textColor)
          .frame(width: 42, alignment: .trailing)
      }
      // Compound name for this state, when present (e.g. a multi-record
      // SDF's per-record title), so docking-pose / library sets can be
      // navigated by name while scrubbing states (issue #203).
      if let name = meta?.title(forState: cur) {
        HStack(spacing: 6) {
          inspectorLabel("Name")
          Text(name)
            .font(.system(size: 10, weight: .medium))
            .foregroundColor(TimelineTheme.accent)
            .lineLimit(1).truncationMode(.middle)
            .textSelection(.enabled)
            .help(name)
          Spacer(minLength: 0)
        }
      }
      // Play/pause + per-object fps — animates THIS object's models (via a
      // Swift timer + `set state`), independent of the movie and other objects.
      HStack(spacing: 6) {
        inspectorLabel("Play")
        Button {
          engine.toggleObjectStates(entry.name)
        } label: {
          Image(systemName: engine.playingObjects.contains(entry.name) ? "pause.fill" : "play.fill")
            .font(.system(size: 13, weight: .semibold))
            .foregroundColor(
              engine.movieTabActive ? PanelTheme.disabledColor : TimelineTheme.accent
            )
            .frame(width: 24, height: 20).contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .disabled(engine.movieTabActive)  // authoring a movie → model playback off
        // A disabled Button swallows its own hover events, so `.help()` on
        // it never fires. Wrap it so the tooltip lives on a non-disabled
        // container that still receives hover while playback is off.
        .allowsHitTesting(!engine.movieTabActive)
        .overlay {
          if engine.movieTabActive {
            Color.clear.contentShape(Rectangle())
              .help(
                "Model playback is off while the movie timeline is open — close it to inspect models"
              )
          }
        }
        .help(engine.playingObjects.contains(entry.name) ? "Pause" : "Play models")
        Spacer(minLength: 0)
        ClassicPicker(
          compact: true,
          selection: Binding(
            get: { Int(engine.objectPlaybackFPS(entry.name)) },
            set: { engine.setObjectFPS(entry.name, Double($0)) }),
          options: [1, 5, 10, 15, 30].map { ($0, "\($0) fps") }
        )
        .help("Playback speed for this object")
      }
      HStack(spacing: 6) {
        inspectorLabel("Overlay all")
        ToggleSetting(value: (meta?.overlayAll ?? false) ? 1 : 0) { on in
          engine.runCommand("set all_states, \(on ? 1 : 0), \(entry.name)", naming: entry.name)
        }
        Spacer(minLength: 0)
        stateActionButton("Fit") {
          engine.runCommand("intra_fit \(entry.name)", naming: entry.name)
        }
        stateActionButton("Split") {
          engine.runCommand("split_states \(entry.name)", naming: entry.name)
        }
        stateActionButton("Sync") {
          engine.runCommand("unset state, \(entry.name)", naming: entry.name)
        }
      }
    }
  }

  /// Drop the object-level overrides of the rep's card settings and color.
  private func resetRep(_ rep: String) {
    guard let spec = RepCatalog.spec(rep) else { return }
    let settings =
      spec.properties.map(\.setting) + (spec.colorSetting.isEmpty ? [] : [spec.colorSetting])
    engine.runCommand(
      Set(settings).sorted().map { "unset \($0), \(entry.name)" }.joined(separator: "; "),
      naming: entry.name)
  }

  private func setState(_ n: Int) {
    engine.runCommand("set state, \(n), \(entry.name)", naming: entry.name)
  }

  private func stateActionButton(_ title: String, _ action: @escaping () -> Void) -> some View {
    Button(title, action: action).buttonStyle(ClassicButtonStyle(compact: true))
  }

  // Visible toggle + delete (X) for the current layer. Shown whether the layer
  // is visible or hidden — hiding keeps it listed (toggle back on to reset),
  // the X removes the layer entirely.
  @ViewBuilder
  private func layerHeader(_ rep: String) -> some View {
    let shown = activeSet.contains(rep)
    HStack(spacing: 6) {
      inspectorLabel("Visible")
      ToggleSetting(value: shown ? 1 : 0) { on in
        if on {
          // Stays in keptHidden: until the next poll reports it shown, that is
          // what keeps this tab listed (and selected) instead of jumping away.
          engine.runCommand("show \(rep), \(entry.name)", naming: entry.name)
        } else {
          engine.runCommand("hide \(rep), \(entry.name)", naming: entry.name)
          engine.keptHidden[entry.name, default: []].insert(rep)
        }
      }
      Spacer(minLength: 0)
      Button("Reset") { resetRep(rep) }
        .buttonStyle(ClassicButtonStyle(compact: true))
        .help("Return this representation's settings and color to the defaults")
      Button {
        engine.runCommand("hide \(rep), \(entry.name)", naming: entry.name)
        engine.keptHidden[entry.name]?.remove(rep)
        if selectedRep == rep { selectedRep = nil }
      } label: {
        Image(systemName: "xmark").font(.system(size: 8, weight: .bold))
      }
      .buttonStyle(ClassicButtonStyle(compact: true))
      .accessibilityLabel("Delete \(RepCatalog.display(rep)) layer")
    }
  }
}

// MARK: - Representation chips

private struct RepChips: View {
  let objName: String
  let listed: [String]  // shown ∪ kept-hidden, in catalog order
  let active: Set<String>  // currently shown (others are dimmed = hidden)
  let current: String?
  var heatmap: Bool = false  // the object has a CSV heatmap (C ▸ by file…)
  let onSelect: (String) -> Void
  @EnvironmentObject var engine: PyMOLEngine

  // The "+" menu offers only reps not already listed as a layer.
  private var inactive: [String] {
    RepCatalog.order.filter { !listed.contains($0) }
  }

  var body: some View {
    ScrollView(.horizontal, showsIndicators: false) {
      HStack(spacing: 4) {
        // First, so it is never scrolled out of view.
        if heatmap {
          let sel = current == kHeatmapChip
          Button("heatmap") { onSelect(kHeatmapChip) }
            .buttonStyle(ClassicButtonStyle(selected: sel, compact: true))
        }
        ForEach(listed, id: \.self) { rep in
          let sel = rep == current
          let shown = active.contains(rep)
          Button(RepCatalog.display(rep)) { onSelect(rep) }
            .buttonStyle(ClassicButtonStyle(selected: sel, compact: true))
            .opacity(shown ? 1.0 : 0.45)  // hidden layers are dimmed
        }
        ClassicMenuTrigger {
          Image(systemName: "plus").font(.system(size: 8, weight: .bold))
            .foregroundColor(Color(white: 0.95))
            .frame(width: 20, height: 17)
            .background(ClassicBevel(fill: Color(white: 0.4)))
        } items: {
          inactive.map { rep in
            .action(RepCatalog.display(rep), destructive: false) {
              engine.runCommand("show \(rep), \(objName)", naming: objName)
              engine.keptHidden[objName]?.remove(rep)
              onSelect(rep)
            }
          }
        }
      }
    }
  }
}

// MARK: - Property grid

private struct RepPropertyGrid: View {
  let objName: String
  let spec: RepSpec
  let state: RepState
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    VStack(spacing: 3) {
      if !spec.colorSetting.isEmpty {
        gridRow("Color") {
          RepColorControl(
            objName: objName, rep: spec.rep,
            colorSetting: spec.colorSetting,
            defaultColor: spec.defaultColor,
            colorState: state.color, effColor: state.effColor)
        }
      }
      ForEach(spec.properties) { p in
        gridRow(p.label) { control(for: p) }
        // Per-atom transparency detail sits directly under the matching
        // transparency slider so it's clear the slider is only a baseline.
        if let at = state.atomTransp, at.setting == p.setting {
          atomTranspRow(at)
        }
      }
      // Fallback: surface the detail even if this rep's spec has no slider
      // for its transparency setting (so the info is never lost).
      if let at = state.atomTransp,
        !spec.properties.contains(where: { $0.setting == at.setting })
      {
        atomTranspRow(at)
      }
    }
    .padding(.top, 2)
  }

  // "per-atom: min–max" readout + a Clear action, shown when atom-level
  // transparency overrides the object-level slider for this rep.
  @ViewBuilder
  private func atomTranspRow(_ at: AtomTransp) -> some View {
    HStack(spacing: 6) {
      Image(systemName: "drop.halffull").font(.system(size: 10))
      Text("per-atom: \(rangeLabel(at))").font(.system(size: 10))
      Spacer(minLength: 4)
      Button("Clear") { clearAtomTransp(at.setting) }
        .buttonStyle(ClassicButtonStyle(compact: true))
    }
    .foregroundColor(PanelTheme.atomTranspColor)
    .padding(.leading, 84)  // align under the control column (label width 78 + gap)
    .help(
      "Some atoms set \(at.setting) individually, so the slider above only sets a baseline. Clear removes the per-atom values and hands control back to the slider."
    )
  }

  private func rangeLabel(_ at: AtomTransp) -> String {
    let lo = fmtTransp(at.min)
    let hi = fmtTransp(at.max)
    return lo == hi ? lo : "\(lo)–\(hi)"
  }

  private func fmtTransp(_ v: Double) -> String {
    var s = String(format: "%.2f", v)
    if s.contains(".") {
      while s.hasSuffix("0") { s.removeLast() }
      if s.hasSuffix(".") { s.removeLast() }
    }
    return s
  }

  // Remove the per-atom overrides for `setting` on this object. The
  // atom-selection form `unset(setting, (obj))` clears atom-level values while
  // keeping the object-level slider value; a rebuild refreshes the baked
  // cartoon/surface geometry, then re-poll so the row disappears promptly.
  private func clearAtomTransp(_ setting: String) {
    engine.runCommand("unset \(setting), (\(objName))", naming: objName)
    engine.runCommand("rebuild \(objName)", naming: objName)
    engine.refreshExpandedDetail()
  }

  @ViewBuilder
  private func control(for p: RepProperty) -> some View {
    let v = state.values[p.setting] ?? 0
    switch p.kind {
    case .slider:
      // commitOnly props skip live updates — they apply once on release
      // (avoids re-running an expensive rebuild on every drag tick).
      LabeledSlider(
        prop: p, value: v,
        onLive: { if !p.commitOnly { engine.setSilently(p.setting, fmt($0), objName) } },
        onCommit: { set(p.setting, $0) })
    case .segmented:
      SegmentedSetting(prop: p, value: v) { set(p.setting, $0) }
    case .toggle:
      ToggleSetting(value: v) { set(p.setting, $0 ? 1 : 0) }
    case .color:
      // A standalone color control bound to an arbitrary color setting
      // (e.g. surface_contour_color). -1 = inherit (here: the surface color).
      SettingColorControl(
        objName: objName, rep: spec.rep, setting: p.setting,
        colorState: state.settingColors[p.setting] ?? "inherit")
    }
  }

  private func set(_ setting: String, _ value: Double) {
    engine.runCommand("set \(setting), \(fmt(value)), \(objName)", naming: objName)
  }

  private func fmt(_ value: Double) -> String {
    (value == value.rounded()) ? String(Int(value)) : String(format: "%.4f", value)
  }

  @ViewBuilder
  private func gridRow<Content: View>(_ label: String, @ViewBuilder _ content: () -> Content)
    -> some View
  {
    HStack(spacing: 6) {
      inspectorLabel(label)
      content()
      Spacer(minLength: 0)
    }
  }
}

// MARK: - Scene (global) card

// Small (?) affordance that reveals a one-line description on click (and on hover).
private struct HelpButton: View {
  let text: String
  @State private var show = false
  var body: some View {
    Button {
      show.toggle()
    } label: {
      Image(systemName: "questionmark.circle")
        .font(.system(size: 11))
        .foregroundColor(PanelTheme.disabledColor)
        .contentShape(Rectangle())
    }
    .buttonStyle(.plain)
    .help(text)
    .popover(isPresented: $show, arrowEdge: .trailing) {
      Text(text)
        .font(.system(size: 12))
        .foregroundColor(PanelTheme.textColor)
        .frame(width: 230, alignment: .leading)
        .fixedSize(horizontal: false, vertical: true)
        .padding(12)
    }
  }
}

// The body of the Viewer Settings tab (ViewerSettingsPanel).
struct SceneCard: View {
  @EnvironmentObject var engine: PyMOLEngine
  @State private var showSettings = false
  // Per-sub-group expand state; the heavier groups start collapsed.
  @State private var openGroups: Set<String> = ["Canvas", "Camera", "Lighting", "Effects"]
  var body: some View {
    VStack(spacing: 3) {
      ForEach(SceneCatalog.groups, id: \.self) { group in
        sceneGroup(group)
      }
      Button {
        showSettings = true
      } label: {
        HStack(spacing: 6) {
          Image(systemName: "slider.horizontal.3")
          Text("All Settings…")
          Spacer()
        }
        .font(.system(size: 12, weight: .medium))
      }
      .buttonStyle(ClassicButtonStyle(prominent: true))
      .padding(.top, 4).padding(.bottom, 6)
    }
    .padding(.horizontal, 8).padding(.vertical, 4)
    .background(PanelTheme.rowAltBackground.opacity(0.6))
    .sheet(isPresented: $showSettings) {
      SettingsSheet().background(ClassicUI.panelBg).classicControls()
    }
  }

  // A collapsible sub-group of Scene settings (Canvas, Camera, Lighting, …).
  @ViewBuilder
  private func sceneGroup(_ group: String) -> some View {
    let isOpen = openGroups.contains(group)
    VStack(spacing: 3) {
      Button {
        if isOpen { openGroups.remove(group) } else { openGroups.insert(group) }
      } label: {
        HStack(spacing: 4) {
          Image(systemName: isOpen ? "chevron.down" : "chevron.right")
            .font(.system(size: 8, weight: .semibold))
            .foregroundColor(PanelTheme.disabledColor)
            .frame(width: 10)
          Text(group.uppercased())
            .font(.system(size: 9.5, weight: .semibold))
            .tracking(0.4)
            .foregroundColor(PanelTheme.headerColor)
          Spacer()
        }
        .padding(.vertical, 2)
        .contentShape(Rectangle())
      }
      .buttonStyle(.plain)
      if isOpen {
        ForEach(SceneCatalog.params.filter { $0.group == group }) { p in
          SceneParamRow(param: p, engine: engine)
        }
        if group == "Effects" { resetEffectsButton }
      }
    }
  }

  // One-tap restore of the Effects group to its built-in defaults — an
  // obvious way back to a neutral render.
  private var resetEffectsButton: some View {
    Button {
      resetEffects()
    } label: {
      HStack(spacing: 6) {
        Image(systemName: "arrow.counterclockwise")
        Text("Reset effects to defaults")
        Spacer()
      }
      .font(.system(size: 12, weight: .medium))
    }
    .buttonStyle(ClassicButtonStyle())
    .padding(.leading, 10).padding(.vertical, 3)
  }

  // Defaults live on the engine (engine.resetEffects). The ~500ms
  // scene-state poll re-syncs the toggles/sliders after the sets land.
  private func resetEffects() { engine.resetEffects() }

}

// MARK: - Scene-setting row

// One scene-setting row (label + control + help) of SceneCard. Self-hides when
// its dependsOn parent toggle is off, so callers render it unconditionally.
struct SceneParamRow: View {
  let param: SceneParam
  @ObservedObject var engine: PyMOLEngine

  var body: some View {
    if let dep = param.dependsOn, (engine.sceneState.values[dep] ?? 0) <= 0.5 {
      EmptyView()
    } else {
      sceneRow(param.label, help: param.help) {
        sceneControl(param)
      }
      .padding(.leading, param.dependsOn != nil ? 12 : 0)
    }
  }

  @ViewBuilder
  private func sceneControl(_ p: SceneParam) -> some View {
    if p.isColor {
      DebouncedColorPicker(
        get: {
          let c = engine.sceneState.outlineColor
          return Color(
            .sRGB, red: c.count > 0 ? c[0] : 0,
            green: c.count > 1 ? c[1] : 0, blue: c.count > 2 ? c[2] : 0)
        },
        apply: { setOutlineColor($0) }
      )
      .frame(width: 28)
    } else {
      let v = engine.sceneState.values[p.setting] ?? 0
      switch p.kind {
      case .toggle:
        if p.setting == "hover_indicator" {
          // Not a PyMOL setting: an app-level preference (persisted via
          // UserDefaults) that gates the hover pre-selection preview.
          ToggleSetting(value: engine.hoverPreviewEnabled ? 1 : 0) { on in
            engine.hoverPreviewEnabled = on
            if !on { engine.clearHoverPreview() }
          }
        } else {
          ToggleSetting(value: v) { on in engine.runCommand("set \(p.setting), \(on ? 1 : 0)") }
        }
      case .segmented:
        SegmentedSetting(
          prop: RepProperty(
            setting: p.setting, label: p.label, kind: .segmented, options: p.options),
          value: v
        ) { engine.runCommand("set \(p.setting), \(Int($0))") }
      case .slider:
        if p.setting == "field_of_view" {
          // Lens control: the slider is a 35mm-equivalent focal length (mm),
          // shown by converting the polled field_of_view, and applied via
          // set_fov (a dolly zoom) so it swaps perspective rather than
          // zooming. Inert in orthoscopic mode.
          let fovDeg = engine.sceneState.values["field_of_view"] ?? 20
          let orthoOn = (engine.sceneState.values["ortho"] ?? 0) > 0.5
          LabeledSlider(
            prop: RepProperty(
              setting: p.setting, label: p.label, kind: .slider,
              min: p.min, max: p.max, step: p.step, decimals: p.decimals),
            value: fovToMM(fovDeg),
            onLive: { engine.runPython("from pymol import cmd as _c\n_c.set_fov(\(mmToFOV($0)))") },
            onCommit: { engine.runCommand("set_fov \(String(format: "%.3f", mmToFOV($0)))") }
          )
          .disabled(orthoOn)
          .opacity(orthoOn ? 0.4 : 1.0)
        } else if p.setting == "zoom" {
          // Zoom: absolute apparent MAGNIFICATION, independent of the Lens.
          // M = scene_radius / (cam_dist * tan(fov/2)) is invariant under
          // the Lens dolly-zoom, so the two sliders don't fight. Dragging
          // dollies the camera to the target distance.
          let camDist = engine.sceneState.values["cam_dist"] ?? 0
          let radius = engine.sceneState.values["scene_radius"] ?? 0
          let fovDeg = engine.sceneState.values["field_of_view"] ?? 20
          let zoomReady = radius > 0 && camDist > 0
          LabeledSlider(
            prop: RepProperty(
              setting: p.setting, label: p.label, kind: .slider,
              min: p.min, max: p.max, step: p.step, decimals: p.decimals),
            value: zoomMag(camDist: camDist, radius: radius, fovDeg: fovDeg),
            onLive: { engine.setZoomMagnification($0, radius: radius, fovDeg: fovDeg) },
            onCommit: { engine.setZoomMagnification($0, radius: radius, fovDeg: fovDeg) }
          )
          .disabled(!zoomReady)
          .opacity(zoomReady ? 1.0 : 0.4)
        } else {
          LabeledSlider(
            prop: RepProperty(
              setting: p.setting, label: p.label, kind: .slider,
              min: p.min, max: p.max, step: p.step, decimals: p.decimals),
            value: v,
            onLive: { engine.setSilently(p.setting, fmtScene($0, p)) },
            onCommit: { engine.runCommand("set \(p.setting), \(fmtScene($0, p))") })
        }
      case .color:
        EmptyView()  // scene colors use p.isColor above, not the .color kind
      }
    }
  }

  private func fmtScene(_ v: Double, _ p: SceneParam) -> String {
    p.decimals == 0 ? String(Int(v.rounded())) : String(format: "%.4f", v)
  }

  // Lens control: map PyMOL's vertical field_of_view (degrees) to a 35mm-
  // equivalent focal length using the full-frame sensor height (24mm):
  //   fov = 2*atan(12/f)   <->   f = 12/tan(fov/2)
  // fovToMM clamps into the slider's mm range so the thumb never leaves bounds.
  private func fovToMM(_ fovDeg: Double) -> Double {
    let f = 12.0 / tan(fovDeg * .pi / 360.0)
    return min(max(f, 12.0), 135.0)
  }
  private func mmToFOV(_ mm: Double) -> Double {
    2.0 * atan(12.0 / mm) * 180.0 / .pi
  }

  // Apparent magnification for the Zoom slider: ~1 when the scene fits the view,
  // higher when zoomed in. Clamped to the slider's [0.5, 8] range.
  private func zoomMag(camDist: Double, radius: Double, fovDeg: Double) -> Double {
    let t = tan(fovDeg * .pi / 360.0)
    guard camDist > 1e-4, t > 1e-6, radius > 0 else { return 1.0 }
    return min(max(radius / (camDist * t), 0.5), 8.0)
  }

  private func setOutlineColor(_ color: Color) {
    engine.runCommand(
      "set_color _outlinecol, \(rgb01List(color))\nset metal_outline_color, _outlinecol")
  }

  @ViewBuilder
  private func sceneRow<Content: View>(
    _ label: String, help: String = "", @ViewBuilder _ content: () -> Content
  ) -> some View {
    HStack(spacing: 6) {
      Text(label)
        .font(.system(size: 10))
        .foregroundColor(PanelTheme.textColor)
        .frame(width: 110, alignment: .leading)
      content()
      Spacer(minLength: 0)
      if !help.isEmpty { HelpButton(text: help) }
    }
  }
}

// MARK: - Object panel JSON payload

/// Decoded shape of the tempfile written by `appkit_inspector.poll_panel`.
///
/// Extracted from the function body so `@testable import mMol` can decode it
/// in unit tests.
struct PanelPayload: Decodable {
  let objects: [String]
  let selections: [String]
  let enabled: [String]
  let sel_counts: [String: Int]
  let nstate: [String: Int]?
  // Group tree (#255). OPTIONAL on purpose: a non-optional field makes
  // the whole decode fail against an older bundled mmol/inspector.py,
  // which would freeze the panel on its last list rather than degrade
  // to a flat one.
  let groups: [String]?
  let maps: [String]?
  let parent: [String: String]?
}

// MARK: - PyMOLEngine extensions for object polling

extension PyMOLEngine {
  /// Parse the object-list JSON (written by appkit_inspector.poll_panel to a
  /// temp file; the feedback line is just the "OBJPANEL:ready" trigger) and
  /// update the objects array. Called by the existing pollFeedback timer.
  /// File-based to avoid the ~1KB feedback-line cap splitting the payload —
  /// inline, the list truncated at ~16 objects (so the panel froze on the stale
  /// list after e.g. `split_states`) and the continuation leaked to the log (#231).
  func parseObjectPanelFeedback(_ line: String) {
    guard line.hasPrefix("OBJPANEL:") else { return }
    // Taken before the latest visibility click: would undo it until the next poll.
    guard objectPollEpoch == objectListEpoch else { return }
    let path = TempChannel.path(TempChannel.Stem.objectPanel)
    guard let data = FileManager.default.contents(atPath: path) else { return }

    guard let payload = try? JSONDecoder().decode(PanelPayload.self, from: data) else {
      return
    }

    let enabledSet = Set(payload.enabled)
    // Enabled *objects* only — a selection's pink markers don't put a row in
    // the sequence strip, so folding them in here would trigger pointless
    // re-fetches every time a selection is (de)selected. See #380.
    let enabledObjects = Set(payload.objects.filter { enabledSet.contains($0) })
    let groupSet = Set(payload.groups ?? [])
    let mapSet = Set(payload.maps ?? [])
    let parentMap = payload.parent ?? [:]
    var entries: [ObjectEntry] = []

    for name in payload.objects {
      entries.append(
        ObjectEntry(
          id: "obj_\(name)",
          name: name,
          isEnabled: enabledSet.contains(name),
          isSelection: false,
          atomCount: nil,
          stateCount: max(payload.nstate?[name] ?? 1, 1),
          isGroup: groupSet.contains(name),
          parent: parentMap[name],
          isMap: mapSet.contains(name)
        ))
    }

    for name in payload.selections {
      entries.append(
        ObjectEntry(
          id: "sel_\(name)",
          name: name,
          isEnabled: enabledSet.contains(name),
          isSelection: true,
          atomCount: payload.sel_counts[name]
        ))
    }

    let epoch = objectListEpoch
    DispatchQueue.main.async {
      // A visibility click that landed after the check above (this block
      // runs later) would be undone by this older snapshot: drop it.
      guard self.objectListEpoch == epoch else { return }
      // Guard: the ~500ms poll usually returns the same object list;
      // re-assigning an equal array still fires @Published and re-renders
      // the panel (resetting open menus). Only assign on real changes.
      if self.objects != entries { self.objects = entries }
      // The sequence strip rows only enabled objects (#380). This is the
      // core's own view of what is enabled, so it catches every path that
      // can change it and cannot race the optimistic checkbox flip. Track
      // the set even while the panel is hidden — `sequenceVisible`'s
      // false->true edge does its own fetch, and a stale tracker would
      // otherwise suppress the next real change.
      if self.lastSequenceEnabled != enabledObjects {
        self.lastSequenceEnabled = enabledObjects
        if self.sequenceVisible { self.fetchSequences() }
      }
      // Every tick, changed or not: the guarded assignments above publish
      // only changes.
      self.panelPolled.send()
    }
  }
}

// MARK: - A/S/H/L menus

// The A/S/H/L/C buttons open the opaque, square desktop-PyMOL menus
// (ClassicMenu.swift). The A button and an object row's right-click open the
// same A menu.

/// The A menu for `entry`'s kind of row (object, selection, group, map).
private func classicActionItems(for entry: ObjectEntry, engine: PyMOLEngine) -> [ClassicMenuItem] {
  classicActionItems(
    actionMenuItems(isSelection: entry.isSelection, isGroup: entry.isGroup, isMap: entry.isMap),
    name: entry.name, engine: engine)
}

/// The A menu (`actionMenuItems`) as classic menu items.
private func classicActionItems(
  _ items: [ActionMenuItem], name: String,
  engine: PyMOLEngine
) -> [ClassicMenuItem] {
  items.map { item -> ClassicMenuItem in
    switch item {
    case .action(let label, let key):
      return .action(label) { runActionCommand(key, name: name, engine: engine) }
    case .separator:
      return .separator
    case .submenu(let label, let children):
      return .submenu(label, classicActionItems(children, name: name, engine: engine))
    case .moveToGroup(let label):
      let groups = engine.objects.filter { $0.isGroup && $0.name != name }
      var kids: [ClassicMenuItem] = groups.map { g in
        .action(g.name) {
          engine.runCommand("group \(g.name), \(name), action=add", naming: g.name, name)
        }
      }
      if !groups.isEmpty { kids.append(.separator) }
      kids.append(.action("New Group…") { engine.pendingGroupFor = name })
      if engine.objects.first(where: { $0.name == name })?.parent != nil {
        kids.append(.separator)
        kids.append(
          .action("Remove from Group") {
            engine.runCommand("ungroup \(name)", naming: name)
          })
      }
      return .submenu(label, kids)
    case .copyToObject(let label):
      var kids: [ClassicMenuItem] = [.action("New Object…") { engine.pendingCopyToNew = name }]
      let targets = engine.objects.filter {
        !$0.isSelection && !$0.isGroup && $0.name != name
      }
      if !targets.isEmpty {
        kids.append(.separator)
        kids += targets.map { t in
          .action(t.name) {
            engine.runCommand(
              copyToObjectCommand(sele: name, target: t.name),
              naming: name, t.name)
          }
        }
      }
      return .submenu(label, kids)
    case .alignToMolecule(let label):
      let targets = engine.objects.filter { !$0.isSelection && $0.name != name }
      return .submenu(
        label,
        targets.isEmpty
          ? [.disabled("(no other molecules)")]
          : targets.map { t in
            .action(t.name) { runAlignCommand(mobile: name, target: t.name, engine: engine) }
          })
    case .alignToSelection(let label):
      let targets = engine.objects.filter { $0.isSelection && $0.name != name }
      return .submenu(
        label,
        targets.isEmpty
          ? [.disabled("(no selections)")]
          : targets.map { t in
            .action(t.name) { runAlignCommand(mobile: name, target: t.name, engine: engine) }
          })
    }
  }
}

extension RepMenuButton {
  fileprivate var showItems: [ClassicMenuItem] {
    let n = name
    var items: [ClassicMenuItem] = [
      .submenu(
        "side chains",
        [
          .action("as sticks") {
            engine.runCommand(
              "show sticks, (\(n)) and \(kSidechainPred) and not hydro; set cartoon_side_chain_helper, 1, \(n)",
              naming: n)
          },
          .action("as lines") {
            engine.runCommand(
              "show lines, (\(n)) and \(kSidechainPred) and not hydro; set cartoon_side_chain_helper, 1, \(n)",
              naming: n)
          },
          .action("as spheres") {
            engine.runCommand("show spheres, (\(n)) and sidechain", naming: n)
          },
        ]),
      .action("hydrogens") { engine.runCommand("show sticks, (\(n)) and hydro", naming: n) },
      .action("all states") { engine.setAllStates(n, true) },
      .separator,
    ]
    for opt in showHideOptions {
      if opt.label == "---" {
        items.append(.separator)
      } else if let rep = opt.rep, rep != "everything" {
        items.append(.action(opt.label) { engine.runCommand("show \(rep), \(n)", naming: n) })
      }
    }
    return items
  }

  fileprivate var hideItems: [ClassicMenuItem] {
    let n = name
    var items: [ClassicMenuItem] = [
      .action("side chains") {
        engine.runCommand(
          "hide sticks, (\(n)) and \(kSidechainPred); hide lines, (\(n)) and \(kSidechainPred)",
          naming: n)
      },
      .submenu(
        "hydrogens",
        [
          .action("all") { engine.runCommand("hide (\(n) and hydro)", naming: n) },
          .action("nonpolar") {
            engine.runCommand("hide (\(n) and hydro and (elem C extend 1))", naming: n)
          },
        ]),
      .action("all states") { engine.setAllStates(n, false) },
      .separator,
    ]
    for opt in showHideOptions {
      if opt.label == "---" {
        items.append(.separator)
      } else if let rep = opt.rep {
        items.append(.action(opt.label) { engine.runCommand("hide \(rep), \(n)", naming: n) })
      }
    }
    return items
  }

  fileprivate var labelItems: [ClassicMenuItem] {
    let n = name
    return labelOptions.compactMap { opt -> ClassicMenuItem? in
      if opt.label == "---" { return .separator }
      guard let expr = opt.expr else { return nil }
      return .action(opt.label) {
        engine.runCommand(expr.isEmpty ? "label \(n)" : "label \(n), \(expr)", naming: n)
      }
    }
  }
}
