// ObjectCommands.swift — Copy / Extract command builders and the name checks
// that guard every command the object panel builds.

import Foundation

// MARK: - Copy / Extract command builders (#461, #468)

/// How far an Extract reaches beyond the selection itself.
///
/// The widened variants copy the bonded neighbourhood into the new object but
/// still remove only the *selected* atoms from the source — that asymmetry is
/// upstream's (`modules/pymol/menu.py:47`), and it is what makes "extend 1"
/// useful for pulling a ligand out with the residues it touches.
enum ExtractScope: String, CaseIterable {
  case object
  case extend1
  case byresExtend1

  /// The selection the new object is built from.
  fileprivate func createSelection(for sele: String) -> String {
    switch self {
    case .object: return "(\(sele))"
    case .extend1: return "((\(sele)) extend 1)"
    case .byresExtend1: return "(byres ((\(sele)) extend 1))"
    }
  }
}

/// Command for "Copy to Object ▸ <existing object>" (#461).
///
/// `copy_to` merges the atoms into an existing object and renames chain/segi/ID
/// so they cannot collide with what is already there; the source object keeps
/// its atoms. `quiet=0` so the " Copied N atoms to object X" line lands in the
/// console feed — the panel gives no other confirmation that anything happened.
///
/// Wrapped in a python block to undo one thing `copy_to` does on its own: it
/// disables every object the selection lives in (`modules/pymol/editing.py:3309`)
/// so the merged result stands out. On the desktop that is a momentary thing in a
/// menu; here the source's row stays in the panel with its checkbox silently
/// clearing, right after a menu item that promises a *copy* — it reads as the
/// copy having eaten the original. So the enabled sources are captured first and
/// switched back on afterwards. Sources that were already hidden stay hidden.
func copyToObjectCommand(sele: String, target: String) -> String {
  let capture =
    "_on = [o for o in cmd.get_object_list(\"(\(sele))\") "
    + "if o in cmd.get_names(\"objects\", enabled_only=1)]"
  let copy = "cmd.copy_to(\"\(target)\", \"(\(sele))\", zoom=0, quiet=0)"
  let restore = "[cmd.enable(o) for o in _on]"
  return "python\n\(capture); \(copy); \(restore)\npython end"
}

/// Command for "Copy to Object ▸ New Object…" (#461).
///
/// Upstream is `cmd.create(None, sele, zoom=0)`, letting PyMOL auto-name the
/// result; here the name comes from the modal instead, so the copy can be
/// labelled at the moment it is made. `zoom=0` keeps the camera still — the new
/// object sits exactly on top of the atoms it was copied from, so framing it
/// would be a jump to nowhere.
func copyToNewObjectCommand(sele: String, name: String) -> String {
  "create \(name), (\(sele)), zoom=0"
}

/// Command for "Extract ▸ …" (#468) — move the selected atoms out of their
/// parent object into a new one.
///
/// Routed through a `python` block rather than the `extract` console keyword
/// because that keyword needs a name and this deliberately passes `None`, which
/// makes PyMOL's own `get_unused_name("obj")` mint `obj01`, `obj02`, … exactly
/// as the desktop menu item does.
///
/// Two notes on the widened variants:
///
/// - they must spell the extract selection out (`extract="<sele>"`) instead of
///   the boolean `extract=1`, because the atoms to remove are a *subset* of the
///   atoms to copy. `cmd.create` calls that form deprecated and says so on the
///   console; upstream's menu has the same wart and the alternative is losing
///   the feature.
/// - the trailing `delete` drops the selection itself. `cmd.extract` only
///   cleans up its own temporary, so the user's named selection would survive
///   with every one of its atoms gone — a 0-atom row in the panel that still
///   looks live. This mirrors what "Remove Atoms" already does.
func extractCommand(sele: String, scope: ExtractScope) -> String {
  let create: String
  switch scope {
  case .object:
    // Selection == extraction target, so `cmd.extract` (which is `create`
    // with the boolean `extract=1`) covers it without the deprecated form.
    create = "cmd.extract(None, \"\(scope.createSelection(for: sele))\", zoom=0)"
  case .extend1, .byresExtend1:
    create =
      "cmd.create(None, \"\(scope.createSelection(for: sele))\", "
      + "extract=\"\(sele)\", zoom=0)"
  }
  return "python\n\(create); cmd.delete(\"\(sele)\")\npython end"
}

/// Mirror of PyMOL's `ExecutiveGetUnusedName(G, "obj")`
/// (`src/pymol/executive/executive.cc:3553` → `ExecutiveMakeUnusedName`, pattern `%02d`
/// starting at 1): the first `objNN` no existing name has taken.
///
/// Computed here rather than asked of the engine so the "New Object…" modal can
/// be prefilled the instant it opens. `existing` should be every name PyMOL
/// knows — objects, selections and groups all share one namespace.
func defaultNewObjectName(existing: [String]) -> String {
  let taken = Set(existing)
  var n = 1
  while true {
    let candidate = String(format: "obj%02d", n)
    if !taken.contains(candidate) { return candidate }
    n += 1
  }
}

/// PyMOL's own legal-name character set: A–Z, a–z, 0–9 and `+ - . ^ _`
/// (`ObjectMakeValidName`, src/pymol/objects/pymol_object.cc).
///
/// Enforced up front rather than left to the engine because the engine's
/// response to an illegal name is to quietly rewrite it — type `my obj!` and you
/// get an object called `my_obj` — so the name in the panel is not the name that
/// was asked for. A comma is worse than cosmetic: it splits the command into a
/// different argument list and `create foo, bar, (sele), zoom=0` throws a Python
/// traceback into the console feed.
func isLegalObjectName(_ name: String) -> Bool {
  !name.isEmpty
    && name.allSatisfy { c in
      guard c.isASCII else { return false }
      return c.isLetter || c.isNumber || "+-.^_".contains(c)
    }
}

/// Whether `token` can be interpolated into a PyMOL command without splitting
/// it (#531), for strings that are NOT object names and so may legitimately
/// fall outside isLegalObjectName: scene names (spaces allowed) and a group
/// name the user is typing. Refuses exactly what can end or split a command: control
/// characters and every line break `cmd.do` splits on (str.splitlines), `;`
/// (PyMOL's command separator), `,` (its argument separator), and `\` -- a
/// line ending in a backslash is a continuation, so PyMOL would merge the
/// NEXT command into this one and lose both.
func isCommandSafeToken(_ token: String) -> Bool {
  !token.isEmpty
    && token.unicodeScalars.allSatisfy { u in
      switch u.value {
      case 0x00...0x1F, 0x7F, 0x85, 0x2028, 0x2029: return false
      case 0x3B, 0x2C, 0x5C: return false  // ; , backslash
      default: return true
      }
    }
}

extension PyMOLEngine {
  /// Run a command built by interpolating NAMES into PyMOL command language
  /// (#531).
  ///
  /// Every object, selection and group name the Inspector shows came from
  /// PyMOL, and with `validate_object_names` on (the default)
  /// ObjectMakeValidName keeps them inside `[A-Za-z0-9+-.^_]`. With it off,
  /// a name can carry a line break, `;` or `,`, and interpolated raw it would
  /// split the command -- running its tail as commands of its own. PyMOL
  /// command language has no general quoting, so rather than escape each of
  /// ~60 builders, refuse the command when a name is outside PyMOL's own
  /// alphabet (`naming:`), or when a non-object token -- a scene name, a
  /// group name the user typed -- contains anything that can split a
  /// command (`tokens:`, see isCommandSafeToken). The console says why
  /// without echoing the name. Every legal name passes, so in the ordinary
  /// case every command runs exactly as it did.
  func runCommand(_ command: String, naming names: String..., tokens: [String] = []) {
    guard names.allSatisfy(isLegalObjectName),
      tokens.allSatisfy(isCommandSafeToken)
    else {
      refuseUnsafeCommand()
      return
    }
    runCommand(command)
  }

  func refuseUnsafeCommand() {
    runCommand(
      "python\nprint(' Inspector: skipped an action -- a name contains "
        + "characters a PyMOL command cannot carry (possible when "
        + "validate_object_names is off). cmd.get_names() lists the objects; "
        + "rename the odd one with cmd.set_name.')"
        + "\npython end")
  }
}

/// Whether `name` can be given to a brand-new object.
///
/// Legal, and not already taken: `create` against a name that already exists
/// does *nothing at all* — no new object, no error, no console line — so an
/// unchecked collision turns "New Object…" into a button that silently fails.
func canNameNewObject(_ name: String, existing: [String]) -> Bool {
  isLegalObjectName(name) && !existing.contains(name)
}
