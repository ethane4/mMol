// MousePanel.swift — the mouse legend under the object panel: current mouse
// mode, button/modifier action grid, and selection mode, as desktop PyMOL draws it.

import SwiftUI

// MARK: - Mouse mode definitions

/// Action code abbreviations matching the Python CODE dict.
private let actionCodes: [String: String] = [
  "rota": "Rota", "move": "Move", "movz": "MovZ", "clip": "Clip",
  "rotz": "RotZ", "clpn": "ClpN", "clpf": "ClpF",
  "lb": " lb ", "mb": " mb ", "rb": " rb ",
  "+lb": "+lb ", "+mb": "+mb ", "+rb": "+rb ",
  "pkat": "PkAt", "pkbd": "PkBd", "rotf": "RotF",
  "torf": "TorF", "movf": "MovF", "orig": "Orig",
  "+lbx": "+lBx", "-lbx": "-lBx", "lbbx": "lbBx",
  "none": " -- ", "cent": "Cent", "pktb": "PkTB",
  "slab": "Slab", "movs": "MovS", "pk1": "Pk1 ",
  "mova": "MovA", "menu": "Menu", "sele": "Sele",
  "+/-": "+/- ", "+box": "+Box", "-box": "-Box",
  "mvsz": "MvSZ", "clik": "Clik", "mvoz": "MvOZ",
  "movo": "MovO", "roto": "RotO", "drgm": "DrgM",
  "rotv": "RotV", "movv": "MovV", "mvvz": "MvVZ",
  "drgo": "DrgO", "mvfz": "MvFZ", "mvaz": "MvAZ",
  "rotl": "RotL", "movl": "MovL", "mvzl": "MvzL",
  "imsz": "IMSZ", "imvz": "IMvZ", "box": " Box",
  "irtz": "IRtZ",
  "rotd": "RotD", "movd": "MovD", "mvdz": "MvDZ",
]

/// Tooltips for the legend's codes and row labels; codes not listed get none.
private let legendHelp: [String: String] = [
  "Rota": "Rotate the view", "Move": "Move the view", "MovZ": "Move closer or farther (zoom)",
  "Clip": "Move the clipping planes", "RotZ": "Rotate about the screen axis",
  "ClpN": "Move the near clipping plane", "ClpF": "Move the far clipping plane",
  "PkAt": "Pick an atom", "PkBd": "Pick a bond", "PkTB": "Pick a torsion bond",
  "Pk1": "Pick atom pk1", "RotF": "Rotate a fragment", "TorF": "Twist a fragment's torsion",
  "MovF": "Move a fragment", "MovA": "Move an atom", "MovO": "Move an object",
  "RotO": "Rotate an object", "MvOZ": "Move an object closer or farther",
  "Orig": "Set the rotation origin", "Cent": "Center on the clicked atom",
  "Slab": "Change the slab thickness", "MovS": "Move the slab",
  "MvSZ": "Move the slab and zoom", "Menu": "Open the context menu", "Sele": "Select",
  "+/-": "Toggle in the selection", "+Box": "Add a box of atoms to the selection",
  "-Box": "Remove a box of atoms from the selection", "--": "No action",
  "Buttons": "Left, middle and right mouse buttons, and the scroll wheel",
  "& Keys": "No modifier key held", "Shft": "With Shift held", "Ctrl": "With Control held",
  "CtSh": "With Control and Shift held", "Clik": "A single click",
]

/// A single mouse button/modifier binding: (button, modifier, action).
private typealias Binding = (btn: String, mod: String, act: String)

/// Mouse mode: key, display name, and bindings list.
private struct MouseModeDefinition {
  let key: String
  let displayName: String
  let bindings: [Binding]
}

/// Index map: (button, modifier) -> slot in a 22-element action array.
/// Slots 0-11: drag actions (L/M/R x none/shft/ctrl/ctsh)
/// Slots 12-15: wheel (none/shft/ctrl/ctsh)
/// Slots 16-18: double-click (L/M/R)
/// Slots 19-21: single-click (L/M/R)
private let buttonModToIndex: [String: Int] = [
  "l:none": 0, "m:none": 1, "r:none": 2,
  "l:shft": 3, "m:shft": 4, "r:shft": 5,
  "l:ctrl": 6, "m:ctrl": 7, "r:ctrl": 8,
  "l:ctsh": 9, "m:ctsh": 10, "r:ctsh": 11,
  "w:none": 12, "w:shft": 13, "w:ctrl": 14, "w:ctsh": 15,
  "double_left:none": 16, "double_middle:none": 17, "double_right:none": 18,
  "single_left:none": 19, "single_middle:none": 20, "single_right:none": 21,
]

/// Build a 22-element action array from a list of bindings.
private func buildModeArray(_ bindings: [Binding]) -> [String] {
  var mode = Array(repeating: "none", count: 22)
  for b in bindings {
    let key = "\(b.btn):\(b.mod)"
    if let idx = buttonModToIndex[key] {
      mode[idx] = b.act.lowercased()
    }
  }
  return mode
}

/// Display code for an action string.
private func codeFor(_ action: String) -> String {
  actionCodes[action.lowercased()] ?? String(action.prefix(4))
}

// MARK: - Mode data (mirrors pymol.controlling.mode_dict)

private let allModes: [MouseModeDefinition] = [
  MouseModeDefinition(
    key: "three_button_viewing", displayName: "3-Button Viewing",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "move"), ("r", "none", "movz"),
      ("l", "shft", "+box"), ("m", "shft", "-box"), ("r", "shft", "clip"),
      ("l", "ctrl", "move"), ("m", "ctrl", "pkat"), ("r", "ctrl", "pk1"),
      ("l", "ctsh", "sele"), ("m", "ctsh", "orig"), ("r", "ctsh", "clip"),
      ("w", "none", "slab"), ("w", "shft", "movs"), ("w", "ctrl", "mvsz"), ("w", "ctsh", "movz"),
      ("double_left", "none", "menu"), ("double_middle", "none", "none"),
      ("double_right", "none", "pkat"),
      ("single_left", "none", "+/-"), ("single_middle", "none", "cent"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "three_button_editing", displayName: "3-Button Editing",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "move"), ("r", "none", "movz"),
      ("l", "shft", "roto"), ("m", "shft", "movo"), ("r", "shft", "mvoz"),
      ("l", "ctrl", "torf"), ("m", "ctrl", "+/-"), ("r", "ctrl", "pktb"),
      ("l", "ctsh", "mova"), ("m", "ctsh", "orig"), ("r", "ctsh", "clip"),
      ("w", "none", "slab"), ("w", "shft", "movs"), ("w", "ctrl", "mvsz"), ("w", "ctsh", "movz"),
      ("double_left", "none", "torf"), ("double_middle", "none", "drgm"),
      ("double_right", "none", "pktb"),
      ("single_left", "none", "pkat"), ("single_middle", "none", "cent"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "three_button_motions", displayName: "3-Button Motions",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "move"), ("r", "none", "movz"),
      ("l", "shft", "rotv"), ("m", "shft", "movv"), ("r", "shft", "mvvz"),
      ("l", "ctrl", "torf"), ("m", "ctrl", "pkat"), ("r", "ctrl", "pktb"),
      ("l", "ctsh", "mova"), ("m", "ctsh", "orig"), ("r", "ctsh", "clip"),
      ("w", "none", "slab"), ("w", "shft", "movs"), ("w", "ctrl", "mvsz"), ("w", "ctsh", "movz"),
      ("double_left", "none", "torf"), ("double_middle", "none", "drgm"),
      ("double_right", "none", "pktb"),
      ("single_left", "none", "pkat"), ("single_middle", "none", "cent"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "three_button_lights", displayName: "3-Button Lights",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "move"), ("r", "none", "movz"),
      ("l", "shft", "rotl"), ("m", "shft", "movl"), ("r", "shft", "mvzl"),
      ("l", "ctrl", "none"), ("m", "ctrl", "none"), ("r", "ctrl", "none"),
      ("l", "ctsh", "none"), ("m", "ctsh", "none"), ("r", "ctsh", "none"),
      ("w", "none", "slab"), ("w", "shft", "movs"), ("w", "ctrl", "mvsz"), ("w", "ctsh", "movz"),
      ("double_left", "none", "none"), ("double_middle", "none", "none"),
      ("double_right", "none", "none"),
      ("single_left", "none", "none"), ("single_middle", "none", "cent"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "three_button_maestro", displayName: "3-Button Maestro",
    bindings: [
      ("l", "none", "box"), ("m", "none", "rota"), ("r", "none", "move"),
      ("l", "shft", "+box"), ("m", "shft", "-box"), ("r", "shft", "clip"),
      ("l", "ctrl", "+/-"), ("m", "ctrl", "irtz"), ("r", "ctrl", "pk1"),
      ("l", "ctsh", "sele"), ("m", "ctsh", "orig"), ("r", "ctsh", "clip"),
      ("w", "none", "imvz"), ("w", "shft", "movs"), ("w", "ctrl", "none"), ("w", "ctsh", "slab"),
      ("double_left", "none", "menu"), ("double_middle", "none", "none"),
      ("double_right", "none", "pkat"),
      ("single_left", "none", "sele"), ("single_middle", "none", "cent"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "two_button_viewing", displayName: "2-Button Viewing",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "none"), ("r", "none", "movz"),
      ("l", "shft", "pk1"), ("m", "shft", "none"), ("r", "shft", "clip"),
      ("l", "ctrl", "move"), ("m", "ctrl", "none"), ("r", "ctrl", "pkat"),
      ("l", "ctsh", "sele"), ("m", "ctsh", "none"), ("r", "ctsh", "cent"),
      ("w", "none", "none"), ("w", "shft", "none"), ("w", "ctrl", "none"), ("w", "ctsh", "none"),
      ("double_left", "none", "menu"), ("double_middle", "none", "none"),
      ("double_right", "none", "pkat"),
      ("single_left", "none", "+/-"), ("single_middle", "none", "none"),
      ("single_right", "none", "menu"),
    ]),
  MouseModeDefinition(
    key: "one_button_viewing", displayName: "1-Button Viewing",
    bindings: [
      ("l", "none", "rota"), ("m", "none", "none"), ("r", "none", "none"),
      ("l", "shft", "none"), ("m", "shft", "none"), ("r", "shft", "none"),
      ("l", "ctrl", "none"), ("m", "ctrl", "none"), ("r", "ctrl", "none"),
      ("l", "ctsh", "none"), ("m", "ctsh", "none"), ("r", "ctsh", "none"),
      ("w", "none", "none"), ("w", "shft", "none"), ("w", "ctrl", "none"), ("w", "ctsh", "none"),
      ("double_left", "none", "menu"), ("double_middle", "none", "none"),
      ("double_right", "none", "none"),
      ("single_left", "none", "cent"), ("single_middle", "none", "none"),
      ("single_right", "none", "menu"),
    ]),
]

/// The default ring of (real PyMOL) modes users cycle through.
private let modeRing = [
  "three_button_viewing", "three_button_editing",
  "three_button_motions", "three_button_lights",
  "three_button_maestro",
]

/// Synthetic mode key for the macOS "Trackpad" gesture mode. Not a PyMOL
/// mouse_mode — selecting it flips engine.trackpadMode (which re-routes the
/// MetalViewport trackpad gestures to rotate/pan/zoom/roll/clip) instead of
/// issuing a `mouse <mode>` command.
private let kTrackpadKey = "trackpad"

/// Gesture-oriented legend shown while Trackpad mode is active (rotate / move /
/// zoom / roll / clip), in place of the L/M/R button matrix.
private struct TrackpadGestureHint: Identifiable {
  let id = UUID()
  let title: String
  let detail: String
}
private let trackpadGestureHints: [TrackpadGestureHint] = [
  .init(title: "Rotate", detail: "Drag"),
  .init(title: "Move", detail: "Two-finger sideways"),
  .init(title: "Slab", detail: "Two-finger up/down"),
  .init(title: "Zoom", detail: "Pinch"),
  .init(title: "Roll", detail: "Two-finger twist"),
  .init(title: "Clip", detail: "Shift + two-finger"),
]

private let selectionModeNames = [
  "Atoms", "Residues", "Chains", "Segments",
  "Objects", "Molecules", "C-alphas",
]

// MARK: - MousePanel View

struct MousePanel: View {
  @EnvironmentObject var engine: PyMOLEngine

  // Persisted across launches. Defaults to 3-Button Viewing (index 0).
  @AppStorage("macMouseModeIndex") private var selectedModeIndex: Int = -1
  @State private var selectionModeIndex: Int = 1  // Residues default

  // The keys shown in the mode picker.
  //
  // Only the two modes that are fully, cleanly implemented in the
  // Metal viewport: 3-Button Viewing (rotate/zoom/clip via drags) and the
  // synthetic Trackpad gesture mode. Every other classic mode is hidden here —
  // on this backend a mode only reprograms the drag rows (clicks are hardcoded),
  // and their distinguishing behavior is either unreachable (the middle-button
  // column — no middle button on most Macs) or lacks its visual. Notably
  // Maestro's signature bare-left box-select DOES select atoms but draws no
  // rubber-band marquee: the Metal path (SceneRenderMetal) forces
  // internal_gui/feedback off and never renders the Ortho overlay CGO that
  // OrthoDrawLoop emits, so the drag has no on-screen feedback. Offering those
  // modes would imply behavior the app doesn't visibly deliver.
  private let displayRing = ["three_button_viewing", kTrackpadKey]

  /// Effective selected index, clamped into the display ring. The persisted
  /// default (-1) resolves to 3-Button Viewing.
  private var effectiveModeIndex: Int {
    displayRing.indices.contains(selectedModeIndex) ? selectedModeIndex : 0
  }

  private func displayName(_ key: String) -> String {
    if key == kTrackpadKey { return "Trackpad" }
    return allModes.first { $0.key == key }?.displayName ?? key
  }

  private var isTrackpad: Bool {
    displayRing.indices.contains(effectiveModeIndex)
      && displayRing[effectiveModeIndex] == kTrackpadKey
  }

  private var currentMode: MouseModeDefinition {
    let ringKey =
      displayRing.indices.contains(effectiveModeIndex)
      ? displayRing[effectiveModeIndex] : modeRing[0]
    return allModes.first { $0.key == ringKey } ?? allModes[0]
  }

  private var modeArray: [String] {
    buildModeArray(currentMode.bindings)
  }

  var body: some View {
    ClassicMouseLegend(
      modeName: displayName(displayRing[effectiveModeIndex]),
      isTrackpad: isTrackpad,
      mode: modeArray,
      code: codeFor,
      selectionMode: selectionModeNames.indices.contains(selectionModeIndex)
        ? selectionModeNames[selectionModeIndex] : "Atoms",
      cycleMode: { selectedModeIndex = (effectiveModeIndex + 1) % displayRing.count },
      cycleSelection: {
        selectionModeIndex = (selectionModeIndex + 1) % selectionModeNames.count
      }
    )
    .onAppear { applyMode(effectiveModeIndex) }
    .onChange(of: selectedModeIndex) { applyMode($0) }
    .onChange(of: selectionModeIndex) { newIndex in
      engine.runCommand("set mouse_selection_mode, \(newIndex)")
    }
  }

  /// Apply the mode at `index`: flip engine.trackpadMode and, for real PyMOL
  /// modes, issue the `mouse <key>` command. The synthetic Trackpad entry only
  /// flips the flag (no PyMOL mouse_mode change).
  private func applyMode(_ index: Int) {
    guard displayRing.indices.contains(index) else { return }
    let key = displayRing[index]
    if key == kTrackpadKey {
      engine.trackpadMode = true
    } else {
      engine.trackpadMode = false
      engine.runCommand("mouse \(key)")
    }
  }
}

// MARK: - Legend

/// Desktop PyMOL's internal-GUI mouse legend: monospaced text in the core's
/// green / red / blue, mode name and selection mode clickable to cycle. The
/// click row is the Metal viewport's fixed Sele / -- / Menu (clicks bypass the
/// per-mode button matrix), not the per-mode values.
private struct ClassicMouseLegend: View {
  let modeName: String
  let isTrackpad: Bool
  let mode: [String]
  let code: (String) -> String
  let selectionMode: String
  let cycleMode: () -> Void
  let cycleSelection: () -> Void
  @EnvironmentObject private var playback: PlaybackState

  private let labelW: CGFloat = 56
  private let cellW: CGFloat = 38

  var body: some View {
    VStack(alignment: .leading, spacing: 0) {
      Button(action: cycleMode) {
        HStack(spacing: 6) {
          t("Mouse Mode", ClassicUI.legendGreen)
          t(modeName, ClassicUI.legendRed)
        }
        .contentShape(Rectangle())
      }
      .buttonStyle(ClassicPressStyle())
      .help("Click to switch mouse mode")

      if isTrackpad {
        ForEach(trackpadGestureHints) { h in
          row(h.title, [h.detail], labelColor: ClassicUI.legendBlue, wide: true)
        }
      } else {
        row(
          "Buttons", ["L", "M", "R", "Wheel"],
          labelColor: ClassicUI.legendRed, cellColor: ClassicUI.legendBlue)
        row(
          "& Keys", [code(mode[0]), code(mode[1]), code(mode[2]), code(mode[12])],
          labelColor: ClassicUI.legendRed)
        row("Shft", [code(mode[3]), code(mode[4]), code(mode[5]), code(mode[13])])
        row("Ctrl", [code(mode[6]), code(mode[7]), code(mode[8]), code(mode[14])])
        row("CtSh", [code(mode[9]), code(mode[10]), code(mode[11]), code(mode[15])])
        row("Clik", ["Sele", "Cent", "Menu"])
      }

      Button(action: cycleSelection) {
        HStack(spacing: 6) {
          t("Selecting", ClassicUI.legendGreen)
          t(selectionMode, ClassicUI.legendRed)
        }
        .contentShape(Rectangle())
      }
      .buttonStyle(ClassicPressStyle())
      .help("Click to change what a click selects")

      HStack(spacing: 6) {
        t("State", ClassicUI.legendGreen)
        t(
          String(format: "%4d/%4d", playback.currentFrame, playback.frameCount),
          ClassicUI.legendText)
      }
    }
    .padding(.horizontal, 4)
    .padding(.vertical, 3)
    .frame(maxWidth: .infinity, alignment: .leading)
    .background(ClassicUI.panelBg)
  }

  private func t(_ s: String, _ c: Color) -> some View {
    Text(s).font(ClassicUI.font).foregroundColor(c).lineLimit(1)
  }

  private func row(
    _ label: String, _ cells: [String],
    labelColor: Color = ClassicUI.legendBlue,
    cellColor: Color = ClassicUI.legendText,
    wide: Bool = false
  ) -> some View {
    HStack(spacing: 0) {
      t(label, labelColor).frame(width: labelW, alignment: .trailing)
        .help(legendHelp[label] ?? "")
      ForEach(Array(cells.enumerated()), id: \.offset) { _, c in
        t(c, cellColor)
          .frame(width: wide ? nil : cellW, alignment: wide ? .leading : .center)
          .padding(.leading, wide ? 6 : 0)
          .help(legendHelp[c.trimmingCharacters(in: .whitespaces)] ?? "")
      }
    }
  }
}

// MARK: - Preview
