// PanelLayout.swift — the persisted panel-layout contract of the layout in
// ContentView (#331 console default, #332 persistence).
//
// Two responsibilities, both deliberately free of SwiftUI:
//
//   1. The UserDefaults key namespace (`mmol.panels.*`). The layout reads
//      these keys through @AppStorage, so pane sizes survive a relaunch.
//   2. The sizing arithmetic. A size the USER chose is stored as a FRACTION of
//      the window dimension, never as absolute points: a height that was sensible
//      in a 1600pt window must not be restored into a 700pt one, where it would
//      leave no viewport. Fractions are clamped on the way in AND on the way out.
//      The untouched DEFAULT is the opposite — absolute points — so a console
//      nobody has resized looks the same on a laptop and on a 6K display.
//
// Everything here is a pure function of its arguments so PanelLayoutTests can
// pin the boundaries without a window.

import CoreGraphics
import Foundation

enum PanelLayout {

  // MARK: - Keys

  private static let ns = "mmol.panels."

  /// The user's top band (console + button grid) height, as a fraction of the window height. Absent
  /// until they resize it, which is what selects the absolute default (#331).
  static let consoleFracKey = ns + "consoleFrac"
  /// The user's inspector (right column) width, as a fraction of the
  /// window width. Absent until they drag the seam (#350), which is what
  /// selects the absolute default width.
  static let inspectorFracKey = ns + "inspectorFrac"
  /// Sequence strip visible. Written by PyMOLEngine, which owns the flag.
  static let sequenceVisibleKey = ns + "sequenceVisible"

  /// Every key this type defines — the namespace/uniqueness check in the tests
  /// runs off this list, so a new key must be added here too.
  static let allKeys: [String] = [consoleFracKey, inspectorFracKey, sequenceVisibleKey]

  // MARK: - Bounds

  /// #331: the top band's height on a launch with nothing stored. ABSOLUTE
  /// points, not a share of the window, so it matches what the app has always
  /// shown on any display. A size the user drags to IS stored as a fraction
  /// (see `consoleFrac`); only this untouched default is fixed.
  static let macDefaultConsoleHeight: CGFloat = 170
  /// Storable band: storing ~0 would restore an unusable sliver on the next
  /// launch.
  static let minStorableFrac: CGFloat = 0.01
  static let maxStorableFrac: CGFloat = 0.95

  /// Usable minimum: the top band must fit the five-row button grid.
  static let macMinConsoleHeight: CGFloat = 150
  /// Floor under the growth ceiling (see `maxConsoleHeight`): in a window too
  /// short to satisfy the viewport minimum, the console may still take this
  /// share rather than collapsing to nothing.
  static let minCeilingFrac: CGFloat = 0.2
  /// The macOS viewport's own minimum. The console may grow until the viewport
  /// is down to this, which is what keeps #317 ("drag it open to read a long
  /// log") working while still guaranteeing a viewport.
  static let macViewportMinHeight: CGFloat = 360

  /// #350: the macOS inspector column's width on a launch with nothing stored.
  /// ABSOLUTE points for the same reason as the console default — an untouched
  /// inspector looks the same on a laptop and on a 6K display.
  static let macDefaultInspectorWidth: CGFloat = 340
  /// Below this an Objects row is nothing but its fixed chrome — gutter (26) +
  /// five 38pt A/S/H/L/C buttons + chevron/indent — with no room left for even
  /// a short name, so shrinking further only breaks the panel.
  static let macMinInspectorWidth: CGFloat = 280
  /// The macOS viewport's own minimum WIDTH; the inspector may grow until the
  /// viewport is down to this. The horizontal sibling of `macViewportMinHeight`.
  static let macViewportMinWidth: CGFloat = 480

  // MARK: - macOS inspector sizing (#350)

  /// How wide the inspector is allowed to get in a window of `windowWidth`,
  /// given a viewport beside it that needs `viewportMin`. Same bracket shape as
  /// `maxConsoleHeight`: never more than half the window (the viewport must
  /// remain the main event), never less than `minCeilingFrac`.
  static func maxInspectorWidth(
    windowWidth: CGFloat,
    viewportMin: CGFloat = macViewportMinWidth
  ) -> CGFloat {
    guard windowWidth.isFinite, windowWidth > 0 else { return 0 }
    return min(
      max(windowWidth - viewportMin, windowWidth * minCeilingFrac),
      windowWidth * 0.5)
  }

  /// The inspector's width in a window of `windowWidth`.
  ///
  /// `frac` is the persisted share and applies only once the user has dragged
  /// the seam themselves; a missing, zero, negative or non-finite value (an
  /// unset UserDefaults Double reads as 0) means "untouched" and yields the
  /// absolute default. The CEILING WINS when it crosses the usable minimum —
  /// in a very narrow window, better a cramped inspector than no viewport.
  static func inspectorWidth(
    frac: CGFloat, windowWidth: CGFloat,
    maxWidth: CGFloat
  ) -> CGFloat {
    guard windowWidth.isFinite, windowWidth > 0 else { return macMinInspectorWidth }
    let w = (frac.isFinite && frac > 0) ? frac * windowWidth : macDefaultInspectorWidth
    return min(max(w, macMinInspectorWidth), maxWidth)
  }

  /// The fraction to persist for an inspector measured at `width`, clamped into
  /// the storable band so nothing unrestorable can ever be written. nil for a
  /// degenerate window — there is nothing meaningful to store, so don't write.
  static func inspectorFrac(width: CGFloat, windowWidth: CGFloat) -> CGFloat? {
    guard windowWidth.isFinite, windowWidth > 0, width.isFinite else { return nil }
    return min(max(width / windowWidth, minStorableFrac), maxStorableFrac)
  }

  // MARK: - Console sizing

  /// How tall the console is allowed to get in a window of `windowHeight`, given
  /// a pane below it that needs `viewportMin`.
  ///
  /// Normally that's "everything except the viewport's minimum" — 324pt of a
  /// 684pt window — which is generous enough to read a long log in. Two guards
  /// bracket it: never more than 85% of the window (the viewport must remain a
  /// viewport), and never less than `minCeilingFrac`.
  static func maxConsoleHeight(
    windowHeight: CGFloat,
    viewportMin: CGFloat = macViewportMinHeight
  ) -> CGFloat {
    guard windowHeight.isFinite, windowHeight > 0 else { return 0 }
    return min(
      max(windowHeight - viewportMin, windowHeight * minCeilingFrac),
      windowHeight * 0.85)
  }

  /// The console's height in a window of `windowHeight`.
  ///
  /// `frac` is the persisted share, and applies only once the user has sized the
  /// console themselves; a missing, zero, negative or non-finite value (an unset
  /// UserDefaults Double reads as 0) means "untouched" and yields `defaultHeight`
  /// — an absolute size, so an untouched console looks the same on a laptop and
  /// on a 6K display.
  ///
  /// `minHeight` is the usable floor; `maxHeight` the ceiling
  /// (`maxConsoleHeight`).
  /// The CEILING WINS when the two cross — in a very short window, better a
  /// cramped console than no viewport.
  static func consoleHeight(
    frac: CGFloat, windowHeight: CGFloat,
    defaultHeight: CGFloat,
    minHeight: CGFloat, maxHeight: CGFloat
  ) -> CGFloat {
    guard windowHeight.isFinite, windowHeight > 0 else { return minHeight }
    let h = (frac.isFinite && frac > 0) ? frac * windowHeight : defaultHeight
    return min(max(h, minHeight), maxHeight)
  }

  /// The fraction to persist for a console measured at `height`, clamped into
  /// the storable band so nothing unrestorable can ever be written. nil for a
  /// degenerate window — there is nothing meaningful to store, so don't write.
  static func consoleFrac(height: CGFloat, windowHeight: CGFloat) -> CGFloat? {
    guard windowHeight.isFinite, windowHeight > 0, height.isFinite else { return nil }
    return min(max(height / windowHeight, minStorableFrac), maxStorableFrac)
  }
}
