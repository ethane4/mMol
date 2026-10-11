// ClassicUI.swift — the "Classic PyMOL Layout": desktop PyMOL 2.5's window
// arrangement and chrome on top of mMol's Metal viewport (macOS only).
//
// The layout itself is assembled in ContentView (classicMacBody); this file holds
// the pieces it is built from that have no mMol counterpart:
//
//   * ClassicButtonGrid  — the Qt GUI's Reset/Zoom/Orient/… grid beside the
//                          console (modules/pmg_qt/pymol_qt_gui.py).
//   * ClassicTransportBar — the internal GUI's |< < ■ ▶ > >| S ▼ F strip under
//                          the object panel (src/pymol/ui/control.cc).
//   * ClassicBevel       — the raised-button look of the internal GUI's
//                          draw_button (src/pymol/executive/executive.cc), colours verbatim.
//
// Everything else (object panel A/S/H/L/C, mouse legend, console) is mMol's own
// view in the classic look.

import AppKit
import SwiftUI

/// sRGB color with plain components (for PyMOL colours and blended chrome).
struct RGBA: Codable, Equatable {
  var r: Double, g: Double, b: Double, a: Double
  init(_ r: Double, _ g: Double, _ b: Double, _ a: Double = 1) {
    self.r = r
    self.g = g
    self.b = b
    self.a = a
  }
  var color: Color { Color(.sRGB, red: r, green: g, blue: b, opacity: a) }
  var nsColor: NSColor { NSColor(srgbRed: r, green: g, blue: b, alpha: a) }
  /// Linear blend toward `other` by `t` (0…1).
  func blended(with other: RGBA, _ t: Double) -> RGBA {
    RGBA(
      r + (other.r - r) * t, g + (other.g - g) * t,
      b + (other.b - b) * t, a + (other.a - a) * t)
  }
}

extension Color {
  /// Resolve a SwiftUI Color to sRGB components.
  var rgba: RGBA {
    let ns = NSColor(self).usingColorSpace(.sRGB) ?? NSColor(self)
    return RGBA(
      Double(ns.redComponent), Double(ns.greenComponent),
      Double(ns.blueComponent), Double(ns.alphaComponent))
  }
}

enum ClassicUI {
  // PyMOL 2.5's Qt chrome: panel grey, near-white text, the "A" toggle's blue
  // as the accent, and the selection red. Console text is Ortho's TextColor.
  static let chromeBg = RGBA(0.2, 0.2, 0.2)
  static let chromeText = RGBA(0.95, 0.95, 0.95)
  static let accentRGBA = RGBA(0.5, 0.5, 1.0)
  static let selectionRGBA = RGBA(1.0, 0.3, 0.3)
  static let accent = accentRGBA.color

  // Internal-GUI palette (src/pymol/executive/executive.cc ExecutiveDrawPanel / draw_button).
  static let nameEnabled = Color(white: 0.5)
  static let nameDisabled = Color(white: 0.25)
  static let lightEdge = Color(white: 0.6)
  static let darkEdge = Color(white: 0.35)
  static let toggleA = Color(.sRGB, red: 0.5, green: 0.5, blue: 1.0, opacity: 1)  // A, L
  static let toggleS = Color(.sRGB, red: 0.6, green: 0.6, blue: 0.8, opacity: 1)  // S
  static let toggleH = Color(.sRGB, red: 0.4, green: 0.4, blue: 0.6, opacity: 1)  // H
  static let toggleLight = Color(.sRGB, red: 0.7, green: 0.7, blue: 0.9, opacity: 1)
  static let toggleDark = Color(.sRGB, red: 0.3, green: 0.3, blue: 0.5, opacity: 1)
  static let panelBg = chromeBg.color
  static let consoleText = RGBA(0.83, 0.83, 1.0)

  // Mouse-legend colours (layer1/Ortho / ButMode).
  static let legendGreen = Color(.sRGB, red: 0.2, green: 1.0, blue: 0.2, opacity: 1)
  static let legendRed = Color(.sRGB, red: 1.0, green: 0.3, blue: 0.3, opacity: 1)
  static let legendBlue = Color(.sRGB, red: 0.5, green: 0.5, blue: 1.0, opacity: 1)
  static let legendText = Color(white: 0.9)

  static let font = Font.system(size: 12, design: .monospaced)

  /// Fill for the A/S/H/L/C toggle with this letter.
  static func toggleFill(_ letter: String) -> Color {
    switch letter {
    case "S": return toggleS
    case "H": return toggleH
    default: return toggleA
    }
  }

}

/// The "C" toggle's fill: draw_button with no colour paints a four-corner blend —
/// green top-left, yellow top-right, red bottom-left, blue bottom-right.
struct ClassicRainbowFill: View {
  private static func c(_ r: Double, _ g: Double, _ b: Double) -> Color {
    Color(.sRGB, red: r, green: g, blue: b, opacity: 1)
  }

  var body: some View {
    ZStack {
      LinearGradient(
        colors: [Self.c(1, 0.1, 0.1), Self.c(0.1, 0.1, 1)],
        startPoint: .leading, endPoint: .trailing)
      LinearGradient(
        colors: [Self.c(0.1, 1, 0.1), Self.c(1, 1, 0.1)],
        startPoint: .leading, endPoint: .trailing
      )
      .mask(
        LinearGradient(
          colors: [.white, .clear],
          startPoint: .top, endPoint: .bottom))
    }
  }
}

/// The internal GUI's raised button: `fill` with a 1px light top/left edge and a
/// dark bottom/right edge, square corners.
struct ClassicBevel<Fill: View>: View {
  var fill: Fill
  var light: Color = ClassicUI.lightEdge
  var dark: Color = ClassicUI.darkEdge

  var body: some View {
    fill
      .overlay(alignment: .top) { light.frame(height: 1) }
      .overlay(alignment: .leading) { light.frame(width: 1) }
      .overlay(alignment: .bottom) { dark.frame(height: 1) }
      .overlay(alignment: .trailing) { dark.frame(width: 1) }
  }
}

/// A flat dark Qt-style push button, as in the 2.5 button grid.
private struct ClassicGridButtonStyle: ButtonStyle {
  func makeBody(configuration: Configuration) -> some View {
    configuration.label
      .font(.system(size: 13))
      .foregroundColor(Color(white: 0.92))
      .lineLimit(1)
      .frame(maxWidth: .infinity, minHeight: 22)
      .background(Color(white: configuration.isPressed ? 0.32 : 0.22))
      .overlay(RoundedRectangle(cornerRadius: 3).strokeBorder(Color(white: 0.38), lineWidth: 1))
      .clipShape(RoundedRectangle(cornerRadius: 3))
      .contentShape(Rectangle())
  }
}

/// The Qt GUI's control grid beside the console. Commands match pmg_qt's
/// (Zoom/Orient animate). Draw/Ray saves the `ray` image (File ▸ Export Image…).
/// Builder and Properties have no counterpart here, so that row holds Rebuild alone.
struct ClassicButtonGrid: View {
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    VStack(spacing: 4) {
      HStack(spacing: 4) {
        button("Reset") { run("reset") }
        button("Zoom") { run("zoom animate=1") }
        button("Orient") { run("orient animate=1") }
        button("Draw/Ray") {
          NotificationCenter.default.post(name: .mmolExportImage, object: nil)
        }
      }
      HStack(spacing: 4) {
        button("Unpick") { run("unpick") }
        button("Deselect") { run("deselect") }
        button("Rock") { run("rock") }
        button("Get View") { run("get_view") }
      }
      HStack(spacing: 4) {
        button("|<") { run("rewind") }
        button("<") { run("backward") }
        button("Stop") { run("mstop") }
        button("Play") { run("mplay") }
        button(">") { run("forward") }
        button(">|") { run("ending") }
        button("MClear") { run("mclear") }
      }
      button("Rebuild") { run("rebuild") }
      HStack(spacing: 4) {
        button("Undo") { run("undo") }
        button("Redo") { run("redo") }
      }
      Spacer(minLength: 0)
    }
    .padding(6)
    .frame(maxHeight: .infinity, alignment: .top)
    .background(ClassicUI.panelBg)
  }

  private func button(_ title: String, action: @escaping () -> Void) -> some View {
    Button(title, action: action).buttonStyle(ClassicGridButtonStyle())
  }

  private func run(_ cmd: String) { engine.runCommand(cmd) }
}

/// The last `ray` image in place of the live scene, drawn the way desktop
/// PyMOL's SceneDrawImageOverlay draws it (src/pymol/scene/scene.cc): the scene area is
/// background colour; an image smaller than the viewport sits centred at its true
/// pixel size inside a 1px border in the background's contrast colour
/// (ColorGetBkrdContColor); transparent pixels show a 16px checkerboard of the
/// background ± 0x30 (show_alpha_checker); an oversize image is fitted.
/// Never takes clicks — interacting with the viewport changes the view, the core
/// drops the image, and the engine clears this overlay.
struct ClassicRayOverlay: View {
  let ray: PyMOLEngine.ClassicRayImage
  @Environment(\.displayScale) private var scale

  var body: some View {
    GeometryReader { geo in
      let imgPt = CGSize(
        width: CGFloat(ray.pixelWidth) / scale,
        height: CGFloat(ray.pixelHeight) / scale)
      let fitsExactly =
        nearly(imgPt, geo.size)
        || nearly(CGSize(width: ray.pixelWidth, height: ray.pixelHeight), geo.size)
      let oversize =
        !fitsExactly
        && (imgPt.width > geo.size.width || imgPt.height > geo.size.height)
      let shown: CGSize =
        fitsExactly
        ? geo.size
        : oversize ? fitted(imgPt, in: geo.size) : imgPt
      let bordered = !fitsExactly && !oversize

      ZStack {
        ray.background.color
        ZStack {
          Canvas { ctx, size in drawChecker(ctx, size) }
          Image(nsImage: ray.image)
            .resizable()
            .interpolation(fitsExactly || bordered ? .none : .high)
        }
        .frame(width: shown.width, height: shown.height)
        .padding(bordered ? 1 / scale : 0)
        .background(bordered ? contrastColor : .clear)
      }
      .frame(width: geo.size.width, height: geo.size.height)
    }
    .allowsHitTesting(false)
  }

  private func nearly(_ a: CGSize, _ b: CGSize) -> Bool {
    abs(a.width - b.width) <= 2 && abs(a.height - b.height) <= 2
  }

  private func fitted(_ s: CGSize, in box: CGSize) -> CGSize {
    let k = min(box.width / s.width, box.height / s.height)
    return CGSize(width: s.width * k, height: s.height * k)
  }

  /// One background channel shifted by TRN_BKG (0x30) toward mid-grey.
  private func shifted(_ c: Double) -> Double {
    let v = Int(255 * c)
    return Double((v & 0x80) != 0 ? v - 0x30 : v + 0x30) / 255
  }

  private func drawChecker(_ ctx: GraphicsContext, _ size: CGSize) {
    let bg = ray.background
    ctx.fill(Path(CGRect(origin: .zero, size: size)), with: .color(bg.color))
    let alt = Color(
      .sRGB, red: shifted(bg.r), green: shifted(bg.g),
      blue: shifted(bg.b), opacity: 1)
    // 16 image pixels per square, in this view's points.
    let sq = 16 * size.width / CGFloat(max(ray.pixelWidth, 1))
    guard sq > 0.5 else { return }
    var path = Path()
    var row = 0
    var y: CGFloat = 0
    while y < size.height {
      var col = 0
      var x: CGFloat = 0
      while x < size.width {
        if (row + col) & 1 == 1 { path.addRect(CGRect(x: x, y: y, width: sq, height: sq)) }
        x += sq
        col += 1
      }
      y += sq
      row += 1
    }
    ctx.fill(path, with: .color(alt))
  }

  /// ColorGetBkrdContColor: white on a dark background, black on a light one,
  /// then any channel still within 0.5 of the background flipped.
  private var contrastColor: Color {
    let bg = ray.background
    let base: Double = (bg.r + bg.g + bg.b) > 0.5 ? 1 : 0
    func ch(_ b: Double) -> Double { abs(b - base) < 0.5 ? 1 - base : base }
    return Color(.sRGB, red: ch(bg.r), green: ch(bg.g), blue: ch(bg.b), opacity: 1)
  }
}

/// The internal GUI's movie/view strip (src/pymol/ui/control.cc): rewind, back, stop,
/// play, forward, end, S (sequence viewer), ▼ (rock), F (full screen).
struct ClassicTransportBar: View {
  @EnvironmentObject var engine: PyMOLEngine

  var body: some View {
    HStack(spacing: 1) {
      cell(Image(systemName: "backward.end.fill")) { engine.runCommand("rewind") }
      cell(Image(systemName: "arrowtriangle.left.fill")) { engine.runCommand("backward") }
      cell(Image(systemName: "stop.fill")) { engine.runCommand("mstop") }
      cell(Image(systemName: "play.fill")) {
        engine.runCommand(engine.playback.isPlaying ? "mstop" : "mplay")
      }
      cell(Image(systemName: "arrowtriangle.right.fill")) { engine.runCommand("forward") }
      cell(Image(systemName: "forward.end.fill")) { engine.runCommand("ending") }
      cell(Text("S").font(ClassicUI.font)) { engine.sequenceVisible.toggle() }
      cell(Image(systemName: "arrowtriangle.down.fill")) { engine.runCommand("rock") }
      cell(Text("F").font(ClassicUI.font)) { NSApp.keyWindow?.toggleFullScreen(nil) }
    }
    .frame(height: 22)
    .padding(.horizontal, 2)
    .background(ClassicUI.panelBg)
  }

  private func cell<L: View>(_ label: L, action: @escaping () -> Void) -> some View {
    Button(action: action) {
      label
        .font(.system(size: 10, weight: .bold))
        .foregroundColor(Color(.sRGB, red: 1.0, green: 0.25, blue: 0.25, opacity: 1))
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .background(ClassicBevel(fill: ClassicUI.nameEnabled))
        .contentShape(Rectangle())
    }
    .buttonStyle(.plain)
  }
}
