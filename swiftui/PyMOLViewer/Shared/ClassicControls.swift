// ClassicControls.swift — PyMOL 2.5-style controls for the side panels (Render
// Options, Settings, Movie, Align, Mutagenesis, plot tabs): the button grid's
// square grey buttons, dropdowns that open the A/S/H/L/C classic menu, and a
// bevelled slider, checkbox and text field in the A button's blue.
//
// `.classicControls()` on a panel restyles every Button and Toggle inside it;
// pickers and sliders are ClassicPicker / ClassicSegmented / ClassicSlider.

import SwiftUI

/// The button grid's flat dark push button. `selected` fills it with the accent
/// (a chosen option); `prominent` is a lighter grey for a panel's main action.
struct ClassicButtonStyle: ButtonStyle {
  var selected = false
  var prominent = false
  /// The object card's small size.
  var compact = false
  @Environment(\.isEnabled) private var enabled

  func makeBody(configuration: Configuration) -> some View {
    configuration.label
      .font(.system(size: compact ? 10 : 12))
      .foregroundColor(selected ? .white : Color(white: enabled ? 0.92 : 0.5))
      .lineLimit(1)
      .padding(.horizontal, compact ? 5 : 6)
      .frame(minHeight: compact ? 17 : 22)
      .background(
        selected
          ? ClassicUI.toggleA.opacity(configuration.isPressed ? 0.8 : 1)
          : Color(white: (prominent ? 0.34 : 0.22) + (configuration.isPressed ? 0.1 : 0))
      )
      .overlay(
        RoundedRectangle(cornerRadius: 3)
          .strokeBorder(
            selected ? ClassicUI.toggleLight : Color(white: prominent ? 0.55 : 0.38), lineWidth: 1)
      )
      .clipShape(RoundedRectangle(cornerRadius: 3))
      .contentShape(Rectangle())
  }
}

/// Keeps a control's own look but dims it on mouse-down, so custom controls answer
/// the press at once (an `onTapGesture` only reacts on release).
struct ClassicPressStyle: ButtonStyle {
  func makeBody(configuration: Configuration) -> some View {
    configuration.label.brightness(configuration.isPressed ? -0.12 : 0)
  }
}

/// A square bevelled checkbox, filled with the accent when on.
struct ClassicToggleStyle: ToggleStyle {
  func makeBody(configuration: Configuration) -> some View {
    Button {
      configuration.isOn.toggle()
    } label: {
      HStack(spacing: 6) {
        ClassicBevel(
          fill: configuration.isOn ? ClassicUI.toggleA : Color(white: 0.16),
          light: configuration.isOn ? ClassicUI.toggleLight : ClassicUI.darkEdge,
          dark: configuration.isOn ? ClassicUI.toggleDark : ClassicUI.lightEdge
        )
        .frame(width: 13, height: 13)
        .overlay {
          if configuration.isOn {
            Image(systemName: "checkmark").font(.system(size: 8, weight: .black))
              .foregroundColor(.white)
          }
        }
        configuration.label.font(.system(size: 12))
        Spacer(minLength: 0)
      }
      .contentShape(Rectangle())
    }
    .buttonStyle(ClassicPressStyle())
  }
}

/// An on/off switch: a sunken rectangular track (accent when on) and a raised
/// square grey knob.
struct ClassicSwitchStyle: ToggleStyle {
  func makeBody(configuration: Configuration) -> some View {
    HStack(spacing: 6) {
      configuration.label.font(.system(size: 12))
      Button {
        configuration.isOn.toggle()
      } label: {
        ZStack(alignment: configuration.isOn ? .trailing : .leading) {
          ClassicBevel(
            fill: configuration.isOn ? ClassicUI.toggleA : Color(white: 0.12),
            light: ClassicUI.darkEdge, dark: ClassicUI.lightEdge
          )
          .frame(width: 26, height: 12)
          ClassicBevel(fill: Color(white: 0.62)).frame(width: 12, height: 14)
        }
        .frame(width: 26, height: 14)
        .animation(.easeOut(duration: 0.1), value: configuration.isOn)
        .contentShape(Rectangle())
      }
      .buttonStyle(ClassicPressStyle())
    }
  }
}

extension View {
  /// Classic buttons, checkboxes and accent for everything inside.
  func classicControls() -> some View {
    buttonStyle(ClassicButtonStyle())
      .toggleStyle(ClassicToggleStyle())
      .tint(ClassicUI.toggleA)
  }

  /// A sunken classic text field (use with .textFieldStyle(.plain)).
  func classicField() -> some View {
    textFieldStyle(.plain)
      .font(.system(size: 12, design: .monospaced))
      .padding(.horizontal, 4).frame(minHeight: 20)
      .background(
        ClassicBevel(fill: Color(white: 0.12), light: ClassicUI.darkEdge, dark: ClassicUI.lightEdge)
      )
  }
}

/// The face of a classic dropdown: grey raised box, text, ▼.
struct ClassicDropdownFace: View {
  let text: String
  var compact = false
  var body: some View {
    HStack(spacing: 4) {
      Text(text).font(.system(size: compact ? 10 : 12)).lineLimit(1)
      Spacer(minLength: 2)
      Image(systemName: "arrowtriangle.down.fill").font(.system(size: compact ? 6 : 7))
    }
    .foregroundColor(Color(white: 0.95))
    .padding(.horizontal, compact ? 4 : 6)
    .frame(minWidth: compact ? 54 : 70, minHeight: compact ? 16 : 20)
    .background(ClassicBevel(fill: Color(white: 0.4)))
    .fixedSize()
  }
}

/// A dropdown: a grey bevelled box showing the choice, opening the classic menu.
struct ClassicPicker<T: Hashable>: View {
  var label: String? = nil
  var compact = false
  @Binding var selection: T
  let options: [(T, String)]

  var body: some View {
    HStack(spacing: 6) {
      if let label {
        Text(label).font(.system(size: 12))
        Spacer(minLength: 4)
      }
      ClassicMenuTrigger {
        ClassicDropdownFace(text: options.first { $0.0 == selection }?.1 ?? "–", compact: compact)
      } items: {
        options.map { value, title in
          .action((value == selection ? "✓ " : "   ") + title, destructive: false) {
            selection = value
          }
        }
      }
    }
  }
}

/// A row of classic buttons, the chosen one in the accent (a segmented picker).
struct ClassicSegmented<T: Hashable>: View {
  @Binding var selection: T
  let options: [(T, String)]

  var body: some View {
    HStack(spacing: 3) {
      ForEach(options.indices, id: \.self) { i in
        Button {
          selection = options[i].0
        } label: {
          Text(options[i].1).minimumScaleFactor(0.75).frame(maxWidth: .infinity)
        }
        .buttonStyle(ClassicButtonStyle(selected: options[i].0 == selection))
      }
    }
  }
}

/// A sunken track filled in the accent with a raised square thumb.
/// `onEditingChanged(false)` fires on release, like Slider's.
struct ClassicSlider: View {
  @Binding var value: Double
  let range: ClosedRange<Double>
  var step: Double? = nil
  var onEditingChanged: (Bool) -> Void = { _ in }
  @State private var dragging = false

  var body: some View {
    GeometryReader { g in
      let thumb: CGFloat = 10
      let span = max(range.upperBound - range.lowerBound, 1e-9)
      let t = CGFloat(
        (min(max(value, range.lowerBound), range.upperBound) - range.lowerBound) / span)
      let x = t * (g.size.width - thumb)
      ZStack(alignment: .leading) {
        ClassicBevel(fill: Color(white: 0.12), light: ClassicUI.darkEdge, dark: ClassicUI.lightEdge)
          .frame(height: 6)
        ClassicUI.toggleA.frame(width: x + thumb / 2, height: 4).padding(.leading, 1)
        ClassicBevel(fill: Color(white: dragging ? 0.75 : 0.6))
          .frame(width: thumb, height: 14)
          .offset(x: x)
      }
      .frame(maxHeight: .infinity)
      .contentShape(Rectangle())
      .gesture(
        DragGesture(minimumDistance: 0)
          .onChanged { d in
            if !dragging {
              dragging = true
              onEditingChanged(true)
            }
            var v =
              range.lowerBound + Double((d.location.x - thumb / 2) / (g.size.width - thumb)) * span
            if let step { v = (v / step).rounded() * step }
            value = min(max(v, range.lowerBound), range.upperBound)
          }
          .onEnded { _ in
            dragging = false
            onEditingChanged(false)
          })
    }
    .frame(height: 16)
  }
}

/// A rectangular bevelled color swatch; a click opens the system color panel,
/// whose changes go to the last-clicked swatch.
struct ClassicColorWell: View {
  @Binding var selection: Color
  @StateObject private var target = ColorPanelTarget()

  var body: some View {
    ClassicBevel(fill: selection, light: ClassicUI.darkEdge, dark: ClassicUI.lightEdge)
      .frame(width: 38, height: 16)
      .contentShape(Rectangle())
      .onTapGesture {
        target.set = { selection = $0 }
        let panel = NSColorPanel.shared
        panel.showsAlpha = false
        panel.setTarget(nil)
        panel.setAction(nil)
        panel.color = NSColor(selection)
        panel.setTarget(target)
        panel.setAction(#selector(ColorPanelTarget.changed(_:)))
        panel.orderFront(nil)
      }
  }
}

final class ColorPanelTarget: NSObject, ObservableObject {
  var set: (Color) -> Void = { _ in }
  @objc func changed(_ sender: NSColorPanel) { set(Color(nsColor: sender.color)) }
}

/// A bordered progress bar centered just above the bottom of the 3D view that
/// fills as a `ray` runs or a surface / mesh / dots builds, reading `progress`
/// (0..1). The work blocks the main thread, so a background thread moves the
/// layer directly (Core Animation transactions committed and flushed off-main).
final class RayProgressBar {
  private let track = CALayer(), fill = CALayer()
  private let lock = NSLock()
  private var running = true

  static func show(progress: @escaping () -> Float = PyMOLBridge_RayProgress)
    -> RayProgressBar?
  {
    guard let content = (NSApp.mainWindow ?? NSApp.windows.first { $0.isVisible })?.contentView,
      let viewport = find(PyMOLMTKView.self, in: content)
    else { return nil }
    content.wantsLayer = true
    guard let root = content.layer else { return nil }
    let bar = RayProgressBar()
    // The ray blocks the main thread right after this, so commit the new
    // layers now rather than in the run loop's implicit transaction.
    CATransaction.begin()
    defer {
      CATransaction.commit()
      CATransaction.flush()
    }
    let vp = viewport.convert(viewport.bounds, to: content)
    let w = min(360, vp.width * 0.45)
    let h: CGFloat = 12
    let inset: CGFloat = 2
    let lift: CGFloat = 18
    let y = content.isFlipped ? vp.maxY - lift - h : vp.minY + lift
    bar.track.frame = CGRect(x: vp.midX - w / 2, y: y, width: w, height: h)
    bar.track.backgroundColor = NSColor(white: 0.1, alpha: 0.9).cgColor
    bar.track.borderColor = NSColor(white: 0.85, alpha: 1).cgColor
    bar.track.borderWidth = 1
    bar.track.zPosition = 1000
    bar.fill.frame = CGRect(x: inset, y: inset, width: 0, height: h - 2 * inset)
    bar.fill.backgroundColor = NSColor(ClassicUI.toggleA).cgColor
    bar.track.addSublayer(bar.fill)
    root.addSublayer(bar.track)
    let span = w - 2 * inset
    Thread.detachNewThread {
      while bar.isRunning {
        let p = CGFloat(max(0, min(1, progress())))
        CATransaction.begin()
        CATransaction.setDisableActions(true)
        bar.fill.frame.size.width = span * p
        CATransaction.commit()
        CATransaction.flush()
        Thread.sleep(forTimeInterval: 0.03)
      }
    }
    return bar
  }

  private static func find<T: NSView>(_ type: T.Type, in view: NSView) -> T? {
    if let v = view as? T { return v }
    for s in view.subviews { if let v = find(type, in: s) { return v } }
    return nil
  }

  private var isRunning: Bool {
    lock.lock()
    defer { lock.unlock() }
    return running
  }

  func hide() {
    lock.lock()
    running = false
    lock.unlock()
    track.removeFromSuperlayer()
  }
}
