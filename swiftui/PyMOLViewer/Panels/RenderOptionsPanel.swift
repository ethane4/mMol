// RenderOptionsPanel.swift — Display ▸ Render Options.
//
// A title-bar tab: the viewport shows a small sample protein (1UBQ) on its own
// while this panel, in the right column, edits the `ray` settings. The list of
// options and their current values come from modules/mmol/render_options.py
// (RENDEROPTS: marker); settings are PyMOL's own, so every change carries over
// to the user's scene. Render ray-traces the sample; Defaults returns to the
// values the startup scripts (~/.pymolrc) left.

import SwiftUI

struct RenderOption: Identifiable {
  let group, name, label, kind: String
  let extra: Any?
  var id: String { name }
}

extension PyMOLEngine {
  /// The tab became active / inactive: show the sample alone, or the scene.
  func showRenderSample(_ on: Bool) {
    if on {
      let path = Bundle.main.path(forResource: "1ubq", ofType: "cif") ?? ""
      runPython("from mmol import render_options as _ro\n_ro.begin(\(pyLiteral(path)))")
    } else {
      runPython("from mmol import render_options as _ro\n_ro.end()")
    }
  }

  func setRenderOption(_ name: String, _ value: Any) {
    renderOptionValues[name] = value
    guard
      let data = try? JSONSerialization.data(withJSONObject: [value], options: .fragmentsAllowed),
      let json = String(data: data, encoding: .utf8)
    else { return }
    runPython(
      "from mmol import render_options as _ro\nimport json as _j\n"
        + "_ro.set_option('\(name)', _j.loads(\(pyLiteral(json)))[0])")
  }

  func resetRenderOptions() {
    runPython("from mmol import render_options as _ro\n_ro.reset()")
  }

  /// RENDEROPTS:<path of a JSON file {spec, values}>
  func parseRenderOptions(_ line: String) {
    guard
      let data = FileManager.default.contents(atPath: String(line.dropFirst("RENDEROPTS:".count))),
      let obj = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
      let spec = obj["spec"] as? [[String: Any]]
    else { return }
    renderOptionSpec = spec.compactMap { d in
      guard let g = d["group"] as? String, let n = d["name"] as? String,
        let l = d["label"] as? String, let k = d["kind"] as? String
      else { return nil }
      return RenderOption(group: g, name: n, label: l, kind: k, extra: d["extra"])
    }
    renderOptionValues = obj["values"] as? [String: Any] ?? [:]
  }
}

struct RenderOptionsPanel: View {
  @EnvironmentObject var engine: PyMOLEngine

  private var groups: [String] {
    var seen: [String] = []
    for o in engine.renderOptionSpec where !seen.contains(o.group) { seen.append(o.group) }
    return seen
  }

  var body: some View {
    VStack(spacing: 8) {
      Button {
        engine.runCommand("ray 1000,1000")
      } label: {
        Label("Render", systemImage: "camera.aperture")
          .font(.system(size: 14, weight: .semibold))
          .frame(maxWidth: .infinity, minHeight: 30)
      }
      .buttonStyle(ClassicButtonStyle(prominent: true))
      .help("Ray-trace the sample at 1000 × 1000 with these settings")
      Button {
        engine.resetRenderOptions()
      } label: {
        Text("Return to Defaults").frame(maxWidth: .infinity)
      }
      .help("Back to the values your ~/.pymolrc set at startup")
      Divider()
      ScrollView {
        VStack(alignment: .leading, spacing: 10) {
          ForEach(groups, id: \.self) { g in
            Text(g.uppercased()).font(.system(size: 10, weight: .bold))
              .foregroundColor(.secondary).padding(.top, 4)
            ForEach(engine.renderOptionSpec.filter { $0.group == g }) { row($0) }
          }
        }
        .padding(.horizontal, 2)
      }
    }
    .padding(10)
    .background(ClassicUI.panelBg)
    .classicControls()
  }

  @ViewBuilder private func row(_ o: RenderOption) -> some View {
    let v = engine.renderOptionValues[o.name]
    switch o.kind {
    case "toggle", "rep":
      Toggle(
        o.label,
        isOn: Binding(
          get: { (v as? Bool) ?? ((v as? NSNumber)?.intValue ?? 0 != 0) },
          set: { engine.setRenderOption(o.name, $0) }))
    case "menu", "cartoon":
      let options = (o.extra as? [[Any]]) ?? []
      ClassicPicker(
        label: o.label,
        selection: Binding(
          get: { "\(v ?? "")" },
          set: { s in
            if let opt = options.first(where: { "\($0[0])" == s }) {
              engine.setRenderOption(o.name, opt[0])
            }
          }), options: options.map { ("\($0[0])", "\($0[1])") })
    case "color":
      let c = (v as? [Double]) ?? [0, 0, 0]
      HStack {
        Text(o.label).font(.system(size: 12))
        Spacer()
        ClassicColorWell(
          selection: Binding(
            get: { Color(red: c[0], green: c[1], blue: c[2]) },
            set: { col in
              guard let n = NSColor(col).usingColorSpace(.sRGB) else { return }
              engine.setRenderOption(
                o.name, [n.redComponent, n.greenComponent, n.blueComponent].map(Double.init))
            }))
      }
    case "vector":
      let vec = (v as? [Double]) ?? [0, 0, -1]
      Text(o.label).font(.system(size: 12))
      ForEach(0..<3, id: \.self) { k in
        slider(["x", "y", "z"][k], vec[k], -1...1) { x in
          var n = vec
          n[k] = x
          engine.setRenderOption(o.name, n)
        }
      }
    default:  // slider
      let r = (o.extra as? [Double]) ?? [0, 1]
      slider(o.label, (v as? NSNumber)?.doubleValue ?? 0, r[0]...r[1]) {
        engine.setRenderOption(o.name, $0)
      }
    }
  }

  /// A labelled slider that applies as it moves (throttled) and on release.
  private func slider(
    _ label: String, _ value: Double, _ range: ClosedRange<Double>,
    set: @escaping (Double) -> Void
  ) -> some View {
    RenderSlider(label: label, value: value, range: range, commit: set)
  }
}

private struct RenderSlider: View {
  let label: String
  let value: Double
  let range: ClosedRange<Double>
  let commit: (Double) -> Void
  @State private var drag: Double? = nil
  @State private var lastLive = Date.distantPast

  var body: some View {
    VStack(alignment: .leading, spacing: 1) {
      HStack {
        Text(label).font(.system(size: 12))
        Spacer()
        Text(String(format: "%.3g", drag ?? value)).font(.system(size: 11).monospacedDigit())
          .foregroundColor(.secondary)
      }
      ClassicSlider(
        value: Binding(
          get: { drag ?? value },
          set: { v in
            drag = v
            // ~20 updates/s so shape changes follow the drag.
            if Date().timeIntervalSince(lastLive) > 0.05 {
              lastLive = Date()
              commit(v)
            }
          }),
        range: range,
        onEditingChanged: { editing in
          if !editing, let d = drag {
            commit(d)
            drag = nil
          }
        })
    }
  }
}
