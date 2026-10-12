// HeatmapPanel.swift — the object card's CSV heatmap editor (C ▸ by file…).

import SwiftUI

// MARK: - Heatmap panel (C ▸ by file…)

/// The chip value that selects an object's heatmap panel rather than a layer.
let kHeatmapChip = "__heatmap"

/// Edits one object's CSV heatmap: colours, value range, value column and
/// index column. Every change recolours the object (Heatmap.swift).
struct HeatmapPanel: View {
  let objName: String
  @EnvironmentObject var engine: PyMOLEngine
  // The range thumbs move freely while dragging and recolour on release.
  @State private var liveRange: ClosedRange<Double>?

  var body: some View {
    if let st = engine.heatmaps[objName] {
      VStack(alignment: .leading, spacing: 6) {
        HStack {
          Text(st.csv.fileName).font(.system(size: 10, weight: .semibold))
            .foregroundColor(PanelTheme.textColor).lineLimit(1).truncationMode(.middle)
          Spacer()
          Text("\(st.residueValues.count) residues")
            .font(.system(size: 9)).foregroundColor(PanelTheme.disabledColor)
        }
        row("Low color") {
          DebouncedColorPicker(
            get: { st.settings.low.color },
            apply: { c in engine.updateHeatmap(objName) { $0.low = c.rgba } })
        }
        row("High color") {
          DebouncedColorPicker(
            get: { st.settings.high.color },
            apply: { c in engine.updateHeatmap(objName) { $0.high = c.rgba } })
        }
        row("White midpoint") {
          Toggle(
            "",
            isOn: Binding(
              get: { st.settings.throughWhite },
              set: { v in engine.updateHeatmap(objName) { $0.throughWhite = v } })
          )
          .toggleStyle(ClassicSwitchStyle()).labelsHidden()
        }
        gradient(st)
        let shown = liveRange ?? (st.settings.minValue...st.settings.maxValue)
        HStack {
          Text("min \(fmt(shown.lowerBound))")
          Spacer()
          Text("max \(fmt(shown.upperBound))")
        }
        .font(.system(size: 9, design: .monospaced)).foregroundColor(PanelTheme.textColor)
        HeatmapRangeSlider(
          bounds: st.dataRange, value: shown,
          onChange: { liveRange = $0 },
          onEnd: { r in
            liveRange = nil
            engine.updateHeatmap(objName) {
              $0.minValue = r.lowerBound
              $0.maxValue = r.upperBound
            }
          }
        )
        .frame(height: 18)
        row("Heatmap column") {
          Picker(
            "",
            selection: Binding(
              get: { st.settings.column },
              set: { v in engine.updateHeatmap(objName) { $0.column = v } })
          ) {
            ForEach(st.csv.numericColumns.filter { $0 != st.settings.indexColumn }, id: \.self) {
              Text($0).tag($0)
            }
          }
          .labelsHidden().controlSize(.small).frame(maxWidth: 150)
        }
        row("Index column") {
          Picker(
            "",
            selection: Binding(
              get: { st.settings.indexColumn ?? "" },
              set: { v in engine.updateHeatmap(objName) { $0.indexColumn = v.isEmpty ? nil : v } })
          ) {
            Text("None (row 1 = residue 1)").tag("")
            ForEach(st.csv.headers, id: \.self) { Text($0).tag($0) }
          }
          .labelsHidden().controlSize(.small).frame(maxWidth: 150)
        }
        HStack {
          Button("Reapply") { engine.applyHeatmap(objName) }
            .help("Recolour from the heatmap again (after colouring the object some other way)")
          Spacer()
          Button("Remove") {
            engine.heatmaps[objName] = nil
          }
          .help("Forget this CSV (the colours stay until you recolour)")
        }
        .controlSize(.small)
      }
      .padding(.vertical, 4)
    }
  }

  private func row<C: View>(_ label: String, @ViewBuilder _ content: () -> C) -> some View {
    HStack {
      Text(label).font(.system(size: 10)).foregroundColor(PanelTheme.textColor)
      Spacer()
      content()
    }
  }

  /// The colour scale itself, low → high.
  private func gradient(_ st: HeatmapState) -> some View {
    let lo = st.settings.minValue
    let hi = st.settings.maxValue
    let stops = (0...8).map { k -> Color in st.color(for: lo + (hi - lo) * Double(k) / 8).color }
    return LinearGradient(colors: stops, startPoint: .leading, endPoint: .trailing)
      .frame(height: 8)
      .clipShape(RoundedRectangle(cornerRadius: 2))
  }

  private func fmt(_ v: Double) -> String {
    abs(v) >= 1000 || (abs(v) < 0.01 && v != 0)
      ? String(format: "%.3g", v) : String(format: "%.2f", v)
  }
}

/// One track, two thumbs: the colour scale's min and max within the data range.
private struct HeatmapRangeSlider: View {
  let bounds: ClosedRange<Double>
  let value: ClosedRange<Double>
  let onChange: (ClosedRange<Double>) -> Void
  let onEnd: (ClosedRange<Double>) -> Void
  @State private var dragging: ClosedRange<Double>?

  var body: some View {
    GeometryReader { geo in
      let w = max(1, geo.size.width - 14)
      let span = max(1e-12, bounds.upperBound - bounds.lowerBound)
      let cur = dragging ?? value
      let x = { (v: Double) in 7 + CGFloat((v - bounds.lowerBound) / span) * w }
      let v = { (px: CGFloat) in bounds.lowerBound + Double(min(max(0, (px - 7) / w), 1)) * span }
      ZStack(alignment: .leading) {
        Capsule().fill(PanelTheme.buttonBackground).frame(height: 4)
        Capsule().fill(PanelTheme.accentColor)
          .frame(width: max(0, x(cur.upperBound) - x(cur.lowerBound)), height: 4)
          .offset(x: x(cur.lowerBound))
        thumb.position(x: x(cur.lowerBound), y: geo.size.height / 2)
          .gesture(
            DragGesture().onChanged { g in
              let lo = min(v(g.location.x), cur.upperBound)
              dragging = lo...cur.upperBound
              onChange(dragging!)
            }.onEnded { _ in
              if let d = dragging { onEnd(d) }
              dragging = nil
            })
        thumb.position(x: x(cur.upperBound), y: geo.size.height / 2)
          .gesture(
            DragGesture().onChanged { g in
              let hi = max(v(g.location.x), cur.lowerBound)
              dragging = cur.lowerBound...hi
              onChange(dragging!)
            }.onEnded { _ in
              if let d = dragging { onEnd(d) }
              dragging = nil
            })
      }
    }
  }

  private var thumb: some View {
    Circle().fill(Color.white).frame(width: 13, height: 13)
      .shadow(color: .black.opacity(0.3), radius: 1, y: 0.5)
  }
}
