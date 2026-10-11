// ContactMapView.swift — A ▸ Contact Map…: a residue × residue Cα contact map
// for one object, in a title-bar tab of the main window, with an Å cutoff slider.
//
// The Cα coordinates are read from PyMOL once (and again on Refresh); every
// pairwise distance is computed here, so moving the slider only re-thresholds
// the cached distances — no round-trip to the core per step.

import AppKit
import SwiftUI

/// One Cα the map has a row/column for.
struct ContactMapResidue: Equatable {
  let chain: String
  let resi: String
  let resn: String
  var label: String { (chain.isEmpty ? "" : chain + ":") + resn + " " + resi }
}

@MainActor
final class ContactMapModel: ObservableObject {
  let objectName: String
  @Published private(set) var residues: [ContactMapResidue] = []
  @Published private(set) var error: String?
  /// Row-major n × n Cα–Cα distances (Å).
  private(set) var distances: [Float] = []

  init(objectName: String) { self.objectName = objectName }

  var count: Int { residues.count }

  func distance(_ i: Int, _ j: Int) -> Float { distances[i * count + j] }

  /// Read the object's Cα atoms (polymer, with coordinates, one per residue
  /// — the first alt-loc) in the current state, then fill the distance table.
  func load(engine: PyMOLEngine) {
    let path = (NSTemporaryDirectory() as NSString)
      .appendingPathComponent("mmol_contactmap_\(UUID().uuidString).json")
    defer { try? FileManager.default.removeItem(atPath: path) }
    let obj = objectName.replacingOccurrences(of: "'", with: "")
    engine.runPython(
      """
      import json
      from pymol import cmd as _c
      _sel = "(%s) and polymer and name CA and present and (alt '' or alt A)" % '\(obj)'
      _res = []
      _xyz = []
      try:
          _c.iterate_state(-1, _sel, '_res.append((chain, resi, resn)); _xyz.extend((x, y, z))',
                           space={'_res': _res, '_xyz': _xyz})
          _out = {'res': _res, 'xyz': _xyz}
      except Exception as _e:
          _out = {'error': str(_e)}
      with open(r'\(path)', 'w') as _f:
          json.dump(_out, _f)
      """)
    struct Payload: Decodable {
      let res: [[String]]?
      let xyz: [Float]?
      let error: String?
    }
    guard let data = FileManager.default.contents(atPath: path),
      let p = try? JSONDecoder().decode(Payload.self, from: data)
    else {
      error = "Couldn't read “\(objectName)” from PyMOL."
      return
    }
    if let e = p.error {
      error = e
      return
    }
    let res = (p.res ?? []).compactMap { r in
      r.count == 3 ? ContactMapResidue(chain: r[0], resi: r[1], resn: r[2]) : nil
    }
    let xyz = p.xyz ?? []
    guard !res.isEmpty, xyz.count == res.count * 3 else {
      error = "“\(objectName)” has no protein Cα atoms to map."
      residues = []
      distances = []
      return
    }
    let n = res.count
    var d = [Float](repeating: 0, count: n * n)
    xyz.withUnsafeBufferPointer { c in
      d.withUnsafeMutableBufferPointer { out in
        for i in 0..<n {
          let xi = c[3 * i]
          let yi = c[3 * i + 1]
          let zi = c[3 * i + 2]
          for j in (i + 1)..<max(i + 1, n) {
            let dx = xi - c[3 * j]
            let dy = yi - c[3 * j + 1]
            let dz = zi - c[3 * j + 2]
            let v = (dx * dx + dy * dy + dz * dz).squareRoot()
            out[i * n + j] = v
            out[j * n + i] = v
          }
        }
      }
    }
    error = nil
    distances = d
    residues = res
  }

  /// Residue indices taking part in a contact (≤ cutoff, not a residue with
  /// itself) inside rows × cols — both partners of each such pair.
  func contactResidues(rows: ClosedRange<Int>, cols: ClosedRange<Int>, cutoff: Float) -> [Int] {
    var hit = Set<Int>()
    for i in rows where i < count {
      for j in cols where j < count && j != i && distances[i * count + j] <= cutoff {
        hit.insert(i)
        hit.insert(j)
      }
    }
    return hit.sorted()
  }

  /// The map at `cutoff` Å: one pixel per residue pair, contacts dark, the
  /// rest white. Also the contact
  /// count (unordered pairs, excluding a residue with itself).
  func image(cutoff: Float) -> (NSImage, Int)? {
    let n = count
    guard n > 0 else { return nil }
    var px = [UInt8](repeating: 255, count: n * n * 4)
    var contacts = 0
    // Raw buffers: this runs n² times on every slider step. (Chain
    // boundaries are drawn by the plot as vector lines, crisp at any zoom.)
    distances.withUnsafeBufferPointer { d in
      px.withUnsafeMutableBufferPointer { out in
        for i in 0..<n {
          for j in 0..<n where d[i &* n &+ j] <= cutoff {
            let k = (i &* n &+ j) &* 4
            out[k] = 20
            out[k &+ 1] = 20
            out[k &+ 2] = 30
            if j > i { contacts &+= 1 }
          }
        }
      }
    }
    guard let provider = CGDataProvider(data: Data(px) as CFData),
      let cg = CGImage(
        width: n, height: n, bitsPerComponent: 8, bitsPerPixel: 32,
        bytesPerRow: n * 4, space: CGColorSpaceCreateDeviceRGB(),
        bitmapInfo: CGBitmapInfo(rawValue: CGImageAlphaInfo.noneSkipLast.rawValue),
        provider: provider, decode: nil, shouldInterpolate: false,
        intent: .defaultIntent)
    else { return nil }
    return (NSImage(cgImage: cg, size: NSSize(width: n, height: n)), contacts)
  }
}

struct ContactMapView: View {
  let objectName: String
  @EnvironmentObject var engine: PyMOLEngine
  @StateObject private var model: ContactMapModel
  @State private var cutoff: Double = 5
  @State private var rendered: (image: NSImage, contacts: Int)?
  @State private var hover: (i: Int, j: Int)?
  @State private var resetToken = 0
  @State private var lastSelection: String?

  init(objectName: String) {
    self.objectName = objectName
    _model = StateObject(wrappedValue: ContactMapModel(objectName: objectName))
  }

  var body: some View {
    VStack(alignment: .leading, spacing: 10) {
      HStack {
        Text(objectName).font(.headline)
        Spacer()
        if let r = rendered {
          Text("\(model.count) residues · \(r.contacts) contacts")
            .font(.caption).foregroundStyle(.secondary)
        }
        Button("Reset Zoom") { resetToken += 1 }
          .help("Show the whole map (scroll or pinch to zoom, ⇧-scroll to pan)")
        Button("Refresh") { reload() }
          .help("Re-read the coordinates (after moving or editing the object)")
      }
      HStack(spacing: 10) {
        Text("Cα cutoff")
        ClassicSlider(value: $cutoff, range: 3...20, step: 0.5)
        Text(String(format: "%.1f Å", cutoff))
          .font(.system(.body, design: .monospaced))
          .frame(width: 60, alignment: .trailing)
      }
      if let e = model.error {
        Text(e).foregroundStyle(.secondary)
          .frame(maxWidth: .infinity, maxHeight: .infinity)
      } else if let r = rendered {
        ContactMapPlot(
          image: r.image, residues: model.residues, resetToken: resetToken,
          onHover: { h in hover = h.map { (i: $0.0, j: $0.1) } },
          onBox: { rows, cols, extend, isClick in
            selectBox(rows: rows, cols: cols, extend: extend, isClick: isClick)
          }
        )
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        HStack {
          Text(hoverText)
          Spacer()
          if let note = lastSelection { Text(note) }
        }
        .font(.system(.caption, design: .monospaced))
        .foregroundStyle(.secondary)
      } else {
        ProgressView().frame(maxWidth: .infinity, maxHeight: .infinity)
      }
    }
    .padding(14)
    .frame(minWidth: 420, minHeight: 480)
    .classicControls()
    .onAppear {
      reload()
    }
    .onChange(of: cutoff) { _ in redraw() }
  }

  private var hoverText: String {
    guard let h = hover, h.i < model.count, h.j < model.count else {
      return "Hover to read a pair · scroll to zoom, ⇧-scroll to pan · drag a box to select"
    }
    let a = model.residues[h.i]
    let b = model.residues[h.j]
    return String(format: "%@  ×  %@   %.1f Å", a.label, b.label, model.distance(h.i, h.j))
  }

  /// Drag-box (or click) on the map: every residue in a contact inside the box
  /// becomes 'sele' in the main window (⇧ adds to it), quietly.
  private func selectBox(
    rows: ClosedRange<Int>, cols: ClosedRange<Int>, extend: Bool,
    isClick: Bool
  ) {
    let idx = model.contactResidues(rows: rows, cols: cols, cutoff: Float(cutoff))
    guard !idx.isEmpty else {
      if isClick && !extend {
        // A click on white space clears the selection.
        engine.runPython(
          "from pymol import cmd as _c\n_c.select('sele', 'none', enable=0, quiet=1)\n")
        engine.syncSequenceSelectionNow()
        lastSelection = "Selection cleared."
      } else {
        lastSelection = "No contacts in that box."
      }
      return
    }
    var byChain: [String: [String]] = [:]
    var order: [String] = []
    for i in idx {
      let r = model.residues[i]
      if byChain[r.chain] == nil { order.append(r.chain) }
      byChain[r.chain, default: []].append(r.resi)
    }
    let obj = objectName.replacingOccurrences(of: "'", with: "")
    let parts = order.map { ch -> String in
      let resis = byChain[ch]!.joined(separator: "+")
      return ch.isEmpty ? "(resi \(resis))" : "(chain \(ch) and resi \(resis))"
    }
    let expr = "(\(obj)) and (\(parts.joined(separator: " or ")))"
    let sel = extend ? "(?sele) or (\(expr))" : expr
    engine.runPython(
      "from pymol import cmd as _c\n_c.select('sele', '\(sel)', enable=1, quiet=1)\n")
    engine.syncSequenceSelectionNow()
    lastSelection = "Selected \(idx.count) residue\(idx.count == 1 ? "" : "s") in contacts → sele"
  }

  private func reload() {
    model.load(engine: engine)
    redraw()
  }

  private func redraw() {
    rendered = model.image(cutoff: Float(cutoff)).map { (image: $0.0, contacts: $0.1) }
  }
}

// MARK: - Zoomable plot

/// The map inside a fixed frame: pinch / scroll-wheel zooms about the pointer,
/// two-finger scroll pans, axis numbers re-space with the zoom, and a drag
/// boxes a region whose contacts become 'sele'.
struct ContactMapPlot: NSViewRepresentable {
  let image: NSImage
  let residues: [ContactMapResidue]
  var resetToken: Int
  var onHover: ((Int, Int)?) -> Void
  var onBox: (ClosedRange<Int>, ClosedRange<Int>, _ extend: Bool, _ isClick: Bool) -> Void

  func makeNSView(context: Context) -> ContactMapPlotView { ContactMapPlotView() }

  func updateNSView(_ v: ContactMapPlotView, context: Context) {
    v.onHover = onHover
    v.onBox = onBox
    v.update(image: image, residues: residues, resetToken: resetToken)
  }
}

final class ContactMapPlotView: NSView {
  var onHover: ((Int, Int)?) -> Void = { _ in }
  var onBox: (ClosedRange<Int>, ClosedRange<Int>, Bool, Bool) -> Void = { _, _, _, _ in }

  private var image: NSImage?
  private var residues: [ContactMapResidue] = []
  private var resetToken = -1
  /// Visible window in matrix units (residue indices): columns x0..<x0+span,
  /// rows y0..<y0+span. Square, so one span.
  private var x0: CGFloat = 0, y0: CGFloat = 0, span: CGFloat = 1
  private var dragStart: NSPoint?, dragNow: NSPoint?

  private let leftMargin: CGFloat = 46, bottomMargin: CGFloat = 22
  private let topMargin: CGFloat = 4, rightMargin: CGFloat = 6
  private var n: CGFloat { CGFloat(residues.count) }

  override var isFlipped: Bool { true }
  override var acceptsFirstResponder: Bool { true }

  func update(image: NSImage, residues: [ContactMapResidue], resetToken: Int) {
    let resized = residues.count != self.residues.count
    self.image = image
    self.residues = residues
    if resized || resetToken != self.resetToken {
      self.resetToken = resetToken
      x0 = 0
      y0 = 0
      span = max(1, n)
    }
    needsDisplay = true
  }

  // MARK: Geometry

  /// The square plot area, inside the axis margins.
  private var plot: NSRect {
    let w = bounds.width - leftMargin - rightMargin
    let h = bounds.height - topMargin - bottomMargin
    let side = max(1, min(w, h))
    // Centred across the width when the area is wider than it is tall.
    return NSRect(x: leftMargin + max(0, (w - side) / 2), y: topMargin, width: side, height: side)
  }

  private var cell: CGFloat { plot.width / span }

  /// Matrix coordinate (column, row) under a view point.
  private func matrixPoint(_ p: NSPoint) -> CGPoint {
    CGPoint(x: x0 + (p.x - plot.minX) / cell, y: y0 + (p.y - plot.minY) / cell)
  }

  private func clampWindow() {
    let minSpan = min(n, 8)
    span = min(max(span, max(1, minSpan)), max(1, n))
    x0 = min(max(0, x0), max(0, n - span))
    y0 = min(max(0, y0), max(0, n - span))
  }

  private func zoom(by factor: CGFloat, about p: NSPoint) {
    guard n > 0 else { return }
    let m = matrixPoint(p)
    span /= factor
    clampWindow()
    x0 = m.x - (p.x - plot.minX) / cell
    y0 = m.y - (p.y - plot.minY) / cell
    clampWindow()
    needsDisplay = true
  }

  // MARK: Drawing

  override func draw(_ dirtyRect: NSRect) {
    NSColor.windowBackgroundColor.setFill()
    bounds.fill()
    guard let image, n > 0 else { return }
    let plot = self.plot
    NSGraphicsContext.current?.imageInterpolation = .none
    // NSImage source rects run bottom-up; the view is flipped.
    let src = NSRect(x: x0, y: n - (y0 + span), width: span, height: span)
    image.draw(
      in: plot, from: src, operation: .copy, fraction: 1, respectFlipped: true,
      hints: [.interpolation: NSImageInterpolation.none.rawValue])
    drawChainGrid(in: plot)
    NSColor.secondaryLabelColor.withAlphaComponent(0.6).setStroke()
    NSBezierPath(rect: plot.insetBy(dx: -0.5, dy: -0.5)).stroke()
    drawAxes(in: plot)
    if let a = dragStart, let b = dragNow {
      let r = NSRect(
        x: min(a.x, b.x), y: min(a.y, b.y),
        width: abs(a.x - b.x), height: abs(a.y - b.y)
      ).intersection(plot)
      NSColor.systemBlue.withAlphaComponent(0.15).setFill()
      r.fill()
      NSColor.systemBlue.setStroke()
      NSBezierPath(rect: r).stroke()
    }
  }

  /// Ticks on residues whose number is a multiple of a "nice" step (1, 2, 5,
  /// 10, 20, 25, 50 …), the step chosen so labels sit ~28 px apart at the
  /// current zoom — zooming in shows finer numbering, zooming out coarser.
  private func drawAxes(in plot: NSRect) {
    let target = 28 / cell
    var step = 1
    for s in [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000] {
      step = s
      if CGFloat(s) >= target { break }
    }
    let attrs: [NSAttributedString.Key: Any] = [
      .font: NSFont.monospacedSystemFont(ofSize: 9, weight: .regular),
      .foregroundColor: NSColor.secondaryLabelColor,
    ]
    NSColor.secondaryLabelColor.setStroke()
    let lo = max(0, Int(floor(min(x0, y0))))
    let hi = min(residues.count, Int(ceil(max(x0, y0) + span)))
    var lastX = -CGFloat.infinity
    var lastY = -CGFloat.infinity
    for i in lo..<hi {
      guard let num = Int(residues[i].resi.filter { $0.isNumber || $0 == "-" }),
        num % step == 0
      else { continue }
      let c = CGFloat(i) + 0.5
      if c >= x0, c <= x0 + span {
        let x = plot.minX + (c - x0) * cell
        let s = residues[i].resi as NSString
        let sz = s.size(withAttributes: attrs)
        if x - sz.width / 2 > lastX + 3 {
          tick(from: NSPoint(x: x, y: plot.maxY), to: NSPoint(x: x, y: plot.maxY + 3))
          s.draw(at: NSPoint(x: x - sz.width / 2, y: plot.maxY + 4), withAttributes: attrs)
          lastX = x + sz.width / 2
        }
      }
      if c >= y0, c <= y0 + span {
        let y = plot.minY + (c - y0) * cell
        let label = (residues[i].chain.isEmpty ? "" : residues[i].chain + " ") + residues[i].resi
        let s = label as NSString
        let sz = s.size(withAttributes: attrs)
        if y - sz.height / 2 > lastY + 1 {
          tick(from: NSPoint(x: plot.minX - 3, y: y), to: NSPoint(x: plot.minX, y: y))
          s.draw(
            at: NSPoint(x: plot.minX - 5 - sz.width, y: y - sz.height / 2), withAttributes: attrs)
          lastY = y + sz.height / 2
        }
      }
    }
  }

  /// Grey lines across the whole map where one chain ends and the next
  /// begins, so each chain × chain block reads as a cell of a grid.
  private func drawChainGrid(in plot: NSRect) {
    let path = NSBezierPath()
    for b in 1..<max(1, residues.count) where residues[b].chain != residues[b - 1].chain {
      let c = CGFloat(b)
      if c > x0, c < x0 + span {
        let x = (plot.minX + (c - x0) * cell).rounded() + 0.5
        path.move(to: NSPoint(x: x, y: plot.minY))
        path.line(to: NSPoint(x: x, y: plot.maxY))
      }
      if c > y0, c < y0 + span {
        let y = (plot.minY + (c - y0) * cell).rounded() + 0.5
        path.move(to: NSPoint(x: plot.minX, y: y))
        path.line(to: NSPoint(x: plot.maxX, y: y))
      }
    }
    path.lineWidth = 1
    NSColor(white: 0.62, alpha: 1).setStroke()
    path.stroke()
  }

  private func tick(from a: NSPoint, to b: NSPoint) {
    let p = NSBezierPath()
    p.move(to: a)
    p.line(to: b)
    p.lineWidth = 1
    p.stroke()
  }

  // MARK: Input

  override func updateTrackingAreas() {
    super.updateTrackingAreas()
    trackingAreas.forEach(removeTrackingArea)
    addTrackingArea(
      NSTrackingArea(
        rect: bounds,
        options: [.mouseMoved, .mouseEnteredAndExited, .activeInKeyWindow, .inVisibleRect],
        owner: self))
  }

  override func mouseMoved(with event: NSEvent) { reportHover(event) }
  override func mouseExited(with event: NSEvent) { onHover(nil) }

  private func reportHover(_ event: NSEvent) {
    let p = convert(event.locationInWindow, from: nil)
    guard plot.contains(p) else {
      onHover(nil)
      return
    }
    let m = matrixPoint(p)
    onHover((min(residues.count - 1, Int(m.y)), min(residues.count - 1, Int(m.x))))
  }

  override func magnify(with event: NSEvent) {
    zoom(by: 1 + event.magnification, about: convert(event.locationInWindow, from: nil))
  }

  override func scrollWheel(with event: NSEvent) {
    let p = convert(event.locationInWindow, from: nil)
    if event.modifierFlags.contains(.shift) {
      // ⇧-scroll pans (a mouse wheel's ⇧ turns vertical into horizontal).
      x0 -= event.scrollingDeltaX / cell
      y0 -= event.scrollingDeltaY / cell
      clampWindow()
      needsDisplay = true
    } else {
      // Scroll (wheel or trackpad) zooms about the pointer; sideways pans.
      let d =
        event.hasPreciseScrollingDeltas ? event.scrollingDeltaY / 50 : event.scrollingDeltaY / 5
      if d != 0 { zoom(by: pow(1.15, d), about: p) }
      if event.hasPreciseScrollingDeltas, event.scrollingDeltaX != 0 {
        x0 -= event.scrollingDeltaX / cell
        clampWindow()
        needsDisplay = true
      }
    }
    reportHover(event)
  }

  override func mouseDown(with event: NSEvent) {
    let p = convert(event.locationInWindow, from: nil)
    guard plot.contains(p) else { return }
    dragStart = p
    dragNow = p
  }

  override func mouseDragged(with event: NSEvent) {
    guard dragStart != nil else { return }
    dragNow = convert(event.locationInWindow, from: nil)
    needsDisplay = true
  }

  override func mouseUp(with event: NSEvent) {
    defer {
      dragStart = nil
      dragNow = nil
      needsDisplay = true
    }
    guard let a = dragStart, let b = dragNow, residues.count > 0 else { return }
    let last = residues.count - 1
    let ma = matrixPoint(a)
    let mb = matrixPoint(b)
    let c0 = max(0, min(last, Int(min(ma.x, mb.x))))
    let c1 = max(0, min(last, Int(max(ma.x, mb.x))))
    let r0 = max(0, min(last, Int(min(ma.y, mb.y))))
    let r1 = max(0, min(last, Int(max(ma.y, mb.y))))
    let isClick = hypot(a.x - b.x, a.y - b.y) < 3
    onBox(r0...r1, c0...c1, event.modifierFlags.contains(.shift), isClick)
  }
}
