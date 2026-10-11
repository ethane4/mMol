// AnalysisTabs.swift — A ▸ Ramachandran Plot (a plot tab, like the contact
// map), and the Align and Mutagenesis tool tabs (right column, live viewport).
// The Python side is modules/mmol/tools.py.

import AppKit
import SwiftUI

extension PyMOLEngine {
  /// Run mmol_tools.<call>(args…, path) and decode the JSON it writes to path.
  func toolCall<T: Decodable>(_ call: String, _ args: [Any] = []) -> T? {
    let path = (NSTemporaryDirectory() as NSString)
      .appendingPathComponent("mmol_tool_\(UUID().uuidString).json")
    defer { try? FileManager.default.removeItem(atPath: path) }
    guard let data = try? JSONSerialization.data(withJSONObject: args + [path]),
      let json = String(data: data, encoding: .utf8)
    else { return nil }
    runPython(
      "from mmol import tools as _mt\nimport json as _j\n"
        + "_mt.\(call)(*_j.loads(r'''\(json)'''))")
    guard let out = FileManager.default.contents(atPath: path) else { return nil }
    return try? JSONDecoder().decode(T.self, from: out)
  }

  func endMutagenesis() {
    runPython("from mmol import tools as _mt\n_mt.mut_end()")
  }

  /// Molecule objects and public selections, for the pickers.
  var moleculeNames: [String] { objects.filter { !$0.isSelection && !$0.isGroup }.map(\.name) }
  var selectionNames: [String] { objects.filter(\.isSelection).map(\.name) }

  /// Make `expr` the 'sele' selection (empty clears it).
  func selectQuietly(_ expr: String?) {
    let code =
      expr.map { "_c.select('sele', '''\($0)''', enable=1, quiet=1)" }
      ?? "_c.select('sele', 'none', enable=0, quiet=1)"
    runPython("from pymol import cmd as _c\n" + code)
    syncSequenceSelectionNow()
  }
}

// MARK: - Ramachandran plot

struct RamaResidue: Decodable {
  let chain, resi, resn: String
  let phi, psi: Double
  var label: String { (chain.isEmpty ? "" : chain + ":") + resn + " " + resi }
  var color: Color {
    resn == "GLY" ? .orange : resn == "PRO" ? .green : Color(red: 0.25, green: 0.45, blue: 0.95)
  }
}

struct RamachandranView: View {
  let objectName: String
  @EnvironmentObject var engine: PyMOLEngine
  @State private var residues: [RamaResidue] = []
  @State private var error: String?
  @State private var hover: Int?
  @State private var box: (CGPoint, CGPoint)?
  @State private var note: String?

  private struct Payload: Decodable {
    let res: [RamaResidue]?
    let error: String?
  }

  var body: some View {
    VStack(alignment: .leading, spacing: 10) {
      HStack(spacing: 14) {
        Text(objectName).font(.headline)
        legend("General", RamaResidue(chain: "", resi: "", resn: "ALA", phi: 0, psi: 0).color)
        legend("Gly", .orange)
        legend("Pro", .green)
        Spacer()
        Text("\(residues.count) residues").font(.caption).foregroundStyle(.secondary)
        Button("Refresh") { load() }
          .help("Re-read the angles (after editing the object)")
      }
      if let error {
        Text(error).foregroundStyle(.secondary)
          .frame(maxWidth: .infinity, maxHeight: .infinity)
      } else {
        GeometryReader { g in
          let side = min(g.size.width, g.size.height) - 44
          plot(side: max(side, 100))
            .frame(maxWidth: .infinity, maxHeight: .infinity)
        }
        HStack {
          Text(hoverText)
          Spacer()
          if let note { Text(note) }
        }
        .font(.system(.caption, design: .monospaced))
        .foregroundStyle(.secondary)
      }
    }
    .padding(14)
    .classicControls()
    .onAppear(perform: load)
  }

  private func legend(_ s: String, _ c: Color) -> some View {
    HStack(spacing: 4) {
      Circle().fill(c).frame(width: 8, height: 8)
      Text(s).font(.caption)
    }
  }

  private var hoverText: String {
    guard let i = hover, i < residues.count else {
      return
        "Hover a point · click to select its residue · drag a box to select · shaded: typical β / α regions (guide only)"
    }
    let r = residues[i]
    return String(format: "%@   φ %.1f°   ψ %.1f°", r.label, r.phi, r.psi)
  }

  private func point(_ r: RamaResidue, _ side: CGFloat) -> CGPoint {
    CGPoint(x: (r.phi + 180) / 360 * side, y: (180 - r.psi) / 360 * side)
  }

  private func plot(side: CGFloat) -> some View {
    Canvas { ctx, size in
      let s = side
      ctx.fill(Path(CGRect(x: 0, y: 0, width: s, height: s)), with: .color(Color(white: 0.97)))
      // Typical regions: β, right-handed α, left-handed α (a guide, not a contour).
      for (phi, psi, rx, ry) in [
        (-120.0, 135.0, 60.0, 45.0), (-65.0, -40.0, 35.0, 30.0), (60.0, 42.0, 20.0, 25.0),
      ] {
        let c = CGPoint(x: (phi + 180) / 360 * s, y: (180 - psi) / 360 * s)
        let rect = CGRect(
          x: c.x - rx / 360 * s, y: c.y - ry / 360 * s,
          width: 2 * rx / 360 * s, height: 2 * ry / 360 * s)
        ctx.fill(Path(ellipseIn: rect), with: .color(Color(red: 0.75, green: 0.85, blue: 1.0)))
      }
      var grid = Path()
      for k in 0...6 {
        let v = CGFloat(k) / 6 * s
        grid.move(to: CGPoint(x: v, y: 0))
        grid.addLine(to: CGPoint(x: v, y: s))
        grid.move(to: CGPoint(x: 0, y: v))
        grid.addLine(to: CGPoint(x: s, y: v))
      }
      ctx.stroke(grid, with: .color(.gray.opacity(0.35)), lineWidth: 0.5)
      var axes = Path()
      axes.move(to: CGPoint(x: s / 2, y: 0))
      axes.addLine(to: CGPoint(x: s / 2, y: s))
      axes.move(to: CGPoint(x: 0, y: s / 2))
      axes.addLine(to: CGPoint(x: s, y: s / 2))
      ctx.stroke(axes, with: .color(.gray.opacity(0.7)), lineWidth: 0.8)
      for (i, r) in residues.enumerated() {
        let p = point(r, s)
        let d: CGFloat = i == hover ? 7 : 4.5
        ctx.fill(
          Path(ellipseIn: CGRect(x: p.x - d / 2, y: p.y - d / 2, width: d, height: d)),
          with: .color(r.color))
      }
      if let (a, b) = box {
        let rect = CGRect(
          x: min(a.x, b.x), y: min(a.y, b.y), width: abs(a.x - b.x), height: abs(a.y - b.y))
        ctx.fill(Path(rect), with: .color(.accentColor.opacity(0.12)))
        ctx.stroke(Path(rect), with: .color(.accentColor), lineWidth: 1)
      }
      ctx.stroke(Path(CGRect(x: 0, y: 0, width: s, height: s)), with: .color(.gray), lineWidth: 1)
    }
    .frame(width: side, height: side)
    .overlay(alignment: .bottom) {  // φ ticks
      HStack {
        ForEach([-180, -120, -60, 0, 60, 120, 180], id: \.self) {
          Text("\($0)")
          if $0 < 180 { Spacer() }
        }
      }
      .font(.system(size: 10)).foregroundStyle(.secondary).offset(y: 16)
    }
    .overlay(alignment: .leading) {  // ψ ticks
      VStack(alignment: .trailing) {
        ForEach([180, 120, 60, 0, -60, -120, -180], id: \.self) {
          Text("\($0)")
          if $0 > -180 { Spacer() }
        }
      }
      .font(.system(size: 10)).foregroundStyle(.secondary).offset(x: -30)
    }
    .overlay(alignment: .bottom) { Text("φ").font(.callout).offset(y: 32) }
    .overlay(alignment: .leading) { Text("ψ").font(.callout).offset(x: -46) }
    .contentShape(Rectangle())
    .onContinuousHover { phase in
      if case .active(let p) = phase { hover = nearest(p, side) } else { hover = nil }
    }
    .gesture(
      DragGesture(minimumDistance: 0)
        .onChanged { v in box = (v.startLocation, v.location) }
        .onEnded { v in
          box = nil
          let extend = NSEvent.modifierFlags.contains(.shift)
          if hypot(v.translation.width, v.translation.height) < 3 {
            select(nearest(v.location, side).map { [$0] } ?? [], extend: extend, click: true)
          } else {
            let rect = CGRect(
              x: min(v.startLocation.x, v.location.x), y: min(v.startLocation.y, v.location.y),
              width: abs(v.translation.width), height: abs(v.translation.height))
            select(
              residues.indices.filter { rect.contains(point(residues[$0], side)) }, extend: extend,
              click: false)
          }
        })
  }

  private func nearest(_ p: CGPoint, _ side: CGFloat) -> Int? {
    var best: (Int, CGFloat)?
    for (i, r) in residues.enumerated() {
      let q = point(r, side)
      let d = hypot(q.x - p.x, q.y - p.y)
      if d < 7, d < (best?.1 ?? .infinity) { best = (i, d) }
    }
    return best?.0
  }

  /// The residues become 'sele' in the main view (⇧ adds to it).
  private func select(_ idx: [Int], extend: Bool, click: Bool) {
    guard !idx.isEmpty else {
      if click && !extend {
        engine.selectQuietly(nil)
        note = "Selection cleared."
      }
      return
    }
    let obj = objectName.replacingOccurrences(of: "'", with: "")
    let parts = Dictionary(grouping: idx.map { residues[$0] }, by: \.chain).map { ch, rs in
      let resis = rs.map(\.resi).joined(separator: "+")
      return ch.isEmpty ? "(resi \(resis))" : "(chain \(ch) and resi \(resis))"
    }
    let expr = "(\(obj)) and (\(parts.joined(separator: " or ")))"
    engine.selectQuietly(extend ? "(?sele) or (\(expr))" : expr)
    note = "Selected \(idx.count) residue\(idx.count == 1 ? "" : "s") → sele"
  }

  private func load() {
    let p: Payload? = engine.toolCall("rama", [objectName])
    error = p == nil ? "Couldn't read “\(objectName)”." : p?.error
    residues = p?.res ?? []
    if error == nil && residues.isEmpty {
      error = "“\(objectName)” has no protein residues with φ/ψ."
    }
  }
}

// MARK: - Align

struct AlignPanel: View {
  @EnvironmentObject var engine: PyMOLEngine
  @State private var mobile = ""
  @State private var target = ""
  @State private var mobileRegion = "all"
  @State private var targetRegion = "all"
  @State private var method = "align"
  @State private var caOnly = true
  @State private var result: AlignResult?

  struct AlignResult: Decodable {
    let rmsd, rmsd_before, total_rmsd: Double?
    let atoms, atoms_before, total_atoms: Int?
    let error: String?
  }

  var body: some View {
    ScrollView {
      VStack(alignment: .leading, spacing: 14) {
        Text("Align Structures").font(.headline)
        side("Move (mobile)", $mobile, $mobileRegion)
        side("Onto (target)", $target, $targetRegion)
        VStack(alignment: .leading, spacing: 6) {
          Text("Method").font(.system(size: 12, weight: .semibold))
          ClassicSegmented(
            selection: $method,
            options: [("align", "align"), ("super", "super"), ("cealign", "cealign")])
          Text(methodBlurb).font(.caption).foregroundStyle(.secondary)
            .fixedSize(horizontal: false, vertical: true)
          Toggle("Fit on Cα atoms only", isOn: $caOnly)
            .disabled(method == "cealign")
        }
        Button(action: run) {
          Label("Align", systemImage: "arrow.triangle.merge")
            .font(.system(size: 14, weight: .semibold))
            .frame(maxWidth: .infinity, minHeight: 30)
        }
        .buttonStyle(ClassicButtonStyle(prominent: true))
        .disabled(mobile.isEmpty || target.isEmpty || mobile == target)
        if let r = result { resultBox(r) }
      }
      .padding(12)
    }
    .background(ClassicUI.panelBg)
    .classicControls()
    .onChange(of: engine.alignMobile) { _ in fillDefaults() }
    .onChange(of: engine.moleculeNames) { _ in fillDefaults() }
    .onAppear(perform: fillDefaults)
  }

  /// Mobile: the object the wizard was opened from, else the second object;
  /// target: the first other object. Re-run as objects come and go.
  private func fillDefaults() {
    let names = engine.moleculeNames
    if let m = engine.alignMobile {
      mobile = m
      engine.alignMobile = nil
    }
    if !names.contains(mobile) { mobile = names.count > 1 ? names[1] : names.first ?? "" }
    if !names.contains(target) || target == mobile { target = names.first { $0 != mobile } ?? "" }
  }

  private var methodBlurb: String {
    switch method {
    case "super":
      return "Structure-based superposition; ignores sequence. Good for low sequence identity."
    case "cealign":
      return "Combinatorial extension on Cα; best for distant homologs. Uses the regions' Cα atoms."
    default:
      return
        "Sequence alignment, then superposition with outlier rejection. Best for close homologs."
    }
  }

  private func side(_ title: String, _ obj: Binding<String>, _ region: Binding<String>) -> some View
  {
    VStack(alignment: .leading, spacing: 6) {
      Text(title).font(.system(size: 12, weight: .semibold))
      ClassicPicker(label: "Object", selection: obj, options: engine.moleculeNames.map { ($0, $0) })
      HStack(spacing: 4) {
        TextField("Region (e.g. resi 10-80, chain A)", text: region)
          .classicField()
        ClassicMenuTrigger {
          Image(systemName: "arrowtriangle.down.fill").font(.system(size: 7))
            .foregroundColor(Color(white: 0.95))
            .frame(width: 20, height: 20)
            .background(ClassicBevel(fill: Color(white: 0.4)))
        } items: {
          (["all"] + engine.selectionNames).map { s in
            .action(s, destructive: false) { region.wrappedValue = s }
          }
        }
        .help("Use a selection as the region")
      }
      Text("Only this region of the object is used for the fit.")
        .font(.caption).foregroundStyle(.secondary)
    }
  }

  private func resultBox(_ r: AlignResult) -> some View {
    VStack(alignment: .leading, spacing: 6) {
      if let e = r.error {
        Text(e).foregroundStyle(.red).font(.caption)
      } else {
        row("Aligned region RMSD", r.rmsd, r.atoms, "atoms")
        if let b = r.rmsd_before, let n = r.atoms_before, n != r.atoms {
          row("  before outlier rejection", b, n, "atoms")
        }
        row("Whole-object Cα RMSD", r.total_rmsd, r.total_atoms, "pairs")
        Text("Whole-object: every sequence-matched Cα pair, measured after the fit.")
          .font(.caption).foregroundStyle(.secondary)
          .fixedSize(horizontal: false, vertical: true)
      }
    }
    .padding(10)
    .frame(maxWidth: .infinity, alignment: .leading)
    .background(RoundedRectangle(cornerRadius: 6).fill(Color.primary.opacity(0.06)))
  }

  private func row(_ label: String, _ v: Double?, _ n: Int?, _ unit: String) -> some View {
    HStack {
      Text(label).font(.system(size: 12))
      Spacer()
      Text(v.map { String(format: "%.3f Å", $0) } ?? "–").font(
        .system(size: 12, weight: .semibold).monospacedDigit())
      Text(n.map { "\($0) \(unit)" } ?? "").font(.caption).foregroundStyle(.secondary)
        .frame(width: 74, alignment: .trailing)
    }
  }

  private func run() {
    result =
      engine.toolCall("align", [mobile, mobileRegion, target, targetRegion, method, caOnly])
      ?? AlignResult(
        rmsd: nil, rmsd_before: nil, total_rmsd: nil, atoms: nil, atoms_before: nil,
        total_atoms: nil, error: "The alignment didn't run.")
    engine.requestViewportRedraw()
  }
}

// MARK: - Mutagenesis

struct MutagenesisPanel: View {
  @EnvironmentObject var engine: PyMOLEngine
  @State private var status: MutStatus?

  struct MutStatus: Decodable {
    let error, residue, mode: String?
    let states, state: Int?
    let scores: [Double]?
  }

  private static let aminoAcids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  ]

  var body: some View {
    ScrollView {
      VStack(alignment: .leading, spacing: 14) {
        Text("Mutagenesis").font(.headline)
        step("1", "Click a residue in the viewer (it becomes sele), then:")
        Button {
          status = engine.toolCall("mut_select", ["sele"])
        } label: {
          Text("Mutate Selected Residue").frame(maxWidth: .infinity)
        }
        if let e = status?.error {
          Text(e).font(.caption).foregroundStyle(.red).fixedSize(horizontal: false, vertical: true)
        }
        if let res = status?.residue {
          Text(res).font(.system(size: 12, design: .monospaced))
          step("2", "Mutate to:")
          LazyVGrid(
            columns: Array(repeating: GridItem(.flexible(), spacing: 4), count: 4), spacing: 4
          ) {
            ForEach(Self.aminoAcids, id: \.self) { aa in
              Button {
                status = engine.toolCall("mut_mode", [aa])
              } label: {
                Text(aa).frame(maxWidth: .infinity)
              }
              .buttonStyle(ClassicButtonStyle(selected: status?.mode == aa))
            }
          }
          Button {
            status = engine.toolCall("mut_mode", ["current"])
          } label: {
            Text("No change (current residue)").frame(maxWidth: .infinity)
          }
          .buttonStyle(ClassicButtonStyle(selected: status?.mode == "current"))
          rotamers
          Button(action: { status = engine.toolCall("mut_apply") }) {
            Label("Apply", systemImage: "checkmark")
              .font(.system(size: 14, weight: .semibold))
              .frame(maxWidth: .infinity, minHeight: 30)
          }
          .buttonStyle(ClassicButtonStyle(prominent: true))
        }
      }
      .padding(12)
    }
    .background(ClassicUI.panelBg)
    .classicControls()
  }

  private func step(_ n: String, _ text: String) -> some View {
    HStack(alignment: .firstTextBaseline, spacing: 6) {
      Text(n).font(.system(size: 12, weight: .bold)).foregroundColor(ClassicUI.accent)
      Text(text).font(.system(size: 12)).fixedSize(horizontal: false, vertical: true)
    }
  }

  @ViewBuilder private var rotamers: some View {
    let n = status?.states ?? 0
    if n > 0 {
      step("3", "Rotamer (lower strain = fewer clashes):")
      VStack(spacing: 2) {
        ForEach(1...n, id: \.self) { i in
          let score = status?.scores.flatMap { i <= $0.count ? $0[i - 1] : nil }
          HStack {
            Text("\(i)").font(.system(size: 12).monospacedDigit())
            Spacer()
            if let score {
              Text(String(format: "strain %.1f", score)).font(.caption.monospacedDigit())
            }
          }
          .padding(.horizontal, 8).padding(.vertical, 3)
          .background(
            RoundedRectangle(cornerRadius: 4)
              .fill(
                status?.state == i ? ClassicUI.accent.opacity(0.35) : Color.primary.opacity(0.05))
          )
          .contentShape(Rectangle())
          .onTapGesture { status = engine.toolCall("mut_state", [i]) }
        }
      }
    }
  }
}
