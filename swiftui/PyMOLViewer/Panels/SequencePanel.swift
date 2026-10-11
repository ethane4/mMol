// SequencePanel.swift — Horizontal scrollable sequence viewer
// SwiftUI reimplementation of PyMOL's seq_view (src/pymol/ui/seeker.cc + src/pymol/ui/seq.cc).
//
// Features mirrored from the original:
//  - One-letter codes per residue, colored by the residue's ACTUAL molecular
//    color (guide-atom color) so spectrum/by-element/custom coloring shows.
//  - Selection highlight synced both ways with the active 'sele' selection.
//  - Click to select a residue; drag to select a range; Shift extends,
//    Ctrl selects-and-centers; click on empty space deselects.
//  - Residue-number ruler above each sequence (every N residues).

import AppKit
import SwiftUI

// MARK: - Data Models

/// A single residue in the sequence display.
struct SequenceResidue: Identifiable, Equatable {
  let id: String  // unique key: "obj/chain/resi/index"
  let objectName: String
  let chain: String
  let oneLetter: String
  let resi: String  // residue index (e.g. "42")
  let resn: String  // three-letter code (e.g. "ALA")
  let color: Color  // real guide-atom color
  var isGap: Bool = false  // alignment gap placeholder (not a real residue)
  var isBreak: Bool = false  // chain-boundary separator bar (not a real residue)
  var isChainLabel: Bool = false  // chain-ID indicator at a chain run's start
  var chainLabel: String = ""  // the chain ID shown when isChainLabel
  var isSpacer: Bool = false  // blank cell between entities (polymer, each ligand…)
  var isHet: Bool = false  // a column of a ligand / ion / water's 3-letter name

  /// Real, selectable residue (excludes gap / break / chain-label / spacer placeholders).
  var isSelectable: Bool { !isGap && !isBreak && !isChainLabel && !isSpacer }

  /// Identity used to match against the active selection set (obj/chain/resi).
  var selKey: String { "\(objectName)/\(chain)/\(resi)" }
}

/// A group of residues belonging to one molecular object.
struct SequenceObject: Identifiable, Equatable {
  let id: String  // object name
  let name: String
  let residues: [SequenceResidue]
}

// MARK: - Three-letter to one-letter mapping

private let aa3to1: [String: String] = [
  "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
  "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
  "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
  "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
  "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
  "DA": "A", "DC": "C", "DG": "G", "DT": "T",
  "A": "A", "C": "C", "G": "G", "U": "U",
]

// MARK: - Look

// Desktop PyMOL's seq_view look (src/pymol/ui/seq.cc) — the scene's black
// background, labels and ruler in the foreground colour (seq_view_label_color
// "front"), residues in their atom colours, 8 x 13 pt character cells, a
// residue number every 5 (seq_view_label_spacing).
/// One object's ruler + residue rows, in points.
let classicSeqRowsHeight: CGFloat = 26
/// Clearance under the rows for the horizontal overlay scroller.
let classicSeqScrollerRoom: CGFloat = 15  // the scroller sits a few points clear of the letters

private let headerColor = Color.white
// Selected-residue highlight uses the selection color (matching the 3D
// selection indicator); text flips to black/white for contrast on that color.
private let selectionBG = ClassicUI.selectionRGBA.color
private var selectionFG: Color {
  let s = ClassicUI.selectionRGBA
  let lum: Double = 0.299 * s.r + 0.587 * s.g + 0.114 * s.b
  return lum > 0.6 ? Color.black : Color.white
}
private let rulerColor = Color.white
private let panelBackground = Color.black
private let cellWidth: CGFloat = 8
private let rulerSpacing = 5  // draw a residue number every N residues
private let rowHeight: CGFloat = 13

// MARK: - SequencePanel View

struct SequencePanel: View {
  @EnvironmentObject var engine: PyMOLEngine
  @State private var lastRange: ClosedRange<Int>? = nil
  @State private var anchorIndex: Int? = nil  // for Shift-click range
  // Optimistic live-highlight during a drag: the selKeys currently under the
  // sweep, unioned into the per-cell highlight test so they light up
  // immediately, without waiting for the async PyMOL selection round-trip
  // (issue #429). Cleared/reconciled against selectedResidueKeys when the drag
  // ends. This is purely visual — the real 'sele' still comes from applyToggle.
  @State private var dragSelKeys: Set<String> = []
  // Desktop PyMOL's drag rule: the residue a drag STARTS on decides whether
  // it selects (start was unselected) or deselects (start was selected). nil
  // between drags. The span is applied against a snapshot of 'sele' taken at
  // the start (`_seqdrag0`), so moving back over residues un-sweeps them.
  @State private var dragAdds: Bool? = nil

  /// Flattened residue list (object order) for range/index mapping.
  private var flat: [SequenceResidue] { engine.sequences.flatMap { $0.residues } }

  /// Fixed gutter width for object-name labels, so every object's residues
  /// start at the same x — required for alignment columns (and aligned rulers)
  /// to line up across the stacked object rows.
  private var labelWidth: CGFloat {
    let longest = engine.sequences.map { $0.name.count }.max() ?? 0
    // Monospaced 12 pt (~7.3 pt a char) and no minimum gutter.
    return min(CGFloat(longest) * 7.3 + 8, 180)
  }

  /// Map residue id -> flat index (for click handling).
  private var idToIndex: [String: Int] {
    var m: [String: Int] = [:]
    for (i, r) in flat.enumerated() { m[r.id] = i }
    return m
  }

  var body: some View {
    // Scroll BOTH axes: horizontally through long sequences, vertically when
    // several objects are loaded (otherwise extra object rows clip against a
    // fixed-height container).
    GeometryReader { geo in
      ScrollView([.horizontal, .vertical], showsIndicators: true) {
        content
          // Pin to the top-left like desktop PyMOL's strip; a two-axis
          // ScrollView otherwise centres content narrower than it. (Less
          // the padding applied below, so nothing scrolls needlessly.)
          .frame(
            minWidth: max(0, geo.size.width - 4),
            minHeight: max(0, geo.size.height - classicSeqScrollerRoom),
            alignment: .topLeading
          )
          .padding(.horizontal, 2)
          // Room for the overlay scroller, which otherwise draws across the
          // bottom row of letters while scrolling.
          .padding(.bottom, classicSeqScrollerRoom)
      }
    }
    .frame(maxWidth: .infinity)
    .background(panelBackground)
    .onAppear { engine.fetchSequences() }
    .onChange(of: engine.objects) { _ in engine.fetchSequences() }
  }

  /// The strip redraws only when its own data or the drag changes: the panel
  /// observes the whole engine, whose other updates (polls, the object card,
  /// playback) would otherwise rebuild every row's ruler and selection mask
  /// — tens of thousands of residues with many objects loaded.
  private struct StripContent: View, Equatable {
    let key: [Int]
    let build: () -> AnyView
    var body: some View { build() }
    static func == (a: Self, b: Self) -> Bool { a.key == b.key }
  }

  private var content: some View {
    StripContent(
      key: [engine.sequenceVersion, dragSelKeys.hashValue, dragAdds.hashValue],
      build: { AnyView(strip) }
    )
    .equatable()
  }

  @ViewBuilder
  private var strip: some View {
    if engine.sequences.isEmpty {
      // Desktop PyMOL's strip is simply empty until something is loaded.
      Color.clear.frame(height: 1)
    } else {
      VStack(alignment: .leading, spacing: 0) {
        let offsets = flatOffsets
        ForEach(engine.sequences) { obj in classicBlock(obj, offset: offsets[obj.id] ?? 0) }
      }
    }
  }

  // MARK: - Strip block

  /// Index of each object's first residue in `flat`.
  private var flatOffsets: [String: Int] {
    var m: [String: Int] = [:]
    var n = 0
    for o in engine.sequences {
      m[o.id] = n
      n += o.residues.count
    }
    return m
  }

  private func classicBlock(_ obj: SequenceObject, offset: Int) -> some View {
    let selected = obj.residues.map { r -> Bool in
      guard r.isSelectable else { return false }
      let inSele = engine.selectedResidueKeys.contains(r.selKey)
      let swept = dragSelKeys.contains(r.selKey)
      return dragAdds == false ? (inSele && !swept) : (inSele || swept)
    }
    return ClassicSequenceRows(
      name: obj.name, residues: obj.residues,
      ruler: rulerCells(obj.residues).map { ClassicRulerCell(text: $0.text, isChain: $0.isChain) },
      selected: selected, labelWidth: labelWidth, cellWidth: cellWidth, rowHeight: rowHeight,
      labelColor: headerColor, rulerColor: rulerColor,
      selectionBG: selectionBG, selectionFG: selectionFG
    )
    .equatable()
    .contentShape(Rectangle())
    .gesture(classicGesture(obj, offset: offset))
  }

  /// Residue slot under local x (nil over the name gutter or past the end).
  private func slot(_ obj: SequenceObject, _ x: CGFloat, clamp: Bool = false) -> Int? {
    let i = Int(floor((x - labelWidth) / cellWidth))
    if clamp { return max(0, min(i, obj.residues.count - 1)) }
    return (x >= labelWidth && i < obj.residues.count) ? i : nil
  }

  /// Desktop PyMOL's strip: a click toggles the residue (or, on the name, adds
  /// the object); a drag sweeps from the residue it started on.
  private func classicGesture(_ obj: SequenceObject, offset: Int) -> some Gesture {
    DragGesture(minimumDistance: 0)
      .onChanged { v in
        guard let start = slot(obj, v.startLocation.x) else { return }
        let moved = abs(v.translation.width) > 3
        guard moved || dragAdds != nil else { return }  // still a click
        let cur = slot(obj, v.location.x, clamp: true) ?? start
        dragStep(from: offset + start, to: offset + cur)
      }
      .onEnded { v in
        if dragAdds != nil {
          dragEnd()
          return
        }
        if v.startLocation.x < labelWidth {
          anchorIndex = nil
          applyToggle(residues: obj.residues, add: true)
        } else if let i = slot(obj, v.startLocation.x), obj.residues[i].isSelectable {
          handleClick(on: obj.residues[i])
        }
      }
  }

  /// Click selection (PyMOL-style, additive/toggle):
  ///  - plain click toggles the residue in/out of 'sele' (add if not
  ///    selected, remove if selected), setting the range anchor;
  ///  - Shift-click adds the range from the anchor to here;
  ///  - Ctrl-click toggles + centers.
  /// Modifiers are read from the current NSEvent on macOS.
  private func handleClick(on r: SequenceResidue) {
    guard let idx = idToIndex[r.id] else { return }
    var shift = false
    var ctrl = false
    let mods = NSEvent.modifierFlags
    shift = mods.contains(.shift)
    ctrl = mods.contains(.control)
    if shift, let a = anchorIndex, !flat.isEmpty {
      // `idx` is resolved against the current `flat`, but `a` was captured on
      // an earlier click; if `flat` shrank since then (an object was removed/
      // replaced), `a` can be >= flat.count. Clamp both ends so the range
      // subscript can never trap out of bounds.
      let n = flat.count
      let lo = max(0, min(min(a, idx), n - 1))
      let hi = max(0, min(max(a, idx), n - 1))
      if lo <= hi {
        applyToggle(residues: Array(flat[lo...hi]), add: true)
      }
    } else {
      anchorIndex = idx
      let isSelected = engine.selectedResidueKeys.contains(r.selKey)
      applyToggle(residues: [r], add: !isSelected, center: ctrl)
    }
  }

  // MARK: - Residue-number ruler

  /// One cell per residue slot, aligned with the residue row below (both use
  /// `cellWidth` slots). Cells carry either a residue-number digit (every
  /// `rulerSpacing`, overflowing into following slots) or — at the start of a
  /// chain run — the chain-ID indicator (`isChain`), which takes priority.
  private func rulerCells(_ residues: [SequenceResidue]) -> [(text: String, isChain: Bool)] {
    var cells = Array(repeating: (text: " ", isChain: false), count: residues.count)
    // Residue numbers (real residues only; numbering restarts per chain).
    // Polymer residues only: ligand and water numbers ran together over
    // the HOH HOH … run ("8880").
    for (i, r) in residues.enumerated() where r.isSelectable && !r.isHet {
      guard let n = Int(r.resi), n % rulerSpacing == 0 else { continue }
      for (j, c) in Array(r.resi).enumerated() where i + j < cells.count {
        cells[i + j] = (text: String(c), isChain: false)
      }
    }
    // Chain IDs sit above each chain run's first residue and win their cell.
    for (i, r) in residues.enumerated() where r.isChainLabel {
      cells[i] = (text: r.chainLabel, isChain: true)
    }
    return cells
  }

  // MARK: - Drag selection

  /// One step of a sweep from flat index `startIdx` to `curIdx`.
  private func dragStep(from startIdx: Int, to curIdx: Int) {
    guard flat.indices.contains(startIdx), flat.indices.contains(curIdx) else { return }
    if dragAdds == nil {
      dragAdds = !engine.selectedResidueKeys.contains(flat[startIdx].selKey)
      engine.runPython(
        "from pymol import cmd as _c\n"
          + "_c.select('_seqdrag0', '?sele', enable=0, quiet=1)\n")
    }
    let range = min(startIdx, curIdx)...max(startIdx, curIdx)
    // Optimistic live highlight: light up the swept residues right now,
    // even before 'sele' is read back. Updated every time the pointer
    // enters a new cell so the visible span tracks the cursor. (issue #429)
    dragSelKeys = Set(flat[range].filter { $0.isSelectable }.map { $0.selKey })
    guard range != lastRange else { return }
    lastRange = range
    applyDragSpan(Array(flat[range]), add: dragAdds ?? true)
    anchorIndex = startIdx
  }

  private func dragEnd() {
    lastRange = nil
    dragAdds = nil
    engine.runPython("from pymol import cmd as _c\n_c.delete('_seqdrag0')\n")
    // Read the real 'sele' back at once, then drop the optimistic overlay
    // on the next turn, after the confirmed keys have landed.
    engine.syncSequenceSelectionNow()
    DispatchQueue.main.async { self.dragSelKeys = [] }
  }

  // MARK: - Selection dispatch

  /// Toggle residues into/out of the active 'sele', mirroring PyMOL's
  /// SeekerSelectionToggle: `add` => `(?sele) or (expr)` (additive),
  /// otherwise `(sele) and not (expr)` (remove). Never replaces — clicks
  /// accumulate, and clicking a selected residue removes just that residue.
  private func applyToggle(residues: [SequenceResidue], add: Bool, center: Bool = false) {
    // Gap/break placeholders carry no real residue — drop them before building
    // a selection expression (an empty resi would produce malformed PyMOL).
    let residues = residues.filter { $0.isSelectable }
    guard !residues.isEmpty else { return }
    let expr = selectionExpression(residues)
    // Quiet Python rather than a console command: desktop PyMOL's sequence
    // viewer doesn't echo `select …` / "Selector: selection defined" on
    // every click.
    var py = Self.modeExpression(expr)
    py +=
      add
      ? "_c.select('sele', '(?sele) or (' + _e + ')', enable=1, quiet=1)\n"
      : "_c.select('sele', '(sele) and not (' + _e + ')', enable=1, quiet=1)\n"
    if center { py += "_c.center('sele', animate=-1)\n" }
    engine.runPython(py)
    // The highlight follows at once: the select above ran synchronously.
    engine.syncSequenceSelectionNow()
  }

  /// One step of a drag: 'sele' becomes the drag-start snapshot plus (or
  /// minus) the residues currently spanned, quietly.
  private func applyDragSpan(_ residues: [SequenceResidue], add: Bool) {
    let residues = residues.filter { $0.isSelectable }
    guard !residues.isEmpty else { return }
    let expr = selectionExpression(residues)
    var py = Self.modeExpression(expr)
    py +=
      add
      ? "_c.select('sele', '(_seqdrag0) or (' + _e + ')', enable=1, quiet=1)\n"
      : "_c.select('sele', '(_seqdrag0) and not (' + _e + ')', enable=1, quiet=1)\n"
    engine.runPython(py)
  }

  /// Python that leaves in `_e` the residue expression widened by the current
  /// selection mode (mouse_selection_mode, the panel's "Selecting …"), so a
  /// click in the strip picks what a click in the viewport would: Chains →
  /// the whole chain, Objects → the object, C-alphas → just the CAs.
  static func modeExpression(_ expr: String) -> String {
    "from pymol import cmd as _c\n"
      + "_e = '\(expr)'\n"
      + "_e = {2: 'bychain (%s)', 3: 'bysegi (%s)', 4: 'byobject (%s)', 5: 'bymolecule (%s)',\n"
      + "      6: '(%s) and name CA'}.get(_c.get_setting_int('mouse_selection_mode'), '%s') % _e\n"
  }

  /// Build a PyMOL selection expression, grouping residues by object+chain
  /// and listing residue numbers (robust to gaps).
  private func selectionExpression(_ residues: [SequenceResidue]) -> String {
    var groups: [String: (obj: String, chain: String, resis: [String])] = [:]
    for r in residues {
      let key = "\(r.objectName)|\(r.chain)"
      if groups[key] == nil {
        groups[key] = (r.objectName, r.chain, [])
      }
      groups[key]!.resis.append(r.resi)
    }
    let parts = groups.values.map { g -> String in
      let resiList = g.resis.joined(separator: "+")
      if g.chain.isEmpty {
        return "(\(g.obj) and resi \(resiList))"
      }
      return "(\(g.obj) and chain \(g.chain) and resi \(resiList))"
    }
    return parts.joined(separator: " or ")
  }
}

// MARK: - PyMOLEngine extensions for sequence polling

extension PyMOLEngine {
  /// Parse the SEQPANEL JSON emitted by fetchSequences() into `sequences`.
  /// Payload: { "objects": [{ "name", "residues": [[chain, resi, resn, colorIdx], ...] }],
  ///            "colors": { "<idx>": [r, g, b] } }
  func parseSequencePanelFeedback(_ line: String) {
    guard line.hasPrefix("SEQPANEL:") else { return }
    let path = TempChannel.path(TempChannel.Stem.sequence)
    guard let data = FileManager.default.contents(atPath: path) else { return }

    struct Payload: Decodable {
      struct Obj: Decodable {
        let name: String
        let residues: [[String]]
      }
      let objects: [Obj]
      let colors: [String: [Double]]
    }

    guard let p = try? JSONDecoder().decode(Payload.self, from: data) else { return }

    var objs: [SequenceObject] = []
    for o in p.objects {
      var residues: [SequenceResidue] = []
      var prevChain: String? = nil  // last real residue's chain (for breaks)
      for (i, t) in o.residues.enumerated() where t.count >= 4 {
        let chain = t[0]
        let resi = t[1]
        let resn = t[2]
        let cidx = t[3]
        // Alignment gap placeholder (emitted by appkit_sequence): a dim,
        // non-selectable dash that keeps aligned residues column-aligned.
        if resn == "-" {
          residues.append(
            SequenceResidue(
              id: "\(o.name)/gap/\(i)",
              objectName: o.name,
              chain: "", oneLetter: "-", resi: "", resn: "-",
              color: Color.gray.opacity(0.35),
              isGap: true))
          continue
        }
        // Chain boundary / first chain: insert a separator bar (between
        // chains only) and a chain-ID label at the start of each chain run
        // (skipped for blank/unnamed chains). Both are non-selectable.
        let chainChanged = (prevChain != chain)
        if chainChanged {
          if prevChain != nil {
            residues.append(
              SequenceResidue(
                id: "\(o.name)/brk/\(i)",
                objectName: o.name,
                chain: "", oneLetter: "", resi: "", resn: "",
                color: rulerColor, isBreak: true))
          }
          if !chain.isEmpty {
            residues.append(
              SequenceResidue(
                id: "\(o.name)/chid/\(i)",
                objectName: o.name,
                chain: chain, oneLetter: "", resi: "", resn: "",
                color: headerColor, isChainLabel: true, chainLabel: chain))
          }
        }
        prevChain = chain
        let rgb = p.colors[cidx] ?? [0.8, 0.8, 0.8]
        let color = Color(
          .sRGB,
          red: rgb.count > 0 ? rgb[0] : 0.8,
          green: rgb.count > 1 ? rgb[1] : 0.8,
          blue: rgb.count > 2 ? rgb[2] : 0.8)
        // HETATM groups (ligands, ions, waters) are tagged 'het' by
        // appkit_sequence. They have no 1-letter code, so spread the
        // 3-letter resn (SY7, EDO, HOH…) across adjacent columns like
        // PyMOL — every cell shares one selKey (obj/chain/resi), so a
        // click on any character selects the whole group (issue #201).
        if t.count >= 5 && t[4] == "het" {
          // One blank column before each group, as desktop PyMOL's
          // viewer separates entities ("…KVGE HEM SO4 HOH HOH").
          residues.append(
            SequenceResidue(
              id: "\(o.name)/sp/\(i)", objectName: o.name,
              chain: "", oneLetter: "", resi: "", resn: "",
              color: .clear, isSpacer: true))
          let letters = resn.isEmpty ? ["?"] : resn.map { String($0) }
          for (j, ch) in letters.enumerated() {
            residues.append(
              SequenceResidue(
                id: "\(o.name)/\(chain)/\(resi)/\(i)/\(j)",
                objectName: o.name,
                chain: chain,
                oneLetter: ch,
                resi: resi,
                resn: resn,
                color: color,
                isHet: true
              ))
          }
          continue
        }
        residues.append(
          SequenceResidue(
            id: "\(o.name)/\(chain)/\(resi)/\(i)",
            objectName: o.name,
            chain: chain,
            oneLetter: aa3to1[resn.uppercased()] ?? "X",
            resi: resi,
            resn: resn,
            color: color
          ))
      }
      if !residues.isEmpty {
        objs.append(SequenceObject(id: o.name, name: o.name, residues: residues))
      }
    }

    DispatchQueue.main.async {
      self.sequences = objs
    }
  }

  /// Parse the SEQSEL JSON (active-selection residue keys) into
  /// `selectedResidueKeys` for highlight sync.
  /// Re-read 'sele' for the strip right now: the Python runs synchronously,
  /// so its file is already written when it returns, and there is no need to
  /// wait for the next feedback poll (up to 100 ms) to pick up the marker.
  func syncSequenceSelectionNow() {
    fetchSequenceSelection()
    parseSequenceSelectionFeedback("SEQSEL:ready")
  }

  func parseSequenceSelectionFeedback(_ line: String) {
    guard line.hasPrefix("SEQSEL:") else { return }
    let path = TempChannel.path(TempChannel.Stem.sequenceSelection)
    guard let data = FileManager.default.contents(atPath: path) else { return }
    guard let keys = try? JSONDecoder().decode([String].self, from: data) else { return }
    let set = Set(keys)
    DispatchQueue.main.async {
      if self.selectedResidueKeys != set {
        self.selectedResidueKeys = set
      }
    }
  }
}

// MARK: - Strip rows

/// One ruler cell (a residue-number digit or a chain ID) above a residue slot.
struct ClassicRulerCell: Equatable {
  let text: String
  let isChain: Bool
}

/// One object's ruler + residue rows, drawn in a single Canvas the way desktop
/// PyMOL's seq_view draws its strip. A view per residue meant thousands of
/// views re-laid-out
/// on every engine change; this redraws only when its own inputs change
/// (`.equatable()`), and hit-testing is arithmetic on `cellWidth`.
struct ClassicSequenceRows: View, Equatable {
  let name: String
  let residues: [SequenceResidue]
  let ruler: [ClassicRulerCell]
  let selected: [Bool]
  let labelWidth: CGFloat
  let cellWidth: CGFloat
  let rowHeight: CGFloat
  let labelColor: Color
  let rulerColor: Color
  let selectionBG: Color
  let selectionFG: Color

  static let font = Font.system(size: 12, design: .monospaced)
  static let rulerFont = Font.system(size: 11, design: .monospaced)

  /// Residues per Canvas. A row is drawn as side-by-side chunks: one Canvas
  /// as wide as a long sequence (3VEY with waters: ~19,000 pt) exceeds the
  /// largest GPU texture, so it was re-rendered while scrolling instead of
  /// just moved. 512 residues stay well under the limit.
  static let chunk = 512

  var body: some View {
    HStack(spacing: 0) {
      ForEach(Array(stride(from: 0, to: max(residues.count, 1), by: Self.chunk)), id: \.self) {
        start in
        chunkCanvas(start, min(start + Self.chunk, residues.count))
      }
    }
    .frame(
      width: labelWidth + CGFloat(residues.count) * cellWidth, height: rowHeight * 2,
      alignment: .leading)
  }

  /// Residues start..<end (and, for the first chunk, the name gutter).
  private func chunkCanvas(_ start: Int, _ end: Int) -> some View {
    let x0 = start == 0 ? labelWidth : 0
    return Canvas { ctx, _ in
      let midRuler = rowHeight / 2
      let midRes = rowHeight * 1.5
      // each letter is laid out once per draw, not once per residue
      var glyphs: [String: GraphicsContext.ResolvedText] = [:]
      func glyph(_ s: String, _ font: Font, _ color: Color) -> GraphicsContext.ResolvedText {
        let key = "\(s)|\(color)|\(font == Self.font)"
        if let g = glyphs[key] { return g }
        let g = ctx.resolve(Text(s).font(font).foregroundColor(color))
        glyphs[key] = g
        return g
      }
      if start == 0 {
        ctx.draw(
          Text(name).font(Self.font).foregroundColor(labelColor),
          at: CGPoint(x: 0, y: midRes), anchor: .leading)
      }
      for i in start..<min(end, ruler.count) where ruler[i].text != " " && !ruler[i].text.isEmpty {
        let cell = ruler[i]
        ctx.draw(
          glyph(cell.text, Self.rulerFont, cell.isChain ? labelColor : rulerColor),
          at: CGPoint(x: x0 + CGFloat(i - start) * cellWidth, y: midRuler), anchor: .leading)
      }
      for i in start..<end {
        let r = residues[i]
        let x = x0 + CGFloat(i - start) * cellWidth
        if r.isBreak {
          ctx.fill(
            Path(
              CGRect(
                x: x + cellWidth / 2, y: rowHeight + 1, width: 1,
                height: rowHeight - 2)),
            with: .color(rulerColor.opacity(0.5)))
          continue
        }
        if r.isChainLabel || r.isSpacer || r.oneLetter.isEmpty { continue }
        let on = i < selected.count && selected[i]
        if on {
          ctx.fill(
            Path(CGRect(x: x, y: rowHeight, width: cellWidth, height: rowHeight)),
            with: .color(selectionBG))
        }
        ctx.draw(
          glyph(r.oneLetter, Self.font, on ? selectionFG : r.color),
          at: CGPoint(x: x, y: midRes), anchor: .leading)
      }
    }
    .frame(width: x0 + CGFloat(end - start) * cellWidth, height: rowHeight * 2)
  }
}
