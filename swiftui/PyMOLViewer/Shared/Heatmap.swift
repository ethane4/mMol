// Heatmap.swift — C ▸ by file…: colour a structure from a column of a CSV.
//
// The CSV is read here (header row, then one row per residue). Each object keeps
// its own CSV and settings; the object's ▸ card shows a "heatmap" chip whose panel
// (HeatmapPanel, ObjectPanel.swift) edits the colours, range, value column and
// index column, and every change recolours the object.
//
// Residue mapping:
//   • no index column: row k (1-based, after the header) is residue number k,
//     in every chain;
//   • an index column: each value names the residue — "45" (residue 45 in every
//     chain) or "A:45" / "A45" (residue 45 of chain A).
// Residues the CSV does not cover are grey.

import AppKit
import Foundation
import SwiftUI

/// A parsed CSV: column names and the raw cell text of every row.
struct HeatmapCSV: Equatable {
  let fileName: String
  let headers: [String]
  let rows: [[String]]

  /// Parse comma- (or, failing that, tab-) separated text with a header row.
  /// Handles quoted fields with embedded commas and doubled quotes.
  static func parse(_ text: String, fileName: String) -> HeatmapCSV? {
    var lines = text.components(separatedBy: .newlines)
      .map { $0.trimmingCharacters(in: CharacterSet(charactersIn: "\r")) }
      .filter { !$0.trimmingCharacters(in: .whitespaces).isEmpty }
    guard !lines.isEmpty else { return nil }
    let sep: Character = (lines[0].contains(",") || !lines[0].contains("\t")) ? "," : "\t"
    let headers = split(lines.removeFirst(), sep)
    guard !headers.isEmpty else { return nil }
    let rows = lines.map { line -> [String] in
      var r = split(line, sep)
      if r.count < headers.count { r += Array(repeating: "", count: headers.count - r.count) }
      return r
    }
    return HeatmapCSV(fileName: fileName, headers: headers, rows: rows)
  }

  private static func split(_ line: String, _ sep: Character) -> [String] {
    var out: [String] = []
    var cur = ""
    var quoted = false
    let it = Array(line)
    var i = 0
    while i < it.count {
      let c = it[i]
      if quoted {
        if c == "\"" {
          if i + 1 < it.count, it[i + 1] == "\"" {
            cur.append("\"")
            i += 1
          } else {
            quoted = false
          }
        } else {
          cur.append(c)
        }
      } else if c == "\"" {
        quoted = true
      } else if c == sep {
        out.append(cur.trimmingCharacters(in: .whitespaces))
        cur = ""
      } else {
        cur.append(c)
      }
      i += 1
    }
    out.append(cur.trimmingCharacters(in: .whitespaces))
    return out
  }

  func values(of column: String) -> [Double?] {
    guard let c = headers.firstIndex(of: column) else { return [] }
    return rows.map { c < $0.count ? Double($0[c]) : nil }
  }

  /// Columns where most cells are numbers — the ones a heatmap can use.
  var numericColumns: [String] {
    headers.filter { h in
      let v = values(of: h)
      return !v.isEmpty && v.compactMap { $0 }.count * 2 >= v.count
    }
  }
}

struct HeatmapSettings: Equatable {
  var column: String
  var indexColumn: String?  // nil: row k is residue k
  var low: RGBA = RGBA(0.13, 0.35, 0.85)
  var high: RGBA = RGBA(0.85, 0.15, 0.15)
  var throughWhite: Bool = true  // diverging: low → white → high
  var minValue: Double
  var maxValue: Double
}

struct HeatmapState: Equatable {
  var csv: HeatmapCSV
  var settings: HeatmapSettings

  /// min / max of the current value column (the range slider's ends).
  var dataRange: ClosedRange<Double> {
    let v = csv.values(of: settings.column).compactMap { $0 }
    guard let lo = v.min(), let hi = v.max() else { return 0...1 }
    return lo...(hi > lo ? hi : lo + 1)
  }

  /// The colour for `value` under the current settings.
  func color(for value: Double) -> RGBA {
    let s = settings
    let span = s.maxValue - s.minValue
    let t = span > 0 ? min(1, max(0, (value - s.minValue) / span)) : 0.5
    if s.throughWhite {
      let white = RGBA(1, 1, 1)
      return t < 0.5
        ? s.low.blended(with: white, t * 2) : white.blended(with: s.high, (t - 0.5) * 2)
    }
    return s.low.blended(with: s.high, t)
  }

  /// Residue selection fragment per row value: [(chain or nil, resi, value)].
  var residueValues: [(chain: String?, resi: String, value: Double)] {
    let vals = csv.values(of: settings.column)
    let idx = settings.indexColumn.flatMap { csv.headers.firstIndex(of: $0) }
    var out: [(String?, String, Double)] = []
    for (k, row) in csv.rows.enumerated() {
      guard k < vals.count, let v = vals[k] else { continue }
      if let idx {
        guard idx < row.count, let key = Self.residueKey(row[idx]) else { continue }
        out.append((key.chain, key.resi, v))
      } else {
        out.append((nil, String(k + 1), v))
      }
    }
    return out
  }

  /// "45" → (nil, "45"); "A:45", "A45", "A_45", "A 45" → ("A", "45").
  static func residueKey(_ raw: String) -> (chain: String?, resi: String)? {
    let s = raw.trimmingCharacters(in: .whitespaces)
    guard !s.isEmpty else { return nil }
    if s.allSatisfy({ $0.isNumber || $0 == "-" }) { return (nil, s) }
    let parts = s.split(whereSeparator: { $0 == ":" || $0 == "_" || $0 == " " })
    if parts.count == 2 { return (String(parts[0]), String(parts[1])) }
    // "A45": leading letters are the chain.
    let letters = s.prefix { $0.isLetter }
    let rest = s.dropFirst(letters.count)
    if !letters.isEmpty, !rest.isEmpty { return (String(letters), String(rest)) }
    return (nil, s)
  }
}

extension PyMOLEngine {
  /// C ▸ by file…: pick a CSV, attach it to `object`, colour by its first
  /// numeric column, and open the object's card on the heatmap chip.
  func loadHeatmapCSV(for object: String) {
    let panel = NSOpenPanel()
    panel.allowedContentTypes = [.commaSeparatedText, .tabSeparatedText, .plainText]
    panel.allowsMultipleSelection = false
    panel.title = "Colour \(object) by file"
    panel.message = "Choose a CSV with a header row and one row per residue."
    guard panel.runModal() == .OK, let url = panel.url else { return }
    attachHeatmap(object, url: url)
  }

  /// Attach the CSV at `url` to `object` and colour it (also the test hook's
  /// entry point, PYMOL_AUTOHEATMAP=<object>|<path>).
  func attachHeatmap(_ object: String, url: URL) {
    guard let text = try? String(contentsOf: url, encoding: .utf8) else {
      logLine(" Heatmap: couldn't read \(url.lastPathComponent).")
      return
    }
    guard let csv = HeatmapCSV.parse(text, fileName: url.lastPathComponent),
      let column = csv.numericColumns.first
    else {
      logLine(" Heatmap: \(url.lastPathComponent) has no numeric column to colour by.")
      return
    }
    // A first column that looks like residue numbers / IDs is the index.
    let first = csv.headers[0].lowercased()
    let index: String? =
      ["resi", "residue", "res", "index", "resid", "position", "pos"]
        .contains(where: { first.hasPrefix($0) }) ? csv.headers[0] : nil
    let valueColumn =
      index != nil
      ? (csv.numericColumns.first { $0 != index } ?? column) : column
    var state = HeatmapState(
      csv: csv,
      settings: HeatmapSettings(
        column: valueColumn, indexColumn: index, minValue: 0, maxValue: 1))
    let r = state.dataRange
    state.settings.minValue = r.lowerBound
    state.settings.maxValue = r.upperBound
    heatmaps[object] = state
    applyHeatmap(object)
    expandedDetail = object
    heatmapChipActive.insert(object)
  }

  /// Update an object's heatmap settings and recolour it.
  func updateHeatmap(_ object: String, _ change: (inout HeatmapSettings) -> Void) {
    guard var st = heatmaps[object] else { return }
    let oldColumn = st.settings.column
    change(&st.settings)
    if st.settings.column != oldColumn {  // new column: range = its data
      heatmaps[object] = st
      let r = st.dataRange
      st.settings.minValue = r.lowerBound
      st.settings.maxValue = r.upperBound
    }
    heatmaps[object] = st
    applyHeatmap(object)
  }

  /// Colour `object` from its heatmap: polymer residues the CSV covers get
  /// their value's colour (quantised to 128 steps, one `color` per step),
  /// the rest grey.
  func applyHeatmap(_ object: String) {
    guard let st = heatmaps[object] else { return }
    let steps = 128
    let s = st.settings
    let span = s.maxValue - s.minValue
    var buckets: [Int: [(String?, String)]] = [:]
    for rv in st.residueValues {
      let t = span > 0 ? min(1, max(0, (rv.value - s.minValue) / span)) : 0.5
      buckets[Int((t * Double(steps - 1)).rounded()), default: []].append((rv.chain, rv.resi))
    }
    let tag = sanitizedName(object)
    let obj = object.replacingOccurrences(of: "'", with: "")
    var py = "from pymol import cmd as _c\n"
    py += "_c.color('gray70', \"(\(obj)) and polymer\", quiet=1)\n"
    for (k, members) in buckets.sorted(by: { $0.key < $1.key }) {
      let c = st.color(for: s.minValue + span * Double(k) / Double(steps - 1))
      let name = "_hm_\(tag)_\(k)"
      var noChain: [String] = []
      var byChain: [String: [String]] = [:]
      for (ch, resi) in members {
        if let ch { byChain[ch, default: []].append(resi) } else { noChain.append(resi) }
      }
      var parts: [String] = []
      if !noChain.isEmpty { parts.append("resi " + noChain.joined(separator: "+")) }
      for (ch, rs) in byChain {
        parts.append("(chain \(ch) and resi \(rs.joined(separator: "+")))")
      }
      py += "_c.set_color('\(name)', [\(c.r), \(c.g), \(c.b)], quiet=1)\n"
      py +=
        "_c.color('\(name)', \"(\(obj)) and polymer and (\(parts.joined(separator: " or ")))\", quiet=1)\n"
    }
    runPython(py)
  }
}
