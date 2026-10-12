// PyMOLColorMenu.swift — the "C" menu's coloring modes and desktop PyMOL's
// tiered named-color table.

import SwiftUI

/// One named PyMOL color and its sRGB value, mirrored from the core color
/// table (src/pymol/settings/color.cc) so the menu swatch matches what `color <name>`
/// actually paints.
struct PyMOLNamedColor: Equatable {
  let name: String
  let red: Double
  let green: Double
  let blue: Double
  var swatch: Color { Color(.sRGB, red: red, green: green, blue: blue) }
}

/// One hue family of desktop PyMOL's two-level color menu (#379).
struct PyMOLColorFamily {
  let label: String
  let shades: [PyMOLNamedColor]
  /// The family entry is tinted with its first (canonical) shade, exactly
  /// like upstream menu.py, which colors the family label with `c_list[0]`.
  var swatch: Color { shades[0].swatch }
}

/// Desktop PyMOL's tiered color menu — `all_colors_list` in
/// modules/pymol/menu.py — surfaced in the SwiftUI "C" menu as one submenu per
/// hue family (#379). Family order and membership mirror upstream (a few
/// shades, e.g. limon and wheat, deliberately appear in two families there);
/// the RGB values are the core's, so no colors are invented here.
enum PyMOLColorMenu {
  static let families: [PyMOLColorFamily] = [
    PyMOLColorFamily(
      label: "reds",
      shades: [
        PyMOLNamedColor(name: "red", red: 1.0, green: 0.0, blue: 0.0),
        PyMOLNamedColor(name: "tv_red", red: 1.0, green: 0.2, blue: 0.2),
        PyMOLNamedColor(name: "raspberry", red: 0.7, green: 0.3, blue: 0.4),
        PyMOLNamedColor(name: "darksalmon", red: 0.73, green: 0.55, blue: 0.52),
        PyMOLNamedColor(name: "salmon", red: 1.0, green: 0.6, blue: 0.6),
        PyMOLNamedColor(name: "deepsalmon", red: 1.0, green: 0.5, blue: 0.5),
        PyMOLNamedColor(name: "warmpink", red: 0.85, green: 0.2, blue: 0.5),
        PyMOLNamedColor(name: "firebrick", red: 0.698, green: 0.13, blue: 0.13),
        PyMOLNamedColor(name: "ruby", red: 0.6, green: 0.2, blue: 0.2),
        PyMOLNamedColor(name: "chocolate", red: 0.555, green: 0.222, blue: 0.111),
        PyMOLNamedColor(name: "brown", red: 0.65, green: 0.32, blue: 0.17),
      ]),
    PyMOLColorFamily(
      label: "greens",
      shades: [
        PyMOLNamedColor(name: "green", red: 0.0, green: 1.0, blue: 0.0),
        PyMOLNamedColor(name: "tv_green", red: 0.2, green: 1.0, blue: 0.2),
        PyMOLNamedColor(name: "chartreuse", red: 0.5, green: 1.0, blue: 0.0),
        PyMOLNamedColor(name: "splitpea", red: 0.52, green: 0.75, blue: 0.0),
        PyMOLNamedColor(name: "smudge", red: 0.55, green: 0.7, blue: 0.4),
        PyMOLNamedColor(name: "palegreen", red: 0.65, green: 0.9, blue: 0.65),
        PyMOLNamedColor(name: "limegreen", red: 0.0, green: 1.0, blue: 0.5),
        PyMOLNamedColor(name: "lime", red: 0.5, green: 1.0, blue: 0.5),
        PyMOLNamedColor(name: "limon", red: 0.75, green: 1.0, blue: 0.25),
        PyMOLNamedColor(name: "forest", red: 0.2, green: 0.6, blue: 0.2),
      ]),
    PyMOLColorFamily(
      label: "blues",
      shades: [
        PyMOLNamedColor(name: "blue", red: 0.0, green: 0.0, blue: 1.0),
        PyMOLNamedColor(name: "tv_blue", red: 0.3, green: 0.3, blue: 1.0),
        PyMOLNamedColor(name: "marine", red: 0.0, green: 0.5, blue: 1.0),
        PyMOLNamedColor(name: "slate", red: 0.5, green: 0.5, blue: 1.0),
        PyMOLNamedColor(name: "lightblue", red: 0.75, green: 0.75, blue: 1.0),
        PyMOLNamedColor(name: "skyblue", red: 0.2, green: 0.5, blue: 0.8),
        PyMOLNamedColor(name: "purpleblue", red: 0.5, green: 0.0, blue: 1.0),
        PyMOLNamedColor(name: "deepblue", red: 0.25, green: 0.25, blue: 0.65),
        PyMOLNamedColor(name: "density", red: 0.1, green: 0.1, blue: 0.6),
      ]),
    PyMOLColorFamily(
      label: "yellows",
      shades: [
        PyMOLNamedColor(name: "yellow", red: 1.0, green: 1.0, blue: 0.0),
        PyMOLNamedColor(name: "tv_yellow", red: 1.0, green: 1.0, blue: 0.2),
        PyMOLNamedColor(name: "paleyellow", red: 1.0, green: 1.0, blue: 0.5),
        PyMOLNamedColor(name: "yelloworange", red: 1.0, green: 0.87, blue: 0.37),
        PyMOLNamedColor(name: "limon", red: 0.75, green: 1.0, blue: 0.25),
        PyMOLNamedColor(name: "wheat", red: 0.99, green: 0.82, blue: 0.65),
        PyMOLNamedColor(name: "sand", red: 0.72, green: 0.55, blue: 0.3),
      ]),
    PyMOLColorFamily(
      label: "magentas",
      shades: [
        PyMOLNamedColor(name: "magenta", red: 1.0, green: 0.0, blue: 1.0),
        PyMOLNamedColor(name: "lightmagenta", red: 1.0, green: 0.2, blue: 0.8),
        PyMOLNamedColor(name: "hotpink", red: 1.0, green: 0.0, blue: 0.5),
        PyMOLNamedColor(name: "pink", red: 1.0, green: 0.65, blue: 0.85),
        PyMOLNamedColor(name: "lightpink", red: 1.0, green: 0.75, blue: 0.87),
        PyMOLNamedColor(name: "dirtyviolet", red: 0.7, green: 0.5, blue: 0.5),
        PyMOLNamedColor(name: "violet", red: 1.0, green: 0.5, blue: 1.0),
        PyMOLNamedColor(name: "violetpurple", red: 0.55, green: 0.25, blue: 0.6),
        PyMOLNamedColor(name: "purple", red: 0.75, green: 0.0, blue: 0.75),
        PyMOLNamedColor(name: "deeppurple", red: 0.6, green: 0.1, blue: 0.6),
      ]),
    PyMOLColorFamily(
      label: "cyans",
      shades: [
        PyMOLNamedColor(name: "cyan", red: 0.0, green: 1.0, blue: 1.0),
        PyMOLNamedColor(name: "palecyan", red: 0.8, green: 1.0, blue: 1.0),
        PyMOLNamedColor(name: "aquamarine", red: 0.5, green: 1.0, blue: 1.0),
        PyMOLNamedColor(name: "greencyan", red: 0.25, green: 1.0, blue: 0.75),
        PyMOLNamedColor(name: "teal", red: 0.0, green: 0.75, blue: 0.75),
        PyMOLNamedColor(name: "deepteal", red: 0.1, green: 0.6, blue: 0.6),
        PyMOLNamedColor(name: "lightteal", red: 0.4, green: 0.7, blue: 0.7),
      ]),
    PyMOLColorFamily(
      label: "oranges",
      shades: [
        PyMOLNamedColor(name: "orange", red: 1.0, green: 0.5, blue: 0.0),
        PyMOLNamedColor(name: "tv_orange", red: 1.0, green: 0.55, blue: 0.15),
        PyMOLNamedColor(name: "brightorange", red: 1.0, green: 0.7, blue: 0.2),
        PyMOLNamedColor(name: "lightorange", red: 1.0, green: 0.8, blue: 0.5),
        PyMOLNamedColor(name: "yelloworange", red: 1.0, green: 0.87, blue: 0.37),
        PyMOLNamedColor(name: "olive", red: 0.77, green: 0.7, blue: 0.0),
        PyMOLNamedColor(name: "deepolive", red: 0.6, green: 0.6, blue: 0.1),
      ]),
    PyMOLColorFamily(
      label: "tints",
      shades: [
        PyMOLNamedColor(name: "wheat", red: 0.99, green: 0.82, blue: 0.65),
        PyMOLNamedColor(name: "palegreen", red: 0.65, green: 0.9, blue: 0.65),
        PyMOLNamedColor(name: "lightblue", red: 0.75, green: 0.75, blue: 1.0),
        PyMOLNamedColor(name: "paleyellow", red: 1.0, green: 1.0, blue: 0.5),
        PyMOLNamedColor(name: "lightpink", red: 1.0, green: 0.75, blue: 0.87),
        PyMOLNamedColor(name: "palecyan", red: 0.8, green: 1.0, blue: 1.0),
        PyMOLNamedColor(name: "lightorange", red: 1.0, green: 0.8, blue: 0.5),
        PyMOLNamedColor(name: "bluewhite", red: 0.85, green: 0.85, blue: 1.0),
      ]),
    PyMOLColorFamily(
      label: "grays",
      shades: [
        PyMOLNamedColor(name: "white", red: 1.0, green: 1.0, blue: 1.0),
        PyMOLNamedColor(name: "gray90", red: 0.909, green: 0.909, blue: 0.909),
        PyMOLNamedColor(name: "gray80", red: 0.808, green: 0.808, blue: 0.808),
        PyMOLNamedColor(name: "gray70", red: 0.707, green: 0.707, blue: 0.707),
        PyMOLNamedColor(name: "gray60", red: 0.606, green: 0.606, blue: 0.606),
        PyMOLNamedColor(name: "gray50", red: 0.505, green: 0.505, blue: 0.505),
        PyMOLNamedColor(name: "gray40", red: 0.404, green: 0.404, blue: 0.404),
        PyMOLNamedColor(name: "gray30", red: 0.303, green: 0.303, blue: 0.303),
        PyMOLNamedColor(name: "gray20", red: 0.202, green: 0.202, blue: 0.202),
        PyMOLNamedColor(name: "gray10", red: 0.101, green: 0.101, blue: 0.101),
        PyMOLNamedColor(name: "black", red: 0.0, green: 0.0, blue: 0.0),
      ]),
  ]

  /// Every distinct named color in the table, in menu order.
  static var allColorNames: [String] {
    var seen = Set<String>()
    return families.flatMap(\.shades).map(\.name).filter { seen.insert($0).inserted }
  }

  /// Look up a named color anywhere in the tiered table.
  static func color(named name: String) -> PyMOLNamedColor? {
    for family in families {
      if let c = family.shades.first(where: { $0.name == name }) { return c }
    }
    return nil
  }

  /// Classic `gray` (0.5) is not in upstream's family table (which lists the
  /// gray10…gray90 ramp), but it is the color people reach for, so it heads
  /// the grays row.
  static let gray = PyMOLNamedColor(name: "gray", red: 0.5, green: 0.5, blue: 0.5)

  /// The rows of the "C" menu's color section (#379). Like desktop PyMOL's
  /// object-panel popup, each row IS a color — clicking "red" paints red —
  /// and the same row expands (arrow) into that hue's named variants. Rows
  /// without a natural base color (tints) are plain submenus.
  static let rows: [PyMOLColorMenuRow] = {
    let base: [String: PyMOLNamedColor?] = [
      "reds": color(named: "red"),
      "greens": color(named: "green"),
      "blues": color(named: "blue"),
      "yellows": color(named: "yellow"),
      "magentas": color(named: "magenta"),
      "cyans": color(named: "cyan"),
      "oranges": color(named: "orange"),
      "tints": nil,
      "grays": gray,
    ]
    return families.map { PyMOLColorMenuRow(family: $0, primary: base[$0.label] ?? nil) }
  }()
}

/// One row of the color section: a hue family plus the color a click on the
/// row itself applies (nil = the row only expands).
struct PyMOLColorMenuRow {
  let family: PyMOLColorFamily
  let primary: PyMOLNamedColor?
  /// What the row reads: the base color's name when clickable, else the
  /// family's name.
  var title: String { primary?.name ?? family.label }
  var swatch: Color { primary?.swatch ?? family.swatch }
}

/// Color options with optional swatch color
struct ColorOption {
  let label: String
  let command: String?
  let swatch: Color?

}

/// The coloring modes at the top of the "C" menu. The colors themselves follow
/// as PyMOLColorMenu.rows (click = base color, expand = its variants), so the
/// old flat red…white list is gone — it duplicated the rows (#379).
let colorOptions: [ColorOption] = [
  ColorOption(label: "by element", command: "util.cnc", swatch: nil),
  ColorOption(label: "by chain", command: "util.cbc", swatch: nil),
  ColorOption(label: "by ss", command: "util.cbss", swatch: nil),
  ColorOption(label: "spectrum", command: "spectrum", swatch: nil),
  // Desktop PyMOL's C ▸ spectrum ▸ rainbow and C ▸ by chain ▸ chainbows
  // (modules/pymol/menu.py): blue at the N-terminus through to red at the C.
  ColorOption(label: "rainbow", command: "rainbow", swatch: nil),
  ColorOption(label: "chainbows", command: "util.chainbow", swatch: nil),
  ColorOption(label: "by b-factor", command: "spectrum_b", swatch: nil),
  // Colour from a CSV column; edited afterwards in the object's ▸ heatmap chip.
  ColorOption(label: "by file…", command: "by_file", swatch: nil),
]
