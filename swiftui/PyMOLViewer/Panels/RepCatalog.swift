// RepCatalog.swift — the inspector card's per-representation controls and
// the Viewer Settings scene parameters.

import Foundation

// MARK: - Representation inspector: control metadata

enum RepControlKind {
  case slider, segmented, toggle, color
}

/// One controllable property row (label + control bound to a PyMOL setting).
struct RepProperty: Identifiable {
  var id: String { setting }
  let setting: String
  let label: String
  let kind: RepControlKind
  var min: Double = 0
  var max: Double = 1
  var step: Double = 0.01
  var decimals: Int = 2
  var options: [(label: String, value: Double)] = []  // for .segmented
  // Apply only on release (not on every live drag tick). For settings whose
  // change forces an expensive rebuild (e.g. solvent_radius re-tessellates the
  // whole surface), live updates would recompute on every drag step.
  var commitOnly: Bool = false
}

/// Static description of a representation: display name, color-override setting
/// (empty = no per-rep color), and the property rows it exposes. Setting names
/// MUST match modules/mmol/inspector.py.
struct RepSpec {
  let rep: String
  let display: String
  let colorSetting: String  // e.g. "surface_color"; "" if none
  let defaultColor: Int  // value meaning "inherit" (-1, labels -6)
  let properties: [RepProperty]
}

enum RepCatalog {
  static let order = [
    "cartoon", "surface", "sticks", "spheres", "ribbon",
    "mesh", "lines", "dots", "nonbonded", "nb_spheres", "labels",
  ]

  static let specs: [String: RepSpec] = [
    "cartoon": RepSpec(
      rep: "cartoon", display: "Cartoon",
      colorSetting: "cartoon_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "cartoon_transparency", label: "Transparency", kind: .slider),
        RepProperty(setting: "cartoon_loop_radius", label: "Loop radius", kind: .slider),
        RepProperty(setting: "cartoon_tube_radius", label: "Tube radius", kind: .slider),
        RepProperty(setting: "cartoon_fancy_helices", label: "Fancy helices", kind: .toggle),
        RepProperty(setting: "cartoon_flat_sheets", label: "Flat sheets", kind: .toggle),
        RepProperty(setting: "cartoon_spline", label: "Spline ribbon", kind: .toggle),
      ]),
    "surface": RepSpec(
      rep: "surface", display: "Surface",
      colorSetting: "surface_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "transparency", label: "Transparency", kind: .slider),
        RepProperty(
          setting: "surface_quality", label: "Quality", kind: .segmented,
          options: [("0", 0), ("1", 1), ("2", 2)]),
        RepProperty(
          setting: "solvent_radius", label: "Solvent radius", kind: .slider, min: 0.5, max: 3,
          step: 0.1, decimals: 1, commitOnly: true),
        RepProperty(
          setting: "surface_clip_front", label: "Clip front", kind: .slider, min: 0, max: 1,
          step: 0.02, decimals: 2),
        RepProperty(
          setting: "surface_clip_back", label: "Clip back", kind: .slider, min: 0, max: 1,
          step: 0.02, decimals: 2),
        RepProperty(setting: "metal_interior_cap", label: "Solid interior", kind: .toggle),
        RepProperty(setting: "surface_contour", label: "Contour", kind: .toggle),
        RepProperty(
          setting: "surface_contour_width", label: "Contour width", kind: .slider, min: 0.5, max: 6,
          step: 0.5, decimals: 1),
        RepProperty(setting: "surface_contour_color", label: "Contour color", kind: .color),
        RepProperty(setting: "surface_contour_opaque", label: "Contour opaque", kind: .toggle),
      ]),
    "sticks": RepSpec(
      rep: "sticks", display: "Sticks",
      colorSetting: "stick_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "stick_transparency", label: "Transparency", kind: .slider),
        RepProperty(setting: "stick_radius", label: "Radius", kind: .slider),
        RepProperty(setting: "stick_h_scale", label: "H scale", kind: .slider),
        RepProperty(setting: "metal_interior_cap", label: "Solid interior", kind: .toggle),
      ]),
    "spheres": RepSpec(
      rep: "spheres", display: "Spheres",
      colorSetting: "sphere_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "sphere_transparency", label: "Transparency", kind: .slider),
        RepProperty(setting: "sphere_scale", label: "Scale", kind: .slider, max: 3, step: 0.05),
        RepProperty(setting: "metal_interior_cap", label: "Solid interior", kind: .toggle),
      ]),
    "ribbon": RepSpec(
      rep: "ribbon", display: "Ribbon",
      colorSetting: "ribbon_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "ribbon_transparency", label: "Transparency", kind: .slider),
        RepProperty(
          setting: "ribbon_width", label: "Width", kind: .slider, max: 20, step: 0.1, decimals: 1),
      ]),
    "mesh": RepSpec(
      rep: "mesh", display: "Mesh",
      colorSetting: "mesh_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "transparency", label: "Transparency", kind: .slider),
        RepProperty(setting: "mesh_width", label: "Width", kind: .slider, max: 2, step: 0.05),
      ]),
    "lines": RepSpec(
      rep: "lines", display: "Lines",
      colorSetting: "line_color", defaultColor: -1,
      properties: [
        RepProperty(
          setting: "line_width", label: "Width", kind: .slider, min: 0.5, max: 10, step: 0.5,
          decimals: 1)
      ]),
    "dots": RepSpec(
      rep: "dots", display: "Dots",
      colorSetting: "dot_color", defaultColor: -1,
      properties: [
        RepProperty(setting: "transparency", label: "Transparency", kind: .slider),
        RepProperty(
          setting: "dot_density", label: "Density", kind: .segmented,
          options: [("0", 0), ("1", 1), ("2", 2), ("3", 3)]),
        RepProperty(setting: "dot_radius", label: "Radius", kind: .slider, max: 1, step: 0.05),
      ]),
    "nonbonded": RepSpec(
      rep: "nonbonded", display: "Nonbonded",
      colorSetting: "", defaultColor: -1,
      properties: [
        RepProperty(setting: "nonbonded_size", label: "Size", kind: .slider, max: 1, step: 0.05)
      ]),
    "nb_spheres": RepSpec(
      rep: "nb_spheres", display: "NB spheres",
      colorSetting: "", defaultColor: -1,
      properties: [
        RepProperty(setting: "nb_spheres_size", label: "Size", kind: .slider, max: 1, step: 0.05)
      ]),
    "labels": RepSpec(
      rep: "labels", display: "Labels",
      colorSetting: "label_color", defaultColor: -6,
      properties: [
        RepProperty(
          setting: "label_size", label: "Size", kind: .slider, min: 5, max: 40, step: 1, decimals: 0
        )
      ]),
  ]

  static func spec(_ rep: String) -> RepSpec? { specs[rep] }
  static func display(_ rep: String) -> String { specs[rep]?.display ?? rep }
}

// MARK: - Scene (global) parameter table

struct SceneParam: Identifiable {
  var id: String { setting }
  let setting: String
  let label: String
  let kind: RepControlKind
  var min: Double = 0
  var max: Double = 1
  var step: Double = 1
  var decimals: Int = 0
  var options: [(label: String, value: Double)] = []
  let group: String
  // When set, this control is a dependent sub-setting: shown indented under
  // its parent toggle and hidden while the parent (a setting key) is off.
  var dependsOn: String? = nil
  // The Outline color row binds to sceneState.outlineColor via a ColorPicker
  // instead of the kind switch.
  var isColor: Bool = false
  // One-line description shown in a (?) popover next to the control.
  var help: String = ""
}

enum SceneCatalog {
  // Ordered sub-groups of the Viewer Settings tab. Background, depth cue and the
  // ray lighting live in Render Options.
  static let groups = [
    "Canvas", "Camera", "Lighting", "Shadows & AO", "Metal optimization", "Effects", "Quality",
  ]

  static let params: [SceneParam] = [
    // --- Canvas: multi-object/state layout ---
    // grid_mode is an int (0=off, 1=by object, 2=by state); the toggle maps
    // off→0 / on→1 and reads on for any non-zero mode.
    SceneParam(
      setting: "grid_mode", label: "Grid", kind: .toggle, group: "Canvas",
      help: "Lay each object out in its own grid cell instead of overlaid in one view."),
    SceneParam(
      setting: "all_states", label: "Overlay all states", kind: .toggle, group: "Canvas",
      help: "Show every coordinate state (NMR models / trajectory frames) at once."),
    // Not a PyMOL setting — bound to engine.hoverPreviewEnabled (see SceneParamRow's
    // special-case) and persisted across launches.
    SceneParam(
      setting: "hover_indicator", label: "Hover indicator", kind: .toggle, group: "Canvas",
      help:
        "Highlight the residue under the cursor (nearest Cα) before you click, in a distinct cyan — without changing the selection. Persists across launches."
    ),

    // --- Camera: viewpoint + lens ---
    SceneParam(
      setting: "field_of_view", label: "Lens (mm)", kind: .slider, min: 12, max: 135, step: 0.5,
      decimals: 0, group: "Camera",
      help:
        "Focal length (35mm-equivalent). Like swapping physical lenses: the camera dollies to keep the subject framed, so a short/wide lens (~12mm) exaggerates perspective (fisheye) and a long lens (~135mm) flattens it (macro/telephoto). Drag for a live, continuous perspective change. Pinch to zoom. No effect in orthoscopic mode."
    ),
    SceneParam(
      setting: "zoom", label: "Zoom (×)", kind: .slider, min: 0.5, max: 8, step: 0.1, decimals: 1,
      group: "Camera",
      help:
        "Apparent magnification — dollies the camera closer/farther. ~1× fits the whole scene; higher zooms in. Independent of the Lens: changing perspective (dolly-zoom) leaves this put, and vice-versa. Pinch to zoom also updates it."
    ),
    SceneParam(
      setting: "ortho", label: "Orthographic", kind: .toggle, group: "Camera",
      help: "Orthographic (parallel) projection — no perspective. Disables the Lens control."),

    // --- Lighting: real-time shading ---
    SceneParam(
      setting: "shininess", label: "Shininess", kind: .slider, min: 0, max: 100, step: 1,
      decimals: 0, group: "Lighting",
      help: "Tightness of specular highlights. Higher = smaller, sharper glints."),

    // --- Shadows & AO: screen-space occlusion ---
    SceneParam(
      setting: "metal_ssao", label: "Ambient occlusion", kind: .toggle, group: "Shadows & AO",
      help: "Screen-space ambient occlusion — darkens crevices and contact points for depth."),

    // --- Metal optimization: GPU quality/perf knobs ---
    SceneParam(
      setting: "metal_msaa", label: "MSAA 4×", kind: .toggle, group: "Metal optimization",
      help: "4× multisample antialiasing — smoother edges at some GPU cost."),

    // --- Effects: stylization + post-processing ---
    SceneParam(
      setting: "metal_outline", label: "Outline", kind: .toggle, group: "Effects",
      help: "Draw a silhouette / toon outline around objects."),
    SceneParam(
      setting: "metal_outline_color", label: "Outline color", kind: .toggle, group: "Effects",
      dependsOn: "metal_outline", isColor: true,
      help: "Color of the outline contour."),
    SceneParam(
      setting: "metal_outline_width", label: "Outline width", kind: .slider, min: 0.5, max: 5.0,
      step: 0.1, decimals: 1, group: "Effects", dependsOn: "metal_outline",
      help: "Thickness of the outline, in pixels."),

    // --- Quality: tessellation ---
    SceneParam(
      setting: "surface_quality", label: "Surface quality", kind: .segmented,
      options: [("0", 0), ("1", 1), ("2", 2)], group: "Quality",
      help: "Surface mesh detail: 0 = coarse/fast, 2 = fine/slow."),
  ]
}
