// PyMOLEngine+ClassicRay.swift — the desktop-PyMOL `ray` image shown over the
// viewport (mmol/classic_ray.py).

import AppKit

extension PyMOLEngine {
  // MARK: - Ray image

  /// The image the last `ray` produced (mmol/classic_ray.py), shown over the
  /// viewport until the core drops it.
  struct ClassicRayImage: Equatable {
    let image: NSImage
    let pixelWidth: Int
    let pixelHeight: Int
    let background: RGBA
  }

  /// RAYIMG:<w>:<h>:<r>,<g>,<b>:<png path>
  func parseClassicRayImage(_ line: String) {
    let parts = line.split(separator: ":", maxSplits: 4, omittingEmptySubsequences: false)
    guard parts.count == 5, let w = Int(parts[1]), let h = Int(parts[2]) else { return }
    let rgb = parts[3].split(separator: ",").compactMap { Double($0) }
    let path = String(parts[4])
    guard rgb.count == 3, let img = NSImage(contentsOfFile: path) else { return }
    try? FileManager.default.removeItem(atPath: path)
    classicRayImage = ClassicRayImage(
      image: img, pixelWidth: w, pixelHeight: h,
      background: RGBA(rgb[0], rgb[1], rgb[2]))
  }

  /// True while a ray image is up — including the moment between the core
  /// getting it and the next feedback poll showing it (a fast `rayg`), when
  /// a hover pick would otherwise make the core drop it unseen.
  var rayImageHeld: Bool {
    if classicRayImage != nil { return true }
    guard let inst = instance else { return false }
    return PyMOLBridge_HasSceneImage(inst) != 0
  }

  /// A click in the viewport dismisses the ray image (desktop PyMOL); merely
  /// moving the pointer must not, so hover previews are suspended while one is up.
  func dismissClassicRay() {
    guard classicRayImage != nil, let inst = instance else { return }
    PyMOLBridge_DropSceneImage(inst)
    classicRayImage = nil
  }

  /// Write the held ray image (what the viewport shows) to `path` — the
  /// export a desktop PyMOL user expects after `ray`. False if none is held.
  @discardableResult
  func saveHeldRayImage(to path: String) -> Bool {
    guard classicRayImage != nil else { return false }
    try? FileManager.default.removeItem(atPath: path)
    runPython("from pymol import cmd as _c\n_c.png(\(pyLiteral(path)), prior=1, quiet=1)\n")
    return FileManager.default.fileExists(atPath: path)
  }

  /// Route `ray` to the desktop-PyMOL display path (mmol/classic_ray.py).
  func enableClassicRay() {
    runPythonQuiet("from mmol import classic_ray\nclassic_ray.enable(True)\n")
  }
}
