// PyMOLEngine+Picking.swift — click picking, the hover preview and the
// right-click hit.

import Foundation

extension PyMOLEngine {
  // Tap-to-select via metal_pick (NDC in [-1,1], aspect = width/height).
  func pick(ndcX: Float, ndcY: Float, aspect: Float) {
    #if DEBUG
      if let viewportInputTap {
        viewportInputTap(.pick(ndcX: ndcX, ndcY: ndcY))
        return
      }
    #endif
    guard let inst = instance else { return }
    PyMOLBridge_Pick(inst, ndcX, ndcY, aspect)
  }

  /// Middle click: center the view (and rotation origin) on the atom at (ndcX, ndcY).
  func centerAt(ndcX: Float, ndcY: Float, aspect: Float) {
    runPython("from mmol import pick as _mp\n_mp.center_at(\(ndcX), \(ndcY), \(aspect))")
  }

  // MARK: - Hover pre-selection preview (issue #165)

  /// Update the hover pre-selection preview to what a click at (ndcX,ndcY)
  /// would select, without committing to 'sele'. No-op while the feature is
  /// disabled. LEADING-EDGE throttled (+ trailing catch-up) so the preview
  /// tracks CONTINUOUSLY while the pointer sweeps — a pure trailing debounce
  /// only fired once the pointer paused, which read as "hover doesn't work".
  func hoverPreview(_ ndcX: Float, _ ndcY: Float, _ aspect: Float) {
    guard isReady, hoverPreviewEnabled else { return }
    if let last = lastHoverNDC {
      let dx = ndcX - last.0
      let dy = ndcY - last.1
      if abs(dx) < kHoverMinNDC && abs(dy) < kHoverMinNDC { return }
    }
    lastHoverNDC = (ndcX, ndcY)
    hoverThrottle.run { [weak self] in
      guard let self = self else { return }
      // A trailing pick scheduled before `ray` must not run over the held
      // ray image (desktop PyMOL has no hover pick).
      if self.rayImageHeld { return }
      self.runPython(
        "from mmol import pick as _mp; "
          + "_mp.hover_preview_at(\(ndcX), \(ndcY), \(aspect))")
    }
  }

  /// Cancel any pending hover pick and empty the preview selection. Cheap and
  /// idempotent — safe to call from mouse-down / drag-start / mouse-exit and
  /// when the user toggles the feature off.
  func clearHoverPreview() {
    hoverThrottle.cancel()
    lastHoverNDC = nil
    guard isReady else { return }
    runPython("from mmol import pick as _mp; _mp.clear_preview()")
  }

  /// Right-click: identify the atom/residue under the click (read-only — does
  /// NOT change the selection) and publish it so ContentView can show the
  /// context menu. pick_info_at writes the hit JSON to the temp dir, which we
  /// read back synchronously (this runs on the main thread from the event).
  func longPressPick(ndcX: Float, ndcY: Float, aspect: Float) {
    guard isReady else { return }
    runPython("from mmol import pick as _mp; _mp.pick_info_at(\(ndcX), \(ndcY), \(aspect))")
    let path = TempChannel.path(TempChannel.Stem.longPress)
    guard let data = FileManager.default.contents(atPath: path),
      let root = try? JSONSerialization.jsonObject(with: data) as? [String: Any]
    else { return }
    let hit = (root["hit"] as? Bool) ?? false
    longPressHit = LongPressHit(
      isEmpty: !hit,
      obj: root["obj"] as? String ?? "",
      chain: root["chain"] as? String ?? "",
      resi: root["resi"] as? String ?? "",
      resn: root["resn"] as? String ?? "",
      name: root["name"] as? String ?? "",
      sel: root["sel"] as? String ?? "")
  }
}
