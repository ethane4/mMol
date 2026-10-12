// PyMOLEngine+Move.swift — Box Select and Move mode (the rigid-body gizmo).

import Foundation

extension PyMOLEngine {
  // MARK: - Box Select (rubber-band selection, #358)

  /// Enter / leave Box Select (Tools ▸ Box Select).
  func toggleBoxSelect() {
    setInteractionMode(interactionMode == .boxSelect ? .viewing : .boxSelect)
  }

  /// Begin a new box: snapshot the selection so that every re-commit of THIS
  /// box composes against what was there before it, not against its own
  /// previous commit. Called when a drag starts a fresh rectangle, not when
  /// one adjusts the existing one.
  func beginBoxSession() {
    boxSessionActive = true
    runPython("from mmol import pick as _mp\n_mp.box_begin('sele')")
  }

  /// Update the rectangle and commit it. There is no preview/accept split: the
  /// committed selection IS the feedback, which is what makes the tool usable
  /// without hunting for a confirm button.
  func setBoxRect(_ rect: BoxRect?) {
    boxRect = rect
    guard let rect = rect, !rect.isDegenerate else { return }
    commitBoxRect(rect)
  }

  private func commitBoxRect(_ rect: BoxRect) {
    let aspect = gizmoAspect
    boxThrottle.run { [weak self] in
      // mode='add': the tool only ever grows the selection it opened on.
      self?.runPython(
        "from mmol import pick as _mp\n"
          + "_mp.box_commit_ndc(\(rect.minX), \(rect.minY), \(rect.maxX), \(rect.maxY), "
          + "\(aspect), name='sele', mode='add')")
    }
  }

  /// Forget the rectangle and close the session. The last commit STAYS — the
  /// tool has no cancel, because with commit-on-drag there is no uncommitted
  /// state left to throw away.
  private func clearBoxState() {
    // Reset the throttle clock too, so the FIRST rectangle of the next box
    // always commits on the leading edge instead of waiting out a window
    // left over from the previous one.
    boxThrottle.cancel(resetClock: true)
    boxDrag = nil
    boxRect = nil
    guard boxSessionActive else { return }
    boxSessionActive = false
    runPython("from mmol import pick as _mp\n_mp.box_finish()")
  }

  /// Leave the tool. The selection the box made is already committed.
  func endBoxSelection() {
    setInteractionMode(.viewing)
  }

  // MARK: - Move mode (rigid-body object gizmo)

  /// Viewport aspect (width / height) for the gizmo projection. Falls back to
  /// 1.0 before the first reshape.
  var gizmoAspect: Float {
    let s = viewportPixelSize
    return s.height > 0 ? Float(s.width / s.height) : 1.0
  }

  func setInteractionMode(_ mode: InteractionMode) {
    interactionMode = mode
    // Leaving Box Select (for viewing OR for another tool) must drop the
    // rectangle and its preview highlight, or the cyan atoms stay lit with no
    // box on screen to explain them.
    if mode != .boxSelect { clearBoxState() }
    // The viewport tools (Move, Box Select) are exclusive with the other
    // modes; .viewing is not a tool and leaves them alone (Esc unwinds them
    // via their own setters).
    if mode == .move || mode == .boxSelect {
      if measureMode != nil { setMeasureMode(nil) }  // mutually exclusive
    }
    if mode == .boxSelect {
      // The box freezes hover picking (MetalViewport.handleMouseMoved), so
      // nothing will refresh or clear the hover highlight while it is up —
      // drop it now rather than leave whatever the pointer happened to be
      // over when the tool was switched on lit up.
      clearHoverPreview()
      // Open with a rectangle already down and already committed. The tool
      // has no other affordance now that the overlay bar is gone, so an
      // empty viewport would just look like nothing happened.
      beginBoxSession()
      setBoxRect(.initial)
    }
    if mode == .move {
      refreshGizmo()
    } else {
      hoveredHandle = nil
      activeMoveObject = nil
      gizmo = nil
      // Clear lastAdjustSent FIRST. adjustFrameToggle and moveShiftHeld
      // both carry `didSet { syncAdjustFrame() }`, and syncAdjustFrame
      // calls readGizmo(), which re-reads the gizmo temp file — still
      // "active" here, because cleanup() below hasn't run yet. That
      // resurrects the gizmo and activeMoveObject just cleared above.
      // Zeroing lastAdjustSent first makes syncAdjustFrame's
      // `want != lastAdjustSent` guard early-return. The skipped
      // set_adjust(0) round-trip is redundant: cleanup() tears the whole
      // gizmo down a few lines later.
      lastAdjustSent = false
      adjustFrameToggle = false
      moveShiftHeld = false
      runPython("from mmol import move as _mm\n_mm.cleanup()")
    }
  }

  /// Push adjust-frame mode to metal_move when it changes (overlay toggle or the
  /// Shift key). set_adjust rebuilds the gizmo (greys it / restores it).
  func syncAdjustFrame() {
    let want = adjustFrameActive
    guard want != lastAdjustSent else { return }
    lastAdjustSent = want
    runPython("from mmol import move as _mm\n_mm.set_adjust(\(want ? 1 : 0))")
    readGizmo()
  }

  /// Clear the active object's custom gizmo frame → snap back to the auto frame.
  /// Distinct from resetActiveMovePosition (which resets the object's position).
  func resetGizmoFrame() {
    runPython("from mmol import move as _mm\n_mm.reset_gizmo()")
    readGizmo()
  }

  /// Grab-what-you-touch: set the active object to whatever is under the point.
  func moveSetActiveAt(ndcX: Float, ndcY: Float, aspect: Float) {
    runPython("from mmol import move as _mm\n_mm.pick_object(\(ndcX), \(ndcY), \(aspect))")
    readGizmo()
  }

  /// Explicitly set the active object (overlay dropdown).
  func setActiveMoveObject(_ name: String) {
    runPython("from mmol import move as _mm\n_mm.set_active(\(pyLiteral(name)), \(gizmoAspect))")
    readGizmo()
  }

  /// Re-emit the gizmo geometry (after the camera changed).
  func refreshGizmo(aspect: Float? = nil) {
    guard interactionMode == .move else { return }
    runPython("from mmol import move as _mm\n_mm.refresh(\(aspect ?? gizmoAspect))")
    readGizmo()
  }

  func gizmoBeginDrag(_ handle: GizmoHandle, ndcX: Float, ndcY: Float, aspect: Float) {
    runPython(
      "from mmol import move as _mm\n_mm.begin_drag('\(handle.pyName)', \(ndcX), \(ndcY), \(aspect))"
    )
    readGizmo()
  }

  func gizmoUpdateDrag(ndcX: Float, ndcY: Float, aspect: Float) {
    // Mid-drag: only issue the move. Skip readGizmo() — metal_move deliberately
    // does NOT re-emit geometry per tick (the gizmo is a 3D CGO that tracks the
    // molecule via its synced TTT, and hit-testing is idle while a handle is
    // grabbed), so there is no file to read back. This keeps the drag as light
    // as a vanilla-GL matrix update; the continuous draw loop repaints from
    // PyMOL's redisplay flag. end_drag re-emits geometry for the next hit-test.
    runPython("from mmol import move as _mm\n_mm.update_drag(\(ndcX), \(ndcY), \(aspect))")
  }

  func gizmoEndDrag() {
    runPython("from mmol import move as _mm\n_mm.end_drag()")
    readGizmo()
  }

  func resetActiveMovePosition() {
    runPython("from mmol import move as _mm\n_mm.reset_active()")
    readGizmo()
  }

  // Read the gizmo geometry metal_move wrote to the temp dir (same-process
  // TMPDIR), synchronously — called right after each metal_move call (the
  // longPressPick pattern). Updates the published gizmo + active object.
  private func readGizmo() {
    let path = TempChannel.path(TempChannel.Stem.gizmo)
    guard let data = FileManager.default.contents(atPath: path),
      let root = try? JSONSerialization.jsonObject(with: data) as? [String: Any]
    else { return }
    if (root["active"] as? Bool) == true, let g = GizmoGeometry(json: root) {
      // Only republish when the projected geometry actually changed. This is
      // called on every hover step (to keep the hit-test current against any
      // view change), so a static view must not fire objectWillChange each
      // move — that would churn every view observing the engine.
      if gizmo != g { gizmo = g }
      if activeMoveObject != g.obj { activeMoveObject = g.obj }
    } else {
      if gizmo != nil { gizmo = nil }
      if activeMoveObject != nil { activeMoveObject = nil }
    }
  }
}
