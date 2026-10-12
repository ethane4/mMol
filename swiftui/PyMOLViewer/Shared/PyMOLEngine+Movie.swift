// PyMOLEngine+Movie.swift — movie playback, per-object model playback, scenes
// and the timeline track.

import Foundation

extension PyMOLEngine {
  // MARK: - Timeline / playback controls
  //
  // All routed through runPython (raw PyRun, no log echo) so high-rate scrub
  // calls don't flood the feedback log. Play/pause uses cmd.mplay/mstop and
  // lets the core's SceneIdle advance frames (the Metal draw loop ticks idle()
  // every frame); we never run a Swift frame-advance timer, which would race
  // the core and double-advance.

  private func movieCmd(_ call: String) {
    runPython("from pymol import cmd as _m\n_m.\(call)")
  }

  func play() {
    movieCmd("mplay()")
    playback.isPlaying = true
  }
  func pause() {
    movieCmd("mstop()")
    playback.isPlaying = false
  }
  func togglePlay() { playback.isPlaying ? pause() : play() }

  // These set the movie frame while paused; force one repaint so the viewport
  // leaves the on-demand gate and shows the new frame (issue #132) — mirrors
  // the per-object `set state` case in runCommandCore.
  func rewindMovie() {
    movieCmd("rewind()")
    requestViewportRedraw()
  }
  func endingMovie() {
    movieCmd("ending()")
    requestViewportRedraw()
  }
  func stepForward() {
    movieCmd("forward()")
    requestViewportRedraw()
  }
  func stepBackward() {
    movieCmd("backward()")
    requestViewportRedraw()
  }

  // MARK: - Per-object model (state) playback — NMR/MD inspection, independent
  // of the movie. Each object animates its OWN coordinate state at its OWN fps
  // via a Swift timer (`set state, k, obj`), so ensembles can play at different
  // rates without touching the global movie frame. The state-row slider follows
  // via the object poll.

  func objectPlaybackFPS(_ name: String) -> Double { objectFPS[name] ?? 15 }

  func setObjectFPS(_ name: String, _ fps: Double) {
    objectFPS[name] = fps
    if playingObjects.contains(name) { startObjectStates(name) }  // re-time
  }

  func toggleObjectStates(_ name: String) {
    playingObjects.contains(name) ? stopObjectStates(name) : startObjectStates(name)
  }

  func startObjectStates(_ name: String) {
    stopObjectStates(name)
    guard let total = objects.first(where: { $0.name == name })?.stateCount, total > 1 else {
      return
    }
    let fps = objectPlaybackFPS(name)
    var k = min(max(objectMeta[name]?.state ?? 1, 1), total)
    playingObjects.insert(name)
    let t = Timer(timeInterval: 1.0 / max(fps, 0.1), repeats: true) { [weak self] _ in
      guard let self else { return }
      k = k >= total ? 1 : k + 1
      self.runCommand("set state, \(k), \(name)", naming: name)  // #531
    }
    RunLoop.main.add(t, forMode: .common)  // keeps ticking during scroll/interaction
    objectStateTimers[name] = t
  }

  func stopObjectStates(_ name: String) {
    objectStateTimers[name]?.invalidate()
    objectStateTimers[name] = nil
    playingObjects.remove(name)
  }

  // Live scrub: clamp, set immediately for snappy UI, throttle the core call.
  func scrub(to frame: Int) {
    let f = max(1, min(frame, max(playback.frameCount, 1)))
    isScrubbing = true
    playback.currentFrame = f
    scrubReleaseWork?.cancel()
    guard f != lastScrubFrame else { return }
    lastScrubFrame = f
    movieCmd("frame(\(f))")
    // Paused-scrub redraw: without an active mplay the on-demand gate can
    // swallow the frame change once a movie exists, so nudge one repaint
    // (issue #132). Harmless while playing (already redrawing every tick).
    requestViewportRedraw()
  }

  // Drag ended: commit the final frame and release the scrub lock after a
  // beat so the ~poll can resume mirroring core state without yanking.
  func endScrub() {
    let work = DispatchWorkItem { [weak self] in
      self?.isScrubbing = false
      self?.lastScrubFrame = -1
    }
    scrubReleaseWork = work
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.25, execute: work)
  }

  func setMovieFPS(_ fps: Double) {
    let f = max(0.1, fps)
    playback.movieFPS = f
    movieCmd("set('movie_fps', \(f))")
  }

  func setMovieLoop(_ on: Bool) {
    playback.movieLoop = on
    movieCmd("set('movie_loop', \(on ? 1 : 0))")
  }

  func setShowFrameRate(_ on: Bool) {
    movieCmd("set('show_frame_rate', \(on ? 1 : 0))")
  }

  // Recall a saved scene (camera + reps + visibility).
  func recallScene(_ name: String) { sceneAction(name, "recall") }

  // Per-scene management for the scene chips' long-press menu. Injection-safe.
  // updateScene overwrites the named scene with the CURRENT view/reps.
  func updateScene(_ name: String) { sceneAction(name, "update") }
  func deleteScene(_ name: String) { sceneAction(name, "delete") }

  private func sceneAction(_ name: String, _ action: String) {
    guard !name.isEmpty else { return }
    runPython("from pymol import cmd as _c\n_c.scene(\(pyLiteral(name)), '\(action)')")
  }

  func renameScene(_ name: String, to newName: String) {
    let n = newName.trimmingCharacters(in: .whitespacesAndNewlines)
    guard !name.isEmpty, !n.isEmpty, n != name else { return }
    runPython(
      "from pymol import cmd as _c\n"
        + "_c.scene(\(pyLiteral(name)), 'rename', new_key=\(pyLiteral(n)))")
  }

  // Move the playhead to `frame` and commit (for tapping a keyframe/ruler in
  // the timeline). Wraps the scrub + release pair used by the transport.
  func seek(to frame: Int) {
    scrub(to: frame)
    endScrub()
  }

  // MARK: - Timeline drag-and-drop (re-time keyframes, place/move scenes)

  // MARK: - Unified timeline track ops
  //
  // The docked TimelinePanel's single lane. Every mutation ends in
  // rebuildMovie(), which recomputes frames from the transition durations and
  // rebuilds the entire core movie from `timelineItems` (Swift is the source
  // of truth). Camera views live in appkit_movie._views, keyed by item UUID.

  // Cumulative 1-based frame of each item, derived from the transition
  // durations: item[0] sits at frame 1 (its inbound transition is ignored);
  // each later item is `transition.seconds` (× fps) past the previous, forced
  // strictly increasing so mview keyframes never collide. [] when empty.
  // Each item's (startFrame, endFrame). Camera/scene are zero-width points
  // (end == start); a states clip occupies its own span (end = start + its
  // duration × fps), so later items begin after it.
  func itemSpans() -> [(start: Int, end: Int)] {
    guard !timelineItems.isEmpty else { return [] }
    let fps = max(playback.movieFPS, 1)
    var out = [(start: Int, end: Int)](repeating: (1, 1), count: timelineItems.count)
    var acc = 0.0
    var lastSeqEnd = 0
    var seenSeq = false
    for (i, it) in timelineItems.enumerated() {
      // Positioned items (a camera dropped at the playhead) sit at their
      // absolute frame and don't participate in the sequential layout.
      if let af = it.atFrame {
        out[i] = (max(1, af), max(1, af))
        continue
      }
      let start: Int
      if !seenSeq {
        start = 1
        seenSeq = true
      } else {
        acc += max(it.transition.seconds, 0.1)
        start = max(lastSeqEnd + 1, 1 + Int((acc * fps).rounded()))
      }
      var end = start
      if case .states(let spec) = it.kind {
        end = start + max(1, Int((spec.durationSeconds * fps).rounded()))
      }
      out[i] = (start, end)
      lastSeqEnd = end
    }
    return out
  }

  func itemFrames() -> [Int] { itemSpans().map { $0.start } }

  var timelineTotalFrames: Int { itemSpans().map { $0.end }.max() ?? 1 }

  // Rebuild the whole core movie from `timelineItems`. Serializes an ordered
  // spec (camera → UUID into _views; scene → base64 name, which
  // appkit_movie.rebuild decodes) with each item's inbound-transition easing.
  // Mirrors playback.frameCount immediately (the ~10/s poll would otherwise lag
  // the ruler a beat).
  func rebuildMovie() {
    guard isReady else { return }
    if timelineItems.isEmpty {
      // "[]" → clears the movie (see appkit_movie.rebuild).
      runPython("from mmol import movie as _am\n_am.rebuild('[]')")
      playback.frameCount = 1
      playback.currentFrame = 1
      return
    }
    let spans = itemSpans()
    let parts: [[String: Any]] = timelineItems.enumerated().map { i, it in
      let (start, end) = spans[i]
      var d: [String: Any] = ["frame": start]
      switch it.kind {
      case .camera:
        d["cam"] = it.id.uuidString
      case .scene(let name):
        d["scene"] = Data(name.utf8).base64EncodedString()
      case .states(let spec):
        d["end"] = end
        d["states"] = 1
        d["objects"] = spec.objects ?? NSNull()
        d["mode"] = spec.mode.rawValue
        d["first"] = spec.firstModel
        d["last"] = spec.lastModel
        return d
      }
      d["power"] = it.transition.linear ? 1.0 : 0.0
      d["linear"] = it.transition.linear ? 1 : 0
      if let ps = it.pinnedState { d["state"] = ps }
      return d
    }
    guard let data = try? JSONSerialization.data(withJSONObject: parts),
      let json = String(data: data, encoding: .utf8)
    else { return }
    runPython("from mmol import movie as _am\n_am.rebuild(\(pyLiteral(json)))")
    playback.frameCount = timelineTotalFrames
    playback.currentFrame = 1
  }

  // Multi-state objects currently loaded (name, state count), and the max — used
  // by the timeline to represent a bare ensemble and to target "Play models".
  func multiStateObjects() -> [(name: String, count: Int)] {
    objects.filter { $0.stateCount > 1 }.map { ($0.name, $0.stateCount) }
  }
  var maxStateCount: Int { objects.map { $0.stateCount }.max() ?? 1 }

  // Append a "Play models" states clip (sweeps a model range over `seconds`).
  func appendStatesClip(
    objects: [String]? = nil, mode: StatesMode = .sweep,
    firstModel: Int = 1, lastModel: Int = 0, seconds: Double = 4.0
  ) {
    guard isReady else { return }
    let spec = StatesSpec(
      objects: objects, mode: mode,
      firstModel: firstModel, lastModel: lastModel, durationSeconds: seconds)
    timelineItems.append(TimelineItem(kind: .states(spec)))
    rebuildMovie()
  }

  // Re-position a frame-anchored camera item (dragged along the lane).
  func setItemAtFrame(_ id: UUID, _ frame: Int) {
    guard let i = timelineItems.firstIndex(where: { $0.id == id }),
      timelineItems[i].atFrame != nil
    else { return }
    timelineItems[i].atFrame = max(1, frame)
    rebuildMovie()
  }

  // True when the lane has a "Play models" clip (camera keyframes then anchor to
  // the playhead so they can be placed within the clip's span).
  var hasStatesClip: Bool {
    timelineItems.contains { if case .states = $0.kind { return true } else { return false } }
  }

  // Update a states clip's spec (mode / duration / objects) in place.
  func updateStatesClip(_ id: UUID, _ spec: StatesSpec) {
    guard let i = timelineItems.firstIndex(where: { $0.id == id }) else { return }
    timelineItems[i].kind = .states(spec)
    rebuildMovie()
  }

  // Append a camera keyframe of the CURRENT view to the end of the lane.
  func captureCameraItem() {
    guard isReady else { return }
    // With a "Play models" clip present, anchor the keyframe at the current
    // playhead so it can be dropped WITHIN the clip's span (camera moves while
    // the models cycle). Otherwise append sequentially as before.
    let at = hasStatesClip ? max(1, playback.currentFrame) : nil
    let item = TimelineItem(kind: .camera, atFrame: at)
    runPython("from mmol import movie as _am\n_am.capture_view('\(item.id.uuidString)')")
    timelineItems.append(item)
    rebuildMovie()
    seekToItem(item.id)
  }

  // Append a scene marker (from the palette) to the end of the lane. The
  // camera flies into the scene and its reps cut in at playback.
  func appendSceneItem(_ name: String) {
    guard isReady, !name.isEmpty else { return }
    let item = TimelineItem(kind: .scene(name: name))
    timelineItems.append(item)
    rebuildMovie()
    seekToItem(item.id)
  }

  // Reorder: move the item at `from` to final index `to` (drag past a
  // neighbor). Ripples all frames + one rebuild.
  func moveItem(from: Int, to: Int) {
    guard timelineItems.indices.contains(from) else { return }
    let clamped = min(max(to, 0), timelineItems.count - 1)
    guard from != clamped else { return }
    let it = timelineItems.remove(at: from)
    timelineItems.insert(it, at: clamped)
    rebuildMovie()
  }

  func deleteItem(_ id: UUID) {
    guard let idx = timelineItems.firstIndex(where: { $0.id == id }) else { return }
    let it = timelineItems.remove(at: idx)
    if case .camera = it.kind {
      runPython("from mmol import movie as _am\n_am.forget_view('\(id.uuidString)')")
    }
    rebuildMovie()
  }

  // Edit the transition INTO `id` (duration + easing) — the long-press config.
  func setTransition(_ id: UUID, seconds: Double, linear: Bool) {
    guard let idx = timelineItems.firstIndex(where: { $0.id == id }) else { return }
    timelineItems[idx].transition = Transition(seconds: seconds, linear: linear)
    rebuildMovie()
  }

  func seekToItem(_ id: UUID) {
    guard let idx = timelineItems.firstIndex(where: { $0.id == id }) else { return }
    let frames = itemFrames()
    if frames.indices.contains(idx) { seek(to: frames[idx]) }
  }

  // Append a camera Roll/Rock as a run of 4 camera items joined by equal
  // transitions that sum to `duration`. Python computes the waypoint views
  // (spin/rock around the current view) under Swift-generated UUIDs.
  func appendCameraTemplate(kind: String, duration: Double, axis: String, angle: Double) {
    guard isReady else { return }
    let count = 4
    let ids = (0..<count).map { _ in UUID() }
    let idsJSON = "[" + ids.map { "\"\($0.uuidString)\"" }.joined(separator: ",") + "]"
    runPython(
      "from mmol import movie as _am\n"
        + "_am.capture_template_views(\(pyLiteral(idsJSON)), "
        + "'\(kind)', axis='\(axis)', angle=\(angle))")
    let per = max(duration / Double(count - 1), 0.2)  // 3 gaps span the duration
    // A roll is one continuous spin: easing each 120° gap would slow the camera
    // to a near-stop at every waypoint (a 6x speed pulse). Rock keeps the ease,
    // since its waypoints are where the motion reverses.
    let linear = kind == "roll"
    for id in ids {
      timelineItems.append(
        TimelineItem(
          id: id, kind: .camera,
          transition: Transition(seconds: per, linear: linear)))
    }
    rebuildMovie()
  }

  // Append a scene marker per saved scene, each `secondsPerScene` apart.
  func appendScenesTemplate(secondsPerScene: Double, scenes: [String]? = nil) {
    guard isReady else { return }
    let names = scenes ?? sceneNames
    guard !names.isEmpty else { return }
    for nm in names {
      timelineItems.append(
        TimelineItem(
          kind: .scene(name: nm),
          transition: Transition(seconds: secondsPerScene, linear: false)))
    }
    rebuildMovie()
  }

  // Empty the unified lane and clear the core movie + stored views.
  func clearMovieItems() {
    guard isReady else { return }
    timelineItems.removeAll()
    runPython("from mmol import movie as _am\n_am.forget_all_views()\n_am.reset_movie()")
    playback.frameCount = 1
    playback.currentFrame = 1
  }
}
