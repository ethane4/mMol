// MetalViewport.swift — MTKView wrapper for SwiftUI (NSViewRepresentable).

import MetalKit
import SwiftUI

// MARK: - Render-loop gate (#396)

/// The decision the render loop makes on every display-link tick, factored out
/// as a pure type so it can be tested (a test cannot drive a CADisplayLink).
///
/// Two things had to change here. The old gate consumed PyMOL's redisplay flag
/// with `GetRedisplay(inst, reset: 1)` *before* the drawable was acquired, so a
/// frame that then failed to get a drawable dropped both the frame and the
/// request for it — the same shape as the earlier "mouse rotate silently
/// dropped" bug. And when the GPU is the bottleneck, the loop would keep
/// handing frames to a queue that could not retire them, so the main thread sat
/// in `currentDrawable` (blocking input, the feedback timer and SwiftUI commits)
/// for about half of every rotation.
///
/// So: the gate PEEKS the flag (`reset: 0`), a separate in-flight cap makes the
/// loop skip a tick rather than pile on, and the flag is consumed only once the
/// frame is actually about to be encoded.
enum RenderGate {

  /// Committed-but-not-completed frames allowed before a tick is skipped.
  ///
  /// `CAMetalLayer.maximumDrawableCount` is 3, so staying below it leaves a
  /// drawable free and keeps the (now late) `currentDrawable` wait short.
  /// Two is still enough in-flight work to keep the GPU fed.
  static let maxFramesInFlight = 2

  enum Decision: Equatable {
    /// Encode a frame. The caller consumes the redisplay flag first.
    case render
    /// Nothing has changed; leave the last presented frame on screen.
    case skip
    /// There IS work to do, but too many frames are already in flight.
    /// The redisplay flag is left set, so a later tick picks this up.
    case throttle
  }

  static func decide(
    forceRedraw: Bool,
    hasRenderedOnce: Bool,
    redisplayPending: Bool,
    framesInFlight: Int
  ) -> Decision {
    // The very first frame always renders (defensive against a blank start
    // before anything has flagged a redisplay), and so does a forced one
    // after a wake/activate — the display can discard the drawable's
    // contents while asleep, and an unchanged scene flags no redisplay.
    let wants = forceRedraw || !hasRenderedOnce || redisplayPending
    guard wants else { return .skip }
    return framesInFlight >= maxFramesInFlight ? .throttle : .render
  }

  /// Whether the next tick must render unconditionally, given whether the
  /// frame just encoded reached the screen. A frame with no drawable rendered
  /// into the offscreen targets but presented nothing, and its redisplay flag
  /// is already consumed — without this the viewport would hold a stale image
  /// until the next interaction.
  static func forceRedrawAfterRender(presented: Bool) -> Bool { !presented }
}

struct MetalViewport: NSViewRepresentable {
  @EnvironmentObject var engine: PyMOLEngine

  func makeNSView(context: Context) -> MTKView {
    let view = PyMOLMTKView(frame: .zero)
    view.device = MTLCreateSystemDefaultDevice()
    view.delegate = context.coordinator
    view.colorPixelFormat = .bgra8Unorm
    view.depthStencilPixelFormat = .depth32Float_stencil8
    // Allow ProMotion (120Hz) on capable displays; the system clamps this to
    // the panel's actual max (e.g. 60 on non-ProMotion). The on-demand gate in
    // draw(in:) keeps the GPU idle on a static scene, so the higher tick only
    // costs a cheap idle poll when nothing is moving. The rate actually in
    // force comes from the view's own display link (viewDidMoveToWindow).
    view.preferredFramesPerSecond = 120
    view.enableSetNeedsDisplay = false
    view.isPaused = false
    context.coordinator.engine = engine
    context.coordinator.mtkView = view
    // Back-reference so the view's NSEvent overrides can reach the
    // coordinator's input handlers. Without this, mouseDown/Dragged/etc.
    // call `coordinator?.handle...` on a nil coordinator and silently
    // no-op — mouse rotate/zoom/pan never reach PyMOL.
    view.coordinator = context.coordinator

    // Repaint when the app re-activates or the system wakes from sleep: the
    // display can discard the drawable's contents while asleep/locked, and
    // the on-demand gate (a static scene flags no redisplay) would otherwise
    // leave the viewport black until the next interaction.
    NotificationCenter.default.addObserver(
      context.coordinator, selector: #selector(Coordinator.handleWake),
      name: NSApplication.didBecomeActiveNotification, object: nil)
    NSWorkspace.shared.notificationCenter.addObserver(
      context.coordinator, selector: #selector(Coordinator.handleWake),
      name: NSWorkspace.didWakeNotification, object: nil)
    // A state/frame change while a movie exists can leave the viewport on its
    // on-demand gate with nothing to repaint it; the engine posts this to
    // force one unconditional frame so the new state shows (issue #132).
    NotificationCenter.default.addObserver(
      context.coordinator, selector: #selector(Coordinator.handleWake),
      name: PyMOLEngine.forceRedrawNotification, object: nil)

    // Trackpad pinch → zoom. Two-finger drag (scrollWheel) → translate;
    // see handleScrollWheel. A real mouse wheel still zooms.
    let magnify = NSMagnificationGestureRecognizer(
      target: context.coordinator,
      action: #selector(Coordinator.handleMagnification(_:)))
    view.addGestureRecognizer(magnify)

    // Trackpad two-finger twist → Z-axis roll. Only acts while the engine's
    // "Trackpad" mouse mode is on (engine.trackpadMode). Otherwise it's a
    // no-op so it never fights the classic per-button mouse modes.
    let rotate = NSRotationGestureRecognizer(
      target: context.coordinator,
      action: #selector(Coordinator.handleRotationGesture(_:)))
    view.addGestureRecognizer(rotate)

    // Click-debug harness (PYMOL_AUTOCLICK="ndcx,ndcy[;ndcx,ndcy...]"): after
    // the scene renders, synthesize real clicks at the given NDC points
    // through the genuine mouse path. Each click's mouse→NDC math and the
    // resulting pick land in PYMOL_PICKDEBUG.
    if let spec = ProcessInfo.processInfo.environment["PYMOL_AUTOCLICK"] {
      let pts: [(CGFloat, CGFloat)] = spec.split(separator: ";").compactMap {
        let c = $0.split(separator: ",").compactMap { Double($0) }
        return c.count == 2 ? (CGFloat(c[0]), CGFloat(c[1])) : nil
      }
      for (i, p) in pts.enumerated() {
        DispatchQueue.main.asyncAfter(deadline: .now() + 3.0 + Double(i) * 1.0) {
          [weak coordinator = context.coordinator] in
          coordinator?.debugClick(ndcX: p.0, ndcY: p.1, in: view)
        }
      }
    }

    return view
  }

  func updateNSView(_ nsView: MTKView, context: Context) {}

  func makeCoordinator() -> Coordinator { Coordinator() }
}

// Custom MTKView subclass to handle mouse/keyboard events on macOS
class PyMOLMTKView: MTKView {
  weak var coordinator: MetalViewport.Coordinator?

  // Decline keyboard first-responder so a click in the viewport does NOT steal
  // focus from the command-line input (issue #73): the command line stays "hot"
  // for typing while the user rotates/picks, matching desktop PyMOL. Mouse events
  // are still delivered to this view via the mouse overrides below + acceptsFirstMouse
  // (they don't require first-responder status).
  // Keyboard shortcuts do NOT depend on this view's responder status: cmd.set_key
  // bindings are dispatched from an app-level NSEvent monitor in ContentView
  // (installPyMOLKeyMonitor, #258), which sees keys regardless of focus.
  override var acceptsFirstResponder: Bool { false }
  override func acceptsFirstMouse(for event: NSEvent?) -> Bool { true }

  // MARK: - Render loop
  //
  // Drive draw(in:) from a display link we own instead of MetalKit's.
  //
  // MTKView creates its own CVDisplayLink in viewDidMoveToWindow, but only
  // starts it for a window that is already visible — and it never retries. A
  // launch that opens a document (Finder double-click on a .pse, `open file.pse`)
  // attaches this view to a window that is not visible yet, so the loop never
  // starts: draw(in:) is then reached ONLY by the handful of manual draw() calls
  // in PyMOLEngine.initialize()'s blank-on-launch guard and in handleWake. The
  // window paints those few frames and is stale from then on. Because PyMOL
  // executes queued mouse input inside a rendered frame (OrthoExecDeferred runs
  // from SceneRenderMetal), every drag is silently swallowed too — the viewport
  // looks frozen while menus, panels and the console stay responsive. Launching
  // the app first and opening the file afterwards hid it, since that view was
  // attached to a window that did become visible in time.
  //
  // Toggling isPaused does NOT fix it: MetalKit will not rebuild a display link
  // it never created. Owning the link sidesteps the ordering entirely. It ticks
  // at the screen's refresh rate and follows the view across displays; the
  // on-demand gate in draw(in:) still skips the GPU work for a static scene, so
  // a still viewport costs the same idle poll it did before.
  private var renderLink: CADisplayLink?

  override func viewDidMoveToWindow() {
    super.viewDidMoveToWindow()
    guard window != nil else {
      // No window: drop our link and hand the loop back to MetalKit, so a
      // view that is re-attached later is never left with nothing driving it.
      renderLink?.invalidate()
      renderLink = nil
      isPaused = false
      return
    }
    guard renderLink == nil else { return }
    let link = displayLink(target: self, selector: #selector(renderTick))
    link.add(to: .main, forMode: .common)
    // A range, not a fixed rate: the system is free to tick slower when the
    // display or thermal state calls for it. The 30 Hz floor keeps a scene
    // that IS animating from being throttled to a crawl.
    link.preferredFrameRateRange = CAFrameRateRange(minimum: 30, maximum: 120, preferred: 120)
    renderLink = link
    isPaused = true  // ours drives now; don't let both loops draw
  }

  @objc private func renderTick() { draw() }

  // Track pointer motion over the viewport so the hover pre-selection preview
  // (issue #165) can update as the mouse moves WITHOUT any button held. A
  // tracking area is required for mouseMoved/mouseExited to fire; recreate it
  // on every layout change so it always spans the current visible bounds.
  private var hoverTrackingArea: NSTrackingArea?
  override func updateTrackingAreas() {
    super.updateTrackingAreas()
    if let existing = hoverTrackingArea {
      removeTrackingArea(existing)
    }
    let area = NSTrackingArea(
      rect: .zero,
      options: [
        .mouseMoved, .mouseEnteredAndExited, .activeInKeyWindow,
        .inVisibleRect,
      ],
      owner: self, userInfo: nil)
    addTrackingArea(area)
    hoverTrackingArea = area
  }

  override func mouseMoved(with event: NSEvent) {
    coordinator?.handleMouseMoved(event, in: self)
  }
  override func mouseExited(with event: NSEvent) {
    coordinator?.handleMouseExited(event, in: self)
  }

  override func mouseDown(with event: NSEvent) {
    coordinator?.dismissClassicRay()
    coordinator?.handleMouseDown(event, in: self)
  }
  override func mouseUp(with event: NSEvent) {
    coordinator?.handleMouseUp(event, in: self)
  }
  override func mouseDragged(with event: NSEvent) {
    coordinator?.handleMouseDragged(event, in: self)
  }
  override func rightMouseDown(with event: NSEvent) {
    coordinator?.dismissClassicRay()
    coordinator?.handleRightMouseDown(event, in: self)
  }
  override func rightMouseUp(with event: NSEvent) {
    coordinator?.handleRightMouseUp(event, in: self)
  }
  override func rightMouseDragged(with event: NSEvent) {
    coordinator?.handleRightMouseDragged(event, in: self)
  }
  override func otherMouseDown(with event: NSEvent) {
    coordinator?.dismissClassicRay()
    coordinator?.handleOtherMouseDown(event, in: self)
  }
  override func otherMouseUp(with event: NSEvent) {
    coordinator?.handleOtherMouseUp(event, in: self)
  }
  override func otherMouseDragged(with event: NSEvent) {
    coordinator?.handleOtherMouseDragged(event, in: self)
  }
  override func scrollWheel(with event: NSEvent) {
    coordinator?.handleScrollWheel(event, in: self)
  }
}

/// PYMOL_FRAMELOG=<path>: per rendered frame, "interval idle render" in ms
/// (main-thread CPU time; profiling aid).
enum FrameLog {
  private static let handle: FileHandle? = ProcessInfo.processInfo.environment["PYMOL_FRAMELOG"]
    .flatMap {
      FileManager.default.createFile(atPath: $0, contents: nil)
      return FileHandle(forWritingAtPath: $0)
    }
  private static var last = CACurrentMediaTime()
  /// PYMOL_HANGLOG=<path>: main-thread stalls over 20 ms ("hang <ms>"), from a
  /// background ping every 5 ms.
  static func startHangLog() {
    guard let path = ProcessInfo.processInfo.environment["PYMOL_HANGLOG"] else { return }
    FileManager.default.createFile(atPath: path, contents: nil)
    guard let h = FileHandle(forWritingAtPath: path) else { return }
    Thread.detachNewThread {
      while true {
        let sent = CACurrentMediaTime()
        let done = DispatchSemaphore(value: 0)
        DispatchQueue.main.async { done.signal() }
        done.wait()
        let ms = (CACurrentMediaTime() - sent) * 1000
        if ms > 20 { h.write(String(format: "hang %.0f\n", ms).data(using: .utf8)!) }
        Thread.sleep(forTimeInterval: 0.005)
      }
    }
  }
  static func note(idle: Double, render: Double) {
    guard let handle else { return }
    let now = CACurrentMediaTime()
    handle.write(
      String(format: "%.1f %.1f %.1f\n", (now - last) * 1000, idle * 1000, render * 1000).data(
        using: .utf8)!)
    last = now
  }
}

// MARK: - Shared Coordinator (MTKViewDelegate + input handling)

extension MetalViewport {
  class Coordinator: NSObject, MTKViewDelegate {
    weak var engine: PyMOLEngine?
    weak var mtkView: MTKView?
    private var viewportSize: CGSize = .zero
    // Set when the app/display wakes (unlock, system wake, re-activate). The
    // next draw(in:) then renders unconditionally, bypassing the on-demand
    // gate, to repaint a drawable whose contents were discarded during sleep.
    fileprivate var forceRedraw = false

    // Wake/activate -> force one repaint. Also kicks an immediate draw() in
    // case the display link hasn't resumed ticking yet.
    @objc fileprivate func handleWake() {
      forceRedraw = true
      DispatchQueue.main.async { [weak self] in self?.mtkView?.draw() }
    }
    deinit {
      NotificationCenter.default.removeObserver(self)
      NSWorkspace.shared.notificationCenter.removeObserver(self)
    }
    // Track the mouse-down point to distinguish a click (pick/select) from
    // a drag (rotate). Point space, view coordinates.
    private var mouseDownLoc: CGPoint = .zero
    private var didDrag = false
    // Move mode: distance (view points) the pointer must travel before a
    // press becomes a DRAG (orbit / gizmo-handle manipulation). Below it, the
    // press stays a candidate TAP → object select. Without this dead-zone a
    // sub-pixel jitter during a tap (macOS emits mouseDragged for <1pt moves)
    // armed didDrag and routed the release to the orbit path, so tap-to-select
    // only worked for a perfectly still click ("doesn't work consistently").
    private let moveDragSlop: CGFloat = 6
    // Move mode: the gizmo handle grabbed at mouse-down (nil = none → camera).
    private var moveHandle: GizmoHandle?
    // Right-button equivalents: a pure right-click raises the viewport
    // context menu (CPU pick → engine.longPressHit), while a right-DRAG
    // still drives PyMOL's clip (slab). The button-down is deferred to the
    // first drag so a bare click never enters clip mode.
    private var rightDownLoc: CGPoint = .zero
    private var rightDidDrag = false
    // Left/right chords (desktop PyMOL behaviour): a left press while the
    // right button is held is ignored for its whole press, and a right press
    // that overlapped a left one never opens the scene menu on release.
    private var leftSuppressed = false
    private var rightChorded = false
    /// Right pressed while left was held: the drag zooms (right) until the
    /// right button comes up, as in desktop PyMOL.
    private var rightOverLeft = false
    /// The left press was interrupted by a right chord: its release is not
    /// a click (no pick).
    private var leftChorded = false

    // Trackpad pinch (NSMagnificationGestureRecognizer) → zoom via an
    // explicit camera dolly (engine.zoomBy). We can't use the scroll-wheel
    // BUTTON path: PyMOL's default three_button_viewing binds the bare wheel
    // to 'slab' (clip), so it would change the slab, not zoom. magnification
    // is cumulative from gesture start; feed the per-callback delta as a
    // zoom fraction (spread = positive = zoom in).
    private var lastMag: CGFloat = 0
    private let kZoomGain: CGFloat = 1.0

    // Trackpad two-finger twist → Z-roll (Trackpad mode only). NSRotationGesture
    // .rotation is cumulative radians; feed the per-callback delta as `turn z`.
    // Negated so a clockwise twist rolls the molecule clockwise on screen.
    private var lastRoll: CGFloat = 0
    private let kRollSign: Float = -1

    // Trackpad two-finger drag (delivered as precise scrollWheel events) →
    // translate. Synthesized as a PyMOL middle-button drag: a MIDDLE-DOWN at
    // the start, drag events that follow an accumulated synthetic cursor, and
    // a MIDDLE-UP when the gesture (incl. momentum) ends. A real mouse wheel
    // (no precise deltas) still zooms.
    private var panActive = false
    private var panCursorX: Int32 = 0
    private var panCursorY: Int32 = 0
    private var panEndDebounce: DispatchWorkItem?
    // Sign so the molecule follows the fingers (grab-and-move). Tunable.
    // Y is negated: macOS scrollingDeltaY is opposite the on-screen pan we
    // want (verified — up/down was inverted before the flip).
    private let kPanSignX: CGFloat = 1
    private let kPanSignY: CGFloat = -1

    // Gesture mode latched at drag START (a mid-drag modifier change can't
    // switch it). Shift held → synthesize a Shift+Right-button drag, which
    // PyMOL's three_button_viewing binds to 'clip' (vertical = move the slab
    // through the scene, horizontal = slab thickness). Otherwise a Middle-drag
    // = translate. Clip is not grab-and-move, so its Y is NOT negated; flip
    // kClipSignY if the up/down direction feels inverted.
    private var gestureButton: Int32 = PYMOL_BUTTON_MIDDLE
    private var gestureMods: Int32 = 0
    private var gestureIsClip = false
    // Trackpad mode: a two-finger gesture picks its action from its first few
    // points of travel. Mostly horizontal moves the view (both axes, until the
    // fingers lift); mostly vertical changes the slab, as the mouse wheel does
    // in 3-button viewing.
    private var axisUndecided = false
    private var axisTravel = CGPoint.zero
    private var gestureIsSlab = false
    private var slabAccum: CGFloat = 0
    private let kAxisDecide: CGFloat = 4  // points of travel before deciding
    private let kSlabStep: CGFloat = 6  // points of vertical travel per wheel step
    private let kClipSignX: CGFloat = 1
    private let kClipSignY: CGFloat = 1

    // Hover pre-selection preview (issue #165): last view-point location fed
    // to a preview, used to skip re-picking when the pointer barely moved.
    // .zero is treated as "no prior move" (any first move re-picks).
    private var lastHoverLoc: CGPoint = .zero

    // MARK: - MTKViewDelegate

    func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
      viewportSize = size
      engine?.viewportPixelSize = size
      engine?.reshape(width: Int(size.width), height: Int(size.height))
    }

    private var wasSuppressed = false
    private var hasRenderedOnce = false
    private var moveSyncCounter = 0
    func draw(in view: MTKView) {
      guard let engine = engine, engine.isReady else { return }
      // A movie export renders frames off the main thread and owns the core
      // exclusively (it reshapes global state per frame). Skip the live
      // render meanwhile so we never race it; the exporter restores the
      // scene + clears this flag when done, and the next tick redraws. (#58 L-59)
      if engine.exportRenderActive { return }
      // Keep the Move gizmo's 2D hit-geometry (+ debug bullseye overlay)
      // continuously in sync with the camera — not just on mouse-move — so a
      // cursor parked after an orbit/zoom still sees current targets and the
      // next grab maps to the visible gizmo. Throttled to ~15 Hz; readGizmo
      // republishes only when the projection actually moved (no churn on a
      // static view). Skipped WHILE dragging a handle: metal_move omits the
      // per-tick emit there and the CGO already tracks via its synced TTT.
      if engine.interactionMode == .move {
        let draggingGizmo = moveHandle != nil
        if !draggingGizmo {
          moveSyncCounter += 1
          if moveSyncCounter >= 8 {
            moveSyncCounter = 0
            engine.refreshGizmo()
          }
        } else {
          moveSyncCounter = 0
        }
      }
      // Panel-resize drag: while suppressed, freeze the drawable size so the
      // renderer doesn't reallocate offscreen targets on every frame (choppy
      // + OOM). On resume, snap the drawable to the current bounds ONCE,
      // which fires drawableSizeWillChange → a single reshape.
      let suppress = engine.suppressDrawableResize
      if suppress != wasSuppressed {
        wasSuppressed = suppress
        view.autoResizeDrawable = !suppress
        if suppress, let layer = view.layer {
          // The frozen frame stays inside the viewport, fitted to it, over the
          // scene's background color.
          var bg: [Float] = [0, 0, 0]
          PyMOLBridge_GetBackground(&bg)
          layer.masksToBounds = true
          layer.contentsGravity = engine.frozenDrawableGravity
          layer.backgroundColor = CGColor(
            red: CGFloat(bg[0]), green: CGFloat(bg[1]), blue: CGFloat(bg[2]), alpha: 1)
        }
        if !suppress {
          view.layer?.contentsGravity = .resize
          let scale = view.window?.backingScaleFactor ?? 2.0
          let target = CGSize(
            width: view.bounds.width * scale,
            height: view.bounds.height * scale)
          if target.width > 0, target.height > 0, target != view.drawableSize {
            view.drawableSize = target  // → drawableSizeWillChange → one reshape
          }
        }
      }
      // Build RendererMetal on the first frame (bridge no-ops thereafter),
      // then run PyMOL's idle work (advances movies/animations and sets the
      // redisplay flag).
      engine.setupMetalRenderer(view: view)
      let t0 = CACurrentMediaTime()
      engine.idle()
      let tIdle = CACurrentMediaTime()

      // On-demand rendering: skip the GPU-expensive frame when nothing
      // needs redrawing — a static structure then costs only a cheap idle
      // poll instead of a full render every tick, the bulk of the
      // battery/thermal win. The last presented frame stays on screen.
      // Mirrors the legacy AppKit loop (main_appkit.mm).
      //
      // The flag is PEEKED here (reset: 0) and consumed below only once
      // this tick has committed to encoding a frame, so a throttled tick
      // leaves the request standing instead of swallowing it (#396).
      let pending = engine.instance.map { PyMOLBridge_GetRedisplay($0, 0) != 0 } ?? false
      switch RenderGate.decide(
        forceRedraw: forceRedraw,
        hasRenderedOnce: hasRenderedOnce,
        redisplayPending: pending,
        framesInFlight: engine.metalFramesInFlight)
      {
      case .skip:
        // Nothing is dirty, so no deferred rep build is waiting on a
        // frame. A synchronous heavy op (image export / Copy Image)
        // leaves the scene clean, and waiting for "build frames" that
        // never render held its overlay up until the 60 s backstop.
        engine.heavyRenderTick(presented: false)
        return
      case .throttle:
        return
      case .render:
        break
      }

      // Consume the flag NOW, before the scene traversal: PyMOL sets it
      // again from inside SceneRenderMetal when queued input (a drag, a
      // click) needs another frame, so clearing it afterwards would eat
      // the next frame of every drag.
      if let inst = engine.instance { _ = PyMOLBridge_GetRedisplay(inst, 1) }

      // The bridge acquires the drawable itself, after the scene is
      // encoded — see PyMOLBridge_RenderMetalFrame. `presented` is false
      // when no drawable arrived: the frame rendered offscreen but nothing
      // reached the screen, so force the next tick to render again.
      let size = view.drawableSize
      let presented = engine.renderMetalFrame(
        view: view,
        width: Int(size.width),
        height: Int(size.height))
      FrameLog.note(idle: tIdle - t0, render: CACurrentMediaTime() - tIdle)
      hasRenderedOnce = hasRenderedOnce || presented
      forceRedraw = RenderGate.forceRedrawAfterRender(presented: presented)
      // This frame built any deferred rep geometry (e.g. a surface mesh);
      // let the engine clear the "Calculating…" overlay once the build
      // frame(s) have completed.
      engine.heavyRenderTick(presented: presented)
    }

    // MARK: - Coordinate conversion

    private func pymolPoint(in view: MTKView, at point: CGPoint) -> (Int32, Int32) {
      let backing = view.convertToBacking(point)
      return (Int32(backing.x), Int32(backing.y))
    }

    private func pymolModifiers(_ flags: UInt) -> Int32 {
      var mods: Int32 = 0
      let nsFlags = NSEvent.ModifierFlags(rawValue: flags)
      if nsFlags.contains(.shift) { mods |= PYMOL_MOD_SHIFT }
      if nsFlags.contains(.control) { mods |= PYMOL_MOD_CTRL }
      if nsFlags.contains(.option) { mods |= PYMOL_MOD_ALT }
      return mods
    }

    // MARK: - Move-mode gizmo hit-testing (shared)

    // A view point -> (ndc_x, ndc_y, aspect) in PyMOL NDC (bottom-left, +y up).
    // View points are already bottom-left.
    private func gizmoNDC(in view: MTKView, at p: CGPoint) -> (Float, Float, Float)? {
      let w = view.bounds.width
      let h = view.bounds.height
      guard w > 0, h > 0 else { return nil }
      let nx = Float(p.x / w) * 2 - 1
      let ny = Float(p.y / h) * 2 - 1
      return (nx, ny, Float(w / h))
    }

    // MARK: - Box Select input (#358)

    // True while Box Select owns the viewport. Every camera gesture below
    // early-returns on this: the issue's "camera controls are disabled" is
    // what makes drawing the box possible at all, since otherwise the same
    // drag would orbit the scene out from under the rectangle.
    var boxSelectActive: Bool { engine?.interactionMode == .boxSelect }

    // Grab band around the box's edges/corners, in NDC-Y units (~1.75% of the
    // viewport height each side, so ~20 pt of target on a 600 pt viewport).
    private let boxGrabSlop: CGFloat = 0.035

    // Press: grab an edge / corner / the interior of the existing box, or —
    // if the press misses it — start drawing a new one.
    //
    // A NEW box opens a session (snapshotting the selection it will compose
    // against); adjusting the existing one keeps the session it already has,
    // so dragging a corner back in takes atoms away again instead of
    // ratcheting them on.
    private func boxBegin(in view: MTKView, at p: CGPoint) {
      guard let engine = engine, let (nx, ny, aspect) = gizmoNDC(in: view, at: p) else { return }
      let ndc = CGPoint(x: CGFloat(nx), y: CGFloat(ny))
      if let existing = engine.boxRect,
        let drag = existing.beginDrag(at: ndc, aspect: CGFloat(aspect), slop: boxGrabSlop)
      {
        engine.boxDrag = drag
      } else {
        engine.boxDrag = BoxRect.beginNewDrag(at: ndc)
        engine.beginBoxSession()
      }
    }

    private func boxUpdate(in view: MTKView, at p: CGPoint) {
      guard let engine = engine, let drag = engine.boxDrag,
        let (nx, ny, _) = gizmoNDC(in: view, at: p)
      else { return }
      engine.setBoxRect(drag.rect(at: CGPoint(x: CGFloat(nx), y: CGFloat(ny))))
    }

    // Release. The drag has been committing all along, so this just lands the
    // final rectangle. A press that never moved is a TAP: on the box it
    // leaves it alone, outside it dismisses the box — WITHOUT committing an
    // empty rectangle, so a stray click can't wipe the selection the box just
    // made.
    private func boxEnd(in view: MTKView, at p: CGPoint) {
      guard let engine = engine else { return }
      defer { engine.boxDrag = nil }
      guard let drag = engine.boxDrag,
        let (nx, ny, _) = gizmoNDC(in: view, at: p)
      else { return }
      let rect = drag.rect(at: CGPoint(x: CGFloat(nx), y: CGFloat(ny)))
      if rect.isDegenerate {
        engine.setBoxRect(drag.translates ? drag.origin : nil)
      } else {
        engine.setBoxRect(rect)
      }
    }

    // The gizmo handle under a view point (Move mode only), or nil.
    private func gizmoHit(in view: MTKView, at p: CGPoint) -> GizmoHandle? {
      guard let g = engine?.gizmo, let (nx, ny, aspect) = gizmoNDC(in: view, at: p) else {
        return nil
      }
      return g.hitTest(ndc: CGPoint(x: CGFloat(nx), y: CGFloat(ny)), aspect: CGFloat(aspect))
    }

    // MARK: - macOS mouse handling

    func handleMouseDown(_ event: NSEvent, in view: MTKView) {
      if NSEvent.pressedMouseButtons & 0b10 != 0 {  // right is held
        leftSuppressed = true
        rightChorded = true
        return
      }
      leftChorded = false
      rightOverLeft = false
      // Don't send PyMOL a button-down yet: a left-click in PyMOL's
      // viewing mode runs SceneClick/Release, whose GL pick is dead on
      // Metal and ends up CLEARING the active selection. We only want
      // PyMOL's mouse handling for an actual drag (rotate), so the
      // button-down is deferred to the first drag event (below). A pure
      // click selects via metal_pick in mouseUp instead.
      // A committing click starts: drop any hover preview so it can't
      // linger under (or fight) the committed selection.
      engine?.clearHoverPreview()
      lastHoverLoc = .zero
      mouseDownLoc = view.convert(event.locationInWindow, from: nil)
      didDrag = false
      // Box Select: the press starts (or grabs) the rectangle; no PyMOL
      // button event is ever sent, so the camera cannot move.
      if boxSelectActive {
        boxBegin(in: view, at: mouseDownLoc)
        return
      }
      // Move mode: remember whether the press landed on a gizmo handle, so
      // the drag manipulates the object; otherwise the drag orbits the camera.
      moveHandle = nil
      if engine?.interactionMode == .move {
        // Shift → adjust-frame mode for this drag (re-anchor the gizmo).
        engine?.moveShiftHeld = event.modifierFlags.contains(.shift)
        // Refresh the hit geometry in case the click arrived with no prior
        // hover move (e.g. right after an external view change), so the grab
        // maps to the currently-visible handle. No-op cost if already fresh.
        if let (_, _, a) = gizmoNDC(in: view, at: mouseDownLoc) {
          engine?.refreshGizmo(aspect: a)
        }
        moveHandle = gizmoHit(in: view, at: mouseDownLoc)
      }
    }

    // Pointer moved over the viewport with no button held → refresh the hover
    // pre-selection preview (issue #165). No-op while measuring or during a
    // drag (didDrag guards the button-held drag path). Computes NDC EXACTLY
    // like handleMouseUp — bottom-left origin, no Y flip — so the preview
    // aligns with what a click would select. A sub-pixel move gate skips the
    // re-pick when the pointer barely moved; the engine additionally
    // debounces the actual Python pick.
    /// Any button press in the viewport dismisses a held `ray` image (desktop
    /// PyMOL behaviour) before the press is handled as usual.
    func dismissClassicRay() { engine?.dismissClassicRay() }

    func handleMouseMoved(_ event: NSEvent, in view: MTKView) {
      guard !didDrag else { return }
      // A held `ray` image stays up while the pointer
      // merely passes over it, as in desktop PyMOL: the hover pick edits a
      // preview selection, which would make the core drop the image.
      guard engine?.rayImageHeld != true else { return }
      // No hover picking while the box tool is on: the box is already
      // driving the selection, and a hover pick would fight it.
      guard !boxSelectActive else { return }
      let loc = view.convert(event.locationInWindow, from: nil)
      if lastHoverLoc != .zero,
        hypot(loc.x - lastHoverLoc.x, loc.y - lastHoverLoc.y) < 2
      {
        return
      }
      lastHoverLoc = loc
      // Move mode: highlight the gizmo handle under the cursor so it's
      // obvious which axis/ring/center a drag will grab. Also reflect Shift
      // so the gizmo greys out (adjust-frame mode) as you hover with it held.
      if engine?.interactionMode == .move {
        engine?.moveShiftHeld = event.modifierFlags.contains(.shift)
        // Re-emit the hit-test geometry for the CURRENT view before testing.
        // The cached 2D projection goes stale on ANY view change — not just
        // the four gesture-end sites that call refreshGizmo, but also the
        // Viewer Settings Lens / Zoom, orient·reset·center, MCP/AI commands,
        // movie playback and window resizes. Hovering the visible ring then
        // hit the OLD cached ring (→ no highlight, wrong grab). Refreshing
        // here — gated by the 2px hover step, cheaper than the per-hover atom
        // pick already run in viewing mode, and republished only when the
        // projection actually moved — makes hover + the subsequent grab
        // always map to the handle the user sees. Uses the hover's own aspect
        // so the emitted geometry and the cursor share one aspect.
        if let (_, _, a) = gizmoNDC(in: view, at: loc) {
          engine?.refreshGizmo(aspect: a)
        }
        let hit = gizmoHit(in: view, at: loc)
        if engine?.hoveredHandle != hit { engine?.hoveredHandle = hit }
        return
      }
      guard engine?.measureMode == nil else { return }
      let w = view.bounds.width
      let h = view.bounds.height
      guard w > 0, h > 0 else { return }
      let ndcX = Float(loc.x / w) * 2 - 1
      let ndcY = Float(loc.y / h) * 2 - 1
      engine?.hoverPreview(ndcX, ndcY, Float(w / h))
    }

    // Pointer left the viewport → clear the preview so it doesn't linger.
    func handleMouseExited(_ event: NSEvent, in view: MTKView) {
      lastHoverLoc = .zero
      guard engine?.rayImageHeld != true else { return }  // see handleMouseMoved
      engine?.clearHoverPreview()
      if engine?.hoveredHandle != nil { engine?.hoveredHandle = nil }
      // Drop Shift adjust-mode when the pointer leaves, so the gizmo doesn't
      // stay greyed if Shift is released outside the viewport.
      if engine?.moveShiftHeld == true { engine?.moveShiftHeld = false }
    }

    func handleMouseUp(_ event: NSEvent, in view: MTKView) {
      if leftSuppressed {
        leftSuppressed = false
        return
      }
      if leftChorded {
        // Left released after a right chord: end any resumed rotate,
        // but never treat it as a click.
        leftChorded = false
        if rightOverLeft {
          // right still held: it keeps zooming via rightMouseDragged
        } else if didDrag {
          let pt = pymolPoint(in: view, at: view.convert(event.locationInWindow, from: nil))
          engine?.button(
            PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1,
            modifiers: pymolModifiers(event.modifierFlags.rawValue))
        }
        didDrag = false
        lastHoverLoc = .zero
        return
      }
      // Clear the drag flag on exit so passive hover (which is gated on
      // !didDrag) resumes immediately after a drag, not only after the next
      // mouse-down.
      defer { didDrag = false }
      let loc = view.convert(event.locationInWindow, from: nil)
      if boxSelectActive {
        boxEnd(in: view, at: loc)
        return
      }
      let mods = pymolModifiers(event.modifierFlags.rawValue)
      let moved = hypot(loc.x - mouseDownLoc.x, loc.y - mouseDownLoc.y)

      if engine?.interactionMode == .move {
        let aspect = Float(view.bounds.width / max(view.bounds.height, 1))
        if moveHandle != nil, didDrag {
          engine?.gizmoEndDrag()
          moveHandle = nil
          return
        }
        moveHandle = nil
        if didDrag {
          // Finish the camera orbit (empty-space drag).
          let pt = pymolPoint(in: view, at: loc)
          engine?.button(
            PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1, modifiers: mods)
          engine?.refreshGizmo(aspect: aspect)
          return
        }
        // Not a drag (the press stayed inside moveDragSlop) → treat as a
        // TAP and grab-what-you-touch. No separate distance gate here: the
        // dead-zone in handleMouseDragged already guarantees didDrag is only
        // set once the pointer really moved, so a shaky tap still selects.
        if let (nx, ny, a) = gizmoNDC(in: view, at: loc) {
          engine?.moveSetActiveAt(ndcX: nx, ndcY: ny, aspect: a)
        }
        return
      }

      if didDrag {
        // Finish the rotate drag.
        let pt = pymolPoint(in: view, at: loc)
        engine?.button(PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1, modifiers: mods)
        // Reset the drag latch so hover resumes after the rotate ends.
        // (didDrag is only cleared on mouseDown; without this, a drag
        // leaves it true and handleMouseMoved's `!didDrag` guard blocks
        // all hover until the next click — hover "stopped after rotating".)
        didDrag = false
        lastHoverLoc = .zero
        return
      }

      // Pure click (no drag, no PyMOL button events sent) → CPU pick.
      // NDC in view-point space, bottom-left origin (macOS views aren't
      // flipped, matching PyMOL's NDC) so no Y flip.
      if moved < 4 {
        let w = view.bounds.width
        let h = view.bounds.height
        if w > 0, h > 0 {
          let ndcX = Float(loc.x / w) * 2 - 1
          let ndcY = Float(loc.y / h) * 2 - 1
          // Pick-debug: record the click in top-down (SwiftUI) points so
          // the overlay crosshair lands exactly where the user clicked.
          if PyMOLEngine.debugPickEnabled {
            engine?.debugClickPoint = CGPoint(x: loc.x, y: h - loc.y)
          }
          Self.pickDbg(
            String(
              format:
                "mouseUp loc=(%.1f,%.1f) bounds=(%.1f,%.1f) backing=%.2f -> ndc=(%.4f,%.4f) aspect=%.4f",
              loc.x, loc.y, w, h, view.window?.backingScaleFactor ?? 0,
              ndcX, ndcY, Float(w / h)))
          if engine?.measureMode != nil {
            engine?.measurePick(ndcX: ndcX, ndcY: ndcY, aspect: Float(w / h))
          } else {
            engine?.pick(ndcX: ndcX, ndcY: ndcY, aspect: Float(w / h))
          }
        }
      }
    }

    // --- Click-debug harness (PYMOL_AUTOCLICK) ---
    // Append a line to PYMOL_PICKDEBUG so the mouse→NDC math is visible
    // alongside pick_at's projection (which logs to the same file).
    static func pickDbg(_ s: String) {
      guard let path = ProcessInfo.processInfo.environment["PYMOL_PICKDEBUG"] else { return }
      if let fh = FileHandle(forWritingAtPath: path)
        ?? {
          FileManager.default.createFile(atPath: path, contents: nil)
          return FileHandle(forWritingAtPath: path)
        }()
      {
        fh.seekToEndOfFile()
        fh.write((s + "\n").data(using: .utf8)!)
        try? fh.close()
      }
    }

    // Synthesize a real left-click at the given NDC by converting NDC → view
    // point → window point and dispatching genuine NSEvents through the
    // view's mouseDown/mouseUp — the EXACT path a user click takes. Lets the
    // debug harness click precise scene positions without Accessibility.
    func debugClick(ndcX: CGFloat, ndcY: CGFloat, in view: MTKView) {
      let w = view.bounds.width
      let h = view.bounds.height
      guard w > 0, h > 0, let win = view.window else { return }
      let vp = CGPoint(x: (ndcX + 1) / 2 * w, y: (ndcY + 1) / 2 * h)  // bottom-left
      let winPt = view.convert(vp, to: nil)
      let ts = ProcessInfo.processInfo.systemUptime
      let mk = { (type: NSEvent.EventType) -> NSEvent? in
        NSEvent.mouseEvent(
          with: type, location: winPt, modifierFlags: [],
          timestamp: ts, windowNumber: win.windowNumber,
          context: nil, eventNumber: 0, clickCount: 1, pressure: 1)
      }
      Self.pickDbg(
        String(
          format: "debugClick ndc=(%.4f,%.4f) -> vpoint=(%.1f,%.1f) winpoint=(%.1f,%.1f)",
          Float(ndcX), Float(ndcY), vp.x, vp.y, winPt.x, winPt.y))
      if let d = mk(.leftMouseDown) { view.mouseDown(with: d) }
      if let u = mk(.leftMouseUp) { view.mouseUp(with: u) }
    }

    func handleMouseDragged(_ event: NSEvent, in view: MTKView) {
      // Left pressed during a right-drag: AppKit now reports the motion
      // as LEFT drags, so keep feeding the right-button (zoom) drag —
      // desktop PyMOL keeps zooming until the right button comes up.
      if leftSuppressed || rightOverLeft {
        if NSEvent.pressedMouseButtons & 0b10 != 0 { handleRightMouseDragged(event, in: view) }
        return
      }
      let loc = view.convert(event.locationInWindow, from: nil)
      if boxSelectActive {
        didDrag = true
        boxUpdate(in: view, at: loc)
        return
      }
      let mods = pymolModifiers(event.modifierFlags.rawValue)

      if engine?.interactionMode == .move {
        guard let (nx, ny, aspect) = gizmoNDC(in: view, at: loc) else { return }
        // Dead-zone: until the pointer has travelled past the slop, keep the
        // press a candidate TAP (don't arm didDrag or start orbit/handle
        // drag). This is what makes tap-to-select fire on a slightly shaky
        // click instead of being swallowed as a tiny camera orbit.
        if !didDrag,
          hypot(loc.x - mouseDownLoc.x, loc.y - mouseDownLoc.y) < moveDragSlop
        {
          return
        }
        if let h = moveHandle {
          // Dragging a gizmo handle → manipulate the active object.
          if !didDrag {
            didDrag = true
            if let (dnx, dny, _) = gizmoNDC(in: view, at: mouseDownLoc) {
              engine?.gizmoBeginDrag(h, ndcX: dnx, ndcY: dny, aspect: aspect)
            }
          }
          engine?.gizmoUpdateDrag(ndcX: nx, ndcY: ny, aspect: aspect)
          return
        }
        // Empty-space drag → orbit the camera. The gizmo is a 3D CGO that
        // re-renders from the new camera automatically, so we do NOT refresh
        // its 2D hit-test geometry per tick — that per-tick Python round-trip
        // + JSON file I/O is what made orbiting laggy (even with no active
        // object). It's refreshed once on mouseUp for the next hover/click.
        if !didDrag {
          didDrag = true
          let down = pymolPoint(in: view, at: mouseDownLoc)
          engine?.button(
            PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_DOWN, x: down.0, y: down.1, modifiers: mods)
        }
        let pt = pymolPoint(in: view, at: loc)
        engine?.drag(x: pt.0, y: pt.1, modifiers: mods)
        return
      }

      if !didDrag {
        // First movement: now send the button-down (at the press point)
        // so PyMOL enters rotate mode for this drag.
        didDrag = true
        // A rotate/drag begins: drop the hover preview (the pointer is no
        // longer just hovering) and reset the move gate.
        engine?.clearHoverPreview()
        lastHoverLoc = .zero
        let down = pymolPoint(in: view, at: mouseDownLoc)
        engine?.button(
          PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_DOWN, x: down.0, y: down.1, modifiers: mods)
      }
      let pt = pymolPoint(in: view, at: loc)
      engine?.drag(x: pt.0, y: pt.1, modifiers: mods)
    }

    func handleRightMouseDown(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      // Defer PyMOL's button-down: raising it here would immediately enter
      // clip mode, and PyMOL's own pop-up menu is never rendered under this
      // Metal backend (internal_gui=0). A bare right-click instead pops the
      // native SwiftUI context menu in handleRightMouseUp; a right-DRAG
      // starts the clip on first movement (handleRightMouseDragged).
      rightDownLoc = view.convert(event.locationInWindow, from: nil)
      rightDidDrag = false
      // Right pressed while left is already down: a chord, not a menu click.
      rightChorded = NSEvent.pressedMouseButtons & 0b1 != 0
      if rightChorded && engine?.interactionMode != .move && moveHandle == nil {
        // Left first, then right: stop the rotate and zoom from here on.
        let mods = pymolModifiers(event.modifierFlags.rawValue)
        let pt = pymolPoint(in: view, at: rightDownLoc)
        if didDrag {
          engine?.button(
            PYMOL_BUTTON_LEFT, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1, modifiers: mods)
        }
        didDrag = true  // a drag, not a click, whatever follows
        leftChorded = true
        rightOverLeft = true
        rightDidDrag = true
        engine?.button(
          PYMOL_BUTTON_RIGHT, state: PYMOL_BUTTON_DOWN, x: pt.0, y: pt.1, modifiers: mods)
      }
    }

    func handleRightMouseUp(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      let loc = view.convert(event.locationInWindow, from: nil)
      let mods = pymolModifiers(event.modifierFlags.rawValue)

      if rightDidDrag {
        // Finish the clip (slab) drag.
        rightChorded = false
        let pt = pymolPoint(in: view, at: loc)
        engine?.button(
          PYMOL_BUTTON_RIGHT, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1, modifiers: mods)
        if rightOverLeft {
          rightOverLeft = false
          // Left still held: further left drags rotate again, from here.
          if NSEvent.pressedMouseButtons & 0b1 != 0 {
            didDrag = false
            mouseDownLoc = loc
          }
        }
        return
      }
      if rightChorded {
        rightChorded = false
        return
      }  // no scene menu after a chord

      // Pure right-click → identify the atom/residue under the cursor and
      // publish it so ContentView shows the viewport context menu. NDC in
      // view-point space, bottom-left origin (macOS views aren't flipped),
      // matching the left-click pick — so no Y flip.
      let w = view.bounds.width
      let h = view.bounds.height
      guard w > 0, h > 0 else { return }
      let ndcX = Float(loc.x / w) * 2 - 1
      let ndcY = Float(loc.y / h) * 2 - 1
      engine?.longPressPick(ndcX: ndcX, ndcY: ndcY, aspect: Float(w / h))
    }

    func handleRightMouseDragged(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      let loc = view.convert(event.locationInWindow, from: nil)
      let mods = pymolModifiers(event.modifierFlags.rawValue)
      let moved = hypot(loc.x - rightDownLoc.x, loc.y - rightDownLoc.y)
      if !rightDidDrag {
        // First real movement past the click threshold: now raise the
        // deferred button-down (at the press point) to begin the clip.
        guard moved >= 4 else { return }
        rightDidDrag = true
        let down = pymolPoint(in: view, at: rightDownLoc)
        engine?.button(
          PYMOL_BUTTON_RIGHT, state: PYMOL_BUTTON_DOWN, x: down.0, y: down.1, modifiers: mods)
      }
      let pt = pymolPoint(in: view, at: loc)
      engine?.drag(x: pt.0, y: pt.1, modifiers: mods)
    }

    // Where the middle button went down: released within 4 pt is a click.
    private var middleDownLoc: CGPoint?

    func handleOtherMouseDown(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      middleDownLoc = view.convert(event.locationInWindow, from: nil)
      let pt = pymolPoint(in: view, at: view.convert(event.locationInWindow, from: nil))
      let mods = pymolModifiers(event.modifierFlags.rawValue)
      engine?.button(
        PYMOL_BUTTON_MIDDLE, state: PYMOL_BUTTON_DOWN, x: pt.0, y: pt.1, modifiers: mods)
    }

    func handleOtherMouseUp(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      let pt = pymolPoint(in: view, at: view.convert(event.locationInWindow, from: nil))
      let mods = pymolModifiers(event.modifierFlags.rawValue)
      engine?.button(PYMOL_BUTTON_MIDDLE, state: PYMOL_BUTTON_UP, x: pt.0, y: pt.1, modifiers: mods)
      // A middle CLICK centers on the atom under it, as in desktop PyMOL. The
      // core's own "Cent" needs GL picking, which the Metal path lacks, so
      // use the CPU picker like a left click does.
      let loc = view.convert(event.locationInWindow, from: nil)
      if let down = middleDownLoc, hypot(loc.x - down.x, loc.y - down.y) < 4, mods == 0 {
        let w = view.bounds.width
        let h = view.bounds.height
        if w > 0, h > 0 {
          engine?.centerAt(
            ndcX: Float(loc.x / w) * 2 - 1, ndcY: Float(loc.y / h) * 2 - 1,
            aspect: Float(w / h))
        }
      }
      middleDownLoc = nil
    }

    func handleOtherMouseDragged(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      let pt = pymolPoint(in: view, at: view.convert(event.locationInWindow, from: nil))
      let mods = pymolModifiers(event.modifierFlags.rawValue)
      engine?.drag(x: pt.0, y: pt.1, modifiers: mods)
    }

    func handleScrollWheel(_ event: NSEvent, in view: MTKView) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      let loc = view.convert(event.locationInWindow, from: nil)
      let pt = pymolPoint(in: view, at: loc)
      let mods = pymolModifiers(event.modifierFlags.rawValue)

      let phase = event.phase
      let momentum = event.momentumPhase

      // A traditional scroll WHEEL has no touch phase (a trackpad / Magic
      // Mouse gesture always sets phase or momentumPhase). Route the wheel
      // to PyMOL's default bare-wheel binding = SLAB (clip), via the scroll
      // button. Touch-surface two-finger scroll falls through to PAN below.
      if phase == [] && momentum == [] {
        let wheel = event.scrollingDeltaY != 0 ? event.scrollingDeltaY : event.deltaY
        guard wheel != 0 else { return }
        let btn: Int32 = wheel > 0 ? PYMOL_BUTTON_SCROLL_FORWARD : PYMOL_BUTTON_SCROLL_REVERSE
        engine?.button(btn, state: PYMOL_BUTTON_DOWN, x: pt.0, y: pt.1, modifiers: mods)
        return
      }

      // Trackpad two-finger drag. Latch the mode at the START: Shift held →
      // CLIP (Shift+Right-drag), else PAN (Middle-drag). Real trackpad
      // gestures begin with phase == .began; synthetic/no-phase precise
      // scrolls (and momentum without a prior .began) start on first delta.
      let scale = view.window?.backingScaleFactor ?? 2.0
      if !panActive && (phase == .began || (phase == [] && momentum != .ended)) {
        gestureIsClip = event.modifierFlags.contains(.shift)
        gestureButton = gestureIsClip ? PYMOL_BUTTON_RIGHT : PYMOL_BUTTON_MIDDLE
        gestureMods = gestureIsClip ? PYMOL_MOD_SHIFT : 0
        gestureIsSlab = false
        panActive = true
        panCursorX = pt.0
        panCursorY = pt.1
        if engine?.trackpadMode == true && !gestureIsClip {
          axisUndecided = true
          axisTravel = .zero
        } else {
          engine?.button(
            gestureButton, state: PYMOL_BUTTON_DOWN,
            x: panCursorX, y: panCursorY, modifiers: gestureMods)
        }
      }

      if axisUndecided {
        axisTravel.x += event.scrollingDeltaX
        axisTravel.y += event.scrollingDeltaY
        guard max(abs(axisTravel.x), abs(axisTravel.y)) >= kAxisDecide else {
          armPanEndDebounce()
          return
        }
        axisUndecided = false
        if abs(axisTravel.y) > abs(axisTravel.x) {
          gestureIsSlab = true
          slabAccum = axisTravel.y
        } else {
          engine?.button(
            gestureButton, state: PYMOL_BUTTON_DOWN,
            x: panCursorX, y: panCursorY, modifiers: gestureMods)
        }
      }

      if gestureIsSlab {
        // Fingers on the pad only: a momentum glide would keep changing it.
        if momentum == [] { slabAccum += event.scrollingDeltaY }
        while abs(slabAccum) >= kSlabStep {
          let forward = slabAccum > 0
          slabAccum -= forward ? kSlabStep : -kSlabStep
          engine?.button(
            forward ? PYMOL_BUTTON_SCROLL_FORWARD : PYMOL_BUTTON_SCROLL_REVERSE,
            state: PYMOL_BUTTON_DOWN, x: pt.0, y: pt.1, modifiers: 0)
        }
        if phase == .cancelled || momentum == .ended {
          endTrackpadPan()
        } else {
          armPanEndDebounce()
        }
        return
      }

      let signX = gestureIsClip ? kClipSignX : kPanSignX
      let signY = gestureIsClip ? kClipSignY : kPanSignY
      let dx = Int32((event.scrollingDeltaX * scale * signX).rounded())
      let dy = Int32((event.scrollingDeltaY * scale * signY).rounded())

      if panActive && (dx != 0 || dy != 0) {
        // macOS views are bottom-left origin (matching PyMOL); a finger
        // moving up has positive scrollingDeltaY, so add directly.
        // For CLIP, hold X fixed so only the near plane (vertical →
        // front clip) moves — horizontal would move the far plane,
        // closing the slab from both sides.
        if !gestureIsClip { panCursorX += dx }
        panCursorY += dy
        engine?.drag(x: panCursorX, y: panCursorY, modifiers: gestureMods)
      }

      // End when the momentum glide finishes (the true end), or on cancel.
      // We deliberately DON'T end at phase == .ended (fingers up): momentum
      // events follow with phase == [] and would re-trigger the start
      // condition, restarting the drag mid-glide. The debounce is the
      // safety net for flicks that produce no momentum and for synthetic
      // no-phase event streams.
      if phase == .cancelled || momentum == .ended {
        endTrackpadPan()
      } else if panActive {
        armPanEndDebounce()
      }
    }

    private func armPanEndDebounce() {
      panEndDebounce?.cancel()
      let work = DispatchWorkItem { [weak self] in self?.endTrackpadPan() }
      panEndDebounce = work
      DispatchQueue.main.asyncAfter(deadline: .now() + 0.12, execute: work)
    }

    private func endTrackpadPan() {
      panEndDebounce?.cancel()
      panEndDebounce = nil
      guard panActive else { return }
      panActive = false
      let pressed = !axisUndecided && !gestureIsSlab
      axisUndecided = false
      gestureIsSlab = false
      guard pressed else { return }
      engine?.button(
        gestureButton, state: PYMOL_BUTTON_UP,
        x: panCursorX, y: panCursorY, modifiers: gestureMods)
    }

    @objc func handleMagnification(_ gesture: NSMagnificationGestureRecognizer) {
      guard !boxSelectActive else { return }  // camera frozen (#358)
      switch gesture.state {
      case .began:
        lastMag = 0
      case .changed:
        // Spreading fingers (magnification increasing) = zoom in.
        let delta = gesture.magnification - lastMag
        lastMag = gesture.magnification
        engine?.zoomBy(Float(delta * kZoomGain))
      case .ended, .cancelled:
        lastMag = 0
        // Zoom changed the projection scale → refresh the gizmo hit-test
        // once (guarded to move mode) so a handle grab after zooming is
        // accurate. Not done per .changed tick — that would relag zoom.
        engine?.refreshGizmo()
      default:
        break
      }
    }

    // Trackpad two-finger twist → Z-axis roll, mirroring the iOS handleRotation
    // gesture. Gated on engine.trackpadMode so it never disturbs the classic
    // per-button mouse modes. runPython (not runCommand) to avoid echoing
    // `turn z` into the feedback log every frame.
    @objc func handleRotationGesture(_ gesture: NSRotationGestureRecognizer) {
      guard engine?.trackpadMode == true, !boxSelectActive else { return }
      switch gesture.state {
      case .began:
        lastRoll = 0
      case .changed:
        let delta = gesture.rotation - lastRoll
        lastRoll = gesture.rotation
        let deg = kRollSign * Float(delta) * 180.0 / .pi
        engine?.runPython("from pymol import cmd as _c; _c.turn('z', \(deg))")
      case .ended, .cancelled:
        lastRoll = 0
      default:
        break
      }
    }
  }
}
