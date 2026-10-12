// PyMOLEngine.swift — Singleton managing the PyMOL instance lifecycle
// Wraps the C bridge API and publishes observable state for SwiftUI views.

import Combine
import Foundation
import MetalKit
import SwiftUI

/// Frequently-changing movie/timeline playback state, kept in its OWN
/// ObservableObject so the ≈10/s frame ticks during playback re-render only the
/// TransportBar (the lone observer) and not the inspector/menus that observe
/// PyMOLEngine. Mirrors the core; updated by PyMOLEngine.parsePlaybackFeedback.
final class PlaybackState: ObservableObject {
  @Published var currentFrame: Int = 1  // 1-based, == cmd.get_frame()
  @Published var frameCount: Int = 1  // == cmd.count_frames()
  @Published var isPlaying: Bool = false  // == cmd.get_movie_playing()
  @Published var movieFPS: Double = 15  // == movie_fps setting
  @Published var movieLoop: Bool = true  // == movie_loop setting
}

enum MeasureKind: String { case distance, angle, dihedral }

/// One PyMOL setting for the Settings panel. type: 1 bool, 2 int, 3 float,
/// 4 float3, 5 color, 6 string (pymol.setting type codes).
struct SettingItem: Identifiable, Equatable, Codable {
  let name: String
  let type: Int
  var val: String
  var id: String { name }
}

/// The atom/residue under a right-click, for the viewport context menu.
/// `isEmpty` means the click landed on background (no atom) → scene-level actions.
struct LongPressHit: Equatable, Identifiable {
  let id = UUID()
  var isEmpty: Bool
  var obj = "", chain = "", resi = "", resn = "", name = "", sel = ""
  /// Menu title, e.g. "ALA 42 · chain A" (or "Scene" for empty space).
  var title: String {
    if isEmpty { return "Scene" }
    var t = resn.isEmpty ? resi : "\(resn) \(resi)"
    if !chain.isEmpty { t += " · chain \(chain)" }
    return t
  }
}

/// The console's lines, observed only by the console view, so output doesn't
/// re-render the object panel and the rest of the window line by line.
final class ConsoleLog: ObservableObject {
  @Published var lines: [String] = [] {
    // Cap at 400, counting what the cap drops. The console diffs with
    // `trimmed` instead of guessing the overlap from the text, which goes
    // wrong on repeated lines (a startup script's identical tracebacks).
    didSet {
      let excess = lines.count - 400
      if excess > 0 {
        trimmed += excess
        lines.removeFirst(excess)
      }
    }
  }
  private(set) var trimmed = 0

  /// Empties the log (the `clear` command).
  func clear() { lines = [] }
}

final class PyMOLEngine: ObservableObject {
  static let shared = PyMOLEngine()

  let console = ConsoleLog()
  var feedbackLog: [String] {
    get { console.lines }
    set { console.lines = newValue }
  }
  /// Last SEQCOLOR fingerprint (see fetchSequenceSelection).
  private var lastSequenceColorDigest: String?
  @Published var objects: [MoleculeObject] = []
  // Fires once per object-panel poll, whether or not anything changed. The guarded
  // @Published lists above emit only on change.
  let panelPolled = PassthroughSubject<Void, Never>()
  @Published var sequences: [SequenceObject] = [] { didSet { sequenceVersion &+= 1 } }
  @Published var selectedResidueKeys: Set<String> = [] { didSet { sequenceVersion &+= 1 } }
  /// Bumped when the strip's data changes, so it can skip the engine's other
  /// updates (SequencePanel.StripContent).
  private(set) var sequenceVersion = 0
  /// A ▸ Contact Map…: the object whose map tab should open (ContentView
  /// opens it and clears this).
  @Published var contactMapRequest: String? = nil
  /// Contact maps open as title-bar tabs (in this order) that replace the 3D
  /// view; nil active = the 3D view ("mMol" tab).
  @Published var contactMapTabs: [String] = []
  @Published var activeContactMap: String? = nil
  /// Open tool tabs and the one showing (ToolTabs.swift); nil = mMol or a plot tab.
  @Published var toolTabs: [ToolTab] = []
  @Published var activeTool: ToolTab? = nil
  /// A ▸ Align ▸ Align Wizard…: the object to preselect as mobile.
  @Published var alignMobile: String? = nil
  @Published var renderOptionSpec: [RenderOption] = []
  @Published var renderOptionValues: [String: Any] = [:]
  /// Edit ▸ New Selection… / all ▸ A ▸ New Selection…: show the selection builder.
  @Published var showSelectionBuilder = false
  /// Scene ▸ Rename…: the scene being renamed (ContentView asks for the name).
  @Published var pendingSceneRename: String? = nil

  func closeContactMapTab(_ name: String) {
    contactMapTabs.removeAll { $0 == name }
    if activeContactMap == name { activeContactMap = nil }
  }
  /// C ▸ by file…: each object's CSV heatmap (Heatmap.swift).
  @Published var heatmaps: [String: HeatmapState] = [:]
  /// Objects whose card shows the heatmap panel rather than a layer. Kept
  /// here, not in the card's @State: the card is rebuilt as the list scrolls
  /// and expands, which reset a local choice right after loading a CSV.
  @Published var heatmapChipActive: Set<String> = []
  // Set when a right-click identifies an atom/residue (or empty space);
  // ContentView observes this to present the context menu, then clears it.
  @Published var longPressHit: LongPressHit?
  @Published var isReady = false
  // Restored from the last launch (#332). A session (.pse) that turns seq_view
  // on still wins — it assigns this property after launch, like any other setter.
  @Published var sequenceVisible = UserDefaults.standard
    .bool(forKey: PanelLayout.sequenceVisibleKey)
  {
    // Showing the strip must (re)fetch the sequence data — toggling the
    // panel on (menu/toolbar) only flipped this bool, so the strip rendered
    // empty until an unrelated load/refresh happened to repopulate it (the
    // "sequence only shows when the console is open" bug). Fetch on the
    // false→true edge; fetchSequences() no-ops until the engine is ready,
    // and the load path re-fetches then.
    didSet {
      UserDefaults.standard.set(sequenceVisible, forKey: PanelLayout.sequenceVisibleKey)
      if sequenceVisible && !oldValue { fetchSequences() }
    }
  }
  // Names the sequence strip was last built from (issue #380): the panel only
  // rows ENABLED objects, so any path that changes what is enabled — the object
  // panel checkbox, "all", a group cascade, `disable` at the prompt, a scene or
  // session recall — has to republish it. Rather than hooking each of those,
  // parseObjectPanelFeedback compares the core's authoritative enabled set
  // against this and re-fetches on a mismatch, which is also immune to the
  // optimistic-UI flip racing the enable/disable command.
  var lastSequenceEnabled: Set<String> = []

  // macOS "Trackpad" mouse mode: when on, the MetalViewport trackpad gesture
  // handlers (scroll/pinch) drive rotate/pan/zoom/roll/clip directly, rather
  // than routing the bare wheel to PyMOL's
  // mouse_mode binding (slab). Toggled by selecting "Trackpad" in MousePanel;
  // the chosen mode is persisted via @AppStorage in MousePanel.
  @Published var trackpadMode: Bool = false

  // Current MTKView drawable size in backing pixels (updated by the viewport
  // coordinator on resize). Used by the Export menu's "current size" render.
  @Published var viewportPixelSize: CGSize = .zero

  // While true (set during a panel-divider drag), the MetalViewport freezes its
  // drawable size so a fast drag doesn't reallocate all offscreen targets every
  // frame (choppy + OOM). Cleared on release → one reshape at the final size.
  var suppressDrawableResize = false
  // How the frozen frame fits the viewport meanwhile: scaled for the console
  // divider (PyMOL shrinks the scene with the viewport's height), centered for
  // the side seam (a narrower viewport just shows less of the sides).
  var frozenDrawableGravity: CALayerContentsGravity = .resizeAspect

  // While a movie export is rendering, it owns the core off the main thread
  // (renderHiResPNG temporarily reshapes global core state + blocks on the
  // GPU). The live draw loop and the feedback poll — the other main-thread
  // core accessors — gate on this flag and skip, so the exporter is the sole
  // core user and the per-frame render no longer beachballs the UI. Set/cleared
  // on the main thread by MovieExporter; the export sheet is modal so no
  // interleaved commands run meanwhile. (#58 L-59)
  var exportRenderActive = false

  // Representation inspector state.
  // Per-object active representations + their current setting values + color
  // override, populated by pollDetails()/parseObjectDetailFeedback().
  @Published var objectDetails: [String: [RepState]] = [:]

  // Global "Scene" parameters (metal_*, depth_cue, fog, fov, surface_quality, bg).
  @Published var sceneState = SceneState()
  // Per-object state metadata (effective current state + overlay-all) for the
  // inspector STATE row; populated alongside objectDetails for expanded cards.
  @Published var objectMeta: [String: ObjStateMeta] = [:]
  // Saved scenes (ordered) + the current scene name, for the Scenes strip.
  @Published var sceneNames: [String] = []
  @Published var currentScene: String = ""
  // Live atom-count preview for the selection builder (nil = not previewing).
  @Published var selectionPreviewCount: Int? = nil
  // Set by the action menu's "Rename" to request a rename modal for this object
  // name (PyMOL's `wizard renaming` has no UI on the Metal/SwiftUI app). The
  // ObjectPanel observes this and presents a name-entry alert.
  @Published var pendingRename: String? = nil
  // Request channel for "New Group…" (#255): the object to put in a new group.
  // ObjectPanel presents the name-entry alert, mirroring pendingRename.
  @Published var pendingGroupFor: String? = nil
  // Request channel for "Copy to Object ▸ New Object…" (#461): the selection to
  // copy into a brand-new object. Same shape as pendingGroupFor — ObjectPanel
  // presents the name-entry alert, prefilled with PyMOL's own objNN default.
  @Published var pendingCopyToNew: String? = nil
  // Pick-debug instrumentation (active when PYMOL_PICKDEBUG is set): the last
  // click point in viewport (top-down, SwiftUI) points, drawn as a crosshair so
  // a screenshot shows click-vs-selection alignment. Set by MetalViewport.
  static let debugPickEnabled = ProcessInfo.processInfo.environment["PYMOL_PICKDEBUG"] != nil
  @Published var debugClickPoint: CGPoint? = nil
  // Interactive measurement: nil = off; otherwise taps pick atoms for a
  // distance/angle/dihedral. measureStatus is the prompt/result shown in the UI.
  @Published var measureMode: MeasureKind? = nil
  @Published var measureStatus: String = ""
  // Move mode: rigid-body object manipulation via an on-screen gizmo (driven by
  // mmol/move.py / TTT matrix). interactionMode gates gesture routing in
  // MetalViewport; gizmo holds the projected handle geometry the overlay draws
  // and the handlers hit-test.
  @Published var interactionMode: InteractionMode = .viewing
  @Published var activeMoveObject: String? = nil
  // Box Select (#358): the rubber-band rectangle, in viewport NDC. nil = the
  // mode is on but no box has been drawn yet. interactionMode == .boxSelect
  // gates gesture routing in MetalViewport (which freezes the camera while it
  // is on); boxDrag is the in-flight grab, live only between press and release.
  @Published var boxRect: BoxRect? = nil
  var boxDrag: BoxDrag? = nil
  // Adjust-frame mode: gizmo controls re-anchor the gizmo's own frame (origin +
  // inclination) instead of moving the structure — for precise custom pivots.
  // Active when the overlay toggle is on OR Shift is held (macOS). The gizmo
  // renders gray + semitransparent while active. Synced to metal_move.set_adjust.
  @Published var adjustFrameToggle = false { didSet { syncAdjustFrame() } }
  // @Published so the overlay toggle reflects Shift; callers must only assign it
  // on an actual change (it's polled on every hover move) to avoid re-render spam.
  @Published var moveShiftHeld = false { didSet { syncAdjustFrame() } }
  private var lastAdjustSent = false
  var adjustFrameActive: Bool { interactionMode == .move && (adjustFrameToggle || moveShiftHeld) }
  // macOS hover highlight — pushes the hovered handle to metal_move so the CGO
  // gizmo emphasizes it. Rebuilds are infrequent (only when the handle changes).
  @Published var hoveredHandle: GizmoHandle? = nil {
    didSet {
      guard interactionMode == .move, oldValue != hoveredHandle else { return }
      runPython(
        "from mmol import move as _mm\n_mm.set_hover('\(hoveredHandle?.pyName ?? "")')")
      // set_hover rebuilt the gizmo CGO with the new handle emphasized, but a
      // static scene's on-demand draw gate can consume the redisplay flag
      // before the next tick and leave the highlight unshown (the "hover does
      // nothing" symptom). Force one repaint so the emphasis actually appears.
      requestViewportRedraw()
    }
  }
  @Published var gizmo: GizmoGeometry? = nil
  // Hover pre-selection preview (issue #165, macOS): when true, moving the
  // pointer over the viewport highlights what a click WOULD select (in a
  // distinct light-cyan) without committing to 'sele'. User-toggleable in
  // Viewer Settings ▸ Canvas. When false, hoverPreview() is a no-op.
  // Persisted across launches via UserDefaults (didSet); default on.
  @Published var hoverPreviewEnabled: Bool =
    (UserDefaults.standard.object(forKey: PyMOLEngine.hoverDefaultsKey) as? Bool ?? true)
  {
    didSet { UserDefaults.standard.set(hoverPreviewEnabled, forKey: PyMOLEngine.hoverDefaultsKey) }
  }
  static let hoverDefaultsKey = "mmol.hoverPreviewEnabled"
  // Full settings catalog for the searchable Settings panel (loaded on demand).
  @Published var settingsCatalog: [SettingItem] = []
  // The single detail view that is currently open (accordion: at most one).
  // nil = none, `sceneDetailKey` = the SCENE card, otherwise an object name.
  // Drives which object the detail poll queries (collapsed = cheap).
  // Poll the just-expanded object's rep detail IMMEDIATELY here (at the source)
  // rather than waiting up to ~500ms for the next pollObjects tick — a heavy
  // surface build can starve that timer and the card lingers on "No
  // representations shown" / renders empty (#107). Doing it in didSet (instead
  // of a per-layout view .onChange) guarantees EVERY expand path fires the
  // refresh: the layout (which had no such observer — the empty-panel bug),
  // session restore, and the PYMOL_AUTOEXPAND test hook.
  @Published var expandedDetail: String? = nil {
    didSet {
      guard expandedDetail != oldValue, expandedDetail != nil else { return }
      refreshExpandedDetail()
    }
  }
  // Sentinel for "the SCENE card is the open detail view" — a control char so
  // it can never collide with a real object name.
  static let sceneDetailKey = "\u{1}scene"

  // Reps the user toggled invisible but wants kept as listed "layers" (per
  // object) so they stay in the inspector and can be toggled back on, instead
  // of vanishing the moment they're hidden. Modified ONLY by explicit user
  // actions (hide keeps, show/delete removes) — never by the detail poll — so
  // there's no re-add race with the ~500ms poll.
  @Published var keptHidden: [String: Set<String>] = [:]
  /// Reps the user turned on per object (by command, menu, or while its card
  /// was open): the card lists these beyond cartoon. lastShownReps spots reps
  /// that appear while a card is open.
  @Published var cardReps: [String: Set<String>] = [:]
  private var lastShownReps: [String: Set<String>] = [:]
  private var detailRefresh: DispatchWorkItem?

  // The .pse session file the user currently has open (Finder open / drag-drop /
  // Open… of a .pse, or the destination of a Save As). nil = never-saved session
  // (or a non-.pse structure was opened). Drives ⌘S overwrite vs Save As, and the
  // macOS window title. Cleared by Clear Session.
  @Published var currentSessionURL: URL? = nil

  // MARK: Timeline / playback (states · trajectories · movies)
  // In PyMOL these are ONE concept: a 1-based movie frame index that maps
  // (via mset) to a coordinate state. The CORE drives frame advance
  // (cmd.mplay → SceneIdle, ticked every Metal frame); Swift only scrubs
  // (cmd.frame) and mirrors core state.
  //
  // The per-frame state lives in a SEPARATE ObservableObject so the ≈10/s
  // frame ticks during playback re-render ONLY the TransportBar (which
  // observes `playback`) and NOT the inspector/menus (which observe `engine`).
  // Re-rendering the inspector on every tick was resetting open menus to the
  // top and dropping in-flight button touches during playback.
  let playback = PlaybackState()
  // Lightweight gate the rest of the UI observes to show/hide the transport
  // bar. Flips only when a multi-state object appears/disappears, so it does
  // NOT cause per-frame re-renders of views that observe `engine`.
  @Published var hasTimeline: Bool = false
  // Camera keyframe frames the user has captured this session (1-based, sorted),
  // used to draw the Timeline's camera track. Manual keyframing only — template
  // builders author their own mview keyframes and clear this list.
  // Session-scoped: not reconstructed from a reloaded .pse (a Phase-2 concern).
  @Published var cameraKeyframes: [Int] = []
  // Scenes placed AS time-positioned markers on the timeline (drag-and-drop),
  // as opposed to the recall palette. Each is an mview keyframe tagged with a
  // scene (camera flies between scenes, reps cut). One keyframe per frame, so
  // sceneMarkers frames are disjoint from cameraKeyframes. Session-scoped too.
  struct SceneMarker: Identifiable, Equatable {
    let frame: Int
    let name: String
    var id: Int { frame }  // one keyframe per frame → frame is a stable id
  }
  @Published var sceneMarkers: [SceneMarker] = []

  // MARK: Unified timeline track
  //
  // The docked TimelinePanel edits ONE ordered lane of "objects" — camera
  // keyframes and scene markers — joined by transitions (duration + easing).
  // Swift is the source of truth for the ORDER, TIMING and EASING (the things
  // the UI edits); the heavy camera view matrices live in Python
  // (appkit_movie._views) keyed by each camera item's UUID. Every edit
  // recomputes frames from the transition durations and rebuilds the whole
  // core movie via appkit_movie.rebuild. Session-scoped, like the legacy
  // tracks above (not reconstructed from a reloaded .pse).
  struct Transition: Equatable {
    var seconds: Double = 2.0  // gap from the previous item to this one
    var linear: Bool = false  // false = smooth (eased) default
  }
  enum StatesMode: String, Equatable { case sweep, loop, lockstep }
  // A "Play models" clip: sweep a model range of one/all multi-state objects.
  struct StatesSpec: Equatable {
    var objects: [String]? = nil  // nil = all multi-state objects
    var mode: StatesMode = .sweep
    var firstModel: Int = 1  // play models firstModel..lastModel
    var lastModel: Int = 0  // 0 = through the object's last model
    var durationSeconds: Double = 4.0
  }
  struct TimelineItem: Identifiable, Equatable {
    enum Kind: Equatable {
      case camera
      case scene(name: String)
      case states(StatesSpec)
    }
    let id: UUID
    var kind: Kind
    var transition: Transition  // transition INTO this item (item[0]'s ignored)
    var pinnedState: Int? = nil  // camera/scene: hold this model (opt out of auto-sweep)
    // Explicit timeline frame. When set, the item sits at this absolute frame
    // (may overlap a states clip) instead of the sequential transition layout —
    // used for camera keyframes dropped at the playhead during a clip.
    var atFrame: Int? = nil
    init(
      id: UUID = UUID(), kind: Kind, transition: Transition = Transition(),
      pinnedState: Int? = nil, atFrame: Int? = nil
    ) {
      self.id = id
      self.kind = kind
      self.transition = transition
      self.pinnedState = pinnedState
      self.atFrame = atFrame
    }
  }
  @Published var timelineItems: [TimelineItem] = []

  // While the user drags the scrubber, the poll must NOT overwrite
  // currentFrame (classic two-way-binding fight). Set on drag start, cleared
  // shortly after release so the next poll can re-sync.
  var isScrubbing: Bool = false
  private var scrubReleaseWork: DispatchWorkItem?
  private var lastScrubFrame: Int = -1

  // The opaque PyMOL instance pointer
  private(set) var instance: PyMOLHandle?

  private var feedbackTimer: Timer?

  // The progress bar of a surface / mesh / dots build, held until this many
  // render frames complete (the build happens during them).
  private var buildBar: RayProgressBar?
  private var pendingHeavyClearFrames = 0
  // Backstop that hides a stuck bar if the render loop ever stalls.
  private var busyBackstop: DispatchWorkItem?

  private init() {
  }

  // MARK: - Lifecycle

  func initialize(resourcePath: String) {
    guard instance == nil else { return }

    instance = PyMOLBridge_New()
    guard let inst = instance else { return }

    // Point the embedded Python's TLS at the bundled CA bundle BEFORE init,
    // so `fetch` (HTTPS to RCSB/PDB) can verify certificates. Must precede
    // Py init.
    if let ca = Bundle.main.path(forResource: "cacert", ofType: "pem", inDirectory: "data") {
      setenv("SSL_CERT_FILE", ca, 1)
    }

    PyMOLBridge_InitPython(inst, resourcePath)
    PyMOLBridge_Start(inst)

    isReady = true

    // `fetch` downloads into fetch_path. Use the temp directory: it's always
    // writable, is the app's own (container) tmp under the sandbox, and is NOT
    // TCC-protected — so fetching never triggers the macOS "access your
    // Documents folder" prompt. Fetched structures are transient cache, not
    // user documents, so a temp location is the right home for them.
    let tmp = NSTemporaryDirectory()
    runPython("from pymol import cmd as _c; _c.set('fetch_path', '\(tmp)')")

    // Test affordance: run arbitrary ;-separated PyMOL commands at launch.
    if let c = ProcessInfo.processInfo.environment["PYMOL_AUTOCMD"] {
      for one in c.split(separator: ";") {
        runCommand(one.trimmingCharacters(in: .whitespaces))
      }
    }

    // Blank-on-launch guard: the on-demand draw gate can leave the very first
    // frame — where deferred rep geometry (cartoon/surface) is still being
    // built — unshown, so the scene stays blank until the user interacts
    // (loading a second object or orbiting "fixed" it). One forced frame builds
    // the geometry; the follow-up forced frames actually paint it. Cheap
    // (a handful of extra frames at startup) and harmless once rendered.
    for delay in [0.15, 0.35, 0.6, 1.0, 1.6] {
      DispatchQueue.main.asyncAfter(deadline: .now() + delay) { [weak self] in
        self?.requestViewportRedraw()
      }
    }

    // Test affordance: open the object card of PYMOL_AUTOEXPAND=<object>.
    if let ex = ProcessInfo.processInfo.environment["PYMOL_AUTOEXPAND"] {
      DispatchQueue.main.asyncAfter(deadline: .now() + 3.0) { [weak self] in
        self?.expandedDetail = ex
      }
    }

    // Poll feedback every 100ms
    feedbackTimer = Timer.scheduledTimer(withTimeInterval: 0.1, repeats: true) { [weak self] _ in
      guard let self = self else { return }
      self.pollFeedback()
      self.drainMCPMainQueue()
      self.pollObjects()
      // While the core is advancing frames, mirror the frame counter at
      // the full 100ms tick so the scrubber tracks playback smoothly.
      // When idle, the cheaper 500ms pollObjects() discovery suffices.
      if self.playback.isPlaying { self.pollPlayback() }
    }
  }

  // Ask the core for the current frame / length / play state. Cheap (a few
  // cmd gets); the PLAYBACK: line is read on the next pollFeedback tick.
  private func pollPlayback() {
    runPythonQuiet("from mmol import movie as _am\n_am.poll()")
  }

  func shutdown() {
    feedbackTimer?.invalidate()
    feedbackTimer = nil
    if let inst = instance {
      PyMOLBridge_Free(inst)
    }
    instance = nil
    isReady = false
  }

  // MARK: - Commands

  /// Toggle an object/selection's enabled state with an OPTIMISTIC UI update.
  /// Flips the matching `objects[idx].isEnabled` immediately on the main thread
  /// (ObjectEntry.isEnabled is a var → the checkbox re-renders this frame, no
  /// waiting on the ~500ms pollObjects round-trip), THEN issues the actual
  /// enable/disable command. The pollObjects/parseObjectPanelFeedback equality
  /// guard reconciles the array later; since we already set the value the poll
  /// will report, the reconcile is a no-op and there is no flicker.
  func setObjectEnabled(_ name: String, _ enabled: Bool) {
    // UI mutation must happen on the main thread (drives @Published).
    if Thread.isMainThread {
      applyEnabledOptimistically(name, enabled)
    } else {
      DispatchQueue.main.async { self.applyEnabledOptimistically(name, enabled) }
    }
    // Guarded like the Inspector's builders (#531): the visibility
    // checkbox is the most-used of them.
    runCommand(enabled ? "enable \(name)" : "disable \(name)", naming: name)
    // The sequence strip lists only ENABLED objects (#380): rebuild it.
    refetchSequencesSoon()
  }

  /// Bumped by every optimistic object-list change (visibility clicks). A panel
  /// poll issued before the latest bump carries a snapshot from before the
  /// click, which would flip the row back until the next poll: it is dropped.
  var objectListEpoch = 0
  /// Last scroll-wheel/trackpad event anywhere in the app (ContentView's monitor).
  var lastScrollAt: CFTimeInterval = 0
  var objectPollEpoch = 0
  private var sequenceRefetchQueued = false

  /// Rebuild the sequence strip after the current click is drawn (it costs tens
  /// of ms with many objects); quick successive clicks share one rebuild.
  func refetchSequencesSoon() {
    guard sequenceVisible, !sequenceRefetchQueued else { return }
    sequenceRefetchQueued = true
    DispatchQueue.main.async {
      self.sequenceRefetchQueued = false
      self.fetchSequences()
    }
  }

  /// Flip `name` and, when it is a group, its whole subtree.
  ///
  /// PyMOL cascades enable/disable down a group already, but the panel would not
  /// see that for up to a poll interval (~500ms), so a group tap would leave its
  /// members' checkboxes visibly stale and any tri-state parent wrong. Mirroring
  /// the cascade locally keeps the tree self-consistent until the poll confirms.
  private func applyEnabledOptimistically(_ name: String, _ enabled: Bool) {
    objectListEpoch += 1
    var targets: Set<String> = [name]
    // Walk down by parent links; cheap (the list is tens of entries) and the
    // seen-set makes a corrupt parent map terminate rather than spin.
    var frontier: Set<String> = [name]
    while !frontier.isEmpty {
      let children = objects.filter { $0.parent.map { frontier.contains($0) } ?? false }
        .map { $0.name }
        .filter { !targets.contains($0) }
      if children.isEmpty { break }
      targets.formUnion(children)
      frontier = Set(children)
    }
    for idx in objects.indices where targets.contains(objects[idx].name) {
      objects[idx].isEnabled = enabled
    }
  }

  func runCommand(_ command: String) {
    guard isReady else { return }
    // A plain `png <file>` (ray=0) wants the RENDERED frame, but PyMOL's
    // ScenePNG reads a GL framebuffer we don't have → it writes nothing.
    // Capture the Metal render ourselves instead. `ray=1` still uses the
    // (working) core CPU ray-trace path.
    if maybeCaptureRenderedPNG(command) { return }
    // Known long ops (surface build, in-place ray-trace, quality bumps) run
    // off-main so the "Calculating…" overlay can render + animate. The
    // overlay's blocking scrim prevents the user from issuing an interleaved
    // command mid-op, so selectively backgrounding these stays correctly
    // ordered. Light/interactive commands stay synchronous (snappy).
    if let label = heavyLabel(for: command) {
      runHeavy(label) { [weak self] in self?.runCommandCore(command) }
    } else {
      runCommandCore(command)
    }
  }

  // Synchronous command body. Safe on the main thread (light commands) or on
  // coreQueue (heavy commands): it calls only direct bridge/runPython helpers,
  // never runCommand, so it can't re-enter the heavy-dispatch path.
  private func runCommandCore(_ command: String) {
    PyMOLBridge_RunCommand(command)
    noteShownReps(command)
    clearRepColors(after: command)
    // The open card reflects the change now, not at the next ~500 ms poll.
    if expandedDetail != nil {
      detailRefresh?.cancel()
      let w = DispatchWorkItem { [weak self] in self?.pollDetails() }
      detailRefresh = w
      DispatchQueue.main.asyncAfter(deadline: .now() + 0.05, execute: w)
    }
    // The sequence strip only refetches when objects come, go or toggle, so
    // `remove resn HOH`, `delete`, recolouring etc. left it stale (waters
    // still listed after A ▸ Remove Waters). Refetch after any command that
    // can change which residues exist or their colours.
    if sequenceVisible && Self.changesSequence(command) { fetchSequences() }
    handleSessionViewport(for: command)
    maybeWidenClipForSurface(for: command)
    // A per-object `set state, N, obj` (and the related all_states / unset
    // controls) changes which model renders, but once a movie exists the
    // redisplay flag can be consumed before the viewport's on-demand gate
    // checks it — leaving the viewer frozen while the state counter advances
    // (issue #132). Force one repaint so the new state shows.
    let l = command.lowercased()
    if l.hasPrefix("set state") || l.hasPrefix("set_state")
      || l.hasPrefix("set all_states") || l.hasPrefix("unset state")
      || l.hasPrefix("unset all_states")
    {
      requestViewportRedraw()
    }
  }

  /// `show <rep>[, sel]` / `as <rep>[, sel]`: list that rep on the cards of the
  /// objects named in sel (all objects when it names none).
  private func noteShownReps(_ command: String) {
    let names = objects.filter { !$0.isSelection }.map(\.name)
    for stmt in command.split(whereSeparator: { $0 == "\n" || $0 == ";" }) {
      let words = stmt.trimmingCharacters(in: .whitespaces)
        .split(separator: ",", maxSplits: 1).map { $0.trimmingCharacters(in: .whitespaces) }
      let head = words.first?.split(separator: " ", maxSplits: 1).map(String.init) ?? []
      guard head.count == 2, head[0] == "show" || head[0] == "as",
        RepCatalog.order.contains(head[1])
      else { continue }
      let sel = words.count > 1 ? words[1] : ""
      let tokens = Set(
        sel.split(whereSeparator: {
          !($0.isLetter || $0.isNumber || $0 == "_" || $0 == "-" || $0 == ".")
        }).map(String.init))
      // Names typed in sel count even before the object list has caught up
      // (e.g. right after a fetch); stray words only add unused keys.
      let hit = names.filter(tokens.contains)
      for n in (hit.isEmpty ? names : hit) + tokens { cardReps[n, default: []].insert(head[1]) }
    }
  }

  /// S/H ▸ all states: draw every state of the objects in `name` (an object,
  /// selection or "all") at once, or just the current one.
  func setAllStates(_ name: String, _ on: Bool) {
    let n = name.replacingOccurrences(of: "'", with: "")
    runPython(
      "from pymol import cmd as _c\nfor _o in _c.get_object_list('(\(n))') or []:\n"
        + "    _c.set('all_states', \(on ? 1 : 0), _o)")
  }

  /// A color command (C menu, `color`, util.*, spectrum) wins over the card's
  /// per-rep colors on the objects it touched: clear those overrides.
  private func clearRepColors(after command: String) {
    for stmt in command.split(whereSeparator: { $0 == "\n" || $0 == ";" }) {
      var t = stmt.trimmingCharacters(in: .whitespaces)
      if t.hasPrefix("cmd.") { t.removeFirst(4) }  // the C menu's `python` blocks
      let verb = String(t.prefix { !($0 == " " || $0 == "(" || $0 == ",") }).lowercased()
      var sel = "all"
      if verb == "color" || verb == "colour" {
        sel = t.split(separator: ",", maxSplits: 1).dropFirst().first.map(String.init) ?? "all"
      } else if verb.hasPrefix("util.c") || verb.hasPrefix("util.ss") || verb == "util.rainbow" {
        sel =
          t.dropFirst(verb.count).trimmingCharacters(in: CharacterSet(charactersIn: " ()'\""))
          .split(separator: ",").first.map(String.init) ?? "all"
      } else if verb != "spectrum" {
        continue
      }
      sel = sel.replacingOccurrences(of: "'", with: "").trimmingCharacters(in: .whitespaces)
      runPython(
        "from mmol import inspector as _ai\n_ai.clear_rep_colors('\(sel.isEmpty ? "all" : sel)')"
      )
    }
  }

  /// True if any statement in `command` can change residues or their colours.
  static func changesSequence(_ command: String) -> Bool {
    let verbs: Set<String> = [
      "remove", "delete", "extract", "create", "alter",
      "alter_state", "color", "colour", "spectrum", "set_color",
      "cbc", "cnc", "fix_chemistry", "fab",
    ]
    for stmt in command.split(whereSeparator: { $0 == "\n" || $0 == ";" }) {
      let s = stmt.trimmingCharacters(in: .whitespaces).lowercased()
      if s.contains("util.") { return true }
      let verb = s.prefix { $0.isLetter || $0 == "_" }
      if verbs.contains(String(verb)) { return true }
    }
    return false
  }

  // Classify a command as a known long (>~2s) operation that deserves the
  // "Calculating…" overlay. Only fire-and-forget ops (no synchronous return
  // value the caller reads) are listed — export/PNG paths drive the busy flag
  // through their own runHeavy completion. nil = light/interactive.
  private func heavyLabel(for command: String) -> String? {
    let l = command.lowercased()
    if l.hasPrefix("ray ") || l == "ray" || l.hasPrefix("rayg ") || l == "rayg" {
      return "Ray tracing…"
    }
    if l.contains("show surface") || l.contains("as surface") { return "Calculating surface…" }
    if l.contains("show mesh") || l.contains("as mesh") { return "Calculating mesh…" }
    if l.contains("show dots") || l.contains("as dots") { return "Calculating dots…" }
    if l.contains("set surface_quality") || l.contains("set solvent_radius")
      || l.contains("set surface_carve")
    {
      return "Recomputing surface…"
    }
    return nil
  }

  // MARK: - Progress bar (heavy-op dispatch)

  // Run a long op with the progress bar at the bottom of the viewport. The
  // embedded core's GIL model (PAutoBlock) is NOT safe to call off the main
  // thread, so the op stays on the main thread and a background thread moves
  // the bar.
  func runHeavy(_ label: String, _ work: @escaping () -> Void) {
    if label.hasPrefix("Ray") {
      DispatchQueue.main.async {
        let bar = RayProgressBar.show()
        work()
        bar?.hide()
      }
      return
    }
    // `show surface` / `mesh` / `dots` (and surface-setting changes) only flag
    // the rep: PyMOL builds it lazily on the next Metal frame, the slow part,
    // reporting PyMOLBridge_BuildProgress. So the bar stays up until that frame
    // has rendered (heavyRenderTick), with a backstop if the loop ever stalls.
    PyMOLBridge_ResetBuildProgress()
    if buildBar == nil { buildBar = RayProgressBar.show(progress: PyMOLBridge_BuildProgress) }
    DispatchQueue.main.async { [weak self] in
      work()
      guard let self else { return }
      self.pendingHeavyClearFrames = 2
      self.busyBackstop?.cancel()
      let bs = DispatchWorkItem { [weak self] in
        self?.pendingHeavyClearFrames = 0
        self?.hideBuildBar()
      }
      self.busyBackstop = bs
      DispatchQueue.main.asyncAfter(deadline: .now() + 60, execute: bs)
    }
  }

  private func hideBuildBar() {
    buildBar?.hide()
    buildBar = nil
  }

  // Called once per completed Metal frame (from MetalViewport.draw, main
  // thread). Hides the build's progress bar after the frame that built the
  // rep. `presented` is false when the frame rendered but no drawable was
  // available, so nothing reached the screen.
  func heavyRenderTick(presented: Bool = true) {
    guard pendingHeavyClearFrames > 0 else { return }
    pendingHeavyClearFrames -= 1
    if pendingHeavyClearFrames == 0 {
      busyBackstop?.cancel()
      busyBackstop = nil
      hideBuildBar()
    }
  }

  // Whole-structure view fits (orient/reset, load/fetch auto_zoom) + showing a
  // surface/mesh/dots set an atom-fit slab; those reps add a ~solvent_radius
  // shell beyond the atoms, so the near plane slices the surface front (exposing
  // the interior). After such a command, re-fit the view with a buffer so the
  // whole surface stays visible (runs synchronously after the command — cmd.do
  // is blocking). .pse loads restore their own exact view → excluded. Bare zoom
  // / center are user-directed framing → not overridden.
  private func maybeWidenClipForSurface(for command: String) {
    let l = command.lowercased()
    if l.contains(".pse") { return }
    let triggers = [
      "orient", "reset", "load ", "fetch ",
      "show surface", "show mesh", "show dots",
    ]
    guard triggers.contains(where: { l.contains($0) }) else { return }
    runPython("from mmol import inspector as _ai\n_ai.widen_clip_for_surface()")
  }

  private func maybeCaptureRenderedPNG(_ command: String) -> Bool {
    let t = command.trimmingCharacters(in: .whitespaces)
    let lower = t.lowercased()
    guard lower == "png" || lower.hasPrefix("png ") else { return false }
    // Ray-traced export → let the core handle it (that path works).
    if lower.replacingOccurrences(of: " ", with: "").contains("ray=1") { return false }

    // Parse `png file[, width[, height[, ...]]]` — both positional and
    // keyword (width=, height=) forms. Filename is the first arg.
    let argStr = String(t.dropFirst(3)).trimmingCharacters(in: .whitespaces)
    let args = argStr.split(separator: ",").map {
      $0.trimmingCharacters(in: CharacterSet(charactersIn: " \"'"))
    }
    guard let first = args.first, !first.isEmpty else { return false }
    let path = (first as NSString).expandingTildeInPath

    // Resolve explicit width/height: positional args 2 & 3, or width=/height=.
    var w = 0
    var h = 0
    for (i, a) in args.enumerated() {
      let kv = a.split(separator: "=", maxSplits: 1).map(String.init)
      if kv.count == 2 {
        let key = kv[0].lowercased().trimmingCharacters(in: .whitespaces)
        let val = Int(Double(kv[1].trimmingCharacters(in: .whitespaces)) ?? 0)
        if key == "width" { w = val }
        if key == "height" { h = val }
      } else if i == 1 {
        w = Int(Double(a) ?? 0)
      } else if i == 2 {
        h = Int(Double(a) ?? 0)
      }
    }

    // Explicit resolution → GPU hi-res offscreen render. Otherwise capture
    // the live Metal frame at window resolution.
    if w > 0 && h > 0 {
      renderHiResPNG(path, width: w, height: h)
    } else {
      capturePNG(path)
    }
    return true
  }

  // Save the next rendered Metal frame to a PNG (captures at window resolution).
  func capturePNG(_ path: String) {
    guard let inst = instance else { return }
    PyMOLBridge_CapturePNG(inst, path)
  }

  // Render the full Metal pipeline offscreen at an arbitrary resolution and
  // write a PNG (Metal-accelerated export — all reps). Synchronous; runs on
  // the main thread.
  func renderHiResPNG(_ path: String, width: Int, height: Int) {
    guard let inst = instance else { return }
    PyMOLBridge_RenderHiResPNG(inst, path, Int32(width), Int32(height))
  }

  // Exporter for cmd.movie_export requests (#581). Owned here, not by a sheet,
  // so a scripted export outlives any UI and has no save panel: it writes
  // straight to the requested path and reports in the log.
  private lazy var scriptedMovieExporter = MovieExporter()

  /// Startup chatter (splash, ~/.pymolrc echoes) stays out of the log; once the
  /// startup script has run, the log starts with the mMol banner.
  private var startupQuiet = true
  private var startupErrors = false

  func finishStartupLog() {
    startupQuiet = false
    feedbackLog =
      AppBuildInfo.logo + ["", AppBuildInfo.consoleBanner]
      + (startupErrors
        ? ["  Note: ~/.pymolrc reported errors at startup (run @~/.pymolrc to see them)."] : [])
  }

  // Append a line to the console log from the main thread.
  func logLine(_ line: String) {
    feedbackLog.append(line)
  }

  // Handle a `MOVIEEXPORT:<json>` request emitted by cmd.movie_export.
  // MUST be called on the main thread.
  func startScriptedMovieExport(_ json: String) {
    let req: MovieExportRequest
    switch MovieExportRequest.decode(json) {
    case .failure(let e):
      logLine(" movie_export: \(e.message)")
      return
    case .success(let r): req = r
    }
    let exporter = scriptedMovieExporter
    guard !exporter.isExporting && !exportRenderActive else {
      logLine(" movie_export: another movie export is already running")
      return
    }
    // From the request, not the playback mirror (which lags by a poll).
    let frames = max(req.frames, 1)
    guard frames > 1 else {
      logLine(" movie_export: there is no movie to export (build one with Movie ▸ Edit Timeline)")
      return
    }
    let first = max(req.first, 1)
    let last = req.last > 0 ? min(req.last, frames) : frames
    guard first <= last else {
      logLine(" movie_export: first frame \(first) is past the last frame \(last)")
      return
    }
    let fps = max(req.options.fpsOverride, 1)
    let dest = URL(fileURLWithPath: req.path)
    let isFrames = req.options.format == .png
    // Check the destination BEFORE rendering, and never delete anything we
    // didn't make: a PNG sequence needs a new or empty folder; a movie file
    // may replace a file but never a folder.
    if let problem = Self.movieDestinationProblem(dest, frames: isFrames) {
      logLine(" movie_export: \(problem)")
      return
    }
    exporter.onFinish = { [weak self] url, error in
      guard let self = self else { return }
      guard let url = url else {
        self.logLine(" movie_export: failed: \(error ?? "unknown error")")
        return
      }
      do {
        let fm = FileManager.default
        try fm.createDirectory(
          at: dest.deletingLastPathComponent(),
          withIntermediateDirectories: true)
        // Re-check: the destination may have changed during the render.
        if let problem = Self.movieDestinationProblem(dest, frames: isFrames) {
          try? fm.removeItem(at: url)
          self.logLine(" movie_export: \(problem)")
          return
        }
        if fm.fileExists(atPath: dest.path) {
          // Only an empty folder (frames) or a file (movie) gets here.
          try fm.removeItem(at: dest)
        }
        try fm.moveItem(at: url, to: dest)
        self.logLine(" movie_export: wrote \(dest.path)")
      } catch {
        self.logLine(" movie_export: could not write \(dest.path): \(error.localizedDescription)")
      }
    }
    logLine(
      " movie_export: rendering \(last - first + 1) frames at "
        + "\(req.options.width)×\(req.options.height)…")
    exporter.start(
      engine: self, options: req.options,
      first: first, last: last, fps: fps)
  }

  // Why `dest` can't receive this export, or nil if it can.
  static func movieDestinationProblem(_ dest: URL, frames: Bool) -> String? {
    let fm = FileManager.default
    var isDir: ObjCBool = false
    guard fm.fileExists(atPath: dest.path, isDirectory: &isDir) else { return nil }
    if frames {
      guard isDir.boolValue else {
        return "\(dest.path) is a file; a PNG sequence needs a new or empty folder"
      }
      let contents = (try? fm.contentsOfDirectory(atPath: dest.path)) ?? ["?"]
      return contents.isEmpty
        ? nil
        : "\(dest.path) is not empty; choose a new or empty folder for the frames"
    }
    return isDir.boolValue ? "\(dest.path) is a folder; give a movie file name" : nil
  }

  // Rebuild dirty object representations for the current frame on the MAIN
  // thread. Movie export calls this before its off-main renderHiResPNG so the
  // rep rebuild — which reaches the Python C-API via the busy-status callback —
  // runs where the embedded interpreter's GIL is held, not on the render queue.
  // MUST be called on the main thread.
  func updateScene() {
    guard let inst = instance else { return }
    PyMOLBridge_UpdateScene(inst)
  }

  // Apply a Metal-renderer letterbox so a loaded .pse reproduces its saved
  // viewport aspect. Loading non-session content (or reinitialize) clears it.
  private func handleSessionViewport(for command: String) {
    let lower = command.lowercased()
    if lower.contains(".pse"), let path = sessionPath(from: command) {
      // grid_mode (objects laid out in side-by-side cells) is now rendered
      // natively by the Metal renderer (per-slot viewport loop in
      // SceneRenderMetal), so the session keeps whatever grid_mode it saved
      // — no reset needed. See GitHub issue: grid-mode Metal support.
      // Read the session's 'main' [W,H] (→ letterbox via SESSIONVP) AND fix
      // the camera: modern .pse files store a 25-float SceneViewType, but
      // our embedded core is 18-float and mis-restores it (front/back
      // flip). Convert the stored col-major 4x4 rotation into the 18-float
      // set_view's 9-float rotation, which is ALSO column-major (same as
      // get_view — do NOT transpose, or the restored camera orientation
      // differs from the saved one). Take the 4x4's upper-left 3x3 in the
      // same column-major order. (18-float sessions restore natively.)
      runPython(
        "import pickle, os\n"
          + "from pymol import cmd as _sc\n"
          + "try:\n"
          + "    _d = pickle.load(open(os.path.expanduser(r'''\(path)'''), 'rb'))\n"
          + "    _m = _d.get('main')\n"
          + "    if _m and _m[0] and _m[1]: print('SESSIONVP:%d,%d' % (int(_m[0]), int(_m[1])))\n"
          + "    _v = _d.get('view')\n"
          + "    if _v and len(_v) >= 25:\n"
          + "        _R = [0.0]*9\n"
          + "        for _i in range(3):\n"
          + "            for _j in range(3):\n"
          + "                _R[_i*3+_j] = _v[_i*4+_j]\n"
          + "        _sc.set_view(_R + list(_v[16:19]) + list(_v[19:22]) + list(_v[22:24]) + [_v[24]])\n"
          + "except Exception:\n"
          + "    pass")
      // A .pse saved while in Move mode bakes in the transient move gizmo
      // (the "_move_gizmo" CGO). interactionMode is NOT persisted, so on
      // cold-launch restore / manual open we'd otherwise show a stray gizmo
      // with no active mode (and its "_"-prefixed name hides it from the
      // object panel, so it can't be removed). Strip it unless we're
      // currently in Move mode — the guard avoids deleting a live gizmo
      // when a file is opened mid-move. cleanup() is idempotent.
      if interactionMode != .move {
        runPython("from mmol import move as _mm\n_mm.cleanup()")
      }
    } else if lower.hasPrefix("load ") || lower.hasPrefix("fetch ")
      || lower.hasPrefix("reinitialize")
    {
      setLetterboxAspect(0)  // new non-session content → fill the window
    }
  }

  /// The console hint for the text at the prompt (modules/mmol/console_hint.py).
  func consoleHint(_ text: String) -> ConsoleHint? {
    guard isReady, !text.isEmpty, let c = PyMOLBridge_Hint(text) else { return nil }
    defer { PyMOLBridge_FreeFeedback(c) }
    return try? JSONDecoder().decode(
      ConsoleHint.self, from: Data(bytesNoCopy: c, count: strlen(c), deallocator: .none))
  }

  // Extract the file path from a `load <path>[, ...]` command (best effort).
  private func sessionPath(from command: String) -> String? {
    let t = command.trimmingCharacters(in: .whitespaces)
    guard let r = t.range(of: "load ", options: [.caseInsensitive]) else { return nil }
    var rest = String(t[r.upperBound...])
    if let comma = rest.firstIndex(of: ",") { rest = String(rest[..<comma]) }
    rest = rest.trimmingCharacters(in: .whitespaces)
    return rest.isEmpty ? nil : rest
  }

  func setLetterboxAspect(_ aspect: Float) {
    guard let inst = instance else { return }
    PyMOLBridge_SetLetterboxAspect(inst, aspect)
  }

  // Debug: run raw Python in the embedded interpreter.
  #if DEBUG
    /// Test seam: observes every Python string this engine emits. Placed BEFORE
    /// the isReady guard so it reports what the engine *intended* to run,
    /// independent of whether the core is initialized — that keeps tests of
    /// teardown behaviour (metal_move.cleanup / appkit_measure.reset, which
    /// leave no Swift-side trace) working under any host configuration.
    /// Compiled out of Release; cost in Debug is one optional-closure call.
    var pythonTap: ((String) -> Void)? = nil
  #endif

  func runPython(_ code: String) {
    #if DEBUG
      pythonTap?(code)
    #endif
    guard isReady else { return }
    PyMOLBridge_RunPython(code)
  }

  /// Read-only Python for the periodic UI polls: runs like `runPython` but
  /// does not leave a redisplay behind. The polls' temporary selections
  /// (count_atoms / count_states) make the core mark the scene changed, which
  /// re-rendered a full frame on every poll — 2-3 frames/s on a static view,
  /// 20-30 ms of GPU each with ray tracing on. Use ONLY for code that changes
  /// nothing visible.
  func runPythonQuiet(_ code: String) {
    #if DEBUG
      pythonTap?(code)
    #endif
    guard isReady, let inst = instance else { return }
    PyMOLBridge_RunPythonQuiet(inst, code)
  }

  /// Two-stage "clear selection" used by both the selection-chip X and the Esc
  /// key (issues #163 + #166). Stage 1: if any public selection is enabled
  /// (its pink markers are showing), `deselect` — hide the highlight but keep
  /// the named selection so a follow-up refinement is still possible. Stage 2:
  /// if nothing is enabled but a `sele` named selection still exists, `delete`
  /// it outright so the object list is clean. The logic runs INLINE in Python
  /// (one runPython round-trip) rather than in Swift so it reads the CURRENT
  /// core state directly — the Swift `objects` mirror is only refreshed by the
  /// ~500ms feedback poll, so a Swift-side branch would race the poll lag.
  func escapeClearSelection() {
    guard isReady else { return }
    var py = "from pymol import cmd as _c\n"
    py += "if _c.get_names('public_selections', enabled_only=1):\n"
    py += "    _c.deselect()\n"
    py += "elif 'sele' in (_c.get_names('public_selections') or []):\n"
    py += "    _c.delete('sele')\n"
    runPython(py)
  }

  /// Write the whole session to `url` (a .pse) via cmd.save (the only path to the
  /// C++ saver is runPython — there's no Swift cmd.save wrapper) and track it as
  /// the open document so a subsequent ⌘S overwrites it with no panel. Raw triple
  /// quotes tolerate spaces/quotes in the path.
  /// NOTE: writes a plain URL — works for the Developer-ID build and for any
  /// Save-As-chosen URL (auto write-granted). A sandboxed/MAS build that silently
  /// overwrites a Finder-opened file would additionally need a security-scoped
  /// bookmark resolved here.
  /// TODO: security-scoped bookmark for sandboxed overwrite.
  func saveSession(to url: URL) {
    runPython("from pymol import cmd as _c\n_c.save(r'''\(url.path)''')")
    currentSessionURL = url
  }

  // MARK: - Classic defaults

  /// The fixed look the app runs with (desktop PyMOL 2.5): ~/.pymolrc owns the
  /// scene (background, colours, how new objects show); this only sets the
  /// selection colours and switches the real-time Metal effects off — even over
  /// a restored session that saved them on — so the viewport rasterises like
  /// desktop PyMOL until `ray`.
  /// A slider's in-between values: applied without echoing to the console
  /// (the release runs the logged `set`, so the log shows the value landed on).
  func setSilently(_ setting: String, _ value: String, _ object: String = "") {
    let obj = object.replacingOccurrences(of: "'", with: "")
    runPython("from pymol import cmd as _c\n_c.set('\(setting)', \(value), '\(obj)', quiet=1)")
  }

  func applyClassicDefaults() {
    guard isReady else { return }
    if let inst = instance {
      PyMOLBridge_SetSelectionColor(inst, 1.0, 0.3, 0.3)
      // Hover pre-selection preview (issue #165): light cyan, distinct
      // from the committed selection colour.
      PyMOLBridge_SetPreselectionColor(inst, 0.40, 0.85, 1.0)
    }
    // active_selections off: turning one selection on doesn't hide the others,
    // so several named selections can show at once.
    runPython(
      "from pymol import cmd as _c\n_c.set('metal_outline', 0)\n_c.set('metal_ssao', 0)\n"
        + "_c.set('active_selections', 0)\n")
  }

  /// True when `path` is a session file (.pse/.psw). A session restores its own
  /// object names, so the coordinate-file load contract (object name) does not
  /// apply to it.
  static func isSessionFile(_ path: String) -> Bool {
    let ext = (path as NSString).pathExtension.lowercased()
    return ext == "pse" || ext == "psw"
  }

  /// The exact command `loadStructure` runs — factored pure for tests.
  /// Sessions load bare: passing an object name to `load foo.pse` is meaningless
  /// (nothing by that name is created, and the selector logs `Invalid selection
  /// name` straight to the console — issue #272).
  static func loadInvocation(path: String, name: String) -> String {
    isSessionFile(path) ? "load \(path)" : "load \(path), \(name)"
  }

  /// Load a structure (or restore a session). Routes all UI/agent/open-with
  /// load paths.
  func loadStructure(path: String, name: String) {
    guard isReady else { return }
    runCommand(Self.loadInvocation(path: path, name: name))
  }

  func fetchStructure(id: String) {
    guard isReady else { return }
    let clean = id.replacingOccurrences(of: "'", with: "")
    runCommand("fetch \(clean), async=0, type=pdb")
  }

  /// Clear the whole session (File ▸ Clear Session): wipe all objects, selections,
  /// and camera via `reinitialize`, then re-apply the classic defaults. The
  /// objects panel refreshes on the next poll tick.
  func clearSession() {
    guard isReady else { return }
    // No document is open after a clear — the next ⌘S becomes a Save As.
    currentSessionURL = nil
    // reinitialize wipes the core movie (mset/mview); drop the camera track's
    // keyframes so the Timeline doesn't show diamonds for frames that no
    // longer exist.
    cameraKeyframes.removeAll()
    sceneMarkers.removeAll()
    runCommand("reinitialize")
    // reinitialize also resets engine settings to defaults — restore the
    // fetch_path that initialize() set (the writable temp dir) so a post-clear
    // `fetch` doesn't fall back to cwd (which can fail or prompt for file
    // access), then re-apply the classic defaults.
    runPython("from pymol import cmd as _c; _c.set('fetch_path', '\(NSTemporaryDirectory())')")
    applyClassicDefaults()
  }

  /// Reset the Effects-group post-processing/stylization settings to their
  /// SettingInfo.h defaults, for the Inspector's "Reset effects" button.
  /// Non-destructive — keeps loaded structures.
  func resetEffects() {
    let defaults: [(String, String)] = [
      ("metal_outline", "0"),
      ("metal_outline_width", "1.4"),
      ("metal_outline_color", "0x000000"),
    ]
    for (k, v) in defaults { runCommand("set \(k), \(v)") }
  }

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
  @Published var playingObjects: Set<String> = []
  @Published var objectFPS: [String: Double] = [:]
  private var objectStateTimers: [String: Timer] = [:]

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

  // Recall a saved scene (camera + reps + visibility). Injection-safe (base64).
  func recallScene(_ name: String) {
    guard !name.isEmpty else { return }
    let b64 = Data(name.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom pymol import cmd as _c\n"
        + "_c.scene(_b64.b64decode('\(b64)').decode('utf-8'), 'recall')")
  }

  // Per-scene management for the scene chips' long-press menu. Injection-safe.
  // updateScene overwrites the named scene with the CURRENT view/reps.
  func updateScene(_ name: String) { sceneAction(name, "update") }
  func deleteScene(_ name: String) { sceneAction(name, "delete") }

  private func sceneAction(_ name: String, _ action: String) {
    guard !name.isEmpty else { return }
    let b64 = Data(name.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom pymol import cmd as _c\n"
        + "_c.scene(_b64.b64decode('\(b64)').decode('utf-8'), '\(action)')")
  }

  func renameScene(_ name: String, to newName: String) {
    let n = newName.trimmingCharacters(in: .whitespacesAndNewlines)
    guard !name.isEmpty, !n.isEmpty, n != name else { return }
    let b = Data(name.utf8).base64EncodedString()
    let nb = Data(n.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom pymol import cmd as _c\n"
        + "_c.scene(_b64.b64decode('\(b)').decode('utf-8'), 'rename', "
        + "new_key=_b64.b64decode('\(nb)').decode('utf-8'))")
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
  // spec (camera → UUID into _views; scene → base64 name) with each item's
  // inbound-transition easing, base64-wrapped so scene names / JSON can't break
  // the Python literal. Mirrors playback.frameCount immediately (the ~10/s poll
  // would otherwise lag the ruler a beat).
  func rebuildMovie() {
    guard isReady else { return }
    if timelineItems.isEmpty {
      // "[]" → clears the movie (see appkit_movie.rebuild).
      runPython("from mmol import movie as _am\n_am.rebuild('W10=')")
      playback.frameCount = 1
      playback.currentFrame = 1
      return
    }
    let spans = itemSpans()
    var parts: [String] = []
    for (i, it) in timelineItems.enumerated() {
      let (start, end) = spans[i]
      let ease =
        it.transition.linear
        ? "\"power\":1.0,\"linear\":1"
        : "\"power\":0.0,\"linear\":0"
      switch it.kind {
      case .camera:
        var s = "{\"frame\":\(start),\"cam\":\"\(it.id.uuidString)\",\(ease)"
        if let ps = it.pinnedState { s += ",\"state\":\(ps)" }
        parts.append(s + "}")
      case .scene(let name):
        let b64 = Data(name.utf8).base64EncodedString()
        var s = "{\"frame\":\(start),\"scene\":\"\(b64)\",\(ease)"
        if let ps = it.pinnedState { s += ",\"state\":\(ps)" }
        parts.append(s + "}")
      case .states(let spec):
        let objs =
          spec.objects.map {
            "[" + $0.map { "\"\($0)\"" }.joined(separator: ",") + "]"
          } ?? "null"
        parts.append(
          "{\"frame\":\(start),\"end\":\(end),\"states\":1,"
            + "\"objects\":\(objs),\"mode\":\"\(spec.mode.rawValue)\","
            + "\"first\":\(spec.firstModel),\"last\":\(spec.lastModel)}")
      }
    }
    let json = "[" + parts.joined(separator: ",") + "]"
    let jb64 = Data(json.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom mmol import movie as _am\n"
        + "_am.rebuild(_b64.b64decode('\(jb64)').decode('utf-8'))")
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

  func frameForItem(_ id: UUID) -> Int? {
    guard let idx = timelineItems.firstIndex(where: { $0.id == id }) else { return nil }
    let frames = itemFrames()
    return frames.indices.contains(idx) ? frames[idx] : nil
  }

  func seekToItem(_ id: UUID) {
    if let f = frameForItem(id) { seek(to: f) }
  }

  // Append a camera Roll/Rock as a run of 4 camera items joined by equal
  // transitions that sum to `duration`. Python computes the waypoint views
  // (spin/rock around the current view) under Swift-generated UUIDs.
  func appendCameraTemplate(kind: String, duration: Double, axis: String, angle: Double) {
    guard isReady else { return }
    let count = 4
    let ids = (0..<count).map { _ in UUID() }
    let idsJSON = "[" + ids.map { "\"\($0.uuidString)\"" }.joined(separator: ",") + "]"
    let idsB64 = Data(idsJSON.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom mmol import movie as _am\n"
        + "_am.capture_template_views(_b64.b64decode('\(idsB64)').decode('utf-8'), "
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

  // MARK: - Selection builder support

  // Preview how many atoms a selection expression would match, without
  // creating a lasting selection. Emits SELPREVIEW:<n> (or :err). The caller
  // debounces. Uses a throwaway '_pv' selection.
  func previewSelection(_ expr: String) {
    // Pass the expression as base64 and decode in Python, so a selection
    // string containing quotes / ''' / backslashes can neither break the
    // literal nor inject code (the old strip-backslashes + r'''…''' was both
    // lossy and injectable via an embedded ''').
    let b64 = Data(expr.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\n"
        + "from pymol import cmd as _c\n"
        + "_e = _b64.b64decode('\(b64)').decode('utf-8')\n"
        + "try:\n"
        + "    _n = _c.select('_pv', _e)\n"
        + "    print('SELPREVIEW:%d' % int(_n))\n"
        + "    _c.delete('_pv')\n"
        + "except Exception:\n"
        + "    print('SELPREVIEW:err')")
  }

  // Create/overwrite a named selection and enable it.
  func createSelection(name: String, expr: String) {
    // base64 the name + expression (see previewSelection): injection-safe.
    let nb = Data(name.utf8).base64EncodedString()
    let eb = Data(expr.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom pymol import cmd as _c\n"
        + "_c.select(_b64.b64decode('\(nb)').decode('utf-8'), "
        + "_b64.b64decode('\(eb)').decode('utf-8'), enable=1)")
  }

  // Rename an object/selection.
  func renameObject(_ old: String, to newName: String) {
    guard !newName.isEmpty else { return }
    // base64 both names (see previewSelection): injection-safe. The C core
    // sanitizes the decoded name via ObjectMakeValidName.
    let ob = Data(old.utf8).base64EncodedString()
    let nb = Data(newName.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom pymol import cmd as _c\n"
        + "_c.set_name(_b64.b64decode('\(ob)').decode('utf-8'), "
        + "_b64.b64decode('\(nb)').decode('utf-8'))")
  }

  // MARK: - Interactive measurement

  func setMeasureMode(_ k: MeasureKind?) {
    measureMode = k
    if let k = k {
      leaveViewportTool()  // mutually exclusive
      runPython("from mmol import measure as _am\n_am.set_mode('\(k.rawValue)')")
    } else {
      runPython("from mmol import measure as _am\n_am.reset()")
    }
  }

  // A tap while in measure mode: accumulate an atom pick (NDC in [-1,1]).
  func measurePick(ndcX: Float, ndcY: Float, aspect: Float) {
    runPython("from mmol import measure as _am\n_am.pick(\(ndcX), \(ndcY), \(aspect))")
  }

  func clearMeasurements() {
    runPython("from mmol import measure as _am\n_am.clear_all()")
  }

  /// Leave the viewport tool a docked-bar mode cannot share the screen with: Move
  /// (its gizmo takes the drags). Used by the Measure setter on entry. Box Select
  /// is left alone.
  private func leaveViewportTool() {
    if interactionMode == .move { setInteractionMode(.viewing) }
  }

  // MARK: - Exclusive interaction modes
  //         (Move / Box Select / Measure)

  /// Leave whichever exclusive interaction mode is active, each through that
  /// mode's OWN existing exit path — so the Esc key (see
  /// ContentView.installEscKeyMonitor) and the overlays' ✕ buttons converge on
  /// identical behavior; Esc must never grow a second, divergent path. (#235)
  ///
  /// The modes are mutually exclusive, so in practice at most one branch runs. They
  /// are checked independently rather than with early returns so a desynchronized
  /// state still unwinds completely instead of leaving one mode stranded.
  ///
  /// EVERY exclusive mode must appear here: Move, Box Select and Measure. A
  /// mode that enters but is not listed is one Esc cannot leave.
  ///
  /// - Returns: whether any mode was actually exited. Callers use this to fall
  ///   through to their next behavior (Esc clears the selection) when Escape
  ///   found no mode to leave.
  @discardableResult
  func exitActiveInteractionMode() -> Bool {
    var exited = false
    if interactionMode == .move {
      setInteractionMode(.viewing)
      exited = exited || interactionMode != .move
    }
    if interactionMode == .boxSelect {
      setInteractionMode(.viewing)
      exited = exited || interactionMode != .boxSelect
    }
    if measureMode != nil {
      setMeasureMode(nil)
      exited = exited || measureMode == nil
    }
    return exited
  }

  // MARK: - Box Select (rubber-band selection, #358)

  // Leading-edge throttle for the live commit, mirroring hoverPreview: a drag
  // produces a new rectangle every frame, and each commit projects EVERY drawn
  // atom, so an unthrottled drag would queue Python faster than it can run.
  private var boxWork: DispatchWorkItem?
  private var lastBoxFire = Date.distantPast
  private let kBoxCommitInterval = 0.045  // <= ~22 commits/s
  // Whether a box session is open (box_begin has run, so the pre-box snapshot
  // exists). Tracked so leaving the mode only pays for the teardown round-trip
  // when there is something to tear down — setInteractionMode runs on every
  // Move/Measure toggle too.
  private var boxSessionActive = false

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
    boxWork?.cancel()
    let fire: () -> Void = { [weak self] in
      guard let self = self else { return }
      self.lastBoxFire = Date()
      // mode='add': the tool only ever grows the selection it opened on.
      self.runPython(
        "from mmol import pick as _mp\n"
          + "_mp.box_commit_ndc(\(rect.minX), \(rect.minY), \(rect.maxX), \(rect.maxY), "
          + "\(aspect), name='sele', mode='add')")
    }
    let elapsed = Date().timeIntervalSince(lastBoxFire)
    if elapsed >= kBoxCommitInterval {
      fire()  // leading edge — track the drag
    } else {
      let work = DispatchWorkItem(block: fire)  // trailing — final rest position
      boxWork = work
      DispatchQueue.main.asyncAfter(
        deadline: .now() + (kBoxCommitInterval - elapsed), execute: work)
    }
  }

  /// Forget the rectangle and close the session. The last commit STAYS — the
  /// tool has no cancel, because with commit-on-drag there is no uncommitted
  /// state left to throw away.
  private func clearBoxState() {
    boxWork?.cancel()
    boxWork = nil
    boxDrag = nil
    boxRect = nil
    // Reset the throttle clock too, so the FIRST rectangle of the next box
    // always commits on the leading edge instead of waiting out a window
    // left over from the previous one.
    lastBoxFire = .distantPast
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
    let nb = Data(name.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom mmol import move as _mm\n"
        + "_mm.set_active(_b64.b64decode('\(nb)').decode('utf-8'), \(gizmoAspect))")
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

  // MARK: - Settings panel

  func loadSettingsCatalog() {
    runPython("from mmol import settings as _as\n_as.catalog()")
  }

  func setSetting(_ name: String, _ value: String) {
    // base64 the name + value (see previewSelection): injection-safe.
    let nb = Data(name.utf8).base64EncodedString()
    let vb = Data(value.utf8).base64EncodedString()
    runPython(
      "import base64 as _b64\nfrom mmol import settings as _as\n"
        + "_as.set_value(_b64.b64decode('\(nb)').decode('utf-8'), "
        + "_b64.b64decode('\(vb)').decode('utf-8'))")
  }

  // Read the settings catalog JSON the bridge wrote to the temp dir.
  private func loadSettingsCatalogFile() {
    let path = TempChannel.path(TempChannel.Stem.settings)
    guard let data = FileManager.default.contents(atPath: path),
      let items = try? JSONDecoder().decode([SettingItem].self, from: data)
    else { return }
    DispatchQueue.main.async { self.settingsCatalog = items }
  }

  // SETVAL:<name>=<val> — update one row after an edit (reflects clamping/parse).
  private func parseSetValFeedback(_ line: String) {
    let body = String(line.dropFirst("SETVAL:".count))
    guard let eq = body.firstIndex(of: "=") else { return }
    let name = String(body[..<eq])
    let val = String(body[body.index(after: eq)...])
    DispatchQueue.main.async {
      if let i = self.settingsCatalog.firstIndex(where: { $0.name == name }) {
        self.settingsCatalog[i].val = val
      }
    }
  }

  // Parse MEASURE:<json> {kind,count,need,value?} → measureStatus text.
  func parseMeasureFeedback(_ line: String) {
    let js = String(line.dropFirst("MEASURE:".count))
    guard let data = js.data(using: .utf8),
      let root = try? JSONSerialization.jsonObject(with: data) as? [String: Any]
    else { return }
    let kind = root["kind"] as? String ?? "distance"
    let count = (root["count"] as? NSNumber)?.intValue ?? 0
    let need = (root["need"] as? NSNumber)?.intValue ?? 2
    let value = (root["value"] as? NSNumber)?.doubleValue
    let unit = kind == "distance" ? " Å" : "°"
    let status: String
    if let v = value {
      let s =
        (v == v.rounded())
        ? String(Int(v)) : String(format: kind == "distance" ? "%.2f" : "%.1f", v)
      status = "\(kind.capitalized): \(s)\(unit)  ·  tap to start another"
    } else if count == 0 {
      status = "Tap \(need) atoms"
    } else {
      status = "Tap \(need - count) more atom\(need - count == 1 ? "" : "s")"
    }
    DispatchQueue.main.async { self.measureStatus = status }
  }

  // Fetch per-object residue sequences (one guide atom per residue) and emit
  // a SEQPANEL: feedback line parsed into `sequences`. Run via raw Python so
  // the command isn't echoed to the log.
  func fetchSequences() {
    guard isReady else { return }
    // Write the (potentially large) per-residue JSON to a temp file rather
    // than the feedback buffer — PyMOL feedback lines are capped (~1KB) and
    // long sequences would be split/truncated. Print only a short marker;
    // Swift reads the file (same process TMPDIR) on seeing it.
    // Each residue carries its guide-atom color index; a color table maps
    // index -> RGB so the panel can reflect the real molecular coloring.
    // The query (incl. BIMO-style alignment gap layout) lives in the bundled
    // mmol.sequence module so it stays readable and testable.
    runPythonQuiet(
      "from mmol import sequence as _sq\n"
        + "_sq.poll()\n"
        + "print('SEQPANEL:ready')"
    )
  }

  // Poll which residues are in the active 'sele' selection so the sequence
  // panel can highlight them (3D -> sequence sync). Writes compact keys
  // ("obj/chain/resi") to a temp file; emits a short SEQSEL marker.
  /// `colorCheck` adds the SEQCOLOR fingerprint (a scan of every drawn atom's
  /// colour). Only the half-second poll asks for it; a click's refresh
  /// shouldn't pay for it.
  func fetchSequenceSelection(colorCheck: Bool = false) {
    guard isReady else { return }
    runPythonQuiet(
      colorCheck ? sequenceSelectionPython + sequenceColorPython : sequenceSelectionPython)
  }

  private var sequenceSelectionPython: String {
    "from mmol import sequence as _sq\n_sq.poll_selection()\n"
  }

  /// Fingerprint of the colours (and count) of the atoms the strip draws: a
  /// recolour from anywhere — an alias, a `run` script, a menu — or atoms
  /// removed by a script changes it, and the strip refetches.
  private var sequenceColorPython: String {
    "from pymol import cmd as _sc\nprint('SEQCOLOR:%d' % _sc._cmd.color_fingerprint(_sc._COb))\n"
  }

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

  // Trailing coalescer for the hover re-pick so a fast mouse sweep collapses to
  // ~one pick every ~33 ms (mirrors MetalViewport.armPanEndDebounce). Plus a
  // last-NDC gate so a barely-moved pointer doesn't re-pick at all — the actual
  // projection over all atoms is the cost we're throttling.
  private var hoverWork: DispatchWorkItem?
  private var lastHoverNDC: (Float, Float)? = nil
  private var lastHoverFire = Date.distantPast
  private let kHoverInterval = 0.045  // leading-edge throttle: ≤~22 picks/s
  private let kHoverMinNDC: Float = 0.004  // sub-pixel move gate (NDC²-ish)

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
    hoverWork?.cancel()
    let fire: () -> Void = { [weak self] in
      guard let self = self else { return }
      // A trailing pick scheduled before `ray` must not run over the held
      // ray image (desktop PyMOL has no hover pick).
      if self.rayImageHeld { return }
      self.lastHoverFire = Date()
      self.runPython(
        "from mmol import pick as _mp; "
          + "_mp.hover_preview_at(\(ndcX), \(ndcY), \(aspect))")
    }
    let elapsed = Date().timeIntervalSince(lastHoverFire)
    if elapsed >= kHoverInterval {
      fire()  // leading edge — update now
    } else {
      let work = DispatchWorkItem(block: fire)  // trailing — final rest position
      hoverWork = work
      DispatchQueue.main.asyncAfter(
        deadline: .now() + (kHoverInterval - elapsed), execute: work)
    }
  }

  /// Cancel any pending hover pick and empty the preview selection. Cheap and
  /// idempotent — safe to call from mouse-down / drag-start / mouse-exit and
  /// when the user toggles the feature off.
  func clearHoverPreview() {
    hoverWork?.cancel()
    hoverWork = nil
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
    let path = (NSTemporaryDirectory() as NSString).appendingPathComponent("pymol_longpress.json")
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

  // MARK: - Render loop hooks (called by MetalViewport)

  func idle() {
    guard let inst = instance else { return }
    _ = PyMOLBridge_Idle(inst)
  }

  // Posted to force the MetalViewport off its on-demand gate for one frame.
  // The coordinator observes this (same handler as wake) and renders once
  // unconditionally. Needed when a command changes what's displayed but the
  // core's redisplay flag doesn't reach the next draw(in:) — notably a
  // per-object `set state` / movie `frame` scrub AFTER a movie exists, where
  // the flag is otherwise consumed before the gate checks it (issue #132).
  static let forceRedrawNotification = Notification.Name("PyMOLEngine.forceRedraw")

  // Request a single unconditional viewport repaint. Cheap (one extra frame);
  // safe during mplay (that path already renders every tick, so an extra
  // forced frame is a no-op there). Main-thread only — the coordinator kicks
  // an immediate draw() which must run on the UI thread.
  func requestViewportRedraw() {
    if Thread.isMainThread {
      NotificationCenter.default.post(name: PyMOLEngine.forceRedrawNotification, object: nil)
    } else {
      DispatchQueue.main.async {
        NotificationCenter.default.post(name: PyMOLEngine.forceRedrawNotification, object: nil)
      }
    }
  }

  func reshape(width: Int, height: Int) {
    guard let inst = instance else { return }
    PyMOLBridge_Reshape(inst, Int32(width), Int32(height))
  }

  #if DEBUG
    /// What the viewport hands the core for a pointer: a PyMOL button or
    /// drag event (the camera path), or a click's atom pick.
    enum ViewportInputEvent: Equatable {
      case button(button: Int32, state: Int32, x: Int32, y: Int32, modifiers: Int32)
      case drag(x: Int32, y: Int32, modifiers: Int32)
      case pick(ndcX: Float, ndcY: Float)
    }

    /// Test seam (#622's routing tests): while set, `button`, `drag` and
    /// `pick` report here INSTEAD of reaching the core, so a test can drive
    /// MetalViewport's handlers and see what reached the camera path without
    /// moving the shared engine's camera (the core runs queued mouse input
    /// in the next rendered frame). Compiled out of Release.
    var viewportInputTap: ((ViewportInputEvent) -> Void)? = nil
  #endif

  func button(_ btn: Int32, state: Int32, x: Int32, y: Int32, modifiers: Int32) {
    #if DEBUG
      if let viewportInputTap {
        viewportInputTap(.button(button: btn, state: state, x: x, y: y, modifiers: modifiers))
        return
      }
    #endif
    guard let inst = instance else { return }
    PyMOLBridge_Button(inst, btn, state, x, y, modifiers)
  }

  func drag(x: Int32, y: Int32, modifiers: Int32) {
    #if DEBUG
      if let viewportInputTap {
        viewportInputTap(.drag(x: x, y: y, modifiers: modifiers))
        return
      }
    #endif
    guard let inst = instance else { return }
    PyMOLBridge_Drag(inst, x, y, modifiers)
  }

  // Explicit camera-dolly zoom: frac > 0 zooms IN, frac < 0 zooms OUT.
  // Used by the pinch/scroll gestures instead of the scroll-wheel BUTTON path,
  // because PyMOL's default three_button_viewing binds the bare wheel to
  // 'slab' (clip), not zoom. The dolly distance is scaled by the current slab
  // depth (front+back)/2 so the feel is scene-size-independent, mirroring
  // cButModeZoomForward. get_view() here is the embedded 18-float layout
  // (front=v[15], back=v[16]). move('z', +d) moves the camera in → zoom in.
  func zoomBy(_ frac: Float) {
    guard isReady, frac.isFinite, abs(frac) > 1e-5 else { return }
    runPython(
      "from pymol import cmd as _zc\n"
        + "_zv = _zc.get_view()\n"
        + "_zc.move('z', (_zv[15] + _zv[16]) * 0.5 * (\(frac)))")
  }

  // Absolute zoom: set the apparent magnification `mag` (independent of the Lens
  // control). Target camera distance = radius / (mag * tan(fov/2)); dolly the
  // camera there via move('z', current - target). `radius`/`fovDeg` come from the
  // scene poll (same values used to display the slider), so the target distance
  // is consistent with the shown magnification. move z + = camera moves in.
  func setZoomMagnification(_ mag: Double, radius: Double, fovDeg: Double) {
    guard isReady, radius > 0, mag > 0 else { return }
    let t = tan(fovDeg * .pi / 360.0)
    guard t > 1e-6 else { return }
    let target = radius / (mag * t)
    guard target.isFinite, target > 0.01 else { return }
    runPython(
      "from pymol import cmd as _zc\n"
        + "_zv = _zc.get_view()\n"
        + "_cur = abs(_zv[11])\n"
        + "_zc.move('z', _cur - \(target))")
  }

  /// Fire a cmd.set_key binding by canonical PyMOL key token, returning
  /// whether one existed. The caller consumes the key event iff this is true,
  /// which is what lets an unbound ⌃M fall through to the Move menu item
  /// while a user-bound ⌃M shadows it (#258).
  ///
  /// Synchronous and in-process (same model as `complete(_:)`), so the answer
  /// is available while the NSEvent monitor still has to decide.
  func invokeKeyBinding(_ token: String) -> Bool {
    guard isReady, !token.isEmpty else { return false }
    return PyMOLBridge_InvokeKey(token) == 1
  }

  // Construct the Metal renderer from the MTKView (idempotent in the bridge).
  func setupMetalRenderer(view: MTKView) {
    guard let inst = instance else { return }
    PyMOLBridge_SetupMetalRenderer(inst, Unmanaged.passUnretained(view).toOpaque())
  }

  // Render one live frame at width x height (backing pixels).
  //
  // The view — not a drawable — is what gets handed over: the bridge acquires
  // the CAMetalDrawable itself, after the scene has been encoded, since only
  // the final post pass writes to it (#396). Returns false when no drawable
  // arrived, i.e. nothing was presented and the caller should render again.
  @discardableResult
  func renderMetalFrame(view: MTKView, width: Int, height: Int) -> Bool {
    guard let inst = instance else { return false }
    let vPtr = Unmanaged.passUnretained(view).toOpaque()
    return PyMOLBridge_RenderMetalFrame(inst, vPtr, Int32(width), Int32(height)) != 0
  }

  /// Live Metal frames committed but not yet completed on the GPU. The render
  /// loop caps this (RenderGate.maxFramesInFlight) so it skips a tick instead
  /// of queueing frames the GPU cannot retire.
  var metalFramesInFlight: Int {
    guard let inst = instance else { return 0 }
    return Int(PyMOLBridge_MetalFramesInFlight(inst))
  }

  // MARK: - Ray image

  /// The image the last `ray` produced (mmol/classic_ray.py), shown over the
  /// viewport until the core drops it.
  struct ClassicRayImage: Equatable {
    let image: NSImage
    let pixelWidth: Int
    let pixelHeight: Int
    let background: RGBA
  }
  @Published var classicRayImage: ClassicRayImage? = nil

  /// RAYIMG:<w>:<h>:<r>,<g>,<b>:<png path>
  private func parseClassicRayImage(_ line: String) {
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
    let quoted = path.replacingOccurrences(of: "\\", with: "\\\\")
      .replacingOccurrences(of: "'", with: "\\'")
    runPython("from pymol import cmd as _c\n_c.png('\(quoted)', prior=1, quiet=1)\n")
    return FileManager.default.fileExists(atPath: path)
  }

  /// Route `ray` to the desktop-PyMOL display path (mmol/classic_ray.py).
  func enableClassicRay() {
    runPythonQuiet("from mmol import classic_ray\nclassic_ray.enable(True)\n")
  }

  // MARK: - MCP view requests (mmol.mcp/tools.py)

  /// HEATMAP:<b64 json> → attach + configure a CSV heatmap (as C ▸ by file…);
  /// CONTACTMAP:<b64 json> → open that object's contact map window.
  private func handleMCPViewRequest(_ line: String) {
    guard let colon = line.firstIndex(of: ":"),
      let data = Data(base64Encoded: String(line[line.index(after: colon)...])),
      let req = (try? JSONSerialization.jsonObject(with: data)) as? [String: Any],
      let object = req["object"] as? String
    else { return }
    if line.hasPrefix("CONTACTMAP:") {
      DispatchQueue.main.async { self.contactMapRequest = object }
      return
    }
    guard let path = req["path"] as? String else { return }
    DispatchQueue.main.async {
      self.attachHeatmap(object, url: URL(fileURLWithPath: path))
      guard self.heatmaps[object] != nil else { return }
      func rgba(_ v: Any?) -> RGBA? {
        guard let a = v as? [Double], a.count == 3 else { return nil }
        return RGBA(a[0], a[1], a[2])
      }
      // Column / index first: a new column resets the range to its data…
      self.updateHeatmap(object) { s in
        if let c = req["column"] as? String { s.column = c }
        if let i = req["index_column"] as? String { s.indexColumn = i.isEmpty ? nil : i }
        if let c = rgba(req["low_color"]) { s.low = c }
        if let c = rgba(req["high_color"]) { s.high = c }
        if let w = req["through_white"] as? Bool { s.throughWhite = w }
      }
      // …then any explicit range.
      if req["min"] != nil || req["max"] != nil {
        self.updateHeatmap(object) { s in
          if let v = req["min"] as? Double { s.minValue = v }
          if let v = req["max"] as? Double { s.maxValue = v }
        }
      }
    }
  }

  // MARK: - Polling

  private func pollFeedback() {
    guard let inst = instance else { return }
    // While a movie export owns the core off the main thread, skip polling —
    // reading feedback here is a core access that would race the exclusive
    // exporter (which reshapes global state). The exporter drives its own
    // progress, so nothing is lost. (#58 L-59)
    if exportRenderActive { return }
    // The core drops a ray image on any change to the scene or view
    // (SceneInvalidateCopy); follow it, as desktop PyMOL does.
    if classicRayImage != nil && PyMOLBridge_HasSceneImage(inst) == 0 {
      classicRayImage = nil
    }
    guard let cStr = PyMOLBridge_GetFeedback(inst) else { return }
    let text = String(cString: cStr)
    PyMOLBridge_FreeFeedback(cStr)
    if !text.isEmpty {
      var printed: [String] = []  // console lines, appended in one update
      defer {
        if !printed.isEmpty {
          DispatchQueue.main.async {
            if self.startupQuiet {
              if printed.contains(where: { $0.contains("Error") || $0.contains("Traceback") }) {
                self.startupErrors = true
              }
              return
            }
            self.feedbackLog.append(contentsOf: printed)
          }
        }
      }
      // Check for ObjectPanel feedback before appending to log
      for line in text.components(separatedBy: "\n") {
        if line.hasPrefix("OBJPANEL:") {
          parseObjectPanelFeedback(line)
        } else if line.hasPrefix("OBJDETAIL:") {
          parseObjectDetailFeedback(line)
        } else if line.hasPrefix("SESSIONVP:") {
          let parts = line.dropFirst("SESSIONVP:".count).split(separator: ",")
          if parts.count == 2, let w = Double(parts[0]), let h = Double(parts[1]), h > 0 {
            let aspect = Float(w / h)
            DispatchQueue.main.async { self.setLetterboxAspect(aspect) }
          }
        } else if line.hasPrefix("SEQPANEL:") {
          parseSequencePanelFeedback(line)
        } else if line.hasPrefix("SEQSEL:") {
          parseSequenceSelectionFeedback(line)
        } else if line.hasPrefix("PLAYBACK:") {
          parsePlaybackFeedback(line)
        } else if line.hasPrefix("PLAYBACK_ERR:") {
          // swallow (don't flood the log with poll errors)
        } else if line.hasPrefix("MOVIEEXPORT:") {
          // cmd.movie_export (#581): a scripted export request.
          let json = String(line.dropFirst("MOVIEEXPORT:".count))
          DispatchQueue.main.async { self.startScriptedMovieExport(json) }
        } else if line.hasPrefix("SELPREVIEW:") {
          let v = String(line.dropFirst("SELPREVIEW:".count))
          let count = Int(v)
          DispatchQueue.main.async { self.selectionPreviewCount = count }
        } else if line.hasPrefix("MEASURE:") {
          parseMeasureFeedback(line)
        } else if line.hasPrefix("MEASURE_ERR:") {
          // swallow
        } else if line.hasPrefix("SETTINGS:ready") {
          loadSettingsCatalogFile()
        } else if line.hasPrefix("SETTINGS:err") {
          // swallow
        } else if line.hasPrefix("SETVAL:") {
          parseSetValFeedback(line)
        } else if line.hasPrefix("SEQCOLOR:") {
          let digest = String(line.dropFirst("SEQCOLOR:".count))
          if digest != lastSequenceColorDigest {
            let first = lastSequenceColorDigest == nil
            lastSequenceColorDigest = digest
            if !first && sequenceVisible { fetchSequences() }
          }
        } else if line.hasPrefix("HEATMAP:") || line.hasPrefix("CONTACTMAP:") {
          // Requests from the MCP tools color_by_csv / contact_map.
          handleMCPViewRequest(line)
        } else if line.hasPrefix("RENDEROPTS:") {
          parseRenderOptions(line)
        } else if line.hasPrefix("RAYIMG:") {
          parseClassicRayImage(line)
        } else if line.hasPrefix("MCP:") {
          let body = String(line.dropFirst("MCP:".count))
          if let colon = body.firstIndex(of: ":") {
            let kind = String(body[..<colon])
            let b64 = String(body[body.index(after: colon)...])
            let detail =
              Data(base64Encoded: b64)
              .flatMap { String(data: $0, encoding: .utf8) } ?? ""
            DispatchQueue.main.async {
              MCPServerManager.shared.handleFeedbackEvent(kind, detail)
            }
          }
        } else if !line.isEmpty {
          printed.append(line)
        }
      }
    }
  }

  // Run MCP tool work queued from the in-process server's HTTP request threads
  // on THIS (main) thread — the only thread where touching the embedded core is
  // safe (off-main cmd calls race SceneRenderMetal / corrupt the interpreter).
  // mmol.mcp.mainthread.run_on_main enqueues from the handler thread and blocks
  // until we drain here, once per ~100ms Timer tick. Gated on the server running
  // so there's no per-tick bridge call when MCP is off. A queued op runs
  // synchronously (it can block this tick — the same on-main invariant runHeavy
  // relies on); when idle the drain is an instant empty-queue check.
  private func drainMCPMainQueue() {
    guard isReady, MCPServerManager.shared.isRunning else { return }
    PyMOLBridge_RunPython("import mmol.mcp.mainthread as _mt; _mt.drain_main_thread_queue()")
  }

  private var objectPollCounter = 0

  private func pollObjects() {
    // Poll every 5th tick (~500ms at 100ms timer) to avoid flooding
    objectPollCounter += 1
    guard objectPollCounter % 5 == 0 else { return }
    // Not mid-scroll: a poll (tens of ms with many objects) would stall the
    // scrolling sequence strip or zoom; it runs on the next quiet tick.
    if CACurrentMediaTime() - lastScrollAt < 0.3 {
      objectPollCounter -= 1
      return
    }

    // Keep the sequence-panel selection highlight in sync with the active
    // selection (3D-view picks/selects reflect in the sequence).
    if sequenceVisible {
      fetchSequenceSelection(colorCheck: true)
    }

    // Run via runPython (raw PyRun), NOT runCommand/cmd.do — cmd.do echoes
    // the whole command block into the feedback log every poll, which floods
    // the log and starves the UI. The print('OBJPANEL:ready') still reaches
    // the feedback buffer (parsed by pollFeedback); only the echo is avoided.
    // poll_panel writes the list to a temp file and prints just that marker,
    // so the payload can't overflow the ~1KB feedback-line cap (#231).
    // Quiet: a poll must not cost a re-render of a static scene.
    objectPollEpoch = objectListEpoch
    // The undo snapshot of a settled change rides on the same quiet poll, so
    // it is taken after the screen has updated (mmol/undo.py).
    runPythonQuiet(
      "from mmol import inspector as _ai, undo as _undo\n_undo.checkpoint()\n_ai.poll_panel()")

    pollDetails()
    // Discover/refresh the timeline length + frame (cheap). Fast updates
    // during playback are handled by pollPlayback() on the 100ms tick.
    pollPlayback()
  }

  // Query active reps + per-rep settings + scene globals for the currently
  // EXPANDED object cards only (collapsed cards cost nothing). Emits an
  // OBJDETAIL feedback line via the bundled appkit_inspector module.
  // Public trigger: poll the expanded object's rep detail RIGHT NOW instead of
  // waiting up to ~500ms for the next pollObjects tick. Called when a card is
  // expanded so the rep list appears immediately (otherwise a heavy surface
  // build can starve the timer and the card lingers on "No representations").
  func refreshExpandedDetail() {
    pollDetails()
  }

  private func pollDetails() {
    // Always run (cheap when nothing expanded — just the ~8 scene gets) so
    // the Scene card stays fresh; per-object rep detail is queried only for
    // expanded cards.
    // At most one object card is open; the SCENE sentinel polls no object.
    var names: [String] = []
    if let d = expandedDetail, d != Self.sceneDetailKey { names = [d] }
    let pyList = names.map { "'\($0.replacingOccurrences(of: "'", with: ""))'" }.joined(
      separator: ", ")
    // The scene extent (Settings > Zoom) is measured only while that tab shows.
    let extent = activeTool == .settings ? "True" : "False"
    runPythonQuiet(
      "from mmol import inspector as _ai\n_ai.poll([\(pyList)], extent=\(extent))")
  }

  /// One object's `objmeta` entry -> ObjStateMeta.
  ///
  /// Static so a test can exercise it without an engine. The DEFAULTS are
  /// the point: they match each setting's own, so a payload from a build that
  /// predates a key reads as "not set" rather than as a value the user chose.
  static func parseObjMeta(_ m: [String: Any]) -> ObjStateMeta {
    ObjStateMeta(
      state: (m["state"] as? NSNumber)?.intValue ?? 1,
      overlayAll: ((m["all"] as? NSNumber)?.intValue ?? 0) != 0,
      titles: (m["titles"] as? [Any])?.map { $0 as? String ?? "" } ?? [])
  }

  // Parse the inspector JSON (written by appkit_inspector.poll to a temp file;
  // the feedback line is just the "OBJDETAIL:ready" trigger) → objectDetails +
  // sceneState. File-based to avoid the ~1KB feedback-line cap splitting the
  // payload and leaking continuation lines into the terminal log.
  func parseObjectDetailFeedback(_ line: String) {
    let path = TempChannel.path(TempChannel.Stem.objectDetail)
    guard let data = FileManager.default.contents(atPath: path),
      let root = try? JSONSerialization.jsonObject(with: data) as? [String: Any]
    else { return }

    var details: [String: [RepState]] = [:]
    if let detail = root["detail"] as? [String: Any] {
      for (obj, repsAny) in detail {
        guard let reps = repsAny as? [[String: Any]] else { continue }
        details[obj] = reps.map { r in
          var values: [String: Double] = [:]
          if let vals = r["vals"] as? [String: Any] {
            for (k, v) in vals { values[k] = (v as? NSNumber)?.doubleValue ?? 0 }
          }
          var settingColors: [String: String] = [:]
          if let cols = r["colors"] as? [String: Any] {
            for (k, v) in cols { settingColors[k] = v as? String ?? "inherit" }
          }
          var atomTransp: AtomTransp? = nil
          if let at = r["atom_transp"] as? [String: Any],
            let setting = at["setting"] as? String
          {
            atomTransp = AtomTransp(
              setting: setting,
              min: (at["min"] as? NSNumber)?.doubleValue ?? 0,
              max: (at["max"] as? NSNumber)?.doubleValue ?? 0)
          }
          return RepState(
            rep: r["rep"] as? String ?? "",
            visible: ((r["vis"] as? NSNumber)?.intValue ?? 1) != 0,
            values: values,
            color: r["color"] as? String ?? "inherit",
            effColor: r["eff"] as? String,
            settingColors: settingColors,
            atomTransp: atomTransp)
        }
      }
    }

    var scene = SceneState()
    if let sc = root["scene"] as? [String: Any] {
      for (k, v) in sc {
        if k == "outline_rgb", let arr = v as? [Any] {
          scene.outlineColor = arr.map { ($0 as? NSNumber)?.doubleValue ?? 0 }
        } else {
          scene.values[k] = (v as? NSNumber)?.doubleValue ?? 0
        }
      }
    }

    var meta: [String: ObjStateMeta] = [:]
    if let om = root["objmeta"] as? [String: Any] {
      for (obj, mAny) in om {
        guard let m = mAny as? [String: Any] else { continue }
        meta[obj] = PyMOLEngine.parseObjMeta(m)
      }
    }

    // Sessions saved by earlier builds carry hidden note bookmarks as
    // scenes; keep them out of the Scene menu and the timeline composer.
    let scenes = ((root["scenes"] as? [String]) ?? [])
      .filter { !$0.hasPrefix("__mmol_note_") }
    let curScene = (root["cur_scene"] as? String) ?? ""

    DispatchQueue.main.async {
      // Equality-guard: the ~500ms poll re-emits identical detail most of
      // the time; re-publishing unchanged values re-renders the inspector
      // and resets any open menu to the top. Only assign on real changes.
      if self.objectDetails != details { self.objectDetails = details }
      // Reps that appeared while the card was open (scripts, MCP…) get a tab.
      for (obj, reps) in details {
        let shown = Set(reps.map(\.rep))
        if let prev = self.lastShownReps[obj], !shown.subtracting(prev).isEmpty {
          self.cardReps[obj, default: []].formUnion(shown.subtracting(prev))
        }
        self.lastShownReps[obj] = shown
      }
      if self.sceneState != scene { self.sceneState = scene }
      if self.objectMeta != meta { self.objectMeta = meta }
      if self.sceneNames != scenes { self.sceneNames = scenes }
      if self.currentScene != curScene { self.currentScene = curScene }
    }
  }

  // Parse PLAYBACK:<json> → currentFrame / frameCount / isPlaying / fps / loop.
  // Mirrors core state; never overrides a user scrub in progress.
  func parsePlaybackFeedback(_ line: String) {
    let js = String(line.dropFirst("PLAYBACK:".count))
    guard let data = js.data(using: .utf8),
      let root = try? JSONSerialization.jsonObject(with: data) as? [String: Any]
    else { return }
    let frame = (root["frame"] as? NSNumber)?.intValue ?? 1
    let count = max((root["count"] as? NSNumber)?.intValue ?? 1, 1)
    let playing = ((root["playing"] as? NSNumber)?.intValue ?? 0) != 0
    let loop = ((root["loop"] as? NSNumber)?.intValue ?? 1) != 0
    let fps = (root["fps"] as? NSNumber)?.doubleValue ?? 15

    DispatchQueue.main.async {
      // Equality-guard every assignment: a @Published set fires the
      // publisher even when the value is unchanged, so unguarded writes
      // here would re-render the transport bar on every poll. Only the
      // genuinely-changing currentFrame ticks during playback.
      let pb = self.playback
      if pb.frameCount != count { pb.frameCount = count }
      if pb.isPlaying != playing { pb.isPlaying = playing }
      if pb.movieLoop != loop { pb.movieLoop = loop }
      if fps > 0 && pb.movieFPS != fps { pb.movieFPS = fps }
      // Don't fight an active drag; otherwise track the core frame.
      if !self.isScrubbing {
        let f = min(max(frame, 1), count)
        if pb.currentFrame != f { pb.currentFrame = f }
      }
      // The transport-bar visibility gate (observed by `engine`) flips
      // only when crossing the 1-frame threshold — not per frame.
      let has = count > 1
      if self.hasTimeline != has { self.hasTimeline = has }
    }
  }

}

// MARK: - Data models
// ObjectEntry is the canonical model, defined in ObjectPanel.swift.
// MoleculeObject is a typealias for backward compatibility.
typealias MoleculeObject = ObjectEntry
