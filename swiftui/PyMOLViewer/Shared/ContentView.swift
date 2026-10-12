// ContentView.swift — Main layout: viewport + side panels

import AppKit
import SwiftUI
import UniformTypeIdentifiers

// File-menu commands (defined on the App scene) post these; ContentView
// observes them and runs the matching open/save/export action, so the native
// menu items share the folder menu's logic + get standard shortcuts.
extension Notification.Name {
  static let mmolOpenFile = Notification.Name("mmol.menu.openFile")
  static let mmolFetch = Notification.Name("mmol.menu.fetch")
  static let mmolClearSession = Notification.Name("mmol.menu.clearSession")
  static let mmolSaveSession = Notification.Name("mmol.menu.saveSession")
  static let mmolSaveSessionAs = Notification.Name("mmol.menu.saveSessionAs")
  static let mmolExportImage = Notification.Name("mmol.menu.exportImage")
  static let mmolToggleTimeline = Notification.Name("mmol.menu.toggleTimeline")
  static let mcpOpenConnectSheet = Notification.Name("mmol.mcp.openConnectSheet")
}

struct ContentView: View {
  @EnvironmentObject var engine: PyMOLEngine

  // The top band's height the USER chose, as a fraction of the window height,
  // so restoring into a differently-sized window can't squeeze the viewport
  // away (#332). 0 means "never resized", which selects the absolute default
  // instead (#331). Keys live in PanelLayout.
  @AppStorage(PanelLayout.consoleFracKey) private var consoleFrac = 0.0
  // The top band's height the current drag started from, in points. nil when
  // no drag is in flight; the committed size lives in `consoleFrac`.
  @State private var macConsoleDragAnchor: CGFloat? = nil
  // The right column's width the USER chose, as a fraction of the window width
  // (#350). 0 means "never resized", which selects the absolute default instead.
  @AppStorage(PanelLayout.inspectorFracKey) private var inspectorFrac = 0.0
  // The right column's width the current seam drag started from, in points.
  @State private var macInspectorDragAnchor: CGFloat? = nil
  /// Sizes while a divider is dragged: layout follows these, and the stored
  /// fractions (@AppStorage, a defaults write and a full re-layout per change)
  /// are written once, on release.
  @State private var liveInspectorW: CGFloat? = nil
  @State private var liveConsoleH: CGFloat? = nil
  /// The Claude terminal's dragged height (0 until the user drags its seam),
  /// and its live height and drag anchor while they do.
  @AppStorage(PanelLayout.claudeHeightKey) private var claudeHeight = 0.0
  @State private var liveClaudeH: CGFloat? = nil
  @State private var claudeDragAnchor: CGFloat? = nil

  // "Fetch from PDB…" alert state (Open File… uses an NSOpenPanel directly).
  @State private var showMacFetch = false
  @State private var macFetchID = ""
  // File ▸ Clear Session confirmation.
  @State private var showClearConfirm = false
  // Scene ▸ Rename… text (the scene is engine.pendingSceneRename).
  @State private var sceneRenameText = ""
  // Drag-and-drop: true while a file is hovered over the viewport (draws a border).
  @State private var isViewportDropTargeted = false
  // Local key-down monitor token for the Esc → clear-selection handler
  // (issues #163 + #166). Installed in macOSLayout.onAppear, removed on
  // .onDisappear. NSEvent.addLocalMonitorForEvents (not .onKeyPress, which is
  // macOS 14+) keeps us on the macOS 13 deployment target.
  @State private var escKeyMonitor: Any?
  // Local key-down monitor token for cmd.set_key dispatch (#258). Same
  // rationale as escKeyMonitor above: MetalViewport declines first-responder
  // status (#73), so the viewport never receives keyDown and a monitor is the
  // only way to see keys at all. Installed/removed alongside it.
  @State private var pymolKeyMonitor: Any?
  // Local key-down monitor token for Tab / Shift+Tab → cycle the selection
  // mode (#319). Same rationale and lifecycle as the two monitors above.
  @State private var selectionModeKeyMonitor: Any?
  @EnvironmentObject private var mcpManager: MCPServerManager
  @State private var showConnectSheet = false

  var body: some View {
    macOSLayout
  }

  // MARK: - Layout

  // Drag handle under the top band (console + button grid) that resizes it,
  // persisting the share (#332). Drawn as a 1pt hairline with a taller
  // invisible hit area and a resize cursor.
  @ViewBuilder
  private func macConsoleDivider(windowHeight: CGFloat) -> some View {
    let maxH = macMaxConsoleHeight(windowHeight: windowHeight)
    Rectangle()
      .fill(hairlineColor)
      .frame(height: 1)
      .padding(.vertical, 2)
      .contentShape(Rectangle())
      .onHover { inside in
        if inside { NSCursor.resizeUpDown.push() } else { NSCursor.pop() }
      }
      .gesture(
        // Window coordinates: the divider moves with the drag, so its own
        // space would shift under the pointer and the size would jitter.
        DragGesture(minimumDistance: 1, coordinateSpace: .global)
          .onChanged { v in
            let start = macConsoleDragAnchor ?? topBandHeight(windowHeight: windowHeight)
            macConsoleDragAnchor = start
            engine.frozenDrawableGravity = .resizeAspect
            engine.suppressDrawableResize = true
            liveConsoleH = min(
              max(
                start + v.translation.height,
                PanelLayout.macMinConsoleHeight), maxH)
          }
          .onEnded { _ in
            if let h = liveConsoleH,
              let f = PanelLayout.consoleFrac(height: h, windowHeight: windowHeight)
            {
              consoleFrac = Double(f)
            }
            liveConsoleH = nil
            macConsoleDragAnchor = nil
            engine.suppressDrawableResize = false
          }
      )
  }

  // The sequence strip's height: two 13 pt rows per object plus the
  // horizontal scroller, up to 200 pt; 0 while hidden.
  private var classicSeqHeight: CGFloat {
    guard engine.sequenceVisible else { return 0 }
    return min(
      CGFloat(max(engine.sequences.count, 1)) * classicSeqRowsHeight
        + classicSeqScrollerRoom + 4, 200)
  }

  // The console's ceiling: what the window leaves above the divider, the
  // sequence strip and a minimum-height viewport.
  private func macMaxConsoleHeight(windowHeight: CGFloat) -> CGFloat {
    PanelLayout.maxConsoleHeight(
      windowHeight: windowHeight,
      viewportMin: PanelLayout.macViewportMinHeight + classicSeqHeight + 5)
  }

  // The console never drops below its minimum, even when the window is too
  // short for the usual ceiling.
  private func topBandHeight(windowHeight: CGFloat) -> CGFloat {
    if let h = liveConsoleH { return h }
    return max(
      PanelLayout.macMinConsoleHeight,
      PanelLayout.consoleHeight(
        frac: CGFloat(consoleFrac),
        windowHeight: windowHeight,
        defaultHeight: PanelLayout.macDefaultConsoleHeight,
        minHeight: PanelLayout.macMinConsoleHeight,
        maxHeight: macMaxConsoleHeight(windowHeight: windowHeight)))
  }

  // The right column's current width: the share the user dragged the seam to,
  // or the absolute default until they have (#350) — clamped every layout pass
  // so a window shrink can't squeeze the viewport below its minimum width.
  private func macInspectorWidth(windowWidth: CGFloat) -> CGFloat {
    if let w = liveInspectorW { return w }
    return PanelLayout.inspectorWidth(
      frac: CGFloat(inspectorFrac),
      windowWidth: windowWidth,
      maxWidth: PanelLayout.maxInspectorWidth(windowWidth: windowWidth))
  }

  // The seam between the viewport and the right column, doubling as the drag
  // handle that resizes it (#350). Same recipe as macConsoleDivider: a 1pt
  // hairline with a wider invisible hit strip, a resize cursor, and the
  // committed size persisted as a fraction of the window width.
  @ViewBuilder
  private func macInspectorDivider(windowWidth: CGFloat) -> some View {
    let maxW = PanelLayout.maxInspectorWidth(windowWidth: windowWidth)
    Rectangle()
      .fill(hairlineColor)
      .frame(width: 1)
      // A 9pt grab band grown entirely on the VIEWPORT side, so the visible
      // hairline still sits against the panel's leading edge. Padding is
      // layout width in the HStack, not an overlay, so it cannot steal
      // clicks from the panel's own controls.
      .padding(.leading, 8)
      .contentShape(Rectangle())
      .onHover { inside in
        if inside { NSCursor.resizeLeftRight.push() } else { NSCursor.pop() }
      }
      .gesture(
        DragGesture(minimumDistance: 1, coordinateSpace: .global)
          .onChanged { v in
            let start = macInspectorDragAnchor ?? macInspectorWidth(windowWidth: windowWidth)
            macInspectorDragAnchor = start
            // Freeze the Metal drawable while the seam moves so the
            // renderer doesn't reallocate its offscreen targets every
            // frame; one reshape snaps the viewport crisp on release.
            engine.frozenDrawableGravity = .center
            engine.suppressDrawableResize = true
            // The column sits on the RIGHT: dragging left grows it.
            liveInspectorW = min(
              max(
                start - v.translation.width,
                PanelLayout.macMinInspectorWidth), maxW)
          }
          .onEnded { _ in
            if let w = liveInspectorW,
              let f = PanelLayout.inspectorFrac(width: w, windowWidth: windowWidth)
            {
              inspectorFrac = Double(f)
            }
            liveInspectorW = nil
            macInspectorDragAnchor = nil
            engine.suppressDrawableResize = false
          }
      )
      // onEnded doesn't fire for a cancelled gesture; release the freeze
      // here too so the drawable can't stay stuck at its mid-drag size.
      .onDisappear {
        liveInspectorW = nil
        macInspectorDragAnchor = nil
        engine.suppressDrawableResize = false
      }
  }

  private var macOSLayoutBase: some View {
    // The window's content size, used for the top band's default and its
    // persisted share (#331/#332). A GeometryReader is the only reliable
    // source here: NSWindow isn't on screen yet at first layout.
    GeometryReader { winGeo in
      VStack(spacing: 0) {
        classicMacBody(windowSize: winGeo.size)
      }
      .alert("Fetch from PDB", isPresented: $showMacFetch) {
        TextField("PDB ID (e.g. 1ubq)", text: $macFetchID)
        Button("Fetch") { macFetch() }
        Button("Cancel", role: .cancel) {}
      } message: {
        Text("Download a structure from the RCSB PDB.")
      }
      .alert("Clear the session?", isPresented: $showClearConfirm) {
        Button("Clear", role: .destructive) { engine.clearSession() }
        Button("Cancel", role: .cancel) {}
      } message: {
        Text("Unsaved changes will be lost.")
      }
      .alert(
        "Rename Scene",
        isPresented: Binding(
          get: { engine.pendingSceneRename != nil },
          set: { if !$0 { engine.pendingSceneRename = nil } })
      ) {
        TextField("Scene name", text: $sceneRenameText)
        Button("Rename") {
          if let old = engine.pendingSceneRename { engine.renameScene(old, to: sceneRenameText) }
          engine.pendingSceneRename = nil
        }
        Button("Cancel", role: .cancel) { engine.pendingSceneRename = nil }
      }
      .onChange(of: engine.pendingSceneRename) { if let n = $0 { sceneRenameText = n } }
      .toolbar {
        macOpenToolbar
        // The "mMol" / Render Options / Settings / contact-map tabs, where the
        // title was. Without the glass capsule macOS 26 draws behind toolbar
        // items: it squeezed the "mMol" label and clipped the map tabs.
        if #available(macOS 26.0, *) {
          ToolbarItem(placement: .navigation) { classicTabStrip }
            .sharedBackgroundVisibility(.hidden)
        } else {
          ToolbarItem(placement: .navigation) { classicTabStrip }
        }
        // The MCP button, and left of it, while Claude drives the session,
        // its indicator. Both draw their own classic face, so they go without
        // the glass capsule, and sit in the toolbar so they never move the
        // window's content.
        if #available(macOS 26.0, *) {
          if mcpManager.activeTool {
            ToolbarItem(placement: .primaryAction) { MCPDrivingIndicator() }
              .sharedBackgroundVisibility(.hidden)
          }
          ToolbarItem(placement: .primaryAction) { MCPStatusView() }
            .sharedBackgroundVisibility(.hidden)
        } else {
          if mcpManager.activeTool {
            ToolbarItem(placement: .primaryAction) { MCPDrivingIndicator() }
          }
          ToolbarItem(placement: .primaryAction) { MCPStatusView() }
        }
      }
      // Native File-menu commands → reuse the same actions as the folder menu.
      .onReceive(NotificationCenter.default.publisher(for: .mmolOpenFile)) { _ in macOpenFile() }
      .onReceive(NotificationCenter.default.publisher(for: .mmolFetch)) { _ in
        macFetchID = ""
        showMacFetch = true
      }
      .onReceive(NotificationCenter.default.publisher(for: .mmolClearSession)) { _ in
        showClearConfirm = true
      }
      .onReceive(NotificationCenter.default.publisher(for: .mmolSaveSession)) { _ in saveSession() }
      .onReceive(NotificationCenter.default.publisher(for: .mmolSaveSessionAs)) { _ in
        saveSessionAs()
      }
      .onReceive(NotificationCenter.default.publisher(for: .mmolExportImage)) { _ in saveImage() }
      .onReceive(NotificationCenter.default.publisher(for: .mmolToggleTimeline)) { _ in
        engine.openTool(.movie)
      }
      .onReceive(NotificationCenter.default.publisher(for: .mcpOpenConnectSheet)) { _ in
        showConnectSheet = true
      }
      .sheet(isPresented: $engine.showSelectionBuilder) {
        SelectionBuilderSheet()
      }
      .sheet(isPresented: $showConnectSheet) {
        MCPConnectSheet().environmentObject(mcpManager)
      }
      .alert(
        "Allow Claude to control mMol?",
        isPresented: Binding(
          get: { mcpManager.pendingApproval },
          set: { if !$0 { mcpManager.pendingApproval = false } })
      ) {
        Button("Stop server", role: .destructive) { mcpManager.denyAndStop() }
        Button("Allow") { mcpManager.approveSession() }
      } message: {
        Text(
          "A local app connected to mMol and can now run commands, "
            + "run Python, and load structures until you stop it.")
      }
    }  // end GeometryReader
    .frame(maxWidth: .infinity, maxHeight: .infinity)
  }

  // MARK: - Title-bar tabs

  // "mMol" (the 3D view), then Render Options / Settings / Movie when open, then
  // one per open contact map, each with a close box. Picking a map tab swaps the
  // viewport area for the map; the others swap the right column (Movie also
  // docks the timeline under the viewport).
  private var classicTabStrip: some View {
    HStack(spacing: 2) {
      classicTab(
        title: "mMol", active: engine.activeContactMap == nil && engine.activeTool == nil,
        select: { engine.showTab() }, close: nil)
      ForEach(engine.toolTabs, id: \.self) { t in
        classicTab(
          title: t.rawValue, active: engine.activeTool == t,
          select: { engine.openTool(t) }, close: { engine.closeTool(t) })
      }
      ForEach(engine.contactMapTabs, id: \.self) { name in
        classicTab(
          title: name.hasPrefix("rama:")
            ? "\(name.dropFirst(5)) Ramachandran" : "\(name) contact map",
          active: engine.activeContactMap == name,
          select: { showContactMapTab(name) },
          close: { engine.closeContactMapTab(name) })
      }
    }
  }

  private func showContactMapTab(_ name: String) {
    engine.showTab(plot: name)
  }

  private func classicTab(
    title: String, active: Bool, select: @escaping () -> Void,
    close: (() -> Void)?
  ) -> some View {
    HStack(spacing: 5) {
      Text(title)
        .font(.system(size: 13, weight: active ? .semibold : .regular))
        .foregroundColor(active ? .primary : .secondary)
        .lineLimit(1)
      if let close {
        Button(action: close) {
          Image(systemName: "xmark").font(.system(size: 10, weight: .bold))
            .foregroundColor(.primary.opacity(0.8))
            .frame(width: 20, height: 20)
            .background(Color.primary.opacity(0.12))
            .clipShape(RoundedRectangle(cornerRadius: 4))
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .help("Close this tab")
      }
    }
    .padding(.horizontal, 10).padding(.vertical, 4)
    .background(
      RoundedRectangle(cornerRadius: 6)
        .fill(active ? Color.primary.opacity(0.14) : Color.clear)
    )
    .contentShape(Rectangle())
    .onTapGesture(perform: select)
  }

  // MARK: - Window body

  // Desktop PyMOL 2.5's window:
  //
  //   ┌ console log + input ─────────────┬ button grid ─┐
  //   ├ (sequence) ──────────────────────┼ object panel ┤
  //   │ Metal viewport                   │ mouse legend │
  //   │                                  │ |< < ■ ▶ … F │
  //   └──────────────────────────────────┴──────────────┘
  //
  // Two columns. The console's divider resizes the console against the
  // viewport only; the right column (button grid over the object panel) keeps
  // its own layout, so the object panel stays put. Render Options, Settings
  // and the other tools take over the object panel's place while their tab is
  // active.
  @ViewBuilder
  private func classicMacBody(windowSize: CGSize) -> some View {
    let rightW = macInspectorWidth(windowWidth: windowSize.width)
    HStack(spacing: 0) {
      VStack(spacing: 0) {
        CommandPanel()
          .frame(height: topBandHeight(windowHeight: windowSize.height))
        macConsoleDivider(windowHeight: windowSize.height)
        ZStack {
          VSplitView {
            if engine.sequenceVisible {
              // Sized to its rows (two 13 pt lines per object, plus the
              // horizontal scroller), like desktop PyMOL's strip.
              // Fixed, not just ideal: a two-axis ScrollView centres
              // content shorter than itself, which read as padding.
              let seqH = classicSeqHeight
              SequencePanel()
                .frame(minHeight: seqH, maxHeight: seqH)
                .id(engine.sequences.count)
            }
            macViewportStack
              // Drag a .pdb/.cif/.pse/etc. onto the viewport to load it
              // (same path as File ▸ Open / Finder "Open With").
              .onDrop(of: [.fileURL], isTargeted: $isViewportDropTargeted) { providers in
                handleViewportDrop(providers)
              }
              .overlay {
                if isViewportDropTargeted {
                  Rectangle()
                    .strokeBorder(ClassicUI.accent, lineWidth: 3)
                    .allowsHitTesting(false)
                }
              }
          }
          .frame(maxWidth: .infinity, maxHeight: .infinity)
          // Kept alive (not removed) while a map tab is up: rebuilding the
          // Metal view on every switch back would be slow.
          .opacity(engine.activeContactMap == nil ? 1 : 0)
          .allowsHitTesting(engine.activeContactMap == nil)
          // Every open map stays alive too, so zoom / cutoff survive a switch.
          ForEach(engine.contactMapTabs, id: \.self) { name in
            let shown = engine.activeContactMap == name
            Group {
              if name.hasPrefix("rama:") {
                RamachandranView(objectName: String(name.dropFirst(5)))
              } else {
                ContactMapView(objectName: name)
              }
            }
            .background(Color(nsColor: .windowBackgroundColor))
            .opacity(shown ? 1 : 0)
            .allowsHitTesting(shown)
          }
        }
      }
      // Never taller than the window, which would push the right column down.
      .frame(height: windowSize.height, alignment: .top)
      .clipped()

      macInspectorDivider(windowWidth: windowSize.width)
      VStack(spacing: 0) {
        // The console's default height, so the two seams line up until the
        // console is resized.
        ClassicButtonGrid()
          .frame(height: PanelLayout.macDefaultConsoleHeight)
        Rectangle().fill(hairlineColor).frame(height: 1).padding(.vertical, 2)
        switch engine.activeTool {
        case .render?: RenderOptionsPanel().frame(width: rightW)
        case .settings?: ViewerSettingsPanel().frame(width: rightW)
        case .movie?:
          ScrollView { MovieExportControls() }
            .frame(width: rightW)
            .background(ClassicUI.panelBg)
            .classicControls()
        case .align?: AlignPanel().frame(width: rightW)
        case .mutagenesis?: MutagenesisPanel().frame(width: rightW)
        case nil, .claude?:
          VStack(spacing: 0) {
            ObjectPanel()
            MousePanel()
            ClassicTransportBar()
          }
          .frame(width: rightW)
          .background(ClassicUI.panelBg)
        }
      }
      .frame(width: rightW)
    }
    .background(ClassicUI.panelBg)
  }

  private var macOSLayout: some View {
    macOSLayoutBase
      .preferredColorScheme(.dark)
      .tint(ClassicUI.accent)
      // Claude drives mMol through the MCP server: switching it off ends the
      // Claude session and closes its tab.
      .onChange(of: mcpManager.isRunning) { running in
        if !running, engine.toolTabs.contains(.claude) { engine.closeTool(.claude) }
      }
      .onChange(of: engine.isReady) { ready in
        if ready {
          engine.applyClassicDefaults()
          loadPymolrcAndAudit()
          engine.enableClassicRay()
        }
      }
      // A ▸ Contact Map…: a tab beside "mMol" in the title bar, shown at once.
      .onChange(of: engine.contactMapRequest) { name in
        guard let name else { return }
        if !engine.contactMapTabs.contains(name) { engine.contactMapTabs.append(name) }
        showContactMapTab(name)
        engine.contactMapRequest = nil
      }
      .onChange(of: engine.renderTabActive) { engine.showRenderSample($0) }
      .onAppear {
        initializeEngine()
        FrameLog.startHangLog()
        // Polls wait while the user scrolls (PyMOLEngine.pollObjects).
        NSEvent.addLocalMonitorForEvents(matching: .scrollWheel) { e in
          engine.lastScrollAt = CACurrentMediaTime()
          return e
        }
        installEscKeyMonitor()
        installPyMOLKeyMonitor()
        installSelectionModeKeyMonitor()
        // The sequence strip is on at every launch, listing whatever gets
        // loaded (Display ▸ Sequence still hides it).
        engine.sequenceVisible = true
      }
      .onDisappear {
        if let token = escKeyMonitor {
          NSEvent.removeMonitor(token)
          escKeyMonitor = nil
        }
        if let token = pymolKeyMonitor {
          NSEvent.removeMonitor(token)
          pymolKeyMonitor = nil
        }
        if let token = selectionModeKeyMonitor {
          NSEvent.removeMonitor(token)
          selectionModeKeyMonitor = nil
        }
      }
  }

  // The viewport column: the 3D Metal view + the Timeline dock. Extracted so
  // classicMacBody stays within the type-checker's size limit — same pattern
  // as macViewport itself.
  @ViewBuilder private var macViewportStack: some View {
    if engine.activeTool == .claude {
      // Measured so the Claude terminal fits inside the window: a fixed height
      // pushed its last rows, where Claude Code's input line is, off the
      // bottom. Only while it is open, so the usual layout is untouched.
      GeometryReader { geo in
        let claudeH =
          liveClaudeH
          ?? PanelLayout.claudePaneHeight(
            stackHeight: geo.size.height, preferred: CGFloat(claudeHeight))
        VStack(spacing: 0) {
          viewportWithRay(
            minHeight: min(
              PanelLayout.macViewportMinHeight, max(0, geo.size.height - claudeH - 5)))
          claudeDivider(stackHeight: geo.size.height, height: claudeH)
          ClaudeTerminalView().frame(height: claudeH)
        }
      }
    } else {
      VStack(spacing: 0) {
        viewportWithRay(minHeight: PanelLayout.macViewportMinHeight)
        // The timeline editor docks here while the Movie tab is active.
        if engine.movieTabActive {
          Divider()
          TimelinePanel()
        }
      }
    }
  }

  /// The seam above the Claude terminal, doubling as the drag handle that
  /// resizes it; the same recipe as macConsoleDivider. The dragged height is
  /// kept in points, so the terminal stays the size the user chose.
  private func claudeDivider(stackHeight: CGFloat, height: CGFloat) -> some View {
    Rectangle()
      .fill(hairlineColor)
      .frame(height: 1)
      .padding(.vertical, 2)
      .contentShape(Rectangle())
      .onHover { inside in
        if inside { NSCursor.resizeUpDown.push() } else { NSCursor.pop() }
      }
      .gesture(
        DragGesture(minimumDistance: 1, coordinateSpace: .global)
          .onChanged { v in
            let start = claudeDragAnchor ?? height
            claudeDragAnchor = start
            engine.frozenDrawableGravity = .resizeAspect
            engine.suppressDrawableResize = true
            // Dragging the seam up makes the terminal taller.
            liveClaudeH = PanelLayout.claudePaneHeight(
              stackHeight: stackHeight, preferred: start - v.translation.height)
          }
          .onEnded { _ in
            if let h = liveClaudeH { claudeHeight = Double(h) }
            liveClaudeH = nil
            claudeDragAnchor = nil
            engine.suppressDrawableResize = false
          }
      )
  }

  /// The viewport with the last `ray` result over it, until the core drops it.
  private func viewportWithRay(minHeight: CGFloat) -> some View {
    macViewport(minHeight: minHeight)
      .overlay {
        if let ray = engine.classicRayImage { ClassicRayOverlay(ray: ray) }
      }
  }

  // Shared by all three local key-down monitors below: a modal, sheet,
  // popover or NSPanel owns the keyboard while it is up, so keys belong to it
  // (pgup/pgdn must not change scenes behind an open Fetch-from-PDB sheet,
  // Esc must dismiss the sheet rather than clear the selection, …).
  //
  // Detect secondary windows STRUCTURALLY: SwiftUI's NSApp.mainWindow is
  // unreliable (often nil in Window-scene apps), so a `keyWindow != mainWindow`
  // test wrongly swallows keys on the main window.
  private func secondaryWindowOwnsKeys() -> Bool {
    if NSApp.modalWindow != nil { return true }
    // The Claude terminal takes every key (Esc, Tab, arrows…) while focused.
    if ClaudeSession.shared.hasFocus { return true }
    if let keyWindow = NSApp.keyWindow {
      // Sheets set isSheet; popovers/NSMenu helpers are NSPanels.
      return keyWindow.isSheet || keyWindow is NSPanel
    }
    return false
  }

  // The two keyboard-focus facts the dispatching monitors need:
  //   focused — an editable field is first responder (drives KeyRouting's
  //     Tier A non-US-keyboard yield for ALT/CTSH combos).
  //   editing — focused AND non-empty (drives Tier B, the arrows/home/end/
  //     ctrl-letter yield while the user is composing).
  // We require the text view to be editable or a field editor: the feedback
  // log is a selectable-but-not-editable NSTextView (ConsoleTextView, #406),
  // and it DOES take first responder when clicked — its `string` is the entire
  // log, so counting it as a focused field would silently disable arrows/
  // home/end for as long as the console held focus.
  private func textFocusFlags() -> (focused: Bool, editing: Bool) {
    let responder = NSApp.keyWindow?.firstResponder
    if let tv = responder as? NSTextView, tv.isEditable || tv.isFieldEditor {
      return (true, !tv.string.isEmpty)
    }
    if let tf = responder as? NSTextField, tf.isEditable {
      return (true, !tf.stringValue.isEmpty)
    }
    return (false, false)
  }

  // Esc → back out one level. A local key-down monitor rather than .onKeyPress
  // because the whole main window must catch Esc even when the viewport isn't
  // the SwiftUI focus — MetalViewport deliberately declines first-responder
  // status (issue #73) so the command line stays hot for typing.
  //
  // The ladder, in order (issues #163 + #166, then #235):
  //   (a) a sheet / panel / popover is up (their window is key, not the main
  //       mMol window) — Esc belongs to it, pass through so it dismisses;
  //   (b) an exclusive interaction mode (Move / Measure) is active —
  //       leave it, then consume;
  //   (c) otherwise — two-stage clear selection, then consume.
  // Non-Esc keys always pass through untouched.
  private func installEscKeyMonitor() {
    guard escKeyMonitor == nil else { return }
    escKeyMonitor = NSEvent.addLocalMonitorForEvents(matching: .keyDown) { event in
      guard event.keyCode == 53 else { return event }  // 53 = Esc
      // (a) A modal/sheet/panel/popover owns the interaction → let it handle
      // Esc (dismiss). See secondaryWindowOwnsKeys above.
      if secondaryWindowOwnsKeys() { return event }
      // NOTE: intentionally does NOT defer to a focused text field — Esc
      // acts regardless of keyboard focus (incl. while the command-line box
      // is focused), per product decision. Only true modal/sheet/panel
      // windows above still get Esc for dismissal.
      //
      // (b) Leaving an interaction mode outranks clearing the selection —
      // #166 anticipated exactly this ("exiting a mode should take priority
      // before the selection stages kick in"). Routed through the engine so
      // Esc reuses each mode's own exit path rather than duplicating it.
      if engine.exitActiveInteractionMode() { return nil }
      // (c) No mode was active → the selection stages.
      engine.escapeClearSelection()
      return nil  // consume — don't beep or propagate
    }
  }

  // cmd.set_key dispatch (#258). A local key-down monitor for the same reason
  // the Esc handler is one: MetalViewport deliberately declines
  // first-responder status (#73), so keyDown never reaches the viewport and
  // the app would otherwise never see a key at all.
  //
  // The consume rule is the whole conflict policy: KeyRouting decides whether
  // an event is even a candidate, then we consume it ONLY if a binding
  // actually fired. So an unbound ⌃M falls through to the Move menu item
  // naturally, while a user-bound ⌃M shadows it — no reserved-key table
  // needed. Esc (keyCode 53) never yields a token, so this monitor and the Esc
  // one never contend.
  private func installPyMOLKeyMonitor() {
    guard pymolKeyMonitor == nil else { return }
    pymolKeyMonitor = NSEvent.addLocalMonitorForEvents(matching: .keyDown) { event in
      // Shared modal/sheet/panel guard (secondaryWindowOwnsKeys above):
      // if a sheet, alert, or popover is the key window, keys belong to
      // it — pgup/pgdn must not change scenes behind an open
      // Fetch-from-PDB sheet, for example.
      if secondaryWindowOwnsKeys() { return event }

      let focus = textFocusFlags()
      guard
        let token = KeyRouting.token(
          keyCode: event.keyCode,
          charactersIgnoringModifiers: event.charactersIgnoringModifiers,
          modifiers: event.modifierFlags,
          textFieldFocused: focus.focused,
          textEditingActive: focus.editing)
      else { return event }
      return engine.invokeKeyBinding(token) ? nil : event
    }
  }

  // Tab / Shift+Tab → cycle the selection mode (#319), i.e. the same
  // atoms/residues/chains/… list SelectionModeMenu shows. A third local
  // monitor for the same reason as the two above (#73), and it cannot contend
  // with them: Esc only looks at keyCode 53, and Tab carries neither CTRL nor
  // ALT so KeyRouting.token always returns nil for it.
  //
  // Which Tab presses are ours lives in KeyRouting.selectionModeStep (pure,
  // unit-tested) — bare Tab and Shift+Tab only, and never while a field is
  // mid-edit, so the command line's Tab completion is untouched. The
  // modal/sheet/panel guard stays here alongside the other monitors', and as
  // with the Esc monitor the event is consumed only when the shortcut
  // actually fires; otherwise it falls through to normal focus navigation.
  private func installSelectionModeKeyMonitor() {
    guard selectionModeKeyMonitor == nil else { return }
    selectionModeKeyMonitor = NSEvent.addLocalMonitorForEvents(matching: .keyDown) { event in
      guard
        let step = KeyRouting.selectionModeStep(
          keyCode: event.keyCode,
          modifiers: event.modifierFlags,
          textEditingActive: textFocusFlags().editing)
      else { return event }
      if secondaryWindowOwnsKeys() { return event }
      SelectionModeMenu.cycle(engine, forward: step > 0)
      return nil  // consume — don't move keyboard focus as well
    }
  }

  // The macOS viewport: the Metal view plus its floating overlays and the
  // right-click context menu. Extracted from macOSLayout's body so the
  // type-checker can resolve each in isolation (the inline chain tripped the
  // "unable to type-check in reasonable time" limit).
  @ViewBuilder
  private func macViewport(minHeight: CGFloat) -> some View {
    MetalViewport()
      .frame(minWidth: 400, minHeight: minHeight)
      .layoutPriority(1)
      .overlay(alignment: .top) {
        if engine.measureMode != nil {
          measureOverlay
        } else if engine.interactionMode == .move {
          moveOverlay
        }
      }
      // Box Select rubber band (#358). Unlike the Move gizmo (a 3D CGO in
      // the Metal scene) the box is pure 2D screen furniture, so SwiftUI
      // draws it; MetalViewport hit-tests the same rectangle in NDC.
      .overlay {
        if engine.interactionMode == .boxSelect {
          BoxSelectRectOverlay(
            rect: engine.boxRect,
            accent: ClassicUI.accent,
            onClose: { engine.endBoxSelection() })
        }
      }
      // (The Move-mode gizmo is a 3D CGO object rendered in the Metal
      // scene by mmol/move.py; no SwiftUI overlay is needed. Input is
      // hit-tested against the projected geometry in MetalViewport.)
      // Pick-debug crosshair: marks exactly where the last click landed,
      // so a screenshot shows click-vs-selection offset.
      .overlay { debugClickMarker }
      // Right-click (or a trackpad double-click) in the viewport picks the
      // atom/residue under the pointer (or empty space) and sets
      // engine.longPressHit; open the classic menu for it at the pointer.
      // PyMOL's own pop-up menu is never drawn under this Metal backend
      // (internal_gui=0).
      .onChange(of: engine.longPressHit) { hit in
        guard let hit else { return }
        engine.longPressHit = nil
        ClassicMenuController.shared.showAtPointer(viewportMenuItems(hit))
      }
  }

  // Allowed import types — the molecular/map/session extension set.
  private var macImportTypes: [UTType] {
    let exts = [
      "pdb", "ent", "cif", "mmcif", "mcif", "sdf", "mol", "mol2",
      "xyz", "pdbqt", "pqr", "mae", "pse", "ccp4", "mrc", "map",
      "dx", "mtz", "fasta", "pir",
    ]
    return exts.compactMap { UTType(filenameExtension: $0) } + [.data]
  }

  // Drag-and-drop: load each dropped file URL through the same path as the Open
  // menu / Finder "Open With" — loadOpenedFile handles security scope, temp copy,
  // name sanitizing, and engine-not-ready retry. A dropped .pse restores a session.
  private func handleViewportDrop(_ providers: [NSItemProvider]) -> Bool {
    let engine = self.engine
    var accepted = false
    for provider in providers where provider.canLoadObject(ofClass: URL.self) {
      accepted = true
      _ = provider.loadObject(ofClass: URL.self) { url, _ in
        guard let url, url.isFileURL else { return }
        Task { @MainActor in loadOpenedFile(url, into: engine) }
      }
    }
    return accepted
  }

  // Open a molecule/session via NSOpenPanel and load it. PyMOL infers the format
  // from the extension; the object name is the filename stem (sanitized).
  private func macOpenFile() {
    let panel = NSOpenPanel()
    panel.allowsMultipleSelection = false
    panel.canChooseDirectories = false
    panel.canChooseFiles = true
    panel.allowedContentTypes = macImportTypes
    panel.title = "Open Structure"
    guard panel.runModal() == .OK, let url = panel.url else { return }
    // A chosen .pse replaces the whole session (no multi-window yet, #29) —
    // same warn/save/cancel stopgap as the Finder open path (issue #349).
    guard confirmReplaceSessionIfNeeded(opening: url, engine: engine) else { return }
    let raw = url.deletingPathExtension().lastPathComponent
    var name = String(raw.map { $0.isLetter || $0.isNumber ? $0 : "_" })
    if name.isEmpty { name = "mol" }
    engine.loadStructure(path: url.path, name: name)
    // Track an opened .pse as the current document so ⌘S overwrites it; a
    // non-.pse structure clears the tracked document.
    engine.currentSessionURL = (url.pathExtension.lowercased() == "pse") ? url : nil
  }

  private func macFetch() {
    let id = macFetchID.trimmingCharacters(in: .whitespaces)
      .replacingOccurrences(of: "'", with: "")
    guard !id.isEmpty else { return }
    engine.fetchStructure(id: id)
  }

  // The viewport's right-click menu. Empty space → scene-level actions; a
  // hit → residue-scoped actions on hit.sel (an obj/chain/resi selector).
  private func viewportMenuItems(_ hit: LongPressHit) -> [ClassicMenuItem] {
    var items: [ClassicMenuItem] = [.disabled(hit.title), .separator]
    if hit.isEmpty {
      items += [
        .action("Reset view") { engine.runCommand("reset") },
        .action("Deselect all") { engine.runCommand("deselect") },
      ]
      return items
    }
    let sel = hit.sel
    let colors: [(String, Color)] = [
      ("red", .red), ("orange", .orange), ("yellow", .yellow), ("green", .green),
      ("cyan", .cyan), ("blue", .blue), ("magenta", Color(red: 1, green: 0, blue: 1)),
      ("white", .white),
    ]
    items += [
      .action("Zoom to residue") { engine.runCommand("zoom (\(sel)), animate=1") },
      .action("Select residue") {
        engine.runCommand("select sele, (?sele) or (\(sel))\nenable sele")
      },
      .action("Label residue") {
        engine.runCommand("label first (\(sel)), '\(hit.resn)\(hit.resi)'")
      },
      .action("Hide residue") { engine.runCommand("hide everything, (\(sel))") },
      .action("Center here") { engine.runCommand("center (\(sel))") },
      .submenu(
        "Color",
        colors.map { name, swatch in
          .action(name.capitalized, swatch: swatch) { engine.runCommand("color \(name), (\(sel))") }
        } + [
          .separator,
          .action("By element") {
            engine.runCommand("python\nfrom pymol import util; util.cnc('(\(sel))')\npython end")
          },
        ]),
    ]
    return items
  }

  // Thin seam between the viewport and the docked panels.
  private var hairlineColor: Color { ClassicUI.chromeText.color.opacity(0.18) }

  // MARK: - Toolbar

  // The folder menu: open, fetch, save and share. ⌘O / ⇧⌘O live in the
  // File menu.
  private var macOpenToolbar: some ToolbarContent {
    ToolbarItem(placement: .navigation) {
      ClassicMenuTrigger {
        HStack(spacing: 3) {
          Image(systemName: "folder").font(.system(size: 14))
          Image(systemName: "chevron.down").font(.system(size: 8, weight: .semibold))
        }
        .foregroundColor(ClassicUI.chromeText.color)
        .padding(.horizontal, 8).frame(height: 24)
      } items: {
        let formats: [(String, String)] = [
          ("PDB (.pdb)", "pdb"), ("mmCIF (.cif)", "cif"), ("SDF (.sdf)", "sdf"),
          ("MOL (.mol)", "mol"), ("MOL2 (.mol2)", "mol2"), ("XYZ (.xyz)", "xyz"),
          ("PQR (.pqr)", "pqr"),
        ]
        return [
          .action("Open File…") { macOpenFile() },
          .action("Fetch from PDB…") {
            macFetchID = ""
            showMacFetch = true
          },
          .separator,
          .submenu(
            "Save Structure",
            formats.map { title, ext in .action(title) { saveStructure(ext: ext) } } + [
              .separator,
              .action("VRML (.wrl)") { saveStructure(ext: "wrl") },
              .action("POV-Ray (.pov)") { saveStructure(ext: "pov") },
            ]),
          .action("Save Session (.pse)…") { saveSession() },
          .submenu("Share", [.action("Session…") { shareSession() }]),
        ]
      }
      .help("Open, fetch, save and share")
    }
  }

  // MARK: - Save / export

  // File ▸ Export Image… saves the image the last `ray` left on screen, at its
  // own size, as desktop PyMOL's does. With no ray image up there is nothing
  // to save: say how to make one rather than rendering something else.
  private func saveImage() {
    guard let ray = engine.classicRayImage else {
      let alert = NSAlert()
      alert.messageText = "No ray-traced image to export"
      alert.informativeText =
        "Run ray first (for example, type ray in the console), "
        + "then choose Export Image… to save the result."
      alert.runModal()
      return
    }
    let panel = NSSavePanel()
    panel.allowedContentTypes = [.png]
    panel.nameFieldStringValue = "render.png"
    panel.canCreateDirectories = true
    panel.title = "Save Ray-Traced Image (\(ray.pixelWidth) × \(ray.pixelHeight))"
    guard panel.runModal() == .OK, let url = panel.url else { return }
    engine.saveHeldRayImage(to: url.path)
  }

  // ⌘S: overwrite the currently-open .pse with no panel. Falls back to Save As
  // when no document is tracked (never-saved session, or a non-.pse was opened).
  private func saveSession() {
    if let url = engine.currentSessionURL {
      engine.saveSession(to: url)
    } else {
      saveSessionAs()
    }
  }

  // ⇧⌘S: always show the Save panel, prefilled from the tracked document when
  // there is one, then save to the chosen URL and make it the open document.
  private func saveSessionAs() {
    let panel = NSSavePanel()
    if let pse = UTType(filenameExtension: "pse") { panel.allowedContentTypes = [pse] }
    if let current = engine.currentSessionURL {
      panel.directoryURL = current.deletingLastPathComponent()
      panel.nameFieldStringValue = current.lastPathComponent
    } else {
      panel.nameFieldStringValue = "session.pse"
    }
    panel.canCreateDirectories = true
    panel.title = "Save Session"
    guard panel.runModal() == .OK, let url = panel.url else { return }
    engine.saveSession(to: url)
  }

  // Save the whole scene to a molecular or 3D file. cmd.save infers the format
  // from the extension; the user types the extension (.pdb/.cif/.mol2/.sdf/.xyz
  // /.mae/.pqr molecular, or .wrl/.pov 3D — glTF/COLLADA/STL aren't available
  // on this libxml-off / NO_OPENGL build).
  private func saveStructure(ext: String) {
    let panel = NSSavePanel()
    if let t = UTType(filenameExtension: ext) { panel.allowedContentTypes = [t] }
    panel.allowsOtherFileTypes = true
    panel.nameFieldStringValue = "structure.\(ext)"
    panel.canCreateDirectories = true
    panel.title = "Save Structure (.\(ext))"
    guard panel.runModal() == .OK, let url = panel.url else { return }
    engine.runPython("from pymol import cmd as _c\n_c.save(r'''\(url.path)''')")
  }

  private func shareSession() {
    let tmp = (NSTemporaryDirectory() as NSString).appendingPathComponent("pymol_share.pse")
    engine.runPython("from pymol import cmd as _c; _c.save(r'''\(tmp)''')")
    let sessionURL = URL(fileURLWithPath: tmp)
    presentShare(items: [sessionURL])
  }

  private func presentShare(items: [URL]) {
    guard !items.isEmpty else { return }
    guard let window = NSApp.keyWindow, let anchor = window.contentView else { return }
    let picker = NSSharingServicePicker(items: items)
    picker.show(relativeTo: .zero, of: anchor, preferredEdge: .minY)
  }

  // MARK: - Measurement overlay (shared)

  // A thin bar over the top of the viewport while measure mode is active:
  // pick the measurement type, see the live prompt/result, clear, or exit.
  // Pick-debug: a cyan crosshair + ring at the exact pixel of the last click,
  // overlaid on the viewport (same top-down coordinate space as the MTKView).
  // Lets a screenshot directly compare where the user clicked vs where the pink
  // selection square rendered. Only present when PYMOL_PICKDEBUG is set.
  @ViewBuilder
  private var debugClickMarker: some View {
    if PyMOLEngine.debugPickEnabled, let p = engine.debugClickPoint {
      ZStack {
        Circle().stroke(Color.cyan, lineWidth: 1.5).frame(width: 22, height: 22)
        Rectangle().fill(Color.cyan).frame(width: 1.5, height: 14)
        Rectangle().fill(Color.cyan).frame(width: 14, height: 1.5)
      }
      .position(p)
      .allowsHitTesting(false)
    }
  }

  private var measureOverlay: some View {
    HStack(spacing: 10) {
      Picker(
        "",
        selection: Binding(
          get: { engine.measureMode ?? .distance },
          set: { engine.setMeasureMode($0) })
      ) {
        Text("Distance").tag(MeasureKind.distance)
        Text("Angle").tag(MeasureKind.angle)
        Text("Dihedral").tag(MeasureKind.dihedral)
      }
      .pickerStyle(.segmented)
      .frame(maxWidth: 240)
      Text(engine.measureStatus)
        .font(.system(size: 12, weight: .medium))
        .foregroundColor(ClassicUI.chromeText.color)
        .lineLimit(1).minimumScaleFactor(0.7)
      Spacer(minLength: 0)
      Button {
        engine.clearMeasurements()
      } label: {
        Image(systemName: "trash").foregroundColor(ClassicUI.chromeText.color)
      }.buttonStyle(.plain).help("Delete all measurements")
      Button {
        engine.setMeasureMode(nil)
      } label: {
        Image(systemName: "xmark.circle.fill")
          .foregroundColor(ClassicUI.chromeText.color.opacity(0.6))
      }.buttonStyle(.plain).accessibilityLabel("Exit measure mode")
    }
    .padding(.horizontal, 12).padding(.vertical, 8)
    .background(ClassicUI.panelBg)
    .tint(ClassicUI.accent)
  }

  // Move-mode overlay bar (mirrors measureOverlay): Move/Rotate tool toggle,
  // active-object dropdown, live readout, reset, exit.
  private var moveOverlay: some View {
    HStack(spacing: 10) {
      Image(systemName: "move.3d")
        .foregroundColor(ClassicUI.accent)
      ClassicMenuTrigger {
        HStack(spacing: 3) {
          Image(systemName: "scope")
          Text(engine.activeMoveObject ?? "Tap an object")
            .lineLimit(1)
          Image(systemName: "chevron.down").font(.system(size: 9))
        }
        .font(.system(size: 12, weight: .medium))
        .foregroundColor(ClassicUI.chromeText.color)
      } items: {
        let names = engine.objects.filter { !$0.isSelection }.map { $0.name }
        if names.isEmpty { return [.disabled("No objects loaded")] }
        return names.map { n in
          .check(n, engine.activeMoveObject == n) { engine.setActiveMoveObject(n) }
        }
      }

      // Adjust-frame controls: the toggle (macOS: also hold Shift) makes the
      // gizmo grey out and its controls re-anchor the frame (origin + tilt)
      // instead of moving the structure; reset snaps the frame back to the
      // automatic molecular center.
      Button {
        engine.adjustFrameToggle.toggle()
      } label: {
        Image(systemName: "gyroscope")
          .foregroundColor(
            engine.adjustFrameActive
              ? ClassicUI.accent
              : ClassicUI.chromeText.color)
      }
      .buttonStyle(.plain)
      .disabled(engine.activeMoveObject == nil)
      .help(
        "Adjust the gizmo frame — drag to set its origin & tilt (or hold Shift). Moves the gizmo, not the structure."
      )

      Button {
        engine.resetGizmoFrame()
      } label: {
        Image(systemName: "arrow.counterclockwise")
          .foregroundColor(ClassicUI.chromeText.color)
      }
      .buttonStyle(.plain)
      .disabled(engine.activeMoveObject == nil)
      .help("Reset the gizmo frame to the automatic center")

      Spacer(minLength: 0)

      Button {
        engine.resetActiveMovePosition()
      } label: {
        Image(systemName: "arrow.uturn.backward")
          .foregroundColor(ClassicUI.chromeText.color)
      }
      .buttonStyle(.plain)
      .disabled(engine.activeMoveObject == nil)
      .help("Reset this object's position")

      Button {
        engine.setInteractionMode(.viewing)
      } label: {
        Image(systemName: "xmark.circle.fill")
          .foregroundColor(ClassicUI.chromeText.color.opacity(0.6))
      }.buttonStyle(.plain).accessibilityLabel("Exit move mode")
    }
    .padding(.horizontal, 12).padding(.vertical, 8)
    .background(ClassicUI.panelBg)
    .tint(ClassicUI.accent)
  }

  // MARK: - Initialization

  private func initializeEngine() {
    guard !engine.isReady else { return }
    let resourcePath = Bundle.main.resourcePath ?? ""
    engine.initialize(resourcePath: resourcePath)
  }

  // Run ~/.pymolrc(.py) after applyClassicDefaults, so a user's startup
  // script can override them (e.g. a custom bg_color) — matching desktop PyMOL,
  // where .pymolrc runs after all built-in defaults are set. Then the
  // render-options snapshot, and the shortcut audit: ~/.pymolrc may bind a key
  // mMol also uses as a menu shortcut, which it says once (#258).
  private func loadPymolrcAndAudit() {
    // Load first in its own call so a stale mmol_keys import can never
    // abort the rc load (#258 review note).
    engine.runPython(
      "from mmol import pymolrc as _mmolrc; "
        + "_mmolrc.load(); "
        + "from mmol import render_options as _ro; _ro.snapshot()")
    // Undo and redo record changes from here on, not the startup script's.
    engine.runPython("from mmol import undo as _undo; _undo.install()")
    engine.runPython(
      "from mmol import keys as _mmol_keys; "
        + "_mmol_keys.audit_shadowed()")
    // The feedback poll is still delivering startup lines; then the banner.
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.8) { engine.finishStartupLog() }
    // Test affordance: PYMOL_AUTOCMD_LATE="cmd1; cmd2" runs once the defaults
    // and startup script are in, as if typed at the prompt (PYMOL_AUTOCMD
    // runs at engine init, before either).
    if let late = ProcessInfo.processInfo.environment["PYMOL_AUTOCMD_LATE"] {
      let delay =
        Double(ProcessInfo.processInfo.environment["PYMOL_AUTOCMD_LATE_DELAY"] ?? "") ?? 2.0
      DispatchQueue.main.asyncAfter(deadline: .now() + delay) {
        for one in late.split(separator: ";") {
          engine.runCommand(one.trimmingCharacters(in: .whitespaces))
        }
        // …PYMOL_AUTOTOOL=<tab title, e.g. Align> opens a tool tab…
        if let t = ProcessInfo.processInfo.environment["PYMOL_AUTOTOOL"].flatMap(ToolTab.init) {
          if t == .claude { engine.openClaude() } else { engine.openTool(t) }
        }
        // …and PYMOL_AUTOCONTACTMAP=<object> (rama:<object>) opens its plot.
        if let obj = ProcessInfo.processInfo.environment["PYMOL_AUTOCONTACTMAP"] {
          engine.contactMapRequest = obj
        }
      }
    }
  }
}

/// The Box Select rubber band (#358): the rectangle, its eight grab handles, a
/// dim wash over everything outside it, and the ✕ that leaves the tool.
///
/// Everything except the ✕ is chrome — `allowsHitTesting(false)` — because the
/// drag routing lives in MetalViewport, which hit-tests the SAME rectangle in
/// NDC (BoxRect.edges). Two hit-testers over one rectangle would be two chances
/// to disagree; here the SwiftUI layer only draws. The ✕ is the one exception,
/// and it sits OUTSIDE the frame so it can't swallow the corner resize handle.
struct BoxSelectRectOverlay: View {
  let rect: BoxRect?
  let accent: Color
  var onClose: () -> Void = {}

  /// How far outside the top-right corner the ✕ sits, in points.
  private let closeInset: CGFloat = 13

  var body: some View {
    GeometryReader { geo in
      if let rect = rect {
        let r = rect.inPoints(geo.size)
        ZStack {
          Group {
            // Even-odd fill of (whole viewport + box) = everything but the box.
            Path { p in
              p.addRect(CGRect(origin: .zero, size: geo.size))
              p.addRect(r)
            }
            .fill(Color.black.opacity(0.22), style: FillStyle(eoFill: true))

            Path { $0.addRect(r) }
              .stroke(accent, style: StrokeStyle(lineWidth: 1.5, dash: [6, 4]))

            ForEach(Array(Self.handleCenters(r).enumerated()), id: \.offset) { _, c in
              Circle()
                .fill(accent)
                .frame(width: 8, height: 8)
                .position(c)
            }
          }
          .allowsHitTesting(false)

          Button(action: onClose) {
            Image(systemName: "xmark.circle.fill")
              .font(.system(size: 18))
              .symbolRenderingMode(.palette)
              .foregroundStyle(Color.white, accent)
          }
          .buttonStyle(.plain)
          .position(x: r.maxX + closeInset, y: r.minY - closeInset)
          .help("Close Box Select — the selection stays (Esc)")
          .accessibilityLabel("Close box select")
        }
      }
    }
  }

  /// The eight resize handles: every corner and edge midpoint (the centre is
  /// the translate area, which needs no dot).
  static func handleCenters(_ r: CGRect) -> [CGPoint] {
    [r.minX, r.midX, r.maxX].flatMap { x in
      [r.minY, r.midY, r.maxY].map { CGPoint(x: x, y: $0) }
    }
    .filter { !($0.x == r.midX && $0.y == r.midY) }
  }
}
