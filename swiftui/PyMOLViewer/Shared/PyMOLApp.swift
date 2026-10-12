// PyMOLApp.swift — SwiftUI entry point for the macOS app

import AppKit
import SwiftUI
import UniformTypeIdentifiers

// macOS delivers files to open through NSApplicationDelegate.application(_:open:)
// as the COMPLETE [URL] array — Terminal `open a.pdb b.pdb c.pdb`, Finder "Open
// With" multi-select, and drag-drop onto the Dock icon all arrive this way.
// SwiftUI's .onOpenURL surfaces only the FIRST of those URLs and silently drops
// the rest, so a multi-file open loaded just one file (issue #222). Implementing
// the delegate method and looping over the whole array is the reliable way to
// receive every file; macOS therefore routes OS opens here instead of .onOpenURL.
final class MMolAppDelegate: NSObject, NSApplicationDelegate {
  func application(_ application: NSApplication, open urls: [URL]) {
    // Called on the main thread; hop to the main actor to reach loadOpenedFile
    // (@MainActor), mirroring ContentView's drag-drop handler.
    Task { @MainActor in handleOpenedURLs(urls, into: PyMOLEngine.shared) }
  }

  func applicationWillTerminate(_ notification: Notification) {
    // The tempfile channels are named after this pid (#399), so nothing will
    // ever reuse them — without this they'd accumulate one set per run.
    TempChannel.removeAll()
  }
}

// MARK: - Keyboard shortcut table (#360, #361)

/// The one table behind the app's discoverable shortcuts: each is registered on
/// a real menu-bar command (a menu key equivalent is what fires reliably on
/// macOS), from these constants.
///
/// Every ⌃-letter entry must be mirrored in modules/mmol/keys.py
/// (APP_SHORTCUTS), the launch-time audit that warns when ~/.pymolrc binds
/// over a menu shortcut. ⌘ entries need no mirror: KeyRouting passes ⌘ events
/// straight through to the menus, so a user binding can never shadow them.
enum AppShortcuts {
  // Tool picker (#360): the exclusive viewport modes, one ⌃-mnemonic each —
  // Move / mEasure — plus Box Select (⌃S, S for "select") on the same scheme.
  // ⌃E is the first letter of "measure" not already spoken for (⌃M is Move).
  static let moveTool = KeyboardShortcut("m", modifiers: .control)
  static let measureTool = KeyboardShortcut("e", modifiers: .control)
  static let boxSelect = KeyboardShortcut("s", modifiers: .control)
  /// Display ▸ Sequence.
  static let sequencePane = KeyboardShortcut("2", modifiers: .command)

  /// Every entry above. The collision test runs off this list, so a new
  /// shortcut must be added here too.
  static let all: [KeyboardShortcut] = [moveTool, measureTool, boxSelect, sequencePane]
}

struct PyMOLApp: App {
  @StateObject private var engine = PyMOLEngine.shared
  @StateObject private var mcp = MCPServerManager.shared
  // Receives OS file-open events (application(_:open:)); see MMolAppDelegate.
  @NSApplicationDelegateAdaptor(MMolAppDelegate.self) private var appDelegate

  var body: some Scene {
    // Single, unique window (`Window`, not `WindowGroup`): mMol's engine is
    // one shared PyMOL session, so a second window would only duplicate the
    // same view and stay in sync. `Window` drops the "New Window" command and
    // re-focuses the existing window instead. Per-window sessions: issue #29.
    Window("mMol", id: "mmol-main") { rootView }
      .windowStyle(.titleBar)
      .windowToolbarStyle(.unifiedCompact)  // shorter bar, smaller corner radius
      .defaultSize(width: 1200, height: 800)
      .commands { macCommands }
  }

  // Content of the main window: engine injection.
  @ViewBuilder private var rootView: some View {
    ContentView()
      .environmentObject(engine)
      .environmentObject(engine.playback)
      // Bring the app/window to the front on launch (a GUI app should
      // foreground itself; also lets it be launched from a terminal).
      .onAppear { NSApplication.shared.activate(ignoringOtherApps: true) }
      // No title text: the title-bar tab strip ("mMol" / Render Options /
      // contact maps) stands in its place.
      .navigationTitle("")
      .environmentObject(mcp)
      .onAppear { mcp.bind(engine: engine) }
  }

  // Native menus. File: Open / Fetch / Save / Export. Connect: MCP server
  // control. Buttons post notifications ContentView observes (reusing the
  // toolbar's open/save/export logic).
  @CommandsBuilder private var macCommands: some Commands {
    CommandGroup(after: .newItem) {
      Button("Open…") {
        NotificationCenter.default.post(name: .mmolOpenFile, object: nil)
      }.keyboardShortcut("o", modifiers: .command)
      Button("Fetch from PDB…") {
        NotificationCenter.default.post(name: .mmolFetch, object: nil)
      }.keyboardShortcut("o", modifiers: [.command, .shift])
      Divider()
      // ⌘S overwrites the open .pse with no panel (Save As if never saved).
      Button("Save Session") {
        NotificationCenter.default.post(name: .mmolSaveSession, object: nil)
      }.keyboardShortcut("s", modifiers: .command)
      // ⇧⌘S always shows the Save panel and updates the tracked document.
      Button("Save Session As…") {
        NotificationCenter.default.post(name: .mmolSaveSessionAs, object: nil)
      }.keyboardShortcut("s", modifiers: [.command, .shift])
      // Saves the image the last `ray` left on screen.
      Button("Export Image…") {
        NotificationCenter.default.post(name: .mmolExportImage, object: nil)
      }.keyboardShortcut("e", modifiers: [.command, .shift])
      Divider()
      Button("Clear Session") {
        NotificationCenter.default.post(name: .mmolClearSession, object: nil)
      }
    }
    // Edit ▸ New Selection… — the selection builder sheet.
    CommandGroup(after: .pasteboard) {
      Divider()
      Button("New Selection…") { engine.showSelectionBuilder = true }
    }
    // TOOLS — every exclusive interaction mode, in one menu.
    //
    // These are mutually exclusive: entering any one leaves the others
    // (see PyMOLEngine's clearing setters and `exitActiveInteractionMode`).
    // One menu makes the exclusion set legible and puts the key
    // equivalents on screen together. Shortcuts come from AppShortcuts (#360).
    CommandMenu("Tools") {
      Button("Align Structures…") { engine.openTool(.align) }
      Button("Mutagenesis…") { engine.openTool(.mutagenesis) }
      Button("Claude Code") { engine.openClaude() }
      Divider()
      // The viewport tools first: these three reinterpret a click or a
      // drag, which is what this menu was called "Mouse" for.
      Button(engine.interactionMode == .move ? "Stop Moving Objects" : "Move Objects") {
        engine.setInteractionMode(engine.interactionMode == .move ? .viewing : .move)
      }
      .keyboardShortcut(AppShortcuts.moveTool)
      Button(engine.measureMode != nil ? "Stop Measuring" : "Measure Distances") {
        engine.setMeasureMode(engine.measureMode == nil ? .distance : nil)
      }
      .keyboardShortcut(AppShortcuts.measureTool)
      // ⌃S — S for "select" (#358).
      Button(engine.interactionMode == .boxSelect ? "Stop Box Select" : "Box Select") {
        engine.toggleBoxSelect()
      }
      .keyboardShortcut(AppShortcuts.boxSelect)
    }
    // @CommandsBuilder takes at most 10 children; past that, wrap menus
    // in a `Group`.
    Group {
      // Display — desktop PyMOL 2.5's menu of that name.
      CommandMenu("Display") {
        Toggle(
          "Sequence",
          isOn: Binding(
            get: { engine.sequenceVisible },
            set: { engine.sequenceVisible = $0 })
        )
        .keyboardShortcut(AppShortcuts.sequencePane)
        Divider()
        Menu("Background") {
          ForEach(
            [("Black", "black"), ("White", "white"), ("Light Grey", "grey90"), ("Grey", "grey50")],
            id: \.1
          ) { label, color in
            Button(label) { engine.runCommand("bg_color \(color)") }
          }
        }
        Button("Render Options…") { engine.openTool(.render) }
        Button("Viewer Settings…") { engine.openTool(.settings) }
      }
      // Scene — desktop PyMOL 2.5's menu of that name. PgUp / PgDn
      // already step through the scenes (PyMOL's own key bindings).
      CommandMenu("Scene") {
        Button("Store New Scene") { engine.runCommand("scene new, store") }
        Button("Update Current Scene") { engine.runCommand("scene auto, update") }
          .disabled(engine.currentScene.isEmpty)
        Button("Delete Current Scene") { engine.runCommand("scene auto, delete") }
          .disabled(engine.currentScene.isEmpty)
        Button("Rename…") { engine.pendingSceneRename = engine.currentScene }
          .disabled(engine.currentScene.isEmpty)
        if !engine.sceneNames.isEmpty {
          Divider()
          ForEach(engine.sceneNames, id: \.self) { name in
            Toggle(
              name,
              isOn: Binding(
                get: { engine.currentScene == name },
                set: { _ in engine.recallScene(name) }))
          }
        }
        Divider()
        Button("Previous") { engine.runCommand("scene auto, previous") }
          .disabled(engine.sceneNames.isEmpty)
        Button("Next") { engine.runCommand("scene auto, next") }
          .disabled(engine.sceneNames.isEmpty)
      }
      CommandMenu("Movie") {
        Button("Movie Builder") {
          NotificationCenter.default.post(name: .mmolToggleTimeline, object: nil)
        }.keyboardShortcut("m", modifiers: [.command, .option])
      }
    }
    CommandMenu("Connect") {
      Toggle(
        "Enable AI control",
        isOn: Binding(
          get: { mcp.isRunning }, set: { _ in mcp.toggle() })
      )
      .keyboardShortcut("m", modifiers: [.control, .command])
      Divider()
      if mcp.isRunning, let port = mcp.port {
        Text("Listening on 127.0.0.1:\(port)")
        Text("Clients: \(mcp.clientCount)")
        Divider()
      }
      Button("Connect an AI app…") {
        NotificationCenter.default.post(name: .mcpOpenConnectSheet, object: nil)
      }
      Divider()
      Button("Copy connection details") {
        if let port = mcp.port {
          let s =
            "URL: http://127.0.0.1:\(port)/mcp\n"
            + "Authorization: Bearer \(mcp.token)"
          NSPasteboard.general.clearContents()
          NSPasteboard.general.setString(s, forType: .string)
        }
      }.disabled(!mcp.isRunning)
    }
  }
}

/// `s` with every character but letters and digits replaced by "_", for
/// object names and the PyMOL color names derived from them.
func sanitizedName(_ s: String) -> String {
  String(s.map { $0.isLetter || $0.isNumber ? $0 : "_" })
}

// Load a file the OS handed to mMol (Finder double-click / "Open With"). The OS
// may deliver the URL before the engine has finished initializing (cold launch
// from a file), so retry on the main queue until the engine is ready (capped so a
// failed init never loops forever). The URL may be security-scoped, so copy it
// into the temp dir before handing the path to PyMOL, which infers the format
// from the extension. The object name is the sanitized stem.
@MainActor
func loadOpenedFile(_ url: URL, into engine: PyMOLEngine, attempt: Int = 0) {
  guard engine.isReady else {
    guard attempt < 40 else { return }  // ~10s cap (40 × 250ms)
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.25) {
      loadOpenedFile(url, into: engine, attempt: attempt + 1)
    }
    return
  }
  guard confirmReplaceSessionIfNeeded(opening: url, engine: engine) else { return }
  let scoped = url.startAccessingSecurityScopedResource()
  defer { if scoped { url.stopAccessingSecurityScopedResource() } }
  let ext = url.pathExtension.isEmpty ? "pdb" : url.pathExtension
  let temp = FileManager.default.temporaryDirectory
    .appendingPathComponent("open_\(UUID().uuidString.prefix(8)).\(ext)")
  try? FileManager.default.removeItem(at: temp)
  let path: String
  if (try? FileManager.default.copyItem(at: url, to: temp)) != nil {
    path = temp.path
  } else {
    path = url.path  // fall back to the original path (e.g. local macOS file)
  }
  var name = sanitizedName(url.deletingPathExtension().lastPathComponent)
  if name.isEmpty { name = "mol" }
  engine.loadStructure(path: path, name: name)
  // Track an opened .pse as the current document so ⌘S overwrites it.
  engine.currentSessionURL = (ext.lowercased() == "pse") ? url : nil
}

// Opening a session file REPLACES the whole current session — PyMOL sessions are
// app-global and there is no multi-window support yet (#29), so a Finder
// double-click on a .pse silently wiped whatever was on screen (issue #349).
// True when that destructive replace is about to happen: a session file is being
// opened while objects are loaded. Pure and engine-free so it's unit-testable
// (same pattern as handleOpenedURLs below).
func openWouldReplaceSession(_ url: URL, hasObjects: Bool) -> Bool {
  hasObjects && PyMOLEngine.isSessionFile(url.path)
}

// Stopgap for issue #349 until sessions can open in their own window/tab: before
// a .pse/.psw replaces a non-empty session, offer to save the current one first,
// replace it outright, or cancel the open. Returns false when the open must not
// proceed (Cancel, or a failed/cancelled save).
@MainActor
func confirmReplaceSessionIfNeeded(opening url: URL, engine: PyMOLEngine) -> Bool {
  guard openWouldReplaceSession(url, hasObjects: !engine.objects.isEmpty)
  else { return true }
  let alert = NSAlert()
  alert.messageText = "Replace the current session?"
  alert.informativeText = """
    Opening “\(url.lastPathComponent)” will replace everything in the current \
    session. mMol can’t open sessions in separate windows yet.
    """
  alert.alertStyle = .warning
  alert.addButton(withTitle: "Save and Replace…")
  alert.addButton(withTitle: "Replace")
  alert.addButton(withTitle: "Cancel")
  switch alert.runModal() {
  case .alertFirstButtonReturn:
    return saveCurrentSessionForReplace(engine: engine)
  case .alertSecondButtonReturn:
    return true
  default:
    return false
  }
}

// Save the outgoing session before it's replaced, mirroring ⌘S semantics: an open
// document is overwritten silently; an untitled session shows the Save panel.
// Ordering is safe without waiting — cmd.save and the subsequent `load` both run
// synchronously through the bridge on the main actor, in issue order. Returns
// false when the user cancels the panel, which aborts the replace.
@MainActor
private func saveCurrentSessionForReplace(engine: PyMOLEngine) -> Bool {
  let dest: URL
  if let current = engine.currentSessionURL {
    dest = current
  } else {
    let panel = NSSavePanel()
    if let pse = UTType(filenameExtension: "pse") { panel.allowedContentTypes = [pse] }
    panel.nameFieldStringValue = "session.pse"
    panel.canCreateDirectories = true
    panel.title = "Save Session"
    guard panel.runModal() == .OK, let url = panel.url else { return false }
    dest = url
  }
  engine.saveSession(to: dest)
  return true
}

// Load EVERY URL the OS handed us in one open (multi-file Terminal `open`, Finder
// "Open With" multi-select, Dock-icon drop). Factored out of MMolAppDelegate so
// the "load them all, not just the first" contract (issue #222) is unit-testable
// without booting the engine: tests pass a spy for `load` (which defaults to
// loadOpenedFile) and assert every URL is forwarded, in order.
@MainActor
func handleOpenedURLs(
  _ urls: [URL], into engine: PyMOLEngine,
  load: @MainActor (URL, PyMOLEngine) -> Void = { loadOpenedFile($0, into: $1) }
) {
  for url in urls { load(url, engine) }
}
