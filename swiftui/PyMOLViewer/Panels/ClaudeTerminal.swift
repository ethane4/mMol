// ClaudeTerminal.swift — `claude` at the mMol> prompt (or Tools ▸ Claude Code):
// the user's own Claude Code in a terminal docked under the viewport, wired to
// this window's MCP server so its mmol tools drive the live session.
//
// The session outlives tab switches (the terminal view is kept here); closing
// the tab ends it. The MCP server is started if needed and this one connection
// is trusted without the approval prompt, since the user launched it.

import SwiftTerm
import SwiftUI

@MainActor
final class ClaudeSession: NSObject, ObservableObject, LocalProcessTerminalViewDelegate {
  static let shared = ClaudeSession()

  /// Bumped for each new session so SwiftUI swaps in its terminal view.
  @Published private(set) var generation = 0
  private(set) var terminal: LocalProcessTerminalView?
  private var running = false
  /// Called when Claude Code exits normally (/exit), to close its tab.
  var onExit: (() -> Void)?

  private static let prompt = """
    You are running inside mMol, a desktop molecular viewer (a PyMOL fork). The mmol MCP \
    tools drive the user's live mMol session in the window around this terminal: \
    run_pymol_command, run_python, get_session_state, capture_viewport, search_pdb, \
    color_by_csv, contact_map. Start with get_session_state; after visual changes, \
    check them with capture_viewport. Keep changes to what the user asks for.
    """

  var hasFocus: Bool {
    guard let t = terminal, let r = t.window?.firstResponder as? NSView else { return false }
    return r === t || r.isDescendant(of: t)
  }

  /// Claude's working folder (its files, e.g. saved images, land here). Its own
  /// folder rather than home, so Claude Code's trust prompt covers only this.
  private static var workspace: String {
    let dir = FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
      .appendingPathComponent("mMol/Claude", isDirectory: true)
    try? FileManager.default.createDirectory(at: dir, withIntermediateDirectories: true)
    return dir.path
  }

  func start() {
    guard !running else { return }
    let mcp = MCPServerManager.shared
    if !mcp.isRunning { mcp.start() }
    waitForServer(tries: 50)
  }

  func stop() {
    if running { terminal?.process.terminate() }
    running = false
    terminal = nil
  }

  private func waitForServer(tries: Int) {
    let mcp = MCPServerManager.shared
    if mcp.isRunning, let port = mcp.port { return launch(port: port, token: mcp.token) }
    guard tries > 0 else {
      newTerminal().feed(text: "Couldn't start mMol's MCP server, so Claude can't connect.\r\n")
      return
    }
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.1) { self.waitForServer(tries: tries - 1) }
  }

  private func newTerminal() -> LocalProcessTerminalView {
    let t = LocalProcessTerminalView(frame: NSRect(x: 0, y: 0, width: 800, height: 300))
    t.processDelegate = self
    t.font = .monospacedSystemFont(ofSize: 12, weight: .regular)
    t.nativeBackgroundColor = NSColor(white: 0.08, alpha: 1)
    t.nativeForegroundColor = NSColor(white: 0.9, alpha: 1)
    terminal = t
    generation += 1
    return t
  }

  private func launch(port: Int, token: String) {
    // The server's address and token, readable only by this user.
    let config = [
      "mcpServers": [
        "mmol": [
          "type": "http", "url": "http://127.0.0.1:\(port)/mcp",
          "headers": ["Authorization": "Bearer \(token)"],
        ]
      ]
    ]
    let path = (NSTemporaryDirectory() as NSString).appendingPathComponent("mmol_claude_mcp.json")
    FileManager.default.createFile(
      atPath: path, contents: try? JSONSerialization.data(withJSONObject: config),
      attributes: [.posixPermissions: 0o600])
    MCPServerManager.shared.noteUserInitiatedConnect()

    // The app's own Python/PyMOL variables must not leak into Claude's shell.
    var env = ProcessInfo.processInfo.environment.filter { k, _ in
      !(k.hasPrefix("PYTHON") || k.hasPrefix("PYMOL") || k.hasPrefix("MMOL")
        || k == "SSL_CERT_FILE")
    }
    env["TERM"] = "xterm-256color"
    env["COLORTERM"] = "truecolor"
    env["MMOL_MCP_CONFIG"] = path
    env["MMOL_PROMPT"] = Self.prompt
    // A login, interactive zsh so PATH matches the user's Terminal.
    let script = """
      command -v claude >/dev/null || { printf 'Claude Code (claude) is not installed or not on your PATH.\\n\
      Install it from https://docs.claude.com/en/docs/claude-code and type claude again.\\n'; exit 127; }
      exec claude --mcp-config "$MMOL_MCP_CONFIG" --allowedTools mcp__mmol --append-system-prompt "$MMOL_PROMPT"
      """
    running = true
    newTerminal().startProcess(
      executable: "/bin/zsh", args: ["-l", "-i", "-c", script],
      environment: env.map { "\($0)=\($1)" }, execName: "zsh",
      currentDirectory: Self.workspace)
  }

  // MARK: LocalProcessTerminalViewDelegate

  func processTerminated(source: TerminalView, exitCode: Int32?) {
    running = false
    // A normal exit closes the tab; anything else (Claude Code missing, a
    // crash) leaves the terminal up so its message can be read.
    if exitCode == 0, let onExit {
      DispatchQueue.main.async { onExit() }
      return
    }
    source.feed(
      text: "\r\n[Claude session ended. Type claude at the mMol> prompt to start a new one.]\r\n")
  }

  func processFailedToStart(source: TerminalView, error: LocalProcessError) {
    running = false
    source.feed(text: "Couldn't start the shell: \(error)\r\n")
  }

  func sizeChanged(source: LocalProcessTerminalView, newCols: Int, newRows: Int) {}
  func setTerminalTitle(source: LocalProcessTerminalView, title: String) {}
  func hostCurrentDirectoryUpdate(source: TerminalView, directory: String?) {}
}

extension PyMOLEngine {
  /// `claude` at the prompt / Tools ▸ Claude Code: show the tab, starting a session if none runs.
  @MainActor func openClaude() {
    openTool(.claude)
    ClaudeSession.shared.onExit = { [weak self] in self?.closeTool(.claude) }
    ClaudeSession.shared.start()
  }
}

/// The docked terminal (shown under the viewport while the Claude tab is active).
struct ClaudeTerminalView: View {
  @ObservedObject private var session = ClaudeSession.shared

  var body: some View {
    if let t = session.terminal {
      TerminalHost(terminal: t).id(session.generation)
    } else {
      Color(white: 0.08)
    }
  }

  private struct TerminalHost: NSViewRepresentable {
    let terminal: LocalProcessTerminalView
    func makeNSView(context: Context) -> LocalProcessTerminalView {
      DispatchQueue.main.async { terminal.window?.makeFirstResponder(terminal) }
      return terminal
    }
    func updateNSView(_ v: LocalProcessTerminalView, context: Context) {}
  }
}
