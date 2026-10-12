// MCPStatusView.swift — the MCP server's toolbar button and menu (macOS),
// drawn like the rest of the classic PyMOL window: a flat Qt-style button that
// opens a classic menu (ClassicMenu.swift) rather than a system popover.
import SwiftUI

struct MCPStatusView: View {
  @EnvironmentObject var mcp: MCPServerManager
  @EnvironmentObject var engine: PyMOLEngine

  private var dotColor: Color {
    if !mcp.isRunning { return Color(white: 0.5) }
    return mcp.clientCount > 0 ? ClassicUI.legendGreen : .yellow
  }

  var body: some View {
    ClassicMenuTrigger {
      HStack(spacing: 5) {
        Circle()
          .fill(dotColor)
          .frame(width: 8, height: 8)
          .opacity(mcp.activeTool ? 0.35 : 1.0)
          .animation(
            mcp.activeTool
              ? .easeInOut(duration: 0.6).repeatForever(autoreverses: true)
              : .default, value: mcp.activeTool)
        Text(mcp.isRunning ? "mMol MCP" : "MCP off")
        if mcp.clientCount > 0 {
          Text("· \(mcp.clientCount)").foregroundColor(Color(white: 0.6))
        }
      }
      .padding(.horizontal, 10)
      .frame(height: 22)
      // Fully rounded ends, echoing the window corner it sits beside.
      .classicButtonFace(cornerRadius: 11)
    } items: {
      menuItems()
    }
    .help("mMol MCP server")
  }

  private func menuItems() -> [ClassicMenuItem] {
    var items: [ClassicMenuItem] = [
      .check("Enable MCP server", mcp.isRunning) { mcp.toggle() },
      .separator,
    ]
    if mcp.isRunning, let port = mcp.port {
      items.append(.disabled("Listening on 127.0.0.1:\(port)"))
      items.append(.disabled("Clients connected: \(mcp.clientCount)"))
      if !mcp.lastAction.isEmpty {
        items.append(.disabled("Last action: \(Self.short(mcp.lastAction))"))
      }
    } else {
      items.append(.disabled("Server is off"))
    }
    if !mcp.activityLog.isEmpty {
      items.append(
        .submenu("Recent activity", mcp.activityLog.suffix(12).map { .disabled(Self.short($0)) }))
    }
    items.append(.separator)
    // One-click connect for the two common clients; both trust the session.
    // Their result goes to the console, where the popover used to show it.
    if mcp.isRunning {
      items.append(
        .action("Connect Claude Code") {
          mcp.connectClaudeCode { msg in DispatchQueue.main.async { engine.logLine(msg) } }
        })
      items.append(
        .action("Connect Claude App") {
          mcp.noteUserInitiatedConnect()
          engine.logLine(MCPDesktopInstaller.installViaConfig().message)
        })
    } else {
      items += [.disabled("Connect Claude Code"), .disabled("Connect Claude App")]
    }
    items.append(
      .action("Connect an AI app…") {
        NotificationCenter.default.post(name: .mcpOpenConnectSheet, object: nil)
      })
    return items
  }

  /// `text` cut to one menu row's worth.
  private static func short(_ text: String) -> String {
    text.count > 60 ? String(text.prefix(59)) + "…" : text
  }
}
