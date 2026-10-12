// MCPDrivingIndicator.swift — "Claude is controlling mMol" (macOS): a
// non-blocking capsule in the toolbar, left of the MCP button. It used to be a
// banner across the top of the window, which shifted the whole layout down
// every time it appeared.
import SwiftUI

struct MCPDrivingIndicator: View {
  @EnvironmentObject var mcp: MCPServerManager

  var body: some View {
    HStack(spacing: 6) {
      Circle().fill(ClassicUI.legendGreen).frame(width: 8, height: 8)
      // Only the action shortens when space is tight, never the label.
      Text("Claude is controlling mMol").fixedSize()
      if !mcp.lastAction.isEmpty {
        Text("· \(mcp.lastAction)")
          .foregroundColor(Color(white: 0.6))
          .truncationMode(.tail)
          .frame(maxWidth: 220, alignment: .leading)
      }
      Button("Stop") { mcp.stop() }
        .buttonStyle(.plain)
        .foregroundColor(ClassicMenuStyle.destructiveText)
        .help("Stop the MCP server")
    }
    .padding(.horizontal, 10)
    .frame(height: 22)
    .classicButtonFace(cornerRadius: 11)
  }
}
