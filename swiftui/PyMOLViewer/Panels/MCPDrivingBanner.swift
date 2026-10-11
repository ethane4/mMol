// MCPDrivingBanner.swift — non-blocking "AI is driving" banner (macOS).
import SwiftUI

struct MCPDrivingBanner: View {
  @EnvironmentObject var mcp: MCPServerManager

  var body: some View {
    Group {
      if mcp.activeTool {
        HStack(spacing: 10) {
          Circle().fill(Color.green).frame(width: 8, height: 8)
          Text("Claude is controlling mMol")
            .font(.caption).fontWeight(.semibold)
          if !mcp.lastAction.isEmpty {
            Text("· \(mcp.lastAction)")
              .font(.caption).foregroundStyle(.secondary).lineLimit(1)
          }
          Spacer()
          Button("Stop") { mcp.stop() }
            .controlSize(.small).tint(.red)
        }
        .padding(.horizontal, 12).padding(.vertical, 6)
        .frame(maxWidth: .infinity)
        .background(Color.green.opacity(0.12))
        .overlay(
          Rectangle().frame(height: 1)
            .foregroundStyle(Color.green.opacity(0.3)), alignment: .bottom
        )
        .transition(.move(edge: .top).combined(with: .opacity))
      }
    }
    .animation(.easeInOut(duration: 0.2), value: mcp.activeTool)
  }
}
