// ToolTabs.swift — the title-bar tabs beside "mMol" that swap the right column
// for a tool (Render Options, Settings, Movie, Align, Mutagenesis) or dock a
// panel under the viewport (Movie's timeline, Claude's terminal). Plot tabs
// (contact maps, Ramachandran) swap the viewport instead: engine.contactMapTabs.

import SwiftUI

enum ToolTab: String {
  case render = "Render Options"
  case settings = "Settings"
  case movie = "Movie"
  case
    align = "Align"
  case mutagenesis = "Mutagenesis"
  case claude = "Claude"
}

extension PyMOLEngine {
  var renderTabActive: Bool { activeTool == .render }
  var movieTabActive: Bool { activeTool == .movie }

  func openTool(_ t: ToolTab) {
    if !toolTabs.contains(t) { toolTabs.append(t) }
    activeContactMap = nil
    activeTool = t
  }

  @MainActor func closeTool(_ t: ToolTab) {
    toolTabs.removeAll { $0 == t }
    if activeTool == t { activeTool = nil }
    switch t {
    case .render: runPython("from mmol import render_options as _ro\n_ro.end(close=True)")
    case .mutagenesis: endMutagenesis()
    case .claude:
      ClaudeSession.shared.stop()
      if MCPServerManager.shared.isRunning { MCPServerManager.shared.stop() }
    default: break
    }
  }

  /// The mMol tab (nil) or a plot tab.
  func showTab(plot: String? = nil) {
    activeTool = nil
    activeContactMap = plot
  }
}
