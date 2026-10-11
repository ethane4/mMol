// ViewerSettingsPanel.swift — Display ▸ Viewer Settings….
//
// A title-bar tab like Render Options, but over the user's own scene: the
// viewport stays live while the right column shows the viewer settings
// (SceneCard) and its "All Settings…" button opens the searchable SettingsSheet.

import SwiftUI

struct ViewerSettingsPanel: View {
  var body: some View {
    ScrollView {
      SceneCard().padding(10)
    }
    .background(ClassicUI.panelBg)
    .classicControls()
  }
}
