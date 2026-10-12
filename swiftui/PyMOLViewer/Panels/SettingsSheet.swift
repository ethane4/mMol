// SettingsSheet.swift — the searchable list of every PyMOL setting.

import SwiftUI

// MARK: - Searchable Settings panel (all ~825 PyMOL settings)

struct SettingsSheet: View {
  @EnvironmentObject var engine: PyMOLEngine
  @Environment(\.dismiss) private var dismiss
  @State private var search = ""
  private var filtered: [SettingItem] {
    let q = search.trimmingCharacters(in: .whitespaces).lowercased()
    if q.isEmpty { return engine.settingsCatalog }
    return engine.settingsCatalog.filter { $0.name.contains(q) }
  }

  var body: some View {
    VStack(spacing: 0) {
      HStack {
        Text("Settings").font(.headline)
        Spacer()
        Button("Done") { dismiss() }
      }.padding(16)

      HStack {
        Image(systemName: "magnifyingglass").foregroundStyle(.secondary)
        TextField("Search \(engine.settingsCatalog.count) settings…", text: $search)
          .textFieldStyle(.plain)
          .autocorrectionDisabled()
        if !search.isEmpty {
          Button {
            search = ""
          } label: {
            Image(systemName: "xmark.circle.fill")
          }
          .buttonStyle(.plain).foregroundStyle(.secondary)
        }
      }
      .padding(8).background(Color.gray.opacity(0.15))
      .clipShape(RoundedRectangle(cornerRadius: 8))
      .padding(.horizontal, 16)

      if engine.settingsCatalog.isEmpty {
        Spacer()
        ProgressView("Loading settings…")
        Spacer()
      } else {
        List(filtered) { item in
          SettingRow(item: item)
        }
        .listStyle(.plain)
      }
    }
    .onAppear {
      if engine.settingsCatalog.isEmpty { engine.loadSettingsCatalog() }
    }
    .frame(width: 460, height: 560)
  }

}

private struct SettingRow: View {
  let item: SettingItem
  @EnvironmentObject var engine: PyMOLEngine
  @State private var text = ""

  private var isBool: Bool { item.type == 1 }
  private var boolOn: Bool {
    if let d = Double(item.val) { return d != 0 }
    return ["on", "true", "yes"].contains(item.val.lowercased())
  }

  var body: some View {
    HStack(spacing: 8) {
      Text(item.name)
        .font(.system(size: 12, design: .monospaced))
        .lineLimit(1).minimumScaleFactor(0.7)
      Spacer(minLength: 8)
      if isBool {
        Toggle(
          "",
          isOn: Binding(
            get: { boolOn },
            set: { engine.setSetting(item.name, $0 ? "1" : "0") })
        )
        .labelsHidden().tint(TimelineTheme.accent)
      } else {
        TextField("", text: $text)
          .multilineTextAlignment(.trailing)
          .classicField()
          .frame(width: 110)
          .autocorrectionDisabled()
          .onSubmit { engine.setSetting(item.name, text) }
          .onAppear { text = item.val }
          .onChange(of: item.val) { text = $0 }
      }
    }
    .font(.system(size: 12))
  }
}
