// SelectionBuilderSheet.swift — Edit ▸ New Selection…

import SwiftUI

// MARK: - Selection builder (named selections + spatial algebra)

struct SelectionBuilderSheet: View {
  @EnvironmentObject var engine: PyMOLEngine
  @Environment(\.dismiss) private var dismiss

  enum Op: String, CaseIterable, Identifiable {
    case none = "base only"
    case within = "within … of"
    case around = "around"
    case
      expand = "expand"
    case extend = "extend (bonds)"
    case and = "and"
    case or = "or"
    case not = "not"
    var id: String { rawValue }
    var needsDist: Bool { self == .within || self == .around || self == .expand }
    var needsCount: Bool { self == .extend }
    var needsOther: Bool { self == .within || self == .and || self == .or }
  }

  @State private var base = "sele"
  @State private var op: Op = .around
  @State private var dist = "5"
  @State private var other = "all"
  @State private var byres = true
  @State private var name = "sel01"
  @State private var previewWork: DispatchWorkItem?

  @State private var showRename = false
  @State private var renameTarget = ""
  @State private var renameText = ""

  private var bases: [String] {
    ["sele", "all", "polymer", "organic", "solvent"] + engine.objects.map { $0.name }
  }
  private var selections: [ObjectEntry] { engine.objects.filter { $0.isSelection } }

  var body: some View {
    VStack(spacing: 0) {
      HStack {
        Text("Selections").font(.headline)
        Spacer()
        Button("Done") { dismiss() }
      }.padding(16)

      ScrollView {
        VStack(alignment: .leading, spacing: 16) {
          builder
          Divider()
          manage
        }.padding(16)
      }
    }
    .onChange(of: expr) { _ in schedulePreview() }
    .onAppear { schedulePreview() }
    .onDisappear { engine.selectionPreviewCount = nil }
    .alert("Rename “\(renameTarget)”", isPresented: $showRename) {
      TextField("New name", text: $renameText)
      Button("Rename") { engine.renameObject(renameTarget, to: renameText) }
      Button("Cancel", role: .cancel) {}
    }
    .frame(width: 440, height: 540)
  }

  // MARK: builder

  private var expr: String {
    let b = "(\(base))"
    var e: String
    switch op {
    case .none: e = b
    case .within: e = "\(b) within \(distNum) of (\(other))"
    case .around: e = "\(b) around \(distNum)"
    case .expand: e = "\(b) expand \(distNum)"
    case .extend: e = "\(b) extend \(Int(distNum.rounded()))"
    case .and: e = "\(b) and (\(other))"
    case .or: e = "\(b) or (\(other))"
    case .not: e = "not \(b)"
    }
    if byres { e = "byres (\(e))" }
    return e
  }
  private var distNum: Double { Double(dist) ?? 5 }

  private var builder: some View {
    VStack(alignment: .leading, spacing: 12) {
      Text("New selection").font(.system(size: 13, weight: .semibold))
      row("From") { picker($base, bases) }
      row("Operator") {
        Picker("", selection: $op) { ForEach(Op.allCases) { Text($0.rawValue).tag($0) } }
          .labelsHidden()
      }
      if op.needsDist || op.needsCount {
        row(op.needsCount ? "Bonds" : "Distance (Å)") {
          TextField("5", text: $dist)
            .frame(width: 70).textFieldStyle(.roundedBorder)
        }
      }
      if op.needsOther {
        row("Of") { picker($other, bases) }
      }
      Toggle("Whole residues (byres)", isOn: $byres).tint(TimelineTheme.accent)

      // Live preview of the composed expression + match count.
      VStack(alignment: .leading, spacing: 4) {
        Text(expr).font(.system(size: 11, design: .monospaced))
          .foregroundStyle(.secondary).textSelection(.enabled)
        Text(previewText).font(.caption).foregroundStyle(.secondary)
      }
      HStack {
        TextField("name", text: $name).frame(width: 120).textFieldStyle(.roundedBorder)
        Spacer()
        Button {
          engine.createSelection(name: name, expr: expr)
          dismiss()
        } label: {
          Label("Create", systemImage: "plus").font(.system(size: 13, weight: .semibold))
        }
        .buttonStyle(.borderedProminent).tint(TimelineTheme.accent)
        .disabled(name.trimmingCharacters(in: .whitespaces).isEmpty)
      }
    }
  }

  private var previewText: String {
    if let c = engine.selectionPreviewCount { return "\(c) atom\(c == 1 ? "" : "s")" }
    return "…"
  }

  private func schedulePreview() {
    previewWork?.cancel()
    let e = expr
    let work = DispatchWorkItem { engine.previewSelection(e) }
    previewWork = work
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.3, execute: work)
  }

  // MARK: manage existing selections

  private var manage: some View {
    VStack(alignment: .leading, spacing: 8) {
      Text("Manage").font(.system(size: 13, weight: .semibold))
      if selections.isEmpty {
        Text("No named selections yet.").font(.caption).foregroundStyle(.secondary)
      } else {
        ForEach(selections) { sel in
          HStack(spacing: 10) {
            Toggle(
              "",
              isOn: Binding(
                get: { sel.isEnabled },
                set: { on in
                  engine.runCommand("\(on ? "enable" : "disable") \(sel.name)", naming: sel.name)
                })
            )
            .labelsHidden()
            Text(sel.name).font(.system(size: 13))
            if let c = sel.atomCount { Text("(\(c))").font(.caption).foregroundStyle(.secondary) }
            Spacer()
            Button {
              renameTarget = sel.name
              renameText = sel.name
              showRename = true
            } label: {
              Image(systemName: "pencil")
            }.buttonStyle(.borderless)
            Button(role: .destructive) {
              engine.runCommand("delete \(sel.name)", naming: sel.name)
            } label: {
              Image(systemName: "trash")
            }.buttonStyle(.borderless)
          }
        }
      }
    }
  }

  // MARK: helpers

  @ViewBuilder
  private func row<C: View>(_ label: String, @ViewBuilder _ content: () -> C) -> some View {
    HStack {
      Text(label).font(.system(size: 12)).frame(width: 90, alignment: .leading)
      content()
      Spacer()
    }
  }

  private func picker(_ sel: Binding<String>, _ opts: [String]) -> some View {
    ClassicPicker(selection: sel, options: opts.map { ($0, $0) })
  }
}
