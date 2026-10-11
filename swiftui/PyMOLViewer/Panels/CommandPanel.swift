// CommandPanel.swift — Log viewer + command input for PyMOL
// Replaces modules/pymol/appkit_command_panel.py with pure SwiftUI.

import SwiftUI

struct CommandPanel: View {
  @EnvironmentObject var engine: PyMOLEngine

  @State private var commandText = ""
  @State private var commandHistory: [String] = []
  @State private var historyIndex = -1
  @State private var hint: ConsoleHint?
  @State private var pick: Int?  // the choice Tab is cycling through

  private var bgColor: Color { ClassicUI.panelBg }
  private var textColor: Color { ClassicUI.consoleText.color }

  var body: some View {
    VStack(spacing: 0) {
      // Scrolling log area
      ConsoleLogView(log: engine.console)

      Divider()
        .background(Color.gray.opacity(0.4))

      // Command input bar
      HStack(spacing: 4) {
        Text("mMol>")
          .font(ClassicUI.font)
          .foregroundColor(textColor)

        CommandTextField(
          text: $commandText,
          textColor: textColor,
          bgColor: bgColor,
          fontSize: 12,
          onSubmit: submitCommand,
          onUpArrow: historyBack,
          onDownArrow: historyForward,
          onComplete: { _ in tabComplete() }
        )
      }
      .padding(.horizontal, 8)
      .padding(.vertical, 6)
      .background(bgColor)
      if let hint, !hint.cmd.isEmpty || !hint.choices.isEmpty {
        hintRow(hint)
      }
    }
    .background(bgColor)
    .onAppear {
      // Test affordance: type PYMOL_CONSOLE_TEXT at the prompt once ready.
      if let t = ProcessInfo.processInfo.environment["PYMOL_CONSOLE_TEXT"] {
        DispatchQueue.main.asyncAfter(deadline: .now() + 8) { commandText = t }
      }
    }
    .onChange(of: commandText) { text in
      // Typing (not Tab cycling) asks for a fresh hint.
      if let h = hint, let p = pick, text == h.applying(h.choices[p], to: text, final: false) {
        return
      }
      pick = nil
      hint = engine.consoleHint(text)
    }
  }

  // MARK: - Hints

  // Fish-style help under the prompt: the signature with the argument being
  // typed in the accent, what it means (or what the command does), the
  // completions with their descriptions (click one; Tab cycles), an example.
  private func hintRow(_ h: ConsoleHint) -> some View {
    VStack(alignment: .leading, spacing: 2) {
      if !h.cmd.isEmpty {
        signature(h).lineLimit(1)
        let info = h.arg >= 0 ? "\(h.params[h.arg].name): \(h.argdoc)" : h.about
        if !info.isEmpty, !info.hasSuffix(": ") {
          Text(info).foregroundColor(textColor.opacity(0.6)).lineLimit(1)
        }
      }
      let shown = h.choices.filter { $0.text != h.word(in: commandText) }
      if !shown.isEmpty {
        LazyVGrid(
          columns: Array(repeating: GridItem(.flexible(), alignment: .leading), count: 3),
          alignment: .leading, spacing: 1
        ) {
          ForEach(Array(shown.prefix(9).enumerated()), id: \.offset) { _, c in
            let i = h.choices.firstIndex { $0.text == c.text }
            (Text(c.text).foregroundColor(i == pick ? ClassicUI.accent : textColor)
              + Text(c.desc.isEmpty ? "" : "  \(c.desc)").foregroundColor(textColor.opacity(0.5)))
              .lineLimit(1)
              .onTapGesture { commandText = h.applying(c, to: commandText, final: true) }
          }
        }
        if shown.count > 9 {
          Text("+\(shown.count - 9) more (Tab cycles)").foregroundColor(textColor.opacity(0.4))
        }
      } else if !h.example.isEmpty {
        Text("e.g.  \(h.example)").foregroundColor(textColor.opacity(0.45)).lineLimit(1)
      }
    }
    .font(ClassicUI.font)
    .foregroundColor(textColor.opacity(0.85))
    .padding(.horizontal, 8).padding(.bottom, 5)
    .frame(maxWidth: .infinity, alignment: .leading)
    .background(bgColor)
  }

  /// `fetch code [name] [state] [type=cif] …`: required plain, optional in
  /// brackets (one level, no nesting), the current one in the accent.
  private func signature(_ h: ConsoleHint) -> Text {
    h.params.enumerated().reduce(Text(h.cmd).bold()) { t, e in
      let (i, p) = e
      if p.minor == true && i != h.arg { return t }
      let def = p.default.isEmpty || p.default.count > 12 ? "" : "=\(p.default)"
      let word = p.req ? p.name : "[\(p.name)\(def)]"
      let piece =
        i == h.arg
        ? Text(word).bold().foregroundColor(ClassicUI.accent)
        : Text(word).foregroundColor(textColor.opacity(p.req ? 0.9 : 0.55))
      return t + Text(" ") + piece
    }
  }

  /// Tab: one choice → take it; several → extend to their common prefix, or,
  /// when that adds nothing, step through them.
  private func tabComplete() -> String? {
    guard let h = hint, !h.choices.isEmpty else { return nil }
    if h.choices.count == 1 { return h.applying(h.choices[0], to: commandText, final: true) }
    let texts = h.choices.map(\.text)
    let prefix = texts.dropFirst().reduce(texts[0]) { String($0.commonPrefix(with: $1)) }
    let typed = h.word(in: commandText)
    if pick == nil, prefix.count > typed.count {
      return h.applying(.init(text: prefix, desc: "", post: ""), to: commandText, final: false)
    }
    pick = ((pick ?? -1) + 1) % h.choices.count
    return h.applying(h.choices[pick!], to: commandText, final: false)
  }

  // MARK: - Actions

  private func submitCommand() {
    let trimmed = commandText.trimmingCharacters(in: .whitespaces)
    guard !trimmed.isEmpty else { return }

    commandHistory.append(trimmed)
    historyIndex = commandHistory.count
    if trimmed == "claude" {  // Claude Code in a terminal tab (ClaudeTerminal.swift)
      engine.feedbackLog.append("mMol>claude")
      engine.openClaude()
      commandText = ""
      hint = nil
      return
    }

    engine.runCommand(trimmed)

    commandText = ""
    hint = nil
  }

  private func historyBack() {
    guard !commandHistory.isEmpty, historyIndex > 0 else { return }
    historyIndex -= 1
    commandText = commandHistory[historyIndex]
  }

  private func historyForward() {
    guard !commandHistory.isEmpty else { return }
    historyIndex += 1
    if historyIndex < commandHistory.count {
      commandText = commandHistory[historyIndex]
    } else {
      historyIndex = commandHistory.count
      commandText = ""
    }
  }
}

// MARK: - Console log text view (macOS)

/// The console log rendered as ONE selectable `NSTextView` instead of a SwiftUI
/// `Text` per line (#406).
///
/// `.textSelection(.enabled)` scopes a selection to a single text view, so with
/// a `Text` per line a drag could never reach past the line it started on and
/// ⌘C could only ever copy that one line; the `LazyVStack` compounded it by
/// never materialising the scrolled-away rows. A single text view gets
/// multi-line drag selection, ⇧-click, ⌘A and ⌘C (Edit ▸ Copy) for free.
///
/// Deliberately non-editable and not a field editor: `ContentView`'s
/// `textFocusFlags()` counts only editable text views as "the user is typing",
/// so a focused log must NOT look like a text field — otherwise arrows, home,
/// end and the ctrl-letter chords would stop reaching PyMOL's key bindings for
/// as long as the console had focus.
private struct ConsoleTextView: NSViewRepresentable {
  let entries: [String]
  /// ConsoleLog.trimmed — lines the cap has dropped so far.
  let trimmed: Int

  /// Monospaced 12 pt in the console colour, on the panel grey.
  private static let attributes: [NSAttributedString.Key: Any] = {
    let base = NSFont.systemFont(ofSize: 12)
    let font =
      base.fontDescriptor.withDesign(.monospaced)
      .flatMap { NSFont(descriptor: $0, size: 12) } ?? base
    return [.font: font, .foregroundColor: ClassicUI.consoleText.nsColor]
  }()

  func makeNSView(context: Context) -> NSScrollView {
    let scrollView = ConsoleLogText.makeScrollView()
    scrollView.backgroundColor = ClassicUI.chromeBg.nsColor
    if let textView = scrollView.documentView as? NSTextView {
      textView.backgroundColor = ClassicUI.chromeBg.nsColor
      textView.typingAttributes = Self.attributes
    }
    return scrollView
  }

  func updateNSView(_ scrollView: NSScrollView, context: Context) {
    guard let textView = scrollView.documentView as? NSTextView else { return }
    let coordinator = context.coordinator

    let update = ConsoleLogText.plan(
      previous: coordinator.rendered, current: entries,
      dropped: trimmed - coordinator.renderedTrimmed)
    guard update != .none else { return }
    // Follow new output only when the user is already parked at the bottom.
    // Scrolling unconditionally (what the SwiftUI log did) would yank the
    // view away from anyone reading — or mid-drag selecting — older output
    // every time a command printed a line.
    let follow = coordinator.rendered.isEmpty || Self.isScrolledToBottom(scrollView)
    // Reading further up while the cap trims the front: measure the trimmed
    // text so the view can be shifted by exactly that much, and what the
    // user is reading stays put instead of sliding up out of view.
    var trimmedHeight: CGFloat = 0
    if !follow, case .trimThenAppend(let drop, _) = update,
      let lm = textView.layoutManager, let tc = textView.textContainer
    {
      let firstKept = NSRange(location: drop, length: 0)
      lm.ensureLayout(forCharacterRange: NSRange(location: 0, length: drop))
      trimmedHeight =
        lm.boundingRect(
          forGlyphRange: lm.glyphRange(
            forCharacterRange: firstKept,
            actualCharacterRange: nil),
          in: tc
        ).minY
    }

    ConsoleLogText.apply(update, to: textView, attributes: Self.attributes)
    coordinator.rendered = entries
    coordinator.renderedTrimmed = trimmed

    if follow {
      // Lay out and scroll in the same pass. A deferred scroll left a frame
      // in which a front-trim had shifted the text up past the old bottom,
      // showing a blank or jumped log until the scroll caught up — the
      // "output appears and disappears" flicker once the log is at its cap.
      if let tc = textView.textContainer { textView.layoutManager?.ensureLayout(for: tc) }
      textView.scrollToEndOfDocument(nil)
      // Belt and braces: once more after layout settles (the old behaviour).
      DispatchQueue.main.async { textView.scrollToEndOfDocument(nil) }
    } else if trimmedHeight > 0 {
      let clip = scrollView.contentView
      clip.scroll(
        to: NSPoint(
          x: clip.bounds.origin.x,
          y: max(0, clip.bounds.origin.y - trimmedHeight)))
      scrollView.reflectScrolledClipView(clip)
    }
  }

  func makeCoordinator() -> Coordinator { Coordinator() }

  /// Holds what the text view currently shows, so each update is a diff
  /// rather than a rebuild. SwiftUI hands `updateNSView` only the new value.
  final class Coordinator {
    var rendered: [String] = []
    var renderedTrimmed = 0
  }

  /// Within a point of the end. Slack because the clip view's bounds and the
  /// document height rarely land on exactly the same fractional value.
  private static func isScrolledToBottom(_ scrollView: NSScrollView) -> Bool {
    guard let documentHeight = scrollView.documentView?.frame.height else { return true }
    return scrollView.contentView.bounds.maxY >= documentHeight - 1
  }
}

// MARK: - Console log text model (macOS)

/// The edit that turns what the console log currently displays into what a new
/// `feedbackLog` value says it should display.
enum ConsoleLogUpdate: Equatable {
  /// Nothing changed.
  case none
  /// Rewrite the whole document — first fill, or a log that was cleared or
  /// rewritten rather than extended.
  case replace(String)
  /// Append at the end. The common case, and the only one that leaves the
  /// user's selection and scroll position completely untouched.
  case append(String)
  /// The log ALSO lost leading lines — `PyMOLEngine` trims the front once the
  /// feedback log passes its 400-line cap. Delete that many UTF-16 units from
  /// the start, then append (which may be empty).
  case trimThenAppend(dropCharacters: Int, append: String)
}

/// Turns `feedbackLog` (an array of lines) into edits for a single text view,
/// and owns the text view those edits land in.
///
/// Split out of `ConsoleTextView` so it can be tested: a test cannot drag-select
/// in a live app, but it CAN plan an edit, apply it to a real text view, select
/// all and copy — which is exactly what #406 is about. The interesting case is
/// the one that is invisible until it bites: a log trimmed at the front must not
/// be mistaken for a log that only grew at the back.
enum ConsoleLogText {

  /// One line per line: the separator IS what makes ⌘C paste as multiple
  /// lines rather than one run-on string.
  static let lineSeparator = "\n"

  static func text(for entries: some Sequence<String>) -> String {
    entries.joined(separator: lineSeparator)
  }

  /// The log's scroll view + text view, configured but unstyled. A factory
  /// rather than inline setup in `makeNSView` so a test can hold the real
  /// thing: an `NSViewRepresentable`'s `Context` cannot be built outside
  /// SwiftUI, which would otherwise put every AppKit-side guarantee here —
  /// selectable, NOT editable, `copy:` reachable — out of reach of a test.
  static func makeScrollView() -> NSScrollView {
    let scrollView = NSTextView.scrollableTextView()
    scrollView.hasVerticalScroller = true
    scrollView.drawsBackground = true

    guard let textView = scrollView.documentView as? NSTextView else { return scrollView }
    textView.isEditable = false  // see ConsoleTextView: keeps textFocusFlags() honest
    textView.isSelectable = true
    textView.isRichText = false
    textView.allowsUndo = false
    textView.drawsBackground = true
    // Matches the .padding(4) the SwiftUI log had.
    textView.textContainerInset = NSSize(width: 4, height: 4)
    // Wrap long lines instead of scrolling sideways, as the Text rows did.
    textView.textContainer?.widthTracksTextView = true
    textView.isHorizontallyResizable = false
    return scrollView
  }

  /// Performs a planned edit on the log's text view.
  static func apply(
    _ update: ConsoleLogUpdate, to textView: NSTextView,
    attributes: [NSAttributedString.Key: Any]
  ) {
    guard let storage = textView.textStorage else { return }
    switch update {
    case .none:
      break
    case .replace(let text):
      storage.setAttributedString(NSAttributedString(string: text, attributes: attributes))
    case .append(let text):
      storage.append(NSAttributedString(string: text, attributes: attributes))
    case .trimThenAppend(let dropCharacters, let text):
      let selection = textView.selectedRange()
      storage.beginEditing()
      storage.deleteCharacters(in: NSRange(location: 0, length: dropCharacters))
      if !text.isEmpty {
        storage.append(NSAttributedString(string: text, attributes: attributes))
      }
      storage.endEditing()
      textView.setSelectedRange(shiftSelection(selection, by: dropCharacters))
    }
  }

  /// The smallest edit that takes the view from `previous` to `current`.
  static func plan(previous: [String], current: [String]) -> ConsoleLogUpdate {
    if previous == current { return .none }
    if previous.isEmpty || current.isEmpty { return .replace(text(for: current)) }

    // Longest suffix of `previous` that is still a prefix of `current`: what
    // survived on screen. Anything before it was trimmed off the front,
    // anything after it is new output. Quadratic in the worst case but
    // linear in practice — a wrong overlap dies on its first comparison.
    for overlap in stride(from: min(previous.count, current.count), through: 1, by: -1) {
      guard previous.suffix(overlap).elementsEqual(current.prefix(overlap)) else { continue }
      let added = current.dropFirst(overlap)
      let appended = added.isEmpty ? "" : lineSeparator + text(for: added)
      let dropped = previous.count - overlap
      if dropped == 0 { return .append(appended) }
      // The trimmed lines plus the newline that separated them from the
      // first surviving one.
      let head = text(for: previous.prefix(dropped)) + lineSeparator
      return .trimThenAppend(dropCharacters: (head as NSString).length, append: appended)
    }
    // No overlap at all — the log was replaced wholesale.
    return .replace(text(for: current))
  }

  /// The edit from `previous` to `current` when the log's owner says exactly
  /// how many lines its cap dropped off the front in between. No guessing:
  /// `plan(previous:current:)` infers that number from the text, and on a log
  /// of repeated lines (a startup script's identical tracebacks) it can pick a
  /// too-long overlap and drop the new lines. Anything inconsistent with the
  /// count — a log cleared or rewritten — is rebuilt wholesale.
  static func plan(previous: [String], current: [String], dropped: Int) -> ConsoleLogUpdate {
    if dropped == 0 && previous == current { return .none }
    let overlap = previous.count - dropped
    guard dropped >= 0, overlap > 0, current.count >= overlap,
      previous.suffix(overlap).elementsEqual(current.prefix(overlap))
    else { return .replace(text(for: current)) }
    let added = current.dropFirst(overlap)
    let appended = added.isEmpty ? "" : lineSeparator + text(for: added)
    if dropped == 0 { return .append(appended) }
    let head = text(for: previous.prefix(dropped)) + lineSeparator
    return .trimThenAppend(dropCharacters: (head as NSString).length, append: appended)
  }

  /// Slides a selection back over a front-trim, giving up whatever part of it
  /// the trim consumed. Without this, hitting the 400-line cap mid-selection
  /// would silently leave the highlight sitting on different text.
  static func shiftSelection(_ range: NSRange, by dropCharacters: Int) -> NSRange {
    let lost = max(0, min(range.length, dropCharacters - range.location))
    return NSRange(
      location: max(0, range.location - dropCharacters),
      length: range.length - lost)
  }
}

// MARK: - Console hint (modules/mmol/console_hint.py)

struct ConsoleHint: Decodable {
  struct Param: Decodable {
    let name: String
    let req: Bool
    let `default`: String
    var minor: Bool? = nil  // left out of the signature unless being typed
  }
  struct Choice: Decodable { let text, desc, post: String }
  let cmd, about, argdoc, example: String
  let params: [Param]
  let arg: Int  // index into params of the argument being typed, or -1
  let choices: [Choice]
  let start: Int  // where the current word starts (unicode scalars)

  func word(in line: String) -> String {
    String(String.UnicodeScalarView(line.unicodeScalars.dropFirst(start)))
  }

  /// `line` with the current word replaced by `choice` (+ its postfix when final).
  func applying(_ choice: Choice, to line: String, final: Bool) -> String {
    String(String.UnicodeScalarView(line.unicodeScalars.prefix(start))) + choice.text
      + (final && !choice.text.hasSuffix("/") ? choice.post : "")
  }
}

/// The log, observing only ConsoleLog (see PyMOLEngine.console).
private struct ConsoleLogView: View {
  @ObservedObject var log: ConsoleLog
  var body: some View { ConsoleTextView(entries: log.lines, trimmed: log.trimmed) }
}

// MARK: - Command Text Field (handles up/down arrow keys)

struct CommandTextField: NSViewRepresentable {
  @Binding var text: String
  var textColor: Color
  var bgColor: Color
  var fontSize: CGFloat
  var onSubmit: () -> Void
  var onUpArrow: () -> Void
  var onDownArrow: () -> Void
  var onComplete: (String) -> String?

  func makeNSView(context: Context) -> NSTextField {
    // Arrow-key history is handled in the delegate's doCommandBy (moveUp:/
    // moveDown:), not via an NSTextField.keyDown override — a focused field's
    // key events are swallowed by the window's field editor, so keyDown never
    // sees the arrows. A plain NSTextField is therefore sufficient.
    let field = NSTextField()
    field.delegate = context.coordinator
    field.font = .monospacedSystemFont(ofSize: fontSize, weight: .regular)
    field.textColor = NSColor(textColor)
    field.backgroundColor = NSColor(bgColor)
    field.isBordered = false
    field.focusRingType = .none
    field.placeholderString = "Enter command..."
    field.cell?.sendsActionOnEndEditing = false
    // Paint our own background, and keep macOS's autocompletion off a
    // command line: either can flash an empty white box over the field
    // as it takes focus.
    field.drawsBackground = true
    field.isAutomaticTextCompletionEnabled = false
    return field
  }

  func updateNSView(_ nsView: NSTextField, context: Context) {
    if nsView.stringValue != text {
      // While editing, the field editor holds the text: write it there.
      if let editor = nsView.currentEditor() as? NSTextView {
        editor.string = text
        editor.setSelectedRange(NSRange(location: (text as NSString).length, length: 0))
      } else {
        nsView.stringValue = text
      }
    }
    context.coordinator.onSubmit = onSubmit
    context.coordinator.onUpArrow = onUpArrow
    context.coordinator.onDownArrow = onDownArrow
    context.coordinator.parent = self
  }

  func makeCoordinator() -> Coordinator {
    Coordinator(self)
  }

  class Coordinator: NSObject, NSTextFieldDelegate {
    var parent: CommandTextField
    var onSubmit: () -> Void
    var onUpArrow: () -> Void
    var onDownArrow: () -> Void

    init(_ parent: CommandTextField) {
      self.parent = parent
      self.onSubmit = parent.onSubmit
      self.onUpArrow = parent.onUpArrow
      self.onDownArrow = parent.onDownArrow
    }

    func controlTextDidChange(_ notification: Notification) {
      guard let field = notification.object as? NSTextField else { return }
      parent.text = field.stringValue
    }

    func control(
      _ control: NSControl, textView: NSTextView,
      doCommandBy selector: Selector
    ) -> Bool {
      if selector == #selector(NSResponder.insertNewline(_:)) {
        onSubmit()
        return true
      }
      // Tab → completion from the console hint (CommandPanel.tabComplete).
      // Always consume Tab (don't shift keyboard focus out of the field).
      if selector == #selector(NSResponder.insertTab(_:)) {
        let current = textView.string
        if let completed = parent.onComplete(current), completed != current {
          textView.string = completed
          parent.text = completed
          textView.setSelectedRange(NSRange(location: (completed as NSString).length, length: 0))
        }
        return true
      }
      // Up/Down arrows → command history. While the field is focused its key
      // events go to the window's shared field editor (an NSTextView), so a
      // focused NSTextField never sees keyDown for the arrows — they arrive
      // here as moveUp:/moveDown: instead (the same delegate path that makes
      // Return and Tab work). Recall via the history closures, then push the
      // recalled text straight into the field editor (setting the field's
      // stringValue while it is being edited is unreliable — mirror the Tab
      // handling above) and put the caret at the end.
      if selector == #selector(NSResponder.moveUp(_:)) {
        onUpArrow()
        let s = parent.text
        textView.string = s
        textView.setSelectedRange(NSRange(location: (s as NSString).length, length: 0))
        return true
      }
      if selector == #selector(NSResponder.moveDown(_:)) {
        onDownArrow()
        let s = parent.text
        textView.string = s
        textView.setSelectedRange(NSRange(location: (s as NSString).length, length: 0))
        return true
      }
      return false
    }
  }
}

// MARK: - Preview
