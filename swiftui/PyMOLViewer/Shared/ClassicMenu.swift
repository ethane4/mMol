// ClassicMenu.swift — desktop-PyMOL-style popup menus for the Classic PyMOL
// Layout's A/S/H/L/C buttons.
//
// The system menu (NSMenu, which SwiftUI's Menu becomes) is drawn translucent with
// rounded corners and cannot be restyled. PyMOL 2.5's menus are opaque, square and
// flat, so the classic layout presents its own: borderless panels holding a
// SwiftUI item list, a submenu panel opening to the right of its row on hover,
// everything closing on a choice, Escape, a click elsewhere or the app resigning.
// The items are the same commands the mMol menus run (see the `classicItems`
// builders in ObjectPanel.swift); only the presentation differs.

import AppKit
import SwiftUI

/// One row of a classic menu.
struct ClassicMenuItem: Identifiable {
  enum Kind {
    case action(() -> Void)
    /// A submenu; `primary`, when set, also makes the row itself clickable
    /// (the colour menu's "click the base colour, hover for its shades").
    case submenu([ClassicMenuItem], primary: (() -> Void)?)
    case separator
    case disabled
  }

  let id = UUID()
  var title: String
  var kind: Kind
  var swatch: Color? = nil
  /// Red text: removes or deletes something.
  var destructive: Bool = false

  static func action(
    _ title: String, swatch: Color? = nil, destructive: Bool? = nil,
    _ run: @escaping () -> Void
  ) -> ClassicMenuItem {
    ClassicMenuItem(
      title: title, kind: .action(run), swatch: swatch,
      destructive: destructive ?? isRemoval(title))
  }
  static func submenu(
    _ title: String, swatch: Color? = nil, primary: (() -> Void)? = nil,
    _ children: [ClassicMenuItem]
  ) -> ClassicMenuItem {
    ClassicMenuItem(title: title, kind: .submenu(children, primary: primary), swatch: swatch)
  }
  static var separator: ClassicMenuItem { ClassicMenuItem(title: "", kind: .separator) }
  static func disabled(_ title: String) -> ClassicMenuItem {
    ClassicMenuItem(title: title, kind: .disabled)
  }

  /// "Remove Waters", "delete", "Delete everything", Hydrogens ▸ remove …
  static func isRemoval(_ title: String) -> Bool {
    let t = title.lowercased()
    return t.hasPrefix("remove") || t.hasPrefix("delete")
  }

  var height: CGFloat {
    if case .separator = kind { return ClassicMenuStyle.separatorHeight }
    return ClassicMenuStyle.rowHeight
  }
}

enum ClassicMenuStyle {
  static let rowHeight: CGFloat = 20
  static let separatorHeight: CGFloat = 7
  static let background = Color(white: 0.18)
  static let border = NSColor(white: 0.42, alpha: 1)
  static let text = Color(white: 0.95)
  static let disabledText = Color(white: 0.5)
  static let destructiveText = Color(.sRGB, red: 1.0, green: 0.36, blue: 0.36, opacity: 1)
  static let highlight = Color(.sRGB, red: 0.22, green: 0.42, blue: 0.85, opacity: 1)
  static let font = Font.system(size: 13)
  static let minWidth: CGFloat = 130
}

/// Owns the open menu panels (one per open level) and the monitors that close them.
@MainActor
final class ClassicMenuController {
  static let shared = ClassicMenuController()

  private var panels: [NSPanel] = []
  private var monitors: [Any] = []
  /// The button the open menu belongs to: clicking it again closes the menu
  /// (rather than the outside-click monitor closing it and the button's own
  /// mouseDown reopening it at once).
  private weak var opener: NSView?

  func isOpen(from view: NSView) -> Bool { !panels.isEmpty && opener === view }
  private var resignObserver: NSObjectProtocol?

  /// Open `items` with the menu's top-left corner at `topLeft` (screen coords).
  func show(_ items: [ClassicMenuItem], topLeft: NSPoint, from view: NSView? = nil) {
    close()
    opener = view
    open(items, level: 0, topLeft: topLeft)
    installMonitors()
  }

  func close() {
    pendingHover?.cancel()
    opener = nil
    panels.forEach { $0.orderOut(nil) }
    panels.removeAll()
    monitors.forEach { NSEvent.removeMonitor($0) }
    monitors.removeAll()
    if let o = resignObserver { NotificationCenter.default.removeObserver(o) }
    resignObserver = nil
  }

  private var pendingHover: DispatchWorkItem?
  /// Where the pointer entered the previous row: the apex of the "safe triangle".
  private var lastHoverPoint: NSPoint?

  /// Hovering a row at `level`. Switching is instant, except while the pointer is
  /// heading for an open submenu (inside the triangle from where it entered the
  /// last row to the submenu's near edge): then crossing sibling rows on the way
  /// waits briefly, and entering the submenu cancels the wait.
  func hover(level: Int, _ apply: @escaping () -> Void) {
    pendingHover?.cancel()
    let point = NSEvent.mouseLocation
    defer { lastHoverPoint = point }
    guard panels.count > level + 1, let apex = lastHoverPoint,
      headingInto(panels[level + 1].frame, from: apex, to: point)
    else { return apply() }
    let work = DispatchWorkItem(block: apply)
    pendingHover = work
    DispatchQueue.main.asyncAfter(deadline: .now() + 0.3, execute: work)
  }

  /// Whether `p` lies in the triangle from `apex` to the submenu's near edge.
  private func headingInto(_ sub: NSRect, from apex: NSPoint, to p: NSPoint) -> Bool {
    let x = apex.x <= sub.minX ? sub.minX : sub.maxX
    let (b, c) = (NSPoint(x: x, y: sub.maxY), NSPoint(x: x, y: sub.minY))
    func side(_ u: NSPoint, _ v: NSPoint) -> CGFloat {
      (v.x - u.x) * (p.y - u.y) - (v.y - u.y) * (p.x - u.x)
    }
    let (d1, d2, d3) = (side(apex, b), side(b, c), side(c, apex))
    return !((d1 < 0 || d2 < 0 || d3 < 0) && (d1 > 0 || d2 > 0 || d3 > 0))
  }

  /// Back on the row whose submenu is open: drop any pending switch.
  func cancelPendingHover() { pendingHover?.cancel() }

  /// Hovering a submenu row: open its children beside it, replacing deeper levels.
  func openSubmenu(_ items: [ClassicMenuItem], level: Int, rowTop: CGFloat) {
    closeLevels(from: level)
    guard level - 1 < panels.count else { return }
    let parent = panels[level - 1].frame
    open(items, level: level, topLeft: NSPoint(x: parent.maxX - 1, y: parent.maxY - rowTop))
  }

  /// Hovering a plain row: anything opened from a sibling goes away.
  func closeLevels(from level: Int) {
    while panels.count > level { panels.removeLast().orderOut(nil) }
  }

  func choose(_ run: @escaping () -> Void) {
    close()
    // After the panels are gone, so a command that opens a sheet or a panel
    // of its own is not fighting a menu still on screen.
    DispatchQueue.main.async(execute: run)
  }

  private func open(_ items: [ClassicMenuItem], level: Int, topLeft: NSPoint) {
    let host = NSHostingView(rootView: ClassicMenuList(items: items, level: level))
    let size = host.fittingSize
    let panel = NSPanel(
      contentRect: NSRect(origin: .zero, size: size),
      styleMask: [.borderless, .nonactivatingPanel],
      backing: .buffered, defer: false)
    panel.isOpaque = true
    panel.backgroundColor = NSColor(ClassicMenuStyle.background)
    panel.hasShadow = true
    panel.level = .popUpMenu
    panel.isReleasedWhenClosed = false
    panel.contentView = host
    host.wantsLayer = true
    host.layer?.borderWidth = 1
    host.layer?.borderColor = ClassicMenuStyle.border.cgColor
    host.layer?.cornerRadius = 0

    // Keep it on the screen: flip a submenu to the left of its parent, or
    // nudge a root menu, when it would run off the edge.
    var origin = NSPoint(x: topLeft.x, y: topLeft.y - size.height)
    if let screen = NSScreen.screens.first(where: { $0.frame.contains(topLeft) }) ?? NSScreen.main {
      let vis = screen.visibleFrame
      if origin.x + size.width > vis.maxX {
        origin.x =
          level > 0 && level - 1 < panels.count
          ? panels[level - 1].frame.minX - size.width + 1
          : vis.maxX - size.width
      }
      origin.y = max(vis.minY, min(origin.y, vis.maxY - size.height))
    }
    panel.setFrameOrigin(origin)
    panel.orderFrontRegardless()
    panels.append(panel)
  }

  private func installMonitors() {
    // A click in any other window of the app closes the menu and still counts.
    if let m = NSEvent.addLocalMonitorForEvents(
      matching: [.leftMouseDown, .rightMouseDown, .otherMouseDown],
      handler: { [weak self] event in
        guard let self else { return event }
        if let w = event.window, self.panels.contains(where: { $0 === w }) { return event }
        // A click on the menu's own button is left to the button, which
        // toggles the menu closed.
        if let o = self.opener, let w = o.window, event.window === w,
          o.bounds.contains(o.convert(event.locationInWindow, from: nil))
        {
          return event
        }
        self.close()
        return event
      })
    {
      monitors.append(m)
    }
    // Escape closes; every other key is left alone.
    if let m = NSEvent.addLocalMonitorForEvents(
      matching: .keyDown,
      handler: { [weak self] event in
        if event.keyCode == 53 {
          self?.close()
          return nil
        }
        return event
      })
    {
      monitors.append(m)
    }
    // Clicks in other apps, and the app losing focus.
    if let m = NSEvent.addGlobalMonitorForEvents(
      matching: [.leftMouseDown, .rightMouseDown],
      handler: { [weak self] _ in self?.close() })
    {
      monitors.append(m)
    }
    resignObserver = NotificationCenter.default.addObserver(
      forName: NSApplication.didResignActiveNotification, object: nil, queue: .main
    ) { [weak self] _ in
      MainActor.assumeIsolated { self?.close() }
    }
  }
}

/// One level of a classic menu.
private struct ClassicMenuList: View {
  let items: [ClassicMenuItem]
  let level: Int
  @State private var hovered: UUID?

  var body: some View {
    VStack(alignment: .leading, spacing: 0) {
      ForEach(Array(items.enumerated()), id: \.element.id) { index, item in
        row(item, top: rowTop(index))
      }
    }
    .padding(.vertical, 3)
    .frame(minWidth: ClassicMenuStyle.minWidth, alignment: .leading)
    .fixedSize()
    .background(ClassicMenuStyle.background)
  }

  /// Distance from the panel's top edge to row `index` (the 3 pt top padding
  /// plus the rows above), so a submenu can open level with its row.
  private func rowTop(_ index: Int) -> CGFloat {
    3 + items.prefix(index).reduce(0) { $0 + $1.height }
  }

  @ViewBuilder
  private func row(_ item: ClassicMenuItem, top: CGFloat) -> some View {
    switch item.kind {
    case .separator:
      Rectangle().fill(Color(white: 0.38)).frame(height: 1)
        .padding(.vertical, 3)
        .frame(maxWidth: .infinity)
    default:
      let isOn = hovered == item.id
      let enabled: Bool = {
        if case .disabled = item.kind { return false }
        return true
      }()
      HStack(spacing: 6) {
        if let swatch = item.swatch {
          Rectangle().fill(swatch).frame(width: 10, height: 10)
            .overlay(Rectangle().strokeBorder(Color(white: 0.55), lineWidth: 0.5))
        } else if items.contains(where: { $0.swatch != nil }) {
          Color.clear.frame(width: 10, height: 10)  // keep labels aligned
        }
        Text(item.title)
          .font(ClassicMenuStyle.font)
          .foregroundColor(
            !enabled
              ? ClassicMenuStyle.disabledText
              : isOn
                ? .white
                : item.destructive
                  ? ClassicMenuStyle.destructiveText
                  : ClassicMenuStyle.text
          )
          .lineLimit(1)
          .fixedSize()
        Spacer(minLength: 16)
        if case .submenu = item.kind {
          Text("▸").font(.system(size: 10))
            .foregroundColor(isOn ? .white : ClassicMenuStyle.text)
        }
      }
      .padding(.horizontal, 10)
      .frame(height: ClassicMenuStyle.rowHeight)
      .frame(maxWidth: .infinity, alignment: .leading)
      .background(isOn && enabled ? ClassicMenuStyle.highlight : Color.clear)
      .contentShape(Rectangle())
      .onHover { inside in
        guard inside, enabled else { return }
        let menu = ClassicMenuController.shared
        if hovered == item.id { return menu.cancelPendingHover() }
        menu.hover(level: level) {
          hovered = item.id
          if case .submenu(let children, _) = item.kind {
            menu.openSubmenu(children, level: level + 1, rowTop: top)
          } else {
            menu.closeLevels(from: level + 1)
          }
        }
      }
      .onTapGesture {
        switch item.kind {
        case .action(let run): ClassicMenuController.shared.choose(run)
        case .submenu(_, let primary?): ClassicMenuController.shared.choose(primary)
        default: break
        }
      }
    }
  }
}

/// An A/S/H/L/C face that opens a classic menu below itself on mouse-down,
/// as the internal GUI's buttons do.
struct ClassicMenuTrigger<Face: View>: View {
  let face: Face
  let items: () -> [ClassicMenuItem]

  init(@ViewBuilder face: () -> Face, items: @escaping () -> [ClassicMenuItem]) {
    self.face = face()
    self.items = items
  }

  var body: some View {
    face.overlay(ClassicMenuHitArea(items: items))
      // The menu opens from the NSView's mouseDown; this claims the same
      // click for SwiftUI so it doesn't also reach the row's tap gesture
      // (which toggles the object on and off).
      .contentShape(Rectangle())
      .onTapGesture {}
  }
}

private struct ClassicMenuHitArea: NSViewRepresentable {
  let items: () -> [ClassicMenuItem]

  func makeNSView(context: Context) -> HitView { HitView() }
  func updateNSView(_ view: HitView, context: Context) { view.items = items }

  final class HitView: NSView {
    var items: () -> [ClassicMenuItem] = { [] }
    override func acceptsFirstMouse(for event: NSEvent?) -> Bool { true }
    override func mouseDown(with event: NSEvent) {
      // Second click on the button while its menu is open: close it.
      let open = MainActor.assumeIsolated { ClassicMenuController.shared.isOpen(from: self) }
      if open {
        MainActor.assumeIsolated { ClassicMenuController.shared.close() }
      } else {
        openMenu()
      }
    }

    func openMenu() {
      guard let window else { return }
      let r = window.convertToScreen(convert(bounds, to: nil))
      MainActor.assumeIsolated {
        ClassicMenuController.shared.show(
          items(), topLeft: NSPoint(x: r.minX, y: r.minY),
          from: self)
      }
    }
  }
}

extension ClassicMenuItem {
  /// A choice row, check-marked when `on`.
  static func check(_ title: String, _ on: Bool, _ run: @escaping () -> Void) -> ClassicMenuItem {
    .action((on ? "✓ " : "   ") + title, destructive: false, run)
  }
}

extension ClassicMenuController {
  /// Open `items` with the menu's top-left corner at the pointer.
  func showAtPointer(_ items: [ClassicMenuItem]) {
    show(items, topLeft: NSEvent.mouseLocation)
  }
}

extension View {
  /// Right-click (or Control-click) opens `items` as a classic menu at the
  /// pointer, in place of the system context menu.
  func classicContextMenu(_ items: @escaping () -> [ClassicMenuItem]) -> some View {
    overlay(ClassicRightClickArea(items: items))
  }
}

private struct ClassicRightClickArea: NSViewRepresentable {
  let items: () -> [ClassicMenuItem]

  func makeNSView(context: Context) -> HitView { HitView() }
  func updateNSView(_ view: HitView, context: Context) { view.items = items }

  final class HitView: NSView {
    var items: () -> [ClassicMenuItem] = { [] }

    // Only right-clicks and Control-clicks land here; every other event
    // (clicks, drags, scrolling, hover) passes through to the view beneath.
    override func hitTest(_ point: NSPoint) -> NSView? {
      guard let e = NSApp.currentEvent else { return nil }
      let menuClick =
        e.type == .rightMouseDown
        || (e.type == .leftMouseDown && e.modifierFlags.contains(.control))
      return menuClick ? super.hitTest(point) : nil
    }
    override func rightMouseDown(with event: NSEvent) { open(event) }
    override func mouseDown(with event: NSEvent) { open(event) }

    // At the click, which is where the pointer is.
    private func open(_ event: NSEvent) {
      guard let window else { return }
      let at = window.convertPoint(toScreen: event.locationInWindow)
      let items = items()
      MainActor.assumeIsolated { ClassicMenuController.shared.show(items, topLeft: at) }
    }
  }
}
