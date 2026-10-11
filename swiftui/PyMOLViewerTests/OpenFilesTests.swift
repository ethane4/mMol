import XCTest

@testable import mMol

// Regression tests for issue #222: `open a.pdb b.pdb c.pdb` from a Terminal (and
// Finder "Open With" multi-select / Dock-icon drop) only loaded the first file.
// macOS delivers all of them via NSApplicationDelegate.application(_:open:) as one
// [URL] array; MMolAppDelegate forwards that array through handleOpenedURLs.
// These lock in that EVERY url is forwarded — a spy stands in for the real
// loadOpenedFile so the engine/PyMOL core never has to boot.
@MainActor
final class OpenFilesTests: XCTestCase {

  func testAllOpenedURLsAreForwardedInOrder() {
    let urls = ["a.pdb", "b.pdb", "c.pdb"].map { URL(fileURLWithPath: "/tmp/\($0)") }
    var seen: [String] = []
    handleOpenedURLs(urls, into: PyMOLEngine.shared) { url, _ in
      seen.append(url.lastPathComponent)
    }
    // The bug loaded only "a.pdb"; the fix must forward all three, in order.
    XCTAssertEqual(seen, ["a.pdb", "b.pdb", "c.pdb"])
  }

  func testSingleOpenedURLStillForwarded() {
    // A plain Finder double-click delivers a one-element array; it must still
    // load (the delegate now handles single- and multi-file opens uniformly).
    let urls = [URL(fileURLWithPath: "/tmp/only.pse")]
    var seen: [String] = []
    handleOpenedURLs(urls, into: PyMOLEngine.shared) { url, _ in
      seen.append(url.lastPathComponent)
    }
    XCTAssertEqual(seen, ["only.pse"])
  }

  func testEmptyURLListIsNoOp() {
    var called = false
    handleOpenedURLs([], into: PyMOLEngine.shared) { _, _ in called = true }
    XCTAssertFalse(called)
  }

  // MARK: - Issue #272: sessions load bare — no object name

  func testSessionLoadHasNoObjectName() {
    // `load foo.pse, name` creates nothing named `name` (the session restores
    // its own objects).
    XCTAssertEqual(
      PyMOLEngine.loadInvocation(
        path: "/tmp/open_E277B7EE.pse",
        name: "top6_candidates"),
      "load /tmp/open_E277B7EE.pse")
  }

  func testSessionExtensionMatchIsCaseInsensitiveAndCoversPsw() {
    for path in ["/tmp/a.PSE", "/tmp/b.psw", "/tmp/c.Psw"] {
      XCTAssertEqual(PyMOLEngine.loadInvocation(path: path, name: "x"), "load \(path)", path)
    }
  }

  func testCoordinateLoadKeepsObjectName() {
    XCTAssertEqual(
      PyMOLEngine.loadInvocation(path: "/tmp/5hbh.pdb", name: "5hbh"),
      "load /tmp/5hbh.pdb, 5hbh")
  }

  // MARK: - Issue #349: a .pse open must not silently wipe a non-empty session

  func testSessionOpenOverNonEmptySessionNeedsConfirmation() {
    let pse = URL(fileURLWithPath: "/tmp/other.pse")
    XCTAssertTrue(openWouldReplaceSession(pse, hasObjects: true))
  }

  func testSessionOpenIntoEmptySessionNeedsNoConfirmation() {
    // Cold launch from a Finder double-click: nothing to lose, no prompt.
    let pse = URL(fileURLWithPath: "/tmp/other.pse")
    XCTAssertFalse(openWouldReplaceSession(pse, hasObjects: false))
  }

  func testCoordinateOpenNeverPrompts() {
    // A .pdb/.cif ADDS an object — it doesn't replace the session.
    let pdb = URL(fileURLWithPath: "/tmp/5hbh.pdb")
    XCTAssertFalse(openWouldReplaceSession(pdb, hasObjects: true))
  }
}
