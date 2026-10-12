import XCTest

@testable import mMol

/// pyLiteral, Throttle and the feedback routing table: small engine helpers
/// testable without a running core.
final class PyLiteralTests: XCTestCase {

  /// The decoded payload of a pyLiteral expression, or nil if it is not one.
  private func decode(_ expr: String) -> String? {
    let head = "__import__('base64').b64decode('"
    let tail = "').decode()"
    guard expr.hasPrefix(head), expr.hasSuffix(tail) else { return nil }
    let b64 = String(expr.dropFirst(head.count).dropLast(tail.count))
    return Data(base64Encoded: b64).flatMap { String(data: $0, encoding: .utf8) }
  }

  func testRoundTripsAwkwardStrings() {
    for s in [
      "", "plain", "it's", "a \"quoted\" path", "x'''); import os; ('''",
      "back\\slash\nnewline", "Ünïcødé ✓ 分子",
    ] {
      let expr = pyLiteral(s)
      XCTAssertEqual(decode(expr), s)
      // Only base64's alphabet sits inside the quotes, so nothing can end the
      // literal early.
      XCTAssertEqual(expr.filter { $0 == "'" }.count, 4, expr)
    }
  }

  func testThrottleFiresLeadingEdgeThenDefers() {
    let t = Throttle(10)
    var fired = 0
    t.run { fired += 1 }
    XCTAssertEqual(fired, 1, "the first call fires at once")
    t.run { fired += 1 }
    XCTAssertEqual(fired, 1, "a call inside the interval waits for the trailing edge")
    t.cancel(resetClock: true)
    t.run { fired += 1 }
    XCTAssertEqual(fired, 2, "a reset clock fires on the leading edge again")
    t.cancel()
  }

  /// pollFeedback takes the first matching prefix, so no prefix may shadow
  /// another, and a handled prefix must not also be swallowed.
  func testFeedbackPrefixesAreUnambiguous() {
    let handled = PyMOLEngine.feedbackHandlers.map(\.prefix)
    let all = handled + PyMOLEngine.ignoredFeedbackPrefixes
    XCTAssertEqual(Set(all).count, all.count)
    for a in all {
      for b in all where a != b {
        XCTAssertFalse(b.hasPrefix(a), "\(a) shadows \(b)")
      }
    }
  }
}
