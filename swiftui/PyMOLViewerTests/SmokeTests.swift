import XCTest

@testable import mMol  // module name = PRODUCT_NAME "mMol"; try `import PyMOLViewer` if this fails

final class SmokeTests: XCTestCase {
  func testHarnessRuns() { XCTAssertEqual(1 + 1, 2) }
}
