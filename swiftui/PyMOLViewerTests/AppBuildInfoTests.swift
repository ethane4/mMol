import XCTest

@testable import mMol

/// Covers the pure half of `AppBuildInfo` — the part that decides what a user
/// sees when asked "which build is this?".
final class AppBuildInfoTests: XCTestCase {

  // MARK: - displayVersion

  func testUsesVersionAndBuild() {
    XCTAssertEqual(
      AppBuildInfo.displayVersion(version: "1.9.1", build: "27"),
      "1.9.1 (27)")
  }

  // Incomplete Info.plist: degrade to whichever half exists rather than
  // rendering "1.9.1 ()" or " (27)" at the user.
  func testMissingHalvesDegradeCleanly() {
    XCTAssertEqual(
      AppBuildInfo.displayVersion(version: "1.9.1", build: ""),
      "1.9.1")
    XCTAssertEqual(
      AppBuildInfo.displayVersion(version: "", build: "27"),
      "(27)")
    XCTAssertEqual(
      AppBuildInfo.displayVersion(version: "", build: ""),
      "")
  }

  // Whatever the bundle says, the footer must have something to show.
  func testBundleDisplayVersionIsNonEmpty() {
    XCTAssertFalse(AppBuildInfo.displayVersion.isEmpty)
  }

  // MARK: - consoleBanner (#377)

  // The console's first line must come from the bundle, not a hand-bumped
  // string: the old marker read "v35" for two months of unrelated releases.
  func testConsoleBannerCarriesTheDisplayVersion() {
    XCTAssertEqual(
      AppBuildInfo.consoleBanner(displayVersion: "1.10.0 (27)"),
      "  mMol 1.10.0 (27) · built on open-source PyMOL (Schrödinger, LLC)")
  }

  // A plist missing both keys must still print something a tester can quote,
  // not a trailing-space "mMol " that looks like a rendering glitch.
  func testConsoleBannerWithoutVersionSaysSo() {
    XCTAssertEqual(
      AppBuildInfo.consoleBanner(displayVersion: ""),
      "  mMol (unknown version) · built on open-source PyMOL (Schrödinger, LLC)")
  }

  // The live property must go through the same formatter as the tested one.
  func testLiveConsoleBannerMatchesDisplayVersion() {
    XCTAssertEqual(
      AppBuildInfo.consoleBanner,
      AppBuildInfo.consoleBanner(displayVersion: AppBuildInfo.displayVersion))
    XCTAssertTrue(AppBuildInfo.consoleBanner.hasPrefix("  mMol "))
  }
}
