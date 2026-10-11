import Foundation

/// The app's own version identity, as a user should read it: the version and
/// the build number together, since the version alone does not say which build
/// is running.
enum AppBuildInfo {

  /// `CFBundleShortVersionString`, e.g. "1.9.1". Empty if somehow absent.
  static var version: String { string(for: "CFBundleShortVersionString") }

  /// `CFBundleVersion`, e.g. "27". Empty if somehow absent.
  static var build: String { string(for: "CFBundleVersion") }

  /// What to show the user, e.g. "1.9.1 (27)".
  static var displayVersion: String {
    Self.displayVersion(version: version, build: build)
  }

  /// The line the in-app console prints first at launch, e.g.
  /// " [build] mMol 1.10.0 (27)". Derived from the bundle so it can never go
  /// stale: it replaced a hand-bumped "v35" marker that sat two months behind
  /// the real version (GitHub issue #377).
  static var consoleBanner: String { consoleBanner(displayVersion: displayVersion) }

  // MARK: - Pure logic (unit-tested; no Bundle access)

  /// The "version (build)" convention; degrades to whichever half exists if
  /// the Info.plist is incomplete rather than rendering stray parentheses.
  static func displayVersion(version: String, build: String) -> String {
    switch (version.isEmpty, build.isEmpty) {
    case (false, false): return "\(version) (\(build))"
    case (false, true): return version
    case (true, false): return "(\(build))"
    case (true, true): return ""
    }
  }

  /// Leading space matches the console's other status lines. An Info.plist
  /// with neither key still prints an honest line rather than "mMol ".
  static func consoleBanner(displayVersion: String) -> String {
    (displayVersion.isEmpty ? "  mMol (unknown version)" : "  mMol \(displayVersion)")
      + " · built on open-source PyMOL (Schrödinger, LLC)"
  }

  /// The startup logo at the top of the console.
  static let logo = [
    "               ██▄   ▄██           ██",
    "   ▄▄▄▄▄▄▄▄▄   ██▀█▄█▀██   ▄▄▄▄▄   ██",
    "   ██  ██  ██  ██  ▀  ██  ██   ██  ██",
    "   ██  ██  ██  ██     ██  ▀█▄▄▄█▀  ▀█▄▄",
  ]

  // MARK: - Bundle access

  private static func string(for key: String) -> String {
    (Bundle.main.object(forInfoDictionaryKey: key) as? String) ?? ""
  }
}
