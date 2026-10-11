#!/bin/bash
# build_mac.sh — build mMol.app (the Xcode project builds the PyMOL core too).
#
#   swiftui/build_mac.sh                 optimized Release build -> build_dd/mMol.app
#   CONFIG=Debug swiftui/build_mac.sh    unoptimized, for debugging
#   NAME=mMol-test swiftui/build_mac.sh  a side copy (e.g. while mMol is open)
#   INSTALL=1 swiftui/build_mac.sh       also copy to /Applications
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DD="$ROOT/build_dd"
CONFIG="${CONFIG:-Release}"
NAME="${NAME:-mMol}"
export DEVELOPER_DIR="${DEVELOPER_DIR:-/Applications/Xcode.app/Contents/Developer}"

test -d "$ROOT/deps_macos/python-standalone/python" || {
  echo "ERROR: run scripts/fetch_macos_python.sh first"
  exit 1
}
(cd "$ROOT/swiftui" && xcodegen generate --quiet)

xcodebuild -project "$ROOT/swiftui/PyMOLViewer.xcodeproj" -scheme PyMOLViewer_macOS \
  -configuration "$CONFIG" -destination 'platform=macOS' -derivedDataPath "$DD" \
  -skipPackagePluginValidation -skipMacroValidation build

DST="$DD/$NAME.app"
rm -rf "$DST"
cp -R "$DD/Build/Products/$CONFIG/mMol.app" "$DST"
if [[ "$NAME" != mMol ]]; then # rename the side copy, then re-seal its signature
  PLIST="$DST/Contents/Info.plist"
  mv "$DST/Contents/MacOS/mMol" "$DST/Contents/MacOS/$NAME"
  for k in CFBundleExecutable CFBundleName CFBundleDisplayName; do
    plutil -replace "$k" -string "$NAME" "$PLIST"
  done
  codesign --force --deep --sign - "$DST"
fi
codesign -v "$DST"

if [[ "${INSTALL:-0}" == 1 ]]; then
  rm -rf "/Applications/$NAME.app"
  cp -R "$DST" "/Applications/$NAME.app"
  echo "Installed /Applications/$NAME.app"
fi
echo "Built $DST"
