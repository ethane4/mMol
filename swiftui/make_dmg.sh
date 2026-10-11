#!/bin/bash
# make_dmg.sh — package mMol.app as build_dmg/mMol-<version>.dmg.
#
#   swiftui/make_dmg.sh     ad-hoc signed: on first launch users right-click
#                           mMol > Open (or allow it in System Settings > Privacy)
#
#   NOTARIZE=1 DEVID="Developer ID Application: NAME (TEAMID)" \
#   NOTARY_PROFILE=mmol-notary swiftui/make_dmg.sh
#                           Developer ID signed + notarized + stapled: opens with
#                           no warning. One-time setup: a Developer ID Application
#                           certificate in the keychain, and
#                           xcrun notarytool store-credentials mmol-notary \
#                             --apple-id <id> --team-id <TEAMID> --password <app-specific>
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SW="$ROOT/swiftui"
OUT="$ROOT/build_dmg"

# A side-copy name, so a running build_dd/mMol.app is never replaced; the DMG
# takes Xcode's own product, which is always named mMol.app.
CONFIG=Release NAME=mMol-dmg bash "$SW/build_mac.sh"
APP="$ROOT/build_dd/Build/Products/Release/mMol.app"
VERSION=$(plutil -extract CFBundleShortVersionString raw "$APP/Contents/Info.plist")

rm -rf "$OUT"
mkdir -p "$OUT/root"
cp -R "$APP" "$OUT/root/"
APP="$OUT/root/mMol.app"

if [[ "${NOTARIZE:-0}" == 1 ]]; then
  : "${DEVID:?Set DEVID to your 'Developer ID Application: NAME (TEAMID)' identity}"
  PROFILE="${NOTARY_PROFILE:-mmol-notary}"
  # every nested Mach-O first (deepest paths first), then the bundle
  find "$APP/Contents" -type f \( -name '*.dylib' -o -name '*.so' -o -perm -u+x \) -print0 |
    xargs -0 file | grep 'Mach-O' | cut -d: -f1 | awk '{ print length, $0 }' | sort -rn | cut -d' ' -f2- |
    while read -r f; do codesign --force --options runtime --timestamp --sign "$DEVID" "$f"; done
  codesign --force --options runtime --timestamp \
    --entitlements "$SW/mMol_DeveloperID.entitlements" --sign "$DEVID" "$APP"
  codesign --verify --deep --strict "$APP"
  ditto -c -k --keepParent "$APP" "$OUT/mMol.zip"
  xcrun notarytool submit "$OUT/mMol.zip" --keychain-profile "$PROFILE" --wait
  xcrun stapler staple "$APP"
  rm "$OUT/mMol.zip"
fi

ln -s /Applications "$OUT/root/Applications"
DMG="$OUT/mMol-$VERSION.dmg"
hdiutil create -volname mMol -srcfolder "$OUT/root" -ov -format UDZO "$DMG" >/dev/null
rm -rf "$OUT/root"

if [[ "${NOTARIZE:-0}" == 1 ]]; then
  codesign --force --timestamp --sign "$DEVID" "$DMG"
  xcrun notarytool submit "$DMG" --keychain-profile "$PROFILE" --wait
  xcrun stapler staple "$DMG"
fi
echo "Built $DMG"
