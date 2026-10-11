#!/bin/bash
# bundle_biopython.sh — bundle Biopython's `Bio` package into the embedded
# Python site-packages of the macOS app (python-build-standalone 3.13), so
# `import Bio` works in the embedded interpreter that runs Raymond's run_python
# tool. The cp313 arm64 C extensions are kept (they load natively).
#
# Idempotent: re-run to refresh. The deps_macos tree is gitignored.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
PYVER="3.13"
ABI="cp313"
BIO_VERSION="${BIO_VERSION:-1.87}"

MAC_SP="$REPO/deps_macos/python-standalone/python/lib/python${PYVER}/site-packages"

TMP="$(mktemp -d /tmp/raymol_bio.XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

echo "=== Downloading biopython==${BIO_VERSION} (${ABI}, macOS arm64 wheel) ==="
python3 -m pip download "biopython==${BIO_VERSION}" \
  --only-binary=:all: --no-deps \
  --python-version "$PYVER" --platform macosx_11_0_arm64 --abi "$ABI" \
  -d "$TMP" >/dev/null
WHL="$(ls "$TMP"/biopython-*.whl)"
echo "  $WHL"

EX="$TMP/extract"
mkdir -p "$EX"
(cd "$EX" && unzip -q "$WHL")

trim() { # drop __pycache__ + bundled tests to keep the bundle lean
  find "$1" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
  find "$1" -type d -name "tests" -exec rm -rf {} + 2>/dev/null || true
}

# --- full Bio (keep arm64 .so; the app has numpy) ---------------------------
if [ -d "$MAC_SP" ]; then
  echo "=== macOS: install Bio -> $MAC_SP ==="
  rm -rf "$MAC_SP/Bio" "$MAC_SP"/biopython-*.dist-info
  cp -R "$EX/Bio" "$MAC_SP/Bio"
  cp -R "$EX"/biopython-*.dist-info "$MAC_SP/"
  trim "$MAC_SP/Bio"
  du -sh "$MAC_SP/Bio"
else
  echo "  SKIP macOS (run scripts/fetch_macos_python.sh first; $MAC_SP missing)"
fi

echo "=== Done. Verify with: <embedded-python> -c 'import Bio; print(Bio.__version__)' ==="
