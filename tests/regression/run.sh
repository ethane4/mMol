#!/bin/bash
# Runs the mMol regression test against a built app.
#
# Usage: tests/regression/run.sh [--update] [APP]
#
# APP defaults to build_dd/mMol.app. The app starts with an empty home folder,
# so no .pymolrc is read, runs regression.py and is stopped when it finishes.
# The exit status is 0 when every check matches golden/. With --update, the
# reference files in golden/ are rewritten from this run instead.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
update=0
if [[ ${1:-} == --update ]]; then
  update=1
  shift
fi
app=${1:-$here/../../build_dd/mMol.app}
exe="$app/Contents/MacOS/$(basename "$app" .app)"
work=$(mktemp -d -t mmol-test)
mkdir "$work/home"
mkdir -p "$here/golden"

HOME="$work/home" MMOL_TEST_ROOT="$here" MMOL_TEST_OUT="$work" \
  MMOL_TEST_UPDATE=$update PYMOL_AUTOCMD_LATE="run $here/regression.py" \
  PYMOL_AUTOCMD_LATE_DELAY=3 "$exe" >"$work/app.log" 2>&1 &
pid=$!
for _ in $(seq 300); do
  if grep -q '^DONE' "$work/report.txt" 2>/dev/null; then
    break
  fi
  if ! kill -0 "$pid" 2>/dev/null; then
    echo "mMol exited early; see $work/app.log" >&2
    exit 2
  fi
  sleep 1
done
kill "$pid" 2>/dev/null || true
wait "$pid" 2>/dev/null || true

if [[ ! -f $work/report.txt ]]; then
  echo "no report after 300 s; see $work/app.log" >&2
  exit 2
fi
if ((update)); then
  echo "golden files updated from $work"
  exit 0
fi
cat "$work/report.txt"
echo "output: $work"
grep -q '^DONE 0 failure' "$work/report.txt"
