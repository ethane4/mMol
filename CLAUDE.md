# CLAUDE.md

This file guides Claude Code (claude.ai/code) when it works on this repository.
mMol is a macOS molecular viewer (SwiftUI and Metal) built on PyMOL: a C++17
core organized by component, a Python API on top, and a native app around both.

## Build and ship

The Xcode project, generated from `swiftui/project.yml` by xcodegen, builds
everything: the PyMOL core (the `pymol_core` and `champ` static-library
targets), then the app.

One-time setup:

```shell
brew install xcodegen glm libpng freetype libomp
scripts/setup_macos_deps.sh
```

`scripts/setup_macos_deps.sh` installs the embedded Python 3.13 with NumPy and
Biopython into `deps_macos/`.

*   **Build:** `swiftui/build_mac.sh` makes an optimized Release build at
    `build_dd/mMol.app`; use `CONFIG=Debug` for debugging. While an mMol is
    running, build a side copy with `NAME=mMol-test` instead of overwriting it.
*   **Ship:** `swiftui/make_dmg.sh` makes `build_dmg/mMol-<version>.dmg`,
    ad-hoc signed. Use `NOTARIZE=1 DEVID=...` for a notarized Developer ID
    DMG.

## Run the tests

The app's Swift unit tests live in `swiftui/PyMOLViewerTests/` and run through
the `UnitTests_macOS` scheme:

```shell
DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer xcodebuild test \
    -project swiftui/PyMOLViewer.xcodeproj -scheme UnitTests_macOS \
    -destination 'platform=macOS' -skipPackagePluginValidation
```

The regression test in `tests/regression/` checks the whole app: it opens a
built mMol with an empty home folder, ray traces four scenes on the GPU and the
CPU, records every setting, and runs core commands (symmetry, alignment,
sessions, file formats, exports, maps) and the `mmol` modules, then compares
everything with `tests/regression/golden/`:

```shell
tests/regression/run.sh build_dd/mMol-test.app
```

Run it before and after every change to the C++ core. A refactor must pass
unchanged. When a change is meant to alter the output, check the new renders
and values, then rewrite the references with `run.sh --update` and commit them
with the change.

## Code formatting

Google style throughout:

*   C, C++ and Objective-C++: `.clang-format`; format every file you change.
*   Python: yapf (`.style.yapf`) with numpydoc docstrings.
*   Swift: `.swift-format`.
*   Shell: `shfmt -i 2 -ci`.
*   Markdown: the Google Markdown style guide.

## Architecture

### C++ core

The core lives in `src/`, one folder per component under `src/pymol/`.
Headers are included by their path from `src/`, for example
`#include "pymol/scene/scene.h"`, and files are snake_case `.cc` and `.h`.

*   `base/`: globals (`PyMOLGlobals`), feedback, errors, strings, utilities.
*   `platform/`: OS, Python and stubbed OpenGL shims, memory debugging.
*   `geometry/`, `fields/`: vectors, matrices, spatial maps, triangulation,
    sphere meshes; scalar fields, isosurfaces and contouring.
*   `graphics/`: the renderer interface, CGOs, GPU buffers; `metal/` holds
    `RendererMetal`, and `gl/` the OpenGL stubs.
*   `settings/`, `color/`: the settings system (built-in defaults live in
    `settings/setting_info.h`); named colors, ramps and color tables.
*   `scene/`, `ray/`, `ui/`: the scene, camera, movies and picking; the CPU
    and Metal ray tracers; PyMOL's internal GUI pieces (ortho, wizards).
*   `objects/`, `reps/`: object types, coordinate sets and atom info; the
    representations (cartoon, surface, sticks and others).
*   `select/`, `executive/`: the selection language; commands and the object
    list.
*   `io/`, `symmetry/`, `chem/`, `align/`, `image/`: file readers and
    writers, crystal symmetry, chemistry, structure alignment, images.
*   `glyph/`, `text/`: the glyph cache and FreeType faces; fonts and label
    drawing.
*   `python/`: the `_cmd` extension (`cmd.cc`, with the bindings split into
    `cmd_*.cc` by the Python module that calls them) and the Python glue.
*   `app/`: the `PyMOL` instance lifecycle.

`src/ov/` is PyMOL's container library. Third-party code is in
`third_party/`.

Large files are split by topic within their component, for example
`executive/executive_load.cc` and `select/selector_parse.cc`. Helpers that
several files of one component share go in its `<name>_internal.h`, which
other components must not include.

Components are layered: one may include headers only from its own folder and
from lower layers, listed in `LAYERS` in `scripts/check_deps.py`, which runs
in every build. When a lower component must call up, for example `color`
finding a ramp object, the higher component registers a function pointer at
startup (`ColorSetRampFns`, `SettingSetColorFns`). The components in `CORE`
still depend on each other; when untangling one, move it into `LAYERS`.

### Python modules

The `modules/` directory holds:

*   `pymol/`: the main Python API, split by domain (`importing.py`,
    `editing.py`, `selecting.py`, `viewing.py`, `exporting.py`, `creating.py`,
    `fitting.py` and others). The `_cmd` C extension bridges to the C++ core.
*   `mmol/`: the app-side Python the Swift app calls (panel queries, picking,
    movies, console hints), and `mmol/mcp/`, the MCP server that lets Claude
    drive the app.
*   `chempy/`: chemistry utilities and data structures.
*   `pymol2/`: classes for running PyMOL instances from Python
    (`pymol2.PyMOL`).

### Other key directories

*   `swiftui/`: the macOS app (SwiftUI).
*   `third_party/`: third-party code (CHAMP chemistry, pocketfft, TNT/JAMA
    and the UIUC VMD file-format plugins).
*   `tests/regression/`: the whole-app regression test, its structures
    (`data/`) and reference outputs (`golden/`).
*   `data/`: files bundled into the app and read at runtime through
    `$PYMOL_DATA`: the label fonts (`fonts/`), the sphere meshes
    (`sphere_meshes.bin`), the space-group tables
    (`pymol/space_groups.json`), chemistry data and setting help. Ship
    tables as data files here rather than as source code.

### App conventions

*   **Menus:** use the PyMOL 2.5-style menus in
    `swiftui/PyMOLViewer/Shared/ClassicMenu.swift`, through
    `ClassicMenuTrigger` for buttons and `.classicContextMenu` for
    right-clicks, not SwiftUI `Menu`, `.contextMenu` or confirmation dialogs.
*   **Long operations:** show `RayProgressBar` at the bottom of the viewport.
    `ray` reports `g_mmol_ray_progress`; surface, mesh and dots builds report
    `g_mmol_busy_progress` (`src/pymol/base/busy.h`).
*   **Hot paths:** UI polls that run every frame or every few hundred
    milliseconds must not build selections over the whole session. Read the
    atoms directly in a small C++ function in `src/pymol/python/cmd.cc` (for
    example `object_rep_colors` and `residue_indices`) and keep a thin Python
    wrapper.

### Metal-only build

The core has no OpenGL, WebGL, OpenVR, Windows or Python-less code paths:
PyMOL's OpenGL calls reduce to stubs (`src/pymol/platform/os_gl.h` and
`src/pymol/graphics/gl/gl_stubs.cc`) and `RendererMetal` draws the CGOs. Don't
add `#ifdef`s for other platforms. Core compiler flags and defines live on the
`pymol_core` target in `swiftui/project.yml`, and the app's bridge flags in
`swiftui/PyMOLBridge.xcconfig`.

## Requirements

*   macOS 14 or newer on Apple silicon, Xcode 16 or newer, and xcodegen.
*   Homebrew glm, libpng, freetype and libomp at build time. The app bundles
    the dylibs it needs.
