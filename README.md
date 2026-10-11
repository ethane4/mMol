# mMol

<img src="docs/icon.png" alt="mMol icon" width="128">

mMol is a native macOS molecular viewer with the look and feel of desktop
PyMOL 2.5. It runs the real PyMOL engine and its full Python API (drawn with
Metal) and its GPU `ray` gives the same image as PyMOL's CPU
ray tracer.

![mMol showing hemoglobin (4HHB) as a cartoon colored by chain](/docs/screenshot.png)

## Features

*   **PyMOL commands and scripts** Every PyMOL command works at the `mMol>`
    prompt, as do `.pml` and `.py` scripts. Your `~/.pymolrc` and
    `~/.pymolrc.py` run at startup.
*   **Help as you type** The console shows each command's arguments, defaults
    and allowed values (for example the `type=` choices of `fetch`).
    `help <command>` prints readable documentation.
*   **Metal rendering** The live view, picking and `ray` all run on the GPU,
    including cartoon, sticks, spheres, surfaces, meshes, labels and
    transparency. Large structures and multi-state objects stay interactive.
*   **Analysis and tools tabs:**
    *   contact maps and Ramachandran plots
    *   an alignment wizard that reports total and aligned-region RMSD
    *   a mutagenesis wizard
    *   a movie builder (scenes, keyframes, export)
    *   render settings with live sliders
*   **Trajectories and maps** Use `load_traj` for trajectory files (DCD, XTC,
    TRR, GRO and others, through the VMD plugins), and contour CCP4, MRC and
    DX maps.
*   **Figure-ready defaults** Out of the box mMol uses:
    *   an orthoscopic view with no specular highlights;
    *   wide ribbons and curved sheets;
    *   60% transparent surfaces and 50% transparent spheres;
    *   outlined, unshadowed `ray` renders at antialias 2, with transparent PNG
        backgrounds.

    Change any of these with `set`, or in your `~/.pymolrc`.
*   **Claude integration**
    *   Type `claude` at the `mMol>` prompt to open a Claude Code session
        docked under the viewer. It can inspect and change the open session.
    *   Other AI apps can connect through the **mMol MCP** button at the top
        right. It runs a local, token-protected MCP server and offers
        one-click setup for Claude Code and Claude Desktop.

## Install mMol

1.  Download `mMol-<version>.dmg` from the Releases page.
2.  Open it and drag **mMol** into **Applications**.
3.  The first time, right-click **mMol** and choose **Open**, then confirm.

mMol needs macOS 14 (Sonoma) or newer on Apple silicon. The app bundles its own
Python 3.13, NumPy and Biopython, so it needs nothing else installed. Typing
`claude` in the console also needs [Claude Code][claude-code].

[claude-code]: https://docs.claude.com/en/docs/claude-code

## Quick start

```text
fetch 4hhb                    # type=assembly for the biological unit
hide everything
show cartoon
util.cbc                      # color by chain
show sticks, hetatm and not resn HOH
orient
ray 2400, 1600                # GPU ray trace
png ~/Desktop/4hhb.png
```

Files also open when you drag them onto the window or choose **File > Open**:
PDB, mmCIF, SDF, MOL, MOL2, PSE sessions, CCP4, MRC and DX maps, and more.

### Example config file

mMol runs `~/.pymolrc` and `~/.pymolrc.py` at launch, like PyMOL. The repo has
an example, [`examples/.pymolrc`](/examples/.pymolrc), with extra colors, a
`get_lighting` command, and aliases for coloring, figures and quick views
(`clrd`, `clrss`, `fig`, `figure`, `fasta`, `listed` and others). To use it:

```shell
cp examples/.pymolrc ~/.pymolrc
```

Its display settings are mMol's built-in defaults. They're repeated in the file
so the same look carries over to desktop PyMOL.

## Build from source

You need:

*   a Mac with Apple silicon
*   Xcode 16 or newer
*   [Homebrew](https://brew.sh)

Set up once, then build:

```shell
git clone <repo-url> mMol && cd mMol
brew install xcodegen glm libpng freetype libomp
scripts/setup_macos_deps.sh   # Python 3.13, NumPy, Biopython -> deps_macos/
swiftui/build_mac.sh          # optimized build -> build_dd/mMol.app
open build_dd/mMol.app
```

`build_mac.sh` generates the Xcode project from `swiftui/project.yml`, then
builds the PyMOL core and the app. It takes these options:

*   `CONFIG=Debug swiftui/build_mac.sh`: an unoptimized build for debugging.
*   `NAME=mMol-test swiftui/build_mac.sh`: a side copy, so a running mMol
    isn't overwritten.
*   `INSTALL=1 swiftui/build_mac.sh`: also copies the app to `/Applications`.

To work in Xcode, run `xcodegen generate` in `swiftui/`, open
`PyMOLViewer.xcodeproj`, and run the `PyMOLViewer_macOS` scheme.

For MacPorts and other setups without Homebrew, see the
[app build notes](/swiftui/README.md).

### Make a DMG

```shell
swiftui/make_dmg.sh   # build_dmg/mMol-<version>.dmg, ad-hoc signed
```

For a release that opens without the right-click step, sign it with an Apple
Developer ID and notarize it. Store your notary credentials once, then build:

```shell
xcrun notarytool store-credentials mmol-notary
NOTARIZE=1 DEVID="Developer ID Application: Your Name (TEAMID)" \
    swiftui/make_dmg.sh
```

### Run the tests

```shell
xcodebuild test -project swiftui/PyMOLViewer.xcodeproj \
    -scheme UnitTests_macOS -destination 'platform=macOS' \
    -derivedDataPath build_dd -skipPackagePluginValidation
tests/regression/run.sh build_dd/mMol.app
```

The regression test opens the app for a few seconds, renders reference scenes
and runs core commands, and compares the results with `tests/regression/golden`.

## Repository layout

*   `src/pymol`: PyMOL's C++ core, one folder per component (`scene`,
    `objects`, `reps`, `select`, `executive`, `ray`, `graphics/metal` and
    others); `src/ov` is PyMOL's small container library.
*   `modules/pymol`: PyMOL's Python API.
*   `modules/mmol`: mMol's app-side Python (panel queries, picking, movies)
    and the MCP server (`mmol/mcp`) that lets AI apps drive mMol.
*   `swiftui`: the macOS app (SwiftUI), the Xcode project spec, and the build
    and DMG scripts.
*   `data`: fonts, chemistry tables and help text bundled with the app.
*   `third_party`: third-party code (VMD file-format plugins, CHAMP,
    pocketfft, TNT/JAMA linear algebra).
*   `tests/regression`: the whole-app regression test, its structures and
    reference outputs.
*   `scripts`: dependency setup, dylib bundling and the component layering
    check (`check_deps.py`).

## Credits and license

mMol builds on [RayMol](https://github.com/javierbq/RayMol) (javierbq), a Metal
port of the open-source PyMOL, Copyright (c)
[Schrödinger, LLC](https://www.schrodinger.com/). It's distributed under
PyMOL's BSD-style license; see [`LICENSE`](/LICENSE). Third-party code in
`third_party/` keeps its own licenses.

PyMOL is a trademark of Schrödinger, LLC. mMol is an independent project and
isn't affiliated with or endorsed by Schrödinger.

## See also

*   [App build notes](/swiftui/README.md): build details and MacPorts.
*   [`CLAUDE.md`](/CLAUDE.md): architecture notes for working on the code.
*   [`examples/.pymolrc`](/examples/.pymolrc): the example startup file.
