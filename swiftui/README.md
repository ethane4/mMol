# mMol app build notes

This directory holds the native SwiftUI and Metal app (`PyMOLViewer.xcodeproj`,
generated from `project.yml` by xcodegen). The app embeds a standalone Python
and a Metal-only build of the PyMOL C++ core, which the same project builds as
the `pymol_core` and `champ` static-library targets.

There's no CI for this build yet; it's built and tested by hand.

## Build for macOS

1.  **Install the build-time dependencies.** These headers and libraries are
    needed only to compile the core and link the app. Python comes from the
    next step. OpenGL, GLEW and GLUT aren't needed: the app draws with Metal,
    and the core has no OpenGL code left (`src/pymol/platform/os_gl.h` only
    declares stubs).

    With Homebrew:

    ```shell
    brew install xcodegen glm libpng freetype libomp
    ```

    With MacPorts:

    ```shell
    sudo port install glm libpng freetype libomp
    ```

    Install xcodegen from its [releases](https://github.com/yonaskolb/XcodeGen)
    or with Homebrew.

2.  **Set up the embedded Python.** Run this once, and again to refresh:

    ```shell
    scripts/setup_macos_deps.sh
    ```

    It downloads a `python-build-standalone` CPython 3.13 release into the
    gitignored `deps_macos/python-standalone/python`, then installs the pinned
    NumPy and Biopython that the app ships with.

3.  **Build the core, then the app:**

    ```shell
    swiftui/build_mac.sh                # Release -> build_dd/mMol.app
    CONFIG=Debug swiftui/build_mac.sh   # unoptimized
    ```

    Or run `xcodegen generate` in `swiftui/` and open `PyMOLViewer.xcodeproj`.

## Use a package manager other than Homebrew

`swiftui/PyMOLBridge.xcconfig` holds the build settings for the core and the
app. It looks for glm, libpng, freetype and libomp under Homebrew's
`/opt/homebrew` prefix, through a single override, `PYMOL_EXTERNAL_PREFIX`.

With Homebrew there's nothing to do. With MacPorts, or any other prefix, create
the gitignored local override `swiftui/PyMOLBridge.local.xcconfig`:

```text
PYMOL_EXTERNAL_PREFIX = /opt/local
PYMOL_LIBOMP_STATIC = -L/opt/local/lib/libomp -lomp
```

`PYMOL_LIBOMP_STATIC` can't be derived from `PYMOL_EXTERNAL_PREFIX` because
Homebrew and MacPorts ship different kinds of `libomp`:

*   Homebrew ships a static `libomp.a` under `opt/libomp/lib`.
*   MacPorts ships only a dynamic `libomp.dylib` under `lib/libomp`.

So the variable takes either a full linker flag (`-L<dir> -lomp`) or, for a
prefix laid out like Homebrew, a bare path such as
`/your/prefix/opt/libomp/lib/libomp.a`.
