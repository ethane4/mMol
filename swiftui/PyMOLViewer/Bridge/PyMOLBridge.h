// PyMOLBridge.h — C bridging header for Swift ↔ PyMOL embedding API
// This header exposes the C functions Swift needs to drive PyMOL.

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

// Opaque PyMOL instance handle (CPyMOL* in the implementation)
typedef void* PyMOLHandle;

// --- Lifecycle ---
PyMOLHandle PyMOLBridge_New(void);
void PyMOLBridge_Free(PyMOLHandle instance);
void PyMOLBridge_InitPython(PyMOLHandle instance, const char* resourcePath);
void PyMOLBridge_Start(PyMOLHandle instance);

// --- Render loop ---
int PyMOLBridge_Idle(PyMOLHandle instance);
void PyMOLBridge_Reshape(PyMOLHandle instance, int width, int height);
int PyMOLBridge_GetRedisplay(PyMOLHandle instance, int reset);

// --- Input ---
void PyMOLBridge_Button(PyMOLHandle instance, int button, int state, int x,
                        int y, int modifiers);
void PyMOLBridge_Drag(PyMOLHandle instance, int x, int y, int modifiers);
void PyMOLBridge_SetLetterboxAspect(PyMOLHandle instance, float aspect);
// RGB (0..1) of the 3D selection indicator squares — set from the active theme.
void PyMOLBridge_SetSelectionColor(PyMOLHandle instance, float r, float g,
                                   float b);
// RGB (0..1) of the transient hover-preview indicator points (issue #165).
void PyMOLBridge_SetPreselectionColor(PyMOLHandle instance, float r, float g,
                                      float b);
void PyMOLBridge_CapturePNG(PyMOLHandle instance, const char* path);
// Hi-res offscreen render → PNG: reshape PyMOL to width×height, render the full
// Metal pipeline (all reps) into offscreen targets at that resolution, write
// the PNG, then restore the window size. Synchronous (blocks until the file is
// written). Runs on the main thread, sequenced with the live draw loop.
void PyMOLBridge_RenderHiResPNG(PyMOLHandle instance, const char* path,
                                int width, int height);
// 1 while the core holds a ray/load_png image in place of the live scene
// (Classic layout overlay).
int PyMOLBridge_HasSceneImage(PyMOLHandle instance);
// Discard that image (a click in the viewport dismisses it, as in desktop
// PyMOL).
void PyMOLBridge_DropSceneImage(PyMOLHandle instance);
/// #601: write PNGs uncompressed and unfiltered (process-wide) -- for a movie
/// export's throwaway per-frame files only; turn it back off after the render.
void PyMOLBridge_SetFastPNG(int on);
// Rebuild dirty object representations on the calling thread. Call on the MAIN
// thread before an off-main renderHiResPNG so the rep rebuild (which touches
// the Python C-API via the busy-status callback) doesn't run on the render
// queue.
void PyMOLBridge_UpdateScene(PyMOLHandle instance);
// Fires a cmd.set_key binding by its canonical PyMOL key token ("left",
// "CTRL-T", "ALT-A", "F1"). Returns 1 if a binding existed and was invoked
// (including when the binding's callable raised — the key is still consumed so
// it isn't ALSO dispatched to a menu). Returns 0 when the key is unbound, and
// also 0 when pymol.internal cannot be imported (shutdown / corrupted state) —
// in that case no binding could have run, so the key falls through to macOS
// menus and mMol's own ⌃M / ⌃D shortcuts keep working (#258).
int PyMOLBridge_InvokeKey(const char* key);

// --- Python execution ---
void PyMOLBridge_RunCommand(const char* command);

// Tab autocomplete: runs PyMOL's own command-line completion
// (cmd._parser.complete) on the current input and returns the completed string
// (extended to the unambiguous prefix; the candidate list, when ambiguous, is
// printed to the feedback log). Returns NULL if there's no completion. Caller
// frees with PyMOLBridge_FreeFeedback.
char* PyMOLBridge_Hint(const char* text);
/// The running `ray`'s progress 0..1, or < 0 when none runs (any thread).
float PyMOLBridge_RayProgress(void);
/// Progress 0..1 of the surface / mesh / dots being built (any thread).
float PyMOLBridge_BuildProgress(void);
/// Restart the build progress at 0 before a new build is requested.
void PyMOLBridge_ResetBuildProgress(void);
/// The scene background color (bg_rgb), 0..1.
void PyMOLBridge_GetBackground(float rgb[3]);
char* PyMOLBridge_GetFeedback(PyMOLHandle instance);
void PyMOLBridge_FreeFeedback(char* str);

// Construct the Metal renderer from the MTKView (idempotent), and hand off the
// per-frame drawable + render-pass descriptor (mtkView/drawable/passDescriptor
// passed as opaque void* so this C header needs no Metal import).
void PyMOLBridge_SetupMetalRenderer(PyMOLHandle instance, void* mtkView);

// Render one live frame at width x height (backing pixels).
//
// The CAMetalDrawable is acquired HERE, from the view, only once the scene has
// been encoded — not by the caller beforehand (#396). Only the final post pass
// writes to the drawable, so asking for it earlier just parked the main thread
// in `currentDrawable` (and stalled input) for as long as the GPU was behind.
//
// Returns 1 if a frame was presented, 0 if the drawable never arrived (the
// offscreen work is still committed; the caller should render again) or if
// there is no renderer yet.
int PyMOLBridge_RenderMetalFrame(PyMOLHandle instance, void* mtkView, int width,
                                 int height);

// Live Metal frames committed but not yet completed on the GPU. The render loop
// skips a tick above a small cap so the late currentDrawable wait stays short.
// 0 when there is no renderer.
int PyMOLBridge_MetalFramesInFlight(PyMOLHandle instance);

// Debug: execute raw Python (PyRun_SimpleString) under the GIL.
void PyMOLBridge_RunPython(const char* code);
// Run read-only Python (the periodic UI polls) WITHOUT leaving a redisplay
// behind: if PyMOL's redisplay flag was clear before, it is cleared again
// afterwards. The polls create/destroy temporary selections (count_atoms,
// count_states, ...), and the core marks the scene changed for those, so an
// idle window was re-rendering a full frame on every poll (~2-3 frames/s
// while nothing moved).
void PyMOLBridge_RunPythonQuiet(PyMOLHandle instance, const char* code);

// Tap-to-select: NDC coords in [-1,1], aspect = width/height. Runs the CPU-side
// metal_pick.pick_at (GL color picking is unavailable on the Metal backend).
void PyMOLBridge_Pick(PyMOLHandle instance, float ndcX, float ndcY,
                      float aspect);

// --- Button/modifier constants (match PyMOL's defines) ---
#define PYMOL_BUTTON_LEFT 0
#define PYMOL_BUTTON_MIDDLE 1
#define PYMOL_BUTTON_RIGHT 2
#define PYMOL_BUTTON_SCROLL_FORWARD 3
#define PYMOL_BUTTON_SCROLL_REVERSE 4
#define PYMOL_BUTTON_DOWN 0
#define PYMOL_BUTTON_UP 1

#define PYMOL_MOD_SHIFT 1
#define PYMOL_MOD_CTRL 2
#define PYMOL_MOD_ALT 4

#ifdef __cplusplus
}
#endif
