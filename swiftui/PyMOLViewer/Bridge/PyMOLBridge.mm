// PyMOLBridge.mm — C bridge implementation connecting Swift to PyMOL's embedding API
// Compiled as Objective-C++ to access both the C++ PyMOL internals and ObjC Metal APIs.
//
// NOTE: PyMOL include paths come from OTHER_CPLUSPLUSFLAGS in
// PyMOLBridge.xcconfig (not HEADER_SEARCH_PATHS, which would join Clang
// module builds of the system frameworks).

#include "pymol/app/pymol.h"
#include "pymol/base/pymol_options.h"
#include "pymol/python/p.h"
#include "pymol/scene/scene.h"
#include "pymol/image/my_png.h"     // MyPNGSetFastWrite (#601)
#include "pymol/scene/scene_def.h"  // CScene::Image / CopyType (classic ray overlay)
#include "pymol/color/color.h"
#include "pymol/settings/setting.h"

#include <cstring>

#import <Foundation/Foundation.h>
#import <Python.h>
#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>

#include "pymol/graphics/metal/renderer_metal.h"

// The bridging header uses PyMOLHandle (void*) to avoid CPyMOL typedef
// conflicts between the bridging header and PyMOL.h. We include PyMOL.h
// directly here (not the bridging header) and cast at function boundaries.
typedef void *PyMOLHandle;
#define INST(h) static_cast<CPyMOL *>(h)

// Forward declarations from PyMOL internals
extern "C" {
PyObject *PyInit__cmd(void);
void init_cmd(void);
PyObject *PyInit__champ(void);  // chempy.champ._champ (statically linked)
}

// C++-linkage forward decl (defined in src/pymol/scene/scene_render.cc). Declared
// OUTSIDE the extern "C" block below so the call resolves to the core's
// C++ (mangled) symbol, not a C one.
void SceneRenderMetal(PyMOLGlobals *G);
void SceneUpdate(PyMOLGlobals *G, int force);  // src/pymol/scene/scene.cc
extern void ImmBatch_SetActiveRenderer(pymol::Renderer *r);

// All PyMOLBridge_* entry points MUST have C linkage to match the Swift
// bridging header (PyMOLBridge.h declares them inside extern "C"); otherwise
// the Swift side references _PyMOLBridge_* while the .mm emits mangled C++
// names and the link fails.
extern "C" {

// --- _mmol_render: Python access to the core's GPU `ray` and PNG writer ---

// GPU port of PyMOL's ray tracer (src/pymol/ray/ray_metal.mm, in the core).
void RayMetalSetEnabled(int on);
int RayMetalEnabled();
const char *RayMetalLastReport();
const char *RayMetalStages();

static PyObject *mmolRender_fastPng(PyObject *, PyObject *args) {
  int on = 0;
  if (!PyArg_ParseTuple(args, "i", &on)) return nullptr;
  MyPNGSetFastWrite(on != 0);
  Py_RETURN_NONE;
}

static PyObject *mmolRender_gpuRayStages(PyObject *, PyObject *) {
  return PyUnicode_FromString(RayMetalStages());
}

static PyObject *mmolRender_setGpuRay(PyObject *, PyObject *args) {
  int on = 0;
  if (!PyArg_ParseTuple(args, "i", &on)) return nullptr;
  int was = RayMetalEnabled();
  RayMetalSetEnabled(on);
  return PyLong_FromLong(was);
}

static PyObject *mmolRender_gpuRayReport(PyObject *, PyObject *) {
  return PyUnicode_FromString(RayMetalLastReport());
}

static PyMethodDef mmolRenderMethods[] = {
    {"set_gpu_ray", mmolRender_setGpuRay, METH_VARARGS,
     "set_gpu_ray(on) -> previous: run PyMOL's `ray` tracer on the GPU when the scene allows"},
    {"fast_png", mmolRender_fastPng, METH_VARARGS,
     "fast_png(on): write PNGs uncompressed (for throwaway files)"},
    {"gpu_ray_stages", mmolRender_gpuRayStages, METH_NOARGS,
     "gpu_ray_stages() -> per-stage times of the last ray render"},
    {"gpu_ray_report", mmolRender_gpuRayReport, METH_NOARGS,
     "gpu_ray_report() -> how the last `ray` was traced (gpu: … / cpu: reason)"},
    {nullptr, nullptr, 0, nullptr}};

static struct PyModuleDef mmolRenderModule = {PyModuleDef_HEAD_INIT, "_mmol_render", nullptr, -1,
                                              mmolRenderMethods};

static PyObject *PyInit__mmol_render(void) { return PyModule_Create(&mmolRenderModule); }

// The ray image PyMOL is currently showing in place of the live scene, if any:
// set by `ray` / `load_png`, dropped by the core (SceneInvalidateCopy) on any
// change to the scene or view. Returns 1 while one is live, else 0.
int PyMOLBridge_HasSceneImage(PyMOLHandle h) {
  if (!h) return 0;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G || !G->Scene) return 0;
  return (G->Scene->Image && G->Scene->CopyType) ? 1 : 0;
}

// Drop the held ray image (a click in the viewport, as in desktop PyMOL).
void PyMOLBridge_DropSceneImage(PyMOLHandle h) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (G && G->Scene && G->Scene->CopyType) SceneInvalidate(G);
}

// --- Lifecycle ---

PyMOLHandle PyMOLBridge_New(void) {
  CPyMOLOptions *options = PyMOLOptions_New();
  if (!options) return nullptr;

  options->show_splash = 1;
  options->internal_gui = 0;
  options->internal_feedback = 1;
  options->external_gui = 0;

  CPyMOL *instance = PyMOL_NewWithOptions(options);
  PyMOLOptions_Free(options);
  return static_cast<PyMOLHandle>(instance);
}

void PyMOLBridge_Free(PyMOLHandle h) {
  if (h) {
    PyMOL_Stop(INST(h));
    PyMOL_Free(INST(h));
  }
}

void PyMOLBridge_InitPython(PyMOLHandle h, const char *resourcePath) {
  if (!h || !resourcePath) return;

  // Register the statically-linked _cmd builtin BEFORE init (top-level name only).
  PyImport_AppendInittab("_cmd", PyInit__cmd);
  // Same for _champ (chemistry/charge assignment for vacuum electrostatics
  // etc.). inittab only takes top-level names; chempy imports it as the
  // submodule chempy.champ._champ, so we pre-seed sys.modules after init.
  PyImport_AppendInittab("_champ", PyInit__champ);
  PyImport_AppendInittab("_mmol_render", PyInit__mmol_render);

  NSString *resPath = [NSString stringWithUTF8String:resourcePath];
  NSString *pythonHome =
      [resPath stringByAppendingPathComponent:@"python"];  // contains lib/python3.13
  NSString *modulesPath = [resPath stringByAppendingPathComponent:@"modules"];
  NSString *dataPath = [resPath stringByAppendingPathComponent:@"data"];

  // Modern PyConfig boot (mirrors layer5/main_appkit.mm). NOT isolated: PyMOL
  // relies on a normally-populated sys.path + site.py. config.home must be the
  // directory CONTAINING lib/python3.13 (BeeWare layout), i.e. <res>/python.
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  config.isolated = 0;
  config.site_import = 1;
  config.write_bytecode = 0;  // signed/read-only bundle: cannot write .pyc
  config.buffered_stdio = 0;
  PyConfig_SetBytesString(&config, &config.program_name, "PyMOL");
  PyConfig_SetBytesString(&config, &config.home, [pythonHome UTF8String]);

  PyStatus status = Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);
  if (PyStatus_Exception(status)) {
    NSLog(@"[PyMOL] Python init failed: %s", status.err_msg ? status.err_msg : "(unknown)");
    return;
  }

  init_cmd();  // register pymol._cmd in sys.modules

  PyObject *sysPath = PySys_GetObject("path");
  if (sysPath) {
    PyObject *p = PyUnicode_FromString([modulesPath UTF8String]);
    PyList_Insert(sysPath, 0, p);
    Py_DECREF(p);
  }

  setenv("PYMOL_PATH", [resPath UTF8String], 1);
  setenv("PYMOL_DATA", [dataPath UTF8String], 1);

  // _champ is a TOP-LEVEL builtin (inittab), but chempy/champ/__init__.py does
  // `from . import _champ` (i.e. imports chempy.champ._champ). Pre-seed
  // sys.modules so that submodule import resolves to the builtin — otherwise
  // charge assignment (util.protein_vacuum_esp etc.) raises ImportError.
  PyRun_SimpleString(
      "import sys, importlib\n"
      "try:\n"
      "    sys.modules.setdefault('chempy.champ._champ', importlib.import_module('_champ'))\n"
      "except Exception as _e:\n"
      "    import os; os.write(2, ('[PyMOL] _champ seed failed: %r\\n' % (_e,)).encode())\n");

  // NOTE: PInit / PyMOL_Start / stage-1 happen in PyMOLBridge_Start, AFTER the
  // C subsystems exist. Calling them here (before PyMOL_Start) dereferences
  // uninitialized button-mode/Setting state and crashes.
}

void PyMOLBridge_Start(PyMOLHandle h) {
  if (!h) return;
  // Mirror the macOS embedding sequence (layer5/main_appkit.mm): PyMOL_Start
  // brings up all C subsystems (Setting/ButMode/Scene/...), THEN PInit wires
  // the Python layer as the global instance, THEN stage 1 enables deferred
  // command processing in PyMOL_Idle. (PyMOL_StartWithPython would PInit with
  // global_instance=false; macOS uses true, so we do the steps explicitly.)
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  // Bind the high-level pymol.cmd singleton to THIS instance's globals BEFORE
  // PInit (mirrors main_appkit.mm:353). PInit wraps &SingletonPyMOLGlobals in a
  // capsule for the Python cmd singleton; without this it has no G and every
  // cmd.* call raises CmdException('G'). Also gates Python stdout->feedback.
  SingletonPyMOLGlobals = G;
  G->HaveGUI = true;
  PyMOL_Start(INST(h));
  PInit(G, true);
  PyMOL_SetPythonInitStage(INST(h), 1);

  // NOTE: ~/.mmolrc(.py) is intentionally NOT loaded here. PyMOLEngine
  // applies the mMol theme defaults (bg color etc.) asynchronously right
  // after isReady flips true (ContentView's `applyPersistedTheme`), which
  // runs AFTER this function returns — loading mmolrc here would have it
  // clobbered by that theme re-assertion moments later. Swift instead calls
  // PyMOLBridge_RunPython("from mmol import pymolrc as mmolrc; mmolrc.load()")
  // at the END of applyPersistedTheme, so a user's startup script always
  // wins over the app's defaults (mMol#225).
}

// --- Render loop ---

int PyMOLBridge_Idle(PyMOLHandle h) { return h ? PyMOL_Idle(INST(h)) : 0; }

void PyMOLBridge_Reshape(PyMOLHandle h, int width, int height) {
  if (h) PyMOL_Reshape(INST(h), width, height, 0);
}

int PyMOLBridge_GetRedisplay(PyMOLHandle h, int reset) {
  return h ? PyMOL_GetRedisplay(INST(h), reset) : 0;
}

// --- Input ---

// Shift incoming mouse coords (full-drawable backing px) into the letterboxed
// scene sub-rect's coordinate frame, so picking/drag line up with the render.
static void lbOffset(PyMOLHandle h, int &x, int &y) {
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  auto *r = static_cast<pymol::RendererMetal *>(G->Renderer);
  if (!r) return;
  x -= r->letterboxOriginX();
  y -= r->letterboxOriginY();
}

void PyMOLBridge_Button(PyMOLHandle h, int button, int state, int x, int y, int modifiers) {
  if (!h) return;
  lbOffset(h, x, y);
  PyMOL_Button(INST(h), button, state, x, y, modifiers);
}

void PyMOLBridge_Drag(PyMOLHandle h, int x, int y, int modifiers) {
  if (!h) return;
  lbOffset(h, x, y);
  PyMOL_Drag(INST(h), x, y, modifiers);
}

void PyMOLBridge_SetLetterboxAspect(PyMOLHandle h, float aspect) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  auto *r = static_cast<pymol::RendererMetal *>(G->Renderer);
  if (r) r->setLetterboxAspect(aspect);
}

void PyMOLBridge_SetSelectionColor(PyMOLHandle h, float r, float g, float b) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G || !G->Renderer) return;
  G->Renderer->setSelectionColor(r, g, b);
}

void PyMOLBridge_SetPreselectionColor(PyMOLHandle h, float r, float g, float b) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G || !G->Renderer) return;
  G->Renderer->setPreselectionColor(r, g, b);
}

void PyMOLBridge_CapturePNG(PyMOLHandle h, const char *path) {
  if (!h || !path) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  auto *r = static_cast<pymol::RendererMetal *>(G->Renderer);
  if (r) r->requestPNGCapture(std::string(path));
}

int PyMOLBridge_InvokeKey(const char *key) {
  if (!key || !key[0]) return 0;
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if (!G) return 0;
  int blk = PAutoBlock(G);
  int fired = 0;
  // pymol.internal._invoke_key(key, quiet=1) — the same routine cmd._special
  // and cmd._ctrl funnel through. Reached via pymol.internal because cmd
  // re-exports only _special/_ctrl/_alt, not _invoke_key. quiet=1 suppresses
  // the "No key mapping" print for the (very common) unbound case.
  PyObject *mod = PyImport_ImportModule("pymol.internal");
  if (mod) {
    PyObject *res = PyObject_CallMethod(mod, "_invoke_key", "si", key, 1);
    if (res) {
      fired = PyObject_IsTrue(res) == 1 ? 1 : 0;
      Py_DECREF(res);
    }
    if (PyErr_Occurred()) {
      // The user's bound callable raised. A binding DID exist, so
      // consume the key (fired = 1) — same policy as any other PyMOL
      // script error. Print so the user sees the traceback.
      PyErr_Print();
      fired = 1;
    }
    Py_DECREF(mod);
  } else {
    // pymol.internal could not be imported (shutdown, corrupted state).
    // No binding could have fired, so leave fired = 0 so the key falls
    // through to the macOS menus. Print the error and clear it.
    PyErr_Print();
  }
  PAutoUnblock(G, blk);
  return fired;
}

// --- Context management ---

// --- Python execution ---

void PyMOLBridge_RunCommand(const char *command) {
  if (!command) return;
  // Route through PyMOL's command interpreter (cmd.do) so input is parsed as
  // PyMOL command language ("load ...", "show cartoon", "python ... python
  // end"), NOT raw Python. Use PyMOL's GIL model (PAutoBlock), NOT
  // PyGILState_Ensure — mixing the two with PyMOL_Idle's manual GIL corrupts
  // the interpreter thread state.
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if (!G) return;
  int blk = PAutoBlock(G);
  PyObject *pymol = PyImport_ImportModule("pymol");
  if (pymol) {
    PyObject *cmd = PyObject_GetAttrString(pymol, "cmd");
    if (cmd) {
      PyObject *res = PyObject_CallMethod(cmd, "do", "s", command);
      Py_XDECREF(res);
      Py_DECREF(cmd);
    }
    Py_DECREF(pymol);
  }
  if (PyErr_Occurred()) PyErr_Print();
  PAutoUnblock(G, blk);
}

extern "C" volatile float g_mmol_ray_progress;
float PyMOLBridge_RayProgress(void) { return g_mmol_ray_progress; }
extern "C" volatile float g_mmol_busy_progress;
float PyMOLBridge_BuildProgress(void) { return g_mmol_busy_progress; }
void PyMOLBridge_ResetBuildProgress(void) { g_mmol_busy_progress = 0.0f; }

void PyMOLBridge_GetBackground(float rgb[3]) {
  rgb[0] = rgb[1] = rgb[2] = 0.0f;
  if (PyMOLGlobals *G = SingletonPyMOLGlobals) {
    const float *c = ColorGet(G, SettingGetGlobal_color(G, cSetting_bg_rgb));
    if (c)
      for (int i = 0; i < 3; ++i) rgb[i] = c[i];
  }
}

char *PyMOLBridge_Hint(const char *text) {
  if (!text) return nullptr;
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if (!G) return nullptr;
  int blk = PAutoBlock(G);
  char *result = nullptr;
  // mmol.console_hint.hint(text): JSON usage line + completions.
  PyObject *mod = PyImport_ImportModule("mmol.console_hint");
  if (mod) {
    PyObject *res = PyObject_CallMethod(mod, "hint", "s", text);
    if (res && PyUnicode_Check(res)) {
      const char *s = PyUnicode_AsUTF8(res);
      if (s) result = strdup(s);
    }
    Py_XDECREF(res);
    Py_DECREF(mod);
  }
  if (PyErr_Occurred()) PyErr_Clear();
  PAutoUnblock(G, blk);
  return result;
}

void PyMOLBridge_RunPython(const char *code) {
  if (!code) return;
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if (!G) return;
  int blk = PAutoBlock(G);
  PyRun_SimpleString(code);
  if (PyErr_Occurred()) PyErr_Print();
  PAutoUnblock(G, blk);
}

void PyMOLBridge_RunPythonQuiet(PyMOLHandle h, const char *code) {
  if (!code) return;
  if (!h) {
    PyMOLBridge_RunPython(code);
    return;
  }
  // Peek (reset=0) so a redisplay that was already pending stays pending.
  int const pendingBefore = PyMOL_GetRedisplay(INST(h), 0);
  // A held `ray` image (Classic PyMOL Layout) survives the poll. The polls'
  // selection churn (count_states temp selections, cmd.select) makes the core
  // invalidate the scene, which drops the image half a second after `ray` —
  // yet by the contract below nothing visible changed, so put it back.
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  CScene *scene = G ? G->Scene : nullptr;
  auto const heldImage = (scene && scene->CopyType) ? scene->Image : nullptr;
  int const heldCopyType = heldImage ? scene->CopyType : 0;
  bool const heldForced = heldImage ? scene->CopyForced : false;
  PyMOLBridge_RunPython(code);
  if (heldImage && !scene->CopyType && !scene->Image) {
    scene->Image = heldImage;
    scene->CopyType = heldCopyType;
    scene->CopyForced = heldForced;
  }
  // Only what the poll itself raised is dropped. The polls are read-only by
  // contract (appkit_inspector.poll_panel / poll, appkit_movie.poll): their
  // temp-selection churn changes nothing visible, so the frame it would
  // trigger is pure waste. Anything a real command changes before or after
  // this call sets the flag on its own and is untouched here.
  if (!pendingBefore) (void)PyMOL_GetRedisplay(INST(h), 1);
}

void PyMOLBridge_Pick(PyMOLHandle h, float ndcX, float ndcY, float aspect) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  char script[256];
  snprintf(script, sizeof(script), "from mmol.pick import pick_at; pick_at(%f, %f, %f)", ndcX, ndcY,
           aspect);
  int blk = PAutoBlock(G);
  PyRun_SimpleString(script);
  if (PyErr_Occurred()) PyErr_Print();
  PAutoUnblock(G, blk);
}

char *PyMOLBridge_GetFeedback(PyMOLHandle h) {
  if (!h) return nullptr;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return nullptr;

  // Use PAutoBlock + PyObject calls. The old PyRun_String(..., PyEval_GetGlobals(),
  // PyEval_GetLocals()) needed an active Python frame, which does not exist when
  // called from the C timer -> 'SystemError: frame does not exist' every tick.
  int blk = PAutoBlock(G);
  char *text = nullptr;
  PyObject *pymol = PyImport_ImportModule("pymol");
  if (pymol) {
    PyObject *cmd = PyObject_GetAttrString(pymol, "cmd");
    if (cmd) {
      PyObject *fb = PyObject_CallMethod(cmd, "_get_feedback", nullptr);
      if (fb) {
        if (PyList_Check(fb) && PyList_Size(fb) > 0) {
          PyObject *sep = PyUnicode_FromString("\n");
          PyObject *joined = PyUnicode_Join(sep, fb);
          Py_XDECREF(sep);
          if (joined) {
            const char *str = PyUnicode_AsUTF8(joined);
            if (str && str[0]) text = strdup(str);
            Py_DECREF(joined);
          }
        }
        Py_DECREF(fb);
      }
      Py_DECREF(cmd);
    }
    Py_DECREF(pymol);
  }
  if (PyErr_Occurred()) PyErr_Clear();
  PAutoUnblock(G, blk);
  return text;
}

void PyMOLBridge_FreeFeedback(char *str) { free(str); }

// --- Metal rendering ---

// Keep the Metal device + queue alive for the process lifetime: libpymol_core
// is built WITHOUT ARC, so RendererMetal::_device/_queue are non-retaining raw
// pointers (main_appkit holds them in long-lived view ivars; we hold them here).
// A locally-created queue would be released by ARC on return, leaving _queue
// dangling and crashing in beginFrame's [_queue commandBuffer].
static id<MTLDevice> g_metalDevice = nil;
static id<MTLCommandQueue> g_metalQueue = nil;

void PyMOLBridge_SetupMetalRenderer(PyMOLHandle h, void *mtkViewPtr) {
  if (!h || !mtkViewPtr) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G || G->Renderer) return;  // idempotent: build the renderer once
  MTKView *v = (__bridge MTKView *)mtkViewPtr;
  // readPixels blits from the drawable texture, so it can't be
  // framebuffer-only. This is benign for the normal render path.
  v.framebufferOnly = NO;
  g_metalDevice = v.device ?: MTLCreateSystemDefaultDevice();
  g_metalQueue = [g_metalDevice newCommandQueue];
  G->HaveGUI = true;  // mirror main_appkit.mm:688
  G->Renderer = new pymol::RendererMetal(g_metalDevice, g_metalQueue);
}

int PyMOLBridge_RenderMetalFrame(PyMOLHandle h, void *mtkViewPtr, int width, int height) {
  if (!h || !mtkViewPtr) return 0;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return 0;
  auto *renderer = static_cast<pymol::RendererMetal *>(G->Renderer);
  if (!renderer) return 0;
  MTKView *view = (__bridge MTKView *)mtkViewPtr;
  if (width < 1 || height < 1) return 0;

  // Keep PyMOL's window/scene-block size in sync with the actual drawable
  // (backing pixels). Without this the scene block stays at the default
  // 640x480, so mouse clicks/drags outside that region miss the Scene block
  // and never rotate/zoom — even though rendering fills the full viewport.
  // G->Option->winX/winY feed SceneRenderMetal's one-time OrthoReshape; set
  // them here (before SceneRenderMetal) so that reshape uses the real size.
  // Letterbox: when a loaded .pse asks for a specific viewport aspect, render
  // the scene into a centered sub-rect of that aspect (PyMOL reshaped to the
  // sub-rect so the projection matches), reproducing the session's saved
  // framing. The surrounding bars are the cleared background. 0 = fill.
  float lbAspect = renderer->letterboxAspect();
  int vpW = width, vpH = height, ox = 0, oy = 0;
  if (lbAspect > 0.0f && width > 0 && height > 0) {
    float winAspect = (float)width / (float)height;
    if (winAspect > lbAspect) {
      vpH = height;
      vpW = (int)(height * lbAspect + 0.5f);
    } else {
      vpW = width;
      vpH = (int)(width / lbAspect + 0.5f);
    }
    ox = (width - vpW) / 2;
    oy = (height - vpH) / 2;
  }

  static int s_lastVpW = 0, s_lastVpH = 0;
  if (vpW != s_lastVpW || vpH != s_lastVpH) {
    G->Option->winX = vpW;
    G->Option->winY = vpH;
    PyMOL_Reshape(INST(h), vpW, vpH, 0);
    s_lastVpW = vpW;
    s_lastVpH = vpH;
  }

  // Mirror main_appkit.mm drawInMTKView (805-869); ordering is load-bearing.
  // beginLiveFrame takes only the SIZE — the drawable is fetched below (#396).
  renderer->beginLiveFrame(width, height);
  renderer->setLetterboxOrigin(ox, oy);
  renderer->viewport(ox, oy, vpW, vpH);
  renderer->beginFrame();
  ImmBatch_SetActiveRenderer(renderer);
  PyMOL_PushValidContext(INST(h));
  SceneRenderMetal(G);
  PyMOL_PopValidContext(INST(h));
  ImmBatch_SetActiveRenderer(nullptr);

  // Now — with the whole scene already encoded into the offscreen targets —
  // ask for the drawable that only the final post pass needs. Any wait the
  // display server imposes has been overlapped with this frame's CPU encoding
  // instead of preceding it. Both properties must come from the same access
  // window; `currentRenderPassDescriptor` is derived from `currentDrawable`,
  // so a nil drawable means we simply present nothing this frame and endFrame
  // commits the offscreen work as-is.
  id<CAMetalDrawable> drawable = view.currentDrawable;
  MTLRenderPassDescriptor *passDesc = drawable ? view.currentRenderPassDescriptor : nil;
  if (drawable && passDesc) renderer->setPresentTarget(drawable, passDesc);
  renderer->endFrame();
  return (drawable && passDesc) ? 1 : 0;
}

int PyMOLBridge_MetalFramesInFlight(PyMOLHandle h) {
  if (!h) return 0;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return 0;
  auto *renderer = static_cast<pymol::RendererMetal *>(G->Renderer);
  return renderer ? renderer->framesInFlight() : 0;
}

// Render one offscreen frame at the current (already-reshaped) size and write
// it to path.
static void renderOneOffscreen(PyMOLHandle h, PyMOLGlobals *G, pymol::RendererMetal *renderer,
                               int width, int height, const std::string &path) {
  renderer->beginOffscreen(width, height, path);
  renderer->beginFrame();
  ImmBatch_SetActiveRenderer(renderer);
  PyMOL_PushValidContext(INST(h));
  SceneRenderMetal(G);
  PyMOL_PopValidContext(INST(h));
  ImmBatch_SetActiveRenderer(nullptr);
  renderer->endOffscreen();  // runs post chain, commits, blocks, writes PNG
}

void PyMOLBridge_SetFastPNG(int on) { MyPNGSetFastWrite(on != 0); }

void PyMOLBridge_RenderHiResPNG(PyMOLHandle h, const char *path, int width, int height) {
  if (!h || !path || width < 1 || height < 1) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  auto *renderer = static_cast<pymol::RendererMetal *>(G->Renderer);
  if (!renderer) return;

  // Save the live window/scene-block size so we can restore it afterward.
  // (The live RenderMetalFrame keys its reshape off these, so leaving them at
  // the hi-res size would desync the on-screen scene block.)
  int savedW = G->Option->winX;
  int savedH = G->Option->winY;

  // Reshape PyMOL (projection + scene block) to the export resolution so the
  // scene is composed for that aspect, then render the full pipeline offscreen.
  G->Option->winX = width;
  G->Option->winY = height;
  PyMOL_Reshape(INST(h), width, height, 0);

  renderOneOffscreen(h, G, renderer, width, height, std::string(path));

  // Restore the window size. The renderer's post targets self-heal on the
  // next setDrawable (drawable size != hi-res size → recreated).
  G->Option->winX = savedW;
  G->Option->winY = savedH;
  PyMOL_Reshape(INST(h), savedW, savedH, 0);
}

// Rebuild any dirty object representations for the current frame/state on the
// CALLING thread. Movie export calls this on the main thread before dispatching
// the GPU render to a background queue: SceneUpdate can trigger
// ObjectMolecule::update() -> OrthoBusyFast -> PLockStatusAttempt, which calls
// into the Python C-API. Under _PYMOL_EMBEDDED, PAutoBlock is a no-op (the main
// thread owns the interpreter's GIL persistently), so that Python call is only
// safe on the main thread — doing it on the render queue corrupts the Python
// heap (SIGSEGV in _PyObject_Malloc). Running the rebuild here leaves the
// subsequent off-main SceneRenderMetal with clean reps and no Python touch.
void PyMOLBridge_UpdateScene(PyMOLHandle h) {
  if (!h) return;
  PyMOLGlobals *G = PyMOL_GetGlobals(INST(h));
  if (!G) return;
  SceneUpdate(G, false);
}

// --- Getters ---

}  // extern "C"
