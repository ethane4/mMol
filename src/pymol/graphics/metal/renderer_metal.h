#pragma once
#include "pymol/graphics/renderer.h"

#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>
#include <simd/simd.h>

#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stack>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pymol {

class RendererMetal : public Renderer {
 public:
  RendererMetal(id<MTLDevice> device, id<MTLCommandQueue> queue);
  ~RendererMetal() override;

  // Live frame setup, deliberately split in two (#396). Only the FINAL post
  // pass writes to the drawable, so the drawable is acquired as late as
  // possible: everything that just needs the target SIZE (MSAA rebuild,
  // offscreen post targets) happens in beginLiveFrame from the drawable's
  // dimensions, and setPresentTarget hands over the actual CAMetalDrawable
  // right before endFrame. Waiting on `currentDrawable` before encoding the
  // scene blocked the main thread (and therefore input) for roughly half of
  // every GPU-bound frame.
  //
  // beginLiveFrame clears any previous drawable, so a frame whose drawable
  // never arrives presents nothing rather than re-presenting a stale one.
  void beginLiveFrame(int drawableW, int drawableH);
  void setPresentTarget(id<CAMetalDrawable> drawable, MTLRenderPassDescriptor* passDesc);

  // Live command buffers that have been committed but whose GPU work has not
  // completed yet. The render loop uses this to skip a tick instead of piling
  // frames onto the queue, which is what keeps the (now late) currentDrawable
  // wait short — see RenderGate in MetalViewport.swift. Offscreen/export
  // frames block until completion and are never counted.
  int framesInFlight() const { return _inFlight ? _inFlight->load() : 0; }

  // Frame lifecycle
  void beginFrame() override;
  void endFrame() override;

  // MSAA: stash the desired sample count (from metal_msaa). Applied at the top
  // of the next beginLiveFrame, before any encoder is open, so a toggle never
  // mismatches an in-flight encoder. n < 1 is clamped to 1.
  void setDesiredSampleCount(int n) override { _desiredSampleCount = (n < 1) ? 1 : (NSUInteger)n; }

  // Viewport and clear
  void viewport(int x, int y, int w, int h) override;
  bool getViewportRect(int& x, int& y, int& w, int& h) const override;
  void clear(bool color, bool depth, bool stencil) override;
  void clearColor(float r, float g, float b, float a) override;
  void scissor(int x, int y, int w, int h) override;

  // State management
  void enable(Capability cap) override;
  void disable(Capability cap) override;
  void blendFunc(BlendFunc src, BlendFunc dst) override;
  void depthFunc(DepthFunc func) override;
  void depthMask(bool write) override;
  void colorMask(bool r, bool g, bool b, bool a) override;
  void lineWidth(float w) override;
  void pointSize(float s) override;
  void setDepthClamp(bool enabled) override;

  // Drawing
  void drawArrays(PrimitiveType mode, int first, int count) override;
  void drawElements(PrimitiveType mode, int count, const void* indices) override;

  // Buffers
  uint32_t createBuffer() override;
  void deleteBuffer(uint32_t id) override;
  void bindBuffer(BufferTarget target, uint32_t id) override;
  void bufferData(BufferTarget target, size_t size, const void* data, BufferUsage usage) override;

  // Vertex attributes
  void vertexAttribPointer(int index, int size, int type, bool normalized, int stride,
                           const void* offset) override;
  void enableVertexAttribArray(int index) override;
  void disableVertexAttribArray(int index) override;

  // Shaders
  void useProgram(uint32_t programId) override;
  void setUniform1i(int location, int v) override;
  void setUniform1f(int location, float v) override;
  void setUniform2f(int location, float v0, float v1) override;
  void setUniform3f(int location, float v0, float v1, float v2) override;
  void setUniform4f(int location, float v0, float v1, float v2, float v3) override;
  void setUniformMatrix4fv(int location, const float* value) override;
  void setUniformMatrix3fv(int location, const float* value) override;

  // Textures
  uint32_t createTexture() override;
  void deleteTexture(uint32_t id) override;
  void bindTexture(TextureTarget target, uint32_t id) override;
  void activeTexture(int unit) override;
  void texParameteri(TextureTarget target, int pname, int param) override;

  // Framebuffers
  uint32_t createFramebuffer() override;
  void deleteFramebuffer(uint32_t id) override;
  void bindFramebuffer(uint32_t id) override;

  // Matrix stack
  void matrixMode(int mode) override;
  void loadIdentity() override;
  void loadMatrixf(const float* m) override;
  void pushMatrix() override;
  void popMatrix() override;
  void translatef(float x, float y, float z) override;
  void scalef(float x, float y, float z) override;
  void multMatrixf(const float* m) override;

  // Immediate mode replacement
  void beginBatch(PrimitiveType mode) override;
  void batchVertex3f(float x, float y, float z) override;
  void batchVertex3fv(const float* v) override;
  void batchVertex2f(float x, float y) override;
  void batchVertex2i(int x, int y) override;
  void batchColor3f(float r, float g, float b) override;
  void batchColor3fv(const float* c) override;
  void batchColor4f(float r, float g, float b, float a) override;
  void batchColor4fv(const float* c) override;
  void batchColor4ub(unsigned char r, unsigned char g, unsigned char b, unsigned char a) override;
  void batchNormal3fv(const float* n) override;
  void endBatch() override;

  // Render readiness
  bool isRenderReady() const override;
  bool hasActiveEncoder() const override;

  // Queries
  void getIntegerv(int pname, int* params) override;
  const char* getString(int name) override;
  int getError() override;

  // Misc
  void flush() override;
  void finish() override;
  void readPixels(int x, int y, int w, int h, int format, int type, void* pixels) override;
  void pixelStorei(int pname, int param) override;

  // VBO rendering
  void drawVBO(PrimitiveType mode, int vertexCount, const void* data, size_t dataSize,
               size_t stride, int posOffset, int normalOffset, int colorOffset, int colorType,
               int interiorCap = 0) override;
  void drawVBOIndexed(PrimitiveType mode, int indexCount, const void* vertexData,
                      size_t vertexDataSize, size_t stride, int posOffset, int normalOffset,
                      int colorOffset, int colorType, const void* indexData, size_t indexDataSize,
                      int interiorCap = 0) override;
  void setInteriorCapColor(float r, float g, float b, bool overrideColor) override;
  void setRepClip(float front, float back) override;
  void setRepContour(bool enabled, const float* rgba, float widthPx) override;
  void setRepScreenAO(bool exempt) override;
  void invalidateVBOCache(uint64_t key) override;
  void invalidateVBOCacheEntry(const void* cpuData) override;
  void drawLabels(const LabelDrawCall& call) override;
  void drawConnectors(const ConnectorDrawCall& call) override;
  void drawSphereImpostors(const SphereImpostorDrawCall& call) override;
  void drawCylinderImpostors(const CylinderImpostorDrawCall& call) override;
  void setPostParams(int fogEnabled, float fogStart, float fogEnd, float bgR, float bgG, float bgB,
                     int aoEnabled, int aaEnabled, int outlineEnabled, float projA, float projB,
                     float outlineR, float outlineG, float outlineB, float outlineWidth,
                     int ortho) override;
  void setLightingParams(float ambient, float direct, float reflect, float specular,
                         float shininess) override;
  void setKeyLightDir(const float* lightv) override;

  // Letterbox: render the scene into a centered sub-rect of the given aspect
  // (W/H) so a loaded .pse reproduces its saved-viewport framing. 0 = fill.
  void setLetterboxAspect(float a) { _letterboxAspect = a; }
  float letterboxAspect() const { return _letterboxAspect; }
  // Letterbox origin in backing px (for mouse-coordinate offsetting).
  int letterboxOriginX() const { return _lbOriginX; }
  int letterboxOriginY() const { return _lbOriginY; }
  void setLetterboxOrigin(int x, int y) {
    _lbOriginX = x;
    _lbOriginY = y;
  }

  // Request a PNG capture of the next RENDERED frame (png ray=0). The Metal
  // app has no GL framebuffer for PyMOL's ScenePNG to read, so we grab the
  // renderer's final offscreen color and write it via MyPNGWrite.
  void requestPNGCapture(const std::string& path) { _capturePath = path; }

  // Hi-res offscreen render → PNG. beginOffscreen sizes the post targets to
  // an arbitrary W×H, points the scene pass at them (no drawable), and arms
  // the PNG capture; the caller then runs the normal beginFrame / scene-draw /
  // endOffscreen sequence. endOffscreen runs the full post chain (skipping the
  // drawable blit), commits, and blocks until the GPU has written the PNG.
  // Targets self-heal to the window size on the next live beginLiveFrame.
  void beginOffscreen(int w, int h, const std::string& path);
  void endOffscreen();
  void beginTransparentOIT() override;
  void endTransparentOIT() override;
  void drawBezierTubes(const void* controlPoints, size_t dataSize, float radius, float r, float g,
                       float b) override;

 private:
  void buildImpostorPipelines();
  // The cylinder VBO layout (stride/offsets/formats) varies with the rep, so
  // the cylinder pipeline is built lazily from the first draw call's layout
  // and rebuilt only if a later call has a different stride.
  void releaseCylinderPipelines();
  void buildCylinderImpostorPipeline(const CylinderImpostorDrawCall& call);
  void buildLabelPipeline();
  // (Re)upload the glyph atlas to an MTLTexture if the generation changed.
  void ensureLabelAtlas(const unsigned char* pixels, int w, int h, uint64_t generation);
  // The connector VBO layout comes from CGOOptimizeConnectors, so the vertex
  // descriptor is built from the first draw call's offsets and rebuilt only if
  // a later call disagrees.
  void buildConnectorPipeline(const ConnectorDrawCall& call);

  // 4x4 matrix stored column-major
  using Mat4 = std::array<float, 16>;

  static Mat4 identityMatrix();
  static Mat4 multiplyMatrices(const Mat4& a, const Mat4& b);
  static Mat4 translationMatrix(float x, float y, float z);
  static Mat4 scaleMatrix(float x, float y, float z);

  void ensureEncoder();
  void applyDepthStencilState();
  // Cached depth-stencil states for the transparent (OIT) and opaque bezier-tube
  // draw paths. Built once on first use and reused across all draws — previously
  // a fresh state (+ descriptor) was created per draw call, leaking one of each
  // every transparent/tube draw under MRC. Released in the destructor.
  id<MTLDepthStencilState> oitDepthState();     // LessEqual, no depth write
  id<MTLDepthStencilState> bezierDepthState();  // LessEqual, depth write
  MTLPrimitiveType toMTL(PrimitiveType t);
  void buildVBOPipelines();

  // Anti-aliased line rendering: scene VBO line segments are expanded on the
  // CPU into feathered screen-space quads (PyMOL "trilines") instead of using
  // 1px MTLPrimitiveTypeLine, so lines/ribbon/mesh-wireframe/dashes are smooth.
  void drawLinesAA(PrimitiveType mode, int vertexCount, const void* data, size_t stride,
                   int posOffset, int colorOffset, int colorType, bool flat);
  void buildLinePipeline();

  // Shared by drawVBO and drawVBOIndexed (ibo nil for the non-indexed draw).
  struct VBOUniforms {
    float modelview[16];
    float projection[16];
    float pointSize;
    float _pad[3];
  };
  id<MTLBuffer> cachedVBOBuffer(const void* data, size_t size);
  // Picks and binds the pipeline for the vertex layout, then the depth state,
  // vertex buffer, uniforms, lighting and per-rep clip. False skips the draw.
  bool bindVBODraw(id<MTLBuffer> vbo, size_t stride, int posOffset, int normalOffset,
                   int colorOffset, int colorType, bool unlit, VBOUniforms& u);
  void encodeVBODraw(PrimitiveType mode, int count, id<MTLBuffer> ibo);
  void drawVBOInteriorCap(PrimitiveType mode, int count, id<MTLBuffer> vbo, id<MTLBuffer> ibo,
                          size_t stride, int posOffset, const VBOUniforms& u);
  void recordVBOCoverage(id<MTLBuffer> vbo, id<MTLBuffer> ibo, int count, size_t stride,
                         int posOffset);

  // Vertex attribute description
  struct VertexAttrib {
    int size = 0;  // component count (1-4)
    int type = 0;  // GL type constant (unused in Metal, always float)
    bool normalized = false;
    int stride = 0;
    uintptr_t offset = 0;
    bool enabled = false;
  };

  // Metal objects
  id<MTLDevice> _device;
  id<MTLCommandQueue> _queue;
  id<MTLCommandBuffer> _cmdBuffer;
  id<MTLRenderCommandEncoder> _encoder;
  MTLRenderPassDescriptor* _passDesc;
  id<CAMetalDrawable> _drawable;
  // In-flight live frame count. Held behind a shared_ptr because the command
  // buffer's completion handler (a background thread) decrements it and may
  // outlive this renderer.
  std::shared_ptr<std::atomic<int>> _inFlight;

  // Buffer pool
  uint32_t _nextBufferId = 1;
  std::unordered_map<uint32_t, id<MTLBuffer>> _buffers;
  uint32_t _boundArrayBuffer = 0;
  uint32_t _boundElementBuffer = 0;

  // Texture pool
  uint32_t _nextTextureId = 1;
  std::unordered_map<uint32_t, id<MTLTexture>> _textures;
  uint32_t _boundTexture = 0;
  int _activeTextureUnit = 0;

  // Framebuffer pool
  uint32_t _nextFBOId = 1;
  std::unordered_map<uint32_t, id<MTLTexture>> _fbColorAttachments;
  std::unordered_map<uint32_t, id<MTLTexture>> _fbDepthAttachments;
  uint32_t _boundFBO = 0;

  // Pipeline state cache
  id<MTLRenderPipelineState> _currentPipeline;
  id<MTLRenderPipelineState> _batchPipeline;     // built-in batch shader pipeline
  id<MTLRenderPipelineState> _vboPipelineUByte;  // VBO with UByte4Norm color
  id<MTLRenderPipelineState> _vboPipelineFloat;  // VBO with Float4 color
  id<MTLFunction> _vboVertexFunc;
  id<MTLFunction> _vboFragmentFunc;
  id<MTLFunction> _vboVertexUnlitFunc;  // flat-color (no normal) for lines/dots
  id<MTLFunction> _vboFragmentUnlitFunc;
  // Unlit, position-ONLY (no per-vertex color attribute): used for uniform-
  // colored line geometry (alignment objects, distance/angle dashes) whose VBO
  // has no a_Color. Reads its color from a uniform (buffer 2) instead of an
  // attribute, so the pipeline can be created (the regular unlit shader requires
  // attribute 2 and would fail to compile a pipeline for such a layout).
  id<MTLFunction> _vboVertexUnlitFlatFunc;
  // Impostor ray-casting (analytic spheres/cylinders). nil-init (MRC).
  id<MTLRenderPipelineState> _sphereImpostorPipeline = nil;
  // Cylinder impostor pipelines are cached PER VERTEX LAYOUT — (stride, a_cap
  // offset) — not in a single slot. a_cap's offset is part of the vertex
  // descriptor, so a stick VBO (per-vertex a_cap) and a CGO VBO (one constant
  // a_cap) need different pipelines even at the same stride; a single slot would
  // recompile the MSL library on every draw in a scene that has both, which is
  // exactly Move mode (issue #441). Values are +1-owned (MRC): release each
  // before erasing, like _vboCache. The two ivars below are NON-OWNING aliases
  // of the entry selected by the last buildCylinderImpostorPipeline().
  struct CylinderPipelines {
    id<MTLRenderPipelineState> opaque = nil;
    id<MTLRenderPipelineState> oit = nil;
  };
  std::map<std::pair<NSUInteger, int>, CylinderPipelines> _cylinderPipelines;
  id<MTLRenderPipelineState> _cylinderImpostorPipeline = nil;  // alias, not owned

  // Post-processing: the scene renders to offscreen color+depth, then
  // fullscreen passes (SSAO, fog/depth-cue, FXAA) composite to the drawable.
  // _passDesc is pointed at _scenePassDesc so existing scene-draw code is
  // unchanged; _screenPassDesc (from Swift) is used only by the final pass.
  id<MTLTexture> _sceneColor = nil;  // resolved (single-sample), post chain reads
  id<MTLTexture> _sceneDepth = nil;  // resolved (single-sample), post chain reads
  id<MTLTexture> _postColor = nil;   // ping-pong target for intermediate passes
  // MSAA: when _sampleCount > 1 the scene renders to these multisampled targets
  // and resolves into _sceneColor/_sceneDepth. Opaque pipelines use
  // _sampleCount; OIT + post pipelines stay single-sample.
  id<MTLTexture> _sceneColorMS = nil;
  id<MTLTexture> _sceneDepthMS = nil;
  NSUInteger _sampleCount = 4;         // 4x MSAA by default (metal_msaa)
  NSUInteger _desiredSampleCount = 4;  // applied at next beginLiveFrame (no encoder)
  void setSampleCount(NSUInteger n);   // rebuilds targets+pipelines on change
  void buildBatchPipeline();
  MTLRenderPassDescriptor* _scenePassDesc = nil;
  MTLRenderPassDescriptor* _screenPassDesc = nil;
  NSUInteger _rtW = 0, _rtH = 0;
  // Live render-target height, captured in beginLiveFrame. Pixel-radius post passes
  // (outline thickness) are authored against the live resolution;
  // offscreen exports run at a higher _rtH, so scaling those radii by
  // _rtH/_liveRefH keeps the effect resolution-relative (WYSIWYG). 0 until the
  // first live frame, where the scale is 1 (live appearance unchanged).
  NSUInteger _liveRefH = 0;
  float pixelRadiusScale() const { return (_liveRefH > 0) ? (float)_rtH / (float)_liveRefH : 1.0f; }
  id<MTLRenderPipelineState> _blitPipeline = nil;
  id<MTLRenderPipelineState> _ssaoPipeline = nil;
  id<MTLRenderPipelineState> _fxaaPipeline = nil;
  id<MTLSamplerState> _postSampler = nil;
  void ensurePostTargets(NSUInteger w, NSUInteger h);
  void buildPostPipelines();
  void runPostChain();

  // Order-independent transparency (weighted-blended). Transparent geometry
  // accumulates here (depth-tested vs opaque _sceneDepth, no write); a resolve
  // pass composites over the opaque color during runPostChain.
  id<MTLTexture> _oitAccum = nil;   // RGBA16Float, additive
  id<MTLTexture> _oitReveal = nil;  // R16Float, revealage (multiplicative)
  MTLRenderPassDescriptor* _oitPassDesc = nil;
  id<MTLRenderPipelineState> _vboOitPipelineUByte = nil;
  id<MTLRenderPipelineState> _vboOitPipelineFloat = nil;
  id<MTLFunction> _vboFragmentOitFunc = nil;  // for one-off OIT pipelines
  // Build a weighted-blended OIT MRT pipeline (vbo_vertex + vbo_fragment_oit)
  // for an arbitrary vertex layout (e.g. the surface's stride-44 layout).
  id<MTLRenderPipelineState> oitPipelineForVD(MTLVertexDescriptor* vd);
  // Build-once cache for one-off VBO pipelines whose vertex layout does not match
  // a prebuilt stride (e.g. the molecular-surface stride-44 layout). Without it,
  // drawVBO/drawVBOIndexed rebuilt a pipeline on EVERY such draw — a per-frame
  // MRC leak plus the (significant) cost of pipeline-state compilation. The cache
  // OWNS each +1 pipeline; callers borrow. Released in setSampleCount + the dtor.
  enum class VBOPipelineVariant { Lit, Unlit, UnlitFlat, Oit };
  id<MTLRenderPipelineState> cachedVBOPipeline(VBOPipelineVariant variant, size_t stride,
                                               int posOffset, int normalOffset, int colorOffset,
                                               int colorType, MTLVertexDescriptor* vd);
  id<MTLRenderPipelineState> _sphereOitPipeline = nil;
  id<MTLRenderPipelineState> _cylinderOitPipeline = nil;  // alias, not owned
  NSUInteger _cylinderOitStride = 0;
  id<MTLRenderPipelineState> _oitResolvePipeline = nil;
  bool _oitActive = false;      // true while the transparent pass is rendering
  bool _oitHasContent = false;  // true if any transparent fragments drew

  // Surface interior-cap (stencil): mark = position-only/no-color/stencil-INVERT,
  // fill = full-screen quad gated on stencil. Built once; mark rebuilt per stride.
  id<MTLRenderPipelineState> _capMarkPipeline = nil;
  id<MTLRenderPipelineState> _capFillPipeline = nil;
  id<MTLDepthStencilState> _capMarkDSS = nil;
  id<MTLDepthStencilState> _capFillDSS = nil;
  size_t _capMarkStride = 0;
  // Interior-cap color (ray_interior_color). Override => use _capColor for caps;
  // else per-primitive default (atom color darkened / surface gray).
  float _capColor[3] = {0.32f, 0.32f, 0.36f};
  bool _capColorOverride = false;
  // Per-rep clip planes (eye-space distances) for the next lit-VBO draw.
  // _repClipFront < 0 => disabled (use the global slab). Set via setRepClip.
  float _repClipFront = -1.0f;
  float _repClipBack = 1e6f;
  // Surface outer-contour outline (per-surface, coverage-boundary). When armed
  // (setRepContour), the next surface draw is stashed; after the scene the
  // stashed geometry is rendered to a coverage mask and a post pass outlines the
  // mask boundary. Works on transparent/clipped surfaces. See drawVBO / runPostChain.
  bool _repContourEnabled = false;        // arm capture for the current draw
  float _contourColor[4] = {0, 0, 0, 1};  // line RGBA (frame param, last writer)
  float _contourWidth = 2.0f;             // px (constant on-screen)
  bool _contourActive = false;            // any contour draw stashed this frame
  id<MTLTexture> _surfaceCoverageTex = nil;
  id<MTLRenderPipelineState> _coveragePipeline = nil;  // stride-keyed
  size_t _coverageStride = 0;
  id<MTLRenderPipelineState> _surfaceContourPipeline = nil;
  id<MTLFunction> _coverageVtxFunc = nil, _coverageFragFunc = nil;
  struct CoverageDraw {
    id<MTLBuffer> vbo;
    id<MTLBuffer> ibo;  // nil => non-indexed
    int count;          // index count (indexed) or vertex count (non-indexed)
    size_t stride;
    int posOffset;
    float modelview[16];
    float projection[16];
    float clipFront;
    float clipBack;
  };
  std::vector<CoverageDraw> _coverageDraws;
  // Per-rep screen-space SSAO exemption (#79). Cartoon/ribbon lit-VBO draws are
  // stashed (when _repAOExempt is armed) and rasterized DEPTH-TESTED against
  // _sceneDepth into _aoExemptMaskTex (R8), so only the front-most cartoon pixels
  // are marked. post_ssao_fog reads that mask and skips the SSAO crease/contour
  // term there, leaving surface pockets untouched.
  bool _repAOExempt = false;                         // armed for the current draw
  std::vector<CoverageDraw> _aoExemptDraws;          // cartoon/ribbon draws this frame
  id<MTLTexture> _aoExemptMaskTex = nil;             // R8 front-most cartoon mask
  id<MTLRenderPipelineState> _aoMaskPipeline = nil;  // stride-keyed (coverage fns)
  size_t _aoMaskStride = 0;
  id<MTLDepthStencilState> _aoMaskDepthState = nil;  // LessEqual, no depth write
  bool renderAOExemptMask();                         // fills _aoExemptMaskTex; see .mm
  id<MTLFunction> _capMarkVtxFunc = nil, _capMarkFragFunc = nil;
  id<MTLFunction> _capFillVtxFunc = nil, _capFillFragFunc = nil;

  // Anti-aliased screen-space line quads. _vboLinePipeline (lazy, sample-count
  // aware) renders CPU-expanded feathered quads. _lineExpand is the per-frame
  // scratch where segments are expanded (9 floats/vert); each draw uploads it
  // into its OWN transient MTLBuffer (see drawVBOLines: a buffer shared across
  // draws is read by all of them at commit time, #462).
  id<MTLRenderPipelineState> _vboLinePipeline = nil;
  id<MTLFunction> _lineAAVtxFunc = nil;
  id<MTLFunction> _lineAAFragFunc = nil;
  std::vector<float> _lineExpand;

  // GPU-tessellated Bezier tube ("tube cartoon") pipeline.
  id<MTLRenderPipelineState> _bezierTubePipeline = nil;
  id<MTLBuffer> _bezierTessFactors = nil;  // MTLQuadTessellationFactorsHalf/patch
  NSUInteger _bezierTessPatchCap = 0;      // patches the factor buffer covers
  void buildBezierTubePipeline();
  // Per-frame post params (fog/depth-cue + SSAO), set by SceneRenderMetal.
  int _postFogEnabled = 0;
  float _fogStart = 0.f, _fogEnd = 1.f;
  float _bgR = 0.f, _bgG = 0.f, _bgB = 0.f;
  int _aoEnabled = 1;       // SSAO (cSetting_metal_ssao)
  int _aaEnabled = 1;       // FXAA (cSetting_antialias_shader != 0)
  int _outlineEnabled = 0;  // silhouette outline (cSetting_metal_outline)
  // Outline contour color + thickness in px (cSetting_metal_outline_color/_width).
  float _outlineR = 0.f, _outlineG = 0.f, _outlineB = 0.f;
  float _outlineWidth = 1.4f;
  id<MTLRenderPipelineState> _outlinePipeline = nil;
  // Offscreen-export-only pass: rewrites the framebuffer alpha from scene depth
  // (background = far -> alpha 0) so a transparent-background PNG can be written
  // on the Metal fast path. Gated on _offscreen && transparent clear (_clearA<0.5).
  id<MTLRenderPipelineState> _exportAlphaPipeline = nil;
  // PyMOL lighting model (cSetting_ambient/direct/reflect/specular/shininess),
  // set per frame by SceneRenderMetal and uploaded into each lit shader's
  // uniform so the Scene-panel lighting sliders actually affect the render.
  // Defaults match the values the shaders previously hard-coded.
  float _lightAmbient = 0.14f, _lightDirect = 0.45f, _lightReflect = 0.481f;
  float _lightSpecular = 0.5f, _lightShininess = 55.0f;
  // Key-light direction TOWARD the light in eye space = -normalize(cSetting_light).
  // Default reproduces the previously hard-coded normalize(0.4,0.4,1.0), which is
  // exactly -normalize(PyMOL's default light). Fed into every lit shader.
  float _keyLightEye[3] = {0.34815531f, 0.34815531f, 0.87038828f};
  float _projA = -1.f, _projB = 0.f;   // projection[10], projection[14]
  float _projOrtho = 0.f;              // 1 = orthographic (linear eye-z recon; #139)
  float _letterboxAspect = 0.f;        // saved-viewport W/H; 0 = fill window
  int _lbOriginX = 0, _lbOriginY = 0;  // letterbox sub-rect origin (backing px)
  std::string _capturePath;            // pending png ray=0 capture (empty = none)
  bool _offscreen = false;             // hi-res offscreen render (no drawable)

  // Label/text rendering (screen-aligned textured glyph quads). Initialized to
  // nil — this is a C++ class under MRC, so id ivars are not zero-initialized.
  id<MTLRenderPipelineState> _labelPipeline = nil;
  id<MTLSamplerState> _labelSampler = nil;
  id<MTLTexture> _labelAtlas = nil;  // glyph atlas, uploaded from CPU copy
  uint64_t _labelAtlasGen = 0;       // generation of the uploaded atlas
  // Label connectors (background box / outline / connector line). The vertex
  // shader amplifies one per-connector record into kConnectorVertsPerInstance
  // vertices, so there is no per-vertex buffer to cache beyond the CGO's own.
  id<MTLRenderPipelineState> _connectorPipeline = nil;
  size_t _connectorStride = 0;  // layout the pipeline was built for
  // The per-connector records are uploaded into a transient per-draw buffer
  // (see drawConnectors). Deliberately NOT the pointer-keyed _vboCache: the
  // connector CGO is rebuilt whenever a label setting changes, and a recycled
  // heap address would serve stale background/connector colors from the cache.
  uint32_t _currentProgram = 0;

  // Depth/stencil state
  id<MTLDepthStencilState> _depthStencilState;
  // Cached transparent/opaque-tube depth-stencil states (see oitDepthState/
  // bezierDepthState). nil-init (MRC: not zero-initialized otherwise).
  id<MTLDepthStencilState> _oitDepthState = nil;
  id<MTLDepthStencilState> _bezierDepthState = nil;
  bool _depthTestEnabled = false;
  bool _depthWriteEnabled = true;
  MTLCompareFunction _depthCompareFunc = MTLCompareFunctionLess;
  bool _depthStencilDirty = true;

  // Blend state (tracked, applied when pipeline is created)
  bool _blendEnabled = false;
  MTLBlendFactor _blendSrcFactor = MTLBlendFactorOne;
  MTLBlendFactor _blendDstFactor = MTLBlendFactorZero;

  // Color mask
  MTLColorWriteMask _colorWriteMask = MTLColorWriteMaskAll;

  // VBO buffer cache — reuse Metal buffers across frames
  std::unordered_map<const void*, id<MTLBuffer>> _vboCache;
  // One-off VBO pipeline cache (see cachedVBOPipeline). Owns each +1 pipeline.
  std::unordered_map<uint64_t, id<MTLRenderPipelineState>> _vboPipelineCache;
  id<MTLBuffer> _batchBuffer;  // reusable buffer for batch drawing

  // Clear values
  float _clearR = 0.0f, _clearG = 0.0f, _clearB = 0.0f, _clearA = 1.0f;

  // Viewport
  MTLViewport _viewport = {0, 0, 1, 1, 0.0, 1.0};
  MTLScissorRect _scissorRect = {0, 0, 1, 1};
  bool _scissorEnabled = false;

  // Capability flags
  bool _cullFaceEnabled = false;
  bool _stencilTestEnabled = false;
  bool _lightingEnabled = false;
  bool _fogEnabled = false;

  // Vertex attributes
  static constexpr int kMaxVertexAttribs = 8;
  VertexAttrib _vertexAttribs[kMaxVertexAttribs];

  // Uniform buffer (generic float storage)
  static constexpr int kMaxUniforms = 64;
  float _uniformData[kMaxUniforms * 4];  // up to 64 vec4s

  // Matrix stack — mode 0 = modelview, mode 1 = projection
  int _matrixMode = 0;  // 0x1700 = GL_MODELVIEW mapped to 0
  Mat4 _modelviewMatrix;
  Mat4 _projectionMatrix;
  std::stack<Mat4> _modelviewStack;
  std::stack<Mat4> _projectionStack;

  // Batch system
  struct BatchVertex {
    float x, y, z;
    float r, g, b, a;
    float nx, ny, nz;
  };

  PrimitiveType _batchMode{};
  std::vector<BatchVertex> _batchVertices;
  float _curR = 1.0f, _curG = 1.0f, _curB = 1.0f, _curA = 1.0f;
  float _curNX = 0.0f, _curNY = 0.0f, _curNZ = 1.0f;

  // Line width / point size
  float _lineWidth = 1.0f;
  float _pointSize = 1.0f;
};

}  // namespace pymol
