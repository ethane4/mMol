// RendererMetal: frames, encoders, render targets, matrices and the batch
// pipeline. Drawing and post-processing live in the renderer_metal_*.mm files.

#include "pymol/graphics/metal/renderer_metal.h"

namespace pymol {

#pragma mark - Static matrix helpers
RendererMetal::Mat4 RendererMetal::identityMatrix() {
  Mat4 m{};
  m[0] = m[5] = m[10] = m[15] = 1.0f;
  return m;
}

RendererMetal::Mat4 RendererMetal::multiplyMatrices(const Mat4& a, const Mat4& b) {
  Mat4 r{};
  for (int col = 0; col < 4; ++col) {
    for (int row = 0; row < 4; ++row) {
      float sum = 0.0f;
      for (int k = 0; k < 4; ++k) {
        sum += a[k * 4 + row] * b[col * 4 + k];
      }
      r[col * 4 + row] = sum;
    }
  }
  return r;
}

RendererMetal::Mat4 RendererMetal::translationMatrix(float x, float y, float z) {
  Mat4 m = identityMatrix();
  m[12] = x;
  m[13] = y;
  m[14] = z;
  return m;
}

RendererMetal::Mat4 RendererMetal::scaleMatrix(float x, float y, float z) {
  Mat4 m{};
  m[0] = x;
  m[5] = y;
  m[10] = z;
  m[15] = 1.0f;
  return m;
}
#pragma mark - Metal primitive type mapping
MTLPrimitiveType RendererMetal::toMTL(PrimitiveType t) {
  switch (t) {
    case PrimitiveType::Points:
      return MTLPrimitiveTypePoint;
    case PrimitiveType::Lines:
      return MTLPrimitiveTypeLine;
    case PrimitiveType::LineStrip:
      return MTLPrimitiveTypeLineStrip;
    case PrimitiveType::Triangles:
      return MTLPrimitiveTypeTriangle;
    case PrimitiveType::TriangleStrip:
      return MTLPrimitiveTypeTriangleStrip;
    case PrimitiveType::TriangleFan:
      // Metal has no triangle fan — callers should convert to triangles
      return MTLPrimitiveTypeTriangle;
    case PrimitiveType::Quads:
      // Metal has no quads — callers should convert to triangles
      return MTLPrimitiveTypeTriangle;
  }
  return MTLPrimitiveTypeTriangle;
}
#pragma mark - Metal blend factor mapping
#pragma mark - Constructor / Destructor
RendererMetal::RendererMetal(id<MTLDevice> device, id<MTLCommandQueue> queue)
    : _device(device),
      _queue(queue),
      _cmdBuffer(nil),
      _encoder(nil),
      _passDesc(nil),
      _drawable(nil),
      _currentPipeline(nil),
      _batchPipeline(nil),
      _vboPipelineUByte(nil),
      _vboPipelineFloat(nil),
      _vboVertexFunc(nil),
      _vboFragmentFunc(nil),
      _vboVertexUnlitFunc(nil),
      _vboFragmentUnlitFunc(nil),
      _vboVertexUnlitFlatFunc(nil),
      _batchBuffer(nil),
      _sphereImpostorPipeline(nil),
      _cylinderImpostorPipeline(nil),
      _depthStencilState(nil),
      _inFlight(std::make_shared<std::atomic<int>>(0)) {
  _modelviewMatrix = identityMatrix();
  _projectionMatrix = identityMatrix();
  std::memset(_uniformData, 0, sizeof(_uniformData));

  // Create initial depth stencil state
  applyDepthStencilState();

  // Build the built-in pipelines (re-callable so setSampleCount can rebuild
  // them at a new MSAA sample count).
  buildBatchPipeline();
  buildVBOPipelines();
}

// Built-in batch pipeline (immediate-mode beginBatch/endBatch + selection
// indicators). Rebuilt by setSampleCount when the MSAA sample count changes.
void RendererMetal::buildBatchPipeline() {
  // MRC: release the previous pipeline before rebuilding (e.g. on MSAA change).
  // On the first call (from the ctor) _batchPipeline is nil, so this is a no-op.
  [_batchPipeline release];
  _batchPipeline = nil;
  NSString* batchSrc = @R"(
using namespace metal;

struct BatchVertexIn {
  float3 position [[attribute(0)]];
  float4 color    [[attribute(1)]];
  float3 normal   [[attribute(2)]];
};

struct BatchUniforms {
  float4x4 modelview;
  float4x4 projection;
  float pointSize;
  float _pad[3];
};

struct BatchVertexOut {
  float4 position [[position]];
  float4 color;
  float pointSize [[point_size]];
};

vertex BatchVertexOut batch_vertex(
    BatchVertexIn in [[stage_in]],
    constant BatchUniforms& uniforms [[buffer(1)]])
{
  BatchVertexOut out;
  out.position = uniforms.projection * uniforms.modelview * float4(in.position, 1.0);
  out.color = in.color;
  out.pointSize = uniforms.pointSize;
  return out;
}

fragment float4 batch_fragment(BatchVertexOut in [[stage_in]],
    constant BatchUniforms& uniforms [[buffer(1)]],
    float2 pc [[point_coord]])
{
  if (uniforms._pad[0] > 0.5) {  /* selection marker (points only) */
    float r = length(pc - 0.5);
    if (r < 0.19) return float4(1.0, 1.0, 1.0, 1.0);
    if (r < 0.31) return float4(0.0, 0.0, 0.0, 1.0);
  }
  return in.color;
}
)";

  NSError* error = nil;
  id<MTLLibrary> batchLib = [_device newLibraryWithSource:batchSrc options:nil error:&error];
  if (batchLib) {
    id<MTLFunction> vertFunc = [batchLib newFunctionWithName:@"batch_vertex"];
    id<MTLFunction> fragFunc = [batchLib newFunctionWithName:@"batch_fragment"];

    if (vertFunc && fragFunc) {
      // Vertex descriptor matching BatchVertex { float3 pos, float4 color, float3 normal }
      MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
      // position: offset 0, float3
      vd.attributes[0].format = MTLVertexFormatFloat3;
      vd.attributes[0].offset = 0;
      vd.attributes[0].bufferIndex = 0;
      // color: offset 12, float4
      vd.attributes[1].format = MTLVertexFormatFloat4;
      vd.attributes[1].offset = 3 * sizeof(float);
      vd.attributes[1].bufferIndex = 0;
      // normal: offset 28, float3
      vd.attributes[2].format = MTLVertexFormatFloat3;
      vd.attributes[2].offset = 7 * sizeof(float);
      vd.attributes[2].bufferIndex = 0;
      // stride = sizeof(BatchVertex) = 10 floats = 40 bytes
      vd.layouts[0].stride = 10 * sizeof(float);
      vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

      MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
      psd.vertexFunction = vertFunc;
      psd.fragmentFunction = fragFunc;
      psd.vertexDescriptor = vd;
      psd.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
      psd.colorAttachments[0].blendingEnabled = YES;
      psd.colorAttachments[0].sourceRGBBlendFactor = MTLBlendFactorSourceAlpha;
      psd.colorAttachments[0].destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
      psd.colorAttachments[0].sourceAlphaBlendFactor = MTLBlendFactorOne;
      psd.colorAttachments[0].destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
      psd.rasterSampleCount = _sampleCount;
      psd.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
      psd.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;

      _batchPipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
      if (!_batchPipeline) {
        NSLog(@"RendererMetal: failed to create batch pipeline: %@", error);
      }
      [psd release];  // MRC: descriptor (+1) consumed by pipeline creation
    }
    // MRC: the functions (+1) and library (+1) aren't needed past pipeline build.
    [vertFunc release];
    [fragFunc release];
    [batchLib release];
  } else {
    NSLog(@"RendererMetal: failed to compile batch shader: %@", error);
  }
}

void RendererMetal::setSampleCount(NSUInteger n) {
  if (n < 1) n = 1;
  // The MSAA scene pass resolves the Depth32Float_Stencil8 depth attachment
  // (StoreAndMultisampleResolve, see ensurePostTargets). A depth-stencil
  // format is only a legal MSAA resolve target on Apple GPU family 5+ (A12 and
  // later) and on Mac GPUs; other GPUs reject it —
  // under Metal API validation (any run from Xcode) that is a hard assert in
  // beginFrame on the very first frame, i.e. a crash on launch. Fall back to
  // single-sample rendering there rather than resolving an unsupported format.
  if (n > 1 && _device) {
    bool depthStencilResolveOK = false;
    if (@available(macOS 10.15, *)) {
      depthStencilResolveOK =
          [_device supportsFamily:MTLGPUFamilyApple5] || [_device supportsFamily:MTLGPUFamilyMac2];
    }
    if (!depthStencilResolveOK) {
      static bool warned = false;
      if (!warned) {
        NSLog(@"RendererMetal: MSAA depth-stencil resolve unsupported on this GPU; "
              @"rendering single-sample (metal_msaa ignored)");
        warned = true;
      }
      n = 1;
    }
  }
  if (n == _sampleCount) return;
  _sampleCount = n;
  // Rebuild every sample-count-dependent (opaque-pass) pipeline. OIT and
  // post-process pipelines stay single-sample. Lazy pipelines rebuild on next
  // use once nil'd; batch + VBO are rebuilt eagerly here.
  // MRC: release the +1 sample-count-dependent pipelines before discarding them.
  // Lazy ones (sphere/cylinder/bezier/label/line) rebuild on next use; the eager
  // batch/VBO pipelines are rebuilt by build*Pipelines() below.
  [_sphereImpostorPipeline release];
  _sphereImpostorPipeline = nil;
  releaseCylinderPipelines();
  [_bezierTubePipeline release];
  _bezierTubePipeline = nil;
  [_labelPipeline release];
  _labelPipeline = nil;
  [_connectorPipeline release];
  _connectorPipeline = nil;
  _connectorStride = 0;
  [_vboLinePipeline release];
  _vboLinePipeline = nil;  // lazy; rebuilt on next line draw
  // The interior-cap MARK pipeline is lazy + sample-count-dependent but is NOT
  // rebuilt by buildVBOPipelines (it's stride-keyed, built in drawVBO*). Drop it
  // so it rebuilds at the new count — otherwise its rasterSampleCount would
  // mismatch the scene pass (leak + Metal validation error on capped surfaces).
  [_capMarkPipeline release];
  _capMarkPipeline = nil;
  _capMarkStride = 0;
  // The layout pipeline cache is keyed by sample count; drop it so one-off
  // pipelines rebuild at the new count (release the +1s the cache owns).
  for (auto& kv : _vboPipelineCache) [kv.second release];
  _vboPipelineCache.clear();
  buildBatchPipeline();
  buildVBOPipelines();
  // Force scene-target recreation (single-sample vs multisampled) next frame.
  // MRC: release the old MS targets now; ensurePostTargets recreates them and its
  // own release of these ivars then sees nil (a safe no-op).
  [_sceneColorMS release];
  _sceneColorMS = nil;
  [_sceneDepthMS release];
  _sceneDepthMS = nil;
  _rtW = 0;
  _rtH = 0;
}

RendererMetal::~RendererMetal() {
  // MRC (no ARC, no autorelease pool): every +1-owned Metal/CF object created
  // via alloc/init or a -new* method must be released here, or the entire GPU
  // resource set leaks each time the renderer is torn down (device/context loss,
  // view reinit). NOT released here: _device/_queue (plain ctor assignment, no
  // retain — owned by the Metal layer), the per-frame _cmdBuffer/_encoder/
  // _drawable (autoreleased), and _passDesc/_screenPassDesc/_currentPipeline
  // (non-owning aliases of owned objects).

  // Pooled Metal objects held by value in STL maps: clear()/erase() does NOT
  // send -release under MRC, so release each entry before clearing.
  for (auto& kv : _vboCache) [kv.second release];
  _vboCache.clear();
  for (auto& kv : _vboPipelineCache) [kv.second release];
  _vboPipelineCache.clear();
  for (auto& kv : _buffers) [kv.second release];
  _buffers.clear();
  for (auto& kv : _textures) [kv.second release];
  _textures.clear();
  for (auto& kv : _fbColorAttachments) [kv.second release];
  _fbColorAttachments.clear();
  for (auto& kv : _fbDepthAttachments) [kv.second release];
  _fbDepthAttachments.clear();

  // Offscreen scene/post/OIT targets + atlas + capture copy
  // (mirror the ensurePostTargets release block, which only covers resize).
  [_sceneColor release];
  [_sceneDepth release];
  [_postColor release];
  [_sceneColorMS release];
  [_sceneDepthMS release];
  [_oitAccum release];
  [_oitReveal release];
  [_labelAtlas release];

  // Render-pass descriptors ([[MTLRenderPassDescriptor alloc] init], +1).
  [_scenePassDesc release];
  [_oitPassDesc release];

  // Pipeline states (newRenderPipelineStateWithDescriptor, +1).
  [_batchPipeline release];
  [_vboPipelineUByte release];
  [_vboPipelineFloat release];
  [_sphereImpostorPipeline release];
  [_vboOitPipelineUByte release];
  [_vboOitPipelineFloat release];
  [_sphereOitPipeline release];
  [_oitResolvePipeline release];
  releaseCylinderPipelines();  // per-layout cylinder pipelines (owner)
  [_capMarkPipeline release];
  [_capFillPipeline release];
  [_coveragePipeline release];
  [_surfaceContourPipeline release];
  [_surfaceCoverageTex release];
  [_aoExemptMaskTex release];
  [_aoMaskPipeline release];
  [_aoMaskDepthState release];
  [_vboLinePipeline release];
  [_bezierTubePipeline release];
  [_blitPipeline release];
  [_ssaoPipeline release];
  [_fxaaPipeline release];
  [_outlinePipeline release];
  [_exportAlphaPipeline release];
  [_labelPipeline release];
  [_connectorPipeline release];

  // Shader functions (newFunctionWithName, +1).
  [_vboVertexFunc release];
  [_vboFragmentFunc release];
  [_vboVertexUnlitFunc release];
  [_vboFragmentUnlitFunc release];
  [_vboVertexUnlitFlatFunc release];
  [_vboFragmentOitFunc release];
  [_capMarkVtxFunc release];
  [_capMarkFragFunc release];
  [_capFillVtxFunc release];
  [_capFillFragFunc release];
  [_lineAAVtxFunc release];
  [_lineAAFragFunc release];

  // Samplers (newSamplerStateWithDescriptor, +1).
  [_postSampler release];
  [_labelSampler release];

  // Depth-stencil states (newDepthStencilStateWithDescriptor, +1).
  [_depthStencilState release];
  [_capMarkDSS release];
  [_capFillDSS release];
  [_oitDepthState release];
  [_bezierDepthState release];

  [_bezierTessFactors release];

  // Reusable batch buffer (+1; grown via newBufferWithLength).
  [_batchBuffer release];
}
#pragma mark - Drawable setup
void RendererMetal::beginLiveFrame(int drawableW, int drawableH) {
  // No drawable yet — it is acquired after the scene is encoded (#396). Clear
  // the previous frame's target so a frame that never gets one cannot present
  // stale contents; setPresentTarget fills these in before endFrame.
  _drawable = nil;
  _screenPassDesc = nil;
  // Apply any pending MSAA sample-count change here, before any encoder is open
  // (endFrame ended the previous one and beginFrame hasn't run yet). This
  // rebuilds the opaque pipelines + forces target recreation below, so the new
  // count is consistent for the whole upcoming frame. Early-returns if unchanged.
  setSampleCount(_desiredSampleCount);
  // Render the scene into offscreen targets sized to the drawable; the existing
  // scene-draw code keys off _passDesc, so point it at the offscreen descriptor.
  // The size comes from the caller (MTKView.drawableSize) rather than
  // drawable.texture, precisely so no drawable has to exist yet.
  if (drawableW < 1) drawableW = 1;
  if (drawableH < 1) drawableH = 1;
  ensurePostTargets((NSUInteger)drawableW, (NSUInteger)drawableH);
  // Remember the live render height so pixel-radius post passes (outline)
  // scale up for hi-res offscreen exports (which size _rtH directly via
  // beginOffscreen and never touch _liveRefH). On the live path this keeps
  // pixelRadiusScale()==1, so the on-screen appearance is unchanged.
  _liveRefH = (NSUInteger)drawableH;
  _passDesc = _scenePassDesc;
}

void RendererMetal::setPresentTarget(id<CAMetalDrawable> drawable,
                                     MTLRenderPassDescriptor* passDesc) {
  _drawable = drawable;
  _screenPassDesc = passDesc;  // the drawable target (final composite pass)
}
#pragma mark - Encoder management
void RendererMetal::ensureEncoder() {
  if (_encoder) return;
  if (!_cmdBuffer || !_passDesc) return;

  // Don't create an encoder on an already-committed command buffer.
  // After endFrame() commits the buffer, its status changes from
  // MTLCommandBufferStatusNotEnqueued/Enqueued to Committed/Completed.
  MTLCommandBufferStatus status = [_cmdBuffer status];
  if (status >= MTLCommandBufferStatusCommitted) {
    return;
  }

  _encoder = [_cmdBuffer renderCommandEncoderWithDescriptor:_passDesc];
  if (!_encoder) return;

  // Restore viewport
  [_encoder setViewport:_viewport];

  // Restore scissor if enabled
  if (_scissorEnabled) {
    [_encoder setScissorRect:_scissorRect];
  }

  // Apply depth stencil state
  if (_depthStencilState) {
    [_encoder setDepthStencilState:_depthStencilState];
  }

  // Set cull mode
  [_encoder setCullMode:_cullFaceEnabled ? MTLCullModeBack : MTLCullModeNone];
}

void RendererMetal::applyDepthStencilState() {
  if (!_depthStencilDirty) return;

  MTLDepthStencilDescriptor* desc = [[MTLDepthStencilDescriptor alloc] init];
  desc.depthCompareFunction = _depthTestEnabled ? _depthCompareFunc : MTLCompareFunctionAlways;
  desc.depthWriteEnabled = _depthWriteEnabled;
  // MRC: this runs whenever the depth mode changes (frequently — per rep).
  // Release the previous state before replacing it and release the descriptor,
  // or both leak on every call. Any state still referenced by an in-flight
  // command buffer is kept alive by Metal's own retain.
  [_depthStencilState release];
  _depthStencilState = [_device newDepthStencilStateWithDescriptor:desc];
  [desc release];
  _depthStencilDirty = false;

  if (_encoder) {
    [_encoder setDepthStencilState:_depthStencilState];
  }
}

void RendererMetal::setInteriorCapColor(float r, float g, float b, bool overrideColor) {
  _capColor[0] = r;
  _capColor[1] = g;
  _capColor[2] = b;
  _capColorOverride = overrideColor;
}

void RendererMetal::setRepClip(float front, float back) {
  _repClipFront = front;  // < 0 disables per-rep clip in the lit fragment
  _repClipBack = back;
}

void RendererMetal::setRepContour(bool enabled, const float* rgba, float widthPx) {
  _repContourEnabled = enabled;  // armed => the next surface draw is stashed
  if (enabled) {
    _contourColor[0] = rgba[0];
    _contourColor[1] = rgba[1];
    _contourColor[2] = rgba[2];
    _contourColor[3] = rgba[3];
    _contourWidth = widthPx;
  }
}

void RendererMetal::setRepScreenAO(bool exempt) {
  _repAOExempt = exempt;  // armed => the next lit-VBO draw is stashed for the mask
}

// Rasterize the stashed cartoon/ribbon draws into _aoExemptMaskTex (R8), DEPTH-
// TESTED (LessEqual, no write) against the resolved _sceneDepth so only the
// front-most cartoon pixels are marked (a cartoon pixel hidden behind a stick or
// surface fails the test, so those foreground pixels keep their AO). Reuses the
// coverage vertex/fragment (position-only, writes 1.0). Returns false (mask not
// bound; caller sets aoExemptEnabled=0) when there is nothing to exempt.
// See post_ssao_fog for how the mask gates the SSAO term.
bool RendererMetal::renderAOExemptMask() {
  if (_aoExemptDraws.empty() || !_aoExemptMaskTex || !_coverageVtxFunc || !_coverageFragFunc ||
      !_sceneDepth || !_cmdBuffer)
    return false;

  // Depth-compare state: test vs stored scene depth, never write it back.
  if (!_aoMaskDepthState) {
    MTLDepthStencilDescriptor* dsd = [[MTLDepthStencilDescriptor alloc] init];
    dsd.depthCompareFunction = MTLCompareFunctionLessEqual;
    dsd.depthWriteEnabled = NO;
    _aoMaskDepthState = [_device newDepthStencilStateWithDescriptor:dsd];
    [dsd release];
  }

  MTLRenderPassDescriptor* mpd = [MTLRenderPassDescriptor renderPassDescriptor];
  mpd.colorAttachments[0].texture = _aoExemptMaskTex;
  mpd.colorAttachments[0].loadAction = MTLLoadActionClear;
  mpd.colorAttachments[0].clearColor = MTLClearColorMake(0, 0, 0, 0);
  mpd.colorAttachments[0].storeAction = MTLStoreActionStore;
  // Load the existing scene depth to test against; STORE it unchanged (no depth
  // write happens) so the following SSAO pass still reads valid _sceneDepth.
  mpd.depthAttachment.texture = _sceneDepth;
  mpd.depthAttachment.loadAction = MTLLoadActionLoad;
  mpd.depthAttachment.storeAction = MTLStoreActionStore;
  mpd.stencilAttachment.texture = _sceneDepth;
  mpd.stencilAttachment.loadAction = MTLLoadActionLoad;
  mpd.stencilAttachment.storeAction = MTLStoreActionStore;
  id<MTLRenderCommandEncoder> me = [_cmdBuffer renderCommandEncoderWithDescriptor:mpd];
  if (_aoMaskDepthState) [me setDepthStencilState:_aoMaskDepthState];
  // Slope-scaled negative depth bias: pull the mask fragments slightly toward the
  // camera so grazing cartoon triangles at folds (whose single-sample re-raster
  // depth can exceed the MSAA-resolved scene depth by a sub-texel amount) still
  // pass the LessEqual test and get masked — otherwise an AO contour line leaks
  // through right at the fold. Slope-scaled so nearly-flat faces get almost no
  // bias (can't sneak a cartoon that's genuinely behind a stick into the mask).
  [me setDepthBias:0.0f slopeScale:-3.0f clamp:0.0f];
  [me setCullMode:MTLCullModeNone];
  for (const auto& cd : _aoExemptDraws) {
    if (!_aoMaskPipeline || _aoMaskStride != cd.stride) {
      MTLVertexDescriptor* vd = [[MTLVertexDescriptor alloc] init];
      vd.attributes[0].format = MTLVertexFormatFloat3;
      vd.attributes[0].offset = (NSUInteger)(cd.posOffset < 0 ? 0 : cd.posOffset);
      vd.attributes[0].bufferIndex = 0;
      vd.layouts[0].stride = (NSUInteger)cd.stride;
      vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;
      MTLRenderPipelineDescriptor* pd2 = [[MTLRenderPipelineDescriptor alloc] init];
      pd2.vertexFunction = _coverageVtxFunc;
      pd2.fragmentFunction = _coverageFragFunc;
      pd2.vertexDescriptor = vd;
      pd2.colorAttachments[0].pixelFormat = MTLPixelFormatR8Unorm;
      // Depth-tested against _sceneDepth (Depth32Float_Stencil8).
      pd2.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
      pd2.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
      pd2.rasterSampleCount = 1;
      NSError* pe = nil;
      [_aoMaskPipeline release];  // MRC: release the previous-stride pipeline
      _aoMaskPipeline = [_device newRenderPipelineStateWithDescriptor:pd2 error:&pe];
      _aoMaskStride = cd.stride;
      if (!_aoMaskPipeline) NSLog(@"RendererMetal: AO-exempt mask pipeline failed: %@", pe);
      [vd release];
      [pd2 release];  // MRC: consumed by pipeline creation
    }
    if (!_aoMaskPipeline) continue;
    [me setRenderPipelineState:_aoMaskPipeline];
    [me setVertexBuffer:cd.vbo offset:0 atIndex:0];
    struct {
      float modelview[16];
      float projection[16];
      float pointSize;
      float pad[3];
    } cu;
    std::memcpy(cu.modelview, cd.modelview, 64);
    std::memcpy(cu.projection, cd.projection, 64);
    cu.pointSize = 1.0f;
    cu.pad[0] = cu.pad[1] = cu.pad[2] = 0.0f;
    [me setVertexBytes:&cu length:sizeof(cu) atIndex:1];
    {
      struct {
        float front, back, enabled, pad;
      } _cl = {cd.clipFront, cd.clipBack, cd.clipFront >= 0.0f ? 1.0f : 0.0f, 0.0f};
      [me setFragmentBytes:&_cl length:sizeof(_cl) atIndex:1];
    }
    if (cd.ibo) {
      [me drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                     indexCount:(NSUInteger)cd.count
                      indexType:MTLIndexTypeUInt32
                    indexBuffer:cd.ibo
              indexBufferOffset:0];
    } else {
      [me drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:(NSUInteger)cd.count];
    }
  }
  [me endEncoding];
  return true;
}
#pragma mark - Frame lifecycle
void RendererMetal::beginFrame() {
  _cmdBuffer = [_queue commandBuffer];
  _encoder = nil;
  _oitActive = false;
  _oitHasContent = false;
  _coverageDraws.clear();  // surface outer-contour: re-stashed during this frame
  _contourActive = false;
  _aoExemptDraws.clear();  // cartoon/ribbon AO-exempt mask: re-stashed this frame

  // Configure clear values on the render pass descriptor
  if (_passDesc) {
    _passDesc.colorAttachments[0].loadAction = MTLLoadActionClear;
    _passDesc.colorAttachments[0].clearColor =
        MTLClearColorMake(_clearR, _clearG, _clearB, _clearA);
    _passDesc.depthAttachment.loadAction = MTLLoadActionClear;
    _passDesc.depthAttachment.clearDepth = 1.0;
    _passDesc.stencilAttachment.loadAction = MTLLoadActionClear;
    _passDesc.stencilAttachment.clearStencil = 0;

    // Create the encoder immediately to ensure the clear executes
    _encoder = [_cmdBuffer renderCommandEncoderWithDescriptor:_passDesc];
  }
}

void RendererMetal::endFrame() {
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }

  // The scene was rendered to _sceneColor/_sceneDepth; run the post-process
  // chain, whose final pass writes to the drawable.
  //
  // No drawable is a legitimate outcome now that it is acquired after the scene
  // is encoded (#396): `currentDrawable` can time out. Nothing is presented and
  // any pending PNG capture stays pending (runPostChain owns _capturePath), so
  // the caller re-renders — see the forced redraw in MetalViewport's RenderGate.
  if (_drawable && _cmdBuffer && _screenPassDesc && _sceneColor) {
    runPostChain();
    [_cmdBuffer presentDrawable:_drawable];
  } else if (_drawable && _cmdBuffer) {
    [_cmdBuffer presentDrawable:_drawable];
  }

  if (_cmdBuffer) {
    // MMOL_GPU_TIMING=<file>: append rendered frames/s and mean GPU time per
    // frame every ~2 s to <file> (dev diagnostics for "GPU usage is high":
    // distinguishes an expensive frame from a view that keeps re-rendering
    // while static). A file, not NSLog: the unified log does not reliably
    // surface this app's NSLog output when it is launched via `open --env`.
    static const char* gpuTimingPath = getenv("MMOL_GPU_TIMING");
    if (gpuTimingPath && *gpuTimingPath) {
      [_cmdBuffer addCompletedHandler:^(id<MTLCommandBuffer> cb) {
        static double accum = 0.0, windowStart = 0.0;
        static int frames = 0;
        double now = CFAbsoluteTimeGetCurrent();
        if (windowStart == 0.0) windowStart = now;
        accum += (cb.GPUEndTime - cb.GPUStartTime) * 1000.0;
        frames++;
        if (now - windowStart >= 2.0) {
          if (FILE* f = fopen(gpuTimingPath, "a")) {
            fprintf(f, "%.0f frames/s=%.1f gpu_ms/frame=%.2f gpu_busy=%.0f%%\n", now,
                    frames / (now - windowStart), accum / frames,
                    100.0 * (accum / 1000.0) / (now - windowStart));
            fclose(f);
          }
          accum = 0.0;
          frames = 0;
          windowStart = now;
        }
      }];
    }
    // Count this frame as in flight until the GPU reports it done, so the
    // render loop can throttle itself instead of queueing frames faster than
    // they retire (#396). The handler runs on a background thread and can fire
    // after this renderer is gone, hence the shared atomic.
    auto inFlight = _inFlight;
    inFlight->fetch_add(1, std::memory_order_relaxed);
    [_cmdBuffer addCompletedHandler:^(id<MTLCommandBuffer>) {
      inFlight->fetch_sub(1, std::memory_order_relaxed);
    }];
    [_cmdBuffer commit];
    _cmdBuffer = nil;
  }

  _drawable = nil;
  _passDesc = nil;
  _screenPassDesc = nil;
}
#pragma mark - Hi-res offscreen render → PNG
// Size the post targets to an arbitrary W×H and point the scene pass at them.
// There is no drawable: the caller runs beginFrame / scene-draw / endOffscreen,
// and runPostChain (gated on !_offscreen) skips the drawable blit, leaving the
// post-processed color in sceneSrc for the PNG capture block to write out.
void RendererMetal::beginOffscreen(int w, int h, const std::string& path) {
  if (w < 1) w = 1;
  if (h < 1) h = 1;
  _offscreen = true;
  _drawable = nil;
  _screenPassDesc = nil;
  // Force target recreation at the requested resolution, then aim the scene
  // pass + viewport at it. (ensurePostTargets early-returns if size matches.)
  ensurePostTargets((NSUInteger)w, (NSUInteger)h);
  _passDesc = _scenePassDesc;
  setLetterboxOrigin(0, 0);
  viewport(0, 0, w, h);
  requestPNGCapture(path);
}

// End the offscreen render: close the scene encoder, run the full post chain
// (which writes the PNG via the capture block), commit, and block until the
// GPU finishes so the file exists on return. Leaves _offscreen=false; the next
// live beginLiveFrame recreates the targets at the window size.
void RendererMetal::endOffscreen() {
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }
  if (_cmdBuffer && _sceneColor) {
    runPostChain();
  }
  if (_cmdBuffer) {
    [_cmdBuffer commit];
    [_cmdBuffer waitUntilCompleted];  // completion handler writes the PNG
    _cmdBuffer = nil;
  }
  _offscreen = false;
  _passDesc = nil;
}
#pragma mark - Post-processing (offscreen scene target + fullscreen passes)
#pragma mark - Viewport and clear
#pragma mark - State management
#pragma mark - Drawing
#pragma mark - Buffers
#pragma mark - Vertex attributes
#pragma mark - Shaders / Programs
#pragma mark - Textures
#pragma mark - Framebuffers
#pragma mark - Matrix stack
#pragma mark - Batch system (immediate mode replacement)
#pragma mark - Render readiness
#pragma mark - Queries
#pragma mark - Misc
#pragma mark - VBO Pipeline
#pragma mark - VBO Drawing
#pragma mark - Impostor Ray-Casting (analytic spheres)
#pragma mark - Impostor Ray-Casting (analytic cylinders)
#pragma mark - GPU-tessellated Bezier tubes ("tube cartoon")
#pragma mark - Label / Text Drawing
#pragma mark - Label Connectors (background box / outline / connector line)

}  // namespace pymol
