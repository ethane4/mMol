
// Drawing labels and label connectors.

#include "pymol/graphics/metal/renderer_metal.h"

namespace pymol {

// Self-contained MSL port of data/shaders/label.vs + label.fs. Screen-aligned
// textured glyph quads; the atlas RGBA already carries the baked label color
// (RGB = color, A = glyph coverage). Picking/fog/background-color effects from
// the full GL label shader are intentionally omitted here.
static NSString* const kLabelShaderSrc = @R"(
using namespace metal;

struct LabelVtxIn {
  float3 worldpos          [[attribute(0)]];
  float3 targetpos         [[attribute(1)]];
  float3 screenoffset      [[attribute(2)]];
  float2 texcoords         [[attribute(3)]];
  float3 screenworldoffset [[attribute(4)]];
  float  relative_mode     [[attribute(5)]];
};

struct LabelVtxOut {
  float4 position [[position]];
  float2 texcoords;
};

struct LabelUniforms {
  float4x4 modelview;
  float4x4 projection;
  float2 screenSize;
  float screenOriginVertexScale;
  float scaleByVertexScale;
  float labelTextureSize;
  float front;
  float clipRange;
};

static float convertNormalZToScreenZ(float normalz, float front,
                                     float clipRange, float4x4 proj) {
  float a_centerN = (normalz + 1.0) / 2.0;
  float z = -(front + clipRange * a_centerN);
  float4 p = proj * float4(0.0, 0.0, z, 1.0);
  return p.z / p.w;
}

vertex LabelVtxOut label_vertex(LabelVtxIn in [[stage_in]],
    constant LabelUniforms& U [[buffer(1)]]) {
  LabelVtxOut out;
  float isScreenCoord = step(2.0, fmod(in.relative_mode, 4.0));
  float isPixelCoord  = step(4.0, fmod(in.relative_mode, 8.0));
  float zTarget       = step(8.0, fmod(in.relative_mode, 16.0));
  float isProjected   = step(isScreenCoord + isPixelCoord, 0.5);

  float3 viewVector = (float4(0.0, 0.0, -1.0, 0.0) * U.modelview).xyz;
  float sovx = U.screenOriginVertexScale;
  float screenVertexScale =
      U.scaleByVertexScale * sovx * U.labelTextureSize +
      (1.0 - U.scaleByVertexScale);

  float4 wpos = float4(in.worldpos, 1.0);
  float4 tpos = float4(in.targetpos, 1.0);
  float4 transformedPosition = U.projection * U.modelview * wpos;
  float4 targetPosition = U.projection * U.modelview * tpos;
  targetPosition.xyz /= targetPosition.w; targetPosition.w = 1.0;
  transformedPosition.xyz /= transformedPosition.w;
  transformedPosition.xy =
      (floor(transformedPosition.xy * U.screenSize + 0.5) + 0.5) / U.screenSize;

  float4 a_center = wpos + in.screenworldoffset.z * float4(viewVector, 0.0);
  float4 tposZ = U.projection * U.modelview * a_center;
  tposZ.xyz /= tposZ.w; tposZ.w = 1.0;
  float2 pixOffset = (2.0 * in.worldpos.xy / U.screenSize) - 1.0;
  transformedPosition = isProjected * transformedPosition +
      isScreenCoord * wpos +
      isPixelCoord * float4(pixOffset.x, pixOffset.y, -0.5, 0.0);
  transformedPosition.xy += in.screenworldoffset.xy / (U.screenSize * sovx);
  transformedPosition.z = (1.0 - zTarget) *
          ((isProjected * tposZ.z) + (1.0 - isProjected) *
              convertNormalZToScreenZ(in.worldpos.z, U.front, U.clipRange,
                  U.projection)) +
      zTarget * targetPosition.z;
  transformedPosition.xy +=
      in.screenoffset.xy * 2.0 / (U.screenSize * screenVertexScale);
  transformedPosition.w = 1.0;

  out.position = transformedPosition;
  // The projection matrix is glm's GL-convention one (NDC z in [-1,1]) but
  // Metal's clip volume is [0,w], so remap exactly as every other pipeline in
  // this file does. Without it the near half of the slab is clipped away
  // (screen-anchored labels sit at z = -0.5 and never draw at all) and the
  // rest depth-tests against the scene at half its true depth.
  out.position.z = 0.5 * (out.position.z + out.position.w);
  out.texcoords = in.texcoords;
  return out;
}

fragment float4 label_fragment(LabelVtxOut in [[stage_in]],
    texture2d<float> atlas [[texture(0)]],
    sampler smp [[sampler(0)]]) {
  float4 c = atlas.sample(smp, in.texcoords);
  if (c.a < 0.05)
    discard_fragment();
  return c;
}
)";

void RendererMetal::buildLabelPipeline() {
  if (_labelPipeline) return;

  NSError* error = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kLabelShaderSrc options:nil error:&error];
  if (!lib) {
    NSLog(@"RendererMetal: failed to compile label shader: %@", error);
    return;
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"label_vertex"];
  id<MTLFunction> ffn = [lib newFunctionWithName:@"label_fragment"];
  if (!vfn || !ffn) {
    NSLog(@"RendererMetal: label shader functions not found");
    return;
  }

  // Vertex descriptor matching the interleaved label VBO. Offsets are filled
  // in per-draw via a stored layout; here we build for the canonical packing
  // (worldpos F3, targetpos F3, screenoffset F3, texcoords F2,
  //  screenworldoffset F3, relative_mode F1) with 4-byte alignment, which is
  // how CGOConvertToLabelShader lays it out.
  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  vd.attributes[0].format = MTLVertexFormatFloat3;
  vd.attributes[0].offset = 0;
  vd.attributes[1].format = MTLVertexFormatFloat3;
  vd.attributes[1].offset = 12;
  vd.attributes[2].format = MTLVertexFormatFloat3;
  vd.attributes[2].offset = 24;
  vd.attributes[3].format = MTLVertexFormatFloat2;
  vd.attributes[3].offset = 36;
  vd.attributes[4].format = MTLVertexFormatFloat3;
  vd.attributes[4].offset = 44;
  vd.attributes[5].format = MTLVertexFormatFloat;
  vd.attributes[5].offset = 56;
  for (int i = 0; i < 6; ++i) vd.attributes[i].bufferIndex = 0;
  vd.layouts[0].stride = 60;
  vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

  MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
  psd.vertexFunction = vfn;
  psd.fragmentFunction = ffn;
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

  _labelPipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
  if (!_labelPipeline) {
    NSLog(@"RendererMetal: failed to create label pipeline: %@", error);
    return;
  }

  MTLSamplerDescriptor* sd = [[MTLSamplerDescriptor alloc] init];
  sd.minFilter = MTLSamplerMinMagFilterNearest;
  sd.magFilter = MTLSamplerMinMagFilterNearest;
  sd.sAddressMode = MTLSamplerAddressModeClampToEdge;
  sd.tAddressMode = MTLSamplerAddressModeClampToEdge;
  _labelSampler = [_device newSamplerStateWithDescriptor:sd];
}

void RendererMetal::ensureLabelAtlas(const unsigned char* pixels, int w, int h,
                                     uint64_t generation) {
  if (!pixels || w <= 0 || h <= 0) return;
  if (_labelAtlas && _labelAtlasGen == generation && (int)_labelAtlas.width == w &&
      (int)_labelAtlas.height == h)
    return;

  if (!_labelAtlas || (int)_labelAtlas.width != w || (int)_labelAtlas.height != h) {
    MTLTextureDescriptor* td =
        [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatRGBA8Unorm
                                                           width:w
                                                          height:h
                                                       mipmapped:NO];
    td.usage = MTLTextureUsageShaderRead;
    td.storageMode = MTLStorageModeShared;
    // MRC: the atlas doubles in size whenever it fills up, so this runs more
    // than once per session; release the old one or each grow leaks a texture.
    [_labelAtlas release];
    _labelAtlas = [_device newTextureWithDescriptor:td];
  }
  if (!_labelAtlas) return;
  [_labelAtlas replaceRegion:MTLRegionMake2D(0, 0, w, h)
                 mipmapLevel:0
                   withBytes:pixels
                 bytesPerRow:static_cast<NSUInteger>(w) * 4];
  _labelAtlasGen = generation;
}

void RendererMetal::drawLabels(const LabelDrawCall& call) {
  if (!call.data || call.dataSize == 0 || call.vertexCount <= 0) return;
  ensureEncoder();
  if (!_encoder) return;
  buildLabelPipeline();
  if (!_labelPipeline || !_labelSampler) return;
  ensureLabelAtlas(call.atlasPixels, call.atlasW, call.atlasH, call.atlasGen);
  if (!_labelAtlas) return;

  // Upload (or reuse cached) the interleaved label vertex data.
  id<MTLBuffer> vbo = nil;
  auto it = _vboCache.find(call.data);
  if (it != _vboCache.end()) {
    vbo = it->second;
  } else {
    vbo = [_device newBufferWithBytes:call.data
                               length:call.dataSize
                              options:MTLResourceStorageModeShared];
    if (!vbo) return;
    _vboCache[call.data] = vbo;
  }

  [_encoder setRenderPipelineState:_labelPipeline];
  applyDepthStencilState();
  if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];

  [_encoder setVertexBuffer:vbo offset:0 atIndex:0];

  struct {
    float modelview[16];
    float projection[16];
    float screenSize[2];
    float screenOriginVertexScale;
    float scaleByVertexScale;
    float labelTextureSize;
    float front;
    float clipRange;
    float _pad;
  } U;
  std::memcpy(U.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(U.projection, _projectionMatrix.data(), 64);
  U.screenSize[0] = call.screenW;
  U.screenSize[1] = call.screenH;
  U.screenOriginVertexScale = call.screenOriginVertexScale;
  U.scaleByVertexScale = call.scaleByVertexScale;
  U.labelTextureSize = call.labelTextureSize;
  U.front = call.front;
  U.clipRange = call.clipRange;
  U._pad = 0.f;
  [_encoder setVertexBytes:&U length:sizeof(U) atIndex:1];

  [_encoder setFragmentTexture:_labelAtlas atIndex:0];
  [_encoder setFragmentSamplerState:_labelSampler atIndex:0];

  [_encoder drawPrimitives:MTLPrimitiveTypeTriangle
               vertexStart:0
               vertexCount:static_cast<NSUInteger>(call.vertexCount)];
}
// Vertices the connector vertex shader amplifies each per-connector record
// into. Slots that the record's draw flags don't ask for collapse to a
// degenerate (off-clip) position, so one draw covers every combination:
//
//    0.. 5  background box     (1 quad)
//    6..29  background outline (4 edge quads, 6 verts each)
//   30..35  connector line, segment A
//   36..41  connector line, segment B (label_connector_mode 2/4 only)
//   42..44  connector elbow wedge     (label_connector_mode 2/4 only)
static constexpr int kConnectorVertsPerInstance = 45;

// Inline MSL port of data/shaders/connector.vs + connector.gs + connector.fs.
// The GL path amplifies one vertex per connector into all of the above with a
// geometry shader; Metal has none, so the same math runs here per output
// vertex, selected by [[vertex_id]], with the per-connector record fetched
// per-instance. The transform helpers are line-for-line ports of
// data/shaders/connector.shared.
static NSString* const kConnectorShaderSrc = @R"(
using namespace metal;

struct ConnVtxIn {
  float3 a_target    [[attribute(0)]];
  float3 a_center    [[attribute(1)]];
  float2 a_indent    [[attribute(2)]];
  float3 a_swOffset  [[attribute(3)]];
  float2 a_textSize  [[attribute(4)]];
  float4 a_color     [[attribute(5)]];
  uchar  a_relMode   [[attribute(6)]];
  uchar  a_drawFlags [[attribute(7)]];
  float4 a_bgColor   [[attribute(8)]];
  float  a_relExtLen [[attribute(9)]];
  float  a_conWidth  [[attribute(10)]];
};

struct ConnVtxOut {
  float4 position [[position]];
  float4 color;
  float  lineEdge;
  float  aaCutoff;
};

struct ConnUniforms {
  float4x4 modelview;
  float4x4 projection;
  float2 screenSize;
  float screenOriginVertexScale;
  float textureToLabelSize;
  float front;
  float clipRange;
  float lineSmooth;
  float _pad;
};

// Quad emitted as a 4-vertex triangle strip in the GL geometry shader; these
// are the two triangles that cover it.
constant uint kStrip[6] = {0u, 1u, 2u, 1u, 3u, 2u};

static float4 connNormalize(float4 p, float2 ss) {
  float4 r = float4(p.xyz / p.w, 1.0);
  // match the label shader's rounding to the middle of a pixel
  r.xy = (floor(r.xy * ss + 0.5) + 0.5) / ss;
  return r;
}

static float connNormalZToScreenZ(float nz, float front, float clipRange,
                                  float4x4 proj) {
  float cn = (nz + 1.0) / 2.0;
  float4 p = proj * float4(0.0, 0.0, -(front + clipRange * cn), 1.0);
  return p.z / p.w;
}

// Perpendicular half-offset of a screen-space line, as connector.gs's
// drawLineAsGeometry computes it.
static float2 connLineDirv(float4 p1, float4 p2, float lineWidth, float2 ss) {
  float3 a = float3(floor(p1.xy * ss) / ss, p1.z);
  float3 b = float3(floor(p2.xy * ss) / ss, p2.z);
  float3 diffV = a - b;
  float2 n = normalize(cross(diffV, float3(0.0, 0.0, 1.0))).xy;
  return lineWidth * n / ss;
}

static float2 connOffset0(float2 dv) {
  float hmid = step(abs(dv.x), 0.5);
  float vmid = step(abs(dv.y), 0.5);
  float right = (1.0 - hmid) * step(0.0, dv.x) + hmid * 0.5;
  float top   = (1.0 - vmid) * step(0.0, dv.y) + vmid * 0.5;
  return float2(2.0 * (right - 0.5), 2.0 * (top - 0.5));
}

static float2 connOffset1(float2 dv) {
  float2 n = normalize(dv);
  float absyx = step(abs(n.y), abs(n.x));
  float notabsyx = 1.0 - absyx;
  float hdir = 2.0 * (step(0.0, n.x) - 0.5);
  float vdir = 2.0 * (step(0.0, n.y) - 0.5);
  return float2((absyx * hdir) + (notabsyx * vdir * (dv.x / dv.y)),
                (notabsyx * vdir) + (absyx * hdir * (dv.y / dv.x)));
}

static float2 connOffset2(float2 dv, thread float2& eoff, float2 ts,
                          float extLength) {
  float rightorig = step(0.0, dv.x);
  float ratio = abs(ts.y / ts.x);
  float right = rightorig + (rightorig - 0.5) *
      clamp(abs(dv.x) - 1.0, 0.0, 2.0 * min(1.0, extLength * ratio));
  float top = step(0.0, dv.y);
  float yoff = 2.0 * (top - 0.5);
  eoff = float2(2.0 * (right - 0.5), yoff);
  return float2(2.0 * (rightorig - 0.5), yoff);
}

static float2 connOffset3(float2 dv) {
  float hmid = step(abs(dv.x), 1.0);
  float vmid = step(abs(dv.y), 1.0);
  float right = (1.0 - hmid) * step(0.0, dv.x) + hmid * ((1.0 + dv.x) / 2.0);
  float top   = (1.0 - vmid) * step(0.0, dv.y) + vmid * ((1.0 + dv.y) / 2.0);
  return float2(2.0 * (right - 0.5), 2.0 * (top - 0.5));
}

static float2 connOffset4(float2 dv, thread float2& eoff, float2 ts,
                          float extLength) {
  float hmid = step(abs(dv.x), 1.0);
  float vmid = step(abs(dv.y), 1.0);
  float rightorig =
      (1.0 - hmid) * step(0.0, dv.x) + hmid * ((1.0 + dv.x) / 2.0);
  float ratio = abs(ts.y / ts.x);
  float right = rightorig + (1.0 - hmid) * (rightorig - 0.5) *
      clamp(abs(dv.x) - 1.0, 0.0, 2.0 * min(1.0, extLength * ratio));
  float top = (1.0 - vmid) * step(0.0, dv.y) + vmid * ((1.0 + dv.y) / 2.0);
  float yoff = 2.0 * (top - 0.5);
  eoff = float2(2.0 * (right - 0.5), yoff);
  return float2(2.0 * (rightorig - 0.5), yoff);
}

vertex ConnVtxOut connector_vertex(ConnVtxIn in [[stage_in]],
    uint vid [[vertex_id]],
    constant ConnUniforms& U [[buffer(1)]]) {
  ConnVtxOut out;
  out.position = float4(0.0, 0.0, -2.0, 1.0); // degenerate: clipped away
  out.color = float4(0.0);
  out.lineEdge = 0.0;
  out.aaCutoff = 0.0;

  // --- getDrawFlags (connector.shared) ---
  float df = float(in.a_drawFlags);
  float cm4 = step(64.0, fmod(df, 128.0));
  float cm3 = step(32.0, fmod(df, 64.0));
  float cm2 = step(16.0, fmod(df, 32.0));
  float cm1 = step(8.0,  fmod(df, 16.0));
  float drawBgOutline = step(4.0, fmod(df, 8.0));
  float drawBg        = step(2.0, fmod(df, 4.0));
  float drawLine      = step(1.0, fmod(df, 2.0));
  float cm0 = 1.0 - step(0.5, cm1 + cm2 + cm3 + cm4);

  bool wantBg      = vid < 6u;
  bool wantOutline = vid >= 6u && vid < 30u;
  bool wantLine    = vid >= 30u;
  if ((wantBg && drawBg < 0.5) || (wantOutline && drawBgOutline < 0.5) ||
      (wantLine && drawLine < 0.5))
    return out;

  // --- label size / extension length (connector.gs main) ---
  float t2l = U.textureToLabelSize;
  float negLabelSize = step(t2l, 0.0);
  float extLength = (1.0 - negLabelSize) * in.a_relExtLen +
      negLabelSize * (-in.a_relExtLen / (U.screenOriginVertexScale * 2.0));
  float2 textSize = max(in.a_textSize / abs(t2l), float2(1e-6));

  // --- calculatePreConnectorInfo (connector.shared) ---
  float4 centerPt = float4(in.a_center, 1.0);
  float4 targetPt = float4(in.a_target, 1.0);
  float3 viewVector = (float4(0.0, 0.0, -1.0, 0.0) * U.modelview).xyz;
  float4x4 mvp = U.projection * U.modelview;
  float4 tCenter = connNormalize(mvp * centerPt, U.screenSize);
  float4 tTarget = connNormalize(mvp * targetPt, U.screenSize);
  float4 tPosZ = connNormalize(
      mvp * (centerPt + in.a_swOffset.z * float4(viewVector, 0.0)),
      U.screenSize);

  float rm = float(in.a_relMode);
  float isScreenCoord = step(2.0, fmod(rm, 4.0));
  float isPixelCoord  = step(4.0, fmod(rm, 8.0));
  float isProjected   = step(isPixelCoord + isScreenCoord, 0.5);
  float zTarget       = step(8.0, fmod(rm, 16.0));

  float zValue = (1.0 - zTarget) * ((isProjected * tPosZ.z) +
      (1.0 - isProjected) * connNormalZToScreenZ(in.a_center.z, U.front,
          U.clipRange, U.projection)) + zTarget * tTarget.z;
  zValue += 1e-4; // keep the background box just behind the glyphs

  // Keep an off-slab connector from being drawn across the scene. Done here,
  // while both depths are still in GL NDC, so the blend is meaningful.
  float ttz = (1.0 - zTarget) * (1.0 - isProjected);
  tTarget.z = ttz * tTarget.z + (1.0 - ttz) * zValue;
  // GL NDC z [-1,1] -> Metal [0,1], matching label_vertex (so the box lands
  // with its glyphs) and every other pipeline in this file. Every position
  // emitted below takes its z from one of these two.
  zValue = 0.5 * (zValue + 1.0);
  tTarget.z = 0.5 * (tTarget.z + 1.0);

  float2 pixOffset = (2.0 * in.a_center.xy / U.screenSize) - 1.0;
  tCenter = isProjected * tCenter + isScreenCoord * centerPt +
      isPixelCoord * float4(pixOffset.x, pixOffset.y, -0.5, 0.0);
  tCenter.w = 1.0;

  float2 tsScreen = textSize / U.screenSize;
  float2 offset = in.a_indent * tsScreen +
      in.a_swOffset.xy / (U.screenSize * U.screenOriginVertexScale);
  float2 dVector = (tTarget.xy - tCenter.xy - offset) / tsScreen;
  float2 limit = 1.0 + 0.5 * in.a_conWidth / textSize;
  float doNotDraw = (1.0 - step(limit.x, abs(dVector.x))) *
                    (1.0 - step(limit.y, abs(dVector.y)));

  // --- calculateConnectorInfo (connector.shared) ---
  float2 eoff = float2(0.0);
  float2 conOff;
  if (cm0 > 0.0)      conOff = connOffset0(dVector);
  else if (cm1 > 0.0) conOff = connOffset1(dVector);
  else if (cm2 > 0.0) conOff = connOffset2(dVector, eoff, textSize, extLength);
  else if (cm3 > 0.0) conOff = connOffset3(dVector);
  else                conOff = connOffset4(dVector, eoff, textSize, extLength);

  float4 endpointOnBBX = tCenter;
  endpointOnBBX.xy += offset + tsScreen * conOff;
  endpointOnBBX.z = zValue;
  float4 endpointExt = tCenter;
  endpointExt.xy += offset + tsScreen * eoff;
  endpointExt.z = zValue;

  // --- label bounding box, in NDC ---
  float4 tp3 = tCenter;
  tp3.xy += offset;
  tp3.z = zValue;
  float4 upperLeft  = tp3; upperLeft.xy  += tsScreen * float2(-1.0,  1.0);
  float4 upperRight = tp3; upperRight.xy += tsScreen * float2( 1.0,  1.0);
  float4 lowerLeft  = tp3; lowerLeft.xy  += tsScreen * float2(-1.0, -1.0);
  float4 lowerRight = tp3; lowerRight.xy += tsScreen * float2( 1.0, -1.0);

  // --- background box (inset by the outline width when outlined) ---
  if (wantBg) {
    float2 lw = drawBgOutline * (in.a_conWidth / U.screenSize);
    float4 q[4] = {upperLeft, upperRight, lowerLeft, lowerRight};
    q[0].xy += float2( lw.x, -lw.y);
    q[1].xy += float2(-lw.x, -lw.y);
    q[2].xy += float2( lw.x,  lw.y);
    q[3].xy += float2(-lw.x,  lw.y);
    out.position = q[kStrip[vid]];
    out.color = in.a_bgColor;
    return out;
  }

  // --- background outline: four mitred edge quads ---
  if (wantOutline) {
    uint seg = (vid - 6u) / 6u;
    float4 p1, p2;
    float topext, botext;
    if (seg == 0u) {
      p1 = upperLeft;  p2 = upperRight; topext = 1.0; botext = 0.0;
    } else if (seg == 1u) {
      p1 = lowerRight; p2 = upperRight; topext = 0.0; botext = 1.0;
    } else if (seg == 2u) {
      p1 = lowerLeft;  p2 = lowerRight; topext = 0.0; botext = 1.0;
    } else {
      p1 = lowerLeft;  p2 = upperLeft;  topext = 1.0; botext = 0.0;
    }
    float3 diffV = p1.xyz - p2.xyz;
    float2 linev = in.a_conWidth * normalize(diffV).xy / U.screenSize;
    float2 dirv =
        in.a_conWidth * normalize(float2(diffV.y, -diffV.x)) / U.screenSize;
    float4 v[4] = {p1, p2, p1, p2};
    v[0].xy = p1.xy + dirv + topext * linev;
    v[1].xy = p2.xy + dirv - topext * linev;
    v[2].xy = p1.xy - dirv + botext * linev;
    v[3].xy = p2.xy - dirv - botext * linev;
    out.position = v[kStrip[(vid - 6u) % 6u]];
    out.color = in.a_color;
    return out;
  }

  // --- connector line ---
  if (doNotDraw > 0.5)
    return out; // target is inside the label box: nothing to connect

  // The GL geometry shader widens the line by `antialiasedLines *
  // aa * width`, but nothing ever sets that uniform, so it stays 0 there too.
  float aa = U.lineSmooth * min(1.0, 2.0 / in.a_conWidth);
  float cm24 = cm2 + cm4;

  if (cm24 > 0.0) {
    float2 dirv2 =
        connLineDirv(endpointExt, endpointOnBBX, in.a_conWidth, U.screenSize);
    float2 dirv1 =
        connLineDirv(tTarget, endpointExt, in.a_conWidth, U.screenSize);
    if (vid < 42u) {
      // two straight runs: target -> elbow -> label bounding box
      bool segB = vid >= 36u;
      float4 p1 = segB ? tTarget : endpointExt;
      float4 p2 = segB ? endpointExt : endpointOnBBX;
      float2 dirv = segB ? dirv1 : dirv2;
      float4 v[4] = {p1, p2, p1, p2};
      v[0].xy = p1.xy + dirv;
      v[1].xy = p2.xy + dirv;
      v[2].xy = p1.xy - dirv;
      v[3].xy = p2.xy - dirv;
      uint k = kStrip[(vid - 30u) % 6u];
      out.position = v[k];
      out.lineEdge = (k < 2u) ? -1.0 : 1.0;
      out.aaCutoff = segB ? aa : 0.0;
    } else {
      // wedge that fills the outside of the elbow between the two runs
      float norder = step(0.0, dirv1.x * dirv2.y - dirv1.y * dirv2.x);
      uint e = vid - 42u;
      float4 p = endpointExt;
      if (e > 0u) {
        out.lineEdge = -1.0;
        bool first = (e == 1u);
        if (norder < 0.5) {
          p.xy += first ? dirv1 : dirv2;
          out.aaCutoff = first ? aa : 0.0;
        } else {
          p.xy -= first ? dirv2 : dirv1;
          out.aaCutoff = first ? 0.0 : aa;
        }
      }
      out.position = p;
    }
    out.color = in.a_color;
    return out;
  }

  if (vid >= 36u)
    return out; // single-segment modes use slots 30..35 only

  float4 p1 = tTarget;
  float4 p2 = endpointOnBBX;
  if (cm1 > 0.0 && drawBg > 0.0) {
    // connector_mode 1 over a background box: the GL geometry shader runs the
    // line to the box centre and clips it against the box edges. tp3 IS that
    // centre, so the crossing reduces to a scale along target -> centre.
    float2 rel = tTarget.xy - tp3.xy;
    float2 t = 1.0 - tsScreen / max(abs(rel), float2(1e-6));
    float tEnter = clamp(max(t.x, t.y), 0.0, 1.0);
    p2.xy = tp3.xy + rel * (1.0 - tEnter);
    p2.z = zValue;
  }
  float2 dirv = connLineDirv(p1, p2, in.a_conWidth, U.screenSize);
  float4 v[4] = {p1, p2, p1, p2};
  v[0].xy = p1.xy + dirv;
  v[1].xy = p2.xy + dirv;
  v[2].xy = p1.xy - dirv;
  v[3].xy = p2.xy - dirv;
  uint k = kStrip[vid - 30u];
  out.position = v[k];
  out.lineEdge = (k < 2u) ? -1.0 : 1.0;
  out.aaCutoff = aa;
  out.color = in.a_color;
  return out;
}

fragment float4 connector_fragment(ConnVtxOut in [[stage_in]]) {
  float alpha = in.color.a;
  if (in.aaCutoff > 0.0)
    alpha *= smoothstep(0.0, in.aaCutoff, 1.0 - abs(in.lineEdge));
  if (alpha <= 0.0)
    discard_fragment();
  return float4(in.color.rgb, alpha);
}
)";

void RendererMetal::buildConnectorPipeline(const ConnectorDrawCall& call) {
  if (_connectorPipeline && _connectorStride == call.stride) return;
  [_connectorPipeline release];
  _connectorPipeline = nil;
  _connectorStride = call.stride;

  NSError* error = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kConnectorShaderSrc options:nil error:&error];
  if (!lib) {
    NSLog(@"RendererMetal: failed to compile connector shader: %@", error);
    return;
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"connector_vertex"];
  id<MTLFunction> ffn = [lib newFunctionWithName:@"connector_fragment"];
  [lib release];  // MRC: the functions keep the library alive
  if (!vfn || !ffn) {
    NSLog(@"RendererMetal: connector shader functions not found");
    [vfn release];
    [ffn release];
    return;
  }

  // One record per connector, stepped per instance; the vertex shader turns it
  // into kConnectorVertsPerInstance vertices. Offsets come from the CGO's own
  // interleaved layout (CGOOptimizeConnectors), all 4-byte aligned.
  struct {
    int attr;
    MTLVertexFormat fmt;
    int off;
  } kAttrs[] = {
      {0, MTLVertexFormatFloat3, call.targetOff},
      {1, MTLVertexFormatFloat3, call.centerOff},
      {2, MTLVertexFormatFloat2, call.indentOff},
      {3, MTLVertexFormatFloat3, call.screenWorldOff},
      {4, MTLVertexFormatFloat2, call.textSizeOff},
      {5, MTLVertexFormatUChar4Normalized, call.colorOff},
      {6, MTLVertexFormatUChar, call.relModeOff},
      {7, MTLVertexFormatUChar, call.drawFlagsOff},
      {8, MTLVertexFormatUChar4Normalized, call.bgColorOff},
      {9, MTLVertexFormatFloat, call.relExtLenOff},
      {10, MTLVertexFormatFloat, call.conWidthOff},
  };
  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  for (const auto& a : kAttrs) {
    if (a.off < 0) {
      NSLog(@"RendererMetal: connector VBO is missing attribute %d", a.attr);
      [vfn release];
      [ffn release];
      return;
    }
    vd.attributes[a.attr].format = a.fmt;
    vd.attributes[a.attr].offset = static_cast<NSUInteger>(a.off);
    vd.attributes[a.attr].bufferIndex = 0;
  }
  vd.layouts[0].stride = call.stride;
  vd.layouts[0].stepFunction = MTLVertexStepFunctionPerInstance;
  vd.layouts[0].stepRate = 1;

  MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
  psd.vertexFunction = vfn;
  psd.fragmentFunction = ffn;
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

  _connectorPipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
  [psd release];
  [vfn release];
  [ffn release];
  if (!_connectorPipeline) NSLog(@"RendererMetal: failed to create connector pipeline: %@", error);
}

void RendererMetal::drawConnectors(const ConnectorDrawCall& call) {
  if (!call.data || call.dataSize == 0 || call.connectorCount <= 0) return;
  ensureEncoder();
  if (!_encoder) return;
  buildConnectorPipeline(call);
  if (!_connectorPipeline) return;

  // Re-upload every draw (no pointer-keyed caching), into a transient buffer:
  // a frame issues one connector draw per labelled object, and a shared buffer
  // rewritten at offset 0 would leave every recorded draw reading the last
  // call's data once the command buffer executes -- the same defect fixed for
  // the line-AA path in #462. The encoder retains the buffer for the frame.
  id<MTLBuffer> connVBO = [_device newBufferWithLength:std::max(call.dataSize, (size_t)4096)
                                               options:MTLResourceStorageModeShared];
  if (!connVBO) return;
  std::memcpy(connVBO.contents, call.data, call.dataSize);

  [_encoder setRenderPipelineState:_connectorPipeline];
  applyDepthStencilState();
  if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];
  // Screen-aligned quads with no reliable winding: don't inherit a cull mode.
  [_encoder setCullMode:MTLCullModeNone];
  [_encoder setVertexBuffer:connVBO offset:0 atIndex:0];
  [connVBO release];  // MRC: the encoder holds its own reference for the frame

  struct {
    float modelview[16];
    float projection[16];
    float screenSize[2];
    float screenOriginVertexScale;
    float textureToLabelSize;
    float front;
    float clipRange;
    float lineSmooth;
    float _pad;
  } U;
  std::memcpy(U.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(U.projection, _projectionMatrix.data(), 64);
  U.screenSize[0] = call.screenW;
  U.screenSize[1] = call.screenH;
  U.screenOriginVertexScale = call.screenOriginVertexScale;
  U.textureToLabelSize = call.textureToLabelSize;
  U.front = call.front;
  U.clipRange = call.clipRange;
  U.lineSmooth = call.lineSmooth ? 1.f : 0.f;
  U._pad = 0.f;
  [_encoder setVertexBytes:&U length:sizeof(U) atIndex:1];

  [_encoder drawPrimitives:MTLPrimitiveTypeTriangle
               vertexStart:0
               vertexCount:kConnectorVertsPerInstance
             instanceCount:static_cast<NSUInteger>(call.connectorCount)];

  [_encoder setCullMode:_cullFaceEnabled ? MTLCullModeBack : MTLCullModeNone];
}

}  // namespace pymol
