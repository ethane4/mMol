
// Drawing vertex buffers: meshes, antialiased lines and their pipelines.

#include "pymol/graphics/metal/renderer_metal.h"

namespace pymol {

void RendererMetal::buildVBOPipelines() {
  // MRC: release the previous build's owned objects before this call rebuilds
  // them (e.g. on an MSAA sample-count change via setSampleCount). On the first
  // call (from the ctor) every ivar is nil, so each release is a no-op.
  [_vboVertexFunc release];
  _vboVertexFunc = nil;
  [_vboFragmentFunc release];
  _vboFragmentFunc = nil;
  [_vboVertexUnlitFunc release];
  _vboVertexUnlitFunc = nil;
  [_vboFragmentUnlitFunc release];
  _vboFragmentUnlitFunc = nil;
  [_vboVertexUnlitFlatFunc release];
  _vboVertexUnlitFlatFunc = nil;
  [_capMarkVtxFunc release];
  _capMarkVtxFunc = nil;
  [_capMarkFragFunc release];
  _capMarkFragFunc = nil;
  [_capFillVtxFunc release];
  _capFillVtxFunc = nil;
  [_capFillFragFunc release];
  _capFillFragFunc = nil;
  [_lineAAVtxFunc release];
  _lineAAVtxFunc = nil;
  [_lineAAFragFunc release];
  _lineAAFragFunc = nil;
  [_vboFragmentOitFunc release];
  _vboFragmentOitFunc = nil;
  [_capMarkDSS release];
  _capMarkDSS = nil;
  [_capFillDSS release];
  _capFillDSS = nil;
  [_capFillPipeline release];
  _capFillPipeline = nil;
  [_vboPipelineUByte release];
  _vboPipelineUByte = nil;
  [_vboPipelineFloat release];
  _vboPipelineFloat = nil;
  [_vboOitPipelineUByte release];
  _vboOitPipelineUByte = nil;
  [_vboOitPipelineFloat release];
  _vboOitPipelineFloat = nil;

  // Metal shader for VBO-based molecular geometry (cartoons, surfaces, etc.)
  // Supports position + normal + color with basic Lambertian lighting.
  NSString* vboSrc = @R"(
using namespace metal;

struct VBOVertexIn {
  float3 position [[attribute(0)]];
  float3 normal   [[attribute(1)]];
  float4 color    [[attribute(2)]];
};

struct VBOUniforms {
  float4x4 modelview;
  float4x4 projection;
  float pointSize;
};

struct VBOVertexOut {
  float4 position [[position]];
  float4 color;
  float3 normalEye;   // eye-space normal, interpolated → per-fragment (Phong)
  float  eyeDist;     // distance from camera (-eyeZ), for per-rep clipping
};

// Per-representation clip planes (eye-space distances from the camera). Lets one
// rep (e.g. the surface) clip tighter than the global slab so the user can peek
// inside while cartoon/sticks stay whole. enabled<0.5 => no per-rep clip.
struct ClipU { float front; float back; float enabled; float _pad; };
static void apply_rep_clip(ClipU clip, float eyeDist) {
  if (clip.enabled > 0.5 && (eyeDist < clip.front || eyeDist > clip.back))
    discard_fragment();
}

// PyMOL two-light model. The ambient/direct/reflect/specular/shininess come
// from the live PyMOL settings (Scene lighting sliders) via LightU, instead of
// hard-coded constants, so the sliders take effect. Two-sided: the interpolated
// normal is flipped to face the viewer so cartoon undersides / surface
// interiors light up instead of going dark.
struct LightU { float ambient, direct, reflect, spec, shininess;
                 float klx, kly, klz; };  // klx/y/z: key-light dir (toward light, eye space)
static float3 vbo_shade(float3 baseColor, float3 nEye, LightU lt) {
  float3 normal = normalize(nEye);
  if (normal.z < 0.0) normal = -normal;
  float ambient    = lt.ambient;
  float direct     = lt.direct;
  float reflect    = lt.reflect;
  float spec_value = lt.spec;
  float shininess  = max(lt.shininess, 1.0);
  const float3 L0 = float3(0.0, 0.0, 1.0);
  const float3 L1 = normalize(float3(lt.klx, lt.kly, lt.klz));  // key light (cSetting_light)
  float intensity = ambient;
  float specular = 0.0;
  float n0 = dot(normal, L0);
  intensity += direct * saturate(n0);
  float n1 = dot(normal, L1);
  intensity += reflect * saturate(n1);
  if (n1 > 0.0) {
    float3 H1 = normalize(L1 + float3(0.0, 0.0, 1.0));
    // Soften specular on the lit-VBO path (cartoon + molecular surface). The
    // fixed view half-vector makes a flat ribbon/sheet face light up uniformly
    // bright, which reads as harsh triangular facets once screen-space AO is
    // removed from cartoons (#79) — the geometry normals are smooth, but the
    // strong highlight on each flat facet exaggerates the low-poly faces. Scale
    // the specular toward the softer ray-tracer/desktop look. Sphere/cylinder
    // impostors have their own shaders and keep full specular (sticks/spheres
    // stay glossy).
    const float kLitVBOSpecScale = 0.4;
    specular += kLitVBOSpecScale * spec_value * pow(max(dot(normal, H1), 0.0), shininess);
  }
  return baseColor * min(intensity, 1.0) + specular;
}

// Unlit input has no normal attribute (lines/ribbon/dots disable lighting).
struct VBOVertexInUnlit {
  float3 position [[attribute(0)]];
  float4 color    [[attribute(2)]];
};

struct VBOVertexOutUnlit {
  float4 position [[position]];
  float4 color;
  float  pointSize [[point_size]];
};

vertex VBOVertexOut vbo_vertex(
    VBOVertexIn in [[stage_in]],
    constant VBOUniforms& uniforms [[buffer(1)]])
{
  VBOVertexOut out;
  float4 eyePos = uniforms.modelview * float4(in.position, 1.0);
  out.position = uniforms.projection * eyePos;
  // PyMOL's projection is GL-convention (clip Z in [-1,1]); the impostors write
  // window depth explicitly as 0.5+0.5*ndcz ([0,1]). Remap the rasterized VBO z
  // the same way so VBO geometry (cartoon/surface) and impostors (spheres/
  // cylinders) share one depth convention and occlude correctly. (Also un-clips
  // the GL near half.)
  out.position.z = 0.5 * (out.position.z + out.position.w);
  // Carry the eye-space normal + raw color; lighting is done PER-FRAGMENT
  // (Phong) in vbo_fragment. Per-vertex (Gouraud) lighting baked the shade
  // into the interpolated color, which faceted the coarse cartoon mesh
  // (visible triangle banding); per-pixel shading is smooth.
  out.normalEye = (uniforms.modelview * float4(in.normal, 0.0)).xyz;
  out.color = in.color;
  // Eye-space distance from the camera (eyePos.z is negative in front), used by
  // the fragment stage to discard fragments outside this rep's clip planes.
  out.eyeDist = -eyePos.z;
  return out;
}

fragment float4 vbo_fragment(VBOVertexOut in [[stage_in]],
    constant LightU& lt [[buffer(0)]],
    constant ClipU& clip [[buffer(1)]])
{
  apply_rep_clip(clip, in.eyeDist);
  return float4(vbo_shade(in.color.rgb, in.normalEye, lt), in.color.a);
}

// Unlit: flat color, no lighting. Used for lines/ribbon (GL_LINES) and dots
// (GL_POINTS), which disable lighting and carry no normal attribute. Emits
// point_size so dots render at the requested screen size.
vertex VBOVertexOutUnlit vbo_vertex_unlit(
    VBOVertexInUnlit in [[stage_in]],
    constant VBOUniforms& uniforms [[buffer(1)]])
{
  VBOVertexOutUnlit out;
  out.position = uniforms.projection * (uniforms.modelview * float4(in.position, 1.0));
  out.position.z = 0.5 * (out.position.z + out.position.w);  // GL z [-1,1] -> Metal [0,1]
  out.color = in.color;
  out.pointSize = max(uniforms.pointSize, 1.0);
  return out;
}

fragment float4 vbo_fragment_unlit(VBOVertexOutUnlit in [[stage_in]])
{
  return in.color;
}

// Unlit, position-only: for uniform-colored line geometry whose VBO carries NO
// per-vertex color (alignment objects, distance/angle dashes). The regular
// unlit shader declares a color attribute [[attribute(2)]]; binding it against
// a VBO layout that has no color attribute makes pipeline creation FAIL ("Vertex
// attribute 2 is not defined in the vertex descriptor"), which on device floods
// every frame and gets the app killed. This variant takes only position and
// reads the color from a uniform (buffer 2), so the pipeline always builds.
struct VBOVertexInPosOnly {
  float3 position [[attribute(0)]];
};
vertex VBOVertexOutUnlit vbo_vertex_unlit_flat(
    VBOVertexInPosOnly in [[stage_in]],
    constant VBOUniforms& uniforms [[buffer(1)]],
    constant float4& flatColor [[buffer(2)]])
{
  VBOVertexOutUnlit out;
  out.position = uniforms.projection * (uniforms.modelview * float4(in.position, 1.0));
  out.position.z = 0.5 * (out.position.z + out.position.w);  // GL z [-1,1] -> Metal [0,1]
  out.color = flatColor;
  out.pointSize = max(uniforms.pointSize, 1.0);
  return out;
}

// --- Interior cap for a clipped closed SURFACE (stencil capping) ---
// MARK: render the surface's NOT-clipped faces position-only (cull none, no
// color/depth write, stencil = INVERT). Each face toggles stencil bit 0, so odd
// parity == the near (slab) plane is inside the solid the surface encloses.
struct CapMarkIn { float3 position [[attribute(0)]]; };
struct CapMarkOut { float4 position [[position]]; float eyeDist; };
vertex CapMarkOut cap_mark_vertex(CapMarkIn in [[stage_in]],
    constant VBOUniforms& u [[buffer(1)]]) {
  CapMarkOut o;
  float4 eye = u.modelview * float4(in.position, 1.0);
  float4 p = u.projection * eye;
  p.z = 0.5 * (p.z + p.w);   // GL [-1,1] -> Metal [0,1] (match the VBO)
  o.position = p;
  o.eyeDist = -eye.z;        // per-rep clip parity (matches vbo_vertex)
  return o;
}
// Discard the per-rep-clipped-away faces so the stencil parity marks the cut at
// the per-rep FRONT plane (not just the global near plane). The bound ClipU has
// back left open so only the viewer-facing cut is capped; when per-rep clip is
// off (enabled=0) this is a no-op and the global-slab/impostor cap is unchanged.
fragment float4 cap_mark_fragment(CapMarkOut in [[stage_in]],
    constant ClipU& clip [[buffer(1)]]) {
  apply_rep_clip(clip, in.eyeDist);
  return float4(0.0);  // color masked off (writeMask none)
}
// FILL: a full-screen triangle at the cut-plane depth (capDepth; the global near
// sentinel 1e-5 when per-rep clip is off); a stencil EQUAL test keeps only the
// interior pixels, filled with a flat interior color (and the pass clears the
// stencil bit so multiple capped surfaces don't interfere).
struct CapFillOut { float4 position [[position]]; };
vertex CapFillOut cap_fill_vertex(uint vid [[vertex_id]],
    constant float& capDepth [[buffer(0)]]) {
  float2 p = float2(float((vid << 1) & 2), float(vid & 2));   // (0,0)(2,0)(0,2)
  CapFillOut o;
  o.position = float4(p * 2.0 - 1.0, capDepth, 1.0);  // cover screen at cut depth
  return o;
}
fragment float4 cap_fill_fragment(constant float4& interiorColor [[buffer(0)]]) {
  // The cut cross-section faces the viewer, so light it with a +Z normal
  // through the shared two-light model (matches desktop's interior_normal
  // {0,0,1}). Caps use the default lighting (a minor flat fill).
  LightU lt = {0.14, 0.45, 0.481, 0.5, 55.0,
               0.34815531, 0.34815531, 0.87038828};  // default key light (+Z-ish)
  return float4(vbo_shade(interiorColor.rgb, float3(0.0, 0.0, 1.0), lt), interiorColor.a);
}

// --- Surface coverage (for the outer-contour outline) ---
// Rendered AFTER the scene into an R8 mask: position-only, applies the per-rep
// clip (so coverage matches the displayed/clipped surface), writes 1 for every
// surviving fragment. post_surface_contour then outlines the mask boundary.
struct CovOut { float4 position [[position]]; float eyeDist; };
vertex CovOut coverage_vertex(CapMarkIn in [[stage_in]],
    constant VBOUniforms& u [[buffer(1)]]) {
  CovOut o;
  float4 eye = u.modelview * float4(in.position, 1.0);
  o.position = u.projection * eye;
  o.position.z = 0.5 * (o.position.z + o.position.w);
  o.eyeDist = -eye.z;
  return o;
}
fragment float4 coverage_fragment(CovOut in [[stage_in]],
    constant ClipU& clip [[buffer(1)]]) {
  apply_rep_clip(clip, in.eyeDist);
  return float4(1.0, 0.0, 0.0, 1.0); // R8 mask = 1 where the surface covers
}

// Weighted-blended OIT output (McGuire/Bavoil). Transparent geometry writes
// premultiplied color*weight to the accum target and its alpha to the reveal
// target; the weight de-emphasizes far/low-alpha fragments. z = window depth.
struct OITFragOut {
  float4 accum  [[color(0)]];
  float  reveal [[color(1)]];
};
static float oit_weight(float a, float z) {
  return clamp(pow(min(1.0, a * 10.0) + 0.01, 3.0) * 1e8 *
               pow(1.0 - z * 0.9, 3.0), 1e-2, 3e3);
}
fragment OITFragOut vbo_fragment_oit(VBOVertexOut in [[stage_in]],
    constant LightU& lt [[buffer(0)]],
    constant ClipU& clip [[buffer(1)]])
{
  apply_rep_clip(clip, in.eyeDist);
  float4 c = float4(vbo_shade(in.color.rgb, in.normalEye, lt), in.color.a);
  float w = oit_weight(c.a, in.position.z);
  OITFragOut o;
  o.accum = float4(c.rgb * c.a, c.a) * w;
  o.reveal = c.a;
  return o;
}

// --- Anti-aliased screen-space line quads ("trilines"). The CPU expands each
// line segment into a feathered quad in clip space (see drawLinesAA); these
// shaders just pass the pre-projected verts through and feather-fade the edges
// by signed distance from the segment center line. ---
struct LineAAIn {
  float4 posClip    [[attribute(0)]];
  float4 color      [[attribute(1)]];
  float  centerDist [[attribute(2)]];
};
struct LineAAOut {
  float4 position [[position]];
  float4 color;
  float  centerDist;
};
struct LineAAU { float halfWidth; float feather; };
// CPU emits pre-projected clip-space quad verts; pass through with the standard
// GL->Metal z remap (matches vbo_vertex) so lines occlude consistently.
vertex LineAAOut line_aa_vertex(LineAAIn in [[stage_in]]) {
  LineAAOut o;
  o.position = in.posClip;
  o.position.z = 0.5 * (o.position.z + o.position.w);
  o.color = in.color;
  o.centerDist = in.centerDist;
  return o;
}
fragment float4 line_aa_fragment(LineAAOut in [[stage_in]],
    constant LineAAU& u [[buffer(0)]]) {
  float cov = 1.0 - clamp((abs(in.centerDist) - u.halfWidth) / max(u.feather, 1e-3), 0.0, 1.0);
  if (cov <= 0.0) discard_fragment();
  return float4(in.color.rgb, in.color.a * cov);
}
)";

  NSError* error = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:vboSrc options:nil error:&error];
  if (!lib) {
    NSLog(@"RendererMetal: failed to compile VBO shader: %@", error);
    return;
  }

  _vboVertexFunc = [lib newFunctionWithName:@"vbo_vertex"];
  _vboFragmentFunc = [lib newFunctionWithName:@"vbo_fragment"];
  _vboVertexUnlitFunc = [lib newFunctionWithName:@"vbo_vertex_unlit"];
  _vboFragmentUnlitFunc = [lib newFunctionWithName:@"vbo_fragment_unlit"];
  _vboVertexUnlitFlatFunc = [lib newFunctionWithName:@"vbo_vertex_unlit_flat"];
  _capMarkVtxFunc = [lib newFunctionWithName:@"cap_mark_vertex"];
  _capMarkFragFunc = [lib newFunctionWithName:@"cap_mark_fragment"];
  _capFillVtxFunc = [lib newFunctionWithName:@"cap_fill_vertex"];
  _capFillFragFunc = [lib newFunctionWithName:@"cap_fill_fragment"];
  _coverageVtxFunc = [lib newFunctionWithName:@"coverage_vertex"];
  _coverageFragFunc = [lib newFunctionWithName:@"coverage_fragment"];
  _lineAAVtxFunc = [lib newFunctionWithName:@"line_aa_vertex"];
  _lineAAFragFunc = [lib newFunctionWithName:@"line_aa_fragment"];
  if (!_vboVertexFunc || !_vboFragmentFunc) {
    NSLog(@"RendererMetal: VBO shader functions not found");
    [lib release];  // MRC: library (+1) consumed
    return;
  }

  // Surface interior-cap states + fill pipeline (mark pipeline is stride-keyed,
  // built lazily in drawVBOIndexed). MARK: stencil INVERT on every not-clipped
  // face (depth ALWAYS, no depth write). FILL: stencil EQUAL 1 -> write cap +
  // clear the bit (ZERO) so multiple capped surfaces don't accumulate.
  {
    MTLStencilDescriptor* sm = [[MTLStencilDescriptor alloc] init];
    sm.stencilCompareFunction = MTLCompareFunctionAlways;
    sm.depthStencilPassOperation = MTLStencilOperationInvert;
    sm.depthFailureOperation = MTLStencilOperationInvert;
    sm.stencilFailureOperation = MTLStencilOperationKeep;
    sm.readMask = 0x1;
    sm.writeMask = 0x1;
    MTLDepthStencilDescriptor* md = [[MTLDepthStencilDescriptor alloc] init];
    md.depthCompareFunction = MTLCompareFunctionAlways;
    md.depthWriteEnabled = NO;
    md.frontFaceStencil = sm;
    md.backFaceStencil = sm;
    _capMarkDSS = [_device newDepthStencilStateWithDescriptor:md];

    MTLStencilDescriptor* sf = [[MTLStencilDescriptor alloc] init];
    sf.stencilCompareFunction = MTLCompareFunctionEqual;     // ref = 1
    sf.depthStencilPassOperation = MTLStencilOperationZero;  // clear the bit
    sf.depthFailureOperation = MTLStencilOperationZero;
    sf.stencilFailureOperation = MTLStencilOperationKeep;
    sf.readMask = 0x1;
    sf.writeMask = 0x1;
    MTLDepthStencilDescriptor* fd = [[MTLDepthStencilDescriptor alloc] init];
    fd.depthCompareFunction = MTLCompareFunctionAlways;
    fd.depthWriteEnabled = YES;
    fd.frontFaceStencil = sf;
    fd.backFaceStencil = sf;
    _capFillDSS = [_device newDepthStencilStateWithDescriptor:fd];
    // MRC: the stencil/depth descriptors (+1) are consumed by the states above.
    [sm release];
    [md release];
    [sf release];
    [fd release];

    if (_capFillVtxFunc && _capFillFragFunc) {
      MTLRenderPipelineDescriptor* fp = [[MTLRenderPipelineDescriptor alloc] init];
      fp.vertexFunction = _capFillVtxFunc;
      fp.fragmentFunction = _capFillFragFunc;
      fp.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
      fp.rasterSampleCount = _sampleCount;
      fp.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
      fp.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
      NSError* ce = nil;
      _capFillPipeline = [_device newRenderPipelineStateWithDescriptor:fp error:&ce];
      if (!_capFillPipeline) NSLog(@"RendererMetal: cap fill pipeline failed: %@", ce);
      [fp release];  // MRC: descriptor (+1) consumed by pipeline creation
    }
  }

  // Create pipeline for UByte4Norm color format:
  // position: Float3 @ offset 0, normal: Float3 @ offset 12, color: UChar4Norm @ offset 24
  // stride = 28
  {
    MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
    vd.attributes[0].format = MTLVertexFormatFloat3;
    vd.attributes[0].offset = 0;
    vd.attributes[0].bufferIndex = 0;
    vd.attributes[1].format = MTLVertexFormatFloat3;
    vd.attributes[1].offset = 12;
    vd.attributes[1].bufferIndex = 0;
    vd.attributes[2].format = MTLVertexFormatUChar4Normalized;
    vd.attributes[2].offset = 24;
    vd.attributes[2].bufferIndex = 0;
    vd.layouts[0].stride = 28;
    vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

    MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
    psd.vertexFunction = _vboVertexFunc;
    psd.fragmentFunction = _vboFragmentFunc;
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

    _vboPipelineUByte = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
    if (!_vboPipelineUByte) {
      NSLog(@"RendererMetal: failed to create VBO UByte pipeline: %@", error);
    }
    [psd release];  // MRC: descriptor (+1) consumed by pipeline creation
  }

  // Create pipeline for Float4 color format:
  // position: Float3 @ offset 0, normal: Float3 @ offset 12, color: Float4 @ offset 24
  // stride = 40
  {
    MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
    vd.attributes[0].format = MTLVertexFormatFloat3;
    vd.attributes[0].offset = 0;
    vd.attributes[0].bufferIndex = 0;
    vd.attributes[1].format = MTLVertexFormatFloat3;
    vd.attributes[1].offset = 12;
    vd.attributes[1].bufferIndex = 0;
    vd.attributes[2].format = MTLVertexFormatFloat4;
    vd.attributes[2].offset = 24;
    vd.attributes[2].bufferIndex = 0;
    vd.layouts[0].stride = 40;
    vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

    MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
    psd.vertexFunction = _vboVertexFunc;
    psd.fragmentFunction = _vboFragmentFunc;
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

    _vboPipelineFloat = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
    if (!_vboPipelineFloat) {
      NSLog(@"RendererMetal: failed to create VBO Float pipeline: %@", error);
    }
    [psd release];  // MRC: descriptor (+1) consumed by pipeline creation
  }

  // Weighted-blended OIT variants: same vertex shader + geometry, but the
  // fragment writes to MRT (accum RGBA16F additive, reveal R16F multiplicative)
  // for order-independent transparency. Prebuild the two common layouts; other
  // layouts (e.g. the surface's stride-44 layout) get a one-off via
  // oitPipelineForVD() at draw time.
  _vboFragmentOitFunc = [lib newFunctionWithName:@"vbo_fragment_oit"];
  if (_vboFragmentOitFunc) {
    auto mkvd = [](MTLVertexFormat colorFmt, NSUInteger strideBytes) -> MTLVertexDescriptor* {
      MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
      vd.attributes[0].format = MTLVertexFormatFloat3;
      vd.attributes[0].offset = 0;
      vd.attributes[0].bufferIndex = 0;
      vd.attributes[1].format = MTLVertexFormatFloat3;
      vd.attributes[1].offset = 12;
      vd.attributes[1].bufferIndex = 0;
      vd.attributes[2].format = colorFmt;
      vd.attributes[2].offset = 24;
      vd.attributes[2].bufferIndex = 0;
      vd.layouts[0].stride = strideBytes;
      vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;
      return vd;
    };
    _vboOitPipelineUByte = oitPipelineForVD(mkvd(MTLVertexFormatUChar4Normalized, 28));
    _vboOitPipelineFloat = oitPipelineForVD(mkvd(MTLVertexFormatFloat4, 40));
  }
  [lib release];  // MRC: library (+1) no longer needed once functions are created
}

// Lazy, sample-count-aware. The CPU-expanded
// quad layout is: posClip Float4 @0, color Float4 @16, centerDist Float @32,
// stride 36, per-vertex. Src-alpha blending lets the feather alpha antialias.
void RendererMetal::buildLinePipeline() {
  if (_vboLinePipeline) return;
  if (!_lineAAVtxFunc || !_lineAAFragFunc) return;

  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  vd.attributes[0].format = MTLVertexFormatFloat4;
  vd.attributes[0].offset = 0;
  vd.attributes[0].bufferIndex = 0;
  vd.attributes[1].format = MTLVertexFormatFloat4;
  vd.attributes[1].offset = 16;
  vd.attributes[1].bufferIndex = 0;
  vd.attributes[2].format = MTLVertexFormatFloat;
  vd.attributes[2].offset = 32;
  vd.attributes[2].bufferIndex = 0;
  vd.layouts[0].stride = 36;
  vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

  MTLRenderPipelineDescriptor* p = [[MTLRenderPipelineDescriptor alloc] init];
  p.vertexFunction = _lineAAVtxFunc;
  p.fragmentFunction = _lineAAFragFunc;
  p.vertexDescriptor = vd;
  p.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
  p.colorAttachments[0].blendingEnabled = YES;
  p.colorAttachments[0].sourceRGBBlendFactor = MTLBlendFactorSourceAlpha;
  p.colorAttachments[0].destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
  p.colorAttachments[0].sourceAlphaBlendFactor = MTLBlendFactorOne;
  p.colorAttachments[0].destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
  p.rasterSampleCount = _sampleCount;
  p.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  p.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  NSError* e = nil;
  _vboLinePipeline = [_device newRenderPipelineStateWithDescriptor:p error:&e];
  if (!_vboLinePipeline) NSLog(@"RendererMetal: line AA pipeline failed: %@", e);
}

void RendererMetal::drawLinesAA(PrimitiveType mode, int vertexCount, const void* data,
                                size_t stride, int posOffset, int colorOffset, int colorType,
                                bool flat) {
  int nSeg = (mode == PrimitiveType::LineStrip) ? (vertexCount - 1) : (vertexCount / 2);
  if (nSeg <= 0) return;
  buildLinePipeline();
  if (!_vboLinePipeline) return;

  // mvp = projection * modelview (both column-major GL). Small column-major
  // mat4*vec4 helper (multiplyMatrices already gives mat4*mat4).
  const Mat4 mvp = multiplyMatrices(_projectionMatrix, _modelviewMatrix);
  auto mulMV4 = [&mvp](float x, float y, float z, float out[4]) {
    for (int r = 0; r < 4; ++r) {
      out[r] = mvp[0 * 4 + r] * x + mvp[1 * 4 + r] * y + mvp[2 * 4 + r] * z + mvp[3 * 4 + r];
    }
  };

  const float lw = std::max(_lineWidth, 1.0f);
  const float feather = 1.5f;
  const float halfExtent = lw * 0.5f + feather;

  const float vpW = (float)_viewport.width;
  const float vpH = (float)_viewport.height;
  if (vpW <= 0.0f || vpH <= 0.0f) return;
  const float halfVpX = 0.5f * vpW;
  const float halfVpY = 0.5f * vpH;

  const auto* base = static_cast<const unsigned char*>(data);
  auto readPos = [&](int idx, float p[3]) {
    const float* fp = reinterpret_cast<const float*>(base + (size_t)idx * stride + posOffset);
    p[0] = fp[0];
    p[1] = fp[1];
    p[2] = fp[2];
  };
  auto readColor = [&](int idx, float c[4]) {
    if (flat || colorOffset < 0) {
      c[0] = c[1] = c[2] = c[3] = 1.0f;
      return;
    }
    const unsigned char* cp = base + (size_t)idx * stride + colorOffset;
    if (colorType == 0) {  // 4x uchar normalized
      c[0] = cp[0] / 255.0f;
      c[1] = cp[1] / 255.0f;
      c[2] = cp[2] / 255.0f;
      c[3] = cp[3] / 255.0f;
    } else {  // 4x float
      const float* fc = reinterpret_cast<const float*>(cp);
      c[0] = fc[0];
      c[1] = fc[1];
      c[2] = fc[2];
      c[3] = fc[3];
    }
  };

  _lineExpand.clear();
  _lineExpand.reserve((size_t)nSeg * 6 * 9);

  // 2 triangles, 6 verts. corner endpoint index {0,0,1,1,0,1}, side {-1,1,-1,-1,1,1}.
  static const int kCornerEnd[6] = {0, 0, 1, 1, 0, 1};
  static const float kCornerSide[6] = {-1.f, 1.f, -1.f, -1.f, 1.f, 1.f};

  for (int s = 0; s < nSeg; ++s) {
    int i0, i1;
    if (mode == PrimitiveType::LineStrip) {
      i0 = s;
      i1 = s + 1;
    } else {
      i0 = 2 * s;
      i1 = 2 * s + 1;
    }

    float p0[3], p1[3], c0[4], c1[4];
    readPos(i0, p0);
    readPos(i1, p1);
    readColor(i0, c0);
    readColor(i1, c1);

    float clip0[4], clip1[4];
    mulMV4(p0[0], p0[1], p0[2], clip0);
    mulMV4(p1[0], p1[1], p1[2], clip1);
    if (clip0[3] == 0.0f || clip1[3] == 0.0f) continue;

    const float ndc0x = clip0[0] / clip0[3], ndc0y = clip0[1] / clip0[3];
    const float ndc1x = clip1[0] / clip1[3], ndc1y = clip1[1] / clip1[3];
    // Screen-space positions (px) for direction; only the direction matters.
    const float sx0 = ndc0x * halfVpX, sy0 = ndc0y * halfVpY;
    const float sx1 = ndc1x * halfVpX, sy1 = ndc1y * halfVpY;
    float dx = sx1 - sx0, dy = sy1 - sy0;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-4f) {
      dx = 1.0f;
      dy = 0.0f;
    } else {
      dx /= len;
      dy /= len;
    }
    const float px = -dy, py = dx;  // perpendicular

    for (int v = 0; v < 6; ++v) {
      const int end = kCornerEnd[v];
      const float side = kCornerSide[v];
      const float* clip = end ? clip1 : clip0;
      const float* col = end ? c1 : c0;
      const float ndcx = end ? ndc1x : ndc0x;
      const float ndcy = end ? ndc1y : ndc0y;
      // perp offset (across) + dir cap extension (along, ± to fill joins)
      const float capSign = end ? 1.0f : -1.0f;
      const float soX = px * side * halfExtent + dx * capSign * halfExtent;
      const float soY = py * side * halfExtent + dy * capSign * halfExtent;
      const float ndcOX = soX / halfVpX;
      const float ndcOY = soY / halfVpY;
      const float w = clip[3];
      // out clip = ((ndc + ndcOffset) * w, clipz, w)
      _lineExpand.push_back((ndcx + ndcOX) * w);
      _lineExpand.push_back((ndcy + ndcOY) * w);
      _lineExpand.push_back(clip[2]);
      _lineExpand.push_back(w);
      _lineExpand.push_back(col[0]);
      _lineExpand.push_back(col[1]);
      _lineExpand.push_back(col[2]);
      _lineExpand.push_back(col[3]);
      _lineExpand.push_back(side * halfExtent);  // centerDist
    }
  }

  const NSUInteger emitVerts = (NSUInteger)(_lineExpand.size() / 9);
  if (emitVerts == 0) return;
  const size_t bytes = _lineExpand.size() * sizeof(float);
  // One frame issues one line draw PER enabled object showing `lines`. The draw
  // is only RECORDED here; the GPU reads the vertex buffer when the command
  // buffer executes at commit. A single shared buffer rewritten at offset 0
  // therefore leaves every draw reading the LAST object's vertices, so only the
  // last object's lines appear and the others' silently vanish (#462; the same
  // trap endBatch documents for the immediate-mode batch VBO). Allocate a fresh
  // transient buffer per draw -- the encoder retains it until the frame
  // completes, so each draw reads its own vertices. This also removes the
  // cross-frame version of the hazard, where the next frame's memcpy landed on
  // a buffer the previous frame's GPU work could still be reading.
  id<MTLBuffer> lineVBO = [_device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
  if (!lineVBO) return;
  std::memcpy(lineVBO.contents, _lineExpand.data(), bytes);

  [_encoder setRenderPipelineState:_vboLinePipeline];
  // Lines depth-test/write like normal scene geometry.
  applyDepthStencilState();
  if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];
  [_encoder setVertexBuffer:lineVBO offset:0 atIndex:0];
  [lineVBO release];  // MRC: the encoder holds its own reference for the frame
  struct LineAAU {
    float halfWidth;
    float feather;
  } u = {lw * 0.5f, feather};
  [_encoder setFragmentBytes:&u length:sizeof(u) atIndex:0];
  [_encoder drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:emitVerts];
}

id<MTLRenderPipelineState> RendererMetal::oitPipelineForVD(MTLVertexDescriptor* vd) {
  if (!_vboVertexFunc || !_vboFragmentOitFunc) return nil;
  MTLRenderPipelineDescriptor* p = [[MTLRenderPipelineDescriptor alloc] init];
  p.vertexFunction = _vboVertexFunc;
  p.fragmentFunction = _vboFragmentOitFunc;
  p.vertexDescriptor = vd;
  p.colorAttachments[0].pixelFormat = MTLPixelFormatRGBA16Float;
  p.colorAttachments[0].blendingEnabled = YES;
  p.colorAttachments[0].sourceRGBBlendFactor = MTLBlendFactorOne;
  p.colorAttachments[0].destinationRGBBlendFactor = MTLBlendFactorOne;
  p.colorAttachments[0].sourceAlphaBlendFactor = MTLBlendFactorOne;
  p.colorAttachments[0].destinationAlphaBlendFactor = MTLBlendFactorOne;
  p.colorAttachments[1].pixelFormat = MTLPixelFormatR16Float;
  p.colorAttachments[1].blendingEnabled = YES;
  p.colorAttachments[1].sourceRGBBlendFactor = MTLBlendFactorZero;
  p.colorAttachments[1].destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceColor;
  p.colorAttachments[1].sourceAlphaBlendFactor = MTLBlendFactorZero;
  p.colorAttachments[1].destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceColor;
  p.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  p.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  NSError* e = nil;
  id<MTLRenderPipelineState> ps = [_device newRenderPipelineStateWithDescriptor:p error:&e];
  if (!ps) NSLog(@"RendererMetal: VBO OIT pipeline failed: %@", e);
  [p release];  // MRC: descriptor (+1) consumed by pipeline creation
  return ps;
}

// Build-once cache for one-off VBO pipelines (non-prebuilt vertex layouts, e.g.
// the molecular-surface stride-44). Returns a BORROWED pipeline owned by
// _vboPipelineCache; callers must NOT release it. Replaces the former per-draw
// oitPipelineForVD/inline-fallback creations, which leaked a
// +1 pipeline every frame (and paid full pipeline compilation each time).
id<MTLRenderPipelineState> RendererMetal::cachedVBOPipeline(VBOPipelineVariant variant,
                                                            size_t stride, int posOffset,
                                                            int normalOffset, int colorOffset,
                                                            int colorType,
                                                            MTLVertexDescriptor* vd) {
  // Stable FNV-1a key over the layout + variant + sample count.
  uint64_t key = 1469598103934665603ULL;
  auto mix = [&](uint64_t x) {
    key ^= x;
    key *= 1099511628211ULL;
  };
  mix((uint64_t)variant);
  mix((uint64_t)stride);
  mix((uint32_t)posOffset);
  mix((uint32_t)normalOffset);
  mix((uint32_t)colorOffset);
  mix((uint64_t)colorType);
  mix((uint64_t)_sampleCount);
  auto it = _vboPipelineCache.find(key);
  if (it != _vboPipelineCache.end()) return it->second;  // borrowed (cache-owned)

  id<MTLRenderPipelineState> ps = nil;
  if (variant == VBOPipelineVariant::Oit) {
    ps = oitPipelineForVD(vd);  // +1
  } else {
    // Opaque pass fallback for a non-prebuilt layout (mirrors the former inline
    // fallback): Lit uses the per-fragment shader, Unlit/UnlitFlat the flat ones.
    id<MTLFunction> vfn =
        (variant == VBOPipelineVariant::UnlitFlat)
            ? _vboVertexUnlitFlatFunc
            : (variant == VBOPipelineVariant::Unlit ? _vboVertexUnlitFunc : _vboVertexFunc);
    id<MTLFunction> ffn =
        (variant == VBOPipelineVariant::Lit) ? _vboFragmentFunc : _vboFragmentUnlitFunc;
    if (vfn && ffn) {
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
      NSError* error = nil;
      ps = [_device newRenderPipelineStateWithDescriptor:psd error:&error];
      if (!ps) NSLog(@"RendererMetal: fallback VBO pipeline failed: %@", error);
      [psd release];  // MRC: descriptor (+1) consumed by pipeline creation
    }
  }
  if (ps) _vboPipelineCache[key] = ps;  // cache takes ownership of the +1
  return ps;                            // borrowed
}
// Window-space depth ([0,1], matching the VBO's 0.5*(z+w) convention) of the
// per-rep front cut plane, for placing the interior cap-fill at the cut instead
// of the global near plane. Returns the near sentinel (1e-5) when per-rep clip
// is off (repClipFront < 0), preserving the original global-slab/impostor cap.
// P is the column-major projection matrix; clip.z/clip.w of eye point (0,0,ez)
// use P[10],P[14] and P[11],P[15].
static float repCapDepth(const float* P, float repClipFront) {
  if (repClipFront < 0.0f) return 1e-5f;
  float ez = -repClipFront;  // eye-space z of the front cut plane (negative)
  float cz = P[10] * ez + P[14];
  float cw = P[11] * ez + P[15];
  if (std::fabs(cw) < 1e-6f) return 1e-5f;
  return 0.5f * (cz / cw + 1.0f);  // GL NDC [-1,1] -> window [0,1]
}

// Reuses the cached Metal buffer for a CPU data pointer, or creates one.
id<MTLBuffer> RendererMetal::cachedVBOBuffer(const void* data, size_t size) {
  auto it = _vboCache.find(data);
  if (it != _vboCache.end()) return it->second;
  id<MTLBuffer> buf = [_device newBufferWithBytes:data
                                           length:size
                                          options:MTLResourceStorageModeShared];
  if (buf) _vboCache[data] = buf;
  return buf;
}

bool RendererMetal::bindVBODraw(id<MTLBuffer> vbo, size_t stride, int posOffset, int normalOffset,
                                int colorOffset, int colorType, bool unlit, VBOUniforms& u) {
  // Build a vertex descriptor matching the VBO layout
  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];

  // Position (attribute 0) — always Float3
  if (posOffset >= 0) {
    vd.attributes[0].format = MTLVertexFormatFloat3;
    vd.attributes[0].offset = static_cast<NSUInteger>(posOffset);
    vd.attributes[0].bufferIndex = 0;
  }

  // Normal (attribute 1) — always Float3
  if (normalOffset >= 0) {
    vd.attributes[1].format = MTLVertexFormatFloat3;
    vd.attributes[1].offset = static_cast<NSUInteger>(normalOffset);
    vd.attributes[1].bufferIndex = 0;
  }

  // Color (attribute 2)
  if (colorOffset >= 0) {
    vd.attributes[2].format =
        (colorType == 0) ? MTLVertexFormatUChar4Normalized : MTLVertexFormatFloat4;
    vd.attributes[2].offset = static_cast<NSUInteger>(colorOffset);
    vd.attributes[2].bufferIndex = 0;
  }

  vd.layouts[0].stride = static_cast<NSUInteger>(stride);
  vd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;

  // Uniform-colored geometry (no per-vertex color attribute): use the
  // position-only flat shader so the pipeline can be created (the normal unlit
  // shader needs attribute 2). Color comes from a uniform (buffer 2) below.
  bool flat = unlit && (colorOffset < 0);

  id<MTLRenderPipelineState> pipeline = nil;

  // Check if layout matches pre-built pipelines
  if (!unlit && posOffset == 0 && normalOffset == 12 && colorOffset == 24) {
    if (_oitActive) {
      // Transparent pass: route lit common layouts to the OIT MRT pipelines.
      if (colorType == 0 && stride == 28)
        pipeline = _vboOitPipelineUByte;
      else if (colorType == 1 && stride == 40)
        pipeline = _vboOitPipelineFloat;
    } else if (colorType == 0 && stride == 28 && _vboPipelineUByte) {
      pipeline = _vboPipelineUByte;
    } else if (colorType == 1 && stride == 40 && _vboPipelineFloat) {
      pipeline = _vboPipelineFloat;
    }
  }
  // In the OIT pass only MRT pipelines can render; build a one-off OIT
  // pipeline for non-standard lit layouts (e.g. surface stride 44). Unlit
  // geometry (lines/dots) is skipped in the transparent pass.
  if (_oitActive && !pipeline) {
    if (unlit) return false;
    pipeline = cachedVBOPipeline(VBOPipelineVariant::Oit, stride, posOffset, normalOffset,
                                 colorOffset, colorType, vd);
    if (!pipeline) return false;
  }

  // Fallback for a non-prebuilt opaque layout: build-once via the layout cache.
  // `flat` (no per-vertex color) → position-only shader; `unlit` (no normal:
  // lines/ribbon/dots) → flat-color shader; else the lit per-fragment shader.
  if (!pipeline) {
    VBOPipelineVariant v = flat ? VBOPipelineVariant::UnlitFlat
                                : (unlit ? VBOPipelineVariant::Unlit : VBOPipelineVariant::Lit);
    pipeline = cachedVBOPipeline(v, stride, posOffset, normalOffset, colorOffset, colorType, vd);
  }

  if (!pipeline) return false;

  [_encoder setRenderPipelineState:pipeline];

  // Apply depth/stencil
  if (_oitActive) {
    [_encoder setDepthStencilState:oitDepthState()];
  } else {
    applyDepthStencilState();
    if (_depthStencilState) {
      [_encoder setDepthStencilState:_depthStencilState];
    }
  }

  // Bind vertex buffer
  [_encoder setVertexBuffer:vbo offset:0 atIndex:0];

  // Upload modelview, projection matrices, and point size as uniforms
  std::memcpy(u.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(u.projection, _projectionMatrix.data(), 64);
  u.pointSize = _pointSize > 0.0f ? _pointSize : 1.0f;
  u._pad[0] = u._pad[1] = u._pad[2] = 0.0f;
  [_encoder setVertexBytes:&u length:sizeof(u) atIndex:1];
  // Lighting (Scene sliders) for the lit vbo_fragment / vbo_fragment_oit at
  // fragment buffer(0). Harmless for the unlit/flat pipelines (they don't read
  // it); without it lit geometry reads an unbound LightU and renders black.
  {
    struct {
      float a, d, r, s, sh, klx, kly, klz;
    } _lt = {_lightAmbient,   _lightDirect,    _lightReflect,   _lightSpecular,
             _lightShininess, _keyLightEye[0], _keyLightEye[1], _keyLightEye[2]};
    [_encoder setFragmentBytes:&_lt length:sizeof(_lt) atIndex:0];
  }
  // Per-rep clip planes for the lit vbo_fragment / vbo_fragment_oit at fragment
  // buffer(1). enabled=0 (front<0) leaves the rep clipped only by the global slab.
  {
    struct {
      float front, back, enabled, pad;
    } _cl = {_repClipFront, _repClipBack, _repClipFront >= 0.0f ? 1.0f : 0.0f, 0.0f};
    [_encoder setFragmentBytes:&_cl length:sizeof(_cl) atIndex:1];
  }

  // Flat (uniform-colored) geometry: supply the color the flat shader reads
  // from buffer 2. No per-vertex color is available here (the GL path would set
  // an a_Color uniform that the Metal backend drops), so render white — visible
  // and, crucially, non-crashing. TODO: thread the CGO's current color through.
  if (flat) {
    const float flatColor[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    [_encoder setVertexBytes:flatColor length:sizeof(flatColor) atIndex:2];
  }
  return true;
}

void RendererMetal::encodeVBODraw(PrimitiveType mode, int count, id<MTLBuffer> ibo) {
  if (ibo) {
    [_encoder drawIndexedPrimitives:toMTL(mode)
                         indexCount:static_cast<NSUInteger>(count)
                          indexType:MTLIndexTypeUInt32
                        indexBuffer:ibo
                  indexBufferOffset:0];
  } else {
    [_encoder drawPrimitives:toMTL(mode) vertexStart:0 vertexCount:static_cast<NSUInteger>(count)];
  }
}

// --- Surface interior cap (stencil) ---
// For a clipped closed surface with metal_interior_cap on: MARK the slab's
// interior cross-section (stencil parity over the not-clipped faces), then FILL
// it with a flat darkened color at the near-plane sentinel depth (which the
// post passes already exclude from SSAO). Opaque pass only.
void RendererMetal::drawVBOInteriorCap(PrimitiveType mode, int count, id<MTLBuffer> vbo,
                                       id<MTLBuffer> ibo, size_t stride, int posOffset,
                                       const VBOUniforms& u) {
  if (_oitActive || !_capFillPipeline || !_capMarkVtxFunc || !_capMarkFragFunc) return;
  if (!_capMarkPipeline || _capMarkStride != stride) {
    MTLVertexDescriptor* mvd = [[MTLVertexDescriptor alloc] init];
    mvd.attributes[0].format = MTLVertexFormatFloat3;
    mvd.attributes[0].offset = static_cast<NSUInteger>(posOffset < 0 ? 0 : posOffset);
    mvd.attributes[0].bufferIndex = 0;
    mvd.layouts[0].stride = static_cast<NSUInteger>(stride);
    mvd.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;
    MTLRenderPipelineDescriptor* mp = [[MTLRenderPipelineDescriptor alloc] init];
    mp.vertexFunction = _capMarkVtxFunc;
    mp.fragmentFunction = _capMarkFragFunc;
    mp.vertexDescriptor = mvd;
    mp.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    mp.colorAttachments[0].writeMask = MTLColorWriteMaskNone;  // stencil only
    mp.rasterSampleCount = _sampleCount;
    mp.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    mp.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    NSError* me = nil;
    [_capMarkPipeline release];  // MRC: release the previous-stride pipeline (+1)
    _capMarkPipeline = [_device newRenderPipelineStateWithDescriptor:mp error:&me];
    _capMarkStride = stride;
    if (!_capMarkPipeline) NSLog(@"RendererMetal: cap mark pipeline failed: %@", me);
    [mvd release];  // MRC: descriptors (+1) consumed by pipeline creation
    [mp release];
  }
  if (!_capMarkPipeline) return;
  [_encoder setRenderPipelineState:_capMarkPipeline];
  [_encoder setDepthStencilState:_capMarkDSS];
  [_encoder setStencilReferenceValue:1];
  [_encoder setCullMode:MTLCullModeNone];
  [_encoder setVertexBuffer:vbo offset:0 atIndex:0];
  [_encoder setVertexBytes:&u length:sizeof(u) atIndex:1];
  // MARK pass: per-rep clip (front-only) so stencil parity marks the cut at
  // the per-rep front plane; no-op (enabled=0) when per-rep clip is off.
  {
    struct {
      float front, back, enabled, pad;
    } _capcl = {_repClipFront, 1e6f, _repClipFront >= 0.0f ? 1.0f : 0.0f, 0.0f};
    [_encoder setFragmentBytes:&_capcl length:sizeof(_capcl) atIndex:1];
  }
  encodeVBODraw(mode, count, ibo);
  [_encoder setRenderPipelineState:_capFillPipeline];
  [_encoder setDepthStencilState:_capFillDSS];
  {
    float capDepth = repCapDepth(_projectionMatrix.data(), _repClipFront);
    [_encoder setVertexBytes:&capDepth length:sizeof(capDepth) atIndex:0];
  }
  const float interior[4] = {_capColorOverride ? _capColor[0] : 0.32f,
                             _capColorOverride ? _capColor[1] : 0.32f,
                             _capColorOverride ? _capColor[2] : 0.36f, 1.0f};  // darkened
  [_encoder setFragmentBytes:interior length:sizeof(interior) atIndex:0];
  [_encoder drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
  // Restore default state for subsequent geometry.
  [_encoder setCullMode:_cullFaceEnabled ? MTLCullModeBack : MTLCullModeNone];
  applyDepthStencilState();
  if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];
}

void RendererMetal::recordVBOCoverage(id<MTLBuffer> vbo, id<MTLBuffer> ibo, int count,
                                      size_t stride, int posOffset) {
  if (!_repContourEnabled && !_repAOExempt) return;
  CoverageDraw cd;
  cd.vbo = vbo;
  cd.ibo = ibo;
  cd.count = count;
  cd.stride = stride;
  cd.posOffset = posOffset;
  std::memcpy(cd.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(cd.projection, _projectionMatrix.data(), 64);
  cd.clipFront = _repClipFront;
  cd.clipBack = _repClipBack;
  // Surface outer-contour: rendered into the coverage mask after the scene
  // (see runPostChain).
  if (_repContourEnabled) {
    _coverageDraws.push_back(cd);
    _contourActive = true;
  }
  // Per-rep AO exemption (#79): cartoon/ribbon draws rasterized into the
  // AO-exempt mask after the scene (see renderAOExemptMask).
  if (_repAOExempt) _aoExemptDraws.push_back(cd);
}

void RendererMetal::drawVBO(PrimitiveType mode, int vertexCount, const void* data, size_t dataSize,
                            size_t stride, int posOffset, int normalOffset, int colorOffset,
                            int colorType, int interiorCap) {
  if (!data || dataSize == 0 || vertexCount <= 0) return;
  ensureEncoder();
  if (!_encoder) return;

  id<MTLBuffer> vbo = cachedVBOBuffer(data, dataSize);
  if (!vbo) return;

  // Anti-aliased lines: in the normal (non-OIT) pass, expand line segments on
  // the CPU into feathered screen-space quads instead of drawing 1px
  // MTLPrimitiveTypeLine. Lines are skipped in the transparent pass.
  if ((mode == PrimitiveType::Lines || mode == PrimitiveType::LineStrip) && !_oitActive) {
    bool flatLine = (colorOffset < 0);
    drawLinesAA(mode, vertexCount, data, stride, posOffset, colorOffset, colorType, flatLine);
    return;
  }

  // Lines/ribbon (no normal) and dots (Points) use the unlit flat-color +
  // point-size shader, even when their interleaved layout matches a prebuilt
  // (lit) pipeline.
  bool unlit = (normalOffset < 0) || (mode == PrimitiveType::Points);
  VBOUniforms u;
  if (!bindVBODraw(vbo, stride, posOffset, normalOffset, colorOffset, colorType, unlit, u)) return;
  encodeVBODraw(mode, vertexCount, nil);
  if (interiorCap) drawVBOInteriorCap(mode, vertexCount, vbo, nil, stride, posOffset, u);
  recordVBOCoverage(vbo, nil, vertexCount, stride, posOffset);
}

void RendererMetal::drawVBOIndexed(PrimitiveType mode, int indexCount, const void* vertexData,
                                   size_t vertexDataSize, size_t stride, int posOffset,
                                   int normalOffset, int colorOffset, int colorType,
                                   const void* indexData, size_t indexDataSize, int interiorCap) {
  if (!vertexData || !indexData || vertexDataSize == 0 || indexDataSize == 0 || indexCount <= 0)
    return;
  ensureEncoder();
  if (!_encoder) return;

  id<MTLBuffer> vbo = cachedVBOBuffer(vertexData, vertexDataSize);
  id<MTLBuffer> ibo = cachedVBOBuffer(indexData, indexDataSize);
  if (!vbo || !ibo) return;

  bool unlit = (normalOffset < 0);
  VBOUniforms u;
  if (!bindVBODraw(vbo, stride, posOffset, normalOffset, colorOffset, colorType, unlit, u)) return;
  encodeVBODraw(mode, indexCount, ibo);
  if (interiorCap) drawVBOInteriorCap(mode, indexCount, vbo, ibo, stride, posOffset, u);
  recordVBOCoverage(vbo, ibo, indexCount, stride, posOffset);
}

void RendererMetal::invalidateVBOCache(uint64_t key) {
  // Clear entire cache — key type changed to pointer-based.
  // MRC: the map holds +1 MTLBuffers (newBufferWithBytes); STL clear() does not
  // release them, so drop each ownership first. Any buffer still referenced by an
  // in-flight command buffer is kept alive by Metal's own retain until that
  // command buffer completes (same reasoning as ensurePostTargets).
  for (auto& kv : _vboCache) [kv.second release];
  _vboCache.clear();
}

void RendererMetal::invalidateVBOCacheEntry(const void* cpuData) {
  auto it = _vboCache.find(cpuData);
  if (it == _vboCache.end()) return;
  // MRC: the map holds a +1 MTLBuffer (newBufferWithBytes). Anything still
  // referenced by an in-flight command buffer stays alive via Metal's retain.
  [it->second release];
  _vboCache.erase(it);
}

}  // namespace pymol
