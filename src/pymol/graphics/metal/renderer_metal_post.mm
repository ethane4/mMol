
// Post-processing (SSAO, FXAA, outlines, fog) and order-independent transparency.

#include "pymol/graphics/metal/renderer_metal.h"
#include "pymol/image/my_png.h"

namespace pymol {

// Read a Metal BGRA8 texture and write it to a PNG via PyMOL's libpng writer.
// Converts BGRA→RGBA and flips to PyMOL's bottom-up row order. Runs off the
// render thread (command-buffer completion handler) after the GPU finishes.
static void rtWriteTextureToPNG(id<MTLTexture> tex, const std::string& path) {
  if (!tex) return;
  NSUInteger w = tex.width, h = tex.height;
  if (w == 0 || h == 0) return;
  NSUInteger bpr = w * 4;
  std::vector<uint8_t> bgra((size_t)bpr * h);
  [tex getBytes:bgra.data() bytesPerRow:bpr fromRegion:MTLRegionMake2D(0, 0, w, h) mipmapLevel:0];
  pymol::Image img((int)w, (int)h);
  uint8_t* dst = img.bits();
  for (NSUInteger y = 0; y < h; ++y) {
    const uint8_t* s = bgra.data() + (size_t)y * bpr;
    uint8_t* d = dst + (size_t)(h - 1 - y) * bpr;  // top-down -> bottom-up
    for (NSUInteger x = 0; x < w; ++x) {
      d[x * 4 + 0] = s[x * 4 + 2];  // R <- B
      d[x * 4 + 1] = s[x * 4 + 1];  // G
      d[x * 4 + 2] = s[x * 4 + 0];  // B <- R
      d[x * 4 + 3] = s[x * 4 + 3];  // A
    }
  }
  MyPNGWrite(path.c_str(), img, 0.0f, cMyPNG_FormatPNG, 1, 1.0f, 1.0f);
}
// Fullscreen-triangle vertex shader + post-process fragment shaders. A single
// library so all post pipelines share the vertex function.
static NSString* const kPostSrc = @R"(
using namespace metal;

// Linear eye distance (positive, toward the scene) from window depth [0,1].
// ortho>0.5: the projection is orthographic, so eye-z is LINEAR in ndc-z
//   (ez = (ndcz - projB)/projA with projA=proj[10], projB=proj[14]); the
//   perspective inverse (-projB/(ndcz+projA)) would divide by a near-zero /
//   wrong denominator and produce garbage distances.
static float post_linear_depth(float d, float projA, float projB,
                               float ortho = 0.0) {
  float ndcz = 2.0 * d - 1.0;
  float ez = (ortho > 0.5) ? ((ndcz - projB) / projA)  // ortho: linear
                           : (-projB / (ndcz + projA)); // persp: inverse
  return -ez;                         // distance from camera (positive)
}

// Screen-space crease occlusion: the SSAO ring test of post_ssao_fog. Samples
// N depths on a disk of radiusPx pixels around uv and counts neighbours that are 1..6% closer to the camera than the centre
// (zc, linear eye distance) — a depth step of that size is a crease / contact /
// occluding silhouette, which is exactly the contact darkening users read as
// depth. Returns the occluded fraction 0..1; callers scale it by an intensity.
// `cell` (x0,y0,x1,y1 in uv) limits the ring to the pixel's own grid_mode cell:
// a sample from the neighbouring cell is another object entirely and must not
// read as an occluder (#478). Default: the whole frame.
static bool post_in_cell(float2 uv, float4 cell) {
  return uv.x >= cell.x && uv.x < cell.z && uv.y >= cell.y && uv.y < cell.w;
}
static float post_ssao_occlusion(depth2d<float> depthTex, sampler s, float2 uv,
                                 float2 invres, float zc, float radiusPx,
                                 float projA, float projB, float ortho,
                                 float4 cell = float4(0.0, 0.0, 1.0, 1.0)) {
  const int N = 12;
  const float TWO_PI = 6.28318530718;
  const float range = 0.06; // ignore occluders farther than 6% of center z
  float occ = 0.0;
  for (int i = 0; i < N; i++) {
    float ang = (float(i) + 0.5) * (TWO_PI / float(N));
    // vary the radius across the ring to cover the disk
    float rr = radiusPx * (0.35 + 0.65 * float((i % 4) + 1) / 4.0);
    float2 off = float2(cos(ang), sin(ang)) * rr * invres;
    if (!post_in_cell(uv + off, cell)) continue; // other cell: not an occluder
    float dn = depthTex.sample(s, uv + off);
    if (dn >= 0.99999) continue; // background neighbor: no occlusion (no halo)
    float zn = post_linear_depth(dn, projA, projB, ortho);
    float diff = zc - zn; // > 0 when neighbor is closer to camera (occluder)
    if (diff > 0.0) {
      float rel = diff / max(zc, 1e-4);
      float w = smoothstep(0.0, 0.01, rel) *
                (1.0 - smoothstep(range * 0.5, range, rel));
      occ += w;
    }
  }
  return occ / float(N);
}

// Per-rep AO exemption (#79): aoMaskTex.r == 1 on front-most cartoon/ribbon
// pixels, where the crease term paints dark contour lines on ribbon silhouettes
// and self-folds. Dilated by two texels so grazing cartoon triangles at folds
// (whose single-sample re-raster depth can just miss the MSAA-resolved scene
// depth) are still covered and no AO line leaks through at the fold.
static bool post_ao_exempt(texture2d<float> aoMaskTex, sampler s, float2 uv,
                           float2 invres) {
  float m = 0.0;
  for (int dy = -2; dy <= 2; dy++)
    for (int dx = -2; dx <= 2; dx++)
      m = max(m, aoMaskTex.sample(s, uv + float2(float(dx), float(dy)) * invres).r);
  return m > 0.5;
}

struct PostVOut { float4 position [[position]]; float2 uv; };

vertex PostVOut post_vertex(uint vid [[vertex_id]]) {
  // Oversized triangle covering the screen: ids 0,1,2 -> (0,0),(2,0),(0,2).
  float2 p = float2((vid << 1) & 2, vid & 2);
  PostVOut o;
  o.position = float4(p * 2.0 - 1.0, 0.0, 1.0);
  o.uv = float2(p.x, 1.0 - p.y); // texture v=0 is top of screen
  return o;
}

fragment float4 post_blit(PostVOut in [[stage_in]],
    texture2d<float> src [[texture(0)]], sampler s [[sampler(0)]]) {
  return src.sample(s, in.uv);
}

// Offscreen transparent export: keep the fully post-processed RGB but rewrite
// alpha from scene depth so the background (depth at the far plane) becomes
// transparent and geometry stays opaque. Runs only for the offscreen PNG path
// when a transparent background was requested (ray_opaque_background 0); the
// live view never uses it. Alpha is binary (no matte fringing toward bg_rgb);
// exporting at 2x and downscaling anti-aliases the cutout.
fragment float4 post_export_alpha(PostVOut in [[stage_in]],
    texture2d<float> src [[texture(0)]],
    depth2d<float> depthTex [[texture(1)]],
    texture2d<float> revealTex [[texture(2)]], sampler s [[sampler(0)]]) {
  float d = depthTex.sample(s, in.uv);
  // Opaque geometry (depth < far) stays fully opaque; the far-plane background is
  // cut out. BUT transparent reps (the molecular surface) render via weighted-
  // blended OIT and do NOT write depth, so a depth-only cutout erases them from
  // the exported alpha wherever the translucent shell overhangs the background —
  // the surface vanishes from a transparent-background PNG while showing fine in
  // the viewport (issue #22). Recover their coverage from the OIT transmittance
  // (cover = 1 - reveal, matching oit_resolve) so the surface survives the matte.
  float cover = 1.0 - revealTex.sample(s, in.uv).r;
  float a = (d < 0.99995) ? 1.0 : cover;   // far plane == cleared background
  return float4(src.sample(s, in.uv).rgb, a);
}

// Weighted-blended OIT resolve: composite accumulated transparent color over
// the (already post-processed) opaque color. reveal = Π(1-α) = transmittance.
fragment float4 oit_resolve(PostVOut in [[stage_in]],
    texture2d<float> opaqueTex [[texture(0)]],
    texture2d<float> accumTex [[texture(1)]],
    texture2d<float> revealTex [[texture(2)]],
    sampler s [[sampler(0)]]) {
  float3 opaque = opaqueTex.sample(s, in.uv).rgb;
  float reveal = revealTex.sample(s, in.uv).r;
  float4 accum = accumTex.sample(s, in.uv);
  float3 avg = accum.rgb / max(accum.a, 1e-5);
  float cover = 1.0 - reveal;            // fraction covered by transparency
  return float4(avg * cover + opaque * (1.0 - cover), 1.0);
}

// FXAA (Timothy Lottes' classic console variant) — edge-aware blur that
// smooths ALL silhouettes, including the impostor ray-cast sphere/cylinder
// edges that geometry MSAA can't touch. Operates on luma of the composited
// scene color.
static float fxaa_luma(float3 c) { return dot(c, float3(0.299, 0.587, 0.114)); }

fragment float4 post_fxaa(PostVOut in [[stage_in]],
    texture2d<float> tex [[texture(0)]], sampler smp [[sampler(0)]]) {
  float2 inv = 1.0 / float2(tex.get_width(), tex.get_height());
  float3 rgbM = tex.sample(smp, in.uv).rgb;
  float lumaM  = fxaa_luma(rgbM);
  float lumaNW = fxaa_luma(tex.sample(smp, in.uv + float2(-1.0, -1.0) * inv).rgb);
  float lumaNE = fxaa_luma(tex.sample(smp, in.uv + float2( 1.0, -1.0) * inv).rgb);
  float lumaSW = fxaa_luma(tex.sample(smp, in.uv + float2(-1.0,  1.0) * inv).rgb);
  float lumaSE = fxaa_luma(tex.sample(smp, in.uv + float2( 1.0,  1.0) * inv).rgb);
  float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
  float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));
  float range = lumaMax - lumaMin;
  // Skip near-flat regions (no visible edge).
  if (range < max(0.0312, lumaMax * 0.125))
    return float4(rgbM, 1.0);

  float2 dir;
  dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
  dir.y =  ((lumaNW + lumaSW) - (lumaNE + lumaSE));
  float dirReduce =
      max((lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 * 0.125), 1.0 / 128.0);
  float rcpDirMin = 1.0 / (min(abs(dir.x), abs(dir.y)) + dirReduce);
  dir = clamp(dir * rcpDirMin, -8.0, 8.0) * inv;

  float3 rgbA = 0.5 * (tex.sample(smp, in.uv + dir * (1.0 / 3.0 - 0.5)).rgb +
                       tex.sample(smp, in.uv + dir * (2.0 / 3.0 - 0.5)).rgb);
  float3 rgbB = rgbA * 0.5 + 0.25 * (tex.sample(smp, in.uv + dir * -0.5).rgb +
                                     tex.sample(smp, in.uv + dir *  0.5).rgb);
  float lumaB = fxaa_luma(rgbB);
  if (lumaB < lumaMin || lumaB > lumaMax)
    return float4(rgbA, 1.0);
  return float4(rgbB, 1.0);
}

// SSAO (depth-difference / crease shading) + depth-cue fog, in one pass.
// Reads the scene color + scene depth, darkens creases where neighbors are
// modestly closer to the camera, then fades toward the background color with
// distance (PyMOL depth_cue). Scale-invariant: occlusion uses RELATIVE depth
// difference so it needs no per-scene world-unit tuning.
struct PostU {
  float projA;       // projection[10]
  float projB;       // projection[14]
  float fogStart;    // eye-space distance
  float fogEnd;
  float bgR, bgG, bgB;
  float fogEnabled;
  float aoEnabled;
  float aoIntensity;
  float aoRadiusPx;  // sample radius in pixels at render resolution
  float aoExemptEnabled; // >0.5: aoMaskTex marks pixels that skip the SSAO term (#79)
  float projOrtho;       // >0.5: orthographic projection (linear eye-z / no foreshortening)
  float pad0;
};

fragment float4 post_ssao_fog(PostVOut in [[stage_in]],
    texture2d<float> colorTex [[texture(0)]],
    depth2d<float> depthTex [[texture(1)]],
    texture2d<float> aoMaskTex [[texture(2)]],
    sampler s [[sampler(0)]],
    constant PostU& u [[buffer(0)]]) {
  float3 color = colorTex.sample(s, in.uv).rgb;
  float d = depthTex.sample(s, in.uv);

  // Interior-cap cross-sections are drawn at a sentinel depth right at the near
  // plane (~1e-5); nothing else sits that close (it would be clipped). Treat
  // those pixels as a flat solid fill: skip AO so the cut face stays
  // clean instead of being shaded by the screen-space post passes.
  bool flatCap = (d <= 0.0015);

  // Per-rep AO exemption (#79): aoMaskTex.r == 1 on front-most cartoon/ribbon
  // pixels. The depth-step SSAO term paints dark contour lines on their
  // silhouettes/self-folds/coil-crossings, so it is SKIPPED there (surface
  // pockets keep their AO). Dilate the mask by two texels so grazing cartoon
  // triangles at folds (whose single-sample re-raster depth can just miss the
  // MSAA-resolved scene depth) are still covered and no AO line leaks through.
  bool aoExempt = false;
  if (u.aoExemptEnabled > 0.5) {
    float2 iv = 1.0 / float2(colorTex.get_width(), colorTex.get_height());
    aoExempt = post_ao_exempt(aoMaskTex, s, in.uv, iv);
  }

  float ao = 1.0;
  if (u.aoEnabled > 0.5 && d < 0.99999 && !flatCap && !aoExempt) {
    float zc = post_linear_depth(d, u.projA, u.projB, u.projOrtho);
    float2 invres = 1.0 / float2(colorTex.get_width(), colorTex.get_height());
    float occ = post_ssao_occlusion(depthTex, s, in.uv, invres, zc, u.aoRadiusPx,
                                    u.projA, u.projB, u.projOrtho);
    ao = clamp(1.0 - occ * u.aoIntensity, 0.0, 1.0);
  }
  color *= ao;

  if (u.fogEnabled > 0.5 && d < 0.99999) {
    float dist = post_linear_depth(d, u.projA, u.projB, u.projOrtho);
    float fog = clamp((u.fogEnd - dist) / max(u.fogEnd - u.fogStart, 1e-4),
                      0.0, 1.0);
    color = mix(float3(u.bgR, u.bgG, u.bgB), color, fog);
  }
  return float4(color, 1.0);
}

// Silhouette / toon outlines: depth-based Sobel edge detection. Strong depth
// discontinuities (object silhouettes + where one part occludes another) get a
// dark contour. Scale-invariant: the gradient is taken relative to the center
// linear depth. Composited over the resolved scene color.
struct OutlineU {
  float projA, projB;     // projection[10]/[14] for linear depth
  float invW, invH;       // 1/resolution
  float colR, colG, colB; // outline color
  float thickness;        // sample step in pixels
};
fragment float4 post_outline(PostVOut in [[stage_in]],
    texture2d<float> colorTex [[texture(0)]],
    depth2d<float> depthTex [[texture(1)]],
    sampler s [[sampler(0)]],
    constant OutlineU& u [[buffer(0)]]) {
  float3 color = colorTex.sample(s, in.uv).rgb;
  float dc = depthTex.sample(s, in.uv);
  float2 px = float2(u.invW, u.invH) * u.thickness;
  // 3x3 linear-depth samples
  float zc = post_linear_depth(dc, u.projA, u.projB);
  float2 o = px;
  float zl = post_linear_depth(depthTex.sample(s, in.uv + float2(-o.x, 0)), u.projA, u.projB);
  float zr = post_linear_depth(depthTex.sample(s, in.uv + float2( o.x, 0)), u.projA, u.projB);
  float zt = post_linear_depth(depthTex.sample(s, in.uv + float2(0, -o.y)), u.projA, u.projB);
  float zb = post_linear_depth(depthTex.sample(s, in.uv + float2(0,  o.y)), u.projA, u.projB);
  // gradient magnitude relative to center depth (scale-invariant)
  float gx = abs(zr - zl);
  float gy = abs(zb - zt);
  float grad = (gx + gy) / max(zc, 1e-3);
  float edge = smoothstep(0.02, 0.08, grad);
  // don't outline pure background (both center far)
  edge *= step(dc, 0.99999);
  color = mix(color, float3(u.colR, u.colG, u.colB), edge);
  return float4(color, 1.0);
}

// Surface OUTER contour: outline the boundary of the surface coverage mask
// (covTex.r = 1 where the surface covers the screen, 0 elsewhere). A pixel is on
// the contour where coverage differs from a neighbor within the line width, so
// the line traces the outer silhouette (and any interior holes) of the surface —
// independent of the surface's transparency or the depth buffer. Composited over
// the scene color.
struct ContourU {
  float invW, invH;             // 1/resolution
  float thickness;              // ring radius in pixels (~half the line width)
  float colR, colG, colB, colA; // line color; colA folds in opaque/transparency
};
fragment float4 post_surface_contour(PostVOut in [[stage_in]],
    texture2d<float> colorTex [[texture(0)]],
    texture2d<float> covTex [[texture(1)]],
    sampler s [[sampler(0)]],
    constant ContourU& u [[buffer(0)]]) {
  float3 color = colorTex.sample(s, in.uv).rgb;
  float c = covTex.sample(s, in.uv).r;
  float2 px = float2(u.invW, u.invH) * max(u.thickness, 1.0);
  float e = 0.0;
  e = max(e, abs(c - covTex.sample(s, in.uv + float2( px.x, 0)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2(-px.x, 0)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2(0,  px.y)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2(0, -px.y)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2( px.x,  px.y)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2(-px.x,  px.y)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2( px.x, -px.y)).r));
  e = max(e, abs(c - covTex.sample(s, in.uv + float2(-px.x, -px.y)).r));
  float a = clamp(e, 0.0, 1.0) * u.colA;
  color = mix(color, float3(u.colR, u.colG, u.colB), a);
  return float4(color, 1.0);
}

)";

void RendererMetal::ensurePostTargets(NSUInteger w, NSUInteger h) {
  if (w == 0 || h == 0) return;
  buildPostPipelines();
  if (_sceneColor && _rtW == w && _rtH == h) return;
  _rtW = w;
  _rtH = h;

  // MRC (no ARC): release the PREVIOUS size's targets before reallocating.
  // Without this, every resize leaks a full target set, and a continuous
  // resize drag leaks dozens of full-res (MSAA) sets per second → the app is
  // OOM-jettisoned on a memory-constrained device ("crashes on resize").
  // The render-pass descriptors are updated below, so any in-flight command
  // buffer keeps its own retain until it completes — this is safe.
  [_sceneColor release];
  [_postColor release];
  [_sceneDepth release];
  [_sceneColorMS release];
  [_sceneDepthMS release];
  [_oitAccum release];
  [_oitReveal release];
  [_surfaceCoverageTex release];
  [_aoExemptMaskTex release];
  _sceneColor = _postColor = _sceneDepth = nil;
  _sceneColorMS = _sceneDepthMS = nil;
  _oitAccum = _oitReveal = nil;
  _surfaceCoverageTex = nil;
  _aoExemptMaskTex = nil;

  MTLTextureDescriptor* cd =
      [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatBGRA8Unorm
                                                         width:w
                                                        height:h
                                                     mipmapped:NO];
  cd.usage = MTLTextureUsageRenderTarget | MTLTextureUsageShaderRead;
  cd.storageMode = MTLStorageModePrivate;
  _sceneColor = [_device newTextureWithDescriptor:cd];
  _postColor = [_device newTextureWithDescriptor:cd];

  // Surface outer-contour coverage mask (single-sample R8): the surface footprint
  // is rendered here after the scene, then post_surface_contour outlines its
  // boundary. Single-sample so the post pass reads it directly (no MSAA resolve).
  MTLTextureDescriptor* covd =
      [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatR8Unorm
                                                         width:w
                                                        height:h
                                                     mipmapped:NO];
  covd.usage = MTLTextureUsageRenderTarget | MTLTextureUsageShaderRead;
  covd.storageMode = MTLStorageModePrivate;
  _surfaceCoverageTex = [_device newTextureWithDescriptor:covd];
  // Per-rep AO-exempt mask (#79): same single-sample R8, rasterized depth-tested
  // against _sceneDepth so only front-most cartoon/ribbon pixels are marked.
  _aoExemptMaskTex = [_device newTextureWithDescriptor:covd];

  MTLTextureDescriptor* dd =
      [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatDepth32Float_Stencil8
                                                         width:w
                                                        height:h
                                                     mipmapped:NO];
  dd.usage = MTLTextureUsageRenderTarget | MTLTextureUsageShaderRead;
  dd.storageMode = MTLStorageModePrivate;
  _sceneDepth = [_device newTextureWithDescriptor:dd];

  if (!_scenePassDesc) _scenePassDesc = [[MTLRenderPassDescriptor alloc] init];
  if (_sampleCount > 1) {
    // MSAA: the scene renders to multisampled targets and resolves into the
    // single-sample _sceneColor/_sceneDepth that the post chain samples.
    // StoreAndMultisampleResolve keeps the MS content so the scene encoder can
    // be resumed (selection indicators) across the OIT pass before the final
    // resolve.
    MTLTextureDescriptor* cms =
        [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatBGRA8Unorm
                                                           width:w
                                                          height:h
                                                       mipmapped:NO];
    cms.textureType = MTLTextureType2DMultisample;
    cms.sampleCount = _sampleCount;
    cms.usage = MTLTextureUsageRenderTarget;
    cms.storageMode = MTLStorageModePrivate;
    _sceneColorMS = [_device newTextureWithDescriptor:cms];
    MTLTextureDescriptor* dms =
        [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatDepth32Float_Stencil8
                                                           width:w
                                                          height:h
                                                       mipmapped:NO];
    dms.textureType = MTLTextureType2DMultisample;
    dms.sampleCount = _sampleCount;
    dms.usage = MTLTextureUsageRenderTarget;
    dms.storageMode = MTLStorageModePrivate;
    _sceneDepthMS = [_device newTextureWithDescriptor:dms];

    _scenePassDesc.colorAttachments[0].texture = _sceneColorMS;
    _scenePassDesc.colorAttachments[0].resolveTexture = _sceneColor;
    _scenePassDesc.colorAttachments[0].storeAction = MTLStoreActionStoreAndMultisampleResolve;
    _scenePassDesc.depthAttachment.texture = _sceneDepthMS;
    _scenePassDesc.depthAttachment.resolveTexture = _sceneDepth;
    _scenePassDesc.depthAttachment.depthResolveFilter = MTLMultisampleDepthResolveFilterSample0;
    _scenePassDesc.depthAttachment.storeAction = MTLStoreActionStoreAndMultisampleResolve;
    _scenePassDesc.stencilAttachment.texture = _sceneDepthMS;
    _scenePassDesc.stencilAttachment.storeAction = MTLStoreActionDontCare;
  } else {
    _sceneColorMS = nil;
    _sceneDepthMS = nil;
    _scenePassDesc.colorAttachments[0].texture = _sceneColor;
    _scenePassDesc.colorAttachments[0].resolveTexture = nil;
    _scenePassDesc.colorAttachments[0].storeAction = MTLStoreActionStore;
    _scenePassDesc.depthAttachment.texture = _sceneDepth;
    _scenePassDesc.depthAttachment.resolveTexture = nil;
    _scenePassDesc.depthAttachment.storeAction = MTLStoreActionStore;
    _scenePassDesc.stencilAttachment.texture = _sceneDepth;
  }

  // OIT accumulation targets (weighted-blended). accum holds Σ premul·weight,
  // reveal holds Π(1-α). They share the opaque depth (tested, not written).
  MTLTextureDescriptor* ad =
      [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatRGBA16Float
                                                         width:w
                                                        height:h
                                                     mipmapped:NO];
  ad.usage = MTLTextureUsageRenderTarget | MTLTextureUsageShaderRead;
  ad.storageMode = MTLStorageModePrivate;
  _oitAccum = [_device newTextureWithDescriptor:ad];
  MTLTextureDescriptor* rd =
      [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatR16Float
                                                         width:w
                                                        height:h
                                                     mipmapped:NO];
  rd.usage = MTLTextureUsageRenderTarget | MTLTextureUsageShaderRead;
  rd.storageMode = MTLStorageModePrivate;
  _oitReveal = [_device newTextureWithDescriptor:rd];

  if (!_oitPassDesc) _oitPassDesc = [[MTLRenderPassDescriptor alloc] init];
  _oitPassDesc.colorAttachments[0].texture = _oitAccum;
  _oitPassDesc.colorAttachments[0].loadAction = MTLLoadActionClear;
  _oitPassDesc.colorAttachments[0].clearColor = MTLClearColorMake(0, 0, 0, 0);
  _oitPassDesc.colorAttachments[0].storeAction = MTLStoreActionStore;
  _oitPassDesc.colorAttachments[1].texture = _oitReveal;
  _oitPassDesc.colorAttachments[1].loadAction = MTLLoadActionClear;
  _oitPassDesc.colorAttachments[1].clearColor = MTLClearColorMake(1, 0, 0, 0);
  _oitPassDesc.colorAttachments[1].storeAction = MTLStoreActionStore;
  _oitPassDesc.depthAttachment.texture = _sceneDepth;
  _oitPassDesc.depthAttachment.loadAction = MTLLoadActionLoad;  // keep opaque depth
  _oitPassDesc.depthAttachment.storeAction = MTLStoreActionStore;
  _oitPassDesc.stencilAttachment.texture = _sceneDepth;
  _oitPassDesc.stencilAttachment.loadAction = MTLLoadActionLoad;
}

void RendererMetal::buildPostPipelines() {
  if (_blitPipeline) return;
  NSError* err = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kPostSrc options:nil error:&err];
  if (!lib) {
    NSLog(@"RendererMetal: post lib compile failed: %@", err);
    return;
  }

  if (!_postSampler) {
    MTLSamplerDescriptor* sd = [[MTLSamplerDescriptor alloc] init];
    sd.minFilter = MTLSamplerMinMagFilterLinear;
    sd.magFilter = MTLSamplerMinMagFilterLinear;
    sd.sAddressMode = MTLSamplerAddressModeClampToEdge;
    sd.tAddressMode = MTLSamplerAddressModeClampToEdge;
    _postSampler = [_device newSamplerStateWithDescriptor:sd];
    [sd release];  // MRC: descriptor (+1) consumed by sampler creation
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"post_vertex"];
  auto mkpipe = [&](NSString* name) -> id<MTLRenderPipelineState> {
    id<MTLFunction> ffn = [lib newFunctionWithName:name];
    if (!vfn || !ffn) {
      [ffn release];
      return nil;
    }
    MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
    psd.vertexFunction = vfn;
    psd.fragmentFunction = ffn;
    psd.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    NSError* e = nil;
    id<MTLRenderPipelineState> p = [_device newRenderPipelineStateWithDescriptor:psd error:&e];
    if (!p) NSLog(@"RendererMetal: post pipeline %@ failed: %@", name, e);
    [psd release];  // MRC: descriptor (+1) consumed by pipeline creation
    [ffn release];  // MRC: fragment function (+1) retained by the pipeline
    return p;
  };
  _blitPipeline = mkpipe(@"post_blit");
  _exportAlphaPipeline = mkpipe(@"post_export_alpha");
  _fxaaPipeline = mkpipe(@"post_fxaa");
  _ssaoPipeline = mkpipe(@"post_ssao_fog");
  _oitResolvePipeline = mkpipe(@"oit_resolve");
  _outlinePipeline = mkpipe(@"post_outline");
  _surfaceContourPipeline = mkpipe(@"post_surface_contour");
  // MRC: the shared vertex function and the source library (+1) are no longer
  // needed once all post pipelines are built.
  [vfn release];
  [lib release];
}

void RendererMetal::setPostParams(int fogEnabled, float fogStart, float fogEnd, float bgR,
                                  float bgG, float bgB, int aoEnabled, int aaEnabled,
                                  int outlineEnabled, float projA, float projB, float outlineR,
                                  float outlineG, float outlineB, float outlineWidth, int ortho) {
  _outlineR = outlineR;
  _outlineG = outlineG;
  _outlineB = outlineB;
  _outlineWidth = outlineWidth;
  _postFogEnabled = fogEnabled;
  _fogStart = fogStart;
  _fogEnd = fogEnd;
  _bgR = bgR;
  _bgG = bgG;
  _bgB = bgB;
  _aoEnabled = aoEnabled;
  _aaEnabled = aaEnabled;
  _outlineEnabled = outlineEnabled;
  _projA = projA;
  _projB = projB;
  _projOrtho = ortho ? 1.0f : 0.0f;
}

void RendererMetal::setLightingParams(float ambient, float direct, float reflect, float specular,
                                      float shininess) {
  _lightAmbient = ambient;
  _lightDirect = direct;
  _lightReflect = reflect;
  _lightSpecular = specular;
  _lightShininess = shininess;
}

void RendererMetal::setKeyLightDir(const float* lightv) {
  // Direction TOWARD the light (used in dot(N,L)) is
  // -normalize(cSetting_light). PyMOL's default light yields the historical
  // normalize(0.4,0.4,1.0), so default behavior is unchanged.
  if (!lightv) return;
  float x = lightv[0], y = lightv[1], z = lightv[2];
  float len = std::sqrt(x * x + y * y + z * z);
  if (len < 1e-6f) {
    // Degenerate light vector: fall back to the historical key-light direction.
    _keyLightEye[0] = 0.34815531f;
    _keyLightEye[1] = 0.34815531f;
    _keyLightEye[2] = 0.87038828f;
    return;
  }
  _keyLightEye[0] = -x / len;
  _keyLightEye[1] = -y / len;
  _keyLightEye[2] = -z / len;
}

// Screen-space AO (metal_ssao) tuning: the darkening strength and the ring
// radius (a fraction of the render height).
static const float kSSAOIntensity = 0.8f;
static const float kSSAORadiusFrac = 0.015f;  // ~1.5% of height

void RendererMetal::runPostChain() {
  if (!_blitPipeline) return;

  static bool noAA = getenv("PYMOL_NO_AA") != nullptr;
  static bool noAO = getenv("PYMOL_NO_AO") != nullptr;
  bool doAO = _ssaoPipeline && _aoEnabled && !noAO;
  bool doFog = _ssaoPipeline && _postFogEnabled;
  id<MTLTexture> sceneSrc = _sceneColor;

  // Pass 1: SSAO + depth-cue/fog (color+depth -> post).
  if ((doAO || doFog) && _postColor) {
    // Per-rep AO exemption (#79): rasterize the stashed cartoon/ribbon draws
    // (depth-tested vs the scene) into _aoExemptMaskTex BEFORE the SSAO pass,
    // so post_ssao_fog can skip the contour term on those pixels. Only
    // meaningful when AO is actually running.
    bool aoMaskReady = doAO ? renderAOExemptMask() : false;

    struct {
      float projA, projB, fogStart, fogEnd;
      float bgR, bgG, bgB, fogEnabled;
      float aoEnabled, aoIntensity, aoRadiusPx, aoExemptEnabled;
      float projOrtho;  // matches MSL PostU: 1 = orthographic (#139)
      float pad0;
    } u;
    u.pad0 = 0.0f;
    u.projA = _projA;
    u.projB = _projB;
    u.fogStart = _fogStart;
    u.fogEnd = _fogEnd;
    u.bgR = _bgR;
    u.bgG = _bgG;
    u.bgB = _bgB;
    u.fogEnabled = doFog ? 1.0f : 0.0f;
    u.aoEnabled = doAO ? 1.0f : 0.0f;
    u.aoIntensity = kSSAOIntensity;
    u.aoRadiusPx = (float)_rtH * kSSAORadiusFrac;
    u.aoExemptEnabled = aoMaskReady ? 1.0f : 0.0f;
    u.projOrtho = _projOrtho;

    MTLRenderPassDescriptor* pd = [MTLRenderPassDescriptor renderPassDescriptor];
    pd.colorAttachments[0].texture = _postColor;
    pd.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    pd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> e1 = [_cmdBuffer renderCommandEncoderWithDescriptor:pd];
    [e1 setRenderPipelineState:_ssaoPipeline];
    [e1 setFragmentTexture:_sceneColor atIndex:0];
    [e1 setFragmentTexture:_sceneDepth atIndex:1];
    // texture(2) = AO-exempt mask; bind a valid 2D texture even when unused
    // (aoExemptEnabled gates the sample), since the shader declares it.
    [e1 setFragmentTexture:(aoMaskReady ? _aoExemptMaskTex : _sceneColor) atIndex:2];
    [e1 setFragmentSamplerState:_postSampler atIndex:0];
    [e1 setFragmentBytes:&u length:sizeof(u) atIndex:0];
    [e1 drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [e1 endEncoding];
    sceneSrc = _postColor;
  }

  // Pass 2: OIT resolve — composite accumulated transparency over the opaque
  // (post-processed) color. Ping-pongs to whichever target isn't the source.
  if (_oitHasContent && _oitResolvePipeline && _oitAccum) {
    id<MTLTexture> dst = (sceneSrc == _sceneColor) ? _postColor : _sceneColor;
    MTLRenderPassDescriptor* pd = [MTLRenderPassDescriptor renderPassDescriptor];
    pd.colorAttachments[0].texture = dst;
    pd.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    pd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> e2 = [_cmdBuffer renderCommandEncoderWithDescriptor:pd];
    [e2 setRenderPipelineState:_oitResolvePipeline];
    [e2 setFragmentTexture:sceneSrc atIndex:0];
    [e2 setFragmentTexture:_oitAccum atIndex:1];
    [e2 setFragmentTexture:_oitReveal atIndex:2];
    [e2 setFragmentSamplerState:_postSampler atIndex:0];
    [e2 drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [e2 endEncoding];
    sceneSrc = dst;
  }

  // Pass 3: silhouette/toon outlines (depth-based edges over the scene color).
  if (_outlineEnabled && _outlinePipeline && _sceneDepth) {
    id<MTLTexture> dst = (sceneSrc == _sceneColor) ? _postColor : _sceneColor;
    struct {
      float projA, projB, invW, invH;
      float colR, colG, colB, thickness;
    } u;
    u.projA = _projA;
    u.projB = _projB;
    u.invW = (_rtW > 0) ? 1.0f / (float)_rtW : 0.0f;
    u.invH = (_rtH > 0) ? 1.0f / (float)_rtH : 0.0f;
    u.colR = _outlineR;
    u.colG = _outlineG;
    u.colB = _outlineB;
    // Resolution-relative thickness so outlines aren't hairline-thin in hi-res
    // exports (pixelRadiusScale()==1 on the live view; see #48).
    u.thickness = _outlineWidth * pixelRadiusScale();
    MTLRenderPassDescriptor* pd = [MTLRenderPassDescriptor renderPassDescriptor];
    pd.colorAttachments[0].texture = dst;
    pd.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    pd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> e3 = [_cmdBuffer renderCommandEncoderWithDescriptor:pd];
    [e3 setRenderPipelineState:_outlinePipeline];
    [e3 setFragmentTexture:sceneSrc atIndex:0];
    [e3 setFragmentTexture:_sceneDepth atIndex:1];
    [e3 setFragmentSamplerState:_postSampler atIndex:0];
    [e3 setFragmentBytes:&u length:sizeof(u) atIndex:0];
    [e3 drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [e3 endEncoding];
    sceneSrc = dst;
  }

  // Pass 3.4: surface OUTER contour. Render the stashed surface footprint into a
  // coverage mask, then outline the mask boundary over the scene color. Works on
  // transparent/clipped surfaces (coverage is independent of depth/OIT).
  if (_contourActive && !_coverageDraws.empty() && _surfaceContourPipeline && _surfaceCoverageTex &&
      _coverageVtxFunc && _coverageFragFunc) {
    // (a) Coverage mask: rasterize the stashed surface draws (position-only,
    // per-rep clip applied) into the R8 mask, cleared to 0.
    MTLRenderPassDescriptor* cpd = [MTLRenderPassDescriptor renderPassDescriptor];
    cpd.colorAttachments[0].texture = _surfaceCoverageTex;
    cpd.colorAttachments[0].loadAction = MTLLoadActionClear;
    cpd.colorAttachments[0].clearColor = MTLClearColorMake(0, 0, 0, 0);
    cpd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> ce = [_cmdBuffer renderCommandEncoderWithDescriptor:cpd];
    for (const auto& cd : _coverageDraws) {
      if (!_coveragePipeline || _coverageStride != cd.stride) {
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
        pd2.rasterSampleCount = 1;
        NSError* pe = nil;
        [_coveragePipeline release];  // MRC: release the previous-stride pipeline
        _coveragePipeline = [_device newRenderPipelineStateWithDescriptor:pd2 error:&pe];
        _coverageStride = cd.stride;
        if (!_coveragePipeline) NSLog(@"RendererMetal: coverage pipeline failed: %@", pe);
        [vd release];
        [pd2 release];  // MRC: consumed by pipeline creation
      }
      if (!_coveragePipeline) continue;
      [ce setRenderPipelineState:_coveragePipeline];
      [ce setCullMode:MTLCullModeNone];
      [ce setVertexBuffer:cd.vbo offset:0 atIndex:0];
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
      [ce setVertexBytes:&cu length:sizeof(cu) atIndex:1];
      // Coverage mask is intentionally CLIP-AGNOSTIC for the per-rep clip: the
      // outer contour must trace the surface's true outer silhouette, not the
      // per-rep "peek inside" clip cross-section. Feeding cd.clipFront/Back here
      // discards the clipped fragments, so the clip cut edge would enter the
      // mask boundary and get outlined (the reported bug). Bind a disabled clip
      // (enabled=0) so the mask is the full surface footprint. The GLOBAL slab
      // still applies (it rides in cd.projection, shared with the live draw), and
      // the visible clipped surface + its interior cap are drawn by separate
      // passes that keep the per-rep clip — only the contour changes.
      {
        struct {
          float front, back, enabled, pad;
        } _cl = {-1.0f, 1e6f, 0.0f, 0.0f};
        [ce setFragmentBytes:&_cl length:sizeof(_cl) atIndex:1];
      }
      if (cd.ibo) {
        [ce drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                       indexCount:(NSUInteger)cd.count
                        indexType:MTLIndexTypeUInt32
                      indexBuffer:cd.ibo
                indexBufferOffset:0];
      } else {
        [ce drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:(NSUInteger)cd.count];
      }
    }
    [ce endEncoding];

    // (b) Outline the coverage boundary over the scene color.
    id<MTLTexture> dst = (sceneSrc == _sceneColor) ? _postColor : _sceneColor;
    struct {
      float invW, invH, thickness, colR, colG, colB, colA, _pad;
    } u;
    u.invW = (_rtW > 0) ? 1.0f / (float)_rtW : 0.0f;
    u.invH = (_rtH > 0) ? 1.0f / (float)_rtH : 0.0f;
    float halfW = _contourWidth * 0.5f;
    u.thickness = (halfW < 0.5f ? 0.5f : halfW) * pixelRadiusScale();
    u.colR = _contourColor[0];
    u.colG = _contourColor[1];
    u.colB = _contourColor[2];
    u.colA = _contourColor[3];
    u._pad = 0.0f;
    MTLRenderPassDescriptor* pd = [MTLRenderPassDescriptor renderPassDescriptor];
    pd.colorAttachments[0].texture = dst;
    pd.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    pd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> ec = [_cmdBuffer renderCommandEncoderWithDescriptor:pd];
    [ec setRenderPipelineState:_surfaceContourPipeline];
    [ec setFragmentTexture:sceneSrc atIndex:0];
    [ec setFragmentTexture:_surfaceCoverageTex atIndex:1];
    [ec setFragmentSamplerState:_postSampler atIndex:0];
    [ec setFragmentBytes:&u length:sizeof(u) atIndex:0];
    [ec drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [ec endEncoding];
    sceneSrc = dst;
  }

  // Transparent-export alpha (offscreen only). When a transparent background was
  // requested (ray_opaque_background 0 -> clear alpha 0), rewrite the captured
  // alpha from scene depth so the background is cut out while geometry stays
  // opaque. The post passes above all emit alpha=1, so this restores a usable
  // matte regardless of which of them ran. Live rendering (to a drawable) never
  // takes this path — its layer is opaque.
  if (_offscreen && _clearA < 0.5f && _exportAlphaPipeline && _sceneDepth) {
    id<MTLTexture> dst = (sceneSrc == _sceneColor) ? _postColor : _sceneColor;
    MTLRenderPassDescriptor* pd = [MTLRenderPassDescriptor renderPassDescriptor];
    pd.colorAttachments[0].texture = dst;
    pd.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    pd.colorAttachments[0].storeAction = MTLStoreActionStore;
    id<MTLRenderCommandEncoder> ea = [_cmdBuffer renderCommandEncoderWithDescriptor:pd];
    [ea setRenderPipelineState:_exportAlphaPipeline];
    [ea setFragmentTexture:sceneSrc atIndex:0];
    [ea setFragmentTexture:_sceneDepth atIndex:1];
    [ea setFragmentTexture:_oitReveal atIndex:2];  // recover transparent (surface) coverage
    [ea setFragmentSamplerState:_postSampler atIndex:0];
    [ea drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [ea endEncoding];
    sceneSrc = dst;
  }

  // Final pass → drawable: FXAA if available, else a 1:1 blit.
  // PYMOL_NO_AA forces the blit (A/B testing + user escape hatch).
  // Skipped entirely in offscreen mode — there is no drawable; the post chain
  // output in sceneSrc feeds the PNG capture below directly. FXAA still ran as
  // a ping-pong pass above if enabled, so the captured image is anti-aliased.
  if (!_offscreen) {
    _screenPassDesc.depthAttachment.texture = nil;
    _screenPassDesc.stencilAttachment.texture = nil;
    _screenPassDesc.colorAttachments[0].loadAction = MTLLoadActionDontCare;
    id<MTLRenderPipelineState> finalPipe =
        (_fxaaPipeline && _aaEnabled && !noAA) ? _fxaaPipeline : _blitPipeline;
    id<MTLRenderCommandEncoder> enc =
        [_cmdBuffer renderCommandEncoderWithDescriptor:_screenPassDesc];
    [enc setRenderPipelineState:finalPipe];
    [enc setFragmentTexture:sceneSrc atIndex:0];
    [enc setFragmentSamplerState:_postSampler atIndex:0];
    [enc drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
    [enc endEncoding];
  }

  // Pending rendered-frame PNG capture (png ray=0): copy the final
  // post-processed color into a CPU-readable texture and write the PNG once the
  // GPU completes. Captures at the current render resolution.
  if (!_capturePath.empty() && sceneSrc) {
    NSUInteger cw = sceneSrc.width, ch = sceneSrc.height;
    // Allocate a FRESH staging texture per capture (not a shared ivar): two
    // captures a frame or two apart would otherwise blit into the same texture
    // while the prior completion handler is still reading it off the GPU thread
    // — a data race producing a torn/garbled PNG.
    MTLTextureDescriptor* d =
        [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatBGRA8Unorm
                                                           width:cw
                                                          height:ch
                                                       mipmapped:NO];
    d.usage = MTLTextureUsageShaderRead;
    d.storageMode = MTLStorageModeShared;
    id<MTLTexture> capTex = [_device newTextureWithDescriptor:d];  // +1 (owned here)
    id<MTLBlitCommandEncoder> be = [_cmdBuffer blitCommandEncoder];
    [be copyFromTexture:sceneSrc
              sourceSlice:0
              sourceLevel:0
             sourceOrigin:MTLOriginMake(0, 0, 0)
               sourceSize:MTLSizeMake(cw, ch, 1)
                toTexture:capTex
         destinationSlice:0
         destinationLevel:0
        destinationOrigin:MTLOriginMake(0, 0, 0)];
    [be endEncoding];
    std::string path = _capturePath;
    _capturePath.clear();
    [_cmdBuffer addCompletedHandler:^(id<MTLCommandBuffer> cb) {
      rtWriteTextureToPNG(capTex, path);
    }];
    // MRC: the copied handler block retained capTex and releases it when the
    // block is destroyed (after it runs), so drop our +1 now — the texture stays
    // alive for exactly this capture's handler and is then freed once.
    [capTex release];
  }
}

id<MTLDepthStencilState> RendererMetal::oitDepthState() {
  // Transparent/OIT draws: depth-test vs opaque (LessEqual), never write depth.
  // Identical for every transparent draw, so build once and cache (MRC: +1,
  // released in the destructor).
  if (!_oitDepthState) {
    MTLDepthStencilDescriptor* dsd = [[MTLDepthStencilDescriptor alloc] init];
    dsd.depthCompareFunction = MTLCompareFunctionLessEqual;
    dsd.depthWriteEnabled = NO;
    _oitDepthState = [_device newDepthStencilStateWithDescriptor:dsd];
    [dsd release];
  }
  return _oitDepthState;
}

id<MTLDepthStencilState> RendererMetal::bezierDepthState() {
  // Opaque bezier-tube draws: standard LessEqual + depth write. Cached once
  // (MRC: +1, released in the destructor).
  if (!_bezierDepthState) {
    MTLDepthStencilDescriptor* dsd = [[MTLDepthStencilDescriptor alloc] init];
    dsd.depthCompareFunction = MTLCompareFunctionLessEqual;
    dsd.depthWriteEnabled = YES;
    _bezierDepthState = [_device newDepthStencilStateWithDescriptor:dsd];
    [dsd release];
  }
  return _bezierDepthState;
}

void RendererMetal::beginTransparentOIT() {
  // OIT requires its targets + pipelines; if any are missing, leave the scene
  // encoder active so transparent draws fall back to normal blending.
  if (!_cmdBuffer || !_oitPassDesc || !_vboOitPipelineUByte || !_oitResolvePipeline) return;
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }

  _passDesc = _oitPassDesc;
  _encoder = [_cmdBuffer renderCommandEncoderWithDescriptor:_oitPassDesc];
  if (!_encoder) {
    _passDesc = _scenePassDesc;
    return;
  }
  [_encoder setViewport:_viewport];
  // Depth-test against opaque depth (LEQUAL), but DO NOT write depth, so
  // transparent fragments occlude/are-occluded by opaque geometry yet never
  // hide each other.
  [_encoder setDepthStencilState:oitDepthState()];
  [_encoder setCullMode:MTLCullModeNone];
  _oitActive = true;
  _oitHasContent = true;
}

void RendererMetal::endTransparentOIT() {
  if (!_oitActive) return;
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }
  _oitActive = false;

  // Resume rendering into the scene color (LOAD, don't clear) so any
  // post-transparent draws (e.g. selection indicators) land on the opaque
  // image. The OIT resolve runs later in runPostChain.
  _passDesc = _scenePassDesc;
  _scenePassDesc.colorAttachments[0].loadAction = MTLLoadActionLoad;
  _scenePassDesc.depthAttachment.loadAction = MTLLoadActionLoad;
  _scenePassDesc.stencilAttachment.loadAction = MTLLoadActionLoad;
  _encoder = [_cmdBuffer renderCommandEncoderWithDescriptor:_scenePassDesc];
  if (_encoder) {
    [_encoder setViewport:_viewport];
    _depthTestEnabled = true;
    _depthWriteEnabled = true;
    _depthStencilDirty = true;
    applyDepthStencilState();
    [_encoder setCullMode:_cullFaceEnabled ? MTLCullModeBack : MTLCullModeNone];
  }
}

}  // namespace pymol
