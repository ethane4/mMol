
// Sphere, cylinder and Bezier-tube impostors.

#include "pymol/graphics/metal/renderer_metal.h"

namespace pymol {

// Inline MSL port of data/shaders/sphere.vs + sphere.fs. Per-pixel ray-sphere
// intersection with [[depth(any)]] output and diffuse + Blinn-Phong specular.
static NSString* const kSphereImpostorSrc = @R"(
using namespace metal;

struct SphereIn {
  float4 vertex_radius [[attribute(0)]];  // center xyz + radius
  float4 color         [[attribute(1)]];  // UByte4Norm -> float4
  float  rightUpFlags  [[attribute(2)]];
};
struct SphereU {
  float4x4 modelview;
  float4x4 projection;
  float sphere_size_scale;
  float ortho;          // 1 = orthographic
  float depthZeroToOne; // 1 if clip Z already [0,1]; else apply 0.5+0.5 remap
  float interiorCap;    // 1 = cap the slab cross-section with a solid interior color
  float4 interiorColor; // rgb cap color; .a > 0.5 => use it (else atom*0.45)
  // PyMOL lighting model (Scene sliders): ambient/direct/reflect/spec/shininess.
  float lAmbient, lDirect, lReflect, lSpecular, lShininess;
  float klx, kly, klz;  // key-light dir (toward light, eye space) = -normalize(cSetting_light)
};
struct SphereVOut {
  float4 position [[position]];
  float4 color;
  float3 sphere_center;  // eye space
  float  radius2;
  float3 point;          // eye-space impostor point
};
struct SphereFOut {
  float4 color [[color(0)]];
  float  depth [[depth(any)]];
};

static float2 outer_tangent_adjustment(float3 center, float radius_sq) {
  float2 xy_dist = float2(length(center.xz), length(center.yz));
  float2 cos_a = clamp(center.z / xy_dist, -1.0, 1.0);
  float2 cos_b = xy_dist / sqrt(radius_sq + (xy_dist * xy_dist));
  float2 cos_ab = cos_a * cos_b + sqrt((1.0 - cos_a*cos_a) * (1.0 - cos_b*cos_b));
  float2 cos_ab_sq = cos_ab * cos_ab;
  float2 tan_ab_sq = (1.0 - cos_ab_sq) / cos_ab_sq;
  return min(sqrt(tan_ab_sq + 1.0), 10.0);
}

vertex SphereVOut sphere_impostor_vertex(SphereIn in [[stage_in]],
    constant SphereU& u [[buffer(1)]]) {
  SphereVOut out;
  float radius = in.vertex_radius.w * u.sphere_size_scale;
  float3 mvcol0 = float3(u.modelview[0].x, u.modelview[0].y, u.modelview[0].z);
  radius /= max(length(mvcol0), 1e-6);
  float right = -1.0 + 2.0 * fmod(in.rightUpFlags, 2.0);
  float up    = -1.0 + 2.0 * floor(fmod(in.rightUpFlags / 2.0, 2.0));
  float4 tmppos = u.modelview * float4(in.vertex_radius.xyz, 1.0);
  out.color = in.color;
  out.radius2 = radius * radius;
  float2 corner = float2(right, up);
  if (u.ortho < 0.5)
    corner *= outer_tangent_adjustment(tmppos.xyz, out.radius2);
  float4 eyePos = tmppos;
  eyePos.xy += radius * corner;
  out.sphere_center = tmppos.xyz / tmppos.w;
  out.point = eyePos.xyz / eyePos.w;
  out.position = u.projection * eyePos;
  // GL->Metal clip-z remap (matches the VBO and the impostor fragment depth) so
  // the near half of the frustum (GL ndcz [-1,0)) isn't clipped — front atoms
  // stay visible instead of being sliced off by Metal's z>=0 clip volume.
  out.position.z = 0.5 * (out.position.z + out.position.w);
  return out;
}

// Weighted-blended OIT output for transparent spheres (accum + reveal + the
// ray-cast depth, so opaque geometry still occludes them).
struct SphereOITOut {
  float4 accum  [[color(0)]];
  float  reveal [[color(1)]];
  float  depth  [[depth(any)]];
};
static float sph_oit_weight(float a, float z) {
  return clamp(pow(min(1.0, a * 10.0) + 0.01, 3.0) * 1e8 *
               pow(1.0 - z * 0.9, 3.0), 1e-2, 3e3);
}

// Shared ray-sphere intersection + PyMOL two-light shading. Discards on miss /
// out-of-range depth. Returns lit rgb, alpha, and window depth.
static void sphere_shade(SphereVOut in, constant SphereU& u,
    thread float3& rgb, thread float& alpha, thread float& depth) {
  float3 ray_origin, ray_dir, sphere_dir;
  if (u.ortho >= 0.5) {
    ray_origin = in.point; ray_dir = float3(0.0,0.0,-1.0);
    sphere_dir = ray_origin - in.sphere_center;
  } else {
    ray_origin = float3(0.0); ray_dir = normalize(in.point);
    sphere_dir = in.sphere_center;
  }
  float b = dot(sphere_dir, ray_dir);
  float position = b*b + in.radius2 - dot(sphere_dir, sphere_dir);
  if (position < 0.0) discard_fragment();
  float nearest = b - sqrt(position);
  float3 ipoint = nearest * ray_dir + ray_origin;
  float4 clip = u.projection * float4(ipoint, 1.0);
  float ndcz = clip.z / clip.w;
  depth = (u.depthZeroToOne >= 0.5) ? ndcz : (0.5 + 0.5 * ndcz);

  if (depth >= 1.0) discard_fragment();   // behind the far plane
  if (depth <= 0.0) {
    // The near (slab) plane cuts in front of this surface point. If the BACK of
    // the sphere is still inside the slab, the plane slices through the sphere:
    // cap the cross-section with a flat interior color, else discard (=see-thru).
    if (u.interiorCap < 0.5) discard_fragment();
    float farthest = b + sqrt(position);
    float3 fpoint = farthest * ray_dir + ray_origin;
    float4 fclip = u.projection * float4(fpoint, 1.0);
    float fdepth = (u.depthZeroToOne >= 0.5) ? (fclip.z / fclip.w)
                                             : (0.5 + 0.5 * fclip.z / fclip.w);
    if (fdepth <= 0.0) discard_fragment();  // whole sphere in front of the slab
    depth = 1e-5;                           // flat disc at the slab plane (front)
    rgb = (u.interiorColor.a > 0.5) ? u.interiorColor.rgb  // ray_interior_color
                                    : in.color.rgb * 0.45; // else atom darkened
    alpha = in.color.a;
    return;
  }
  float3 normal = normalize(ipoint - in.sphere_center);

  // PyMOL two-light model, driven by the live ambient/direct/reflect/specular/
  // shininess settings (Scene sliders) uploaded in SphereU.
  float ambient    = u.lAmbient;
  float direct     = u.lDirect;
  float reflect    = u.lReflect;
  float spec_value = u.lSpecular;
  float shininess  = max(u.lShininess, 1.0);
  const float3 L0 = float3(0.0, 0.0, 1.0);
  const float3 L1 = normalize(float3(u.klx, u.kly, u.klz));  // key light (cSetting_light)
  float intensity = ambient;
  float specular = 0.0;
  float n0 = dot(normal, L0);
  intensity += direct * saturate(n0);
  float n1 = dot(normal, L1);
  intensity += reflect * saturate(n1);
  if (n1 > 0.0) {
    float3 H1 = normalize(L1 + float3(0.0, 0.0, 1.0));
    specular += spec_value * pow(max(dot(normal, H1), 0.0), shininess);
  }
  rgb = in.color.rgb * min(intensity, 1.0) + specular;
  alpha = in.color.a;
}

fragment SphereFOut sphere_impostor_fragment(SphereVOut in [[stage_in]],
    constant SphereU& u [[buffer(1)]]) {
  float3 rgb; float a; float depth;
  sphere_shade(in, u, rgb, a, depth);
  SphereFOut out;
  out.color = float4(rgb, a);
  out.depth = depth;
  return out;
}

fragment SphereOITOut sphere_impostor_fragment_oit(SphereVOut in [[stage_in]],
    constant SphereU& u [[buffer(1)]]) {
  float3 rgb; float a; float depth;
  sphere_shade(in, u, rgb, a, depth);
  float w = sph_oit_weight(a, depth);
  SphereOITOut out;
  out.accum = float4(rgb * a, a) * w;
  out.reveal = a;
  out.depth = depth;
  return out;
}

)";

void RendererMetal::buildImpostorPipelines() {
  if (_sphereImpostorPipeline) return;
  NSError* err = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kSphereImpostorSrc options:nil error:&err];
  if (!lib) {
    NSLog(@"RendererMetal: sphere impostor compile failed: %@", err);
    return;
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"sphere_impostor_vertex"];
  id<MTLFunction> ffn = [lib newFunctionWithName:@"sphere_impostor_fragment"];
  if (!vfn || !ffn) {
    NSLog(@"RendererMetal: sphere impostor funcs missing");
    return;
  }

  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  vd.attributes[0].format = MTLVertexFormatFloat4;  // a_vertex_radius
  vd.attributes[0].offset = 0;
  vd.attributes[0].bufferIndex = 0;
  vd.attributes[1].format = MTLVertexFormatUChar4Normalized;  // a_Color
  vd.attributes[1].offset = 16;
  vd.attributes[1].bufferIndex = 0;
  vd.attributes[2].format = MTLVertexFormatFloat;  // a_rightUpFlags
  vd.attributes[2].offset = 20;
  vd.attributes[2].bufferIndex = 0;
  vd.layouts[0].stride = 24;
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
  _sphereImpostorPipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&err];
  if (!_sphereImpostorPipeline) NSLog(@"RendererMetal: sphere impostor pipeline failed: %@", err);

  // Transparent sphere OIT variant: same vertex shader + geometry, MRT
  // accum/reveal output, ray-cast depth retained for occlusion.
  id<MTLFunction> offn = [lib newFunctionWithName:@"sphere_impostor_fragment_oit"];
  if (offn) {
    MTLRenderPipelineDescriptor* op = [[MTLRenderPipelineDescriptor alloc] init];
    op.vertexFunction = vfn;
    op.fragmentFunction = offn;
    op.vertexDescriptor = vd;
    op.colorAttachments[0].pixelFormat = MTLPixelFormatRGBA16Float;
    op.colorAttachments[0].blendingEnabled = YES;
    op.colorAttachments[0].sourceRGBBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].destinationRGBBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].sourceAlphaBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].destinationAlphaBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[1].pixelFormat = MTLPixelFormatR16Float;
    op.colorAttachments[1].blendingEnabled = YES;
    op.colorAttachments[1].sourceRGBBlendFactor = MTLBlendFactorZero;
    op.colorAttachments[1].destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceColor;
    op.colorAttachments[1].sourceAlphaBlendFactor = MTLBlendFactorZero;
    op.colorAttachments[1].destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceColor;
    op.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    op.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    _sphereOitPipeline = [_device newRenderPipelineStateWithDescriptor:op error:&err];
    if (!_sphereOitPipeline) NSLog(@"RendererMetal: sphere OIT pipeline failed: %@", err);
  }
}

void RendererMetal::drawSphereImpostors(const SphereImpostorDrawCall& call) {
  if (!call.data || call.dataSize == 0 || call.sphereCount <= 0) return;
  ensureEncoder();
  if (!_encoder) return;
  buildImpostorPipelines();
  if (!_sphereImpostorPipeline) return;
  if (_oitActive && !_sphereOitPipeline) return;  // no OIT variant: skip

  // Only the canonical packing (pos@0, color@16, rightUp@20 Float, stride 24)
  // is handled by the prebuilt pipeline. Log and bail otherwise (revisit if hit).
  bool common = (call.posRadiusOff == 0 && call.colorOff == 16 && call.rightUpOff == 20 &&
                 call.stride == 24 && call.rightUpIsFloat);
  if (!common) {
    NSLog(@"RendererMetal: sphere impostor non-standard packing pos=%d col=%d ru=%d(float=%d) "
          @"stride=%zu",
          call.posRadiusOff, call.colorOff, call.rightUpOff, call.rightUpIsFloat, call.stride);
    return;
  }

  id<MTLBuffer> vbo = nil;
  auto it = _vboCache.find(call.data);
  if (it != _vboCache.end())
    vbo = it->second;
  else {
    vbo = [_device newBufferWithBytes:call.data
                               length:call.dataSize
                              options:MTLResourceStorageModeShared];
    if (!vbo) return;
    _vboCache[call.data] = vbo;
  }

  // The impostor VBO is a plain triangle list (our build emits 6 verts/sphere
  // with corner flags {0,1,3,3,2,0}); draw all verts directly, no index buffer.
  NSUInteger vertexCount = call.stride ? (call.dataSize / call.stride) : 0;
  if (vertexCount < 3) return;

  if (_oitActive) {
    [_encoder setRenderPipelineState:_sphereOitPipeline];
    [_encoder setDepthStencilState:oitDepthState()];
  } else {
    [_encoder setRenderPipelineState:_sphereImpostorPipeline];
    applyDepthStencilState();
    if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];
  }
  [_encoder setVertexBuffer:vbo offset:0 atIndex:0];

  struct {
    float modelview[16];
    float projection[16];
    float sphere_size_scale;
    float ortho;
    float depthZeroToOne;
    float interiorCap;
    float interiorColor[4];
    float lAmbient, lDirect, lReflect, lSpecular, lShininess;
    float klx, kly, klz;  // key-light dir (toward light, eye space)
  } u;
  std::memcpy(u.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(u.projection, _projectionMatrix.data(), 64);
  u.lAmbient = _lightAmbient;
  u.lDirect = _lightDirect;
  u.lReflect = _lightReflect;
  u.lSpecular = _lightSpecular;
  u.lShininess = _lightShininess;
  u.klx = _keyLightEye[0];
  u.kly = _keyLightEye[1];
  u.klz = _keyLightEye[2];
  u.sphere_size_scale = call.sphereSizeScale;
  u.ortho = (float)call.ortho;
  // Projection is GL-convention ([-1,1] clip Z): remap to [0,1] for Metal's
  // fragment depth (matches the GL sphere.fs `0.5 + 0.5 * clipZ/clipW`).
  u.depthZeroToOne = 0.0f;
  // Cap the slab cross-section only in the opaque pass (not OIT).
  u.interiorCap = (call.interiorCap && !_oitActive) ? 1.0f : 0.0f;
  // ray_interior_color: .a>0.5 => use this rgb for the cap (else atom*0.45).
  u.interiorColor[0] = _capColor[0];
  u.interiorColor[1] = _capColor[1];
  u.interiorColor[2] = _capColor[2];
  u.interiorColor[3] = _capColorOverride ? 1.0f : 0.0f;
  [_encoder setVertexBytes:&u length:sizeof(u) atIndex:1];
  // The fragment shader also reads `u` (projection for depth, ortho flag) at
  // buffer index 1 — it must be bound to the fragment stage too, otherwise it
  // reads zero and clip.z/clip.w = 0/0 = NaN (all fragments fail the depth
  // range test / produce garbage depth).
  [_encoder setFragmentBytes:&u length:sizeof(u) atIndex:1];

  [_encoder drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:vertexCount];
}
// Inline MSL port of data/shaders/cylinder.vs + cylinder.fs. An 8-vertex box
// impostor (36 indices) with per-pixel ray-cylinder intersection, flat/round
// caps, two-color interpolation along the bond, [[depth(any)]] output, and the
// same PyMOL two-light shading as the sphere impostor. `a_cap` arrives EITHER as
// vertex attribute 6 (when the CGO baked the cap/interp bits per cylinder — the
// stick rep always does) OR, when the CGO emitted one constant for every
// cylinder, as the cap_const uniform; cap_const < 0 selects the per-vertex path.
// Both are live — see CylU::cap_const and GitHub issue #441.
static NSString* const kCylinderImpostorSrc = @R"(
using namespace metal;

struct CylIn {
  float3 vertex1 [[attribute(0)]];
  float3 vertex2 [[attribute(1)]];
  float4 color   [[attribute(2)]];
  float4 color2  [[attribute(3)]];
  float  radius  [[attribute(4)]];
  uchar  flags   [[attribute(5)]];
  uchar  cap     [[attribute(6)]];
};
struct CylU {
  float4x4 modelview;
  float4x4 projection;
  float uni_radius;
  float ortho;
  float depthZeroToOne;
  float no_flat_caps;
  float cap_const;      // >= 0: the CGO supplied one constant a_cap for every
                        // cylinder, use it. < 0: a_cap is per-vertex (attribute 6).
  float half_bond;
  float inv_height;
  float interiorCap;    // 1 = cap the slab cross-section with a solid interior color
  float4 interiorColor; // rgb cap color; .a > 0.5 => use it (else bond*0.45)
  // PyMOL lighting model (Scene sliders): ambient/direct/reflect/spec/shininess.
  float lAmbient, lDirect, lReflect, lSpecular, lShininess;
  float klx, kly, klz;  // key-light dir (toward light, eye space) = -normalize(cSetting_light)
};
struct CylVOut {
  float4 position [[position]];
  float3 surface_point;
  float3 axis;
  float3 base;
  float3 end_cyl;
  float3 U;
  float3 V;
  float radius;
  float inv_sqr_height;
  float cap;            // resolved a_cap bits: per-vertex when the CGO baked
                        // them per cylinder, else the CGO's constant (cap_const)
  float4 color1;
  float4 color2;
};
struct CylFOut { float4 color [[color(0)]]; float depth [[depth(any)]]; };

static float cyl_bit(thread float& bits) {
  float bit = fmod(bits, 2.0);
  bits = (bits - bit) / 2.0;
  return step(0.5, bit);
}

vertex CylVOut cyl_impostor_vertex(CylIn in [[stage_in]],
    constant CylU& u [[buffer(1)]]) {
  CylVOut o;
  // Normal matrix == upper-left 3x3 of the (rotation+uniform-scale) modelview.
  float3x3 N = float3x3(u.modelview[0].xyz, u.modelview[1].xyz, u.modelview[2].xyz);
  float uniformglscale = length(N[0]);
  float radius = (u.uni_radius != 0.0) ? (u.uni_radius * in.radius) : in.radius;
  o.color1 = in.color;
  o.color2 = in.color2;
  float3 attr_axis = in.vertex2 - in.vertex1;
  float ish = length(attr_axis) / uniformglscale;
  ish *= ish;
  o.inv_sqr_height = 1.0 / ish;
  float3 h = normalize(attr_axis);
  o.axis = normalize(N * h);
  float3 uu = cross(h, float3(1.0, 0.0, 0.0));
  if (dot(uu, uu) < 0.001) uu = cross(h, float3(0.0, 1.0, 0.0));
  uu = normalize(uu);
  float3 vv = normalize(cross(uu, h));
  o.U = normalize(N * uu);
  o.V = normalize(N * vv);
  float4 base4 = u.modelview * float4(in.vertex1, 1.0); o.base = base4.xyz;
  float4 end4  = u.modelview * float4(in.vertex2, 1.0); o.end_cyl = end4.xyz;

  float4 vert = float4(in.vertex1, 1.0);
  float packed = float(in.flags);
  float out_v   = cyl_bit(packed);
  float up_v    = cyl_bit(packed);
  float right_v = cyl_bit(packed);
  vert.xyz += up_v * attr_axis;
  vert.xyz += (2.0 * right_v - 1.0) * radius * uu;
  vert.xyz += (2.0 * out_v   - 1.0) * radius * vv;
  vert.xyz += (2.0 * up_v    - 1.0) * radius * h;

  float4 tvertex = u.modelview * vert;
  o.surface_point = tvertex.xyz;
  float4 pos = u.projection * tvertex;

  // Clamp z on the front clipping plane if the impostor box would be clipped
  // (we want to clip on the per-pixel depth, not the box face). See cylinder.vs.
  if (pos.z / pos.w < -1.0) {
    float diff = abs(base4.z - end4.z) + radius * 3.5;
    float4 inset = u.modelview * vert;
    inset.z -= diff;
    inset = u.projection * inset;
    if (inset.z / inset.w > -1.0) pos.z = -pos.w;
  }
  // GL->Metal clip-z remap, AFTER the GL-convention near clamp above: the GL
  // near plane (ndcz -1) maps to Metal z 0, so the near half of the box isn't
  // clipped away (front sticks stay visible). Matches the VBO + fragment depth.
  pos.z = 0.5 * (pos.z + pos.w);
  o.position = pos;
  o.radius = radius / uniformglscale;
  // Resolve a_cap here, once: CGOConvertShaderCylindersToCylinderShader emits it
  // as a constant only when every cylinder in the CGO shares one value, and bakes
  // it per-vertex otherwise (sticks, whose caps vary per bond). cap_const < 0
  // means "the VBO carries it" — see CylU::cap_const.
  o.cap = (u.cap_const >= 0.0) ? u.cap_const : float(in.cap);
  return o;
}

struct CylOITOut {
  float4 accum  [[color(0)]];
  float  reveal [[color(1)]];
  float  depth  [[depth(any)]];
};
static float cyl_oit_weight(float a, float z) {
  return clamp(pow(min(1.0, a * 10.0) + 0.01, 3.0) * 1e8 *
               pow(1.0 - z * 0.9, 3.0), 1e-2, 3e3);
}

// Shared ray-cylinder intersection (caps + two-color interp) + PyMOL shading.
// Discards on miss; returns lit rgb, alpha, window depth.
static void cyl_shade(CylVOut in, constant CylU& u,
    thread float3& rgb, thread float& alpha, thread float& depth) {
  float3 ray_target = in.surface_point;
  float3 ray_origin, ray_dir;
  if (u.ortho >= 0.5) { ray_origin = in.surface_point; ray_dir = float3(0.0,0.0,1.0); }
  else { ray_origin = float3(0.0); ray_dir = normalize(-ray_target); }

  float3x3 basis = float3x3(in.U, in.V, in.axis);
  float2 P = ((ray_target - in.base) * basis).xy;
  float2 D = (ray_dir * basis).xy;
  float radius2 = in.radius * in.radius;
  float a0 = P.x*P.x + P.y*P.y - radius2;
  float a1 = P.x*D.x + P.y*D.y;
  float a2 = D.x*D.x + D.y*D.y;
  float d = a1*a1 - a0*a2;
  if (d < 0.0) discard_fragment();
  float dist = (-a1 + sqrt(d)) / a2;
  float3 new_point = ray_target + dist * ray_dir;
  float3 tmp_point = new_point - in.base;
  float3 normal = normalize(tmp_point - in.axis * dot(tmp_point, in.axis));

  // cap bits: 0 frontcap, 1 endcap, 2 frontcapround, 3 endcapround, 4 interp
  float fcap = in.cap + 0.001;
  bool frontcap      = cyl_bit(fcap) > 0.5;
  bool endcap        = cyl_bit(fcap) > 0.5;
  bool frontcapround = (cyl_bit(fcap) > 0.5) && (u.no_flat_caps > 0.5);
  bool endcapround   = (cyl_bit(fcap) > 0.5) && (u.no_flat_caps > 0.5);
  bool nocolorinterp = !(cyl_bit(fcap) > 0.5);

  float ratio = dot(new_point - in.base, in.end_cyl - in.base) * in.inv_sqr_height;
  if (nocolorinterp) {
    float dp = clamp(-u.half_bond * new_point.z * u.inv_height, 0.0, 0.5);
    ratio = smoothstep(0.5 - dp, 0.5 + dp, ratio);
  } else {
    ratio = clamp(ratio, 0.0, 1.0);
  }
  float4 color = mix(in.color1, in.color2, ratio);

  bool cap_test_base = 0.0 > dot(new_point - in.base, in.axis);
  bool cap_test_end  = 0.0 < dot(new_point - in.end_cyl, in.axis);
  if (cap_test_base || cap_test_end) {
    float3 thisaxis = -in.axis;
    float3 thisbase = in.base;
    if (cap_test_end) {
      thisaxis = in.axis; thisbase = in.end_cyl;
      frontcap = endcap; frontcapround = endcapround;
    }
    if (!frontcap) discard_fragment();
    if (frontcapround) {
      float3 sd = thisbase - ray_origin;
      float b = dot(sd, ray_dir);
      float pos = b*b + radius2 - dot(sd, sd);
      if (pos < 0.0) discard_fragment();
      float nr = sqrt(pos) + b;
      new_point = nr * ray_dir + ray_origin;
      normal = normalize(new_point - thisbase);
    } else {
      float dNV = dot(thisaxis, ray_dir);
      if (dNV < 0.0) discard_fragment();
      float nr = dot(thisaxis, thisbase - ray_origin) / dNV;
      new_point = ray_dir * nr + ray_origin;
      if (dot(new_point - thisbase, new_point - thisbase) > radius2) discard_fragment();
      normal = thisaxis;
    }
  }

  float4 clip = u.projection * float4(new_point, 1.0);
  float ndcz = clip.z / clip.w;
  depth = (u.depthZeroToOne >= 0.5) ? ndcz : (0.5 + 0.5 * ndcz);
  if (depth <= 0.0) {
    // Surface point is in front of the slab. If the cylinder body still straddles
    // the slab (its far wall is inside), cap the cross-section with a flat
    // interior color instead of discarding (= see-through). Else discard.
    if (u.interiorCap < 0.5) discard_fragment();
    float dist2 = (-a1 - sqrt(d)) / a2;            // far body intersection
    float3 p2 = ray_target + dist2 * ray_dir;
    float4 c2 = u.projection * float4(p2, 1.0);
    float d2 = (u.depthZeroToOne >= 0.5) ? (c2.z / c2.w) : (0.5 + 0.5 * c2.z / c2.w);
    if (d2 <= 0.0) discard_fragment();             // whole stick in front of slab
    depth = 1e-5;                                  // flat disc at the slab plane
    rgb = (u.interiorColor.a > 0.5) ? u.interiorColor.rgb  // ray_interior_color
                                    : color.rgb * 0.45;    // else bond darkened
    alpha = color.a;
    return;
  }

  // PyMOL two-light model, driven by the live lighting settings in CylU.
  float ambient    = u.lAmbient;
  float direct     = u.lDirect;
  float reflect    = u.lReflect;
  float spec_value = u.lSpecular;
  float shininess  = max(u.lShininess, 1.0);
  const float3 L0 = float3(0.0, 0.0, 1.0);
  const float3 L1 = normalize(float3(u.klx, u.kly, u.klz));  // key light (cSetting_light)
  float intensity = ambient;
  float specular = 0.0;
  float n0 = dot(normal, L0);
  intensity += direct * saturate(n0);
  float n1 = dot(normal, L1);
  intensity += reflect * saturate(n1);
  if (n1 > 0.0) {
    float3 H1 = normalize(L1 + float3(0.0, 0.0, 1.0));
    specular += spec_value * pow(max(dot(normal, H1), 0.0), shininess);
  }
  rgb = color.rgb * min(intensity, 1.0) + specular;
  alpha = color.a;
}

fragment CylFOut cyl_impostor_fragment(CylVOut in [[stage_in]],
    constant CylU& u [[buffer(1)]]) {
  float3 rgb; float a; float depth;
  cyl_shade(in, u, rgb, a, depth);
  CylFOut o;
  o.color = float4(rgb, a);
  o.depth = depth;
  return o;
}

fragment CylOITOut cyl_impostor_fragment_oit(CylVOut in [[stage_in]],
    constant CylU& u [[buffer(1)]]) {
  float3 rgb; float a; float depth;
  cyl_shade(in, u, rgb, a, depth);
  float w = cyl_oit_weight(a, depth);
  CylOITOut o;
  o.accum = float4(rgb * a, a) * w;
  o.reveal = a;
  o.depth = depth;
  return o;
}

)";

void RendererMetal::buildCylinderImpostorPipeline(const CylinderImpostorDrawCall& call) {
  // Cache per vertex layout. a_cap's offset is part of the descriptor, so a
  // stick VBO (per-vertex a_cap) and a CGO VBO (constant a_cap) need different
  // pipelines even at the same stride — and Move mode draws both every frame, so
  // a single slot would recompile the MSL library twice per frame and (MRC) leak
  // the displaced pipelines. Point the ivars at this layout's entry instead.
  const auto layout = std::make_pair(static_cast<NSUInteger>(call.stride), call.capOff);
  {
    auto it = _cylinderPipelines.find(layout);
    if (it != _cylinderPipelines.end()) {
      _cylinderImpostorPipeline = it->second.opaque;
      _cylinderOitPipeline = it->second.oit;
      return;
    }
  }
  _cylinderImpostorPipeline = nil;
  _cylinderOitPipeline = nil;

  NSError* err = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kCylinderImpostorSrc options:nil error:&err];
  if (!lib) {
    NSLog(@"RendererMetal: cyl impostor compile failed: %@", err);
    return;
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"cyl_impostor_vertex"];
  id<MTLFunction> ffn = [lib newFunctionWithName:@"cyl_impostor_fragment"];
  if (!vfn || !ffn) {
    NSLog(@"RendererMetal: cyl impostor funcs missing");
    return;
  }

  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  vd.attributes[0].format = MTLVertexFormatFloat3;  // attr_vertex1
  vd.attributes[0].offset = call.v1Off;
  vd.attributes[0].bufferIndex = 0;
  vd.attributes[1].format = MTLVertexFormatFloat3;  // attr_vertex2
  vd.attributes[1].offset = call.v2Off;
  vd.attributes[1].bufferIndex = 0;
  MTLVertexFormat colFmt =
      call.colorIsFloat ? MTLVertexFormatFloat4 : MTLVertexFormatUChar4Normalized;
  vd.attributes[2].format = colFmt;  // a_Color
  vd.attributes[2].offset = call.colorOff;
  vd.attributes[2].bufferIndex = 0;
  vd.attributes[3].format = colFmt;  // a_Color2
  vd.attributes[3].offset = call.color2Off;
  vd.attributes[3].bufferIndex = 0;
  vd.attributes[4].format = MTLVertexFormatFloat;  // attr_radius
  vd.attributes[4].offset = call.radiusOff;
  vd.attributes[4].bufferIndex = 0;
  vd.attributes[5].format = MTLVertexFormatUChar;  // attr_flags (UByte)
  vd.attributes[5].offset = call.flagsOff;
  vd.attributes[5].bufferIndex = 0;
  // a_cap (UByte). Metal requires every attribute the MSL declares to exist in
  // the descriptor, so when the CGO supplied a CONSTANT a_cap (capOff < 0) we
  // still bind slot 6 — to attr_flags, a byte we know is in range — and the
  // shader ignores it because cap_const >= 0. See CylU::cap_const.
  vd.attributes[6].format = MTLVertexFormatUChar;  // a_cap (UByte)
  vd.attributes[6].offset = (call.capOff >= 0) ? call.capOff : call.flagsOff;
  vd.attributes[6].bufferIndex = 0;
  vd.layouts[0].stride = call.stride;
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
  _cylinderImpostorPipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&err];
  if (!_cylinderImpostorPipeline) NSLog(@"RendererMetal: cyl impostor pipeline failed: %@", err);

  // Transparent cylinder OIT variant (MRT accum/reveal, ray-cast depth kept).
  id<MTLFunction> offn = [lib newFunctionWithName:@"cyl_impostor_fragment_oit"];
  if (offn) {
    MTLRenderPipelineDescriptor* op = [[MTLRenderPipelineDescriptor alloc] init];
    op.vertexFunction = vfn;
    op.fragmentFunction = offn;
    op.vertexDescriptor = vd;
    op.colorAttachments[0].pixelFormat = MTLPixelFormatRGBA16Float;
    op.colorAttachments[0].blendingEnabled = YES;
    op.colorAttachments[0].sourceRGBBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].destinationRGBBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].sourceAlphaBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[0].destinationAlphaBlendFactor = MTLBlendFactorOne;
    op.colorAttachments[1].pixelFormat = MTLPixelFormatR16Float;
    op.colorAttachments[1].blendingEnabled = YES;
    op.colorAttachments[1].sourceRGBBlendFactor = MTLBlendFactorZero;
    op.colorAttachments[1].destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceColor;
    op.colorAttachments[1].sourceAlphaBlendFactor = MTLBlendFactorZero;
    op.colorAttachments[1].destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceColor;
    op.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    op.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
    _cylinderOitPipeline = [_device newRenderPipelineStateWithDescriptor:op error:&err];
    if (!_cylinderOitPipeline) NSLog(@"RendererMetal: cyl OIT pipeline failed: %@", err);
  }

  // The map takes ownership of the +1 pipelines; the ivars stay as aliases.
  // Only cache a layout whose opaque pipeline compiled, so a transient failure
  // is retried rather than cached forever.
  if (_cylinderImpostorPipeline) {
    _cylinderPipelines[layout] = CylinderPipelines{_cylinderImpostorPipeline, _cylinderOitPipeline};
  } else {
    [_cylinderOitPipeline release];
    _cylinderOitPipeline = nil;
  }
}

// Release every cached cylinder pipeline (MRC: clear() does not send -release)
// and drop the aliases pointing into the cache.
void RendererMetal::releaseCylinderPipelines() {
  for (auto& kv : _cylinderPipelines) {
    [kv.second.opaque release];
    [kv.second.oit release];
  }
  _cylinderPipelines.clear();
  _cylinderImpostorPipeline = nil;
  _cylinderOitPipeline = nil;
}

void RendererMetal::drawCylinderImpostors(const CylinderImpostorDrawCall& call) {
  if (!call.vdata || !call.idata || call.indexCount <= 0) return;
  if (call.v1Off < 0 || call.v2Off < 0 || call.radiusOff < 0 || call.flagsOff < 0 ||
      call.colorOff < 0 || call.color2Off < 0)
    return;
  if (call.flagsIsFloat) {
    // The MSL declares attr_flags as `uchar`; a Float-format flags VBO would
    // mismatch. This layout isn't produced by the current rep; bail loudly.
    NSLog(@"RendererMetal: cyl impostor unexpected Float attr_flags; skipping");
    return;
  }
  ensureEncoder();
  if (!_encoder) return;
  buildCylinderImpostorPipeline(call);
  if (!_cylinderImpostorPipeline) return;
  if (_oitActive && !_cylinderOitPipeline) return;  // no OIT variant: skip

  id<MTLBuffer> vbo = nil, ibo = nil;
  {
    auto it = _vboCache.find(call.vdata);
    if (it != _vboCache.end())
      vbo = it->second;
    else {
      vbo = [_device newBufferWithBytes:call.vdata
                                 length:call.vdataSize
                                options:MTLResourceStorageModeShared];
      if (vbo) _vboCache[call.vdata] = vbo;
    }
  }
  {
    auto it = _vboCache.find(call.idata);
    if (it != _vboCache.end())
      ibo = it->second;
    else {
      ibo = [_device newBufferWithBytes:call.idata
                                 length:call.idataSize
                                options:MTLResourceStorageModeShared];
      if (ibo) _vboCache[call.idata] = ibo;
    }
  }
  if (!vbo || !ibo) return;

  if (_oitActive) {
    [_encoder setRenderPipelineState:_cylinderOitPipeline];
    [_encoder setDepthStencilState:oitDepthState()];
  } else {
    [_encoder setRenderPipelineState:_cylinderImpostorPipeline];
    applyDepthStencilState();
    if (_depthStencilState) [_encoder setDepthStencilState:_depthStencilState];
  }
  [_encoder setVertexBuffer:vbo offset:0 atIndex:0];

  struct {
    float modelview[16];
    float projection[16];
    float uni_radius;
    float ortho;
    float depthZeroToOne;
    float no_flat_caps;
    float cap_const;
    float half_bond;
    float inv_height;
    float interiorCap;
    float interiorColor[4];
    float lAmbient, lDirect, lReflect, lSpecular, lShininess;
    float klx, kly, klz;  // key-light dir (toward light, eye space)
  } u;
  std::memcpy(u.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(u.projection, _projectionMatrix.data(), 64);
  u.lAmbient = _lightAmbient;
  u.lDirect = _lightDirect;
  u.lReflect = _lightReflect;
  u.lSpecular = _lightSpecular;
  u.lShininess = _lightShininess;
  u.klx = _keyLightEye[0];
  u.kly = _keyLightEye[1];
  u.klz = _keyLightEye[2];
  u.uni_radius = call.uniRadius;
  u.ortho = (float)call.ortho;
  u.depthZeroToOne = 0.0f;  // GL-convention clip Z (matches the sphere path)
  u.no_flat_caps = (float)call.noFlatCaps;
  // Negative tells the shader to read a_cap per-vertex (attribute 6) instead.
  // Sticks always land here: their caps vary per bond, so the CGO bakes a_cap
  // per-vertex and never emits a constant (GitHub issue #441).
  u.cap_const = (call.capOff >= 0) ? -1.0f : call.capConst;
  u.half_bond = 0.0f;   // smooth_half_bonds default off
  u.inv_height = 1.0f;  // only used when half_bond != 0
  // Cap the slab cross-section only in the opaque pass (not OIT).
  u.interiorCap = (call.interiorCap && !_oitActive) ? 1.0f : 0.0f;
  // ray_interior_color: .a>0.5 => use this rgb for the cap (else bond*0.45).
  u.interiorColor[0] = _capColor[0];
  u.interiorColor[1] = _capColor[1];
  u.interiorColor[2] = _capColor[2];
  u.interiorColor[3] = _capColorOverride ? 1.0f : 0.0f;
  [_encoder setVertexBytes:&u length:sizeof(u) atIndex:1];
  [_encoder setFragmentBytes:&u length:sizeof(u) atIndex:1];

  [_encoder drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                       indexCount:call.indexCount
                        indexType:MTLIndexTypeUInt32
                      indexBuffer:ibo
                indexBufferOffset:0];
}
// Metal tessellation pipeline: each cubic Bezier patch (4 control points) is
// tessellated over a QUAD domain — u runs along the curve, v around a ring —
// to extrude a smooth tube. (Metal has no isoline tessellation, so a tube via
// quad tessellation is both the natural fit and the showcase.)
static NSString* const kBezierTubeSrc = @R"(
using namespace metal;

struct CP { float3 pos [[attribute(0)]]; };
struct TubeU {
  float4x4 modelview;
  float4x4 projection;
  float radius;
  float _p0, _p1, _p2;
  float4 color;
};
struct TubeOut {
  float4 position [[position]];
  float3 eyeNormal;
  float4 color;
};

[[patch(quad, 4)]]
vertex TubeOut bezier_tube_vertex(
    patch_control_point<CP> cp [[stage_in]],
    float2 uv [[position_in_patch]],
    constant TubeU& U [[buffer(1)]]) {
  float t = uv.x;
  float omt = 1.0 - t;
  float3 p0 = cp[0].pos, p1 = cp[1].pos, p2 = cp[2].pos, p3 = cp[3].pos;
  // cubic Bezier position + tangent (derivative)
  float3 P = omt*omt*omt*p0 + 3.0*omt*omt*t*p1 + 3.0*omt*t*t*p2 + t*t*t*p3;
  float3 T = 3.0*omt*omt*(p1-p0) + 6.0*omt*t*(p2-p1) + 3.0*t*t*(p3-p2);
  // Robust tangent: fall back to the chord if the derivative degenerates
  // (coincident control points at sharp turns) to avoid collapsed/flipped
  // rings and gaps at joints.
  if (dot(T, T) < 1e-6) T = p3 - p0;
  T = normalize(T);
  // a stable frame around the tangent
  float3 up = (abs(T.y) < 0.99) ? float3(0.0,1.0,0.0) : float3(1.0,0.0,0.0);
  float3 N = normalize(cross(up, T));
  float3 B = normalize(cross(T, N));
  float ang = uv.y * 6.28318530718;
  float3 radial = cos(ang) * N + sin(ang) * B;
  float3 world = P + radial * U.radius;
  float4 eye = U.modelview * float4(world, 1.0);
  TubeOut o;
  o.position = U.projection * eye;
  o.position.z = 0.5 * (o.position.z + o.position.w);  // GL z [-1,1] -> Metal [0,1]
  o.eyeNormal = normalize((U.modelview * float4(radial, 0.0)).xyz);
  o.color = U.color;
  return o;
}

struct LightU { float ambient, direct, reflect, spec, shininess;
                 float klx, kly, klz; };  // klx/y/z: key-light dir (toward light, eye space)
fragment float4 bezier_tube_fragment(TubeOut in [[stage_in]],
    constant LightU& lt [[buffer(0)]]) {
  // PyMOL two-light model driven by the live lighting settings (Scene sliders).
  // Two-sided via the view vector: orient the normal toward the camera
  // (eye-space +z) so both the outer surface and the open tube ends/interior
  // are lit, never black.
  float3 nrm = normalize(in.eyeNormal);
  if (nrm.z < 0.0) nrm = -nrm;
  float ambient = lt.ambient, direct = lt.direct, reflectv = lt.reflect;
  float spec_value = lt.spec, shininess = max(lt.shininess, 1.0);
  const float3 L0 = float3(0.0,0.0,1.0);
  const float3 L1 = normalize(float3(lt.klx, lt.kly, lt.klz));  // key light (cSetting_light)
  float intensity = ambient, specular = 0.0;
  float n0 = dot(nrm, L0);
  intensity += direct * saturate(n0);
  float n1 = dot(nrm, L1);
  intensity += reflectv * saturate(n1);
  if (n1 > 0.0) {
    float3 H1 = normalize(L1 + float3(0.0,0.0,1.0));
    specular += spec_value * pow(max(dot(nrm, H1), 0.0), shininess);
  }
  float3 rgb = in.color.rgb * min(intensity, 1.0) + specular;
  return float4(rgb, in.color.a);
}
)";

void RendererMetal::buildBezierTubePipeline() {
  if (_bezierTubePipeline) return;
  NSError* err = nil;
  id<MTLLibrary> lib = [_device newLibraryWithSource:kBezierTubeSrc options:nil error:&err];
  if (!lib) {
    NSLog(@"RendererMetal: bezier tube compile failed: %@", err);
    return;
  }
  id<MTLFunction> vfn = [lib newFunctionWithName:@"bezier_tube_vertex"];
  id<MTLFunction> ffn = [lib newFunctionWithName:@"bezier_tube_fragment"];
  if (!vfn || !ffn) {
    NSLog(@"RendererMetal: bezier tube funcs missing");
    return;
  }

  MTLVertexDescriptor* vd = [MTLVertexDescriptor vertexDescriptor];
  vd.attributes[0].format = MTLVertexFormatFloat3;  // control point position
  vd.attributes[0].offset = 0;
  vd.attributes[0].bufferIndex = 0;
  vd.layouts[0].stride = 12;
  vd.layouts[0].stepFunction = MTLVertexStepFunctionPerPatchControlPoint;

  MTLRenderPipelineDescriptor* psd = [[MTLRenderPipelineDescriptor alloc] init];
  psd.vertexFunction = vfn;
  psd.fragmentFunction = ffn;
  psd.vertexDescriptor = vd;
  psd.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
  psd.rasterSampleCount = _sampleCount;
  psd.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  psd.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float_Stencil8;
  psd.maxTessellationFactor = 64;
  psd.tessellationFactorStepFunction = MTLTessellationFactorStepFunctionConstant;
  psd.tessellationFactorFormat = MTLTessellationFactorFormatHalf;
  psd.tessellationOutputWindingOrder = MTLWindingClockwise;
  psd.tessellationPartitionMode = MTLTessellationPartitionModeInteger;
  _bezierTubePipeline = [_device newRenderPipelineStateWithDescriptor:psd error:&err];
  if (!_bezierTubePipeline) NSLog(@"RendererMetal: bezier tube pipeline failed: %@", err);
}

// half-float bit pattern for a tessellation factor
static inline uint16_t f16(float v) {
  _Float16 h = (_Float16)v;
  uint16_t bits;
  std::memcpy(&bits, &h, sizeof(bits));
  return bits;
}

void RendererMetal::drawBezierTubes(const void* cp, size_t dataSize, float radius, float r, float g,
                                    float b) {
  if (!cp || dataSize < 48) return;  // need at least one 4-point patch
  ensureEncoder();
  if (!_encoder) return;
  buildBezierTubePipeline();
  if (!_bezierTubePipeline) return;

  NSUInteger numPatches = dataSize / 48;  // 4 control points * 3 floats * 4 B
  if (numPatches == 0) return;

  // Control-point buffer (cache by pointer like the other VBOs).
  id<MTLBuffer> cpBuf = nil;
  auto it = _vboCache.find(cp);
  if (it != _vboCache.end())
    cpBuf = it->second;
  else {
    cpBuf = [_device newBufferWithBytes:cp length:dataSize options:MTLResourceStorageModeShared];
    if (!cpBuf) return;
    _vboCache[cp] = cpBuf;
  }

  // Constant per-patch tessellation factors (u = along curve, v = around ring).
  const float segF = 24.0f;   // subdivisions along the curve
  const float ringF = 14.0f;  // subdivisions around the tube
  struct QuadFactors {
    uint16_t edge[4];
    uint16_t inside[2];
  };
  if (!_bezierTessFactors || _bezierTessPatchCap < numPatches) {
    // MRC: drop the smaller buffer before growing (an in-flight command
    // buffer that still uses it keeps it alive through Metal's own retain).
    [_bezierTessFactors release];
    _bezierTessFactors = [_device newBufferWithLength:numPatches * sizeof(QuadFactors)
                                              options:MTLResourceStorageModeShared];
    if (!_bezierTessFactors) return;
    QuadFactors qf;
    qf.edge[0] = f16(ringF);
    qf.edge[1] = f16(segF);
    qf.edge[2] = f16(ringF);
    qf.edge[3] = f16(segF);
    qf.inside[0] = f16(segF);
    qf.inside[1] = f16(ringF);
    QuadFactors* dst = (QuadFactors*)_bezierTessFactors.contents;
    for (NSUInteger i = 0; i < numPatches; ++i) dst[i] = qf;
    _bezierTessPatchCap = numPatches;
  }

  [_encoder setRenderPipelineState:_bezierTubePipeline];
  // Opaque geometry: standard depth test + write.
  [_encoder setDepthStencilState:bezierDepthState()];
  [_encoder setCullMode:MTLCullModeNone];
  [_encoder setVertexBuffer:cpBuf offset:0 atIndex:0];

  struct {
    float modelview[16];
    float projection[16];
    float radius;
    float _p0, _p1, _p2;
    float color[4];
  } u;
  std::memcpy(u.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(u.projection, _projectionMatrix.data(), 64);
  u.radius = radius;
  u._p0 = u._p1 = u._p2 = 0.0f;
  u.color[0] = r;
  u.color[1] = g;
  u.color[2] = b;
  u.color[3] = 1.0f;
  [_encoder setVertexBytes:&u length:sizeof(u) atIndex:1];
  // Lighting (Scene sliders) for bezier_tube_fragment at fragment buffer(0).
  {
    struct {
      float a, d, r, s, sh, klx, kly, klz;
    } _lt = {_lightAmbient,   _lightDirect,    _lightReflect,   _lightSpecular,
             _lightShininess, _keyLightEye[0], _keyLightEye[1], _keyLightEye[2]};
    [_encoder setFragmentBytes:&_lt length:sizeof(_lt) atIndex:0];
  }

  [_encoder setTessellationFactorBuffer:_bezierTessFactors offset:0 instanceStride:0];
  [_encoder drawPatches:4
                  patchStart:0
                  patchCount:numPatches
            patchIndexBuffer:nil
      patchIndexBufferOffset:0
               instanceCount:1
                baseInstance:0];
}

}  // namespace pymol
