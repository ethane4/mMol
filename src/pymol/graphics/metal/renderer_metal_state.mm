
// The GL-style state, buffer, texture, matrix and immediate-mode API.

#include "pymol/graphics/metal/renderer_metal.h"

namespace pymol {

static MTLBlendFactor toMTLBlend(BlendFunc f) {
  switch (f) {
    case BlendFunc::Zero:
      return MTLBlendFactorZero;
    case BlendFunc::One:
      return MTLBlendFactorOne;
    case BlendFunc::SrcAlpha:
      return MTLBlendFactorSourceAlpha;
    case BlendFunc::OneMinusSrcAlpha:
      return MTLBlendFactorOneMinusSourceAlpha;
    case BlendFunc::DstAlpha:
      return MTLBlendFactorDestinationAlpha;
    case BlendFunc::OneMinusDstAlpha:
      return MTLBlendFactorOneMinusDestinationAlpha;
    case BlendFunc::SrcColor:
      return MTLBlendFactorSourceColor;
    case BlendFunc::OneMinusSrcColor:
      return MTLBlendFactorOneMinusSourceColor;
  }
  return MTLBlendFactorOne;
}

static MTLCompareFunction toMTLCompare(DepthFunc f) {
  switch (f) {
    case DepthFunc::Never:
      return MTLCompareFunctionNever;
    case DepthFunc::Less:
      return MTLCompareFunctionLess;
    case DepthFunc::Equal:
      return MTLCompareFunctionEqual;
    case DepthFunc::LessEqual:
      return MTLCompareFunctionLessEqual;
    case DepthFunc::Greater:
      return MTLCompareFunctionGreater;
    case DepthFunc::NotEqual:
      return MTLCompareFunctionNotEqual;
    case DepthFunc::GreaterEqual:
      return MTLCompareFunctionGreaterEqual;
    case DepthFunc::Always:
      return MTLCompareFunctionAlways;
  }
  return MTLCompareFunctionLess;
}
void RendererMetal::viewport(int x, int y, int w, int h) {
  _viewport = {(double)x, (double)y, (double)w, (double)h, 0.0, 1.0};

  // Also update scissor to match if scissor is not explicitly set
  _scissorRect = {static_cast<NSUInteger>(x), static_cast<NSUInteger>(y),
                  static_cast<NSUInteger>(w), static_cast<NSUInteger>(h)};

  if (_encoder) {
    [_encoder setViewport:_viewport];
    if (!_scissorEnabled) {
      [_encoder setScissorRect:_scissorRect];
    }
  }
}

bool RendererMetal::getViewportRect(int& x, int& y, int& w, int& h) const {
  x = static_cast<int>(_viewport.originX);
  y = static_cast<int>(_viewport.originY);
  w = static_cast<int>(_viewport.width);
  h = static_cast<int>(_viewport.height);
  return true;
}

void RendererMetal::clear(bool color, bool depth, bool stencil) {
  // In Metal, clears happen via the render pass descriptor's load actions.
  // If we already have an encoder, we need to end it and start a new pass.
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }

  if (_passDesc) {
    if (color) {
      _passDesc.colorAttachments[0].loadAction = MTLLoadActionClear;
      _passDesc.colorAttachments[0].clearColor =
          MTLClearColorMake(_clearR, _clearG, _clearB, _clearA);
    } else {
      _passDesc.colorAttachments[0].loadAction = MTLLoadActionLoad;
    }

    if (depth) {
      _passDesc.depthAttachment.loadAction = MTLLoadActionClear;
      _passDesc.depthAttachment.clearDepth = 1.0;
    } else {
      _passDesc.depthAttachment.loadAction = MTLLoadActionLoad;
    }

    if (stencil) {
      _passDesc.stencilAttachment.loadAction = MTLLoadActionClear;
      _passDesc.stencilAttachment.clearStencil = 0;
    } else {
      _passDesc.stencilAttachment.loadAction = MTLLoadActionLoad;
    }
  }

  // The next ensureEncoder() call will create a new encoder with these
  // load actions, effectively performing the clear.
}

void RendererMetal::clearColor(float r, float g, float b, float a) {
  _clearR = r;
  _clearG = g;
  _clearB = b;
  _clearA = a;
}

void RendererMetal::scissor(int x, int y, int w, int h) {
  _scissorRect = {static_cast<NSUInteger>(x), static_cast<NSUInteger>(y),
                  static_cast<NSUInteger>(w), static_cast<NSUInteger>(h)};

  if (_encoder && _scissorEnabled) {
    [_encoder setScissorRect:_scissorRect];
  }
}
void RendererMetal::enable(Capability cap) {
  switch (cap) {
    case Capability::DepthTest:
      _depthTestEnabled = true;
      _depthStencilDirty = true;
      applyDepthStencilState();
      break;
    case Capability::Blend:
      _blendEnabled = true;
      // Blend is baked into pipeline state — will apply on next draw
      break;
    case Capability::CullFace:
      _cullFaceEnabled = true;
      if (_encoder) [_encoder setCullMode:MTLCullModeBack];
      break;
    case Capability::ScissorTest:
      _scissorEnabled = true;
      if (_encoder) [_encoder setScissorRect:_scissorRect];
      break;
    case Capability::StencilTest:
      _stencilTestEnabled = true;
      break;
    case Capability::Lighting:
      _lightingEnabled = true;
      break;
    case Capability::Fog:
      _fogEnabled = true;
      break;
    // These have no direct Metal equivalent — tracked as flags for shaders
    case Capability::LineSmooth:
    case Capability::Normalize:
    case Capability::ColorMaterial:
    case Capability::AlphaTest:
    case Capability::PolygonOffset:
    case Capability::Texture2D:
      break;
  }
}

void RendererMetal::disable(Capability cap) {
  switch (cap) {
    case Capability::DepthTest:
      _depthTestEnabled = false;
      _depthStencilDirty = true;
      applyDepthStencilState();
      break;
    case Capability::Blend:
      _blendEnabled = false;
      break;
    case Capability::CullFace:
      _cullFaceEnabled = false;
      if (_encoder) [_encoder setCullMode:MTLCullModeNone];
      break;
    case Capability::ScissorTest:
      _scissorEnabled = false;
      // Reset scissor to viewport
      if (_encoder) {
        MTLScissorRect fullRect = {
            static_cast<NSUInteger>(_viewport.originX), static_cast<NSUInteger>(_viewport.originY),
            static_cast<NSUInteger>(_viewport.width), static_cast<NSUInteger>(_viewport.height)};
        [_encoder setScissorRect:fullRect];
      }
      break;
    case Capability::StencilTest:
      _stencilTestEnabled = false;
      break;
    case Capability::Lighting:
      _lightingEnabled = false;
      break;
    case Capability::Fog:
      _fogEnabled = false;
      break;
    case Capability::LineSmooth:
    case Capability::Normalize:
    case Capability::ColorMaterial:
    case Capability::AlphaTest:
    case Capability::PolygonOffset:
    case Capability::Texture2D:
      break;
  }
}

void RendererMetal::blendFunc(BlendFunc src, BlendFunc dst) {
  _blendSrcFactor = toMTLBlend(src);
  _blendDstFactor = toMTLBlend(dst);
}

void RendererMetal::depthFunc(DepthFunc func) {
  _depthCompareFunc = toMTLCompare(func);
  _depthStencilDirty = true;
  applyDepthStencilState();
}

void RendererMetal::depthMask(bool write) {
  _depthWriteEnabled = write;
  _depthStencilDirty = true;
  applyDepthStencilState();
}

void RendererMetal::colorMask(bool r, bool g, bool b, bool a) {
  _colorWriteMask = MTLColorWriteMaskNone;
  if (r) _colorWriteMask |= MTLColorWriteMaskRed;
  if (g) _colorWriteMask |= MTLColorWriteMaskGreen;
  if (b) _colorWriteMask |= MTLColorWriteMaskBlue;
  if (a) _colorWriteMask |= MTLColorWriteMaskAlpha;
  // Color mask is baked into pipeline state
}

void RendererMetal::lineWidth(float w) {
  _lineWidth = w;
  // Metal does not support variable line width — tracked for compatibility
}

void RendererMetal::pointSize(float s) {
  _pointSize = s;
  // Point size is set via vertex shader in Metal
}

void RendererMetal::setDepthClamp(bool enabled) {
  // MTLDepthClipModeClamp rasterizes primitives outside the near/far planes with
  // clamped depth instead of clipping them. Used for selection-indicator points
  // so a selected atom's marker stays on screen even when the atom falls outside
  // the clip slab. Encoder state — scope it around the draw and reset after.
  if (_encoder) [_encoder setDepthClipMode:enabled ? MTLDepthClipModeClamp : MTLDepthClipModeClip];
}
void RendererMetal::drawArrays(PrimitiveType mode, int first, int count) {
  ensureEncoder();
  if (!_encoder) return;

  // Bind the current array buffer at index 0 if one is bound
  auto it = _buffers.find(_boundArrayBuffer);
  if (it != _buffers.end()) {
    [_encoder setVertexBuffer:it->second offset:0 atIndex:0];
  }

  // Set current pipeline if available
  if (_currentPipeline) {
    [_encoder setRenderPipelineState:_currentPipeline];
  }

  [_encoder drawPrimitives:toMTL(mode)
               vertexStart:static_cast<NSUInteger>(first)
               vertexCount:static_cast<NSUInteger>(count)];
}

void RendererMetal::drawElements(PrimitiveType mode, int count, const void* indices) {
  ensureEncoder();
  if (!_encoder) return;

  // Bind vertex buffer
  auto vit = _buffers.find(_boundArrayBuffer);
  if (vit != _buffers.end()) {
    [_encoder setVertexBuffer:vit->second offset:0 atIndex:0];
  }

  if (_currentPipeline) {
    [_encoder setRenderPipelineState:_currentPipeline];
  }

  // If an element array buffer is bound, use it
  auto eit = _buffers.find(_boundElementBuffer);
  if (eit != _buffers.end()) {
    [_encoder drawIndexedPrimitives:toMTL(mode)
                         indexCount:static_cast<NSUInteger>(count)
                          indexType:MTLIndexTypeUInt32
                        indexBuffer:eit->second
                  indexBufferOffset:0];
  } else if (indices) {
    // Create a temporary index buffer from the provided pointer
    NSUInteger byteSize = static_cast<NSUInteger>(count) * sizeof(uint32_t);
    id<MTLBuffer> tmpIdx = [_device newBufferWithBytes:indices
                                                length:byteSize
                                               options:MTLResourceStorageModeShared];
    [_encoder drawIndexedPrimitives:toMTL(mode)
                         indexCount:static_cast<NSUInteger>(count)
                          indexType:MTLIndexTypeUInt32
                        indexBuffer:tmpIdx
                  indexBufferOffset:0];
    // MRC: the encoder retains the buffer until the command buffer completes, so
    // this +1 transient can be released now (no caching for this client path).
    [tmpIdx release];
  }
}
uint32_t RendererMetal::createBuffer() {
  uint32_t id = _nextBufferId++;
  // No Metal buffer allocated yet — created on bufferData()
  return id;
}

void RendererMetal::deleteBuffer(uint32_t id) {
  _buffers.erase(id);
  if (_boundArrayBuffer == id) _boundArrayBuffer = 0;
  if (_boundElementBuffer == id) _boundElementBuffer = 0;
}

void RendererMetal::bindBuffer(BufferTarget target, uint32_t id) {
  switch (target) {
    case BufferTarget::Array:
      _boundArrayBuffer = id;
      break;
    case BufferTarget::ElementArray:
      _boundElementBuffer = id;
      break;
  }
}

void RendererMetal::bufferData(BufferTarget target, size_t size, const void* data,
                               BufferUsage /*usage*/) {
  uint32_t boundId = (target == BufferTarget::Array) ? _boundArrayBuffer : _boundElementBuffer;
  if (boundId == 0 || size == 0) return;

  id<MTLBuffer> mtlBuf;
  if (data) {
    mtlBuf = [_device newBufferWithBytes:data length:size options:MTLResourceStorageModeShared];
  } else {
    mtlBuf = [_device newBufferWithLength:size options:MTLResourceStorageModeShared];
  }

  _buffers[boundId] = mtlBuf;
}
void RendererMetal::vertexAttribPointer(int index, int size, int type, bool normalized, int stride,
                                        const void* offset) {
  if (index < 0 || index >= kMaxVertexAttribs) return;

  auto& attr = _vertexAttribs[index];
  attr.size = size;
  attr.type = type;
  attr.normalized = normalized;
  attr.stride = stride;
  attr.offset = reinterpret_cast<uintptr_t>(offset);
}

void RendererMetal::enableVertexAttribArray(int index) {
  if (index >= 0 && index < kMaxVertexAttribs) {
    _vertexAttribs[index].enabled = true;
  }
}

void RendererMetal::disableVertexAttribArray(int index) {
  if (index >= 0 && index < kMaxVertexAttribs) {
    _vertexAttribs[index].enabled = false;
  }
}
void RendererMetal::useProgram(uint32_t programId) {
  _currentProgram = programId;
  // Pipeline state lookup will be done by ShaderMgr or a future pipeline
  // cache that maps programId → MTLRenderPipelineState.
}

void RendererMetal::setUniform1i(int location, int v) {
  if (location < 0 || location >= kMaxUniforms) return;
  _uniformData[location * 4] = static_cast<float>(v);
}

void RendererMetal::setUniform1f(int location, float v) {
  if (location < 0 || location >= kMaxUniforms) return;
  _uniformData[location * 4] = v;
}

void RendererMetal::setUniform2f(int location, float v0, float v1) {
  if (location < 0 || location >= kMaxUniforms) return;
  _uniformData[location * 4 + 0] = v0;
  _uniformData[location * 4 + 1] = v1;
}

void RendererMetal::setUniform3f(int location, float v0, float v1, float v2) {
  if (location < 0 || location >= kMaxUniforms) return;
  _uniformData[location * 4 + 0] = v0;
  _uniformData[location * 4 + 1] = v1;
  _uniformData[location * 4 + 2] = v2;
}

void RendererMetal::setUniform4f(int location, float v0, float v1, float v2, float v3) {
  if (location < 0 || location >= kMaxUniforms) return;
  _uniformData[location * 4 + 0] = v0;
  _uniformData[location * 4 + 1] = v1;
  _uniformData[location * 4 + 2] = v2;
  _uniformData[location * 4 + 3] = v3;
}

void RendererMetal::setUniformMatrix4fv(int location, const float* value) {
  if (location < 0 || !value) return;
  // Store 16 floats starting at this location's slot
  // Each location occupies 4 float slots, a 4x4 matrix needs 4 locations
  if (location + 3 >= kMaxUniforms) return;
  std::memcpy(&_uniformData[location * 4], value, 16 * sizeof(float));
}

void RendererMetal::setUniformMatrix3fv(int location, const float* value) {
  if (location < 0 || !value) return;
  // 3x3 matrix = 9 floats, stored in 3 vec4 slots (padded)
  if (location + 2 >= kMaxUniforms) return;
  for (int col = 0; col < 3; ++col) {
    _uniformData[(location + col) * 4 + 0] = value[col * 3 + 0];
    _uniformData[(location + col) * 4 + 1] = value[col * 3 + 1];
    _uniformData[(location + col) * 4 + 2] = value[col * 3 + 2];
    _uniformData[(location + col) * 4 + 3] = 0.0f;
  }
}
uint32_t RendererMetal::createTexture() {
  uint32_t id = _nextTextureId++;
  // Texture is allocated on first texImage call (not part of Renderer yet)
  return id;
}

void RendererMetal::deleteTexture(uint32_t id) {
  _textures.erase(id);
  if (_boundTexture == id) _boundTexture = 0;
}

void RendererMetal::bindTexture(TextureTarget /*target*/, uint32_t id) { _boundTexture = id; }

void RendererMetal::activeTexture(int unit) { _activeTextureUnit = unit; }

void RendererMetal::texParameteri(TextureTarget /*target*/, int /*pname*/, int /*param*/) {
  // Texture parameters are set via MTLSamplerDescriptor in Metal.
  // Will be handled when we build the sampler state cache.
}
uint32_t RendererMetal::createFramebuffer() {
  uint32_t id = _nextFBOId++;
  // Actual textures created when attachments are configured
  return id;
}

void RendererMetal::deleteFramebuffer(uint32_t id) {
  _fbColorAttachments.erase(id);
  _fbDepthAttachments.erase(id);
  if (_boundFBO == id) _boundFBO = 0;
}

void RendererMetal::bindFramebuffer(uint32_t id) {
  _boundFBO = id;

  // When binding a non-default FBO, we need to end the current encoder
  // and start a new render pass targeting the FBO's textures.
  if (id != 0 && _encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }

  // For FBO 0, the pass descriptor points to the drawable (default).
  // For custom FBOs, we'd configure _passDesc to point to their textures.
  // Full implementation deferred to framebuffer attachment methods.
}
void RendererMetal::matrixMode(int mode) {
  // GL_MODELVIEW = 0x1700, GL_PROJECTION = 0x1701
  _matrixMode = (mode == 0x1701) ? 1 : 0;
}

void RendererMetal::loadIdentity() {
  if (_matrixMode == 0) {
    _modelviewMatrix = identityMatrix();
  } else {
    _projectionMatrix = identityMatrix();
  }
}

void RendererMetal::loadMatrixf(const float* m) {
  if (!m) return;
  Mat4& mat = (_matrixMode == 0) ? _modelviewMatrix : _projectionMatrix;
  std::memcpy(mat.data(), m, 16 * sizeof(float));
}

void RendererMetal::pushMatrix() {
  if (_matrixMode == 0) {
    _modelviewStack.push(_modelviewMatrix);
  } else {
    _projectionStack.push(_projectionMatrix);
  }
}

void RendererMetal::popMatrix() {
  if (_matrixMode == 0) {
    if (!_modelviewStack.empty()) {
      _modelviewMatrix = _modelviewStack.top();
      _modelviewStack.pop();
    }
  } else {
    if (!_projectionStack.empty()) {
      _projectionMatrix = _projectionStack.top();
      _projectionStack.pop();
    }
  }
}

void RendererMetal::translatef(float x, float y, float z) {
  Mat4& mat = (_matrixMode == 0) ? _modelviewMatrix : _projectionMatrix;
  mat = multiplyMatrices(mat, translationMatrix(x, y, z));
}

void RendererMetal::scalef(float x, float y, float z) {
  Mat4& mat = (_matrixMode == 0) ? _modelviewMatrix : _projectionMatrix;
  mat = multiplyMatrices(mat, scaleMatrix(x, y, z));
}

void RendererMetal::multMatrixf(const float* m) {
  if (!m) return;
  Mat4& mat = (_matrixMode == 0) ? _modelviewMatrix : _projectionMatrix;
  Mat4 other;
  std::memcpy(other.data(), m, 16 * sizeof(float));
  mat = multiplyMatrices(mat, other);
}
void RendererMetal::beginBatch(PrimitiveType mode) {
  _batchMode = mode;
  _batchVertices.clear();
}

void RendererMetal::batchVertex3f(float x, float y, float z) {
  _batchVertices.push_back({x, y, z, _curR, _curG, _curB, _curA, _curNX, _curNY, _curNZ});
}

void RendererMetal::batchVertex3fv(const float* v) { batchVertex3f(v[0], v[1], v[2]); }

void RendererMetal::batchVertex2f(float x, float y) { batchVertex3f(x, y, 0.0f); }

void RendererMetal::batchVertex2i(int x, int y) {
  batchVertex3f(static_cast<float>(x), static_cast<float>(y), 0.0f);
}

void RendererMetal::batchColor3f(float r, float g, float b) {
  _curR = r;
  _curG = g;
  _curB = b;
  _curA = 1.0f;
}

void RendererMetal::batchColor3fv(const float* c) { batchColor3f(c[0], c[1], c[2]); }

void RendererMetal::batchColor4f(float r, float g, float b, float a) {
  _curR = r;
  _curG = g;
  _curB = b;
  _curA = a;
}

void RendererMetal::batchColor4fv(const float* c) { batchColor4f(c[0], c[1], c[2], c[3]); }

void RendererMetal::batchColor4ub(unsigned char r, unsigned char g, unsigned char b,
                                  unsigned char a) {
  batchColor4f(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
}

void RendererMetal::batchNormal3fv(const float* n) {
  _curNX = n[0];
  _curNY = n[1];
  _curNZ = n[2];
}

void RendererMetal::endBatch() {
  if (_batchVertices.empty()) return;

  ensureEncoder();
  if (!_encoder || !_batchPipeline) return;

  // Handle triangle fan conversion: fan vertices [0,1,2,3,...,N] become
  // triangles [0,1,2], [0,2,3], [0,3,4], ..., [0,N-1,N]
  std::vector<BatchVertex>* drawVertices = &_batchVertices;
  std::vector<BatchVertex> converted;

  if (_batchMode == PrimitiveType::TriangleFan && _batchVertices.size() >= 3) {
    converted.reserve((_batchVertices.size() - 2) * 3);
    for (size_t i = 1; i + 1 < _batchVertices.size(); ++i) {
      converted.push_back(_batchVertices[0]);
      converted.push_back(_batchVertices[i]);
      converted.push_back(_batchVertices[i + 1]);
    }
    drawVertices = &converted;
  } else if (_batchMode == PrimitiveType::Quads && _batchVertices.size() >= 4) {
    // Convert quads to triangles: [0,1,2,3] → [0,1,2], [0,2,3]
    converted.reserve((_batchVertices.size() / 4) * 6);
    for (size_t i = 0; i + 3 < _batchVertices.size(); i += 4) {
      converted.push_back(_batchVertices[i + 0]);
      converted.push_back(_batchVertices[i + 1]);
      converted.push_back(_batchVertices[i + 2]);
      converted.push_back(_batchVertices[i + 0]);
      converted.push_back(_batchVertices[i + 2]);
      converted.push_back(_batchVertices[i + 3]);
    }
    drawVertices = &converted;
  }

  NSUInteger byteSize = drawVertices->size() * sizeof(BatchVertex);
  // A single frame can issue MULTIPLE batches (e.g. the hover-preview points
  // then the committed-selection points). All draws are queued into one command
  // buffer and execute together at commit, so a SHARED buffer reused at offset 0
  // would leave every draw reading the LAST batch's vertices — the earlier
  // batches' points silently vanish (this is why the hover preview disappeared
  // whenever a committed selection was enabled). Allocate a fresh transient
  // buffer per batch; the encoder retains it until the frame completes, so each
  // draw reads its own vertices. (CGO uses its own batch object, and the
  // renderer's immediate path issues only a handful of batches per frame, so the
  // per-batch allocation is negligible.)
  id<MTLBuffer> tmpVBO = [_device newBufferWithLength:std::max(byteSize, (NSUInteger)16)
                                              options:MTLResourceStorageModeShared];
  memcpy(tmpVBO.contents, drawVertices->data(), byteSize);

  // Bind batch VBO at buffer index 0
  [_encoder setVertexBuffer:tmpVBO offset:0 atIndex:0];
  [tmpVBO release];  // MRC: the encoder holds its own reference for the frame

  // Upload modelview, projection matrices, and point size as uniforms
  struct {
    float modelview[16];
    float projection[16];
    float pointSize;
    float _pad[3];
  } uniforms;
  std::memcpy(uniforms.modelview, _modelviewMatrix.data(), 64);
  std::memcpy(uniforms.projection, _projectionMatrix.data(), 64);
  uniforms.pointSize = _pointSize > 0.0f ? _pointSize : 1.0f;
  uniforms._pad[0] = uniforms._pad[1] = uniforms._pad[2] = 0.0f;
  if (selectionMarkers && _batchMode == PrimitiveType::Points) uniforms._pad[0] = 1.0f;
  [_encoder setVertexBytes:&uniforms length:sizeof(uniforms) atIndex:1];
  [_encoder setFragmentBytes:&uniforms length:sizeof(uniforms) atIndex:1];
  // Lighting (Scene sliders) for the lit vbo_fragment / vbo_fragment_oit at
  // fragment buffer(0). Harmless for the unlit/flat pipelines (they don't read
  // it). Scoped so the local doesn't clash across multiple draw paths.
  {
    struct {
      float a, d, r, s, sh, klx, kly, klz;
    } _lt = {_lightAmbient,   _lightDirect,    _lightReflect,   _lightSpecular,
             _lightShininess, _keyLightEye[0], _keyLightEye[1], _keyLightEye[2]};
    [_encoder setFragmentBytes:&_lt length:sizeof(_lt) atIndex:0];
  }

  // Always use the built-in batch pipeline for batch rendering.
  // _currentPipeline is for VBO draws with a different vertex layout.
  [_encoder setRenderPipelineState:_batchPipeline];

  MTLPrimitiveType mtlPrim =
      (_batchMode == PrimitiveType::TriangleFan || _batchMode == PrimitiveType::Quads)
          ? MTLPrimitiveTypeTriangle
          : toMTL(_batchMode);

  [_encoder drawPrimitives:mtlPrim
               vertexStart:0
               vertexCount:static_cast<NSUInteger>(drawVertices->size())];

  _batchVertices.clear();
}
bool RendererMetal::isRenderReady() const { return _cmdBuffer && _passDesc && _batchPipeline; }

bool RendererMetal::hasActiveEncoder() const {
  // The encoder must exist AND the command buffer must not yet be committed.
  // After endFrame() commits the buffer, _cmdBuffer is nilled. But if
  // SceneRenderAll ends/restarts the encoder, _encoder may be nil while
  // _cmdBuffer is still valid — that's fine, ensureEncoder() can recreate it.
  // We only reject the case where _cmdBuffer itself is gone or committed.
  if (!_cmdBuffer) return false;
  MTLCommandBufferStatus status = [_cmdBuffer status];
  return status < MTLCommandBufferStatusCommitted;
}
void RendererMetal::getIntegerv(int /*pname*/, int* params) {
  // Stub — most GL queries don't map to Metal.
  // Specific queries will be added as needed.
  if (params) *params = 0;
}

const char* RendererMetal::getString(int /*name*/) { return "Metal"; }

int RendererMetal::getError() {
  return 0;  // Metal uses NSError, not error codes
}
void RendererMetal::flush() {
  // Metal command buffers are committed in endFrame().
  // No equivalent of glFlush needed.
}

void RendererMetal::finish() {
  if (_cmdBuffer) {
    [_cmdBuffer waitUntilCompleted];
  }
}

void RendererMetal::readPixels(int x, int y, int w, int h, int /*format*/, int /*type*/,
                               void* pixels) {
  if (!pixels || !_drawable) return;

  // Reading pixels from the drawable texture requires a blit encoder.
  // This is used for ray tracing output and screenshot capture.
  id<MTLTexture> tex = _drawable.texture;
  if (!tex) return;

  // End current render encoder if active
  if (_encoder) {
    [_encoder endEncoding];
    _encoder = nil;
  }

  MTLRegion region = MTLRegionMake2D(static_cast<NSUInteger>(x), static_cast<NSUInteger>(y),
                                     static_cast<NSUInteger>(w), static_cast<NSUInteger>(h));
  NSUInteger bytesPerRow = static_cast<NSUInteger>(w) * 4;  // RGBA8

  [tex getBytes:pixels bytesPerRow:bytesPerRow fromRegion:region mipmapLevel:0];
}

void RendererMetal::pixelStorei(int /*pname*/, int /*param*/) {
  // Metal doesn't have pixel store state — alignment is handled
  // explicitly when creating textures and reading pixels.
}

}  // namespace pymol
