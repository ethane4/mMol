#include "pymol/graphics/post_process.h"

#include "pymol/platform/gl_def.h"
#include "pymol/graphics/generic_buffer.h"
#include "pymol/platform/os_gl.h"
#include "pymol/memory.h"

void PostProcess::activateRTAsTexture(std::size_t idx, GLuint textureUnit) {
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  auto& rt = m_renderTargets[idx];
  if (rt && rt->textures().front()) {
    rt->textures().front()->bind();
  }
}

void PostProcess::bindRT(std::size_t idx, bool clear) {
  if (idx < m_renderTargets.size()) {
    auto& t = m_renderTargets[idx];
    if (t) {
      t->bind(clear);
    }
  }
}

void PostProcess::bindFBORBO(std::size_t idx) {
  if (idx < m_renderTargets.size()) {
    auto& t = m_renderTargets[idx];
    if (t) {
      t->bindFBORBO();
    }
  }
}

const RenderTargetGL::shape_type PostProcess::size(
    std::size_t idx) const noexcept {
  return m_renderTargets[idx]->size();
}

OIT_PostProcess::OIT_PostProcess(int width, int height, RenderbufferGL* rbo) {
  if (TM3_IS_ONEBUF) {
    auto rt0 = std::make_unique<RenderTargetGL>(width, height);
    rt0->layout({{4, rt_layout_t::FLOAT}}, rbo);
    m_renderTargets.push_back(std::move(rt0));

    auto rt1 = std::make_unique<RenderTargetGL>(width, height);
    rt1->layout({{1, rt_layout_t::FLOAT}}, rt0->rbo());
    m_renderTargets.push_back(std::move(rt1));
  } else {
    std::vector<rt_layout_t> layouts;
    layouts.emplace_back(4, rt_layout_t::FLOAT);
    if (GLEW_VERSION_3_0) {
      layouts.emplace_back(1, rt_layout_t::FLOAT);
    } else {
      layouts.emplace_back(2, rt_layout_t::FLOAT);
    }

    auto rt0 = std::make_unique<RenderTargetGL>(width, height);
    rt0->layout(std::move(layouts), rbo);
    m_renderTargets.push_back(std::move(rt0));
  }
}

void OIT_PostProcess::bindRT(std::size_t idx, bool clear) {}

void OIT_PostProcess::activateRTAsTexture(std::size_t idx, GLuint textureUnit) {
  glActiveTexture(GL_TEXTURE0 + textureUnit);
  if (TM3_IS_ONEBUF) {
    auto& rt = m_renderTargets[idx];
    auto& tex = rt->textures().front();
    if (tex) {
      tex->bind();
    }
  } else {
    auto& rt = m_renderTargets.front();
    if (rt) {
      rt->textures()[idx]->bind();
    }
  }
}
