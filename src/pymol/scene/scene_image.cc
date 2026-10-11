/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* sc
-*
-*
Z* -------------------------------------------------------------------
*/

// Capturing, sizing and saving scene images, and drawing image overlays.

#include "pymol/utility.h"
#include "pymol/scene/scene.h"
#include "pymol/ray/scene_ray.h"
#include "pymol/scene/movie.h"
#include "pymol/image/my_png.h"
#include "pymol/python/p.h"
#include "pymol/executive/executive.h"
#include "pymol/text/text.h"
#include "pymol/base/feedback.h"
#include "pymol/base/util2.h"
#include "pymol/scene/scene_internal.h"

int SceneCaptureWindow(PyMOLGlobals* G) {
  CScene* I = G->Scene;
  int ok = true;

  /* check assumptions */

  if (ok && G->HaveGUI && G->ValidContext) {
    auto drawBuffer = SceneDrawBothGetConfig(G);

    ScenePurgeImage(G);

    SceneCopy(G, drawBuffer, true, true);

    if (!I->Image) ok = false;

    if (ok && I->Image) {
      I->DirtyFlag = false;
      I->CopyType = 2; /* suppresses display of copied image */
      if (SettingGetGlobal_b(G, cSetting_opaque_background))
        I->Image->m_needs_alpha_reset = true;
    }
  } else {
    ok = false;
  }
  return ok;
}

/**
 * Sets the Scene's cached copy image
 * @param image the image to set
 * @param dirty whether the image is dirty (needed?)
 * @param copyForced whether the image was forced to be copied
 */
static void SceneSetCopyImage(PyMOLGlobals* G, pymol::Image image, bool dirty,
                              bool copyForced) {
  auto I = G->Scene;
  I->Image = std::make_shared<pymol::Image>(std::move(image));
  I->CopyType = true;
  I->CopyForced = copyForced;
  I->DirtyFlag = dirty;
}

/**
 * Get the maximum dimensions of the viewport
 * @return the maximum dimensions of OpenGL viewport
 */
Extent2D SceneGLGetMaxDimensions(PyMOLGlobals* G) {
  // TODO: Move this function to OpenGL Context manager
  GLint dims[2];
  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, reinterpret_cast<GLint*>(&dims));
  return Extent2D{static_cast<std::uint32_t>(dims[0]),
                  static_cast<std::uint32_t>(dims[1])};
}

Extent2D ExtentClampByAspectRatio(Extent2D extent, const Extent2D& maxDim) {
  float extentAspect = static_cast<float>(extent.width) / extent.height;
  if (extent.width > maxDim.width) {
    extent.height = static_cast<std::uint32_t>(maxDim.width / extentAspect);
    extent.width = maxDim.width;
  }
  if (extent.height > maxDim.height) {
    extent.width = static_cast<std::uint32_t>(maxDim.height * extentAspect);
    extent.height = maxDim.height;
  }
  return extent;
}

/**
 * @brief Returns supersampled image
 * @param src source image
 * @param factor supersampling factor
 * @param shift shift value
 * @return supersampled image
 */
static pymol::Image SceneSupersampleImage(
    const pymol::Image& src, const UpscaledExtentInfo& upscaledExtent) {
  auto factor = upscaledExtent.factor;
  auto shift = upscaledExtent.shift;
  auto oriWidth = src.getWidth();
  auto oriHeight = src.getHeight();
  unsigned int src_row_bytes = oriWidth * pymol::Image::getPixelSize();

  auto width = oriWidth / factor;
  auto height = oriHeight / factor;

  auto* p = src.bits();
  pymol::Image newImg(width, height);
  unsigned char* q = newImg.bits();
  unsigned int factor_col_bytes = factor * pymol::Image::getPixelSize();
  unsigned int factor_row_bytes = factor * src_row_bytes;

  // TODO: Rename variables and cleanup -- This was pulled almost 1:1
  for (int b = 0; b < height; b++) { /* rows */
    auto* pp = p;
    for (int a = 0; a < width; a++) { /* cols */
      unsigned int c1{};
      unsigned int c2{};
      unsigned int c3{};
      unsigned int c4{};
      auto* ppp = pp;
      for (int d = 0; d < factor; d++) { /* box rows */
        auto* pppp = ppp;
        for (int c = 0; c < factor; c++) { /* box cols */
          unsigned int alpha = pppp[3];
          c4 += alpha;
          c1 += *(pppp++) * alpha;
          c2 += *(pppp++) * alpha;
          c3 += *(pppp++) * alpha;
          pppp++;
        }
        ppp += src_row_bytes;
      }
      if (c4) { /* divide out alpha channel & average */
        c1 = c1 / c4;
        c2 = c2 / c4;
        c3 = c3 / c4;
      } else { /* alpha zero! so compute average RGB */
        c1 = c2 = c3 = 0;
        ppp = pp;
        for (int d = 0; d < factor; d++) { /* box rows */
          auto* pppp = ppp;
          for (int c = 0; c < factor; c++) { /* box cols */
            c1 += *(pppp++);
            c2 += *(pppp++);
            c3 += *(pppp++);
            pppp++;
          }
          ppp += src_row_bytes;
        }
        c1 = c1 >> shift;
        c2 = c2 >> shift;
        c3 = c3 >> shift;
      }
      *(q++) = c1;
      *(q++) = c2;
      *(q++) = c3;
      *(q++) = c4 >> shift;
      pp += factor_col_bytes;
    }
    p += factor_row_bytes;
  }
  return newImg;
}

pymol::Image GLImageToPyMOLImage(PyMOLGlobals* G,
                                 const GLFramebufferConfig& config,
                                 const Rect2D& srcRect) {
  auto imgData = G->ShaderMgr->readPixelsFrom(G, srcRect, config);
  pymol::Image img(srcRect.extent.width, srcRect.extent.height);
  if (!imgData.empty()) {
    img.setVecData(std::move(imgData));
  }
  return img;
}

void PyMOLImageCopy(const pymol::Image& srcImage, pymol::Image& dstImage,
                    const Rect2D& srcRect, const Rect2D& dstRect) {
  auto srcPx = srcImage.pixels();
  auto dstPx =
      dstImage.pixels() + (dstRect.offset.x * dstRect.extent.width) +
      (dstRect.offset.y * dstRect.extent.height) * srcRect.extent.width;
  int y_limit;
  int x_limit;

  if (((dstRect.offset.y + 1) * dstRect.extent.height) > srcRect.extent.height)
    y_limit =
        srcRect.extent.height - (dstRect.offset.y * dstRect.extent.height);
  else
    y_limit = dstRect.extent.height;
  if (((dstRect.offset.x + 1) * dstRect.extent.width) > srcRect.extent.width)
    x_limit = srcRect.extent.width - (dstRect.offset.x * dstRect.extent.width);
  else
    x_limit = dstRect.extent.width;
  for (int a = 0; a < y_limit; a++) {
    std::copy_n(srcPx, x_limit, dstPx);
    srcPx += srcRect.extent.width;
    dstPx += dstRect.extent.width;
  }
}

UpscaledExtentInfo ExtentGetUpscaleInfo(PyMOLGlobals* G, Extent2D extent,
                                        const Extent2D& maxExtent,
                                        int antialias) {
  int factor = 0;
  int shift = 0;
  if (antialias == 1) {
    factor = 2;
    shift = 2;
  }
  if (antialias >= 2) {
    factor = 4;
    shift = 4;
  }
  while (factor > 1) {
    if (((extent.width * factor) < maxExtent.width) &&
        ((extent.height * factor) < maxExtent.height)) {
      extent.width *= factor;
      extent.height *= factor;
      break;
    } else {
      factor >>= 1;
      shift -= 2;
      if (factor < 2) {
        G->Feedback->autoAdd(
            FB_Scene, FB_Blather,
            "Scene-Warning: Maximum OpenGL viewport exceeded. Antialiasing "
            "disabled.");
        break;
      }
    }
  }
  if (factor < 2) {
    factor = 0;
  }

  return UpscaledExtentInfo{extent = extent, factor = factor, shift = shift};
}

pymol::Result<> SceneMakeSizedImage(PyMOLGlobals* G, Extent2D extent,
                                    int antialias, bool excludeSelections,
                                    SceneRenderWhich renderWhich) {
  CScene* I = G->Scene;

  float sceneAspectRatio = SceneGetAspectRatio(G);

  // Calculate from implicit extents if needed
  if (extent.width == 0 && extent.height == 0) {
    extent = SceneGetExtent(G);
  } else if (extent.width != 0 && extent.height == 0) {
    extent.height = static_cast<std::uint32_t>(extent.width / sceneAspectRatio);
  } else if (extent.height != 0 && extent.width == 0) {
    extent.width = static_cast<std::uint32_t>(extent.height * sceneAspectRatio);
  }

  std::optional<Extent2D> saveExtent;
  if (!((extent.width > 0) && (extent.height > 0) && (I->Width > 0) &&
        (I->Height > 0))) {
    if (saveExtent) {
      SceneSetExtent(G, *saveExtent);
    }
    return pymol::make_error(
        "SceneMakeSizedImage-Error: invalid image dimensions");
  }

  if (!(G->HaveGUI && G->ValidContext)) {
    if (saveExtent) {
      SceneSetExtent(G, *saveExtent);
    }
    return {};
  }

  auto maxDim = SceneGLGetMaxDimensions(G);
  auto clampedExtents = ExtentClampByAspectRatio(extent, maxDim);
  auto upscaledExtentInfo =
      ExtentGetUpscaleInfo(G, clampedExtents, maxDim, antialias);
  extent = upscaledExtentInfo.extent;

  if (!saveExtent) {
    saveExtent = SceneGetExtent(G);
  }

  SceneSetExtent(G, extent);
  G->ShaderMgr->bindOffscreenSizedImage(extent, true);

  int nXStep = (extent.width / (I->Width + 1)) + 1;
  int nYStep = (extent.height / (I->Height + 1)) + 1;

  pymol::Image final_image(extent.width, extent.height);

  int total_steps = nXStep * nYStep;

  OrthoBusyPrime(G);

  /* so the trick here is that we need to move the camera around
      so that we get a pixel-perfect mosaic */
  for (int y = 0; y < nYStep; y++) {
    int y_offset = -(I->Height * y);

    for (int x = 0; x < nXStep; x++) {
      int x_offset = -(I->Width * x);

      OrthoBusyFast(G, y * nXStep + x, total_steps);

      auto offscreenFBO = G->ShaderMgr->bindOffscreenSizedImage(extent, true);

      Rect2D viewport{};
      // No offset since this is now done offscreen. In the future,
      // if there's ever a desire to do this in-place, then the
      // offset should be applied.
      viewport.extent = extent;
      SceneRenderInfo renderInfo{};
      renderInfo.mousePos = Offset2D{x_offset, y_offset};
      renderInfo.viewportOverride = viewport;
      renderInfo.excludeSelections = excludeSelections;
      renderInfo.renderWhich = renderWhich;
      renderInfo.offscreenConfig = offscreenFBO;
      G->ShaderMgr->setDrawBuffer(offscreenFBO);

      // Top-level configs are used by PP shaders.
      // Temporarily swap them out and restore.
      // TODO: This ought to be done in a stack-manner
      // in place of the GPUMgr's framegraph.
      auto oldTLConfig = G->ShaderMgr->topLevelConfig;
      G->ShaderMgr->topLevelConfig = offscreenFBO;

      // JJ: SceneSetViewport in SceneRender will extract the glViewport
      // rather than the Scene extent. Unsure why this is done, so for now
      // just preset the viewport.
      SceneSetViewport(G, viewport);
      SceneRender(G, renderInfo);

      G->ShaderMgr->topLevelConfig = oldTLConfig;

      Rect2D srcRect{};
      srcRect.extent = extent;
      auto image = GLImageToPyMOLImage(G, offscreenFBO, srcRect);

      if (!image.empty()) {
        Rect2D dstRect{{x, y}, SceneGetExtent(G)};
        PyMOLImageCopy(image, final_image, srcRect, dstRect);
      }
    }
  }

  if (!OrthoDeferredWaiting(G)) {
    if (SettingGet<int>(G, cSetting_draw_mode) == -2) {
      ExecutiveSetSettingFromString(G, cSetting_draw_mode, "-1", "", -1, true,
                                    true);
      SceneUpdate(G, false);
    }
  }

  SceneInvalidateCopy(G, true);

  if (upscaledExtentInfo.factor != 0) { /* are we oversampling? */
    final_image = SceneSupersampleImage(final_image, upscaledExtentInfo);
  }

  ScenePurgeImage(G);

  SceneSetCopyImage(G, std::move(final_image), false, true);

  if (SettingGet<bool>(G, cSetting_opaque_background)) {
    I->Image->m_needs_alpha_reset = true;
  }

  if (saveExtent) {
    SceneSetExtent(G, *saveExtent);
  }

  return {};
}

/**
 * Prepares the scene image for PNG export.
 *
 * If there is no rendered image (`G->Scene->CopyType` is false) then read the
 * current image from the OpenGL back buffer.
 *
 * Sets the alpha channel to opaque if the `opaque_background` setting is on.
 *
 * Modifies `G->Scene->Image`.
 *
 * @param prior_only Return nullptr if there is no prior image
 * (`G->Scene->Image` is nullptr)
 * @return The scene image
 */
pymol::Image* SceneImagePrepare(PyMOLGlobals* G, bool prior_only) {
  CScene* I = G->Scene;
  pymol::Image* image = nullptr;
  int save_stereo = (I->StereoMode == 1);

  if (!(I->CopyType || prior_only)) {
    if (G->HaveGUI && G->ValidContext) {
      ScenePurgeImage(G);
      I->Image =
          std::make_shared<pymol::Image>(I->Width, I->Height, save_stereo);
      image = I->Image.get();

      PyMOLReadPixels(I->rect.left, I->rect.bottom, I->Width, I->Height,
                      GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*)(image->bits()));
      I->Image->m_needs_alpha_reset = true;
    }
  } else if (I->Image) {
    image = I->Image.get();
  }
  if (image) {
    int opaque_back = SettingGetGlobal_b(G, cSetting_opaque_background);
    if (opaque_back && I->Image->m_needs_alpha_reset) {
      int i, s = image->getSizeInBytes() * (image->isStereo() ? 2 : 1);
      for (i = pymol::Image::Channel::ALPHA; i < s;
           i += pymol::Image::getPixelSize())
        image->bits()[i] = 0xFF;
      I->Image->m_needs_alpha_reset = false;
    }
  }
  return image;
}

bool ScenePNG(PyMOLGlobals* G, pymol::zstring_view png, float dpi, int quiet,
              int prior_only, int format, png_outbuf_t* outbuf) {
  CScene* I = G->Scene;
  SceneImagePrepare(G, prior_only);
  if (I->Image) {
    int width, height;
    std::tie(width, height) = I->Image->getSize();
    auto saveImage = I->Image;
    if (I->Image->isStereo()) {
      saveImage = std::make_shared<pymol::Image>();
      *(saveImage) = I->Image->interlace();
    }
    if (dpi < 0.0F) dpi = SettingGetGlobal_f(G, cSetting_image_dots_per_inch);
    auto screen_gamma = SettingGetGlobal_f(G, cSetting_png_screen_gamma);
    auto file_gamma = SettingGetGlobal_f(G, cSetting_png_file_gamma);
    // outbuf != nullptr means we're capturing to an in-memory buffer (e.g.
    // scene thumbnails) rather than a file, so there is no path to report and
    // the "check directory" hint would be misleading (issue #113).
    const bool to_buffer = (outbuf != nullptr);
    if (MyPNGWrite(png, *saveImage, dpi, format, quiet, screen_gamma,
                   file_gamma, outbuf)) {
      if (!quiet && !to_buffer) {
        PRINTFB(G, FB_Scene, FB_Actions)
        " %s: wrote %dx%d pixel image to file \"%s\".\n", __func__, width,
            I->Image->getHeight(), png.c_str() ENDFB(G);
      }
    } else if (!to_buffer) {
      PRINTFB(G, FB_Scene, FB_Errors)
      " %s-Error: error writing \"%s\"! Please check directory...\n", __func__,
          png.c_str() ENDFB(G);
    } else if (!quiet) {
      // In-memory capture failed (e.g. build without libpng). Callers that want
      // this silent (scene thumbnails) pass quiet, so don't spam the console.
      PRINTFB(G, FB_Scene, FB_Errors)
      " %s-Error: error encoding %dx%d pixel image to memory buffer.\n",
          __func__, width, I->Image->getHeight() ENDFB(G);
    }
  }
  return I->Image.get() != nullptr;
}
int SceneValidateImageMode(PyMOLGlobals* G, int mode, bool defaultdraw) {
  switch (mode) {
    case cSceneImage_Normal:
    case cSceneImage_Draw:
    case cSceneImage_Ray:
      return mode;
  }

  if (mode != cSceneImage_Default) {
    PRINTFB(G, FB_Scene, FB_Warnings)
    " %s-Warning: invalid mode %d\n", __FUNCTION__, mode ENDFB(G);
  }

  if (!G->HaveGUI || SettingGetGlobal_b(G, cSetting_ray_trace_frames)) {
    return cSceneImage_Ray;
  }

  if (defaultdraw || SettingGetGlobal_b(G, cSetting_draw_frames)) {
    return cSceneImage_Draw;
  }

  return cSceneImage_Normal;
}
int SceneMakeMovieImage(PyMOLGlobals* G, int show_timing, int validate,
                        int mode, int width, int height) {
  CScene* I = G->Scene;
  auto requestedExtent = Extent2D{static_cast<std::uint32_t>(width),
                                  static_cast<std::uint32_t>(height)};
  //  float *v;
  int valid = true;
  PRINTFB(G, FB_Scene, FB_Blather)
  " Scene: Making movie image.\n" ENDFB(G);

  // PYMOL-3209 objects inside hidden groups become visible
  ExecutiveUpdateSceneMembers(G);

  mode = SceneValidateImageMode(G, mode, width || height);

  I->DirtyFlag = false;
  switch (mode) {
    case cSceneImage_Ray:
      SceneRay(G, width, height,
               SettingGetGlobal_i(G, cSetting_ray_default_renderer), nullptr,
               nullptr, 0.0F, 0.0F, false, nullptr, show_timing, -1);
      break;
    case cSceneImage_Draw:
      SceneMakeSizedImage(G, requestedExtent,
                          SettingGet<int>(G, cSetting_antialias),
                          /*excludeSelections*/ false);
      break;
    case cSceneImage_Normal: {
      auto drawBuffer = SceneDrawBothGetConfig(G);
      if (G->HaveGUI && G->ValidContext) {
        G->ShaderMgr->setDrawBuffer(drawBuffer);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        /* insert OpenGL context validation code here? */
        SceneRenderInfo renderInfo{};
        SceneRender(G, renderInfo);
        SceneGLClearColor(0.0, 0.0, 0.0, 1.0);
        SceneCopy(G, drawBuffer, true, false);
        /* insert OpenGL context validation code here? */
      }
    } break;
  }
  MovieSetImage(G,
                MovieFrameToImage(G, SettingGetGlobal_i(G, cSetting_frame) - 1),
                I->Image);
  if (I->Image) I->CopyType = true;
  return valid;
}
int SceneLoadPNG(PyMOLGlobals* G, const char* fname, int movie_flag, int stereo,
                 int quiet) {
  CScene* I = G->Scene;
  int ok = false;
  if (I->Image) {
    ScenePurgeImage(G);
    I->CopyType = false;
    OrthoInvalidateDoDraw(
        G);  // right now, need to invalidate since text could be shown
  }
  I->Image = MyPNGRead(fname);
  if (I->Image) {
    if (!quiet) {
      PRINTFB(G, FB_Scene, FB_Details)
      " Scene: loaded image from '%s'.\n", fname ENDFB(G);
    }
    if ((stereo > 0) ||
        ((stereo < 0) && (I->Image->getWidth() == 2 * I->Width) &&
         (I->Image->getHeight() == I->Height))) {
      *(I->Image) = I->Image->deinterlace(stereo == 2);
    }

    I->CopyType = true;
    I->CopyForced = true;
    OrthoRemoveSplash(G);
    SettingSetGlobal_b(G, cSetting_text, 0);
    if (movie_flag && I->Image && !I->Image->empty()) {
      MovieSetImage(
          G, MovieFrameToImage(G, SettingGetGlobal_i(G, cSetting_frame) - 1),
          I->Image);
      I->MovieFrameFlag = true;
    } else {
      I->DirtyFlag = false; /* make sure we don't overwrite image */
    }
    OrthoDirty(G);
    ok = true;
  } else {
    if (!quiet) {
      PRINTFB(G, FB_Scene, FB_Errors)
      " Scene: unable to load image from '%s'.\n", fname ENDFB(G);
    }
  }
  return (ok);
}

// TODO: Replace with ShaderMgr::drawPixelsTo
static void RendererWritePixelsTo(PyMOLGlobals* G, const Rect2D& rect,
                                  unsigned char* buffer) {
  PyMOLDrawPixels(rect.extent.width, rect.extent.height, GL_RGBA,
                  GL_UNSIGNED_BYTE, buffer);
}

static bool SceneOverlayOversize(PyMOLGlobals* G, unsigned char* data,
                                 CGO* orthoCGO) {
  bool drawn = false;
  auto I = G->Scene;
  auto show_alpha = SettingGet<bool>(G, cSetting_show_alpha_checker);
  const float* bg_color =
      ColorGet(G, SettingGet<int>(G, nullptr, nullptr, cSetting_bg_rgb));
  unsigned int bg_rr, bg_r = (unsigned int)(255 * bg_color[0]);
  unsigned int bg_gg, bg_g = (unsigned int)(255 * bg_color[1]);
  unsigned int bg_bb, bg_b = (unsigned int)(255 * bg_color[2]);

  int factor = 1;
  int shift = 0;
  int tmp_height = I->Image->getHeight();
  int tmp_width = I->Image->getWidth();
  int src_row_bytes = I->Image->getWidth() * pymol::Image::getPixelSize();
  unsigned int color_word;
  float rgba[4] = {0.0F, 0.0F, 0.0F, 1.0F};

  ColorGetBkrdContColor(G, rgba, false);
  color_word = ColorGet32BitWord(G, rgba);

  while (tmp_height && tmp_width &&
         ((tmp_height > (I->Height - 3)) || (tmp_width > (I->Width - 3)))) {
    tmp_height = (tmp_height >> 1);
    tmp_width = (tmp_width >> 1);
    factor = (factor << 1);
    shift++;
  }
  tmp_width += 2;
  tmp_height += 2;

  if (tmp_height && tmp_width) {
    unsigned int buffer_size = tmp_height * tmp_width * 4;
    std::vector<unsigned char> buffer_vec(buffer_size);
    auto* buffer = buffer_vec.data();

    if (!buffer_vec.empty() && data) {
      unsigned char* p = data;
      unsigned char* q = buffer;
      unsigned char *pp, *ppp, *pppp;
      unsigned int c1, c2, c3, c4, tot, bg;
      unsigned int factor_col_bytes = factor * 4;
      unsigned int factor_row_bytes = factor * src_row_bytes;

      shift = shift + shift;

      for (int a = 0; a < tmp_width; a++) { /* border, first row */
        *((unsigned int*)(q)) = color_word;
        q += 4;
      }
      for (int b = 1; b < tmp_height - 1; b++) { /* rows */
        pp = p;
        *((unsigned int*)(q)) = color_word; /* border */
        q += 4;
        for (int a = 1; a < tmp_width - 1; a++) { /* cols */
          ppp = pp;

          c1 = c2 = c3 = c4 = tot = 0;

          if (show_alpha &&
              (((a >> 4) + (b >> 4)) & 0x1)) { /* introduce checkerboard */
            bg_rr = ((bg_r & 0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
            bg_gg = ((bg_g & 0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
            bg_bb = ((bg_b & 0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
          } else {
            bg_rr = bg_r;
            bg_gg = bg_g;
            bg_bb = bg_b;
          }

          for (int d = 0; d < factor; d++) { /* box rows */
            pppp = ppp;
            for (int c = 0; c < factor; c++) { /* box cols */
              unsigned char alpha = pppp[3];
              c1 += *(pppp++) * alpha;
              c2 += *(pppp++) * alpha;
              c3 += *(pppp++) * alpha;
              pppp++;
              c4 += alpha;
              tot += 0xFF;
            }
            ppp += src_row_bytes;
          }
          if (c4) {
            bg = tot - c4;
            *(q++) = (c1 + bg_rr * bg) / tot;
            *(q++) = (c2 + bg_gg * bg) / tot;
            *(q++) = (c3 + bg_bb * bg) / tot;
            *(q++) = 0xFF;
          } else {
            *(q++) = bg_rr;
            *(q++) = bg_gg;
            *(q++) = bg_bb;
            *(q++) = 0xFF;
          }
          pp += factor_col_bytes;
        }
        *((unsigned int*)(q)) = color_word; /* border */
        q += 4;
        p += factor_row_bytes;
      }
      for (int a = 0; a < tmp_width; a++) { /* border, last row */
        *((unsigned int*)(q)) = color_word;
        q += 4;
      }
      Rect2D rect{{(int)((I->Width - tmp_width) / 2 + I->rect.left),
                   (int)((I->Height - tmp_height) / 2 + I->rect.bottom)},
                  {tmp_width, tmp_height}};
      RendererWritePixelsTo(G, rect, buffer);
      drawn = true;
    }
  }
  int text_pos = (I->Height - tmp_height) / 2 - 15;
  int x_pos, y_pos;
  if (text_pos < 0) {
    text_pos = (I->Height - tmp_height) / 2 + 3;
    x_pos = (I->Width - tmp_width) / 2 + 3;
    y_pos = text_pos;
  } else {
    x_pos = (I->Width - tmp_width) / 2;
    y_pos = text_pos;
  }

  auto buffer = pymol::join_to_string("Image size = ", I->Image->getWidth(),
                                      " x ", I->Image->getHeight());

  TextSetColor3f(G, rgba[0], rgba[1], rgba[2]);
  TextDrawStrAt(G, buffer.c_str(), x_pos + I->rect.left, y_pos + I->rect.bottom,
                orthoCGO);
  return drawn;
}

static bool SceneOverlayOversizeBorder(PyMOLGlobals* G, int width, int height,
                                       unsigned char* data) {
  /* but a border around image */
  auto I = G->Scene;
  bool drawn = false;

  auto show_alpha = SettingGet<bool>(G, cSetting_show_alpha_checker);
  const float* bg_color =
      ColorGet(G, SettingGet<int>(G, nullptr, nullptr, cSetting_bg_rgb));
  unsigned int bg_rr, bg_r = (unsigned int)(255 * bg_color[0]);
  unsigned int bg_gg, bg_g = (unsigned int)(255 * bg_color[1]);
  unsigned int bg_bb, bg_b = (unsigned int)(255 * bg_color[2]);

  unsigned int color_word;
  float rgba[4] = {0.0F, 0.0F, 0.0F, 1.0F};
  unsigned int tmp_height = height + 2;
  unsigned int tmp_width = width + 2;
  unsigned int border = 1;
  unsigned int upscale = 1;

  // Upscale for Retina/4K
  if (DIP2PIXEL(height) == I->Height && DIP2PIXEL(width) == I->Width) {
    upscale = DIP2PIXEL(1);
    tmp_height = DIP2PIXEL(height);
    tmp_width = DIP2PIXEL(width);
    border = 0;
  }

  unsigned int n_word = tmp_height * tmp_width;
  std::vector<unsigned int> tmp_buffer_vec(n_word);
  ColorGetBkrdContColor(G, rgba, false);
  color_word = ColorGet32BitWord(G, rgba);

  if (!tmp_buffer_vec.empty()) {
    auto* tmp_buffer = tmp_buffer_vec.data();
    unsigned int a, b;
    unsigned int* p = (unsigned int*)data;
    unsigned int* q = tmp_buffer;

    // top border
    for (a = 0; a < border; ++a) {
      for (b = 0; b < tmp_width; b++) *(q++) = color_word;
    }

    for (a = border; a < tmp_height - border; a++) {
      // left border
      for (b = 0; b < border; ++b) {
        *(q++) = color_word;
      }

      for (b = border; b < tmp_width - border; b++) {
        unsigned char* qq = (unsigned char*)q;
        unsigned char* pp = (unsigned char*)p;
        unsigned char bg;
        if (show_alpha &&
            (((a >> 4) + (b >> 4)) & 0x1)) { /* introduce checkerboard */
          bg_rr = ((bg_r & 0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
          bg_gg = ((bg_g & 0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
          bg_bb = ((bg_b & 0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
        } else {
          bg_rr = bg_r;
          bg_gg = bg_g;
          bg_bb = bg_b;
        }
        if (pp[3]) {
          bg = 0xFF - pp[3];
          *(qq++) = (pp[0] * pp[3] + bg_rr * bg) / 0xFF;
          *(qq++) = (pp[1] * pp[3] + bg_gg * bg) / 0xFF;
          *(qq++) = (pp[2] * pp[3] + bg_bb * bg) / 0xFF;
          *(qq++) = 0xFF;
        } else {
          *(qq++) = bg_rr;
          *(qq++) = bg_gg;
          *(qq++) = bg_bb;
          *(qq++) = 0xFF;
        }
        q++;

        if ((b + 1 - border) % upscale == 0) {
          p++;
        }
      }

      if ((a + 1 - border) % upscale != 0) {
        // read row again
        p -= width;
      }

      // right border
      for (b = 0; b < border; ++b) {
        *(q++) = color_word;
      }
    }

    // bottom border
    for (a = 0; a < border; ++a) {
      for (b = 0; b < tmp_width; b++) *(q++) = color_word;
    }

    Rect2D rect{{(int)((I->Width - tmp_width) / 2 + I->rect.left),
                 (int)((I->Height - tmp_height) / 2 + I->rect.bottom)},
                {tmp_width, tmp_height}};
    RendererWritePixelsTo(G, rect, (unsigned char*)tmp_buffer);
    drawn = true;
  }
  return drawn;
}

static bool SceneOverlayExactFitNoAlpha(PyMOLGlobals* G, int width, int height,
                                        unsigned char* data) {
  auto I = G->Scene;
  Rect2D rect{{(int)((I->Width - width) / 2 + I->rect.left),
               (int)((I->Height - height) / 2 + I->rect.bottom)},
              {width, height}};
  RendererWritePixelsTo(G, rect, data);
  return true;
}

static bool SceneOverlayExactFit(PyMOLGlobals* G, int width, int height,
                                 unsigned char* data) {
  auto I = G->Scene;
  float rgba[4] = {0.0F, 0.0F, 0.0F, 1.0F};
  unsigned int n_word = height * width;
  std::vector<unsigned int> tmp_buffer_vec(n_word);
  ColorGetBkrdContColor(G, rgba, false);

  auto show_alpha = SettingGet<bool>(G, cSetting_show_alpha_checker);
  const float* bg_color =
      ColorGet(G, SettingGet<int>(G, nullptr, nullptr, cSetting_bg_rgb));
  unsigned int bg_rr, bg_r = (unsigned int)(255 * bg_color[0]);
  unsigned int bg_gg, bg_g = (unsigned int)(255 * bg_color[1]);
  unsigned int bg_bb, bg_b = (unsigned int)(255 * bg_color[2]);

  if (tmp_buffer_vec.empty()) {
    return false;
  }
  auto* tmp_buffer = tmp_buffer_vec.data();
  unsigned int a, b;
  unsigned int* p = (unsigned int*)data;
  unsigned int* q = tmp_buffer;
  for (a = 0; a < (unsigned int)height; a++) {
    for (b = 0; b < (unsigned int)width; b++) {
      unsigned char* qq = (unsigned char*)q;
      unsigned char* pp = (unsigned char*)p;
      unsigned char bg;
      if (show_alpha &&
          (((a >> 4) + (b >> 4)) & 0x1)) { /* introduce checkerboard */
        bg_rr = ((bg_r & 0x80) ? bg_r - TRN_BKG : bg_r + TRN_BKG);
        bg_gg = ((bg_g & 0x80) ? bg_g - TRN_BKG : bg_g + TRN_BKG);
        bg_bb = ((bg_b & 0x80) ? bg_b - TRN_BKG : bg_b + TRN_BKG);
      } else {
        bg_rr = bg_r;
        bg_gg = bg_g;
        bg_bb = bg_b;
      }
      if (pp[3]) {
        bg = 0xFF - pp[3];
        *(qq++) = (pp[0] * pp[3] + bg_rr * bg) / 0xFF;
        *(qq++) = (pp[1] * pp[3] + bg_gg * bg) / 0xFF;
        *(qq++) = (pp[2] * pp[3] + bg_bb * bg) / 0xFF;
        *(qq++) = 0xFF;
      } else {
        *(qq++) = bg_rr;
        *(qq++) = bg_gg;
        *(qq++) = bg_bb;
        *(qq++) = 0xFF;
      }
      q++;
      p++;
    }
  }
  Rect2D rect{{(int)((I->Width - width) / 2 + I->rect.left),
               (int)((I->Height - height) / 2 + I->rect.bottom)},
              {width, height}};
  RendererWritePixelsTo(G, rect, (unsigned char*)tmp_buffer);
  return true;
}

int SceneDrawImageOverlay(PyMOLGlobals* G, int override, CGO* orthoCGO) {
  CScene* I = G->Scene;
  int drawn = false;
  int text = SettingGet<bool>(G, cSetting_text);
  /* is the text/overlay (ESC) on? */
  int overlay = OrthoGetOverlayStatus(G);
  bool draw_overlay = (!text || overlay) && (override || I->CopyType == true) &&
                      I->Image && !I->Image->empty();

  if (!draw_overlay) {
    return drawn;
  }

  int width = I->Image->getWidth();
  int height = I->Image->getHeight();
  unsigned char* data = I->Image->bits();

  if ((height > I->Height) || (width > I->Width)) { /* image is oversize */
    drawn = SceneOverlayOversize(G, data, orthoCGO);
  } else if (((width < I->Width) || (height < I->Height)) &&
             ((I->Width - width) > 2) && ((I->Height - height) > 2)) {
    drawn = SceneOverlayOversizeBorder(G, width, height, data);
  } else if (I->CopyForced) { /* near-exact fit */
    drawn = SceneOverlayExactFit(G, width, height, data);
  } else { /* not a forced copy, so don't show/blend alpha */
    drawn = SceneOverlayExactFitNoAlpha(G, width, height, data);
  }

  I->LastRender = UtilGetSeconds(G);
  return drawn;
}

/**
 * Will call cmd.raw_image_callback(img) with the current RGBA image, copied
 * to a WxHx4 numpy array. Return false if no callback is defined
 * (cmd.raw_image_callback == None).
 */
static bool call_raw_image_callback(PyMOLGlobals* G) {
  bool done = false;

  int blocked = PAutoBlock(G);
  auto raw_image_callback =
      PyObject_GetAttrString(G->P_inst->cmd, "raw_image_callback");

  if (raw_image_callback != Py_None) {
    PRINTFB(G, FB_Scene, FB_Errors)
    " raw_image_callback-Error: no numpy support\n" ENDFB(G);
  }

  Py_XDECREF(raw_image_callback);
  PAutoUnblock(G, blocked);

  return done;
}

/**
 * Creates an image of the current scene.
 *
 * @param extent requested extent
 * @param antialias antialias mode
 * @param dpi dots per inch
 * @param format ??
 * @param quiet if true, don't print messages
 * @param out_img if not nullptr, store the image data here
 * @param filename if not empty, store the image to this file as PNG
 */

static void SceneImage(PyMOLGlobals* G, const Extent2D& extent, int antialias,
                       float dpi, int format, bool quiet, pymol::Image* out_img,
                       const std::string& filename) {
  auto allButGizmos_i = pymol::to_underlying(SceneRenderWhich::All) &
                        ~pymol::to_underlying(SceneRenderWhich::Gizmos);
  auto allButGizmos = static_cast<SceneRenderWhich>(allButGizmos_i);
  SceneMakeSizedImage(G, extent, antialias, /*excludeSelecions*/ true,
                      allButGizmos);
  if (!filename.empty()) {
    ScenePNG(G, filename.c_str(), dpi, quiet, false, format, nullptr);
  } else if (out_img) {
    png_outbuf_t outbuf;
    ScenePNG(G, "", dpi, quiet, false, format, &outbuf);
    out_img->setVecData(std::move(outbuf));
  } else if (call_raw_image_callback(G)) {
  } else if (G->HaveGUI && SettingGetGlobal_b(G, cSetting_auto_copy_images)) {
  }
}

bool SceneDeferImage(PyMOLGlobals* G, const Extent2D& extent,
                     const char* filename, int antialias, float dpi, int format,
                     int quiet, pymol::Image* out_img) {
  std::string filenameStr = filename ? filename : "";
  std::function<void()> deferred = [=]() {
    SceneImage(G, extent, antialias, dpi, format, quiet, out_img, filenameStr);
  };

  if (G->ValidContext) {
    deferred();
    return false;
  }

  OrthoDefer(G, std::move(deferred));
  return true;
}
void SceneApplyImageGamma(PyMOLGlobals* G, unsigned int* buffer, int width,
                          int height) {
  float gamma = SettingGetGlobal_f(G, cSetting_gamma);

  if (gamma > R_SMALL4)
    gamma = 1.0F / gamma;
  else
    gamma = 1.0F;

  if (buffer && height && width) {
    float _inv3 = 1 / (255 * 3.0F);
    float _1 = 1 / 3.0F;
    unsigned char* p;
    int x, y;
    float c1, c2, c3, inp, sig;
    unsigned int i1, i2, i3;
    p = (unsigned char*)buffer;
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        c1 = p[0];
        c2 = p[1];
        c3 = p[2];
        inp = (c1 + c2 + c3) * _inv3;
        if (inp < R_SMALL4)
          sig = _1;
        else
          sig = (float)(pow(inp, gamma) / inp);
        i1 = (unsigned int)(sig * c1);
        i2 = (unsigned int)(sig * c2);
        i3 = (unsigned int)(sig * c3);
        if (i1 > 255) i1 = 255;
        if (i2 > 255) i2 = 255;
        if (i3 > 255) i3 = 255;
        p[0] = i1;
        p[1] = i2;
        p[2] = i3;
        p += 4;
      }
    }
  }
}
