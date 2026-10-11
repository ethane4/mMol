
/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_Font
#define _H_Font

#include "pymol/base/pymol_globals.h"
#include "pymol/base/base.h"

#define cFontNULL 0
#define cFontGLUT 1
#define cFontBitmap 2
#define cFontVector 3
#define cFontExtrude 4
#define cFontPixmap 5
#define cFontTexture 6

#define cFontModePixels 0
#define cFontModeSpacial 1

#include "pymol/platform/os_gl.h"

struct CFont {
  PyMOLGlobals* G = nullptr;
  int TextID = 0;
  virtual const char* RenderOpenGL(const RenderInfo* info, const char* text,
                                   float size, const float* rpos, bool needSize,
                                   short relativeMode, bool shouldRender,
                                   CGO* shaderCGO) {
    return nullptr;
  }
  virtual const char* RenderOpenGLFlat(const RenderInfo* info, const char* text,
                                       float size, const float* rpos,
                                       bool needSize, short relativeMode,
                                       bool shouldRender, CGO* shaderCGO) {
    return nullptr;
  }
  virtual const char* RenderRay(CRay* ray, const char* text, float size,
                                const float* rpos, bool needSize,
                                short relativeMode) {
    return nullptr;
  }
  CFont(PyMOLGlobals* G) : G(G) {};
  virtual ~CFont() = 0;
};

/// Glyph drawing shared by the bitmap and FreeType fonts.
void CharacterRenderOpenGLPrime(PyMOLGlobals* G, const RenderInfo* info);
short CharacterRenderOpenGL(PyMOLGlobals* G, const RenderInfo* info, int id,
                            short isworldlabel, short relativeMode,
                            CGO* shaderCGO);
void CharacterRenderOpenGLDone(PyMOLGlobals* G, const RenderInfo* info);

#endif
