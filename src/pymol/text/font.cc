

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

#include "pymol/base/util.h"
#include "pymol/text/font.h"

#include <glm/gtc/type_ptr.hpp>

#include "pymol/geometry/vector.h"
#include "pymol/glyph/character.h"
#include "pymol/graphics/cgo.h"
#include "pymol/graphics/texture.h"
#include "pymol/platform/os_gl.h"
#include "pymol/settings/setting.h"
#include "pymol/text/text.h"

CFont::~CFont() = default;

void CharacterRenderOpenGLPrime(PyMOLGlobals* G, const RenderInfo* info) {
  if (G->HaveGUI && G->ValidContext) {
    if ((info && !info->use_shaders) ||
        (!info && !SettingGetGlobal_b(G, cSetting_use_shaders))) {
      glEnable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    }
    /*    glEnable(GL_BLEND);
       glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); */
  }
}

void CharacterRenderOpenGLDone(PyMOLGlobals* G, const RenderInfo* info) {
  if (G->HaveGUI && G->ValidContext) {
    if ((info && !info->use_shaders) ||
        (!info && !SettingGetGlobal_b(G, cSetting_use_shaders))) {
      glDisable(GL_TEXTURE_2D);
    }
    /*    glDisable(GL_BLEND); */
  }
}

short CharacterRenderOpenGL(PyMOLGlobals* G, const RenderInfo* info, int id,
                            short isworldlabel, short relativeMode,
                            CGO* shaderCGO)

/* need orientation matrix */
{
  CCharacter* I = G->Character;
  CharRec* rec = I->Char + id;
  short success = 1;
  auto isTextured = TextureIsCharTextured(G, id, rec->extent);
  float sampling = 1.0F;

  if (G->HaveGUI && G->ValidContext && isTextured) {
    if (info) sampling = (float)info->sampling;
    if (isTextured) {
      /*    if(glIsTexture(texture_id)) -- BAD -- impacts performance */
      float *v, v0[3];
      float v1[3];
      if (!shaderCGO) {
        TextureBindTexture(G);
      }
      v = TextGetPos(G);
      copy3f(v, v0);
      v0[0] -= rec->XOrig / sampling;
      v0[1] -= rec->YOrig / sampling;
      copy3f(v0, v1);
      v1[0] += rec->Width / sampling;
      v1[1] += rec->Height / sampling;
      /*      glColor4f(0.5F,0.5F,0.5F,1.0F); */
      if (shaderCGO) {
        float* worldPos = TextGetWorldPos(G);
        if (isworldlabel) {
          float* targetPos = TextGetTargetPos(G);
          float* screenWorldOffset = TextGetScreenWorldOffset(G);
          shaderCGO->add<cgo::draw::label>(
              glm::make_vec3(worldPos), glm::make_vec3(screenWorldOffset),
              glm::make_vec3(v0), glm::make_vec3(v1),
              glm::make_vec4(rec->extent), relativeMode,
              glm::make_vec3(targetPos));
        } else {
          CGODrawTexture(shaderCGO, worldPos, v0, v1, rec->extent);
        }
      } else {
      }
    }
    TextAdvance(G, rec->Advance / sampling);
  } else {
    if (!isTextured) success = 0;
  }
  return success;
}
