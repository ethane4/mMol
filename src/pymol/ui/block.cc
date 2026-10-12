

/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/
#include "pymol/platform/os_python.h"
#include "pymol/platform/os_predef.h"
#include "pymol/platform/os_gl.h"

#include "pymol/ui/block.h"
#include "pymol/app/main.h"
#include "pymol/graphics/cgo.h"
#include "pymol/graphics/immediate_helper.h"

int Block::getWidth() const { return rect.right - rect.left; }
int Block::getHeight() const { return rect.top - rect.bottom; }
void Block::fill(CGO* orthoCGO) {
  if (m_G->HaveGUI && m_G->ValidContext) {
    if (orthoCGO) {
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.right, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.bottom, 0.f);
      CGOEnd(orthoCGO);
    } else {
      ImmBatch batch;
      batch.begin(GL_POLYGON);
      batch.vertex2i(rect.right, rect.top);
      batch.vertex2i(rect.right, rect.bottom);
      batch.vertex2i(rect.left, rect.bottom);
      batch.vertex2i(rect.left, rect.top);
      batch.end();
    }
  }
}
void Block::drawLeftEdge(CGO* orthoCGO) {
  if (m_G->HaveGUI && m_G->ValidContext) {
    if (orthoCGO) {
      CGOColor(orthoCGO, .3f, .3f, .3f);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left + 1.f, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left + 1.f, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      if (m_G->HaveGUI && m_G->ValidContext) {
        ImmBatch batch;
        batch.begin(GL_LINES);
        batch.color3f(0.3, 0.3, 0.3);
        batch.vertex2i(rect.left, rect.bottom);
        batch.vertex2i(rect.left, rect.top);
        batch.end();
      }
    }
  }
}
void Block::drawTopEdge() {}
void Block::setMargin(int t, int l, int b, int r) {
  margin.top = t;
  margin.left = l;
  margin.bottom = b;
  margin.right = r;
}
void Block::reshape(int width, int height) {
  rect.top = (height - margin.top);
  rect.left = margin.left;
  rect.bottom = margin.bottom;
  rect.right = (width - margin.right);
}
void Block::translate(int dx, int dy) {
  rect.top += dy;
  rect.left += dx;
  rect.bottom += dy;
  rect.right += dx;
}
void Block::recursiveDraw(CGO* orthoCGO) {
  if (active) {
    draw(orthoCGO);
  }
}

Block* Block::recursiveFind(int x, int y) {
  if (!active || !rectXYInside(x, y)) {
    return nullptr;
  }
  return this;
}

bool Block::rectXYInside(int x, int y) const {
  return ((y <= rect.top) && (y >= rect.bottom) && (x <= rect.right) &&
          (x >= rect.left));
}
