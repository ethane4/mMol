#pragma once

#include "pymol/zstring_view.h"

struct PyMOLGlobals;
class CGO;
struct RenderInfo;
struct CSetting;
struct Rep;
struct PickContext;

void CGORenderGL(CGO* I, const float* color, CSetting* set1, CSetting* set2,
                 RenderInfo* info, Rep* rep);
void CGORenderGLAlpha(CGO* I, RenderInfo* info, bool calcDepth);
