#pragma once

#include <string>

#include "pymol/base/rect.h"

struct SceneElem {
  std::string name;
  pymol::Rect<int> rect{};
  bool drawn;
  SceneElem(std::string elem_, bool drawn);
};
