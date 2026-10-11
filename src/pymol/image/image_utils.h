#pragma once

#include "pymol/image/image.h"
#include "pymol/base/result.h"

struct PyMOLGlobals;

/**
 * @brief Composites an image with an overlay
 * @param image the base image
 * @param overlay the overlay image
 * @return the composited image
 */
pymol::Result<pymol::Image> PyMOLImageComposite(PyMOLGlobals* G,
                                                const pymol::Image& image,
                                                const pymol::Image& overlay);
