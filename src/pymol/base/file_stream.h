/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO with C++ streams.
 */

#pragma once

#include <fstream>
#include <string>

#include "pymol/zstring_view.h"

namespace pymol {
/**
 * Reads entire file into a string
 * @param filename Path in native filesystem encoding or UTF-8
 * @throw ... If file cannot be opened
 */
std::string file_get_contents(pymol::zstring_view filename);

/**
 * File stream open wrapper with UTF-8 support on Windows. On Unix, this simply
 * calls `stream.open`.
 *
 * @tparam T std::fstream, std::ifstream, std::ofstream, or a subclass
 * @param stream Stream to open
 * @param filename Path in native filesystem encoding or UTF-8
 * @param mode see std::fstream::open
 * @throw ... If file cannot be opened
 */
template <typename T>
void fstream_open(T& stream, zstring_view filename,
                  std::ios_base::openmode mode) {
  try {
    stream.open(filename.data(), mode);
  } catch (const std::ios_base::failure&) {
    throw;
  }
}
}  // namespace pymol
