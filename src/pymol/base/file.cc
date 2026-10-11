/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO.
 */

#include <stdio.h>
#include <stdlib.h>

#include "pymol/base/file.h"
#include "pymol/base/file_stream.h"
#include "pymol/platform/memory_debug.h"

#include <fstream>
#include <string>

#include "pymol/zstring_view.h"

/**
 * Get the size from the current file pointer to the end of the file
 */
static long fgetsize(FILE* fp) {
  long filesize, current = ftell(fp);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  fseek(fp, current, SEEK_SET);
  return filesize;
}

/**
 * Allocate memory and read the entire file from the given file pointer.
 * The file size is stored into the size pointer if not nullptr.
 */
static char* fgetcontents(FILE* fp, long* size) {
  long filesize = fgetsize(fp);

  char* contents = pymol::malloc<char>(filesize + 255);
  if (!contents) return nullptr;

  if (1 != fread(contents, filesize, 1, fp)) {
    mfree(contents);
    return nullptr;
  }

  if (size) *size = filesize;

  contents[filesize] = '\0';
  return contents;
}

/**
 * Allocate memory and read the entire file for the given filename.
 * The file size is stored into the size pointer if not nullptr.
 */
char* FileGetContents(const char* filename, long* size) {
  char* contents;
  FILE* fp = pymol_fopen(filename, "rb");

  if (!fp) return nullptr;

  contents = fgetcontents(fp, size);
  fclose(fp);
  return contents;
}

namespace pymol {

static long istream_get_size(std::istream& file) {
  auto current = file.tellg();
  file.seekg(0, std::ios::end);
  auto filesize = file.tellg();
  file.seekg(current);
  return filesize;
}

static std::string istream_get_contents(std::istream& file) {
  auto filesize = istream_get_size(file);
  std::string contents(filesize, '\0');
  file.read(&contents[0], filesize);
  return contents;
}

std::string file_get_contents(pymol::zstring_view filename) {
  std::ifstream file;
  fstream_open(file, filename, std::ios::in | std::ios::binary);
  return istream_get_contents(file);
}

}  // namespace pymol
