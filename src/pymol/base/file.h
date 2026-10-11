/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO.
 */

#ifndef _H_File
#define _H_File

#define pymol_fopen fopen

char* FileGetContents(const char* filename, long* size);

#endif
