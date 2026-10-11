

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
#ifndef _H_os_predef
#define _H_os_predef

/* Macros used by Fortify source in GCC 4.1.x are incompatible with
   PyMOL's Feedback system... */

#ifdef _FORTIFY_SOURCE
#undef _FORTIFY_SOURCE
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */

#define PATH_SEP "/"

/* commercial product */

/* collaboration product (placarded) */

/* educational product (placarded) */

/* evaluation product (placarded) */

/* END PROPRIETARY CODE SEGMENT */

#include <stddef.h>

#include "ov/ov_types.h"

// alternative to std::swap if references are not allowed (e.g. bit fields)
#define SWAP_NOREF(a, b) \
  {                      \
    auto _t = (a);       \
    (a) = (b);           \
    (b) = _t;            \
  }

#endif
