
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
#ifndef _H_main
#define _H_main

int MainSavingUnderWhileIdle(void);

#include "pymol/platform/os_python.h"
#include "pymol/base/pymol_globals.h"

PyObject* MainAsPyList(PyMOLGlobals* G);
int MainFromPyList(PyMOLGlobals* G, PyObject* list);

#define PYMOL_MAX_OPT_STR 1025

int main_shared(int);

#endif /* _H_main header */
