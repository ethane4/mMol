#include "pymol/platform/main_thread.h"

#include "pymol/platform/os_python.h"

static long s_main_thread_id = -1;

void PSetGlutThread() { s_main_thread_id = PyThread_get_thread_ident(); }

int PIsGlutThread() { return PyThread_get_thread_ident() == s_main_thread_id; }
