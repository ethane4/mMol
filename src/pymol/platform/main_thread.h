// Which thread is PyMOL's main (GUI) thread.

#pragma once

/// Records the calling thread as the main thread (at startup).
void PSetGlutThread();

/// True on the main thread.
int PIsGlutThread();
