#ifndef _H_OVHeap
#define _H_OVHeap

#include "ov/ov_port.h"
#include "ov/ov_types.h"

/* OVHeap is a macro wrapper around stdlib's memory management routines */

typedef int* OVHeap;

#define OVHeap_New() ((void*)1)
#define OVHeap_Del(I)

#define OVHeap_Malloc(I, size) malloc(size)
#define OVHeap_Calloc(I, num, size) calloc(num, size)
#define OVHeap_Realloc(I, ptr, size) realloc(ptr, size)
#define OVHeap_Free(I, ptr) free(ptr)

#define OVHeap_Dump(I, log_to_files)

/* convenience macros */

#define OVHeap_ALLOC(I, type) (type*)OVHeap_Calloc(I, 1, sizeof(type))

#define OVHeap_MALLOC(I, type, num) \
  (type*)OVHeap_Malloc(I, sizeof(type) * (num))
#define OVHeap_CALLOC(I, type, num) (type*)OVHeap_Calloc(I, num, sizeof(type))
#define OVHeap_REALLOC(I, ptr, type, num) \
  OVHeap_Realloc(I, ptr, sizeof(type) * (num))
#define OVHeap_FREE_AUTO_NULL(I, ptr) \
  {                                   \
    if (ptr) {                        \
      OVHeap_Free(I, ptr);            \
      ptr = NULL;                     \
    }                                 \
  }

#endif
