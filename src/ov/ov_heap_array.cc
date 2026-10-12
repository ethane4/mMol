#include "ov/ov_heap_array.h"

#include "ov/ov_utility.h"

void* _OVHeapArray_Check(void* ptr, ov_size index) { /* array-growth routine */
  _OVHeapArray* vla = &(((_OVHeapArray*)ptr)[-1]);
  if (index >= vla->size) {
    ov_size new_size = (index + (index >> 1) + 1);
    _OVHeapArray* new_vla;
    new_vla = (_OVHeapArray*)OVHeap_Realloc(
        vla->heap, vla, sizeof(_OVHeapArray) + (vla->unit_size * new_size));
    if (!new_vla) {
      fprintf(stderr, "_OVHeapArray_Check-Error: realloc failed\n");
    } else {
      vla = new_vla;
      if (vla->auto_zero) {
        char* start =
            ((char*)vla) + sizeof(_OVHeapArray) + (vla->unit_size * vla->size);
        char* stop =
            ((char*)vla) + sizeof(_OVHeapArray) + (vla->unit_size * new_size);
        ov_utility_zero_range(start, stop);
      }
      vla->size = new_size;
    }
  }
  return ((void*)&(vla[1]));
}

void* _OVHeapArray_Alloc(OVHeap* heap, ov_size unit_size, ov_size size,
                         int zero) {
  _OVHeapArray* vla;
  if (zero) {
    vla = (_OVHeapArray*)OVHeap_Calloc(
        heap, 1, sizeof(_OVHeapArray) + (unit_size * size));
  } else {
    vla = (_OVHeapArray*)OVHeap_Malloc(
        heap, sizeof(_OVHeapArray) + (unit_size * size));
  }
  if (!vla) {
    fprintf(stderr, "_OVHeapArray: realloc failed\n");
    return OV_NULL;
  } else {
    vla->heap = heap;
    vla->size = size;
    vla->unit_size = unit_size;
    vla->auto_zero = zero;
    return ((void*)&(vla[1]));
  }
}

void _OVHeapArray_Free(void* ptr) {
  if (!ptr) {
    fprintf(stderr, "_OVHeapArray_Free-Error: tried to free NULL pointer!\n");
  } else {
    _OVHeapArray* vla = &(((_OVHeapArray*)ptr)[-1]);
    OVHeap_Free(vla->heap, vla);
  }
}

ov_size OVHeapArray_GetSize(void* ptr) {
  _OVHeapArray* vla;
  vla = &((_OVHeapArray*)ptr)[-1];
  return (vla->size);
}

void* _OVHeapArray_SetSize(void* ptr, ov_size new_size) {
  _OVHeapArray *vla, *new_vla;
  vla = &((_OVHeapArray*)ptr)[-1];
  new_vla = (_OVHeapArray*)OVHeap_Realloc(
      vla->heap, vla, (vla->unit_size * new_size) + sizeof(_OVHeapArray));
  if (!new_vla) {
    fprintf(stderr, "VLASetSize-ERR: realloc failed.\n");
  } else {
    vla = new_vla;
    if (new_size > vla->size && vla->auto_zero) {
      char* start =
          ((char*)vla) + sizeof(_OVHeapArray) + (vla->unit_size * vla->size);
      char* stop =
          ((char*)vla) + sizeof(_OVHeapArray) + (vla->unit_size * new_size);
      ov_utility_zero_range(start, stop);
    }
    vla->size = new_size;
  }
  return ((void*)&(vla[1]));
}
