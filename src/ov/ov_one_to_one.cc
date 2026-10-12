
#include "ov/ov_one_to_one.h"
#include "ov/ov_heap_array.h"
#include "ov/ov_utility.h"

#define HASH(value, mask) \
  (((value ^ (value >> 24)) ^ ((value >> 8) ^ (value >> 16))) & mask)

/* FYI: "up" stands for UniquePair -- a precursor to OneToOne */

namespace {
typedef struct {
  int active;
  ov_word forward_value, reverse_value;
  ov_size forward_next, reverse_next;
} up_element;
}  // namespace

struct _OVOneToOne {
  OVHeap* heap;
  ov_uword mask;
  ov_size size, n_inactive;
  ov_word next_inactive;
  up_element* elem;
  ov_word* forward;
  ov_word* reverse;
};

OVstatus OVOneToOne_Init(OVOneToOne* up, OVHeap* heap) {
  ov_utility_zero_range(up, up + 1);
  up->heap = heap;
  return_OVstatus_SUCCESS;
}

OVOneToOne* OVOneToOne_New(OVHeap* heap) {
  OVOneToOne* up;
  up = OVHeap_ALLOC(heap, OVOneToOne);
  up->heap = heap;
  return up;
}

void OVOneToOne_Purge(OVOneToOne* up) {
  if (up) {
    OVHeapArray_FREE_AUTO_NULL(up->elem);
    OVHeap_FREE_AUTO_NULL(up->heap, up->forward);
    OVHeap_FREE_AUTO_NULL(up->heap, up->reverse);
  }
}

void OVOneToOne_Del(OVOneToOne* up) {
  if (up) {
    OVOneToOne_Purge(up);
    OVHeap_FREE_AUTO_NULL(up->heap, up);
  }
}

void OVOneToOne_Reset(OVOneToOne* up) {
  OVOneToOne_Purge(up);
  OVOneToOne_Init(up, up->heap);
}

static void Reload(OVOneToOne* up) { /* assumes hash tables are clean and
                                        initialized to zero */
  ov_uword mask = up->mask;

  if (up->elem && mask) {
    {
      up_element* elem = up->elem;
      ov_uword a;
      for (a = 0; a < up->size; a++) {
        if (elem->active) {
          elem->forward_next = 0; /* 0 is the sentinel for end of list */
          elem->reverse_next = 0;
        }
        elem++;
      }
    }

    {
      ov_uword a;
      ov_word* forward = up->forward;
      ov_word* reverse = up->reverse;
      up_element* elem = up->elem;
      {
        ov_word fwd, rev;
        ov_word fwd_val;
        ov_word rev_val;
        for (a = 0; a < up->size; a++) {
          if (elem->active) {
            fwd_val = elem->forward_value;
            rev_val = elem->reverse_value;
            fwd = HASH(fwd_val, mask);
            rev = HASH(rev_val, mask);
            elem->forward_next = forward[fwd];
            forward[fwd] = a + 1; /* NOTE: 1 based indices */
            elem->reverse_next = reverse[rev];
            reverse[rev] = a + 1;
          }
          elem++;
        }
      }
    }
  }
}

OVreturn_word OVOneToOne_GetForward(OVOneToOne* up, ov_word forward_value) {
  if (!up) {
    OVreturn_word result = {OVstatus_NULL_PTR};
    return result;
  } else {
    ov_uword mask = up->mask;
    if (mask) {
      ov_word hash = HASH(forward_value, mask);
      up_element* elem = up->elem;
      ov_word index = up->forward[hash];
      up_element* cur_elem = elem + (index - 1);
      while (index) {
        if (cur_elem->forward_value == forward_value) {
          OVreturn_word result = {OVstatus_SUCCESS};
          result.word = cur_elem->reverse_value;
          return result;
        }
        index = cur_elem->forward_next;
        cur_elem = elem + (index - 1);
      }
    }
    {
      OVreturn_word result = {OVstatus_NOT_FOUND};
      return result;
    }
  }
}

static OVstatus Recondition(OVOneToOne* up, ov_uword size, int force) {
  if (!up) {
    return_OVstatus_NULL_PTR;
  } else {
    ov_uword mask = up->mask;
    if ((size > mask) || ((size << 2) < mask) || force) {
      while ((size << 2) < mask) {
        mask = mask >> 1;
        if (mask < 2) break;
      }

      while (size > mask) {
        mask = (mask << 1) + 1;
      }

      {
        if (!up->elem) {
          up->elem = OVHeapArray_CALLOC(up->heap, up_element, size);
          if (!up->elem) {
            return_OVstatus_OUT_OF_MEMORY;
          }
        }
        if (mask != up->mask) {
          ov_word* tmp_forward = OVHeap_CALLOC(up->heap, ov_word, mask + 1);
          ov_word* tmp_reverse = OVHeap_CALLOC(up->heap, ov_word, mask + 1);
          if (!(tmp_forward && tmp_reverse)) { /* validate */
            OVHeap_FREE_AUTO_NULL(up->heap, tmp_forward);
            OVHeap_FREE_AUTO_NULL(up->heap, tmp_reverse);
            /* being unable to condition is not an error */
          } else {
            /* impossible to fail after here... */
            OVHeap_FREE_AUTO_NULL(up->heap, up->forward);
            OVHeap_FREE_AUTO_NULL(up->heap, up->reverse);
            up->forward = tmp_forward;
            up->reverse = tmp_reverse;
            up->mask = mask;
          }
        } else {
          ov_utility_zero_range(up->forward, up->forward + (up->mask + 1));
          ov_utility_zero_range(up->reverse, up->reverse + (up->mask + 1));
        }
        Reload(up);
      }
    }
  }
  return_OVstatus_SUCCESS;
}

OVstatus OVOneToOne_Pack(OVOneToOne* up) {
  if (!up) {
    return_OVstatus_NULL_PTR;
  } else {
    if (up->n_inactive && up->elem) {
      ov_uword new_size = 0;
      up_element *src = up->elem, *dst = up->elem;
      ov_uword a;

      for (a = 0; a < up->size; a++) {
        if (src->active) {
          if (src > dst) {
            *dst = *src;
          }
          dst++;
          new_size++;
        }
        src++;
      }
      up->n_inactive = 0;
      up->next_inactive = 0;
      if (new_size > 0 && new_size < up->size) {
        if (!OVHeapArray_SET_SIZE(up->elem, up_element, new_size))
          ov_utility_zero_range(up->elem + new_size, up->elem + up->size);
      }
      up->size = new_size;
      return Recondition(up, new_size, OV_TRUE);
    }
    return_OVstatus_SUCCESS;
  }
}

OVstatus OVOneToOne_DelReverse(OVOneToOne* up, ov_word reverse_value) {
  if (!up) {
    return_OVstatus_NULL_PTR;
  } else {
    ov_word mask = up->mask;
    if (mask) {
      ov_word rev_hash = HASH(reverse_value, mask);
      ov_word rev = up->reverse[rev_hash];

      if (!rev) {
        return_OVstatus_NOT_FOUND;
      } else {
        up_element* rev_elem = NULL;
        up_element* elem = up->elem;
        ov_word rev_last = 0;

        while (rev) {
          rev_elem = elem + (rev - 1);
          if (rev_elem->reverse_value == reverse_value) break;
          rev_last = rev;
          rev = rev_elem->reverse_next;
        }

        if (rev_elem) {
          ov_word rev_fwd_val = rev_elem->forward_value;
          ov_word fwd_hash = HASH(rev_fwd_val, mask);
          ov_word fwd = up->forward[fwd_hash];
          ov_word fwd_last = 0;
          up_element* fwd_elem = NULL;

          while (fwd) {
            fwd_elem = elem + (fwd - 1);
            if (fwd_elem == rev_elem) break;
            fwd_last = fwd;
            fwd = fwd_elem->forward_next;
          }

          if (rev && (rev == fwd)) {
            /* excise element */

            if (rev_last)
              up->elem[rev_last - 1].reverse_next = rev_elem->reverse_next;
            else
              up->reverse[rev_hash] = rev_elem->reverse_next;

            if (fwd_last)
              up->elem[fwd_last - 1].forward_next = fwd_elem->forward_next;
            else
              up->forward[fwd_hash] = fwd_elem->forward_next;

            /* store as inactive */

            rev_elem->active = OV_FALSE;
            rev_elem->forward_next = up->next_inactive;
            up->next_inactive = rev;
            up->n_inactive++;
            if (up->n_inactive >
                (up->size >> 1)) /* over half of bits are inactive */
              OVOneToOne_Pack(up);
            return_OVstatus_SUCCESS;
          }
        }
      }
    }
    return_OVstatus_NOT_FOUND;
  }
}

OVstatus OVOneToOne_Set(OVOneToOne* up, ov_word forward_value,
                        ov_word reverse_value) {
  if (!up) {
    return_OVstatus_NULL_PTR;
  } else {
    ov_word mask = up->mask;
    ov_word fwd_hash = HASH(forward_value, mask);
    ov_word rev_hash = HASH(reverse_value, mask);
    up_element* fwd_elem = NULL;
    up_element* rev_elem = NULL;
    ov_word fwd;
    ov_word rev;
    if (!mask) {
      fwd = 0;
      rev = 0;
    } else {
      fwd = up->forward[fwd_hash];
      rev = up->reverse[rev_hash];

      { /* find elements if they exist, and detect erroneous conditions */

        up_element* elem = up->elem;

        while (fwd) {
          fwd_elem = elem + (fwd - 1);
          if (fwd_elem->forward_value == forward_value) break;
          fwd = fwd_elem->forward_next;
        }
        while (rev) {
          rev_elem = elem + (rev - 1);
          if (rev_elem->reverse_value == reverse_value) break;
          rev = rev_elem->reverse_next;
        }
      }

      if ((fwd && (!rev)) || (rev && (!fwd))) {
        return_OVstatus_DUPLICATE;
      }
    }
    if (!(fwd || rev)) {
      ov_size new_index;
      /* new pair */
      if (up->n_inactive) {
        new_index = up->next_inactive;
        up->next_inactive = up->elem[new_index - 1].forward_next;
        up->n_inactive--;
      } else {
        if (up->elem && (!OVHeapArray_CHECK(up->elem, up_element, up->size))) {
          return_OVstatus_OUT_OF_MEMORY;
        } else {
          OVstatus result;
          if (OVreturn_IS_ERROR(result =
                                    Recondition(up, up->size + 1, OV_FALSE))) {
            return result;
          } else {
            /* guaranteed to succeed past this point, so we can increase size */
            new_index = ++up->size;
          }
        }
      }
      {
        up_element* elem = up->elem + (new_index - 1);
        elem->forward_value = forward_value;
        elem->reverse_value = reverse_value;
        elem->active = OV_TRUE;

        /* regenerate new hashes */
        mask = up->mask;
        fwd_hash = HASH(forward_value, mask);
        rev_hash = HASH(reverse_value, mask);

        {
          ov_word* forward_start_index = up->forward + fwd_hash;
          ov_word* reverse_start_index = up->reverse + rev_hash;

          elem->forward_next = *forward_start_index;
          *forward_start_index = new_index; /* note the +1 offset */
          elem->reverse_next = *reverse_start_index;
          *reverse_start_index = new_index; /* note the +1 offset */
        }
      }

    } else if (fwd_elem != rev_elem) {
      return_OVstatus_MISMATCH;
    } else {
      return_OVstatus_NO_EFFECT;
      /* exists and matched, so do nothing */
    }
  }
  return_OVstatus_SUCCESS;
}
