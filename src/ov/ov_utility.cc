#include "ov/ov_port.h"
#include "ov/ov_utility.h"

void ov_utility_zero_range(void* start, void* stop) {
  char* p = (char*)start;
  char* q = (char*)stop;
  if (q - p) memset(p, 0, q - p);
}
