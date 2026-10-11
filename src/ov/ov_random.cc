
#include "ov/ov_random.h"

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

/* Modified by Warren DeLano for use with the ObVec library
 * (1) added object orientation (eliminated globals & statics,
 *                            & included instance pointer)
 * (2) renamed methods
 */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

struct _OVRandom {
  OVHeap* heap;
  ov_uint32 mt[N]; /* the array for the state vector  */
  ov_int32 mti;    /* mti==N+1 means mt[N] is not initialized */
  ov_uint32 mag01[2];
};

void OVRandom_Del(OVRandom* I) { OVHeap_FREE_AUTO_NULL(I->heap, I); }

/* initializes mt[N] with a seed */

/*void init_genrand(unsigned long s)*/
OVRandom* OVRandom_NewBySeed(OVHeap* heap, ov_uint32 seed) {
  OVRandom* I = OVHeap_ALLOC(heap, OVRandom);

  if (I) {
    ov_uint32* mt = I->mt;
    ov_int32 mti = N + 1;

    I->heap = heap;
    mt[0] = seed & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
    I->mag01[0] = 0x0UL;
    I->mag01[1] = MATRIX_A;
    I->mti = mti;
  }
  return (I);
}

/* initialize by an array with array-length */

/* init_key is the array for initializing keys */

/* key_length is its length */

/* slight change for C++, 2004/2/26 */
OVRandom* OVRandom_NewByArray(OVHeap* heap, ov_uint32 init_key[],
                              ov_int32 key_length) {
  OVRandom* I = OVRandom_NewBySeed(heap, 19650218UL);
  if (I) {
    int i, j, k;
    ov_uint32* mt = I->mt;
    i = 1;
    j = 0;
    k = (N > key_length ? N : key_length);
    for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) +
              init_key[j] + j; /* non linear */
      mt[i] &= 0xffffffffUL;   /* for WORDSIZE > 32 machines */
      i++;
      j++;
      if (i >= N) {
        mt[0] = mt[N - 1];
        i = 1;
      }
      if (j >= key_length) j = 0;
    }
    for (k = N - 1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) -
              i;             /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i >= N) {
        mt[0] = mt[N - 1];
        i = 1;
      }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
  }
  return I;
}

/* generates a random number on [0,0xffffffff]-interval */
ov_uint32 OVRandom_Get_int32(OVRandom* I) {
  ov_uint32 y;
  ov_uint32* mt = I->mt;
  ov_uint32* mag01 = I->mag01;
  ov_uint32 mti = I->mti;

  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mti >= N) { /* generate N words at one time */
    int kk;

    for (kk = 0; kk < N - M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
  }

  y = mt[mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  I->mti = mti;
  return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
ov_int32 OVRandom_Get_int31(OVRandom* I) {
  return (ov_int32)(OVRandom_Get_int32(I) >> 1);
}

/* generates a random number on [0,1]-real-interval */
ov_float64 OVRandom_Get_float64_inc1(OVRandom* I) {
  return OVRandom_Get_int32(I) * (1.0 / 4294967295.0);
  /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
ov_float64 OVRandom_Get_float64_exc1(OVRandom* I) {
  return OVRandom_Get_int32(I) * (1.0 / 4294967296.0);
  /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
ov_float64 OVRandom_Get_float64_exc01(OVRandom* I) {
  return (((double)OVRandom_Get_int32(I)) + 0.5) * (1.0 / 4294967296.0);
  /* divided by 2^32 */
}
