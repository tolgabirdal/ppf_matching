//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

// Note - The x86 and x64 versions do _not_ produce the same results, as the
// algorithms are optimized for their respective platforms. You can still
// compile and run any of them on any platform, but your performance with the
// non-native version will be less than optimal.

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

// Microsoft Visual Studio

#if defined(_MSC_VER)

#define FORCE_INLINE    __forceinline

#include <stdlib.h>
#include "TBaseInclude.h"

#define ROTL32(x,y)     _rotl(x,y)
#define ROTL64(x,y)     _rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else   // defined(_MSC_VER)

#define FORCE_INLINE __attribute__((always_inline))

inline int rotl32 ( int x, int8_t r )
{
  return (x << r) | (x >> (32 - r));
}

inline TLong rotl64 ( TLong x, int8_t r )
{
  return (x << r) | (x >> (64 - r));
}

#define ROTL32(x,y)     rotl32(x,y)
#define ROTL64(x,y)     rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

FORCE_INLINE int getblock ( const int * p, int i )
{
  return p[i];
}

FORCE_INLINE TLong getblock64 ( const TLong * p, int i )
{
  return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

FORCE_INLINE int fmix ( int h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

//----------

FORCE_INLINE TLong fmix64 ( TLong k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}

//-----------------------------------------------------------------------------

FORCE_INLINE void MurmurHash3_x86_32 ( const void * key, int len,
                          int seed, void * out )
{
  const unsigned char * data = (const unsigned char*)key;
  const int nblocks = len / 4;

  int h1 = seed;

  int c1 = 0xcc9e2d51;
  int c2 = 0x1b873593;
  int i;

  //----------
  // body

  const int * blocks = (const int *)(data + nblocks*4);
	unsigned char * tail;
  int k1 = 0;

  for(i = -nblocks; i; i++)
  {
    //int k1 = getblock(blocks,i);
	  int k1 = blocks[i];

    k1 *= c1;
    k1 = ROTL32(k1,15);
    k1 *= c2;
    
    h1 ^= k1;
    h1 = ROTL32(h1,13); 
    h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  tail = (unsigned char*)(data + nblocks*4);

  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix(h1);

  *(int*)out = h1;
} 
