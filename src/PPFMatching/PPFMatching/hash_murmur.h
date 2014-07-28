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

#define ROTL32(x,y)     _rotl(x,y)
#define ROTL64(x,y)     _rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else   // defined(_MSC_VER)

#define long long TLong

#define FORCE_INLINE __attribute__((always_inline))

inline int rotl32 ( int x, int8_t r )
{
	return (x << r) | (x >> (32 - r));
}

inline long long rotl64 ( long long x, int8_t r )
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

FORCE_INLINE long long getblock64 ( const long long * p, int i )
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


//-----------------------------------------------------------------------------

FORCE_INLINE void MurmurHash3_x86_32 ( const void * key, int len, int seed, void * out )
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
/*
inline void bmix32 ( unsigned int & h1, unsigned int & k1, unsigned int & c1, unsigned int & c2 )
{
k1 *= c1; 
k1  = _rotl(k1,11); 
k1 *= c2;
h1 ^= k1;

h1 = h1*3+0x52dce729;

c1 = c1*5+0x7b7d159c;
c2 = c2*5+0x6bce6396;
}*/

//
//FORCE_INLINE long fmix64 ( long k )
//{
//	k ^= k >> 33;
//	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
//	k ^= k >> 33;
//	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
//	k ^= k >> 33;
//
//	return k;
//}
//
//inline void bmix32 ( unsigned int & h1, unsigned int & h2, unsigned int & k1, unsigned int & k2, unsigned int & c1, unsigned int & c2 )
//{
//	k1 *= c1; 
//	k1  = _rotl(k1,11); 
//	k1 *= c2;
//	h1 ^= k1;
//	h1 += h2;
//
//	h2 = _rotl(h2,17);
//
//	k2 *= c2; 
//	k2  = _rotl(k2,11);
//	k2 *= c1;
//	h2 ^= k2;
//	h2 += h1;
//
//	h1 = h1*3+0x52dce729;
//	h2 = h2*3+0x38495ab5;
//
//	c1 = c1*5+0x7b7d159c;
//	c2 = c2*5+0x6bce6396;
//}
//
//inline unsigned int fmix32 ( unsigned int h )
//{
//	h ^= h >> 16;
//	h *= 0x85ebca6b;
//	h ^= h >> 13;
//	h *= 0xc2b2ae35;
//	h ^= h >> 16;
//
//	return h;
//}
//
//FORCE_INLINE void MurmurHash3_x86_64 ( const void * key, const int len, const unsigned int seed, void * out )
//{
//	const unsigned char * data = (const unsigned char*)key;
//	const int nblocks = len / 8;
//
//	unsigned int h1 = 0x8de1c3ac ^ seed;
//	unsigned int h2 = 0xbab98226 ^ seed;
//
//	unsigned int c1 = 0x95543787;
//	unsigned int c2 = 0x2ad7eb25;
//
//	//----------
//	// body
//
//	const unsigned int * blocks = (const unsigned int *)(data + nblocks*8);
//
//	for(int i = -nblocks; i; i++)
//	{
//		unsigned int k1 = blocks[i*2];
//		unsigned int k2 = blocks[i*2+1];
//
//		bmix32(h1,h2,k1,k2,c1,c2);
//	}
//
//	//----------
//	// tail
//
//	const unsigned char * tail = (const unsigned char*)(data + nblocks*8);
//
//	unsigned int k1 = 0;
//	unsigned int k2 = 0;
//
//	switch(len & 7)
//	{
//	case 7: k2 ^= tail[6] << 16;
//	case 6: k2 ^= tail[5] << 8;
//	case 5: k2 ^= tail[4] << 0;
//	case 4: k1 ^= tail[3] << 24;
//	case 3: k1 ^= tail[2] << 16;
//	case 2: k1 ^= tail[1] << 8;
//	case 1: k1 ^= tail[0] << 0;
//		bmix32(h1,h2,k1,k2,c1,c2);
//	};
//
//	//----------
//	// finalization
//
//	h2 ^= len;
//
//	h1 += h2;
//	h2 += h1;
//
//	h1 = fmix32(h1);
//	h2 = fmix32(h2);
//
//	h1 += h2;
//	h2 += h1;
//
//	((unsigned int*)out)[0] = h1;
//	((unsigned int*)out)[1] = h2;
//}
//
//
//
////----------
//// Block mix - combine the key bits with the hash bits and scramble everything
//
//inline void bmix64 ( unsigned long & h1, unsigned long & h2, unsigned long & k1, unsigned long & k2, unsigned long & c1, unsigned long & c2 )
//{
//	k1 *= c1; 
//	k1  = _rotl64(k1,23); 
//	k1 *= c2;
//	h1 ^= k1;
//	h1 += h2;
//
//	h2 = _rotl64(h2,41);
//
//	k2 *= c2; 
//	k2  = _rotl64(k2,23);
//	k2 *= c1;
//	h2 ^= k2;
//	h2 += h1;
//
//	h1 = h1*3+0x52dce729;
//	h2 = h2*3+0x38495ab5;
//
//	c1 = c1*5+0x7b7d159c;
//	c2 = c2*5+0x6bce6396;
//}
//
////----------
//// Finalization mix - avalanches all bits to within 0.05% bias
//
//inline unsigned long fmix64 ( unsigned long k )
//{
//	k ^= k >> 33;
//	k *= 0xff51afd7ed558ccd;
//	k ^= k >> 33;
//	k *= 0xc4ceb9fe1a85ec53;
//	k ^= k >> 33;
//
//	return k;
//}
//
////----------
//
//void MurmurHash3_x64_128 ( const void * key, const int len, const unsigned int seed, void * out )
//{
//	const unsigned char * data = (const unsigned char*)key;
//	const int nblocks = len / 16;
//
//	unsigned long h1 = 0x9368e53c2f6af274 ^ seed;
//	unsigned long h2 = 0x586dcd208f7cd3fd ^ seed;
//
//	unsigned long c1 = 0x87c37b91114253d5;
//	unsigned long c2 = 0x4cf5ad432745937f;
//
//	//----------
//	// body
//
//	const unsigned long * blocks = (const unsigned long *)(data);
//
//	for(int i = 0; i < nblocks; i++)
//	{
//		unsigned long k1 = blocks[i*2];
//		unsigned long k2 = blocks[i*2+1];
//
//		bmix64(h1,h2,k1,k2,c1,c2);
//	}
//
//	//----------
//	// tail
//
//	const unsigned char * tail = (const unsigned char*)(data + nblocks*16);
//
//	unsigned long k1 = 0;
//	unsigned long k2 = 0;
//
//	switch(len & 15)
//	{
//	case 15: k2 ^= (unsigned long)(tail[14]) << 48;
//	case 14: k2 ^= (unsigned long)(tail[13]) << 40;
//	case 13: k2 ^= (unsigned long)(tail[12]) << 32;
//	case 12: k2 ^= (unsigned long)(tail[11]) << 24;
//	case 11: k2 ^= (unsigned long)(tail[10]) << 16;
//	case 10: k2 ^= (unsigned long)(tail[ 9]) << 8;
//	case  9: k2 ^= (unsigned long)(tail[ 8]) << 0;
//
//	case  8: k1 ^= (unsigned long)(tail[ 7]) << 56;
//	case  7: k1 ^= (unsigned long)(tail[ 6]) << 48;
//	case  6: k1 ^= (unsigned long)(tail[ 5]) << 40;
//	case  5: k1 ^= (unsigned long)(tail[ 4]) << 32;
//	case  4: k1 ^= (unsigned long)(tail[ 3]) << 24;
//	case  3: k1 ^= (unsigned long)(tail[ 2]) << 16;
//	case  2: k1 ^= (unsigned long)(tail[ 1]) << 8;
//	case  1: k1 ^= (unsigned long)(tail[ 0]) << 0;
//		bmix64(h1,h2,k1,k2,c1,c2);
//	};
//
//	//----------
//	// finalization
//
//	h2 ^= len;
//
//	h1 += h2;
//	h2 += h1;
//
//	h1 = fmix64(h1);
//	h2 = fmix64(h2);
//
//	h1 += h2;
//	h2 += h1;
//
//	((unsigned long*)out)[0] = h1;
//	((unsigned long*)out)[1] = h2;
//}
//
////-----------------------------------------------------------------------------
//// If we need a smaller hash value, it's faster to just use a portion of the 
//// 128-bit hash
//
//void MurmurHash3_x64_32 ( const void * key, int len, unsigned int seed, void * out )
//{
//	unsigned int temp[4];
//
//	MurmurHash3_x64_128(key,len,seed,temp);
//
//	*(unsigned int*)out = temp[0];
//}
//
////----------
//
//void MurmurHash3_x64_64 ( const void * key, int len, unsigned int seed, void * out )
//{
//	unsigned long temp[2];
//
//	MurmurHash3_x64_128(key,len,seed,temp);
//
//	*(unsigned long*)out = temp[0];
//} 
