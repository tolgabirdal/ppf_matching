
#ifndef THASHINT_H
#define THASHINT_H

#include<stdlib.h>

typedef struct hashnode_i {
	unsigned int key;
	void *data;
	struct hashnode_i *next;
} hashnode_i ;

typedef struct HSHTBL_i {
	size_t size;
	struct hashnode_i **nodes;
	size_t (*hashfunc)(unsigned int);
} hashtable_int;


static __inline unsigned int next_power_of_two(unsigned int value)
{
	/* Round up to the next highest power of 2 */
	/* from http://www-graphics.stanford.edu/~seander/bithacks.html */

	--value;
	value |= value >> 1;
	value |= value >> 2;
	value |= value >> 4;
	value |= value >> 8;
	value |= value >> 16;
	++value;

	return value;
}

#if defined (__cplusplus)
extern "C" {
#endif

	hashtable_int *hashtable_int_create(size_t size, size_t (*hashfunc)(unsigned int));
	void hashtable_int_destroy(hashtable_int *hashtbl);
	unsigned int hashtable_int_insert(hashtable_int *hashtbl, unsigned int key, void *data);
	unsigned int hashtable_int_insert_hashed(hashtable_int *hashtbl, unsigned int key, void *data);
	unsigned int hashtable_int_remove(hashtable_int *hashtbl, unsigned int key);
	void *hashtable_int_get(hashtable_int *hashtbl, unsigned int key);
	hashnode_i* hashtable_int_get_bucket_hashed(hashtable_int *hashtbl, unsigned int key);
	unsigned int hashtable_int_resize(hashtable_int *hashtbl, size_t size);
	hashtable_int *hashtable_int_clone(hashtable_int *hashtbl);
	hashtable_int *hashtable_int_read(FILE* f);
	int hashtable_int_write(const hashtable_int * hashtbl, const size_t dataSize, FILE* f);
	void hashtable_int_print(hashtable_int *hashtbl);

#if defined (__cplusplus)
}
#endif

#endif

