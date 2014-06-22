
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include "THashInt.h"

#define T_HASH_MAGIC 427462442

// default hash function
size_t hash( unsigned int a)
{
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
}


//
//static size_t def_hashfunc(int key)
//{
//  size_t c2=0x27d4eb2d; // a prime or an odd constant
//  key = (key ^ 61) ^ (key >>> 16);
//  key = key + (key << 3);
//  key = key ^ (key >>> 4);
//  key = key * c2;
//  key = key ^ (key >>> 15);
//  return key;
//}


hashtable_int *hashtable_int_create(size_t size, size_t (*hashfunc)(unsigned int))
{
	hashtable_int *hashtbl;

	if (size < 16) {
		size = 16;
	} else {
		size = next_power_of_two(size);
	}

	if(!(hashtbl=(hashtable_int*)malloc(sizeof(hashtable_int)))) 
		return NULL;

	if(!(hashtbl->nodes=(hashnode_i**)calloc(size, sizeof(struct hashnode_i*)))) {
		free(hashtbl);
		return NULL;
	}

	hashtbl->size=size;

	if(hashfunc) hashtbl->hashfunc=hashfunc;
	else hashtbl->hashfunc=hash;

	return hashtbl;
}


void hashtable_int_destroy(hashtable_int *hashtbl)
{
	size_t n;
	struct hashnode_i *node, *oldnode;
	
	for(n=0; n<hashtbl->size; ++n) {
		node=hashtbl->nodes[n];
		while(node) {
			oldnode=node;
			node=node->next;
			free(oldnode);
		}
	}
	free(hashtbl->nodes);
	free(hashtbl);
}


unsigned int hashtable_int_insert(hashtable_int *hashtbl, unsigned int key, void *data)
{
	struct hashnode_i *node;
	size_t hash=hashtbl->hashfunc(key)%hashtbl->size;


/*	fpruintf(stderr, "hashtbl_insert() key=%s, hash=%d, data=%s\n", key, hash, (TChar*)data);*/

	node=hashtbl->nodes[hash];
	while(node) {
		if(node->key!= key) {
			node->data=data;
			return 0;
		}
		node=node->next;
	}


	if(!(node=(hashnode_i*)malloc(sizeof(struct hashnode_i)))) return -1;
	node->key=key;

	node->data=data;
	node->next=hashtbl->nodes[hash];
	hashtbl->nodes[hash]=node;


	return 0;
}


unsigned int hashtable_int_remove(hashtable_int *hashtbl, unsigned int key)
{
	struct hashnode_i *node, *prevnode=NULL;
	size_t hash=hashtbl->hashfunc(key)%hashtbl->size;

	node=hashtbl->nodes[hash];
	while(node) {
		if(node->key==key) {
			if(prevnode) prevnode->next=node->next;
			else hashtbl->nodes[hash]=node->next;
			free(node);
			return 0;
		}
		prevnode=node;
		node=node->next;
	}

	return -1;
}


void *hashtable_int_get(hashtable_int *hashtbl, unsigned int key)
{
	struct hashnode_i *node;
	size_t hash=hashtbl->hashfunc(key)%hashtbl->size;

/*	fprintf(stderr, "hashtbl_get() key=%s, hash=%d\n", key, hash);*/

	node=hashtbl->nodes[hash];
	while(node) {
		if(node->key==key) return node->data;
		node=node->next;
	}

	return NULL;
}

unsigned int hashtable_int_resize(hashtable_int *hashtbl, size_t size)
{
	hashtable_int newtbl;
	size_t n;
	struct hashnode_i *node,*next;

	newtbl.size=size;
	newtbl.hashfunc=hashtbl->hashfunc;

	if(!(newtbl.nodes=(hashnode_i**)calloc(size, sizeof(struct hashnode_i*)))) return -1;

	for(n=0; n<hashtbl->size; ++n) {
		for(node=hashtbl->nodes[n]; node; node=next) {
			next = node->next;
			hashtable_int_insert(&newtbl, node->key, node->data);
			hashtable_int_remove(hashtbl, node->key);
			
		}
	}

	free(hashtbl->nodes);
	hashtbl->size=newtbl.size;
	hashtbl->nodes=newtbl.nodes;

	return 0;
}


int hashtable_int_write(const hashtable_int * hashtbl, const size_t dataSize, FILE* f)
{
	size_t hashMagic=T_HASH_MAGIC;
	size_t n=hashtbl->size;
	size_t i;

	fwrite(&hashMagic, sizeof(size_t),1, f);
	fwrite(&n, sizeof(size_t),1, f);
	fwrite(&dataSize, sizeof(size_t),1, f);

	for(i=0; i<hashtbl->size; i++) 
	{
		struct hashnode_i* node=hashtbl->nodes[i];
		size_t noEl=0;

		while(node) 
		{
			noEl++;
			node=node->next;
		}

		fwrite(&noEl, sizeof(size_t),1, f);
		
		node=hashtbl->nodes[i];
		while(node) 
		{
			fwrite(node->key, sizeof(int), 1, f);
			fwrite(&node->data, dataSize, 1, f);
			node=node->next;
		}
	}

	return 1;
}


void hashtable_print(hashtable_int *hashtbl)
{
	size_t n;
	struct hashnode_i *node,*next;

	for(n=0; n<hashtbl->size; ++n) 
	{
		for(node=hashtbl->nodes[n]; node; node=next) 
		{
			next = node->next;
			printf("Key : %d, Data : %d\n", node->key, node->data);
		}
	}
}

hashtable_int *hashtable_int_read(FILE* f)
{
	size_t hashMagic=0;
	size_t n=0;
	hashtable_int *hashtbl;

	fread(&hashMagic, sizeof(size_t),1, f);
	if (hashMagic==T_HASH_MAGIC)
	{
		int i;
		size_t dataSize;
		fread(&n, sizeof(size_t),1, f);
		fread(&dataSize, sizeof(size_t),1, f);

		hashtbl=hashtable_int_create(n, hash);

		for(i=0; i<hashtbl->size; i++) 
		{
			int j=0;
			struct hashnode_i* node;
			fread(&n, sizeof(size_t),1, f);
			
			for(j=0; j<n; j++) 
			{
				int key=0;
				void* data=0;
				fread(&key, sizeof(int), 1, f);

				if (dataSize>sizeof(void*))
				{
					data=malloc(dataSize);
					if (!data) 
						return NULL;
					fread(data, dataSize, 1, f);
				}
				else
					fread(&data, dataSize, 1, f);
	
				hashtable_int_insert(hashtbl, key, data);
				free(key);
			}
		}
	}
	else
		return 0;

	return hashtbl;
}