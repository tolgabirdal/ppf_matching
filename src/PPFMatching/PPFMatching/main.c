
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include <tommy.h>

typedef struct THash {
    int id;
    void* data;
	tommy_node node;
} THash;

typedef struct
{
	int magic;
	double maxDist, angleStep, distStep;
	//Tobject inputPC;
	double* alpha_m;
	int n;
	THash* hashTable;
	//TPoseList* Poses;
}TPPFModelPC;


#define T_PPF_LENGTH 5

void compute_ppf_features(const double p1[4], const double n1[4],
						  const double p2[4], const double n2[4],
						  double f[4])
{
	double delta[4] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2], 0};
	double f1,f2,f3,f4;

	f[3] = sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);

	if (f[3])
	{
		delta[0] /= f[3];
		delta[1] /= f[3];
		delta[2] /= f[3];
	}

	f[0] = n1[0] * delta[0] + n1[1] * delta[1] + n1[2] * delta[2];
	f[1] = n2[0] * delta[0] + n2[1] * delta[1] + n2[2] * delta[2];
	f[2] = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
}

void 

 int compare(const void* arg, const void* obj)
 {
     return *(const int*)arg != ((THash*)obj)->id;
 }

// test hash table
int main()
{
	int value_to_find = 227;
	THash* objFound ;
	THash* hashNode = (THash*)malloc(sizeof(THash));

	tommy_hashtable* hashTable = (tommy_hashtable*)malloc(sizeof(tommy_hashtable));

	tommy_hashtable_init(hashTable, 100000);

	hashNode->id = 227;
	hashNode->data = 5;
	tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

	hashNode = (THash*)malloc(sizeof(THash));
	hashNode->id = 112;
	hashNode->data = 2;
	tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

	objFound = tommy_hashtable_search(hashTable, compare, &value_to_find, tommy_inthash_u32(value_to_find));

	if (objFound)
		printf("%d\n", objFound->id);
	else
		printf("Nof found\n");
	

	return 0;
}