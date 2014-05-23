
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*#include "opencv2/features2d.hpp"
#include "opencv2/nonfree.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc.hpp"*/
#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/rgbd.hpp"
#include "helpers.h"
#include "c_utils.h"
#include "hash_murmur.h"
#include <tommy.h>

#include <Eigen/Core>
#include "opencv2/core/eigen.hpp"


using namespace cv;

typedef struct THash {
	int id;
	void* data;
	tommy_node node;
} THash;

typedef struct
{
	int magic;
	double maxDist, angleStep, distStep;
	Mat inputPC;
	flann::Index pcTree;
	Mat alpha_m;
	int n;
	THash* hashTable;
}TPPFModelPC;


#define T_PPF_LENGTH 5

// compute per point PPF as in paper
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

// quantize ppf and hash it for proper indexing
int hash_ppf(const double f[4], const double AngleStep, const double DistanceStep)
{
	const int d1 = (int) (floor ((double)f[0] / (double)AngleStep));
	const int d2 = (int) (floor ((double)f[1] / (double)AngleStep));
	const int d3 = (int) (floor ((double)f[2] / (double)AngleStep));
	const int d4 = (int) (floor ((double)f[3] / (double)DistanceStep));
	int key[4]={d1,d2,d3,d4};
	int hashKey=0;
	MurmurHash3_x86_32(key, 4, 42, &hashKey);
	return hashKey;
}


int t_create_match_model_pc(const Mat PC, const double RelSamplingStep, const double RelativeAngleStep, const double RelativeDistanceStep, TPPFModelPC** Model3D)
{


	return 0;
}

int compare(const void* arg, const void* obj)
{
	return *(const int*)arg != ((THash*)obj)->id;
}

// test bounding box
int main()
{
	int useNormals = 1;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);

	TPPFModelPC* Model3D = 0;

	float xRange[2], yRange[2], zRange[2];
	compute_obb(pc, xRange, yRange, zRange);
	

//	t_create_match_model_pc(pc, const double AngleStep, const double DistanceStep, &Model3D);

}

// test point cloud reading & visualization
// for now uniform sampling looks better
int main_ply()
{
	int useNormals = 1;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";

	Mat pc = load_ply_simple(fn, numVert, useNormals);
	
	Mat spc1 = sample_pc_uniform(pc, numVert/350);
	Mat spc2 = sample_pc_random(pc, 350);
	visualize_pc(pc, useNormals, "PC1");
	visualize_pc(spc1, useNormals, "Uniform Sampled");
	visualize_pc(spc2, useNormals, "Random Sampled");

	return 0;
}
    
// test hash table
int main_hash_test()
{
	int value_to_find = 227;
	THash* objFound ;
	THash* hashNode = (THash*)malloc(sizeof(THash));

	tommy_hashtable* hashTable = (tommy_hashtable*)malloc(sizeof(tommy_hashtable));

	// 262144 = 2^18
	int size = next_power_of_two(100000);
	tommy_hashtable_init(hashTable, size);

	hashNode->id = 227;
	hashNode->data = (void*)5;
	tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

	hashNode = (THash*)malloc(sizeof(THash));
	hashNode->id = 112;
	hashNode->data = (void*)2;
	tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

	objFound = (THash*)tommy_hashtable_search(hashTable, compare, &value_to_find, tommy_inthash_u32(value_to_find));

	if (objFound)
		printf("%d\n", objFound->id);
	else
		printf("Nof found\n");


	return 0;
}