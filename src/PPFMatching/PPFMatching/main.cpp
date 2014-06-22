
//<<<<<<< HEAD
/*#include "opencv2/features2d.hpp"
#include "opencv2/nonfree.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc.hpp"*/
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/rgbd.hpp"
#include "opencv2/flann/flann.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//=======
//#include <opencv2/opencv.hpp>
//>>>>>>> 37dd8f1b41c87a7848dc42b1965e88fa5a6115f9
#include "helpers.h"
#include "visualize_win.h"
#include "c_utils.h"
#include "THashInt.h"
#include "hash_murmur.h"
#include <tommy.h>

using namespace cv;

typedef struct THash {
	int id;
	int i, j, ppfInd;
	tommy_node node;
} THash;

class PPFPose
{
public:
	PPFPose(){alpha=0; modelIndex=0; numVotes=0;};
	PPFPose(double Alpha =0, unsigned int ModelIndex=0, unsigned int NumVotes=0)
	{
		alpha = Alpha;
		modelIndex = ModelIndex;
		numVotes = NumVotes;
	};
	~PPFPose(){};

	double alpha;
	unsigned int modelIndex;
	unsigned int numVotes;
};

class TPPFModelPC
{
public:

	TPPFModelPC() {} ;
	~TPPFModelPC() {};

	int magic;
	double maxDist, angleStep, distStep;
	Mat inputPC, PPF;
	cvflann::Index<Distance_32F>* flannIndex;
	//Mat alpha_m;
	int n, numRefPoints, sampledStep, ppfStep;
	tommy_hashtable* hashTable;
//	hashtable_int* hashTable;
};


#define T_PPF_LENGTH 5

// compute per point PPF as in paper
void compute_ppf_features(const double p1[4], const double n1[4],
						  const double p2[4], const double n2[4],
						  double f[4])
{
	/* 
		Vectors will be defined as of length 4 instead of 3, because of:
		- Further SIMD vectorization
		- Cache alignment
	*/

	double d[4] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2], 0};
	double c[4];

	double norm = TNorm3(d); 
	f[3] = norm;

	if (norm)
	{
		d[0] /= f[3];
		d[1] /= f[3];
		d[2] /= f[3];
	}
	else
	{
		// TODO: Handle this
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		return ;
	}

	/*
		Issues of numerical stability is of concern here.
		Bertram's suggestion: atan2(a dot b, |axb|)
		My correction : 
		I guess it should be: angle = atan2(norm(cross(a,b)), dot(a,b))
		The macro is implemented accordingly.
	*/

	TAngle3(n1, d, c, f[0]);
	TAngle3(n2, d, c, f[1]);
	TAngle3(n1, n1, c, f[2]);
}

// simple hashing
int hash_ppf_simple(const double f[4], const double AngleStep, const double DistanceStep)
{
	const int d1 = (int) (floor ((double)f[0] / (double)AngleStep));
	const int d2 = (int) (floor ((double)f[1] / (double)AngleStep));
	const int d3 = (int) (floor ((double)f[2] / (double)AngleStep));
	const int d4 = (int) (floor ((double)f[3] / (double)DistanceStep));
	
	return (d1 | (d2<<8) | (d3<<16) | (d4<<24));
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

size_t hash_murmur(unsigned int key)
{
	size_t hashKey=0;
	MurmurHash3_x86_32((void*)&key, 4, 42, &hashKey);
	return hashKey;
}

// tommy's compare function
int compare(const void* arg, const void* obj)
{
	return *(const int*)arg != ((THash*)obj)->id;
}

// TODO: An initial attempt. I will double check this
double compute_alpha(const double p1[4], const double n1[4], const double p2[4])
{
	double Tmg[3], mpt[3], row2[3], row3[3], alpha;

	compute_transform_rt_yz(p1, n1, row2, row3, Tmg);

	mpt[1] = Tmg[1] + row2[0] * p2[0] + row2[1] * p2[1] + row2[2] * p2[2];
    mpt[2] = Tmg[2] + row3[0] * p2[0] + row3[1] * p2[1] + row3[2] * p2[2];

	alpha=atan2(-mpt[2], mpt[1]);

	if ( alpha != alpha)
    {
		printf("NaN value!\n");
		return 0;
	}

	if (sin(alpha)*mpt[2]<0.0)
		alpha=-alpha;

	return (-alpha);
}


Mat compute_ppf_pc_train(const Mat PC, const double distanceStep, const double angleStep)
{
	Mat PPFMat = Mat(PC.rows*PC.rows, T_PPF_LENGTH, CV_32FC1);

	for (int i=0; i<PC.rows; i++)
	{
		for (int j=0; j<PC.rows; j++)
		{
			// cannnot compute the ppf with myself
			if (i!=j)
			{
				float* f1 = (float*)(&PC.data[i * PC.step]);
				float* f2 = (float*)(&PC.data[j * PC.step]);
				const double p1[4] = {f1[0], f1[1], f1[2], 1};
				const double p2[4] = {f2[0], f1[1], f1[2], 1};
				const double n1[4] = {f1[3], f1[4], f1[5], 1};
				const double n2[4] = {f2[3], f1[4], f1[5], 1};

				double f[4]={0};
				compute_ppf_features(p1, n1, p2, n2, f);
				double alpha = compute_alpha(p1, n1, p2);

				int corrInd = i*PC.rows+j;
				PPFMat.data[ corrInd ] = f[0];
				PPFMat.data[ corrInd + 1 ] = f[1];
				PPFMat.data[ corrInd + 2 ] = f[2];
				PPFMat.data[ corrInd + 3 ] = f[3];
				PPFMat.data[ corrInd + 4 ] = (float)alpha;
			}
		}
	}

	return PPFMat;
}

// TODO: Check all step sizes to be positive
Mat train_pc_ppf(const Mat PC, const double sampling_step_relative, const double distance_step_relative, const double angle_step_relative, TPPFModelPC** Model3D)
{
	const int numPoints = PC.rows;

	// compute bbox
	float xRange[2], yRange[2], zRange[2];
	compute_obb(PC, xRange, yRange, zRange);

	// compute sampling step from diameter of bbox
	float dx = xRange[1] - xRange[0];
	float dy = yRange[1] - yRange[0];
	float dz = zRange[1] - zRange[0];
	float diameter = sqrt ( dx * dx + dy * dy + dz * dz );
	float distanceStep = diameter * sampling_step_relative;

	Mat sampled = sample_pc_octree(PC, xRange, yRange, zRange, sampling_step_relative);

    float angleStepRadians = (360/angle_step_relative)*M_PI/180;

	tommy_hashtable* hashTable = (tommy_hashtable*)malloc(sizeof(tommy_hashtable));
	// 262144 = 2^18
	int size = next_power_of_two(sampled.rows*sampled.rows);
	tommy_hashtable_init(hashTable, size);

	//hashtable_int* hashTable = hashtable_int_create(sampled.rows*sampled.rows, NULL);
	
	*Model3D = new TPPFModelPC();

	int numPPF = sampled.rows*sampled.rows;
	(*Model3D)->PPF = Mat(numPPF, T_PPF_LENGTH, CV_32FC1);
	int ppfStep = (*Model3D)->PPF.step;
	int sampledStep = sampled.step;
	
	// TODO: Maybe I could sample 1/5th of them here. Check the performance later.
	int numRefPoints = sampled.rows;
	for (int i=0; i<numRefPoints; i++)
	{
		for (int j=0; j<sampled.rows; j++)
		{
			// cannnot compute the ppf with myself
			if (i!=j)
			{
				float* f1 = (float*)(&sampled.data[i * sampledStep]);
				float* f2 = (float*)(&sampled.data[j * sampledStep]);
				const double p1[4] = {f1[0], f1[1], f1[2], 0};
				const double p2[4] = {f2[0], f1[1], f1[2], 0};
				const double n1[4] = {f1[3], f1[4], f1[5], 0};
				const double n2[4] = {f2[3], f1[4], f1[5], 0};

				double f[4]={0};
				compute_ppf_features(p1, n1, p2, n2, f);
				unsigned int hash = hash_ppf_simple(f, angleStepRadians, distanceStep);
				double alpha = compute_alpha(p1, n1, p2);
				unsigned int corrInd = i*numRefPoints+j;
				unsigned int ppfInd = corrInd*ppfStep;

				//hashtable_int_insert(hashTable, hash, (void*)corrInd);
				//printf("%f %f %f %f \n", f[0], f[1], f[2], f[3]);
				//printf("%d\n", hash);

				THash* hashNode = (THash*)calloc(1, sizeof(THash));
				hashNode->id = hash;
				//hashNode->data = (void*)corrInd;
				hashNode->i = i;
				hashNode->j = j;
				hashNode->ppfInd = ppfInd;
				tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

				float* ppfRow = (float*)(&(*Model3D)->PPF.data[ ppfInd ]);
				ppfRow[0] = f[0];
				ppfRow[1] = f[1];
				ppfRow[2] = f[2];
				ppfRow[3] = f[3];
				ppfRow[4] = (float)alpha;
			}
		}
	}

	//*Model3D = (TPPFModelPC*)calloc(1, sizeof(TPPFModelPC));

	(*Model3D)->angleStep = angleStepRadians;
	(*Model3D)->distStep = distanceStep;
	(*Model3D)->hashTable = hashTable;
	(*Model3D)->sampledStep = sampledStep;
	(*Model3D)->ppfStep = ppfStep;
	(*Model3D)->numRefPoints = numRefPoints;
	//(*Model3D)->PPF = PPFMat;

	//(*Model3D)->magic=T_MAGIC_VAL_PC_MODEL;



	//return compute_ppf_pc_train(pc, distanceStep, angleStepRadians);
	
	return Mat();
}

Mat t_load_ppf_model(const char* FileName)
{
	Mat ppf = Mat();

	return ppf;
}

void t_match_pc_ppf(Mat pc, float SearchRadius, int SampleStep, TPPFModelPC* ppfModel)
{
	cvflann::Matrix<float> data;
	int i;
	int numNeighbors = 1000;
	int numAngles = (int) (floor (2 * M_PI / ppfModel->angleStep));
	int max_votes_i = 0, max_votes_j = 0;
	int max_votes = 0;
	unsigned int* accumulator;
	cvflann::Index<Distance_32F>* flannIndex;
	float angleStepRadians = ppfModel->angleStep;
	float distanceStep = ppfModel->distStep;
	int sampledStep = ppfModel->sampledStep;
	int ppfStep = ppfModel->ppfStep;
	int numRefPoints = ppfModel->numRefPoints;
	unsigned int n = numRefPoints;

	// allocate the accumulator
	accumulator = (unsigned int*)calloc(numAngles*n, sizeof(unsigned int));
	
	
	// obtain the tree representation for fast search
	flannIndex  = (cvflann::Index<Distance_32F>*)index_pc_flann(pc, data);

	cv::Mat1i ind(pc.rows, numNeighbors);
	cvflann::Matrix<int> indices((int*) ind.data, ind.rows, ind.cols);
	cvflann::Matrix<float> dists(new float[pc.rows*numNeighbors], pc.rows, numNeighbors);

	// TODO: Can be parallelized!
#pragma omp parallel for
	for (i = 0; i < pc.rows; i += 15)
	{
		int j;

		float* f1 = (float*)(&pc.data[i * pc.step]);
		const double p1[4] = {f1[0], f1[1], f1[2], 0};
		const double n1[4] = {f1[3], f1[4], f1[5], 0};
		double p1t[4];
		double row1[3]={0}, row2[3]={0}, row3[3]={0}, Tsg[3]={0};

		compute_transform_rt_yz(p1, n1, row2, row3, Tsg);

		// invert Tsg : We will need this
		
		
		//compute_transform_rt(psr, nsr, row1, row2, row3, Tsg);

		// This is a later issue: We might want to look into a local neighborhood only
		// flannIndex->radiusSearch(data, indices, dists, radius, searchParams);

		for (j = 0; j < pc.rows; j ++)
		{
			if (i!=j)
			{
				float* f2 = (float*)(&pc.data[j * pc.step]);
				const double p2[4] = {f2[0], f2[1], f2[2], 0};
				const double n2[4] = {f2[3], f2[4], f2[5], 0};
				double p2t[4], alpha_scene;
				
				double f[4]={0};
				compute_ppf_features(p1, n1, p2, n2, f);
				unsigned int hash = hash_ppf_simple(f, angleStepRadians, distanceStep);

				// we don't need to call this here, as we already estimate the Tsg from scene reference point
				//double alpha = compute_alpha(p1, n1, p2);
				p2t[1] = Tsg[1] + row2[0] * p2[0] + row2[1] * p2[1] + row2[2] * p2[2];
				p2t[2] = Tsg[2] + row3[0] * p2[0] + row3[1] * p2[1] + row3[2] * p2[2];

				alpha_scene=atan2(-p2t[2], p2t[1]);

				if ( alpha_scene != alpha_scene)
				{
					printf("NaN value!\n");
					//return ;
				}

				if (sin(alpha_scene)*p2t[2]<0.0)
					alpha_scene=-alpha_scene;

				//hashtable_int_insert(hashTable, hash, (void*)corrInd);
				//printf("%f %f %f %f \n", f[0], f[1], f[2], f[3]);
				//printf("%d\n", hash);

				tommy_hashtable_node* node = tommy_hashtable_bucket(ppfModel->hashTable, tommy_inthash_u32(hash));

				int numNodes = 0;
				while (node)
				{
					THash* tData = (THash*) node->data;
					int corrI = (int)tData->i;
					//int corrJ = (int)tData->j;
					int ppfInd = (int)tData->ppfInd;
					//int corrInd = (int)tData->i*sampledStep+j;
					float* ppfCorrScene = (float*)(&ppfModel->PPF.data[ppfInd]);
					double alpha_model = (double)ppfCorrScene[4];
					double alpha = alpha_scene - alpha_model;
					//unsigned int hashInd = corrI*numRefPoints + corrJ;
					

					/*  Map alpha to the indices:
						atan2 generates results in (-pi pi]
						That's why alpha should be in range [-2pi 2pi]
						So the quantization would be :
						numAngles * (alpha+2pi)/(4pi)
					*/
					
					int alpha_index = (int)(numAngles*(alpha + 2*PI) / (4*PI));

					unsigned int accIndex = corrI * numAngles + alpha_index;
#pragma omp atomic
					accumulator[accIndex]++;

					node = node->next;
					numNodes++;
				}
				//tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashNode->id));

				//printf("%d\n", numNodes);
			}
		}

		// now maximize the accumulator

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < numAngles; ++j)
			{
				int accInd = i*numAngles + j;
				const int accVal = accumulator[ accInd ];
				if (accVal > max_votes)
				{
					max_votes = accVal;
					max_votes_i = i;
					max_votes_j = j;
				}
				// Reset accumulator_array for the next set of iterations with a new scene reference point
				accumulator[accInd ] = 0;
			}
		}

		printf("Model Reference: %d, Alpha Index: %d\n", max_votes_i, max_votes_j);

		// TODO : Compute pose
		float* fMax = (float*)(&pc.data[max_votes_i * pc.step]);
		const double pMax[4] = {f1[0], f1[1], f1[2], 0};
		const double nMax[4] = {f1[3], f1[4], f1[5], 0};
		double Tmg[3], pose[4][4];

		// convert alpha_index to alpha
		int alpha_index = max_votes_j;


		//compute_transform_rt_yz(p1, n1, row2, row3, Tmg);




		//PPFPose
	}
}

int main()
{
	int useNormals = 1;
	int withBbox = 1;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);
	//Mat pc = Mat(100,100,CV_32FC1);

	TPPFModelPC* ppfModel = 0;
	Mat PPFMAt = train_pc_ppf(pc, 0.05, 0.05, 30, &ppfModel);

	t_match_pc_ppf(pc, 15, 5, ppfModel);

	//compute_ppf_pc(pc, PPFMAt, const double RelSamplingStep, const double RelativeAngleStep, const double RelativeDistanceStep, TPPFModelPC** Model3D)

	return 0;
}


// test octree point cloud sampling
int main_octree_sampling()
{
	int useNormals = 1;
	int withBbox = 1;
	int withOctree = 0;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);

	float xRange[2], yRange[2], zRange[2];
	compute_obb(pc, xRange, yRange, zRange);

	Mat sampled = sample_pc_octree(pc, xRange, yRange, zRange, 0.05);

#ifdef _MSV_VER
	visualize_pc(sampled, 0, 1, 0, "Point Cloud");
#endif

	return 0;
}

// test octree visualization
int main_octree_vis()
{
	int useNormals = 1;
	int withBbox = 1;
	int withOctree = 0;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);

#ifdef _MSV_VER
	visualize_pc(pc, useNormals, withBbox, withOctree, "Point Cloud");
#endif

	return 0;
}

bool is_in_bbox(const float *point, const float range[2]) 
{
	return 
		point[0] >= range[0] &&
		point[1] >= range[0] &&
		point[2] >= range[0] &&
		point[0] <= range[1] &&
		point[1] <= range[1] &&
		point[2] <= range[1];
}

// test octree
int main_octree()
{
	TOctreeNode oc;
	const float dim = 2;
	float root[3]={0,0,0};
	float dimHalf[3]={dim/2,dim/2,dim/2};

	float range[2]={-0.05, 0.05};

	t_octree_init(&oc, root, dimHalf);

	std::vector<float*> points;
	std::vector<float*> resultsGT, resultsOT;

	// Create a bunch of random points
	const int nPoints = 1 * 1000 * 1000; 
	for(int i=0; i<nPoints; ++i) 
	{
		float px = (0-dim/2) + (dim*rand()) * ((dim/2) / RAND_MAX);
		float py = (0-dim/2) + (dim*rand()) * ((dim/2) / RAND_MAX);
		float pz = (0-dim/2) + (dim*rand()) * ((dim/2) / RAND_MAX);

		float* pt =new float[3];
		pt[0]=px;
		pt[1]=py;
		pt[2]=pz;
		points.push_back(pt);

		bool gt = is_in_bbox(pt, range) ;
		resultsGT.push_back(pt);

		if (gt)
			printf("%f %f %f\n", pt[0], pt[1], pt[2]);

		t_octree_insert(&oc, pt);
	}

	// now compute octree results

	printf("\n----------------------------------------------------------------\n");

	t_octree_query_in_bbox(&oc, range, range, range, resultsOT);

	for(unsigned int i=0; i<resultsOT.size(); ++i) 
	{
		float* rot = resultsOT[i];
		printf("%f %f %f\n", rot[0], rot[1], rot[2]);
	}



	return 0;
}

// test bounding box  _bbox
int main_bbox()
{
	int useNormals = 1;
	int withBbox = 1;
	int numVert = 176920;
	const char* fn = "../../../data/cheff2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);

	float xRange[2], yRange[2], zRange[2];
	compute_obb(pc, xRange, yRange, zRange);

	printf("Bounding box -- x: (%f, %f), y: (%f, %f), z: (%f, %f)\n", xRange[0], xRange[1], yRange[0], yRange[1], zRange[0], zRange[1]);

#ifdef _MSV_VER
	visualize_pc(pc, useNormals, withBbox, 0, "Point Cloud");
#endif

	return 0;
}

// test point cloud reading & visualization
// for now uniform sampling looks better
int main_ply()
{
	int useNormals = 1;
	int withBbox = 0;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";

	Mat pc = load_ply_simple(fn, numVert, useNormals);
	
	Mat spc1 = sample_pc_uniform(pc, numVert/350);
	Mat spc2 = sample_pc_random(pc, 350);

#ifdef _MSV_VER
	visualize_pc(pc, useNormals, withBbox, 0, "PC1");
	visualize_pc(spc1, useNormals, withBbox, 0, "Uniform Sampled");
	visualize_pc(spc2, useNormals, withBbox, 0, "Random Sampled");
#endif

	return 0;
}
    
// test hash table
/*int main_hash_test()
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
}*/
