/*#include "opencv2/features2d.hpp"
#include "opencv2/nonfree.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc.hpp"*/
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
//#include <opencv2/rgbd/rgbd.hpp>
#include <opencv2/flann/flann.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "helpers.h"
#include "visualize_win.h"
#include "c_utils.h"
#include "hash_murmur.h"
#include "THashInt.h"

using namespace cv;

#define USE_TOMMY_HASHTABLE

#if defined( USE_TOMMY_HASHTABLE )
#include <tommy.h>
#endif

#if defined (T_OPENMP)
#include <omp.h>
#endif

typedef struct THash {
	int id;
	int i, ppfInd;
	#if defined( USE_TOMMY_HASHTABLE )
	tommy_node node;
	#endif
} THash;

class PPFPose
{
public:
	PPFPose()
	{
		alpha=0; 
		modelIndex=0; 
		numVotes=0;

		for (int i=0; i<16; i++)
			Pose[i]=0;
	};

	PPFPose(double Alpha =0, unsigned int ModelIndex=0, unsigned int NumVotes=0)
	{
		alpha = Alpha;
		modelIndex = ModelIndex;
		numVotes = NumVotes;

		for (int i=0; i<16; i++)
			Pose[i]=0;
	};

	void update_pose(double NewPose[16])
	{
		double R[9];

		for (int i=0; i<16; i++)
			Pose[i]=NewPose[i];

		R[0] = Pose[0];	R[1] = Pose[1]; R[2] = Pose[2];
		R[3] = Pose[4];	R[4] = Pose[5]; R[5] = Pose[6];
		R[6] = Pose[8];	R[7] = Pose[9]; R[8] = Pose[10];

		t[0]=Pose[3]; t[1]=Pose[7]; t[2]=Pose[11];

		// compute the angle
		const double trace = R[0] + R[4] + R[8];

		if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
		else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
		else							 {angle = ( acos((trace - 1)/2) ); }

		// compute the quaternion
		matrix_to_quaternion(R, q);
	}

	void update_pose(double NewR[9], double NewT[3])
	{
		//for (int i=0; i<16; i++)
		//	Pose[i]=NewPose[i];
		Pose[0]=NewR[0]; Pose[1]=NewR[1]; Pose[2]=NewR[2]; Pose[3]=NewT[0];
		Pose[4]=NewR[3]; Pose[5]=NewR[4]; Pose[6]=NewR[5]; Pose[7]=NewT[1];
		Pose[8]=NewR[6]; Pose[9]=NewR[7]; Pose[10]=NewR[8]; Pose[11]=NewT[2];
		Pose[12]=0;		 Pose[13]=0;	  Pose[14]=0;	   Pose[15]=1;

		// compute the angle
		const double trace = NewR[0] + NewR[4] + NewR[8];

		if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
		else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
		else							 {angle = ( acos((trace - 1)/2) ); }

		// compute the quaternion
		matrix_to_quaternion(NewR, q);
	}

	void update_pose_quat(double Q[4], double NewT[3])
	{
		double NewR[9];

		quaternion_to_matrix(Q, NewR);
		q[0]=Q[0]; q[1]=Q[1]; q[2]=Q[2]; q[3]=Q[3]; 

		//for (int i=0; i<16; i++)
		//	Pose[i]=NewPose[i];
		Pose[0]=NewR[0]; Pose[1]=NewR[1]; Pose[2]=NewR[2]; Pose[3]=NewT[0];
		Pose[4]=NewR[3]; Pose[5]=NewR[4]; Pose[6]=NewR[5]; Pose[7]=NewT[1];
		Pose[8]=NewR[6]; Pose[9]=NewR[7]; Pose[10]=NewR[8]; Pose[11]=NewT[2];
		Pose[12]=0;		 Pose[13]=0;	  Pose[14]=0;	   Pose[15]=1;

		// compute the angle
		const double trace = NewR[0] + NewR[4] + NewR[8];

		if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
		else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
		else							 {angle = ( acos((trace - 1)/2) ); }
	}

	PPFPose* clone()
	{
		PPFPose* pose = new PPFPose(alpha, modelIndex, numVotes);
		for (int i=0; i<16; i++)
			pose->Pose[i]= Pose[i];

		pose->q[0]=q[0];
		pose->q[1]=q[1];
		pose->q[2]=q[2];
		pose->q[3]=q[3];

		pose->t[0]=t[0];
		pose->t[1]=t[1];
		pose->t[2]=t[2];

		pose->angle=angle;

		return pose;
	}

	int write_pose(FILE* f)
	{
		int POSE_MAGIC = 7673;
		fwrite(&POSE_MAGIC, sizeof(int), 1, f);
		fwrite(&angle, sizeof(double), 1, f);
		fwrite(&numVotes, sizeof(int), 1, f);
		fwrite(&modelIndex, sizeof(int), 1, f);
		fwrite(Pose, sizeof(double)*16, 1, f);
		fwrite(t, sizeof(double)*3, 1, f);
		fwrite(q, sizeof(double)*4, 1, f);
		return 0;
	}

	int read_pose(FILE* f)
	{
		int POSE_MAGIC = 7673, magic;
		
		fread(&magic, sizeof(int), 1, f);
		if (magic == POSE_MAGIC)
		{
			fread(&angle, sizeof(double), 1, f);
			fread(&numVotes, sizeof(int), 1, f);
			fread(&modelIndex, sizeof(int), 1, f);
			fread(Pose, sizeof(double)*16, 1, f);
			fread(t, sizeof(double)*3, 1, f);
			fread(q, sizeof(double)*4, 1, f);
			return 0;
		}
		
		return -1;
	}

	int write_pose(const char* FileName)
	{
		FILE* f = fopen(FileName, "wb");

		if (!f)
			return -1;

		int status = write_pose(f);

		fclose(f);
		return status;
	}

	int read_pose(const char* FileName)
	{
		FILE* f = fopen(FileName, "rb");

		if (!f)
			return -1;

		int status = read_pose(f);
		
		fclose(f);
		return status;
	}
	
	~PPFPose(){};

	double alpha;
	unsigned int modelIndex;
	unsigned int numVotes;
	double Pose[16], angle, t[3], q[4];
};

class PoseCluster
{
public:
	PoseCluster() 
	{
		poseList.clear();
		numVotes=0;
		id=0;
	};
	
	PoseCluster(PPFPose* newPose) 
	{
		poseList.clear();
		poseList.push_back(newPose);
		numVotes=newPose->numVotes;
		id=0;
	};

	PoseCluster(PPFPose* newPose, int id) 
	{
		poseList.clear();
		poseList.push_back(newPose);
		this->numVotes = newPose->numVotes;
		this->id = id;
	};

	~PoseCluster()
	{
		numVotes=0;
		id=0;
		poseList.clear();
	};

	void add_pose(PPFPose* newPose) 
	{
		poseList.push_back(newPose);
		this->numVotes += newPose->numVotes;
	};

	vector < PPFPose* > poseList;
	int numVotes;
	int id;
};

class TPPFModelPC
{
public:

	TPPFModelPC() {} ;
	~TPPFModelPC() {};

	int magic;
	double maxDist, angleStep, distStep;
	double sampling_step_relative;
	Mat inputPC, sampledPC, PPF;
	cvflann::Index<Distance_32F>* flannIndex;
	int n, numRefPoints, sampledStep, ppfStep;
#if defined (USE_TOMMY_HASHTABLE)
	tommy_hashtable* hashTable;
#else
	hashtable_int* hashTable;
#endif
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
	TAngle3(n1, n2, c, f[2]);

	// We need angles in range [0;pi]
	// TAngle3 already provides that.
	
	// If not, revert to range [0;pi]
	/*if (f[0]<0)
		f[0]+=M_PI;

	if (f[1]<0)
		f[1]+=M_PI;

	if (f[2]<0)
		f[2]+=M_PI;*/
}


// simple hashing
int hash_ppf_simple(const double f[4], const double AngleStep, const double DistanceStep)
{
	const int d1 = (int) (floor ((double)f[0] / (double)AngleStep));
	const int d2 = (int) (floor ((double)f[1] / (double)AngleStep));
	const int d3 = (int) (floor ((double)f[2] / (double)AngleStep));
	const int d4 = (int) (floor ((double)f[3] / (double)DistanceStep));

	//printf("%d, %d, %d,%d\n",d1,d2,d3,d4);
	
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
	MurmurHash3_x86_32(key, 4*sizeof(int), 42, &hashKey);
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

	// checked row2, row3: They are correct
	
	mpt[1] = Tmg[1] + row2[0] * p2[0] + row2[1] * p2[1] + row2[2] * p2[2];
    mpt[2] = Tmg[2] + row3[0] * p2[0] + row3[1] * p2[1] + row3[2] * p2[2];

	alpha=atan2(-mpt[2], mpt[1]);

	if ( alpha != alpha)
    {
		//printf("NaN value!\n");
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
	compute_bbox_std(PC, xRange, yRange, zRange);

	// compute sampling step from diameter of bbox
	float dx = xRange[1] - xRange[0];
	float dy = yRange[1] - yRange[0];
	float dz = zRange[1] - zRange[0];
	float diameter = sqrt ( dx * dx + dy * dy + dz * dz );

	//float diameter = compute_diameter(pc);

	float distanceStep = diameter * sampling_step_relative;

	//Mat sampled = sample_pc_octree(PC, xRange, yRange, zRange, sampling_step_relative);
	Mat sampled = sample_pc_by_quantization(PC, xRange, yRange, zRange, sampling_step_relative);
	
	//visualize_pc(sampled, 1, 1, 0, "reg");

    double angleStepRadians = (360.0/angle_step_relative)*PI/180.0;

	int size = next_power_of_two(sampled.rows*sampled.rows);

#if defined (USE_TOMMY_HASHTABLE)
	tommy_hashtable* hashTable = (tommy_hashtable*)malloc(sizeof(tommy_hashtable));
	// 262144 = 2^18	
	tommy_hashtable_init(hashTable, size);
#else
	hashtable_int* hashTable = hashtable_int_create(size, NULL);
#endif
	
	*Model3D = new TPPFModelPC();

	int numPPF = sampled.rows*sampled.rows;
	(*Model3D)->PPF = Mat(numPPF, T_PPF_LENGTH, CV_32FC1);
	int ppfStep = (*Model3D)->PPF.step;
	int sampledStep = sampled.step;
	
	// TODO: Maybe I could sample 1/5th of them here. Check the performance later.
	int numRefPoints = sampled.rows;
	for (int i=0; i<numRefPoints; i++)
	{
		float* f1 = (float*)(&sampled.data[i * sampledStep]);
		const double p1[4] = {f1[0], f1[1], f1[2], 0};
		const double n1[4] = {f1[3], f1[4], f1[5], 0};

		//printf("///////////////////// NEW REFERENCE ////////////////////////\n");
		for (int j=0; j<numRefPoints; j++)
		{
			// cannnot compute the ppf with myself
			if (i!=j)
			{
				float* f2 = (float*)(&sampled.data[j * sampledStep]);
				const double p2[4] = {f2[0], f2[1], f2[2], 0};
				const double n2[4] = {f2[3], f2[4], f2[5], 0};

				double f[4]={0};
				compute_ppf_features(p1, n1, p2, n2, f);
				//unsigned int hashValue = hash_ppf_simple(f, angleStepRadians, distanceStep);
				unsigned int hashValue = hash_ppf(f, angleStepRadians, distanceStep);
				double alpha = compute_alpha(p1, n1, p2);
				unsigned int corrInd = i*numRefPoints+j;
				unsigned int ppfInd = corrInd*ppfStep;

				THash* hashNode = (THash*)calloc(1, sizeof(THash));
				hashNode->id = hashValue;
				//hashNode->data = (void*)corrInd;
				hashNode->i = i;
				//hashNode->j = j;
				hashNode->ppfInd = ppfInd;

#if defined(USE_TOMMY_HASHTABLE)
				//tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, (hashValue));
				tommy_hashtable_insert(hashTable, &hashNode->node, hashNode, tommy_inthash_u32(hashValue));
#else
				hashtable_int_insert_hashed(hashTable, hashValue, (void*)hashNode);
#endif

				float* ppfRow = (float*)(&(*Model3D)->PPF.data[ ppfInd ]);
				ppfRow[0] = f[0];
				ppfRow[1] = f[1];
				ppfRow[2] = f[2];
				ppfRow[3] = f[3];
				ppfRow[4] = (float)alpha;
			}
		}

		//printf("Training reference : %d\n", i);
	}

	//*Model3D = (TPPFModelPC*)calloc(1, sizeof(TPPFModelPC));

	(*Model3D)->angleStep = angleStepRadians;
	(*Model3D)->distStep = distanceStep;
	(*Model3D)->hashTable = hashTable;
	(*Model3D)->sampledStep = sampledStep;
	(*Model3D)->ppfStep = ppfStep;
	(*Model3D)->numRefPoints = numRefPoints;
	(*Model3D)->sampling_step_relative = sampling_step_relative;
	(*Model3D)->sampledPC = sampled;
	
	return Mat();
}

Mat t_load_ppf_model(const char* FileName)
{
	Mat ppf = Mat();

	return ppf;
}

bool match_pose(const PPFPose sourcePose, const PPFPose targetPose, const double PositionThrehsold, const double RotationThreshold)
{
	// translational difference
	const double* Pose = sourcePose.Pose;
	const double* PoseT = targetPose.Pose;
	double dv[3] = {targetPose.t[0]-sourcePose.t[0], targetPose.t[1]-sourcePose.t[1], targetPose.t[2]-sourcePose.t[2]};
	double dNorm = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
	
	const double phi = fabs ( sourcePose.angle - targetPose.angle );

	return (phi<RotationThreshold && dNorm < PositionThrehsold);
}

int qsort_pose_cmp (const void * a, const void * b)
{
	PPFPose* pose1 = *(PPFPose**)a;
	PPFPose* pose2 = *(PPFPose**)b;
   return ( pose2->numVotes - pose1->numVotes );
}

int sort_pose_clusters (const PoseCluster* a, const PoseCluster* b)
{
   return ( a->numVotes > b->numVotes );
}

int cluster_poses(PPFPose** poseList, const int numPoses, const double PositionThreshold, const double RotationThreshold, const double MinMatchScore, vector < PPFPose* >& finalPoses)
{
	vector<PoseCluster*> poseClusters;
	poseClusters.clear();

	finalPoses.clear();

	// sort the poses for stability
	qsort(poseList, numPoses, sizeof(PPFPose*), qsort_pose_cmp);

	for (int i=0; i<numPoses; i++)
	{
		PPFPose* pose = poseList[i];
		bool assigned = false;

		// search all clusters
		for (int j=0; j<poseClusters.size() && !assigned; j++)
		{
			const PPFPose* poseCenter = poseClusters[j]->poseList[0];
			if (match_pose(*pose, *poseCenter, PositionThreshold, RotationThreshold))
			{
				poseClusters[j]->add_pose(pose);
				assigned = true;
			}
		}

		if (!assigned)
		{
			poseClusters.push_back ( new PoseCluster(pose));
		}
	}

	// sort the clusters so that we could output multiple hypothesis
	std::sort (poseClusters.begin(), poseClusters.end(), sort_pose_clusters);

	finalPoses.resize(poseClusters.size());

	// TODO: Use MinMatchScore

	for (int i=0; i<poseClusters.size(); i++)
	{
		// We could only average the quaternions. So I will make use of them here		
		double qAvg[4]={0}, tAvg[3]={0}, R[9]={0}, Pose[16]={0};

		// Perform the final averaging
		PoseCluster* curCluster = poseClusters[i];
		vector<PPFPose*> curPoses = curCluster->poseList;
		const int curSize = curPoses.size();

		for (int j=0; j<curSize; j++)
		{
			qAvg[0]+= curPoses[j]->q[0];
			qAvg[1]+= curPoses[j]->q[1];
			qAvg[2]+= curPoses[j]->q[2];
			qAvg[3]+= curPoses[j]->q[3];

			tAvg[0]+= curPoses[j]->t[0];
			tAvg[1]+= curPoses[j]->t[1];
			tAvg[2]+= curPoses[j]->t[2];
		}

		tAvg[0]/=(double)curSize;
		tAvg[1]/=(double)curSize;
		tAvg[2]/=(double)curSize;

		qAvg[0]/=(double)curSize;
		qAvg[1]/=(double)curSize;
		qAvg[2]/=(double)curSize;
		qAvg[3]/=(double)curSize;

		curPoses[0]->update_pose_quat(qAvg, tAvg);
		curPoses[0]->numVotes=curCluster->numVotes;

		//finalPoses.push_back(curPoses[0]);
		finalPoses[i]=curPoses[0]->clone();

		// we won't need this
		delete poseClusters[i];
	}

	poseClusters.clear();

	return 0;
}

void t_match_pc_ppf(Mat pc, const float SearchRadius, const int SampleStep, const float sampling_step_relative, const TPPFModelPC* ppfModel, vector < PPFPose* >& results)
{
	cvflann::Matrix<float> data;
	int i;
	int numNeighbors = 1000;
	int numAngles = (int) (floor (2 * M_PI / ppfModel->angleStep));
	cvflann::Index<Distance_32F>* flannIndex;
	float angleStepRadians = ppfModel->angleStep;
	float distanceStep = ppfModel->distStep;
	int sampledStep = ppfModel->sampledStep;
	int ppfStep = ppfModel->ppfStep;
	int numRefPoints = ppfModel->numRefPoints;
	unsigned int n = numRefPoints;
	PPFPose** poseList;
	int sceneSamplingStep = SampleStep, c = 0;

	// compute bbox
	float xRange[2], yRange[2], zRange[2];
	compute_bbox_std(pc, xRange, yRange, zRange);
	//float diameter = compute_diameter(pc);
	// sample the point cloud
	//float sampling_step_relative = (float)ppfModel->sampling_step_relative;
	float dx = xRange[1] - xRange[0];
	float dy = yRange[1] - yRange[0];
	float dz = zRange[1] - zRange[0];
	float diameter = sqrt ( dx * dx + dy * dy + dz * dz );
	float distanceSampleStep = diameter * sampling_step_relative;
	//Mat sampled = sample_pc_octree(pc, xRange, yRange, zRange, sampling_step_relative);
	Mat sampled = sample_pc_by_quantization(pc, xRange, yRange, zRange, sampling_step_relative);
	//Mat sampled = pc.clone();

	//visualize_pc(sampled, 1, 1, 0, "reg");

	// allocate the accumulator
#if !defined (T_OPENMP)
	unsigned int* accumulator = (unsigned int*)calloc(numAngles*n, sizeof(unsigned int));
#endif

	poseList = (PPFPose**)calloc((sampled.rows/sceneSamplingStep)+4, sizeof(PPFPose*));

	// obtain the tree representation for fast search
	//flannIndex  = (cvflann::Index<Distance_32F>*)index_pc_flann(sampled, data);

	//cv::Mat1i ind(sampled.rows, numNeighbors);
	//cvflann::Matrix<int> indices((int*) ind.data, ind.rows, ind.cols);
	//cvflann::Matrix<float> dists(new float[sampled.rows*numNeighbors], sampled.rows, numNeighbors);

	// TODO: Can be parallelized!
#if defined T_OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < sampled.rows; i += sceneSamplingStep)
	{
		unsigned int refIndMax = 0, alphaIndMax = 0;
		unsigned int maxVotes = 0;

		int j;

		float* f1 = (float*)(&sampled.data[i * sampled.step]);
		const double p1[4] = {f1[0], f1[1], f1[2], 0};
		const double n1[4] = {f1[3], f1[4], f1[5], 0};
		double p1t[4];
		double *row1, *row2, *row3, tsg[3]={0}, Rsg[9]={0}, RInv[9]={0};

#if defined (T_OPENMP)
		unsigned int* accumulator = (unsigned int*)calloc(numAngles*n, sizeof(unsigned int));
		//printf("In thread: %d\n", omp_get_thread_num());
#endif

		//compute_transform_rt_yz(p1, n1, row2, row3, tsg);
		compute_transform_rt(p1, n1, Rsg, tsg);
		row1=&Rsg[0]; row2=&Rsg[3]; row3=&Rsg[6];

		// This is a later issue: We might want to look into a local neighborhood only
		// flannIndex->radiusSearch(data, indices, dists, radius, searchParams);

		for (j = 0; j < sampled.rows; j ++)
		{
			if (i!=j)
			{
				float* f2 = (float*)(&sampled.data[j * sampled.step]);
				const double p2[4] = {f2[0], f2[1], f2[2], 0};
				const double n2[4] = {f2[3], f2[4], f2[5], 0};
				double p2t[4], alpha_scene;
				
				double f[4]={0};
				compute_ppf_features(p1, n1, p2, n2, f);
				//unsigned int hashValue = hash_ppf_simple(f, angleStepRadians, distanceStep);
				unsigned int hashValue = hash_ppf(f, angleStepRadians, distanceStep);

				// we don't need to call this here, as we already estimate the tsg from scene reference point
				//double alpha = compute_alpha(p1, n1, p2);
				p2t[1] = tsg[1] + row2[0] * p2[0] + row2[1] * p2[1] + row2[2] * p2[2];
				p2t[2] = tsg[2] + row3[0] * p2[0] + row3[1] * p2[1] + row3[2] * p2[2];

				alpha_scene=atan2(-p2t[2], p2t[1]);

				if ( alpha_scene != alpha_scene)
				{
					continue;
				}

				if (sin(alpha_scene)*p2t[2]<0.0)
					alpha_scene=-alpha_scene;

				alpha_scene=-alpha_scene;


#if defined (USE_TOMMY_HASHTABLE)
				//tommy_hashtable_node* node = tommy_hashtable_bucket(ppfModel->hashTable, (hashValue));
				tommy_hashtable_node* node = tommy_hashtable_bucket(ppfModel->hashTable, tommy_inthash_u32(hashValue));
#else
				hashnode_i* node = hashtable_int_get_bucket_hashed(ppfModel->hashTable, (hashValue));
#endif

				//int numNodes = 0;
				while (node)
				{
					THash* tData = (THash*) node->data;
					int corrI = (int)tData->i;
					int ppfInd = (int)tData->ppfInd;
					float* ppfCorrScene = (float*)(&ppfModel->PPF.data[ppfInd]);
					double alpha_model = (double)ppfCorrScene[T_PPF_LENGTH-1];
					double alpha = alpha_model - alpha_scene;
					

					/*  Map alpha to the indices:
						atan2 generates results in (-pi pi]
						That's why alpha should be in range [-2pi 2pi]
						So the quantization would be :
						numAngles * (alpha+2pi)/(4pi)
					*/

					//printf("%f\n", alpha);					
					int alpha_index = (int)(numAngles*(alpha + 2*PI) / (4*PI));

					unsigned int accIndex = corrI * numAngles + alpha_index;

					accumulator[accIndex]++;
					node = node->next;

					//numNodes++;
				}
				//printf("%d\n", numNodes);
			}
		}

		// Maximize the accumulator
		for (int k = 0; k < n; k++)
		{
			for (int j = 0; j < numAngles; j++)
			{
				const unsigned int accInd = k*numAngles + j;
				const unsigned int accVal = accumulator[ accInd ];
				if (accVal > maxVotes)
				{
					maxVotes = accVal;
					refIndMax = k;
					alphaIndMax = j;
				}
				
#if !defined (T_OPENMP)
				accumulator[accInd ] = 0;
#endif
			}
		}

		// invert Tsg : Luckily rotation is orthogonal: Inverse = Transpose.
		// We are not required to invert.
		double tInv[3], tmg[3], Rmg[9], Ralpha[9];
		matrix_transpose33(Rsg, RInv);
		matrix_product331(RInv, tsg, tInv);

		double TsgInv[16] =	{	RInv[0], RInv[1], RInv[2], -tInv[0],
								RInv[3], RInv[4], RInv[5], -tInv[1],
								RInv[6], RInv[7], RInv[8], -tInv[2],
								0, 0, 0, 1
								};

		// TODO : Compute pose
		const float* fMax = (float*)(&ppfModel->sampledPC.data[refIndMax * ppfModel->sampledPC.step]);
		const double pMax[4] = {fMax[0], fMax[1], fMax[2], 1};
		const double nMax[4] = {fMax[3], fMax[4], fMax[5], 1};
		double pose[4][4];

		compute_transform_rt(pMax, nMax, Rmg, tmg);
		row1=&Rsg[0]; row2=&Rsg[3]; row3=&Rsg[6];

		double Tmg[16] =	{	Rmg[0], Rmg[1], Rmg[2], tmg[0],
								Rmg[3], Rmg[4], Rmg[5], tmg[1],
								Rmg[6], Rmg[7], Rmg[8], tmg[2],
								0, 0, 0, 1
							};

		// convert alpha_index to alpha
		int alpha_index = alphaIndMax;
		double alpha = (alpha_index*(4*PI))/numAngles-2*PI;
		
		// Equation 2:
		double Talpha[16]={0};
		get_unit_x_rotation_44(alpha, Talpha);
		
		double Temp[16]={0};
		double Pose[16]={0};
		matrix_product44(Talpha, Tmg, Temp);
		matrix_product44(TsgInv, Temp, Pose);

		PPFPose *ppf = new PPFPose(alpha, refIndMax, maxVotes);
		
		ppf->update_pose(Pose);

		poseList[i/sceneSamplingStep] = ppf;

//		printf("Model Reference: %d, Alpha Index: %d, Alpha: %f\n", refIndMax, alphaIndMax, alpha);

#if defined (T_OPENMP)
		free(accumulator);
#endif
	}

	// TODO : Make the parameters relative if not arguments.
	double RotationThreshold = ((360/ppfModel->angleStep) / 180.0 * M_PI);
	double PositionThreshold = ppfModel->sampling_step_relative;
	double MinMatchScore = 0.5;
	
	int numPosesAdded = sampled.rows/sceneSamplingStep;
	cluster_poses(poseList, numPosesAdded, PositionThreshold, RotationThreshold, MinMatchScore, results);

	// free up the used space
	sampled.release();

	for (int i=0; i<numPosesAdded; i++)
	{
		PPFPose* pose = poseList[i];
		delete pose;
		poseList[i]=0;
	}

	free(poseList);
#if !defined (T_OPENMP)
		free(accumulator);
#endif
}


// real data
int main()
{
	int useNormals = 1;
	int withBbox = 1;
	int numVert = 16399;
	const char* fn = "../../../data/cheff_simple.ply";
	//int numVert = 135142;
	//const char* fn = "../../../data/chicken_high.ply";
	//int numVert = 28291;
	//const char* fn = "../../../data/parasaurolophus_low_normals2.ply";
	//int numVert = 6700;
	//const char* fn = "../../../data/parasaurolophus_6700.ply";
	//int numVert = 51954;
	//const char* fn = "../../../data/SpaceTime/Scena1/scene1-model1_0_ascii.ply";
	//int numVert = 33368;
	//const char* fn = "../../../data/chicken2.ply";
	//int numVert = 13550;
	//const char* fn = "../../../data/chicken3.ply";
	
	Mat pc = load_ply_simple(fn, numVert, useNormals);
	//Mat pc = Mat(100,100,CV_32FC1);

	TPPFModelPC* ppfModel = 0;
	printf("Training...");
	Mat PPFMAt = train_pc_ppf(pc, 0.05, 0.05, 30, &ppfModel);
	printf("\nTraining complete. Loading model...");

	//numVert = 122503;
	//fn = "../../../data/SpaceTime/Scena1/scene1-scene4_0_ascii.ply";
	//numVert = 113732;
	//fn = "../../../data/Retrieval/rs22_proc2.ply";
	//numVert = 114373;
	//fn = "../../../data/rs1_normals.ply";
	numVert = 135985;
	fn = "../../../data/Retrieval/rs8_proc.ply";
	//numVert = 12345;
	//fn = "../../../data/rs1_normals2.ply";
	Mat pcTest = load_ply_simple(fn, numVert, useNormals);
	printf("\nStarting matching...");

	int64 tick1 = cv::getTickCount();
	vector < PPFPose* > results;
	t_match_pc_ppf(pcTest, 15, 5, 0.05, ppfModel, results);
	int64 tick2 = cv::getTickCount();
	printf("Elapsed Time %f sec\n", (double)(tick2-tick1)/ cv::getTickFrequency());

	printf("Estimated Poses:\n");

	// debug first five poses
	for (int i=0; i<MIN(5, results.size()); i++)
	{
		PPFPose* pose = results[i];

		// Print the pose
		printf("Pose %d : Voted by %d, Alpha is %f\n", i, pose->numVotes, pose->alpha);
		for (int j=0; j<4; j++)
		{
			for (int k=0; k<4; k++)
			{
				printf("%f ", pose->Pose[j*4+k]);
			}
			printf("\n");
		}
		printf("\n");

		// Visualize registration
		Mat pct = transform_pc_pose(pc, pose->Pose);
#ifdef _MSC_VER
		visualize_registration(pcTest, pct, "Registration");
#endif
	}

	return 0;
}

//
int main_synthetic()
{
	int useNormals = 1;
	int withBbox = 1;
	//int numVert = 6700;
	//const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	int numVert = 34834;
	const char* fn = "../../../data/Stanford/bunny/bun_zipper_ascii.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);
	//Mat pc = Mat(100,100,CV_32FC1);

	TPPFModelPC* ppfModel = 0;
	printf("Training...");
	Mat PPFMAt = train_pc_ppf(pc, 0.05, 0.03, 30, &ppfModel);
	printf("\nTraining complete. Starting matching...\n");

	// generate 10 poses and compare
	for (int numTrials =0 ; numTrials <10; numTrials++)
	{
		// Make a sample pose:
		double Pose[16]={0};
		//generate_random_pose(Pose, 0.5);
		get_random_pose(Pose);
		printf("Random Pose (Ground Truth):\n");
		matrix_print(Pose, 4,4);	

		// add noise and transform
		Mat pcPerturb = add_noise_pc(pc, 0.000015);
		//Mat pcPerturb = pc.clone();

		// observe a partial object: Occlusion simulation
		Mat pcPartial = pcPerturb.rowRange(pcPerturb.rows/3, pcPerturb.rows);

		// transform with arbitrary pose
		Mat pcPerturbTrans = transform_pc_pose(pcPartial, Pose );
		//Mat pcPerturb = pc.clone();

		int64 tick1 = cv::getTickCount();
		vector < PPFPose* > results;
		t_match_pc_ppf(pcPerturbTrans, 15, 5, 0.05, ppfModel, results);
		int64 tick2 = cv::getTickCount();
		printf("Elapsed Time %f sec\n", (double)(tick2-tick1)/ cv::getTickFrequency());

		printf("Estimated Poses (Ground Truth):\n");

		// debug first five poses
		for (int i=0; i<MIN(5, results.size()); i++)
		{
			PPFPose* pose = results[i];

			// Print the pose
			printf("Pose %d : Voted by %d, Alpha is %f\n", i, pose->numVotes, pose->alpha);
			for (int j=0; j<4; j++)
			{
				for (int k=0; k<4; k++)
				{
					printf("%f ", pose->Pose[j*4+k]);
				}
				printf("\n");
			}
			printf("\n");

			// Visualize registration
			Mat pct = transform_pc_pose(pc, pose->Pose);
#ifdef _MSC_VER
			visualize_registration(pcPerturbTrans, pct, "Registration");
#endif
		}

	}
	/*for (int i=0; i<MIN(3, results.size()); i++)
	{
		PPFPose* pose = results[i];

		printf("Pose %d : Voted by %d, Alpha is %f\n", i, pose->numVotes, pose->alpha);
		for (int j=0; j<4; j++)
		{
			for (int k=0; k<4; k++)
			{
				printf("%f ", pose->Pose[j*4+k]);
			}
			printf("\n");
		}
		printf("\n");
	}*/

	//compute_ppf_pc(pc, PPFMAt, const double RelSamplingStep, const double RelativeAngleStep, const double RelativeDistanceStep, TPPFModelPC** Model3D)

	return 0;
}

// rotation test
int main_rt()
{
	double r[3]={1,0,0};
	double angle = 0.319;
	double R[16]={0};
	get_unit_x_rotation_44(angle, R);
	for (int i=0; i<16; i++)
		printf("%f,", R[i]);

	printf("\n\n");
	memset(R,0,sizeof(double)*16);
	//axis_angle_to_matrix(r, angle, R);
	aa_to_R(angle, r, R);
	for (int i=0; i<9; i++)
		printf("%f,", R[i]);

	return 0;
}


// test different point cloud sampling
// _sampling
int main_sampling()
{
	int useNormals = 1;
	int withBbox = 1;
	int withOctree = 0;
	int numVert = 6700;
	const char* fn = "../../../data/parasaurolophus_6700_2.ply";
	Mat pc = load_ply_simple(fn, numVert, useNormals);

	float xRange[2], yRange[2], zRange[2];
	compute_bbox_std(pc, xRange, yRange, zRange);

	float dx = xRange[1] - xRange[0];
	float dy = yRange[1] - yRange[0];
	float dz = zRange[1] - zRange[0];
	float diameter = sqrt ( dx * dx + dy * dy + dz * dz );

#if !defined (T_OPENMP)
	omp_set_num_threads(8);
#endif

	int64 t1 = cv::getTickCount();
	Mat sampled2 = sample_pc_octree(pc, xRange, yRange, zRange, 0.01);
	int64 t2 = cv::getTickCount();
	printf("Elapsed : %f\n", (t2-t1)/cv::getTickFrequency());

	t1 = cv::getTickCount();
	Mat sampled = sample_pc_by_quantization(pc, xRange, yRange, zRange, 0.01);
	t2 = cv::getTickCount();
	printf("Elapsed : %f\n", (t2-t1)/cv::getTickFrequency());
	//Mat sampled = sample_pc_kd_tree(pc, diameter*0.025, 15);

#ifdef _MSC_VER
	visualize_pc(sampled, 0, 1, 0, "Point Cloud");
#endif

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
	compute_bbox_std(pc, xRange, yRange, zRange);

	Mat sampled = sample_pc_octree(pc, xRange, yRange, zRange, 0.05);

#ifdef _MSC_VER
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

#ifdef _MSC_VER
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
	compute_bbox_std(pc, xRange, yRange, zRange);

	printf("Bounding box -- x: (%f, %f), y: (%f, %f), z: (%f, %f)\n", xRange[0], xRange[1], yRange[0], yRange[1], zRange[0], zRange[1]);

#ifdef _MSC_VER
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

#ifdef _MSC_VER
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

