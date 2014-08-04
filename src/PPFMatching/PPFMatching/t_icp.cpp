
// Author: Tolga Birdal

#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <math.h>
#include "helpers.h"
#include "c_utils.h"
#include "t_icp.h"
#include "THashInt.h"
#include "visualize_win.h"

using namespace std;
using namespace cv;

void subtract_columns(Mat SrcPC, double mean[3])
{
	int height = SrcPC.rows;
	int width = SrcPC.cols;

	for (int i=0; i<height; i++)
	{
		float *row = (float*)(&SrcPC.data[i*SrcPC.step]);
		{
			row[0]-=(float)mean[0];
			row[1]-=(float)mean[1];
			row[2]-=(float)mean[2];
		}
	}
}


// as in PCA
void compute_mean_cols(Mat SrcPC, double mean[3])
{
	int height = SrcPC.rows;
	int width = SrcPC.cols;

	double mean1=0, mean2 = 0, mean3 = 0;

	for (int i=0; i<height; i++)
	{
		const float *row = (float*)(&SrcPC.data[i*SrcPC.step]);
		{
			mean1 += (double)row[0];
			mean2 += (double)row[1];
			mean3 += (double)row[2];
		}
	}

	mean1/=(double)height;
	mean2/=(double)height;
	mean3/=(double)height;

	mean[0] = mean1;
	mean[1] = mean2;
	mean[2] = mean3;
}

// as in PCA
void subtract_mean_from_columns(Mat SrcPC, double mean[3])
{
	compute_mean_cols(SrcPC, mean);
	subtract_columns(SrcPC, mean);
}

double compute_dist_to_origin(Mat SrcPC)
{
	int height = SrcPC.rows;
	int width = SrcPC.cols;
	double dist = 0;

	for (int i=0; i<height; i++)
	{
		const float *row = (float*)(&SrcPC.data[i*SrcPC.step]);
		dist += sqrt(row[0]*row[0]+row[1]*row[1]+row[2]*row[2]);
	}

	return dist;
}


#define ELEM_SWAP_F(a,b) { float temp=(a);(a)=(b);(b)=temp; }
static float median_F(float arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) >>1;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP_F(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) >>1;
    if (arr[middle] > arr[high])    ELEM_SWAP_F(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP_F(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP_F(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP_F(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP_F(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP_F(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

float madsigma(float* r, int m)
{
	float* t=(float*)calloc(m, sizeof(float));
	int i=0;
	float s=0, medR;

	memcpy(t, r, m*sizeof(float));
	medR=median_F(t, m);

	for (i=0; i<m; i++)
		t[i] = fabs(r[i]-medR);

	s= median_F(t, m)/ 0.6745;

	free(t);
	return s;
}

float get_rejection_threshold(float* r, int m, float outlierScale)
{
	float* t=(float*)calloc(m, sizeof(float));
	int i=0;
	float s=0, medR, threshold;

	memcpy(t, r, m*sizeof(float));
	medR=median_F(t, m);

	for (i=0; i<m; i++)
		t[i] = (float)fabs((double)r[i]-(double)medR);

	s = 1.48257968 * median_F(t, m);

	threshold = (outlierScale*s+medR);

	free(t);
	return threshold;
}

// Kok Lim Low's linearization 
void minimize_point_to_plane_metric(Mat Src, Mat Dst, Mat& X)
{
	//Mat sub = Dst - Src;
	Mat A=Mat(Src.rows, 6, CV_64F);
	Mat b=Mat(Src.rows, 1, CV_64F);

#if defined T_OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<Src.rows; i++)
	{
		const double *srcPt = (double*)&Src.data[i*Src.step];
		const double *dstPt = (double*)&Dst.data[i*Dst.step];
		const double *normals = &dstPt[3];
		double *bVal = (double*)&b.data[i*b.step];
		double *aRow = (double*)&A.data[i*A.step];

		const double sub[3]={dstPt[0]-srcPt[0], dstPt[1]-srcPt[1], dstPt[2]-srcPt[2]};

		*bVal = TDot3(sub, normals);
		TCross(srcPt, normals, aRow);

		aRow[3] = normals[0];
		aRow[4] = normals[1];
		aRow[5] = normals[2];
	}

	cv::solve(A, b, X, DECOMP_SVD);
}


void get_transform_mat(Mat X, double Pose[16])
{
	Mat DCM;
	double *r1, *r2, *r3;
	double* x = (double*)X.data;
	
	const double sx = sin(x[0]);
	const double cx = cos(x[0]);
	const double sy = sin(x[1]);
	const double cy = cos(x[1]);
	const double sz = sin(x[2]);
	const double cz = cos(x[2]);

	Mat R1 = Mat::eye(3,3, CV_64F);
	Mat R2 = Mat::eye(3,3, CV_64F);
	Mat R3 = Mat::eye(3,3, CV_64F);

	r1= (double*)R1.data;
	r2= (double*)R2.data;
	r3= (double*)R3.data;

	r1[4]= cx; r1[5]= -sx;
	r1[7]= sx; r1[8]= cx;

	r2[0]= cy; r2[2]= sy;
	r2[6]= -sy; r2[8]= cy;

	r3[0]= cz; r3[1]= -sz;
	r3[3]= sz; r3[4]= cz;

	DCM = R1*(R2*R3);
	
	Pose[0] = DCM.at<double>(0,0);
	Pose[1] = DCM.at<double>(0,1);
	Pose[2] = DCM.at<double>(0,2);
	Pose[4] = DCM.at<double>(1,0);
	Pose[5] = DCM.at<double>(1,1);
	Pose[6] = DCM.at<double>(1,2);
	Pose[8] = DCM.at<double>(2,0);
	Pose[9] = DCM.at<double>(2,1);
	Pose[10] = DCM.at<double>(2,2);
	Pose[3]=x[3];
	Pose[7]=x[4];
	Pose[11]=x[5];
	Pose[15]=1;
}

/* Fast way to look up the duplicates
   duplicates is pre-allocated
   make sure that the max element in array will not exceed maxElement
*/
hashtable_int* get_hashtable(int* data, int length, int numMaxElement)
{
	hashtable_int* hashtable = hashtable_int_create(numMaxElement*2, 0);
	for(int i = 0; i < length; i++)
	{
		const int key = data[i];

		hashtable_int_insert_hashed(hashtable, key+1, (void*)(i+1));
	}
	
	return hashtable;
}

// source point clouds are assumed to contain their normals
int t_icp_register(const Mat SrcPC, const Mat DstPC, const float Tolerence, const int MaxIterations, const float RejectionScale, const int NumNeighborsCorr, const int NumLevels, const int SampleType, const T_ICP_CALLBACK RegistrationVisFunction, float* Residual, double Pose[16])
{
	int n = SrcPC.rows;
	//double PoseInit[16];

	bool UseRobustReject = RejectionScale>0;

	Mat SrcTemp = SrcPC.clone();
	Mat DstTemp = DstPC.clone();
	double meanSrc[3], meanDst[3];
	compute_mean_cols(SrcTemp, meanSrc);
	compute_mean_cols(DstTemp, meanDst);
	double meanAvg[3]={0.5*(meanSrc[0]+meanDst[0]), 0.5*(meanSrc[1]+meanDst[1]), 0.5*(meanSrc[2]+meanDst[2])};
	subtract_columns(SrcTemp, meanAvg);
	subtract_columns(DstTemp, meanAvg);

	double distSrc = compute_dist_to_origin(SrcTemp);
	double distDst = compute_dist_to_origin(DstTemp);

	double scale = (double)n / ((distSrc + distDst)*0.5);

	SrcTemp(cv::Range(0, SrcTemp.rows), cv::Range(0,3)) = (SrcTemp(cv::Range(0, SrcTemp.rows), cv::Range(0,3)) )  * scale;
	DstTemp(cv::Range(0, DstTemp.rows), cv::Range(0,3)) = (DstTemp(cv::Range(0, DstTemp.rows), cv::Range(0,3)) )  * scale;

	Mat SrcPC0 = SrcTemp;
	Mat DstPC0 = DstTemp;
	
	// initialize pose
	matrix_ident(4, Pose);

	void* flann = index_pc_flann(DstPC0);
	Mat M = Mat::eye(4,4,CV_64F);
	
	// walk the pyramid
	for (int level = NumLevels-1; level >=0; level--)
	{
		const double impact = 2;
		double div = pow((double)2, (double)level);
		double div2 = div*div;
		const int numSamples = round((double)(n/(div)));
		//const double TolP = Tolerence*div2;
		const double TolP = Tolerence*(double)(level+1)*(level+1);
		const int MaxIterationsPyr = round((double)MaxIterations/(level+1));

		// Obtain the sampled point clouds for this level: Also rotates the normals
		Mat SrcPC = transform_pc_pose(SrcPC0, Pose);
		//Mat SrcPC;
		//move_points(SrcPC0, Pose, SrcPC);

		const int sampleStep = round((double)n/(double)numSamples);
		std::vector<int> srcSampleInd;

		/*
			Note by Tolga Birdal
			Downsample the model point clouds. If more optimization is required,
			one could also downsample the scene points, but I think this might 
			decrease the accuracy. That's why I won't be implementing it at this
			moment.

			Also note that you have to compute a KD-tree for each level.
		*/
		SrcPC = sample_pc_uniform_ind(SrcPC, sampleStep, srcSampleInd);

		double fval_old=9999999999;
		double fval_perc=0;
		double fval_min=9999999999;
		Mat Src_Moved = SrcPC.clone();
		
		int i=0;

		int numElSrc = Src_Moved.rows;
		int sizesResult[2] = {numElSrc, 1};
		float* distances = new float[numElSrc];
		int* indices = new int[numElSrc];

		Mat Indices(2, sizesResult, CV_32S, indices, 0);
		Mat Distances(2, sizesResult, CV_32F, distances, 0);

		// use robust weighting for outlier treatment
		int* indicesModel = new int[numElSrc];
		int* indicesScene = new int[numElSrc];

		int* newI = new int[numElSrc];
		int* newJ = new int[numElSrc];

		double PoseX[16]={0};
		matrix_ident(4, PoseX);

		while( (!(fval_perc<(1+TolP) && fval_perc>(1-TolP))) && i<MaxIterationsPyr)
		{
			int selInd = 0, di=0;
			
			query_pc_flann(flann, Src_Moved, Indices, Distances);

			for (di=0; di<numElSrc; di++)
			{
				newI[di]=di;
				newJ[di]=indices[di];
			}
			
			if (UseRobustReject)
			{
				int numInliers = 0;
				float threshold = get_rejection_threshold(distances, Distances.rows, RejectionScale);
				Mat acceptInd = Distances<threshold;
				
				uchar *accPtr = (uchar*)acceptInd.data;
				for (int l=0; l<acceptInd.rows; l++)
				{
					if (accPtr[l])
					{
						newI[numInliers] = l;
						newJ[numInliers] = indices[l];
						numInliers++;
					}
				}	
				numElSrc=numInliers;
			}

			// Step 2: Picky ICP
			// Among the resulting corresponding pairs, if more than one scene point p_i
			// is assigned to the same model point m_j, then select p_i that corresponds
			// to the minimum distance

			hashtable_int* duplicateTable = get_hashtable(newJ, numElSrc, DstPC0.rows);

			for(di=0; di<duplicateTable->size; di++) {
				hashnode_i *node = duplicateTable->nodes[di];

				if (node) {
					// select the first node
					int idx = (int)node->data-1, dn=0;
					int dup = (int)node->key-1;
					int minIdxD = idx;
					float minDist = distances[idx];

					while( node ) {
						idx = (int)node->data-1;

						if (distances[idx] < minDist) {
							minDist = distances[idx];
							minIdxD = idx;
						}

						node = node->next;
						dn++;
					}

					indicesModel[ selInd ] = newI[ minIdxD ];
					indicesScene[ selInd ] = dup ;
					selInd++;
				}
			}

			hashtable_int_destroy(duplicateTable);

			if (selInd)
			{

				Mat Src_Match = Mat(selInd, SrcPC.cols, CV_64F);
				Mat Dst_Match = Mat(selInd, SrcPC.cols, CV_64F);
				
				for(di=0; di<selInd; di++) 
				{
					const int indModel = indicesModel[di];
					const int indScene = indicesScene[di];
					const float *srcPt = (float*)&SrcPC.data[indModel*SrcPC.step];
					const float *dstPt = (float*)&DstPC0.data[indScene*DstPC0.step];
					double *srcMatchPt = (double*)&Src_Match.data[di*Src_Match.step];
					double *dstMatchPt = (double*)&Dst_Match.data[di*Dst_Match.step];
					int ci=0;

					for (ci=0; ci<SrcPC.cols; ci++)
					{
						srcMatchPt[ci] = (double)srcPt[ci];
						dstMatchPt[ci] = (double)dstPt[ci];
					}
				}			

				Mat X;
				minimize_point_to_plane_metric(Src_Match, Dst_Match, X);
				
				get_transform_mat(X, PoseX);
				Src_Moved = transform_pc_pose(SrcPC, PoseX);

				double fval = cv::norm(Src_Match, Dst_Match)/(double)(Src_Moved.rows);

				// Calculate change in error between itterations
				fval_perc=fval/fval_old;

				// Store error value
				fval_old=fval;

				if (fval < fval_min)
					fval_min = fval;
			}
			else
				break;

			//RegistrationVisFunction(PoseX, (void*)SrcPC);

			i++;

			// visualize on demand
#if defined (_MSC_VER)
			//Src_Moved = transform_pc_pose(SrcPC, PoseX);
			//visualize_registration(DstPC0, Src_Moved, "Registration", 150);
#endif
		}

		double TempPose[16];
		matrix_product44(PoseX, Pose, TempPose);

		// no need to copy the last 4 rows	
		for (int c=0; c<12; c++)
			Pose[c] = TempPose[c];

		delete[] newI;
		delete[] newJ;
		delete[] indicesModel;
		delete[] indicesScene;
		delete[] distances;
		delete[] indices;
	}

	// Pose(1:3, 4) = Pose(1:3, 4)./scale;
	Pose[3] = Pose[3]/scale + meanAvg[0];
	Pose[7] = Pose[7]/scale + meanAvg[1];
	Pose[11] = Pose[11]/scale + meanAvg[2];

	//Pose(1:3, 4) = Pose(1:3, 4)./scale + meanAvg' - Pose(1:3, 1:3)*meanAvg';
	double Rpose[9], Cpose[3];
	pose_to_r(Pose, Rpose);
	matrix_product331(Rpose, meanAvg, Cpose);
	Pose[3] -= Cpose[0];
	Pose[7] -= Cpose[1];
	Pose[11] -= Cpose[2];

	destroy_flann(flann); flann = 0;
	return 0;
}