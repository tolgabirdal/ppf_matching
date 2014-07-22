
#include "helpers.h"
#include "c_utils.h"
#include "gdiam.hpp"
#include <time.h>
#include <fstream>
#include <vector>
#include <iostream>

#include "flann/flann.h"
#include "dsyevh3.h"

#if defined (T_OPENMP)
#include <omp.h>
#endif

using namespace std;
using namespace cv;

#define FLANN_ALGORITHM_KD FLANN_INDEX_KDTREE_SINGLE
#define FLANN_CHECKS FLANN_CHECKS_UNLIMITED
#define FLANN_NUM_TREES 8
#define FLANN_NUM_ITERATIONS 12 // not used
#define FLANN_NUM_BRANCHINGS 32 // not used
#define FLANN_TARGET_PRECISION -1
#define FLANN_MAX_LEAF_SIZE 12

Mat load_ply_simple(const char* fileName, int numVertices, int withNormals)
{
	Mat cloud;

	if (withNormals)
		cloud=Mat(numVertices, 6, CV_32FC1);
	else
		cloud=Mat(numVertices, 3, CV_32FC1);

	ifstream ifs(fileName);

	if (!ifs.is_open())
	{
		printf("Cannot open file...\n");
		return Mat();
	}

	string str;
	while (str.substr(0, 10) !="end_header")
		getline(ifs, str);

	float dummy =  0;
	for(size_t i = 0; i < numVertices; i++)
	{
		float* data = (float*)(&cloud.data[i*cloud.step[0]]);
		if (withNormals)
		{
			ifs >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];

			// normalize to unit norm
			double norm = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);
			if (norm>0.00001)
			{
				data[3]/=norm;
				data[4]/=norm;
				data[5]/=norm;
			}
		}
		else
		{
			ifs >> data[0] >> data[1] >> data[2];
		}
	}

	//cloud *= 5.0f;
	return cloud;
}


void write_ply(Mat PC, const char* FileName)
{
    ofstream outFile( FileName );
 
    if ( !outFile )
    {
        //cerr << "Error opening output file: " << FileName << "!" << endl;
		printf("Error opening output file: %s!\n", FileName); 
        exit( 1 );
    }
 
    ////
    // Header
    ////
 
    const int pointNum    = ( int ) PC.rows;
    const int vertNum    = ( int ) PC.cols;
 
    outFile << "ply" << endl;
    outFile << "format ascii 1.0" << endl;
    outFile << "element vertex " << pointNum << endl;
    outFile << "property float x" << endl;
    outFile << "property float y" << endl;
    outFile << "property float z" << endl;
	if(vertNum==6)
	{
		outFile << "property float nx" << endl;
		outFile << "property float ny" << endl;
		outFile << "property float nz" << endl;
	}
    outFile << "end_header" << endl;
 
    ////
    // Points
    ////
 
    for ( int pi = 0; pi < pointNum; ++pi )
    {
        const float* point = (float*)(&PC.data[ pi*PC.step ]);
 
		outFile << point[0] << " "<<point[1]<<" "<<point[2];

		if (vertNum==6)
		{
			outFile<<" " << point[3] << " "<<point[4]<<" "<<point[5];
		}
         
        outFile << endl;
    }
 
    return;
}

TOctreeNode* Mat2Octree(Mat pc)
{
	float xRange[2], yRange[2], zRange[2];
	compute_bbox_std(pc, xRange, yRange, zRange);

	float cx = (xRange[1] + xRange[0])*0.5f;
	float cy = (yRange[1] + yRange[0])*0.5f;
	float cz = (zRange[1] + zRange[0])*0.5f;

	float maxDim = MAX( MAX(xRange[1], yRange[1]), zRange[1]);
	float minDim = MIN( MIN(xRange[1], yRange[1]), zRange[1]);

	float root[3]={cx, cy, cz};
	float half_dim[3]={(xRange[1]-xRange[0])*0.5, (yRange[1]-yRange[0])*0.5, (zRange[1]-zRange[0])*0.5};

	TOctreeNode * oc = new TOctreeNode();
	t_octree_init(oc, root, half_dim);

	int nPoints = pc.rows;
	for(int i=0; i<nPoints; ++i) 
	{
		float* data = (float*)(&pc.data[i*pc.step[0]]);
		t_octree_insert(oc, data);
	}

	return oc;
}

Mat sample_pc_uniform(Mat PC, int sampleStep)
{
	int numRows = PC.rows/sampleStep;
	Mat sampledPC = Mat(numRows, PC.cols, PC.type());

	int c=0;
	for (int i=0; i<PC.rows && c<numRows; i+=sampleStep)
	{
		PC.row(i).copyTo(sampledPC.row(c++));
	}

	return sampledPC;
}

Mat sample_pc_uniform_ind(Mat PC, int sampleStep, vector<int> &indices)
{
	int numRows = round((double)PC.rows/(double)sampleStep);
	indices.resize(numRows);
	Mat sampledPC = Mat(numRows, PC.cols, PC.type());

	int c=0;
	for (int i=0; i<PC.rows && c<numRows; i+=sampleStep)
	{
		indices[c] = i;
		PC.row(i).copyTo(sampledPC.row(c++));
	}

	return sampledPC;
}

// not yet implemented
Mat sample_pc_perfect_uniform(Mat PC, int sampleStep)
{
	// make a symmetric square matrix for sampling
	int numPoints = PC.rows;
	int numPointsSampled = numPoints/sampleStep;

	int n = sqrt((double)numPoints);
	int numTotalPerfect = n*n;

	
	return Mat();
}

/*void* index_pc_flann(Mat pc, cvflann::Matrix<float>& data)
{	
	cvflann::LinearIndexParams params;
	
	data = cvflann::Matrix<float>( (float*)pc.data, pc.rows, 3, pc.step/sizeof(float));
	
	cvflann::Index < Distance_32F>* flannIndex = new cvflann::Index< Distance_32F >(data, params);
	flannIndex->buildIndex();

	return (void*)flannIndex;
}*/

void* index_pc_flann(Mat pc)
{	
	FLANNParameters p;
	FLANN_INDEX flannIndex;
	float speedup=0;
	int i;
	float* dataset;
	bool isCont = pc.isContinuous();

	p.log_level = FLANN_LOG_NONE;
	//p.log_destination = NULL;
	p.algorithm = FLANN_ALGORITHM_KD;
	p.checks = FLANN_CHECKS;
	p.trees = FLANN_NUM_TREES;
	p.branching = FLANN_NUM_BRANCHINGS;
	p.iterations = FLANN_NUM_ITERATIONS;
	p.target_precision = FLANN_TARGET_PRECISION;
	p.leaf_max_size = FLANN_MAX_LEAF_SIZE;
	p.eps = 0;

	flann_set_distance_type(flann::FLANN_DIST_EUCLIDEAN, 0);

#if defined (T_OPENMP)
		p.cores = omp_get_num_threads();
#else
		p.cores = 1;
#endif

	if (isCont && pc.rows==3)
	{
		dataset=(float*)pc.data;
	}
	else
	{
		dataset = new float[pc.rows*3];
		float* temp=dataset;
		for (i=0; i<pc.rows; i++)
		{
			const float* src = (float*)(&pc.data[i*pc.step]);
			float* dst = (float*)(&dataset[i*3]);

			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
		}
	}

 	flannIndex = flann_build_index(dataset, pc.rows, 3, &speedup, &p);

	if (!isCont)
		delete[] dataset;
	
	return (void*)flannIndex;
}

void destroy_flann(void* flannIndex)
{
	FLANNParameters p;
	p.log_level = FLANN_LOG_NONE;
	p.algorithm = FLANN_ALGORITHM_KD;
	p.checks = FLANN_CHECKS;
	p.trees = FLANN_NUM_TREES;
	p.branching = FLANN_NUM_BRANCHINGS;
	p.iterations = FLANN_NUM_ITERATIONS;
	p.target_precision = FLANN_TARGET_PRECISION;
	p.leaf_max_size = FLANN_MAX_LEAF_SIZE;
	p.eps = 0;
	
	flann_free_index(flannIndex, &p);
}

// For speed purposes this function assumes that PC, Indices and Distances are created with continuous structures
void query_pc_flann(void* flannIndex, Mat PC, Mat& Indices, Mat& Distances)
{	
	FLANNParameters p;
	float speedup=0;
	int i;
	
	const int numNeighbors = Indices.cols;

	p.log_level = FLANN_LOG_NONE;
	//p.log_destination = NULL;
	p.algorithm = FLANN_ALGORITHM_KD;
	p.checks = FLANN_CHECKS;
	p.trees = FLANN_NUM_TREES;
	p.branching = FLANN_NUM_BRANCHINGS;
	p.iterations = FLANN_NUM_ITERATIONS;
	p.target_precision = FLANN_TARGET_PRECISION;
	p.leaf_max_size = FLANN_MAX_LEAF_SIZE;
	p.cores=1;
	p.eps = 0;
	#if defined (T_OPENMP)
	p.cores=8;
	omp_set_num_threads(8);
#endif

	float* dataset;
	if (PC.isContinuous() && PC.rows==3)
	{
		dataset = (float*)PC.data;
	}
	else
	{
		dataset = new float[PC.rows*3];
		for (i=0; i<PC.rows; i++)
		{
			const float* src = (float*)(&PC.data[i*PC.step]);
			float* dst = (float*)(&dataset[i*3]);
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
		}
	}
/*
	float* dataset = new float[PC.rows*3];
	float* distances = new float[PC.rows*numNeighbors];
	int* indices = new int[PC.rows*numNeighbors];
	for (i=0; i<PC.rows; i++)
	{
		const float* src = (float*)(&PC.data[i*PC.step]);
		float* dst = (float*)(&dataset[i*3]);
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
	}*/

	flann_find_nearest_neighbors_index_float(flannIndex, dataset, PC.rows, (int*)Indices.data, (float*)Distances.data, numNeighbors, &p);

	// copy to opencv matrices
	/*for (i=0; i<PC.rows; i++)
	{
		const int *indicesRow = &indices[i*numNeighbors];
		const float *distancesRow = &distances[i*numNeighbors];
		int *IndicesRow = (float*)(&Indices.data[i*Indices.step]);
		float *DistancesRow = (float*)&(Distances.data[i*Distances.step]);
		for (j=0; j<numNeighbors; j++)
		{
			IndicesRow[j] = indicesRow[j];
			DistancesRow[j] = distancesRow[j];
		}
	}

	delete[] indices;
	delete[] distances;
	delete[] dataset;*/
	
	if (PC.isContinuous() && PC.rows==3)
		delete[] dataset;

}


// not yet complete
Mat sample_pc_kd_tree(Mat pc, float radius, int numNeighbors)
{
	void* flannIndex = index_pc_flann(pc);

	int interpNormals = (pc.cols==6);

	int* index = new int[numNeighbors];
	float* dist = new float[numNeighbors];

	for (int i=0; i<pc.rows; i++)
	{
		float* pcRow = (float*)(&pc.data[i*pc.step]);
		
		/*cvflann::Matrix<float> queryPt(pcRow, 1, 3);
		cvflann::Matrix<int> indexPt(index, 1,numNeighbors );
		cvflann::Matrix<float> distPt(dist, 1,numNeighbors );
 		flannIndex->radiusSearch(queryPt, indexPt, distPt, radius, searchParams);*/

		//int* rowInd = indices[i*indices.stride];
		//for (int j=0; j<indices.cols; j++)
		{


			/*if (!interpNormals)
			{
				int numPts = numNeighbors;
				// average the points
				for (j=0; j<numPts; j++)
				{
					float* pcRow = (float*)(&pc.data[rowInd[j]*pc.step]);
					
					px += (double)pcRow[0];
					py += (double)pcRow[1];
					pz += (double)pcRow[2];
				}

				px/=(double)numPts;
				py/=(double)numPts;
				pz/=(double)numPts;

				pcData[0]=(float)px;
				pcData[1]=(float)py;
				pcData[2]=(float)pz;
			}
			else
			{
				for (j=0; j<results.size(); j++)
				{
					px += (double)results[j][0];
					py += (double)results[j][1];
					pz += (double)results[j][2];
					nx += (double)results[j][3];
					ny += (double)results[j][4];
					nz += (double)results[j][5];
				}

				px/=(double)results.size();
				py/=(double)results.size();
				pz/=(double)results.size();
				nx/=(double)results.size();
				ny/=(double)results.size();
				nz/=(double)results.size();

				pcData[0]=(float)px;
				pcData[1]=(float)py;
				pcData[2]=(float)pz;

				// normalize the normals
				double norm = sqrt(nx*nx+ny*ny+nz*nz);

				if (norm>EPS)
				{
					pcData[3]=(float)(nx/norm);
					pcData[4]=(float)(ny/norm);
					pcData[5]=(float)(nz/norm);
				}
			}*/
		}
	}

	return Mat();	
}

// uses a volume instead of an octree
// TODO: Right now normals are required. 
// This is much faster than sample_pc_octree
Mat sample_pc_by_quantization(Mat pc, float xrange[2], float yrange[2], float zrange[2], float sampleStep)
{
	vector < vector<int> > map;

	int numSamplesDim = (int)(1.0/sampleStep);

	float xr = xrange[1] - xrange[0];
	float yr = yrange[1] - yrange[0];
	float zr = zrange[1] - zrange[0];

	int numPoints = 0;

	map.resize((numSamplesDim+1)*(numSamplesDim+1)*(numSamplesDim+1));

//#pragma omp parallel for
	for (int i=0; i<pc.rows; i++)
	{
		const float* point = (float*)(&pc.data[i * pc.step]);

		// quantize a point
		const int xCell =(int) ((float)numSamplesDim*(point[0]-xrange[0])/xr);
		const int yCell =(int) ((float)numSamplesDim*(point[1]-yrange[0])/yr);
		const int zCell =(int) ((float)numSamplesDim*(point[2]-zrange[0])/zr);
		const int index = xCell*numSamplesDim*numSamplesDim+yCell*numSamplesDim+zCell;

/*#pragma omp critical 
		{*/
			map[index].push_back(i);
//		}
	}

	for (int i=0; i<map.size(); i++)
	{
		numPoints += (map[i].size()>0);
	}

	Mat pcSampled = Mat(numPoints, pc.cols, CV_32F);
	int c = 0;

//#pragma omp parallel for
	for (int i=0; i<map.size(); i++)
	{
		double px=0, py=0, pz=0;
		double nx=0, ny=0, nz=0;

		vector<int> curCell = map[i];
		const int cn = curCell.size();
		if (cn>0)
		{
			for (int j=0; j<cn; j++)
			{
				const int ptInd = curCell[j];
				float* point = (float*)(&pc.data[ptInd * pc.step]);

				px += (double)point[0];
				py += (double)point[1];
				pz += (double)point[2];
				nx += (double)point[3];
				ny += (double)point[4];
				nz += (double)point[5];
			}						

			px/=(double)cn;
			py/=(double)cn;
			pz/=(double)cn;
			nx/=(double)cn;
			ny/=(double)cn;
			nz/=(double)cn;

			float *pcData = (float*)(&pcSampled.data[c*pcSampled.step[0]]);
			pcData[0]=(float)px;
			pcData[1]=(float)py;
			pcData[2]=(float)pz;

			// normalize the normals
			double norm = sqrt(nx*nx+ny*ny+nz*nz);

			if (norm>EPS)
			{
				pcData[3]=(float)(nx/norm);
				pcData[4]=(float)(ny/norm);
				pcData[5]=(float)(nz/norm);
			}
//#pragma omp atomic
			c++;

			curCell.clear();
		}
	}

	map.clear();
	return pcSampled;
}

// TODO : This uses a recursive octree. The recursion can sometimes get really deep 
// (if the points are too many). Instead I might use KD-Tree sampling but still
// we could find a more efficient way to do this.
Mat sample_pc_octree(Mat pc, float xrange[2], float yrange[2], float zrange[2], float sampleStep)
{
	TOctreeNode *oc = Mat2Octree(pc);

	// modified to sample within the cube

	float xstep = (xrange[1]-xrange[0]) * sampleStep;
	float ystep = (yrange[1]-yrange[0]) * sampleStep;
	float zstep = (zrange[1]-zrange[0]) * sampleStep;
	/*float xstep = sampleStep;
	float ystep = sampleStep;
	float zstep = sampleStep;*/

	float pdx = xrange[0], pdy=yrange[0], pdz=zrange[0];
	float dx=pdx+xstep, dy=pdy+ystep, dz=pdz+zstep;

	int numPoints = 0, c=0;
	int interpNormals = (pc.cols==6);		

	// count the number of points
	while (pdx<=xrange[1])
	{
		pdy=yrange[0]; 
		while (pdy<=yrange[1])
		{
			pdz=zrange[0];
			while (pdz<=zrange[1])
			{
				numPoints++;
				pdz+=zstep;
			}
			pdy+=ystep;
		}
		pdx+=xstep;
	}

	Mat pcSampled = Mat(numPoints, pc.cols, CV_32FC1);

	pdx = xrange[0]; 
	dx=pdx+xstep;  

	//while (dx<xrange[1] && dy<yrange[1] && dz<zrange[1])
	// query discrete bounding boxes over the octree
	while (pdx<=xrange[1])
	{
		float xbox[2] = {pdx, dx};
		pdy=yrange[0]; 
		dy=pdy+ystep;
		while (pdy<=yrange[1])
		{
			float ybox[2] = {pdy, dy};
			pdz=zrange[0];
			dz=pdz+zstep;
			while (pdz<=zrange[1])
			{
				float zbox[2] = {pdz, dz};
				int j;
				double px=0, py=0, pz=0;
				double nx=0, ny=0, nz=0;
				std::vector<float*> results;
				float *pcData = (float*)(&pcSampled.data[c*pcSampled.step[0]]);

				t_octree_query_in_bbox ( oc, xbox, ybox, zbox, results );

				if (results.size())
				{
					if (!interpNormals)
					{
						for (j=0; j<results.size(); j++)
						{
							px += (double)results[j][0];
							py += (double)results[j][1];
							pz += (double)results[j][2];
						}

						px/=(double)results.size();
						py/=(double)results.size();
						pz/=(double)results.size();

						pcData[0]=(float)px;
						pcData[1]=(float)py;
						pcData[2]=(float)pz;
					}
					else
					{
						for (j=0; j<results.size(); j++)
						{
							px += (double)results[j][0];
							py += (double)results[j][1];
							pz += (double)results[j][2];
							nx += (double)results[j][3];
							ny += (double)results[j][4];
							nz += (double)results[j][5];
						}

						px/=(double)results.size();
						py/=(double)results.size();
						pz/=(double)results.size();
						nx/=(double)results.size();
						ny/=(double)results.size();
						nz/=(double)results.size();

						pcData[0]=(float)px;
						pcData[1]=(float)py;
						pcData[2]=(float)pz;

						// normalize the normals
						double norm = sqrt(nx*nx+ny*ny+nz*nz);
						
						if (norm>EPS)
						{
							pcData[3]=(float)(nx/norm);
							pcData[4]=(float)(ny/norm);
							pcData[5]=(float)(nz/norm);
						}
					}

					c++;
				}

				results.clear();

				pdz=dz;
				dz+=zstep;
			}
			pdy=dy;
			dy+=ystep; 
		}
		pdx=dx;
		dx+=xstep; 
	}


	t_octree_destroy(oc);

	pcSampled=pcSampled.rowRange(0, c);

	return pcSampled;
}

void shuffle(int *array, size_t n)
{
	size_t i;
	for (i = 0; i < n - 1; i++) 
	{
		size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
		int t = array[j];
		array[j] = array[i];
		array[i] = t;
	}
}

Mat sample_pc_random(Mat PC, int numPoints)
{
	Mat sampledPC = Mat(numPoints, PC.cols, PC.type());

	// This is a slight inefficient way of doing it. Will optimize in final version.
	srand(time(0));
	int* randInd = new int[PC.rows];
	for (int i=0; i<PC.rows; i++)
		randInd[i]=i;
	shuffle(randInd, PC.rows);

	for (int i=0; i<numPoints; i++)
	{
		PC.row(randInd[i]).copyTo(sampledPC.row(i));
	}

	delete[] (randInd);

	return sampledPC;
}

// compute the standard bounding box
void compute_bbox_std(Mat pc, float xRange[2], float yRange[2], float zRange[2])
{
	Mat pcPts = pc.colRange(0, 3);
	int num = pcPts.rows;

	float* points = (float*)pcPts.data;
    GPointPair   pair;

    //printf( "Axis parallel bounding box\n" );
    GBBox   bbx;
    bbx.init();
    for  ( int  ind = 0; ind < num; ind++ )
        bbx.bound( (float*)(pcPts.data + (ind * pcPts.step)) );
    //bbx.dump();

	xRange[0]=bbx.min_coord(0);
	yRange[0]=bbx.min_coord(1);
	zRange[0]=bbx.min_coord(2);

	xRange[1]=bbx.max_coord(0);
	yRange[1]=bbx.max_coord(1);
	zRange[1]=bbx.max_coord(2);
}

// compute the oriented bounding box
double compute_diameter(Mat pc)
{
	Mat pcPts = pc.colRange(0, 3);
	int num = pcPts.rows;

	float* points = (float*)pcPts.data;
	gdiam_real  * gPoints;

	gPoints = (gdiam_point)malloc( sizeof( gdiam_point_t ) * num );

	// Initialize points in vector
	for  ( int  ind = 0; ind < num; ind++ ) 
	{
		const float* row = (float*)(pcPts.data + (ind * pcPts.step));
		gPoints[ ind * 3 + 0 ] = row[0];
		gPoints[ ind * 3 + 1 ] = row[1];
		gPoints[ ind * 3 + 2 ] = row[2];
	}

	GPointPair   pair;

	pair = gdiam_approx_diam_pair( (gdiam_real *)points, num, 0.0 );

	free(gPoints);
	
	return pair.distance;
}
/*
void compute_obb(Mat pc, float xRange[2], float yRange[2], float zRange[2])
{
	Mat pcPts = pc.colRange(0, 3);
	int num = pcPts.rows;

	float* points = (float*)pcPts.data;
	gdiam_real  * gPoints;

	gPoints = (gdiam_point)malloc( sizeof( gdiam_point_t ) * num );

	// Initialize points in vector
	for  ( int  ind = 0; ind < num; ind++ ) 
	{
		const float* row = (float*)(pcPts.data + (ind * pcPts.step));
		gPoints[ ind * 3 + 0 ] = row[0];
		gPoints[ ind * 3 + 1 ] = row[1];
		gPoints[ ind * 3 + 2 ] = row[2];
	}

	gdiam_point  * pnt_arr;
    gdiam_bbox   bb;

    pnt_arr = gdiam_convert( (gdiam_real *)points, num );
    bb = gdiam_approx_mvbb_grid_sample( pnt_arr, num, 5, 400 );

	free(pnt_arr);
	free(gPoints);
}*/



Mat normalize_pc(Mat pc, float scale)
{
	double minVal=0, maxVal=0;

	Mat x,y,z, pcn;
	pc.col(0).copyTo(x);
	pc.col(1).copyTo(y);
	pc.col(2).copyTo(z);

	float cx = cv::mean(x).val[0];
	float cy = cv::mean(y).val[0];
	float cz = cv::mean(z).val[0];

	cv::minMaxIdx(pc, &minVal, &maxVal);

	x=x-cx;
	y=y-cy;
	z=z-cz;
	pcn.create(pc.rows, 3, CV_32FC1);
	x.copyTo(pcn.col(0));
	y.copyTo(pcn.col(1));
	z.copyTo(pcn.col(2));

	cv::minMaxIdx(pcn, &minVal, &maxVal);
	pcn=(float)scale*(pcn)/((float)maxVal-(float)minVal);
	
	return pcn;
}

Mat normalize_pc_coeff(Mat pc, float scale, float* Cx, float* Cy, float* Cz, float* MinVal, float* MaxVal)
{
	double minVal=0, maxVal=0;

	Mat x,y,z, pcn;
	pc.col(0).copyTo(x);
	pc.col(1).copyTo(y);
	pc.col(2).copyTo(z);

	float cx = cv::mean(x).val[0];
	float cy = cv::mean(y).val[0];
	float cz = cv::mean(z).val[0];

	cv::minMaxIdx(pc, &minVal, &maxVal);

	x=x-cx;
	y=y-cy;
	z=z-cz;
	pcn.create(pc.rows, 3, CV_32FC1);
	x.copyTo(pcn.col(0));
	y.copyTo(pcn.col(1));
	z.copyTo(pcn.col(2));

	cv::minMaxIdx(pcn, &minVal, &maxVal);
	pcn=(float)scale*(pcn)/((float)maxVal-(float)minVal);

	*MinVal=minVal;
	*MaxVal=maxVal;
	*Cx=cx;
	*Cy=cy;
	*Cz=cz;
	
	return pcn;
}

Mat trans_pc_coeff(Mat pc, float scale, float Cx, float Cy, float Cz, float MinVal, float MaxVal)
{
	Mat x,y,z, pcn;
	pc.col(0).copyTo(x);
	pc.col(1).copyTo(y);
	pc.col(2).copyTo(z);

	x=x-Cx;
	y=y-Cy;
	z=z-Cz;
	pcn.create(pc.rows, 3, CV_32FC1);
	x.copyTo(pcn.col(0));
	y.copyTo(pcn.col(1));
	z.copyTo(pcn.col(2));

	pcn=(float)scale*(pcn)/((float)MaxVal-(float)MinVal);
	
	return pcn;
}

Mat transform_pc_pose(Mat pc, double Pose[16])
{
	Mat pct = Mat(pc.rows, pc.cols, CV_32F);

	double R[9], t[3];
	pose_to_rt(Pose, R, t); 

#if defined T_OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<pc.rows; i++)
	{
		const float *pcData = (float*)(&pc.data[i*pc.step]);
		float *pcDataT = (float*)(&pct.data[i*pct.step]);
		const float *n1 = &pcData[3];
		float *nT = &pcDataT[3];

		double p[4] = {(double)pcData[0], (double)pcData[1], (double)pcData[2], 1};
		double p2[4];
		
		matrix_product441(Pose, p, p2);

		// p2[3] should normally be 1
		if (fabs(p2[3])>EPS)
		{
			pcDataT[0] = (float)(p2[0]/p2[3]);
			pcDataT[1] = (float)(p2[1]/p2[3]);
			pcDataT[2] = (float)(p2[2]/p2[3]);
		}

		// Rotate the normals, too
		double n[3] = {(double)n1[0], (double)n1[1], (double)n1[2]}, n2[3];
		//matrix_product441(Pose, n, n2);
		
		matrix_product331(R, n, n2);
		double nNorm = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);

		if (nNorm>EPS)
		{
			nT[0]=(float)(n2[0]/nNorm);
			nT[1]=(float)(n2[1]/nNorm);
			nT[2]=(float)(n2[2]/nNorm);
		}
	}

	return pct;
}

Mat gen_random_mat(int rows, int cols, double mean, double stddev, int type)
{
	cv::Mat meanMat = mean*cv::Mat::ones(1,1,type);
	cv::Mat sigmaMat= stddev*cv::Mat::ones(1,1,type);
	cv::RNG rng(time(0));
	cv::Mat matr(rows, cols,type);
	rng.fill(matr, cv::RNG::NORMAL, meanMat, sigmaMat);

	return matr;
}

void get_rand_quat(double q[4])
{
	q[0] = (float)rand()/(float)(RAND_MAX);
	q[1] = (float)rand()/(float)(RAND_MAX);
	q[2] = (float)rand()/(float)(RAND_MAX);
	q[3] = (float)rand()/(float)(RAND_MAX);

	double n = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
	q[0]/=n;
	q[1]/=n;
	q[2]/=n;
	q[3]/=n;

	q[0]=fabs(q[0]);
}

void get_random_rotation(double R[9])
{
	double q[4];
	get_rand_quat(q);
	quaternion_to_matrix(q, R);
}

void get_random_pose(double Pose[16])
{
	double R[9], t[3];

	srand(time(0));
	get_random_rotation(R);

	t[0] = (float)rand()/(float)(RAND_MAX);
	t[1] = (float)rand()/(float)(RAND_MAX);
	t[2] = (float)rand()/(float)(RAND_MAX);

	rt_to_pose(R,t,Pose);
}

// this is not completely correct. Use get_random_pose instead
void generate_random_pose(double Pose[16], double scale)
{
	Mat S, U, V, R;
	Mat randR = gen_random_mat(3,3,scale,scale,CV_64FC1);

	cv::SVDecomp(randR, S, U, V);
	S = Mat::eye(3, 3, CV_64F);
	R = U*S*V;

	Mat randT = gen_random_mat(3,1,scale,scale,CV_64FC1);
	
	Pose[0]=R.at<double>(0,0);
	Pose[1]=R.at<double>(0,1);
	Pose[2]=R.at<double>(0,2);
	Pose[4]=R.at<double>(1,0);
	Pose[5]=R.at<double>(1,1);
	Pose[6]=R.at<double>(1,2);
	Pose[8]=R.at<double>(2,0);
	Pose[9]=R.at<double>(2,1);
	Pose[10]=R.at<double>(2,2);

	Pose[3]=randT.at<double>(0)-0.5;
	Pose[7]=randT.at<double>(1)-0.5;
	Pose[11]=randT.at<double>(2)-0.5;

	Pose[15] = 1;
}

Mat add_noise_pc(Mat pc, double scale)
{
	Mat randT = gen_random_mat(pc.rows,pc.cols,0,scale,CV_32FC1);
	return randT + pc;
}

//////////// COMPUTE NORMALS OF POINT CLOUD //////////////////////////////////
/*
	The routines below use the eigenvectors of the local covariance matrix
	to compute the normals of a point cloud.
	The algorithm uses FLANN and Joachim Kopp's fast 3x3 eigenvector computations
	to improve accuracy and increase speed
	Also, view point flipping as in point cloud library is implemented
*/

void mean_cov_local_pc(const float* pc, const int ws, const int point_count, double CovMat[3][3], double Mean[4])
{
	int i;
	double accu[16];	

	// For each point in the cloud
	for (i = 0; i < point_count; ++i)
	{
		const float* cloud = &pc[i*ws];
		accu [0] += cloud[0] * cloud[0];
		accu [1] += cloud[0] * cloud[1];
		accu [2] += cloud[0] * cloud[2];
		accu [3] += cloud[1] * cloud[1]; // 4
		accu [4] += cloud[1] * cloud[2]; // 5
		accu [5] += cloud[2] * cloud[2]; // 8
		accu [6] += cloud[0];
		accu [7] += cloud[1];
		accu [8] += cloud[2];
	}

	for (i = 0; i < 9; ++i)
		accu[i]/=(double)point_count;

	Mean[0] = accu[6]; 
	Mean[1] = accu[7]; 
	Mean[2] = accu[8];
	Mean[3] = 0;
	CovMat[0][0] = accu [0] - accu [6] * accu [6];
	CovMat[0][1] = accu [1] - accu [6] * accu [7];
	CovMat[0][2] = accu [2] - accu [6] * accu [8];
	CovMat[1][1] = accu [3] - accu [7] * accu [7];
	CovMat[1][2] = accu [4] - accu [7] * accu [8];
	CovMat[2][2] = accu [5] - accu [8] * accu [8];
	CovMat[1][0] = CovMat[0][1];
	CovMat[2][0] = CovMat[0][2];
	CovMat[2][1] = CovMat[1][2];

}

void mean_cov_local_pc_ind(const float* pc, const int* Indices, const int ws, const int point_count, double CovMat[3][3], double Mean[4])
{
	int i;
	double accu[16]={0};	

	for (i = 0; i < point_count; ++i)
	{
		const float* cloud = &pc[ Indices[i] * ws ];
		accu [0] += cloud[0] * cloud[0];
		accu [1] += cloud[0] * cloud[1];
		accu [2] += cloud[0] * cloud[2];
		accu [3] += cloud[1] * cloud[1]; // 4
		accu [4] += cloud[1] * cloud[2]; // 5
		accu [5] += cloud[2] * cloud[2]; // 8
		accu [6] += cloud[0];
		accu [7] += cloud[1];
		accu [8] += cloud[2];
	}

	for (i = 0; i < 9; ++i)
		accu[i]/=(double)point_count;

	Mean[0] = accu[6]; 
	Mean[1] = accu[7]; 
	Mean[2] = accu[8];
	Mean[3] = 0;
	CovMat[0][0] = accu [0] - accu [6] * accu [6];
	CovMat[0][1] = accu [1] - accu [6] * accu [7];
	CovMat[0][2] = accu [2] - accu [6] * accu [8];
	CovMat[1][1] = accu [3] - accu [7] * accu [7];
	CovMat[1][2] = accu [4] - accu [7] * accu [8];
	CovMat[2][2] = accu [5] - accu [8] * accu [8];
	CovMat[1][0] = CovMat[0][1];
	CovMat[2][0] = CovMat[0][2];
	CovMat[2][1] = CovMat[1][2];

}

void flip_normal_viewpoint(const float* point, double vp_x, double vp_y, double vp_z, double *nx, double *ny, double *nz)
{
	double cos_theta;

	// See if we need to flip any plane normals
	vp_x -= (double)point[0];
	vp_y -= (double)point[1];
	vp_z -= (double)point[2];

	// Dot product between the (viewpoint - point) and the plane normal
	cos_theta = (vp_x * (*nx) + vp_y * (*ny) + vp_z * (*nz));

	// Flip the plane normal
	if (cos_theta < 0)
	{
		(*nx) *= -1;
		(*ny) *= -1;
		(*nz) *= -1;
	}
}

int compute_normals_pc_3d(const Mat PC, Mat& PCNormals, const int NumNeighbors, const bool FlipViewpoint, const double viewpoint[3])
{
	int i, j, ws;
	
	if(PC.cols!=3 && PC.cols!=6) // 3d data is expected
	{
		return -1;
	}

	//size_t steps[2] = {1, 3};
	int sizes[2] = {PC.rows, 3};
	int sizesResult[2] = {PC.rows, NumNeighbors};
	float* dataset = new float[PC.rows*3];
	float* distances = new float[PC.rows*NumNeighbors];
	int* indices = new int[PC.rows*NumNeighbors];

	for (i=0; i<PC.rows; i++)
	{
		const float* src = (float*)(&PC.data[i*PC.step]);
		float* dst = (float*)(&dataset[i*3]);

		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
	}

	Mat PCInput(2, sizes, CV_32F, dataset, 0);

	void* flannIndex = index_pc_flann(PCInput);

	Mat Indices(2, sizesResult, CV_32S, indices, 0);
	Mat Distances(2, sizesResult, CV_32F, distances, 0);

	query_pc_flann(flannIndex, PCInput, Indices, Distances);
	destroy_flann(flannIndex); flannIndex = 0;

	PCNormals = Mat(PC.rows, 6, CV_32F);
	
	for (i=0; i<PC.rows; i++)
	{		
		double C[3][3], mu[4];
		double Q[3][3], w[3], curvature=0;
		double trace ;
		int k;
		const float* pci = &dataset[i*3];
		float* pcr = (float*)(&PCNormals.data[i*PCNormals.step]);
		int minEigID = 2;
		double nr[3];

		int* indLocal = &indices[i*NumNeighbors];

		// compute covariance matrix
		mean_cov_local_pc_ind(dataset, indLocal, 3, NumNeighbors, C, mu);

		// eigenvectors of covariance matrix
		dsyevh3(C, Q, w);

		// find min eigenvalue
		if (w[0]<w[1])
		{
			if (w[0]<w[2])
				minEigID = 0; // w0 is min
		}
		else
		{
			if (w[1]<w[2])
				minEigID = 1;
		}

		//printf("%f %f %f, %f\n", w[0], w[1], w[2], w[minEigID]);

		trace = C[0][0] + C[1][1] + C[2][2];
		
		if (trace>0)
		{
			curvature = fabs(w[minEigID] / trace);
		}

		pcr[0] = pci[0];
		pcr[1] = pci[1];
		pcr[2] = pci[2];

		nr[0] = Q[0][minEigID];
		nr[1] = Q[1][minEigID];
		nr[2] = Q[2][minEigID];

		if (FlipViewpoint)
		{
			flip_normal_viewpoint(pci, viewpoint[0], viewpoint[1], viewpoint[2], &nr[0], &nr[1], &nr[2]);
		}

		pcr[3] = nr[0];
		pcr[4] = nr[1];
		pcr[5] = nr[2];
	}

	delete[] indices;
	delete[] distances;
	delete[] dataset;
	
	return 1;
}

//////////////////////////////////////// END OF NORMAL COMPUTATIONS ///////////////////////////////////

////////////// POINT SELECTION FOR ICP //////////////////

//
//void prune_scene_for_icp(Mat Model, Mat Scene, Mat& ScenePrune, Mat& ModelPrune)
//{
//	int sizesResult[2] = {Model.rows, 1};
//	float* distances = new float[Model.rows];
//	int* indices = new int[Model.rows];
//
//	void* flannIndex = index_pc_flann(Scene);
//
//	Mat Indices(2, sizesResult, CV_32S, indices, 0);
//	Mat Distances(2, sizesResult, CV_32F, distances, 0);
//
//	query_pc_flann(flannIndex, Model, Indices, Distances);
//	destroy_flann(flannIndex); flannIndex = 0;
//
//	// use robust weighting for outlier treatment
//	int* indicesModel = new int[Model.rows];
//	int* indicesScene = new int[Model.rows];
//	int numInliers = 0;
//
//	float sigma = madsigma(distances, Model.rows);
//
//	for (int i=0; i<Model.rows;	i++)
//	{
//		if (distances[i] < 1.4826 * sigma)
//		{
//			indicesModel[numInliers] = i;
//			indicesScene[numInliers] = indices[i];
//			numInliers++;
//		}
//	}
//
//	ScenePrune = Mat(numInliers, Model.cols, CV_32F);
//	ModelPrune = Mat(numInliers, Model.cols, CV_32F);
//
//	delete[] indicesModel;	
//	delete[] indicesScene;	
//}