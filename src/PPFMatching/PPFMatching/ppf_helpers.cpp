//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                          License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2014, OpenCV Foundation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
// Author: Tolga Birdal

#include "precomp.hpp"
#include "ppf_helpers.hpp"
#include "c_utils.hpp"
#include <time.h>
#include <fstream>
#include <vector>
#include <iostream>
#include "flann/flann.h"

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

namespace cv
{
namespace ppf_match_3d
{

Mat loadPLYSimple(const char* fileName, int numVertices, int withNormals)
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
    {
        if ( str.substr(0, 14) == "element vertex" )
        {
        
        }
        getline(ifs, str);
    }
    
    float dummy =  0;
    for (int i = 0; i < numVertices; i++)
    {
        float* data = (float*)(&cloud.data[i*cloud.step[0]]);
        if (withNormals)
        {
            ifs >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
            
            // normalize to unit norm
            double norm = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);
            if (norm>0.00001)
            {
                data[3]/=(float)norm;
                data[4]/=(float)norm;
                data[5]/=(float)norm;
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


void writePLY(Mat PC, const char* FileName)
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
    if (vertNum==6)
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

Mat samplePCUniform(Mat PC, int sampleStep)
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

Mat samplePCUniformInd(Mat PC, int sampleStep, vector<int> &indices)
{
    int numRows = cvRound((double)PC.rows/(double)sampleStep);
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

void* indexPCFlann(Mat pc)
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

void destroyFlann(void* flannIndex)
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
void queryPCFlann(void* flannIndex, Mat PC, Mat& Indices, Mat& Distances)
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
    
    flann_find_nearest_neighbors_index_float(flannIndex, dataset, PC.rows, (int*)Indices.data, (float*)Distances.data, numNeighbors, &p);
    
    if (PC.isContinuous() && PC.rows==3)
        delete[] dataset;
}

// uses a volume instead of an octree
// TODO: Right now normals are required.
// This is much faster than sample_pc_octree
Mat samplePCByQuantization(Mat pc, float xrange[2], float yrange[2], float zrange[2], float sampleStep, int weightByCenter)
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
        //  }
    }
    
    for (unsigned int i=0; i<map.size(); i++)
    {
        numPoints += (map[i].size()>0);
    }
    
    Mat pcSampled = Mat(numPoints, pc.cols, CV_32F);
    int c = 0;
    
    for (unsigned int i=0; i<map.size(); i++)
    {
        double px=0, py=0, pz=0;
        double nx=0, ny=0, nz=0;
        
        vector<int> curCell = map[i];
        const int cn = curCell.size();
        if (cn>0)
        {
            if (weightByCenter)
            {
                int xCell, yCell, zCell;
                double xc, yc, zc;
                double weightSum = 0 ;
                zCell = i % numSamplesDim;
                yCell = ((i-zCell)/numSamplesDim) % numSamplesDim;
                xCell = ((i-zCell-yCell*numSamplesDim)/(numSamplesDim*numSamplesDim));
                
                xc = ((double)xCell+0.5) * (double)xr/numSamplesDim + (double)xrange[0];
                yc = ((double)yCell+0.5) * (double)yr/numSamplesDim + (double)yrange[0];
                zc = ((double)zCell+0.5) * (double)zr/numSamplesDim + (double)zrange[0];
                
                for (int j=0; j<cn; j++)
                {
                    const int ptInd = curCell[j];
                    float* point = (float*)(&pc.data[ptInd * pc.step]);
                    const double dx = point[0]-xc;
                    const double dy = point[1]-yc;
                    const double dz = point[2]-zc;
                    const double d = sqrt(dx*dx+dy*dy+dz*dz);
                    double w = 0;
                    // exp( - (distance/h)**2 )
                    //const double w = exp(-d*d);
                    
                    if (d>EPS)
                        w = 1.0/d;
                        
                    //float weights[3]={1,1,1};
                    px += w*(double)point[0];
                    py += w*(double)point[1];
                    pz += w*(double)point[2];
                    nx += w*(double)point[3];
                    ny += w*(double)point[4];
                    nz += w*(double)point[5];
                    
                    weightSum+=w;
                }
                px/=(double)weightSum;
                py/=(double)weightSum;
                pz/=(double)weightSum;
                nx/=(double)weightSum;
                ny/=(double)weightSum;
                nz/=(double)weightSum;
            }
            else
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
                
            }
            
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

// compute the standard bounding box
void computeBboxStd(Mat pc, float xRange[2], float yRange[2], float zRange[2])
{
    Mat pcPts = pc.colRange(0, 3);
    int num = pcPts.rows;
    
    float* points = (float*)pcPts.data;
    
    xRange[0] = points[0];
    xRange[1] = points[0];
    yRange[0] = points[1];
    yRange[1] = points[1];
    zRange[0] = points[2];
    zRange[1] = points[2];
    
    for  ( int  ind = 0; ind < num; ind++ )
    {
        const float* row = (float*)(pcPts.data + (ind * pcPts.step));
        const float x = row[0];
        const float y = row[1];
        const float z = row[2];
        
        if (x<xRange[0])
            xRange[0]=x;
        if (x>xRange[1])
            xRange[1]=x;
            
        if (y<yRange[0])
            yRange[0]=y;
        if (y>yRange[1])
            yRange[1]=y;
            
        if (z<zRange[0])
            zRange[0]=z;
        if (z>zRange[1])
            zRange[1]=z;
    }
}

Mat normalizePCCoeff(Mat pc, float scale, float* Cx, float* Cy, float* Cz, float* MinVal, float* MaxVal)
{
    double minVal=0, maxVal=0;
    
    Mat x,y,z, pcn;
    pc.col(0).copyTo(x);
    pc.col(1).copyTo(y);
    pc.col(2).copyTo(z);
    
    float cx = (float) cv::mean(x).val[0];
    float cy = (float) cv::mean(y).val[0];
    float cz = (float) cv::mean(z).val[0];
    
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
    
    *MinVal=(float)minVal;
    *MaxVal=(float)maxVal;
    *Cx=(float)cx;
    *Cy=(float)cy;
    *Cz=(float)cz;
    
    return pcn;
}

Mat transPCCoeff(Mat pc, float scale, float Cx, float Cy, float Cz, float MinVal, float MaxVal)
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

Mat transformPCPose(Mat pc, double Pose[16])
{
    Mat pct = Mat(pc.rows, pc.cols, CV_32F);
    
    double R[9], t[3];
    poseToRT(Pose, R, t);
    
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
        
        matrixProduct441(Pose, p, p2);
        
        // p2[3] should normally be 1
        if (fabs(p2[3])>EPS)
        {
            pcDataT[0] = (float)(p2[0]/p2[3]);
            pcDataT[1] = (float)(p2[1]/p2[3]);
            pcDataT[2] = (float)(p2[2]/p2[3]);
        }
        
        // Rotate the normals, too
        double n[3] = {(double)n1[0], (double)n1[1], (double)n1[2]}, n2[3];
        //matrixProduct441(Pose, n, n2);
        
        matrixProduct331(R, n, n2);
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

Mat genRandomMat(int rows, int cols, double mean, double stddev, int type)
{
    cv::Mat meanMat = mean*cv::Mat::ones(1,1,type);
    cv::Mat sigmaMat= stddev*cv::Mat::ones(1,1,type);
    cv::RNG rng(time(0));
    cv::Mat matr(rows, cols,type);
    rng.fill(matr, cv::RNG::NORMAL, meanMat, sigmaMat);
    
    return matr;
}

void getRandQuat(double q[4])
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

void getRandomRotation(double R[9])
{
    double q[4];
    getRandQuat(q);
    quatToDCM(q, R);
}

void getRandomPose(double Pose[16])
{
    double R[9], t[3];
    
    srand((unsigned int)time(0));
    getRandomRotation(R);
    
    t[0] = (float)rand()/(float)(RAND_MAX);
    t[1] = (float)rand()/(float)(RAND_MAX);
    t[2] = (float)rand()/(float)(RAND_MAX);
    
    rtToPose(R,t,Pose);
}

Mat addNoisePC(Mat pc, double scale)
{
    Mat randT = genRandomMat(pc.rows,pc.cols,0,scale,CV_32FC1);
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

void meanCovLocalPC(const float* pc, const int ws, const int point_count, double CovMat[3][3], double Mean[4])
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

void meanCovLocalPCInd(const float* pc, const int* Indices, const int ws, const int point_count, double CovMat[3][3], double Mean[4])
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


CV_EXPORTS int computeNormalsPC3d(const Mat PC, Mat& PCNormals, const int NumNeighbors, const bool FlipViewpoint, const double viewpoint[3])
{
    int i;
    
    if (PC.cols!=3 && PC.cols!=6) // 3d data is expected
    {
        //return -1;
        CV_Error(cv::Error::BadImageSize, "PC should have 3 or 6 elements in its columns");
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
    
    void* flannIndex = indexPCFlann(PCInput);
    
    Mat Indices(2, sizesResult, CV_32S, indices, 0);
    Mat Distances(2, sizesResult, CV_32F, distances, 0);
    
    queryPCFlann(flannIndex, PCInput, Indices, Distances);
    destroyFlann(flannIndex);
    flannIndex = 0;
    
    PCNormals = Mat(PC.rows, 6, CV_32F);
    
    for (i=0; i<PC.rows; i++)
    {
        double C[3][3], mu[4];
        double w[3], curvature=0;
        double trace;
        const float* pci = &dataset[i*3];
        float* pcr = (float*)(&PCNormals.data[i*PCNormals.step]);
        int minEigID = 2;
        double nr[3];
        
        int* indLocal = &indices[i*NumNeighbors];
        
        // compute covariance matrix
        meanCovLocalPCInd(dataset, indLocal, 3, NumNeighbors, C, mu);
        
        // eigenvectors of covariance matrix
        eigenLowest33(C, nr);
        
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
        
        trace = C[0][0] + C[1][1] + C[2][2];
        
        if (trace>0)
        {
            curvature = fabs(w[minEigID] / trace);
        }
        
        pcr[0] = pci[0];
        pcr[1] = pci[1];
        pcr[2] = pci[2];
        
        if (FlipViewpoint)
        {
            flipNormalViewpoint(pci, viewpoint[0], viewpoint[1], viewpoint[2], &nr[0], &nr[1], &nr[2]);
        }
        
        pcr[3] = (float)nr[0];
        pcr[4] = (float)nr[1];
        pcr[5] = (float)nr[2];
    }
    
    delete[] indices;
    delete[] distances;
    delete[] dataset;
    
    return 1;
}

//////////////////////////////////////// END OF NORMAL COMPUTATIONS ///////////////////////////////////
}

}