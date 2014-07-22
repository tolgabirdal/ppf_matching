
#ifndef __T_ICP_H_
#define __T_ICP_H_

#include <opencv2/core/core.hpp>
//#include <opencv2/flann/flann.hpp>
#include "t_octree.h"

//using namespace cv;

typedef int (*T_ICP_CALLBACK)(double Pose[4][4], void* UserData);

#if defined (__cplusplus)
extern "C" {
#endif

	int t_icp_register(const cv::Mat SrcPC, const cv::Mat DstPC, const float Tolerence, const int MaxIterations, const float RejectionScale, const int NumNeighborsCorr, const int NumLevels, const int SampleType, const T_ICP_CALLBACK RegistrationVisFunction, float* Residual, double Pose[16]);

#if defined (__cplusplus)
}
#endif

#endif
