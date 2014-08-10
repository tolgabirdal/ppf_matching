
#ifndef __OPENCV_ICP_HPP__
#define __OPENCV_ICP_HPP__

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION > 2
#include <opencv2/core/utility.hpp>
#else
#include <opencv2/core/core.hpp>
#endif

#include "Pose3D.hpp"
#include "c_utils.h"
#include <vector>

using namespace std;

namespace cv 
{
	namespace ppf_match_3d 
	{
		class CV_EXPORTS ICP
		{
		public:

			enum ICP_SAMPLING_TYPE
			{
				ICP_SAMPLING_TYPE_UNIFORM, ICP_SAMPLING_TYPE_GELFAND
			};

			ICP() 
			{
				Tolerence = 0.05;
				RejectionScale = 2.5;
				MaxIterations = 250;
				NumLevels = 6;
				SampleType = 0;
				NumNeighborsCorr = 1;
			}

			~ICP() { }

			ICP(const int iterations, const float tolerence=0.05, const float rejectionScale=2.5, const int numLevels=6, const int sampleType=0, const int numMaxCorr=1) 
			{
				Tolerence = tolerence;
				NumNeighborsCorr = numMaxCorr;
				RejectionScale = rejectionScale;
				MaxIterations = iterations;
				NumLevels = numLevels;
				SampleType = sampleType;
			};

			void operator()(InputArray SrcPC, InputArray DstPC, double& Residual, double Pose[16]) const;
			int registerModelToScene(const cv::Mat& SrcPC, const cv::Mat& DstPC, double& Residual, double Pose[16]);
			int registerModelToScene(const cv::Mat& SrcPC, const cv::Mat& DstPC, std::vector<Pose3D*>& Poses);

		private:
			float Tolerence;
			int MaxIterations;
			float RejectionScale;
			int NumNeighborsCorr;
			int NumLevels;
			int SampleType;


			

		};
	}
}

#endif