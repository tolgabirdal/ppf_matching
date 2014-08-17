
#ifndef __OPENCV_ICP_HPP__
#define __OPENCV_ICP_HPP__

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION > 2
#include <opencv2/core/utility.hpp>
#else
#include <opencv2/core/core.hpp>
#endif

#include "Pose3D.hpp"
#include "c_utils.hpp"
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

			/*!
			Default constructor
    `		*/
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

			/**
			 *  \brief Brief
			 *  \param [in] tolerence Tolerence parameter controls the accuracy of registration at each iteration of ICP.
			 *  \param [in] rejectionScale Robust outlier rejection is applied for robustness. This value actually corresponds to the standard deviation coefficient. Points with rejectionScale * \sigma are ignored during registration.
			 *  \param [in] numLevels Number of pyramid levels to proceed. Deep pyramids increase speed but decrease accuracy. Too coarse pyramids might have computational overhead on top of the inaccurate registrtaion. This parameter should be chosen to optimize a balance. Typical values range from 4 to 10.
			 *  \param [in] sampleType Currently this parameter is ignored and only uniform sampling is applied. Leave it as 0.
			 *  \param [in] numMaxCorr Currently this parameter is ignored and only PickyICP is applied. Leave it as 1.
			 *  \return
			 *  
			 *  \details Constructor
			 */
			ICP(const int iterations, const float tolerence=0.05, const float rejectionScale=2.5, const int numLevels=6, const int sampleType=0, const int numMaxCorr=1) 
			{
				Tolerence = tolerence;
				NumNeighborsCorr = numMaxCorr;
				RejectionScale = rejectionScale;
				MaxIterations = iterations;
				NumLevels = numLevels;
				SampleType = sampleType;
			};

			/**
			 *  \brief Brief
			 *  
			 *  \return Return_Description
			 *  
			 *  \details Details
			 */void operator()(InputArray SrcPC, InputArray DstPC, double& Residual, double Pose[16]) const;
			 
			/**
			 *  \brief Brief
			 *  
			 *  \param [in] SrcPC The input point cloud for the model. Expected to have the normals (Nx6). Currently,
			 *  CV_32F is the only supported data type.
			 *  \param [in] DstPC The input point cloud for the scene. It is assumed that the model is registered on the scene. Scene remains static. Expected to have the normals (Nx6). Currently, CV_32F is the only supported data type.
			 *  \param [out] Residual The output registration error. 
			 *  \return On successful termination, the function returns 0.
			 *  
			 *  \details It is assumed that the model is registered on the scene. Scene remains static, while the model transforms. The output poses transform the models onto the scene. Because of the point to plane minimization, the scene is expected to have the normals available. Expected to have the normals (Nx6).
			 */ 
			int registerModelToScene(const Mat& SrcPC, const Mat& DstPC, double& Residual, double Pose[16]);
			
			/**
			 *  \brief Brief
			 *  
			 *  \param [in] SrcPC The input point cloud for the model. Expected to have the normals (Nx6). Currently,
			 *  CV_32F is the only supported data type.
			 *  \param [in] DstPC The input point cloud for the scene. Currently, CV_32F is the only supported data type.
			 *  \param [out] Poses List output of poses. For more detailed information check out Pose3D.
			 *  \return On successful termination, the function returns 0.
			 *  
			 *  \details It is assumed that the model is registered on the scene. Scene remains static, while the model transforms. The output poses transform the models onto the scene. Because of the point to plane minimization, the scene is expected to have the normals available. Expected to have the normals (Nx6).
			 */ 
			int registerModelToScene(const Mat& SrcPC, const Mat& DstPC, std::vector<Pose3D*>& Poses);

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