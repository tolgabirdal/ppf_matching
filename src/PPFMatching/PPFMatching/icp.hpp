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

#ifndef __OPENCV_ICP_HPP__
#define __OPENCV_ICP_HPP__

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION > 2
#include <opencv2/core/utility.hpp>
#else
#include <opencv2/core/core.hpp>
#endif

#include "pose_3d.hpp"
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
				m_tolerence = 0.05;
				m_rejectionScale = 2.5;
				m_maxItereations = 250;
				m_numLevels = 6;
				m_sampleType = 0;
				m_numNeighborsCorr = 1;
			}

			virtual ~ICP() { }

			/**
			 *  \brief Brief
			 *  \param [in] tolerence tolerence parameter controls the accuracy of registration at each iteration of ICP.
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
				m_tolerence = tolerence;
				m_numNeighborsCorr = numMaxCorr;
				m_rejectionScale = rejectionScale;
				m_maxItereations = iterations;
				m_numLevels = numLevels;
				m_sampleType = sampleType;
			};

			/**
			 *  \brief Perform registration
			 *  
			 *  \param [in] srcPC The input point cloud for the model. Expected to have the normals (Nx6). Currently,
			 *  CV_32F is the only supported data type.
			 *  \param [in] dstPC The input point cloud for the scene. It is assumed that the model is registered on the scene. Scene remains static. Expected to have the normals (Nx6). Currently, CV_32F is the only supported data type.
			 *  \param [out] Residual The output registration error. 
			 *  \return On successful termination, the function returns 0.
			 *  
			 *  \details It is assumed that the model is registered on the scene. Scene remains static, while the model transforms. The output poses transform the models onto the scene. Because of the point to plane minimization, the scene is expected to have the normals available. Expected to have the normals (Nx6).
			 */ 
			int registerModelToScene(const Mat& srcPC, const Mat& dstPC, double& Residual, double Pose[16]);
			
			/**
			 *  \brief Perform registration with multiple initial poses
			 *  
			 *  \param [in] srcPC The input point cloud for the model. Expected to have the normals (Nx6). Currently,
			 *  CV_32F is the only supported data type.
			 *  \param [in] dstPC The input point cloud for the scene. Currently, CV_32F is the only supported data type.
			 *  \param [out] Poses List output of poses. For more detailed information check out Pose3D.
			 *  \return On successful termination, the function returns 0.
			 *  
			 *  \details It is assumed that the model is registered on the scene. Scene remains static, while the model transforms. The output poses transform the models onto the scene. Because of the point to plane minimization, the scene is expected to have the normals available. Expected to have the normals (Nx6).
			 */ 
			int registerModelToScene(const Mat& srcPC, const Mat& dstPC, std::vector<Pose3D*>& Poses);

		private:
			float m_tolerence;
			int m_maxItereations;
			float m_rejectionScale;
			int m_numNeighborsCorr;
			int m_numLevels;
			int m_sampleType;

		};
	}
}

#endif