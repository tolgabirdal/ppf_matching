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

#ifndef _OPENCV_POSE3D_HPP_
#define _OPENCV_POSE3D_HPP_

#include <vector>
#include "c_utils.hpp"

namespace cv 
{
	namespace ppf_match_3d 
	{
		class CV_EXPORTS Pose3D
		{
		public:
			Pose3D()
			{
				alpha=0; 
				modelIndex=0; 
				numVotes=0;
				residual = 0;

				for (int i=0; i<16; i++)
					Pose[i]=0;
			};

			Pose3D(double Alpha, unsigned int ModelIndex=0, unsigned int NumVotes=0)
			{
				alpha = Alpha;
				modelIndex = ModelIndex;
				numVotes = NumVotes;
				residual=0;

				for (int i=0; i<16; i++)
					Pose[i]=0;
			};

			void updatePose(double NewPose[16]);
			void updatePose(double NewR[9], double NewT[3]);
			void updatePoseQuat(double Q[4], double NewT[3]);
			void appendPose(double IncrementalPose[16]);
			void printPose();

			Pose3D* clone();

			int writePose(FILE* f);
			int readPose(FILE* f);
			int writePose(const std::string& FileName);
			int readPose(const std::string& FileName);

			virtual ~Pose3D(){};

			double alpha, residual;
			unsigned int modelIndex;
			unsigned int numVotes;
			double Pose[16], angle, t[3], q[4];
		};

		class CV_EXPORTS PoseCluster3D
		{
		public:
			PoseCluster3D() 
			{
				//poseList.clear();
				numVotes=0;
				id=0;
			};

			PoseCluster3D(Pose3D* newPose) 
			{
				//poseList.clear();
				poseList.push_back(newPose);
				numVotes=newPose->numVotes;
				id=0;
			};

			PoseCluster3D(Pose3D* newPose, int id) 
			{
				//poseList.clear();
				poseList.push_back(newPose);
				this->numVotes = newPose->numVotes;
				this->id = id;
			};

			virtual ~PoseCluster3D()
			{
				numVotes=0;
				id=0;
				//poseList.clear();
			};

			void addPose(Pose3D* newPose) ;
			int writePoseCluster(FILE* f);
			int readPoseCluster(FILE* f);
			int writePoseCluster(const std::string& FileName);
			int readPoseCluster(const std::string& FileName);

			std::vector < Pose3D* > poseList;
			int numVotes;
			int id;
		};
	}
}

#endif
