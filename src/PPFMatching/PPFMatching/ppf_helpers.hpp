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

#ifndef __OPENCV_PPF_3D_HELPERS_HPP__
#define __OPENCV_PPF_3D_HELPERS_HPP__

#include <opencv2/core/core.hpp>

namespace cv
{
namespace ppf_match_3d
{

cv::Mat loadPLYSimple(const char* fileName, int numVertices, int withNormals);
void writePLY(cv::Mat PC, const char* FileName);

cv::Mat samplePCUniform(cv::Mat PC, int sampleStep);
cv::Mat samplePCUniformInd(cv::Mat PC, int sampleStep, std::vector<int>& indices);
cv::Mat samplePCByQuantization(cv::Mat pc, float xrange[2], float yrange[2], float zrange[2], float sample_step_relative, int weightByCenter=0);

void computeBboxStd(cv::Mat pc, float xRange[2], float yRange[2], float zRange[2]);

void* indexPCFlann(cv::Mat pc);
void destroyFlann(void* flannIndex);
void queryPCFlann(void* flannIndex, cv::Mat PC, cv::Mat& Indices, cv::Mat& Distances);

cv::Mat normalize_pc(cv::Mat pc, float scale);
cv::Mat normalizePCCoeff(cv::Mat pc, float scale, float* Cx, float* Cy, float* Cz, float* MinVal, float* MaxVal);
cv::Mat transPCCoeff(cv::Mat pc, float scale, float Cx, float Cy, float Cz, float MinVal, float MaxVal);
cv::Mat transformPCPose(cv::Mat pc, double Pose[16]);

void getRandomPose(double Pose[16]);
cv::Mat addNoisePC(cv::Mat pc, double scale);

/**
 *  \brief Compute the normals of an arbitrary point cloud
 *
 *  \param [in] PC Input point cloud to compute the normals for.
 *  \param [in] PCNormals Output point cloud
 *  \param [in] NumNeighbors Number of neighbors to take into account in a local region
 *  \param [in] FlipViewpoint Should normals be flipped to a viewing direction?
 *  \return Returns 0 on success
 *
 *  \details computeNormalsPC3d uses a plane fitting approach to smoothly compute
 *  local normals. Normals are obtained through the eigenvector of the covariance
 *  matrix, corresponding to the smallest eigen value.
 *  If PCNormals is provided to be an Nx6 matrix, then no new allocation
 *  is made, instead the existing memory is overwritten.
 */
int computeNormalsPC3d(const cv::Mat PC, cv::Mat& PCNormals, const int NumNeighbors, const bool FlipViewpoint, const double viewpoint[3]);
}
}

#endif