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
/**
** ppf_match_3d : Interfaces for matching 3d surfaces in 3d scenes. This module implements the algorithm from Bertram Drost and Slobodan Ilic.
** Use: Read a 3D model, load a 3D scene and match the model to the scene
**
**
**  Creation - 2014
**      Author: Tolga Birdal (tbirdal@gmail.com)
**
** Refer to the following research paper for more information:
**  B. Drost, Markus Ulrich, N. Navab, S. Ilic
Model Globally, Match Locally: Efficient and Robust 3D Object Recognition
IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR), San Francisco, California (USA), June 2010.
***/
// Author: Tolga Birdal


#ifndef __OPENCV_OBJDETECT_PPF_3D_HPP__
#define __OPENCV_OBJDETECT_PPF_3D_HPP__

#include "t_hash_int.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION > 2
#include <opencv2/core/utility.hpp>
#else
#include <opencv2/core/core.hpp>
#endif

#include "ppf_helpers.hpp"
#include "c_utils.hpp"
#include "pose_3d.hpp"
#include <vector>

using namespace std;

namespace cv
{
namespace ppf_match_3d
{

#define T_PPF_LENGTH 5

/**
* @struct THash
* @brief Struct, holding a node in the hashtable
*/
typedef struct THash
{
    int id;
    int i, ppfInd;
} THash;

/**
* @class PPF3DDetector
* @brief Class, allowing the load and matching 3D models.
* Typical Use:
*
* // Train a model
*   ppf_match_3d::PPF3DDetector detector(0.05, 0.05);
* detector.trainModel(pc);
* // Search the model in a given scene
*   vector < Pose3D* > results;
* detector.match(pcTest, results, 1.0/5.0,0.05);
*
*/
class CV_EXPORTS PPF3DDetector
{
    public:
    
        /**
            * \brief Empty constructor. Sets default arguments
            */
        PPF3DDetector();
        
        /**
        * Constructor with arguments
        * @param [in] RelativeSamplingStep Set the sampling distance for the pose refinement relative to the object's diameter. Decreasing this value leads to a more accurate pose refinement but a larger model and slower model generation and refinement. Increasing the value leads to a less accurate pose refinement but a smaller model and faster model generation and matching. Beware of the memory consumption when using large values.
        * @param [in] RelativeDistanceStep Set the discretization distance of the point pair distance relative to the model's diameter. This value should default to the value of RelativeSamplingStep. For noisy scenes, the value can be increased to improve the robustness of the matching against noisy points.
        * @param [in] Set the discretization of the point pair orientation as the number of subdivisions of the angle. Increasing the value increases the precision of the matching but decreases the robustness against incorrect normal directions. Decreasing the value decreases the precision of the matching but increases the robustness against incorrect normal directions. For very noisy scenes where the normal directions can not be computed accurately, the value can be set to 25 or 20.
        */
        PPF3DDetector(const double RelativeSamplingStep, const double RelativeDistanceStep=0.05, const double NumAngles=30);
        
        virtual ~PPF3DDetector();
        
        /**
        *  Set the parameters for the search    *
        *  @param [in] numPoses The maximum number of poses to return
        *  @param [in] positionThreshold Position threshold controlling the similarity of translations. Depends on the units of calibration.
        *  @param [in] rotationThreshold Position threshold controlling the similarity of rotations. This parameter can be perceived as a threshold over the difference of angles
        *  @param [in] minMatchScore Not used at the moment
        *  @param [in] useWeightedClustering The algorithm by default clusters the poses without waiting. A non-zero value would indicate that the pose clustering should take into account the number of votes as the weights and perform a weighted averaging instead of a simple one.
        *  \return No return value is available.
        *
        */
        void SetSearchParams(const int numPoses=5, const double positionThreshold=-1, const double rotationThreshold=-1, const double minMatchScore=0.5, const bool useWeightedClustering=false);
        
        /**
        *  \brief Trains a new model.
        *
        *  @param [in] Model The input point cloud with normals (Nx6)
        *  \return Returns 0 on success.
        *
        *  \details Uses the parameters set in the constructor to downsample and learn a new model
        */
        int trainModel(const Mat& Model);
        
        /**
        *  \brief Brief
        *
        *  @param [in] Scene Point cloud for the scene
        *  @param [out] results List of output poses
        *  @param [in] RelativeSceneSampleStep The ratio of scene points to be used for the matching. For example, if this value is set to 1.0/5.0, every 5th point from the scene is used for pose refinement. This parameter allows an easy tradeoff between speed and accuracy of the matching. Increasing the value leads to less points being used and in turn to a faster but less accurate pose computation. Decreasing the value has the inverse effect.
        *  @param [in] RelativeSceneDistance Set the distance threshold relative to the diameter of the model. Only scene points that are closer to the object than this distance are used for the optimization. Scene points further away are ignored.
        *  \return Return_Description
        *
        *  \details Details
        */
        void match(const Mat& Scene, vector < Pose3D* >& results, const double RelativeSceneSampleStep=1.0/5.0, const double RelativeSceneDistance=0.03);
        
        void read(const FileNode& fn);
        void write(FileStorage& fs) const;
        
    protected:
    
        int magic;
        double maxDist, angle_step, angleStepRadians, distance_step;
        double samplingStepRelative, angleStepRelative, distanceStepRelative;
        Mat inputPC, sampledPC, PPF;
        int n, num_ref_points, sampled_step, ppf_step;
        hashtable_int* hash_table;
        THash* hash_nodes;
        
        int NumPoses;
        double PositionThreshold, RotationThreshold, MinMatchScore;
        bool UseWeightedAvg;
        
        float sampleStepSearch;
        int SceneSampleStep;
        
        void clearTrainingModels();
        
    private:
        void computePPFFeatures( const double p1[4], const double n1[4],
                                 const double p2[4], const double n2[4],
                                 double f[4]);
                                 
        bool matchPose(const Pose3D& sourcePose, const Pose3D& targetPose);
        
        int clusterPoses(Pose3D** poseList, int numPoses, vector < Pose3D* >& finalPoses);
        
        bool trained;
};
};
};


#endif