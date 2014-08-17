
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

		typedef struct THash {
			int id;
			int i, ppfInd;
		} THash;

		class CV_EXPORTS PPF3DDetector
		{
		public:

			/**
			* \brief Empty constructor, initialize with read().
			*/
			PPF3DDetector()
			{
				samplingStepRelative = 0.05;
				distanceStepRelative = 0.05;
				SceneSampleStep = 1/0.04;
				angleStepRelative = 30;
				angleStepRadians = (360.0/angleStepRelative)*PI/180.0;
				angle_step = angleStepRadians;
				trained = false;

				SetSearchParams();
			}
			
			PPF3DDetector(const double RelativeSamplingStep, const double RelativeDistanceStep=0.05, const double NumAngles=30)
			{
				samplingStepRelative = RelativeSamplingStep;
				distanceStepRelative = RelativeDistanceStep;
				angleStepRelative = NumAngles;
				angleStepRadians = (360.0/angleStepRelative)*PI/180.0;
				//SceneSampleStep = 1.0/RelativeSceneSampleStep;
				angle_step = angleStepRadians;
				trained = false;

				SetSearchParams();
			};

			void SetSearchParams(const int numPoses=5, const double positionThreshold=-1, const double rotationThreshold=-1, const double minMatchScore=0.5, const bool useWeightedClustering=false)
			{
				NumPoses=numPoses;

				if (positionThreshold<0)
					PositionThreshold = samplingStepRelative;
				else
					PositionThreshold = positionThreshold;

				if (rotationThreshold<0)
					RotationThreshold = ((360/angle_step) / 180.0 * M_PI);
				else
					RotationThreshold = rotationThreshold;

				UseWeightedAvg = useWeightedClustering;
				MinMatchScore = minMatchScore;
			}

			~PPF3DDetector();

			int trainModel(const Mat& Model);

			//int matchModel(const Mat Scene, const double samplingStepRelative, const double distanceStepRelative, const double angleStepRelative);
			void match(const Mat& Scene, vector < Pose3D* >& results, const double RelativeSceneSampleStep=1.0/5.0, const double RelativeSceneDistance=0.03);

			void read(const FileNode& fn);
			void write(FileStorage& fs) const;

			String readClass(const FileNode& fn, const String &class_id_override = "");
			void writeClass(const String& class_id, FileStorage& fs) const;

			void readClasses(	const std::vector<String>& class_ids,
								const String& format = "templates_%s.yml.gz");

			void writeClasses(const String& format = "templates_%s.yml.gz") const;

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

			static int qsort_pose_cmp (const void * a, const void * b)
			{
				Pose3D* pose1 = *(Pose3D**)a;
				Pose3D* pose2 = *(Pose3D**)b;
				return ( pose2->numVotes - pose1->numVotes );
			}

			static int sort_pose_clusters (const PoseCluster3D* a, const PoseCluster3D* b)
			{
				return ( a->numVotes > b->numVotes );
			}

			void clear_training_models();

		private:
			void compute_ppf_features(	const double p1[4], const double n1[4],
										const double p2[4], const double n2[4],
										double f[4]);
			
			bool match_pose(const Pose3D& sourcePose, const Pose3D& targetPose);
			//int qsort_pose_cmp (const void * a, const void * b);
			//int sort_pose_clusters (const PoseCluster3D* a, const PoseCluster3D* b);

			int cluster_poses(Pose3D** poseList, int numPoses, vector < Pose3D* >& finalPoses);

			bool trained;
		};
	};
};

#endif