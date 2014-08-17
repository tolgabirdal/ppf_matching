
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

			void update_pose(double NewPose[16]);
			void update_pose(double NewR[9], double NewT[3]);
			void update_pose_quat(double Q[4], double NewT[3]);
			void append_pose(double IncrementalPose[16]);
			void print_pose();

			Pose3D* clone();

			int write_pose(FILE* f);
			int read_pose(FILE* f);
			int write_pose(const std::string& FileName);
			int read_pose(const std::string& FileName);

			~Pose3D(){};

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

			~PoseCluster3D()
			{
				numVotes=0;
				id=0;
				//poseList.clear();
			};

			void add_pose(Pose3D* newPose) ;
			int write_pose_cluster(FILE* f);
			int read_pose_cluster(FILE* f);
			int write_pose_cluster(const std::string& FileName);
			int read_pose_cluster(const std::string& FileName);

			std::vector < Pose3D* > poseList;
			int numVotes;
			int id;
		};
	}
}

#endif
