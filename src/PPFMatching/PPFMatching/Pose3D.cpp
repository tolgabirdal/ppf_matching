

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION > 2
  #include <opencv2/core/utility.hpp>
#else
  #include <opencv2/core/core.hpp>
#endif

#include "Pose3D.hpp"

using namespace std;

namespace cv 
{
	namespace ppf_match_3d 
	{

		void Pose3D::update_pose(double NewPose[16])
		{
			double R[9];

			for (int i=0; i<16; i++)
				Pose[i]=NewPose[i];

			R[0] = Pose[0];	R[1] = Pose[1]; R[2] = Pose[2];
			R[3] = Pose[4];	R[4] = Pose[5]; R[5] = Pose[6];
			R[6] = Pose[8];	R[7] = Pose[9]; R[8] = Pose[10];

			t[0]=Pose[3]; t[1]=Pose[7]; t[2]=Pose[11];

			// compute the angle
			const double trace = R[0] + R[4] + R[8];

			if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
			else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
			else							 {angle = ( acos((trace - 1)/2) ); }

			// compute the quaternion
			matrix_to_quaternion(R, q);
		}

		void Pose3D::append_pose(double IncrementalPose[16])
		{
			double R[9], PoseFull[16]={0};

			matrix_product44(IncrementalPose, this->Pose, PoseFull);

			R[0] = PoseFull[0];	R[1] = PoseFull[1]; R[2] = PoseFull[2];
			R[3] = PoseFull[4];	R[4] = PoseFull[5]; R[5] = PoseFull[6];
			R[6] = PoseFull[8];	R[7] = PoseFull[9]; R[8] = PoseFull[10];

			t[0]=PoseFull[3]; t[1]=PoseFull[7]; t[2]=PoseFull[11];

			// compute the angle
			const double trace = R[0] + R[4] + R[8];

			if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
			else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
			else							 {angle = ( acos((trace - 1)/2) ); }

			// compute the quaternion
			matrix_to_quaternion(R, q);
		}

		void Pose3D::update_pose(double NewR[9], double NewT[3])
		{
			//for (int i=0; i<16; i++)
			//	Pose[i]=NewPose[i];
			Pose[0]=NewR[0]; Pose[1]=NewR[1]; Pose[2]=NewR[2]; Pose[3]=NewT[0];
			Pose[4]=NewR[3]; Pose[5]=NewR[4]; Pose[6]=NewR[5]; Pose[7]=NewT[1];
			Pose[8]=NewR[6]; Pose[9]=NewR[7]; Pose[10]=NewR[8]; Pose[11]=NewT[2];
			Pose[12]=0;		 Pose[13]=0;	  Pose[14]=0;	   Pose[15]=1;

			// compute the angle
			const double trace = NewR[0] + NewR[4] + NewR[8];

			if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
			else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
			else							 {angle = ( acos((trace - 1)/2) ); }

			// compute the quaternion
			matrix_to_quaternion(NewR, q);
		}

		void Pose3D::update_pose_quat(double Q[4], double NewT[3])
		{
			double NewR[9];

			quaternion_to_matrix(Q, NewR);
			q[0]=Q[0]; q[1]=Q[1]; q[2]=Q[2]; q[3]=Q[3]; 

			//for (int i=0; i<16; i++)
			//	Pose[i]=NewPose[i];
			Pose[0]=NewR[0]; Pose[1]=NewR[1]; Pose[2]=NewR[2]; Pose[3]=NewT[0];
			Pose[4]=NewR[3]; Pose[5]=NewR[4]; Pose[6]=NewR[5]; Pose[7]=NewT[1];
			Pose[8]=NewR[6]; Pose[9]=NewR[7]; Pose[10]=NewR[8]; Pose[11]=NewT[2];
			Pose[12]=0;		 Pose[13]=0;	  Pose[14]=0;	   Pose[15]=1;

			// compute the angle
			const double trace = NewR[0] + NewR[4] + NewR[8];

			if (fabs(trace - 3) <= EPS)		 { angle = 0;	}
			else if (fabs(trace + 1) <= EPS) { angle = M_PI;	}
			else							 {angle = ( acos((trace - 1)/2) ); }
		}

		Pose3D* Pose3D::clone()
		{
			Pose3D* pose = new Pose3D(alpha, modelIndex, numVotes);
			for (int i=0; i<16; i++)
				pose->Pose[i]= Pose[i];

			pose->q[0]=q[0];
			pose->q[1]=q[1];
			pose->q[2]=q[2];
			pose->q[3]=q[3];

			pose->t[0]=t[0];
			pose->t[1]=t[1];
			pose->t[2]=t[2];

			pose->angle=angle;

			return pose;
		}

		void Pose3D::print_pose()
		{
			printf("\n-- Pose to Model Index %d: NumVotes = %d, Residual = %f\n", this->modelIndex, this->numVotes, this->residual);
			for (int j=0; j<4; j++)
			{
				for (int k=0; k<4; k++)
				{
					printf("%f ", this->Pose[j*4+k]);
				}
				printf("\n");
			}
			printf("\n");
		}

		int Pose3D::write_pose(FILE* f)
		{
			int POSE_MAGIC = 7673;
			fwrite(&POSE_MAGIC, sizeof(int), 1, f);
			fwrite(&angle, sizeof(double), 1, f);
			fwrite(&numVotes, sizeof(int), 1, f);
			fwrite(&modelIndex, sizeof(int), 1, f);
			fwrite(Pose, sizeof(double)*16, 1, f);
			fwrite(t, sizeof(double)*3, 1, f);
			fwrite(q, sizeof(double)*4, 1, f);
			fwrite(&residual, sizeof(double), 1, f);
			return 0;
		}

		int Pose3D::read_pose(FILE* f)
		{
			int POSE_MAGIC = 7673, magic;

			fread(&magic, sizeof(int), 1, f);
			if (magic == POSE_MAGIC)
			{
				fread(&angle, sizeof(double), 1, f);
				fread(&numVotes, sizeof(int), 1, f);
				fread(&modelIndex, sizeof(int), 1, f);
				fread(Pose, sizeof(double)*16, 1, f);
				fread(t, sizeof(double)*3, 1, f);
				fread(q, sizeof(double)*4, 1, f);
				fread(&residual, sizeof(double), 1, f);
				return 0;
			}

			return -1;
		}

		int Pose3D::write_pose(const std::string& FileName)
		{
			FILE* f = fopen(FileName.c_str(), "wb");

			if (!f)
				return -1;

			int status = write_pose(f);

			fclose(f);
			return status;
		}

		int Pose3D::read_pose(const std::string& FileName)
		{
			FILE* f = fopen(FileName.c_str(), "rb");

			if (!f)
				return -1;

			int status = read_pose(f);

			fclose(f);
			return status;
		}


		void PoseCluster3D::add_pose(Pose3D* newPose) 
		{
			poseList.push_back(newPose);
			this->numVotes += newPose->numVotes;
		};

		int PoseCluster3D::write_pose_cluster(FILE* f)
		{
			int POSE_CLUSTER_MAGIC_IO = 8462597;
			fwrite(&POSE_CLUSTER_MAGIC_IO, sizeof(int), 1, f);
			fwrite(&id, sizeof(int), 1, f);
			fwrite(&numVotes, sizeof(int), 1, f);

			int numPoses = poseList.size();
			fwrite(&numPoses, sizeof(int), 1, f);

			for (int i=0; i<numPoses; i++)
				poseList[i]->write_pose(f);

			return 0;
		}

		int PoseCluster3D::read_pose_cluster(FILE* f)
		{
			int POSE_CLUSTER_MAGIC_IO = 8462597;
			int magic=0, numPoses=0;
			fwrite(&magic, sizeof(int), 1, f);

			if (magic==POSE_CLUSTER_MAGIC_IO)
				return -1;

			fread(&id, sizeof(int), 1, f);
			fread(&numVotes, sizeof(int), 1, f);
			fread(&numPoses, sizeof(int), 1, f);

			poseList.clear();
			poseList.resize(numPoses);
			for (int i=0; i<poseList.size(); i++)
			{
				poseList[i] = new Pose3D();
				poseList[i]->read_pose(f);
			}
		}

		int PoseCluster3D::write_pose_cluster(const std::string& FileName)
		{
			FILE* f = fopen(FileName.c_str(), "wb");

			if (!f)
				return -1;

			int status = write_pose_cluster(f);

			fclose(f);
			return status;
		}

		int PoseCluster3D::read_pose_cluster(const std::string& FileName)
		{
			FILE* f = fopen(FileName.c_str(), "rb");

			if (!f)
				return -1;

			int status = read_pose_cluster(f);

			fclose(f);
			return status;
		}
	}
}