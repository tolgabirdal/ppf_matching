


#include "tmesh.h"
#include "ppf_match_3d.hpp"
#include "icp.hpp"
#include "visualize_win.h"

using namespace std;
using namespace cv;
using namespace ppf_match_3d;

// real data
int main()
{
	int numVert ;

	//const char* fn = "../../../data/kinect/model/Frog_ascii2.ply";
	//const char* fn = "../../../data/SpaceTime/Scena1/scene1-model1_0_ascii.ply";
	const char* fn = "/data/parasaurolophus_low_normals2.ply";

	// output filenames
	const char* outputResultFile = "../../../data/out/PPFICPOutputFrog.ply";
	const char* scenePCFile = "../../../data/out/PPFICPSceneFrog.ply";

	// Read the model as a mesh
	TMesh* mesh = 0;
	read_mesh_ply(&mesh, fn);
	Mat pc = Mat(mesh->NumVertices, 6, CV_32F);
	get_mesh_vertices(mesh, pc.data, pc.step, 1, 1);

	// Now train the model
	printf("Training...");
	int64 tick1 = cv::getTickCount();
	ppf_match_3d::PPF3DDetector detector(0.05, 0.05);
	detector.trainModel(pc);	
	int64 tick2 = cv::getTickCount();
	printf("\nTraining complete in %f ms. Loading model...", (double)(tick2-tick1)/ cv::getTickFrequency());

	// Read the scene

	//numVert = 264310;
	//fn = "../../../data/kinect/scene/frog_scene_5_ascii.ply";
	//numVert = 131834;
	//fn = "../../../data/SpaceTime/Scena1/scene1-scene1_0_ascii.ply";
	numVert = 113732;
	fn = "../../../data/Retrieval/rs22_proc2.ply";
	Mat pcTest = load_ply_simple(fn, numVert, 1);
	printf("\nStarting matching...");

	// Match the model to the scene and get the pose
	tick1 = cv::getTickCount();
	vector < Pose3D* > results;
	detector.match(pcTest, results, 1.0/5.0,0.05);
	//t_match_pc_ppf(pcTest, 15, 5, 0.03, ppfModel, results);
	tick2 = cv::getTickCount();
	printf("\nPPF Elapsed Time %f sec\n", (double)(tick2-tick1)/ cv::getTickFrequency());

	// Create an instance of ICP
	ICP icp(200, 0.005, 3, 8);

	int64 t1 = cv::getTickCount();
	float residualOutput = 0;

	// Get only first N results
	int N=1;
	vector<Pose3D*>::const_iterator first = results.begin();
	vector<Pose3D*>::const_iterator last = results.begin() + N;
	vector<Pose3D*> resultsSub(first, last);

	// Register for all selected poses
	icp.registerModelToScene(pc, pcTest, resultsSub);
	int64 t2 = cv::getTickCount();

	printf("Elapsed Time: %f\n\n", (double)(t2-t1)/cv::getTickFrequency());

	printf("Estimated Poses:\n");
	// debug first five poses
	for (int i=0; i<resultsSub.size(); i++)
	{
		Pose3D* pose = resultsSub[i];

		printf("Pose Result %d:\n", i);
		pose->print_pose();

		Mat pct = transform_pc_pose(pc, pose->Pose);

		// write the mesh if desired
		//TMesh* MeshTr = transform_mesh_new(mesh, pose->Pose);
		//write_ply(pct, outputResultFile);
		//t_write_mesh_ply(MeshTr, outputResultFile);
		//write_ply(pcTest, scenePCFile);
		//destroy_mesh(&MeshTr);

		// Visualize registration
#if defined (_MSC_VER)
		visualize_registration(pcTest, pct, "Registration");
#endif		
	}

	return 0;
}
