
//#include "t_icp.h"
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
	const char* fn = "../../../data/SpaceTime/Scena1/scene1-model1_0_ascii.ply";

	TMesh* mesh = 0;
	read_mesh_ply(&mesh, fn);
	Mat pc = Mat(mesh->NumVertices, 6, CV_32F);
	get_mesh_vertices(mesh, pc.data, pc.step, 1, 1);

	const char* outputResultFile = "../../../data/out/PPFICPOutputFrog.ply";
	const char* scenePCFile = "../../../data/out/PPFICPSceneFrog.ply";

	printf("Training...");
	int64 tick1 = cv::getTickCount();
//	Mat PPFMAt = train_pc_ppf(pc, 0.05, 0.03, 30, &ppfModel);
	ppf_match_3d::PPF3DDetector detector(0.05, 0.03);
	detector.trainModel(pc);
	
	int64 tick2 = cv::getTickCount();
	printf("\nTraining complete in %f ms. Loading model...", (double)(tick2-tick1)/ cv::getTickFrequency());

	numVert = 264310;
	fn = "../../../data/kinect/scene/frog_scene_5_ascii.ply";
	Mat pcTest = load_ply_simple(fn, numVert, 1);
	printf("\nStarting matching...");

	tick1 = cv::getTickCount();
	vector < Pose3D* > results;
	detector.match(pcTest, results, 1.0/5.0,0.03);
	//t_match_pc_ppf(pcTest, 15, 5, 0.03, ppfModel, results);
	tick2 = cv::getTickCount();
	printf("\nPPF Elapsed Time %f sec\n", (double)(tick2-tick1)/ cv::getTickFrequency());

	printf("Estimated Poses:\n");

	ICP icp;

	int64 t1 = cv::getTickCount();
//	icp.registerModelToScene(pc, pcTest, results);
	float residualOutput = 0;
	icp.registerModelToScene(pc, pcTest, results);
	int64 t2 = cv::getTickCount();

	printf("Elapsed Time: %f\n\n", (double)(t2-t1)/cv::getTickFrequency());

	// debug first five poses
	for (int i=0; i<std::min(static_cast<size_t>(10), results.size()); i++)
	{
		Pose3D* pose = results[i];

		printf("Pose Result %d:\n");
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
