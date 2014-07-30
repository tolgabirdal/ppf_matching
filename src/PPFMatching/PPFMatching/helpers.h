
#ifndef __HELPERS_H_
#define __HELPERS_H_

#include <opencv2/core/core.hpp>
//#include <opencv2/flann/flann.hpp>
#include "t_octree.h"

//using namespace cv;

//typedef cvflann::L2<float> Distance_32F;

#if defined (__cplusplus)
extern "C" {
#endif

	cv::Mat load_ply_simple(const char* fileName, int numVertices, int withNormals);
	void write_ply(cv::Mat PC, const char* FileName);
	
	cv::Mat sample_pc_uniform(cv::Mat PC, int sampleStep);
	cv::Mat sample_pc_uniform_ind(cv::Mat PC, int sampleStep, std::vector<int>& indices);
	cv::Mat sample_pc_perfect_uniform(cv::Mat PC, int sampleStep);
	cv::Mat sample_pc_random(cv::Mat PC, int numPoints);
	cv::Mat sample_pc_octree(cv::Mat pc, float xrange[2], float yrange[2], float zrange[2], float resolution);
	cv::Mat sample_pc_by_quantization(cv::Mat pc, float xrange[2], float yrange[2], float zrange[2], float sample_step_relative, int weightByCenter=0);
	cv::Mat sample_pc_kd_tree(cv::Mat pc, float radius, int numNeighbors);

	void compute_bbox_std(cv::Mat pc, float xRange[2], float yRange[2], float zRange[2]);
	//double compute_diameter(cv::Mat pc);
	//void compute_obb(cv::Mat pc, float xRange[2], float yRange[2], float zRange[2]);
	
	TOctreeNode* Mat2Octree(cv::Mat pc);
	void* index_pc_flann(cv::Mat pc);
	void destroy_flann(void* flannIndex);
	void query_pc_flann(void* flannIndex, cv::Mat PC, cv::Mat& Indices, cv::Mat& Distances);
	
	cv::Mat normalize_pc(cv::Mat pc, float scale);
	cv::Mat normalize_pc_coeff(cv::Mat pc, float scale, float* Cx, float* Cy, float* Cz, float* MinVal, float* MaxVal);
	cv::Mat trans_pc_coeff(cv::Mat pc, float scale, float Cx, float Cy, float Cz, float MinVal, float MaxVal);
	cv::Mat transform_pc_pose(cv::Mat pc, double Pose[16]);

	void get_random_pose(double Pose[16]);
	void generate_random_pose(double Pose[16], double scale);
	cv::Mat add_noise_pc(cv::Mat pc, double scale);

	int compute_normals_pc_3d(const cv::Mat PC, cv::Mat& PCNormals, const int NumNeighbors, const bool FlipViewpoint, const double viewpoint[3]);


#if defined (__cplusplus)
}
#endif

#endif