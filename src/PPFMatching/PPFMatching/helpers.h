
#ifndef __HELPERS_H_
#define __HELPERS_H_

#include <opencv2/core/core.hpp>
#include <opencv2/flann/flann.hpp>
#include "t_octree.h"

using namespace cv;

typedef cvflann::L2<float> Distance_32F;

#if defined (__cplusplus)
extern "C" {
#endif

	Mat load_ply_simple(const char* fileName, int numVertices, int withNormals);
	
	Mat sample_pc_uniform(Mat PC, int sampleStep);
	Mat sample_pc_perfect_uniform(Mat PC, int sampleStep);
	Mat sample_pc_random(Mat PC, int numPoints);
	Mat sample_pc_octree(Mat pc, float xrange[2], float yrange[2], float zrange[2], float resolution);

	void compute_obb(Mat pc, float xRange[2], float yRange[2], float zRange[2]);
	
	TOctreeNode* Mat2Octree(Mat pc);
	void* index_pc_flann(Mat pc, cvflann::Matrix<float>& data);
	
	Mat normalize_pc(Mat pc, float scale);
	Mat transform_pc_pose(Mat pc, double Pose[16]);

	void generate_random_pose(double Pose[16]);
	Mat add_noise_pc(Mat pc);


#if defined (__cplusplus)
}
#endif

#endif