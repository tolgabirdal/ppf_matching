
#ifndef __HELPERS_H_
#define __HELPERS_H_

#include <opencv2/core.hpp>

using namespace cv;

#if defined (__cplusplus)
extern "C" {
#endif

	Mat load_ply_simple(const char* fileName, int numVertices, int withNormals);
	void* visualize_pc(Mat cloud, int withNormals, int withBbox, char* Title);
	
	Mat sample_pc_uniform(Mat PC, int sampleStep);
	Mat sample_pc_perfect_uniform(Mat PC, int sampleStep);
	Mat sample_pc_random(Mat PC, int numPoints);

	void compute_obb(Mat pc, float xRange[2], float yRange[2], float zRange[2]);

#if defined (__cplusplus)
}
#endif

#endif