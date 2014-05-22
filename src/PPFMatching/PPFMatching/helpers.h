
#ifndef __HELPERS_H_
#define __HELPERS_H_

#include <opencv2/core.hpp>

using namespace cv;

#if defined (__cplusplus)
extern "C" {
#endif

	Mat load_ply_simple(const char* fileName, int numVertices, int withNormals);
	void* visualize_pc(Mat cloud, int camera_pov);
	

#if defined (__cplusplus)
}
#endif

#endif