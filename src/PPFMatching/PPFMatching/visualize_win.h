#ifndef __VISUALIZE_WIN_H_
#define __VISUALIZE_WIN_H_

//#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "t_octree.h"

using namespace cv;

#if defined (__cplusplus)
extern "C" {
#endif

	void* visualize_pc(Mat cloud, int withNormals, int withBbox, int withOctree, char* Title);

#if defined (__cplusplus)
}
#endif

#endif
