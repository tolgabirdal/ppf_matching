#ifndef __VISUALIZE_WIN_H_
#define __VISUALIZE_WIN_H_

//#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "t_octree.h"

#if defined (__cplusplus)
extern "C" {
#endif

	void* visualize_pc(cv::Mat cloud, int withNormals, int withBbox, int withOctree, char* Title);
	void* visualize_registration(cv::Mat pc1, cv::Mat pc2, char* Title, int waitMS=-1);

#if defined (__cplusplus)
}
#endif

#endif
