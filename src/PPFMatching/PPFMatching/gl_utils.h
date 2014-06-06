
#ifndef __GL_UTILSH_
#define __GL_UTILSH_

//#include <opencv2/core.hpp>
#include "t_octree.h"

#if defined(__cplusplus)
extern "C" {
#endif 

	void draw_cube(double center_x, double center_y, double center_z, double size);
	void draw_bbox(float xrange[2], float yrange[2], float zrange[2]);
	int draw_octree(const TOctreeNode* tree, float center_x, float center_y, float center_z, float cube_sizex, float cube_sizey, float cube_sizez);

	//int t_set_camera_gl(const Mat P, const int width, const int height, const float zNear, const float zFar);

#if defined(__cplusplus)
}
#endif 

#endif