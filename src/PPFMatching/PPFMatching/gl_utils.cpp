
#if defined WIN32
#include "windows.h"
#endif

#include "gl_utils.h"
#include <gl/gl.h>
//#include <gl/glu.h>

void draw_cube(float center_x, float center_y, float center_z, float size)
{
    float half_size = size / 2.0;
    float front     = center_z - half_size;
    float back      = center_z + half_size;
    float left      = center_x - half_size;
    float right     = center_x + half_size;
    float bottom    = center_y - half_size;
    float top       = center_y + half_size;

    glPushMatrix();
	glPolygonMode(GL_BACK, GL_LINE); 
	glPolygonMode(GL_FRONT, GL_LINE); 

    // red side - front
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 0.0, 0.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(left,  top,    front);
    glVertex3f(left,  bottom, front);
    glEnd();

    // green side - back
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 1.0, 0.0);
    glVertex3f(right, bottom, back);
    glVertex3f(right, top,    back);
    glVertex3f(left,  top,    back);
    glVertex3f(left,  bottom, back);
    glEnd();

    // blue side - right
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(right, top,    back);
    glVertex3f(right, bottom, back);
    glEnd();

    // yellow side - left
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 1.0, 0.0);
    glVertex3f(left, bottom, back);
    glVertex3f(left, top,    back);
    glVertex3f(left, top,    front);
    glVertex3f(left, bottom, front);
    glEnd();

    // magneta side - top
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 0.0, 1.0);
    glVertex3f(right, top, back);
    glVertex3f(right, top, front);
    glVertex3f(left,  top, front);
    glVertex3f(left,  top, back);
    glEnd();

    // cyan side - bottom
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 1.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, back);
    glVertex3f(left,  bottom, back);
    glVertex3f(left,  bottom, front);
    glEnd();

    glPopMatrix();
}

void draw_prism(float center_x, float center_y, float center_z, float sizex, float sizey, float sizez)
{
    float half_sizex = sizex / 2.0;
    float half_sizey = sizey / 2.0;
    float half_sizez = sizez / 2.0;
    float front     = center_z - half_sizez;
    float back      = center_z + half_sizez;
    float left      = center_x - half_sizex;
    float right     = center_x + half_sizex;
    float bottom    = center_y - half_sizey;
    float top       = center_y + half_sizey;

    glPushMatrix();
	glPolygonMode(GL_BACK, GL_LINE); 
	glPolygonMode(GL_FRONT, GL_LINE); 

    // red side - front
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 0.0, 0.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(left,  top,    front);
    glVertex3f(left,  bottom, front);
    glEnd();

    // green side - back
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 1.0, 0.0);
    glVertex3f(right, bottom, back);
    glVertex3f(right, top,    back);
    glVertex3f(left,  top,    back);
    glVertex3f(left,  bottom, back);
    glEnd();

    // blue side - right
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(right, top,    back);
    glVertex3f(right, bottom, back);
    glEnd();

    // yellow side - left
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 1.0, 0.0);
    glVertex3f(left, bottom, back);
    glVertex3f(left, top,    back);
    glVertex3f(left, top,    front);
    glVertex3f(left, bottom, front);
    glEnd();

    // magneta side - top
    glBegin(GL_POLYGON);
    //glColor3f(1.0, 0.0, 1.0);
    glVertex3f(right, top, back);
    glVertex3f(right, top, front);
    glVertex3f(left,  top, front);
    glVertex3f(left,  top, back);
    glEnd();

    // cyan side - bottom
    glBegin(GL_POLYGON);
    //glColor3f(0.0, 1.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, back);
    glVertex3f(left,  bottom, back);
    glVertex3f(left,  bottom, front);
    glEnd();

    glPopMatrix();
}

void draw_bbox(float xRange[2], float yRange[2], float zRange[2])
{

	glBegin(GL_LINES);

	glVertex3f(xRange[0], yRange[0], zRange[0]);
	glVertex3f(xRange[1], yRange[0], zRange[0]);

	glVertex3f(xRange[0], yRange[0], zRange[0]);
	glVertex3f(xRange[0], yRange[1], zRange[0]);

	glVertex3f(xRange[0], yRange[0], zRange[0]);
	glVertex3f(xRange[0], yRange[0], zRange[1]);

	glVertex3f(xRange[1], yRange[1], zRange[1]);
	glVertex3f(xRange[1], yRange[1], zRange[0]);

	glVertex3f(xRange[1], yRange[1], zRange[1]);
	glVertex3f(xRange[1], yRange[0], zRange[1]);

	glVertex3f(xRange[1], yRange[1], zRange[1]);
	glVertex3f(xRange[0], yRange[1], zRange[1]);

	glVertex3f(xRange[1], yRange[0], zRange[0]);
	glVertex3f(xRange[1], yRange[0], zRange[1]);

	glVertex3f(xRange[1], yRange[0], zRange[0]);
	glVertex3f(xRange[1], yRange[1], zRange[0]);

	glVertex3f(xRange[0], yRange[1], zRange[0]);
	glVertex3f(xRange[0], yRange[1], zRange[1]);

	glVertex3f(xRange[0], yRange[1], zRange[0]);
	glVertex3f(xRange[1], yRange[1], zRange[0]);

	glVertex3f(xRange[0], yRange[0], zRange[1]);
	glVertex3f(xRange[0], yRange[1], zRange[1]);

	glVertex3f(xRange[0], yRange[0], zRange[1]);
	glVertex3f(xRange[1], yRange[0], zRange[1]);

	glEnd();
}

int draw_octree(const TOctreeNode* tree, float center_x, float center_y, float center_z, float cube_sizex, float cube_sizey, float cube_sizez)
{
    if (!tree)
    {
        //draw_cube(center_x, center_y, center_z, cube_size);
		draw_prism(center_x, center_y, center_z, cube_sizex, cube_sizey, cube_sizez);
        return 0;
    }
    else if (tree->children[0]==NULL)
    {
        //draw_cube(center_x, center_y, center_z, cube_size);
		draw_prism(center_x, center_y, center_z, cube_sizex, cube_sizey, cube_sizez);
        return 0;
    }
	else
	{
		int i;

        draw_prism(center_x, center_y, center_z, cube_sizex, cube_sizey, cube_sizez);
		for (i=0; i!=8; i++)
		{
			//if (tree->children[i]->data)
			{
				//float new_cube_size = cube_size / 2.0;
				//float bias          = new_cube_size / 2.0;
				float left          = tree->children[i]->center[0] - tree->children[i]->halfDim[0];
				float right         = tree->children[i]->center[0] + tree->children[i]->halfDim[0];
				float bottom        = tree->children[i]->center[1] - tree->children[i]->halfDim[1];
				float top           = tree->children[i]->center[1] + tree->children[i]->halfDim[1];
				float front         = tree->children[i]->center[2] - tree->children[i]->halfDim[2];
				float back          = tree->children[i]->center[2] + tree->children[i]->halfDim[2];

				draw_octree(tree->children[i], left,  bottom, back,  2*tree->children[i]->halfDim[0],  2*tree->children[i]->halfDim[1],  2*tree->children[i]->halfDim[2]);
			}
		}
        return 0;
    }
}

//
//int t_set_camera_gl(const Mat P, const int width, const int height, const float zNear, const float zFar)
//{
//#define T_ACCESS_ELEM(p,y,x) p[y*w+x]
//
//	float mat[16];
//
//	float* p;
//	int w,h;
//
//	p=(float*)P.data;
//	w=P.cols;
//	h=P.rows;
//
//	if (T_ACCESS_ELEM(p,2,3) < 0.0)
//	{
//		//t_mul_image_scalar(P, &P2, -1);
//		Mat P2=-P;
//		p=(float*)P2.data;
//	}
//
//	mat[0] = T_ACCESS_ELEM(p,0,0);
//	mat[4] = T_ACCESS_ELEM(p,0,1);
//	mat[8] = T_ACCESS_ELEM(p,0,2);
//	mat[12] = T_ACCESS_ELEM(p,0,3);
//	mat[1] = T_ACCESS_ELEM(p,1,0);
//	mat[5] = T_ACCESS_ELEM(p,1,1);
//	mat[9] = T_ACCESS_ELEM(p,1,2);
//	mat[13] = T_ACCESS_ELEM(p,1,3);
//	mat[2] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,0);
//	mat[6] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,1);
//	mat[10] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,2);
//	mat[14] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,3)+zNear*zFar;
//	mat[3] = T_ACCESS_ELEM(p,2,0);
//	mat[7] = T_ACCESS_ELEM(p,2,1);
//	mat[11] = T_ACCESS_ELEM(p,2,2);
//	mat[15] = T_ACCESS_ELEM(p,2,3);
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	glOrtho(0.0, width, 0.0, height, zNear, zFar);
//	glTranslatef(0,height,0);
//	glScalef(1.0,-1.0,1.0);
//	glMultMatrixf(mat);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//
//	return 0;
//}