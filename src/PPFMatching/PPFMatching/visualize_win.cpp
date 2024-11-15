
#include "WindowGL.h"
#include "visualize_win.h"
#include "ppf_helpers.hpp"
//#include "gl_utils.h"

using namespace cv;
using namespace cv::ppf_match_3d;

typedef struct
{
	Mat PC;
	Mat PC2;
//	TOctreeNode* octree;
	TWindowGL* window;
	int withNormals, withBbox, withOctree;
}TWindowData;

static float light_diffuse[] = {100.0f, 100.0f, 100.0f, 100.0f}; 
static float light_position[] = {-100.0f, 100.0f, 1.0f, 0.0f};
static float amb[] =  {0.24f, 0.24f, 0.24f, 0.0f};
static float dif[] =  {1.0f, 1.0f, 1.0f, 0.0f};


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
    float half_sizex = sizex / 2.0f;
    float half_sizey = sizey / 2.0f;
    float half_sizez = sizez / 2.0f;
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

int display(void* UserData)
{
	TWindowData* wd = (TWindowData*)UserData;

	Mat pc = wd->PC;
	TWindowGL* window = wd->window;
	int withNormals = wd->withNormals;
	int withBbox = wd->withBbox;
	int withOctree = wd->withOctree;

	double minVal = 0, maxVal = 0;
	double diam=5;
	
	Mat x,y,z, pcn;
	pc.col(0).copyTo(x);
	pc.col(1).copyTo(y);
	pc.col(2).copyTo(z);

	float cx = (float) cv::mean(x).val[0];
	float cy = (float) cv::mean(y).val[0];
	float cz = (float) cv::mean(z).val[0];

	cv::minMaxIdx(pc, &minVal, &maxVal);

	x=x-cx;
	y=y-cy;
	z=z-cz;
	pcn.create(pc.rows, 3, CV_32FC1);
	x.copyTo(pcn.col(0));
	y.copyTo(pcn.col(1));
	z.copyTo(pcn.col(2));

	cv::minMaxIdx(pcn, &minVal, &maxVal);
	pcn=(float)diam*(pcn)/((float)maxVal-(float)minVal);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1.0, 1, diam*4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	gluLookAt (-diam, diam, diam/2,
                0,0, 0,
                0.0, 1.0, 0.0);

	if (window->tracking)
		glRotatef(window->angle, window->rx, window->ry, window->rz);

	glRotatef(window->gangle, window->grx, window->gry, window->grz);
	//else

	/*glRotatef(80,0,0,1);
	glRotatef(40,0,1,0);*/

	// set lights
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);

	//glColor4f(1,1,1,1);
	int glSize = 2;
	glLineWidth(glSize);
	glPointSize(glSize);

	/*glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);*/
	glDisable(GL_LIGHTING);

	glBegin(GL_POINTS);

	for (int i=0; i!=pcn.rows; i++)
	{
		float* dataPC = (float*)(&pc.data[i*pc.step[0]]);
		float* data = (float*)(&pcn.data[i*pcn.step[0]]);
		if (withNormals)
		{
			glColor4f(dataPC[3],dataPC[4],dataPC[5],1);
			glNormal3f(dataPC[3], dataPC[4], dataPC[5]);
		}

		glVertex3f(data[0], data[1], data[2]);
	}

	glEnd();

	float xRange[2], yRange[2], zRange[2];
	if (withBbox || withOctree)
	{
		//compute_obb(pcn, xRange, yRange, zRange);
		computeBboxStd(pcn, xRange, yRange, zRange);
	}

	if (withBbox)
	{
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHTING);
		glColor4f(0,0.9,0.1,0.75);		
		
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

	if (withOctree)
	{
	/*	float cx = (xRange[1] + xRange[0])*0.5f;
		float cy = (yRange[1] + yRange[0])*0.5f;
		float cz = (zRange[1] + zRange[0])*0.5f;
		wd->octree = Mat2Octree(pcn);
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHTING);
		glSize = 1;
		glLineWidth(glSize);
		glPointSize(glSize);
		glColor4f(1,1,1,1);
		draw_octree(wd->octree, cx, cy, cz, xRange[1] - xRange[0], yRange[1] - yRange[0], zRange[1] - zRange[0]);
		t_octree_destroy(wd->octree);*/
	}

	//SwapBuffers(window->hDC);

	return 0;
}

int display_registration(void* UserData)
{
	TWindowData* wd = (TWindowData*)UserData;

	Mat pc = wd->PC;
	Mat pc2 = wd->PC2;
	TWindowGL* window = wd->window;
	double diam=5;
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1.0, 1, diam*4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	gluLookAt (-diam, diam, diam/2,
                0,0, 0,
                0.0, 1.0, 0.0);

	if (window->tracking)
		glRotatef(window->angle, window->rx, window->ry, window->rz);

	glRotatef(window->gangle, window->grx, window->gry, window->grz);
	//else

	/*glRotatef(80,0,0,1);
	glRotatef(40,0,1,0);*/

	// set lights
    glEnable(GL_LIGHTING);
    //glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);
    glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);

	//glColor4f(1,1,1,1);
	int glSize = 2;
	glLineWidth(glSize);
	glPointSize(glSize);

	/*glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);*/
	glDisable(GL_LIGHTING);

	glColor4f(0.65,0.65,0.65,0.65);

	glBegin(GL_POINTS);
	for (int i=0; i!=pc.rows; i++)
	{
		float* data = (float*)(&pc.data[i*pc.step[0]]);
		glNormal3f(data[3], data[4], data[5]);
		glVertex3f(data[0], data[1], data[2]);
	}
	glEnd();

	glColor4f(0,0.8,0.2,0.45);
	glBegin(GL_POINTS);
	for (int i=0; i!=pc2.rows; i++)
	{
		float* data = (float*)(&pc2.data[i*pc2.step[0]]);
		glNormal3f(data[3], data[4], data[5]);
		glVertex3f(data[0], data[1], data[2]);
	}
	glEnd();

	return 0;
}

void* visualize_pc(Mat pc, int withNormals, int withBbox, int withOctree, char* Title)
{
	int width = 1024;
	int height = 1024;

	TWindowGL* window=(TWindowGL*)calloc(1, sizeof(TWindowGL));
	int status=CreateGLWindow(window, Title, 300, 300, width, height, 24);

	MoveGLWindow(window, 300, 300);
	update_window(window);

	glEnable3D(45, 1, 5, width, height);

	TWindowData* wd = new TWindowData();
	wd->PC = pc;
	wd->window = window;
	wd->withNormals=withNormals;
	wd->withBbox=withBbox;
	wd->withOctree = withOctree;
//	wd->octree=0;

	if (withOctree)
	{
//		wd->octree = Mat2Octree(pc);
	}

	//draw_custom_gl_scene(window, display, wd);
	register_custom_gl_scene(window, display, wd);

	//update_window(window);
	wait_window(window);

	close_window(window);

	return (void*)window;
}


TWindowGL* window=0;
void* visualize_registration(Mat pc1, Mat pc2, char* Title, int waitMS)
{
	int width = 1024;
	int height = 1024;

	if (!window)
	{
		window=(TWindowGL*)calloc(1, sizeof(TWindowGL));
		int status=CreateGLWindow(window, Title, 300, 50, width, height, 24);
		MoveGLWindow(window, 300, 50);
		update_window(window);
	}


	glEnable3D(45, 1, 5, width, height);

	TWindowData* wd = new TWindowData();
	
	float cx=0, cy=0, cz=0, minv=0, maxv=0;
	wd->PC = normalizePCCoeff(pc1, 5, &cx, &cy, &cz, &minv, &maxv);
	wd->PC2 = transPCCoeff(pc2, 5, cx, cy, cz, minv, maxv);
	//wd->PC2 = normalize_pc(pc2, 5);

	wd->window = window;
	
	//draw_custom_gl_scene(window, display, wd);
	register_custom_gl_scene(window, display_registration, wd);

	//update_window(window);
	if (waitMS <=0)
		wait_window(window);
	else
		wait_window_ms(window, waitMS);

	//close_window(window);

	

	return (void*)window;
}


