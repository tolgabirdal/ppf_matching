
#include "WindowGL.h"
#include "helpers.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;

Mat load_ply_simple(const char* fileName, int numVertices, int withNormals)
{
	Mat cloud;

	if (withNormals)
		cloud=Mat(numVertices, 6, CV_32FC1);
	else
		cloud=Mat(numVertices, 3, CV_32FC1);

	ifstream ifs(fileName);

	string str;
	while (str!="end_header")
		getline(ifs, str);

	float dummy =  0;
	for(size_t i = 0; i < numVertices; i++)
	{
		float* data = (float*)(&cloud.data[i*cloud.step[0]]);
		if (withNormals)
		{
			ifs >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
		}
		else
		{
			ifs >> data[0] >> data[1] >> data[2];
		}
	}

	//cloud *= 5.0f;
	return cloud;
}


int t_set_camera_gl(const Mat P, const int width, const int height, const float zNear, const float zFar)
{
#define T_ACCESS_ELEM(p,y,x) p[y*w+x]

	float mat[16];

	float* p;
	int w,h;

	p=(float*)P.data;
	w=P.cols;
	h=P.rows;

	if (T_ACCESS_ELEM(p,2,3) < 0.0)
	{
		//t_mul_image_scalar(P, &P2, -1);
		Mat P2=-P;
		p=(float*)P2.data;
	}

	mat[0] = T_ACCESS_ELEM(p,0,0);
	mat[4] = T_ACCESS_ELEM(p,0,1);
	mat[8] = T_ACCESS_ELEM(p,0,2);
	mat[12] = T_ACCESS_ELEM(p,0,3);
	mat[1] = T_ACCESS_ELEM(p,1,0);
	mat[5] = T_ACCESS_ELEM(p,1,1);
	mat[9] = T_ACCESS_ELEM(p,1,2);
	mat[13] = T_ACCESS_ELEM(p,1,3);
	mat[2] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,0);
	mat[6] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,1);
	mat[10] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,2);
	mat[14] = -(zNear+zFar)*T_ACCESS_ELEM(p,2,3)+zNear*zFar;
	mat[3] = T_ACCESS_ELEM(p,2,0);
	mat[7] = T_ACCESS_ELEM(p,2,1);
	mat[11] = T_ACCESS_ELEM(p,2,2);
	mat[15] = T_ACCESS_ELEM(p,2,3);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, width, 0.0, height, zNear, zFar);
	glTranslatef(0,height,0);
	glScalef(1.0,-1.0,1.0);
	glMultMatrixf(mat);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	return 0;
}

typedef struct
{
	Mat PC;
	TWindowGL* window;
	int withNormals;
}TWindowData;

static float light_diffuse[] = {100, 0, 100, 100.0f}; 
static float light_position[] = {-100.0, 100.0, 1.0, 0.0};
static float amb[] =  {0.4, 0.4, 0.4, 0.0};
static float dif[] =  {1.0, 1.0, 1.0, 0.0};

int display(void* UserData)
{
	TWindowData* wd = (TWindowData*)UserData;

	Mat pc = wd->PC;
	TWindowGL* window = wd->window;
	int withNormals = wd->withNormals;

	double minVal = 0, maxVal = 0;
	double diam=5;
	
	Mat x,y,z, pcn;
	pc.col(0).copyTo(x);
	pc.col(1).copyTo(y);
	pc.col(2).copyTo(z);

	float cx = cv::mean(x).val[0];
	float cy = cv::mean(y).val[0];
	float cz = cv::mean(z).val[0];

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

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1.0, 1, diam*2);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	gluLookAt (-diam, diam, diam/2,
                0,0, 0,
                0.0, 1.0, 0.0);
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

	glColor4f(1,1,1,1);
	glLineWidth(0.5);

	glBegin(GL_POINTS);

	for (int i=0; i!=pcn.rows; i++)
	{
		float* dataPC = (float*)(&pc.data[i*pc.step[0]]);
		float* data = (float*)(&pcn.data[i*pcn.step[0]]);
		if (withNormals)
			glNormal3f(dataPC[3], dataPC[4], dataPC[5]);

		glVertex3f(data[0], data[1], data[2]);
	}

	glEnd();

	return 0;
}


void* visualize_pc(Mat pc, int withNormals)
{
	TWindowGL* window=(TWindowGL*)calloc(1, sizeof(TWindowGL));
	int status=CreateGLWindow(window, "Point Cloud", 300, 300, 512, 512, 24);

	MoveGLWindow(window, 300, 300);
	update_window(window);

	glEnable3D(45, 1, 3, 512, 512);

	TWindowData* wd = new TWindowData();
	wd->PC = pc;
	wd->window = window;
	wd->withNormals=withNormals;
	draw_custom_gl_scene(window, display, wd);

	//update_window(window);
	wait_window(window);

	return (void*)window;
}