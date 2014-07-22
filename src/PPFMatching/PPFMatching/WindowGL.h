#pragma once

#ifndef __WINDOWGL_H
#define __WINDOWGL_H

#include <windows.h>
// gl.h
#include <gl/gl.h>
//#include <gl/glaux.h>
#include <gl/glu.h>
//#include "TMesh3D.h"

typedef struct DisplayNode
{
	int data;
	struct DisplayNode *next; /* pointer to next element in list */
} DisplayNode;

typedef struct {
	float red;
	float green;
	float blue;
	float dummy;
} TWinColor;

enum ImageDrawMode {DrawModeSTRETCH, DrawModeNORMAL};
enum RegionDrawMode {RegionDrawModeFILLED, RegionDrawModeMARGIN, RegionDrawModeWIREFRAME};

typedef int (__cdecl *TGLOpCallback)(void* UserData);

typedef struct TWindowGL
{
	int magicVal;
	char* title;
	HDC			hDC;									// Private GDI Device Context
	HGLRC		hRC;									// Permanent Rendering Context
	HWND		hWnd;									// Holds Our Window Handle
	HINSTANCE	hInstance;							// Holds The Instance Of The Application
	DisplayNode *dispList;
	int length;
	int		keys[256];									// Array Used For The Keyboard Routine
	int		active;								// Window Active Flag Set To TRUE By Default
	int width, height;
	int contentWidth;
	int contentHeight;
	short   line_width; /* width for linedrawing */
	int     crs_x;            /* Textcursor-Position X */
	int     crs_y;            /* Textcursor-Position Y */
	float R,G,B;
	HFONT font;
	GLuint _font_base;
	GLuint texture;
	int xOffset, yOffset;
	int ImageDrawMode;
	int RegionDrawMode;
	TWinColor color;

	struct TWindowGL* next;
	int displaySize;
	void* preservedData;

	float rx, ry, rz, angle;
	float grx, gry, grz, gangle;
	float sx, sy, sz;
	float tx, ty, tz;
	float rot[4], globalRot[4];

	int tracking;

	TGLOpCallback PaintCallback;

} TWindowGL;

__inline static DisplayNode *list_add(DisplayNode **p, int i)
{
	DisplayNode *n = (DisplayNode*)malloc(sizeof(DisplayNode));
	if (n == NULL)
		return NULL;

	n->next = *p; /* the previous element (*p) now becomes the "next" element */
	*p = n;       /* add new empty element to the front (head) of the list */
	n->data = i;

	return *p;
}

__inline static void list_remove(DisplayNode **p) /* remove head */
{
	if (*p != NULL)
	{
		struct DisplayNode *n = *p;
		*p = (*p)->next;
		free(n);
	}
}

__inline static void clear_list(DisplayNode **p)
{
	DisplayNode* temp=*p;
	while(temp->next!=NULL)
	{
		DisplayNode* next=temp->next;
		list_remove(&temp);
		temp=0;
		temp=next;
	}
	//list_remove(&temp);
}

#if defined (__cplusplus)
extern "C" {
#endif

	int wait_window(TWindowGL* window);
	int wait_window_ms(TWindowGL* window, int miliSecs);
	int update_window(TWindowGL* window);

	static LRESULT	CALLBACK WindowProc(HWND, UINT, WPARAM, LPARAM);	// Declaration For WndProc
	BOOL CreateGLWindow(TWindowGL*, char* title, int X, int Y, int width, int height, int bits);
	BOOL MoveGLWindow(TWindowGL* window, const int x, const int y);
	GLvoid KillGLWindow(TWindowGL* window);
	int DrawGLSceneDefault(GLvoid)	;
	void disp_points(TWindowGL* window, float* x, float* y, const int length);
	void disp_rectangle(TWindowGL* window, float x1, float y1, float x2, float y2);
	void disp_line(TWindowGL* window, float x1, float y1, float x2, float y2);
	void disp_circle(TWindowGL* window, const float x, const float y, const float radius);
	void disp_polygon(TWindowGL* window, float* x, float* y, const int length);
	void disp_arrow(TWindowGL* window, float x1, float y1, float x2, float y2, float size);
	void disp_image(TWindowGL* window, const int width, const int height, const unsigned char* data);
	//void disp_mesh(TWindowGL* window, TMesh3D *mesh, int useColors);
	void t_write_string_gl_t(TWindowGL* window, const char* string, const int x, const int y);

	void draw_custom_gl_scene(TWindowGL* window, TGLOpCallback Callback, void* Data);
	void register_custom_gl_scene(TWindowGL* window, TGLOpCallback Callback, void* Data);

	//void draw_mesh_3D_gl(TMesh3D *mesh, int useColors);

	void set_color(TWindowGL* win, unsigned char R, unsigned char G, unsigned char B);
	void set_color_string(TWindowGL* win, const char* color);
	void set_line_width(TWindowGL* win, int lineWidth);
	void set_font(TWindowGL* window, const char* fontString);
	void set_font_default(TWindowGL* window);
	void set_transform(TWindowGL* window, float sx, float sy, float sz, float rx, float ry, float rz);
	void set_rotation(TWindowGL* window, float rx, float ry, float rz);
	void set_translation(TWindowGL* window, float tx, float ty, float tz);
	void set_scale(TWindowGL* window, float sx, float sy, float sz);

	void clear_window(TWindowGL* window);
	int close_window(TWindowGL* window);
	float GetOpenGLVersion ( void );
	void get_pixel_buffer(unsigned char* dst, const int x, const int y, const int w, const int h, const int ws);

	void glEnable3D( float FOVy, float zNear, float zFar, int width, int height );
	void glDisable3D( void );

	void disp_points_3d(TWindowGL* window, float* x, float* y, float* z, const int length);

#if defined (__cplusplus)
}
#endif

#endif