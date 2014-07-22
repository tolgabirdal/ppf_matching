

#include <math.h>
#include <time.h>
#include "WindowGL.h"
#include "trackball.h"
#include "TRgbGL.h"

#ifndef GL_TEXTURE_RECTANGLE_ARB
#define GL_TEXTURE_RECTANGLE_ARB 0x84F5
#endif

#ifndef GL_BGR
#define GL_BGR 0x80E0
#endif

#ifndef PI
#ifdef  M_PI
#define PI               M_PI
#define PI_F             3.1415926535897932384626433832795f
#else
#define PI               3.1415926535897932384626433832795
#define PI_F             3.1415926535897932384626433832795f
#endif
#endif

static TWindowGL** tWindows;
TWindowGL* currentWindow;

static DWORD g_dwGL_VersionNumber = 0;

float GetOpenGLVersion ( void )
{
	const GLubyte	*szGLVer ;
	if ( g_dwGL_VersionNumber )
		return g_dwGL_VersionNumber;

	szGLVer = glGetString( GL_VERSION );

	if ( !( szGLVer==0 ) )
	{
		DWORD	dwVer = 0;
		INT		iNumDigits = 0;
		const GLubyte	*pCurr = szGLVer;

		// Get the Major number.
		while ( *pCurr != '.' )
		{
			if ( isdigit( *pCurr ) )
			{
				dwVer *= 10;
				dwVer += *pCurr - '0';
			
				*pCurr++;
			}
			else
			{
				return 0UL;
			}
		}

		*pCurr++;

		// Get the Minor number.
		while ( ( *pCurr != '.' ) && ( *pCurr != 0 ) )
		{
			if ( isdigit( *pCurr ) )
			{
				dwVer *= 10;
				dwVer += *pCurr - '0';
				iNumDigits++;
				
				*pCurr++;
			}
			else
			{
				break;
				//return 0UL;
			}
			
		}
		
		// Add padding zeros.
		if ( iNumDigits == 1 )
			dwVer *= 100;
		else if ( iNumDigits == 2 )
			dwVer *= 10;

		g_dwGL_VersionNumber = dwVer;

		return (float)dwVer/(float)1000;
	}
	
	return 0UL;

} // glfwGetOpenGLVersion()

void glEnable3D(float FOVy, float zNear, float zFar, int width, int height )
{
	GLint iViewport[4];

	// Get a copy of the viewport
	glGetIntegerv( GL_VIEWPORT, iViewport );

	// Save a copy of the projection matrix so that we can restore it 
	// when it's time to do 3D rendering again.
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();

	// Set up the orthographic projection
	gluPerspective((GLdouble)FOVy,(GLfloat)width/(GLfloat)height,(GLdouble)zNear, (GLdouble)zFar);
	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadIdentity();

	// Make sure depth testing and lighting are disabled for 2D rendering until
	// we are finished rendering in 2D
	glPushAttrib( GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT );
	//glEnable( GL_DEPTH_TEST );
	//glEnable( GL_LIGHTING );
	//glDisable(GL_DEPTH_TEST);
}

void glDisable3D( void )
{
	glPopAttrib();
	glMatrixMode( GL_PROJECTION );
	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
	glPopMatrix();
}

int InitGL(const int w, const int h)										// All Setup For OpenGL Goes Here
{
	#define GL_TEXTURE_RECTANGLE_ARB 0x84F5

	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);				// Background
	glClearDepth(1.0f);									// Depth Buffer Setup
	glClearStencil(0);									// Clear The Stencil Buffer To 0
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations
	
	glEnable(GL_TEXTURE_RECTANGLE_ARB);							// Enable 2D Texture Mapping
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	glGenTextures(1, &currentWindow->texture);
	glBindTexture(GL_TEXTURE_RECTANGLE_ARB, currentWindow->texture);
	glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho(0, w, h, 0, 0, 1);
	glMatrixMode (GL_MODELVIEW);

	return TRUE;										// Initialization Went OK
}

GLvoid ReSizeGLScene(GLsizei width, GLsizei height)
{
	if (height==0)			// Prevent A Divide By Zero By
		height=1;			// Making Height Equal One

	glViewport(0,0,width,height);	// Reset The Current Viewport

	glMatrixMode(GL_PROJECTION);	// Select The Projection Matrix
	glLoadIdentity();		// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	//gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);
	glOrtho(0, width, height, 0, 0, 1);

	glMatrixMode(GL_MODELVIEW);	// Select The Modelview Matrix
	glLoadIdentity();		// Reset The Modelview Matrix
}

void DrawString(HWND windowHandle, HDC hdc, HFONT font, GLuint font_base, const char* text, const char* color, const int x, const int y)
{
	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);

	glRasterPos2f(x,y);   

	glPushAttrib(GL_LIST_BIT);                          // Pushes The Display List Bits   
	glListBase(font_base - 32);                        // Sets The font base Character to 32   
	glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);  // Draws The Display List Text   
	glPopAttrib();// Pops The Display List Bits   

}

void drawLine(const float x1, const float y1, const float x2, const float y2)
{
	glPushMatrix();

	//glLoadIdentity();
	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);

	glBegin(GL_LINES);
	glVertex2f(x1, y1); 
	glVertex2f(x2, y2);
	glEnd();
	glPopMatrix();
}

void drawRectangle(const float x1, const float y1, const float x2, const float y2)
{
	//glLoadIdentity();
	glPushMatrix();

	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);

	glBegin(GL_QUADS);
	glVertex2f(x1, y1); 
	glVertex2f(x2, y1);
	glVertex2f(x2, y2);
	glVertex2f(x1, y2);
	glEnd();
	
	glPopMatrix();
}

void drawCircle(const int xCenter, const int yCenter, const int r)
{	
	int x=0,y=r;
	int d=3-(2*r);

	//glLoadIdentity();
	glPushMatrix();
	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);	

	glBegin(GL_POINTS);
	while(x<=y){

		glVertex2i(xCenter+x,yCenter+y);
		glVertex2i(xCenter+y,yCenter+x);
		glVertex2i(xCenter-x,yCenter+y);
		glVertex2i(xCenter+y,yCenter-x);
		glVertex2i(xCenter-x,yCenter-y);
		glVertex2i(xCenter-y,yCenter-x);
		glVertex2i(xCenter+x,yCenter-y);
		glVertex2i(xCenter-y,yCenter+x);

		if (d<0)
			d += (4*x)+6;
		else
		{
			d += (4*(x-y))+10;
			y -= 1;
		}
		x++;
	}
	glEnd();

	glPopMatrix();
}


void drawPolygon(const float* x, const float* y, const int length)
{
	int i=0;

	//glLoadIdentity();
	glPushMatrix();

	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);

	// Example of polygon
	glBegin(GL_POLYGON);

	for (; i!=length; i++)
		glVertex2f(x[i], y[i]);

	glEnd();
	glPopMatrix();
}

void drawArrow(const float x1, const float y1, const float x2, const float y2, const float size)
{
	float x1l=x1;
	float x2l=x2;
	float y1l=y1;
	float y2l=y2;

	float angle = atan2( (float) y1l - y2l, (float) x1l - x2l );
	float hypotenuse = sqrt((float)( (y1l- y2l)*(y1l - y2l) + (x1l - x2l)*(x1l - x2l) ));

	/* Here we lengthen the arrow by a factor of three. */
	x2l = (float) (x1l -  hypotenuse * cos(angle));
	y2l = (float) (y1l - hypotenuse * sin(angle));

	drawLine(x1l, y1l, x2l, y2l);
	//cvLine( frameToDraw, p, q, line_color, line_thickness, CV_AA, 0 );
	/* Now draw the tips of the arrow.  I do some scaling so that the
	* tips look proportional to the main line of the arrow.
	*/			
	x1l = (float) (x2l + size * cos(angle + PI / 4));
	y1l = (float) (y2l + size * sin(angle + PI / 4));
	drawLine(x1l, y1l, x2l, y2l);
	x1l= (float) (x2l + size * cos(angle - PI / 4));
	y1l = (float) (y2l + size * sin(angle - PI / 4));
	drawLine(x1l, y1l, x2l, y2l);
}

void drawPoints(const float* x, const float* y, const int length)
{
	int i=0;

	//glLoadIdentity();
	glPushMatrix();

	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);

	// Example of polygon
	glBegin(GL_POINTS);

	for (; i!=length; i++)
		glVertex2f(x[i], y[i]);

	glEnd();

	glPopMatrix();
}

// has to be in 3D mode
void drawPoints_3d(const float* x, const float* y, const float* z, const int length)
{
	int i=0;

	//glLoadIdentity();
	//glPushAttrib (GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	//glLoadIdentity();

	//glMatrixMode(GL_MODELVIEW);
	if (currentWindow->sx==0)
		currentWindow->sx=1;
	if (currentWindow->sy==0)
		currentWindow->sy=1;
	if (currentWindow->sz==0)
		currentWindow->sz=1;
	glScalef(currentWindow->sx, currentWindow->sy, currentWindow->sz);
	
	glRotatef(currentWindow->angle, currentWindow->rx, currentWindow->ry, currentWindow->rz);
	//glRotatef(currentWindow->gangle, currentWindow->grx, currentWindow->gry, currentWindow->grz);

	/*glRotatef(currentWindow->rx, 1,0,0);
	glRotatef(currentWindow->ry, 0,1,0);
	glRotatef(currentWindow->rz, 0,0,1);*/
	
	//glTranslatef(currentWindow->tx, currentWindow->ty, currentWindow->tz);

	glColor4f(currentWindow->color.red,currentWindow->color.green,currentWindow->color.blue,1);
	glLineWidth(currentWindow->line_width);

	// Example of polygon
	glBegin(GL_POINTS);

	for (; i!=length; i++)
		glVertex3f(x[i], y[i], z[i]);

	glEnd();

	glPopMatrix();
	//glPopAttrib();
}

int draw_image(const int width, const int height, const unsigned char* data)
{
	glPushMatrix();
	glDisable(GL_DEPTH_TEST);
	glEnable( GL_TEXTURE_RECTANGLE_ARB );
	glTexImage2D( GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGB, width, height, 0, GL_BGR, GL_UNSIGNED_BYTE, data);

	if (currentWindow->ImageDrawMode=DrawModeSTRETCH)
	{
		glBegin( GL_QUADS );
		glTexCoord2d(0.0,0.0); glVertex2d(0.0,0.0);
		glTexCoord2d(width,0.0); glVertex2d(currentWindow->width,0.0);
		glTexCoord2d(width,height); glVertex2d(currentWindow->width,currentWindow->height);
		glTexCoord2d(0.0,height); glVertex2d(0.0,currentWindow->height);
		glEnd();
	}
	else if (currentWindow->ImageDrawMode=DrawModeNORMAL)
	{
		glBegin( GL_QUADS );
		glTexCoord2d(0.0,0.0); glVertex2d(0.0,0.0);
		glTexCoord2d(width,0.0); glVertex2d(width,0.0);
		glTexCoord2d(width,height); glVertex2d(width,height);
		glTexCoord2d(0.0,height); glVertex2d(0.0,height);
		glEnd();
	}

	glDisable(GL_TEXTURE_RECTANGLE_ARB);
	glEnable(GL_DEPTH_TEST);

	glPopMatrix();

	return 0;
}

void disp_image(TWindowGL* window, const int width, const int height, const unsigned char* data)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	draw_image(width, height, data);

	glEndList();
	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_points(TWindowGL* window, float* x, float* y, const int length)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawPoints(x,y,length);

	glEndList();
	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_points_3d(TWindowGL* window, float* x, float* y, float* z, const int length)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawPoints_3d(x,y,z,length);

	glEndList();
	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_line(TWindowGL* window, const float x1, const float y1, const float x2, const float y2)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	//glColor3f(1.0,0.0,0.0);
	drawLine(x1,y1,x2,y2);
	/*glBegin(GL_TRIANGLES);								// Drawing Using Triangles
		glVertex3f( 0.0f, 1.0f, 0.0f);					// Top
		glVertex3f(-1.0f,-1.0f, 0.0f);					// Bottom Left
		glVertex3f( 1.0f,-1.0f, 0.0f);					// Bottom Right
	glEnd();*/

	glEndList();

	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_rectangle(TWindowGL* window, const float x1, const float y1, const float x2, const float y2)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawRectangle(x1,y1,x2,y2);

	glEndList();

	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_circle(TWindowGL* window, const float x, const float y, const float radius)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawCircle(x,y,radius);

	glEndList();

	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_polygon(TWindowGL* window, float* x, float* y, const int length)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawPolygon(x,y,length);

	glEndList();

	list_add(&window->dispList,index);
	window->displaySize++;
}

void disp_arrow(TWindowGL* window, float x1, float y1, float x2, float y2, float size)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	drawArrow(x1,y1,x2,y2,size);
	//drawLine(0.1,0.1,0.5,0.5);

	glEndList();

	list_add(&window->dispList,index);
	window->displaySize++;
}

void t_write_string_gl_t(TWindowGL* window, const char* string, const int x, const int y)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	DrawString(window->hWnd, window->hDC, window->font, window->_font_base, string, "red", x, y);

	glEndList();
	list_add(&window->dispList,index);
	window->displaySize++;
}

void draw_custom_gl_scene(TWindowGL* window, TGLOpCallback Callback, void* Data)
{
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);

	Callback(Data);

	glEndList();
	list_add(&window->dispList,index);
	window->displaySize++;
}

void register_custom_gl_scene(TWindowGL* window, TGLOpCallback Callback, void* Data)
{
	window->PaintCallback=Callback;
	window->preservedData = Data;
}

int DrawGLSceneDefault(void* Data)		// Here's Where We Do All The Drawing
{
	DisplayNode* temp=currentWindow->dispList;
	int i=currentWindow->displaySize-1;
	// Clear Screen And Depth Buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity ();
	
	if (i>=0)
	{
		while (i>=0)
		{
			/*glPushAttrib (GL_ALL_ATTRIB_BITS);
			glPushMatrix();
			//glLoadIdentity();

			//glMatrixMode(GL_MODELVIEW);
			glRotatef(currentWindow->rx, 1,0,0);
			glRotatef(currentWindow->ry, 0,1,0);
			glRotatef(currentWindow->rz, 0,0,1);
			if (currentWindow->sx==0)
				currentWindow->sx=1;
			if (currentWindow->sy==0)
				currentWindow->sy=1;
			if (currentWindow->sz==0)
				currentWindow->sz=1;

			glScalef(currentWindow->sx, currentWindow->sy, currentWindow->sz);
			glTranslatef(currentWindow->tx, currentWindow->ty, currentWindow->tz);
*/
			glCallList(temp->data); 
/*
			glPopMatrix();
			glPopAttrib();*/

			temp=temp->next;
			i--;
		}
	}

	return TRUE;			// Keep Going
}

int DrawGLScene(GLvoid)
{
	DrawGLSceneDefault(0);
	return 0;
}

 GLvoid KillGLWindow(TWindowGL* window)								// Properly Kill The Window
{
	WCHAR windowClass[80];

	if (window->hRC)											// Do We Have A Rendering Context?
	{
		if (!wglMakeCurrent(NULL,NULL))					// Are We Able To Release The DC And RC Contexts?
		{
			MessageBox(NULL,(LPCWSTR)("Release Of DC And RC Failed."),(LPCWSTR)("SHUTDOWN ERROR"),MB_OK | MB_ICONINFORMATION);
		}

		if (!wglDeleteContext(window->hRC))						// Are We Able To Delete The RC?
		{
			MessageBox(NULL,(LPCWSTR)("Release Rendering Context Failed."),(LPCWSTR)("SHUTDOWN ERROR"),MB_OK | MB_ICONINFORMATION);
		}
		window->hRC=NULL;										// Set RC To NULL
	}

	if (window->hDC && !ReleaseDC(window->hWnd,window->hDC))					// Are We Able To Release The DC
	{
		MessageBox(NULL,(LPCWSTR)("Release Device Context Failed."),(LPCWSTR)("SHUTDOWN ERROR"),MB_OK | MB_ICONINFORMATION);
		window->hDC=NULL;										// Set DC To NULL
	}

	if (window->hWnd && !DestroyWindow(window->hWnd))					// Are We Able To Destroy The Window?
	{
		MessageBox(NULL,(LPCWSTR)("Could Not Release hWnd."),(LPCWSTR)("SHUTDOWN ERROR"),MB_OK | MB_ICONINFORMATION);
		window->hWnd=NULL;										// Set hWnd To NULL
	}

	MultiByteToWideChar(CP_ACP,0,window->title, -1,windowClass,80); 
	if (!UnregisterClass((LPCWSTR)(windowClass),window->hInstance))			// Are We Able To Unregister Class
	{
		MessageBox(NULL,(LPCWSTR)TEXT("Could Not Unregister Class."),(LPCWSTR)TEXT("SHUTDOWN ERROR"),MB_OK | MB_ICONINFORMATION);
		window->hInstance=NULL;									// Set hInstance To NULL
	}
}

BOOL MoveGLWindow(TWindowGL* window, const int x, const int y)
{
	MoveWindow(window->hWnd, x,y, window->width,window->height, 1);
	return TRUE;
}

// bits ignored. default: 24
BOOL CreateGLWindow(TWindowGL* window, char* title, int X, int Y, int width, int height, const int bits)
{
	//T_LAST_OPNAME();
	GLuint		PixelFormat;							// Holds The Results After Searching For A Match
	WNDCLASS	wc;										// Windows Class Structure
	DWORD		dwExStyle;								// Window Extended Style
	DWORD		dwStyle;								// Window Style
	const char* t=title;
	WCHAR pTitle[80];
	WCHAR ErrorText[80];
	WCHAR windowErr[16];
	WCHAR fontName[16];
	HFONT   oldfont;
	int m_iYOffset, m_iXOffset;
	float version;

	static	PIXELFORMATDESCRIPTOR pfd=					// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),					// Size Of This Pixel Format Descriptor
		1,												// Version Number
		PFD_DRAW_TO_WINDOW |							// Format Must Support Window
		PFD_SUPPORT_OPENGL |							// Format Must Support OpenGL
	//	PFD_SUPPORT_GDI |
		PFD_DOUBLEBUFFER,								// Must Support Double Buffering
		PFD_TYPE_RGBA,									// Request An RGBA Format
		24,											// Select Our Color Depth
		0, 0, 0, 0, 0, 0,								// Color Bits Ignored
		0,												// No Alpha Buffer
		0,												// Shift Bit Ignored
		0,												// No Accumulation Buffer
		0, 0, 0, 0,										// Accumulation Bits Ignored
		16,												// 16Bit Z-Buffer (Depth Buffer)  
		1,												// Use Stencil Buffer ( * Important * )
		0,												// No Auxiliary Buffer
		PFD_MAIN_PLANE,									// Main Drawing Layer
		0,												// Reserved
		0, 0, 0											// Layer Masks Ignored
	};

	if(strlen(title)>=80)
	{
		printf( "Window Name should be less than 80 characters");
		return -1;
	}

	MultiByteToWideChar(CP_ACP,0,"TWindow Error",-1,windowErr,80); 
	MultiByteToWideChar(CP_ACP,0,t,-1,pTitle,80); 
	MultiByteToWideChar(CP_ACP,0,"Courier New",-1,fontName,16); 

	window->width=width;
	window->height=height;
	window->title=title;
	window->dispList=(DisplayNode*)(malloc(sizeof(DisplayNode)));
	window->line_width=1;
	window->dispList->next=NULL;
	window->displaySize=0;
	window->ImageDrawMode=DrawModeNORMAL;

	window->_font_base = glGenLists(100);
	window->font=(HFONT)GetStockObject(SYSTEM_FIXED_FONT);

	window->hInstance			= GetModuleHandle(NULL);		// Grab An Instance For Our Window
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc		= (WNDPROC) WindowProc;			// WndProc Handles Messages
	wc.cbClsExtra		= 0;							// No Extra Window Data
	wc.cbWndExtra		= 0;							// No Extra Window Data
	wc.hInstance		= window->hInstance;					// Set The Instance
	wc.hIcon			= LoadIcon(NULL, IDI_WINLOGO);	// Load The Default Icon
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);	// Load The Arrow Pointer
	wc.hbrBackground	= NULL;							// No Background Required For GL
	wc.lpszMenuName		= NULL;							// We Don't Want A Menu
	wc.lpszClassName	= (LPCWSTR)pTitle;						// Set The Class Name

	if (!RegisterClass(&wc))							// Attempt To Register The Window Class
	{
		//WCHAR errorStr[80];
		//MultiByteToWideChar(CP_ACP,0,"Failed To Register The Window Class.",-1,errorStr,80); 
		//MessageBox(NULL,errorStr,windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("Failed To register the window class");
		return -1;
	}

	dwExStyle=  WS_EX_TOOLWINDOW | WS_EX_APPWINDOW ;//WS_EX_TOOLWINDOW;//WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;//0;//WS_EX_APPWINDOW | WS_EX_WINDOWEDGE; //;WS_EX_TOOLWINDOW;//WS_EX_APPWINDOW ;//| WS_EX_WINDOWEDGE;	// Window Extended Style
	dwStyle= WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS;//| WS_CHILD | WS_OVERLAPPED;//;WS_OVERLAPPED | WS_SIZEBOX | WS_CLIPCHILDREN | WS_CLIPSIBLINGS;// WS_CAPTION;//WS_VISIBLE  ;//| WS_CLIPCHILDREN | WS_CLIPSIBLINGS ;	// Windows Style

	//dwExStyle=WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
	//dwStyle=WS_OVERLAPPEDWINDOW;							// Windows Style


	// Create The Window

	m_iYOffset = GetSystemMetrics(SM_CYMENU) + GetSystemMetrics(SM_CYBORDER)*2+1;
	m_iXOffset =GetSystemMetrics(SM_CXBORDER)*2+1;
	window->yOffset=m_iYOffset;
	window->xOffset=m_iXOffset;

	window->contentWidth=width;
	window->contentHeight=height;

	if (!(window->hWnd=CreateWindowEx(	dwExStyle,				// Extended Style For The Window
		(LPCWSTR)(pTitle),				// Class Name
		pTitle,//(LPCWSTR)(_T(window->title)),					// Window Title
		dwStyle,				// Window Style
		X,Y,					// Window Position
		width+window->xOffset+GetSystemMetrics(SM_CXFIXEDFRAME), height+m_iYOffset+GetSystemMetrics(SM_CYFIXEDFRAME),//+GetSystemMetrics(SM_CYSIZEFRAME),			// Selected Width And Height
		NULL,					// No Parent Window
		NULL,					// No Menu
		window->hInstance,				// Instance
		NULL)))					// Dont Pass Anything To WM_CREATE
	{
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,(LPCWSTR)("Window Creation Error."),windowErr,MB_OK|MB_ICONEXCLAMATION);
		//printf("Window creation error");
		printf("Window creation error");
		return -1;
		//return FALSE;									// Return FALSE
	}

	if (!(window->hDC=GetDC(window->hWnd)))								// Did We Get A Device Context?
	{
		//MultiByteToWideChar(CP_ACP,0,"Can't Activate The GL Rendering Context.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("Window DC cannot be acquired");
		return -1;
		//return FALSE;									// Return FALSE
	}

	if (!(PixelFormat=ChoosePixelFormat(window->hDC,&pfd)))		// Did Windows Find A Matching Pixel Format?
	{
		//MultiByteToWideChar(CP_ACP,0,"Can't Find A Suitable Pixel Format.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("Can't Find A Suitable Pixel Format. Check that 24 bit GL pixel format is supported by the system");
		return -1;
		//return FALSE;									// Return FALSE
	}

	if(!SetPixelFormat(window->hDC,PixelFormat,&pfd))			// Are We Able To Set The Pixel Format?
	{
		//MultiByteToWideChar(CP_ACP,0,"Can't Set The PixelFormat.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("Can't Set The PixelFormat");
		return -1;
		//return FALSE;									// Return FALSE
	}

	if (!(window->hRC=wglCreateContext(window->hDC)))					// Are We Able To Get A Rendering Context?
	{
		//MultiByteToWideChar(CP_ACP,0,"Can't Create A GL Rendering Context.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("PixelFormat cannot be set");
		return -1;
		//return FALSE;									// Return FALSE
	}

	if(!wglMakeCurrent(window->hDC, window->hRC))						// Try To Activate The Rendering Context
	{
		//MultiByteToWideChar(CP_ACP,0,"Can't Activate The GL Rendering Context.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		//MessageBox(NULL,(LPCWSTR)("Can't Activate The GL Rendering Context."),windowErr,MB_OK|MB_ICONEXCLAMATION);
		printf("GL rendering context couldn't be activated");
		return -1;
	}

	version=GetOpenGLVersion();
	if (version<2)
	{
		MultiByteToWideChar(CP_ACP,0,"TWindowGL at least requires OpenGL 2.0 to function properly.",-1,ErrorText,80); 
		KillGLWindow(window);
		MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		return -1;
	}

	ShowWindow(window->hWnd,SW_SHOW);							// Show The Window
	SetForegroundWindow(window->hWnd);							// Slightly Higher Priority
	SetFocus(window->hWnd);										// Sets Keyboard Focus To The Window

	ReSizeGLScene(width, height);						// Set Up Our Perspective GL Screen

	currentWindow=window;

	if (!InitGL(width,height))										// Initialize Our Newly Created GL Window
	{
		MultiByteToWideChar(CP_ACP,0,"Initialization Failed.",-1,ErrorText,80); 
		KillGLWindow(window);									// Reset The Display
		MessageBox(NULL,ErrorText,windowErr,MB_OK|MB_ICONEXCLAMATION);
		return -1;									// Return FALSE
	}

	oldfont = (HFONT)SelectObject(window->hDC, window->font);           // Selects The Font We Want   
	SelectObject(window->hDC, oldfont);                         // Selects The Font We Want   
	wglUseFontBitmaps(window->hDC, 32, 100, window->_font_base);         // Builds 96 Characters Starting At Character 32   
	DeleteObject(window->font); // Delete The Font 
	window->font=oldfont;

	window->active=TRUE;

	if (tWindows==NULL)
		tWindows=(TWindowGL**)(calloc(sizeof(TWindowGL*),1));//new TWindowGL*();

	tWindows[0]=currentWindow;
	currentWindow->preservedData=0;
	currentWindow->PaintCallback=DrawGLSceneDefault;

	// start trackball: Not thread safe
	//startTrackball (width/2, height/2, 0, 0, width, height);
	

	return 0;										// Success
}

static LRESULT CALLBACK WindowProc(	HWND hWnd,	// Handle For This Window
						 UINT	uMsg,	// Message For This Window
						 WPARAM	wParam,  // Additional Message Information
						 LPARAM	lParam)	// Additional Message Information
{
	POINT mousePos;
	RECT imagePos;

	if (currentWindow==0)
		return DefWindowProc(hWnd,uMsg,wParam,lParam);

	switch (uMsg)			// Check For Windows Messages
	{
	case WM_ACTIVATE:	// Watch For Window Activate Message
		{
			if (!HIWORD(wParam))// Check Minimization State
			{
				currentWindow->active=TRUE;// Program Is Active
			}
			else
			{
				currentWindow->active=FALSE;// Program Is No Longer Active
			}

			return 0;	// Return To The Message Loop
		}

	case WM_SYSCOMMAND:	// Intercept System Commands
		{
			switch (wParam)	// Check System Calls
			{
				// Screensaver Trying To Start?
			case SC_SCREENSAVE:	
				// Monitor Trying To Enter Powersave?
			case SC_MONITORPOWER:
				return 0;	// Prevent From Happening
			}
			break;			// Exit
		}

	case WM_CLOSE:	// Did We Receive A Close Message?
		{
			PostQuitMessage(0);// Send A Quit Message
			return 0;		// Jump Back
		}

	case WM_KEYDOWN:		// Is A Key Being Held Down?
		{
			currentWindow->keys[wParam] = TRUE;// If So, Mark It As TRUE
			return 0;		// Jump Back
		}

	case WM_MOUSEMOVE:		// Is mouse moving ?
		{
			//currentWindow->keys[wParam] = TRUE;// If So, Mark It As TRUE

			if (currentWindow->tracking && WM_LBUTTONDOWN)
			{
				POINT mousePos;
				RECT imagePos;

				mousePos.x = LOWORD(lParam); 
				mousePos.y = HIWORD(lParam);

				/*GetWindowRect(GetDlgItem(hwnd, IDC_BANNERSITE), &imagePos);
				if(PtInRect(&imagePos, mousePos))*/
				rollToTrackball((long)mousePos.x, (long)mousePos.y, currentWindow->rot); // rot is output rotation angle

				currentWindow->angle = currentWindow->rot[0];
				currentWindow->rx = currentWindow->rot[1];
				currentWindow->ry = currentWindow->rot[2];
				currentWindow->rz = currentWindow->rot[3];
			}
			
			return 0;		// Jump Back
		}

	case WM_LBUTTONDOWN:
		
		mousePos.x = LOWORD(lParam); 
		mousePos.y = HIWORD(lParam);

		currentWindow->rot[0]=0;
		currentWindow->rot[1]=0;
		currentWindow->rot[2]=0;
		currentWindow->rot[3]=0;
		startTrackball (mousePos.x, mousePos.y, 0, 0, currentWindow->width, currentWindow->height);
		currentWindow->tracking = 1;
		rollToTrackball((long)mousePos.x, (long)mousePos.y, currentWindow->rot);

		return 0;
	case WM_LBUTTONUP:
		mousePos.x = LOWORD(lParam); 
		mousePos.y = HIWORD(lParam);
		rollToTrackball((long)mousePos.x, (long)mousePos.y, currentWindow->rot); // rot is output rotation angle
		addToRotationTrackball(currentWindow->rot, currentWindow->globalRot);
		currentWindow->gangle = currentWindow->globalRot[0];
		currentWindow->grx = currentWindow->globalRot[1];
		currentWindow->gry = currentWindow->globalRot[2];
		currentWindow->grz = currentWindow->globalRot[3];
		currentWindow->rot[0]=0;
		currentWindow->rot[1]=0;
		currentWindow->rot[2]=0;
		currentWindow->rot[3]=0;
		currentWindow->tracking = 0;
		return 0;
		
	case WM_KEYUP:		// Has A Key Been Released?
		{
			currentWindow->keys[wParam] = FALSE;// If So, Mark It As FALSE
			return 0;			// Jump Back
		}

	case WM_SIZE:		// Resize The OpenGL Window
		{
			// LoWord=Width, HiWord=Height
			ReSizeGLScene(LOWORD(lParam),HIWORD(lParam)); 
			return 0;			// Jump Back
		}
	case WM_MOVING:
		//SwapBuffers(currentWindow->hDC);
		break;
	 case WM_ENTERSIZEMOVE:
		/* ValidateRect(currentWindow->hWnd,NULL);
		 SwapBuffers(currentWindow->hDC);
		 ShowWindow(currentWindow->hWnd, SW_SHOW);*/

		 InvalidateRect(currentWindow->hWnd,NULL,1);
		 break;
	 case WM_PAINT:
		 //DrawGLScene();
		 SwapBuffers(currentWindow->hDC);
		 //ValidateRect(currentWindow->hWnd,NULL);
		 break;
	}

	//SwapBuffers(currentWindow->hDC);
	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

int update_window(TWindowGL* window)
{
	int is_processed = 0;
	MSG msg;

	if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))		// Is There A Message Waiting?
	{
		if (msg.message==WM_QUIT)					// Have We Received A Quit Message?
		{
			return 0;								// If So done=TRUE
		}
		else										// If Not, Deal With Window Messages
		{
			window->PaintCallback(window->preservedData);
			SwapBuffers(window->hDC);
			TranslateMessage(&msg);					// Translate The Message
			DispatchMessage(&msg);					// Dispatch The Message
		}
	}
	else											// If There Are No Messages
	{
		// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
		if (window->active)									// Program Active?
		{
			if (window->keys[VK_ESCAPE])					// Was Escape Pressed?
			{
				return 0;							// ESC Signalled A Quit
			}
			else									// Not Time To Quit, Update Screen
			{
				window->PaintCallback(window->preservedData);
				//DrawGLScene(currentWindow->preservedData);						// Draw The Scene
				SwapBuffers(window->hDC);			// Swap Buffers (Double Buffering)
			}
		}
	}

	return 0;
}

 int wait_window(TWindowGL* window)
{
	int is_processed = 0;
	MSG msg;

	while(1)										// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))		// Is There A Message Waiting?
		{
			if (msg.message==WM_QUIT)					// Have We Received A Quit Message?
			{
				break;								// If So done=TRUE
			}
			else										// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);					// Translate The Message
				DispatchMessage(&msg);					// Dispatch The Message
			}
		}
		else											// If There Are No Messages
		{
			// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			if (window->active)									// Program Active?
			{
				if (window->keys[VK_ESCAPE])					// Was Escape Pressed?
				{
					break;							// ESC Signalled A Quit
				}
				else									// Not Time To Quit, Update Screen
				{					
					window->PaintCallback(window->preservedData);
					//DrawGLScene();						// Draw The Scene
					SwapBuffers(window->hDC);			// Swap Buffers (Double Buffering)
				}
			}
		}
	}

	return 0;

}

int wait_window_ms(TWindowGL* window, int miliSecs)
{
	int is_processed = 0;
	MSG msg;
	clock_t t1=time(0);
	clock_t t2=t1;

	while((double)(t2-t1)*1000.0 < miliSecs)										// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))		// Is There A Message Waiting?
		{
			if (msg.message==WM_QUIT)					// Have We Received A Quit Message?
			{
				break;								// If So done=TRUE
			}
			else										// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);					// Translate The Message
				DispatchMessage(&msg);					// Dispatch The Message
			}
		}
		else											// If There Are No Messages
		{
			// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			if (window->active)									// Program Active?
			{
				if (window->keys[VK_ESCAPE])					// Was Escape Pressed?
				{
					break;							// ESC Signalled A Quit
				}
				else									// Not Time To Quit, Update Screen
				{					
					window->PaintCallback(window->preservedData);
					//DrawGLScene();						// Draw The Scene
					SwapBuffers(window->hDC);			// Swap Buffers (Double Buffering)
				}
			}
		}
	
		t2=time(0);

	}

	return 0;

}

void clear_window(TWindowGL* window)
{
	clear_list(&window->dispList);
	currentWindow->displaySize=0;
	currentWindow->dispList=0;
	glDeleteLists(0,1000);
}

int close_window(TWindowGL* window)
{
	clear_window(window);
	KillGLWindow(window);
	return 0;
}

void set_color(TWindowGL* window, unsigned char R, unsigned char G, unsigned char B)
{
	currentWindow->color.red=(float)R/255.0;
	currentWindow->color.green=(float)G/255.0;
	currentWindow->color.blue=(float)B/255.0;
}

void set_rotation(TWindowGL* window, float rx, float ry, float rz)
{
	currentWindow->rx = rx;
	currentWindow->ry = ry;
	currentWindow->rz = rz;
}

void set_translation(TWindowGL* window, float tx, float ty, float tz)
{
	currentWindow->tx = tx;
	currentWindow->ty = ty;
	currentWindow->tz = tz;
}

void set_scale(TWindowGL* window, float sx, float sy, float sz)
{
	currentWindow->sx = sx;
	currentWindow->sy = sy;
	currentWindow->sz = sz;
}

TColorGL search_binary_string(const char* colorString)
{
	int i=0;
	int j=RGBSIZE;

	for (;;)
	{
		int k = ((j-i)>>1)+i;
		const int res = strcmp(colorString,rgbName[k]);

		if (res<0)
			j = k;
		else
			if (res>0)
				i = k+1;
			else
			{
				return rgbValue[k];
			}

		if (i==k && j==k || i==RGBSIZE)
		{
			//Color not found
			TColorGL retColor=COLORblack;
			return retColor;
		}
	}
}

void set_color_string(TWindowGL* win, const char* color)
{
	TColorGL c=search_binary_string(color);

	win->color.red=(float)c.red;
	win->color.green=(float)c.green;
	win->color.blue=(float)c.blue;
}

void set_line_width(TWindowGL* window, int lineWidth)
{
	window->line_width=lineWidth;
}

void set_draw_mode(TWindowGL* window, int DrawMode)
{
	window->RegionDrawMode=DrawMode;
}

void set_window(TWindowGL* window)
{
	currentWindow=window;
}

// -FontName-Height-Width-Italic-Underlined-Strikeout-[Bold-][CharSet-]
void set_font(TWindowGL* window, const char* fontString)
{
	WCHAR fontName[16];
	GLYPHMETRICSFLOAT gmf[256];
	HFONT   oldfont;
	char* fontNamePtr;
	int italic=0, underline=0, strikeout=0;
	int heightFont=0, widthFont=0;
	int weight=0; //FW_BOLD
	int i=0;
	int exitLoop=0;
	char * pch;
	//char* fntStr=(char*)fontString;
	char fntStr[100];
	int result;
	strcpy (fntStr,fontString);
	
	pch = strtok(fntStr,"-");
	while (pch != NULL && !exitLoop)
	{
		switch (i)
		{
		case 0:
			if (strcmp(pch,"*"))
				MultiByteToWideChar(CP_ACP,0,pch,-1,fontName,16);
			else
				MultiByteToWideChar(CP_ACP,0,"SYSTEM",-1,fontName,16);
			break;
		case 1:
			if (strcmp(pch,"*"))
				heightFont=atoi(pch);
			break;
		case 2:
			if (strcmp(pch,"*"))
				widthFont=atoi(pch);
			break;
		case 3:
			if (strcmp(pch,"*"))
				italic=atoi(pch);
			break;
		case 4:
			if (strcmp(pch,"*"))
				underline=atoi(pch);
			break;
		case 5:
			if (strcmp(pch,"*"))
				strikeout=atoi(pch);
			break;
		case 6:
			if (strcmp(pch,"*"))
				if (atoi(pch))
					weight=FW_BOLD;
			break;
		default:
			exitLoop=1;
			break;
		}

		pch = strtok (NULL, "-");	
		i++;
	}
	
	DeleteObject(window->font);
	window->font= CreateFont(-heightFont,                            // Height Of Font   
		-widthFont,                              // Width Of Font   
		0,                              // Angle Of Escapement   
		0,                              // Orientation Angle   
		weight,                       // Font Weight   
		italic,                          // Italic   
		underline,                          // Underline   
		strikeout,                          // Strikeout   
		ANSI_CHARSET,                   // Character Set Identifier   
		OUT_TT_PRECIS,                  // Output Precision   
		CLIP_DEFAULT_PRECIS,            // Clipping Precision   
		ANTIALIASED_QUALITY,            // Output Quality   
		FF_DONTCARE|DEFAULT_PITCH,      // Family And Pitch   
		fontName);                 // Font Name   

	oldfont = (HFONT)SelectObject(window->hDC, window->font);           // Selects The Font We Want   
	
	window->_font_base = glGenLists(100);
	wglUseFontBitmaps(window->hDC, 32, 100, window->_font_base);         // Builds 96 Characters Starting At Character 32   
	/*wglUseFontOutlines(	window->hDC,				// Select The Current DC
				32,				// Starting Character
				100,				// Number Of Display Lists To Build
				window->_font_base,				// Starting Display Lists
				0.0f,				// Deviation From The True Outlines
				0.2f,				// Font Thickness In The Z Direction
				WGL_FONT_POLYGONS,		// Use Polygons, Not Lines
				gmf);				// Address Of Buffer To Recieve Data*/
	SelectObject(window->hDC, oldfont);                         // Selects The Font We Want   
	DeleteObject(window->font); // Delete The Font 
	window->font=oldfont;
}

void set_font_default(TWindowGL* window)
{
	HFONT   oldfont;
	GLYPHMETRICSFLOAT gmf[256];

	DeleteObject(window->font);
	window->font=(HFONT)GetStockObject(SYSTEM_FIXED_FONT);
	oldfont = (HFONT)SelectObject(window->hDC, window->font);           // Selects The Font We Want   
	window->_font_base = glGenLists(100);
	wglUseFontBitmaps(window->hDC, 32, 100, window->_font_base);         // Builds 96 Characters Starting At Character 32   
	/*wglUseFontOutlines(	window->hDC,				// Select The Current DC
				32,				// Starting Character
				100,				// Number Of Display Lists To Build
				window->_font_base,				// Starting Display Lists
				0.0f,				// Deviation From The True Outlines
				0.2f,				// Font Thickness In The Z Direction
				WGL_FONT_POLYGONS,		// Use Polygons, Not Lines
				gmf);				// Address Of Buffer To Recieve Data*/


	SelectObject(window->hDC, oldfont);                         // Selects The Font We Want   
	DeleteObject(window->font); // Delete The Font 
	window->font=oldfont;
}

void get_pixel_buffer(unsigned char* dst, const int x, const int y, const int w, const int h, const int ws)
{
	//glScalef(1,-1,1);
	DrawGLScene();
	glReadPixels(0,0,w,h, GL_BGR, GL_UNSIGNED_BYTE,dst);	
}
