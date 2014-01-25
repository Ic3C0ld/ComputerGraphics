#include "controls.h"

///// GLOBAL VARS ONLY TO THIS FILE /////////////////
//Controls

static  bool isElapsedTime;
static double mouseLoc[2];
static  bool isMousePressed[4];
bool isKeyPressed[256];

//Camera
static Matrix eye, center, up, aux, forward, left; //// ALL VECTORS == Matrix(3,1) 
static double r, theta, phi;

static double xScroll, yScroll;

static bool FPSmode;










//// CONTROLS ////////////////////////////////////////////////////////////////////////////////
void initControlVars()
{
	for (int i = 0; i < 256; i++){ isKeyPressed[i] = false; }
	for (int i = 0; i < 2; i++){ isMousePressed[i] = false; }

	//Mouse loc to be Initialized in other way hopefully by glGet* or something
	mouseLoc[0] = 450;// what ever initWindowSize says divided by 2
	mouseLoc[1] = 450;

}
void keyboardDown(unsigned char key, int x, int y)
{
	isKeyPressed[key] = true;
	keyboardFunc(key);

}
void keyboardUp(unsigned char key, int x, int y)
{
	isKeyPressed[key] = false;
}
void keyboardFunc(unsigned char key)
{
	if (isKeyPressed['w']){ eye = eye + forward; }
	if (isKeyPressed['a']){ eye = eye + left; }
	if (isKeyPressed['s']){ eye =eye - forward; }
	if (isKeyPressed['d']){ eye = eye - left; }


}
void specialKeyFunc(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_HOME:
		 FPSmode = ! FPSmode;
		if (FPSmode == true)
		{
			glutFullScreen();
			glutSetCursor(GLUT_CURSOR_NONE);
		}
		else
		{
			glutReshapeWindow(900, 900);
			glutPositionWindow(100, 100);
			glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
		}
		break;


	default:
		break;
	}
}

//////////MOUSE//
void  mousePassiveMotion(int x, int y)
{
	if ( FPSmode == true)
	{

		int width = glutGet(GLUT_WINDOW_WIDTH);
		int height = glutGet(GLUT_WINDOW_HEIGHT);

		if (x >= width - 1)
		{
			x = 0 + 1;
			mouseLoc[0] = x;
			glutWarpPointer(x, y);

		}
		else if (x <= 0)
		{
			x = width - 2;
			mouseLoc[0] = x;
			glutWarpPointer(x, y);

		}
		if (y >= height - 1)
		{
			y = 0 + 1;
			mouseLoc[1] = y;
			glutWarpPointer(x, y);

		}
		else if (y <= 0)
		{
			y = height - 2;
			mouseLoc[1] = y;
			glutWarpPointer(x, y);

		}

		 phi -= (x - mouseLoc[0]) / 400;
		 theta += (y - mouseLoc[1]) / 400;
	}

	mouseLoc[0] = x;
	mouseLoc[1] = y;

}
void  mouseFunc(int button, int state, int x, int y)
{
	//LEFT BUTTON
	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
	{
		isMousePressed[0] = true;

	}
	if (state == GLUT_UP && button == GLUT_LEFT_BUTTON)
	{
		isMousePressed[0] = false;
	}
	//MIDDLE BUTTON
	if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
	{
		glutSetCursor(GLUT_CURSOR_NONE);
		isMousePressed[2] = true;

	}
	if (state == GLUT_UP && button == GLUT_MIDDLE_BUTTON)
	{
		glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
		isMousePressed[2] = false;
	}

}
void  mouseMotion(int x, int y)
{
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);

	if (x > width - 2)
	{
		x = 0 + 1;
		mouseLoc[0] = x;
		glutWarpPointer(x, y);

	}
	else if (x <1)
	{
		x = width - 2;
		mouseLoc[0] = x;
		glutWarpPointer(x, y);

	}
	if (y > height - 2)
	{
		y = 2;
		mouseLoc[1] = y;
		glutWarpPointer(x, y);
	}
	else if (y <2)
	{
		y = height - 2;
		mouseLoc[1] = y;
		glutWarpPointer(x, y);

	}

	//////TO BE ADDED
	//Mouse event handling based and isMousePressed && isKeyPressed
	if (isMousePressed[2] == true)
	{
		 theta += (y - mouseLoc[1]) / 400;
		 phi -= (x - mouseLoc[0]) / 400;
	}


	//HUD
	if (isMousePressed[0] == true)
	{
		 yScroll -= (y - mouseLoc[1]);
		if ( yScroll < 0)  yScroll = 0; //Dont scroll above the first page
		//camera.phi -= (x - mouseLoc[0]) / 400;
	}

	/////
	mouseLoc[0] = x;
	mouseLoc[1] = y;
}





//////////////// CAMERA ///////////////////////////
void initCameraVars(double eyex, double eyey, double eyez, double R, double THETA, double PHI)
{

	up = Matrix(3, 1);
	aux = Matrix(3, 1);
	forward = Matrix(3, 1);
	left = Matrix(3, 1);

	double eyeVector[] = { eyex, eyey, eyez };
	double centerVector[] = { 0, 0, 0 };

	eye = Matrix(3, 1, eyeVector);
	center = Matrix(3, 1, centerVector);

	r = eye.norm();

	theta = acos(-eye.mat[2] / r);				//// theta = acos(-eye.z / r);
	phi = atan2(-eye.mat[1], -eye.mat[0]);		//// phi = atan2(-eye.y, -eye.x);


	///HUD
	xScroll = 0;
	yScroll = 0;
}
void LookAt()
{

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//these actions consider the familiar xyz system , not opengl's, SEE END of function
	forward.mat[0] = r*sin(theta)*cos(phi);		//x
	forward.mat[1] = r*sin(theta)*sin(phi);   //y
	forward.mat[2] = r*cos(theta);			//z
	normalize(forward);

	up.mat[0] = r*sin(theta - M_PI_2)*cos(phi);		//x
	up.mat[1] = r*sin(theta - M_PI_2)*sin(phi);		//y
	up.mat[2] = r*cos(theta - M_PI_2);				//z
	normalize(up);

	center = eye + forward;
	left = cross(up, forward);
	normalize(left);


	//!!!!!!!!!!!!!!!!!!!//REAL WORLD Coordinate system change to openGL, swap xyz to y'z'x'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*gluLookAt(eye.y, eye.z, eye.x,
		center.y, center.z, center.x,
		up.y, up.z, up.x);*/

	gluLookAt(eye.mat[1], eye.mat[2], eye.mat[0],
		center.mat[1], center.mat[2], center.mat[0],
		up.mat[1], up.mat[2], up.mat[0]); 


}

