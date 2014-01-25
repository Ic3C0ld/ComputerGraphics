#ifndef CONTROLS_H
#define CONTROLS_H

#include "matrix.h"
#include "globals.h"

class Camera;

//void initControls;
//
//void keyboardDown(unsigned char key, int x, int y);
//void keyboardDown(unsigned char key, int x, int y);
//void keyboardUp(unsigned char key, int x, int y);
//void keyboardFunc(unsigned char key);
//void specialKeyFunc(int key, int x, int y);
//
////////////MOUSE//
//void mousePassiveMotion(int x, int y);
//void mouseFunc(int button, int state, int x, int y);
//void mouseMotion(int x, int y);

	
	 void initControlVars();

	 void keyboardDown(unsigned char key, int x, int y);
	 void keyboardUp(unsigned char key, int x, int y);
	 void keyboardFunc(unsigned char key);
	 void specialKeyFunc(int Key, int x, int y);


	 void mouseFunc(int button, int state, int x, int y);
	 void mouseMotion(int x, int y);
	 void mousePassiveMotion(int x, int y);


	////void menuAll(int choice, Container&);
	//void initMenus();
	//int idMenu1, idMenu2, idMenu3;


	



/////////////////////////// CAMERA ///////////////////////
	 void initCameraVars(double eyex = -50, double eyey = -50, double eyez = 50, double r = 1, double theta = M_PI_2, double phi = M_PI);
		 void LookAt();

	



#endif