#ifndef VISUALS_H
#define VISUALS_H	

#include "matrix.h"
#include "quaternion.h"
#include "glut.h"
#include "controls.h"
#include "rigid.h"
#include "simulation.h"

//-------- Functions --------------------------------

void Render();
// The function responsible for drawing everything in the 
// OpenGL context associated to a window. 

void Resize(int w, int h);
// Handle the window size changes and define the world coordinate 
// system and projection type

void Setup();
// Set up the OpenGL state machine and create a light source

void Idle();

void lights();

void projectionReset();



void drawFloor();

#endif