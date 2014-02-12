#include "visuals.h"
#include "globals.h"

float roty = 0.0;

Simulation simulation(/*BoxSize*/20,/*Spheres*/0,/*Particles*/1, /*Springs*/0);
double targetdt=0.005;



int count = 5;
double mass = 10;
double radius = 1;
double Pxyz[] = { 0, 20, 0, 0 };
double Vxyz[] = { 3, 1, 5, 0 };
double w[] = { 0, 20, 0, 0 };
double color3[] = {1,0,1};
//Particle particle(count, radius, mass, Pxyz, Vxyz, w,color3);


void testVariableSetup()
{
	
	
}


void Render()
{
	///////////////// Preparation /////////////////////
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	projectionReset();
	LookAt();

	//draw necessities
	lights();
	drawFloor();




	////////////// DRAWING /////////////////////////////

	simulation.draw();

	//// TESTING   ///////////////////////////////////////////////
	
	


	/*particle.draw();
	particle.update(targetdt/10);*/



	glutSwapBuffers();          
}

//-----------------------------------------------------------

void Idle()
{
	roty += 0.002;
	simulation.update(targetdt);
	glutPostRedisplay();
}


void Setup()  // TOUCH IT !! 
{
	testVariableSetup();


	initCameraVars();
	initControlVars();

	glEnable(GL_BLEND);

	//Parameter handling
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);  //renders a fragment if its z value is less or equal of the stored value
	glClearDepth(1);

	// polygon rendering mode
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	

	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);


	// Black background
	glClearColor(0.05f, 0.05f, 0.05f, 1.0f);

}



void Resize(int w, int h)
{
	// define the visible area of the window ( in pixels )
	if (h == 0) h = 1;
	glViewport(0, 0, w, h);

	// Setup viewing volume

	

}

void projectionReset()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);

	gluPerspective(65.0, (float)w / (float)h, 1.0, 800.0);
}

void lights()
{
	//Set up light source
	GLfloat light_position[] = { 0.0, 30.0, 50.0, 0.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	GLfloat ambientLight[] = { 0.1, 0.1, 0.1, 1.0 };
	GLfloat diffuseLight[] = { 0.9, 0.9, 0.9, 1.0 };
	GLfloat specularLight[] = { 1.0, 1.0, 1.0, 1.0 };


	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);

	glEnable(GL_LIGHT0);


	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(1, 1, 1);
	glPushMatrix();
			glTranslatef(0.0, 30, 50.0);
		glutSolidSphere(1, 10, 10);
	glPopMatrix();
	glPopAttrib();

}


void drawFloor()
{
	double bx = -100;
	double bz = -100;
	double dx = 8;
	double dz = 8;

	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 0.5, 0.5);
	glPolygonMode(GL_FRONT, GL_LINE);

	for (int j = 0; j <30; j++)
	{
		glBegin(GL_TRIANGLE_STRIP);


		for (int i = 0; i < 30; i++)
		{

			glVertex3f(bx, 0, bz);
			glVertex3f(bx, 0, bz + dz);
			glVertex3f(bx + dx, 0, bz);
			glVertex3f(bx + dx, 0, bz + dz);

			bx += dx;
		}

		glEnd();

		bx = -100;
		bz += dz;

	}
	glPolygonMode(GL_FRONT, GL_FILL);

	glPopAttrib();
}