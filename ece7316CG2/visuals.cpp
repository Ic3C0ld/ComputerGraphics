#include "visuals.h"
#include "globals.h"


float rotx = 0.0;

void Render()
{
	///////////////// Preparation /////////////////////
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	projectionReset();
	LookAt();
	lights();




	////////////// DRAWING /////////////////////////////

	double n[] = { 0, 1, 0 };
	Matrix normal(3, 1,n );

	double p[] = { 0, 0, 0 };
	Matrix Point(3, 1, p);


	double A = normal.mat[0];
	double B = normal.mat[1];
	double C = normal.mat[2];

	Matrix temp(normal.traspose());



	Matrix temp2(normal.traspose()*Point);

	double D = (normal*Point).getMat()[0];




	//// TESTING   ///////////////////////////////////////////////
		
	//glColor3f(1, 1, 1);
	double axis[] = {1, 1, 1};
	double point[] = { 0, 0, 0 };
	//glTranslatef(0, 0, -100);
	
	

	double PlanePoints[] = {-0.5,	 -0.5,	0.5,	0.5,
							 0,		 0,		0,		0,
							-0.5,	 0.5,	0.5,	 -0.5,
							1,		1,		1,		1	};

	double refNormal[] = { 0, 1, 0 };

	Matrix RefPlane(4, 4, PlanePoints);
	Matrix RefNormal(3, 1, refNormal);

	//glScalef(20, 20, 20);



	//glMultMatrixd(rotTdegrees((rotx), axis, point).getMat());


	/*glBegin(GL_QUAD_STRIP);
	glVertex3f(-0.5,-0.5,0);
	glVertex3f(0.5,-0.5,0);
	glVertex3f(-0.5, 0.5, 0);
	glVertex3f(0.5, -0.5, 0);
	glEnd();*/

	
	glScalef(10, 10, 10);

	glBegin(GL_QUADS);
	glVertex3f(RefPlane.getColumn(0).mat[0], RefPlane.getColumn(0).mat[1], RefPlane.getColumn(0).mat[2] );
	glVertex3f(RefPlane.getColumn(1).mat[0], RefPlane.getColumn(1).mat[1], RefPlane.getColumn(1).mat[2]);
	glVertex3f(RefPlane.getColumn(2).mat[0], RefPlane.getColumn(2).mat[1], RefPlane.getColumn(2).mat[2]);
	glVertex3f(RefPlane.getColumn(3).mat[0], RefPlane.getColumn(3).mat[1], RefPlane.getColumn(3).mat[2]);
	glEnd();

	double axisX[] = { 1, 0, 0 };

	
	
	for (int i = 0; i < 100; i++)
	{
		//RefPlane = rot(M_PI/100 , axisX, point)*RefPlane;
		//RefPlane = translate(0, 2, 0)*rotX(180 / 100)*RefPlane;

		//RefPlane = translate(0, 2, 0)*RefPlane;

		

		glRotated(180 / 100, 1, 0, 0);
		glTranslated(0, 2, 0);



		glBegin(GL_QUADS);
		glVertex3f(RefPlane.mat[0], RefPlane.mat[4], RefPlane.mat[8]);
		glVertex3f(RefPlane.mat[1], RefPlane.mat[5], RefPlane.mat[9]);
		glVertex3f(RefPlane.mat[2], RefPlane.mat[6], RefPlane.mat[10]);
		glVertex3f(RefPlane.mat[3], RefPlane.mat[7], RefPlane.mat[11]);
		glEnd();


	}
	
	



	//glMultMatrixd(rotTdegrees((rotx), axis, point).getMat());




	//glutSolidSphere(10, 1000, 1000);
	//glutSolidTeapot(5.0);

	glutSwapBuffers();          
}

//-----------------------------------------------------------



void Idle()
{
	rotx += 0.2;

	glutPostRedisplay();
}


void Setup()  // TOUCH IT !! 
{
	initCameraVars();
	initControlVars();



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

	/*glEnable(GL_CULL_FACE);
	glFrontFace(GL_CW);*/


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

	GLfloat ambientLight[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat specularLight[] = { 1.0, 1.0, 1.0, 1.0 };


	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);

	glEnable(GL_LIGHT0);



	glPushMatrix();

	glTranslatef(0.0, 30, 50.0);
	glutSolidSphere(1, 3, 3);
	glPopMatrix();

}

