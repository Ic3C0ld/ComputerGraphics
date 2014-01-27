#include "visuals.h"
#include "globals.h"


float roty = 0.0;


Matrix R(4,4);
Matrix w(4, 1);

Matrix dR_dt(4, 4);

double particle[] = { 0, 1, 2, 3, 4, -1, -2, -3, -4,
					  0, 0, 0, 0, 0, 0, 0, 0, 0,
					  0, 0, 0, 0, 0, 0, 0, 0, 0,
					  1, 1, 1, 1, 1, 1, 1, 1, 1, };

Matrix particle9(4,9);
Matrix particles2(4,9);


double mat[100 * 100];

void testVariableSetup()
{
	setIdentity(R);

	double wMat[] = {0,0.1,0,0};
	w = Matrix(4,1,wMat);

	double dR_mat[] = {		  0,	-w.mat[2],		w.mat[1],		0,
					   w.mat[2],			0,		-w.mat[0],		0,
					  -w.mat[1],	 w.mat[0],			   0,		0,
							  0,			0,			   0,		1	};


	for (int i = 0; i < 100 * 100; i++)
	{
		mat[i] = 0;
	}
	dR_dt = Matrix(4, 4, dR_mat);
	particle9 = Matrix(4, 9, particle);
	particles2 = Matrix(4, 9, particle);

}


void Render()
{
	///////////////// Preparation /////////////////////
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	projectionReset();
	LookAt();
	lights();




	////////////// DRAWING /////////////////////////////


	//// TESTING   ///////////////////////////////////////////////
	

//	R = dR_dt*R;

	R = R+dR_dt*R;
	particles2 = R*particles2;
	
	
	for (int i = 0; i < 10; i++)
	{
		Matrix(100, 100,mat);
	}
	
	/*for (int i = 0; i < particles2.m_columns; i++)
	{
		glPushMatrix();
		glTranslatef(particles2.mat[0 * particles2.m_columns + i], particles2.mat[1 * particles2.m_columns + i], particles2.mat[2 * particles2.m_columns + i]);
		glutSolidSphere(1, 20, 20);
		glPopMatrix();

	}*/
	
	//glMultMatrixd(R.transpose().getMat());
    //	glutSolidTeapot(10);




	

	glutSwapBuffers();          
}

//-----------------------------------------------------------



void Idle()
{
	roty += 0.2;

	glutPostRedisplay();
}


void Setup()  // TOUCH IT !! 
{
	testVariableSetup();


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

