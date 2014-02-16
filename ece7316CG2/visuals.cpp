#include "visuals.h"
#include "globals.h"

clock_t simTimeClocks = 0;
int boxSize = 30;
int sphereCount = 5;
int particleCount = 2;
int springCount =3;
extern Simulation* simulation;
double targetdt=0.005;

double colorG[] = { 1, 0.3, 0.15 };
Plane *top = new Plane(10, colorG);

double x1[] = { -40, -200, 0, 1 };
double x2[] = { 40, 40, 0, 1 };
double v[] = { 0, 0, 0, 0 };

SpringSystem test(top, 35, 65, 10, 10, 5, 4, 4, 1, x1, x2, v, v, colorG);





void testVariableSetup()
{
	
	
	top->m_plane = translate(0, 50, 0)*top->m_plane;
	test.x0 = translate(0, 50, 0)*test.x0;

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

	simulation->draw();


	old_draw_HUD();
	//// TESTING   ///////////////////////////////////////////////
	/*test.draw();
	top->draw();
*/
	glutSwapBuffers();          
}

//-----------------------------------------------------------

void Idle()
{ //baa, DEN katafera na to isorrophsw na kanei aneksarthta tou idleDT

	clock_t timeNow=clock();

	double deltaT = (timeNow - simTimeClocks) / CLOCKS_PER_SEC;
	if (deltaT > targetdt)
	{
		simulation->update(targetdt);
		simTimeClocks -= targetdt*CLOCKS_PER_SEC;
		
	}
	
	glutPostRedisplay();

	
}


void Setup()  // TOUCH IT !! 
{
	testVariableSetup();
	initMenus();
	simTimeClocks = clock();
	initCameraVars();
	initControlVars();
	simulation = new  Simulation(boxSize, sphereCount, particleCount, springCount);
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

	/*glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);*/

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