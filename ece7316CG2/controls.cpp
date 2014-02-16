#include "controls.h"

///// GLOBAL VARS ONLY TO THIS FILE /////////////////
//Controls

Simulation* simulation;
extern int boxSize;
extern int sphereCount;
extern int particleCount;
extern int springCount;


static  bool isElapsedTime;
static double mouseLoc[2];
static  bool isMousePressed[4];
bool isKeyPressed[256];

//OTHER
MenuSystem graphicMenu;

//classic menu system
#define TOGGLE_STATS					00
#define B2B_COLLISIONS_NO_SPEED_CHANGE	10
#define B2B_COLLISIONS_SPEED_CHANGE		11
#define USE_MATERIAL_PROPERTIES			20
#define PARTICLE_COLLISION_NO_ROTATION  21
#define PARTICLE_COLLISION_ROTATION		22
#define SPRING_COLLISION				23
#define FOLLOW_UP_CAMERA				24
#define RESTART							25




int menu ,subMenuA, subMenuB;;


//Camera
static Matrix eye, center, up, aux, forward, left; //// ALL VECTORS == Matrix(3,1) 
static double r, theta, phi;

static double xScroll, yScroll;

static bool FPSmode;

//////// Simulation Info for HUD //////////////////////////////

#define DEFAULT_SCREEN			=00;
#define MAIN_SCREEN				=10;
#define PART_SELECTION_SCREEN	=20;
#define STATS_SCREEN			=30;

bool default_screen = true;
bool part_selection_screen=false;
bool main_screen = false;
bool stats_screen = false;

////
bool hudScreen[40]; //will be used for exampleas  hud[MAIN_SCREEN]==true or else;.//not all values will have meaning


//// collected stats ////
double totalKinetic;

extern int sphere_2_wallCollisions;
extern int particle_2_wallCollisions ;
extern int spring_2_wallCollisions ;

extern int collisions_2_everyone_else;

extern int followUpObjID ;//the first 0-5 elements are the planes and they dont have v_t,x_t so it will break the program

/////  bool flags that control the simulation properties ////////////////////
bool restartSimulation = false;
bool preRenderflag = true;
double preRenderSeconds = 30;

///PartA 
extern bool pA_collide_no_speedchange ;
extern bool pA_collide_with_speedchange ;
extern bool pA_collisionstats ;

//PartB
extern bool pB1_use_material_properties ;
extern bool pB2_use_particles_no_rotation ;
extern bool pB2_use_particles_with_rotation;

extern bool pB3_use_followUp_camera;

extern bool pb4_use_springs_with_speedchange ;




//// CONTROLS ////////////////////////////////////////////////////////////////////////////////
void initControlVars()
{
	for (int i = 0; i < 256; i++){ isKeyPressed[i] = false; }
	for (int i = 0; i < 2; i++){ isMousePressed[i] = false; }

	//Mouse loc to be Initialized in other way hopefully by glGet* or something
	mouseLoc[0] = 450;// what ever initWindowSize says divided by 2
	mouseLoc[1] = 450;

	totalKinetic = 0;


	for (int i = 0; i < 40; i++)
	{
		hudScreen[i] = false;
	}


	//menu = MenuSystem();





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
	if (FPSmode == true || isKeyPressed['x'] == true)
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
	if ((state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)  )
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
	if (isMousePressed[2] == true )
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

	if (pB3_use_followUp_camera == false)
	{
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

	else
	{
		Matrix x, v ;
		double distance;
		
		switch (simulation->objects[followUpObjID]->m_id)
		{
		case RIGID_SPHERE:
		{
			Sphere * s = static_cast<Sphere*>(simulation->objects[followUpObjID]);

			 v = s->V_t;
			
			 x = s->X_t;

			distance = 4 * s->m_radius;
		}
		default:
			break;
		}
		

		center = x;
		forward = v;		//z
		normalize(forward);

		up.mat[0] = 0;		//x
		up.mat[1] = 1;		//y
		up.mat[2] = 0;		//z
	

		eye = center - distance*forward;

		



		gluLookAt(eye.mat[0], eye.mat[1], eye.mat[2],
			center.mat[0], center.mat[1], center.mat[2],
			up.mat[0], up.mat[1], up.mat[2]);
	}


	
	


}












////// HUD STUFF  ////////////////////////////////////////////////////////

//class MenuSystem

MenuSystem::MenuSystem()
{
	go_part_selection_screen = new MenuEntry(part_selection_screen, "Part Selection", -20, +50, 20, 5);





}

////CLASS MENUENTRY

MenuEntry::MenuEntry(bool& flagToControl,std::string entryName, double posX, double posY, double width, double height)
{
	xMin = posX;
	yMin = posY;
	xMax = posX + width;
	yMax = posY + height;
	clicked = false;
	text = entryName;
	flagPointer = &flagToControl;
}

void MenuEntry::drawRect()
{
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(0.5, 0.5, 0.5);
	glLineWidth(4);

	glPushMatrix();

	glBegin(GL_LINE_STRIP);
		glVertex3f(xMin, yMin, 0);
		glVertex3f(xMax, yMin, 0);
		glVertex3f(xMax, yMax, 0);
		glVertex3f(xMin, yMax, 0);
		glVertex3f(xMin, yMin, 0);
	glEnd();
	
	glPopMatrix();
	glPopAttrib();

}

void MenuEntry::drawText()
{
	DrawText(text,0.002, xMin, yMin);

}

void MenuEntry::draw()
{
	drawRect();
	drawText();
}

bool MenuEntry::isEntryClicked(double x,double y)
{
	if (( xMin < x < xMax) && (yMin < y < yMax))
		return true;
	
	return false;	

}

void MenuEntry::onClickAction()
{
	(*flagPointer) = !(*flagPointer);//toggle the flag
	
}

//////////////////////// OLD STUFF ////////////////////
void DrawText(std::string& str, double size, double x, double y)
{
	glPushMatrix();
	
	glTranslatef(x, y, 0);
	glScalef(size, size, size);

	for (int i = 0; i<str.size(); i++)
		glutStrokeCharacter(GLUT_STROKE_ROMAN, str[i]);
	glPopMatrix();


}



void old_draw_HUD()
{
	glDisable(GL_DEPTH_TEST);//So that HUD doesnt disappear behind other objects
	glDepthMask(GL_FALSE); //Just because internet said so!  :P



	//SAVE MATRICES PROJECTION AND MODELVIEW and SET their state for HUD drawing
	glMatrixMode(GL_PROJECTION);//Save projection Matrix
	glPushMatrix();
	glLoadIdentity();

	//set up projection
	double orthoSize = 100;
	glOrtho(-orthoSize, orthoSize, -orthoSize, orthoSize, orthoSize, -orthoSize);


	glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(0, 0, -1);//draw everything slightly further away
	glTranslatef(0.002*orthoSize*xScroll, 0.002*orthoSize*yScroll, 0);

	/////////////////////	START drawing	///////////////////////////////////////////////////////

	glColor3f(1, 1, 1);

	//DRAW quick Tips
	std::stringstream quickTip(" Left Click & scrolldown : STATS");
	DrawText(quickTip.str(), 0.0002*orthoSize, -orthoSize, 0.96*orthoSize);

	std::stringstream quickTip1(" Right Click : MENU");
	DrawText(quickTip1.str(), 0.0002*orthoSize,  -orthoSize, 0.92 * orthoSize);

	std::stringstream quickTip2(" W,A,S,D && middle mouse or HOME key : MOVE");
	DrawText(quickTip2.str(), 0.0002*orthoSize, - orthoSize, 0.88* orthoSize);

	std::stringstream quickTip3(" For other objects plz adjust count at visuals.cpp");
	DrawText(quickTip3.str(), 0.0002*orthoSize, - orthoSize,0.84 * orthoSize);


	//DrawData
	glPushMatrix();
	glTranslatef(-0.95*orthoSize,-0.95*orthoSize,0);
	drawAxis(orthoSize, data_sphere_2_wallCollisions, data_time, std::string("ball2wall collisions"));
	glTranslatef(+0.95*orthoSize, 0, 0);
	drawAxis(orthoSize, data_collisions_2_everyone_else, data_time, std::string("any Other collisions"));

	glPopMatrix();


	////////////////////	END drawing 	///////////////////////////////////////////////////////
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();


	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
}



void drawAxis(double size, std::vector<double> &data, std::vector<double> &time, std::string &str)
{
	//Draw the x,y axi
	glPushMatrix();

	glLineWidth(2);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0, 0, 0);
	glVertex3f(0.75*size, 0, 0);
	glVertex3f(0.73*size, 0.005*size, 0);
	glVertex3f(0.73*size, -0.005*size, 0);
	glVertex3f(0.75*size, 0, 0);
	glEnd();


	glRotatef(90, 0, 0, 1);

	glLineWidth(2);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0, 0, 0);
	glVertex3f(0.55*size, 0, 0);
	glVertex3f(0.53*size, 0.005*size, 0);
	glVertex3f(0.53*size, -0.005*size, 0);
	glVertex3f(0.55*size, 0, 0);
	glEnd();

	glPopMatrix();


	//Draw the Data points
	float y=0,axisMax = 0,t=0,hits = 0, maxY = 0,currentData=0;
	float maxT = 0;
	glPushMatrix();

	glPointSize(2);
	glBegin(GL_POINTS);

	
	for (int i = 0; i < data.size(); i++)
	{
		 maxT = 5>data_time[data_time.size() - 1] ? 5 : data_time[data_time.size() - 1];
		if (i < data_time.size() - 1)
		{
			t = data_time[i] / maxT;
		}
		else
		{
			t = maxT;
		}

		 maxY = 100>data[data.size() - 1] ? 100 : data[data.size() - 1];
		 y = data[i] / maxY;

		glVertex3f(0.7*size*t, 0.5*size*y, 0);
	}
	glPointSize(1);

	glEnd();

	//Draw Data name
	std::stringstream name("");
	name << "#" << str << std::endl;
	glPushMatrix();
	DrawText(name.str(), 0.0002*size, 0.2*size, 0.6*size);
	glPopMatrix();

	//Draw the #HITS
	std::stringstream temp2("");
	temp2 <<"#"<< maxY <<std::endl;
	glPushMatrix();
	DrawText(temp2.str(), 0.0002*size, -0.05*size, 0.6*size);
	glPopMatrix();


	std::stringstream temp("");
	temp << maxT << " sec "<<std::endl;
	glPushMatrix();
	DrawText(temp.str(), 0.0002*size, 0.75*size, 0.02*size);
	glPopMatrix();


	glPopMatrix();



}




void initMenus()
{
	
	

	subMenuA = glutCreateMenu(menuHandler);
	glutAddMenuEntry("TOGGLE Sphere collisions", B2B_COLLISIONS_NO_SPEED_CHANGE);
	glutAddMenuEntry("TOGGLE Sphere ELASTIC collisions", B2B_COLLISIONS_SPEED_CHANGE);

	subMenuB = glutCreateMenu(menuHandler);
	glutAddMenuEntry("TOGGLE particle ", PARTICLE_COLLISION_NO_ROTATION);
	glutAddMenuEntry("TOGGLE particle collisions WITH rotation", PARTICLE_COLLISION_ROTATION);
	glutAddMenuEntry("TOGGLE springs", SPRING_COLLISION);
	


	 glutCreateMenu(menuHandler);
	glutAddSubMenu("Part A", subMenuA);
	glutAddSubMenu("Part B", subMenuB);
	glutAddMenuEntry("TOGGLE STATS", TOGGLE_STATS);
	glutAddMenuEntry("TOGGLE use Material properties", USE_MATERIAL_PROPERTIES);
	glutAddMenuEntry("TOGGLE follow up camera", FOLLOW_UP_CAMERA);
	glutAddMenuEntry("RESTART simulation", RESTART);

	glutAttachMenu(GLUT_RIGHT_BUTTON);



}

void menuHandler(int id)
{
	switch (id)
	{
	case TOGGLE_STATS:
	{
		pA_collisionstats = !pA_collisionstats;
		break;
	}
	case B2B_COLLISIONS_NO_SPEED_CHANGE:
	{
		pA_collide_no_speedchange = !pA_collide_no_speedchange;
		break;
	}
	case B2B_COLLISIONS_SPEED_CHANGE:
	{
		pA_collide_with_speedchange = !pA_collide_with_speedchange;
		break;
	}
	case USE_MATERIAL_PROPERTIES:
	{
		pB1_use_material_properties = !pB1_use_material_properties;
		break;
	}
	case PARTICLE_COLLISION_NO_ROTATION:
	{
		pB2_use_particles_no_rotation = !pB2_use_particles_no_rotation;
		break;
	}
	case PARTICLE_COLLISION_ROTATION:
	{
		pB2_use_particles_with_rotation = !pB2_use_particles_with_rotation;
		break;
	}
	case SPRING_COLLISION:
	{
		pb4_use_springs_with_speedchange = !pb4_use_springs_with_speedchange;
		break;
	}
	case FOLLOW_UP_CAMERA:
	{
		pB3_use_followUp_camera = !pB3_use_followUp_camera;
		break;
	}
	case RESTART:
	{
		delete(simulation);
		simulation = new  Simulation(boxSize, sphereCount, particleCount, springCount);
		break;
	}
	default:
		break;
	}
}



