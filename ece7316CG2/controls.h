#ifndef CONTROLS_H
#define CONTROLS_H

#include <iostream>
#include <string>
#include <sstream>
#include "simulation.h"
#include "rigid.h"
#include "matrix.h"
#include "globals.h"

class Camera;
class MenuSystem;

	 void initControlVars();

	 void keyboardDown(unsigned char key, int x, int y);
	 void keyboardUp(unsigned char key, int x, int y);
	 void keyboardFunc(unsigned char key);
	 void specialKeyFunc(int Key, int x, int y);


	 void mouseFunc(int button, int state, int x, int y);
	 void mouseMotion(int x, int y);
	 void mousePassiveMotion(int x, int y);


	 void initMenus();
	 void menuHandler(int id);


	 void old_draw_HUD();

	
	 void hud();
	 void hudMain();

	 void mainScreen();
	 void partSelection();
	 
	 void statsScreen();

	 //// just so that i can know if a menu entry is pressed
	 class MenuEntry
	 {
	 public:
		 MenuEntry(bool& flagToControl,std::string entryName = std::string(""), double posX = 0, double posY = 0, double width = 0, double height = 0);

		 double xMax, xMin,yMax,yMin;
		 bool clicked;
		 bool* flagPointer;

		 void onClickAction();
		 bool isEntryClicked(double x, double y);
		 void drawRect();
		 void drawText();
		 void draw();
		 std::string text;
		
	 };

	 void DrawText(std::string&,double size, double x, double y);
	 void drawAxis(double size, std::vector<double> &data, std::vector<double> &time, std::string &str);

	 class MenuSystem
	 {
	 public:
		 MenuSystem();
		 //hudMain
		 MenuEntry* go_part_selection_screen,*preRender,*restart,*exit;

		//hudpartSelection
		 MenuEntry *go_stats_screen,
			 *a_collision_no_speedchange, *a_collide_with_speedchange,
			 *b1_use_material_properties, *b2_use_particles_no_rotation, *b2_use_particles_with_rotation,
			 *b3_use_followUp_camera,
			 *b4_use_springs_no_speedchange, *b4_use_springs_with_speedchange;
		 int current_screen;
	
		 void hud();//this will be called only,it will do the rest with the functions below

		 void hudMain();//will show help info ,keyboard shortcuts etc
		 void mainScreen();//will show part_selection_button ,etc.
		 void partSelection();//will show the partselection screen 
		 void statsScreen();//will show the stats of collisions


	 };

/////////////////////////// CAMERA ///////////////////////
	 void initCameraVars(double eyex = -50, double eyey = -50, double eyez = 50, double r = 1, double theta = M_PI_2, double phi = M_PI);
	 void LookAt();

	 



#endif