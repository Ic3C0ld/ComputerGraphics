#include "simulation.h"



//Class Simulation//


Simulation::Simulation(double boxSize, int spheres, int particles, int springs)
{
	sim_dt = target_dt = 0.0001; //for now
	
	
//// BOX: Create planes  ////////////////////////

//TOP green
	double colorG[] = { 0.15, 0.3,0.15};
	Plane *top = new Plane(boxSize,colorG);
	top->m_plane = translate(0, boxSize / 2, 0)*top->m_plane;

//FRONT dark grey
	double colorDGrey[] = { 0.3, 0.3, 0.3 };
	Plane *front = new Plane(boxSize, colorDGrey);
	front->m_plane = rotX(M_PI_2)* top->m_plane;
 
//BOTTOM white
	double colorR[] = { 1, 1, 1 };
	Plane *bottom = new Plane(boxSize,colorR);
	bottom->m_plane = rotX(M_PI_2)*front->m_plane;

//BACK dark grey
	Plane *back = new Plane(boxSize, colorDGrey);
	back->m_plane = rotX(M_PI_2)*bottom->m_plane;

//LEFT light grey
	double colorGrey[] = { 0.5, 0.5, 0.5 };
	Plane *left = new Plane(boxSize, colorGrey);
	left->m_plane = rotY(M_PI_2)*front->m_plane;

//Right light grey
	Plane *right = new Plane(boxSize, colorGrey);
	right->m_plane = rotY(-M_PI_2)*front->m_plane;


	objects.push_back(top);
	objects.push_back(bottom);
	objects.push_back(front);
	objects.push_back(back);
	objects.push_back(left);
	objects.push_back(right);

	//translate the whole box to sit on the ground, not to have 0,0,0 as center
	for (int i = 0; i < 6; i++)
	{
		objects[i]->p_plane->m_plane = translate(0, boxSize / 2, 0)*objects[i]->p_plane->m_plane;
	}

//////////////////////////// BOX READY //////////////




}



void Simulation::draw()
{
	for (int i = 0; i < objects.size(); i++)
	{
		objects[i]->draw();
	}
}