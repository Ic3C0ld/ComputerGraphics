#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <ctime>

#include "rigid.h"

class Simulation
{
public:
	Simulation(double boxSize=100,int sphereCount = 20, int particleCount = 0,int springCount=0);



	void update(double dt);
	void draw();


	double sim_dt;
	double target_dt;
	std::vector<Rigid*> objects; // will always add 6 planes for a complete box
	


};






#endif