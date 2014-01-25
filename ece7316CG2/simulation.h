#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

#include "rigid.h"


class Simulation
{
public:
	Simulation(int ObjectCount=20);


	std::vector<Rigid> objects; // will always add 4 planes

	void ode(Rigid&, double deltaT); //this will only call each Rigid's "calc_ddt" and "setStates"


};






#endif