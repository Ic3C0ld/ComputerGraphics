#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

#include "rigid.h"


class Simulation
{
public:
	Simulation(int ObjectCount=20);



	void update();
	std::vector<Rigid*> objects; // will always add 4 planes
	


};






#endif