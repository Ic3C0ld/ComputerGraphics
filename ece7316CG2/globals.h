#ifndef GLOBALS_H
#define GLOBALS_H

#include "rigid.h"
#include "controls.h"
#include "main.h"
#include "simulation.h"
#include "visuals.h"





extern int boxSize;
extern int sphereCount;
extern int particleCount;
extern int springCount;
extern double simTime ;


extern std::vector<double>data_sphere_2_wallCollisions;
extern std::vector<double>data_particle_2_wallCollisions;
extern std::vector<double>data_spring_2_wallCollisions;
extern std::vector<double>data_collisions_2_everyone_else;
extern std::vector<double>data_time;


///PartA 
extern bool pA_collide_no_speedchange ;
extern bool pA_collide_with_speedchange ;
extern bool pA_collisionstats ;

//PartB
extern bool pB1_use_material_properties ;
extern bool pB2_use_particles_no_rotation ;
extern bool pB2_use_particles_with_rotation ;

extern bool pB3_use_followUp_camera ;

extern bool pb4_use_springs_with_speedchange ;




#endif