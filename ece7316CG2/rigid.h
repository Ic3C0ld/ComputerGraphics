#ifndef RIGID_H
#define RIGID_H	

#include <vector>
#include "matrix.h"
#include "quaternion.h"

#define RIGID_BASE 0
#define RIGID_PLANE 1
#define RIGID_SPHERE 2
#define RIGID_PARTICLE 3
#define RIGID_SPRING 4

class Rigid;
class Sphere;
class Particle;
class Plane;
class SpringSystem;


void calcCollision(Plane*, Sphere*);
void calcCollision(Plane*, Particle*);
void calcCollision(Plane*, SpringSystem*);
void calcCollision(Plane*, Sphere*);

void calcCollision(Sphere*, Sphere*);
void calcCollision(Sphere*, Particle*);
void calcCollision(Sphere*, SpringSystem*);

void calcCollision(Particle*, SpringSystem*);



class Rigid
{
public:

	Rigid();



	////// State Variables //////////////////
	//
	//Each subClass will have its own, because of variety 
	//it doesnt make sense to define so many to cover them all
	//
	///// Unnecessary :/  ////////////////////
	bool m_justCollided;
	double m_color[3];

	////Hopefully these will make casts unnecessary
	////These pointers will either be NULL or pointer == *this if ->id == this->id
	int m_id;
	int getWhatAmI();				//// will return defines suppose Plane==1,Sphere==2 etc...

	Plane* p_plane;
	Sphere* p_sphere;
	Particle* p_particle;
	SpringSystem* p_spring;

	/////////////////////////////////////
			
	virtual void draw()=0;

	virtual void checkCollision(Rigid*)=0;
	virtual void applyCollisionResponse()=0;
	virtual void update(double dt)=0;

	
};


class Plane : public Rigid
{
public:
	Plane(double size=1);
	Plane(double size, double color[]);

	//			m_plane N==normal of plane and all rest are the plane's points minimum 3 
	////	[	nx	p0x p1x p2x	p3x	...	]
	////	[	ny	p0y	p2y	p2x	p3x	...	]
	////	[	nz	p0z	p2z	p2x	p3x	...	]
	////	[	0	1	1	1	1	...	]
	////
	Matrix m_plane;

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);


};


class Sphere :public Rigid
{
public:
	Sphere(double radius, double mass, double Pxyz[], double Vxyz[], double color3[]);
	
	double m_radius, m_mass;
	
	Matrix X_t,V_t;

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);

	//
	
	Matrix applyV_t;//this will be calculated in a collision,

};


class Particle : public Rigid
{
public:
	Particle(int partCount, double radius, double totalMass, double Pxyz[], double Vxyz[], double Wxyz[]);//Homogeneous points and vectors

	///// ALL INFO here /////////////////////
	Matrix points,	I,I_1,	x_t,v_t,	w_t,R_t;
	Quaternion rotation;

	double m_radius, m_mass;

	Matrix currentPoints;
	///// ALL INFO here /////////////////////



	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);

};


class StringSystem : public Rigid
{
public:
	StringSystem();

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);

};








#endif