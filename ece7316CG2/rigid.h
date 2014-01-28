#ifndef RIGID_H
#define RIGID_H	

#include <vector>
#include "matrix.h"


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
	virtual void calcCollisionResponse(Rigid*) = 0;
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
	virtual void calcCollisionResponse(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);


};


class Sphere :public Rigid
{
public:
	Sphere(double size, double mass, double Pxyz[], double Vxyz[]);
	

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void calcCollisionResponse(Rigid*) ;
	virtual void applyCollisionResponse();
	virtual void update(double dt);

};


class Particle : public Rigid
{
public:
	Particle(int partCount,double totalMass,double Pxyz[],double Vxyz[],double Wxyz );

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void calcCollisionResponse(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);

};


class StringSystem : public Rigid
{
public:
	StringSystem();

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void calcCollisionResponse(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);

};








#endif