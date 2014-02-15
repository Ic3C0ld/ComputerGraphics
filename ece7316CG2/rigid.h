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
void calcCollision(Plane* , Particle* ,// Matrix cp/*ContactPoint*/,
													 Matrix&  rp/*RelativePosFromCenterMass*/,
													// Matrix vp /*VelocityofCp*/,
													 Matrix& n /*collision normal */,
													 double vRel/*Relative Velocity on the normal*/);
void calcCollision(Plane*, SpringSystem*);


void calcCollision(Sphere* s1, Sphere* s2, Matrix normal12, double vrel);
void calcCollision(Sphere*, Particle*, Matrix& r1,
										Matrix& r2,
										Matrix& n,
										double vRel);
void calcCollision(Sphere*, SpringSystem*);

void calcCollision(Particle* p1, Particle* p2, Matrix& r1, Matrix& r2, Matrix& n, double  vRel);
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
	virtual double getKinetik()=0;
	virtual Matrix getMomentum()=0;

	
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
	virtual double getKinetik();
	virtual Matrix getMomentum();

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
	virtual double getKinetik();
	virtual Matrix getMomentum();
	//
	
	Matrix Jtotal;
	

};


class Particle : public Rigid
{
public:
	Particle(int partCount, double radius, double totalMass, double Pxyz[], double Vxyz[], double Wxyz[],double color3[]);//Homogeneous points and vectors

	///// ALL INFO here /////////////////////
	Matrix points,	I,I_1,	x_t,v_t,	w_t,R_t;
	Quaternion rotation;

	double m_radius, m_mass;
	
	Matrix currentPoints,cI,cI_1;

	std::vector<Matrix> J;
	std::vector<Matrix> rP;

	///// ALL INFO here /////////////////////



	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);
	virtual double getKinetik() ;
	virtual Matrix getMomentum() ;
};


class SpringSystem : public Rigid
{
public:
	SpringSystem(Plane* plane2beAttachedto,
				double xPercent, double yPercent,
				double k1, double k2,
				double mass1,double mass2,
				double radius1,double radius2,
				double x1t_4_1[], double x2t_4_1[],
				double v1t_4_1[],double v2t_4_1[],
				double color[]);

	///// ALL INFO here /////////////////////
	Matrix x0; //where the spring is attached

	Matrix	x1_t, x2_t; //// positions w.r.t global frame
	Matrix v0, v1_t, v2_t; //// velocities 
	Matrix a0, a1_t, a2_t; ////accelerations

	double m_k1, m_k2;
	double m_mass1, m_mass2;
	double m_radius1, m_radius2;
	
	double xPer, yPer;
	Plane* plane;

	Matrix J1,J2;

	///// ALL INFO here /////////////////////

	virtual void draw();

	virtual void checkCollision(Rigid*);
	virtual void applyCollisionResponse();
	virtual void update(double dt);
	virtual double getKinetik();
	virtual Matrix getMomentum();
};








#endif