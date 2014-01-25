#ifndef RIGID_H
#define RIGID_H	

#include <vector>
#include "matrix.h"

class Rigid;
class Sphere;
class Plane;
class StringSystem;


class Rigid
{
public:

	Rigid();



	////// State Variables //////////////////
	Matrix  X_t ;		//3x1	position of center mass
	Matrix  R_t ;		//3x3	rotation matrix
	Matrix  P_t ;		//3x1	m*u
	Matrix  L_t ;		//3x1	I*w

	////// State Derivatives ///////////////
	Matrix  dX_dt;		//3x1	position
	Matrix  dR_dt;		//3x3	rotation matrix
	Matrix  dP_dt;		//3x1	m*u
	Matrix  dL_dt;		//3x1	torque

	///// Extra info //////////////////////
	Matrix I;			//3x3 Inertia
	Matrix Iinv ;		//3x3 Inverse Inertia
	Matrix v ;			//3x3 velocity
	Matrix w ;			//3x3 angular velocity

	bool justCollided;
	///// Unnecessary ////////////////////
	double color[3];
	int whatAmI;
	/////////////////////////////////////
	bool Questionflags[10]; ////to implement the questions in order




	virtual void calculate_forces();
	virtual void calculate_ddt();

	virtual void setState(std::vector<double> &);
	virtual std::vector<double> getState();

	virtual void draw();

	virtual void checkCollision(Rigid&);
	virtual void collide(Rigid&);
	virtual int getWhatAmI();				//// will return defines suppose Sphere==1,Particle==2 etc...


};




class Sphere :public Rigid
{
public:
	Sphere(double size, double mass, double Pxyz[], double Vxyz[]);
	virtual void calculate_forces();
	virtual void calculate_ddt();

	virtual void setState(std::vector<double> &);
	virtual std::vector<double> getState();

	virtual void draw();

	virtual void checkCollision(Rigid&);
	virtual void collide(Rigid&);
	virtual int getWhatAmI();				//// will return defines suppose Sphere==1,Particle==2 etc...


};


class Particle : public Rigid
{
public:
	Particle(int partCount,double totalMass,double Pxyz[],double Vxyz[],double Wxyz );

	virtual void calculate_forces();
	virtual void calculate_ddt();

	virtual void setState(std::vector<double> &);
	virtual std::vector<double> getState();

	virtual void draw();

	virtual void checkCollision(Rigid&);
	virtual void collide(Rigid&);
	virtual int getWhatAmI();				//// will return defines suppose Sphere==1,Particle==2 etc...


};

class Plane : public Rigid
{
public:
	Plane(double size);
	virtual void calculate_forces();
	virtual void calculate_ddt();

	virtual void setState(std::vector<double> &);
	virtual std::vector<double> getState();

	virtual void draw();

	virtual void checkCollision(Rigid&);
	virtual void collide(Rigid&);
	virtual int getWhatAmI();				//// will return defines suppose Sphere==1,Particle==2 etc...


	Matrix normal;
	double A,B,C,D;         //in order to formulate
};

class StringSystem : public Rigid
{
public:
	StringSystem();

	virtual void calculate_forces();
	virtual void calculate_ddt();

	virtual void setState(std::vector<double> &);
	virtual std::vector<double> getState();

	virtual void draw();

	virtual void checkCollision(Rigid&);
	virtual void collide(Rigid&);
	virtual int getWhatAmI();				//// will return defines suppose Sphere==1,Particle==2 etc...


};








#endif