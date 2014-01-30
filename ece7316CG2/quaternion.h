#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>
#include "matrix.h"

////FRIEND FUNCTIONS////////////

 class Quaternion;

 Quaternion operator*(Quaternion& q1, Quaternion& q2);
 Quaternion operator*(Quaternion& q1, double k);
 Quaternion operator*(double k, Quaternion& q1);

 Quaternion operator/(Quaternion& q1, double k);

 Quaternion operator+(Quaternion& q1, Quaternion& q2);
 Quaternion operator-(Quaternion& q1, Quaternion& q2);

////
 Matrix Quaternion2RotMatrix(Quaternion& q);		//R(4,4)
 Matrix Quaternion2MatrixVector(Quaternion& q);	//->matrix v(4,1) 

 Quaternion Matrix2Quaternion(Matrix& v);			//matrix v(4,1) only

  void QRotate(double theta, double Axisx, double Axisy, double Axisz, double& Px, double& Py, double& Pz);
 
  void normalize(Quaternion&);

  Matrix QRotate(double theta, double* axis);
  Matrix QRotate(double theta, Matrix& Axis);

  Matrix QRotate(double theta, double* Axis, double* Pxyz);
  Matrix QRotate(double theta, Matrix& Axis, Matrix& Point); //Axis/Point 4x1, or 3x1 doesnt matter 

   Matrix QRotate(double theta, double* axis, double* point, bool useDegrees);
   Matrix QRotate(double theta, Matrix& axis, Matrix& point, bool useDegrees);
//////////////////////////////////////

class Quaternion
{
public:
	Quaternion();
	Quaternion(double w, double x, double y, double z);
	Quaternion(double* wxyz);
	Quaternion(const Quaternion& );
	// QUATERNION  Q={w, x,y,z}
	double q[4];
	//
	Quaternion conjugate();
	Quaternion inverse();

	double norm();
	//
	Quaternion operator=(Quaternion& other);

	friend Quaternion operator*(Quaternion& q1, Quaternion& q2);
	friend Quaternion operator*(Quaternion& q1, double k);
	friend Quaternion operator*(double k, Quaternion& q1);

	friend Quaternion operator/(Quaternion& q1, double k);

	friend Quaternion operator+(Quaternion& q1, Quaternion& q2);
	friend Quaternion operator-(Quaternion& q1, Quaternion& q2);

	////
	friend Matrix Quaternion2RotMatrix(Quaternion& q);		//R(4,4)
	friend Matrix Quaternion2MatrixVector(Quaternion& q);	//->matrix v(4,1) 

	friend Quaternion Matrix2Quaternion(Matrix& v);			//matrix v(4,1) only

	friend void normalize(Quaternion&);


	friend void QRotate(double theta, double Axisx, double Axisy, double Axisz, double& Px, double& Py, double& Pz);
	
	friend Matrix QRotate(double theta,double* axis);
	friend Matrix QRotate(double theta, Matrix& Axis);

	friend Matrix QRotate(double theta, double* axis, double* point);
	friend Matrix QRotate(double theta, Matrix& axis, Matrix& point);

	friend Matrix QRotate(double theta, double* axis, double* point, bool useDegrees);
	friend Matrix QRotate(double theta, Matrix& axis, Matrix& point, bool useDegrees);

	friend Matrix QRotate(double theta, double* axis, bool useDegrees);
	friend Matrix QRotate(double theta, Matrix& Axis, bool useDegrees);
};





#endif