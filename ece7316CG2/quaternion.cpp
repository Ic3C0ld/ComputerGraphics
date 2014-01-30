#include "quaternion.h"

/////	CLASS QUATERNION	//////////

////	Constructors	//////////////


Quaternion::Quaternion()
{
	for (int i = 0; i < 4; i++)
	{
		q[i] = 0;
	}
}
Quaternion::Quaternion(double w, double x, double y, double z)
{
	q[0] = w;
	q[1] = x;
	q[2] = y;
	q[3] = z;
}
Quaternion::Quaternion(double* wxyz)
{
	for (int i = 0; i < 4; i++)
	{
		q[i] = wxyz[i];
	}
}
Quaternion::Quaternion(const Quaternion& other)
{
	for (int i = 0; i < 4; i++)
	{
		q[i] = other.q[i];
	}
}

////////////////////////////////
////////////////////////////////

////	Friend Functions	/////////////
Quaternion operator+(Quaternion& q1, Quaternion& q2)
{
	Quaternion sum;

	for (int i = 0; i < 4; i++)
	{
		sum.q[i] = q1.q[i] + q2.q[i];
	}

	return sum;
}
Quaternion operator-(Quaternion& q1, Quaternion& q2)
{
	Quaternion temp;

	for (int i = 0; i < 4; i++)
	{
		temp.q[i] = q1.q[i] + q2.q[2];
	}

	return temp;
}
Quaternion operator*(double k, Quaternion& q1)
{
	Quaternion temp;

	for (int i = 0; i < 4; i++)
	{
		temp.q[i] = q1.q[i] * k;
	}

	return temp;
}
Quaternion operator*(Quaternion& q1, double k)
{
	Quaternion temp;

	for (int i = 0; i < 4; i++)
	{
		temp.q[i] = q1.q[i] * k;
	}

	return temp;
}
Quaternion operator*(Quaternion& q1, Quaternion& q2)
{
	////q1*q2==(w1*w2-V1*V2 , w1*V2+w2*V1 + V1 x V2 )   //Vi={x,y,z}, * dot, x Cross
	
	double w = q1.q[0] * q2.q[0] - 
			 ( q1.q[1] * q2.q[1] +
			   q1.q[2] * q2.q[2] + 
			   q1.q[3] * q2.q[3]);

	double x =  q1.q[0] * q2.q[1] + 
				q1.q[1] * q2.q[0] +
				q1.q[2] * q2.q[3] -
				q1.q[3] * q2.q[2] ;
	
	double y =  q1.q[0] * q2.q[2] + 
				q1.q[2] * q2.q[0] +
				q1.q[3] * q2.q[1] -
				q1.q[1] * q2.q[3] ;

	double z =  q1.q[0] * q2.q[3] + 
				q1.q[3] * q2.q[0] +
				q1.q[1] * q2.q[2] -
				q1.q[2] * q2.q[1] ;

	return Quaternion(w, x, y, z);
}
Quaternion operator/(Quaternion& q1, double k)
{
	Quaternion temp;

	for (int i = 0; i < 4; i++)
	{
		temp.q[i] = q1.q[i]/k;
	}

	return temp;
}

Quaternion Quaternion::operator=(Quaternion& other)
{
	for (int i = 0; i < 4; i++)
	{
		q[i] = other.q[i];
	}

	return *this;
}
Matrix Quaternion2RotMatrix(Quaternion& q)
{
	//just to make the actions more readable
	double w = q.q[0];
	double x = q.q[1];
	double y = q.q[2];
	double z = q.q[3];
		
	//Simplified version : assumes norm()==1 , i believe;
	double mat[] = { 1-2*(y*y+z*z)	,		2*(x*y-w*z)    ,		2*(x*z+w*y)  ,		0,
					 2 * (x*y + w*z),	  1-2 * (x*x + z*z),		2 * (y*z-w*x),		0,
					 2 * (x*z - w*y),		2*(y*z+w*x) ,			1-2*(x*x+y*y),		0,
					 0				,			0		,				0		 ,		1};
					 
	return Matrix(4,4,mat);
};

Matrix QRotate(double th, double* axis, bool useDegrees)
{
	double theta = th;
	if (useDegrees == true)		theta *= M_PI / 180;

	double norm = sqrtl(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
	axis[0] /= norm;
	axis[1] /= norm;
	axis[2] /= norm;

	double w = cos(theta / 2);
	double x = sin(theta / 2)*axis[0];
	double y = sin(theta / 2)*axis[1];
	double z = sin(theta / 2)*axis[2];

	Quaternion qRot(w, x, y, z);
	//normalize(qRot);

	return Quaternion2RotMatrix(qRot);
}

Matrix QRotate(double th, Matrix& Axis, bool useDegrees)
{
	double theta = th;
	if (useDegrees == true)		theta *= M_PI / 180;
	

	double w = cos(theta / 2);
	double x = sin(theta / 2)*Axis.mat[0] / Axis.norm();
	double y = sin(theta / 2)*Axis.mat[1] / Axis.norm();
	double z = sin(theta / 2)*Axis.mat[2] / Axis.norm();

	Quaternion qRot(w, x, y, z);
	//normalize(qRot);

	return Quaternion2RotMatrix(qRot);
}

Matrix QRotate(double theta, double* axis)
{	
	double norm = sqrtl(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]  );
	axis[0] /= norm;
	axis[1] /= norm;
	axis[2] /= norm;

	double w = cos(theta/2);
	double x = sin(theta /2)*axis[0];
	double y = sin(theta /2)*axis[1];
	double z = sin(theta /2)*axis[2];

	Quaternion qRot(w, x, y, z);
	//normalize(qRot);

	return Quaternion2RotMatrix(qRot);
}
Matrix QRotate(double theta, Matrix& Axis)
{
	double w = cos(theta / 2);
	double x = sin(theta/2)*Axis.mat[0] / Axis.norm();
	double y = sin(theta / 2)*Axis.mat[1] / Axis.norm();
	double z = sin(theta / 2)*Axis.mat[2] / Axis.norm();

	Quaternion qRot(w, x, y, z);
	//normalize(qRot);

	return Quaternion2RotMatrix(qRot);
}


Matrix QRotate(double theta, double* axis, double* point)
{
	double pMinus[] = { -point[0], -point[1], -point[2] };

	return (translate(point)*QRotate(theta, axis)*translate( pMinus));
}


Matrix QRotate(double th, double* axis, double* point, bool useDegrees)
{
	double theta = th;
	if (useDegrees == true)		 theta *= M_PI / 180;


	double pMinus[] = { -point[0], -point[1], -point[2] };

	return (translate(point)*QRotate(theta, axis)*translate(pMinus));

}







void normalize(Quaternion& q)
{
	double norm = q.norm();

	for (int i = 0; i < 4; i++)
	{
		q.q[i] /= norm;
	}

}
///////

Quaternion  Quaternion::conjugate()
{
	// q_={w, -xyz}
	return Quaternion(q[0], -q[1], -q[2], -q[3]);
}
double Quaternion::norm()
{
	return sqrtl((q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]  ));
}
Quaternion  Quaternion::inverse()
{
	return ( conjugate() / norm() );
}


