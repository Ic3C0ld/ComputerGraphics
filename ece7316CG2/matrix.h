#ifndef MATRIX_H
#define MATRIX_H

//////INCLUDES////////////////////////
#include <cmath>
#include <vector>
#include "glut.h"
////
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616




class Matrix;

//////Matrix Method declarations/////


Matrix operator+(Matrix&, Matrix&);
Matrix operator+(Matrix&, double);
Matrix operator+(double, Matrix&);

Matrix operator-(Matrix&, Matrix&);
Matrix operator-(Matrix&, double);
Matrix operator-(double, Matrix&);

Matrix operator*(Matrix&, Matrix&);
Matrix operator*(Matrix&, double);
Matrix operator*(double, Matrix&);

Matrix operator/(Matrix&, double);
Matrix cross(Matrix&, Matrix&);

void transpose(Matrix& toBeTransposed);
void invert(Matrix& toBeInverted);
void normalize(Matrix& toBeNormalized);
double determinant(Matrix&);
/////OTHER///////
//Matrix vector3DtoMatrix(Vector3D&);
//Vector3D MatrixtoVector3D(Matrix&);
//common Matrices ,each return 
Matrix Identity();
Matrix Zeros();
Matrix Ones();
//
//Rotation // Translation

Matrix translate(double dx, double dy, double dz);
Matrix translate(double* dXdYdZ);
Matrix rot(double theta, double axis[], double point[]);
Matrix rotT(double w, double axis[], double point[]);
Matrix scale(double x, double y, double z);

//
Matrix rotX(double theta, bool useDegrees = false);
Matrix rotY(double theta, bool useDegrees = false);
Matrix rotZ(double theta, bool useDegrees = false);



//////DECLARATIONS////////////////
class Matrix
{
public:
	Matrix();															//default constructor of 1x1 Matrix equal to 0;
	Matrix(unsigned rows, unsigned columns, double* MatTable);	//
	Matrix(unsigned Rows, unsigned Columns);						//eg.  as above where dimensions=Sizes.size()
	Matrix(const Matrix&);												//

	//~Matrix();


////ALL INFORMATION IS HERE////////
//// Actual Matrix ////// N=rows M= colums
	// eg. 4x4 Matrix //
	// [ 0  1  2  3 ]   
	// [ 4  5  6  7 ]	
	// [ 8  9  10 11]
	// [ 12 13 14 15]
	std::vector<double> mat; //row major
	unsigned m_rows, m_columns;	//this sets the N,M
	unsigned m_size;			//N*M

	//
//	std::vector<double*> allocatedArrays;//from   double* getMat();
/////////////////////////////////////


////METHODS///////
	Matrix& operator=(const Matrix& );
	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const Matrix&);

	double norm();//has meaning if Matrix Nx1 or 1xM...
	
	Matrix inverse();
	Matrix transpose();
	double det();

	void DeleteRow(unsigned);
	void DeleteColumn(unsigned);
	Matrix getRow(unsigned whichRow);
	Matrix getColumn(unsigned whichColumn);
	void addRow(Matrix& rowVector,unsigned where2put);
	void addColumn(Matrix&, unsigned where2put);
	double* getMat();

	bool isVector();




	friend Matrix operator+(Matrix&, Matrix&);
	friend Matrix operator-(Matrix&, Matrix&);
	
	friend Matrix operator*(Matrix&, Matrix&);
	friend Matrix operator*(Matrix&, double);
	friend Matrix operator*(double, Matrix&);

	friend Matrix operator/(Matrix&, double);
	friend Matrix cross(Matrix&, Matrix&);



	friend void normalize(Matrix&);//has meaning if Matrix Nx1 or 1xM...
	friend void transpose(Matrix&); //implementation for 2-dimensions only ,
	friend void invert(Matrix&);	  //so is this implementation
	friend double determinant(Matrix&);

	
	//common Matrices ,each  
	friend void setIdentity(Matrix&); //  2-dimensions only 
	friend void setZeros(Matrix&);
	friend void setOnes(Matrix&);
	//
	
	//
	//// Translation / Rotation etc //// 4x4 Matrices  ///////////////

	friend Matrix translate(double dx, double dy, double dz);
	friend Matrix translate(double*);
	friend Matrix rot(double theta, double axis[], double point[]);
	friend Matrix rotTdegrees(double w, double axis[], double point[]);
	friend Matrix scale(double x, double y, double z);
	//////OTHER
	//friend Matrix vector3DtoMatrix(Vector3D&);
	//friend Vector3D MatrixtoVector3D(Matrix&);
	//
private:
	
};






#endif