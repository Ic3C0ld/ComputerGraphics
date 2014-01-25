#include "matrix.h"

////	MATRIX Constructors	 ////////////////////////// 

Matrix::Matrix()
{
	m_size = 1;
	m_rows = 1; m_columns = 1; // x2 so that N*M= 1*1 == m_size

	mat.push_back((double)0);

}											
Matrix::Matrix(unsigned rows, unsigned columns, double* MatTable)
{
	//ensure they are not zero size
	m_rows = (rows > 1) ? rows : 1;
	m_columns = (columns>1) ? columns : 1;
	m_size = m_rows* m_columns;

	for (int i = 0; i < m_size; i++)
	{
		mat.push_back( MatTable[i] );
	}

}	

Matrix::Matrix(unsigned Rows,unsigned Columns)	
{
	//ensure rows/columns >=1;
	m_rows = Rows > 1 ? Rows : 1;
	m_columns = Columns > 1 ? Columns : 1;
	m_size = m_rows*m_columns;

	//Initialize mat to 0s
	mat.assign(m_size, (double)0);

}			
Matrix::Matrix(const Matrix& other)
{
	m_size = other.m_size;
	m_rows = other.m_rows;
	m_columns = other.m_columns;

	for (int i = 0; i < m_size; i++)
	{
		mat.push_back(other.mat[i]);
	}

}										

//Matrix::~Matrix()
//{
//	double* pointer;
//	for (int i = 0; i < allocatedArrays.size(); i++)
//	{
//		pointer = allocatedArrays[i];
//		delete[] pointer;
//	}
//}
/////////////////  MATRIX  ////////////////////

Matrix& Matrix::operator=(const Matrix& other)
{
	m_size = other.m_size;
	m_rows = other.m_rows;
	m_columns = other.m_columns;

	//SELF ASSIGNMENT PROTECTION :) if A=A ,i cannot delete before i copy
	std::vector<double>tempAssign;

	for (int i = 0; i < m_size; i++)
	{
		tempAssign.push_back(other.mat[i]);
	}
	
	mat.clear();

	for (int i = 0; i < m_size; i++)
	{
		mat.push_back(tempAssign[i]);
	}

	return *this;
}

///////////////// MATRIX's FRIENDS :) ////////////////
Matrix operator+(Matrix& A, Matrix& B)
{
	Matrix temp(A.m_rows,A.m_columns);

	if ( (A.m_rows == B.m_rows) &&  (A.m_columns == B.m_columns) )
	{
		for (int i = 0; i < A.m_size; i++)
		{
			temp.mat[i] = A.mat[i] + B.mat[i];
		}
	}
	
	return temp;

}

Matrix operator-(Matrix& A, Matrix& B)
{
	Matrix temp(A.m_rows, A.m_columns);

	if ((A.m_rows == B.m_rows) && (A.m_columns == B.m_columns))
	{
		for (int i = 0; i < A.m_size; i++)
		{
			temp.mat[i] = A.mat[i] - B.mat[i];
		}
	}

	return temp;
}

Matrix operator*(Matrix& A, Matrix& B)
{
	//create temp Matrix of size A.rows,B.columns
	Matrix temp( A.m_rows,B.m_columns );

	if (A.m_columns == B.m_rows)
	{
		for (int i = 0; i < A.m_rows; i++) //trasverse through A.rows 
		{
			for (int j = 0; j < B.m_columns; j++) //trasverse through columns ,could be j<B.colCount in general
			{
				//now we can say that temp.mat[i][j]=temp.mat[B.m_columns*i+j]
				double tempSum = 0;//

				for (int k = 0; k < A.m_columns; k++) //trasverse through cols of A and rows of B ,could be k<A.colCount same as k<B.rowCount in general
				{
					tempSum += A.mat[A.m_columns * i + k] * B.mat[j + B.m_columns * k];
				}

				temp.mat[temp.m_columns* i + j] = tempSum;
			}

		}
	}

	return temp;

}
Matrix operator*(Matrix& A, double k)
{
	Matrix temp( A.m_rows, A.m_columns );

	for (int i = 0; i < A.m_size; i++)
	{
		temp.mat[i] = k*A.mat[i];
	}

	return temp;
}
Matrix operator*(double k, Matrix& A)
{
	Matrix temp( A.m_rows, A.m_columns );

	for (int i = 0; i < A.m_size; i++)
	{
		temp.mat[i] = k*A.mat[i];
	}

	return temp;
}

Matrix operator/(Matrix& A, double k)
{
	Matrix temp( A.m_rows, A.m_columns );

	for (int i = 0; i < A.m_size; i++)
	{
		temp.mat[i] = A.mat[i]/k;
	}

	return temp;
}

Matrix cross(Matrix& v1, Matrix& v2)
{
	double x1 = v1.mat[0]; double y1 = v1.mat[1]; double z1 = v1.mat[2];
	double x2 = v2.mat[0]; double y2 = v2.mat[1]; double z2 = v2.mat[2];


	double mat[] = { y1*z2-z1*y2,	
					 z1*x2-x1*z2,
					 x1*y2-y1*x2 };

	return Matrix(3, 1, mat);
}

void transpose(Matrix& A )
{
	// Á [Col0 ..															[Col0  Col1  Col2]
	//   [Col0 ..		-> transpose -> [Col0 Col1 Col2]  -> append ->      [ ..	..	  .. ]
	//	 [Col0 ..															[ ..    ..    .. ]
	

	// SLOW IMPLEMENTATION
	//Matrix temp=A.getColumn(0);
	//temp.m_columns = temp.m_rows;//
	//temp.m_rows = 1;			// Hack to traspose first column into a row

	//for (int i = 1; i < A.m_columns; i++)
	//{
	//	//Yet again need to transpose the vector
	//	Matrix tempVector = A.getColumn(i);
	//	tempVector.m_columns = tempVector.m_rows;
	//	tempVector.m_rows = 1;

	//	temp.addRow(tempVector,i);
	//}
	//A = temp;



	Matrix T(A.m_columns,A.m_rows);

	for (int i = 0; i < T.m_rows; i++)
	{
		for (int j = 0; j < T.m_columns; j++)
		{
			T.mat[i*T.m_columns + j] = A.mat[j*T.m_rows + i];
		}
	}

	A = T;

}
void invert(Matrix& A)
{
	if (determinant(A) == 0 || (A.m_rows!=A.m_columns)) //Actually need to check too if max(Coffactors[])* 1/det ->OVERFLOW
	{
		
	}
	else
	{
		Matrix Cofactors(A.m_rows, A.m_columns);

		//Compute the Coffactor elementes
		for (int i = 0; i < Cofactors.m_columns; i++)
		{
			for (int j = 0; j < Cofactors.m_rows; j++)
			{
				Matrix temp = A;
				double coefficient = pow(-1.0, i + j) ;

				temp.DeleteColumn(i);
				temp.DeleteRow(j);

				Cofactors.mat[j*A.m_columns + i] = coefficient*determinant(temp);

			}
		}

		Matrix Inverse = Cofactors.traspose() / determinant(A);

		A = Inverse;
	}
	

}
double determinant(Matrix& A)
{
	if ((A.m_rows != A.m_columns) || A.m_rows<1)//no need to check columns too for <1
	{
		return 0;
	}
	else if (A.m_rows == 1){ return A.mat[0]; }
	//else if (A.m_rows == 2)
	//{
	//	//		 A
	//	//		[ 0  1 ]
	//	//		[ 2  3 ]
	//	//
	//	return (A.mat[0] * A.mat[3] - A.mat[1] * A.mat[2]);
	//}
	else
	{
		double det=0;

		for (int i = 0; i < A.m_columns; i++)
		{
			Matrix temp = A;
			double coefficient = pow(-1.0, i)*A.mat[i];

			temp.DeleteColumn(i);
			temp.DeleteRow(0);

			det += coefficient*determinant(temp);//recursion here :D
		}

		return det;
	}

}

void setIdentity(Matrix& A)
{
	unsigned least = A.m_rows < A.m_columns ? A.m_rows : A.m_columns;

	for (int i = 0; i < least; i++)
	{
		A.mat[i] = 1;
	}


}
void setZeros(Matrix& A)
{
	for (int i = 0; i < A.m_size; i++)
	{
		A.mat[i] = 1;
	}
}
void setOnes(Matrix& A)
{
	for (int i = 0; i < A.m_size; i++)
	{
		A.mat[i] = 1;
	}

}

void Matrix::DeleteRow(unsigned row)
{
	unsigned start = row*m_columns;//size of dimension 0= rows, how many cells it holds
	unsigned finish = start + m_columns;
	
	mat.erase(mat.begin()+start, mat.begin()+finish);

	//update Matrix information
	m_size = mat.size();
	m_rows -= 1;
}
void Matrix::DeleteColumn(unsigned col)
{
	for (int i = m_rows-1; i >=0; i--)//inverse deletion of cells because eg. 1st. delete(9th)= X(9), 2nd. delete(9th)= X(10) etc...
	{
		mat.erase(mat.begin() + (col + i*m_columns) );
	}
	//update information
	m_columns -= 1;
	m_size = mat.size();

}
void Matrix::addRow(Matrix& rowVector, unsigned where2put)
{
	if ( (rowVector.isVector() == true) && (rowVector.m_columns == m_columns))
	{
		for (int i = 0; i < m_columns; i++)
		{
			mat.insert(mat.begin() + where2put*m_columns+i, rowVector.mat[i]);
		}
		m_rows += 1;
		m_size = mat.size();
	}
}
void Matrix::addColumn(Matrix& columnVector, unsigned where2put)
{
	if (columnVector.isVector() == true && columnVector.m_rows == m_rows)
	{
		for (int i = 0; i < m_rows; i++)
		{
			mat.insert(mat.begin() + i*(m_columns+1) + where2put, columnVector.mat[i]);
		}

		m_columns += 1;
		m_size = mat.size();
	}

}
Matrix Matrix::getRow(unsigned whichRow)
{
	Matrix temp(1, m_columns);
	
	for (int i = 0; i < m_columns; i++)
	{
		temp.mat[i] = mat[whichRow*m_columns+i];
	}
	
	return temp;
}
Matrix Matrix::getColumn(unsigned whichColumn)
{
	Matrix temp(m_rows, 1);

	for (int i = 0; i < m_rows; i++)
	{
		temp.mat[i] = mat[i*m_columns + whichColumn];
	}

	return temp;
}
bool Matrix::isVector()
{
	if ( (m_rows == 1) || (m_columns == 1))
	{
		return true;
	}
	else
	{
		return false;
	}
}
Matrix Matrix::traspose()
{
	/*Matrix transposed = (*this);
	transpose(transposed);
	return transposed;*/

	if (m_rows == 1 || m_columns == 1)
	{
		Matrix T(*this);
		T.m_rows = m_columns;
		T.m_columns = m_rows;
		return T;	
	}

	Matrix T(m_columns, m_rows);

	for (int i = 0; i < T.m_rows; i++)
	{
		for (int j = 0; j < T.m_columns; j++)
		{
			T.mat[i*T.m_columns + j] = mat[j*T.m_rows + i];
		}
	}

	return T;
}
Matrix Matrix::inverse()
{
	if (determinant(*this) == 0 || (m_rows != m_columns)) //Actually need to check too if max(Coffactors[])* 1/det ->OVERFLOW
	{
		return *this;//do nothing
	}
	else
	{
		Matrix Cofactors(m_rows, m_columns);

		//Compute the Coffactor elementes
		for (int i = 0; i < Cofactors.m_columns; i++)
		{
			for (int j = 0; j < Cofactors.m_rows; j++)
			{
				Matrix temp = *this;
				double coefficient = pow(-1.0, i + j)  /* *A.mat[j*A.m_columns + i]*/;

				temp.DeleteColumn(i);
				temp.DeleteRow(j);

				Cofactors.mat[j*m_columns + i] = coefficient*determinant(temp);

			}
		}

		Matrix Inverse = Cofactors.traspose() / determinant(*this);

		return Inverse;
	}



}
double Matrix::det()
{
	return determinant(*this);
}

double* Matrix::getMat()
{
	////??? WHEN DO I NEED TO FREE THIS///////
	double* table;
	table = new double[m_size];

	for (int i = 0; i < m_size; i++)
	{
		table[i] = mat[i];
	}

	//allocatedArrays.push_back(table); //After testing allocation i dont seem to have any problem ,so...

	return table;
}

//// SPECIAL Vector Methods///////////////////////

double Matrix::norm()
{
	double norm = 0;
	for (int i = 0; i < m_size; i++)
	{
		norm += mat[i]*mat[i];
	}

	return sqrtl(norm);
}
void normalize(Matrix& toBeNormalized)
{
	double dividor = toBeNormalized.norm();
	for (int i = 0; i < toBeNormalized.m_size; i++)
	{
		toBeNormalized.mat[i] /= dividor;
	}
}


//// Translation / Rotation etc /////////////////


Matrix translate(double dx, double dy, double dz)
{
	double mat[] = { 1, 0, 0, dx,
					 0, 1, 0, dy,
					 0, 0, 1, dz,
					 0, 0, 0, 1};

	Matrix T(4, 4, mat);
	
	return T;
}

Matrix translate(double* dXdYdZ)
{
	double mat[]= { 1, 0, 0, dXdYdZ[0],
					0, 1, 0, dXdYdZ[1],
					0, 0, 1, dXdYdZ[2],
					0, 0, 0,	 1		};

	Matrix T(4, 4, mat);

	return T;
}


Matrix rotX(double w, bool useDegrees)
{
	if (useDegrees==true)
	{
		w = w*(M_PI) / 180.0;
	}
	
	double mat[] = {1,		0,		0,			 0,
					0,  cos(w),  -sin(w),		 0,
					0,  sin(w),   cos(w),		 0,
					0,		0,		0,			 1};

	Matrix Rx(4, 4, mat);

	return Rx;

}

Matrix rotY(double w, bool useDegrees)
{
	if (useDegrees == true)
	{
		w = w*(M_PI) / 180.0;
	}

	double mat[] = { cos(w),	 0,		sin(w),		0,
					 0,			 1,		0,			0,
					 -sin(w),	 0,		cos(w),		0,
					 0,			 0,		 0,			1 };

	Matrix Rx(4, 4, mat);

	return Rx;

}

Matrix rotZ(double w, bool useDegrees)
{
	if (useDegrees == true)
	{
		w = w*(M_PI) / 180.0;
	}

	double mat[] = { cos(w),  -sin(w),	0,		0,
					 sin(w),   cos(w),  0,		0,
					 0,			0,		1,		0,
					 0,			0,		0,		1};
					 
	Matrix Rx(4, 4, mat);

	return Rx;

}

Matrix rotTdegrees(double w, double axis[], double point[])//rotation Transposed and argument w in degrees 
{
	////
	double matrix[16];

	glPushMatrix();
	glLoadIdentity();

	glTranslated(point[0], point[1], point[2]);
	glRotated(w, axis[0], axis[1], axis[2]);
	glTranslated(-point[0], -point[1], -point[2]);

	glGetDoublev(GL_MODELVIEW_MATRIX, matrix);////this might create issues on other machines but here GLdouble is  " typedef double GLDouble " 

	glPopMatrix();

	Matrix C(4, 4, matrix);

	return C;

}


Matrix rot(double w, double axis[], double point[])
{
	////
	double theta = w * 180 / M_PI; //needed as glRotate* takes degrees
	double matrix[16];

	glPushMatrix();
	glLoadIdentity();

	glTranslated(point[0], point[1], point[2]);
	glRotated(theta, axis[0], axis[1], axis[2]);
	glTranslated(-point[0], -point[1], -point[2]);

	glGetDoublev(GL_MODELVIEW_MATRIX, matrix);

	glPopMatrix();

	Matrix C(4, 4, matrix);

	return C.traspose();

	
}


//// OTHER //////////////
//Matrix vector3DtoMatrix(Vector3D& v)
//{
//	double tempV[] = { v.x, v.y, v.z ,0};
//	Matrix A(4, 1, tempV);
//
//	return A;
//}
//Vector3D MatrixtoVector3D(Matrix& A)
//{
//	Vector3D v(A.mat[0], A.mat[1], A.mat[2]);
//	return v;
//}