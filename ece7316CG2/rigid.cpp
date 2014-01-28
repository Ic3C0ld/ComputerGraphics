#include "rigid.h"

/////////	Class RIGID		////////////////////////

Rigid::Rigid()
{
	m_justCollided = false;
	m_color[0] = m_color[1] = m_color[2] = 1;

	m_id = RIGID_BASE;

	p_plane = NULL;
	p_sphere = NULL;
	p_particle = NULL;
	p_spring = NULL;

}


//////////////////////////////// END RIGID ////////


/////////	Class PLANE		////////////////////////
Plane::Plane(double size)
		:Rigid()
{
	m_id = RIGID_PLANE;
	m_color[0] =m_color[1]=m_color[2] =1;
	m_justCollided = false;
	p_plane = this;

	//Create sample 1x1 square at center Normal= {0,1,0} as a vector it has w=0,
	//while the actual points have w=1 in homogeneous representation
	double planeMat[] = {	
							0,	0,	0.5,	 0.5,	-0.5,	-0.5,	0.5,
							1,	0,	  0,	   0,	   0,	   0,	0,
							0,	0,	0.5,	-0.5,	-0.5,	 0.5,	0.5,
							0,	1,	  1,	   1,	   1,	   1,	1	};

	//the point {0,0,0} is so that i can draw with GL_TRIANGLE_FAN

	m_plane = Matrix(4, 7, planeMat);
	m_plane = scale(size, size, size)*m_plane; //// scale(size)=size*Identity4x4

	//hack normalize of normal vector
	m_plane.mat[0 + 1 * m_plane.m_columns] /= size;
}
Plane::Plane(double size, double color[]) 
		:Rigid()
{
	m_id = RIGID_PLANE;
	m_color[0] = color[0];	m_color[1] = color[1];	m_color[2] = color[2];
	m_justCollided = false;
	p_plane = this;

	//Create sample 1x1 square at center Normal= {0,1,0} as a vector it has w=0,
	//while the actual points have w=1 in homogeneous representation
	double planeMat[] = {	
							0,	0,	0.5,	 0.5,	-0.5,	-0.5,	0.5,
							1,	0,	  0,	   0,	   0,	   0,	0,
							0,	0,	0.5,	-0.5,	-0.5,	 0.5,	0.5,
							0,	1,	  1,	   1,	   1,	   1,	1	};

	//the point {0,0,0} is so that i can draw with GL_TRIANGLE_FAN
	
	m_plane = Matrix(4, 7, planeMat);
	m_plane = scale(size, size, size)*m_plane; //// scale(size)=size*Identity4x4

	//hack normalize of normal vector
	m_plane.mat[0 + 1 * m_plane.m_columns] /= size; //{ 0,1,0	,0} *scale(size)  /size


	
}

void Plane::draw()
{
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(m_color[0], m_color[1], m_color[2]);

	glPolygonMode(GL_FRONT, GL_LINE);

	glPushMatrix();
		
		glBegin(GL_TRIANGLE_FAN);
		for (int i = 1; i < m_plane.m_columns; i++)
		{
			glVertex3f(m_plane.mat[i], m_plane.mat[i+m_plane.m_columns], m_plane.mat[i + 2 * m_plane.m_columns]);
		}
		glEnd();

	glPopMatrix();

	glPolygonMode(GL_FRONT, GL_FILL);

	glPopAttrib();
}

void Plane::checkCollision(Rigid*)
{

}
void Plane::calcCollisionResponse(Rigid*)
{
	 
}
void Plane::applyCollisionResponse()
{

}
void Plane::update(double dt)
{

}


//////////////////////////////// END PLANE ////////
