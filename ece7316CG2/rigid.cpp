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

	glDisable(GL_CULL_FACE);
	glLineWidth(5);
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
	glLineWidth(1);
	glEnable(GL_CULL_FACE);

	glPopAttrib();
}
void Plane::checkCollision(Rigid* obj)
{
	switch (obj->m_id)
	{
	case RIGID_SPHERE:
	{
		Matrix vector = m_plane.getColumn(1);
		vector = vector - obj->p_sphere->X_t;
		double dotProduct = (vector.transpose()*m_plane.getColumn(0)).mat[0];

		if (dotProduct < obj->p_sphere->m_radius)
		{
			calcCollision(this,obj->p_sphere);
		}
		break;

		
	}
	default:
		break;
	}

}
void Plane::applyCollisionResponse()
{

}
void Plane::update(double dt)
{

}


//////////////////////////////// END PLANE ////////

/////////	Class SPHERE	////////////////////////

Sphere::Sphere(double mass, double radius, double Pxyz[], double Vxyz[],double color3[])
		:Rigid()
{
	p_sphere = this;
	m_id = RIGID_SPHERE;

	m_radius = radius;
	m_mass = mass;
	
	m_color[0] = color3[0];	m_color[1] = color3[1];	m_color[2] = color3[2];

	double xMat[] = { Pxyz[0], Pxyz[1], Pxyz[2] ,1};
	double vMat[] = { Vxyz[0], Vxyz[1], Vxyz[2], 0};

	X_t = Matrix(4, 1, xMat);
	V_t = Matrix(4, 1, vMat);

	applyV_t = Matrix(4, 1);//zero-ed
}

void Sphere::draw()
{
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(m_color[0], m_color[1], m_color[2]);

	glPushMatrix();
		glTranslatef(X_t.mat[0], X_t.mat[1], X_t.mat[2]);
		glutSolidSphere(m_radius, 5+m_radius, 5+m_radius);
	glPopMatrix();

	glPopAttrib();

}
void Sphere::checkCollision(Rigid* obj)
{
	switch (obj->m_id)
	{
		case RIGID_PLANE:
		{
			Matrix vector = obj->p_plane->m_plane.getColumn(1);
			vector = vector - X_t;
			double dotProduct = (vector.transpose()*obj->p_plane->m_plane.getColumn(0)).mat[0];

			if (dotProduct < m_radius)
			{
				calcCollision(obj->p_plane, this);
			}
			break;
		}
		case RIGID_SPHERE:
		{
			Matrix r12 = X_t - obj->p_sphere->X_t;

			if (r12.norm() < (m_radius + obj->p_sphere->m_radius))
			{
				calcCollision(this, obj->p_sphere);
			}

			break;
		}
		default:
			break;
	}
}
void Sphere::applyCollisionResponse()
{
	V_t = V_t + applyV_t;

	setZeros(applyV_t);

}
void Sphere::update(double dt)
{
	X_t = X_t+ V_t*dt;
}

//////////////////////////////// END SPHERE ////////



/////////	Class PARTICLE	////////////////////////
Particle::Particle(int PartCount,double radius,double totalMass,double cmPosition[],double cmVel[],double W[])
:Rigid()
{
	m_id = RIGID_PARTICLE;
	m_mass = totalMass;
	m_radius = radius; //all spheres will have same size and mass

	x_t =Matrix(4, 1, cmPosition);
	v_t=Matrix(4, 1, cmVel);
	
	w_t=Matrix(4,1,W);
	
	rotation = Quaternion(1, 0, 0, 0);
	R_t = Quaternion2RotMatrix(rotation);


	if ((PartCount % 2)==0)
	{
		double point0[] = { -radius, 0, 0, 0 };
		double point1[] = { +radius, 0, 0, 0 };
		
		points = Matrix(4, 1, point0);
		points.addColumn(Matrix(4, 1, point1), 1);

		for (int i = 2; i < PartCount; i++)
		{
			point0[0] = pow(-1,i)*(i * 2 * radius -radius);
			points.addColumn(Matrix(4, 1, point0), 0);
			
		}

	}
	else
	{
		double point0[] = { 0, 0, 0, 0 };
		points = Matrix(4, 1, point0);

		for (int i = 1; i < 0.5*PartCount; i++)
		{
			point0[0] = -i * 2 * radius;
			points.addColumn(Matrix(4, 1, point0),0);
			point0[0] = +i * 2 * radius;
			points.addColumn(Matrix(4, 1, point0), 0);
		}

	}

	I = I_1 =  Matrix(4, 4);
	setIdentity(I);
	setIdentity(I_1);

	currentPoints = points;

}

void Particle::draw()
{
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(m_color[0], m_color[1], m_color[2]);


	double x,y,z;

		for (int i = 0; i < currentPoints.m_columns; i++)
		{
			glPushMatrix();
			x = currentPoints.mat[i + currentPoints.m_columns * 0];
			y = currentPoints.mat[i + currentPoints.m_columns * 1];
			z = currentPoints.mat[i + currentPoints.m_columns * 2];

			glTranslatef(x, y, z);
			glutSolidSphere(m_radius, 10, 10);
				
			glPopMatrix();

		}

	
	glPopAttrib();

}

void Particle::update(double dt)
{

}
void Particle::checkCollision(Rigid* obj)
{

}
void Particle::applyCollisionResponse()
{

}


//////////////////////////////// END PARTICLE ////////



///////////////////// GLOBAL ///////////////////////////////////////////

void calcCollision(Plane* plane, Sphere* sphere)
{
	
	double nSpeed = (sphere->V_t.transpose()*plane->m_plane.getColumn(0)).mat[0];

	//sphere->applyV_t = sphere->applyV_t - 2 * nSpeed* plane->m_plane.getColumn(0);


	Matrix planeNormal = plane->m_plane.getColumn(0);

	sphere->V_t = sphere->V_t -2 * nSpeed* plane->m_plane.getColumn(0);
}
void calcCollision(Sphere* s1, Sphere* s2)
{
	Matrix r12 = s2->X_t - s1->X_t;
	double m1 = s1->m_mass;
	double m2 = s2->m_mass;

	double error = s1->m_radius + s2->m_radius - r12.norm();

	normalize(r12);

	
	s1->X_t = s1->X_t -(r12*error*(s2->m_radius) / (s1->m_radius + s2->m_radius));
	s2->X_t = s2->X_t +(r12*error*(s1->m_radius) / (s1->m_radius + s2->m_radius));



	double v1prev = (s1->V_t.transpose()*r12).mat[0];
	double v2prev = (s2->V_t.transpose()*r12).mat[0];

	
	double v1after = (v1prev*(m1 - m2) + 2 * m2*v2prev) / (m1 + m2);
	double v2after = (v2prev*(m2 - m1) + 2 * m1*v1prev) / (m1 + m2);


	s1->V_t = s1->V_t + (v1after - v1prev)*r12;
	s2->V_t = s2->V_t + (v2after - v2prev)*r12;


}