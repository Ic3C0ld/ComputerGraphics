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
	glColor4f(m_color[0], m_color[1], m_color[2],0.2);

	glDisable(GL_CULL_FACE);
	glLineWidth(4);
	glPolygonMode(GL_FRONT, GL_LINE);

	glPushMatrix();
		
		glBegin(GL_TRIANGLE_FAN);
		for (int i = 1; i < m_plane.m_columns; i++)
		{
			glVertex3f(m_plane.mat[i], m_plane.mat[i+m_plane.m_columns], m_plane.mat[i + 2 * m_plane.m_columns]);
		}
		glEnd();

		Matrix p1 = m_plane.getColumn(1);
		Matrix p2 =p1 + 10*m_plane.getColumn(0);//show the normals direction //just to be sure 

		glColor3f(1, 1, 1);
		glLineWidth(2);

		glBegin(GL_LINES);
		glVertex3f(p1.mat[0], p1.mat[1], p1.mat[2]);
		glVertex3f(p2.mat[0], p2.mat[1], p2.mat[2]);
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
		Sphere *s = static_cast<Sphere*>(obj);

		Matrix n = -1*m_plane.getColumn(0);//normally the planes normal points out of the box.So now it points towards the sphere

		double vdotProduct = (s->V_t.transpose()*n).mat[0];

		if (vdotProduct < 0)//if the sphere heads towards the plane
		{
			Matrix vector = m_plane.getColumn(1);
			vector = vector - s->X_t;

			double dotProduct = (vector.transpose()*m_plane.getColumn(0)).mat[0];

			if (dotProduct <s->m_radius)//checking for actual overlaps
			{
				calcCollision(this, s);
			}
		}
		break;
	}
	case RIGID_PARTICLE:
	{
		//Last try
		Particle *p = static_cast<Particle*>(obj);
		Matrix n = -1 * m_plane.getColumn(0); //so that it points inwards the box











		//// TRY ONE
		//Particle *p = static_cast<Particle*>(obj);

		//Matrix n = -1 * m_plane.getColumn(0);					//get the plane normal facing inwards the box
		//Matrix p1 = m_plane.getColumn(1);						//get one plane points
		//Matrix r12 = ((p->x_t - p1).transpose()*n).mat[0] * n; //r12 = ((x_t-p1)*n)*n == projection on the n vector

		//int N = p->currentPoints.m_columns;



		//std::vector<Matrix> J, RP;

		//for (int i = 0; i < N; i++)
		//{
		//	Matrix testP = p->currentPoints.getColumn(i); //testing all spheres of the particle

		//	//check for collision
		//	Matrix v = testP - p1;
		//	double dotError = (v.transpose()*n).mat[0];

		//	

		//	if (dotError < p->m_radius) //actual intersection
		//	{
		//		Matrix cP = testP - n*dotError;
		//		Matrix rP = testP - p->x_t;

		//		Matrix vRot = cross(p->w_t, rP);
		//		Matrix vRel = -1 * (p->v_t + vRot);

		//		double vrel = (vRel.transpose()*n).mat[0];

		//		if (vrel>0)// the point is going towards the plane,
		//		{

		//			double num = -2 * vrel;
		//			double den1 = 1.0 / p->m_mass;
		//			double den2 = (n.transpose()*p->cI_1*(cross(cross(rP, n), rP))).mat[0];


		//			Matrix t = cross(cross(n, vRel), n);

		//			normalize(t);
		//			Matrix nt = n + 0.0*t;
		//			//normalize(nt);
		//			Matrix  Jtemp = nt* num / (den1 + den2);
		//			J.push_back(Jtemp);
		//			RP.push_back(rP);
 
		//		}
		//	}


		//}

		//for (int i = 0; i < J.size(); i++)
		//{
		//p->v_t = p->v_t - J[i] / p->m_mass/J.size();
		//p->w_t = p->w_t - p->cI_1*cross(RP[i], J[i] / J.size());
		//}

		////TRY TWO
		//Particle *p = static_cast<Particle*>(obj);

		//Matrix n = -1 * m_plane.getColumn(0);					//get the plane normal facing inwards the box
		//Matrix p1 = m_plane.getColumn(1);						//get one plane points
		//Matrix r12 = ((p->x_t - p1).transpose()*n).mat[0] * n; //r12 = ((x_t-p1)*n)*n == projection on the n vector

		//int N = p->currentPoints.m_columns;



		//Matrix rPavg(4,1);
		//for (int i = 0; i < N; i++)
		//{
		//	Matrix testP = p->currentPoints.getColumn(i); //testing all spheres of the particle

		//	//check for collision
		//	Matrix v = testP - p1;
		//	double dotError = (v.transpose()*n).mat[0];

		//	if (dotError < p->m_radius) //actual intersection
		//	{
		//		Matrix cP = testP - n*dotError;
		//		Matrix rP = testP - p->x_t;

		//		Matrix vRot = cross(p->w_t, rP);
		//		Matrix vRel = -1 * (p->v_t + vRot);

		//		double vrel = (vRel.transpose()*n).mat[0];

		//		if (vrel<0)// the point is going towards the plane,
		//		{
		//			rPavg = rPavg + rP;
		//								
		//		}
		//	}
		//}

		//Matrix vRot = cross(p->w_t, rPavg);
		//Matrix vRel = -1 * (p->v_t + vRot);

		//double num = -2 * ((vRel.transpose()*n).mat[0]);
		//double den1 = 1.0 / obj->p_particle->m_mass;
		//double den2 = (n.transpose()*obj->p_particle->cI_1*cross(cross(rPavg, n), rPavg)).mat[0];
		//
		//Matrix t = cross(cross( n,vRel), n);
		//normalize(t);
		//Matrix nt = n + 0.0*t;
		////normalize(nt);
		//Matrix  J = nt* num / (den1 + den2);
		//
		//p->v_t = p->v_t - J/p->m_mass ;
		//p->w_t = p->w_t + p->cI_1*cross(rPavg, -1*J);
		//

//		//TRY 0
//		//since all spheres in particle have the same radius		    |<--->|
//		//the bounding sphere of it all is maximum  pCount*m_radius    ооооо
//		Matrix ppoint = m_plane.getColumn(1);
//		Matrix vector = ppoint - obj->p_particle->x_t;
//		Matrix normal = m_plane.getColumn(0);
//		
//		double dotProduct = (vector.transpose()*normal).mat[0];
//		double pCount = obj->p_particle->currentPoints.m_columns;
//
//		
//
//		/*if (dotProduct < obj->p_particle->m_radius*pCount)
//		{*/
//			bool collided = false;
//
//			//determine which of the spheres is the one potentially colliding
//			for (int i = 0; (i < pCount); i++)
//			{
//				vector = ppoint - obj->p_particle->currentPoints.getColumn(i);
//				dotProduct = (vector.transpose()*normal).mat[0];
//
//				Matrix n = -1 * normal;
//				if (dotProduct < obj->p_particle->m_radius)
//				{
//
//					Matrix cP = obj->p_particle->currentPoints.getColumn(i) + normal*obj->p_particle->m_radius;
//
//					Matrix rP = cP - obj->p_particle->x_t;
//
//					Matrix vRot = cross(obj->p_particle->w_t, rP);
//					Matrix vRel = -1*(obj->p_particle->v_t + vRot); //V = Vlinear + W x  R;
//
//					double dot = (vRel.transpose()*n).mat[0];
//
//					if (dot > 0)
//					{
//						double num = -2 * ((vRel.transpose()*n).mat[0]);
//						double den1 = 1.0 / obj->p_particle->m_mass;
//						double den2 = (n.transpose()*obj->p_particle->cI_1*cross(cross(rP, n), rP)).mat[0];
//
//						Matrix t = cross(cross( n,vRel), n);
//						normalize(t);
//						Matrix nt = n + 0.3*t;
//						//normalize(nt);
//						Matrix  J = nt* num / (den1 + den2);
//
//						obj->p_particle->v_t = obj->p_particle->v_t - J / obj->p_particle->m_mass; /*- J / obj->p_particle->m_mass;*/;
//
//						Matrix dw = obj->p_particle->cI_1*cross(rP, -1*J);
//						obj->p_particle->w_t = obj->p_particle->w_t + dw;
//
//
//						/*Matrix T = n*dotProduct;
//						obj->p_particle->currentPoints = translate(T.mat[0], T.mat[1], T.mat[2])*obj->p_particle->currentPoints;
//*/
//
//						i=pCount;
//					}
//
//
//					
//					
//					////////////////////////////////////////////////////////////////
//					//Matrix cP = obj->p_particle->x_t + dotProduct*normal;
//					////TODO compute necessary variables for collision
//					//Matrix rP = cP - obj->p_particle->x_t;
//					//Matrix vP = obj->p_particle->v_t + cross(obj->p_particle->w_t, rP); //V = Vlinear + W x  R;
//
//					//double num = -2*((vP*n).mat[0]);
//					//double den1 = 1.0 / obj->p_particle->m_mass;
//					//double den2 = (n.transpose()*obj->p_particle->cI_1*cross(cross(rP, n), rP)).mat[0];
//
//					//Matrix  J = num / (den1 + den2) *n;
//
//					//obj->p_particle->v_t = obj->p_particle->v_t - J / obj->p_particle->m_mass;
//					//obj->p_particle->w_t = obj->p_particle->w_t - obj->p_particle->cI_1*cross(cP, J);
//
//					 ////do not allow more collisions on this body //might fix it later somehow
//				}
//				
//			
//			}

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
double Plane::getKinetik()
{
	return 0; //return zero kinetic
}
Matrix Plane::getMomentum()
{
	return Matrix(4, 1); //return 0 vector
}
//////////////////////////////// END PLANE ////////

/////////	Class SPHERE	////////////////////////

Sphere::Sphere(double radius, double mass, double Pxyz[], double Vxyz[],double color3[])
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
		glutSolidSphere(m_radius, 15+m_radius, 15+m_radius);
	glPopMatrix();

	glPopAttrib();

}
void Sphere::checkCollision(Rigid* obj)
{
	switch (obj->m_id)
	{
		case RIGID_PLANE:
		{
			Plane *p = static_cast<Plane*>(obj);
				
			Matrix n = -1 * p->m_plane.getColumn(0);//normally the planes normal points out of the box.So now it points towards the sphere

			double vdotProduct = (V_t.transpose()*n).mat[0];

			if (vdotProduct < 0)//if the sphere heads towards the plane
			{
				Matrix vector = p->m_plane.getColumn(1);
				vector = vector - X_t;

				double dotProduct = (vector.transpose()*p->m_plane.getColumn(0)).mat[0];

				if (dotProduct < m_radius)//checking for actual overlaps
				{
					calcCollision(p, this);
				}
			}

			break;
		}
		case RIGID_SPHERE:
		{
			Sphere* s = static_cast<Sphere*>(obj);
			Matrix n12 = s->X_t - X_t;

			if (n12.norm() < (m_radius + s->m_radius))
			{
				normalize(n12);
				Matrix vRelative = s->V_t - V_t;
				double vrel = (vRelative.transpose()*n12).mat[0];
				
				if (vrel < 0)
				calcCollision(this, s, n12,vrel);
			}

			break;
		}
		default:
			break;
	}
}
void Sphere::applyCollisionResponse()
{
	V_t = V_t + applyV_t;//to be updated with J impulse

	setZeros(applyV_t);

}
void Sphere::update(double dt)
{
	X_t = X_t+ V_t*dt;
}
double Sphere::getKinetik()
{
	return 0.5*m_mass*V_t.norm();
}
Matrix Sphere::getMomentum()
{
	return (m_mass*V_t);
}

//////////////////////////////// END SPHERE ////////



/////////	Class PARTICLE	////////////////////////
Particle::Particle(int PartCount,double radius,double totalMass,double cmPosition[],double cmVel[],double W[],double color3[])
:Rigid()
{
	m_id = RIGID_PARTICLE;
	p_particle = this;
	m_mass = totalMass;
	m_radius = radius; //all spheres will have same size and mass

	x_t =Matrix(4, 1, cmPosition);
	v_t=Matrix(4, 1, cmVel);
	
	w_t=Matrix(4,1,W);
	
	rotation = Quaternion(1, 0, 0, 0);
	R_t = Quaternion2RotMatrix(rotation);

	dv = dw = Matrix(4, 1);

	m_color[0] = color3[0]; m_color[15] = color3[1]; m_color[2] = color3[2];




	//radius *=3;

	double point0[] = { radius, radius, radius, 1 };
	double point1[] = { -radius, radius, radius, 1 };

	double point2[] = { -radius, -radius, radius, 1 };
	double point3[] = {  radius, -radius, radius, 1 };	
	
	double point4[] = { radius, radius, -radius, 1 };
	double point5[] = { -radius, radius, -radius, 1 };

	double point6[] = { -radius, -radius, -radius, 1 };
	double point7[] = { radius, -radius, -radius, 1 };


	points = Matrix(4, 1, point0);

	points.addColumn(Matrix(4, 1, point1), 0);
	points.addColumn(Matrix(4, 1, point2), 0);
	points.addColumn(Matrix(4, 1, point3), 0);
	points.addColumn(Matrix(4, 1, point4), 0);
	points.addColumn(Matrix(4, 1, point5), 0);
	points.addColumn(Matrix(4, 1, point6), 0);
	points.addColumn(Matrix(4, 1, point7), 0);


	//radius /= 3;


	/*if ((PartCount % 2)==0)
	{
		double point0[] = { -radius, 0, 0, 1 };
		double point1[] = { +radius, 0, 0, 1 };
		
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
		double point0[] = { 0, 0, 0, 1 };
		points = Matrix(4, 1, point0);

		for (int i = 1; i < 0.5*PartCount; i++)
		{
			point0[0] = -i * 2 * radius;
			points.addColumn(Matrix(4, 1, point0),0);
			point0[0] = +i * 2 * radius;
			points.addColumn(Matrix(4, 1, point0), 0);
		}

	}*/
	

	//// Calculate Inertia and Inverse /////////////
	double Ixx=0, Iyy=0, Izz=0, Ixy=0, Ixz=0, Iyz=0;
	double x , y , z ;
	for (int i = 0; i < points.m_columns; i++)
	{
		x = points.mat[i];
		y = points.mat[i + points.m_columns];
		z = points.mat[i + 2 * points.m_columns];

		Ixx += (y*y + z*z)*m_mass / points.m_columns;
		Iyy += (x*x + z*z)*m_mass / points.m_columns;
		Izz += (x*x + y*y)*m_mass / points.m_columns;

		Ixy += (x*y)*m_mass / points.m_columns;
		Ixz += (x*z)*m_mass / points.m_columns;
		Iyz += (y*z)*m_mass / points.m_columns;

	}

	I = Matrix(4, 4);
	I.mat[0] = Ixx; 	I.mat[1] = -Ixy;	I.mat[2] = -Ixz;	I.mat[3] = 0;
	I.mat[4] = -Ixy;	I.mat[5] =  Iyy;	I.mat[6] = -Iyz;	I.mat[7] = 0;
	I.mat[8] = -Ixz;	I.mat[9] = -Iyz;	I.mat[10] = Izz;	I.mat[11] = 0;
	I.mat[12] = 0;		I.mat[13] = 0;		I.mat[14] = 0;		I.mat[15] = 1;

	I_1 = I.inverse();
	//// End Inertia //////////
	double m = totalMass / 5.0;

	//I.mat[0] = 2.0* m*radius*radius;
	//I.mat[5] = m*radius*radius + 12 * m*radius;
	//I.mat[10] = m*radius*radius + 12 * m*radius;
	//testing shit
	/*setIdentity(I);
	setIdentity(I_1);
	

	I_1 = I.inverse();

	I_1 = I_1 ;*/
	currentPoints = translate(cmPosition)*points;
	cI = I;
	cI_1 = I_1;

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
		glPushMatrix();
		
		//Draw Ы 
		Matrix p2 = x_t + (2*m_radius)*w_t;
		glBegin(GL_LINES);
		glVertex3f(x_t.mat[0], x_t.mat[1], x_t.mat[2] );
		glVertex3f(p2.mat[0], p2.mat[1], p2.mat[2]);

		glEnd();
		glPopMatrix();
	
	glPopAttrib();

}

void Particle::update(double dt)
{
	//Quaternions : q(n+1) = q(n)+dq -> dq/dt= 0.5*w*q -> q=q+(dq/dt)*dt

	//update quaternion Rotation and R_t
	Quaternion qw(0, w_t.mat[0], w_t.mat[1], w_t.mat[2]);
	Quaternion dq = 0.5*qw*rotation*dt;
	rotation = rotation + dq;
	normalize(rotation);
	R_t = Quaternion2RotMatrix(rotation);

	//update the Inertia matrices
	cI = R_t*I*R_t.transpose();
	cI_1 = R_t.transpose()*I_1*R_t;

	//update center mass and points: possition and Velocity, 
	x_t = x_t + v_t*dt;
	currentPoints = translate(x_t.mat[0] , x_t.mat[1] , x_t.mat[2] )*R_t*points;
	
}

void Particle::checkCollision(Rigid* obj)
{
	switch (obj->m_id)
	{
	case RIGID_PLANE:
	{

	}
	default:
		break;
	}
}
void Particle::applyCollisionResponse()
{
	v_t = v_t + dv;
	w_t = w_t + dw;

	setZeros(dv);
	setZeros(dw);

}

double  Particle::getKinetik()
{
	double sumKi = 0;
	for (int i = 0; i < points.m_columns; i++)
	{
		Matrix ri = currentPoints.getColumn(i) - x_t;
		Matrix v = v_t + cross(w_t, ri);
		sumKi += (v.transpose()*v).mat[0];

	}

	sumKi *= 0.5*m_mass/points.m_columns;//common multiplier 1/2 *m of each piece
	return sumKi;
}

Matrix Particle::getMomentum()
{
	Matrix momentum(4, 1);
	for (int i = 0; i < points.m_columns; i++)
	{
		Matrix ri = currentPoints.getColumn(i) - x_t;
		Matrix v = v_t + cross(w_t, ri);
		momentum = momentum + v;
	}
	return momentum*m_mass / points.m_columns;
}

//////////////////////////////// END PARTICLE ////////


///////////////////// GLOBAL ///////////////////////////////////////////

void calcCollision(Plane* plane, Sphere* sphere)
{	
	double nSpeed = (sphere->V_t.transpose()*plane->m_plane.getColumn(0)).mat[0];

	if (nSpeed>0)
	{
		Matrix planeNormal = plane->m_plane.getColumn(0);
		sphere->V_t = sphere->V_t - 2 * nSpeed* plane->m_plane.getColumn(0);
	}
	

}
void calcCollision(Sphere* s1, Sphere* s2,Matrix n12,double vrel)
{	
	double kinetic1 = s1->getKinetik()+s2->getKinetik();

	Matrix J = (-2.0*vrel / (1.0 / s1->m_mass + 1.0 / s2->m_mass)) * n12;

	s1->V_t = s1->V_t -J/s1->m_mass; 
	s2->V_t = s2->V_t +J /s2->m_mass;

	//double error = s1->m_radius + s2->m_radius - r12.norm();
	double kinetic2 = s1->getKinetik() + s2->getKinetik();

	//normalize(r12);

	////s1->X_t = s1->X_t - error*s2->m_radius / (s1->m_radius + s2->m_radius)*r12;
	////s2->X_t = s2->X_t + error*s1->m_radius / (s1->m_radius + s2->m_radius)*r12;

	//double v1prev = (s1->V_t.transpose()*r12).mat[0];
	//double v2prev = (s2->V_t.transpose()*r12).mat[0];

	//
	//	double m1 = s1->m_mass;
	//	double m2 = s2->m_mass;
	//	double invSumMass = 1.0 / (m1 + m2);

	//	double v1after = invSumMass*(v1prev*(m1 - m2) + 2 * m2*v2prev) ;
	//	double v2after = invSumMass*(v2prev*(m2 - m1) + 2 * m1*v1prev) ;

	//	s1->V_t = s1->V_t + (v1after - v1prev)*r12;
	//	s2->V_t = s2->V_t + (v2after - v2prev)*r12;

	
	
}


void calcCollision(Plane* plane, Particle* particle, Matrix cp/*ContactPoint*/,
													Matrix  rp/*RelativePosFromCenterMass*/,
													Matrix vp /*VelocityofCp*/,
													Matrix n /*collision normal */,
													double vr/*Relative Velocity on the normal*/)
{
	//	J= - (vp*n) / (1/m1  + n*(invInertia1*(rp x n) x rp	)


}


