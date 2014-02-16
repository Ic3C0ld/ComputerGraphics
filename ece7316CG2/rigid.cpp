#include "rigid.h"

////extern variables and flags
int sphere_2_wallCollisions = 0;
int particle_2_wallCollisions = 0;
int spring_2_wallCollisions = 0;
double simTime = 0;

std::vector<double>data_sphere_2_wallCollisions;
std::vector<double>data_particle_2_wallCollisions;
std::vector<double>data_spring_2_wallCollisions;
std::vector<double>data_time;



int collisions_2_everyone_else = 0;
std::vector<double>data_collisions_2_everyone_else;


int followUpObjID = 6;



 bool pA_collide_no_speedchange = true;
 bool pA_collide_with_speedchange = true;
 bool pA_collisionstats = true;

//PartB
 bool pB1_use_material_properties = true;
 bool pB2_use_particles_no_rotation = false;
 bool pB2_use_particles_with_rotation = true;

 bool pB3_use_followUp_camera = false;

 bool pb4_use_springs_with_speedchange = true;



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
		Matrix n = -1 * m_plane.getColumn(0); //"normal" , so that it points inwards the box
		Matrix p0 = m_plane.getColumn(1);		//"point0", get one plane point
		int count = p->currentPoints.m_columns;// "#parts" , consisting the particle body

		
		Matrix rPavg(4,1);
		int colPoints = 0;
		for (int i = 0; i < count; i++)
		{
			Matrix tP = p->currentPoints.getColumn(i); //"testPoint" , current body part to test

			//check for plane collision
			double distance = ((tP - p0).transpose()*n).mat[0];

			if (distance < p->m_radius) 
			{//we have intersection or have gone completely out of box
				
				Matrix cP = tP - p->m_radius*n;			   //"collisionPoint" , current collision point 
				Matrix rP = cP - p->x_t;					//collision point with respect to the center mass of the body
				Matrix vP = p->v_t + cross(p->w_t, rP);		// vTotal= vLinear + w_t x r;


				double  vRel = -1*(vP.transpose()*n).mat[0];	//vRelative = vPlane - Vpoint = 0 -vP -> dot with normal

				if ( vRel > 0)
				{//we have a collision
					colPoints++;
					rPavg = (rPavg + rP);
					//calcCollision(this, p,  rP,  n, vRel);
					
				}


			}	
		}

		if (colPoints>0)
		{
			rPavg = rPavg /colPoints;

			Matrix vP = p->v_t + cross(p->w_t, rPavg);
			double  vRel = -1 * (vP.transpose()*n).mat[0];
			calcCollision(this, p, rPavg, n, vRel);

		}

		break;
	}
	case RIGID_SPRING:
	{
		Plane* p= this;
		SpringSystem* s = static_cast<SpringSystem*>(obj);


		Matrix pPoint = p->m_plane.getColumn(1);
		Matrix n = -1 * p->m_plane.getColumn(0);//vector looking inwards

		Matrix r01 = s->x1_t - pPoint;
		Matrix r02 = s->x2_t - pPoint;

		if ((r01.transpose()*n).mat[0] < s->m_radius1)
		{//intersection with sphere 1

			double vRel = ((-1 * s->v1_t).transpose()*n).mat[0];

			if (vRel>0)
			{//collision with sphere 1

				Matrix j = (-2 * vRel) / (1 / s->m_mass1)*n;
				s->J1 = s->J1 - j;
			}

		}

		if ((r02.transpose()*n).mat[0] < s->m_radius2)
		{//intersection with sphere 2

			double vRel = ((-1 * s->v2_t).transpose()*n).mat[0];

			if (vRel>0)
			{//collision with sphere 2

				Matrix j = (-2 * vRel) / (1 / s->m_mass2)*n;
				s->J2 = s->J2 - j;
			}

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

	Jtotal = Matrix(4, 1);//zero-ed
}
void Sphere::draw()
{
	if (pB1_use_material_properties == true)
	{
		glDisable(GL_COLOR_MATERIAL);
		GLfloat material_diffuse[] = { m_color[0], m_color[1], m_color[2], 1 };
		GLfloat material_specular[] = { m_color[0], m_color[1], m_color[2], 1 };
		GLfloat material_shininess[] = {30*( m_color[0]+ m_color[1]+ m_color[2] )};
		glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, material_shininess);

	}
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(m_color[0], m_color[1], m_color[2]);

	glPushMatrix();
		glTranslatef(X_t.mat[0], X_t.mat[1], X_t.mat[2]);
		glutSolidSphere(m_radius, 15+m_radius, 15+m_radius);
	glPopMatrix();

	glPopAttrib();

	if (pB1_use_material_properties == true)		glEnable(GL_COLOR_MATERIAL);


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
			Sphere* s1 = this;
			Sphere* s2 = static_cast<Sphere*>(obj);
			Matrix n12 = s2->X_t -s1-> X_t;

			if (n12.norm() < (s1->m_radius + s2->m_radius))
			{
				normalize(n12);
				Matrix vRelative = s1->V_t - s2->V_t;
				double vrel = (vRelative.transpose()*n12).mat[0];

				if (vrel > 0)
				{
					if (pA_collide_with_speedchange == true)
					{
						calcCollision(s1, s2, n12, vrel);
					}
					else if (pA_collide_no_speedchange)
					{
						double v1dot = (s1->V_t.transpose()*n12).mat[0];
						double v2dot = (s2->V_t.transpose()*n12).mat[0];

						s1->V_t = s1->V_t - 2 * v1dot*n12;
						s2->V_t = s2->V_t - 2 * v2dot*n12;

					}

				}
			}
				

			
			break;
		}
		case RIGID_PARTICLE:
		{
			Particle* p = static_cast<Particle*>(obj);
			Sphere* s = this;

			int pPoints = p->currentPoints.m_columns;

			int colPoints = 0;
			Matrix r1avg(4, 1);
			Matrix r2avg(4, 1);
			Matrix navg(4, 1);

			for (int i = 0; i < pPoints; i++)
			{
				Matrix tp = p->currentPoints.getColumn(i);

				Matrix r12 = tp - s->X_t;
				Matrix n = r12 / r12.norm();

				if (r12.norm() < s->m_radius + p->m_radius)
				{//intersection

					Matrix cP = s->X_t + r12*(s->m_radius) / (s->m_radius + p->m_radius);

					Matrix r1 = cP - s->X_t;
					Matrix r2 = cP - p->x_t;

					Matrix v1 = s->V_t;
					Matrix v2 = p->v_t + cross(p->w_t, r2);
					
					double vRel = ((v1-v2).transpose()*n).mat[0];

					if (vRel > 0)
					{//collision

						r1avg = r1avg + r1;
						r2avg = r2avg + r2;
						navg = navg + n;

						colPoints++;

						i = pPoints;//break after the first hit

					}//END if collision


				}//END IF intersection

			}//END for (i:pPoints)

			if (colPoints > 0)
			{
				
					r1avg = r1avg / colPoints;
					r2avg = r2avg / colPoints;
					navg = navg / colPoints;
					normalize(navg);

				//v1 is condidered the same as i dont allow any rotation for the spheres
				Matrix v2 = p->v_t + cross(p->w_t, r2avg);

				double vRel = ((s->V_t - v2).transpose()*navg).mat[0];

				calcCollision(s, p,r1avg,r2avg, navg, vRel);
			}
		}
			break;
		case RIGID_SPRING:
		{
			SpringSystem* spring = static_cast<SpringSystem*>(obj);
			Sphere* sphere = this;

			Matrix r01 = spring->x1_t - sphere->X_t;
			Matrix r02 = spring->x2_t - sphere->X_t;

			Matrix n01 = r01 / r01.norm();
			Matrix n02 = r02 / r02.norm();


			if (r01.norm() < sphere->m_radius + spring->m_radius1)
			{//intersection with first sphere of spring

				double vRel1 = ((sphere->V_t - spring->v1_t).transpose()*n01).mat[0];

				if (vRel1>0)
				{//collision

					Matrix j = (-2 * vRel1) / (1 / sphere->m_mass + 1 / spring->m_mass1)*n01;
					sphere->Jtotal = sphere->Jtotal + j;
					spring->J1 = spring->J1 - j;

					data_collisions_2_everyone_else.push_back(++collisions_2_everyone_else);

				}
			}

			if (r02.norm() < sphere->m_radius + spring->m_radius2)
			{//intersection with second sphere of spring

				double vRel2 = ((sphere->V_t - spring->v2_t).transpose()*n01).mat[0];

				if (vRel2>0)
				{//collision

					Matrix j = (-2 * vRel2) / (1 / sphere->m_mass + 1 / spring->m_mass2)*n01;
					sphere->Jtotal = sphere->Jtotal + j;
					spring->J2 = spring->J2 - j;

					data_collisions_2_everyone_else.push_back(++collisions_2_everyone_else);

				}
			}

			break;

		}
		default:
			break;
	}
}
void Sphere::applyCollisionResponse()
{
	V_t = V_t + Jtotal/m_mass;//to be updated with J impulse

	setZeros(Jtotal);

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

	m_color[0] = color3[0]; m_color[15] = color3[1]; m_color[2] = color3[2];


	radius *=0.5;

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


	radius *=2;


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
	//// spheres' masses  are considered as point masses ////
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
		
		//Draw ù 
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
		//this method will not be called since will planes are first on the object 
		break;
	}
	case RIGID_SPHERE:
	{
		break;
	}
	case RIGID_PARTICLE:
	{
		Particle *p1 = this; ////for legibility
		Particle *p2 = static_cast<Particle*>(obj);

		int p1Points = p1->currentPoints.m_columns;
		int p2Points = p2->currentPoints.m_columns;


		int colPoints=0;
		Matrix rP1avg(4,1);
		Matrix rP2avg(4,1);
		Matrix navg(4,1);

		for (int i = 0; i < p1Points; i++)
		{
			for (int j = 0; j < p2Points; j++)
			{
				Matrix tP1 = p1->currentPoints.getColumn(i);
				Matrix tP2 = p2->currentPoints.getColumn(j);

				Matrix r12 = tP2 - tP1;
				Matrix n = r12 / r12.norm();
				if (r12.norm() < (p1->m_radius + p2->m_radius))
				{//intersection
					Matrix cP = tP1+r12*(p1->m_radius) / (p1->m_radius + p1->m_radius); ////not completely accurate 

					Matrix rp1 = cP - p1->x_t;
					Matrix rp2 = cP - p2->x_t;

					Matrix v1 = p1->v_t + cross(p1->w_t, rp1);
					Matrix v2 = p2->v_t + cross(p2->w_t, rp2);

					double vRel = ((v1 - v2).transpose()*n).mat[0];

					if (vRel > 0)
					{//collision
						rP1avg = rP1avg + rp1;
						rP2avg = rP2avg + rp2;
						navg = navg + n;
						colPoints++;
						
						//i = j = 100; // break after the first collision

					}
				}
			}//FOR j:p2Points
		}//FOR i:p1Points

		if (colPoints > 0)
		{
	
			rP1avg = rP1avg / colPoints;
			rP2avg = rP2avg / colPoints;
			normalize(navg);

			Matrix vP1 = p1->v_t + cross(p1->w_t, rP1avg);
			Matrix vP2 = p2->v_t + cross(p2->w_t, rP2avg);
			double vRel = ((vP1 - vP2).transpose()*navg).mat[0];

			calcCollision(p1, p2, rP1avg, rP2avg, navg, vRel);

		}

		break;
	}
	case RIGID_SPRING:
	{
		Particle* p = this;
		SpringSystem* s = static_cast<SpringSystem*>(obj);


		int pPoints = p->currentPoints.m_columns;

		int col1Points = 0;
		int col2Points = 0;

		Matrix r12avg(4, 1);	//Spring rAvgs	
		Matrix r22avg(4, 1);
		Matrix rb1avg(4, 1);	//particle rAvg
		Matrix rb2avg(4, 1);

		Matrix navg1(4,1), navg2(4,1);

		for (int i = 0; i < pPoints; i++)
		{
			Matrix tp = p->currentPoints.getColumn(i);

			Matrix r12 = tp-s->x1_t ;
			Matrix r22 = tp-s->x2_t ;

			Matrix n1 = r12 / r12.norm();
			Matrix n2 = r22 / r22.norm();

			if (r12.norm() < s->m_radius1 + p->m_radius)
			{//intersection with the first spring ball

				Matrix cP = s->x1_t + n1*s->m_radius1;

				Matrix r1 = cP - s->x1_t;
				Matrix r2 = cP - p->x_t;

				Matrix v1 = s->v1_t;
				Matrix v2 = p->v_t + cross(p->w_t, r2);

				double vRel = ((v1 - v2).transpose()*n1).mat[0];

				if (vRel > 0)
				{//collision with the first spring ball

					r12avg = r12avg + r1;
					rb1avg = rb1avg + r2;
					navg1 = navg1 + n1;

					col1Points++;
					//i = pPoints;//dont allow other collisions ,break the loop
				}//endif collision
			}//ENDIF intersection

			if (r22.norm() < s->m_radius1 + p->m_radius)
			{//intersection with the first spring ball

				Matrix cP = s->x2_t + n2*s->m_radius1;

				Matrix r1 = cP - s->x2_t;
				Matrix r2 = cP - p->x_t;

				Matrix v1 = s->v2_t;
				Matrix v2 = p->v_t + cross(p->w_t, r2);

				double vRel = ((v1 - v2).transpose()*n2).mat[0];

				if (vRel > 0)
				{//collision with the first spring ball

					r22avg = r22avg + r1;
					rb2avg = rb2avg + r2;
					navg2 = navg2 + n2;

					col2Points++;
					//i = pPoints;//dont allow other collisions ,break the loop

				}//endif vRel1>0

			}//ENDIF intersection

		}//END for(i:particle.points)


		if (col1Points > 0)
		{
			r12avg = r12avg / col1Points;
			rb1avg = rb1avg / col1Points;
			navg1 = navg1 / col1Points;
			normalize(navg1);

			//v1=s->v1_t;
			Matrix v2 = p->v_t + cross(p->w_t, rb1avg);

			double vRel = ((s->v1_t - v2).transpose()*navg1).mat[0];

			Matrix j = -2 * vRel*navg1 / (1 / s->m_mass1 + 1 / p->m_mass + (navg1.transpose()*p->cI_1*cross(cross(rb1avg, navg1), rb1avg)).mat[0]);

			s->J1 = s->J1 + j;

			p->J.push_back(-1 * j);
			p->rP.push_back(rb1avg);

		}
		if (col2Points > 0)
		{
			r22avg = r22avg / col2Points;
			rb2avg = rb2avg / col2Points;
			navg2 = navg2 / col2Points;
			normalize(navg2);

			//v1=s->v2_t;
			Matrix v2 = p->v_t + cross(p->w_t, rb2avg);

			double vRel = ((s->v2_t - v2).transpose()*navg2).mat[0];

			Matrix j = -2 * vRel*navg2 / (1 / s->m_mass2 + 1 / p->m_mass + (navg2.transpose()*p->cI_1*cross(cross(rb2avg, navg2), rb2avg)).mat[0]);

			s->J2 = s->J2 + j;

			p->J.push_back(-1 * j);
			p->rP.push_back(rb1avg);

		}


	}
	default:
		break;
	}
}
void Particle::applyCollisionResponse()
{
	for (int i = 0; i < J.size(); i++)
	{
		v_t = v_t + J[i] / m_mass;

		if (pB2_use_particles_with_rotation==true)
		w_t = w_t + cI_1*cross(rP[i], J[i]);
	}
	
	
	J.clear();
	rP.clear();

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






/////////	Class SPRING	////////////////////////

SpringSystem::SpringSystem(Plane* p, double xPercent,double yPercent,
									 double k1, double k2,
									 double mass_1, double mass_2,
									 double radius_1, double radius_2,
									 double x1t_4_1[], double x2t_4_1[],
									 double v1t_4_1[], double v2t_4_1[],
									 double color[])
:Rigid()
{
	m_id = RIGID_SPRING;
	p_spring = this; 

	m_color[0] = color[0];	m_color[1] = color[1];	m_color[2] = color[2];
	m_k1 = k1;				m_k2 = k2;
	m_mass1 = mass_1;		m_mass2 = mass_2;
	m_radius1 = radius_1;	m_radius2 = radius_2;

	x1_t = Matrix(4, 1, x1t_4_1);	x2_t = Matrix(4, 1, x2t_4_1);
	v1_t = Matrix(4, 1, v1t_4_1);	v2_t = Matrix(4, 1, v2t_4_1);

	Matrix xV = p->m_plane.getColumn(3) - p->m_plane.getColumn(2);  //collumn(0)=normalVector,collumn(1)=centerPoint,the next ones are the side points
	Matrix yV = p->m_plane.getColumn(4) - p->m_plane.getColumn(3);

	x0 = p->m_plane.getColumn(2) + 0.01*xPercent*xV + 0.01*yPercent*yV;

	v1_t = v2_t = a1_t=a2_t=Matrix(4, 1); //zero

	double downMat[] = { 0, -1, 0, 0 };
	Matrix down(4, 1,downMat);
	x1_t = x0 - 2 * m_radius1*down;
	x1_t = x1_t - 2 * m_radius2*down;

	J1 = J2 = Matrix(4, 1);


	xPer = xPercent;
	yPer = yPercent;

	plane = p;
	
}


void SpringSystem::draw()
{
	glPushAttrib(GL_COLOR_BUFFER_BIT);
	glColor3f(m_color[0], m_color[1], m_color[2]);
	glPushMatrix();

	//draw Spheres

	glPushMatrix();
	glTranslatef(x1_t.mat[0], x1_t.mat[1], x1_t.mat[2]);
	glutSolidSphere(m_radius1,20,20);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(x2_t.mat[0], x2_t.mat[1], x2_t.mat[2]);
	glutSolidSphere(m_radius2, 20, 20);
	glPopMatrix();

	//draw lines for Springs
	
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	glVertex3f(x0.mat[0], x0.mat[1], x0.mat[2]);
	glColor3f(1,0,0);

	glVertex3f(x1_t.mat[0], x1_t.mat[1], x1_t.mat[2]);
	glColor3f(0,1,0);

	glVertex3f(x2_t.mat[0], x2_t.mat[1], x2_t.mat[2]);
	glEnd();
	glLineWidth(1);


	glPopMatrix();
	glPopMatrix();

}


void SpringSystem::update(double dt)
{
	
	//// calculate forces 
	// gravity and vectors
	Matrix g(4, 1);
	g.mat[1] = -9.81*10;//hack to make gravity feel like one ,,dont know why it behaves like the moon instead of the earth ,maybe the NOT real time calculations;? FIX PENDING

	Matrix r01 = x1_t - x0;
	Matrix n01 = r01 / r01.norm();

	Matrix r12 = x2_t - x1_t;
	Matrix n12 = r12 / r12.norm();

	//Spring Force1
	Matrix F21 = (m_k2*(2 * m_radius1 - r12.norm()) /*+ 0.1*m_k2*pow((2 * m_radius1 - r12.norm()),3)*/)*n12;
	Matrix F10 = (m_k1*(2 * m_radius2 - r01.norm()) /*+ 0.1*m_k1*pow(2 * m_radius2 - r01.norm(),3)*/)*n01;


	////CheckCollision with itself
	if (r12.norm() < m_radius1 + m_radius2)
	{//intersection

		double vRel = ((v1_t - v2_t).transpose()*n12).mat[0];
		if (vRel>0)
		{//colliding

			Matrix j = (-2 * vRel)/(1/m_mass1+1/m_mass2)*n12 ;
			J1 = J1 + j;
			J2 = J2 - j;
		}
	}

	v1_t = v1_t + J1 / m_mass1;
	v2_t = v2_t + J2 / m_mass2;


	Matrix vrel01 = v1_t;
	Matrix vrel12 = v2_t-v1_t;


	////Update variables

	x1_t = x1_t + v1_t*dt;
	x2_t = x2_t + v2_t*dt;

	v1_t = v1_t + a1_t*dt;
	v2_t = v2_t + a2_t*dt;

	a1_t = ((m_mass1 + m_mass2)*g + F10 - F21 - 3* vrel01) / m_mass1;
	a2_t = (m_mass2*g + F21 - 3* vrel12) / m_mass2;



	setZeros(J1);
	setZeros(J2);


}

void SpringSystem::checkCollision(Rigid* obj)
{
	switch (obj->m_id)
	{
	case RIGID_PLANE:
	{
		SpringSystem* s = this;
		Plane* p = static_cast<Plane*>(obj);

		Matrix pPoint = p->m_plane.getColumn(1);
		Matrix n = -1 * p->m_plane.getColumn(0);//vector looking inwards

		Matrix r01 = s->x1_t - pPoint;
		Matrix r02 = s->x2_t - pPoint;

		if (r01.norm() < s->m_radius1)
		{//intersection with sphere 1

			double vRel = ((-1 * s->v1_t).transpose()*n).mat[0];
			
			if (vRel>0)
			{//collision with sphere 1

				Matrix j = (-2 * vRel) / (1 / s->m_mass1)*n;
				s->J1 = s->J1 - j;
			}

		}

		if (r02.norm() < s->m_radius2)
		{//intersection with sphere 2

			double vRel = ((-1 * s->v2_t).transpose()*n).mat[0];

			if (vRel>0)
			{//collision with sphere 2

				Matrix j = (-2 * vRel) / (1 / s->m_mass2)*n;
				s->J2 = s->J2 - j;
			}

		}



		break;
	}
	case RIGID_SPHERE:
	{

		Sphere* sphere = static_cast<Sphere*>(obj);
		SpringSystem* spring = this;

		Matrix r01 = spring->x1_t - sphere->X_t;
		Matrix r02 = spring->x2_t - sphere->X_t;

		Matrix n01 = r01 / r01.norm();
		Matrix n02 = r02 / r02.norm();


		if (r01.norm() < sphere->m_radius + spring->m_radius1)
		{//intersection with first sphere of spring

			double vRel1 = ((sphere->V_t - spring->v1_t).transpose()*n01).mat[0];

			if (vRel1>0)
			{//collision

				Matrix j = (-2 * vRel1) / (1 / sphere->m_mass + 1 / spring->m_mass1)*n01;
				sphere->Jtotal = sphere->Jtotal + j;
				spring->J1 = spring->J1 - j;
			}
		}

		if (r02.norm() < sphere->m_radius + spring->m_radius2)
		{//intersection with second sphere of spring

			double vRel2 = ((sphere->V_t - spring->v2_t).transpose()*n01).mat[0];

			if (vRel2>0)
			{//collision

				Matrix j = (-2 * vRel2) / (1 / sphere->m_mass + 1 / spring->m_mass2)*n01;
				sphere->Jtotal = sphere->Jtotal + j;
				spring->J2 = spring->J2 - j;
			}
		}

		break;
	}

	case RIGID_PARTICLE:
	{
		Particle* p = static_cast<Particle*>(obj);
		SpringSystem* s = this;


		int pPoints = p->currentPoints.m_columns;

		int col1Points = 0;
		int col2Points = 0;

		Matrix r12avg(4,1);	//Spring rAvgs	
		Matrix r22avg(4, 1);
		Matrix rb1avg(4, 1);	//particle rAvg
		Matrix rb2avg(4, 1);

		Matrix navg1, navg2;

		for (int i = 0; i < pPoints; i++)
		{
			Matrix tp = p->currentPoints.getColumn(i);

			Matrix r12 = s->x1_t - tp;
			Matrix r22 = s->x2_t - tp;

			Matrix n1 = r12 / r12.norm();
			Matrix n2 = r22 / r22.norm();

			if (r12.norm() < s->m_radius1 + p->m_radius)
			{//intersection with the first spring ball

				Matrix cP = s->x1_t + n1*s->m_radius1;

				Matrix r1 = cP - s->x1_t;
				Matrix r2 = cP - p->x_t;

				Matrix v1 = s->v1_t;
				Matrix v2 = p->v_t + cross(p->w_t, r2);

				double vRel = ((v1 - v2).transpose()*n1).mat[0];

				if (vRel > 0)
				{//collision with the first spring ball

					r12avg = r12avg + r1;
					rb1avg = rb1avg + r2;
					navg1 = navg1 + n1;

					col1Points++;
				}//endif collision
			}//ENDIF intersection

			if (r22.norm() < s->m_radius1 + p->m_radius)
			{//intersection with the first spring ball

				Matrix cP = s->x2_t + n2*s->m_radius1;

				Matrix r1 = cP - s->x2_t;
				Matrix r2 = cP - p->x_t;

				Matrix v1 = s->v2_t;
				Matrix v2 = p->v_t + cross(p->w_t, r2);

				double vRel = ((v1 - v2).transpose()*n2).mat[0];

				if (vRel > 0)
				{//collision with the first spring ball

					r22avg = r22avg + r1;
					rb2avg = rb2avg + r2;
					navg2 = navg2 + n2;

					col2Points++;
				}//endif vRel1>0

			}//ENDIF intersection

		}//END for(i:particle.points)


		if (col1Points > 0)
		{
			r12avg = r12avg / col1Points;
			rb1avg = rb1avg / col1Points;
			navg1 = navg1 / col1Points;
			normalize(navg1);

			//v1=s->v1_t;
			Matrix v2 = p->v_t + cross(p->w_t, rb1avg);

			double vRel = ((s->v1_t - v2).transpose()*navg1).mat[0];

			Matrix j = -2 * vRel*navg1 / (1 / s->m_mass1 + 1 / p->m_mass + (navg1.transpose()*p->cI_1*cross(cross(rb1avg, navg1), rb1avg)).mat[0]);

			s->J1 = s->J1 + j;

			p->J.push_back(-1 * j);
			p->rP.push_back(rb1avg);
			
		}
		if (col2Points > 0)
		{
			r22avg = r22avg / col1Points;
			rb2avg = rb2avg / col1Points;
			navg2 = navg2 / col2Points;
			normalize(navg2);

			//v1=s->v2_t;
			Matrix v2 = p->v_t + cross(p->w_t, rb2avg);

			double vRel = ((s->v2_t - v2).transpose()*navg2).mat[0];

			Matrix j = -2 * vRel*navg2 / (1 / s->m_mass2 + 1 / p->m_mass + (navg2.transpose()*p->cI_1*cross(cross(rb2avg, navg2), rb2avg)).mat[0]);

			s->J2 = s->J2 + j;

			p->J.push_back(-1 * j);
			p->rP.push_back(rb1avg);

		}


	}
	case RIGID_SPRING:
	{
		SpringSystem* s1=static_cast<SpringSystem*>(obj);
		SpringSystem* s2 = this;

		Matrix r11 = s2->x1_t - s1->x1_t;
		Matrix r12 = s2->x2_t - s1->x1_t;

		Matrix r21 = s2->x1_t - s1->x2_t;
		Matrix r22 = s2->x2_t - s1->x2_t;

		Matrix n11 = r11 / r11.norm();
		Matrix n12 = r12 / r12.norm();
		Matrix n21 = r21 / r21.norm();
		Matrix n22 = r22 / r22.norm();


		if (r11.norm() < s1->m_radius1 + s2->m_radius1)
		{//intersection with first sphere of spring
			double vRel1 = ((s1->v1_t - s2->v1_t).transpose()*n11).mat[0];

			if (vRel1>0)
			{//collision
				Matrix j = (-2 * vRel1) / (1 / m_mass1 + 1 / m_mass2)*n11;
				s1->J1 = s1->J1 + j;
				s2->J1 = s2->J1 - j;
			}
		}

		if (r12.norm() < s1->m_radius1 + s2->m_radius1)
		{//intersection with first sphere of spring

			double vRel12 = ((s1->v1_t - s2->v2_t).transpose()*n12).mat[0];

			if (vRel12>0)
			{//collision

				Matrix j = (-2 * vRel12) / (1 / m_mass1 + 1 / m_mass2)*n12;
				s1->J1 = s1->J1 + j;
				s2->J2 = s2->J2 - j;

			}
		}

		if (r21.norm() < s1->m_radius1 + s2->m_radius1)
		{//intersection with first sphere of spring

			double vRel21 = ((s1->v2_t - s2->v1_t).transpose()*n21).mat[0];

			if (vRel21>0)
			{//collision

				Matrix j = (-2 * vRel21) / (1 / m_mass1 + 1 / m_mass2)*n21;
				s1->J2 = s1->J2 + j;
				s2->J1 = s2->J1 - j;

			}
		}

		if (r22.norm() < s1->m_radius1 + s2->m_radius1)
		{//intersection with first sphere of spring

			double vRel22 = ((s1->v2_t - s2->v2_t).transpose()*n22).mat[0];

			if (vRel22>0)
			{//collision

				Matrix j = (-2 * vRel22) / (1 / m_mass1 + 1 / m_mass2)*n22;
				s1->J2 = s1->J2 + j;
				s2->J2 = s2->J2 - j;

			}
		}


		break;
	}
	default:
		break;
	}
}
void SpringSystem::applyCollisionResponse(){}

double  SpringSystem::getKinetik(){ return 0; }

Matrix SpringSystem::getMomentum(){ return Matrix(4, 1); }
//////////////////////////////// END SPRING ////////

///////////////////// GLOBAL ///////////////////////////////////////////

void calcCollision(Plane* plane, Sphere* sphere)
{	
	double nSpeed = (sphere->V_t.transpose()*plane->m_plane.getColumn(0)).mat[0];

	if (nSpeed>0)
	{
		Matrix planeNormal = plane->m_plane.getColumn(0);
		sphere->V_t = sphere->V_t - 2 * nSpeed* plane->m_plane.getColumn(0);
	}
	
	
	data_sphere_2_wallCollisions.push_back(++sphere_2_wallCollisions);

}
void calcCollision(Sphere* s1, Sphere* s2,Matrix n12,double vrel)
{	
	double kinetic1 = s1->getKinetik()+s2->getKinetik();

	Matrix J = (-2.0*vrel / (1.0 / s1->m_mass + 1.0 / s2->m_mass)) * n12;


	
	s1->Jtotal = s1->Jtotal + J;
	s2->Jtotal = s2->Jtotal - J;

	/*s1->V_t = s1->V_t +J/s1->m_mass; 
	s2->V_t = s2->V_t -J /s2->m_mass;*/

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

	data_collisions_2_everyone_else.push_back(++collisions_2_everyone_else);
	
}


void calcCollision(Plane* , Particle* p, //Matrix cp/*ContactPoint*/,
													Matrix&  rp/*RelativePosFromCenterMass*/,
												//	Matrix vp /*VelocityofCp*/,
													Matrix& n /*collision normal */,
													double  vRel/*Relative Velocity on the normal*/)
{
	//	J= -2*(vp*n) / (1/m1  + n*(invInertia1*(rp x n) x rp	)


	double J = -2 * vRel / (1 / p->m_mass + (n.transpose()*(p->cI_1*cross(cross(rp, n), rp))).mat[0]);

	p->J.push_back(-1 * J*n );
	p->rP.push_back(rp);

	data_particle_2_wallCollisions.push_back(++particle_2_wallCollisions);
}


void calcCollision(Particle* p1, Particle* p2, Matrix& r1, Matrix& r2, Matrix& n, double vRel)
{
	double nom = -2 * vRel;
	double den1 = 1 / p1->m_mass + 1/p2->m_mass;
	double den2 = (n.transpose()*(p1->cI_1*cross(cross(r1, n), r1))).mat[0];
	double den3 = (n.transpose()*(p2->cI_1*cross(cross(r2, n), r2))).mat[0];


	double J = nom / (den1 + den2 + den3);

	p1->J.push_back(J*n);
	p2->J.push_back(-1 * J*n);

	p1->rP.push_back(r1);
	p2->rP.push_back(r2);


	data_collisions_2_everyone_else.push_back(++collisions_2_everyone_else);

}


void calcCollision(Sphere* s, Particle* p, Matrix& r1, Matrix& r2, Matrix& n, double vRel)
{
	double nom = -2 * vRel;

	double den1 = 1 / s->m_mass + 1 / p->m_mass;
	double den2 = (n.transpose()*(p->cI_1*cross(cross(r2, n), r2))).mat[0];

	double J = nom / (den1 + den2);

	p->J.push_back(-1*J*n);
	p->rP.push_back(r2);

	s->Jtotal = s->Jtotal + J*n;

	data_collisions_2_everyone_else.push_back(++collisions_2_everyone_else);

}