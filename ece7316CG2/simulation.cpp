#include "simulation.h"


extern double totalKinetic;

//Class Simulation//

Simulation::Simulation(double boxSize, int spheres, int particles, int springs)
{
	sim_dt = target_dt = 0.001; //for now
	
	
//// BOX: Create planes  ////////////////////////

//TOP green
	double colorG[] = { 0.15, 0.3,0.15};
	Plane *top = new Plane(boxSize,colorG);
	top->m_plane = translate(0, boxSize / 2, 0)*top->m_plane;

//FRONT dark grey
	double colorDGrey[] = { 0.3, 0.3, 0.3 };
	Plane *front = new Plane(boxSize, colorDGrey);
	front->m_plane = rotX(M_PI_2)* top->m_plane;
 
//BOTTOM white
	double colorR[] = { 1, 1, 1 };
	Plane *bottom = new Plane(boxSize,colorR);
	bottom->m_plane = rotX(M_PI_2)*front->m_plane;

//BACK dark grey
	Plane *back = new Plane(boxSize, colorDGrey);
	back->m_plane = rotX(M_PI_2)*bottom->m_plane;

//LEFT light grey
	double colorGrey[] = { 0.5, 0.5, 0.5 };
	Plane *left = new Plane(boxSize, colorGrey);
	left->m_plane = rotY(M_PI_2)*front->m_plane;

//Right light grey
	Plane *right = new Plane(boxSize, colorGrey);
	right->m_plane = rotY(-M_PI_2)*front->m_plane;


	objects.push_back(top);
	objects.push_back(bottom);
	objects.push_back(front);
	objects.push_back(back);
	objects.push_back(left);
	objects.push_back(right);

	//translate the whole box to sit on the ground, not to have 0,0,0 as center
	for (int i = 0; i < 6; i++)
	{
		objects[i]->p_plane->m_plane = translate(0, boxSize, 0)*rotX(30,true)*objects[i]->p_plane->m_plane;
	}

//////////////////////////// BOX READY //////////////

	srand(time(NULL));
/////	SPHERES		////////////////
	for (int i = 0; i < spheres;i++)
	{


		double radius = 0.02*(rand() % 100) + 1;
		double mass = 10*pow(radius, 3);

		double x_t[] = { (rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius), (rand() % (int)(boxSize - 2*radius))+radius, (rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius) };
		double v_t[] = { (rand() % 200) - 10, (rand() % 200) - 10, (rand() % 200) - 10 };

		double color3[] = { 0.1*(rand() % 10), 0.1*(rand() % 10), 0.1*(rand() % 10) };

		Sphere *temp = new Sphere(radius, mass, x_t, v_t, color3);
		
		
		objects.push_back(temp);
	}



	/*for (int i = 0; i < 6; i++)
	{
		objects[i]->p_plane->m_plane = translate(0,10,0)*rotX(45, true)* objects[i]->p_plane->m_plane;
	}*/
	//////////////////////////// Spheres READY //////////////


	/////	PARTICLES		////////////////
	for (int i = 0; i < particles; i++)
	{
		double radius = 0.02*(rand() % 100) + 0.3;
		double mass = 10 * pow(radius, 3);


		double x_t[] = { (rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius),
									(rand() % (int)(boxSize - 2 * radius)) + radius	,
						(rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius)	,
																				1 };
		double v_t[] = { (rand() % 200) - 10, (rand() % 200) - 10, (rand() % 200) - 10 ,0};
		double w_t[] = { 0.1*(rand() % 200) - 10, 0.1*(rand() % 200) - 10, 0.1*(rand() % 200) - 10 ,0};
		/*double v_t[] = { -30,0,0,0};

		double w_t[] = {0, 10, 0,0};

		double x_t[] = { (rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius),
									15,
						(rand() % (int)(boxSize - radius)) - 0.5*(boxSize - radius)	,
																				1 };*/
		double color3[] = { 0.1*(rand() % 10), 0.1*(rand() % 10), 0.1*(rand() % 10) };

		Particle *p = new Particle(5, radius, mass, x_t, v_t, w_t,color3);
		objects.push_back(p);

	}
	//////////////////////////// Particles READY //////////////


	//// SPRINGS //////////////////////


	for (int i = 0; i < springs; i++)
	{
		double xPercent = (rand() % 80) + 10;
		double yPercent = (rand() % 80) + 10;

		double color3[] = { 0.1*(rand() % 10), 0.1*(rand() % 10), 0.1*(rand() % 10) };


		double radius1 = 0.03*(rand() % 100) + 0.15;
		double mass1 = 4* pow(radius1, 3);
		double radius2 = 0.03*(rand() % 100) + 0.15;
		double mass2 = 4 * pow(radius2, 3);

		double k1 = rand() %100 + 60;
		double k2= rand() % 80 + 40;

	    double x_t[] = { (rand() % (int)(boxSize - radius1)) - 0.5*(boxSize - radius1),
									(rand() % (int)(boxSize - 2 * radius1)) + radius1	,
						(rand() % (int)(boxSize - radius1)) - 0.5*(boxSize - radius1)	,
																				1 };
		 double x2_t[] = { (rand() % (int)(boxSize - radius1)) - 0.5*(boxSize - radius1),
									(rand() % (int)(boxSize - 2 * radius1)) + radius1	,
						(rand() % (int)(boxSize - radius1)) - 0.5*(boxSize - radius1)	,
																				1 };
		 double v[] = { 0, 0, 0, 0 };

		 SpringSystem* s=new SpringSystem(top, xPercent,yPercent,k1,k2,mass1,mass2,radius1,radius2,x_t,x2_t,v,v,color3);
		 objects.push_back(s);
	
	}


}



void Simulation::draw()
{
	for (int i = 0; i < objects.size(); i++)
	{
		objects[i]->draw();
	}

}
void Simulation::update(double simdt)
{
	for (int i = 0; i < objects.size(); i++)
	{
		for (int j = i+1; j < objects.size(); j++)
		{
			objects[i]->checkCollision(objects[j]);
		}

	}

	for (int i = 0; i < objects.size(); i++)
	{
		objects[i]->applyCollisionResponse();
		objects[i]->update(simdt);
		
	}

	double totalKinetic = 0;
	for (int i = 0; i < objects.size(); i++)
	{
		totalKinetic += objects[i]->getKinetik();
	}

	//std::cout << "TotalKinetic== " << totalKinetic << std::endl;
	//system("CLS");

	totalKinetic = 0;

}
