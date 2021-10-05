/* TODO */
/*
			- Univerzalni funkci pro generovani linek castic
			- Varovani, pokud delka steny nekoresponduje se zadanymi rozmery a dp
*/


#include "SPH_particle_system.h"



 void Create_Block_of_Fluid(Particle_system  &particles, double x_0, double y_0, double x_1, double y_1)
 {

 	real dp = particles.data_const.dp;
 	realvec r, v, a;

 	//gets length of wall line and number of particles for the wall
 	real len_x = LEN(x_0, y_0, x_1, y_0);
 	real len_y = LEN(x_0, y_0, x_0, y_1);


 	int npart_x = (int)(len_x/dp);
 	int npart_y = (int)(len_y/dp);

 	/* Debug */
 	std::cout << "CREATE -> len_x: " << len_x << " len_y: " << len_y << std::endl;
 	std::cout << "CREATE -> npart_x: " << npart_x << " npart_y: " << npart_y << std::endl;

 	for(int i = 0; i < npart_x+1; i ++){

 		r.x = x_0 + i*dp;

 		for(int j = 0; j < npart_y+1; j++){

 			r.y = y_0 + j*dp;

 			v.x = 0; v.y = 0;
 			a.x = 0; a.y = 0;

 			particles.Add_particle(fluid, 0, particles.data_const.rho0, r, v, a);

 		}
 	}

 }

void Create_Wall_x_bt(Particle_system & particles, double x_0, double x_1, double y, int orientation, realvec n)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;
	realvec gn;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = 1; //+1 tempfor asdasd
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.y = y + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.x = x_0 + i*dp;
			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;
			particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

			particles.special.n.push_back(n);
		}
	}

}

void Create_Wall_y_bt(Particle_system & particles, double x, double y_0, double y_1, int orientation, realvec n)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;
	realvec gn;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = 1;
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.x = x + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.y = y_0 + i*dp;
			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;

			particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

			particles.special.n.push_back(n);
		}
	}

}

void Create_corner_bt(Particle_system & particles, int orientationX, int orientationY, realvec xyv, realvec n)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;
	realvec gn;

	// number of virtual layers
	int nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.y = xyv.y + dp*l*orientationY;

		for(int i = 0; i < nvl; i++)
		{
			r.x = xyv.x + dp*i*orientationX;

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;
			particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

			particles.special.n.push_back(n);

		}
	}

}
