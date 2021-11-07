/* TODO */
/*
			- Univerzalni funkci pro generovani linek castic
			- Varovani, pokud delka steny nekoresponduje se zadanymi rozmery a dp
*/


#include "SPH_particle_system.h"
#include "SPH_pressure.h"


void Create_Wall(Particle_system & particles, double x_0, double y_0, double x_1, double y_1)
{
	real dp = particles.data_const.dp;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y_0, x_1, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	for(int i = 0; i < npart+1; i++)
	{
		r.x = 0;
		r.y = 0;

		v.x = 0; v.y = 0;
		a.x = 0; a.y = 0;

		particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

	}


}

void Create_Wall_x(Particle_system & particles, double x_0, double x_1, double y, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp) +1; //+1 tempfor asdasd
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.y = y + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.x = x_0 + i*dp;

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;

			particles.ALESPH_Add_particle(wall, 0, particles.data_const.rho0, r, v, a);
		}
	}

}

void Create_Wall_y(Particle_system & particles, double x, double y_0, double y_1, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp) ;
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
		}
	}

}

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

			particles.ALESPH_Add_particle(fluid, 0, particles.data_const.rho0, r, v, a);

		}
	}

}

void Create_Block_of_Fluid_with_speed(Particle_system  &particles, double x_0, double y_0, double x_1, double y_1, realvec vel)
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

			v.x = vel.x; v.y = vel.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(fluid, 0, particles.data_const.rho0, r, v, a);

		}
	}

}

/*** *** *** *** *** *** *** ***/
void Create_Wall_x_mov(Particle_system & particles, double x, double y_0, double y_1, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp);
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.x = x + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.y = y_0 + i*dp;

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;

			particles.Add_particle(mov_a, 0, particles.data_const.rho0, r, v, a);
		}
	}

}

void Create_Wall_y_mov(Particle_system & particles, double x, double y_0, double y_1, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp);
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.x = x + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.y = y_0 + i*dp;

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;

			particles.Add_particle(mov_a, 0, particles.data_const.rho0, r, v, a);
		}
	}

}

/*** *** *** *** *** *** *** ***/

void Inlet_x_direction(Particle_system &particles, double x_0, double x_1, double y, int orientation, int nbl, realvec v0)
{

	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = nbl; //number of virtual layers = number of buffer layers
	//nvl = 1;

	//for(int l = 0; l < nvl; l++)
	for(int l = 0; l < nbl; l++)
	{

		r.y = y + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.x = x_0 + i*dp;

			v.x = v0.x; v.y = v0.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(inlet, 0, particles.data_const.rho0, r, v, a);
		}
	}


}

void Inlet_y_direction(Particle_system &particles, double x, double y_0, double y_1 , int orientation, int nbl, realvec v0)
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp);
	//nvl = 1;


	//for(int l = 0; l < nvl; l++)
	for(int l = 0; l < nbl; l++)
	{

		r.x = x + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.y = y_0 + i*dp;

			v.x = v0.x; v.y = v0.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(inlet, 0, particles.data_const.rho0, r, v, a);
		}
	}

}

void Outlet_y_direction(Particle_system &particles, double x, double y_0, double y_1 , int orientation, int nbl, realvec v0)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp);
	//nvl = 1;


	//for(int l = 0; l < nvl; l++)
	for(int l = 0; l < nbl; l++)
	{

		r.x = x + dp*l*orientation;

		for(int i = 0; i < npart+1; i++)
		{
			r.y = y_0 + i*dp;

			v.x = v0.x; v.y = v0.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(outlet, 0, particles.data_const.rho0, r, v, a);
		}
	}

}

/*** *** *** *** *** *** *** ***/

void Create_Wall_x_mDBC(Particle_system & particles, double x_0, double x_1, double y, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp);
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
		}
	}

}

// INLET EXPERIMENT
struct Inlet_buffer_x
{

	real xb_0, yb_0, xb_m, yb_m; //size of buffer (filled with particles)
	realvec v_inl; //inlet velocity
	idx npb; //number of particles in buffer

	real b_len; //buffer length
	real b_hig; //buffer hight

	//(real _x_0, real _x_m, real _y_0, real _y_m ) : x_0(_x_0), x_m(_x_m), y_0(_y_0), y_m(_y_m) {}
	Inlet_buffer_x
	(Particle_system &particles, real _xb_0, real _yb_0, real _xb_m, real _yb_m, real dp, realvec _v_inl)
	: xb_0(_xb_0), yb_0(_yb_0), xb_m(_xb_m), yb_m(_yb_m), v_inl(_v_inl)
	{

		b_len = xb_m - xb_0;
		b_hig = yb_m - yb_0;

		int npart_x = (int)(b_len/dp);
		int npart_y = (int)(b_hig/dp);

		realvec r, v, a;

		for(int i = 0; i < npart_x+1; i ++)
		{

			r.x = xb_0 + i*dp;

			for(int j = 0; j < npart_y+1; j++)
			{

				r.y = yb_0 + j*dp;

				v.x = v_inl.x; v.y = v_inl.y;
				a.x = 0; a.y = 0;

				particles.Add_particle(fluid, 0, particles.data_const.rho0, r, v, a);

			}

		}

	}

};

/*** NEW AND BETTER CREATE FUNCTIONS ***/
void DRAW_CreateInletY(Particle_system &particles, double x, double y_0, double y_1 , int orientation, int nbl, realvec v0, int type)
{

	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x, y_0, x, y_1);
	unsigned int npart = std::ceil(len/dp) + 1;

	//create particles
	for(int l = 0; l < nbl; l++)
	{
		r.x = x + dp*l*orientation;

		for(int h = 0; h < npart; h++)
		{
			r.y = y_0 + h*dp;

			v.x = v0.x; v.y = v0.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(type, 0, particles.data_const.rho0, r, v, a);
		}
	}

	//write buffer information
	//particles.zones.push_back(InletBuffer(x, y_0, x + nbl*dp*orientation, y_1, nbl, v0));
	particles.zones.push_back(BufferZone(x, x + nbl*dp*orientation, y_0, y_1, nbl, v0, {(real)orientation, 0}, type));

} // function

void DRAW_CreateInletX(Particle_system &particles, double y, double x_0, double x_1 , int orientation, int nbl, realvec v0, int type)
{

	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	unsigned int npart = std::ceil(len/dp) + 1;

	//create particles
	for(int l = 0; l < nbl; l++)
	{
		r.y = y + dp*l*orientation;

		for(int h = 0; h < npart; h++)
		{
			r.x = x_0 + h*dp;

			v.x = v0.x; v.y = v0.y;
			a.x = 0; a.y = 0;

			particles.Add_particle(type, 0, particles.data_const.rho0, r, v, a);
		}
	}

	//write buffer information
	//particles.zones.push_back(InletBuffer(y, y_0, y + nbl*dp*orientation, y_1, nbl, v0));
	particles.zones.push_back(BufferZone(x_0, x_1, y, y + nbl*dp*orientation, nbl, v0, {0, (real)orientation}, type));

} // function

/*** New better functions ***/

/*** NEW AND BETTER CREATE FUNCTIONS ***/
void DRAW_CreateOutletY(Particle_system &particles, double x, double y_0, double y_1 , int orientation, int nbl, realvec v0, int type, bool empty)
{

		real dp = particles.data_const.dp;
		real kh = particles.data_const.h * particles.data_const.kap;
		realvec r, v, a;

	if(empty == false){

		//gets length of wall line and number of particles for the wall
		real len = LEN(x, y_0, x, y_1);
		unsigned int npart = std::ceil(len/dp) + 1;

		//create particles
		for(int l = 0; l < nbl; l++)
		{
			r.x = x + dp*l*orientation;

			for(int h = 0; h < npart; h++)
			{
				r.y = y_0 + h*dp;

				v.x = v0.x; v.y = v0.y;
				a.x = 0; a.y = 0;

				particles.Add_particle(inlet, 0, particles.data_const.rho0, r, v, a);
			}
		}
	}

	//write buffer information
	//particles.zones.push_back(InletBuffer(x, y_0, x + nbl*dp*orientation, y_1, nbl, v0));
	particles.zones.push_back(BufferZone(x, x + nbl*dp*orientation, y_0, y_1, nbl, v0, {(real)orientation, 0}, type));

} // function

///===================================================================================================///

void Create_Wall_x_mdbc(Particle_system & particles, double x_0, double x_1, double y, int orientation, double y_w)
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
	int nvl = std::ceil(kh/dp) +1; //+1 tempfor asdasd
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

			gn.x = r.x;
			gn.y = y_w + fabs(r.y - y_w);
			particles.special.gnr.push_back(gn);
		}
	}

}

void Create_Wall_y_mdbc(Particle_system & particles, double x, double y_0, double y_1, int orientation, double x_w)
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
	int nvl = std::ceil(kh/dp) +1;
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

			gn.x = x_w + fabs(r.x - x_w);
			gn.y = r.y;
			particles.special.gnr.push_back(gn);
		}
	}

}

void Create_corner_mdbc(Particle_system & particles, int orientationX, int orientationY, realvec xyv)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;
	realvec gn;

	// number of virtual layers
	int nvl = std::ceil(kh/dp) +1; //+1 tempfor asdasd

	for(int l = 0; l < nvl; l++)
	{

		r.y = xyv.y + dp*l*orientationY + dp/2*orientationY;

		for(int i = 0; i < nvl; i++)
		{
			r.x = xyv.x + dp*i*orientationX + dp/2*orientationX;

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;
			particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

			gn.x = xyv.x + fabs(r.x - xyv.x);
			gn.y = xyv.y + fabs(r.y - xyv.y);
			particles.special.gnr.push_back(gn);

		}
	}

}

void Create_Wall_x_lattice(Particle_system & particles, double x_0, double x_1, double y, int orientation)
// orientation +1: -1:
{
	real dp = particles.data_const.dp;
	real kh = particles.data_const.h * particles.data_const.kap;
	realvec r, v, a;

	//gets length of wall line and number of particles for the wall
	real len = LEN(x_0, y, x_1, y);
	int npart = (int)(len/dp);
	/*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

	// number of virtual layers
	int nvl = std::ceil(kh/dp) +1; //+1 tempfor asdasd
	//nvl = 1;

	for(int l = 0; l < nvl; l++)
	{

		r.y = y + dp*l*orientation/2;

		for(int i = 0; i < npart+1; i++)
		{
			if(l%2 == 0){r.x = x_0 + i*dp + dp/2;}
			else{r.x = x_0 + i*dp;}

			v.x = 0; v.y = 0;
			a.x = 0; a.y = 0;

			particles.Add_particle(wall, 0, particles.data_const.rho0, r, v, a);
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
			particles.ALESPH_Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

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

			particles.ALESPH_Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

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
			particles.ALESPH_Add_particle(wall, 0, particles.data_const.rho0, r, v, a);

			particles.special.n.push_back(n);

		}
	}

}

void Create_Block_of_FluidWP(Particle_system  &particles, double x_0, double y_0, double x_1, double y_1)
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

			real rhoin = Pressure_to_density(9.81*(y_1 - r.y)*1000, 1000., particles.data_const.cs);

			//particles.ALESPH_Add_particle(fluid, 0, particles.data_const.rho0, r, v, a);
			particles.ALESPH_Add_particle(fluid, 9.81*(y_1 - r.y)*1000, rhoin, r, v, a);

		}
	}

}
