#pragma once
#include "SPH_defs.h"

// Active valuces
struct SPH_data
{

	std::vector<idx> index; //id of particle
	std::vector<unsigned short int> part_type; //type of particle (fluid, wall, virtual, buffer...)

	std::vector<real> p;

	std::vector<real> rho;
	std::vector<real> drho;

	std::vector<realvec> r;
	std::vector<realvec> v;
	std::vector<realvec> a;

	//TEST
	std::vector<realvec> r_o;
	std::vector<realvec> v_o;
	std::vector<realvec> v_oo;
	std::vector<real> rho_o;
	std::vector<real> rho_oo;

};

//Mostly const. parameters
struct SPH_constants
{

	real h; //smoothig length
	real kap; //kernel support coef.
	real dp; //initial particle distance

	real cs; //numerical speed of sound
	real rho0; //reference density

	realvec graviy; //grav acc.
	real avisc; //alpha parameter of artificial viscosity
	real delta; //density diffusion term constant

	real m; //particle mass

};

struct SPH_special
{

	std::vector<realvec> gnr; //mDBC, ghost node position for walls
	std::vector<realvec> n; //BT, boundary particles normals

	std::vector<real> omega; //ALESPH particle volume
	std::vector<real> omegarho;
	std::vector<realvec> omegav;

	std::vector<real> omega_o; //ALESPH particle volume
	std::vector<real> omegarho_o;
	std::vector<realvec> omegav_o;
	std::vector<real> omega_oo; //ALESPH particle volume
	std::vector<real> omegarho_oo;
	std::vector<realvec> omegav_oo;

	std::vector<real> domega; //ALESPH particle volume change
	std::vector<real> domegarho; //ALESPH density and omega change
	std::vector<realvec> domegav; //ALESPH acceleration and omega

	std::vector<realvec> gradrho; //ALESPH MUSCL
	std::vector<realvec> gradvx; //ALESPH MUSCL
	std::vector<realvec> gradvy; //ALESPH MUSCL

	//EXPERIMENT
	std::vector<realvec> gradrhoF; //ALESPH MUSCL
	std::vector<realvec> gradvxF; //ALESPH MUSCL
	std::vector<realvec> gradvyF; //ALESPH MUSCL
	std::vector<realvec> gradrhoB; //ALESPH MUSCL
	std::vector<realvec> gradvxB; //ALESPH MUSCL
	std::vector<realvec> gradvyB; //ALESPH MUSCL

	real dtcv_temp=0;

};

//Pairs etc. informations
struct Interaction_data
{

	unsigned int ncx, ncy; //number of cells

};

//Cell structure for pairs searching
struct Cell
{

	idx id;
	//realvec r;

	unsigned int np; //number of particles
	std::vector<idx> cp; //contained particles

	Cell
		//(idx _id, realvec _r) : id(_id), r(_r) {}
		(idx _id) : id(_id) { np = 0; }

};

//Domain data
struct Simulation_data
{

	real x_0, x_m, y_0, y_m; //domain size
	real hh; //smoothed length obtained from dp
	real kh; //effective radius

	/* Figure out better. */
	int nvl; //number of layer of virtual particles in walls

	Simulation_data
		(real _x_0, real _x_m, real _y_0, real _y_m ) : x_0(_x_0), x_m(_x_m), y_0(_y_0), y_m(_y_m) {}


};
