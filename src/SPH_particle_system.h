#pragma once

#include "SPH_defs.h"
#include "SPH_data.h"
#include "SPH_inlet_outlet_data.h"

struct Particle_system
{

	unsigned int np; //number of particles

	SPH_data data; //main SPH data (pos, vel, p, rho...)
	SPH_constants data_const; //constants of system
	SPH_special special; //special variables (ghost node pos., ...)
	Interaction_data pairs; //information about interaction in actual time step

	/* Temp, output testing */
	SPH_data particles_out; //list of particles leaving the simulation
	unsigned int np_out; // number of particles leaving the simulation

	std::vector<Cell> cells; //grid over domain to work with particle pairs
	std::vector<BufferZone> zones; //information about inlet, outlets

	//Function to create new particles
	void Add_particle
	(short unsigned int _type, real _p, real _rho, realvec _r, realvec _v, realvec _a)
	{

		np++;
		data.part_type.push_back(_type);
		data.index.push_back(np - 1); //particles are labeled from 0

		data.p.push_back(_p);
		data.rho.push_back(_rho);
		data.drho.push_back(0.); //init drho/dt with 0

		data.r.push_back(_r);
		data.v.push_back(_v);
		data.a.push_back(_a);

		//TEST
		data.r_o.push_back(_r);
		data.rho_o.push_back(_rho);
		data.rho_oo.push_back(_rho);
		data.v_o.push_back(_v);
		data.v_oo.push_back(_v);

	}

	//Function to create new particles
	void ALESPH_Add_particle
	(short unsigned int _type, real _p, real _rho, realvec _r, realvec _v, realvec _a)
	{

		np++;
		data.part_type.push_back(_type);
		data.index.push_back(np - 1); //particles are labeled from 0

		data.p.push_back(_p);
		data.rho.push_back(_rho);
		data.drho.push_back(0.); //init drho/dt with 0

		data.r.push_back(_r);
		data.v.push_back(_v);
		data.a.push_back(_a);


		//ALESPH variables
		special.omega.push_back(data_const.dp * data_const.dp);
		special.domega.push_back(0.);
		special.domegarho.push_back(0.);
		special.omegaa.push_back({0., 0.});
		special.omegav.push_back(_v*(data_const.dp * data_const.dp));
		special.omegarho.push_back(_rho*(data_const.dp * data_const.dp));

		//ALESPH MUSCL variables
		special.gradrho.push_back({0., 0.});
		special.gradvx.push_back({0., 0.});
		special.gradvy.push_back({0., 0.});

	}

	//Function to create new particles
	void Remove_particle
	(idx p) /* There also should be type, to confirm p match right part. */
	{

		//pList.erase(pList.begin()+i);
		np_out++;
		particles_out.part_type.push_back(data.part_type[p]);
		particles_out.index.push_back(data.index[p]);

		particles_out.r.push_back(data.r[p]);
		particles_out.v.push_back(data.v[p]);

		particles_out.p.push_back(data.p[p]);
		particles_out.rho.push_back(data.rho[p]);

		np--;
		data.part_type.erase(data.part_type.begin() + p);
		data.index.erase(data.index.begin() + p); //particles are labeled from 0

		data.p.erase(data.p.begin() + p);
		data.rho.erase(data.rho.begin() + p);
		data.drho.erase(data.drho.begin() + p); //init drho/dt with 0

		data.r.erase(data.r.begin() + p);
		data.v.erase(data.v.begin() + p);
		data.a.erase(data.a.begin() + p);


	}

	//Computation on particles
	void Compute_Acceleration();


	//Constructor
	Particle_system
	/*(int _np, double _h, double _kap, double _cs, double _rho0, double _dp) : np(_np)*/
	/*(int _np, double _h, double _kap, double _cs, double _rho0, double _dp) */
	(double _h, double _kap, double _cs, double _rho0, double _visco, double _delta, double _dp)
	{

	np = 0; //system is initialized without particles
	np_out = 0; //system is initialized without particles

 data_const.h=_h;
 data_const.kap=_kap;
 data_const.cs=_cs;
 data_const.rho0=_rho0;
	data_const.dp=_dp;
	data_const.graviy = {0., 9.81};
	//data_const.graviy = {0., 0.};
	data_const.m = _rho0 / (1/_dp * 1/_dp);
	data_const.avisc = _visco;
	data_const.delta = _delta;

	}

};



