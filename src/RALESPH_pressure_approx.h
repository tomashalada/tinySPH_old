#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"

#include "ALESPH_riemann_solver.h"
//Include state equation:
#include "SPH_pressure.h"


std::array<real, 4> RALESPH_Boundary_pressure_approximation
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	//Sum of variables to point ar
	realvec VEL_APPROX = {0. , 0.};
	real RHO_APPROX = 0;
	real P_APPROX = 0;
	std::array<real, 4> Variables_APPROX;

	//Temp variables
	real W_sum = 0;
	real rhoW_sum = 0;
	real pW_sum = 0;
	realvec vW_sum = {0., 0.};

	// Get cell and nb's cells for point of aprotixmation.
	const real kh = particles.data_const.kap*particles.data_const.h;
	const real of = particles.data_const.h*0.273;
	const real h = particles.data_const.h;
	const real m = particles.data_const.m;
	const real c0 = particles.data_const.cs;

	unsigned int npx, npy;
	unsigned int ncx, ncy;

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

	double *kernel;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		idx ni = particles.cells[cl].cp[n]; //actual particle index

		//Load data of neighbour particle
		const realvec nr = particles.data.r[ni]; //position of ative particle
		const realvec nv = particles.data.v[ni];; //position of ative particle
		const real nrho = particles.data.rho[ni];; //actual particle densiy
		const real np = particles.data.p[ni];
		const real no = particles.special.omega[ni];

		//Interation variables
		realvec dr; //position
		real drs; //dr size

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));

		//Assign kernel values
		kernel = Wendland_kernel(drs, h);

		//Boundary Riemann problem
		realvec vw = {0., 0.};
		vw = vw  - nv;
		realvec nb = {1., 0.};
		real pE = np - nrho*c0*(vw.x * nb.x + vw.y*nb.y);

		//Add to sum
		W_sum += kernel[0];
		//rhoW_sum += nrho*kernel[0];
		//pW_sum += + np*kernel[0];
		//vW_sum.x += nv.x*kernel[0];
		//vW_sum.y += nv.y*kernel[0];
		pW_sum += 2*pE*no*kernel[0];

		} // cycle over particles in neighbour cells
	} // cycle over neighbour cells

	//In case of empty neighbour return zero values, otherwise CSPH
	if(W_sum == 0)
	{
		RHO_APPROX = 0.;
		P_APPROX = 0.;
		VEL_APPROX = {0. ,0.};
	}
	else
	{
		RHO_APPROX += 0;
		P_APPROX +=  pW_sum;
		VEL_APPROX.x +=  vW_sum.x / W_sum;
		VEL_APPROX.y +=  vW_sum.y / W_sum;
	}

	//Return approximated values
	Variables_APPROX[0] = RHO_APPROX;
	Variables_APPROX[1] = P_APPROX;
	Variables_APPROX[2] = VEL_APPROX.x;
	Variables_APPROX[3] = VEL_APPROX.y;

	return Variables_APPROX;
}


