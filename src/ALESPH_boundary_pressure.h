#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

#include "ALESPH_riemann_solver.h"


/* Approximate velocity with kernel aproximation to certain point */
real Boundary_Pressure_ALESPH
(Particle_system &particles, Simulation_data simulation_data, realvec ar, realvec bn)
{

	//realvec VEL_APROX = {0. , 0.};
	real wp = 0; //wall pressure

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;

	real rho0 = particles.data_const.rho0;
	real c0 = particles.data_const.cs;

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

	realvec vWsum = {0., 0.};
	real Wsum = 0;
	realvec vvv = {0., 0.};
	real mWsum = 0;
	//realvec vWsum;
	//real Wsum ;

	real pw_temp = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		realvec nr = {0., 0.};
		realvec nv = {0., 0.};
		real nrho = 0;
		real np = 0;
		real nomega = 0;

		realvec dr = {0., 0.};
		real drs = 0;

		real ps;

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];
		np = particles.data.p[particles.cells[cl].cp[n]];
		nomega = particles.special.omega[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		vWsum = vWsum + nv*kernel[0]; // += not overloaded yet
		vvv = vvv + nv*kernel[0]*(m/nrho);
		mWsum += m*kernel[0]/nrho;

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;
		ps = np - nrho * cs(nrho, rho0, c0) * (bn.x * (0 - nv.x) + bn.y * (0 - nv.y));
		pw_temp += nomega * kernel[0] * 2 * ps;


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	wp = pw_temp;
	if(wp > 0){std::cout << "[BOUNDARY PRESSURE] Wall pessure: " << wp << std::endl;}


	return wp;
}

void Wall_pressure
(Particle_system &particles, Simulation_data simulation_data)
{

	realvec r;
	realvec bn;

	#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] == wall)
		{

			r = particles.data.r[p];
			bn = particles.special.n[p];

			particles.data.p[p] = Boundary_Pressure_ALESPH(particles, simulation_data, r, bn);

		} // if function - check part type

	} // cycle over particles

}


