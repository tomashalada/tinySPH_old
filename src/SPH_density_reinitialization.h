#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"

void Density_reinitialization_sumation
(Particle_system &particles)
{

	//Prepare cells to find interacting particles
	unsigned int ncx, ncy;
	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	//#pragma omp parallel for
	#pragma omp parallel for schedule(dynamic, 3)
	for(int c = 0; c < particles.cells.size(); c++)
	{

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

		//Temp. data, temp. variable
		realvec ac_sum = {0.,0.}; //acceleration sum

		real rho_sum = 0; //density sum
		real WV_sum = 0;

		const real h = particles.data_const.h;
		const real m = particles.data_const.m;
		const real c0 = particles.data_const.cs;
		const real delta = particles.data_const.delta;
		const real avisco = particles.data_const.avisc;
		const real rho0 = particles.data_const.rho0;
		const real dp = particles.data_const.dp;
		const realvec gravity = particles.data_const.graviy;

		double *kernel;

		real dt_visco = 0;
		real dtacc = 0;

		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i]; //actual particle index

			const realvec ar = particles.data.r[ai];; //position of ative particle
			const real arho = particles.data.rho[ai]; //actual particle densiy


			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

				//Load data of neighbour particle
				idx ni = particles.cells[cl].cp[n]; //actual particle index
				if(ai == ni){continue;}

				const realvec nr = particles.data.r[ni]; //position of ative particle
				const real nrho = particles.data.rho[ni];; //actual particle densiy

				//Interation variables
				realvec dr; //position, velocity difference
				real drs; //dr size
				real W; //kernel value

				//Position and velocity difference
				dr = ar - nr;
				drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));

				//Assign kernel values
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];

				// ========== Density_sumation ========== //
				rho_sum += m*W;
				WV_sum += m*W/nrho;

				//rho_sum += nrho*W;
				//WV_sum += W;
				// ====================================== //

				} // cycle over particles in neighbour cells
			} // cycle over neighbour cells

			if(WV_sum != 0)
			particles.special.rhoTemp[ai] = rho_sum/WV_sum;

			//Reset temp. variables
			rho_sum = 0;
			WV_sum = 0;

		} // cycle over particles in active cell
	} // cycle over cells

	#pragma omp barrier

} // function

void Density_reinitialization_overwrite
(Particle_system &particles)
{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] != fluid){continue;}

			if(particles.data.rho[p] < 990 || particles.data.rho[p] > 1010)
			std::cerr << "Particle density after reinitialisation too high: rho = " << particles.data.rho[p] << " rhoTemp = " << particles.special.rhoTemp[p] << std::endl;

			particles.data.rho[p] = particles.special.rhoTemp[p];

		}
} // function

