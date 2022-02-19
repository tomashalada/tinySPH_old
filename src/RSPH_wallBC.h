#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

#include "SPH_kernel_approx.h"
#include "SPH_kernel.h"

void wallRSPH_precompute_normals
(Particle_system &particles, Simulation_data simulation_data)
{
	#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{
		if(particles.data.part_type[p] == wall)
		{

			realvec PHI_sum = {0., 0.};
			const realvec ar = particles.data.r[p];

			/* Get cell and nb's cells for point of aprotixmation. */
			const real kh = particles.data_const.kap*particles.data_const.h;
			const real of = particles.data_const.h*0.273;
			const real h = particles.data_const.h;
			const real m = particles.data_const.m;

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

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				//experiment for(int n = 0; n < particles.cells[cl].np; n++)
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

					idx ni = particles.cells[cl].cp[n]; //neighbour particle index
					if(particles.data.part_type[ni] != wall){continue;}

					const realvec nr = particles.data.r[ni]; //position of ative particle
					const real nrho = particles.data.rho[ni];; //actual particle densiy
					const real np = particles.data.p[ni];

					//Interation variables
					realvec dr; //position, velocity difference
					real drs; //dr size
					realvec dW; //kernel gradient
					real drdv; //dot product of velocity dif. and position dif.
					real dvdW; //dot product of velocity dif. and smoothing function gradient
					real drdW; //dot product of position dif. and smoothing function gradient
					double *kernel;

					//Position and velocity difference
					dr = ar - nr;
					drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));

					//Assign kernel values
					kernel = Wendland_kernel(drs, h);
					dW = dr * kernel[1];
					drdW = dr.x*dW.x + dr.y*dW.y;

					PHI_sum.x += - m/nrho * dW.x;
					PHI_sum.y += - m/nrho * dW.y;

				} // if function - check part type
			} // cycle over nb cells

			particles.special.n[p] = PHI_sum/NORM(PHI_sum.x, PHI_sum.y); // if function - check part type
			PHI_sum = {0., 0.};

		} // if function - check part type
	} // cycle over particles
} // function

