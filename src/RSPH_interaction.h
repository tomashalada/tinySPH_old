#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"

#include "ALESPH_riemann_solver.h"

void Compute_Forces
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

		real drho_sum = 0; //density sum
		real gamma_sum = 0; //normalization factor

		const real h = particles.data_const.h;
		const real m = particles.data_const.m;
		const real c0 = particles.data_const.cs;
		const real delta = particles.data_const.delta;
		const real avisco = particles.data_const.avisc;
		const real rho0 = particles.data_const.rho0;
		const realvec gravity = particles.data_const.graviy;
		const real dp = particles.data_const.dp;

		double *kernel;

		real dt_visco = 0;
		real dtacc = 0;

		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i]; //actual particle index

			const realvec ar = particles.data.r[ai];; //position of ative particle
			const realvec av = particles.data.v[ai]; //position of ative particle
			const real arho = particles.data.rho[ai]; //actual particle densiy
			const real ap = particles.data.p[ai]; //actual particle density


			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

				//Load data of neighbour particle
				idx ni = particles.cells[cl].cp[n]; //actual particle index
				if(ai == ni){continue;}

				const realvec nr = particles.data.r[ni]; //position of ative particle
				const realvec nv = particles.data.v[ni];; //position of ative particle
				const real nrho = particles.data.rho[ni];; //actual particle densiy
				const real np = particles.data.p[ni];

				//Interation variables
				realvec dr, dv; //position, velocity difference
				real drs; //dr size
				real W; //kernel
				realvec dW; //kernel gradient
				real drdv; //dot product of velocity dif. and position dif.
				real dvdW; //dot product of velocity dif. and smoothing function gradient
				real drdW; //dot product of position dif. and smoothing function gradient

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));
				drdv = dr.x*dv.x + dr.y*dv.y;

				//Assign kernel values
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];
				dW = dr * kernel[1];
				dvdW = dv.x*dW.x + dv.y*dW.y;
				drdW = dr.x*dW.x + dr.y*dW.y;

				// ========== SOLVE RIEMANN PROBLEM ========== //
				realvec drn = dr*(1/drs)*(-1); //pair normal vector
				real vL = av.x * drn.x + av.y * drn.y; //v-left state
				real vR = nv.x * drn.x + nv.y * drn.y; //v-right state

				real vss = velRiemannLinearizedwithPressure(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				//real ps = pRiemannLinearized(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				real ps = pRiemannLinearizedWithLimiter(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5); //reconstruct

				// ========== MOMENTUM EQN ========== //
				real pR = 0; //wall pressure

				if(particles.data.part_type[ni] == wall)
				{
					realvec nb = particles.special.n[ni];
				 real	drn = dr.x*nb.x + dr.y*nb.y;

					realvec rwf = dr;

					vL = (-1)*( av.x * nb.x + av.y * nb.y );
					vR = -vL;
					pR = ap + arho*(-1)*(particles.data_const.graviy.y*rwf.y);

					ps = pRiemannLinearized(arho, nrho, vL, vR, ap, pR, 0.5*(arho + nrho), c0);
				 //ps = pRiemannLinearizedWithLimiter(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);

					ac_sum.x -= 2.*dW.x*(ps/(arho*nrho))*m;
					ac_sum.y -= 2.*dW.y*(ps/(arho*nrho))*m;
				}
				else
				{
					ac_sum.x -= 2.*dW.x*(ps/(arho*nrho))*m;
					ac_sum.y -= 2.*dW.y*(ps/(arho*nrho))*m;
				}

				// ========== CONTINUITY EQN ========== //
				realvec dvs = av - vs;

				if(particles.data.part_type[ni] == wall)
				{
					realvec nb = particles.special.n[ni];
					real dvn = dvs.x*nb.x + dvs.y*nb.y;

					drho_sum -=  2.*nrho*dvn*W*dp;
				}
				else
				{
					real dvsdW = dvs.x*dW.x + dvs.y*dW.y;

					drho_sum += 2.*dvsdW*m;
				}

				// ==================================== //

				gamma_sum += W*m/nrho;

				} // cycle over particles in neighbour cells
			} // cycle over neighbour cells

			//Assign acceleration to particle
			particles.data.a[ai] = ac_sum - particles.data_const.graviy;
			//particles.data.rho[ai] = drho_sum*arho + diff_sum;
			particles.data.drho[ai] = drho_sum;


			//Reset temp. variables
			ac_sum = {0., 0.};
			drho_sum = 0;
			gamma_sum = 0;

		} // cycle over particles in active cell
	} // cycle over cells

	#pragma omp barrier

} // function


