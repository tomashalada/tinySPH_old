#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

#include "ALESPH_riemann_solver.h"


void Kernel_gradient_approx
(Particle_system &particles)
{

	//Prepare cells to find interacting particles
	unsigned int ncx, ncy;
	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	//Indices of neigbour cells

	#pragma omp parallel for
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

		/* Debug */
		//std::cout << "ACCELERATION -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

		//Temp. data, temp. variable
		real p_temp = 0;
		real visco = 0;
		real gamma = 0;

		realvec gradrhoF_sum = {0., 0.};
		realvec gradvxF_sum = {0., 0.};
		realvec gradvyF_sum = {0., 0.};

		realvec gradrhoB_sum = {0., 0.};
		realvec gradvxB_sum = {0., 0.};
		realvec gradvyB_sum = {0., 0.};


		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real dp = particles.data_const.dp;
		real c0 = particles.data_const.cs;
		real rho0 = particles.data_const.rho0;

		double *kernel;

		//exprimet for(int i = 0; i < particles.cells[zz].np; i++)
		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}

			realvec ar; //position of ative particle
			realvec av; //position of ative particle
			real arho; //actual particle densiy
			real aomega; //ALE

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i];
			ar = particles.data.r[ai];
			av = particles.data.v[ai];
			arho = particles.data.rho[ai];
			aomega = particles.special.omega[ai];


			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				//experiment for(int n = 0; n < particles.cells[cl].np; n++)
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

				if(particles.cells[zz].cp[i] == particles.cells[cl].cp[n]){continue;}

				realvec nr; //position of neighbour particle
				realvec dr; //position diference
				realvec nv; //position of neighbour particle
				realvec dv; //velocity diference
				real nrho; //neighbour density
				real nomega; //ALE

				//Load data of neighbour particle
				idx ni = particles.cells[cl].cp[n];
				nr = particles.data.r[ni];
				nv = particles.data.v[ni];
				nrho = particles.data.rho[ni];
				nomega = particles.special.omega[ni];

				real drs; //dr size
				realvec dW; //smoothing function gradient
				real W;

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));

				/* get kernel values
				double *Wendland_kernel(double r, double h) */
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];
				dW = dr * kernel[1];


				//ALE equations
				if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				{

					realvec nb = {0., 0.};
					nb = particles.special.n[particles.cells[cl].cp[n]]*(-1);

					if(nr.x>ar.x)
					{
					gradrhoF_sum.x += (nrho - arho)*nb.x*W*dp;
					gradrhoF_sum.y += (nrho - arho)*nb.y*W*dp;
					gradvxF_sum.x += (nv.x - av.x)*nb.x*W*dp;
					gradvxF_sum.y += (nv.x - av.x)*nb.y*W*dp;
					gradvyF_sum.x += (nv.y - av.y)*nb.x*W*dp;
					gradvyF_sum.y += (nv.y - av.y)*nb.y*W*dp;
					}

					if(nr.x<ar.x)
					{
					gradrhoB_sum.x += (nrho - arho)*nb.x*W*dp;
					gradrhoB_sum.y += (nrho - arho)*nb.y*W*dp;
					gradvxB_sum.x += (nv.x - av.x)*nb.x*W*dp;
					gradvxB_sum.y += (nv.x - av.x)*nb.y*W*dp;
					gradvyB_sum.x += (nv.y - av.y)*nb.x*W*dp;
					gradvyB_sum.y += (nv.y - av.y)*nb.y*W*dp;
					}

				}
				else
				{

					if(nr.x>ar.x)
					{
					gradrhoF_sum.x += (nrho - arho)*nomega*dW.x;
					gradrhoF_sum.y += (nrho - arho)*nomega*dW.y;
					gradvxF_sum.x += (nv.x - av.x)*nomega*dW.x;
					gradvxF_sum.y += (nv.x - av.x)*nomega*dW.y;
					gradvyF_sum.x += (nv.y - av.y)*nomega*dW.x;
					gradvyF_sum.y += (nv.y - av.y)*nomega*dW.y;
					}

					if(nr.x<ar.x)
					{
					gradrhoB_sum.x += (nrho - arho)*nomega*dW.x;
					gradrhoB_sum.y += (nrho - arho)*nomega*dW.y;
					gradvxB_sum.x += (nv.x - av.x)*nomega*dW.x;
					gradvxB_sum.y += (nv.x - av.x)*nomega*dW.y;
					gradvyB_sum.x += (nv.y - av.y)*nomega*dW.x;
					gradvyB_sum.y += (nv.y - av.y)*nomega*dW.y;
					}

				}

				gamma += W*m/nrho;
				//gamma += W;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			if(gamma == 0)
			{

				particles.special.gradrhoF[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvxF[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvyF[particles.cells[zz].cp[i]] = {0., 0.};

				particles.special.gradrhoB[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvxB[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvyB[particles.cells[zz].cp[i]] = {0., 0.};

			}
			else
			{

				//particles.special.gradrho[particles.cells[zz].cp[i]] = gradvx_temp*aomega;
				//particles.special.gradvx[particles.cells[zz].cp[i]] = gradvx_temp*aomega;
				//particles.special.gradvy[particles.cells[zz].cp[i]] = gradvx_temp*aomega;
				particles.special.gradrhoF[particles.cells[zz].cp[i]] = gradrhoF_sum*2.;
				particles.special.gradvxF[particles.cells[zz].cp[i]] = gradvxF_sum*2.;
				particles.special.gradvyF[particles.cells[zz].cp[i]] = gradvyF_sum*2.;

				particles.special.gradrhoB[particles.cells[zz].cp[i]] = gradrhoB_sum*2.;
				particles.special.gradvxB[particles.cells[zz].cp[i]] = gradvxB_sum*2.;
				particles.special.gradvyB[particles.cells[zz].cp[i]] = gradvyB_sum*2.;

			}

			gamma = 0;
			gradrhoF_sum = {0., 0.};
			gradvxF_sum = {0., 0.};
			gradvyF_sum = {0., 0.};

			gradrhoB_sum = {0., 0.};
			gradvxB_sum = {0., 0.};
			gradvyB_sum = {0., 0.};

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function
