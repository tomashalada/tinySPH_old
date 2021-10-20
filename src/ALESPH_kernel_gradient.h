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
		realvec gradrho_temp = {0., 0.};
		realvec gradvx_temp = {0., 0.};
		realvec gradvy_temp = {0., 0.};


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
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];
			aomega = particles.special.omega[particles.cells[zz].cp[i]];


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
				nr = particles.data.r[particles.cells[cl].cp[n]];
				nv = particles.data.v[particles.cells[cl].cp[n]];
				nrho = particles.data.rho[particles.cells[cl].cp[n]];
				nomega = particles.special.omega[particles.cells[cl].cp[n]];

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

					gradrho_temp.x += (nrho - arho)*nb.x*W*dp;
					gradrho_temp.y += (nrho - arho)*nb.y*W*dp;

					gradvx_temp.x += (nv.x - av.x)*nb.x*W*dp;
					gradvx_temp.y += (nv.x - av.x)*nb.y*W*dp;

					gradvy_temp.x += (nv.y - av.y)*nb.x*W*dp;
					gradvy_temp.y += (nv.y - av.y)*nb.y*W*dp;

				}
				else
				{

					gradrho_temp.x += (nrho - arho)*nomega*dW.x;
					gradrho_temp.y += (nrho - arho)*nomega*dW.y;

					gradvx_temp.x += (nv.x - av.x)*nomega*dW.x;
					gradvx_temp.y += (nv.x - av.x)*nomega*dW.y;

					gradvy_temp.x += (nv.y - av.y)*nomega*dW.x;
					gradvy_temp.y += (nv.y - av.y)*nomega*dW.y;

				}

				gamma += W*m/nrho;
				//gamma += W;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			if(gamma == 0)
			{

				particles.special.gradrho[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvx[particles.cells[zz].cp[i]] = {0., 0.};
				particles.special.gradvy[particles.cells[zz].cp[i]] = {0., 0.};

			}
			else
			{

				particles.special.gradrho[particles.cells[zz].cp[i]] = gradvx_temp;
				particles.special.gradvx[particles.cells[zz].cp[i]] = gradvx_temp;
				particles.special.gradvy[particles.cells[zz].cp[i]] = gradvx_temp;

			}

			gamma = 0;
			gradrho_temp = {0., 0.};
			gradvx_temp = {0., 0.};
			gradvy_temp = {0., 0.};

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function
