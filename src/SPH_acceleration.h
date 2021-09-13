#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"

void Compute_Acceleration
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
		realvec ac_temp = {0.,0.};
		real p_temp = 0;
		real visco = 0;

		real drs; //dr size
		real drdv; //vect(v) \cdot \nabla W
		realvec dW; //smoothing function gradient

		real h = particles.data_const.h;
		real m = particles.data_const.m;

		realvec ar; //position of ative particle
		realvec nr; //position of neighbour particle
		realvec dr; //position diference

		realvec av; //position of ative particle
		realvec nv; //position of neighbour particle
		realvec dv; //velocity diference

		real nrho; //neighbour density
		real arho; //actual particle densiy

		real np; //neighbour pressure
		real ap; //actual particle density

		double *kernel;

		//exprimet for(int i = 0; i < particles.cells[zz].np; i++)
		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}
			/*
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}

			if(particles.data.part_type[particles.cells[zz].cp[i]] == wall || particles.data.part_type[particles.cells[zz].cp[i]] == virt) {continue;}
			*/


			//Load data of actual particle
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];
			ap = particles.data.p[particles.cells[zz].cp[i]];

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				//experiment for(int n = 0; n < particles.cells[cl].np; n++)
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

					/*
				if(particles.data.part_type[particles.cells[cl].cp[n]] == outletf){continue;}
				*/

				//Load data of neighbour particle
				nr = particles.data.r[particles.cells[cl].cp[n]];
				nv = particles.data.v[particles.cells[cl].cp[n]];
				nrho = particles.data.rho[particles.cells[cl].cp[n]];
				np = particles.data.p[particles.cells[cl].cp[n]];

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
				drdv = dr.x*dv.x + dr.y*dv.y;

				/* get kernel values
				double *Wendland_kernel(double r, double h) */
				kernel = Wendland_kernel(drs, h);
				dW = dr * kernel[1]; // <--- check this, if its ok

				/* pressure term */
				p_temp = ap/pow(arho, 2) + np/pow(nrho,2);
				//p_temp = (ap + np)/(arho*nrho);

				/* real Artificial_Viscosity
				(real h, real drs, real drdv, real rho0, real c0, real alpha) */
				visco = Artificial_Viscosity(h, drs, drdv, particles.data_const.rho0, particles.data_const.cs, particles.data_const.avisc);

				// <--- check this, if its ok
				ac_temp.x = ac_temp.x - dW.x*(p_temp + visco)*m;
				ac_temp.y = ac_temp.y - dW.y*(p_temp + visco)*m;

				//ac_temp = ac_temp - dW*(p_temp + visco)*m; //+= operator is not overloaded yet

				/* Temp, step */
				p_temp = 0; visco = 0;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			//Assign acceleration to particle
			particles.data.a[particles.cells[zz].cp[i]] = ac_temp - particles.data_const.graviy;
			ac_temp = {0., 0.};

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


