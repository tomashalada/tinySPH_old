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

	real dtcv = 0;

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

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real c0 = particles.data_const.cs;
		real rho0 = particles.data_const.rho0;

		double *kernel;

		real dt_visco = 0;
		real dtacc = 0;

		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}
			/*
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}
			*/

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i]; //actual particle index

			realvec ar = particles.data.r[ai];; //position of ative particle
			realvec av = particles.data.v[ai]; //position of ative particle
			real arho = particles.data.rho[ai]; //actual particle densiy
			real ap = particles.data.p[ai]; //actual particle density

			realvec nr; //position of neighbour particle
			realvec dr; //position diference


			realvec dv; //velocity diference

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
				idx ni = particles.cells[cl].cp[n]; //actual particle index
				//if(ai == ni){continue;}

				realvec nr = particles.data.r[ni]; //position of ative particle
				realvec nv = particles.data.v[ni];; //position of ative particle
				real nrho = particles.data.rho[ni];; //actual particle densiy
				real np = particles.data.p[ni];

				//Interation variables
				realvec dr, dv; //position, velocity difference
				real drs; //dr size
				realvec dW; //kernel gradient
				real drdv; //dot product of velocity dif. and position dif.

				real p_temp = 0; //pressure term
				real visco = 0; //viscosity term

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));
				drdv = dr.x*dv.x + dr.y*dv.y;

				//Assign kernel values
				//double *Wendland_kernel(double r, double h)
				kernel = Wendland_kernel(drs, h);
				dW = dr * kernel[1];

				//Assign acceleration terms
				//p_temp = ap/pow(arho, 2) + np/pow(nrho,2);
				p_temp = (ap + np)/(arho*nrho);
				real alpha = particles.data_const.avisc;
				//if(particles.data.p[ni] == wall){alpha = alpha*(-1)*100;}
				visco = Artificial_Viscosity(h, drs, drdv, 0.5*(arho + nrho), c0, alpha);

				real dt_visco_i = 0;
				dt_visco_i = fabs(h*drdv/(drs*drs));
				//std::cout << "dt_visco_i: " << dt_visco_i << std::endl;
				if(dt_visco_i > dt_visco){ dt_visco =  dt_visco_i; }

				//Assign sum to particle
				ac_temp.x -= dW.x*(p_temp + visco)*m;
				ac_temp.y -= dW.y*(p_temp + visco)*m;

				p_temp = 0;
				visco = 0;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells
			real dtcv_tp = h/(c0 + dt_visco);
			if(dtcv_tp > dtcv)
			{
					#pragma critical
					{
					if(dtcv_tp > dtcv) dtcv = dtcv_tp;
					}
			}

			//Assign acceleration to particle
			particles.data.a[ai] = ac_temp - particles.data_const.graviy;

			//Reset temp. variables
			ac_temp = {0., 0.};

		} // cycle over particles in active cell



	} // cycle over cells
	particles.special.dtcv_temp = dtcv;

	#pragma omp barrier

} // function


