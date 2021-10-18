#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

void Compute_Density
(Particle_system &particles)
{

	//Prepare cells to find interacting particles
	unsigned int ncx, ncy;
	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

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
		//std::cout << "DENSITY -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

		//Temp. data, temp. variable
		real drho_temp = 0;
		real dt_temp = 0;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real delta = particles.data_const.delta;
		real c0 = particles.data_const.cs;
		real rho0 = particles.data_const.rho0;

		double *kernel;

		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			//if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == wall){continue;}

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i]; //actual particle index

			realvec ar = particles.data.r[ai];; //position of ative particle
			realvec av = particles.data.v[ai];; //position of ative particle
			real arho = particles.data.rho[ai];; //actual particle densiy


			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}

				for(int n = 0; n < particles.cells[cl].np; n++)
				{

				/*
				if(particles.data.part_type[particles.cells[cl].cp[n]] == outletf){continue;}

				if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
							particles.data.part_type[particles.cells[cl].cp[n]] == inlet) {continue;}

				if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
							particles.data.part_type[particles.cells[cl].cp[n]] == outlet) {continue;}

				if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
							particles.data.part_type[particles.cells[cl].cp[n]] == outletf) {continue;}

					*/
				//if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
				//			particles.data.part_type[particles.cells[cl].cp[n]] == wall) {continue;}


				//Load data of neighbour particle
				idx ni = particles.cells[cl].cp[n]; //actual particle index
				//if(ai == ni){continue;}

				realvec nr = particles.data.r[ni]; //position of ative particle
				realvec nv = particles.data.v[ni];; //position of ative particle
				real nrho = particles.data.rho[ni];; //actual particle densiy

				//Interation variables
				realvec dr, dv; //position, velocity difference
				real drs; //dr size
				realvec dW; //kernel gradient
				real dvdW; //dot product of velocity smoothing function gradient

				//Diffusion term variables
				realvec Psi;
				real PsidW;

				//Assign values
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));

				//Assign kernel values
				//double *Wendland_kernel(double r, double h)
				kernel = Wendland_kernel(drs, h);
				dW = dr * kernel[1];
				dvdW = dv.x*dW.x + dv.y*dW.y;

				//Assign diffusive term valuse
				//Psi.x =  (dr.x / (pow(drs,2) + eps*h*h))* 2 * (arho - nrho);
				//Psi.y =  (dr.y / (pow(drs,2) + eps*h*h))* 2 * (arho - nrho);

				 real dpH = rho0* dr.y * particles.data_const.graviy.y;
				 real cb = c0*c0*rho0 / 7.;
				 //real drhoH = rho0 * (pow(((dpH + 1)/cb),1/7) - 1);
				 real drhoH = rho0 * (pow(((dpH + 1)/cb),1/7) - 1);

			  Psi.x =  (dr.x / (pow(drs,2) + eps*h*h))* 2 * ((arho - nrho) - drhoH);
			  Psi.y =  (dr.y / (pow(drs,2) + eps*h*h))* 2 * ((arho - nrho) - drhoH);

				PsidW = Psi.x*dW.x + Psi.y*dW.y;

				//Add pair increment
				drho_temp += dvdW*m/nrho;

				//drho_temp += dvdW*m;
				//if(particles.data.part_type[ni] != wall){
				dt_temp += h*delta*c0*PsidW * m/ nrho;
				//}

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			//Assign sum to particle
			//particles.data.drho[ai] = drho_temp;
			particles.data.drho[ai] = drho_temp*arho + dt_temp;

			//Reset temp. variables
			drho_temp = 0;
			dt_temp = 0;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


