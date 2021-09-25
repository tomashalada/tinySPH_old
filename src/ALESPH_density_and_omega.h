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

		//Temp. data, temp. variable
		real drho_temp = 0;
		real domega_temp = 0;
		real dt_temp = 0;
		real gamma = 0;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real delta = particles.data_const.delta;
		real c0 = particles.data_const.cs;
		real dp = particles.data_const.dp;


		double *kernel;

		for(int i = 0; i < particles.cells[zz].np; i++)
		{

			/* Debug */
			//std::cout << "DENSITY -> PARICLE AC LOOP: Particle idx: " << particles.cells[zz].cp[i] << std::endl;

			if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}
			if(particles.data.part_type[particles.cells[zz].cp[i]] == wall){continue;}

			realvec ar; //position of ative particle
			realvec av; //position of ative particle
			real arho; //actual particle densiy
			real aomega; //actual particle volume

			//Load data of actual particle
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}

				/* Debug */
				//std::cout << "DENSITY -> NEIGHBOUR CELLS LOOP: Number of particles in actual cell: " << particles.cells[cl].np << std::endl;

				for(int n = 0; n < particles.cells[cl].np; n++)
				{



				//if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
				//			particles.data.part_type[particles.cells[cl].cp[n]] == wall) {continue;}


				/* Debug */
				//std::cout << "DENSITY -> PARTICLES IN NEIGHBOUR CELL: Particle ID: " << particles.cells[cl].cp[n] << std::endl;

				realvec nr; //position of ative particle
				realvec nv; //position of ative particle
				real nrho; //actual particle densiy
				real nomega; //neighbour particle

				real drs; //dr size
				realvec dW; //smoothing function gradient
				real dvdW; //vect(v) \cdot \nabla W
				realvec dr, dv;

				real W = 0;

				//dif. term
				realvec Psi;
				real PsidW;

				//Load data of neighbour particle
				nr = particles.data.r[particles.cells[cl].cp[n]];
				nv = particles.data.v[particles.cells[cl].cp[n]];
				nrho = particles.data.rho[particles.cells[cl].cp[n]];
				nomega = particles.special.omega[particles.cells[cl].cp[n]];

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));

				/* get kernel values
				double *Wendland_kernel(double r, double h) */
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];
				dW = dr * kernel[1]; // <--- check this, if its ok
				dvdW = dv.x*dW.x + dv.y*dW.y;

				//std::cout << "Normal loading... " << std::endl;

				if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				{

					realvec nb = {0., 0.};
					nb = particles.special.n[particles.cells[cl].cp[n]];

					real dvn = 0;
					dvn = dv.x*nb.x + dv.y*nb.y;

					drho_temp -=  nrho*dvn*W*dp;
					domega_temp -= dvn*W*dp;

				}
				else
				{

					drho_temp += dvdW*m;
					domega_temp += dvdW*nomega;

				}
				//std::cout << "Normal loaded. DONE. n = [ " << nb.x << "," << nb.y << " ]drho_temp: " << drho_temp << " nrhp: " << nrho << std::endl;

				gamma += W*m/nrho;



				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			//Assign density derivative to particle

			if(gamma == 0)
			{

				//particles.data.drho[particles.cells[zz].cp[i]] = 1000.;
				particles.data.drho[particles.cells[zz].cp[i]] = 0.;
				particles.special.domega[particles.cells[zz].cp[i]] = 0.;

			}
			else
			{

				particles.data.drho[particles.cells[zz].cp[i]] = drho_temp/gamma; //normalize
				particles.special.domega[particles.cells[zz].cp[i]] = domega_temp*aomega;

			}
			drho_temp = 0;
			dt_temp = 0;
			gamma = 0;
			domega_temp = 0;

			//std::cout << "..assigned " << std::endl;

			/* Debug */
			//std::cout << "DENSITY -> Assign drho_temp was SUCCESSFUL. " << drho_temp << std::endl;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


