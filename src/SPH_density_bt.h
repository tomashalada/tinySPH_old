#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"


void Compute_Density_BT
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
		real gamma = 0;

		/*
		real drs; //dr size
		realvec dW; //smoothing function gradient
		real dvdW; //vect(v) \cdot \nabla W
		*/

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real delta = particles.data_const.delta;
		real c0 = particles.data_const.cs;
		real dp = particles.data_const.dp;

		/*
		realvec ar; //position of ative particle
		realvec nr; //position of neighbour particle
		realvec dr; //position diference

		realvec av; //position of ative particle
		realvec nv; //position of neighbour particle
		realvec dv; //velocity diference

		real nrho; //neighbour density
		real arho; //actual particle densiy
		*/

		double *kernel;

		/* Debug */
		//std::cout << "DENSITY -> CELLS LOOP: Number of particles in actual cell: " << particles.cells[zz].np << std::endl;

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
					//dvn = nv.x*nb.x + nv.y*nb.y;
					dvn = dv.x*nb.x + dv.y*nb.y;

					//drho_temp -= arho*dvn*W*dp;
					drho_temp -=  nrho*dvn*W*dp;
					//std::cout << "Normal loaded. DONE. n = [ " << nb.x << "," << nb.y << " ] dv = [ " << dv.x << "," << dv.y << " ] drho_temp: " << drho_temp << " nrhp: " << nrho << std::endl;

				}
				else
				{

					drho_temp += dvdW*m;

				}
				//std::cout << "Normal loaded. DONE. n = [ " << nb.x << "," << nb.y << " ]drho_temp: " << drho_temp << " nrhp: " << nrho << std::endl;

				gamma += W*m/nrho;


				/*  compute diffusive term */
				//Psi.x =  (dr.x / (pow(drs,2) + eps*h*h))* 2 * (nrho - arho);
				//Psi.y =  (dr.y / (pow(drs,2) + eps*h*h))* 2 * (nrho - arho);
				//PsidW = Psi.x*dW.x + Psi.y*dW.y;

				/* Debug */
				//if(abs(kernel[1]) > 0.0001)
				//{

				//	std::cout << "DENSITY -> kernel dW: " << kernel[1] << std::endl;

				//}
				if(dv.x*dv.x + dv.y*dv.y > 0.0001)
				{

					//std::cout << "DENSITY -> dv dot dv: " << dv.x*dv.x + dv.y*dv.y << " dvx: " << dv.x << " dvy: " << dv.y << std::endl;

				}
				//if(dvdW > 0.00001)
				//{

				//	std::cout << "DENSITY -> dvdW: " << dvdW << std::endl;

				//}

				//drho_temp += dvdW*m/nrho;
				//dt_temp += h*delta*c0*PsidW * m/ nrho;


				/* Debug */
				//std::cout << "DENSITY -> m: " << m << std::endl;
				//std::cout << "DENSITY -> PARTICLES IN NEIGHBOUR dr: [" << dr.x << ", " << dr.y << "] dv: [" << dv.x << ", " << dv.y << "]" << std::endl;
				//std::cout << "DENSITY -> PARTICLES IN NEIGHBOUR dW: [" << dW.x << ", " << dW.y << "], [" << kernel[1] << "] drs: " << drs << " q: " << fabs(drs/h) <<", dvdW: " << dvdW << " nrho: " << nrho << std::endl;


				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			//Assign density derivative to particle

			/* Debug */
			//std::cout << "DENSITY -> Assign drho_temp to array drho_temp: " << drho_temp << " arho: " << arho << " Idx to assign drho_temp: "<< particles.cells[zz].cp[i] << std::endl;

			/* Debug */
			//if(abs(drho_temp) > 0.0001)
			//{

			//	std::cout << "DENSITY -> drho non-zero: " << drho_temp << std::endl;
			//	std::cout << "DENSITY -> Assign drho_temp was SUCCESSFUL. drho_temp: " << drho_temp <<  " dt_temp: " << dt_temp << std::endl;

			//}

			/* Debug */
			//std::cout << "DENSITY -> Assign drho_temp was SUCCESSFUL. drho_temp: " << drho_temp <<  " dt_temp: " << dt_temp << std::endl;

			//particles.data.drho[particles.cells[zz].cp[i]] = drho_temp*arho;
			//particles.data.drho[particles.cells[zz].cp[i]] = drho_temp - dt_temp;
			//std::cout << "Assing = " << gamma << std::endl;

			if(gamma == 0)
			{

				//particles.data.drho[particles.cells[zz].cp[i]] = 1000.;
				particles.data.drho[particles.cells[zz].cp[i]] = 0.;

			}
			else
			{

				particles.data.drho[particles.cells[zz].cp[i]] = drho_temp/gamma; //normalize

			}
			drho_temp = 0;
			dt_temp = 0;
			gamma = 0;

			//std::cout << "..assigned " << std::endl;

			/* Debug */
			//std::cout << "DENSITY -> Assign drho_temp was SUCCESSFUL. " << drho_temp << std::endl;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


