#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "ALESPH_riemann_solver.h"


void RSPHCompute_Density
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
		real dt_temp = 0;
		real gamma = 0;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real delta = particles.data_const.delta;
		real c0 = particles.data_const.cs;
		real dp = particles.data_const.dp;

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
			real ap;

			//Load data of actual particle
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];
			ap = particles.data.p[particles.cells[zz].cp[i]];

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}


				for(int n = 0; n < particles.cells[cl].np; n++)
				{


				if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
							particles.data.part_type[particles.cells[cl].cp[n]] == wall) {continue;}

					if(particles.cells[zz].cp[i] == particles.cells[cl].cp[n]){continue;}

				/* Debug */
				//std::cout << "DENSITY -> PARTICLES IN NEIGHBOUR CELL: Particle ID: " << particles.cells[cl].cp[n] << std::endl;

				realvec nr; //position of ative particle
				realvec nv; //position of ative particle
				real nrho; //actual particle densiy
				real np;

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
				np = particles.data.p[particles.cells[cl].cp[n]];

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
				//RSPH
				realvec drn = dr*(1/drs)*(-1);
				real vL = av.x * drn.x + av.y * drn.y;
				real vR = nv.x * drn.x + nv.y * drn.y;
				real vss = velRiemannLinearizedwithPressure(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				real ps = pRiemannLinearized(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5);

				realvec dvs = av - vs;
				real dvsdW = dvs.x*dW.x + dvs.y*dW.y;



				 if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				 {

				 	realvec nb = {0., 0.};
				 	nb = particles.special.n[particles.cells[cl].cp[n]];
				 	//nb = {0., 1.};

						vL = (-1)*( av.x * nb.x + av.y * nb.y );
						vR = -vL;

						realvec rwf;
						realvec rw;
						rw.x = nr.x;
						rw.y = 0.;

						//rwf = dr;
						rwf = ar -rw;

						real pR = ap + arho*(-1)*(particles.data_const.graviy.y*rwf.y);

				 	//dvn = nv.x*nb.x + nv.y*nb.y;

						real vss = velRiemannLinearizedwithPressure(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
						real ps = pRiemannLinearized(arho, nrho, vL, vR, ap, pR, 0.5*(arho + nrho), c0);
						realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5);

						realvec dvs = av - vs;
						real dvsdW = dvs.x*dW.x + dvs.y*dW.y;

				 	real dvn = 0;
				 	//dvn = dv.x*nb.x + dv.y*nb.y; //toto je ok
				 	dvn = dvs.x*nb.x + dvs.y*nb.y;

				 	//drho_temp -= arho*dvn*W*dp;
				 	//drho_temp -=  nrho*dvn*W*dp;
						drho_temp += dvsdW*m;
				 	//std::cout << "Normal loaded. DONE. n = [ " << nb.x << "," << nb.y << " ] dv = [ " << dv.x << "," << dv.y << " ] drho_temp: " << drho_temp << " nrhp: " << nrho << std::endl;

				 }
				 else
				 {

					//drho_temp += dvdW*m;
					drho_temp += dvsdW*m;

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

			gamma = 1;
				particles.data.drho[particles.cells[zz].cp[i]] = 2.0*drho_temp/gamma; //normalize

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


