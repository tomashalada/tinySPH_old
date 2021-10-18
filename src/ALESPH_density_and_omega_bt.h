#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "ALESPH_riemann_solver.h"


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
		real domega_temp = 0;
		real domegarho_temp = 0;
		real gamma = 0;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real delta = particles.data_const.delta;
		real rho0 = particles.data_const.rho0;
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

			//if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}

			realvec ar; //position of ative particle
			realvec av; //position of ative particle
			real arho; //actual particle densiy
			real aomega; //actual particle volume
			real ap;

			// *** MUSCL **** //
			realvec gradarho;
			realvec gradavx;
			realvec gradavy;

			//Load data of actual particle
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];
			aomega = particles.special.omega[particles.cells[zz].cp[i]];
			ap = particles.data.p[particles.cells[zz].cp[i]];

			// *** MUSCL **** //
			gradarho = particles.special.gradrho[particles.cells[zz].cp[i]];
			gradavx = particles.special.gradvx[particles.cells[zz].cp[i]];
			gradavy =particles.special.gradvy[particles.cells[zz].cp[i]];

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}

				/* Debug */
				//std::cout << "DENSITY -> NEIGHBOUR CELLS LOOP: Number of particles in actual cell: " << particles.cells[cl].np << std::endl;

				for(int n = 0; n < particles.cells[cl].np; n++)
				{



				if(particles.data.part_type[particles.cells[zz].cp[i]] == wall &&
							particles.data.part_type[particles.cells[cl].cp[n]] == wall) {continue;}
				if(particles.cells[zz].cp[i] == particles.cells[cl].cp[n]){continue;}


				/* Debug */
				//std::cout << "DENSITY -> PARTICLES IN NEIGHBOUR CELL: Particle ID: " << particles.cells[cl].cp[n] << std::endl;

				//Lead nb particle data
				realvec nr; //position of ative particle
				realvec nv; //position of ative particle
				real nrho; //actual particle densiy
				real nomega; //neighbour particle
				real np; //neighbour pressure

				// *** MUSCL **** //
				realvec gradnrho;
				realvec gradnvx;
				realvec gradnvy;

				//Load data of neighbour particle
				nr = particles.data.r[particles.cells[cl].cp[n]];
				nv = particles.data.v[particles.cells[cl].cp[n]];
				nrho = particles.data.rho[particles.cells[cl].cp[n]];
				nomega = particles.special.omega[particles.cells[cl].cp[n]];
				np = particles.data.p[particles.cells[cl].cp[n]];

				// *** MUSCL **** //
				gradnrho = particles.special.gradrho[particles.cells[cl].cp[n]];
				gradnvx = particles.special.gradvx[particles.cells[cl].cp[n]];
				gradnvy =particles.special.gradvy[particles.cells[cl].cp[n]] ;

				real drs; //dr size
				realvec dW; //smoothing function gradient
				real dvdW; //vect(v) \cdot \nabla W
				realvec dr, dv;
				real W;

				//ALESPH
				real vr, vl; //left and right state
				realvec drn; //unit pair vector
				real rhos; //rho^star, Riemann problem solution
				real vss; //v^star, Riemann problem solution
				realvec vs; //v^star, Riemann problem solution
				real vsdW;

				real arhoL;
				real nrhoR;
				realvec avL;
				realvec nvR;

				//dif. term
				realvec Psi;
				real PsidW;

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				//dv = nv - av;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
				drn = dr*(1/drs)*(-1); //unit pair vector

				// Get kernels values
				//double *Wendland_kernel(double r, double h)
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];
				dW = dr * kernel[1];

				// *** MUSCL **** //
				real param_r = 1;
				real gradn;

				arhoL = MUSCLleft(arho, gradarho, dr, param_r);
				nrhoR = MUSCLright(nrho, gradnrho, dr, param_r);

				avL.x = MUSCLleft(av.x, gradavx, dr, param_r);
				avL.y = MUSCLleft(av.y, gradavy, dr, param_r);

				nvR.x = MUSCLright(nv.x, gradnvx, dr, param_r);
				nvR.y = MUSCLright(nv.y, gradnvy, dr, param_r);

				//Riemann problem
				vl = avL.x * drn.x + avL.y * drn.y;
				vr = nvR.x * drn.x + nvR.y * drn.y;

				//real avgc = cs(arho, rho0, c0) + cs(nrho, rho0, c0);
				real avgc = c0;
				rhos = densRiemannLinearized(arhoL, nrhoR, vl, vr, 0.5*(arhoL + nrhoR), avgc);
				//rhos = densRiemannLinearizedWithLimiter(arhoL, nrhoR, vl, vr, 0.5*(arhoL + nrhoR), avgc);
				vss = velRiemannLinearized(arhoL, nrhoR, vl, vr, 0.5*(arhoL + nrhoR), avgc);

				// //Riemann problem
				// vl = av.x * drn.x + av.y * drn.y;
				// vr = nv.x * drn.x + nv.y * drn.y;

				// //real avgc = cs(arho, rho0, c0) + cs(nrho, rho0, c0);
				// real avgc = c0;
				// rhos = densRiemannLinearized(arho, nrho, vl, vr, 0.5*(arho + nrho), avgc);
				// vss = velRiemannLinearized(arho, nrho, vl, vr, 0.5*(arho + nrho), avgc);

				//vss = velRiemannLinearizedwithPressure(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), avgc);

				//***experiment***
				//real ps = pRiemannLinearizedWithLimiter(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), avgc);
				//rhos = ps/(c0*c0) + rho0;
				//***experiment***

				vs = drn * vss + ((av + nv)*0.5 - drn*(vl + vr)*0.5);

				//realvec vz = (av + nv)*0.5;
				//realvec vz = (av + nv)/2;
				realvec vz = av;

				vsdW = (vs.x - vz.x)*dW.x + (vs.y - vz.y)*dW.y;

				//std::cout << "[DENSITY] >> aV: [" << av.x << "," << av.y << "] nV:" << av.x << "," << av.y << "] vl:" << vl << " vr: " << vr << " Vs:" << av.x << "," << av.y << "] Vss: " << vss << " aRho: "<< arho << " nRho: " << nrho << " Rhos:" << rhos << std::endl;

				//ALE equations
				if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				{

					realvec nb = {0., 0.};
					nb = particles.special.n[particles.cells[cl].cp[n]]*(-1);

					real dvn = 0;
					real vsn = 0;
					dvn = dv.x*nb.x + dv.y*nb.y;
					vsn = dv.x*nb.x + dv.y*nb.y;

					domega_temp -= dvn*W*dp;
					//domegarho_temp -= 2 * rhos * vsn * W * dp; //zero due BC

				}
				else
				{

					domega_temp -= dvdW*nomega; //tdy muze byt +
					domegarho_temp -= 2 * rhos * nomega * vsdW;

				}

				gamma += W*m/nrho;
				//gamma += W;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			//Assign density derivative to particle
			if(gamma == 0)
			{

				particles.special.domega[particles.cells[zz].cp[i]] = 0.;
				particles.special.domegarho[particles.cells[zz].cp[i]] = 0.;

			}
			else
			{

				particles.special.domega[particles.cells[zz].cp[i]] = domega_temp*aomega;
				particles.special.domegarho[particles.cells[zz].cp[i]] = domegarho_temp*aomega;

			}

				//std::cout << "DENSITY >> Write/ aomega: " << aomega <<" domega: " << domega_temp <<  " domegarho_temp: " << domegarho_temp <<  std::endl;

				//std::cout << "[DENSITY2] >> aOmega: " << aomega << " dOmega: " << domega_temp << " dOmegaRho: " << domegarho_temp << std::endl;

					//av.y << "] nV:" << av.x << "," << av.y << "] vl:" << vl << " vr: " << vr << " Vs:" << av.x << "," << av.y << "] Vss: " << vss << " aRho: "<< arho << " nRho: " << nrho << " Rhos:" << rhos << std::endl;

			gamma = 0;
			domega_temp = 0;
			domegarho_temp = 0;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


