#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"
#include "ALESPH_riemann_solver.h"

void RSPHCompute_Acceleration
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
		realvec ac_temp_visco = {0.,0.};
		real p_temp = 0;
		real visco = 0;
		real gamma = 0;

		real drs; //dr size
		real drdv; //vect(v) \cdot \nabla W
		realvec dW; //smoothing function gradient
		real W;
		real ndr;
		real rdW;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real dp = particles.data_const.dp;
		real c0 = particles.data_const.cs;

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
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}

			/*
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

				if(particles.cells[zz].cp[i] == particles.cells[cl].cp[n]){continue;}

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
				W = kernel[0];

				/* pressure term */
				p_temp = ap/pow(arho, 2) + np/pow(nrho,2);
				//p_temp = (ap + np)/(arho*nrho);

				/* real Artificial_Viscosity
				(real h, real drs, real drdv, real rho0, real c0, real alpha) */
				visco = Artificial_Viscosity(h, drs, drdv, particles.data_const.rho0, particles.data_const.cs, particles.data_const.avisc);


				//RSPH
				realvec drn = dr*(1/drs)*(-1);
				real vL = av.x * drn.x + av.y * drn.y;
				real vR = nv.x * drn.x + nv.y * drn.y;
				real vss = velRiemannLinearizedwithPressure(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				//real ps = pRiemannLinearized(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				real ps = pRiemannLinearizedWithLimiter(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5);

				realvec dvs = av - vs;
				real dvsdW = dvs.x*dW.x + dvs.y*dW.y;
				// <--- check this, if its ok

				//ac_temp = ac_temp - dW*(p_temp + visco)*m; //+= operator is not overloaded yet

				 if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				 {

				 	realvec nb = {0., 0.};
				 	nb = particles.special.n[particles.cells[cl].cp[n]];
				 	//nb = {0., 1.};

				 	ndr = dr.x*nb.x + dr.y*nb.y;

						vL = (-1)*( av.x * nb.x + av.y * nb.y );
						vR = -vL;

						realvec rwf;
						rwf = dr;
						real pR = ap + arho*(-1)*(particles.data_const.graviy.y*rwf.y);

						real vss = velRiemannLinearizedwithPressure(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
						//real ps = pRiemannLinearized(arho, nrho, vL, vR, ap, pR, 0.5*(arho + nrho), c0);
						real ps = pRiemannLinearizedWithLimiter(arho, nrho, vL, vR, ap, np, 0.5*(arho + nrho), c0);
						realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5);

						realvec dvs = av - vs;
						real dvsdW = dvs.x*dW.x + dvs.y*dW.y;

				 	real dvn = 0;
				 	//dvn = dv.x*nb.x + dv.y*nb.y; //toto je ok
				 	dvn = dvs.x*nb.x + dvs.y*nb.y;

						ac_temp.x = ac_temp.x - dW.x*2*(ps/(arho*nrho))*m;
						ac_temp.y = ac_temp.y - dW.y*2*(ps/(arho*nrho))*m;

				 	//ac_temp.x += p_temp * nb.x * W * dp * nrho;
				 	//ac_temp.y += p_temp * nb.y * W * dp * nrho;

				 	//// ac_temp.x += (ps/(arho*nrho)) * nb.x * W * dp * nrho;
				 	//// ac_temp.y += (ps/(arho*nrho)) * nb.y * W * dp * nrho;

				 	//ac_temp_visco.x -= dv.x * ndr * W * dp / (drs*drs);
				 	//ac_temp_visco.y -= dv.y * ndr * W * dp / (drs*drs);

				 	//ac_temp.x += (p_temp + visco) * nb.x * W * dp * nrho;
				 	//ac_temp.y += (p_temp + visco) * nb.y * W * dp * nrho;

				 }
				 else
				 {

					rdW = dr.x * dW.x + dr.y * dW.y;

					//ac_temp.x = ac_temp.x - dW.x*(p_temp + visco)*m;
					//ac_temp.y = ac_temp.y - dW.y*(p_temp + visco)*m;

					ac_temp.x = ac_temp.x - dW.x*2*(ps/(arho*nrho))*m;
					ac_temp.y = ac_temp.y - dW.y*2*(ps/(arho*nrho))*m;

					//ac_temp.x = ac_temp.x - dW.x*(p_temp )*m;
					//ac_temp.y = ac_temp.y - dW.y*(p_temp )*m;

					//ac_temp_visco.x = ac_temp_visco.x + dv.x * rdW * m / (nrho * drs * drs);
					//ac_temp_visco.y = ac_temp_visco.y + dv.y * rdW * m / (nrho * drs * drs);

				 }

				gamma += W*m/nrho;

				/* Temp, step */
				p_temp = 0; visco = 0;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			if(gamma == 0)
			{

				realvec zero = {0., 0.};
				particles.data.a[particles.cells[zz].cp[i]] = zero - particles.data_const.graviy;

			}
			else
			{

				gamma=1;
				particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma +  ac_temp_visco*2*0.001/(gamma * arho) - particles.data_const.graviy;
				//particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma +  ac_temp_visco/gamma - particles.data_const.graviy;
				//particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma   - particles.data_const.graviy;

			}
			ac_temp = {0., 0.};
			ac_temp_visco = {0., 0.};
			gamma = 0;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function


