#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"

#include "ALESPH_riemann_solver.h"
#include "_testMUSCL_schemes.h"
//Include state equation:
#include "SPH_pressure.h"

void Compute_Forces
(Particle_system &particles, real dt)
{

	//Prepare cells to find interacting particles
	unsigned int ncx, ncy;
	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	//#pragma omp parallel for
	#pragma omp parallel for schedule(dynamic, 3)
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
		realvec ac_sum = {0.,0.}; //acceleration sum

		real domega_sum = 0; //d(omega) sum
		real domegarho_sum = 0; //d(omega*rho) sum
		realvec domegav_sum = {0., 0.}; //d(omega*v) sum
		real gamma_sum = 0; //normalization factor

		//realvec gradrho_sum = {0., 0.};
		//realvec gradvx_sum = {0., 0.};
		//realvec gradvy_sum = {0., 0.};

		const real h = particles.data_const.h;
		const real m = particles.data_const.m;
		const real c0 = particles.data_const.cs;
		const real delta = particles.data_const.delta;
		const real avisco = particles.data_const.avisc;
		const real rho0 = particles.data_const.rho0;
		const realvec gravity = particles.data_const.graviy;
		const real dp = particles.data_const.dp;

		double *kernel;

		real dt_visco = 0;
		real dtacc = 0;

		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}

			//Load data of actual particle
			idx ai = particles.cells[zz].cp[i]; //actual particle index

			const realvec ar = particles.data.r[ai];; //position of ative particle
			const realvec av = particles.data.v[ai]; //position of ative particle
			const real arho = particles.data.rho[ai]; //actual particle densiy
			const real ap = particles.data.p[ai]; //actual particle density

			const real ao = particles.special.omega[ai]; //actual particle omega
			const real aorho = particles.special.omegarho[ai]; //actual particle omega*v
			const realvec aov = particles.special.omegav[ai]; //actual particle omega*v

			//MUSCL
			const realvec gradarhoF = particles.special.gradrhoF[ai];
			const realvec gradavxF = particles.special.gradvxF[ai];
			const realvec gradavyF =particles.special.gradvyF[ai];

			const realvec gradarhoB = particles.special.gradrhoB[ai];
			const realvec gradavxB = particles.special.gradvxB[ai];
			const realvec gradavyB =particles.special.gradvyB[ai];

			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

				//Load data of neighbour particle
				idx ni = particles.cells[cl].cp[n]; //actual particle index
				if(ai == ni){continue;}

				const realvec nr = particles.data.r[ni]; //position of ative particle
				const realvec nv = particles.data.v[ni];; //position of ative particle
				const real nrho = particles.data.rho[ni];; //actual particle densiy
				const real np = particles.data.p[ni];

				const real no = particles.special.omega[ni]; //actual particle omega
				const real norho = particles.special.omegarho[ni]; //actual particle omega*v
				const realvec nov = particles.special.omegav[ni]; //actual particle omega*v

				//MUSCL
				const realvec gradnrhoF = particles.special.gradrhoF[ni];
				const realvec gradnvxF = particles.special.gradvxF[ni];
				const realvec gradnvyF = particles.special.gradvyF[ni] ;

				const realvec gradnrhoB = particles.special.gradrhoB[ni];
				const realvec gradnvxB = particles.special.gradvxB[ni];
				const realvec gradnvyB = particles.special.gradvyB[ni] ;

				//Interation variables
				realvec dr, dv; //position, velocity difference
				real drs; //dr size
				real W; //kernel
				realvec dW; //kernel gradient
				real drdv; //dot product of velocity dif. and position dif.
				real dvdW; //dot product of velocity dif. and smoothing function gradient
				real drdW; //dot product of position dif. and smoothing function gradient

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));
				drdv = dr.x*dv.x + dr.y*dv.y;

				//Assign kernel values
				kernel = Wendland_kernel(drs, h);
				W = kernel[0];
				dW = dr * kernel[1];
				dvdW = dv.x*dW.x + dv.y*dW.y;
				drdW = dr.x*dW.x + dr.y*dW.y;

				// ========== MUSCL APPROX. ========== //

				// ========== SOLVE RIEMANN PROBLEM ========== //
				realvec drn = dr*(1/drs)*(-1); //pair normal vector

				//----- 1st order ----- //
				//1st order: real vL = av.x * drn.x + av.y * drn.y; //v-left state
				//1st order: real vR = nv.x * drn.x + nv.y * drn.y; //v-right state

				//1st order: real rhoL = arho;
				//1st order: real rhoR = nrho;

				//----- 2nd order, MUSCL ----- //
				//2nd order: //with limiter: realvec avL, nvR;
				//2nd order: //with limiter: real rhoL = MUSCLleft(arho, gradarho, dr, minmod_r(arho, nrho, gradarho, dr));
				//2nd order: //with limiter: real rhoR = MUSCLright(nrho, gradnrho, dr, minmod_r(arho, nrho, gradnrho, dr));
				//2nd order: //with limiter: avL.x = MUSCLleft(av.x, gradavx, dr, minmod_r(avL.x, nvR.x, gradavx, dr));
				//2nd order: //with limiter: avL.y = MUSCLleft(av.y, gradavy, dr, minmod_r(avL.x, nvR.x, gradnvx, dr));
				//2nd order: //with limiter: nvR.x = MUSCLright(nv.x, gradnvx, dr, minmod_r(avL.y, nvR.y, gradavy, dr));
				//2nd order: //with limiter: nvR.y = MUSCLright(nv.y, gradnvy, dr, minmod_r(avL.y, nvR.y, gradnvy, dr));

				realvec _avL, _avR, _nvL, _nvR, avL, nvR;

				real arhoL = MUSCL_superBee_L(arho, gradarhoF, gradarhoB, dr);
				real arhoR = MUSCL_superBee_R(arho, gradarhoF, gradarhoB, dr);
				real nrhoL = MUSCL_superBee_L(nrho, gradnrhoF, gradnrhoB, dr);
				real nrhoR = MUSCL_superBee_R(nrho, gradnrhoF, gradnrhoB, dr);

				_avL.x = MUSCL_superBee_L(av.x, gradavxF, gradavxB, dr);
				_avL.y = MUSCL_superBee_L(av.y, gradavyF, gradavyB, dr);
				_avR.x = MUSCL_superBee_R(av.x, gradavxF, gradavxB, dr);
				_avR.y = MUSCL_superBee_R(av.y, gradavyF, gradavyB, dr);

				_nvL.x = MUSCL_superBee_L(nv.x, gradnvxF, gradnvxB, dr);
				_nvL.y = MUSCL_superBee_L(nv.y, gradnvyF, gradnvyB, dr);
				_nvR.x = MUSCL_superBee_R(nv.x, gradnvxF, gradnvxB, dr);
				_nvR.y = MUSCL_superBee_R(nv.y, gradnvyF, gradnvyB, dr);

				std::array<real,3> _aW = {arho, arho*av.x, arho*av.y};
				std::array<real,3> _nW = {nrho, nrho*av.x, nrho*av.y};
				std::array<real,3> _aWL = {arhoL, arhoL*_avL.x, arhoL*_avL.y};
				std::array<real,3> _aWR = {arhoR, arhoR*_avR.x, arhoR*_avR.y};
				std::array<real,3> _nWL = {nrhoL, nrhoL*_nvL.x, nrhoL*_nvL.y};
				std::array<real,3> _nWR = {nrhoR, nrhoR*_nvR.x, nrhoR*_nvR.y};

				std::array<real,3> aWL = MUSCL_final(_aWL, _aWR, _aW, dt, h, rho0, c0);
				std::array<real,3> nWR = MUSCL_final(_nWL, _nWR, _nW, dt, h, rho0, c0);

				real rhoL = aWL[0];
				real rhoR = nWR[0];
				avL.x = aWL[1]/aWL[0];
				avL.y = aWL[2]/aWL[0];

				nvR.x = nWR[1]/nWR[0];
				nvR.y = nWR[2]/nWR[0];

				//real rhoL = MUSCL_superBee_L(arho, gradarhoF, gradarhoB, dr);
				//real rhoR = MUSCL_superBee_R(nrho, gradnrhoF, gradnrhoB, dr);

				//avL.x = MUSCL_superBee_L(av.x, gradavxF, gradavxB, dr);
				//avL.y = MUSCL_superBee_L(av.y, gradavyF, gradavyB, dr);

				//nvR.x = MUSCL_superBee_R(nv.x, gradnvxF, gradnvxB, dr);
				//nvR.y = MUSCL_superBee_R(nv.y, gradnvyF, gradnvyB, dr);

				real vL = avL.x * drn.x + avL.y * drn.y;
				real vR = nvR.x * drn.x + nvR.y * drn.y;

				//real avgc = cs(arho, rho0, c0) + cs(nrho, rho0, c0);
				real avgc = c0;
				real avgrho = 0.5*(arho + nrho);

				//----- VFRoe ----- //
				//VFRoe: real beta0 = (np-ap)/(avgrho*avgc*avgc) + NORM(nv.x-av.x, nv.y-av.x)/avgc;
				//VFRoe: real Ma = NORM(av.x, av.y)/avgc;
				//VFRoe: real beta = MIN(1, MAX(Ma,beta0));

				//VFRoe: real vss = vVFRoeTurkel(rhoL, rhoR, vL, vR, avgrho, avgc, beta);
				//VFRoe: real rhos = densVFRoeTurkel(rhoL, rhoR, vL, vR, avgrho, avgc, beta);
				//VFRoe: real ps = Compute_Pressure2(rhos,  rho0,  c0);
				//VFRoe: realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5); //reconstruct

				//----- Roe with pressure ----- //
				//RoeWithPress: //real vss = velRiemannLinearizedwithPressure(rhoL, rhoR, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				//RoeWithPress: real vss = velRiemannLinearized(rhoL, rhoR, vL, vR, avgrho, avgc);
				//RoeWithPress: real ps = pRiemannLinearized(rhoL, rhoR, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				//RoeWithPress: //real ps = pRiemannLinearizedWithLimiter(rhoL, rhoR, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				//RoeWithPress: real rhos = Pressure_to_density(ps, rho0, c0);
				//RoeWithPress: realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5); //reconstruct

				//----- Roe with density ----- //
				//real vss = velRiemannLinearizedwithPressure(rhoL, rhoR, vL, vR, ap, np, 0.5*(arho + nrho), c0);
				real vss = velRiemannLinearized(rhoL, rhoR, vL, vR, avgrho, avgc);
				real rhos = densRiemannLinearized(rhoL, rhoR, vL, vR, avgrho, avgc);
				//real rhos = densRiemannLinearizedWithLimiter(arho, nrho, vL, vR, 0.5*(arho + nrho), avgc);
				real ps = Compute_Pressure(rhos,  rho0,  c0);
				realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vL + vR)*0.5); //reconstruct

				// ========== ALESPH HELP VARIABLES ========== //

				//realvec vz = (av + nv)*0.5;
				realvec vz = av;
				real vsdW = (vs.x - vz.x)*dW.x + (vs.y - vz.y)*dW.y;

				// ========== MOMENTUM EQN ========== //
				real pR = 0; //wall pressure

				if(particles.data.part_type[ni] == wall)
				{
					realvec nb = particles.special.n[ni]*(-1);
				 real	drn = dr.x*nb.x + dr.y*nb.y;

					realvec rwf = dr;

					rhos = arho - ((0 - av.x)*nb.x + (0 - av.y)*nb.y)*arho/c0;
					ps = Compute_Pressure(rhos, rho0, c0)*1;

					domegav_sum.x -= 2*W*dp*( ps * nb.x);
					domegav_sum.y -= 2*W*dp*( ps * nb.y);
				}
				else
				{
					domegav_sum.x -= 2*no*((rhos*vs.x*(vs.x - vz.x) + ps) * dW.x + (rhos*vs.x*(vs.y - vz.y)) * dW.y);
					domegav_sum.y -= 2*no*((rhos*vs.y*(vs.x - vz.x)) * dW.x + (rhos*vs.y*(vs.y - vz.y) + ps) * dW.y);
				}

				// ========== CONTINUITY EQN ========== //
				realvec dvs = av - vs;

				if(particles.data.part_type[ni] == wall)
				{
					realvec nb = particles.special.n[ni]*(-1);
					real dvn = dv.x*nb.x + dv.y*nb.y;

					domega_sum -= dvn*W*dp;
				}
				else
				{
					vz = av;
					vsdW = (vs.x - vz.x)*dW.x + (vs.y - vz.y)*dW.y;

					domega_sum -= dvdW*no;
					domegarho_sum -= 2 * rhos * no * vsdW;
				}

				// ==================================== //

				gamma_sum += W*m/nrho;

				} // cycle over particles in neighbour cells
			} // cycle over neighbour cells

			//Assign acceleration to particle
			if(gamma_sum == 0)
			{
			particles.special.domega[ai] = 0.;
			particles.special.domegarho[ai] = 0.;
			particles.special.domegav[ai] = {0., 0.};
			}
			else
			{
			particles.special.domega[ai] = domega_sum*ao;
			particles.special.domegarho[ai] = domegarho_sum*ao;
			particles.special.domegav[ai] = domegav_sum*ao - particles.data_const.graviy*ao*arho;
			}


			//Reset temp. variables
			domega_sum = 0;
			domegarho_sum = 0;
			domegav_sum = {0., 0.};

		} // cycle over particles in active cell
	} // cycle over cells

	#pragma omp barrier
	std::cout << "bezimJa" << std::endl;

} // function


