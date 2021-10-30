#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_pressure.h"


/* Completely wrong, copletely bad, copletely ... */
void Integrate_LeapFrog_partOne_ALESPH
(Particle_system &particles, double dt)
{

	#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] == fluid)
		{

			//particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt * 0.5;
			particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt * 0.5;
			particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);

			//std::cout << "[INTEGRATE_PART_ONE] >> V: [" << particles.data.v[p].x << "," << particles.data.v[p].y  << "] OmegaV: [" << particles.special.omegav[p].x << "," << particles.special.omegav[p].y << "] Omega: " << particles.special.omega[p] << " dOmegaA: [" << particles.special.omegaa[p].x << "," << particles.special.omegaa[p].x << "] dOmega: " << particles.special.domega[p] << std::endl;

		} // if function - check part type


	} // cycle over particles

} // function

void Integrate_LeapFrog_partTwo_ALESPH
(Particle_system &particles, double dt, int step)
{

	if(step == 1)
	{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] == fluid)
			{


				//particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt * 0.5;
				particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt * 0.5;
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);

				//particles.data.v[p] = particles.data.v[p] + particles.data.a[p] * dt * 0.5;
				particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt * 0.5;

			} // if function - check part type

		} // loop over particles

	} // if function - check step
	else
	{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] == fluid)
			{

				//particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt;
				particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt;
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p] * particles.special.omega[p]);

			//	particles.data.v[p] = particles.data.v[p] + particles.data.a[p] * dt;
				particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;

			} // if function - check part type

		} // loop over particles

	} // if function - check step

} // function


void Integrate_density_compute_pressure
(Particle_system &particles, double dt)
{

		//#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			/* Test for outlet */
			if(particles.data.part_type[p] == outlet){continue;}

			/* Test for outlet */
			if(particles.data.part_type[p] == outletf){continue;}

			/* ALESPHTEST */
			if(particles.data.part_type[p] != fluid){continue;}

		 particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt;
			particles.special.omegarho[p] = particles.special.omegarho[p] + particles.special.domegarho[p] * dt;
			//	particles.data.v[p] = particles.special.omegav[p] / particles.special.omega[p];

			//particles.data.rho[p] += particles.data.drho[p]*dt;
			particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);

			if(particles.data.rho[p] > 1500. || particles.data.rho[p] < 0.)
			{std::cout << " ** DEBUG **  Particle on position rho = [" << particles.data.rho[p] << "]  was excluded. " << std::endl;}

			//std::cout << " rho,pINTEGRATE ... rho: " << particles.data.rho[p] << " omegarho: " << particles.special.omegarho[p] << " omega: " <<  particles.special.omega[p] << " domegarho: " << particles.special.domegarho[p] << std::endl;

			/* Debug */
			//if(abs(particles.data.drho[p]) > 0.0005)
			//{

			//	std::cout << "SIMULATION -> RUN -> PRESSURE: drho non-zero: " << particles.data.drho[p] << std::endl;

			//}

			particles.data.p[p] = Compute_Pressure2(particles.data.rho[p], particles.data_const.rho0, particles.data_const.cs);
			//std::cout << "[DENSITY&PRESSURE] >> Rho: " << particles.data.rho[p] << " dOmegaRho: " << particles.special.domegarho[p] << " Omega:" << particles.special.omega[p] << std::endl;


			/* Debug */
			//std::cout << "SIMULATION -> RUN -> PRESSURE: p: " << particles.data.p[p] << std::endl;
			if(particles.data.p[p] < 0){
				particles.data.p[p] = 0;
			}

		}

} // function

void Integrate_ALESPH
(Particle_system &particles, double dt, int step)
{


		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] == fluid)
			{

				particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt;
				particles.special.omegarho[p] = particles.special.omegarho[p] + particles.special.domegarho[p] * dt;
				particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt;

				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p] * particles.special.omega[p]);
				particles.data.rho[p] = particles.special.omegarho[p] / (particles.special.omega[p]);
				particles.data.p[p] = Compute_Pressure(particles.data.rho[p], particles.data_const.rho0, particles.data_const.cs);
				particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;

			}

		}

} // function

