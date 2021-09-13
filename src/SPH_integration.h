#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_pressure.h"

#include "SPH_acceleration.h" /* MOVE OUT, its because outlet */

/* Completely wrong, copletely bad, copletely ... */
void Integrate_LeapFrog_partOne
(Particle_system &particles, double dt)
{

	#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] == fluid)
		{

			particles.data.v[p] = particles.data.v[p] + particles.data.a[p] * dt * 0.5;

		} // if function - check part type


	} // cycle over particles

} // function

void Integrate_LeapFrog_partTwo
(Particle_system &particles, double dt, int step)
{

	if(step == 1)
	{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] == fluid)
			{

				particles.data.v[p] = particles.data.v[p] + particles.data.a[p] * dt * 0.5;
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

			particles.data.v[p] = particles.data.v[p] + particles.data.a[p] * dt;
			particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;

			} // if function - check part type

		} // loop over particles

	} // if function - check step

} // function


void Integrate_density_compute_pressure
(Particle_system &particles, double dt)
{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			/* Test for outlet */
			if(particles.data.part_type[p] == outlet){continue;}

			/* Test for outlet */
			if(particles.data.part_type[p] == outletf){continue;}



			particles.data.rho[p] += particles.data.drho[p]*dt;

			/* Debug */
			//if(abs(particles.data.drho[p]) > 0.0005)
			//{

			//	std::cout << "SIMULATION -> RUN -> PRESSURE: drho non-zero: " << particles.data.drho[p] << std::endl;

			//}

			particles.data.p[p] = Compute_Pressure(particles.data.rho[p], particles.data_const.rho0, particles.data_const.cs);


			/* Debug */
			//std::cout << "SIMULATION -> RUN -> PRESSURE: p: " << particles.data.p[p] << std::endl;
			if(particles.data.p[p] < 0){
				//particles.data.p[p] = 0;
			}

		}

} // function

