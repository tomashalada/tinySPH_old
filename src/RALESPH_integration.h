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
			particles.special.omega_o[p] = particles.special.omega[p];
			particles.special.omegarho_o[p] = particles.special.omegarho[p];
			particles.special.omegav_o[p] = particles.special.omegav[p];

		 particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt * 0.5;
			particles.special.omegarho[p] = particles.special.omegarho[p] + particles.special.domegarho[p] * dt * 0.5;
			particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt * 0.5;

			//particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);
			//particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);

			particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);
			particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);
		}

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

				particles.special.omega_o[p] = particles.special.omega[p];
				particles.special.omegarho_o[p] = particles.special.omegarho[p];
				particles.special.omegav_o[p] = particles.special.omegav[p];

		 	particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt * 0.5;
				particles.special.omegarho[p] = particles.special.omegarho[p] + particles.special.domegarho[p] * dt * 0.5;
				particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt * 0.5;

				//particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);
				//particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);

				particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);

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

		 	//to tu bylo: particles.special.omega[p] = particles.special.omega_o[p] + particles.special.domega[p] * dt * 0.5;
				//to tu bylo: particles.special.omegarho[p] = particles.special.omegarho_o[p] + particles.special.domegarho[p] * dt * 0.5;
				//to tu bylo: particles.special.omegav[p] = particles.special.omegav_o[p] + particles.special.domegav[p] * dt * 0.5;

				//to tu bylo: particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);
				//to tu bylo: particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);

				//to tu bylo: particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;

				particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;

		 	particles.special.omega[p] = particles.special.omega_o[p] + particles.special.domega[p] * dt;
				particles.special.omegarho[p] = particles.special.omegarho_o[p] + particles.special.domegarho[p] * dt;
				particles.special.omegav[p] = particles.special.omegav_o[p] + particles.special.domegav[p] * dt;

				//particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);
				//particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);

				particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);


			} // if function - check part type
		} // loop over particles
	} // if function - check step
} // function


void Integrate_Verlet
 (Particle_system &particles, double dt)
 {

 	#pragma omp parallel for schedule(static)
 	for(int p = 0; p < particles.np; p++)
 	{

 		//if(particles.data.part_type[p] == fluid)
 		{

		 	particles.special.omega[p] = particles.special.omega_oo[p] + particles.special.domega[p] * dt * 2;
				particles.special.omegarho[p] = particles.special.omegarho_oo[p] + particles.special.domegarho[p] * dt * 2;
				particles.special.omegav[p] = particles.special.omegav_oo[p] + particles.special.domegav[p] * dt * 2;

 			//particles.data.v[p] = particles.data.v_oo[p] + particles.data.a[p] * dt * 2 ;
				particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);

 			particles.data.r[p] = particles.data.r[p] + particles.data.v_o[p] * dt + particles.special.domegav[p] * dt *dt * 0.5/particles.special.omega[p];



 			//Move time layers
 			particles.data.v_oo[p] = particles.data.v_o[p];
 			particles.data.v_o[p] = particles.data.v[p];

 			particles.data.rho_oo[p] = particles.data.rho_o[p];
 			particles.data.rho_o[p] = particles.data.rho[p];

				particles.special.omega_oo[p] = particles.special.omega_o[p];
				particles.special.omegarho_oo[p] = particles.special.omegarho_o[p];
				particles.special.omegav_oo[p] = particles.special.omegav_o[p];

				particles.special.omega_o[p] = particles.special.omega[p];
				particles.special.omegarho_o[p] = particles.special.omegarho[p];
				particles.special.omegav_o[p] = particles.special.omegav[p];

 		} // if function - check part type
 	} // cycle over particles
} // function

 void Integrate_Euler
 (Particle_system &particles, double dt)
 {

 	#pragma omp parallel for schedule(static)
 	for(int p = 0; p < particles.np; p++)
 	{

 		//if(particles.data.part_type[p] == fluid)
 		{

		 	particles.special.omega[p] = particles.special.omega[p] + particles.special.domega[p] * dt;
				particles.special.omegarho[p] = particles.special.omegarho[p] + particles.special.domegarho[p] * dt;
				particles.special.omegav[p] = particles.special.omegav[p] + particles.special.domegav[p] * dt ;

				particles.data.rho[p] = particles.special.omegarho[p] /( particles.special.omega[p]);
				particles.data.v[p] = particles.special.omegav[p] / (particles.data.rho[p]*particles.special.omega[p]);

 			particles.data.r[p] = particles.data.r[p] + particles.data.v_o[p] * dt + particles.special.domegav[p] * dt *dt * 0.5/particles.special.omega[p];


 			//Move time layers
 			particles.data.v_oo[p] = particles.data.v_o[p];
 			particles.data.v_o[p] = particles.data.v[p];

 			particles.data.rho_oo[p] = particles.data.rho_o[p];
 			particles.data.rho_o[p] = particles.data.rho[p];

				particles.special.omega_oo[p] = particles.special.omega_o[p];
				particles.special.omegarho_oo[p] = particles.special.omegarho_o[p];
				particles.special.omegav_oo[p] = particles.special.omegav_o[p];

				particles.special.omega_o[p] = particles.special.omega[p];
				particles.special.omegarho_o[p] = particles.special.omegarho[p];
				particles.special.omegav_o[p] = particles.special.omegav[p];


 		} // if function - check part type
 	} // cycle over particles
 } // function
