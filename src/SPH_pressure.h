#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

//Free surface flow state eq. ver. A
real Compute_Pressure
(real rho, real rho0, real c0)
{

	real p;

	const real gamma = 7.; //Weakly compressible state eq. Poisson constant
	const real b_const = c0*c0*rho0/gamma; //Weakly compressible state eq. constant
	const real rho_r = rho/rho0; //Weakly compressible state eq. relative density

	//Weakly compressible state equation
	p = b_const*( pow(rho_r, gamma) - 1 );

	return p;
}

//Free surface flow state eq. ver. B
real Compute_Pressure2
(real rho, real rho0, real c0)
{

	real p;

	//Weakly compressible state equation
	p =c0*c0*( rho - rho0 );

	return p;
}

//Free surface flow state eq. ver. A
real Pressure_to_density
(real p, real rho0, real c0)
{

	real rho;

	const real gamma = 7.; //Weakly compressible state eq. Poisson constant
	const real b_const = c0*c0*rho0/gamma;; //Weakly compressible state eq. constant
	real rho_r; //Weakly compressible state eq. relative density

	//Weakly compressible state equation
	rho = (pow(p / b_const + 1, 1./7) )*rho0;

	return rho;
}

real Pressure_to_density2
(real p, real rho0, real c0)
{

	real rho;

	//Weakly compressible state equation
	rho = p / (c0*c0) + rho0;

	return rho;
}

void Density_to_pressure
(Particle_system &particles)
{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] == outlet){continue;}
			if(particles.data.part_type[p] == inlet){continue;}

			if(particles.data.rho[p] > 1500. || particles.data.rho[p] < 0.)
			{
				std::cout << "=== ERROR ===  Particle density too high, rho: " << particles.data.rho[p] << " position r: [" << particles.data.r[p].x << "," << particles.data.r[p].y << "]" << std::endl;
				//particles.data.rho[p] = 1000.;
			}


			particles.data.p[p] = Compute_Pressure(particles.data.rho[p], particles.data_const.rho0, particles.data_const.cs);

			/* Test */
			//if(particles.data.p[p] < 0){ particles.data.p[p] = 0; }

		}
} // function


