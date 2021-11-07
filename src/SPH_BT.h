#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

#include "SPH_kernel_approx.h"

//mDBC density aproximation for boundary particle
void BT_compute_density_BTGeo
(Particle_system &particles, Simulation_data simulation_data)
{


	#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] == wall)
		{

			realvec r;
			realvec gn = {0., 0.};
			realvec gnd = {0., 0.};

			r = particles.data.r[p];


			/* DEBUG */
			//if((r.x == 0.70 && r.y == -0.02) || (r.x == 0.44 && r.y == -0.04))
			//if((r.x > 0.8-eps && r.x < 0.9+eps) && (r.y == 0.02))

			//particles.data.rho[p] = Kernel_density_approximation_mDBC(particles, simulation_data, gn);
			//if(particles.data.rho[p] < particles.data_const.rho0)
			//particles.data.rho[p] = particles.data_const.rho0;

			//particles.data.rho[p] = Kernel_density_approximation_MATRIX_mDBC(particles, simulation_data, gn, gnd);
			//particles.data.rho[p] = Kernel_density_approximation_BT(particles, simulation_data, gn, gnd);
			particles.data.rho[p] = Kernel_density_approximation_BT_update(particles, simulation_data, r, gnd);
			if(particles.data.rho[p] < 1000.){particles.data.rho[p] = 1000.;}

			if((r.x > 0.0 - eps) && (r.x < 0.0 + eps) && (r.y > 0.7 - eps) && (r.y < 0.7 + eps))
			{
				std::cout << "WALL PARTICLE POSITION: wr: [" << r.x << "," << r.y << "] rho: " << particles.data.rho[p] << std::endl;
			}

		} // if function - check part type

	} // cycle over particles

}

