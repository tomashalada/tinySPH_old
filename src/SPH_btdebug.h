#include "SPH_defs.h"
#include "SPH_particle_system.h"

void RemoveParticlesOutOfDomain
(Particle_system &particles, Simulation_data simulation_data)
{


	real x0 = simulation_data.x_0;
	real xm = simulation_data.x_m;
	real y0 = simulation_data.y_0;
	real ym = simulation_data.y_m;
	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] != fluid){continue;}


			realvec r = particles.data.r[p];
			//if(particles.data.r[p].y < simulation_data.y_0) // change fluid particle to buffer particle
			if((r.y < (y0 - eps)) || (r.y > ym + eps) || (r.x < x0 - eps) || (r.x > xm - eps)) // change fluid particle to buffer particle
			{

				std::cout << " ** DEBUG **  Particle on position r = [" << r.x << "," << r.y << "]  was excluded. " << std::endl;
				particles.Remove_particle(p); //remove particle outgoing from buffer


			}
	}
}

void RemoveParticlesOutOfDomain_ALE
(Particle_system &particles, Simulation_data simulation_data)
{

	real x0 = simulation_data.x_0;
	real xm = simulation_data.x_m;
	real y0 = simulation_data.y_0;
	real ym = simulation_data.y_m;
	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] != fluid){continue;}


			realvec r = particles.data.r[p];
			//if(particles.data.r[p].y < simulation_data.y_0) // change fluid particle to buffer particle
			if((r.y < (y0 - eps)) || (r.y > ym + eps) || (r.x < x0 - eps) || (r.x > xm - eps)) // change fluid particle to buffer particle
			{

				std::cout << " ** DEBUG **  Particle on position r = [" << r.x << "," << r.y << "]  was excluded. " << std::endl;
				particles.ALESPH_Remove_particle(p); //remove particle outgoing from buffer


			}
	}
}
