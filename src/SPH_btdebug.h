#include "SPH_defs.h"
#include "SPH_particle_system.h"

void RemoveParticlesOutOfDomain
(Particle_system &particles, Simulation_data simulation_data)
{
	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] != fluid){continue;}


			realvec r = particles.data.r[p];
			//if(particles.data.r[p].y < simulation_data.y_0) // change fluid particle to buffer particle
			if((r.y < (0. - eps)) || (r.y > 1.8) || (r.x < 0. - eps) || (r.x > 1.6)) // change fluid particle to buffer particle
			{

				std::cout << " ** DEBUG **  Particle on position r = [" << r.x << "," << r.y << "]  was excluded. " << std::endl;
				particles.Remove_particle(p); //remove particle outgoing from buffer


			}
	}
}

