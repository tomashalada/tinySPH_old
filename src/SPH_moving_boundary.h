#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

void Moving_boundary
(Particle_system &particles, real dp)
{

		#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] != mov_a){continue;}
			particles.data.r[p].x += dp*0.025;

		}

}
