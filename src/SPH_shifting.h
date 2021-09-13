#include "SPH_defs.h"
#include "SPH_particle_system.h"

#include "SPH_kernel_approx.h"

void ShiftParticles
(Particle_system &particles, Simulation_data simulation_data, real dt)
{

	const real Afsm = 2.;
	const real Afst = 1.5;
	const real A = 2;

	const real h = particles.data_const.h;


	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] != fluid){continue;}


		real Afsc = 0;
		realvec deltar = {0., 0.};
		realvec dC = {0., 0.};
		real divr = 0;

		realvec ar = particles.data.r[p];
		realvec av = particles.data.v[p];
		real vs = sqrt(pow(av.x, 2) + pow(av.y, 2));

		divr = Kernel_positionVec_divergence(particles, simulation_data, ar);
	 dC = Kernel_particles_concentration_gradient(particles, simulation_data, ar);
		Afsc = (divr - Afst)/(Afsm - Afst);

		if(divr - Afst < 0)
		{


			deltar = dC * ( Afsc * A * h * vs * dt);

		}
		else if((divr - Afst) == 0)
		{

			deltar = dC * ( A * h * vs * dt);

		}
		else
		{

			deltar = {0., 0.};

		}

		particles.data.r[p] = particles.data.r[p] + deltar;
		/* Debug */
		//std::cout << "[ SHIFITNG ] delta r = [" << deltar.x << "," << deltar.y << "]" << std::endl;

	}

}

