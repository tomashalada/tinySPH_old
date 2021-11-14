#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "testing.h"

real Pressure_in_point
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{
	real P;

	P = Kernel_velocity_approximation_TEST(particles, simulation_data, ar)[1];

	return P;
}
