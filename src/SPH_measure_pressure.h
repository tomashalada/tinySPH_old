#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "testing.h"

real Pressure_in_point
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{
	real P;

	P = Kernel_velocity_approximation_TEST(particles, simulation_data, ar)[1];
	//P = Kernel_velocity_approximation_TEST_noWall(particles, simulation_data, ar)[1];

	return P;
}

void Output_hydrostatic_pressure
(Particle_system &particles, std::string file_name)
{
	std::ofstream filePH;
	filePH.open(file_name);

	real H_coef = 1./0.5;
	real p_coef = 1./(9.81*1000.*0.5);

	for(int i = 0; i < particles.np; i++)
	{
		if(particles.data.part_type[i] != fluid){continue;}
		//filePH << particles.data.r[i].y*H_coef << " " << particles.data.p[i]*p_coef << std::endl;
		filePH << particles.data.p[i]*p_coef << " " << particles.data.r[i].y*H_coef << std::endl;
	}

	filePH.close();
}
