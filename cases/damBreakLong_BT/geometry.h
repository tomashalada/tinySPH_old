#include "SPH_defs.h"
#include "SPH_particle_system.h"


		// Create_Wall_x_mdbc(particles, 0.+dp , 0.8 , 0. , -1, dp/2);
		// //Create_Wall_x_lattice(particles, 0.+dp , 0.8 - dp , 0. , -1);
		// Create_Wall_y_mdbc(particles, 0. , 0 + dp, 1. , -1, dp/2);
		// Create_Wall_y_mdbc(particles, 0.8 , 0 + dp, 1. , +1, 0.8-dp/2);

		// Create_corner_mdbc(particles, -1, -1, {0. + dp/2, 0. + dp/2});
		// Create_corner_mdbc(particles, +1, -1, {0.8 - dp/2, 0. + dp/2});

		// //Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.2, 0.5);
		// Create_Block_of_FluidWP(particles, 0 + dp, 0 + dp, 0.2, 0.5);

		Create_Wall_x_bt(particles, 0.+dp , 1.61 -dp, 0. , -1, {0., 1.});
		Create_Wall_y_bt(particles, 0. , 0. + dp, 1. , -1, {1., 0.});
		Create_Wall_y_bt(particles, 1.61 , 0. + dp , 1. , +1, {-1. , 0.});

		Create_corner_bt(particles, -1, -1, {0., 0.}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, +1, -1, {1.61, 0.}, {- sqrt(2), sqrt(2)});

		//Create_Block_of_Fluid(particles, 0.3 + dp, 0 + dp, 0.5, 0.5);
		Create_Block_of_FluidWP(particles, 0.0 + dp, 0 + dp, 0.6, 0.3);
