#include "SPH_defs.h"
#include "SPH_particle_system.h"

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
						//disabled for mDBC Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);

						//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
						//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.3 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);
