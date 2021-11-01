#include "SPH_defs.h"
#include "SPH_particle_system.h"


		//Create_Wall_x_mdbc(particles, 0.0, 0.4-3*dp, 0.1-dp/3 , -1,0.1-dp/3 + dp/2);
		Create_Wall_x_mdbc(particles, 0.0, 0.4-3*dp, 0.1 , -1,0.1 + dp/2);
		Create_Wall_y_mdbc(particles, 0.4 , 0 + dp, 0.1- 3*dp , -1, 0.4 + dp/2);
		Create_Wall_x_mdbc(particles, 0.4+dp , 1.2 , 0. , -1, 0 + dp/2);


		Create_corner_mdbc(particles, -1, -1, {0.4 + dp/2, 0. + dp/2});

		Create_corner_mdbc(particles, -1, -1, {0.4 + dp/2, 0.1 + dp/2});

		DRAW_CreateInletY(particles, 0.1 , 0.1 + dp, 0.2, -1, 4, {0.75, 0.}, inlet1);
		//DRAW_CreateInletY(particles, 0.6  + dp, 0.0 + dp, 0.1, -1, 4, {0.5, 0.}, inlet1);
		DRAW_CreateOutletY(particles, 1.0, 0.0 + dp, 0.5, +1, 4, {0.5, 0.}, outlet1, true);


		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.2, 0.5);
		//Create_Block_of_FluidWP(particles, 0 + dp, 0 + dp, 0.2, 0.5);

