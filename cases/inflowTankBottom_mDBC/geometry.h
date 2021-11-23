#include "SPH_defs.h"
#include "SPH_particle_system.h"

		//BOX
		Create_Wall_x_mdbc(particles, 0.2+dp , 1.2-dp, 0. , -1, dp/2);
		Create_Wall_y_mdbc(particles, 0.2 , 0 + dp, 0.4 , -1, dp/2);
		Create_Wall_y_mdbc(particles, 1.2 , 0 + dp, 0.8 , +1, 1.2-dp/2);

		Create_corner_mdbc(particles, -1, -1, {0.2 + dp/2, 0. + dp/2});
		Create_corner_mdbc(particles, +1, -1, {1.2 - dp/2, 0. + dp/2});

		//INFLOW PIPE
		Create_Wall_x_mdbc(particles, 0.0+dp , 0.4, 0.4+3*dp , -1, dp/2);

		//BLOCK OF WATER
		Create_Block_of_FluidWP(particles, 0.2 + dp, 0 + dp, 1.2, 0.15);

		//INLET
		DRAW_CreateInletY(particles, 0.1 , 0.4 + 4*dp, 0.5, -1, 4, {1.5, 0.}, inlet1);

		//: //DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.}, inlet1);
		//: DRAW_CreateInletX(particles, 0.2  + 4*dp, 0.8 + dp, 1., +1, 4, {0.0, -1.5}, inlet1);
		//: DRAW_CreateOutletY(particles, 1.6  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.}, outlet1, true);
		//: //DRAW_CreateOutletY(particles, 0.5  - 4*dp, 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet2, true);
		//: DRAW_CreateOutletY(particles, 0.2 , 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet2, true);

		//: //DRAW_CreateInletY(particles, 0.9  - 4*dp, 0.0 + dp, 0.1, +1, 4, {-4.0, 0.}, inlet2);
		//: //DRAW_CreateOutletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet, true);
