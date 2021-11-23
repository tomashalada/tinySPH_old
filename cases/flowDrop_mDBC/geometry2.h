#include "SPH_defs.h"
#include "SPH_particle_system.h"

		//Create_Wall_x_bt(particles, 0. - dp*2, 2. + dp*2, 0. , -1, {0., 1.});
		Create_Wall_x_mdbc(particles, 0.+dp , 1.61 , 0. , -1, dp/2);

		//DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.1 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);


		//DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.}, inlet1);
		DRAW_CreateInletX(particles, 0.3  + 4*dp, 0.8 + dp, 1., +1, 4, {0.0, -3.0}, inlet1);
		DRAW_CreateOutletY(particles, 1.3  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.}, outlet1, true);
		//DRAW_CreateOutletY(particles, 0.5  - 4*dp, 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet2, true);
		DRAW_CreateOutletY(particles, 0.5 , 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet2, true);

		//DRAW_CreateInletY(particles, 0.9  - 4*dp, 0.0 + dp, 0.1, +1, 4, {-4.0, 0.}, inlet2);
		//DRAW_CreateOutletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet, true);
